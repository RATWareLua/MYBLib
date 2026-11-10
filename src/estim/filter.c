/* SPDX-License-Identifier: BSD-3-Clause */
/* The extended Kalman filter on the state, the parameters of a force model and the systematic
 * shifts of the measurements; the fixed-interval smoother over what the filter recorded. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/estim.h"
#include "myb/num.h"
#include "myb/orbit.h"
#include "myb/uncert.h"

static const double pi = 3.141592653589793;

struct myb_filter {
  myb_allocator allocator;
  myb_force_model* model;
  myb_orbit* orbit;
  size_t np;
  size_t nb;
  size_t n; /* 6 + np + nb */
  double q;
  int started;
  int running; /* the propagator is started at (t, y) */
  double t;
  double t_end;
  double y[6];
  double* bias; /* nb */
  double* p;    /* n x n */
  double* f;    /* n x n: transition, then I - K H */
  double* next; /* n x n */
  double* pht;  /* n x MYB_MEASURE_MAX */
  double* gain; /* n x MYB_MEASURE_MAX */
  /* the record for the smoother: one node per instant */
  size_t capacity;
  size_t nodes;
  int smoothed;
  double* node_t;  /* capacity */
  double* x_minus; /* capacity x n: before the observations of the instant */
  double* x_plus;  /* capacity x n: after them; the smoothed values after myb_filter_smooth */
  double* p_minus; /* capacity x n x n */
  double* p_plus;  /* capacity x n x n */
  double* step;    /* capacity x n x n: the transition from the node before */
  double* work;    /* 3 n x n + n */
};

void myb_filter_destroy(myb_filter* filter) {
  if (filter == NULL) {
    return;
  }
  const myb_allocator a = filter->allocator;
  myb_orbit_destroy(filter->orbit);
  myb_internal_release(&a, filter->p);
  myb_internal_release(&a, filter->node_t);
  myb_internal_release(&a, filter);
}

myb_status myb_filter_create(myb_force_model* model, const myb_filter_options* options,
                             const myb_allocator* allocator, myb_filter** filter) {
  MYB_CHECK(model != NULL && filter != NULL);
  const myb_filter_options none = {0};
  const myb_filter_options* opt = options != NULL ? options : &none;
  if (!(opt->epsilon >= 0.0) || !(opt->process_noise >= 0.0) || !isfinite(opt->epsilon) ||
      !isfinite(opt->process_noise) || opt->bias_count > MYB_FIT_MAX_BIASES) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_filter), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_filter* f = block;
  memset(f, 0, sizeof *f);
  f->allocator = a;
  f->model = model;
  f->np = myb_force_model_parameter_count(model);
  f->nb = opt->bias_count;
  f->n = 6 + f->np + f->nb;
  f->q = opt->process_noise;
  f->capacity = opt->history;
  const size_t n = f->n;
  status = myb_internal_alloc_array(&a, 3 * n * n + 2 * n * MYB_MEASURE_MAX + f->nb + 1,
                                    sizeof(double), &block);
  if (status == MYB_OK) {
    f->p = block;
    f->f = f->p + n * n;
    f->next = f->f + n * n;
    f->pht = f->next + n * n;
    f->gain = f->pht + n * MYB_MEASURE_MAX;
    f->bias = f->gain + n * MYB_MEASURE_MAX;
    memset(f->bias, 0, (f->nb + 1) * sizeof(double));
    const myb_orbit_options orbit_options = {.epsilon = opt->epsilon, .variational = 1};
    status = myb_orbit_create(model, &orbit_options, allocator, &f->orbit);
  }
  if (status == MYB_OK && f->capacity > 0) {
    /* per node: the time, two vectors and three matrices; then the work of the smoother */
    const size_t per_node = 1 + 2 * n + 3 * n * n;
    if (f->capacity > (SIZE_MAX / sizeof(double) - (3 * n * n + n)) / per_node) {
      status = MYB_E_OVERFLOW;
    } else {
      status = myb_internal_alloc_array(&a, f->capacity * per_node + 3 * n * n + n,
                                        sizeof(double), &block);
    }
    if (status == MYB_OK) {
      f->node_t = block;
      f->x_minus = f->node_t + f->capacity;
      f->x_plus = f->x_minus + f->capacity * n;
      f->p_minus = f->x_plus + f->capacity * n;
      f->p_plus = f->p_minus + f->capacity * n * n;
      f->step = f->p_plus + f->capacity * n * n;
      f->work = f->step + f->capacity * n * n;
    }
  }
  if (status != MYB_OK) {
    myb_filter_destroy(f);
    return status;
  }
  *filter = f;
  return MYB_OK;
}

/* Starts the propagator at the filter's time and state, unless the end is reached. */
static myb_status restart(myb_filter* f) {
  f->running = 0;
  if (f->t == f->t_end) {
    return MYB_OK;
  }
  const myb_status status = myb_orbit_start(f->orbit, f->t, f->y, f->t_end);
  if (status == MYB_OK) {
    f->running = 1;
  }
  return status;
}

/* The unknowns as one vector: the state, the parameters of the model, the shifts. */
static void pack(const myb_filter* f, double* x) {
  memcpy(x, f->y, sizeof f->y);
  for (size_t k = 0; k < f->np; ++k) {
    x[6 + k] = *myb_force_model_parameter(f->model, k);
  }
  for (size_t k = 0; k < f->nb; ++k) {
    x[6 + f->np + k] = f->bias[k];
  }
}

/* Writes the filter's present estimate into the node of the present instant. */
static void record_plus(myb_filter* f) {
  if (f->capacity == 0) {
    return;
  }
  const size_t n = f->n;
  const size_t k = f->nodes - 1;
  pack(f, f->x_plus + k * n);
  memcpy(f->p_plus + k * n * n, f->p, n * n * sizeof(double));
}

myb_status myb_filter_start(myb_filter* filter, double t0, const double y0[6],
                            const double* covariance, double t_end) {
  MYB_CHECK(filter != NULL && y0 != NULL && covariance != NULL);
  myb_filter* f = filter;
  f->started = 0;
  if (!isfinite(t0) || !isfinite(t_end) || !(t_end > t0)) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < f->n * f->n; ++i) {
    if (!isfinite(covariance[i])) {
      return MYB_E_DOMAIN;
    }
  }
  for (int i = 0; i < 6; ++i) {
    if (!isfinite(y0[i])) {
      return MYB_E_DOMAIN;
    }
  }
  memcpy(f->y, y0, sizeof f->y);
  memcpy(f->p, covariance, f->n * f->n * sizeof(double));
  memset(f->bias, 0, f->nb * sizeof(double));
  f->t = t0;
  f->t_end = t_end;
  f->nodes = 0;
  f->smoothed = 0;
  if (f->capacity > 0) {
    const size_t n = f->n;
    f->nodes = 1;
    f->node_t[0] = t0;
    pack(f, f->x_minus);
    memcpy(f->p_minus, f->p, n * n * sizeof(double));
    memset(f->step, 0, n * n * sizeof(double));
    for (size_t i = 0; i < n; ++i) {
      f->step[i * n + i] = 1.0;
    }
    record_plus(f);
  }
  const myb_status status = restart(f);
  if (status == MYB_OK) {
    f->started = 1;
  }
  return status;
}

myb_status myb_filter_set_bias(myb_filter* filter, const double* bias) {
  MYB_CHECK(filter != NULL && filter->started && (bias != NULL || filter->nb == 0));
  for (size_t k = 0; k < filter->nb; ++k) {
    if (!isfinite(bias[k])) {
      return MYB_E_DOMAIN;
    }
  }
  if (filter->nb > 0) {
    memcpy(filter->bias, bias, filter->nb * sizeof(double));
  }
  if (filter->capacity > 0 && filter->nodes == 1) {
    pack(filter, filter->x_minus);
    record_plus(filter);
  }
  return MYB_OK;
}

void myb_filter_bias(const myb_filter* filter, double* bias) {
  MYB_CHECK(filter != NULL && filter->started && (bias != NULL || filter->nb == 0));
  if (filter->nb > 0) {
    memcpy(bias, filter->bias, filter->nb * sizeof(double));
  }
}

myb_status myb_filter_predict(myb_filter* filter, double t) {
  MYB_CHECK(filter != NULL && filter->started);
  myb_filter* f = filter;
  if (!isfinite(t) || t < f->t || t > f->t_end) {
    return MYB_E_DOMAIN;
  }
  if (t == f->t) {
    return MYB_OK;
  }
  if (f->capacity > 0 && f->nodes == f->capacity) {
    return MYB_E_OVERFLOW;
  }
  const size_t n = f->n;
  const size_t np = f->np;
  double y[6];
  double phi[36];
  double sens[6 * MYB_MODEL_MAX_PARAMETERS];
  myb_status status = myb_orbit_partials(f->orbit, t, y, phi, np > 0 ? sens : NULL);
  if (status != MYB_OK) {
    return status;
  }
  /* F = | Phi S 0 |
   *     |  0  1 0 |
   *     |  0  0 1 | */
  memset(f->f, 0, n * n * sizeof(double));
  for (size_t i = 0; i < 6; ++i) {
    for (size_t j = 0; j < 6; ++j) {
      f->f[i * n + j] = phi[6 * i + j];
    }
    for (size_t k = 0; k < np; ++k) {
      f->f[i * n + 6 + k] = sens[i * np + k];
    }
  }
  for (size_t k = 6; k < n; ++k) {
    f->f[k * n + k] = 1.0;
  }
  myb_covariance_transform(n, f->f, f->p, f->next);
  const double dt = t - f->t;
  for (size_t i = 0; i < 3; ++i) {
    f->next[i * n + i] += f->q * dt * dt * dt / 3.0;
    f->next[i * n + 3 + i] += f->q * dt * dt / 2.0;
    f->next[(3 + i) * n + i] += f->q * dt * dt / 2.0;
    f->next[(3 + i) * n + 3 + i] += f->q * dt;
  }
  memcpy(f->p, f->next, n * n * sizeof(double));
  memcpy(f->y, y, sizeof f->y);
  f->t = t;
  if (f->capacity > 0) {
    const size_t k = f->nodes++;
    f->smoothed = 0;
    f->node_t[k] = t;
    memcpy(f->step + k * n * n, f->f, n * n * sizeof(double));
    pack(f, f->x_minus + k * n);
    memcpy(f->p_minus + k * n * n, f->p, n * n * sizeof(double));
    record_plus(f);
  }
  return restart(f);
}

myb_status myb_filter_update(myb_filter* filter, const myb_observation* observation,
                             double* innovation) {
  MYB_CHECK(filter != NULL && filter->started && observation != NULL);
  myb_filter* f = filter;
  const myb_observation* o = observation;
  if (o->measure == NULL || o->count < 1 || o->count > MYB_MEASURE_MAX || o->bias < 0 ||
      (o->bias > 0 && (size_t)o->bias - 1 + o->count > f->nb)) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < o->count; ++i) {
    if (!isfinite(o->value[i]) || !(o->sigma[i] > 0.0) || !isfinite(o->sigma[i])) {
      return MYB_E_DOMAIN;
    }
  }
  myb_status status = myb_filter_predict(f, o->t);
  if (status != MYB_OK) {
    return status;
  }
  const size_t n = f->n;
  const size_t m = o->count;
  /* the column of the unknowns that holds the shift of value i, or n for none */
  const size_t shift = o->bias > 0 ? 6 + f->np + (size_t)o->bias - 1 : n;
  double computed[MYB_MEASURE_MAX];
  double h[6 * MYB_MEASURE_MAX];
  status = o->measure(o->context, o->t, f->y, computed, h);
  if (status != MYB_OK) {
    return status;
  }
  double r[MYB_MEASURE_MAX];
  for (size_t i = 0; i < m; ++i) {
    if (shift < n) {
      computed[i] += f->bias[(size_t)o->bias - 1 + i];
    }
    r[i] = o->value[i] - computed[i];
    if ((o->periodic >> i) & 1u) {
      r[i] -= 2.0 * pi * floor((r[i] + pi) / (2.0 * pi));
    }
    if (!isfinite(r[i])) {
      return MYB_E_DOMAIN;
    }
    if (innovation != NULL) {
      innovation[i] = r[i];
    }
  }
  /* P H^T (the measurement sees the state and its own shift) and S = H P H^T + R */
  for (size_t a = 0; a < n; ++a) {
    for (size_t i = 0; i < m; ++i) {
      double sum = 0.0;
      for (size_t j = 0; j < 6; ++j) {
        sum += f->p[a * n + j] * h[6 * i + j];
      }
      if (shift < n) {
        sum += f->p[a * n + shift + i];
      }
      f->pht[a * MYB_MEASURE_MAX + i] = sum;
    }
  }
  double s[MYB_MEASURE_MAX * MYB_MEASURE_MAX];
  for (size_t i = 0; i < m; ++i) {
    for (size_t k = 0; k < m; ++k) {
      double sum = 0.0;
      for (size_t j = 0; j < 6; ++j) {
        sum += h[6 * i + j] * f->pht[j * MYB_MEASURE_MAX + k];
      }
      if (shift < n) {
        sum += f->pht[(shift + i) * MYB_MEASURE_MAX + k];
      }
      s[i * m + k] = sum;
    }
    s[i * m + i] += o->sigma[i] * o->sigma[i];
  }
  /* symmetric to rounding: the factorisation reads the lower triangle */
  for (size_t i = 0; i < m; ++i) {
    for (size_t k = 0; k < i; ++k) {
      s[i * m + k] = 0.5 * (s[i * m + k] + s[k * m + i]);
    }
  }
  status = myb_cholesky(m, s);
  if (status != MYB_OK) {
    return status;
  }
  /* K = P H^T S^-1, row by row */
  for (size_t a = 0; a < n; ++a) {
    double row[MYB_MEASURE_MAX];
    for (size_t i = 0; i < m; ++i) {
      row[i] = f->pht[a * MYB_MEASURE_MAX + i];
    }
    myb_cholesky_solve(m, s, row);
    for (size_t i = 0; i < m; ++i) {
      f->gain[a * MYB_MEASURE_MAX + i] = row[i];
    }
  }
  /* the state, the parameters and the shifts */
  for (size_t a = 0; a < n; ++a) {
    double step = 0.0;
    for (size_t i = 0; i < m; ++i) {
      step += f->gain[a * MYB_MEASURE_MAX + i] * r[i];
    }
    if (a < 6) {
      f->y[a] += step;
    } else if (a < 6 + f->np) {
      *myb_force_model_parameter(f->model, a - 6) += step;
    } else {
      f->bias[a - 6 - f->np] += step;
    }
  }
  /* Joseph's form: P = (1 - K H) P (1 - K H)^T + K R K^T */
  memset(f->f, 0, n * n * sizeof(double));
  for (size_t a = 0; a < n; ++a) {
    f->f[a * n + a] = 1.0;
    for (size_t j = 0; j < 6; ++j) {
      double sum = 0.0;
      for (size_t i = 0; i < m; ++i) {
        sum += f->gain[a * MYB_MEASURE_MAX + i] * h[6 * i + j];
      }
      f->f[a * n + j] -= sum;
    }
    if (shift < n) {
      for (size_t i = 0; i < m; ++i) {
        f->f[a * n + shift + i] -= f->gain[a * MYB_MEASURE_MAX + i];
      }
    }
  }
  myb_covariance_transform(n, f->f, f->p, f->next);
  for (size_t a = 0; a < n; ++a) {
    for (size_t b = 0; b < n; ++b) {
      double sum = 0.0;
      for (size_t i = 0; i < m; ++i) {
        sum += f->gain[a * MYB_MEASURE_MAX + i] * o->sigma[i] * o->sigma[i] *
               f->gain[b * MYB_MEASURE_MAX + i];
      }
      f->next[a * n + b] += sum;
    }
  }
  for (size_t a = 0; a < n; ++a) {
    for (size_t b = 0; b <= a; ++b) {
      const double mean = 0.5 * (f->next[a * n + b] + f->next[b * n + a]);
      f->p[a * n + b] = mean;
      f->p[b * n + a] = mean;
    }
  }
  f->smoothed = 0;
  record_plus(f);
  return restart(f);
}

void myb_filter_get(const myb_filter* filter, double* t, double y[6], double* covariance) {
  MYB_CHECK(filter != NULL && filter->started);
  if (t != NULL) {
    *t = filter->t;
  }
  if (y != NULL) {
    memcpy(y, filter->y, sizeof filter->y);
  }
  if (covariance != NULL) {
    memcpy(covariance, filter->p, filter->n * filter->n * sizeof(double));
  }
}

/* ---- The smoother ----------------------------------------------------------------------------
 * H. E. Rauch, F. Tung, C. T. Striebel, AIAA J. 3, 1445 (1965): from the last node backwards
 *   C_k = P+_k F_{k+1}^T (P-_{k+1})^-1,
 *   x^s_k = x+_k + C_k (x^s_{k+1} - x-_{k+1}),
 *   P^s_k = P+_k + C_k (P^s_{k+1} - P-_{k+1}) C_k^T. */

size_t myb_filter_steps(const myb_filter* filter) {
  MYB_CHECK(filter != NULL);
  return filter->nodes;
}

myb_status myb_filter_smooth(myb_filter* filter) {
  MYB_CHECK(filter != NULL && filter->started);
  myb_filter* f = filter;
  if (f->capacity == 0) {
    return MYB_E_DOMAIN;
  }
  if (f->smoothed) {
    return MYB_OK;
  }
  const size_t n = f->n;
  double* factor = f->work;          /* the factor of P-_{k+1} */
  double* c = f->work + n * n;       /* C_k */
  double* d = f->work + 2 * n * n;   /* P^s_{k+1} - P-_{k+1}, then C D C^T */
  double* column = f->work + 3 * n * n;
  for (size_t k = f->nodes - 1; k-- > 0;) {
    const double* p_plus = f->p_plus + k * n * n;
    const double* p_next_minus = f->p_minus + (k + 1) * n * n;
    const double* step = f->step + (k + 1) * n * n;
    memcpy(factor, p_next_minus, n * n * sizeof(double));
    const myb_status status = myb_cholesky(n, factor);
    if (status != MYB_OK) {
      return status;
    }
    /* row a of C: (P-)^-1 applied to row a of P+ F^T */
    for (size_t a = 0; a < n; ++a) {
      for (size_t b = 0; b < n; ++b) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
          sum += p_plus[a * n + j] * step[b * n + j];
        }
        column[b] = sum;
      }
      myb_cholesky_solve(n, factor, column);
      memcpy(c + a * n, column, n * sizeof(double));
    }
    double* x = f->x_plus + k * n;
    const double* x_next = f->x_plus + (k + 1) * n;
    const double* x_next_minus = f->x_minus + (k + 1) * n;
    for (size_t a = 0; a < n; ++a) {
      double sum = 0.0;
      for (size_t b = 0; b < n; ++b) {
        sum += c[a * n + b] * (x_next[b] - x_next_minus[b]);
      }
      x[a] += sum;
    }
    const double* p_next = f->p_plus + (k + 1) * n * n;
    for (size_t i = 0; i < n * n; ++i) {
      d[i] = p_next[i] - p_next_minus[i];
    }
    myb_covariance_transform(n, c, d, factor);
    double* p = f->p_plus + k * n * n;
    for (size_t a = 0; a < n; ++a) {
      for (size_t b = 0; b <= a; ++b) {
        const double mean = p[a * n + b] + 0.5 * (factor[a * n + b] + factor[b * n + a]);
        p[a * n + b] = mean;
        p[b * n + a] = mean;
      }
    }
  }
  f->smoothed = 1;
  return MYB_OK;
}

void myb_filter_step(const myb_filter* filter, size_t step, double* t, double* unknowns,
                     double* covariance) {
  MYB_CHECK(filter != NULL && step < filter->nodes);
  const size_t n = filter->n;
  if (t != NULL) {
    *t = filter->node_t[step];
  }
  if (unknowns != NULL) {
    memcpy(unknowns, filter->x_plus + step * n, n * sizeof(double));
  }
  if (covariance != NULL) {
    memcpy(covariance, filter->p_plus + step * n * n, n * n * sizeof(double));
  }
}
