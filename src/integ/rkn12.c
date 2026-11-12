/* SPDX-License-Identifier: BSD-3-Clause */
/* RKN12(10): J. R. Dormand, M. E. A. El-Mikkawy, P. J. Prince, IMA J. Numer. Anal. 7, 423
 * (1987), for x'' = a(t, x), with the coefficients of integ/rkn12_coefficients.h. Step-size
 * control as in E. Hairer, S. P. Norsett, G. Wanner, Solving ODE I (1993), §II.4, with the
 * exponent 1/11 of the order-10 error estimate; the estimate covers positions and velocities.
 * The pair has no interpolant of its own: the state inside the last step is one more step of the
 * pair from its start (16 evaluations), or, with options.interpolant, comes from one built on
 * the first request in the step: the acceleration as a polynomial of degree 12 through its
 * values at the 13 Chebyshev-Lobatto nodes of the step, integrated twice from the start (a
 * Chebyshev series integrated term by term), the values at the 11 inner nodes taken at the
 * places of the previous polynomial — Picard's iteration started from the quintic of Hermite
 * of the ends — until the places at the nodes move less than a thousandth of the tolerance,
 * then made to end exactly at the step's end by a term in theta^2, theta^3 of the size of the
 * mismatch. Should it not settle in 12 rounds (a step too long for the degree), the state is
 * one more step of the pair from the start of the step (16 evaluations a request). Positions,
 * velocities and time are updated with compensated summation. */
#include <math.h> /* sqrt, fabs */
#include <stdckdint.h>
#include <stdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "integ/event.h"
#include "integ/rkn12_coefficients.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "num/eft.h"

enum {
  MAX_REJECTIONS = 100,
  DEGREE = 12,           /* of the polynomial of the acceleration over a step */
  NODES = DEGREE + 1,    /* Chebyshev-Lobatto, both ends included */
  ROUNDS = 12            /* of Picard's iteration before falling back to a step of the pair */
};

static const int64_t default_max_steps = 1000000000;
static const double safety = 0.9;
static const double min_factor = 0.2;
static const double max_factor = 5.0;

struct myb_rkn12 {
  myb_allocator allocator;
  size_t n;
  myb_rkn12_options options;
  myb_position_accel_function f;
  void* context;
  double t;
  double t_comp;
  double t_old;
  double t_comp_old;
  double h_last;
  double h; /* next step (signed); 0 = automatic */
  int started;
  int have_step;
  myb_rkn12_stats stats;
  double* block;
  double* y;        /* 2n: x, v */
  double* y_comp;   /* 2n */
  double* y_old;    /* 2n */
  double* comp_old; /* 2n */
  double* y_new;    /* 2n */
  double* comp_new; /* 2n */
  double* error;    /* 2n: error estimate of the last attempt */
  double* scratch;  /* 2n: events */
  double* a_now;    /* n: a(t, x) */
  double* a_old;    /* n: a at the start of the last accepted step */
  double* stage_x;  /* n */
  double* k;        /* MYB_RKN12_STAGES x n */
  /* the interpolant of the last step */
  int dense_state;          /* 0 not built, 1 built, -1 did not settle */
  double node[NODES];       /* theta of the nodes, from 0 to 1 */
  double dct[NODES][NODES]; /* Chebyshev coefficients from the values at the nodes */
  double* node_x;           /* NODES x n: places at the nodes */
  double* node_a;           /* NODES x n: accelerations there */
  double* series_v;         /* (DEGREE + 2) x n: V(tau) = int_-1^tau a */
  double* series_x;         /* (DEGREE + 3) x n: X(tau) = int_-1^tau V */
  double* fix;              /* 2 x n: the terms alpha theta^2 + beta theta^3 */
};

static int finite(double x) {
  return x - x == 0.0;
}

static double power(double x, double p) { /* x > 0 */
  return myb_exp(myb_log(x) * p);
}

static double* stage(const myb_rkn12* s, int i) {
  return s->k + (size_t)i * s->n;
}

myb_status myb_rkn12_create(size_t n, const myb_rkn12_options* options,
                            const myb_allocator* allocator, myb_rkn12** integrator) {
  MYB_CHECK(options != NULL && integrator != NULL);
  myb_rkn12_options o = *options;
  if (n == 0 || !(o.rtol > 0.0) || !finite(o.rtol) || !(o.atol > 0.0) || !finite(o.atol) ||
      !(o.h_initial >= 0.0) || !finite(o.h_initial) || !(o.h_max >= 0.0) || !(o.h_min >= 0.0) ||
      o.max_steps < 0) {
    return MYB_E_DOMAIN;
  }
  if (o.max_steps == 0) {
    o.max_steps = default_max_steps;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_rkn12), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_rkn12* s = block;
  *s = (myb_rkn12){.allocator = a, .n = n, .options = o};
  /* 8 vectors of 2n, 3 of n, 17 stages, the interpolant: 2 NODES + 2 DEGREE + 7 */
  size_t count = 0;
  if (ckd_mul(&count, (size_t)(8 * 2 + 3 + MYB_RKN12_STAGES + 2 * NODES + 2 * DEGREE + 7), n)) {
    status = MYB_E_OVERFLOW;
    goto fail;
  }
  status = myb_internal_alloc_array(&a, count, sizeof(double), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  double* p = block;
  memset(p, 0, count * sizeof(double));
  s->block = p;
  double** wide[8] = {&s->y,     &s->y_comp,   &s->y_old, &s->comp_old,
                      &s->y_new, &s->comp_new, &s->error, &s->scratch};
  for (int i = 0; i < 8; ++i) {
    *wide[i] = p;
    p += 2 * n;
  }
  double** narrow[3] = {&s->a_now, &s->a_old, &s->stage_x};
  for (int i = 0; i < 3; ++i) {
    *narrow[i] = p;
    p += n;
  }
  s->k = p;
  p += (size_t)MYB_RKN12_STAGES * n;
  s->node_x = p;
  p += (size_t)NODES * n;
  s->node_a = p;
  p += (size_t)NODES * n;
  s->series_v = p;
  p += (size_t)(DEGREE + 2) * n;
  s->series_x = p;
  p += (size_t)(DEGREE + 3) * n;
  s->fix = p;
  /* the nodes theta_j = (1 - cos(j pi / M)) / 2 and the discrete cosine transform that gives the
   * coefficients of sum c_k T_k(tau), tau = 2 theta - 1, T_k(tau_j) = (-1)^k cos(k j pi / M) */
  const double pi = 3.14159265358979323846;
  for (int j = 0; j < NODES; ++j) {
    s->node[j] = j == DEGREE ? 1.0 : 0.5 * (1.0 - myb_cos(pi * (double)j / (double)DEGREE));
  }
  for (int k = 0; k < NODES; ++k) {
    for (int j = 0; j < NODES; ++j) {
      const int phase = (k * j) % (2 * DEGREE);
      double c = myb_cos(pi * (double)phase / (double)DEGREE);
      c *= (k % 2 == 0 ? 1.0 : -1.0) * 2.0 / (double)DEGREE;
      if (j == 0 || j == DEGREE) {
        c *= 0.5;
      }
      if (k == 0 || k == DEGREE) {
        c *= 0.5;
      }
      s->dct[k][j] = c;
    }
  }
  *integrator = s;
  return MYB_OK;

fail:
  myb_internal_release(&a, s);
  return status;
}

void myb_rkn12_destroy(myb_rkn12* integrator) {
  if (integrator == NULL) {
    return;
  }
  const myb_allocator a = integrator->allocator;
  myb_internal_release(&a, integrator->block);
  myb_internal_release(&a, integrator);
}

static myb_status evaluate(myb_rkn12* s, double t, const double* x, double* a) {
  s->stats.evaluations++;
  return s->f(s->context, t, x, a);
}

static double rms(size_t m, const double* v) {
  double sum = 0.0;
  for (size_t i = 0; i < m; ++i) {
    sum += v[i] * v[i];
  }
  return sqrt(sum / (double)m);
}

/* Hairer's starting step for the first-order form y = (x, v), y' = (v, a), order 12. */
static myb_status initial_step(myb_rkn12* s, double direction, double* h_out) {
  const size_t n = s->n;
  double* scaled = s->scratch;
  double* scale = s->error;
  for (size_t i = 0; i < 2 * n; ++i) {
    scale[i] = s->options.atol + s->options.rtol * fabs(s->y[i]);
  }
  for (size_t i = 0; i < 2 * n; ++i) {
    scaled[i] = s->y[i] / scale[i];
  }
  const double d0 = rms(2 * n, scaled);
  for (size_t i = 0; i < n; ++i) {
    scaled[i] = s->y[n + i] / scale[i];
    scaled[n + i] = s->a_now[i] / scale[n + i];
  }
  const double d1 = rms(2 * n, scaled);
  const double h0 = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1;
  /* Euler step of length h0 and the change of the derivative over it */
  for (size_t i = 0; i < n; ++i) {
    s->stage_x[i] = s->y[i] + direction * h0 * s->y[n + i];
  }
  double* a1 = stage(s, 1);
  const myb_status status = evaluate(s, s->t + direction * h0, s->stage_x, a1);
  if (status != MYB_OK) {
    return status;
  }
  for (size_t i = 0; i < n; ++i) {
    scaled[i] = direction * h0 * s->a_now[i] / scale[i]; /* change of v over h0 */
    scaled[n + i] = (a1[i] - s->a_now[i]) / scale[n + i];
  }
  const double d2 = rms(2 * n, scaled) / h0;
  const double larger = d1 > d2 ? d1 : d2;
  const double h1 = larger <= 1e-15 ? (h0 * 1e-3 > 1e-6 ? h0 * 1e-3 : 1e-6)
                                    : power(0.01 / larger, 1.0 / 13.0);
  *h_out = direction * (100.0 * h0 < h1 ? 100.0 * h0 : h1);
  return MYB_OK;
}

myb_status myb_rkn12_start(myb_rkn12* integrator, myb_position_accel_function a, void* context,
                           double t0, const double* y0) {
  MYB_CHECK(integrator != NULL && a != NULL && y0 != NULL);
  myb_rkn12* s = integrator;
  if (!finite(t0)) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < 2 * s->n; ++i) {
    if (!finite(y0[i])) {
      return MYB_E_DOMAIN;
    }
    s->y[i] = y0[i];
    s->y_comp[i] = 0.0;
  }
  s->f = a;
  s->context = context;
  s->t = t0;
  s->t_comp = 0.0;
  s->h = s->options.h_initial;
  s->have_step = 0;
  s->dense_state = 0;
  s->stats = (myb_rkn12_stats){};
  const myb_status status = evaluate(s, t0, s->y, s->a_now);
  s->started = status == MYB_OK;
  return status;
}

/* One step of length h from (t + t_comp, y, comp) with a0 = a(t, x): stages, the order-12
 * state into out/out_comp (compensated) and, if error != NULL, the error estimate. */
static myb_status pair_step(myb_rkn12* s, double t, double t_comp, const double* y,
                            const double* comp, const double* a0, double h, double* out,
                            double* out_comp, double* error) {
  const size_t n = s->n;
  memcpy(stage(s, 0), a0, n * sizeof(double));
  for (int st = 1; st < MYB_RKN12_STAGES; ++st) {
    const double* a = myb_rkn12_a[st];
    const double ch = myb_rkn12_c[st] * h;
    for (size_t i = 0; i < n; ++i) {
      double sum = 0.0;
      for (int j = 0; j < st; ++j) {
        sum += a[j] * stage(s, j)[i];
      }
      s->stage_x[i] = y[i] + (comp[i] + (ch * y[n + i] + h * h * sum));
    }
    const myb_status status = evaluate(s, t + (t_comp + ch), s->stage_x, stage(s, st));
    if (status != MYB_OK) {
      return status;
    }
  }
  for (size_t i = 0; i < n; ++i) {
    double sb = 0.0;
    double sbp = 0.0;
    double sbt = 0.0;
    double sbpt = 0.0;
    for (int j = 0; j < MYB_RKN12_STAGES; ++j) {
      const double kj = stage(s, j)[i];
      sb += myb_rkn12_b[j] * kj;
      sbp += myb_rkn12_bp[j] * kj;
      sbt += myb_rkn12_bt[j] * kj;
      sbpt += myb_rkn12_bpt[j] * kj;
    }
    const double dx = comp[i] + h * (y[n + i] + h * sb);
    const double dv = comp[n + i] + h * sbp;
    if (out_comp != NULL) {
      out[i] = myb_eft_two_sum(y[i], dx, &out_comp[i]);
      out[n + i] = myb_eft_two_sum(y[n + i], dv, &out_comp[n + i]);
    } else {
      out[i] = y[i] + dx;
      out[n + i] = y[n + i] + dv;
    }
    if (error != NULL) {
      error[i] = h * h * sbt;
      error[n + i] = h * sbpt;
    }
  }
  return MYB_OK;
}

static double error_norm(const myb_rkn12* s) {
  const size_t m = 2 * s->n;
  double sum = 0.0;
  for (size_t i = 0; i < m; ++i) {
    const double y_abs = fabs(s->y[i]) > fabs(s->y_new[i]) ? fabs(s->y[i]) : fabs(s->y_new[i]);
    const double q = s->error[i] / (s->options.atol + s->options.rtol * y_abs);
    sum += q * q;
  }
  return sqrt(sum / (double)m);
}

static double clamp_step(const myb_rkn12* s, double h) {
  if (s->options.h_max > 0.0 && fabs(h) > s->options.h_max) {
    return h > 0.0 ? s->options.h_max : -s->options.h_max;
  }
  return h;
}

myb_status myb_rkn12_step(myb_rkn12* integrator, double t_limit) {
  MYB_CHECK(integrator != NULL && integrator->started);
  myb_rkn12* s = integrator;
  if (!finite(t_limit)) {
    return MYB_E_DOMAIN;
  }
  const double remaining = (t_limit - s->t) - s->t_comp;
  if (remaining == 0.0) {
    return MYB_OK;
  }
  const double direction = remaining > 0.0 ? 1.0 : -1.0;
  if (s->h == 0.0 || s->h * direction < 0.0) {
    double h0 = 0.0;
    const myb_status status = initial_step(s, direction, &h0);
    if (status != MYB_OK) {
      return status;
    }
    s->h = h0;
  }
  double h = clamp_step(s, s->h);
  int previous_rejected = 0;
  /* bound: MAX_REJECTIONS attempts; every rejection shrinks h by a factor <= 0.9 */
  for (int tries = 0; tries < MAX_REJECTIONS; ++tries) {
    const int last = 1.01 * fabs(h) >= fabs(remaining);
    const double h_planned = h;
    if (last) {
      h = remaining;
    }
    if (!last && s->options.h_min > 0.0 && fabs(h) < s->options.h_min) {
      return MYB_E_STEP_MIN;
    }
    if (!last && s->t + h == s->t) {
      return MYB_E_STEP_MIN;
    }
    myb_status status = pair_step(s, s->t, s->t_comp, s->y, s->y_comp, s->a_now, h, s->y_new,
                                  s->comp_new, s->error);
    if (status != MYB_OK) {
      return status;
    }
    const double norm = error_norm(s);
    if (!finite(norm)) {
      return MYB_E_STEP_MIN;
    }
    if (norm <= 1.0) {
      double factor = norm == 0.0 ? max_factor : safety * power(norm, -1.0 / 11.0);
      factor = factor < max_factor ? factor : max_factor;
      if (previous_rejected && factor > 1.0) {
        factor = 1.0;
      }
      /* The acceleration at the new state: the first stage of the next step. */
      memcpy(s->a_old, s->a_now, s->n * sizeof(double));
      double err = 0.0;
      const double t_new = myb_eft_two_sum(s->t, h + s->t_comp, &err);
      status = evaluate(s, last ? t_limit : t_new + err, s->y_new, s->a_now);
      if (status != MYB_OK) {
        memcpy(s->a_now, s->a_old, s->n * sizeof(double));
        return status;
      }
      s->t_old = s->t;
      s->t_comp_old = s->t_comp;
      memcpy(s->y_old, s->y, 2 * s->n * sizeof(double));
      memcpy(s->comp_old, s->y_comp, 2 * s->n * sizeof(double));
      memcpy(s->y, s->y_new, 2 * s->n * sizeof(double));
      memcpy(s->y_comp, s->comp_new, 2 * s->n * sizeof(double));
      s->t = last ? t_limit : t_new;
      s->t_comp = last ? 0.0 : err;
      s->h_last = h;
      const double h_next = h * factor;
      s->h = clamp_step(s, last && fabs(h_planned) > fabs(h_next) ? h_planned : h_next);
      s->have_step = 1;
      s->dense_state = 0;
      s->stats.accepted++;
      return MYB_OK;
    }
    s->stats.rejected++;
    double factor = safety * power(norm, -1.0 / 11.0);
    factor = factor > min_factor ? factor : min_factor;
    h *= factor;
    previous_rejected = 1;
  }
  return MYB_E_NOCONV;
}

myb_status myb_rkn12_integrate(myb_rkn12* integrator, double t_end) {
  MYB_CHECK(integrator != NULL && integrator->started);
  /* bound: max_steps accepted steps */
  for (int64_t i = 0; i < integrator->options.max_steps; ++i) {
    if (integrator->t == t_end && integrator->t_comp == 0.0) {
      return MYB_OK;
    }
    const myb_status status = myb_rkn12_step(integrator, t_end);
    if (status != MYB_OK) {
      return status;
    }
  }
  return integrator->t == t_end && integrator->t_comp == 0.0 ? MYB_OK : MYB_E_NOCONV;
}

double myb_rkn12_time(const myb_rkn12* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->t + integrator->t_comp;
}

const double* myb_rkn12_state(const myb_rkn12* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->y;
}

myb_rkn12_stats myb_rkn12_statistics(const myb_rkn12* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->stats;
}

/* ---- The interpolant of a step ---------------------------------------------------------- */

/* The sum of a Chebyshev series of `count` coefficients (stride n) at tau (Clenshaw). */
static double chebyshev(const double* c, size_t stride, int count, double tau) {
  double b1 = 0.0;
  double b2 = 0.0;
  for (int k = count - 1; k >= 1; --k) {
    const double b0 = c[(size_t)k * stride] + 2.0 * tau * b1 - b2;
    b2 = b1;
    b1 = b0;
  }
  return c[0] + tau * b1 - b2;
}

/* The coefficients of the integral from -1 of a series of `count` terms: `count + 1` terms. */
static void integrate_series(const double* c, double* out, size_t stride, int count) {
  for (int k = 1; k <= count; ++k) {
    const double before = k == 1 ? 2.0 * c[0] : c[(size_t)(k - 1) * stride];
    const double after = k + 1 < count ? c[(size_t)(k + 1) * stride] : 0.0;
    out[(size_t)k * stride] = (before - after) / (2.0 * (double)k);
  }
  double at_minus_one = 0.0;
  for (int k = 1; k <= count; ++k) {
    at_minus_one += (k % 2 == 0 ? 1.0 : -1.0) * out[(size_t)k * stride];
  }
  out[0] = -at_minus_one;
}

/* The place and velocity of component i of the interpolant at theta. */
static void interpolant_at(const myb_rkn12* s, size_t i, double theta, double* x, double* v) {
  const size_t n = s->n;
  const double h = s->h_last;
  const double tau = 2.0 * theta - 1.0;
  const double x0 = s->y_old[i] + s->comp_old[i];
  const double v0 = s->y_old[n + i] + s->comp_old[n + i];
  const double alpha = s->fix[i];
  const double beta = s->fix[n + i];
  const double big_x = chebyshev(s->series_x + i, n, DEGREE + 3, tau);
  const double big_v = chebyshev(s->series_v + i, n, DEGREE + 2, tau);
  *x = x0 + (h * theta * v0 + (0.25 * h * h * big_x + theta * theta * (alpha + beta * theta)));
  *v = v0 + (0.5 * h * big_v + theta * (2.0 * alpha + 3.0 * beta * theta) / h);
}

/* The series from the accelerations at the nodes, and the terms that make it end at the end
 * of the step. */
static void make_series(myb_rkn12* s) {
  const size_t n = s->n;
  const double h = s->h_last;
  double c[NODES];
  double v_series[DEGREE + 2];
  double x_series[DEGREE + 3];
  for (size_t i = 0; i < n; ++i) {
    for (int k = 0; k < NODES; ++k) {
      double sum = 0.0;
      for (int j = 0; j < NODES; ++j) {
        sum += s->dct[k][j] * s->node_a[(size_t)j * n + i];
      }
      c[k] = sum;
    }
    integrate_series(c, v_series, 1, NODES);
    integrate_series(v_series, x_series, 1, DEGREE + 2);
    for (int k = 0; k < DEGREE + 2; ++k) {
      s->series_v[(size_t)k * n + i] = v_series[k];
    }
    for (int k = 0; k < DEGREE + 3; ++k) {
      s->series_x[(size_t)k * n + i] = x_series[k];
    }
    /* the mismatch at the end (tau = 1: the sums of the coefficients) */
    double end_x = 0.0;
    double end_v = 0.0;
    for (int k = 0; k < DEGREE + 3; ++k) {
      end_x += x_series[k];
    }
    for (int k = 0; k < DEGREE + 2; ++k) {
      end_v += v_series[k];
    }
    const double x0 = s->y_old[i] + s->comp_old[i];
    const double v0 = s->y_old[n + i] + s->comp_old[n + i];
    const double dx = ((s->y[i] - x0) + s->y_comp[i]) - (h * v0 + 0.25 * h * h * end_x);
    const double dv = ((s->y[n + i] - v0) + s->y_comp[n + i]) - 0.5 * h * end_v;
    s->fix[i] = 3.0 * dx - h * dv;
    s->fix[n + i] = h * dv - 2.0 * dx;
  }
}

/* Builds the interpolant of the last step: 1 if it settled. */
static myb_status build_interpolant(myb_rkn12* s, int* settled) {
  const size_t n = s->n;
  const double h = s->h_last;
  *settled = 0;
  /* the quintic of Hermite through x, h v, h^2 a at both ends */
  for (int j = 0; j < NODES; ++j) {
    const double q = s->node[j];
    const double q2 = q * q;
    const double q3 = q2 * q;
    const double q4 = q3 * q;
    const double q5 = q4 * q;
    const double h00 = 1.0 - 10.0 * q3 + 15.0 * q4 - 6.0 * q5;
    const double h10 = q - 6.0 * q3 + 8.0 * q4 - 3.0 * q5;
    const double h20 = 0.5 * (q2 - 3.0 * q3 + 3.0 * q4 - q5);
    const double h01 = 10.0 * q3 - 15.0 * q4 + 6.0 * q5;
    const double h11 = -4.0 * q3 + 7.0 * q4 - 3.0 * q5;
    const double h21 = 0.5 * (q3 - 2.0 * q4 + q5);
    for (size_t i = 0; i < n; ++i) {
      const double x0 = s->y_old[i] + s->comp_old[i];
      const double x1 = s->y[i] + s->y_comp[i];
      s->node_x[(size_t)j * n + i] =
          x0 * h00 + h * s->y_old[n + i] * h10 + h * h * s->a_old[i] * h20 + x1 * h01 +
          h * s->y[n + i] * h11 + h * h * s->a_now[i] * h21;
    }
  }
  memcpy(s->node_a, s->a_old, n * sizeof(double));
  memcpy(s->node_a + (size_t)DEGREE * n, s->a_now, n * sizeof(double));
  /* Picard's iteration */
  for (int round = 0; round < ROUNDS && !*settled; ++round) { /* bound: ROUNDS */
    for (int j = 1; j < DEGREE; ++j) {
      const double t = s->t_old + (s->t_comp_old + s->node[j] * h);
      const myb_status status =
          evaluate(s, t, s->node_x + (size_t)j * n, s->node_a + (size_t)j * n);
      if (status != MYB_OK) {
        return status;
      }
    }
    make_series(s);
    double worst = 0.0;
    int quiet = 1;
    for (int j = 1; j < DEGREE; ++j) {
      for (size_t i = 0; i < n; ++i) {
        double x = 0.0;
        double v = 0.0;
        interpolant_at(s, i, s->node[j], &x, &v);
        double* old = s->node_x + (size_t)j * n + i;
        const double change = fabs(x - *old);
        const double scale = s->options.atol + s->options.rtol * fabs(x);
        worst = fmax(worst, change / scale);
        quiet = quiet && change <= 4e-16 * fabs(x);
        *old = x;
      }
    }
    *settled = worst < 1e-3 || quiet;
  }
  return MYB_OK;
}

myb_status myb_rkn12_dense(myb_rkn12* integrator, double t, double* y) {
  MYB_CHECK(integrator != NULL && y != NULL);
  myb_rkn12* s = integrator;
  if (!s->have_step) {
    return MYB_E_RANGE;
  }
  const double lo = s->h_last > 0.0 ? s->t_old : s->t;
  const double hi = s->h_last > 0.0 ? s->t : s->t_old;
  if (!(t >= lo && t <= hi)) {
    return MYB_E_RANGE;
  }
  const double h = (t - s->t_old) - s->t_comp_old;
  if (h == 0.0) {
    for (size_t i = 0; i < 2 * s->n; ++i) {
      y[i] = s->y_old[i] + s->comp_old[i];
    }
    return MYB_OK;
  }
  if (!s->options.interpolant) {
    return pair_step(s, s->t_old, s->t_comp_old, s->y_old, s->comp_old, s->a_old, h, y, NULL,
                     NULL);
  }
  if (t == s->t + s->t_comp) {
    for (size_t i = 0; i < 2 * s->n; ++i) {
      y[i] = s->y[i] + s->y_comp[i];
    }
    return MYB_OK;
  }
  if (s->dense_state == 0) {
    int settled = 0;
    const myb_status status = build_interpolant(s, &settled);
    if (status != MYB_OK) {
      return status;
    }
    s->dense_state = settled ? 1 : -1;
  }
  if (s->dense_state < 0) {
    return pair_step(s, s->t_old, s->t_comp_old, s->y_old, s->comp_old, s->a_old, h, y, NULL,
                     NULL);
  }
  const double theta = fmin(fmax(h / s->h_last, 0.0), 1.0);
  for (size_t i = 0; i < s->n; ++i) {
    interpolant_at(s, i, theta, &y[i], &y[s->n + i]);
  }
  return MYB_OK;
}

static myb_status dense_view(void* integrator, double t, double* y) {
  return myb_rkn12_dense(integrator, t, y);
}

myb_status myb_rkn12_locate(myb_rkn12* integrator, myb_event_function g, void* context,
                            int* found, double* t_event, double* y_event) {
  MYB_CHECK(integrator != NULL && g != NULL && found != NULL && t_event != NULL);
  myb_rkn12* s = integrator;
  *found = 0;
  if (!s->have_step) {
    return MYB_E_RANGE;
  }
  const myb_internal_step_view view = {.dense = dense_view,
                                       .integrator = s,
                                       .m = 2 * s->n,
                                       .t0 = s->t_old,
                                       .y0 = s->y_old,
                                       .t1 = s->t,
                                       .y1 = s->y,
                                       .scratch = s->scratch};
  return myb_internal_locate_event(&view, g, context, found, t_event, y_event);
}
