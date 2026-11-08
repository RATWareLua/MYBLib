/* SPDX-License-Identifier: BSD-3-Clause */
/* DOP853: E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I,
 * 2nd ed. (1993), §II.5 (method, error estimate), §II.6 (dense output), §II.4 (step-size
 * control and starting step). The structure of the error norm and of the dense output follows
 * SciPy's DOP853 (scipy/integrate/_ivp/rk.py, BSD-3-Clause), which reproduces Hairer's code. */
#include <math.h> /* sqrt, fabs */
#include <stdckdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "integ/dop853_coefficients.h"
#include "integ/event.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "num/eft.h"

enum { MAX_REJECTIONS = 100, DEFAULT_MAX_STEPS = 1000000000 };

static const double safety = 0.9;
static const double min_factor = 0.2;
static const double max_factor = 10.0;

struct myb_dop853 {
  myb_allocator allocator;
  size_t n;
  myb_dop853_options options;
  myb_ode_function f;
  void* context;
  double t;      /* time ...       */
  double t_comp; /* ... + its compensation */
  double t_old;
  double h_last;    /* last accepted step (signed) */
  double h;         /* next step to try (signed); 0 = choose automatically */
  int started;
  int have_step;
  int dense_ready;
  myb_dop853_stats stats;
  double e3[13];
  double* y;
  double* y_comp;
  double* y_old;
  double* f_now;  /* f(t, y) */
  double* y_new;
  double* comp_new;
  double* increment; /* exact increment of the last accepted step */
  double* scratch;
  double* k;     /* MYB_DOP853_STAGES_EXTENDED x n */
  double* dense; /* MYB_DOP853_INTERPOLATOR x n */
};

static int finite(double x) {
  return x - x == 0.0;
}

static double power(double x, double p) { /* x > 0 */
  return myb_exp(myb_log(x) * p);
}

static double* stage(myb_dop853* s, int i) {
  return s->k + (size_t)i * s->n;
}

myb_status myb_dop853_create(size_t n, const myb_dop853_options* options,
                             const myb_allocator* allocator, myb_dop853** integrator) {
  MYB_CHECK(options != NULL && integrator != NULL);
  const myb_dop853_options o = *options;
  if (n == 0 || !(o.rtol > 0.0) || !finite(o.rtol) || !(o.atol > 0.0) || !finite(o.atol) ||
      !(o.h_initial >= 0.0) || !(o.h_max >= 0.0) || !(o.h_min >= 0.0) || o.max_steps < 0) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_dop853), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_dop853* s = block;
  *s = (myb_dop853){.allocator = a, .n = n, .options = o};
  if (s->options.max_steps == 0) {
    s->options.max_steps = DEFAULT_MAX_STEPS;
  }
  /* 8 vectors + 16 stages + 7 interpolation rows. */
  const size_t vectors = 8 + MYB_DOP853_STAGES_EXTENDED + MYB_DOP853_INTERPOLATOR;
  size_t count = 0;
  if (ckd_mul(&count, vectors, n)) {
    status = MYB_E_OVERFLOW;
    goto fail;
  }
  status = myb_internal_alloc_array(&a, count, sizeof(double), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  double* p = block;
  memset(p, 0, count * sizeof(double));
  double** slots[8] = {&s->y,        &s->y_comp,   &s->y_old,     &s->f_now,
                       &s->y_new,    &s->comp_new, &s->increment, &s->scratch};
  for (int i = 0; i < 8; ++i) {
    *slots[i] = p + (size_t)i * n;
  }
  s->k = p + 8 * n;
  s->dense = s->k + (size_t)MYB_DOP853_STAGES_EXTENDED * n;
  /* E3 = B - (0.2440..., 0.7338..., 0.02205...) at stages 0, 8, 11; E3[12] = 0. */
  for (int j = 0; j < 12; ++j) {
    s->e3[j] = myb_dop853_a[12][j];
  }
  s->e3[12] = 0.0;
  for (int i = 0; i < 3; ++i) {
    const int j = (int)myb_dop853_e3_sub[i][0];
    s->e3[j] -= myb_dop853_e3_sub[i][1];
  }
  *integrator = s;
  return MYB_OK;

fail:
  myb_internal_release(&a, s);
  return status;
}

void myb_dop853_destroy(myb_dop853* integrator) {
  if (integrator == NULL) {
    return;
  }
  const myb_allocator a = integrator->allocator;
  myb_internal_release(&a, integrator->y); /* the single block of all vectors */
  myb_internal_release(&a, integrator);
}

static myb_status evaluate(myb_dop853* s, double t, const double* y, double* dydt) {
  s->stats.evaluations++;
  return s->f(s->context, t, y, dydt);
}

myb_status myb_dop853_start(myb_dop853* integrator, myb_ode_function f, void* context, double t0,
                            const double* y0) {
  MYB_CHECK(integrator != NULL && f != NULL && y0 != NULL);
  myb_dop853* s = integrator;
  if (!finite(t0)) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < s->n; ++i) {
    if (!finite(y0[i])) {
      return MYB_E_DOMAIN;
    }
    s->y[i] = y0[i];
    s->y_comp[i] = 0.0;
  }
  s->f = f;
  s->context = context;
  s->t = t0;
  s->t_comp = 0.0;
  s->h = s->options.h_initial;
  s->have_step = 0;
  s->dense_ready = 0;
  s->stats = (myb_dop853_stats){};
  const myb_status status = evaluate(s, t0, s->y, s->f_now);
  s->started = status == MYB_OK;
  return status;
}

static double rms_scaled(const myb_dop853* s, const double* v, const double* scale) {
  double sum = 0.0;
  for (size_t i = 0; i < s->n; ++i) {
    const double q = v[i] / scale[i];
    sum += q * q;
  }
  return sqrt(sum / (double)s->n);
}

/* Hairer's starting step (Solving ODE I, §II.4). */
static myb_status initial_step(myb_dop853* s, double direction, double* h_out) {
  double* scale = s->scratch;
  for (size_t i = 0; i < s->n; ++i) {
    scale[i] = s->options.atol + s->options.rtol * fabs(s->y[i]);
  }
  const double d0 = rms_scaled(s, s->y, scale);
  const double d1 = rms_scaled(s, s->f_now, scale);
  const double h0 = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1;
  double* y1 = s->y_new;
  for (size_t i = 0; i < s->n; ++i) {
    y1[i] = s->y[i] + direction * h0 * s->f_now[i];
  }
  double* f1 = stage(s, 1);
  const myb_status status = evaluate(s, s->t + direction * h0, y1, f1);
  if (status != MYB_OK) {
    return status;
  }
  for (size_t i = 0; i < s->n; ++i) {
    f1[i] -= s->f_now[i];
  }
  const double d2 = rms_scaled(s, f1, scale) / h0;
  const double larger = d1 > d2 ? d1 : d2;
  const double h1 = larger <= 1e-15 ? (h0 * 1e-3 > 1e-6 ? h0 * 1e-3 : 1e-6)
                                    : power(0.01 / larger, 1.0 / 8.0);
  *h_out = direction * (100.0 * h0 < h1 ? 100.0 * h0 : h1);
  return MYB_OK;
}

/* One attempt with step h: stages, compensated new state, error norm. */
static myb_status attempt(myb_dop853* s, double h, double* error_norm) {
  const size_t n = s->n;
  memcpy(stage(s, 0), s->f_now, n * sizeof(double));
  for (int st = 1; st <= MYB_DOP853_STAGES; ++st) {
    const double* a = myb_dop853_a[st];
    double* target = st < MYB_DOP853_STAGES ? s->scratch : s->increment;
    for (size_t i = 0; i < n; ++i) {
      double sum = 0.0;
      for (int j = 0; j < st; ++j) {
        sum += a[j] * stage(s, j)[i];
      }
      target[i] = h * sum;
    }
    if (st < MYB_DOP853_STAGES) {
      for (size_t i = 0; i < n; ++i) {
        s->scratch[i] = s->y[i] + s->scratch[i];
      }
      const myb_status status = evaluate(s, s->t + myb_dop853_c[st] * h, s->scratch, stage(s, st));
      if (status != MYB_OK) {
        return status;
      }
    }
  }
  /* y_new = y + increment with compensation: (y_new + comp_new) carries the exact sum. */
  for (size_t i = 0; i < n; ++i) {
    double err;
    s->y_new[i] = myb_eft_two_sum(s->y[i], s->increment[i] + s->y_comp[i], &err);
    s->comp_new[i] = err;
  }
  const myb_status status = evaluate(s, s->t + h, s->y_new, stage(s, MYB_DOP853_STAGES));
  if (status != MYB_OK) {
    return status;
  }
  double e5_sum = 0.0;
  double e3_sum = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double e5 = 0.0;
    double e3 = 0.0;
    for (int j = 0; j <= MYB_DOP853_STAGES; ++j) {
      e5 += myb_dop853_e5[j] * stage(s, j)[i];
      e3 += s->e3[j] * stage(s, j)[i];
    }
    const double y_abs = fabs(s->y[i]) > fabs(s->y_new[i]) ? fabs(s->y[i]) : fabs(s->y_new[i]);
    const double scale = s->options.atol + s->options.rtol * y_abs;
    e5_sum += (e5 / scale) * (e5 / scale);
    e3_sum += (e3 / scale) * (e3 / scale);
  }
  if (e5_sum == 0.0 && e3_sum == 0.0) {
    *error_norm = 0.0;
  } else {
    *error_norm = fabs(h) * e5_sum / sqrt((e5_sum + 0.01 * e3_sum) * (double)n);
  }
  return finite(*error_norm) ? MYB_OK : MYB_E_STEP_MIN;
}

static double clamp_step(const myb_dop853* s, double h) {
  if (s->options.h_max > 0.0 && fabs(h) > s->options.h_max) {
    return h > 0.0 ? s->options.h_max : -s->options.h_max;
  }
  return h;
}

myb_status myb_dop853_step(myb_dop853* integrator, double t_limit) {
  MYB_CHECK(integrator != NULL && integrator->started);
  myb_dop853* s = integrator;
  if (!finite(t_limit)) {
    return MYB_E_DOMAIN;
  }
  /* Remaining interval, using the compensated time. */
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
    /* As in Hairer's code: a step within 1% of the remaining interval is stretched to its end,
     * so no sliver of a step (below the resolution of t) is ever left over. */
    const int last = 1.01 * fabs(h) >= fabs(remaining);
    const double h_planned = h;
    if (last) {
      h = remaining;
    }
    if (!last && s->options.h_min > 0.0 && fabs(h) < s->options.h_min) {
      return MYB_E_STEP_MIN;
    }
    if (!last && s->t + h == s->t) {
      return MYB_E_STEP_MIN; /* the step collapsed below the resolution of t */
    }
    double norm = 0.0;
    const myb_status status = attempt(s, h, &norm);
    if (status != MYB_OK) {
      return status;
    }
    if (norm < 1.0) {
      double factor = norm == 0.0 ? max_factor : safety * power(norm, -1.0 / 8.0);
      factor = factor < max_factor ? factor : max_factor;
      if (previous_rejected && factor > 1.0) {
        factor = 1.0;
      }
      s->t_old = s->t;
      memcpy(s->y_old, s->y, s->n * sizeof(double));
      double err;
      const double t_new = myb_eft_two_sum(s->t, h + s->t_comp, &err);
      s->t = last ? t_limit : t_new;
      s->t_comp = last ? 0.0 : err;
      memcpy(s->y, s->y_new, s->n * sizeof(double));
      memcpy(s->y_comp, s->comp_new, s->n * sizeof(double));
      memcpy(s->f_now, stage(s, MYB_DOP853_STAGES), s->n * sizeof(double));
      s->h_last = h;
      /* A step shortened to land on t_limit says little about the next one: keep the plan. */
      const double h_next = h * factor;
      s->h = clamp_step(s, last && fabs(h_planned) > fabs(h_next) ? h_planned : h_next);
      s->have_step = 1;
      s->dense_ready = 0;
      s->stats.accepted++;
      return MYB_OK;
    }
    s->stats.rejected++;
    double factor = safety * power(norm, -1.0 / 8.0);
    factor = factor > min_factor ? factor : min_factor;
    h *= factor;
    previous_rejected = 1;
  }
  return MYB_E_NOCONV;
}

myb_status myb_dop853_integrate(myb_dop853* integrator, double t_end) {
  MYB_CHECK(integrator != NULL && integrator->started);
  /* bound: max_steps accepted steps */
  for (int64_t i = 0; i < integrator->options.max_steps; ++i) {
    if (integrator->t == t_end && integrator->t_comp == 0.0) {
      return MYB_OK;
    }
    const myb_status status = myb_dop853_step(integrator, t_end);
    if (status != MYB_OK) {
      return status;
    }
  }
  return integrator->t == t_end ? MYB_OK : MYB_E_NOCONV;
}

double myb_dop853_time(const myb_dop853* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->t + integrator->t_comp;
}

const double* myb_dop853_state(const myb_dop853* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->y;
}

myb_dop853_stats myb_dop853_statistics(const myb_dop853* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->stats;
}

/* Extra stages and the seven interpolation rows of the last step. */
static myb_status prepare_dense(myb_dop853* s) {
  const size_t n = s->n;
  const double h = s->h_last;
  for (int st = MYB_DOP853_STAGES + 1; st < MYB_DOP853_STAGES_EXTENDED; ++st) {
    const double* a = myb_dop853_a[st];
    for (size_t i = 0; i < n; ++i) {
      double sum = 0.0;
      for (int j = 0; j < st; ++j) {
        sum += a[j] * stage(s, j)[i];
      }
      s->scratch[i] = s->y_old[i] + h * sum;
    }
    const myb_status status =
        evaluate(s, s->t_old + myb_dop853_c[st] * h, s->scratch, stage(s, st));
    if (status != MYB_OK) {
      return status;
    }
  }
  const double* f_old = stage(s, 0);
  const double* f_new = stage(s, MYB_DOP853_STAGES);
  double* row[MYB_DOP853_INTERPOLATOR];
  for (int r = 0; r < MYB_DOP853_INTERPOLATOR; ++r) {
    row[r] = s->dense + (size_t)r * n;
  }
  for (size_t i = 0; i < n; ++i) {
    const double dy = s->increment[i];
    row[0][i] = dy;
    row[1][i] = h * f_old[i] - dy;
    row[2][i] = 2.0 * dy - h * (f_new[i] + f_old[i]);
    for (int r = 0; r < 4; ++r) {
      double sum = 0.0;
      for (int j = 0; j < MYB_DOP853_STAGES_EXTENDED; ++j) {
        sum += myb_dop853_d[r][j] * stage(s, j)[i];
      }
      row[3 + r][i] = h * sum;
    }
  }
  s->dense_ready = 1;
  return MYB_OK;
}

myb_status myb_dop853_dense(myb_dop853* integrator, double t, double* y) {
  MYB_CHECK(integrator != NULL && y != NULL);
  myb_dop853* s = integrator;
  if (!s->have_step) {
    return MYB_E_RANGE;
  }
  const double lo = s->h_last > 0.0 ? s->t_old : s->t;
  const double hi = s->h_last > 0.0 ? s->t : s->t_old;
  if (!(t >= lo && t <= hi)) {
    return MYB_E_RANGE;
  }
  if (!s->dense_ready) {
    const myb_status status = prepare_dense(s);
    if (status != MYB_OK) {
      return status;
    }
  }
  const double x = (t - s->t_old) / s->h_last;
  for (size_t i = 0; i < s->n; ++i) {
    double v = 0.0;
    /* Hairer's nested form: rows 6, 5, ..., 0 multiplied alternately by x and 1 - x. */
    for (int k = 0; k < MYB_DOP853_INTERPOLATOR; ++k) {
      v += s->dense[(size_t)(MYB_DOP853_INTERPOLATOR - 1 - k) * s->n + i];
      v *= (k % 2 == 0) ? x : 1.0 - x;
    }
    y[i] = s->y_old[i] + v;
  }
  return MYB_OK;
}

static myb_status dense_view(void* integrator, double t, double* y) {
  return myb_dop853_dense(integrator, t, y);
}

myb_status myb_dop853_locate(myb_dop853* integrator, myb_event_function g, void* context,
                             int* found, double* t_event, double* y_event) {
  MYB_CHECK(integrator != NULL && g != NULL && found != NULL && t_event != NULL);
  myb_dop853* s = integrator;
  *found = 0;
  if (!s->have_step) {
    return MYB_E_RANGE;
  }
  const myb_internal_step_view view = {.dense = dense_view,
                                       .integrator = s,
                                       .m = s->n,
                                       .t0 = s->t_old,
                                       .y0 = s->y_old,
                                       .t1 = s->t,
                                       .y1 = s->y,
                                       .scratch = s->scratch};
  return myb_internal_locate_event(&view, g, context, found, t_event, y_event);
}
