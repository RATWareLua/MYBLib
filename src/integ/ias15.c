/* SPDX-License-Identifier: BSD-3-Clause */
/* IAS15: H. Rein, D. S. Spiegel, MNRAS 446, 1424 (2015) — Gauss-Radau integration of
 * x'' = a(t, x, v) in Everhart's form (E. Everhart, in Dynamics of Comets, 1985), with the
 * step-size criterion of D. Pham, H. Rein, D. S. Spiegel, Open J. Astrophys. 7 (2024).
 * Implemented from the papers; the constants are derived in integ/ias15_coefficients.h.
 *
 * Within a step t = t0 + h dt, h in [0, 1]:
 *   a(h) = a0 + b0 h + b1 h^2 + ... + b6 h^7,
 *   v(h) = v0 + dt h (a0 + b0 h / 2 + b1 h^2 / 3 + ... + b6 h^7 / 8),
 *   x(h) = x0 + dt h (v0 + dt h (a0 / 2 + b0 h / 6 + b1 h^2 / 12 + ... + b6 h^7 / 72)).
 * A sweep evaluates a at the nodes h_1..h_7 with the current b, updates the divided differences
 * g and with them b. Prediction for the next step: the polynomial of the last step re-expanded
 * around its end (q = new step / old step),
 *   e_m = q^(m+1) sum_{k>=m} binom(k+1, m+1) b_k,
 * plus the last step's correction (b - e) scaled by q^(m+1). */
#include <math.h> /* sqrt, fabs */
#include <stdckdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "integ/event.h"
#include "integ/ias15_internal.h"
#include "integ/ias15_coefficients.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "num/eft.h"

enum {
  ORDER = 7,           /* b_0..b_6 */
  MAX_SWEEPS = 12,     /* predictor-corrector sweeps per attempt */
  MAX_REJECTIONS = 100,
  DEFAULT_MAX_STEPS = 1000000000,
};

static const double default_epsilon = 1e-9;
static const double safety = 0.25;          /* step change bound, both ways */
static const double converged_change = 1e-16; /* |delta b6| / |a| at which the sweeps stop */

/* 1 / ((k + 2)(k + 3)) and 1 / (k + 2): integration of h^(k+1) twice and once. */
static const double x_weight[ORDER] = {1.0 / 6.0,  1.0 / 12.0, 1.0 / 20.0, 1.0 / 30.0,
                                       1.0 / 42.0, 1.0 / 56.0, 1.0 / 72.0};
static const double v_weight[ORDER] = {1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0,
                                       1.0 / 6.0, 1.0 / 7.0, 1.0 / 8.0};

struct myb_ias15 {
  myb_allocator allocator;
  size_t n;
  myb_ias15_options options;
  double eta; /* (7! epsilon)^(1/7) */
  myb_accel_function f;
  void* context;
  double t;
  double t_comp;
  double t_old;
  double t_comp_old;
  double h_last; /* last accepted step (signed) */
  double h;      /* next step (signed); 0 = automatic */
  int started;
  int have_step;
  int a0_valid;         /* a0 holds a(t, y) */
  int prediction;       /* e holds a genuine prediction for the current step */
  int last_prediction;  /* e_last held one for the last step */
  myb_ias15_stats stats;
  double* y;        /* 2n: x, v */
  double* y_comp;   /* 2n */
  double* y_old;    /* 2n: start of the last accepted step */
  double* comp_old; /* 2n */
  double* y_new;    /* 2n */
  double* comp_new; /* 2n */
  double* y_sub;    /* 2n: state at a node */
  double* scratch;  /* 2n */
  double* a0;       /* n: a(t, y) */
  double* a_old;    /* n: a at the start of the last accepted step */
  double* a_sub;    /* n */
  double* g;        /* 7n, component-major: g[7 i + k] */
  double* b;        /* 7n */
  double* e;        /* 7n: prediction the current step started from */
  double* b_last;   /* 7n: converged b of the last accepted step */
  double* e_last;   /* 7n */
  /* a step withdrawn (myb_internal_ias15_retake): its polynomial predicts the next step */
  int hint;
  double hint_h;  /* the length of the step withdrawn (signed) */
  double hint_s0; /* where the next step starts in it, as a fraction */
  double* b_hint; /* 7n */
  /* the groups of the options: where each ends (one group ending at n without them) */
  size_t group_count;
  size_t* group_end;
};

static int finite(double x) {
  return x - x == 0.0;
}

static double larger_of(double a, double b) {
  return a > b ? a : b;
}

myb_status myb_ias15_create(size_t n, const myb_ias15_options* options,
                            const myb_allocator* allocator, myb_ias15** integrator) {
  MYB_CHECK(options != NULL && integrator != NULL);
  myb_ias15_options o = *options;
  if (o.epsilon == 0.0) {
    o.epsilon = default_epsilon;
  }
  if (n == 0 || !(o.epsilon > 0.0) || !finite(o.epsilon) || !(o.h_initial >= 0.0) ||
      !finite(o.h_initial) || !(o.h_max >= 0.0) || !(o.h_min >= 0.0) || o.max_steps < 0 ||
      (o.block != 0 && n % o.block != 0)) {
    return MYB_E_DOMAIN;
  }
  if (o.max_steps == 0) {
    o.max_steps = DEFAULT_MAX_STEPS;
  }
  if (o.block == 0) {
    o.block = n;
  }
  if ((o.group_count == 0) != (o.groups == NULL)) {
    return MYB_E_DOMAIN;
  }
  size_t covered = 0;
  for (size_t k = 0; k < o.group_count; ++k) {
    if (o.groups[k] == 0 || o.groups[k] % o.block != 0 || o.groups[k] > n - covered) {
      return MYB_E_DOMAIN;
    }
    covered += o.groups[k];
  }
  if (o.group_count != 0 && covered != n) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_ias15), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_ias15* s = block;
  *s = (myb_ias15){.allocator = a, .n = n, .options = o};
  s->options.groups = NULL; /* the caller's array is not kept */
  s->eta = myb_exp(myb_log(5040.0 * o.epsilon) / 7.0);
  s->group_count = o.group_count != 0 ? o.group_count : 1;
  status = myb_internal_alloc_array(&a, s->group_count, sizeof(size_t), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  s->group_end = block;
  s->group_end[0] = n;
  covered = 0;
  for (size_t k = 0; k < o.group_count; ++k) {
    covered += o.groups[k];
    s->group_end[k] = covered;
  }
  /* 8 vectors of 2n, 3 of n, 6 of 7n: 61 n values. */
  size_t count = 0;
  if (ckd_mul(&count, (size_t)(8 * 2 + 3 + 6 * ORDER), n)) {
    status = MYB_E_OVERFLOW;
    goto fail;
  }
  status = myb_internal_alloc_array(&a, count, sizeof(double), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  double* p = block;
  memset(p, 0, count * sizeof(double));
  double** wide[8] = {&s->y,     &s->y_comp,   &s->y_old, &s->comp_old,
                      &s->y_new, &s->comp_new, &s->y_sub, &s->scratch};
  for (int i = 0; i < 8; ++i) {
    *wide[i] = p;
    p += 2 * n;
  }
  double** narrow[3] = {&s->a0, &s->a_old, &s->a_sub};
  for (int i = 0; i < 3; ++i) {
    *narrow[i] = p;
    p += n;
  }
  double** series[6] = {&s->g, &s->b, &s->e, &s->b_last, &s->e_last, &s->b_hint};
  for (int i = 0; i < 6; ++i) {
    *series[i] = p;
    p += (size_t)ORDER * n;
  }
  *integrator = s;
  return MYB_OK;

fail:
  myb_internal_release(&a, s->group_end);
  myb_internal_release(&a, s);
  return status;
}

void myb_ias15_destroy(myb_ias15* integrator) {
  if (integrator == NULL) {
    return;
  }
  const myb_allocator a = integrator->allocator;
  myb_internal_release(&a, integrator->y); /* the single block of all vectors */
  myb_internal_release(&a, integrator->group_end);
  myb_internal_release(&a, integrator);
}

static myb_status evaluate(myb_ias15* s, double t, const double* y, double* a) {
  s->stats.evaluations++;
  return s->f(s->context, t, y, y + s->n, a);
}

/* Norm of a block of components: sqrt of the sum of squares. */
static double block_norm(const double* v, size_t count) {
  double sum = 0.0;
  for (size_t i = 0; i < count; ++i) {
    sum += v[i] * v[i];
  }
  return sqrt(sum);
}

static double automatic_step(const myb_ias15* s) {
  const size_t n = s->n;
  const size_t m = s->options.block;
  double tau = (double)INFINITY;
  for (size_t first = 0; first < n; first += m) {
    const double x = block_norm(s->y + first, m);
    const double v = block_norm(s->y + n + first, m);
    const double a = block_norm(s->a0 + first, m);
    if (!(a > 0.0) || !finite(a)) {
      continue;
    }
    const double from_x = sqrt(x / a);
    const double from_v = v / a;
    if (from_x > 0.0 && from_x < tau) {
      tau = from_x;
    }
    if (from_v > 0.0 && from_v < tau) {
      tau = from_v;
    }
  }
  return s->eta * tau; /* infinity: no information, the remaining interval is taken */
}

myb_status myb_ias15_start(myb_ias15* integrator, myb_accel_function a, void* context, double t0,
                           const double* y0) {
  MYB_CHECK(integrator != NULL && a != NULL && y0 != NULL);
  myb_ias15* s = integrator;
  if (!finite(t0)) {
    return MYB_E_DOMAIN;
  }
  const size_t n = s->n;
  for (size_t i = 0; i < 2 * n; ++i) {
    if (!finite(y0[i])) {
      return MYB_E_DOMAIN;
    }
    s->y[i] = y0[i];
    s->y_comp[i] = 0.0;
  }
  memset(s->b, 0, (size_t)ORDER * n * sizeof(double));
  memset(s->e, 0, (size_t)ORDER * n * sizeof(double));
  s->f = a;
  s->context = context;
  s->t = t0;
  s->t_comp = 0.0;
  s->hint = 0;
  s->have_step = 0;
  s->prediction = 0;
  s->last_prediction = 0;
  s->stats = (myb_ias15_stats){};
  const myb_status status = evaluate(s, t0, s->y, s->a0);
  s->a0_valid = status == MYB_OK;
  s->started = status == MYB_OK;
  if (status == MYB_OK) {
    s->h = s->options.h_initial > 0.0 ? s->options.h_initial : automatic_step(s);
  }
  return status;
}

/* Position and velocity polynomials at h for component i (Horner, highest term first). */
static void node_state(const double* b, double a0, double h, double* p, double* q) {
  double px = b[ORDER - 1] * x_weight[ORDER - 1];
  double pv = b[ORDER - 1] * v_weight[ORDER - 1];
  for (int k = ORDER - 2; k >= 0; --k) {
    px = px * h + b[k] * x_weight[k];
    pv = pv * h + b[k] * v_weight[k];
  }
  *p = px * h + 0.5 * a0;
  *q = pv * h + a0;
}

/* State at fraction h of a step dt from (y, comp) with accelerations a0 and series b. The
 * compensation enters before the rounding of the sum, as in the update of the state. */
static void state_at(size_t n, const double* y, const double* comp, const double* a0,
                     const double* b, double dt, double h, double* out, double* out_comp) {
  const double dth = dt * h;
  for (size_t i = 0; i < n; ++i) {
    double p = 0.0;
    double q = 0.0;
    node_state(b + (size_t)ORDER * i, a0[i], h, &p, &q);
    const double dx = comp[i] + dth * (y[n + i] + dth * p);
    const double dv = comp[n + i] + dth * q;
    if (out_comp != NULL) {
      out[i] = myb_eft_two_sum(y[i], dx, &out_comp[i]);
      out[n + i] = myb_eft_two_sum(y[n + i], dv, &out_comp[n + i]);
    } else {
      out[i] = y[i] + dx;
      out[n + i] = y[n + i] + dv;
    }
  }
}

/* g from b: g_j = sum_k D[j][k] b_k. */
static void g_from_b(myb_ias15* s) {
  for (size_t i = 0; i < s->n; ++i) {
    const double* b = s->b + (size_t)ORDER * i;
    double* g = s->g + (size_t)ORDER * i;
    for (int j = 0; j < ORDER; ++j) {
      double sum = 0.0;
      for (int k = ORDER - 1; k >= j; --k) {
        sum += myb_ias15_d[j][k] * b[k];
      }
      g[j] = sum;
    }
  }
}

/* Predictor-corrector sweeps for the step dt from the current b. */
static myb_status sweeps(myb_ias15* s, double dt) {
  const size_t n = s->n;
  g_from_b(s);
  double change_last = (double)INFINITY;
  /* bound: MAX_SWEEPS */
  for (int sweep = 0; sweep < MAX_SWEEPS; ++sweep) {
    s->stats.iterations++;
    /* the relative change of b6, the worst over the groups; NaN once a group gives one */
    double relative = 0.0;
    for (int node = 1; node <= ORDER; ++node) {
      const double h = myb_ias15_h[node];
      state_at(n, s->y, s->y_comp, s->a0, s->b, dt, h, s->y_sub, NULL);
      const myb_status status = evaluate(s, s->t + (s->t_comp + h * dt), s->y_sub, s->a_sub);
      if (status != MYB_OK) {
        return status;
      }
      const double* r = myb_ias15_r[node];
      const int j = node - 1; /* the divided difference this node completes */
      double change = 0.0;
      double scale = 0.0;
      size_t group = 0;
      for (size_t i = 0; i < n; ++i) {
        double* g = s->g + (size_t)ORDER * i;
        double* b = s->b + (size_t)ORDER * i;
        double d = (s->a_sub[i] - s->a0[i]) * r[0];
        for (int k = 1; k < node; ++k) {
          d = (d - g[k - 1]) * r[k];
        }
        const double delta = d - g[j];
        g[j] = d;
        for (int k = 0; k <= j; ++k) {
          b[k] += myb_ias15_c[k][j] * delta;
        }
        if (node == ORDER) {
          change = larger_of(change, fabs(delta));
          scale = larger_of(scale, fabs(s->a_sub[i]));
          if (i + 1 == s->group_end[group]) {
            /* Relative change of b6 (C[6][6] = 1, so delta b6 = delta g6) of the group. */
            const double of_group =
                scale > 0.0 ? change / scale : (change > 0.0 ? change : 0.0);
            if (of_group != of_group || relative != relative) {
              relative = (double)NAN;
            } else {
              relative = larger_of(relative, of_group);
            }
            change = 0.0;
            scale = 0.0;
            ++group;
          }
        }
      }
    }
    /* NaN compares false. */
    if (relative != relative) {
      return MYB_OK; /* NaN: the series diverged; the step is rejected by its proposal */
    }
    if (relative < converged_change || (sweep >= 2 && !(relative < change_last))) {
      return MYB_OK; /* converged, or the corrections stopped shrinking (round-off level) */
    }
    change_last = relative;
  }
  s->stats.unconverged++;
  return MYB_OK;
}

/* Proposed next step after an attempt with dt; NaN if the series is not finite. The norm of
 * each derivative over the whole state is the maximum over blocks of the Euclidean norm of the
 * block (Pham, Rein & Spiegel 2024, eq. 7), taken separately for a, a' and a'' before eq. 16:
 * a block whose acceleration passes through zero (the middle body of a symmetric collinear
 * configuration) does not shrink the step, as a minimum of per-block timescales would. */
static double proposed_of_group(const myb_ias15* s, double dt, size_t begin, size_t end) {
  const size_t m = s->options.block;
  double y2 = 0.0; /* max over blocks of |a|^2, |dt a'|^2, |dt^2 a''|^2 */
  double y3 = 0.0;
  double y4 = 0.0;
  for (size_t first = begin; first < end; first += m) {
    double b2 = 0.0;
    double b3 = 0.0;
    double b4 = 0.0;
    for (size_t i = first; i < first + m; ++i) {
      const double* b = s->b + (size_t)ORDER * i;
      /* a, dt a', dt^2 a'' at h = 1 */
      const double a = s->a0[i] + (b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6]);
      const double j =
          b[0] + 2.0 * b[1] + 3.0 * b[2] + 4.0 * b[3] + 5.0 * b[4] + 6.0 * b[5] + 7.0 * b[6];
      const double k =
          2.0 * b[1] + 6.0 * b[2] + 12.0 * b[3] + 20.0 * b[4] + 30.0 * b[5] + 42.0 * b[6];
      b2 += a * a;
      b3 += j * j;
      b4 += k * k;
    }
    if (!finite(b2) || !finite(b3) || !finite(b4)) {
      return (double)NAN;
    }
    y2 = b2 > y2 ? b2 : y2;
    y3 = b3 > y3 ? b3 : y3;
    y4 = b4 > y4 ? b4 : y4;
  }
  if (!(y2 > 0.0)) {
    return (double)INFINITY; /* no acceleration anywhere: no timescale */
  }
  /* tau^2 / dt^2 = 2 / (y3 / y2 + sqrt(y4 / y2)), free of overflow in y4 y2 */
  const double denominator = y3 / y2 + sqrt(y4 / y2);
  if (!(denominator > 0.0)) {
    return (double)INFINITY;
  }
  return s->eta * sqrt(2.0 / denominator) * dt;
}

/* The shortest proposal of the groups (in magnitude; dt carries the direction). */
static double proposed_step(const myb_ias15* s, double dt) {
  double proposed = (double)INFINITY;
  size_t begin = 0;
  for (size_t k = 0; k < s->group_count; ++k) {
    const double of_group = proposed_of_group(s, dt, begin, s->group_end[k]);
    if (of_group != of_group) {
      return of_group;
    }
    if (fabs(of_group) < fabs(proposed)) {
      proposed = of_group;
    }
    begin = s->group_end[k];
  }
  return proposed;
}

static double power_int(double q, int k) {
  double r = 1.0;
  for (int i = 0; i < k; ++i) {
    r *= q;
  }
  return r;
}

/* binom(k + 1, m + 1) for m <= k */
static const double binomial[ORDER][ORDER] = {
    {1, 2, 3, 4, 5, 6, 7},   {0, 1, 3, 6, 10, 15, 21}, {0, 0, 1, 4, 10, 20, 35},
    {0, 0, 0, 1, 5, 15, 35}, {0, 0, 0, 0, 1, 6, 21},   {0, 0, 0, 0, 0, 1, 7},
    {0, 0, 0, 0, 0, 0, 1},
};

/* b and e for a step of q times the step withdrawn, starting at the fraction s0 of it: the
 * polynomial sum b_k s^(k+1) of that step in the variable s = s0 + q u. */
static void predict_from_hint(myb_ias15* s, double q) {
  const size_t n = s->n;
  double qp[ORDER];
  double sp[ORDER];
  for (int m = 0; m < ORDER; ++m) {
    qp[m] = power_int(q, m + 1);
    sp[m] = power_int(s->hint_s0, m);
  }
  for (size_t i = 0; i < n; ++i) {
    const double* bh = s->b_hint + (size_t)ORDER * i;
    double* e = s->e + (size_t)ORDER * i;
    double* b = s->b + (size_t)ORDER * i;
    for (int m = 0; m < ORDER; ++m) {
      double sum = 0.0;
      for (int k = ORDER - 1; k >= m; --k) {
        sum += binomial[m][k] * sp[k - m] * bh[k];
      }
      e[m] = qp[m] * sum;
      b[m] = e[m];
    }
  }
  s->prediction = 1;
}

/* b and e for a step of q times the last accepted one, from b_last and e_last. */
static void predict(myb_ias15* s, double q) {
  const size_t n = s->n;
  if (!s->have_step) {
    memset(s->b, 0, (size_t)ORDER * n * sizeof(double));
    memset(s->e, 0, (size_t)ORDER * n * sizeof(double));
    s->prediction = 0;
    return;
  }
  double qp[ORDER];
  for (int m = 0; m < ORDER; ++m) {
    qp[m] = power_int(q, m + 1);
  }
  for (size_t i = 0; i < n; ++i) {
    const double* bl = s->b_last + (size_t)ORDER * i;
    const double* el = s->e_last + (size_t)ORDER * i;
    double* e = s->e + (size_t)ORDER * i;
    double* b = s->b + (size_t)ORDER * i;
    for (int m = 0; m < ORDER; ++m) {
      double sum = 0.0;
      for (int k = ORDER - 1; k >= m; --k) {
        sum += binomial[m][k] * bl[k];
      }
      e[m] = qp[m] * sum;
      b[m] = s->last_prediction ? e[m] + qp[m] * (bl[m] - el[m]) : e[m];
    }
  }
  s->prediction = 1;
}

static double clamp_step(const myb_ias15* s, double h) {
  if (s->options.h_max > 0.0 && fabs(h) > s->options.h_max) {
    return h > 0.0 ? s->options.h_max : -s->options.h_max;
  }
  return h;
}

static void accept(myb_ias15* s, double dt, int last, double t_limit) {
  const size_t n = s->n;
  s->t_old = s->t;
  s->t_comp_old = s->t_comp;
  memcpy(s->y_old, s->y, 2 * n * sizeof(double));
  memcpy(s->comp_old, s->y_comp, 2 * n * sizeof(double));
  memcpy(s->a_old, s->a0, n * sizeof(double));
  memcpy(s->y, s->y_new, 2 * n * sizeof(double));
  memcpy(s->y_comp, s->comp_new, 2 * n * sizeof(double));
  double err = 0.0;
  const double t_new = myb_eft_two_sum(s->t, dt + s->t_comp, &err);
  s->t = last ? t_limit : t_new;
  s->t_comp = last ? 0.0 : err;
  memcpy(s->b_last, s->b, (size_t)ORDER * n * sizeof(double));
  memcpy(s->e_last, s->e, (size_t)ORDER * n * sizeof(double));
  s->last_prediction = s->prediction;
  s->h_last = dt;
  s->have_step = 1;
  s->a0_valid = 0;
  s->hint = 0;
  s->stats.accepted++;
}

myb_status myb_ias15_step(myb_ias15* integrator, double t_limit) {
  MYB_CHECK(integrator != NULL && integrator->started);
  myb_ias15* s = integrator;
  if (!finite(t_limit)) {
    return MYB_E_DOMAIN;
  }
  const double remaining = (t_limit - s->t) - s->t_comp;
  if (remaining == 0.0) {
    return MYB_OK;
  }
  if (!s->a0_valid) {
    const myb_status status = evaluate(s, s->t + s->t_comp, s->y, s->a0);
    if (status != MYB_OK) {
      return status;
    }
    s->a0_valid = 1;
  }
  const double direction = remaining > 0.0 ? 1.0 : -1.0;
  double dt = clamp_step(s, direction * fabs(s->h));
  /* bound: MAX_REJECTIONS; every rejection shortens the step at least 4 times */
  for (int tries = 0; tries < MAX_REJECTIONS; ++tries) {
    const int last = !(1.01 * fabs(dt) < fabs(remaining)); /* also for an infinite dt */
    const double dt_planned = dt;
    if (last) {
      dt = remaining;
    }
    if (!last && s->options.h_min > 0.0 && fabs(dt) < s->options.h_min) {
      return MYB_E_STEP_MIN;
    }
    if (!last && s->t + dt == s->t) {
      return MYB_E_STEP_MIN;
    }
    if (s->hint) {
      predict_from_hint(s, dt / s->hint_h);
    } else {
      predict(s, s->have_step ? dt / s->h_last : 0.0); /* for the actual length of the attempt */
    }
    const myb_status status = sweeps(s, dt);
    if (status != MYB_OK) {
      return status;
    }
    const double proposed = proposed_step(s, dt);
    if (!(fabs(proposed) >= safety * fabs(dt))) { /* also NaN: the series diverged */
      s->stats.rejected++;
      dt = finite(proposed) ? proposed : safety * dt;
      continue;
    }
    state_at(s->n, s->y, s->y_comp, s->a0, s->b, dt, 1.0, s->y_new, s->comp_new);
    for (size_t i = 0; i < 2 * s->n; ++i) {
      if (!finite(s->y_new[i])) {
        return MYB_E_STEP_MIN; /* the solution left the doubles */
      }
    }
    accept(s, dt, last, t_limit);
    /* Growth is bounded relative to the planned step: a cut step says nothing about it. */
    const double reference = last ? larger_of(fabs(dt_planned), fabs(dt)) : fabs(dt);
    const double next = fabs(proposed) < reference / safety ? fabs(proposed) : reference / safety;
    s->h = clamp_step(s, direction * next);
    return MYB_OK;
  }
  return MYB_E_NOCONV;
}

double myb_internal_ias15_planned_step(const myb_ias15* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->started && integrator->have_step ? fabs(integrator->h) : 0.0;
}

void myb_internal_ias15_plan_step(myb_ias15* integrator, double h) {
  MYB_CHECK(integrator != NULL && integrator->started && h > 0.0);
  integrator->h = h;
}

myb_status myb_internal_ias15_retake(myb_ias15* integrator, double t_limit) {
  MYB_CHECK(integrator != NULL && integrator->started);
  myb_ias15* s = integrator;
  const size_t n = s->n;
  const double withdrawn = s->h_last;
  const double dt = (t_limit - s->t_old) - s->t_comp_old;
  if (!s->have_step || !finite(t_limit) || !(dt * withdrawn > 0.0) ||
      !(fabs(dt) < fabs(withdrawn))) {
    return MYB_E_DOMAIN;
  }
  /* back to the beginning of the step, with its polynomial as the prediction */
  memcpy(s->b_hint, s->b_last, (size_t)ORDER * n * sizeof(double));
  memcpy(s->y, s->y_old, 2 * n * sizeof(double));
  memcpy(s->y_comp, s->comp_old, 2 * n * sizeof(double));
  memcpy(s->a0, s->a_old, n * sizeof(double));
  s->t = s->t_old;
  s->t_comp = s->t_comp_old;
  s->a0_valid = 1;
  s->have_step = 0; /* until the step is taken again there is none */
  s->stats.accepted--;
  s->stats.rejected++;
  s->hint = 1;
  s->hint_h = withdrawn;
  s->hint_s0 = 0.0;
  predict_from_hint(s, dt / withdrawn);
  const myb_status status = sweeps(s, dt);
  if (status != MYB_OK) {
    return status;
  }
  const double proposed = proposed_step(s, dt);
  state_at(n, s->y, s->y_comp, s->a0, s->b, dt, 1.0, s->y_new, s->comp_new);
  for (size_t i = 0; i < 2 * n; ++i) {
    if (!finite(s->y_new[i])) {
      return MYB_E_STEP_MIN;
    }
  }
  accept(s, dt, 1, t_limit);
  s->hint = 1;
  s->hint_s0 = dt / withdrawn;
  double next = fabs(withdrawn);
  if (finite(proposed) && fabs(proposed) < next) {
    next = fabs(proposed);
  }
  s->h = clamp_step(s, withdrawn > 0.0 ? next : -next);
  return MYB_OK;
}

myb_status myb_ias15_integrate(myb_ias15* integrator, double t_end) {
  MYB_CHECK(integrator != NULL && integrator->started);
  /* bound: max_steps accepted steps */
  for (int64_t i = 0; i < integrator->options.max_steps; ++i) {
    if (integrator->t == t_end && integrator->t_comp == 0.0) {
      return MYB_OK;
    }
    const myb_status status = myb_ias15_step(integrator, t_end);
    if (status != MYB_OK) {
      return status;
    }
  }
  return integrator->t == t_end && integrator->t_comp == 0.0 ? MYB_OK : MYB_E_NOCONV;
}

double myb_ias15_time(const myb_ias15* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->t + integrator->t_comp;
}

const double* myb_ias15_state(const myb_ias15* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->y;
}

myb_ias15_stats myb_ias15_statistics(const myb_ias15* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->stats;
}

myb_status myb_ias15_dense(myb_ias15* integrator, double t, double* y) {
  MYB_CHECK(integrator != NULL && y != NULL);
  const myb_ias15* s = integrator;
  if (!s->have_step) {
    return MYB_E_RANGE;
  }
  const double lo = s->h_last > 0.0 ? s->t_old : s->t;
  const double hi = s->h_last > 0.0 ? s->t : s->t_old;
  if (!(t >= lo && t <= hi)) {
    return MYB_E_RANGE;
  }
  const double h = ((t - s->t_old) - s->t_comp_old) / s->h_last;
  state_at(s->n, s->y_old, s->comp_old, s->a_old, s->b_last, s->h_last, h, y, NULL);
  return MYB_OK;
}

static myb_status dense_view(void* integrator, double t, double* y) {
  return myb_ias15_dense(integrator, t, y);
}

myb_status myb_internal_ias15_position_box(const myb_ias15* integrator, size_t coordinate,
                                           double* lo, double* hi) {
  MYB_CHECK(integrator != NULL && lo != NULL && hi != NULL && coordinate < integrator->n);
  const myb_ias15* s = integrator;
  if (!s->have_step) {
    return MYB_E_RANGE;
  }
  enum { DEGREE = ORDER + 2 };
  /* x(h) = x0 + (dt v0) h + (dt^2 a0 / 2) h^2 + sum_k dt^2 b_k w_k h^(k+3) */
  const size_t n = s->n;
  const double dt = s->h_last;
  const double* b = s->b_last + (size_t)ORDER * coordinate;
  double c[DEGREE + 1];
  c[0] = s->y_old[coordinate] + s->comp_old[coordinate];
  c[1] = dt * s->y_old[n + coordinate];
  c[2] = dt * dt * 0.5 * s->a_old[coordinate];
  for (int k = 0; k < ORDER; ++k) {
    c[k + 3] = dt * dt * b[k] * x_weight[k];
  }
  /* Bernstein coefficients: beta_j = sum_{i <= j} C(j, i) / C(9, i) c_i */
  static const double choose9[DEGREE + 1] = {1, 9, 36, 84, 126, 126, 84, 36, 9, 1};
  double low = c[0];
  double high = c[0];
  double size = fabs(c[0]);
  for (int j = 1; j <= DEGREE; ++j) {
    double beta = 0.0;
    double choose_j = 1.0; /* C(j, i) */
    for (int i = 0; i <= j; ++i) {
      beta += choose_j / choose9[i] * c[i];
      choose_j = choose_j * (double)(j - i) / (double)(i + 1);
    }
    low = fmin(low, beta);
    high = fmax(high, beta);
    size = fmax(size, fabs(beta));
  }
  for (int i = 1; i <= DEGREE; ++i) {
    size = fmax(size, fabs(c[i]) * 10.0);
  }
  const double margin = 64.0 * 2.220446049250313e-16 * size;
  *lo = low - margin;
  *hi = high + margin;
  return MYB_OK;
}

myb_status myb_ias15_locate(myb_ias15* integrator, myb_event_function g, void* context,
                            int* found, double* t_event, double* y_event) {
  MYB_CHECK(integrator != NULL && g != NULL && found != NULL && t_event != NULL);
  myb_ias15* s = integrator;
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
