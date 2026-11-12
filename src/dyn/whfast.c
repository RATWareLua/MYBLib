/* SPDX-License-Identifier: BSD-3-Clause */
/* Wisdom-Holman splitting in Jacobi coordinates (Wisdom & Holman 1991; Rein & Tamayo 2015).
 *
 * Jacobi coordinates of the massive bodies 0 .. m-1 with eta_i = GM_0 + ... + GM_i:
 *   q_0 = (sum_k GM_k x_k) / eta_(m-1)                (centre of mass)
 *   q_i = x_i - R_(i-1),  R_(i-1) = (sum_(k<i) GM_k x_k) / eta_(i-1),   i >= 1,
 * and test bodies q_t = x_t - q_0; velocities likewise. The inverse runs outwards in:
 *   R_(m-1) = q_0,  R_(i-1) = R_i - (GM_i / eta_i) q_i,  x_i = R_(i-1) + q_i,  x_0 = R_0.
 * With masses m_i = GM_i / G the kinetic energy is diagonal in these coordinates (Jacobi
 * masses m_i eta_(i-1) / eta_i), and H = sum_i H_Kepler,i + H_interaction, where body i >= 1
 * moves on a Kepler orbit of GM eta_i (test bodies: eta_(m-1)) and the centre of mass drifts.
 * The interaction is the potential less the Kepler terms, so the Jacobi acceleration of the
 * kick is the linear transform of the inertial N-body accelerations plus the Kepler term:
 *   a'_i = a_i - (sum_(k<i) GM_k a_k) / eta_(i-1) + eta_i q_i / |q_i|^3,
 *   a'_0 = (sum_k GM_k a_k) / eta_(m-1)    (exactly 0 for a reciprocal graph).
 * The Jacobi state is accumulated with compensated summation: the drift returns the increment
 * of the state (kepler/drift.h), the kick the increment of the velocity. */
#include <math.h> /* sqrt, fabs */
#include <stdckdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "dyn/nbody_internal.h"
#include "dyn/wh_coefficients.h"
#include "kepler/drift.h"
#include "myb/dyn.h"
#include "num/eft.h"

/* Direction of the corrector: the kernel state is C^-1 applied to the state, C(+b) maps the
 * kernel back (fixed by the order tests of tests/unit/test_whfast.c). */
static const double corrector_sign = 1.0;

static const double max_steps = 0x1p53;

struct myb_wh {
  myb_allocator allocator;
  myb_nbody* system;
  size_t n; /* bodies */
  size_t m; /* massive bodies */
  const myb_wh_table* table;
  int stages; /* corrector: Z pairs per half, 0 for none */
  double h;
  double t0;
  int64_t k; /* node index */
  int started;
  int ja_valid; /* ja holds the kick acceleration at the current positions */
  int reciprocal;
  myb_wh_stats stats;
  double* eta; /* m */
  double* q;   /* 3n each: Jacobi positions and velocities with their compensations */
  double* q_lo;
  double* p;
  double* p_lo;
  double* save; /* 12n: q, q_lo, p, p_lo */
  double* x;    /* 3n: inertial positions (kick, output) */
  double* a;    /* 3n: inertial accelerations */
  double* ja;   /* 3n: Jacobi accelerations */
};

static int finite(double x) {
  return x - x == 0.0;
}

myb_status myb_wh_create(myb_nbody* system, const myb_wh_options* options,
                         const myb_allocator* allocator, myb_wh** wh) {
  MYB_CHECK(system != NULL && options != NULL && wh != NULL);
  const myb_wh_options o = *options;
  if (!finite(o.h) || o.h == 0.0 || (int)o.scheme < 0 || (int)o.scheme >= MYB_WH_SCHEMES) {
    return MYB_E_DOMAIN;
  }
  int stages = 0;
  if (o.corrector != 0) {
    if (o.scheme != MYB_WH_DKD || o.corrector < 3 || o.corrector > 17 || o.corrector % 2 == 0) {
      return MYB_E_DOMAIN;
    }
    stages = (o.corrector - 1) / 2;
  }
  const size_t m = myb_nbody_massive(system);
  const size_t n = myb_nbody_count(system);
  if (!(myb_nbody_gm(system, 0) > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator al = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&al, 1, sizeof(myb_wh), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_wh* w = block;
  *w = (myb_wh){.allocator = al,
                .system = system,
                .n = n,
                .m = m,
                .table = &myb_wh_tables[o.scheme],
                .stages = stages,
                .h = o.h};
  size_t count = 0;
  if (ckd_mul(&count, n, (size_t)33) || ckd_add(&count, count, m)) {
    status = MYB_E_OVERFLOW;
    goto fail;
  }
  status = myb_internal_alloc_array(&al, count, sizeof(double), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  double* d = block;
  memset(d, 0, count * sizeof(double));
  const size_t v = 3 * n;
  w->q = d;
  w->q_lo = d + v;
  w->p = d + 2 * v;
  w->p_lo = d + 3 * v;
  w->save = d + 4 * v;
  w->x = d + 8 * v;
  w->a = d + 9 * v;
  w->ja = d + 10 * v;
  w->eta = d + 11 * v;
  double eta = 0.0;
  for (size_t i = 0; i < m; ++i) {
    eta += myb_nbody_gm(system, i);
    w->eta[i] = eta;
  }
  *wh = w;
  return MYB_OK;

fail:
  myb_internal_release(&al, w);
  return status;
}

void myb_wh_destroy(myb_wh* wh) {
  if (wh == NULL) {
    return;
  }
  const myb_allocator a = wh->allocator;
  myb_internal_release(&a, wh->q);
  myb_internal_release(&a, wh);
}

/* ---- Coordinates --------------------------------------------------------------------------- */

/* Inertial (x, 3 per body) to Jacobi (q); the same transform serves velocities. */
static void to_jacobi(const myb_wh* w, const double* x, double* q) {
  double s[3];
  const double g0 = myb_nbody_gm(w->system, 0);
  for (int c = 0; c < 3; ++c) {
    s[c] = g0 * x[c];
  }
  for (size_t i = 1; i < w->m; ++i) {
    const double g = myb_nbody_gm(w->system, i);
    for (int c = 0; c < 3; ++c) {
      const size_t at = 3 * i + (size_t)c;
      q[at] = x[at] - s[c] / w->eta[i - 1];
      s[c] += g * x[at];
    }
  }
  const double total = w->eta[w->m - 1];
  for (int c = 0; c < 3; ++c) {
    q[c] = s[c] / total;
  }
  for (size_t i = w->m; i < w->n; ++i) {
    for (int c = 0; c < 3; ++c) {
      q[3 * i + (size_t)c] = x[3 * i + (size_t)c] - q[c];
    }
  }
}

/* Jacobi (q) to inertial (x). */
static void from_jacobi(const myb_wh* w, const double* q, double* x) {
  double r[3] = {q[0], q[1], q[2]};
  for (size_t i = w->m; i < w->n; ++i) {
    for (int c = 0; c < 3; ++c) {
      x[3 * i + (size_t)c] = r[c] + q[3 * i + (size_t)c];
    }
  }
  for (size_t i = w->m - 1; i >= 1; --i) {
    const double ratio = myb_nbody_gm(w->system, i) / w->eta[i];
    for (int c = 0; c < 3; ++c) {
      const size_t at = 3 * i + (size_t)c;
      r[c] -= ratio * q[at];
      x[at] = r[c] + q[at];
    }
  }
  for (int c = 0; c < 3; ++c) {
    x[c] = r[c];
  }
}

static void add_compensated(double* y, double* lo, size_t i, double delta) {
  y[i] = myb_eft_two_sum(y[i], delta + lo[i], &lo[i]);
}

/* ---- The two flows -------------------------------------------------------------------------- */

/* Kepler drifts of every Jacobi body, free drift of the centre of mass. */
static myb_status drift(myb_wh* w, double dt) {
  w->stats.drifts++;
  w->ja_valid = 0;
  for (int c = 0; c < 3; ++c) {
    add_compensated(w->q, w->q_lo, (size_t)c, dt * w->p[c]);
  }
  const double total = w->eta[w->m - 1];
  for (size_t i = 1; i < w->n; ++i) {
    const double mu = i < w->m ? w->eta[i] : total;
    myb_pv pv0;
    for (int c = 0; c < 3; ++c) {
      pv0.r[c] = w->q[3 * i + (size_t)c];
      pv0.v[c] = w->p[3 * i + (size_t)c];
    }
    myb_pv delta;
    const myb_status status = myb_internal_kepler_delta(mu, &pv0, dt, &delta);
    if (status != MYB_OK) {
      return status;
    }
    for (int c = 0; c < 3; ++c) {
      add_compensated(w->q, w->q_lo, 3 * i + (size_t)c, delta.r[c]);
      add_compensated(w->p, w->p_lo, 3 * i + (size_t)c, delta.v[c]);
    }
  }
  return MYB_OK;
}

static void add_kepler_term(double* ja, const double* q, size_t i, double mu) {
  const double* r = q + 3 * i;
  const double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
  const double k = mu / (r2 * sqrt(r2));
  for (int c = 0; c < 3; ++c) {
    ja[3 * i + (size_t)c] += k * r[c];
  }
}

/* Jacobi accelerations of the kick at the current positions (time t). */
static myb_status kick_acceleration(myb_wh* w, double t) {
  from_jacobi(w, w->q, w->x);
  w->stats.evaluations++;
  const myb_status status = myb_nbody_acceleration(w->system, t, w->x, w->a);
  if (status != MYB_OK) {
    return status;
  }
  const double* a = w->a;
  double* ja = w->ja;
  double s[3];
  const double g0 = myb_nbody_gm(w->system, 0);
  for (int c = 0; c < 3; ++c) {
    s[c] = g0 * a[c];
  }
  for (size_t i = 1; i < w->m; ++i) {
    const double g = myb_nbody_gm(w->system, i);
    for (int c = 0; c < 3; ++c) {
      const size_t at = 3 * i + (size_t)c;
      ja[at] = a[at] - s[c] / w->eta[i - 1];
      s[c] += g * a[at];
    }
    add_kepler_term(ja, w->q, i, w->eta[i]);
  }
  const double total = w->eta[w->m - 1];
  double center[3] = {0.0, 0.0, 0.0};
  if (!w->reciprocal) {
    for (int c = 0; c < 3; ++c) {
      center[c] = s[c] / total;
    }
  }
  for (int c = 0; c < 3; ++c) {
    ja[c] = center[c];
  }
  for (size_t i = w->m; i < w->n; ++i) {
    for (int c = 0; c < 3; ++c) {
      const size_t at = 3 * i + (size_t)c;
      ja[at] = a[at] - center[c];
    }
    add_kepler_term(ja, w->q, i, total);
  }
  for (size_t i = 0; i < 3 * w->n; ++i) {
    if (!finite(ja[i])) {
      return MYB_E_DOMAIN;
    }
  }
  w->ja_valid = 1;
  return MYB_OK;
}

static myb_status kick(myb_wh* w, double t, double dt) {
  if (!w->ja_valid) {
    const myb_status status = kick_acceleration(w, t);
    if (status != MYB_OK) {
      return status;
    }
  }
  for (size_t i = 0; i < 3 * w->n; ++i) {
    add_compensated(w->p, w->p_lo, i, dt * w->ja[i]);
  }
  return MYB_OK;
}

/* One step of the scheme of length dt from time t. */
static myb_status scheme_step(myb_wh* w, double t, double dt) {
  const myb_wh_table* s = w->table;
  const int drifts = s->kick_first ? s->kicks - 1 : s->kicks + 1;
  const int operations = s->kicks + drifts;
  double elapsed = 0.0; /* fraction of the step done by the drifts */
  int next_drift = 0;
  int next_kick = 0;
  for (int op = 0; op < operations; ++op) {
    const int is_kick = (op % 2 == 0) == (s->kick_first != 0);
    myb_status status = MYB_OK;
    if (is_kick) {
      status = kick(w, t + elapsed * dt, s->kick[next_kick] * dt);
      ++next_kick;
    } else {
      status = drift(w, s->drift[next_drift] * dt);
      elapsed += s->drift[next_drift];
      ++next_drift;
    }
    if (status != MYB_OK) {
      return status;
    }
  }
  return MYB_OK;
}

/* Z(a, b) = A(a) B(b) A(-2a) B(-b) A(a), fractions of hc. */
static myb_status corrector_z(myb_wh* w, double t, double hc, double a, double b) {
  myb_status status = drift(w, a * hc);
  if (status == MYB_OK) {
    status = kick(w, t, b * hc);
  }
  if (status == MYB_OK) {
    status = drift(w, -2.0 * a * hc);
  }
  if (status == MYB_OK) {
    status = kick(w, t, -b * hc);
  }
  if (status == MYB_OK) {
    status = drift(w, a * hc);
  }
  return status;
}

/* C(sign): Z(a_m, s b_m) ... Z(a_1, s b_1) Z(-a_1, -s b_1) ... Z(-a_m, -s b_m) for step hc;
 * C(-sign) is its inverse. */
static myb_status corrector(myb_wh* w, double t, double hc, double sign) {
  const int m = w->stages;
  if (m == 0) {
    return MYB_OK;
  }
  const double* b = myb_wh_corrector_b[m - 1];
  const double* a = myb_wh_corrector_a;
  for (int i = m - 1; i >= 0; --i) {
    const myb_status status = corrector_z(w, t, hc, a[i], sign * b[i]);
    if (status != MYB_OK) {
      return status;
    }
  }
  for (int i = 0; i < m; ++i) {
    const myb_status status = corrector_z(w, t, hc, -a[i], -sign * b[i]);
    if (status != MYB_OK) {
      return status;
    }
  }
  return MYB_OK;
}

static void save_state(myb_wh* w) {
  memcpy(w->save, w->q, 12 * w->n * sizeof(double)); /* q, q_lo, p, p_lo are contiguous */
}

static void restore_state(myb_wh* w) {
  memcpy(w->q, w->save, 12 * w->n * sizeof(double));
  w->ja_valid = 0;
}

/* ---- Public operations ---------------------------------------------------------------------- */

myb_status myb_wh_start(myb_wh* wh, double t0, const double* x, const double* v) {
  MYB_CHECK(wh != NULL && x != NULL && v != NULL);
  myb_wh* w = wh;
  w->started = 0;
  /* the splitting needs accelerations that depend on the positions only, and exact: the
   * error of an approximate sum is relative to the whole force, the central body's included,
   * which the splitting subtracts */
  if (!finite(t0) || myb_internal_nbody_velocity_dependent(w->system) ||
      myb_nbody_approximation(w->system, NULL) != 0) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < 3 * w->n; ++i) {
    if (!finite(x[i]) || !finite(v[i])) {
      return MYB_E_DOMAIN;
    }
  }
  w->reciprocal = myb_internal_nbody_reciprocal(w->system);
  to_jacobi(w, x, w->q);
  to_jacobi(w, v, w->p);
  memset(w->q_lo, 0, 3 * w->n * sizeof(double));
  memset(w->p_lo, 0, 3 * w->n * sizeof(double));
  w->t0 = t0;
  w->k = 0;
  w->stats = (myb_wh_stats){};
  w->ja_valid = 0;
  const myb_status status = corrector(w, t0, w->h, -corrector_sign);
  if (status != MYB_OK) {
    return status;
  }
  w->started = 1;
  return MYB_OK;
}

double myb_wh_time(const myb_wh* wh) {
  MYB_CHECK(wh != NULL);
  return wh->t0 + (double)wh->k * wh->h;
}

/* One kernel step of direction (+1 or -1); the state is unchanged on error. */
static myb_status kernel_step(myb_wh* w, int direction) {
  const double dt = direction > 0 ? w->h : -w->h;
  save_state(w);
  const myb_status status = scheme_step(w, myb_wh_time(w), dt);
  if (status != MYB_OK) {
    restore_state(w);
    return status;
  }
  w->k += direction;
  w->stats.steps++;
  return MYB_OK;
}

myb_status myb_wh_step(myb_wh* wh) {
  MYB_CHECK(wh != NULL && wh->started);
  if (fabs((double)wh->k + 1.0) > max_steps) {
    return MYB_E_DOMAIN;
  }
  return kernel_step(wh, 1);
}

myb_status myb_wh_integrate(myb_wh* wh, double t_end) {
  MYB_CHECK(wh != NULL && wh->started);
  myb_wh* w = wh;
  if (!finite(t_end)) {
    return MYB_E_DOMAIN;
  }
  const double target = __builtin_roundeven((t_end - w->t0) / w->h);
  if (!(fabs(target) <= max_steps)) {
    return MYB_E_DOMAIN;
  }
  const int64_t k_end = (int64_t)target;
  const int64_t steps = k_end > w->k ? k_end - w->k : w->k - k_end;
  const int direction = k_end > w->k ? 1 : -1;
  /* bound: |k_end - k| <= 2^54 steps, each moving k one node towards k_end */
  for (int64_t i = 0; i < steps; ++i) {
    const myb_status status = kernel_step(w, direction);
    if (status != MYB_OK) {
      return status;
    }
  }
  return MYB_OK;
}

/* Corrected state at the node, then one step of length dt (corrected as well) when dt != 0;
 * the result to x, v; the kernel state is restored. */
static myb_status output(myb_wh* w, double dt, double* x, double* v) {
  const double t = myb_wh_time(w);
  save_state(w);
  myb_status status = corrector(w, t, w->h, corrector_sign);
  if (status == MYB_OK && dt != 0.0) {
    status = corrector(w, t, dt, -corrector_sign);
    if (status == MYB_OK) {
      status = scheme_step(w, t, dt);
    }
    if (status == MYB_OK) {
      status = corrector(w, t + dt, dt, corrector_sign);
    }
  }
  if (status == MYB_OK) {
    from_jacobi(w, w->q, x);
    from_jacobi(w, w->p, v);
  }
  restore_state(w);
  return status;
}

myb_status myb_wh_state(myb_wh* wh, double* x, double* v) {
  MYB_CHECK(wh != NULL && wh->started && x != NULL && v != NULL);
  return output(wh, 0.0, x, v);
}

myb_status myb_wh_dense(myb_wh* wh, double t, double* x, double* v) {
  MYB_CHECK(wh != NULL && wh->started && x != NULL && v != NULL);
  if (!finite(t)) {
    return MYB_E_DOMAIN;
  }
  const double dt = t - myb_wh_time(wh);
  if (!(fabs(dt) <= fabs(wh->h))) {
    return MYB_E_RANGE;
  }
  return output(wh, dt, x, v);
}

myb_wh_stats myb_wh_statistics(const myb_wh* wh) {
  MYB_CHECK(wh != NULL);
  return wh->stats;
}
