/* SPDX-License-Identifier: BSD-3-Clause */
/* The direct method of Sims and Flanagan (myb/optimal.h): impulses in the middle of equal
 * segments, conics between them — or the arcs of a force model — the two halves of the
 * transfer meeting in the middle. */
#include <math.h> /* sqrt, fabs */
#include <string.h>

#include "core/alloc.h"
#include "core/check.h"
#include "myb/cmath.h"
#include "myb/kepler.h"
#include "myb/optimal.h"
#include "num/vec3.h"
#include "transfer/leg_internal.h"

enum {
  MAX_SEGMENTS = MYB_SIMS_FLANAGAN_MAX_SEGMENTS,
  COAST = 42,           /* a state and its transition matrix */
  COAST_STEPS = 1000000 /* of the integrator, for one coast */
};

/* The tolerance of a coast in a field, relative and absolute in the units of the solver. */
static const double coast_tolerance = 1e-13;

/* ---- Coasts ---------------------------------------------------------------------------------- */

/* The field at (t, r, v) in the units of the solver, and when wanted its derivatives by r and
 * v (row-major). */
static myb_status field_values(const myb_internal_coast* coast, double t, const double r[3],
                               const double v[3], double g[3], double dgdr[9], double dgdv[9]) {
  const double speed = coast->length / coast->time;
  const double x[3] = {r[0] * coast->length, r[1] * coast->length, r[2] * coast->length};
  const double w[3] = {v[0] * speed, v[1] * speed, v[2] * speed};
  const double at = coast->origin + t * coast->time;
  double a[3];
  myb_status st = myb_force_model_acceleration(coast->field, at, x, w, a);
  if (st != MYB_OK) {
    return st;
  }
  const double scale = coast->time * coast->time / coast->length;
  for (int i = 0; i < 3; ++i) {
    g[i] = a[i] * scale;
  }
  if (dgdr != NULL) {
    st = myb_force_model_jacobian(coast->field, at, x, w, dgdr, dgdv, NULL);
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 9; ++i) {
      dgdr[i] *= coast->time * coast->time;
      dgdv[i] *= coast->time;
    }
  }
  return MYB_OK;
}

/* The state and its variational equations: Phi' = A Phi, A = [0 I; dg/dr dg/dv]. */
static myb_status coast_rhs(void* context, double t, const double* z, double* dz) {
  const myb_internal_coast* coast = context;
  double g[3];
  double dgdr[9];
  double dgdv[9];
  const myb_status st = field_values(coast, t, z, z + 3, g, dgdr, dgdv);
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 3; ++i) {
    dz[i] = z[3 + i];
    dz[3 + i] = g[i];
  }
  const double* phi = z + 6;
  double* rate = dz + 6;
  for (int j = 0; j < 6; ++j) {
    for (int i = 0; i < 3; ++i) {
      rate[6 * i + j] = phi[6 * (3 + i) + j];
      double sum = 0.0;
      for (int k = 0; k < 3; ++k) {
        sum += dgdr[3 * i + k] * phi[6 * k + j] + dgdv[3 * i + k] * phi[6 * (3 + k) + j];
      }
      rate[6 * (3 + i) + j] = sum;
    }
  }
  return MYB_OK;
}

myb_status myb_internal_coast_open(myb_internal_coast* coast, myb_force_model* field,
                                   double length, double time, double origin,
                                   const myb_allocator* allocator) {
  memset(coast, 0, sizeof *coast);
  coast->field = field;
  coast->length = length;
  coast->time = time;
  coast->origin = origin;
  if (field == NULL) {
    return MYB_OK;
  }
  const myb_dop853_options options = {.rtol = coast_tolerance, .atol = coast_tolerance,
                                      .max_steps = COAST_STEPS};
  return myb_dop853_create(COAST, &options, allocator, &coast->integrator);
}

void myb_internal_coast_close(myb_internal_coast* coast) {
  myb_dop853_destroy(coast->integrator);
  coast->integrator = NULL;
}

myb_status myb_internal_coast_fly(myb_internal_coast* coast, double t, const double y0[6],
                                  double dt, double y[6], double* phi) {
  if (coast->field == NULL) {
    myb_pv from;
    myb_pv to;
    memcpy(from.r, y0, sizeof from.r);
    memcpy(from.v, y0 + 3, sizeof from.v);
    const myb_status st = phi != NULL ? myb_kepler_stm(1.0, &from, dt, &to, phi)
                                      : myb_kepler_propagate(1.0, &from, dt, &to);
    if (st == MYB_OK) {
      memcpy(y, to.r, sizeof to.r);
      memcpy(y + 3, to.v, sizeof to.v);
    }
    return st;
  }
  /* the transition matrix always goes along: the steps, and with them the state, are then
   * the same whether it is asked for or not */
  double z[COAST];
  memcpy(z, y0, 6 * sizeof *z);
  memset(z + 6, 0, 36 * sizeof *z);
  for (int i = 0; i < 6; ++i) {
    z[6 + 7 * i] = 1.0;
  }
  myb_status st = myb_dop853_start(coast->integrator, coast_rhs, coast, t, z);
  if (st == MYB_OK) {
    st = myb_dop853_integrate(coast->integrator, t + dt);
  }
  if (st != MYB_OK) {
    return st;
  }
  const double* end = myb_dop853_state(coast->integrator);
  memcpy(y, end, 6 * sizeof *y);
  if (phi != NULL) {
    memcpy(phi, end + 6, 36 * sizeof *phi);
  }
  return MYB_OK;
}

myb_status myb_internal_coast_rate(myb_internal_coast* coast, double t, const double y[6],
                                   double rate[6]) {
  memcpy(rate, y + 3, 3 * sizeof *rate);
  if (coast->field == NULL) {
    const double d = myb_vec3_norm(y);
    if (!(d > 0.0) || !isfinite(d)) {
      return MYB_E_DOMAIN;
    }
    const double k = -1.0 / (d * d * d);
    for (int i = 0; i < 3; ++i) {
      rate[3 + i] = k * y[i];
    }
    return MYB_OK;
  }
  return field_values(coast, t, y, y + 3, rate + 3, NULL, NULL);
}

/* ---- The sweep of a leg ---------------------------------------------------------------------- */

size_t myb_internal_sweep_work(size_t n) {
  return n * (2 * 36 + 6);
}

/* p = p m, 6 x 6. */
static void right_multiply(double p[36], const double m[36]) {
  double out[36];
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double sum = 0.0;
      for (int k = 0; k < 6; ++k) {
        sum += p[6 * i + k] * m[6 * k + j];
      }
      out[6 * i + j] = sum;
    }
  }
  memcpy(p, out, sizeof out);
}

/* y = m x, 6 x 6 by 6. */
static void apply(const double m[36], const double x[6], double y[6]) {
  for (int i = 0; i < 6; ++i) {
    double sum = 0.0;
    for (int k = 0; k < 6; ++k) {
      sum += m[6 * i + k] * x[k];
    }
    y[i] = sum;
  }
}

/* The block of segment i of g (6 x 3 n): the columns of the velocity of p. */
static void plant_columns(const double p[36], double* g, size_t columns, size_t i) {
  for (size_t r = 0; r < 6; ++r) {
    for (size_t c = 0; c < 3; ++c) {
      g[r * columns + 3 * i + c] = p[6 * r + 3 + c];
    }
  }
}

/* The impulse dv added to the velocity of y (times sign) at the time t, and the rate of the
 * state before it minus the rate after it: what the impulse made later does, through the
 * transition matrix from there. */
static myb_status jump(myb_internal_coast* coast, double t, double y[6], const double dv[3],
                       double sign, double out[6]) {
  double before[6];
  double after[6];
  myb_status st = myb_internal_coast_rate(coast, t, y, before);
  for (int k = 0; k < 3; ++k) {
    y[3 + k] += sign * dv[k];
  }
  if (st == MYB_OK) {
    st = myb_internal_coast_rate(coast, t, y, after);
  }
  for (int k = 0; k < 6 && st == MYB_OK; ++k) {
    out[k] = before[k] - after[k];
  }
  return st;
}

/* The identity, 6 x 6. */
static void identity(double p[36]) {
  memset(p, 0, 36 * sizeof *p);
  for (int i = 0; i < 6; ++i) {
    p[7 * i] = 1.0;
  }
}

/* The two halves flown, c = forward minus backward; the transition matrices of the coasts in
 * phis (2 per segment, in the order of time) and the jumps of the rates at the impulses
 * (forward: before minus after; backward: after minus before) when asked for. */
static myb_status fly_halves(myb_internal_coast* coast, const myb_internal_leg* leg,
                             double* phis, double* jumps, double ahead[6], double behind[6]) {
  const double half = 0.5 * leg->dt;
  myb_status st = MYB_OK;
  memcpy(ahead, leg->start, 6 * sizeof *ahead);
  /* bound: the segments flown forward */
  for (size_t i = 0; i < leg->forward && st == MYB_OK; ++i) {
    const double from = leg->t + (double)i * leg->dt;
    const double middle = leg->t + ((double)i + 0.5) * leg->dt;
    st = myb_internal_coast_fly(coast, from, ahead, half, ahead,
                                phis != NULL ? phis + 72 * i : NULL);
    if (st == MYB_OK && jumps != NULL) {
      st = jump(coast, middle, ahead, leg->impulses + 3 * i, 1.0, jumps + 6 * i);
    } else if (st == MYB_OK) {
      for (size_t k = 0; k < 3; ++k) {
        ahead[3 + k] += leg->impulses[3 * i + k];
      }
    }
    if (st == MYB_OK) {
      st = myb_internal_coast_fly(coast, middle, ahead, half, ahead,
                                  phis != NULL ? phis + 72 * i + 36 : NULL);
    }
  }
  memcpy(behind, leg->end, 6 * sizeof *behind);
  /* bound: the segments flown backward */
  for (size_t i = leg->n; i-- > leg->forward && st == MYB_OK;) {
    const double to = leg->t + ((double)i + 1.0) * leg->dt;
    const double middle = leg->t + ((double)i + 0.5) * leg->dt;
    st = myb_internal_coast_fly(coast, to, behind, -half, behind,
                                phis != NULL ? phis + 72 * i + 36 : NULL);
    if (st == MYB_OK && jumps != NULL) {
      st = jump(coast, middle, behind, leg->impulses + 3 * i, -1.0, jumps + 6 * i);
    } else if (st == MYB_OK) {
      for (size_t k = 0; k < 3; ++k) {
        behind[3 + k] -= leg->impulses[3 * i + k];
      }
    }
    if (st == MYB_OK) {
      st = myb_internal_coast_fly(coast, middle, behind, -half, behind,
                                  phis != NULL ? phis + 72 * i : NULL);
    }
  }
  return st;
}

myb_status myb_internal_leg_sweep(myb_internal_coast* coast, const myb_internal_leg* leg,
                                  double* g, int full, double* work, myb_internal_sweep* out) {
  const size_t n = leg->n;
  const size_t forward = leg->forward;
  const size_t columns = 3 * n;
  full = full && g != NULL;
  double* phis = g != NULL ? work : NULL;
  double* jumps = full ? work + 72 * n : NULL;
  double ahead[6];
  double behind[6];
  myb_status st = fly_halves(coast, leg, phis, jumps, ahead, behind);
  if (st != MYB_OK) {
    return st;
  }
  for (int k = 0; k < 6; ++k) {
    out->c[k] = ahead[k] - behind[k];
  }
  if (g == NULL) {
    return MYB_OK;
  }
  /* the rates at the ends and at the meeting point, for the derivatives by the times */
  const double meet = leg->t + (double)forward * leg->dt;
  double rate_start[6] = {0};
  double rate_end[6] = {0};
  double rate_ahead[6] = {0};
  double rate_behind[6] = {0};
  if (full) {
    st = myb_internal_coast_rate(coast, leg->t, leg->start, rate_start);
    if (st == MYB_OK) {
      st = myb_internal_coast_rate(coast, leg->t + (double)n * leg->dt, leg->end, rate_end);
    }
    if (st == MYB_OK) {
      st = myb_internal_coast_rate(coast, meet, ahead, rate_ahead);
    }
    if (st == MYB_OK) {
      st = myb_internal_coast_rate(coast, meet, behind, rate_behind);
    }
    if (st != MYB_OK) {
      return st;
    }
  }
  /* from the meeting point back to the start: p = Phi(meet, the time reached); an impulse made
   * later moves the meeting point by p (rate before - rate after) */
  double p[36];
  double moved[6];
  double step_ahead[6] = {0};
  double shift_ahead[6] = {0};
  identity(p);
  /* bound: the segments flown forward */
  for (size_t i = forward; i-- > 0;) {
    right_multiply(p, phis + 72 * i + 36);
    plant_columns(p, g, columns, i);
    if (full) {
      apply(p, jumps + 6 * i, moved);
      for (int k = 0; k < 6; ++k) {
        step_ahead[k] += ((double)i + 0.5) * moved[k];
        shift_ahead[k] += moved[k];
      }
    }
    right_multiply(p, phis + 72 * i);
  }
  memcpy(out->by_start, p, sizeof out->by_start);
  if (full) {
    apply(p, rate_start, moved);
    for (int k = 0; k < 6; ++k) {
      step_ahead[k] += (double)forward * rate_ahead[k];
      shift_ahead[k] += rate_ahead[k] - moved[k];
    }
  }
  /* from the meeting point on to the end: the backward half took each impulse off, so the
   * mismatch grows with it by p (velocity columns), and moves by -p (rate after - before)
   * when it is made later */
  double step_behind[6] = {0};
  double shift_behind[6] = {0};
  identity(p);
  /* bound: the segments flown backward */
  for (size_t i = forward; i < n; ++i) {
    right_multiply(p, phis + 72 * i);
    plant_columns(p, g, columns, i);
    if (full) {
      apply(p, jumps + 6 * i, moved);
      for (int k = 0; k < 6; ++k) {
        step_behind[k] += ((double)i + 0.5) * moved[k];
        shift_behind[k] += moved[k];
      }
    }
    right_multiply(p, phis + 72 * i + 36);
  }
  for (int k = 0; k < 36; ++k) {
    out->by_end[k] = -p[k];
  }
  if (full) {
    apply(p, rate_end, moved);
    for (int k = 0; k < 6; ++k) {
      step_behind[k] += (double)forward * rate_behind[k] - (double)n * moved[k];
      shift_behind[k] += rate_behind[k] - moved[k];
      out->by_step[k] = step_ahead[k] - step_behind[k];
      /* a field that does not depend on time: the leg moved in time is the same leg */
      out->by_shift[k] = coast->field != NULL ? shift_ahead[k] - shift_behind[k] : 0.0;
    }
  }
  return MYB_OK;
}

/* ---- A transfer ------------------------------------------------------------------------------ */

/* A transfer in the units of the solver. */
typedef struct leg {
  size_t n;       /* segments */
  size_t forward; /* of them flown from the start */
  double length;  /* m */
  double time;    /* s */
  double mass;    /* kg */
  double dt;      /* of a segment */
  double impulse; /* thrust dt: the velocity a full segment gives a unit of mass */
  double rate;    /* thrust dt / exhaust: the mass a full segment spends */
  myb_pv start;
  myb_pv target;
  myb_force_model* field; /* NULL: the conics */
  double t0;              /* s */
} leg;

static myb_status leg_open(leg* l, const myb_sims_flanagan* p) {
  double sum = p->mu + p->thrust + p->exhaust + p->t0 + p->mass0 + p->t1;
  for (int i = 0; i < 6; ++i) {
    sum += p->start[i] + p->target[i];
  }
  if (!isfinite(sum) || !(p->mu > 0.0) || !(p->thrust > 0.0) || !(p->exhaust > 0.0) ||
      !(p->mass0 > 0.0) || !(p->t1 > p->t0) || p->segments < 2 || p->segments > MAX_SEGMENTS ||
      !(myb_vec3_norm(p->start) > 0.0) || !(myb_vec3_norm(p->target) > 0.0)) {
    return MYB_E_DOMAIN;
  }
  l->n = p->segments;
  l->forward = (p->segments + 1) / 2;
  l->field = p->field;
  l->t0 = p->t0;
  l->length = myb_vec3_norm(p->start);
  l->time = sqrt(l->length * l->length * l->length / p->mu);
  l->mass = p->mass0;
  l->dt = (p->t1 - p->t0) / l->time / (double)p->segments;
  l->impulse = p->thrust * l->time * l->time / (l->length * l->mass) * l->dt;
  l->rate = l->impulse / (p->exhaust * l->time / l->length);
  for (int i = 0; i < 3; ++i) {
    l->start.r[i] = p->start[i] / l->length;
    l->start.v[i] = p->start[3 + i] * l->time / l->length;
    l->target.r[i] = p->target[i] / l->length;
    l->target.v[i] = p->target[3 + i] * l->time / l->length;
  }
  return MYB_OK;
}

/* The 6 x 3 block of the segment in g (6 x columns): multiplied from the left by m. */
static void carry(const double m[36], double* g, size_t columns, size_t segment) {
  double block[18];
  for (size_t r = 0; r < 6; ++r) {
    for (size_t c = 0; c < 3; ++c) {
      double sum = 0.0;
      for (size_t k = 0; k < 6; ++k) {
        sum += m[6 * r + k] * g[k * columns + 3 * segment + c];
      }
      block[3 * r + c] = sum;
    }
  }
  for (size_t r = 0; r < 6; ++r) {
    for (size_t c = 0; c < 3; ++c) {
      g[r * columns + 3 * segment + c] = block[3 * r + c];
    }
  }
}

/* The same block set to what an impulse does to the state: nothing to the position, itself
 * to the velocity. */
static void plant(double* g, size_t columns, size_t segment) {
  for (size_t r = 0; r < 6; ++r) {
    for (size_t c = 0; c < 3; ++c) {
      g[r * columns + 3 * segment + c] = r == 3 + c ? 1.0 : 0.0;
    }
  }
}

/* Half a segment of the conic from `from`; the transition matrix when wanted. */
static myb_status coast(const myb_pv* from, double dt, myb_pv* to, double* phi) {
  return phi != NULL ? myb_kepler_stm(1.0, from, dt, to, phi)
                     : myb_kepler_propagate(1.0, from, dt, to);
}

/* The mismatch c of the two halves for the sizes u and the unit directions d (3 per segment)
 * of the throttles; g (6 x 3 n; may be NULL): its derivatives by the impulses. */
static myb_status meet(const leg* l, const double* u, const double* d, double c[6], double* g) {
  const size_t n = l->n;
  const size_t columns = 3 * n;
  double total = 0.0;
  for (size_t i = 0; i < n; ++i) {
    total += u[i];
  }
  if (!(1.0 - l->rate * total > 0.0)) {
    return MYB_E_DOMAIN; /* more propellant than mass */
  }
  if (g != NULL) {
    memset(g, 0, 6 * columns * sizeof *g);
  }
  double phi[36];
  double* stm = g != NULL ? phi : NULL;
  myb_pv s = l->start;
  myb_pv half;
  double before = 0.0;
  for (size_t i = 0; i < l->forward; ++i) {
    const double dv = u[i] * l->impulse / (1.0 - l->rate * (before + 0.5 * u[i]));
    myb_status st = coast(&s, 0.5 * l->dt, &half, stm);
    if (st != MYB_OK) {
      return st;
    }
    for (size_t j = 0; g != NULL && j < i; ++j) {
      carry(phi, g, columns, j);
    }
    for (size_t k = 0; k < 3; ++k) {
      half.v[k] += dv * d[3 * i + k];
    }
    st = coast(&half, 0.5 * l->dt, &s, stm);
    if (st != MYB_OK) {
      return st;
    }
    if (g != NULL) {
      plant(g, columns, i);
      for (size_t j = 0; j <= i; ++j) {
        carry(phi, g, columns, j);
      }
    }
    before += u[i];
  }
  const myb_pv ahead = s;
  s = l->target;
  before = total;
  for (size_t i = n; i-- > l->forward;) {
    before -= u[i];
    const double dv = u[i] * l->impulse / (1.0 - l->rate * (before + 0.5 * u[i]));
    myb_status st = coast(&s, -0.5 * l->dt, &half, stm);
    if (st != MYB_OK) {
      return st;
    }
    for (size_t j = i + 1; g != NULL && j < n; ++j) {
      carry(phi, g, columns, j);
    }
    for (size_t k = 0; k < 3; ++k) {
      half.v[k] -= dv * d[3 * i + k];
    }
    st = coast(&half, -0.5 * l->dt, &s, stm);
    if (st != MYB_OK) {
      return st;
    }
    if (g != NULL) {
      plant(g, columns, i);
      for (size_t j = i; j < n; ++j) {
        carry(phi, g, columns, j);
      }
    }
  }
  for (size_t k = 0; k < 3; ++k) {
    c[k] = ahead.r[k] - s.r[k];
    c[3 + k] = ahead.v[k] - s.v[k];
  }
  return MYB_OK;
}

/* As meet, with the coasts in the field of a model: the sweep of leg_internal.h, the same
 * impulses. work: myb_internal_sweep_work(n) doubles, used with g only. */
static myb_status meet_field(const leg* l, myb_internal_coast* coast, double* work,
                             const double* u, const double* d, double c[6], double* g) {
  const size_t n = l->n;
  double total = 0.0;
  for (size_t i = 0; i < n; ++i) {
    total += u[i];
  }
  if (!(1.0 - l->rate * total > 0.0)) {
    return MYB_E_DOMAIN; /* more propellant than mass */
  }
  double impulses[3 * MAX_SEGMENTS];
  double before = 0.0;
  for (size_t i = 0; i < l->forward; ++i) {
    const double dv = u[i] * l->impulse / (1.0 - l->rate * (before + 0.5 * u[i]));
    for (size_t k = 0; k < 3; ++k) {
      impulses[3 * i + k] = dv * d[3 * i + k];
    }
    before += u[i];
  }
  before = total;
  for (size_t i = n; i-- > l->forward;) {
    before -= u[i];
    const double dv = u[i] * l->impulse / (1.0 - l->rate * (before + 0.5 * u[i]));
    for (size_t k = 0; k < 3; ++k) {
      impulses[3 * i + k] = dv * d[3 * i + k];
    }
  }
  myb_internal_leg half = {.n = n, .forward = l->forward, .t = 0.0, .dt = l->dt,
                           .impulses = impulses};
  memcpy(half.start, l->start.r, 3 * sizeof *half.start);
  memcpy(half.start + 3, l->start.v, 3 * sizeof *half.start);
  memcpy(half.end, l->target.r, 3 * sizeof *half.end);
  memcpy(half.end + 3, l->target.v, 3 * sizeof *half.end);
  myb_internal_sweep out;
  const myb_status st = myb_internal_leg_sweep(coast, &half, g, 0, work, &out);
  if (st == MYB_OK) {
    memcpy(c, out.c, sizeof out.c);
  }
  return st;
}

myb_status myb_sims_flanagan_mismatch(const myb_sims_flanagan* problem, const double* throttles,
                                      double mismatch[6], double* mass1, double* jacobian) {
  MYB_CHECK(problem != NULL && throttles != NULL && mismatch != NULL);
  leg l;
  myb_status st = leg_open(&l, problem);
  if (st != MYB_OK) {
    return st;
  }
  const size_t n = l.n;
  double u[MAX_SEGMENTS];
  double d[3 * MAX_SEGMENTS];
  double total = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double size = myb_vec3_norm(throttles + 3 * i);
    if (!(size <= 1.0 + 1e-9)) {
      return MYB_E_DOMAIN;
    }
    u[i] = size > 1.0 ? 1.0 : size;
    for (size_t k = 0; k < 3; ++k) {
      d[3 * i + k] = size > 0.0 ? throttles[3 * i + k] / size : 0.0;
    }
    total += u[i];
  }
  double c[6];
  if (l.field == NULL) {
    st = meet(&l, u, d, c, jacobian);
  } else {
    myb_internal_coast coast;
    const myb_allocator a = myb_internal_allocator(NULL);
    void* work = NULL;
    st = myb_internal_coast_open(&coast, l.field, l.length, l.time, l.t0, NULL);
    if (st == MYB_OK && jacobian != NULL) {
      st = myb_internal_alloc_array(&a, myb_internal_sweep_work(n), sizeof(double), &work);
    }
    if (st == MYB_OK) {
      st = meet_field(&l, &coast, work, u, d, c, jacobian);
    }
    myb_internal_release(&a, work);
    myb_internal_coast_close(&coast);
  }
  if (st != MYB_OK) {
    return st;
  }
  const double speed = l.length / l.time;
  for (size_t k = 0; k < 3; ++k) {
    mismatch[k] = c[k] * l.length;
    mismatch[3 + k] = c[3 + k] * speed;
  }
  if (mass1 != NULL) {
    *mass1 = l.mass * (1.0 - l.rate * total);
  }
  if (jacobian != NULL) {
    /* from the impulses to the throttles: dv_i = tau_i impulse / m_i, and the masses of the
     * segment and of all the later ones fall with |tau_i| */
    const size_t columns = 3 * n;
    double later[6] = {0};
    double before = total;
    for (size_t i = n; i-- > 0;) {
      before -= u[i];
      const double mid = 1.0 - l.rate * (before + 0.5 * u[i]);
      const double q = l.impulse / mid;
      for (size_t r = 0; r < 6; ++r) {
        double* row = jacobian + r * columns + 3 * i;
        const double along = row[0] * d[3 * i] + row[1] * d[3 * i + 1] + row[2] * d[3 * i + 2];
        const double by_size = along * q * u[i] * l.rate / (2.0 * mid) + later[r];
        const double unit = r < 3 ? l.length : speed;
        for (size_t k = 0; k < 3; ++k) {
          row[k] = (row[k] * q + by_size * d[3 * i + k]) * unit;
        }
        later[r] += along * u[i] * q * l.rate / mid;
      }
    }
  }
  return MYB_OK;
}

/* The program handed to myb_nlp_solve: x = (sizes, azimuths, elevations) of the throttles. */
typedef struct program {
  leg l;
  double* g;         /* 6 x 3 n */
  double* direction; /* 3 n */
  myb_internal_coast coast; /* with a field */
  double* work;             /* with a field: of the sweep */
} program;

static myb_status program_function(void* context, const double* x, double* f, double* c,
                                   double* gradient, double* jacobian) {
  program* p = context;
  const size_t n = p->l.n;
  const size_t columns = 3 * n;
  const double* u = x;
  const double* azimuth = x + n;
  const double* elevation = x + 2 * n;
  double total = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double ce = myb_cos(elevation[i]);
    p->direction[3 * i] = ce * myb_cos(azimuth[i]);
    p->direction[3 * i + 1] = ce * myb_sin(azimuth[i]);
    p->direction[3 * i + 2] = myb_sin(elevation[i]);
    total += u[i];
  }
  double* g = jacobian != NULL ? p->g : NULL;
  const myb_status st = p->l.field == NULL
                            ? meet(&p->l, u, p->direction, c, g)
                            : meet_field(&p->l, &p->coast, p->work, u, p->direction, c, g);
  if (st != MYB_OK) {
    return st;
  }
  *f = p->l.rate * total; /* the part of the mass that is spent */
  if (jacobian == NULL) {
    return MYB_OK;
  }
  for (size_t i = 0; i < columns; ++i) {
    gradient[i] = i < n ? p->l.rate : 0.0;
  }
  double later[6] = {0};
  double before = total;
  for (size_t i = n; i-- > 0;) {
    before -= u[i];
    const double mid = 1.0 - p->l.rate * (before + 0.5 * u[i]);
    const double q = p->l.impulse / mid;
    const double* d = p->direction + 3 * i;
    const double ca = myb_cos(azimuth[i]);
    const double sa = myb_sin(azimuth[i]);
    const double ce = myb_cos(elevation[i]);
    const double se = myb_sin(elevation[i]);
    const double by_azimuth[3] = {-ce * sa, ce * ca, 0.0};
    const double by_elevation[3] = {-se * ca, -se * sa, ce};
    for (size_t r = 0; r < 6; ++r) {
      const double* row = p->g + r * columns + 3 * i;
      const double along = row[0] * d[0] + row[1] * d[1] + row[2] * d[2];
      jacobian[r * columns + i] = along * q * (1.0 + u[i] * p->l.rate / (2.0 * mid)) + later[r];
      jacobian[r * columns + n + i] =
          u[i] * q * (row[0] * by_azimuth[0] + row[1] * by_azimuth[1] + row[2] * by_azimuth[2]);
      jacobian[r * columns + 2 * n + i] = u[i] * q * (row[0] * by_elevation[0] +
                                                      row[1] * by_elevation[1] +
                                                      row[2] * by_elevation[2]);
      later[r] += along * u[i] * q * p->l.rate / mid;
    }
  }
  return MYB_OK;
}

myb_status myb_sims_flanagan_solve(const myb_sims_flanagan* problem,
                                   const myb_nlp_options* options, double* throttles,
                                   const myb_allocator* allocator,
                                   myb_sims_flanagan_result* result) {
  MYB_CHECK(problem != NULL && throttles != NULL && result != NULL);
  memset(result, 0, sizeof *result);
  program p;
  myb_status st = leg_open(&p.l, problem);
  if (st != MYB_OK) {
    return st;
  }
  const size_t n = p.l.n;
  const size_t columns = 3 * n;
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  st = myb_internal_coast_open(&p.coast, p.l.field, p.l.length, p.l.time, p.l.t0, allocator);
  if (st == MYB_OK) {
    const size_t sweep = p.l.field != NULL ? myb_internal_sweep_work(n) : 0;
    st = myb_internal_alloc_array(&a, 10 * columns + sweep, sizeof(double), &block);
  }
  if (st != MYB_OK) {
    myb_internal_coast_close(&p.coast);
    return st;
  }
  p.g = block;
  p.direction = p.g + 6 * columns;
  double* x = p.direction + columns;
  double* lower = x + columns;
  double* upper = lower + columns;
  p.work = upper + columns;
  int given = 0;
  for (size_t i = 0; i < columns; ++i) {
    given = given || throttles[i] != 0.0;
    if (!isfinite(throttles[i])) {
      st = MYB_E_DOMAIN;
    }
  }
  for (size_t i = 0; i < n && st == MYB_OK; ++i) {
    /* the direction of the motion on the coast from the nearer end: where a throttle has
     * none of its own */
    myb_pv at;
    st = i < p.l.forward
             ? myb_kepler_propagate(1.0, &p.l.start, ((double)i + 0.5) * p.l.dt, &at)
             : myb_kepler_propagate(1.0, &p.l.target, -((double)(n - i) - 0.5) * p.l.dt, &at);
    if (st != MYB_OK) {
      break;
    }
    const double* tau = throttles + 3 * i;
    const double size = myb_vec3_norm(tau);
    const double speed = myb_vec3_norm(at.v);
    double d[3];
    for (size_t k = 0; k < 3; ++k) {
      d[k] = size > 0.0 ? tau[k] / size : (speed > 0.0 ? at.v[k] / speed : (k == 0 ? 1.0 : 0.0));
    }
    x[i] = given ? (size > 1.0 ? 1.0 : size) : 0.5;
    x[n + i] = myb_atan2(d[1], d[0]);
    x[2 * n + i] = myb_asin(d[2] > 1.0 ? 1.0 : (d[2] < -1.0 ? -1.0 : d[2]));
    lower[i] = 0.0;
    upper[i] = 1.0;
    lower[n + i] = lower[2 * n + i] = -HUGE_VAL;
    upper[n + i] = upper[2 * n + i] = HUGE_VAL;
  }
  if (st == MYB_OK) {
    const myb_nlp_problem nlp = {.n = columns, .equalities = 6, .lower = lower, .upper = upper,
                                 .function = program_function, .context = &p};
    const myb_nlp_options own = {.feasibility = 1e-9, .optimality = 1e-6};
    st = myb_nlp_solve(&nlp, options != NULL ? options : &own, x, NULL, allocator, &result->nlp);
  }
  if (st == MYB_OK || st == MYB_E_NOCONV) {
    double total = 0.0;
    double c[6];
    double ignored = 0.0;
    /* the directions of the solution are those of its last evaluation only if that was the
     * last point: compute them again */
    const myb_status again = program_function(&p, x, &ignored, c, NULL, NULL);
    if (again == MYB_OK) {
      for (size_t i = 0; i < n; ++i) {
        for (size_t k = 0; k < 3; ++k) {
          throttles[3 * i + k] = x[i] * p.direction[3 * i + k];
        }
        total += x[i];
      }
      const double speed = p.l.length / p.l.time;
      for (size_t k = 0; k < 3; ++k) {
        result->mismatch[k] = c[k] * p.l.length;
        result->mismatch[3 + k] = c[3 + k] * speed;
      }
      result->mass1 = p.l.mass * (1.0 - p.l.rate * total);
    } else {
      st = again;
    }
  }
  myb_internal_release(&a, block);
  myb_internal_coast_close(&p.coast);
  return st;
}
