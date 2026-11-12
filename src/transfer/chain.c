/* SPDX-License-Identifier: BSD-3-Clause */
/* Chains of legs of Sims and Flanagan joined by flybys (myb/optimal.h): one nonlinear program
 * whose unknowns are the times, the excess velocities at the bodies and the throttles of every
 * leg. The legs are swept by transfer/leg_internal.h, which gives the derivatives of each
 * mismatch at fixed impulses; here they are carried to the unknowns: the impulses grow with
 * the length of their segment and fall with the mass, the mass carries from leg to leg, and
 * the ends of a leg move with the bodies. */
#include <math.h> /* sqrt, fabs */
#include <string.h>

#include "core/alloc.h"
#include "core/check.h"
#include "myb/cmath.h"
#include "myb/kepler.h"
#include "myb/lambert.h"
#include "myb/optimal.h"
#include "num/vec3.h"
#include "transfer/leg_internal.h"

enum {
  MAX_LEGS = MYB_CHAIN_MAX_LEGS,
  MAX_SEGMENTS = MYB_SIMS_FLANAGAN_MAX_SEGMENTS,
  PER_SEGMENT = 4,  /* columns of the core for a segment: its size through the masses, then its
                       impulse at fixed masses (3) */
  FIRST_OUTER = 10, /* updates of the multipliers in the first pass */
  FIRST_INNER = 300 /* and steps between two of them */
};

/* The step of the central differences of the velocity of a body (units of the solver). */
static const double rate_step = 1e-4;

/* The size of the throttles the solver starts from when none is given. */
static const double start_size = 0.5;

/* The first weight of the constraints in both passes: the propellant should not pull the first
 * steps far from them, from a guess that nearly meets them or from a pass that did. */
static const double first_weight = 1e4;

/* A chain in the units of the solver, with its work. The solver takes the throttles in two
 * passes (myb/optimal.h tells why). In the first each is a vector tau, the impulse being q tau
 * (q = thrust dt / mass), with a size u apart from it that the masses go with, |tau|^2 <= u^2:
 * the propellant is linear in the sizes, at the least it spends u = |tau|, and the direction
 * of a throttle that is off is not lost; but the cone is degenerate where the engine is off —
 * its constraint has no gradient there — and the solver crawls. So the first pass is cut short
 * at 1e-6, and the second takes the throttles as a size and two angles from where it left
 * them. On the chain of test_chain, coplanar and with its dates fixed, the first pass alone —
 * 50 updates of the multipliers, 14 000 steps — left the constraints at 1e-8; the second,
 * after a short first, took 4 updates and 284 steps to 1e-10. */
typedef struct state {
  const myb_chain* chain;
  size_t legs;
  size_t total;                /* segments in all */
  size_t offset[MAX_LEGS + 1]; /* the first segment of each leg */
  size_t largest;              /* segments of the longest leg */
  size_t n;                    /* unknowns of the chain, and of the second pass */
  size_t m;                    /* constraints of the chain */
  size_t solver;               /* unknowns of the first pass: 4 per segment */
  size_t width;                /* columns of the core: 7 legs + 1 + 4 total */
  double length;               /* m */
  double time;                 /* s */
  double speed;                /* m/s */
  double origin;               /* s: the earliest departure, the time 0 of the solver */
  double thrust;               /* the acceleration of the first mass */
  double exhaust;
  double start_mass[MAX_LEGS + 1]; /* at each body, of the last evaluation */
  myb_internal_coast coast;
  myb_allocator allocator;
  void* block;
  double* core;      /* m x width */
  double* gradient;  /* width: of the objective */
  double* g;         /* 6 x 3 largest */
  double* sweep;     /* the work of a sweep */
  double* impulses;  /* 3 largest */
  double* masses;    /* largest: in the middle of each segment */
  double* sizes;     /* total: of the throttles, for the masses */
  double* throttles; /* 3 total: the throttles, for the impulses */
  double* turns;     /* 9 total: their directions and its derivatives by the two angles */
  double* z;         /* solver: the unknowns of the first pass */
  double* lower;     /* solver */
  double* upper;     /* solver */
  double* y;         /* n: the unknowns of the second pass */
  double* low;       /* n */
  double* high;      /* n */
} state;

/* ---- The problem ----------------------------------------------------------------------------- */

static int chain_valid(const myb_chain* c, size_t* total) {
  *total = 0;
  if (c->legs < 1 || c->legs > MAX_LEGS) {
    return 0;
  }
  double sum = c->mu + c->thrust + c->exhaust + c->mass0 + c->departure[0] + c->departure[1] +
               c->launch + c->arrival + c->epoch.f;
  int valid = c->mu > 0.0 && c->thrust > 0.0 && c->exhaust > 0.0 && c->mass0 > 0.0 &&
              c->departure[0] <= c->departure[1] && c->launch >= 0.0 && c->arrival >= 0.0;
  size_t segments = 0;
  for (size_t k = 0; k < c->legs; ++k) {
    sum += c->flight[k][0] + c->flight[k][1];
    valid = valid && c->segments[k] >= 2 && c->segments[k] <= MAX_SEGMENTS &&
            c->flight[k][0] > 0.0 && c->flight[k][0] <= c->flight[k][1];
    segments += c->segments[k] <= MAX_SEGMENTS ? c->segments[k] : 0;
  }
  for (size_t b = 0; b <= c->legs; ++b) {
    const myb_chain_body* body = &c->bodies[b];
    valid = valid && body->state != NULL;
    if (b > 0 && b < c->legs) {
      sum += body->gm + body->periapsis;
      valid = valid && body->gm > 0.0 && body->periapsis > 0.0;
    }
  }
  *total = segments;
  return valid && isfinite(sum) && 7 * c->legs + 1 + 4 * segments <= MYB_NLP_MAX_UNKNOWNS &&
         8 * c->legs + segments <= MYB_NLP_MAX_CONSTRAINTS;
}

myb_status myb_chain_size(const myb_chain* chain, size_t* unknowns, size_t* equalities,
                          size_t* inequalities) {
  MYB_CHECK(chain != NULL);
  size_t total = 0;
  if (!chain_valid(chain, &total)) {
    return MYB_E_DOMAIN;
  }
  if (unknowns != NULL) {
    *unknowns = 7 * chain->legs + 1 + 3 * total;
  }
  if (equalities != NULL) {
    *equalities = 7 * chain->legs - 1;
  }
  if (inequalities != NULL) {
    *inequalities = chain->legs + 1;
  }
  return MYB_OK;
}

/* The state of body b at the time t (s of the chain), in SI. */
static myb_status body_state(const myb_chain* chain, size_t b, double t, double y[6]) {
  const myb_chain_body* body = &chain->bodies[b];
  const myb_status st = body->state(body->context, myb_time_add(chain->epoch, t), y);
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 6; ++i) {
    if (!isfinite(y[i])) {
      return MYB_E_DOMAIN;
    }
  }
  return MYB_OK;
}

static void state_close(state* s) {
  myb_internal_coast_close(&s->coast);
  myb_internal_release(&s->allocator, s->block);
  s->block = NULL;
}

static myb_status state_open(state* s, const myb_chain* chain, const myb_allocator* allocator) {
  memset(s, 0, sizeof *s);
  s->allocator = myb_internal_allocator(allocator);
  if (!chain_valid(chain, &s->total)) {
    return MYB_E_DOMAIN;
  }
  s->chain = chain;
  s->legs = chain->legs;
  for (size_t k = 0; k < s->legs; ++k) {
    s->offset[k + 1] = s->offset[k] + chain->segments[k];
    s->largest = chain->segments[k] > s->largest ? chain->segments[k] : s->largest;
  }
  s->n = 7 * s->legs + 1 + 3 * s->total;
  s->m = 8 * s->legs;
  s->solver = 7 * s->legs + 1 + 4 * s->total;
  s->width = 7 * s->legs + 1 + PER_SEGMENT * s->total;
  /* the units: the distance of the first body at the earliest departure */
  double y[6];
  myb_status st = body_state(chain, 0, chain->departure[0], y);
  if (st != MYB_OK) {
    return st;
  }
  s->length = myb_vec3_norm(y);
  if (!(s->length > 0.0)) {
    return MYB_E_DOMAIN;
  }
  s->time = sqrt(s->length * s->length * s->length / chain->mu);
  s->speed = s->length / s->time;
  s->origin = chain->departure[0];
  s->thrust = chain->thrust * s->time * s->time / (s->length * chain->mass0);
  s->exhaust = chain->exhaust * s->time / s->length;
  st = myb_internal_coast_open(&s->coast, chain->field, s->length, s->time, s->origin, allocator);
  if (st != MYB_OK) {
    return st;
  }
  const size_t sweep = myb_internal_sweep_work(s->largest);
  const size_t doubles = s->m * s->width + s->width + 18 * s->largest + sweep + 3 * s->largest +
                         s->largest + 13 * s->total + 3 * s->solver + 3 * s->n;
  st = myb_internal_alloc_array(&s->allocator, doubles, sizeof(double), &s->block);
  if (st != MYB_OK) {
    myb_internal_coast_close(&s->coast);
    return st;
  }
  s->core = s->block;
  s->gradient = s->core + s->m * s->width;
  s->g = s->gradient + s->width;
  s->sweep = s->g + 18 * s->largest;
  s->impulses = s->sweep + sweep;
  s->masses = s->impulses + 3 * s->largest;
  s->sizes = s->masses + s->largest;
  s->throttles = s->sizes + s->total;
  s->turns = s->throttles + 3 * s->total;
  s->z = s->turns + 9 * s->total;
  s->lower = s->z + s->solver;
  s->upper = s->lower + s->solver;
  s->y = s->upper + s->solver;
  s->low = s->y + s->n;
  s->high = s->low + s->n;
  return MYB_OK;
}

/* The state of body b at the time t of the solver and, when wanted, its rate (v, a) — the
 * acceleration by central differences of the velocity. */
static myb_status body_at(const state* s, size_t b, double t, double y[6], double* rate) {
  double w[6];
  myb_status st = body_state(s->chain, b, s->origin + t * s->time, w);
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 3; ++i) {
    y[i] = w[i] / s->length;
    y[3 + i] = w[3 + i] / s->speed;
  }
  if (rate == NULL) {
    return MYB_OK;
  }
  double later[6];
  double earlier[6];
  st = body_state(s->chain, b, s->origin + (t + rate_step) * s->time, later);
  if (st == MYB_OK) {
    st = body_state(s->chain, b, s->origin + (t - rate_step) * s->time, earlier);
  }
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 3; ++i) {
    rate[i] = y[3 + i];
    rate[3 + i] = (later[3 + i] - earlier[3 + i]) / s->speed / (2.0 * rate_step);
  }
  return MYB_OK;
}

/* ---- The core: everything in the units of the solver ----------------------------------------- */

/* The columns of the core: the departure, the flights, the excess velocities (vector v: out of
 * body k at 2 k, into body k + 1 at 2 k + 1), then PER_SEGMENT per segment. */
static size_t column_excess(const state* s, size_t vector) {
  return 1 + s->legs + 3 * vector;
}

static size_t column_segment(const state* s, size_t segment) {
  return 7 * s->legs + 1 + PER_SEGMENT * segment;
}

/* What leg k does: its impulses from the throttles, its mismatch c, the part of the first mass
 * a unit of throttle spends in it (*spent) and the sum of its throttles (*sum); with
 * derivatives its rows of the core and its part of the gradient of the objective. at: the
 * times of the bodies; mass: at its start; spent, sums: of the earlier legs. */
static myb_status leg_values(state* s, size_t k, const double* at, double flight, double mass,
                             const double* bodies, const double* rates, const double* excess,
                             const double* spent, const double* sums, int derivatives, double c[6],
                             double* spent_k, double* sum_k) {
  const myb_chain* chain = s->chain;
  const size_t n = chain->segments[k];
  const size_t first = s->offset[k];
  const size_t width = s->width;
  const double* u = s->sizes + first;
  const double* tau = s->throttles + 3 * first;
  const double dt = flight / (double)n;
  const double ratio = s->thrust / s->exhaust;
  const double rho = ratio * dt;
  double total = 0.0;
  for (size_t i = 0; i < n; ++i) {
    total += u[i];
  }
  if (!(mass - rho * total > 0.0)) {
    return MYB_E_DOMAIN; /* more propellant than mass */
  }
  double before = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double m = mass - rho * (before + 0.5 * u[i]);
    const double q = s->thrust * dt / m;
    s->masses[i] = m;
    for (size_t j = 0; j < 3; ++j) {
      s->impulses[3 * i + j] = q * tau[3 * i + j];
    }
    before += u[i];
  }
  myb_internal_leg leg = {
      .n = n, .forward = (n + 1) / 2, .t = at[k], .dt = dt, .impulses = s->impulses};
  for (size_t i = 0; i < 6; ++i) {
    leg.start[i] = bodies[6 * k + i] + (i < 3 ? 0.0 : excess[6 * k + i - 3]);
    leg.end[i] = bodies[6 * (k + 1) + i] + (i < 3 ? 0.0 : excess[6 * k + i]);
  }
  myb_internal_sweep out;
  const myb_status st =
      myb_internal_leg_sweep(&s->coast, &leg, derivatives ? s->g : NULL, 1, s->sweep, &out);
  if (st != MYB_OK) {
    return st;
  }
  memcpy(c, out.c, sizeof out.c);
  *spent_k = rho;
  *sum_k = total;
  if (!derivatives) {
    return MYB_OK;
  }
  /* the impulses: at fixed masses q g (q = thrust dt / m); through the masses, the size of a
   * throttle lowers the mass of its own segment by rho / 2 and of every later one by rho; the
   * mass at the start lowers them all, and dt raises each impulse dv by dv M / (dt m) */
  const size_t columns = 3 * n;
  double* rows = s->core + 6 * k * width;
  double by_mass[6] = {0};
  double grow[6];
  double later[6] = {0};
  memcpy(grow, out.by_step, sizeof grow);
  /* bound: the segments of the leg */
  for (size_t i = n; i-- > 0;) {
    const double m = s->masses[i];
    const double q = s->thrust * dt / m;
    const double* dv = s->impulses + 3 * i;
    double* segment = rows + column_segment(s, first + i);
    for (size_t r = 0; r < 6; ++r) {
      const double* gr = s->g + r * columns + 3 * i;
      const double moved = gr[0] * dv[0] + gr[1] * dv[1] + gr[2] * dv[2]; /* by the impulse */
      by_mass[r] -= moved / m;
      grow[r] += moved * mass / (dt * m);
      segment[r * width] = moved * rho / (2.0 * m) + later[r];
      for (size_t j = 0; j < 3; ++j) {
        segment[r * width + 1 + j] = q * gr[j];
      }
      later[r] += moved * rho / m;
    }
  }
  /* the times: the leg moved as a whole with its bodies, and its segments made longer (the end
   * and its body move n times as fast) */
  const double* start_rate = rates + 6 * k;
  const double* end_rate = rates + 6 * (k + 1);
  for (size_t r = 0; r < 6; ++r) {
    double shift = out.by_shift[r];
    double end = 0.0;
    for (size_t j = 0; j < 6; ++j) {
      shift += out.by_start[6 * r + j] * start_rate[j];
      end += out.by_end[6 * r + j] * end_rate[j];
    }
    shift += end;
    double* row = rows + r * width;
    row[0] = shift;
    for (size_t j = 0; j < k; ++j) {
      /* an earlier flight: later here, and less mass at the start */
      row[1 + j] = shift - by_mass[r] * ratio * sums[j] / (double)chain->segments[j];
    }
    row[1 + k] = (grow[r] + (double)n * end) / (double)n;
    for (size_t j = 0; j < 3; ++j) {
      row[column_excess(s, 2 * k) + j] = out.by_start[6 * r + 3 + j];
      row[column_excess(s, 2 * k + 1) + j] = out.by_end[6 * r + 3 + j];
    }
    for (size_t j = 0; j < k; ++j) {
      for (size_t i = s->offset[j]; i < s->offset[j + 1]; ++i) {
        row[column_segment(s, i)] = -by_mass[r] * spent[j];
      }
    }
  }
  s->gradient[1 + k] = ratio * total / (double)n;
  for (size_t i = 0; i < n; ++i) {
    s->gradient[column_segment(s, first + i)] = rho;
  }
  return MYB_OK;
}

/* The objective (the part of the first mass spent), the constraints (equalities first) and,
 * with derivatives, the core and its gradient. times: the departure (from the earliest), then
 * the flights; excess: 2 legs vectors; the throttles in s->throttles, their sizes for the
 * masses in s->sizes. */
static myb_status core(state* s, const double* times, const double* excess, double* f, double* c,
                       int derivatives) {
  const myb_chain* chain = s->chain;
  const size_t legs = s->legs;
  const size_t width = s->width;
  if (derivatives) {
    memset(s->core, 0, s->m * width * sizeof *s->core);
    memset(s->gradient, 0, width * sizeof *s->gradient);
  }
  double at[MAX_LEGS + 1];
  at[0] = times[0];
  for (size_t k = 0; k < legs; ++k) {
    if (!(times[1 + k] > 0.0)) {
      return MYB_E_DOMAIN; /* the bodies out of order */
    }
    at[k + 1] = at[k] + times[1 + k];
  }
  double bodies[6 * (MAX_LEGS + 1)];
  double rates[6 * (MAX_LEGS + 1)];
  myb_status st = MYB_OK;
  for (size_t b = 0; b <= legs && st == MYB_OK; ++b) {
    st = body_at(s, b, at[b], bodies + 6 * b, derivatives ? rates + 6 * b : NULL);
  }
  double spent[MAX_LEGS];
  double sums[MAX_LEGS];
  double mass = 1.0;
  double objective = 0.0;
  s->start_mass[0] = 1.0;
  for (size_t k = 0; k < legs && st == MYB_OK; ++k) {
    st = leg_values(s, k, at, times[1 + k], mass, bodies, rates, excess, spent, sums, derivatives,
                    c + 6 * k, &spent[k], &sums[k]);
    objective += spent[k] * sums[k];
    mass -= spent[k] * sums[k];
    s->start_mass[k + 1] = mass;
  }
  if (st != MYB_OK) {
    return st;
  }
  *f = objective;
  /* the flybys: the excess speed kept, and the turn within the least distance; the speeds
   * themselves and not their squares, whose gradients would be as small as the speeds */
  size_t row = 6 * legs;
  for (size_t j = 1; j < legs; ++j, ++row) {
    const double* in = excess + 3 * (2 * j - 1);
    const double* out = excess + 3 * (2 * j);
    const double a = myb_vec3_norm(in);
    const double b = myb_vec3_norm(out);
    if (!(a > 0.0) || !(b > 0.0)) {
      return MYB_E_DOMAIN; /* a flyby without an excess speed has no turn */
    }
    c[row] = a - b;
    if (derivatives) {
      for (size_t i = 0; i < 3; ++i) {
        s->core[row * width + column_excess(s, 2 * j - 1) + i] = in[i] / a;
        s->core[row * width + column_excess(s, 2 * j) + i] = -out[i] / b;
      }
    }
  }
  for (size_t j = 1; j < legs; ++j, ++row) {
    const double* in = excess + 3 * (2 * j - 1);
    const double* out = excess + 3 * (2 * j);
    const double a = myb_vec3_norm(in);
    const double b = myb_vec3_norm(out);
    const myb_chain_body* body = &chain->bodies[j];
    const double kappa = body->periapsis / s->length / (body->gm / chain->mu);
    const double e = 1.0 + kappa * a * b;
    const double p = myb_vec3_dot(in, out);
    c[row] = 1.0 - 2.0 / (e * e) - p / (a * b);
    if (derivatives) {
      const double tighter = 4.0 * kappa / (e * e * e);
      for (size_t i = 0; i < 3; ++i) {
        s->core[row * width + column_excess(s, 2 * j - 1) + i] =
            tighter * b * in[i] / a - out[i] / (a * b) + p * in[i] / (a * a * a * b);
        s->core[row * width + column_excess(s, 2 * j) + i] =
            tighter * a * out[i] / b - in[i] / (a * b) + p * out[i] / (a * b * b * b);
      }
    }
  }
  /* the excess speeds at the two ends (no derivative at a speed of zero, where a limit of zero
   * holds the excess velocity) */
  const double* ends[2] = {excess, excess + 3 * (2 * legs - 1)};
  const double limits[2] = {chain->launch / s->speed, chain->arrival / s->speed};
  const size_t vectors[2] = {0, 2 * legs - 1};
  for (size_t e = 0; e < 2; ++e, ++row) {
    const double speed = myb_vec3_norm(ends[e]);
    c[row] = speed - limits[e];
    for (size_t i = 0; i < 3 && derivatives && speed > 0.0; ++i) {
      s->core[row * width + column_excess(s, vectors[e]) + i] = ends[e][i] / speed;
    }
  }
  return MYB_OK;
}

/* The scale of constraint r in SI: the mismatches in m and m/s, the speeds in m/s, the turns
 * without a unit. */
static double row_unit(const state* s, size_t r) {
  const size_t legs = s->legs;
  if (r < 6 * legs) {
    return r % 6 < 3 ? s->length : s->speed;
  }
  if (r < 7 * legs - 1) {
    return s->speed;
  }
  if (r < 8 * legs - 2) {
    return 1.0;
  }
  return s->speed;
}

/* The times and the excess velocities of the unknowns x (SI) in those of the solver (z). */
static void scale_in(const state* s, const double* x, double* z) {
  z[0] = (x[0] - s->origin) / s->time;
  for (size_t k = 0; k < s->legs; ++k) {
    z[1 + k] = x[1 + k] / s->time;
  }
  for (size_t i = 0; i < 6 * s->legs; ++i) {
    z[1 + s->legs + i] = x[1 + s->legs + i] / s->speed;
  }
}

/* The throttles of the chain (3 per segment): each its own size, one longer than 1 cut to 1;
 * MYB_E_DOMAIN for one longer than 1 + 1e-9. */
static myb_status take_throttles(state* s, const double* throttles) {
  for (size_t i = 0; i < s->total; ++i) {
    const double size = myb_vec3_norm(throttles + 3 * i);
    if (!(size <= 1.0 + 1e-9)) {
      return MYB_E_DOMAIN;
    }
    s->sizes[i] = size > 1.0 ? 1.0 : size;
    for (size_t j = 0; j < 3; ++j) {
      s->throttles[3 * i + j] = size > 1.0 ? throttles[3 * i + j] / size : throttles[3 * i + j];
    }
  }
  return MYB_OK;
}

myb_status myb_chain_evaluate(const myb_chain* chain, const double* x, double* propellant,
                              double* constraints, double* gradient, double* jacobian,
                              const myb_allocator* allocator) {
  MYB_CHECK(chain != NULL && x != NULL);
  state s;
  myb_status st = state_open(&s, chain, allocator);
  if (st != MYB_OK) {
    return st;
  }
  const size_t base = 7 * s.legs + 1;
  const int derivatives = gradient != NULL || jacobian != NULL;
  double f = 0.0;
  double* c = s.lower; /* m <= n values of work */
  for (size_t i = 0; i < s.n && st == MYB_OK; ++i) {
    if (!isfinite(x[i])) {
      st = MYB_E_DOMAIN;
    }
  }
  if (st == MYB_OK) {
    scale_in(&s, x, s.z);
    st = take_throttles(&s, x + base);
  }
  if (st == MYB_OK) {
    st = core(&s, s.z, s.z + 1 + s.legs, &f, c, derivatives);
  }
  if (st == MYB_OK) {
    if (propellant != NULL) {
      *propellant = f * chain->mass0;
    }
    for (size_t r = 0; r < s.m && constraints != NULL; ++r) {
      constraints[r] = c[r] * row_unit(&s, r);
    }
    /* by the unknowns in their units: the times per second, the excess velocities per m/s, a
     * throttle tau at a fixed mass and through its size |tau| */
    for (size_t r = 0; r <= s.m && derivatives; ++r) {
      const int objective = r == s.m;
      if ((objective && gradient == NULL) || (!objective && jacobian == NULL)) {
        continue;
      }
      const double* from = objective ? s.gradient : s.core + r * s.width;
      double* to = objective ? gradient : jacobian + r * s.n;
      const double unit = objective ? chain->mass0 : row_unit(&s, r);
      for (size_t j = 0; j < 1 + s.legs; ++j) {
        to[j] = from[j] * unit / s.time;
      }
      for (size_t j = 1 + s.legs; j < base; ++j) {
        to[j] = from[j] * unit / s.speed;
      }
      for (size_t i = 0; i < s.total; ++i) {
        const double* segment = from + column_segment(&s, i);
        const double* tau = s.throttles + 3 * i;
        const double size = s.sizes[i];
        for (size_t j = 0; j < 3; ++j) {
          const double along = size > 0.0 ? tau[j] / size : 0.0;
          to[base + 3 * i + j] = (segment[1 + j] + segment[0] * along) * unit;
        }
      }
    }
  }
  state_close(&s);
  return st;
}

/* ---- The first guess ------------------------------------------------------------------------- */

myb_status myb_chain_guess(const myb_chain* chain, double* x) {
  MYB_CHECK(chain != NULL && x != NULL);
  size_t unknowns = 0;
  myb_status st = myb_chain_size(chain, &unknowns, NULL, NULL);
  if (st != MYB_OK) {
    return st;
  }
  const size_t legs = chain->legs;
  memset(x, 0, unknowns * sizeof *x);
  x[0] = 0.5 * (chain->departure[0] + chain->departure[1]);
  double at = x[0];
  double here[6];
  st = body_state(chain, 0, at, here);
  if (st != MYB_OK) {
    return st;
  }
  /* the sense of the motion of the first body */
  const int retrograde = here[0] * here[4] - here[1] * here[3] < 0.0;
  for (size_t k = 0; k < legs; ++k) {
    x[1 + k] = 0.5 * (chain->flight[k][0] + chain->flight[k][1]);
    double there[6];
    st = body_state(chain, k + 1, at + x[1 + k], there);
    if (st != MYB_OK) {
      return st;
    }
    myb_lambert_solution arc;
    size_t count = 0;
    st = myb_lambert(chain->mu, here, there, x[1 + k], retrograde, 0, &arc, 1, &count);
    if (st != MYB_OK) {
      return st;
    }
    if (count == 0) {
      return MYB_E_NOCONV;
    }
    double* out = x + 1 + legs + 6 * k;
    for (size_t i = 0; i < 3; ++i) {
      out[i] = arc.v1[i] - here[3 + i];
      out[3 + i] = arc.v2[i] - there[3 + i];
    }
    at += x[1 + k];
    memcpy(here, there, sizeof here);
  }
  /* cut down to what the launcher gives and what may be left at the arrival */
  double* ends[2] = {x + 1 + legs, x + 1 + legs + 3 * (2 * legs - 1)};
  const double limits[2] = {chain->launch, chain->arrival};
  for (size_t e = 0; e < 2; ++e) {
    const double speed = myb_vec3_norm(ends[e]);
    if (speed > limits[e]) {
      for (size_t i = 0; i < 3; ++i) {
        ends[e][i] *= limits[e] / speed;
      }
    }
  }
  return MYB_OK;
}

/* ---- The program ----------------------------------------------------------------------------- */

/* The unknowns of the solver: the departure, the flights, the excess velocities (all in its
 * units), then the throttles tau (3 per segment) and their sizes u; its constraints those of
 * the chain, then |tau|^2 - u^2 <= 0 for each segment. */
static myb_status program_function(void* context, const double* z, double* f, double* c,
                                   double* gradient, double* jacobian) {
  state* s = context;
  const size_t base = 7 * s->legs + 1;
  const size_t total = s->total;
  const size_t n = s->solver;
  const double* tau = z + base;
  const double* u = z + base + 3 * total;
  memcpy(s->throttles, tau, 3 * total * sizeof *tau);
  memcpy(s->sizes, u, total * sizeof *u);
  const myb_status st = core(s, z, z + 1 + s->legs, f, c, jacobian != NULL);
  if (st != MYB_OK) {
    return st;
  }
  for (size_t i = 0; i < total; ++i) {
    c[s->m + i] = myb_vec3_dot(tau + 3 * i, tau + 3 * i) - u[i] * u[i];
  }
  if (jacobian == NULL) {
    return MYB_OK;
  }
  for (size_t r = 0; r <= s->m; ++r) {
    const int objective = r == s->m;
    const double* from = objective ? s->gradient : s->core + r * s->width;
    double* to = objective ? gradient : jacobian + r * n;
    memcpy(to, from, base * sizeof *to);
    for (size_t i = 0; i < total; ++i) {
      const double* segment = from + column_segment(s, i);
      for (size_t j = 0; j < 3; ++j) {
        to[base + 3 * i + j] = segment[1 + j];
      }
      to[base + 3 * total + i] = segment[0];
    }
  }
  for (size_t i = 0; i < total; ++i) {
    double* row = jacobian + (s->m + i) * n;
    memset(row, 0, n * sizeof *row);
    for (size_t j = 0; j < 3; ++j) {
      row[base + 3 * i + j] = 2.0 * tau[3 * i + j];
    }
    row[base + 3 * total + i] = -2.0 * u[i];
  }
  return MYB_OK;
}

/* The second pass: the departure, the flights, the excess velocities, then the sizes, the
 * azimuths and the elevations of the throttles; its constraints those of the chain. */
static myb_status program_angles(void* context, const double* z, double* f, double* c,
                                 double* gradient, double* jacobian) {
  state* s = context;
  const size_t base = 7 * s->legs + 1;
  const size_t total = s->total;
  const size_t n = s->n;
  for (size_t i = 0; i < total; ++i) {
    const double u = z[base + i];
    const double ca = myb_cos(z[base + total + i]);
    const double sa = myb_sin(z[base + total + i]);
    const double ce = myb_cos(z[base + 2 * total + i]);
    const double se = myb_sin(z[base + 2 * total + i]);
    double* turn = s->turns + 9 * i;
    turn[0] = ce * ca; /* the direction */
    turn[1] = ce * sa;
    turn[2] = se;
    turn[3] = -ce * sa; /* by the azimuth */
    turn[4] = ce * ca;
    turn[5] = 0.0;
    turn[6] = -se * ca; /* by the elevation */
    turn[7] = -se * sa;
    turn[8] = ce;
    s->sizes[i] = u;
    for (size_t j = 0; j < 3; ++j) {
      s->throttles[3 * i + j] = u * turn[j];
    }
  }
  const myb_status st = core(s, z, z + 1 + s->legs, f, c, jacobian != NULL);
  if (st != MYB_OK || jacobian == NULL) {
    return st;
  }
  for (size_t r = 0; r <= s->m; ++r) {
    const int objective = r == s->m;
    const double* from = objective ? s->gradient : s->core + r * s->width;
    double* to = objective ? gradient : jacobian + r * n;
    memcpy(to, from, base * sizeof *to);
    for (size_t i = 0; i < total; ++i) {
      const double* segment = from + column_segment(s, i);
      const double* turn = s->turns + 9 * i;
      const double u = s->sizes[i];
      to[base + i] =
          segment[0] + segment[1] * turn[0] + segment[2] * turn[1] + segment[3] * turn[2];
      to[base + total + i] =
          u * (segment[1] * turn[3] + segment[2] * turn[4] + segment[3] * turn[5]);
      to[base + 2 * total + i] =
          u * (segment[1] * turn[6] + segment[2] * turn[7] + segment[3] * turn[8]);
    }
  }
  return MYB_OK;
}

/* The direction of a throttle that has none: that of the motion in the middle of its segment
 * on the coast from the nearer end of its leg. */
static myb_status along_motion(const state* s, size_t k, size_t i, const double* z, double d[3]) {
  const myb_chain* chain = s->chain;
  const size_t n = chain->segments[k];
  double at = z[0];
  for (size_t j = 0; j < k; ++j) {
    at += z[1 + j];
  }
  const double dt = z[1 + k] / (double)n;
  const int forward = i < (n + 1) / 2;
  double y[6];
  myb_status st = body_at(s, forward ? k : k + 1, forward ? at : at + z[1 + k], y, NULL);
  if (st != MYB_OK) {
    return st;
  }
  const double* excess = z + 1 + s->legs + 3 * (forward ? 2 * k : 2 * k + 1);
  myb_pv from;
  myb_pv to;
  for (size_t j = 0; j < 3; ++j) {
    from.r[j] = y[j];
    from.v[j] = y[3 + j] + excess[j];
  }
  const double dt_coast = forward ? ((double)i + 0.5) * dt : -((double)(n - i) - 0.5) * dt;
  st = myb_kepler_propagate(1.0, &from, dt_coast, &to);
  if (st != MYB_OK) {
    return st;
  }
  const double speed = myb_vec3_norm(to.v);
  for (size_t j = 0; j < 3; ++j) {
    d[j] = speed > 0.0 ? to.v[j] / speed : (j == 0 ? 1.0 : 0.0);
  }
  return MYB_OK;
}

/* The bounds of the times and of the excess velocities, the same in both passes: the times
 * between theirs, the excess velocity at an end held at zero where nothing is allowed. */
static void bounds_of_times(const state* s, double* lower, double* upper) {
  const myb_chain* chain = s->chain;
  const size_t legs = s->legs;
  const size_t base = 7 * legs + 1;
  lower[0] = 0.0;
  upper[0] = (chain->departure[1] - chain->departure[0]) / s->time;
  for (size_t k = 0; k < legs; ++k) {
    lower[1 + k] = chain->flight[k][0] / s->time;
    upper[1 + k] = chain->flight[k][1] / s->time;
  }
  for (size_t i = 1 + legs; i < base; ++i) {
    lower[i] = -HUGE_VAL;
    upper[i] = HUGE_VAL;
  }
  for (size_t i = 0; i < 3; ++i) {
    if (chain->launch == 0.0) {
      lower[1 + legs + i] = upper[1 + legs + i] = 0.0;
    }
    if (chain->arrival == 0.0) {
      lower[base - 3 + i] = upper[base - 3 + i] = 0.0;
    }
  }
}

myb_status myb_chain_solve(const myb_chain* chain, const myb_nlp_options* options, double* x,
                           const myb_allocator* allocator, myb_chain_result* result) {
  MYB_CHECK(chain != NULL && x != NULL && result != NULL);
  memset(result, 0, sizeof *result);
  state s;
  myb_status st = state_open(&s, chain, allocator);
  if (st != MYB_OK) {
    return st;
  }
  const size_t legs = s.legs;
  const size_t base = 7 * legs + 1;
  const size_t total = s.total;
  int given = 0;
  for (size_t i = 0; i < s.n; ++i) {
    given = given || x[i] != 0.0;
    if (!isfinite(x[i])) {
      st = MYB_E_DOMAIN;
    }
  }
  if (st == MYB_OK && !given) {
    st = myb_chain_guess(chain, x);
  }
  if (st == MYB_OK) {
    scale_in(&s, x, s.z);
    bounds_of_times(&s, s.lower, s.upper);
  }
  /* the throttles: those given, each its own size; none given, half thrust along the motion,
   * as Sims and Flanagan alone start */
  int throttled = 0;
  for (size_t i = base; i < s.n; ++i) {
    throttled = throttled || x[i] != 0.0;
  }
  for (size_t k = 0; k < legs && st == MYB_OK; ++k) {
    for (size_t i = 0; i < chain->segments[k] && st == MYB_OK; ++i) {
      const size_t at = s.offset[k] + i;
      const double* tau = x + base + 3 * at;
      const double size = myb_vec3_norm(tau);
      double* to = s.z + base + 3 * at;
      if (throttled) {
        for (size_t j = 0; j < 3; ++j) {
          to[j] = size > 1.0 ? tau[j] / size : tau[j];
        }
        s.z[base + 3 * total + at] = size > 1.0 ? 1.0 : size;
      } else {
        st = along_motion(&s, k, i, s.z, to);
        for (size_t j = 0; j < 3; ++j) {
          to[j] *= start_size;
        }
        s.z[base + 3 * total + at] = start_size;
      }
      for (size_t j = 0; j < 3; ++j) {
        s.lower[base + 3 * at + j] = -1.0;
        s.upper[base + 3 * at + j] = 1.0;
      }
      s.lower[base + 3 * total + at] = 0.0;
      s.upper[base + 3 * total + at] = 1.0;
    }
  }
  const myb_nlp_options asked = options != NULL ? *options : (myb_nlp_options){0};
  /* the first pass: the throttles as vectors, cut short */
  if (st == MYB_OK) {
    const myb_nlp_problem first = {.n = s.solver,
                                   .equalities = 7 * legs - 1,
                                   .inequalities = legs + 1 + total,
                                   .lower = s.lower,
                                   .upper = s.upper,
                                   .function = program_function,
                                   .context = &s};
    const myb_nlp_options rough = {.feasibility = 1e-6,
                                   .optimality = 1e-4,
                                   .max_outer = FIRST_OUTER,
                                   .max_inner = FIRST_INNER,
                                   .weight = asked.weight > 0.0 ? asked.weight : first_weight};
    st = myb_nlp_solve(&first, &rough, s.z, NULL, allocator, &result->first);
    if (st == MYB_E_NOCONV) {
      st = MYB_OK; /* where it got to is all the second pass needs */
    }
  }
  /* the second: a size and two angles each, from there */
  if (st == MYB_OK) {
    memcpy(s.y, s.z, base * sizeof *s.y);
    memcpy(s.low, s.lower, base * sizeof *s.low);
    memcpy(s.high, s.upper, base * sizeof *s.high);
    for (size_t i = 0; i < total; ++i) {
      const double* tau = s.z + base + 3 * i;
      const double size = myb_vec3_norm(tau);
      const double up = size > 0.0 ? tau[2] / size : 0.0;
      s.y[base + i] = s.z[base + 3 * total + i];
      s.y[base + total + i] = myb_atan2(tau[1], tau[0]);
      s.y[base + 2 * total + i] = myb_asin(up > 1.0 ? 1.0 : (up < -1.0 ? -1.0 : up));
      s.low[base + i] = 0.0;
      s.high[base + i] = 1.0;
      s.low[base + total + i] = s.low[base + 2 * total + i] = -HUGE_VAL;
      s.high[base + total + i] = s.high[base + 2 * total + i] = HUGE_VAL;
    }
    const myb_nlp_problem second = {.n = s.n,
                                    .equalities = 7 * legs - 1,
                                    .inequalities = legs + 1,
                                    .lower = s.low,
                                    .upper = s.high,
                                    .function = program_angles,
                                    .context = &s};
    myb_nlp_options fine = {.feasibility = 1e-9, .optimality = 1e-6, .weight = first_weight};
    if (options != NULL) {
      fine = asked;
      fine.weight = asked.weight > 0.0 ? asked.weight : first_weight;
    }
    st = myb_nlp_solve(&second, &fine, s.y, NULL, allocator, &result->nlp);
  }
  if (st == MYB_OK || st == MYB_E_NOCONV) {
    /* back to the unknowns of the chain; what the solution looks like is evaluated as
     * myb_chain_evaluate has it */
    const double* y = s.y;
    double f = 0.0;
    double* c = s.core; /* m values of work */
    myb_status again = program_angles(&s, y, &f, c, NULL, NULL);
    if (again == MYB_OK) {
      x[0] = s.origin + y[0] * s.time;
      for (size_t k = 0; k < legs; ++k) {
        x[1 + k] = y[1 + k] * s.time;
      }
      for (size_t i = 1 + legs; i < base; ++i) {
        x[i] = y[i] * s.speed;
      }
      memcpy(x + base, s.throttles, 3 * total * sizeof *x);
      again = take_throttles(&s, x + base);
    }
    if (again == MYB_OK) {
      again = core(&s, y, y + 1 + legs, &f, c, 0);
    }
    if (again == MYB_OK) {
      for (size_t b = 0; b <= legs; ++b) {
        result->mass[b] = s.start_mass[b] * chain->mass0;
        const double* in = b > 0 ? y + 1 + legs + 3 * (2 * b - 1) : NULL;
        const double* out = b < legs ? y + 1 + legs + 3 * (2 * b) : NULL;
        result->excess[b][0] = in != NULL ? myb_vec3_norm(in) * s.speed : 0.0;
        result->excess[b][1] = out != NULL ? myb_vec3_norm(out) * s.speed : 0.0;
        if (in != NULL && out != NULL) {
          double normal[3];
          myb_vec3_cross(in, out, normal);
          const double turn = myb_atan2(myb_vec3_norm(normal), myb_vec3_dot(in, out));
          const double v2 = result->excess[b][0] * result->excess[b][1];
          const double half = myb_sin(0.5 * turn);
          result->turn[b] = turn;
          result->periapsis[b] =
              half > 0.0 && v2 > 0.0 ? chain->bodies[b].gm / v2 * (1.0 / half - 1.0) : HUGE_VAL;
        }
      }
      for (size_t k = 0; k < legs; ++k) {
        for (size_t r = 0; r < 6; ++r) {
          result->mismatch[k][r] = c[6 * k + r] * row_unit(&s, 6 * k + r);
        }
      }
    } else {
      st = again;
    }
  }
  state_close(&s);
  return st;
}
