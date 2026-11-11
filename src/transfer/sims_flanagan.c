/* SPDX-License-Identifier: BSD-3-Clause */
/* The direct method of Sims and Flanagan (myb/optimal.h): impulses in the middle of equal
 * segments, conics between them, the two halves of the transfer meeting in the middle. */
#include <math.h> /* sqrt, fabs */
#include <string.h>

#include "core/alloc.h"
#include "core/check.h"
#include "myb/cmath.h"
#include "myb/kepler.h"
#include "myb/optimal.h"
#include "num/vec3.h"

enum { MAX_SEGMENTS = MYB_SIMS_FLANAGAN_MAX_SEGMENTS };

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
  st = meet(&l, u, d, c, jacobian);
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
  const myb_status st = meet(&p->l, u, p->direction, c, jacobian != NULL ? p->g : NULL);
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
  st = myb_internal_alloc_array(&a, 10 * columns, sizeof(double), &block);
  if (st != MYB_OK) {
    return st;
  }
  p.g = block;
  p.direction = p.g + 6 * columns;
  double* x = p.direction + columns;
  double* lower = x + columns;
  double* upper = lower + columns;
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
  return st;
}
