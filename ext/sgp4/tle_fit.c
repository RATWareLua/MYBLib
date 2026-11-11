/* SPDX-License-Identifier: BSD-3-Clause */
/* Fitting an element set of SGP4 to positions (myb_sgp4.h): the mean elements of the theory
 * that carry its orbit through an ephemeris, by Gauss-Newton steps on the residuals of the
 * positions with the derivatives taken by central differences of the theory itself.
 *
 * The unknowns are free of the singularities of the classical angles (a circular or an
 * equatorial orbit):
 *   n                              mean motion (rev/day),
 *   af = e cos(w + node), ag = e sin(w + node),
 *   chi = tan(i/2) cos(node), psi = tan(i/2) sin(node),
 *   L = M + w + node               (rad),
 * and, if asked for, B*. */
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "core/alloc.h"
#include "core/check.h"
#include "myb/cmath.h"
#include "myb/kepler.h"
#include "myb/num.h"
#include "myb_sgp4.h"

enum { MAX_UNKNOWNS = 7, HALVINGS = 10 };

static const double pi = 3.14159265358979323846264338327950288;

static double degrees_in_turn(double radians) {
  double d = radians * 180.0 / pi;
  d -= 360.0 * floor(d / 360.0);
  return d >= 360.0 ? 0.0 : d;
}

/* The unknowns of an element set. */
static void unknowns_of(const myb_tle* tle, double x[MAX_UNKNOWNS]) {
  const double rad = pi / 180.0;
  const double node = tle->node * rad;
  const double along = (tle->argument_of_perigee + tle->node) * rad;
  const double half = myb_tan(0.5 * tle->inclination * rad);
  x[0] = tle->mean_motion;
  x[1] = tle->eccentricity * myb_cos(along);
  x[2] = tle->eccentricity * myb_sin(along);
  x[3] = half * myb_cos(node);
  x[4] = half * myb_sin(node);
  x[5] = (tle->mean_anomaly + tle->argument_of_perigee + tle->node) * rad;
  x[6] = tle->bstar;
}

/* The element set of the unknowns (the rest of *tle is kept). 0 if they are no orbit. */
static int tle_of(const double x[MAX_UNKNOWNS], myb_tle* tle) {
  const double e = sqrt(x[1] * x[1] + x[2] * x[2]);
  const double half = sqrt(x[3] * x[3] + x[4] * x[4]);
  if (!(x[0] > 0.0) || !(e < 1.0) || !isfinite(x[0] + e + half + x[5] + x[6])) {
    return 0;
  }
  const double along = e > 0.0 ? myb_atan2(x[2], x[1]) : 0.0;
  const double node = half > 0.0 ? myb_atan2(x[4], x[3]) : 0.0;
  tle->mean_motion = x[0];
  tle->eccentricity = e;
  tle->inclination = 2.0 * myb_atan(half) * 180.0 / pi;
  tle->node = degrees_in_turn(node);
  tle->argument_of_perigee = degrees_in_turn(along - node);
  tle->mean_anomaly = degrees_in_turn(x[5] - along);
  tle->bstar = x[6];
  return 1;
}

typedef struct fit {
  const double* minutes;
  const double* positions;
  size_t count;
  myb_sgp4_gravity gravity;
  myb_sgp4_mode mode;
  myb_tle shell; /* what the unknowns do not touch */
} fit;

/* The residuals (computed less given, 3 per time) of the unknowns x; their sum of squares. */
static myb_status residuals(const fit* f, const double x[MAX_UNKNOWNS], double* r,
                            double* squares) {
  myb_tle tle = f->shell;
  if (!tle_of(x, &tle)) {
    return MYB_E_DOMAIN;
  }
  myb_sgp4 theory;
  myb_status st = myb_sgp4_init(&tle, f->gravity, f->mode, &theory);
  if (st != MYB_OK) {
    return st;
  }
  double sum = 0.0;
  for (size_t k = 0; k < f->count; ++k) {
    double p[3];
    double v[3];
    st = myb_sgp4_propagate(&theory, f->minutes[k], p, v, NULL);
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 3; ++i) {
      const double d = p[i] - f->positions[3 * k + (size_t)i];
      r[3 * k + (size_t)i] = d;
      sum += d * d;
    }
  }
  *squares = sum;
  return MYB_OK;
}

myb_status myb_tle_fit(const double* minutes, const double* positions, size_t count,
                       const myb_tle_fit_options* options, myb_tle* tle, double* rms,
                       int* iterations) {
  MYB_CHECK(minutes != NULL && positions != NULL && tle != NULL);
  const myb_tle_fit_options none = {.gravity = MYB_SGP4_WGS72};
  const myb_tle_fit_options* o = options != NULL ? options : &none;
  const size_t unknowns = o->bstar ? 7 : 6;
  if (count < 3 || !(o->tolerance >= 0.0) || o->max_iterations < 0 || count > SIZE_MAX / 64) {
    return MYB_E_DOMAIN;
  }
  const double tolerance = o->tolerance > 0.0 ? o->tolerance : 1e-9;
  const int budget = o->max_iterations > 0 ? o->max_iterations : 40;
  const size_t m = 3 * count;
  const fit f = {.minutes = minutes, .positions = positions, .count = count,
                 .gravity = o->gravity, .mode = o->mode, .shell = *tle};
  const myb_allocator a = myb_internal_allocator(NULL);
  void* block = NULL;
  /* the Jacobian, the residuals, two more columns of residuals for the differences */
  myb_status st = myb_internal_alloc_array(&a, m * (unknowns + 4), sizeof(double), &block);
  if (st != MYB_OK) {
    return st;
  }
  double* jacobian = block;
  double* r = jacobian + m * unknowns;
  double* plus = r + m;
  double* minus = plus + m;
  double* rhs = minus + m;
  double x[MAX_UNKNOWNS];
  unknowns_of(tle, x);
  double squares = 0.0;
  st = residuals(&f, x, r, &squares);
  int pass = 0;
  int settled = 0;
  /* bound: the budget of iterations */
  for (pass = 0; pass < budget && st == MYB_OK && !settled; ++pass) {
    for (size_t j = 0; j < unknowns && st == MYB_OK; ++j) {
      /* steps that the theory resolves: a part in 10^7 of the mean motion, 10^-7 of the
       * others, a thousandth of B* or of its usual size */
      const double h = j == 0 ? 1e-7 * x[0] : (j == 6 ? 1e-3 * fmax(fabs(x[6]), 1e-5) : 1e-7);
      double xp[MAX_UNKNOWNS];
      double xm[MAX_UNKNOWNS];
      memcpy(xp, x, sizeof xp);
      memcpy(xm, x, sizeof xm);
      xp[j] += h;
      xm[j] -= h;
      double unused = 0.0;
      st = residuals(&f, xp, plus, &unused);
      if (st == MYB_OK) {
        st = residuals(&f, xm, minus, &unused);
      }
      for (size_t i = 0; i < m && st == MYB_OK; ++i) {
        jacobian[i * unknowns + j] = (plus[i] - minus[i]) / (xp[j] - xm[j]);
      }
    }
    if (st != MYB_OK) {
      break;
    }
    for (size_t i = 0; i < m; ++i) {
      rhs[i] = -r[i];
    }
    st = myb_least_squares(m, unknowns, jacobian, rhs, NULL);
    if (st != MYB_OK) {
      break;
    }
    /* the step, halved while the residuals do not fall */
    double scale = 1.0;
    double trial[MAX_UNKNOWNS];
    double trial_squares = squares;
    int improved = 0;
    for (int halving = 0; halving < HALVINGS && !improved; ++halving) {
      memcpy(trial, x, sizeof trial);
      for (size_t j = 0; j < unknowns; ++j) {
        trial[j] += scale * rhs[j];
      }
      const myb_status tried = residuals(&f, trial, plus, &trial_squares);
      improved = tried == MYB_OK && trial_squares < squares;
      scale *= 0.5;
    }
    if (!improved) {
      settled = 1; /* no step lowers the residuals: at the minimum the differences resolve */
      break;
    }
    const double before = sqrt(squares / (double)count);
    const double after = sqrt(trial_squares / (double)count);
    memcpy(x, trial, sizeof x);
    memcpy(r, plus, m * sizeof *r);
    squares = trial_squares;
    settled = before - after <= tolerance * before;
  }
  myb_internal_release(&a, block);
  if (st != MYB_OK) {
    return st;
  }
  if (!tle_of(x, tle)) {
    return MYB_E_DOMAIN;
  }
  if (rms != NULL) {
    *rms = sqrt(squares / (double)count);
  }
  if (iterations != NULL) {
    *iterations = pass;
  }
  return settled ? MYB_OK : MYB_E_NOCONV;
}

myb_status myb_tle_from_state(const double r[3], const double v[3], myb_sgp4_gravity gravity,
                              myb_tle* tle) {
  MYB_CHECK(r != NULL && v != NULL && tle != NULL);
  myb_sgp4_constants c;
  myb_sgp4_constants_of(gravity, &c);
  /* GM in m^3/s^2 from xke (Earth radii^1.5 per minute) */
  const double radius = c.radius * 1000.0;
  const double mu = c.xke * c.xke * radius * radius * radius / 3600.0;
  myb_pv pv;
  memcpy(pv.r, r, sizeof pv.r);
  memcpy(pv.v, v, sizeof pv.v);
  myb_elements el;
  myb_status st = myb_elements_from_pv(mu, &pv, &el);
  if (st != MYB_OK) {
    return st;
  }
  if (!(el.e < 1.0)) {
    return MYB_E_DOMAIN;
  }
  double eccentric = 0.0;
  double mean = 0.0;
  st = myb_anomaly_eccentric_from_true(el.nu, el.e, &eccentric);
  if (st == MYB_OK) {
    st = myb_anomaly_mean_from_eccentric(eccentric, el.e, &mean);
  }
  if (st != MYB_OK) {
    return st;
  }
  const double a = el.p / (1.0 - el.e * el.e);
  tle->mean_motion = sqrt(mu / (a * a * a)) * 86400.0 / (2.0 * pi);
  tle->eccentricity = el.e;
  tle->inclination = el.i * 180.0 / pi;
  tle->node = degrees_in_turn(el.raan);
  tle->argument_of_perigee = degrees_in_turn(el.argp);
  tle->mean_anomaly = degrees_in_turn(mean);
  return MYB_OK;
}
