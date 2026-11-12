/* SPDX-License-Identifier: BSD-3-Clause */
/* Random impacts of micrometeoroids (myb/meteoroid.h): the flux of Grun et al. (1985), the
 * Poisson process of the impacts, and the impulses applied to an orbit. */
#include "myb/meteoroid.h"

#include <math.h> /* sqrt, ceil */
#include <string.h>

#include "core/base.h"
#include "myb/cmath.h"
#include "num/root.h"
#include "num/vec3.h"

enum { MASS_ITERATIONS = 200 }; /* bracketed Newton: more than the bits of a double */

static const double AU = 149597870700.0;          /* m, IAU 2012 B2 */
static const double LN_GRAMS = 6.907755278982137; /* ln 1000: kg -> g */
static const double DEFAULT_SPEED = 2.0e4;        /* m/s */
static const double TWO_PI = 6.283185307179586;

/* ---- The flux ---------------------------------------------------------------------------------
 * F of Grun et al. (1985, eq. A3) as a function of s = ln m, m in grams, with dF/ds. A term
 * that underflows (a mass far beyond the range of the fit) is left out of the derivative,
 * where it would be 0 times infinity. */
static double grun(double s, double* slope) {
  const double m = myb_exp(s);
  const double m2 = m * m;
  const double small = myb_exp(0.306 * s); /* m^g4 */
  const double g1 = 2.2e3 * small + 15.0;
  const double t1 = myb_exp(-4.38 * myb_log(g1));
  const double g2 = m + 1.0e11 * m2 + 1.0e27 * m2 * m2;
  const double t2 = 1.3e-9 * myb_exp(-0.36 * myb_log(g2));
  const double g3 = m + 1.0e6 * m2;
  const double t3 = 1.3e-16 * myb_exp(-0.85 * myb_log(g3));
  if (slope != NULL) {
    double d = 0.0;
    if (t1 > 0.0) {
      d += -4.38 * t1 / g1 * (2.2e3 * 0.306 * small);
    }
    if (t2 > 0.0) {
      d += -0.36 * t2 / g2 * (m + 2.0e11 * m2 + 4.0e27 * m2 * m2);
    }
    if (t3 > 0.0) {
      d += -0.85 * t3 / g3 * (m + 2.0e6 * m2);
    }
    *slope = d;
  }
  return t1 + t2 + t3;
}

static int finite_positive(double x) {
  return x > 0.0 && x - x == 0.0;
}

double myb_meteoroid_flux(double mass) {
  if (!finite_positive(mass)) {
    return (double)NAN;
  }
  return grun(myb_log(mass) + LN_GRAMS, NULL);
}

static int valid_options(const myb_meteoroid_options* o) {
  return finite_positive(o->cross_section) && finite_positive(o->mass_min) &&
         finite_positive(o->mass_max) && o->mass_max > o->mass_min && o->distance >= 0.0 &&
         o->distance - o->distance == 0.0 && o->exponent - o->exponent == 0.0 && o->speed >= 0.0 &&
         o->speed - o->speed == 0.0 && o->enhancement >= 0.0 &&
         o->enhancement - o->enhancement == 0.0 && finite_positive(o->object_mass);
}

myb_status myb_meteoroid_rate(const myb_meteoroid_options* options, double* rate) {
  MYB_CHECK(options != NULL && rate != NULL);
  if (!valid_options(options)) {
    return MYB_E_DOMAIN;
  }
  const double band = myb_meteoroid_flux(options->mass_min) - myb_meteoroid_flux(options->mass_max);
  double scale = 1.0;
  if (options->distance > 0.0) {
    scale = myb_exp(-options->exponent * myb_log(options->distance / AU));
  }
  const double r = 4.0 * options->cross_section * band * scale;
  if (!finite_positive(r)) {
    return MYB_E_DOMAIN;
  }
  *rate = r;
  return MYB_OK;
}

/* ---- The mass of an impact ------------------------------------------------------------------- */

typedef struct mass_target {
  double level; /* F(m) looked for */
} mass_target;

/* level - F(s), increasing in s = ln m (grams). */
static void mass_eval(const void* context, double s, double* f, double* df) {
  const mass_target* target = context;
  double slope = 0.0;
  *f = target->level - grun(s, &slope);
  *df = -slope;
}

myb_status myb_meteoroid_mass(const myb_meteoroid_options* options, double u, double* mass) {
  MYB_CHECK(options != NULL && mass != NULL);
  if (!valid_options(options) || !(u >= 0.0 && u < 1.0)) {
    return MYB_E_DOMAIN;
  }
  if (u == 0.0) {
    *mass = options->mass_min;
    return MYB_OK;
  }
  const double lo = myb_log(options->mass_min) + LN_GRAMS;
  const double hi = myb_log(options->mass_max) + LN_GRAMS;
  const double top = grun(lo, NULL);
  const double bottom = grun(hi, NULL);
  const mass_target target = {top - u * (top - bottom)};
  const myb_internal_root_function function = {mass_eval, &target};
  double s = 0.0;
  const myb_status st =
      myb_internal_newton_bracketed(&function, lo, hi, 0.5 * (lo + hi), MASS_ITERATIONS, &s);
  if (st != MYB_OK) {
    return st;
  }
  const double m = myb_exp(s - LN_GRAMS);
  *mass =
      m < options->mass_min ? options->mass_min : (m > options->mass_max ? options->mass_max : m);
  return MYB_OK;
}

/* ---- The impacts ----------------------------------------------------------------------------- */

myb_status myb_meteoroid_next(const myb_meteoroid_options* options, myb_rng* rng, double t,
                              myb_meteoroid_impact* impact) {
  MYB_CHECK(options != NULL && rng != NULL && impact != NULL);
  if (t - t != 0.0) {
    return MYB_E_DOMAIN;
  }
  double rate = 0.0;
  myb_status st = myb_meteoroid_rate(options, &rate);
  if (st != MYB_OK) {
    return st;
  }
  /* the wait: u in [0, 1), so 1 - u in (0, 1] and the logarithm is finite */
  const double wait = -myb_log1p(-myb_rng_uniform(rng)) / rate;
  double mass = 0.0;
  st = myb_meteoroid_mass(options, myb_rng_uniform(rng), &mass);
  if (st != MYB_OK) {
    return st;
  }
  double speed = options->speed > 0.0 ? options->speed : DEFAULT_SPEED;
  const double when = t + wait;
  if (options->speed_function != NULL) {
    st = options->speed_function(options->speed_context, rng, when, &speed);
    if (st != MYB_OK) {
      return st;
    }
    if (!finite_positive(speed)) {
      return MYB_E_DOMAIN;
    }
  }
  double direction[3];
  if (options->direction_function != NULL) {
    st = options->direction_function(options->direction_context, rng, when, direction);
    if (st != MYB_OK) {
      return st;
    }
    const double length = myb_vec3_norm(direction);
    if (!finite_positive(length)) {
      return MYB_E_DOMAIN;
    }
    for (int i = 0; i < 3; ++i) {
      direction[i] /= length;
    }
  } else {
    const double z = 2.0 * myb_rng_uniform(rng) - 1.0;
    const double phi = TWO_PI * myb_rng_uniform(rng);
    const double across = sqrt((1.0 - z) * (1.0 + z));
    direction[0] = across * myb_cos(phi);
    direction[1] = across * myb_sin(phi);
    direction[2] = z;
  }
  const double beta = options->enhancement > 0.0 ? options->enhancement : 1.0;
  const double momentum = beta * mass * speed;
  impact->t = when;
  impact->mass = mass;
  impact->speed = speed;
  for (int i = 0; i < 3; ++i) {
    impact->direction[i] = direction[i];
    impact->impulse[i] = momentum * direction[i];
    impact->dv[i] = impact->impulse[i] / options->object_mass;
  }
  return MYB_OK;
}

myb_status myb_meteoroid_propagate(myb_orbit* orbit, const myb_meteoroid_options* options,
                                   myb_rng* rng, double t_end, myb_meteoroid_impact* impacts,
                                   size_t capacity, size_t* count) {
  MYB_CHECK(orbit != NULL && options != NULL && rng != NULL && count != NULL &&
            (impacts != NULL || capacity == 0));
  *count = 0;
  double t = myb_orbit_time(orbit);
  double rate = 0.0;
  myb_status st = myb_meteoroid_rate(options, &rate);
  if (st != MYB_OK) {
    return st;
  }
  if (!(t_end >= t) || t_end - t_end != 0.0) {
    return MYB_E_DOMAIN;
  }
  /* far beyond the Poisson tail: the count exceeds lambda T by 20 of its sigma with a
   * probability below 1e-80 */
  const double expected = rate * (t_end - t);
  const double budget = 1000.0 + expected + 20.0 * sqrt(expected);
  if (!(budget < 9.0e15)) {
    return MYB_E_DOMAIN;
  }
  const size_t limit = (size_t)ceil(budget);
  double y[6];
  /* bound: limit impacts */
  for (size_t k = 0; k <= limit; ++k) {
    if (k == limit) {
      return MYB_E_NOCONV;
    }
    myb_meteoroid_impact impact;
    st = myb_meteoroid_next(options, rng, t, &impact);
    if (st != MYB_OK) {
      return st;
    }
    if (!(impact.t < t_end)) {
      break;
    }
    st = myb_orbit_state(orbit, impact.t, y);
    if (st == MYB_OK) {
      st = myb_orbit_kick(orbit, impact.dv);
    }
    if (st != MYB_OK) {
      return st;
    }
    if (*count < capacity) {
      impacts[*count] = impact;
    }
    ++*count;
    t = impact.t;
  }
  return myb_orbit_state(orbit, t_end, y);
}
