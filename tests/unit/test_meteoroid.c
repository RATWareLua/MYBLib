/* SPDX-License-Identifier: BSD-3-Clause */
/* Random impacts of micrometeoroids (myb/meteoroid.h):
 *  - the flux of Grun et al. (1985) against the formula evaluated here with the C library's
 *    pow, and its scaling with the distance from the Sun;
 *  - the number of impacts over a long time against lambda T, and the counts of many runs
 *    against Poisson's statistics (variance = mean);
 *  - the masses against the distribution F gives, by the test of Kolmogorov and Smirnov, and the
 *    inverse of the distribution against F itself;
 *  - the same impacts from the same seed, bit for bit, and others from another;
 *  - the total impulse: its sum against lambda T beta v E[m], each component against zero mean
 *    and the variance lambda T (beta v)^2 E[m^2] / 3 of an isotropic flux, E[m] and E[m^2] by
 *    quadrature of F;
 *  - the impacts applied to an orbit with myb_orbit_kick against the free motion they give;
 *  - the functions of the caller for the speed and the direction; what is refused. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/meteoroid.h"
#include "myb/model.h"
#include "myb/orbit.h"
#include "myb/uncert.h"

/* F of Grun et al. with the C library's pow, m in kg. */
static double grun_reference(double mass) {
  const double m = 1.0e3 * mass;
  return pow(2.2e3 * pow(m, 0.306) + 15.0, -4.38) +
         1.3e-9 * pow(m + 1.0e11 * m * m + 1.0e27 * pow(m, 4.0), -0.36) +
         1.3e-16 * pow(m + 1.0e6 * m * m, -0.85);
}

static myb_meteoroid_options base_options(void) {
  return (myb_meteoroid_options){
      .cross_section = 1.0, .mass_min = 1.0e-15, .mass_max = 1.0e-9, .object_mass = 100.0};
}

/* The fraction of the impacts below m: (F(m_min) - F(m)) / (F(m_min) - F(m_max)). */
static double cumulative(const myb_meteoroid_options* o, double m) {
  const double top = myb_meteoroid_flux(o->mass_min);
  const double bottom = myb_meteoroid_flux(o->mass_max);
  return (top - myb_meteoroid_flux(m)) / (top - bottom);
}

/* E[m^k] for k = 1, 2: m_min^k + int k m^(k-1) (1 - G(m)) dm, Simpson in ln m. */
static double mass_moment(const myb_meteoroid_options* o, int k) {
  const int n = 4000;
  const double a = log(o->mass_min);
  const double b = log(o->mass_max);
  const double h = (b - a) / n;
  double sum = 0.0;
  for (int i = 0; i <= n; ++i) {
    const double m = exp(a + h * i);
    const double f = k * pow(m, (double)k) * (1.0 - cumulative(o, m));
    const double w = (i == 0 || i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
    sum += w * f;
  }
  return pow(o->mass_min, (double)k) + sum * h / 3.0;
}

static int compare_doubles(const void* a, const void* b) {
  const double x = *(const double*)a;
  const double y = *(const double*)b;
  return (x > y) - (x < y);
}

/* ---- The flux -------------------------------------------------------------------------------- */

static void test_flux(void) {
  double worst = 0.0;
  for (int k = 0; k <= 200; ++k) {
    const double m = pow(10.0, -21.0 + 20.0 * k / 200.0); /* 1e-21 .. 0.1 kg */
    const double reference = grun_reference(m);
    worst = fmax(worst, fabs(myb_meteoroid_flux(m) - reference) / reference);
  }
  const double year = 365.25 * 86400.0;
  printf(
      "flux of Grun et al.: against pow %.1e (1e-21 .. 0.1 kg); above 1e-15 kg %.3e, "
      "1e-12 kg %.3e, 1e-9 kg %.3e, 1e-6 kg %.3e per m^2 per year\n",
      worst, myb_meteoroid_flux(1e-15) * year, myb_meteoroid_flux(1e-12) * year,
      myb_meteoroid_flux(1e-9) * year, myb_meteoroid_flux(1e-6) * year);
  MT_CHECK(worst < 1e-13);
  MT_CHECK(isnan(myb_meteoroid_flux(0.0)) && isnan(myb_meteoroid_flux(-1.0)) &&
           isnan(myb_meteoroid_flux((double)INFINITY)));
  /* decreasing */
  int monotone = 1;
  for (int k = 1; k <= 400; ++k) {
    monotone = monotone && myb_meteoroid_flux(pow(10.0, -21.0 + 0.05 * k)) <
                               myb_meteoroid_flux(pow(10.0, -21.0 + 0.05 * (k - 1)));
  }
  MT_CHECK(monotone);
  /* the rate and its scaling with the distance */
  myb_meteoroid_options o = base_options();
  double at_1au = 0.0;
  double at_2au = 0.0;
  MT_CHECK(myb_meteoroid_rate(&o, &at_1au) == MYB_OK);
  MT_CHECK(fabs(at_1au - 4.0 * (grun_reference(1e-15) - grun_reference(1e-9))) < 1e-13 * at_1au);
  o.distance = 2.0 * 149597870700.0;
  o.exponent = 1.3;
  MT_CHECK(myb_meteoroid_rate(&o, &at_2au) == MYB_OK);
  MT_CHECK(fabs(at_2au / at_1au - pow(2.0, -1.3)) < 1e-15);
}

/* ---- The number of impacts ------------------------------------------------------------------- */

static void test_counts(void) {
  const myb_meteoroid_options o = base_options();
  double rate = 0.0;
  MT_CHECK(myb_meteoroid_rate(&o, &rate) == MYB_OK);
  myb_rng rng;
  myb_rng_seed(&rng, 2026);
  /* one long run of some 50 000 impacts */
  const double span = 5.0e4 / rate;
  long count = 0;
  double t = 0.0;
  myb_meteoroid_impact impact;
  /* bound: the impacts of the span, well within 10^6 */
  for (long k = 0; k < 1000000; ++k) {
    MT_CHECK(myb_meteoroid_next(&o, &rng, t, &impact) == MYB_OK);
    if (impact.t >= span) {
      break;
    }
    t = impact.t;
    ++count;
  }
  const double expected = rate * span;
  const double sigmas = ((double)count - expected) / sqrt(expected);
  /* 2000 runs of a mean of 20: the mean and the variance of the counts */
  const int runs = 2000;
  const double short_span = 20.0 / rate;
  double sum = 0.0;
  double sum2 = 0.0;
  for (int run = 0; run < runs; ++run) {
    double tr = 0.0;
    int n = 0;
    /* bound: far beyond a mean of 20 */
    for (int k = 0; k < 10000; ++k) {
      MT_CHECK(myb_meteoroid_next(&o, &rng, tr, &impact) == MYB_OK);
      if (impact.t >= short_span) {
        break;
      }
      tr = impact.t;
      ++n;
    }
    sum += n;
    sum2 += (double)n * n;
  }
  const double mean = sum / runs;
  const double variance = (sum2 - sum * sum / runs) / (runs - 1);
  const double mean_sigmas = (mean - 20.0) / sqrt(20.0 / runs);
  /* the variance of the sample variance of Poisson counts: (lambda + 2 lambda^2) / runs */
  const double variance_sigmas = (variance - 20.0) / sqrt((20.0 + 2.0 * 400.0) / runs);
  printf(
      "impacts: %.4e per s on 1 m^2 between 1e-15 and 1e-9 kg; %ld in %.4e s against "
      "%.1f (%.2f sigma); %d runs of mean 20: mean %.4f (%.2f sigma), variance %.3f "
      "(%.2f sigma)\n",
      rate, count, span, expected, sigmas, runs, mean, mean_sigmas, variance, variance_sigmas);
  MT_CHECK(fabs(sigmas) < 4.0 && fabs(mean_sigmas) < 4.0 && fabs(variance_sigmas) < 4.0);
}

/* ---- The masses ------------------------------------------------------------------------------ */

static void test_masses(void) {
  const myb_meteoroid_options o = base_options();
  /* the inverse of the distribution */
  double worst = 0.0;
  for (int k = 0; k < 1000; ++k) {
    const double u = k / 1000.0;
    double m = 0.0;
    MT_CHECK(myb_meteoroid_mass(&o, u, &m) == MYB_OK);
    MT_CHECK(m >= o.mass_min && m <= o.mass_max);
    worst = fmax(worst, fabs(cumulative(&o, m) - u));
  }
  double m0 = 0.0;
  MT_CHECK(myb_meteoroid_mass(&o, 0.0, &m0) == MYB_OK);

  /* the masses drawn, by Kolmogorov and Smirnov */
  enum { SAMPLES = 50000 };

  static double g[SAMPLES];
  myb_rng rng;
  myb_rng_seed(&rng, 7);
  double t = 0.0;
  double smallest = o.mass_max;
  for (int k = 0; k < SAMPLES; ++k) {
    myb_meteoroid_impact impact;
    MT_CHECK(myb_meteoroid_next(&o, &rng, t, &impact) == MYB_OK);
    t = impact.t;
    smallest = fmin(smallest, impact.mass);
    g[k] = cumulative(&o, impact.mass);
  }
  qsort(g, SAMPLES, sizeof g[0], compare_doubles);
  double d = 0.0;
  for (int k = 0; k < SAMPLES; ++k) {
    d = fmax(d, fmax((k + 1.0) / (double)SAMPLES - g[k], g[k] - (double)k / (double)SAMPLES));
  }
  const double statistic = sqrt((double)SAMPLES) * d;
  /* the median mass */
  double median = 0.0;
  MT_CHECK(myb_meteoroid_mass(&o, 0.5, &median) == MYB_OK);
  printf(
      "masses: the inverse of the distribution to %.1e; %d drawn: Kolmogorov-Smirnov "
      "sqrt(n) D = %.3f (1.95 at 0.001); the median %.4e kg, the smallest drawn %.4e kg\n",
      worst, SAMPLES, statistic, median, smallest);
  MT_CHECK(worst < 1e-12 && m0 == o.mass_min);
  MT_CHECK(statistic < 1.95);
}

/* ---- Reproducibility ------------------------------------------------------------------------- */

static void test_seed(void) {
  const myb_meteoroid_options o = base_options();
  myb_rng a;
  myb_rng b;
  myb_rng c;
  myb_rng_seed(&a, 99);
  myb_rng_seed(&b, 99);
  myb_rng_seed(&c, 100);
  double ta = 0.0;
  double tb = 0.0;
  double tc = 0.0;
  int same = 1;
  int differ = 0;
  for (int k = 0; k < 1000; ++k) {
    myb_meteoroid_impact ia;
    myb_meteoroid_impact ib;
    myb_meteoroid_impact ic;
    MT_CHECK(myb_meteoroid_next(&o, &a, ta, &ia) == MYB_OK);
    MT_CHECK(myb_meteoroid_next(&o, &b, tb, &ib) == MYB_OK);
    MT_CHECK(myb_meteoroid_next(&o, &c, tc, &ic) == MYB_OK);
    same = same && memcmp(&ia, &ib, sizeof ia) == 0;
    differ += ia.t != ic.t && ia.mass != ic.mass;
    ta = ia.t;
    tb = ib.t;
    tc = ic.t;
  }
  /* a fixed value of the sequence of seed 99 (the platform does not enter) */
  myb_rng_seed(&a, 99);
  myb_meteoroid_impact first;
  MT_CHECK(myb_meteoroid_next(&o, &a, 0.0, &first) == MYB_OK);
  printf(
      "seed: 1000 impacts the same bit for bit from one seed, %d of 1000 differ from another; "
      "the first of seed 99: t = %a s, m = %a kg, direction (%a, %a, %a)\n",
      differ, first.t, first.mass, first.direction[0], first.direction[1], first.direction[2]);
  MT_CHECK(same && differ > 990);
}

/* ---- The total impulse ----------------------------------------------------------------------- */

static void test_impulse(void) {
  myb_meteoroid_options o = base_options();
  o.mass_min = 1.0e-12;
  o.mass_max = 1.0e-11;
  o.enhancement = 2.0;
  double rate = 0.0;
  MT_CHECK(myb_meteoroid_rate(&o, &rate) == MYB_OK);
  const double beta_v = 2.0 * 2.0e4; /* the default speed */
  const double mean_mass = mass_moment(&o, 1);
  const double mean_mass2 = mass_moment(&o, 2);
  const int runs = 1000;
  const double span = 100.0 / rate;
  myb_rng rng;
  myb_rng_seed(&rng, 31415);
  double scalar = 0.0;
  double sums[3] = {0.0, 0.0, 0.0};
  double squares = 0.0;
  long impacts = 0;
  double worst_unit = 0.0;
  for (int run = 0; run < runs; ++run) {
    double j[3] = {0.0, 0.0, 0.0};
    double t = 0.0;
    /* bound: far beyond a mean of 100 */
    for (int k = 0; k < 100000; ++k) {
      myb_meteoroid_impact impact;
      MT_CHECK(myb_meteoroid_next(&o, &rng, t, &impact) == MYB_OK);
      if (impact.t >= span) {
        break;
      }
      t = impact.t;
      ++impacts;
      scalar += beta_v * impact.mass;
      worst_unit = fmax(worst_unit, fabs(sqrt(impact.direction[0] * impact.direction[0] +
                                              impact.direction[1] * impact.direction[1] +
                                              impact.direction[2] * impact.direction[2]) -
                                         1.0));
      for (int i = 0; i < 3; ++i) {
        j[i] += impact.impulse[i];
        MT_CHECK(impact.dv[i] == impact.impulse[i] / o.object_mass);
      }
    }
    for (int i = 0; i < 3; ++i) {
      sums[i] += j[i];
      squares += j[i] * j[i];
    }
  }
  const double total_time = runs * span;
  const double expected_scalar = rate * total_time * beta_v * mean_mass;
  const double scalar_sigma = sqrt(rate * total_time) * beta_v * sqrt(mean_mass2);
  const double component_variance = rate * span * beta_v * beta_v * mean_mass2 / 3.0;
  double worst_mean = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst_mean = fmax(worst_mean, fabs(sums[i] / runs) / sqrt(component_variance / runs));
  }
  const double variance = squares / (3.0 * runs);
  /* the relative spread of a variance of 3 runs nearly normal components */
  const double variance_sigmas = (variance / component_variance - 1.0) / sqrt(2.0 / (3.0 * runs));
  printf(
      "impulse (beta 2, 20 km/s, 1e-12 .. 1e-11 kg, E[m] %.4e kg): %ld impacts; the sum of "
      "|J| against lambda T beta v E[m] %.2f sigma; the mean of a component of a run "
      "%.2f sigma from 0; its variance %.4e against %.4e N^2 s^2 (%.2f sigma)\n",
      mean_mass, impacts, (scalar - expected_scalar) / scalar_sigma, worst_mean, variance,
      component_variance, variance_sigmas);
  MT_CHECK(fabs(scalar - expected_scalar) < 4.0 * scalar_sigma);
  MT_CHECK(worst_mean < 4.0 && fabs(variance_sigmas) < 4.0);
  MT_CHECK(worst_unit < 1e-15);
}

/* ---- The impacts on an orbit ----------------------------------------------------------------- */

static myb_status no_force(void* context, double t, const double x[3], const double v[3],
                           double a[3]) {
  (void)context;
  (void)t;
  (void)x;
  (void)v;
  a[0] = a[1] = a[2] = 0.0;
  return MYB_OK;
}

static void test_orbit(void) {
  myb_meteoroid_options o = base_options();
  o.object_mass = 1.0e-3; /* a gram: the impulses show */
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "nothing", no_force, myb_term_no_jacobian, NULL, NULL) ==
           MYB_OK);
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, NULL, NULL, &orbit) == MYB_OK);
  double rate = 0.0;
  MT_CHECK(myb_meteoroid_rate(&o, &rate) == MYB_OK);
  const double t_end = 300.0 / rate;
  const double y0[6] = {1.0e3, -2.0e3, 5.0e2, 3.0, -1.0, 2.0};
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, 2.0 * t_end) == MYB_OK);

  enum { CAPACITY = 1000 };

  static myb_meteoroid_impact impacts[CAPACITY];
  size_t count = 0;
  myb_rng rng;
  myb_rng_seed(&rng, 5);
  MT_CHECK(myb_meteoroid_propagate(orbit, &o, &rng, t_end, impacts, CAPACITY, &count) == MYB_OK);
  MT_CHECK(myb_orbit_time(orbit) == t_end && count <= CAPACITY);
  double y[6];
  MT_CHECK(myb_orbit_state(orbit, t_end, y) == MYB_OK);
  /* free motion with the kicks: x = x0 + v0 T + sum dv_i (T - t_i), v = v0 + sum dv_i */
  double x[3];
  double v[3];
  double dv_total = 0.0;
  for (int i = 0; i < 3; ++i) {
    x[i] = y0[i] + y0[3 + i] * t_end;
    v[i] = y0[3 + i];
  }
  for (size_t k = 0; k < count; ++k) {
    for (int i = 0; i < 3; ++i) {
      x[i] += impacts[k].dv[i] * (t_end - impacts[k].t);
      v[i] += impacts[k].dv[i];
    }
    dv_total += sqrt(impacts[k].dv[0] * impacts[k].dv[0] + impacts[k].dv[1] * impacts[k].dv[1] +
                     impacts[k].dv[2] * impacts[k].dv[2]);
  }
  double dx = 0.0;
  double dvv = 0.0;
  for (int i = 0; i < 3; ++i) {
    dx = fmax(dx, fabs(y[i] - x[i]) / fabs(x[i]));
    dvv = fmax(dvv, fabs(y[3 + i] - v[i]));
  }
  printf(
      "on an orbit: %zu impacts in %.4e s, sum of |dv| %.4e m/s; the free motion with the "
      "kicks to %.1e (position, relative) and %.1e m/s\n",
      count, t_end, dv_total, dx, dvv);
  MT_CHECK(count > 200 && dx < 1e-13 && dvv < 1e-15);
  /* behind the time reached */
  MT_CHECK(myb_meteoroid_propagate(orbit, &o, &rng, 0.5 * t_end, NULL, 0, &count) == MYB_E_DOMAIN);
  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);
}

/* ---- The caller's distributions; what is refused --------------------------------------------- */

static myb_status uniform_speed(void* context, myb_rng* rng, double t, double* speed) {
  (void)context;
  (void)t;
  *speed = 1.0e4 + 2.0e4 * myb_rng_uniform(rng);
  return MYB_OK;
}

static myb_status from_apex(void* context, myb_rng* rng, double t, double direction[3]) {
  (void)rng;
  (void)t;
  memcpy(direction, context, 3 * sizeof *direction);
  return MYB_OK;
}

static myb_status nowhere(void* context, myb_rng* rng, double t, double direction[3]) {
  (void)context;
  (void)rng;
  (void)t;
  direction[0] = direction[1] = direction[2] = 0.0;
  return MYB_OK;
}

static void test_caller(void) {
  myb_meteoroid_options o = base_options();
  double apex[3] = {0.0, -3.0, 4.0}; /* any length */
  o.speed_function = uniform_speed;
  o.direction_function = from_apex;
  o.direction_context = apex;
  o.enhancement = 3.0;
  myb_rng rng;
  myb_rng_seed(&rng, 11);
  double slowest = 1e300;
  double fastest = 0.0;
  double worst = 0.0;
  double t = 0.0;
  for (int k = 0; k < 2000; ++k) {
    myb_meteoroid_impact impact;
    MT_CHECK(myb_meteoroid_next(&o, &rng, t, &impact) == MYB_OK);
    t = impact.t;
    slowest = fmin(slowest, impact.speed);
    fastest = fmax(fastest, impact.speed);
    const double p = 3.0 * impact.mass * impact.speed;
    worst = fmax(worst, fabs(impact.impulse[0]) + fabs(impact.impulse[1] - p * -0.6) +
                            fabs(impact.impulse[2] - p * 0.8));
    worst = fmax(worst, fabs(impact.direction[1] + 0.6) + fabs(impact.direction[2] - 0.8));
  }
  MT_CHECK(slowest >= 1.0e4 && fastest < 3.0e4 && fastest - slowest > 1.9e4);
  MT_CHECK(worst < 1e-15);
  /* refused */
  myb_meteoroid_impact impact;
  o.direction_function = nowhere;
  MT_CHECK(myb_meteoroid_next(&o, &rng, 0.0, &impact) == MYB_E_DOMAIN);
  const myb_meteoroid_options good = base_options();
  myb_meteoroid_options bad = good;
  double value = 0.0;
  bad.mass_max = bad.mass_min;
  MT_CHECK(myb_meteoroid_rate(&bad, &value) == MYB_E_DOMAIN);
  bad = good;
  bad.cross_section = 0.0;
  MT_CHECK(myb_meteoroid_rate(&bad, &value) == MYB_E_DOMAIN);
  bad = good;
  bad.object_mass = -1.0;
  MT_CHECK(myb_meteoroid_next(&bad, &rng, 0.0, &impact) == MYB_E_DOMAIN);
  bad = good;
  bad.speed = (double)NAN;
  MT_CHECK(myb_meteoroid_rate(&bad, &value) == MYB_E_DOMAIN);
  MT_CHECK(myb_meteoroid_mass(&good, 1.0, &value) == MYB_E_DOMAIN);
  MT_CHECK(myb_meteoroid_mass(&good, -0.1, &value) == MYB_E_DOMAIN);
  MT_CHECK(myb_meteoroid_next(&good, &rng, (double)INFINITY, &impact) == MYB_E_DOMAIN);
}

int main(void) {
  test_flux();
  test_counts();
  test_masses();
  test_seed();
  test_impulse();
  test_orbit();
  test_caller();
  return mt_finish("test_meteoroid");
}
