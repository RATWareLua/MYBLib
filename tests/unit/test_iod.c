/* SPDX-License-Identifier: BSD-3-Clause */
/* Initial orbit determination (myb/iod.h):
 *  - Gibbs: three points of Kepler orbits of every kind give the velocity of the orbit; the
 *    example 7-3 of Vallado; points off one plane are measured;
 *  - Herrick-Gibbs: the error against the spacing (fourth power); the example 7-4 of Vallado;
 *  - Gauss: directions from three places of a turning Earth to objects on low, eccentric and
 *    high orbits, over short and long arcs — the orbit comes back to rounding once refined;
 *    how far the series of Gauss alone is; what is refused;
 *  - Laplace: the same arcs;
 *  - the double-r iteration: arcs of a third to two thirds of a turn from rough distances;
 *  - Gooding: the same long arcs, and arcs of one and two whole turns and more with the
 *    branches of Lambert's problem, from rough ranges. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/iod.h"
#include "myb/kepler.h"

static const double MU = 3.986004418e14;
static const double PI = 3.141592653589793;

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static double worst3(const double a[3], const double b[3]) {
  double w = 0.0;
  for (int i = 0; i < 3; ++i) {
    w = fmax(w, fabs(a[i] - b[i]));
  }
  return w / norm3(b);
}

static myb_pv orbit(double a, double e, double i, double raan, double argp, double nu) {
  const myb_elements el = {.p = a * (1.0 - e * e), .e = e, .i = i, .raan = raan, .argp = argp,
                           .nu = nu};
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(MU, &el, &pv) == MYB_OK);
  return pv;
}

static void test_gibbs(void) {
  /* points of orbits from circular to hyperbolic, a few tens of degrees apart */
  const double e[5] = {0.0, 0.2, 0.7, 0.95, 1.4};
  double worst = 0.0;
  double worst_plane = 0.0;
  for (int k = 0; k < 5; ++k) {
    const double a = e[k] < 1.0 ? 2.0e7 : -2.0e7;
    const double p = a * (1.0 - e[k] * e[k]);
    const myb_elements el = {.p = p, .e = e[k], .i = 0.9, .raan = 2.0, .argp = 1.0, .nu = -0.5};
    myb_pv pv[3];
    for (int j = 0; j < 3; ++j) {
      myb_elements at = el;
      at.nu = -0.5 + 0.45 * j;
      MT_CHECK(myb_pv_from_elements(MU, &at, &pv[j]) == MYB_OK);
    }
    double v2[3];
    double plane = 1.0;
    MT_CHECK(myb_iod_gibbs(MU, pv[0].r, pv[1].r, pv[2].r, v2, &plane) == MYB_OK);
    worst = fmax(worst, worst3(v2, pv[1].v));
    worst_plane = fmax(worst_plane, fabs(plane));
  }
  printf("Gibbs, e from 0 to 1.4: the velocity %.2e, out of plane %.1e rad\n", worst,
         worst_plane);
  MT_CHECK(worst < 2e-14 && worst_plane < 1e-15);
  /* Vallado, example 7-3 (km, km/s) */
  const double mu_km = 398600.4418;
  const double r1[3] = {0.0, 0.0, 6378.137};
  const double r2[3] = {0.0, -4464.696, -5102.509};
  const double r3[3] = {0.0, 5740.323, 3189.068};
  double v2[3];
  MT_CHECK(myb_iod_gibbs(mu_km, r1, r2, r3, v2, NULL) == MYB_OK);
  printf("  Vallado 7-3: v2 = (%.6f, %.6f, %.6f) km/s (0, 5.531148, -5.191806)\n", v2[0], v2[1],
         v2[2]);
  MT_CHECK(fabs(v2[0]) < 1e-12 && fabs(v2[1] - 5.531148) < 2e-6 && fabs(v2[2] + 5.191806) < 2e-6);
  /* a point lifted out of the plane: the angle says so */
  const myb_pv a = orbit(1.0e7, 0.1, 0.5, 1.0, 2.0, 0.0);
  const myb_pv b = orbit(1.0e7, 0.1, 0.5, 1.0, 2.0, 0.6);
  const myb_pv c = orbit(1.0e7, 0.1, 0.5, 1.0, 2.0, 1.3);
  double h[3] = {b.r[1] * c.r[2] - b.r[2] * c.r[1], b.r[2] * c.r[0] - b.r[0] * c.r[2],
                 b.r[0] * c.r[1] - b.r[1] * c.r[0]};
  const double hn = norm3(h);
  double lifted[3];
  for (int i = 0; i < 3; ++i) {
    lifted[i] = a.r[i] + 1.0e4 * h[i] / hn;
  }
  double plane = 0.0;
  MT_CHECK(myb_iod_gibbs(MU, lifted, b.r, c.r, v2, &plane) == MYB_OK);
  MT_CHECK(fabs(plane - 1.0e4 / norm3(lifted)) < 1e-9);
  /* refused */
  const double zero[3] = {0.0, 0.0, 0.0};
  const double twice[3] = {2.0 * a.r[0], 2.0 * a.r[1], 2.0 * a.r[2]};
  const double thrice[3] = {3.0 * a.r[0], 3.0 * a.r[1], 3.0 * a.r[2]};
  MT_CHECK(myb_iod_gibbs(MU, zero, b.r, c.r, v2, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_iod_gibbs(MU, a.r, twice, thrice, v2, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_iod_gibbs(-1.0, a.r, b.r, c.r, v2, NULL) == MYB_E_DOMAIN);
}

static void test_herrick_gibbs(void) {
  /* the same orbit sampled closer and closer: the error falls as the fourth power */
  const myb_pv middle = orbit(7.2e6, 0.05, 1.1, 0.3, 0.8, 0.7);
  const double period = 2.0 * PI * sqrt(7.2e6 * 7.2e6 * 7.2e6 / MU);
  double errors[3];
  const double spacing[3] = {period / 40.0, period / 80.0, period / 160.0};
  for (int k = 0; k < 3; ++k) {
    myb_pv before;
    myb_pv after;
    MT_CHECK(myb_kepler_propagate(MU, &middle, -spacing[k], &before) == MYB_OK);
    MT_CHECK(myb_kepler_propagate(MU, &middle, spacing[k], &after) == MYB_OK);
    double v2[3];
    MT_CHECK(myb_iod_herrick_gibbs(MU, before.r, -spacing[k], middle.r, 0.0, after.r, spacing[k],
                                   v2) == MYB_OK);
    errors[k] = worst3(v2, middle.v);
  }
  printf("Herrick-Gibbs, spacing 1/40, 1/80, 1/160 of the period: %.2e, %.2e, %.2e (ratios "
         "%.1f, %.1f)\n",
         errors[0], errors[1], errors[2], errors[0] / errors[1], errors[1] / errors[2]);
  MT_CHECK(errors[0] < 3e-5 && errors[2] < 1e-7);
  MT_CHECK(fabs(errors[0] / errors[1] - 16.0) < 1.0 && fabs(errors[1] / errors[2] - 16.0) < 1.0);
  /* Vallado, example 7-4 (km, s, km/s) */
  const double r1[3] = {3419.85564, 6019.82602, 2784.60022};
  const double r2[3] = {2935.91195, 6326.18324, 2660.59584};
  const double r3[3] = {2434.95202, 6597.38674, 2521.52311};
  double v2[3];
  MT_CHECK(myb_iod_herrick_gibbs(398600.4418, r1, 0.0, r2, 76.48, r3, 153.04, v2) == MYB_OK);
  printf("  Vallado 7-4: v2 = (%.6f, %.6f, %.6f) km/s (-6.441557, 3.777559, -1.720567)\n", v2[0],
         v2[1], v2[2]);
  MT_CHECK(fabs(v2[0] + 6.441557) < 2e-6 && fabs(v2[1] - 3.777559) < 2e-6 &&
           fabs(v2[2] + 1.720567) < 2e-6);
  MT_CHECK(myb_iod_herrick_gibbs(MU, r1, 0.0, r2, 0.0, r3, 1.0, v2) == MYB_E_DOMAIN);
}

/* A place on the surface of a sphere turning about z: where it is at t. */
static void place_at(double latitude, double longitude, double t, double out[3]) {
  const double radius = 6378137.0;
  const double angle = longitude + 7.292115e-5 * t;
  out[0] = radius * cos(latitude) * cos(angle);
  out[1] = radius * cos(latitude) * sin(angle);
  out[2] = radius * sin(latitude);
}

/* Directions from the place to an object of the orbit `middle` (at t = 0) at three times;
 * the orbit that Gauss finds nearest to the true one against it. */
static double gauss_case(const myb_pv* middle, double span, double* unrefined, int* count,
                         int* passes) {
  const double t[3] = {-0.45 * span, 0.0, 0.55 * span};
  double observer[3][3];
  double sight[3][3];
  for (int k = 0; k < 3; ++k) {
    myb_pv at;
    MT_CHECK(myb_kepler_propagate(MU, middle, t[k], &at) == MYB_OK);
    place_at(0.7, -1.9, t[k], observer[k]);
    double d[3];
    for (int i = 0; i < 3; ++i) {
      d[i] = at.r[i] - observer[k][i];
    }
    const double n = norm3(d);
    for (int i = 0; i < 3; ++i) {
      sight[k][i] = d[i] / n;
    }
  }
  myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS];
  size_t found = 0;
  const myb_status st = myb_iod_gauss(MU, t, observer, sight, solutions, &found);
  *count = (int)found;
  if (st != MYB_OK) {
    return 1.0;
  }
  double best = 1.0;
  for (size_t k = 0; k < found; ++k) {
    const double error = fmax(worst3(solutions[k].r, middle->r), worst3(solutions[k].v, middle->v));
    if (error < best) {
      best = error;
      *passes = solutions[k].iterations;
    }
    MT_CHECK(solutions[k].range[0] > 0.0 && solutions[k].range[1] > 0.0 &&
             solutions[k].range[2] > 0.0);
  }
  /* what the series alone would give: the first step of the same method, from the true
   * distance — the velocity from the three true positions and the truncated coefficients */
  const double x3 = pow(norm3(middle->r), 3.0);
  myb_pv ends[2];
  MT_CHECK(myb_kepler_propagate(MU, middle, t[0], &ends[0]) == MYB_OK);
  MT_CHECK(myb_kepler_propagate(MU, middle, t[2], &ends[1]) == MYB_OK);
  const double f1 = 1.0 - 0.5 * MU * t[0] * t[0] / x3;
  const double f3 = 1.0 - 0.5 * MU * t[2] * t[2] / x3;
  const double g1 = t[0] - MU * t[0] * t[0] * t[0] / (6.0 * x3);
  const double g3 = t[2] - MU * t[2] * t[2] * t[2] / (6.0 * x3);
  double v[3];
  for (int i = 0; i < 3; ++i) {
    v[i] = (-f3 * ends[0].r[i] + f1 * ends[1].r[i]) / (f1 * g3 - f3 * g1);
  }
  *unrefined = worst3(v, middle->v);
  return best;
}

/* The sights of an orbit (`middle` at t = 0) from the turning place at three times. */
static void sights_of(const myb_pv* middle, const double t[3], double observer[3][3],
                      double sight[3][3], double range[3]) {
  for (int k = 0; k < 3; ++k) {
    myb_pv at;
    MT_CHECK(myb_kepler_propagate(MU, middle, t[k], &at) == MYB_OK);
    place_at(0.7, -1.9, t[k], observer[k]);
    double d[3];
    for (int i = 0; i < 3; ++i) {
      d[i] = at.r[i] - observer[k][i];
    }
    range[k] = norm3(d);
    for (int i = 0; i < 3; ++i) {
      sight[k][i] = d[i] / range[k];
    }
  }
}

static double error_of(const myb_iod_solution* s, const myb_pv* middle) {
  return fmax(worst3(s->r, middle->r), worst3(s->v, middle->v));
}

/* The largest angle (rad) between a sight and the direction to the orbit of a solution: an
 * orbit through the three sights, true or not, has it at rounding. */
static double off_sights(const myb_iod_solution* s, const double t[3],
                         const double observer[3][3], const double sight[3][3]) {
  myb_pv middle;
  memcpy(middle.r, s->r, sizeof middle.r);
  memcpy(middle.v, s->v, sizeof middle.v);
  double worst = 0.0;
  for (int k = 0; k < 3; ++k) {
    myb_pv at;
    MT_CHECK(myb_kepler_propagate(MU, &middle, t[k] - t[1], &at) == MYB_OK);
    double d[3];
    for (int i = 0; i < 3; ++i) {
      d[i] = at.r[i] - observer[k][i];
    }
    const double cross[3] = {d[1] * sight[k][2] - d[2] * sight[k][1],
                             d[2] * sight[k][0] - d[0] * sight[k][2],
                             d[0] * sight[k][1] - d[1] * sight[k][0]};
    worst = fmax(worst, norm3(cross) / norm3(d));
  }
  return worst;
}

static void test_laplace(void) {
  const struct {
    const char* name;
    double a, e, i, nu, span;
  } cases[5] = {{"low orbit, 4 minutes", 7.0e6, 0.001, 0.9, 0.5, 240.0},
                {"low orbit, 12 minutes", 7.0e6, 0.001, 0.9, 0.5, 720.0},
                {"eccentric, 30 minutes", 2.6e7, 0.7, 1.1, 2.4, 1800.0},
                {"navigation orbit, 2 hours", 2.656e7, 0.01, 0.96, 1.0, 7200.0},
                {"geostationary height, 5 hours", 4.2164e7, 0.0003, 0.1, 0.3, 18000.0}};
  for (int k = 0; k < 5; ++k) {
    const myb_pv middle = orbit(cases[k].a, cases[k].e, cases[k].i, 1.3, 0.4, cases[k].nu);
    const double t[3] = {-0.45 * cases[k].span, 0.0, 0.55 * cases[k].span};
    double observer[3][3];
    double sight[3][3];
    double range[3];
    sights_of(&middle, t, observer, sight, range);
    myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS];
    size_t found = 0;
    MT_CHECK(myb_iod_laplace(MU, t, observer, sight, solutions, &found) == MYB_OK);
    double best = 1.0;
    int passes = 0;
    for (size_t s = 0; s < found; ++s) {
      if (error_of(&solutions[s], &middle) < best) {
        best = error_of(&solutions[s], &middle);
        passes = solutions[s].iterations;
      }
    }
    printf("Laplace, %-30s: %d candidate(s), the orbit to %.2e after %d passes\n",
           cases[k].name, (int)found, best, passes);
    MT_CHECK(found >= 1 && best < 1e-8);
  }
  /* sights that do not turn */
  const double t[3] = {0.0, 100.0, 200.0};
  double observer[3][3];
  for (int k = 0; k < 3; ++k) {
    place_at(0.7, -1.9, t[k], observer[k]);
  }
  const double still[3][3] = {{0.6, 0.8, 0.0}, {0.6, 0.8, 0.0}, {0.6, 0.8, 0.0}};
  myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS];
  size_t found = 0;
  MT_CHECK(myb_iod_laplace(MU, t, observer, still, solutions, &found) == MYB_E_DOMAIN);
}

/* Long arcs: the double-r iteration from distances 5 % off, Gooding from ranges 5 % off; what
 * Gauss makes of the same sights. Three sights may hold more than one orbit: over two thirds
 * of a turn both methods settle from 5 % off on another orbit through the same sights, and
 * on the true one from 1 % off. */
static void test_long_arcs(void) {
  const struct {
    const char* name;
    double a, e, i, nu, span; /* span in periods */
  } cases[4] = {{"low orbit, a third of a turn", 7.0e6, 0.001, 0.9, 0.5, 0.33},
                {"low orbit, two thirds of a turn", 7.0e6, 0.02, 1.7, 0.5, 0.66},
                {"eccentric, half a turn", 2.6e7, 0.6, 1.1, 2.4, 0.5},
                {"navigation orbit, 8 hours", 2.656e7, 0.01, 0.96, 1.0, 0.63}};
  for (int k = 0; k < 4; ++k) {
    const myb_pv middle = orbit(cases[k].a, cases[k].e, cases[k].i, 1.3, 0.4, cases[k].nu);
    const double period = 2.0 * PI * sqrt(pow(cases[k].a, 3.0) / MU);
    const double span = cases[k].span * period;
    const double t[3] = {-0.45 * span, 0.0, 0.55 * span};
    double observer[3][3];
    double sight[3][3];
    double range[3];
    sights_of(&middle, t, observer, sight, range);
    myb_pv first;
    MT_CHECK(myb_kepler_propagate(MU, &middle, t[0], &first) == MYB_OK);
    double error_r = 1.0;
    double error_g = 1.0;
    double other = 0.0;    /* how far another orbit through the sights is */
    double crossing = 0.0; /* how far from the sights the orbits found pass */
    for (int pass = 0; pass < 2 && (error_r > 1e-8 || error_g > 1e-8); ++pass) {
      const double off = pass == 0 ? 0.05 : 0.01;
      const double distance[2] = {(1.0 + off) * norm3(first.r), (1.0 - off) * norm3(middle.r)};
      const double ranges[2] = {(1.0 + off) * range[0], (1.0 - off) * range[2]};
      myb_iod_solution solution;
      MT_CHECK(myb_iod_double_r(MU, t, observer, sight, distance, &solution) == MYB_OK);
      error_r = error_of(&solution, &middle);
      crossing = fmax(crossing, off_sights(&solution, t, observer, sight));
      MT_CHECK(myb_iod_gooding(MU, t, observer, sight, ranges, 0, 0, &solution) == MYB_OK);
      error_g = error_of(&solution, &middle);
      crossing = fmax(crossing, off_sights(&solution, t, observer, sight));
      if (pass == 0 && error_r > 1e-8) {
        other = error_r;
      }
    }
    myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS];
    size_t found = 0;
    double error_gauss = 1.0;
    if (myb_iod_gauss(MU, t, observer, sight, solutions, &found) == MYB_OK) {
      for (size_t s = 0; s < found; ++s) {
        error_gauss = fmin(error_gauss, error_of(&solutions[s], &middle));
      }
    }
    printf("%-34s: double-r %.2e, Gooding %.2e; off the sights %.1e rad (Gauss: %d "
           "candidate(s), the nearest %.1e)\n",
           cases[k].name, error_r, error_g, crossing, (int)found, error_gauss);
    if (other > 0.0) {
      printf("  from 5 %% off another orbit through the same sights, %.2e from the true one; "
             "from 1 %% off the true one\n",
             other);
    }
    MT_CHECK(error_r < 1e-8 && error_g < 1e-8 && crossing < 1e-12);
  }
}

/* Whole turns between the first and the third sight: only Gooding. Each branch of Lambert's
 * problem with the turns of the truth is tried; the true orbit is one of them. */
static void test_gooding_turns(void) {
  const struct {
    const char* name;
    double a, e, i, nu, turns; /* turns between the first and the third time */
    int revolutions;
  } cases[3] = {{"low orbit, 1.85 turns", 7.0e6, 0.001, 0.9, 0.5, 1.85, 1},
                {"eccentric, 2.7 turns", 1.5e7, 0.4, 1.1, 2.4, 2.7, 2},
                {"navigation orbit, 3.3 turns", 2.656e7, 0.01, 0.96, 1.0, 3.3, 3}};
  for (int k = 0; k < 3; ++k) {
    const myb_pv middle = orbit(cases[k].a, cases[k].e, cases[k].i, 1.3, 0.4, cases[k].nu);
    const double period = 2.0 * PI * sqrt(pow(cases[k].a, 3.0) / MU);
    const double span = cases[k].turns * period;
    const double t[3] = {-0.3 * span, 0.0, 0.7 * span};
    double observer[3][3];
    double sight[3][3];
    double range[3];
    sights_of(&middle, t, observer, sight, range);
    const double ranges[2] = {1.05 * range[0], 0.95 * range[2]};
    double best = 1.0;
    int which = 0;
    for (int branch = -1; branch <= 1; branch += 2) {
      myb_iod_solution solution;
      const myb_status st = myb_iod_gooding(MU, t, observer, sight, ranges, cases[k].revolutions,
                                            branch, &solution);
      if (st == MYB_OK) {
        /* the other branch, where it settles, is another orbit through the same sights */
        MT_CHECK(off_sights(&solution, t, observer, sight) < 1e-12);
      }
      if (st == MYB_OK && error_of(&solution, &middle) < best) {
        best = error_of(&solution, &middle);
        which = branch;
      }
    }
    printf("Gooding, %-28s: %d turn(s), branch %+d gives the orbit to %.2e\n", cases[k].name,
           cases[k].revolutions, which, best);
    MT_CHECK(best < 1e-8);
  }
  /* refused */
  const myb_pv middle = orbit(7.0e6, 0.001, 0.9, 1.3, 0.4, 0.5);
  const double t[3] = {-1000.0, 0.0, 1000.0};
  double observer[3][3];
  double sight[3][3];
  double range[3];
  sights_of(&middle, t, observer, sight, range);
  myb_iod_solution solution;
  const double ranges[2] = {range[0], range[2]};
  MT_CHECK(myb_iod_gooding(MU, t, observer, sight, ranges, 9, 1, &solution) == MYB_E_DOMAIN);
  MT_CHECK(myb_iod_gooding(MU, t, observer, sight, ranges, 1, 0, &solution) == MYB_E_DOMAIN);
  const double negative[2] = {-range[0], range[2]};
  MT_CHECK(myb_iod_gooding(MU, t, observer, sight, negative, 0, 0, &solution) == MYB_E_DOMAIN);
  const double none[2] = {0.0, 7.0e6};
  MT_CHECK(myb_iod_double_r(MU, t, observer, sight, none, &solution) == MYB_E_DOMAIN);
}

static void test_gauss(void) {
  const struct {
    const char* name;
    double a, e, i, nu, span;
  } cases[5] = {{"low orbit, 4 minutes", 7.0e6, 0.001, 0.9, 0.5, 240.0},
                {"low orbit, 12 minutes", 7.0e6, 0.001, 0.9, 0.5, 720.0},
                {"eccentric, 30 minutes", 2.6e7, 0.7, 1.1, 2.4, 1800.0},
                {"navigation orbit, 2 hours", 2.656e7, 0.01, 0.96, 1.0, 7200.0},
                {"geostationary height, 5 hours", 4.2164e7, 0.0003, 0.1, 0.3, 18000.0}};
  for (int k = 0; k < 5; ++k) {
    const myb_pv middle = orbit(cases[k].a, cases[k].e, cases[k].i, 1.3, 0.4, cases[k].nu);
    double unrefined = 0.0;
    int count = 0;
    int passes = 0;
    const double error = gauss_case(&middle, cases[k].span, &unrefined, &count, &passes);
    printf("Gauss, %-30s: %d candidate(s), the orbit to %.2e after %d passes (the series "
           "alone: %.1e)\n",
           cases[k].name, count, error, passes, unrefined);
    MT_CHECK(count >= 1);
    MT_CHECK(error < 1e-8);
    MT_CHECK(error < 1e-3 * unrefined || unrefined < 1e-9);
  }
  /* refused: times out of order, a sight that is not a unit vector, three sights in a plane */
  const double t[3] = {0.0, 100.0, 200.0};
  const double back[3] = {0.0, 100.0, 50.0};
  double observer[3][3];
  double sight[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  for (int k = 0; k < 3; ++k) {
    place_at(0.7, -1.9, t[k], observer[k]);
  }
  myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS];
  size_t found = 7;
  MT_CHECK(myb_iod_gauss(MU, back, observer, sight, solutions, &found) == MYB_E_DOMAIN);
  MT_CHECK(found == 0);
  sight[0][0] = 1.001;
  MT_CHECK(myb_iod_gauss(MU, t, observer, sight, solutions, &found) == MYB_E_DOMAIN);
  const double flat[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.6, 0.8, 0.0}};
  MT_CHECK(myb_iod_gauss(MU, t, observer, flat, solutions, &found) == MYB_E_DOMAIN);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_gibbs();
  test_herrick_gibbs();
  test_gauss();
  test_laplace();
  test_long_arcs();
  test_gooding_turns();
  return mt_finish("test_iod");
}
