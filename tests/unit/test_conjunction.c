/* SPDX-License-Identifier: BSD-3-Clause */
/* Closest approaches and the probability of collision (myb/conjunction.h):
 *  - two circular orbits of one radius in planes 0.3 rad apart meet at their node: the time
 *    and the distance in closed form; coplanar circles: one approach a synodic period;
 *  - two propagated orbits against the same two as providers;
 *  - the probability: a centred round Gaussian (closed form), an offset round one (the series
 *    of Rice, which is Chan's), unequal sigmas against a sum over a fine grid, far tails, the
 *    plane of an encounter from the vectors, the largest probability over the size of the
 *    covariance. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/conjunction.h"
#include "myb/env.h"
#include "myb/kepler.h"
#include "myb/model.h"
#include "myb/uncert.h"

static const double MU = 3.986004418e14;
static const double PI = 3.141592653589793;

typedef struct kepler_orbit {
  myb_pv start; /* at t = 0 */
} kepler_orbit;

static myb_status kepler_state(void* context, myb_time t, double state[6]) {
  const kepler_orbit* o = context;
  myb_pv pv;
  const myb_status st = myb_kepler_propagate(MU, &o->start, myb_time_seconds(t), &pv);
  memcpy(state, pv.r, 3 * sizeof(double));
  memcpy(state + 3, pv.v, 3 * sizeof(double));
  return st;
}

static myb_pv orbit(double a, double e, double i, double raan, double argp, double nu) {
  const myb_elements el = {.p = a * (1.0 - e * e), .e = e, .i = i, .raan = raan, .argp = argp,
                           .nu = nu};
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(MU, &el, &pv) == MYB_OK);
  return pv;
}

static void test_approaches(void) {
  const double a = 7.2e6;
  const double n = sqrt(MU / (a * a * a));
  const double period = 2.0 * PI / n;
  /* the same circle in two planes through the x axis, both objects 20 degrees before the
   * node, the second 2 km behind the first along its orbit */
  const double lag = 2000.0 / a;
  kepler_orbit first = {orbit(a, 0.0, 0.0, 0.0, 0.0, -0.35)};
  kepler_orbit second = {orbit(a, 0.0, 0.3, 0.0, 0.0, -0.35 - lag)};
  myb_approach found[8];
  size_t count = 0;
  MT_CHECK(myb_approach_find(kepler_state, &first, kepler_state, &second, 0.0, 1.2 * period,
                             30.0, found, 8, &count) == MYB_OK);
  /* r (cos u, sin u, 0) and r (cos w, sin w cos i, sin w sin i), w = u - lag: the distance is
   * least where the two straddle the node, u + w = 2 k pi, and is 2 r sin(lag / 2) there to
   * the first order of the tilt */
  const double t_node = (0.35 + 0.5 * lag) / n;
  printf("two circles crossing at a node: %zu approaches in 1.2 periods; the first at %.6f s "
         "(%.6f), %.3f m apart; the second half a period later: %.6f s\n",
         count, found[0].t, t_node, found[0].distance, found[1].t - found[0].t - 0.5 * period);
  MT_CHECK(count == 3);
  MT_CHECK(fabs(found[0].t - t_node) < 1e-3);
  MT_CHECK(fabs(found[0].distance - 2000.0) < 30.0);
  MT_CHECK(fabs(found[1].t - found[0].t - 0.5 * period) < 1e-3);
  for (size_t k = 0; k < count; ++k) {
    /* at a minimum the range rate vanishes: to the resolution of the time, 1e-12 s */
    const double rate = (found[k].position[0] * found[k].velocity[0] +
                         found[k].position[1] * found[k].velocity[1] +
                         found[k].position[2] * found[k].velocity[2]) / found[k].distance;
    MT_CHECK(fabs(rate) < 1e-6);
  }
  /* coplanar circles: closest when in line, once a synodic period, by the difference of radii */
  kepler_orbit low = {orbit(7.0e6, 0.0, 0.5, 1.0, 0.0, 0.4)};
  kepler_orbit high = {orbit(7.3e6, 0.0, 0.5, 1.0, 0.0, 0.0)};
  const double n1 = sqrt(MU / (7.0e6 * 7.0e6 * 7.0e6));
  const double n2 = sqrt(MU / (7.3e6 * 7.3e6 * 7.3e6));
  const double synodic = 2.0 * PI / (n1 - n2);
  MT_CHECK(myb_approach_find(kepler_state, &low, kepler_state, &high, 0.0, 2.5 * synodic, 60.0,
                             found, 8, &count) == MYB_OK);
  printf("coplanar circles: %zu approaches; %.9f m apart (300 km); a synodic period between "
         "them to %.1e s\n",
         count, found[0].distance, fabs(found[1].t - found[0].t - synodic));
  MT_CHECK(count == 2); /* at 0.94 and 1.94 synodic periods */
  MT_CHECK(fabs(found[0].distance - 3.0e5) < 1e-6);
  MT_CHECK(fabs(found[1].t - found[0].t - synodic) < 1e-6);
  MT_CHECK(fabs(found[0].t - (2.0 * PI - 0.4) / (n1 - n2)) < 1e-6);
  /* more than the capacity: counted all the same */
  MT_CHECK(myb_approach_find(kepler_state, &low, kepler_state, &high, 0.0, 2.5 * synodic, 60.0,
                             found, 1, &count) == MYB_OK);
  MT_CHECK(count == 2);
  MT_CHECK(myb_approach_find(kepler_state, &low, kepler_state, &high, 10.0, 0.0, 60.0, found, 8,
                             &count) == MYB_E_DOMAIN);
  MT_CHECK(myb_approach_find(kepler_state, &low, kepler_state, &high, 0.0, 10.0, 0.0, found, 8,
                             &count) == MYB_E_DOMAIN);

  /* the same pair of crossing orbits propagated: against the providers */
  myb_point_mass_term centre = {MU};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "centre", myb_term_point_mass, myb_term_point_mass_jacobian,
                               &centre, NULL) == MYB_OK);
  myb_orbit* o1 = NULL;
  myb_orbit* o2 = NULL;
  const myb_orbit_options options = {.epsilon = 1e-11};
  MT_CHECK(myb_orbit_create(model, &options, NULL, &o1) == MYB_OK);
  MT_CHECK(myb_orbit_create(model, &options, NULL, &o2) == MYB_OK);
  double y1[6];
  double y2[6];
  memcpy(y1, first.start.r, 3 * sizeof(double));
  memcpy(y1 + 3, first.start.v, 3 * sizeof(double));
  memcpy(y2, second.start.r, 3 * sizeof(double));
  memcpy(y2 + 3, second.start.v, 3 * sizeof(double));
  MT_CHECK(myb_orbit_start(o1, 0.0, y1, 2.0 * period) == MYB_OK);
  MT_CHECK(myb_orbit_start(o2, 0.0, y2, 2.0 * period) == MYB_OK);
  myb_approach marched[8];
  size_t marched_count = 0;
  MT_CHECK(myb_approach_orbits(o1, o2, 1.2 * period, 10.0, marched, 8, &marched_count) ==
           MYB_OK);
  MT_CHECK(myb_approach_find(kepler_state, &first, kepler_state, &second, 0.0, 1.2 * period,
                             30.0, found, 8, &count) == MYB_OK);
  double worst_t = 0.0;
  double worst_d = 0.0;
  for (size_t k = 0; k < count && k < marched_count; ++k) {
    worst_t = fmax(worst_t, fabs(marched[k].t - found[k].t));
    worst_d = fmax(worst_d, fabs(marched[k].distance - found[k].distance));
  }
  printf("two propagated orbits in steps of 10 s against the providers: %zu approaches, time "
         "%.1e s, distance %.1e m\n",
         marched_count, worst_t, worst_d);
  MT_CHECK(marched_count == count);
  MT_CHECK(worst_t < 1e-5 && worst_d < 1e-3);
  MT_CHECK(myb_orbit_time(o1) == 1.2 * period && myb_orbit_time(o2) == 1.2 * period);
  MT_CHECK(myb_approach_orbits(o1, o2, period, 10.0, marched, 8, &marched_count) ==
           MYB_E_DOMAIN); /* behind them */
  myb_orbit_destroy(o1);
  myb_orbit_destroy(o2);
  myb_force_model_destroy(model);
}

/* The sum of the Gaussian over a fine grid on the disk. */
static double grid_probability(const myb_encounter* e, int points) {
  const double h = 2.0 * e->radius / points;
  double sum = 0.0;
  for (int i = 0; i < points; ++i) {
    const double x = -e->radius + (i + 0.5) * h;
    for (int j = 0; j < points; ++j) {
      const double y = -e->radius + (j + 0.5) * h;
      if (x * x + y * y <= e->radius * e->radius) {
        const double dx = (x - e->miss[0]) / e->sigma[0];
        const double dy = (y - e->miss[1]) / e->sigma[1];
        sum += exp(-0.5 * (dx * dx + dy * dy));
      }
    }
  }
  return sum * h * h / (2.0 * PI * e->sigma[0] * e->sigma[1]);
}

static void test_probability(void) {
  double p = 0.0;
  /* a round Gaussian centred on the disk: 1 - exp(-R^2 / 2 sigma^2) */
  double worst = 0.0;
  const double radii[5] = {0.01, 0.3, 1.0, 3.0, 8.0};
  for (int k = 0; k < 5; ++k) {
    const myb_encounter e = {.miss = {0.0, 0.0}, .sigma = {1.0, 1.0}, .radius = radii[k]};
    MT_CHECK(myb_collision_probability(&e, &p) == MYB_OK);
    const double exact = -expm1(-0.5 * radii[k] * radii[k]);
    worst = fmax(worst, fabs(p / exact - 1.0));
  }
  printf("probability, centred round Gaussian: against 1 - exp(-R^2 / 2 sigma^2) %.1e\n", worst);
  MT_CHECK(worst < 1e-13);
  /* offset: the series of Rice, which Chan's is when the sigmas are equal */
  worst = 0.0;
  double smallest = 1.0;
  const double offsets[5] = {0.5, 2.0, 5.0, 12.0, 30.0};
  for (int k = 0; k < 5; ++k) {
    const myb_encounter e = {.miss = {0.6 * offsets[k], 0.8 * offsets[k]},
                             .sigma = {1.0, 1.0}, .radius = 0.7};
    MT_CHECK(myb_collision_probability(&e, &p) == MYB_OK);
    const double series = myb_collision_probability_chan(&e);
    worst = fmax(worst, fabs(p / series - 1.0));
    smallest = fmin(smallest, p);
  }
  printf("  offset round Gaussian, misses of 0.5 to 30 sigma (down to %.1e): against the "
         "series of Rice %.1e\n",
         smallest, worst);
  MT_CHECK(worst < 1e-11 && smallest < 1e-180);
  /* unequal sigmas: a sum over a grid of 4000 x 4000 */
  const myb_encounter flat = {.miss = {1.3, -0.4}, .sigma = {3.0, 0.5}, .radius = 0.9};
  MT_CHECK(myb_collision_probability(&flat, &p) == MYB_OK);
  const double summed = grid_probability(&flat, 4000);
  const double chan = myb_collision_probability_chan(&flat);
  printf("  sigmas 3 and 0.5: %.12f, a grid gives %.12f; Chan's circle of equal area %.6f\n", p,
         summed, chan);
  MT_CHECK(fabs(p / summed - 1.0) < 1e-5);
  MT_CHECK(fabs(chan / p - 1.0) < 0.2); /* an approximation there */
  /* the sign of the miss does not matter, nor which way it is turned by half a turn */
  const myb_encounter mirrored = {.miss = {-1.3, 0.4}, .sigma = {3.0, 0.5}, .radius = 0.9};
  double q = 0.0;
  MT_CHECK(myb_collision_probability(&mirrored, &q) == MYB_OK);
  MT_CHECK(fabs(q / p - 1.0) < 1e-13);
  /* a disk much larger than the Gaussian that sits inside it: certain */
  const myb_encounter inside = {.miss = {1.0, 1.0}, .sigma = {0.5, 0.2}, .radius = 20.0};
  MT_CHECK(myb_collision_probability(&inside, &p) == MYB_OK);
  MT_CHECK(fabs(p - 1.0) < 1e-13);
  const myb_encounter bad = {.miss = {0.0, 0.0}, .sigma = {1.0, 0.0}, .radius = 1.0};
  MT_CHECK(myb_collision_probability(&bad, &p) == MYB_E_DOMAIN);
  MT_CHECK(isnan(myb_collision_probability_chan(&bad)));

  /* the plane of an encounter: a miss of 300 m across a relative velocity along (1, 2, 2),
   * a covariance with the axes 400, 100, 50 m, the long one across the velocity */
  const double v[3] = {1000.0, 2000.0, 2000.0};
  const double across[3] = {2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0}; /* unit, across v */
  const double third[3] = {2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0};  /* across both */
  const double along[3] = {1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0};
  double covariance[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      covariance[3 * i + j] = 400.0 * 400.0 * across[i] * across[j] +
                              100.0 * 100.0 * third[i] * third[j] +
                              50.0 * 50.0 * along[i] * along[j];
    }
  }
  /* the miss 300 m along `across` and 40 m along `third`, plus a part along the velocity
   * that the plane does not see */
  double r[3];
  for (int i = 0; i < 3; ++i) {
    r[i] = 300.0 * across[i] + 40.0 * third[i] + 77.0 * along[i];
  }
  myb_encounter e;
  MT_CHECK(myb_encounter_plane(r, v, covariance, 12.0, &e) == MYB_OK);
  printf("  the plane of an encounter: sigmas %.9f and %.9f m, miss (%.9f, %.9f) m\n",
         e.sigma[0], e.sigma[1], fabs(e.miss[0]), fabs(e.miss[1]));
  MT_CHECK(fabs(e.sigma[0] - 400.0) < 1e-9 && fabs(e.sigma[1] - 100.0) < 1e-9);
  MT_CHECK(fabs(fabs(e.miss[0]) - 300.0) < 1e-9 && fabs(fabs(e.miss[1]) - 40.0) < 1e-9);
  MT_CHECK(e.radius == 12.0);
  MT_CHECK(myb_collision_probability(&e, &p) == MYB_OK);
  /* a disk small against the sigmas: the density at the miss times its area */
  const double density = exp(-0.5 * (0.75 * 0.75 + 0.4 * 0.4)) / (2.0 * PI * 400.0 * 100.0);
  MT_CHECK(fabs(p / (density * PI * 144.0) - 1.0) < 5e-3); /* to (R / sigma)^2 */
  const double zero[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_encounter_plane(r, zero, covariance, 12.0, &e) == MYB_E_DOMAIN);
  covariance[1] += 1.0;
  MT_CHECK(myb_encounter_plane(r, v, covariance, 12.0, &e) == MYB_E_DOMAIN);

  /* the largest probability over the size of a round covariance, for a small disk:
   * R^2 / (e d^2), at sigma = d / sqrt 2 */
  const myb_encounter far = {.miss = {3000.0, 4000.0}, .sigma = {10.0, 10.0}, .radius = 5.0};
  double factor = 0.0;
  MT_CHECK(myb_collision_probability_max(&far, &p, &factor) == MYB_OK);
  const double expected = 25.0 / (exp(1.0) * 5000.0 * 5000.0);
  printf("  the largest probability over the size of the covariance: %.6e (R^2 / e d^2 = "
         "%.6e) at sigma = %.3f m (d / sqrt 2 = %.3f)\n",
         p, expected, factor * 10.0, 5000.0 / sqrt(2.0));
  MT_CHECK(fabs(p / expected - 1.0) < 1e-5);
  MT_CHECK(fabs(factor * 10.0 / (5000.0 / sqrt(2.0)) - 1.0) < 1e-4);
  const myb_encounter through = {.miss = {1.0, 1.0}, .sigma = {10.0, 10.0}, .radius = 5.0};
  MT_CHECK(myb_collision_probability_max(&through, &p, &factor) == MYB_OK);
  MT_CHECK(p == 1.0 && factor == 0.0);
}

/* The covariances carried to the approach by the orbits themselves, against the closed
 * transition matrix of the two-body motion. */
static void test_carried(void) {
  const double a = 7.2e6;
  const double period = 2.0 * PI * sqrt(a * a * a / MU);
  const double lag = 150.0 / a;
  const myb_pv start1 = orbit(a, 0.001, 0.2, 0.3, 0.5, -0.9);
  const myb_pv start2 = orbit(a, 0.001, 1.4, 0.3, 0.5, -0.9 - lag);
  /* 60 m and 6 cm/s along every axis for the first, twice that for the second, with a
   * correlation of the position and the velocity along x */
  double p1[36] = {0};
  double p2[36] = {0};
  for (int i = 0; i < 3; ++i) {
    p1[7 * i] = 60.0 * 60.0;
    p1[7 * (3 + i)] = 0.06 * 0.06;
    p2[7 * i] = 120.0 * 120.0;
    p2[7 * (3 + i)] = 0.12 * 0.12;
  }
  p1[3] = p1[18] = 0.5 * 60.0 * 0.06;
  myb_point_mass_term centre = {MU};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "centre", myb_term_point_mass, myb_term_point_mass_jacobian,
                               &centre, NULL) == MYB_OK);
  const myb_orbit_options options = {.epsilon = 1e-11, .variational = 1};
  myb_orbit* o1 = NULL;
  myb_orbit* o2 = NULL;
  MT_CHECK(myb_orbit_create(model, &options, NULL, &o1) == MYB_OK);
  MT_CHECK(myb_orbit_create(model, &options, NULL, &o2) == MYB_OK);
  double y1[6];
  double y2[6];
  memcpy(y1, start1.r, 3 * sizeof(double));
  memcpy(y1 + 3, start1.v, 3 * sizeof(double));
  memcpy(y2, start2.r, 3 * sizeof(double));
  memcpy(y2 + 3, start2.v, 3 * sizeof(double));
  MT_CHECK(myb_orbit_start(o1, 0.0, y1, 2.0 * period) == MYB_OK);
  MT_CHECK(myb_orbit_start(o2, 0.0, y2, 2.0 * period) == MYB_OK);
  myb_conjunction found[4];
  size_t count = 0;
  MT_CHECK(myb_conjunction_find(o1, p1, o2, p2, 15.0, 1.6 * period, 20.0, found, 4, &count) ==
           MYB_OK);
  MT_CHECK(count >= 2 && count <= 4);
  double worst_covariance = 0.0;
  double worst_probability = 0.0;
  for (size_t k = 0; k < count; ++k) {
    /* the two covariances at that instant in closed form */
    myb_pv pv;
    double phi[36];
    double carried[36];
    double sum[9] = {0};
    MT_CHECK(myb_kepler_stm(MU, &start1, found[k].approach.t, &pv, phi) == MYB_OK);
    myb_covariance_transform(6, phi, p1, carried);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        sum[3 * i + j] += carried[6 * i + j];
      }
    }
    MT_CHECK(myb_kepler_stm(MU, &start2, found[k].approach.t, &pv, phi) == MYB_OK);
    myb_covariance_transform(6, phi, p2, carried);
    double scale = 0.0;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        sum[3 * i + j] += carried[6 * i + j];
        scale = fmax(scale, fabs(sum[3 * i + j]));
      }
    }
    for (int i = 0; i < 9; ++i) {
      worst_covariance = fmax(worst_covariance, fabs(found[k].covariance[i] - sum[i]) / scale);
    }
    myb_encounter plane;
    double probability = 0.0;
    MT_CHECK(myb_encounter_plane(found[k].approach.position, found[k].approach.velocity, sum,
                                 15.0, &plane) == MYB_OK);
    MT_CHECK(myb_collision_probability(&plane, &probability) == MYB_OK);
    worst_probability =
        fmax(worst_probability, fabs(found[k].probability / probability - 1.0));
    MT_CHECK(found[k].probability > 0.0 && found[k].probability < 1.0);
  }
  printf("covariances carried by the orbits to %zu approaches (steps of 20 s): the first "
         "%.1f m apart at %.1f s, sigmas %.1f and %.1f m, probability %.4e; against the closed "
         "transition matrix: covariance %.1e, probability %.1e\n",
         count, found[0].approach.distance, found[0].approach.t, found[0].encounter.sigma[0],
         found[0].encounter.sigma[1], found[0].probability, worst_covariance,
         worst_probability);
  MT_CHECK(worst_covariance < 1e-6);
  MT_CHECK(worst_probability < 1e-5);
  /* the uncertainty grows along the orbit: the later approach is the less certain */
  MT_CHECK(found[1].encounter.sigma[0] > found[0].encounter.sigma[0]);
  myb_orbit_destroy(o1);
  myb_orbit_destroy(o2);
  /* orbits without the transition matrix are refused */
  const myb_orbit_options bare = {.epsilon = 1e-11};
  MT_CHECK(myb_orbit_create(model, &bare, NULL, &o1) == MYB_OK);
  MT_CHECK(myb_orbit_create(model, &bare, NULL, &o2) == MYB_OK);
  MT_CHECK(myb_orbit_start(o1, 0.0, y1, 2.0 * period) == MYB_OK);
  MT_CHECK(myb_orbit_start(o2, 0.0, y2, 2.0 * period) == MYB_OK);
  MT_CHECK(myb_conjunction_find(o1, p1, o2, p2, 15.0, 1.6 * period, 20.0, found, 4, &count) ==
           MYB_E_DOMAIN);
  myb_orbit_destroy(o1);
  myb_orbit_destroy(o2);
  myb_force_model_destroy(model);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_approaches();
  test_probability();
  test_carried();
  return mt_finish("test_conjunction");
}
