/* SPDX-License-Identifier: BSD-3-Clause */
/* Propagation of one object (myb/orbit.h):
 *  - the two-body problem against the universal-variable propagator, forward and backward;
 *  - the state at a time does not depend on what else was asked: bit for bit with and without
 *    500 other samples;
 *  - events: the apsides of an ellipse (roots of x.v) against Kepler's equation, the direction
 *    filter, a search that starts at a root, the limit;
 *  - partials: the state transition matrix against its closed form (myb_kepler_stm);
 *  - partials in a rotating field: the matrix from the analytic gradient is symplectic;
 *  - decay of a circular orbit in a thin atmosphere against da = -2 pi cd A/m rho a^2 per
 *    revolution;
 *  - the order of the times, the step budget, errors. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/myb.h"

static const double mu = 3.986004418e14;
static const double two_pi = 6.283185307179586;

static double apsis(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[0] * y[3] + y[1] * y[4] + y[2] * y[5];
}

static double difference(const double a[6], const double b[6], double* dv) {
  const double dx = hypot(hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2]);
  *dv = hypot(hypot(a[3] - b[3], a[4] - b[4]), a[5] - b[5]);
  return dx;
}

typedef struct problem {
  myb_point_mass_term central;
  myb_force_model* model;
  myb_pv pv0;
  double y0[6];
  double period;
  double mean_motion;
  double mean_anomaly0;
} problem;

static void problem_init(problem* p) {
  p->central.gm = mu;
  MT_CHECK(myb_force_model_create(NULL, &p->model) == MYB_OK);
  MT_CHECK(myb_force_model_add(p->model, "central", myb_term_point_mass,
                               myb_term_point_mass_jacobian, &p->central, NULL) == MYB_OK);
  const double e = 0.3;
  const double a = 1.0e7;
  const myb_elements el = {.p = a * (1.0 - e * e), .e = e, .i = 0.9, .raan = 0.4, .argp = 1.1,
                           .nu = 1.0};
  MT_CHECK(myb_pv_from_elements(mu, &el, &p->pv0) == MYB_OK);
  memcpy(p->y0, p->pv0.r, 3 * sizeof(double));
  memcpy(p->y0 + 3, p->pv0.v, 3 * sizeof(double));
  p->mean_motion = sqrt(mu / (a * a * a));
  p->period = two_pi / p->mean_motion;
  double ecc = 0.0;
  MT_CHECK(myb_anomaly_eccentric_from_true(el.nu, e, &ecc) == MYB_OK);
  MT_CHECK(myb_anomaly_mean_from_eccentric(ecc, e, &p->mean_anomaly0) == MYB_OK);
}

static void kepler(const problem* p, double t, double y[6]) {
  myb_pv pv;
  MT_CHECK(myb_kepler_propagate(mu, &p->pv0, t, &pv) == MYB_OK);
  memcpy(y, pv.r, 3 * sizeof(double));
  memcpy(y + 3, pv.v, 3 * sizeof(double));
}

static void test_two_body(const problem* p) {
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(p->model, NULL, NULL, &orbit) == MYB_OK);
  const double t_end = 5.0 * p->period;
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_time(orbit) == 0.0);
  double y[6];
  double ref[6];
  MT_CHECK(myb_orbit_state(orbit, 0.0, y) == MYB_OK);
  for (int i = 0; i < 6; ++i) {
    MT_CHECK_SAME(y[i], p->y0[i]);
  }
  double worst_x = 0.0;
  double worst_v = 0.0;
  for (int k = 1; k <= 200; ++k) {
    const double t = t_end * k / 200.0;
    MT_CHECK(myb_orbit_state(orbit, t, y) == MYB_OK);
    MT_CHECK(myb_orbit_time(orbit) == t);
    kepler(p, t, ref);
    double dv = 0.0;
    worst_x = fmax(worst_x, difference(y, ref, &dv));
    worst_v = fmax(worst_v, dv);
  }
  printf("two bodies, 5 revolutions (e = 0.3): %.1e m, %.1e m/s from the closed form, %lld steps\n",
         worst_x, worst_v, (long long)myb_orbit_steps(orbit));
  MT_CHECK(worst_x < 1e-4);
  MT_CHECK(worst_v < 1e-7);
  MT_CHECK(myb_orbit_steps(orbit) > 10);

  /* backward */
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, -2.0 * p->period) == MYB_OK);
  MT_CHECK(myb_orbit_steps(orbit) == 0);
  for (int k = 1; k <= 20; ++k) {
    const double t = -2.0 * p->period * k / 20.0;
    MT_CHECK(myb_orbit_state(orbit, t, y) == MYB_OK);
    kepler(p, t, ref);
    double dv = 0.0;
    MT_CHECK(difference(y, ref, &dv) < 1e-4);
    MT_CHECK(dv < 1e-7);
  }

  /* the order of the times */
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_state(orbit, 1000.0, y) == MYB_OK);
  MT_CHECK(myb_orbit_state(orbit, 999.0, y) == MYB_E_DOMAIN);
  MT_CHECK(myb_orbit_state(orbit, t_end + 1.0, y) == MYB_E_DOMAIN);
  MT_CHECK(myb_orbit_state(orbit, (double)NAN, y) == MYB_E_DOMAIN);
  MT_CHECK(myb_orbit_time(orbit) == 1000.0);
  MT_CHECK(myb_orbit_state(orbit, 1000.0, y) == MYB_OK); /* the same time again */
  MT_CHECK(myb_orbit_state(orbit, t_end, y) == MYB_OK);
  MT_CHECK(myb_orbit_partials(orbit, t_end, y, NULL, NULL) == MYB_E_DOMAIN); /* not variational */
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, 0.0) == MYB_E_DOMAIN);
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, (double)INFINITY) == MYB_E_DOMAIN);
  double bad[6];
  memcpy(bad, p->y0, sizeof bad);
  bad[2] = (double)NAN;
  MT_CHECK(myb_orbit_start(orbit, 0.0, bad, t_end) == MYB_E_DOMAIN);
  myb_orbit_destroy(orbit);

  /* the step budget */
  const myb_orbit_options few = {.max_steps = 3};
  MT_CHECK(myb_orbit_create(p->model, &few, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_state(orbit, t_end, y) == MYB_E_NOCONV);
  myb_orbit_destroy(orbit);
  const myb_orbit_options negative = {.epsilon = -1.0};
  MT_CHECK(myb_orbit_create(p->model, &negative, NULL, &orbit) == MYB_E_DOMAIN);
  myb_orbit_destroy(NULL);
}

static void test_sampling(const problem* p) {
  const double t_end = 3.0 * p->period;
  const double probe = 2.37 * p->period;
  myb_orbit* direct = NULL;
  myb_orbit* sampled = NULL;
  MT_CHECK(myb_orbit_create(p->model, NULL, NULL, &direct) == MYB_OK);
  MT_CHECK(myb_orbit_create(p->model, NULL, NULL, &sampled) == MYB_OK);
  MT_CHECK(myb_orbit_start(direct, 0.0, p->y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_start(sampled, 0.0, p->y0, t_end) == MYB_OK);
  double a[6];
  double b[6];
  MT_CHECK(myb_orbit_state(direct, probe, a) == MYB_OK);
  static double times[500];
  static double states[6 * 500];
  for (int k = 0; k < 500; ++k) {
    times[k] = probe * (k + 1) / 501.0;
  }
  MT_CHECK(myb_orbit_states(sampled, times, 500, states) == MYB_OK);
  /* events on the way do not change it either */
  int found = 0;
  double te = 0.0;
  MT_CHECK(myb_orbit_event(sampled, apsis, NULL, probe, 0, &found, &te, NULL) == MYB_OK);
  MT_CHECK(myb_orbit_state(sampled, probe, b) == MYB_OK);
  for (int i = 0; i < 6; ++i) {
    MT_CHECK_SAME(b[i], a[i]);
  }
  MT_CHECK(myb_orbit_steps(direct) == myb_orbit_steps(sampled));
  /* myb_orbit_states is myb_orbit_state at each time */
  MT_CHECK(myb_orbit_start(direct, 0.0, p->y0, t_end) == MYB_OK);
  for (int k = 0; k < 500; k += 37) {
    MT_CHECK(myb_orbit_state(direct, times[k], a) == MYB_OK);
    for (int i = 0; i < 6; ++i) {
      MT_CHECK_SAME(states[6 * k + i], a[i]);
    }
  }
  MT_CHECK(myb_orbit_states(direct, times, 3, states) == MYB_E_DOMAIN); /* behind */
  MT_CHECK(myb_orbit_states(direct, NULL, 0, NULL) == MYB_OK);
  myb_orbit_destroy(direct);
  myb_orbit_destroy(sampled);
}

static void test_events(const problem* p) {
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(p->model, NULL, NULL, &orbit) == MYB_OK);
  const double t_end = 4.0 * p->period;
  /* x.v < 0 before the periapsis and > 0 after it: periapses are roots where g grows */
  const double periapsis = (two_pi - p->mean_anomaly0) / p->mean_motion;
  const double apoapsis = periapsis + 0.5 * p->period - p->period; /* the one before it */
  int found = 0;
  double te = 0.0;
  double ye[6];
  double worst = 0.0;

  /* both kinds in order: apoapsis, periapsis, apoapsis, ... */
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  for (int k = 0; k < 6; ++k) {
    MT_CHECK(myb_orbit_event(orbit, apsis, NULL, t_end, 0, &found, &te, ye) == MYB_OK);
    MT_CHECK(found);
    const double expected = apoapsis + 0.5 * p->period * k;
    worst = fmax(worst, fabs(te - expected));
    MT_CHECK(myb_orbit_time(orbit) == te);
    const double r = hypot(hypot(ye[0], ye[1]), ye[2]);
    const double speed = hypot(hypot(ye[3], ye[4]), ye[5]);
    MT_CHECK(fabs(apsis(NULL, te, ye)) < 1e-9 * r * speed);
  }
  /* only the periapses */
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  for (int k = 0; k < 3; ++k) {
    MT_CHECK(myb_orbit_event(orbit, apsis, NULL, t_end, 1, &found, &te, NULL) == MYB_OK);
    MT_CHECK(found);
    worst = fmax(worst, fabs(te - (periapsis + p->period * k)));
  }
  /* only the apoapses, with states asked in between */
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  for (int k = 0; k < 3; ++k) {
    MT_CHECK(myb_orbit_event(orbit, apsis, NULL, t_end, -1, &found, &te, NULL) == MYB_OK);
    MT_CHECK(found);
    worst = fmax(worst, fabs(te - (apoapsis + p->period * k)));
    MT_CHECK(myb_orbit_state(orbit, te + 100.0, ye) == MYB_OK);
  }
  printf("events: apsides of 4 revolutions against Kepler's equation %.1e s\n", worst);
  MT_CHECK(worst < 1e-7);

  /* the limit: nothing before it, then the root just behind it */
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_event(orbit, apsis, NULL, apoapsis - 1.0, 0, &found, &te, NULL) == MYB_OK);
  MT_CHECK(!found);
  MT_CHECK(myb_orbit_time(orbit) == apoapsis - 1.0);
  MT_CHECK(myb_orbit_event(orbit, apsis, NULL, apoapsis + 1.0, 0, &found, &te, NULL) == MYB_OK);
  MT_CHECK(found && fabs(te - apoapsis) < 1e-7);
  MT_CHECK(myb_orbit_event(orbit, apsis, NULL, apoapsis - 0.5, 0, &found, &te, NULL) ==
           MYB_E_DOMAIN);

  /* backward: the periapsis before the start, where g still grows with t */
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, -p->period) == MYB_OK);
  MT_CHECK(myb_orbit_event(orbit, apsis, NULL, -p->period, 1, &found, &te, NULL) == MYB_OK);
  MT_CHECK(found && fabs(te - (periapsis - p->period)) < 1e-7);
  MT_CHECK(myb_orbit_event(orbit, apsis, NULL, -p->period, 1, &found, &te, NULL) == MYB_OK);
  MT_CHECK(!found);
  myb_orbit_destroy(orbit);
}

static void test_partials(const problem* p) {
  const myb_orbit_options options = {.variational = 1, .epsilon = 1e-11};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(p->model, &options, NULL, &orbit) == MYB_OK);
  const double t_end = 2.0 * p->period;
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  double y[6];
  double phi[36];
  MT_CHECK(myb_orbit_partials(orbit, 0.0, y, phi, NULL) == MYB_OK);
  for (int i = 0; i < 36; ++i) {
    MT_CHECK_SAME(phi[i], i % 7 == 0 ? 1.0 : 0.0);
  }
  double worst = 0.0;
  for (int k = 1; k <= 8; ++k) {
    const double t = t_end * k / 8.0;
    MT_CHECK(myb_orbit_partials(orbit, t, y, phi, NULL) == MYB_OK);
    myb_pv pv;
    double exact[36];
    MT_CHECK(myb_kepler_stm(mu, &p->pv0, t, &pv, exact) == MYB_OK);
    double scale = 0.0;
    double error = 0.0;
    for (int i = 0; i < 36; ++i) {
      scale = fmax(scale, fabs(exact[i]));
      error = fmax(error, fabs(phi[i] - exact[i]));
    }
    worst = fmax(worst, error / scale);
    double ys[6];
    MT_CHECK(myb_orbit_state(orbit, t, ys) == MYB_OK); /* states of a variational object too */
    for (int i = 0; i < 6; ++i) {
      MT_CHECK_SAME(ys[i], y[i]);
    }
  }
  printf("partials: state transition matrix over 2 revolutions against the closed form %.1e\n",
         worst);
  MT_CHECK(worst < 1e-11);
  myb_orbit_destroy(orbit);
}

/* A circular orbit in an atmosphere of constant scale height on a planet that does not rotate:
 * to first order in the drag the semi-major axis falls by 2 pi cd A/m rho a^2 per revolution
 * (D. King-Hele, Satellite Orbits in an Atmosphere, 1987, eq. 4.9 for e = 0). */
static void test_decay(void) {
  const myb_time epoch = myb_time_from_seconds(8.0e8);
  const myb_env_options options = {.epoch = epoch};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  size_t planet = 0;
  const myb_env_body b_planet = {.gm = mu};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  myb_exponential_atmosphere air = {.ellipsoid = {.a = 6378137.0, .f = 0.0},
                                    .count = 1,
                                    .base = {400000.0},
                                    .density = {4.0e-12},
                                    .scale = {60000.0}};
  myb_object object = {.mass = 1000.0, .area = 1.0, .cr = 1.0, .cd = 2.2};
  myb_point_mass_term central = {.gm = mu};
  myb_drag_term drag = {.env = env, .body = planet, .density = myb_density_exponential,
                        .density_context = &air, .object = &object};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "central", myb_term_point_mass,
                               myb_term_point_mass_jacobian, &central, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "drag", myb_term_drag, NULL, &drag, NULL) == MYB_OK);
  const double a0 = 6378137.0 + 400000.0;
  const double speed = sqrt(mu / a0);
  const double y0[6] = {a0, 0.0, 0.0, 0.0, speed, 0.0};
  const double period = two_pi * sqrt(a0 * a0 * a0 / mu);
  const int revolutions = 10;
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, NULL, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, revolutions * period) == MYB_OK);
  double y[6];
  MT_CHECK(myb_orbit_state(orbit, revolutions * period, y) == MYB_OK);
  const double r = hypot(hypot(y[0], y[1]), y[2]);
  const double v2 = y[3] * y[3] + y[4] * y[4] + y[5] * y[5];
  const double a1 = 1.0 / (2.0 / r - v2 / mu);
  const double expected =
      -two_pi * object.cd * object.area / object.mass * air.density[0] * a0 * a0 * revolutions;
  printf("drag: decay of a circular orbit over %d revolutions %.4f m, first order %.4f m\n",
         revolutions, a1 - a0, expected);
  MT_CHECK(fabs((a1 - a0) / expected - 1.0) < 2e-3);
  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);
  myb_env_destroy(env);
}

typedef struct turning {
  double rate;
} turning;

static myb_status turning_orientation(void* context, myb_time t, myb_orientation* o) {
  const turning* s = context;
  const double angles[3] = {0.3, 0.2, 1.0 + s->rate * myb_time_seconds(t)};
  const double rates[3] = {0.0, 0.0, s->rate};
  myb_rotation_euler313(angles, rates, &o->rotation);
  o->ut = t;
  o->xp = 0.0;
  o->yp = 0.0;
  return MYB_OK;
}

/* The state transition matrix in the field of a rotating oblate body, from the analytic
 * gradient of the field: symplectic (the motion is Hamiltonian, also with the rotation),
 * equal to differences of whole propagations, and far closer to both than the matrix from a
 * Jacobian by differences. */
static void test_field_partials(void) {
  const myb_env_options env_options = {.epoch = myb_time_from_seconds(0.0)};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&env_options, NULL, &env) == MYB_OK);
  turning spin = {.rate = 7.292115e-5};
  size_t planet = 0;
  const myb_env_body body = {.gm = mu, .orientation = turning_orientation,
                             .orientation_context = &spin};
  MT_CHECK(myb_env_add_body(env, &body, &planet) == MYB_OK);
  double c[45] = {0.0};
  double sn[45] = {0.0};
  c[0] = 1.0;
  mt_rng rng = {.state = 5};
  for (int n = 2; n <= 8; ++n) {
    for (int m = 0; m <= n; ++m) {
      const double size = (n == 2 && m == 0 ? 4.8e-4 : 2e-6) / (n - 1);
      c[n * (n + 1) / 2 + m] = size * ((double)(mt_next(&rng) >> 11) / 4503599627370496.0 - 1.0);
      sn[n * (n + 1) / 2 + m] =
          m == 0 ? 0.0 : size * ((double)(mt_next(&rng) >> 11) / 4503599627370496.0 - 1.0);
    }
  }
  myb_gravity_field* field = NULL;
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_field_create(mu, 6378136.3, 8, c, sn, NULL, &field) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(8, NULL, &work) == MYB_OK);
  myb_field_term gravity = {.env = env, .body = planet, .field = field, .work = work,
                            .n_min = 0, .n_max = -1};
  myb_force_model* analytic = NULL;
  myb_force_model* numeric = NULL;
  MT_CHECK(myb_force_model_create(NULL, &analytic) == MYB_OK);
  MT_CHECK(myb_force_model_create(NULL, &numeric) == MYB_OK);
  MT_CHECK(myb_force_model_add(analytic, "field", myb_term_field, myb_term_field_jacobian,
                               &gravity, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(numeric, "field", myb_term_field, NULL, &gravity, NULL) == MYB_OK);
  const double a = 7.2e6;
  const myb_elements el = {.p = a * (1.0 - 0.01 * 0.01), .e = 0.01, .i = 1.1, .raan = 0.4,
                           .argp = 1.1, .nu = 0.3};
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(mu, &el, &pv) == MYB_OK);
  double y0[6];
  memcpy(y0, pv.r, 3 * sizeof(double));
  memcpy(y0 + 3, pv.v, 3 * sizeof(double));
  const double n_mean = sqrt(mu / (a * a * a));
  const double t_end = 2.0 * two_pi / n_mean;
  const myb_orbit_options options = {.variational = 1, .epsilon = 1e-11};
  double y[6];
  double phi[36];
  double phi_numeric[36];
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(analytic, &options, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_partials(orbit, t_end, y, phi, NULL) == MYB_OK);
  myb_orbit_destroy(orbit);
  MT_CHECK(myb_orbit_create(numeric, &options, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_partials(orbit, t_end, y, phi_numeric, NULL) == MYB_OK);
  myb_orbit_destroy(orbit);
  /* in the units of the orbit (a, 1 / n) every block of Phi is a pure number */
  double scaled[36];
  double scaled_numeric[36];
  double largest = 0.0;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      const double unit = (i < 3 ? 1.0 : 1.0 / n_mean) / (j < 3 ? 1.0 : 1.0 / n_mean);
      scaled[6 * i + j] = phi[6 * i + j] * unit;
      scaled_numeric[6 * i + j] = phi_numeric[6 * i + j] * unit;
      largest = fmax(largest, fabs(scaled[6 * i + j]));
    }
  }
  /* Phi^T J Phi = J */
  double defect = 0.0;
  double defect_numeric = 0.0;
  for (int pass = 0; pass < 2; ++pass) {
    const double* f = pass == 0 ? scaled : scaled_numeric;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        double sum = 0.0;
        for (int k = 0; k < 3; ++k) {
          sum += f[6 * k + i] * f[6 * (k + 3) + j] - f[6 * (k + 3) + i] * f[6 * k + j];
        }
        const double expected = (j == i + 3) ? 1.0 : (i == j + 3 ? -1.0 : 0.0);
        const double d = fabs(sum - expected) / (largest * largest);
        if (pass == 0) {
          defect = fmax(defect, d);
        } else {
          defect_numeric = fmax(defect_numeric, d);
        }
      }
    }
  }
  /* differences of whole propagations */
  myb_orbit_options plain = {.epsilon = 1e-11};
  MT_CHECK(myb_orbit_create(analytic, &plain, NULL, &orbit) == MYB_OK);
  double worst = 0.0;
  for (int j = 0; j < 6; ++j) {
    const double h = j < 3 ? 1.0 : 1.0e-3;
    double plus[6];
    double minus[6];
    double yp[6];
    double ym[6];
    memcpy(plus, y0, sizeof plus);
    memcpy(minus, y0, sizeof minus);
    plus[j] += h;
    minus[j] -= h;
    MT_CHECK(myb_orbit_start(orbit, 0.0, plus, t_end) == MYB_OK);
    MT_CHECK(myb_orbit_state(orbit, t_end, yp) == MYB_OK);
    MT_CHECK(myb_orbit_start(orbit, 0.0, minus, t_end) == MYB_OK);
    MT_CHECK(myb_orbit_state(orbit, t_end, ym) == MYB_OK);
    for (int i = 0; i < 6; ++i) {
      const double unit = (i < 3 ? 1.0 : 1.0 / n_mean) / (j < 3 ? 1.0 : 1.0 / n_mean);
      const double numeric_column = (yp[i] - ym[i]) / (2.0 * h) * unit;
      worst = fmax(worst, fabs(scaled[6 * i + j] - numeric_column) / largest);
    }
  }
  myb_orbit_destroy(orbit);
  printf("partials in a rotating field of degree 8, 2 revolutions: symplectic to %.1e (%.1e "
         "with the Jacobian by differences); against differences of propagations %.1e\n",
         defect, defect_numeric, worst);
  MT_CHECK(defect < 1e-12);
  MT_CHECK(defect_numeric > 10.0 * defect);
  MT_CHECK(worst < 1e-6);
  myb_force_model_destroy(analytic);
  myb_force_model_destroy(numeric);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
  myb_env_destroy(env);
}

/* The other integrator: the same answers through the same interface. */
static void test_dop853(const problem* p) {
  const myb_orbit_options options = {.integrator = MYB_ORBIT_DOP853};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(p->model, &options, NULL, &orbit) == MYB_OK);
  const double t_end = 5.0 * p->period;
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  double y[6];
  double ref[6];
  double worst_x = 0.0;
  double worst_v = 0.0;
  for (int k = 1; k <= 200; ++k) {
    const double t = t_end * k / 200.0;
    MT_CHECK(myb_orbit_state(orbit, t, y) == MYB_OK);
    kepler(p, t, ref);
    double dv = 0.0;
    worst_x = fmax(worst_x, difference(y, ref, &dv));
    worst_v = fmax(worst_v, dv);
  }
  printf("DOP853, two bodies, 5 revolutions: %.1e m, %.1e m/s from the closed form, %lld "
         "steps\n",
         worst_x, worst_v, (long long)myb_orbit_steps(orbit));
  MT_CHECK(worst_x < 1e-2 && worst_v < 1e-5);
  /* the state at a time is the same whatever else was asked */
  double alone[6];
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_state(orbit, 0.37 * t_end, alone) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  for (int k = 1; k <= 300; ++k) {
    MT_CHECK(myb_orbit_state(orbit, 0.37 * t_end * k / 300.0, y) == MYB_OK);
  }
  MT_CHECK(memcmp(y, alone, sizeof y) == 0);
  /* events: the next pericentre */
  int found = 0;
  double t_event = 0.0;
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_event(orbit, apsis, NULL, t_end, +1, &found, &t_event, y) == MYB_OK);
  const double expected = (two_pi - p->mean_anomaly0) / p->mean_motion;
  MT_CHECK(found && fabs(t_event - expected) < 1e-6);
  /* an impulse, then on */
  const double dv[3] = {0.0, 5.0, 0.0};
  MT_CHECK(myb_orbit_kick(orbit, dv) == MYB_OK);
  double kicked[6];
  MT_CHECK(myb_orbit_state(orbit, t_event, kicked) == MYB_OK);
  MT_CHECK(kicked[4] == y[4] + 5.0 && kicked[0] == y[0]);
  MT_CHECK(myb_orbit_state(orbit, t_end, y) == MYB_OK);
  /* no variational equations with it */
  MT_CHECK(myb_orbit_partials(orbit, t_end, y, NULL, NULL) == MYB_E_DOMAIN);
  myb_orbit_destroy(orbit);
  const myb_orbit_options both = {.integrator = MYB_ORBIT_DOP853, .variational = 1};
  MT_CHECK(myb_orbit_create(p->model, &both, NULL, &orbit) == MYB_E_DOMAIN);
  const myb_orbit_options unknown = {.integrator = (myb_orbit_integrator)7};
  MT_CHECK(myb_orbit_create(p->model, &unknown, NULL, &orbit) == MYB_E_DOMAIN);
  /* a looser tolerance: fewer steps, a larger error */
  const myb_orbit_options loose = {.integrator = MYB_ORBIT_DOP853, .rtol = 1e-8, .atol = 1e-5};
  MT_CHECK(myb_orbit_create(p->model, &loose, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, p->y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_state(orbit, t_end, y) == MYB_OK);
  kepler(p, t_end, ref);
  double dv_loose = 0.0;
  const double dx_loose = difference(y, ref, &dv_loose);
  MT_CHECK(dx_loose > worst_x && dx_loose < 100.0);
  myb_orbit_destroy(orbit);
}

int main(void) {
  problem p;
  problem_init(&p);
  test_two_body(&p);
  test_sampling(&p);
  test_events(&p);
  test_partials(&p);
  test_field_partials();
  test_decay();
  test_dop853(&p);
  myb_force_model_destroy(p.model);
  return mt_finish("test_orbit");
}
