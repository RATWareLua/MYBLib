/* SPDX-License-Identifier: BSD-3-Clause */
/* Events and breaks (myb/events.h, myb_force_model_add_break):
 *  - breaks at the edges of a shadow: with radiation pressure and the Earth's shadow the
 *    final state changes with the precision of the integrator as long as steps straddle the
 *    edges, and stops doing so once the edges are declared;
 *  - the transition matrix and the sensitivities carried across the breaks, against
 *    differences of whole propagations;
 *  - ready-made events on a Kepler orbit against closed forms: apsides, a node, a distance;
 *    a close approach to a moving body; an altitude and an elevation against the geodetic
 *    and the station routines; the impact on a shape;
 *  - eclipses against CSPICE N0067 (gfoclt_c on an SPK written by this test; spiceypy,
 *    tests/data/events/cspice.ref), and with the same kernels the surface geometry of
 *    myb/geometry.h against sincpt_c, subpnt_c and ilumin_c.
 * MYB_EVENTS_DUMP=<path> writes the SPK the oracle reads. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/myb.h"

static const double mu = 3.986004418e14;
static const double pi = 3.141592653589793;
static const double two_pi = 6.283185307179586;
static const double au = 1.495978707e11;
static const double earth_radius = 6378136.6;
static const double sun_radius = 6.957e8;

/* ---- The world: the Earth at the origin, turning; the Sun on a circle around it ----------- */

typedef struct mover {
  double radius;
  double rate;
  double phase;
  double inclination;
} mover;

static myb_status mover_state(void* context, myb_time t, double y[6]) {
  const mover* m = context;
  const double u = m->phase + m->rate * myb_time_seconds(t);
  const double cu = myb_cos(u);
  const double su = myb_sin(u);
  const double ci = myb_cos(m->inclination);
  const double si = myb_sin(m->inclination);
  y[0] = m->radius * cu;
  y[1] = m->radius * su * ci;
  y[2] = m->radius * su * si;
  y[3] = -m->radius * m->rate * su;
  y[4] = m->radius * m->rate * cu * ci;
  y[5] = m->radius * m->rate * cu * si;
  return MYB_OK;
}

/* Body-fixed = R3(angle0 + rate t) inertial: what a PCK with the pole along z gives. */
typedef struct spinner {
  double angle0;
  double rate;
} spinner;

static void spin_matrix(const spinner* s, double t, myb_rotation* r) {
  const double w = s->angle0 + s->rate * t;
  const double c = myb_cos(w);
  const double sn = myb_sin(w);
  const double m[9] = {c, sn, 0.0, -sn, c, 0.0, 0.0, 0.0, 1.0};
  const double dm[9] = {-sn * s->rate, c * s->rate, 0.0, -c * s->rate, -sn * s->rate, 0.0,
                        0.0,           0.0,         0.0};
  memcpy(r->m, m, sizeof m);
  memcpy(r->dm, dm, sizeof dm);
}

static myb_status spinner_orientation(void* context, myb_time t, myb_orientation* o) {
  spin_matrix(context, myb_time_seconds(t), &o->rotation);
  o->ut = t;
  o->xp = 0.0;
  o->yp = 0.0;
  return MYB_OK;
}

typedef struct world {
  mover sun_motion;
  mover moon_motion;
  spinner spin;
  myb_env* env;
  size_t earth;
  size_t sun;
  size_t moon;
  myb_object object;
  myb_point_mass_term central;
  myb_radiation_term light;
  myb_eclipse_event edges[2];
  myb_force_model* model; /* the centre and the light */
  double y0[6];
  double period;
} world;

static void world_init(world* w, int with_light, int with_breaks) {
  memset(w, 0, sizeof *w);
  w->sun_motion = (mover){.radius = au, .rate = 1.99e-7, .phase = 0.3, .inclination = 0.4};
  w->moon_motion = (mover){.radius = 3.844e8, .rate = 2.66e-6, .phase = 2.0, .inclination = 0.3};
  /* W = 30 deg + 360.9856235 deg/day, and the 90 deg of the IAU convention */
  w->spin = (spinner){.angle0 = (90.0 + 30.0) * pi / 180.0,
                      .rate = 360.9856235 * pi / 180.0 / 86400.0};
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &w->env) == MYB_OK);
  const myb_env_body earth = {.id = 399, .gm = mu, .orientation = spinner_orientation,
                              .orientation_context = &w->spin, .radius = earth_radius};
  const myb_env_body sun = {.id = 10, .gm = 1.32712440018e20, .state = mover_state,
                            .state_context = &w->sun_motion, .radius = sun_radius,
                            .luminosity = 3.828e26};
  const myb_env_body moon = {.id = 301, .gm = 4.9028e12, .state = mover_state,
                             .state_context = &w->moon_motion, .radius = 1.7374e6};
  MT_CHECK(myb_env_add_body(w->env, &earth, &w->earth) == MYB_OK);
  MT_CHECK(myb_env_add_body(w->env, &sun, &w->sun) == MYB_OK);
  MT_CHECK(myb_env_add_body(w->env, &moon, &w->moon) == MYB_OK);
  w->object = (myb_object){.mass = 100.0, .area = 20.0, .cr = 1.3, .cd = 2.2};
  w->central.gm = mu;
  w->light = (myb_radiation_term){.env = w->env, .source = w->sun, .object = &w->object,
                                  .c = 299792458.0, .occulters = {w->earth},
                                  .occulter_count = 1};
  MT_CHECK(myb_force_model_create(NULL, &w->model) == MYB_OK);
  MT_CHECK(myb_force_model_add(w->model, "central", myb_term_point_mass,
                               myb_term_point_mass_jacobian, &w->central, NULL) == MYB_OK);
  if (with_light) {
    MT_CHECK(myb_force_model_add(w->model, "light", myb_term_radiation,
                                 myb_term_radiation_jacobian, &w->light, NULL) == MYB_OK);
  }
  if (with_breaks) {
    MT_CHECK(myb_term_radiation_breaks(&w->light, w->model, w->edges) == MYB_OK);
    MT_CHECK(myb_force_model_break_count(w->model) == 2);
  }
  const double a = 1.2e7;
  const myb_elements el = {.p = a * (1.0 - 0.01 * 0.01), .e = 0.01, .i = 0.2, .raan = 0.3,
                           .argp = 1.0, .nu = 0.5};
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(mu, &el, &pv) == MYB_OK);
  memcpy(w->y0, pv.r, 3 * sizeof(double));
  memcpy(w->y0 + 3, pv.v, 3 * sizeof(double));
  w->period = two_pi * sqrt(a * a * a / mu);
}

static void world_destroy(world* w) {
  myb_force_model_destroy(w->model);
  myb_env_destroy(w->env);
}

static double distance6(const double a[6], const double b[6]) {
  return hypot(hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2]);
}

/* The state after `span` at the precision epsilon. */
static void run(world* w, double epsilon, double span, double y[6], int64_t* steps,
                int64_t* breaks) {
  const myb_orbit_options options = {.epsilon = epsilon};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(w->model, &options, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, w->y0, span) == MYB_OK);
  MT_CHECK(myb_orbit_state(orbit, span, y) == MYB_OK);
  if (steps != NULL) {
    *steps = myb_orbit_steps(orbit);
  }
  if (breaks != NULL) {
    *breaks = myb_orbit_breaks(orbit);
  }
  myb_orbit_destroy(orbit);
}

static void test_breaks(void) {
  static const double epsilons[4] = {1e-9, 1e-10, 1e-11, 1e-12};
  const double span = 2.0 * 86400.0;
  double spread[2] = {0.0, 0.0};
  double truth[6];
  for (int with = 1; with >= 0; --with) {
    world w;
    world_init(&w, 1, with);
    if (with) {
      int64_t breaks = 0;
      run(&w, 1e-13, span, truth, NULL, &breaks);
      /* two revolutions in three hours and a half: four edges each */
      const double revolutions = span / w.period;
      printf("%.1f revolutions, %lld breaks\n", revolutions, (long long)breaks);
      MT_CHECK(breaks >= 4 * (int64_t)revolutions - 4 && breaks <= 4 * (int64_t)revolutions + 4);
    }
    for (int k = 0; k < 4; ++k) {
      double y[6];
      int64_t steps = 0;
      run(&w, epsilons[k], span, y, &steps, NULL);
      const double d = distance6(y, truth);
      printf("  %s breaks, epsilon %.0e: %lld steps, %.3e m from the reference\n",
             with ? "with" : "without", epsilons[k], (long long)steps, d);
      spread[with] = fmax(spread[with], d);
    }
    world_destroy(&w);
  }
  printf("the final state over the precisions: %.2e m without breaks, %.2e m with\n", spread[0],
         spread[1]);
  /* without the breaks the error does not fall with epsilon: it is that of the steps over
   * the edges; with them it does */
  MT_CHECK(spread[0] > 1.0e-4);
  MT_CHECK(spread[1] < 2.0e-6);
  MT_CHECK(spread[1] < spread[0] / 100.0);
}

/* The same edges with the other integrator. */
static void test_breaks_dop853(void) {
  const double span = 2.0 * 86400.0;
  double reference[6];
  double spread[2] = {0.0, 0.0};
  for (int with = 1; with >= 0; --with) {
    world w;
    world_init(&w, 1, with);
    if (with) {
      run(&w, 1e-13, span, reference, NULL, NULL); /* IAS15 */
    }
    static const double tolerances[3] = {1e-10, 1e-11, 1e-12};
    for (int k = 0; k < 3; ++k) {
      const myb_orbit_options options = {.integrator = MYB_ORBIT_DOP853, .rtol = tolerances[k],
                                         .atol = 1e-3 * tolerances[k]};
      myb_orbit* orbit = NULL;
      double y[6];
      MT_CHECK(myb_orbit_create(w.model, &options, NULL, &orbit) == MYB_OK);
      MT_CHECK(myb_orbit_start(orbit, 0.0, w.y0, span) == MYB_OK);
      MT_CHECK(myb_orbit_state(orbit, span, y) == MYB_OK);
      const double d = distance6(y, reference);
      spread[with] = fmax(spread[with], d);
      if (with && k == 2) {
        printf("DOP853 through the same shadow: %lld steps, %lld breaks, %.2e m from IAS15\n",
               (long long)myb_orbit_steps(orbit), (long long)myb_orbit_breaks(orbit), d);
        MT_CHECK(myb_orbit_breaks(orbit) >= 48 && myb_orbit_breaks(orbit) <= 56);
      }
      myb_orbit_destroy(orbit);
    }
    world_destroy(&w);
  }
  printf("  over its tolerances: %.2e m without breaks, %.2e m with\n", spread[0], spread[1]);
  /* millimetres are the accuracy of DOP853 at these tolerances over two days */
  MT_CHECK(spread[1] < 1e-2);
  MT_CHECK(spread[1] < spread[0] / 1000.0);
}

/* Phi and the sensitivity to cr across the breaks against differences of propagations. */
static void test_partials(void) {
  world w;
  world_init(&w, 1, 1);
  MT_CHECK(myb_force_model_add_parameter(w.model, "cr", &w.object.cr, 1e-3, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_bind_parameter(w.model, 0, 1) == MYB_OK);
  const double span = 1.5 * w.period; /* six edges */
  const myb_orbit_options options = {.epsilon = 1e-11, .variational = 1};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(w.model, &options, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, w.y0, span) == MYB_OK);
  double y[6];
  double phi[36];
  double sens[6];
  /* a time inside a later stretch, then the end */
  double y_mid[6];
  double phi_mid[36];
  MT_CHECK(myb_orbit_partials(orbit, 0.7 * span, y_mid, phi_mid, NULL) == MYB_OK);
  MT_CHECK(myb_orbit_partials(orbit, span, y, phi, sens) == MYB_OK);
  MT_CHECK(myb_orbit_breaks(orbit) >= 5);
  const int64_t steps_with = myb_orbit_steps(orbit);
  const int64_t breaks_met = myb_orbit_breaks(orbit);
  myb_orbit_destroy(orbit);
  /* the same without the breaks: what they cost in steps */
  {
    world plain_world;
    world_init(&plain_world, 1, 0);
    MT_CHECK(myb_force_model_add_parameter(plain_world.model, "cr", &plain_world.object.cr, 1e-3,
                                           NULL) == MYB_OK);
    MT_CHECK(myb_force_model_bind_parameter(plain_world.model, 0, 1) == MYB_OK);
    MT_CHECK(myb_orbit_create(plain_world.model, &options, NULL, &orbit) == MYB_OK);
    MT_CHECK(myb_orbit_start(orbit, 0.0, plain_world.y0, span) == MYB_OK);
    double unused[6];
    MT_CHECK(myb_orbit_partials(orbit, span, unused, NULL, NULL) == MYB_OK);
    printf("partials: %lld breaks in 1.5 revolutions; %lld steps with them, %lld without\n",
           (long long)breaks_met, (long long)steps_with, (long long)myb_orbit_steps(orbit));
    MT_CHECK(steps_with < 2 * myb_orbit_steps(orbit));
    myb_orbit_destroy(orbit);
    world_destroy(&plain_world);
  }
  /* the state is that of the plain propagation with the same breaks */
  double plain[6];
  run(&w, 1e-11, span, plain, NULL, NULL);
  MT_CHECK(distance6(plain, y) < 1e-6);
  /* central differences of the initial state */
  double worst = 0.0;
  double worst_mid = 0.0;
  for (int j = 0; j < 6; ++j) {
    const double h = j < 3 ? 10.0 : 1.0e-2;
    double plus[6];
    double minus[6];
    double plus_mid[6];
    double minus_mid[6];
    const double saved = w.y0[j];
    w.y0[j] = saved + h;
    run(&w, 1e-12, span, plus, NULL, NULL);
    run(&w, 1e-12, 0.7 * span, plus_mid, NULL, NULL);
    w.y0[j] = saved - h;
    run(&w, 1e-12, span, minus, NULL, NULL);
    run(&w, 1e-12, 0.7 * span, minus_mid, NULL, NULL);
    w.y0[j] = saved;
    for (int i = 0; i < 6; ++i) {
      /* rows of positions in m per unit, of velocities in m/s per unit: compare each with
       * the size of its row and column */
      const double scale = (i < 3 ? 1.0 : 1.0e-3) * (j < 3 ? 1.0 : 1.0e3);
      worst = fmax(worst, fabs((plus[i] - minus[i]) / (2.0 * h) - phi[6 * i + j]) / scale);
      worst_mid = fmax(worst_mid, fabs((plus_mid[i] - minus_mid[i]) / (2.0 * h) -
                                       phi_mid[6 * i + j]) / scale);
    }
  }
  double largest = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      largest = fmax(largest, fabs(phi[6 * i + j]));
    }
  }
  /* the sensitivity to cr */
  const double cr = w.object.cr;
  double plus[6];
  double minus[6];
  w.object.cr = cr + 0.05;
  run(&w, 1e-12, span, plus, NULL, NULL);
  w.object.cr = cr - 0.05;
  run(&w, 1e-12, span, minus, NULL, NULL);
  w.object.cr = cr;
  double worst_s = 0.0;
  double largest_s = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst_s = fmax(worst_s, fabs((plus[i] - minus[i]) / 0.1 - sens[i]));
    largest_s = fmax(largest_s, fabs(sens[i]));
  }
  printf("across the breaks: Phi within %.1e of differences (largest element %.1f), at 0.7 of "
         "the arc %.1e; d x / d cr within %.1e m of %.1f m\n",
         worst, largest, worst_mid, worst_s, largest_s);
  MT_CHECK(worst < 2.0e-6 * largest && worst_mid < 2.0e-6 * largest);
  MT_CHECK(worst_s < 1.0e-5 * largest_s);
  world_destroy(&w);
}

/* ---- Events on a Kepler orbit --------------------------------------------------------------- */

typedef struct kepler_orbit {
  myb_force_model* model;
  myb_point_mass_term central;
  myb_elements el;
  double a;
  double n;
  double m0;
  double y0[6];
} kepler_orbit;

static void kepler_init(kepler_orbit* k) {
  k->central.gm = mu;
  MT_CHECK(myb_force_model_create(NULL, &k->model) == MYB_OK);
  MT_CHECK(myb_force_model_add(k->model, "central", myb_term_point_mass,
                               myb_term_point_mass_jacobian, &k->central, NULL) == MYB_OK);
  k->a = 1.0e7;
  k->el = (myb_elements){.p = k->a * (1.0 - 0.3 * 0.3), .e = 0.3, .i = 0.9, .raan = 0.4,
                         .argp = 1.1, .nu = 1.0};
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(mu, &k->el, &pv) == MYB_OK);
  memcpy(k->y0, pv.r, 3 * sizeof(double));
  memcpy(k->y0 + 3, pv.v, 3 * sizeof(double));
  k->n = sqrt(mu / (k->a * k->a * k->a));
  double ecc = 0.0;
  MT_CHECK(myb_anomaly_eccentric_from_true(k->el.nu, k->el.e, &ecc) == MYB_OK);
  MT_CHECK(myb_anomaly_mean_from_eccentric(ecc, k->el.e, &k->m0) == MYB_OK);
}

/* The first time after 0 at which the true anomaly is nu. */
static double time_of(const kepler_orbit* k, double nu) {
  double ecc = 0.0;
  double mean = 0.0;
  MT_CHECK(myb_anomaly_eccentric_from_true(nu, k->el.e, &ecc) == MYB_OK);
  MT_CHECK(myb_anomaly_mean_from_eccentric(ecc, k->el.e, &mean) == MYB_OK);
  double t = (mean - k->m0) / k->n;
  const double period = two_pi / k->n;
  /* bound: a few periods */
  while (t <= 0.0) {
    t += period;
  }
  while (t > period) {
    t -= period;
  }
  return t;
}

static double first_event(const kepler_orbit* k, myb_event_function g, void* context,
                          int direction, double limit, double y[6]) {
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(k->model, NULL, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, k->y0, limit) == MYB_OK);
  int found = 0;
  double t = (double)NAN;
  MT_CHECK(myb_orbit_event(orbit, g, context, limit, direction, &found, &t, y) == MYB_OK);
  MT_CHECK(found);
  myb_orbit_destroy(orbit);
  return found ? t : (double)NAN;
}

static void test_kepler_events(void) {
  kepler_orbit k;
  kepler_init(&k);
  const double period = two_pi / k.n;
  double y[6];
  /* apsides: the rate of the distance to the origin */
  myb_range_event centre = {.body = MYB_EVENT_ORIGIN};
  const double t_peri = first_event(&k, myb_event_range_rate, &centre, +1, 2.0 * period, y);
  const double t_apo = first_event(&k, myb_event_range_rate, &centre, -1, 2.0 * period, y);
  printf("pericentre %.2e s, apocentre %.2e s from Kepler's equation\n",
         fabs(t_peri - time_of(&k, 0.0)), fabs(t_apo - time_of(&k, pi)));
  MT_CHECK(fabs(t_peri - time_of(&k, 0.0)) < 1e-8);
  MT_CHECK(fabs(t_apo - time_of(&k, pi)) < 1e-8);
  MT_CHECK(fabs(hypot(hypot(y[0], y[1]), y[2]) - k.a * 1.3) < 1e-5);
  /* the ascending node on the plane z = 0: the argument of latitude is 0 */
  myb_plane_event equator = {.body = MYB_EVENT_ORIGIN, .normal = {0.0, 0.0, 2.0}};
  const double t_node = first_event(&k, myb_event_plane, &equator, +1, 2.0 * period, y);
  const double t_down = first_event(&k, myb_event_plane, &equator, -1, 2.0 * period, y);
  printf("ascending node %.2e s, descending %.2e s\n",
         fabs(t_node - time_of(&k, two_pi - k.el.argp)),
         fabs(t_down - time_of(&k, pi - k.el.argp)));
  MT_CHECK(fabs(t_node - time_of(&k, two_pi - k.el.argp)) < 1e-8);
  MT_CHECK(fabs(t_down - time_of(&k, pi - k.el.argp)) < 1e-8);
  /* a distance: r = p / (1 + e cos nu) */
  myb_range_event sphere = {.body = MYB_EVENT_ORIGIN, .range = 1.1e7};
  const double nu_out = acos((k.el.p / sphere.range - 1.0) / k.el.e);
  const double t_out = first_event(&k, myb_event_range, &sphere, +1, 2.0 * period, y);
  const double t_in = first_event(&k, myb_event_range, &sphere, -1, 2.0 * period, y);
  printf("leaving the sphere %.2e s, entering %.2e s\n", fabs(t_out - time_of(&k, nu_out)),
         fabs(t_in - time_of(&k, two_pi - nu_out)));
  MT_CHECK(fabs(t_out - time_of(&k, nu_out)) < 1e-8);
  MT_CHECK(fabs(t_in - time_of(&k, two_pi - nu_out)) < 1e-8);
  myb_force_model_destroy(k.model);
}

/* Events that need the bodies of an environment. */
static void test_world_events(void) {
  world w;
  world_init(&w, 0, 0);
  /* an inclined eccentric orbit that dips below 500 km: for the altitude and the shape */
  const double a = 1.0e7;
  const myb_elements el = {.p = a * (1.0 - 0.32 * 0.32), .e = 0.32, .i = 0.9, .raan = 0.4,
                           .argp = 1.1, .nu = 2.5};
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(mu, &el, &pv) == MYB_OK);
  double y0[6];
  memcpy(y0, pv.r, 3 * sizeof(double));
  memcpy(y0 + 3, pv.v, 3 * sizeof(double));
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(w.model, NULL, NULL, &orbit) == MYB_OK);
  const double day = 86400.0;
  int found = 0;
  double t = 0.0;
  double y[6];

  /* the closest approach to the Moon: the distance has a minimum there */
  myb_range_event moon = {.env = w.env, .body = w.moon};
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, 10.0 * day) == MYB_OK);
  MT_CHECK(myb_orbit_event(orbit, myb_event_range_rate, &moon, 10.0 * day, +1, &found, &t, y) ==
           MYB_OK);
  MT_CHECK(found);
  if (found) {
    double at[3];
    double around[2];
    MT_CHECK(myb_orbit_start(orbit, 0.0, y0, 10.0 * day) == MYB_OK);
    double ya[6];
    MT_CHECK(myb_orbit_state(orbit, t - 5.0, ya) == MYB_OK);
    around[0] = myb_event_range(&moon, t - 5.0, ya);
    MT_CHECK(myb_orbit_state(orbit, t, ya) == MYB_OK);
    at[0] = myb_event_range(&moon, t, ya);
    MT_CHECK(fabs(myb_event_range_rate(&moon, t, ya)) < 1e-6);
    MT_CHECK(myb_orbit_state(orbit, t + 5.0, ya) == MYB_OK);
    around[1] = myb_event_range(&moon, t + 5.0, ya);
    printf("closest approach to the Moon at %.3f s: %.3f km, %.2e m more 5 s before and after\n",
           t, at[0] / 1e3, 0.5 * (around[0] + around[1]) - at[0]);
    MT_CHECK(around[0] > at[0] && around[1] > at[0]);
  }

  /* an altitude over the ellipsoid: the geodetic height at the event */
  myb_altitude_event level = {.env = w.env, .body = w.earth,
                              .ellipsoid = {.a = 6378137.0, .f = 1.0 / 298.257223563},
                              .altitude = 1.0e6};
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, day) == MYB_OK);
  MT_CHECK(myb_orbit_event(orbit, myb_event_altitude, &level, day, -1, &found, &t, y) == MYB_OK);
  MT_CHECK(found);
  if (found) {
    myb_rotation r;
    spin_matrix(&w.spin, t, &r);
    double fixed[3];
    myb_rotation_apply_position(&r, y, fixed);
    double lon = 0.0;
    double lat = 0.0;
    double height = 0.0;
    MT_CHECK(myb_geodetic_from_cartesian(&level.ellipsoid, fixed, &lon, &lat, &height) == MYB_OK);
    printf("descent through 1000 km at %.3f s: geodetic height off by %.2e m\n", t,
           height - level.altitude);
    MT_CHECK(fabs(height - level.altitude) < 1e-5);
    MT_CHECK(myb_event_altitude(&level, t - 1.0, y) != myb_event_altitude(&level, t, y));
  }

  /* rise over a station: the elevation of the measurement at the event */
  myb_station station;
  MT_CHECK(myb_station_init(&station, w.env, w.earth, &level.ellipsoid, 0.5, 0.6, 250.0) ==
           MYB_OK);
  myb_elevation_event horizon = {.station = &station, .elevation = 10.0 * pi / 180.0};
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, 2.0 * day) == MYB_OK);
  int passes = 0;
  double worst_elevation = 0.0;
  /* bound: the passes of two days */
  for (int k = 0; k < 40; ++k) {
    MT_CHECK(myb_orbit_event(orbit, myb_event_elevation, &horizon, 2.0 * day, 0, &found, &t, y) ==
             MYB_OK);
    if (!found) {
      break;
    }
    double azel[2];
    MT_CHECK(myb_measure_azel(&station, t, y, azel, NULL) == MYB_OK);
    worst_elevation = fmax(worst_elevation, fabs(azel[1] - horizon.elevation));
    ++passes;
  }
  printf("station: %d crossings of 10 degrees in two days, elevation there within %.1e rad\n",
         passes, worst_elevation);
  MT_CHECK(passes >= 4 && passes % 2 == 0);
  MT_CHECK(worst_elevation < 1e-11);

  /* the impact on a shape: an octahedron with its vertices 12 000 km out on the axes, its
   * faces 6928 km from the centre: the orbit starts outside it and its pericentre is inside */
  const double s = 1.2e7;
  const double vertices[18] = {s, 0, 0, -s, 0, 0, 0, s, 0, 0, -s, 0, 0, 0, s, 0, 0, -s};
  const size_t faces[24] = {0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4,
                            2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5};
  myb_polyhedron* shape = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, 6, faces, 8, 1.0, NULL, &shape) == MYB_OK);
  myb_shape_event surface = {.env = w.env, .body = w.earth, .shape = shape};
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, day) == MYB_OK);
  MT_CHECK(myb_event_shape(&surface, 0.0, y0) > 0.0);
  MT_CHECK(myb_orbit_event(orbit, myb_event_shape, &surface, day, -1, &found, &t, y) == MYB_OK);
  MT_CHECK(found);
  if (found) {
    myb_rotation r;
    spin_matrix(&w.spin, t, &r);
    double fixed[3];
    myb_rotation_apply_position(&r, y, fixed);
    /* the faces of this octahedron: |x| + |y| + |z| = s */
    const double off = (fabs(fixed[0]) + fabs(fixed[1]) + fabs(fixed[2]) - s) / sqrt(3.0);
    printf("impact on the octahedron at %.6f s, %.2e m from its face\n", t, off);
    MT_CHECK(fabs(off) < 1e-4);
  }
  myb_polyhedron_destroy(shape);

  /* a body that fails: no event, no error */
  myb_range_event origin_only = {.env = w.env, .body = w.earth, .range = 1.0e7};
  MT_CHECK(fabs(myb_event_range(&origin_only, 0.0, y0) -
                (hypot(hypot(y0[0], y0[1]), y0[2]) - 1.0e7)) < 1e-6);
  myb_orbit_destroy(orbit);
  world_destroy(&w);
}

/* ---- Eclipses and surface geometry against CSPICE ------------------------------------------- */

enum { MAX_WINDOWS = 64, MAX_GEOMETRY = 32, DUMP_STATES = 4321 };

typedef struct window {
  double start;
  double stop;
} window;

typedef struct geometry_case {
  double et;
  double v[15]; /* intercept (3), nearest (3), altitude, centre intercept (3), phase,
                   incidence, emission, at the nearest point: incidence, emission */
} geometry_case;

static const double triaxial[3] = {6378136.6, 6350000.0, 6300000.0};
static const double boresight_tilt[3] = {0.3, -0.5, 0.8};

static int collect(world* w, int inner, int direction, double span, double* times, int max) {
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(w->model, NULL, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, w->y0, span) == MYB_OK);
  myb_eclipse_event edge = {.env = w->env, .source = w->sun, .occulter = w->earth,
                            .inner = inner};
  int count = 0;
  /* bound: max events */
  while (count < max) {
    int found = 0;
    double t = 0.0;
    MT_CHECK(myb_orbit_event(orbit, myb_event_eclipse, &edge, span, direction, &found, &t,
                             NULL) == MYB_OK);
    if (!found) {
      break;
    }
    times[count++] = t;
  }
  myb_orbit_destroy(orbit);
  return count;
}

static myb_status rest_state(void* context, myb_time t, double y[6]) {
  (void)context;
  (void)t;
  memset(y, 0, 6 * sizeof(double));
  return MYB_OK;
}

static void dump_kernel(world* w, double span, const char* path) {
  static double epochs[DUMP_STATES];
  static double states[6 * DUMP_STATES];
  for (int k = 0; k < DUMP_STATES; ++k) {
    epochs[k] = span * k / (DUMP_STATES - 1);
  }
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(w->model, NULL, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, w->y0, span) == MYB_OK);
  MT_CHECK(myb_orbit_states(orbit, epochs, DUMP_STATES, states) == MYB_OK);
  myb_orbit_destroy(orbit);
  myb_spk_writer* writer = NULL;
  MT_CHECK(myb_spk_writer_create("test_events", NULL, &writer) == MYB_OK);
  const myb_spk_segment_header object = {.target = -999, .center = 399, .frame = 1,
                                         .name = "object"};
  MT_CHECK(myb_spk_writer_add_states(writer, &object, 13, 11, epochs, states, DUMP_STATES) ==
           MYB_OK);
  const myb_spk_segment_header sun = {.target = 10, .center = 399, .frame = 1, .name = "sun"};
  double fit = 0.0;
  MT_CHECK(myb_spk_writer_add_chebyshev(writer, &sun, 3, mover_state, &w->sun_motion, -100.0,
                                        span + 100.0, 48, 12, &fit) == MYB_OK);
  printf("dump: the Sun's polynomials within %.1e m\n", fit);
  /* CSPICE wants every body tied to the barycentre of the solar system: the Earth rests there */
  const myb_spk_segment_header earth = {.target = 399, .center = 0, .frame = 1, .name = "earth"};
  MT_CHECK(myb_spk_writer_add_chebyshev(writer, &earth, 3, rest_state, NULL, -100.0,
                                        span + 100.0, 1, 2, NULL) == MYB_OK);
  MT_CHECK(myb_spk_writer_write_file(writer, path) == MYB_OK);
  myb_spk_writer_destroy(writer);
}

static void test_cspice(void) {
  world w;
  world_init(&w, 0, 0);
  const double span = 86400.0;
  const char* dump = getenv("MYB_EVENTS_DUMP");
  if (dump != NULL) {
    dump_kernel(&w, span, dump);
  }
  FILE* f = fopen(MYB_TEST_DATA_DIR "/events/cspice.ref", "r");
  MT_CHECK(f != NULL || dump != NULL);
  if (f == NULL) {
    world_destroy(&w);
    return;
  }
  static window any[MAX_WINDOWS];
  static window full[MAX_WINDOWS];
  static geometry_case cases[MAX_GEOMETRY];
  int any_count = 0;
  int full_count = 0;
  int case_count = 0;
  char line[1024];
  /* bound: the lines of the file */
  while (fgets(line, sizeof line, f) != NULL) {
    window win;
    if (sscanf(line, "ANY %lf %lf", &win.start, &win.stop) == 2 && any_count < MAX_WINDOWS) {
      any[any_count++] = win;
    } else if (sscanf(line, "FULL %lf %lf", &win.start, &win.stop) == 2 &&
               full_count < MAX_WINDOWS) {
      full[full_count++] = win;
    } else if (strncmp(line, "GEOM ", 5) == 0 && case_count < MAX_GEOMETRY) {
      geometry_case* c = &cases[case_count];
      int used = 0;
      const char* p = line + 5;
      int ok = sscanf(p, "%lf%n", &c->et, &used) == 1;
      for (int k = 0; k < 15 && ok; ++k) {
        p += used;
        ok = sscanf(p, "%lf%n", &c->v[k], &used) == 1;
      }
      MT_CHECK(ok);
      case_count += ok;
    }
  }
  fclose(f);
  MT_CHECK(any_count >= 5 && full_count >= 5 && case_count >= 10);

  /* eclipses: the four contacts of every passage through the shadow */
  static double t_in[2][MAX_WINDOWS];
  static double t_out[2][MAX_WINDOWS];
  double worst = 0.0;
  for (int inner = 0; inner < 2; ++inner) {
    const window* ref = inner ? full : any;
    const int expected = inner ? full_count : any_count;
    const int entries = collect(&w, inner, -1, span, t_in[inner], MAX_WINDOWS);
    const int exits = collect(&w, inner, +1, span, t_out[inner], MAX_WINDOWS);
    /* a window that the span cuts has one end only */
    static double ref_in[MAX_WINDOWS];
    static double ref_out[MAX_WINDOWS];
    int n_in = 0;
    int n_out = 0;
    for (int k = 0; k < expected; ++k) {
      if (ref[k].start > 0.0) {
        ref_in[n_in++] = ref[k].start;
      }
      if (ref[k].stop < span) {
        ref_out[n_out++] = ref[k].stop;
      }
    }
    MT_CHECK(entries == n_in && exits == n_out);
    int compared = 0;
    for (int k = 0; k < entries && k < n_in; ++k) {
      worst = fmax(worst, fabs(t_in[inner][k] - ref_in[k]));
      ++compared;
    }
    for (int k = 0; k < exits && k < n_out; ++k) {
      worst = fmax(worst, fabs(t_out[inner][k] - ref_out[k]));
      ++compared;
    }
    printf("%s: %d entries and %d exits, %d contacts compared with gfoclt\n",
           inner ? "umbra" : "penumbra", entries, exits, compared);
    MT_CHECK(compared >= 2 * expected - 2);
  }
  printf("eclipse contacts against CSPICE gfoclt: within %.2e s\n", worst);
  MT_CHECK(worst < 2.0e-7);

  /* the surface: intercept of a ray, the point under the object, the angles of illumination */
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(w.model, NULL, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, w.y0, span) == MYB_OK);
  double worst_point = 0.0;
  double worst_angle = 0.0;
  double worst_altitude = 0.0;
  for (int k = 0; k < case_count; ++k) {
    const geometry_case* c = &cases[k];
    double y[6];
    MT_CHECK(myb_orbit_state(orbit, c->et, y) == MYB_OK);
    double sun_state[6];
    MT_CHECK(myb_env_state(w.env, c->et, w.sun, sun_state) == MYB_OK);
    myb_orientation o;
    MT_CHECK(myb_env_orientation(w.env, c->et, w.earth, &o) == MYB_OK);
    double observer[3];
    double sun[3];
    myb_rotation_apply_position(&o.rotation, y, observer);
    myb_rotation_apply_position(&o.rotation, sun_state, sun);
    /* the boresight: off the nadir by the fixed tilt */
    const double r = hypot(hypot(y[0], y[1]), y[2]);
    const double tilt = hypot(hypot(boresight_tilt[0], boresight_tilt[1]), boresight_tilt[2]);
    double look_inertial[3];
    double look[3];
    for (int i = 0; i < 3; ++i) {
      look_inertial[i] = -y[i] / r + 0.2 * boresight_tilt[i] / tilt;
    }
    myb_rotation_apply_position(&o.rotation, look_inertial, look);
    double hit[3];
    double distance = 0.0;
    MT_CHECK(myb_ray_ellipsoid(triaxial, observer, look, &distance, hit) == 1);
    double near[3];
    double altitude = 0.0;
    MT_CHECK(myb_ellipsoid_nearest(triaxial, observer, near, &altitude) == MYB_OK);
    const double down[3] = {-observer[0], -observer[1], -observer[2]};
    double under[3];
    MT_CHECK(myb_ray_ellipsoid(triaxial, observer, down, NULL, under) == 1);
    double normal[3];
    double angles[5];
    myb_ellipsoid_normal(triaxial, hit, normal);
    MT_CHECK(myb_illumination_angles(hit, normal, sun, observer, &angles[0], &angles[1],
                                     &angles[2]) == MYB_OK);
    myb_ellipsoid_normal(triaxial, near, normal);
    MT_CHECK(myb_illumination_angles(near, normal, sun, observer, NULL, &angles[3],
                                     &angles[4]) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      worst_point = fmax(worst_point, fabs(hit[i] - c->v[i]));
      worst_point = fmax(worst_point, fabs(near[i] - c->v[3 + i]));
      worst_point = fmax(worst_point, fabs(under[i] - c->v[7 + i]));
    }
    worst_altitude = fmax(worst_altitude, fabs(altitude - c->v[6]));
    for (int i = 0; i < 5; ++i) {
      worst_angle = fmax(worst_angle, fabs(angles[i] - c->v[10 + i]));
    }
    /* the point under the object sees it at the zenith */
    MT_CHECK(angles[4] < 1e-9);
  }
  myb_orbit_destroy(orbit);
  printf("surface against CSPICE sincpt, subpnt, ilumin (%d epochs): points %.2e m, altitude "
         "%.2e m, angles %.2e rad\n",
         case_count, worst_point, worst_altitude, worst_angle);
  MT_CHECK(worst_point < 1e-6 && worst_altitude < 1e-6 && worst_angle < 1e-12);
  world_destroy(&w);
}

/* Reads up to `max` numbers from a line. */
static int numbers(const char* line, double* out, int max) {
  int count = 0;
  int used = 0;
  while (count < max && sscanf(line, "%lf%n", &out[count], &used) == 1) {
    line += used;
    ++count;
  }
  return count;
}

/* With light time, the limb and the terminator: sincpt_c with CN, spkpos_c with CN and CN+S,
 * edlimb_c, edterm_c. */
static void test_cspice_views(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/events/cspice.ref", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  world w;
  world_init(&w, 0, 0);
  const double span = 86400.0;
  const double c = 299792458.0;
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(w.model, NULL, NULL, &orbit) == MYB_OK);
  char line[2048];
  double worst_point = 0.0;
  double worst_time = 0.0;
  double worst_sun = 0.0;
  double worst_sun_time = 0.0;
  double worst_aberration = 0.0;
  double worst_limb = 0.0;
  double worst_terminator = 0.0;
  double worst_own = 0.0;
  double least_moved = (double)INFINITY;
  int counts[4] = {0, 0, 0, 0};
  /* bound: the lines of the file */
  while (fgets(line, sizeof line, f) != NULL) {
    double v[16];
    const int limb = strncmp(line, "LIMB ", 5) == 0;
    const int sincn = strncmp(line, "SINCN ", 6) == 0;
    const int suncn = strncmp(line, "SUNCN ", 6) == 0;
    const int term = strncmp(line, "TERM ", 5) == 0;
    if (!limb && !sincn && !suncn && !term) {
      continue;
    }
    int penumbral = 0;
    const char* rest = line + (limb || term ? 5 : 6);
    if (term) {
      penumbral = strncmp(rest, "PENUMBRAL", 9) == 0;
      rest += penumbral ? 10 : 7;
    }
    const int n = numbers(rest, v, 16);
    const double et = v[0];
    /* the object and the bodies at et, in the body-fixed frame of the Earth */
    double y[6];
    MT_CHECK(myb_orbit_start(orbit, 0.0, w.y0, span) == MYB_OK);
    MT_CHECK(myb_orbit_state(orbit, et, y) == MYB_OK);
    myb_orientation o;
    MT_CHECK(myb_env_orientation(w.env, et, w.earth, &o) == MYB_OK);
    double sun_state[6];
    MT_CHECK(myb_env_state(w.env, et, w.sun, sun_state) == MYB_OK);
    double observer[3];
    double sun[3];
    myb_rotation_apply_position(&o.rotation, y, observer);
    myb_rotation_apply_position(&o.rotation, sun_state, sun);
    if (sincn && n == 5) {
      const double r = hypot(hypot(y[0], y[1]), y[2]);
      const double tilt = hypot(hypot(boresight_tilt[0], boresight_tilt[1]), boresight_tilt[2]);
      double look[3];
      for (int i = 0; i < 3; ++i) {
        look[i] = -y[i] / r + 0.2 * boresight_tilt[i] / tilt;
      }
      int hit = 0;
      double point[3];
      double lt = 0.0;
      MT_CHECK(myb_view_intercept(w.env, et, w.earth, triaxial, y, look, c, &hit, point, &lt) ==
               MYB_OK);
      MT_CHECK(hit);
      for (int i = 0; i < 3; ++i) {
        worst_point = fmax(worst_point, fabs(point[i] - v[1 + i]));
      }
      worst_time = fmax(worst_time, fabs(lt - v[4]));
      /* the Earth turned on while the light travelled: the point of the body as it is now
       * is metres away */
      double look_now[3];
      double still[3];
      myb_rotation_apply_position(&o.rotation, look, look_now);
      MT_CHECK(myb_ray_ellipsoid(triaxial, observer, look_now, NULL, still) == 1);
      const double moved =
          hypot(hypot(still[0] - point[0], still[1] - point[1]), still[2] - point[2]);
      least_moved = fmin(least_moved, moved);
      ++counts[0];
    } else if (suncn && n == 8) {
      double position[3];
      double lt = 0.0;
      MT_CHECK(myb_view_body(w.env, et, w.sun, y, c, 0, position, NULL, &lt) == MYB_OK);
      double seen[3];
      MT_CHECK(myb_view_body(w.env, et, w.sun, y, c, 1, seen, NULL, NULL) == MYB_OK);
      for (int i = 0; i < 3; ++i) {
        worst_sun = fmax(worst_sun, fabs(position[i] - v[1 + i]));
        worst_aberration = fmax(worst_aberration, fabs(seen[i] - v[5 + i]));
      }
      worst_sun_time = fmax(worst_sun_time, fabs(lt - v[4]));
      /* the aberration itself: a million metres at the distance of the Sun */
      MT_CHECK(hypot(hypot(seen[0] - position[0], seen[1] - position[1]),
                     seen[2] - position[2]) > 1.0e5);
      ++counts[1];
    } else if (limb && n == 10) {
      double centre[3];
      double major[3];
      double minor[3];
      MT_CHECK(myb_ellipsoid_limb(triaxial, observer, centre, major, minor) == MYB_OK);
      /* the axes of an ellipse are given up to their signs */
      const double sign_major = major[0] * v[4] + major[1] * v[5] + major[2] * v[6] < 0.0 ? -1.0 : 1.0;
      const double sign_minor = minor[0] * v[7] + minor[1] * v[8] + minor[2] * v[9] < 0.0 ? -1.0 : 1.0;
      for (int i = 0; i < 3; ++i) {
        worst_limb = fmax(worst_limb, fabs(centre[i] - v[1 + i]));
        worst_limb = fmax(worst_limb, fabs(sign_major * major[i] - v[4 + i]));
        worst_limb = fmax(worst_limb, fabs(sign_minor * minor[i] - v[7 + i]));
      }
      /* a point of it: on the surface, the line of sight in the tangent plane */
      const double p[3] = {centre[0] + 0.6 * major[0] + 0.8 * minor[0],
                           centre[1] + 0.6 * major[1] + 0.8 * minor[1],
                           centre[2] + 0.6 * major[2] + 0.8 * minor[2]};
      double normal[3];
      myb_ellipsoid_normal(triaxial, p, normal);
      const double sight[3] = {p[0] - observer[0], p[1] - observer[1], p[2] - observer[2]};
      MT_CHECK(fabs(normal[0] * sight[0] + normal[1] * sight[1] + normal[2] * sight[2]) <
               1e-8 * hypot(hypot(sight[0], sight[1]), sight[2]));
      ++counts[2];
    } else if (term && n == 13) {
      /* the centre of the Sun over the tangent plane: -R at the edge of the dark, +R at the
       * edge of the full light */
      const double level = penumbral ? sun_radius : -sun_radius;
      for (int k = 0; k < 4; ++k) {
        const double* p = v + 1 + 3 * k;
        double normal[3];
        myb_ellipsoid_normal(triaxial, p, normal);
        const double over = normal[0] * (sun[0] - p[0]) + normal[1] * (sun[1] - p[1]) +
                            normal[2] * (sun[2] - p[2]);
        worst_terminator = fmax(worst_terminator, fabs(over - level));
      }
      double points[12];
      MT_CHECK(myb_ellipsoid_terminator(triaxial, sun, sun_radius, penumbral, 4, points) ==
               MYB_OK);
      for (int k = 0; k < 4; ++k) {
        const double* p = points + 3 * k;
        double normal[3];
        myb_ellipsoid_normal(triaxial, p, normal);
        const double over = normal[0] * (sun[0] - p[0]) + normal[1] * (sun[1] - p[1]) +
                            normal[2] * (sun[2] - p[2]);
        worst_own = fmax(worst_own, fabs(over - level));
        double on = 0.0;
        for (int i = 0; i < 3; ++i) {
          on += (p[i] / triaxial[i]) * (p[i] / triaxial[i]);
        }
        MT_CHECK(fabs(on - 1.0) < 1e-14);
      }
      ++counts[3];
    }
  }
  fclose(f);
  myb_orbit_destroy(orbit);
  printf("with light time against CSPICE: the ray (CN, %d epochs) %.2e m, %.1e s; the Sun "
         "from the object (%d) %.2e m, %.1e s, with aberration %.1f m (CSPICE's is the "
         "classical one: (v/c)^2 of 1 au is 25 m)\n",
         counts[0], worst_point, worst_time, counts[1], worst_sun, worst_sun_time,
         worst_aberration);
  printf("limb against edlimb (%d epochs): %.2e m; terminator (%d): the points of edterm "
         "within %.1e m of the horizon condition, those of the library %.1e m\n",
         counts[2], worst_limb, counts[3], worst_terminator, worst_own);
  MT_CHECK(counts[0] == 20 && counts[1] == 20 && counts[2] == 20 && counts[3] == 12);
  /* CSPICE gives the epoch of the point, not the light time: the difference of two epochs
   * of 8e4 s is known to 1.5e-11 s */
  MT_CHECK(worst_point < 1e-5 && worst_time < 5e-11);
  printf("the point seen and the point of the body as it is now: at least %.2f m apart\n",
         least_moved);
  MT_CHECK(least_moved > 1.0);
  MT_CHECK(worst_sun < 1e-2 && worst_sun_time < 1e-10);
  MT_CHECK(worst_aberration < 100.0);
  MT_CHECK(worst_limb < 1e-6);
  MT_CHECK(worst_terminator < 1.0 && worst_own < 1e-3);
  /* errors */
  double centre[3];
  double major[3];
  double minor[3];
  const double inside[3] = {1.0e6, 0.0, 0.0};
  MT_CHECK(myb_ellipsoid_limb(triaxial, inside, centre, major, minor) == MYB_E_DOMAIN);
  double points[12];
  MT_CHECK(myb_ellipsoid_terminator(triaxial, inside, 1.0, 0, 4, points) == MYB_E_DOMAIN);
  const double far[3] = {1.0e11, 2.0e10, -3.0e10};
  MT_CHECK(myb_ellipsoid_terminator(triaxial, far, -1.0, 0, 4, points) == MYB_E_DOMAIN);
  /* a point source: both edges are where it is on the horizon */
  double dark[12];
  double lit[12];
  MT_CHECK(myb_ellipsoid_terminator(triaxial, far, 0.0, 0, 4, dark) == MYB_OK);
  MT_CHECK(myb_ellipsoid_terminator(triaxial, far, 0.0, 1, 4, lit) == MYB_OK);
  MT_CHECK(memcmp(dark, lit, sizeof dark) == 0);
  /* a line of sight that misses */
  int hit = 1;
  const double away[3] = {w.y0[0], w.y0[1], w.y0[2]};
  MT_CHECK(myb_view_intercept(w.env, 100.0, w.earth, triaxial, w.y0, away, c, &hit, NULL,
                              NULL) == MYB_OK);
  MT_CHECK(hit == 0);
  MT_CHECK(myb_view_body(w.env, 100.0, w.sun, w.y0, -1.0, 0, NULL, NULL, NULL) == MYB_E_DOMAIN);
  world_destroy(&w);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_breaks();
  test_breaks_dop853();
  test_partials();
  test_kepler_events();
  test_world_events();
  test_cspice();
  test_cspice_views();
  return mt_finish("test_events");
}
