/* SPDX-License-Identifier: BSD-3-Clause */
/* The environment of a force model and its ready-made terms (myb/env.h), with providers made
 * for the test (bodies on circles, uniformly rotating frames), so nothing is downloaded:
 *  - scales: the tabulated difference of two scales against the exact conversion (TDB - TT
 *    from the series of Fairhead & Bretagnon, TCB - TT with its rate), and the same bits
 *    whatever was asked before;
 *  - bodies: what the providers are asked and how often (once per instant), the origin, errors;
 *  - every term against the same functions of myb/grav.h, myb/force.h and myb/rel.h called by
 *    hand with the providers' values: bit for bit, the terms being nothing but that
 *    composition; analytic Jacobians against central differences; the field of a body away
 *    from the origin against the point-mass perturbation; drag against v - omega x r;
 *  - the layered atmosphere against its formula. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/myb.h"
#include "num/vec3.h"

static const double c_light = 299792458.0;

/* ---- Providers of the test ----------------------------------------------------------------- */

typedef struct mover {
  double radius;
  double rate;
  double phase;
  double inclination;
  int calls;
  myb_time last;
  myb_status fail;
  int poison;
} mover;

static myb_status mover_state(void* context, myb_time t, double y[6]) {
  mover* m = context;
  m->calls++;
  m->last = t;
  if (m->fail != MYB_OK) {
    return m->fail;
  }
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
  if (m->poison) {
    y[4] = (double)NAN;
  }
  return MYB_OK;
}

typedef struct spinner {
  double angles[3];
  double rate;
  double ut_shift;
  double xp;
  double yp;
  int calls;
  myb_time last;
} spinner;

static myb_status spinner_orientation(void* context, myb_time t, myb_orientation* o) {
  spinner* s = context;
  s->calls++;
  s->last = t;
  const double angles[3] = {s->angles[0], s->angles[1],
                            s->angles[2] + s->rate * myb_time_seconds(t)};
  const double rates[3] = {0.0, 0.0, s->rate};
  myb_rotation_euler313(angles, rates, &o->rotation);
  o->ut = myb_time_add(t, s->ut_shift);
  o->xp = s->xp;
  o->yp = s->yp;
  return MYB_OK;
}

static double rel3(const double a[3], const double b[3]) {
  const double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  return myb_vec3_norm(d) / myb_vec3_norm(b);
}

static void check_same3(const double a[3], const double b[3]) {
  for (int i = 0; i < 3; ++i) {
    MT_CHECK_SAME(a[i], b[i]);
  }
}

/* ---- Scales ---------------------------------------------------------------------------------- */

static void test_scales(void) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/time-scales.tk", &diag) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/tdb-tt-fb1990.tk", &diag) == MYB_OK);
  myb_timesys* sys = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_scale tt;
  myb_scale tdb;
  myb_scale tcb;
  MT_CHECK(myb_timesys_scale(sys, "TT", &tt) == MYB_OK);
  MT_CHECK(myb_timesys_scale(sys, "TDB", &tdb) == MYB_OK);
  MT_CHECK(myb_timesys_scale(sys, "TCB", &tcb) == MYB_OK);
  const myb_calendar date = {2025, 6, 1, 0, 0, 0.0};
  myb_time epoch;
  MT_CHECK(myb_calendar_to_time(&date, &epoch) == MYB_OK);

  mover sun = {.radius = 1.5e11, .rate = 2e-7};
  mover moon = {.radius = 3.8e8, .rate = 2.7e-6};
  const myb_env_options options = {.epoch = epoch, .timesys = sys, .scale = "TT"};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  size_t i_sun = 99;
  size_t i_moon = 99;
  const myb_env_body b_sun = {.id = 10, .gm = 1.3e20, .state = mover_state,
                              .state_context = &sun, .state_scale = "TDB"};
  const myb_env_body b_moon = {.id = 301, .gm = 4.9e12, .state = mover_state,
                               .state_context = &moon, .state_scale = "TCB"};
  MT_CHECK(myb_env_add_body(env, &b_sun, &i_sun) == MYB_OK && i_sun == 0);
  MT_CHECK(myb_env_add_body(env, &b_moon, &i_moon) == MYB_OK && i_moon == 1);

  /* the table against the exact conversion, 40 days either way */
  mt_rng rng = {.state = 12};
  double worst_tdb = 0.0;
  double worst_tcb = 0.0;
  for (int k = 0; k < 300; ++k) {
    const double t = 40.0 * 86400.0 * ((double)(mt_next(&rng) >> 11) / 4503599627370496.0 - 1.0);
    const myb_time model = myb_time_add(epoch, t);
    myb_time exact;
    myb_time got;
    MT_CHECK(myb_timesys_convert(sys, model, tt, tdb, &exact) == MYB_OK);
    MT_CHECK(myb_env_time(env, t, "TDB", &got) == MYB_OK);
    worst_tdb = fmax(worst_tdb, fabs(myb_dd_to_double(myb_time_diff(got, exact))));
    MT_CHECK(myb_timesys_convert(sys, model, tt, tcb, &exact) == MYB_OK);
    MT_CHECK(myb_env_time(env, t, "TCB", &got) == MYB_OK);
    worst_tcb = fmax(worst_tcb, fabs(myb_dd_to_double(myb_time_diff(got, exact))));
    MT_CHECK(myb_env_time(env, t, NULL, &got) == MYB_OK);
    MT_CHECK(myb_time_compare(got, model) == 0);
    MT_CHECK(myb_env_time(env, t, "TT", &got) == MYB_OK);
    MT_CHECK(myb_time_compare(got, model) == 0);
  }
  printf("scales: table of one hour against the exact conversion: TDB %.2e s, TCB %.2e s\n",
         worst_tdb, worst_tcb);
  MT_CHECK(worst_tdb < 1e-13);
  MT_CHECK(worst_tcb < 1e-13);

  /* a node of the table is the exact value */
  myb_time exact;
  myb_time got;
  MT_CHECK(myb_timesys_convert(sys, myb_time_add(epoch, -7200.0), tt, tdb, &exact) == MYB_OK);
  MT_CHECK(myb_env_time(env, -7200.0, "TDB", &got) == MYB_OK);
  MT_CHECK(fabs(myb_dd_to_double(myb_time_diff(got, exact))) < 2.5e-16);

  /* the same bits whatever was asked before: a fresh environment against a used one */
  myb_env* fresh = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &fresh) == MYB_OK);
  MT_CHECK(myb_env_add_body(fresh, &b_sun, NULL) == MYB_OK);
  const double probes[4] = {123456.789, -3599.5, 3600.0, 2.5e6};
  for (int k = 0; k < 4; ++k) {
    myb_time a;
    myb_time b;
    MT_CHECK(myb_env_time(env, probes[k], "TDB", &a) == MYB_OK);
    MT_CHECK(myb_env_time(fresh, probes[k], "TDB", &b) == MYB_OK);
    MT_CHECK(a.s == b.s);
    MT_CHECK_SAME(a.f, b.f);
  }
  myb_env_destroy(fresh);

  /* the providers get the instant on their own scale */
  double state[6];
  const double t = 5000.25;
  MT_CHECK(myb_env_state(env, t, i_sun, state) == MYB_OK);
  MT_CHECK(myb_env_time(env, t, "TDB", &got) == MYB_OK);
  MT_CHECK(myb_time_compare(sun.last, got) == 0);
  MT_CHECK(myb_env_state(env, t, i_moon, state) == MYB_OK);
  MT_CHECK(myb_env_time(env, t, "TCB", &got) == MYB_OK);
  MT_CHECK(myb_time_compare(moon.last, got) == 0);

  /* errors */
  MT_CHECK(myb_env_time(env, t, "GPS", &got) == MYB_E_NOT_FOUND); /* no provider uses it */
  MT_CHECK(myb_env_time(env, t, "XYZ", &got) == MYB_E_NOT_FOUND);
  MT_CHECK(myb_env_time(env, (double)NAN, "TDB", &got) == MYB_E_DOMAIN);
  MT_CHECK(myb_env_time(env, 1e18, "TDB", &got) == MYB_E_RANGE);
  const myb_env_body bad_scale = {.gm = 1.0, .state = mover_state, .state_context = &sun,
                                  .state_scale = "XYZ"};
  MT_CHECK(myb_env_add_body(env, &bad_scale, NULL) == MYB_E_NOT_FOUND);
  MT_CHECK(myb_env_body_count(env) == 2);
  myb_env_destroy(env);

  myb_env* none = NULL;
  const myb_env_options no_scale = {.epoch = epoch, .timesys = sys};
  MT_CHECK(myb_env_create(&no_scale, NULL, &none) == MYB_E_DOMAIN);
  const myb_env_options wrong_scale = {.epoch = epoch, .timesys = sys, .scale = "XYZ"};
  MT_CHECK(myb_env_create(&wrong_scale, NULL, &none) == MYB_E_NOT_FOUND);
  const myb_env_options bad_step = {.epoch = epoch, .table_step = -1.0};
  MT_CHECK(myb_env_create(&bad_step, NULL, &none) == MYB_E_DOMAIN);
  MT_CHECK(none == NULL);
  myb_timesys_destroy(sys);
}

/* ---- Bodies ---------------------------------------------------------------------------------- */

static void test_bodies(void) {
  const myb_time epoch = myb_time_from_seconds(8.0e8);
  const myb_env_options options = {.epoch = epoch};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  mover moon = {.radius = 3.8e8, .rate = 2.7e-6, .phase = 0.4, .inclination = 0.3};
  spinner spin = {.angles = {0.2, 0.4, 1.0}, .rate = 7.3e-5, .ut_shift = 0.25, .xp = 1e-6,
                  .yp = -2e-6};
  size_t planet = 99;
  size_t i_moon = 99;
  const myb_env_body b_planet = {.id = 399, .gm = 3.986e14, .orientation = spinner_orientation,
                                 .orientation_context = &spin, .radius = 6.4e6};
  const myb_env_body b_moon = {.id = 301, .gm = 4.9e12, .state = mover_state,
                               .state_context = &moon};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_moon, &i_moon) == MYB_OK);
  MT_CHECK(myb_env_body_count(env) == 2);
  MT_CHECK_SAME(myb_env_gm(env, planet), 3.986e14);
  MT_CHECK(myb_time_compare(myb_env_epoch(env), epoch) == 0);
  size_t found = 99;
  MT_CHECK(myb_env_find(env, 301, &found) == MYB_OK && found == i_moon);
  MT_CHECK(myb_env_find(env, 5, &found) == MYB_E_NOT_FOUND);

  /* once per instant, at epoch + t */
  double state[6];
  double expected[6];
  const double t = 1234.5;
  MT_CHECK(myb_env_state(env, t, i_moon, state) == MYB_OK);
  MT_CHECK(myb_env_state(env, t, i_moon, state) == MYB_OK);
  MT_CHECK(moon.calls == 1);
  MT_CHECK(myb_time_compare(moon.last, myb_time_add(epoch, t)) == 0);
  mover twin = moon;
  MT_CHECK(mover_state(&twin, myb_time_add(epoch, t), expected) == MYB_OK);
  for (int i = 0; i < 6; ++i) {
    MT_CHECK_SAME(state[i], expected[i]);
  }
  MT_CHECK(myb_env_state(env, t + 1.0, i_moon, state) == MYB_OK);
  MT_CHECK(myb_env_state(env, t, i_moon, state) == MYB_OK);
  MT_CHECK(moon.calls == 3);
  MT_CHECK(myb_env_state(env, t, planet, state) == MYB_OK);
  for (int i = 0; i < 6; ++i) {
    MT_CHECK_SAME(state[i], 0.0);
  }

  myb_orientation o;
  myb_orientation oe;
  MT_CHECK(myb_env_orientation(env, t, planet, &o) == MYB_OK);
  MT_CHECK(myb_env_orientation(env, t, planet, &o) == MYB_OK);
  MT_CHECK(spin.calls == 1);
  spinner twin_spin = spin;
  MT_CHECK(spinner_orientation(&twin_spin, myb_time_add(epoch, t), &oe) == MYB_OK);
  MT_CHECK(memcmp(&o, &oe, sizeof o) == 0);
  /* a body without a provider keeps the inertial axes */
  MT_CHECK(myb_env_orientation(env, t, i_moon, &o) == MYB_OK);
  MT_CHECK(o.rotation.m[0] == 1.0 && o.rotation.m[4] == 1.0 && o.rotation.m[8] == 1.0 &&
           o.rotation.m[1] == 0.0 && o.rotation.dm[0] == 0.0 && o.xp == 0.0);
  MT_CHECK(myb_time_compare(o.ut, myb_time_add(epoch, t)) == 0);

  /* errors of the providers come through; non-finite results are refused */
  moon.fail = MYB_E_RANGE;
  MT_CHECK(myb_env_state(env, 77.0, i_moon, state) == MYB_E_RANGE);
  moon.fail = MYB_OK;
  moon.poison = 1;
  MT_CHECK(myb_env_state(env, 78.0, i_moon, state) == MYB_E_DOMAIN);
  moon.poison = 0;
  MT_CHECK(myb_env_state(env, 78.0, i_moon, state) == MYB_OK);
  MT_CHECK(myb_env_state(env, (double)INFINITY, i_moon, state) == MYB_E_DOMAIN);

  /* invalid bodies */
  myb_env_body bad = b_planet;
  MT_CHECK(myb_env_add_body(env, &bad, NULL) == MYB_E_DOMAIN); /* a second origin */
  bad = b_moon;
  bad.gm = -1.0;
  MT_CHECK(myb_env_add_body(env, &bad, NULL) == MYB_E_DOMAIN);
  bad = b_moon;
  bad.radius = (double)NAN;
  MT_CHECK(myb_env_add_body(env, &bad, NULL) == MYB_E_DOMAIN);
  bad = b_moon;
  bad.state_scale = "TDB"; /* no time system */
  MT_CHECK(myb_env_add_body(env, &bad, NULL) == MYB_E_DOMAIN);
  for (int k = 2; k < MYB_ENV_MAX_BODIES; ++k) {
    MT_CHECK(myb_env_add_body(env, &b_moon, NULL) == MYB_OK);
  }
  MT_CHECK(myb_env_add_body(env, &b_moon, NULL) == MYB_E_OVERFLOW);
  myb_env_destroy(env);
  myb_env_destroy(NULL);
}

/* ---- Point masses ------------------------------------------------------------------------------ */

static void test_point_masses(void) {
  const double x[3] = {7.0e6, -1.2e6, 2.5e6};
  const double v[3] = {1.0e3, 7.2e3, -0.5e3};
  double a[3];
  double jx[9];
  double jv[9];
  myb_point_mass_term central = {.gm = 3.986004418e14};
  MT_CHECK(myb_term_point_mass(&central, 0.0, x, v, a) == MYB_OK);
  const double r = myb_vec3_norm(x);
  const double k = -central.gm / (r * r * r);
  const double expected[3] = {k * x[0], k * x[1], k * x[2]};
  check_same3(a, expected);
  double jref[9];
  myb_point_mass_jacobian(central.gm, x, jref);
  MT_CHECK(myb_term_point_mass_jacobian(&central, 0.0, x, v, jx, jv) == MYB_OK);
  for (int i = 0; i < 9; ++i) {
    MT_CHECK_SAME(jx[i], jref[i]);
    MT_CHECK_SAME(jv[i], 0.0);
  }
  const double zero[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_term_point_mass(&central, 0.0, zero, v, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_term_point_mass_jacobian(&central, 0.0, zero, v, jx, jv) == MYB_E_DOMAIN);

  /* a body at the origin: direct minus indirect */
  const myb_time epoch = myb_time_from_seconds(8.0e8);
  const myb_env_options options = {.epoch = epoch};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  mover moon = {.radius = 3.8e8, .rate = 2.7e-6, .phase = 0.4, .inclination = 0.3};
  size_t planet = 0;
  size_t i_moon = 0;
  const myb_env_body b_planet = {.gm = 3.986e14};
  const myb_env_body b_moon = {.gm = 4.9e12, .state = mover_state, .state_context = &moon};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_moon, &i_moon) == MYB_OK);
  myb_third_body_term third = {.env = env, .body = i_moon};
  const double t = 4321.0;
  double s[6];
  mover twin = moon;
  MT_CHECK(mover_state(&twin, myb_time_add(epoch, t), s) == MYB_OK);
  double by_hand[3] = {0.0, 0.0, 0.0};
  myb_force_third_body(b_moon.gm, s, x, by_hand);
  MT_CHECK(myb_term_third_body(&third, t, x, v, a) == MYB_OK);
  check_same3(a, by_hand);
  myb_third_body_term self = {.env = env, .body = planet};
  MT_CHECK(myb_term_third_body(&self, t, x, v, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_term_third_body(&third, t, s, v, a) == MYB_E_DOMAIN); /* at the body's centre */

  /* analytic Jacobian against central differences of the same term */
  myb_force_model* analytic = NULL;
  myb_force_model* numeric = NULL;
  MT_CHECK(myb_force_model_create(NULL, &analytic) == MYB_OK);
  MT_CHECK(myb_force_model_create(NULL, &numeric) == MYB_OK);
  MT_CHECK(myb_force_model_add(analytic, "moon", myb_term_third_body,
                               myb_term_third_body_jacobian, &third, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(numeric, "moon", myb_term_third_body, NULL, &third, NULL) ==
           MYB_OK);
  double nx[9];
  double nv[9];
  MT_CHECK(myb_force_model_jacobian(analytic, t, x, v, jx, jv, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_jacobian(numeric, t, x, v, nx, nv, NULL) == MYB_OK);
  double scale = 0.0;
  double worst = 0.0;
  for (int i = 0; i < 9; ++i) {
    scale = fmax(scale, fabs(jx[i]));
    worst = fmax(worst, fabs(jx[i] - nx[i]));
    MT_CHECK_SAME(jv[i], 0.0);
  }
  printf("third body: analytic Jacobian against central differences %.1e\n", worst / scale);
  MT_CHECK(worst < 1e-8 * scale);
  myb_force_model_destroy(analytic);
  myb_force_model_destroy(numeric);
  myb_env_destroy(env);

  /* an inertial origin: the direct attraction only */
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_moon, &i_moon) == MYB_OK);
  third.env = env;
  third.body = i_moon;
  MT_CHECK(myb_term_third_body(&third, t, x, v, a) == MYB_OK);
  const double d[3] = {s[0] - x[0], s[1] - x[1], s[2] - x[2]};
  const double dd = myb_vec3_norm(d);
  const double q = b_moon.gm / (dd * dd * dd);
  const double direct[3] = {q * d[0], q * d[1], q * d[2]};
  check_same3(a, direct);
  myb_env_destroy(env);
}

/* ---- Gravity field ----------------------------------------------------------------------------- */

enum { DEGREE = 4, COEFFICIENTS = (DEGREE + 1) * (DEGREE + 2) / 2 };

static void make_field(double gm, myb_gravity_field** field) {
  double c[COEFFICIENTS] = {0.0};
  double s[COEFFICIENTS] = {0.0};
  c[0] = 1.0;
  c[3] = -4.84165e-4; /* C20 */
  c[4] = 1.3e-9;      /* C21 */
  s[4] = -2.1e-9;
  c[5] = 2.439e-6; /* C22 */
  s[5] = -1.400e-6;
  c[6] = 9.57e-7;  /* C30 */
  c[7] = 2.03e-6;  /* C31 */
  s[7] = 2.48e-7;
  c[10] = 5.4e-7;  /* C40 */
  c[12] = 3.5e-7;  /* C42 */
  s[12] = 6.6e-7;
  MT_CHECK(myb_gravity_field_create(gm, 6378136.3, DEGREE, c, s, NULL, field) == MYB_OK);
}

static void test_field(void) {
  const myb_time epoch = myb_time_from_seconds(8.0e8);
  const myb_env_options options = {.epoch = epoch};
  const double gm = 3.986004415e14;
  const double x[3] = {7.0e6, -1.2e6, 2.5e6};
  const double v[3] = {1.0e3, 7.2e3, -0.5e3};
  const double t = 9876.5;
  const myb_time when = myb_time_add(epoch, t);
  spinner spin = {.angles = {0.2, 0.4, 1.0}, .rate = 7.3e-5, .ut_shift = -0.3, .xp = 6e-7,
                  .yp = 1.6e-6};
  mover sun = {.radius = 1.5e11, .rate = 2e-7, .phase = 1.0, .inclination = 0.4};
  mover moon = {.radius = 3.8e8, .rate = 2.7e-6, .phase = 0.4, .inclination = 0.3};
  myb_gravity_field* field = NULL;
  myb_gravity_field* reference = NULL;
  myb_gravity_work* work = NULL;
  make_field(gm, &field);
  MT_CHECK(myb_gravity_field_clone(field, NULL, &reference) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(DEGREE, NULL, &work) == MYB_OK);

  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  size_t planet = 0;
  size_t i_sun = 0;
  size_t i_moon = 0;
  const myb_env_body b_planet = {.gm = gm, .orientation = spinner_orientation,
                                 .orientation_context = &spin};
  const myb_env_body b_sun = {.gm = 1.32712440018e20, .state = mover_state,
                              .state_context = &sun};
  const myb_env_body b_moon = {.gm = 4.9028e12, .state = mover_state, .state_context = &moon};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_sun, &i_sun) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_moon, &i_moon) == MYB_OK);

  /* the body at the origin: rotate, evaluate, rotate back */
  myb_field_term term = {.env = env, .body = planet, .field = field, .work = work, .n_min = 0,
                         .n_max = -1};
  double a[3];
  MT_CHECK(myb_term_field(&term, t, x, v, a) == MYB_OK);
  myb_orientation o;
  spinner twin = spin;
  MT_CHECK(spinner_orientation(&twin, when, &o) == MYB_OK);
  myb_rotation inverse;
  myb_rotation_invert(&o.rotation, &inverse);
  double r_fixed[3];
  double g_fixed[3];
  double by_hand[3];
  myb_rotation_apply_position(&o.rotation, x, r_fixed);
  MT_CHECK(myb_gravity_field_eval(reference, work, r_fixed, 0, -1, NULL, g_fixed) == MYB_OK);
  myb_rotation_apply_position(&inverse, g_fixed, by_hand);
  check_same3(a, by_hand);
  /* a part of the degrees */
  term.n_min = 2;
  term.n_max = 3;
  MT_CHECK(myb_term_field(&term, t, x, v, a) == MYB_OK);
  MT_CHECK(myb_gravity_field_eval(reference, work, r_fixed, 2, 3, NULL, g_fixed) == MYB_OK);
  myb_rotation_apply_position(&inverse, g_fixed, by_hand);
  check_same3(a, by_hand);

  /* tides: set the time, apply, evaluate — by hand on the reference field */
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/grav/iers2010-tides.tk", &diag) == MYB_OK);
  myb_tides* tides = NULL;
  MT_CHECK(myb_tides_create(pool, "IERS10_TIDE", NULL, &tides, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  const int flags = MYB_TIDES_SOLID | MYB_TIDES_FREQUENCY | MYB_TIDES_POLE_SOLID;
  myb_field_term tidal = {.env = env, .body = planet, .field = field, .work = work, .n_min = 0,
                          .n_max = -1, .update = 1, .tides = tides, .tide_flags = flags,
                          .permanent = MYB_PERMANENT_TIDE_KEEP, .raisers = {i_sun, i_moon},
                          .raiser_count = 2};
  MT_CHECK(myb_term_field(&tidal, t, x, v, a) == MYB_OK);
  double s_sun[6];
  double s_moon[6];
  mover twin_sun = sun;
  mover twin_moon = moon;
  MT_CHECK(mover_state(&twin_sun, when, s_sun) == MYB_OK);
  MT_CHECK(mover_state(&twin_moon, when, s_moon) == MYB_OK);
  myb_tide_raiser raisers[2] = {{.gm = b_sun.gm}, {.gm = b_moon.gm}};
  myb_rotation_apply_position(&o.rotation, s_sun, raisers[0].r);
  myb_rotation_apply_position(&o.rotation, s_moon, raisers[1].r);
  const myb_tides_input input = {.raisers = raisers, .count = 2, .t = o.ut, .xp = o.xp,
                                 .yp = o.yp, .flags = flags,
                                 .permanent = MYB_PERMANENT_TIDE_KEEP};
  myb_gravity_field_set_time(reference, when);
  MT_CHECK(myb_tides_apply(tides, reference, &input) == MYB_OK);
  MT_CHECK(myb_gravity_field_eval(reference, work, r_fixed, 0, -1, NULL, g_fixed) == MYB_OK);
  myb_rotation_apply_position(&inverse, g_fixed, by_hand);
  check_same3(a, by_hand);
  double static_part[3];
  myb_gravity_field_set_time(reference, when);
  MT_CHECK(myb_gravity_field_eval(reference, work, r_fixed, 0, -1, NULL, g_fixed) == MYB_OK);
  myb_rotation_apply_position(&inverse, g_fixed, static_part);
  printf("field: the tides change the acceleration by %.1e of it\n", rel3(a, static_part));
  MT_CHECK(rel3(a, static_part) > 1e-9);
  /* the same instant again, and after another one: the same bits */
  double again[3];
  MT_CHECK(myb_term_field(&tidal, t, x, v, again) == MYB_OK);
  check_same3(again, a);
  MT_CHECK(myb_term_field(&tidal, t + 60.0, x, v, again) == MYB_OK);
  MT_CHECK(again[0] != a[0]);
  MT_CHECK(myb_term_field(&tidal, t, x, v, again) == MYB_OK);
  check_same3(again, a);
  myb_env_destroy(env);

  /* a body away from the origin: its field at the object minus its field at the origin */
  myb_gravity_field_set_time(field, when);
  myb_gravity_field_set_time(reference, when);
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  mover near = {.radius = 2.4e7, .rate = 3e-5, .phase = 2.0, .inclination = 0.5};
  size_t other = 0;
  const myb_env_body b_central = {.gm = 1.0e15};
  const myb_env_body b_other = {.gm = gm, .state = mover_state, .state_context = &near,
                                .orientation = spinner_orientation,
                                .orientation_context = &spin};
  MT_CHECK(myb_env_add_body(env, &b_central, NULL) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_other, &other) == MYB_OK);
  myb_field_term away = {.env = env, .body = other, .field = field, .work = work, .n_min = 0,
                         .n_max = -1};
  MT_CHECK(myb_term_field(&away, t, x, v, a) == MYB_OK);
  double s[6];
  mover twin_near = near;
  MT_CHECK(mover_state(&twin_near, when, s) == MYB_OK);
  const double r[3] = {x[0] - s[0], x[1] - s[1], x[2] - s[2]};
  const double r0[3] = {-s[0], -s[1], -s[2]};
  double direct[3];
  double indirect[3];
  myb_rotation_apply_position(&o.rotation, r, r_fixed);
  MT_CHECK(myb_gravity_field_eval(reference, work, r_fixed, 0, -1, NULL, g_fixed) == MYB_OK);
  myb_rotation_apply_position(&inverse, g_fixed, direct);
  myb_rotation_apply_position(&o.rotation, r0, r_fixed);
  MT_CHECK(myb_gravity_field_eval(reference, work, r_fixed, 0, -1, NULL, g_fixed) == MYB_OK);
  myb_rotation_apply_position(&inverse, g_fixed, indirect);
  for (int i = 0; i < 3; ++i) {
    by_hand[i] = direct[i] - indirect[i];
  }
  check_same3(a, by_hand);
  /* degree 0 alone is the point-mass perturbation */
  away.n_max = 0;
  MT_CHECK(myb_term_field(&away, t, x, v, a) == MYB_OK);
  double third[3] = {0.0, 0.0, 0.0};
  myb_force_third_body(gm, s, x, third);
  printf("field: degree 0 of a body away from the origin against the third-body form %.1e\n",
         rel3(a, third));
  MT_CHECK(rel3(a, third) < 1e-14);
  myb_env_destroy(env);

  /* an inertial origin: no indirect term */
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_other, &other) == MYB_OK);
  away.env = env;
  away.body = other;
  away.n_max = -1;
  MT_CHECK(myb_term_field(&away, t, x, v, a) == MYB_OK);
  check_same3(a, direct);
  myb_env_destroy(env);

  myb_tides_destroy(tides);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(reference);
  myb_gravity_field_destroy(field);
}

/* ---- Radiation --------------------------------------------------------------------------------- */

static void test_radiation(void) {
  const myb_time epoch = myb_time_from_seconds(8.0e8);
  const myb_env_options options = {.epoch = epoch};
  const double t = 500.0;
  const myb_time when = myb_time_add(epoch, t);
  const double au = 149597870700.0;
  const double pi = 3.141592653589793;
  mover sun = {.radius = au}; /* at rest on the x axis */
  mover moon = {.radius = 3.8e8, .phase = 2.0, .inclination = 0.3};
  spinner spin = {.angles = {0.2, 0.4, 1.0}, .rate = 7.3e-5};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  size_t planet = 0;
  size_t i_sun = 0;
  size_t i_moon = 0;
  const myb_env_body b_planet = {.gm = 3.986e14, .orientation = spinner_orientation,
                                 .orientation_context = &spin, .radius = 6378137.0};
  const myb_env_body b_sun = {.gm = 1.327e20, .state = mover_state, .state_context = &sun,
                              .radius = 6.957e8, .luminosity = 4.0 * pi * au * au * 1361.0,
                              .limb = {0.6, 0.0}};
  const myb_env_body b_moon = {.gm = 4.9e12, .state = mover_state, .state_context = &moon,
                               .radius = 1737400.0};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_sun, &i_sun) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_moon, &i_moon) == MYB_OK);
  myb_object object = {.mass = 406.965, .area = pi * 0.3 * 0.3, .cr = 1.13, .cd = 2.2};
  myb_radiation_term term = {.env = env, .source = i_sun, .object = &object, .c = c_light,
                             .occulters = {planet, i_moon}, .occulter_count = 2};
  double s_sun[6];
  double s_moon[6];
  mover twin_sun = sun;
  mover twin_moon = moon;
  MT_CHECK(mover_state(&twin_sun, when, s_sun) == MYB_OK);
  MT_CHECK(mover_state(&twin_moon, when, s_moon) == MYB_OK);
  const double v[3] = {0.0, 5.7e3, 0.0};
  /* the Sun is on +x: lit on the day side, the umbra behind the planet, the penumbra at its
   * edge (12 270 km behind the centre the umbra ends 6 322 km from the axis and the penumbra
   * 6 436 km from it) */
  const double places[4][3] = {{1.227e7, 1.0e5, 2.0e5},
                               {-1.227e7, 1.0e5, 2.0e5},
                               {-1.227e7, 6.38e6, 0.0},
                               {-1.227e7, 6.9e6, 0.0}};
  double fractions[4];
  for (int k = 0; k < 4; ++k) {
    const double* x = places[k];
    double a[3];
    MT_CHECK(myb_term_radiation(&term, t, x, v, a) == MYB_OK);
    const double d[3] = {x[0] - s_sun[0], x[1] - s_sun[1], x[2] - s_sun[2]};
    const double dn = myb_vec3_norm(d);
    const double u[3] = {d[0] / dn, d[1] / dn, d[2] / dn};
    const double to_source[3] = {-d[0], -d[1], -d[2]};
    const double to_planet[3] = {0.0 - x[0], 0.0 - x[1], 0.0 - x[2]};
    const double to_moon[3] = {s_moon[0] - x[0], s_moon[1] - x[1], s_moon[2] - x[2]};
    double visible = 1.0;
    visible *= myb_shadow_fraction(to_source, b_sun.radius, to_planet, b_planet.radius, b_sun.limb);
    visible *= myb_shadow_fraction(to_source, b_sun.radius, to_moon, b_moon.radius, b_sun.limb);
    fractions[k] = visible;
    const double pressure = myb_radiation_flux(b_sun.luminosity, dn) / c_light;
    double by_hand[3] = {0.0, 0.0, 0.0};
    myb_force_radiation_sphere(visible * pressure, object.cr, object.area / object.mass, u,
                               by_hand);
    check_same3(a, by_hand);
  }
  printf("radiation: visible fractions %.3f (lit), %.3f (umbra), %.3f (penumbra), %.3f (past it)\n",
         fractions[0], fractions[1], fractions[2], fractions[3]);
  MT_CHECK(fractions[0] == 1.0 && fractions[1] == 0.0 && fractions[3] == 1.0);
  MT_CHECK(fractions[2] > 0.01 && fractions[2] < 0.99);
  /* size: 1361 W/m^2 on LAGEOS */
  double a[3];
  MT_CHECK(myb_term_radiation(&term, t, places[0], v, a) == MYB_OK);
  const double expected = 1361.0 / c_light * 1.13 * object.area / object.mass;
  MT_CHECK(fabs(myb_vec3_norm(a) / expected - 1.0) < 2e-4); /* 12 270 km closer than 1 au */
  /* without occulters */
  myb_radiation_term open = term;
  open.occulter_count = 0;
  MT_CHECK(myb_term_radiation(&open, t, places[1], v, a) == MYB_OK);
  MT_CHECK(fabs(myb_vec3_norm(a) / expected - 1.0) < 2e-4);
  /* errors */
  myb_radiation_term bad = term;
  bad.occulters[1] = i_sun;
  MT_CHECK(myb_term_radiation(&bad, t, places[0], v, a) == MYB_E_DOMAIN);
  bad = term;
  bad.c = 0.0;
  MT_CHECK(myb_term_radiation(&bad, t, places[0], v, a) == MYB_E_DOMAIN);
  object.mass = 0.0;
  MT_CHECK(myb_term_radiation(&term, t, places[0], v, a) == MYB_E_DOMAIN);
  object.mass = 406.965;

  /* the planet's own radiation, in its body-fixed frame */
  static const myb_planet_radiation knocke = {.radius = 6378137.0,
                                              .albedo = {0.34, 0.0, 0.10, 0.0, 0.29},
                                              .emissivity = {0.68, 0.0, -0.07, 0.0, -0.18},
                                              .omega = 1.99106e-7};
  myb_planet_radiation_term planet_term = {.env = env, .planet = planet, .source = i_sun,
                                           .model = &knocke, .object = &object, .phase0 = 0.7,
                                           .rings = 6};
  const double* x = places[0];
  MT_CHECK(myb_term_planet_radiation(&planet_term, t, x, v, a) == MYB_OK);
  myb_orientation o;
  spinner twin = spin;
  MT_CHECK(spinner_orientation(&twin, when, &o) == MYB_OK);
  myb_rotation inverse;
  myb_rotation_invert(&o.rotation, &inverse);
  double r_fixed[3];
  double sun_fixed[3];
  myb_rotation_apply_position(&o.rotation, x, r_fixed);
  myb_rotation_apply_position(&o.rotation, s_sun, sun_fixed);
  const double flux = myb_radiation_flux(b_sun.luminosity, myb_vec3_norm(s_sun));
  double a_fixed[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_planet_radiation(&knocke, r_fixed, sun_fixed, flux, 0.7 + knocke.omega * t,
                                      object.cr * (object.area / object.mass), 6, a_fixed) ==
           MYB_OK);
  double by_hand[3];
  myb_rotation_apply_position(&inverse, a_fixed, by_hand);
  check_same3(a, by_hand);
  MT_CHECK(myb_vec3_norm(a) > 1e-11 && myb_vec3_dot(a, x) > 0.0); /* away from the planet */
  planet_term.rings = 0;
  MT_CHECK(myb_term_planet_radiation(&planet_term, t, x, v, a) == MYB_E_DOMAIN);
  myb_env_destroy(env);

  /* occulters need the disk of the source */
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  myb_env_body point_sun = b_sun;
  point_sun.radius = 0.0;
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &point_sun, &i_sun) == MYB_OK);
  myb_radiation_term point = {.env = env, .source = i_sun, .object = &object, .c = c_light,
                              .occulters = {planet}, .occulter_count = 1};
  MT_CHECK(myb_term_radiation(&point, t, places[0], v, a) == MYB_E_DOMAIN);
  point.occulter_count = 0;
  MT_CHECK(myb_term_radiation(&point, t, places[0], v, a) == MYB_OK);
  myb_env_destroy(env);
}

/* ---- Atmosphere --------------------------------------------------------------------------------- */

static void test_drag(void) {
  /* the layered atmosphere */
  myb_exponential_atmosphere air = {.ellipsoid = {.a = 6378137.0, .f = 0.0},
                                          .count = 3,
                                          .base = {0.0, 25000.0, 500000.0},
                                          .density = {1.225, 3.899e-2, 6.967e-13},
                                          .scale = {7249.0, 6349.0, 63822.0}};
  const myb_time any = myb_time_from_seconds(0.0);
  const double heights[6] = {-100.0, 0.0, 12000.0, 25000.0, 499999.0, 800000.0};
  const int layers[6] = {0, 0, 0, 1, 1, 2};
  for (int k = 0; k < 6; ++k) {
    const double r[3] = {0.0, air.ellipsoid.a + heights[k], 0.0};
    double h = 0.0;
    double lon = 0.0;
    double lat = 0.0;
    MT_CHECK(myb_geodetic_from_cartesian(&air.ellipsoid, r, &lon, &lat, &h) == MYB_OK);
    double rho = -1.0;
    MT_CHECK(myb_density_exponential(&air, any, r, &rho) == MYB_OK);
    const int n = layers[k];
    MT_CHECK_SAME(rho, air.density[n] * myb_exp(-(h - air.base[n]) / air.scale[n]));
  }
  double rho = 0.0;
  const double place[3] = {7.0e6, 0.0, 0.0};
  myb_exponential_atmosphere broken = air;
  broken.count = 0;
  MT_CHECK(myb_density_exponential(&broken, any, place, &rho) == MYB_E_DOMAIN);
  broken = air;
  broken.base[2] = 20000.0;
  MT_CHECK(myb_density_exponential(&broken, any, place, &rho) == MYB_E_DOMAIN);
  broken = air;
  broken.scale[1] = 0.0;
  MT_CHECK(myb_density_exponential(&broken, any, place, &rho) == MYB_E_DOMAIN);
  broken = air;
  broken.ellipsoid.a = 0.0;
  MT_CHECK(myb_density_exponential(&broken, any, place, &rho) == MYB_E_DOMAIN);

  /* drag in the co-rotating atmosphere of a planet spinning about z */
  const myb_time epoch = myb_time_from_seconds(8.0e8);
  const myb_env_options options = {.epoch = epoch};
  const double t = 250.0;
  const myb_time when = myb_time_add(epoch, t);
  spinner spin = {.angles = {0.0, 0.0, 0.3}, .rate = 7.292115e-5};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  size_t planet = 0;
  const myb_env_body b_planet = {.gm = 3.986e14, .orientation = spinner_orientation,
                                 .orientation_context = &spin};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  myb_object object = {.mass = 500.0, .area = 2.0, .cr = 1.3, .cd = 2.2};
  myb_exponential_atmosphere copy = air;
  myb_drag_term term = {.env = env, .body = planet, .density = myb_density_exponential,
                        .density_context = &copy, .object = &object};
  const double x[3] = {4.0e6, 5.0e6, 2.6e6};
  const double v[3] = {-5.5e3, 3.4e3, 3.9e3};
  double a[3];
  MT_CHECK(myb_term_drag(&term, t, x, v, a) == MYB_OK);
  /* by hand with the functions of the library */
  myb_orientation o;
  spinner twin = spin;
  MT_CHECK(spinner_orientation(&twin, when, &o) == MYB_OK);
  myb_rotation inverse;
  myb_rotation_invert(&o.rotation, &inverse);
  const double y[6] = {x[0], x[1], x[2], v[0], v[1], v[2]};
  double fixed[6];
  myb_rotation_apply_state(&o.rotation, y, fixed);
  MT_CHECK(myb_density_exponential(&copy, when, fixed, &rho) == MYB_OK);
  double w[3];
  myb_rotation_apply_position(&inverse, fixed + 3, w);
  double by_hand[3] = {0.0, 0.0, 0.0};
  myb_force_drag(rho, object.cd * (object.area / object.mass), w, by_hand);
  check_same3(a, by_hand);
  /* by the formula: v - omega x r */
  const double omega[3] = {0.0, 0.0, spin.rate};
  double wxr[3];
  myb_vec3_cross(omega, x, wxr);
  const double relative[3] = {v[0] - wxr[0], v[1] - wxr[1], v[2] - wxr[2]};
  const double speed = myb_vec3_norm(relative);
  const double k = -0.5 * rho * object.cd * object.area / object.mass * speed;
  const double formula[3] = {k * relative[0], k * relative[1], k * relative[2]};
  printf("drag: against -1/2 rho cd A/m |w| w with w = v - omega x r: %.1e\n",
         rel3(a, formula));
  MT_CHECK(rel3(a, formula) < 1e-14);
  MT_CHECK(rho > 0.0 && myb_vec3_dot(a, relative) < 0.0);
  object.cd = (double)NAN;
  MT_CHECK(myb_term_drag(&term, t, x, v, a) == MYB_E_DOMAIN);
  object.cd = 2.2;
  copy.count = 40;
  MT_CHECK(myb_term_drag(&term, t, x, v, a) == MYB_E_DOMAIN);
  myb_env_destroy(env);
}

/* ---- Relativity and empirical accelerations ------------------------------------------------------ */

static void test_relativity(void) {
  const myb_time epoch = myb_time_from_seconds(8.0e8);
  const myb_env_options options = {.epoch = epoch};
  const double t = 777.0;
  const myb_time when = myb_time_add(epoch, t);
  const myb_ppn ppn = {c_light, 1.0, 1.0};
  mover sun = {.radius = 1.5e11, .rate = 2e-7, .phase = 1.0, .inclination = 0.4};
  spinner spin = {.angles = {0.2, 0.4, 1.0}, .rate = 7.3e-5};
  const double x[3] = {7.0e6, -1.2e6, 2.5e6};
  const double v[3] = {1.0e3, 7.2e3, -0.5e3};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  size_t planet = 0;
  size_t i_sun = 0;
  const myb_env_body b_planet = {.gm = 3.986004418e14, .orientation = spinner_orientation,
                                 .orientation_context = &spin, .spin = {0.0, 0.0, 9.8e8}};
  const myb_env_body b_sun = {.gm = 1.32712440018e20, .state = mover_state,
                              .state_context = &sun};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_sun, &i_sun) == MYB_OK);
  myb_relativity_term term = {.env = env, .ppn = ppn, .central = planet, .external = i_sun};
  double a[3];
  double by_hand[3];
  MT_CHECK(myb_term_schwarzschild(&term, t, x, v, a) == MYB_OK);
  MT_CHECK(myb_rel_schwarzschild(&ppn, b_planet.gm, x, v, by_hand) == MYB_OK);
  check_same3(a, by_hand);
  myb_orientation o;
  spinner twin = spin;
  MT_CHECK(spinner_orientation(&twin, when, &o) == MYB_OK);
  myb_rotation inverse;
  myb_rotation_invert(&o.rotation, &inverse);
  double j[3];
  myb_rotation_apply_position(&inverse, b_planet.spin, j);
  MT_CHECK(myb_term_lense_thirring(&term, t, x, v, a) == MYB_OK);
  MT_CHECK(myb_rel_lense_thirring(&ppn, b_planet.gm, j, x, v, by_hand) == MYB_OK);
  check_same3(a, by_hand);
  double s[6];
  mover twin_sun = sun;
  MT_CHECK(mover_state(&twin_sun, when, s) == MYB_OK);
  const double big_r[3] = {0.0 - s[0], 0.0 - s[1], 0.0 - s[2]};
  const double big_v[3] = {0.0 - s[3], 0.0 - s[4], 0.0 - s[5]};
  MT_CHECK(myb_term_de_sitter(&term, t, x, v, a) == MYB_OK);
  MT_CHECK(myb_rel_de_sitter(&ppn, b_sun.gm, big_r, big_v, v, by_hand) == MYB_OK);
  check_same3(a, by_hand);
  term.external = planet;
  MT_CHECK(myb_term_de_sitter(&term, t, x, v, a) == MYB_E_DOMAIN);
  myb_env_destroy(env);

  /* a central body that moves: the terms use the state relative to it */
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  mover wander = {.radius = 9.0e8, .rate = 1e-6, .phase = 0.3, .inclination = 0.2};
  myb_env_body moving = b_planet;
  moving.state = mover_state;
  moving.state_context = &wander;
  MT_CHECK(myb_env_add_body(env, &moving, &planet) == MYB_OK);
  term.env = env;
  term.central = planet;
  mover twin_wander = wander;
  MT_CHECK(mover_state(&twin_wander, when, s) == MYB_OK);
  const double xb[3] = {s[0] + x[0], s[1] + x[1], s[2] + x[2]};
  const double vb[3] = {s[3] + v[0], s[4] + v[1], s[5] + v[2]};
  const double xr[3] = {xb[0] - s[0], xb[1] - s[1], xb[2] - s[2]};
  const double vr[3] = {vb[0] - s[3], vb[1] - s[4], vb[2] - s[5]};
  MT_CHECK(myb_term_schwarzschild(&term, t, xb, vb, a) == MYB_OK);
  MT_CHECK(myb_rel_schwarzschild(&ppn, b_planet.gm, xr, vr, by_hand) == MYB_OK);
  check_same3(a, by_hand);
  myb_env_destroy(env);
}

static void test_empirical(void) {
  myb_empirical_term term = {.constant = {1e-10, -2e-10, 3e-11},
                                   .cosine = {4e-11, 5e-10, -6e-10},
                                   .sine = {-7e-11, 8e-11, 9e-10}};
  const double x[3] = {7.0e6, -1.2e6, 2.5e6};
  const double v[3] = {1.0e3, 7.2e3, -0.5e3};
  double a[3];
  MT_CHECK(myb_term_empirical(&term, 0.0, x, v, a) == MYB_OK);
  /* u from the elements of the orbit */
  myb_pv pv;
  memcpy(pv.r, x, sizeof pv.r);
  memcpy(pv.v, v, sizeof pv.v);
  myb_elements el;
  MT_CHECK(myb_elements_from_pv(3.986004418e14, &pv, &el) == MYB_OK);
  const double u = el.argp + el.nu;
  double rtn[3];
  for (int i = 0; i < 3; ++i) {
    rtn[i] = term.constant[i] + term.cosine[i] * cos(u) + term.sine[i] * sin(u);
  }
  double expected[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_empirical_rtn(x, v, rtn, expected) == MYB_OK);
  printf("empirical: against the argument of latitude of the elements %.1e\n", rel3(a, expected));
  MT_CHECK(rel3(a, expected) < 1e-14);
  /* an orbit in the xy plane counts u from the x axis, in the sense of the motion */
  myb_empirical_term once = {.cosine = {1.0, 0.0, 0.0}, .sine = {0.0, 1.0, 0.0}};
  const double xe[3] = {0.0, 8.0e6, 0.0};
  const double prograde[3] = {-7.0e3, 0.0, 0.0};
  const double retrograde[3] = {7.0e3, 0.0, 0.0};
  MT_CHECK(myb_term_empirical(&once, 0.0, xe, prograde, a) == MYB_OK);
  /* u = 90 deg: cos u = 0, sin u = 1: a = 1 along the transverse direction -x */
  MT_CHECK(fabs(a[0] + 1.0) < 1e-15 && fabs(a[1]) < 1e-15 && a[2] == 0.0);
  MT_CHECK(myb_term_empirical(&once, 0.0, xe, retrograde, a) == MYB_OK);
  /* retrograde: u = -90 deg counted along the motion: sin u = -1, transverse +x */
  MT_CHECK(fabs(a[0] + 1.0) < 1e-15 && fabs(a[1]) < 1e-15 && a[2] == 0.0);
  MT_CHECK(myb_term_empirical(&once, 0.0, xe, xe, a) == MYB_E_DOMAIN);
}

/* ---- Analytic Jacobians --------------------------------------------------------------------------- */

/* The Jacobian a term gives against central differences of the term (the force model computes
 * both); returns the largest difference of d a / d x and of d a / d v relative to the largest
 * element of each. */
static void jacobian_errors(myb_term_function function, myb_term_jacobian jacobian, void* term,
                            double t, const double x[3], const double v[3], double* ex,
                            double* ev) {
  myb_force_model* analytic = NULL;
  myb_force_model* numeric = NULL;
  MT_CHECK(myb_force_model_create(NULL, &analytic) == MYB_OK);
  MT_CHECK(myb_force_model_create(NULL, &numeric) == MYB_OK);
  MT_CHECK(myb_force_model_add(analytic, "term", function, jacobian, term, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(numeric, "term", function, NULL, term, NULL) == MYB_OK);
  double ax[9];
  double av[9];
  double nx[9];
  double nv[9];
  MT_CHECK(myb_force_model_jacobian(analytic, t, x, v, ax, av, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_jacobian(numeric, t, x, v, nx, nv, NULL) == MYB_OK);
  double sx = 0.0;
  double sv = 0.0;
  double dx = 0.0;
  double dv = 0.0;
  for (int i = 0; i < 9; ++i) {
    sx = fmax(sx, fabs(nx[i]));
    sv = fmax(sv, fabs(nv[i]));
    dx = fmax(dx, fabs(ax[i] - nx[i]));
    dv = fmax(dv, fabs(av[i] - nv[i]));
  }
  *ex = sx > 0.0 ? dx / sx : dx;
  *ev = sv > 0.0 ? dv / sv : dv;
  myb_force_model_destroy(analytic);
  myb_force_model_destroy(numeric);
}

static void test_jacobians(void) {
  const myb_time epoch = myb_time_from_seconds(8.0e8);
  const myb_env_options options = {.epoch = epoch};
  const double au = 149597870700.0;
  const double pi = 3.141592653589793;
  const double t = 640.0;
  const double x[3] = {4.0e6, 5.0e6, 2.6e6};
  const double v[3] = {-5.5e3, 3.4e3, 3.9e3};
  spinner spin = {.angles = {0.2, 0.4, 1.0}, .rate = 7.292115e-5, .ut_shift = -0.3, .xp = 6e-7,
                  .yp = 1.6e-6};
  mover sun = {.radius = au};
  mover moon = {.radius = 3.8e8, .rate = 2.7e-6, .phase = 2.0, .inclination = 0.3};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  size_t planet = 0;
  size_t i_sun = 0;
  size_t i_moon = 0;
  const myb_env_body b_planet = {.gm = 3.986004415e14, .orientation = spinner_orientation,
                                 .orientation_context = &spin, .radius = 6378137.0,
                                 .spin = {1.0e8, -2.0e8, 9.8e8}};
  const myb_env_body b_sun = {.gm = 1.32712440018e20, .state = mover_state,
                              .state_context = &sun, .radius = 6.957e8,
                              .luminosity = 4.0 * pi * au * au * 1361.0, .limb = {0.6, 0.0}};
  const myb_env_body b_moon = {.gm = 4.9028e12, .state = mover_state, .state_context = &moon,
                               .radius = 1737400.0};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_sun, &i_sun) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_moon, &i_moon) == MYB_OK);
  double ex = 0.0;
  double ev = 0.0;

  /* the field with tides */
  myb_gravity_field* field = NULL;
  myb_gravity_work* work = NULL;
  make_field(b_planet.gm, &field);
  MT_CHECK(myb_gravity_work_create(DEGREE, NULL, &work) == MYB_OK);
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/grav/iers2010-tides.tk", &diag) == MYB_OK);
  myb_tides* tides = NULL;
  MT_CHECK(myb_tides_create(pool, "IERS10_TIDE", NULL, &tides, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_field_term gravity = {.env = env, .body = planet, .field = field, .work = work,
                            .n_min = 0, .n_max = -1, .update = 1, .tides = tides,
                            .tide_flags = MYB_TIDES_SOLID | MYB_TIDES_FREQUENCY,
                            .permanent = MYB_PERMANENT_TIDE_KEEP, .raisers = {i_sun, i_moon},
                            .raiser_count = 2};
  jacobian_errors(myb_term_field, myb_term_field_jacobian, &gravity, t, x, v, &ex, &ev);
  printf("jacobians against central differences: field %.1e", ex);
  MT_CHECK(ex < 1e-9 && ev == 0.0);
  /* the Jacobian alone brings the field to its instant: asked first at a new time */
  double jx[9];
  double jv[9];
  double a[3];
  MT_CHECK(myb_term_field_jacobian(&gravity, t + 5.0, x, v, jx, jv) == MYB_OK);
  MT_CHECK(myb_term_field(&gravity, t + 5.0, x, v, a) == MYB_OK);
  double again[9];
  MT_CHECK(myb_term_field_jacobian(&gravity, t + 5.0, x, v, again, jv) == MYB_OK);
  for (int i = 0; i < 9; ++i) {
    MT_CHECK_SAME(again[i], jx[i]);
  }
  /* symmetric and free of trace in any axes */
  double size = 0.0;
  for (int i = 0; i < 9; ++i) {
    size = fmax(size, fabs(jx[i]));
  }
  MT_CHECK(fabs(jx[1] - jx[3]) < 1e-14 * size && fabs(jx[2] - jx[6]) < 1e-14 * size &&
           fabs(jx[5] - jx[7]) < 1e-14 * size);
  MT_CHECK(fabs(jx[0] + jx[4] + jx[8]) < 1e-14 * size);

  /* radiation: lit, umbra, penumbra */
  myb_object object = {.mass = 406.965, .area = pi * 0.3 * 0.3, .cr = 1.13, .cd = 2.2};
  myb_radiation_term light = {.env = env, .source = i_sun, .object = &object, .c = c_light,
                              .occulters = {planet, i_moon}, .occulter_count = 2};
  const double lit[3] = {1.227e7, 1.0e5, 2.0e5};
  const double dark[3] = {-1.227e7, 1.0e5, 2.0e5};
  const double edge[3] = {-1.227e7, 6.38e6, 0.0};
  /* the pressure changes over the distance of the Sun: the differences of the force model
   * (steps of 74 m) resolve its gradient to 1e-6 only; steps of 1000 km do to 1e-10 */
  jacobian_errors(myb_term_radiation, myb_term_radiation_jacobian, &light, t, lit, v, &ex, &ev);
  MT_CHECK(ex < 1e-5 && ev == 0.0);
  MT_CHECK(myb_term_radiation_jacobian(&light, t, lit, v, jx, jv) == MYB_OK);
  myb_radiation_term open_light = light;
  open_light.occulter_count = 0;
  ex = 0.0;
  for (int j = 0; j < 3; ++j) {
    double plus[3] = {lit[0], lit[1], lit[2]};
    double minus[3] = {lit[0], lit[1], lit[2]};
    plus[j] += 1.0e6;
    minus[j] -= 1.0e6;
    double ap[3];
    double am[3];
    MT_CHECK(myb_term_radiation(&open_light, t, plus, v, ap) == MYB_OK);
    MT_CHECK(myb_term_radiation(&open_light, t, minus, v, am) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      ex = fmax(ex, fabs(jx[3 * i + j] - (ap[i] - am[i]) / 2.0e6));
    }
  }
  ex /= fabs(jx[0]);
  printf(", light %.1e", ex);
  MT_CHECK(ex < 1e-9);
  MT_CHECK(myb_term_radiation_jacobian(&light, t, dark, v, jx, jv) == MYB_OK);
  for (int i = 0; i < 9; ++i) {
    MT_CHECK(jx[i] == 0.0 && jv[i] == 0.0);
  }
  jacobian_errors(myb_term_radiation, myb_term_radiation_jacobian, &light, t, edge, v, &ex, &ev);
  printf(" (penumbra %.1e)", ex);
  MT_CHECK(ex < 1e-4 && ev == 0.0);
  MT_CHECK(myb_term_radiation_jacobian(&light, t, edge, v, jx, jv) == MYB_OK);
  MT_CHECK(fabs(jx[1]) > 1e-15); /* the shadow's edge: 3e-14 s^-2 against 5e-20 in full light */

  /* drag */
  myb_exponential_atmosphere air = {.ellipsoid = {.a = 6378137.0, .f = 1.0 / 298.257223563},
                                    .count = 2,
                                    .base = {200000.0, 500000.0},
                                    .density = {2.5e-10, 6.967e-13},
                                    .scale = {37105.0, 63822.0}};
  myb_drag_term drag = {.env = env, .body = planet, .density = myb_density_exponential,
                        .density_context = &air, .object = &object};
  jacobian_errors(myb_term_drag, myb_term_drag_jacobian, &drag, t, x, v, &ex, &ev);
  printf(", drag %.1e and %.1e", ex, ev);
  MT_CHECK(ex < 1e-6 && ev < 1e-8);

  /* relativity */
  mover star = {.radius = au, .rate = 2e-7, .phase = 1.0, .inclination = 0.4};
  size_t i_star = 0;
  const myb_env_body b_star = {.gm = 1.32712440018e20, .state = mover_state,
                               .state_context = &star};
  MT_CHECK(myb_env_add_body(env, &b_star, &i_star) == MYB_OK);
  myb_relativity_term rel = {.env = env, .ppn = {c_light, 1.1, 0.9}, .central = planet,
                             .external = i_star};
  jacobian_errors(myb_term_schwarzschild, myb_term_schwarzschild_jacobian, &rel, t, x, v, &ex,
                  &ev);
  printf(", Schwarzschild %.1e and %.1e", ex, ev);
  MT_CHECK(ex < 1e-8 && ev < 1e-8);
  jacobian_errors(myb_term_lense_thirring, myb_term_lense_thirring_jacobian, &rel, t, x, v, &ex,
                  &ev);
  printf(", Lense-Thirring %.1e and %.1e", ex, ev);
  MT_CHECK(ex < 1e-8 && ev < 1e-8);
  jacobian_errors(myb_term_de_sitter, myb_term_de_sitter_jacobian, &rel, t, x, v, &ex, &ev);
  printf(", de Sitter %.1e\n", ev);
  MT_CHECK(ex == 0.0 && ev < 1e-8);

  /* a parameter bound to its term: the same partial from one term instead of the model */
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  size_t i_field = 0;
  size_t i_light = 0;
  size_t i_cr = 0;
  MT_CHECK(myb_force_model_add(model, "field", myb_term_field, myb_term_field_jacobian,
                               &gravity, &i_field) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "light", myb_term_radiation, myb_term_radiation_jacobian,
                               &light, &i_light) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(model, "cr", &object.cr, 0.01, &i_cr) == MYB_OK);
  double free_partial[3];
  double bound_partial[3];
  MT_CHECK(myb_force_model_jacobian(model, t, lit, v, jx, jv, free_partial) == MYB_OK);
  MT_CHECK(myb_force_model_bind_parameter(model, i_cr, i_light) == MYB_OK);
  MT_CHECK(myb_force_model_jacobian(model, t, lit, v, jx, jv, bound_partial) == MYB_OK);
  MT_CHECK(myb_term_radiation(&light, t, lit, v, a) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    /* the acceleration is linear in cr: the partial is a / cr */
    MT_CHECK(fabs(bound_partial[i] - a[i] / object.cr) < 1e-13 * fabs(a[0] / object.cr));
    /* differences of the whole model lose the digits the field takes: 1e-5 here */
    MT_CHECK(fabs(free_partial[i] - a[i] / object.cr) < 1e-3 * fabs(a[0] / object.cr));
  }
  myb_force_model_set_enabled(model, i_light, 0);
  MT_CHECK(myb_force_model_jacobian(model, t, lit, v, jx, jv, bound_partial) == MYB_OK);
  MT_CHECK(bound_partial[0] == 0.0 && bound_partial[1] == 0.0 && bound_partial[2] == 0.0);
  MT_CHECK(myb_force_model_bind_parameter(model, i_cr, 7) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_model_bind_parameter(model, 3, i_light) == MYB_E_DOMAIN);
  MT_CHECK(myb_term_no_jacobian(NULL, t, x, v, jx, jv) == MYB_OK);
  for (int i = 0; i < 9; ++i) {
    MT_CHECK(jx[i] == 0.0 && jv[i] == 0.0);
  }
  myb_force_model_destroy(model);
  myb_tides_destroy(tides);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
  myb_env_destroy(env);
}

/* The SPK body and the fixed orientation as providers. */
static void test_adapters(void) {
  myb_rotation tilt;
  const double angles[3] = {0.3, 0.2, 0.1};
  const double rates[3] = {0.0, 0.0, 0.0};
  myb_rotation_euler313(angles, rates, &tilt);
  myb_orientation o;
  const myb_time when = myb_time_from_seconds(123.0);
  MT_CHECK(myb_orientation_fixed(&tilt, when, &o) == MYB_OK);
  MT_CHECK(memcmp(&o.rotation, &tilt, sizeof tilt) == 0);
  MT_CHECK(myb_time_compare(o.ut, when) == 0 && o.xp == 0.0 && o.yp == 0.0);

  myb_spk* spk = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  MT_CHECK(myb_spk_load_file(spk, MYB_TEST_DATA_DIR "/spk/t02_03.bsp", &diag) == MYB_OK);
  myb_spk_segment_info info;
  myb_spk_segment(spk, 0, &info);
  const myb_time middle = myb_time_from_seconds(0.5 * (info.start + info.stop));
  double expected[6];
  double state[6];
  int32_t frame = 0;
  MT_CHECK(myb_spk_state(spk, info.target, info.center, middle, expected, &frame) == MYB_OK);
  myb_spk_body body = {.spk = spk, .target = info.target, .center = info.center};
  MT_CHECK(myb_spk_body_state(&body, middle, state) == MYB_OK);
  for (int i = 0; i < 6; ++i) {
    MT_CHECK_SAME(state[i], expected[i]);
  }
  body.frame = frame;
  MT_CHECK(myb_spk_body_state(&body, middle, state) == MYB_OK);
  body.frame = frame + 1;
  MT_CHECK(myb_spk_body_state(&body, middle, state) == MYB_E_FORMAT);
  body.frame = 0;
  MT_CHECK(myb_spk_body_state(&body, myb_time_from_seconds(info.stop + 1e9), state) ==
           MYB_E_RANGE);
  myb_spk_destroy(spk);
}

/* A copy of an environment gives the same values and keeps its own. */
static void test_clone(void) {
  const myb_time epoch = myb_time_from_seconds(1000.0);
  const myb_env_options options = {.epoch = epoch};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  mover moon = {.radius = 3.8e8, .rate = 2.7e-6, .phase = 0.4, .inclination = 0.3};
  spinner spin = {.angles = {0.2, 0.4, 1.0}, .rate = 7.3e-5};
  size_t planet = 0;
  size_t i_moon = 0;
  const myb_env_body b_planet = {.id = 399, .gm = 3.986e14, .orientation = spinner_orientation,
                                 .orientation_context = &spin, .radius = 6.4e6};
  const myb_env_body b_moon = {.id = 301, .gm = 4.9e12, .state = mover_state,
                               .state_context = &moon};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_moon, &i_moon) == MYB_OK);
  double first[6];
  MT_CHECK(myb_env_state(env, 500.0, i_moon, first) == MYB_OK);
  myb_env* copy = NULL;
  MT_CHECK(myb_env_clone(env, NULL, &copy) == MYB_OK);
  MT_CHECK(myb_env_body_count(copy) == 2 && myb_env_gm(copy, i_moon) == 4.9e12);
  size_t found = 9;
  MT_CHECK(myb_env_find(copy, 301, &found) == MYB_OK && found == i_moon);
  /* the copy starts with nothing computed: asking it calls the provider again */
  const int calls = moon.calls;
  double second[6];
  MT_CHECK(myb_env_state(copy, 500.0, i_moon, second) == MYB_OK);
  MT_CHECK(moon.calls == calls + 1);
  MT_CHECK(memcmp(first, second, sizeof first) == 0);
  /* each keeps its own instant */
  double other[6];
  MT_CHECK(myb_env_state(copy, 900.0, i_moon, other) == MYB_OK);
  MT_CHECK(myb_env_state(env, 500.0, i_moon, second) == MYB_OK);
  MT_CHECK(moon.calls == calls + 2); /* the original still had 500 s */
  MT_CHECK(memcmp(first, second, sizeof first) == 0 && other[0] != first[0]);
  myb_orientation o1;
  myb_orientation o2;
  MT_CHECK(myb_env_orientation(env, 700.0, planet, &o1) == MYB_OK);
  MT_CHECK(myb_env_orientation(copy, 700.0, planet, &o2) == MYB_OK);
  MT_CHECK(memcmp(o1.rotation.m, o2.rotation.m, sizeof o1.rotation.m) == 0);
  myb_env_destroy(copy);
  myb_env_destroy(env);
}

int main(void) {
  test_scales();
  test_bodies();
  test_point_masses();
  test_field();
  test_radiation();
  test_drag();
  test_relativity();
  test_empirical();
  test_jacobians();
  test_adapters();
  test_clone();
  return mt_finish("test_env");
}
