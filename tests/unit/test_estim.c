/* SPDX-License-Identifier: BSD-3-Clause */
/* Estimation of an orbit (myb/estim.h) on a made-up planet (oblate field, uniform rotation,
 * a layer of atmosphere) with three stations, so nothing is downloaded:
 *  - measurements: derivatives against central differences; the station's state;
 *  - least squares: observations without noise bring a displaced guess back to the truth
 *    (state and drag coefficient); with noise the estimate lies within its own covariance and
 *    the weighted RMS is 1; fixed unknowns, a priori information, an arc backward in time,
 *    observations that do not determine the unknowns;
 *  - the filter: the same noisy observations one by one from a wide covariance end within the
 *    filter's covariance of the truth; process noise; errors. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/myb.h"

static const double gm = 3.986004415e14;
static const double two_pi = 6.283185307179586;

typedef struct spinner {
  double rate;
} spinner;

static myb_status spinner_orientation(void* context, myb_time t, myb_orientation* o) {
  const spinner* s = context;
  const double angles[3] = {0.3, 0.2, 1.0 + s->rate * myb_time_seconds(t)};
  const double rates[3] = {0.0, 0.0, s->rate};
  myb_rotation_euler313(angles, rates, &o->rotation);
  o->ut = t;
  o->xp = 0.0;
  o->yp = 0.0;
  return MYB_OK;
}

/* The scenario: environment, model with the drag coefficient as its parameter, stations. */
typedef struct scenario {
  spinner spin;
  myb_env* env;
  size_t planet;
  myb_gravity_field* field;
  myb_gravity_work* work;
  myb_exponential_atmosphere air;
  myb_object object;
  myb_field_term gravity;
  myb_drag_term drag;
  myb_force_model* model;
  myb_station station[3];
  double y0[6];
  double period;
} scenario;

static void scenario_init(scenario* s) {
  memset(s, 0, sizeof *s);
  s->spin.rate = 7.292115e-5;
  const myb_env_options options = {.epoch = myb_time_from_seconds(0.0)};
  MT_CHECK(myb_env_create(&options, NULL, &s->env) == MYB_OK);
  const myb_env_body planet = {.gm = gm, .orientation = spinner_orientation,
                               .orientation_context = &s->spin};
  MT_CHECK(myb_env_add_body(s->env, &planet, &s->planet) == MYB_OK);
  double c[15] = {0.0};
  double sn[15] = {0.0};
  c[0] = 1.0;
  c[3] = -4.84165e-4;
  c[5] = 2.439e-6;
  sn[5] = -1.400e-6;
  c[6] = 9.57e-7;
  c[10] = 5.4e-7;
  MT_CHECK(myb_gravity_field_create(gm, 6378136.3, 4, c, sn, NULL, &s->field) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(4, NULL, &s->work) == MYB_OK);
  s->air = (myb_exponential_atmosphere){.ellipsoid = {6378137.0, 1.0 / 298.257223563},
                                        .count = 1,
                                        .base = {400000.0},
                                        .density = {3.0e-12},
                                        .scale = {58000.0}};
  s->object = (myb_object){.mass = 200.0, .area = 4.0, .cr = 1.3, .cd = 2.2};
  s->gravity = (myb_field_term){.env = s->env, .body = s->planet, .field = s->field,
                                .work = s->work, .n_min = 0, .n_max = -1};
  s->drag = (myb_drag_term){.env = s->env, .body = s->planet,
                            .density = myb_density_exponential, .density_context = &s->air,
                            .object = &s->object};
  MT_CHECK(myb_force_model_create(NULL, &s->model) == MYB_OK);
  size_t i_drag = 0;
  size_t i_cd = 0;
  MT_CHECK(myb_force_model_add(s->model, "field", myb_term_field, myb_term_field_jacobian,
                               &s->gravity, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(s->model, "drag", myb_term_drag, myb_term_drag_jacobian,
                               &s->drag, &i_drag) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(s->model, "cd", &s->object.cd, 0.01, &i_cd) == MYB_OK);
  MT_CHECK(myb_force_model_bind_parameter(s->model, i_cd, i_drag) == MYB_OK);
  const double sites[3][3] = {{0.3, 0.9, 150.0}, {2.4, -0.5, 2300.0}, {-1.9, 0.2, 10.0}};
  for (int k = 0; k < 3; ++k) {
    MT_CHECK(myb_station_init(&s->station[k], s->env, s->planet, &s->air.ellipsoid, sites[k][0],
                              sites[k][1], sites[k][2]) == MYB_OK);
  }
  const double a = 6378137.0 + 420000.0;
  const myb_elements el = {.p = a * (1.0 - 1e-6), .e = 0.001, .i = 0.9, .raan = 0.4,
                           .argp = 1.1, .nu = 0.3};
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(gm, &el, &pv) == MYB_OK);
  memcpy(s->y0, pv.r, 3 * sizeof(double));
  memcpy(s->y0 + 3, pv.v, 3 * sizeof(double));
  s->period = two_pi * sqrt(a * a * a / gm);
}

static void scenario_destroy(scenario* s) {
  myb_force_model_destroy(s->model);
  myb_gravity_work_destroy(s->work);
  myb_gravity_field_destroy(s->field);
  myb_env_destroy(s->env);
}

/* ---- Measurements ------------------------------------------------------------------------------ */

static void test_measurements(scenario* s) {
  const double t = 1234.5;
  const double y[6] = {5.1e6, 3.9e6, 2.3e6, -4.9e3, 3.2e3, 5.0e3};
  myb_station* st = &s->station[0];
  /* the station: on the ellipsoid, moving with the planet */
  double site[6];
  MT_CHECK(myb_station_state(st, t, site) == MYB_OK);
  myb_orientation o;
  MT_CHECK(myb_env_orientation(s->env, t, s->planet, &o) == MYB_OK);
  double fixed[6];
  myb_rotation_apply_state(&o.rotation, site, fixed);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(fixed[i] - st->position[i]) < 1e-8);
    MT_CHECK(fabs(fixed[3 + i]) < 1e-11); /* at rest in the body-fixed frame */
  }
  const double speed = sqrt(site[3] * site[3] + site[4] * site[4] + site[5] * site[5]);
  MT_CHECK(speed > 100.0 && speed < 500.0);

  const myb_measure_function measures[4] = {myb_measure_range, myb_measure_range_rate,
                                            myb_measure_radec, myb_measure_azel};
  const size_t counts[4] = {1, 1, 2, 2};
  const char* names[4] = {"range", "range rate", "radec", "azel"};
  printf("partials of the measurements against central differences:");
  for (int k = 0; k < 4; ++k) {
    double value[2];
    double partial[12];
    MT_CHECK(measures[k](st, t, y, value, partial) == MYB_OK);
    double worst = 0.0;
    for (int j = 0; j < 6; ++j) {
      const double h = j < 3 ? 10.0 : 0.01;
      double plus[6];
      double minus[6];
      memcpy(plus, y, sizeof plus);
      memcpy(minus, y, sizeof minus);
      plus[j] += h;
      minus[j] -= h;
      double vp[2];
      double vm[2];
      MT_CHECK(measures[k](st, t, plus, vp, NULL) == MYB_OK);
      MT_CHECK(measures[k](st, t, minus, vm, NULL) == MYB_OK);
      for (size_t i = 0; i < counts[k]; ++i) {
        const double numeric = (vp[i] - vm[i]) / (2.0 * h);
        const double* row = partial + 6 * i;
        double scale = 0.0;
        for (int q = j < 3 ? 0 : 3; q < (j < 3 ? 3 : 6); ++q) {
          scale = fmax(scale, fabs(row[q]));
        }
        if (scale > 0.0) {
          worst = fmax(worst, fabs(row[j] - numeric) / scale);
        } else {
          MT_CHECK(numeric == 0.0 && row[j] == 0.0);
        }
      }
    }
    printf(" %s %.1e", names[k], worst);
    MT_CHECK(worst < 1e-7);
  }
  printf("\n");
  /* what the angles are: the object straight above the station has elevation 90 deg, to its
   * north azimuth 0 */
  double up[3];
  double north[3];
  myb_rotation inverse;
  myb_rotation_invert(&o.rotation, &inverse);
  myb_rotation_apply_position(&inverse, st->enu + 6, up);
  myb_rotation_apply_position(&inverse, st->enu + 3, north);
  double above[6] = {0};
  for (int i = 0; i < 3; ++i) {
    above[i] = site[i] + 1.0e5 * up[i] + 1.0e3 * north[i];
  }
  double angles[2];
  MT_CHECK(myb_measure_azel(st, t, above, angles, NULL) == MYB_OK);
  MT_CHECK(fabs(angles[1] - (1.5707963267948966 - 0.01)) < 1e-5);
  MT_CHECK(angles[0] < 1e-9 || angles[0] > two_pi - 1e-9);
  double range[1];
  MT_CHECK(myb_measure_range(st, t, above, range, NULL) == MYB_OK);
  MT_CHECK(fabs(range[0] - sqrt(1.0e10 + 1.0e6)) < 1e-6);
  /* errors */
  MT_CHECK(myb_measure_range(st, t, site, range, NULL) == MYB_E_DOMAIN);
  double straight[6] = {0};
  for (int i = 0; i < 3; ++i) {
    straight[i] = site[i] + 1.0e5 * up[i];
  }
  MT_CHECK(myb_measure_azel(st, t, straight, angles, NULL) == MYB_E_DOMAIN ||
           myb_measure_azel(st, t, straight, angles, NULL) == MYB_OK); /* the zenith: rounding */
  double pos[3];
  double state[6];
  double ident[36];
  MT_CHECK(myb_measure_position(NULL, t, y, pos, ident) == MYB_OK);
  MT_CHECK(pos[2] == y[2] && ident[0] == 1.0 && ident[7] == 1.0 && ident[14] == 1.0 &&
           ident[1] == 0.0 && ident[3] == 0.0);
  MT_CHECK(myb_measure_state(NULL, t, y, state, ident) == MYB_OK);
  MT_CHECK(state[5] == y[5] && ident[35] == 1.0 && ident[34] == 0.0);
}

/* ---- Observations of the truth ------------------------------------------------------------------ */

enum { MAX_OBS = 2000 };

/* Observations every `step` seconds from t0 towards t_end, the kinds and the stations taken
 * in turn, with noise of `noise` sigmas from rng (NULL: none). */
static size_t observe(scenario* s, double t0, double t_end, double step, myb_rng* rng,
                      double noise, myb_observation* obs) {
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(s->model, NULL, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, t0, s->y0, t_end) == MYB_OK);
  const myb_measure_function measures[4] = {myb_measure_range, myb_measure_range_rate,
                                            myb_measure_radec, myb_measure_azel};
  const size_t counts[4] = {1, 1, 2, 2};
  const double sigmas[4] = {1.0, 1e-3, 2e-6, 2e-6};
  const unsigned periodic[4] = {0u, 0u, 1u, 1u};
  size_t n = 0;
  const double direction = t_end > t0 ? 1.0 : -1.0;
  for (int k = 1; n < MAX_OBS; ++k) {
    const double t = t0 + direction * step * k;
    if ((t_end - t) * direction < 0.0) {
      break;
    }
    double y[6];
    MT_CHECK(myb_orbit_state(orbit, t, y) == MYB_OK);
    const int kind = k % 4;
    myb_observation* o = &obs[n++];
    memset(o, 0, sizeof *o);
    o->t = t;
    o->measure = measures[kind];
    o->context = &s->station[k % 3];
    o->count = counts[kind];
    o->periodic = periodic[kind];
    MT_CHECK(o->measure(o->context, t, y, o->value, NULL) == MYB_OK);
    for (size_t i = 0; i < o->count; ++i) {
      o->sigma[i] = sigmas[kind];
      if (rng != NULL) {
        o->value[i] += noise * sigmas[kind] * myb_rng_normal(rng);
      }
    }
  }
  myb_orbit_destroy(orbit);
  return n;
}

static void displaced(const scenario* s, double y[6]) {
  const double shift[6] = {150.0, -220.0, 90.0, 0.15, 0.08, -0.21};
  for (int i = 0; i < 6; ++i) {
    y[i] = s->y0[i] + shift[i];
  }
}

static double distance(const double a[6], const double b[6], double* dv) {
  *dv = hypot(hypot(a[3] - b[3], a[4] - b[4]), a[5] - b[5]);
  return hypot(hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2]);
}

/* ---- Least squares ------------------------------------------------------------------------------ */

static void test_fit(scenario* s) {
  static myb_observation obs[MAX_OBS];
  static double residuals[2 * MAX_OBS];
  const double arc = 3.0 * s->period;
  const double cd_true = s->object.cd;

  /* without noise: back to the truth */
  size_t n = observe(s, 0.0, arc, 60.0, NULL, 0.0, obs);
  MT_CHECK(n > 200);
  double y[6];
  displaced(s, y);
  s->object.cd = 1.7;
  myb_fit_options options = {.max_iterations = 8, .tolerance = 1e-12};
  myb_fit_result result;
  double covariance[49];
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, obs, n, &options, NULL, &result, covariance,
                         residuals) == MYB_OK);
  double dv = 0.0;
  const double dx = distance(y, s->y0, &dv);
  printf("least squares without noise: %zu observations over 3 revolutions, weighted RMS "
         "%.2e -> %.2e in %d corrections; state %.1e m, %.1e m/s, cd %.1e from the truth\n",
         n, result.rms_before, result.rms, result.iterations, dx, dv,
         fabs(s->object.cd - cd_true));
  MT_CHECK(result.rms_before > 100.0);
  MT_CHECK(result.rms < 1e-5);
  MT_CHECK(dx < 1e-5 && dv < 1e-8);
  MT_CHECK(fabs(s->object.cd - cd_true) < 1e-7);
  MT_CHECK(result.values > n);
  double worst = 0.0;
  for (size_t k = 0; k < result.values; ++k) {
    worst = fmax(worst, fabs(residuals[k]));
  }
  MT_CHECK(worst < 1e-4);

  /* with noise: within the covariance, weighted RMS 1 */
  myb_rng rng;
  myb_rng_seed(&rng, 20261110);
  s->object.cd = cd_true;
  n = observe(s, 0.0, arc, 60.0, &rng, 1.0, obs);
  displaced(s, y);
  s->object.cd = 1.7;
  options = (myb_fit_options){0};
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, obs, n, &options, NULL, &result, covariance, NULL) ==
           MYB_OK);
  MT_CHECK(result.converged);
  double largest = 0.0;
  for (int i = 0; i < 6; ++i) {
    largest = fmax(largest, fabs(y[i] - s->y0[i]) / sqrt(covariance[7 * i + i]));
  }
  const double cd_sigma = sqrt(covariance[48]);
  largest = fmax(largest, fabs(s->object.cd - cd_true) / cd_sigma);
  printf("least squares with noise: weighted RMS %.3f in %d corrections; position sigma %.2f m, "
         "cd %.4f +- %.4f; largest error %.1f sigma\n",
         result.rms, result.iterations, sqrt(covariance[0] + covariance[8] + covariance[16]),
         s->object.cd, cd_sigma, largest);
  MT_CHECK(result.rms > 0.9 && result.rms < 1.1);
  MT_CHECK(largest < 4.0);
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < i; ++j) {
      MT_CHECK(fabs(covariance[7 * i + j] - covariance[7 * j + i]) <=
               1e-12 * sqrt(covariance[8 * i] * covariance[8 * j]));
    }
  }

  /* a fixed unknown keeps its value; its rows of the covariance are zero */
  displaced(s, y);
  s->object.cd = 2.5;
  const unsigned char fixed[7] = {0, 0, 0, 0, 0, 0, 1};
  options = (myb_fit_options){.fixed = fixed};
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, obs, n, &options, NULL, &result, covariance, NULL) ==
           MYB_OK);
  MT_CHECK(s->object.cd == 2.5);
  MT_CHECK(covariance[48] == 0.0 && covariance[6] == 0.0 && covariance[42] == 0.0);
  MT_CHECK(result.rms > 1.05); /* the wrong coefficient shows in the residuals */
  /* a priori information holds an unknown near its first value */
  displaced(s, y);
  s->object.cd = 2.5;
  const double prior[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-4};
  options = (myb_fit_options){.prior_sigma = prior};
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, obs, n, &options, NULL, &result, covariance, NULL) ==
           MYB_OK);
  MT_CHECK(fabs(s->object.cd - 2.5) < 3e-4 && s->object.cd < 2.5);
  MT_CHECK(sqrt(covariance[48]) < 1e-4);
  s->object.cd = cd_true;

  /* the position of the object itself as the observations, on an arc backward in time */
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(s->model, NULL, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, s->y0, -arc) == MYB_OK);
  n = 0;
  for (int k = 0; k <= 150; ++k) {
    myb_observation* o = &obs[n++];
    memset(o, 0, sizeof *o);
    o->t = -arc * k / 150.0;
    o->measure = myb_measure_position;
    o->count = 3;
    double state[6];
    MT_CHECK(myb_orbit_state(orbit, o->t, state) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      o->value[i] = state[i];
      o->sigma[i] = 1.0;
    }
  }
  myb_orbit_destroy(orbit);
  displaced(s, y);
  s->object.cd = 1.9;
  options = (myb_fit_options){.max_iterations = 8, .tolerance = 1e-12};
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, obs, n, &options, NULL, &result, NULL, NULL) ==
           MYB_OK);
  const double back = distance(y, s->y0, &dv);
  printf("least squares on positions, backward: state %.1e m, cd %.1e from the truth\n", back,
         fabs(s->object.cd - cd_true));
  MT_CHECK(back < 1e-5 && fabs(s->object.cd - cd_true) < 1e-7);
  s->object.cd = cd_true;

  /* errors */
  memcpy(y, s->y0, sizeof y);
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, obs, 0, NULL, NULL, &result, NULL, NULL) ==
           MYB_E_DOMAIN);
  myb_observation swapped[3] = {obs[1], obs[3], obs[2]};
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, swapped, 3, NULL, NULL, &result, NULL, NULL) ==
           MYB_E_DOMAIN); /* out of order */
  MT_CHECK(myb_fit_orbit(s->model, -2.0 * arc, y, obs, n, NULL, NULL, &result, NULL, NULL) ==
           MYB_E_DOMAIN); /* ordered towards the epoch instead of away from it */
  myb_observation bad = obs[1];
  bad.sigma[1] = 0.0;
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, &bad, 1, NULL, NULL, &result, NULL, NULL) ==
           MYB_E_DOMAIN);
  bad = obs[1];
  bad.count = 7;
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, &bad, 1, NULL, NULL, &result, NULL, NULL) ==
           MYB_E_DOMAIN);
  /* one position does not determine seven unknowns */
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, &obs[1], 1, NULL, NULL, &result, NULL, NULL) ==
           MYB_E_DOMAIN);
  /* nothing to estimate: the residuals only */
  const unsigned char all_fixed[7] = {1, 1, 1, 1, 1, 1, 1};
  options = (myb_fit_options){.fixed = all_fixed};
  memcpy(y, s->y0, sizeof y);
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, obs, n, &options, NULL, &result, NULL, NULL) ==
           MYB_OK);
  MT_CHECK(result.converged && result.iterations == 0 && result.rms < 1e-5);
}

/* ---- The filter ---------------------------------------------------------------------------------- */

static void test_filter(scenario* s) {
  static myb_observation obs[MAX_OBS];
  const double arc = 2.0 * s->period;
  const double cd_true = s->object.cd;
  myb_rng rng;
  myb_rng_seed(&rng, 7);
  const size_t n = observe(s, 0.0, arc, 30.0, &rng, 1.0, obs);
  /* the truth at the end */
  double truth[6];
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(s->model, NULL, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, s->y0, arc) == MYB_OK);
  MT_CHECK(myb_orbit_state(orbit, obs[n - 1].t, truth) == MYB_OK);
  myb_orbit_destroy(orbit);

  double y[6];
  displaced(s, y);
  s->object.cd = 1.9;
  double p0[49] = {0.0};
  for (int i = 0; i < 3; ++i) {
    p0[8 * i] = 1.0e6;         /* (1 km)^2 */
    p0[8 * (3 + i)] = 1.0;     /* (1 m/s)^2 */
  }
  p0[48] = 0.25;
  myb_filter* filter = NULL;
  MT_CHECK(myb_filter_create(s->model, NULL, NULL, &filter) == MYB_OK);
  MT_CHECK(myb_filter_start(filter, 0.0, y, p0, arc) == MYB_OK);
  double first_innovation = 0.0;
  double last_innovation = 0.0;
  for (size_t k = 0; k < n; ++k) {
    double innovation[MYB_MEASURE_MAX];
    MT_CHECK(myb_filter_update(filter, &obs[k], innovation) == MYB_OK);
    const double scaled = fabs(innovation[0]) / obs[k].sigma[0];
    if (k == 0) {
      first_innovation = scaled;
    }
    if (k + 20 >= n) {
      last_innovation = fmax(last_innovation, scaled);
    }
  }
  double t = 0.0;
  double covariance[49];
  myb_filter_get(filter, &t, y, covariance);
  MT_CHECK(t == obs[n - 1].t);
  double largest = 0.0;
  for (int i = 0; i < 6; ++i) {
    largest = fmax(largest, fabs(y[i] - truth[i]) / sqrt(covariance[8 * i]));
  }
  largest = fmax(largest, fabs(s->object.cd - cd_true) / sqrt(covariance[48]));
  double dv = 0.0;
  const double dx = distance(y, truth, &dv);
  printf("filter: %zu observations over 2 revolutions; first innovation %.0f sigma, the last 20 "
         "below %.1f; end %.2f m, %.1e m/s from the truth (sigma %.2f m), cd %.4f +- %.4f; "
         "largest error %.1f sigma\n",
         n, first_innovation, last_innovation, dx, dv,
         sqrt(covariance[0] + covariance[8] + covariance[16]), s->object.cd, sqrt(covariance[48]),
         largest);
  MT_CHECK(first_innovation > 50.0);
  MT_CHECK(last_innovation < 5.0);
  MT_CHECK(dx < 5.0);
  MT_CHECK(largest < 4.0);
  for (int i = 0; i < 7; ++i) {
    MT_CHECK(covariance[8 * i] > 0.0);
    for (int j = 0; j < i; ++j) {
      MT_CHECK(covariance[7 * i + j] == covariance[7 * j + i]);
    }
  }
  /* errors of the calls */
  MT_CHECK(myb_filter_update(filter, &obs[0], NULL) == MYB_E_DOMAIN); /* behind the filter */
  MT_CHECK(myb_filter_predict(filter, arc + 1.0) == MYB_E_DOMAIN);
  MT_CHECK(myb_filter_predict(filter, arc) == MYB_OK);
  myb_filter_get(filter, &t, NULL, NULL);
  MT_CHECK(t == arc);
  MT_CHECK(myb_filter_predict(filter, arc) == MYB_OK); /* the end reached: nothing to do */
  myb_filter_destroy(filter);
  s->object.cd = cd_true;

  /* prediction alone: the covariance grows, faster with process noise */
  const double p_small[49] = {[0] = 1.0, [8] = 1.0, [16] = 1.0, [24] = 1e-6, [32] = 1e-6,
                              [40] = 1e-6, [48] = 1e-6};
  double quiet[49];
  double noisy[49];
  MT_CHECK(myb_filter_create(s->model, NULL, NULL, &filter) == MYB_OK);
  MT_CHECK(myb_filter_start(filter, 0.0, s->y0, p_small, 600.0) == MYB_OK);
  MT_CHECK(myb_filter_predict(filter, 600.0) == MYB_OK);
  myb_filter_get(filter, NULL, y, quiet);
  myb_filter_destroy(filter);
  const myb_filter_options with_noise = {.process_noise = 1e-8};
  MT_CHECK(myb_filter_create(s->model, &with_noise, NULL, &filter) == MYB_OK);
  MT_CHECK(myb_filter_start(filter, 0.0, s->y0, p_small, 600.0) == MYB_OK);
  MT_CHECK(myb_filter_predict(filter, 300.0) == MYB_OK);
  MT_CHECK(myb_filter_predict(filter, 600.0) == MYB_OK);
  myb_filter_get(filter, NULL, y, noisy);
  const double trace_quiet = quiet[0] + quiet[8] + quiet[16];
  const double trace_noisy = noisy[0] + noisy[8] + noisy[16];
  /* q dt^3 / 3 per axis for a free particle: 3 * 1e-8 * 600^3 / 3 = 2.16 m^2 (gravity bends
   * it a little) */
  printf("filter: position variance after 600 s %.3f m^2, with process noise %.3f m^2 more "
         "(q dt^3: 2.16)\n",
         trace_quiet, trace_noisy - trace_quiet);
  MT_CHECK(trace_quiet > 3.0);
  MT_CHECK(fabs((trace_noisy - trace_quiet) / 2.16 - 1.0) < 0.1);
  MT_CHECK(myb_filter_start(filter, 0.0, s->y0, p_small, -10.0) == MYB_E_DOMAIN);
  myb_filter_destroy(filter);
  myb_filter_destroy(NULL);
  const myb_filter_options negative = {.process_noise = -1.0};
  MT_CHECK(myb_filter_create(s->model, &negative, NULL, &filter) == MYB_E_DOMAIN);
}

/* Systematic shifts of the measurements as unknowns: the ranges of the three stations each
 * carry a delay, the range rates of the first a common offset. */
static void test_biases(scenario* s) {
  static myb_observation obs[MAX_OBS];
  const double arc = 3.0 * s->period;
  const double cd_true = s->object.cd;
  const double truth[4] = {35.0, -12.0, 4.5, 2.0e-3};
  myb_rng rng;
  myb_rng_seed(&rng, 31);
  const size_t n = observe(s, 0.0, arc, 40.0, &rng, 1.0, obs);
  size_t shifted = 0;
  for (size_t k = 0; k < n; ++k) {
    int station = 0;
    for (int j = 0; j < 3; ++j) {
      if (obs[k].context == &s->station[j]) {
        station = j;
      }
    }
    if (obs[k].measure == myb_measure_range) {
      obs[k].bias = 1 + station;
    } else if (obs[k].measure == myb_measure_range_rate && station == 0) {
      obs[k].bias = 4;
    }
    if (obs[k].bias > 0) {
      obs[k].value[0] += truth[obs[k].bias - 1];
      ++shifted;
    }
  }
  MT_CHECK(shifted > 100);
  double y[6];
  displaced(s, y);
  s->object.cd = 1.8;
  double bias[4] = {0.0, 0.0, 0.0, 0.0};
  myb_fit_options options = {.bias_count = 4, .bias = bias};
  myb_fit_result result;
  static double covariance[11 * 11];
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, obs, n, &options, NULL, &result, covariance, NULL) ==
           MYB_OK);
  double largest = 0.0;
  for (int k = 0; k < 4; ++k) {
    largest = fmax(largest, fabs(bias[k] - truth[k]) / sqrt(covariance[12 * (7 + k)]));
  }
  for (int i = 0; i < 6; ++i) {
    largest = fmax(largest, fabs(y[i] - s->y0[i]) / sqrt(covariance[12 * i]));
  }
  printf("least squares with shifts of the measurements: range delays %.2f, %.2f, %.2f m "
         "(sigma %.2f), a range-rate offset %.2e m/s; weighted RMS %.3f; largest error %.1f "
         "sigma\n",
         bias[0], bias[1], bias[2], sqrt(covariance[12 * 7]), bias[3], result.rms, largest);
  MT_CHECK(result.converged && result.rms > 0.9 && result.rms < 1.1);
  MT_CHECK(largest < 4.0);
  MT_CHECK(fabs(bias[0] - 35.0) < 2.0 && fabs(bias[3] - 2.0e-3) < 5.0e-4);
  /* without the shifts as unknowns the residuals show them */
  displaced(s, y);
  s->object.cd = 1.8;
  myb_fit_options plain = {0};
  const myb_status ignored =
      myb_fit_orbit(s->model, 0.0, y, obs, n, &plain, NULL, &result, NULL, NULL);
  MT_CHECK(ignored == MYB_E_DOMAIN); /* an observation names a shift the fit does not have */
  /* a shift held at its value, another kept near its first value by a priori information */
  displaced(s, y);
  s->object.cd = 1.8;
  double held[4] = {35.0, 0.0, 0.0, 0.0};
  const unsigned char fixed[11] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
  const double prior[11] = {[8] = 1.0e-3};
  options = (myb_fit_options){.bias_count = 4, .bias = held, .fixed = fixed,
                              .prior_sigma = prior};
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y, obs, n, &options, NULL, &result, covariance, NULL) ==
           MYB_OK);
  MT_CHECK(held[0] == 35.0 && covariance[12 * 7] == 0.0);
  MT_CHECK(fabs(held[1]) < 0.01); /* the truth is -12 */
  MT_CHECK(fabs(held[2] - 4.5) < 2.0);
  s->object.cd = cd_true;

  /* the same shifts through the filter */
  const myb_filter_options filter_options = {.bias_count = 4};
  myb_filter* filter = NULL;
  MT_CHECK(myb_filter_create(s->model, &filter_options, NULL, &filter) == MYB_OK);
  static double p0[11 * 11];
  memset(p0, 0, sizeof p0);
  for (int i = 0; i < 3; ++i) {
    p0[12 * i] = 1.0e6;
    p0[12 * (3 + i)] = 1.0;
    p0[12 * (7 + i)] = 1.0e4;
  }
  p0[12 * 6] = 0.25;
  p0[12 * 10] = 1.0e-4;
  displaced(s, y);
  s->object.cd = 1.8;
  MT_CHECK(myb_filter_start(filter, 0.0, y, p0, arc) == MYB_OK);
  for (size_t k = 0; k < n; ++k) {
    MT_CHECK(myb_filter_update(filter, &obs[k], NULL) == MYB_OK);
  }
  double filtered[4];
  myb_filter_bias(filter, filtered);
  myb_filter_get(filter, NULL, NULL, covariance);
  largest = 0.0;
  for (int k = 0; k < 4; ++k) {
    largest = fmax(largest, fabs(filtered[k] - truth[k]) / sqrt(covariance[12 * (7 + k)]));
  }
  printf("filter with the same shifts: %.2f, %.2f, %.2f m, %.2e m/s; largest error %.1f sigma\n",
         filtered[0], filtered[1], filtered[2], filtered[3], largest);
  MT_CHECK(largest < 4.0 && fabs(filtered[0] - 35.0) < 3.0);
  myb_filter_destroy(filter);
  s->object.cd = cd_true;
}

/* The smoother against the fit over the arc: with the same observations and the same a priori
 * information, and no process noise, the smoothed estimate at the first instant is the
 * least-squares orbit. */
static void test_smoother(scenario* s) {
  static myb_observation obs[MAX_OBS];
  const double arc = 1.5 * s->period;
  const double cd_true = s->object.cd;
  myb_rng rng;
  myb_rng_seed(&rng, 99);
  const size_t n = observe(s, 0.0, arc, 60.0, &rng, 1.0, obs);
  const double sigma0[7] = {300.0, 300.0, 300.0, 0.3, 0.3, 0.3, 0.3};
  /* the fit, started close so that it is the linear problem the filter solves */
  double y_fit[6];
  double first[6];
  for (int i = 0; i < 6; ++i) {
    first[i] = s->y0[i] + (i < 3 ? 20.0 : 0.02) * (i % 2 ? 1.0 : -1.0);
  }
  memcpy(y_fit, first, sizeof first);
  s->object.cd = 2.25;
  myb_fit_options options = {.prior_sigma = sigma0, .max_iterations = 8, .tolerance = 1e-10};
  myb_fit_result result;
  double cov_fit[49];
  MT_CHECK(myb_fit_orbit(s->model, 0.0, y_fit, obs, n, &options, NULL, &result, cov_fit,
                         NULL) == MYB_OK);
  const double cd_fit = s->object.cd;

  const myb_filter_options filter_options = {.history = n + 1};
  myb_filter* filter = NULL;
  MT_CHECK(myb_filter_create(s->model, &filter_options, NULL, &filter) == MYB_OK);
  double p0[49] = {0.0};
  for (int i = 0; i < 7; ++i) {
    p0[8 * i] = sigma0[i] * sigma0[i];
  }
  s->object.cd = 2.25;
  MT_CHECK(myb_filter_start(filter, 0.0, first, p0, arc) == MYB_OK);
  MT_CHECK(myb_filter_steps(filter) == 1);
  for (size_t k = 0; k < n; ++k) {
    MT_CHECK(myb_filter_update(filter, &obs[k], NULL) == MYB_OK);
  }
  MT_CHECK(myb_filter_steps(filter) == n + 1);
  double filtered_first[7];
  double p_first[49];
  myb_filter_step(filter, 0, NULL, filtered_first, p_first);
  MT_CHECK(filtered_first[0] == first[0]); /* before the pass: what the filter knew then */
  MT_CHECK(myb_filter_smooth(filter) == MYB_OK);
  MT_CHECK(myb_filter_smooth(filter) == MYB_OK); /* again: nothing to do */
  double t0 = -1.0;
  double x[7];
  double cov[49];
  myb_filter_step(filter, 0, &t0, x, cov);
  MT_CHECK(t0 == 0.0);
  double worst = 0.0;
  double worst_sigma = 0.0;
  for (int i = 0; i < 7; ++i) {
    const double estimate = i < 6 ? y_fit[i] : cd_fit;
    worst = fmax(worst, fabs(x[i] - estimate) / sqrt(cov_fit[8 * i]));
    worst_sigma = fmax(worst_sigma, fabs(sqrt(cov[8 * i]) / sqrt(cov_fit[8 * i]) - 1.0));
  }
  printf("smoother: the first of %zu instants against the fit over the arc — estimates within "
         "%.3f sigma, sigmas within %.1e; position sigma %.1f m before the pass, %.2f m "
         "after\n",
         n + 1, worst, worst_sigma, sqrt(p_first[0] + p_first[8] + p_first[16]),
         sqrt(cov[0] + cov[8] + cov[16]));
  MT_CHECK(worst < 0.05);
  MT_CHECK(worst_sigma < 0.02);
  /* the last instant is the filter's own estimate; in between the smoothed covariance is
   * smaller than the filter's */
  double y_end[6];
  double cov_end[49];
  double x_end[7];
  double smoothed_end[49];
  myb_filter_get(filter, NULL, y_end, cov_end);
  myb_filter_step(filter, n, NULL, x_end, smoothed_end);
  MT_CHECK(memcmp(y_end, x_end, sizeof y_end) == 0 &&
           memcmp(cov_end, smoothed_end, sizeof cov_end) == 0);
  /* against the truth in the middle of the arc */
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(s->model, NULL, NULL, &orbit) == MYB_OK);
  s->object.cd = cd_true;
  MT_CHECK(myb_orbit_start(orbit, 0.0, s->y0, arc) == MYB_OK);
  double middle[6];
  double t_middle = 0.0;
  myb_filter_step(filter, n / 2, &t_middle, x, cov);
  MT_CHECK(myb_orbit_state(orbit, t_middle, middle) == MYB_OK);
  myb_orbit_destroy(orbit);
  double largest = 0.0;
  for (int i = 0; i < 6; ++i) {
    largest = fmax(largest, fabs(x[i] - middle[i]) / sqrt(cov[8 * i]));
  }
  MT_CHECK(largest < 4.0);
  /* the record is full: one more instant does not fit */
  MT_CHECK(myb_filter_predict(filter, arc) == (obs[n - 1].t == arc ? MYB_OK : MYB_E_OVERFLOW));
  myb_filter_destroy(filter);
  /* a filter without history has nothing to smooth */
  MT_CHECK(myb_filter_create(s->model, NULL, NULL, &filter) == MYB_OK);
  MT_CHECK(myb_filter_start(filter, 0.0, first, p0, arc) == MYB_OK);
  MT_CHECK(myb_filter_smooth(filter) == MYB_E_DOMAIN && myb_filter_steps(filter) == 0);
  myb_filter_destroy(filter);
  s->object.cd = cd_true;
}

/* A point of the object away from its centre of mass. */
static void test_offset(scenario* s) {
  myb_orbital_attitude orbital = {.env = s->env, .body = s->planet};
  /* 2 m along the radius, 1 m along the normal of the orbit */
  myb_offset_measure reflector = {.measure = myb_measure_state, .attitude = myb_attitude_orbital,
                                  .attitude_context = &orbital, .offset = {2.0, 0.0, 1.0}};
  const double* y = s->y0;
  double point[6];
  double partials[36];
  MT_CHECK(myb_measure_offset(&reflector, 50.0, y, point, partials) == MYB_OK);
  const double r = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
  const double h[3] = {y[1] * y[5] - y[2] * y[4], y[2] * y[3] - y[0] * y[5],
                       y[0] * y[4] - y[1] * y[3]};
  const double hn = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
  double worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(point[i] - (y[i] + 2.0 * y[i] / r + 1.0 * h[i] / hn)));
  }
  /* the offset turns with the orbit: h / r^2 about the normal, on the 2 m along the radius */
  const double rate = hn / (r * r);
  const double turning[3] = {(h[1] * y[2] - h[2] * y[1]) / (hn * r),
                             (h[2] * y[0] - h[0] * y[2]) / (hn * r),
                             (h[0] * y[1] - h[1] * y[0]) / (hn * r)};
  double worst_v = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst_v = fmax(worst_v, fabs(point[3 + i] - (y[3 + i] + 2.0 * rate * turning[i])));
  }
  printf("a point 2 m up and 1 m across from the centre of mass: position within %.1e m, "
         "velocity %.1e m/s (it moves %.2e m/s against the centre)\n",
         worst, worst_v, 2.0 * rate);
  MT_CHECK(worst < 1e-8);
  MT_CHECK(worst_v < 1e-3 * 2.0 * rate);
  MT_CHECK(partials[0] == 1.0 && partials[7] == 1.0 && partials[1] == 0.0);
  /* the range of a station to that point is the range to where it is */
  myb_offset_measure ranged = reflector;
  ranged.measure = myb_measure_range;
  ranged.context = &s->station[0];
  double to_point = 0.0;
  double direct = 0.0;
  MT_CHECK(myb_measure_offset(&ranged, 50.0, y, &to_point, NULL) == MYB_OK);
  MT_CHECK(myb_measure_range(&s->station[0], 50.0, point, &direct, NULL) == MYB_OK);
  MT_CHECK(to_point == direct);
  double centre = 0.0;
  MT_CHECK(myb_measure_range(&s->station[0], 50.0, y, &centre, NULL) == MYB_OK);
  MT_CHECK(fabs(to_point - centre) > 0.1 && fabs(to_point - centre) < 2.3);
}

int main(void) {
  static scenario s;
  scenario_init(&s);
  test_measurements(&s);
  test_fit(&s);
  test_filter(&s);
  test_biases(&s);
  test_smoother(&s);
  test_offset(&s);
  scenario_destroy(&s);
  return mt_finish("test_estim");
}
