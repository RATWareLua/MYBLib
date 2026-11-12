/* SPDX-License-Identifier: BSD-3-Clause */
/* Measurements with the light time (myb/estim.h) against the light-time equation solved apart
 * (myb_light_time and myb_round_trip of myb/observe.h on the exact two-body motion of the
 * object and the state of the station from the environment), their derivatives against
 * differences, an orbit fitted to two-way ranges, and the derivatives of the measurement of
 * an offset point through an attitude that follows the orbit. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/env.h"
#include "myb/estim.h"
#include "myb/kepler.h"
#include "myb/observe.h"
#include "myb/surface.h"

static const double gm = 3.986004418e14;
static const double c_light = 299792458.0;

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

/* The object on its two-body orbit through y at t0, as a provider of states. */
typedef struct orbit {
  double t0;
  double y[6];
} orbit;

static myb_status orbit_state(void* context, myb_time t, double state[6]) {
  const orbit* o = context;
  myb_pv from;
  myb_pv to;
  memcpy(from.r, o->y, sizeof from.r);
  memcpy(from.v, o->y + 3, sizeof from.v);
  const myb_status status = myb_kepler_propagate(gm, &from, myb_time_seconds(t) - o->t0, &to);
  memcpy(state, to.r, sizeof to.r);
  memcpy(state + 3, to.v, sizeof to.v);
  return status;
}

static myb_status station_provider(void* context, myb_time t, double state[6]) {
  return myb_station_state(context, myb_time_seconds(t), state);
}

typedef struct scene {
  spinner spin;
  myb_env* env;
  myb_station station;
} scene;

static void scene_init(scene* s) {
  s->spin.rate = 7.292115e-5;
  const myb_env_options options = {.epoch = myb_time_from_seconds(0.0)};
  MT_CHECK(myb_env_create(&options, NULL, &s->env) == MYB_OK);
  size_t planet = 0;
  const myb_env_body body = {.gm = gm, .orientation = spinner_orientation,
                             .orientation_context = &s->spin};
  MT_CHECK(myb_env_add_body(s->env, &body, &planet) == MYB_OK);
  const myb_ellipsoid wgs84 = {6378137.0, 1.0 / 298.257223563};
  MT_CHECK(myb_station_init(&s->station, s->env, planet, &wgs84, 0.3, 0.7, 100.0) == MYB_OK);
}

/* The range of the equation solved apart: one way, and half the round trip. */
static double reference_range(scene* s, const orbit* o, double t, int two_way) {
  const myb_light_model model = {.c = c_light, .gamma = 1.0};
  orbit copy = *o;
  if (two_way) {
    double total = 0.0;
    MT_CHECK(myb_round_trip(&model, station_provider, &s->station, orbit_state, &copy,
                            myb_time_from_seconds(t), &total, NULL, NULL) == MYB_OK);
    return 0.5 * c_light * total;
  }
  double site[6];
  MT_CHECK(myb_station_state(&s->station, t, site) == MYB_OK);
  myb_light_path path;
  MT_CHECK(myb_light_time(&model, orbit_state, &copy, myb_time_from_seconds(t), site,
                          MYB_LIGHT_RECEPTION, &path) == MYB_OK);
  return path.range;
}

/* Derivatives of a measurement of `count` values against central differences. A difference
 * is as good as the value allows: 2e-6 of the largest derivative of its kind (place or
 * velocity) and the rounding of the value over the step are allowed, and the result is the
 * largest excess over that allowance (below 1: within it). */
static double check_partials(myb_measure_function measure, void* context, double t,
                             const double y[6], int count, unsigned periodic) {
  double values[MYB_MEASURE_MAX];
  double partials[6 * MYB_MEASURE_MAX];
  MT_CHECK(measure(context, t, y, values, partials) == MYB_OK);
  double worst = 0.0;
  for (int j = 0; j < 6; ++j) {
    const double h = j < 3 ? 20.0 : 0.02;
    double q[6];
    double plus[MYB_MEASURE_MAX];
    double minus[MYB_MEASURE_MAX];
    memcpy(q, y, sizeof q);
    q[j] = y[j] + h;
    MT_CHECK(measure(context, t, q, plus, NULL) == MYB_OK);
    q[j] = y[j] - h;
    MT_CHECK(measure(context, t, q, minus, NULL) == MYB_OK);
    for (int i = 0; i < count; ++i) {
      double difference = plus[i] - minus[i];
      if ((periodic >> i) & 1u) {
        difference = remainder(difference, 6.283185307179586);
      }
      double scale = 0.0;
      for (int k = 0; k < 3; ++k) {
        scale = fmax(scale, fabs(partials[6 * i + (j < 3 ? k : 3 + k)]));
      }
      const double allowed = 2e-6 * scale + 8.0 * 2.2e-16 * fmax(fabs(values[i]), 1.0) / h;
      const double e = fabs(partials[6 * i + j] - difference / (2.0 * h));
      worst = fmax(worst, e / allowed);
    }
  }
  return worst;
}

static void test_ranges(scene* s) {
  /* a low orbit over the station, a high one, and an object at the distance of the Moon */
  static const struct {
    const char* name;
    double y[6];
    double limit;
    double h; /* s: the step of the difference of ranges */
  } objects[3] = {
      {"800 km", {4.3e6, 3.6e6, 4.4e6, -5.1e3, 5.3e3, 0.7e3}, 2e-8, 0.02},
      {"36000 km", {3.2e7, 2.2e7, 1.6e7, -1.7e3, 2.4e3, 0.3e3}, 5e-8, 0.2},
      {"380000 km", {2.9e8, 2.2e8, 1.1e8, -0.6e3, 0.8e3, 0.1e3}, 5e-7, 1.0},
  };
  const double t = 1000.0;
  for (int k = 0; k < 3; ++k) {
    orbit o = {.t0 = t};
    memcpy(o.y, objects[k].y, sizeof o.y);
    double geometric = 0.0;
    MT_CHECK(myb_measure_range(&s->station, t, o.y, &geometric, NULL) == MYB_OK);
    for (int two_way = 0; two_way < 2; ++two_way) {
      myb_light_measure light = {.station = &s->station, .c = c_light, .gm = gm,
                                 .two_way = two_way};
      double range = 0.0;
      double rate = 0.0;
      MT_CHECK(myb_measure_light_range(&light, t, o.y, &range, NULL) == MYB_OK);
      MT_CHECK(myb_measure_light_range_rate(&light, t, o.y, &rate, NULL) == MYB_OK);
      const double expected = reference_range(s, &o, t, two_way);
      /* the rate in the time of reception: the range of the equation at times a step
       * apart, the object carried along its orbit */
      const double h = objects[k].h;
      double r[4];
      for (int i = 0; i < 4; ++i) {
        static const double at[4] = {-2.0, -1.0, 1.0, 2.0};
        orbit moved = {.t0 = t + at[i] * h};
        myb_pv from;
        myb_pv to;
        memcpy(from.r, o.y, sizeof from.r);
        memcpy(from.v, o.y + 3, sizeof from.v);
        MT_CHECK(myb_kepler_propagate(gm, &from, at[i] * h, &to) == MYB_OK);
        memcpy(moved.y, to.r, sizeof to.r);
        memcpy(moved.y + 3, to.v, sizeof to.v);
        r[i] = reference_range(s, &moved, t + at[i] * h, two_way);
      }
      const double expected_rate = (r[0] - 8.0 * r[1] + 8.0 * r[2] - r[3]) / (12.0 * h);
      const double worst_range =
          check_partials(myb_measure_light_range, &light, t, o.y, 1, 0);
      const double worst_rate =
          check_partials(myb_measure_light_range_rate, &light, t, o.y, 1, 0);
      printf("%s, %s: range %.4f m (%.4f m from the geometric one), from the equation "
             "solved apart %.1e m; rate %.6f m/s, from the difference of ranges %.1e m/s; "
             "derivatives against differences, of what a difference allows: %.2f, %.2f\n",
             objects[k].name, two_way ? "two-way" : "one-way", range, range - geometric,
             fabs(range - expected), rate, fabs(rate - expected_rate), worst_range, worst_rate);
      MT_CHECK(fabs(range - expected) < objects[k].limit);
      MT_CHECK(fabs(rate - expected_rate) < 1e-6);
      MT_CHECK(worst_range < 1.0 && worst_rate < 1.0);
      MT_CHECK(fabs(range - geometric) > 1e-3);
    }
  }
  /* without the mass at the origin the object goes straight over the light time: the same
   * to the square of the light time */
  myb_light_measure straight = {.station = &s->station, .c = c_light};
  myb_light_measure curved = straight;
  curved.gm = gm;
  double a = 0.0;
  double b = 0.0;
  MT_CHECK(myb_measure_light_range(&straight, t, objects[0].y, &a, NULL) == MYB_OK);
  MT_CHECK(myb_measure_light_range(&curved, t, objects[0].y, &b, NULL) == MYB_OK);
  printf("the low orbit taken as a straight line over the light time: %.2e m off\n",
         fabs(a - b));
  MT_CHECK(fabs(a - b) > 1e-6 && fabs(a - b) < 1e-3);
  /* refused */
  myb_light_measure wrong = curved;
  wrong.c = 0.0;
  MT_CHECK(myb_measure_light_range(&wrong, t, objects[0].y, &a, NULL) == MYB_E_DOMAIN);
  wrong = curved;
  wrong.gm = -1.0;
  MT_CHECK(myb_measure_light_range(&wrong, t, objects[0].y, &a, NULL) == MYB_E_DOMAIN);
  double site[6];
  MT_CHECK(myb_station_state(&s->station, t, site) == MYB_OK);
  MT_CHECK(myb_measure_light_range(&curved, t, site, &a, NULL) == MYB_E_DOMAIN);
  /* an end faster than light: no light time */
  wrong = curved;
  wrong.c = 1000.0;
  MT_CHECK(myb_measure_light_range(&wrong, t, objects[0].y, &a, NULL) != MYB_OK);
}

static void test_directions(scene* s) {
  const double t = 2500.0;
  orbit o = {.t0 = t, .y = {4.3e6, 3.6e6, 4.4e6, -5.1e3, 5.3e3, 0.7e3}};
  const myb_light_model model = {.c = c_light, .gamma = 1.0};
  double site[6];
  MT_CHECK(myb_station_state(&s->station, t, site) == MYB_OK);
  myb_light_path path;
  MT_CHECK(myb_light_time(&model, orbit_state, &o, myb_time_from_seconds(t), site,
                          MYB_LIGHT_RECEPTION, &path) == MYB_OK);
  const double d[3] = {path.target[0] - site[0], path.target[1] - site[1],
                       path.target[2] - site[2]};
  const double length = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  const double unit[3] = {d[0] / length, d[1] / length, d[2] / length};
  double seen[3];
  myb_aberration(unit, site + 3, c_light, seen);
  for (int aberration = 0; aberration < 2; ++aberration) {
    myb_light_measure light = {.station = &s->station, .c = c_light, .gm = gm,
                               .aberration = aberration};
    double radec[2];
    double azel[2];
    MT_CHECK(myb_measure_light_radec(&light, t, o.y, radec, NULL) == MYB_OK);
    MT_CHECK(myb_measure_light_azel(&light, t, o.y, azel, NULL) == MYB_OK);
    double ra = 0.0;
    double dec = 0.0;
    myb_direction_radec(aberration ? seen : unit, &ra, &dec);
    /* the same direction in the axes of the station by the geometric measurement of a point
     * placed along it */
    const double* p = aberration ? seen : unit;
    const double along[6] = {site[0] + 1e6 * p[0], site[1] + 1e6 * p[1], site[2] + 1e6 * p[2],
                             0.0, 0.0, 0.0};
    double expected_azel[2];
    MT_CHECK(myb_measure_azel(&s->station, t, along, expected_azel, NULL) == MYB_OK);
    double geometric[2];
    MT_CHECK(myb_measure_radec(&s->station, t, o.y, geometric, NULL) == MYB_OK);
    const double worst_radec = check_partials(myb_measure_light_radec, &light, t, o.y, 2, 1);
    const double worst_azel = check_partials(myb_measure_light_azel, &light, t, o.y, 2, 1);
    const double e = fmax(fmax(fabs(radec[0] - ra), fabs(radec[1] - dec)),
                          fmax(fabs(azel[0] - expected_azel[0]), fabs(azel[1] - expected_azel[1])));
    printf("the direction %s: %.2f\" from the geometric one; from the equation solved apart "
           "%.1e rad; derivatives of what a difference allows: %.2f, %.2f\n",
           aberration ? "the station sees" : "the light comes from",
           206264.8 * hypot((radec[0] - geometric[0]) * cos(dec), radec[1] - geometric[1]), e,
           worst_radec, worst_azel);
    MT_CHECK(e < 5e-13);
    MT_CHECK(worst_radec < 1.0 && worst_azel < 1.0);
  }
}

/* An orbit fitted to two-way ranges with the light time: from a displaced state back to the
 * one the ranges were made with. */
static void test_fit(scene* s) {
  enum { COUNT = 90 };
  static myb_observation observations[COUNT];
  myb_point_mass_term centre = {.gm = gm};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "centre", myb_term_point_mass,
                               myb_term_point_mass_jacobian, &centre, NULL) == MYB_OK);
  myb_light_measure light = {.station = &s->station, .c = c_light, .gm = gm, .two_way = 1};
  const orbit truth = {.t0 = 0.0, .y = {4.3e6, 3.6e6, 4.4e6, -5.1e3, 5.3e3, 0.7e3}};
  orbit copy = truth;
  for (int k = 0; k < COUNT; ++k) {
    const double t = 200.0 * k;
    double y[6];
    MT_CHECK(orbit_state(&copy, myb_time_from_seconds(t), y) == MYB_OK);
    myb_observation* o = &observations[k];
    *o = (myb_observation){.t = t, .measure = myb_measure_light_range, .context = &light,
                           .count = 1, .sigma = {0.01}};
    MT_CHECK(myb_measure_light_range(&light, t, y, o->value, NULL) == MYB_OK);
  }
  double y0[6];
  memcpy(y0, truth.y, sizeof y0);
  y0[0] += 300.0;
  y0[1] -= 200.0;
  y0[4] += 0.2;
  const myb_fit_options options = {.max_iterations = 8, .tolerance = 1e-12};
  myb_fit_result result;
  MT_CHECK(myb_fit_orbit(model, 0.0, y0, observations, COUNT, &options, NULL, &result, NULL,
                         NULL) == MYB_OK);
  double dx = 0.0;
  for (int i = 0; i < 3; ++i) {
    dx = fmax(dx, fabs(y0[i] - truth.y[i]));
  }
  printf("an orbit fitted to %d two-way ranges with the light time: weighted RMS %.1e -> "
         "%.1e in %d corrections, the state %.1e m from the truth\n",
         (int)COUNT, result.rms_before, result.rms, result.iterations, dx);
  MT_CHECK(result.rms < 1e-3 && dx < 1e-3);
  /* the same ranges taken as geometric ones put the orbit elsewhere: the light time is
   * metres of range */
  for (int k = 0; k < COUNT; ++k) {
    observations[k].measure = myb_measure_range;
    observations[k].context = &s->station;
  }
  memcpy(y0, truth.y, sizeof y0);
  MT_CHECK(myb_fit_orbit(model, 0.0, y0, observations, COUNT, &options, NULL, &result, NULL,
                         NULL) == MYB_OK);
  double off = 0.0;
  for (int i = 0; i < 3; ++i) {
    off = fmax(off, fabs(y0[i] - truth.y[i]));
  }
  printf("  the same ranges read as geometric: the orbit %.1f m away\n", off);
  MT_CHECK(off > 1.0);
  myb_force_model_destroy(model);
}

/* An attitude law of the caller's own: the body axis x along the place, y along the
 * angular momentum (differentiated by differences of the law). */
static myb_status radial_first(void* context, double t, const double x[3], const double v[3],
                              double m[9]) {
  (void)context;
  (void)t;
  const double reach = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  const double h[3] = {x[1] * v[2] - x[2] * v[1], x[2] * v[0] - x[0] * v[2],
                       x[0] * v[1] - x[1] * v[0]};
  const double hn = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
  if (!(reach > 0.0) || !(hn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    m[i] = x[i] / reach;
    m[3 + i] = h[i] / hn;
  }
  m[6] = m[1] * m[5] - m[2] * m[4];
  m[7] = m[2] * m[3] - m[0] * m[5];
  m[8] = m[0] * m[4] - m[1] * m[3];
  return MYB_OK;
}

/* A law that follows the velocity: x along it, z along the angular momentum. */
static myb_status along_velocity(void* context, double t, const double x[3], const double v[3],
                                 double m[9]) {
  (void)context;
  (void)t;
  const double speed = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  const double h[3] = {x[1] * v[2] - x[2] * v[1], x[2] * v[0] - x[0] * v[2],
                       x[0] * v[1] - x[1] * v[0]};
  const double hn = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
  if (!(speed > 0.0) || !(hn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    m[i] = v[i] / speed;
    m[6 + i] = h[i] / hn;
  }
  m[3] = m[7] * m[2] - m[8] * m[1];
  m[4] = m[8] * m[0] - m[6] * m[2];
  m[5] = m[6] * m[1] - m[7] * m[0];
  return MYB_OK;
}

/* The point of the arm on the Kepler orbit through y at t0, a time dt later. */
static void point_on_orbit(const myb_offset_measure* arm, double t0, const double y[6],
                           double dt, double point[3]) {
  myb_pv from;
  myb_pv to;
  memcpy(from.r, y, sizeof from.r);
  memcpy(from.v, y + 3, sizeof from.v);
  MT_CHECK(myb_kepler_propagate(gm, &from, dt, &to) == MYB_OK);
  double m[9];
  MT_CHECK(arm->attitude(arm->attitude_context, t0 + dt, to.r, to.v, m) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    point[i] = to.r[i] + m[i] * arm->offset[0] + m[3 + i] * arm->offset[1] +
               m[6 + i] * arm->offset[2];
  }
}

/* The velocity of the point of an arm whose attitude follows the velocity: against the
 * difference of the point on the true orbit; with gm the turning of the velocity is in it,
 * without gm it is not. */
static void test_offset_velocity(void) {
  const double t = 700.0;
  const double y[6] = {4.3e6, 3.6e6, 4.4e6, -5.1e3, 5.3e3, 0.7e3};
  myb_offset_measure arm = {.measure = myb_measure_state, .attitude = along_velocity,
                            .offset = {30.0, -20.0, 10.0}, .gm = gm};
  double values[6];
  MT_CHECK(myb_measure_offset(&arm, t, y, values, NULL) == MYB_OK);
  /* the velocity of the point: a difference of fourth order over 0.25 s */
  const double h = 0.25;
  double p[4][3];
  static const double at[4] = {-2.0, -1.0, 1.0, 2.0};
  for (int k = 0; k < 4; ++k) {
    point_on_orbit(&arm, t, y, at[k] * h, p[k]);
  }
  double with_gm = 0.0;
  double turning = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double v = (p[0][i] - 8.0 * p[1][i] + 8.0 * p[2][i] - p[3][i]) / (12.0 * h);
    with_gm = fmax(with_gm, fabs(values[3 + i] - v));
    turning = fmax(turning, fabs(v - y[3 + i]));
  }
  arm.gm = 0.0;
  MT_CHECK(myb_measure_offset(&arm, t, y, values, NULL) == MYB_OK);
  double without_gm = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double v = (p[0][i] - 8.0 * p[1][i] + 8.0 * p[2][i] - p[3][i]) / (12.0 * h);
    without_gm = fmax(without_gm, fabs(values[3 + i] - v));
  }
  arm.gm = gm;
  const double worst_partials = check_partials(myb_measure_offset, &arm, t, y, 6, 0);
  printf("a point 37 m from the centre of mass, the attitude along the velocity: it moves "
         "%.2e m/s against the centre; its velocity against the difference of its places on "
         "the orbit %.1e m/s with gm, %.1e m/s with the velocity kept; derivatives against "
         "differences, of what a difference allows: %.2f\n",
         turning, with_gm, without_gm, worst_partials);
  MT_CHECK(with_gm < 1e-8 && without_gm > 1e-3 && worst_partials < 1.0);
  arm.gm = -1.0;
  MT_CHECK(myb_measure_offset(&arm, t, y, values, NULL) == MYB_E_DOMAIN);
}

/* The measurement of a point off the centre of mass: its derivatives with the attitude that
 * follows the orbit. A long arm makes the part of the attitude visible. */
static void test_offset(scene* s) {
  const double t = 700.0;
  const double y[6] = {4.3e6, 3.6e6, 4.4e6, -5.1e3, 5.3e3, 0.7e3};
  myb_orbital_attitude orbital = {.env = s->env, .body = 0};
  static const struct {
    const char* name;
    myb_attitude_function law;
  } laws[2] = {{"the orbital axes", myb_attitude_orbital}, {"a law of the caller", radial_first}};
  for (int k = 0; k < 2; ++k) {
    myb_offset_measure arm = {.measure = myb_measure_state, .attitude = laws[k].law,
                              .attitude_context = k == 0 ? (void*)&orbital : NULL,
                              .offset = {30.0, -20.0, 10.0}};
    const double worst_state = check_partials(myb_measure_offset, &arm, t, y, 6, 0);
    double partials[36];
    double point[6];
    MT_CHECK(myb_measure_offset(&arm, t, y, point, partials) == MYB_OK);
    /* what the attitude adds to the unit matrix */
    double added = 0.0;
    for (int i = 0; i < 36; ++i) {
      added = fmax(added, fabs(partials[i] - (i % 7 == 0 ? 1.0 : 0.0)));
    }
    myb_light_measure light = {.station = &s->station, .c = c_light, .gm = gm, .two_way = 1};
    arm.measure = myb_measure_light_range;
    arm.context = &light;
    const double worst_range = check_partials(myb_measure_offset, &arm, t, y, 1, 0);
    arm.measure = myb_measure_azel;
    arm.context = &s->station;
    const double worst_azel = check_partials(myb_measure_offset, &arm, t, y, 2, 1);
    printf("a point 37 m from the centre of mass, %s: the attitude adds up to %.1e to the "
           "derivatives of its state; against differences, of what a difference allows: "
           "state %.2f, two-way range %.2f, azimuth and elevation %.2f\n",
           laws[k].name, added, worst_state, worst_range, worst_azel);
    MT_CHECK(added > 1e-6 && added < 1e-2);
    MT_CHECK(worst_state < 1.0 && worst_range < 1.0 && worst_azel < 1.0);
    /* without derivatives asked the values are the same */
    double again[6];
    arm.measure = myb_measure_state;
    arm.context = NULL;
    MT_CHECK(myb_measure_offset(&arm, t, y, again, NULL) == MYB_OK);
    MT_CHECK(memcmp(again, point, sizeof point) == 0);
  }
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  static scene s;
  scene_init(&s);
  test_ranges(&s);
  test_directions(&s);
  test_fit(&s);
  test_offset(&s);
  test_offset_velocity();
  myb_env_destroy(s.env);
  return mt_finish("test_light_measure");
}
