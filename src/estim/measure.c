/* SPDX-License-Identifier: BSD-3-Clause */
/* Ready-made measurements of myb/estim.h: the object itself, and range, range rate and
 * directions from a station on a body of an environment. */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "env/env_internal.h"
#include "myb/cmath.h"
#include "myb/estim.h"
#include "num/vec3.h"

static const double two_pi = 6.283185307179586;

myb_status myb_measure_position(void* context, double t, const double y[6], double* values,
                                double* partials) {
  MYB_CHECK(y != NULL && values != NULL);
  (void)context;
  (void)t;
  for (int i = 0; i < 3; ++i) {
    values[i] = y[i];
  }
  if (partials != NULL) {
    memset(partials, 0, 18 * sizeof(double));
    partials[0] = 1.0;
    partials[7] = 1.0;
    partials[14] = 1.0;
  }
  return MYB_OK;
}

myb_status myb_measure_state(void* context, double t, const double y[6], double* values,
                             double* partials) {
  MYB_CHECK(y != NULL && values != NULL);
  (void)context;
  (void)t;
  for (int i = 0; i < 6; ++i) {
    values[i] = y[i];
  }
  if (partials != NULL) {
    memset(partials, 0, 36 * sizeof(double));
    for (int i = 0; i < 6; ++i) {
      partials[7 * i] = 1.0;
    }
  }
  return MYB_OK;
}

myb_status myb_station_init(myb_station* station, myb_env* env, size_t body,
                            const myb_ellipsoid* ellipsoid, double longitude, double latitude,
                            double height) {
  MYB_CHECK(station != NULL && env != NULL && ellipsoid != NULL &&
            body < myb_env_body_count(env));
  const myb_status status =
      myb_geodetic_to_cartesian(ellipsoid, longitude, latitude, height, station->position);
  if (status != MYB_OK) {
    return status;
  }
  myb_rotation enu;
  myb_topocentric_enu(longitude, latitude, &enu);
  memcpy(station->enu, enu.m, sizeof station->enu);
  station->env = env;
  station->body = body;
  return MYB_OK;
}

/* State of the station and the rotation inertial axes -> body-fixed frame at t. */
static myb_status station_at(const myb_station* st, double t, double state[6],
                             const myb_orientation** orientation) {
  const double* s = NULL;
  myb_status status = myb_internal_env_state(st->env, t, st->body, &s);
  if (status != MYB_OK) {
    return status;
  }
  const myb_rotation* inverse = NULL;
  status = myb_internal_env_orientation(st->env, t, st->body, orientation, &inverse);
  if (status != MYB_OK) {
    return status;
  }
  const double fixed[6] = {st->position[0], st->position[1], st->position[2], 0.0, 0.0, 0.0};
  double inertial[6];
  myb_rotation_apply_state(inverse, fixed, inertial);
  for (int i = 0; i < 6; ++i) {
    state[i] = s[i] + inertial[i];
  }
  return MYB_OK;
}

myb_status myb_station_state(const myb_station* station, double t, double state[6]) {
  MYB_CHECK(station != NULL && station->env != NULL && state != NULL);
  const myb_orientation* o = NULL;
  return station_at(station, t, state, &o);
}

/* The vector from the station to the object, its rate and its length. */
static myb_status line_of_sight(const myb_station* st, double t, const double y[6], double d[3],
                                double w[3], double* range, const myb_orientation** o) {
  double site[6];
  const myb_status status = station_at(st, t, site, o);
  if (status != MYB_OK) {
    return status;
  }
  for (int i = 0; i < 3; ++i) {
    d[i] = y[i] - site[i];
    w[i] = y[3 + i] - site[3 + i];
  }
  *range = myb_vec3_norm(d);
  if (!(*range > 0.0) || !isfinite(*range)) {
    return MYB_E_DOMAIN;
  }
  return MYB_OK;
}

myb_status myb_measure_range(void* station, double t, const double y[6], double* values,
                             double* partials) {
  const myb_station* st = station;
  MYB_CHECK(st != NULL && st->env != NULL && y != NULL && values != NULL);
  double d[3];
  double w[3];
  double range = 0.0;
  const myb_orientation* o = NULL;
  const myb_status status = line_of_sight(st, t, y, d, w, &range, &o);
  if (status != MYB_OK) {
    return status;
  }
  values[0] = range;
  if (partials != NULL) {
    for (int i = 0; i < 3; ++i) {
      partials[i] = d[i] / range;
      partials[3 + i] = 0.0;
    }
  }
  return MYB_OK;
}

myb_status myb_measure_range_rate(void* station, double t, const double y[6], double* values,
                                  double* partials) {
  const myb_station* st = station;
  MYB_CHECK(st != NULL && st->env != NULL && y != NULL && values != NULL);
  double d[3];
  double w[3];
  double range = 0.0;
  const myb_orientation* o = NULL;
  const myb_status status = line_of_sight(st, t, y, d, w, &range, &o);
  if (status != MYB_OK) {
    return status;
  }
  const double rate = myb_vec3_dot(d, w) / range;
  values[0] = rate;
  if (partials != NULL) {
    /* d (d.w / |d|) / d d = (w - rate d / |d|) / |d|, d / d w = d / |d| */
    for (int i = 0; i < 3; ++i) {
      partials[i] = (w[i] - rate * d[i] / range) / range;
      partials[3 + i] = d[i] / range;
    }
  }
  return MYB_OK;
}

/* Longitude in [0, 2 pi) and latitude of the vector p in its own axes, with their gradients
 * with respect to p. MYB_E_DOMAIN on the pole of the axes. */
static myb_status direction_angles(const double p[3], double* longitude, double* latitude,
                                   double dlongitude[3], double dlatitude[3]) {
  const double h2 = p[0] * p[0] + p[1] * p[1];
  const double h = sqrt(h2);
  const double r2 = h2 + p[2] * p[2];
  if (!(h > 0.0)) {
    return MYB_E_DOMAIN;
  }
  double lon = myb_atan2(p[1], p[0]);
  if (lon < 0.0) {
    lon += two_pi;
  }
  *longitude = lon;
  *latitude = myb_atan2(p[2], h);
  dlongitude[0] = -p[1] / h2;
  dlongitude[1] = p[0] / h2;
  dlongitude[2] = 0.0;
  dlatitude[0] = -p[0] * p[2] / (r2 * h);
  dlatitude[1] = -p[1] * p[2] / (r2 * h);
  dlatitude[2] = h / r2;
  return MYB_OK;
}

myb_status myb_measure_radec(void* station, double t, const double y[6], double* values,
                             double* partials) {
  const myb_station* st = station;
  MYB_CHECK(st != NULL && st->env != NULL && y != NULL && values != NULL);
  double d[3];
  double w[3];
  double range = 0.0;
  const myb_orientation* o = NULL;
  myb_status status = line_of_sight(st, t, y, d, w, &range, &o);
  if (status != MYB_OK) {
    return status;
  }
  double da[3];
  double dd[3];
  status = direction_angles(d, &values[0], &values[1], da, dd);
  if (status != MYB_OK) {
    return status;
  }
  if (partials != NULL) {
    for (int i = 0; i < 3; ++i) {
      partials[i] = da[i];
      partials[3 + i] = 0.0;
      partials[6 + i] = dd[i];
      partials[9 + i] = 0.0;
    }
  }
  return MYB_OK;
}

myb_status myb_measure_azel(void* station, double t, const double y[6], double* values,
                            double* partials) {
  const myb_station* st = station;
  MYB_CHECK(st != NULL && st->env != NULL && y != NULL && values != NULL);
  double d[3];
  double w[3];
  double range = 0.0;
  const myb_orientation* o = NULL;
  myb_status status = line_of_sight(st, t, y, d, w, &range, &o);
  if (status != MYB_OK) {
    return status;
  }
  /* local = E M d with M inertial -> body-fixed and E body-fixed -> east, north, up; the
   * angles are those of (north, east, up): azimuth from north towards east */
  double fixed[3];
  myb_rotation_apply_position(&o->rotation, d, fixed);
  const double* e = st->enu;
  const double east = e[0] * fixed[0] + e[1] * fixed[1] + e[2] * fixed[2];
  const double north = e[3] * fixed[0] + e[4] * fixed[1] + e[5] * fixed[2];
  const double up = e[6] * fixed[0] + e[7] * fixed[1] + e[8] * fixed[2];
  const double p[3] = {north, east, up};
  double da[3];
  double de[3];
  status = direction_angles(p, &values[0], &values[1], da, de);
  if (status != MYB_OK) {
    return status;
  }
  if (partials != NULL) {
    /* gradients with respect to (north, east, up) -> body-fixed -> inertial */
    const double* m = o->rotation.m;
    for (int k = 0; k < 2; ++k) {
      const double* g = k == 0 ? da : de;
      double g_fixed[3];
      for (int j = 0; j < 3; ++j) {
        g_fixed[j] = g[0] * e[3 + j] + g[1] * e[j] + g[2] * e[6 + j];
      }
      for (int j = 0; j < 3; ++j) {
        partials[6 * k + j] = g_fixed[0] * m[j] + g_fixed[1] * m[3 + j] + g_fixed[2] * m[6 + j];
        partials[6 * k + 3 + j] = 0.0;
      }
    }
  }
  return MYB_OK;
}

/* The offset of the observed point in the inertial axes at t for the state (x, v). */
static myb_status offset_at(const myb_offset_measure* o, double t, const double x[3],
                            const double v[3], double out[3]) {
  double m[9];
  const myb_status status = o->attitude(o->attitude_context, t, x, v, m);
  if (status != MYB_OK) {
    return status;
  }
  for (int i = 0; i < 3; ++i) {
    out[i] = m[i] * o->offset[0] + m[3 + i] * o->offset[1] + m[6 + i] * o->offset[2];
  }
  return myb_vec3_finite(out) ? MYB_OK : MYB_E_DOMAIN;
}

myb_status myb_measure_offset(void* offset, double t, const double y[6], double* values,
                              double* partials) {
  const myb_offset_measure* o = offset;
  MYB_CHECK(o != NULL && o->measure != NULL && o->attitude != NULL && y != NULL &&
            values != NULL);
  static const double half = 0.5; /* s: the attitude a second apart gives its rate */
  double now[3];
  double before[3];
  double after[3];
  const double x_before[3] = {y[0] - half * y[3], y[1] - half * y[4], y[2] - half * y[5]};
  const double x_after[3] = {y[0] + half * y[3], y[1] + half * y[4], y[2] + half * y[5]};
  myb_status status = offset_at(o, t, y, y + 3, now);
  if (status == MYB_OK) {
    status = offset_at(o, t - half, x_before, y + 3, before);
  }
  if (status == MYB_OK) {
    status = offset_at(o, t + half, x_after, y + 3, after);
  }
  if (status != MYB_OK) {
    return status;
  }
  double point[6];
  for (int i = 0; i < 3; ++i) {
    point[i] = y[i] + now[i];
    point[3 + i] = y[3 + i] + (after[i] - before[i]) / (2.0 * half);
  }
  return o->measure(o->context, t, point, values, partials);
}
