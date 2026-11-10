/* SPDX-License-Identifier: BSD-3-Clause */
/* Ready-made event functions (myb/events.h). */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "env/env_internal.h"
#include "myb/cmath.h"
#include "myb/events.h"
#include "myb/force.h"
#include "myb/geometry.h"
#include "num/vec3.h"

static const double two_pi = 6.283185307179586;

/* The position and the velocity of the object relative to a body; 0 if the environment
 * fails. */
static int relative(myb_env* env, size_t body, double t, const double* y, double d[3],
                    double w[3]) {
  if (body == MYB_EVENT_ORIGIN) {
    memcpy(d, y, 3 * sizeof(double));
    memcpy(w, y + 3, 3 * sizeof(double));
    return 1;
  }
  MYB_CHECK(env != NULL);
  const double* s = NULL;
  if (myb_internal_env_state(env, t, body, &s) != MYB_OK) {
    return 0;
  }
  for (int i = 0; i < 3; ++i) {
    d[i] = y[i] - s[i];
    w[i] = y[3 + i] - s[3 + i];
  }
  return 1;
}

/* The position of the object in the body-fixed frame of a body; 0 if the environment fails. */
static int body_fixed(myb_env* env, size_t body, double t, const double* y, double fixed[3]) {
  MYB_CHECK(env != NULL && body != MYB_EVENT_ORIGIN);
  double d[3];
  double w[3];
  if (!relative(env, body, t, y, d, w)) {
    return 0;
  }
  const myb_orientation* o = NULL;
  if (myb_internal_env_orientation(env, t, body, &o, NULL) != MYB_OK) {
    return 0;
  }
  myb_rotation_apply_position(&o->rotation, d, fixed);
  return 1;
}

double myb_event_range(void* event, double t, const double* y) {
  const myb_range_event* e = event;
  MYB_CHECK(e != NULL && y != NULL);
  double d[3];
  double w[3];
  if (!relative(e->env, e->body, t, y, d, w)) {
    return (double)NAN;
  }
  return myb_vec3_norm(d) - e->range;
}

double myb_event_range_rate(void* event, double t, const double* y) {
  const myb_range_event* e = event;
  MYB_CHECK(e != NULL && y != NULL);
  double d[3];
  double w[3];
  if (!relative(e->env, e->body, t, y, d, w)) {
    return (double)NAN;
  }
  const double range = myb_vec3_norm(d);
  return range > 0.0 ? myb_vec3_dot(d, w) / range : (double)NAN;
}

double myb_event_plane(void* event, double t, const double* y) {
  const myb_plane_event* e = event;
  MYB_CHECK(e != NULL && y != NULL);
  double d[3];
  double w[3];
  if (!relative(e->env, e->body, t, y, d, w)) {
    return (double)NAN;
  }
  return myb_vec3_dot(e->normal, d);
}

double myb_event_altitude(void* event, double t, const double* y) {
  const myb_altitude_event* e = event;
  MYB_CHECK(e != NULL && y != NULL);
  double fixed[3];
  if (!body_fixed(e->env, e->body, t, y, fixed)) {
    return (double)NAN;
  }
  double longitude = 0.0;
  double latitude = 0.0;
  double height = 0.0;
  if (myb_geodetic_from_cartesian(&e->ellipsoid, fixed, &longitude, &latitude, &height) !=
      MYB_OK) {
    return (double)NAN;
  }
  return height - e->altitude;
}

double myb_event_elevation(void* event, double t, const double* y) {
  const myb_elevation_event* e = event;
  MYB_CHECK(e != NULL && e->station != NULL && y != NULL);
  const myb_station* st = e->station;
  double site[6];
  if (myb_station_state(st, t, site) != MYB_OK) {
    return (double)NAN;
  }
  const myb_orientation* o = NULL;
  if (myb_internal_env_orientation(st->env, t, st->body, &o, NULL) != MYB_OK) {
    return (double)NAN;
  }
  const double d[3] = {y[0] - site[0], y[1] - site[1], y[2] - site[2]};
  double fixed[3];
  myb_rotation_apply_position(&o->rotation, d, fixed);
  const double up = st->enu[6] * fixed[0] + st->enu[7] * fixed[1] + st->enu[8] * fixed[2];
  const double range = myb_vec3_norm(fixed);
  if (!(range > 0.0)) {
    return (double)NAN;
  }
  /* the angle from the horizon: atan2 keeps its accuracy at the zenith */
  const double east = st->enu[0] * fixed[0] + st->enu[1] * fixed[1] + st->enu[2] * fixed[2];
  const double north = st->enu[3] * fixed[0] + st->enu[4] * fixed[1] + st->enu[5] * fixed[2];
  return myb_atan2(up, sqrt(east * east + north * north)) - e->elevation;
}

double myb_event_eclipse(void* event, double t, const double* y) {
  const myb_eclipse_event* e = event;
  MYB_CHECK(e != NULL && e->env != NULL && y != NULL);
  const myb_env_body* source = myb_internal_env_body(e->env, e->source);
  const myb_env_body* occulter = myb_internal_env_body(e->env, e->occulter);
  const double* s = NULL;
  if (myb_internal_env_state(e->env, t, e->source, &s) != MYB_OK) {
    return (double)NAN;
  }
  const double to_source[3] = {s[0] - y[0], s[1] - y[1], s[2] - y[2]};
  if (myb_internal_env_state(e->env, t, e->occulter, &s) != MYB_OK) {
    return (double)NAN;
  }
  const double to_occulter[3] = {s[0] - y[0], s[1] - y[1], s[2] - y[2]};
  return myb_shadow_margin(to_source, source->radius, to_occulter, occulter->radius, e->inner);
}

myb_status myb_term_radiation_breaks(const myb_radiation_term* term, myb_force_model* model,
                                     myb_eclipse_event* events) {
  MYB_CHECK(term != NULL && term->env != NULL && model != NULL &&
            (events != NULL || term->occulter_count == 0) &&
            term->occulter_count <= MYB_ENV_MAX_OCCULTERS);
  for (size_t k = 0; k < term->occulter_count; ++k) {
    for (int inner = 0; inner < 2; ++inner) {
      myb_eclipse_event* e = &events[2 * k + (size_t)inner];
      *e = (myb_eclipse_event){.env = term->env,
                               .source = term->source,
                               .occulter = term->occulters[k],
                               .inner = inner};
      const myb_status status = myb_force_model_add_break(model, myb_event_eclipse, e);
      if (status != MYB_OK) {
        return status;
      }
    }
  }
  return MYB_OK;
}

double myb_event_shape(void* event, double t, const double* y) {
  const myb_shape_event* e = event;
  MYB_CHECK(e != NULL && e->shape != NULL && y != NULL);
  double fixed[3];
  if (!body_fixed(e->env, e->body, t, y, fixed)) {
    return (double)NAN;
  }
  return two_pi - myb_polyhedron_solid_angle(e->shape, fixed);
}
