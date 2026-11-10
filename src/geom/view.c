/* SPDX-License-Identifier: BSD-3-Clause */
/* What an observer sees of a body of an environment: the body where it was when the light
 * left it (myb/geometry.h). */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "env/env_internal.h"
#include "myb/geometry.h"
#include "myb/observe.h"
#include "num/vec3.h"

enum { LIGHT_ITERATIONS = 30 };

/* The state of the body at t - tau and the vector from the observer to it. */
static myb_status towards(myb_env* env, double t, double tau, size_t body, const double y[6],
                          double d[3]) {
  const double* s = NULL;
  const myb_status status = myb_internal_env_state(env, t - tau, body, &s);
  if (status != MYB_OK) {
    return status;
  }
  for (int i = 0; i < 3; ++i) {
    d[i] = s[i] - y[i];
  }
  return MYB_OK;
}

myb_status myb_view_body(myb_env* env, double t, size_t body, const double y[6], double c,
                         int aberration, double position[3], myb_orientation* orientation,
                         double* light_time) {
  MYB_CHECK(env != NULL && y != NULL && body < myb_env_body_count(env));
  if (!(c > 0.0) || !isfinite(c) || !isfinite(t)) {
    return MYB_E_DOMAIN;
  }
  double tau = 0.0;
  double d[3] = {0.0, 0.0, 0.0};
  /* bound: LIGHT_ITERATIONS; each gains a factor c / v */
  for (int k = 0; k < LIGHT_ITERATIONS; ++k) {
    const myb_status status = towards(env, t, tau, body, y, d);
    if (status != MYB_OK) {
      return status;
    }
    const double next = myb_vec3_norm(d) / c;
    if (!isfinite(next)) {
      return MYB_E_DOMAIN;
    }
    const int settled = fabs(next - tau) <= 1e-16 * next;
    tau = next;
    if (settled) {
      break;
    }
  }
  myb_status status = towards(env, t, tau, body, y, d);
  if (status != MYB_OK) {
    return status;
  }
  if (position != NULL) {
    const double distance = myb_vec3_norm(d);
    if (aberration && distance > 0.0) {
      const double unit[3] = {d[0] / distance, d[1] / distance, d[2] / distance};
      double seen[3];
      myb_aberration(unit, y + 3, c, seen);
      for (int i = 0; i < 3; ++i) {
        position[i] = distance * seen[i];
      }
    } else {
      memcpy(position, d, sizeof d);
    }
  }
  if (orientation != NULL) {
    status = myb_env_orientation(env, t - tau, body, orientation);
    if (status != MYB_OK) {
      return status;
    }
  }
  if (light_time != NULL) {
    *light_time = tau;
  }
  return MYB_OK;
}

myb_status myb_view_intercept(myb_env* env, double t, size_t body, const double radii[3],
                              const double y[6], const double direction[3], double c, int* hit,
                              double point[3], double* light_time) {
  MYB_CHECK(env != NULL && radii != NULL && y != NULL && direction != NULL && hit != NULL &&
            body < myb_env_body_count(env));
  *hit = 0;
  double tau = 0.0;
  myb_status status = myb_view_body(env, t, body, y, c, 0, NULL, NULL, &tau);
  if (status != MYB_OK) {
    return status;
  }
  double found[3] = {0.0, 0.0, 0.0};
  /* the light time is that of the point met, which depends on where the body was and how it
   * was turned then; bound: LIGHT_ITERATIONS */
  for (int k = 0; k < LIGHT_ITERATIONS; ++k) {
    double d[3];
    status = towards(env, t, tau, body, y, d);
    if (status != MYB_OK) {
      return status;
    }
    const myb_orientation* o = NULL;
    status = myb_internal_env_orientation(env, t - tau, body, &o, NULL);
    if (status != MYB_OK) {
      return status;
    }
    const double from_body[3] = {-d[0], -d[1], -d[2]};
    double observer[3];
    double look[3];
    myb_rotation_apply_position(&o->rotation, from_body, observer);
    myb_rotation_apply_position(&o->rotation, direction, look);
    double distance = 0.0;
    if (!myb_ray_ellipsoid(radii, observer, look, &distance, found)) {
      return MYB_OK;
    }
    const double next = distance / c;
    const int settled = fabs(next - tau) <= 1e-16 * next;
    tau = next;
    if (settled) {
      break;
    }
  }
  *hit = 1;
  if (point != NULL) {
    memcpy(point, found, sizeof found);
  }
  if (light_time != NULL) {
    *light_time = tau;
  }
  return MYB_OK;
}
