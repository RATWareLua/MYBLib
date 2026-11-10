/* SPDX-License-Identifier: BSD-3-Clause */
/* A satellite in a low orbit about an oblate planet with an atmosphere.
 *
 * The force model is assembled from ready-made terms (myb/env.h) — the field of the planet
 * with its J2 and the drag of a layer of gas — and propagated for a day (myb/orbit.h). The
 * program finds the ascending nodes of the day and prints the drift of the node and the
 * shortening of the period against the first-order theories.
 *
 * The library knows no body: every number below is the example's own (close to the Earth). */
#include <math.h>
#include <stdio.h>

#include "myb/myb.h"

static const double pi = 3.141592653589793;

/* Zero when the object crosses the equator; grows with time at the ascending node. */
static double equator(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[2];
}

static int fail(const char* what, myb_status status) {
  fprintf(stderr, "%s: %s\n", what, myb_status_name(status));
  return 1;
}

int main(void) {
  /* the planet */
  const double gm = 3.986004418e14; /* m^3/s^2 */
  const double radius = 6378137.0;  /* m */
  const double j2 = 1.08263e-3;
  /* the fully normalised coefficient: C20 = -J2 / sqrt(5); indices n (n + 1) / 2 + m */
  double c[6] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double s[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  c[3] = -j2 / sqrt(5.0);
  myb_gravity_field* field = NULL;
  myb_gravity_work* work = NULL;
  myb_status status = myb_gravity_field_create(gm, radius, 2, c, s, NULL, &field);
  if (status == MYB_OK) {
    status = myb_gravity_work_create(2, NULL, &work);
  }
  if (status != MYB_OK) {
    return fail("field", status);
  }

  /* the environment: one body at the origin; its axes are the inertial ones (a field of
   * zonal terms does not notice the rotation) */
  myb_env_options options = {0};
  options.epoch = myb_time_from_seconds(0.0);
  myb_env* env = NULL;
  status = myb_env_create(&options, NULL, &env);
  myb_env_body planet_body = {0};
  planet_body.id = 1;
  planet_body.gm = gm;
  size_t planet = 0;
  if (status == MYB_OK) {
    status = myb_env_add_body(env, &planet_body, &planet);
  }
  if (status != MYB_OK) {
    return fail("environment", status);
  }

  /* the terms and the object */
  myb_field_term gravity = {0};
  gravity.env = env;
  gravity.body = planet;
  gravity.field = field;
  gravity.work = work;
  gravity.n_min = 0;
  gravity.n_max = -1;
  myb_exponential_atmosphere air = {0};
  air.ellipsoid.a = radius;
  air.count = 1;
  air.base[0] = 400000.0;   /* m */
  air.density[0] = 3.0e-12; /* kg/m^3 at 400 km */
  air.scale[0] = 58000.0;   /* m */
  myb_object object = {0};
  object.mass = 500.0; /* kg */
  object.area = 4.0;   /* m^2 */
  object.cr = 1.3;
  object.cd = 2.2;
  myb_drag_term drag = {0};
  drag.env = env;
  drag.body = planet;
  drag.density = myb_density_exponential;
  drag.density_context = &air;
  drag.object = &object;
  myb_force_model* model = NULL;
  status = myb_force_model_create(NULL, &model);
  if (status == MYB_OK) {
    status = myb_force_model_add(model, "gravity", myb_term_field, NULL, &gravity, NULL);
  }
  if (status == MYB_OK) {
    status = myb_force_model_add(model, "drag", myb_term_drag, NULL, &drag, NULL);
  }
  if (status != MYB_OK) {
    return fail("model", status);
  }

  /* the orbit: 400 km, nearly circular, inclined by 51.6 degrees */
  const double a0 = radius + 400000.0;
  myb_elements elements = {0};
  elements.e = 0.001;
  elements.p = a0 * (1.0 - elements.e * elements.e);
  elements.i = 51.6 * pi / 180.0;
  elements.raan = 0.5;
  elements.argp = 1.0;
  elements.nu = 0.0;
  myb_pv pv;
  status = myb_pv_from_elements(gm, &elements, &pv);
  if (status != MYB_OK) {
    return fail("elements", status);
  }
  const double y0[6] = {pv.r[0], pv.r[1], pv.r[2], pv.v[0], pv.v[1], pv.v[2]};
  const double day = 86400.0;
  myb_orbit* orbit = NULL;
  status = myb_orbit_create(model, NULL, NULL, &orbit);
  if (status == MYB_OK) {
    status = myb_orbit_start(orbit, 0.0, y0, day);
  }
  if (status != MYB_OK) {
    return fail("propagator", status);
  }

  /* the ascending nodes of the day */
  enum { MAX_NODES = 32 };
  double nodes[MAX_NODES];
  int count = 0;
  for (int k = 0; k < MAX_NODES; ++k) {
    int found = 0;
    double t_node = 0.0;
    status = myb_orbit_event(orbit, equator, NULL, day, 1, &found, &t_node, NULL);
    if (status != MYB_OK) {
      return fail("event", status);
    }
    if (!found) {
      break;
    }
    nodes[count++] = t_node;
  }
  if (count < 3) {
    fprintf(stderr, "only %d ascending nodes in a day\n", count);
    return 1;
  }
  const double first_period = nodes[1] - nodes[0];
  const double last_period = nodes[count - 1] - nodes[count - 2];
  printf("%d ascending nodes; the first at %.3f, %.3f, %.3f s\n", count, nodes[0], nodes[1],
         nodes[2]);

  /* the state a day later (the search above stopped there) */
  double y[6];
  status = myb_orbit_state(orbit, day, y);
  if (status != MYB_OK) {
    return fail("propagation", status);
  }
  for (int i = 0; i < 3; ++i) {
    pv.r[i] = y[i];
    pv.v[i] = y[3 + i];
  }
  myb_elements after;
  status = myb_elements_from_pv(gm, &pv, &after);
  if (status != MYB_OK) {
    return fail("elements", status);
  }
  const double n = sqrt(gm / (a0 * a0 * a0));
  const double theory = -1.5 * n * j2 * (radius / elements.p) * (radius / elements.p) *
                        cos(elements.i); /* rad/s */
  const double drift = (after.raan - elements.raan) / day;
  /* drag: da = -2 pi cd A/m rho a^2 per revolution, and dP / P = 3/2 da / a */
  const double decay = 2.0 * pi * object.cd * object.area / object.mass * air.density[0] * a0 *
                       a0 * (count - 2);
  printf("after one day, %lld steps:\n", (long long)myb_orbit_steps(orbit));
  printf("  node drift %.4f deg/day (first order in J2: %.4f)\n", drift * day * 180.0 / pi,
         theory * day * 180.0 / pi);
  printf("  nodal period %.3f s, shorter by %.3f s at the end (drag, first order: %.3f)\n",
         first_period, first_period - last_period, 1.5 * decay / a0 * first_period);

  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);
  myb_env_destroy(env);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
  /* the osculating node oscillates about its mean drift by half a percent of a day's drift */
  return fabs(drift / theory - 1.0) < 0.02 && last_period < first_period ? 0 : 1;
}
