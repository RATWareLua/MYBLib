/* SPDX-License-Identifier: BSD-3-Clause */
/* The C reference of the test of the Go binding (bindings/go): the calls of
 * examples/low_orbit.c, one for one — the field with J2 and the drag of a layer of gas, a day,
 * the ascending nodes and the state at the end — with the bits of every number printed, so
 * that the Go test compares the orbit it builds from its types with these, bit for bit. The
 * inputs computed here (C20, p, i) are printed too: the Go side takes them as they are and
 * computes nothing.
 *
 * Output, one record a line, every double as the 16 hexadecimal digits of its bits:
 *   input <name> <bits>          the inputs of the example
 *   node <t> <y0> ... <y5>       an ascending node: its time and the state there
 *   final <y0> ... <y5>          the state a day later
 *   steps <count>                the accepted steps of the propagation */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "myb/myb.h"

static const double pi = 3.141592653589793;

static double equator(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[2];
}

static void print_bits(double x) {
  uint64_t bits;
  memcpy(&bits, &x, sizeof bits);
  printf(" %016" PRIx64, bits);
}

static void print_input(const char* name, double x) {
  printf("input %s", name);
  print_bits(x);
  printf("\n");
}

static int fail(const char* what, myb_status status) {
  fprintf(stderr, "%s: %s\n", what, myb_status_name(status));
  return 1;
}

int main(void) {
  const double gm = 3.986004418e14;
  const double radius = 6378137.0;
  const double j2 = 1.08263e-3;
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
  air.base[0] = 400000.0;
  air.density[0] = 3.0e-12;
  air.scale[0] = 58000.0;
  myb_object object = {0};
  object.mass = 500.0;
  object.area = 4.0;
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
  const double a0 = radius + 400000.0;
  myb_elements elements = {0};
  elements.e = 0.001;
  elements.p = a0 * (1.0 - elements.e * elements.e);
  elements.i = 51.6 * pi / 180.0;
  elements.raan = 0.5;
  elements.argp = 1.0;
  elements.nu = 0.0;
  print_input("c20", c[3]);
  print_input("p", elements.p);
  print_input("e", elements.e);
  print_input("i", elements.i);
  print_input("raan", elements.raan);
  print_input("argp", elements.argp);
  print_input("nu", elements.nu);
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

  enum { MAX_NODES = 32 };

  int count = 0;
  for (int k = 0; k < MAX_NODES; ++k) {
    int found = 0;
    double t_node = 0.0;
    double y_node[6];
    status = myb_orbit_event(orbit, equator, NULL, day, 1, &found, &t_node, y_node);
    if (status != MYB_OK) {
      return fail("event", status);
    }
    if (!found) {
      break;
    }
    printf("node");
    print_bits(t_node);
    for (int i = 0; i < 6; ++i) {
      print_bits(y_node[i]);
    }
    printf("\n");
    ++count;
  }
  double y[6];
  status = myb_orbit_state(orbit, day, y);
  if (status != MYB_OK) {
    return fail("propagation", status);
  }
  printf("final");
  for (int i = 0; i < 6; ++i) {
    print_bits(y[i]);
  }
  printf("\nsteps %lld\n", (long long)myb_orbit_steps(orbit));
  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);
  myb_env_destroy(env);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
  return count >= 14 ? 0 : 1;
}
