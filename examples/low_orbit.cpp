// SPDX-License-Identifier: BSD-3-Clause
// The satellite of low_orbit.c through the C++ face of the library (myb/myb.hpp): a low orbit
// about an oblate planet with an atmosphere, the ascending nodes of a day, the drift of the
// node against the first-order theory.
//
// The library knows no body: every number below is the example's own (close to the Earth).
#include <cmath>
#include <cstdio>
#include <vector>

#include "myb/myb.hpp"

namespace {

const double pi = 3.141592653589793;

int fail(const char* what, myb::status code, const std::string& message = {}) {
  std::fprintf(stderr, "%s: %s %s\n", what, myb::status_name(code), message.c_str());
  return 1;
}

} // namespace

int main() {
  const double gm = 3.986004418e14; // m^3/s^2
  const double radius = 6378137.0;  // m
  const double j2 = 1.08263e-3;
  const double day = 86400.0;

  // the field: the fully normalised C20 = -J2 / sqrt(5) at the index n (n + 1) / 2 + m
  auto field = myb::share(myb::gravity_field::create(
      gm, radius, 2, {1.0, 0.0, 0.0, -j2 / std::sqrt(5.0), 0.0, 0.0}, std::vector<double>(6, 0.0)));
  // the environment: one body at the origin
  auto env = myb::share(myb::environment::create());
  auto model = myb::share(myb::force_model::create());
  if (!field || !env || !model) {
    return fail("setup", !field ? field.code() : !env ? env.code() : model.code());
  }
  ::myb_env_body planet_body{};
  planet_body.id = 1;
  planet_body.gm = gm;
  const auto planet = (*env)->add_body(planet_body);
  if (!planet) {
    return fail("environment", planet.code());
  }

  // the terms: the model keeps the field, the environment and the structures it is given
  ::myb_exponential_atmosphere air{};
  air.ellipsoid.a = radius;
  air.count = 1;
  air.base[0] = 400000.0;   // m
  air.density[0] = 3.0e-12; // kg/m^3 at 400 km
  air.scale[0] = 58000.0;   // m
  ::myb_object object{};
  object.mass = 500.0; // kg
  object.area = 4.0;   // m^2
  object.cr = 1.3;
  object.cd = 2.2;
  const auto gravity = (*model)->add_field("gravity", *env, *planet, *field);
  const auto drag = (*model)->add_drag("drag", *env, *planet, air, object);
  if (!gravity || !drag) {
    return fail("model", !gravity ? gravity.code() : drag.code());
  }

  // the orbit: 400 km, nearly circular, inclined by 51.6 degrees
  const double a0 = radius + 400000.0;
  ::myb_elements elements{};
  elements.e = 0.001;
  elements.p = a0 * (1.0 - elements.e * elements.e);
  elements.i = 51.6 * pi / 180.0;
  elements.raan = 0.5;
  elements.argp = 1.0;
  const auto y0 = myb::state_from_elements(gm, elements);
  auto orbit = myb::orbit::create(*model);
  if (!y0 || !orbit) {
    return fail("propagator", !y0 ? y0.code() : orbit.code());
  }
  if (const auto started = orbit->start(0.0, *y0, day); !started) {
    return fail("start", started.code());
  }

  // the ascending nodes of the day: the equator crossed upwards
  std::vector<double> nodes;
  for (int k = 0; k < 32; ++k) { // bound: the nodes of a day
    const auto node = orbit->next_event([](double, const myb::state& y) { return y[2]; }, day, 1);
    if (!node) {
      return fail("event", node.code());
    }
    if (!node->has_value()) {
      break;
    }
    nodes.push_back((*node)->t);
  }
  const auto last = orbit->at(day);
  if (!last || nodes.size() < 3) {
    return fail("propagation", last.code());
  }
  const auto after = myb::elements_from_state(gm, *last);
  if (!after) {
    return fail("elements", after.code());
  }
  const double n = std::sqrt(gm / (a0 * a0 * a0));
  const double theory =
      -1.5 * n * j2 * (radius / elements.p) * (radius / elements.p) * std::cos(elements.i);
  const double drift = (after->raan - elements.raan) / day;
  const double first_period = nodes[1] - nodes[0];
  const double last_period = nodes[nodes.size() - 1] - nodes[nodes.size() - 2];
  std::printf("%zu ascending nodes, %lld steps\n", nodes.size(),
              static_cast<long long>(orbit->steps()));
  std::printf("  node drift %.4f deg/day (first order in J2: %.4f)\n", drift * day * 180.0 / pi,
              theory * day * 180.0 / pi);
  std::printf("  nodal period %.3f s, shorter by %.3f s at the end\n", first_period,
              first_period - last_period);
  return std::fabs(drift / theory - 1.0) < 0.02 && last_period < first_period ? 0 : 1;
}
