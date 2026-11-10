// SPDX-License-Identifier: BSD-3-Clause
// The C++ face of the library (bindings/cpp/include/myb/myb.hpp):
//  - the low orbit of examples/low_orbit.c built from the classes gives the nodes and the
//    final state of the same calls made in C, bit for bit;
//  - a failed outcome carries the status (and the message of a parser) and no value;
//  - terms and events as lambdas; a term that fails stops the evaluation with its status;
//  - no exception handling: the test is compiled with exceptions switched off;
//  - a model keeps alive the environment and the field it was given;
//  - kernels and time scales, an SPK set and a body of it as a provider, messages written
//    and read in both notations.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "myb/myb.hpp"

namespace {

long checks = 0;
long failures = 0;

void report(const char* file, int line, const char* text) {
  ++failures;
  std::fprintf(stderr, "%s:%d: check failed: %s\n", file, line, text);
}

#define CHECK(condition)                    \
  do {                                      \
    ++checks;                               \
    if (!(condition)) {                     \
      report(__FILE__, __LINE__, #condition); \
    }                                       \
  } while (0)

bool same_bits(const myb::state& a, const myb::state& b) {
  return std::memcmp(a.data(), b.data(), sizeof(double) * 6) == 0;
}

const double pi = 3.141592653589793;
const double gm = 3.986004418e14;
const double radius = 6378137.0;
const double j2 = 1.08263e-3;
const double day = 86400.0;

::myb_exponential_atmosphere make_air() {
  ::myb_exponential_atmosphere air{};
  air.ellipsoid.a = radius;
  air.count = 1;
  air.base[0] = 400000.0;
  air.density[0] = 3.0e-12;
  air.scale[0] = 58000.0;
  return air;
}

::myb_object make_object() {
  ::myb_object object{};
  object.mass = 500.0;
  object.area = 4.0;
  object.cr = 1.3;
  object.cd = 2.2;
  return object;
}

myb::state initial_state() {
  ::myb_elements elements{};
  elements.e = 0.001;
  elements.p = (radius + 400000.0) * (1.0 - elements.e * elements.e);
  elements.i = 51.6 * pi / 180.0;
  elements.raan = 0.5;
  elements.argp = 1.0;
  elements.nu = 0.0;
  return myb::state_from_elements(gm, elements).value();
}

double equator(void*, double, const double* y) {
  return y[2];
}

struct run {
  std::vector<double> nodes;
  myb::state final{};
  long long steps = 0;
};

// The example, call for call in C.
run low_orbit_c() {
  run out;
  double c[6] = {1.0, 0.0, 0.0, -j2 / std::sqrt(5.0), 0.0, 0.0};
  double s[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  ::myb_gravity_field* field = nullptr;
  ::myb_gravity_work* work = nullptr;
  CHECK(::myb_gravity_field_create(gm, radius, 2, c, s, nullptr, &field) == MYB_OK);
  CHECK(::myb_gravity_work_create(2, nullptr, &work) == MYB_OK);
  ::myb_env_options options{};
  ::myb_env* env = nullptr;
  CHECK(::myb_env_create(&options, nullptr, &env) == MYB_OK);
  ::myb_env_body planet_body{};
  planet_body.id = 1;
  planet_body.gm = gm;
  std::size_t planet = 0;
  CHECK(::myb_env_add_body(env, &planet_body, &planet) == MYB_OK);
  ::myb_field_term gravity{};
  gravity.env = env;
  gravity.body = planet;
  gravity.field = field;
  gravity.work = work;
  gravity.n_max = -1;
  ::myb_exponential_atmosphere air = make_air();
  ::myb_object object = make_object();
  ::myb_drag_term drag{};
  drag.env = env;
  drag.body = planet;
  drag.density = ::myb_density_exponential;
  drag.density_context = &air;
  drag.object = &object;
  ::myb_force_model* model = nullptr;
  CHECK(::myb_force_model_create(nullptr, &model) == MYB_OK);
  CHECK(::myb_force_model_add(model, "gravity", ::myb_term_field, ::myb_term_field_jacobian,
                              &gravity, nullptr) == MYB_OK);
  CHECK(::myb_force_model_add(model, "drag", ::myb_term_drag, ::myb_term_drag_jacobian, &drag,
                              nullptr) == MYB_OK);
  const myb::state y0 = initial_state();
  ::myb_orbit* orbit = nullptr;
  CHECK(::myb_orbit_create(model, nullptr, nullptr, &orbit) == MYB_OK);
  CHECK(::myb_orbit_start(orbit, 0.0, y0.data(), day) == MYB_OK);
  for (int k = 0; k < 32; ++k) { // bound: the nodes of a day
    int found = 0;
    double t_node = 0.0;
    CHECK(::myb_orbit_event(orbit, equator, nullptr, day, 1, &found, &t_node, nullptr) == MYB_OK);
    if (!found) {
      break;
    }
    out.nodes.push_back(t_node);
  }
  CHECK(::myb_orbit_state(orbit, day, out.final.data()) == MYB_OK);
  out.steps = ::myb_orbit_steps(orbit);
  ::myb_orbit_destroy(orbit);
  ::myb_force_model_destroy(model);
  ::myb_env_destroy(env);
  ::myb_gravity_work_destroy(work);
  ::myb_gravity_field_destroy(field);
  return out;
}

// The same through the classes.
run low_orbit_cpp() {
  run out;
  auto field = myb::share(myb::gravity_field::create(
      gm, radius, 2, {1.0, 0.0, 0.0, -j2 / std::sqrt(5.0), 0.0, 0.0}, std::vector<double>(6, 0.0)));
  auto env = myb::share(myb::environment::create());
  auto model = myb::share(myb::force_model::create());
  CHECK(field.ok() && env.ok() && model.ok());
  if (!field || !env || !model) {
    return out;
  }
  ::myb_env_body planet_body{};
  planet_body.id = 1;
  planet_body.gm = gm;
  const auto planet = (*env)->add_body(planet_body);
  CHECK(planet.ok() && *planet == 0);
  CHECK((*model)->add_field("gravity", *env, *planet, *field).ok());
  CHECK((*model)->add_drag("drag", *env, *planet, make_air(), make_object()).ok());
  CHECK((*model)->term_count() == 2 && (*model)->term_name(1) == "drag");
  // the model holds what it needs: the caller's own references may go
  std::shared_ptr<myb::force_model> forces = *model;
  field = MYB_E_NOMEM;
  env = MYB_E_NOMEM;
  model = MYB_E_NOMEM;
  auto orbit = myb::orbit::create(forces);
  CHECK(orbit.ok());
  forces.reset(); // and the orbit holds the model
  CHECK(orbit->start(0.0, initial_state(), day).ok());
  for (int k = 0; k < 32; ++k) { // bound: the nodes of a day
    const auto node =
        orbit->next_event([](double, const myb::state& y) { return y[2]; }, day, 1);
    CHECK(node.ok());
    if (!node.ok() || !node->has_value()) {
      break;
    }
    out.nodes.push_back((*node)->t);
  }
  const auto last = orbit->at(day);
  CHECK(last.ok());
  out.final = last.value_or(myb::state{});
  out.steps = orbit->steps();
  return out;
}

void test_low_orbit() {
  const run c = low_orbit_c();
  const run cpp = low_orbit_cpp();
  CHECK(c.nodes.size() >= 14 && c.nodes.size() == cpp.nodes.size());
  bool nodes_same = c.nodes.size() == cpp.nodes.size();
  for (std::size_t k = 0; nodes_same && k < c.nodes.size(); ++k) {
    nodes_same = std::memcmp(&c.nodes[k], &cpp.nodes[k], sizeof(double)) == 0;
  }
  CHECK(nodes_same);
  CHECK(same_bits(c.final, cpp.final));
  CHECK(c.steps == cpp.steps && c.steps > 100);
  std::printf("low orbit through the classes: %zu nodes, %lld steps, the final state %s that "
              "of the C calls\n",
              cpp.nodes.size(), cpp.steps, same_bits(c.final, cpp.final) ? "bit for bit" : "NOT");
}

void test_outcomes() {
  ::myb_elements wrong{};
  wrong.p = -1.0;
  const auto failed = myb::state_from_elements(gm, wrong);
  CHECK(!failed.ok() && !failed && failed.code() == MYB_E_DOMAIN);
  CHECK(std::string(myb::status_name(failed.code())).find("DOMAIN") != std::string::npos);
  const myb::state fallback{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  CHECK(failed.value_or(fallback) == fallback);
  // two bodies: there and back
  const myb::state y0 = initial_state();
  const auto later = myb::kepler_propagate(gm, y0, 1234.5);
  CHECK(later.ok());
  const auto back = myb::kepler_propagate(gm, *later, -1234.5);
  CHECK(back.ok());
  double worst = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    worst = std::fmax(worst, std::fabs((*back)[i] - y0[i]));
  }
  CHECK(worst < 1e-6);
  const auto elements = myb::elements_from_state(gm, y0);
  CHECK(elements.ok() && std::fabs(elements->e - 0.001) < 1e-12);
  // a parser's message comes with the status
  auto pool = myb::kernel_pool::create();
  CHECK(pool.ok());
  const auto missing = pool->load_file(MYB_DATA_DIR "/no-such-kernel.tk");
  CHECK(!missing.ok() && missing.code() == MYB_E_IO);
  const auto garbage = myb::ephemeris_message::parse("CCSDS_OEM_VERS = 2.0\nMETA_START\n");
  CHECK(!garbage.ok() && garbage.code() == MYB_E_FORMAT && !garbage.message().empty());
  std::printf("a message that is not one: %s\n", garbage.message().c_str());
  const auto no_field = myb::gravity_field::create(gm, radius, 4, {1.0}, {0.0});
  CHECK(!no_field.ok() && no_field.code() == MYB_E_DOMAIN);
}

void test_callables() {
  auto model = myb::share(myb::force_model::create());
  CHECK(model.ok());
  if (!model) {
    return;
  }
  const myb::vec3 push{1.0e-6, -2.0e-6, 3.0e-6};
  CHECK((*model)->add_point_mass("planet", gm).ok());
  CHECK((*model)
            ->add_term("push",
                       [push](double, const myb::vec3&, const myb::vec3&, myb::vec3& a) {
                         a = push;
                         return MYB_OK;
                       })
            .ok());
  const myb::vec3 x{7.0e6, 1.0e5, -2.0e5};
  const myb::vec3 v{0.0, 7.5e3, 1.0e3};
  const auto both = (*model)->acceleration(0.0, x, v);
  (*model)->set_enabled(1, false);
  CHECK(!(*model)->enabled(1));
  const auto planet = (*model)->acceleration(0.0, x, v);
  (*model)->set_enabled(1, true);
  CHECK(both.ok() && planet.ok());
  for (std::size_t i = 0; i < 3; ++i) {
    CHECK(std::fabs((*both)[i] - (*planet)[i] - push[i]) < 1e-14); // the sum is rounded at 8 m/s^2
  }
  // an orbit with the push: the event and the transition matrix
  ::myb_orbit_options options{};
  options.variational = 1;
  auto orbit = myb::orbit::create(*model, options);
  CHECK(orbit.ok());
  const myb::state y0 = initial_state();
  CHECK(orbit->start(0.0, y0, 6000.0).ok());
  const auto crossing =
      orbit->next_event([](double, const myb::state& y) { return y[2]; }, 6000.0, -1);
  CHECK(crossing.ok() && crossing->has_value());
  if (crossing.ok() && crossing->has_value()) {
    CHECK(std::fabs((*crossing)->y[2]) < 1e-6 && (*crossing)->y[5] < 0.0);
  }
  const auto partials = orbit->partials_at(5000.0);
  CHECK(partials.ok());
  const auto same_time = orbit->at(5000.0);
  CHECK(same_time.ok() && same_bits(partials->y, *same_time));
  CHECK(std::fabs(partials->phi[0]) > 0.0 && std::isfinite(partials->phi[35]));
  const auto many = orbit->at(std::vector<double>{5100.0, 5200.0, 5300.0});
  CHECK(many.ok() && many->size() == 3);
  CHECK(orbit->kick({0.0, 1.0, 0.0}).ok());
  CHECK(orbit->time() == 5300.0 && orbit->steps() > 0);

  // a term that fails: its status comes back from the evaluation
  auto failing = myb::share(myb::force_model::create());
  CHECK(failing.ok());
  CHECK((*failing)
            ->add_term("refuses", [](double, const myb::vec3&, const myb::vec3&, myb::vec3&) {
              return MYB_E_RANGE;
            })
            .ok());
  CHECK((*failing)->acceleration(0.0, x, v).code() == MYB_E_RANGE);
  // a model that was moved works where it is now
  auto first = myb::force_model::create();
  CHECK(first.ok() && first->add_point_mass("planet", gm).ok());
  myb::force_model second = std::move(first).value();
  CHECK(second.acceleration(0.0, x, v).ok() && second.term_count() == 1);
}

void test_time_and_ephemeris() {
  auto pool = myb::kernel_pool::create();
  CHECK(pool.ok());
  for (const char* kernel : {MYB_DATA_DIR "/time/time-scales.tk",
                             MYB_DATA_DIR "/time/tdb-tt-fb1990.tk",
                             MYB_DATA_DIR "/time/utc-tai.tk"}) {
    CHECK(pool->load_file(kernel).ok());
  }
  const auto system = myb::time_system::create(*pool);
  CHECK(system.ok());
  if (!system) {
    return;
  }
  const auto tt = system->scale("TT");
  const auto tai = system->scale("TAI");
  CHECK(tt.ok() && tai.ok() && !system->scale("NO SUCH SCALE").ok());
  const ::myb_calendar label{2026, 11, 11, 12, 30, 15.25};
  const auto instant = system->from_utc(label, *tt);
  CHECK(instant.ok());
  const auto back = system->to_utc(*instant, *tt);
  CHECK(back.ok() && back->year == 2026 && back->minute == 30 &&
        std::fabs(back->second - 15.25) < 1e-9);
  const auto in_tai = system->convert(*instant, *tt, *tai);
  CHECK(in_tai.ok());
  CHECK(std::fabs(::myb_dd_to_double(::myb_time_diff(*instant, *in_tai)) - 32.184) < 1e-12);

  // an SPK set: the state through the class and through C; a body of it in an environment
  auto loaded = myb::ephemeris::create();
  CHECK(loaded.ok());
  CHECK(loaded->load_file(MYB_TEST_DATA_DIR "/spk/t02_03.bsp").ok());
  CHECK(!loaded->load_file(MYB_TEST_DATA_DIR "/spk/missing.bsp").ok());
  auto eph = std::make_shared<myb::ephemeris>(std::move(loaded).value());
  ::myb_spk_segment_info info{};
  ::myb_spk_segment(eph->get(), 0, &info);
  const ::myb_time middle = ::myb_time_from_seconds(0.5 * (info.start + info.stop));
  const auto there = eph->at(info.target, info.center, middle);
  myb::state expected{};
  CHECK(::myb_spk_state(eph->get(), info.target, info.center, middle, expected.data(),
                        nullptr) == MYB_OK);
  CHECK(there.ok() && same_bits(*there, expected));
  CHECK(!eph->at(info.target, info.center, ::myb_time_from_seconds(info.stop + 1.0e9)).ok());
  ::myb_env_options options{};
  options.epoch = middle;
  auto env = myb::share(myb::environment::create(options));
  CHECK(env.ok());
  ::myb_env_body origin{};
  origin.id = info.center;
  origin.gm = gm;
  ::myb_env_body moving{};
  moving.id = info.target;
  moving.gm = 1.0;
  moving.state = ::myb_spk_body_state;
  moving.state_context = eph->body(info.target, info.center);
  CHECK((*env)->add_body(origin).ok());
  const auto index = (*env)->add_body(moving);
  CHECK(index.ok() && (*env)->body_count() == 2 && (*env)->find(info.target).value_or(99) == *index);
  (*env)->keep(eph);
  eph.reset(); // the environment keeps the ephemeris
  const auto now = (*env)->at(0.0, *index);
  CHECK(now.ok() && same_bits(*now, expected));
}

void test_messages() {
  const ::myb_calendar epochs[2] = {{2026, 3, 1, 12, 0, 0.0}, {2026, 3, 1, 12, 1, 0.5}};
  const double states[12] = {7.0e6, 1.0, -2.5, 0.25, 7.5e3, -1.0e3,
                             7.0e6 + 15.0, 4.5e5, -6.0e4, -480.0, 7.4e3, -990.0};
  ::myb_oem_segment_data segment{};
  segment.object_name = "A & B";
  segment.object_id = "2026-001A";
  segment.center_name = "EARTH";
  segment.ref_frame = "EME2000";
  segment.time_system = "UTC";
  segment.count = 2;
  segment.epochs = epochs;
  segment.states = states;
  ::myb_oem_data data{};
  data.originator = "myb";
  data.creation_date = {2026, 3, 2, 0, 0, 0.0};
  data.segment_count = 1;
  data.segments = &segment;
  for (const ::myb_odm_notation notation : {MYB_ODM_KVN, MYB_ODM_XML}) {
    const auto text = myb::to_text(data, notation);
    CHECK(text.ok() && (text->front() == '<') == (notation == MYB_ODM_XML));
    const auto message = myb::ephemeris_message::parse(*text);
    CHECK(message.ok());
    if (!message) {
      continue;
    }
    CHECK(message->segment_count() == 1 && message->state_count(0) == 2);
    CHECK(message->meta(0, "OBJECT_NAME").value_or("") == "A & B");
    CHECK(!message->meta(0, "INTERPOLATION").has_value());
    const auto point = message->at(0, 1);
    CHECK(point.epoch.minute == 1 && point.epoch.second == 0.5);
    CHECK(std::memcmp(point.y.data(), states + 6, sizeof(double) * 6) == 0);
  }
  ::myb_oem_data empty{};
  CHECK(myb::to_text(empty).code() == MYB_E_DOMAIN);
}

} // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  test_low_orbit();
  test_outcomes();
  test_callables();
  test_time_and_ephemeris();
  test_messages();
  std::printf("test_cpp_binding: %ld checks, %ld failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
