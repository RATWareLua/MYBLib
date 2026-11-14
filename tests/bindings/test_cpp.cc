// SPDX-License-Identifier: BSD-3-Clause
// The C++ face of the library (bindings/cpp/include/myb/myb.hpp):
//  - the low orbit of examples/low_orbit.c built from the classes gives the nodes and the
//    final state of the same calls made in C, bit for bit;
//  - a failed outcome carries the status (and the message of a parser) and no value;
//  - terms and events as lambdas; a term that fails stops the evaluation with its status;
//  - no exception handling: the test is compiled with exceptions switched off;
//  - a model keeps alive the environment and the field it was given;
//  - kernels and time scales, an SPK set and a body of it as a provider, messages written
//    and read in both notations;
//  - the upper layers of the rest of the interface, area by area — N bodies and systems of
//    bodies, estimation, the CR3BP, targeting and optimal control, attitude, observables, the
//    files, the array forms: in each the classes against the same C calls, bit for bit, and
//    failures that come back as results with the status of the library;
//  - the two-line elements of ext/sgp4 when the build links it (MYB_CPP_TEST_SGP4).
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#if defined(MYB_CPP_TEST_SGP4)
#include "myb_sgp4.h" // before myb.hpp: the header then wraps the extension too
#endif
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

// ---- The upper layers of the rest of the interface ---------------------------------------------

bool bits(const void* a, const void* b, std::size_t size) {
  return std::memcmp(a, b, size) == 0;
}

bool same_doubles(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && (a.empty() || bits(a.data(), b.data(), sizeof(double) * a.size()));
}

template <std::size_t N>
bool same_array(const std::array<double, N>& a, const double* b) {
  return bits(a.data(), b, sizeof(double) * N);
}

bool same_time(const ::myb_time& a, const ::myb_time& b) {
  return a.s == b.s && bits(&a.f, &b.f, sizeof(double));
}

bool same_rotation(const ::myb_rotation& a, const ::myb_rotation& b) {
  return bits(a.m, b.m, sizeof a.m) && bits(a.dm, b.dm, sizeof a.dm);
}

::myb_rotation identity_rotation() {
  ::myb_rotation r{};
  r.m[0] = 1.0;
  r.m[4] = 1.0;
  r.m[8] = 1.0;
  return r;
}

// The two-body motion of a state as a provider of the C interface.
struct kepler_motion {
  double mu;
  myb::state y0;
};

myb_status kepler_state(void* context, ::myb_time t, double y[6]) {
  const auto* motion = static_cast<const kepler_motion*>(context);
  const ::myb_pv pv0 = myb::to_pv(motion->y0);
  ::myb_pv pv{};
  const myb_status code = ::myb_kepler_propagate(motion->mu, &pv0, ::myb_time_seconds(t), &pv);
  if (code != MYB_OK) {
    return code;
  }
  for (int k = 0; k < 3; ++k) {
    y[k] = pv.r[k];
    y[k + 3] = pv.v[k];
  }
  return MYB_OK;
}

// The same as a lambda.
myb::provider kepler_provider(const myb::state& y0) {
  return myb::provider([y0](::myb_time t, myb::state& y) {
    const auto later = myb::kepler_propagate(gm, y0, ::myb_time_seconds(t));
    if (!later) {
      return later.code();
    }
    y = *later;
    return MYB_OK;
  });
}

myb::state tilted_state() {
  ::myb_elements elements{};
  elements.e = 0.002;
  elements.p = (radius + 400200.0) * (1.0 - elements.e * elements.e);
  elements.i = 51.0 * pi / 180.0;
  elements.raan = 0.5;
  elements.argp = 1.0;
  elements.nu = 0.0;
  return myb::state_from_elements(gm, elements).value();
}

// ---- N bodies

myb_status push_action(void* context, double, const ::myb_system_view*, double acceleration[3],
                       double[3]) {
  const auto* push = static_cast<const double*>(context);
  acceleration[0] = push[0];
  acceleration[1] = push[1];
  acceleration[2] = push[2];
  return MYB_OK;
}

void test_nbody() {
  const std::vector<double> gms = {1.32712440018e20, 1.26686534e17, 3.7931187e16};
  const std::vector<double> x = {0.0, 0.0, 0.0, 7.78e11, 0.0, 0.0, -1.43e12, 1.0e10, 0.0,
                                 1.0e11, 1.0e11, 0.0};
  const std::vector<double> v = {0.0, 0.0, 0.0, 0.0, 13070.0, 0.0, 0.0, -9690.0, 100.0,
                                 -2.0e4, 2.0e4, 0.0};
  ::myb_nbody* raw = nullptr;
  CHECK(::myb_nbody_create(3, 1, gms.data(), nullptr, &raw) == MYB_OK);
  std::vector<double> a_c(12);
  CHECK(::myb_nbody_acceleration(raw, 0.0, x.data(), a_c.data()) == MYB_OK);
  ::myb_nbody_invariants invariants_c{};
  CHECK(::myb_nbody_invariants_of(raw, x.data(), v.data(), &invariants_c) == MYB_OK);

  auto made = myb::nbody::create(3, 1, gms);
  CHECK(made.ok());
  if (!made) {
    ::myb_nbody_destroy(raw);
    return;
  }
  myb::nbody bodies = *made;
  CHECK(bodies.count() == 4 && bodies.massive() == 3 && bodies.gm(3) == 0.0);
  const auto a = bodies.acceleration(0.0, x);
  CHECK(a.ok() && same_doubles(*a, a_c));
  const auto invariants = bodies.invariants_of(x, v);
  CHECK(invariants.ok() && bits(&*invariants, &invariants_c, sizeof invariants_c));
  // threads give the same bits; the pool is kept by the system
  auto threads = myb::share(myb::thread_pool::create(2));
  CHECK(threads.ok() && (*threads)->threads() == 2);
  if (threads) {
    bodies.set_pool(*threads);
    threads = MYB_E_NOMEM;
    const auto a_threads = bodies.acceleration(0.0, x);
    CHECK(a_threads.ok() && same_doubles(*a_threads, a_c));
    bodies.set_pool(nullptr);
  }
  CHECK(bodies.set_approximation(8, 0.5).ok());
  CHECK(bodies.approximation().order == 8 && bodies.approximation().theta == 0.5);
  CHECK(bodies.set_approximation(0, 0.0).ok() && bodies.approximation().order == 0);
  CHECK(bodies.set_approximation(99, 0.5).code() == MYB_E_DOMAIN);
  CHECK(bodies.set_interaction(1, 1, MYB_INTERACTION_NONE).code() == MYB_E_DOMAIN);
  // the 1PN terms: only the accelerations with velocities, and the integrators refuse them
  const ::myb_ppn ppn{299792458.0, 1.0, 1.0};
  CHECK(bodies.set_relativity(ppn).ok() && ::myb_nbody_set_relativity(raw, &ppn) == MYB_OK);
  CHECK(bodies.acceleration(0.0, x).code() == MYB_E_DOMAIN);
  std::vector<double> a_pn_c(12);
  CHECK(::myb_nbody_acceleration_xv(raw, 0.0, x.data(), v.data(), a_pn_c.data()) == MYB_OK);
  const auto a_pn = bodies.acceleration_xv(0.0, x, v);
  CHECK(a_pn.ok() && same_doubles(*a_pn, a_pn_c));
  ::myb_nbody_pn_invariants pn_c{};
  CHECK(::myb_nbody_pn_invariants_of(raw, x.data(), v.data(), &pn_c) == MYB_OK);
  const auto pn = bodies.pn_invariants_of(x, v);
  CHECK(pn.ok() && bits(&*pn, &pn_c, sizeof pn_c));
  ::myb_wh_options options{};
  options.h = 8.64e5;
  options.scheme = MYB_WH_SABA2;
  auto refused = myb::wh::create(bodies, options);
  CHECK(!refused.ok() || refused->start(0.0, x, v).code() == MYB_E_DOMAIN);
  CHECK(bodies.clear_relativity().ok() && ::myb_nbody_set_relativity(raw, nullptr) == MYB_OK);
  // Wisdom-Holman: 100 steps of ten days
  ::myb_wh* wh_c = nullptr;
  CHECK(::myb_wh_create(raw, &options, nullptr, &wh_c) == MYB_OK);
  CHECK(::myb_wh_start(wh_c, 0.0, x.data(), v.data()) == MYB_OK);
  CHECK(::myb_wh_integrate(wh_c, 100.0 * options.h) == MYB_OK);
  std::vector<double> x_c(12);
  std::vector<double> v_c(12);
  CHECK(::myb_wh_state(wh_c, x_c.data(), v_c.data()) == MYB_OK);
  auto wh = myb::wh::create(bodies, options);
  CHECK(wh.ok());
  if (wh) {
    CHECK(wh->start(0.0, x, v).ok() && wh->integrate(100.0 * options.h).ok());
    const auto now = wh->current();
    CHECK(now.ok() && same_doubles(now->x, x_c) && same_doubles(now->v, v_c));
    CHECK(wh->time() == ::myb_wh_time(wh_c) && wh->statistics().steps == 100);
    CHECK(wh->start(0.0, x, std::vector<double>(3)).code() == MYB_E_DOMAIN);
    std::printf("Wisdom-Holman through the class: %lld steps, the state %s that of the C calls\n",
                static_cast<long long>(wh->statistics().steps),
                now.ok() && same_doubles(now->x, x_c) ? "bit for bit" : "NOT");
  }
  ::myb_wh_destroy(wh_c);
  ::myb_nbody_destroy(raw);
  CHECK(myb::nbody::create(3, 1, {1.0, 2.0}).code() == MYB_E_DOMAIN);
  CHECK(myb::nbody::create(2, 0, {1.0, -1.0}).code() == MYB_E_DOMAIN);
}

// Two spheres let go at rest (tests/unit/test_system.c): contact, an elastic rebound, contact
// again and the merger, through the class and in C.
struct encounter_run {
  double t_contact = 0.0;
  double t_second = 0.0;
  double energy = 0.0;
  std::vector<double> final_state;
};

encounter_run encounters_c() {
  encounter_run out;
  const double gms[2] = {4.0e5, 1.0e5};
  ::myb_system* s = nullptr;
  CHECK(::myb_system_create(2, gms, nullptr, &s) == MYB_OK);
  CHECK(::myb_system_set_radius(s, 0, 60.0) == MYB_OK);
  CHECK(::myb_system_set_radius(s, 1, 40.0) == MYB_OK);
  const double y[12] = {-1000.0, 0.0, 0.0, 4000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  ::myb_ias15_options options{};
  options.epsilon = 1e-10;
  ::myb_ias15* in = nullptr;
  CHECK(::myb_system_integrator(s, &options, nullptr, &in) == MYB_OK);
  CHECK(::myb_ias15_start(in, ::myb_system_acceleration, s, 0.0, y) == MYB_OK);
  ::myb_system_event event{};
  CHECK(::myb_system_advance(s, in, 1.0e6, 0.0, &event) == MYB_OK);
  out.t_contact = event.t;
  CHECK(::myb_system_bounce(s, in, 0, 1, 1.0) == MYB_OK);
  CHECK(::myb_system_advance(s, in, 2.0 * out.t_contact, 0.0, &event) == MYB_OK);
  CHECK(::myb_system_advance(s, in, 1.0e6, 0.0, &event) == MYB_OK);
  out.t_second = event.t;
  ::myb_system_merger merger{};
  CHECK(::myb_system_merge(s, in, 0, 1, &merger) == MYB_OK);
  out.energy = merger.energy;
  CHECK(::myb_system_advance(s, in, event.t + 100.0, 0.0, &event) == MYB_OK);
  out.final_state.assign(::myb_ias15_state(in), ::myb_ias15_state(in) + 12);
  ::myb_ias15_destroy(in);
  ::myb_system_destroy(s);
  return out;
}

encounter_run encounters_cpp() {
  encounter_run out;
  auto s = myb::system::create({4.0e5, 1.0e5});
  CHECK(s.ok());
  if (!s) {
    return out;
  }
  CHECK(s->set_radius(0, 60.0).ok() && s->set_radius(1, 40.0).ok());
  CHECK(s->radius(0) == 60.0 && s->bodies() == 2 && s->coordinates() == 6);
  const std::vector<double> y = {-1000.0, 0.0, 0.0, 4000.0, 0.0, 0.0,
                                 0.0,     0.0, 0.0, 0.0,    0.0, 0.0};
  ::myb_ias15_options options{};
  options.epsilon = 1e-10;
  auto in = s->integrator(&options);
  CHECK(in.ok() && in->coordinates() == 6);
  if (!in) {
    return out;
  }
  CHECK(s->start(*in, 0.0, y).ok());
  const auto first = s->advance(*in, 1.0e6);
  CHECK(first.ok() && first->kind == MYB_SYSTEM_CONTACT && first->a == 0 && first->b == 1);
  out.t_contact = first.ok() ? first->t : 0.0;
  CHECK(s->bounce(*in, 0, 1, 1.0).ok());
  const auto back = s->advance(*in, 2.0 * out.t_contact);
  CHECK(back.ok() && back->kind == MYB_SYSTEM_ARRIVED);
  const auto second = s->advance(*in, 1.0e6);
  CHECK(second.ok() && second->kind == MYB_SYSTEM_CONTACT);
  out.t_second = second.ok() ? second->t : 0.0;
  const auto merger = s->merge(*in, 0, 1);
  CHECK(merger.ok() && s->present(0) && !s->present(1));
  CHECK(s->nbody().gm(1) == 0.0 && s->nbody().gm(0) == 5.0e5);
  out.energy = merger.ok() ? merger->energy : 0.0;
  CHECK(s->advance(*in, out.t_second + 100.0).ok());
  out.final_state = in->current();
  // failures: the body is gone; a watch that is negative
  CHECK(s->merge(*in, 0, 1).code() == MYB_E_DOMAIN);
  CHECK(s->advance(*in, 0.0, -1.0).code() == MYB_E_DOMAIN);
  CHECK(s->start(*in, 0.0, std::vector<double>(3)).code() == MYB_E_DOMAIN);
  return out;
}

void test_system() {
  const encounter_run c = encounters_c();
  const encounter_run cpp = encounters_cpp();
  CHECK(c.t_contact > 0.0 && bits(&c.t_contact, &cpp.t_contact, sizeof(double)));
  CHECK(bits(&c.t_second, &cpp.t_second, sizeof(double)) && bits(&c.energy, &cpp.energy, sizeof(double)));
  CHECK(same_doubles(c.final_state, cpp.final_state));
  std::printf("two spheres through the class: contact at %.9f s, merger at %.9f s, %s the C calls\n",
              cpp.t_contact, cpp.t_second,
              same_doubles(c.final_state, cpp.final_state) ? "bit for bit" : "NOT");

  // a planet with a field and a frame from outside, a satellite pushed by an action
  const double c_coeff[6] = {1.0, 0.0, 0.0, -j2 / std::sqrt(5.0), 0.0, 0.0};
  const double s_coeff[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const myb::state satellite = initial_state();
  std::vector<double> y0(12, 0.0);
  for (std::size_t k = 0; k < 3; ++k) {
    y0[3 + k] = satellite[k];
    y0[9 + k] = satellite[3 + k];
  }
  double push[3] = {1.0e-6, -2.0e-6, 5.0e-7};
  ::myb_rotation identity = identity_rotation();
  ::myb_ias15_options options{};
  options.epsilon = 1e-10;
  // C
  ::myb_gravity_field* field_c = nullptr;
  CHECK(::myb_gravity_field_create(gm, radius, 2, c_coeff, s_coeff, nullptr, &field_c) == MYB_OK);
  const double gms[2] = {gm, 1.0};
  ::myb_system* s_c = nullptr;
  CHECK(::myb_system_create(2, gms, nullptr, &s_c) == MYB_OK);
  CHECK(::myb_system_set_orientation(s_c, 0, ::myb_orientation_fixed, &identity) == MYB_OK);
  CHECK(::myb_system_set_field(s_c, 0, field_c, 0) == MYB_OK);
  CHECK(::myb_system_add_field(s_c, 0, 1) == MYB_OK);
  CHECK(::myb_system_add_action(s_c, 1, push_action, push) == MYB_OK);
  std::vector<double> a_c(6);
  CHECK(::myb_system_acceleration(s_c, 0.0, y0.data(), y0.data() + 6, a_c.data()) == MYB_OK);
  ::myb_ias15* in_c = nullptr;
  CHECK(::myb_system_integrator(s_c, &options, nullptr, &in_c) == MYB_OK);
  CHECK(::myb_ias15_start(in_c, ::myb_system_acceleration, s_c, 0.0, y0.data()) == MYB_OK);
  CHECK(::myb_ias15_integrate(in_c, 3000.0) == MYB_OK);
  const std::vector<double> end_c(::myb_ias15_state(in_c), ::myb_ias15_state(in_c) + 12);
  const std::size_t which_c[1] = {1};
  std::vector<double> y1_c(12);
  std::vector<double> phi_c(144);
  std::vector<double> sensitivity_c(12);
  CHECK(::myb_system_partials(s_c, &options, 0.0, y0.data(), 300.0, which_c, 1, y1_c.data(),
                              phi_c.data(), sensitivity_c.data()) == MYB_OK);
  ::myb_ias15_destroy(in_c);
  ::myb_system_destroy(s_c);
  ::myb_gravity_field_destroy(field_c);
  // C++
  auto field = myb::share(myb::gravity_field::create(gm, radius, 2,
                                                     {1.0, 0.0, 0.0, -j2 / std::sqrt(5.0), 0.0, 0.0},
                                                     std::vector<double>(6, 0.0)));
  auto s = myb::system::create({gm, 1.0});
  CHECK(field.ok() && s.ok());
  if (!field || !s) {
    return;
  }
  CHECK(s->set_orientation(0, myb::fixed_orientation(identity_rotation())).ok());
  CHECK(s->set_field(0, *field).ok());
  field = MYB_E_NOMEM; // the system keeps the field
  CHECK(s->add_field(0, 1).ok());
  const myb::vec3 push_cpp{push[0], push[1], push[2]};
  CHECK(s->add_action(1, [push_cpp](double, const ::myb_system_view&, myb::vec3& a, myb::vec3&) {
          a = push_cpp;
          return MYB_OK;
        }).ok());
  const auto a = s->acceleration(0.0, std::vector<double>(y0.begin(), y0.begin() + 6),
                                 std::vector<double>(y0.begin() + 6, y0.end()));
  CHECK(a.ok() && same_doubles(*a, a_c));
  auto in = s->integrator(&options);
  CHECK(in.ok());
  if (in) {
    CHECK(s->start(*in, 0.0, y0).ok() && in->integrate(3000.0).ok());
    CHECK(same_doubles(in->current(), end_c) && in->time() == 3000.0);
    const auto dense = in->dense(2999.0);
    CHECK(dense.ok() && dense->size() == 12);
  }
  const auto derivatives = s->partials(0.0, y0, 300.0, {1}, &options);
  CHECK(derivatives.ok() && same_doubles(derivatives->y1, y1_c) &&
        same_doubles(derivatives->phi, phi_c) && same_doubles(derivatives->sensitivity, sensitivity_c));
  CHECK(s->pack_rotation(0, {1.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 1e-4}, std::nullopt, y0).code() ==
        MYB_E_DOMAIN); // the rotation of the body is not integrated
  CHECK(s->nbody().count() == 2 && s->set_radius(5, 1.0).code() == MYB_E_DOMAIN);
  std::printf("a planet with a field and a pushed satellite through the class: the state and the "
              "partials %s the C calls\n",
              in.ok() && same_doubles(in->current(), end_c) && derivatives.ok() &&
                      same_doubles(derivatives->phi, phi_c)
                  ? "bit for bit with"
                  : "NOT as");

  // an integrated rotation into the state and out of it
  ::myb_system_rotation rotation{};
  rotation.moments[0] = 1.0e12;
  rotation.moments[1] = 1.1e12;
  rotation.moments[2] = 1.2e12;
  ::myb_system* r_c = nullptr;
  CHECK(::myb_system_create(2, gms, nullptr, &r_c) == MYB_OK);
  CHECK(::myb_system_set_rotation(r_c, 0, &rotation) == MYB_OK);
  const double q[4] = {0.9, 0.1, -0.3, 0.2};
  const double omega[3] = {1.0e-5, -2.0e-5, 7.0e-5};
  std::vector<double> yr_c(2 * ::myb_system_coordinates(r_c), 0.0);
  yr_c[3] = 7.0e6; // the second body away from the first
  CHECK(::myb_system_pack_rotation(r_c, 0, q, omega, nullptr, yr_c.data()) == MYB_OK);
  double inertia_c[9] = {};
  CHECK(::myb_system_inertia(r_c, 0, 0.0, yr_c.data(), inertia_c) == MYB_OK);
  double q_c[4] = {};
  double omega_c[3] = {};
  double core_c[3] = {};
  CHECK(::myb_system_unpack_rotation(r_c, 0, yr_c.data(), q_c, omega_c, core_c) == MYB_OK);
  ::myb_system_destroy(r_c);
  auto r = myb::system::create({gm, 1.0});
  CHECK(r.ok() && r->set_rotation(0, rotation).ok());
  if (r) {
    std::vector<double> yr(2 * r->coordinates(), 0.0);
    CHECK(yr.size() == 30 && yr.size() == yr_c.size());
    yr[3] = 7.0e6;
    CHECK(r->pack_rotation(0, {q[0], q[1], q[2], q[3]}, {omega[0], omega[1], omega[2]},
                           std::nullopt, yr)
              .ok());
    CHECK(same_doubles(yr, yr_c));
    const auto inertia = r->inertia(0, 0.0, yr);
    CHECK(inertia.ok() && same_array(*inertia, inertia_c));
    const auto back = r->unpack_rotation(0, yr);
    CHECK(back.ok() && same_array(back->q, q_c) && same_array(back->omega, omega_c) &&
          same_array(back->omega_core, core_c));
    CHECK(r->unpack_spins(0, yr).code() == MYB_E_DOMAIN); // no pair with spins
  }
}

// ---- Estimation

myb_status point_mass_force(::myb_force_model** model, ::myb_point_mass_term* term) {
  const myb_status code = ::myb_force_model_create(nullptr, model);
  if (code != MYB_OK) {
    return code;
  }
  return ::myb_force_model_add(*model, "planet", ::myb_term_point_mass,
                               ::myb_term_point_mass_jacobian, term, nullptr);
}

void test_estimation() {
  ::myb_rotation identity = identity_rotation();
  auto env = myb::share(myb::environment::create());
  CHECK(env.ok());
  if (!env) {
    return;
  }
  ::myb_env_body planet{};
  planet.id = 399;
  planet.gm = gm;
  planet.orientation = ::myb_orientation_fixed;
  planet.orientation_context = &identity;
  CHECK((*env)->add_body(planet).ok());
  const ::myb_ellipsoid wgs84{6378137.0, 1.0 / 298.257223563};
  ::myb_station station_c{};
  CHECK(::myb_station_init(&station_c, (*env)->get(), 0, &wgs84, 0.3, 0.8, 100.0) == MYB_OK);
  auto station = myb::station::create(*env, 0, wgs84, 0.3, 0.8, 100.0);
  CHECK(station.ok());
  if (!station) {
    return;
  }
  const myb::state y0 = initial_state();
  double values_c[2] = {};
  double partials_c[12] = {};
  CHECK(::myb_measure_azel(&station_c, 10.0, y0.data(), values_c, partials_c) == MYB_OK);
  const auto azel = station->measure_azel().evaluate(10.0, y0);
  CHECK(azel.ok() && azel->values.size() == 2 && bits(azel->values.data(), values_c, sizeof values_c) &&
        bits(azel->partials.data(), partials_c, sizeof partials_c));
  double light_c[1] = {};
  ::myb_light_measure light{};
  light.station = &station_c;
  light.c = 299792458.0;
  light.gm = gm;
  CHECK(::myb_measure_light_range(&light, 10.0, y0.data(), light_c, nullptr) == MYB_OK);
  const auto light_range = station->measure_light_range(299792458.0, gm, false).evaluate(10.0, y0);
  CHECK(light_range.ok() && bits(light_range->values.data(), light_c, sizeof light_c));
  double place_c[6] = {};
  CHECK(::myb_station_state(&station_c, 10.0, place_c) == MYB_OK);
  const auto place = station->at(10.0);
  CHECK(place.ok() && same_array(*place, place_c));

  // a fit of positions and ranges about a point mass, the guess 100 m off
  std::vector<::myb_observation> raw;
  std::vector<myb::observation> observations;
  const myb::measure ranges = station->measure_range();
  for (int k = 1; k <= 5; ++k) {
    const double t = 600.0 * k;
    const auto truth = myb::kepler_propagate(gm, y0, t);
    const auto later = myb::kepler_propagate(gm, y0, t + 300.0);
    CHECK(truth.ok() && later.ok());
    ::myb_observation position{};
    position.t = t;
    position.measure = ::myb_measure_position;
    position.count = 3;
    for (std::size_t i = 0; i < 3; ++i) {
      position.value[i] = (*truth)[i];
      position.sigma[i] = 1.0;
    }
    raw.push_back(position);
    const auto o = myb::observation::create(myb::measure_position(), t,
                                            {(*truth)[0], (*truth)[1], (*truth)[2]},
                                            {1.0, 1.0, 1.0});
    CHECK(o.ok());
    observations.push_back(*o);
    double range = 0.0;
    CHECK(::myb_measure_range(&station_c, t + 300.0, later->data(), &range, nullptr) == MYB_OK);
    ::myb_observation range_c{};
    range_c.t = t + 300.0;
    range_c.measure = ::myb_measure_range;
    range_c.context = &station_c;
    range_c.count = 1;
    range_c.value[0] = range;
    range_c.sigma[0] = 0.5;
    raw.push_back(range_c);
    const auto r = myb::observation::create(ranges, t + 300.0, {range}, {0.5});
    CHECK(r.ok());
    observations.push_back(*r);
  }
  myb::state guess = y0;
  guess[0] += 100.0;
  ::myb_point_mass_term term{gm};
  ::myb_force_model* model_c = nullptr;
  CHECK(point_mass_force(&model_c, &term) == MYB_OK);
  myb::state fitted_c = guess;
  ::myb_fit_result result_c{};
  std::vector<double> covariance_c(36);
  std::vector<double> residuals_c(20);
  CHECK(::myb_fit_orbit(model_c, 0.0, fitted_c.data(), raw.data(), raw.size(), nullptr, nullptr,
                        &result_c, covariance_c.data(), residuals_c.data()) == MYB_OK);
  auto model = myb::share(myb::force_model::create());
  CHECK(model.ok() && (*model)->add_point_mass("planet", gm).ok());
  if (!model) {
    ::myb_force_model_destroy(model_c);
    return;
  }
  const auto fitted = myb::fit_orbit(**model, 0.0, guess, observations);
  CHECK(fitted.ok());
  if (fitted) {
    CHECK(bits(fitted->y0.data(), fitted_c.data(), sizeof(double) * 6));
    CHECK(same_doubles(fitted->covariance, covariance_c) && same_doubles(fitted->residuals, residuals_c));
    CHECK(fitted->summary.iterations == result_c.iterations && fitted->summary.values == 20 &&
          bits(&fitted->summary.rms, &result_c.rms, sizeof(double)));
    std::printf("a fit of 5 positions and 5 ranges through the class: %d iterations, rms %.3g, %s\n",
                fitted->summary.iterations, fitted->summary.rms,
                bits(fitted->y0.data(), fitted_c.data(), sizeof(double) * 6)
                    ? "bit for bit with the C calls"
                    : "NOT as the C calls");
  }
  CHECK(myb::fit_orbit(**model, 0.0, guess, {}).code() == MYB_E_DOMAIN);
  myb::fit_options wrong{};
  wrong.fixed.assign(3, 0);
  CHECK(myb::fit_orbit(**model, 0.0, guess, observations, wrong).code() == MYB_E_DOMAIN);
  CHECK(myb::observation::create(myb::measure_position(), 1.0, {1.0, 2.0}, {1.0, 1.0}).code() ==
        MYB_E_DOMAIN);

  // the filter and its smoother over the same observations
  std::vector<double> prior(36, 0.0);
  for (std::size_t k = 0; k < 6; ++k) {
    prior[7 * k] = k < 3 ? 1.0e4 : 1.0;
  }
  ::myb_filter_options filter_options{};
  filter_options.history = 12;
  ::myb_filter* filter_c = nullptr;
  CHECK(::myb_filter_create(model_c, &filter_options, nullptr, &filter_c) == MYB_OK);
  CHECK(::myb_filter_start(filter_c, 0.0, guess.data(), prior.data(), 4000.0) == MYB_OK);
  std::vector<double> innovation_c(3);
  for (const ::myb_observation& o : raw) {
    CHECK(::myb_filter_update(filter_c, &o, innovation_c.data()) == MYB_OK);
  }
  CHECK(::myb_filter_predict(filter_c, 3900.0) == MYB_OK);
  double t_c = 0.0;
  double y_c[6] = {};
  std::vector<double> p_c(36);
  ::myb_filter_get(filter_c, &t_c, y_c, p_c.data());
  CHECK(::myb_filter_smooth(filter_c) == MYB_OK);
  std::vector<double> first_c(6);
  std::vector<double> first_p_c(36);
  double first_t_c = -1.0;
  ::myb_filter_step(filter_c, 0, &first_t_c, first_c.data(), first_p_c.data());
  const std::size_t steps_c = ::myb_filter_steps(filter_c);
  ::myb_filter_destroy(filter_c);
  auto filter = myb::filter::create(*model, filter_options);
  CHECK(filter.ok());
  if (filter) {
    CHECK(filter->start(0.0, guess, std::vector<double>(35), 4000.0).code() == MYB_E_DOMAIN);
    CHECK(filter->start(0.0, guess, prior, 4000.0).ok());
    bool updated = true;
    std::vector<double> last_innovation;
    for (const myb::observation& o : observations) {
      const auto innovation = filter->update(o);
      updated = updated && innovation.ok();
      if (innovation.ok()) {
        last_innovation = *innovation;
      }
    }
    CHECK(updated && last_innovation.size() == 1 && bits(last_innovation.data(), innovation_c.data(), sizeof(double)));
    CHECK(filter->predict(3900.0).ok());
    const myb::filter::estimate now = filter->current();
    CHECK(now.t == t_c && same_array(now.y, y_c) && same_doubles(now.covariance, p_c));
    CHECK(filter->smooth().ok() && filter->steps() == steps_c && filter->unknown_count() == 6);
    const myb::filter::kept_estimate first = filter->step(0);
    CHECK(first.t == first_t_c && same_doubles(first.unknowns, first_c) &&
          same_doubles(first.covariance, first_p_c));
    CHECK(filter->bias().empty() && filter->set_bias({1.0}).code() == MYB_E_DOMAIN);
    CHECK(filter->predict(5000.0).code() == MYB_E_DOMAIN); // beyond t_end
  }

  // the variational equations of the model
  ::myb_variational* variational_c = nullptr;
  CHECK(::myb_variational_create(model_c, nullptr, nullptr, &variational_c) == MYB_OK);
  CHECK(::myb_variational_start(variational_c, 0.0, y0.data()) == MYB_OK);
  CHECK(::myb_variational_integrate(variational_c, 1000.0) == MYB_OK);
  double yv_c[6] = {};
  double phi_c[36] = {};
  ::myb_variational_get(variational_c, yv_c, phi_c, nullptr);
  ::myb_variational_destroy(variational_c);
  auto variational = myb::variational::create(*model);
  CHECK(variational.ok());
  if (variational) {
    CHECK(variational->start(0.0, y0).ok() && variational->integrate(1000.0).ok());
    const auto now = variational->current();
    CHECK(same_array(now.y, yv_c) && same_array(now.phi, phi_c) && now.sensitivity.empty());
    CHECK(variational->dense(5000.0).code() == MYB_E_RANGE);
  }
  ::myb_force_model_destroy(model_c);

  // covariances and random numbers
  ::myb_rng rng_c{};
  ::myb_rng_seed(&rng_c, 42u);
  myb::rng generator(42u);
  CHECK(generator.next() == ::myb_rng_next(&rng_c));
  const double normal_c = ::myb_rng_normal(&rng_c);
  const double normal = generator.normal();
  CHECK(bits(&normal, &normal_c, sizeof(double)));
  myb::mat6 p{};
  for (std::size_t k = 0; k < 6; ++k) {
    p[7 * k] = 1.0 + static_cast<double>(k);
  }
  double rtn_c[36] = {};
  CHECK(::myb_covariance_rtn(y0.data(), y0.data() + 3, p.data(), rtn_c) == MYB_OK);
  const auto rtn = myb::covariance_rtn({y0[0], y0[1], y0[2]}, {y0[3], y0[4], y0[5]}, p);
  CHECK(rtn.ok() && same_array(*rtn, rtn_c));
  CHECK(myb::covariance_rtn({1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, p).code() == MYB_E_DOMAIN);
  CHECK(myb::sample_statistics(1, 2, {1.0, 2.0}).code() == MYB_E_DOMAIN);
}

// ---- The CR3BP

void test_cr3bp() {
  const double mu = 1.215058560962404e-02;
  double points_c[5][3] = {};
  CHECK(::myb_cr3bp_lagrange(mu, points_c) == MYB_OK);
  const auto points = myb::cr3bp_lagrange(mu);
  CHECK(points.ok() && bits(points->data(), points_c, sizeof points_c));
  CHECK(myb::cr3bp_lagrange(0.7).code() == MYB_E_DOMAIN);
  CHECK(myb::cr3bp_richardson(0.0121, 3, 0.01, 1).code() == MYB_E_DOMAIN);
  // the L1 Lyapunov orbit of the JPL catalogue
  double y_c[6] = {8.2522726142834313e-01, 0.0, 0.0, 0.0, 1.0748883813645686e-01, 0.0};
  const myb::state guess{y_c[0], y_c[1], y_c[2], y_c[3], y_c[4], y_c[5]};
  ::myb_cr3bp_correction correction{};
  correction.symmetry = MYB_CR3BP_SYMMETRY_XZ;
  correction.fixed = 0;
  double period_c = 0.0;
  CHECK(::myb_cr3bp_correct(mu, &correction, y_c, &period_c) == MYB_OK);
  double mono_c[36] = {};
  CHECK(::myb_cr3bp_monodromy(mu, y_c, period_c, mono_c) == MYB_OK);
  double re_c[6] = {};
  double im_c[6] = {};
  double index_c = 0.0;
  CHECK(::myb_cr3bp_stability(mono_c, re_c, im_c, &index_c) == MYB_OK);
  const auto orbit = myb::cr3bp_correct(mu, correction, guess);
  CHECK(orbit.ok() && same_array(orbit->y0, y_c) && orbit->period == period_c);
  if (!orbit) {
    return;
  }
  const auto mono = myb::cr3bp_monodromy(mu, orbit->y0, orbit->period);
  CHECK(mono.ok() && same_array(*mono, mono_c));
  const auto floquet = myb::cr3bp_stability(*mono);
  CHECK(floquet.ok() && same_array(floquet->re, re_c) && same_array(floquet->im, im_c) &&
        floquet->index == index_c);
  double branch_c[6] = {};
  CHECK(::myb_cr3bp_manifold(mu, y_c, period_c, 0.25, 1, 1, 1e-6, branch_c) == MYB_OK);
  const auto branch = myb::cr3bp_manifold(mu, orbit->y0, orbit->period, 0.25, true, 1, 1e-6);
  CHECK(branch.ok() && same_array(*branch, branch_c));
  CHECK(myb::cr3bp_manifold(mu, orbit->y0, orbit->period, 1.5, true, 1, 1e-6).code() ==
        MYB_E_DOMAIN);
  // to the section of the Moon
  ::myb_cr3bp_section section{};
  section.axis = 0;
  section.value = 1.0 - mu;
  section.direction = 0;
  double end_c[6] = {};
  double elapsed_c = 0.0;
  int crossed_c = 0;
  CHECK(::myb_cr3bp_propagate(mu, branch_c, 10.0, &section, end_c, &elapsed_c, &crossed_c) ==
        MYB_OK);
  const auto arrival = myb::cr3bp_propagate(mu, *branch, 10.0, section);
  CHECK(arrival.ok() && same_array(arrival->y, end_c) && arrival->elapsed == elapsed_c &&
        arrival->crossed == (crossed_c != 0));
  CHECK(myb::cr3bp_jacobi(mu, orbit->y0) == ::myb_cr3bp_jacobi(mu, y_c));
  // multiple shooting and a step of the family
  ::myb_cr3bp_orbit patched_c{};
  CHECK(::myb_cr3bp_orbit_from_state(mu, y_c, period_c, 4, &patched_c) == MYB_OK);
  CHECK(::myb_cr3bp_orbit_correct(mu, &patched_c, 0.0, 0) == MYB_OK);
  CHECK(::myb_cr3bp_orbit_continue(mu, &patched_c, 1e-3, 0.0, 0) == MYB_OK);
  auto patched = myb::cr3bp_orbit_from_state(mu, orbit->y0, orbit->period, 4);
  CHECK(patched.ok());
  if (patched) {
    CHECK(myb::cr3bp_orbit_correct(mu, *patched).ok() &&
          myb::cr3bp_orbit_continue(mu, *patched, 1e-3).ok());
    CHECK(bits(patched->states, patched_c.states, sizeof patched_c.states) &&
          patched->period == patched_c.period);
  }
  CHECK(myb::cr3bp_orbit_from_state(mu, orbit->y0, orbit->period, 11).code() == MYB_E_DOMAIN);
  // a Lissajous torus about it, its spectrum and a branch of its manifold
  ::myb_cr3bp_torus torus_c{};
  CHECK(::myb_cr3bp_torus_from_orbit(mu, y_c, period_c, 25, 0, 1e-3, &torus_c) == MYB_OK);
  CHECK(::myb_cr3bp_torus_continue(mu, &torus_c, MYB_CR3BP_TORUS_JACOBI, 0.0, 0.0, 0) == MYB_OK);
  auto spectrum_c = std::make_unique<::myb_cr3bp_torus_spectrum>();
  CHECK(::myb_cr3bp_torus_stability(mu, &torus_c, 0.0, spectrum_c.get()) == MYB_OK);
  auto torus = myb::cr3bp_torus_from_orbit(mu, orbit->y0, orbit->period, 25, 0, 1e-3);
  CHECK(torus.ok());
  if (!torus) {
    return;
  }
  CHECK(myb::cr3bp_torus_continue(mu, *torus, MYB_CR3BP_TORUS_JACOBI, 0.0).ok());
  CHECK(bits(torus->states, torus_c.states, sizeof torus_c.states) &&
        torus->period == torus_c.period && torus->rotation == torus_c.rotation);
  const auto spectrum = myb::cr3bp_torus_stability(mu, *torus);
  CHECK(spectrum.ok());
  if (spectrum) {
    CHECK(spectrum->count == spectrum_c->count &&
          bits(spectrum->re, spectrum_c->re, sizeof spectrum_c->re) &&
          bits(spectrum->unstable_direction, spectrum_c->unstable_direction,
               sizeof spectrum_c->unstable_direction) &&
          spectrum->unstable == spectrum_c->unstable &&
          spectrum->normal_rotation == spectrum_c->normal_rotation);
    double leaving_c[6] = {};
    CHECK(::myb_cr3bp_torus_manifold(&torus_c, spectrum_c.get(), 0.5, 1, 1, 1e-6, leaving_c) ==
          MYB_OK);
    const auto leaving = myb::cr3bp_torus_manifold(*torus, *spectrum, 0.5, true, 1, 1e-6);
    CHECK(leaving.ok() && same_array(*leaving, leaving_c));
    std::printf("a Lissajous torus of 25 points through the class: lambda %.6f, %zu hyperbolic, "
                "%s the C calls\n",
                spectrum->unstable, spectrum->hyperbolic,
                leaving.ok() && same_array(*leaving, leaving_c) ? "bit for bit with" : "NOT as");
  }
  double point_c[6] = {};
  ::myb_cr3bp_torus_state(&torus_c, 1.0, point_c);
  CHECK(same_array(myb::cr3bp_torus_state(*torus, 1.0), point_c));
  // the term of a force model
  auto model = myb::force_model::create();
  CHECK(model.ok() && model->add_cr3bp("cr3bp", mu).ok());
  if (model) {
    const myb::vec3 x{0.8, 0.01, 0.0};
    const myb::vec3 v{0.0, 0.1, 0.0};
    const auto a = model->acceleration(0.0, x, v);
    CHECK(a.ok() && *a == myb::cr3bp_acceleration(mu, x, v));
  }
}

// ---- Targeting, optimal control, low thrust, conjunctions

myb_status square_root_of_two(void*, const double* u, double* f, double* jacobian) {
  f[0] = u[0] * u[0] - 2.0;
  if (jacobian != nullptr) {
    jacobian[0] = 2.0 * u[0];
  }
  return MYB_OK;
}

myb_status nearest_on_line(void*, const double* x, double* f, double* c, double* gradient,
                           double* jacobian) {
  *f = (x[0] - 1.0) * (x[0] - 1.0) + (x[1] - 2.0) * (x[1] - 2.0);
  c[0] = x[0] + x[1] - 1.0;
  if (gradient != nullptr) {
    gradient[0] = 2.0 * (x[0] - 1.0);
    gradient[1] = 2.0 * (x[1] - 2.0);
    jacobian[0] = 1.0;
    jacobian[1] = 1.0;
  }
  return MYB_OK;
}

myb_status constant_torque(void* context, double, const double[4], double torque[3]) {
  const auto* n = static_cast<const double*>(context);
  torque[0] = n[0];
  torque[1] = n[1];
  torque[2] = n[2];
  return MYB_OK;
}

void test_targeting() {
  const myb::state y0 = initial_state();
  const auto later = myb::kepler_propagate(gm, y0, 1800.0);
  CHECK(later.ok());
  if (!later) {
    return;
  }
  const myb::vec3 r1{y0[0], y0[1], y0[2]};
  const myb::vec3 r2{(*later)[0], (*later)[1], (*later)[2]};
  ::myb_lambert_solution solutions_c[3] = {};
  std::size_t count_c = 0;
  CHECK(::myb_lambert(gm, r1.data(), r2.data(), 1800.0, 0, 1, solutions_c, 3, &count_c) ==
        MYB_OK);
  const auto solutions = myb::lambert(gm, r1, r2, 1800.0, false, 1);
  CHECK(solutions.ok() && solutions->size() == count_c && count_c >= 1);
  if (solutions.ok() && !solutions->empty()) {
    CHECK(bits((*solutions)[0].v1, solutions_c[0].v1, sizeof solutions_c[0].v1) &&
          bits((*solutions)[0].v2, solutions_c[0].v2, sizeof solutions_c[0].v2));
  }
  CHECK(myb::lambert(gm, {0.0, 0.0, 0.0}, r2, 1800.0).code() == MYB_E_DOMAIN);
  // a first orbit from three positions
  const auto third = myb::kepler_propagate(gm, y0, 3600.0);
  CHECK(third.ok());
  const myb::vec3 r3{(*third)[0], (*third)[1], (*third)[2]};
  double v2_c[3] = {};
  double tilt_c = 0.0;
  CHECK(::myb_iod_gibbs(gm, r1.data(), r2.data(), r3.data(), v2_c, &tilt_c) == MYB_OK);
  const auto gibbs = myb::iod_gibbs(gm, r1, r2, r3);
  CHECK(gibbs.ok() && same_array(gibbs->v2, v2_c) && gibbs->out_of_plane == tilt_c);
  CHECK(myb::iod_herrick_gibbs(gm, r1, 2.0, r2, 1.0, r3, 3.0).code() == MYB_E_DOMAIN);
  // relative motion and the B-plane
  const myb::state deputy = tilted_state();
  double relative_c[6] = {};
  const ::myb_pv chief_pv = myb::to_pv(y0);
  const ::myb_pv deputy_pv = myb::to_pv(deputy);
  CHECK(::myb_relative_from_inertial(&chief_pv, &deputy_pv, relative_c) == MYB_OK);
  const auto relative = myb::relative_from_inertial(y0, deputy);
  CHECK(relative.ok() && same_array(*relative, relative_c));
  double stm_c[36] = {};
  CHECK(::myb_relative_stm(gm, &chief_pv, 600.0, stm_c) == MYB_OK);
  const auto stm = myb::relative_stm(gm, y0, 600.0);
  CHECK(stm.ok() && same_array(*stm, stm_c));
  const double hyperbola[6] = {7.0e6, 0.0, 0.0, 0.0, 12000.0, 1000.0};
  const double pole[3] = {0.0, 0.0, 1.0};
  ::myb_bplane bplane_c{};
  CHECK(::myb_bplane_compute(gm, hyperbola, hyperbola + 3, pole, &bplane_c) == MYB_OK);
  const auto bplane = myb::bplane_compute(gm, {7.0e6, 0.0, 0.0}, {0.0, 12000.0, 1000.0},
                                          {0.0, 0.0, 1.0});
  CHECK(bplane.ok() && bits(&*bplane, &bplane_c, sizeof bplane_c));
  CHECK(myb::bplane_compute(gm, r1, {y0[3], y0[4], y0[5]}, {0.0, 0.0, 1.0}).code() ==
        MYB_E_DOMAIN); // an ellipse
  // shooting and a nonlinear program, the callables as lambdas
  std::vector<double> u_c = {1.0};
  ::myb_shoot_options shoot_options{1e-14, 0};
  std::vector<double> work(::myb_shoot_work_size(1, 1));
  int iterations_c = 0;
  CHECK(::myb_shoot(square_root_of_two, nullptr, 1, 1, u_c.data(), &shoot_options, work.data(),
                    &iterations_c, nullptr) == MYB_OK);
  std::vector<double> u = {1.0};
  const auto shot = myb::shoot(
      [](const double* at, double* f, double* jacobian) {
        f[0] = at[0] * at[0] - 2.0;
        if (jacobian != nullptr) {
          jacobian[0] = 2.0 * at[0];
        }
        return MYB_OK;
      },
      1, u, shoot_options);
  CHECK(shot.ok() && same_doubles(u, u_c) && shot->iterations == iterations_c);
  ::myb_nlp_problem problem_c{};
  problem_c.n = 2;
  problem_c.equalities = 1;
  problem_c.function = nearest_on_line;
  std::vector<double> x_c = {0.0, 0.0};
  ::myb_nlp_result report_c{};
  CHECK(::myb_nlp_solve(&problem_c, nullptr, x_c.data(), nullptr, nullptr, &report_c) == MYB_OK);
  myb::nlp_problem problem;
  problem.n = 2;
  problem.equalities = 1;
  problem.function = [](const double* at, double& f, double* c, double* gradient,
                        double* jacobian) { return nearest_on_line(nullptr, at, &f, c, gradient, jacobian); };
  std::vector<double> x = {0.0, 0.0};
  const auto solved = myb::nlp_solve(problem, x);
  CHECK(solved.ok() && same_doubles(x, x_c) &&
        bits(&solved->report.objective, &report_c.objective, sizeof(double)) &&
        solved->multipliers.size() == 1);
  std::vector<double> short_x = {0.0};
  CHECK(myb::nlp_solve(problem, short_x).code() == MYB_E_DOMAIN);
  // a leg of Sims and Flanagan with the engine off, and a powered flyby
  ::myb_sims_flanagan leg{};
  leg.mu = gm;
  leg.thrust = 0.1;
  leg.exhaust = 3.0e4;
  leg.t0 = 0.0;
  leg.mass0 = 500.0;
  leg.t1 = 1800.0;
  leg.segments = 4;
  for (std::size_t k = 0; k < 6; ++k) {
    leg.start[k] = y0[k];
    leg.target[k] = (*later)[k];
  }
  const std::vector<double> throttles(12, 0.0);
  double mismatch_c[6] = {};
  double mass1_c = 0.0;
  std::vector<double> jacobian_c(72);
  CHECK(::myb_sims_flanagan_mismatch(&leg, throttles.data(), mismatch_c, &mass1_c,
                                     jacobian_c.data()) == MYB_OK);
  const auto mismatch = myb::sims_flanagan_mismatch(leg, throttles);
  CHECK(mismatch.ok() && same_array(mismatch->mismatch, mismatch_c) && mismatch->mass1 == mass1_c &&
        same_doubles(mismatch->jacobian, jacobian_c));
  CHECK(myb::sims_flanagan_mismatch(leg, std::vector<double>(3)).code() == MYB_E_DOMAIN);
  double dv_c = 0.0;
  double turn_c = 0.0;
  double flyby_partials_c[6] = {};
  CHECK(::myb_chain_flyby(3.24859e14, 3000.0, 3100.0, 7.0e6, &dv_c, &turn_c, flyby_partials_c) ==
        MYB_OK);
  const auto flyby = myb::chain_flyby(3.24859e14, 3000.0, 3100.0, 7.0e6);
  CHECK(flyby.ok() && flyby->dv == dv_c && flyby->turn == turn_c &&
        same_array(flyby->partials, flyby_partials_c));
  CHECK(myb::chain_flyby(3.24859e14, -1.0, 3100.0, 7.0e6).code() == MYB_E_DOMAIN);
  // low thrust: an arc, Edelbaum, the Q-law
  ::myb_thrust_arc arc{};
  arc.start = 0.0;
  arc.stop = 1000.0;
  arc.thrust = 0.5;
  arc.exhaust = 2.0e4;
  arc.steering = MYB_STEER_VELOCITY;
  ::myb_thrust thrust_c{};
  CHECK(::myb_thrust_init(&thrust_c, 100.0, &arc, 1) == MYB_OK);
  const auto thrust = myb::thrust_init(100.0, {arc});
  CHECK(thrust.ok() && myb::thrust_mass(*thrust, 500.0) == ::myb_thrust_mass(&thrust_c, 500.0));
  double rtn_c[3] = {};
  CHECK(::myb_thrust_rtn(&thrust_c, 500.0, y0.data(), y0.data() + 3, rtn_c) == MYB_OK);
  const auto rtn = myb::thrust_rtn(*thrust, 500.0, r1, {y0[3], y0[4], y0[5]});
  CHECK(rtn.ok() && same_array(*rtn, rtn_c));
  CHECK(myb::thrust_init(-1.0, {arc}).code() == MYB_E_DOMAIN);
  double edelbaum_dv = 0.0;
  double beta0 = 0.0;
  CHECK(::myb_edelbaum(7700.0, 3070.0, 0.5, &edelbaum_dv, &beta0) == MYB_OK);
  const auto edelbaum = myb::edelbaum(7700.0, 3070.0, 0.5);
  CHECK(edelbaum.ok() && edelbaum->dv == edelbaum_dv && edelbaum->beta0 == beta0);
  CHECK(myb::edelbaum(0.0, 3070.0, 0.5).code() == MYB_E_DOMAIN);
  ::myb_qlaw law{};
  law.mu = gm;
  law.target[0] = 7.5e6;
  law.target[1] = 0.01;
  law.target[2] = 0.9;
  law.weight[0] = 1.0;
  law.weight[1] = 1.0;
  law.weight[2] = 1.0;
  const auto elements = myb::elements_from_state(gm, y0);
  CHECK(elements.ok());
  double direction_c[3] = {};
  double q_c = 0.0;
  CHECK(::myb_qlaw_direction(&law, &*elements, direction_c, &q_c, nullptr, nullptr) == MYB_OK);
  const auto steering = myb::qlaw_direction(law, *elements);
  CHECK(steering.ok() && same_array(steering->direction, direction_c) && steering->q == q_c);
  ::myb_qlaw bad_law{};
  CHECK(myb::qlaw_direction(bad_law, *elements).code() == MYB_E_DOMAIN);

  // two orbits that pass each other: approaches, the encounter and its probability
  const kepler_motion first_motion{gm, y0};
  const kepler_motion second_motion{gm, deputy};
  ::myb_approach approaches_c[16] = {};
  std::size_t found_c = 0;
  CHECK(::myb_approach_find(kepler_state, const_cast<kepler_motion*>(&first_motion), kepler_state,
                            const_cast<kepler_motion*>(&second_motion), 0.0, 6000.0, 60.0,
                            approaches_c, 16, &found_c) == MYB_OK);
  const auto approaches =
      myb::approach_find(kepler_provider(y0), kepler_provider(deputy), 0.0, 6000.0, 60.0, 16);
  CHECK(approaches.ok() && approaches->found == found_c && found_c >= 1);
  if (approaches.ok() && !approaches->approaches.empty()) {
    const ::myb_approach& closest = approaches->approaches[0];
    CHECK(bits(&closest, &approaches_c[0], sizeof closest));
    const myb::mat3 covariance{100.0, 0.0, 0.0, 0.0, 400.0, 0.0, 0.0, 0.0, 900.0};
    ::myb_encounter encounter_c{};
    CHECK(::myb_encounter_plane(closest.position, closest.velocity, covariance.data(), 10.0,
                                &encounter_c) == MYB_OK);
    const auto encounter = myb::encounter_plane(
        {closest.position[0], closest.position[1], closest.position[2]},
        {closest.velocity[0], closest.velocity[1], closest.velocity[2]}, covariance, 10.0);
    CHECK(encounter.ok() && bits(&*encounter, &encounter_c, sizeof encounter_c));
    double probability_c = 0.0;
    CHECK(::myb_collision_probability(&encounter_c, &probability_c) == MYB_OK);
    const auto probability = myb::collision_probability(*encounter);
    CHECK(probability.ok() && *probability == probability_c);
    std::printf("two orbits through the class: %zu approaches, the closest %.1f m, probability "
                "%.3g, %s the C calls\n",
                approaches->found, closest.distance, probability.value_or(-1.0),
                probability.ok() && *probability == probability_c ? "bit for bit with" : "NOT as");
    CHECK(myb::encounter_plane({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, covariance, 10.0).code() ==
          MYB_E_DOMAIN);
  }
  CHECK(myb::approach_find(kepler_provider(y0), kepler_provider(deputy), 10.0, 0.0, 60.0).code() ==
        MYB_E_DOMAIN);
  // a C provider in the class: the same as the lambda
  const myb::provider c_provider(kepler_state, const_cast<kepler_motion*>(&first_motion));
  const auto from_c = c_provider.at(::myb_time_from_seconds(100.0));
  const auto from_lambda = kepler_provider(y0).at(::myb_time_from_seconds(100.0));
  CHECK(from_c.ok() && from_lambda.ok() && same_bits(*from_c, *from_lambda));

  // attitude
  const myb::quat a{0.9, 0.1, -0.3, 0.2};
  const myb::quat b{0.5, -0.5, 0.5, 0.5};
  double product_c[4] = {};
  ::myb_quat_multiply(a.data(), b.data(), product_c);
  CHECK(same_array(myb::quat_multiply(a, b), product_c));
  double matrix_c[9] = {};
  ::myb_quat_to_matrix(a.data(), matrix_c);
  CHECK(same_array(myb::quat_to_matrix(a), matrix_c));
  myb::quat unit = a;
  double unit_c[4] = {a[0], a[1], a[2], a[3]};
  CHECK(myb::quat_normalize(unit) == ::myb_quat_normalize(unit_c) && same_array(unit, unit_c));
  const std::array<int, 3> sequence{3, 1, 3};
  const int sequence_c[3] = {3, 1, 3};
  double angles_c[3] = {};
  CHECK(::myb_euler_from_matrix(sequence_c, matrix_c, angles_c) == MYB_OK);
  const auto angles = myb::euler_from_matrix(sequence, myb::quat_to_matrix(a));
  CHECK(angles.ok() && same_array(*angles, angles_c));
  CHECK(myb::euler_to_matrix({1, 1, 2}, {0.1, 0.2, 0.3}).code() == MYB_E_DOMAIN);
  const myb::vec3 moments{1.0, 2.0, 3.0};
  const myb::vec3 l{0.1, 0.2, 2.0};
  ::myb_free_rotation free_c{};
  CHECK(::myb_free_rotation_init(moments.data(), unit.data(), l.data(), &free_c) == MYB_OK);
  const auto free_rotation = myb::free_rotation::init(moments, unit, l);
  CHECK(free_rotation.ok());
  if (free_rotation) {
    double q_at_c[4] = {};
    double l_at_c[3] = {};
    ::myb_free_rotation_state(&free_c, 7.5, q_at_c, l_at_c);
    const myb::attitude_state at = free_rotation->at(7.5);
    CHECK(same_array(at.q, q_at_c) && same_array(at.l, l_at_c));
  }
  CHECK(myb::free_rotation::init({0.0, 1.0, 2.0}, unit, l).code() == MYB_E_DOMAIN);
  double torque[3] = {1.0e-3, 0.0, -2.0e-3};
  double q_split_c[4] = {unit[0], unit[1], unit[2], unit[3]};
  double l_split_c[3] = {l[0], l[1], l[2]};
  CHECK(::myb_rotation_split(moments.data(), constant_torque, torque, 4, 0.0, 0.01, 200, q_split_c,
                             l_split_c) == MYB_OK);
  const myb::vec3 torque_cpp{torque[0], torque[1], torque[2]};
  const auto split = myb::rotation_split(
      moments,
      [torque_cpp](double, const myb::quat&, myb::vec3& n) {
        n = torque_cpp;
        return MYB_OK;
      },
      4, 0.0, 0.01, 200, myb::attitude_state{unit, l});
  CHECK(split.ok() && same_array(split->q, q_split_c) && same_array(split->l, l_split_c));
  CHECK(myb::rotation_split(moments, [](double, const myb::quat&, myb::vec3&) { return MYB_OK; },
                            3, 0.0, 0.01, 1, myb::attitude_state{unit, l})
            .code() == MYB_E_DOMAIN);
}

// ---- Observables

myb_status still_station(void*, ::myb_time, double y[6]) {
  y[0] = 6.4e6;
  y[1] = 0.0;
  y[2] = 0.0;
  y[3] = 0.0;
  y[4] = 0.0;
  y[5] = 0.0;
  return MYB_OK;
}

myb_status at_origin(void*, ::myb_time, double y[6]) {
  for (int k = 0; k < 6; ++k) {
    y[k] = 0.0;
  }
  return MYB_OK;
}

void test_observables() {
  const myb::state y0 = initial_state();
  kepler_motion motion{gm, y0};
  ::myb_light_body body_c{gm, at_origin, nullptr};
  ::myb_light_model model_c{299792458.0, 1.0, &body_c, 1};
  const double observer[6] = {6.4e6, 0.0, 0.0, 0.0, 0.0, 0.0};
  const ::myb_time t = ::myb_time_from_seconds(500.0);
  ::myb_light_path path_c{};
  CHECK(::myb_light_time(&model_c, kepler_state, &motion, t, observer, MYB_LIGHT_RECEPTION,
                         &path_c) == MYB_OK);
  myb::light_model model(299792458.0, 1.0);
  model.add_body(gm, myb::provider([](::myb_time, myb::state& y) {
    y = myb::state{};
    return MYB_OK;
  }));
  const auto path = myb::light_time(model, kepler_provider(y0), t,
                                    {observer[0], 0.0, 0.0, 0.0, 0.0, 0.0});
  CHECK(path.ok() && same_time(path->t_target, path_c.t_target) &&
        bits(path->target, path_c.target, sizeof path_c.target) &&
        path->light_time == path_c.light_time && path->delay == path_c.delay);
  double total_c = 0.0;
  ::myb_time bounce_c{};
  ::myb_time transmit_c{};
  CHECK(::myb_round_trip(&model_c, still_station, nullptr, kepler_state, &motion, t, &total_c,
                         &bounce_c, &transmit_c) == MYB_OK);
  const auto trip = myb::round_trip(model, myb::provider(still_station, nullptr),
                                    kepler_provider(y0), t);
  CHECK(trip.ok() && trip->total == total_c && same_time(trip->t_bounce, bounce_c) &&
        same_time(trip->t_transmit, transmit_c));
  std::printf("light time through the class: %.12f s, the round trip %.12f s, %s the C calls\n",
              path.ok() ? path->light_time : 0.0, trip.ok() ? trip->total : 0.0,
              trip.ok() && trip->total == total_c ? "bit for bit with" : "NOT as");
  myb::light_model bad(-1.0, 1.0);
  CHECK(myb::light_time(bad, kepler_provider(y0), t, myb::state{}).code() == MYB_E_DOMAIN);
  const myb::vec3 p{0.6, 0.0, 0.8};
  double aberrated_c[3] = {};
  ::myb_aberration(p.data(), y0.data() + 3, 299792458.0, aberrated_c);
  CHECK(same_array(myb::aberration(p, {y0[3], y0[4], y0[5]}, 299792458.0), aberrated_c));
  double ra_c = 0.0;
  double dec_c = 0.0;
  ::myb_direction_radec(p.data(), &ra_c, &dec_c);
  CHECK(myb::direction_radec(p).ra == ra_c && myb::direction_radec(p).dec == dec_c);
  ::myb_plasma_profile corona{};
  corona.radius = 6.96e8;
  corona.count = 1;
  corona.density[0] = 1.0e14;
  corona.power[0] = 2.0;
  const myb::vec3 p1{1.5e11, 1.0e10, 0.0};
  const myb::vec3 p2{-1.5e11, 2.0e10, 0.0};
  double column_c = 0.0;
  CHECK(::myb_plasma_column(&corona, p1.data(), p2.data(), &column_c) == MYB_OK);
  const auto column = myb::plasma_column(corona, p1, p2);
  CHECK(column.ok() && *column == column_c);
  CHECK(myb::plasma_column(corona, {1.0e11, 0.0, 0.0}, {-1.0e11, 0.0, 0.0}).code() ==
        MYB_E_DOMAIN); // through the centre
}

// ---- Files

void test_files() {
  // a table
  const std::string csv = "t,x\n0,1.5\n1,2.25\n2,-3\n";
  const auto parsed = myb::table::parse(csv);
  CHECK(parsed.ok() && parsed->rows() == 3 && parsed->columns() == 2 && parsed->name(1) == "x");
  if (parsed) {
    const std::vector<double> row = parsed->row(1);
    CHECK(row.size() == 2 && row[0] == 1.0 && row[1] == 2.25);
  }
  const auto broken = myb::table::parse("t,x\n0,1\n1,y\n");
  CHECK(!broken.ok() && broken.code() == MYB_E_FORMAT && !broken.message().empty());
  const std::vector<double> values = {0.0, 1.5, 1.0, 2.25};
  const char* names_c[2] = {"t", "x"};
  char written_c[256] = {};
  CHECK(::myb_table_write(names_c, 2, values.data(), 2, written_c, sizeof written_c, nullptr) ==
        MYB_OK);
  const auto written = myb::table_write({"t", "x"}, 2, values);
  CHECK(written.ok() && *written == written_c);
  CHECK(myb::table_write({"t"}, 2, values).code() == MYB_E_DOMAIN);

  // an SPK file of states and one of Chebyshev series, written in memory and read back
  const myb::state y0 = initial_state();
  std::vector<double> epochs;
  std::vector<myb::state> states;
  for (int k = 0; k < 20; ++k) {
    const double t = 120.0 * k;
    const auto y = myb::kepler_propagate(gm, y0, t);
    CHECK(y.ok());
    epochs.push_back(t);
    states.push_back(y.value_or(myb::state{}));
  }
  ::myb_spk_segment_header header{-1000, 399, 1, "test orbit"};
  ::myb_spk_writer* writer_c = nullptr;
  CHECK(::myb_spk_writer_create("myb test", nullptr, &writer_c) == MYB_OK);
  CHECK(::myb_spk_writer_add_states(writer_c, &header, 13, 5, epochs.data(), states.front().data(),
                                    states.size()) == MYB_OK);
  kepler_motion motion{gm, y0};
  ::myb_spk_segment_header chebyshev_header{-1001, 399, 1, "chebyshev"};
  CHECK(::myb_spk_writer_add_chebyshev(writer_c, &chebyshev_header, 3, kepler_state, &motion, 0.0,
                                       2280.0, 4, 11, nullptr) == MYB_OK);
  std::vector<unsigned char> bytes_c(::myb_spk_writer_size(writer_c));
  CHECK(::myb_spk_writer_write(writer_c, bytes_c.data(), bytes_c.size()) == MYB_OK);
  ::myb_spk_writer_destroy(writer_c);
  auto writer = myb::spk_writer::create("myb test");
  CHECK(writer.ok());
  if (!writer) {
    return;
  }
  CHECK(writer->add_states(header, 13, 5, epochs, states).ok());
  const auto max_error = writer->add_chebyshev_measured(chebyshev_header, 3, kepler_provider(y0),
                                                        0.0, 2280.0, 4, 11);
  CHECK(max_error.ok() && *max_error >= 0.0 && *max_error < 1.0);
  CHECK(writer->add_states(header, 13, 5, epochs, {}).code() == MYB_E_DOMAIN);
  CHECK(writer->add_states(header, 13, 0, epochs, states).code() == MYB_E_DOMAIN);
  const auto bytes = writer->write();
  CHECK(bytes.ok() && *bytes == bytes_c);
  auto loaded = myb::ephemeris::create();
  CHECK(loaded.ok() && bytes.ok() && loaded->load_memory(*bytes).ok());
  CHECK(loaded.ok() && loaded->segment_count() == 2 && loaded->segment(1).target == -1001);
  if (loaded.ok()) {
    const auto back = loaded->at(-1000, 399, ::myb_time_from_seconds(600.0));
    CHECK(back.ok() && std::fabs((*back)[0] - states[5][0]) < 1e-6);
    std::printf("an SPK file of types 13 and 3 through the class: %zu bytes, %s the C writer, the "
                "Chebyshev fit within %.2e m\n",
                bytes.ok() ? bytes->size() : 0, bytes.ok() && *bytes == bytes_c ? "the bytes of" : "NOT as",
                max_error.value_or(-1.0));
  }

  // SP3
  const char* sp3_path = MYB_TEST_DATA_DIR "/sp3/ilrsa.lageos1.250607.excerpt.sp3";
  ::myb_sp3* sp3_c = nullptr;
  CHECK(::myb_sp3_load_file(sp3_path, nullptr, &sp3_c, nullptr) == MYB_OK);
  const auto orbits = myb::sp3::load_file(sp3_path);
  CHECK(orbits.ok());
  if (orbits && sp3_c != nullptr) {
    CHECK(orbits->satellite_count() == ::myb_sp3_satellite_count(sp3_c) &&
          orbits->satellite(0) == "L51" && orbits->frame() == ::myb_sp3_frame(sp3_c));
    double position_c[3] = {};
    double velocity_c[3] = {};
    const int given_c = ::myb_sp3_state(sp3_c, 0, 3, position_c, velocity_c);
    const myb::sp3::record record = orbits->at(0, 3);
    CHECK(record.given == given_c && same_array(record.position, position_c));
    ::myb_calendar epoch_c{};
    ::myb_sp3_epoch(sp3_c, 3, &epoch_c);
    const ::myb_calendar epoch = orbits->epoch(3);
    CHECK(bits(&epoch, &epoch_c, sizeof epoch_c));
  }
  ::myb_sp3_destroy(sp3_c);
  CHECK(myb::sp3::load_file(MYB_TEST_DATA_DIR "/sp3/missing.sp3").code() == MYB_E_IO);

  // a spacecraft clock and an orientation file
  auto pool = myb::kernel_pool::create();
  CHECK(pool.ok() && pool->load_file(MYB_TEST_DATA_DIR "/ck/test.tsc").ok());
  if (!pool) {
    return;
  }
  ::myb_sclk* clock_c = nullptr;
  CHECK(::myb_sclk_create(pool->get(), 9099, nullptr, &clock_c, nullptr) == MYB_OK);
  const auto clock = myb::share(myb::sclk::create(*pool, 9099));
  CHECK(clock.ok() && myb::sclk::create(*pool, 9097).code() == MYB_E_NOT_FOUND);
  if (clock && clock_c != nullptr) {
    ::myb_time t_c{};
    CHECK(::myb_sclk_to_time(clock_c, 6.1e6, &t_c) == MYB_OK);
    const auto t = (*clock)->to_time(6.1e6);
    CHECK(t.ok() && same_time(*t, t_c) && (*clock)->partitions() == ::myb_sclk_partitions(clock_c));
    const auto text = (*clock)->decode(6.1e6);
    char text_c[64] = {};
    CHECK(::myb_sclk_decode(clock_c, 6.1e6, text_c, sizeof text_c, nullptr) == MYB_OK);
    CHECK(text.ok() && *text == text_c);
    const auto ticks = (*clock)->encode(text.value_or(""));
    CHECK(ticks.ok() && *ticks == 6.1e6);
    CHECK((*clock)->encode("not a clock").code() == MYB_E_FORMAT);
  }
  ::myb_sclk_destroy(clock_c);
  const std::vector<double> ticks = {1000.0, 2000.0, 3000.0, 4000.0};
  std::vector<myb::quat> quaternions;
  for (int k = 0; k < 4; ++k) {
    const double angle = 0.1 * k;
    quaternions.push_back({std::cos(angle), 0.0, 0.0, std::sin(angle)});
  }
  ::myb_ck_segment_header ck_header{-77001, 1, "pointings"};
  ::myb_ck_writer* ck_writer_c = nullptr;
  CHECK(::myb_ck_writer_create("myb test", nullptr, &ck_writer_c) == MYB_OK);
  CHECK(::myb_ck_writer_add_pointings(ck_writer_c, &ck_header, 3, ticks.data(),
                                      quaternions.front().data(), nullptr, ticks.size(), nullptr,
                                      0) == MYB_OK);
  std::vector<unsigned char> ck_bytes_c(::myb_ck_writer_size(ck_writer_c));
  CHECK(::myb_ck_writer_write(ck_writer_c, ck_bytes_c.data(), ck_bytes_c.size()) == MYB_OK);
  ::myb_ck_writer_destroy(ck_writer_c);
  auto ck_writer = myb::ck_writer::create("myb test");
  CHECK(ck_writer.ok() && ck_writer->add_pointings(ck_header, 3, ticks, quaternions).ok());
  if (!ck_writer) {
    return;
  }
  CHECK(ck_writer->add_pointings(ck_header, 3, ticks, {}).code() == MYB_E_DOMAIN);
  const auto ck_bytes = ck_writer->write();
  CHECK(ck_bytes.ok() && *ck_bytes == ck_bytes_c);
  auto files = myb::share(myb::ck::create());
  CHECK(files.ok() && ck_bytes.ok() && (*files)->load_memory(*ck_bytes).ok());
  if (files && ck_bytes.ok()) {
    const auto pointing = (*files)->pointing_at(-77001, 2500.0, 0.0);
    CHECK(pointing.ok() && (*files)->segment_count() == 1 && (*files)->segment(0).type == 3);
    CHECK((*files)->pointing_at(-1, 2500.0, 0.0).code() == MYB_E_RANGE);
  }

  // an OPM written and read back, in both notations
  ::myb_opm opm{};
  std::snprintf(opm.originator, sizeof opm.originator, "%s", "myb");
  std::snprintf(opm.object_name, sizeof opm.object_name, "%s", "SAT");
  std::snprintf(opm.object_id, sizeof opm.object_id, "%s", "2026-001A");
  std::snprintf(opm.center_name, sizeof opm.center_name, "%s", "EARTH");
  std::snprintf(opm.ref_frame, sizeof opm.ref_frame, "%s", "EME2000");
  std::snprintf(opm.time_system, sizeof opm.time_system, "%s", "UTC");
  opm.creation_date = {2026, 3, 2, 0, 0, 0.0};
  opm.epoch = {2026, 3, 1, 12, 0, 0.0};
  for (std::size_t k = 0; k < 6; ++k) {
    opm.state[k] = y0[k];
  }
  for (const ::myb_odm_notation notation : {MYB_ODM_KVN, MYB_ODM_XML}) {
    const auto text = myb::to_text(opm, notation);
    CHECK(text.ok());
    const auto back = myb::parse_opm(text.value_or(""));
    CHECK(back.ok() && std::fabs(back->state[0] - y0[0]) < 1e-6 &&
          std::string(back->object_name) == "SAT");
  }
  const auto not_opm = myb::parse_opm("CCSDS_OPM_VERS = 2.0\n");
  CHECK(!not_opm.ok() && not_opm.code() == MYB_E_FORMAT && !not_opm.message().empty());
  const double line_values[2] = {1.5, -2.25};
  char line_c[128] = {};
  CHECK(::myb_ocm_format_line(nullptr, 30.0, line_values, 2, line_c, sizeof line_c, nullptr) ==
        MYB_OK);
  const auto line = myb::ocm_format_line(std::nullopt, 30.0, {1.5, -2.25});
  CHECK(line.ok() && *line == line_c);

  // numbers in text and the variables of a kernel pool
  double number_c = 0.0;
  std::size_t consumed_c = 0;
  CHECK(::myb_parse_double("6.02214076D+23 rest", 19, &number_c, &consumed_c) == MYB_OK);
  const auto number = myb::parse_double("6.02214076D+23 rest");
  CHECK(number.ok() && number->value == number_c && number->consumed == consumed_c);
  CHECK(myb::parse_double("x1").code() == MYB_E_FORMAT);
  char shortest_c[MYB_DOUBLE_CHARS] = {};
  const std::size_t shortest_length = ::myb_format_double(0.1 + 0.2, shortest_c);
  CHECK(myb::format_double(0.1 + 0.2) == std::string(shortest_c, shortest_length));
  const auto fixed = myb::format_fixed(2.5, 0);
  const auto scientific = myb::format_scientific(123456.0, 4);
  CHECK(fixed.ok() && *fixed == "2" && scientific.ok() && *scientific == "1.235e+05");
  CHECK(myb::format_fixed(1.0, 99).code() == MYB_E_DOMAIN);
  auto kernels = myb::kernel_pool::create();
  CHECK(kernels.ok() &&
        kernels->load_text("\\begindata\nVALUES = ( 1.5 2.5 )\nNAME = 'myb'\n\\begintext\n").ok());
  if (kernels) {
    const auto described = kernels->describe("VALUES");
    CHECK(described.ok() && described->type == MYB_KVAR_NUMBERS && described->count == 2);
    const auto values_read = kernels->numbers("VALUES");
    CHECK(values_read.ok() && *values_read == std::vector<double>({1.5, 2.5}));
    const auto name = kernels->string("NAME");
    CHECK(name.ok() && *name == "myb");
    CHECK(kernels->numbers("NAME").code() == MYB_E_FORMAT &&
          kernels->describe("MISSING").code() == MYB_E_NOT_FOUND);
    const auto bad = kernels->load_text("\\begindata\nX = (\n");
    CHECK(!bad.ok() && bad.code() == MYB_E_FORMAT && !bad.message().empty());
  }
}

// ---- The array forms

myb_status turning_orientation(void*, ::myb_time t, ::myb_orientation* orientation) {
  const double angles[3] = {0.0, 0.0, 1.0e-4 * ::myb_time_seconds(t)};
  const double rates[3] = {0.0, 0.0, 1.0e-4};
  ::myb_rotation_euler313(angles, rates, &orientation->rotation);
  orientation->ut = t;
  orientation->xp = 0.0;
  orientation->yp = 0.0;
  return MYB_OK;
}

void test_array_forms() {
  // instants
  const std::vector<double> seconds = {0.5, 7.0e8 + 0.25, -3.75, 812345678.125};
  const std::vector<::myb_time> instants = myb::time_from_seconds(seconds);
  bool same = instants.size() == seconds.size();
  for (std::size_t k = 0; same && k < seconds.size(); ++k) {
    same = same_time(instants[k], ::myb_time_from_seconds(seconds[k]));
  }
  CHECK(same && myb::time_seconds(instants).size() == 4);
  const std::vector<::myb_calendar> labels = {{2026, 1, 1, 0, 0, 0.0}, {2026, 2, 28, 23, 59, 59.5}};
  const auto converted = myb::calendar_to_time(labels);
  ::myb_time second_c{};
  CHECK(::myb_calendar_to_time(&labels[1], &second_c) == MYB_OK);
  CHECK(converted.ok() && same_time((*converted)[1], second_c));
  const std::vector<::myb_calendar> wrong = {labels[0], labels[1], {2026, 2, 30, 0, 0, 0.0}};
  const auto failed = myb::calendar_to_time(wrong);
  CHECK(!failed.ok() && failed.code() == MYB_E_DOMAIN && failed.message() == "element 2");
  std::printf("an array form that stops: %s (%s)\n", myb::status_name(failed.code()),
              failed.message().c_str());
  const auto calendars = myb::time_to_calendar(instants);
  ::myb_calendar calendar_c{};
  ::myb_time_to_calendar(instants[1], &calendar_c);
  CHECK(calendars.size() == 4 && bits(&calendars[1], &calendar_c, sizeof calendar_c));
  // time scales
  auto pool = myb::kernel_pool::create();
  CHECK(pool.ok());
  if (!pool) {
    return;
  }
  for (const char* kernel :
       {MYB_DATA_DIR "/time/time-scales.tk", MYB_DATA_DIR "/time/tdb-tt-fb1990.tk",
        MYB_DATA_DIR "/time/utc-tai.tk", MYB_DATA_DIR "/frames/iau1980.tk",
        MYB_DATA_DIR "/frames/iau2006-cip.tk", MYB_TEST_DATA_DIR "/wgccre/pck00011.tpc"}) {
    CHECK(pool->load_file(kernel).ok());
  }
  auto system = myb::share(myb::time_system::create(*pool));
  CHECK(system.ok());
  if (!system) {
    return;
  }
  const auto tt = (*system)->scale("TT");
  const auto tdb = (*system)->scale("TDB");
  CHECK(tt.ok() && tdb.ok());
  const auto in_tdb = (*system)->convert(instants, *tt, *tdb);
  ::myb_time tdb_c{};
  CHECK(::myb_timesys_convert((*system)->get(), instants[3], *tt, *tdb, &tdb_c) == MYB_OK);
  CHECK(in_tdb.ok() && same_time((*in_tdb)[3], tdb_c));
  const auto utc = (*system)->to_utc(instants, *tt);
  CHECK(utc.ok() && utc->size() == 4);
  const auto from_utc = (*system)->from_utc(labels, *tt);
  CHECK(from_utc.ok() && from_utc->size() == 2);
  // two bodies
  const myb::state y0 = initial_state();
  const std::vector<double> dt = {0.0, 100.0, -2500.0, 86400.0};
  const auto propagated = myb::kepler_propagate(gm, y0, dt);
  const auto single = myb::kepler_propagate(gm, y0, 86400.0);
  CHECK(propagated.ok() && single.ok() && same_bits((*propagated)[3], *single));
  const auto transitions = myb::kepler_stm(gm, y0, dt);
  const auto one = myb::kepler_stm(gm, y0, -2500.0);
  CHECK(transitions.ok() && one.ok() && same_bits((*transitions)[2].y, one->y) &&
        (*transitions)[2].stm == one->stm);
  const myb::state still{7.0e6, 0.0, 0.0, 0.0, 0.0, 0.0};
  const auto elements = myb::elements_from_state(gm, std::vector<myb::state>{y0, still, y0});
  CHECK(!elements.ok() && elements.code() == MYB_E_DOMAIN && elements.message() == "element 1");
  const auto equinoctial = myb::equinoctial_from_state(gm, std::vector<myb::state>{y0});
  const auto equinoctial_one = myb::equinoctial_from_state(gm, y0);
  CHECK(equinoctial.ok() && equinoctial_one.ok() &&
        bits(&(*equinoctial)[0], &*equinoctial_one, sizeof(::myb_equinoctial)));
  const auto round = myb::state_from_equinoctial(gm, *equinoctial);
  CHECK(round.ok() && round->size() == 1);
  // ephemerides at a series of instants
  auto eph = myb::share(myb::ephemeris::create());
  CHECK(eph.ok() && (*eph)->load_file(MYB_TEST_DATA_DIR "/spk/t02_03.bsp").ok());
  if (eph) {
    const ::myb_spk_segment_info info = (*eph)->segment(0);
    std::vector<::myb_time> times;
    for (int k = 0; k < 5; ++k) {
      times.push_back(::myb_time_from_seconds(info.start + (info.stop - info.start) * (0.1 + 0.2 * k)));
    }
    const auto series = (*eph)->at(info.target, info.center, times);
    const auto alone = (*eph)->at(info.target, info.center, times[4]);
    CHECK(series.ok() && alone.ok() && same_bits((*series)[4], *alone));
    times.push_back(::myb_time_from_seconds(info.stop + 1.0e9));
    const auto beyond = (*eph)->at(info.target, info.center, times);
    CHECK(!beyond.ok() && beyond.code() == MYB_E_RANGE && beyond.message() == "element 5");
    const myb::provider body = myb::spk_provider(*eph, info.target, info.center);
    const auto from_provider = body.at(times[4]);
    CHECK(from_provider.ok() && same_bits(*from_provider, *alone));
  }
  // rotations
  const std::vector<::myb_rotation> rotations = {identity_rotation(), identity_rotation()};
  const auto turned = myb::rotation_apply_state(rotations, std::vector<myb::state>{y0, y0});
  CHECK(turned.ok() && same_bits((*turned)[1], myb::rotation_apply_state(rotations[1], y0)));
  CHECK(myb::rotation_apply_state(rotations, std::vector<myb::state>{y0}).code() == MYB_E_DOMAIN);
  const auto equinox = myb::equinox_model::create(*pool);
  CHECK(equinox.ok());
  if (equinox) {
    const std::vector<::myb_rotation> tod = equinox->rotation(MYB_EQUINOX_TOD, instants);
    ::myb_rotation tod_c{};
    ::myb_equinox_rotation(equinox->get(), MYB_EQUINOX_TOD, instants[2], &tod_c);
    CHECK(tod.size() == 4 && same_rotation(tod[2], tod_c));
    const std::vector<::myb_rotation> pef = equinox->teme_pef(instants);
    CHECK(pef.size() == 4 && same_rotation(pef[1], equinox->teme_pef(instants[1])));
  }
  const auto cip = myb::cip_model::create(*pool);
  CHECK(cip.ok());
  if (cip) {
    const myb::cip_model::poles poles = cip->xys(instants);
    double xys_c[3] = {};
    double rates_c[3] = {};
    ::myb_cip_xys(cip->get(), instants[1], xys_c, rates_c);
    CHECK(poles.xys.size() == 4 && same_array(poles.xys[1], xys_c) && same_array(poles.rates[1], rates_c));
    const auto terrestrial = cip->terrestrial(instants, instants);
    CHECK(terrestrial.ok() && same_rotation((*terrestrial)[0], cip->terrestrial(instants[0], instants[0])));
    CHECK(cip->terrestrial(instants, std::vector<::myb_time>(1)).code() == MYB_E_DOMAIN);
  }
  CHECK(myb::cip_model::create(*pool, "NO_SUCH_MODEL").code() == MYB_E_NOT_FOUND);
  const auto mars = myb::share(myb::wgccre_model::create(*pool, 499));
  CHECK(mars.ok());
  if (mars) {
    const std::vector<::myb_rotation> series = (*mars)->rotation(instants);
    ::myb_rotation mars_c{};
    ::myb_wgccre_rotation((*mars)->get(), instants[3], &mars_c);
    CHECK(series.size() == 4 && same_rotation(series[3], mars_c) && (*mars)->frame() == 1);
    const auto orientations = myb::orientation_many(myb::wgccre_orientation(*mars), instants);
    CHECK(orientations.ok() && same_rotation((*orientations)[3].rotation, mars_c));
  }
  // EOP over a series
  auto eop = myb::eop_series::create();
  CHECK(eop.ok() && eop->load_c04_file(MYB_TEST_DATA_DIR "/eop/eopc04-excerpt.txt").ok());
  if (eop) {
    CHECK(eop->rows().rows > 100);
    const auto june = (*system)->from_utc(std::vector<::myb_calendar>{{2018, 6, 1, 0, 0, 0.0},
                                                                      {2018, 6, 2, 6, 0, 0.0}},
                                          *tt);
    CHECK(june.ok());
    if (june) {
      const auto series = eop->at(**system, *june, *tt);
      const auto alone = eop->at(**system, (*june)[1], *tt);
      CHECK(series.ok() && alone.ok() &&
            bits(&series->values[1], &alone->values, sizeof(::myb_eop_values)) &&
            same_time(series->ut1[1], alone->ut1));
    }
    CHECK(eop->load_c04("not a table").code() == MYB_E_FORMAT);
  }
  // a graph of frames: a frame that turns, from a lambda and from C
  auto graph = myb::frame_graph::create();
  CHECK(graph.ok());
  if (graph) {
    CHECK(graph->add(10, 1, myb::orientation_provider([](::myb_time t, ::myb_orientation& o) {
            return turning_orientation(nullptr, t, &o);
          })).ok());
    CHECK(graph->add(20, 1, myb::orientation_provider(turning_orientation, nullptr)).ok());
    CHECK(graph->add_fixed(30, 10, identity_rotation()).ok());
    CHECK(graph->add(10, 1, myb::fixed_orientation(identity_rotation())).code() == MYB_E_DOMAIN);
    const auto series = graph->rotation(30, 20, instants);
    const auto alone = graph->rotation(30, 20, instants[3]);
    CHECK(series.ok() && alone.ok() && same_rotation((*series)[3], *alone));
    CHECK(graph->rotation(30, 99, instants).code() == MYB_E_NOT_FOUND);
    const auto through_provider = graph->orientation(10, 1).at(instants[2]);
    const auto direct = graph->rotation(10, 1, instants[2]);
    CHECK(through_provider.ok() && direct.ok() && same_rotation(through_provider->rotation, *direct));
  }
}

// ---- Contact, breaks, the timescale of IAS15, the orders of a PN pair

myb_status moon_gravity(void*, double, const double[3], const double[3], double a[3]) {
  a[0] = 0.0;
  a[1] = 0.0;
  a[2] = -1.62;
  return MYB_OK;
}

bool same_contact(const ::myb_contact_result& a, const ::myb_contact_result& b) {
  return bits(a.force, b.force, sizeof a.force) && bits(a.torque, b.torque, sizeof a.torque) &&
         bits(&a.energy, &b.energy, sizeof(double)) && bits(&a.gap, &b.gap, sizeof(double)) &&
         a.touching == b.touching;
}

void test_contact_and_breaks() {
  ::myb_contact_surface ground{};
  ground.shape = MYB_CONTACT_PLANE;
  ground.normal[2] = 1.0;
  ::myb_contact_law law{};
  law.stiffness = 1.0e6;
  law.exponent = 1.5;
  law.damping = 2.0;
  // a ball dropped on the plane by IAS15 with the timescale of the jerk, its steps aimed at the
  // touches (tests/unit/test_contact.c)
  ::myb_orbit_options options{};
  options.epsilon = 1e-10;
  options.timescale = MYB_IAS15_TIMESCALE_JERK;
  const double y0_c[6] = {0.0, 0.0, 0.3, 0.0, 0.0, 0.0};
  ::myb_contact_term term_c{};
  term_c.surface = ground;
  term_c.law = law;
  term_c.radius = 0.1;
  ::myb_force_model* model_c = nullptr;
  CHECK(::myb_force_model_create(nullptr, &model_c) == MYB_OK);
  CHECK(::myb_force_model_add(model_c, "gravity", moon_gravity, nullptr, nullptr, nullptr) ==
        MYB_OK);
  CHECK(::myb_force_model_add(model_c, "contact", ::myb_term_contact, nullptr, &term_c,
                              nullptr) == MYB_OK);
  CHECK(::myb_force_model_add_break(model_c, ::myb_event_contact, &term_c) == MYB_OK);
  ::myb_force_model_set_break_scheme(model_c, MYB_BREAKS_PREDICT);
  ::myb_orbit* orbit_c = nullptr;
  CHECK(::myb_orbit_create(model_c, &options, nullptr, &orbit_c) == MYB_OK);
  CHECK(::myb_orbit_start(orbit_c, 0.0, y0_c, 2.0) == MYB_OK);
  double end_c[6] = {};
  CHECK(::myb_orbit_state(orbit_c, 2.0, end_c) == MYB_OK);
  const std::int64_t steps_c = ::myb_orbit_steps(orbit_c);
  const std::int64_t breaks_c = ::myb_orbit_breaks(orbit_c);
  ::myb_orbit_destroy(orbit_c);
  ::myb_force_model_destroy(model_c);
  auto model = myb::share(myb::force_model::create());
  CHECK(model.ok());
  if (!model) {
    return;
  }
  CHECK((*model)->break_scheme() == MYB_BREAKS_RETAKE);
  CHECK((*model)
            ->add_term("gravity",
                       [](double, const myb::vec3&, const myb::vec3&, myb::vec3& a) {
                         a = {0.0, 0.0, -1.62};
                         return MYB_OK;
                       })
            .ok());
  CHECK((*model)->add_contact("contact", nullptr, 0, ground, law, 0.1, true).ok());
  CHECK((*model)->break_count() == 1);
  (*model)->set_break_scheme(MYB_BREAKS_PREDICT);
  CHECK((*model)->break_scheme() == MYB_BREAKS_PREDICT);
  const myb::state y0{0.0, 0.0, 0.3, 0.0, 0.0, 0.0};
  auto orbit = myb::orbit::create(*model, options);
  CHECK(orbit.ok());
  if (orbit) {
    CHECK(orbit->start(0.0, y0, 2.0).ok());
    const auto end = orbit->at(2.0);
    CHECK(end.ok() && same_array(*end, end_c) && orbit->steps() == steps_c &&
          orbit->breaks() == breaks_c && breaks_c > 0);
    std::printf("a ball dropped on a plane through the class (IAS15, the jerk, aimed breaks): "
                "%lld steps, %lld breaks, %s the C calls\n",
                static_cast<long long>(orbit->steps()), static_cast<long long>(orbit->breaks()),
                end.ok() && same_array(*end, end_c) ? "bit for bit with" : "NOT as");
  }
  // the timescale of a alone stalls at the touch; an unknown timescale is refused
  ::myb_orbit_options plain{};
  plain.epsilon = 1e-10;
  auto stalled = myb::orbit::create(*model, plain);
  CHECK(stalled.ok() && stalled->start(0.0, y0, 2.0).ok() &&
        stalled->at(2.0).code() == MYB_E_STEP_MIN);
  ::myb_orbit_options unknown = options;
  const int seven = 7;
  static_assert(sizeof unknown.timescale == sizeof seven, "an enum of C is an int here");
  std::memcpy(&unknown.timescale, &seven, sizeof seven);
  CHECK(myb::orbit::create(*model, unknown).code() == MYB_E_DOMAIN);
  // the term with a friction that sticks needs its state: the model refuses it
  ::myb_contact_law stick = law;
  stick.friction = 0.6;
  stick.kind = MYB_CONTACT_FRICTION_STICK;
  stick.tangential_stiffness = 1.0e4;
  stick.tangential_damping = 100.0;
  auto refusing = myb::force_model::create();
  CHECK(refusing.ok() && refusing->add_contact("contact", nullptr, 0, ground, stick, 0.1).ok());
  if (refusing) {
    CHECK(refusing->acceleration(0.0, {0.0, 0.0, 0.099}, {0.0, 0.0, 0.0}).code() == MYB_E_DOMAIN);
  }

  // the force of a surface, smooth and sticking
  const std::vector<::myb_contact_point> points = {{{0.0, 0.0, 0.0}, 0.1}};
  const myb::vec3 x{0.0, 0.0, 0.099};
  const myb::vec3 v{0.01, 0.0, -0.001};
  ::myb_contact_law smooth = law;
  smooth.friction = 0.5;
  smooth.slip = 1e-3;
  ::myb_contact_result smooth_c{};
  CHECK(::myb_contact_force(&ground, &smooth, points.data(), 1, x.data(), v.data(), nullptr,
                            nullptr, &smooth_c) == MYB_OK);
  const auto smooth_force = myb::contact_force(ground, smooth, points, x, v);
  CHECK(smooth_force.ok() && same_contact(*smooth_force, smooth_c) && smooth_c.touching == 1);
  CHECK(myb::contact_force(ground, stick, points, x, v).code() == MYB_E_DOMAIN);
  const std::vector<double> deflection = {1.0e-5, -2.0e-6, 0.0};
  double rate_c[3] = {};
  ::myb_contact_result stick_c{};
  CHECK(::myb_contact_force_stick(&ground, &stick, points.data(), 1, x.data(), v.data(), nullptr,
                                  nullptr, deflection.data(), rate_c, &stick_c) == MYB_OK);
  const auto stick_force = myb::contact_force_stick(ground, stick, points, x, v, deflection);
  CHECK(stick_force.ok() && same_contact(stick_force->result, stick_c) &&
        bits(stick_force->deflection_rate.data(), rate_c, sizeof rate_c));
  CHECK(myb::contact_force_stick(ground, stick, points, x, v, {1.0}).code() == MYB_E_DOMAIN);
  CHECK(myb::contact_force_stick(ground, stick, points, x, v, {}).code() == MYB_E_DOMAIN);
  double restitution_c = 0.0;
  CHECK(::myb_contact_restitution(2.0, 1.5, &restitution_c) == MYB_OK);
  const auto restitution = myb::contact_restitution(2.0, 1.5);
  CHECK(restitution.ok() && *restitution == restitution_c);
  CHECK(myb::contact_restitution(-1.0, 1.0).code() == MYB_E_DOMAIN);
  const auto nearest = myb::contact_surface_nearest(ground, {1.0, 2.0, 0.5});
  CHECK(nearest.ok() && nearest->height == 0.5);

  // a system with a contact that sticks: the deflections in the state
  const double gms[2] = {0.0, 1.0};
  ::myb_system* s_c = nullptr;
  CHECK(::myb_system_create(2, gms, nullptr, &s_c) == MYB_OK);
  const ::myb_contact_point legs[2] = {{{0.0, 0.0, -0.5}, 0.1}, {{0.3, 0.0, -0.5}, 0.1}};
  ::myb_system_contact contact_c{};
  contact_c.surface = ground;
  contact_c.law = stick;
  contact_c.points = legs;
  contact_c.point_count = 2;
  CHECK(::myb_system_add_contact(s_c, 1, 0, &contact_c) == MYB_OK);
  const std::size_t n = ::myb_system_coordinates(s_c);
  std::vector<double> y_c(2 * n, 0.0);
  y_c[5] = 0.59;
  y_c[n + 3] = 0.02;
  const std::vector<double> springs = {1.0e-5, 0.0, 0.0, -1.0e-5, 2.0e-6, 0.0};
  CHECK(::myb_system_pack_contact(s_c, 0, springs.data(), y_c.data()) == MYB_OK);
  std::vector<double> a_c(n);
  CHECK(::myb_system_acceleration(s_c, 0.0, y_c.data(), y_c.data() + n, a_c.data()) == MYB_OK);
  ::myb_system_destroy(s_c);
  auto s = myb::system::create({0.0, 1.0});
  CHECK(s.ok() && s->add_contact(1, 0, ground, stick, {legs[0], legs[1]}).ok());
  if (s) {
    CHECK(s->coordinates() == n && n == 12);
    std::vector<double> y(2 * n, 0.0);
    y[5] = 0.59;
    y[n + 3] = 0.02;
    CHECK(s->pack_contact(0, springs, y).ok() && same_doubles(y, y_c));
    const auto back = s->unpack_contact(0, y);
    CHECK(back.ok() && same_doubles(*back, springs));
    const auto a = s->acceleration(0.0, std::vector<double>(y.begin(), y.begin() + 12),
                                   std::vector<double>(y.begin() + 12, y.end()));
    CHECK(a.ok() && same_doubles(*a, a_c));
    CHECK(s->pack_contact(1, springs, y).code() == MYB_E_DOMAIN);
    CHECK(s->pack_contact(0, {1.0, 2.0, 3.0}, y).code() == MYB_E_DOMAIN);
    CHECK(s->unpack_contact(0, std::vector<double>(3)).code() == MYB_E_DOMAIN);
  }
  auto smooth_system = myb::system::create({0.0, 1.0});
  CHECK(smooth_system.ok() &&
        smooth_system->add_contact(1, 0, ground, smooth, {legs[0], legs[1]}).ok());
  if (smooth_system) {
    std::vector<double> y(2 * smooth_system->coordinates(), 0.0);
    CHECK(smooth_system->coordinates() == 6);
    CHECK(smooth_system->pack_contact(0, springs, y).code() == MYB_E_DOMAIN); // no springs
  }

  // the integrator of a system with the timescale of the jerk, and an unknown one
  ::myb_ias15_options jerk{};
  jerk.epsilon = 1e-10;
  jerk.timescale = MYB_IAS15_TIMESCALE_JERK;
  const double pair_gm[2] = {4.0e5, 1.0e5};
  const double start_c[12] = {-1000.0, 0.0, 0.0, 4000.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, -0.04, 0.0};
  ::myb_system* p_c = nullptr;
  CHECK(::myb_system_create(2, pair_gm, nullptr, &p_c) == MYB_OK);
  ::myb_ias15* in_c = nullptr;
  CHECK(::myb_system_integrator(p_c, &jerk, nullptr, &in_c) == MYB_OK);
  CHECK(::myb_ias15_start(in_c, ::myb_system_acceleration, p_c, 0.0, start_c) == MYB_OK);
  CHECK(::myb_ias15_integrate(in_c, 300.0) == MYB_OK);
  const std::vector<double> after_c(::myb_ias15_state(in_c), ::myb_ias15_state(in_c) + 12);
  ::myb_ias15_destroy(in_c);
  ::myb_system_destroy(p_c);
  auto p = myb::system::create({4.0e5, 1.0e5});
  CHECK(p.ok());
  if (p) {
    auto in = p->integrator(&jerk);
    CHECK(in.ok());
    if (in) {
      CHECK(p->start(*in, 0.0, std::vector<double>(start_c, start_c + 12)).ok() &&
            in->integrate(300.0).ok() && same_doubles(in->current(), after_c));
    }
    ::myb_ias15_options wrong = jerk;
    std::memcpy(&wrong.timescale, &seven, sizeof seven);
    CHECK(p->integrator(&wrong).code() == MYB_E_DOMAIN);
  }

  // a pair beyond 1PN with the spin-spin terms of the next order: its spins in the state
  const double binary_gm[2] = {1.3e20, 1.0e20};
  ::myb_system_pn_pair pair{};
  pair.orders = MYB_PN_2 | MYB_PN_SO | MYB_PN_SS3;
  pair.c = 299792458.0;
  pair.kappa1 = 1.0;
  pair.kappa2 = 1.0;
  const double chi1[3] = {0.0, 0.0, 0.5};
  const double chi2[3] = {0.0, 0.3, 0.0};
  ::myb_system* b_c = nullptr;
  CHECK(::myb_system_create(2, binary_gm, nullptr, &b_c) == MYB_OK);
  CHECK(::myb_system_add_pn_pair(b_c, 0, 1, &pair) == MYB_OK);
  const std::size_t nb = ::myb_system_coordinates(b_c);
  std::vector<double> yb_c(2 * nb, 0.0);
  yb_c[0] = -4.0e8;
  yb_c[3] = 6.0e8;
  yb_c[nb + 1] = -2.0e5;
  yb_c[nb + 4] = 2.6e5;
  CHECK(::myb_system_pack_spins(b_c, 0, chi1, chi2, yb_c.data()) == MYB_OK);
  std::vector<double> ab_c(nb);
  CHECK(::myb_system_acceleration(b_c, 0.0, yb_c.data(), yb_c.data() + nb, ab_c.data()) == MYB_OK);
  ::myb_system_destroy(b_c);
  auto binary = myb::system::create({1.3e20, 1.0e20});
  CHECK(binary.ok() && binary->add_pn_pair(0, 1, pair).ok());
  if (binary) {
    CHECK(binary->coordinates() == nb && nb == 12);
    std::vector<double> yb(2 * nb, 0.0);
    yb[0] = -4.0e8;
    yb[3] = 6.0e8;
    yb[nb + 1] = -2.0e5;
    yb[nb + 4] = 2.6e5;
    CHECK(binary->pack_spins(0, {chi1[0], chi1[1], chi1[2]}, {chi2[0], chi2[1], chi2[2]}, yb).ok());
    const auto ab = binary->acceleration(0.0, std::vector<double>(yb.begin(), yb.begin() + 12),
                                         std::vector<double>(yb.begin() + 12, yb.end()));
    CHECK(ab.ok() && same_doubles(*ab, ab_c));
    ::myb_system_pn_pair unknown_order = pair;
    unknown_order.orders = 1u << 20;
    CHECK(binary->add_pn_pair(0, 1, unknown_order).code() == MYB_E_DOMAIN);
    std::printf("a pair with MYB_PN_SS3, a contact that sticks, IAS15 by the jerk through the "
                "class: %s the C calls\n",
                ab.ok() && same_doubles(*ab, ab_c) ? "bit for bit with" : "NOT as");
  }
}

#if defined(MYB_CPP_TEST_SGP4)
void test_two_line_elements() {
  const std::string line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
  const std::string line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
  ::myb_tle tle_c{};
  CHECK(::myb_tle_parse(line1.c_str(), line2.c_str(), 0, &tle_c, nullptr) == MYB_OK);
  const auto tle = myb::tle_parse(line1, line2);
  CHECK(tle.ok() && tle->catalog_number == 5 && tle->mean_motion == tle_c.mean_motion);
  if (!tle) {
    return;
  }
  const auto lines = myb::tle_format(*tle);
  CHECK(lines.ok() && lines->line1 == line1 && lines->line2 == line2);
  ::myb_sgp4 record_c{};
  CHECK(::myb_sgp4_init(&tle_c, MYB_SGP4_WGS72, MYB_SGP4_IMPROVED, &record_c) == MYB_OK);
  double r_c[3] = {};
  double v_c[3] = {};
  CHECK(::myb_sgp4_propagate(&record_c, 360.0, r_c, v_c, nullptr) == MYB_OK);
  const auto theory = myb::sgp4::init(*tle);
  CHECK(theory.ok());
  if (theory) {
    const auto y = theory->at(360.0);
    CHECK(y.ok() && bits(y->data(), r_c, sizeof r_c) && bits(y->data() + 3, v_c, sizeof v_c));
    std::printf("SGP4 through the class: %s the C calls\n",
                y.ok() && bits(y->data(), r_c, sizeof r_c) ? "bit for bit with" : "NOT as");
  }
  const auto omm = myb::omm_from_tle(*tle);
  CHECK(omm.ok());
  if (omm) {
    const auto again = myb::tle_from_omm(*omm);
    CHECK(again.ok() && again->mean_motion == tle->mean_motion);
  }
  const auto broken = myb::tle_parse(line1, "2 garbage");
  CHECK(!broken.ok() && broken.code() == MYB_E_FORMAT && !broken.message().empty());
}
#endif

} // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  test_low_orbit();
  test_outcomes();
  test_callables();
  test_time_and_ephemeris();
  test_messages();
  test_nbody();
  test_system();
  test_estimation();
  test_cr3bp();
  test_targeting();
  test_observables();
  test_files();
  test_array_forms();
  test_contact_and_breaks();
#if defined(MYB_CPP_TEST_SGP4)
  test_two_line_elements();
#endif
  std::printf("test_cpp_binding: %ld checks, %ld failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
