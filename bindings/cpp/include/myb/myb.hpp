// SPDX-License-Identifier: BSD-3-Clause
// myb/myb.hpp — myb for C++17.
//
// The C headers (myb/myb.h) are C++ already: every function and structure of the library is
// usable from C++ as it is. This header adds what C cannot give:
//   - ownership: the objects of the library are classes that destroy what they hold; a force
//     model keeps alive the environment, the fields and the term structures it points to;
//   - outcomes: a function that can fail returns myb::result<T> — the status and, when it is
//     MYB_OK, the value. Asking a failed result for its value is a broken contract and stops
//     the program with a message, as MYB_CHECK does in the library;
//   - callables: terms of a model and events take lambdas.
// There is no exception handling here: nothing is thrown and nothing is caught, and the
// header compiles with exceptions switched off. A callable handed to the library reports a
// failure through its status; it must not throw — it is called from C, through a noexcept
// function, so an exception leaving it ends the program (std::terminate).
// It wraps the path most programs take — kernels and time, ephemerides, the environment and
// the force model, propagation, events, the text messages. For everything else call the C
// functions on what get() returns: the two mix freely.
//
// No arithmetic is done here: a result obtained through this header has the bits of the same
// calls made in C.
#ifndef MYB_MYB_HPP_
#define MYB_MYB_HPP_

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "myb/myb.h"

namespace myb {

using status = ::myb_status;
using vec3 = std::array<double, 3>;
using state = std::array<double, 6>; // position (m), velocity (m/s)

inline const char* status_name(status code) noexcept {
  return ::myb_status_name(code);
}

namespace detail {

[[noreturn]] inline void broken_contract(const char* what, status code) noexcept {
  std::fprintf(stderr, "myb: %s (%s)\n", what, ::myb_status_name(code));
  std::abort();
}

} // namespace detail

// ---- Outcomes ---------------------------------------------------------------------------------

// The outcome of a call: a status, the message of the library if it gave one (parsers name
// the line), and the value when the status is MYB_OK.
template <class T>
class [[nodiscard]] result {
 public:
  result(T value) : value_(std::move(value)) {} // NOLINT: a value is a successful outcome
  result(status code, std::string message = {}) // NOLINT: a status is a failed one
      : status_(code), message_(std::move(message)) {
    if (code == MYB_OK) {
      detail::broken_contract("a successful result needs its value", code);
    }
  }

  bool ok() const noexcept { return status_ == MYB_OK; }
  explicit operator bool() const noexcept { return ok(); }
  status code() const noexcept { return status_; }
  const std::string& message() const noexcept { return message_; }

  T& value() & {
    require();
    return *value_;
  }
  const T& value() const& {
    require();
    return *value_;
  }
  T&& value() && {
    require();
    return std::move(*value_);
  }
  T& operator*() & { return value(); }
  const T& operator*() const& { return value(); }
  T&& operator*() && { return std::move(*this).value(); }
  T* operator->() { return &value(); }
  const T* operator->() const { return &value(); }
  T value_or(T fallback) const& { return ok() ? *value_ : std::move(fallback); }

 private:
  void require() const noexcept {
    if (!ok()) {
      detail::broken_contract("the value of a failed result was asked for", status_);
    }
  }

  status status_ = MYB_OK;
  std::string message_;
  std::optional<T> value_;
};

template <>
class [[nodiscard]] result<void> {
 public:
  result() = default;
  result(status code, std::string message = {}) // NOLINT
      : status_(code), message_(std::move(message)) {}

  bool ok() const noexcept { return status_ == MYB_OK; }
  explicit operator bool() const noexcept { return ok(); }
  status code() const noexcept { return status_; }
  const std::string& message() const noexcept { return message_; }

 private:
  status status_ = MYB_OK;
  std::string message_;
};

// A result that owns its object through a shared pointer: what a force model and an orbit
// take, so that they can keep it alive.
template <class T>
result<std::shared_ptr<T>> share(result<T>&& made) {
  if (!made.ok()) {
    return {made.code(), made.message()};
  }
  return std::make_shared<T>(std::move(made).value());
}

namespace detail {

inline std::string describe(const ::myb_diag& diag) {
  std::string text = diag.message;
  if (diag.line > 0) {
    text = "line " + std::to_string(diag.line) + ": " + text;
  }
  return text;
}

inline result<void> outcome(status code, const ::myb_diag& diag) {
  if (code == MYB_OK) {
    return {};
  }
  return {code, describe(diag)};
}

inline result<void> outcome(status code) {
  if (code == MYB_OK) {
    return {};
  }
  return {code};
}

struct kpool_deleter {
  void operator()(::myb_kpool* p) const noexcept { ::myb_kpool_destroy(p); }
};
struct timesys_deleter {
  void operator()(::myb_timesys* p) const noexcept { ::myb_timesys_destroy(p); }
};
struct spk_deleter {
  void operator()(::myb_spk* p) const noexcept { ::myb_spk_destroy(p); }
};
struct field_deleter {
  void operator()(::myb_gravity_field* p) const noexcept { ::myb_gravity_field_destroy(p); }
};
struct work_deleter {
  void operator()(::myb_gravity_work* p) const noexcept { ::myb_gravity_work_destroy(p); }
};
struct env_deleter {
  void operator()(::myb_env* p) const noexcept { ::myb_env_destroy(p); }
};
struct model_deleter {
  void operator()(::myb_force_model* p) const noexcept { ::myb_force_model_destroy(p); }
};
struct orbit_deleter {
  void operator()(::myb_orbit* p) const noexcept { ::myb_orbit_destroy(p); }
};
struct oem_deleter {
  void operator()(::myb_oem* p) const noexcept { ::myb_oem_destroy(p); }
};

// A text of the library's writers: they tell the length when the buffer is too small.
template <class Write>
result<std::string> text(Write&& write) {
  std::size_t length = 0;
  status code = write(static_cast<char*>(nullptr), std::size_t{0}, &length);
  if (code != MYB_OK && code != MYB_E_OVERFLOW) {
    return code;
  }
  std::string out(length + 1, '\0');
  code = write(out.data(), out.size(), &length);
  if (code != MYB_OK) {
    return code;
  }
  out.resize(length);
  return out;
}

} // namespace detail

// ---- Two bodies -------------------------------------------------------------------------------

inline state to_state(const ::myb_pv& pv) noexcept {
  return {pv.r[0], pv.r[1], pv.r[2], pv.v[0], pv.v[1], pv.v[2]};
}

inline ::myb_pv to_pv(const state& y) noexcept {
  return {{y[0], y[1], y[2]}, {y[3], y[4], y[5]}};
}

inline result<state> state_from_elements(double mu, const ::myb_elements& elements) {
  ::myb_pv pv{};
  const status code = ::myb_pv_from_elements(mu, &elements, &pv);
  if (code != MYB_OK) {
    return code;
  }
  return to_state(pv);
}

inline result<::myb_elements> elements_from_state(double mu, const state& y) {
  const ::myb_pv pv = to_pv(y);
  ::myb_elements elements{};
  const status code = ::myb_elements_from_pv(mu, &pv, &elements);
  if (code != MYB_OK) {
    return code;
  }
  return elements;
}

// The state a time dt later on the orbit of two bodies.
inline result<state> kepler_propagate(double mu, const state& y0, double dt) {
  const ::myb_pv pv0 = to_pv(y0);
  ::myb_pv pv{};
  const status code = ::myb_kepler_propagate(mu, &pv0, dt, &pv);
  if (code != MYB_OK) {
    return code;
  }
  return to_state(pv);
}

// ---- Kernels and time -------------------------------------------------------------------------

// A pool of text kernels (myb/io.h).
class kernel_pool {
 public:
  static result<kernel_pool> create() {
    ::myb_kpool* pool = nullptr;
    const status code = ::myb_kpool_create(nullptr, &pool);
    if (code != MYB_OK) {
      return code;
    }
    return kernel_pool(pool);
  }

  result<void> load_file(const std::string& path) {
    ::myb_diag diag{};
    return detail::outcome(::myb_kpool_load_file(pool_.get(), path.c_str(), &diag), diag);
  }

  ::myb_kpool* get() const noexcept { return pool_.get(); }

 private:
  explicit kernel_pool(::myb_kpool* pool) noexcept : pool_(pool) {}
  std::unique_ptr<::myb_kpool, detail::kpool_deleter> pool_;
};

// The time scales of a pool (myb/time.h). It keeps nothing of the pool.
class time_system {
 public:
  static result<time_system> create(const kernel_pool& pool) {
    ::myb_timesys* system = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_timesys_create(pool.get(), nullptr, &system, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return time_system(system);
  }

  result<::myb_scale> scale(const std::string& name) const {
    ::myb_scale found{};
    const status code = ::myb_timesys_scale(system_.get(), name.c_str(), &found);
    if (code != MYB_OK) {
      return code;
    }
    return found;
  }

  result<::myb_time> convert(::myb_time t, ::myb_scale from, ::myb_scale to) const {
    ::myb_time out{};
    const status code = ::myb_timesys_convert(system_.get(), t, from, to, &out);
    if (code != MYB_OK) {
      return code;
    }
    return out;
  }

  result<::myb_time> from_utc(const ::myb_calendar& utc, ::myb_scale to) const {
    ::myb_time out{};
    const status code = ::myb_timesys_from_utc(system_.get(), &utc, to, &out);
    if (code != MYB_OK) {
      return code;
    }
    return out;
  }

  result<::myb_calendar> to_utc(::myb_time t, ::myb_scale from) const {
    ::myb_calendar utc{};
    const status code = ::myb_timesys_to_utc(system_.get(), t, from, &utc);
    if (code != MYB_OK) {
      return code;
    }
    return utc;
  }

  ::myb_timesys* get() const noexcept { return system_.get(); }

 private:
  explicit time_system(::myb_timesys* system) noexcept : system_(system) {}
  std::unique_ptr<::myb_timesys, detail::timesys_deleter> system_;
};

// ---- Ephemerides ------------------------------------------------------------------------------

// A set of SPK files (myb/kin.h).
class ephemeris {
 public:
  static result<ephemeris> create() {
    ::myb_spk* spk = nullptr;
    const status code = ::myb_spk_create(nullptr, &spk);
    if (code != MYB_OK) {
      return code;
    }
    return ephemeris(spk);
  }

  result<void> load_file(const std::string& path) {
    ::myb_diag diag{};
    return detail::outcome(::myb_spk_load_file(spk_.get(), path.c_str(), &diag), diag);
  }

  // The state of target relative to center at t on the time scale of the files.
  result<state> at(std::int32_t target, std::int32_t center, ::myb_time t) const {
    state y{};
    const status code = ::myb_spk_state(spk_.get(), target, center, t, y.data(), nullptr);
    if (code != MYB_OK) {
      return code;
    }
    return y;
  }

  // A body of the set as a provider of an environment: state = myb_spk_body_state, context =
  // the pointer returned, which lives as long as this object. frame 0 accepts any frame.
  ::myb_spk_body* body(std::int32_t target, std::int32_t center, std::int32_t frame = 0) {
    bodies_.push_back(std::make_unique<::myb_spk_body>(
        ::myb_spk_body{spk_.get(), target, center, frame}));
    return bodies_.back().get();
  }

  ::myb_spk* get() const noexcept { return spk_.get(); }

 private:
  explicit ephemeris(::myb_spk* spk) noexcept : spk_(spk) {}
  std::unique_ptr<::myb_spk, detail::spk_deleter> spk_;
  std::vector<std::unique_ptr<::myb_spk_body>> bodies_;
};

// ---- Fields and the environment ----------------------------------------------------------------

// The field of a body in spherical harmonics (myb/grav.h).
class gravity_field {
 public:
  // c, s: fully normalised coefficients at the index n (n + 1) / 2 + m, to the degree given.
  static result<gravity_field> create(double gm, double radius, int degree,
                                      const std::vector<double>& c,
                                      const std::vector<double>& s) {
    const std::size_t needed =
        degree < 0 ? 0 : static_cast<std::size_t>(degree + 1) * static_cast<std::size_t>(degree + 2) / 2;
    if (degree < 0 || c.size() < needed || s.size() < needed) {
      return MYB_E_DOMAIN;
    }
    ::myb_gravity_field* field = nullptr;
    const status code =
        ::myb_gravity_field_create(gm, radius, degree, c.data(), s.data(), nullptr, &field);
    if (code != MYB_OK) {
      return code;
    }
    return gravity_field(field);
  }

  // An ICGEM file (gfc), to max_degree (a negative value: all of it).
  static result<gravity_field> load_icgem_file(const std::string& path, int max_degree) {
    ::myb_gravity_field* field = nullptr;
    ::myb_diag diag{};
    const status code =
        ::myb_gravity_field_load_icgem_file(path.c_str(), max_degree, nullptr, &field, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return gravity_field(field);
  }

  int degree() const noexcept { return ::myb_gravity_field_degree(field_.get()); }
  double gm() const noexcept { return ::myb_gravity_field_gm(field_.get()); }
  double radius() const noexcept { return ::myb_gravity_field_radius(field_.get()); }
  ::myb_gravity_field* get() const noexcept { return field_.get(); }

 private:
  explicit gravity_field(::myb_gravity_field* field) noexcept : field_(field) {}
  std::unique_ptr<::myb_gravity_field, detail::field_deleter> field_;
};

// The bodies of a scenario (myb/env.h): what the terms of a force model ask for.
class environment {
 public:
  static result<environment> create(const ::myb_env_options& options = {}) {
    ::myb_env* env = nullptr;
    const status code = ::myb_env_create(&options, nullptr, &env);
    if (code != MYB_OK) {
      return code;
    }
    return environment(env);
  }

  // The index of the body in the environment. The contexts of its providers must live as
  // long as the environment: hand their owners to keep().
  result<std::size_t> add_body(const ::myb_env_body& body) {
    std::size_t index = 0;
    const status code = ::myb_env_add_body(env_.get(), &body, &index);
    if (code != MYB_OK) {
      return code;
    }
    return index;
  }

  // Ties the life of something the environment points to (an ephemeris, a context of a
  // provider) to the environment.
  void keep(std::shared_ptr<void> owner) { kept_.push_back(std::move(owner)); }

  result<std::size_t> find(std::int64_t id) const {
    std::size_t index = 0;
    const status code = ::myb_env_find(env_.get(), id, &index);
    if (code != MYB_OK) {
      return code;
    }
    return index;
  }

  // The state of a body relative to the origin at the model time t.
  result<state> at(double t, std::size_t body) {
    state y{};
    const status code = ::myb_env_state(env_.get(), t, body, y.data());
    if (code != MYB_OK) {
      return code;
    }
    return y;
  }

  std::size_t body_count() const noexcept { return ::myb_env_body_count(env_.get()); }
  double gm(std::size_t body) const noexcept { return ::myb_env_gm(env_.get(), body); }
  ::myb_env* get() const noexcept { return env_.get(); }

 private:
  explicit environment(::myb_env* env) noexcept : env_(env) {}
  std::unique_ptr<::myb_env, detail::env_deleter> env_;
  std::vector<std::shared_ptr<void>> kept_;
};

// ---- The force model ----------------------------------------------------------------------------

// A force on the object, as a callable: it writes into a (zeros when it is called) the
// acceleration of the term at the model time t and the state (x, v) relative to the origin of
// the model, and returns MYB_OK — or the status that stops the evaluation. It must not throw.
using term_function = std::function<status(double t, const vec3& x, const vec3& v, vec3& a)>;

// The forces on one object (myb/model.h) with the ready-made terms of myb/env.h. The model
// owns the structures of its terms and keeps what they point to alive.
class force_model {
 public:
  static result<force_model> create() {
    ::myb_force_model* model = nullptr;
    const status code = ::myb_force_model_create(nullptr, &model);
    if (code != MYB_OK) {
      return code;
    }
    return force_model(model);
  }

  // Stores a value for the life of the model: the reference stays valid (and at the same
  // address) until the model is destroyed. For the contexts of terms added with add_raw.
  template <class T>
  T& keep(T value) {
    auto holder = std::make_shared<T>(std::move(value));
    T& reference = *holder;
    kept_.push_back(std::move(holder));
    return reference;
  }

  // Any term of the C interface; the context must live as long as the model (see keep).
  result<std::size_t> add_raw(const std::string& name, ::myb_term_function function,
                              ::myb_term_jacobian jacobian, void* context) {
    std::size_t index = 0;
    const status code =
        ::myb_force_model_add(model_.get(), name.c_str(), function, jacobian, context, &index);
    if (code != MYB_OK) {
      return code;
    }
    return index;
  }

  // The attraction of a point mass at the origin.
  result<std::size_t> add_point_mass(const std::string& name, double gm) {
    auto& term = keep(::myb_point_mass_term{gm});
    return add_raw(name, ::myb_term_point_mass, ::myb_term_point_mass_jacobian, &term);
  }

  // The attraction of a body of the environment away from the origin, with its indirect part.
  result<std::size_t> add_third_body(const std::string& name, std::shared_ptr<environment> env,
                                     std::size_t body) {
    auto& term = keep(::myb_third_body_term{env->get(), body});
    kept_.push_back(std::move(env));
    return add_raw(name, ::myb_term_third_body, ::myb_term_third_body_jacobian, &term);
  }

  // The field of a body of the environment in harmonics of the degrees n_min .. n_max
  // (n_max < 0: all the field has).
  result<std::size_t> add_field(const std::string& name, std::shared_ptr<environment> env,
                                std::size_t body, std::shared_ptr<gravity_field> field,
                                int n_min = 0, int n_max = -1) {
    ::myb_gravity_work* work = nullptr;
    const status code = ::myb_gravity_work_create(field->degree(), nullptr, &work);
    if (code != MYB_OK) {
      return code;
    }
    kept_.push_back(std::shared_ptr<::myb_gravity_work>(work, detail::work_deleter{}));
    ::myb_field_term made{};
    made.env = env->get();
    made.body = body;
    made.field = field->get();
    made.work = work;
    made.n_min = n_min;
    made.n_max = n_max;
    auto& term = keep(made);
    kept_.push_back(std::move(env));
    kept_.push_back(std::move(field));
    return add_raw(name, ::myb_term_field, ::myb_term_field_jacobian, &term);
  }

  // The pressure of the light of a body of the environment on the object, in the shadows of
  // the occulting bodies given.
  result<std::size_t> add_radiation(const std::string& name, std::shared_ptr<environment> env,
                                    std::size_t source, const ::myb_object& object, double c,
                                    const std::vector<std::size_t>& occulters = {}) {
    if (occulters.size() > MYB_ENV_MAX_OCCULTERS) {
      return MYB_E_DOMAIN;
    }
    ::myb_radiation_term made{};
    made.env = env->get();
    made.source = source;
    made.object = &keep(object);
    made.c = c;
    for (std::size_t k = 0; k < occulters.size(); ++k) {
      made.occulters[k] = occulters[k];
    }
    made.occulter_count = occulters.size();
    auto& term = keep(made);
    kept_.push_back(std::move(env));
    return add_raw(name, ::myb_term_radiation, ::myb_term_radiation_jacobian, &term);
  }

  // The drag of a layered exponential atmosphere that turns with a body of the environment.
  result<std::size_t> add_drag(const std::string& name, std::shared_ptr<environment> env,
                               std::size_t body, const ::myb_exponential_atmosphere& air,
                               const ::myb_object& object) {
    ::myb_drag_term made{};
    made.env = env->get();
    made.body = body;
    made.density = ::myb_density_exponential;
    made.density_context = &keep(air);
    made.object = &keep(object);
    auto& term = keep(made);
    kept_.push_back(std::move(env));
    return add_raw(name, ::myb_term_drag, ::myb_term_drag_jacobian, &term);
  }

  // The post-Newtonian correction of a central body of the environment (Schwarzschild).
  result<std::size_t> add_schwarzschild(const std::string& name,
                                        std::shared_ptr<environment> env, std::size_t central,
                                        double c, double beta = 1.0, double gamma = 1.0) {
    ::myb_relativity_term made{};
    made.env = env->get();
    made.ppn = ::myb_ppn{c, beta, gamma};
    made.central = central;
    made.external = central;
    auto& term = keep(made);
    kept_.push_back(std::move(env));
    return add_raw(name, ::myb_term_schwarzschild, ::myb_term_schwarzschild_jacobian, &term);
  }

  // A term of the caller's own; its derivatives are taken by differences of the model.
  result<std::size_t> add_term(const std::string& name, term_function function) {
    auto& stored = keep(std::move(function));
    return add_raw(name, &force_model::call_term, nullptr, &stored);
  }

  // The acceleration of the object at the model time t and the state (x, v).
  result<vec3> acceleration(double t, const vec3& x, const vec3& v) {
    vec3 a{};
    const status code =
        ::myb_force_model_acceleration(model_.get(), t, x.data(), v.data(), a.data());
    if (code != MYB_OK) {
      return code;
    }
    return a;
  }

  std::size_t term_count() const noexcept { return ::myb_force_model_count(model_.get()); }
  std::string term_name(std::size_t index) const {
    return ::myb_force_model_name(model_.get(), index);
  }
  bool enabled(std::size_t index) const noexcept {
    return ::myb_force_model_enabled(model_.get(), index) != 0;
  }
  void set_enabled(std::size_t index, bool on) noexcept {
    ::myb_force_model_set_enabled(model_.get(), index, on ? 1 : 0);
  }
  ::myb_force_model* get() const noexcept { return model_.get(); }

 private:
  explicit force_model(::myb_force_model* model) noexcept : model_(model) {}

  static status call_term(void* context, double t, const double x[3], const double v[3],
                          double a[3]) noexcept {
    const auto& function = *static_cast<const term_function*>(context);
    const vec3 position{x[0], x[1], x[2]};
    const vec3 velocity{v[0], v[1], v[2]};
    vec3 acceleration{};
    const status code = function(t, position, velocity, acceleration);
    if (code != MYB_OK) {
      return code;
    }
    a[0] = acceleration[0];
    a[1] = acceleration[1];
    a[2] = acceleration[2];
    return MYB_OK;
  }

  // kept_ is declared first and so destroyed last: the model goes before what it points to
  std::vector<std::shared_ptr<void>> kept_;
  std::unique_ptr<::myb_force_model, detail::model_deleter> model_;
};

// ---- Propagation --------------------------------------------------------------------------------

// A function of the state whose zeros are looked for: g(t, y). It must not throw.
using event_function = std::function<double(double t, const state& y)>;

// One object under a force model (myb/orbit.h): states at the times asked for, events.
class orbit {
 public:
  struct crossing {
    double t;
    state y;
  };
  struct partials {
    state y;
    std::array<double, 36> phi; // the state transition matrix from the start, by rows
  };

  static result<orbit> create(std::shared_ptr<force_model> model,
                              const ::myb_orbit_options& options = {}) {
    ::myb_orbit* made = nullptr;
    const status code = ::myb_orbit_create(model->get(), &options, nullptr, &made);
    if (code != MYB_OK) {
      return code;
    }
    return orbit(made, std::move(model));
  }

  // The object is at y0 at t0 and will be asked about up to t_end (before t0: backwards).
  result<void> start(double t0, const state& y0, double t_end) {
    return detail::outcome(::myb_orbit_start(orbit_.get(), t0, y0.data(), t_end));
  }

  result<state> at(double t) {
    state y{};
    const status code = ::myb_orbit_state(orbit_.get(), t, y.data());
    if (code != MYB_OK) {
      return code;
    }
    return y;
  }

  // The states at times in the order of the propagation.
  result<std::vector<state>> at(const std::vector<double>& times) {
    std::vector<state> states(times.size());
    static_assert(sizeof(state) == 6 * sizeof(double), "a state is six doubles in a row");
    const status code = ::myb_orbit_states(orbit_.get(), times.data(), times.size(),
                                           states.empty() ? nullptr : states.front().data());
    if (code != MYB_OK) {
      return code;
    }
    return states;
  }

  // The next zero of g up to t_limit: direction > 0 for a rising g, < 0 for a falling one,
  // 0 for both. No value: none before t_limit.
  result<std::optional<crossing>> next_event(const event_function& g, double t_limit,
                                             int direction = 0) {
    int found = 0;
    crossing hit{};
    const status code =
        ::myb_orbit_event(orbit_.get(), &orbit::call_event,
                          const_cast<void*>(static_cast<const void*>(&g)), t_limit, direction,
                          &found, &hit.t, hit.y.data());
    if (code != MYB_OK) {
      return code;
    }
    if (!found) {
      return std::optional<crossing>{};
    }
    return std::optional<crossing>{hit};
  }

  // An impulse: the velocity changes by dv where the object is now.
  result<void> kick(const vec3& dv) {
    return detail::outcome(::myb_orbit_kick(orbit_.get(), dv.data()));
  }

  // The state and the transition matrix (options.variational must be set).
  result<partials> partials_at(double t) {
    partials out{};
    const status code =
        ::myb_orbit_partials(orbit_.get(), t, out.y.data(), out.phi.data(), nullptr);
    if (code != MYB_OK) {
      return code;
    }
    return out;
  }

  double time() const noexcept { return ::myb_orbit_time(orbit_.get()); }
  std::int64_t steps() const noexcept { return ::myb_orbit_steps(orbit_.get()); }
  std::int64_t breaks() const noexcept { return ::myb_orbit_breaks(orbit_.get()); }
  ::myb_orbit* get() const noexcept { return orbit_.get(); }

 private:
  orbit(::myb_orbit* made, std::shared_ptr<force_model> model) noexcept
      : model_(std::move(model)), orbit_(made) {}

  static double call_event(void* context, double t, const double* y) noexcept {
    const auto& g = *static_cast<const event_function*>(context);
    const state now{y[0], y[1], y[2], y[3], y[4], y[5]};
    return g(t, now);
  }

  std::shared_ptr<force_model> model_; // outlives the propagator
  std::unique_ptr<::myb_orbit, detail::orbit_deleter> orbit_;
};

// ---- Messages -----------------------------------------------------------------------------------

// An orbit ephemeris message as text, in the notation given.
inline result<std::string> to_text(const ::myb_oem_data& data,
                                   ::myb_odm_notation notation = MYB_ODM_KVN) {
  return detail::text([&](char* buffer, std::size_t size, std::size_t* length) {
    return ::myb_oem_write_as(&data, notation, buffer, size, length);
  });
}

inline result<std::string> to_text(const ::myb_opm& opm,
                                   ::myb_odm_notation notation = MYB_ODM_KVN) {
  return detail::text([&](char* buffer, std::size_t size, std::size_t* length) {
    return ::myb_opm_write_as(&opm, notation, buffer, size, length);
  });
}

inline result<std::string> to_text(const ::myb_omm& omm,
                                   ::myb_odm_notation notation = MYB_ODM_KVN) {
  return detail::text([&](char* buffer, std::size_t size, std::size_t* length) {
    return ::myb_omm_write_as(&omm, notation, buffer, size, length);
  });
}

inline result<std::string> to_text(const ::myb_ocm_data& data,
                                   ::myb_odm_notation notation = MYB_ODM_KVN) {
  return detail::text([&](char* buffer, std::size_t size, std::size_t* length) {
    return ::myb_ocm_write(&data, notation, buffer, size, length);
  });
}

// An orbit ephemeris message read from text in either notation (myb/ccsds.h).
class ephemeris_message {
 public:
  struct point {
    ::myb_calendar epoch;
    state y;
  };

  static result<ephemeris_message> parse(const std::string& text) {
    ::myb_oem* oem = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_oem_load(text.data(), text.size(), nullptr, &oem, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return ephemeris_message(oem);
  }

  std::size_t segment_count() const noexcept { return ::myb_oem_segment_count(oem_.get()); }
  std::size_t state_count(std::size_t segment) const noexcept {
    return ::myb_oem_state_count(oem_.get(), segment);
  }

  // A keyword of the metadata of a segment; no value if the message does not have it.
  std::optional<std::string> meta(std::size_t segment, const std::string& key) const {
    const char* value = ::myb_oem_meta(oem_.get(), segment, key.c_str());
    return value != nullptr ? std::optional<std::string>{value} : std::nullopt;
  }

  point at(std::size_t segment, std::size_t index) const noexcept {
    point out{};
    (void)::myb_oem_state(oem_.get(), segment, index, &out.epoch, out.y.data(), nullptr);
    return out;
  }

  ::myb_oem* get() const noexcept { return oem_.get(); }

 private:
  explicit ephemeris_message(::myb_oem* oem) noexcept : oem_(oem) {}
  std::unique_ptr<::myb_oem, detail::oem_deleter> oem_;
};

} // namespace myb

#endif // MYB_MYB_HPP_
