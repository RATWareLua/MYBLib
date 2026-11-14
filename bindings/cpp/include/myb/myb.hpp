// SPDX-License-Identifier: BSD-3-Clause
// myb/myb.hpp — myb for C++17.
//
// The C headers (myb/myb.h) are C++ already: every function and structure of the library is
// usable from C++ as it is. This header adds what C cannot give:
//   - ownership: the objects of the library are classes that destroy what they hold; an object
//     keeps alive what it points to — a force model its environment, fields and the structures
//     of its terms, a propagator its model, a system of bodies its providers and fields;
//   - outcomes: a function that can fail returns myb::result<T> — the status and, when it is
//     MYB_OK, the value. Asking a failed result for its value is a broken contract and stops
//     the program with a message, as MYB_CHECK does in the library. An array form that stops at
//     an element names it in the message ("element 3");
//   - callables: terms of a model, providers of states and orientations, measurements, actions
//     on bodies, torques, events and the functions of the solvers take lambdas.
// There is no exception handling here: nothing is thrown and nothing is caught, and the
// header compiles with exceptions switched off. A callable handed to the library reports a
// failure through its status; it must not throw — it is called from C, through a noexcept
// function, so an exception leaving it ends the program (std::terminate).
//
// The areas: kernels and time, frames and ephemerides, the environment and the force model,
// propagation and events, contact, N bodies and systems of bodies, estimation, uncertainty, the
// circular restricted three-body problem, transfers and optimal control, low thrust,
// conjunctions, attitude, observables, the readers and writers of files (with the two-line
// elements of ext/sgp4 when its header is included before this one), and the array forms of
// time, ephemerides, rotations and two bodies over vectors. Structures of the C interface that
// are plain data (elements, options, problems, results) are used as they are. For what is not
// here call the C functions on what get() returns: the two mix freely.
//
// No arithmetic is done here: a result obtained through this header has the bits of the same
// calls made in C. (Sizes of buffers are counted here; numbers are not.)
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
#include <type_traits>
#include <utility>
#include <vector>

#include "myb/myb.h"

namespace myb {

using status = ::myb_status;
using vec3 = std::array<double, 3>;
using state = std::array<double, 6>; // position (m), velocity (m/s)
using quat = std::array<double, 4>;  // scalar first (myb/attitude.h)
using mat3 = std::array<double, 9>;  // by rows
using mat6 = std::array<double, 36>; // by rows

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

// The value of a call that returned `code` (taken after the call: the status first, always).
template <class T>
result<T> give(status code, T value) {
  if (code != MYB_OK) {
    return code;
  }
  return value;
}

template <class T>
result<T> give(status code, const ::myb_diag& diag, T value) {
  if (code != MYB_OK) {
    return {code, describe(diag)};
  }
  return value;
}

// The value of an array form, which stops at the element that fails.
template <class T>
result<T> give_many(status code, std::size_t failed, T value) {
  if (code != MYB_OK) {
    return {code, "element " + std::to_string(failed)};
  }
  return value;
}

template <std::size_t N>
std::array<double, N> to_array(const double* values) noexcept {
  std::array<double, N> out{};
  for (std::size_t k = 0; k < N; ++k) {
    out[k] = values[k];
  }
  return out;
}

template <std::size_t N>
void copy_out(const std::array<double, N>& in, double* out) noexcept {
  for (std::size_t k = 0; k < N; ++k) {
    out[k] = in[k];
  }
}

// A vector of fixed-size arrays as the contiguous doubles the C functions take.
template <std::size_t N>
double* flat(std::vector<std::array<double, N>>& v) noexcept {
  static_assert(sizeof(std::array<double, N>) == N * sizeof(double), "arrays of doubles in a row");
  return v.empty() ? nullptr : v.front().data();
}
template <std::size_t N>
const double* flat(const std::vector<std::array<double, N>>& v) noexcept {
  static_assert(sizeof(std::array<double, N>) == N * sizeof(double), "arrays of doubles in a row");
  return v.empty() ? nullptr : v.front().data();
}

template <class T>
const T* data_or_null(const std::vector<T>& v) noexcept {
  return v.empty() ? nullptr : v.data();
}

// An object of the library with what it points to: the object goes first, then what it kept.
// Shared, so that a view into a larger object (the point masses of a system) keeps it alive.
template <class T, class Deleter>
struct core {
  core() = default;
  core(const core&) = delete;
  core& operator=(const core&) = delete;
  ~core() { Deleter{}(object); }
  std::vector<std::shared_ptr<const void>> kept;
  T* object = nullptr;
};

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
struct ocm_deleter {
  void operator()(::myb_ocm* p) const noexcept { ::myb_ocm_destroy(p); }
};
struct cip_deleter {
  void operator()(::myb_cip* p) const noexcept { ::myb_cip_destroy(p); }
};
struct eop_deleter {
  void operator()(::myb_eop* p) const noexcept { ::myb_eop_destroy(p); }
};
struct equinox_deleter {
  void operator()(::myb_equinox* p) const noexcept { ::myb_equinox_destroy(p); }
};
struct wgccre_deleter {
  void operator()(::myb_wgccre* p) const noexcept { ::myb_wgccre_destroy(p); }
};
struct pck_deleter {
  void operator()(::myb_pck* p) const noexcept { ::myb_pck_destroy(p); }
};
struct graph_deleter {
  void operator()(::myb_frame_graph* p) const noexcept { ::myb_frame_graph_destroy(p); }
};
struct pool_deleter {
  void operator()(::myb_pool* p) const noexcept { ::myb_pool_destroy(p); }
};
struct nbody_deleter {
  void operator()(::myb_nbody* p) const noexcept { ::myb_nbody_destroy(p); }
};
struct wh_deleter {
  void operator()(::myb_wh* p) const noexcept { ::myb_wh_destroy(p); }
};
struct system_deleter {
  void operator()(::myb_system* p) const noexcept { ::myb_system_destroy(p); }
};
struct ias15_deleter {
  void operator()(::myb_ias15* p) const noexcept { ::myb_ias15_destroy(p); }
};
struct filter_deleter {
  void operator()(::myb_filter* p) const noexcept { ::myb_filter_destroy(p); }
};
struct variational_deleter {
  void operator()(::myb_variational* p) const noexcept { ::myb_variational_destroy(p); }
};
struct table_deleter {
  void operator()(::myb_table* p) const noexcept { ::myb_table_destroy(p); }
};
struct spk_writer_deleter {
  void operator()(::myb_spk_writer* p) const noexcept { ::myb_spk_writer_destroy(p); }
};
struct sp3_deleter {
  void operator()(::myb_sp3* p) const noexcept { ::myb_sp3_destroy(p); }
};
struct sclk_deleter {
  void operator()(::myb_sclk* p) const noexcept { ::myb_sclk_destroy(p); }
};
struct ck_deleter {
  void operator()(::myb_ck* p) const noexcept { ::myb_ck_destroy(p); }
};
struct ck_writer_deleter {
  void operator()(::myb_ck_writer* p) const noexcept { ::myb_ck_writer_destroy(p); }
};

// A text of the library's writers: they tell the length when the buffer is too small. The
// first buffer is a real one (some writers take no null), large enough for the short texts.
template <class Write>
result<std::string> text(Write&& write) {
  std::size_t length = 0;
  char first[128] = {};
  status code = write(first, sizeof first, &length);
  if (code == MYB_OK) {
    return std::string(first, length);
  }
  if (code != MYB_E_OVERFLOW) {
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

// The bytes of a writer of binary files: its size, then the file.
template <class Writer>
result<std::vector<unsigned char>> bytes(std::size_t size, Writer&& write) {
  std::vector<unsigned char> out(size);
  const status code = write(out.empty() ? nullptr : out.data(), out.size());
  if (code != MYB_OK) {
    return code;
  }
  return out;
}

inline void* context_of(const void* callable) noexcept {
  return const_cast<void*>(callable);
}

} // namespace detail

// ---- Providers ----------------------------------------------------------------------------------

// The state of a point at an instant, for what takes a provider (myb_state_function of
// myb/kin.h): it writes the position (m) and the velocity (m/s) into y and returns MYB_OK — or
// the status that stops the caller. It must not throw.
using state_function = std::function<status(::myb_time t, state& y)>;

// The orientation of a body at an instant (myb_orientation_function of myb/frames.h): it fills
// the rotation from the inertial axes to the body (with its rate) and, for a body whose rotation
// is observed, its rotation time and pole; it returns MYB_OK or the status that stops the caller.
using orientation_function = std::function<status(::myb_time t, ::myb_orientation& orientation)>;

namespace detail {
// What C calls for them; noexcept, as call_term below.
inline status call_state(void* context, ::myb_time t, double y[6]) noexcept {
  const auto& function = *static_cast<const state_function*>(context);
  state out{};
  const status code = function(t, out);
  if (code != MYB_OK) {
    return code;
  }
  copy_out(out, y);
  return MYB_OK;
}
static_assert(noexcept(call_state(nullptr, ::myb_time{}, nullptr)),
              "a function called from C must not let an exception through");

inline status call_orientation(void* context, ::myb_time t,
                               ::myb_orientation* orientation) noexcept {
  const auto& function = *static_cast<const orientation_function*>(context);
  ::myb_orientation out{};
  const status code = function(t, out);
  if (code != MYB_OK) {
    return code;
  }
  *orientation = out;
  return MYB_OK;
}
static_assert(noexcept(call_orientation(nullptr, ::myb_time{}, nullptr)),
              "a function called from C must not let an exception through");
} // namespace detail

// A provider of states for the C structures and calls that take one (a function and its
// context): a callable of the caller's, or a provider of the C interface with its context and,
// optionally, what keeps that context alive. Copies share the callable.
class provider {
 public:
  template <class F, class = std::enable_if_t<std::is_constructible<state_function, F>::value &&
                                              !std::is_same<std::decay_t<F>, provider>::value>>
  provider(F&& function) // NOLINT: a callable is a provider
      : function_(&detail::call_state) {
    auto stored = std::make_shared<state_function>(std::forward<F>(function));
    context_ = stored.get();
    owner_ = std::move(stored);
  }
  provider(::myb_state_function function, void* context,
           std::shared_ptr<const void> owner = {}) noexcept
      : function_(function), context_(context), owner_(std::move(owner)) {}

  ::myb_state_function function() const noexcept { return function_; }
  void* context() const noexcept { return context_; }
  const std::shared_ptr<const void>& owner() const noexcept { return owner_; }

  // The state at t, as the users of the provider see it.
  result<state> at(::myb_time t) const {
    state y{};
    const status code = function_(context_, t, y.data());
    return detail::give(code, y);
  }

 private:
  ::myb_state_function function_;
  void* context_ = nullptr;
  std::shared_ptr<const void> owner_;
};

// The same for orientations (myb_orientation_function).
class orientation_provider {
 public:
  template <class F,
            class = std::enable_if_t<std::is_constructible<orientation_function, F>::value &&
                                     !std::is_same<std::decay_t<F>, orientation_provider>::value>>
  orientation_provider(F&& function) // NOLINT: a callable is a provider
      : function_(&detail::call_orientation) {
    auto stored = std::make_shared<orientation_function>(std::forward<F>(function));
    context_ = stored.get();
    owner_ = std::move(stored);
  }
  orientation_provider(::myb_orientation_function function, void* context,
                       std::shared_ptr<const void> owner = {}) noexcept
      : function_(function), context_(context), owner_(std::move(owner)) {}

  ::myb_orientation_function function() const noexcept { return function_; }
  void* context() const noexcept { return context_; }
  const std::shared_ptr<const void>& owner() const noexcept { return owner_; }

  result<::myb_orientation> at(::myb_time t) const {
    ::myb_orientation out{};
    const status code = function_(context_, t, &out);
    return detail::give(code, out);
  }

 private:
  ::myb_orientation_function function_;
  void* context_ = nullptr;
  std::shared_ptr<const void> owner_;
};

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

inline result<state> state_from_equinoctial(double mu, const ::myb_equinoctial& elements) {
  ::myb_pv pv{};
  const status code = ::myb_pv_from_equinoctial(mu, &elements, &pv);
  if (code != MYB_OK) {
    return code;
  }
  return to_state(pv);
}

inline result<::myb_equinoctial> equinoctial_from_state(double mu, const state& y) {
  const ::myb_pv pv = to_pv(y);
  ::myb_equinoctial elements{};
  const status code = ::myb_equinoctial_from_pv(mu, &pv, &elements);
  return detail::give(code, elements);
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

// The state a time dt later and its transition matrix.
struct kepler_transition {
  state y;
  mat6 stm; // d y_i / d y0_j, by rows
};

inline result<kepler_transition> kepler_stm(double mu, const state& y0, double dt) {
  const ::myb_pv pv0 = to_pv(y0);
  ::myb_pv pv{};
  kepler_transition out{};
  const status code = ::myb_kepler_stm(mu, &pv0, dt, &pv, out.stm.data());
  out.y = to_state(pv);
  return detail::give(code, out);
}

// The array forms (myb/kepler.h): element k of the output from element k of the input.
namespace detail {
inline std::vector<::myb_pv> to_pvs(const std::vector<state>& states) {
  std::vector<::myb_pv> pv(states.size());
  for (std::size_t k = 0; k < states.size(); ++k) {
    pv[k] = to_pv(states[k]);
  }
  return pv;
}
inline std::vector<state> to_states(const std::vector<::myb_pv>& pv) {
  std::vector<state> states(pv.size());
  for (std::size_t k = 0; k < pv.size(); ++k) {
    states[k] = to_state(pv[k]);
  }
  return states;
}
} // namespace detail

inline result<std::vector<::myb_elements>> elements_from_state(double mu,
                                                               const std::vector<state>& states) {
  const std::vector<::myb_pv> pv = detail::to_pvs(states);
  std::vector<::myb_elements> elements(states.size());
  std::size_t failed = 0;
  const status code = ::myb_elements_from_pv_many(mu, detail::data_or_null(pv), pv.size(),
                                                  elements.data(), &failed);
  return detail::give_many(code, failed, std::move(elements));
}

inline result<std::vector<state>> state_from_elements(double mu,
                                                      const std::vector<::myb_elements>& elements) {
  std::vector<::myb_pv> pv(elements.size());
  std::size_t failed = 0;
  const status code = ::myb_pv_from_elements_many(mu, detail::data_or_null(elements),
                                                  elements.size(), pv.data(), &failed);
  return detail::give_many(code, failed, detail::to_states(pv));
}

inline result<std::vector<::myb_equinoctial>> equinoctial_from_state(
    double mu, const std::vector<state>& states) {
  const std::vector<::myb_pv> pv = detail::to_pvs(states);
  std::vector<::myb_equinoctial> elements(states.size());
  std::size_t failed = 0;
  const status code = ::myb_equinoctial_from_pv_many(mu, detail::data_or_null(pv), pv.size(),
                                                     elements.data(), &failed);
  return detail::give_many(code, failed, std::move(elements));
}

inline result<std::vector<state>> state_from_equinoctial(
    double mu, const std::vector<::myb_equinoctial>& elements) {
  std::vector<::myb_pv> pv(elements.size());
  std::size_t failed = 0;
  const status code = ::myb_pv_from_equinoctial_many(mu, detail::data_or_null(elements),
                                                     elements.size(), pv.data(), &failed);
  return detail::give_many(code, failed, detail::to_states(pv));
}

// From one state to the time offsets dt[k].
inline result<std::vector<state>> kepler_propagate(double mu, const state& y0,
                                                   const std::vector<double>& dt) {
  const ::myb_pv pv0 = to_pv(y0);
  std::vector<::myb_pv> pv(dt.size());
  std::size_t failed = 0;
  const status code = ::myb_kepler_propagate_many(mu, &pv0, detail::data_or_null(dt), dt.size(),
                                                  pv.data(), &failed);
  return detail::give_many(code, failed, detail::to_states(pv));
}

inline result<std::vector<kepler_transition>> kepler_stm(double mu, const state& y0,
                                                         const std::vector<double>& dt) {
  const ::myb_pv pv0 = to_pv(y0);
  std::vector<::myb_pv> pv(dt.size());
  std::vector<double> stm(36 * dt.size());
  std::size_t failed = 0;
  const status code = ::myb_kepler_stm_many(mu, &pv0, detail::data_or_null(dt), dt.size(),
                                            pv.data(), stm.data(), &failed);
  std::vector<kepler_transition> out(dt.size());
  for (std::size_t k = 0; code == MYB_OK && k < dt.size(); ++k) {
    out[k].y = to_state(pv[k]);
    out[k].stm = detail::to_array<36>(stm.data() + 36 * k);
  }
  return detail::give_many(code, failed, std::move(out));
}

// ---- Numbers in text (myb/io.h) -----------------------------------------------------------------

// The correctly rounded number at the start of the text and the characters it took.
struct parsed_number {
  double value;
  std::size_t consumed;
};
inline result<parsed_number> parse_double(const std::string& text) {
  parsed_number out{};
  const status code = ::myb_parse_double(text.data(), text.size(), &out.value, &out.consumed);
  return detail::give(code, out);
}

// The shortest decimal that reads back to the value.
inline std::string format_double(double value) {
  char out[MYB_DOUBLE_CHARS] = {};
  const std::size_t length = ::myb_format_double(value, out);
  return std::string(out, length);
}

// `decimals` digits after the point ("%.Nf"), and `digits` significant digits ("%.(N-1)e").
inline result<std::string> format_fixed(double value, int decimals) {
  return detail::text([&](char* buffer, std::size_t size, std::size_t* length) {
    return ::myb_format_fixed(value, decimals, buffer, size, length);
  });
}
inline result<std::string> format_scientific(double value, int digits) {
  return detail::text([&](char* buffer, std::size_t size, std::size_t* length) {
    return ::myb_format_scientific(value, digits, buffer, size, length);
  });
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

  result<void> load_text(const std::string& text) {
    ::myb_diag diag{};
    return detail::outcome(::myb_kpool_load_text(pool_.get(), text.data(), text.size(), &diag),
                           diag);
  }

  struct variable {
    ::myb_kvar_type type;
    std::size_t count;
  };
  result<variable> describe(const std::string& name) const {
    variable out{};
    const status code = ::myb_kpool_describe(pool_.get(), name.c_str(), &out.type, &out.count);
    return detail::give(code, out);
  }
  result<std::vector<double>> numbers(const std::string& name) const {
    const double* values = nullptr;
    std::size_t count = 0;
    const status code = ::myb_kpool_numbers(pool_.get(), name.c_str(), &values, &count);
    if (code != MYB_OK) {
      return code;
    }
    return std::vector<double>(values, values + count);
  }
  result<std::string> string(const std::string& name, std::size_t index = 0) const {
    const char* value = nullptr;
    std::size_t length = 0;
    const status code = ::myb_kpool_string(pool_.get(), name.c_str(), index, &value, &length);
    if (code != MYB_OK) {
      return code;
    }
    return std::string(value, length);
  }

  ::myb_kpool* get() const noexcept { return pool_.get(); }

 private:
  explicit kernel_pool(::myb_kpool* pool) noexcept : pool_(pool) {}
  std::unique_ptr<::myb_kpool, detail::kpool_deleter> pool_;
};

// The array forms of the instants (myb/time.h).
inline std::vector<::myb_time> time_from_seconds(const std::vector<double>& seconds) {
  std::vector<::myb_time> t(seconds.size());
  ::myb_time_from_seconds_many(detail::data_or_null(seconds), seconds.size(), t.data());
  return t;
}

inline std::vector<double> time_seconds(const std::vector<::myb_time>& t) {
  std::vector<double> seconds(t.size());
  ::myb_time_seconds_many(detail::data_or_null(t), t.size(), seconds.data());
  return seconds;
}

inline result<::myb_time> calendar_to_time(const ::myb_calendar& calendar) {
  ::myb_time t{};
  const status code = ::myb_calendar_to_time(&calendar, &t);
  return detail::give(code, t);
}

inline result<std::vector<::myb_time>> calendar_to_time(
    const std::vector<::myb_calendar>& calendar) {
  std::vector<::myb_time> t(calendar.size());
  std::size_t failed = 0;
  const status code = ::myb_calendar_to_time_many(detail::data_or_null(calendar), calendar.size(),
                                                  t.data(), &failed);
  return detail::give_many(code, failed, std::move(t));
}

inline std::vector<::myb_calendar> time_to_calendar(const std::vector<::myb_time>& t) {
  std::vector<::myb_calendar> calendar(t.size());
  ::myb_time_to_calendar_many(detail::data_or_null(t), t.size(), calendar.data());
  return calendar;
}

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

  result<std::vector<::myb_time>> convert(const std::vector<::myb_time>& t, ::myb_scale from,
                                          ::myb_scale to) const {
    std::vector<::myb_time> out(t.size());
    std::size_t failed = 0;
    const status code = ::myb_timesys_convert_many(system_.get(), detail::data_or_null(t),
                                                   t.size(), from, to, out.data(), &failed);
    return detail::give_many(code, failed, std::move(out));
  }

  result<::myb_time> from_utc(const ::myb_calendar& utc, ::myb_scale to) const {
    ::myb_time out{};
    const status code = ::myb_timesys_from_utc(system_.get(), &utc, to, &out);
    if (code != MYB_OK) {
      return code;
    }
    return out;
  }

  result<std::vector<::myb_time>> from_utc(const std::vector<::myb_calendar>& utc,
                                           ::myb_scale to) const {
    std::vector<::myb_time> out(utc.size());
    std::size_t failed = 0;
    const status code = ::myb_timesys_from_utc_many(system_.get(), detail::data_or_null(utc),
                                                    utc.size(), to, out.data(), &failed);
    return detail::give_many(code, failed, std::move(out));
  }

  result<::myb_calendar> to_utc(::myb_time t, ::myb_scale from) const {
    ::myb_calendar utc{};
    const status code = ::myb_timesys_to_utc(system_.get(), t, from, &utc);
    if (code != MYB_OK) {
      return code;
    }
    return utc;
  }

  result<std::vector<::myb_calendar>> to_utc(const std::vector<::myb_time>& t,
                                             ::myb_scale from) const {
    std::vector<::myb_calendar> utc(t.size());
    std::size_t failed = 0;
    const status code = ::myb_timesys_to_utc_many(system_.get(), detail::data_or_null(t),
                                                  t.size(), from, utc.data(), &failed);
    return detail::give_many(code, failed, std::move(utc));
  }

  ::myb_timesys* get() const noexcept { return system_.get(); }

 private:
  explicit time_system(::myb_timesys* system) noexcept : system_(system) {}
  std::unique_ptr<::myb_timesys, detail::timesys_deleter> system_;
};

// ---- Frames -------------------------------------------------------------------------------------

struct angle_and_rate {
  double angle; // rad
  double rate;  // rad/s
};

struct angles_and_rates {
  vec3 angles; // rad
  vec3 rates;  // rad/s
};

inline vec3 rotation_apply_position(const ::myb_rotation& rotation, const vec3& position) noexcept {
  vec3 out{};
  ::myb_rotation_apply_position(&rotation, position.data(), out.data());
  return out;
}

inline state rotation_apply_state(const ::myb_rotation& rotation, const state& y) noexcept {
  state out{};
  ::myb_rotation_apply_state(&rotation, y.data(), out.data());
  return out;
}

// One rotation per position or state (array forms of myb/frames.h).
inline result<std::vector<vec3>> rotation_apply_position(
    const std::vector<::myb_rotation>& rotations, const std::vector<vec3>& positions) {
  if (rotations.size() != positions.size()) {
    return MYB_E_DOMAIN;
  }
  std::vector<vec3> out(positions.size());
  ::myb_rotation_apply_position_many(detail::data_or_null(rotations), detail::flat(positions),
                                     positions.size(), detail::flat(out));
  return out;
}

inline result<std::vector<state>> rotation_apply_state(
    const std::vector<::myb_rotation>& rotations, const std::vector<state>& states) {
  if (rotations.size() != states.size()) {
    return MYB_E_DOMAIN;
  }
  std::vector<state> out(states.size());
  ::myb_rotation_apply_state_many(detail::data_or_null(rotations), detail::flat(states),
                                  states.size(), detail::flat(out));
  return out;
}

// Any orientation provider at a series of instants.
inline result<std::vector<::myb_orientation>> orientation_many(
    const orientation_provider& orientation, const std::vector<::myb_time>& t) {
  std::vector<::myb_orientation> out(t.size());
  std::size_t failed = 0;
  const status code = ::myb_orientation_many(orientation.function(), orientation.context(),
                                             detail::data_or_null(t), t.size(), out.data(),
                                             &failed);
  return detail::give_many(code, failed, std::move(out));
}

// A constant orientation (myb_orientation_fixed).
inline orientation_provider fixed_orientation(const ::myb_rotation& rotation) {
  auto kept = std::make_shared<::myb_rotation>(rotation);
  ::myb_rotation* raw = kept.get();
  return orientation_provider(::myb_orientation_fixed, raw, std::move(kept));
}

// The CIO-based transformation of a body (myb/frames.h): X, Y, s, the rotation angle, the
// celestial -> terrestrial rotation.
class cip_model {
 public:
  static result<cip_model> create(const kernel_pool& pool, const std::string& prefix = "IAU2006A") {
    ::myb_cip* cip = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_cip_create(pool.get(), prefix.c_str(), nullptr, &cip, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return cip_model(cip);
  }

  struct pole {
    vec3 xys;   // X, Y, s (rad)
    vec3 rates; // rad per TT second
  };
  pole xys(::myb_time tt) const noexcept {
    pole out{};
    ::myb_cip_xys(cip_.get(), tt, out.xys.data(), out.rates.data());
    return out;
  }
  struct poles {
    std::vector<vec3> xys;
    std::vector<vec3> rates;
  };
  poles xys(const std::vector<::myb_time>& tt) const {
    poles out{std::vector<vec3>(tt.size()), std::vector<vec3>(tt.size())};
    ::myb_cip_xys_many(cip_.get(), detail::data_or_null(tt), tt.size(), detail::flat(out.xys),
                       detail::flat(out.rates));
    return out;
  }

  result<void> tabulate(::myb_time tt0, ::myb_time tt1, double step) {
    return detail::outcome(::myb_cip_tabulate(cip_.get(), tt0, tt1, step));
  }

  angle_and_rate era(::myb_time ut1) const noexcept {
    angle_and_rate out{};
    out.angle = ::myb_cip_era(cip_.get(), ut1, &out.rate);
    return out;
  }
  double sprime(::myb_time tt) const noexcept { return ::myb_cip_sprime(cip_.get(), tt); }

  // Celestial -> terrestrial at tt and ut1, with the EOP (none: all zero).
  ::myb_rotation terrestrial(::myb_time tt, ::myb_time ut1,
                             const std::optional<::myb_eop_values>& eop = std::nullopt) const
      noexcept {
    ::myb_rotation out{};
    ::myb_cip_terrestrial(cip_.get(), tt, ut1, eop ? &*eop : nullptr, &out);
    return out;
  }
  // The same over series; eop empty: all zero at every instant.
  result<std::vector<::myb_rotation>> terrestrial(
      const std::vector<::myb_time>& tt, const std::vector<::myb_time>& ut1,
      const std::vector<::myb_eop_values>& eop = {}) const {
    if (ut1.size() != tt.size() || (!eop.empty() && eop.size() != tt.size())) {
      return MYB_E_DOMAIN;
    }
    std::vector<::myb_rotation> out(tt.size());
    ::myb_cip_terrestrial_many(cip_.get(), detail::data_or_null(tt), detail::data_or_null(ut1),
                               detail::data_or_null(eop), tt.size(), out.data());
    return out;
  }

  ::myb_cip* get() const noexcept { return cip_.get(); }

 private:
  explicit cip_model(::myb_cip* cip) noexcept : cip_(cip) {}
  std::unique_ptr<::myb_cip, detail::cip_deleter> cip_;
};

// C (celestial -> intermediate) from X, Y, s and their rates.
inline ::myb_rotation cip_c2i(const vec3& xys, const vec3& rates) noexcept {
  ::myb_rotation out{};
  ::myb_cip_c2i(xys.data(), rates.data(), &out);
  return out;
}

// Earth orientation parameters of the IERS (myb/frames.h).
class eop_series {
 public:
  static result<eop_series> create() {
    ::myb_eop* eop = nullptr;
    const status code = ::myb_eop_create(nullptr, &eop);
    if (code != MYB_OK) {
      return code;
    }
    return eop_series(eop);
  }

  result<void> load_c04(const std::string& text) {
    ::myb_diag diag{};
    return detail::outcome(::myb_eop_load_c04(eop_.get(), text.data(), text.size(), &diag), diag);
  }
  result<void> load_c04_file(const std::string& path) {
    ::myb_diag diag{};
    return detail::outcome(::myb_eop_load_c04_file(eop_.get(), path.c_str(), &diag), diag);
  }
  result<void> load_finals(const std::string& text) {
    ::myb_diag diag{};
    return detail::outcome(::myb_eop_load_finals(eop_.get(), text.data(), text.size(), &diag),
                           diag);
  }
  result<void> load_finals_file(const std::string& path) {
    ::myb_diag diag{};
    return detail::outcome(::myb_eop_load_finals_file(eop_.get(), path.c_str(), &diag), diag);
  }
  result<void> set_tides(const kernel_pool& pool, const std::string& prefix = "IERS10_EOPT") {
    ::myb_diag diag{};
    return detail::outcome(::myb_eop_set_tides(eop_.get(), pool.get(), prefix.c_str(), &diag),
                           diag);
  }

  struct coverage {
    std::size_t rows;
    double first; // MJD (UTC)
    double last;
  };
  coverage rows() const noexcept {
    coverage out{};
    out.rows = ::myb_eop_rows(eop_.get(), &out.first, &out.last);
    return out;
  }

  struct values_at {
    ::myb_eop_values values;
    ::myb_time ut1;
  };
  result<values_at> at(const time_system& system, ::myb_time t, ::myb_scale scale) const {
    values_at out{};
    const status code = ::myb_eop_at(eop_.get(), system.get(), t, scale, &out.values, &out.ut1);
    return detail::give(code, out);
  }
  struct series {
    std::vector<::myb_eop_values> values;
    std::vector<::myb_time> ut1;
  };
  result<series> at(const time_system& system, const std::vector<::myb_time>& t,
                    ::myb_scale scale) const {
    series out{std::vector<::myb_eop_values>(t.size()), std::vector<::myb_time>(t.size())};
    std::size_t failed = 0;
    const status code =
        ::myb_eop_at_many(eop_.get(), system.get(), detail::data_or_null(t), t.size(), scale,
                          out.values.data(), out.ut1.data(), &failed);
    return detail::give_many(code, failed, std::move(out));
  }

  ::myb_eop* get() const noexcept { return eop_.get(); }

 private:
  explicit eop_series(::myb_eop* eop) noexcept : eop_(eop) {}
  std::unique_ptr<::myb_eop, detail::eop_deleter> eop_;
};

// The equinox-based systems of the Earth (myb/frames.h): MOD, TOD, TEME, PEF.
class equinox_model {
 public:
  static result<equinox_model> create(const kernel_pool& pool,
                                      const std::string& prefix = "IAU1980") {
    ::myb_equinox* equinox = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_equinox_create(pool.get(), prefix.c_str(), nullptr, &equinox, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return equinox_model(equinox);
  }

  angles_and_rates precession(::myb_time tt) const noexcept {
    angles_and_rates out{};
    ::myb_equinox_precession(equinox_.get(), tt, out.angles.data(), out.rates.data());
    return out;
  }
  angle_and_rate obliquity(::myb_time tt) const noexcept {
    angle_and_rate out{};
    out.angle = ::myb_equinox_obliquity(equinox_.get(), tt, &out.rate);
    return out;
  }
  struct nutation_angles {
    std::array<double, 2> nutation; // Delta psi, Delta eps (rad)
    std::array<double, 2> rates;
  };
  nutation_angles nutation(::myb_time tt) const noexcept {
    nutation_angles out{};
    ::myb_equinox_nutation(equinox_.get(), tt, out.nutation.data(), out.rates.data());
    return out;
  }

  ::myb_rotation rotation(::myb_equinox_frame frame, ::myb_time tt) const noexcept {
    ::myb_rotation out{};
    ::myb_equinox_rotation(equinox_.get(), frame, tt, &out);
    return out;
  }
  std::vector<::myb_rotation> rotation(::myb_equinox_frame frame,
                                       const std::vector<::myb_time>& tt) const {
    std::vector<::myb_rotation> out(tt.size());
    ::myb_equinox_rotation_many(equinox_.get(), frame, detail::data_or_null(tt), tt.size(),
                                out.data());
    return out;
  }

  angle_and_rate gmst(::myb_time ut1) const noexcept {
    angle_and_rate out{};
    out.angle = ::myb_equinox_gmst(equinox_.get(), ut1, &out.rate);
    return out;
  }
  ::myb_rotation teme_pef(::myb_time ut1) const noexcept {
    ::myb_rotation out{};
    ::myb_equinox_teme_pef(equinox_.get(), ut1, &out);
    return out;
  }
  std::vector<::myb_rotation> teme_pef(const std::vector<::myb_time>& ut1) const {
    std::vector<::myb_rotation> out(ut1.size());
    ::myb_equinox_teme_pef_many(equinox_.get(), detail::data_or_null(ut1), ut1.size(), out.data());
    return out;
  }

  ::myb_equinox* get() const noexcept { return equinox_.get(); }

 private:
  explicit equinox_model(::myb_equinox* equinox) noexcept : equinox_(equinox) {}
  std::unique_ptr<::myb_equinox, detail::equinox_deleter> equinox_;
};

// A rotation model of the WGCCRE (myb/frames.h).
class wgccre_model {
 public:
  static result<wgccre_model> create(const kernel_pool& pool, std::int64_t body) {
    ::myb_wgccre* model = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_wgccre_create(pool.get(), body, nullptr, &model, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return wgccre_model(model);
  }

  std::int32_t frame() const noexcept { return ::myb_wgccre_frame(model_.get()); }
  angles_and_rates euler(::myb_time t) const noexcept {
    angles_and_rates out{};
    ::myb_wgccre_euler(model_.get(), t, out.angles.data(), out.rates.data());
    return out;
  }
  ::myb_rotation rotation(::myb_time t) const noexcept {
    ::myb_rotation out{};
    ::myb_wgccre_rotation(model_.get(), t, &out);
    return out;
  }
  std::vector<::myb_rotation> rotation(const std::vector<::myb_time>& t) const {
    std::vector<::myb_rotation> out(t.size());
    ::myb_wgccre_rotation_many(model_.get(), detail::data_or_null(t), t.size(), out.data());
    return out;
  }

  ::myb_wgccre* get() const noexcept { return model_.get(); }

 private:
  explicit wgccre_model(::myb_wgccre* model) noexcept : model_(model) {}
  std::unique_ptr<::myb_wgccre, detail::wgccre_deleter> model_;
};

// The model as a provider (myb_orientation_wgccre); base: from the inertial axes in use to the
// reference frame of the model (none: they coincide).
inline orientation_provider wgccre_orientation(std::shared_ptr<const wgccre_model> model,
                                               const std::optional<::myb_rotation>& base =
                                                   std::nullopt) {
  struct holder {
    std::shared_ptr<const wgccre_model> model;
    ::myb_rotation base;
    ::myb_wgccre_orientation context;
  };
  auto kept = std::make_shared<holder>();
  kept->model = std::move(model);
  kept->base = base.value_or(::myb_rotation{});
  kept->context.model = kept->model->get();
  kept->context.base = base ? &kept->base : nullptr;
  ::myb_wgccre_orientation* raw = &kept->context;
  return orientation_provider(::myb_orientation_wgccre, raw, std::move(kept));
}

// A set of binary PCK files (myb/frames.h).
class binary_pck {
 public:
  static result<binary_pck> create() {
    ::myb_pck* pck = nullptr;
    const status code = ::myb_pck_create(nullptr, &pck);
    if (code != MYB_OK) {
      return code;
    }
    return binary_pck(pck);
  }

  result<void> load_file(const std::string& path) {
    ::myb_diag diag{};
    return detail::outcome(::myb_pck_load_file(pck_.get(), path.c_str(), &diag), diag);
  }
  // An image in memory, kept by the set.
  result<void> load_memory(std::vector<unsigned char> image) {
    auto kept = std::make_shared<std::vector<unsigned char>>(std::move(image));
    ::myb_diag diag{};
    const status code =
        ::myb_pck_load_memory(pck_.get(), kept->data(), kept->size(), &diag);
    if (code == MYB_OK) {
      images_.push_back(std::move(kept));
    }
    return detail::outcome(code, diag);
  }

  struct euler_angles {
    angles_and_rates angles;
    std::int32_t base;
  };
  result<euler_angles> euler(std::int32_t id, ::myb_time t) const {
    euler_angles out{};
    const status code = ::myb_pck_euler(pck_.get(), id, t, out.angles.angles.data(),
                                        out.angles.rates.data(), &out.base);
    return detail::give(code, out);
  }
  struct based_rotation {
    ::myb_rotation rotation; // base frame -> frame of class id
    std::int32_t base;
  };
  result<based_rotation> rotation(std::int32_t id, ::myb_time t) const {
    based_rotation out{};
    const status code = ::myb_pck_rotation(pck_.get(), id, t, &out.rotation, &out.base);
    return detail::give(code, out);
  }
  struct based_rotations {
    std::vector<::myb_rotation> rotations;
    std::vector<std::int32_t> bases;
  };
  result<based_rotations> rotation(std::int32_t id, const std::vector<::myb_time>& t) const {
    based_rotations out{std::vector<::myb_rotation>(t.size()), std::vector<std::int32_t>(t.size())};
    std::size_t failed = 0;
    const status code = ::myb_pck_rotation_many(pck_.get(), id, detail::data_or_null(t), t.size(),
                                                out.rotations.data(), out.bases.data(), &failed);
    return detail::give_many(code, failed, std::move(out));
  }

  ::myb_pck* get() const noexcept { return pck_.get(); }

 private:
  explicit binary_pck(::myb_pck* pck) noexcept : pck_(pck) {}
  // images_ first: destroyed after the set that reads them
  std::vector<std::shared_ptr<std::vector<unsigned char>>> images_;
  std::unique_ptr<::myb_pck, detail::pck_deleter> pck_;
};

// A frame of a PCK set as a provider (myb_orientation_pck); base as for wgccre_orientation.
inline orientation_provider pck_orientation(std::shared_ptr<const binary_pck> pck, std::int32_t id,
                                            const std::optional<::myb_rotation>& base =
                                                std::nullopt) {
  struct holder {
    std::shared_ptr<const binary_pck> pck;
    ::myb_rotation base;
    ::myb_pck_orientation context;
  };
  auto kept = std::make_shared<holder>();
  kept->pck = std::move(pck);
  kept->base = base.value_or(::myb_rotation{});
  kept->context.pck = kept->pck->get();
  kept->context.id = id;
  kept->context.base = base ? &kept->base : nullptr;
  ::myb_pck_orientation* raw = &kept->context;
  return orientation_provider(::myb_orientation_pck, raw, std::move(kept));
}

// The graph of frames (myb/frames.h). A handle: copies share one graph, which keeps alive the
// providers it was given and its time system.
class frame_graph {
 public:
  // One time scale: every provider gets the instant asked.
  static result<frame_graph> create() { return make(nullptr, nullptr); }
  // Instants on the scale named, providers converted to theirs by the time system.
  static result<frame_graph> create(std::shared_ptr<const time_system> timesys,
                                    const std::string& scale) {
    return make(std::move(timesys), &scale);
  }

  // Frame id given in frame base by the provider; scale: that of its argument (empty: the
  // graph's).
  result<void> add(std::int32_t id, std::int32_t base, const orientation_provider& orientation,
                   const std::string& scale = {}) {
    const status code =
        ::myb_frame_graph_add(core_->object, id, base, orientation.function(),
                              orientation.context(), scale.empty() ? nullptr : scale.c_str());
    if (code == MYB_OK) {
      core_->kept.push_back(orientation.owner());
    }
    return detail::outcome(code);
  }
  result<void> add_fixed(std::int32_t id, std::int32_t base, const ::myb_rotation& rotation) {
    return detail::outcome(::myb_frame_graph_add_fixed(core_->object, id, base, &rotation));
  }
  result<void> add_tkframe(const kernel_pool& pool, std::int32_t id) {
    ::myb_diag diag{};
    return detail::outcome(::myb_frame_graph_add_tkframe(core_->object, pool.get(), id, &diag),
                           diag);
  }

  result<::myb_rotation> rotation(std::int32_t from, std::int32_t to, ::myb_time t) const {
    ::myb_rotation out{};
    const status code = ::myb_frame_graph_rotation(core_->object, from, to, t, &out);
    return detail::give(code, out);
  }
  result<std::vector<::myb_rotation>> rotation(std::int32_t from, std::int32_t to,
                                               const std::vector<::myb_time>& t) const {
    std::vector<::myb_rotation> out(t.size());
    std::size_t failed = 0;
    const status code = ::myb_frame_graph_rotation_many(core_->object, from, to,
                                                        detail::data_or_null(t), t.size(),
                                                        out.data(), &failed);
    return detail::give_many(code, failed, std::move(out));
  }

  // The rotation from -> to as a provider (myb_orientation_graph); it keeps the graph alive.
  orientation_provider orientation(std::int32_t from, std::int32_t to) const {
    struct holder {
      std::shared_ptr<graph_core> graph;
      ::myb_graph_orientation context;
    };
    auto kept = std::make_shared<holder>();
    kept->graph = core_;
    kept->context = ::myb_graph_orientation{core_->object, from, to};
    ::myb_graph_orientation* raw = &kept->context;
    return orientation_provider(::myb_orientation_graph, raw, std::move(kept));
  }

  ::myb_frame_graph* get() const noexcept { return core_->object; }

 private:
  using graph_core = detail::core<::myb_frame_graph, detail::graph_deleter>;

  static result<frame_graph> make(std::shared_ptr<const time_system> timesys,
                                  const std::string* scale) {
    auto made = std::make_shared<graph_core>();
    ::myb_frame_graph_options options{};
    options.timesys = timesys ? timesys->get() : nullptr;
    options.scale = scale != nullptr ? scale->c_str() : nullptr;
    const status code = ::myb_frame_graph_create(&options, nullptr, &made->object);
    if (code != MYB_OK) {
      return code;
    }
    made->kept.push_back(std::move(timesys));
    return frame_graph(std::move(made));
  }
  explicit frame_graph(std::shared_ptr<graph_core> made) noexcept : core_(std::move(made)) {}
  std::shared_ptr<graph_core> core_;
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

  // An image in memory (what an SPK writer gives), kept by the set.
  result<void> load_memory(std::vector<unsigned char> image) {
    auto kept = std::make_shared<std::vector<unsigned char>>(std::move(image));
    ::myb_diag diag{};
    const status code = ::myb_spk_load_memory(spk_.get(), kept->data(), kept->size(), &diag);
    if (code == MYB_OK) {
      images_.push_back(std::move(kept));
    }
    return detail::outcome(code, diag);
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

  // The same at a series of instants.
  result<std::vector<state>> at(std::int32_t target, std::int32_t center,
                                const std::vector<::myb_time>& t) const {
    std::vector<state> states(t.size());
    std::size_t failed = 0;
    const status code =
        ::myb_spk_state_many(spk_.get(), target, center, detail::data_or_null(t), t.size(),
                             detail::flat(states), nullptr, &failed);
    return detail::give_many(code, failed, std::move(states));
  }

  // The state in the frame given, the segments turned into it by the graph.
  result<state> at_in(const frame_graph& graph, std::int32_t target, std::int32_t center,
                      ::myb_time t, std::int32_t frame) const {
    state y{};
    const status code =
        ::myb_spk_state_in(spk_.get(), graph.get(), target, center, t, frame, y.data());
    return detail::give(code, y);
  }
  result<std::vector<state>> at_in(const frame_graph& graph, std::int32_t target,
                                   std::int32_t center, const std::vector<::myb_time>& t,
                                   std::int32_t frame) const {
    std::vector<state> states(t.size());
    std::size_t failed = 0;
    const status code = ::myb_spk_state_in_many(spk_.get(), graph.get(), target, center,
                                                detail::data_or_null(t), t.size(), frame,
                                                detail::flat(states), &failed);
    return detail::give_many(code, failed, std::move(states));
  }

  std::size_t segment_count() const noexcept { return ::myb_spk_segment_count(spk_.get()); }
  ::myb_spk_segment_info segment(std::size_t index) const noexcept {
    ::myb_spk_segment_info info{};
    ::myb_spk_segment(spk_.get(), index, &info);
    return info;
  }

  // A body of the set as a provider of an environment: state = myb_spk_body_state, context =
  // the pointer returned, which lives as long as this object. frame 0 accepts any frame.
  ::myb_spk_body* body(std::int32_t target, std::int32_t center, std::int32_t frame = 0) {
    bodies_.push_back(std::make_unique<::myb_spk_body>(
        ::myb_spk_body{spk_.get(), target, center, frame, nullptr}));
    return bodies_.back().get();
  }

  ::myb_spk* get() const noexcept { return spk_.get(); }

 private:
  explicit ephemeris(::myb_spk* spk) noexcept : spk_(spk) {}
  std::vector<std::shared_ptr<std::vector<unsigned char>>> images_;
  std::unique_ptr<::myb_spk, detail::spk_deleter> spk_;
  std::vector<std::unique_ptr<::myb_spk_body>> bodies_;
};

// A body of an SPK set as a provider (myb_spk_body_state) that keeps the set alive; with a graph
// the state is given in `frame` whatever frames the segments are in.
inline provider spk_provider(std::shared_ptr<const ephemeris> eph, std::int32_t target,
                             std::int32_t center, std::int32_t frame = 0,
                             std::shared_ptr<const frame_graph> graph = nullptr) {
  struct holder {
    std::shared_ptr<const ephemeris> eph;
    std::shared_ptr<const frame_graph> graph;
    ::myb_spk_body body;
  };
  auto kept = std::make_shared<holder>();
  kept->body = ::myb_spk_body{eph->get(), target, center, frame,
                              graph ? graph->get() : nullptr};
  kept->eph = std::move(eph);
  kept->graph = std::move(graph);
  ::myb_spk_body* raw = &kept->body;
  return provider(::myb_spk_body_state, raw, std::move(kept));
}

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

// A function of the state whose zeros are looked for (an event of an orbit, a break of a force
// model): g(t, y). It must not throw.
using event_function = std::function<double(double t, const state& y)>;

namespace detail {
// What C calls for it; noexcept for the same reason as call_term.
inline double call_event(void* context, double t, const double* y) noexcept {
  const auto& g = *static_cast<const event_function*>(context);
  const state now{y[0], y[1], y[2], y[3], y[4], y[5]};
  return g(t, now);
}
static_assert(noexcept(call_event(nullptr, 0.0, nullptr)),
              "a function called from C must not let an exception through");
} // namespace detail

// A force on the object, as a callable: it writes into a (zeros when it is called) the
// acceleration of the term at the model time t and the state (x, v) relative to the origin of
// the model, and returns MYB_OK — or the status that stops the evaluation. It must not throw.
using term_function = std::function<status(double t, const vec3& x, const vec3& v, vec3& a)>;

namespace detail {
// What C calls for such a term. It is noexcept, and the assertion keeps it so: an exception
// leaving the callable ends the program here and never unwinds through the frames of C.
inline status call_term(void* context, double t, const double x[3], const double v[3],
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
static_assert(noexcept(call_term(nullptr, 0.0, nullptr, nullptr, nullptr)),
              "a function called from C must not let an exception through");
} // namespace detail

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

  // The rotating frame of the circular restricted three-body problem (myb/cr3bp.h).
  result<std::size_t> add_cr3bp(const std::string& name, double mu) {
    auto& term = keep(mu);
    return add_raw(name, ::myb_cr3bp_term, ::myb_cr3bp_term_jacobian, &term);
  }

  // Burn arcs of constant thrust (myb/thrust.h; the structure as myb_thrust_init made it).
  result<std::size_t> add_thrust(const std::string& name, const ::myb_thrust& thrust) {
    auto& term = keep(thrust);
    return add_raw(name, ::myb_thrust_acceleration, nullptr, &term);
  }

  // The engine of a Q-law (myb_term_qlaw); the law is copied and its copy pointed to.
  result<std::size_t> add_qlaw(const std::string& name, const ::myb_qlaw& law,
                               const ::myb_qlaw_term& engine) {
    ::myb_qlaw_term made = engine;
    made.law = &keep(law);
    auto& term = keep(made);
    return add_raw(name, ::myb_term_qlaw, nullptr, &term);
  }

  // A term of the caller's own; its derivatives are taken by differences of the model.
  result<std::size_t> add_term(const std::string& name, term_function function) {
    auto& stored = keep(std::move(function));
    return add_raw(name, &detail::call_term, nullptr, &stored);
  }

  // A parameter to estimate or differentiate: *value lives in the context of a term (see keep).
  result<std::size_t> add_parameter(const std::string& name, double* value, double step) {
    std::size_t index = 0;
    const status code =
        ::myb_force_model_add_parameter(model_.get(), name.c_str(), value, step, &index);
    return detail::give(code, index);
  }
  std::size_t parameter_count() const noexcept {
    return ::myb_force_model_parameter_count(model_.get());
  }

  // A sphere of `radius` on the surface of a body of the environment (env empty: a surface
  // fixed at the origin, in the axes of the model): myb_term_contact, the structure copied;
  // with touches_as_breaks its touches and partings are breaks of the model (myb_event_contact).
  result<std::size_t> add_contact(const std::string& name, std::shared_ptr<environment> env,
                                  std::size_t body, const ::myb_contact_surface& surface,
                                  const ::myb_contact_law& law, double radius,
                                  bool touches_as_breaks = false) {
    ::myb_contact_term made{};
    made.env = env ? env->get() : nullptr;
    made.body = body;
    made.surface = surface;
    made.law = law;
    made.radius = radius;
    auto& term = keep(made);
    kept_.push_back(std::move(env));
    const result<std::size_t> index = add_raw(name, ::myb_term_contact, nullptr, &term);
    if (!index || !touches_as_breaks) {
      return index;
    }
    const result<void> added = add_break_raw(::myb_event_contact, &term);
    if (!added) {
      return added.code();
    }
    return index;
  }

  // Where the force is continuous but not smooth: a function whose sign changes there.
  result<void> add_break(event_function g) {
    auto& stored = keep(std::move(g));
    return add_break_raw(&detail::call_event, &stored);
  }
  // Any break function of the C interface; the context must live as long as the model.
  result<void> add_break_raw(::myb_break_function g, void* context) {
    return detail::outcome(::myb_force_model_add_break(model_.get(), g, context));
  }
  std::size_t break_count() const noexcept { return ::myb_force_model_break_count(model_.get()); }

  // How an orbit by IAS15 ends its steps at the breaks (from its next start).
  void set_break_scheme(::myb_break_scheme scheme) noexcept {
    ::myb_force_model_set_break_scheme(model_.get(), scheme);
  }
  ::myb_break_scheme break_scheme() const noexcept {
    return ::myb_force_model_break_scheme(model_.get());
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

  // kept_ is declared first and so destroyed last: the model goes before what it points to
  std::vector<std::shared_ptr<void>> kept_;
  std::unique_ptr<::myb_force_model, detail::model_deleter> model_;
};

// ---- Propagation --------------------------------------------------------------------------------

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
        ::myb_orbit_event(orbit_.get(), &detail::call_event,
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

// An orbit parameter message and a mean-elements message read from text in either notation.
inline result<::myb_opm> parse_opm(const std::string& text) {
  ::myb_opm opm{};
  ::myb_diag diag{};
  const status code = ::myb_opm_parse(text.data(), text.size(), &opm, &diag);
  return detail::give(code, diag, opm);
}

inline result<::myb_omm> parse_omm(const std::string& text) {
  ::myb_omm omm{};
  ::myb_diag diag{};
  const status code = ::myb_omm_parse(text.data(), text.size(), &omm, &diag);
  return detail::give(code, diag, omm);
}

// An orbit ephemeris message read from text in either notation (myb/ccsds.h).
class ephemeris_message {
 public:
  struct point {
    ::myb_calendar epoch;
    state y;
    vec3 acceleration;     // zeros where the line has none
    bool has_acceleration;
  };
  struct covariance_record {
    ::myb_calendar epoch;
    std::optional<std::string> ref_frame; // none: that of the segment
    mat6 covariance;                      // m^2, m^2/s, m^2/s^2
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

  static result<ephemeris_message> load_file(const std::string& path) {
    ::myb_oem* oem = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_oem_load_file(path.c_str(), nullptr, &oem, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return ephemeris_message(oem);
  }

  // A keyword of the header; no value if the message does not have it.
  std::optional<std::string> header(const std::string& key) const {
    const char* value = ::myb_oem_header(oem_.get(), key.c_str());
    return value != nullptr ? std::optional<std::string>{value} : std::nullopt;
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
    const int given = ::myb_oem_state(oem_.get(), segment, index, &out.epoch, out.y.data(),
                                      out.acceleration.data());
    out.has_acceleration = given == 3;
    return out;
  }

  std::size_t covariance_count(std::size_t segment) const noexcept {
    return ::myb_oem_covariance_count(oem_.get(), segment);
  }
  covariance_record covariance(std::size_t segment, std::size_t index) const {
    covariance_record out{};
    const char* frame = nullptr;
    ::myb_oem_covariance(oem_.get(), segment, index, &out.epoch, &frame, out.covariance.data());
    if (frame != nullptr) {
      out.ref_frame = std::string(frame);
    }
    return out;
  }

  ::myb_oem* get() const noexcept { return oem_.get(); }

 private:
  explicit ephemeris_message(::myb_oem* oem) noexcept : oem_(oem) {}
  std::unique_ptr<::myb_oem, detail::oem_deleter> oem_;
};

// A comprehensive message (OCM, myb/ccsds.h): keywords by name, lines as text and numbers.
class comprehensive_message {
 public:
  static result<comprehensive_message> parse(const std::string& text) {
    ::myb_ocm* ocm = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_ocm_load(text.data(), text.size(), nullptr, &ocm, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return comprehensive_message(ocm);
  }

  static result<comprehensive_message> load_file(const std::string& path) {
    ::myb_ocm* ocm = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_ocm_load_file(path.c_str(), nullptr, &ocm, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return comprehensive_message(ocm);
  }

  std::optional<std::string> header(const std::string& key) const {
    return text_or_none(::myb_ocm_header(ocm_.get(), key.c_str()));
  }
  std::optional<std::string> meta(const std::string& key) const {
    return text_or_none(::myb_ocm_meta(ocm_.get(), key.c_str()));
  }
  std::size_t block_count(::myb_ocm_kind kind) const noexcept {
    return ::myb_ocm_block_count(ocm_.get(), kind);
  }
  std::optional<std::string> value(::myb_ocm_kind kind, std::size_t block,
                                   const std::string& key) const {
    return text_or_none(::myb_ocm_value(ocm_.get(), kind, block, key.c_str()));
  }
  std::size_t line_count(::myb_ocm_kind kind, std::size_t block) const noexcept {
    return ::myb_ocm_line_count(ocm_.get(), kind, block);
  }
  std::string line_text(::myb_ocm_kind kind, std::size_t block, std::size_t line) const {
    const char* text = ::myb_ocm_line_text(ocm_.get(), kind, block, line);
    return text != nullptr ? std::string(text) : std::string();
  }

  // A line as numbers: its time (absolute: the label; otherwise seconds from EPOCH_TZERO) and
  // the values after it in the units of the message.
  struct numbers {
    bool absolute;
    ::myb_calendar epoch;
    double offset;
    std::vector<double> values;
  };
  result<numbers> line(::myb_ocm_kind kind, std::size_t block, std::size_t index) const {
    numbers out{};
    int absolute = 0;
    std::size_t count = 0;
    status code = ::myb_ocm_line(ocm_.get(), kind, block, index, &absolute, &out.epoch,
                                 &out.offset, nullptr, 0, &count);
    if (code != MYB_OK) {
      return code;
    }
    out.values.resize(count);
    code = ::myb_ocm_line(ocm_.get(), kind, block, index, &absolute, &out.epoch, &out.offset,
                          out.values.data(), out.values.size(), &count);
    out.absolute = absolute != 0;
    return detail::give(code, std::move(out));
  }

  // A line of a Cartesian trajectory in SI (count: 3, 6 or 9 values of state).
  struct trajectory_point {
    bool absolute;
    ::myb_calendar epoch;
    double offset;
    std::array<double, 9> state;
    int count;
  };
  result<trajectory_point> trajectory(std::size_t block, std::size_t index) const {
    trajectory_point out{};
    int absolute = 0;
    const status code = ::myb_ocm_state(ocm_.get(), block, index, &absolute, &out.epoch,
                                        &out.offset, out.state.data(), &out.count);
    out.absolute = absolute != 0;
    return detail::give(code, out);
  }

  ::myb_ocm* get() const noexcept { return ocm_.get(); }

 private:
  explicit comprehensive_message(::myb_ocm* ocm) noexcept : ocm_(ocm) {}
  static std::optional<std::string> text_or_none(const char* value) {
    return value != nullptr ? std::optional<std::string>{value} : std::nullopt;
  }
  std::unique_ptr<::myb_ocm, detail::ocm_deleter> ocm_;
};

// A line of data of an OCM: the label (or, without one, seconds from EPOCH_TZERO) and values.
inline result<std::string> ocm_format_line(const std::optional<::myb_calendar>& epoch,
                                           double offset, const std::vector<double>& values) {
  return detail::text([&](char* buffer, std::size_t size, std::size_t* length) {
    return ::myb_ocm_format_line(epoch ? &*epoch : nullptr, offset, detail::data_or_null(values),
                                 values.size(), buffer, size, length);
  });
}

// ---- Contact --------------------------------------------------------------------------------------

// The point of a contact surface nearest to a position, the outward normal there and the height
// of the position over it (myb/contact.h).
struct surface_point {
  vec3 point;
  vec3 normal;
  double height; // m, negative inside
};
inline result<surface_point> contact_surface_nearest(const ::myb_contact_surface& surface,
                                                     const vec3& position) {
  surface_point out{};
  const status code = ::myb_contact_surface_nearest(&surface, position.data(), out.point.data(),
                                                    out.normal.data(), &out.height);
  return detail::give(code, out);
}

// The force and torque of the surface on a body with its contact points: x, v of its centre in
// the frame of the surface, attitude from the axes of the surface to the body (none: the same),
// omega its angular velocity relative to the surface in the body axes (none: zero).
inline result<::myb_contact_result> contact_force(const ::myb_contact_surface& surface,
                                                  const ::myb_contact_law& law,
                                                  const std::vector<::myb_contact_point>& points,
                                                  const vec3& x, const vec3& v,
                                                  const std::optional<mat3>& attitude = std::nullopt,
                                                  const std::optional<vec3>& omega = std::nullopt) {
  ::myb_contact_result out{};
  const status code = ::myb_contact_force(&surface, &law, detail::data_or_null(points),
                                          points.size(), x.data(), v.data(),
                                          attitude ? attitude->data() : nullptr,
                                          omega ? omega->data() : nullptr, &out);
  return detail::give(code, out);
}

// The same with the state of a friction that sticks: the deflections of the springs, 3 a point
// in the axes of the surface (empty: none, for the smooth friction), and their rates.
struct contact_stick_result {
  ::myb_contact_result result;
  std::vector<double> deflection_rate; // 3 a point; empty without deflections
};
inline result<contact_stick_result> contact_force_stick(
    const ::myb_contact_surface& surface, const ::myb_contact_law& law,
    const std::vector<::myb_contact_point>& points, const vec3& x, const vec3& v,
    const std::vector<double>& deflection, const std::optional<mat3>& attitude = std::nullopt,
    const std::optional<vec3>& omega = std::nullopt) {
  if (!deflection.empty() && deflection.size() != 3 * points.size()) {
    return MYB_E_DOMAIN;
  }
  contact_stick_result out{::myb_contact_result{}, std::vector<double>(deflection.size())};
  const status code = ::myb_contact_force_stick(
      &surface, &law, detail::data_or_null(points), points.size(), x.data(), v.data(),
      attitude ? attitude->data() : nullptr, omega ? omega->data() : nullptr,
      detail::data_or_null(deflection), out.deflection_rate.empty() ? nullptr : out.deflection_rate.data(),
      &out.result);
  return detail::give(code, std::move(out));
}

inline result<double> contact_restitution(double damping, double speed) {
  double restitution = 0.0;
  const status code = ::myb_contact_restitution(damping, speed, &restitution);
  return detail::give(code, restitution);
}
inline result<double> contact_damping(double restitution, double speed) {
  double damping = 0.0;
  const status code = ::myb_contact_damping(restitution, speed, &damping);
  return detail::give(code, damping);
}

// ---- N bodies -----------------------------------------------------------------------------------

// Threads for the accelerations of the point masses (myb/pool.h).
class thread_pool {
 public:
  static result<thread_pool> create(int threads) {
    ::myb_pool* made = nullptr;
    const status code = ::myb_pool_create(threads, nullptr, &made);
    if (code != MYB_OK) {
      return code;
    }
    return thread_pool(made);
  }
  int threads() const noexcept { return ::myb_pool_threads(pool_.get()); }
  ::myb_pool* get() const noexcept { return pool_.get(); }

 private:
  explicit thread_pool(::myb_pool* made) noexcept : pool_(made) {}
  std::unique_ptr<::myb_pool, detail::pool_deleter> pool_;
};

// Positions and velocities of all the bodies, 3 values each (m, m/s).
struct nbody_state {
  std::vector<double> x;
  std::vector<double> v;
};

// The point masses of an N-body problem (myb/dyn.h): massive bodies first, then test bodies. A
// handle: copies share one system, and the point masses of a myb::system are a view into it.
class nbody {
 public:
  static result<nbody> create(std::size_t massive_bodies, std::size_t test_bodies,
                              const std::vector<double>& gms) {
    if (gms.size() != massive_bodies) {
      return MYB_E_DOMAIN;
    }
    auto made = std::make_shared<nbody_core>();
    const status code = ::myb_nbody_create(massive_bodies, test_bodies, gms.data(), nullptr,
                                           &made->object);
    if (code != MYB_OK) {
      return code;
    }
    ::myb_nbody* raw = made->object;
    std::vector<std::shared_ptr<const void>>* kept = &made->kept;
    return nbody(std::move(made), raw, kept);
  }

  std::size_t count() const noexcept { return ::myb_nbody_count(nbody_); }
  std::size_t massive() const noexcept { return ::myb_nbody_massive(nbody_); }
  double gm(std::size_t body) const noexcept { return ::myb_nbody_gm(nbody_, body); }

  // Threads for the accelerations (none: the calling thread); the pool is kept alive.
  void set_pool(std::shared_ptr<thread_pool> threads) {
    ::myb_nbody_set_pool(nbody_, threads ? threads->get() : nullptr);
    kept_->push_back(std::move(threads));
  }

  result<void> set_interaction(std::size_t receiver, std::size_t source,
                               ::myb_interaction level) {
    return detail::outcome(::myb_nbody_set_interaction(nbody_, receiver, source, level));
  }
  // The 1PN terms (EIH) as the default of every pair, and back to Newton.
  result<void> set_relativity(const ::myb_ppn& ppn) {
    return detail::outcome(::myb_nbody_set_relativity(nbody_, &ppn));
  }
  result<void> clear_relativity() {
    return detail::outcome(::myb_nbody_set_relativity(nbody_, nullptr));
  }

  // The accelerations at the positions x (3 per body), and with the velocities (1PN terms).
  result<std::vector<double>> acceleration(double t, const std::vector<double>& x) const {
    if (x.size() != 3 * count()) {
      return MYB_E_DOMAIN;
    }
    std::vector<double> a(x.size());
    const status code = ::myb_nbody_acceleration(nbody_, t, x.data(), a.data());
    return detail::give(code, std::move(a));
  }
  result<std::vector<double>> acceleration_xv(double t, const std::vector<double>& x,
                                              const std::vector<double>& v) const {
    if (x.size() != 3 * count() || v.size() != x.size()) {
      return MYB_E_DOMAIN;
    }
    std::vector<double> a(x.size());
    const status code = ::myb_nbody_acceleration_xv(nbody_, t, x.data(), v.data(), a.data());
    return detail::give(code, std::move(a));
  }

  result<::myb_nbody_invariants> invariants_of(const std::vector<double>& x,
                                               const std::vector<double>& v) const {
    if (x.size() != 3 * count() || v.size() != x.size()) {
      return MYB_E_DOMAIN;
    }
    ::myb_nbody_invariants out{};
    const status code = ::myb_nbody_invariants_of(nbody_, x.data(), v.data(), &out);
    return detail::give(code, out);
  }
  result<::myb_nbody_pn_invariants> pn_invariants_of(const std::vector<double>& x,
                                                     const std::vector<double>& v) const {
    if (x.size() != 3 * count() || v.size() != x.size()) {
      return MYB_E_DOMAIN;
    }
    ::myb_nbody_pn_invariants out{};
    const status code = ::myb_nbody_pn_invariants_of(nbody_, x.data(), v.data(), &out);
    return detail::give(code, out);
  }

  // The fast multipole method (order 0: exact direct sums again).
  result<void> set_approximation(int order, double theta) {
    return detail::outcome(::myb_nbody_set_approximation(nbody_, order, theta));
  }
  struct approximation_setting {
    int order; // 0: exact
    double theta;
  };
  approximation_setting approximation() const noexcept {
    approximation_setting out{};
    out.order = ::myb_nbody_approximation(nbody_, &out.theta);
    return out;
  }

  ::myb_nbody* get() const noexcept { return nbody_; }

 private:
  friend class system;
  using nbody_core = detail::core<::myb_nbody, detail::nbody_deleter>;

  nbody(std::shared_ptr<const void> anchor, ::myb_nbody* raw,
        std::vector<std::shared_ptr<const void>>* kept) noexcept
      : anchor_(std::move(anchor)), nbody_(raw), kept_(kept) {}

  std::shared_ptr<const void> anchor_; // what owns the point masses and what they point to
  ::myb_nbody* nbody_;
  std::vector<std::shared_ptr<const void>>* kept_;
};

// The Wisdom-Holman integrators of a system of point masses (myb/dyn.h). It keeps the system.
class wh {
 public:
  static result<wh> create(const nbody& system, const ::myb_wh_options& options) {
    ::myb_wh* made = nullptr;
    const status code = ::myb_wh_create(system.get(), &options, nullptr, &made);
    if (code != MYB_OK) {
      return code;
    }
    return wh(system, made);
  }

  result<void> start(double t0, const std::vector<double>& x, const std::vector<double>& v) {
    if (x.size() != 3 * system_.count() || v.size() != x.size()) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(::myb_wh_start(wh_.get(), t0, x.data(), v.data()));
  }
  result<void> step() { return detail::outcome(::myb_wh_step(wh_.get())); }
  result<void> integrate(double t_end) {
    return detail::outcome(::myb_wh_integrate(wh_.get(), t_end));
  }
  double time() const noexcept { return ::myb_wh_time(wh_.get()); }

  // The corrected state at the current node, and at any t within a step of it.
  result<nbody_state> current() {
    nbody_state out{std::vector<double>(3 * system_.count()),
                    std::vector<double>(3 * system_.count())};
    const status code = ::myb_wh_state(wh_.get(), out.x.data(), out.v.data());
    return detail::give(code, std::move(out));
  }
  result<nbody_state> dense(double t) {
    nbody_state out{std::vector<double>(3 * system_.count()),
                    std::vector<double>(3 * system_.count())};
    const status code = ::myb_wh_dense(wh_.get(), t, out.x.data(), out.v.data());
    return detail::give(code, std::move(out));
  }
  ::myb_wh_stats statistics() const noexcept { return ::myb_wh_statistics(wh_.get()); }
  ::myb_wh* get() const noexcept { return wh_.get(); }

 private:
  wh(const nbody& system, ::myb_wh* made) : system_(system), wh_(made) {}
  nbody system_; // outlives the integrator
  std::unique_ptr<::myb_wh, detail::wh_deleter> wh_;
};

// An IAS15 integrator of a system of bodies (myb_system_integrator): n coordinates, a state of
// 2 n values. It keeps the system whose accelerations it integrates.
class ias15 {
 public:
  double time() const noexcept { return ::myb_ias15_time(ias15_.get()); }
  std::vector<double> current() const {
    const double* y = ::myb_ias15_state(ias15_.get());
    return std::vector<double>(y, y + 2 * n_);
  }
  ::myb_ias15_stats statistics() const noexcept { return ::myb_ias15_statistics(ias15_.get()); }
  result<void> step(double t_limit) {
    return detail::outcome(::myb_ias15_step(ias15_.get(), t_limit));
  }
  result<void> integrate(double t_end) {
    return detail::outcome(::myb_ias15_integrate(ias15_.get(), t_end));
  }
  // The state at t within the last step (dense output).
  result<std::vector<double>> dense(double t) {
    std::vector<double> y(2 * n_);
    const status code = ::myb_ias15_dense(ias15_.get(), t, y.data());
    return detail::give(code, std::move(y));
  }
  std::size_t coordinates() const noexcept { return n_; }
  ::myb_ias15* get() const noexcept { return ias15_.get(); }

 private:
  friend class system;
  ias15(std::shared_ptr<const void> anchor, ::myb_ias15* made, std::size_t n) noexcept
      : anchor_(std::move(anchor)), ias15_(made), n_(n) {}
  std::shared_ptr<const void> anchor_;
  std::unique_ptr<::myb_ias15, detail::ias15_deleter> ias15_;
  std::size_t n_;
};

// What acts on a body beyond the levels of a system (myb_system_action): it writes the
// acceleration of the body (inertial axes) and the torque on it per unit of its mass (its own
// frame), both zero when it is called, and returns MYB_OK or the status that stops the
// evaluation. It must not throw.
using action_function = std::function<status(double t, const ::myb_system_view& view,
                                             vec3& acceleration, vec3& torque)>;

namespace detail {
inline status call_action(void* context, double t, const ::myb_system_view* view,
                          double acceleration[3], double torque[3]) noexcept {
  const auto& function = *static_cast<const action_function*>(context);
  vec3 a{};
  vec3 n{};
  const status code = function(t, *view, a, n);
  if (code != MYB_OK) {
    return code;
  }
  copy_out(a, acceleration);
  copy_out(n, torque);
  return MYB_OK;
}
static_assert(noexcept(call_action(nullptr, 0.0, nullptr, nullptr, nullptr)),
              "a function called from C must not let an exception through");
} // namespace detail

// A system of bodies with their models (myb/system.h). A handle: copies share one system, which
// keeps alive its providers, fields, actions and contact points.
class system {
 public:
  static result<system> create(const std::vector<double>& gms) {
    auto made = std::make_shared<system_core>();
    const status code = ::myb_system_create(gms.size(), gms.data(), nullptr, &made->object);
    if (code != MYB_OK) {
      return code;
    }
    return system(std::move(made));
  }

  std::size_t bodies() const noexcept { return ::myb_system_bodies(core_->object); }

  // Ties the life of something the system points to (a mesh of a light action or of a contact
  // surface) to the system.
  void keep(std::shared_ptr<const void> owner) { core_->kept.push_back(std::move(owner)); }

  // ---- the model of a body
  result<void> set_orientation(std::size_t body, const orientation_provider& orientation) {
    const status code = ::myb_system_set_orientation(core_->object, body, orientation.function(),
                                                     orientation.context());
    return kept_if(code, orientation.owner());
  }
  result<void> set_field(std::size_t body, std::shared_ptr<const gravity_field> field,
                         int degree = 0) {
    const status code = ::myb_system_set_field(core_->object, body, field->get(), degree);
    return kept_if(code, std::move(field));
  }
  result<void> set_rotation(std::size_t body, const ::myb_system_rotation& rotation) {
    return detail::outcome(::myb_system_set_rotation(core_->object, body, &rotation));
  }
  result<void> set_tide(std::size_t body, const ::myb_system_tide& tide) {
    return detail::outcome(::myb_system_set_tide(core_->object, body, &tide));
  }
  result<void> set_mass_rate(std::size_t body, double rate, double epoch) {
    return detail::outcome(::myb_system_set_mass_rate(core_->object, body, rate, epoch));
  }

  // ---- the levels of the graph
  result<void> add_field(std::size_t body, std::size_t point) {
    return detail::outcome(::myb_system_add_field(core_->object, body, point));
  }
  result<void> add_figures(std::size_t a, std::size_t b) {
    return detail::outcome(::myb_system_add_figures(core_->object, a, b));
  }
  result<void> add_tide(std::size_t body, std::size_t raiser, std::size_t point) {
    return detail::outcome(::myb_system_add_tide(core_->object, body, raiser, point));
  }
  result<void> add_distortion(std::size_t body, std::size_t raiser) {
    return detail::outcome(::myb_system_add_distortion(core_->object, body, raiser));
  }
  result<void> add_pn_pair(std::size_t a, std::size_t b, const ::myb_system_pn_pair& pair) {
    return detail::outcome(::myb_system_add_pn_pair(core_->object, a, b, &pair));
  }
  result<void> add_action(std::size_t body, action_function action) {
    auto stored = std::make_shared<action_function>(std::move(action));
    void* context = stored.get();
    return add_action_raw(body, &detail::call_action, context, std::move(stored));
  }
  // Any action of the C interface; owner (may be empty) keeps its context alive.
  result<void> add_action_raw(std::size_t body, ::myb_system_action action, void* context,
                              std::shared_ptr<const void> owner = {}) {
    const status code = ::myb_system_add_action(core_->object, body, action, context);
    return kept_if(code, std::move(owner));
  }
  // Light on the mesh of a body (myb_system_light_action, the structure copied); the shape must
  // outlive the system: keep() it.
  result<void> add_light(std::size_t body, const ::myb_system_light& light) {
    auto stored = std::make_shared<::myb_system_light>(light);
    void* context = stored.get();
    return add_action_raw(body, ::myb_system_light_action, context, std::move(stored));
  }
  result<void> add_frame_dragging(std::size_t body, std::size_t point, const vec3& spin) {
    return detail::outcome(
        ::myb_system_add_frame_dragging(core_->object, body, point, spin.data()));
  }
  // A body moved from outside (an ephemeris): its number among the guided bodies.
  result<std::size_t> add_guided(double gm, const provider& motion) {
    std::size_t index = 0;
    const status code = ::myb_system_add_guided(core_->object, gm, motion.function(),
                                                motion.context(), &index);
    if (code == MYB_OK) {
      core_->kept.push_back(motion.owner());
    }
    return detail::give(code, index);
  }
  std::size_t guided() const noexcept { return ::myb_system_guided(core_->object); }
  // A body on its contact points touches the ground (the points are copied and kept).
  result<void> add_contact(std::size_t body, std::size_t ground,
                           const ::myb_contact_surface& surface, const ::myb_contact_law& law,
                           std::vector<::myb_contact_point> points) {
    auto stored = std::make_shared<std::vector<::myb_contact_point>>(std::move(points));
    ::myb_system_contact contact{};
    contact.surface = surface;
    contact.law = law;
    contact.points = detail::data_or_null(*stored);
    contact.point_count = stored->size();
    const status code = ::myb_system_add_contact(core_->object, body, ground, &contact);
    if (code == MYB_OK) {
      contacts_->push_back(contact_record{body, ground, stored->size()});
    }
    return kept_if(code, std::move(stored));
  }

  // ---- the state
  std::size_t coordinates() const noexcept { return ::myb_system_coordinates(core_->object); }

  struct rotation_state {
    quat q;           // the attitude
    vec3 omega;       // of the mantle, in its frame (rad/s)
    vec3 omega_core;  // of the core, in the mantle frame (rad/s)
  };
  // The rotation of a body into a state y (2 n values) and out of it.
  result<void> pack_rotation(std::size_t body, const quat& q, const vec3& omega,
                             const std::optional<vec3>& omega_core,
                             std::vector<double>& y) const {
    if (y.size() != 2 * coordinates()) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(::myb_system_pack_rotation(core_->object, body, q.data(),
                                                      omega.data(),
                                                      omega_core ? omega_core->data() : nullptr,
                                                      y.data()));
  }
  result<rotation_state> unpack_rotation(std::size_t body, const std::vector<double>& y) const {
    if (y.size() != 2 * coordinates()) {
      return MYB_E_DOMAIN;
    }
    rotation_state out{};
    const status code =
        ::myb_system_unpack_rotation(core_->object, body, y.data(), out.q.data(),
                                     out.omega.data(), out.omega_core.data());
    return detail::give(code, out);
  }
  struct spin_pair {
    vec3 chi1; // of body a of the pair
    vec3 chi2; // of body b
  };
  result<void> pack_spins(std::size_t pair, const vec3& chi1, const vec3& chi2,
                          std::vector<double>& y) const {
    if (y.size() != 2 * coordinates()) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(
        ::myb_system_pack_spins(core_->object, pair, chi1.data(), chi2.data(), y.data()));
  }
  result<spin_pair> unpack_spins(std::size_t pair, const std::vector<double>& y) const {
    if (y.size() != 2 * coordinates()) {
      return MYB_E_DOMAIN;
    }
    spin_pair out{};
    const status code =
        ::myb_system_unpack_spins(core_->object, pair, y.data(), out.chi1.data(), out.chi2.data());
    return detail::give(code, out);
  }
  // The deflections of the springs of a contact whose friction sticks (the contacts in place,
  // in the order they were added): 3 a point, in the axes of its surface (m).
  result<void> pack_contact(std::size_t contact, const std::vector<double>& deflection,
                            std::vector<double>& y) const {
    if (contact >= contacts_->size() || deflection.size() != 3 * (*contacts_)[contact].points ||
        y.size() != 2 * coordinates()) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(
        ::myb_system_pack_contact(core_->object, contact, deflection.data(), y.data()));
  }
  result<std::vector<double>> unpack_contact(std::size_t contact,
                                             const std::vector<double>& y) const {
    if (contact >= contacts_->size() || y.size() != 2 * coordinates()) {
      return MYB_E_DOMAIN;
    }
    std::vector<double> deflection(3 * (*contacts_)[contact].points);
    const status code =
        ::myb_system_unpack_contact(core_->object, contact, y.data(), deflection.data());
    return detail::give(code, std::move(deflection));
  }

  // ---- integration
  // An IAS15 integrator for the system as it stands (of the options the precision, the step
  // limits, the budget and the timescale are used); start it with start().
  result<ias15> integrator(const ::myb_ias15_options* options = nullptr) const {
    ::myb_ias15* made = nullptr;
    const status code = ::myb_system_integrator(core_->object, options, nullptr, &made);
    if (code != MYB_OK) {
      return code;
    }
    return ias15(core_, made, coordinates());
  }
  result<void> start(ias15& ias, double t0, const std::vector<double>& y0) {
    if (y0.size() != 2 * ias.coordinates()) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(::myb_ias15_start(ias.get(), ::myb_system_acceleration,
                                             core_->object, t0, y0.data()));
  }
  // The second derivatives of the coordinates (n values each).
  result<std::vector<double>> acceleration(double t, const std::vector<double>& x,
                                           const std::vector<double>& v) {
    if (x.size() != coordinates() || v.size() != x.size()) {
      return MYB_E_DOMAIN;
    }
    std::vector<double> a(x.size());
    const status code =
        ::myb_system_acceleration(core_->object, t, x.data(), v.data(), a.data());
    return detail::give(code, std::move(a));
  }
  result<mat3> inertia(std::size_t body, double t, const std::vector<double>& y) {
    if (y.size() != 2 * coordinates()) {
      return MYB_E_DOMAIN;
    }
    mat3 out{};
    const status code = ::myb_system_inertia(core_->object, body, t, y.data(), out.data());
    return detail::give(code, out);
  }
  // The derivatives of the state at t1 by the state at t0 and by the GM of the bodies given.
  struct derivatives {
    std::vector<double> y1;          // 2 n
    std::vector<double> phi;         // 2 n x 2 n, by rows
    std::vector<double> sensitivity; // 2 n x the bodies given
  };
  result<derivatives> partials(double t0, const std::vector<double>& y0, double t1,
                               const std::vector<std::size_t>& which = {},
                               const ::myb_ias15_options* options = nullptr) {
    const std::size_t size = 2 * coordinates();
    if (y0.size() != size) {
      return MYB_E_DOMAIN;
    }
    derivatives out{std::vector<double>(size), std::vector<double>(size * size),
                    std::vector<double>(size * which.size())};
    const status code = ::myb_system_partials(
        core_->object, options, t0, y0.data(), t1, detail::data_or_null(which), which.size(),
        out.y1.data(), out.phi.data(), which.empty() ? nullptr : out.sensitivity.data());
    return detail::give(code, std::move(out));
  }

  // ---- encounters
  result<void> set_radius(std::size_t body, double r) {
    return detail::outcome(::myb_system_set_radius(core_->object, body, r));
  }
  double radius(std::size_t body) const noexcept {
    return ::myb_system_radius(core_->object, body);
  }
  bool present(std::size_t body) const noexcept {
    return ::myb_system_present(core_->object, body) != 0;
  }
  // Steps the integrator towards t_end and stops at the first contact (or closest approach
  // within watch).
  result<::myb_system_event> advance(ias15& ias, double t_end, double watch = 0.0) {
    ::myb_system_event event{};
    const status code = ::myb_system_advance(core_->object, ias.get(), t_end, watch, &event);
    return detail::give(code, event);
  }
  result<::myb_system_merger> merge(ias15& ias, std::size_t a, std::size_t b) {
    ::myb_system_merger merger{};
    const status code = ::myb_system_merge(core_->object, ias.get(), a, b, &merger);
    if (code == MYB_OK) { // the system drops every contact that names b, the others keep order
      std::vector<contact_record> kept;
      for (const contact_record& record : *contacts_) {
        if (record.body != b && record.ground != b) {
          kept.push_back(record);
        }
      }
      *contacts_ = std::move(kept);
    }
    return detail::give(code, merger);
  }
  result<void> bounce(ias15& ias, std::size_t a, std::size_t b, double restitution) {
    return detail::outcome(::myb_system_bounce(core_->object, ias.get(), a, b, restitution));
  }

  // The point masses (myb_system_nbody): the graph, the relativity, the threads are set there.
  myb::nbody nbody() const {
    return myb::nbody(core_, ::myb_system_nbody(core_->object), &core_->kept);
  }

  ::myb_system* get() const noexcept { return core_->object; }

 private:
  using system_core = detail::core<::myb_system, detail::system_deleter>;
  // A contact in place: what sizes the deflections of its points.
  struct contact_record {
    std::size_t body;
    std::size_t ground;
    std::size_t points;
  };
  explicit system(std::shared_ptr<system_core> made)
      : core_(std::move(made)), contacts_(std::make_shared<std::vector<contact_record>>()) {}

  result<void> kept_if(status code, std::shared_ptr<const void> owner) {
    if (code == MYB_OK) {
      core_->kept.push_back(std::move(owner));
    }
    return detail::outcome(code);
  }

  std::shared_ptr<system_core> core_;
  std::shared_ptr<std::vector<contact_record>> contacts_; // shared by the copies, as the system
};

// ---- Estimation ---------------------------------------------------------------------------------

// Values computed from a state and their derivatives (myb/estim.h).
struct measurement {
  std::vector<double> values;
  std::vector<double> partials; // d value_i / d y_j at 6 i + j
};

// A measurement as a callable: it writes the values (as many as the measurement has) and, when
// partials is not null, their derivatives by the state (6 a value), and returns MYB_OK or the
// status that stops the estimation. It must not throw.
using measure_function =
    std::function<status(double t, const state& y, double* values, double* partials)>;

namespace detail {
inline status call_measure(void* context, double t, const double y[6], double* values,
                           double* partials) noexcept {
  const auto& function = *static_cast<const measure_function*>(context);
  const state now{y[0], y[1], y[2], y[3], y[4], y[5]};
  return function(t, now, values, partials);
}
static_assert(noexcept(call_measure(nullptr, 0.0, nullptr, nullptr, nullptr)),
              "a function called from C must not let an exception through");
} // namespace detail

// A measurement of `count` values: a callable of the caller's or a measurement of the C
// interface with its context (and what keeps it alive). Copies share the callable.
class measure {
 public:
  template <class F, class = std::enable_if_t<std::is_constructible<measure_function, F>::value>>
  measure(std::size_t count, F&& function) : count_(count), function_(&detail::call_measure) {
    auto stored = std::make_shared<measure_function>(std::forward<F>(function));
    context_ = stored.get();
    owner_ = std::move(stored);
  }
  measure(std::size_t count, ::myb_measure_function function, void* context,
          std::shared_ptr<const void> owner = {}) noexcept
      : count_(count), function_(function), context_(context), owner_(std::move(owner)) {}

  std::size_t count() const noexcept { return count_; }
  ::myb_measure_function function() const noexcept { return function_; }
  void* context() const noexcept { return context_; }
  const std::shared_ptr<const void>& owner() const noexcept { return owner_; }

  // The values at the state y at the model time t, with their derivatives.
  result<measurement> evaluate(double t, const state& y) const {
    measurement out{std::vector<double>(count_), std::vector<double>(6 * count_)};
    const status code = function_(context_, t, y.data(), out.values.data(), out.partials.data());
    return detail::give(code, std::move(out));
  }

 private:
  std::size_t count_;
  ::myb_measure_function function_;
  void* context_ = nullptr;
  std::shared_ptr<const void> owner_;
};

// The position (3 values) and the whole state (6) of the object itself.
inline measure measure_position() {
  return measure(3, ::myb_measure_position, nullptr);
}
inline measure measure_state() {
  return measure(6, ::myb_measure_state, nullptr);
}

// A measurement taken at a point of the object offset from its centre of mass and turned by its
// attitude (myb_measure_offset); attitude_owner keeps the context of the attitude alive.
inline measure measure_offset(const measure& inner, ::myb_attitude_function attitude,
                              void* attitude_context, const vec3& offset, double gm,
                              std::shared_ptr<const void> attitude_owner = {}) {
  struct holder {
    measure inner;
    std::shared_ptr<const void> attitude_owner;
    ::myb_offset_measure offset;
  };
  auto kept = std::make_shared<holder>(holder{inner, std::move(attitude_owner), {}});
  kept->offset.measure = inner.function();
  kept->offset.context = inner.context();
  kept->offset.attitude = attitude;
  kept->offset.attitude_context = attitude_context;
  detail::copy_out(offset, kept->offset.offset);
  kept->offset.gm = gm;
  ::myb_offset_measure* raw = &kept->offset;
  return measure(inner.count(), ::myb_measure_offset, raw, std::move(kept));
}

// An observation: the values observed at the model time t by a measurement, with their
// standard deviations (myb_observation); it keeps the measurement alive.
class observation {
 public:
  // periodic: bit i for a value that is an angle; bias: 0, or the first systematic shift + 1.
  static result<observation> create(const measure& by, double t, const std::vector<double>& value,
                                    const std::vector<double>& sigma, unsigned periodic = 0,
                                    int bias = 0) {
    if (by.count() == 0 || by.count() > MYB_MEASURE_MAX || value.size() != by.count() ||
        sigma.size() != by.count()) {
      return MYB_E_DOMAIN;
    }
    ::myb_observation made{};
    made.t = t;
    made.measure = by.function();
    made.context = by.context();
    made.count = by.count();
    for (std::size_t k = 0; k < by.count(); ++k) {
      made.value[k] = value[k];
      made.sigma[k] = sigma[k];
    }
    made.periodic = periodic;
    made.bias = bias;
    return observation(by, made);
  }

  const ::myb_observation& get() const noexcept { return observation_; }
  double time() const noexcept { return observation_.t; }
  std::size_t count() const noexcept { return observation_.count; }

 private:
  observation(const measure& by, const ::myb_observation& made) : measure_(by), observation_(made) {}
  measure measure_;
  ::myb_observation observation_;
};

namespace detail {
struct station_core {
  ::myb_station station{};
  ::myb_station_motion motion{};
  std::shared_ptr<environment> env;
  std::shared_ptr<const void> motion_owner;
};
} // namespace detail

// An observer on a body of an environment (myb_station), and its measurements. A handle: the
// measurements made from it keep it, and it keeps its environment.
class station {
 public:
  static result<station> create(std::shared_ptr<environment> env, std::size_t body,
                                const ::myb_ellipsoid& ellipsoid, double longitude,
                                double latitude, double height) {
    auto made = std::make_shared<detail::station_core>();
    const status code = ::myb_station_init(&made->station, env->get(), body, &ellipsoid,
                                           longitude, latitude, height);
    if (code != MYB_OK) {
      return code;
    }
    made->env = std::move(env);
    return station(std::move(made));
  }

  // What moves the station on its body (copied; owner keeps the tides and loads it points to).
  void set_motion(const ::myb_station_motion& motion, std::shared_ptr<const void> owner = {}) {
    core_->motion = motion;
    core_->motion_owner = std::move(owner);
    core_->station.motion = &core_->motion;
  }

  // The state of the station at the model time t (origin and axes of the environment).
  result<state> at(double t) const {
    state y{};
    const status code = ::myb_station_state(&core_->station, t, y.data());
    return detail::give(code, y);
  }

  measure measure_range() const { return geometric(1, ::myb_measure_range); }
  measure measure_range_rate() const { return geometric(1, ::myb_measure_range_rate); }
  measure measure_radec() const { return geometric(2, ::myb_measure_radec); }
  measure measure_azel() const { return geometric(2, ::myb_measure_azel); }

  // With the light time (myb_light_measure): c the speed of light, gm of the body at the origin.
  measure measure_light_range(double c, double gm, bool two_way) const {
    return with_light(1, ::myb_measure_light_range, c, gm, two_way, false);
  }
  measure measure_light_range_rate(double c, double gm, bool two_way) const {
    return with_light(1, ::myb_measure_light_range_rate, c, gm, two_way, false);
  }
  measure measure_light_radec(double c, double gm, bool two_way, bool aberration) const {
    return with_light(2, ::myb_measure_light_radec, c, gm, two_way, aberration);
  }
  measure measure_light_azel(double c, double gm, bool two_way, bool aberration) const {
    return with_light(2, ::myb_measure_light_azel, c, gm, two_way, aberration);
  }

  const ::myb_station* get() const noexcept { return &core_->station; }

 private:
  explicit station(std::shared_ptr<detail::station_core> made) noexcept : core_(std::move(made)) {}

  measure geometric(std::size_t count, ::myb_measure_function function) const {
    return measure(count, function, &core_->station, core_);
  }
  measure with_light(std::size_t count, ::myb_measure_function function, double c, double gm,
                     bool two_way, bool aberration) const {
    struct holder {
      std::shared_ptr<detail::station_core> core;
      ::myb_light_measure light;
    };
    auto kept = std::make_shared<holder>();
    kept->core = core_;
    kept->light.station = &core_->station;
    kept->light.c = c;
    kept->light.gm = gm;
    kept->light.two_way = two_way ? 1 : 0;
    kept->light.aberration = aberration ? 1 : 0;
    ::myb_light_measure* raw = &kept->light;
    return measure(count, function, raw, std::move(kept));
  }

  std::shared_ptr<detail::station_core> core_;
};

// Least squares over an arc (myb_fit_orbit).
struct fit_options {
  int max_iterations = 0;            // 0 = 10
  double tolerance = 0.0;            // 0 = 1e-6
  double epsilon = 0.0;              // 0 = 1e-10
  std::vector<unsigned char> fixed;  // empty: all estimated; else 6 + np + biases flags
  std::vector<double> prior_sigma;   // empty: none; else 6 + np + biases deviations
  std::vector<double> bias;          // the systematic shifts: their first guesses
};

struct fit_outcome {
  state y0;                       // the estimate at t0
  ::myb_fit_result summary;
  std::vector<double> covariance; // (6 + np + biases)^2
  std::vector<double> residuals;  // observed - computed of every value, at the estimate
  std::vector<double> bias;       // the estimated shifts
};

// Fits y0 and the parameters of the model (in place) to the observations, ordered in time away
// from t0.
inline result<fit_outcome> fit_orbit(force_model& model, double t0, const state& y0,
                                     const std::vector<observation>& observations,
                                     const fit_options& options = {}) {
  const std::size_t unknowns = 6 + model.parameter_count() + options.bias.size();
  if ((!options.fixed.empty() && options.fixed.size() != unknowns) ||
      (!options.prior_sigma.empty() && options.prior_sigma.size() != unknowns)) {
    return MYB_E_DOMAIN;
  }
  std::vector<::myb_observation> raw(observations.size());
  std::size_t values = 0;
  for (std::size_t k = 0; k < observations.size(); ++k) {
    raw[k] = observations[k].get();
    values += observations[k].count();
  }
  fit_outcome out{};
  out.y0 = y0;
  out.bias = options.bias;
  out.covariance.resize(unknowns * unknowns);
  out.residuals.resize(values);
  ::myb_fit_options settings{};
  settings.max_iterations = options.max_iterations;
  settings.tolerance = options.tolerance;
  settings.epsilon = options.epsilon;
  settings.fixed = detail::data_or_null(options.fixed);
  settings.prior_sigma = detail::data_or_null(options.prior_sigma);
  settings.bias_count = out.bias.size();
  settings.bias = out.bias.empty() ? nullptr : out.bias.data();
  const status code = ::myb_fit_orbit(model.get(), t0, out.y0.data(), detail::data_or_null(raw),
                                      raw.size(), &settings, nullptr, &out.summary,
                                      out.covariance.data(),
                                      out.residuals.empty() ? nullptr : out.residuals.data());
  return detail::give(code, std::move(out));
}

// The sequential filter with its smoother (myb_filter). It keeps its model.
class filter {
 public:
  static result<filter> create(std::shared_ptr<force_model> model,
                               const ::myb_filter_options& options = {}) {
    ::myb_filter* made = nullptr;
    const status code = ::myb_filter_create(model->get(), &options, nullptr, &made);
    if (code != MYB_OK) {
      return code;
    }
    const std::size_t parameters = model->parameter_count();
    return filter(std::move(model), made, parameters, options.bias_count);
  }

  // Starts at t0 with y0 and the covariance of the 6 + np unknowns (by rows).
  result<void> start(double t0, const state& y0, const std::vector<double>& covariance,
                     double t_end) {
    const std::size_t side = 6 + parameters_;
    if (covariance.size() != side * side) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(
        ::myb_filter_start(filter_.get(), t0, y0.data(), covariance.data(), t_end));
  }
  // Propagates to the observation and applies it: observed - predicted before the correction.
  result<std::vector<double>> update(const observation& o) {
    std::vector<double> innovation(o.count());
    const status code = ::myb_filter_update(filter_.get(), &o.get(), innovation.data());
    return detail::give(code, std::move(innovation));
  }
  result<void> predict(double t) { return detail::outcome(::myb_filter_predict(filter_.get(), t)); }

  struct estimate {
    double t;
    state y;
    std::vector<double> covariance; // (6 + np + biases)^2
  };
  estimate current() const {
    estimate out{0.0, state{}, std::vector<double>(unknown_count() * unknown_count())};
    ::myb_filter_get(filter_.get(), &out.t, out.y.data(), out.covariance.data());
    return out;
  }

  result<void> set_bias(const std::vector<double>& bias) {
    if (bias.size() != biases_) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(::myb_filter_set_bias(filter_.get(), detail::data_or_null(bias)));
  }
  std::vector<double> bias() const {
    std::vector<double> out(biases_);
    if (!out.empty()) {
      ::myb_filter_bias(filter_.get(), out.data());
    }
    return out;
  }

  std::size_t steps() const noexcept { return ::myb_filter_steps(filter_.get()); }
  result<void> smooth() { return detail::outcome(::myb_filter_smooth(filter_.get())); }
  struct kept_estimate {
    double t;
    std::vector<double> unknowns;   // state, parameters, shifts
    std::vector<double> covariance;
  };
  kept_estimate step(std::size_t index) const {
    kept_estimate out{0.0, std::vector<double>(unknown_count()),
                      std::vector<double>(unknown_count() * unknown_count())};
    ::myb_filter_step(filter_.get(), index, &out.t, out.unknowns.data(), out.covariance.data());
    return out;
  }

  std::size_t unknown_count() const noexcept { return 6 + parameters_ + biases_; }
  ::myb_filter* get() const noexcept { return filter_.get(); }

 private:
  filter(std::shared_ptr<force_model> model, ::myb_filter* made, std::size_t parameters,
         std::size_t biases) noexcept
      : model_(std::move(model)), filter_(made), parameters_(parameters), biases_(biases) {}
  std::shared_ptr<force_model> model_;
  std::unique_ptr<::myb_filter, detail::filter_deleter> filter_;
  std::size_t parameters_;
  std::size_t biases_;
};

// ---- Uncertainty --------------------------------------------------------------------------------

// The variational equations of a force model (myb/uncert.h). It keeps its model.
class variational {
 public:
  static result<variational> create(std::shared_ptr<force_model> model,
                                    const ::myb_variational_options& options = {}) {
    ::myb_variational* made = nullptr;
    const status code = ::myb_variational_create(model->get(), &options, nullptr, &made);
    if (code != MYB_OK) {
      return code;
    }
    const std::size_t parameters = model->parameter_count();
    return variational(std::move(model), made, parameters);
  }

  result<void> start(double t0, const state& y0) {
    return detail::outcome(::myb_variational_start(variational_.get(), t0, y0.data()));
  }
  result<void> step(double t_limit) {
    return detail::outcome(::myb_variational_step(variational_.get(), t_limit));
  }
  result<void> integrate(double t_end) {
    return detail::outcome(::myb_variational_integrate(variational_.get(), t_end));
  }
  double time() const noexcept { return ::myb_variational_time(variational_.get()); }

  struct transition {
    state y;
    mat6 phi;                        // d y / d y0, by rows
    std::vector<double> sensitivity; // d y / d p: 6 x np, by rows
  };
  transition current() const {
    transition out{state{}, mat6{}, std::vector<double>(6 * parameters_)};
    ::myb_variational_get(variational_.get(), out.y.data(), out.phi.data(),
                          out.sensitivity.empty() ? nullptr : out.sensitivity.data());
    return out;
  }
  result<transition> dense(double t) {
    transition out{state{}, mat6{}, std::vector<double>(6 * parameters_)};
    const status code =
        ::myb_variational_dense(variational_.get(), t, out.y.data(), out.phi.data(),
                                out.sensitivity.empty() ? nullptr : out.sensitivity.data());
    return detail::give(code, std::move(out));
  }
  std::size_t parameter_count() const noexcept { return parameters_; }
  ::myb_variational* get() const noexcept { return variational_.get(); }

 private:
  variational(std::shared_ptr<force_model> model, ::myb_variational* made,
              std::size_t parameters) noexcept
      : model_(std::move(model)), variational_(made), parameters_(parameters) {}
  std::shared_ptr<force_model> model_;
  std::unique_ptr<::myb_variational, detail::variational_deleter> variational_;
  std::size_t parameters_;
};

// M P M^T of n x n matrices.
inline result<std::vector<double>> covariance_transform(std::size_t n, const std::vector<double>& m,
                                                        const std::vector<double>& p) {
  if (m.size() != n * n || p.size() != n * n) {
    return MYB_E_DOMAIN;
  }
  std::vector<double> out(n * n);
  ::myb_covariance_transform(n, detail::data_or_null(m), detail::data_or_null(p),
                             out.empty() ? nullptr : out.data());
  return out;
}

// A state covariance in the radial, transverse, normal axes of (r, v).
inline result<mat6> covariance_rtn(const vec3& r, const vec3& v, const mat6& p) {
  mat6 out{};
  const status code = ::myb_covariance_rtn(r.data(), v.data(), p.data(), out.data());
  return detail::give(code, out);
}

// The deterministic random numbers of the library (xoshiro256**).
class rng {
 public:
  explicit rng(std::uint64_t seed) noexcept { ::myb_rng_seed(&rng_, seed); }
  std::uint64_t next() noexcept { return ::myb_rng_next(&rng_); }
  double uniform() noexcept { return ::myb_rng_uniform(&rng_); }
  double normal() noexcept { return ::myb_rng_normal(&rng_); }
  // A sample of N(mean, L L^T), L lower triangular (n x n by rows).
  result<std::vector<double>> gaussian(const std::vector<double>& mean,
                                       const std::vector<double>& l) {
    if (l.size() != mean.size() * mean.size()) {
      return MYB_E_DOMAIN;
    }
    std::vector<double> out(mean.size());
    ::myb_rng_gaussian(&rng_, mean.size(), detail::data_or_null(mean), detail::data_or_null(l),
                       out.empty() ? nullptr : out.data());
    return out;
  }
  ::myb_rng* get() noexcept { return &rng_; }

 private:
  ::myb_rng rng_{};
};

// The mean and covariance of count samples of dimension n, one after the other.
struct sample_moments {
  std::vector<double> mean;
  std::vector<double> covariance;
};
inline result<sample_moments> sample_statistics(std::size_t count, std::size_t n,
                                                const std::vector<double>& samples) {
  if (samples.size() != count * n) {
    return MYB_E_DOMAIN;
  }
  sample_moments out{std::vector<double>(n), std::vector<double>(n * n)};
  const status code = ::myb_sample_statistics(count, n, detail::data_or_null(samples),
                                              out.mean.empty() ? nullptr : out.mean.data(),
                                              out.covariance.empty() ? nullptr
                                                                     : out.covariance.data());
  return detail::give(code, std::move(out));
}

// ---- The circular restricted three-body problem -----------------------------------------------

inline vec3 cr3bp_acceleration(double mu, const vec3& x, const vec3& v) noexcept {
  vec3 a{};
  ::myb_cr3bp_acceleration(mu, x.data(), v.data(), a.data());
  return a;
}

struct jacobian3 {
  mat3 dadx; // by rows
  mat3 dadv;
};
inline jacobian3 cr3bp_jacobian(double mu, const vec3& x) noexcept {
  jacobian3 out{};
  ::myb_cr3bp_jacobian(mu, x.data(), out.dadx.data(), out.dadv.data());
  return out;
}

inline double cr3bp_jacobi(double mu, const state& y) noexcept {
  return ::myb_cr3bp_jacobi(mu, y.data());
}

// L1 .. L5.
inline result<std::array<vec3, 5>> cr3bp_lagrange(double mu) {
  double points[5][3] = {};
  const status code = ::myb_cr3bp_lagrange(mu, points);
  std::array<vec3, 5> out{};
  for (std::size_t k = 0; k < 5; ++k) {
    out[k] = detail::to_array<3>(points[k]);
  }
  return detail::give(code, out);
}

struct periodic_orbit {
  state y0;
  double period;
};

// y0 corrected to a periodic orbit of the symmetry asked for, with its period.
inline result<periodic_orbit> cr3bp_correct(double mu, const ::myb_cr3bp_correction& correction,
                                            const state& y0) {
  periodic_orbit out{y0, 0.0};
  const status code = ::myb_cr3bp_correct(mu, &correction, out.y0.data(), &out.period);
  return detail::give(code, out);
}

inline result<mat6> cr3bp_monodromy(double mu, const state& y0, double period) {
  mat6 out{};
  const status code = ::myb_cr3bp_monodromy(mu, y0.data(), period, out.data());
  return detail::give(code, out);
}

struct floquet {
  std::array<double, 6> re;
  std::array<double, 6> im;
  double index;
};
inline result<floquet> cr3bp_stability(const mat6& monodromy) {
  floquet out{};
  const status code =
      ::myb_cr3bp_stability(monodromy.data(), out.re.data(), out.im.data(), &out.index);
  return detail::give(code, out);
}

inline result<state> cr3bp_manifold(double mu, const state& y0, double period, double phase,
                                    bool unstable, int sign, double epsilon) {
  state y{};
  const status code = ::myb_cr3bp_manifold(mu, y0.data(), period, phase, unstable ? 1 : 0, sign,
                                           epsilon, y.data());
  return detail::give(code, y);
}

inline result<periodic_orbit> cr3bp_richardson(double mu, int point, double az, int branch) {
  periodic_orbit out{};
  const status code = ::myb_cr3bp_richardson(mu, point, az, branch, out.y0.data(), &out.period);
  return detail::give(code, out);
}

// Between the pulsating rotating frame of two real bodies and their common frame.
inline result<state> cr3bp_from_ephemeris(const state& body1, const state& body2, double gm1,
                                          double gm2, const state& y) {
  state out{};
  const status code =
      ::myb_cr3bp_from_ephemeris(body1.data(), body2.data(), gm1, gm2, y.data(), out.data());
  return detail::give(code, out);
}
inline result<state> cr3bp_to_ephemeris(const state& body1, const state& body2, double gm1,
                                        double gm2, const state& rotating) {
  state out{};
  const status code = ::myb_cr3bp_to_ephemeris(body1.data(), body2.data(), gm1, gm2,
                                               rotating.data(), out.data());
  return detail::give(code, out);
}

// Periodic orbits of any shape by multiple shooting, and their families. The orbit is corrected
// and continued in place (on MYB_E_NOCONV of a continuation it is left as it was).
inline result<::myb_cr3bp_orbit> cr3bp_orbit_from_state(double mu, const state& y0, double period,
                                                        std::size_t patches) {
  ::myb_cr3bp_orbit out{};
  const status code = ::myb_cr3bp_orbit_from_state(mu, y0.data(), period, patches, &out);
  return detail::give(code, out);
}
inline result<void> cr3bp_orbit_correct(double mu, ::myb_cr3bp_orbit& orbit,
                                        double tolerance = 0.0, int max_iterations = 0) {
  return detail::outcome(::myb_cr3bp_orbit_correct(mu, &orbit, tolerance, max_iterations));
}
inline result<void> cr3bp_orbit_continue(double mu, ::myb_cr3bp_orbit& orbit, double step,
                                         double tolerance = 0.0, int max_iterations = 0) {
  return detail::outcome(
      ::myb_cr3bp_orbit_continue(mu, &orbit, step, tolerance, max_iterations));
}

// Quasi-periodic orbits: invariant circles of the flow, their families, stability, manifolds.
inline result<::myb_cr3bp_torus> cr3bp_torus_from_orbit(double mu, const state& y0, double period,
                                                        std::size_t points, int pair,
                                                        double size) {
  ::myb_cr3bp_torus out{};
  const status code =
      ::myb_cr3bp_torus_from_orbit(mu, y0.data(), period, points, pair, size, &out);
  return detail::give(code, out);
}
inline result<void> cr3bp_torus_continue(double mu, ::myb_cr3bp_torus& torus,
                                         ::myb_cr3bp_torus_family family, double step,
                                         double tolerance = 0.0, int max_iterations = 0) {
  return detail::outcome(
      ::myb_cr3bp_torus_continue(mu, &torus, family, step, tolerance, max_iterations));
}
inline state cr3bp_torus_state(const ::myb_cr3bp_torus& torus, double angle) noexcept {
  state y{};
  ::myb_cr3bp_torus_state(&torus, angle, y.data());
  return y;
}
inline result<::myb_cr3bp_torus_spectrum> cr3bp_torus_stability(double mu,
                                                                const ::myb_cr3bp_torus& torus,
                                                                double tolerance = 0.0) {
  auto out = std::make_unique<::myb_cr3bp_torus_spectrum>();
  const status code = ::myb_cr3bp_torus_stability(mu, &torus, tolerance, out.get());
  if (code != MYB_OK) {
    return code;
  }
  return *out;
}
inline result<state> cr3bp_torus_manifold(const ::myb_cr3bp_torus& torus,
                                          const ::myb_cr3bp_torus_spectrum& spectrum,
                                          double angle, bool unstable, int sign, double epsilon) {
  state y{};
  const status code = ::myb_cr3bp_torus_manifold(&torus, &spectrum, angle, unstable ? 1 : 0, sign,
                                                 epsilon, y.data());
  return detail::give(code, y);
}

// The flow for a time, or to the first crossing of a section.
struct cr3bp_arrival {
  state y;
  double elapsed;
  bool crossed;
};
inline result<cr3bp_arrival> cr3bp_propagate(double mu, const state& y0, double duration,
                                             const std::optional<::myb_cr3bp_section>& section =
                                                 std::nullopt) {
  cr3bp_arrival out{};
  int crossed = 0;
  const status code = ::myb_cr3bp_propagate(mu, y0.data(), duration, section ? &*section : nullptr,
                                            out.y.data(), &out.elapsed, &crossed);
  out.crossed = crossed != 0;
  return detail::give(code, out);
}

// ---- Transfers ------------------------------------------------------------------------------------

inline result<::myb_bplane> bplane_compute(double gm, const vec3& r, const vec3& v,
                                           const vec3& pole) {
  ::myb_bplane out{};
  const status code = ::myb_bplane_compute(gm, r.data(), v.data(), pole.data(), &out);
  return detail::give(code, out);
}

inline double sphere_of_influence(double a, double mass_ratio) noexcept {
  return ::myb_sphere_of_influence(a, mass_ratio);
}
inline double hill_radius(double a, double e, double mass_ratio) noexcept {
  return ::myb_hill_radius(a, e, mass_ratio);
}
inline double flyby_turn(double gm, double v_infinity, double periapsis) noexcept {
  return ::myb_flyby_turn(gm, v_infinity, periapsis);
}
inline double departure_impulse(double gm, double radius, double v_infinity) noexcept {
  return ::myb_departure_impulse(gm, radius, v_infinity);
}

// The residuals of a shooting problem: f (m values) at u, and when jacobian is not null the
// m x n Jacobian by rows. It returns MYB_OK or the status that stops the solver; it must not
// throw.
using shoot_function = std::function<status(const double* u, double* f, double* jacobian)>;

namespace detail {
inline status call_shoot(void* context, const double* u, double* f, double* jacobian) noexcept {
  const auto& function = *static_cast<const shoot_function*>(context);
  return function(u, f, jacobian);
}
static_assert(noexcept(call_shoot(nullptr, nullptr, nullptr, nullptr)),
              "a function called from C must not let an exception through");
} // namespace detail

struct shoot_report {
  int iterations;
  double residual; // max |f_i|
};
// Newton's method for f(u) = 0 with u.size() unknowns and `equations` equations; u in place.
inline result<shoot_report> shoot(const shoot_function& function, std::size_t equations,
                                  std::vector<double>& u, const ::myb_shoot_options& options) {
  std::vector<double> work(::myb_shoot_work_size(u.size(), equations));
  shoot_report out{};
  const status code = ::myb_shoot(&detail::call_shoot, detail::context_of(&function), u.size(),
                                  equations, u.empty() ? nullptr : u.data(), &options,
                                  work.empty() ? nullptr : work.data(), &out.iterations,
                                  &out.residual);
  return detail::give(code, out);
}

// Every conic with 0 .. max_revolutions whole turns between r1 and r2 in the time tof.
inline result<std::vector<::myb_lambert_solution>> lambert(double mu, const vec3& r1,
                                                           const vec3& r2, double tof,
                                                           bool retrograde = false,
                                                           int max_revolutions = 0) {
  const std::size_t capacity =
      1 + 2 * static_cast<std::size_t>(max_revolutions > 0 ? max_revolutions : 0);
  std::vector<::myb_lambert_solution> out(capacity);
  std::size_t count = 0;
  const status code = ::myb_lambert(mu, r1.data(), r2.data(), tof, retrograde ? 1 : 0,
                                    max_revolutions, out.data(), out.size(), &count);
  out.resize(count < capacity ? count : capacity);
  return detail::give(code, std::move(out));
}

// ---- A first orbit from a few measurements --------------------------------------------------------

struct gibbs_velocity {
  vec3 v2;
  double out_of_plane; // rad
};
inline result<gibbs_velocity> iod_gibbs(double mu, const vec3& r1, const vec3& r2,
                                        const vec3& r3) {
  gibbs_velocity out{};
  const status code =
      ::myb_iod_gibbs(mu, r1.data(), r2.data(), r3.data(), out.v2.data(), &out.out_of_plane);
  return detail::give(code, out);
}

inline result<vec3> iod_herrick_gibbs(double mu, const vec3& r1, double t1, const vec3& r2,
                                      double t2, const vec3& r3, double t3) {
  vec3 v2{};
  const status code =
      ::myb_iod_herrick_gibbs(mu, r1.data(), t1, r2.data(), t2, r3.data(), t3, v2.data());
  return detail::give(code, v2);
}

namespace detail {
struct sights {
  double t[3];
  double observer[3][3];
  double sight[3][3];
};
inline sights to_sights(const std::array<double, 3>& t, const std::array<vec3, 3>& observer,
                        const std::array<vec3, 3>& sight) noexcept {
  sights out{};
  for (std::size_t k = 0; k < 3; ++k) {
    out.t[k] = t[k];
    copy_out(observer[k], out.observer[k]);
    copy_out(sight[k], out.sight[k]);
  }
  return out;
}
} // namespace detail

// Three unit sights from three observers at t[0] < t[1] < t[2] (positions relative to the
// centre): the candidate orbits, each refined to all three.
inline result<std::vector<::myb_iod_solution>> iod_gauss(double mu, const std::array<double, 3>& t,
                                                         const std::array<vec3, 3>& observer,
                                                         const std::array<vec3, 3>& sight) {
  const detail::sights in = detail::to_sights(t, observer, sight);
  ::myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS] = {};
  std::size_t count = 0;
  const status code = ::myb_iod_gauss(mu, in.t, in.observer, in.sight, solutions, &count);
  return detail::give(code, std::vector<::myb_iod_solution>(solutions, solutions + count));
}
inline result<std::vector<::myb_iod_solution>> iod_laplace(double mu,
                                                           const std::array<double, 3>& t,
                                                           const std::array<vec3, 3>& observer,
                                                           const std::array<vec3, 3>& sight) {
  const detail::sights in = detail::to_sights(t, observer, sight);
  ::myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS] = {};
  std::size_t count = 0;
  const status code = ::myb_iod_laplace(mu, in.t, in.observer, in.sight, solutions, &count);
  return detail::give(code, std::vector<::myb_iod_solution>(solutions, solutions + count));
}
inline result<::myb_iod_solution> iod_double_r(double mu, const std::array<double, 3>& t,
                                               const std::array<vec3, 3>& observer,
                                               const std::array<vec3, 3>& sight,
                                               const std::array<double, 2>& distance) {
  const detail::sights in = detail::to_sights(t, observer, sight);
  ::myb_iod_solution out{};
  const status code =
      ::myb_iod_double_r(mu, in.t, in.observer, in.sight, distance.data(), &out);
  return detail::give(code, out);
}
inline result<::myb_iod_solution> iod_gooding(double mu, const std::array<double, 3>& t,
                                              const std::array<vec3, 3>& observer,
                                              const std::array<vec3, 3>& sight,
                                              const std::array<double, 2>& range,
                                              int revolutions = 0, int branch = 0) {
  const detail::sights in = detail::to_sights(t, observer, sight);
  ::myb_iod_solution out{};
  const status code = ::myb_iod_gooding(mu, in.t, in.observer, in.sight, range.data(),
                                        revolutions, branch, &out);
  return detail::give(code, out);
}

// ---- Relative motion ------------------------------------------------------------------------------

// The relative state of a deputy about a chief in the turning frame (R, T, N), and back.
inline result<state> relative_from_inertial(const state& chief, const state& deputy) {
  const ::myb_pv c = to_pv(chief);
  const ::myb_pv d = to_pv(deputy);
  state out{};
  const status code = ::myb_relative_from_inertial(&c, &d, out.data());
  return detail::give(code, out);
}
inline result<state> relative_to_inertial(const state& chief, const state& relative) {
  const ::myb_pv c = to_pv(chief);
  ::myb_pv d{};
  const status code = ::myb_relative_to_inertial(&c, relative.data(), &d);
  return detail::give(code, to_state(d));
}
inline result<mat6> relative_stm(double mu, const state& chief, double dt) {
  const ::myb_pv c = to_pv(chief);
  mat6 out{};
  const status code = ::myb_relative_stm(mu, &c, dt, out.data());
  return detail::give(code, out);
}
inline mat6 cw_stm(double n, double dt) noexcept {
  mat6 out{};
  ::myb_cw_stm(n, dt, out.data());
  return out;
}
struct rendezvous_impulses {
  vec3 dv0;
  vec3 dv1;
};
// target_velocity: none, at rest with respect to the chief.
inline result<rendezvous_impulses> relative_rendezvous(double mu, const state& chief, double dt,
                                                       const state& relative, const vec3& target,
                                                       const std::optional<vec3>& target_velocity =
                                                           std::nullopt) {
  const ::myb_pv c = to_pv(chief);
  rendezvous_impulses out{};
  const status code = ::myb_relative_rendezvous(
      mu, &c, dt, relative.data(), target.data(),
      target_velocity ? target_velocity->data() : nullptr, out.dv0.data(), out.dv1.data());
  return detail::give(code, out);
}
inline result<std::array<double, 6>> roe_from_states(double mu, const state& chief,
                                                     const state& deputy) {
  const ::myb_pv c = to_pv(chief);
  const ::myb_pv d = to_pv(deputy);
  std::array<double, 6> out{};
  const status code = ::myb_roe_from_states(mu, &c, &d, out.data());
  return detail::give(code, out);
}
inline vec3 roe_to_position(double a, double u, const std::array<double, 6>& roe) noexcept {
  vec3 out{};
  ::myb_roe_to_position(a, u, roe.data(), out.data());
  return out;
}

// ---- Optimal transfers ----------------------------------------------------------------------------

// The indirect method (myb/optimal.h); guess: the eight costates of a solution nearby. The
// field of the problem is the caller's (its get()).
inline result<::myb_indirect_solution> indirect_solve(
    const ::myb_indirect_problem& problem,
    const std::optional<std::array<double, 8>>& guess = std::nullopt) {
  ::myb_indirect_solution out{};
  const status code =
      ::myb_indirect_solve(&problem, guess ? guess->data() : nullptr, nullptr, &out);
  return detail::give(code, out);
}
inline result<std::vector<::myb_indirect_point>> indirect_sample(
    const ::myb_indirect_problem& problem, const ::myb_indirect_solution& solution,
    const std::vector<double>& times) {
  std::vector<::myb_indirect_point> out(times.size());
  const status code = ::myb_indirect_sample(&problem, &solution, detail::data_or_null(times),
                                            times.size(), nullptr, out.data());
  return detail::give(code, std::move(out));
}

// The objective and constraints of a nonlinear program at x (myb_nlp_function): *f, c, and when
// gradient is not null also the gradient and the Jacobian of c. It must not throw.
using nlp_function = std::function<status(const double* x, double& f, double* c,
                                          double* gradient, double* jacobian)>;

namespace detail {
inline status call_nlp(void* context, const double* x, double* f, double* c, double* gradient,
                       double* jacobian) noexcept {
  const auto& function = *static_cast<const nlp_function*>(context);
  return function(x, *f, c, gradient, jacobian);
}
static_assert(noexcept(call_nlp(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr)),
              "a function called from C must not let an exception through");
} // namespace detail

struct nlp_problem {
  std::size_t n = 0;            // unknowns
  std::size_t equalities = 0;   // constraints c_i = 0, first
  std::size_t inequalities = 0; // then c_i <= 0
  std::vector<double> lower;    // empty: no bounds; else n values
  std::vector<double> upper;
  nlp_function function;
};
struct nlp_solution {
  ::myb_nlp_result report;
  std::vector<double> multipliers; // of the constraints
};
// x: the first guess, the solution on return (on MYB_E_NOCONV what was reached).
inline result<nlp_solution> nlp_solve(const nlp_problem& problem, std::vector<double>& x,
                                      const ::myb_nlp_options* options = nullptr) {
  if (x.size() != problem.n || (!problem.lower.empty() && problem.lower.size() != problem.n) ||
      (!problem.upper.empty() && problem.upper.size() != problem.n)) {
    return MYB_E_DOMAIN;
  }
  ::myb_nlp_problem made{};
  made.n = problem.n;
  made.equalities = problem.equalities;
  made.inequalities = problem.inequalities;
  made.lower = detail::data_or_null(problem.lower);
  made.upper = detail::data_or_null(problem.upper);
  made.function = &detail::call_nlp;
  made.context = detail::context_of(&problem.function);
  nlp_solution out{::myb_nlp_result{},
                   std::vector<double>(problem.equalities + problem.inequalities)};
  const status code = ::myb_nlp_solve(&made, options, x.empty() ? nullptr : x.data(),
                                      out.multipliers.empty() ? nullptr : out.multipliers.data(),
                                      nullptr, &out.report);
  return detail::give(code, std::move(out));
}

// The direct method of Sims and Flanagan: throttles 3 a segment.
struct leg_mismatch {
  std::array<double, 6> mismatch; // forward - backward at the meeting point (m, m/s)
  double mass1;                   // at the arrival (kg)
  std::vector<double> jacobian;   // 6 x 3 segments, by rows
};
inline result<leg_mismatch> sims_flanagan_mismatch(const ::myb_sims_flanagan& leg,
                                                   const std::vector<double>& throttles) {
  if (throttles.size() != 3 * leg.segments) {
    return MYB_E_DOMAIN;
  }
  leg_mismatch out{std::array<double, 6>{}, 0.0, std::vector<double>(6 * throttles.size())};
  const status code = ::myb_sims_flanagan_mismatch(&leg, throttles.data(), out.mismatch.data(),
                                                   &out.mass1, out.jacobian.data());
  return detail::give(code, std::move(out));
}
// throttles: the first guess (all zeros: the one built in), the solution on return.
inline result<::myb_sims_flanagan_result> sims_flanagan_solve(
    const ::myb_sims_flanagan& leg, std::vector<double>& throttles,
    const ::myb_nlp_options* options = nullptr) {
  if (throttles.size() != 3 * leg.segments) {
    return MYB_E_DOMAIN;
  }
  ::myb_sims_flanagan_result out{};
  const status code =
      ::myb_sims_flanagan_solve(&leg, options, throttles.data(), nullptr, &out);
  return detail::give(code, out);
}

// Chains of legs with flybys. A body of a chain from a provider (which must outlive the chain).
inline ::myb_chain_body chain_body(const provider& body, double gm = 0.0,
                                   double periapsis = 0.0) noexcept {
  ::myb_chain_body out{};
  out.state = body.function();
  out.context = body.context();
  out.gm = gm;
  out.periapsis = periapsis;
  return out;
}
// A powered flyby: the impulse at the periapsis and the turn of the excess velocity, with their
// derivatives (d dv / d(v_in, v_out, periapsis), then d turn / d(the same)).
struct flyby_impulse {
  double dv;
  double turn; // rad
  std::array<double, 6> partials;
};
inline result<flyby_impulse> chain_flyby(double gm, double v_in, double v_out, double periapsis) {
  flyby_impulse out{};
  const status code = ::myb_chain_flyby(gm, v_in, v_out, periapsis, &out.dv, &out.turn,
                                        out.partials.data());
  return detail::give(code, out);
}
struct chain_sizes {
  std::size_t unknowns;
  std::size_t equalities;
  std::size_t inequalities;
};
inline result<chain_sizes> chain_size(const ::myb_chain& chain) {
  chain_sizes out{};
  const status code =
      ::myb_chain_size(&chain, &out.unknowns, &out.equalities, &out.inequalities);
  return detail::give(code, out);
}
inline result<std::vector<double>> chain_guess(const ::myb_chain& chain) {
  const result<chain_sizes> sizes = chain_size(chain);
  if (!sizes) {
    return sizes.code();
  }
  std::vector<double> x(sizes->unknowns);
  const status code = ::myb_chain_guess(&chain, x.data());
  return detail::give(code, std::move(x));
}
struct chain_value {
  double propellant;               // kg
  std::vector<double> constraints;
  std::vector<double> gradient;    // empty without derivatives
  std::vector<double> jacobian;    // constraints x unknowns, by rows
};
inline result<chain_value> chain_evaluate(const ::myb_chain& chain, const std::vector<double>& x,
                                          bool derivatives = true) {
  const result<chain_sizes> sizes = chain_size(chain);
  if (!sizes) {
    return sizes.code();
  }
  if (x.size() != sizes->unknowns) {
    return MYB_E_DOMAIN;
  }
  const std::size_t constraints = sizes->equalities + sizes->inequalities;
  chain_value out{0.0, std::vector<double>(constraints),
                  std::vector<double>(derivatives ? x.size() : 0),
                  std::vector<double>(derivatives ? constraints * x.size() : 0)};
  const status code = ::myb_chain_evaluate(
      &chain, x.data(), &out.propellant, out.constraints.data(),
      derivatives ? out.gradient.data() : nullptr, derivatives ? out.jacobian.data() : nullptr,
      nullptr);
  return detail::give(code, std::move(out));
}
// x: the first guess (all zeros: that of chain_guess), the solution on return.
inline result<::myb_chain_result> chain_solve(const ::myb_chain& chain, std::vector<double>& x,
                                              const ::myb_nlp_options* options = nullptr) {
  const result<chain_sizes> sizes = chain_size(chain);
  if (!sizes) {
    return sizes.code();
  }
  if (x.size() != sizes->unknowns) {
    return MYB_E_DOMAIN;
  }
  auto out = std::make_unique<::myb_chain_result>();
  const status code = ::myb_chain_solve(&chain, options, x.data(), nullptr, out.get());
  if (code != MYB_OK) {
    return code;
  }
  return *out;
}

// ---- Low thrust -----------------------------------------------------------------------------------

inline result<::myb_thrust> thrust_init(double mass0, const std::vector<::myb_thrust_arc>& arcs) {
  ::myb_thrust out{};
  const status code = ::myb_thrust_init(&out, mass0, detail::data_or_null(arcs), arcs.size());
  return detail::give(code, out);
}
inline double thrust_mass(const ::myb_thrust& thrust, double t) noexcept {
  return ::myb_thrust_mass(&thrust, t);
}
inline result<vec3> thrust_acceleration(const ::myb_thrust& thrust, double t, const vec3& x,
                                        const vec3& v) {
  vec3 a{};
  const status code = ::myb_thrust_acceleration(detail::context_of(&thrust), t, x.data(),
                                                v.data(), a.data());
  return detail::give(code, a);
}
inline result<vec3> thrust_rtn(const ::myb_thrust& thrust, double t, const vec3& x,
                               const vec3& v) {
  vec3 a{};
  const status code = ::myb_thrust_rtn(&thrust, t, x.data(), v.data(), a.data());
  return detail::give(code, a);
}
inline std::array<double, 6> equinoctial_rates(double mu, const ::myb_equinoctial& elements,
                                               const vec3& a_rtn) noexcept {
  std::array<double, 6> out{};
  ::myb_equinoctial_rates(mu, &elements, a_rtn.data(), out.data());
  return out;
}
struct edelbaum_transfer {
  double dv;    // m/s
  double beta0; // rad
};
inline result<edelbaum_transfer> edelbaum(double v0, double v1, double di) {
  edelbaum_transfer out{};
  const status code = ::myb_edelbaum(v0, v1, di, &out.dv, &out.beta0);
  return detail::give(code, out);
}
struct edelbaum_point {
  double v;
  double di;
};
inline edelbaum_point edelbaum_state(double v0, double beta0, double dv) noexcept {
  edelbaum_point out{};
  ::myb_edelbaum_state(v0, beta0, dv, &out.v, &out.di);
  return out;
}

// The Q-law.
struct qlaw_steering {
  vec3 direction; // radial, transverse, normal
  double q;
  double q_rate;
  double effectivity;
};
inline result<qlaw_steering> qlaw_direction(const ::myb_qlaw& law, const ::myb_elements& now) {
  qlaw_steering out{};
  const status code = ::myb_qlaw_direction(&law, &now, out.direction.data(), &out.q,
                                           &out.q_rate, &out.effectivity);
  return detail::give(code, out);
}
inline bool qlaw_arrived(const ::myb_qlaw& law, const ::myb_elements& now,
                         const std::array<double, 5>& tolerance) noexcept {
  return ::myb_qlaw_arrived(&law, &now, tolerance.data()) != 0;
}
// The right-hand side with the mass as a seventh component (the flight is the caller's).
inline result<std::array<double, 7>> qlaw_flight_rhs(const ::myb_qlaw_flight& flight, double t,
                                                     const std::array<double, 7>& y) {
  std::array<double, 7> dydt{};
  const status code =
      ::myb_qlaw_flight_rhs(detail::context_of(&flight), t, y.data(), dydt.data());
  return detail::give(code, dydt);
}
inline result<::myb_qlaw_arrival> qlaw_fly(const ::myb_qlaw_flight& flight, double t0,
                                           const std::array<double, 7>& y0, double t_limit,
                                           double step) {
  ::myb_qlaw_arrival out{};
  const status code = ::myb_qlaw_fly(&flight, t0, y0.data(), t_limit, step, nullptr, &out);
  return detail::give(code, out);
}

// ---- Conjunctions -------------------------------------------------------------------------------

struct approaches_found {
  std::vector<::myb_approach> approaches; // at most the capacity, in the order of the time
  std::size_t found;                      // how many there were
};
inline result<approaches_found> approach_find(const provider& first, const provider& second,
                                              double t0, double t1, double step,
                                              std::size_t capacity = 64) {
  approaches_found out{std::vector<::myb_approach>(capacity), 0};
  const status code =
      ::myb_approach_find(first.function(), first.context(), second.function(), second.context(),
                          t0, t1, step, out.approaches.data(), capacity, &out.found);
  out.approaches.resize(out.found < capacity ? out.found : capacity);
  return detail::give(code, std::move(out));
}
// Two orbits started and standing at one time, taken together to t1.
inline result<approaches_found> approach_orbits(orbit& first, orbit& second, double t1,
                                                double step, std::size_t capacity = 64) {
  approaches_found out{std::vector<::myb_approach>(capacity), 0};
  const status code = ::myb_approach_orbits(first.get(), second.get(), t1, step,
                                            out.approaches.data(), capacity, &out.found);
  out.approaches.resize(out.found < capacity ? out.found : capacity);
  return detail::give(code, std::move(out));
}
inline result<::myb_encounter> encounter_plane(const vec3& position, const vec3& velocity,
                                               const mat3& covariance, double radius) {
  ::myb_encounter out{};
  const status code =
      ::myb_encounter_plane(position.data(), velocity.data(), covariance.data(), radius, &out);
  return detail::give(code, out);
}
inline result<double> collision_probability(const ::myb_encounter& encounter) {
  double probability = 0.0;
  const status code = ::myb_collision_probability(&encounter, &probability);
  return detail::give(code, probability);
}
inline double collision_probability_chan(const ::myb_encounter& encounter) noexcept {
  return ::myb_collision_probability_chan(&encounter);
}
struct probability_bound {
  double probability;
  double factor;
};
inline result<probability_bound> collision_probability_max(const ::myb_encounter& encounter) {
  probability_bound out{};
  const status code =
      ::myb_collision_probability_max(&encounter, &out.probability, &out.factor);
  return detail::give(code, out);
}
struct conjunctions_found {
  std::vector<::myb_conjunction> conjunctions;
  std::size_t found;
};
// Both orbits created with `variational`; the covariances of their states at the start.
inline result<conjunctions_found> conjunction_find(orbit& first, const mat6& first_covariance,
                                                   orbit& second, const mat6& second_covariance,
                                                   double radius, double t1, double step,
                                                   std::size_t capacity = 64) {
  conjunctions_found out{std::vector<::myb_conjunction>(capacity), 0};
  const status code = ::myb_conjunction_find(first.get(), first_covariance.data(), second.get(),
                                             second_covariance.data(), radius, t1, step,
                                             out.conjunctions.data(), capacity, &out.found);
  out.conjunctions.resize(out.found < capacity ? out.found : capacity);
  return detail::give(code, std::move(out));
}
inline result<double> collision_rate(const state& mean, const mat6& covariance, double radius) {
  double rate = 0.0;
  const status code = ::myb_collision_rate(mean.data(), covariance.data(), radius, &rate);
  return detail::give(code, rate);
}

// The relative state as a Gaussian at t: its mean and covariance. It must not throw.
using gaussian_function = std::function<status(double t, state& mean, mat6& covariance)>;

namespace detail {
inline status call_gaussian(void* context, double t, double mean[6],
                            double covariance[36]) noexcept {
  const auto& function = *static_cast<const gaussian_function*>(context);
  state m{};
  mat6 c{};
  const status code = function(t, m, c);
  if (code != MYB_OK) {
    return code;
  }
  copy_out(m, mean);
  copy_out(c, covariance);
  return MYB_OK;
}
static_assert(noexcept(call_gaussian(nullptr, 0.0, nullptr, nullptr)),
              "a function called from C must not let an exception through");
} // namespace detail

struct probability_estimate {
  double probability;
  double error;
};
inline result<probability_estimate> collision_probability_long(const gaussian_function& relative,
                                                               double t0, double t1, double step,
                                                               double radius) {
  probability_estimate out{};
  const status code = ::myb_collision_probability_long(&detail::call_gaussian,
                                                       detail::context_of(&relative), t0, t1,
                                                       step, radius, &out.probability, &out.error);
  return detail::give(code, out);
}
inline result<probability_estimate> conjunction_long(orbit& first, const mat6& first_covariance,
                                                     orbit& second, const mat6& second_covariance,
                                                     double radius, double t1, double step) {
  probability_estimate out{};
  const status code = ::myb_conjunction_long(first.get(), first_covariance.data(), second.get(),
                                             second_covariance.data(), radius, t1, step,
                                             &out.probability, &out.error);
  return detail::give(code, out);
}

// ---- Attitude -------------------------------------------------------------------------------------

inline quat quat_multiply(const quat& a, const quat& b) noexcept {
  quat out{};
  ::myb_quat_multiply(a.data(), b.data(), out.data());
  return out;
}
inline quat quat_conjugate(const quat& q) noexcept {
  quat out{};
  ::myb_quat_conjugate(q.data(), out.data());
  return out;
}
// q made unit in place; the former norm.
inline double quat_normalize(quat& q) noexcept {
  return ::myb_quat_normalize(q.data());
}
inline mat3 quat_to_matrix(const quat& q) noexcept {
  mat3 out{};
  ::myb_quat_to_matrix(q.data(), out.data());
  return out;
}
inline quat quat_from_matrix(const mat3& m) noexcept {
  quat out{};
  ::myb_quat_from_matrix(m.data(), out.data());
  return out;
}
inline vec3 quat_rotate(const quat& q, const vec3& v) noexcept {
  vec3 out{};
  ::myb_quat_rotate(q.data(), v.data(), out.data());
  return out;
}
inline quat quat_from_rotation_vector(const vec3& phi) noexcept {
  quat out{};
  ::myb_quat_from_rotation_vector(phi.data(), out.data());
  return out;
}
inline vec3 quat_to_rotation_vector(const quat& q) noexcept {
  vec3 out{};
  ::myb_quat_to_rotation_vector(q.data(), out.data());
  return out;
}
inline quat quat_slerp(const quat& a, const quat& b, double t) noexcept {
  quat out{};
  ::myb_quat_slerp(a.data(), b.data(), t, out.data());
  return out;
}
inline quat quat_derivative(const quat& q, const vec3& omega) noexcept {
  quat out{};
  ::myb_quat_derivative(q.data(), omega.data(), out.data());
  return out;
}

inline result<mat3> euler_to_matrix(const std::array<int, 3>& axes, const vec3& angles) {
  mat3 out{};
  const status code = ::myb_euler_to_matrix(axes.data(), angles.data(), out.data());
  return detail::give(code, out);
}
inline result<vec3> euler_from_matrix(const std::array<int, 3>& axes, const mat3& m) {
  vec3 out{};
  const status code = ::myb_euler_from_matrix(axes.data(), m.data(), out.data());
  return detail::give(code, out);
}

struct principal {
  vec3 moments; // ascending
  mat3 axes;    // rows: the principal axes
};
inline result<principal> principal_axes(const mat3& inertia) {
  principal out{};
  const status code = ::myb_principal_axes(inertia.data(), out.moments.data(), out.axes.data());
  return detail::give(code, out);
}

struct rigid_rates {
  quat dq;
  vec3 dl;
};
inline rigid_rates rigid_body_rates(const vec3& moments, const quat& q, const vec3& l,
                                    const std::optional<vec3>& torque = std::nullopt) noexcept {
  rigid_rates out{};
  ::myb_rigid_body_rates(moments.data(), q.data(), l.data(), torque ? torque->data() : nullptr,
                         out.dq.data(), out.dl.data());
  return out;
}

struct core_mantle {
  vec3 rate_mantle;
  vec3 rate_core;
};
inline result<core_mantle> core_mantle_rates(const mat3& inertia_mantle,
                                             const std::optional<mat3>& inertia_rate,
                                             const std::array<double, 2>& core, double kv,
                                             const vec3& omega_mantle, const vec3& omega_core,
                                             const vec3& torque) {
  core_mantle out{};
  const status code = ::myb_core_mantle_rates(
      inertia_mantle.data(), inertia_rate ? inertia_rate->data() : nullptr, core.data(), kv,
      omega_mantle.data(), omega_core.data(), torque.data(), out.rate_mantle.data(),
      out.rate_core.data());
  return detail::give(code, out);
}

inline vec3 torque_point_mass(double gm, const vec3& r, const vec3& moments) noexcept {
  vec3 out{};
  ::myb_torque_point_mass(gm, r.data(), moments.data(), out.data());
  return out;
}
inline vec3 torque_field(double mass_ratio, const vec3& r, const vec3& a) noexcept {
  vec3 out{};
  ::myb_torque_field(mass_ratio, r.data(), a.data(), out.data());
  return out;
}
inline result<vec3> torque_figure_figure(double gm_a, const mat3& gq_a, const mat3& gq_b,
                                         const vec3& r) {
  vec3 out{};
  const status code =
      ::myb_torque_figure_figure(gm_a, gq_a.data(), gq_b.data(), r.data(), out.data());
  return detail::give(code, out);
}

struct attitude_state {
  quat q; // the attitude
  vec3 l; // the body angular momentum (principal axes)
};

// Torque-free rotation in closed form.
class free_rotation {
 public:
  static result<free_rotation> init(const vec3& moments, const quat& q, const vec3& l) {
    ::myb_free_rotation made{};
    const status code = ::myb_free_rotation_init(moments.data(), q.data(), l.data(), &made);
    if (code != MYB_OK) {
      return code;
    }
    return free_rotation(made);
  }
  attitude_state at(double t) const noexcept {
    attitude_state out{};
    ::myb_free_rotation_state(&rotation_, t, out.q.data(), out.l.data());
    return out;
  }
  const ::myb_free_rotation& get() const noexcept { return rotation_; }

 private:
  explicit free_rotation(const ::myb_free_rotation& made) noexcept : rotation_(made) {}
  ::myb_free_rotation rotation_;
};

// The torque on a body at t in the attitude q (body frame), for the splitting integrator. It
// must not throw.
using torque_function = std::function<status(double t, const quat& q, vec3& torque)>;

namespace detail {
inline status call_torque(void* context, double t, const double q[4], double torque[3]) noexcept {
  const auto& function = *static_cast<const torque_function*>(context);
  const quat attitude{q[0], q[1], q[2], q[3]};
  vec3 out{};
  const status code = function(t, attitude, out);
  if (code != MYB_OK) {
    return code;
  }
  copy_out(out, torque);
  return MYB_OK;
}
static_assert(noexcept(call_torque(nullptr, 0.0, nullptr, nullptr)),
              "a function called from C must not let an exception through");
} // namespace detail

// (q, l) advanced from t0 by `steps` steps of h (orders 2, 4, 6, 8).
inline result<attitude_state> rotation_split(const vec3& moments, const torque_function& torque,
                                             int order, double t0, double h, std::int64_t steps,
                                             const attitude_state& start) {
  attitude_state out = start;
  const status code =
      ::myb_rotation_split(moments.data(), &detail::call_torque, detail::context_of(&torque),
                           order, t0, h, steps, out.q.data(), out.l.data());
  return detail::give(code, out);
}

// ---- Observables ------------------------------------------------------------------------------

// The light-time model (myb_light_model): the speed of light, PPN gamma, the bodies that delay
// light (each a GM and a provider, kept here).
class light_model {
 public:
  light_model(double c, double gamma) noexcept : c_(c), gamma_(gamma) {}

  void add_body(double gm, const provider& body) {
    providers_.push_back(body);
    ::myb_light_body made{};
    made.gm = gm;
    made.state = body.function();
    made.context = body.context();
    bodies_.push_back(made);
  }

  ::myb_light_model model() const noexcept {
    ::myb_light_model out{};
    out.c = c_;
    out.gamma = gamma_;
    out.bodies = detail::data_or_null(bodies_);
    out.count = bodies_.size();
    return out;
  }

 private:
  double c_;
  double gamma_;
  std::vector<provider> providers_;
  std::vector<::myb_light_body> bodies_;
};

inline result<::myb_light_path> light_time(const light_model& model, const provider& target,
                                           ::myb_time t, const state& observer,
                                           ::myb_light_direction direction = MYB_LIGHT_RECEPTION) {
  const ::myb_light_model made = model.model();
  ::myb_light_path out{};
  const status code = ::myb_light_time(&made, target.function(), target.context(), t,
                                       observer.data(), direction, &out);
  return detail::give(code, out);
}

inline double shapiro_delay(double gm, double c, double gamma, const vec3& x1, const vec3& x2,
                            const vec3& xb) noexcept {
  return ::myb_shapiro_delay(gm, c, gamma, x1.data(), x2.data(), xb.data());
}

struct round_trip_times {
  double total;         // s
  ::myb_time t_bounce;
  ::myb_time t_transmit;
};
inline result<round_trip_times> round_trip(const light_model& model, const provider& station,
                                           const provider& target, ::myb_time t) {
  const ::myb_light_model made = model.model();
  round_trip_times out{};
  const status code =
      ::myb_round_trip(&made, station.function(), station.context(), target.function(),
                       target.context(), t, &out.total, &out.t_bounce, &out.t_transmit);
  return detail::give(code, out);
}

inline vec3 aberration(const vec3& p, const vec3& v, double c) noexcept {
  vec3 out{};
  ::myb_aberration(p.data(), v.data(), c, out.data());
  return out;
}
inline vec3 light_deflection(double gm, double c, double gamma, const vec3& p, const vec3& q,
                             const vec3& e, double em, double dlim) noexcept {
  vec3 out{};
  ::myb_light_deflection(gm, c, gamma, p.data(), q.data(), e.data(), em, dlim, out.data());
  return out;
}
struct radec {
  double ra;
  double dec;
};
inline radec direction_radec(const vec3& p) noexcept {
  radec out{};
  ::myb_direction_radec(p.data(), &out.ra, &out.dec);
  return out;
}
inline double plasma_group_delay(double column, double frequency, double e, double m_e,
                                 double eps0) noexcept {
  return ::myb_plasma_group_delay(column, frequency, e, m_e, eps0);
}
inline result<double> plasma_column(const ::myb_plasma_profile& profile, const vec3& p1,
                                    const vec3& p2) {
  double column = 0.0;
  const status code = ::myb_plasma_column(&profile, p1.data(), p2.data(), &column);
  return detail::give(code, column);
}

// ---- Files ------------------------------------------------------------------------------------------

// A table of numbers (CSV and the like, myb/io.h).
class table {
 public:
  static result<table> parse(const std::string& text) {
    ::myb_table* made = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_table_load(text.data(), text.size(), nullptr, &made, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return table(made);
  }
  static result<table> load_file(const std::string& path) {
    ::myb_table* made = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_table_load_file(path.c_str(), nullptr, &made, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return table(made);
  }
  std::size_t rows() const noexcept { return ::myb_table_rows(table_.get()); }
  std::size_t columns() const noexcept { return ::myb_table_columns(table_.get()); }
  std::string name(std::size_t column) const { return ::myb_table_name(table_.get(), column); }
  std::vector<double> row(std::size_t index) const {
    const double* values = ::myb_table_row(table_.get(), index);
    return std::vector<double>(values, values + columns());
  }
  ::myb_table* get() const noexcept { return table_.get(); }

 private:
  explicit table(::myb_table* made) noexcept : table_(made) {}
  std::unique_ptr<::myb_table, detail::table_deleter> table_;
};

// Rows of `columns` values (row by row) as comma-separated text; names empty: no row of names.
inline result<std::string> table_write(const std::vector<std::string>& names, std::size_t columns,
                                       const std::vector<double>& values) {
  if ((!names.empty() && names.size() != columns) ||
      (columns != 0 && values.size() % columns != 0)) {
    return MYB_E_DOMAIN;
  }
  std::vector<const char*> pointers(names.size());
  for (std::size_t k = 0; k < names.size(); ++k) {
    pointers[k] = names[k].c_str();
  }
  const std::size_t rows = columns != 0 ? values.size() / columns : 0;
  return detail::text([&](char* buffer, std::size_t size, std::size_t* length) {
    return ::myb_table_write(pointers.empty() ? nullptr : pointers.data(), columns,
                             detail::data_or_null(values), rows, buffer, size, length);
  });
}

// A writer of SPK files (myb/kin.h): types 9 and 13 from states, 2 and 3 from a provider, any
// type from its words.
class spk_writer {
 public:
  static result<spk_writer> create(const std::string& internal_name = "myb") {
    ::myb_spk_writer* made = nullptr;
    const status code = ::myb_spk_writer_create(internal_name.c_str(), nullptr, &made);
    if (code != MYB_OK) {
      return code;
    }
    return spk_writer(made);
  }

  result<void> add_states(const ::myb_spk_segment_header& header, int type, int degree,
                          const std::vector<double>& epochs, const std::vector<state>& states) {
    if (epochs.size() != states.size()) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(::myb_spk_writer_add_states(writer_.get(), &header, type, degree,
                                                       detail::data_or_null(epochs),
                                                       detail::flat(states), states.size()));
  }
  result<void> add_chebyshev(const ::myb_spk_segment_header& header, int type,
                             const provider& source, double start, double stop,
                             std::size_t intervals, int degree) {
    return detail::outcome(::myb_spk_writer_add_chebyshev(writer_.get(), &header, type,
                                                          source.function(), source.context(),
                                                          start, stop, intervals, degree,
                                                          nullptr));
  }
  // The same, with the largest distance (m) between the provider and the polynomials.
  result<double> add_chebyshev_measured(const ::myb_spk_segment_header& header, int type,
                                        const provider& source, double start, double stop,
                                        std::size_t intervals, int degree) {
    double max_error = 0.0;
    const status code = ::myb_spk_writer_add_chebyshev(writer_.get(), &header, type,
                                                       source.function(), source.context(),
                                                       start, stop, intervals, degree,
                                                       &max_error);
    return detail::give(code, max_error);
  }
  result<void> add_raw(const ::myb_spk_segment_header& header, int type, double start,
                       double stop, const std::vector<double>& words) {
    return detail::outcome(::myb_spk_writer_add_raw(writer_.get(), &header, type, start, stop,
                                                    detail::data_or_null(words), words.size()));
  }

  std::size_t size() const noexcept { return ::myb_spk_writer_size(writer_.get()); }
  result<std::vector<unsigned char>> write() const {
    return detail::bytes(size(), [&](void* buffer, std::size_t length) {
      return ::myb_spk_writer_write(writer_.get(), buffer, length);
    });
  }
  result<void> write_file(const std::string& path) const {
    return detail::outcome(::myb_spk_writer_write_file(writer_.get(), path.c_str()));
  }
  ::myb_spk_writer* get() const noexcept { return writer_.get(); }

 private:
  explicit spk_writer(::myb_spk_writer* made) noexcept : writer_(made) {}
  std::unique_ptr<::myb_spk_writer, detail::spk_writer_deleter> writer_;
};

// SP3 precise orbits (myb/kin.h).
class sp3 {
 public:
  static result<sp3> parse(const std::string& text) {
    ::myb_sp3* made = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_sp3_load(text.data(), text.size(), nullptr, &made, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return sp3(made);
  }
  static result<sp3> load_file(const std::string& path) {
    ::myb_sp3* made = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_sp3_load_file(path.c_str(), nullptr, &made, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return sp3(made);
  }

  std::size_t satellite_count() const noexcept { return ::myb_sp3_satellite_count(sp3_.get()); }
  std::string satellite(std::size_t index) const { return ::myb_sp3_satellite(sp3_.get(), index); }
  std::size_t epoch_count() const noexcept { return ::myb_sp3_epoch_count(sp3_.get()); }
  std::string time_system() const { return ::myb_sp3_time_system(sp3_.get()); }
  std::string frame() const { return ::myb_sp3_frame(sp3_.get()); }
  ::myb_calendar epoch(std::size_t index) const noexcept {
    ::myb_calendar out{};
    ::myb_sp3_epoch(sp3_.get(), index, &out);
    return out;
  }

  struct record {
    int given; // 1: the position, 2: the velocity
    vec3 position;
    vec3 velocity;
  };
  record at(std::size_t sat, std::size_t epoch_index) const noexcept {
    record out{};
    out.given = ::myb_sp3_state(sp3_.get(), sat, epoch_index, out.position.data(),
                                out.velocity.data());
    return out;
  }
  struct clock_record {
    int given; // 1: the clock, 2: its rate
    double correction; // s
    double rate;       // s/s
  };
  clock_record clock(std::size_t sat, std::size_t epoch_index) const noexcept {
    clock_record out{};
    out.given = ::myb_sp3_clock(sp3_.get(), sat, epoch_index, &out.correction, &out.rate);
    return out;
  }
  std::array<double, 8> sigma(std::size_t sat, std::size_t epoch_index) const noexcept {
    std::array<double, 8> out{};
    ::myb_sp3_sigma(sp3_.get(), sat, epoch_index, out.data());
    return out;
  }
  double accuracy(std::size_t sat) const noexcept { return ::myb_sp3_accuracy(sp3_.get(), sat); }
  unsigned flags(std::size_t sat, std::size_t epoch_index) const noexcept {
    return ::myb_sp3_flags(sp3_.get(), sat, epoch_index);
  }
  ::myb_sp3* get() const noexcept { return sp3_.get(); }

 private:
  explicit sp3(::myb_sp3* made) noexcept : sp3_(made) {}
  std::unique_ptr<::myb_sp3, detail::sp3_deleter> sp3_;
};

inline result<std::string> sp3_write(const ::myb_sp3_data& data) {
  return detail::text([&](char* buffer, std::size_t size, std::size_t* length) {
    return ::myb_sp3_write(&data, buffer, size, length);
  });
}

// A spacecraft clock of type 1 (myb/ck.h).
class sclk {
 public:
  static result<sclk> create(const kernel_pool& pool, std::int64_t id) {
    ::myb_sclk* made = nullptr;
    ::myb_diag diag{};
    const status code = ::myb_sclk_create(pool.get(), id, nullptr, &made, &diag);
    if (code != MYB_OK) {
      return {code, detail::describe(diag)};
    }
    return sclk(made);
  }
  int time_system() const noexcept { return ::myb_sclk_time_system(sclk_.get()); }
  std::size_t partitions() const noexcept { return ::myb_sclk_partitions(sclk_.get()); }
  result<::myb_time> to_time(double ticks) const {
    ::myb_time t{};
    const status code = ::myb_sclk_to_time(sclk_.get(), ticks, &t);
    return detail::give(code, t);
  }
  result<double> from_time(::myb_time t) const {
    double ticks = 0.0;
    const status code = ::myb_sclk_from_time(sclk_.get(), t, &ticks);
    return detail::give(code, ticks);
  }
  result<double> encode(const std::string& text) const {
    double ticks = 0.0;
    const status code = ::myb_sclk_encode(sclk_.get(), text.data(), text.size(), &ticks);
    return detail::give(code, ticks);
  }
  result<std::string> decode(double ticks) const {
    return detail::text([&](char* buffer, std::size_t size, std::size_t* length) {
      return ::myb_sclk_decode(sclk_.get(), ticks, buffer, size, length);
    });
  }
  ::myb_sclk* get() const noexcept { return sclk_.get(); }

 private:
  explicit sclk(::myb_sclk* made) noexcept : sclk_(made) {}
  std::unique_ptr<::myb_sclk, detail::sclk_deleter> sclk_;
};

// A set of CK files (myb/ck.h).
class ck {
 public:
  static result<ck> create() {
    ::myb_ck* made = nullptr;
    const status code = ::myb_ck_create(nullptr, &made);
    if (code != MYB_OK) {
      return code;
    }
    return ck(made);
  }
  result<void> load_file(const std::string& path) {
    ::myb_diag diag{};
    return detail::outcome(::myb_ck_load_file(ck_.get(), path.c_str(), &diag), diag);
  }
  // An image in memory (what a CK writer gives), kept by the set.
  result<void> load_memory(std::vector<unsigned char> image) {
    auto kept = std::make_shared<std::vector<unsigned char>>(std::move(image));
    ::myb_diag diag{};
    const status code = ::myb_ck_load_memory(ck_.get(), kept->data(), kept->size(), &diag);
    if (code == MYB_OK) {
      images_.push_back(std::move(kept));
    }
    return detail::outcome(code, diag);
  }
  std::size_t segment_count() const noexcept { return ::myb_ck_segment_count(ck_.get()); }
  ::myb_ck_segment_info segment(std::size_t index) const noexcept {
    ::myb_ck_segment_info info{};
    ::myb_ck_segment(ck_.get(), index, &info);
    return info;
  }
  // The orientation of an instrument at `ticks` of its clock, or the nearest within tolerance.
  result<::myb_ck_pointing> pointing_at(std::int32_t instrument, double ticks, double tolerance,
                                        bool need_rate = false) const {
    ::myb_ck_pointing out{};
    const status code =
        ::myb_ck_pointing_at(ck_.get(), instrument, ticks, tolerance, need_rate ? 1 : 0, &out);
    return detail::give(code, out);
  }
  ::myb_ck* get() const noexcept { return ck_.get(); }

 private:
  explicit ck(::myb_ck* made) noexcept : ck_(made) {}
  std::vector<std::shared_ptr<std::vector<unsigned char>>> images_;
  std::unique_ptr<::myb_ck, detail::ck_deleter> ck_;
};

// An instrument frame as a provider (myb_orientation_ck) that keeps the set and the clock.
inline orientation_provider ck_orientation(std::shared_ptr<const ck> files,
                                           std::shared_ptr<const sclk> clock,
                                           std::int32_t instrument, std::int32_t base = 0,
                                           double tolerance = 0.0) {
  struct holder {
    std::shared_ptr<const ck> files;
    std::shared_ptr<const sclk> clock;
    ::myb_ck_orientation context;
  };
  auto kept = std::make_shared<holder>();
  kept->context = ::myb_ck_orientation{files->get(), clock->get(), instrument, base, tolerance};
  kept->files = std::move(files);
  kept->clock = std::move(clock);
  ::myb_ck_orientation* raw = &kept->context;
  return orientation_provider(::myb_orientation_ck, raw, std::move(kept));
}

// A writer of CK files (myb/ck.h): types 1 and 3 (pointings), 2 (intervals), 4 (Chebyshev),
// 5 (packets) and 6 (mini-segments).
class ck_writer {
 public:
  static result<ck_writer> create(const std::string& internal_name = "myb") {
    ::myb_ck_writer* made = nullptr;
    const status code = ::myb_ck_writer_create(internal_name.c_str(), nullptr, &made);
    if (code != MYB_OK) {
      return code;
    }
    return ck_writer(made);
  }

  // rates empty: none; interval_starts empty: one interval of them all.
  result<void> add_pointings(const ::myb_ck_segment_header& header, int type,
                             const std::vector<double>& ticks, const std::vector<quat>& quaternions,
                             const std::vector<vec3>& rates = {},
                             const std::vector<double>& interval_starts = {}) {
    if (quaternions.size() != ticks.size() || (!rates.empty() && rates.size() != ticks.size())) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(::myb_ck_writer_add_pointings(
        writer_.get(), &header, type, detail::data_or_null(ticks), detail::flat(quaternions),
        detail::flat(rates), ticks.size(), detail::data_or_null(interval_starts),
        interval_starts.size()));
  }
  result<void> add_intervals(const ::myb_ck_segment_header& header,
                             const std::vector<double>& starts, const std::vector<double>& stops,
                             const std::vector<quat>& quaternions, const std::vector<vec3>& rates,
                             const std::vector<double>& seconds_per_tick) {
    const std::size_t count = starts.size();
    if (stops.size() != count || quaternions.size() != count || rates.size() != count ||
        seconds_per_tick.size() != count) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(::myb_ck_writer_add_intervals(
        writer_.get(), &header, detail::data_or_null(starts), detail::data_or_null(stops),
        detail::flat(quaternions), detail::flat(rates), detail::data_or_null(seconds_per_tick),
        count));
  }
  // counts: 7 per interval; coefficients: the series one after another.
  result<void> add_chebyshev(const ::myb_ck_segment_header& header,
                             const std::vector<double>& starts, const std::vector<double>& stops,
                             const std::vector<int>& counts,
                             const std::vector<double>& coefficients, bool has_rates) {
    const std::size_t count = starts.size();
    if (stops.size() != count || counts.size() != 7 * count) {
      return MYB_E_DOMAIN;
    }
    std::size_t needed = 0;
    for (const int n : counts) {
      needed += n > 0 ? static_cast<std::size_t>(n) : 0;
    }
    if (coefficients.size() < needed) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(::myb_ck_writer_add_chebyshev(
        writer_.get(), &header, detail::data_or_null(starts), detail::data_or_null(stops),
        detail::data_or_null(counts), detail::data_or_null(coefficients), count,
        has_rates ? 1 : 0));
  }
  result<void> add_packets(const ::myb_ck_segment_header& header, const ::myb_ck_packets& packets,
                           bool has_rates, double start, double stop,
                           const std::vector<double>& interval_starts = {}) {
    return detail::outcome(::myb_ck_writer_add_packets(
        writer_.get(), &header, &packets, has_rates ? 1 : 0, start, stop,
        detail::data_or_null(interval_starts), interval_starts.size()));
  }
  result<void> add_minisegments(const ::myb_ck_segment_header& header,
                                const std::vector<::myb_ck_packets>& minisegments,
                                const std::vector<double>& boundaries, bool has_rates,
                                bool select_last, double start, double stop) {
    if (boundaries.size() != minisegments.size() + 1) {
      return MYB_E_DOMAIN;
    }
    return detail::outcome(::myb_ck_writer_add_minisegments(
        writer_.get(), &header, detail::data_or_null(minisegments), boundaries.data(),
        minisegments.size(), has_rates ? 1 : 0, select_last ? 1 : 0, start, stop));
  }

  std::size_t size() const noexcept { return ::myb_ck_writer_size(writer_.get()); }
  result<std::vector<unsigned char>> write() const {
    return detail::bytes(size(), [&](void* buffer, std::size_t length) {
      return ::myb_ck_writer_write(writer_.get(), buffer, length);
    });
  }
  result<void> write_file(const std::string& path) const {
    return detail::outcome(::myb_ck_writer_write_file(writer_.get(), path.c_str()));
  }
  ::myb_ck_writer* get() const noexcept { return writer_.get(); }

 private:
  explicit ck_writer(::myb_ck_writer* made) noexcept : writer_(made) {}
  std::unique_ptr<::myb_ck_writer, detail::ck_writer_deleter> writer_;
};

// ---- Two-line elements (ext/sgp4: include myb_sgp4.h before this header, link myb::sgp4) -----
#if defined(MYB_SGP4_H_)

inline result<::myb_tle> tle_parse(const std::string& line1, const std::string& line2,
                                   unsigned flags = 0) {
  ::myb_tle out{};
  ::myb_diag diag{};
  const status code = ::myb_tle_parse(line1.c_str(), line2.c_str(), flags, &out, &diag);
  return detail::give(code, diag, out);
}
inline ::myb_calendar tle_epoch(const ::myb_tle& tle) noexcept {
  ::myb_calendar out{};
  ::myb_tle_epoch(&tle, &out);
  return out;
}
struct tle_lines {
  std::string line1;
  std::string line2;
};
inline result<tle_lines> tle_format(const ::myb_tle& tle) {
  char line1[70] = {};
  char line2[70] = {};
  const status code = ::myb_tle_format(&tle, line1, line2);
  return detail::give(code, tle_lines{line1, line2});
}
inline result<::myb_tle> tle_from_omm(const ::myb_omm& omm) {
  ::myb_tle out{};
  const status code = ::myb_tle_from_omm(&omm, &out);
  return detail::give(code, out);
}
inline result<::myb_omm> omm_from_tle(const ::myb_tle& tle) {
  ::myb_omm out{};
  const status code = ::myb_omm_from_tle(&tle, &out);
  return detail::give(code, out);
}

// The theory of an element set: TEME states (m, m/s) at minutes from its epoch.
class sgp4 {
 public:
  static result<sgp4> init(const ::myb_tle& tle, ::myb_sgp4_gravity gravity = MYB_SGP4_WGS72,
                           ::myb_sgp4_mode mode = MYB_SGP4_IMPROVED) {
    ::myb_sgp4 made{};
    const status code = ::myb_sgp4_init(&tle, gravity, mode, &made);
    if (code != MYB_OK) {
      return code;
    }
    return sgp4(made);
  }
  result<state> at(double minutes) const {
    state y{};
    const status code = ::myb_sgp4_propagate(&record_, minutes, y.data(), y.data() + 3, nullptr);
    return detail::give(code, y);
  }
  const ::myb_sgp4& get() const noexcept { return record_; }

 private:
  explicit sgp4(const ::myb_sgp4& made) noexcept : record_(made) {}
  ::myb_sgp4 record_;
};

struct tle_fit_report {
  double rms; // m
  int iterations;
};
// tle: the first guess, the fitted set on return (on MYB_E_NOCONV the best so far).
inline result<tle_fit_report> tle_fit(const std::vector<double>& minutes,
                                      const std::vector<vec3>& positions, ::myb_tle& tle,
                                      const ::myb_tle_fit_options* options = nullptr) {
  if (positions.size() != minutes.size()) {
    return MYB_E_DOMAIN;
  }
  tle_fit_report out{};
  const status code = ::myb_tle_fit(detail::data_or_null(minutes), detail::flat(positions),
                                    minutes.size(), options, &tle, &out.rms, &out.iterations);
  return detail::give(code, out);
}
inline result<void> tle_from_state(const vec3& r, const vec3& v, ::myb_sgp4_gravity gravity,
                                   ::myb_tle& tle) {
  return detail::outcome(::myb_tle_from_state(r.data(), v.data(), gravity, &tle));
}

#endif // MYB_SGP4_H_

} // namespace myb

#endif // MYB_MYB_HPP_
