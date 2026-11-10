# SPDX-License-Identifier: BSD-3-Clause
"""myb for Python: the motion of objects in space.

Two layers over the shared library of myb:

  * `myb.lib` and `myb.ffi` (cffi) reach every function and structure of the C interface,
    by the names of the C headers: `myb.lib.myb_kepler_propagate(...)`. The declarations
    are generated from the headers, so nothing is missing and nothing can differ.
  * The classes and functions of this module wrap the path most programs take — kernels
    and time, ephemerides, the environment and the force model, propagation, events, the
    text messages — with ownership (objects are released when Python drops them), errors
    as exceptions (`myb.Error` carries the status of the library) and Python callables for
    terms, providers and events.

No arithmetic is done here: a number obtained through this package has the bits of the same
calls made in C. States are sequences of six floats (position in m, velocity in m/s); times
of a model are seconds from the epoch of its environment.

The optional libraries of ext/ are loaded on request: `myb.extension("sgp4")`.
"""
from collections import namedtuple

from ._native import extension, ffi, lib, library_path

__all__ = [
    "Elements", "Environment", "Ephemeris", "Error", "ForceModel", "GravityField",
    "KernelPool", "Orbit", "TimeSystem", "elements_from_state", "extension", "ffi",
    "kepler_propagate", "lib", "library_path", "read_oem", "seconds", "state_from_elements",
    "time_from_seconds", "version", "write_oem",
]

version = (lib.MYB_VERSION_MAJOR, lib.MYB_VERSION_MINOR, lib.MYB_VERSION_PATCH,
           lib.MYB_VERSION_TWEAK)
__version__ = ".".join(str(part) for part in version)

Elements = namedtuple("Elements", "p e i raan argp nu")


class Error(Exception):
    """A call of the library failed: `status` is its code (lib.MYB_E_DOMAIN, ...), `name`
    the name of the code, `message` what a parser said (with the line), if anything."""

    def __init__(self, status, message=""):
        self.status = int(status)
        self.name = ffi.string(lib.myb_status_name(status)).decode("ascii")
        self.message = message
        super().__init__("%s: %s" % (self.name, message) if message else self.name)


def _check(status, diag=None):
    if status == lib.MYB_OK:
        return
    message = ""
    if diag is not None:
        message = ffi.string(diag.message).decode("utf-8", "replace")
        if diag.line > 0:
            message = "line %d: %s" % (diag.line, message)
    raise Error(status, message)


def _doubles(values, count, what):
    values = [float(v) for v in values]
    if len(values) != count:
        raise ValueError("%s: %d numbers expected, %d given" % (what, count, len(values)))
    return ffi.new("double[]", values)


def _calendar(label):
    """(year, month, day, hour, minute, second) -> myb_calendar."""
    year, month, day, hour, minute, second = label
    return ffi.new("myb_calendar*", {"year": int(year), "month": int(month), "day": int(day),
                                     "hour": int(hour), "minute": int(minute),
                                     "second": float(second)})


def _label(calendar):
    return (calendar.year, calendar.month, calendar.day, calendar.hour, calendar.minute,
            calendar.second)


class _Owner:
    """An object of the library: released when dropped or closed, usable in `with`."""

    _destroy = None

    def __init__(self, pointer):
        self._pointer = ffi.gc(pointer, self._destroy)
        self._kept = []  # what the object points to must live as long as it does

    @property
    def pointer(self):
        """The object of the C interface, for the functions of `myb.lib`."""
        if self._pointer is None:
            raise ValueError("%s is closed" % type(self).__name__)
        return self._pointer

    def close(self):
        if self._pointer is not None:
            ffi.release(self._pointer)
            self._pointer = None
            self._kept = []

    def __enter__(self):
        return self

    def __exit__(self, *unused):
        self.close()
        return False


# ---- Time ---------------------------------------------------------------------------------------

def time_from_seconds(seconds_past_epoch):
    """An instant (myb_time) from seconds past J2000 of its time scale."""
    return lib.myb_time_from_seconds(float(seconds_past_epoch))


def seconds(instant):
    """The seconds past J2000 of an instant, as one float."""
    return lib.myb_time_seconds(instant)


class KernelPool(_Owner):
    """A pool of text kernels (myb/io.h)."""

    _destroy = staticmethod(lib.myb_kpool_destroy)

    def __init__(self, *paths):
        out = ffi.new("myb_kpool**")
        _check(lib.myb_kpool_create(ffi.NULL, out))
        super().__init__(out[0])
        for path in paths:
            self.load_file(path)

    def load_file(self, path):
        diag = ffi.new("myb_diag*")
        _check(lib.myb_kpool_load_file(self.pointer, str(path).encode("utf-8"), diag), diag)


class TimeSystem(_Owner):
    """The time scales of a pool (myb/time.h). Scales are named: "TT", "TAI", "TDB", ...;
    instants are myb_time (see time_from_seconds, seconds); UTC is a calendar label
    (year, month, day, hour, minute, second)."""

    _destroy = staticmethod(lib.myb_timesys_destroy)

    def __init__(self, pool):
        out = ffi.new("myb_timesys**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_timesys_create(pool.pointer, ffi.NULL, out, diag), diag)
        super().__init__(out[0])

    def _scale(self, name):
        scale = ffi.new("myb_scale*")
        _check(lib.myb_timesys_scale(self.pointer, name.encode("ascii"), scale))
        return scale[0]

    def convert(self, instant, from_scale, to_scale):
        out = ffi.new("myb_time*")
        _check(lib.myb_timesys_convert(self.pointer, instant, self._scale(from_scale),
                                       self._scale(to_scale), out))
        return out[0]

    def from_utc(self, label, to_scale):
        out = ffi.new("myb_time*")
        _check(lib.myb_timesys_from_utc(self.pointer, _calendar(label), self._scale(to_scale),
                                        out))
        return out[0]

    def to_utc(self, instant, from_scale):
        out = ffi.new("myb_calendar*")
        _check(lib.myb_timesys_to_utc(self.pointer, instant, self._scale(from_scale), out))
        return _label(out)


# ---- Two bodies ---------------------------------------------------------------------------------

def _pv(state):
    y = [float(v) for v in state]
    if len(y) != 6:
        raise ValueError("a state is six numbers, %d given" % len(y))
    return ffi.new("myb_pv*", {"r": y[:3], "v": y[3:]})


def _state_of(pv):
    return (pv.r[0], pv.r[1], pv.r[2], pv.v[0], pv.v[1], pv.v[2])


def state_from_elements(mu, elements):
    """The state of the classical elements (p, e, i, raan, argp, nu): m and rad."""
    given = ffi.new("myb_elements*", dict(zip(Elements._fields, (float(v) for v in elements))))
    pv = ffi.new("myb_pv*")
    _check(lib.myb_pv_from_elements(float(mu), given, pv))
    return _state_of(pv)


def elements_from_state(mu, state):
    out = ffi.new("myb_elements*")
    _check(lib.myb_elements_from_pv(float(mu), _pv(state), out))
    return Elements(out.p, out.e, out.i, out.raan, out.argp, out.nu)


def kepler_propagate(mu, state, dt):
    """The state a time dt later on the orbit of two bodies."""
    out = ffi.new("myb_pv*")
    _check(lib.myb_kepler_propagate(float(mu), _pv(state), float(dt), out))
    return _state_of(out)


# ---- Ephemerides and fields ---------------------------------------------------------------------

class Ephemeris(_Owner):
    """A set of SPK files (myb/kin.h)."""

    _destroy = staticmethod(lib.myb_spk_destroy)

    def __init__(self, *paths):
        out = ffi.new("myb_spk**")
        _check(lib.myb_spk_create(ffi.NULL, out))
        super().__init__(out[0])
        for path in paths:
            self.load_file(path)

    def load_file(self, path):
        diag = ffi.new("myb_diag*")
        _check(lib.myb_spk_load_file(self.pointer, str(path).encode("utf-8"), diag), diag)

    def state(self, target, center, instant):
        """Target relative to center at an instant of the files' time scale."""
        y = ffi.new("double[6]")
        _check(lib.myb_spk_state(self.pointer, int(target), int(center), instant, y, ffi.NULL))
        return tuple(y)

    def segments(self):
        """The segments of the set: target, center, frame, type, start, stop, name."""
        info = ffi.new("myb_spk_segment_info*")
        found = []
        for index in range(lib.myb_spk_segment_count(self.pointer)):
            lib.myb_spk_segment(self.pointer, index, info)
            found.append({"target": info.target, "center": info.center, "frame": info.frame,
                          "type": info.type, "start": info.start, "stop": info.stop,
                          "name": ffi.string(info.name).decode("ascii", "replace")})
        return found

    def body(self, target, center, frame=0):
        """A body of the set as the state of a body of an Environment."""
        return _SpkBody(self, target, center, frame)


class _SpkBody:
    def __init__(self, ephemeris, target, center, frame):
        self.ephemeris = ephemeris
        self.context = ffi.new("myb_spk_body*", {"spk": ephemeris.pointer, "target": int(target),
                                                 "center": int(center), "frame": int(frame)})


class GravityField(_Owner):
    """The field of a body in spherical harmonics (myb/grav.h)."""

    _destroy = staticmethod(lib.myb_gravity_field_destroy)

    def __init__(self, gm, radius, degree, c, s):
        """c, s: fully normalised coefficients at the index n (n + 1) / 2 + m."""
        count = (degree + 1) * (degree + 2) // 2
        out = ffi.new("myb_gravity_field**")
        _check(lib.myb_gravity_field_create(float(gm), float(radius), int(degree),
                                            _doubles(c, count, "c"), _doubles(s, count, "s"),
                                            ffi.NULL, out))
        super().__init__(out[0])

    @classmethod
    def load_icgem_file(cls, path, max_degree=-1):
        out = ffi.new("myb_gravity_field**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_gravity_field_load_icgem_file(str(path).encode("utf-8"), int(max_degree),
                                                     ffi.NULL, out, diag), diag)
        field = cls.__new__(cls)
        _Owner.__init__(field, out[0])
        return field

    @property
    def degree(self):
        return lib.myb_gravity_field_degree(self.pointer)

    @property
    def gm(self):
        return lib.myb_gravity_field_gm(self.pointer)

    @property
    def radius(self):
        return lib.myb_gravity_field_radius(self.pointer)


# ---- The environment ----------------------------------------------------------------------------

class _Calls:
    """Python callables behind the C callbacks of an object. An exception raised by one
    cannot cross the C code: it is kept, the callback reports a failure, and the call of the
    library that was going on raises it again."""

    def __init__(self):
        self.pending = None

    def failed(self, error):
        if self.pending is None:
            self.pending = error
        return error.status if isinstance(error, Error) else lib.MYB_E_DOMAIN

    def check(self, status, diag=None):
        pending, self.pending = self.pending, None
        if pending is not None:
            raise pending
        _check(status, diag)


class Environment(_Owner):
    """The bodies of a scenario (myb/env.h): what the terms of a force model ask for.

    epoch: the instant (myb_time) of the model time 0; timesys and scale: the time system
    and the scale of the model time, needed when providers work on other scales."""

    _destroy = staticmethod(lib.myb_env_destroy)

    def __init__(self, epoch=None, timesys=None, scale=None, table_step=0.0):
        options = ffi.new("myb_env_options*")
        if epoch is not None:
            options.epoch = epoch
        scale_text = None
        if timesys is not None:
            options.timesys = timesys.pointer
        if scale is not None:
            scale_text = ffi.new("char[]", scale.encode("ascii"))
            options.scale = scale_text
        options.table_step = float(table_step)
        out = ffi.new("myb_env**")
        _check(lib.myb_env_create(options, ffi.NULL, out))
        super().__init__(out[0])
        self._kept += [timesys, scale_text]
        self._calls = _Calls()

    def add_body(self, id, gm, state=None, radius=0.0, luminosity=0.0, limb=(0.0, 0.0)):
        """Adds a body and returns its index. state: None for the body at the origin; a
        body of an Ephemeris (Ephemeris.body); or a callable f(seconds) -> six numbers, the
        state relative to the origin at an instant given in seconds past J2000 of the
        model's scale."""
        body = ffi.new("myb_env_body*")
        body.id = int(id)
        body.gm = float(gm)
        body.radius = float(radius)
        body.luminosity = float(luminosity)
        body.limb = [float(limb[0]), float(limb[1])]
        if isinstance(state, _SpkBody):
            body.state = lib.myb_spk_body_state
            body.state_context = state.context
            self._kept.append(state)
        elif state is not None:
            calls = self._calls

            @ffi.callback("myb_status(void*, myb_time, double[6])")
            def provider(context, instant, out):
                try:
                    values = [float(v) for v in state(lib.myb_time_seconds(instant))]
                    if len(values) != 6:
                        raise ValueError("a state is six numbers, %d given" % len(values))
                    out[0:6] = values
                    return lib.MYB_OK
                except BaseException as error:  # noqa: B902 - nothing may cross the C code
                    return calls.failed(error)

            body.state = provider
            self._kept.append(provider)
        index = ffi.new("size_t*")
        _check(lib.myb_env_add_body(self.pointer, body, index))
        return index[0]

    def find(self, id):
        index = ffi.new("size_t*")
        _check(lib.myb_env_find(self.pointer, int(id), index))
        return index[0]

    def state(self, t, body):
        """The state of a body relative to the origin at the model time t."""
        y = ffi.new("double[6]")
        self._calls.check(lib.myb_env_state(self.pointer, float(t), int(body), y))
        return tuple(y)

    @property
    def body_count(self):
        return lib.myb_env_body_count(self.pointer)


# ---- The force model ----------------------------------------------------------------------------

class ForceModel(_Owner):
    """The forces on one object (myb/model.h) with the ready-made terms of myb/env.h. The
    model keeps alive the environments, fields and structures its terms point to."""

    _destroy = staticmethod(lib.myb_force_model_destroy)

    def __init__(self):
        out = ffi.new("myb_force_model**")
        _check(lib.myb_force_model_create(ffi.NULL, out))
        super().__init__(out[0])
        self._calls = _Calls()

    def _add(self, name, function, jacobian, context, *kept):
        index = ffi.new("size_t*")
        _check(lib.myb_force_model_add(self.pointer, name.encode("ascii"), function, jacobian,
                                       context, index))
        self._kept += [context] + list(kept)
        return index[0]

    def add_raw(self, name, function, jacobian, context, *kept):
        """Any term of the C interface: a function and a Jacobian of `myb.lib` (or
        ffi.NULL) and their context (cdata). The context and whatever else is given are
        kept alive with the model."""
        return self._add(name, function, jacobian, context, *kept)

    def add_point_mass(self, name, gm):
        """The attraction of a point mass at the origin."""
        term = ffi.new("myb_point_mass_term*", {"gm": float(gm)})
        return self._add(name, lib.myb_term_point_mass, lib.myb_term_point_mass_jacobian, term)

    def add_third_body(self, name, env, body):
        """The attraction of a body of the environment away from the origin."""
        term = ffi.new("myb_third_body_term*", {"env": env.pointer, "body": int(body)})
        return self._add(name, lib.myb_term_third_body, lib.myb_term_third_body_jacobian, term,
                         env)

    def add_field(self, name, env, body, field, n_min=0, n_max=-1):
        """The field of a body in harmonics of the degrees n_min .. n_max (n_max < 0: all)."""
        out = ffi.new("myb_gravity_work**")
        _check(lib.myb_gravity_work_create(field.degree, ffi.NULL, out))
        work = ffi.gc(out[0], lib.myb_gravity_work_destroy)
        term = ffi.new("myb_field_term*")
        term.env = env.pointer
        term.body = int(body)
        term.field = field.pointer
        term.work = work
        term.n_min = int(n_min)
        term.n_max = int(n_max)
        return self._add(name, lib.myb_term_field, lib.myb_term_field_jacobian, term, env, field,
                         work)

    @staticmethod
    def _object(mass, area, cr=0.0, cd=0.0):
        return ffi.new("myb_object*", {"mass": float(mass), "area": float(area), "cr": float(cr),
                                       "cd": float(cd)})

    def add_radiation(self, name, env, source, mass, area, cr, c, occulters=()):
        """The pressure of the light of a body on an object of the given mass (kg),
        cross-section (m^2) and coefficient cr, in the shadows of the occulting bodies."""
        occulters = [int(k) for k in occulters]
        if len(occulters) > lib.MYB_ENV_MAX_OCCULTERS:
            raise ValueError("at most %d occulters" % lib.MYB_ENV_MAX_OCCULTERS)
        thing = self._object(mass, area, cr=cr)
        term = ffi.new("myb_radiation_term*")
        term.env = env.pointer
        term.source = int(source)
        term.object = thing
        term.c = float(c)
        term.occulters[0:len(occulters)] = occulters
        term.occulter_count = len(occulters)
        return self._add(name, lib.myb_term_radiation, lib.myb_term_radiation_jacobian, term, env,
                         thing)

    def add_drag(self, name, env, body, mass, area, cd, layers, equatorial_radius, flattening=0.0):
        """The drag of an exponential atmosphere that turns with a body. layers: one
        (base height m, density kg/m^3 at the base, scale height m) or several, by height,
        over the ellipsoid given."""
        layers = [tuple(float(v) for v in layer) for layer in layers]
        if not 1 <= len(layers) <= lib.MYB_ATMOSPHERE_LAYERS:
            raise ValueError("1 .. %d layers" % lib.MYB_ATMOSPHERE_LAYERS)
        air = ffi.new("myb_exponential_atmosphere*")
        air.ellipsoid.a = float(equatorial_radius)
        air.ellipsoid.f = float(flattening)
        air.count = len(layers)
        for k, (base, density, scale) in enumerate(layers):
            air.base[k] = base
            air.density[k] = density
            air.scale[k] = scale
        thing = self._object(mass, area, cd=cd)
        term = ffi.new("myb_drag_term*")
        term.env = env.pointer
        term.body = int(body)
        term.density = lib.myb_density_exponential
        term.density_context = air
        term.object = thing
        return self._add(name, lib.myb_term_drag, lib.myb_term_drag_jacobian, term, env, air, thing)

    def add_schwarzschild(self, name, env, central, c, beta=1.0, gamma=1.0):
        """The post-Newtonian correction of a central body of the environment."""
        term = ffi.new("myb_relativity_term*")
        term.env = env.pointer
        term.ppn = {"c": float(c), "beta": float(beta), "gamma": float(gamma)}
        term.central = int(central)
        term.external = int(central)
        return self._add(name, lib.myb_term_schwarzschild, lib.myb_term_schwarzschild_jacobian,
                         term, env)

    def add_term(self, name, function):
        """A term of the caller's own: function(t, x, v) -> the three components of its
        acceleration at the model time t and the state (x, v), each a tuple of three. Its
        derivatives are taken by differences of the model."""
        calls = self._calls

        @ffi.callback("myb_status(void*, double, const double[3], const double[3], double[3])")
        def term(context, t, x, v, a):
            try:
                values = [float(w) for w in function(t, (x[0], x[1], x[2]), (v[0], v[1], v[2]))]
                if len(values) != 3:
                    raise ValueError("an acceleration is three numbers, %d given" % len(values))
                a[0:3] = values
                return lib.MYB_OK
            except BaseException as error:  # noqa: B902 - nothing may cross the C code
                return calls.failed(error)

        return self._add(name, term, ffi.NULL, ffi.NULL, term, function)

    def acceleration(self, t, x, v):
        a = ffi.new("double[3]")
        self._check(lib.myb_force_model_acceleration(self.pointer, float(t),
                                                     _doubles(x, 3, "x"), _doubles(v, 3, "v"), a))
        return tuple(a)

    def _check(self, status):
        """Raises what a callable of the model, or of one of its environments, raised."""
        for kept in self._kept:
            calls = getattr(kept, "_calls", None)
            if calls is not None and calls.pending is not None:
                self._calls.pending, calls.pending = self._calls.pending or calls.pending, None
        self._calls.check(status)

    @property
    def term_names(self):
        return [ffi.string(lib.myb_force_model_name(self.pointer, k)).decode("ascii")
                for k in range(lib.myb_force_model_count(self.pointer))]

    def set_enabled(self, index, enabled):
        lib.myb_force_model_set_enabled(self.pointer, int(index), 1 if enabled else 0)

    def enabled(self, index):
        return lib.myb_force_model_enabled(self.pointer, int(index)) != 0


# ---- Propagation --------------------------------------------------------------------------------

class Orbit(_Owner):
    """One object under a force model (myb/orbit.h): states at the times asked for, events.

    epsilon: the precision of IAS15 (0: 1e-10); variational: also integrate the transition
    matrix, for partials(); integrator: "ias15" or "dop853" (with rtol, atol)."""

    _destroy = staticmethod(lib.myb_orbit_destroy)
    _integrators = {"ias15": lib.MYB_ORBIT_IAS15, "dop853": lib.MYB_ORBIT_DOP853}

    def __init__(self, model, epsilon=0.0, h_max=0.0, max_steps=0, variational=False,
                 integrator="ias15", rtol=0.0, atol=0.0):
        options = ffi.new("myb_orbit_options*")
        options.epsilon = float(epsilon)
        options.h_max = float(h_max)
        options.max_steps = int(max_steps)
        options.variational = 1 if variational else 0
        options.integrator = self._integrators[integrator]
        options.rtol = float(rtol)
        options.atol = float(atol)
        out = ffi.new("myb_orbit**")
        _check(lib.myb_orbit_create(model.pointer, options, ffi.NULL, out))
        super().__init__(out[0])
        self._model = model
        self._kept.append(model)

    def start(self, t0, state, t_end):
        """The object is at `state` at t0 and will be asked about up to t_end."""
        self._model._check(lib.myb_orbit_start(self.pointer, float(t0), _doubles(state, 6, "state"),
                                               float(t_end)))

    def state(self, t):
        y = ffi.new("double[6]")
        self._model._check(lib.myb_orbit_state(self.pointer, float(t), y))
        return tuple(y)

    def states(self, times):
        """The states at times in the order of the propagation."""
        times = [float(t) for t in times]
        out = ffi.new("double[]", 6 * len(times))
        self._model._check(lib.myb_orbit_states(self.pointer, ffi.new("double[]", times),
                                                len(times), out))
        return [tuple(out[6 * k:6 * k + 6]) for k in range(len(times))]

    def next_event(self, g, t_limit, direction=0):
        """The next zero of g(t, state) up to t_limit: direction > 0 for a rising g, < 0 for
        a falling one, 0 for both. Returns (t, state) or None if there is none."""
        calls = self._model._calls

        @ffi.callback("double(void*, double, const double*)")
        def event(context, t, y):
            try:
                return float(g(t, tuple(y[0:6])))
            except BaseException as error:  # noqa: B902 - nothing may cross the C code
                calls.failed(error)
                return float("nan")

        found = ffi.new("int*")
        t_event = ffi.new("double*")
        y_event = ffi.new("double[6]")
        status = lib.myb_orbit_event(self.pointer, event, ffi.NULL, float(t_limit),
                                     int(direction), found, t_event, y_event)
        self._model._check(status)
        return (t_event[0], tuple(y_event)) if found[0] else None

    def kick(self, dv):
        """An impulse: the velocity changes by dv where the object is now."""
        self._model._check(lib.myb_orbit_kick(self.pointer, _doubles(dv, 3, "dv")))

    def partials(self, t):
        """The state at t and the transition matrix from the start, as six rows."""
        y = ffi.new("double[6]")
        phi = ffi.new("double[36]")
        self._model._check(lib.myb_orbit_partials(self.pointer, float(t), y, phi, ffi.NULL))
        return tuple(y), [tuple(phi[6 * row:6 * row + 6]) for row in range(6)]

    @property
    def time(self):
        return lib.myb_orbit_time(self.pointer)

    @property
    def steps(self):
        return lib.myb_orbit_steps(self.pointer)


# ---- Messages -----------------------------------------------------------------------------------

_notations = {"kvn": lib.MYB_ODM_KVN, "xml": lib.MYB_ODM_XML}


def _text(write):
    """A text of the library's writers: they tell the length when the buffer is too small."""
    length = ffi.new("size_t*")
    status = write(ffi.NULL, 0, length)
    if status not in (lib.MYB_OK, lib.MYB_E_OVERFLOW):
        _check(status)
    buffer = ffi.new("char[]", length[0] + 1)
    _check(write(buffer, length[0] + 1, length))
    return ffi.string(buffer, length[0]).decode("utf-8")


def write_oem(originator, creation_date, segments, notation="kvn"):
    """An orbit ephemeris message of the CCSDS as text. creation_date and the epochs are
    calendar labels (year, month, day, hour, minute, second); each segment is a mapping
    with object_name, object_id, center_name, ref_frame, time_system, epochs and states
    (six numbers each, m and m/s), and optionally interpolation with interpolation_degree."""
    kept = []

    def text(value):
        kept.append(ffi.new("char[]", str(value).encode("utf-8")))
        return kept[-1]

    made = ffi.new("myb_oem_segment_data[]", len(segments))
    for k, segment in enumerate(segments):
        epochs = list(segment["epochs"])
        states = [float(v) for state in segment["states"] for v in state]
        if len(states) != 6 * len(epochs):
            raise ValueError("segment %d: six numbers a state, one state an epoch" % k)
        labels = ffi.new("myb_calendar[]", [_calendar(e)[0] for e in epochs])
        values = ffi.new("double[]", states)
        kept += [labels, values]
        for key in ("object_name", "object_id", "center_name", "ref_frame", "time_system"):
            setattr(made[k], key, text(segment[key]))
        if segment.get("interpolation") is not None:
            made[k].interpolation = text(segment["interpolation"])
            made[k].interpolation_degree = int(segment.get("interpolation_degree", 0))
        made[k].count = len(epochs)
        made[k].epochs = labels
        made[k].states = values
    data = ffi.new("myb_oem_data*")
    data.originator = text(originator)
    data.creation_date = _calendar(creation_date)[0]
    data.segment_count = len(segments)
    data.segments = made
    return _text(lambda buffer, size, length: lib.myb_oem_write_as(
        data, _notations[notation], buffer, size, length))


def read_oem(text):
    """An orbit ephemeris message in either notation: a list of segments, each a mapping
    with "meta" (the keywords of its metadata) and "states" (pairs of a calendar label and
    six numbers); the keywords asked for most are also keys of the segment itself."""
    raw = text.encode("utf-8")
    out = ffi.new("myb_oem**")
    diag = ffi.new("myb_diag*")
    _check(lib.myb_oem_load(raw, len(raw), ffi.NULL, out, diag), diag)
    oem = ffi.gc(out[0], lib.myb_oem_destroy)
    keys = ("OBJECT_NAME", "OBJECT_ID", "CENTER_NAME", "REF_FRAME", "TIME_SYSTEM", "START_TIME",
            "STOP_TIME", "INTERPOLATION", "INTERPOLATION_DEGREE")
    segments = []
    epoch = ffi.new("myb_calendar*")
    y = ffi.new("double[6]")
    for s in range(lib.myb_oem_segment_count(oem)):
        meta = {}
        for key in keys:
            value = lib.myb_oem_meta(oem, s, key.encode("ascii"))
            if value != ffi.NULL:
                meta[key] = ffi.string(value).decode("utf-8")
        states = []
        for k in range(lib.myb_oem_state_count(oem, s)):
            lib.myb_oem_state(oem, s, k, epoch, y, ffi.NULL)
            states.append((_label(epoch), tuple(y)))
        segments.append({"meta": meta, "states": states})
    return segments
