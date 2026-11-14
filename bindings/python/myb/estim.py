# SPDX-License-Identifier: BSD-3-Clause
"""Estimation of an orbit from observations (myb/estim.h): stations on bodies of an
environment, measurements with their partials, least squares over an arc and the sequential
filter with its smoother.

A measurement is a Measurement: the function of the C interface (or a Python callable) and
its context, with the number of values it gives. An Observation holds the time, the
measurement, the observed values and their sigmas. The unknowns are the state at the epoch,
the parameters of the force model (ForceModel.add_parameter) and the systematic shifts, in
that order; matrices over them come out as lists of rows."""
from collections import namedtuple

from ._core import (_attitude_provider, _Calls, _check, _doubles, _flat, _Owner, _rows,
                    _settle, _tuple)
from ._native import ffi, lib

__all__ = [
    "Filter", "FitResult", "Measurement", "Observation", "Station", "fit_orbit",
    "measure_azel", "measure_light_azel", "measure_light_radec", "measure_light_range",
    "measure_light_range_rate", "measure_offset", "measure_position", "measure_radec",
    "measure_range", "measure_range_rate", "measure_state",
]

FitResult = namedtuple("FitResult", "state iterations converged values rms_before rms bias "
                                    "covariance residuals")


class Station:
    """An observer on a body of an environment (myb_station): at geodetic longitude,
    latitude (rad) and height (m) on the ellipsoid (a, f) of the body (myb_station_init), or
    at a body-fixed position with the axes east, north, up (rows) given (Station.at).
    motion: a myb_station_motion* of the C interface (cdata), what moves the station on its
    body (myb/site.h), with what it points to kept alive by `keep`."""

    def __init__(self, env, body, longitude, latitude, height, a, f=0.0, motion=None, keep=()):
        self.pointer = ffi.new("myb_station*")
        ellipsoid = ffi.new("myb_ellipsoid*", {"a": float(a), "f": float(f)})
        _check(lib.myb_station_init(self.pointer, env.pointer, int(body), ellipsoid,
                                    float(longitude), float(latitude), float(height)))
        self._set_motion(motion)
        self._kept = [env, motion] + list(keep)

    @classmethod
    def at(cls, env, body, position, enu=(1, 0, 0, 0, 1, 0, 0, 0, 1), motion=None, keep=()):
        made = cls.__new__(cls)
        made.pointer = ffi.new("myb_station*")
        made.pointer.env = env.pointer
        made.pointer.body = int(body)
        made.pointer.position = list(_doubles(position, 3, "position"))
        made.pointer.enu = list(_doubles(enu, 9, "enu"))
        made._set_motion(motion)
        made._kept = [env, motion] + list(keep)
        return made

    def _set_motion(self, motion):
        if motion is not None:
            self.pointer.motion = motion

    @property
    def position(self):
        return tuple(self.pointer.position)

    @property
    def enu(self):
        return tuple(self.pointer.enu)

    def state(self, t):
        """The state of the station relative to the origin of the environment at the model
        time t (myb_station_state)."""
        out = ffi.new("double[6]")
        _settle(lib.myb_station_state(self.pointer, float(t), out), self)
        return tuple(out)


class Measurement:
    """count values computed from the state of the object: a function of the C interface
    (myb_measure_function) and its context, or a Python callable (Measurement.custom)."""

    def __init__(self, function, context, count, *kept):
        self.function = function
        self.context = context
        self.count = int(count)
        self._kept = [context] + list(kept)

    @classmethod
    def custom(cls, function, count):
        """function(t, y) -> (values, partials): `count` values and their derivatives by the
        state (count rows of six; None when only the values are asked for — the second
        argument of the call, `want_partials`, says which): function(t, y, want_partials)."""
        calls = _Calls()

        @ffi.callback("myb_status(void*, double, const double[6], double*, double*)")
        def measure(context, t, y, values, partials):
            try:
                want = partials != ffi.NULL
                v, p = function(t, tuple(y[0:6]), want)
                values[0:count] = list(_doubles(v, count, "values"))
                if want:
                    partials[0:6 * count] = list(_doubles(p, 6 * count, "partials"))
                return lib.MYB_OK
            except BaseException as error:  # noqa: B902 - nothing may cross the C code
                return calls.failed(error)

        made = cls(measure, ffi.NULL, count, function)
        made._calls = calls
        return made

    def __call__(self, t, y, partials=True):
        """The values at the model time t and the state y, and their partials (count rows of
        six) or None."""
        values = ffi.new("double[]", self.count)
        rows = ffi.new("double[]", 6 * self.count) if partials else ffi.NULL
        _settle(self.function(self.context, float(t), _doubles(y, 6, "y"), values, rows), self)
        return _tuple(values, self.count), (_rows(rows, self.count, 6) if partials else None)


def measure_position():
    """The position of the object itself (3 values)."""
    return Measurement(lib.myb_measure_position, ffi.NULL, 3)


def measure_state():
    """The whole state of the object (6 values)."""
    return Measurement(lib.myb_measure_state, ffi.NULL, 6)


def measure_range(station):
    return Measurement(lib.myb_measure_range, station.pointer, 1, station)


def measure_range_rate(station):
    return Measurement(lib.myb_measure_range_rate, station.pointer, 1, station)


def measure_radec(station):
    """Right ascension and declination from a station (2 values; the first periodic)."""
    return Measurement(lib.myb_measure_radec, station.pointer, 2, station)


def measure_azel(station):
    """Azimuth and elevation at a station (2 values; the first periodic)."""
    return Measurement(lib.myb_measure_azel, station.pointer, 2, station)


def _light(station, c, gm, two_way, aberration):
    return ffi.new("myb_light_measure*", {"station": station.pointer, "c": float(c),
                                          "gm": float(gm), "two_way": 1 if two_way else 0,
                                          "aberration": 1 if aberration else 0})


def measure_light_range(station, c, gm=0.0, two_way=False):
    """The range with the light time (myb_measure_light_range)."""
    return Measurement(lib.myb_measure_light_range, _light(station, c, gm, two_way, False), 1,
                       station)


def measure_light_range_rate(station, c, gm=0.0, two_way=False):
    return Measurement(lib.myb_measure_light_range_rate, _light(station, c, gm, two_way, False),
                       1, station)


def measure_light_radec(station, c, gm=0.0, aberration=False):
    return Measurement(lib.myb_measure_light_radec, _light(station, c, gm, False, aberration), 2,
                       station)


def measure_light_azel(station, c, gm=0.0, aberration=False):
    return Measurement(lib.myb_measure_light_azel, _light(station, c, gm, False, aberration), 2,
                       station)


def measure_offset(measurement, attitude, offset, gm=0.0):
    """A measurement taken at a point offset (m, body frame) from the centre of mass, the
    body turned by an attitude law: a Provider of myb/surface.h (myb_attitude_orbital, ...)
    or a callable f(t, x, v) -> nine numbers (myb_measure_offset)."""
    calls = _Calls()
    law = _attitude_provider(attitude, calls)
    context = ffi.new("myb_offset_measure*")
    context.measure = measurement.function
    context.context = measurement.context
    context.attitude = law.function
    context.attitude_context = law.context
    context.offset = list(_doubles(offset, 3, "offset"))
    context.gm = float(gm)
    made = Measurement(lib.myb_measure_offset, context, measurement.count, measurement, law)
    made._calls = calls
    return made


class Observation:
    """An observation at the model time t: the measurement, the observed values, their
    sigmas, which values are angles (bits of `periodic`) and the first systematic shift
    (`bias`: 0 none, k the shifts k - 1 ...)."""

    def __init__(self, t, measurement, value, sigma, periodic=0, bias=0):
        count = measurement.count
        self.measurement = measurement
        self.raw = ffi.new("myb_observation*")
        self.raw.t = float(t)
        self.raw.measure = measurement.function
        self.raw.context = measurement.context
        self.raw.count = count
        self.raw.value[0:count] = list(_doubles(value, count, "value"))
        self.raw.sigma[0:count] = list(_doubles(sigma, count, "sigma"))
        self.raw.periodic = int(periodic)
        self.raw.bias = int(bias)
        self._kept = [measurement]

    @property
    def t(self):
        return self.raw.t


def _unknowns(model, bias_count):
    return 6 + lib.myb_force_model_parameter_count(model.pointer) + int(bias_count)


def fit_orbit(model, t0, y0, observations, max_iterations=0, tolerance=0.0, epsilon=0.0,
              fixed=None, prior_sigma=None, bias=None, covariance=True, residuals=True):
    """Least squares over an arc (myb_fit_orbit): the state at t0 from a first guess y0, the
    parameters of the model (corrected in place) and the systematic shifts (`bias`: their
    first guesses, one per shift). A FitResult: the state, the counts, the RMS before and
    after, the shifts, the covariance (rows) and the residuals. On a failure Error.partial is
    the state of the last iteration."""
    observations = list(observations)
    raw = ffi.new("myb_observation[]", max(len(observations), 1))
    for k, observation in enumerate(observations):
        raw[k] = observation.raw[0]
    bias = None if bias is None else _flat(bias)
    bias_count = 0 if bias is None else len(bias)
    n = _unknowns(model, bias_count)
    options = ffi.new("myb_fit_options*")
    options.max_iterations = int(max_iterations)
    options.tolerance = float(tolerance)
    options.epsilon = float(epsilon)
    kept = []
    if fixed is not None:
        flags = ffi.new("unsigned char[]", [1 if f else 0 for f in fixed])
        if len(flags) != n:
            raise ValueError("fixed: %d flags expected" % n)
        options.fixed = flags
        kept.append(flags)
    if prior_sigma is not None:
        sigma = _doubles(prior_sigma, n, "prior_sigma")
        options.prior_sigma = sigma
        kept.append(sigma)
    shifts = ffi.NULL
    if bias_count:
        shifts = ffi.new("double[]", bias)
        options.bias_count = bias_count
        options.bias = shifts
    state = _doubles(y0, 6, "y0")
    result = ffi.new("myb_fit_result*")
    matrix = ffi.new("double[]", n * n) if covariance else ffi.NULL
    values = sum(o.measurement.count for o in observations)
    left = ffi.new("double[]", max(values, 1)) if residuals else ffi.NULL
    status = lib.myb_fit_orbit(model.pointer, float(t0), state, raw, len(observations), options,
                               ffi.NULL, result, matrix, left)
    _settle(status, model, *observations, partial=lambda: tuple(state))
    return FitResult(tuple(state), result.iterations, result.converged != 0, result.values,
                     result.rms_before, result.rms, _tuple(shifts, bias_count),
                     _rows(matrix, n, n) if covariance else None,
                     _tuple(left, result.values) if residuals else None)


class Filter(_Owner):
    """The extended Kalman filter on a force model (myb_filter): epsilon of IAS15, the
    spectral density q of an unmodelled acceleration, bias_count systematic shifts, history
    instants kept for the smoother."""

    _destroy = staticmethod(lib.myb_filter_destroy)

    def __init__(self, model, epsilon=0.0, process_noise=0.0, bias_count=0, history=0):
        options = ffi.new("myb_filter_options*", {"epsilon": float(epsilon),
                                                  "process_noise": float(process_noise),
                                                  "bias_count": int(bias_count),
                                                  "history": int(history)})
        out = ffi.new("myb_filter**")
        _check(lib.myb_filter_create(model.pointer, options, ffi.NULL, out))
        super().__init__(out[0])
        self._model = model
        self._kept.append(model)
        self._bias_count = int(bias_count)
        self._n = _unknowns(model, bias_count)

    @property
    def unknowns(self):
        """6 + parameters + shifts."""
        return self._n

    def start(self, t0, y0, covariance, t_end):
        """Starts at t0 with y0 and the covariance (rows) of all the unknowns — 6 + np +
        bias_count a side: myb_filter_start reads the rows and columns of the shifts too."""
        _settle(lib.myb_filter_start(self.pointer, float(t0), _doubles(y0, 6, "y0"),
                                     _doubles(covariance, self._n * self._n, "covariance"),
                                     float(t_end)), self._model)

    def update(self, observation):
        """Propagates to the observation and applies it: the innovation (observed -
        predicted)."""
        innovation = ffi.new("double[]", observation.measurement.count)
        _settle(lib.myb_filter_update(self.pointer, observation.raw, innovation), self._model,
                observation)
        self._kept.append(observation)
        return _tuple(innovation, observation.measurement.count)

    def predict(self, t):
        _settle(lib.myb_filter_predict(self.pointer, float(t)), self._model)

    def get(self):
        """(t, state, covariance rows)."""
        t = ffi.new("double*")
        y = ffi.new("double[6]")
        p = ffi.new("double[]", self._n * self._n)
        lib.myb_filter_get(self.pointer, t, y, p)
        return t[0], tuple(y), _rows(p, self._n, self._n)

    def set_bias(self, bias):
        _check(lib.myb_filter_set_bias(self.pointer, _doubles(bias, self._bias_count, "bias")))

    @property
    def bias(self):
        out = ffi.new("double[]", max(self._bias_count, 1))
        lib.myb_filter_bias(self.pointer, out)
        return _tuple(out, self._bias_count)

    @property
    def steps(self):
        """The instants kept for the smoother."""
        return lib.myb_filter_steps(self.pointer)

    def smooth(self):
        _check(lib.myb_filter_smooth(self.pointer))

    def step(self, index):
        """(t, unknowns, covariance rows) at an instant kept."""
        t = ffi.new("double*")
        z = ffi.new("double[]", self._n)
        p = ffi.new("double[]", self._n * self._n)
        lib.myb_filter_step(self.pointer, int(index), t, z, p)
        return t[0], _tuple(z, self._n), _rows(p, self._n, self._n)

