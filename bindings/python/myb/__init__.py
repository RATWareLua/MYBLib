# SPDX-License-Identifier: BSD-3-Clause
"""myb for Python: the motion of objects in space.

Two layers over the shared library of myb:

  * `myb.lib` and `myb.ffi` (cffi) reach every function and structure of the C interface,
    by the names of the C headers: `myb.lib.myb_kepler_propagate(...)`. The declarations
    are generated from the headers, so nothing is missing and nothing can differ.
  * The classes and functions of this package wrap the interface with ownership (objects
    are released when Python drops them), errors as exceptions (`myb.Error` carries the
    status of the library) and Python callables for terms, providers, events, measurements
    and the functions of the solvers. This module holds kernels and time, ephemerides, the
    environment and the force model, propagation, events, the messages of the CCSDS and the
    array forms of time and Kepler; the modules beside it the rest — `myb.nbody` (N bodies,
    the system of bodies), `myb.estim` (stations, measurements, the fit and the filter),
    `myb.cr3bp`, `myb.targeting` (transfers, Lambert, optimal control, thrust, initial orbits,
    relative motion, conjunctions, uncertainty), `myb.attitude`, `myb.observe`, `myb.files`
    (tables, SPK, SP3, CK and SCLK, OPM, OMM, OCM, TLE) and `myb.frames` (the providers of
    orientation and their array forms) — whose names are also names of this module.

No arithmetic is done here: a number obtained through this package has the bits of the same
calls made in C. States are sequences of six floats (position in m, velocity in m/s); times
of a model are seconds from the epoch of its environment. Series (the array forms, the
states of many bodies) come in as sequences — of numbers, of tuples, or an array.array('d')
taken without a copy — and go out as lists: of floats, of tuples, of instants.

The optional libraries of ext/ are loaded on request: `myb.extension("sgp4")`.
"""
from ._core import (Elements, Equinoctial, Error, Orientation, Provider, Rotation, _calendar,
                    _check, _Calls, _doubles, _flat, _instant, _instants, _label,
                    _orientation_provider, _Owner, _pv, _rows, _series, _settle, _state_of,
                    _state_provider, _string, _text, _times, _tuple)
from ._native import extension, ffi, lib, library_path

version = (lib.MYB_VERSION_MAJOR, lib.MYB_VERSION_MINOR, lib.MYB_VERSION_PATCH,
           lib.MYB_VERSION_TWEAK)
__version__ = ".".join(str(part) for part in version)


# ---- Time ---------------------------------------------------------------------------------------

def time_from_seconds(seconds_past_epoch):
    """An instant (myb_time) from seconds past J2000 of its time scale."""
    return lib.myb_time_from_seconds(float(seconds_past_epoch))


def seconds(instant):
    """The seconds past J2000 of an instant, as one float."""
    return lib.myb_time_seconds(instant)


def time_from_seconds_many(values):
    """Instants from seconds past J2000 (myb_time_from_seconds_many): a list of myb_time."""
    given, count = _series(values, 1, "seconds")
    out = ffi.new("myb_time[]", max(count, 1))
    lib.myb_time_from_seconds_many(given, count, out)
    return _instants(out, count)


def seconds_many(instants):
    """The seconds past J2000 of instants (myb_time_seconds_many): a list of floats."""
    given, count = _times(instants)
    out = ffi.new("double[]", max(count, 1))
    lib.myb_time_seconds_many(given, count, out)
    return list(_tuple(out, count))


def calendar_to_time(label):
    """The instant of a calendar label (year, month, day, hour, minute, second) of a
    continuous scale (myb_calendar_to_time)."""
    out = ffi.new("myb_time*")
    _check(lib.myb_calendar_to_time(_calendar(label), out))
    return out[0]


def time_to_calendar(instant):
    out = ffi.new("myb_calendar*")
    lib.myb_time_to_calendar(instant, out)
    return _label(out)


def _labels(labels):
    labels = list(labels)
    return ffi.new("myb_calendar[]", [_calendar(label)[0] for label in labels] or 1), len(labels)


def _failed_check(status, failed):
    _check(status, index=None if status == lib.MYB_OK else failed[0])


def calendar_to_time_many(labels):
    """The instants of calendar labels (myb_calendar_to_time_many): a list of myb_time. On a
    failure Error.index names the label."""
    given, count = _labels(labels)
    out = ffi.new("myb_time[]", max(count, 1))
    failed = ffi.new("size_t*")
    _failed_check(lib.myb_calendar_to_time_many(given, count, out, failed), failed)
    return _instants(out, count)


def time_to_calendar_many(instants):
    """The calendar labels of instants (myb_time_to_calendar_many): a list of labels."""
    given, count = _times(instants)
    out = ffi.new("myb_calendar[]", max(count, 1))
    lib.myb_time_to_calendar_many(given, count, out)
    return [_label(out[k]) for k in range(count)]


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

    def load_text(self, text):
        raw = text.encode("utf-8")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_kpool_load_text(self.pointer, raw, len(raw), diag), diag)

    def describe(self, name):
        """("numbers" or "strings", the number of values) of a variable."""
        kind = ffi.new("myb_kvar_type*")
        count = ffi.new("size_t*")
        _check(lib.myb_kpool_describe(self.pointer, name.encode("ascii"), kind, count))
        return ("numbers" if kind[0] == lib.MYB_KVAR_NUMBERS else "strings", count[0])

    def numbers(self, name):
        values = ffi.new("double**")
        count = ffi.new("size_t*")
        _check(lib.myb_kpool_numbers(self.pointer, name.encode("ascii"), values, count))
        return _tuple(values[0], count[0])

    def string(self, name, index=0):
        value = ffi.new("char**")
        length = ffi.new("size_t*")
        _check(lib.myb_kpool_string(self.pointer, name.encode("ascii"), int(index), value,
                                    length))
        return ffi.string(value[0], length[0]).decode("utf-8")

    def strings(self, name):
        kind, count = self.describe(name)
        return [self.string(name, k) for k in range(count)]


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

    def convert_many(self, instants, from_scale, to_scale):
        """convert for a series (myb_timesys_convert_many): a list of myb_time."""
        given, count = _times(instants)
        out = ffi.new("myb_time[]", max(count, 1))
        failed = ffi.new("size_t*")
        _failed_check(lib.myb_timesys_convert_many(self.pointer, given, count,
                                                   self._scale(from_scale),
                                                   self._scale(to_scale), out, failed), failed)
        return _instants(out, count)

    def from_utc_many(self, labels, to_scale):
        given, count = _labels(labels)
        out = ffi.new("myb_time[]", max(count, 1))
        failed = ffi.new("size_t*")
        _failed_check(lib.myb_timesys_from_utc_many(self.pointer, given, count,
                                                    self._scale(to_scale), out, failed), failed)
        return _instants(out, count)

    def to_utc_many(self, instants, from_scale):
        given, count = _times(instants)
        out = ffi.new("myb_calendar[]", max(count, 1))
        failed = ffi.new("size_t*")
        _failed_check(lib.myb_timesys_to_utc_many(self.pointer, given, count,
                                                  self._scale(from_scale), out, failed), failed)
        return [_label(out[k]) for k in range(count)]


# ---- Two bodies ---------------------------------------------------------------------------------

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


def kepler_stm(mu, state, dt):
    """The state a time dt later and its transition matrix (six rows)."""
    out = ffi.new("myb_pv*")
    stm = ffi.new("double[36]")
    _check(lib.myb_kepler_stm(float(mu), _pv(state), float(dt), out, stm))
    return _state_of(out), _rows(stm, 6, 6)


def _pvs(states):
    flat = _flat(states)
    if len(flat) % 6:
        raise ValueError("states: six numbers each, %d numbers given" % len(flat))
    count = len(flat) // 6
    pvs = ffi.new("myb_pv[]", max(count, 1))
    for k in range(count):
        pvs[k].r = flat[6 * k:6 * k + 3]
        pvs[k].v = flat[6 * k + 3:6 * k + 6]
    return pvs, count


def _elements_array(elements, ctype, fields):
    elements = list(elements)
    made = ffi.new(ctype + "[]", max(len(elements), 1))
    for k, e in enumerate(elements):
        values = _flat(e)
        if len(values) != 6:
            raise ValueError("a set of elements is six numbers, %d given" % len(values))
        for name, value in zip(fields, values):
            setattr(made[k], name, value)
    return made, len(elements)


def elements_from_state_many(mu, states):
    """elements_from_state for a series (myb_elements_from_pv_many): a list of Elements."""
    pvs, count = _pvs(states)
    out = ffi.new("myb_elements[]", max(count, 1))
    failed = ffi.new("size_t*")
    _failed_check(lib.myb_elements_from_pv_many(float(mu), pvs, count, out, failed), failed)
    return [Elements(e.p, e.e, e.i, e.raan, e.argp, e.nu) for e in (out[k] for k in range(count))]


def state_from_elements_many(mu, elements):
    """state_from_elements for a series (myb_pv_from_elements_many): a list of states."""
    given, count = _elements_array(elements, "myb_elements", Elements._fields)
    out = ffi.new("myb_pv[]", max(count, 1))
    failed = ffi.new("size_t*")
    _failed_check(lib.myb_pv_from_elements_many(float(mu), given, count, out, failed), failed)
    return [_state_of(out[k]) for k in range(count)]


def equinoctial_from_state(mu, state):
    """The modified equinoctial elements (p, f, g, h, k, L) of a state."""
    out = ffi.new("myb_equinoctial*")
    _check(lib.myb_equinoctial_from_pv(float(mu), _pv(state), out))
    return Equinoctial(out.p, out.f, out.g, out.h, out.k, out.l)


def state_from_equinoctial(mu, elements):
    given = ffi.new("myb_equinoctial*", dict(zip(Equinoctial._fields,
                                                 (float(v) for v in elements))))
    out = ffi.new("myb_pv*")
    _check(lib.myb_pv_from_equinoctial(float(mu), given, out))
    return _state_of(out)


def equinoctial_from_state_many(mu, states):
    pvs, count = _pvs(states)
    out = ffi.new("myb_equinoctial[]", max(count, 1))
    failed = ffi.new("size_t*")
    _failed_check(lib.myb_equinoctial_from_pv_many(float(mu), pvs, count, out, failed), failed)
    return [Equinoctial(e.p, e.f, e.g, e.h, e.k, e.l) for e in (out[k] for k in range(count))]


def state_from_equinoctial_many(mu, elements):
    given, count = _elements_array(elements, "myb_equinoctial", Equinoctial._fields)
    out = ffi.new("myb_pv[]", max(count, 1))
    failed = ffi.new("size_t*")
    _failed_check(lib.myb_pv_from_equinoctial_many(float(mu), given, count, out, failed), failed)
    return [_state_of(out[k]) for k in range(count)]


def kepler_propagate_many(mu, state, dts):
    """The states at the offsets dts from one state (myb_kepler_propagate_many): a list."""
    given, count = _series(dts, 1, "dt")
    out = ffi.new("myb_pv[]", max(count, 1))
    failed = ffi.new("size_t*")
    _failed_check(lib.myb_kepler_propagate_many(float(mu), _pv(state), given, count, out,
                                                failed), failed)
    return [_state_of(out[k]) for k in range(count)]


def kepler_stm_many(mu, state, dts):
    """The states at the offsets dts and their transition matrices (myb_kepler_stm_many):
    a list of states and a list of matrices (six rows each)."""
    given, count = _series(dts, 1, "dt")
    out = ffi.new("myb_pv[]", max(count, 1))
    stm = ffi.new("double[]", 36 * max(count, 1))
    failed = ffi.new("size_t*")
    _failed_check(lib.myb_kepler_stm_many(float(mu), _pv(state), given, count, out, stm, failed),
                  failed)
    flat = _rows(stm, 6 * count, 6)
    return ([_state_of(out[k]) for k in range(count)],
            [flat[6 * k:6 * k + 6] for k in range(count)])


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

    def load_memory(self, data):
        """Adds an SPK image (bytes); the set keeps it."""
        image = ffi.from_buffer(bytes(data))
        diag = ffi.new("myb_diag*")
        _check(lib.myb_spk_load_memory(self.pointer, image, len(image), diag), diag)
        self._kept.append(image)

    def state(self, target, center, instant):
        """Target relative to center at an instant (or seconds past J2000) of the files'
        time scale."""
        y = ffi.new("double[6]")
        _check(lib.myb_spk_state(self.pointer, int(target), int(center), _instant(instant), y,
                                 ffi.NULL))
        return tuple(y)

    def state_in(self, graph, target, center, instant, frame):
        """The state in `frame` of a graph of frames (myb_spk_state_in)."""
        y = ffi.new("double[6]")
        _settle(lib.myb_spk_state_in(self.pointer, graph.pointer, int(target), int(center),
                                     _instant(instant), int(frame), y), graph)
        return tuple(y)

    def state_many(self, target, center, instants, frames=False):
        """The states at many instants (myb_spk_state_many): a list of states, and with
        frames=True also the list of the frame codes."""
        given, count = _times(instants)
        out = ffi.new("double[]", 6 * max(count, 1))
        codes = ffi.new("int32_t[]", max(count, 1)) if frames else ffi.NULL
        failed = ffi.new("size_t*")
        _failed_check(lib.myb_spk_state_many(self.pointer, int(target), int(center), given,
                                             count, out, codes, failed), failed)
        states = _rows(out, count, 6)
        return (states, list(_tuple(codes, count))) if frames else states

    def state_in_many(self, graph, target, center, instants, frame):
        """state_in at many instants (myb_spk_state_in_many): a list of states."""
        given, count = _times(instants)
        out = ffi.new("double[]", 6 * max(count, 1))
        failed = ffi.new("size_t*")
        status = lib.myb_spk_state_in_many(self.pointer, graph.pointer, int(target), int(center),
                                           given, count, int(frame), out, failed)
        _settle(status, graph, index=None if status == lib.MYB_OK else failed[0])
        return _rows(out, count, 6)

    def segments(self):
        """The segments of the set: target, center, frame, type, start, stop, file, name."""
        info = ffi.new("myb_spk_segment_info*")
        found = []
        for index in range(lib.myb_spk_segment_count(self.pointer)):
            lib.myb_spk_segment(self.pointer, index, info)
            found.append({"target": info.target, "center": info.center, "frame": info.frame,
                          "type": info.type, "start": info.start, "stop": info.stop,
                          "file": info.file,
                          "name": ffi.string(info.name).decode("ascii", "replace")})
        return found

    def segment_size(self, index):
        return lib.myb_spk_segment_size(self.pointer, int(index))

    def segment_data(self, index, offset=0, count=None):
        """Double words of a segment (myb_spk_segment_data): a list of floats."""
        if count is None:
            count = self.segment_size(index) - int(offset)
        out = ffi.new("double[]", max(int(count), 1))
        _check(lib.myb_spk_segment_data(self.pointer, int(index), int(offset), int(count), out))
        return list(_tuple(out, int(count)))

    def body(self, target, center, frame=0, graph=None):
        """A body of the set as a provider of states (myb_spk_body_state): for a body of an
        Environment, a System, a chain of legs, the observables."""
        return _SpkBody(self, target, center, frame, graph)


class _SpkBody:
    def __init__(self, ephemeris, target, center, frame, graph=None):
        self.ephemeris = ephemeris
        self.graph = graph
        self.context = ffi.new("myb_spk_body*", {"spk": ephemeris.pointer, "target": int(target),
                                                 "center": int(center), "frame": int(frame)})
        if graph is not None:
            self.context.graph = graph.pointer
        self._kept = [ephemeris, graph]


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
        return cls._adopt(out[0])

    @property
    def degree(self):
        return lib.myb_gravity_field_degree(self.pointer)

    @property
    def gm(self):
        return lib.myb_gravity_field_gm(self.pointer)

    @property
    def radius(self):
        return lib.myb_gravity_field_radius(self.pointer)


class Polyhedron(_Owner):
    """A closed triangular mesh (myb/grav.h): the shape of a body for light, contact and
    its field. vertices: three numbers each (m); faces: three vertex indices each,
    counter-clockwise seen from outside; g_rho: G times the density (s^-2)."""

    _destroy = staticmethod(lib.myb_polyhedron_destroy)

    def __init__(self, vertices, faces, g_rho=0.0):
        points, count = _series(vertices, 3, "vertices")
        indices = [int(k) for k in _flat(faces)]
        if len(indices) % 3:
            raise ValueError("faces: three indices each")
        made = ffi.new("size_t[]", indices or 1)
        out = ffi.new("myb_polyhedron**")
        _check(lib.myb_polyhedron_create(points, count, made, len(indices) // 3, float(g_rho),
                                         ffi.NULL, out))
        super().__init__(out[0])

    def volume(self):
        """(volume m^3, centre of mass)."""
        center = ffi.new("double[3]")
        volume = lib.myb_polyhedron_volume(self.pointer, center)
        return volume, tuple(center)


# ---- The environment ----------------------------------------------------------------------------

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

    def add_body(self, id, gm, state=None, radius=0.0, luminosity=0.0, limb=(0.0, 0.0),
                 orientation=None, spin=(0.0, 0.0, 0.0)):
        """Adds a body and returns its index. state: None for the body at the origin; a
        body of an Ephemeris (Ephemeris.body); a Provider; or a callable f(seconds) -> six
        numbers, the state relative to the origin at an instant given in seconds past J2000
        of the model's scale. orientation (None: the inertial axes): a Provider of
        myb/frames.h or a callable f(seconds) -> Rotation or Orientation."""
        body = ffi.new("myb_env_body*")
        body.id = int(id)
        body.gm = float(gm)
        body.radius = float(radius)
        body.luminosity = float(luminosity)
        body.limb = [float(limb[0]), float(limb[1])]
        body.spin = list(_doubles(spin, 3, "spin"))
        if state is not None:
            provider = _state_provider(state, self._calls)
            body.state = provider.function
            body.state_context = provider.context
            self._kept.append(provider)
        if orientation is not None:
            provider = _orientation_provider(orientation, self._calls)
            body.orientation = provider.function
            body.orientation_context = provider.context
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
    """The forces on one object (myb/model.h) with the ready-made terms of myb/env.h and
    the terms of the CR3BP, of thrust arcs and of the Q-law. The model keeps alive the
    environments, fields and structures its terms point to."""

    _destroy = staticmethod(lib.myb_force_model_destroy)

    def __init__(self):
        out = ffi.new("myb_force_model**")
        _check(lib.myb_force_model_create(ffi.NULL, out))
        super().__init__(out[0])
        self._calls = _Calls()
        self._contexts = {}

    def _add(self, name, function, jacobian, context, *kept):
        index = ffi.new("size_t*")
        _check(lib.myb_force_model_add(self.pointer, name.encode("ascii"), function, jacobian,
                                       context, index))
        self._kept += [context] + list(kept)
        self._contexts[index[0]] = context
        return index[0]

    def add_raw(self, name, function, jacobian, context, *kept):
        """Any term of the C interface: a function and a Jacobian of `myb.lib` (or
        ffi.NULL) and their context (cdata). The context and whatever else is given are
        kept alive with the model."""
        return self._add(name, function, jacobian, context, *kept)

    def term_context(self, index):
        """The context (cdata) of a term added through this class — for a parameter
        (add_parameter with ffi.addressof) or to change a value in place."""
        return self._contexts[int(index)]

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

    def add_cr3bp(self, name, mu):
        """The circular restricted three-body problem in its rotating frame and units
        (myb_cr3bp_term, with its Jacobian)."""
        term = ffi.new("double*", float(mu))
        return self._add(name, lib.myb_cr3bp_term, lib.myb_cr3bp_term_jacobian, term)

    def add_thrust(self, name, thrust):
        """The burn arcs of a myb.targeting.Thrust (myb_thrust_acceleration)."""
        return self._add(name, lib.myb_thrust_acceleration, ffi.NULL, thrust.pointer, thrust)

    def add_qlaw(self, name, engine):
        """The engine of a Q-law, a myb.targeting.QLawTerm (myb_term_qlaw)."""
        return self._add(name, lib.myb_term_qlaw, ffi.NULL, engine.pointer, engine)

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
        """Raises what a callable of the model, or of what it keeps, raised."""
        _settle(status, self)

    @property
    def term_names(self):
        return [ffi.string(lib.myb_force_model_name(self.pointer, k)).decode("ascii")
                for k in range(lib.myb_force_model_count(self.pointer))]

    def set_enabled(self, index, enabled):
        lib.myb_force_model_set_enabled(self.pointer, int(index), 1 if enabled else 0)

    def enabled(self, index):
        return lib.myb_force_model_enabled(self.pointer, int(index)) != 0

    def add_parameter(self, name, value, step, *kept):
        """A parameter for the variational equations and the estimators: value is a double*
        (cdata: ffi.new("double*", x), or ffi.addressof(model.term_context(k), "member"))
        that a term reads; step that of its central differences. Returns its index."""
        index = ffi.new("size_t*")
        _check(lib.myb_force_model_add_parameter(self.pointer, name.encode("ascii"), value,
                                                 float(step), index))
        self._kept += [value] + list(kept)
        return index[0]

    def bind_parameter(self, parameter, term):
        _check(lib.myb_force_model_bind_parameter(self.pointer, int(parameter), int(term)))

    @property
    def parameter_names(self):
        return [ffi.string(lib.myb_force_model_parameter_name(self.pointer, k)).decode("ascii")
                for k in range(lib.myb_force_model_parameter_count(self.pointer))]

    def parameter(self, index):
        """The present value of a parameter (the estimators correct it in place)."""
        return lib.myb_force_model_parameter(self.pointer, int(index))[0]


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


def write_oem(originator, creation_date, segments, notation=None):
    """An orbit ephemeris message of the CCSDS as text. creation_date and the epochs are
    calendar labels (year, month, day, hour, minute, second); each segment is a mapping
    with object_name, object_id, center_name, ref_frame, time_system, epochs and states
    (six numbers each, m and m/s), and optionally interpolation with interpolation_degree,
    and covariances: a list of (epoch, 6 x 6 matrix) with covariance_frame. notation: "kvn"
    or "xml" (myb_oem_write_as), None for the default of myb_oem_write (KVN)."""
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
        covariances = list(segment.get("covariances", ()))
        if covariances:
            cov_epochs = ffi.new("myb_calendar[]", [_calendar(e)[0] for e, _ in covariances])
            numbers = [_flat(m) for _, m in covariances]
            if any(len(m) != 36 for m in numbers):
                raise ValueError("segment %d: a covariance is 36 numbers" % k)
            matrices = ffi.new("double[]", [v for m in numbers for v in m])
            kept += [cov_epochs, matrices]
            made[k].covariance_count = len(covariances)
            made[k].covariance_epochs = cov_epochs
            made[k].covariances = matrices
            if segment.get("covariance_frame") is not None:
                made[k].covariance_frame = text(segment["covariance_frame"])
    data = ffi.new("myb_oem_data*")
    data.originator = text(originator)
    data.creation_date = _calendar(creation_date)[0]
    data.segment_count = len(segments)
    data.segments = made
    if notation is None:
        return _text(lambda buffer, size, length: lib.myb_oem_write(data, buffer, size, length))
    return _text(lambda buffer, size, length: lib.myb_oem_write_as(
        data, _notations[notation], buffer, size, length))


def _read_oem(oem):
    keys = ("OBJECT_NAME", "OBJECT_ID", "CENTER_NAME", "REF_FRAME", "TIME_SYSTEM", "START_TIME",
            "STOP_TIME", "INTERPOLATION", "INTERPOLATION_DEGREE")
    segments = []
    epoch = ffi.new("myb_calendar*")
    y = ffi.new("double[6]")
    a = ffi.new("double[3]")
    frame = ffi.new("char**")
    matrix = ffi.new("double[36]")
    for s in range(lib.myb_oem_segment_count(oem)):
        meta = {}
        for key in keys:
            value = lib.myb_oem_meta(oem, s, key.encode("ascii"))
            if value != ffi.NULL:
                meta[key] = ffi.string(value).decode("utf-8")
        states = []
        accelerations = []
        for k in range(lib.myb_oem_state_count(oem, s)):
            given = lib.myb_oem_state(oem, s, k, epoch, y, a)
            states.append((_label(epoch), tuple(y)))
            accelerations.append(tuple(a) if given == 3 else None)
        covariances = []
        for k in range(lib.myb_oem_covariance_count(oem, s)):
            lib.myb_oem_covariance(oem, s, k, epoch, frame, matrix)
            covariances.append((_label(epoch), _string(frame[0]), _rows(matrix, 6, 6)))
        segment = {"meta": meta, "states": states, "covariances": covariances}
        if any(acceleration is not None for acceleration in accelerations):
            segment["accelerations"] = accelerations
        segments.append(segment)
    return segments


def read_oem(text):
    """An orbit ephemeris message in either notation: a list of segments, each a mapping
    with "meta" (the keywords of its metadata), "states" (pairs of a calendar label and
    six numbers), "covariances" (triples of a label, the frame or None and six rows) and,
    where the lines give them, "accelerations" (three numbers or None for each state)."""
    raw = text.encode("utf-8")
    out = ffi.new("myb_oem**")
    diag = ffi.new("myb_diag*")
    _check(lib.myb_oem_load(raw, len(raw), ffi.NULL, out, diag), diag)
    return _read_oem(ffi.gc(out[0], lib.myb_oem_destroy))


def read_oem_file(path):
    out = ffi.new("myb_oem**")
    diag = ffi.new("myb_diag*")
    _check(lib.myb_oem_load_file(str(path).encode("utf-8"), ffi.NULL, out, diag), diag)
    return _read_oem(ffi.gc(out[0], lib.myb_oem_destroy))


def oem_header(text, key):
    """The value of a header keyword of an ephemeris message (myb_oem_header), or None."""
    raw = text.encode("utf-8")
    out = ffi.new("myb_oem**")
    diag = ffi.new("myb_diag*")
    _check(lib.myb_oem_load(raw, len(raw), ffi.NULL, out, diag), diag)
    oem = ffi.gc(out[0], lib.myb_oem_destroy)
    return _string(lib.myb_oem_header(oem, key.encode("ascii")))


# ---- The rest of the interface -------------------------------------------------------------------

from . import attitude, cr3bp, estim, files, frames, nbody, observe, targeting  # noqa: E402
from .attitude import *  # noqa: E402,F401,F403
from .cr3bp import *  # noqa: E402,F401,F403
from .estim import *  # noqa: E402,F401,F403
from .files import *  # noqa: E402,F401,F403
from .frames import *  # noqa: E402,F401,F403
from .nbody import *  # noqa: E402,F401,F403
from .observe import *  # noqa: E402,F401,F403
from .targeting import *  # noqa: E402,F401,F403

__all__ = [
    "Elements", "Environment", "Ephemeris", "Equinoctial", "Error", "ForceModel",
    "GravityField", "KernelPool", "Orbit", "Orientation", "Polyhedron", "Provider", "Rotation",
    "TimeSystem", "calendar_to_time", "calendar_to_time_many", "elements_from_state",
    "elements_from_state_many", "equinoctial_from_state", "equinoctial_from_state_many",
    "extension", "ffi", "kepler_propagate", "kepler_propagate_many", "kepler_stm",
    "kepler_stm_many", "lib", "library_path", "oem_header", "read_oem", "read_oem_file",
    "seconds", "seconds_many", "state_from_elements", "state_from_elements_many",
    "state_from_equinoctial", "state_from_equinoctial_many", "time_from_seconds",
    "time_from_seconds_many", "time_to_calendar", "time_to_calendar_many", "version",
    "write_oem",
] + attitude.__all__ + cr3bp.__all__ + estim.__all__ + files.__all__ + frames.__all__ \
    + nbody.__all__ + observe.__all__ + targeting.__all__
