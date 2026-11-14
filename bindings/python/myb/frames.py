# SPDX-License-Identifier: BSD-3-Clause
"""Frames (myb/frames.h) as far as the array forms need them: the models of the rotation of
the Earth and of the bodies, binary PCK files, the graph of frames, orientation providers —
each with its call at one instant and its array form at many (one crossing into C for the
whole series).

A rotation is a Rotation (m, dm): the nine numbers of the matrix (rows) and of its
derivative per second. An orientation is an Orientation (rotation, ut, xp, yp). Instants are
myb_time; the array forms also take seconds past J2000 (turned into instants by
myb_time_from_seconds_many)."""
from collections import namedtuple

from ._core import (Provider, _check, _instant, _Owner, _orientation_of,
                    _orientation_provider, _rotation_in, _rotation_of, _rows, _series, _settle,
                    _times, _Calls, _instants, _tuple)
from ._native import ffi, lib

__all__ = [
    "Cip", "Eop", "EopValues", "Equinox", "FrameGraph", "Pck", "PckWriter", "Wgccre",
    "orientation_fixed",
    "orientation_many", "rotation_apply_position_many", "rotation_apply_state_many",
]

EopValues = namedtuple("EopValues", "xp yp dut1 lod dx dy")

_equinox_frames = {"mod": lib.MYB_EQUINOX_MOD, "tod": lib.MYB_EQUINOX_TOD,
                   "teme": lib.MYB_EQUINOX_TEME}


def _rotations(values):
    values = list(values)
    made = ffi.new("myb_rotation[]", max(len(values), 1))
    for k, rotation in enumerate(values):
        _rotation_in(rotation, ffi.addressof(made, k))
    return made, len(values)


def _rotation_list(buffer, count):
    return [_rotation_of(buffer[k]) for k in range(count)]


def _failed(status, failed, *sources):
    _settle(status, *sources, index=None if status == lib.MYB_OK else failed[0])


def rotation_apply_position_many(rotations, positions):
    """Each position carried by its rotation (myb_rotation_apply_position_many): a list."""
    made, count = _rotations(rotations)
    given, n = _series(positions, 3, "positions")
    if n != count:
        raise ValueError("one rotation a position: %d and %d" % (count, n))
    out = ffi.new("double[]", 3 * max(count, 1))
    lib.myb_rotation_apply_position_many(made, given, count, out)
    return _rows(out, count, 3)


def rotation_apply_state_many(rotations, states):
    """Each state carried by its rotation, the turning of the frame included
    (myb_rotation_apply_state_many): a list."""
    made, count = _rotations(rotations)
    given, n = _series(states, 6, "states")
    if n != count:
        raise ValueError("one rotation a state: %d and %d" % (count, n))
    out = ffi.new("double[]", 6 * max(count, 1))
    lib.myb_rotation_apply_state_many(made, given, count, out)
    return _rows(out, count, 6)


def orientation_fixed(rotation):
    """A constant orientation as a provider (myb_orientation_fixed)."""
    return Provider(lib.myb_orientation_fixed, _rotation_in(rotation))


def orientation_many(provider, instants):
    """What a provider of orientations gives at many instants (myb_orientation_many): a
    list of Orientation. provider: a Provider (orientation_fixed, Cip.orientation, ...,
    Ck.orientation) or a callable f(seconds) -> Rotation or Orientation."""
    calls = _Calls()
    provider = _orientation_provider(provider, calls)
    given, count = _times(instants)
    out = ffi.new("myb_orientation[]", max(count, 1))
    failed = ffi.new("size_t*")
    status = lib.myb_orientation_many(provider.function, provider.context, given, count, out,
                                      failed)
    _failed(status, failed, calls, provider)
    return [_orientation_of(out[k]) for k in range(count)]


class Cip(_Owner):
    """The CIO-based rotation of the Earth (myb_cip): the model <prefix>_* of a pool
    (data/frames/iau2006-cip.tk: "IAU2006A")."""

    _destroy = staticmethod(lib.myb_cip_destroy)

    def __init__(self, pool, prefix="IAU2006A"):
        out = ffi.new("myb_cip**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_cip_create(pool.pointer, prefix.encode("ascii"), ffi.NULL, out, diag), diag)
        super().__init__(out[0])

    def xys(self, tt):
        """X, Y, s and their rates at the TT instant tt."""
        xys = ffi.new("double[3]")
        rates = ffi.new("double[3]")
        lib.myb_cip_xys(self.pointer, _instant(tt), xys, rates)
        return tuple(xys), tuple(rates)

    def xys_many(self, tt):
        """xys at many instants (myb_cip_xys_many): a list of X, Y, s and a list of rates."""
        given, count = _times(tt)
        xys = ffi.new("double[]", 3 * max(count, 1))
        rates = ffi.new("double[]", 3 * max(count, 1))
        lib.myb_cip_xys_many(self.pointer, given, count, xys, rates)
        return _rows(xys, count, 3), _rows(rates, count, 3)

    def tabulate(self, tt0, tt1, step=0.0):
        _check(lib.myb_cip_tabulate(self.pointer, _instant(tt0), _instant(tt1), float(step)))

    def era(self, ut1):
        """(the rotation angle, its rate) at the UT1 instant ut1."""
        rate = ffi.new("double*")
        angle = lib.myb_cip_era(self.pointer, _instant(ut1), rate)
        return angle, rate[0]

    def sprime(self, tt):
        return lib.myb_cip_sprime(self.pointer, _instant(tt))

    def terrestrial(self, tt, ut1, eop=None):
        """The rotation celestial -> terrestrial (myb_cip_terrestrial); eop: EopValues or
        None (all zero)."""
        out = ffi.new("myb_rotation*")
        lib.myb_cip_terrestrial(self.pointer, _instant(tt), _instant(ut1), _eop_values(eop),
                                out)
        return _rotation_of(out)

    def terrestrial_many(self, tt, ut1, eop=None):
        """terrestrial at many instants (myb_cip_terrestrial_many): eop a list of EopValues
        (one an instant) or None. A list of Rotation."""
        given_tt, count = _times(tt)
        given_ut1, n = _times(ut1)
        if n != count:
            raise ValueError("as many UT1 instants as TT instants")
        values = ffi.NULL
        if eop is not None:
            eop = list(eop)
            if len(eop) != count:
                raise ValueError("one set of EOP an instant")
            values = ffi.new("myb_eop_values[]", [_eop_values(e)[0] for e in eop] or 1)
        out = ffi.new("myb_rotation[]", max(count, 1))
        lib.myb_cip_terrestrial_many(self.pointer, given_tt, given_ut1, values, count, out)
        return _rotation_list(out, count)

    def orientation(self, eop, timesys, tt="TT"):
        """The rotation of the Earth with observed EOP as a provider (myb_orientation_cip)."""
        context = ffi.new("myb_cip_orientation*", {"cip": self.pointer, "eop": eop.pointer,
                                                   "timesys": timesys.pointer,
                                                   "tt": timesys._scale(tt)})
        return Provider(lib.myb_orientation_cip, context, self, eop, timesys)


def _eop_values(eop):
    if eop is None:
        return ffi.NULL
    return ffi.new("myb_eop_values*", dict(zip(EopValues._fields, (float(v) for v in eop))))


class Eop(_Owner):
    """A table of Earth orientation parameters (myb_eop): files of the C04 or finals2000A
    form, optionally the diurnal and semidiurnal tides."""

    _destroy = staticmethod(lib.myb_eop_destroy)

    def __init__(self, c04=None, finals=None):
        out = ffi.new("myb_eop**")
        _check(lib.myb_eop_create(ffi.NULL, out))
        super().__init__(out[0])
        if c04 is not None:
            self.load_c04_file(c04)
        if finals is not None:
            self.load_finals_file(finals)

    def _load(self, function, text):
        raw = text.encode("utf-8")
        diag = ffi.new("myb_diag*")
        _check(function(self.pointer, raw, len(raw), diag), diag)

    def _load_file(self, function, path):
        diag = ffi.new("myb_diag*")
        _check(function(self.pointer, str(path).encode("utf-8"), diag), diag)

    def load_c04(self, text):
        self._load(lib.myb_eop_load_c04, text)

    def load_c04_file(self, path):
        self._load_file(lib.myb_eop_load_c04_file, path)

    def load_finals(self, text):
        self._load(lib.myb_eop_load_finals, text)

    def load_finals_file(self, path):
        self._load_file(lib.myb_eop_load_finals_file, path)

    def set_tides(self, pool, prefix="IERS10_EOPT"):
        diag = ffi.new("myb_diag*")
        _check(lib.myb_eop_set_tides(self.pointer, pool.pointer, prefix.encode("ascii"), diag),
               diag)

    def rows(self):
        """(rows, first MJD, last MJD)."""
        first = ffi.new("double*")
        last = ffi.new("double*")
        count = lib.myb_eop_rows(self.pointer, first, last)
        return count, first[0], last[0]

    def at(self, timesys, instant, scale="TT"):
        """(EopValues, the instant on UT1) at an instant of a scale of the time system."""
        values = ffi.new("myb_eop_values*")
        ut1 = ffi.new("myb_time*")
        _check(lib.myb_eop_at(self.pointer, timesys.pointer, _instant(instant),
                              timesys._scale(scale), values, ut1))
        return EopValues(values.xp, values.yp, values.dut1, values.lod, values.dx,
                         values.dy), ut1[0]

    def at_many(self, timesys, instants, scale="TT"):
        """at for a series (myb_eop_at_many): a list of EopValues and a list of UT1
        instants."""
        given, count = _times(instants)
        values = ffi.new("myb_eop_values[]", max(count, 1))
        ut1 = ffi.new("myb_time[]", max(count, 1))
        failed = ffi.new("size_t*")
        _failed(lib.myb_eop_at_many(self.pointer, timesys.pointer, given, count,
                                    timesys._scale(scale), values, ut1, failed), failed)
        return ([EopValues(v.xp, v.yp, v.dut1, v.lod, v.dx, v.dy)
                 for v in (values[k] for k in range(count))], _instants(ut1, count))


class Equinox(_Owner):
    """The equinox-based systems of the Earth (myb_equinox): mean of date, true of date,
    TEME; the model <prefix>_* of a pool."""

    _destroy = staticmethod(lib.myb_equinox_destroy)

    def __init__(self, pool, prefix="IAU1980"):
        out = ffi.new("myb_equinox**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_equinox_create(pool.pointer, prefix.encode("ascii"), ffi.NULL, out, diag),
               diag)
        super().__init__(out[0])

    def rotation(self, frame, tt):
        """The rotation J2000 -> frame ("mod", "tod", "teme") at tt."""
        out = ffi.new("myb_rotation*")
        lib.myb_equinox_rotation(self.pointer, _equinox_frames[frame], _instant(tt), out)
        return _rotation_of(out)

    def rotation_many(self, frame, tt):
        given, count = _times(tt)
        out = ffi.new("myb_rotation[]", max(count, 1))
        lib.myb_equinox_rotation_many(self.pointer, _equinox_frames[frame], given, count, out)
        return _rotation_list(out, count)

    def gmst(self, ut1):
        rate = ffi.new("double*")
        angle = lib.myb_equinox_gmst(self.pointer, _instant(ut1), rate)
        return angle, rate[0]

    def teme_pef(self, ut1):
        out = ffi.new("myb_rotation*")
        lib.myb_equinox_teme_pef(self.pointer, _instant(ut1), out)
        return _rotation_of(out)

    def teme_pef_many(self, ut1):
        given, count = _times(ut1)
        out = ffi.new("myb_rotation[]", max(count, 1))
        lib.myb_equinox_teme_pef_many(self.pointer, given, count, out)
        return _rotation_list(out, count)


class Wgccre(_Owner):
    """The rotation model of a body by the IAU WGCCRE (myb_wgccre) from a pool
    (pck00011.tpc)."""

    _destroy = staticmethod(lib.myb_wgccre_destroy)

    def __init__(self, pool, body):
        out = ffi.new("myb_wgccre**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_wgccre_create(pool.pointer, int(body), ffi.NULL, out, diag), diag)
        super().__init__(out[0])

    @property
    def frame(self):
        return lib.myb_wgccre_frame(self.pointer)

    def euler(self, t):
        angles = ffi.new("double[3]")
        rates = ffi.new("double[3]")
        lib.myb_wgccre_euler(self.pointer, _instant(t), angles, rates)
        return tuple(angles), tuple(rates)

    def rotation(self, t):
        out = ffi.new("myb_rotation*")
        lib.myb_wgccre_rotation(self.pointer, _instant(t), out)
        return _rotation_of(out)

    def rotation_many(self, t):
        given, count = _times(t)
        out = ffi.new("myb_rotation[]", max(count, 1))
        lib.myb_wgccre_rotation_many(self.pointer, given, count, out)
        return _rotation_list(out, count)

    def orientation(self, base=None):
        """The model as a provider (myb_orientation_wgccre); base: a Rotation from the
        inertial axes in use to the model's reference frame, None when they coincide."""
        made = ffi.NULL if base is None else _rotation_in(base)
        context = ffi.new("myb_wgccre_orientation*", {"model": self.pointer, "base": made})
        return Provider(lib.myb_orientation_wgccre, context, self, made)


class Pck(_Owner):
    """A set of binary PCK files (myb_pck)."""

    _destroy = staticmethod(lib.myb_pck_destroy)

    def __init__(self, *paths):
        out = ffi.new("myb_pck**")
        _check(lib.myb_pck_create(ffi.NULL, out))
        super().__init__(out[0])
        for path in paths:
            self.load_file(path)

    def load_file(self, path):
        diag = ffi.new("myb_diag*")
        _check(lib.myb_pck_load_file(self.pointer, str(path).encode("utf-8"), diag), diag)

    def load_memory(self, data):
        image = ffi.from_buffer(bytes(data))
        diag = ffi.new("myb_diag*")
        _check(lib.myb_pck_load_memory(self.pointer, image, len(image), diag), diag)
        self._kept.append(image)

    def euler(self, id, t):
        """(angles, rates, base frame) of frame class id at t."""
        angles = ffi.new("double[3]")
        rates = ffi.new("double[3]")
        base = ffi.new("int32_t*")
        _check(lib.myb_pck_euler(self.pointer, int(id), _instant(t), angles, rates, base))
        return tuple(angles), tuple(rates), base[0]

    def rotation(self, id, t):
        """(Rotation base -> id, base frame) at t."""
        out = ffi.new("myb_rotation*")
        base = ffi.new("int32_t*")
        _check(lib.myb_pck_rotation(self.pointer, int(id), _instant(t), out, base))
        return _rotation_of(out), base[0]

    def rotation_many(self, id, t):
        """rotation at many instants (myb_pck_rotation_many): a list of Rotation and a list
        of base frames."""
        given, count = _times(t)
        out = ffi.new("myb_rotation[]", max(count, 1))
        bases = ffi.new("int32_t[]", max(count, 1))
        failed = ffi.new("size_t*")
        _failed(lib.myb_pck_rotation_many(self.pointer, int(id), given, count, out, bases,
                                          failed), failed)
        return _rotation_list(out, count), list(_tuple(bases, count))

    def orientation(self, id, base=None):
        """A frame of the set as a provider (myb_orientation_pck)."""
        made = ffi.NULL if base is None else _rotation_in(base)
        context = ffi.new("myb_pck_orientation*", {"pck": self.pointer, "id": int(id),
                                                   "base": made})
        return Provider(lib.myb_orientation_pck, context, self, made)


class PckWriter(_Owner):
    """Segments collected in memory and laid out as one binary PCK file (myb_pck_writer):
    Chebyshev polynomials of the Euler angles fitted to a provider of orientations, or any
    type from its words."""

    _destroy = staticmethod(lib.myb_pck_writer_destroy)

    def __init__(self, internal_name=None):
        out = ffi.new("myb_pck_writer**")
        name = ffi.NULL if internal_name is None else ffi.new("char[]",
                                                              internal_name.encode("ascii"))
        _check(lib.myb_pck_writer_create(name, ffi.NULL, out))
        super().__init__(out[0])

    def add_chebyshev(self, frame, base, name, orientation, start, stop, intervals, degree,
                      max_error=False):
        """A segment of type 2 for frame class `frame` in the inertial frame `base`, fitted
        to a provider (a Provider or a callable f(seconds) -> Rotation or Orientation); with
        max_error=True the largest angle between the provider and the polynomials."""
        calls = _Calls()
        provider = _orientation_provider(orientation, calls)
        error = ffi.new("double*") if max_error else ffi.NULL
        _settle(lib.myb_pck_writer_add_chebyshev(self.pointer, int(frame), int(base),
                                                 name.encode("ascii"), provider.function,
                                                 provider.context, float(start), float(stop),
                                                 int(intervals), int(degree), error),
                calls, provider)
        return error[0] if max_error else None

    def add_raw(self, frame, base, name, type, start, stop, words):
        values, count = _series(words, 1, "words")
        _check(lib.myb_pck_writer_add_raw(self.pointer, int(frame), int(base),
                                          name.encode("ascii"), int(type), float(start),
                                          float(stop), values, count))

    @property
    def size(self):
        return lib.myb_pck_writer_size(self.pointer)

    def write(self):
        """The file as bytes."""
        size = self.size
        buffer = ffi.new("char[]", max(size, 1))
        _check(lib.myb_pck_writer_write(self.pointer, buffer, size))
        return ffi.buffer(buffer, size)[:]

    def write_file(self, path):
        _check(lib.myb_pck_writer_write_file(self.pointer, str(path).encode("utf-8")))


class FrameGraph(_Owner):
    """A graph of frames (myb_frame_graph): frames given in others by providers of
    orientation or by constant rotations; the rotation between any two of one tree."""

    _destroy = staticmethod(lib.myb_frame_graph_destroy)

    def __init__(self, timesys=None, scale=None):
        options = ffi.new("myb_frame_graph_options*")
        scale_text = None
        if timesys is not None:
            options.timesys = timesys.pointer
        if scale is not None:
            scale_text = ffi.new("char[]", scale.encode("ascii"))
            options.scale = scale_text
        out = ffi.new("myb_frame_graph**")
        _check(lib.myb_frame_graph_create(options, ffi.NULL, out))
        super().__init__(out[0])
        self._kept += [timesys, scale_text]
        self._calls = _Calls()

    def add(self, id, base, orientation, scale=None):
        """Frame id given in base by a provider (a Provider or a callable f(seconds) ->
        Rotation or Orientation, the rotation base -> id); scale: that of its argument."""
        provider = _orientation_provider(orientation, self._calls)
        scale_text = ffi.NULL if scale is None else ffi.new("char[]", scale.encode("ascii"))
        _check(lib.myb_frame_graph_add(self.pointer, int(id), int(base), provider.function,
                                       provider.context, scale_text))
        self._kept += [provider, scale_text]

    def add_fixed(self, id, base, rotation):
        _check(lib.myb_frame_graph_add_fixed(self.pointer, int(id), int(base),
                                             _rotation_in(rotation)))

    def add_tkframe(self, pool, id):
        diag = ffi.new("myb_diag*")
        _check(lib.myb_frame_graph_add_tkframe(self.pointer, pool.pointer, int(id), diag), diag)

    def rotation(self, from_frame, to_frame, t):
        out = ffi.new("myb_rotation*")
        _settle(lib.myb_frame_graph_rotation(self.pointer, int(from_frame), int(to_frame),
                                              _instant(t), out), self)
        return _rotation_of(out)

    def rotation_many(self, from_frame, to_frame, t):
        """rotation at many instants (myb_frame_graph_rotation_many): a list of Rotation."""
        given, count = _times(t)
        out = ffi.new("myb_rotation[]", max(count, 1))
        failed = ffi.new("size_t*")
        _failed(lib.myb_frame_graph_rotation_many(self.pointer, int(from_frame), int(to_frame),
                                                  given, count, out, failed), failed, self)
        return _rotation_list(out, count)

    def orientation(self, from_frame, to_frame):
        """Two frames of the graph as a provider (myb_orientation_graph)."""
        context = ffi.new("myb_graph_orientation*", {"graph": self.pointer,
                                                     "from": int(from_frame),
                                                     "to": int(to_frame)})
        return Provider(lib.myb_orientation_graph, context, self)
