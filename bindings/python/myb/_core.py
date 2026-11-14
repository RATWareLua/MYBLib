# SPDX-License-Identifier: BSD-3-Clause
"""What the modules of the package share: the error, ownership, the crossing of numbers and
of Python callables into C. Nothing here computes: values pass through as they are."""
import array
from collections import namedtuple

from ._native import ffi, lib

Elements = namedtuple("Elements", "p e i raan argp nu")
Equinoctial = namedtuple("Equinoctial", "p f g h k l")
Rotation = namedtuple("Rotation", "m dm")
Orientation = namedtuple("Orientation", "rotation ut xp yp")


class Error(Exception):
    """A call of the library failed: `status` is its code (lib.MYB_E_DOMAIN, ...), `name`
    the name of the code, `message` what a parser said (with the line), if anything.

    `index`: for the array forms, the element that failed (the outputs before it were
    filled). `partial`: where the library leaves what it reached in its outputs on a failure
    (MYB_E_NOCONV of the solvers, the state of a fit on its last iteration), that — else
    None."""

    def __init__(self, status, message="", index=None, partial=None):
        self.status = int(status)
        self.name = ffi.string(lib.myb_status_name(status)).decode("ascii")
        self.message = message
        self.index = index
        self.partial = partial
        text = "%s: %s" % (self.name, message) if message else self.name
        if index is not None:
            text += " (element %d)" % index
        super().__init__(text)


def _check(status, diag=None, index=None, partial=None):
    if status == lib.MYB_OK:
        return
    message = ""
    if diag is not None:
        message = ffi.string(diag.message).decode("utf-8", "replace")
        if diag.line > 0:
            message = "line %d: %s" % (diag.line, message)
    raise Error(status, message, index=index, partial=partial() if partial else None)


# ---- Numbers in and out --------------------------------------------------------------------------

def _flat(values):
    """Numbers, or sequences of them (nested at any depth), as one list of floats."""
    if isinstance(values, array.array) and values.typecode == "d":
        return values.tolist()
    out = []
    for value in values:
        try:
            out.append(float(value))
        except TypeError:
            out.extend(_flat(value))
    return out


def _doubles(values, count, what):
    """double[count] from `count` numbers (or nested sequences of as many)."""
    values = _flat(values)
    if len(values) != count:
        raise ValueError("%s: %d numbers expected, %d given" % (what, count, len(values)))
    return ffi.new("double[]", values)


def _optional(values, count, what):
    return ffi.NULL if values is None else _doubles(values, count, what)


def _series(values, width, what):
    """A series of elements of `width` numbers each — a sequence of sequences, or flat — as
    double[] and the number of elements."""
    if isinstance(values, array.array) and values.typecode == "d":
        flat = values
    else:
        flat = _flat(values)
    if len(flat) % width:
        raise ValueError("%s: a multiple of %d numbers expected, %d given"
                         % (what, width, len(flat)))
    count = len(flat) // width
    if isinstance(flat, array.array):
        return ffi.from_buffer("double[]", flat), count
    return ffi.new("double[]", flat if flat else 1), count


def _tuple(buffer, count):
    return tuple(ffi.unpack(buffer, count)) if count else ()


def _rows(buffer, rows, columns):
    """A matrix or a series out of a flat double[]: `rows` tuples of `columns` numbers."""
    if rows == 0 or columns == 0:
        return [() for _ in range(rows)]
    flat = ffi.unpack(buffer, rows * columns)
    return [tuple(flat[k:k + columns]) for k in range(0, rows * columns, columns)]


def _square(values, what):
    """A square matrix given as rows or flat: double[n * n] and n."""
    flat = _flat(values)
    n = int(round(len(flat) ** 0.5))
    if n * n != len(flat) or n == 0:
        raise ValueError("%s: a square matrix expected, %d numbers given" % (what, len(flat)))
    return ffi.new("double[]", flat), n


def _text_in(value):
    """bytes for a char* argument (kept alive by the caller)."""
    return ffi.new("char[]", str(value).encode("utf-8"))


def _string(pointer):
    return None if pointer == ffi.NULL else ffi.string(pointer).decode("utf-8", "replace")


def _owned(ctype, value):
    """A struct copied out of an array: an element of a cffi array does not keep it alive."""
    return ffi.new(ctype + "*", value)[0]


def _calendar(label):
    """(year, month, day, hour, minute, second) -> myb_calendar*."""
    year, month, day, hour, minute, second = label
    return ffi.new("myb_calendar*", {"year": int(year), "month": int(month), "day": int(day),
                                     "hour": int(hour), "minute": int(minute),
                                     "second": float(second)})


def _label(calendar):
    return (calendar.year, calendar.month, calendar.day, calendar.hour, calendar.minute,
            calendar.second)


def _pv(state):
    y = _flat(state)
    if len(y) != 6:
        raise ValueError("a state is six numbers, %d given" % len(y))
    return ffi.new("myb_pv*", {"r": y[:3], "v": y[3:]})


def _state_of(pv):
    return (pv.r[0], pv.r[1], pv.r[2], pv.v[0], pv.v[1], pv.v[2])


def _rotation_of(rotation):
    return Rotation(tuple(rotation.m), tuple(rotation.dm))


def _rotation_in(rotation, out=None):
    """A Rotation (m, dm: nine numbers each, rows) or a bare matrix m into a myb_rotation."""
    out = ffi.new("myb_rotation*") if out is None else out
    if isinstance(rotation, Rotation) or (len(rotation) == 2 and not _is_number(rotation[0])):
        m, dm = rotation
    else:
        m, dm = rotation, None
    out.m = list(_doubles(m, 9, "a rotation"))
    out.dm = [0.0] * 9 if dm is None else list(_doubles(dm, 9, "its derivative"))
    return out


def _is_number(value):
    try:
        float(value)
        return True
    except TypeError:
        return False


def _text(write):
    """A text of the library's writers: they tell the length when the buffer is too small."""
    length = ffi.new("size_t*")
    status = write(ffi.NULL, 0, length)
    if status not in (lib.MYB_OK, lib.MYB_E_OVERFLOW):
        _check(status)
    buffer = ffi.new("char[]", length[0] + 1)
    _check(write(buffer, length[0] + 1, length))
    return ffi.string(buffer, length[0]).decode("utf-8")


def _times(instants):
    """myb_time[] from instants (myb_time) or from seconds past J2000 — these turned into
    instants by myb_time_from_seconds_many, as myb.time_from_seconds would one by one."""
    instants = list(instants)
    if instants and _is_number(instants[0]):
        seconds = ffi.new("double[]", [float(s) for s in instants])
        out = ffi.new("myb_time[]", len(instants))
        lib.myb_time_from_seconds_many(seconds, len(instants), out)
        return out, len(instants)
    return ffi.new("myb_time[]", instants if instants else 1), len(instants)


def _instant(instant):
    """An instant (myb_time), or seconds past J2000 made one by myb_time_from_seconds."""
    return lib.myb_time_from_seconds(float(instant)) if _is_number(instant) else instant


def _instants(buffer, count):
    return [_owned("myb_time", buffer[k]) for k in range(count)]


# ---- Objects of the library ----------------------------------------------------------------------

class _Owner:
    """An object of the library: released when dropped or closed, usable in `with`."""

    _destroy = None

    def __init__(self, pointer):
        self._pointer = ffi.gc(pointer, self._destroy)
        self._kept = []  # what the object points to must live as long as it does

    @classmethod
    def _adopt(cls, pointer):
        """An object of the class around a pointer that a function of the library made."""
        made = cls.__new__(cls)
        _Owner.__init__(made, pointer)
        return made

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


def _settle(status, *sources, diag=None, index=None, partial=None):
    """Raises what a Python callable behind any of the sources raised (an object with
    `_calls`, and what it keeps, at any depth), else the status of the call."""
    found = []
    seen = set()
    stack = list(sources)
    while stack:
        source = stack.pop()
        if source is None or id(source) in seen:
            continue
        seen.add(id(source))
        calls = source if isinstance(source, _Calls) else getattr(source, "_calls", None)
        if isinstance(calls, _Calls) and calls.pending is not None:
            found.append(calls.pending)
            calls.pending = None
        kept = getattr(source, "_kept", None)
        if isinstance(kept, list):
            stack.extend(kept)
    if found:
        raise found[0]
    _check(status, diag, index=index, partial=partial)


class Provider:
    """A function of the C interface with its context — a provider of states
    (myb_state_function), of orientations (myb_orientation_function), an attitude law — and
    whatever they need kept alive: Provider(lib.myb_orientation_fixed, rotation)."""

    def __init__(self, function, context, *kept):
        self.function = function
        self.context = context
        self._kept = [context] + list(kept)


def _state_provider(state, calls):
    """A provider of states (myb_state_function) from a Provider, a body of an Ephemeris
    (Ephemeris.body) or a Python callable f(seconds) -> six numbers, the argument the instant
    in seconds past J2000 of its scale."""
    if isinstance(state, Provider):
        return state
    if hasattr(state, "context") and hasattr(state, "ephemeris"):
        return Provider(lib.myb_spk_body_state, state.context, state)

    @ffi.callback("myb_status(void*, myb_time, double[6])")
    def provider(context, instant, out):
        try:
            values = _flat(state(lib.myb_time_seconds(instant)))
            if len(values) != 6:
                raise ValueError("a state is six numbers, %d given" % len(values))
            out[0:6] = values
            return lib.MYB_OK
        except BaseException as error:  # noqa: B902 - nothing may cross the C code
            return calls.failed(error)

    return Provider(provider, ffi.NULL, state)


def _orientation_provider(orientation, calls):
    """A provider of orientations (myb_orientation_function) from a Provider or a Python
    callable f(seconds) -> a Rotation (m, dm) or an Orientation (rotation, ut, xp, yp), ut in
    seconds past J2000 of the rotation time."""
    if isinstance(orientation, Provider):
        return orientation

    @ffi.callback("myb_status(void*, myb_time, myb_orientation*)")
    def provider(context, instant, out):
        try:
            value = orientation(lib.myb_time_seconds(instant))
            if isinstance(value, Orientation):
                _rotation_in(value.rotation, ffi.addressof(out, "rotation"))
                out.ut = lib.myb_time_from_seconds(float(value.ut))
                out.xp = float(value.xp)
                out.yp = float(value.yp)
            else:
                _rotation_in(value, ffi.addressof(out, "rotation"))
                out.ut = instant
                out.xp = 0.0
                out.yp = 0.0
            return lib.MYB_OK
        except BaseException as error:  # noqa: B902 - nothing may cross the C code
            return calls.failed(error)

    return Provider(provider, ffi.NULL, orientation)


def _orientation_of(orientation):
    return Orientation(_rotation_of(orientation.rotation), _owned("myb_time", orientation.ut),
                       orientation.xp, orientation.yp)


def _attitude_provider(attitude, calls):
    """An attitude law (myb_attitude_function of myb/surface.h) from a Provider or a Python
    callable f(t, x, v) -> the nine numbers of the matrix inertial -> body (rows)."""
    if isinstance(attitude, Provider):
        return attitude

    @ffi.callback("myb_status(void*, double, const double[3], const double[3], double[9])")
    def law(context, t, x, v, m):
        try:
            values = _flat(attitude(t, (x[0], x[1], x[2]), (v[0], v[1], v[2])))
            if len(values) != 9:
                raise ValueError("a matrix is nine numbers, %d given" % len(values))
            m[0:9] = values
            return lib.MYB_OK
        except BaseException as error:  # noqa: B902 - nothing may cross the C code
            return calls.failed(error)

    return Provider(law, ffi.NULL, attitude)


def _nlp_options(feasibility=0.0, optimality=0.0, max_outer=0, max_inner=0, weight=0.0):
    return ffi.new("myb_nlp_options*", {"feasibility": float(feasibility),
                                        "optimality": float(optimality),
                                        "max_outer": int(max_outer), "max_inner": int(max_inner),
                                        "weight": float(weight)})


NlpResult = namedtuple("NlpResult", "objective violation stationarity outer inner evaluations")


def _nlp_result(result):
    return NlpResult(result.objective, result.violation, result.stationarity, result.outer,
                     result.inner, result.evaluations)
