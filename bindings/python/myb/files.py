# SPDX-License-Identifier: BSD-3-Clause
"""Files: numbers and tables of numbers (myb/io.h), the SPK writer and SP3 precise orbits
(myb/kin.h), the orbit parameter, mean-elements and comprehensive messages of the CCSDS
(myb/ccsds.h; the ephemeris message is myb.write_oem and myb.read_oem), CK orientation files
and SCLK clocks with the CK writer (myb/ck.h), two-line elements (ext/sgp4).

Epochs are calendar labels (year, month, day, hour, minute, second). An OPM, an OMM and a TLE
are mappings with the names of the members of the C structures (myb_opm, myb_omm, myb_tle);
a member that is left out is zero."""
from collections import namedtuple

from ._core import (Provider, _calendar, _Calls, _check, _doubles, _flat, _label, _Owner,
                    _rotation_of, _rows, _series, _settle, _state_provider, _string, _text,
                    _tuple)
from ._native import extension, ffi, lib

__all__ = [
    "Ck", "CkPointing", "CkWriter", "Ocm", "Sclk", "Sp3", "SpkWriter", "Table",
    "format_double", "format_fixed", "format_scientific", "omm_from_tle", "parse_double",
    "read_omm", "read_opm", "tle_epoch", "tle_fit", "tle_format", "tle_from_omm",
    "tle_from_state", "tle_parse", "write_ocm", "write_ocm_line", "write_omm", "write_opm",
    "write_sp3", "write_table",
]

CkPointing = namedtuple("CkPointing", "q rate has_rate ticks frame rotation")

_notations = {"kvn": lib.MYB_ODM_KVN, "xml": lib.MYB_ODM_XML}
_ocm_kinds = {"traj": lib.MYB_OCM_TRAJ, "phys": lib.MYB_OCM_PHYS, "cov": lib.MYB_OCM_COV,
              "man": lib.MYB_OCM_MAN, "pert": lib.MYB_OCM_PERT, "od": lib.MYB_OCM_OD,
              "user": lib.MYB_OCM_USER}


def _kind(kind):
    return _ocm_kinds[kind.lower()] if isinstance(kind, str) else int(kind)


def _notation(notation):
    return _notations[notation] if isinstance(notation, str) else int(notation)


# ---- Numbers -----------------------------------------------------------------------------------

def parse_double(text):
    """(value, characters used) of the number at the start of text (myb_parse_double)."""
    raw = text.encode("ascii")
    value = ffi.new("double*")
    consumed = ffi.new("size_t*")
    _check(lib.myb_parse_double(raw, len(raw), value, consumed))
    return value[0], consumed[0]


def format_double(value):
    """The shortest decimal that reads back to the value (myb_format_double)."""
    out = ffi.new("char[]", lib.MYB_DOUBLE_CHARS)
    length = lib.myb_format_double(float(value), out)
    return ffi.string(out, length).decode("ascii")


def format_fixed(value, decimals):
    return _text(lambda buffer, size, length: lib.myb_format_fixed(
        float(value), int(decimals), buffer, size, length))


def format_scientific(value, digits):
    return _text(lambda buffer, size, length: lib.myb_format_scientific(
        float(value), int(digits), buffer, size, length))


class Table(_Owner):
    """A table of numbers read from CSV-like text (myb_table): Table.load(text),
    Table.load_file(path)."""

    _destroy = staticmethod(lib.myb_table_destroy)

    @classmethod
    def load(cls, text):
        raw = text.encode("utf-8")
        out = ffi.new("myb_table**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_table_load(raw, len(raw), ffi.NULL, out, diag), diag)
        return cls._adopt(out[0])

    @classmethod
    def load_file(cls, path):
        out = ffi.new("myb_table**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_table_load_file(str(path).encode("utf-8"), ffi.NULL, out, diag), diag)
        return cls._adopt(out[0])

    @property
    def rows(self):
        return lib.myb_table_rows(self.pointer)

    @property
    def columns(self):
        return lib.myb_table_columns(self.pointer)

    @property
    def names(self):
        return [ffi.string(lib.myb_table_name(self.pointer, k)).decode("utf-8")
                for k in range(self.columns)]

    def row(self, index):
        return _tuple(lib.myb_table_row(self.pointer, int(index)), self.columns)

    @property
    def values(self):
        """Every row: a list of tuples."""
        return [self.row(k) for k in range(self.rows)]


def write_table(values, names=None):
    """Rows of numbers as comma-separated text (myb_table_write); names: a first row."""
    rows = [_flat(row) for row in values]
    columns = len(names) if names is not None else (len(rows[0]) if rows else 0)
    if any(len(row) != columns for row in rows):
        raise ValueError("every row of %d numbers" % columns)
    flat = ffi.new("double[]", [v for row in rows for v in row] or 1)
    titles = ffi.NULL
    kept = []
    if names is not None:
        kept = [ffi.new("char[]", str(name).encode("utf-8")) for name in names]
        titles = ffi.new("char*[]", kept)
    return _text(lambda buffer, size, length: lib.myb_table_write(
        titles, columns, flat, len(rows), buffer, size, length))


# ---- SPK and SP3 -------------------------------------------------------------------------------

def _segment_header(target, center, frame, name):
    text = ffi.new("char[]", name.encode("ascii"))
    header = ffi.new("myb_spk_segment_header*", {"target": int(target), "center": int(center),
                                                 "frame": int(frame), "name": text})
    return header, text


class SpkWriter(_Owner):
    """Segments collected in memory and laid out as one SPK file (myb_spk_writer): types 9
    and 13 from states, 2 and 3 fitted to a provider, any type from its words."""

    _destroy = staticmethod(lib.myb_spk_writer_destroy)

    def __init__(self, internal_name=None):
        out = ffi.new("myb_spk_writer**")
        name = ffi.NULL if internal_name is None else ffi.new("char[]",
                                                              internal_name.encode("ascii"))
        _check(lib.myb_spk_writer_create(name, ffi.NULL, out))
        super().__init__(out[0])

    def add_states(self, target, center, frame, name, type, degree, epochs, states):
        """A segment of type 9 (Lagrange) or 13 (Hermite) through states at epochs (seconds
        past J2000 of TDB)."""
        header, text = _segment_header(target, center, frame, name)
        times, count = _series(epochs, 1, "epochs")
        values, n = _series(states, 6, "states")
        if n != count:
            raise ValueError("one state an epoch")
        _check(lib.myb_spk_writer_add_states(self.pointer, header, int(type), int(degree), times,
                                             values, count))

    def add_chebyshev(self, target, center, frame, name, type, state, start, stop, intervals,
                      degree, max_error=False):
        """A segment of type 2 or 3 fitted to a provider (a body of an Ephemeris, a Provider
        or a callable f(seconds) -> six numbers); with max_error=True the largest distance
        between the provider and the polynomials is returned."""
        calls = _Calls()
        provider = _state_provider(state, calls)
        header, text = _segment_header(target, center, frame, name)
        error = ffi.new("double*") if max_error else ffi.NULL
        _settle(lib.myb_spk_writer_add_chebyshev(self.pointer, header, int(type),
                                                 provider.function, provider.context,
                                                 float(start), float(stop), int(intervals),
                                                 int(degree), error), calls, provider)
        return error[0] if max_error else None

    def add_raw(self, target, center, frame, name, type, start, stop, words):
        """A segment of any type from its double words (Ephemeris.segment_data)."""
        header, text = _segment_header(target, center, frame, name)
        values, count = _series(words, 1, "words")
        _check(lib.myb_spk_writer_add_raw(self.pointer, header, int(type), float(start),
                                          float(stop), values, count))

    @property
    def size(self):
        return lib.myb_spk_writer_size(self.pointer)

    def write(self):
        """The file as bytes."""
        size = self.size
        buffer = ffi.new("char[]", max(size, 1))
        _check(lib.myb_spk_writer_write(self.pointer, buffer, size))
        return ffi.buffer(buffer, size)[:]

    def write_file(self, path):
        _check(lib.myb_spk_writer_write_file(self.pointer, str(path).encode("utf-8")))


class Sp3(_Owner):
    """An SP3-c or SP3-d file of precise orbits (myb_sp3): Sp3.load(text),
    Sp3.load_file(path)."""

    _destroy = staticmethod(lib.myb_sp3_destroy)

    @classmethod
    def load(cls, text):
        raw = text.encode("utf-8")
        out = ffi.new("myb_sp3**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_sp3_load(raw, len(raw), ffi.NULL, out, diag), diag)
        return cls._adopt(out[0])

    @classmethod
    def load_file(cls, path):
        out = ffi.new("myb_sp3**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_sp3_load_file(str(path).encode("utf-8"), ffi.NULL, out, diag), diag)
        return cls._adopt(out[0])

    @property
    def satellites(self):
        return [ffi.string(lib.myb_sp3_satellite(self.pointer, k)).decode("ascii")
                for k in range(lib.myb_sp3_satellite_count(self.pointer))]

    @property
    def epochs(self):
        label = ffi.new("myb_calendar*")
        found = []
        for k in range(lib.myb_sp3_epoch_count(self.pointer)):
            lib.myb_sp3_epoch(self.pointer, k, label)
            found.append(_label(label))
        return found

    @property
    def time_system(self):
        return _string(lib.myb_sp3_time_system(self.pointer))

    @property
    def frame(self):
        return _string(lib.myb_sp3_frame(self.pointer))

    def state(self, satellite, epoch):
        """(position or None, velocity or None) of a satellite at an epoch."""
        position = ffi.new("double[3]")
        velocity = ffi.new("double[3]")
        given = lib.myb_sp3_state(self.pointer, int(satellite), int(epoch), position, velocity)
        return (tuple(position) if given & 1 else None, tuple(velocity) if given & 2 else None)

    def clock(self, satellite, epoch):
        """(clock or None, rate or None)."""
        clock = ffi.new("double*")
        rate = ffi.new("double*")
        given = lib.myb_sp3_clock(self.pointer, int(satellite), int(epoch), clock, rate)
        return (clock[0] if given & 1 else None, rate[0] if given & 2 else None)

    def sigma(self, satellite, epoch):
        out = ffi.new("double[8]")
        lib.myb_sp3_sigma(self.pointer, int(satellite), int(epoch), out)
        return tuple(out)

    def accuracy(self, satellite):
        return lib.myb_sp3_accuracy(self.pointer, int(satellite))

    def flags(self, satellite, epoch):
        return lib.myb_sp3_flags(self.pointer, int(satellite), int(epoch))


def write_sp3(satellites, epochs, interval, positions, velocities=None, version="d",
              file_type="G", time_system="GPS", frame="IGS20", orbit_type="FIT", agency="MYB",
              data_used="ORBIT", clocks=None, clock_rates=None, accuracy=None, sigmas=None,
              flags=None):
    """An SP3 file as text (myb_sp3_write). positions[epoch][satellite]: three numbers (m;
    None or non-finite numbers for a missing one); velocities likewise or None; clocks,
    clock_rates [epoch][satellite]; accuracy [satellite]; sigmas [epoch][satellite]: eight
    numbers; flags [epoch][satellite]."""
    satellites = list(satellites)
    epochs = list(epochs)
    count = len(satellites) * len(epochs)
    nan = float("nan")
    kept = []

    def text(value):
        kept.append(ffi.new("char[]", str(value).encode("ascii")))
        return kept[-1]

    def per_record(values, width, what):
        flat = []
        for row in values:
            for value in row:
                flat.extend([nan] * width if value is None else _flat([value]))
        if len(flat) != width * count:
            raise ValueError("%s: %d numbers a record" % (what, width))
        made = ffi.new("double[]", flat or 1)
        kept.append(made)
        return made

    data = ffi.new("myb_sp3_data*")
    data.version = version.encode("ascii")
    data.file_type = file_type.encode("ascii")
    data.time_system = text(time_system)
    data.frame = text(frame)
    data.orbit_type = text(orbit_type)
    data.agency = text(agency)
    data.data_used = text(data_used)
    names = [text(s) for s in satellites]
    data.satellite_count = len(satellites)
    # the array itself is kept: what a field of the structure gives back is a bare pointer,
    # which does not keep the array alive (reading data.satellites here let it be collected)
    satellite_names = ffi.new("char*[]", names or [ffi.NULL])
    kept.append(satellite_names)
    data.satellites = satellite_names
    labels = ffi.new("myb_calendar[]", [_calendar(e)[0] for e in epochs] or 1)
    data.epoch_count = len(epochs)
    data.epochs = labels
    data.interval = float(interval)
    data.positions = per_record(positions, 3, "positions")
    if velocities is not None:
        data.velocities = per_record(velocities, 3, "velocities")
    if clocks is not None:
        data.clocks = per_record(clocks, 1, "clocks")
    if clock_rates is not None:
        data.clock_rates = per_record(clock_rates, 1, "clock_rates")
    if accuracy is not None:
        made = _doubles(accuracy, len(satellites), "accuracy")
        kept.append(made)
        data.accuracy = made
    if sigmas is not None:
        data.sigmas = per_record(sigmas, 8, "sigmas")
    if flags is not None:
        bits = [int(f) for row in flags for f in row]
        if len(bits) != count:
            raise ValueError("flags: one a record")
        made = ffi.new("unsigned char[]", bits or 1)
        kept.append(made)
        data.flags = made
    return _text(lambda buffer, size, length: lib.myb_sp3_write(data, buffer, size, length))


# ---- Messages of the CCSDS ---------------------------------------------------------------------

def _from_struct(value, ctype):
    """A structure of the C interface as a mapping of its members."""
    out = {}
    for name, field in ffi.typeof(ctype).fields:
        member = getattr(value, name)
        kind = field.type
        if kind.kind == "array" and kind.item.cname == "char":
            out[name] = ffi.string(member).decode("utf-8", "replace")
        elif kind.cname == "char":
            out[name] = member.decode("latin-1")
        elif kind.cname == "myb_calendar":
            out[name] = _label(member)
        elif kind.kind == "array" and kind.item.kind == "struct":
            out[name] = [_from_struct(member[k], kind.item.cname) for k in range(kind.length)]
        elif kind.kind == "array":
            out[name] = tuple(member)
        else:
            out[name] = member
    return out


def _to_struct(values, ctype, out=None):
    """A mapping of members into a structure of the C interface (others stay zero)."""
    out = ffi.new(ctype + "*") if out is None else out
    fields = dict(ffi.typeof(ctype).fields)
    for name, value in values.items():
        if name not in fields:
            raise KeyError("%s has no member %r" % (ctype, name))
        kind = fields[name].type
        if kind.kind == "array" and kind.item.cname == "char":
            raw = str(value).encode("utf-8")
            if len(raw) >= kind.length:
                raise ValueError("%s: at most %d characters" % (name, kind.length - 1))
            setattr(out, name, raw)
        elif kind.cname == "char":
            setattr(out, name, str(value).encode("latin-1")[:1])
        elif kind.cname == "myb_calendar":
            setattr(out, name, _calendar(value)[0])
        elif kind.kind == "array" and kind.item.kind == "struct":
            items = list(value)
            for k, item in enumerate(items):
                _to_struct(item, kind.item.cname, ffi.addressof(getattr(out, name), k))
        elif kind.kind == "array":
            setattr(out, name, _flat(value))
        else:
            setattr(out, name, value)
    return out


def _opm_of(opm):
    made = _from_struct(opm, "myb_opm")
    made["maneuvers"] = made["maneuvers"][:opm.maneuver_count]
    made["covariance"] = _rows(opm.covariance, 6, 6)
    return made


def read_opm(text):
    """An orbit parameter message in either notation (myb_opm_parse): a mapping."""
    raw = text.encode("utf-8")
    out = ffi.new("myb_opm*")
    diag = ffi.new("myb_diag*")
    _check(lib.myb_opm_parse(raw, len(raw), out, diag), diag)
    return _opm_of(out)


def write_opm(opm, notation=None):
    """An orbit parameter message from a mapping: notation "kvn" or "xml"
    (myb_opm_write_as), None for the default of myb_opm_write (KVN)."""
    values = dict(opm)
    if "maneuvers" in values and "maneuver_count" not in values:
        values["maneuver_count"] = len(values["maneuvers"])
    made = _to_struct(values, "myb_opm")
    if notation is None:
        return _text(lambda buffer, size, length: lib.myb_opm_write(made, buffer, size, length))
    return _text(lambda buffer, size, length: lib.myb_opm_write_as(
        made, _notation(notation), buffer, size, length))


def _omm_of(omm):
    made = _from_struct(omm, "myb_omm")
    made["covariance"] = _rows(omm.covariance, 6, 6)
    return made


def read_omm(text):
    """A mean-elements message in either notation (myb_omm_parse): a mapping."""
    raw = text.encode("utf-8")
    out = ffi.new("myb_omm*")
    diag = ffi.new("myb_diag*")
    _check(lib.myb_omm_parse(raw, len(raw), out, diag), diag)
    return _omm_of(out)


def write_omm(omm, notation=None):
    """A mean-elements message from a mapping, as write_opm."""
    made = _to_struct(dict(omm), "myb_omm")
    if notation is None:
        return _text(lambda buffer, size, length: lib.myb_omm_write(made, buffer, size, length))
    return _text(lambda buffer, size, length: lib.myb_omm_write_as(
        made, _notation(notation), buffer, size, length))


class Ocm(_Owner):
    """A comprehensive message (myb_ocm): Ocm.load(text), Ocm.load_file(path). Kinds of
    blocks: "traj", "phys", "cov", "man", "pert", "od", "user"."""

    _destroy = staticmethod(lib.myb_ocm_destroy)

    @classmethod
    def load(cls, text):
        raw = text.encode("utf-8")
        out = ffi.new("myb_ocm**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_ocm_load(raw, len(raw), ffi.NULL, out, diag), diag)
        return cls._adopt(out[0])

    @classmethod
    def load_file(cls, path):
        out = ffi.new("myb_ocm**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_ocm_load_file(str(path).encode("utf-8"), ffi.NULL, out, diag), diag)
        return cls._adopt(out[0])

    def header(self, key):
        return _string(lib.myb_ocm_header(self.pointer, key.encode("ascii")))

    def meta(self, key):
        return _string(lib.myb_ocm_meta(self.pointer, key.encode("ascii")))

    def block_count(self, kind):
        return lib.myb_ocm_block_count(self.pointer, _kind(kind))

    def value(self, kind, index, key):
        return _string(lib.myb_ocm_value(self.pointer, _kind(kind), int(index),
                                         key.encode("ascii")))

    def line_count(self, kind, index):
        return lib.myb_ocm_line_count(self.pointer, _kind(kind), int(index))

    def line_text(self, kind, index, line):
        return _string(lib.myb_ocm_line_text(self.pointer, _kind(kind), int(index), int(line)))

    def line(self, kind, index, line):
        """(time, values): the time a calendar label or seconds from EPOCH_TZERO, the values
        in the units of the message."""
        absolute = ffi.new("int*")
        epoch = ffi.new("myb_calendar*")
        offset = ffi.new("double*")
        count = ffi.new("size_t*")
        _check(lib.myb_ocm_line(self.pointer, _kind(kind), int(index), int(line), absolute,
                                epoch, offset, ffi.NULL, 0, count))
        values = ffi.new("double[]", max(count[0], 1))
        _check(lib.myb_ocm_line(self.pointer, _kind(kind), int(index), int(line), absolute,
                                epoch, offset, values, count[0], count))
        return (_label(epoch) if absolute[0] else offset[0]), _tuple(values, count[0])

    def state(self, index, line):
        """(time, state) of a Cartesian trajectory line in m, m/s, m/s^2 (3, 6 or 9
        numbers)."""
        absolute = ffi.new("int*")
        epoch = ffi.new("myb_calendar*")
        offset = ffi.new("double*")
        state = ffi.new("double[9]")
        count = ffi.new("int*")
        _check(lib.myb_ocm_state(self.pointer, int(index), int(line), absolute, epoch, offset,
                                 state, count))
        return (_label(epoch) if absolute[0] else offset[0]), _tuple(state, count[0])


def _pairs(items, kept):
    items = list(items.items()) if hasattr(items, "items") else list(items)
    made = ffi.new("myb_ocm_pair[]", max(len(items), 1))
    for k, (key, value) in enumerate(items):
        key_text = ffi.new("char[]", str(key).encode("utf-8"))
        value_text = ffi.new("char[]", str(value).encode("utf-8"))
        kept += [key_text, value_text]
        made[k].key = key_text
        made[k].value = value_text
    kept.append(made)
    return made, len(items)


def write_ocm(originator, creation_date, meta, blocks, notation="kvn"):
    """A comprehensive message (myb_ocm_write): meta a mapping (or pairs) with TIME_SYSTEM
    and EPOCH_TZERO; blocks a list of mappings with "kind", "pairs" and "lines" (text;
    write_ocm_line makes one from numbers)."""
    kept = []
    blocks = list(blocks)
    made = ffi.new("myb_ocm_block_data[]", max(len(blocks), 1))
    for k, block in enumerate(blocks):
        made[k].kind = _kind(block["kind"])
        pairs, count = _pairs(block.get("pairs", ()), kept)
        made[k].pair_count = count
        made[k].pairs = pairs
        lines = [ffi.new("char[]", str(line).encode("utf-8")) for line in block.get("lines", ())]
        kept += lines
        if lines:
            table = ffi.new("char*[]", lines)
            kept.append(table)
            made[k].lines = table
        made[k].line_count = len(lines)
    data = ffi.new("myb_ocm_data*")
    name = ffi.new("char[]", str(originator).encode("utf-8"))
    data.originator = name
    data.creation_date = _calendar(creation_date)[0]
    pairs, count = _pairs(meta, kept)
    data.meta_count = count
    data.meta = pairs
    data.block_count = len(blocks)
    data.blocks = made
    return _text(lambda buffer, size, length: lib.myb_ocm_write(
        data, _notation(notation), buffer, size, length))


def write_ocm_line(time, values):
    """A line of data (myb_ocm_format_line): time a calendar label or seconds from
    EPOCH_TZERO."""
    numbers, count = _series(values, 1, "values")
    if isinstance(time, (int, float)):
        epoch, offset = ffi.NULL, float(time)
    else:
        epoch, offset = _calendar(time), 0.0
    return _text(lambda buffer, size, length: lib.myb_ocm_format_line(
        epoch, offset, numbers, count, buffer, size, length))


# ---- CK and SCLK -------------------------------------------------------------------------------

class Sclk(_Owner):
    """A spacecraft clock of type 1 (myb_sclk) from the variables _<id> of a pool."""

    _destroy = staticmethod(lib.myb_sclk_destroy)

    def __init__(self, pool, id):
        out = ffi.new("myb_sclk**")
        diag = ffi.new("myb_diag*")
        _check(lib.myb_sclk_create(pool.pointer, int(id), ffi.NULL, out, diag), diag)
        super().__init__(out[0])

    @property
    def time_system(self):
        """1 = TDB, 2 = TT."""
        return lib.myb_sclk_time_system(self.pointer)

    @property
    def partitions(self):
        return lib.myb_sclk_partitions(self.pointer)

    def to_time(self, ticks):
        out = ffi.new("myb_time*")
        _check(lib.myb_sclk_to_time(self.pointer, float(ticks), out))
        return out[0]

    def from_time(self, instant):
        out = ffi.new("double*")
        _check(lib.myb_sclk_from_time(self.pointer, instant, out))
        return out[0]

    def encode(self, text):
        raw = text.encode("ascii")
        out = ffi.new("double*")
        _check(lib.myb_sclk_encode(self.pointer, raw, len(raw), out))
        return out[0]

    def decode(self, ticks):
        """The clock string of encoded ticks (myb_sclk_decode)."""
        length = ffi.new("size_t*")
        size = 64
        while True:
            text = ffi.new("char[]", size)
            status = lib.myb_sclk_decode(self.pointer, float(ticks), text, size, length)
            if status != lib.MYB_E_OVERFLOW or length[0] < size:
                break
            size = length[0] + 1
        _check(status)
        return ffi.string(text, length[0]).decode("ascii")


class Ck(_Owner):
    """A set of CK files (myb_ck)."""

    _destroy = staticmethod(lib.myb_ck_destroy)

    def __init__(self, *paths):
        out = ffi.new("myb_ck**")
        _check(lib.myb_ck_create(ffi.NULL, out))
        super().__init__(out[0])
        for path in paths:
            self.load_file(path)

    def load_file(self, path):
        diag = ffi.new("myb_diag*")
        _check(lib.myb_ck_load_file(self.pointer, str(path).encode("utf-8"), diag), diag)

    def load_memory(self, data):
        image = ffi.from_buffer(bytes(data))
        diag = ffi.new("myb_diag*")
        _check(lib.myb_ck_load_memory(self.pointer, image, len(image), diag), diag)
        self._kept.append(image)

    def segments(self):
        info = ffi.new("myb_ck_segment_info*")
        found = []
        for index in range(lib.myb_ck_segment_count(self.pointer)):
            lib.myb_ck_segment(self.pointer, index, info)
            found.append({"instrument": info.instrument, "frame": info.frame, "type": info.type,
                          "has_rates": info.has_rates != 0, "start": info.start,
                          "stop": info.stop, "file": info.file,
                          "name": ffi.string(info.name).decode("ascii", "replace")})
        return found

    def pointing(self, instrument, ticks, tolerance=0.0, need_rate=False):
        """The orientation of an instrument at ticks of its clock (myb_ck_pointing_at)."""
        out = ffi.new("myb_ck_pointing*")
        _check(lib.myb_ck_pointing_at(self.pointer, int(instrument), float(ticks),
                                      float(tolerance), 1 if need_rate else 0, out))
        return CkPointing(tuple(out.q), tuple(out.rate), out.has_rate != 0, out.ticks, out.frame,
                          _rotation_of(out.rotation))

    def orientation(self, sclk, instrument, base=0, tolerance=0.0):
        """An instrument frame as a provider of orientations (myb_orientation_ck)."""
        context = ffi.new("myb_ck_orientation*", {"ck": self.pointer, "sclk": sclk.pointer,
                                                  "instrument": int(instrument),
                                                  "base": int(base),
                                                  "tolerance": float(tolerance)})
        return Provider(lib.myb_orientation_ck, context, self, sclk)


def _ck_header(instrument, frame, name):
    text = ffi.new("char[]", name.encode("ascii"))
    return ffi.new("myb_ck_segment_header*", {"instrument": int(instrument),
                                              "frame": int(frame), "name": text}), text


_packet_sizes = {0: 8, 1: 4, 2: 14, 3: 7}


def _ck_packets(subtype, window, seconds_per_tick, ticks, packets, kept):
    times, count = _series(ticks, 1, "ticks")
    flat = _flat(packets)
    size = _packet_sizes.get(int(subtype))
    if size is not None and len(flat) != size * count:
        raise ValueError("subtype %d: %d numbers a packet, one packet a tick" % (subtype, size))
    values = ffi.new("double[]", flat or 1)
    kept += [times, values]
    return ffi.new("myb_ck_packets*", {"subtype": int(subtype), "window": int(window),
                                       "seconds_per_tick": float(seconds_per_tick),
                                       "ticks": times, "packets": values, "count": count})


class CkWriter(_Owner):
    """Segments collected in memory and laid out as one CK file (myb_ck_writer): types 1 and
    3 (pointings), 2 (intervals), 4 (Chebyshev), 5 (packets), 6 (mini-segments)."""

    _destroy = staticmethod(lib.myb_ck_writer_destroy)

    def __init__(self, internal_name=None):
        out = ffi.new("myb_ck_writer**")
        name = ffi.NULL if internal_name is None else ffi.new("char[]",
                                                              internal_name.encode("ascii"))
        _check(lib.myb_ck_writer_create(name, ffi.NULL, out))
        super().__init__(out[0])

    def add_pointings(self, instrument, frame, name, type, ticks, quaternions, rates=None,
                      interval_starts=None):
        header, text = _ck_header(instrument, frame, name)
        times, count = _series(ticks, 1, "ticks")
        q = _doubles(quaternions, 4 * count, "quaternions")
        w = ffi.NULL if rates is None else _doubles(rates, 3 * count, "rates")
        starts, n = (ffi.NULL, 0) if interval_starts is None else _series(interval_starts, 1,
                                                                          "interval_starts")
        _check(lib.myb_ck_writer_add_pointings(self.pointer, header, int(type), times, q, w,
                                               count, starts, n))

    def add_intervals(self, instrument, frame, name, starts, stops, quaternions, rates,
                      seconds_per_tick):
        header, text = _ck_header(instrument, frame, name)
        begin, count = _series(starts, 1, "starts")
        _check(lib.myb_ck_writer_add_intervals(
            self.pointer, header, begin, _doubles(stops, count, "stops"),
            _doubles(quaternions, 4 * count, "quaternions"), _doubles(rates, 3 * count, "rates"),
            _doubles(seconds_per_tick, count, "seconds_per_tick"), count))

    def add_chebyshev(self, instrument, frame, name, starts, stops, counts, coefficients,
                      has_rates):
        """counts: seven a interval (q0 .. q3 and the three rates)."""
        header, text = _ck_header(instrument, frame, name)
        begin, count = _series(starts, 1, "starts")
        numbers = [int(c) for c in _flat(counts)]
        if len(numbers) != 7 * count:
            raise ValueError("counts: seven an interval")
        values = _flat(coefficients)
        if len(values) != sum(numbers):
            raise ValueError("coefficients: %d numbers, as the counts say" % sum(numbers))
        _check(lib.myb_ck_writer_add_chebyshev(
            self.pointer, header, begin, _doubles(stops, count, "stops"),
            ffi.new("int[]", numbers or 1), ffi.new("double[]", values or 1), count,
            1 if has_rates else 0))

    def add_packets(self, instrument, frame, name, subtype, window, seconds_per_tick, ticks,
                    packets, has_rates, start, stop, interval_starts=None):
        header, text = _ck_header(instrument, frame, name)
        kept = []
        made = _ck_packets(subtype, window, seconds_per_tick, ticks, packets, kept)
        starts, n = (ffi.NULL, 0) if interval_starts is None else _series(interval_starts, 1,
                                                                          "interval_starts")
        _check(lib.myb_ck_writer_add_packets(self.pointer, header, made, 1 if has_rates else 0,
                                             float(start), float(stop), starts, n))

    def add_minisegments(self, instrument, frame, name, minisegments, boundaries, has_rates,
                         select_last, start, stop):
        """minisegments: mappings with subtype, window, seconds_per_tick, ticks, packets."""
        header, text = _ck_header(instrument, frame, name)
        kept = []
        minisegments = list(minisegments)
        made = ffi.new("myb_ck_packets[]", max(len(minisegments), 1))
        for k, m in enumerate(minisegments):
            made[k] = _ck_packets(m["subtype"], m["window"], m["seconds_per_tick"], m["ticks"],
                                  m["packets"], kept)[0]
        edges = _doubles(boundaries, len(minisegments) + 1, "boundaries")
        _check(lib.myb_ck_writer_add_minisegments(self.pointer, header, made, edges,
                                                  len(minisegments), 1 if has_rates else 0,
                                                  1 if select_last else 0, float(start),
                                                  float(stop)))

    @property
    def size(self):
        return lib.myb_ck_writer_size(self.pointer)

    def write(self):
        size = self.size
        buffer = ffi.new("char[]", max(size, 1))
        _check(lib.myb_ck_writer_write(self.pointer, buffer, size))
        return ffi.buffer(buffer, size)[:]

    def write_file(self, path):
        _check(lib.myb_ck_writer_write_file(self.pointer, str(path).encode("utf-8")))


# ---- Two-line elements (ext/sgp4) --------------------------------------------------------------

def _sgp4():
    return extension("sgp4")


def tle_parse(line1, line2, checksum=True):
    """The two lines of a TLE as a mapping of the members of myb_tle (myb_tle_parse)."""
    sgp4 = _sgp4()
    out = ffi.new("myb_tle*")
    diag = ffi.new("myb_diag*")
    flags = 0 if checksum else sgp4.MYB_TLE_NO_CHECKSUM
    _check(sgp4.myb_tle_parse(line1.encode("ascii"), line2.encode("ascii"), flags, out, diag),
           diag)
    return _from_struct(out, "myb_tle")


def tle_format(tle):
    """(line1, line2) of a TLE given as a mapping (myb_tle_format)."""
    sgp4 = _sgp4()
    made = _to_struct(dict(tle), "myb_tle")
    one = ffi.new("char[70]")
    two = ffi.new("char[70]")
    _check(sgp4.myb_tle_format(made, one, two))
    return ffi.string(one).decode("ascii"), ffi.string(two).decode("ascii")


def tle_epoch(tle):
    """The epoch of a TLE as a UTC label (myb_tle_epoch)."""
    sgp4 = _sgp4()
    out = ffi.new("myb_calendar*")
    sgp4.myb_tle_epoch(_to_struct(dict(tle), "myb_tle"), out)
    return _label(out)


def tle_from_omm(omm):
    sgp4 = _sgp4()
    out = ffi.new("myb_tle*")
    _check(sgp4.myb_tle_from_omm(_to_struct(dict(omm), "myb_omm"), out))
    return _from_struct(out, "myb_tle")


def omm_from_tle(tle):
    sgp4 = _sgp4()
    out = ffi.new("myb_omm*")
    _check(sgp4.myb_omm_from_tle(_to_struct(dict(tle), "myb_tle"), out))
    return _omm_of(out)


_gravities = {"wgs72old": 0, "wgs72": 1, "wgs84": 2}


def tle_from_state(r, v, gravity="wgs72"):
    """A TLE whose mean elements are the osculating ones of (r, v) in km and km/s, TEME
    (myb_tle_from_state)."""
    sgp4 = _sgp4()
    out = ffi.new("myb_tle*")
    grav = _gravities[gravity] if isinstance(gravity, str) else int(gravity)
    _check(sgp4.myb_tle_from_state(_doubles(r, 3, "r"), _doubles(v, 3, "v"), grav, out))
    return _from_struct(out, "myb_tle")


def tle_fit(minutes, positions, tle, bstar=False, max_iterations=0, tolerance=0.0,
            gravity="wgs72", mode=0):
    """A TLE fitted by SGP4 to positions (km, TEME) at minutes from its epoch, from the
    first guess `tle` (myb_tle_fit): (tle, rms, iterations)."""
    sgp4 = _sgp4()
    times, count = _series(minutes, 1, "minutes")
    points = _doubles(positions, 3 * count, "positions")
    options = ffi.new("myb_tle_fit_options*")
    options.bstar = 1 if bstar else 0
    options.max_iterations = int(max_iterations)
    options.tolerance = float(tolerance)
    options.gravity = _gravities[gravity] if isinstance(gravity, str) else int(gravity)
    options.mode = int(mode)
    made = _to_struct(dict(tle), "myb_tle")
    rms = ffi.new("double*")
    iterations = ffi.new("int*")
    _check(sgp4.myb_tle_fit(times, points, count, options, made, rms, iterations))
    return _from_struct(made, "myb_tle"), rms[0], iterations[0]
