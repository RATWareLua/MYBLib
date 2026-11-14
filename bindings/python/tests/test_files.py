# SPDX-License-Identifier: BSD-3-Clause
"""Files through the package against the same C calls through cffi, bit for bit (the bytes
of a file, the numbers read): numbers and tables, kernels, the SPK writer, SP3, the OPM,
OMM and OCM messages, CK and SCLK with the CK writer, two-line elements; failures as
myb.Error with the status and the line of a parser."""
import math
import os

import myb
from myb import ffi, lib

from common import TEST_DATA, Case, doubles, ok

OPM = """CCSDS_OPM_VERS = 2.0
CREATION_DATE = 2000-06-03T05:33:00.123
ORIGINATOR = GSOC
OBJECT_NAME = EUTELSAT W4
OBJECT_ID = 2000-028A
CENTER_NAME = EARTH
REF_FRAME = TOD
TIME_SYSTEM = UTC
EPOCH = 2006-06-03T00:00:00.000
X = 6655.9942 [km]
Y = -40218.5751 [km]
Z = -82.9177 [km]
X_DOT = 3.11548208 [km/s]
Y_DOT = 0.47042605 [km/s]
Z_DOT = -0.00101495 [km/s]
MASS = 1913.000 [kg]
SOLAR_RAD_AREA = 10.000 [m**2]
SOLAR_RAD_COEFF = 1.300
DRAG_AREA = 10.000 [m**2]
DRAG_COEFF = 2.300
MAN_EPOCH_IGNITION = 2000-06-03T09:00:34.1
MAN_DURATION = 132.60 [s]
MAN_DELTA_MASS = -18.418 [kg]
MAN_REF_FRAME = EME2000
MAN_DV_1 = -0.02325700 [km/s]
MAN_DV_2 = 0.01683160 [km/s]
MAN_DV_3 = -0.00893444 [km/s]
"""


def written(write):
    """The text of a writer of the C interface (the length first, then the text)."""
    length = ffi.new("size_t*")
    write(ffi.NULL, 0, length)
    buffer = ffi.new("char[]", length[0] + 1)
    ok(write(buffer, length[0] + 1, length))
    return ffi.string(buffer, length[0]).decode("utf-8")


class Numbers(Case):
    def test_numbers_tables_and_kernels(self):
        value, used = myb.parse_double("6.02214076e23 mol")
        raw = ffi.new("double*")
        consumed = ffi.new("size_t*")
        ok(lib.myb_parse_double(b"6.02214076e23 mol", 17, raw, consumed))
        self.assertSame((value, used), (raw[0], consumed[0]))
        self.assertStatus(lib.MYB_E_FORMAT, myb.parse_double, "x1")
        self.assertEqual(myb.format_double(0.1), "0.1")
        self.assertEqual(myb.format_fixed(math.pi, 5), "3.14159")
        self.assertEqual(myb.format_scientific(123456.0, 3), "1.23e+05")
        self.assertStatus(lib.MYB_E_DOMAIN, myb.format_fixed, 1.0, 99)
        text = "t,x,y\n0,1.5,2\n1,,3.25\n"
        table = myb.Table.load(text)
        made = ffi.new("myb_table**")
        diag = ffi.new("myb_diag*")
        ok(lib.myb_table_load(text.encode(), len(text), ffi.NULL, made, diag))
        self.assertEqual(table.names, ["t", "x", "y"])
        self.assertSame(table.values[0], doubles(lib.myb_table_row(made[0], 0), 3))
        self.assertTrue(math.isnan(table.values[1][1]))
        lib.myb_table_destroy(made[0])
        out = myb.write_table([(0.0, 1.5, 2.0), (1.0, float("nan"), 3.25)], ["t", "x", "y"])
        names = [ffi.new("char[]", n) for n in (b"t", b"x", b"y")]
        values = ffi.new("double[6]", [0.0, 1.5, 2.0, 1.0, float("nan"), 3.25])
        self.assertEqual(out, written(lambda b, s, n: lib.myb_table_write(
            ffi.new("char*[]", names), 3, values, 2, b, s, n)))
        error = self.assertStatus(lib.MYB_E_FORMAT, myb.Table.load, "a,b\n1,2\n3,x\n")
        self.assertIn("line 3", str(error))
        pool = myb.KernelPool()
        pool.load_text("\\begindata\nA = ( 1.5 2.5 )\nB = ( 'one' 'two' )\n\\begintext\n")
        self.assertEqual(pool.describe("A"), ("numbers", 2))
        self.assertEqual(pool.numbers("A"), (1.5, 2.5))
        self.assertEqual(pool.strings("B"), ["one", "two"])
        self.assertStatus(lib.MYB_E_NOT_FOUND, pool.numbers, "C")
        self.assertStatus(lib.MYB_E_FORMAT, pool.numbers, "B")


class Spk(Case):
    def test_the_writer(self):
        epochs = [6.0e8 + 60.0 * k for k in range(8)]
        states = [(7.0e6 * math.cos(1e-3 * 60 * k), 7.0e6 * math.sin(1e-3 * 60 * k), 1.0e3 * k,
                   -7.0e3 * math.sin(1e-3 * 60 * k), 7.0e3 * math.cos(1e-3 * 60 * k), 1.0)
                  for k in range(8)]
        writer = myb.SpkWriter("python")
        writer.add_states(-77, 399, 1, "states", 13, 7, epochs, states)
        image = writer.write()
        raw = ffi.new("myb_spk_writer**")
        ok(lib.myb_spk_writer_create(b"python", ffi.NULL, raw))
        name = ffi.new("char[]", b"states")
        header = ffi.new("myb_spk_segment_header*", {"target": -77, "center": 399, "frame": 1,
                                                     "name": name})
        ok(lib.myb_spk_writer_add_states(raw[0], header, 13, 7, ffi.new("double[]", epochs),
                                         ffi.new("double[]", [v for s in states for v in s]), 8))
        size = lib.myb_spk_writer_size(raw[0])
        buffer = ffi.new("char[]", size)
        ok(lib.myb_spk_writer_write(raw[0], buffer, size))
        self.assertEqual(image, ffi.buffer(buffer, size)[:])
        lib.myb_spk_writer_destroy(raw[0])
        eph = myb.Ephemeris()
        eph.load_memory(image)
        back = eph.state_many(-77, 399, epochs)
        self.assertLess(max(abs(a - b) for s, t in zip(back, states) for a, b in zip(s, t)),
                        1e-8)
        # a provider in Python fitted by Chebyshev polynomials, words copied from another file
        error = writer.add_chebyshev(-78, 399, 1, "fitted", 3, lambda s: eph.state(
            -77, 399, myb.time_from_seconds(s)), epochs[0], epochs[-1], 2, 9, max_error=True)
        self.assertGreaterEqual(error, 0.0)
        other = myb.Ephemeris(os.path.join(TEST_DATA, "spk", "t02_03.bsp"))
        first = other.segments()[0]
        writer.add_raw(first["target"], first["center"], first["frame"], "copy", first["type"],
                       first["start"], first["stop"], other.segment_data(0))
        copy = myb.Ephemeris()
        copy.load_memory(writer.write())
        middle = myb.time_from_seconds(0.5 * (first["start"] + first["stop"]))
        self.assertSame(copy.state(first["target"], first["center"], middle),
                        other.state(first["target"], first["center"], middle))
        self.assertStatus(lib.MYB_E_DOMAIN, writer.add_states, -77, 399, 1, "bad", 13, 7,
                          list(reversed(epochs)), states)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.SpkWriter().write)


class Sp3(Case):
    def test_reading_and_writing(self):
        path = os.path.join(TEST_DATA, "sp3", "ilrsa.lageos1.250607.excerpt.sp3")
        sp3 = myb.Sp3.load_file(path)
        made = ffi.new("myb_sp3**")
        diag = ffi.new("myb_diag*")
        ok(lib.myb_sp3_load_file(path.encode(), ffi.NULL, made, diag))
        r = ffi.new("double[3]")
        v = ffi.new("double[3]")
        self.assertEqual(lib.myb_sp3_state(made[0], 0, 3, r, v), 3)
        self.assertSame(sp3.state(0, 3), (doubles(r, 3), doubles(v, 3)))
        self.assertEqual((sp3.satellites, sp3.time_system, sp3.frame), (["L51"], "UTC", "SLR20"))
        epochs = sp3.epochs
        positions = [[sp3.state(0, e)[0]] for e in range(len(epochs))]
        velocities = [[sp3.state(0, e)[1]] for e in range(len(epochs))]
        text = myb.write_sp3(sp3.satellites, epochs, 360.0, positions, velocities, version="c",
                             file_type="L", time_system="UTC", frame="SLR20")
        again = myb.Sp3.load(text)
        self.assertEqual(again.epochs, epochs)
        self.assertSame([again.state(0, e) for e in range(len(epochs))],
                        [sp3.state(0, e) for e in range(len(epochs))])
        lib.myb_sp3_destroy(made[0])
        # clocks, accuracies, sigmas and flags written and read back
        text = myb.write_sp3(["L51"], epochs[:3], 360.0, positions[:3],
                             clocks=[[1e-6], [2e-6], [None]], accuracy=[0.002],
                             sigmas=[[(0.01,) * 8]] * 3,
                             flags=[[0], [lib.MYB_SP3_MANEUVER], [0]], file_type="L",
                             time_system="UTC")
        back = myb.Sp3.load(text)
        self.assertEqual((back.clock(0, 0), back.clock(0, 2)), ((1e-6, None), (None, None)))
        self.assertEqual((back.accuracy(0), back.flags(0, 1)), (0.002, lib.MYB_SP3_MANEUVER))
        self.assertEqual(len(back.sigma(0, 1)), 8)
        self.assertStatus(lib.MYB_E_IO, myb.Sp3.load_file, path + ".none")
        self.assertStatus(lib.MYB_E_DOMAIN, myb.write_sp3, ["L51"], epochs[:2], -1.0,
                          positions[:2])


class Messages(Case):
    def test_an_ephemeris_message_with_covariances(self):
        matrix = [[float(r == c) * (r + 1) for c in range(6)] for r in range(6)]
        segment = {"object_name": "A", "object_id": "1", "center_name": "EARTH",
                   "ref_frame": "EME2000", "time_system": "UTC",
                   "epochs": [(2026, 1, 1, 0, 0, 0.0)], "states": [(7e6, 0.0, 0.0, 0.0, 7.5e3,
                                                                    0.0)],
                   "covariances": [((2026, 1, 1, 0, 0, 0.0), matrix)], "covariance_frame": "RTN"}
        text = myb.write_oem("myb", (2026, 1, 2, 0, 0, 0.0), [segment])
        self.assertEqual(myb.oem_header(text, "ORIGINATOR"), "myb")
        self.assertEqual(myb.read_oem(text)[0]["covariances"],
                         [((2026, 1, 1, 0, 0, 0.0), "RTN", [tuple(row) for row in matrix])])

    def test_opm_omm_ocm(self):
        opm = myb.read_opm(OPM)
        made = ffi.new("myb_opm*")
        diag = ffi.new("myb_diag*")
        ok(lib.myb_opm_parse(OPM.encode(), len(OPM), made, diag))
        self.assertSame(opm["state"], tuple(made.state))
        self.assertEqual(opm["object_name"], "EUTELSAT W4")
        self.assertSame(opm["maneuvers"][0]["dv"], tuple(made.maneuvers[0].dv))
        for notation in ("kvn", "xml"):
            text = myb.write_opm(opm, notation)
            self.assertEqual(text, written(lambda b, s, n: lib.myb_opm_write_as(
                made, lib.MYB_ODM_XML if notation == "xml" else lib.MYB_ODM_KVN, b, s, n)))
            self.assertSame(myb.read_opm(text)["state"], opm["state"])
        error = self.assertStatus(lib.MYB_E_FORMAT, myb.read_opm, "CCSDS_OPM_VERS = 2.0\nX = 1\n")
        self.assertIn("CREATION_DATE", error.message)
        omm = {"originator": "MYB", "object_name": "SAT", "object_id": "2026-001A",
               "center_name": "EARTH", "ref_frame": "TEME", "time_system": "UTC",
               "mean_element_theory": "SGP4", "creation_date": (2026, 11, 14, 0, 0, 0.0),
               "epoch": (2026, 11, 13, 12, 0, 0.5), "has_mean_motion": 1, "mean_motion": 15.5,
               "eccentricity": 0.001, "inclination": 51.6, "ra_of_asc_node": 10.0,
               "arg_of_pericenter": 20.0, "mean_anomaly": 30.0}
        text = myb.write_omm(omm)
        back = myb.read_omm(text)
        self.assertEqual(back["object_id"], "2026-001A")
        self.assertSame(back["mean_motion"], 15.5)
        line = myb.write_ocm_line(60.0, (7000.0, 0.5, 0.25, 1.0, 7.5, 0.0))
        ocm_text = myb.write_ocm("MYB", (2026, 11, 14, 0, 0, 0.0),
                                 {"TIME_SYSTEM": "UTC", "EPOCH_TZERO": "2026-11-14T00:00:00"},
                                 [{"kind": "traj", "pairs": {"CENTER_NAME": "EARTH",
                                                             "TRAJ_REF_FRAME": "GCRF",
                                                             "TRAJ_TYPE": "CARTPV"},
                                   "lines": [line]}])
        ocm = myb.Ocm.load(ocm_text)
        self.assertEqual((ocm.meta("TIME_SYSTEM"), ocm.block_count("traj"),
                          ocm.line_count("traj", 0)), ("UTC", 1, 1))
        self.assertEqual(ocm.line("traj", 0, 0), (60.0, (7000.0, 0.5, 0.25, 1.0, 7.5, 0.0)))
        self.assertSame(ocm.state(0, 0), (60.0, (7.0e6, 500.0, 250.0, 1000.0, 7500.0, 0.0)))
        self.assertEqual(ocm.value("traj", 0, "TRAJ_TYPE"), "CARTPV")
        self.assertIsNone(ocm.header("NO_SUCH_KEY"))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.write_ocm, "MYB", (2026, 11, 14, 0, 0, 0.0),
                          {"TIME_SYSTEM": "UTC"}, [])


class Orientation(Case):
    def test_ck_sclk_and_the_writer(self):
        pool = myb.KernelPool(os.path.join(TEST_DATA, "ck", "test.tsc"))
        clock = myb.Sclk(pool, 9099)
        self.assertEqual((clock.time_system, clock.partitions), (1, 2))
        self.assertSame(myb.seconds(clock.to_time(0.0)), 599996400.0)
        self.assertEqual(clock.decode(0.0), "1/0000000:000:0")
        self.assertEqual(clock.encode("1/0000000:000:0"), 0.0)
        self.assertSame(clock.from_time(clock.to_time(4096000.0)), 4096000.0)
        self.assertStatus(lib.MYB_E_NOT_FOUND, myb.Sclk, pool, 9097)
        ck = myb.Ck(os.path.join(TEST_DATA, "ck", "test.bc"))
        self.assertEqual(len(ck.segments()), 7)
        pointing = ck.pointing(-9099002, 6.2e6)
        raw_ck = ffi.new("myb_ck**")
        diag = ffi.new("myb_diag*")
        ok(lib.myb_ck_create(ffi.NULL, raw_ck))
        ok(lib.myb_ck_load_file(raw_ck[0], os.path.join(TEST_DATA, "ck", "test.bc").encode(),
                                diag))
        made = ffi.new("myb_ck_pointing*")
        ok(lib.myb_ck_pointing_at(raw_ck[0], -9099002, 6.2e6, 0.0, 0, made))
        self.assertSame(pointing, myb.CkPointing(tuple(made.q), tuple(made.rate),
                                                 made.has_rate != 0, made.ticks, made.frame,
                                                 myb.Rotation(tuple(made.rotation.m),
                                                              tuple(made.rotation.dm))))
        lib.myb_ck_destroy(raw_ck[0])
        self.assertStatus(lib.MYB_E_RANGE, ck.pointing, -1, 6.1e6)
        self.assertStatus(lib.MYB_E_DOMAIN, ck.pointing, -9099001, 6.1e6, -1.0)
        provider = ck.orientation(clock, -9099001, base=1, tolerance=1.0e9)
        seen = myb.orientation_many(provider, [clock.to_time(6.1e6)])
        self.assertEqual(len(seen[0].rotation.m), 9)
        # the writer: pointings of type 3, the same bytes as the C calls give
        ticks = [1000.0 * k for k in range(6)]
        quaternions = [(math.cos(0.05 * k), 0.0, 0.0, math.sin(0.05 * k)) for k in range(6)]
        writer = myb.CkWriter()
        writer.add_pointings(-77001, 1, "turning", 3, ticks, quaternions,
                             rates=[(0.0, 0.0, 1e-4)] * 6)
        image = writer.write()
        raw = ffi.new("myb_ck_writer**")
        ok(lib.myb_ck_writer_create(ffi.NULL, ffi.NULL, raw))
        name = ffi.new("char[]", b"turning")
        header = ffi.new("myb_ck_segment_header*", {"instrument": -77001, "frame": 1,
                                                    "name": name})
        ok(lib.myb_ck_writer_add_pointings(raw[0], header, 3, ffi.new("double[]", ticks),
                                           ffi.new("double[]", [v for q in quaternions
                                                                for v in q]),
                                           ffi.new("double[]", [0.0, 0.0, 1e-4] * 6), 6,
                                           ffi.NULL, 0))
        size = lib.myb_ck_writer_size(raw[0])
        buffer = ffi.new("char[]", size)
        ok(lib.myb_ck_writer_write(raw[0], buffer, size))
        self.assertEqual(image, ffi.buffer(buffer, size)[:])
        lib.myb_ck_writer_destroy(raw[0])
        back = myb.Ck()
        back.load_memory(image)
        self.assertEqual(back.segments()[0]["instrument"], -77001)
        self.assertLess(max(abs(a - b) for a, b in zip(back.pointing(-77001, 2000.0).q,
                                                       quaternions[2])), 1e-15)
        self.assertStatus(lib.MYB_E_DOMAIN, writer.add_pointings, -77001, 1, "bad", 3,
                          list(reversed(ticks)), quaternions)
        # the other types of segments
        more = myb.CkWriter("more")
        more.add_intervals(-1, 1, "intervals", (0.0, 100.0), (50.0, 150.0), quaternions[:2],
                           [(0.0, 0.0, 1e-3)] * 2, (1.0, 1.0))
        more.add_chebyshev(-2, 1, "chebyshev", (0.0,), (100.0,), (1, 1, 1, 1, 0, 0, 0),
                           (1.0, 0.0, 0.0, 0.0), False)
        packets = [(math.cos(0.01 * t), 0.0, 0.0, math.sin(0.01 * t)) for t in ticks[:4]]
        more.add_packets(-3, 1, "packets", 1, 2, 1.0, ticks[:4], packets, False, 0.0, 3000.0)
        more.add_minisegments(-4, 1, "mini", [{"subtype": 1, "window": 2,
                                               "seconds_per_tick": 1.0, "ticks": ticks[:4],
                                               "packets": packets}],
                              (0.0, 3000.0), False, False, 0.0, 3000.0)
        with self.assertRaises(ValueError):
            more.add_packets(-3, 1, "short", 0, 2, 1.0, ticks[:4], packets, False, 0.0, 3000.0)
        path = os.path.join(os.environ.get("TEMP", os.environ.get("TMPDIR", "/tmp")),
                            "myb-python-test-%d.bc" % os.getpid())
        more.write_file(path)
        try:
            self.assertEqual([s["type"] for s in myb.Ck(path).segments()], [2, 4, 5, 6])
        finally:
            os.remove(path)


class TwoLines(Case):
    LINE1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    LINE2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"

    def test_two_line_elements(self):
        try:
            sgp4 = myb.extension("sgp4")
        except ImportError:
            self.skipTest("libmyb_sgp4 is not next to the library")
        tle = myb.tle_parse(self.LINE1, self.LINE2)
        made = ffi.new("myb_tle*")
        diag = ffi.new("myb_diag*")
        ok(sgp4.myb_tle_parse(self.LINE1.encode(), self.LINE2.encode(), 0, made, diag))
        self.assertSame((tle["inclination"], tle["mean_motion"], tle["catalog_number"]),
                        (made.inclination, made.mean_motion, made.catalog_number))
        self.assertEqual(myb.tle_format(tle), (self.LINE1, self.LINE2))
        epoch = myb.tle_epoch(tle)
        self.assertEqual(epoch[:3], (2000, 6, 27))
        omm = myb.omm_from_tle(tle)
        self.assertSame(myb.tle_from_omm(omm)["mean_motion"], tle["mean_motion"])
        self.assertStatus(lib.MYB_E_FORMAT, myb.tle_parse, self.LINE1, self.LINE2[:-1] + "0")
        self.assertEqual(myb.tle_parse(self.LINE1, self.LINE2[:-1] + "0",
                                       checksum=False)["catalog_number"], 5)

