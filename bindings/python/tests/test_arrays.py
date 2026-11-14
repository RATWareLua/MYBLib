# SPDX-License-Identifier: BSD-3-Clause
"""The array forms over sequences against the scalar C calls made one by one through cffi:
time, the time scales, SPK states, Kepler and the elements, the rotation of the Earth and of
the bodies, PCK, orientation providers, the graph of frames — bit for bit; a failing element
is named by Error.index."""
import array
import math
import os

import myb
from myb import ffi, lib

from common import ELEMENTS, GM, TEST_DATA, Case, kernels, ok

SECONDS = [6.32e8 + 2.5e6 * k + 0.375 * k for k in range(12)]


def earth_pool():
    return myb.KernelPool(*kernels("time/time-scales.tk", "time/tdb-tt-fb1990.tk",
                                   "time/utc-tai.tk", "frames/iau2006-cip.tk",
                                   "frames/iau1980.tk", "frames/eop-tides.tk"))


def raw_rotation(rotation):
    return myb.Rotation(tuple(rotation.m), tuple(rotation.dm))


class Time(Case):
    def test_instants_and_seconds(self):
        many = myb.time_from_seconds_many(array.array("d", SECONDS))
        self.assertSame(many, [lib.myb_time_from_seconds(s) for s in SECONDS])
        self.assertSame(myb.seconds_many(many), [lib.myb_time_seconds(t) for t in many])
        # seconds where instants are asked for take the same way in
        self.assertSame(myb.seconds_many(SECONDS), myb.seconds_many(many))
        self.assertEqual(myb.time_from_seconds_many([]), [])

    def test_calendars(self):
        labels = [(2024, 2, 29, 23, 59, 59.75), (1999, 12, 31, 0, 0, 0.125),
                  (2030, 7, 1, 12, 30, 1.0)]
        many = myb.calendar_to_time_many(labels)
        one = ffi.new("myb_time*")
        for label, instant in zip(labels, many):
            ok(lib.myb_calendar_to_time(ffi.new("myb_calendar*", dict(zip(
                ("year", "month", "day", "hour", "minute", "second"), label))), one))
            self.assertSame(instant, one[0])
        self.assertEqual(myb.time_to_calendar_many(many), [myb.time_to_calendar(t) for t in many])
        error = self.assertStatus(lib.MYB_E_DOMAIN, myb.calendar_to_time_many,
                                  labels[:2] + [(2030, 13, 1, 0, 0, 0.0)])
        self.assertEqual(error.index, 2)

    def test_time_scales(self):
        pool = myb.KernelPool(*kernels("time/time-scales.tk", "time/tdb-tt-fb1990.tk",
                                       "time/utc-tai.tk"))
        system = myb.TimeSystem(pool)
        instants = myb.time_from_seconds_many(SECONDS)
        tdb = system.convert_many(instants, "TT", "TDB")
        tt = ffi.new("myb_scale*")
        to = ffi.new("myb_scale*")
        ok(lib.myb_timesys_scale(system.pointer, b"TT", tt))
        ok(lib.myb_timesys_scale(system.pointer, b"TDB", to))
        out = ffi.new("myb_time*")
        for instant, converted in zip(instants, tdb):
            ok(lib.myb_timesys_convert(system.pointer, instant, tt[0], to[0], out))
            self.assertSame(converted, out[0])
        labels = system.to_utc_many(instants, "TT")
        self.assertEqual(labels, [system.to_utc(t, "TT") for t in instants])
        self.assertSame(system.from_utc_many(labels, "TT"),
                        [system.from_utc(label, "TT") for label in labels])
        self.assertStatus(lib.MYB_E_NOT_FOUND, system.convert_many, instants, "TT", "NONE")


class Ephemerides(Case):
    def test_states_at_many_epochs(self):
        eph = myb.Ephemeris(os.path.join(TEST_DATA, "spk", "t02_03.bsp"))
        first = eph.segments()[0]
        times = [first["start"] + (first["stop"] - first["start"]) * (k + 0.5) / 12
                 for k in range(12)]
        instants = myb.time_from_seconds_many(times)
        states, frames = eph.state_many(first["target"], first["center"], instants, frames=True)
        y = ffi.new("double[6]")
        frame = ffi.new("int32_t*")
        for instant, state, code in zip(instants, states, frames):
            ok(lib.myb_spk_state(eph.pointer, first["target"], first["center"], instant, y,
                                 frame))
            self.assertSame(state, tuple(y))
            self.assertEqual(code, frame[0])
        self.assertSame(eph.state_many(first["target"], first["center"], times), states)
        beyond = list(times)
        beyond[5] = first["stop"] + 1.0e9
        error = self.assertStatus(lib.MYB_E_RANGE, eph.state_many, first["target"],
                                  first["center"], beyond)
        self.assertEqual(error.index, 5)


class Kepler(Case):
    def test_elements_and_propagation(self):
        states = myb.state_from_elements_many(GM, [ELEMENTS._replace(nu=0.4 * k)
                                                   for k in range(6)])
        pv = ffi.new("myb_pv*")
        for k, state in enumerate(states):
            given = ffi.new("myb_elements*", ELEMENTS._replace(nu=0.4 * k)._asdict())
            ok(lib.myb_pv_from_elements(GM, given, pv))
            self.assertSame(state, tuple(pv.r) + tuple(pv.v))
        elements = myb.elements_from_state_many(GM, states)
        self.assertSame(elements, [myb.elements_from_state(GM, s) for s in states])
        equinoctial = myb.equinoctial_from_state_many(GM, states)
        self.assertSame(equinoctial, [myb.equinoctial_from_state(GM, s) for s in states])
        self.assertSame(myb.state_from_equinoctial_many(GM, equinoctial),
                        [myb.state_from_equinoctial(GM, e) for e in equinoctial])
        dts = [600.0 * k + 0.125 for k in range(10)]
        propagated = myb.kepler_propagate_many(GM, states[0], dts)
        self.assertSame(propagated, [myb.kepler_propagate(GM, states[0], dt) for dt in dts])
        later, stms = myb.kepler_stm_many(GM, states[0], dts)
        self.assertSame(later, propagated)
        stm = ffi.new("double[36]")
        for dt, matrix in zip(dts, stms):
            ok(lib.myb_kepler_stm(GM, myb._core._pv(states[0]), dt, pv, stm))
            self.assertSame(matrix, [tuple(stm[6 * r:6 * r + 6]) for r in range(6)])
        bad = [ELEMENTS, ELEMENTS._replace(p=-1.0), ELEMENTS]
        error = self.assertStatus(lib.MYB_E_DOMAIN, myb.state_from_elements_many, GM, bad)
        self.assertEqual(error.index, 1)


class Frames(Case):
    def test_the_rotation_of_the_earth(self):
        pool = earth_pool()
        system = myb.TimeSystem(pool)
        cip = myb.Cip(pool)
        equinox = myb.Equinox(pool)
        eop = myb.Eop(c04=os.path.join(TEST_DATA, "eop", "eopc04-excerpt.txt"))
        eop.set_tides(pool)
        tt = myb.time_from_seconds_many(SECONDS)
        xys, rates = cip.xys_many(tt)
        raw_xys = ffi.new("double[3]")
        raw_rates = ffi.new("double[3]")
        for instant, a, b in zip(tt, xys, rates):
            lib.myb_cip_xys(cip.pointer, instant, raw_xys, raw_rates)
            self.assertSame((a, b), (tuple(raw_xys), tuple(raw_rates)))
        values, ut1 = eop.at_many(system, tt)
        self.assertSame(list(zip(values, ut1)), [eop.at(system, t) for t in tt])
        rotations = cip.terrestrial_many(tt, ut1, values)
        out = ffi.new("myb_rotation*")
        for k, instant in enumerate(tt):
            given = ffi.new("myb_eop_values*", values[k]._asdict())
            lib.myb_cip_terrestrial(cip.pointer, instant, ut1[k], given, out)
            self.assertSame(rotations[k], raw_rotation(out))
        self.assertSame(cip.terrestrial_many(tt, ut1), [cip.terrestrial(a, b)
                                                         for a, b in zip(tt, ut1)])
        for frame in ("mod", "tod", "teme"):
            self.assertSame(equinox.rotation_many(frame, tt),
                            [equinox.rotation(frame, t) for t in tt])
        self.assertSame(equinox.teme_pef_many(ut1), [equinox.teme_pef(t) for t in ut1])
        # a series of states carried by them
        states = [(7.0e6 + k, 1.0e5, -2.0e5, 10.0, 7.5e3, 1.0e3) for k in range(len(tt))]
        carried = myb.rotation_apply_state_many(rotations, states)
        made = ffi.new("myb_rotation*")
        y = ffi.new("double[6]")
        for rotation, state, result in zip(rotations, states, carried):
            made.m = list(rotation.m)
            made.dm = list(rotation.dm)
            lib.myb_rotation_apply_state(made, ffi.new("double[6]", state), y)
            self.assertSame(result, tuple(y))
        x = ffi.new("double[3]")
        for rotation, state, result in zip(rotations, states, myb.rotation_apply_position_many(
                rotations, [s[:3] for s in states])):
            made.m = list(rotation.m)
            lib.myb_rotation_apply_position(made, ffi.new("double[3]", state[:3]), x)
            self.assertSame(result, tuple(x))
        late = list(tt)
        late[4] = myb.time_from_seconds(4.0e9)
        error = self.assertStatus(lib.MYB_E_RANGE, eop.at_many, system, late)
        self.assertEqual(error.index, 4)
        # the rotation with the observed EOP as a provider: in a graph, alone
        provider = cip.orientation(eop, system)
        graph = myb.FrameGraph()
        graph.add(3000, 1, provider)
        through = myb.orientation_many(provider, tt)
        self.assertSame([o.rotation for o in through], graph.rotation_many(1, 3000, tt))
        self.assertSame(through[0].ut, ut1[0])
        rate = ffi.new("double*")
        self.assertSame(cip.era(ut1[0]), (lib.myb_cip_era(cip.pointer, ut1[0], rate), rate[0]))
        self.assertSame(equinox.gmst(ut1[0]), (lib.myb_equinox_gmst(equinox.pointer, ut1[0],
                                                                     rate), rate[0]))
        self.assertStatus(lib.MYB_E_NOT_FOUND, myb.Cip, pool, "NOSUCH")

    def test_bodies_pck_providers_and_the_graph(self):
        mars = myb.Wgccre(myb.KernelPool(os.path.join(TEST_DATA, "wgccre", "pck00011.tpc")), 499)
        t = myb.time_from_seconds_many([-3.0e8 + 5.1e7 * k + 0.0625 * k for k in range(12)])
        out = ffi.new("myb_rotation*")
        for instant, rotation in zip(t, mars.rotation_many(t)):
            lib.myb_wgccre_rotation(mars.pointer, instant, out)
            self.assertSame(rotation, raw_rotation(out))
        provider = mars.orientation()
        orientations = myb.orientation_many(provider, t)
        o = ffi.new("myb_orientation*")
        for instant, orientation in zip(t, orientations):
            ok(lib.myb_orientation_wgccre(provider.context, instant, o))
            self.assertSame(orientation, myb.Orientation(raw_rotation(o.rotation), o.ut, o.xp,
                                                         o.yp))
        pck = myb.Pck(os.path.join(TEST_DATA, "frames", "graph.bpc"))
        t = [600002748.0 + 160.0 * k + 0.25 * k for k in range(12)]
        rotations, bases = pck.rotation_many(31099, t)
        self.assertSame(list(zip(rotations, bases)),
                        [pck.rotation(31099, myb.time_from_seconds(s)) for s in t])
        # the frame written again by the PCK writer from the set as a provider
        writer = myb.PckWriter()
        error = writer.add_chebyshev(31099, bases[0], "again", pck.orientation(31099), t[0],
                                     t[-1], 2, 12, max_error=True)
        self.assertLess(error, 1e-9)
        again = myb.Pck()
        again.load_memory(writer.write())
        self.assertLess(max(abs(a - b) for a, b in zip(again.rotation(31099, t[3])[0].m,
                                                       rotations[3].m)), 1e-9)
        late = list(t)
        late[6] = 600002748.0 + 1.0e10
        self.assertEqual(self.assertStatus(lib.MYB_E_RANGE, pck.rotation_many, 31099,
                                           late).index, 6)
        # a graph with a fixed frame and a frame turned by a Python provider
        graph = myb.FrameGraph()
        tilt = myb.Rotation((1.0, 0.0, 0.0, 0.0, 0.6, 0.8, 0.0, -0.8, 0.6), (0.0,) * 9)
        graph.add_fixed(2, 1, tilt)
        graph.add(3, 2, lambda seconds: myb.Rotation(
            (math.cos(1e-5 * seconds), math.sin(1e-5 * seconds), 0,
             -math.sin(1e-5 * seconds), math.cos(1e-5 * seconds), 0, 0, 0, 1), None))
        many = graph.rotation_many(1, 3, t)
        self.assertSame(many, [graph.rotation(1, 3, myb.time_from_seconds(s)) for s in t])
        self.assertSame(myb.orientation_many(graph.orientation(1, 3), t)[0].rotation, many[0])
        self.assertSame(myb.orientation_many(myb.orientation_fixed(tilt), t[:2])[1].rotation,
                        tilt)
        self.assertStatus(lib.MYB_E_NOT_FOUND, graph.rotation_many, 1, 99, t)

        def broken(seconds):
            raise ArithmeticError("no frame today")

        graph.add(4, 1, broken)
        with self.assertRaises(ArithmeticError):
            graph.rotation_many(1, 4, t)

