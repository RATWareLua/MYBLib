# SPDX-License-Identifier: BSD-3-Clause
"""The Python package against the C interface it wraps.

  - the low orbit of examples/low_orbit.c built from the classes gives the nodes and the
    final state of the same calls made through myb.lib — the C functions themselves —
    exactly;
  - errors are exceptions that carry the status of the library and the message of a parser;
  - terms, providers and events as Python callables; an exception raised inside one comes
    out of the call of the library that was going on;
  - objects keep alive what they point to, and can be closed;
  - kernels and time scales, an SPK set and a body of it in an environment, messages
    written and read in both notations, an optional library of ext/;
  - the declarations of the package are those the headers give now.

Run with the library named by MYB_LIBRARY; MYB_DATA_DIR and MYB_TEST_DATA_DIR point at data/
and tests/data/ of the repository (found from this file if they are not set).
"""
import gc
import math
import os
import shutil
import subprocess
import sys
import unittest

import myb
from myb import ffi, lib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
DATA = os.environ.get("MYB_DATA_DIR", os.path.join(ROOT, "data"))
TEST_DATA = os.environ.get("MYB_TEST_DATA_DIR", os.path.join(ROOT, "tests", "data"))

GM = 3.986004418e14
RADIUS = 6378137.0
J2 = 1.08263e-3
DAY = 86400.0
C20 = -J2 / math.sqrt(5.0)
ELEMENTS = myb.Elements(p=(RADIUS + 400000.0) * (1.0 - 0.001 ** 2), e=0.001,
                        i=51.6 * math.pi / 180.0, raan=0.5, argp=1.0, nu=0.0)


def low_orbit_raw():
    """The example call for call through the C functions."""
    c = ffi.new("double[6]", [1.0, 0.0, 0.0, C20, 0.0, 0.0])
    s = ffi.new("double[6]")
    field = ffi.new("myb_gravity_field**")
    work = ffi.new("myb_gravity_work**")
    assert lib.myb_gravity_field_create(GM, RADIUS, 2, c, s, ffi.NULL, field) == lib.MYB_OK
    assert lib.myb_gravity_work_create(2, ffi.NULL, work) == lib.MYB_OK
    options = ffi.new("myb_env_options*")
    env = ffi.new("myb_env**")
    assert lib.myb_env_create(options, ffi.NULL, env) == lib.MYB_OK
    body = ffi.new("myb_env_body*", {"id": 1, "gm": GM})
    planet = ffi.new("size_t*")
    assert lib.myb_env_add_body(env[0], body, planet) == lib.MYB_OK
    gravity = ffi.new("myb_field_term*")
    gravity.env = env[0]
    gravity.body = planet[0]
    gravity.field = field[0]
    gravity.work = work[0]
    gravity.n_max = -1
    air = ffi.new("myb_exponential_atmosphere*")
    air.ellipsoid.a = RADIUS
    air.count = 1
    air.base[0] = 400000.0
    air.density[0] = 3.0e-12
    air.scale[0] = 58000.0
    thing = ffi.new("myb_object*", {"mass": 500.0, "area": 4.0, "cr": 0.0, "cd": 2.2})
    drag = ffi.new("myb_drag_term*")
    drag.env = env[0]
    drag.body = planet[0]
    drag.density = lib.myb_density_exponential
    drag.density_context = air
    drag.object = thing
    model = ffi.new("myb_force_model**")
    assert lib.myb_force_model_create(ffi.NULL, model) == lib.MYB_OK
    assert lib.myb_force_model_add(model[0], b"gravity", lib.myb_term_field,
                                   lib.myb_term_field_jacobian, gravity, ffi.NULL) == lib.MYB_OK
    assert lib.myb_force_model_add(model[0], b"drag", lib.myb_term_drag,
                                   lib.myb_term_drag_jacobian, drag, ffi.NULL) == lib.MYB_OK
    elements = ffi.new("myb_elements*", ELEMENTS._asdict())
    pv = ffi.new("myb_pv*")
    assert lib.myb_pv_from_elements(GM, elements, pv) == lib.MYB_OK
    y0 = ffi.new("double[6]", list(pv.r) + list(pv.v))
    orbit = ffi.new("myb_orbit**")
    assert lib.myb_orbit_create(model[0], ffi.NULL, ffi.NULL, orbit) == lib.MYB_OK
    assert lib.myb_orbit_start(orbit[0], 0.0, y0, DAY) == lib.MYB_OK

    @ffi.callback("double(void*, double, const double*)")
    def equator(context, t, y):
        return y[2]

    nodes = []
    found = ffi.new("int*")
    t_node = ffi.new("double*")
    for _ in range(32):
        assert lib.myb_orbit_event(orbit[0], equator, ffi.NULL, DAY, 1, found, t_node,
                                   ffi.NULL) == lib.MYB_OK
        if not found[0]:
            break
        nodes.append(t_node[0])
    final = ffi.new("double[6]")
    assert lib.myb_orbit_state(orbit[0], DAY, final) == lib.MYB_OK
    steps = lib.myb_orbit_steps(orbit[0])
    lib.myb_orbit_destroy(orbit[0])
    lib.myb_force_model_destroy(model[0])
    lib.myb_env_destroy(env[0])
    lib.myb_gravity_work_destroy(work[0])
    lib.myb_gravity_field_destroy(field[0])
    return nodes, tuple(final), steps


def low_orbit_model():
    """The model of the example from the classes; what it is built from is dropped."""
    field = myb.GravityField(GM, RADIUS, 2, [1.0, 0.0, 0.0, C20, 0.0, 0.0], [0.0] * 6)
    env = myb.Environment()
    planet = env.add_body(1, GM)
    model = myb.ForceModel()
    model.add_field("gravity", env, planet, field)
    model.add_drag("drag", env, planet, mass=500.0, area=4.0, cd=2.2,
                   layers=[(400000.0, 3.0e-12, 58000.0)], equatorial_radius=RADIUS)
    return model


class LowOrbit(unittest.TestCase):
    def test_classes_give_the_bits_of_the_c_calls(self):
        nodes_raw, final_raw, steps_raw = low_orbit_raw()
        model = low_orbit_model()
        gc.collect()  # the model keeps the field and the environment
        orbit = myb.Orbit(model)
        del model
        gc.collect()  # and the orbit keeps the model
        orbit.start(0.0, myb.state_from_elements(GM, ELEMENTS), DAY)
        nodes = []
        for _ in range(32):
            node = orbit.next_event(lambda t, y: y[2], DAY, direction=1)
            if node is None:
                break
            nodes.append(node[0])
        final = orbit.state(DAY)
        self.assertGreaterEqual(len(nodes), 14)
        self.assertEqual(nodes, nodes_raw)
        self.assertEqual(final, final_raw)
        self.assertEqual(orbit.steps, steps_raw)
        self.assertEqual(orbit.time, DAY)
        # the node drifts as the first-order theory says
        a0 = RADIUS + 400000.0
        theory = -1.5 * math.sqrt(GM / a0 ** 3) * J2 * (RADIUS / ELEMENTS.p) ** 2 * math.cos(ELEMENTS.i)
        after = myb.elements_from_state(GM, final)
        self.assertAlmostEqual((after.raan - ELEMENTS.raan) / DAY / theory, 1.0, delta=0.02)

    def test_states_at_many_times_and_partials(self):
        model = myb.ForceModel()
        model.add_point_mass("planet", GM)
        orbit = myb.Orbit(model, variational=True)
        y0 = myb.state_from_elements(GM, ELEMENTS)
        orbit.start(0.0, y0, 6000.0)
        many = orbit.states([1000.0, 2000.0, 3000.0])
        self.assertEqual(len(many), 3)
        for t, y in zip((1000.0, 2000.0, 3000.0), many):
            kepler = myb.kepler_propagate(GM, y0, t)
            self.assertLess(max(abs(a - b) for a, b in zip(y[:3], kepler[:3])), 1e-6)
        y, phi = orbit.partials(4000.0)
        self.assertEqual(y, orbit.state(4000.0))
        self.assertEqual((len(phi), len(phi[0])), (6, 6))
        self.assertNotEqual(phi[0][0], 0.0)
        orbit.kick((0.0, 1.0, 0.0))
        other = myb.Orbit(model, integrator="dop853", rtol=1e-11, atol=1e-6)
        other.start(0.0, y0, 3000.0)
        self.assertLess(max(abs(a - b) for a, b in zip(other.state(3000.0)[:3], many[2][:3])), 1e-2)
        orbit.close()
        with self.assertRaises(ValueError):
            orbit.state(4500.0)


class Errors(unittest.TestCase):
    def test_a_status_becomes_an_exception(self):
        with self.assertRaises(myb.Error) as caught:
            myb.state_from_elements(GM, myb.Elements(-1.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        self.assertEqual(caught.exception.status, lib.MYB_E_DOMAIN)
        self.assertEqual(caught.exception.name, "MYB_E_DOMAIN")
        with self.assertRaises(ValueError):
            myb.kepler_propagate(GM, (1.0, 2.0, 3.0), 10.0)

    def test_a_parser_names_the_line(self):
        with self.assertRaises(myb.Error) as caught:
            myb.read_oem("CCSDS_OEM_VERS = 2.0\nMETA_START\n")
        self.assertEqual(caught.exception.status, lib.MYB_E_FORMAT)
        self.assertIn("line 2", str(caught.exception))
        with self.assertRaises(myb.Error) as caught:
            myb.KernelPool(os.path.join(DATA, "no-such-kernel.tk"))
        self.assertEqual(caught.exception.status, lib.MYB_E_IO)


class Callables(unittest.TestCase):
    X = (7.0e6, 1.0e5, -2.0e5)
    V = (0.0, 7.5e3, 1.0e3)

    def test_a_term_of_the_caller(self):
        model = myb.ForceModel()
        model.add_point_mass("planet", GM)
        push = (1.0e-6, -2.0e-6, 3.0e-6)
        index = model.add_term("push", lambda t, x, v: push)
        self.assertEqual(model.term_names, ["planet", "push"])
        both = model.acceleration(0.0, self.X, self.V)
        model.set_enabled(index, False)
        self.assertFalse(model.enabled(index))
        planet = model.acceleration(0.0, self.X, self.V)
        for a, b, p in zip(both, planet, push):
            self.assertAlmostEqual(a - b, p, delta=1e-14)

    def test_an_exception_inside_a_term_comes_out(self):
        class Refused(Exception):
            pass

        def refuses(t, x, v):
            raise Refused("not today")

        model = myb.ForceModel()
        model.add_term("refuses", refuses)
        with self.assertRaises(Refused):
            model.acceleration(0.0, self.X, self.V)
        # and out of a propagation, which stops there
        model.add_point_mass("planet", GM)
        orbit = myb.Orbit(model)
        with self.assertRaises(Refused):
            orbit.start(0.0, myb.state_from_elements(GM, ELEMENTS), 600.0)
            orbit.state(600.0)
        # a term that reports a status of the library
        failing = myb.ForceModel()

        def out_of_range(t, x, v):
            raise myb.Error(lib.MYB_E_RANGE)

        failing.add_term("range", out_of_range)
        with self.assertRaises(myb.Error) as caught:
            failing.acceleration(0.0, self.X, self.V)
        self.assertEqual(caught.exception.status, lib.MYB_E_RANGE)
        # a wrong shape of the answer
        wrong = myb.ForceModel()
        wrong.add_term("two", lambda t, x, v: (1.0, 2.0))
        with self.assertRaises(ValueError):
            wrong.acceleration(0.0, self.X, self.V)

    def test_an_event_that_raises(self):
        model = myb.ForceModel()
        model.add_point_mass("planet", GM)
        orbit = myb.Orbit(model)
        orbit.start(0.0, myb.state_from_elements(GM, ELEMENTS), 6000.0)

        def broken(t, y):
            raise KeyError("no such thing")

        with self.assertRaises(KeyError):
            orbit.next_event(broken, 6000.0)

    def test_a_provider_in_python(self):
        rate = 2.0e-7
        distance = 1.5e11

        def sun(seconds):
            angle = rate * seconds
            return (distance * math.cos(angle), distance * math.sin(angle), 0.0,
                    -distance * rate * math.sin(angle), distance * rate * math.cos(angle), 0.0)

        env = myb.Environment()
        env.add_body(399, GM)
        index = env.add_body(10, 1.32712440018e20, state=sun)
        self.assertEqual(env.body_count, 2)
        self.assertEqual(env.find(10), index)
        self.assertEqual(env.state(1000.0, index), sun(1000.0))
        model = myb.ForceModel()
        model.add_point_mass("planet", GM)
        model.add_third_body("sun", env, index)
        del env
        gc.collect()
        with_sun = model.acceleration(0.0, self.X, self.V)
        model.set_enabled(1, False)
        alone = model.acceleration(0.0, self.X, self.V)
        tide = math.sqrt(sum((a - b) ** 2 for a, b in zip(with_sun, alone)))
        self.assertTrue(1e-7 < tide < 1e-5)  # 2 GM r / d^3 at 7000 km is 5e-7 m/s^2


class TimeAndEphemeris(unittest.TestCase):
    def test_time_scales(self):
        kernels = [os.path.join(DATA, "time", name)
                   for name in ("time-scales.tk", "tdb-tt-fb1990.tk", "utc-tai.tk")]
        with myb.KernelPool(*kernels) as pool:
            system = myb.TimeSystem(pool)
        label = (2026, 11, 11, 12, 30, 15.25)
        instant = system.from_utc(label, "TT")
        back = system.to_utc(instant, "TT")
        self.assertEqual(back[:5], label[:5])
        self.assertAlmostEqual(back[5], label[5], delta=1e-9)
        tai = system.convert(instant, "TT", "TAI")
        self.assertAlmostEqual(myb.seconds(instant) - myb.seconds(tai), 32.184, delta=1e-6)
        with self.assertRaises(myb.Error) as caught:
            system.convert(instant, "TT", "NO SUCH SCALE")
        self.assertEqual(caught.exception.status, lib.MYB_E_NOT_FOUND)

    def test_an_spk_set_and_a_body_of_it(self):
        eph = myb.Ephemeris(os.path.join(TEST_DATA, "spk", "t02_03.bsp"))
        first = eph.segments()[0]
        middle = myb.time_from_seconds(0.5 * (first["start"] + first["stop"]))
        there = eph.state(first["target"], first["center"], middle)
        raw = ffi.new("double[6]")
        self.assertEqual(lib.myb_spk_state(eph.pointer, first["target"], first["center"], middle,
                                           raw, ffi.NULL), lib.MYB_OK)
        self.assertEqual(there, tuple(raw))
        with self.assertRaises(myb.Error) as caught:
            eph.state(first["target"], first["center"],
                      myb.time_from_seconds(first["stop"] + 1.0e9))
        self.assertEqual(caught.exception.status, lib.MYB_E_RANGE)
        env = myb.Environment(epoch=middle)
        env.add_body(first["center"], GM)
        index = env.add_body(first["target"], 1.0,
                             state=eph.body(first["target"], first["center"]))
        del eph
        gc.collect()  # the environment keeps the ephemeris
        self.assertEqual(env.state(0.0, index), there)


class Messages(unittest.TestCase):
    def test_an_ephemeris_message_in_both_notations(self):
        epochs = [(2026, 3, 1, 12, 0, 0.0), (2026, 3, 1, 12, 1, 0.5)]
        states = [(7.0e6, 1.0, -2.5, 0.25, 7.5e3, -1.0e3),
                  (7.0e6 + 15.0, 4.5e5, -6.0e4, -480.0, 7.4e3, -990.0)]
        segment = {"object_name": "A & B", "object_id": "2026-001A", "center_name": "EARTH",
                   "ref_frame": "EME2000", "time_system": "UTC", "epochs": epochs,
                   "states": states, "interpolation": "HERMITE", "interpolation_degree": 7}
        for notation in ("kvn", "xml"):
            text = myb.write_oem("myb", (2026, 3, 2, 0, 0, 0.0), [segment], notation)
            self.assertEqual(text.startswith("<"), notation == "xml")
            read = myb.read_oem(text)
            self.assertEqual(len(read), 1)
            self.assertEqual(read[0]["meta"]["OBJECT_NAME"], "A & B")
            self.assertEqual(read[0]["meta"]["INTERPOLATION_DEGREE"], "7")
            self.assertEqual([label for label, _ in read[0]["states"]], epochs)
            self.assertEqual([state for _, state in read[0]["states"]], states)
        with self.assertRaises(myb.Error):
            myb.write_oem("myb", (2026, 3, 2, 0, 0, 0.0), [])


class Everything(unittest.TestCase):
    def test_the_whole_interface_is_there(self):
        # a function no class wraps: the shadow of two bodies at once
        source = ffi.new("double[3]", [-1.5e11, 0.0, 0.0])
        bodies = ffi.new("double[6]", [-4.0e8, 7.3e6, 0.0, -1.0e8, -1.9e6, 0.0])
        radii = ffi.new("double[2]", [6.4e6, 1.7e6])
        limb = ffi.new("double[2]", [0.4, 0.2])
        seen = lib.myb_shadow_fraction_many(source, 6.957e8, bodies, radii, 2, limb)
        self.assertTrue(0.0 < seen < 1.0)
        self.assertEqual(lib.myb_version_number(), lib.MYB_VERSION_NUMBER)
        self.assertEqual(myb.__version__, "%d.%d.%d.%d" % myb.version)
        self.assertTrue(os.path.exists(myb.library_path) or os.path.dirname(myb.library_path) == "")

    def test_an_optional_library(self):
        try:
            sgp4 = myb.extension("sgp4")
        except ImportError:
            self.skipTest("libmyb_sgp4 is not next to the library")
        line1 = b"1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
        line2 = b"2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"
        tle = ffi.new("myb_tle*")
        diag = ffi.new("myb_diag*")
        self.assertEqual(sgp4.myb_tle_parse(line1, line2, 0, tle, diag), lib.MYB_OK)
        self.assertEqual(tle.catalog_number, 5)
        self.assertAlmostEqual(tle.inclination, 34.2682, delta=1e-12)
        with self.assertRaises(ValueError):
            myb.extension("no-such-module")

    def test_the_declarations_are_those_of_the_headers(self):
        clang = os.environ.get("MYB_CLANG") or shutil.which("clang")
        if not clang or not os.path.exists(os.path.join(ROOT, "include", "myb", "myb.h")):
            self.skipTest("no clang or no headers to regenerate the declarations from")
        tool = os.path.join(ROOT, "bindings", "python", "tools", "make_cdef.py")
        done = subprocess.run([sys.executable, tool, "--clang", clang, "--check", ROOT],
                              capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)


if __name__ == "__main__":
    unittest.main()
