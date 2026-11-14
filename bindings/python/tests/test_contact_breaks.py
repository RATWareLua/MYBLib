# SPDX-License-Identifier: BSD-3-Clause
"""The last parts of the interface through the package against the same C calls through
cffi, bit for bit: the schemes of the breaks of a force model, the timescale of the steps of
IAS15 in an Orbit and in the integrator of a System, the friction that sticks (the force with
its springs, their place in the state of a System), the PN orders by name with the new bits;
failures as myb.Error with the status."""
import myb
from myb import ffi, lib

from common import ELEMENTS, GM, Case, doubles, ok

PLANE = {"shape": "plane", "point": (0.0, 0.0, 0.0), "normal": (0.0, 0.0, 1.0)}
STICK = {"stiffness": 1.0e4, "exponent": 1.5, "damping": 0.1, "friction": 0.5,
         "kind": "stick", "tangential_stiffness": 1.0e3, "tangential_damping": 10.0}
SMOOTH = {"stiffness": 1.0e4, "exponent": 1.5, "damping": 0.1, "friction": 0.5, "slip": 0.01}
LEGS = [((0.3, 0.0, -0.5), 0.1), ((-0.15, 0.26, -0.5), 0.1), ((-0.15, -0.26, -0.5), 0.1)]


def push(z):
    """A force that is continuous but not smooth where z changes sign."""
    return 1.0e-12 * max(z, 0.0)


def raw_contact(law):
    surface = ffi.new("myb_contact_surface*", {"shape": lib.MYB_CONTACT_PLANE,
                                               "normal": [0.0, 0.0, 1.0]})
    made = ffi.new("myb_contact_law*", {"stiffness": law["stiffness"],
                                        "exponent": law["exponent"],
                                        "damping": law["damping"], "friction": law["friction"],
                                        "slip": law.get("slip", 0.0)})
    if law.get("kind") == "stick":
        made.kind = lib.MYB_CONTACT_FRICTION_STICK
        made.tangential_stiffness = law["tangential_stiffness"]
        made.tangential_damping = law["tangential_damping"]
    points = ffi.new("myb_contact_point[]", [{"offset": list(o), "radius": r} for o, r in LEGS])
    return surface, made, points


class Breaks(Case):
    def test_schemes_of_breaks_and_the_timescale(self):
        y0 = myb.state_from_elements(GM, ELEMENTS)
        model = myb.ForceModel()
        model.add_point_mass("earth", GM)
        model.add_term("push", lambda t, x, v: (0.0, 0.0, push(x[2])))
        model.add_break(lambda t, y: y[2])
        self.assertEqual((model.break_count, model.break_scheme),
                         (1, lib.MYB_BREAKS_RETAKE))
        model.set_break_scheme("predict")
        self.assertEqual(model.break_scheme, lib.MYB_BREAKS_PREDICT)
        orbit = myb.Orbit(model, epsilon=1e-11, timescale="jerk")
        orbit.start(0.0, y0, 12000.0)
        final = orbit.state(12000.0)
        # the same through the C functions
        raw = ffi.new("myb_force_model**")
        ok(lib.myb_force_model_create(ffi.NULL, raw))
        term = ffi.new("myb_point_mass_term*", [GM])
        ok(lib.myb_force_model_add(raw[0], b"earth", lib.myb_term_point_mass,
                                   lib.myb_term_point_mass_jacobian, term, ffi.NULL))

        @ffi.callback("myb_status(void*, double, const double[3], const double[3], double[3])")
        def raw_push(context, t, x, v, a):
            a[0:3] = [0.0, 0.0, push(x[2])]
            return lib.MYB_OK

        @ffi.callback("double(void*, double, const double*)")
        def raw_break(context, t, y):
            return y[2]

        ok(lib.myb_force_model_add(raw[0], b"push", raw_push, ffi.NULL, ffi.NULL, ffi.NULL))
        ok(lib.myb_force_model_add_break(raw[0], raw_break, ffi.NULL))
        lib.myb_force_model_set_break_scheme(raw[0], lib.MYB_BREAKS_PREDICT)
        options = ffi.new("myb_orbit_options*", {"epsilon": 1e-11,
                                                 "timescale": lib.MYB_IAS15_TIMESCALE_JERK})
        raw_orbit = ffi.new("myb_orbit**")
        ok(lib.myb_orbit_create(raw[0], options, ffi.NULL, raw_orbit))
        ok(lib.myb_orbit_start(raw_orbit[0], 0.0, ffi.new("double[6]", y0), 12000.0))
        y = ffi.new("double[6]")
        ok(lib.myb_orbit_state(raw_orbit[0], 12000.0, y))
        self.assertSame(final, doubles(y, 6))
        self.assertEqual(orbit.steps, lib.myb_orbit_steps(raw_orbit[0]))
        lib.myb_orbit_destroy(raw_orbit[0])
        lib.myb_force_model_destroy(raw[0])
        # failures
        self.assertStatus(lib.MYB_E_DOMAIN, myb.Orbit, model, timescale=7)
        with self.assertRaises(KeyError):
            model.set_break_scheme("sometimes")
        with self.assertRaises(ValueError):
            model.set_break_scheme(5)

        def broken(t, y):
            raise TimeoutError("no edge today")

        failing = myb.ForceModel()
        failing.add_point_mass("earth", GM)
        failing.add_break(broken)
        lost = myb.Orbit(failing)
        with self.assertRaises(TimeoutError):
            lost.start(0.0, y0, 600.0)
            lost.state(600.0)


class Contact(Case):
    X = (0.01, 0.0, 0.55)
    V = (0.1, 0.02, -0.01)

    def test_the_force_of_a_friction_that_sticks(self):
        deflection = (1e-4, 0.0, 0.0, 0.0, 2e-5, 0.0, 0.0, 0.0, 0.0)
        result, rate = myb.contact_force_stick(PLANE, STICK, LEGS, self.X, self.V, deflection,
                                               omega=(0.0, 0.0, 0.1))
        surface, law, points = raw_contact(STICK)
        made = ffi.new("myb_contact_result*")
        z = ffi.new("double[9]")
        ok(lib.myb_contact_force_stick(surface, law, points, 3, ffi.new("double[3]", self.X),
                                       ffi.new("double[3]", self.V), ffi.NULL,
                                       ffi.new("double[3]", [0.0, 0.0, 0.1]),
                                       ffi.new("double[9]", deflection), z, made))
        self.assertSame((result, rate), (myb.ContactResult(tuple(made.force), tuple(made.torque),
                                                           made.energy, made.gap, made.touching),
                                         doubles(z, 9)))
        self.assertEqual(result.touching, 3)
        smooth = myb.contact_force(PLANE, SMOOTH, LEGS, self.X, self.V)
        surface, law, points = raw_contact(SMOOTH)
        ok(lib.myb_contact_force(surface, law, points, 3, ffi.new("double[3]", self.X),
                                 ffi.new("double[3]", self.V), ffi.NULL, ffi.NULL, made))
        self.assertSame(smooth.force, tuple(made.force))
        # the smooth friction through the function of the one that sticks: the same force
        self.assertSame(myb.contact_force_stick(PLANE, SMOOTH, LEGS, self.X, self.V,
                                                [0.0] * 9)[0], smooth)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.contact_force, PLANE, STICK, LEGS, self.X,
                          self.V)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.contact_force_stick, PLANE, STICK, LEGS, self.X,
                          self.V, [float("nan")] * 9)
        with self.assertRaises(ValueError):
            myb.contact_force_stick(PLANE, STICK, LEGS, self.X, self.V, [0.0] * 6)

    def test_the_springs_in_the_state_of_a_system(self):
        ground = 4.9028e12
        system = myb.System([ground, 0.0])
        system.add_contact(1, 0, PLANE, STICK, LEGS)
        n = system.coordinates
        self.assertEqual(n, 6 + 9)
        deflection = [1e-5 * k for k in range(9)]
        y = [0.0] * (2 * n)
        y[3:6] = (0.0, 0.0, 1.0e6)
        y = system.pack_contact(y, 0, deflection)
        self.assertSame(system.unpack_contact(y, 0), tuple(deflection))
        a = system.acceleration(0.0, y[:n], y[n:])
        # the same through the C functions
        s = ffi.new("myb_system**")
        ok(lib.myb_system_create(2, ffi.new("double[2]", [ground, 0.0]), ffi.NULL, s))
        surface, law, points = raw_contact(STICK)
        contact = ffi.new("myb_system_contact*", {"surface": surface[0], "law": law[0],
                                                  "points": points, "point_count": 3})
        ok(lib.myb_system_add_contact(s[0], 1, 0, contact))
        raw_y = ffi.new("double[30]")
        raw_y[5] = 1.0e6
        ok(lib.myb_system_pack_contact(s[0], 0, ffi.new("double[9]", deflection), raw_y))
        self.assertSame(y, doubles(raw_y, 30))
        out = ffi.new("double[15]")
        ok(lib.myb_system_acceleration(s[0], 0.0, raw_y, raw_y + 15, out))
        self.assertSame(a, doubles(out, 15))
        lib.myb_system_destroy(s[0])
        # the integrator of a system with the timescale of the jerk
        integrator = system.integrator(timescale="jerk")
        integrator.start(0.0, y)
        integrator.integrate(10.0)
        self.assertEqual(integrator.time, 10.0)
        self.assertStatus(lib.MYB_E_DOMAIN, system.integrator, timescale=7)
        smooth = myb.System([ground, 0.0])
        smooth.add_contact(1, 0, PLANE, SMOOTH, LEGS)
        self.assertStatus(lib.MYB_E_DOMAIN, smooth.unpack_contact, [0.0] * 12, 0)
        with self.assertRaises(IndexError):
            system.pack_contact(y, 3, deflection)

    def test_the_timescale_of_the_integrator_of_a_system(self):
        gm = (4.0e5, 1.0e5)
        y0 = [-1000.0, 0.0, 0.0, 4000.0, 0.0, 0.0] + [0.0] * 6
        system = myb.System(gm)
        system.set_radius(0, 60.0)
        system.set_radius(1, 40.0)
        integrator = system.integrator(epsilon=1e-10, timescale="jerk")
        integrator.start(0.0, y0)
        event = system.advance(integrator, 1.0e6)
        s = ffi.new("myb_system**")
        ok(lib.myb_system_create(2, ffi.new("double[2]", gm), ffi.NULL, s))
        ok(lib.myb_system_set_radius(s[0], 0, 60.0))
        ok(lib.myb_system_set_radius(s[0], 1, 40.0))
        ias = ffi.new("myb_ias15**")
        options = ffi.new("myb_ias15_options*", {"epsilon": 1e-10,
                                                 "timescale": lib.MYB_IAS15_TIMESCALE_JERK})
        ok(lib.myb_system_integrator(s[0], options, ffi.NULL, ias))
        ok(lib.myb_ias15_start(ias[0], lib.myb_system_acceleration, s[0], 0.0,
                               ffi.new("double[12]", y0)))
        made = ffi.new("myb_system_event*")
        ok(lib.myb_system_advance(s[0], ias[0], 1.0e6, 0.0, made))
        self.assertSame(tuple(event), (made.kind, made.a, made.b, made.t, made.distance))
        self.assertSame(integrator.state(), doubles(lib.myb_ias15_state(ias[0]), 12))
        self.assertEqual(integrator.statistics.accepted,
                         lib.myb_ias15_statistics(ias[0]).accepted)
        lib.myb_ias15_destroy(ias[0])
        lib.myb_system_destroy(s[0])


class PostNewtonian(Case):
    def test_orders_by_name(self):
        self.assertEqual(myb.pn_orders(["SS3", "SORR5", "SOTAIL", "TAILX"]),
                         lib.MYB_PN_SS3 | lib.MYB_PN_SORR5 | lib.MYB_PN_SOTAIL | lib.MYB_PN_TAILX)
        self.assertEqual(myb.pn_orders(["2", lib.MYB_PN_SO, "4rr"]),
                         lib.MYB_PN_2 | lib.MYB_PN_SO | lib.MYB_PN_4RR)
        with self.assertRaises(ValueError):
            myb.pn_orders(["5"])
        gm = (1.0e20, 6.0e19)
        system = myb.System(gm)
        system.add_pn_pair(0, 1, ["2", "SO", "SS", "SS3"], c=299792458.0)
        n = system.coordinates
        y = [0.0] * (2 * n)
        y[3] = 1.0e9
        y[n + 4] = 1.0e5
        y = system.pack_spins(y, 0, (0.1, 0.0, 0.5), (0.0, 0.2, 0.1))
        a = system.acceleration(0.0, y[:n], y[n:])
        s = ffi.new("myb_system**")
        ok(lib.myb_system_create(2, ffi.new("double[2]", gm), ffi.NULL, s))
        pair = ffi.new("myb_system_pn_pair*", {"orders": lib.MYB_PN_2 | lib.MYB_PN_SO |
                                               lib.MYB_PN_SS | lib.MYB_PN_SS3,
                                               "c": 299792458.0})
        ok(lib.myb_system_add_pn_pair(s[0], 0, 1, pair))
        out = ffi.new("double[]", n)
        given = ffi.new("double[]", y)
        ok(lib.myb_system_acceleration(s[0], 0.0, given, given + n, out))
        self.assertSame(a, doubles(out, n))
        lib.myb_system_destroy(s[0])
        self.assertStatus(lib.MYB_E_DOMAIN, system.add_pn_pair, 0, 1, ["4"], c=299792458.0)
