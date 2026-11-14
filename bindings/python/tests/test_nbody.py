# SPDX-License-Identifier: BSD-3-Clause
"""N bodies through the classes against the same C calls through cffi, bit for bit: the
accelerations and invariants of an NBody, the fast multipole approximation, a Wisdom-Holman
integration, a System of two spheres that fall, touch, rebound and merge; the levels of a
System with a Python action; failures as myb.Error with the status."""
import gc

import myb
from myb import ffi, lib

from common import Case, doubles, ok, rows

SUN = 1.32712440018e20
PLANETS = [SUN, 1.26686534e17, 3.7931187e16]
X = [(0.0, 0.0, 0.0), (7.78e11, 0.0, 1.0e9), (0.0, 1.43e12, -2.0e9)]
V = [(0.0, 0.0, 0.0), (0.0, 13.07e3, 0.0), (-9.68e3, 0.0, 100.0)]


def flat(values):
    return [v for row in values for v in row]


class NBodies(Case):
    def test_accelerations_and_invariants(self):
        nbody = myb.NBody(PLANETS, test=1)
        self.assertEqual((nbody.count, nbody.massive), (4, 3))
        self.assertEqual(nbody.gm(3), 0.0)
        x = X + [(1.0e11, 1.0e11, 0.0)]
        v = V + [(0.0, 3.0e4, 0.0)]
        a = nbody.acceleration(0.0, x)
        raw = ffi.new("myb_nbody**")
        ok(lib.myb_nbody_create(3, 1, ffi.new("double[]", PLANETS), ffi.NULL, raw))
        out = ffi.new("double[12]")
        ok(lib.myb_nbody_acceleration(raw[0], 0.0, ffi.new("double[]", flat(x)), out))
        self.assertSame(a, rows(out, 4, 3))
        invariants = nbody.invariants(x, v)
        made = ffi.new("myb_nbody_invariants*")
        ok(lib.myb_nbody_invariants_of(raw[0], ffi.new("double[]", flat(x)),
                                       ffi.new("double[]", flat(v)), made))
        self.assertSame(invariants.energy, (made.energy.hi, made.energy.lo))
        self.assertSame(invariants.angular_momentum, tuple(made.angular_momentum))
        # the graph, the approximation, the threads
        nbody.set_interaction(1, 2, "none")
        ok(lib.myb_nbody_set_interaction(raw[0], 1, 2, lib.MYB_INTERACTION_NONE))
        nbody.set_approximation(8, 0.5)
        ok(lib.myb_nbody_set_approximation(raw[0], 8, 0.5))
        self.assertEqual(nbody.approximation, (8, 0.5))
        with myb.Pool(2) as pool:
            nbody.set_pool(pool)
            self.assertEqual(pool.threads, 2)
            ok(lib.myb_nbody_acceleration(raw[0], 0.0, ffi.new("double[]", flat(x)), out))
            self.assertSame(nbody.acceleration(0.0, x), rows(out, 4, 3))
            nbody.set_pool(None)
        nbody.set_approximation(0)
        self.assertEqual(nbody.approximation, (0, 0.0))
        # the 1PN terms: the accelerations need the velocities
        nbody.set_relativity(299792458.0)
        ok(lib.myb_nbody_set_approximation(raw[0], 0, 0.0))
        ok(lib.myb_nbody_set_relativity(raw[0], ffi.new("myb_ppn*", [299792458.0, 1.0, 1.0])))
        ok(lib.myb_nbody_acceleration_xv(raw[0], 0.0, ffi.new("double[]", flat(x)),
                                         ffi.new("double[]", flat(v)), out))
        self.assertSame(nbody.acceleration_xv(0.0, x, v), rows(out, 4, 3))
        pn = nbody.pn_invariants(x, v)
        self.assertEqual(len(pn.center), 3)
        lib.myb_nbody_destroy(raw[0])
        self.assertStatus(lib.MYB_E_DOMAIN, nbody.acceleration, 0.0, x)
        self.assertStatus(lib.MYB_E_DOMAIN, nbody.set_approximation, 30, 0.5)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.NBody, [-1.0])
        coincident = myb.NBody([1.0, 1.0])
        self.assertStatus(lib.MYB_E_DOMAIN, coincident.acceleration, 0.0, [(1, 2, 3), (1, 2, 3)])
        with self.assertRaises(ValueError):
            coincident.acceleration(0.0, [(1, 2, 3)])

    def test_wisdom_holman(self):
        nbody = myb.NBody(PLANETS)
        wh = myb.WisdomHolman(nbody, 86400.0 * 10, scheme="dkd", corrector=11)
        del nbody
        gc.collect()  # the integrator keeps the system
        wh.start(0.0, X, V)
        wh.integrate(86400.0 * 3650)
        x, v = wh.state()
        dense = wh.dense(wh.time + 86400.0)
        raw_nbody = ffi.new("myb_nbody**")
        ok(lib.myb_nbody_create(3, 0, ffi.new("double[]", PLANETS), ffi.NULL, raw_nbody))
        options = ffi.new("myb_wh_options*", {"h": 86400.0 * 10, "scheme": lib.MYB_WH_DKD,
                                              "corrector": 11})
        raw = ffi.new("myb_wh**")
        ok(lib.myb_wh_create(raw_nbody[0], options, ffi.NULL, raw))
        ok(lib.myb_wh_start(raw[0], 0.0, ffi.new("double[]", flat(X)),
                            ffi.new("double[]", flat(V))))
        ok(lib.myb_wh_integrate(raw[0], 86400.0 * 3650))
        rx = ffi.new("double[9]")
        rv = ffi.new("double[9]")
        ok(lib.myb_wh_state(raw[0], rx, rv))
        self.assertSame((x, v), (rows(rx, 3, 3), rows(rv, 3, 3)))
        ok(lib.myb_wh_dense(raw[0], lib.myb_wh_time(raw[0]) + 86400.0, rx, rv))
        self.assertSame(dense, (rows(rx, 3, 3), rows(rv, 3, 3)))
        self.assertEqual(wh.statistics.steps, lib.myb_wh_statistics(raw[0]).steps)
        lib.myb_wh_destroy(raw[0])
        lib.myb_nbody_destroy(raw_nbody[0])
        self.assertStatus(lib.MYB_E_RANGE, wh.dense, wh.time + 1.0e9)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.WisdomHolman, myb.NBody(PLANETS), 0.0)


class Systems(Case):
    GM = (4.0e5, 1.0e5)
    D = 5000.0

    def start(self):
        total = sum(self.GM)
        return [-self.GM[1] / total * self.D, 0.0, 0.0, self.GM[0] / total * self.D, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    def test_spheres_that_touch_rebound_and_merge(self):
        system = myb.System(self.GM)
        system.set_radius(0, 60.0)
        system.set_radius(1, 40.0)
        integrator = system.integrator(epsilon=1e-10)
        integrator.start(0.0, self.start())
        events = [system.advance(integrator, 1.0e6)]
        speed = integrator.state()
        system.bounce(integrator, 0, 1, 1.0)
        events.append(system.advance(integrator, 2.0 * events[0].t))
        events.append(system.advance(integrator, 1.0e6))
        merger = system.merge(integrator, 0, 1)
        final = integrator.state()
        # the same through the C functions
        s = ffi.new("myb_system**")
        ok(lib.myb_system_create(2, ffi.new("double[]", self.GM), ffi.NULL, s))
        ok(lib.myb_system_set_radius(s[0], 0, 60.0))
        ok(lib.myb_system_set_radius(s[0], 1, 40.0))
        ias = ffi.new("myb_ias15**")
        ok(lib.myb_system_integrator(s[0], ffi.new("myb_ias15_options*", {"epsilon": 1e-10}),
                                     ffi.NULL, ias))
        ok(lib.myb_ias15_start(ias[0], lib.myb_system_acceleration, s[0], 0.0,
                               ffi.new("double[]", self.start())))
        event = ffi.new("myb_system_event*")
        raw_events = []
        ok(lib.myb_system_advance(s[0], ias[0], 1.0e6, 0.0, event))
        raw_events.append((event.kind, event.a, event.b, event.t, event.distance))
        raw_speed = tuple(ffi.unpack(lib.myb_ias15_state(ias[0]), 12))
        ok(lib.myb_system_bounce(s[0], ias[0], 0, 1, 1.0))
        ok(lib.myb_system_advance(s[0], ias[0], 2.0 * raw_events[0][3], 0.0, event))
        raw_events.append((event.kind, event.a, event.b, event.t, event.distance))
        ok(lib.myb_system_advance(s[0], ias[0], 1.0e6, 0.0, event))
        raw_events.append((event.kind, event.a, event.b, event.t, event.distance))
        made = ffi.new("myb_system_merger*")
        ok(lib.myb_system_merge(s[0], ias[0], 0, 1, made))
        self.assertSame([tuple(e) for e in events], raw_events)
        self.assertSame(speed, raw_speed)
        self.assertSame(final, tuple(ffi.unpack(lib.myb_ias15_state(ias[0]), 12)))
        self.assertSame(merger, myb.SystemMerger(made.energy, tuple(made.angular_momentum),
                                                 tuple(made.spin)))
        self.assertEqual([e.kind for e in events], [lib.MYB_SYSTEM_CONTACT,
                                                    lib.MYB_SYSTEM_ARRIVED,
                                                    lib.MYB_SYSTEM_CONTACT])
        self.assertTrue(system.present(0) and not system.present(1))
        self.assertEqual(system.radius(0), lib.myb_system_radius(s[0], 0))
        self.assertEqual(system.nbody.gm(1), 0.0)
        lib.myb_ias15_destroy(ias[0])
        lib.myb_system_destroy(s[0])
        self.assertStatus(lib.MYB_E_DOMAIN, system.merge, integrator, 0, 1)  # b is absent
        self.assertStatus(lib.MYB_E_DOMAIN, system.set_radius, 7, 1.0)
        self.assertStatus(lib.MYB_E_DOMAIN, system.bounce, integrator, 0, 1, 2.0)

    def test_levels_actions_and_partials(self):
        system = myb.System([GM_EARTH, 4.9028e12])
        system.set_rotation(1, (0.39, 0.391, 0.392))
        n = system.coordinates
        self.assertEqual(n, 6 + 9)
        y = [0.0] * (2 * n)
        y[3:6] = (3.844e8, 0.0, 0.0)
        y[n + 4] = 1.02e3
        y = system.pack_rotation(y, 1, (1.0, 0.0, 0.0, 0.0), (0.0, 0.0, 2.66e-6))
        q, omega, core = system.unpack_rotation(y, 1)
        self.assertSame(q, (1.0, 0.0, 0.0, 0.0))
        self.assertSame(omega, (0.0, 0.0, 2.66e-6))
        pushes = []

        def push(t, view):
            pushes.append(view)
            return (0.0, 1.0e-9, 0.0), (0.0, 0.0, 1.0e-12)

        system.add_action(1, push)
        a = system.acceleration(0.0, y[:n], y[n:])
        self.assertEqual(pushes[-1].body, 1)
        self.assertEqual(pushes[-1].x[1], (3.844e8, 0.0, 0.0))
        self.assertEqual(len(pushes[-1].frame), 9)
        self.assertEqual(len(a), n)
        inertia = system.inertia(1, 0.0, y)
        self.assertSame(inertia[0], 0.39)
        integrator = system.integrator()
        integrator.start(0.0, y)
        integrator.integrate(3600.0)
        self.assertEqual(integrator.time, 3600.0)
        y1, phi, sensitivity = system.partials(0.0, y, 600.0, bodies=[0])
        self.assertEqual((len(phi), len(phi[0]), len(sensitivity[0])), (2 * n, 2 * n, 1))
        # an exception in the action comes out of the call that was evaluating
        failing = myb.System([GM_EARTH, 1.0])

        def broken(t, view):
            raise LookupError("no thrust")

        failing.add_action(1, broken)
        with self.assertRaises(LookupError):
            failing.acceleration(0.0, [0, 0, 0, 7.0e6, 0, 0], [0.0] * 6)
        # a guided body and the point masses with a provider of states
        guided = myb.System([GM_EARTH, 0.0])
        index = guided.add_guided(4.9e12, lambda seconds: (3.844e8, 0.0, 0.0, 0.0, 1.0e3, 0.0))
        self.assertEqual((index, guided.guided), (0, 1))
        alone = myb.System([GM_EARTH, 0.0])
        x = [0.0, 0.0, 0.0, 7.0e6, 0.0, 0.0]
        pulled = guided.acceleration(0.0, x, [0.0] * 6)
        self.assertNotEqual(pulled[3], alone.acceleration(0.0, x, [0.0] * 6)[3])
        self.assertStatus(lib.MYB_E_DOMAIN, system.unpack_rotation, y, 0)
        self.assertStatus(lib.MYB_E_DOMAIN, system.add_figures, 0, 0)
        nbody = system.nbody
        system.close()
        with self.assertRaises(ValueError):
            nbody.count

    def test_the_models_of_bodies_against_the_c_calls(self):
        identity = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
        c = [1.0, 0.0, 0.0, -4.8e-4, 0.0, 0.0]
        field = myb.GravityField(GM_EARTH, 6378137.0, 2, c, [0.0] * 6)
        system = myb.System([GM_EARTH, 4.9e12, 1.0e3])
        system.set_field(0, field, 2)
        system.set_orientation(0, myb.orientation_fixed(identity))
        system.add_field(0, 1)
        system.set_tide(0, 6378137.0, (0.3, 0.3, 0.3), (0.0, 600.0, 600.0), (0.0, 0.0, 0.0),
                        7.29e-5)
        system.add_tide(0, 1, 1)
        system.nbody.set_relativity(299792458.0)
        system.add_frame_dragging(0, 2, (0.0, 0.0, 1.0e10))
        system.set_mass_rate(0, -1e-15, 0.0)
        n = system.coordinates
        x = [0.0, 0.0, 0.0, 3.844e8, 1.0e6, 2.0e6, 7.0e6, 0.0, 0.0]
        v = [0.0, 0.0, 0.0, 0.0, 1.02e3, 0.0, 0.0, 7.5e3, 0.0]
        a = system.acceleration(10.0, x, v)
        # the same through the C functions
        s = ffi.new("myb_system**")
        ok(lib.myb_system_create(3, ffi.new("double[]", [GM_EARTH, 4.9e12, 1.0e3]), ffi.NULL, s))
        ok(lib.myb_system_set_field(s[0], 0, field.pointer, 2))
        rotation = ffi.new("myb_rotation*", {"m": identity})
        ok(lib.myb_system_set_orientation(s[0], 0, lib.myb_orientation_fixed, rotation))
        ok(lib.myb_system_add_field(s[0], 0, 1))
        tide = ffi.new("myb_system_tide*", {"radius": 6378137.0, "k2": [0.3, 0.3, 0.3],
                                            "orbital_lag": [0.0, 600.0, 600.0],
                                            "rotation_rate": 7.29e-5})
        ok(lib.myb_system_set_tide(s[0], 0, tide))
        ok(lib.myb_system_add_tide(s[0], 0, 1, 1))
        ok(lib.myb_nbody_set_relativity(lib.myb_system_nbody(s[0]),
                                        ffi.new("myb_ppn*", [299792458.0, 1.0, 1.0])))
        ok(lib.myb_system_add_frame_dragging(s[0], 0, 2, ffi.new("double[3]", [0.0, 0.0, 1.0e10])))
        ok(lib.myb_system_set_mass_rate(s[0], 0, -1e-15, 0.0))
        out = ffi.new("double[9]")
        ok(lib.myb_system_acceleration(s[0], 10.0, ffi.new("double[9]", x), ffi.new("double[9]", v),
                                       out))
        self.assertSame(a, doubles(out, n))
        lib.myb_system_destroy(s[0])
        # a lit mesh and a contact on a body with a frame, two spinning bodies beyond 1PN
        shape = myb.Polyhedron([(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1),
                                (0, 0, -1)],
                               [(0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4), (2, 0, 5), (1, 2, 5),
                                (3, 1, 5), (0, 3, 5)])
        self.assertEqual(shape.volume()[0], 4.0 / 3.0)
        lander = myb.System([GM_EARTH, 1.0])
        lander.set_orientation(1, lambda seconds: myb.Rotation(identity, None))
        lander.add_light_action(1, 0, shape, 1361.0, 1.5e11, 299792458.0, 1.0, specular=0.1)
        lander.add_contact(1, 0, {"shape": "ellipsoid", "radii": (6.4e6,) * 3},
                           {"stiffness": 1e3, "exponent": 1.5, "damping": 0.1},
                           [((0.0, 0.0, -1.0), 0.1)])
        self.assertEqual(len(lander.acceleration(0.0, [0, 0, 0, 7.0e6, 0, 0], [0.0] * 6)), 6)
        self.assertStatus(lib.MYB_E_DOMAIN, lander.add_contact, 1, 1, {"shape": "plane"}, {},
                          [((0.0, 0.0, 0.0), 0.1)])
        pair = myb.System([1.0e20, 1.0e20])
        pair.add_pn_pair(0, 1, [lib.MYB_PN_2, lib.MYB_PN_SO], c=299792458.0)
        m = pair.coordinates
        y = pair.pack_spins([0.0] * (2 * m), 0, (0.0, 0.0, 0.5), (0.0, 0.0, 0.1))
        self.assertEqual(pair.unpack_spins(y, 0), ((0.0, 0.0, 0.5), (0.0, 0.0, 0.1)))
        self.assertStatus(lib.MYB_E_DOMAIN, pair.add_pn_pair, 0, 1, 1 << 20)


GM_EARTH = 3.986004418e14
