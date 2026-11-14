# SPDX-License-Identifier: BSD-3-Clause
"""Attitude and observables through the package against the same C calls through cffi, bit
for bit: quaternions, Euler angles, principal axes, the rates, the torques, the exact free
rotation and the splitting integrator with a torque in Python; the light time with a body
that delays it, a round trip, aberration, deflection, the plasma; failures as myb.Error."""
import math

import myb
from myb import ffi, lib

from common import Case, doubles, ok

Q = (0.9, 0.1, -0.3, 0.2)
MOMENTS = (1.0, 2.0, 3.0)


class Attitude(Case):
    def test_quaternions_and_angles(self):
        out = ffi.new("double[9]")
        lib.myb_quat_to_matrix(ffi.new("double[4]", Q), out)
        self.assertSame(myb.quat_to_matrix(Q), doubles(out, 9))
        q = ffi.new("double[4]")
        lib.myb_quat_from_matrix(out, q)
        self.assertSame(myb.quat_from_matrix(doubles(out, 9)), doubles(q, 4))
        lib.myb_quat_multiply(ffi.new("double[4]", Q), q, q)
        self.assertSame(myb.quat_multiply(Q, myb.quat_from_matrix(myb.quat_to_matrix(Q))),
                        doubles(q, 4))
        unit, norm = myb.quat_normalize(Q)
        made = ffi.new("double[4]", Q)
        self.assertSame(norm, lib.myb_quat_normalize(made))
        self.assertSame(unit, doubles(made, 4))
        phi = myb.quat_to_rotation_vector(unit)
        self.assertLess(max(abs(a - b) for a, b in zip(myb.quat_from_rotation_vector(phi), unit)),
                        1e-15)
        self.assertEqual(len(myb.quat_slerp(unit, (1.0, 0.0, 0.0, 0.0), 0.5)), 4)
        self.assertEqual(len(myb.quat_derivative(unit, (0.0, 0.0, 1.0))), 4)
        self.assertEqual(len(myb.quat_rotate(unit, (1.0, 2.0, 3.0))), 3)
        self.assertEqual(myb.quat_conjugate(Q), (0.9, -0.1, 0.3, -0.2))
        angles = myb.euler_from_matrix((3, 1, 3), myb.quat_to_matrix(unit))
        a = ffi.new("double[3]")
        ok(lib.myb_euler_from_matrix(ffi.new("int[3]", [3, 1, 3]),
                                     ffi.new("double[9]", myb.quat_to_matrix(unit)), a))
        self.assertSame(angles, doubles(a, 3))
        self.assertEqual(len(myb.euler_to_matrix((3, 2, 1), angles)), 9)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.euler_to_matrix, (3, 3, 1), angles)
        moments, axes = myb.principal_axes((2.0, 0.1, 0.0, 0.1, 3.0, 0.0, 0.0, 0.0, 1.0))
        self.assertEqual(len(axes), 9)
        self.assertLess(moments[0], moments[2])
        self.assertStatus(lib.MYB_E_DOMAIN, myb.principal_axes, (1.0, 2.0, 0, 0, 1, 0, 0, 0, 1))

    def test_dynamics_and_torques(self):
        unit = myb.quat_normalize(Q)[0]
        l = (0.1, 0.5, 2.0)
        dq, dl = myb.rigid_body_rates(MOMENTS, unit, l, (0.0, 0.0, 1e-3))
        q = ffi.new("double[4]")
        m = ffi.new("double[3]")
        lib.myb_rigid_body_rates(ffi.new("double[3]", MOMENTS), ffi.new("double[4]", unit),
                                 ffi.new("double[3]", l), ffi.new("double[3]", [0.0, 0.0, 1e-3]),
                                 q, m)
        self.assertSame((dq, dl), (doubles(q, 4), doubles(m, 3)))
        free = myb.FreeRotation(MOMENTS, unit, l)
        raw = ffi.new("myb_free_rotation*")
        ok(lib.myb_free_rotation_init(ffi.new("double[3]", MOMENTS), ffi.new("double[4]", unit),
                                      ffi.new("double[3]", l), raw))
        lib.myb_free_rotation_state(raw, 100.0, q, m)
        self.assertSame(free.state(100.0), (doubles(q, 4), doubles(m, 3)))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.FreeRotation, (0.0, 1.0, 1.0), unit, l)
        # the splitting integrator with a torque in Python
        def none(t, q):
            return (0.0, 0.0, 0.0)

        split = myb.rotation_split(MOMENTS, none, 4, 0.0, 0.1, 1000, unit, l)

        @ffi.callback("myb_status(void*, double, const double[4], double[3])")
        def raw_none(context, t, q, out):
            out[0:3] = [0.0, 0.0, 0.0]
            return lib.MYB_OK

        attitude = ffi.new("double[4]", unit)
        momentum = ffi.new("double[3]", l)
        ok(lib.myb_rotation_split(ffi.new("double[3]", MOMENTS), raw_none, ffi.NULL, 4, 0.0, 0.1,
                                  1000, attitude, momentum))
        self.assertSame(split, (doubles(attitude, 4), doubles(momentum, 3)))
        gravity = myb.rotation_split(MOMENTS, lambda t, q: myb.torque_point_mass(
            1e-6, myb.quat_rotate(q, (1.0, 0.0, 0.0)), MOMENTS), 2, 0.0, 0.1, 100, unit, l)
        self.assertNotEqual(gravity, myb.rotation_split(MOMENTS, none, 2, 0.0, 0.1, 100, unit, l))

        def broken(t, q):
            raise FloatingPointError("no torque")

        with self.assertRaises(FloatingPointError):
            myb.rotation_split(MOMENTS, broken, 2, 0.0, 0.1, 10, unit, l)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.rotation_split, MOMENTS, none, 3, 0.0, 0.1, 10,
                          unit, l)
        torque = ffi.new("double[3]")
        lib.myb_torque_field(0.1, ffi.new("double[3]", [1.0, 2.0, 3.0]),
                             ffi.new("double[3]", [0.0, -1.0, 0.5]), torque)
        self.assertSame(myb.torque_field(0.1, (1.0, 2.0, 3.0), (0.0, -1.0, 0.5)),
                        doubles(torque, 3))
        rates = myb.core_mantle_rates((0.39, 0, 0, 0, 0.392, 0, 0, 0, 0.394), None,
                                      (0.001, 0.0012), 1e-9, (0.0, 0.0, 2.66e-6),
                                      (0.0, 1e-8, 2.66e-6), (0.0, 0.0, 0.0))
        self.assertEqual(len(rates), 2)
        self.assertEqual(len(myb.torque_figure_figure(1.0, (1.0,) * 9, (2.0,) * 9,
                                                      (1e8, 0.0, 0.0))), 3)


class Observables(Case):
    C = 299792458.0

    def test_light_time_and_directions(self):
        sun = 1.32712440018e20

        def target(seconds):
            return (1.5e11, 1.0e9 + 1.0e3 * seconds, 0.0, 0.0, 1.0e3, 0.0)

        model = myb.LightModel(self.C, bodies=[(sun, lambda s: (0.0,) * 6)])
        t = myb.time_from_seconds(1000.0)
        observer = (1.0e7, 0.0, 0.0, 0.0, 3.0e4, 0.0)
        path = model.light_time(target, t, observer)

        @ffi.callback("myb_status(void*, myb_time, double[6])")
        def raw_target(context, instant, y):
            y[0:6] = list(target(lib.myb_time_seconds(instant)))
            return lib.MYB_OK

        @ffi.callback("myb_status(void*, myb_time, double[6])")
        def raw_sun(context, instant, y):
            y[0:6] = [0.0] * 6
            return lib.MYB_OK

        bodies = ffi.new("myb_light_body[1]", [{"gm": sun, "state": raw_sun}])
        raw = ffi.new("myb_light_model*", {"c": self.C, "gamma": 1.0, "bodies": bodies,
                                           "count": 1})
        made = ffi.new("myb_light_path*")
        ok(lib.myb_light_time(raw, raw_target, ffi.NULL, t, ffi.new("double[6]", observer),
                              lib.MYB_LIGHT_RECEPTION, made))
        self.assertSame(path, myb.LightPath(made.t_target, tuple(made.target), made.light_time,
                                            made.range, made.delay, made.range_rate,
                                            made.iterations))
        self.assertGreater(path.delay, 0.0)
        total, bounce, transmit = model.round_trip(lambda s: observer, target, t)
        self.assertLess(abs(total - 2.0 * path.light_time), 1e-3)

        def broken(seconds):
            raise KeyError("no such body")

        with self.assertRaises(KeyError):
            model.light_time(broken, t, observer)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.LightModel(-1.0).light_time, target, t,
                          observer)
        x1, x2 = (1.5e11, 0.0, 0.0), (-1.0e11, 1.0e9, 0.0)
        self.assertSame(myb.shapiro_delay(sun, self.C, 1.0, x1, x2, (0.0, 0.0, 0.0)),
                        lib.myb_shapiro_delay(sun, self.C, 1.0, ffi.new("double[3]", x1),
                                              ffi.new("double[3]", x2),
                                              ffi.new("double[3]", [0.0, 0.0, 0.0])))
        p = (0.6, 0.8, 0.0)
        out = ffi.new("double[3]")
        lib.myb_aberration(ffi.new("double[3]", p), ffi.new("double[3]", [3e4, 0.0, 0.0]),
                           self.C, out)
        self.assertSame(myb.aberration(p, (3e4, 0.0, 0.0), self.C), doubles(out, 3))
        self.assertEqual(len(myb.light_deflection(sun, self.C, 1.0, p, (1.0, 0.0, 0.0),
                                                  (0.0, 1.0, 0.0), 1.5e11, 1e-9)), 3)
        ra, dec = myb.direction_radec((0.0, 1.0, 1.0))
        self.assertSame((ra, dec), (math.pi / 2, math.pi / 4))
        self.assertSame(myb.plasma_group_delay(1e17, 1.5e9, 1.602176634e-19, 9.1093837015e-31,
                                               8.8541878128e-12),
                        lib.myb_plasma_group_delay(1e17, 1.5e9, 1.602176634e-19,
                                                   9.1093837015e-31, 8.8541878128e-12))
        column = myb.plasma_column(6.957e8, [(1e14, 6.0), (1e12, 2.0)], (1.5e11, 0.0, 0.0),
                                   (-1.0e11, 3.0e9, 0.0))
        self.assertGreater(column, 0.0)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.plasma_column, 6.957e8, [(1e14, 6.0)],
                          (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0))
