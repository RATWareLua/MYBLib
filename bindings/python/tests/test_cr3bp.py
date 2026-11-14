# SPDX-License-Identifier: BSD-3-Clause
"""The CR3BP through the package against the same C calls through cffi, bit for bit: the
Lagrange points, a Lyapunov orbit corrected by symmetry and by multiple shooting, its
monodromy, stability and manifold, a section, a torus about it with its spectrum and
manifold, the frames of two bodies, the term of a model; failures as myb.Error."""
import myb
from myb import ffi, lib

from common import Case, doubles, ok, rows

MU = 1.215058560962404e-02
LYAPUNOV = (8.2522726142834313e-01, 0.0, 0.0, 0.0, 1.0748883813645686e-01, 0.0)


def raw_correct(y0):
    correction = ffi.new("myb_cr3bp_correction*", {"symmetry": lib.MYB_CR3BP_SYMMETRY_XZ})
    y = ffi.new("double[6]", y0)
    period = ffi.new("double*")
    ok(lib.myb_cr3bp_correct(MU, correction, y, period))
    return tuple(y), period[0]


class Problem(Case):
    def test_points_orbits_and_manifolds(self):
        points = ffi.new("double[5][3]")
        ok(lib.myb_cr3bp_lagrange(MU, points))
        self.assertSame(myb.cr3bp_lagrange(MU), [tuple(points[k]) for k in range(5)])
        y0, period = myb.cr3bp_correct(MU, LYAPUNOV)
        self.assertSame((y0, period), raw_correct(LYAPUNOV))
        self.assertSame(myb.cr3bp_jacobi(MU, y0), lib.myb_cr3bp_jacobi(MU, ffi.new("double[6]",
                                                                                   y0)))
        monodromy = myb.cr3bp_monodromy(MU, y0, period)
        raw = ffi.new("double[36]")
        ok(lib.myb_cr3bp_monodromy(MU, ffi.new("double[6]", y0), period, raw))
        self.assertSame(monodromy, rows(raw, 6, 6))
        stability = myb.cr3bp_stability(monodromy)
        re = ffi.new("double[6]")
        im = ffi.new("double[6]")
        index = ffi.new("double*")
        ok(lib.myb_cr3bp_stability(raw, re, im, index))
        self.assertSame(tuple(stability), (doubles(re, 6), doubles(im, 6), index[0]))
        self.assertGreater(stability.index, 1.0)
        on = myb.cr3bp_manifold(MU, y0, period, 0.25, True, 1, 1e-6)
        y = ffi.new("double[6]")
        ok(lib.myb_cr3bp_manifold(MU, ffi.new("double[6]", y0), period, 0.25, 1, 1, 1e-6, y))
        self.assertSame(on, doubles(y, 6))
        # to the plane y = 0, crossed with y falling
        crossing = myb.cr3bp_propagate(MU, y0, period, (1, 0.0, -1))
        elapsed = ffi.new("double*")
        crossed = ffi.new("int*")
        ok(lib.myb_cr3bp_propagate(MU, ffi.new("double[6]", y0), period,
                                   ffi.new("myb_cr3bp_section*", [1, 0.0, -1]), y, elapsed,
                                   crossed))
        self.assertSame(crossing, (doubles(y, 6), elapsed[0], crossed[0] != 0))
        self.assertTrue(crossing[2])
        halo = myb.cr3bp_richardson(MU, 1, 0.02, 1)
        ok(lib.myb_cr3bp_richardson(MU, 1, 0.02, 1, y, elapsed))
        self.assertSame(halo, (doubles(y, 6), elapsed[0]))
        a = myb.cr3bp_acceleration(MU, y0[:3], y0[3:])
        dadx, dadv = myb.cr3bp_jacobian(MU, y0[:3])
        self.assertEqual((len(a), len(dadx), len(dadv)), (3, 3, 3))
        model = myb.ForceModel()
        model.add_cr3bp("cr3bp", MU)
        term = ffi.new("double[3]")
        ok(lib.myb_cr3bp_term(ffi.new("double*", MU), 0.0, ffi.new("double[3]", y0[:3]),
                              ffi.new("double[3]", y0[3:]), term))
        self.assertSame(a, doubles(term, 3))
        self.assertEqual(model.acceleration(0.0, y0[:3], y0[3:]), a)  # the sum makes -0 +0
        self.assertStatus(lib.MYB_E_DOMAIN, myb.cr3bp_lagrange, 0.7)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.cr3bp_correct, MU, (0.82, 0.1, 0, 0, 0.16, 0))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.cr3bp_propagate, MU, y0, 1.0, (7, 0.0, 0))

    def test_families_and_tori(self):
        y0, period = raw_correct(LYAPUNOV)
        orbit = myb.Cr3bpOrbit.from_state(MU, y0, period, 4)
        orbit.correct()
        orbit.continue_(1e-3)
        raw = ffi.new("myb_cr3bp_orbit*")
        ok(lib.myb_cr3bp_orbit_from_state(MU, ffi.new("double[6]", y0), period, 4, raw))
        ok(lib.myb_cr3bp_orbit_correct(MU, raw, 0.0, 0))
        ok(lib.myb_cr3bp_orbit_continue(MU, raw, 1e-3, 0.0, 0))
        self.assertSame((orbit.states, orbit.period, orbit.tangent),
                        (rows(raw.states, 4, 6), raw.period, doubles(raw.tangent, 25)))
        again = myb.Cr3bpOrbit(MU, orbit.states, orbit.period)
        self.assertFalse(again.has_tangent)
        self.assertStatus(lib.MYB_E_DOMAIN, again.continue_, 1e-3)
        # a torus about the orbit, its spectrum and a state on its manifold
        torus = myb.Cr3bpTorus.from_orbit(MU, y0, period, 9, 0, 1e-3)
        torus.continue_("jacobi", 0.0)
        spectrum = torus.stability()
        manifold = torus.manifold(spectrum, 0.5, True, 1, 1e-6)
        made = ffi.new("myb_cr3bp_torus*")
        ok(lib.myb_cr3bp_torus_from_orbit(MU, ffi.new("double[6]", y0), period, 9, 0, 1e-3,
                                          made))
        ok(lib.myb_cr3bp_torus_continue(MU, made, lib.MYB_CR3BP_TORUS_JACOBI, 0.0, 0.0, 0))
        s = ffi.new("myb_cr3bp_torus_spectrum*")
        ok(lib.myb_cr3bp_torus_stability(MU, made, 0.0, s))
        y = ffi.new("double[6]")
        ok(lib.myb_cr3bp_torus_manifold(made, s, 0.5, 1, 1, 1e-6, y))
        self.assertSame((torus.states, torus.period, torus.rotation, torus.defect),
                        (rows(made.states, 9, 6), made.period, made.rotation, made.defect))
        self.assertSame((spectrum.re, spectrum.im, spectrum.index, spectrum.unstable,
                         spectrum.unstable_direction),
                        (doubles(s.re, s.count), doubles(s.im, s.count), s.index, s.unstable,
                         doubles(s.unstable_direction, 54)))
        self.assertSame((spectrum.normal_rotation, spectrum.normal_eigenvalue,
                         spectrum.normal_re), (s.normal_rotation, tuple(s.normal_eigenvalue),
                                               doubles(s.normal_re, 54)))
        self.assertSame(manifold, doubles(y, 6))
        self.assertEqual(spectrum.count, 54)
        lib.myb_cr3bp_torus_state(made, 1.0, y)
        self.assertSame(torus.state(1.0), doubles(y, 6))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.Cr3bpTorus.from_orbit, MU, y0, period, 4, 0,
                          1e-3)

    def test_frames_of_two_bodies(self):
        earth = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        moon = (3.844e8, 0.0, 0.0, 0.0, 1.02e3, 0.0)
        state = (3.0e8, 1.0e7, 2.0e6, 10.0, 900.0, 1.0)
        rotating = myb.cr3bp_from_ephemeris(earth, moon, 3.986e14, 4.9e12, state)
        out = ffi.new("double[6]")
        ok(lib.myb_cr3bp_from_ephemeris(ffi.new("double[6]", earth), ffi.new("double[6]", moon),
                                        3.986e14, 4.9e12, ffi.new("double[6]", state), out))
        self.assertSame(rotating, doubles(out, 6))
        back = myb.cr3bp_to_ephemeris(earth, moon, 3.986e14, 4.9e12, rotating)
        self.assertLess(max(abs(a - b) for a, b in zip(back, state)), 1e-3)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.cr3bp_from_ephemeris, earth, earth, 3.986e14,
                          4.9e12, state)
