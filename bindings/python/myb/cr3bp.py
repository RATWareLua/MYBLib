# SPDX-License-Identifier: BSD-3-Clause
"""The circular restricted three-body problem (myb/cr3bp.h), in its rotating frame and
units: the distance of the primaries, their mass and their mean motion are 1, mu the mass of
the smaller over the sum. Lagrange points, the Jacobi constant, periodic orbits (corrected by
symmetry or by multiple shooting, continued along their families), monodromy, stability and
manifolds, quasi-periodic tori with their stability and manifolds, propagation to sections,
the frames of two real bodies. The term of a force model is ForceModel.add_cr3bp.

Matrices come out as lists of rows (the monodromy: six rows of six)."""
from collections import namedtuple

from ._core import _check, _doubles, _flat, _rows, _tuple
from ._native import ffi, lib

__all__ = [
    "Cr3bpOrbit", "Cr3bpTorus", "Cr3bpTorusSpectrum", "cr3bp_acceleration", "cr3bp_correct",
    "cr3bp_from_ephemeris", "cr3bp_jacobi", "cr3bp_jacobian", "cr3bp_lagrange",
    "cr3bp_manifold", "cr3bp_monodromy", "cr3bp_propagate", "cr3bp_richardson",
    "cr3bp_stability", "cr3bp_to_ephemeris",
]

_symmetries = {"xz": lib.MYB_CR3BP_SYMMETRY_XZ, "x": lib.MYB_CR3BP_SYMMETRY_X}
_families = {"jacobi": lib.MYB_CR3BP_TORUS_JACOBI, "period": lib.MYB_CR3BP_TORUS_PERIOD,
             "rotation": lib.MYB_CR3BP_TORUS_ROTATION}

Cr3bpStability = namedtuple("Cr3bpStability", "re im index")


def cr3bp_acceleration(mu, x, v):
    a = ffi.new("double[3]")
    lib.myb_cr3bp_acceleration(float(mu), _doubles(x, 3, "x"), _doubles(v, 3, "v"), a)
    return tuple(a)


def cr3bp_jacobian(mu, x):
    """(d a / d x, d a / d v), three rows each."""
    dadx = ffi.new("double[9]")
    dadv = ffi.new("double[9]")
    lib.myb_cr3bp_jacobian(float(mu), _doubles(x, 3, "x"), dadx, dadv)
    return _rows(dadx, 3, 3), _rows(dadv, 3, 3)


def cr3bp_jacobi(mu, y):
    return lib.myb_cr3bp_jacobi(float(mu), _doubles(y, 6, "y"))


def cr3bp_lagrange(mu):
    """The five Lagrange points L1 .. L5: a list of five positions."""
    points = ffi.new("double[5][3]")
    _check(lib.myb_cr3bp_lagrange(float(mu), points))
    return [tuple(points[k]) for k in range(5)]


def cr3bp_correct(mu, y0, symmetry="xz", fixed=0, crossing=1, tolerance=0.0, max_iterations=0):
    """A periodic orbit of the symmetry ("xz" or "x") through a first guess y0 by single
    shooting on the half period: (y0 corrected, period)."""
    correction = ffi.new("myb_cr3bp_correction*")
    correction.symmetry = _symmetries[symmetry] if isinstance(symmetry, str) else int(symmetry)
    correction.fixed = int(fixed)
    correction.crossing = int(crossing)
    correction.tolerance = float(tolerance)
    correction.max_iterations = int(max_iterations)
    y = _doubles(y0, 6, "y0")
    period = ffi.new("double*")
    _check(lib.myb_cr3bp_correct(float(mu), correction, y, period))
    return tuple(y), period[0]


def cr3bp_monodromy(mu, y0, period):
    """The state transition matrix over one period: six rows."""
    out = ffi.new("double[36]")
    _check(lib.myb_cr3bp_monodromy(float(mu), _doubles(y0, 6, "y0"), float(period), out))
    return _rows(out, 6, 6)


def cr3bp_stability(monodromy):
    """The eigenvalues of a monodromy matrix (re, im: six each) and the stability index."""
    re = ffi.new("double[6]")
    im = ffi.new("double[6]")
    index = ffi.new("double*")
    _check(lib.myb_cr3bp_stability(_doubles(monodromy, 36, "monodromy"), re, im, index))
    return Cr3bpStability(tuple(re), tuple(im), index[0])


def cr3bp_manifold(mu, y0, period, phase, unstable, sign, epsilon):
    """A state on the unstable (unstable true) or stable manifold of a periodic orbit."""
    y = ffi.new("double[6]")
    _check(lib.myb_cr3bp_manifold(float(mu), _doubles(y0, 6, "y0"), float(period), float(phase),
                                  1 if unstable else 0, int(sign), float(epsilon), y))
    return tuple(y)


def cr3bp_richardson(mu, point, az, branch=1):
    """The third-order halo orbit of Richardson about L1 or L2: (y0, period)."""
    y = ffi.new("double[6]")
    period = ffi.new("double*")
    _check(lib.myb_cr3bp_richardson(float(mu), int(point), float(az), int(branch), y, period))
    return tuple(y), period[0]


def cr3bp_from_ephemeris(body1, body2, gm1, gm2, state):
    """A state (m, m/s) in the rotating frame and units of two real bodies."""
    out = ffi.new("double[6]")
    _check(lib.myb_cr3bp_from_ephemeris(_doubles(body1, 6, "body1"), _doubles(body2, 6, "body2"),
                                        float(gm1), float(gm2), _doubles(state, 6, "state"), out))
    return tuple(out)


def cr3bp_to_ephemeris(body1, body2, gm1, gm2, rotating):
    out = ffi.new("double[6]")
    _check(lib.myb_cr3bp_to_ephemeris(_doubles(body1, 6, "body1"), _doubles(body2, 6, "body2"),
                                      float(gm1), float(gm2), _doubles(rotating, 6, "rotating"),
                                      out))
    return tuple(out)


def cr3bp_propagate(mu, y0, duration, section=None):
    """The flow from y0 for `duration`, stopping at the first crossing of a section (axis,
    value, direction) when one is given: (y, elapsed, crossed)."""
    plane = ffi.NULL
    if section is not None:
        axis, value, direction = section
        plane = ffi.new("myb_cr3bp_section*", {"axis": int(axis), "value": float(value),
                                               "direction": int(direction)})
    y = ffi.new("double[6]")
    elapsed = ffi.new("double*")
    crossed = ffi.new("int*")
    _check(lib.myb_cr3bp_propagate(float(mu), _doubles(y0, 6, "y0"), float(duration), plane, y,
                                   elapsed, crossed))
    return tuple(y), elapsed[0], crossed[0] != 0


class Cr3bpOrbit:
    """A periodic orbit of any shape cut into patches (myb_cr3bp_orbit), corrected by
    multiple shooting and continued along its family by pseudo-arclength. The methods change
    the orbit in place, as the C functions do."""

    def __init__(self, mu, states, period, tangent=None):
        states = _flat(states)
        if len(states) % 6 or not states:
            raise ValueError("states: six numbers a patch")
        self.mu = float(mu)
        self.raw = ffi.new("myb_cr3bp_orbit*")
        self.raw.patches = len(states) // 6
        self.raw.states[0:len(states)] = states
        self.raw.period = float(period)
        if tangent is not None:
            tangent = _flat(tangent)
            self.raw.tangent[0:len(tangent)] = tangent
            self.raw.has_tangent = 1

    @classmethod
    def from_state(cls, mu, y0, period, patches):
        """The cuts of the orbit through y0 with the given period (myb_cr3bp_orbit_from_state)."""
        made = cls.__new__(cls)
        made.mu = float(mu)
        made.raw = ffi.new("myb_cr3bp_orbit*")
        _check(lib.myb_cr3bp_orbit_from_state(made.mu, _doubles(y0, 6, "y0"), float(period),
                                              int(patches), made.raw))
        return made

    def correct(self, tolerance=0.0, max_iterations=0):
        _check(lib.myb_cr3bp_orbit_correct(self.mu, self.raw, float(tolerance),
                                           int(max_iterations)))

    def continue_(self, step, tolerance=0.0, max_iterations=0):
        """The next orbit of the family a distance `step` along the tangent
        (myb_cr3bp_orbit_continue); on MYB_E_NOCONV the orbit is left as it was."""
        _check(lib.myb_cr3bp_orbit_continue(self.mu, self.raw, float(step), float(tolerance),
                                            int(max_iterations)))

    @property
    def patches(self):
        return self.raw.patches

    @property
    def states(self):
        return _rows(self.raw.states, self.raw.patches, 6)

    @property
    def period(self):
        return self.raw.period

    @property
    def has_tangent(self):
        return self.raw.has_tangent != 0

    @property
    def tangent(self):
        return _tuple(self.raw.tangent, 6 * self.raw.patches + 1)


class Cr3bpTorusSpectrum:
    """The stability of a torus (myb_cr3bp_torus_spectrum): every eigenvalue of the map of
    its circle, their classes, the hyperbolic directions and the normal pair of an elliptic
    torus. raw: the structure of the C interface."""

    def __init__(self, raw, points):
        self.raw = raw
        self._points = points

    def __getattr__(self, name):
        if name in ("count", "trivial", "elliptic", "hyperbolic", "unresolved", "index",
                    "unstable", "stable", "residual", "tail", "normal_rotation",
                    "normal_frequency", "normal_residual"):
            return getattr(self.raw, name)
        raise AttributeError(name)

    @property
    def re(self):
        return _tuple(self.raw.re, self.raw.count)

    @property
    def im(self):
        return _tuple(self.raw.im, self.raw.count)

    @property
    def unstable_direction(self):
        return _tuple(self.raw.unstable_direction, 6 * self._points)

    @property
    def stable_direction(self):
        return _tuple(self.raw.stable_direction, 6 * self._points)

    @property
    def normal_eigenvalue(self):
        """(re, im) of exp(i sigma) of the normal pair."""
        return tuple(self.raw.normal_eigenvalue)

    @property
    def normal_re(self):
        return _tuple(self.raw.normal_re, 6 * self._points)

    @property
    def normal_im(self):
        return _tuple(self.raw.normal_im, 6 * self._points)


class Cr3bpTorus:
    """A quasi-periodic torus as an invariant circle of `points` points (myb_cr3bp_torus):
    found about a periodic orbit, continued along a family ("jacobi", "period" or
    "rotation" kept), its points at any angle, its stability and its manifolds. The methods
    change the torus in place, as the C functions do."""

    def __init__(self, mu):
        self.mu = float(mu)
        self.raw = ffi.new("myb_cr3bp_torus*")

    @classmethod
    def from_orbit(cls, mu, y0, period, points, pair=0, size=1e-3):
        """The smallest torus about the periodic orbit through y0
        (myb_cr3bp_torus_from_orbit)."""
        made = cls(mu)
        _check(lib.myb_cr3bp_torus_from_orbit(made.mu, _doubles(y0, 6, "y0"), float(period),
                                              int(points), int(pair), float(size), made.raw))
        return made

    def continue_(self, family, step, tolerance=0.0, max_iterations=0):
        """The torus of the family a distance `step` along its tangent (0: corrects it)."""
        family = _families[family] if isinstance(family, str) else int(family)
        _check(lib.myb_cr3bp_torus_continue(self.mu, self.raw, family, float(step),
                                            float(tolerance), int(max_iterations)))

    def state(self, angle):
        """The point of the circle at an angle (rad)."""
        y = ffi.new("double[6]")
        lib.myb_cr3bp_torus_state(self.raw, float(angle), y)
        return tuple(y)

    def stability(self, tolerance=0.0):
        """The spectrum of the torus (myb_cr3bp_torus_stability)."""
        spectrum = ffi.new("myb_cr3bp_torus_spectrum*")
        _check(lib.myb_cr3bp_torus_stability(self.mu, self.raw, float(tolerance), spectrum))
        return Cr3bpTorusSpectrum(spectrum, self.raw.points)

    def manifold(self, spectrum, angle, unstable, sign, epsilon):
        """A state on the unstable or stable manifold of the torus."""
        y = ffi.new("double[6]")
        _check(lib.myb_cr3bp_torus_manifold(self.raw, spectrum.raw, float(angle),
                                            1 if unstable else 0, int(sign), float(epsilon), y))
        return tuple(y)

    @property
    def points(self):
        return self.raw.points

    @property
    def states(self):
        return _rows(self.raw.states, self.raw.points, 6)

    @property
    def period(self):
        return self.raw.period

    @property
    def rotation(self):
        return self.raw.rotation

    @property
    def defect(self):
        return self.raw.defect

    @property
    def has_tangent(self):
        return self.raw.has_tangent != 0

    @property
    def tangent(self):
        return _tuple(self.raw.tangent, 6 * self.raw.points + 2)
