# SPDX-License-Identifier: BSD-3-Clause
"""Observables (myb/observe.h): the light time with the relativistic delay, round trips,
aberration, light deflection, right ascension and declination, the delay in a plasma.

The participants are providers of states — bodies of an Ephemeris, Providers or callables
f(seconds) -> six numbers — in one inertial frame on the time of that frame."""
from collections import namedtuple

from ._core import _Calls, _check, _doubles, _owned, _settle, _state_provider
from ._native import ffi, lib

__all__ = [
    "LightModel", "LightPath", "aberration", "direction_radec", "light_deflection",
    "plasma_column", "plasma_group_delay", "shapiro_delay",
]

LightPath = namedtuple("LightPath", "t_target target light_time range delay range_rate "
                                    "iterations")

_directions = {"reception": lib.MYB_LIGHT_RECEPTION, "transmission": lib.MYB_LIGHT_TRANSMISSION}


class LightModel:
    """How light travels (myb_light_model): the speed c, the PPN gamma and the bodies that
    delay it, pairs (gm, provider of states); no bodies: Newtonian light time."""

    def __init__(self, c, gamma=1.0, bodies=()):
        self._calls = _Calls()
        bodies = list(bodies)
        self._bodies = ffi.new("myb_light_body[]", max(len(bodies), 1))
        self._kept = []
        for k, (gm, state) in enumerate(bodies):
            provider = _state_provider(state, self._calls)
            self._bodies[k].gm = float(gm)
            self._bodies[k].state = provider.function
            self._bodies[k].context = provider.context
            self._kept.append(provider)
        self.raw = ffi.new("myb_light_model*", {"c": float(c), "gamma": float(gamma),
                                                "count": len(bodies)})
        if bodies:
            self.raw.bodies = self._bodies

    def light_time(self, target, t, observer, direction="reception"):
        """The one-way light time between the observer (state at the instant t) and the
        target (myb_light_time): a LightPath."""
        provider = _state_provider(target, self._calls)
        path = ffi.new("myb_light_path*")
        direction = _directions[direction] if isinstance(direction, str) else int(direction)
        _settle(lib.myb_light_time(self.raw, provider.function, provider.context, t,
                                   _doubles(observer, 6, "observer"), direction, path),
                self, provider)
        return LightPath(_owned("myb_time", path.t_target), tuple(path.target), path.light_time,
                         path.range, path.delay, path.range_rate, path.iterations)

    def round_trip(self, station, target, t):
        """A round trip received by the station at the instant t (myb_round_trip):
        (total time, instant of the bounce, instant of the transmission)."""
        one = _state_provider(station, self._calls)
        two = _state_provider(target, self._calls)
        total = ffi.new("double*")
        bounce = ffi.new("myb_time*")
        transmit = ffi.new("myb_time*")
        _settle(lib.myb_round_trip(self.raw, one.function, one.context, two.function,
                                   two.context, t, total, bounce, transmit), self, one, two)
        return total[0], bounce[0], transmit[0]


def shapiro_delay(gm, c, gamma, x1, x2, xb):
    return lib.myb_shapiro_delay(float(gm), float(c), float(gamma), _doubles(x1, 3, "x1"),
                                 _doubles(x2, 3, "x2"), _doubles(xb, 3, "xb"))


def aberration(p, v, c):
    out = ffi.new("double[3]")
    lib.myb_aberration(_doubles(p, 3, "p"), _doubles(v, 3, "v"), float(c), out)
    return tuple(out)


def light_deflection(gm, c, gamma, p, q, e, em, dlim):
    out = ffi.new("double[3]")
    lib.myb_light_deflection(float(gm), float(c), float(gamma), _doubles(p, 3, "p"),
                             _doubles(q, 3, "q"), _doubles(e, 3, "e"), float(em), float(dlim),
                             out)
    return tuple(out)


def direction_radec(p):
    """(right ascension in [0, 2 pi), declination) of a direction."""
    ra = ffi.new("double*")
    dec = ffi.new("double*")
    lib.myb_direction_radec(_doubles(p, 3, "p"), ra, dec)
    return ra[0], dec[0]


def plasma_group_delay(column, frequency, e, m_e, eps0):
    return lib.myb_plasma_group_delay(float(column), float(frequency), float(e), float(m_e),
                                      float(eps0))


def plasma_column(radius, terms, p1, p2):
    """The electron column (m^-2) of the profile sum density (radius / r)^power along the
    segment p1 -> p2; terms: pairs (density, power)."""
    terms = [(float(d), float(p)) for d, p in terms]
    if not 1 <= len(terms) <= lib.MYB_PLASMA_TERMS:
        raise ValueError("1 .. %d terms" % lib.MYB_PLASMA_TERMS)
    profile = ffi.new("myb_plasma_profile*")
    profile.radius = float(radius)
    profile.count = len(terms)
    for k, (density, power) in enumerate(terms):
        profile.density[k] = density
        profile.power[k] = power
    out = ffi.new("double*")
    _check(lib.myb_plasma_column(profile, _doubles(p1, 3, "p1"), _doubles(p2, 3, "p2"), out))
    return out[0]
