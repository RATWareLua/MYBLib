# SPDX-License-Identifier: BSD-3-Clause
"""Targeting and what goes with it: the B-plane, patched conics and shooting
(myb/transfer.h), Lambert's problem (myb/lambert.h), optimal transfers with low thrust — the
indirect method, the solver of nonlinear programs, Sims and Flanagan, chains of legs with
flybys (myb/optimal.h) —, thrust arcs, Edelbaum and the Q-law (myb/thrust.h), initial orbits
(myb/iod.h), relative motion (myb/relative.h), closest approaches and the probability of
collision (myb/conjunction.h), the variational equations, covariances and the deterministic
random numbers (myb/uncert.h).

States are six numbers (m, m/s); matrices come out as lists of rows. The functions that
the solvers call — residuals, objectives, providers of states — are Python callables; an
exception raised in one comes out of the call. On MYB_E_NOCONV the solvers leave what they
reached: Error.partial holds it."""
from collections import namedtuple

from ._core import (Elements, NlpResult, _Calls, _check, _doubles, _flat, _nlp_options,
                    _nlp_result, _optional, _Owner, _pv, _rows, _series, _settle,
                    _state_of, _state_provider, _square, _tuple)
from ._native import ffi, lib

__all__ = [
    "Approach", "BPlane", "Chain", "ChainResult", "Conjunction", "Encounter", "IndirectPoint",
    "IndirectProblem", "IndirectSolution", "IodSolution", "LambertSolution", "NlpResult",
    "QLaw", "QLawArrival", "QLawFlight", "QLawTerm", "Rng", "SimsFlanagan",
    "SimsFlanaganResult", "Thrust", "Variational", "approach_find", "approach_orbits",
    "bplane", "chain_flyby", "collision_probability", "collision_probability_chan",
    "collision_probability_long", "collision_probability_max", "collision_rate",
    "conjunction_find", "conjunction_long", "covariance_rtn", "covariance_transform",
    "cw_stm", "departure_impulse", "edelbaum", "edelbaum_state", "encounter_plane",
    "equinoctial_rates", "flyby_turn", "hill_radius", "iod_double_r", "iod_gauss",
    "iod_gibbs", "iod_gooding", "iod_herrick_gibbs", "iod_laplace", "lambert", "nlp_solve",
    "relative_from_inertial", "relative_rendezvous", "relative_stm", "relative_to_inertial",
    "roe_from_states", "roe_to_position", "sample_statistics", "shoot", "sphere_of_influence",
]

BPlane = namedtuple("BPlane", "v_infinity c3 s t r b b_t b_r b_magnitude theta periapsis "
                              "eccentricity time_to_periapsis")
LambertSolution = namedtuple("LambertSolution", "v1 v2 revolutions branch")
IndirectSolution = namedtuple("IndirectSolution", "costate t1 mass1 smoothing residual tries "
                                                  "stages phase")
IndirectPoint = namedtuple("IndirectPoint", "y direction throttle switching hamiltonian costate")
SimsFlanaganResult = namedtuple("SimsFlanaganResult", "throttles mass1 mismatch nlp")
ChainResult = namedtuple("ChainResult", "x mass excess turn periapsis impulse mismatch first nlp")
QLawArrival = namedtuple("QLawArrival", "t y arrived")
IodSolution = namedtuple("IodSolution", "r v range iterations")
Approach = namedtuple("Approach", "t distance position velocity")
Encounter = namedtuple("Encounter", "miss sigma radius")
Conjunction = namedtuple("Conjunction", "approach covariance encounter probability")


# ---- Transfers ---------------------------------------------------------------------------------

def bplane(gm, r, v, pole=(0.0, 0.0, 1.0)):
    """The B-plane of the hyperbola through (r, v) about a body of gm (myb_bplane_compute)."""
    out = ffi.new("myb_bplane*")
    _check(lib.myb_bplane_compute(float(gm), _doubles(r, 3, "r"), _doubles(v, 3, "v"),
                                  _doubles(pole, 3, "pole"), out))
    return BPlane(out.v_infinity, out.c3, tuple(out.s), tuple(out.t), tuple(out.r), tuple(out.b),
                  out.b_t, out.b_r, out.b_magnitude, out.theta, out.periapsis, out.eccentricity,
                  out.time_to_periapsis)


def sphere_of_influence(a, mass_ratio):
    return lib.myb_sphere_of_influence(float(a), float(mass_ratio))


def hill_radius(a, e, mass_ratio):
    return lib.myb_hill_radius(float(a), float(e), float(mass_ratio))


def flyby_turn(gm, v_infinity, periapsis):
    return lib.myb_flyby_turn(float(gm), float(v_infinity), float(periapsis))


def departure_impulse(gm, radius, v_infinity):
    return lib.myb_departure_impulse(float(gm), float(radius), float(v_infinity))


def shoot(function, u, m, tolerance, max_iterations=0, jacobian=False):
    """Newton's method for f(u) = 0, n = len(u) unknowns, m equations (myb_shoot).
    function(u) -> m values; with jacobian=True function(u, want_jacobian) -> (values,
    m rows of n or None). Returns (u, iterations, residual)."""
    calls = _Calls()
    x = _doubles(u, len(_flat(u)), "u")
    n = len(x)
    m = int(m)

    @ffi.callback("myb_status(void*, const double*, double*, double*)")
    def residual(context, given, f, rows):
        try:
            point = _tuple(given, n)
            if jacobian:
                values, matrix = function(point, rows != ffi.NULL)
                if rows != ffi.NULL:
                    rows[0:m * n] = list(_doubles(matrix, m * n, "the Jacobian"))
            else:
                values = function(point)
            f[0:m] = list(_doubles(values, m, "the residuals"))
            return lib.MYB_OK
        except BaseException as error:  # noqa: B902 - nothing may cross the C code
            return calls.failed(error)

    options = ffi.new("myb_shoot_options*", {"tolerance": float(tolerance),
                                             "max_iterations": int(max_iterations)})
    work = ffi.new("double[]", max(lib.myb_shoot_work_size(n, m), 1))
    iterations = ffi.new("int*")
    left = ffi.new("double*")
    status = lib.myb_shoot(residual, ffi.NULL, n, m, x, options, work, iterations, left)
    _settle(status, calls, partial=lambda: (_tuple(x, n), iterations[0], left[0]))
    return _tuple(x, n), iterations[0], left[0]


def lambert(mu, r1, r2, tof, retrograde=False, max_revolutions=0):
    """Every arc of 0 .. max_revolutions revolutions joining r1 and r2 in tof (myb_lambert):
    a list of LambertSolution (v1, v2, revolutions, branch)."""
    capacity = 1 + 2 * max(int(max_revolutions), 0)
    out = ffi.new("myb_lambert_solution[]", capacity)
    count = ffi.new("size_t*")
    _check(lib.myb_lambert(float(mu), _doubles(r1, 3, "r1"), _doubles(r2, 3, "r2"), float(tof),
                           1 if retrograde else 0, int(max_revolutions), out, capacity, count))
    return [LambertSolution(tuple(s.v1), tuple(s.v2), s.revolutions, s.branch)
            for s in (out[k] for k in range(count[0]))]


# ---- Optimal transfers -------------------------------------------------------------------------

_goals = {"fuel": lib.MYB_INDIRECT_FUEL, "time": lib.MYB_INDIRECT_TIME}
_arrivals = {"rendezvous": lib.MYB_INDIRECT_RENDEZVOUS, "orbit": lib.MYB_INDIRECT_ORBIT}


class IndirectProblem:
    """A transfer by the indirect method (myb_indirect_problem): goal "fuel" or "time",
    arrival "rendezvous" or "orbit"; field: a ForceModel of everything but the engine (None:
    the point mass mu)."""

    def __init__(self, mu, thrust, exhaust, t0, start, mass0, t1, target, goal="fuel",
                 smoothing=0.0, tolerance=0.0, tries=0, arrival="rendezvous", phase=0.0,
                 field=None):
        self.raw = ffi.new("myb_indirect_problem*")
        raw = self.raw
        raw.mu = float(mu)
        raw.field = ffi.NULL if field is None else field.pointer
        raw.thrust = float(thrust)
        raw.exhaust = float(exhaust)
        raw.t0 = float(t0)
        raw.start = list(_doubles(start, 6, "start"))
        raw.mass0 = float(mass0)
        raw.t1 = float(t1)
        raw.target = list(_doubles(target, 6, "target"))
        raw.goal = _goals[goal] if isinstance(goal, str) else int(goal)
        raw.smoothing = float(smoothing)
        raw.tolerance = float(tolerance)
        raw.tries = int(tries)
        raw.arrival = _arrivals[arrival] if isinstance(arrival, str) else int(arrival)
        raw.phase = float(phase)
        self._kept = [field]

    def solve(self, guess=None):
        """The solution (myb_indirect_solve) from the eight costates of a guess, or from
        first guesses of a fixed sequence."""
        out = ffi.new("myb_indirect_solution*")
        status = lib.myb_indirect_solve(self.raw, _optional(guess, 8, "guess"), ffi.NULL, out)
        _settle(status, self, partial=lambda: _solution_of(out))
        return _solution_of(out)

    def sample(self, solution, times):
        """The solution at times in increasing order: a list of IndirectPoint."""
        made = ffi.new("myb_indirect_solution*")
        made.costate = list(_doubles(solution.costate, 8, "costate"))
        for name in ("t1", "mass1", "smoothing", "residual", "tries", "stages", "phase"):
            setattr(made, name, getattr(solution, name))
        given, count = _series(times, 1, "times")
        out = ffi.new("myb_indirect_point[]", max(count, 1))
        _settle(lib.myb_indirect_sample(self.raw, made, given, count, ffi.NULL, out), self)
        return [IndirectPoint(tuple(p.y), tuple(p.direction), p.throttle, p.switching,
                              p.hamiltonian, tuple(p.costate))
                for p in (out[k] for k in range(count))]


def _solution_of(out):
    return IndirectSolution(tuple(out.costate), out.t1, out.mass1, out.smoothing, out.residual,
                            out.tries, out.stages, out.phase)


def nlp_solve(function, x, equalities=0, inequalities=0, lower=None, upper=None,
              feasibility=0.0, optimality=0.0, max_outer=0, max_inner=0, weight=0.0):
    """Minimise f(x) with c_i(x) = 0 (the first `equalities`), c_i(x) <= 0 (the
    `inequalities` after them) and lower <= x <= upper (myb_nlp_solve). function(x,
    want_gradient) -> (f, c) or, when want_gradient, (f, c, gradient, Jacobian rows).
    Returns (x, multipliers, NlpResult)."""
    calls = _Calls()
    point = _doubles(x, len(_flat(x)), "x")
    n = len(point)
    k = int(equalities) + int(inequalities)

    @ffi.callback("myb_status(void*, const double*, double*, double*, double*, double*)")
    def objective(context, given, f, c, gradient, jacobian):
        try:
            want = gradient != ffi.NULL
            answer = function(_tuple(given, n), want)
            f[0] = float(answer[0])
            if k:
                c[0:k] = list(_doubles(answer[1], k, "the constraints"))
            if want:
                gradient[0:n] = list(_doubles(answer[2], n, "the gradient"))
                if k:
                    jacobian[0:k * n] = list(_doubles(answer[3], k * n, "the Jacobian"))
            return lib.MYB_OK
        except BaseException as error:  # noqa: B902 - nothing may cross the C code
            return calls.failed(error)

    problem = ffi.new("myb_nlp_problem*")
    problem.n = n
    problem.equalities = int(equalities)
    problem.inequalities = int(inequalities)
    low = _optional(lower, n, "lower")
    high = _optional(upper, n, "upper")
    problem.lower = low
    problem.upper = high
    problem.function = objective
    multipliers = ffi.new("double[]", max(k, 1))
    result = ffi.new("myb_nlp_result*")
    status = lib.myb_nlp_solve(problem, _nlp_options(feasibility, optimality, max_outer,
                                                     max_inner, weight), point, multipliers,
                               ffi.NULL, result)
    _settle(status, calls, partial=lambda: (_tuple(point, n), _tuple(multipliers, k),
                                            _nlp_result(result)))
    return _tuple(point, n), _tuple(multipliers, k), _nlp_result(result)


class SimsFlanagan:
    """A leg of Sims and Flanagan (myb_sims_flanagan): `segments` impulses between two
    states; field: a ForceModel between the impulses (None: the conics of mu)."""

    def __init__(self, mu, thrust, exhaust, t0, start, mass0, t1, target, segments, field=None):
        self.raw = ffi.new("myb_sims_flanagan*")
        raw = self.raw
        raw.mu = float(mu)
        raw.thrust = float(thrust)
        raw.exhaust = float(exhaust)
        raw.t0 = float(t0)
        raw.start = list(_doubles(start, 6, "start"))
        raw.mass0 = float(mass0)
        raw.t1 = float(t1)
        raw.target = list(_doubles(target, 6, "target"))
        raw.segments = int(segments)
        raw.field = ffi.NULL if field is None else field.pointer
        self._kept = [field]

    def mismatch(self, throttles, jacobian=False):
        """(mismatch, mass1, Jacobian rows or None) for the throttles (three a segment)."""
        n = 3 * self.raw.segments
        out = ffi.new("double[6]")
        mass1 = ffi.new("double*")
        rows = ffi.new("double[]", 6 * n) if jacobian else ffi.NULL
        _settle(lib.myb_sims_flanagan_mismatch(self.raw, _doubles(throttles, n, "throttles"), out,
                                               mass1, rows), self)
        return tuple(out), mass1[0], (_rows(rows, 6, n) if jacobian else None)

    def solve(self, throttles=None, feasibility=0.0, optimality=0.0, max_outer=0, max_inner=0,
              weight=0.0):
        """The throttles of the least propellant (myb_sims_flanagan_solve); None: the guess
        built in."""
        n = 3 * self.raw.segments
        given = ffi.new("double[]", n) if throttles is None else _doubles(throttles, n,
                                                                          "throttles")
        result = ffi.new("myb_sims_flanagan_result*")
        status = lib.myb_sims_flanagan_solve(self.raw, _nlp_options(feasibility, optimality,
                                                                    max_outer, max_inner, weight),
                                             given, ffi.NULL, result)

        def reached():
            return SimsFlanaganResult(_tuple(given, n), result.mass1, tuple(result.mismatch),
                                      _nlp_result(result.nlp))

        _settle(status, self, partial=reached)
        return reached()


def chain_flyby(gm, v_in, v_out, periapsis):
    """A powered flyby (myb_chain_flyby): (dv, turn, partials), partials d dv / d(v_in,
    v_out, periapsis) then d turn / d(the same)."""
    dv = ffi.new("double*")
    turn = ffi.new("double*")
    partials = ffi.new("double[6]")
    _check(lib.myb_chain_flyby(float(gm), float(v_in), float(v_out), float(periapsis), dv, turn,
                               partials))
    return dv[0], turn[0], tuple(partials)


class Chain:
    """Legs of Sims and Flanagan from body to body joined by flybys (myb_chain). bodies: legs
    + 1 mappings with "state" (a body of an Ephemeris, a Provider or a callable f(seconds) ->
    six numbers, relative to the centre at epoch + t), and for those flown by "gm" and
    "periapsis", optionally "impulse_exhaust" (> 0: a powered flyby). segments: one number a
    leg; departure: (earliest, latest); flight: (shortest, longest) a leg."""

    def __init__(self, mu, epoch, thrust, exhaust, mass0, bodies, segments, departure, flight,
                 launch, arrival=0.0, field=None):
        bodies = list(bodies)
        segments = [int(s) for s in segments]
        flight = [tuple(float(v) for v in f) for f in flight]
        legs = len(segments)
        if len(bodies) != legs + 1 or len(flight) != legs:
            raise ValueError("legs + 1 bodies and one flight a leg")
        if legs > lib.MYB_CHAIN_MAX_LEGS:
            raise ValueError("at most %d legs" % lib.MYB_CHAIN_MAX_LEGS)
        self._calls = _Calls()
        self.raw = ffi.new("myb_chain*")
        raw = self.raw
        raw.mu = float(mu)
        raw.field = ffi.NULL if field is None else field.pointer
        raw.epoch = epoch
        raw.thrust = float(thrust)
        raw.exhaust = float(exhaust)
        raw.mass0 = float(mass0)
        raw.legs = legs
        raw.segments[0:legs] = segments
        self._kept = [field]
        for k, body in enumerate(bodies):
            provider = _state_provider(body["state"], self._calls)
            raw.bodies[k].state = provider.function
            raw.bodies[k].context = provider.context
            raw.bodies[k].gm = float(body.get("gm", 0.0))
            raw.bodies[k].periapsis = float(body.get("periapsis", 0.0))
            raw.impulse_exhaust[k] = float(body.get("impulse_exhaust", 0.0))
            self._kept.append(provider)
        raw.departure = [float(departure[0]), float(departure[1])]
        for k, (low, high) in enumerate(flight):
            raw.flight[k] = [low, high]
        raw.launch = float(launch)
        raw.arrival = float(arrival)

    def size(self):
        """(unknowns, equalities, inequalities) (myb_chain_size)."""
        unknowns = ffi.new("size_t*")
        equalities = ffi.new("size_t*")
        inequalities = ffi.new("size_t*")
        _check(lib.myb_chain_size(self.raw, unknowns, equalities, inequalities))
        return unknowns[0], equalities[0], inequalities[0]

    def guess(self):
        """The first guess built in (myb_chain_guess)."""
        n = self.size()[0]
        x = ffi.new("double[]", n)
        _settle(lib.myb_chain_guess(self.raw, x), self)
        return _tuple(x, n)

    def evaluate(self, x, gradient=False):
        """(propellant, constraints, gradient or None, Jacobian rows or None) at x."""
        n, equalities, inequalities = self.size()
        k = equalities + inequalities
        propellant = ffi.new("double*")
        constraints = ffi.new("double[]", max(k, 1))
        grad = ffi.new("double[]", n) if gradient else ffi.NULL
        jacobian = ffi.new("double[]", max(k * n, 1)) if gradient else ffi.NULL
        _settle(lib.myb_chain_evaluate(self.raw, _doubles(x, n, "x"), propellant, constraints,
                                       grad, jacobian, ffi.NULL), self)
        return (propellant[0], _tuple(constraints, k), _tuple(grad, n) if gradient else None,
                _rows(jacobian, k, n) if gradient else None)

    def solve(self, x=None, feasibility=0.0, optimality=0.0, max_outer=0, max_inner=0,
              weight=0.0):
        """The chain of the least propellant (myb_chain_solve) from x (None: all zeros, the
        guess built in): a ChainResult."""
        n = self.size()[0]
        given = ffi.new("double[]", n) if x is None else _doubles(x, n, "x")
        result = ffi.new("myb_chain_result*")
        status = lib.myb_chain_solve(self.raw, _nlp_options(feasibility, optimality, max_outer,
                                                            max_inner, weight), given, ffi.NULL,
                                     result)
        legs = self.raw.legs

        def reached():
            return ChainResult(_tuple(given, n), _tuple(result.mass, legs + 1),
                               [tuple(result.excess[k]) for k in range(legs + 1)],
                               _tuple(result.turn, legs + 1), _tuple(result.periapsis, legs + 1),
                               _tuple(result.impulse, legs + 1),
                               [tuple(result.mismatch[k]) for k in range(legs)],
                               _nlp_result(result.first), _nlp_result(result.nlp))

        _settle(status, self, partial=reached)
        return reached()


# ---- Thrust ------------------------------------------------------------------------------------

_steering = {"velocity": lib.MYB_STEER_VELOCITY, "anti_velocity": lib.MYB_STEER_ANTI_VELOCITY,
             "inertial": lib.MYB_STEER_INERTIAL, "rtn": lib.MYB_STEER_RTN,
             "edelbaum": lib.MYB_STEER_EDELBAUM}


class Thrust:
    """Burn arcs of constant thrust (myb_thrust): mass0 at the start of the first arc; each
    arc a mapping with start, stop, thrust (N), exhaust (m/s), steering ("velocity",
    "anti_velocity", "inertial", "rtn", "edelbaum"), direction, edelbaum_v0, edelbaum_beta0,
    edelbaum_sense. A term of a model: ForceModel.add_thrust."""

    def __init__(self, mass0, arcs):
        arcs = list(arcs)
        made = ffi.new("myb_thrust_arc[]", max(len(arcs), 1))
        for k, arc in enumerate(arcs):
            for key in ("start", "stop", "thrust", "exhaust", "edelbaum_v0", "edelbaum_beta0",
                        "edelbaum_sense"):
                setattr(made[k], key, float(arc.get(key, 0.0)))
            steering = arc.get("steering", "velocity")
            made[k].steering = _steering[steering] if isinstance(steering, str) else int(steering)
            made[k].direction = list(_doubles(arc.get("direction", (0.0, 0.0, 0.0)), 3,
                                              "direction"))
        self._raw = ffi.new("myb_thrust*")
        _check(lib.myb_thrust_init(self._raw, float(mass0), made, len(arcs)))

    @property
    def pointer(self):
        return self._raw

    def mass(self, t):
        return lib.myb_thrust_mass(self._raw, float(t))

    def mass_at(self, t):
        """The same through the mass function of an object (myb_thrust_mass_at)."""
        return lib.myb_thrust_mass_at(self._raw, float(t))

    def acceleration(self, t, x, v):
        a = ffi.new("double[3]")
        _check(lib.myb_thrust_acceleration(self._raw, float(t), _doubles(x, 3, "x"),
                                           _doubles(v, 3, "v"), a))
        return tuple(a)

    def rtn(self, t, x, v):
        """The acceleration in the radial, transverse, normal frame of (x, v)."""
        a = ffi.new("double[3]")
        _check(lib.myb_thrust_rtn(self._raw, float(t), _doubles(x, 3, "x"), _doubles(v, 3, "v"),
                                  a))
        return tuple(a)


def equinoctial_rates(mu, elements, a_rtn):
    """The Gauss equations in modified equinoctial elements (p, f, g, h, k, L)."""
    given = ffi.new("myb_equinoctial*", dict(zip("pfghkl", (float(v) for v in elements))))
    rates = ffi.new("double[6]")
    lib.myb_equinoctial_rates(float(mu), given, _doubles(a_rtn, 3, "a_rtn"), rates)
    return tuple(rates)


def edelbaum(v0, v1, di):
    """Edelbaum's transfer between circles: (dv, beta0)."""
    dv = ffi.new("double*")
    beta0 = ffi.new("double*")
    _check(lib.myb_edelbaum(float(v0), float(v1), float(di), dv, beta0))
    return dv[0], beta0[0]


def edelbaum_state(v0, beta0, dv):
    """(speed, inclination change) after dv along Edelbaum's law."""
    v = ffi.new("double*")
    di = ffi.new("double*")
    lib.myb_edelbaum_state(float(v0), float(beta0), float(dv), v, di)
    return v[0], di[0]


def _elements(now):
    return ffi.new("myb_elements*", dict(zip(Elements._fields, (float(v) for v in now))))


class QLaw:
    """The Q-law of Petropoulos (myb_qlaw): the target a, e, i, raan, argp and the weights
    of the elements; rp_min the lowest periapsis wanted."""

    def __init__(self, mu, target, weight, rp_min=0.0, rp_weight=0.0):
        self._raw = ffi.new("myb_qlaw*")
        self._raw.mu = float(mu)
        self._raw.target = list(_doubles(target, 5, "target"))
        self._raw.weight = list(_doubles(weight, 5, "weight"))
        self._raw.rp_min = float(rp_min)
        self._raw.rp_weight = float(rp_weight)

    @property
    def pointer(self):
        return self._raw

    def direction(self, now):
        """(direction RTN, Q, its rate, the effectivity) for the orbit and the place given
        by the elements `now` (p, e, i, raan, argp, nu)."""
        direction = ffi.new("double[3]")
        q = ffi.new("double*")
        rate = ffi.new("double*")
        effectivity = ffi.new("double*")
        _check(lib.myb_qlaw_direction(self._raw, _elements(now), direction, q, rate,
                                      effectivity))
        return tuple(direction), q[0], rate[0], effectivity[0]

    def arrived(self, now, tolerance):
        return lib.myb_qlaw_arrived(self._raw, _elements(now),
                                    _doubles(tolerance, 5, "tolerance")) != 0


class QLawTerm:
    """The engine that a Q-law points (myb_qlaw_term): a term of a model
    (ForceModel.add_qlaw) and the engine of a QLawFlight."""

    def __init__(self, law, thrust, exhaust, mass0, t0, tolerance, effectivity=0.0):
        self._raw = ffi.new("myb_qlaw_term*")
        self._raw.law = law.pointer
        self._raw.thrust = float(thrust)
        self._raw.exhaust = float(exhaust)
        self._raw.mass0 = float(mass0)
        self._raw.t0 = float(t0)
        self._raw.tolerance = list(_doubles(tolerance, 5, "tolerance"))
        self._raw.effectivity = float(effectivity)
        self._kept = [law]

    @property
    def pointer(self):
        return self._raw

    def acceleration(self, t, x, v):
        """The term at (t, x, v) (myb_term_qlaw)."""
        a = ffi.new("double[3]")
        _check(lib.myb_term_qlaw(self._raw, float(t), _doubles(x, 3, "x"), _doubles(v, 3, "v"),
                                 a))
        return tuple(a)


class QLawFlight:
    """The engine of a QLawTerm with the mass as a seventh component (myb_qlaw_flight);
    field: a ForceModel of everything but the engine (None: the point mass of the law)."""

    def __init__(self, engine, field=None, rtol=0.0, atol=0.0):
        self._raw = ffi.new("myb_qlaw_flight*")
        self._raw.engine = engine.pointer
        self._raw.field = ffi.NULL if field is None else field.pointer
        self._raw.rtol = float(rtol)
        self._raw.atol = float(atol)
        self._kept = [engine, field]

    def rhs(self, t, y):
        """The right-hand side for (x, v, m) (myb_qlaw_flight_rhs)."""
        dydt = ffi.new("double[7]")
        _settle(lib.myb_qlaw_flight_rhs(self._raw, float(t), _doubles(y, 7, "y"), dydt), self)
        return tuple(dydt)

    def fly(self, t0, y0, t_limit, step):
        """Flies the law from y0 (x, v, m) until the orbit is within the tolerances or
        t_limit (myb_qlaw_fly): a QLawArrival; on a failure Error.partial is the state
        reached."""
        arrival = ffi.new("myb_qlaw_arrival*")
        status = lib.myb_qlaw_fly(self._raw, float(t0), _doubles(y0, 7, "y0"), float(t_limit),
                                  float(step), ffi.NULL, arrival)

        def reached():
            return QLawArrival(arrival.t, tuple(arrival.y), arrival.arrived != 0)

        _settle(status, self, partial=reached)
        return reached()


# ---- Initial orbits ----------------------------------------------------------------------------

def iod_gibbs(mu, r1, r2, r3):
    """The velocity at r2 of the conic through three positions: (v2, out_of_plane)."""
    v2 = ffi.new("double[3]")
    angle = ffi.new("double*")
    _check(lib.myb_iod_gibbs(float(mu), _doubles(r1, 3, "r1"), _doubles(r2, 3, "r2"),
                             _doubles(r3, 3, "r3"), v2, angle))
    return tuple(v2), angle[0]


def iod_herrick_gibbs(mu, r1, t1, r2, t2, r3, t3):
    v2 = ffi.new("double[3]")
    _check(lib.myb_iod_herrick_gibbs(float(mu), _doubles(r1, 3, "r1"), float(t1),
                                     _doubles(r2, 3, "r2"), float(t2), _doubles(r3, 3, "r3"),
                                     float(t3), v2))
    return tuple(v2)


def _iod_of(s):
    return IodSolution(tuple(s.r), tuple(s.v), tuple(s.range), s.iterations)


def iod_gauss(mu, t, observer, sight):
    """Orbits from three directions `sight` seen from `observer` at the times t (Gauss's
    method, refined): a list of IodSolution."""
    times = _doubles(t, 3, "t")
    observers = _doubles(observer, 9, "observer")
    sights = _doubles(sight, 9, "sight")
    out = ffi.new("myb_iod_solution[]", lib.MYB_IOD_MAX_SOLUTIONS)
    count = ffi.new("size_t*")
    _check(lib.myb_iod_gauss(float(mu), times, ffi.cast("double(*)[3]", observers),
                             ffi.cast("double(*)[3]", sights), out, count))
    return [_iod_of(out[k]) for k in range(count[0])]


def iod_laplace(mu, t, observer, sight):
    """The same by Laplace's method."""
    times = _doubles(t, 3, "t")
    observers = _doubles(observer, 9, "observer")
    sights = _doubles(sight, 9, "sight")
    out = ffi.new("myb_iod_solution[]", lib.MYB_IOD_MAX_SOLUTIONS)
    count = ffi.new("size_t*")
    _check(lib.myb_iod_laplace(float(mu), times, ffi.cast("double(*)[3]", observers),
                               ffi.cast("double(*)[3]", sights), out, count))
    return [_iod_of(out[k]) for k in range(count[0])]


def iod_double_r(mu, t, observer, sight, distance):
    """The double-r iteration from guesses of the distances at the first two times."""
    times = _doubles(t, 3, "t")
    observers = _doubles(observer, 9, "observer")
    sights = _doubles(sight, 9, "sight")
    out = ffi.new("myb_iod_solution*")
    _check(lib.myb_iod_double_r(float(mu), times, ffi.cast("double(*)[3]", observers),
                                ffi.cast("double(*)[3]", sights),
                                _doubles(distance, 2, "distance"), out))
    return _iod_of(out)


def iod_gooding(mu, t, observer, sight, range, revolutions=0, branch=0):
    """Gooding's method from guesses of the ranges along the first and the third sight."""
    times = _doubles(t, 3, "t")
    observers = _doubles(observer, 9, "observer")
    sights = _doubles(sight, 9, "sight")
    out = ffi.new("myb_iod_solution*")
    _check(lib.myb_iod_gooding(float(mu), times, ffi.cast("double(*)[3]", observers),
                               ffi.cast("double(*)[3]", sights), _doubles(range, 2, "range"),
                               int(revolutions), int(branch), out))
    return _iod_of(out)


# ---- Relative motion ---------------------------------------------------------------------------

def relative_from_inertial(chief, deputy):
    """The relative state (R, T, N) of the deputy from the inertial states of both."""
    out = ffi.new("double[6]")
    _check(lib.myb_relative_from_inertial(_pv(chief), _pv(deputy), out))
    return tuple(out)


def relative_to_inertial(chief, relative):
    """The inertial state of the deputy."""
    out = ffi.new("myb_pv*")
    _check(lib.myb_relative_to_inertial(_pv(chief), _doubles(relative, 6, "relative"), out))
    return _state_of(out)


def relative_stm(mu, chief, dt):
    """The linear map of the relative state over dt: six rows."""
    phi = ffi.new("double[36]")
    _check(lib.myb_relative_stm(float(mu), _pv(chief), float(dt), phi))
    return _rows(phi, 6, 6)


def cw_stm(n, dt):
    """The map of Clohessy and Wiltshire: six rows."""
    phi = ffi.new("double[36]")
    lib.myb_cw_stm(float(n), float(dt), phi)
    return _rows(phi, 6, 6)


def relative_rendezvous(mu, chief, dt, relative, target, target_velocity=None):
    """The two impulses of a rendezvous in the linear motion: (dv0, dv1)."""
    dv0 = ffi.new("double[3]")
    dv1 = ffi.new("double[3]")
    _check(lib.myb_relative_rendezvous(float(mu), _pv(chief), float(dt),
                                       _doubles(relative, 6, "relative"),
                                       _doubles(target, 3, "target"),
                                       _optional(target_velocity, 3, "target_velocity"), dv0,
                                       dv1))
    return tuple(dv0), tuple(dv1)


def roe_from_states(mu, chief, deputy):
    """The relative orbital elements (da, dlambda, dex, dey, dix, diy)."""
    out = ffi.new("double[6]")
    _check(lib.myb_roe_from_states(float(mu), _pv(chief), _pv(deputy), out))
    return tuple(out)


def roe_to_position(a, u, roe):
    out = ffi.new("double[3]")
    lib.myb_roe_to_position(float(a), float(u), _doubles(roe, 6, "roe"), out)
    return tuple(out)


# ---- Conjunctions ------------------------------------------------------------------------------

def _approach_of(a):
    return Approach(a.t, a.distance, tuple(a.position), tuple(a.velocity))


def approach_find(first, second, t0, t1, step, capacity=64):
    """Every closest approach of two objects given by providers (bodies of an Ephemeris,
    Providers or callables f(seconds) -> six numbers) in [t0, t1] (myb_approach_find): a
    list of Approach."""
    calls = _Calls()
    one = _state_provider(first, calls)
    two = _state_provider(second, calls)
    count = ffi.new("size_t*")
    while True:
        out = ffi.new("myb_approach[]", max(int(capacity), 1))
        _settle(lib.myb_approach_find(one.function, one.context, two.function, two.context,
                                      float(t0), float(t1), float(step), out, int(capacity),
                                      count), calls, one, two)
        if count[0] <= capacity:
            return [_approach_of(out[k]) for k in range(count[0])]
        capacity = count[0]


def approach_orbits(first, second, t1, step, capacity=1024):
    """The closest approaches of two Orbits taken together to t1 (myb_approach_orbits): a
    list of at most `capacity` Approach (the orbits move: a larger capacity is the way to
    more)."""
    out = ffi.new("myb_approach[]", int(capacity))
    count = ffi.new("size_t*")
    _settle(lib.myb_approach_orbits(first.pointer, second.pointer, float(t1), float(step), out,
                                    int(capacity), count), first, second)
    return [_approach_of(out[k]) for k in range(min(count[0], int(capacity)))]


def _encounter_in(encounter):
    return ffi.new("myb_encounter*", {"miss": list(encounter.miss),
                                      "sigma": list(encounter.sigma),
                                      "radius": float(encounter.radius)})


def _encounter_of(e):
    return Encounter(tuple(e.miss), tuple(e.sigma), e.radius)


def encounter_plane(position, velocity, covariance, radius):
    """The encounter in the plane across the relative velocity (myb_encounter_plane)."""
    out = ffi.new("myb_encounter*")
    _check(lib.myb_encounter_plane(_doubles(position, 3, "position"),
                                   _doubles(velocity, 3, "velocity"),
                                   _doubles(covariance, 9, "covariance"), float(radius), out))
    return _encounter_of(out)


def collision_probability(encounter):
    out = ffi.new("double*")
    _check(lib.myb_collision_probability(_encounter_in(encounter), out))
    return out[0]


def collision_probability_chan(encounter):
    return lib.myb_collision_probability_chan(_encounter_in(encounter))


def collision_probability_max(encounter):
    """(the largest probability, the factor of the sigmas that gives it)."""
    probability = ffi.new("double*")
    factor = ffi.new("double*")
    _check(lib.myb_collision_probability_max(_encounter_in(encounter), probability, factor))
    return probability[0], factor[0]


def conjunction_find(first, first_covariance, second, second_covariance, radius, t1, step,
                     capacity=1024):
    """The approaches of two Orbits created with variational=True, with their covariances
    carried there (myb_conjunction_find): a list of Conjunction."""
    out = ffi.new("myb_conjunction[]", int(capacity))
    count = ffi.new("size_t*")
    _settle(lib.myb_conjunction_find(first.pointer, _doubles(first_covariance, 36, "covariance"),
                                     second.pointer,
                                     _doubles(second_covariance, 36, "covariance"),
                                     float(radius), float(t1), float(step), out, int(capacity),
                                     count), first, second)
    return [Conjunction(_approach_of(c.approach), _rows(c.covariance, 3, 3),
                        _encounter_of(c.encounter), c.probability)
            for c in (out[k] for k in range(min(count[0], int(capacity))))]


def collision_rate(mean, covariance, radius):
    out = ffi.new("double*")
    _check(lib.myb_collision_rate(_doubles(mean, 6, "mean"), _doubles(covariance, 36,
                                                                      "covariance"),
                                  float(radius), out))
    return out[0]


def collision_probability_long(function, t0, t1, step, radius):
    """The expected number of entries over [t0, t1] (myb_collision_probability_long);
    function(t) -> (mean six numbers, covariance 6 x 6). Returns (probability, error)."""
    calls = _Calls()

    @ffi.callback("myb_status(void*, double, double[6], double[36])")
    def gaussian(context, t, mean, covariance):
        try:
            m, p = function(t)
            mean[0:6] = list(_doubles(m, 6, "mean"))
            covariance[0:36] = list(_doubles(p, 36, "covariance"))
            return lib.MYB_OK
        except BaseException as error:  # noqa: B902 - nothing may cross the C code
            return calls.failed(error)

    probability = ffi.new("double*")
    error = ffi.new("double*")
    _settle(lib.myb_collision_probability_long(gaussian, ffi.NULL, float(t0), float(t1),
                                               float(step), float(radius), probability, error),
            calls)
    return probability[0], error[0]


def conjunction_long(first, first_covariance, second, second_covariance, radius, t1, step):
    """The same for two Orbits with variational=True (myb_conjunction_long):
    (probability, error)."""
    probability = ffi.new("double*")
    error = ffi.new("double*")
    _settle(lib.myb_conjunction_long(first.pointer, _doubles(first_covariance, 36, "covariance"),
                                     second.pointer,
                                     _doubles(second_covariance, 36, "covariance"),
                                     float(radius), float(t1), float(step), probability, error),
            first, second)
    return probability[0], error[0]


# ---- Uncertainty -------------------------------------------------------------------------------

class Variational(_Owner):
    """The variational equations of a force model (myb_variational): the state with its
    transition matrix and its sensitivities to the parameters of the model."""

    _destroy = staticmethod(lib.myb_variational_destroy)

    def __init__(self, model, epsilon=0.0, h_max=0.0, max_steps=0):
        options = ffi.new("myb_variational_options*", {"epsilon": float(epsilon),
                                                       "h_max": float(h_max),
                                                       "max_steps": int(max_steps)})
        out = ffi.new("myb_variational**")
        _check(lib.myb_variational_create(model.pointer, options, ffi.NULL, out))
        super().__init__(out[0])
        self._model = model
        self._kept.append(model)
        self._np = lib.myb_force_model_parameter_count(model.pointer)

    def start(self, t0, y0):
        _settle(lib.myb_variational_start(self.pointer, float(t0), _doubles(y0, 6, "y0")),
                self._model)

    def step(self, t_limit):
        _settle(lib.myb_variational_step(self.pointer, float(t_limit)), self._model)

    def integrate(self, t_end):
        _settle(lib.myb_variational_integrate(self.pointer, float(t_end)), self._model)

    @property
    def time(self):
        return lib.myb_variational_time(self.pointer)

    def _outputs(self):
        return (ffi.new("double[6]"), ffi.new("double[36]"),
                ffi.new("double[]", max(6 * self._np, 1)))

    def get(self):
        """(y, phi six rows, sensitivity six rows of np) now."""
        y, phi, s = self._outputs()
        lib.myb_variational_get(self.pointer, y, phi, s)
        return tuple(y), _rows(phi, 6, 6), _rows(s, 6, self._np)

    def dense(self, t):
        y, phi, s = self._outputs()
        _check(lib.myb_variational_dense(self.pointer, float(t), y, phi, s))
        return tuple(y), _rows(phi, 6, 6), _rows(s, 6, self._np)


def covariance_transform(m, p):
    """M P M^T for square matrices (rows or flat)."""
    mm, n = _square(m, "m")
    pp, n2 = _square(p, "p")
    if n != n2:
        raise ValueError("m and p: matrices of one size")
    out = ffi.new("double[]", n * n)
    lib.myb_covariance_transform(n, mm, pp, out)
    return _rows(out, n, n)


def covariance_rtn(r, v, p):
    """A 6 x 6 covariance in the radial, transverse, normal axes of (r, v)."""
    out = ffi.new("double[36]")
    _check(lib.myb_covariance_rtn(_doubles(r, 3, "r"), _doubles(v, 3, "v"),
                                  _doubles(p, 36, "p"), out))
    return _rows(out, 6, 6)


class Rng:
    """The deterministic random numbers of the library (myb_rng: xoshiro256**)."""

    def __init__(self, seed):
        self._raw = ffi.new("myb_rng*")
        lib.myb_rng_seed(self._raw, int(seed))

    def next(self):
        """64 random bits as an int."""
        return int(lib.myb_rng_next(self._raw))

    def uniform(self):
        return lib.myb_rng_uniform(self._raw)

    def normal(self):
        return lib.myb_rng_normal(self._raw)

    def gaussian(self, mean, l):
        """A sample of N(mean, L L^T), L lower triangular (rows or flat)."""
        mean = _flat(mean)
        n = len(mean)
        out = ffi.new("double[]", max(n, 1))
        lib.myb_rng_gaussian(self._raw, n, ffi.new("double[]", mean or 1),
                             _doubles(l, n * n, "l"), out)
        return _tuple(out, n)


def sample_statistics(samples):
    """(mean, covariance rows) of samples of one dimension (myb_sample_statistics)."""
    samples = [_flat(s) for s in samples]
    n = len(samples[0]) if samples else 0
    if any(len(s) != n for s in samples):
        raise ValueError("samples of one dimension")
    flat = ffi.new("double[]", [v for s in samples for v in s] or 1)
    mean = ffi.new("double[]", max(n, 1))
    covariance = ffi.new("double[]", max(n * n, 1))
    _check(lib.myb_sample_statistics(len(samples), n, flat, mean, covariance))
    return _tuple(mean, n), _rows(covariance, n, n)
