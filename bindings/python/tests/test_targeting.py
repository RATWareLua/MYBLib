# SPDX-License-Identifier: BSD-3-Clause
"""Targeting through the package against the same C calls through cffi, bit for bit:
Lambert, the B-plane, shooting and the solver of nonlinear programs with Python functions,
the indirect method, a leg of Sims and Flanagan, a chain of legs and a powered flyby, thrust
arcs and the Q-law, initial orbits, relative motion, approaches and conjunctions of two
orbits, the variational equations, covariances and random numbers; failures as myb.Error."""
import math

import myb
from myb import ffi, lib

from common import ELEMENTS, GM, Case, doubles, ok, rows

DAY = 86400.0


def circle(radius, angle, tilt=0.0):
    v = math.sqrt(GM / radius)
    return (radius * math.cos(angle), radius * math.sin(angle) * math.cos(tilt),
            radius * math.sin(angle) * math.sin(tilt), -v * math.sin(angle),
            v * math.cos(angle) * math.cos(tilt), v * math.cos(angle) * math.sin(tilt))


class Transfers(Case):
    def test_lambert_bplane_and_patched_conics(self):
        r1, r2 = (7.0e6, 0.0, 0.0), (0.0, 8.0e6, 1.0e5)
        solutions = myb.lambert(GM, r1, r2, 3.0 * 3600.0, max_revolutions=1)
        out = ffi.new("myb_lambert_solution[3]")
        count = ffi.new("size_t*")
        ok(lib.myb_lambert(GM, ffi.new("double[3]", r1), ffi.new("double[3]", r2), 3.0 * 3600.0,
                           0, 1, out, 3, count))
        self.assertSame(solutions, [myb.LambertSolution(tuple(s.v1), tuple(s.v2), s.revolutions,
                                                        s.branch)
                                    for s in (out[k] for k in range(count[0]))])
        self.assertStatus(lib.MYB_E_DOMAIN, myb.lambert, GM, r1, r2, -1.0)
        r, v = (-1.0e8, 2.0e6, 0.0), (6.0e3, 0.0, 100.0)
        plane = myb.bplane(4.28e13, r, v)
        made = ffi.new("myb_bplane*")
        ok(lib.myb_bplane_compute(4.28e13, ffi.new("double[3]", r), ffi.new("double[3]", v),
                                  ffi.new("double[3]", [0.0, 0.0, 1.0]), made))
        self.assertSame((plane.b_t, plane.b_r, plane.b, plane.time_to_periapsis),
                        (made.b_t, made.b_r, tuple(made.b), made.time_to_periapsis))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.bplane, 4.28e13, r, (1.0, 0.0, 0.0))
        self.assertSame(myb.flyby_turn(4.28e13, 5e3, 4.0e6),
                        lib.myb_flyby_turn(4.28e13, 5e3, 4.0e6))
        self.assertSame((myb.sphere_of_influence(1.5e11, 3e-6), myb.hill_radius(1.5e11, 0.1, 3e-6),
                         myb.departure_impulse(GM, 7e6, 3e3)),
                        (lib.myb_sphere_of_influence(1.5e11, 3e-6),
                         lib.myb_hill_radius(1.5e11, 0.1, 3e-6),
                         lib.myb_departure_impulse(GM, 7e6, 3e3)))

    def test_shooting_and_a_nonlinear_program(self):
        def f(u, want):
            return ((u[0] * u[0] - 2.0, u[1] - u[0] * u[0] * u[0]),
                    (2.0 * u[0], 0.0, -3.0 * u[0] * u[0], 1.0) if want else None)

        u, iterations, left = myb.shoot(f, (1.0, 1.0), 2, 1e-14, jacobian=True)

        @ffi.callback("myb_status(void*, const double*, double*, double*)")
        def raw_f(context, x, out, jacobian):
            values, matrix = f((x[0], x[1]), jacobian != ffi.NULL)
            out[0:2] = list(values)
            if jacobian != ffi.NULL:
                jacobian[0:4] = list(matrix)
            return lib.MYB_OK

        x = ffi.new("double[2]", [1.0, 1.0])
        work = ffi.new("double[]", lib.myb_shoot_work_size(2, 2))
        n = ffi.new("int*")
        r = ffi.new("double*")
        ok(lib.myb_shoot(raw_f, ffi.NULL, 2, 2, x, ffi.new("myb_shoot_options*", [1e-14, 0]), work,
                         n, r))
        self.assertSame((u, iterations, left), (doubles(x, 2), n[0], r[0]))
        self.assertLess(abs(u[0] - math.sqrt(2.0)), 1e-14)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.shoot, lambda u: u, (1.0,), 1, -1.0)

        def program(x, want):
            value = (x[0] - 1.0) ** 2 + (x[1] - 2.0) ** 2
            constraints = (x[0] + x[1] - 2.0,)
            if not want:
                return value, constraints
            return value, constraints, (2.0 * (x[0] - 1.0), 2.0 * (x[1] - 2.0)), (1.0, 1.0)

        point, multipliers, result = myb.nlp_solve(program, (0.0, 0.0), equalities=1,
                                                   lower=(-5.0, -5.0), upper=(5.0, 5.0))

        @ffi.callback("myb_status(void*, const double*, double*, double*, double*, double*)")
        def raw_program(context, x, value, c, gradient, jacobian):
            answer = program((x[0], x[1]), gradient != ffi.NULL)
            value[0] = answer[0]
            c[0] = answer[1][0]
            if gradient != ffi.NULL:
                gradient[0:2] = list(answer[2])
                jacobian[0:2] = list(answer[3])
            return lib.MYB_OK

        problem = ffi.new("myb_nlp_problem*", {"n": 2, "equalities": 1, "inequalities": 0})
        low = ffi.new("double[2]", [-5.0, -5.0])
        high = ffi.new("double[2]", [5.0, 5.0])
        problem.lower = low
        problem.upper = high
        problem.function = raw_program
        x = ffi.new("double[2]")
        lam = ffi.new("double[1]")
        made = ffi.new("myb_nlp_result*")
        ok(lib.myb_nlp_solve(problem, ffi.NULL, x, lam, ffi.NULL, made))
        self.assertSame((point, multipliers, result.objective, result.evaluations),
                        (doubles(x, 2), doubles(lam, 1), made.objective, made.evaluations))
        self.assertLess(abs(point[0] - 0.5), 1e-6)

        def wrong(x, want):
            raise ZeroDivisionError("nothing to minimise")

        with self.assertRaises(ZeroDivisionError):
            myb.nlp_solve(wrong, (0.0,))


class Optimal(Case):
    def test_indirect_least_time(self):
        start, target = circle(7.0e6, 0.0), circle(8.0e6, 0.5)
        problem = myb.IndirectProblem(GM, 300.0, 3000.0, 0.0, start, 1000.0, 2500.0, target,
                                      goal="time")
        solution = problem.solve()
        raw = ffi.new("myb_indirect_problem*", {"mu": GM, "thrust": 300.0, "exhaust": 3000.0,
                                                "t0": 0.0, "start": start, "mass0": 1000.0,
                                                "t1": 2500.0, "target": target,
                                                "goal": lib.MYB_INDIRECT_TIME})
        made = ffi.new("myb_indirect_solution*")
        ok(lib.myb_indirect_solve(raw, ffi.NULL, ffi.NULL, made))
        self.assertSame(solution, myb.IndirectSolution(tuple(made.costate), made.t1, made.mass1,
                                                       made.smoothing, made.residual, made.tries,
                                                       made.stages, made.phase))
        times = [solution.t1 * k / 4 for k in range(5)]
        points = problem.sample(solution, times)
        out = ffi.new("myb_indirect_point[5]")
        ok(lib.myb_indirect_sample(raw, made, ffi.new("double[]", times), 5, ffi.NULL, out))
        self.assertSame([p.y for p in points], [tuple(out[k].y) for k in range(5)])
        self.assertSame(points[2].throttle, out[2].throttle)
        bad = myb.IndirectProblem(GM, 0.0, 3000.0, 0.0, start, 1000.0, 2500.0, target,
                                  goal="time")
        self.assertStatus(lib.MYB_E_DOMAIN, bad.solve)

    def test_sims_flanagan_and_chains(self):
        start, target = circle(7.0e6, 0.0), circle(8.0e6, 2.5)
        leg = myb.SimsFlanagan(GM, 300.0, 3000.0, 0.0, start, 1000.0, 4500.0, target, 6)
        throttles = [0.3 * ((k % 3) - 1) for k in range(18)]
        mismatch, mass1, jacobian = leg.mismatch(throttles, jacobian=True)
        raw = ffi.new("myb_sims_flanagan*", {"mu": GM, "thrust": 300.0, "exhaust": 3000.0,
                                             "t0": 0.0, "start": start, "mass0": 1000.0,
                                             "t1": 4500.0, "target": target, "segments": 6})
        out = ffi.new("double[6]")
        m = ffi.new("double*")
        rows_out = ffi.new("double[108]")
        ok(lib.myb_sims_flanagan_mismatch(raw, ffi.new("double[]", throttles), out, m, rows_out))
        self.assertSame((mismatch, mass1, jacobian), (doubles(out, 6), m[0],
                                                      rows(rows_out, 6, 18)))
        self.assertStatus(lib.MYB_E_DOMAIN, leg.mismatch, [2.0] * 18)
        # a powered flyby
        dv, turn, partials = myb.chain_flyby(3.2e14, 3.0e3, 3.5e3, 7.0e6)
        a = ffi.new("double*")
        b = ffi.new("double*")
        p = ffi.new("double[6]")
        ok(lib.myb_chain_flyby(3.2e14, 3.0e3, 3.5e3, 7.0e6, a, b, p))
        self.assertSame((dv, turn, partials), (a[0], b[0], doubles(p, 6)))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.chain_flyby, 3.2e14, -1.0, 3.5e3, 7.0e6)
        # a chain of one leg between two bodies on circles about the Sun
        sun = 1.32712440018e20

        def planet(radius, phase):
            rate = math.sqrt(sun / radius ** 3)

            def state(seconds):
                angle = phase + rate * seconds
                return (radius * math.cos(angle), radius * math.sin(angle), 0.0,
                        -radius * rate * math.sin(angle), radius * rate * math.cos(angle), 0.0)
            return state

        earth, mars = planet(1.496e11, 0.0), planet(2.279e11, 0.8)
        chain = myb.Chain(sun, myb.time_from_seconds(0.0), 1.0, 30000.0, 1000.0,
                          [{"state": earth}, {"state": mars}], [5], (0.0, 20 * DAY),
                          [(200 * DAY, 260 * DAY)], 2500.0)
        size = chain.size()
        guess = chain.guess()
        propellant, constraints, gradient, jacobian = chain.evaluate(guess, gradient=True)

        def callback(function):
            @ffi.callback("myb_status(void*, myb_time, double[6])")
            def provider(context, instant, y):
                y[0:6] = list(function(lib.myb_time_seconds(instant)))
                return lib.MYB_OK
            return provider

        providers = [callback(earth), callback(mars)]
        made = ffi.new("myb_chain*")
        made.mu = sun
        made.epoch = myb.time_from_seconds(0.0)
        made.thrust = 1.0
        made.exhaust = 30000.0
        made.mass0 = 1000.0
        made.legs = 1
        made.segments[0] = 5
        made.bodies[0].state = providers[0]
        made.bodies[1].state = providers[1]
        made.departure = [0.0, 20 * DAY]
        made.flight[0] = [200 * DAY, 260 * DAY]
        made.launch = 2500.0
        counts = ffi.new("size_t[3]")
        ok(lib.myb_chain_size(made, counts, counts + 1, counts + 2))
        self.assertEqual(size, tuple(counts))
        x = ffi.new("double[]", size[0])
        ok(lib.myb_chain_guess(made, x))
        self.assertSame(guess, doubles(x, size[0]))
        k = size[1] + size[2]
        f = ffi.new("double*")
        c = ffi.new("double[]", k)
        g = ffi.new("double[]", size[0])
        jac = ffi.new("double[]", k * size[0])
        ok(lib.myb_chain_evaluate(made, x, f, c, g, jac, ffi.NULL))
        self.assertSame((propellant, constraints, gradient, jacobian),
                        (f[0], doubles(c, k), doubles(g, size[0]), rows(jac, k, size[0])))

        def broken(seconds):
            raise RuntimeError("the planet is gone")

        lost = myb.Chain(sun, myb.time_from_seconds(0.0), 1.0, 30000.0, 1000.0,
                         [{"state": earth}, {"state": broken}], [5], (0.0, 20 * DAY),
                         [(200 * DAY, 260 * DAY)], 2500.0)
        with self.assertRaises(RuntimeError):
            lost.guess()


class Thrusting(Case):
    def test_arcs_edelbaum_and_the_qlaw(self):
        arcs = [{"start": 0.0, "stop": 600.0, "thrust": 1.0, "exhaust": 30000.0},
                {"start": 1200.0, "stop": 1800.0, "thrust": 0.5, "exhaust": 20000.0,
                 "steering": "rtn", "direction": (0.0, 1.0, 0.0)}]
        thrust = myb.Thrust(500.0, arcs)
        made = ffi.new("myb_thrust_arc[2]")
        made[0] = {"start": 0.0, "stop": 600.0, "thrust": 1.0, "exhaust": 30000.0}
        made[1] = {"start": 1200.0, "stop": 1800.0, "thrust": 0.5, "exhaust": 20000.0,
                   "steering": lib.MYB_STEER_RTN, "direction": [0.0, 1.0, 0.0]}
        raw = ffi.new("myb_thrust*")
        ok(lib.myb_thrust_init(raw, 500.0, made, 2))
        y = myb.state_from_elements(GM, ELEMENTS)
        a = ffi.new("double[3]")
        ok(lib.myb_thrust_acceleration(raw, 1500.0, ffi.new("double[3]", y[:3]),
                                       ffi.new("double[3]", y[3:]), a))
        self.assertSame(thrust.acceleration(1500.0, y[:3], y[3:]), doubles(a, 3))
        self.assertSame(thrust.mass(1500.0), lib.myb_thrust_mass(raw, 1500.0))
        model = myb.ForceModel()
        model.add_thrust("engine", thrust)
        self.assertEqual(model.acceleration(1500.0, y[:3], y[3:]), doubles(a, 3))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.Thrust, -1.0, arcs)
        dv, beta0 = myb.edelbaum(7.5e3, 3.07e3, 0.5)
        b = ffi.new("double*")
        c = ffi.new("double*")
        ok(lib.myb_edelbaum(7.5e3, 3.07e3, 0.5, b, c))
        self.assertSame((dv, beta0), (b[0], c[0]))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.edelbaum, -1.0, 3.0e3, 0.0)
        law = myb.QLaw(GM, (1.0e7, 0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 0.0, 0.0, 0.0))
        low = myb.Elements(7.0e6, 0.0, 0.1, 0.0, 0.0, 0.3)
        direction = law.direction(low)
        q = ffi.new("myb_qlaw*", {"mu": GM, "target": [1.0e7, 0.0, 0.0, 0.0, 0.0],
                                  "weight": [1.0, 1.0, 0.0, 0.0, 0.0]})
        d = ffi.new("double[3]")
        values = ffi.new("double[3]")
        ok(lib.myb_qlaw_direction(q, ffi.new("myb_elements*", low._asdict()), d, values,
                                  values + 1, values + 2))
        self.assertSame(direction, (doubles(d, 3),) + doubles(values, 3))
        engine = myb.QLawTerm(law, 1.0, 30000.0, 1000.0, 0.0, (2.0e3, 1e-2, 2e-4, 1e-3, 1e-3))
        flight = myb.QLawFlight(engine)
        start = myb.state_from_elements(GM, low) + (1000.0,)
        arrival = flight.fly(0.0, start, 0.5 * DAY, 300.0)
        term = ffi.new("myb_qlaw_term*", {"law": q, "thrust": 1.0, "exhaust": 30000.0,
                                          "mass0": 1000.0,
                                          "tolerance": [2.0e3, 1e-2, 2e-4, 1e-3, 1e-3]})
        raw_flight = ffi.new("myb_qlaw_flight*", {"engine": term})
        out = ffi.new("myb_qlaw_arrival*")
        ok(lib.myb_qlaw_fly(raw_flight, 0.0, ffi.new("double[7]", start), 0.5 * DAY, 300.0,
                            ffi.NULL, out))
        self.assertSame(arrival, myb.QLawArrival(out.t, tuple(out.y), out.arrived != 0))
        self.assertFalse(arrival.arrived)
        self.assertEqual(len(flight.rhs(0.0, start)), 7)
        self.assertEqual(len(engine.acceleration(0.0, start[:3], start[3:6])), 3)
        self.assertStatus(lib.MYB_E_DOMAIN, flight.fly, 0.0, start, 0.0, 300.0)


class InitialOrbits(Case):
    def test_gibbs_and_gauss(self):
        states = myb.kepler_propagate_many(GM, myb.state_from_elements(GM, ELEMENTS),
                                           [0.0, 300.0, 600.0])
        r = [s[:3] for s in states]
        v2, angle = myb.iod_gibbs(GM, *r)
        out = ffi.new("double[3]")
        a = ffi.new("double*")
        ok(lib.myb_iod_gibbs(GM, ffi.new("double[3]", r[0]), ffi.new("double[3]", r[1]),
                             ffi.new("double[3]", r[2]), out, a))
        self.assertSame((v2, angle), (doubles(out, 3), a[0]))
        self.assertSame(myb.iod_herrick_gibbs(GM, r[0], 0.0, r[1], 300.0, r[2], 600.0),
                        myb.iod_herrick_gibbs(GM, r[0], 0.0, r[1], 300.0, r[2], 600.0))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.iod_gibbs, GM, r[0], r[0], r[0])
        observer = [(6.0e6, 0.0, 1.0e6), (5.9e6, 1.0e6, 1.0e6), (5.6e6, 2.0e6, 1.0e6)]
        sight = []
        for o, s in zip(observer, states):
            d = [a - b for a, b in zip(s[:3], o)]
            n = math.sqrt(sum(c * c for c in d))
            sight.append([c / n for c in d])
        found = myb.iod_gauss(GM, (0.0, 300.0, 600.0), observer, sight)
        made = ffi.new("myb_iod_solution[3]")
        count = ffi.new("size_t*")
        ok(lib.myb_iod_gauss(GM, ffi.new("double[3]", [0.0, 300.0, 600.0]),
                             ffi.new("double[3][3]", observer), ffi.new("double[3][3]", sight),
                             made, count))
        self.assertSame(found, [myb.IodSolution(tuple(s.r), tuple(s.v), tuple(s.range),
                                                s.iterations)
                                for s in (made[k] for k in range(count[0]))])
        self.assertLess(max(abs(a - b) for a, b in zip(found[0].v, states[1][3:])), 1e-3)
        self.assertEqual(len(myb.iod_laplace(GM, (0.0, 300.0, 600.0), observer, sight)[0].r), 3)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.iod_gauss, GM, (0.0, 0.0, 600.0), observer,
                          sight)


class Relative(Case):
    def test_relative_motion(self):
        chief = myb.state_from_elements(GM, ELEMENTS)
        deputy = myb.state_from_elements(GM, ELEMENTS._replace(nu=ELEMENTS.nu + 1e-4))
        relative = myb.relative_from_inertial(chief, deputy)
        out = ffi.new("double[6]")
        pv = lambda s: ffi.new("myb_pv*", {"r": s[:3], "v": s[3:]})  # noqa: E731
        ok(lib.myb_relative_from_inertial(pv(chief), pv(deputy), out))
        self.assertSame(relative, doubles(out, 6))
        phi = ffi.new("double[36]")
        ok(lib.myb_relative_stm(GM, pv(chief), 600.0, phi))
        self.assertSame(myb.relative_stm(GM, chief, 600.0), rows(phi, 6, 6))
        lib.myb_cw_stm(1e-3, 600.0, phi)
        self.assertSame(myb.cw_stm(1e-3, 600.0), rows(phi, 6, 6))
        dv0, dv1 = myb.relative_rendezvous(GM, chief, 1200.0, relative, (0.0, 0.0, 0.0))
        self.assertEqual(len(dv0) + len(dv1), 6)
        roe = myb.roe_from_states(GM, chief, deputy)
        ok(lib.myb_roe_from_states(GM, pv(chief), pv(deputy), out))
        self.assertSame(roe, doubles(out, 6))
        self.assertEqual(len(myb.roe_to_position(7.0e6, 0.3, roe)), 3)
        back = myb.relative_to_inertial(chief, relative)
        self.assertLess(max(abs(a - b) for a, b in zip(back, deputy)), 1e-6)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.relative_from_inertial, (0.0,) * 6, deputy)


class Conjunctions(Case):
    def orbits(self, variational=False):
        model = myb.ForceModel()
        model.add_point_mass("earth", GM)
        one = myb.Orbit(model, variational=variational)
        two = myb.Orbit(model, variational=variational)
        a = myb.state_from_elements(GM, ELEMENTS)
        b = myb.state_from_elements(GM, ELEMENTS._replace(i=ELEMENTS.i + 0.01, nu=ELEMENTS.nu))
        one.start(0.0, a, 7200.0)
        two.start(0.0, b, 7200.0)
        return one, two

    def test_approaches_and_probabilities(self):
        one, two = self.orbits()
        found = myb.approach_orbits(one, two, 7200.0, 10.0)
        three, four = self.orbits()
        out = ffi.new("myb_approach[64]")
        count = ffi.new("size_t*")
        ok(lib.myb_approach_orbits(three.pointer, four.pointer, 7200.0, 10.0, out, 64, count))
        self.assertSame(found, [myb.Approach(a.t, a.distance, tuple(a.position),
                                             tuple(a.velocity))
                                for a in (out[k] for k in range(count[0]))])
        self.assertGreater(len(found), 0)
        covariance = [[100.0 if r == c and r < 3 else 1e-4 if r == c else 0.0 for c in range(6)]
                      for r in range(6)]
        one, two = self.orbits(variational=True)
        conjunctions = myb.conjunction_find(one, covariance, two, covariance, 10.0, 7200.0, 10.0)
        three, four = self.orbits(variational=True)
        made = ffi.new("myb_conjunction[64]")
        flat = ffi.new("double[36]", [v for row in covariance for v in row])
        ok(lib.myb_conjunction_find(three.pointer, flat, four.pointer, flat, 10.0, 7200.0, 10.0,
                                    made, 64, count))
        self.assertSame([c.probability for c in conjunctions],
                        [made[k].probability for k in range(count[0])])
        self.assertSame(conjunctions[0].encounter.sigma, tuple(made[0].encounter.sigma))
        e = conjunctions[0].encounter
        p = ffi.new("double*")
        ok(lib.myb_collision_probability(ffi.addressof(made[0], "encounter"), p))
        self.assertSame(myb.collision_probability(e), p[0])
        self.assertSame(myb.collision_probability_chan(e),
                        lib.myb_collision_probability_chan(ffi.addressof(made[0], "encounter")))
        self.assertEqual(len(myb.collision_probability_max(e)), 2)
        plane = myb.encounter_plane((10.0, 0.0, 0.0), (0.0, 1.0e4, 0.0),
                                    (4.0, 0, 0, 0, 4.0, 0, 0, 0, 4.0), 5.0)
        self.assertSame(plane.sigma, (2.0, 2.0))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.encounter_plane, (10.0, 0.0, 0.0),
                          (0.0, 0.0, 0.0), (4.0, 0, 0, 0, 4.0, 0, 0, 0, 4.0), 5.0)
        # a slow encounter: the uncertainty of the velocity half the speed
        wide = [[100.0 if r == c and r < 3 else 0.25 if r == c else 0.0 for c in range(6)]
                for r in range(6)]
        mean = (10.0, 0.0, 0.0, 0.0, 1.0, 0.0)
        rate = myb.collision_rate(mean, wide, 5.0)
        ok(lib.myb_collision_rate(ffi.new("double[6]", mean),
                                  ffi.new("double[36]", [v for row in wide for v in row]), 5.0,
                                  p))
        self.assertSame(rate, p[0])
        probability, error = myb.collision_probability_long(
            lambda t: ((10.0, -100.0 + t, 0.0, 0.0, 1.0, 0.0), wide), 0.0, 200.0, 20.0, 5.0)
        self.assertTrue(0.0 < probability < 1.0)

        def broken(t):
            raise ValueError("lost track")

        with self.assertRaises(ValueError):
            myb.collision_probability_long(broken, 0.0, 200.0, 20.0, 5.0)
        one, two = self.orbits(variational=True)
        self.assertEqual(len(myb.conjunction_long(one, covariance, two, covariance, 10.0, 600.0,
                                                  60.0)), 2)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.approach_orbits, one, two, 600.0, 10.0)


class Uncertainty(Case):
    def test_variational_covariances_and_random_numbers(self):
        model = myb.ForceModel()
        model.add_point_mass("earth", GM)
        y0 = myb.state_from_elements(GM, ELEMENTS)
        variational = myb.Variational(model)
        variational.start(0.0, y0)
        variational.integrate(3000.0)
        now = variational.get()
        raw_model = ffi.new("myb_force_model**")
        ok(lib.myb_force_model_create(ffi.NULL, raw_model))
        term = ffi.new("myb_point_mass_term*", [GM])
        ok(lib.myb_force_model_add(raw_model[0], b"earth", lib.myb_term_point_mass,
                                   lib.myb_term_point_mass_jacobian, term, ffi.NULL))
        raw = ffi.new("myb_variational**")
        ok(lib.myb_variational_create(raw_model[0], ffi.NULL, ffi.NULL, raw))
        ok(lib.myb_variational_start(raw[0], 0.0, ffi.new("double[6]", y0)))
        ok(lib.myb_variational_integrate(raw[0], 3000.0))
        y = ffi.new("double[6]")
        phi = ffi.new("double[36]")
        lib.myb_variational_get(raw[0], y, phi, ffi.NULL)
        self.assertSame(now[:2], (doubles(y, 6), rows(phi, 6, 6)))
        lib.myb_variational_destroy(raw[0])
        lib.myb_force_model_destroy(raw_model[0])
        self.assertStatus(lib.MYB_E_RANGE, variational.dense, -1.0e6)
        p = [[float(r == c) * (r + 1) for c in range(6)] for r in range(6)]
        transformed = myb.covariance_transform(now[1], p)
        out = ffi.new("double[36]")
        lib.myb_covariance_transform(6, phi, ffi.new("double[36]", [v for row in p for v in row]),
                                     out)
        self.assertSame(transformed, rows(out, 6, 6))
        self.assertEqual(len(myb.covariance_rtn(y0[:3], y0[3:], p)), 6)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.covariance_rtn, (1.0, 0, 0), (2.0, 0, 0), p)
        rng = myb.Rng(42)
        raw_rng = ffi.new("myb_rng*")
        lib.myb_rng_seed(raw_rng, 42)
        self.assertEqual([rng.next(), rng.uniform(), rng.normal()],
                         [lib.myb_rng_next(raw_rng), lib.myb_rng_uniform(raw_rng),
                          lib.myb_rng_normal(raw_rng)])
        samples = [rng.gaussian((0.0, 1.0), (1.0, 0.0, 0.5, 1.0)) for _ in range(200)]
        mean, covariance = myb.sample_statistics(samples)
        self.assertLess(abs(mean[1] - 1.0), 0.3)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.sample_statistics, samples[:1])
