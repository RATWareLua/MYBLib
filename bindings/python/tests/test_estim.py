# SPDX-License-Identifier: BSD-3-Clause
"""Estimation through the classes against the same C calls through cffi, bit for bit: the
measurements of a station with their partials, a fit of an orbit to positions, the filter and
its smoother; a parameter of the model estimated with the orbit; measurements in Python;
failures as myb.Error with the status."""
import myb
from myb import ffi, lib

from common import ELEMENTS, GM, RADIUS, Case, doubles, ok, rows

TIMES = [60.0 * k for k in range(1, 31)]


def truth():
    y0 = myb.state_from_elements(GM, ELEMENTS)
    return y0, myb.kepler_propagate_many(GM, y0, TIMES)


def raw_model():
    model = ffi.new("myb_force_model**")
    ok(lib.myb_force_model_create(ffi.NULL, model))
    term = ffi.new("myb_point_mass_term*", {"gm": GM})
    ok(lib.myb_force_model_add(model[0], b"earth", lib.myb_term_point_mass,
                               lib.myb_term_point_mass_jacobian, term, ffi.NULL))
    return model[0], term


class Stations(Case):
    def test_measurements_of_a_station(self):
        env = myb.Environment()
        earth = env.add_body(399, GM)
        station = myb.Station(env, earth, 0.4, 0.8, 120.0, RADIUS, 1.0 / 298.257223563)
        raw = ffi.new("myb_station*")
        ok(lib.myb_station_init(raw, env.pointer, earth,
                                ffi.new("myb_ellipsoid*", [RADIUS, 1.0 / 298.257223563]), 0.4,
                                0.8, 120.0))
        self.assertSame(station.position, tuple(raw.position))
        state = ffi.new("double[6]")
        ok(lib.myb_station_state(raw, 100.0, state))
        self.assertSame(station.state(100.0), tuple(state))
        y = myb.state_from_elements(GM, ELEMENTS)
        values = ffi.new("double[2]")
        partials = ffi.new("double[12]")
        light = ffi.new("myb_light_measure*", {"station": raw, "c": 299792458.0, "gm": GM,
                                               "two_way": 1})
        for made, function, context in (
                (myb.measure_range(station), lib.myb_measure_range, raw),
                (myb.measure_range_rate(station), lib.myb_measure_range_rate, raw),
                (myb.measure_radec(station), lib.myb_measure_radec, raw),
                (myb.measure_azel(station), lib.myb_measure_azel, raw),
                (myb.measure_light_range(station, 299792458.0, GM, two_way=True),
                 lib.myb_measure_light_range, light)):
            ok(function(context, 100.0, ffi.new("double[6]", y), values, partials))
            self.assertSame(made(100.0, y), (doubles(values, made.count),
                                             rows(partials, made.count, 6)))
        self.assertEqual(myb.measure_light_azel(station, 299792458.0, aberration=True).count, 2)
        # a measurement at an offset point with an attitude law in Python
        offset = myb.measure_offset(myb.measure_range(station),
                                    lambda t, x, v: (1, 0, 0, 0, 1, 0, 0, 0, 1),
                                    (1.0, 0.0, 0.0))
        moved = offset(100.0, y)
        self.assertNotEqual(moved[0], made(100.0, y)[0])
        at = myb.Station.at(env, earth, station.position, station.enu)
        self.assertSame(myb.measure_range(at)(10.0, y), myb.measure_range(station)(10.0, y))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.measure_range(station), 0.0,
                          station.state(0.0)[:3] + (0.0, 0.0, 0.0))
        self.assertStatus(lib.MYB_E_DOMAIN, myb.Station, env, earth, 0.0, 0.0, 0.0, -1.0)


class Fit(Case):
    def test_a_fit_of_positions(self):
        y0, states = truth()
        guess = [v + d for v, d in zip(y0, (120.0, -80.0, 40.0, 0.1, -0.05, 0.02))]
        model = myb.ForceModel()
        model.add_point_mass("earth", GM)
        observations = [myb.Observation(t, myb.measure_position(), s[:3], (1.0, 1.0, 1.0))
                        for t, s in zip(TIMES, states)]
        result = myb.fit_orbit(model, 0.0, guess, observations)
        # the same through the C functions
        raw, term = raw_model()
        made = ffi.new("myb_observation[]", len(TIMES))
        for k, (t, s) in enumerate(zip(TIMES, states)):
            made[k].t = t
            made[k].measure = lib.myb_measure_position
            made[k].count = 3
            made[k].value[0:3] = list(s[:3])
            made[k].sigma[0:3] = [1.0, 1.0, 1.0]
        state = ffi.new("double[6]", guess)
        out = ffi.new("myb_fit_result*")
        covariance = ffi.new("double[36]")
        residuals = ffi.new("double[]", 3 * len(TIMES))
        ok(lib.myb_fit_orbit(raw, 0.0, state, made, len(TIMES), ffi.NULL, ffi.NULL, out,
                             covariance, residuals))
        self.assertSame(result.state, tuple(state))
        self.assertSame(result.covariance, rows(covariance, 6, 6))
        self.assertSame(result.residuals, doubles(residuals, 3 * len(TIMES)))
        self.assertSame((result.iterations, result.converged, result.values, result.rms_before,
                         result.rms), (out.iterations, out.converged != 0, out.values,
                                       out.rms_before, out.rms))
        self.assertTrue(result.converged)
        self.assertLess(max(abs(a - b) for a, b in zip(result.state[:3], y0[:3])), 1e-3)
        lib.myb_force_model_destroy(raw)
        self.assertStatus(lib.MYB_E_DOMAIN, myb.fit_orbit, model, 0.0, guess, [])

    def test_a_parameter_and_a_measurement_in_python(self):
        y0, states = truth()
        push = ffi.new("double*", 0.0)
        model = myb.ForceModel()
        model.add_point_mass("earth", GM)
        term = model.add_term("push", lambda t, x, v: (push[0], 0.0, 0.0))
        index = model.add_parameter("push", push, 1e-9)
        model.bind_parameter(index, term)
        self.assertEqual(model.parameter_names, ["push"])
        calls = []

        def position(t, y, want):
            calls.append(want)
            return y[:3], ([1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
                           if want else None)

        mine = myb.Measurement.custom(position, 3)
        observations = [myb.Observation(t, mine, s[:3], (1.0, 1.0, 1.0))
                        for t, s in zip(TIMES, states)]
        result = myb.fit_orbit(model, 0.0, y0, observations,
                               prior_sigma=[0.0] * 6 + [1e-6])
        self.assertEqual(len(result.covariance), 7)
        self.assertLess(abs(model.parameter(index)), 1e-9)
        self.assertTrue(any(calls))

        def refuses(t, y, want):
            raise OverflowError("the sky is closed")

        broken = [myb.Observation(60.0, myb.Measurement.custom(refuses, 1), (1.0,), (1.0,))]
        with self.assertRaises(OverflowError):
            myb.fit_orbit(model, 0.0, y0, broken)

    def test_systematic_shifts(self):
        env = myb.Environment()
        earth = env.add_body(399, GM)
        station = myb.Station(env, earth, 0.4, 0.8, 120.0, RADIUS)
        y0, states = truth()
        rng = myb.measure_range(station)
        observations = [myb.Observation(t, rng, (rng(t, s, False)[0][0] + 5.0,), (1.0,), bias=1)
                        for t, s in zip(TIMES, states)]
        model = myb.ForceModel()
        model.add_point_mass("earth", GM)
        result = myb.fit_orbit(model, 0.0, y0, observations, bias=[0.0], fixed=[1] * 6 + [0])
        raw, term = raw_model()
        made = ffi.new("myb_observation[]", len(observations))
        for k, observation in enumerate(observations):
            made[k] = observation.raw[0]
        shifts = ffi.new("double[1]")
        fixed = ffi.new("unsigned char[7]", [1] * 6 + [0])
        options = ffi.new("myb_fit_options*", {"fixed": fixed, "bias_count": 1, "bias": shifts})
        state = ffi.new("double[6]", y0)
        out = ffi.new("myb_fit_result*")
        covariance = ffi.new("double[49]")
        ok(lib.myb_fit_orbit(raw, 0.0, state, made, len(observations), options, ffi.NULL, out,
                             covariance, ffi.NULL))
        self.assertSame((result.bias, result.state, result.covariance),
                        (doubles(shifts, 1), tuple(state), rows(covariance, 7, 7)))
        self.assertLess(abs(result.bias[0] - 5.0), 1e-6)
        lib.myb_force_model_destroy(raw)
        kalman = myb.Filter(model, bias_count=1)
        self.assertEqual(kalman.unknowns, 7)
        kalman.start(0.0, y0, [[float(r == c) for c in range(7)] for r in range(7)], 3000.0)
        kalman.set_bias([5.0])
        self.assertEqual(kalman.bias, (5.0,))
        kalman.update(observations[0])
        self.assertStatus(lib.MYB_E_DOMAIN, kalman.set_bias, [float("nan")])


class Filtering(Case):
    def test_the_filter_and_the_smoother(self):
        y0, states = truth()
        p0 = [[100.0 if r == c and r < 3 else 0.01 if r == c else 0.0 for c in range(6)]
              for r in range(6)]
        model = myb.ForceModel()
        model.add_point_mass("earth", GM)
        kalman = myb.Filter(model, process_noise=1e-12, history=64)
        guess = [v + d for v, d in zip(y0, (5.0, -3.0, 2.0, 0.01, 0.0, -0.01))]
        kalman.start(0.0, guess, p0, TIMES[-1])
        innovations = []
        for t, s in zip(TIMES[:10], states):
            innovations.append(kalman.update(myb.Observation(t, myb.measure_position(), s[:3],
                                                             (1.0, 1.0, 1.0))))
        kalman.predict(TIMES[12])
        now = kalman.get()
        steps = kalman.steps
        kalman.smooth()
        first = kalman.step(0)
        # the same through the C functions
        raw, term = raw_model()
        f = ffi.new("myb_filter**")
        ok(lib.myb_filter_create(raw, ffi.new("myb_filter_options*", {"process_noise": 1e-12,
                                                                     "history": 64}),
                                 ffi.NULL, f))
        ok(lib.myb_filter_start(f[0], 0.0, ffi.new("double[6]", guess),
                                ffi.new("double[36]", [v for row in p0 for v in row]),
                                TIMES[-1]))
        observation = ffi.new("myb_observation*")
        innovation = ffi.new("double[3]")
        for k, (t, s) in enumerate(zip(TIMES[:10], states)):
            observation.t = t
            observation.measure = lib.myb_measure_position
            observation.count = 3
            observation.value[0:3] = list(s[:3])
            observation.sigma[0:3] = [1.0, 1.0, 1.0]
            ok(lib.myb_filter_update(f[0], observation, innovation))
            self.assertSame(innovations[k], doubles(innovation, 3))
        ok(lib.myb_filter_predict(f[0], TIMES[12]))
        t = ffi.new("double*")
        y = ffi.new("double[6]")
        p = ffi.new("double[36]")
        lib.myb_filter_get(f[0], t, y, p)
        self.assertSame(now, (t[0], tuple(y), rows(p, 6, 6)))
        self.assertEqual(steps, lib.myb_filter_steps(f[0]))
        ok(lib.myb_filter_smooth(f[0]))
        lib.myb_filter_step(f[0], 0, t, y, p)
        self.assertSame(first, (t[0], tuple(y), rows(p, 6, 6)))
        lib.myb_filter_destroy(f[0])
        lib.myb_force_model_destroy(raw)
        self.assertEqual(kalman.bias, ())
        self.assertStatus(lib.MYB_E_DOMAIN, kalman.predict, 0.0)  # behind the filter
        self.assertStatus(lib.MYB_E_DOMAIN, myb.Filter, model, process_noise=-1.0)
        unsmoothed = myb.Filter(model)
        unsmoothed.start(0.0, y0, p0, 100.0)
        self.assertStatus(lib.MYB_E_DOMAIN, unsmoothed.smooth)

