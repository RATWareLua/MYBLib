# Classical integrators: DOP853, RKN12(10), Gauss–Jackson

**English** | [Русский](classical-integrators.ru.md) · [All passports](index.en.md)

## What it is

Three integrators of ordinary differential equations that know nothing about physics.
DOP853 is the explicit Runge–Kutta pair 8(5,3) of Dormand and Prince with step control and a
dense output of order 7, for any smooth first-order system y′ = f(t, y) (an orbit with its
variational equations, say). RKN12(10) is an explicit Runge–Kutta–Nyström pair of orders 12
and 10 for x″ = a(t, x), forces that do not depend on the velocity, made for very tight
tolerances. Gauss–Jackson is a fixed-step multistep predictor–corrector of order 8 in summed
ordinate form for x″ = a(t, x, v). All three update the state and the time with compensated
summation, allocate memory once at creation, give the state between steps and locate events
on it. DOP853 and RKN12(10) can also drive `myb_orbit` in place of the default
[IAS15](ias15.en.md); the fixed-step symplectic methods have their own passport
([symplectic integrators](symplectic-integrators.en.md)).

## Source

- DOP853: E. Hairer, S. P. Nørsett, G. Wanner, *Solving Ordinary Differential Equations I*,
  2nd ed., 1993, §II.5–II.6 (the pair, Hairer's combined error estimate, the dense output of
  order 7). The coefficients are copied as text from SciPy
  (`scipy/integrate/_ivp/dop853_coefficients.py`, BSD-3, retrieved 2026-11-07); the structure
  of the error estimate and of the dense output follows `scipy/integrate/_ivp/rk.py` and
  Hairer; no code is copied. One deliberate difference from Hairer's code: a step cut short at
  `t_limit` does not shorten the next one (the planned step is kept).
- RKN12(10): J. R. Dormand, M. E. A. El-Mikkawy, P. J. Prince, *IMA J. Numer. Anal.* 7, 423
  (1987), 17 stages. The coefficients are the exact fractions of `DPRKN12Tableau` of
  OrdinaryDiffEq.jl (MIT, revision `6a6d63e6`, retrieved 2026-11-07), rounded to double; no
  code is copied. The pair has no interpolant of its own; the optional interpolant of a step
  (phase 22) is the library's construction: the acceleration as a polynomial of degree 12 on
  the Chebyshev–Lobatto nodes of the step, integrated twice, its inner values refined by
  Picard's iteration from the quintic Hermite polynomial of the ends.
- Gauss–Jackson: M. M. Berry, L. M. Healy, *J. Astronaut. Sci.* 52, 331 (2004). The
  coefficients were derived anew from the operator identities hv − s = P₁(D)h²a,
  x − S = P₂(D)h²a (Bernoulli numbers) on the interpolant of degree 8 over the window −4…4, as
  exact rational numbers; the predictor row equals the published one.
- Licences of the borrowed coefficients: SciPy BSD-3, OrdinaryDiffEq.jl MIT (texts in
  [`THIRD_PARTY_NOTICES.md`](../../THIRD_PARTY_NOTICES.md)).

## Interface

Header [`myb/integ.h`](../../include/myb/integ.h):

- DOP853: functions `myb_dop853_create`, `myb_dop853_destroy`, `myb_dop853_start`,
  `myb_dop853_step`, `myb_dop853_integrate`, `myb_dop853_time`, `myb_dop853_state`,
  `myb_dop853_statistics`, `myb_dop853_dense`, `myb_dop853_locate`; structures
  `myb_dop853_options`, `myb_dop853_stats`; the right-hand side `myb_ode_function`;
- RKN12(10): the same set as `myb_rkn12_*`; structures `myb_rkn12_options` (with
  `interpolant`), `myb_rkn12_stats`; the acceleration `myb_position_accel_function`;
- Gauss–Jackson: the same set as `myb_gauss_jackson_*`; structures
  `myb_gauss_jackson_options`, `myb_gauss_jackson_stats`; the acceleration
  `myb_accel_function`;
- events: `myb_event_function`.

In [`myb/orbit.h`](../../include/myb/orbit.h): `MYB_ORBIT_DOP853` and `MYB_ORBIT_RKN12` of
`myb_orbit_integrator` (the field `integrator` of `myb_orbit_options`); for RKN12 every term
of the model is declared a function of the place alone with `myb_force_model_set_position_only`
of [`myb/model.h`](../../include/myb/model.h).

## Inputs and outputs

- Units are the caller's; t is a double counted from an epoch the caller chooses (keep |t|
  moderate for full resolution).
- DOP853: a state of n components and f(t, y) → dy/dt. RKN12(10) and Gauss–Jackson: the state
  y = (x[0..n−1], v[0..n−1]), 2n values; the acceleration of RKN12 has no velocity argument,
  that of Gauss–Jackson may depend on v.
- Error norm of DOP853 and RKN12: err = √(mean_i (e_i / s_i)²),
  s_i = atol + rtol·max(|y_i|, |y_new_i|). DOP853: over the n components, controller
  0.9·err^(−1/8) bounded to [0.2, 10]. RKN12: over the 2n components of (x, v), controller
  0.9·err^(−1/11) bounded to [0.2, 5].
- Gauss–Jackson: a fixed step h, its sign the direction; nodes on the grid t0 + k h; the
  corrector iterates to a relative change of a of `tolerance` (10⁻¹⁵ by default), at most
  `max_corrections` (10) times. There is no error control: `max_difference` of the statistics
  (the largest |x_corrected − x_predicted| / |x| of a step) estimates the local error, to choose
  h; `unconverged` counts steps whose corrector stopped at the limit.
- Dense output: DOP853 — order 7 inside the last accepted step; RKN12 — one more step of the
  pair from the start of the last step, as accurate as the steps, or with `interpolant` the
  interpolant of the step, built on the first request inside it; Gauss–Jackson — the exact
  integral of the window's interpolant of the accelerations (degree 8) anywhere within 4 steps
  of the window's centre, which covers the last 4 steps.
- Events: a sign change of g(t, y) between the ends of the last step, located on the dense
  output by the Illinois variant of regula falsi to the resolution of t; `*found = 0` without
  one.
- Statistics: accepted and rejected steps, evaluations (`int64_t`, so that 10⁹ revolutions do
  not overflow a 32-bit `long`).
- In `myb_orbit`: rtol 10⁻¹² and atol 10⁻⁹ (m, m/s) by default; with DOP853 the variational
  equations are part of one first-order system under the same tolerances; RKN12 carries none.
- A fixed order of operations: one input gives the same bits; all three are in the golden
  reproducibility workload.

## Domain of validity

- DOP853: smooth first-order systems of any dimension, forces of any dependence on the
  velocity.
- RKN12(10): forces that depend on time and place only. `myb_orbit` refuses with
  `MYB_E_DOMAIN` a model with a term not declared position-only (at the start) and the
  variational equations (at creation).
- Gauss–Jackson: the start iterates the method's own formulas over the nodes t0 + j h,
  j = −4…4, so the acceleration must be defined 4 steps before t0 as well. The grid runs one
  way: a `t_end` behind the current node is refused (`MYB_E_DOMAIN`). A step too long either
  stops the start (`MYB_E_NOCONV`) or gives a predictor–corrector difference of the size of the
  orbit (checked with h = P/3); a stalled iteration counts as converged only below 10⁻¹⁰.
- Errors at creation: `MYB_E_DOMAIN` (invalid options, n = 0; for Gauss–Jackson h = 0 or not
  finite, negative options), `MYB_E_NOMEM`, `MYB_E_OVERFLOW`. While stepping:
  `MYB_E_STEP_MIN` (the step falls below `h_min`, e.g. at the singularity of y = 1/(1 − t)),
  `MYB_E_NOCONV` (DOP853: the step budget; RKN12: 100 rejections in a row), `MYB_E_RANGE`
  (dense output outside the last step or the window), and any status of the callback, returned
  unchanged.

## Accuracy

DOP853 — [`test_dop853`](../../tests/unit/test_dop853.c) (two-body orbits against
`myb_kepler_propagate`):

| What | Reference | Result |
|---|---|---|
| Order conditions | all trees to order 5, quadratures to order 8 | satisfied (test bound 10⁻¹⁴) |
| Observed order, fixed steps, e = 0.5 | halving the step | 7.4 → 7.8, converging to 8 |
| 10 revolutions, e = 0.7, rtol 10⁻¹³ | Kepler motion | 1.3·10⁻¹⁰ |
| 10⁴ revolutions, e = 0.3, rtol 10⁻¹³ | Kepler motion | 2·10⁻⁵ (linear energy drift) |
| Dense output inside the steps (e = 0.6, 20 revolutions) | Kepler motion from the start of the step | 2.4·10⁻¹² |
| 40 apsides as events | r·v = 0 on the analytic arc | 7.7·10⁻¹⁰ s |
| Transition matrix from the variational equations, rtol 10⁻¹³ | closed form `myb_kepler_stm` | 2·10⁻¹³ (test bound 10⁻¹⁰) |
| 3 revolutions forward and back, rtol 10⁻¹³ | the start | below the test bound 10⁻¹⁰ |

RKN12(10) — [`test_rkn12`](../../tests/unit/test_rkn12.c),
[`test_mpfr_rkn12`](../../tests/mpfr/test_mpfr_rkn12.c):

| What | Reference | Result |
|---|---|---|
| Every coefficient; order conditions | the exact fractions in GMP rationals | each double the nearest to its fraction; quadratures of order 12, b = b′(1 − c), moments of the order-10 solution, row sums c²/2 exact; the 50-digit rows 14–17 to 10⁻⁵³ |
| Observed order, fixed steps, e = 0.5 | halving the step | 12.4 |
| 10 revolutions, e = 0.9, rtol 10⁻¹⁴ | Kepler motion | 4.8·10⁻¹⁵ of the period, 758 steps |
| Apsides as events | r·v = 0 on the analytic orbit | 1.5·10⁻¹⁰ s |
| Interpolant inside the steps, e = 0.6, rtol 10⁻¹⁴ | Kepler motion | 1.9·10⁻¹¹ s inside against 2.2·10⁻¹¹ s at the ends |
| 3 revolutions forward and back | the start | below the test bound 10⁻¹² |
| Pythagorean three-body problem, t = 10 ([`test_nbody`](../../tests/unit/test_nbody.c)) | IAS15 | 4·10⁻¹² |

Gauss–Jackson — [`test_gauss_jackson`](../../tests/unit/test_gauss_jackson.c),
[`test_mpfr_gauss_jackson`](../../tests/mpfr/test_mpfr_gauss_jackson.c):

| What | Reference | Result |
|---|---|---|
| All 261 coefficients | recomputed in MPFR (256 bits) | bit for bit; predictor row = Berry & Healy 2004 |
| Observed order of the positions, e = 0.1 | halving the step | ~10 (the window is exact for an acceleration of degree 8) |
| 100 revolutions, e = 0.1, 200 steps a revolution | Kepler motion | 8.9·10⁻¹³ of the period |
| Charged particle in uniform fields (force depending on v), h = 0.01, 20 time units | closed form | 6·10⁻¹⁷ |
| Dense output and 20 apsides (e = 0.3, 300 steps a revolution) | Kepler motion | within the test bounds: 10⁻⁸ s between nodes, 5·10⁻¹⁴ at nodes, 10⁻⁷ s for events |
| 3 revolutions backward | Kepler motion | within the test bound 10⁻⁸ s |

In `myb_orbit` — [`test_orbit`](../../tests/unit/test_orbit.c),
[`test_events`](../../tests/unit/test_events.c):

| What | Reference | Result |
|---|---|---|
| DOP853, two bodies, 5 revolutions | closed form | 7·10⁻⁴ m, 289 steps |
| DOP853 with the variational equations, 5 revolutions | `myb_kepler_stm` | 3·10⁻¹¹ |
| RKN12 with its interpolant, two bodies, 5 revolutions | closed form | 3.8·10⁻⁷ m; five pericentres to 5.1·10⁻¹¹ s |
| Two days through the edges of the shadow, with breaks | IAS15 at ε = 10⁻¹³ | RKN12 1.08 mm, DOP853 0.81 mm |

## Cost

Evaluation counts (headers and the plan); timings not measured.

- DOP853: the dense output costs 3 more evaluations at its first request after a step.
- RKN12(10): 17 evaluations per accepted step; the dense output of the pair 16 per request,
  also for every iteration of an event search; the interpolant 44–49 per step (three requests
  of the other kind), after that every request in the step is free.
- Gauss–Jackson: one evaluation for the predictor and one per corrector iteration — 2.55 per
  step over 100 revolutions at 200 steps a revolution.

## Limitations

- DOP853 is not symplectic: the energy drifts linearly and the phase error grows
  quadratically — 2·10⁻⁵ after 10⁴ revolutions, where IAS15 keeps 2.8·10⁻¹¹ of the period. In
  the measurement of Ṗ of PSR B1913+16 it gives a spurious Ṗ = −3.5·10⁻¹⁵ (0.15 % of the
  effect); the measurement uses IAS15. Long arcs belong to IAS15.
- RKN12(10) takes no force that depends on the velocity (drag, relativistic terms); in
  `myb_orbit` it carries no variational equations. The interpolant falls back to the step of
  the pair when a step is too long for degree 12 or Picard's iteration does not settle in 12
  rounds. It was checked against exact Kepler motion, not by order conditions (it has none to
  solve).
- Gauss–Jackson has no step control and no variable step, runs one way, needs the force 4
  steps before t0 and is not offered by `myb_orbit`.
- `*_locate` sees only a sign change between the ends of the last step: two roots inside one
  step cancel. A force that is continuous but not smooth (the edge of a shadow) costs the order
  of a step that straddles it; the bare integrators know nothing about breaks — `myb_orbit`
  ends its DOP853 and RKN12 steps at the breaks of the model.

## Data

None at run time: the coefficients are constants of the code (`src/integ/dop853_coefficients.h`
from SciPy, BSD-3; `src/integ/rkn12_coefficients.h` from OrdinaryDiffEq.jl, MIT;
`src/integ/gj_coefficients.h` derived in the library). `test_mpfr_rkn12` reads the exact
fractions from `tests/data/rkn12/dprkn12.txt` (OrdinaryDiffEq.jl, MIT).
