# IAS15

**English** | [Русский](ias15.ru.md) · [All passports](index.en.md)

## What it is

An implicit integrator of order 15 with adaptive steps for second-order systems
x″ = a(t, x, v): Gauss–Radau quadrature on 8 nodes in Everhart's formulation, the implicit
equations solved by predictor–corrector iteration until the correction reaches round-off, the
step set by the timescale of the acceleration. With the default precision the truncation error
stays below round-off, so the energy error of a conservative system grows like a random walk
(Brouwer's law). It is the default integrator of `myb_orbit`, of the estimators and of the
[system of bodies](system-of-bodies.en.md). Around the published method the library adds a
dense output from the polynomial of the step, events, blocks and groups of coordinates judged
each by itself, withdrawing and retaking a step at a break of the force (the edge of a shadow),
and a restart after an impulse that inherits the step.

## Source

- H. Rein, D. S. Spiegel, *MNRAS* 446, 1424 (2015); E. Everhart (1985) for the Gauss–Radau
  formulation. Implemented from the papers, not from the REBOUND code (GPL).
- The step criterion PRS23: D. Pham, H. Rein, D. S. Spiegel, *Open J. Astrophys.* 7 (2024),
  eqs. 15–17, τ² = 2|a|² / (|a′|² + |a||a″|), dt = (7! ε)^(1/7)·min τ, with a, a′, a″ at the end
  of the step from the step's own polynomial, each measured by the norm of their eq. 7 (the
  maximum over blocks of the Euclidean norm of a block). Until phase 4 the library took the
  minimum of τ over blocks; the figure-eight stalled on it and the paper's norm was restored.
- The prediction of the next step: the polynomial of the step re-expanded, plus the last
  correction scaled by q^(k+1).
- The Gauss–Radau nodes and the C, D, R tables are computed by the library in MPFR at 512 bits;
  no third-party code or data.

## Interface

Header [`myb/integ.h`](../../include/myb/integ.h):

- functions `myb_ias15_create`, `myb_ias15_destroy`, `myb_ias15_start`, `myb_ias15_step`,
  `myb_ias15_integrate`, `myb_ias15_time`, `myb_ias15_state`, `myb_ias15_statistics`,
  `myb_ias15_dense`, `myb_ias15_locate`;
- structures `myb_ias15_options` (`epsilon`, `h_initial`, `h_max`, `h_min`, `max_steps`,
  `block`, `groups`, `group_count`), `myb_ias15_stats` (`accepted`, `rejected`, `evaluations`,
  `iterations`, `unconverged`);
- callbacks `myb_accel_function`, `myb_event_function`.

Withdrawing a step at a break and the restart after an impulse are reached through
`myb_orbit` ([`myb/orbit.h`](../../include/myb/orbit.h): `myb_orbit_kick`, `myb_orbit_breaks`)
and the breaks of the force model (`myb_force_model_add_break` of
[`myb/model.h`](../../include/myb/model.h)); `myb_system_integrator` of
[`myb/system.h`](../../include/myb/system.h) creates an IAS15 with the groups of a system.

## Inputs and outputs

- The state y = (x[0..n−1], v[0..n−1]), 2n values; the acceleration may depend on v. Units are
  the caller's; t is a double counted from an epoch the caller chooses.
- `epsilon`: 10⁻⁹ by default (`myb_orbit` and the estimators use 10⁻¹⁰ by default).
- Blocks: `block` components per block (0: one block of all n; it must divide n) — 3 for the
  Cartesian vector of a body. The norms of a, a′, a″ are maxima over blocks of the Euclidean
  norm of a block, so a body whose acceleration passes through zero (the middle body of a
  symmetric collinear configuration) does not stall the step.
- Groups: `group_count` sizes, each a multiple of `block`, their sum n (the sizes are copied).
  Without groups one group measures the others in its own units (a quaternion beside metres
  disappears); with them each group sets its own timescale and its own convergence, the
  shortest timescale and the worst convergence counting. `myb_system_integrator` makes the
  positions of all bodies, each rotation and each pair of spins groups of their own. Without
  groups the bits are those of before.
- Step: the next step is at most 4 times the last; a step whose successor would be below a
  quarter of it is rejected and redone. The automatic first step is (7! ε)^(1/7) times the
  smallest of √(|x|/|a|) and |v|/|a| over the blocks (zero terms ignored); a first step that is
  too long is rejected and redone.
- The iteration stops at 10⁻¹⁶ or when the correction stops decreasing, after at most 12
  sweeps of 7 evaluations each; steps stopped at the limit are counted in `unconverged`.
- Dense output: (x, v) at any t inside the last accepted step from the step's polynomial,
  without evaluations.
- Events: as `myb_dop853_locate` — a sign change of g between the ends of the last step,
  located on the dense output by the Illinois variant of regula falsi to the resolution of t.
  `myb_orbit_event` examines the sign at the ends of four equal parts of every step.
- Breaks (in `myb_orbit` and the estimators): after a step the break functions are examined on
  its dense output (four parts of the step, the root on the step's polynomial); a step that
  contains a sign change is withdrawn and taken again from its beginning exactly up to the
  root, with the polynomial of the withdrawn step on the part kept as the prediction (the second
  attempt converges in a sweep or two); the step after the edge starts from the same polynomial
  continued past it. The integrator is not restarted: the state, Φ and S continue. The sign of
  every function is remembered, so an edge is not found twice.
- Restart after an impulse (`myb_orbit_kick`): the velocity jumps at the time reached, the
  integrator starts again with the step inherited from the interrupted integration (one that
  is too long is rejected as usual); the variational equations continue from the same values.
- Statistics are `int64_t`; the order of operations is fixed and the state and time are
  updated with compensated summation: one input gives the same bits.

## Domain of validity

- Second-order systems whose acceleration does not vanish over a whole group at once. τ is
  relative to the norm of a and shrinks to zero where the acceleration of the whole state (of a
  whole group) vanishes: a one-dimensional oscillator stalls with `MYB_E_STEP_MIN` — use
  [DOP853](classical-integrators.en.md) there. A system of gravitating bodies does not get
  there.
- Contact: a contact that holds a body against gravity passes the acceleration through zero at
  every touch, and a ball dropped on a plane stops IAS15 with `MYB_E_STEP_MIN` within its first
  contact; contact problems go with DOP853 ([contact and landing](contact-landing.en.md)).
- Close approaches of point masses are integrated through (the Pythagorean problem, closest
  approach 4.1·10⁻⁴); the contact of two spheres in the system is made the end of a step
  ([encounters in the system](system-encounters.en.md)).
- Errors: `MYB_E_DOMAIN` at creation (n = 0, ε not positive or not finite, negative step
  options, a block that does not divide n, groups that are not multiples of the block or do not
  sum to n) and at the start (a state that is not finite); `MYB_E_NOMEM`, `MYB_E_OVERFLOW`;
  `MYB_E_STEP_MIN` (the step falls below `h_min`); `MYB_E_NOCONV` (100 rejections in a row);
  `MYB_E_RANGE` (dense output outside the last step, or before any step); a status of the
  acceleration is returned unchanged.

## Accuracy

The method — [`test_ias15`](../../tests/unit/test_ias15.c),
[`test_mpfr_radau`](../../tests/mpfr/test_mpfr_radau.c) (two-body orbits against
`myb_kepler_propagate`, ε = 10⁻⁹):

| What | Reference | Result |
|---|---|---|
| Gauss–Radau nodes, C, D, R tables | recomputed in MPFR at 512 bits | bit for bit |
| Local order, one step | shortening the step | 15.6 |
| 10 revolutions, e = 0, 0.5, 0.9, 0.99 | Kepler motion | 2·10⁻¹⁵ … 1.6·10⁻¹³ of the period, no unconverged step |
| 10⁴ revolutions, e = 0.3 | Kepler motion | 2.8·10⁻¹¹ of the period (DOP853 on the same arc: 2·10⁻⁵) |
| Brouwer's law: 12 orbits of e = 0.5, 10–1000 revolutions | random walk, slope 1/2 | slope 0.44; r.m.s. relative energy error 4.6·10⁻¹⁵ after 1000 revolutions |
| Charged particle in uniform fields (force depending on v) | closed form | 3·10⁻¹⁶ |
| Hyperbola e = 2, 3 h either side of pericentre | Kepler motion | within the test bound 10⁻⁹ s |
| Dense output and 40 apsides (e = 0.6, 20 revolutions) | Kepler motion from the start of the step | within the test bounds: 10⁻⁹ s inside steps, 2 ulp of t at their ends, 10⁻⁸ s for events |
| Two orbits of periods 1 : 64 in one state, blocks of 3 or one block | Kepler motion | both within the test bound 10⁻¹⁰ s |
| 3 revolutions forward and back | the start | within the test bound 10⁻¹³ |

Systems of bodies — [`test_nbody`](../../tests/unit/test_nbody.c),
[`test_system`](../../tests/unit/test_system.c):

| What | Reference | Result |
|---|---|---|
| Lagrange equilateral solution (GM 1, 10⁻³, 10⁻⁶), 20 periods | rigid rotation | 1·10⁻¹³ of the side |
| Figure-eight, 10 periods | return to the start | 2·10⁻¹⁴; energy 3·10⁻¹⁶ |
| Pythagorean problem (Burrau) | RKN12(10) | 4·10⁻¹² at t = 10; energy kept to 10⁻¹⁰ |
| Free top over 30 days, rotation and motion as one system with groups | exact free rotation | 3·10⁻¹³; \|q\| − 1 = 5·10⁻¹⁵ (without groups the top diverged, \|q\| − 1 = 10³) |

Breaks and impulses — [`test_events`](../../tests/unit/test_events.c),
[`test_lageos`](../../tests/unit/test_lageos.c),
[`test_nongrav`](../../tests/unit/test_nongrav.c):

| What | Reference | Result |
|---|---|---|
| 12 000 km orbit, 0.2 m²/kg, two days, 54 edges of the shadow | spread of the final position over ε = 10⁻⁹ … 10⁻¹² | 6.4 m without breaks, 8·10⁻⁷ m with them |
| Φ across the edges of 1.5 revolutions; ∂x/∂c_r | central differences of whole propagations | 7·10⁻⁸ with a largest element of 25; 2·10⁻⁷ m of 49 m |
| LAGEOS, residual of the first guess over a week | the same at ε = 10⁻¹⁰, 10⁻¹¹, 10⁻¹² | 25.8887 m at all three with breaks (25.8821 m at 10⁻¹⁰ and 25.8882 m at 10⁻¹² without) |
| Hohmann transfer by two impulses (`myb_orbit_kick`) | the computed final orbit | 1.5·10⁻⁸ m, e = 6·10⁻¹⁶ |
| Φ across both impulses | differences of the whole manoeuvre | 2·10⁻⁹ of the largest element |

## Cost

- A sweep of the iteration costs 7 evaluations, at most 12 sweeps a step; the dense output
  costs none.
- Withdrawing the step at breaks costs +26 % of the time of the LAGEOS fit (a cold restart of
  the integrator at every edge cost +70 %, with 352 restarts per propagation); the week of
  LAGEOS with 176 edges takes 2290 steps instead of 2177.

## Limitations

- The published criterion is relative to the norm of a: it stalls where the acceleration of a
  whole group vanishes (a one-dimensional oscillator, contact at every touch); such problems go
  with DOP853.
- On an exactly circular orbit with a constant step the rounding errors of neighbouring steps
  are identical and do not average out: the energy drifts ~10⁻¹⁷ per revolution (3·10⁻¹⁴ over
  3000 revolutions), independently of ε — rounding, not truncation. Asymmetric systems follow
  Brouwer's law.
- Components of different kinds in one state need groups; without them the step is chosen by
  the group with the largest numbers alone.
- The automatic first step is made for a fresh problem with exact zeros: a block whose velocity
  only passes through zero (10⁻¹⁹ instead of 0) gives a timescale of 10⁻¹¹ s. Restarts in
  `myb_orbit` therefore inherit the step; a bare `myb_ias15_start` with the automatic first
  step has the same weakness (`h_initial` sets the first step by hand).
- `myb_ias15_locate` sees only a sign change between the ends of the last step;
  `myb_orbit_event` looks at four parts of each step, so roots closer together than a quarter
  of a step may be missed.
- The bare integrator knows nothing about breaks; withdrawing the step is done by `myb_orbit`,
  the estimators and the system. Ending a step at the edge predicted from the polynomial,
  without withdrawing, was tried and lost: 7·10⁻⁴ m at ε = 10⁻⁹ and 5·10⁻⁶ m at 10⁻¹² against
  8·10⁻⁷ m for withdrawing, since the predicted root is not the true one and the kink stays
  inside the next step.
- Phase 21 adds: the edge of the shadow without withdrawing the step, and contact under IAS15.

## Data

None: the Gauss–Radau nodes and the C, D, R tables are constants of the code
(`src/integ/ias15_coefficients.h`), computed by the library in MPFR.
