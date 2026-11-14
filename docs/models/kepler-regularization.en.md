# Two-body motion and regularisation

**English** | [Русский](kepler-regularization.ru.md) · [All passports](index.en.md)

## What it is

The mathematics of one central body: Kepler's equation for the ellipse, the hyperbola and the
parabola written without cancellation, the conversions between anomalies, classical and
modified equinoctial elements, the propagation of a state over any conic by universal
variables, its state transition matrix in closed form, and mean elements found by averaging
a sampled perturbed orbit. For perturbed motion near a pericentre the library also gives
regularised equations: Kustaanheimo–Stiefel (KS), KS with the time element, and the Sundman
transformation with any exponent, integrated with DOP853 up to a given physical time.

## Source

- Two-body formulas: R. H. Battin, *An Introduction to the Mathematics and Methods of
  Astrodynamics* (1999); D. A. Vallado, *Fundamentals of Astrodynamics and Applications*,
  5th ed. (2022); O. Montenbruck, E. Gill, *Satellite Orbits* (2000).
- Modified equinoctial elements: M. J. H. Walker, B. Ireland, J. Owens, *Celest. Mech.* 36
  (1985); R. A. Broucke, P. J. Cefola, *Celest. Mech.* 5 (1972).
- State transition matrix: S. W. Shepperd, *Celest. Mech.* 35, 129 (1985); Battin (1999),
  §9.7 — the Lagrange coefficients differentiated through the universal Kepler equation.
- KS: P. Kustaanheimo, E. Stiefel, *J. Reine Angew. Math.* 218, 204 (1965); E. Stiefel,
  G. Scheifele, *Linear and Regular Celestial Mechanics* (1971); the time element from their
  sec. 18.
- Sundman: K. F. Sundman, *Acta Math.* 36, 105 (1913); the exponent 3/2 is the intermediate
  anomaly of P. Nacozy, *Celest. Mech.* 16, 309 (1977).
- Implemented from the books and papers; no third-party code, no data files. All roots are
  found by one bracketed Newton solver (the rtsafe rule) of the library.

## Interface

Header [`myb/kepler.h`](../../include/myb/kepler.h):

- Kepler's equation: `myb_kepler_elliptic`, `myb_kepler_hyperbolic`, `myb_kepler_parabolic`;
- anomalies: `myb_anomaly_true_from_eccentric`, `myb_anomaly_eccentric_from_true`,
  `myb_anomaly_mean_from_eccentric`, `myb_anomaly_true_from_hyperbolic`,
  `myb_anomaly_hyperbolic_from_true`, `myb_anomaly_mean_from_hyperbolic`;
- elements: `myb_elements_from_pv`, `myb_pv_from_elements`, `myb_equinoctial_from_pv`,
  `myb_pv_from_equinoctial`, `myb_mean_elements_average`;
- propagation: `myb_kepler_propagate`, `myb_kepler_stm`;
- structures `myb_pv`, `myb_elements`, `myb_equinoctial`, `myb_mean_elements`.

Header [`myb/regular.h`](../../include/myb/regular.h):

- `myb_ks_from_cartesian`, `myb_ks_to_cartesian`, `myb_ks_rhs`;
- `myb_ks_element_from_cartesian`, `myb_ks_element_to_cartesian`, `myb_ks_element_rhs`;
- `myb_sundman_rhs`, `myb_regular_propagate`;
- structure `myb_regular_problem`, enumeration `myb_regular_kind` (`MYB_REGULAR_KS`,
  `MYB_REGULAR_SUNDMAN`, `MYB_REGULAR_KS_ELEMENT`), callback `myb_perturbation_function`;
- the options `myb_dop853_options` of [`myb/integ.h`](../../include/myb/integ.h)
  (see [classical integrators](classical-integrators.en.md)).

## Inputs and outputs

- μ in m³/s², positions in m, velocities in m/s, times in s, angles in radians; vectors in
  any inertial frame centred on the body.
- Classical elements carry the semi-latus rectum p = h²/μ instead of the semi-major axis, so
  that the parabola is included. Ranges: i ∈ [0, π], Ω and ω ∈ [0, 2π), ν ∈ (−π, π].
  Undefined angles: for e = 0 exactly ω = 0 and ν is the argument of latitude; for i = 0 or π
  exactly Ω = 0 and ω is measured from the x axis; for both ν is the true longitude.
- Modified equinoctial elements: f = e cos(ω + Ω), g = e sin(ω + Ω), h = tan(i/2) cos Ω,
  k = tan(i/2) sin Ω, L = Ω + ω + ν ∈ (−π, π].
- Kepler's equations: E − e sin E = M (E in the same revolution as M); e sinh H − H = M;
  Barker D + D³/3 = M with D = tan(ν/2), M = 2 √(μ/p³) (t − T_periapsis).
- `myb_kepler_stm` returns `stm[6 i + j]` = ∂y_i/∂y0_j, y = (r, v), row-major, and the
  same state as `myb_kepler_propagate` bit for bit.
- Mean elements: a (m), the equinoctial f, g, h, k, the mean longitude λ at `epoch` (the mean
  of the times) in [0, 2π) and its rate (rad/s). Input: at least 3 states of 6 values at
  increasing times, neighbours less than a revolution apart, spanning whole revolutions.
- KS state y = (u, u′, h, t), 10 values: u ∈ R⁴, the fictitious time s with dt = r ds, h =
  μ/r − v²/2 the binding energy (positive on a bound orbit). Of the circle of u that gives x
  the one with u₄ = 0 for x₁ ≥ 0, u₃ = 0 otherwise. With the time element y = (u, u′, h, τ),
  τ = t + (u·u′)/h. Sundman state y = (x, v, t), 7 values, dt = r^α ds.
- The perturbation P(t, x, v) is the acceleration without the central term (m/s²).
- `myb_regular_propagate` returns the state at t₁ exactly (the event t(s) = t₁ located on the
  dense output) and, optionally, the number of right-hand sides.

## Domain of validity

- `myb_kepler_elliptic`: 0 ≤ e < 1, any finite M; `myb_kepler_hyperbolic`: e > 1;
  otherwise `MYB_E_DOMAIN`. `MYB_E_NOCONV` is never expected (the root stays in a bracket).
  Tested for e up to 1 − 2⁻⁵² and down to 2⁻⁶⁰ (ellipse), from 1 + 2⁻⁵² (hyperbola).
- Anomalies: `MYB_E_DOMAIN` for an eccentricity outside the range of the conic or a true
  anomaly beyond the asymptotes of a hyperbola.
- Elements: `MYB_E_DOMAIN` for μ ≤ 0, non-finite input, r = 0 or zero angular momentum
  (rectilinear motion has no plane); for p ≤ 0, e < 0 or 1 + e cos ν ≤ 0 on the way back;
  the equinoctial elements also refuse i = π.
- Propagation: every conic, rectilinear motion too while r > 0, either sign of dt.
  `MYB_E_DOMAIN` for μ ≤ 0, r₀ = 0, non-finite input or a state that overflows (a hyperbola
  propagated far); `MYB_E_NOCONV` (iteration budget) is not expected.
- Mean elements: `MYB_E_DOMAIN` for μ ≤ 0, fewer than 3 states, times not increasing, a
  state that is not an ellipse.
- Time element: bound orbits only (`MYB_E_DOMAIN` for h ≤ 0); an orbit that a perturbation
  unbinds on the way ends with `MYB_E_DOMAIN` or `MYB_E_STEP_MIN` — use `MYB_REGULAR_KS`
  there.
- Sundman: any finite α ≥ 0; `MYB_E_DOMAIN` for r = 0 or α negative or not finite.
- `myb_regular_propagate`: `MYB_E_DOMAIN` for μ ≤ 0 or invalid input, and the statuses of
  DOP853 and of the perturbation.

## Accuracy

Two-body routines in `test_kepler` ([tests/unit/test_kepler.c](../../tests/unit/test_kepler.c))
and `test_mpfr_kepler` ([tests/mpfr/test_mpfr_kepler.c](../../tests/mpfr/test_mpfr_kepler.c)):

| What | Reference | Result |
|---|---|---|
| Roots E (e up to 1 − 2⁻⁵²), H, D | the exact root of the same double input, Newton at 300 bits (MPFR) | ≤ 2 ulp (test bound; the plan gives ≤ 1.2 ulp for D) |
| Propagation: ellipses to e = 0.95, hyperbolas, near-parabolic, near-circular, to about a period | an independent universal-variable propagation at 400 bits (MPFR) | ≤ 4.7·10⁻¹⁵ relative (test bound 10⁻¹⁴) |
| Ellipses over 1–30 periods | the same | ≤ 2·10⁻¹⁵ (test bound 10⁻¹⁴) |
| Vallado, Example 2-4 | the printed state | within half of the last printed digit |
| 20 000 ellipses | the route elements → M + n Δt → Kepler's equation → state | < 10⁻¹¹ relative |
| 30 000 propagations of all kinds | energy and angular momentum; back and forth; t₁ + t₂ in one or two calls | < 10⁻¹²; < 10⁻⁹; < 10⁻⁹ |
| 100 000 classical and equinoctial round trips | the input state | < 10⁻¹² relative |
| Transition matrix, all conics, also over 1–30 periods | central differences of the 400-bit propagation | ≤ 3·10⁻¹⁵ in the scales r₀, v₀ (test bound 10⁻¹⁴) |
| Symplecticity of the matrix | ΦᵀJΦ = J | 3·10⁻¹⁴ (test bound 10⁻¹³) |
| The matrix against DOP853 with variational equations, rtol 10⁻¹³ (`test_dop853`) | — | 2·10⁻¹³ (test bound 10⁻¹⁰) |

Single values and ≤ are measured (`docs/PLAN.en.md`, phases 2–4, 9, 19 and the old debts);
values written with < are the bounds the tests assert.

Mean elements in `test_relative` ([tests/unit/test_relative.c](../../tests/unit/test_relative.c)),
an orbit of a = 7200 km, e = 0.01, i = 0.9 rad under J₂, three windows of six revolutions of
180 samples: the osculating semi-major axis swings by 11.5 km, the mean one keeps to 0.75 m
and e to 4·10⁻⁷; the drift of the node agrees with the first-order theory to 7·10⁻⁴, the
rate of the mean longitude to 4·10⁻⁷ (test bounds 1 m, 3·10⁻⁶, 3·10⁻³, 3·10⁻⁶).

Regularisation in `test_regular` ([tests/unit/test_regular.c](../../tests/unit/test_regular.c)):

| What | Reference | Result |
|---|---|---|
| KS transformation and back | the input | x, v < 5·10⁻¹⁶; bilinear relation < 2·10⁻¹⁶ |
| e = 0.999, 20 revolutions, DOP853 10⁻¹² | `myb_kepler_propagate` | KS 1.2·10⁻¹², Cowell 1.3·10⁻⁸, Sundman α = 1 1.8·10⁻⁸ (relative) |
| J₂, e = 0.7, 3 days (7 revolutions) | IAS15 at ε = 10⁻¹² | KS and Sundman agree to 1 mm at 23 000 km (test bound 3 mm) |
| τ′ on a Kepler orbit | μ/(2h) | < 10⁻¹⁵ |
| e = 0.95, 20 revolutions, KS → KS with the time element | `myb_kepler_propagate` | 10⁻⁶: 1.93·10⁻⁵ → 5.27·10⁻⁶; 10⁻⁹: 1.13·10⁻⁸ → 1.48·10⁻⁹; 10⁻¹²: 5.99·10⁻¹² → 1.05·10⁻¹² |
| Time element with J₂ (e = 0.7, 3 days) | plain KS | 7·10⁻⁶ m |
| Rectilinear fall through r = 0 and back | x = 1 | 1 + 2.3·10⁻¹³ |
| Sundman, α = 0 … 2.5, e = 0.9, 5 revolutions | `myb_kepler_propagate` | < 10⁻⁸ for every α |

## Cost

Right-hand sides of DOP853 measured in `test_regular` (from `docs/PLAN.en.md`, phase 9 and the
old debts):

| Case | Evaluations |
|---|---|
| e = 0.999, 20 revolutions, 10⁻¹² | KS 4181, Cowell 75 278, Sundman α = 1 40 289 |
| e = 0.95, 20 revolutions, 10⁻⁶ / 10⁻⁹ / 10⁻¹², KS → with the time element | 785 → 773, 1793 → 1781, 4217 → 4169 |
| e = 0.9, 5 revolutions, Sundman α = 0 / 0.5 / 1 / 1.25 / 1.5 / 2 / 2.5 | 9317 / 7409 / 6053 / 5249 / 5945 / 5621 / 7337 (Cartesian 9374) |

The work of KS does not depend on the eccentricity. The two-body routines are not timed.

## Limitations

- Two bodies only; perturbations enter only the regularised equations, through the caller's
  `myb_perturbation_function`.
- The classical angles are ill-conditioned near e = 0 and i = 0 (only their sums are well
  defined); the equinoctial elements are singular at i = π. Classical and equinoctial
  elements of a near-parabolic orbit far from the pericentre are ill-conditioned: an ulp of
  ν moves r by (1 − e)⁻¹ ulp — a property of the parametrisation; the tests keep such orbits
  within 2.5 rad of the pericentre.
- Mean elements are a numerical average, first order in the oscillation, with nothing assumed
  about the forces; a part of a revolution left over biases them. They are not the mean
  elements of an analytic theory (for SGP4 see [SGP4 and TLE](sgp4-tle.en.md)).
- `myb_regular_propagate` integrates with DOP853 only; the right-hand sides have the
  signature of `myb_ode_function` for other uses. The time element serves bound orbits only.
- Lambert's problem, B-plane and relative motion built on these routines are in
  [transfers and relative motion](transfers-relative-motion.en.md).

## Data

None.
