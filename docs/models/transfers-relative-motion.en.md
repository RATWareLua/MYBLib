# Transfers and relative motion

**English** | [Русский](transfers-relative-motion.ru.md) · [All passports](index.en.md)

## What it is

Tools built on the two-body routines of [two-body motion](kepler-regularization.en.md).
Lambert's problem gives every conic arc that joins two positions in a given time, up to a
chosen number of revolutions. The transfer tools give the B-plane of an approach hyperbola,
the patched-conic quantities (sphere of influence, Hill radius, turn of a flyby, departure
impulse) and a general Newton shooting solver for targeting. Relative motion describes a
deputy about a chief in the frame that turns with the chief: the exact transformation, the
linear map of the relative state for any conic of the chief, the Clohessy–Wiltshire closed
form, a two-impulse rendezvous and the relative orbital elements of D'Amico.

## Source

- Lambert: D. Izzo, "Revisiting Lambert's problem", *Celest. Mech. Dyn. Astron.* 121, 1–15
  (2015) — the Lancaster–Blanchard time equation T(x) with Battin's hypergeometric series
  near x = 1. Izzo's Householder iteration is replaced by Newton's method in a guaranteed
  bracket (T decreases monotonically for zero revolutions; for M ≥ 1 its single minimum
  separates the left and right branches). Gooding's method (1990), named in the plan as the
  second method, is not implemented: each solution is checked by propagation instead.
- B-plane: W. Kizner, JPL EPD-457 (1959); R. Sergeyevsky, G. Snyder, R. Cunniff, JPL Publ.
  82-43 (1983); D. Vallado, *Fundamentals of Astrodynamics* (2013), §12.6.
- Patched conics: the sphere of influence of Laplace and the Hill radius, in their textbook
  forms.
- Relative motion: W. H. Clohessy, R. S. Wiltshire, *J. Aerospace Sci.* 27, 653 (1960);
  G. W. Hill (1878); K. Yamanaka, F. Ankersen, *J. Guid. Control Dyn.* 25, 60 (2002) — the map
  is built here for any conic from the two-body transition matrix (`myb_kepler_stm`) taken
  into the turning frame at both ends, the same linear map without a formula per kind of
  orbit. Relative elements: S. D'Amico, PhD thesis, TU Delft (2010).
- Implemented from the papers; no third-party code, no data files.

## Interface

Header [`myb/lambert.h`](../../include/myb/lambert.h): function `myb_lambert`, structure
`myb_lambert_solution`.

Header [`myb/transfer.h`](../../include/myb/transfer.h):

- B-plane: `myb_bplane_compute`, structure `myb_bplane`;
- patched conics: `myb_sphere_of_influence`, `myb_hill_radius`, `myb_flyby_turn`,
  `myb_departure_impulse`;
- shooting: `myb_shoot`, `myb_shoot_work_size`, structure `myb_shoot_options`, callback
  `myb_shoot_function`.

Header [`myb/relative.h`](../../include/myb/relative.h): `myb_relative_from_inertial`,
`myb_relative_to_inertial`, `myb_relative_stm`, `myb_cw_stm`, `myb_relative_rendezvous`,
`myb_roe_from_states`, `myb_roe_to_position`; the state `myb_pv` of
[`myb/kepler.h`](../../include/myb/kepler.h).

## Inputs and outputs

- μ (gm) in m³/s², positions in m, velocities in m/s, times in s, angles in radians.
- Lambert: tof > 0; `retrograde == 0` selects motion counter-clockwise about +z (angular
  momentum with a positive z component). At most 1 + 2 `max_revolutions` solutions, each with
  v₁, v₂, the number of complete revolutions and the branch (0 for zero revolutions, −1 left,
  +1 right for M ≥ 1); multi-revolution solutions that do not exist for this tof are absent.
- B-plane: S along the incoming asymptote, T = S × k / |S × k| for the caller's pole k (the
  ecliptic or the equator of the target), R = S × T; B from the target to the asymptote,
  perpendicular to S. Fields: v∞, C3 = v∞², S, T, R, B, B·T, B·R, |B|, θ = atan2(B·R, B·T),
  the periapsis radius, the eccentricity and the two-body time to periapsis (negative after
  it).
- Patched conics: sphere of influence a (m/M)^(2/5); Hill radius a (1 − e) (m/(3M))^(1/3);
  turn of a flyby 2 asin(1 / (1 + r_p v∞²/gm)); impulse from a circular orbit of radius r to
  a hyperbola of excess v∞, √(v∞² + 2 gm/r) − √(gm/r).
- Shooting: f(u) = 0 with n unknowns and m equations; the minimum-norm step for m < n, least
  squares for m > n, the step halved at most 10 times while |f| does not decrease; the
  tolerance is on max |f_i|, 20 iterations by default. The residual function fills f and, if
  asked, the m × n Jacobian (row-major) — e.g. from a state transition matrix. Work memory of
  `myb_shoot_work_size(n, m)` doubles comes from the caller (no allocation).
- Relative state: R along the position of the chief, N along its angular momentum, T = N × R
  (the Hill frame, RTN of [frames and geodesy](frames-geodesy.en.md)); six numbers — the
  position of the deputy from the chief (m) and its velocity as seen in the turning frame
  (m/s). `phi` is 6 × 6, row-major: relative(t₀ + Δt) = φ relative(t₀). `myb_cw_stm` takes
  the mean motion n (rad/s, n ≠ 0).
- Rendezvous: Δv₀ at t₀ takes the deputy to the relative position `target` at t₀ + Δt, Δv₁
  there leaves it with `target_velocity` (NULL: at rest with the chief); both in the turning
  frame at their times.
- Relative elements: (δa, δλ, δe_x, δe_y, δi_x, δi_y) with δa = (a_d − a)/a,
  δλ = (u_d − u) + (Ω_d − Ω) cos i, δe_x, δe_y the differences of e cos ω, e sin ω,
  δi_x = i_d − i, δi_y = (Ω_d − Ω) sin i, u the mean argument of latitude.
  `myb_roe_to_position` gives R, T, N to first order in the elements and the eccentricity.

## Domain of validity

- Lambert: `MYB_E_DOMAIN` for non-finite input, r₁ or r₂ zero, r₁ and r₂ collinear (the plane
  of the transfer is undefined); the test also has tof ≤ 0 and μ ≤ 0 refused. `MYB_E_NOCONV`
  (iteration budget) is not expected. The random test skips geometries with
  |r₁ × r₂| < 0.02 |r₁| |r₂|.
- B-plane: `MYB_E_DOMAIN` for an orbit that is not hyperbolic or S parallel to the pole.
- Shooting: `MYB_E_DOMAIN` for n or m outside 1 … 64 or invalid options (a tolerance that is
  not positive); `MYB_E_NOCONV`; the statuses of the residual function. Newton's method is
  local: convergence depends on the first guess.
- Relative motion: `MYB_E_DOMAIN` for a chief at the centre or on a line through it, or
  non-finite input; `myb_relative_stm` adds the errors of `myb_kepler_stm`. The linear map is
  first order: its error is of the order of (separation / distance to the centre) × the
  separation, growing with Δt.
- Rendezvous: `MYB_E_DOMAIN` for a Δt at which the position does not depend on the initial
  velocity in every direction (as after a whole revolution of a circle).
- Relative elements: both orbits must be ellipses (`MYB_E_DOMAIN` otherwise); they are
  singular for a chief in the equator of the frame; the position formula assumes a nearly
  circular chief.

## Accuracy

Lambert in `test_lambert` ([tests/unit/test_lambert.c](../../tests/unit/test_lambert.c)),
transfers in `test_transfer` ([tests/unit/test_transfer.c](../../tests/unit/test_transfer.c)),
relative motion in `test_relative` ([tests/unit/test_relative.c](../../tests/unit/test_relative.c)):

| What | Reference | Result |
|---|---|---|
| Curtis, Example 5.2 | the printed v₁, v₂ | within 0.06 m/s per component |
| Lambert, 20 000 random geometries (radii 6 600–39 600 km, 300 s to 3.8 days, up to 3 revolutions, both senses) | (r₁, v₁) propagated over tof with `myb_kepler_propagate` must reach r₂ with v₂ | < 10⁻¹⁰ relative in r and v (test bound; the README gives a worst of 7·10⁻¹²); sense of motion and number of revolutions checked |
| Lambert by shooting with the closed-form transition matrix | `myb_lambert` | 4·10⁻¹² m/s (test bound 10⁻⁸) |
| B-plane of a hyperbola (r_p 7000 km, v∞ 5 km/s), from 2 days before periapsis | closed forms of \|B\|, e, v∞, C3, time to periapsis | \|B\| < 10⁻⁹, r_p < 10⁻⁸, e, v∞, C3 < 10⁻¹² relative; time < 10⁻⁶ s |
| B-plane geometry | the asymptote 10⁵ days in the past (it converges as 1/r) | S < 10⁻⁶, B·T and B·R < 10⁻⁵ of \|B\| |
| The same hyperbola before and after periapsis | each other | B·T, B·R < 10⁻⁶ m |
| B-plane targeting of an Earth flyby perturbed by the Moon (held fixed), variational transition matrix | the target B·T, B·R | 3 iterations to a miss of 2·10⁻⁸ m (test bound: ≤ 8 iterations, 1 m) |
| Sphere of influence and Hill radius of the Earth | 925 000 km, 1 471 000 km | within 2000 km and 5000 km |
| Minimum-norm step | Jᵀ(JJᵀ)⁻¹ f | < 10⁻¹⁴ |
| Relative state there and back | the input | < 10⁻⁹ relative; the separation keeps its length to < 10⁻¹² |
| Linear map, separations 10 m, 100 m, 1 km: a circle over 0.8 period, an ellipse e = 0.7 through the pericentre, a hyperbola e = 1.5 | the exact motion of both | the error grows 100 times per decade of separation (second order); 1 km leaves < 100 m |
| Linear map on a circle | Clohessy–Wiltshire closed form | 3·10⁻¹⁴ (test bound 10⁻¹⁰) |
| Rendezvous from 5 km behind and 300 m below (chief e = 0.2) to 100 m behind at rest in 40 min, flown exactly | the target position | misses by 11.6 m, the error of the linearisation (test bound 30 m) |
| Relative elements over a revolution (chief a = 7000 km, e = 0.001) | constant but δλ, which drifts by −3/2 δa per radian | 2·10⁻¹⁰ (test bound 10⁻⁹); the position they give within 10 m of the exact one |

Single values are measured (`docs/PLAN.en.md`, phases 8 and 19); values with < are the bounds
the tests assert.

## Cost

Not measured.

## Limitations

- Lambert's problem is the two-body one; r₁ and r₂ collinear (transfer angle 0 or π) are
  refused. A transfer in a perturbed field is found with `myb_shoot` and a transition matrix
  (see [propagation and uncertainty](propagation-uncertainty.en.md)).
- The patched-conic functions are the formulas only; the B-plane and its time to periapsis
  are those of the two-body hyperbola through the given state.
- `myb_shoot` handles at most 64 unknowns and 64 equations; it is a local Newton method with
  step halving, not a global search.
- The relative-motion map is linear in the separation and Keplerian: no J₂, drag or other
  perturbation of the relative motion. The relative elements are singular for an equatorial
  chief, and the position they give is first order in e (a nearly circular chief).

## Data

None.
