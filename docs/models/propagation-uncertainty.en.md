# Force model, propagation and uncertainty

**English** | [Русский](propagation-uncertainty.ru.md) · [All passports](index.en.md)

## What it is

The layer that turns the physics of the library into trajectories. An environment holds the
bodies around the object — each a set of optional numbers and providers (GM, the state from an
ephemeris, the orientation, a radius, a luminosity) evaluated once per instant — and the time
scales of the providers in a table. Ready-made terms compose the physical models with those
values. A force model sums named terms, differentiates them (analytically or by central
differences), registers parameters and declares the breaks where the force is not smooth. An
orbit is propagated by IAS15, DOP853 or RKN12, with states that do not depend on the times asked
for, events, impulsive manoeuvres and, on request, the state transition matrix and the
sensitivities. The uncertainty layer gives the variational equations, covariance transforms, a
deterministic random generator and sample statistics; the budget measures how much each term
moves the object, by full propagations or to first order.

## Source

- IAS15: H. Rein, D. S. Spiegel (2015); the Gauss–Radau constants are computed by the library
  itself (`tests/mpfr/test_mpfr_radau.c`). DOP853: the coefficients of E. Hairer, S. P. Nørsett,
  G. Wanner, copied from SciPy (BSD 3-Clause). RKN12(10): J. R. Dormand, M. E. A. El-Mikkawy,
  P. J. Prince, *IMA J. Numer. Anal.* 7, 423 (1987), the coefficients from OrdinaryDiffEq.jl
  (MIT). See `THIRD_PARTY_NOTICES.md` and the [classical integrators](classical-integrators.en.md)
  and [IAS15](ias15.en.md) passports.
- xoshiro256**: D. Blackman, S. Vigna, *ACM TOMS* 47, 36 (2021), seeded by SplitMix64; normal
  deviates by the Box–Muller transform on the library's correctly rounded log, sin, cos.
- The variational equations, the covariance transforms and the budget are implemented by the
  library; the ready-made terms are compositions of the functions of `myb/grav.h`,
  `myb/force.h` and `myb/rel.h` (their physics is in the passports linked below).

## Interface

- [`myb/env.h`](../../include/myb/env.h): `myb_env_create`, `myb_env_destroy`, `myb_env_clone`,
  `myb_env_add_body`, `myb_env_body_count`, `myb_env_gm`, `myb_env_epoch`, `myb_env_find`,
  `myb_env_time`, `myb_env_state`, `myb_env_orientation`; structures `myb_env_options`,
  `myb_env_body`, `myb_object` (with `myb_mass_function`); the ready-made terms below.
- [`myb/model.h`](../../include/myb/model.h): `myb_force_model_create`, `_destroy`, `_add`,
  `_find`, `_count`, `_name`, `_enabled`, `_set_enabled`, `_set_position_only`,
  `_position_only`, `_add_parameter`, `_parameter_count`, `_parameter_name`, `_parameter`,
  `_bind_parameter`, `_add_break`, `_break_count`, `_break`, `_set_scales`, `_acceleration`,
  `_term`, `_jacobian`, `_accel`; `myb_term_function`, `myb_term_jacobian`,
  `myb_break_function`; `myb_term_no_jacobian`, `myb_point_mass_jacobian`.
- [`myb/orbit.h`](../../include/myb/orbit.h): `myb_orbit_create`, `_destroy`, `_start`, `_time`,
  `_state`, `_states`, `_partials`, `_event`, `_kick`, `_steps`, `_breaks`; `myb_orbit_options`,
  `myb_orbit_integrator` (`MYB_ORBIT_IAS15`, `MYB_ORBIT_DOP853`, `MYB_ORBIT_RKN12`).
- [`myb/uncert.h`](../../include/myb/uncert.h): `myb_variational_create`, `_destroy`, `_start`,
  `_step`, `_integrate`, `_time`, `_get`, `_dense`; `myb_covariance_transform`,
  `myb_covariance_rtn`; `myb_rng` with `myb_rng_seed`, `_next`, `_uniform`, `_normal`,
  `_gaussian`; `myb_sample_statistics`.
- [`myb/budget.h`](../../include/myb/budget.h): `myb_budget`, `myb_budget_linear`;
  `myb_budget_entry`, `myb_budget_options`; flags `MYB_BUDGET_OK`, `MYB_BUDGET_MISSING`,
  `MYB_BUDGET_NEGLIGIBLE`.

The ready-made terms of `myb/env.h`, each a context structure and a `myb_term_function` with its
Jacobian:

| Term | Functions | Jacobian | Physics |
|---|---|---|---|
| Point mass at the origin | `myb_term_point_mass` | analytic | — |
| A body of the environment (direct, minus indirect with a central body) | `myb_term_third_body` | analytic | — |
| Gravity field of a body with its time terms and tides; the indirect term of the figure of the central body | `myb_term_field`, `myb_term_field_indirect` | analytic (the gradient); none for the indirect term | [spherical harmonics](spherical-harmonics.en.md), [IERS tides](tides-iers.en.md) |
| Radiation pressure of a source with the shadows of bodies, optional Poynting–Robertson and refraction | `myb_term_radiation` | analytic; in a penumbra differences of the term | [radiation pressure](radiation-pressure.en.md) |
| Albedo and infrared of a planet (Knocke et al.) | `myb_term_planet_radiation` | analytic | [radiation pressure](radiation-pressure.en.md) |
| Drag in an atmosphere rotating with the body, optional wind | `myb_term_drag` | analytic in v and the rotation; density and wind by differences | [aerodynamics](aerodynamics.en.md), [thermosphere models](thermosphere-models.en.md) |
| Schwarzschild, Lense–Thirring, de Sitter (IERS 2010 eq. 10.12) | `myb_term_schwarzschild`, `myb_term_lense_thirring`, `myb_term_de_sitter` | analytic | [relativity](relativity-1pn.en.md) |
| Non-gravitational acceleration of a small body | `myb_term_nongrav` | analytic; with a delay differences | [comets](comets-outgassing.en.md) |
| Yarkovsky, diurnal and seasonal | `myb_term_yarkovsky` | analytic | [Yarkovsky and YORP](yarkovsky-yorp.en.md) |
| The lagging tide the object raises | `myb_term_tidal_lag` | analytic | [lagged tides](tides-lagged.en.md) |
| Empirical R, T, N: constant and once per revolution | `myb_term_empirical` | analytic | — |

Terms in other headers: `myb_term_surface_radiation`, `myb_term_surface_drag`,
`myb_term_recoil` (`myb/surface.h`, [radiation on surfaces](radiation-surfaces.en.md)),
`myb_thrust_acceleration` and `myb_term_qlaw` ([low thrust](low-thrust-qlaw.en.md)),
`myb_term_contact` ([contact](contact-landing.en.md)), `myb_cr3bp_term`
([three-body problem](cr3bp-periodic.en.md)); the breaks at the edges of shadows
`myb_term_radiation_breaks` ([events](events-geometry.en.md)).

## Inputs and outputs

- Time: the model time t is seconds from `epoch` on a named scale; each provider is called with
  the same instant on its own scale (TDB for the JPL ephemerides, TT for the rotation of the
  Earth — see [time scales](time-scales.en.md)). The difference of two scales comes from a table
  every `table_step` seconds (3600 by default) at epoch + k step, interpolated on four nodes: a
  pure function of t, the same bits whatever was asked before.
- Frames: all bodies and the object share one set of inertial axes. The origin is the body
  without a state provider (its acceleration by the others enters as indirect terms) or, when
  every body has one, the inertial origin of the providers (a barycentre).
- A body: an identifier, GM (m³/s²), a state provider ([SPK ephemerides](spk-ephemerides.en.md)),
  an orientation provider (inertial → body-fixed, [Earth orientation](earth-orientation.en.md)),
  a radius (m), a luminosity (W), limb darkening u1, u2, a spin (m²/s in the body frame); zero
  means absent. The object: mass (kg) or a mass function of t, cross-section (m²), c_r, c_d.
- A term writes its acceleration (it does not add); x and v are relative to the origin. The
  model sums the enabled terms in the order they were added with compensated summation per
  component. Jacobians: dadx[3i + j] = ∂a_i/∂x_j, likewise ∂a/∂v and ∂a/∂p; central differences
  with steps ∛ε·scale (scale = |x| or |v|, at least the length and velocity scales, 1 by
  default; about 10⁻¹⁰ relative). A parameter is a double inside a context, registered by
  address with its own step, optionally bound to the one term it enters.
- A break g(t, y) changes sign where the force is continuous but not smooth; `myb_orbit` and the
  estimators end a step at each of its roots.
- `myb_orbit`: IAS15 (ε 10⁻¹⁰ by default), DOP853 or RKN12 (rtol 10⁻¹², atol 10⁻⁹ m and m/s);
  every state from the dense output of the step that contains it; the steps do not depend on the
  times asked for (only the last is cut at t_end). Times asked in the order of the propagation,
  forward or backward. `myb_orbit_partials`: Φ[6i + j] = ∂y_i/∂y0_j and ∂y_i/∂p_k. Events: the
  sign of g at the ends of four equal parts of every step, a root located on the dense output to
  the resolution of t, a direction filter. `myb_orbit_kick`: the velocity jumps by Δv (m/s,
  inertial axes) at the time reached; Φ and S continue as those of the trajectory with the
  impulse taken as given.
- Variational equations: the columns of Φ = ∂(x, v)/∂(x0, v0) and S = ∂(x, v)/∂p integrated with
  the state as one IAS15 system of 3 (7 + np) coordinates, with dense output (with DOP853 in
  `myb_orbit`, as one first-order system under the same tolerances).
- Covariances: M P Mᵀ; a 6 × 6 covariance in the R, T, N axes of (r, v), without the rotation
  rate of the axes.
- Random numbers: 64-bit words, uniform in [0, 1) with 53 random bits, standard normal, samples
  of N(mean, L Lᵀ); sample mean and covariance (divisor count − 1) with compensated sums.
- The budget: for each term the largest |r_switched − r_reference| over `samples` (100) equal
  times with the term switched, and the largest acceleration of the term; MISSING (off, effect
  above the tolerance), NEGLIGIBLE (on, effect below). `myb_budget_linear`: one propagation,
  every term driving its own variation d_k(t) = ∫ Φ(t, s) a_k(s) ds.

## Domain of validity

- An environment: up to 32 bodies, 4 scales besides the model's, 8 tide raisers, 8 occulters;
  |t| up to 10¹⁷ s or 2⁵² table steps (`MYB_E_RANGE`). An environment and its terms serve one
  thread at a time (`myb_env_clone` for another).
- A model: up to 64 terms, 16 parameters, 32 breaks; names of at most 31 characters.
- RKN12 only for models whose every enabled term is declared a function of the place alone
  (`myb_force_model_set_position_only`; the declaration is not checked), and without variational
  equations. The variational equations go with IAS15 and DOP853.
- Events: roots closer together than a quarter of a step may be missed.
- The linear budget errs to the second order of the effect: right for the small terms, the order
  of magnitude only for those that change the orbit itself.
- Refused with `MYB_E_DOMAIN`: invalid values in a body or a term, the object at the centre of a
  body, times or states not finite, t_end equal to t0, a time behind the time reached or beyond
  t_end, negative options, an unknown integrator, RKN12 with variational equations or a term not
  declared, partials of an orbit created without them, a budget with t1 = t0 or a tolerance ≤ 0,
  statistics of fewer than 2 samples. `MYB_E_NOCONV`: the step budget (10⁹ by default)
  exhausted. `MYB_E_OVERFLOW`: too many bodies, scales, terms, parameters or breaks.
  `MYB_E_NOT_FOUND`: a scale or a body identifier.

## Accuracy

`test_env` ([tests/unit/test_env.c](../../tests/unit/test_env.c)), with providers made for the
test:

| What | Reference | Result |
|---|---|---|
| Every ready-made term | the same functions of grav/force/rel called by hand with the providers' values | bit for bit |
| Jacobian of the third body | central differences | 2·10⁻¹¹ |
| Field of a body away from the origin, degree 0 | Battin's third-body form | 5.5·10⁻¹⁶ |
| Drag | −½ ρ c_d A/m ‖w‖ w, w = v − ω × r | 1.4·10⁻¹⁶ |
| Visible fraction of the source | lit, umbra, penumbra, past it | 1, 0, 0.517, 1 |
| Table of scales, step 1 h, ±40 days | the exact conversion | TDB−TT 3.5·10⁻¹⁶ s; TCB−TT (20 s with a rate) 3.1·10⁻¹⁵ s |

`test_orbit` ([tests/unit/test_orbit.c](../../tests/unit/test_orbit.c)), `test_uncert`
([tests/unit/test_uncert.c](../../tests/unit/test_uncert.c)), `test_events`
([tests/unit/test_events.c](../../tests/unit/test_events.c)), `test_nongrav`
([tests/unit/test_nongrav.c](../../tests/unit/test_nongrav.c)):

| What | Reference | Result |
|---|---|---|
| IAS15, two bodies, 5 revolutions, e = 0.3 | the closed form | 4·10⁻⁸ m |
| Apsides of 4 revolutions, forward and backward | Kepler's equation | 3.6·10⁻¹¹ s |
| The state at a time with and without 500 other samples | itself | bit for bit |
| Transition matrix from `myb_orbit_partials` | the closed form (`myb_kepler_stm`) | 2.6·10⁻¹⁵ |
| Variational equations over 4 revolutions | the closed form | 5.5·10⁻¹⁵ with the analytic Jacobian, 2.3·10⁻⁹ with differences |
| Sensitivity to a pressure parameter; a column of Φ with drag | differences of whole propagations | 3·10⁻⁸ and 9·10⁻⁹ |
| Transition matrix in a rotating field of degree 8, 2 revolutions | symplecticity; differences of propagations | symplectic to 6·10⁻¹⁷ (10⁻¹² with the Jacobian by differences); 5·10⁻¹⁰ |
| Decay of a circular orbit in an atmosphere, 10 revolutions | first order, −25.403 m | −25.408 m |
| DOP853, two bodies, 5 revolutions | the closed form | 7·10⁻⁴ m in 289 steps |
| DOP853 with the variational equations, 5 revolutions | the closed form | 3·10⁻¹¹ |
| RKN12, two bodies, 5 revolutions; five pericentres | the closed form | 3.8·10⁻⁷ m; 5.1·10⁻¹¹ s |
| Breaks: 12 000 km, 0.2 m²/kg, 2 days, 54 edges of the shadow | spread of the final position over ε = 10⁻⁹ … 10⁻¹² | 8·10⁻⁷ m with breaks, 6.4 m without |
| Φ across 8 edges; ∂x/∂c_r | differences of whole propagations | 7·10⁻⁸ (largest element 25); 2·10⁻⁷ m of 49 m |
| DOP853 and RKN12 through the same shadows, two days | IAS15 | 0.81 mm and 1.08 mm |
| A Hohmann transfer by two kicks | the computed orbit | 1.5·10⁻⁸ m, e = 6·10⁻¹⁶; Φ across both kicks 2·10⁻⁹ of the largest element |
| A covariance propagated linearly | 4000 Monte Carlo samples | within 3 % (sampling error 1.1 %) |
| xoshiro256** from the state {1, 2, 3, 4} | the authors' reference values | exact |

`test_lageos` ([tests/unit/test_lageos.c](../../tests/unit/test_lageos.c)): the LAGEOS model of
ready-made terms against the same model assembled by hand — 3.5·10⁻¹⁶ of the acceleration at 42
states; its first-guess residual does not depend on the precision of IAS15 once the edges of the
shadows are breaks (25.8887 m at ε = 10⁻¹⁰, 10⁻¹¹, 10⁻¹²).

`test_transfer` ([tests/unit/test_transfer.c](../../tests/unit/test_transfer.c)), the budget of a
400 km orbit over a day: J2 1100 km, the drag left out flagged (28 km), a term of 10⁻¹⁴ m/s²
negligible (4·10⁻⁶ m), the effect of a small force linear (×2.0000). The linear budget against the
full one: the pressure 36.02082 m against 36.02081 m, the drag 27.76 km against 27.82 km
(0.2 %), J2 1108 km against 1113 km (0.4 %), the same flags; the 10⁻¹⁴ m/s² term 4.2253·10⁻⁶ m,
repeating to 4·10⁻⁷ when another term changes, where two full budgets gave 4.28 and 4.76·10⁻⁶ m.

## Cost

- The environment computes the values of an instant once; the series of TDB−TT (787 sines) is
  evaluated once per table step instead of at every evaluation.
- Breaks: the LAGEOS week takes 2290 steps with its 176 edges against 2177 without; withdrawing
  a step at an edge costs +26 % of the time of the LAGEOS fit (a cold restart of the integrator
  cost +70 %).
- RKN12 with its interpolant: 44–49 evaluations a step, after which every request is free (the
  older dense output: 16 evaluations a request).
- `myb_budget_linear`: one propagation with the Jacobian of the model; a term differentiated by
  differences costs 13 values per value of the model — in the test model J2 was computed
  334 863 times against 80 335 for the six full propagations.

## Limitations

- The library knows no body: every number and every provider is the caller's.
- Terms without an analytic Jacobian are differentiated by central differences of the term
  (about 10⁻¹⁰ relative); the radiation term in a penumbra and the density and the wind of the
  drag term are such.
- RKN12 trusts the caller's declaration that the forces do not depend on the velocity.
- The uncertainty is carried linearly; anything non-linear is the caller's Monte Carlo with
  `myb_rng`. `myb_covariance_rtn` leaves out the rotation rate of the axes.
- An impulse is taken as given: there are no partials with respect to Δv.
- The full budget differences two propagations and has their floor of rounding; the linear one
  is of the first order.
- The variational equations of a whole [system of bodies](system-of-bodies.en.md) are not here.

## Data

None of its own. The time scales come from the caller's `myb_timesys` (the kernels kept in the
repository: `data/time/time-scales.tk`, `data/time/utc-tai.tk`, `data/time/tdb-tt-fb1990.tk`),
the providers from the caller's ephemerides and orientation data. `test_env`, `test_orbit`,
`test_uncert`, `test_transfer` and `test_events` need no downloads; `test_lageos` needs the
downloaded `de440.bsp`, `header.440`, `eopc04.1962-now`, `GOCO06s.gfc` and `fes2004_Cnm-Snm.dat`.
