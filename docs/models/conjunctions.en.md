# Close approaches and the probability of collision

**English** | [Русский](conjunctions.ru.md) · [All passports](index.en.md)

## What it is

How close two objects come, and how likely they are to touch. The closest approaches are the
minima of their distance — the instants where the range rate turns from negative to positive —
for two objects given by providers or for two orbits being propagated. The probability of a
short encounter is the integral of the Gaussian of the relative position over a disk in the
plane across the relative velocity; the covariances can be carried to every approach by the
variational equations of the orbits themselves. Where the two move slowly against each other,
the probability of a long encounter counts the entries of the relative state into the sphere of
the summed radii, integrated in time.

## Source

- The short encounter: J. L. Foster, H. S. Estes, NASA/JSC-25898 (1992); S. Alfano,
  *J. Astronaut. Sci.* 53, 103 (2005). The integral is the library's own reduction to one
  integral of a difference of error functions.
- F. K. Chan, *Spacecraft Collision Probability* (2008), ch. 4 — summed with positive terms only
  (as usually written, the series loses every digit at large misses).
- The largest probability over the size of the covariance: S. Alfano, AAS 03-548.
- The long encounter: V. T. Coppola, AAS 12-247 (2012); D. T. Hall et al., AAS 18-252 (2018).
- Implemented from the publications; no third-party code, no data files.

## Interface

Header [`myb/conjunction.h`](../../include/myb/conjunction.h):

- approaches: `myb_approach_find` (providers `myb_state_function` of
  [`myb/kin.h`](../../include/myb/kin.h)), `myb_approach_orbits` (two `myb_orbit` of
  [`myb/orbit.h`](../../include/myb/orbit.h)); structure `myb_approach`;
- the short encounter: `myb_encounter_plane`, `myb_collision_probability`,
  `myb_collision_probability_chan`, `myb_collision_probability_max`; structure `myb_encounter`;
- with carried covariances: `myb_conjunction_find`; structure `myb_conjunction`;
- the long encounter: `myb_collision_rate`, `myb_collision_probability_long` (callback
  `myb_relative_gaussian_function`), `myb_conjunction_long`.

## Inputs and outputs

- An approach: its time (s), the distance (m), the position (m) and velocity (m/s) of the second
  object relative to the first. Providers are called at `myb_time_from_seconds(t)` and must give
  both objects in one frame; two orbits must be started and standing at one time, and both are
  left at t1.
- The search steps through [t0, t1] in equal steps; an approach is located to the resolution of
  the time (providers) or on the cubic through the relative states at the ends of the step
  (orbits: within step⁴/384 of the fourth derivative of the relative position — a quarter of a
  millimetre for 10 s in a low orbit, a third of a metre for 60 s).
- An encounter: the miss vector (m) and the standard deviations σ0 ≥ σ1 > 0 (m) in the principal
  axes of the covariance projected onto the plane across the relative velocity, and the radius —
  the sum of the radii of the two (m). Its input covariance is that of the relative position,
  3 × 3 row-major in m², the sum of the two covariances, in the frame of the vectors.
- `myb_collision_probability`: one integral across the disk of a difference of two error
  functions (through erfc in the tails, with no subtraction of close numbers), Gauss–Legendre
  panels in an angle that removes the roots at the rim, doubled until the value settles to
  10⁻¹³ relative (or 10⁻³⁰⁰ absolute).
- `myb_conjunction_find`: the 6 × 6 covariances of the two states (m, m/s) at the start, both
  orbits created with `variational`; at every approach the covariance of the relative position
  is Φ_r P Φ_rᵀ of the first plus that of the second, Φ_r the rows of the positions taken on the
  same cubic as the states.
- The long encounter: the rate of entry (1/s) is the flux of the probability density of the
  relative state through the sphere of radius R, the velocity conditioned on the place:
  R² ∫ p_r(R n) E[(−n·v)⁺ | r = R n] dΩ. The sphere is integrated by Gauss–Legendre panels of
  8 nodes in the cosine of the latitude and equal steps in longitude, doubled until the rate
  settles to 10⁻¹². Its integral over time is the expected number of entries: for a function,
  Gauss–Kronrod panels of 15 points no longer than `step`, halved until the estimates agree to
  10⁻¹⁰ of the total (at most 4000 pieces); for two orbits, the same panels not halved (the
  orbits go forward only), with the mean the difference of the states and the covariance
  Φ₁ P₁ Φ₁ᵀ + Φ₂ P₂ Φ₂ᵀ. `error` receives the sum of the disagreements of the two rules.

## Domain of validity

- Two approaches closer together than a step may be missed.
- The short encounter assumes a straight-line pass, quick against the time in which the
  uncertainty changes, so that only the plane across the relative velocity matters.
- Chan's series is exact for equal sigmas; for unequal ones it is an approximation, good while
  the disk is small against both. It returns NaN for an invalid encounter and where
  R²/(2 σx σy) or half the squared miss in sigmas reaches 700.
- The long encounter: the expected number of entries is the probability while it is small (a
  second entry needs an exit first); being inside the sphere at the start is not counted. A
  straight line with a certain velocity and a constant covariance gives exactly the short
  encounter.
- Refused with `MYB_E_DOMAIN`: t1 ≤ t0, a step ≤ 0, non-finite input, orbits at different
  times or t1 not ahead of them, a zero relative velocity, a covariance that is not symmetric
  or not positive in the plane (or, for the rate, a position covariance not positive definite),
  a radius ≤ 0.
- `MYB_E_NOCONV`: the disk integral does not settle; the sphere does not settle (a radius far
  beyond the uncertainty); the time pieces exhausted.

## Accuracy

All in `test_conjunction` ([tests/unit/test_conjunction.c](../../tests/unit/test_conjunction.c)):

| What | Reference | Result |
|---|---|---|
| Two circles crossing at a node; coplanar circles (one approach a synodic period) | closed forms | agree (the times and the distance of the coplanar case to the test bound of 10⁻⁶) |
| Two propagated orbits in steps of 10 s | the same pair as providers | 2·10⁻¹¹ s and 10⁻⁸ m |
| A centred round Gaussian, R/σ from 0.01 to 8 | 1 − exp(−R²/2σ²) | 2·10⁻¹⁶ |
| An offset round Gaussian, misses 0.5 … 30σ (probabilities down to 8·10⁻¹⁹⁰) | the series of Rice (Chan's for equal σ) | 4·10⁻¹⁴ |
| σ = 3 and 0.5 | a sum over a grid of 4000 × 4000 | 2·10⁻⁶ |
| The plane of an encounter from the vectors (axes 400, 100, 50 m) | the construction | σ and the miss to the test bound of 10⁻⁹ m |
| The largest probability over the size of the covariance | R²/(e d²) at σ = d/√2 | 10⁻⁶ |
| Covariances carried by the orbits, 4 approaches, steps of 20 s | the closed transition matrix of the two-body motion | the covariance 2.5·10⁻⁹, the probability 2.6·10⁻⁹ |
| Long encounter: a straight line with a certain velocity | the probability of the short encounter | 3.7·10⁻¹⁵ |
| Long encounter: two orbits crossing at 8.4 km/s | their short encounter (`myb_conjunction_find`) | 4.5·10⁻⁸ |
| Long encounter: 5 cm/s, the velocity as uncertain as itself | Monte Carlo of 2·10⁶ lines: 1.572·10⁻² ± 8.8·10⁻⁵ | 1.584·10⁻² (1.3σ); the short formula gives 8.65·10⁻³ |

## Cost

Not measured.

## Limitations

- Pairs only: no screening of a catalogue.
- The objects are spheres of the summed radius: no shapes, no attitude.
- The uncertainty is Gaussian and carried linearly (the transition matrices); no non-Gaussian
  or non-linear propagation of the covariance.
- `myb_collision_probability_max` assumes the covariance known in shape only — every sigma
  multiplied by one unknown factor.
- In `myb_conjunction_long` the panels are not refined: the step must be well below the time it
  takes to cross the uncertainty, and `error` says where it is not.

## Data

None.
