# Initial orbit determination

**English** | [Русский](initial-orbit-determination.ru.md) · [All passports](index.en.md)

## What it is

A first orbit from a few measurements, in two-body motion about a centre of gravitational
parameter μ. From three positions the velocity at the middle one (Gibbs for points far apart,
Herrick–Gibbs for points close together); from three lines of sight seen from three places the
whole state (Gauss, Laplace, the double-r iteration, Gooding). Every method that starts from
sights ends in the same refinement: Newton's method on the six conditions that the object lies
on each sight, with the transition matrix of the two-body motion, so the orbit returned passes
through all three sights in exact Kepler motion. What comes out is a state to start an
[estimator](estimation.en.md) with, not an estimate: the measurements are taken as exact.

## Source

- Gibbs: J. W. Gibbs (1889), in the form of D. A. Vallado, *Fundamentals of Astrodynamics and
  Applications*, 4th ed., algorithm 54. Herrick–Gibbs: S. Herrick; Vallado, algorithm 55.
- Gauss: Vallado, algorithm 52; H. D. Curtis, *Orbital Mechanics for Engineering Students*,
  algorithms 5.5 and 5.6. The refinement by Newton's method on the six conditions is the
  library's own: the simple iteration of Curtis's scheme diverges in the geometries of the tests
  (a multiplier near −3).
- Laplace: P.-S. Laplace (1780); Vallado, algorithm 51; P. R. Escobal, *Methods of Orbit
  Determination* (1965), ch. 7.
- The double-r iteration: Escobal (1965), §7.3; Vallado, algorithm 53.
- Gooding: R. H. Gooding, *Celest. Mech. Dyn. Astron.* 66, 387 (1997), with Lambert's problem
  of `myb/lambert.h` (see [Kepler motion](kepler-regularization.en.md)).
- Implemented from the publications; no third-party code, no data files.

## Interface

Header [`myb/iod.h`](../../include/myb/iod.h):

- functions `myb_iod_gibbs`, `myb_iod_herrick_gibbs`, `myb_iod_gauss`, `myb_iod_laplace`,
  `myb_iod_double_r`, `myb_iod_gooding`;
- structure `myb_iod_solution`; constant `MYB_IOD_MAX_SOLUTIONS` (3).

## Inputs and outputs

- Positions in m, velocities in m/s, times in s, μ in m³/s², all vectors in one inertial frame
  centred on the body (any consistent units work: the examples of Vallado are run in km, km/s
  and km³/s²).
- Gibbs: three positions in the order of the motion → the velocity at the second; the times are
  not used. `out_of_plane` (may be NULL) receives the angle (rad) between r1 and the plane of r2
  and r3 — how far the three are from one orbit.
- Herrick–Gibbs: three positions at t1 < t2 < t3 → the velocity at the second.
- Gauss and Laplace: the times t[0] < t[1] < t[2], the positions of the three observers relative
  to the centre and three unit sights. The distance at the middle time is a root of
  x⁸ + a x⁶ + b x³ + c = 0; each root that gives positive ranges is a candidate. Up to three
  solutions, in the order of the roots, each with r and v at the middle time, the three ranges
  (all > 0) and the number of Newton steps.
- Double-r: guesses of the distances from the centre at the first two times (e.g. both the
  radius of the orbit expected). Gooding: guesses of the ranges along the first and the third
  sight, the number of whole turns (0 … 8) and, for one or more, the branch of Lambert's problem
  (−1 or +1); the sense of the motion is the one that misses the second sight less at the guess.
- The refinement stops when the step is at the rounding of the state (at most 50 steps).
  Candidates that do not settle, or settle behind an observer, are dropped; two roots that
  settle on one orbit are given once.

## Domain of validity

- Gibbs is exact for any three points of a conic (circle to hyperbola). Points a few degrees
  apart lose digits in the cross products: Herrick–Gibbs is the method there.
- Herrick–Gibbs is a fourth-order Taylor series: the error falls as the fourth power of the
  spacing — three parts in a million for a tenth of a radian between neighbours, seven in 10⁸
  for a twenty-fifth (header).
- Gauss and Laplace are methods of short arcs (the series of Gauss, the parabola through three
  sights); the refinement carries them over arcs the series alone cannot span.
- Double-r: less than a turn between the first and the second time and between the second and
  the third; not for nearly parabolic orbits. Gooding: long arcs and up to 8 whole turns.
- Three sights may lie on more than one orbit. Over two thirds of a turn, from guesses 5 % off,
  both the double-r iteration and Gooding settled on another orbit through the same sights; from
  1 % off, on the true one. Every solution passes through the three sights; which one it is
  depends on the guess.
- Refused with `MYB_E_DOMAIN`: a zero or non-finite vector, μ ≤ 0, positions on one line
  (Gibbs); times not increasing; a sight that is not a unit vector to 10⁻⁹; coplanar sights;
  sights that do not turn (L, L′, L″ in a plane: Laplace); distances or ranges not positive;
  turns outside 0 … 8 or a branch other than ±1.
- `MYB_E_NOCONV`: no candidate settled (Gauss, Laplace), no orbit from these guesses (double-r,
  Gooding).

## Accuracy

All in `test_iod` ([tests/unit/test_iod.c](../../tests/unit/test_iod.c)); the errors are the
largest component of r and v relative to its norm, against the exact Kepler orbit the sights
were made from (an observer on a sphere turning with the Earth):

| What | Reference | Result |
|---|---|---|
| Gibbs, e = 0, 0.2, 0.7, 0.95, 1.4 | the velocity of the orbit | 1.4·10⁻¹⁵ |
| Gibbs, Vallado's example 7-3 | (0, 5.531148, −5.191806) km/s | to the six printed digits |
| Herrick–Gibbs, spacing 1/40, 1/80, 1/160 of a period | the velocity of the orbit | error ratios 16.0 and 16.0 (fourth order) |
| Herrick–Gibbs, Vallado's example 7-4 | (−6.441557, 3.777559, −1.720567) km/s | to the six printed digits |
| Gauss: a low orbit over 4 and 12 min, e = 0.7 over 30 min, a navigation orbit over 2 h, geostationary height over 5 h | the orbit | 5·10⁻¹⁵ … 3·10⁻¹⁶ after 4–5 Newton steps (the series of Gauss alone: 3·10⁻⁵ … 5·10⁻³) |
| Laplace, the same five arcs | the orbit | 4·10⁻¹⁶ … 1.6·10⁻¹⁴ |
| Double-r and Gooding, arcs of a third, a half and two thirds of a turn | the orbit | 10⁻¹⁶ … 3·10⁻¹⁵ |
| Gooding over 1.85, 2.7 and 3.3 turns (1, 2, 3 revolutions of Lambert, both branches tried) | the orbit | 2·10⁻¹⁶ … 2·10⁻¹⁵ |

The test bounds for Gauss, Laplace, the double-r iteration and Gooding are 10⁻⁸ for the orbit
and 10⁻¹² rad for the distance of a solution from its sights.

## Cost

The refinement of Gauss's candidates takes 4–5 Newton steps on the arcs of the test (at most
50). No timings were measured.

## Limitations

- Two-body motion only: no perturbations, no light time, no aberration — a sight is a geometric
  direction at the time given.
- The measurements are exact: no weights and no noise model; a least-squares fit
  ([estimation](estimation.en.md)) is the next step.
- Gauss and Laplace may return up to three candidates; choosing among them (by a fourth sight,
  say) is the caller's.
- The double-r iteration and Gooding need guesses, and a guess decides which orbit through the
  sights is found.
- Herrick–Gibbs keeps the error of its series; Gibbs ignores the times.

## Data

None.
