# Restricted three-body problem: quasi-periodic tori

**English** | [Русский](cr3bp-tori.ru.md) · [All passports](index.en.md)

## What it is

Quasi-periodic orbits of the circular restricted three-body problem — Lissajous orbits about
planar and vertical Lyapunov orbits, quasi-halo orbits, quasi-periodic distant retrograde orbits
(the families measured) — as two-dimensional invariant tori about a periodic orbit whose
monodromy matrix has a pair of eigenvalues exp(±iρ) on the unit circle. A torus is found as an
invariant circle of the flow over a time T: the curve that the flow carries onto itself turned
by the angle ρ. Tori are followed along one-parameter families; their linear stability is the
spectrum of the linearised map of the circle, and the stable and unstable manifolds of a
hyperbolic torus can be started from it and carried to a surface of section; an elliptic torus
has a normal rotation σ and a normal bundle that the map turns by it. The periodic orbits themselves are
in [their passport](cr3bp-periodic.en.md).

## Source

- The invariant circle: G. Gómez, J. M. Mondelo, *Physica D* 157, 283 (2001); Z. P. Olikara,
  D. J. Scheeres, AAS 12-361 (2012).
- The stability of an invariant curve (the Floquet problem): A. Jorba, *Nonlinearity* 14, 943
  (2001).
- The handling of the two conditions the invariance equation holds too few (the flow keeps the
  Jacobi integral and the action ∮p·dq, so two multipliers at their gradients take their
  place; the matrix of the symplectic form in x, v is derived from the equations of motion), the
  grouping of the eigenvalues into circles and the rule that picks the normal rotation (the
  member of the circle whose eigenfunction has the least root-mean-square frequency) are the
  library's own.
- Implemented from the publications; no third-party code, no data files.

## Interface

Header [`myb/cr3bp.h`](../../include/myb/cr3bp.h):

- `myb_cr3bp_torus_from_orbit`, `myb_cr3bp_torus_continue`, `myb_cr3bp_torus_state`;
  structure `myb_cr3bp_torus`, enum `myb_cr3bp_torus_family` (`MYB_CR3BP_TORUS_JACOBI`,
  `MYB_CR3BP_TORUS_PERIOD`, `MYB_CR3BP_TORUS_ROTATION`);
- `myb_cr3bp_torus_stability`, `myb_cr3bp_torus_manifold`; structure
  `myb_cr3bp_torus_spectrum`, with the normal pair of an elliptic torus in `normal_rotation`,
  `normal_eigenvalue`, `normal_frequency`, `normal_residual`, `normal_re`, `normal_im`;
- `myb_cr3bp_propagate` with a `myb_cr3bp_section`;
- constants `MYB_CR3BP_MAX_CIRCLE` (35), `MYB_CR3BP_MAX_SPECTRUM`.

## Inputs and outputs

- The units and frame of the [periodic orbits](cr3bp-periodic.en.md).
- A torus: the circle at an odd number of points (5 … 35) at the angles 2πk/n, the time T of one
  turn along the orbit, the turn ρ (rad) across it in that time, the `defect` — the larger of
  the two multipliers left, which vanish at a solution and measure how far the points resolve
  the circle — and the tangent of its family in (states, T, ρ).
- `myb_cr3bp_torus_from_orbit`: the circle of the linear motion y0 + size (a cos θ − b sin θ),
  a + ib the eigenvector of the monodromy matrix for exp(iρ), scaled so that `size` is the root
  mean square distance of the circle from y0 (positions and velocities together, units of the
  problem); `pair` picks among the pairs on the unit circle in the order of ρ from 0 to π. The
  circle is exact to the order of size².
- `myb_cr3bp_torus_continue`: a step along the tangent (the root mean square displacement of the
  circle with the changes of T and ρ), corrected across it to `tolerance` (10⁻¹⁰) in at most
  20 Newton steps; a step of 0 corrects the torus where it is. The family keeps the mean Jacobi
  constant of the circle, T or ρ at their values on entry. Two phase conditions hold the points
  along the circle and along the flow.
- `myb_cr3bp_torus_state`: the point of the circle at any angle by its trigonometric
  interpolant.
- `myb_cr3bp_torus_stability`: the map A = R(−ρ) diag(Φ1 … Φn) of a displacement of the points
  (Φj the transition matrix of point j over T) and all its 6n eigenvalues, ordered by modulus;
  each eigenvalue of the torus is a circle of n eigenvalues of A that are turns of each other
  by ±ρ. A member counts as resolved when its circle holds more than half the points; the
  classes are trivial (the circles joined to 1), elliptic, hyperbolic and unresolved. Out also:
  the index (|λ| + 1/|λ|)/2, the real unstable and stable eigenvalues and their eigenfunctions at
  the points (unit root mean square, the mean x component non-negative), the residual
  |Av − λv|/(|λ||v|) and the share of their energy above half the highest frequency.
  `tolerance` (10⁻⁴) is the distance within which an eigenvalue counts as on the unit circle, at
  1, or a turn of another.
- The normal pair of an elliptic torus: two circles exp(i(σ − kρ)) and their conjugates, so σ is
  defined modulo ρ. Of the resolved elliptic eigenvalues in the upper half plane the one whose
  complex eigenvector (by inverse iteration on A in complex arithmetic, its Fourier coefficients
  c_k over k = −(n − 1)/2 … (n − 1)/2) has the least root-mean-square frequency
  (Σ k²|c_k|² / Σ |c_k|²)^(1/2) is taken: its eigenfunction tends to the constant eigenvector of
  the orbit as the torus shrinks, and σ to the angle of the orbit's elliptic pair. Out: σ in
  (0, π] (0: no elliptic pair), the eigenvalue, that frequency, the residual |Aψ − λψ|/(|λ||ψ|)
  and ψ at the points (unit root mean square, the largest component of its mean real and
  positive); with its conjugate ψ spans the normal bundle, Φ(θ) ψ(θ) = exp(iσ) ψ(θ + ρ).
- `myb_cr3bp_torus_from_orbit` takes the centre pair of a vertical Lyapunov orbit as of any
  other.
- `myb_cr3bp_torus_manifold`: the point of the circle at an angle displaced by ε (length units)
  along the position part of the eigenfunction there; the side +1 is that of growing x on the
  mean (towards the Moon from a torus about L1). The unstable manifold leaves forward in time,
  growing by |λ| in T while it is near; the stable one backward.

## Domain of validity

- Tori of two frequencies about a periodic orbit with a pair of multipliers on the unit circle.
- The points resolve the members of a circle up to some frequency only; beyond it eigenvalues
  leave the unit circle and are reported unresolved. The Jordan blocks of the trivial circles
  split by the square root of the error of A — 2.6·10⁻⁶ radially for the tori measured — so the
  tolerance must stay well above that and below the turn itself (|1 − exp(iρ)| > 10 tolerance).
- ε of a manifold well above the circle's own error (10⁻¹¹ for the tori measured): 10⁻⁶ … 10⁻⁴
  when the manifold is followed far.
- Refused with `MYB_E_DOMAIN`: μ, the period, the number of points, size ≤ 0, no such pair; an
  invalid torus or one without a tangent; an invalid tolerance or ρ too near 0; a spectrum of
  another number of points, or no hyperbolic pair, for a manifold.
- `MYB_E_NOCONV`: the eigenvector; a continuation step that does not converge (the torus is
  left as it was: a shorter step or more points); the eigenvalues, or a resolved circle off the
  unit circle with no real member (more points, or another tolerance). `MYB_E_NOMEM`.

## Accuracy

`test_cr3bp` ([tests/unit/test_cr3bp.c](../../tests/unit/test_cr3bp.c)), the circle carried by
another integrator (DOP853) than the solver's:

| What | Reference | Result |
|---|---|---|
| A distant retrograde orbit of the JPL catalogue, 21 points, size 10⁻³ | the circle at kρ after k turns | returns onto itself over 20 turns to 3.7·10⁻¹¹; multipliers 10⁻¹⁵; the Jacobi constant over the circle 4·10⁻¹⁵ |
| Four steps along its family to a size of 0.02 | the same | 1.3·10⁻¹⁰ |
| A step along the family and back | the starting circle | 1.1·10⁻¹¹ |
| The two pairs of that orbit | — | two tori: one in its plane, one out of it |
| A Lissajous torus about a small L1 Lyapunov orbit (38 km) | ρ of the linear motion, 0.1764459 | 0.1764423 |
| About the L1 Lyapunov orbit of 4500 km (unstable 2400-fold in a turn), 25 points | the circle after a turn | reaches 3617 km out of the plane; returns to 1.6·10⁻¹¹ |

`test_torus_stability` ([tests/unit/test_torus_stability.c](../../tests/unit/test_torus_stability.c)):

| What | Reference | Result |
|---|---|---|
| Lissajous tori about the L1 Lyapunov orbit (λ = 2446.324), 25 points, sizes 10⁻³ (384 km), 2·10⁻³, 4·10⁻³ | λ of the orbit | 2446.191, 2445.794, 2444.210: −5.4·10⁻⁵, −2.2·10⁻⁴, −8.6·10⁻⁴, ratios 3.997 and 3.988 (the size squared) |
| The same spectra | symplectic pairing; A v = λ v | all 50 hyperbolic members resolved; λ·(1/λ) = 1 to 4·10⁻¹⁰; A v = λ v to 4·10⁻⁹ |
| The eigenvectors | the map by DOP853 and finite differences | 3.7·10⁻⁶ |
| 35 points against 25 | λ | 1.7·10⁻¹⁴ |
| Quasi-halo tori of size 10⁻³ about two L2 halo orbits (indices 12.17 and 74.86) | λ of the orbit | 3.5·10⁻⁶ and 1.8·10⁻⁶ |
| A torus about a distant retrograde orbit, 21 points | — | all elliptic (42 of 126 eigenvalues; the rest trivial or unresolved); at a size of 0.02 spurious eigenvalues 5.6 % off the unit circle reported unresolved, the torus still elliptic |
| Manifolds, ε = 10⁻⁹ | growth by λ in a turn | to 3.4·10⁻⁶, forward on the unstable one, backward on the stable one |
| Manifolds of a Lissajous torus carried to x = 1 − μ, ε 10⁻⁶ … 10⁻⁴ | the manifolds of its Lyapunov orbit | they cross where those do, off by the size out of the plane and by its square in it (4·10⁻⁴ in y, vy at size 10⁻³, ratio 3.994); Jacobi constant kept to 4·10⁻¹⁵; the stable manifold the mirror image of the unstable one |
| The normal rotation: tori of 21 points on either pair of a distant retrograde orbit (pairs at 1.487246569 and 2.226143226 rad), sizes 10⁻³, 2·10⁻³, 4·10⁻³ | the angle of the other pair of the orbit | σ off by 6.86·10⁻⁶, 2.74·10⁻⁵, 1.09·10⁻⁴ — the size squared (ratios 3.997 and 3.987); on the other pair −6.4·10⁻⁷ (ratios 4.02 and 4.06) |
| The normal eigenfunction | the eigenvector of the orbit | 2.6·10⁻³ at the worst point (as the size: ratios 2.000 and 2.001), its mean 2.0·10⁻⁵ (as the square: 4.000) |
| The same | A ψ = λψ; \|λ\| = 1; between the points, by the transition matrices of interpolated points, Φ ψ = λ ψ(θ + ρ) | 2·10⁻¹⁵; 3·10⁻¹⁵; 7.5·10⁻¹³ of ψ |
| 25 and 35 points against 21 | σ | 5·10⁻¹⁵ |
| Tori of 25 points about the vertical Lyapunov orbits of L1 and L2 reaching 11 500 km out of the plane (λ 3231 and 1937, centre pairs 0.2078 and 0.2850 rad), sizes 10⁻³, 2·10⁻³, 4·10⁻³ | the circle after a turn (DOP853) | returns to 1.3·10⁻¹¹; the Jacobi constant one over the circle to 10⁻¹⁴ |
| The same | the centre pair and the unstable eigenvalue of the orbit | ρ off by −6.0·10⁻⁶ at 10⁻³, ratios 4.000 and 4.000 at both points; the unstable eigenvalue −2.0·10⁻⁵ and −1.7·10⁻⁵, ratios 4.000; every hyperbolic member resolved, no normal rotation |

## Cost

The spectrum costs as (6n)³: 35 points give a matrix of 210, 0.2 s (header). The normal rotation
of an elliptic torus adds one complex LU of A per candidate (about n of them): the spectrum takes
0.1 s with 21 points, 0.35 s with 35.

## Limitations

- Two-frequency tori only; the three-frequency tori about an elliptic torus are not computed.
- At most 35 points on a circle; a large torus needs more points than that may resolve, and
  its spurious eigenvalues are reported as unresolved, not as instability.
- A manifold state is the linear displacement along the eigenfunction; the circle's error grows
  with it.
- A continuation step that does not converge must be shortened by the caller.
- The normal rotation is a choice among the members of a circle (σ is defined modulo ρ): the
  rule of the smoothest eigenfunction continues the orbit's own angle as the torus shrinks, and
  other rules would give σ − kρ. A hyperbolic torus has none.
- The vertical Lyapunov orbits of the test come from the second-order guess of the figure eight
  corrected by the symmetric corrector; the tori about them are hyperbolic.

## Data

None. The tests start from orbits of the JPL catalogue whose states are written into them.
