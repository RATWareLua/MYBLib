# Restricted three-body problem: periodic orbits

**English** | [Русский](cr3bp-periodic.ru.md) · [All passports](index.en.md)

## What it is

The circular restricted three-body problem in its rotating frame: the equations of motion, the
Jacobi constant, the analytic Jacobian, the Lagrange points; periodic orbits symmetric to the
xz-plane or the x-axis corrected by single shooting on the half period, and periodic orbits of
any shape by multiple shooting, with their families followed by pseudo-arclength continuation;
the monodromy matrix, the stability index and the invariant manifolds; Richardson's third-order
halo orbits as first guesses; the pulsating frame of two real bodies; propagation to a time or
to a surface of section. The tori about these orbits are in [their own passport](cr3bp-tori.en.md).

## Source

- Equations and the Jacobi constant: V. Szebehely, *Theory of Orbits* (1967).
- Symmetric single shooting: K. C. Howell, *Celest. Mech.* 32, 53 (1984).
- Halo orbits of the third order: D. L. Richardson, *Celest. Mech.* 22, 241 (1980).
- Manifolds from the real eigenvectors of the monodromy matrix: W. S. Koon, M. W. Lo,
  J. E. Marsden, S. D. Ross, *Dynamical Systems, the Three-Body Problem and Space Mission
  Design* (2011).
- Pseudo-arclength continuation: H. B. Keller (1977); E. J. Doedel et al.
- Implemented from the publications; no third-party code. The NASA/JPL Three-Body Periodic
  Orbits catalogue serves as the reference of the tests.

## Interface

Header [`myb/cr3bp.h`](../../include/myb/cr3bp.h):

- `myb_cr3bp_acceleration`, `myb_cr3bp_jacobian`, `myb_cr3bp_jacobi`; terms for a force model
  `myb_cr3bp_term`, `myb_cr3bp_term_jacobian`; `myb_cr3bp_lagrange`;
- symmetric orbits: `myb_cr3bp_correct` (`myb_cr3bp_correction`, `myb_cr3bp_symmetry`),
  `myb_cr3bp_monodromy`, `myb_cr3bp_stability`, `myb_cr3bp_manifold`, `myb_cr3bp_richardson`;
- frames: `myb_cr3bp_from_ephemeris`, `myb_cr3bp_to_ephemeris`;
- orbits of any shape and families: `myb_cr3bp_orbit_from_state`, `myb_cr3bp_orbit_correct`,
  `myb_cr3bp_orbit_continue` (`myb_cr3bp_orbit`, `MYB_CR3BP_MAX_PATCHES` = 10);
- `myb_cr3bp_propagate` (`myb_cr3bp_section`).

## Inputs and outputs

- Units: the distance of the primaries, their total mass and their mean motion are 1; the
  origin is the barycentre, the larger primary at (−μ, 0, 0), the smaller at (1 − μ, 0, 0),
  μ = m2/(m1 + m2). Ω = (x² + y²)/2 + (1 − μ)/r1 + μ/r2; x″ − 2y′ = ∂Ω/∂x,
  y″ + 2x′ = ∂Ω/∂y, z″ = ∂Ω/∂z; C = 2Ω − |v|².
- Symmetry XZ: y0 = (x, 0, z, 0, vy, 0), at T/2 y = 0, vx = vz = 0 (halo, Lyapunov, distant
  retrograde); symmetry X: y0 = (x, 0, 0, 0, vy, vz), at T/2 y = z = vx = 0 (vertical, axial).
  The caller picks the component held fixed and the crossing of y = 0 that ends the half period;
  tolerance 10⁻¹², at most 30 iterations. The state transition matrix comes from the variational
  equations with the analytic Jacobian.
- The stability index is (|λmax| + 1/|λmax|)/2, as in the JPL catalogue (1 for a linearly
  stable orbit).
- A manifold state: the orbit's state at a fraction of the period displaced by ε (length units)
  along the position part of the eigenvector carried to that phase, on either side.
- Richardson: L1 or L2, the out-of-plane amplitude Az in units of the distance of the primaries,
  the northern or the southern branch → y0 at the xz-plane crossing and the period.
- The pulsating frame: origin at the barycentre of the two bodies, x along body 1 → body 2, z
  along their orbital angular momentum, lengths in units of their instantaneous distance, time
  in units of 1/n, n = √((gm1 + gm2)/d³); states in m and m/s relative to any common origin.
- Multiple shooting: the orbit cut into up to 10 parts equal in time; the unknowns are the states
  at the cuts and the period; the conditions are that the parts join — all components but the
  velocity component largest at the first cut, which the Jacobi integral gives — and a phase
  condition. Minimum-norm Newton steps (`myb_shoot`), tolerance 10⁻¹¹, at most 30 steps; no
  symmetry assumed. The tangent of the family is the null vector of the Jacobian of the
  conditions; the first one points towards a longer period.
- `myb_cr3bp_propagate`: IAS15 with ε = 10⁻¹² for a duration (negative: backward), stopping at
  the first crossing of y[axis] = value in a given direction; the crossing located on the dense
  output by the event locator of the integrators.

## Domain of validity

- μ in (0, 0.5].
- Refused with `MYB_E_DOMAIN`: y0 not of the symmetric form, invalid settings, the manifold of an
  orbit that is not unstable, an amplitude with no halo orbit, a number of patches outside
  1 … 10, a period ≤ 0, continuation of an orbit without a tangent, an invalid section.
- `MYB_E_NOCONV`: no convergence, or no crossing within 50 time units (`myb_cr3bp_correct`); a
  continuation step that does not converge leaves the orbit as it was (take a shorter step);
  more than 10⁶ steps in `myb_cr3bp_propagate`. `MYB_E_STEP_MIN` at a collision with a primary.

## Accuracy

`test_cr3bp` ([tests/unit/test_cr3bp.c](../../tests/unit/test_cr3bp.c)):

| What | Reference | Result |
|---|---|---|
| 56 orbits of 8 families of the JPL catalogue (Earth–Moon L1/L2 halo, L1/L2 Lyapunov, vertical, axial, distant retrograde; Sun–Earth L2 halo) | the catalogue | Jacobi constant 3.6·10⁻¹⁵; the catalogue state a fixed point of the corrector to 9·10⁻¹², period 6·10⁻¹² |
| Stability index, the same orbits | the catalogue | to 10⁻⁵, except two L2 Lyapunov orbits passing 500 km above the Moon (1.1·10⁻⁴) |
| Lagrange points | the catalogue | 10⁻¹⁴ |
| Richardson's halo as a first guess | the converged orbit | 1–2 % of γ; period to 0.15 % |
| Manifolds of the L1 Lyapunov orbit (λ = 2446) | growth over one period | by λ, forward on the unstable one, backward on the stable one (test bound 5 %) |
| Multiple shooting: the L1 Lyapunov orbit pushed off by 2·10⁻⁵ at every cut | a periodic orbit | closes again to 2·10⁻¹⁰ |
| 25 steps of continuation along its family | closure; the symmetric corrector | every member closes to 4·10⁻¹³; state 7·10⁻¹⁷, period 3·10⁻¹⁴ |
| The pulsating frame | round trips, L4 at rest | test bound 10⁻¹⁴ |

`test_torus_stability` ([tests/unit/test_torus_stability.c](../../tests/unit/test_torus_stability.c)):
`myb_cr3bp_propagate` on the L1 Lyapunov orbit of the catalogue crosses y = 0 at T/2 to 10⁻¹⁵
and at T to 2·10⁻¹², the state returning to 2·10⁻¹²; against DOP853 over 0.7 T, 9·10⁻¹³.

`test_spacecraft` ([tests/unit/test_spacecraft.c](../../tests/unit/test_spacecraft.c)): JWST
2023-01 … 03 from JPL Horizons, taken into the pulsating frame of the Sun and the Earth–Moon
barycentre, stays at x 1.0079 … 1.0111 (L2 at 1.0101), |y| < 0.005, |z| < 0.003, with the Jacobi
constant 3.000779 just below C(L2) = 3.000894, as libration orbits have.

## Cost

Not measured.

## Limitations

- Circular and restricted: the primaries move on circles, no eccentricity, no third body; in the
  pulsating frame the rotation of the orbital plane is neglected in the velocities.
- Single shooting needs the symmetry; multiple shooting takes at most 10 patches.
- The stability index is sensitive where the orbit passes close to a primary (the two L2
  Lyapunov orbits 500 km above the Moon: 1.1·10⁻⁴ from the catalogue).
- A manifold state is the linear displacement ε along the eigenvector.
- Richardson's expansion is a first guess only (1–2 % in γ).

## Data

None for the library. `test_cr3bp` reads `tests/data/cr3bp/jpl_periodic_orbits.txt` (NASA/JPL
Three-Body Periodic Orbits API, obtained 2026-11-10; the queries are recorded in the file);
`test_spacecraft` reads `tests/data/horizons/*.txt` and the downloaded `de440.bsp` and
`header.440`.
