# Symplectic integrators

**English** | [Русский](symplectic-integrators.ru.md) · [All passports](index.en.md)

## What it is

Fixed-step integrators whose map is symplectic and time-reversible: the energy error stays
bounded over any number of steps instead of drifting, and the phase error grows linearly; the
price is a fixed step without error control. Three families: compositions of the
kick–drift–kick leapfrog for x'' = a(t, x) (Yoshida's triple jump of orders 2, 4, 6, 8 and two
schemes of orders 6 and 8 with smaller error constants); the Wisdom–Holman splitting of an
[N-body system](nbody.en.md) into Kepler motions in Jacobi coordinates and interaction kicks
(WHFast: DKD = SABA1, SABA2–4, SBAB1–4, symplectic correctors of odd orders 3…17); and, briefly,
the splitting of the rotation of a rigid body, torque kick — exact free rotation — kick.

## Source

- Compositions: H. Yoshida, *Phys. Lett. A* 150, 262 (1990) — the triple jump and the 7-stage
  solution A of order 6; the 15-stage method p8 s15 of order 8 from E. Hairer, C. Lubich,
  G. Wanner, *Geometric Numerical Integration*, 2nd ed., 2006, V.3, after W. Kahan, R.-C. Li,
  *Math. Comp.* 66, 1089 (1997) and R. I. McLachlan, *SIAM J. Sci. Comput.* 16, 151 (1995).
- Wisdom–Holman: J. Wisdom, M. Holman, *AJ* 102, 1528 (1991); WHFast — H. Rein, D. Tamayo,
  *MNRAS* 452, 376 (2015); SABA and SBAB — J. Laskar, P. Robutel, *Celest. Mech. Dyn. Astron.*
  80, 39 (2001); symplectic correctors — J. Wisdom, M. Holman, J. Touma, *Fields Inst. Commun.*
  10, 217 (1996).
- Rotation: J. Touma, J. Wisdom, *AJ* 107, 1189 (1994); Celledoni et al. 2008.
- Implemented from the papers; no third-party code. The triple-jump weights are products of
  z1 = 1/(2 − 2^(1/(p+1))), z0 = 1 − 2 z1 computed in MPFR at 512 bits and rounded once; the 7-
  and 15-stage weights are the published values refined by Newton's method on the order
  conditions at 512 bits; SABA/SBAB coefficients come from closed forms. The correctors were
  derived in the library (the coefficients a_i = i·√(7/40), the choice of Wisdom et al., and b
  from a linear system solved in MPFR); for order 3 the derivation gives
  b₁ = 1/(96·√(7/40)) = 0.0249005960277998675, the published value.

## Interface

Header [`myb/integ.h`](../../include/myb/integ.h) — compositions:

- functions `myb_symplectic_create`, `myb_symplectic_destroy`, `myb_symplectic_start`,
  `myb_symplectic_step`, `myb_symplectic_integrate`, `myb_symplectic_time`,
  `myb_symplectic_state`, `myb_symplectic_statistics`;
- structure `myb_symplectic_stats`; callback `myb_position_accel_function`;
- the orders 2, 4, 6, 8 and the constants `MYB_SYMPLECTIC_6_S7`, `MYB_SYMPLECTIC_8_S15`.

Header [`myb/dyn.h`](../../include/myb/dyn.h) — Wisdom–Holman:

- functions `myb_wh_create`, `myb_wh_destroy`, `myb_wh_start`, `myb_wh_step`,
  `myb_wh_integrate`, `myb_wh_time`, `myb_wh_state`, `myb_wh_dense`, `myb_wh_statistics`;
- structures `myb_wh_options`, `myb_wh_stats`; enumeration `myb_wh_scheme` (`MYB_WH_DKD`,
  `MYB_WH_SABA2` … `MYB_WH_SABA4`, `MYB_WH_SBAB1` … `MYB_WH_SBAB4`);
- the system `myb_nbody` (see [the N-body system](nbody.en.md)).

Header [`myb/attitude.h`](../../include/myb/attitude.h) — rotation: `myb_rotation_split`,
callback `myb_torque_function` (the free rotation it builds on and the coupled rotation are in
[figures and their rotation](figure-figure-rotation.en.md)).

## Inputs and outputs

- Compositions: the state is y = (x[0..n−1], v[0..n−1]) for n position components, in any
  consistent units. The acceleration a(t, x) may depend on the time (the time is a coordinate
  of the extended phase space, so the map stays symplectic) but not on the velocity.
  `myb_nbody_acceleration` has the signature of the callback, so an N-body system is integrated
  directly. `myb_symplectic_step` takes a step of either sign; `myb_symplectic_integrate` makes
  equal steps of at most h_max to exactly t_end. Positions, velocities and time are updated
  with compensated summation.
- Wisdom–Holman: inertial positions and velocities, 3 per body (m, m/s), of a `myb_nbody`.
  Body 0 is the central body (GM > 0); the massive bodies follow in order, best from the inside
  out; test bodies are referred to the centre of mass of all massive bodies. Jacobi body i ≥ 1
  follows a Kepler problem with GM = η_i = GM_0 + … + GM_i (the drift is the universal-variable
  propagator, all conics); the kick is the N-body acceleration of the graph less those Kepler
  terms, which are subtracted whatever the graph says. The Jacobi state is accumulated with
  compensated summation.
- The step h has a sign (the forward direction); the time of the k-th node is t0 + k h,
  computed from k. `myb_wh_state` returns the corrected state at the node; `myb_wh_dense` makes
  one step of length t − node of the same scheme and corrector from the corrected node state on
  scratch memory, leaving the integration unchanged.
- A corrector is applied only to the states reported; the kernel steps uncorrected.
- Rotation: principal moments, the attitude quaternion q and the body angular momentum l in
  principal axes, advanced by `steps` steps of h from t0; the torque depends on the time and the
  attitude, not on l (a potential).

## Domain of validity

- Compositions: `MYB_E_DOMAIN` for another order or n = 0 (`myb_symplectic_create`), h = 0
  (`myb_symplectic_step`), h_max ≤ 0 (`myb_symplectic_integrate`). A status of the
  acceleration is returned unchanged and leaves the state in the middle of the step: start
  again.
- Wisdom–Holman: systems dominated by one central body, with a small ratio e of the
  perturbation to the Kepler part. Close encounters between bodies break the splitting: use
  [IAS15](ias15.en.md).
- `myb_wh_create`: `MYB_E_DOMAIN` for h zero or not finite, an unknown scheme, a corrector
  other than 0 or an odd order 3…17, a corrector with a scheme other than DKD, GM_0 not > 0;
  `MYB_E_OVERFLOW`, `MYB_E_NOMEM`.
- `myb_wh_start`: `MYB_E_DOMAIN` for non-finite input, coincident bodies, a system with
  velocity-dependent forces (the 1PN terms of `myb_nbody_set_relativity`) or with an
  approximation (`myb_nbody_set_approximation`: its error is relative to the whole force, the
  central body's included, which the splitting subtracts).
- `myb_wh_integrate`: `MYB_E_DOMAIN` for a non-finite t_end or more than 2^53 steps;
  `myb_wh_dense`: `MYB_E_RANGE` for |t − node| > |h|.
- Rotation: `MYB_E_DOMAIN` for an order other than 2, 4, 6, 8, steps < 0, invalid moments, a
  free rotation on the separatrix; statuses of the torque are returned.

## Accuracy

Compositions, in `test_symplectic` ([tests/unit/test_symplectic.c](../../tests/unit/test_symplectic.c)),
`test_mpfr_yoshida` ([tests/mpfr/test_mpfr_yoshida.c](../../tests/mpfr/test_mpfr_yoshida.c)) and
`test_mpfr_compositions` ([tests/mpfr/test_mpfr_compositions.c](../../tests/mpfr/test_mpfr_compositions.c)):

| What | Reference | Result |
|---|---|---|
| Triple-jump weights | products of z1, z0 in MPFR, 512 bits | bit for bit |
| 7- and 15-stage weights | order conditions from the Baker–Campbell–Hausdorff formula (55 words to weight 8), refined at 512 bits | conditions to 2⁻⁵⁰⁹ and 2⁻⁵¹⁰; doubles bit for bit |
| Observed order, one Kepler orbit of e = 0.3 | `myb_kepler_propagate` | 2.00 / 4.00 / 6.05 / 7.98 (triple jump); 6.00 and 8.00 (7 and 15 stages) |
| Equal work, orders 6 and 8 | the triple jump of the same order | its error is 1.75 and 196 times larger |
| Energy, order 6, e = 0.5, 200 steps per orbit, 5000 orbits | initial energy | the same largest error over orbits 1–100 and 4901–5000 (no drift) |
| 10 orbits of e = 0.7 forward and back, order 8 | the start | 1.7·10⁻¹² |
| Driven oscillator x'' = −x + cos 2t, order 6, h = 0.01, to t = 10 | closed-form solution | within the test bound 10⁻¹² |

Wisdom–Holman, in `test_whfast` ([tests/unit/test_whfast.c](../../tests/unit/test_whfast.c))
and `test_mpfr_wh` ([tests/mpfr/test_mpfr_wh.c](../../tests/mpfr/test_mpfr_wh.c)); positions
in units of the inner semi-major axis:

| What | Reference | Result |
|---|---|---|
| SABA/SBAB coefficients, corrector a_i and b | closed forms and the corrector system in MPFR, 256 bits | nearest doubles; quadratures exact to degree 2n − 1 |
| Two bodies (a planet and a test body), 2000 DKD steps, 100 orbits | `myb_kepler_propagate` | 1·10⁻¹¹ a (rounding of the drifts) |
| Observed orders, two planets with GM ratio ε = 10⁻⁷ | IAS15 | DKD 2.00; correctors 3/5 — 4.00/5.81; SABA2/SBAB2 — 3.93/3.96; orders 6–8 reach rounding at 16 steps per orbit |
| Jupiter–Saturn, 100 Saturn orbits, 50 steps per Jupiter orbit | IAS15 | DKD 9·10⁻⁵; corrector 17 — 1·10⁻⁵ (energy 1.9·10⁻⁷ → 1.9·10⁻¹⁰); SABA4/SBAB4 — 9·10⁻⁷ |
| Jupiter–Saturn, corrector 17, 5000 Jupiter orbits | initial energy | bounded, ≤ 7·10⁻¹⁰ |
| Earth–Mars–Jupiter and a test body, corrector 11, 3000 steps forward and back | the start | 1.3·10⁻¹³ |
| Dense output at 0.37 h inside a step | IAS15 | the accuracy of the nodes |
| 161 bodies, SABA2, with a pool of 4 threads | the same without a pool | the same bits |

The corrector removes the O(e h²) part of the energy error; the e² h² part of the position
error, which no corrector touches, remains (the Jupiter–Saturn row).

Rotation, in `test_attitude` ([tests/unit/test_attitude.c](../../tests/unit/test_attitude.c)):
observed orders 2.00 / 3.99 / 5.98 / 7.95 for the gravity-gradient torque of a point mass,
against DOP853 at rtol 10⁻¹⁴; without torque the result is the exact free rotation (test bound
10⁻¹² rad).

## Cost

Evaluation counts (no timings measured):

- compositions: 1, 3, 9, 27 acceleration evaluations per step for orders 2, 4, 6, 8 of the
  triple jump, 7 and 15 for `MYB_SYMPLECTIC_6_S7` and `MYB_SYMPLECTIC_8_S15`; for one accuracy
  the latter need 10 % and 50 % less work than the triple jump of the same order (header);
- DKD: one N-body evaluation and two drifts per step (`test_whfast`); SBAB reuses the
  acceleration of the last kick in the first kick of the next step; a corrector of order k costs
  2 (k − 1) extra kicks per reported state only.

## Limitations

- A fixed step, no error control, no events; the compositions have no dense output.
- Forces that depend on the velocity are not covered (drag, the 1PN terms): use IAS15.
- Wisdom–Holman needs a dominant central body; close encounters break the splitting. It
  refuses a system with 1PN terms or with the fast multipole approximation.
- The kick subtracts the Kepler terms whatever the interaction graph says; the centre of mass
  moves freely (zero acceleration) only for a mutual graph.
- Correctors exist only for DKD and do not reduce the e² h² error.
- `myb_rotation_split` uses only the triple-jump weights (orders 2, 4, 6, 8), not the 7- and
  15-stage schemes, and takes torques that do not depend on the angular momentum.

## Data

None: the coefficients are constants of the code (`src/integ/symplectic_coefficients.h`,
`src/dyn/wh_coefficients.h`).
