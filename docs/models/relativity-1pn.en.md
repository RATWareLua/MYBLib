# Relativity at 1PN

**English** | [Русский](relativity-1pn.ru.md) · [All passports](index.en.md)

## What it is

General relativity at the first post-Newtonian order, with the PPN parameters β and γ: the
relativistic accelerations of a satellite about its central body (IERS 2010 eq. 10.12 —
Schwarzschild, Lense–Thirring, de Sitter) with their Jacobians and ready-made force terms; the
Einstein–Infeld–Hoffmann equations of N bodies as a level of the interaction graph, with the
1PN conserved quantities; the transformation between the global barycentric system (BCRS) and
the local system of a body (GCRS-like); the rate of proper time; and the Lense–Thirring drag of
a spinning body on the others in a system of bodies.

## Source

- IERS Conventions 2010, IERS Technical Note 36, chapter 10 (eqs. 10.6–10.12); IAU 2000
  Resolution B1.3 and M. Soffel et al., *AJ* 126, 2687 (2003); C. M. Will, *Theory and
  Experiment in Gravitational Physics*, 2nd ed. (2018).
- EIH: T. D. Moyer, DESCANSO Monograph 2 (2003), eq. 4-61 — the equations the JPL ephemerides
  integrate (W. M. Folkner et al., IPN PR 42-196, 2014, eq. 1). The 1PN conserved quantities and
  the independent check: the global PPN Lagrangian of C. Bernus et al., arXiv:2301.01186,
  app. B, its Euler–Lagrange equations derived with sympy and evaluated with 50 digits in mpmath
  (`tests/data/rel/eih_reference.txt`).
- Frame dragging: J. Lense, H. Thirring, 1918.
- Implemented from the sources; no third-party code.

## Interface

- Header [`myb/rel.h`](../../include/myb/rel.h): `myb_ppn`, `myb_ppn_validate`;
  `myb_rel_schwarzschild`, `myb_rel_lense_thirring`, `myb_rel_de_sitter`,
  `myb_rel_schwarzschild_jacobian`, `myb_rel_lense_thirring_jacobian`,
  `myb_rel_de_sitter_jacobian`; `myb_rel_local_frame`, `myb_rel_global_to_local`,
  `myb_rel_local_to_global`, `myb_rel_local_time_shift`; `myb_rel_proper_time_rate`.
- Header [`myb/dyn.h`](../../include/myb/dyn.h): `myb_nbody_set_relativity`,
  `MYB_INTERACTION_EIH` with `myb_nbody_set_interaction`, `myb_nbody_acceleration_xv`,
  `myb_nbody_pn_invariants`, `myb_nbody_pn_invariants_of`.
- Header [`myb/env.h`](../../include/myb/env.h): `myb_relativity_term` with
  `myb_term_schwarzschild`, `myb_term_lense_thirring`, `myb_term_de_sitter` and their
  `_jacobian` functions.
- Header [`myb/system.h`](../../include/myb/system.h): `myb_system_add_frame_dragging` — see
  [system of bodies](system-of-bodies.en.md).

Beyond 1PN: [post-Newtonian binaries](pn-binary.en.md). The N-body machinery:
[N bodies](nbody.en.md). TCB − TCG: [time scales](time-scales.en.md). Light time and the
Shapiro delay: [observables and light time](observables-light-time.en.md).

## Inputs and outputs

- `myb_ppn`: c (m/s, e.g. `MYB_CONST_C` of `data/constants/codata-2022.tk`), β, γ (1 and 1 in
  general relativity).
- The satellite: r, v relative to the central body (m, m/s), its gm (m³/s²); a in m/s²:
  - Schwarzschild: gm/(c² r³) {[2 (β + γ) gm/r − γ v²] r + 2 (1 + γ)(r·v) v};
  - Lense–Thirring: (1 + γ) gm/(c² r³) [(3/r²)(r × v)(r·j) + v × j], j the angular momentum of
    the central body per unit mass (m²/s) in the axes of r;
  - de Sitter: (1 + 2γ) [V × (−gm_ext R/(c² R³))] × v, R and V of the central body relative to
    the external one.
  The Jacobians give ∂a/∂r and ∂a/∂v row-major; the de Sitter term does not depend on r.
- The terms: the state relative to the body `central`; its spin is the body's `spin`
  (body-fixed) turned by its orientation; `external` for de Sitter.
- EIH (Moyer 4-61): for a receiver i and a source j the Newtonian term is joined by the 1PN
  terms in U and the Newtonian acceleration a of each body from its own sources (by the
  graph). The accelerations depend on the velocities: `myb_nbody_acceleration_xv`. The
  Newtonian and the 1PN parts are summed per receiver in a fixed order (Sum2): the bits do not
  depend on the number of threads.
- 1PN invariants, with GM in place of masses (every quantity carries one factor G): the energy
  in double-word arithmetic, the momentum, the relativistic barycentre (the origin of the BCRS
  of the JPL ephemerides) and its velocity.
- The local frame of a body E with barycentric xe, ve, ae at TCB t and w_ext the Newtonian
  potential of the other bodies at xe (m²/s²):
  X = r + c⁻² [½ ve (ve·r) + w_ext r + r (ae·r) − ½ ae r²], r = x − xe;
  T = t − c⁻² [A(t) + ve·r], A(t) = TCB − TCG at the geocentre from the time module. The
  coordinates are in the TCB- and TCG-compatible units of the Resolution (rescale for TDB, TT).
  `myb_rel_local_time_shift` gives −(ve·r)/c² (s).
- Proper time: dτ/dt = 1 − (v²/2 + U)/c², U positive, with the tidal terms of external bodies
  (eq. 10.6) or without (eq. 10.7).
- Frame dragging in a system: the Lense–Thirring acceleration on the point, j copied in the
  axes of the system (a pole that stays where it is), γ and c of the point masses' relativity;
  the body is pushed back by the ratio of the masses, so the momentum is kept.

## Domain of validity

- `myb_ppn_validate`: c > 0 and finite, β and γ finite; `MYB_E_DOMAIN` otherwise.
- The satellite functions and terms: `MYB_E_DOMAIN` for r = 0 or input that is not finite; de
  Sitter with the central body as the external one.
- EIH: `myb_nbody_acceleration` (without velocities) fails with `MYB_E_DOMAIN`, and the
  Wisdom–Holman integrators reject the system; `MYB_INTERACTION_EIH` without
  `myb_nbody_set_relativity`, and the 1PN invariants without it, are `MYB_E_DOMAIN`;
  separations within 10⁻¹⁵⁰ … 10¹⁵⁰ m.
- Frame dragging: `MYB_E_DOMAIN` without the relativity of the point masses, for the same body
  twice or a spin that is not finite; `MYB_E_OVERFLOW` above 64 pairs.
- The local frame: terms of O(c⁻⁴) are left out; the inverse is one fixed-point step.

## Accuracy

Tests: [`test_rel`](../../tests/unit/test_rel.c), [`test_env`](../../tests/unit/test_env.c),
[`test_de440_rel`](../../tests/unit/test_de440_rel.c),
[`test_lageos`](../../tests/unit/test_lageos.c), [`test_system`](../../tests/unit/test_system.c).

| What | Reference | Result |
|---|---|---|
| EIH, 48 cases of 2–5 bodies, (β, γ) = (1, 1), (1.3, 0.7), (0.8, 1.25), c = 10 and 1000 | Euler–Lagrange equations of the PPN Lagrangian, 50 digits (`test_rel`) | 2.7·10⁻¹⁰ of the 1PN terms (the rounding of doubles); every path and pool the same bits |
| 1PN energy, three bodies, v/c ~ 10⁻², 60 time units (IAS15) | conservation (`test_rel`) | 6·10⁻⁹ while the Newtonian energy changes by 1.6·10⁻⁵; 12 times smaller with c doubled |
| Mercury's perihelion, Sun and Mercury, 60 revolutions | 6πGM/(c² a (1 − e²)) (`test_rel`) | 42.9804″ a century, 5·10⁻⁹; the factor (2 + 2γ − β)/3 at (1.5, 0.7) and (0.6, 1.3) to 1.3·10⁻⁷ |
| LAGEOS node over 600 revolutions | 2GMj/(c² a³ (1 − e²)^{3/2}) and (3/2) GM☉ n/(c² a☉) (`test_rel`) | Lense–Thirring 30.631 mas/yr, de Sitter 19.188 mas/yr; the test bound is 10⁻³ |
| Perigee of an orbit with e = 0.3, 200 revolutions, two (β, γ) | the formula with (2 + 2γ − β)/3 (`test_rel`) | 10⁻³ (bound) |
| BCRS → GCRS at 6800 km without field and acceleration | the Lorentz transformation (`test_rel`) | 10⁻⁹ m; the correction on such an orbit is 66 mm |
| Proper time of a GPS clock | IERS 2010 §10.2 (`test_rel`) | faster than TT by 4.4647·10⁻¹⁰, the GPS frequency offset |
| Jacobians, β = 1.1, γ = 0.9 | central differences (`test_env`) | up to 7·10⁻¹¹ (the limit of the differences) |
| The terms of `myb/env.h`, a central body at rest and moving | the functions of `myb/rel.h` by hand (`test_env`) | bit for bit |
| DE440 from J2000, point masses with EIH, one year | DE440 (`test_de440_rel`) | Mercury 106 m (Newton 57 652 m), Venus 10 m, Mars 68 m; after 10 years Mercury 0.95 km (Newton 1823 km) |
| LAGEOS-1 over a week without relativity | the ILRS orbit (`test_lageos`) | 2.06 cm against 2.0 cm (radially 5.3 mm against 2.6) |
| Frame dragging in a system | `myb_rel_lense_thirring` (`test_system`) | 10⁻¹⁷ m/s² (bound) |

For the Sun the drag is 1.7·10⁻¹² of the attraction at Mercury; with it, the loss of mass of the
Sun and the small bodies the re-integration of DE440 holds Mercury to 2.1 cm, Venus to 1.0 cm
and Mars to 1.3 cm after five years ([system of bodies](system-of-bodies.en.md)).

## Cost

Not measured.

## Limitations

- The first order only; the higher orders of a pair are in [PN binaries](pn-binary.en.md).
- EIH is for point masses; figures and spins of the bodies enter through other levels (fields,
  frame dragging).
- Of the PPN parameters only β and γ.
- The local frame: O(c⁻⁴) left out, the inverse by one fixed-point step; T − t gives only the
  part that depends on the position, A(t) is the time module's.
- In the satellite terms j is given by the caller; in a system the spin is copied once (the
  pole does not move).
- Wisdom–Holman integrators cannot take the velocity-dependent 1PN forces.

## Data

None: c, GM, spins and states are the caller's. The tests use `header.440` and `de440.bsp`
(downloaded by `cmake/FetchData.cmake`) and `tests/data/rel/eih_reference.txt`, computed for
the tests with sympy and mpmath.
