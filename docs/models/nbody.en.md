# The N-body system

**English** | [Русский](nbody.ru.md) · [All passports](index.en.md)

## What it is

A system of point masses given by GM: massive bodies that attract each other and test bodies
that feel the massive ones and attract nothing. The accelerations are exact direct sums over all
pairs with compensated summation (Sum2), the same bits whatever the number of threads and the
vector width; an interaction graph sets the level of every directed pair (none, Newtonian, or
Newtonian with the first post-Newtonian EIH terms); the conserved quantities are computed in
double-double. As an approximate level, off by default, the Newtonian sums can come from a fast
multipole method. The system is what IAS15, RKN12(10), Gauss–Jackson, the symplectic
compositions and Wisdom–Holman integrate, and the point-mass core of
[the system of bodies](system-of-bodies.en.md).

## Source

- The sum a_i = Σ_j GM_j (x_j − x_i) / |x_j − x_i|³, terms in increasing j added with Sum2:
  T. Ogita, S. M. Rump, S. Oishi, *SIAM J. Sci. Comput.* 26, 1955 (2005).
- The EIH level: T. D. Moyer, DESCANSO Monograph 2 (2003), eq. 4-61, the equations of the JPL
  ephemerides (W. M. Folkner et al., IPN PR 42-196, 2014, eq. 1), with the PPN parameters β, γ;
  the 1PN integrals from C. Bernus et al., arXiv:2301.01186, app. B. The physics and its checks
  are in [1PN relativity](relativity-1pn.en.md).
- The fast multipole method: L. Greengard, V. Rokhlin, *J. Comput. Phys.* 73, 325 (1987); the
  dual tree walk of W. Dehnen, *J. Comput. Phys.* 179, 27 (2002) and *Comput. Astrophys.
  Cosmol.* 1, 1 (2014); solid-harmonic expansions in the normalisation of White and Head-Gordon
  and of Dehnen, built by recurrences without libm.
- Implemented from the papers; no third-party code, no data files.

## Interface

Header [`myb/dyn.h`](../../include/myb/dyn.h):

- functions `myb_nbody_create`, `myb_nbody_destroy`, `myb_nbody_count`, `myb_nbody_massive`,
  `myb_nbody_gm`, `myb_nbody_set_pool`, `myb_nbody_set_interaction`,
  `myb_nbody_set_relativity`, `myb_nbody_acceleration`, `myb_nbody_acceleration_xv`,
  `myb_nbody_invariants_of`, `myb_nbody_pn_invariants_of`, `myb_nbody_set_approximation`,
  `myb_nbody_approximation`;
- structures `myb_nbody_invariants`, `myb_nbody_pn_invariants`; enumeration `myb_interaction`
  (`MYB_INTERACTION_NONE`, `MYB_INTERACTION_NEWTON`, `MYB_INTERACTION_EIH`); the constant
  `MYB_NBODY_FMM_MAX_ORDER`;
- the PPN parameters `myb_ppn` of [`myb/rel.h`](../../include/myb/rel.h).

Header [`myb/pool.h`](../../include/myb/pool.h): `myb_pool_create`, `myb_pool_destroy`,
`myb_pool_threads`.

## Inputs and outputs

- GM in m³ s⁻²; states are arrays of 3 values per body (x, y, z of body 0, then body 1, …) in
  m, m/s, m/s², in an inertial frame. Bodies 0 … massive − 1 attract; the test bodies follow.
- `myb_nbody_acceleration(nbody, t, x, a)` has the signature of `myb_position_accel_function`
  (RKN12, Gauss–Jackson, compositions); `myb_nbody_acceleration_xv` that of
  `myb_accel_function` (IAS15, Gauss–Jackson). t is not used by Newtonian gravity.
- The graph: a level per directed edge receiver ← source; the source must be massive and differ
  from the receiver. The default is `MYB_INTERACTION_NEWTON`, or `MYB_INTERACTION_EIH` for every
  pair once `myb_nbody_set_relativity` has been called; pairs set explicitly keep their level.
  The potentials U_i and the Newtonian accelerations a_j inside the 1PN terms follow the graph.
- Invariants of the massive bodies, every one carrying a factor G (GM in place of masses):
  G E (m⁵ s⁻⁴) as a double-double (energy to ~10⁻³⁰ relative), G P (m⁴ s⁻³), G L (m⁵ s⁻³),
  the centre of mass and its velocity. Every pair is counted (the graph is not consulted). The
  1PN versions give the 1PN energy, momentum and the relativistic barycentre (the origin of the
  BCRS of the JPL ephemerides).
- Threads (`myb_nbody_set_pool`; NULL — the calling thread only): every output element is
  computed by one thread with the same operations in the same order, so the bits are identical
  for any number of threads; the per-receiver path takes four receivers per vector block,
  without FMA. No memory is allocated by an evaluation of exact sums. The FMM runs on the same
  pool from 2048 bodies on, with the bits of its evaluation without a pool.
- The approximation: order p (0 — exact sums, or 1 … `MYB_NBODY_FMM_MAX_ORDER` = 20) and the
  opening angle θ of the criterion r_A + r_B < θ |z_A − z_B|.

## Domain of validity

- Separations within 10⁻¹⁵⁰ … 10¹⁵⁰ m. The plan targets up to ~10⁴ massive bodies for the
  exact sums.
- `myb_nbody_create`: `MYB_E_DOMAIN` for massive = 0 or a GM negative or not finite;
  `MYB_E_OVERFLOW`, `MYB_E_NOMEM`.
- `myb_nbody_set_interaction`: `MYB_E_DOMAIN` for bad indices or level, or
  `MYB_INTERACTION_EIH` without `myb_nbody_set_relativity`; `myb_nbody_set_relativity`:
  `MYB_E_DOMAIN` for invalid parameters.
- `myb_nbody_acceleration`, `myb_nbody_acceleration_xv`: `MYB_E_DOMAIN` for two interacting
  bodies at the same place or a non-finite position or result (a is then unspecified);
  `myb_nbody_acceleration` also fails with `MYB_E_DOMAIN` once the 1PN terms are set (they need
  the velocities). With an approximation the tree may grow: `MYB_E_NOMEM`.
- Invariants: `MYB_E_DOMAIN` for coincident massive bodies, a zero sum of GM, non-finite input;
  the 1PN invariants also without relativity parameters.
- `myb_nbody_set_approximation`: `MYB_E_DOMAIN` for an order outside 0 … 20 or θ outside
  (0, 1); `MYB_E_NOMEM`, `MYB_E_OVERFLOW` (the former setting is kept).
- `myb_pool_create`: 0 (one per logical processor) … 1024 threads, else `MYB_E_DOMAIN`;
  `MYB_E_NOMEM`, `MYB_E_IO` (the system refused a thread).

## Accuracy

Exact sums and the classical solutions, in `test_nbody`
([tests/unit/test_nbody.c](../../tests/unit/test_nbody.c)):

| What | Reference | Result |
|---|---|---|
| Accelerations, 400 + 50 bodies, GM 10⁵ … 10²⁰ | every term and the sum in double-double | ≤ 3.5·10⁻¹⁶ of Σ\|terms\| (plain summation of the same terms 2.1·10⁻¹⁵) |
| Symmetric and per-receiver paths, vector blocks, 1–17 massive + 0–6 test bodies, pools of 1–13 threads | each other; the library's Sum2 written out by hand for the graph | bit for bit |
| A whole IAS15 integration, 300 bodies, with a pool of 6 threads | the same without a pool | the same bits |
| Invariants of a circular pair | closed form | within the test bound 10⁻¹⁵ |
| Lagrange's equilateral solution (GM 1, 10⁻³, 10⁻⁶), 20 periods with IAS15 | rigid rotation | 1·10⁻¹³ of the side |
| Figure eight, start refined by Gauss–Newton shooting from the 6 digits of Šuvakov & Dmitrašinović (2013) | one period | residual 1.9·10⁻¹⁵; agrees with Simó's start to 10⁻⁹ (period 3·10⁻⁸) |
| Figure eight, 10 periods | the start; initial energy | 2·10⁻¹⁴; 3·10⁻¹⁶ |
| Burrau's Pythagorean problem to t = 70 | Szebehely & Peters (1967) | the same outcome (GM 4 and 5 bound, GM 3 escapes); closest approach 4.1·10⁻⁴ at t = 15.83 |
| The same, IAS15 against RKN12(10) | each other | 4·10⁻¹² at t = 10 (chaos amplifies rounding to 2·10⁻² by t = 70); energy to 10⁻¹⁰ |

The EIH level, in `test_rel` ([tests/unit/test_rel.c](../../tests/unit/test_rel.c)): against the
Euler–Lagrange equations of the global PPN Lagrangian solved with 50 digits, 48 cases of 2–5
bodies and several (β, γ) — 2.7·10⁻¹⁰ of the size of the 1PN terms (double rounding), the same
bits on every path and pool; the 1PN energy of a triple system with v/c ~ 10⁻² over 60 time
units changes by 6·10⁻⁹ while the Newtonian one changes by 1.6·10⁻⁵. More in
[1PN relativity](relativity-1pn.en.md).

The fast multipole method, in `test_fmm` ([tests/unit/test_fmm.c](../../tests/unit/test_fmm.c)):
relative force error, rms over 10⁴ bodies of |a − a_direct| / |a_direct| (the largest over the
bodies is 10 to 50 times more):

| p | cube θ 0.3 | cube θ 0.5 | cube θ 0.7 | Plummer θ 0.3 | Plummer θ 0.5 | Plummer θ 0.7 |
|---|---|---|---|---|---|---|
| 4 | 2.5·10⁻⁵ | 4.3·10⁻⁴ | 2.1·10⁻³ | 1.3·10⁻⁴ | 2.7·10⁻³ | 1.1·10⁻² |
| 6 | 8.5·10⁻⁷ | 4.5·10⁻⁵ | 5.1·10⁻⁴ | 2.9·10⁻⁶ | 3.3·10⁻⁴ | 2.9·10⁻³ |
| 8 | 1.6·10⁻⁸ | 4.3·10⁻⁶ | 1.1·10⁻⁴ | 7.9·10⁻⁸ | 4.3·10⁻⁵ | 8.8·10⁻⁴ |
| 10 | 8.1·10⁻¹⁰ | 7.1·10⁻⁷ | 3.6·10⁻⁵ | 2.2·10⁻⁹ | 5.1·10⁻⁶ | 2.3·10⁻⁴ |
| 12 | 3.8·10⁻¹¹ | 9.7·10⁻⁸ | 1.1·10⁻⁵ | 5.3·10⁻¹³ | 5.8·10⁻⁷ | 7.0·10⁻⁵ |
| 16 | 1.2·10⁻¹⁴ | 1.4·10⁻⁹ | 1.3·10⁻⁶ | 3.2·10⁻¹⁵ | 3.0·10⁻¹⁰ | 1.2·10⁻⁶ |

- Nearly the same for 100, 1000 and 20 000 bodies (p 8, θ 0.5: 6·10⁻⁶ … 4·10⁻⁵); the floor is
  the rounding of the direct sum, 2·10⁻¹⁵. A Barnes–Hut tree (10⁻³ … 10⁻⁴) is p = 2 … 4 here.
- Momentum: |Σ GM_i a_i| / Σ GM_i |a_i| ≤ 10⁻¹⁶ over the table (every cell interaction is
  mutual); along an integration the momentum changes by 2·10⁻¹⁷ of Σ GM |v|.
- Energy injected by the expansions along six dynamical times of a Plummer sphere of 512
  bodies, θ 0.5: about 3·10⁻⁵ (p 4), 4·10⁻⁷ (p 8), 4·10⁻⁸ (p 12) of |E|.
- Order 0 gives back the exact bits; the same input gives the same bits; receivers with an
  excluded source keep their exact sums to the bit; test bodies, the 1PN terms and a
  `myb_system` of 2000 bodies on top of p 8, θ 0.5 stay within the test bound 4·10⁻⁵.
- Over threads (`test_fmm_threads`, from 2048 bodies): the same bits as without a pool for any
  number of threads.

## Cost

- Exact sums, 2000 bodies, one thread, per interacting pair: symmetric path 15.8 ns, vector
  per-receiver path 12.2 ns on x86-64-v2 and 6.2 ns on v3 (although each pair is computed
  twice), scalar per-receiver path 17 ns; 733 bodies on 8 threads: 0.60–0.90 ms per evaluation
  against 3.8–4 ms for the symmetric path (plan, phase 4).
- FMM, one thread, p 8, θ 0.5 (an i5-9300H shared with other work, ±30 %): 10 and 16 µs per
  body for 10⁴ bodies (cube, Plummer) against 55–65 µs for the vectorised direct sum; 17 and
  25 µs for 10⁵ against about 600 µs. From 10⁴ to 10⁵ bodies the time per body grows 1.6–1.7
  times (the direct sum's 10 times); below a few thousand bodies the gain is 2 to 3 times.
- FMM over threads (`test_fmm_threads`, p 8, θ 0.5, the same i5-9300H of four cores shared with
  six other jobs, three runs), µs per body:

  | | serial | 2 threads | 4 threads | direct sum, 4 threads |
  |---|---|---|---|---|
  | cube, 10⁴ | 10.6–15.3 | 7.5–9.9 | 5.4–6.7 | 20–28 |
  | cube, 10⁵ | 21–24 | 12–17 | 9.4–11 | 200–280 (scaled as N²) |
  | Plummer, 10⁴ | 18–22 | 16–18 | 9.5–14 | 21–24 |
  | Plummer, 10⁵ | 30–46 | 21–35 | 15–21 | 210–240 (scaled) |

  0.9–1.8 times faster on two threads, 1.3–2.5 on four (about 2 at the median), for 5000 bodies
  still 1.5 and 2; the FMM on four threads is 1.7–4 times cheaper than the direct sum on four
  threads at 10⁴ bodies, 10–25 times at 10⁵.

## Limitations

- Point masses only: figures, tides and the other levels are those of
  [the system of bodies](system-of-bodies.en.md); there is no softening, and coincident
  interacting bodies are refused.
- Excluding only one direction of a pair breaks the conservation of momentum — the caller's
  choice. The invariants count every pair whatever the graph says.
- With the EIH terms the accelerations depend on the velocities: only
  `myb_nbody_acceleration_xv` computes them, and the Wisdom–Holman integrators refuse the system.
- A system is used by one thread at a time; a pool runs one job at a time, and objects sharing a
  pool must not be used from several threads concurrently (a violation terminates the process).
- The FMM is approximate: the force is not the gradient of one potential and jumps where the
  walk changes, so energy is not conserved exactly. Below 2048 bodies its evaluation is serial;
  from 2048 on it runs on the pool, each pair across a cut of the near field computed by both
  sides — on an idle machine that bounds the gain of four threads near 2.6 for 10⁴ bodies and
  3.2 for 10⁵. It allocates when the tree outgrows its memory, and with threads for the lists of
  the walk and the scratch of each thread (once, kept). With the 1PN terms the potentials U_i and
  the 1PN terms stay direct O(N²) sums; only the Newtonian a_j of the sources come from the FMM.
  The Wisdom–Holman integrators refuse a system with an approximation.
- No Barnes–Hut tree (its 10⁻³ … 10⁻⁴ force error kept it out of the library).

## Data

None.
