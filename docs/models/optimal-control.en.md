# Optimal low-thrust transfers

**English** | [Русский](optimal-control.ru.md) · [All passports](index.en.md)

## What it is

Transfers with low thrust that spend the least propellant or take the least time. The indirect
method solves the conditions of Pontryagin's principle by shooting on the costates, going from
the least "energy" to the least propellant by a homotopy; it also gives the least time and the
insertion into an orbit with the phase free. The direct method of Sims and Flanagan replaces the
thrust of each segment by an impulse in its middle and solves for the throttles with the
library's solver of nonlinear programs (an augmented Lagrangian with a damped BFGS matrix and
bounds); its legs, in the field of the point mass or of a force model, can be joined into chains
by flybys, unpowered or powered (an impulse at the periapsis). What comes out is a local
optimum: the one the first guess leads to.

## Source

- Pontryagin's principle with the energy-to-fuel homotopy: R. Bertrand, R. Epenoy, *Optim.
  Control Appl. Meth.* 23, 171 (2002); the costates and the multiplier of the cost made a unit
  vector, so that a random first guess has a fair chance: F. Jiang, H. Baoyin, J. Li, *J. Guid.
  Control Dyn.* 35, 245 (2012).
- The solver: the method of multipliers (M. R. Hestenes, M. J. D. Powell, 1969;
  R. T. Rockafellar, 1973, for the inequalities) in the scheme of E. G. Birgin, J. M. Martínez,
  *Practical Augmented Lagrangian Methods* (2014); Powell's damped BFGS (1978); the projection
  onto the bounds of D. P. Bertsekas, *SIAM J. Control Optim.* 20, 221 (1982); the approximate
  Wolfe conditions of W. W. Hager, H. Zhang, *SIAM J. Optim.* 16, 170 (2005), where the fall of
  the function is below the resolution of its values.
- J. A. Sims, S. N. Flanagan, AAS 99-338 (1999), for the legs and the flybys. A powered flyby:
  the impulse along the motion at the common periapsis of the two hyperbolas (where the Oberth
  effect makes it least), each turning the excess velocity by its half.
- Implemented from the publications; no third-party code, no data files.

## Interface

Header [`myb/optimal.h`](../../include/myb/optimal.h):

- the indirect method: `myb_indirect_solve`, `myb_indirect_sample`; structures
  `myb_indirect_problem`, `myb_indirect_solution`, `myb_indirect_point`; enums
  `myb_indirect_goal` (`MYB_INDIRECT_FUEL`, `MYB_INDIRECT_TIME`), `myb_indirect_arrival`
  (`MYB_INDIRECT_RENDEZVOUS`, `MYB_INDIRECT_ORBIT`);
- the solver: `myb_nlp_solve`; `myb_nlp_problem`, `myb_nlp_options`, `myb_nlp_result`,
  callback `myb_nlp_function`; `MYB_NLP_MAX_UNKNOWNS`, `MYB_NLP_MAX_CONSTRAINTS` (256 each);
- Sims–Flanagan: `myb_sims_flanagan_mismatch`, `myb_sims_flanagan_solve`; `myb_sims_flanagan`,
  `myb_sims_flanagan_result`; `MYB_SIMS_FLANAGAN_MAX_SEGMENTS` (64);
- chains: `myb_chain_size`, `myb_chain_guess`, `myb_chain_evaluate`, `myb_chain_solve`;
  `myb_chain` (with `impulse_exhaust`), `myb_chain_body`, `myb_chain_result` (with `impulse`);
  `MYB_CHAIN_MAX_LEGS` (8); the powered flyby alone: `myb_chain_flyby`.

## Inputs and outputs

- The dynamics x″ = g(t, x, x′) + (T u/m) d, m′ = −T u/c, 0 ≤ u ≤ 1, |d| = 1; SI units (m, m/s,
  kg, s, N). The field g is the point mass μ or a force model holding everything but the engine,
  the centre included (the costates and the transition matrices use its Jacobian).
- Indirect: the unknowns are l_0, l_x, l_v, l_m at t0 as one unit vector (l_0 > 0). The thrust
  points against l_v; for the least propellant u = 0 for S > ε, 1 for S < −ε, (ε − S)/(2ε)
  between, S = 1 − (c|l_v|/m + l_m)/l_0. ε goes from 1 to `smoothing` (10⁻⁴ by default), each
  solution the guess of the next, predicted along the secant in log ε. The conditions: the state
  at the arrival and l_m = 0 there; for the least time also H at the arrival for a target moving
  on its conic. `MYB_INDIRECT_ORBIT`: the arrival anywhere on the conic through `target` about μ,
  the phase one more unknown and the transversality condition l_x·v_T + l_v·a_T = 0 in its place
  (for the least time H = 0). Inside, lengths are in the first distance from the centre, times
  in √(L³/μ), masses in the first mass; tolerance 10⁻⁹; up to 100 first guesses from a fixed
  sequence (the same result everywhere); DOP853 and a Jacobian of the conditions by central
  differences. Out: the costates, t1, the mass at the arrival, ε reached, the residual, the
  guesses and stages used, the phase. `myb_indirect_sample` gives the state and mass, thrust
  direction, throttle, switching function, Hamiltonian and costates (−l_v is the primer vector).
- The solver: minimise f(x) with c_i = 0 (equalities first), c_i ≤ 0, lower ≤ x ≤ upper; the
  caller gives f, c and, on request, their gradients. Feasibility 10⁻⁸, optimality 10⁻⁶ (the
  largest component of x − P(x − ∇L)), at most 50 updates of the multipliers and 1000 steps
  between them; the unknowns and the functions should be of the order of one. A step is halved
  along its projection into the box until the function falls enough — or, where the fall is
  below 10⁻¹² of the values (near a solution with a large weight), until the slope along the
  step has risen enough and not overshot.
- Sims–Flanagan: 2 … 64 equal segments, a throttle vector |τ| ≤ 1 for each (3 numbers),
  Δv_i = τ_i T Δt / m_i with the mass in the middle of the segment; the first half flown forward,
  the second backward, the mismatch (6: m, m/s) at the meeting point with its Jacobian. The
  solver's unknowns are the size of each throttle in [0, 1] and two angles. With a field each
  half segment is integrated by DOP853 (tolerance 10⁻¹³ in the units of the solver) with its
  transition matrix.
- Chains: up to 8 legs from body to body (providers of the state relative to the centre, at
  epoch + t), the gm and the least periapsis of every body flown by. The unknowns: the departure
  time, the time of flight of each leg, the excess velocities, the throttles; the constraints:
  the mismatch of each leg, the equal excess speeds of each flyby, and the inequalities of the
  turn (cos δ ≥ cos δmax, sin(δmax/2) = 1/e, e = 1 + r_p |w⁻||w⁺|/gm), of the launcher and of
  the arrival. The propellant (kg) is minimised; the derivatives are analytic but for the
  acceleration of a body, taken by central differences of its provider's velocity. Two passes:
  throttles as vectors with a size of their own, then a size and two angles; the weight 10⁴.
- A powered flyby (`impulse_exhaust` of its body > 0, m/s: the exhaust velocity of the engine
  of the impulse; only the bodies flown by): dv = |√(|w⁺|² + 2gm/r_p) − √(|w⁻|² + 2gm/r_p)|,
  the turn δ = asin(1/e⁻) + asin(1/e⁺), e = 1 + r_p|w|²/gm, the mass m → m exp(−dv/c). Two more
  unknowns, r_p (no smaller than the least distance) and dv ≥ 0; the turn is an equality in
  place of the equal speeds, and ±(√(|w⁺|² + 2gm/r_p) − √(|w⁻|² + 2gm/r_p)) − dv ≤ 0, so that
  the solver meets no corner at dv = 0 and drives an impulse that does not pay onto its bound.
  The first pass holds the impulses at 0, the second frees them. `myb_chain_flyby` gives dv and
  δ with their derivatives for given excess speeds and r_p.

## Domain of validity

- A local optimum only. The indirect target orbit must be an ellipse about μ, whatever the
  field; in too short a time there is no solution (the least time to the 8000 km circle is
  3633.451 s; 3 % less has none).
- The solver: up to 256 unknowns and 256 constraints, dense matrices — no replacement for large
  sparse solvers.
- Refused with `MYB_E_DOMAIN`: invalid numbers, a target orbit that is no ellipse; sizes, bounds
  that cross, negative options, a weight above 10¹², a function that is not finite; a throttle
  longer than 1 + 10⁻⁹, more propellant than mass; times out of the order of the legs, a flyby
  without excess speed, more than 256 unknowns or constraints in the first pass of a chain.
- `MYB_E_NOCONV`: no first guess led anywhere (the transfer may be impossible in the time
  given), or the way to the ε asked broke off — the solution then holds the last one found with
  its ε; the solver did not converge — x and the result hold what was reached.

## Accuracy

`test_optimal` ([tests/unit/test_optimal.c](../../tests/unit/test_optimal.c)); the transfer is
7000 → 8000 km with 0.06 rad of plane in 4500 s at 0.3 m/s²:

| What | Reference | Result |
|---|---|---|
| Indirect, least energy and least propellant | the conditions of the principle | the Hamiltonian constant to 10⁻¹²…10⁻¹⁶, the miss 10⁻³ m; least energy 762.689 kg left, least propellant 780.076 kg (ε = 10⁻⁴, 14 steps), the engine on 49 % of the time, off its limits in 0.2 % of the samples |
| Indirect at 1 m/s² over Hohmann's half ellipse | Hohmann's two impulses, 486.82 m/s | 486.84 m/s |
| Least time to a target on the 8000 km circle | the least propellant in that time | 3633.451 s, the switching function negative throughout; the least propellant in that time spends the same 363.345 kg |
| Insertion into the circle with the phase free | rendezvous at fixed places: 780.076 kg at 4.4 rad, 744.5–777.8 kg within 0.3 rad | 783.920 kg at 4.3546 rad; a rendezvous there gives the same mass to 2·10⁻⁷ kg; transversality 3·10⁻¹² of its terms; the Hamiltonian to 2·10⁻¹¹ |
| Least time into the circle | the moving target, 3633.451 s | 3011.759 s, H = 0 to 2·10⁻¹⁵ all the way |
| Hock and Schittkowski, problem 71 | f = 17.0140173 | x to 7·10⁻⁹, multipliers 0.16147 and 0.55229 |
| Rosenbrock's valley behind a wall; a linear program | known solutions | exactly on the wall; in the corner |
| Sims–Flanagan, Jacobian of the mismatch (conics) | central differences | 9·10⁻⁸ |
| Sims–Flanagan, 30 and 60 segments | the indirect method, 219.924 kg | 220.176 kg (+0.11 %; 13 segments at full thrust, 13 off) and 219.961 kg (+0.02 %) |
| Sims–Flanagan with the point mass as a force model | the conics | mismatch 8·10⁻¹⁶ of the orbit; the same least propellant (9·10⁻⁷ kg apart) |
| Sims–Flanagan with the J2 of the Earth | — | 219.892 kg; the halves meet to 4·10⁻³ m and 7·10⁻⁶ m/s; the Jacobian against differences 1.5·10⁻⁷ |

`test_chain` ([tests/unit/test_chain.c](../../tests/unit/test_chain.c)): the Earth, Venus and
Mars on circles (Mars 1.85° off the plane), 1 N on 1000 kg at 30 km/s, 2.5 km/s from the
launcher, a rendezvous with Mars, 10 segments a leg:

| What | Reference | Result |
|---|---|---|
| Derivatives of the program | central differences | 1.4·10⁻⁸ of their sizes (1.5·10⁻⁸ in a field that changes with time), the propellant 4.7·10⁻¹¹ |
| The Sun as a force model | the conics | 3.5·10⁻¹⁴ |
| Earth–Venus–Mars | the same chain kept 10⁶ km from Venus: 731.350 kg | 834.470 kg left: the flyby saves 103 kg; the legs meet to 30 m and 6·10⁻⁶ m/s; excess speeds at Venus 3·10⁻⁶ m/s apart; the turn 1.2112 of the 1.2863 rad allowed (7201 km); the launcher's 2500 m/s used to 7·10⁻⁸ m/s |
| The chain kept 10⁶ km from Venus (a turn of 0.018 rad at most) | from the guess and from the chain by Venus | 731.350 kg (9386 m/s) from both, to 1.3·10⁻⁷ kg, each in 10 updates (before the slope decided, the criticality stopped at 2·10⁻⁶ … 2·10⁻⁴) |
| A powered flyby | its closed form at 50 digits | the impulse and the turn 2.2·10⁻¹⁶, their derivatives too (and against differences 3·10⁻¹⁰), also at excess speeds 3·10⁻⁶ m/s apart |
| The chain with a powered flyby | central differences | derivatives by r_p and dv 1.4·10⁻⁸ |
| By Venus, 1 N, an impulse by an engine of 3000 m/s | the unpowered chain | dv driven to 0; the same chain to 4·10⁻⁷ kg |
| 0.25 N | the unpowered chain does not meet Mars (the flights at their longest, the constraints 0.025 au off) | powered: 380.797 m/s at 26 530 km from Venus, 587.501 kg left, from the guess and from where the unpowered chain stopped alike (6·10⁻⁷ kg); at 0.2 N 469.385 kg with 1146.879 m/s |

## Cost

Counts measured in the tests: the indirect homotopy to ε = 10⁻⁴ in 14 steps; problem 71 in 7
updates of the multipliers and 108 evaluations; the chain by Venus in 10 + 5 updates and 14 000
evaluations. No timings.

## Limitations

- Local optima only. The direct Earth–Mars transfer in the windows of the test is hard for the
  solver: coplanar it gave 712.8 kg at the middle of the windows and 763.7 kg with the dates free;
  with Mars off the plane it ended in a local optimum that spends four fifths of the mass.
- A coplanar transfer between close circles is nearly degenerate (what is spent hardly depends on
  where the engine works); the secant prediction was added after the homotopy stalled there.
- The Jacobian of the indirect conditions is taken by central differences.
- Sims–Flanagan replaces the thrust by impulses: about 0.1 % of the propellant at 30 segments.
- The solver is dense, with at most 256 unknowns and constraints.
- The flybys of the chains are instantaneous, in the field of the body alone; a powered one has
  a single impulse along the motion at the periapsis. Where an impulse pays little the solver
  finds local optima: one of 30 000 m/s pays at 0.3 N by 0.03 … 0.05 kg (two local optima from
  two starts), one of 300 000 m/s at 1 N by 0.100 kg.

## Data

None.
