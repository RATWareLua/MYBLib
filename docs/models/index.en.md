# Model passports

**English** | [Русский](index.ru.md)

A passport is the page to read before trusting a model of `myb`: what it computes, where it comes
from, where it is valid, how accurate it is and how that was measured. Every passport has the
same nine parts — *What it is*, *Source*, *Interface*, *Inputs and outputs*, *Domain of
validity*, *Accuracy*, *Cost*, *Limitations*, *Data* — and every number in it is traceable to a
comment of a public header, to a test (named, with the reference it is compared with) or to the
plan ([`docs/PLAN.en.md`](../PLAN.en.md)). Where nothing was measured, the passport says so.

The column *Validation* gives one headline check; the passport has the rest.

## Time, frames and files

| Model | Passport | Header | Validation |
|---|---|---|---|
| Time scales and their conversions, TDB−TT | [time-scales](time-scales.en.md) | [time.h](../../include/myb/time.h) | against ERFA ≤ 3.2·10⁻¹¹ s; TDB−TT of DE440t — its rate to 1.4·10⁻¹⁶ (`test_erfa_time`, `test_de440_time`) |
| Earth orientation: IAU 2006/2000A, CIO, EOP, nutation table, equinox-based frames | [earth-orientation](earth-orientation.en.md) | [frames.h](../../include/myb/frames.h) | X, Y against ERFA 1.7·10⁻¹⁷ rad; sites on the Earth against Horizons 0.19 m — a model difference (0.077 m through the ITRF93 PCK) (`test_erfa_frames`, `test_earth_frames`) |
| SPK ephemerides (reader of 15 types, writer), SP3 | [spk-ephemerides](spk-ephemerides.en.md) | [kin.h](../../include/myb/kin.h) | `testpo.440` 7.1·10⁻¹⁵ au; CSPICE positions ≤ 6·10⁻¹⁶ (`test_testpo`, `test_spk`) |
| Binary PCK, CK and SCLK (read and written) | [pck-ck-sclk](pck-ck-sclk.en.md) | [frames.h](../../include/myb/frames.h), [ck.h](../../include/myb/ck.h) | CK against CSPICE: the same answers, the matrix 7·10⁻¹⁶ (types 1–3), 4.4·10⁻¹⁶ (types 4–6) (`test_ck`) |
| Rotations, WGCCRE, TK frames, the frame graph, geodesy, orbital frames | [frames-geodesy](frames-geodesy.en.md) | [frames.h](../../include/myb/frames.h) | frame graph against CSPICE `sxform` 2·10⁻¹⁶; geodetic angles 3.3·10⁻¹⁶ rad (`test_frame_graph`, `test_erfa_frames`) |
| Numerical foundations and text formats (correctly rounded functions, double-double, elliptic functions, linear algebra, numbers, CSV, CCSDS) | [numerics-formats](numerics-formats.en.md) | [cmath.h](../../include/myb/cmath.h), [num.h](../../include/myb/num.h), [io.h](../../include/myb/io.h), [ccsds.h](../../include/myb/ccsds.h) | printing of numbers against CPython — no difference in 6124 values; elliptic functions against mpmath ≤ 5.6 ulp (`test_format`, `test_elliptic`) |
| Events and geometry: eclipses and the edges of the shadow (withdrawn or predicted), ray–shape intersection | [events-geometry](events-geometry.en.md) | [events.h](../../include/myb/events.h), [geometry.h](../../include/myb/geometry.h) | eclipse contacts against CSPICE `gfoclt` 4.4·10⁻¹⁰ s; ray and ellipsoid against `surfpt` 5·10⁻¹⁶ (`test_events`, `test_geometry`) |

## Two bodies and integrators

| Model | Passport | Header | Validation |
|---|---|---|---|
| Kepler's equation, elements, universal propagator, two-body STM, mean elements, KS and Sundman | [kepler-regularization](kepler-regularization.en.md) | [kepler.h](../../include/myb/kepler.h), [regular.h](../../include/myb/regular.h) | propagator against MPFR (400 bits) ≤ 4.7·10⁻¹⁵; KS at e = 0.999 — 1.2·10⁻¹² (`test_mpfr_kepler`, `test_regular`) |
| Lambert, B-plane, shooting, relative motion | [transfers-relative-motion](transfers-relative-motion.en.md) | [lambert.h](../../include/myb/lambert.h), [transfer.h](../../include/myb/transfer.h), [relative.h](../../include/myb/relative.h) | Lambert by shooting 4·10⁻¹² m/s; Clohessy–Wiltshire 3·10⁻¹⁴ (`test_transfer`, `test_relative`) |
| DOP853, RKN12(10) with its interpolant, Gauss–Jackson | [classical-integrators](classical-integrators.en.md) | [integ.h](../../include/myb/integ.h) | RKN12(10): 10 revolutions at e = 0.9 — 4.8·10⁻¹⁵ of the period (`test_rkn12`) |
| IAS15 (with the timescale of the jerk for contact) | [ias15](ias15.en.md) | [integ.h](../../include/myb/integ.h) | 10⁴ revolutions at e = 0.3 — 2.8·10⁻¹¹ of the period; Brouwer's law (`test_ias15`) |
| Symplectic compositions (incl. 7 and 15 stages), WH/WHFast with correctors | [symplectic-integrators](symplectic-integrators.en.md) | [integ.h](../../include/myb/integ.h), [dyn.h](../../include/myb/dyn.h) | coefficients bit for bit with MPFR; observed orders 2 … 8 (`test_mpfr_compositions`, `test_whfast`) |

## N bodies

| Model | Passport | Header | Validation |
|---|---|---|---|
| Direct N-body sums, the interaction graph, threads; the FMM as an approximate level, over threads from 2048 bodies | [nbody](nbody.en.md) | [dyn.h](../../include/myb/dyn.h) | Sum2 ≤ 3.5·10⁻¹⁶ of Σ\|terms\|, bits independent of threads; FMM to the rounding of the direct sum at p = 16, its bits independent of threads too (`test_nbody`, `test_fmm`, `test_fmm_threads`) |
| The system of bodies: levels, integrated rotation, tides by band, PN pairs, derivatives by variational equations | [system-of-bodies](system-of-bodies.en.md) | [system.h](../../include/myb/system.h) | re-integration of DE440 over 5 years: the Moon 4.9 cm, Mercury 2.1 cm; derivatives of two bodies against the closed form 2.3·10⁻¹⁶ (`test_de440_system`, `test_system`) |
| Encounters, collisions, mergers, neighbour search | [system-encounters](system-encounters.en.md) | [system.h](../../include/myb/system.h) | contact of two spheres 10⁻¹² of the closed form; 300 spheres — the billiard to 5.6·10⁻¹³ s (`test_system`) |

## Gravity and tides

| Model | Passport | Header | Validation |
|---|---|---|---|
| Spherical harmonics up to degree 2700, gradient | [spherical-harmonics](spherical-harmonics.en.md) | [grav.h](../../include/myb/grav.h) | against MPFR 256 bits: 7.5·10⁻¹⁶ (degree 120), 4.8·10⁻¹⁶ (degree 2190); gradient against mpmath 2.1·10⁻¹⁵ (`test_mpfr_grav`, `test_gradient`) |
| Polyhedron, mascons, ellipsoidal harmonics | [irregular-bodies-gravity](irregular-bodies-gravity.en.md) | [grav.h](../../include/myb/grav.h) | polyhedron against the prism of Nagy et al. 7.6·10⁻¹⁶; mascons 1.2·10⁻¹⁵ (`test_grav`, `test_ellipsoidal`) |
| Figure–figure, rotational dynamics, the torque-free top, lunar librations | [figure-figure-rotation](figure-figure-rotation.en.md) | [grav.h](../../include/myb/grav.h), [attitude.h](../../include/myb/attitude.h) | free top 2·10⁻¹⁴ rad; librations against DE440 0.040″ over a year (`test_attitude`, `test_moon_libration`) |
| IERS 2010 solid tides, FES2004 ocean tides, pole tides, secondary waves | [tides-iers](tides-iers.en.md) | [grav.h](../../include/myb/grav.h) | step 1 against the addition theorem 2·10⁻¹⁶; LAGEOS against ILRS 2.05 → 1.57 cm (`test_tides`, `test_lageos`) |
| Lagged tides: constant time lag, tides by band | [tides-lagged](tides-lagged.en.md) | [force.h](../../include/myb/force.h), [system.h](../../include/myb/system.h) | lunar recession 3.79 cm a year from Mignard's formula, ten revolutions integrated to 6 digits of it; the Moon of DE440 19 cm after 5 years (`test_nongrav`, `test_de440_system`) |

## Relativity

| Model | Passport | Header | Validation |
|---|---|---|---|
| 1PN: Schwarzschild, Lense–Thirring, de Sitter, EIH with β, γ, local frames, proper time | [relativity-1pn](relativity-1pn.en.md) | [rel.h](../../include/myb/rel.h), [dyn.h](../../include/myb/dyn.h) | Mercury 42.9804″ per century (5·10⁻⁹ of the formula); EIH against an independent Lagrangian 2.7·10⁻¹⁰ (`test_rel`) |
| PN binary: 1PN–3.5PN, 4PN with the tail (on any conic), spin terms to the next order, radiation reaction | [pn-binary](pn-binary.en.md) | [rel.h](../../include/myb/rel.h) | coefficients against Blanchet 4·10⁻¹⁶; PSR B1913+16 Ṗ to 3·10⁻⁵ (`test_pn_higher`, `test_pn_binary`) |

## Non-gravitational forces

| Model | Passport | Header | Validation |
|---|---|---|---|
| Radiation pressure, penumbra, shadows of several bodies, refraction, Poynting–Robertson, albedo and IR | [radiation-pressure](radiation-pressure.en.md) | [force.h](../../include/myb/force.h) | shadow fraction against closed forms and quadrature 1.1·10⁻¹⁴; Poynting–Robertson — 7 digits of the averaged theory (`test_force`, `test_nongrav`) |
| Light on surfaces: attitude laws, meshes with self-shadowing, thermal recoil and antenna | [radiation-surfaces](radiation-surfaces.en.md) | [surface.h](../../include/myb/surface.h) | lit area 9.000000000000 m²; 63 W on 241 kg — 8.72·10⁻¹⁰ m/s² (`test_surface`, `test_nongrav`) |
| Yarkovsky and YORP | [yarkovsky-yorp](yarkovsky-yorp.en.md) | [force.h](../../include/myb/force.h), [surface.h](../../include/myb/surface.h) | response against mpmath 5·10⁻¹³; seasonal drift 4·10⁻¹³; YORP changes sign for the mirror body (`test_nongrav`, `test_surface`) |
| Comets and outgassing | [comets-outgassing](comets-outgassing.en.md) | [force.h](../../include/myb/force.h) | A2 drift 5·10⁻⁷ from the formula; a comet 5·10⁻⁶ from Gauss's equation averaged (`test_nongrav`) |
| Drag, gas–surface interaction, continuum, entry heating and ablation | [aerodynamics](aerodynamics.en.md) | [force.h](../../include/myb/force.h), [surface.h](../../include/myb/surface.h) | Sentman's plates over a sphere — the closed form to 5·10⁻¹⁵; ablation against the closed solution 1.9·10⁻¹⁴ (`test_surface`, `test_entry`) |
| Magnetic fields (IGRF), Lorentz force, charging | [magnetic-field-charging](magnetic-field-charging.en.md) | [force.h](../../include/myb/force.h) | IGRF-14 against ppigrf 1.4·10⁻¹⁵; OML against mpmath 6·10⁻¹⁶ (`test_charge`, `test_charging`) |
| Contact and landing: smooth friction or friction that sticks, under DOP853 or IAS15 | [contact-landing](contact-landing.en.md) | [contact.h](../../include/myb/contact.h) | restitution against the closed form 3.5·10⁻¹³; a lander with the friction that sticks stands on 15° to 9.3·10⁻¹⁴ m/s (`test_contact`) |
| Micrometeoroid impacts | [micrometeoroids](micrometeoroids.en.md) | [meteoroid.h](../../include/myb/meteoroid.h) | number of impacts 0.4σ from λT; masses — Kolmogorov–Smirnov 1.17 (`test_meteoroid`) |

## The environment of the Earth

| Model | Passport | Header | Validation |
|---|---|---|---|
| NRLMSISE-00, the input of upper-atmosphere models, space weather; user copies of DTM2020, NRLMSIS 2.1, HWM14 (not in the repository) | [thermosphere-models](thermosphere-models.en.md) | [myb_nrlmsise00.h](../../ext/nrlmsise00/myb_nrlmsise00.h), [weather.h](../../include/myb/weather.h) | 17 reference cases to the 7 printed digits (`test_nrlmsise00`) |
| JB2008 from the papers | [jb2008](jb2008.en.md) | [myb_jb2008.h](../../ext/jb2008/myb_jb2008.h) | against the authors' program: T∞ 1.1·10⁻⁶ K, log₁₀ ρ 8.9·10⁻⁴ rms (`test_jb2008`) |
| Troposphere: IERS ch. 9 models, VMF3/GPT3, gradients | [troposphere](troposphere.en.md) | [myb_troposphere.h](../../ext/troposphere/myb_troposphere.h) | GPT3 6·10⁻¹³ and VMF3 1.3·10⁻¹⁵ of the TU Wien routines (`test_troposphere`) |
| Ionosphere (IONEX 2-D and 3-D, higher orders) and plasma delay | [ionosphere-plasma](ionosphere-plasma.en.md) | [myb_ionosphere.h](../../ext/ionosphere/myb_ionosphere.h), [observe.h](../../include/myb/observe.h) | ray integrals against direct sums 1.3·10⁻⁹; plasma column against mpmath 6·10⁻¹⁶ (`test_ionosphere`, `test_plasma`) |
| SGP4, TLE and its fitting | [sgp4-tle](sgp4-tle.en.md) | [myb_sgp4.h](../../ext/sgp4/myb_sgp4.h) | `tcppver.out` 5.0·10⁻⁹ km (half the last printed digit); SPK type 10 against CSPICE 2.7·10⁻¹² km (`test_sgp4`) |

## Ground station and measurements

| Model | Passport | Header | Validation |
|---|---|---|---|
| Site displacements (IERS 2010 ch. 7) | [site-displacement](site-displacement.en.md) | [site.h](../../include/myb/site.h) | solid tide 1.1·10⁻⁷ m from DEHANTTIDEINEL; ocean loading 5·10⁻⁷ m from HARDISP (`test_site`) |
| Observables and light time | [observables-light-time](observables-light-time.en.md) | [observe.h](../../include/myb/observe.h), [estim.h](../../include/myb/estim.h) | light time against CSPICE 4·10⁻¹² s; aberration bit for bit with ERFA (`test_observe`, `test_light_measure`) |

## Applications

| Model | Passport | Header | Validation |
|---|---|---|---|
| Initial orbit determination: Gibbs, Herrick–Gibbs, Gauss, Laplace, double-r, Gooding | [initial-orbit-determination](initial-orbit-determination.en.md) | [iod.h](../../include/myb/iod.h) | Vallado's examples 7-3 and 7-4 to six digits; Gauss 5·10⁻¹⁵ … 3·10⁻¹⁶ (`test_iod`) |
| Orbit estimation: least squares, filter, smoother | [estimation](estimation.en.md) | [estim.h](../../include/myb/estim.h) | LAGEOS against ILRS 1.57 cm over a week; Roadster 306 m over 7 years, OSIRIS-REx 12.6 m (`test_lageos`, `test_spacecraft`) |
| Conjunctions and the probability of collision (short and long encounters) | [conjunctions](conjunctions.en.md) | [conjunction.h](../../include/myb/conjunction.h) | round Gaussian 2·10⁻¹⁶; long encounter against Monte Carlo 1.3σ (`test_conjunction`) |
| CR3BP: Lagrange points, periodic orbits, families, manifolds | [cr3bp-periodic](cr3bp-periodic.en.md) | [cr3bp.h](../../include/myb/cr3bp.h) | JPL catalogue (56 orbits): Jacobi 3.6·10⁻¹⁵, fixed point 9·10⁻¹² (`test_cr3bp`) |
| CR3BP tori: quasi-periodic orbits, their stability, normal rotation and manifolds | [cr3bp-tori](cr3bp-tori.en.md) | [cr3bp.h](../../include/myb/cr3bp.h) | a circle returns onto itself over 20 turns to 3.7·10⁻¹¹ (`test_cr3bp`, `test_torus_stability`) |
| Low thrust: thrust arcs, Gauss equations, Edelbaum, Q-law | [low-thrust-qlaw](low-thrust-qlaw.en.md) | [thrust.h](../../include/myb/thrust.h) | Q-law 1234.3 m/s against Edelbaum's 1232.6 (`test_thrust`) |
| Optimal control: indirect method, NLP, Sims–Flanagan, chains with flybys (powered too) | [optimal-control](optimal-control.en.md) | [optimal.h](../../include/myb/optimal.h) | 486.84 m/s against Hohmann's 486.82; Sims–Flanagan +0.11 % of the indirect method (`test_optimal`, `test_chain`) |

## Propagation and uncertainty

| Model | Passport | Header | Validation |
|---|---|---|---|
| Environment, force model, ready-made terms, propagation, variational equations, covariance, effect budget | [propagation-uncertainty](propagation-uncertainty.en.md) | [env.h](../../include/myb/env.h), [model.h](../../include/myb/model.h), [orbit.h](../../include/myb/orbit.h), [uncert.h](../../include/myb/uncert.h), [budget.h](../../include/myb/budget.h) | ready-made terms against the assembly by hand — bit for bit; STM against the closed form 5.5·10⁻¹⁵ (`test_env`, `test_uncert`) |
