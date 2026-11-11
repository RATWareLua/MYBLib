# myb — library plan

**English** | [Русский](PLAN.ru.md)

A library for computing the motion of objects in space: N bodies, orbits around planets, interplanetary
transfers, free motion. The goal is models based on the real laws of physics, without simplifications
where the effect is measurable, with results reproducible to the bit and with an error estimate.

## 1. Decisions

| Question | Decision |
|---|---|
| Language | C23 internally (Clang ≥ 19); public headers — the common subset of C17 and C++17 |
| Compiler | **Clang only, CMake only**; Windows: target `x86_64-w64-windows-gnu` (MSYS2 MINGW64 sysroot, lld), Linux: system Clang |
| Verification platforms | Windows 10 (host), Gentoo (VM) |
| License | BSD-3-Clause; only code under a compatible license is ported (§9) |
| Prefix | `myb_` / `MYB_` |
| N-body | exact direct summation only; approximate trees are not part of the core |
| Accuracy | model error below rounding noise wherever possible; where `double` is not enough — double-double; result = value + error estimate |
| Reproducibility | bitwise identical result on Windows and Linux, at any `-O`, thread count and SIMD kernel |

## 2. Principles

Large errors come not from the formulas but from their “bookkeeping”:

| Pitfall | Cost |
|---|---|
| time as a single `double` JD | step ≈ 40 µs → ~0.3 m in LEO, ~1.2 m for the Earth |
| TT instead of TDB (±1.66 ms) | ~13 m along-track in LEO, ~50 m in the Earth's position |
| UTC instead of UT1 (up to 0.9 s) | up to ~400 m at the surface |
| without polar motion | ~10 m |
| “J2000” instead of ICRF (23 mas) | ~17 km at 1 au |
| `G*M` instead of `GM` (2·10⁻⁵ vs 2·10⁻⁹) | of order km/d at GEO |
| `-ffast-math`, FMA contraction on one platform but not on another | different results; in chaotic N-body — macroscopically |
| TLE mean elements fed into Kepler's equation | kilometers |

Cross-cutting rules:
1. **GM, not G and M.** Mass is allowed only as a user parameter (then GM = G·M with σ from G).
2. **Every value carries a time scale, a reference frame and a center.** Time scales are distinct C types; the frame is a checked tag.
3. **Time has two parts:** whole seconds `int64` + fraction `double`. UTC is only a calendar form.
4. **SI** internally and in the API.
5. **GR in the post-Newtonian approximation** wherever the effect is measurable; Newton is a special case.
6. **Model passport** (`docs/models/`): formula → source with equation number → domain of validity → neglected effects and their order → reference test with tolerance. The equation number also goes into the code comment.
7. **No global mutable state**; read-only data can be shared between threads.

## 3. The core knows no specific body

The core is laws and models in general form. The Sun is a body with a luminosity; the Earth is a body with a field, rotation, an atmosphere and tides. Everything specific comes from data.

| Where | What |
|---|---|
| `libmyb` | laws, models in general form, file formats, integrators |
| `data/` | IAU/IERS conventions, body catalogs — text with the source and σ of every number; large files are downloaded |
| `ext/` | empirical models of specific bodies that cannot be reduced to data (NRLMSIS 2.1, JB2008, DTM2020, HWM14) — separate targets |
| `tests/`, `examples/` | assemble the “Solar System” from `data/` |

Body model (everything optional, `NULL` = no effect, every number has a σ and a source): GM; kind of motion (dynamical / kinematic with a provider); shape (sphere, triaxial ellipsoid, polyhedron); field (spherical and ellipsoidal harmonics with time variation, polyhedron, mascons); rotation; inertia tensor, core/mantle; tidal response (Love numbers k_nm(ω), lag/Q); atmosphere; radiation (luminosity L(t), surface albedo and IR); magnetic field; spacecraft surface (plates, optics, gas–surface, thermal).

General forms of “terrestrial” concepts:
- **Time:** coordinate time of the global reference system, coordinate time of the local reference system of any body (a generalization of TCG), proper time of a clock on any world line. TT and TDB are affine conventions from data; UTC is a leap-second table; UT1 is a correction table.
- **Body rotation:** constant; polynomials + series (WGCCRE form); CIP/CIO — X, Y, s as Poisson series + rotation angle + tabulated corrections (this is how IAU 2006/2000A + EOP are expressed); tabulated (Chebyshev from binary PCK); dynamical (Euler equations with torques).
- **Reference frames:** global, local for a body (a generalization of GCRS after Soffel et al. 2003), body-fixed, orbital (RTN, VNC, LVLH), user-defined.
- **Ephemerides:** SPK reader — a provider of “segment: ID A relative to ID B”; the meaning of the IDs is defined by the user.

## 4. Everything that affects motion

Inclusion criterion: an effect is modeled if its contribution over the computation interval exceeds rounding noise. The **effect budget** (`budget`) estimates the contribution of each term in a specific scenario and warns if a disabled effect exceeds the tolerance.

- **Gravity of bodies:** point masses; spherical harmonics; ellipsoidal harmonics (the spherical series diverges inside the Brillouin sphere); homogeneous polyhedron; mascons. *Werner & Scheeres 1997; Garmier & Barriot 2001; Holmes & Featherstone 2002.*
- **Field variability:** solid-body tides (k_nm(ω)), ocean/atmospheric tides as sets of harmonics, pole tide, secular and tabulated series. *IERS TN36 ch. 6; Wahr 1981.*
- **Extended bodies:** “figure–figure” mutual potential, torques, coupling of translational and rotational motion, fluid core. *Borderies 1978; Hou, Scheeres, Xin 2017; Williams et al. 2001.*
- **Tidal dissipation.** *Mignard 1979–80; Efroimsky & Lainey 2007; Williams & Boggs 2016.*
- **GR 1PN:** EIH with PPN β, γ; Lense–Thirring from body spins; geodetic precession of local reference systems; 1PN multipoles (DSX theory). *Will 2018; Soffel et al. 2003; Damour, Soffel, Xu 1991–94.*
- **Higher PN** (2PN, 2.5PN, spin terms) — for compact objects, later. *Blanchet 2014; Poisson & Will 2014.*
- **Source radiation:** penumbra with angular size and limb darkening, refraction in the atmosphere of the occulting body, Doppler/aberration (Poynting–Robertson). *Vokrouhlický, Farinella, Mignard 1993–94; Robertson 2015; Burns, Lamy, Soter 1979.*
- **Object surface:** plates (specular/diffuse/absorption), self-shadowing. *Milani, Nobili, Farinella 1987; Ziebart 2004; Montenbruck et al. 2015.*
- **Radiation from bodies:** albedo and IR. *Knocke et al. 1988; Rodriguez-Solano et al. 2012.*
- **Self-generated heat:** spacecraft thermal recoil, RTG; Yarkovsky and YORP; antenna recoil (P/c). *Turyshev et al. 2012; Vokrouhlický 1998–99; Bottke et al. 2006; Farnocchia et al. 2013.*
- **Plasma and fields:** solar wind pressure; Lorentz force. *Horányi 1996.*
- **Atmosphere:** drag and lift, atmospheric rotation, winds; gas–surface interaction; free-molecular and transitional regimes. *Sentman 1961; Moe & Moe 2005; Doornbos 2012.*
- **Self-generated forces:** thrust, variable mass, center-of-mass offset; outgassing; comets (A1, A2, A3, g(r)). *Marsden, Sekanina, Yeomans 1973.*
- **Non-inertial frames:** indirect terms, inertial forces of rotating frames.
- **Events:** collision with a body's shape, merger conserving momentum and mass.

Outside the first version (now phase 22, the user's decision of 2026-11-11): contact and landing, continuum hypersonic aerodynamics (user tables only), ablation, stochastic micrometeoroid impacts. Later: observables (light time, Shapiro delay, aberration, light deflection).

**Result uncertainty** (`uncert`): variational equations with respect to state and parameters, covariance propagation through the STM, Monte Carlo with a deterministic RNG.

## 5. N bodies

- Exact direct summation only: O(N²) over massive bodies + O(N·M) for test particles; target — up to ~10⁴ massive bodies.
- Interaction graph: the level is set for each pair (point, harmonics, figure–figure, 1PN, tides); the default is point + 1PN.
- Threads — through the platform layer, partitioning over i, summation over j in a fixed order: the result does not depend on the thread count and matches the single-threaded symmetric version.
- SIMD — vectorization over i, not over j; no FMA in the kernels; bit for bit with the scalar code.
- Barnes–Hut gives a force error of ~10⁻³–10⁻⁴ — not included. FMM is a separate module, explicitly marked as approximate (phase 22).

## 6. Architecture

```
L0  core, platform (files, threads, CPUID)
L1  num, cmath (own functions), io (DAF/SPK/PCK, .gfc, IERS series, CCSDS, CSV)
L2  time (general), body (body model, parameters ± σ)
L3  frames, kin (kinematics providers), kepler, integ
L4  grav (fields), rel (1PN, local reference systems, proper time), force (non-gravitational)
L5  dyn (bodies, test particles, rotation, interaction graph), env (environment and ready-made terms of a force model), uncert, budget
L6  orbit, transfer, cr3bp
ext/   empirical models of specific bodies
data/  conventions and catalogs — data only
```
A module uses only the layers below it; `integ` knows nothing about physics.

Key methods by module (sources — §8):
- **num:** double-double [Dekker 1971; Joldes, Muller, Popescu 2017; Muller & Rideau 2022], Sum2/Dot2 [Ogita, Rump, Oishi 2005], Clenshaw, root finding with a budget, Legendre, Stumpff, small linear algebra.
- **time:** IAU resolutions 2000 B1.9, 2006 B3; TDB−TT from the Fairhead & Bretagnon 1990 series or from the ephemeris; leap-second table with an expiry date (`MYB_E_RANGE` beyond it).
- **frames:** IAU 2006/2000A CIO-based scheme [Capitaine et al. 2003; IERS TN36 ch. 5], EOP with interpolation, WGCCRE [Archinal et al. 2018], Moon — binary PCK (PA frame), geodetic coordinates [Vermeille 2011], the matrix and its time derivative; nutation — grid + interpolation with a proven error bound.
- **kin/ephem:** own reader for NAIF DAF/SPK (types 2, 3, 13, 21) and binary PCK; DE440/DE441 [Park et al. 2021]; cross-check against `testpo`.
- **kepler:** elements (classical, equinoctial, modified equinoctial), Kepler's equation [Markley 1995; Raposo-Pulido & Peláez 2017; Conway 1986], universal variables and STM [Shepperd 1985]; SGP4 — only for the sake of TLE [Hoots & Roehrich 1980; Vallado et al. 2006].
- **integ:** IAS15 [Rein & Spiegel 2015], RKN12(10) [Dormand, El-Mikkawy, Prince 1987], DOP853 [Hairer et al.], Gauss–Jackson [Berry & Healy 2004], WH/WHFast [Wisdom & Holman 1991; Rein & Tamayo 2015], Yoshida, SABA/SBAB; dense output, events, compensated update, Brouwer's law.
- **transfer:** Lambert [Izzo 2015; Gooding 1990], patched conics, B-plane [Sergeyevsky et al. 1983], STM-based correction.
- **cr3bp:** Lagrange points, Jacobi, halo [Richardson 1980; Howell 1984], manifolds [Koon, Lo, Marsden, Ross 2011].
- **env, orbit:** the bodies of a scenario with providers of state and orientation, their values computed once per instant; the scales of the providers from a table of differences; ready-made terms of a force model; propagation by IAS15 with its own steps and dense output — the state at an instant does not depend on what else is asked.

## 7. Floating point and reproducibility

Contract (verified by the build and by tests):
- binary64, `FLT_EVAL_METHOD == 0`, `-ffp-contract=off`, no fast-math, `-fno-math-errno`, baseline x86-64-v2 (v3 allowed); the flags are appended last and are not overridden by user CFLAGS; `-ffast-math` → compilation error.
- No libm functions with inexact rounding: `sin`, `cos`, etc. are correctly rounded (CORE-MATH), `fma` is our own integer implementation (the mingw-w64 `fma()` is wrong on roughly 1 input in 6). The `test_no_libm_imports` test checks the library's symbol table.
- `long double` and `float` are forbidden in computations.
- Numbers from files — our own correctly rounding parser, independent of the locale (Eisel–Lemire; Lemire 2021).
- The result hash (`test_repro`) is identical on all builds; a difference = at least one bit has changed.

## 8. Sources

Books: Montenbruck & Gill “Satellite Orbits” (2000); Vallado “Fundamentals of Astrodynamics and Applications”, 5th ed. (2022); Battin (1999); Danby (1988); Murray & Dermott (1999); Szebehely (1967); Roy (2005); Explanatory Supplement, 3rd ed. (2013); Kopeikin, Efroimsky, Kaplan (2011); Soffel & Langhans (2013); Moyer, DESCANSO Monograph 2 (2003); Kaula (1966); Heiskanen & Moritz; Hairer, Nørsett, Wanner (1993); Hairer, Lubich, Wanner (2006); Aarseth (2003); Higham (2002); Muller et al. “Handbook of Floating-Point Arithmetic” (2018); Schaub & Junkins (2018); Markley & Crassidis (2014); Koon, Lo, Marsden, Ross (2011).

Standards: IERS Conventions 2010 (TN36) with updates; IAU resolutions 2000 B1.3–B1.9, 2006 B1–B3, 2009 B2, 2012 B2, 2015 B3; WGCCRE report; NAIF Required Reading; CCSDS 502.0-B, 500.0-G; ECSS-E-ST-10-04C.

Data: JPL (DE440/441, testpo, NAIF generic kernels); IERS (EOP 20 C04, finals2000A, Leap_Second.dat, series tables); ICGEM, PDS Geosciences; CelesTrak, GFZ Kp/ap; JPL Horizons, IGS/CDDIS, ILRS.

## 9. Licenses

| Source | License | Action |
|---|---|---|
| ERFA | BSD-3 | port with notice |
| CORE-MATH | MIT | port with notice (done: sin, cos) |
| fast_float | MIT / Apache / BSL | algorithm — under MIT |
| Ryu | Apache-2.0 / BSL-1.0 | under BSL-1.0 |
| NRLMSISE-00 (C port) | public domain | in `ext/` |
| fdlibm (FreeBSD msun) | SunPro (permissive, retain the notice) | ported: erf, erfc (`src/num/erf.c`) |
| IERS Conventions software | custom, not BSD | do not port; equations — from the papers, tables — as data, examples — for cross-checking |
| NRLMSIS 2.1, JB2008, DTM2020, HWM14 | various | check before porting |
| Orekit, GMAT | Apache-2.0 | do not copy; algorithms — from the papers, the projects — for cross-checking |
| REBOUND | GPL-3 | cross-checking only |
| SOFA, CSPICE | custom | ERFA instead of SOFA; own SPK reader |
| DE440, EGM2008, IERS | free | not in the repository; downloaded |

## 10. Validation

1. Invariants: invertibility of transformations, orthogonality, norms.
2. Integrators: order of convergence, time reversibility, Brouwer's law.
3. Analytical cases: the two-body problem on all conics; “figure-eight”; Lagrange; the Pythagorean problem; torque-free top; Jacobi integral; Mercury's perihelion ≈ 43″/century.
4. References: ERFA, `testpo.440`, Horizons, Vallado, GMAT/Orekit, IGS SP3, ILRS.
5. Re-integration of DE440.

## 11. Development rules

- Expected errors — `myb_status`; contract violation — `MYB_CHECK` (one line to stderr and immediate termination).
- `goto` — only a jump to a cleanup/error label at the end of the function (`fail:`, `out:`): C has no RAII, and a single label releases all resources on any error path. Any other use of `goto` is forbidden.
- Recursion is forbidden (explicit stack); every loop has a visible bound, iterative algorithms have a budget and `MYB_E_NOCONV`.
- System calls — only in `src/platform/`; no memory allocation in an integrator step.
- Code and tests go together; new code compiles without warnings under `-Werror`.

## 12. Roadmap

| Phase | Content | Completion criterion | Status |
|---|---|---|---|
| 0 | CMake skeleton (Clang only), `core`, `num`, `cmath` (fma, sin, cos), reproducibility tests | everything below is green on Windows and Gentoo | Windows, Gentoo: done |
| 1 | `time` in general form, `body`; IAU conventions as data; own number parser, NAIF text kernels | ERFA ~1 ns | Windows: done (ERFA ≤ 3.2·10⁻¹¹ s); Gentoo: done |
| 2 | `kepler`, Lambert; remaining CORE-MATH functions | Vallado, invertibility | Windows, Gentoo: done |
| 3 | `integ` | orders, Brouwer, 10⁴ revolutions | Windows: done (DOP853, RKN12(10), IAS15, Gauss–Jackson, Yoshida compositions, events, dense output, two-body transition matrix); Gentoo: done |
| 4 | `dyn` (Newton), WHFast, interaction graph | “figure-eight”, Pythagorean problem, hash independent of threads | Windows: done (Sum2 direct summation, graph, threads, vector blocks, double-double invariants, WH/WHFast with correctors 3–17, SABA1–4, SBAB1–4); Gentoo: done |
| 5 | `kin`/SPK, `frames` | testpo, Horizons | Windows: done (SPK 14 types, binary and text PCK, CIP/CIO IAU 2006/2000A, EOP with daily corrections, NAIF TK frames, geodesy, orbital frames; testpo — 7·10⁻¹⁵ au, Horizons — Eros 1.8·10⁻¹², sites on the Earth 0.08/0.19 m); Gentoo: done |
| 6 | `rel` (EIH, gravitomagnetism, local reference systems) | Mercury, DE440 | Windows: done (EIH with β, γ, 1PN integrals, Schwarzschild/Lense–Thirring/de Sitter, BCRS↔GCRS, proper time; Mercury 42.9804″/century, DE440: Mercury 106 m/yr vs 58 km for Newton); Gentoo: done |
| 7 | `grav` (polyhedron, ellipsoidal harmonics, figure–figure, tides), `force`, `ext/` atmospheres | ERFA, GMAT/Orekit, IGS/LAGEOS | Windows: done (harmonics up to degree 2700 — MPFR 7.5·10⁻¹⁶; polyhedron, mascons, ellipsoidal harmonics, figure–figure; IERS 2010 and FES2004 tides; radiation with penumbra, albedo, drag; NRLMSISE-00; LAGEOS vs ILRS — 2.0 cm over a week; DE440 with figures — Moon 0.65 m/yr); Gentoo: done (7 presets, 48/48, same hash) |
| 8 | `transfer`, `cr3bp`, `uncert`, `budget` | Horizons for real spacecraft | Windows: done (variational equations — 5.5·10⁻¹⁵ relative to the closed-form expression; JPL periodic orbit catalog, 56 orbits; B-plane and shooting; effect budget; Roadster 7 years — 306 m, OSIRIS-REx through a gravity assist — 12.6 m, JWST near L2); Gentoo: `linux-release` and `linux-asan` 53/53, same hash (the remaining sanitizers — at the end) |
| 9 | attitude, regularization, SGP4, low thrust | own references | Windows: done (SGP4 — `tcppver.out` to printed precision, SPK type 10 vs CSPICE 2.7·10⁻¹² km; equinox-based frames vs ERFA 7·10⁻¹⁷ rad; elliptic functions vs mpmath ≤ 7.4 ulp; torque-free top 2·10⁻¹⁴ rad; lunar librations vs DE440 — 0.04″ over a year; KS at e = 0.999 — 1.2·10⁻¹² over 4 thousand evaluations; low thrust vs closed-form solutions 3·10⁻¹³; stack — frames ≤ 4 KB); Gentoo: `linux-release` and `linux-asan`, same hash |
| 10 | higher PN, observables, charge and plasma | own references | Windows: done (2PN — circular-orbit invariant (27/2 − 7ν); PSR B1913+16 — Ṗ vs Peters–Mathews 3·10⁻⁵; light time vs CSPICE 4·10⁻¹² s; aberration bit for bit with ERFA; IGRF-14 vs ppigrf 1.4·10⁻¹⁵); Gentoo: all 7 presets, same hash |
| 11 | 3PN and 3.5PN, spins, plasma and troposphere, ionosphere, charging | own references | Windows: done (coefficients vs machine parsing of Blanchet 4·10⁻¹⁶; K at 3PN 0.06%; J0737−3039B 5.0735°/yr; IERS ch. 9 10⁻¹⁵; OML vs mpmath 6·10⁻¹⁶); Linux — CI |
| 12 | environment of a force model, ready-made terms, propagation (`env`, `orbit`) | LAGEOS from the ready-made terms with the same result | Windows: done (ready-made terms against the assembly by hand — 3.5·10⁻¹⁶ of the acceleration; LAGEOS — the same 2.0 cm; table of scales against the exact conversion 4·10⁻¹⁶ s; a state does not depend on the sampling — bit for bit); Gentoo: done (`linux-release` 78/78, the new tests under `linux-asan`, same hash) |
| 13 | packaging: installation, CMake package, pkg-config, shared library, `examples/`, CHANGELOG | an outside project builds through `find_package(myb)`; the library exports only `myb_*` | Windows: done (`examples/` against the installation — static and shared; `libmyb.dll` and `libmyb.so` — 442 symbols each, all of the API); Gentoo: done |
| 14 | analytic Jacobians (gradient of the harmonic field, radiation, drag, relativity); orbit estimation: least squares over an arc, a filter, partials of the observables | LAGEOS and the spacecraft of phase 8 through the library's estimator; state transition matrix with a field ≤ 10⁻¹² | Windows and Gentoo: done (gradient of the field against mpmath 2·10⁻¹⁵, at the poles 6·10⁻¹⁶; state transition matrix in a field of degree 8 symplectic to 6·10⁻¹⁷; LAGEOS by the library's estimator — 2.05 cm, the spacecraft of phase 8 — the same 306 m, 12.6 m, 260 m) |
| 15 | output: printing of numbers, writers of SPK (types 9, 13, 2, 3), CCSDS OEM/OPM/OMM (reading and writing), SP3, TLE, CSV | the files are read by CSPICE; writing and reading — bit for bit | Windows and Gentoo: done (printing against CPython — 0 differences on 6124 values; an SPK of 34 segments is read by CSPICE N0067, positions 3.9·10⁻¹⁶; the ILRS SP3 excerpt comes back character for character; 33 TLE sets — bit for bit; OEM: numbers exact, 1.5 % within the 1 ulp of the km ↔ m conversion) |
| 16 | events and geometry: eclipses, apsides, nodes, visibility, close approaches; breaking the integrator's step at the edges of a shadow; ray–shape intersection, sub-point, illumination angles; collision with a shape | CSPICE references (`gfoclt`, `sincpt`); LAGEOS independent of the integrator's precision | Windows and Gentoo: done (eclipse contacts against `gfoclt` 4·10⁻¹⁰ s; ray and ellipsoid against `surfpt` 5·10⁻¹⁶, nearest point against `nearpt` 4·10⁻¹⁶; `sincpt`, `subpnt`, `ilumin` — 10⁻⁸ m and 7·10⁻¹⁶ rad; LAGEOS: residual of the first guess 25.9587 m at ε = 10⁻¹⁰, 10⁻¹¹, 10⁻¹², fit 2.05 cm) |
| 17 | forces of §4: Poynting–Robertson, Yarkovsky and YORP, thermal recoil and antenna, comets and outgassing, lift and gas–surface interaction, winds, plates with attitude and self-shadowing, refraction in the shadow, tidal dissipation, variable mass, impulsive manoeuvres | own references from the papers | Windows and Gentoo: done (Poynting–Robertson: the fall of a and e over 20 revolutions as in the averaged theory, 7 digits; Yarkovsky A2 — 5·10⁻⁷ from the formula; a comet — 5·10⁻⁶ from Gauss's equation averaged; the lagging tide — 2·10⁻⁷, 3.79 cm a year for the Moon; Sentman's plates over a sphere — 5·10⁻¹⁵ from the closed form; self-shadowing: lit area 9.000000000000 m²; YORP changes sign for the mirror body; impulses: Φ across them 2·10⁻⁹). Refraction in the shadow, the Yarkovsky effect of a conducting body and YORP with conduction were not done in the phase itself: they were closed in "Debts of phases 12–17" |
| 18 | N bodies and the body model: levels of the graph (harmonics, figure–figure, tides — the lagged tides of the Earth and the Moon among them), coupled translation and rotation (with the derivative of the rotational deformation), close encounters, collisions and mergers, PN binaries among N bodies | re-integration of DE440 by the library's own means (the Moon with its librations) | Windows and Gentoo: done (DE440 again with the librations of the Moon, from JD 2440400.5: the Moon 1.1 cm after 30 days, 8.8 cm after a year, **19 cm after 5 years** (54 m without the tides of the Earth); the attitude of the Moon 0.003″, 0.041″ and 0.16″; Mercury 0.77 m and Venus 0.56 m after a year; the rest is the small bodies of DE440 that are missing; a free top 3·10⁻¹³; the contact of two spheres 10⁻¹² of the closed form) |
| 19 | applications: initial orbit determination, relative motion, the approach of two objects and the probability of collision, mean elements and TLE fitting, multiple shooting, pseudo-arclength continuation of families and quasi-periodic orbits, variational equations with DOP853 in `myb_orbit`, optimal low thrust | Vallado, the JPL catalogue | Windows and Gentoo: done (orbit determination — the examples 7-3 and 7-4 of Vallado to six digits, Gauss — the orbit to 10⁻¹⁴; relative motion; closest approaches and the probability of collision, the covariances carried by the orbits themselves — 2.5·10⁻⁹ of the closed transition matrix; mean elements and TLE fitting — elements to 10⁻¹⁵; multiple shooting and families — a member closes to 4·10⁻¹³; **quasi-periodic orbits** — the circle of a torus returns onto itself over 20 turns to 3.7·10⁻¹¹; the variational equations with DOP853; low thrust: the Q-law — 1234.3 m/s against Edelbaum's 1232.6, the **indirect method** — the Hamiltonian constant to 10⁻¹², at high thrust 486.84 m/s against Hohmann's 486.82, **Sims–Flanagan** — 0.11 % of the propellant from the indirect method with 30 segments) |
| 20 | files and frames: CK and SCLK (read and written), a graph of frames (SPK across frames), SPK type 15, nutation grid, TDB−TT from the ephemeris, space weather; bodies of a system whose motion is given (a provider) and the small bodies of DE440 from their SPK; the ground station: station displacements (IERS ch. 7), the admittance of the secondary waves of the ocean tide and the atmospheric S1/S2, VMF3/GPT3 and the APG gradients, three-dimensional ionosphere maps and the higher orders of the ionosphere, light time in the measurements; `ext/`: NRLMSISE-00 adapter, JB2008, DTM2020, NRLMSIS 2.1, HWM14 (after the licence checks) | CSPICE; references of the models | — |
| 21 | engineering: fuzzing of the parsers, benchmarks, coverage, model passports `docs/models/`, ARM64 and macOS in CI, consumers built with MSVC; bindings: the upper layers for the rest of the interface, array forms of the calls in C, Go and Haskell, tests under the sanitizers; the edge of a shadow without withdrawing the step (a speed-up, judged by measurement); vcpkg, Conan, PyPI and crates.io packages (they need a tagged release) | CI green on every target; the orbit of the example from every binding bit for bit with the C calls | — |
| 22 | the tails of the models and what §4 and §5 kept outside the first version: optimal low thrust beyond the rendezvous — insertion into an orbit with a free phase and the field of a model in Sims–Flanagan, chains of legs with flybys; the stability and the manifolds of tori; 4PN with the hereditary tail, higher-order spin-orbit terms and the radiation reaction of spins; secondary emission and the limit a ≫ λ_D in charging; compositions with smaller error constants (7 and 15 stages) and a dense output of RKN12(10) for `myb_orbit`; orbit determination by the double-r iteration and the methods of Laplace and Gooding; the probability of a long encounter; atmospheric entry — continuum hypersonic aerodynamics and ablation; contact and landing; random impacts of micrometeoroids; FMM as an approximate level of the graph and a neighbour search for the encounters of many bodies; in the system of bodies — mergers of rotating bodies, derivatives with respect to the state and the parameters, the Lense–Thirring effect and the loss of mass of the Sun; thermal inertia in the torque of light at the moment | own references from the papers | — |

### Completion after 1.1 (phases 12–22)

The review of 2026-11-10 found the physics of the core covered and the layer above it missing.
Phases 12–21 go from what the rest cannot be built without (ready-made terms and propagation,
packaging) to estimation, output, events, the missing forces of §4, the full N-body problem,
applications, files and engineering. The C++, Python and Rust bindings follow phase 13 (they
need the installation and the shared library) and rest on phase 12: a force model is assembled
from data, not from callbacks.

What phases 12–17 left under "Not done" was closed in a separate pass before phase 18 (the
section "Debts of phases 12–17" below); what moved to later phases is named there with its
reason.

The rule for what any phase leaves (2026-11-11): what can be closed without undue trouble is
closed at once; what is harder now gets the phase where it will cost least or give most (or
the final phase 22); an item is dropped only by the user's decision. The leftovers of phases
3–11 that belonged to no phase were sorted by this rule too: the section "Old debts (phases
3–11) and the bindings" below.

### Phase 19 decisions

The phase is closed: everything is done and checked on Windows and Gentoo.

- **Initial orbit determination** (`myb/iod.h`). Gibbs — three positions, exact for any
  conic: the velocity to 1.4·10⁻¹⁵ for e from 0 to 1.4; Vallado's example 7-3 gives (0,
  5.531148, −5.191806) km/s. Herrick–Gibbs — close positions, the error falling as the fourth
  power of the spacing (ratios 16.0 and 16.0); example 7-4 to six digits. **Gauss** — three
  directions from three places: the roots of the polynomial of the eighth degree give the
  candidates, and each is refined by Newton's method on the six conditions (the object lies
  on each sight) with the transition matrix of the two-body motion. A low orbit over 4 and 12
  minutes, an eccentric one over half an hour, a navigation orbit over 2 hours, geostationary
  height over 5 hours: the orbit comes back to **5·10⁻¹⁵…3·10⁻¹⁶** in 4–5 steps (the series
  of Gauss alone gives 3·10⁻⁵…5·10⁻³). The refinement by simple iteration of Curtis's scheme
  diverges in these geometries (a multiplier near −3) — hence Newton.
- **Relative motion** (`myb/relative.h`): the state of a deputy in the R, T, N frame of the
  chief, exactly, both ways; the linear map of the relative state for **any conic** — from
  the transition matrix of the two-body motion taken into the turning frame at both ends
  (what Yamanaka–Ankersen gives on an ellipse and Clohessy–Wiltshire on a circle, without a
  formula per kind of orbit). The error of the map is of the second order of the separation:
  at 10 m, 100 m and 1 km the ratios are 100 exactly on a circle, an ellipse of e = 0.7 and a
  hyperbola; on a circle it agrees with the closed form of Clohessy and Wiltshire to
  3·10⁻¹⁴. A two-impulse rendezvous (`myb_relative_rendezvous`): from 5 km in 40 minutes,
  flown exactly, it misses by 11.6 m — what the linearisation gives. The relative elements of
  D'Amico (`myb_roe_from_states`, `myb_roe_to_position`): constant to 2·10⁻¹⁰ over a
  revolution, but for the drift of the second.
- **Closest approaches and the probability of collision** (`myb/conjunction.h`). The closest
  approaches of two objects are the minima of their distance: for providers
  (`myb_approach_find`, to the resolution of the time) and for two orbits being propagated
  (`myb_approach_orbits`, a cubic over the ends of a step; against the providers 2·10⁻¹¹ s
  and 10⁻⁸ m with steps of 10 s). Two circles crossing at a node and coplanar circles are the
  closed forms. The probability of a short encounter (`myb_encounter_plane`,
  `myb_collision_probability`): the integral of the Gaussian over the disk is reduced to one
  integral of a difference of error functions (through erfc in the tails, with no
  subtraction of close numbers), Gauss–Legendre panels doubled to 10⁻¹³. A centred round
  Gaussian against 1 − exp(−R²/2σ²): 2·10⁻¹⁶; an offset one (to 30σ, probabilities down to
  8·10⁻¹⁹⁰) against the series of Rice: 4·10⁻¹⁴; unequal σ against a sum over a grid of
  4000 × 4000: 2·10⁻⁶. Chan's series (`myb_collision_probability_chan`) is summed with
  positive terms (as usually written it loses every digit at large misses); the largest
  probability over the size of the covariance (`myb_collision_probability_max`) against
  R²/(e d²): 10⁻⁶.
- **Mean elements and TLE fitting**. Averaging (`myb_mean_elements_average`, `myb/kepler.h`):
  under J₂ the osculating semi-major axis swings by 11.5 km and the mean one keeps to
  0.75 m, e to 4·10⁻⁷; the drift of the node agrees with the first-order theory to 7·10⁻⁴,
  the rate of the mean longitude to 4·10⁻⁷. An element set of SGP4 fitted to positions
  (`myb_tle_fit`, `myb_tle_from_state` in `ext/sgp4`): Gauss–Newton in unknowns without the
  singularities of the classical angles, the derivatives by differences of the theory
  itself. From positions of the theory the set comes back: a low orbit over two days with
  B* — the elements to 10⁻¹⁵, B* to six digits, residuals of 5·10⁻⁸ m (8 steps from the
  osculating elements); an eccentric one over three days and a 12-hour orbit (the deep-space
  branch) over four likewise.
- **Multiple shooting and continuation of families** (`myb_cr3bp_orbit_from_state`,
  `myb_cr3bp_orbit_correct`, `myb_cr3bp_orbit_continue` in `myb/cr3bp.h`). The orbit is cut
  into parts equal in time; the unknowns are the states at the cuts and the period; the
  conditions are that the parts join (one component of the closing is given by the Jacobi
  integral) and a phase. No symmetry is assumed; the steps are of minimum norm
  (`myb_shoot`). Continuation is by pseudo-arclength: the tangent of the family is the null
  vector of the Jacobian of the conditions. The L1 Lyapunov orbit of the JPL catalogue is a
  fixed point; pushed off by 2·10⁻⁵ at every cut it returns to a closure of 2·10⁻¹⁰; 25
  steps along the family: every member closes to **4·10⁻¹³** and agrees with the corrector
  that assumes the symmetry to 7·10⁻¹⁷ in the state and 3·10⁻¹⁴ in the period.
- **The variational equations with DOP853 in `myb_orbit`** (left by the debts of 12–17): the
  variations are part of one system of the first order under the same tolerances; the
  transition matrix over 5 revolutions is 3·10⁻¹¹ of the closed form; an impulse leaves it
  as it is.
- **Low thrust: the Q-law** (`myb_qlaw`, `myb_qlaw_direction`, `myb_qlaw_arrived`,
  `myb_term_qlaw` in `myb/thrust.h`; Petropoulos 2004). Q measures the distance to the target
  orbit (the distance in each element over the largest rate at which it can be changed
  anywhere on the orbit); the thrust points where Q falls fastest — the gradient of Q carried
  through the Gauss equations. Neither an optimisation nor a first guess is needed. Raising
  7000 → 10 000 km: **1234.3 m/s against Edelbaum's 1232.6** (0.14 %) in 14.00 days; with
  the plane turned by 0.1 rad 1746 m/s against the optimum of 1641 (6 % above); from an
  ellipse of e = 0.725 to a circle at 42 164 km in 26.6 days, the periapsis never below
  6738 km with the limit at 6600. Found on the way: an engine that always burns **dithers
  near the target and stops advancing** — 100 km from it, half a million integrator steps a
  day: the gradient made a unit vector jumps where it passes through zero. The remedy is
  Petropoulos' effectivity threshold (`myb_qlaw_term.effectivity`): the engine waits where
  thrust at that place of the orbit gains less than a given part of the best. A term has no
  memory of when it burned, so with a threshold it overestimates the mass spent. **The mass
  in the state** (`myb_qlaw_flight_rhs`, `myb_qlaw_fly`): a seventh component that falls
  only while the engine burns. The same transfer to 42 164 km: 26.5 days, the engine on for
  22.6 of them, **2020 m/s against the 2389** a term without memory would count; without
  pauses the mass agrees with the rate of the engine to the step of the look. The field
  without the engine is the point mass or any force model.
- **Approaches with the covariances carried by the orbits themselves**
  (`myb_conjunction_find`). Both orbits are created with the variational equations, the
  covariances of the states are given at the start; at every approach the covariance of the
  relative position is Φ_r P Φ_rᵀ of the first plus that of the second. An approach lies
  inside a step and an orbit does not go back, so the rows of the positions of Φ are taken on
  the same cubic as the states: the rows of the velocities of Φ are their rates. Whatever the
  force models hold acts on the uncertainty as it does on the motion. Against the closed
  transition matrix of the two-body motion (4 approaches, steps of 20 s): the covariance
  **2.5·10⁻⁹**, the probability 2.6·10⁻⁹.
- **The indirect method** (`myb_indirect_solve`, `myb_indirect_sample` in the new
  `myb/optimal.h`): the conditions of Pontryagin's principle, solved by shooting on the
  costates at the start. The least propellant in a given time by a homotopy from the least
  "energy" (ε = 1, a smooth problem) to the least propellant (ε → 0, the engine full on or
  off; Bertrand, Epenoy 2002), each solution the first guess of the next, predicted along the
  secant in the logarithm of ε. The unknowns are the seven costates and the multiplier of the
  cost, together a unit vector: a random point of the sphere is then a fair first guess
  (Jiang, Baoyin, Li 2012); the guesses come from a fixed sequence, so the result is the same
  everywhere. The least time has the engine always on, the time of the arrival unknown, and
  the transversality condition of a moving target. The field is the point mass or a force
  model (the costates go through its Jacobian). Measured: 7000 → 8000 km with the plane
  turned by 0.06 rad in 4500 s at 0.3 m/s² — the Hamiltonian constant to **10⁻¹²…10⁻¹⁶**, the
  miss at the target 10⁻³ m; the least energy leaves 762.689 kg, the least propellant
  **780.076 kg** (ε = 10⁻⁴, 14 steps), the engine on 49 % of the time and off its limits in
  0.2 % of the samples only. At 1 m/s² over Hohmann's half ellipse **486.84 m/s against the
  486.82** of the two impulses. The least time to a target on the circle of 8000 km is
  3633.451 s, the switching function negative throughout; the least propellant in that very
  time spends the same 363.345 kg, and in 3 % less time there is no solution. Found on the
  way: (1) a coplanar transfer between close circles is nearly degenerate — what is spent
  hardly depends on where the engine works (486.8 m/s for Hohmann's, 487.4 for the spiral),
  and the homotopy stopped at ε = 4·10⁻⁴ on it; the secant prediction was added, and a path
  that breaks off returns the last solution found with its ε; (2) the first guesses that lead
  to a solution start slowly — a residual of 0.2…0.9 after ten Newton steps — so they cannot
  be rejected early.
- **A solver of nonlinear programs** (`myb_nlp_solve`): equalities, inequalities and bounds
  on the unknowns, up to 256 of them, dense matrices. The method of multipliers (the
  augmented Lagrangian; the scheme of Birgin and Martínez), inside it a quasi-Newton step
  with Powell's damped BFGS matrix on the free variables and the projection onto the bounds
  (Bertsekas 1982). Problem 71 of Hock and Schittkowski: f = 17.0140173, x to 7·10⁻⁹, the
  multipliers 0.16147 and 0.55229; 7 updates of the multipliers, 108 evaluations.
  Rosenbrock's valley behind a wall ends exactly on the wall; a linear program in its corner.
- **Sims–Flanagan** (`myb_sims_flanagan_mismatch`, `myb_sims_flanagan_solve`): the thrust of
  a segment replaced by an impulse in its middle, conics between the impulses; the first half
  flown from the start, the second backward from the arrival, and the two must meet. The
  mass is taken in the middle of a segment: what is spent is then linear in the throttles and
  the condition on the mass at the meeting point disappears — six are left. The unknowns of
  the solver are the size of a throttle between 0 and 1 and two angles of its direction: the
  engine being off is a bound and not a corner of the function. The Jacobian of the mismatch
  against differences: 9·10⁻⁸. The transfer of the indirect method: 30 segments — **220.176 kg
  against 219.924** (+0.11 %), 13 segments at full thrust and 13 off; 60 segments —
  219.961 kg (+0.02 %). The earlier note that the indirect methods need a solver with
  inequalities was wrong: only Sims–Flanagan does.
- **Quasi-periodic orbits** (`myb_cr3bp_torus_from_orbit`, `myb_cr3bp_torus_continue`,
  `myb_cr3bp_torus_state` in `myb/cr3bp.h`): a torus as an invariant circle of the flow over
  a time T — the curve that the flow carries onto itself turned by an angle ρ (Gómez, Mondelo
  2001; Olikara, Scheeres 2012); the turn through the trigonometric interpolant on an odd
  number of points. The equation of invariance holds two conditions fewer than it seems: the
  flow keeps the Jacobi integral and the action of a closed curve, ∮p·dq. Two multipliers at
  the gradients of these invariants take their place (the matrix of the symplectic form in
  the coordinates x, v is derived from the equations of motion); they vanish at a solution,
  and what is left of them tells whether the points are enough. Two phases hold the points
  along the circle and along the flow; one condition picks a family (the Jacobi constant, T
  or ρ kept), followed by pseudo-arclength. The first guess is the eigenvector of the
  monodromy matrix for a pair on the unit circle. Measured: about a stable distant retrograde
  orbit of the JPL catalogue a circle of 21 points returns onto itself over 20 turns to
  **3.7·10⁻¹¹** (checked with another integrator, DOP853), the multipliers 10⁻¹⁵, the Jacobi
  constant over the circle 4·10⁻¹⁵; four steps along the family to a size of 0.02 —
  1.3·10⁻¹⁰; a step along the family and back — 1.1·10⁻¹¹; both pairs of the orbit give tori,
  out of its plane and in it. Lissajous: about a small Lyapunov orbit of L1 (38 km)
  ρ = 0.1764423 against the 0.1764459 of the linear motion; about one of 4500 km (unstable
  2400-fold in a turn) the torus reaches **3617 km** out of the plane and returns after a
  turn to 1.6·10⁻¹¹. For systems of more than 64 unknowns the LU solver has an internal form
  without the limit on the order.
- **Moved on**: orbit determination by the double-r iteration and the methods of Laplace and
  Gooding — phase 22 (the same inputs as Gauss; for arcs whose sights are nearly coplanar);
  the probability of a long encounter (non-linear relative motion) — phase 22; in optimal low
  thrust — insertion into an orbit with a free phase (a rendezvous with a given state only,
  now), the field of a force model in Sims–Flanagan (the point mass, now), chains of legs
  with flybys — phase 22; the stability of the tori and their manifolds — phase 22; the upper
  layers of the bindings for the new interface — phase 21.
- **Verified** (2026-11-11): Windows — `windows-release` 100 tests of 100, the fourteen
  affected ones at -O0 (with the stack check and the ban on libm), twelve under UBSan; Gentoo
  (Clang 21.1.8) — `linux-release` 100 of 100, a build without warnings, the same hash and
  the same digits of the solvers; the thirteen affected tests under `linux-asan` — nothing
  reported.
- **Reproducibility**: the workload of the whole phase — a flight by the Q-law through a
  pause, the carried covariances, both problems of the indirect method, Sims–Flanagan,
  problem 71 and a torus with a step along its family among it; the reference is
  `0xb22481187ba3923a`, the same at -O0 and -O2; `MYB_REPRO_PREFIX=18` checks
  `0xed85437b94d0e024`. The tests of the solvers take the sine and the cosine of their input
  from the library: a solver follows its input to the last bit, and at -O0 the function of
  the system differed from the one the compiler folds.

### Phase 18 decisions

- **A system of bodies with models** (`myb/system.h`, `myb_system`). The point masses stay
  with the `myb_nbody` inside — its graph, its 1PN terms, its threads and the same bits; on
  top of them every level is asked for by name:
  - a *field* (`myb_system_add_field`): the figure of a body acts on another body as on a
    point, which acts back (momentum is kept) and turns the body. Against the same sum made
    by hand: the same bits; momentum 2·10⁻¹⁶;
  - *figures* (`myb_system_add_figures`): the quadrupoles of two bodies on each other, with
    the torques;
  - a *tide* (`myb_system_add_tide`, `myb_system_set_tide`): the tide one body raises on
    another acts on a third, or on the raiser itself;
  - a *distortion* (`myb_system_add_distortion`): a body whose rotation is integrated is
    distorted by the tide of another and by its own spin;
  - a *PN pair* (`myb_system_add_pn_pair`) and an *action* (`myb_system_add_action`) — below.
- **The tide by band, with lags** — the model of the JPL ephemerides (Folkner et al. 2014,
  eqs. 31–32): the Love numbers k₂ⱼ of the long-period, diurnal and semi-diurnal bands, each
  with its orbital lag (where the raiser was) and its rotational lag (the angle the body has
  turned by). The potential is written band by band, U ∝ Σ k₂ⱼ Wⱼ(r, r*ⱼ)/|r*ⱼ|⁵, which is
  k₂ G m R⁵ P₂/(r³ r*³) for equal k₂ⱼ; the acceleration is its gradient, the reaction goes to
  the body. An earlier position is taken from the motion of the moment: a Taylor series to
  the third order (the acceleration of the point masses and the jerk of the pair) — no
  history is kept. Against the gradient by differences 7·10⁻¹¹; the earlier position against
  Kepler motion 3·10⁻¹¹; the lagging part against Mignard's first-order formula
  (`myb_force_tidal_lag`) 5·10⁻³ at a lag angle of 4·10⁻³.
- **A rotation integrated with the motion** (`myb_system_set_rotation`): a mantle that may
  hold a fluid core — the model of the Moon of the same ephemerides. The quaternion is taken
  as a second-order coordinate, q″ = (q′ω + qω′)/2 with ω = 2 q̄ q′, and the angular velocity
  of the core as the rate of an auxiliary coordinate, so that **IAS15 integrates the motion
  and the rotation as one system**, with no first-order integrator beside it. The inertia
  tensor of the mantle is distorted by the tide (the position of the body and the frame of
  the mantle both lagged) and by the spin; the degree 2 of the field of such a body comes
  from its inertia tensor. A free top over 30 days: **3·10⁻¹³** of the exact motion,
  |q| − 1 = 5·10⁻¹⁵; a body turned by a point mass: orbit and spin exchange 4·10⁻⁵ of the
  angular momentum and keep their sum to **4·10⁻¹⁶**.
- **The rate of the rotational deformation** (left by phase 9). The rate of the spin part of
  the tensor is linear in the ω′ looked for: d/dt[f(ωωᵀ − ω²/3)]ω = f(ω²ω′ + ω(ω·ω′)/3), so the
  equation of the mantle is solved as (I + K)ω′ = …. Measured on the Moon: the attitude changes
  by 0.0003″ in 5 years, the position by 1 cm; the estimate "orders below the torques" held,
  and the term stays in the model.
- **Groups of coordinates in IAS15** (`myb_ias15_options.groups`, `myb_system_integrator`).
  The criterion of the step and the test of convergence take maxima over the blocks of the
  **whole** state; beside accelerations in m/s² the second derivatives of a quaternion
  (1/s²) are not seen at all — the orbit alone chose the step, and the top diverged
  (|q| − 1 = 10³). With groups (the positions of all bodies; each rotation; each pair of
  spins) each is judged by itself, the shortest step and the worst convergence counting.
  Without groups the bits are the former ones (the reproducibility hash did not change).
- **DE440 again** (`test_de440_system`) — the criterion of the phase. Only the state at JD
  2440400.5 is taken from the ephemeris (positions from `de440.bsp`, the angles and angular
  velocities of the Moon and of its core from `header.440`); the system does the rest: EIH;
  J₂ of the Sun on every body; the zonals J₂…J₅ of the Earth and the field of the Moon to
  degree 6 on the Earth or the Moon, the Sun, Mercury, Venus, Mars and Jupiter; the
  quadrupoles of the Earth and the Moon; the tides of the Earth from the Moon and the Sun on
  the Moon with the five lags of DE440; the rotation of the Moon with its core. Differences
  from DE440:

  | | 30 days | 1 year | 5 years | 5 years without the tides of the Earth |
  |---|---|---|---|---|
  | Moon (geocentric) | 0.011 m | 0.088 m | **0.19 m** | 54.3 m |
  | Mercury | 0.13 m | 0.77 m | 23.7 m | 23.7 m |
  | Venus | 0.61 m | 0.56 m | 62.3 m | 62.3 m |
  | Earth | 1.7 m | 13.6 m | 32.8 m | 32.6 m |
  | Mars | 2.1 m | 349 m | 4.0 km | 4.0 km |
  | attitude of the Moon | 0.003″ | 0.041″ | 0.16″ (1.4 m at the surface) | — |

  The planets are bounded by what the model lacks: the 343 asteroids and 30 trans-Neptunian
  objects of DE440 move the barycentre (the Sun by 130 m in a year) and Mars. The pole of the
  Earth is the CIP (IAU 2006/2000A); DE440 has a simpler model. Five years of both models
  take 8 s (1828 steps, 77 498 evaluations).
- **Encounters, collisions, mergers** (`myb_system_set_radius`, `myb_system_advance`,
  `myb_system_merge`, `myb_system_bounce`). The integrator steps to the first event: the
  contact of two spheres — the instant at which their distance falls to the sum of the radii,
  found on the polynomial of the step and made the end of the step (IAS15 withdraws the
  step) — or the closest approach of a pair within a watched distance. A pass with its
  minimum inside a step is not missed. Two spheres falling together: contact at
  554.689657182 s, **10⁻¹²** of the closed form, the centres 100 m − 9·10⁻¹² m apart; an
  elastic rebound returns them to 5000.000000000 m; a merger keeps the mass and the momentum
  and says what it took (the energy of the relative motion — 9 digits — and its angular
  momentum); the radius is that of both volumes. A fly-by: the closest approach against the
  pericentre of the hyperbola, 10⁻¹¹ in time and 10⁻¹³ in distance. The search goes over all
  pairs (N² a step).
- **PN pairs among N bodies** (`myb_system_add_pn_pair`): the terms beyond 1PN of two bodies
  (`myb_pn_binary_evaluate`: 2PN, 2.5PN, 3PN, 3.5PN, spin-orbit, spin-spin) are added to the
  point masses; the relative acceleration is shared by the masses. The terms themselves are
  taken, not the difference of two nearly equal accelerations
  (`myb_internal_pn_binary_terms`). Against the binary 9·10⁻¹⁷; the spins are part of the
  state.
- **Actions and light on a mesh** (`myb_system_add_action`, `myb_system_light_action`; left
  by phase 17): whatever else pushes a body or turns it, as a function. The ready-made action
  is light on the mesh of a body with the torque **of the moment**, the momentum of the light
  absorbed included, which the mean of `myb_force_yorp` leaves out. The torque passed through
  the equations of rotation against the mesh itself: 8·10⁻¹⁵. Why it is needed, measured:
  over a turn with the Sun at one latitude the re-emission gives −6.572·10⁻⁴ pR³ about the
  axis (as `myb_force_yorp`) and the light absorbed −1.026·10⁻³: it vanishes only after the
  average over the orbit.
- **Not done, and where it went**: bodies of a system whose motion is given (a provider),
  and with them the small bodies of DE440 from their SPK — phase 20 (reading files and the
  graph of frames); mergers of bodies whose rotation is integrated or that belong to a PN
  pair (`MYB_E_DOMAIN` now), derivatives of the system with respect to the state and the
  parameters, the Lense–Thirring effect and the loss of mass of the Sun in a system, a
  neighbour search in place of the pass over all pairs, thermal inertia in the torque of
  light at the moment — phase 22; the upper layers of the bindings for the system — phase 21.
- **Reproducibility**: the workload of the phase — every level in one evaluation, three days
  of motion with a rotation and spins, a chain of encounters with rebounds and a merger; the
  reference is `0xed85437b94d0e024`; `MYB_REPRO_PREFIX=17.2` checks the former
  `0x1366eda4f7ba62af`.
- **Verified** (2026-11-11, the working tree, not committed): Windows — `windows-release` 99
  tests of 99, the eleven affected ones at -O0 (with the stack check and the ban on libm);
  Gentoo (Clang 21.1.8) — `linux-release` 99 of 99, a build without warnings, the same hash,
  and the same figures of `test_de440_system` to the last digit (the Moon 0.0108 / 0.0884 /
  0.1880 m, its attitude 0.0031″ / 0.0407″ / 0.1643″); the 13 affected tests under
  `linux-asan` — nothing reported.

### Check on Gentoo (phases 14–17)

2026-11-10, Clang 21.1.8, glibc; the working tree, not committed. `linux-release`: **89 tests
of 89**, LAGEOS (71 s) and the spacecraft of phase 8 among them; the reproducibility hash
`0x99a2ea1ff9d038a9` — the same as on Windows at -O0, -O2 and under UBSan. `linux-asan`: the
15 new and affected tests (printing of numbers, tables, CCSDS, SP3, the SPK writer, SGP4,
events, geometry, small forces, the surface, propagation, estimation, the environment, the
gradient, reproducibility) — no reports. On Windows the same 15 tests under UBSan — no
reports.

### Bindings: C++, Python, Rust

Made after the debts of phases 12–17 (the directory `bindings/`). One scheme for the three
languages, two layers:

- **the lower one is the whole C interface**, with no hand work: in C++ the headers
  themselves (they are C++17 already); in Python the declarations that the preprocessor
  gives from those headers; in Rust the declarations of bindgen. No function is left out,
  and a difference from the headers is caught by a test;
- **the upper one is in the manner of the language**, over the path most programs take:
  kernels and time scales, SPK ephemerides, a gravity field, the environment, the force
  model with its ready-made terms, propagation, events, impulses, the transition matrix, OEM
  messages in both notations. It gives ownership (an object of the library is released by
  itself; a force model keeps alive the environment, the fields and the structures of its
  terms), errors in the manner of the language and callables of the language in place of
  function pointers. The rest (N bodies, estimation, CR3BP, targeting, attitude,
  observables, the SPK/PCK/SP3 writers, …) is taken through the lower layer: the two mix
  freely.

The bindings do **no arithmetic**: a number obtained through them has the bits of the same
calls made in C. That is what is tested: the orbit of the example `examples/low_orbit.c` (a
field with J2 and drag, a day, 15 ascending nodes, 776 steps), built from the classes of each
binding, gives the times of the nodes and the final state **bit for bit** with the C calls.

- **C++** (`bindings/cpp/include/myb/myb.hpp`, CMake target `myb::cpp`, C++17, header only).
  No exceptions: a function that can fail returns `myb::result<T>` — the status, the message
  of a parser and the value; asking a failed result for its value is a broken contract and
  stops the program with a message, as `MYB_CHECK` does in the library. Classes
  `kernel_pool`, `time_system`, `ephemeris`, `gravity_field`, `environment`, `force_model`,
  `orbit`, `ephemeris_message`; terms of a model and events are lambdas. The header has no
  exception handling at all — nothing is thrown and nothing is caught, and it compiles with
  `-fno-exceptions`: a lambda reports a failure by its status, and an exception that leaves
  it ends the program (it is called from C through a `noexcept` function; both such
  functions stand in `myb::detail` next to a `static_assert(noexcept(…))`, so that `noexcept`
  cannot be dropped unnoticed — the header would stop compiling). The test
  `test_cpp_binding` holds 123 checks under the warnings of the library itself; the example
  `examples/low_orbit.cpp` builds against an installation, static and shared
  (`test_install`).
- **Python** (`bindings/python`, the package `myb`, cffi in ABI mode: no compiled part, the
  package loads the shared library). The declarations `_cdef.h` (2171 lines) are written by
  `tools/make_cdef.py`: the preprocessor of Clang over the headers with empty `stddef.h`
  and `stdint.h` — only the library's own comes out, the same text on any system; a test
  regenerates and compares. `myb.lib` and `myb.ffi` are the whole interface; the classes
  `KernelPool`, `TimeSystem`, `Ephemeris`, `GravityField`, `Environment`, `ForceModel`,
  `Orbit` and functions. An error is the exception `myb.Error` with the status of the
  library. An exception raised in a Python function that the C code called (a term, a
  provider, an event) is kept, the callback reports a failure, and the exception is raised
  from the call of the library that was going on. The libraries of `ext/` come on request:
  `myb.extension("sgp4")`. The library is looked for by `MYB_LIBRARY`, in `_lib/` of the
  package and where the system looks. The test `test_python_binding` holds 14 tests.
- **Rust** (`bindings/rust`: the crates `myb-sys` and `myb`; neither depends on another
  crate). `myb-sys` holds the declarations of bindgen 0.72, kept in the crate: building it
  needs no libclang. They are generated for one fixed target, so that the text does not
  depend on the machine; what differs between platforms comes through `c_long` and `usize`,
  which Rust resolves where the crate is built. `build.rs` finds the library by
  `MYB_LIB_DIR` or through pkg-config; static by default, shared only with the MSVC
  toolchain (the static library is built for the GNU ABI). `myb` has types with `Drop`,
  `Result<T, myb::Error>` with the status, `Fn` closures for terms, providers and events; a
  closure reports a failure through `Err`, and a panic is not caught — in a function called
  from C it ends the program; shared objects go through `Rc`, the types are not `Send` (an object of the library has working memory for
  one thread). The modules of `ext/` are behind features of the same names. The test
  `test_rust_binding`: 7 integration tests, 2 tests of `myb-sys` (with the library of
  `ext/sgp4`), the example of the documentation; clippy has nothing to say. One of the seven
  runs itself again as a second process in which the closure of a term, and then of an
  event, panics: the process ends with "panic in a function that cannot unwind" and not with
  the exit code of a test that failed by unwinding.
- **The tests of the bindings** are registered with ctest when there is something to run
  them with: a C++ compiler, a Python with cffi (`MYB_PYTHON` names the interpreter), cargo.
  Under the sanitizers only the C++ one runs.
- **Checked** (2026-11-11, together with the Jacobians of the surface terms and of the
  Yarkovsky term): Windows — Clang 23 with libstdc++ 15, Python 3.13 with cffi 2.0, Rust 1.90
  (the GNU toolchain), `windows-release` 94 tests of 94, the affected tests at -O0 (with the
  stack check) and under UBSan; Gentoo — Clang 21, Python 3.14 with cffi 2.1.1 (from
  Portage), Rust 1.93, `linux-release` 94 of 94, the 11 affected tests (`test_cpp_binding`
  among them) under `linux-asan`; reproducibility hash `0xb76153a8dca45827` on both. On both systems the orbit built from the classes of each
  binding agrees with the C calls bit for bit, and `test_install` builds the C++ example
  against an installation too.

### Bindings: what is ahead

- **Haskell** (asked for by the user on 2026-11-11; by the user's decision at the end, like
  Go: phase 21). A complete binding
  on the scheme of the three that exist. The lower layer is the whole C interface: the
  `foreign import` declarations and the `Storable` instances of the structures, generated
  from the headers (hsc2hs takes the offsets of the fields from the headers themselves;
  nothing by hand). The upper one is in the manner of the language over the same path:
  objects with finalizers (`ForeignPtr`), a failure as a value (`Either` with the status of
  the library), terms of a model, providers and events as Haskell functions through
  `foreign import ccall "wrapper"`; a model keeps alive what its terms point to. No exception
  handling, as in the other bindings: a function reports a failure by its status. No
  arithmetic — the orbit of `examples/low_orbit.c` must agree with the C calls bit for bit. A
  cabal package over the static library (GHC on Windows is a mingw-w64 toolchain too; whether
  its C runtime is the one the library is built for will be found out when it is done, the
  shared library being the fallback), a test in ctest when GHC and cabal are there; on Gentoo
  they come from Portage and stay. By phase 21 phases 18–20 have stopped changing the
  interface, and the upper layer is written once.
- **Go** (the verdict of 2026-11-11, phase 21). Plain cgo, the static `libmyb.a` through
  pkg-config, the coarse layer only: one crossing into C is one task (propagate an orbit,
  find the events, fit, read or write a file). The cost of a crossing (tens of nanoseconds)
  does not depend on the work inside, so a model is assembled from data and everything "per
  epoch" and "per point" goes as one call over a slice. The array forms that are missing
  (ephemerides, rotations, time, Kepler) are added to the library itself in C — Python gains
  from them as well. Go functions as terms and events are a separate, plainly slow path.
  There is no complete lower layer: the types of cgo do not leave their package. purego does
  not serve (structures by value not on every platform, callbacks cannot return a `double`),
  nor does a port to Go (the compiler may fuse `x*y + z` into an FMA).
- **The upper layers for the rest of the interface** (phase 21): N bodies, estimation,
  CR3BP, targeting, attitude, observables, the writers of files — in all the bindings at
  once, when phases 18–20 have finished changing those parts.
- **Publication** (phase 21, after the tag): PyPI and crates.io, with vcpkg and Conan; they
  need the user's accounts.
- **What the bindings do not check yet** (phase 21): consumers built with MSVC — the C++
  header under `cl` and the MSVC toolchain of Rust through the DLL (it needs the import
  library `myb.lib`; the development machine has no MSVC linker); macOS; the tests of Python
  and Rust under the sanitizers.

### Old debts (phases 3–11) and the bindings

The review of 2026-11-11. The "Not done" lists of phases 3–11 still held items that belonged
to no phase of the roadmap; as many had gathered around the bindings. Every item is now
closed or has a phase; nothing is dropped (the user's decision of 2026-11-11).

**Closed at once**

- **The budget of effects to first order** (`myb_budget_linear`; phase 8). One propagation
  instead of one per term: along the trajectory of the model every term drives its own
  variation d″ = (∂a/∂x) d + (∂a/∂v) d′ + a_k, that is d(t) = ∫Φ(t, s) a_k(s) ds. A 400 km
  orbit over a day against the full differences of `myb_budget`: the pressure 36.02082 m
  against 36.02081 m, the drag that was switched off 27.76 km against 27.82 km (0.2 %), J2
  1108 km against 1113 km (0.4 %); the same flags. The variation has no floor of round-off:
  a term of 10⁻¹⁴ m/s² gives 4.2253·10⁻⁶ m and repeats to 4·10⁻⁷ when another term is
  changed, while the difference of two propagations gave 4.28·10⁻⁶ and 4.76·10⁻⁶ m in two
  runs. The work grows as the number of terms and not as its square — but only when the
  terms bring their Jacobians (all the ready-made ones do): a term differentiated by
  differences is computed 13 times per value of the model, and in a model of five such
  terms J2 was computed 334 863 times against 80 335 for the six propagations. The
  description of the function says so.
- **The time element of KS** (`MYB_REGULAR_KS_ELEMENT`, `myb_ks_element_from_cartesian`,
  `myb_ks_element_to_cartesian`, `myb_ks_element_rhs`; phase 9). The time is replaced by the
  element τ = t + (u·u′)/h, whose derivative is constant in Kepler motion (μ/2h — to the last
  bit in the test), so the error of integrating the time is gone. An orbit of e = 0.95, 20
  revolutions, DOP853: at the tolerance 10⁻⁶ the error is 1.93·10⁻⁵ → 5.27·10⁻⁶, at 10⁻⁹
  1.13·10⁻⁸ → 1.48·10⁻⁹, at 10⁻¹² 5.99·10⁻¹² → 1.05·10⁻¹², with no more right-hand sides
  (785 → 773, 1793 → 1781, 4217 → 4169). With J₂ (e = 0.7, three days) it differs from the
  plain form by 7·10⁻⁶ m. Bound orbits only: an unbound one is refused (`MYB_E_DOMAIN`); if
  a perturbation takes the energy through zero on the way, the run ends with an error — the
  plain KS form stays for that.
- **Any Sundman exponent** (`myb_regular_problem.alpha` ≥ 0; phase 9): 1, 3/2 and 2 by
  multiplications and a square root as before (the same bits), the others as exp(α ln r)
  with the library's own functions. e = 0.9, 5 revolutions, right-hand sides: α = 0 — 9317
  (the Cartesian equations 9374), 0.5 — 7409, 1 — 6053, 1.25 — 5249, 1.5 — 5945, 2 — 5621,
  2.5 — 7337. The description in the header was corrected on the way: steps uniform in the
  true anomaly belong to α = 2 (r² dθ ∝ dt), not to 3/2 as it said; 3/2 is the intermediate
  anomaly of Nacozy.
- **XML: CDATA sections and a document type with an internal subset** (debts of phases
  12–17). The value of a keyword or a line of data may be given as one CDATA section — taken
  as written, with no entities; numbers in sections inside a line of data are read. A
  declaration `<!DOCTYPE … [ … ]>` is passed over whole — with its comments, instructions
  and quoted values that hold `]` and `>`; the entities it declares are not known to the
  reader (`&name;` stays as written). A section mixed with other text, two sections in a row
  and a section that does not end are refused with the line number. `test_ccsds` holds 563
  checks (536 before).
- **The bindings: what had not been checked.** C++ — see above (`static_assert`). Rust — the
  test with a panic in a closure. CI gained steps that install cffi and the Rust toolchain
  of the GNU ABI and then check that the tests of the bindings are registered (they would be
  skipped silently otherwise); CI itself has not run with these steps yet — it runs on a tag
  or on request.

**Assigned to phases**

| Item | From | Phase | Why there |
|---|---|---|---|
| Lagged tides of the Earth and the Moon in the N-body problem | 7 | 18 | a level of the graph of interactions, which that phase introduces |
| The derivative of the rotational deformation | 9 | 18 | the model of the Moon's rotation moves from a test into the library; the remaining 0.04″ against DE440 is the criterion of the phase |
| PN binaries among N bodies | 10, 11 | 18 | already there |
| The torque of the light absorbed in YORP (`myb_force_yorp` leaves it out: over a rotation and an orbit it averages to zero) | debts 12–17 | 18 | coupled translation and rotation needs the torque of the moment, not the mean |
| Quasi-periodic (Lissajous) orbits, pseudo-arclength continuation | 8 | 19 | built by the multiple shooting of the same phase |
| Variational equations with DOP853 in `myb_orbit` | debts 12–17 | 19 | multiple shooting needs the transition matrix of an arc from any integrator |
| Optimal low thrust (Q-law, indirect methods, Sims–Flanagan) | 9 | 19 | already there |
| Admittance of the secondary waves of the ocean tide, atmospheric S1/S2 | 7 | 20 | the ocean loading of the station displacements uses the same admittance over the same tables of arguments; S1/S2 gives the pressure loading there |
| VMF3, GPT3, the APG gradients | 11 | 20 | the model of a ground station; GPT3 is a grid file and reading files is the subject of the phase; its gradients replace APG |
| 3D ionosphere maps and the higher orders | 11 | 20 | the second and third orders need the integrals ∫N B cos θ dl and ∫N² dl along the ray, which the three-dimensional maps give (the field is IGRF, already there); a two-dimensional map would need an empirical profile, to be rewritten later |
| DTM2020 | 7 | 20 | next to JB2008 and NRLMSIS 2.1: the same space-weather indices as input |
| CK reader, nutation grid, SPK across frames, SPK type 15, TDB−TT | 5, 7, 9 | 20 | already there |
| The derivative of a measurement of an offset point through the attitude that follows the orbit (a part offset/distance, `myb_measure_offset`) | debts 12–17 | 20 | the measurement model is rewritten there for the light time; the attitude laws are differentiated already (`src/force/surface.c`) |
| Upper layers of the bindings for the rest of the interface, Go, Haskell, publication of the packages, MSVC, macOS, sanitizers | bindings | 21 | see "Bindings: what is ahead" |
| The edge of a shadow without withdrawing the step | debts 12–17 | 21 | a speed-up, judged by measurement, and the benchmarks come with this phase; see below |
| 4PN with the tail, higher-order spin-orbit terms, the radiation reaction of spins | 11 | 22 | no phase opens this part again; the formulas are long (hundreds of terms), the reference is the coefficients of Blanchet's review |
| Secondary emission and the limit a ≫ λ_D in charging | 11 | 22 | the same: nothing else touches charging; a model of the electron yield of materials has to be chosen |
| Compositions with smaller error constants (7 and 15 stages) | 3 | 22 | the coefficients are published, but the rule of the project is to check the order conditions in MPFR and to refine them there |
| RKN12(10) in `myb_orbit` | debts 12–17 | 22 | it needs a dense output of its own order — the order conditions are solved in MPFR, as for the compositions; see below |
| Continuum hypersonic aerodynamics, ablation | §4 | 22 | the Mach number and the pressure coefficient need the temperature and composition of the atmosphere from the models of phase 20; see below |
| Contact and landing | §4 | 22 | rests on the coupled translation and rotation of phase 18 |
| Random impacts of micrometeoroids | §4 | 22 | no phase makes them cheaper: the generator and the impulses are there already |
| FMM | §5 | 22 | attaches to the levels of the graph that phase 18 introduces |

**The user's decision of 2026-11-11: nothing is dropped**

I had proposed to drop five items — for the balance of work and use, and because the plan
itself had put some of them outside the first version from the start. None of them is
impossible; each has a phase now (the rows of the table above), and here is what each needs.

- **The edge of a shadow without withdrawing the step** (phase 21). The first variant —
  ending the step at the predicted edge — lost in accuracy: 7·10⁻⁴ m at ε = 10⁻⁹ against
  8·10⁻⁷ m for withdrawing the step, because the predicted root is not the true one and the
  kink stays inside the next step. Withdrawing the step is exact but costs +26 % of the time
  of the LAGEOS fit. To be done otherwise: bring the end of the step to the edge itself —
  predict, step, refine the edge by the root on the dense output and end the step there,
  without withdrawing a long one. The gain is work, not accuracy, so measurement decides.
- **RKN12(10) in `myb_orbit`** (phase 22). Events and states between the nodes need a dense
  output, and no output of order 10–12 is published for this method: it has to be built with
  additional stages, the order conditions solved and checked in MPFR. The limit of the
  method itself stays — forces that do not depend on the velocity; a model with drag or
  relativity will be refused by `myb_orbit` with it (`MYB_E_DOMAIN`).
- **4PN and the higher spin terms** (phase 22) — without reservations.
- **Continuum hypersonic aerodynamics** (phase 22). The modified Newtonian theory (the
  pressure coefficient C_p,max sin²θ on the windward elements) on the same plates and mesh
  with self-shadowing that serve the free-molecular regime already, and the transition
  between the regimes by the Knudsen number: `myb_gas_bridge` is there, but the continuum
  regime in it is one number now, `continuum_cd`. C_p,max depends on the Mach number and the
  adiabatic index, that is on the temperature and composition of the air — the atmosphere
  models of phase 20 give them.
- **Ablation** (phase 22, with the hypersonics): the heating at the stagnation point and the
  loss of mass; the variable mass of an object is there already
  (`myb_object.mass_function`).
- **Contact and landing** (phase 22): the forces of contact with the surface of a shape — a
  normal one, elastic with damping, and friction; the rebound. The event of a collision with
  a shape is there; the rotation of the object together with its motion comes with phase 18.
- **Random impacts of micrometeoroids** (phase 22): a flux of particles as a Poisson process
  of impulses; the reproducible generator and `myb_orbit_kick` are there.
- **FMM** (phase 22): an approximate level of the graph of interactions, marked as such;
  exact direct summation stays the default rule.

**Reproducibility**: the new computations form a workload of their own (the time element,
the Sundman exponents, the linear budget, a message with CDATA sections); the reference is
`0x1366eda4f7ba62af` at -O0 and -O2. `MYB_REPRO_PREFIX=17.1` checks the previous
`0xb76153a8dca45827` and `17` checks `0x99a2ea1ff9d038a9`: no earlier bit changed, those of
the plain KS form and of the exponents 1, 3/2 and 2 included.

**Verified** (2026-11-11, the working tree, not committed): Windows — `windows-release` 94
tests of 94, the nine affected tests at -O0 (with the stack check and the ban on libm), five
under UBSan; Gentoo (Clang 21.1.8) — `linux-release` 94 of 94, a build without warnings, the
same hash `0x1366eda4f7ba62af`; the seven affected tests under `linux-asan` (regularisation,
the budget, the messages, propagation, the variational equations, reproducibility, the C++
binding) — nothing reported. CI has not run with the new steps for the bindings.

### Check on Gentoo (debts of phases 12–17)

2026-11-11, Clang 21.1.8, glibc; the working tree, not committed. `linux-release`: **91 tests
of 91**, LAGEOS and the spacecraft of phase 8 among them; reproducibility hash
`0xa762f98ab1f66694` — the same as on Windows at -O0, -O2 and under UBSan. The shared
libraries are linked with their version maps: `libmyb.so` has 607 symbols with the versions
`MYB_1.0`, `MYB_1.1`, `MYB_1.2`, `MYB_PRIVATE`; the libraries of `ext/` refer to them
(`myb_spk_segment_data@MYB_1.0`, `myb_internal_check_failed@MYB_PRIVATE`); `test_exports`
compares the version of every symbol with the map. `linux-asan`: the 16 new and affected
tests (events, geometry, SP3, the SPK and PCK writers, CCSDS, small forces, shadows, surface,
propagation, environment, estimation, SGP4, tables, the field, reproducibility) — nothing
reported. On Windows: `windows-release` 91 of 91, the same 16 at -O0 (with the stack check)
and under UBSan.

### Debts of phases 12–17

Each of phases 12–17 left a "Not done" list. It was closed before phase 18; below is what was
done, how it was checked and what moved on, with the reason.

**Estimation (phase 14)**

- **Biases of the measurements as unknowns** (`myb_fit_options.bias_count`,
  `myb_observation.bias`; up to 32): bias number k − 1 + i is added to the i-th computed
  value of an observation; the unknowns are ordered state, parameters of the model, biases.
  Three range delays of 35, −12 and 4.7 m and a range-rate offset of 2.1 mm/s are recovered
  with σ = 0.19 m, the largest error 1.1 σ, weighted RMS 0.963. The same in the filter
  (`myb_filter_set_bias`, `myb_filter_bias`).
- **A smoother** (Rauch–Tung–Striebel; `myb_filter_options.history`, `myb_filter_smooth`,
  `myb_filter_step`): the filter records its instants and a backward pass gives the estimate
  of each from the whole arc. The first of 140 instants against a least-squares fit over the
  same arc: estimates **0.000 σ**, sigmas 7·10⁻⁶; position σ 520 m before the pass, 0.94 m
  after.
- **A measured point away from the centre of mass** (`myb_measure_offset`): an antenna or a
  reflector at an offset in the body axes, with an attitude law; the position of the point
  is exact, its velocity within 3·10⁻¹⁰ m/s.
- **Jacobians of the radiation of a planet and of the empirical accelerations** are analytic
  (`myb_term_planet_radiation_jacobian`: a quadrature of the derivative of the integrand,
  `myb_force_planet_radiation_gradient`; `myb_term_empirical_jacobian`): 6·10⁻¹⁰ and
  4.5·10⁻⁴ from differences (the second number is the accuracy of the differences
  themselves over the grid of rings).
- Moved to phase 20: station displacements (IERS ch. 7) and light time in the measurements —
  they need the loading models and the graph of frames of that phase.

**Forces (phases 12 and 17)**

- **Analytic Jacobians of the new terms** — of all of them, by forward differentiation in six
  variables (`src/num/dual.h`): the formula of a term is written once more on numbers that
  carry their derivatives.
  - Comets and asteroids (`myb_term_nongrav_jacobian`, 1.4·10⁻¹⁰ and 5·10⁻¹¹ from
    differences), the lagging tide (`myb_term_tidal_lag_jacobian`, 2·10⁻⁹).
  - **Yarkovsky** (`myb_term_yarkovsky_jacobian`, `myb_force_yarkovsky_gradient`,
    `myb_force_yarkovsky_season_gradient`). The response of the body depends on the flux only
    through the thermal parameter θ, and 1/R is linear in θ, so dR/dθ = −R(1 − R)/θ — a
    closed form for a body of any size. The seasonal part is the series differentiated in the
    mean anomaly. Against differences: 1.3·10⁻¹⁰, with the season 5·10⁻¹⁰. (In the first pass
    this Jacobian was left to differences on the grounds that no closed derivative exists:
    that was wrong.)
  - **The terms that turn with the object** — light on plates and on a mesh, a gas on plates,
    the recoil of an antenna (`myb_term_surface_radiation_jacobian`,
    `myb_term_surface_drag_jacobian`, `myb_term_recoil_jacobian`): a = Mᵀ f(M w, …), and
    everything is differentiated — the force in the body axes (Sentman's formulas with exp
    and erf, the transitional regime), the direction, the flux, the density and the attitude
    M itself. The four attitude laws of the library are differentiated as they are written
    (they are known by the address of the function); a law of the caller's has no derivative
    to give, and central differences of **the law itself** are taken — 12 calls of it, not 12
    evaluations of the term. A mesh: the derivative of the force in the direction for the
    parts that are lit now (which are lit changes in steps). Against differences of the
    terms with the leading error removed: light on plates 1·10⁻¹⁰…7·10⁻¹⁰ with all five laws,
    a mesh 1·10⁻¹⁰…5·10⁻¹⁰ (with and without ray casting), the gas 2.4·10⁻⁹ (the limit is the
    gradient of the density by differences of its function), the recoil 6·10⁻¹²…3·10⁻¹¹. In a
    penumbra the Jacobian of the light is differences of the term, as for
    `myb_term_radiation`.
  - **The gradient of a wind** is now in the Jacobian of the drag (`myb_term_drag_jacobian`
    used to leave it out): a wind that changes from place to place enters by differences of
    its function, as the density does; 2.4·10⁻⁹.
  No ready-made term is left to the differences of the model.
- **Delayed outgassing** (`myb_force_nongrav_delayed`, `myb_nongrav_term.delay`): g is taken
  at the point of the Keplerian orbit `delay` seconds earlier (Yeomans & Chodas 1989).
  2·10⁶ s after the perihelion the acceleration is that of the perihelion without delay
  (2.579·10⁻⁸ m/s²).
- **The Yarkovsky effect of a conducting body** (`myb_yarkovsky.conductivity`,
  `myb_yarkovsky_response`): the answer of a sphere to a periodic heating, 1/(1 + λζ(z)),
  ζ = z j₁′(z)/j₁(z) — the solution of linear heat conduction; for a radius below 0.3 of the
  depth of the thermal wave a series to z¹² (the closed form loses eight digits there).
  Against mpmath 5·10⁻¹³ over the whole range of sizes. A grain of 2 mm: the diurnal drift
  vanishes, as it must. The diurnal force anywhere on any orbit is `myb_force_yarkovsky`,
  the term `myb_term_yarkovsky`.
- **The seasonal effect on an eccentric orbit** (`myb_yarkovsky_season`,
  `myb_force_yarkovsky_season`, the member `season` of the term): the heating along the spin
  axis is a Fourier series in the mean anomaly, and each harmonic k is answered with the
  response of the body at the frequency kn (the linear theory of Vokrouhlický and Farinella
  1999). As many harmonics as bring the terms down to 10⁻¹⁴ of the first: 2 on a circle,
  38 at e = 0.3, 111 at 0.6, 551 at 0.85; the limit is e ≤ 0.95. Against mpmath (25 digits,
  another number of quadrature nodes, the response from the Bessel function): accelerations
  **8·10⁻¹⁴**, mean drift **4·10⁻¹³** on orbits up to e = 0.85; on a circle the drift agrees
  with the closed form to 5·10⁻¹⁵.
- **YORP with conduction** (`myb_force_yorp`): the torque averaged over a rotation, the
  temperature of the faces lagging by the thermal parameter Θ. Its axial part **does not
  depend on Θ** (−2.1924·10⁻³ pR³ at Θ = 0, 1 and 30) — the result of Čapek and Vokrouhlický
  2004; the transverse parts change severalfold.
- **The torque of the tide** (`myb_torque_tidal_lag`): the torque on the body from the same
  lagging tide, in the units of `myb/attitude.h` — for the change of its rotation.
- **The indirect term of the oblateness** (`myb_term_field_indirect`): the acceleration of
  the origin (the central body) from the attraction of the perturbing bodies by its figure.
  The Earth and the Moon: −1.4834·10⁻¹¹ m/s², as J2 in closed form.
- **The shadow of several bodies at once** (`myb_shadow_fraction_many`, up to 8): every circle
  of the disk of the source loses the **union** of the arcs the bodies cover, cut at their
  limbs and at the intersections of their circles. Against a count of points of the disk
  1.7·10⁻⁶; the product of the fractions, which was taken before, is off by 7·10⁻² (and
  gives 0.15 instead of zero when two bodies cover the disk between them). In the terms of a
  model one body on the disk is the product as before (the same bits), two or more the
  union.
- **Refraction in the atmosphere of the occulting body** (`myb_refraction`,
  `myb_refraction_bending`, `myb_shadow_fraction_refracted`; the member `atmospheres` of the
  radiation terms): every element of the disk is seen at an angle moved by the bending of
  the ray, α(h) = (n − 1)√(2π(R + h)/H), with the weight of the compression of the image and
  the extinction along the path. A ray that grazes the Earth is bent by 65.7′ (two
  horizontal refractions of 34.5′); with next to no air the fraction is the airless one to
  3·10⁻¹⁰; without extinction the light only moves: the flux missing behind the Earth, in
  units of its cross-section, is 1.000559 against 1.000554 for the rays that meet the
  ground. 6000 km from the axis, inside the airless umbra, clear air leaves 2.2 % of the
  light.
- **The transitional flow regime** (`myb_gas_knudsen`, `myb_gas_bridge`; the members
  `length`, `molecule`, `continuum_cd` of the surface drag term): the free-molecular share
  sin²(π(3 + lg Kn)/8) between Kn = 10⁻³ and 10. Air at 100 km: free path 0.14 m, a body of
  2 m has Kn = 0.069, share 0.44.
- **Other integrators in `myb_orbit`** (`myb_orbit_options.integrator`): DOP853, with the same
  states at the times asked for, events and breaks of the step (two bodies over
  5 revolutions 7·10⁻⁴ m, 289 steps). RKN12(10) is left out: it has no dense
  output of its own order and needs forces that do not depend on the velocity (a dense output
  for it is phase 22 — see "Old debts"). The variational equations go with IAS15
  only (`MYB_E_DOMAIN` otherwise); with DOP853 — phase 19.
- **A copy of an environment for a thread** (`myb_env_clone`): the same bodies and providers,
  its own tables and cache of an instant.

**Events and geometry (phase 16)**

- **Predicting the edge of a shadow before the step: the first variant was tried and lost to
  withdrawing the step; it is taken further in phase 21.** The polynomial of
  an IAS15 step, carried beyond its end, predicts the edge; the step was ended short of it.
  Measured over 13 revolutions with 54 edges: the error of the final state is 7·10⁻⁴ m at
  ε = 10⁻⁹ and 5·10⁻⁶ m at 10⁻¹² — against 8·10⁻⁷ m for withdrawing the step, with the same
  number of steps. The predicted root differs from the true one and the kink stays inside
  the next step; to remove it the step has to be withdrawn anyway. Withdrawing the step
  stays.
- **A tree of bounding boxes for rays on a mesh**: the faces are halved along the long side
  of their box down to four in a leaf; it is built with the polyhedron and walked without
  recursion. A mesh of 16 128 faces, 3000 rays: **1.0 µs a ray** against 640 µs over every
  face, no difference in the distance or in the face. Self-shadowing and YORP go through it.
- **The limb and the terminator of an ellipsoid** (`myb_ellipsoid_limb`,
  `myb_ellipsoid_terminator`): the limb is an ellipse (centre and semi-axes), against CSPICE
  `edlimb` 2·10⁻⁷ m; points of the terminator of an extended source (the edge of the full
  shadow and of the full light): the horizon condition holds to 5·10⁻⁵ m, for the points of
  `edterm` to 1.4·10⁻⁴ m.
- **Geometry corrected for light time** (`myb_view_body`, `myb_view_intercept`): the body is
  taken where it was, and turned as it was, when the light left it. A ray against `sincpt`
  with `CN`: 7·10⁻⁹ m and 5·10⁻¹² s; the Sun from the object against `spkpos` with `CN`:
  9·10⁻⁵ m. With aberration the difference from `CN+S` is 13 m at 1 au: CSPICE takes the
  classical formula, the library the exact Lorentz transformation (they differ by (v/c)²,
  25 m). Without the correction a point on the Earth is 8 m off.
- Moved to phase 19: the approach of two propagated objects — with the probability of
  collision, which needs it.

**Files (phase 15)**

- **The XML notation of the CCSDS messages** (505.0-B): the readers tell it by the first
  character and parse it **with the same line parser** — the XML scanner hands over
  keywords, sections and lines of data as the lines of KVN, copying nothing and allocating
  nothing; written by `myb_oem_write_as`, `myb_opm_write_as`, `myb_omm_write_as`. Unit
  attributes, comments, the five standard entities; CDATA sections and a document type with
  an internal subset (added later — see "Old debts").
- **Covariances** of OEM (`myb_oem_covariance`, several matrices a segment, both ways of
  writing the lower triangle), OPM and OMM (`has_covariance`, `cov_ref_frame`, `covariance`):
  written and read in both notations the numbers come back to 1 ulp of the conversion
  km² ↔ m².
- **OCM** (`myb_ocm`, version 3 of the standard): a header, metadata and blocks of seven
  kinds — trajectories, physical properties, covariances, manoeuvres, perturbations, orbit
  determination, user-defined. The standard defines some two hundred keywords and dozens of
  kinds of trajectories, so the message is given as it is: keywords by name, lines as text
  and as numbers in the units of the message; Cartesian trajectories also in SI
  (`myb_ocm_state`). Written by `myb_ocm_write` in both notations. `test_ccsds` holds 536
  checks in all.
- **SPK segments written word for word** (`myb_spk_writer_add_raw`): copies, parts and merges
  of files of any type. 44 segments of 14 types (1, 2, 3, 5, 8, 9, 12, 13, 14, 17–21) were
  copied: the library reads the copies bit for bit as the originals, **and so does CSPICE**
  (308 states). Types 1 and 21 — the difference arrays of JPL's integrator — cannot be
  written otherwise: only it produces them.
- **SPK type 10 from TLE** (`myb_sgp4_spk10_packet`, `myb_sgp4_spk10_words` in `ext/sgp4`):
  the words of the generic segment are those of `SPKW10` for 235 element sets (with the
  directory of epochs at 230 sets); the packet of a two-line set holds the numbers of
  `getelm`, the epoch within 19 µs (SPICE's TDB − TT has one term).
- **A binary PCK writer** (`myb_pck_writer`): Chebyshev polynomials for the Euler angles of
  an orientation provider, the angles followed from node to node. A uniform rotation comes
  back to 7·10⁻¹⁵; a librating frame at degrees 5, 9 and 13 to 8·10⁻⁷, 3·10⁻¹¹ and
  1.5·10⁻¹⁵ rad; **CSPICE reads the file** (`pxform` through a frame of class 2) to
  7·10⁻¹⁶.
- **Clocks and accuracy codes of SP3** (`myb_sp3_clock`, `myb_sp3_sigma`, `myb_sp3_accuracy`,
  `myb_sp3_flags`; the members `clocks`, `clock_rates`, `accuracy`, `sigmas`, `flags` of the
  writer): a file in the layout of the format description, written from what was read, is
  the same in all its 37 lines.
- Moved to phase 20: a CK writer — with the CK and SCLK reader, without which nothing in the
  library can check it.

**Packaging (phase 13)**

- **ELF symbol versions** (`cmake/myb.map`, `ext/<name>/myb_<name>.map`): a function carries
  the version of the release that added it — `MYB_1.0` (394), `MYB_1.1` (5), `MYB_1.2` (202);
  the six support functions for `ext/` are `MYB_PRIVATE`. A program built against a newer
  release is refused when it is loaded, not at the first call. `test_exports` compares the
  version of every symbol with the map; a name the library does not have is refused by the
  linker.
- **pkg-config files for the libraries of `ext/`**: `myb_sgp4.pc`, `myb_nrlmsise00.pc`,
  `myb_troposphere.pc`, `myb_ionosphere.pc` (`Requires: myb`).
- Moved to phase 21: vcpkg and Conan packages — a port of a registry refers to the archive
  of a release by its hash, and 1.2 has no tag yet; a recipe that its own manager has not
  built would be unverified. macOS likewise.

**Reproducibility**: a workload of the debts (forces, shadows, a mesh, geometry, DOP853,
files and messages in both notations, the Jacobians of the surface terms and of the
Yarkovsky term); golden value `0xb76153a8dca45827` at -O0 and -O2;
`MYB_REPRO_PREFIX=17` checks the former `0x99a2ea1ff9d038a9` — no earlier bit changed.

### Phase 17 decisions

- **Light on a moving object** (`myb_force_radiation_moving`; `myb_radiation_term.moving`):
  the terms of the first order in v/c — aberration and the Doppler shift (Burns, Lamy, Soter
  1979, eq. 2); the part along −v is the Poynting–Robertson drag. A grain with β = 0.2 on an
  orbit of e = 0.3 over 20 revolutions: a falls by 7.272665·10⁻³ au, e by 2.190050·10⁻³; the
  averaged equations of Wyatt and Whipple integrated over the same time give **the same seven
  digits**. The Jacobian of the moving term is analytic (2·10⁻¹⁰ and 10⁻⁹ from differences).
- **Comets and asteroids** (`myb_nongrav`, `myb_term_nongrav`): a = g(r)(A1 r̂ + A2 t̂ + A3 n̂)
  with g(r) = α (r/r0)⁻ᵐ (1 + (r/r0)ⁿ)⁻ᵏ. Ready-made: the sublimation of water ice of Marsden,
  Sekanina and Yeomans (g(1 au) = 1 + 2·10⁻⁹ with the ten-digit α of the JPL catalogue) and a
  power law — the Yarkovsky effect as it is fitted to asteroids (A2, exponent 2). The drift
  with A2 = 10⁻¹⁰ m/s²: 9.681642·10⁻⁴ m/s against 2A2(1−e²)(r0/p)²/n = 9.681647·10⁻⁴
  (5·10⁻⁷); a comet with A1 and A2 — 5·10⁻⁶ from Gauss's equation averaged over the orbit by a
  quadrature.
- **The Yarkovsky drift from thermal properties** (`myb_yarkovsky_drift`): Vokrouhlický's
  linear theory for a body much larger than the depth of the thermal wave (κ₁ = κ₂ = κ₃ = ½):
  the diurnal and the seasonal terms. A kilometre of rock at 2.5 au — 2.2·10⁻⁴ au/Myr,
  outwards for a prograde rotation; the signs, the limits, the scaling with the radius and the
  maximum at a thermal parameter of √2 are checked.
- **The recoil of radiation** (`myb_force_emission`, `myb_term_recoil`): the heat of an RTG
  and of the electronics, the beam of an antenna — a power in a direction fixed in the body of
  the object, with an exponential decay (a radioactive source, outgassing that dies away).
  63 W on 241 kg — 8.72·10⁻¹⁰ m/s² (the Pioneer anomaly: 8.74·10⁻¹⁰).
- **The lagging tide** (`myb_force_tidal_lag`, `myb_term_tidal_lag`): Mignard's model with a
  constant delay — the acceleration of a satellite relative to the body from the tide it
  raises itself; the rotation of the body comes from its orientation.
  da/dt = 6 k₂ Δt (m/M)(R/a)⁵ n a (Ω − n): for the Moon (k₂ = 0.3, 600 s) **3.79 cm a year**;
  ten revolutions integrated — 28.348358 m against 28.348364 m of the formula; with a
  retrograde spin the orbit shrinks; a synchronous orbit — nothing.
- **The surface of the object** (`myb/surface.h`):
  - an **attitude law** — a function giving the matrix "inertial axes → body of the object"
    from the time and the state; ready-made: constant, orbital (radius, transverse, normal),
    the yaw steering of a navigation satellite (z to the body, y across the Sun), a uniform
    rotation;
  - **free-molecular flow** (`myb_force_gas_plates`): Sentman's formulas as Doornbos writes
    them — drag and lift of every plate, those turned away from the flow included, with the
    accommodation of energy and the wall temperature. Plates laid over a sphere by a
    quadrature give **the closed form of the sphere to 5·10⁻¹⁵** and no lift (to 2·10⁻¹⁵) at
    speed ratios from 0.3 to 31 — the two formulas check each other (their agreement to 10⁻³⁰
    was verified in mpmath before coding); C_D of a sphere at 7.6 km/s in atomic oxygen at
    900 K — 2.41;
  - **light on a mesh with self-shadowing** (`myb_force_radiation_mesh`): every face is a
    plate; with `shadows = k` a face is divided into k² equal triangles, and a part is lit
    only if the ray from its centre to the source does not meet the mesh
    (`myb_polyhedron_ray`, phase 16). A cube before a slab of 3×3 m: the lit area is
    **9.000000000000 m²** — exactly the silhouette; a cube as a mesh and as six plates —
    10⁻¹⁶; a faceted sphere with Lambertian re-emission — 0.995 of (1 + 4/9)πR²p with 1280
    facets;
  - **YORP**: the same function with a purely thermal surface is Rubincam's
    zero-conductivity model, and it returns the torque too. The torque along the axis averaged
    over a rotation: 6·10⁻¹⁷ for a symmetric body, −6.2190·10⁻³ p R³ for a body with "vanes",
    **+6.2190·10⁻³** for its mirror image — the sign changes, the size agrees to 10⁻¹²;
  - terms of a model: light on plates or on a mesh with the attitude and the shadows of the
    bodies (`myb_term_surface_radiation`: a black plate facing the Sun = the sphere with
    c_r = 1 to 4·10⁻¹⁶), the gas on plates with a wind (`myb_term_surface_drag`).
- **Winds** (`myb_wind_function`, the field `wind` of the drag term): the velocity of the gas
  relative to the rotating body; a wind "to the east in proportion to the distance from the
  axis" equals a faster rotation of the atmosphere to 10⁻¹⁴.
- **A mass that changes** (`myb_object.mass_function`): the terms of the surface take the mass
  at the instant t; `myb_thrust_mass_at` ties it to the thrust arcs — the light pushes harder
  what is left after the burn (the ratio of the accelerations = the ratio of the masses to
  10⁻¹⁴).
- **Impulsive manoeuvres** (`myb_orbit_kick`): the velocity jumps at the time reached and the
  integrator starts again. The variational equations are linear in Φ and S, so they continue
  from the same values (`myb_internal_variational_restart`), with nothing to compose. A
  Hohmann transfer by two impulses: the final orbit 1.5·10⁻⁸ m from the computed one,
  e = 6·10⁻¹⁶; Φ across both impulses against differences of the whole manoeuvre — 2·10⁻⁹ of
  the largest element.
  - Found: the automatic first step of IAS15 is made for a fresh problem, where the
    components that are zero are exactly zero. Continuing with a given Φ, a block whose
    velocity merely passes through zero (10⁻¹⁹ instead of 0) gives a "timescale" of 10⁻¹¹ s
    and `MYB_E_STEP_MIN`. On such a restart the step is now inherited from the integration
    that was interrupted (one that is too long is rejected as usual).
- **Reproducibility**: the workload of phase 17 (every new function and term, an orbit with
  two impulses and Φ); golden value `0x99a2ea1ff9d038a9`; `MYB_REPRO_PREFIX=16` checks
  `0x66d97d643b636e28`. No earlier bit changed (with a constant mass and no wind included).
- **Not done**: the refraction of light in the atmosphere of the occulting body (Vokrouhlický
  et al. 1993); the Yarkovsky theory for a finite depth of the thermal wave and on eccentric
  orbits; YORP with conduction; outgassing of a comet delayed with respect to the perihelion;
  the transitional flow regime; analytic Jacobians of the new terms (they come from
  differences of the model); the change of the rotation of the body under the torque of the
  tide; an offset of the centre of mass; mergers and collisions of bodies — in phase 18. Closed later: see "Debts of phases 12–17".
  Mergers and collisions stay with phase 18.

### Phase 16 decisions

- **Breaks of the step** (`myb_force_model_add_break`, `myb/model.h`). At the edge of a shadow
  the force is continuous but not smooth, and a step of IAS15 that straddles the edge loses
  its order: the error stops falling with ε. A break is a function g(t, y) of the same form as
  an event; it is declared to the **force model**, not to the propagator, so `myb_orbit`, the
  estimator and the filter all see it. After a step the break functions are examined on its
  dense output (as events are: four parts of the step, the root on the step's polynomial); a
  step that contains a sign change is **withdrawn and taken again from its beginning exactly
  up to the root** (`myb_internal_ias15_retake`): the prediction is the polynomial of the
  withdrawn step on the part kept, so the second attempt converges in a sweep or two, and the
  step after the edge starts from the same polynomial continued past it. The integrator is not
  restarted: the state, Φ and S simply continue. The sign of every function is remembered, so
  an edge is not found twice.
  - The first version — a cold restart of the integrator at every edge, with Φ and S composed
    across — gave the same accuracy but cost +70 % of the time of the LAGEOS fit (352 restarts
    per propagation, each with rejected first steps); withdrawing the step costs +26 %.
  - Check: an orbit of 12 000 km with 0.2 m²/kg, two days, 54 edges of the shadow. The spread
    of the final position over ε = 10⁻⁹…10⁻¹²: **6.4 m without breaks, 8·10⁻⁷ m with them**.
    Φ across eight edges against differences of whole propagations — 7·10⁻⁸ with a largest
    element of 25; ∂x/∂c_r — 2·10⁻⁷ m of 49 m.
- **LAGEOS does not depend on the integrator's precision**: the edges of the shadows of the
  Earth and the Moon are declared to the model (`myb_term_radiation_breaks`), 176 edges in the
  week, 2290 steps instead of 2177. The residual of the first guess — **25.9587 m at
  ε = 10⁻¹⁰, 10⁻¹¹ and 10⁻¹²** (without breaks 25.9520 m at 10⁻¹⁰ and 25.9582 m at 10⁻¹²);
  after the fit 2.05 cm at any ε. The finding of phase 14 is closed.
- **Ready-made events** (`myb/events.h`) — functions for `myb_orbit_event` with a context
  structure: the rate of the distance to a body (apsides, close approaches), a distance, a
  plane (nodes), the altitude over the ellipsoid of a body, the elevation over a station, the
  edge of a shadow (`myb_event_eclipse` on `myb_shadow_margin`: the disks touching from outside
  and from inside — the very instants where the force of light stops being smooth), the
  surface of a polyhedron. Bodies come from the environment by index, `MYB_EVENT_ORIGIN` is
  the origin itself. A function cannot return a status: when a provider fails it gives NaN,
  which is no sign change.
  - Against closed forms on an ellipse of e = 0.3: apsides 2·10⁻¹¹ s, nodes and a sphere
    10⁻¹² s; the altitude at the event — 10⁻⁹ m from the threshold, the elevation — 2·10⁻¹⁴ rad
    (against `myb_measure_azel`); the impact on an octahedron — 2·10⁻⁹ m from its face.
  - **Eclipses against CSPICE N0067**: the SPK is written by the test itself
    (`myb_spk_writer`, phase 15), `gfoclt_c` finds the windows of the occultation of the Sun by
    the Earth from the object (convergence 10⁻⁹ s): 28 contacts of the penumbra and the umbra
    in a day — **within 4.4·10⁻¹⁰ s**.
- **The geometry of a surface** (`myb/geometry.h`), in the body-fixed frame:
  - a ray and a triaxial ellipsoid (`myb_ray_ellipsoid`): in coordinates divided by the
    semi-axes, through the point of the line nearest the centre — a discriminant without the
    subtraction of close numbers; from inside, the point where the ray leaves. Against
    `surfpt_c` on five shapes, 300 rays: **5·10⁻¹⁶ of the distance**;
  - the nearest point and the altitude (`myb_ellipsoid_nearest`): the largest root of
    Eberly's equation by Newton's method from a lower bound — the function is convex and
    falling, the iteration is monotone and never passes the root; inside the evolute, in a
    plane of symmetry, there are two nearest points (or a circle) — the one on the positive
    side of the shortest axis is given. Against `nearpt_c`, 300 points, interior ones
    included: **4·10⁻¹⁶ of the size, the altitude 3·10⁻¹⁶**;
  - the normal, the angles of illumination (phase, incidence, emission);
  - a polyhedron: a ray (Möller–Trumbore over every face, faces widened by 10⁻¹² so that a ray
    through a shared edge finds no gap) and the solid angle (4π inside, 0 outside). Against an
    independent solution (planes of the faces and areas) on a convex and a dented icosahedron:
    8000 rays, not one disagreement, distances 2·10⁻¹⁵.
  - Through kernels against `sincpt_c`, `subpnt_c` (both kinds) and `ilumin_c` at 20 epochs:
    points 10⁻⁸ m (the limit of the SPK interpolation), altitude 3·10⁻⁹ m, angles
    **7·10⁻¹⁶ rad**.
- **Reproducibility**: the workload of phase 16 (geometry, an orbit through a shadow with
  breaks and partials, every ready-made event); golden value `0x66d97d643b636e28`;
  `MYB_REPRO_PREFIX=15` checks `0x9128edc3c0c93fb3`. Adding the withdrawal of a step to IAS15
  changed no earlier bit.
- **Not done**: predicting an edge before the step (each edge now costs one withdrawn step);
  an accelerating structure for rays on a mesh (every face is examined); limb and terminator
  points; geometry corrected for light time; the approach of two propagated objects — in
  phase 19. Closed later: see "Debts of phases 12–17". Predicting the edge was tried and rejected; the approach
  of two objects stays with phase 19.

### Phase 15 decisions

- **Printing of numbers** (`myb_format_double`, `myb_format_fixed`, `myb_format_scientific`;
  `myb/io.h`) — the inverse of `myb_parse_double` and as exact: every digit comes from integer
  arithmetic on the exact value of the number (big integers, digit generation after Steele and
  White), without `printf`, the C runtime or the locale. Three forms: the shortest decimal
  that `myb_parse_double` reads back to the same number (as `repr` in Python); a given number
  of decimals (up to 60); a given number of significant digits with an exponent (up to 60) —
  the exact value rounded to nearest, ties to even. Against CPython on 6124 values — **not one
  difference in any form**; 200 000 random numbers read back to the same bits.
- **Writing SPK** (`myb_spk_writer`, `myb/kin.h`): segments of types 9 and 13 from states at
  epochs (Lagrange and Hermite), of types 2 and 3 — Chebyshev polynomials through the states
  of a provider at the Chebyshev nodes of equal parts of the interval, with a report of the
  largest deviation between the nodes (to choose the degree and the number of parts). The DAF
  file: file record, a chain of summary and name records (25 segments each), data;
  little-endian; into a buffer or a file (`myb_platform_write_file` — the only new call of the
  OS). A file of 34 segments is read by **CSPICE N0067: positions 3.9·10⁻¹⁶, velocities
  1.1·10⁻¹³** from the library's own reader; the states of types 9 and 13 come back at their
  epochs within 1.8 ulp of the length of the vector; type 2 of degree 12 over 1000 s —
  6·10⁻⁸ m from the orbit.
- **CCSDS messages** (`myb/ccsds.h`; CCSDS 502.0-B-2/B-3, the keyword = value notation): OEM —
  reading (header, segments with metadata, states and accelerations, both forms of the epoch;
  covariance sections are passed over; mandatory keywords are checked; an error names the
  line) and writing; OPM and OMM — parsing into a structure and writing. In the file km, km/s,
  degrees, km³/s², in the API — SI; the mean elements of an OMM are kept in the units of the
  catalogue, as written: they mean something only with their theory. **A number is written as
  the shortest decimal in the file's unit that the reader turns back into the same SI value**:
  of the quotient and its two neighbours the shortest is taken among those that the
  multiplication returns to the original number. So a message read and written keeps its
  numbers as they were ("7000" stays "7000", not "7000.000000000001"), and one written and
  read gives the states back exactly — except those that no number of kilometres times 1000
  reaches: 183 of 12 000 random ones (1.5 %), which differ by 1 ulp. Seconds of epochs get as
  many decimals as they need, three at least; epochs come back bit for bit.
- **SP3** (`myb_sp3_write`, `myb/kin.h`): SP3-c and SP3-d from a description (satellites,
  epochs, positions, velocities, gaps); the GPS week, the seconds of the week and the MJD of
  the header — from the first epoch. The ILRS excerpt of LAGEOS-1 read and written: **258 720
  characters of epoch, position and velocity lines — character for character**, as are the
  lines of the header (lists of satellites, accuracies, types); the reader gives the same bits
  from the written file.
- **TLE** (`myb_tle_format`, `myb_tle_from_omm`, `myb_omm_from_tle`; `ext/sgp4`): the two
  lines with their checksums, Alpha-5 above 99 999; a value that does not fit its field —
  `MYB_E_DOMAIN`. All 33 sets of Vallado's verification file come back to the same values bit
  for bit, 31 character for character (in two the file itself departs from the catalogue's
  notation: no ephemeris type; a zero B* written as `00000-0`). OMM ⇄ TLE: the numbers of the
  catalogue pass unconverted, the epoch — calendar date ⇄ day of the year.
- **Tables of numbers** (`myb_table_*`, `myb/io.h`): reading of CSV and its relatives (comma,
  semicolon, tab, blanks; a row of names; `#` comments; an empty field is a NaN) and writing
  with the shortest decimals: 3500 random numbers of every magnitude — bit for bit.
- **The writers do not allocate**: the text goes into the caller's buffer, the full length is
  always reported (`MYB_E_OVERFLOW` — call again with a buffer of that size); the internal
  writer is `src/io/text.[ch]`. Only the SPK writer and the readers take memory.
- **Reproducibility**: the workload of phase 15 hashes the bytes themselves — texts of
  numbers, the image of an SPK, messages, SP3, TLE lines, a table — and the numbers read back;
  golden value `0x9128edc3c0c93fb3`; `MYB_REPRO_PREFIX=14` checks `0xb110642ae757e6e0`.
- **Not done**: the XML notation of the CCSDS messages; covariances of OEM and OPM; OCM;
  writing SPK of other types (1, 21, 10), CK and PCK; clocks and accuracy codes in SP3. Closed later: see "Debts of phases 12–17".
  The CK writer moved to phase 20.

### Phase 14 decisions

- **Gradient of a field of spherical harmonics** (`myb_gravity_field_gradient`, `myb/grav.h`):
  the second derivatives of the potential in the body-fixed frame. The second derivatives of
  a solid harmonic of degree n are solid harmonics of degree n + 2 (Cunningham 1970;
  Montenbruck & Gill, eq. 3.33): the operators D± = ∂x ± i∂y and Dz take (n, m) to
  (n + 1, m ± 1) and (n + 1, m). The five independent components (V_zz, V_xz, V_yz,
  V_xx − V_yy, V_xy; V_xx + V_yy = −V_zz) are series of the same form as the potential with
  coefficients of degree n + 2 and orders m, m ± 1, m ± 2; the factors of the fully normalised
  functions are √((2n+1)/(2n+5)) times the root of a product of four numbers n ∓ m + k. The
  series are summed by the same column recursion of Holmes and Featherstone as the potential:
  **no division by cos φ, the pole is an ordinary point**, and the accuracy is that of the
  potential (the recursion tables of a field reach two degrees further; the degree limit of
  the gradient is 2698). Against differentiation of the series in mpmath (60 digits; a field
  of degree 8, 32 cases at all latitudes and at the poles themselves, several ranges of
  degrees): **2.1·10⁻¹⁵ of the largest component, 6·10⁻¹⁶ at the poles**; against central
  differences of the acceleration in a field of degree 60 — 2·10⁻¹⁰.
- **Analytic Jacobians of the terms** (`myb/env.h`, `myb/rel.h`): the field (the gradient
  turned into the inertial axes; the term brings the coefficients to the instant itself),
  radiation pressure (exact in the light and in the umbra, differences of the term itself in a
  penumbra), drag (exact in the velocity and the rotation of the atmosphere, the gradient of
  the density by differences of the density function), Schwarzschild, Lense–Thirring, de
  Sitter (`myb_rel_*_jacobian`). Against central differences: field 1·10⁻¹⁰, light 6·10⁻¹¹,
  drag 4·10⁻⁸ and 2·10⁻¹¹ (position and velocity), relativistic terms up to 7·10⁻¹¹ (the limit
  of the differences themselves). `myb_term_no_jacobian` is a Jacobian of zeros for terms
  whose gradient is left out on purpose; `myb_force_model_bind_parameter` binds a parameter to
  one term: its partial is the difference of that term, not of the whole model (10⁻¹³ instead
  of the 10⁻⁵ of a difference of the sum with the field, at a fraction of the cost).
- **State transition matrix with a field**: in a rotating field of degree 8 over two
  revolutions the matrix from the analytic gradient is **symplectic to 6·10⁻¹⁷** (10⁻¹² with
  the Jacobian by differences) and agrees with differences of whole propagations to 5·10⁻¹⁰.
- **Estimation** (`myb/estim.h`). The unknowns are the state at the epoch and the parameters
  of the force model. A measurement is a function of the state of the object at the time of
  the observation with its derivatives; ready-made: the position and the state of the object
  itself, range, range rate, right ascension and declination, azimuth and elevation from a
  station on a body of the environment (`myb_station`: geodetic coordinates on an ellipsoid,
  the axes east–north–up). Derivatives of the measurements against differences — up to
  3·10⁻¹⁰.
  - **Least squares over an arc** (`myb_fit_orbit`): Gauss–Newton on the variational
    equations, normal equations of the free unknowns scaled by their diagonal and factorised
    by Cholesky; fixed unknowns, a priori sigmas, covariance, residuals; an arc forward or
    backward from the epoch. Without noise a displaced first guess (hundreds of metres, c_d
    off by 0.5) returns to the truth to **2·10⁻⁹ m and 6·10⁻¹¹ in c_d**; with noise (σ: 1 m,
    1 mm/s, 2 µrad) the weighted RMS is 1.035 and the largest error 1.7 σ of its own
    covariance.
  - **Filter** (`myb_filter`): the extended Kalman filter on the same unknowns, the
    covariance carried by Φ and S, process noise as a white acceleration (q·Δt³/3 — checked:
    2.166 m² against 2.16), the update in Joseph's form. The same noisy observations one by
    one from a covariance of a kilometre: 0.68 m from the truth at the end with σ 0.43 m, the
    largest error 1.9 σ.
- **LAGEOS by the library's estimator**: the test's own Gauss–Newton with partials by
  differences and twelve parallel propagations is removed; `myb_fit_orbit` on the variational
  equations (analytic Jacobians of the field with tides, the bodies, relativity and light;
  the gradients of the Earth's radiation and of the empirical terms, 10⁻⁹ of the field's, left
  out) gives **2.05 cm** in two corrections, 21 s instead of the former ~35 s of fitting (the
  whole test 33 s against 44 s). **The spacecraft of phase 8** (`test_spacecraft`) through the
  same estimator — the same numbers: Roadster 306 m, OSIRIS-REx 12.6 m, JWST 260.0 m.
- **Found: the integrator steps over the penumbra.** Before the fit the residual of LAGEOS is
  25.952 m (variational integration, ε = 10⁻¹⁰), 25.954 m (plain), 25.959 / 25.958 / 25.959 m
  for ε = 10⁻¹¹…10⁻¹³; with the shadows switched off the result does not depend on ε
  (25.417 m in both cases). The entry into the penumbra lasts seconds, a step minutes: the
  step criterion built on a, a′, a″ does not see the kink, and a few millimetres along the
  track accumulate in a week. The cure — breaking the step at the edges of a shadow — was
  made in phase 16 (see above): with it the residual is 25.9587 m at any ε.
- **Reproducibility**: workload of phase 14 (gradient, Jacobians of relativity and of a model,
  a fit and the filter on noisy observations); golden value `0xb110642ae757e6e0`;
  `MYB_REPRO_PREFIX=13` checks `0x6e21db689ceda9e6` — no earlier bit changed.
- **Not done**: station displacements (solid Earth tides, ocean loading, pole tide — IERS
  ch. 7) and the light time in the measurements — moved to phase 20; biases of the
  measurements as unknowns; a smoother; analytic Jacobians of the radiation of a planet and of
  the empirical terms. Closed later: see "Debts of phases 12–17".

### Phase 13 decisions

- **Installation** (`cmake --install <build> --prefix <where>`; `cmake/MybInstall.cmake`,
  option `MYB_INSTALL`): the headers `include/myb/` and `include/myb_<ext>.h`, `libmyb.a` and
  the libraries of `ext/`, the data kernels `share/myb/data/` (without `download`), licences
  and notices. CMake package: `find_package(myb 1.2)` gives the targets `myb::myb`,
  `myb::myb_shared`, `myb::sgp4`, `myb::nrlmsise00`, `myb::troposphere`, `myb::ionosphere`
  (and their `_shared`) and the variable `myb_DATA_DIR`; compatible within a major version.
  `myb.pc`: the prefix is computed from the place of the file itself (`${pcfiledir}`), so an
  installation can be moved; a static installation has `-DMYB_STATIC` in `Cflags` (on Windows
  the declarations would otherwise be `dllimport`).
- **Shared library** (`MYB_BUILD_SHARED`, on in the `*-release` presets): the same sources are
  compiled a second time with hidden visibility, and only what is marked `MYB_API` leaves the
  binary. The static library stays the main one: the tests reach into internal functions. The
  libraries of `ext/` are separate shared libraries on top of `libmyb` (`MYB_EXT_API`,
  `MYB_EXT_BUILDING`); the six support functions they are built on (allocator, diagnostics,
  `MYB_CHECK`) are exported for them and are not part of the API. The CORE-MATH objects have
  hidden visibility and serve both libraries.
- **Examples** (`examples/`, public headers only, C17): `low_orbit` — a low orbit about an
  oblate planet with an atmosphere from the ready-made terms: 15 ascending nodes in a day, node
  drift −5.013°/day against −5.002 of the first order in J₂, the period shorter by 0.239 s
  against 0.243 of the first order in the drag; `time_scales` — one instant on TAI, TT, TDB,
  GPS from the kernels of `myb_DATA_DIR`. Built in the tree (`MYB_BUILD_EXAMPLES`) or as a
  separate project against an installation — with any C17 compiler (only the library itself
  needs Clang).
- **Checks**: `test_install` and `test_install_shared` install the build into a scratch
  prefix, configure `examples/` through `find_package(myb)`, build and run them (each example
  returns 0 when its own check holds); `test_exports` and `test_exports_<ext>` read the export
  table (`llvm-readobj --coff-exports` for a DLL, `llvm-nm -D` for ELF): `libmyb.dll` — 442
  symbols, all `myb_`, of the internal ones only the six allowed; the libraries of `ext/` —
  1 to 9. Under the sanitizers the installation check is not registered (an outside project
  would have to link their runtime).
- **CHANGELOG** (`CHANGELOG.md`, `CHANGELOG.ru.md`): versions 1.0.0.0, 1.1.0.0 and the
  current one.
- **Gentoo** (2026-11-10, Clang 21): `linux-release` with the static and the shared libraries
  — 78/78 (5 tests skipped: the data sets were not downloaded); `libmyb.so.1.2.0.0` — the same
  442 symbols, the libraries of `ext/` — 1, 9, 8, 5; the examples against the installation
  print the same digits as on Windows. Found: the library directory there is `lib64` — the
  installation check gets the install directories from the build instead of assuming `lib`.
- **Not done**: vcpkg and Conan packages; ELF symbol versions; pkg-config files for the
  libraries of `ext/`; macOS (phase 21). Closed later: see "Debts of phases 12–17". The vcpkg and Conan packages
  moved to phase 21.

### Phase 12 decisions

- **Why.** A force model (`myb/model.h`) summed terms written by the user: `test_lageos` and
  `test_spacecraft` carried 700 lines of glue each (orientation, third bodies, radiation). From
  a binding that would be a call into the interpreter for every evaluation of the
  acceleration. A model is now assembled from structures of numbers and pointers to providers.
- **Providers** — the common form in which ephemerides and rotation enter the models:
  `myb_state_function` (moved from `observe.h` to `myb/kin.h`; `myb_spk_body` — a body of an
  SPK set) and `myb_orientation_function` (`myb/frames.h`: the rotation inertial axes →
  body-fixed frame and, for a body whose rotation is observed, its rotation time UT1 and the
  coordinates of its pole, which the tides need; ready-made: CIP + EOP, binary PCK, WGCCRE, a
  constant rotation).
- **Environment** (`myb/env.h`, `myb_env`): the bodies of a scenario, each with a GM, a state
  provider, an orientation provider, a radius, a luminosity, limb darkening, a spin; all
  optional (the body model of §3 as far as the forces need it). The body without a state
  provider is the origin (the central body: third-body terms get their indirect part); when
  every body has one, the origin is inertial (a barycentre) and there are no indirect terms.
  The values at an instant t are computed once and kept until t changes: an integrator
  evaluates every substep several times, and six terms need the Sun.
- **Time scales.** The model time is seconds from an epoch on a named scale; each provider
  gets the same instant on its own scale (TDB for the JPL ephemerides, TT for the rotation of
  the Earth). The difference of two scales is read from a table with a step of the model time
  (one hour by default) at the nodes epoch + k·step, Lagrange interpolation on four nodes
  written about the node before t (a large constant part of the difference does not cancel
  against itself). The value is a pure function of t: the same bits whatever was asked
  before. Against the exact conversion over ±40 days: **TDB−TT 3.5·10⁻¹⁶ s, TCB−TT (20 s with
  a rate) 3.1·10⁻¹⁵ s**; the series of 787 sines is evaluated once an hour instead of at every
  evaluation (the "accelerated TDB−TT series" of the phase 7 list).
- **Ready-made terms** — a context structure and a function of the form `myb_term_function`:
  point mass and third body (with analytic Jacobians), the field of a body in its frame with
  its time dependence and tides (for a body away from the origin, minus the same field at the
  origin), radiation pressure of a source with the shadows of bodies, radiation of a planet,
  drag in an atmosphere rotating with the body, Schwarzschild, Lense–Thirring and de Sitter,
  empirical RTN (constant and once per revolution). The object is a `myb_object` (mass, area,
  cr, cd): the parameters to estimate are its fields and those of the terms' structures. An
  atmosphere enters as a density function `myb_density_function`; the core has the layered
  exponential one (`myb_exponential_atmosphere`), the numbers of the layers being data.
- **Propagation** (`myb/orbit.h`, `myb_orbit`): IAS15 from t0 towards t_end with its own
  steps, every state from the dense output of the step that contains its time. The steps do
  not depend on the times asked for (only the last is cut, at t_end), so **the state at a time
  t is the same bits whether 500 other states and events were asked on the way or none**.
  Events are roots of g(t, y): the sign is examined at the ends of the quarters of a step, a
  root is located on the dense output and moved to the first representable time where g has
  changed sign, so that the next search from the root does not find it again; a direction
  filter. On request the state transition matrix and the sensitivities (the variational
  equations of `myb/uncert.h`).
- **Checks.**
  - `test_env` (2774 checks, no downloaded data): every term against the same functions of
    `grav`/`force`/`rel` called by hand with the providers' values — **bit for bit** (a term
    is no more than that composition); Jacobian of the third body against differences
    2·10⁻¹¹; the field of a body away from the origin, degree 0, against Battin's form
    5.5·10⁻¹⁶; drag against −½ρ c_d A/m |w| w, w = v − ω × r — 1.4·10⁻¹⁶; fraction of the disk
    in the shadow: 1 / 0 / 0.517 / 1 (day, umbra, penumbra, past it); a provider is called
    once per instant; errors.
  - `test_orbit`: two bodies, 5 revolutions with e = 0.3 — 4·10⁻⁸ m from the closed form;
    apsides of 4 revolutions against Kepler's equation — 3.6·10⁻¹¹ s, forward and backward,
    with the filter; state transition matrix against the closed form — 2.6·10⁻¹⁵; decay of a
    circular orbit in an atmosphere over 10 revolutions −25.408 m against −25.403 m of the
    first order.
  - **LAGEOS** (`test_lageos` moved to the ready-made terms and `myb_orbit`; the former
    assembly by hand is kept as the reference): acceleration of the model of ready-made terms
    against the one by hand — **3.5·10⁻¹⁶** at 42 states along the orbit; before the fit the
    same 25.954 m and 2228 steps, after it **the same 2.01 cm** (R 2.5 mm, T 1.83 cm,
    N 8.0 mm); the test takes 44 s against 42 s.
- **Reproducibility**: workload of phase 12 (table of scales, a model of every term with its
  Jacobian and a parameter, states, an event, the state transition matrix backward in time);
  golden value `0x6e21db689ceda9e6`; `MYB_REPRO_PREFIX=11` checks `0xc9b659a710498026` — no bit
  of phases 0–11 changed. Version 1.2.0.0 (no tag).
- **Verification**: Windows (Clang 23) — `windows-release`, `windows-dev` (the new and the
  affected tests, the same hash at -O0 and -O2), `windows-ubsan` (`test_env`, `test_orbit`,
  `test_repro`); Gentoo (Clang 21) — `linux-release` 78/78, `test_env`, `test_orbit`,
  `test_repro`, `test_sgp4` under `linux-asan` (ASan + UBSan + LeakSanitizer) — clean, same
  hash. `test_lageos` was not run on Linux (it needs the downloaded data) — left to CI.
- **API changes**: `myb_state_function` is declared in `myb/kin.h` (`observe.h` includes it —
  user sources do not change); `myb/force.h` includes `myb/frames.h`.
- **Not done**: plates with the attitude of the spacecraft and analytic Jacobians of the
  field, radiation, drag and relativity (phases 17 and 14); the indirect term of the
  oblateness of the central body; a copy of an environment for threads (now: one environment
  per thread); other integrators in `myb_orbit`; the exact shadow when two bodies cover the
  disk at once (now: the product of the fractions); an adapter of NRLMSISE-00 to
  `myb_density_function` (phase 20). Closed later: see "Debts of phases 12–17".

### Phase 11 decisions (after 1.0.0.0)

- **3PN and 3.5PN** (`myb_pn_binary_evaluate`, `myb/rel.h`): coefficients A, B of the relative
  acceleration from the Blanchet review (LRR, 2024 revision, eqs. calAcoeff, calBcoeff); the
  gauge length r0′ in the 3PN logarithms is a parameter (`myb_pn_binary.r0`). The old function
  `myb_pn_binary_acceleration` is retained (without 3PN and spins). Checks:
  - **coefficient transcription** — against machine parsing of the review's LaTeX source
    (arXiv:1310.1528, sympy, mpmath 40 digits, 40 states): all orders 1PN…3.5PN to 4·10⁻¹⁶ — a typo
    made in manual transcription would have shown up;
  - **generalized Kepler law** Ω²(r) to 3PN with ln(r/r0′) (eq. keplerlaw): the 3PN coefficient by
    Richardson extrapolation to 10⁻⁶ at ν = 0.25 and 0.1, ln = 0 and 2;
  - **radiation reaction on a circular orbit** (eq. aieom without the 4PN tail): the correction
    −743/336 − 11/4 ν — to 4·10⁻⁴ (residual O(γ));
  - **inspiral rate** ṙ (eq. rOmdot): time for the orbit to shrink 200 → 190 Gm/c² against the
    formula with the 1PN correction −1751/336 − 7/4 ν — 98.7% of the correction (the residual is the
    2PN term);
  - **invariant K** at 3PN (Damour–Jaranowski–Schäfer): 135/2 + (−649/4 + 123/32 π²)ν + 7ν² —
    to 0.06% at ν = 0.25 and 0.1; independent of r0′ (×1000) to 10⁻⁴ (whereas the Kepler law changes
    by ~38 x³) — all logarithmic terms are consistent.
- **Spins** (Kidder 1995; covariant spin supplementary condition): spin-orbit acceleration (1.5PN),
  spin-spin (2PN), rotation-induced quadrupoles (Poisson 1998; κ = 1 for Kerr), and
  spin precession; the spins are dimensionless χ = cS/(Gm²) and part of the state (for IAS15 — as
  “velocities” of auxiliary coordinates). Checks: test-particle limit — the Lense–Thirring
  acceleration (`myb_rel_lense_thirring`) to 4·10⁻¹⁶; spin-spin + quadrupoles — a Hamiltonian
  system: E and J are conserved to 10⁻¹⁷; spin-orbit: Kidder's E_N + E_SO is conserved up to terms
  of second order in spin (×1/10 spin → ×1/100), L_N + L_SO + S — with a residual from the precession
  inside L_SO that falls off as 1/r (beyond the order of the equations); aligned spins do not precess
  (exact zero); **PSR J0737−3039B**: geodetic precession 5.0735°/yr against the Barker–O'Connell
  formula 1.3·10⁻⁵ and the prediction 5.0734°/yr (Breton et al. 2008).
- **Plasma** (`myb/observe.h`): electron column of a power-law profile Σ n_k (R/r)^p_k along a
  segment (adaptive Gauss–Kronrod in the angle seen from the center; radial segments — analytically)
  and group delay 40.308 N/f² (K = e²/(8π²ε₀mₑ) from CODATA 2022, new constants in
  `codata-2022.tk`). Against mpmath — 6·10⁻¹⁶, against the closed forms for r⁻², r⁻⁴ — 10⁻¹⁵.
- **Troposphere** (`ext/troposphere`, target `myb::troposphere`): IERS 2010 Ch. 9 models from the
  published equations — Saastamoinen (ZHD, ZWD), Mendes–Pavlis (optical), FCULa, VMF1,
  GMF, GPT, Chen–Herring gradients; GMF/GPT tables — data from the IERS software (SHA-256). Against
  the test cases of the IERS software: FCULa, VMF1 (with and without height), GMF, GPT — 10⁻¹⁶…10⁻¹⁵.
  Found: the published FCUL_ZD_HPA example is 1.96·10⁻⁶ smaller in both components — a common
  factor f(φ, H), as if the example had been computed with a height 7 m lower; executing the
  program's own operations (with its single-precision literals 0.00266, 3.759) gives the same ratio —
  the example does not match its own input data (4 µm of delay).
- **Ionosphere** (`ext/ionosphere`, target `myb::ionosphere`): IONEX maps (2D, reading per the
  fixed-field format, skipping RMS and height maps, 9999 → no data), interpolation per the IONEX
  description (bilinear, in time — between maps rotated with the Sun), single-layer model
  (SLM and CODE MSLM) with pierce point. Against an independent Python implementation —
  2·10⁻¹¹ (MJD resolution in double, 0.6 µs); pierce point against the intersection of the ray with
  the shell in Cartesian coordinates — 9·10⁻¹⁶.
- **Charging** (`myb/force.h`): OML currents to a sphere for any particle species with drift
  (Whipple 1981, closed forms via erf/erfc; for M < 10⁻³ — a series in M²), photoemission, floating
  potential (monotonic current: bracket expansion and bisection down to adjacent doubles). Against
  integration of the OML cross-section over a shifted Maxwellian (mpmath) — 6·10⁻¹⁶ for M from 0 to
  8; equilibrium in an electron–proton plasma — 1 − W(e√(m_p/m_e)) = −2.50385 kT/e to 10⁻¹⁵. erf and
  erfc — the fdlibm algorithm (FreeBSD msun, SunPro license), within 2 ulp, identical bits everywhere
  (libm erf is forbidden by `test_no_libm_imports`).
- **Reproducibility**: phase 11 workload (3.5PN binaries with spins, columns and delays,
  charging currents and potentials); golden `0xc9b659a710498026`; `MYB_REPRO_PREFIX=10` verifies
  `0xd1bce261642c87ac` — not a single bit of phases 0–10 changed. Version 1.1.0.0 (not tagged).
- **Check**: Windows — new and affected tests in `windows-release`, `windows-dev`,
  `windows-ubsan`; full suite — CI (run 42748579502, all 10 jobs green).
- **CI — on tags only** (and manually): push to main, pull request and the weekly run are
  disabled; `linux-asan-debug` is now in every run, the live EOP series are refreshed in each one.
- **Not done**: 4PN (with the hereditary tail) and higher-order spin-orbit terms, radiation
  reaction due to spins; binaries in the N-body problem; APG gradients, VMF3/GPT3, 3D ionosphere
  maps, higher-order ionospheric terms; secondary emission and a thick Debye sheath in
  charging. Assigned to phases or closed: see "Old debts (phases 3–11) and the bindings".

### Phase 10 decisions

- **Compact binaries beyond 1PN** (`myb/rel.h`, `myb_pn_binary_acceleration`):
  relative acceleration in harmonic coordinates of the center-of-mass frame up to 2.5PN
  (Blanchet 2014, eqs. 203–204; Kidder 1995): 1PN, 2PN, 2.5PN radiation reaction, orders
  selected by flags; the Peters–Mathews formula for Ṗ. Checks:
  - **circular-orbit invariant** K = Ω_φ/Ω_r = 1 + 3x + (27/2 − 7ν)x² (Damour, Jaranowski,
    Schäfer 2000; test-particle limit — Schwarzschild's (1 − 6x)^(−1/2)): (K − 1 − 3x)/x² measured
    on nearly circular orbits of the 1PN + 2PN equations differs from 27/2 − 7ν by 88·x at
    ν = 0.25 and 0.01, x = 10⁻³ and 3·10⁻⁴ — the residual is linear in x, the 2PN coefficient is
    confirmed;
  - **circular inspiral** — dΩ/dt against Peters 1.0006 at x = 5·10⁻⁵;
  - **PSR B1913+16** (masses 1.438 and 1.390, P_b = 0.322997 d, e = 0.617): Ṗ from the times of
    200 periastron passages with IAS15 — −2.40234·10⁻¹², Peters–Mathews for the eccentricity of the
    integrated orbit — ratio 1.000031; periastron advance 4.2256°/yr (observed 4.2266). Found: DOP853
    gives a spurious Ṗ = −3.5·10⁻¹⁵ (energy drift, 0.15% of the effect) — the measurement uses IAS15
    (no secular drift, 10⁻⁵); Newtonian initial conditions at pericenter change the eccentricity of
    the 1PN orbit by 4·10⁻⁵ ((v/c)² v²/|E| ≈ 10⁻⁴), which f(e) turns into 3.5·10⁻⁴ — the formula is
    fed the eccentricity derived from the radii of the integrated orbit.
- **Observables** (`myb/observe.h`): light-time equation (reception and transmission) with the
  Shapiro delay per Moyer 8-38 ((1 + γ) GM/c² terms inside the logarithm, body at the time the ray
  passes, no loss of significance near the limb), range, range rate, two-way
  measurements, aberration via the Lorentz transformation, light deflection (PPN form of eraLd),
  RA/Dec. Checks: against CSPICE N0067 on DE440 (`spkezr` CN and XCN, Moon, Sun, Mars, Jupiter from
  Earth, 1974–2026, 208 rows) — light time 4·10⁻¹² s, position 1.7 mm at 10¹² m, velocity
  5·10⁻¹⁰ m/s; SPICE first-order stellar aberration — 2.5·10⁻⁹ rad from the exact one; against ERFA:
  aberration bit for bit, deflection 9·10⁻⁹ (ERFA SRS constant); delay against quadrature of the
  1PN light ray — 9·10⁻¹⁵; at the solar limb, 1 → 9.5 au — 141.6 µs; range rate against
  differences — 3·10⁻¹¹; two-way Earth–Moon range — 2|x_M − x_E|/c + 2 delays to
  2.5·10⁻⁸ s (order (v/c)²).
- **Charge and plasma** (`myb/force.h`): magnetic fields of bodies via Gauss coefficients
  (Schmidt semi-normalization) through the gravity-field spherical-harmonic evaluator
  (P_n^m = P̄_nm/√(2n + 1), the same workspace), IGRF-format tables (epochs, secular variation),
  Lorentz force, corotation field −(Ω × r) × B, plasma flow pressure (solar wind). Checks: IGRF-14
  (NOAA table) against ppigrf — 1.4·10⁻¹⁵ |B| at 160 points 2000–2025; dipole g₁₀ against the
  formula — 8·10⁻¹⁶; helix in a uniform field over 50 turns — 5·10⁻¹⁵; E × B drift — exact; in a
  dipole the speed and the canonical momentum p_φ are conserved to 3·10⁻¹⁶.
- **Reproducibility**: phase 10 workload (binaries, light time and delays, two-way
  measurements, aberration and deflection, magnetic field, Lorentz, corotation, plasma flow); golden
  `0xd1bce261642c87ac`; `MYB_REPRO_PREFIX=9` verifies `0x827497389931b8a9` — not a single bit of
  phases 0–9 changed.
- **Final check** (2026-11-10): Windows (Clang 23) — `windows-dev`, `windows-release` 64/64,
  `windows-ubsan` 63/63; Gentoo (Clang 21) — all 7 presets: `linux-release`, `linux-dev` 64/64,
  `linux-asan`, `linux-ubsan`, `linux-lsan`, `linux-tsan`, `linux-asan-debug` 63/63; TSan stress —
  200 runs of the multithreaded tests with no reports; hash `0xd1bce261642c87ac` is the same
  everywhere. The binaries test was switched to correctly rounded functions (it used libm `cbrt`,
  `atan2`): its numbers are now identical on both systems to the last digit.
- **Version 1.0.0.0** (`MYB_VERSION_TWEAK` added; the number is major·10⁶ + minor·10⁴ + patch·10² +
  tweak). **CI** enabled: push to main, pull request, weekly run with a refresh of the live series.
- **Not done**: 3PN and higher, spin-orbit and spin-spin terms; binaries in the N-body problem
  (only the relative motion of two bodies); phase and group delay in plasma and troposphere;
  spacecraft charging models (the charge is set by the user). Assigned to phases or closed: see "Old debts (phases 3–11) and the bindings".

### Phase 9 decisions

- **Stack** (before the phase began): five frames of 17–34 KB eliminated — the Poisson-series
  phasor table is packed (≤ 128 entries, 2 KB; powers and number of series from the data, 5 and 5),
  `myb_covariance_transform` computes row by row with a single row of n, `myb_eigenvector_real` and
  `myb_shoot` take workspace from the caller (`myb_shoot_work_size`; the shooting limit raised to
  64 unknowns), `myb_timesys_create` without the 3.5 KB compound literal. The new test
  `test_stack_budget` parses the `-fstack-usage` output (flag `MYB_STACK_USAGE`, on by default outside
  sanitizers): every frame is static and at most 4 KB (8 KB at -O0); the largest is currently 3.4 KB
  (SPK loading). The golden bits did not change.
- **SGP4** (`ext/sgp4`, target `myb::sgp4`; the model is tied to Earth — WGS-72 constants,
  lunisolar terms with fixed coefficients — hence outside the core, like NRLMSISE-00):
  port of Vallado's 2006/2020 code via the MIT-licensed python-sgp4 translation, preserving the order
  of operations, with correctly rounded functions and an exact fmod (one FMA); propagation is a pure
  function (the resonance-integrator cache is not needed: the 720 min steps are fixed, the values are
  the same), the record is shared across threads. TLE parsing (Alpha-5, checksums — flag
  `MYB_TLE_NO_CHECKSUM` for edited sets), WGS-72/72 old/84, improved and AFSPC modes. Against
  `tcppver.out` (33 SGP4-VER sets, 700 rows): **5.0·10⁻⁹ km — half of the last printed digit**;
  3.5 years of deep-space SL-12 R/B (e = 0.79) — 1.15·10⁻⁷ km, exactly as python-sgp4 itself;
  the error codes [1, 1, 6, 6, 4, 3, 6] match. Found along the way: the reference's epoch is computed
  via a Julian date of ~2.45·10⁶ (rounding 40 µs), and for highly elongated orbits this gives up to
  4 mm at t = 0 — the epoch computation path is reproduced exactly.
- **Equinox-based frames** (`myb/frames.h`, kernel `data/frames/iau1980.tk`): IAU 1976
  precession, IAU 1980 obliquity and nutation (106 terms, Poisson series from data), MOD, TOD, TEME
  (TOD pole and the mean equinox projected onto the true equator — the SPICE definition; differs from
  R3(Δψ cos ε) by 10⁻⁹ rad), GMST 1982, TEME→PEF; all with derivatives.
  Against ERFA over 1800–2200: precession angles 10⁻¹⁷ rad, nutation 7·10⁻¹⁷ rad, matrices
  ≤ 5.5·10⁻¹⁶, GMST82 3.9·10⁻¹³ rad (both ERFA and we accumulate up to 1.7·10⁷ s in a single double);
  derivatives against differences 4·10⁻¹⁷ /s.
- **SPK type 10** (`myb_sgp4_spk10_state` in `ext/sgp4` on top of the new `myb_spk_segment_data` —
  raw segment words; the core does not evaluate type 10): generic segment metadata, the pair of
  element sets around t, SGP4 in AFSPC mode with the epoch via a UTC date in a single double, blending
  W = (1 + cos π(t − t1)/(t2 − t1))/2 with its derivative, rotation TEME → J2000. Against CSPICE N0067
  (`spkw10_c`, `spkgeo_c`, `evsgp4_c`): **2.7·10⁻¹² km** in near-Earth orbit and 1.9·10⁻⁹ km in
  deep space (there the epoch matters: SPICE's TDB−TT has a single term, ours is the
  Fairhead–Bretagnon series).
- **Elliptic integrals and Jacobi functions** (`myb/num.h`): R_F, R_D, R_J, R_C by Carlson
  duplication (5th-order series, r = 2⁻⁶⁴), complete and incomplete K, E, Π for any φ
  (quasi-periodicity), am, sn, cn, dn via AGM after reducing u modulo 2K. Found and eliminated
  four losses of significance: R_J at p ≪ x, y, z (92 ulp → the Carlson 1988 form with R_C(α², β²)),
  E(m) at m → 1 (48 ulp → 2R_G with the smallest argument last), 1 − m sin²φ (89 ulp →
  cos²φ + (1 − m) sin²φ), Π at n < 0 (550 ulp → the A&S 17.7.15 transformation to N ∈ (m, 1)).
  Against mpmath (60 digits, 1650 values): **≤ 5.6 ulp**, amplitude ≤ 7.4 ulp (accounting for |u|).
- **Attitude** (`myb/attitude.h`): quaternions (scalar first, Hamilton product,
  M(q) v = q̄ v q — the same “frame rotation” convention as `myb_rotation`), all 12 Euler-angle
  sequences with gimbal-lock handling, principal axes (cyclic Jacobi), Euler's equations
  for the body's angular momentum, torques: of a point mass on an inertia tensor, of a body's full
  field (−(m/M) r × a), quadrupole–quadrupole U₂₂ (via the derivative with respect to rotation;
  angular-momentum balance between the bodies and the orbit — 6.6·10⁻¹⁶).
  - **Exact torque-free top** (Landau–Lifshitz §37; Celledoni et al. 2008): body-frame angular
    momentum via Jacobi functions, precession about the constant angular momentum via Π(n; am u, k²)
    with n depending only on the moments of inertia; the axes are ordered so that the motion circulates
    about axis 3; symmetric bodies, stable rotations and spheres are special cases, the separatrix is
    `MYB_E_DOMAIN`. Against DOP853 with tolerance 10⁻¹⁴: **2·10⁻¹⁴ rad** (circulation about the axis
    of largest and of smallest moment, symmetric bodies in both orderings, stable rotation); symmetric
    top against the closed-form formula — 7·10⁻¹⁵ rad; the angular momentum in the inertial frame is
    constant to 4·10⁻¹⁶.
  - **Splitting** (Touma & Wisdom 1994): torque kick — exact free flow — kick, Yoshida
    weights; observed orders 2.00/3.99/5.98/7.95, without torque — the exact solution.
  - **Body with a fluid core** (`myb_core_mantle_rates`, JPL ephemeris model — Folkner et al. 2014
    §IV): mantle with a time-varying inertia tensor and its derivative, core with friction k_v and
    boundary oblateness.
- **Lunar librations against DE440** (`test_moon_libration`, all constants from `header.440`, initial
  conditions of the DE440 integration 1969-06-28): mantle with moments of inertia from β, γ, J̃₂,
  tidal deformation by Earth (k₂, time lag τ: Earth's position and the mantle frame at t − τ) and
  rotational deformation, core (IFAC, COBLAT, KVC), torques from Earth, Sun, Mercury, Venus,
  Mars, Jupiter on the lunar field up to degree 6 (degree 2 — from the distorted tensor), Earth's J₂ on
  the lunar quadrupole. **0.040″ over the first year, 0.071″ over two (0.6 m at the surface)**,
  angular velocity — 6.6·10⁻⁸ rad/d. Without tidal deformation — 0.64″, without its rate — 0.97″,
  without core friction — 0.11″, without U₂₂ — 0.085″, without rotational deformation — 0.077″; tide
  in the mantle frame at time t instead of t − τ — 29″: every part of the model works and is needed.
- **Regularization** (`myb/regular.h`): Kustaanheimo–Stiefel (u, u′, h, t; Stiefel & Scheifele
  1971) and Sundman (dt = r^α ds, α = 1, 3/2, 2) as right-hand sides for DOP853,
  `myb_regular_propagate` to a given physical time (event t(s) = t₁ on dense output). Two-body problem,
  20 revolutions, DOP853 10⁻¹²: e = 0.999 — **KS 1.2·10⁻¹² in 4 181 evaluations**, Cowell 1.3·10⁻⁸ in
  75 278, Sundman α = 1 — 1.8·10⁻⁸ in 40 289; the KS work does not depend on eccentricity. J₂,
  e = 0.7, 3 days: KS, Sundman and IAS15 agree to 1 mm at 23 000 km. A rectilinear orbit passes
  through r = 0 and comes back (x = 1 + 2.3·10⁻¹³).
- **Low thrust** (`myb/thrust.h`): constant-thrust arcs with linear mass flow (mass is a
  function of time, the integrator state remains (x, v)), steering along velocity, inertial,
  constant in RTN, and the Edelbaum yaw law with sign by cos u; Gauss equations in
  modified equinoctial elements (Walker et al. 1985); Edelbaum transfer (ΔV, β₀,
  state as a function of ΔV). Tsiolkovsky equation — 1.7·10⁻¹⁴; radial thrust 0.9 μ/(8r₀²) against
  the closed-form apocenter (μ − √(μ² − 8aμr₀²))/(4ar₀) — 3·10⁻¹³; 10 days of RTN thrust — Gauss
  against Cowell (IAS15) 0.44 mm at 8 617 km; transfer 7000 → 7200 km with a 1° plane change —
  a = 7200.001 km, i = 27.5025° (theory 27.5°; approximation of order thrust/gravity = 2.4·10⁻⁴).
- **Reproducibility**: phase 9 workload (equinox-based frames, elliptic functions,
  attitude, torque-free top, splitting, core–mantle, KS and Sundman, thrust, SGP4 and SPK type 10
  from `ext/sgp4`); golden `0x827497389931b8a9` is identical at -O0 and -O2; `MYB_REPRO_PREFIX=8`
  verifies the phase 8 golden `0xb1b77fd24243268f` — not a single earlier bit changed.
- **Windows** (Clang 23): `windows-dev`, `windows-release` — 61/61, `windows-ubsan` — 60/60 (no
  stack budget under the sanitizer); the hash is identical.
- **Gentoo** (2026-11-10, Clang 21): `linux-release` and `linux-asan` — 61/61 and 60/60, hash
  0x827497389931b8a9, all of the phase's numbers (SGP4, SPK type 10, librations, top, KS, ERFA) — the
  same to the last digit. Found: a check in `test_elliptic` compared against glibc `tanh`/`cosh`,
  which differ by an ulp — now against the correctly rounded `myb_tanh`/`myb_cosh`.
- **API changes**: `myb_eigenvector_real(n, a, λ, v, work)`, `myb_shoot(..., options, work,
  iterations, residual)` (workspace of n² and of `myb_shoot_work_size(n, m)` from the caller).
- **Not done in this phase**: SPICE CK reader (spacecraft attitude from files); coupled
  integration of the lunar orbit and rotation (positions are taken from DE440); derivative of the
  rotational deformation (estimated to be orders of magnitude smaller than the torques); optimal
  low-thrust control (Q-law, indirect methods, Sims–Flanagan); time element in KS; arbitrary Sundman
  exponent. Assigned to phases or closed: see "Old debts (phases 3–11) and the bindings".

### Phase 8 decisions

- **Small linear algebra** (`myb/num.h`): LU with pivoting, Cholesky, least squares by Householder
  reflections, eigenvalues of nonsymmetric matrices (Parlett–Reinsch balancing by powers of 2,
  reduction to Hessenberg form, implicit Francis double-shift QR following Golub & Van Loan,
  exceptional shifts at the 10th and 20th iterations, budget of 30 n), eigenvectors by inverse
  iteration; fixed order of operations. Check: roots of companion matrices, rotations, spectrum
  1500 … 1/1500 (as for the halo monodromy), trace and determinant of 50 random matrices, a Jordan
  block.
- **Object force model** (`myb/model.h`): a list of named terms a(t, x, v), the sum of the enabled
  ones in order of addition with per-component compensation; a term can be disabled, evaluated
  separately, differentiated — analytically (if the term provides a Jacobian;
  `myb_point_mass_jacobian`) or by central differences with step ∛ε·scale; parameters are
  addresses of doubles inside the term contexts (Cr, A1, empirical terms) with their own
  difference step. This is the common basis of the budget, covariances, estimation and targeting.
- **Uncertainty** (`myb/uncert.h`): variational equations — the columns of the state transition
  matrix and of the parameter sensitivities are integrated together with the state as a single
  IAS15 system (3 (7 + n_p) coordinates, dense output for everything). Against the closed-form
  two-body state transition matrix (`myb_kepler_stm`) over 4 revolutions: 5.5·10⁻¹⁵ with the
  analytic Jacobian, 2.3·10⁻⁹ with the finite-difference one; the sensitivity to the pressure
  parameter and a Φ column with drag against differences of complete propagations — 3·10⁻⁸ and
  9·10⁻⁹. Covariances (M P Mᵀ, RTN axes), the xoshiro256** RNG (the authors' reference values)
  with Box–Muller normal deviates on correctly rounded log, sin, cos — the same sequence
  everywhere; linear covariance against 4000 Monte Carlo samples — within 3% (sampling error
  1.1%).
- **Restricted three-body problem** (`myb/cr3bp.h`): equations, Jacobi constant, analytic
  Jacobian, Lagrange points (bisection of a monotonic function down to double resolution),
  correction of symmetric periodic orbits by single shooting over a half-period (Howell 1984;
  symmetry about the xz plane or the x axis, choice of the fixed variable and of the index of the
  y = 0 crossing), monodromy, stability index, manifolds from real eigenvectors, Richardson (1980)
  third-order halo, pulsating rotating frame of two real bodies. Check — the JPL periodic orbit
  catalog (56 orbits of eight families: L1/L2 halo, L1/L2 Lyapunov, vertical, axial, Earth–Moon
  DRO, Sun–Earth L2 halo): Jacobi constant 3.6·10⁻¹⁵, the catalog state is a fixed point of the
  corrector to 9·10⁻¹², period 6·10⁻¹², stability index to 10⁻⁵, except for two L2 Lyapunov
  orbits passing 500 km above the Moon (1.1·10⁻⁴; our index changes by ~10⁻⁵ when the integration
  is tightened); catalog Lagrange points to 10⁻¹⁴; Richardson's first approximation — within 1–2%
  of γ of the converged orbit, period to 0.15%; the unstable manifold of an L1 Lyapunov orbit
  diverges over one period by exactly a factor λ = 2446, the stable one backward in time —
  likewise.
- **Transfers** (`myb/transfer.h`): B-plane (Kizner; Sergeyevsky et al.; Vallado) — B·T, B·R, θ,
  v∞, periapsis, time to periapsis; sphere of influence, Hill radius, turn angle, departure
  impulse; Newton shooting with a minimum-norm step (m < n), least squares (m > n), step halving.
  Check: |B| and the asymptote geometry from the distant past of a hyperbola; Lambert's problem by
  shooting with the closed-form state transition matrix — agrees with `myb_lambert` to
  4·10⁻¹² m/s; B-plane targeting of an Earth flyby with lunar perturbation using the state
  transition matrix from the variational equations — 3 iterations to a miss of 2·10⁻⁸ m.
- **Effects budget** (`myb/budget.h`): each term is toggled relative to a reference propagation;
  the effect is the largest position difference over the interval (the full nonlinear effect),
  plus the largest acceleration of the term; flags “omitted” (disabled, but the effect exceeds the
  tolerance) and “negligible”. Check: on a 400 km orbit over one day J2 — 1100 km, disabled drag
  is flagged (28 km), a 10⁻¹⁴ m/s² term is negligible (4·10⁻⁶ m), the effect of a small force is
  linear (×2.0000).
- **Real spacecraft against Horizons** (`test_spacecraft`, phase criterion; estimation of the
  initial state and parameters by Gauss–Newton with partial derivatives from the variational
  equations, forces from DE440):
  - **Tesla Roadster** 2018-03 … 2025-03 (7 years, a 0.05 au close approach to Mars in 2020): Sun
    with the Schwarzschild term, planets and Moon, radial A1 (1 au/r)² of JPL solution #11 —
    **RMS 306 m**, maximum 1.3 km at distances of 2·10¹¹ m; A1 = 5.95·10⁻⁸ m/s²;
  - **OSIRIS-REx** through the Earth gravity assist of 2017-09-22: Earth with J2–J4 about the
    CIP, Sun, Moon, planets, pressure with Earth penumbra — **RMS 12.6 m** over 6 days; perigee
    (17 214 km above the Earth) — 54 m in radius and 0.17 s in time from Horizons; B·T, B·R one
    day before perigee — 32 and 30 m;
  - **JWST** 2023-01 … 03 in the pulsating Sun–(Earth–Moon barycenter) frame: x from 1.0079 to
    1.0111 (L2 — 1.0101), |y| < 0.005, |z| < 0.003, Jacobi constant 3.000779 — slightly below
    C(L2) = 3.000894, as for libration orbits; the week 2023-02-01…08 in the ephemeris model with
    pressure of constant direction relative to the Sun line — 260 m RMS (pressure
    2.0·10⁻⁷ m/s²).
- **Reproducibility**: phase 8 workload (LU, eigenvalues, RNG, variational equations with a
  parameter, budget, CR3BP from Richardson to monodromy, B-plane, shooting); golden
  `0xb1b77fd24243268f` is identical at -O0, -O2 and under UBSan; `MYB_REPRO_PREFIX=7` checks the
  previous phase 7 golden.
- **Gentoo** (2026-11-10, Clang 21): `linux-release` and `linux-asan` — 53/53, hash
  0xb1b77fd24243268f, the Roadster, OSIRIS-REx and JWST numbers are the same to the last digit.
  Found: `hypot` from libm in the stability index (test `test_no_libm_imports`) — replaced with
  the correctly rounded `myb_hypot`; the golden bits did not change.
- **Not done in this phase**: multiple shooting and quasi-periodic (Lissajous) orbits; family
  continuation with pseudo-arclength; observables for estimation (range, Doppler, angles) —
  phase 10; a budget from the linear estimate ∫Φ a_k dt (currently — full propagations). Assigned to phases or closed: see "Old debts (phases 3–11) and the bindings".

### Phase 7 decisions

- **Spherical harmonics** (`myb/grav.h`): potential and acceleration by the Holmes–Featherstone
  (2002) column recursion with seeds P̄ₘₘ scaled by 2⁹³⁰ (degree up to 2700 at any latitude); the
  gradient is written without division by cos φ (the P̄/u columns are recursed from their own
  seeds) — the pole is an ordinary point. Coefficients are data: an ICGEM reader (formats 1.0 and
  2.0: `gfct`, `trnd`, `acos`, `asin`, epochs `yyyymmdd[.hhmm]`), truncation by degree, the
  file's tide system; `set_time` recomputes the time-variable terms, `add_correction` adds
  corrections (tides) until the next `set_time`; `clone` — one field per thread. Check: MPFR
  256 bits (potential by a high-precision recursion, acceleration by differences of that
  potential, independently of the derivative formulas): degree 120 at the pole, next to it, at
  the equator — 4.7·10⁻¹⁶ / 7.5·10⁻¹⁶; degree 2190 at latitude 62° (unscaled seeds would
  underflow) — 4.8·10⁻¹⁶ / 2.5·10⁻¹⁶. Point mass, J2 and the split by degree — in closed form.
- **Polyhedron** (Werner & Scheeres 1997): exact field of a homogeneous closed polyhedron inside
  and outside; Laplacian via solid angles; on the surface (edges, vertices, face diagonals) the
  potential and acceleration take their limiting values — no error. Check — the closed-form
  formula for a rectangular prism (Nagy et al. 2000) in 40 digits: 3·10⁻¹⁶ / 7.6·10⁻¹⁶ ×
  max(1, (r/size)²) (the terms grow as r², the field falls off as 1/r — this is written into the
  contract). The mesh is checked for closedness and orientation.
- **Mascons**: a sum of point masses with compensation in array order; against their own
  harmonic series to degree 40 (independent Legendre recursion) — 1.2·10⁻¹⁵.
- **Figure–figure**: the l = l′ = 2 term of the mutual potential in STF form (Damour, Soffel, Xu
  1991): U₂₂ = G/(36 r⁵)[105 (n·Q_A·n)(n·Q_B·n) − 60 n·Q_A Q_B·n + 6 tr(Q_A Q_B)], forces on both
  centers of mass; the “field of one at the center of the other” terms — via harmonics of any
  degree. Check: two clouds of point masses without odd moments, exact double sum: U₂₂ reduces
  the mutual-force error by a factor of 29 at 12 sizes and by a factor of 95 at 24; the remainder
  falls off as r⁻⁸. Q from the normalized C₂ₘ, S₂ₘ and its rotation are verified against degree 2
  of the field.
- **Ellipsoidal harmonics** (Garmier & Barriot 2001): Lamé functions of the four classes K, L,
  M, N from three-term recurrences (derived by substitution into the Lamé equation; the
  eigenvalues are real and distinct — roots of a continuant, then the roots of the polynomial Q in
  (0, k²) via Stieltjes). The exterior term F(λ)E(μ)E(ν) = R(λ)·H(x, y, z), where H is Niven's
  Cartesian harmonic polynomial (a product of the quadrics (t−h²)(t−k²)(x²−t) + y²t(t−k²) +
  z²t(t−h²)), and R(λ) is an integral in a form with a bounded integrand: a gradient free of the
  singular sets of ellipsoidal coordinates. Check: harmonicity of each basis term up to degree 16;
  the Niven form against E(μ)E(ν) in coordinates from the cubic — 3·10⁻¹²; a 4×2.4×1.6 block
  (exact field — polyhedron) inside the Brillouin sphere: the ellipsoidal series 8.9·10⁻³ →
  6.5·10⁻⁵ from degree 4 to 16, the spherical series of the same block (exact moments) diverges
  3.8·10⁻² → 0.26.
- **Tides** (IERS 2010 ch. 6, all numbers are data in `data/grav/iers2010-tides.tk`): step 1 with
  complex Love numbers and k⁽⁺⁾ for degree 4; step 2 — the 71 rows of Tables 6.5a–c as Poisson
  series in χ = GMST + π and the Delaunay arguments (the conversion script cross-checked, in every
  row, the Doodson number against the multipliers and the Doodson multipliers against the Delaunay
  multipliers); step 3 — the permanent tide for zero-tide fields (selected by the file's
  `tide_system`, no guessing for mean-tide); pole tides, solid and ocean, relative to the secular
  pole (ch. 7, eq. 21); ocean tides from files of the `fes2004_Cnm-Snm.dat` kind (eq. 6.15).
  Check: step 1 against the addition theorem (k GMⱼR²ⁿ⁺¹/(rⁿ⁺¹rⱼⁿ⁺¹) Pₙ(cos ψ)) — 2·10⁻¹⁶; step 2
  against the K1 example from the IERS text with ERFA GMST — 1.3·10⁻¹⁶; ocean against eq. 6.15
  with ERFA arguments — 3·10⁻¹⁶ (the remainder is GMST 1982 versus 2006, 0.056″).
- **Forces** (`myb/force.h`, pure functions, accelerations are added): third body in Battin's
  form (against 50 digits, 4·10⁻¹⁶ for |r|/|s| down to 10⁻¹¹); radiation pressure on a sphere and
  on plates (specular, diffuse, absorption with thermal re-emission; a sphere of 320 thousand
  plates gives Cr = 1 + 4ρ/9, and 1 for a mirror); the fraction of light in eclipse — exact
  geometry on the celestial sphere (the disk of the occulting body is a spherical cap) and
  quadratic limb darkening, an integral over annuli with exact arcs (against the closed-form
  formula for the intersection of caps and mpmath quadrature — 1.1·10⁻¹⁴: penumbra at 420 km, at
  LAGEOS, from the distance of the Moon, annular eclipse); planetary albedo and IR per
  Knocke–Ries–Tapley (uniform Lambertian emitter — to 0.4/rings²); drag; empirical RTN.
- **Atmosphere** (`ext/nrlmsise00`, separate target `myb::nrlmsise00`): NRLMSISE-00 — a port of
  Brodowski's C version (public domain; the model and coefficients are U.S. Government material,
  marked per 17 U.S.C. 403 in `NOTICE`): the shared static variables became state passed to the
  functions (thread safety, no global mutable state), tables are `static const`, exp/log/sin/cos
  are myb's correctly rounded functions (reproducibility), ghp7 with its unbounded iteration is
  excluded. 17 reference cases (0–1000 km) — to the printed precision of 7 digits (4.8·10⁻⁷).
- **SP3** (`myb/kin.h`): SP3-c/d reader (positions, velocities, time system, coordinate
  system).
- **LAGEOS-1 against ILRS** (`test_lageos`; phase criterion): the week 2025-06-01…07, ILRS orbit
  in SLRF2020; a model made of myb parts — GOCO06s to degree 30 with trends and annual terms,
  IERS 2010 tides (steps 1–3), FES2004 to degree 20, pole tides, Sun/Moon/planets from DE440, GR
  (eq. 10.12), solar pressure with Earth and Moon penumbra, albedo and IR, empirical terms as used
  by the ILRS centers (constant along-track and 1/rev; Cr = 1.13 fixed, as they do: over a week a
  force of constant direction is indistinguishable from the 1/rev terms). After fitting the
  initial state: **RMS 2.0 cm** (R 2.6 mm, T 1.8 cm, N 8 mm); constant along-track −1.4 pm/s²
  (LAGEOS thermal drag). Contributions (RMS without a part of the model): solid tides 21 cm,
  ocean 6.3 cm, GR 2.06 cm (radially 5.3 mm versus 2.6), shadows 2.1 cm.
- **DE440 with figures** (`test_de440_rel`, third model): EIH + Earth J2–J5 about the CIP, lunar
  harmonics to degree 4 in the DE440 principal-axis frame (binary PCK), solar J2, Earth–Moon
  U₂₂, all from `header.440`. After one year: **Moon 0.65 m** (without figures 23 km), Mercury
  14 m (without solar J2 106 m), Earth 124 m (270); after 10 years the Moon 239 m — exactly the
  DE440 secular tidal acceleration (25.8″/cy² gives 240 m in 10 years); the model has no tides.
- **Performance**: the TT→TDB conversion performed four evaluations of the 787-term series; two
  are enough (after the first the error is < 2·10⁻¹² s, after the second < 2·10⁻²¹ s at a
  resolution of 10⁻¹⁶ s; golden bits unchanged); albedo — azimuth phasors by rotation instead of
  sin/cos for every element; ocean tides write their corrections directly into the coefficients.
- **Reproducibility**: phase 7 workload (harmonics, time-variable terms, tides, polyhedron,
  mascons, figure–figure, ellipsoidal harmonics, forces, SP3); golden `0xbd6bc51c0c5d4fda` is
  identical at -O0, -O2 and under UBSan. `MYB_REPRO_PREFIX=5|6` computes the hash of only the
  phases up to the given one and checks it against that phase's golden (`0x38f060f6206c4bf6`,
  `0x9cb768b2ccddf18b`): adding phase 7 did not change a single bit of the earlier ones. 48 tests.
- **Not done in this phase**: admittance (secondary waves of Table 6.7) and atmospheric S1/S2 for
  ocean tides; Earth and Moon tides with time lags in the N-body problem (DE440: k₂ and τ for
  three bands); rotational dynamics (figure–figure torques) — phase 9; JB2008/DTM2020 in `ext/`;
  an accelerated TDB−TT series (787 sines per call; the LAGEOS test takes it from an hourly
  table). Assigned to phases or closed: see "Old debts (phases 3–11) and the bindings".

### Phase 6 decisions

- **EIH in the N-body problem** (`MYB_INTERACTION_EIH`, `myb_nbody_set_relativity`, `myb/dyn.h`):
  the Einstein–Infeld–Hoffmann equations with the PPN parameters β, γ in Moyer's form (DESCANSO 2,
  eq. 4-61) — the same ones the JPL ephemerides integrate. The interaction graph stores exclusions
  (NONE) and explicit NEWTON/EIH levels; by default, after `set_relativity` all pairs are EIH. The
  potentials U_i and the Newtonian accelerations a_j in the 1PN terms are taken according to the
  graph. The Newtonian pass is unchanged (symmetric, vectorized, parallel — the same bits: the old
  reproducibility golden without the new workloads did not shift by a single bit); after it come
  a potentials pass and a pass of the 1PN terms over receivers in increasing source order with its
  own Sum2. The accelerations depend on velocities: `myb_nbody_acceleration` without velocities
  and the Wisdom–Holman integrators reject such a system (`MYB_E_DOMAIN`).
- **Independent EIH check**: the Euler–Lagrange equations of the global PPN Lagrangian
  (Bernus et al. 2023, arXiv:2301.01186, App. B), derived symbolically (sympy) and solved to first
  order in 1/c² in 50-digit arithmetic (mpmath): 48 cases, 2–5 bodies, (β, γ) = (1, 1),
  (1.3, 0.7), (0.8, 1.25), c = 10 and 1000 — error 2.7·10⁻¹⁰ of the magnitude of the 1PN terms
  (they are ≥ 1.6·10⁻⁶ of the Newtonian ones: i.e. double rounding), all paths and pools — the
  same bits.
- **1PN integrals** (`myb_nbody_pn_invariants_of`): energy, momentum and the relativistic
  barycenter (the BCRS origin of the JPL ephemerides) from the same Lagrangian. Triple system,
  v/c ~ 10⁻², 60 time units (IAS15): the 1PN energy changes by 6·10⁻⁹ while the Newtonian one
  changes by 1.6·10⁻⁵; with c doubled — 12 times smaller (the expected order c⁻⁴ plus integration
  noise).
- **Mercury's perihelion** (Sun and Mercury, DE440 GM, 60 revolutions, perihelia found by the
  event r·v = 0): 42.9804″ per century against 6πGM/(c²a(1−e²)) — discrepancy 5·10⁻⁹; at
  (β, γ) = (1.5, 0.7) and (0.6, 1.3) the factor (2 + 2γ − β)/3 is reproduced to 1.3·10⁻⁷
  (corrections for Mercury's finite mass, mass ratio 1.7·10⁻⁷).
- **Satellite in the local frame** (`myb/rel.h`, IERS 2010 eq. 10.12): the Schwarzschild,
  Lense–Thirring (spin of the central body) and de Sitter (external body) terms — as separate
  functions. LAGEOS orbit, 600 revolutions with IAS15: node precession from Lense–Thirring
  30.631 mas/yr and from de Sitter 19.188 mas/yr — agree with the analytic
  2GMj/(c²a³(1−e²)^{3/2}) and (3/2)GM☉n/(c²a☉) to rounding; the perigee of an eccentric orbit
  (e = 0.3) over 200 revolutions — in agreement with the formula and the factor (2 + 2γ − β)/3.
- **Local frame** (`myb_rel_global_to_local`, IAU 2000 Resolution B1.3, Soffel et al. 2003):
  X = r + c⁻²[½v_E(v_E·r) + w_ext r + r(a_E·r) − ½a_E r²] and the inverse to O(c⁻⁴); in the limit
  without field and acceleration it coincides with the Lorentz transformation (to 10⁻⁹ m at
  6800 km); the correction on a 6800 km orbit is 66 mm. Time offset −v_E·r/c²; TCB − TCG at the
  center — the time module.
  **Proper time** (eqs. 10.6–10.9): GPS clocks run faster than TT by 4.4647·10⁻¹⁰ — exactly the
  amount compensated by the GPS frequency offset (IERS 2010, §10.2).
- **DE440** (`test_de440_rel`, requires `data/download`): the Sun, planets (system barycenters),
  Earth and Moon as point masses with GM and c from `header.440`, IAS15 from DE440 states at
  J2000:

  | after 1 year, m | EIH | Newton |
  |---|---|---|
  | Mercury (heliocentric) | 106 | 57 652 |
  | Venus | 10 | 98 803 |
  | Mars | 68 | 39 723 |
  | Jupiter / Saturn (barycentric) | 93 / 22 | 633 / 93 |
  | Moon (geocentric) | 22 999 | 19 426 |

  After 10 years Mercury — 0.95 km (Newton — 1823 km), Venus — 0.17 km. The remainder is explained
  by what is absent from the model: the figures of the Earth and Moon and tides (the Moon, and
  through it the Earth — 1/82 of the Moon's error), the ~370 small bodies of DE440 (their momentum
  shifts the barycenter by ~270 m/yr — this is visible in the heliocentric outer planets; Mars).
  These are phase 7 tasks.
- **Reproducibility**: EIH accelerations (1 and 4 threads — bitwise identical), 1PN integrals,
  the 10.12 terms and the BCRS→GCRS transformation have been added to the golden; golden
  0x9cb768b2ccddf18b is identical at -O0/-O2 and under UBSan; 41 tests.

### Gentoo check (phases 6–7)

2026-11-10, the same VM (Clang 21.1.8, glibc, CMake 4.3, Ninja, 8 logical processors); the tree
is the working copy (phases 6 and 7 not yet committed), data downloaded with
`cmake -P cmake/FetchData.cmake` (the new GOCO06s and FES2004 — SHA-256 matched). All seven
presets — **48/48**, the build has not a single warning under `-Werror`:

| preset | test time | notes |
|---|---|---|
| `linux-release` | 63 s | hash 0xbd6bc51c0c5d4fda — the same as on Windows; LAGEOS (2.01 cm) and DE440 with figures — the same numbers |
| `linux-dev` | 190 s | |
| `linux-asan` | 208 s | ASan + UBSan + LeakSanitizer |
| `linux-asan-debug` | 978 s | `-O0`, pointer-pair checking |
| `linux-ubsan` | 112 s | |
| `linux-lsan` | 66 s | |
| `linux-tsan` | 326 s | then 20 repetitions of `test_pool`, `test_nbody`, `test_whfast`, `test_repro`, `test_rel` (100 runs) and 2 repetitions of `test_lageos` (12 propagations in the thread pool) — not a single report |

Found: LeakSanitizer — a leak in the `test_grav` test itself (the point-mass field was not
freed); the library is clean. On Windows the same 48 tests pass in `windows-dev`,
`windows-release` and `windows-ubsan`.

### Gentoo check (phases 0–5)

2026-11-09, Gentoo Linux (VMware VM, kernel 6.18, i5-9300H), Clang 21.1.8, glibc, ld.bfd,
CMake 4.3, Ninja; the tree is an archive of the phase 5 commit, data downloaded with
`cmake -P cmake/FetchData.cmake` (SHA-256 matched). Presets `linux-dev`, `linux-release`,
`linux-asan` (ASan + UBSan + LeakSanitizer, `-O2`), `linux-asan-debug` (the same at `-O0` with
pointer-pair checking), `linux-ubsan`, `linux-lsan`, `linux-tsan`: 39/39 everywhere,
`test_no_libm_imports` passes, reproducibility hash 0x38f060f6206c4bf6 is the same as on Windows
(Clang 23, mingw-w64) — bit for bit with a different compiler, C library and linker. Positive
control on the same VM: LSan, ASan and UBSan catch a planted leak, an out-of-bounds array access
and an overflow.

**ThreadSanitizer** (`linux-tsan`, `-O2 -fsanitize=thread`, 153 translation units,
`TSAN_OPTIONS=halt_on_error=1:abort_on_error=1:second_deadlock_stack=1:report_signal_unsafe=1:history_size=7`),
2026-11-09, the same VM (8 logical processors): all 39 tests — 39/39 without a single report
(57 s); then the multithreaded tests — `test_pool` (pools of 1–13 threads, arbitrary range
partitions), `test_nbody` (accelerations on 1–13 threads, all paths, vector blocks, integration
with a pool), `test_whfast` (WH with and without a pool), `test_repro` (N bodies on 1 and 4
threads) — were repeated 30 times each (`ctest --repeat until-fail:30`, 120 runs, 208 s): not a
single data race, not a single lock-order inversion. The positive control (two threads
incrementing a shared counter without synchronization) is caught by TSan (“data race”, code 66).
What this confirms: the shared state of the thread team (`generation`, `pending`, `busy`, `job`)
is changed only under the mutex, worker threads write disjoint parts of the result, and
everything they have written is visible to the calling thread after `myb_platform_team_run` (a
happens-before relation through the mutex and the condition variable). Not verified: the Win32
implementation (SRWLOCK) — TSan does not exist for mingw-w64; it is built the same way and gives
the same bits (`test_repro`).

Found and fixed (did not manifest on Windows):
- `platform_posix.c` compiled for the first time: in strict C23 glibc does not declare
  `O_CLOEXEC` — `_POSIX_C_SOURCE 200809L` added;
- `fegetround`/`feraiseexcept` (CORE-MATH) live in libm in glibc — `myb` links against `m` on
  non-Windows (the ban on inexact libm functions is checked by `test_no_libm_imports`);
- `-Wpsabi`: in the SysV ABI, Clang 32-byte vectors passed by value without AVX change the calling
  convention — the vector blocks of `dyn/nbody.c` are now passed by pointer (bits unchanged);
- ASan pointer-pair checking at `-O2` (Clang 21) fires on comparisons created by the optimizer
  (pointer against −9): at `-O0` the same 39 tests are clean, so `pointer-compare` and
  `pointer-subtract` are enabled only in the ASan Debug build.

### Phase 5 decisions

- **Files are memory-mapped** (`myb_platform_map_file`: `CreateFileMappingW`/`MapViewOfFile`,
  POSIX `mmap`): DE441 (3 GB) is read page by page with no system call per query. The DAF reader
  decodes every number from bytes (`LTL-IEEE` and `BIG-IEEE` on any host), checks the FTP string
  (corruption from a text-mode transfer) and walks the summary record list with protection
  against cycles.
- **SPK** (`myb_spk`, `myb/kin.h`): types 1, 2, 3, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 21;
  the layout of every segment is validated at load time (lengths, monotonic epochs); queries
  change nothing (thread-safe). Priorities follow the SPICE rule (later-loaded file, later-written
  segment). Body chain: from the target and the center up to the common node, only the required
  segments, Sum2 addition. The argument is `myb_time` (two parts): t − t_mid is exact, so the
  normalized Chebyshev argument does not lose digits as it would with a single double.
- **SPK check.** For each type, files were written with CSPICE N0067 (through spiceypy, as a
  reference only; no CSPICE code is used) by `spkw*_c`/`spkw19_`/`spkw01_`, and 7 thousand states were computed at
  random epochs, nodes and midpoints (window selection at boundaries). Agreement: positions
  ≤ 6·10⁻¹⁶, velocities ≤ 10⁻¹⁴. Hermite velocities (types 12, 13, 19) differ by up to 1.1·10⁻¹⁴ —
  this is CSPICE's error: against the exact interpolant (mpmath, 300 bits) ours is ≤ 1.6·10⁻¹⁶,
  CSPICE's is up to 2.3·10⁻¹⁴ (its interpolation works on absolute epochs). Found along the way:
  the SPICE writer writes the type 19 interval selection flag as −1, not 0 as in the
  documentation; the length of a type 21 difference line in the segment trailer is MAXDIM, a
  record is 4·MAXDIM + 11 words (established on a Horizons file); the type 17 velocity is the full
  time derivative of position (numerical check), derived analytically.
- **testpo** (`test_testpo`, requires the downloaded `de440.bsp`, `testpo.440`,
  `moon_pa_de440_200625.bpc`): 11 354 positions and velocities — 7.1·10⁻¹⁵ au (≈1 mm) and
  2·10⁻¹⁷ au/d against the JPL tolerance of 10⁻¹³; 930 lunar libration values — 1 ulp (the angle ψ
  accumulates up to 3.5·10⁴ rad); nutations (14) are absent from the SPK version of DE440.
- **Horizons, body**: a type 21 SPK for 433 Eros from the Horizons API against Horizons vectors —
  1.8·10⁻¹² (0.3 m). CSPICE on the same file gives the same difference: Horizons prints its own
  integration, not the exported difference lines.
- **CIP/CIO IAU 2006/2000A** (`myb_cip`): all coefficients are data
  (`data/frames/iau2006-cip.tk`, IERS 2010 tables 5.2a, 5.2b, 5.2d with the SHA-256 of the
  sources: X 1600, Y 1275, s + XY/2 66 terms, 14 fundamental arguments). A generic Poisson series
  evaluator (`num/poisson.c`): the unique argument vectors (1311 for 2941 terms) are collected at
  load time; for a given date the phasors e^{imF_k} are computed (206 correctly rounded sin/cos),
  and the phasor of an argument is the product of the phasors of its nonzero multipliers. 47 µs
  for X, Y, s with derivatives (ERFA without derivatives — 76 µs). The matrix C = R3(−s) Bᵀ in the
  form of eq. 5.10 with an analytic derivative; W(x_p, y_p, s′), R3(ERA) with the rate
  ERA·(1 − LOD/86400).
- **Against ERFA** (`test_erfa_frames`, 3000 dates 1800–2200): X, Y — 1.7·10⁻¹⁷ rad, s —
  1.5·10⁻¹⁹, s′ — 2·10⁻²⁵; ERA and the matrix — 3.6·10⁻¹³: this is ERFA's rounding (eraEra00 adds
  ~200 revolutions in a single double). Our ERA keeps the whole days separate from the fraction,
  with the product in double-double: error ≤ 1.3·10⁻¹⁵ rad against MPFR (`test_mpfr_frames`).
  Found during checks: the first version rounded Tu = days + fraction before the multiplication
  and lost 10⁻¹³ rad.
- **EOP** (`myb_eop`): EOP 20 C04 and finals2000A (Bulletin A columns; empty LOD, dX, dY in the
  predictions — “no data” = 0 in interpolation, missing x, y, UT1 is an error). Lagrange
  interpolation over 4 rows, as in the IERS interp.f, but in the user's time scale (continuous),
  with UT1 taken as UT1 minus that scale: a leap second inside the window causes no jump (test:
  UT1 steps across 2016-12-31T23:59:60 — 2·10⁻¹¹ s). Sub-daily corrections
  (`data/frames/eop-tides.tk`: 71 ocean terms of tables 8.2–8.3, 10 + 11 libration terms of
  tables 5.1a/b, from interp.f, PMSDNUT2.F, UTLIBR.F) use the same series evaluator; LOD comes
  from the derivative of the UT1 series. The argument χ = GMST + π is taken from UT1 (IERS takes
  it from MJD UTC, which repeats during a leap second). Against an independent Python
  implementation from the original Fortran tables: pole 2·10⁻¹⁴″, UT1 1.7·10⁻¹⁴ s, LOD
  8·10⁻¹¹ s (the finite-difference limit of the reference).
- **WGCCRE** (`myb_wgccre`): RA, Dec, W with precession–nutation terms, phase angles of degree
  `BODYb_MAX_PHASE_DEGREE` (Mars in pck00011 — 2), its own epoch `JED_EPOCH`, its own base frame.
  W = w₀ + w₁d in double-double with reduction modulo 360°. Against CSPICE (18 bodies):
  ≤ 1.1·10⁻¹⁰, and the discrepancy grows with the rotation rate (Venus 10⁻¹³, Earth 2.5·10⁻¹¹,
  Vesta 1.1·10⁻¹⁰) — this is CSPICE's W in a single double; our W is accurate to 10⁻¹⁵ rad (MPFR).
  Found: the derivatives of the nutation terms were initially divided by a degree one extra time.
- **NAIF TK frames** (`myb_tkframe`: ANGLES, MATRIX, QUATERNION) and
  `data/frames/naif-inertial.tk` — the SPICE built-in inertial frames as data (B1950, FK4,
  GALACTIC, ECLIPJ2000, ECLIPB1950): they agree with CSPICE to 1.1·10⁻¹⁶. This was needed
  because the ITRF93 binary PCK is defined relative to ECLIPJ2000 (code 17), not J2000.
- **Geodesy** — the closed form of Vermeille (2011), including the interior of the ellipsoid:
  angles 3.3·10⁻¹⁶ rad, heights 5.6·10⁻¹⁶ of a + |h| at heights from −6000 km to 4·10⁵ km. ERFA
  (eraGc2gde, Fukushima 2006) at heights of 1.3–1.9·10⁷ m loses accuracy up to 2.6·10⁻¹¹ rad, so
  the reference is the original coordinates, and the comparison with ERFA is made only near the
  surface (6·10⁻¹⁶).
- **Horizons, sites on Earth** (`test_earth_frames`): the Moon from the Earth's center minus the
  Moon from two sites (Pulkovo and the southern hemisphere), 10 dates in 2020 — the site vector in
  ICRF. The myb chain (WGS-84 → ITRS → EOP C04 + tides → IAU 2006/2000A → GCRS, TT and UT1 from
  the time scales): 0.19 m and 1.2·10⁻⁵ m/s; there is no radial component, the difference is a
  rotation by 3–11 mas. Via the NAIF ITRF93 binary PCK (`earth_latest_high_prec.bpc` +
  ECLIPJ2000) — 0.077 m: Horizons orients the Earth with the ITRF93 model (IAU 1976/1980 with
  nutation corrections from EOP, J2000 = ICRF), and 0.19 m is the difference between that model
  and IERS 2010, not an error of the chain.
- **Orbital frames** RTN, VNC, LVLH (CCSDS) with derivatives from the acceleration:
  orthonormality 4·10⁻¹⁶, derivatives against central differences along a conic section
  5·10⁻¹¹ s⁻¹.
- **Data**: the large files go via `scripts/fetch-data.cmd` (`cmake -P cmake/FetchData.cmake`) into
  `data/download/`: fixed releases with SHA-256 (DE440, testpo, header, Moon PA, pck00011),
  live series (EOP C04, finals2000A, ITRF93 PCK) — on `refresh`. The repository holds small
  test samples (`tests/data/eop`, `tests/data/spk`, `tests/data/frames`, pck00011.tpc).
- **Reproducibility**: all SPK types, CIP, ERA, EOP with tides, WGCCRE, geodesy and orbital
  frames have been added to the golden hash; the golden hash without the new workload is the
  previous 0xde3d79ea60027aae (not a single old bit changed); the new 0x38f060f6206c4bf6 is
  identical at -O0/-O2, x86-64-v2/v3 and under UBSan.
- **Not done in this phase**: a grid and interpolation of the nutation series (plan §6; currently
  47 µs per date), SPK chains through different reference frames (currently an `MYB_E_FORMAT`
  error; the rotation is done manually via `myb_tkframe`), SPK types 10 (TLE — phase 9) and 15,
  equinox-based frames (TOD, TEME), local GR frames (GCRS↔BCRS — phase 6). Assigned to phases or closed: see "Old debts (phases 3–11) and the bindings".

### Phase 4 decisions

- **N-body system** (`myb_nbody`, `myb/dyn.h`): bodies are specified by GM only; massive (they
  attract) and test bodies (they feel the massive ones but do not attract themselves). The
  acceleration of a receiver is a direct sum over all sources in increasing order of j with Sum2
  addition (TwoSum on every term, errors kept separately; Ogita, Rump, Oishi 2005). Verified
  against a double-double sum: error ≤ 3.5·10⁻¹⁶ of Σ|terms| (plain summation of the same terms —
  2.1·10⁻¹⁵), 400 + 50 bodies with GM from 10⁵ to 10²⁰.
- **Independence from threads and paths.** The symmetric path (each pair once, the term for j is
  −(GM_i q dx), bitwise equal to GM_i q (x_i − x_j)) and the “per-receiver” path give the same
  bits: the order of the terms for each receiver is the same. The per-receiver path is
  parallelized over i via `myb_pool` (platform layer: Win32 SRWLOCK/CONDITION_VARIABLE, POSIX
  pthread; threads are created once and wait for work) and processes receivers in blocks of 4 in
  Clang vectors (no intrinsics: SSE pairs on x86-64-v2, AVX on v3, NEON on ARM) — each lane
  performs exactly the scalar operations; the self-pair is suppressed by a mask, and adding an
  exact +0 to Sum2 does not change a single bit (neither s nor c is ever −0). The test compares
  all paths, all tails (1–17 massive + 0–6 test), pools of 1–13 threads and arbitrary part
  boundaries — bit for bit. Measurement (2000 bodies, one thread, per interacting pair):
  symmetric 15.8 ns, vector per-receiver 12.2 ns on v2 and 6.2 ns on v3 (even though each pair is
  computed twice), scalar per-receiver 17 ns; 733 bodies on 8 threads — 0.60–0.90 ms versus
  3.8–4 ms for the symmetric path.
- **Interaction graph**: a level for each directed pair (currently `NONE` and `NEWTON`; levels for
  harmonics, 1PN and tides will be added by later phases); exclusions are stored sorted, and the
  receiver sums the gaps between excluded sources without branches in the inner loop. Found
  during checks: the bound of the loop over the gaps was computed from a changing counter and cut
  off the tail of the sources for a receiver with many exclusions — fixed, and a check “all
  sources visited, all exclusions accounted for” (`MYB_CHECK`) was added.
- **Invariants** (energy, momentum, angular momentum, center of mass) — in double-double; energy
  as G·E (bodies are specified by GM): integration errors of 10⁻¹⁶ are resolved.
- **Classical solutions** (test `test_nbody`): the Lagrange solution (GM 1, 10⁻³, 10⁻⁶, stable
  by the Routh criterion) — 20 periods with IAS15, deviation from rigid rotation 1·10⁻¹³ of the
  side; the “figure-eight” — the initial condition refined by Gauss–Newton shooting from the
  6-digit values of Šuvakov & Dmitrašinović 2013: p₁ = 0.3471168881189257,
  p₂ = 0.5327249453880308, T = 6.3259139829262061, residual over one period 1.9·10⁻¹⁵; it agrees
  with Simó's condition (8 digits) to 10⁻⁹ and to 3·10⁻⁸ in the period; 10 periods — return to the
  start 2·10⁻¹⁴, energy 3·10⁻¹⁶; Burrau's Pythagorean problem — a close encounter of 4.1·10⁻⁴ at
  t = 15.83, by t = 70 the bodies with GM 4 and 5 form a binary and the body with GM 3 escapes (as
  in Szebehely & Peters 1967); IAS15 and RKN12(10) agree to 4·10⁻¹² at t = 10, and by t = 70 chaos
  amplifies rounding to 2·10⁻² with the same outcome; energy is conserved to 10⁻¹⁰.
- **IAS15 corrected per the paper**: the “figure-eight” revealed that taking the minimum of τ over
  blocks stalls the step when the acceleration of one body passes through zero (the middle body
  of a symmetric collinear configuration, `MYB_E_NOCONV` at the first step). In the paper (Pham,
  Rein & Spiegel 2024, eq. 7) the norms are taken over the whole 3N-vector: the maximum over bodies
  of the Euclidean norms, separately for a, a′, a″ — this is what was implemented; the phase 3
  tests passed without changing the thresholds.
- **Kepler drift in increment form** (`kepler/drift.h`): Δr = (f − 1) r₀ + g v₀,
  Δv = ḟ r₀ + (ġ − 1) v₀ with f − 1 = −U₂/r₀, ġ − 1 = −U₂/r without subtraction from one;
  `myb_kepler_propagate` now returns r₀ + Δ (one rounding instead of three; error against
  MPFR ≤ 4.7·10⁻¹⁵, over 1–30 revolutions ≤ 2·10⁻¹⁵), and the transition matrix gives the same
  state bit for bit.
- **WH/WHFast** (`myb_wh`): Wisdom–Holman splitting in canonical Jacobi coordinates
  (Keplerian part with GM = η_i, drift — the universal propagator for all conic sections;
  kick — the full N-body acceleration converted to Jacobi coordinates, plus the Keplerian term;
  the center of mass drifts freely, its acceleration is exactly 0 for a mutual graph). The Jacobi
  state is accumulated with compensation. Schemes: DKD (= SABA1), SABA2–4 (Gauss–Legendre nodes
  and weights), SBAB1–4 (Gauss–Lobatto; the acceleration of the last kick is reused by the first
  kick of the next step). The coefficients come from closed forms, verified in MPFR bit for bit;
  the quadratures are exact exactly up to degree 2n − 1 (`test_mpfr_wh`).
- **Symplectic correctors** (orders 3, 5, …, 17) are derived here: to first order in the
  perturbation DKD = exp(hA + φ(h ad_A) hB), φ(t) = (t/2)/sinh(t/2); the corrector C = exp(χ),
  χ = −h ((1 − φ(t))/t) B is an odd series with Bernoulli numbers; C is assembled from
  Z(a, b) = A(a)B(b)A(−2a)B(−b)A(a) as Z(a_m, b_m)…Z(a_1, b_1)Z(−a_1, −b_1)…Z(−a_m, −b_m)
  (C⁻¹ is the same with −b), a_i = i·√(7/40) (the choice of Wisdom et al.), b is the solution of a
  linear system in MPFR. For order 3, b₁ = 1/(96·√(7/40)) = 0.0249005960277998675 — the published
  value, i.e. the derivation reproduces Wisdom, Holman & Touma 1996. The corrector is applied only
  to output states; the kernel steps uncorrected.
- **WH check** (`test_whfast`): two bodies — the Kepler flow, 2000 steps — 1·10⁻¹¹ a (rounding
  noise of the drifts: each is accurate to fractions of an ulp in velocity, the energy
  random-walks as with IAS15); observed orders at a small perturbation (ε = 10⁻⁷): DKD 2.00,
  correctors 3/5 — 4.00/5.81, SABA2/SBAB2 — 3.93/3.96, schemes of orders 6–8 reach rounding
  already at 16 steps per revolution; Jupiter–Saturn, 100 Saturn revolutions, 50 steps per
  Jupiter revolution against IAS15: DKD 9·10⁻⁵, corrector 17 — 1·10⁻⁵ (in position the ε²h² term
  remains, which the corrector does not touch; energy 1.9·10⁻⁷ → 1.9·10⁻¹⁰), SABA4/SBAB4 — 9·10⁻⁷;
  5000 revolutions with corrector 17 — the energy is bounded (≤ 7·10⁻¹⁰); 3000 steps forward and
  back — 1.3·10⁻¹³; dense output inside a step has the same accuracy as the nodes; a test body;
  identical bits with and without the thread pool.
- **Reproducibility**: added to the golden hash: N-body accelerations (with exclusion, 1 and 4
  threads — compared bitwise), invariants, the IAS15 “figure-eight”, WH schemes with corrector and
  dense output; the new golden 0xde3d79ea60027aae is identical at -O0/-O2/-O3, x86-64-v2/v3 and
  under UBSan.

### Phase 3 decisions

- **DOP853** — the coefficients were taken as text from SciPy (BSD-3), the structure of the error
  estimate and of the dense output follows Hairer; the state and time are updated with
  compensation (TwoSum), memory is allocated once at creation. The 7th-order dense output is built
  lazily (3 extra right-hand-side evaluations). A step shortened to the `t_limit` boundary does
  not reduce the next one: the planned step is kept (in Hairer's code the next step was computed
  from the shortened one). Verified: order conditions (all trees up to order 5, quadratures up to
  order 8), observed order 7.4 → 7.8 (converging to 8), 10 revolutions e = 0.7 at rtol 10⁻¹³ —
  1.3·10⁻¹⁰, dense output 2.4·10⁻¹², events 7.7·10⁻¹⁰ s, the transition matrix from the
  variational equations against central differences — 3·10⁻⁹ (the limit of the differences). Over
  10⁴ revolutions the error is 2·10⁻⁵: the method is not symplectic, the energy drift is linear —
  long arcs need IAS15.
- **IAS15** is implemented from the papers (Rein & Spiegel 2015; Everhart 1985), not from the
  REBOUND code (GPL). The Gauss–Radau nodes and the C, D, R tables are computed in MPFR (512 bits);
  the test `test_mpfr_radau` recomputes each element and compares bit for bit. The step criterion
  is PRS23 (Pham, Rein & Spiegel 2024, eqs. 15–17): τ² = 2|a|²/(|a′|² + |a||a″|),
  dt = (7!ε)^{1/7}·τ; the norms of a, a′, a″ are the maximum over component blocks (body =
  3 components) of the block's Euclidean norm, as in eq. 7 of the paper (before phase 4 there was
  a minimum of τ over blocks here — see “Phase 4 decisions”). The prediction of b for the next step
  is the re-expansion of the step polynomial plus the previous step's correction scaled by
  q^{k+1}. Iterations stop at 10⁻¹⁶ or when the correction stops decreasing (≤ 12 passes);
  non-converged ones are counted in the statistics. Verified: local order 15.6; 10 revolutions at
  e = 0…0.99 — from 2·10⁻¹⁵ to 1.6·10⁻¹³ of the period; 10⁴ revolutions e = 0.3 — 2.8·10⁻¹¹ of the
  period (DOP853 on the same arc — 2·10⁻⁵); a hyperbola; Brouwer's law (ensemble, e = 0.5): slope
  0.44, 4.6·10⁻¹⁵ per 1000 revolutions; a velocity-dependent force (charge in fields) — 3·10⁻¹⁶.
- **IAS15 limitations** found during checks: (1) τ from PRS23 is relative to the norm of a and
  tends to zero where the acceleration of the whole state vanishes — a one-dimensional oscillator
  ends with `MYB_E_STEP_MIN` (for such problems use DOP853); a system of bodies does not encounter this. (2) On an exactly circular orbit with a constant step the rounding errors of neighboring steps
  are identical and do not average out: the energy drifts ~10⁻¹⁷ per revolution (3·10⁻¹⁴ over
  3000 revolutions; independent of ε — this is rounding, not truncation). In asymmetric systems
  Brouwer's law holds.
- **Event location** is shared by the integrators: Illinois on the dense output down to adjacent doubles.
- **Symplectic compositions** (`myb_symplectic`): the “kick–drift–kick” leapfrog and Yoshida's
  triple jump (orders 2, 4, 6, 8; 1, 3, 9, 27 evaluations per step). Time is a coordinate of the
  extended phase space, so the map is symplectic for a time-dependent potential as well;
  velocity-dependent forces are not covered here (for them — IAS15). The weights are computed in
  MPFR and verified bit for bit (`test_mpfr_yoshida`). Verified: observed orders
  2.00 / 4.00 / 6.05 / 7.98; the energy error is bounded — the same over revolutions 1–100 and
  4901–5000 (for DOP853 over the same lengths it grows linearly); 10 revolutions forward and back
  — 1.7·10⁻¹². Yoshida's solutions with smaller error constants (7 and 15 stages) require solving
  the order equations from scratch — left for later. Assigned to phases or closed: see "Old debts (phases 3–11) and the bindings".
- **8th-order Gauss–Jackson** (`myb_gauss_jackson`) in summed ordinate form.
  The coefficients were re-derived from the operator identities hv − s = P₁(D)h²a,
  x − S = P₂(D)h²a with P₁ = 1/D − 1/(eᴰ − 1), P₂ = 1/D² − eᴰ/(eᴰ − 1)² (Bernoulli numbers) on
  the 8th-degree interpolant of the window u = −4…4 — as exact rational numbers; the predictor row
  matched the published one (Berry & Healy 2004: 3250433/53222400, …, 103798439/159667200), and
  the test `test_mpfr_gauss_jackson` recomputes all 261 numbers bit for bit. The sums use
  compensation. Startup is an iteration of the formulas themselves over the nodes −4…4 from a
  Taylor approximation (the force must be defined 4 steps back as well). Dense output is the exact
  integral of the window interpolant, events use the shared search. There is no step control; the
  statistics report the predictor–corrector difference as an estimate of the local error. Found
  during checks: an iteration “stall” cannot be treated as convergence if the correction is large
  (with too long a step the startup otherwise “converged” to garbage) — a stall is now accepted
  only below 10⁻¹⁰. Verified: observed order in positions ~10 (the window is exact for an
  8th-degree acceleration — positions are exact up to degree 10), 100 revolutions at 200 steps per
  revolution — 8.9·10⁻¹³ of the period at 2.55 evaluations per step, Lorentz force — 6·10⁻¹⁷,
  backward integration.
- **RKN12(10)** (Dormand, El-Mikkawy & Prince 1987) for x″ = a(t, x): the coefficients are the
  exact fractions of `DPRKN12Tableau` from OrdinaryDiffEq.jl (MIT); the test `test_mpfr_rkn12`, in
  GMP rational arithmetic, verifies that each double is the nearest to its fraction and evaluates
  the order conditions anew: 12th-order quadratures, b = b′(1 − c), moments of the embedded
  10th-order solution, row sums c²/2 (rows 14–17 in the source are 50-digit values; they agree to
  10⁻⁵³). Step control is 0.9·err^(−1/11) on positions and velocities. The pair has no
  interpolant of its own, so dense output is one more step of the pair from the start of the last
  step (16 evaluations), as accurate as the steps themselves. Verified: observed order 12.4;
  10 revolutions e = 0.9 at rtol 10⁻¹⁴ — 4.8·10⁻¹⁵ of the period in 758 steps; events
  1.5·10⁻¹⁰ s.
- **Statistics counters and step budgets** are `int64_t`: `long` on Windows is 32-bit and
  would overflow on long runs (10⁹ revolutions).
- **Two-body transition matrix** (`myb_kepler_stm`) in closed form: the Lagrange coefficients
  differentiated through the universal Kepler equation (derivatives of U_n with respect to α at
  fixed χ: −(χU_{n+1} − nU_{n+2})/2; c₄, c₅ are needed — series for |z| < 4). The same paths as
  in the propagator: removed whole periods give an explicit secular term (the dependence of the
  period on energy), a hyperbola far from pericenter is a composition of two arcs from
  pericenter, and the inverse of the symplectic matrix is obtained exactly (by transposing
  blocks). The reference is central differences of a 400-bit propagation in MPFR: ≤ 3·10⁻¹⁵ (in
  units of r₀, v₀), including over 1–30 revolutions; symplecticity 3·10⁻¹⁴; DOP853 with
  variational equations (rtol 10⁻¹³) agrees with it to 2·10⁻¹³.
- **Found while checking the matrix and fixed in the phase 2 propagator**: α = 2/r₀ − v²/μ and the
  period were computed in double; cancellation of terms (up to 2a/r₀ ≈ 39 at e = 0.95) left a few
  ulp in the period, and each removed period added a phase error — 2.8·10⁻¹¹ after 30
  revolutions. Now |r₀|², |v₀|², α and the period are in double-double: 1–30 revolutions
  ≤ 1.9·10⁻¹⁵, one revolution ≤ 4.6·10⁻¹⁵ for all conic sections (was ≤ 1.1·10⁻¹⁴). The MPFR test
  has been extended with multi-revolution cases, and the thresholds tightened to 10⁻¹⁴.
- **Reproducibility**: added to the golden workload: DOP853 (J2) and IAS15 (three bodies) with
  dense output and events, the two-body transition matrix, compositions of orders 2–8,
  Gauss–Jackson, RKN12(10); the new golden 0x49a0d3f351af2744 is identical at -O0/-O2/-O3,
  x86-64-v2/v3 and under UBSan.

### Phase 2 decisions

- **19 CORE-MATH correctly rounded functions** (tan, atan, atan2, asin, acos, exp, expm1,
  log, log1p, sinh, cosh, tanh, asinh, acosh, atanh, cbrt, hypot + sin, cos) are imported by the
  script `third_party/core-math/import.cmake` with three mechanical edits; the “impossible hard
  case” path in atan2 (upstream: printf + exit) became a fail-fast check.
  `pow` is deferred (it needs additional headers and is not required in phase 2).
- **Kepler equation without cancellation**: the residual `(1−e)E + e·(E − sin E) − M` with a
  Stumpff series; roots within 2 ulp of the exact one (MPFR) for e up to 1 − 2⁻⁵². Parabolic —
  closed form + one Newton correction (≤ 1.2 ulp).
- **One guaranteed solver** (`num/root.c`): bracketed Newton with the rtsafe rule. The first
  version returned the midpoint of the bracket when the Newton step was smaller than an ulp
  (`next == lo`); found by cross-checking against the path through elements, fixed, covered by
  tests.
- **Propagator**: universal variables; ellipses — removal of whole periods in double-double;
  hyperbolas far from pericenter (|H₀| > 1) — recomputation from pericenter (otherwise the terms
  of the equation are ~10⁴ times larger than the result). Near-parabolic — the direct path (there
  α is ill-conditioned, and the direct path hardly needs it). Cross-check against MPFR at 400 bits:
  ≤ 2·10⁻¹⁴.
- **Lambert** (Izzo 2015) with the guaranteed solver instead of Householder; an independent check
  serves as the second method (Gooding): each solution is propagated to r₂ by the propagator.
- **Not done in phase 2**: the analytic two-body state transition matrix (Shepperd 1985) —
  moved to phase 3 together with the variational equations of the integrators (done in phase 3).
- **Conditioning of the elements**: classical and equinoctial elements of near-parabolic orbits
  far from pericenter are ill-conditioned (a shift of ν by one ulp changes r by a factor of
  (1 − e)⁻¹); this is a property of the parametrization; in the tests this region is excluded and
  documented.

### Phase 1 decisions

- **Data format — NAIF text kernels.** A documented format in which JPL already distributes its
  kernels; one loader for our conventions and for `gm_de440.tpc` etc.
- **Time scales are run-time values, not C types.** The plan assumed a separate type per scale,
  but scales are defined by data (the core does not know the word “TT”), so an instant is
  `myb_time`, and a scale is a `myb_scale` descriptor of a `myb_timesys` system; mixing scales is
  caught by an explicit conversion through the tree, not by the compiler.
- **UTC is calendar labels only.** UTC is not a continuous scale; it is not in the tree, and
  conversions go through a table with an expiration date (after it — `MYB_E_RANGE`).
- **TDB−TT — the Fairhead & Bretagnon series** (a few ns); the exact TCB−TCG relation from the
  ephemeris will come with the `rel` module (phase 6). Topocentric terms belong to the observer's
  proper time and are not part of the scale convention.
- **Own number parsing**: Clinger's fast path + a big-integer check relative to midpoints;
  `strtod` is not used (it depends on the locale and the CRT).
- **Number output (Ryu)** is deferred until the phase in which writing of results appears (CCSDS
  OEM).
