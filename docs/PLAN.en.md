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

Outside the first version: contact and landing, continuum hypersonic aerodynamics (user tables only), ablation, stochastic micrometeoroid impacts. Later: observables (light time, Shapiro delay, aberration, light deflection).

**Result uncertainty** (`uncert`): variational equations with respect to state and parameters, covariance propagation through the STM, Monte Carlo with a deterministic RNG.

## 5. N bodies

- Exact direct summation only: O(N²) over massive bodies + O(N·M) for test particles; target — up to ~10⁴ massive bodies.
- Interaction graph: the level is set for each pair (point, harmonics, figure–figure, 1PN, tides); the default is point + 1PN.
- Threads — through the platform layer, partitioning over i, summation over j in a fixed order: the result does not depend on the thread count and matches the single-threaded symmetric version.
- SIMD — vectorization over i, not over j; no FMA in the kernels; bit for bit with the scalar code.
- Barnes–Hut gives a force error of ~10⁻³–10⁻⁴ — not included. FMM is a possible separate module, explicitly marked as approximate.

## 6. Architecture

```
L0  core, platform (files, threads, CPUID)
L1  num, cmath (own functions), io (DAF/SPK/PCK, .gfc, IERS series, CCSDS, CSV)
L2  time (general), body (body model, parameters ± σ)
L3  frames, kin (kinematics providers), kepler, integ
L4  grav (fields), rel (1PN, local reference systems, proper time), force (non-gravitational)
L5  dyn (bodies, test particles, rotation, interaction graph), uncert, budget
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
  charging.

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
  spacecraft charging models (the charge is set by the user).

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
  exponent.

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
  phase 10; a budget from the linear estimate ∫Φ a_k dt (currently — full propagations).

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
  table).

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
  equinox-based frames (TOD, TEME), local GR frames (GCRS↔BCRS — phase 6).

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
  the order equations from scratch — left for later.
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
