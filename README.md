# MYBLib

**English** | [Русский](README.ru.md)

A library for simulating and computing the motion of objects in space.

`myb` is a C23 library for the motion of objects in space: N bodies, orbits about planets,
interplanetary transfers, free motion. The models follow the real laws of physics (including
general relativity in the post-Newtonian approximation wherever the effect is measurable), and
the results are reproducible to the bit on Windows and Linux. The core knows no particular
body: the Sun, the planets and spacecraft are given as data.

Status: **version 1.1.0.0, phases 0–11** — numerical foundation, correctly rounded functions,
data reading, time scales, the body model, the two-body problem, Lambert's problem, integrators
(DOP853, RKN12(10), IAS15, Gauss–Jackson, symplectic compositions), the N-body problem (direct
summation, interaction graph, threads, WH/WHFast, SABA, SBAB), SPK ephemerides and reference
frames (IAU 2006/2000A, EOP, WGCCRE, PCK, geodesy, orbital frames), general relativity in the
post-Newtonian approximation (EIH with PPN, Lense–Thirring, de Sitter, local frames), gravity
fields of extended bodies (spherical and ellipsoidal harmonics, polyhedra, mascons,
figure–figure), tides, non-gravitational forces, the NRLMSISE-00 atmosphere, the force model of
an object, variational equations and covariances, the circular restricted three-body problem,
the B-plane and targeting, the effect budget, attitude and rotational dynamics (exact free
rigid body, fluid core, lunar librations), Kustaanheimo–Stiefel and Sundman regularization, low
thrust, equinox-based frames and TEME, SGP4 and SPK type 10, compact binaries to 3.5PN with
spins, observables (light time, Shapiro delay, two-way measurements, aberration, light
deflection, delay in plasma, troposphere and ionosphere), magnetic fields, the Lorentz force,
plasma and the charging of bodies. Detailed plan, sources and roadmap:
[`docs/PLAN.en.md`](docs/PLAN.en.md); provenance of the data: [`data/README.md`](data/README.md).

## What is there

| Module | Contents | Verified |
|---|---|---|
| `myb/core.h` | error codes, version, memory allocator, diagnostics | |
| `myb/num.h` | TwoSum, TwoProd, double-double (+ − × ÷ √), Sum2, Dot2 | accuracy against MPFR: every operation within its proven bounds |
| `myb/cmath.h` | `myb_fma`, `myb_sin`, `myb_cos` — correctly rounded | `fma`: 12·10⁶ inputs bit for bit with hardware FMA; `sin`/`cos`: 1.27·10⁶ inputs (including the hard cases of CORE-MATH) bit for bit with MPFR |
| `myb/io.h` | correctly rounded parsing of numbers (independent of the locale and the C runtime); NAIF text kernels | 4.2·10⁵ strings bit for bit with MPFR, including exact midpoints between doubles and strings longer than 800 digits |
| `myb/time.h` | instants (integer seconds + fraction, resolution < 10⁻¹⁶ s), calendar, Julian dates; time scales TAI, TT, TCG, TDB, TCB, GPS, GST, BDT and UTC from data | against IAU SOFA (ERFA): calendar — 0, TDB−TT — 8·10⁻¹⁶ s, TT↔TCG and TDB↔TCB — 1.3·10⁻¹¹ s, UTC→TAI — 3.2·10⁻¹¹ s (the limit is SOFA's two-part JD resolution) |
| `myb/body.h` | parameters with uncertainties, GM from mass, the body (GM, kind of motion) | |
| `myb/cmath.h` | 17 more CORE-MATH functions: tan, atan, atan2, asin, acos, exp, expm1, log, log1p, sinh, cosh, tanh, asinh, acosh, atanh, cbrt, hypot | 1.44·10⁶ inputs bit for bit with MPFR, including samples of hard cases |
| `myb/kepler.h` | Kepler's equation (ellipse, hyperbola, parabola), anomalies, classical and modified equinoctial elements, universal-variable propagator, closed-form state transition matrix | roots ≤ 2 ulp from the exact ones (MPFR); propagator ≤ 5·10⁻¹⁵ from a 400-bit MPFR reference (also over 30 revolutions); transition matrix ≤ 3·10⁻¹⁵ from MPFR differences, symplectic to 3·10⁻¹⁴; Vallado's example 2-4 in all digits |
| `myb/lambert.h` | Lambert's problem (Izzo 2015), multi-revolution solutions | Curtis' example 5.2 in all digits; 5·10⁴ solutions checked by propagation to r₂: ≤ 7·10⁻¹² |
| `myb/integ.h` | DOP853 (8th order, 7th-order dense output, events) for first-order systems; RKN12(10) with step control for x″ = a(t, x); IAS15 (15th order, Gauss–Radau, PRS23 step criterion) for x″ = a(t, x, v); 8th-order Gauss–Jackson (summed form, dense output, events) for x″ = a(t, x, v) with constant step; Yoshida's symplectic compositions of orders 2–8 for x″ = a(t, x); compensated update of state and time, no allocation within a step | order conditions; observed orders 7.8, 15.6 and 2/4/6/8; IAS15 over 10 revolutions at e ≤ 0.99 — ≤ 1.6·10⁻¹³ of the period, Brouwer's law (energy error growing as ~t^0.44); bounded energy error of the compositions over 5000 revolutions; Gauss–Jackson: 100 revolutions — 9·10⁻¹³ of the period; RKN12: order 12.4, 10 revolutions at e = 0.9 — 5·10⁻¹⁵ of the period; constants of Gauss–Radau, Gauss–Jackson, RKN12 and Yoshida's weights bit for bit with MPFR/GMP |
| `myb/pool.h` | thread pool for parallel loops (created once; the result does not depend on the number of threads) | |
| `myb/dyn.h` | N-body system (GM only; massive and test bodies): direct summation of accelerations with Sum2, interaction graph (directed exclusions), parallel and vector paths with identical bits, invariants in double-double; Wisdom–Holman integrators in Jacobi coordinates: DKD with symplectic correctors of orders 3–17, SABA2–4, SBAB1–4, dense output | accelerations ≤ 3.5·10⁻¹⁶ of Σ\|terms\| (reference: double-double); bits identical for 1–13 threads and all paths; the figure-eight refined to a residual of 1.9·10⁻¹⁵ per period (T = 6.3259139829262); the Pythagorean problem — the outcome of Szebehely & Peters; Lagrange's solution — 10⁻¹³ over 20 periods; WH: observed orders 2/4/6, Jupiter–Saturn against IAS15, coefficients of schemes and correctors bit for bit with MPFR (the 3rd-order corrector matches Wisdom et al. 1996) |
| `myb/kin.h` | own SPK reader (NAIF DAF): types 1, 2, 3, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 21; several files with SPICE priorities, segment chains with compensated summation; files mapped into memory (DE441 — 3 GB) | all types against CSPICE N0067 (7 thousand states): positions ≤ 6·10⁻¹⁶, velocities ≤ 10⁻¹⁴ (the Hermite types are more accurate than CSPICE: 1.6·10⁻¹⁶ against its 2.3·10⁻¹⁴ from the exact interpolant); `testpo.440` of DE440 — 7·10⁻¹⁵ au with JPL's tolerance of 10⁻¹³; asteroid Eros (type 21 from Horizons) — Horizons vectors to 1.8·10⁻¹² |
| `myb/frames.h` | rotations with derivatives; the IAU 2006/2000A CIP/CIO model from the IERS tables (Poisson series as data), ERA, s′, the GCRS→ITRS matrix; EOP (EOP 20 C04, finals2000A) with interpolation across leap seconds and sub-daily corrections (ocean tides, libration); binary PCK (Moon DE440, ITRF93); WGCCRE rotation models from text PCK; fixed FK frames (TK frames) and the built-in inertial frames of NAIF; geodetic coordinates (Vermeille 2011); ENU, RTN, VNC, LVLH | X, Y, s against ERFA ≤ 1.7·10⁻¹⁷ rad, ERA accurate to 1.3·10⁻¹⁵ rad (MPFR; ERFA's to 3·10⁻¹³); EOP against an independent implementation of IERS (interp.f, PMSDNUT2, UTLIBR) — pole 2·10⁻¹⁴″, UT1 1.7·10⁻¹⁴ s; lunar librations of `testpo` — 1 ulp; WGCCRE, 18 bodies, against CSPICE ≤ 1.1·10⁻¹⁰ (CSPICE's limit; our W is accurate to 10⁻¹⁵ rad); geodesy — 3·10⁻¹⁶ rad at all heights (ERFA loses 2.6·10⁻¹¹ at 1.5·10⁷ m); sites on the Earth against Horizons — 0.08 m through the ITRF93 PCK, 0.19 m through IERS 2010 (frame difference 3–11 mas) |
| `myb/rel.h`, `myb/dyn.h` (EIH) | Einstein–Infeld–Hoffmann N-body equations with PPN β, γ (as in the JPL ephemerides) and their integrals (energy, momentum, relativistic barycenter); for a satellite — the Schwarzschild, Lense–Thirring and de Sitter terms (IERS 2010 eq. 10.12); BCRS↔local frame transformation (IAU 2000 B1.3); proper time of clocks | EIH against the Euler–Lagrange equations of the PPN Lagrangian (sympy, 50 digits) — 2.7·10⁻¹⁰ of the 1PN terms, bits independent of threads; 1PN energy conserved to O(c⁻⁴); Mercury's perihelion 42.9804″/century (theory — to 5·10⁻⁹); Lense–Thirring of LAGEOS 30.63 mas/yr, de Sitter 19.19 mas/yr; GPS clocks +4.4647·10⁻¹⁰; DE440 over a year: Mercury 106 m, Venus 10 m, Mars 68 m (Newton — 58, 99, 40 km) |
| `myb/grav.h` | spherical harmonics (ICGEM 1.0/2.0 with trends and periodic terms, degree up to 2700, the pole an ordinary point); IERS 2010 tides (Love numbers, frequency-dependent corrections, permanent tide, pole tides) and ocean tides (FES2004); homogeneous polyhedron (Werner–Scheeres); mascons; ellipsoidal harmonics (Lamé functions, Niven's form); the figure–figure term U₂₂ | MPFR 256 bits: degree 120 near the pole and 2190 at 62° — 7.5·10⁻¹⁶; prism against Nagy's formula (40 digits) — 7.6·10⁻¹⁶; tides against the addition theorem 2·10⁻¹⁶ and the IERS K1 example; the ellipsoidal series converges inside the Brillouin sphere (bar: 6.5·10⁻⁵ at degree 16, where the spherical series diverges); DE440 with the figures of the Earth, the Moon and the Sun's J2 — Moon 0.65 m over a year |
| `myb/force.h` | third body (Battin), radiation pressure on a sphere and plates, penumbral light fraction (exact geometry, limb darkening), planetary albedo and IR (Knocke), drag, empirical RTN | third body against 50 digits 4·10⁻¹⁶; shadow against the intersection of spherical caps 1.1·10⁻¹⁴; LAGEOS-1 against the ILRS orbit: **2.0 cm** RMS over a week |
| `myb/kin.h` (SP3) | SP3-c/d reader (IGS, ILRS) | |
| `myb/num.h` (linear algebra) | LU, Cholesky, least squares (Householder), eigenvalues (Francis QR), eigenvectors | roots of companion matrices, spectrum 1500 … 1/1500, trace and determinant of random matrices |
| `myb/model.h`, `myb/uncert.h` | force model of an object (named terms, Jacobians, parameters); variational equations (transition and sensitivity matrices), covariances, deterministic random number generator, Monte Carlo | transition matrix against the closed two-body formula — 5.5·10⁻¹⁵; linear covariance against 4000 samples — 3% |
| `myb/cr3bp.h` | circular restricted three-body problem: Lagrange points, Jacobi constant, correction of periodic orbits, monodromy and stability, manifolds, Richardson's halo orbits, pulsating system of real bodies | 56 orbits of the JPL catalog (8 families): state 9·10⁻¹², period 6·10⁻¹²; JWST near L2 against Horizons |
| `myb/transfer.h`, `myb/budget.h` | B-plane, patched conics, shooting (Newton, minimum norm, least squares); effect budget of a scenario | Lambert by shooting — to 4·10⁻¹² m/s of `myb_lambert`; Roadster over 7 years — 306 m, OSIRIS-REx through a gravity assist — 12.6 m against Horizons |
| `myb/num.h` (elliptic) | Carlson's integrals R_F, R_D, R_J, R_C; complete and incomplete K, E, Π at any φ; Jacobi functions am, sn, cn, dn | against mpmath (60 digits): ≤ 7.4 ulp, without loss of significance for m → 1 and n < 0 |
| `myb/frames.h` (equinox-based) | IAU 1976 precession, IAU 1980 nutation (106 terms from data), MOD, TOD, TEME, GMST 1982, TEME→PEF | against ERFA: nutation 7·10⁻¹⁷ rad, matrices ≤ 5.5·10⁻¹⁶, GMST82 3.9·10⁻¹³ rad (ERFA's rounding) |
| `myb/attitude.h` | quaternions, 12 Euler angle sequences, principal axes; Euler's equations, torques (point mass, field, figure–figure); exact free rigid body (Jacobi functions and Π); symplectic splitting of orders 2–8; body with a fluid core | free body against DOP853 (10⁻¹⁴) — 2·10⁻¹⁴ rad; symmetric — 7·10⁻¹⁵; orders 2.00/3.99/5.98/7.95; angular momentum balance 6.6·10⁻¹⁶; **lunar librations against DE440 — 0.04″ over a year** (0.6 m at the surface over 2 years) |
| `myb/regular.h` | Kustaanheimo–Stiefel regularization and Sundman transformation for DOP853, propagation to a given time | e = 0.999, 20 revolutions: KS 1.2·10⁻¹² in 4 200 evaluations against 1.3·10⁻⁸ in 75 000 for Cowell; passage through r = 0 |
| `myb/thrust.h` | constant-thrust arcs with mass flow, steering laws (velocity, inertial, RTN, Edelbaum), Gauss's equations in equinoctial elements, Edelbaum transfer | Tsiolkovsky's equation 1.7·10⁻¹⁴; radial thrust against the closed formula 3·10⁻¹³; Gauss against Cowell 5·10⁻¹¹; Edelbaum transfer: a — 1.4·10⁻⁷, i — 0.0025° |
| `myb/rel.h` (binaries) | relative motion of two bodies to 3.5PN in harmonic coordinates (1PN, 2PN, 3PN with the gauge length, radiation reaction at 2.5PN and 3.5PN), spin-orbit and spin-spin terms, spin-induced quadrupoles, spin precession; Peters–Mathews formula | coefficients against a machine transcription of the LaTeX of Blanchet's review 4·10⁻¹⁶; invariant K: 2PN 27/2 − 7ν, 3PN — 0.06%, independent of the gauge; PSR B1913+16 — Ṗ 3·10⁻⁵; PSR J0737−3039B — geodetic precession 5.0735°/yr (prediction 5.0734) |
| `myb/observe.h` | light time (reception, transmission) with the Shapiro delay (Moyer), range and range rate, two-way measurements, aberration (Lorentz), light deflection (PPN), RA/Dec | against CSPICE on DE440: 4·10⁻¹² s; aberration bit for bit with ERFA; delay against quadrature 9·10⁻¹⁵ |
| `myb/force.h` (charge, plasma) | magnetic fields from Gauss coefficients (IGRF tables), Lorentz force, co-rotation field, plasma stream | IGRF-14 against ppigrf 1.4·10⁻¹⁵ of \|B\|; helix 5·10⁻¹⁵; canonical momentum in a dipole 3·10⁻¹⁶ |
| `myb/force.h` (charging), `myb/observe.h` (plasma) | OML currents to a sphere in a drifting plasma, photoemission, floating potential; electron column of power-law profiles (corona, wind) and dispersive delay | currents against integration over the Maxwellian (mpmath) 6·10⁻¹⁶; potential against the solution through Lambert's W 10⁻¹⁵; column against mpmath 6·10⁻¹⁶ |
| `ext/nrlmsise00` | NRLMSISE-00 (separate target `myb::nrlmsise00`): thread-safe, reproducible port | 17 reference cases — to the 7 printed digits |
| `ext/troposphere` | tropospheric delay per IERS 2010 ch. 9 (target `myb::troposphere`): Saastamoinen, Mendes–Pavlis, FCULa, VMF1, GMF, GPT, gradients | test cases of the IERS routines — 10⁻¹⁶…10⁻¹⁵ |
| `ext/ionosphere` | IONEX maps and the single-layer model (target `myb::ionosphere`) | independent implementation of the IONEX interpolation 2·10⁻¹¹; pierce point 9·10⁻¹⁶ |
| `ext/sgp4` | SGP4/SDP4 (Vallado 2006/2020; target `myb::sgp4`): TLE parsing (Alpha-5, checksums), WGS-72/84, AFSPC and improved modes, pure thread-safe propagation; SPK type 10 | Vallado's `tcppver.out`: 700 lines to the printed precision (5·10⁻⁹ km), error codes identical; SPK type 10 against CSPICE N0067 — 2.7·10⁻¹² km (near-Earth), 1.9·10⁻⁹ km (deep space) |

Data (`data/`): the tree of time scales from the IAU resolutions, the TDB−TT series of Fairhead &
Bretagnon, the TAI−UTC table since 1960 (valid to 2027-06-28 per IERS Bulletin C 72), CODATA
2022, the IAU 2006/2000A CIP series and the sub-daily EOP corrections from the IERS 2010 tables,
the built-in inertial frames of NAIF, the IERS 2010 tides (ch. 6). Large sets (DE440, `testpo`,
PCK, EOP, GOCO06s, FES2004, IGRF-14) are downloaded by `scripts/fetch-data.cmd` (Linux:
`cmake -P cmake/FetchData.cmake`) into `data/download/` with SHA-256 checks; without them the
corresponding tests are skipped.

## Building

Requires Clang ≥ 19 and CMake ≥ 3.28. Clang only.

**Windows.** Clang (LLVM for Windows) and MSYS2 with the MINGW64 environment
(`C:\msys64\mingw64`: headers, runtime, `mingw32-make`; for the accuracy tests the packages
`mpfr`, `gmp`). Visual Studio is not needed. Another MSYS2 location: the environment variable
`MYB_MINGW_ROOT`.

```
cmake --preset windows-dev
cmake --build --preset windows-dev --parallel
ctest --preset windows-dev
```

**Linux (Gentoo and others).** `clang`, `ninja`; for the accuracy tests `mpfr`.

```
cmake --preset linux-dev
cmake --build --preset linux-dev --parallel
ctest --preset linux-dev
```

Presets: `*-dev` (Debug, `-Werror`), `*-release` (`-O2`, contract checks kept), `windows-ubsan`
(UBSan in trap mode); on Linux, sanitizers that stop at the first error: `linux-asan` (ASan +
UBSan + LeakSanitizer, `-O2`), `linux-asan-debug` (the same at `-O0` plus the pointer-pair checks
`pointer-compare`/`pointer-subtract`), `linux-ubsan` (UBSan with report and stack), `linux-lsan`
(standalone LeakSanitizer), `linux-tsan` (ThreadSanitizer: data races and lock order of the
thread pool). Results of the Gentoo checks: `docs/PLAN.en.md`, §12.

The test `test_repro` compares a hash of the results with its golden value: it must match on
every platform, at every optimization level and ISA (x86-64-v2/v3).

**CI** (`.github/workflows/ci.yml`, GitHub Actions): the same presets — Linux (Ubuntu 24.04,
Clang 20): `linux-dev`, `linux-release`, `linux-asan`, `linux-asan-debug`, `linux-ubsan`,
`linux-lsan`, `linux-tsan`; Windows (Clang + MSYS2 MINGW64): `windows-dev`, `windows-release`,
`windows-ubsan`. Data are downloaded by `FetchData.cmake` and cached; the live EOP series are
refreshed on every run. It runs only when a tag is pushed (a release) and on request
(`workflow_dispatch`), not on every push.

## Numerical contract

IEEE-754 binary64, round to nearest; the library never changes the floating-point environment.
The build forces off the contraction of `a*b+c` into FMA and fast-math, and a `-ffast-math` from
the user stops the compilation. Details: `docs/PLAN.en.md`, §7.

## License

BSD-3-Clause (`LICENSE`). Third-party code: `THIRD_PARTY_NOTICES.md`.
