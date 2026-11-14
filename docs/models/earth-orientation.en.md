# Earth orientation

**English** | [Русский](earth-orientation.ru.md) · [All passports](index.en.md)

## What it is

The rotation between the celestial frame (GCRS) and the terrestrial frame (ITRS) of the Earth,
with its time derivative. The modern path is the CIO-based IAU 2006/2000A scheme of the IERS
Conventions 2010: the CIP coordinates X, Y and the CIO locator s from Poisson series, the Earth
rotation angle, the TIO locator s′ and polar motion — all coefficients read from kernels. The
observed Earth orientation parameters (EOP 20 C04 or finals2000A) are interpolated as the IERS
does and completed by the diurnal and semidiurnal ocean-tide and libration terms. The series
can be tabulated for speed. The classical equinox-based frames (IAU 1976 precession, IAU 1980
nutation, MOD, TOD, TEME, GMST 1982, TEME → PEF) serve TLE work and old data.

## Source

- IERS Conventions (2010), IERS Technical Note 36 (G. Petit, B. Luzum, eds.), ch. 5: eq. 5.3
  (W), 5.10 (C), 5.13 (s′), 5.15 (ERA), 5.16 (series form), 5.43–5.44 (fundamental arguments);
  N. Capitaine et al. (2003). The coefficients are Tables 5.2a, 5.2b, 5.2d of the IERS
  Conventions Centre, copied without change (SHA-256 of the sources in the kernel).
- EOP: the IERS EOP 20 C04 series and Bulletin A (finals2000A); the interpolation follows the
  IERS routine `interp.f`; the sub-daily terms are Tables 8.2–8.3 (ocean tides, 71 terms, from
  `interp.f`) and 5.1a/b (libration, 10 + 11 terms, from `PMSDNUT2.F`, `UTLIBR.F`) — only the
  numbers are used, no IERS code (its licence is not BSD-compatible).
- Equinox-based: precession of Lieske et al. (1977), obliquity and 1980 nutation (Seidelmann
  1982, theory of Wahr 1981), GMST of Aoki et al. (1982); numbers from ERFA `nut80.c`,
  `prec76.c`, `obl80.c`, `gmst82.c` (BSD-3). TEME as SPICE defines it for SPK type 10.
- The code is the library's own; the Poisson series go through one generic evaluator.

## Interface

Header [`myb/frames.h`](../../include/myb/frames.h):

- CIP/CIO: `myb_cip`, `myb_cip_create`, `myb_cip_destroy`, `myb_cip_xys`, `myb_cip_tabulate`,
  `myb_cip_era`, `myb_cip_sprime`, `myb_cip_c2i`, `myb_cip_terrestrial`;
- EOP: `myb_eop_values`, `myb_eop`, `myb_eop_create`, `myb_eop_destroy`, `myb_eop_load_c04`,
  `myb_eop_load_c04_file`, `myb_eop_load_finals`, `myb_eop_load_finals_file`,
  `myb_eop_set_tides`, `myb_eop_rows`, `myb_eop_at`;
- the orientation provider: `myb_cip_orientation`, `myb_orientation_cip` (the form an
  environment of a force model and the frame graph take);
- equinox-based: `myb_equinox`, `myb_equinox_create`, `myb_equinox_destroy`,
  `myb_equinox_precession`, `myb_equinox_obliquity`, `myb_equinox_nutation`,
  `myb_equinox_frame` (`MYB_EQUINOX_MOD`, `MYB_EQUINOX_TOD`, `MYB_EQUINOX_TEME`),
  `myb_equinox_rotation`, `myb_equinox_gmst`, `myb_equinox_teme_pef`.

Rotations (`myb_rotation`) and the frame graph are described in
[frames and geodesy](frames-geodesy.en.md); the time scales in
[time scales](time-scales.en.md).

## Inputs and outputs

- [terrestrial] = W(t) R3(ERA) C(t) [celestial], C = R3(−s) Bᵀ (eq. 5.10),
  W = R1(−y_p) R2(−x_p) R3(s′). A rotation maps coordinates: x_B = M x_A, with dM/dt in s⁻¹.
- `myb_cip_xys`: X, Y, s in rad at a TT instant, rates in rad per TT second.
- `myb_cip_era`: the angle in [0, 2π) at a UT1 instant, the fraction of the day kept apart;
  rate in rad per UT1 second.
- `myb_eop_values`: x_p, y_p (rad), UT1 − UTC (s), LOD (s), dX, dY (rad). dX, dY must refer to
  the precession-nutation model in use (IAU 2006/2000A for `iau2006-cip.tk`).
- `myb_eop_at`: the EOP at an instant of a continuous scale (e.g. TT) that the UTC table of the
  time system converts, and the same instant on UT1. Four-point Lagrange on the rows around the
  date in the caller's scale; UT1 is interpolated as UT1 minus that scale, so a leap second in
  the window makes no jump. A quantity missing at any node counts as zero (dX, dY, LOD of
  predictions). The sub-daily terms, if set, are added; their LOD follows from the UT1 series.
- `myb_cip_terrestrial`: dX, dY added to X, Y; the derivative is with respect to TT and
  includes the rates of X, Y, s and ERA rate · (1 − LOD/86400 s).
- `myb_orientation_cip`: t is TT (an instant of the scale `tt` of `timesys`); returns the
  celestial → terrestrial rotation, `ut` = UT1, and x_p, y_p for the tides.
- `myb_cip_tabulate`: values and rates of X, Y, s at nodes every `step` s of TT (0: an hour);
  a cubic Hermite between the two nodes around the instant; outside the table the series.
- Equinox-based: MOD = P [J2000], TOD = N [MOD], TEME (the TOD pole, x the mean equinox of date
  projected on the true equator), PEF = R3(GMST82) [TEME]. Times are TT (these conventions do
  not tell TT from TDB) except UT1 for the sidereal time; derivatives per second.

## Domain of validity

- The tests cover 1800–2200 for both schemes; nothing in the code refuses other dates.
- `myb_eop_at`: `MYB_E_RANGE` outside the EOP table or the UTC table, fewer than two rows on
  either side, or a missing x, y or UT1 − UTC; `MYB_E_NOT_FOUND` without a UTC table. The UTC
  table expires (see [time scales](time-scales.en.md)), and so does EOP coverage.
- Loading: `MYB_E_FORMAT` (diag gives the line; rows must have increasing MJD), `MYB_E_IO`,
  `MYB_E_NOMEM`; a failed load leaves the table unchanged.
- `myb_cip_create`, `myb_equinox_create`, `myb_eop_set_tides`: `MYB_E_NOT_FOUND`,
  `MYB_E_FORMAT` (diag names the variable), `MYB_E_NOMEM`; tide prefixes at most 12
  characters.
- `myb_cip_tabulate`: `MYB_E_DOMAIN` for tt1 not after tt0, a negative or non-finite step, or
  more than 2²⁶ nodes.

## Accuracy

| What | Reference | Result |
|---|---|---|
| X, Y, 3000 dates of 1800–2200 | `eraXy06` | 1.7·10⁻¹⁷ rad |
| s; s′ | `eraS06`; `eraSp00` | 1.5·10⁻¹⁹ rad; 2·10⁻²⁵ rad |
| ERA and the matrix | `eraEra00`, `eraC2tcio` | 3.6·10⁻¹³ rad — ERFA's rounding (it sums ~200 turns in one double) |
| ERA, 20 000 instants of 1700–2300 | the same expression in MPFR, 256 bits | ≤ 1.3·10⁻¹⁵ rad |
| Rates of X, Y, s and of the matrix | central differences | within 10⁻¹⁸ rad/s and 10⁻¹³ s⁻¹ (tolerances) |
| Nutation table, step 1 h (6 h) | the series | 4·10⁻¹⁶ rad, rates 4·10⁻¹⁹ rad/s (5·10⁻¹³ rad) |
| EOP with sub-daily terms | independent evaluation the IERS way (Python from the Fortran tables) | pole 2·10⁻¹⁴″, UT1 1.7·10⁻¹⁴ s, LOD 8·10⁻¹¹ s (finite-difference limit of the reference) |
| UT1 across 2016-12-31T23:59:60 | continuity | 2·10⁻¹¹ s |
| Site vectors in the ICRF, 2 sites × 10 dates of 2020 | JPL Horizons | 0.19 m, 1.2·10⁻⁵ m/s — a rotation by 3–11 mas, the difference between Horizons' ITRF93 model and IERS 2010 (through the ITRF93 binary PCK: 0.077 m) |
| Precession angles; nutation; matrices, 1800–2200 | `eraPrec76`, `eraNut80`, `eraPmat76`, `eraNutm80` | 10⁻¹⁷ rad; 7·10⁻¹⁷ rad; ≤ 5.5·10⁻¹⁶ |
| GMST82 | `eraGmst82` | 3.9·10⁻¹³ rad (both sum up to 1.7·10⁷ s in one double) |
| TEME axes | their definition, atan(tan Δψ cos(ε + Δε)) | < 10⁻¹⁵ rad (tolerance); the pole equals that of TOD exactly |
| Equinox-based derivatives | central differences | 4·10⁻¹⁷ s⁻¹ |

Tests: [`test_erfa_frames`](../../tests/erfa/test_erfa_frames.c),
[`test_erfa_equinox`](../../tests/erfa/test_erfa_equinox.c),
[`test_mpfr_frames`](../../tests/mpfr/test_mpfr_frames.c),
[`test_eop`](../../tests/unit/test_eop.c) (EOP, the nutation table),
[`test_earth_frames`](../../tests/unit/test_earth_frames.c) (the whole chain against Horizons).

## Cost

- X, Y, s with their rates from the series: 47 µs a date (phase 5; ERFA without the rates —
  76 µs); 37–39 µs in the measurement of phase 20.
- Through the hourly table: 0.05 µs a call — 700 times faster (phase 20; `test_eop` prints it).

## Limitations

- The EOP are observed values: the rotation is only as good as the table covering the date;
  predictions of finals2000A lack dX, dY, LOD, which count as zero.
- The argument χ = GMST + π of the sub-daily terms is taken from UT1; the IERS routines take
  it from the UTC MJD, which repeats during a leap second.
- The equinox-based frames are of the FK5 era (J2000 differs from the GCRS by the 23 mas frame
  bias); TEME follows SPICE and differs from R3(Δψ cos ε_A) [TOD] by the second-order
  Δψ Δε sin ε, ~10⁻⁹ rad.
- Horizons orients the Earth with ITRF93 (IAU 1976/1980 with EOP nutation corrections), not
  IERS 2010: comparisons with it show 0.19 m at the surface that is a model difference.
- GCRS ↔ BCRS and the local frames belong to [1PN relativity](relativity-1pn.en.md); the
  ITRF93 binary PCK path to [PCK, CK and SCLK](pck-ck-sclk.en.md).
- `myb_cip_tabulate` must not be called while another thread uses the model.

## Data

In the repository (NAIF text-kernel format):

- `data/frames/iau2006-cip.tk` — X 1600, Y 1275, s + XY/2 66 terms (2941 in all), 14
  fundamental arguments, ERA, s′ (IERS Conventions 2010 Tables 5.2a, 5.2b, 5.2d; facts).
- `data/frames/eop-tides.tk` — 71 ocean-tide and 10 + 11 libration terms (IERS Tables 8.2–8.3,
  5.1a/b; facts).
- `data/frames/iau1980.tk` — IAU 1976 precession, IAU 1980 obliquity and nutation (106 terms),
  GMST 1982 (from ERFA, BSD-3).

Downloaded (living series, refreshed with `-DREFRESH=ON`): `eopc04.1962-now` (IERS EOP 20 C04,
weekly) and `finals2000A.all` (IERS Rapid Service, daily). Excerpts for the tests are in
`tests/data/eop/`.
