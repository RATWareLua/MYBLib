# Site displacement (IERS 2010, chapter 7)

**English** | [Русский](site-displacement.ru.md) · [All passports](index.en.md)

## What it is

What the tides and the loads add to the place a point of the crust has in the body-fixed
frame: the tide of the solid body (degrees 2 and 3, steps 1 and 2) and its permanent part, the
pole tide, the load of the ocean tides carried to 342 waves by the admittance, the load of the
atmospheric tides S1 and S2, and the load of the ocean pole tide. The formulas are those of the
IERS Conventions 2010, chapter 7, written from the equations; every number comes from a kernel
(`data/frames/iers2010-site.tk` for the Earth). A station of [estimation](estimation.en.md) can
carry this motion, so that its measurements see the moved place.

## Source

- IERS Conventions 2010, chapter 7 (G. Petit, B. Luzum (eds.), IERS TN 36, version of
  2018-02-01): §7.1.1 eq. 5, 6 (degrees 2 and 3, Love and Shida numbers depending on the
  latitude), eq. 8, 9 (l⁽¹⁾ terms), eq. 10, 11 (out-of-phase parts), eq. 12, 13 with Tables
  7.3a, 7.3b (frequency dependence), eq. 14 (permanent part); §7.1.2 (ocean loading by the
  admittance); eq. 19 (S1, S2); §7.1.4 eq. 21, 26 (secular pole, pole tide); §7.1.5 eq. 29
  (ocean pole tide load).
- The kernel holds numbers only, with the SHA-256 of their sources: the frequency-dependence
  rows of the reference routines `STEP2DIU.F` and `STEP2LON.F` (31 diurnal rows, 29 of them not
  zero, and 5 zonal — more than the 11 and 5 the tables of the text print; K1 out of phase is
  −0.80 mm in the routine and −0.78 in the text), the 342 waves of Cartwright, Tayler and Edden
  from `ADMINT.F` (HARDISP), copied without change. The IERS software is not ported (its
  licence is not BSD); its test cases are the references.
- S1/S2 coefficients: the model of R. D. Ray and R. M. Ponte (*Ann. Geophys.* 21, 1897, 2003)
  as given by the loading calculator of the IERS (T. van Dam) or by TU Wien (the VMF data
  server); the ocean pole tide coefficients of `opoleloadcoefcmcor.txt` (Desai 2002).

## Interface

Header [`myb/site.h`](../../include/myb/site.h):

- the model `myb_site_tides_create`, `myb_site_tides_destroy`, `myb_site_tides_waves`,
  `myb_site_tides_given_waves`;
- ocean loading `myb_ocean_loading_create`, `myb_ocean_loading_load`,
  `myb_ocean_loading_load_file`, `myb_ocean_loading_destroy`, `myb_ocean_loading_at`;
- atmospheric loading `myb_site_atmosphere_load`, `myb_site_atmosphere_load_file`,
  `myb_site_atmosphere_grid_load`, `myb_site_atmosphere_grid_load_file`;
- grids of coefficients `myb_site_grid_load`, `myb_site_grid_load_file`,
  `myb_site_grid_destroy`, `myb_site_grid_values`, `myb_site_grid_at`;
- the displacement `myb_site_displacement` with `myb_site_input`, `myb_site_loads` and the
  flags `MYB_SITE_SOLID`, `_FREQUENCY`, `_POLE`, `_OCEAN`, `_ATMOSPHERE`, `_OCEAN_POLE`, `_ALL`,
  `_MEAN_TIDE`.

A moving station: `myb_station.motion` with `myb_station_motion` in
[`myb/estim.h`](../../include/myb/estim.h).

## Inputs and outputs

- `myb_site_input`: GM of the body (m³ s⁻²); the raising bodies (GM and place in the
  body-fixed frame, `myb_tide_raiser` of [`myb/grav.h`](../../include/myb/grav.h)); `t` — the
  argument time on UT1 (χ = GMST + π, the day of S1 and S2, the date of the secular pole);
  `dynamical` = TT − UT1 (s) — the arguments of the Moon and the Sun are taken on TT; polar
  motion xp, yp (rad); the flags.
- Output: a vector in the body-fixed frame, m. The solid tide and the pole tide are along the
  radius and across it (geocentric latitude); the loads are given along east, north, up, so
  they need the axes `enu` (rows, body-fixed). The place is "conventional tide free" (as in
  the ITRF) unless `MYB_SITE_MEAN_TIDE` takes the permanent part out.
- Ocean loading: the amplitudes (m) and phase lags (degrees) of the given waves (M2 S2 N2 K2
  K1 O1 P1 Q1 Mf Mm Ssa, up to 16) in three rows — up, west, south — from a BLQ record or
  numbers; the admittance A e^−iφ / |H| is interpolated in frequency within each band
  (long-period, diurnal, semidiurnal) by a cubic spline whose end slopes are those of the
  parabola through the three end points, constant outside the given waves.
  `myb_ocean_loading_at` returns up, west, south.
- S1/S2: twelve coefficients (up, east, north; cos and sin of one and two cycles a day of UT1,
  m), from a record of the IERS loading calculator (mm: radial, north–south, east–west), a line
  of the TU Wien station files (18 numbers in mm, S3 not read) or the TU Wien 1° grid.
- Ocean pole tide: six coefficients u_r^R, u_r^I, u_n^R, u_n^I, u_e^R, u_e^I of the point,
  bilinear on the grid; the wobble m₁ = xp − xs, m₂ = −(yp − ys) about the secular pole of
  eq. 21 serves both pole tides.
- The moving station: the raising bodies, UT1, the pole and TT come from the environment; with
  `rate` the velocity of the displacement too (up to 0.02 mm/s), from the displacements 30 s
  before and after.

## Domain of validity

- `myb_site_displacement` refuses with `MYB_E_DOMAIN`: a point or a raiser at the centre, a
  point on the axis of rotation (no longitude), gm not positive, input not finite, loads
  without the axes.
- `myb_site_tides_create`: `MYB_E_NOT_FOUND`, `MYB_E_FORMAT` (the variable in `diag`),
  `MYB_E_DOMAIN` (prefix longer than 20 characters), `MYB_E_NOMEM`. Ocean loading from
  numbers: `MYB_E_DOMAIN` for a negative or non-finite value.
- Readers: `MYB_E_NOT_FOUND` (no such point; names compared without case and surrounding
  blanks), `MYB_E_FORMAT` (the line), `MYB_E_NOMEM`, `MYB_E_IO`.
- `myb_site_grid_at`: `MYB_E_RANGE` outside the grid (the TU Wien S1/S2 grid runs from 89.5° to
  −89.5°: a point nearer a pole is outside), `MYB_E_DOMAIN` not finite; a grid that goes round
  the body continues across its first meridian.

## Accuracy

| What | Reference | Result |
|---|---|---|
| Solid tide, steps 1 and 2, the four cases of `DEHANTTIDEINEL.F` | the routine's arithmetic, done apart in double without its precession addition | **1.1·10⁻⁷ m** (up to 18 m in the fourth case, with the Sun at 0.1 au) |
| The same, cases 1–3 | the values printed in the routine's header | 1.7–2.8·10⁻⁵ m (see Limitations) |
| Permanent part | eq. 14 (printed to 0.1 mm) | below 10⁻⁴ m |
| Pole tide | eq. 26–28 with the secular pole of eq. 21 | 10⁻¹⁸ m |
| Ocean loading, Onsala and Reykjavik, 24 h, up to 38 mm | `HARDISP.F` (its six decimals) with its lunar time | **5·10⁻⁷ m**; with the argument of the Conventions 8·10⁻⁶ m |
| S1/S2 load | eq. 19 | 10⁻¹⁷ m |
| TU Wien S1/S2 grid at WETTZELL, ALGOPARK, AIRA, AGGO | the station coefficients of TU Wien | 4.9·10⁻⁷ m of up to 8.3·10⁻⁴ m |
| Ocean pole tide load, six dates at 43.75° S, 232.25° E | the test file of the Conventions (seven digits) | 5·10⁻⁷ of the value |
| Velocity of the displacement of a moving station | the difference of places a second apart | 4·10⁻¹⁰ m/s |

Also checked: the moving station against the displacement computed by hand from the same
bodies, UT1 and pole (test bound 5·10⁻⁹ m); M2 alone, 10 mm up, swings by 10.0–10.5 mm (the
admittance falls to nothing at N2 and S2); a range to a point above the station changes by the
up displacement to 10⁻⁶ m. The arguments of the Moon and the Sun on TT and the angle of
rotation on UT1 matter: on UT1 alone the phase of the lunar waves lags by 4·10⁻⁴ rad.

Test: `test_site` ([tests/unit/test_site.c](../../tests/unit/test_site.c)); the TU Wien grid
row runs only with `s1_s2_s3_cm_noib_grid.dat` downloaded.

## Cost

Not measured. A moving station with `rate` takes two more evaluations of the bodies and of the
orientation per instant (the header).

## Limitations

- The library follows the chapter where the reference routines depart from it.
  `DEHANTTIDEINEL` adds the precession since J2000 to the mean longitude of the Moon after
  taking the lunar time from it, which turns every diurnal wave by its multiple of 1.4° a
  century (for K1, 0.03 mm of its 12 mm) — hence the 1.7–2.8·10⁻⁵ m from its printed values;
  the fourth case of its header repeats the numbers of the third. `HARDISP` takes the lunar
  time as 360° UT − D, the mean longitude of the Sun without the aberration (20.5″) and on TT:
  in 2009 its angle is 23.6″ (1.57 s of rotation) ahead of χ.
- The admittance is interpolated within a band and held constant outside the given waves; S3
  is not read (eq. 19 has S1 and S2).
- Only the conventional terms of the chapter are modelled; non-tidal loading is not part of
  the model.
- The ocean pole tide was checked on nine nodes of the grid (an excerpt), not the whole file.

## Data

- `data/frames/iers2010-site.tk` (in the repository): the numbers of chapter 7 with the
  SHA-256 of the chapter and of the routines they come from.
- The user's: BLQ files of the Onsala ocean loading service; S1/S2 coefficients (IERS loading
  calculator or TU Wien station files); `opoleloadcoefcmcor.txt` of the Conventions (not in
  FetchData).
- `s1_s2_s3_cm_noib_grid.dat` — the TU Wien S1/S2/S3 grid, 1°, centre of mass, in
  `cmake/FetchData.cmake` with its SHA-256.
- Test data: `tests/data/site/loading.blq` (the two records of the HARDISP case),
  `opole-excerpt.txt` (nine nodes), `atmosphere.atl` and `atmosphere_grid.dat` (made for the
  tests), `atmosphere_sites.dat` (lines in the form of the TU Wien station files).
