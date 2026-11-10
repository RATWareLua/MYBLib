# myb data

**English** | [Русский](README.ru.md)

The library contains no constant of any particular body or convention: everything comes from
files in the format of NAIF text kernels (`myb/io.h`). Each file describes its own sources and
accuracy; this is a summary.

| File | Contents | Source | License |
|---|---|---|---|
| `time/time-scales.tk` | tree of the time scales TAI, TT, TCG, TDB, TCB, GPS, GST, BDT | IAU resolutions 1991 A4, 2000 B1.9, 2006 B3; IS-GPS-200 | definitions (facts) |
| `time/tdb-tt-fb1990.tk` | TDB−TT series (787 terms + 5 JPL corrections) | Fairhead & Bretagnon 1990 as edited by SOFA, numbers from ERFA `dtdb.c` | BSD-3 (ERFA) |
| `time/utc-tai.tk` | TAI−UTC since 1960 (drift before 1972), validity date | USNO tai-utc.dat through ERFA `dat.c`; IERS `Leap_Second.dat` (Bulletin C 72, to 2027-06-28) | BSD-3 (ERFA); IERS data |
| `constants/codata-2022.tk` | G, mₑ, m_p, ε₀ (CODATA 2022), c, e, k (SI definitions) | NIST CODATA 2022; definition of the metre | facts |
| `frames/iau2006-cip.tk` | X, Y, s + XY/2 series of the IAU 2006/2000A CIP model (2941 terms), 14 fundamental arguments, ERA, s′ | IERS Conventions 2010, tables 5.2a, 5.2b, 5.2d (iers-conventions.obspm.fr; SHA-256 in the file) | IERS conventions (facts) |
| `frames/eop-tides.tk` | diurnal and semidiurnal EOP variations: ocean tides (71 terms) and libration (10 + 11) | IERS Conventions 2010, tables 8.2–8.3 and 5.1a/b, tables of the programs interp.f, PMSDNUT2.F, UTLIBR.F (SHA-256 in the file) | IERS conventions (facts) |
| `grav/iers2010-tides.tk` | tides of the Earth's field: Love numbers (table 6.3), step-2 frequency-dependent corrections (tables 6.5a–c, 71 terms, Poisson series), permanent tide, pole tides, secular pole | IERS Conventions 2010, ch. 6 and 7 (2018 revision; SHA-256 in the file) | IERS conventions (facts) |
| `frames/iau1980.tk` | IAU 1976 precession, IAU 1980 obliquity and nutation (106 terms), GMST 1982 | ERFA `nut80.c`, `prec76.c`, `obl80.c`, `gmst82.c` (SHA-256 in the file) | BSD-3 (ERFA) |
| `frames/naif-inertial.tk` | built-in inertial frames of SPICE (B1950, FK4, GALACTIC, ECLIPJ2000, ECLIPB1950) as TK frames | NAIF Frames Required Reading; Lieske 1977/1979; Explanatory Supplement 1992 | facts |

## Updating

- **Leap seconds**: when IERS Bulletin C announces one, add a row to `time/utc-tai.tk` and move
  `MYB_UTC_VALID_UNTIL_MJD` to the new expiry date from `Leap_Second.dat`. After that date the
  UTC conversions return `MYB_E_RANGE` instead of guessing.
- **Large sets** are not kept in the repository: `scripts/fetch-data.cmd` (on Linux
  `cmake -P cmake/FetchData.cmake`) downloads them into `download/` (not versioned). Fixed
  releases are checked against the SHA-256 in `cmake/FetchData.cmake`: DE440 (`de440.bsp`,
  NAIF), `testpo.440` and `header.440` (JPL SSD), the lunar orientation
  `moon_pa_de440_200625.bpc`, `pck00011.tpc` (NAIF), the GOCO06s field model (ICGEM, CC BY 4.0),
  the FES2004 ocean tides as Stokes coefficients (IERS), the IGRF-14 Gauss coefficients
  `igrf14coeffs.txt` (IAGA, NOAA NCEI). The live series are refreshed by
  `fetch-data.cmd refresh`: EOP 20 C04 (`eopc04.1962-now`, IERS EOP Centre), `finals2000A.all`
  (IERS Rapid Service), `earth_latest_high_prec.bpc` (NAIF). Tests that need them are skipped
  when the files are missing.

JPL files (`gm_de440.tpc`, `pck00011.tpc`, ...) are read by the same loader without conversion.
