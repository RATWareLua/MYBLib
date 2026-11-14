# SGP4/SDP4 and two-line element sets

**English** | [Русский](sgp4-tle.ru.md) · [All passports](index.en.md)

## What it is

The SGP4/SDP4 theory of the two-line element sets (TLE) of the US Space Force catalogue, in an
optional module: parsing and writing the two lines (Alpha-5 numbers, checksums), the
conversion between a TLE and a CCSDS mean-elements message (OMM), propagation of an element
set to a TEME state, SPK segments of type 10 (evaluated and written word for word as SPICE
writes them), and the fit of an element set to positions known otherwise. TLE are mean
elements, meaningful only with the theory that produced them; their accuracy is that of the
elements — about a kilometre at the epoch, growing by kilometres a day (the header). For
precise work fit a numerical orbit to the states.

## Source

- F. R. Hoots, R. L. Roehrich, Spacetrack Report #3 (1980), in the revision of D. A. Vallado,
  P. Crawford, R. Hujsak, T. S. Kelso, "Revisiting Spacetrack Report #3", AIAA 2006-6753, and
  its code of 12 March 2020 (distributed by CelesTrak for free use).
- `sgp4.c` is a C translation made from the Python translation of that code in python-sgp4
  by Brandon Rhodes (revision `8e234f81…`, MIT licence); the module is under the same MIT
  licence ([`ext/sgp4/NOTICE`](../../ext/sgp4/NOTICE)). The port follows the code operation
  by operation, with correctly rounded elementary functions and an exact `fmod` (one FMA);
  the deep-space integrator cache is not kept (the 720-minute steps are fixed, the values are
  the same), so propagation is a pure function and a record serves any number of threads.
- SPK type 10: NAIF SPK Required Reading; the words as `SPKW10` lays them out. The element
  set fit (`tle_fit.c`) is myb's own.

## Interface

Header [`myb_sgp4.h`](../../ext/sgp4/myb_sgp4.h) (target `myb::sgp4`, not part of libmyb):

- TLE `myb_tle_parse`, `myb_tle_epoch`, `myb_tle_format`, `myb_tle_from_omm`,
  `myb_omm_from_tle` (the OMM of [`myb/ccsds.h`](../../include/myb/ccsds.h)); the flag
  `MYB_TLE_NO_CHECKSUM`;
- the theory `myb_sgp4_constants_of`, `myb_sgp4_init`, `myb_sgp4_init_elements`,
  `myb_sgp4_propagate`, `myb_sgp4_propagate_km`, `myb_sgp4_gmst`; enumerations
  `myb_sgp4_gravity` (WGS-72 old, WGS-72, WGS-84), `myb_sgp4_mode` (improved, AFSPC),
  `myb_sgp4_error`;
- the fit `myb_tle_fit`, `myb_tle_from_state`, options `myb_tle_fit_options`;
- SPK type 10 `myb_sgp4_spk10_state`, `myb_sgp4_spk10_constants`, `myb_sgp4_spk10_packet`,
  `myb_sgp4_spk10_size`, `myb_sgp4_spk10_words` (the words go to `myb_spk_writer_add_raw` of
  [`myb/kin.h`](../../include/myb/kin.h));
- structures `myb_tle`, `myb_sgp4`, `myb_sgp4_constants`.

## Inputs and outputs

- `myb_tle`: catalogue number as written (5 characters, Alpha-5: A = 10 … Z = 33 without I and
  O) and decoded; epoch year (two-digit 57…99 → 19xx, 00…56 → 20xx) and day of the year
  (1.0 = 1 January 0h UTC); ṅ/2 (rev/day²), n̈/6 (rev/day³), B* (1/Earth radii); angles in
  degrees; the Kozai mean motion in rev/day.
- States are in TEME (true equator, mean equinox of the epoch of the state; see
  [frames](frames-geodesy.en.md), `myb_equinox_teme`) in m and m/s; `myb_sgp4_propagate_km`
  gives km and km/s as the verification program prints them. Time: minutes since the epoch of
  the record. `myb_sgp4_init_elements` takes elements in the units of the theory (epoch in days
  since 1949 December 31 0h UTC, rad, rad/min).
- Constants: WGS-72 old (xke = 0.0743669161), WGS-72 (the catalogue's, default), WGS-84;
  modes improved (Vallado's corrections) and AFSPC. The deep-space branch is taken for periods
  of 225 min and more.
- OMM ⇄ TLE: the numbers pass unconverted, in the units of the catalogue; the epoch is a
  calendar date in the message and a day of the year in the set, exact to the rounding of one
  division. From a set the message gets centre EARTH, frame TEME, time UTC, theory SGP4.
- SPK type 10: eight constants (J2, J3, J4, KE, QO, SO, ER, AE), then per set a packet of 14
  words (NDT20, NDD60, BSTAR, INCL, NODE0, ECC, OMEGA, MO, NO in rad and rad/min, the epoch in
  TDB seconds past J2000, four nutation values). The state at et blends the SGP4 states of the
  two sets around it with W = (1 + cos π(et − t₁)/(t₂ − t₁)) / 2 (with W′ in the velocity), in
  AFSPC mode with each epoch as a UTC Julian date in one double (as SPICE: 40 µs), turned from
  TEME of et to J2000 (et used as TT, as SPICE does); outside the epochs the nearest set alone.
- The fit: `count` ≥ 3 TEME positions (m) at minutes from the epoch of the set; Gauss–Newton
  in unknowns without the singularities of the classical angles, derivatives by differences
  of the theory; B* optionally; 40 steps and a tolerance of 10⁻⁹ of the RMS by default.
  `myb_tle_from_state` gives the first guess: the osculating elements written as mean ones
  (off by the short-period terms of the oblateness, some kilometres).

## Domain of validity

- `myb_tle_parse`: lines of at least 68 characters; a digit in column 69 must equal the
  modulo-10 checksum (digits, `-` = 1) unless `MYB_TLE_NO_CHECKSUM`; `MYB_E_FORMAT` with the
  line and column in `diag` (a wrong checksum, lines swapped, different catalogue numbers).
- `myb_tle_format`: `MYB_E_DOMAIN` for a value that does not fit its field — an eccentricity
  outside [0, 1) (also one that rounds to 1), an angle outside [0, 1000), |ṅ/2| of 1 or more, a
  B* whose decimal exponent has two digits, a number above the width of its field (a
  catalogue number above 339 999, a revolution above 99 999), values not finite.
- `myb_sgp4_init`: `MYB_E_DOMAIN` for an eccentricity outside [0, 1), a mean motion ≤ 0, a
  perigee below the surface. `myb_sgp4_propagate`: the codes of the theory — 1 (mean
  eccentricity outside [−0.001, 1)), 2 (mean motion ≤ 0), 3 (perturbed eccentricity outside
  [0, 1]), 4 (semi-latus rectum < 0) return `MYB_E_DOMAIN`; 6 (below the Earth's radius)
  returns `MYB_OK` with the state.
- `myb_tle_from_omm`: `MYB_E_DOMAIN` for a message without the mean motion or the TLE
  parameters, or a catalogue number above 339 999.
- `myb_tle_fit`: `MYB_E_DOMAIN` (fewer than three positions, invalid options, elements that
  are no orbit of the theory, unknowns the positions do not determine), `MYB_E_NOCONV` (the
  budget; the best set so far is kept). `myb_tle_from_state`: `MYB_E_DOMAIN` if not an
  ellipse.
- `myb_sgp4_spk10_state`: `MYB_E_FORMAT` (not type 10, inconsistent layout), `MYB_E_RANGE`
  (outside the UTC table). `myb_sgp4_spk10_words`: `MYB_E_DOMAIN` (no sets, epochs not
  increasing), `MYB_E_OVERFLOW` (buffer too small).

## Accuracy

| What | Reference | Result |
|---|---|---|
| 33 sets of SGP4-VER.TLE, 700 lines (WGS-72, improved) | `tcppver.out` of Vallado et al. (positions to 10⁻⁸ km) | **5.0·10⁻⁹ km** — half of the last printed digit |
| 3.5 years of the deep-space SL-12 R/B (e = 0.79) | the same | 1.15·10⁻⁷ km, as python-sgp4 itself |
| Error codes along the file | [1, 1, 6, 6, 4, 3, 6] | the same |
| SPK type 10, near-Earth / deep space | CSPICE N0067 (`spkw10_c`, `spkgeo_c`) | **2.7·10⁻¹² km** / 1.9·10⁻⁹ km; velocities 9·10⁻¹⁴ / 10⁻¹² km/s |
| SGP4 in AFSPC mode, near-Earth / deep space | `evsgp4_c` | 4·10⁻¹¹ / 4.4·10⁻⁹ km |
| SPK type 10 written: 3 segments, 235 sets (a directory of epochs at 230) | the words of `SPKW10` | the same words; the copies read as the files, bit for bit |
| The packet of a two-line set | `getelm_c` | elements within 10⁻¹⁴ (relative); the epoch within 19 µs |
| 33 sets written back | the same values | bit for bit; 31 character for character |
| An OMM turned into a set and back | the same message | the same lines; the same propagation, bit for bit |
| Fit to the theory's own positions: low orbit 2 days with B*; eccentric 3 days; 12-hour (deep space) 4 days | the set that made them | elements to 10⁻¹⁵, B* to six digits, residuals 5·10⁻⁸ m (8 steps from the osculating elements) |

In deep space the epoch matters: SPICE converts TDB → UTC with its one-term TDB − TT, myb with
the Fairhead–Bretagnon series; both round the Julian date to 40 µs. The reference's own epoch
is computed through a Julian date of ~2.45·10⁶ (rounding 40 µs), which gives up to 4 mm at
t = 0 for very elongated orbits; the port reproduces that path exactly.

Test: `test_sgp4` ([tests/ext/test_sgp4.c](../../tests/ext/test_sgp4.c)).

## Cost

Not measured.

## Limitations

- The theory itself is the limit: about a kilometre at the epoch, kilometres a day after.
  It is specific to the Earth (WGS-72 constants, lunar and solar terms with fixed
  coefficients).
- The library's own SPK reader does not evaluate type 10; only `myb_sgp4_spk10_state` does.
  The nutation words of a packet are not read (as SPICE since N0066).
- Two sets of the verification file differ from the catalogue's notation (no ephemeris type;
  a zero B* written as `00000-0`), so they come back as the same values but not the same
  characters.
- The fit uses derivatives by differences of the theory; it was checked on positions of the
  theory itself, not on an independent ephemeris.

## Data

- None for the theory: the constants are in the code; the time kernels
  (`data/time/*.tk`) and `data/frames/iau1980.tk` for SPK type 10.
- Tests: `tests/data/sgp4/SGP4-VER.TLE`, `tcppver.out` (Vallado et al., as distributed with
  python-sgp4); `tests/data/spk/t10.bsp`, `t10_many.bsp` (written by CSPICE `spkw10_c`) and
  `t10_reference.txt` (CSPICE outputs).
