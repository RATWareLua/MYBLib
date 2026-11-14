# Thermosphere: NRLMSISE-00, space weather and the common input

**English** | [Русский](thermosphere-models.ru.md) · [All passports](index.en.md)

## What it is

The empirical model NRLMSISE-00 of the Earth's neutral atmosphere (number densities of nine
species, the total mass density, the exospheric and the local temperature, from the ground to
the exosphere) as an optional module, with an adapter that makes it the atmosphere of a force
model; a reader of the space-weather table (F10.7, Kp, ap); and the common input of the
upper-atmosphere models — the UTC day and time, the geodetic place, the local solar time and the
indices at an instant and a body-fixed place. Through that input the user's own copies of
DTM2020, NRLMSIS 2.1 and HWM14 plug in as a `myb_density_function` or a `myb_wind_function`.
**These three models are not in the repository**: the library supplies only their input.
JB2008 has a passport of its own: [JB2008](jb2008.en.md).

## Source

- NRLMSISE-00: J. M. Picone, A. E. Hedin, D. P. Drob, *J. Geophys. Res.* 107, 1468 (2002).
  The code is a port of the C release 20041227 by Dominik Brodowski (revision
  `bc9a2feb…`, 2019-07-09), placed by him in the public domain; the model, its algorithms and
  its coefficient tables are U.S. Government material not subject to copyright (the notice
  under 17 U.S.C. 403 is in [`ext/nrlmsise00/NOTICE`](../../ext/nrlmsise00/NOTICE)). The
  changes: the shared static variables became a state passed to the functions (thread
  safety), the tables are `static const`, exp, log, sin, cos are myb's correctly rounded
  functions, `pow` is exact repeated squaring for integer exponents and `sqrt` for ½, `ghp7`
  (the height of a given pressure) is left out; the computations are otherwise unchanged.
  myb's changes are under BSD-3; the module is not endorsed by, nor affiliated with, the NRL.
- The space-weather table: the CSV layout of CelesTrak (`SW-All.csv`, `SW-Last5Years.csv`;
  T. S. Kelso, "Space Weather Data Documentation"), read by the names of its columns.
- DTM2020, NRLMSIS 2.1, HWM14: their licences forbid porting
  ([PLAN §9](../PLAN.en.md#9-licenses), checked on 2026-11-11). NRLMSIS 2.1 — the NRL "Open
  Source Academic Research License" (derived works go back to NRL, commercial use is
  separate); DTM2020 — CNES, academic non-commercial use only, no modification and no
  redistribution; HWM14 — no licence at all. By the user's decision (2026-11-12) the three
  were translated only on the user's machine, outside git, for a local check; their code is
  not in the repository.

## Interface

Header [`myb_nrlmsise00.h`](../../ext/nrlmsise00/myb_nrlmsise00.h) (target
`myb::nrlmsise00`, not part of libmyb):

- functions `myb_nrlmsise00`, `myb_nrlmsise00_density`;
- structures `myb_nrlmsise00_input`, `myb_nrlmsise00_output`, `myb_nrlmsise00_atmosphere`.

Header [`myb/weather.h`](../../include/myb/weather.h) (in libmyb):

- functions `myb_space_weather_load`, `myb_space_weather_load_file`,
  `myb_space_weather_destroy`, `myb_space_weather_days`, `myb_space_weather_day_at`,
  `myb_space_weather_inputs_at`, `myb_thermosphere_point_at`;
- structures `myb_space_weather_day`, `myb_space_weather_inputs`, `myb_thermosphere_setting`,
  `myb_thermosphere_point`; the marks `MYB_SPACE_WEATHER_OBSERVED`, `_INTERPOLATED`,
  `_PREDICTED`;
- the callbacks `myb_density_function` and `myb_wind_function` of
  [`myb/force.h`](../../include/myb/force.h), through which any model enters a force model.

## Inputs and outputs

- `myb_nrlmsise00`: the day of the year (1…366), the seconds of the UT day, the geodetic
  altitude in km, the geodetic latitude and the longitude in degrees, the local solar time in
  hours (the model's convention: seconds / 3600 + longitude / 15), F10.7 of the day before,
  its 81-day mean centred on the day, the daily Ap, and optionally the seven values of the ap
  history (used when switch 9 is −1). 24 switches as in the model (NULL: the standard model);
  switch 0 is ignored — the output is always SI.
- Output: number densities of He, O, N₂, O₂, Ar, H, N and anomalous O in m⁻³ (indices 0…4,
  6…8), the total mass density in kg m⁻³ (index 5), the exospheric temperature and the
  temperature at the altitude in K. `anomalous_oxygen ≠ 0` (GTD7D) counts the anomalous
  oxygen in the mass density, as appropriate for drag above 500 km.
- `myb_nrlmsise00_density`: the body-fixed place r (m) and the instant t on `scale` → the
  density in kg m⁻³. The date and the time of day are those of UTC, the height and the
  latitude geodetic on the given ellipsoid, the local solar time the model's own. The indices
  come from the weather table (the observed flux; switch 9 = −1 with the ap history where the
  table has the 57 hours) or are the three fixed numbers.
- The table: F10.7 in solar flux units (10⁻²² W m⁻² Hz⁻¹), observed at the Earth and adjusted
  to 1 au, with the 81-day means centred on the day and ending with it; Kp (0…9; the file
  gives tenths) and ap of the eight three-hour spans, the daily Ap, the kind of the flux
  (observed, interpolated, predicted). Days are UTC days (MJD); a value the file does not give
  is a NaN.
- `myb_space_weather_inputs_at` (UTC day and seconds): F10.7 of the day before, its centred
  81-day mean, the daily Ap, the seven-value ap history of A. E. Hedin's models (daily Ap; ap
  of the span and of 3, 6, 9 hours before; the means of the eight values 12–33 and 36–57 hours
  before), Kp of the span. `adjusted ≠ 0` takes the flux adjusted to 1 au (the MSIS models
  were built on the observed one). Where a value of the 57 hours is missing, `has_history = 0`
  and the history holds the daily Ap seven times.
- `myb_thermosphere_point_at`: MJD, day of the year, UTC seconds (below 86400; a leap second
  is the day's end), geodetic height in **m**, latitude and longitude in **rad**, local solar
  time in hours (time of day plus longitude), the indices. The NRLMSISE-00 adapter converts to
  the km and degrees of the model.

## Domain of validity

- NRLMSISE-00 is checked from 0 to 1000 km. Refused with `MYB_E_DOMAIN`: input that is not
  finite, an altitude below 0 or above 20 000 km, a day outside 1…366, a latitude outside
  [−90°, 90°], negative indices, switch 9 = −1 without an ap history.
- The adapter and `myb_thermosphere_point_at` return the errors of `myb_timesys_to_utc` and
  `myb_geodetic_from_cartesian`, and `MYB_E_RANGE` for an instant outside the table; the
  adapter refuses a place below the ellipsoid with `MYB_E_DOMAIN`.
- The reader: `MYB_E_FORMAT` (no DATE column, a date or a number that is not one, days out of
  order; `diag` gives the line), `MYB_E_NOMEM`, `MYB_E_OVERFLOW`, `MYB_E_IO`. A day outside the
  table or in a gap — `MYB_E_RANGE`; so is an instant whose day before is missing or whose
  flux or daily Ap the table does not give (the monthly predictions give the flux alone).
  Seconds outside 0…86401 — `MYB_E_DOMAIN`.

## Accuracy

| What | Reference | Result |
|---|---|---|
| NRLMSISE-00, 17 reference cases, 0–1000 km | section 7 of the C release's DOCUMENTATION (the NRL test output), 7 digits | 4.8·10⁻⁷ — the printed precision (test bound 6·10⁻⁷) |
| GTD7D: the anomalous oxygen in the mass density | 16 amu × its number density | within 10⁻³ of the density |
| The adapter, 420 km over 51.6° N, 2026-03-05 13:00 UTC, table indices with history; fixed indices with anomalous O | `myb_nrlmsise00` called by hand with the same inputs | 10⁻⁹ |
| The table: values of a day, the lags, the ap history | an excerpt in the CelesTrak layout made for the test | exact; the means of the history to 10⁻¹⁵ |
| The thermosphere point at the same instant and place | day 64, 46 800 s, 420 km, 51.6° N, 30° E, 15 h | 10⁻⁶ s, 10⁻⁶ m, 10⁻¹² rad, 10⁻⁹ h |
| DTM2020, NRLMSIS 2.1, HWM14 (the user's local copies) | the original codes | bit for bit in double precision with a correctly rounded libm — checked on the user's machine only, no test in the repository |

Tests: `test_nrlmsise00` ([tests/ext/test_nrlmsise00.c](../../tests/ext/test_nrlmsise00.c)),
`test_weather` ([tests/unit/test_weather.c](../../tests/unit/test_weather.c)).

## Cost

Not measured.

## Limitations

- Only NRLMSISE-00 (and [JB2008](jb2008.en.md) in its own module) is in the repository.
  DTM2020 and NRLMSIS 2.1 cannot be reproduced without their coefficient files, which share
  their licences; HWM14 has no licence. The user plugs a copy in through
  `myb_thermosphere_point_at`.
- NRLMSISE-00 is specific to the Earth; `ghp7` (pressure levels) is not available.
- The local solar time is the model's convention (time of day plus longitude), not the true
  solar time; the adapter always takes the observed flux.
- The reader was checked on an excerpt made in the CelesTrak layout, not on a file of
  CelesTrak: the site was unreachable from the development machine.
- No wind model is in the repository; HWM14 enters only as the user's `myb_wind_function`.

## Data

- None for NRLMSISE-00: the coefficient tables are constants of the code.
- The adapter and the point need a time system with UTC: `data/time/time-scales.tk`,
  `tdb-tt-fb1990.tk`, `utc-tai.tk` (in the repository).
- `SW-Last5Years.csv` (CelesTrak, a living series without a SHA-256) is in the second list of
  `cmake/FetchData.cmake` (`-DLARGE=ON` or `-DONLY=SW-Last5Years.csv`); no test uses it — the
  tests read `tests/data/weather/sw-excerpt.csv`, 15 days made for them. The repository states
  no licence for the CelesTrak files.
