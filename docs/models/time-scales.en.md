# Time scales

**English** | [Русский](time-scales.ru.md) · [All passports](index.en.md)

## What it is

Instants and the conversions between time scales. An instant is whole seconds (int64) plus a
fraction in [0, 1), counted from J2000.0 of its own scale. The scales are not built in: a time
system is a tree built from kernel variables — TAI at the root, TT, TCG, TDB, TCB, GPS, GST,
BDT defined by their relations to a parent (an offset, a rate, or a series). UTC is not a
scale of the tree but calendar labels of TAI read through a leap-second table with an expiry
date. TDB − TT comes from the Fairhead–Bretagnon series by default; it can be given by a
function instead (the time ephemeris of DE440t), or tabulated for speed. An environment of a
force model keeps its own table of the scales its providers use.

## Source

- Definitions: IAU 1991 Resolution A4, IAU 2000 Resolution B1.9 (TT, TCG, L_G), IAU 2006
  Resolution B3 (TDB, TCB, L_B, TDB0); IS-GPS-200 for GPS; the GNSS times as nominal offsets
  from TAI (`data/time/time-scales.tk`).
- TDB − TT: L. Fairhead, P. Bretagnon, *A&A* 229, 240 (1990), the full 787 terms as distributed
  with IAU SOFA (`iau_DTDB`) plus the five JPL-mass adjustments of `eraDtdb`; the numbers are
  copied verbatim from ERFA `src/dtdb.c` (BSD-3) into `data/time/tdb-tt-fb1990.tk`.
- UTC: before 1972 the USNO `tai-utc.dat` rows as tabulated in ERFA `dat.c` (BSD-3); from 1972
  the IERS `Leap_Second.dat` (Bulletin C 72, July 2026), every row checked against ERFA.
- TDB − TT from an ephemeris: the JPL file `de440t.bsp`, which holds TT − TDB at the geocentre
  as the first coordinate of "body" 1000000001 relative to 1000000000.
- The library's code is its own; ERFA is compiled only into the tests as the SOFA reference.

## Interface

Header [`myb/time.h`](../../include/myb/time.h):

- instants: `myb_time`, `myb_time_from_seconds`, `myb_time_from_dd`, `myb_time_add`,
  `myb_time_add_dd`, `myb_time_diff`, `myb_time_seconds`, `myb_time_compare`,
  `myb_time_from_jd`, `myb_time_to_jd`;
- calendar: `myb_calendar`, `myb_calendar_to_time`, `myb_time_to_calendar`;
- time systems: `myb_timesys`, `myb_scale`, `myb_timesys_create`, `myb_timesys_destroy`,
  `myb_timesys_scale`, `myb_timesys_scale_name`, `myb_timesys_convert`,
  `myb_timesys_from_utc`, `myb_timesys_to_utc`;
- the series relation replaced or tabulated: `myb_scale_difference`,
  `myb_timesys_set_difference`, `myb_timesys_tabulate`.

TDB − TT from an ephemeris: `myb_spk_time`, `myb_spk_time_difference` of
[`myb/kin.h`](../../include/myb/kin.h). The table of scales of an environment: the members
`timesys`, `scale`, `table_step` of `myb_env_options` and `myb_env_time` of
[`myb/env.h`](../../include/myb/env.h). The kernels are read by `myb_kpool` of
[`myb/io.h`](../../include/myb/io.h) (see [numerics and formats](numerics-formats.en.md)).

## Inputs and outputs

- An instant: `s` whole seconds and `f` ∈ [0, 1) since J2000.0 (JD 2451545.0) of its scale;
  resolution better than 10⁻¹⁶ s over ±2.9·10¹¹ years. Differences come as double-double
  (`myb_time_diff`, relative error < 2⁻¹⁰³).
- Julian dates are two-part: on output jd1 = 2451545.0 + whole days, 0 ≤ jd2 < 1.
- Calendars are proleptic Gregorian with astronomical years (year 0 = 1 BC); the second is in
  [0, 60) on a continuous scale, [0, 60 + leap) in a UTC label.
- Relation kinds of the kernel: OFFSET (X = FROM + OFFSET), RATE and RATE_INVERSE (an affine
  relation with a rate about T0 given as an exact two-part JD), SERIES (X = FROM + S(X),
  solved by fixed-point iteration).
- `myb_timesys_convert`: offsets add no error beyond the final rounding of the fraction, rates
  are evaluated in double-word arithmetic, series are iterated to below 10⁻¹⁵ s.
- A function given by `myb_timesys_set_difference` returns scale − parent in seconds at an
  instant of the scale; `myb_spk_time_difference` returns TDB − TT at an instant of TDB.
- `myb_timesys_tabulate`: nodes every `step` s (0: an hour) between t0 and t1 of the scale, a
  cubic through the four nodes around the instant; outside the table the relation as before;
  values depend on the nodes only.
- An environment: model time t in seconds from `epoch` on `scale`; every provider gets the
  same instant on the scale it was registered with, read from a table every `table_step` s
  (0: 3600), Lagrange on four nodes — a pure function of t.

## Domain of validity

- `myb_time` arithmetic: finite arguments with |seconds| < 2⁶²; calendars: |year| ≤ 10⁹
  (`MYB_E_DOMAIN` otherwise, and for an invalid date or a second ≥ 60).
- The UTC table runs from 1960-01-01 to its expiry, 2027-06-28 (exclusive):
  `myb_timesys_from_utc` gives `MYB_E_RANGE` before the first row and from the expiry on,
  `MYB_E_DOMAIN` for an invalid label (23:59:60 on a day without a leap second; the 59.95 s
  last minute of 1961-07-31), `MYB_E_NOT_FOUND` without a UTC table.
- The series is valid to a few ns over 1950–2050 (ERFA note: ±3 ns against time ephemerides
  integrated with DE405); the tests compare it with ERFA over 1550–2450.
- DE440t: outside the file `myb_spk_time_difference` gives `MYB_E_RANGE`, and so does a
  conversion through it.
- `myb_timesys_create`: `MYB_E_NOT_FOUND` (missing variable), `MYB_E_FORMAT` (a cycle, an
  unknown kind, an unordered UTC table), `MYB_E_NOMEM`; diag names the variable.
- `myb_timesys_set_difference`, `myb_timesys_tabulate`: `MYB_E_DOMAIN` for a scale whose
  relation is not a series, t1 not after t0, a negative or non-finite step, more than 2²⁶
  nodes.
- Environment: at most `MYB_ENV_MAX_SCALES` = 4 scales besides the model's; `MYB_E_DOMAIN`
  for a step negative, not finite or above 10¹⁵ s; `MYB_E_RANGE` for |t| beyond 10¹⁷ s or
  2⁵² table steps.

## Accuracy

| What | Reference | Result |
|---|---|---|
| Calendar → instant, 200 000 dates of −4700 … 3000 | `eraCal2jd` | exact (tolerance 0) |
| TDB − TT, 20 000 instants of 1550 … 2450 | `eraDtdb` without topocentric terms | ≤ 10⁻¹⁴ s (tolerance; ERFA uses the platform's `sin`) |
| TT ↔ TCG, TDB ↔ TCB, TAI → TT, 50 000 instants each, ±300 years | `eraTttcg`, `eraTcgtt`, `eraTdbtcb`, `eraTcbtdb`, `eraTaitt` | ≤ 3.2·10⁻¹¹ s over all ERFA comparisons (the resolution of ERFA's two-part JD); tolerance 10⁻¹⁰ s |
| UTC → TAI, 100 000 labels of 1960–2026 and 1500 inside leap seconds and pre-1972 steps | `eraUtctai` | within the 10⁻⁹ s tolerance |
| Every pair of the 8 scales there and back | itself | < 10⁻¹⁵ s |
| TAI → UTC → TAI, 200 000 instants of 1960–2027 | itself | < 10⁻¹² s; leap-second labels exact |
| Table of TDB − TT, step of an hour | the series | 4·10⁻¹⁶ s |
| Rate of TT − TDB in DE440t, 400 dates of 1952–2048 | L_C − (v²/2 + U)/c² from the same ephemeris | 1.4·10⁻¹⁶ (the rate itself 3.4·10⁻¹⁰) |
| TDB − TT of DE440t | Fairhead–Bretagnon series | within 39 ns (15 ns on average) |
| DE440t as the relation of TDB, there and back | the file | < 10⁻¹⁵ s |
| Environment table (step 1 h), ±40 days | exact conversion | TDB − TT 3.5·10⁻¹⁶ s; TCB − TT (20 s with a rate) 3.1·10⁻¹⁵ s |

Tests: [`test_erfa_time`](../../tests/erfa/test_erfa_time.c) (ERFA),
[`test_time`](../../tests/unit/test_time.c) (exact properties, UTC, tables),
[`test_de440_time`](../../tests/unit/test_de440_time.c) (DE440t; skipped without the file),
[`test_env`](../../tests/unit/test_env.c) (table of an environment; a node equals the exact
value to 2.5·10⁻¹⁶ s, the same bits whatever was asked before),
[`test_kpool`](../../tests/unit/test_kpool.c) (the kernel reader).

## Cost

- A TT → TDB conversion evaluates the 787-term series twice (after the first evaluation the
  error is < 2·10⁻¹² s, after the second < 2·10⁻²¹ s); it takes 92 µs (phase 20).
- Through the hourly table a conversion takes a fraction of a microsecond (`test_time` prints
  both); an environment evaluates the series once an hour of model time.

## Limitations

- The scales are run-time values, not C types: mixing two scales is caught only by an explicit
  conversion through the tree, not by the compiler.
- The series is a model of the relativistic TCB − TCG relation good to a few ns; only the
  geocentric part is given. Topocentric terms belong to the observer's proper time — see
  [1PN relativity](relativity-1pn.en.md).
- The GNSS times are nominal (TAI − 19 s, TAI − 33 s): their realisations differ by
  nanoseconds, which the library does not model.
- UT1 is not a scale of the tree; it comes from EOP — see
  [Earth orientation](earth-orientation.en.md).
- The leap-second table is a file: a new leap second needs a new row and a new expiry date;
  until then conversions past the expiry fail rather than guess.
- CSPICE's TDB − TT has a single term: epochs compared with SPICE differ on that account (the
  TLE epochs of phase 15 agree with `getelm` within 19 µs).
- `myb_timesys_set_difference` and `myb_timesys_tabulate` must not be called while another
  thread converts.

## Data

In the repository (NAIF text-kernel format):

- `data/time/time-scales.tk` — the tree of TAI, TT, TCG, TDB, TCB, GPS, GST, BDT (IAU and IERS
  constants; definitions are facts).
- `data/time/tdb-tt-fb1990.tk` — the 787 + 5 terms of TDB − TT (from ERFA `dtdb.c`, BSD-3).
- `data/time/utc-tai.tk` — TAI − UTC since 1960, expiry MJD 61584 = 2027-06-28 (ERFA `dat.c`,
  BSD-3; IERS `Leap_Second.dat`).

Optional download (`-DLARGE=ON`): `de440t.bsp` (JPL, 152 MB, SHA-256 in
`cmake/FetchData.cmake`) for TDB − TT from the ephemeris.
