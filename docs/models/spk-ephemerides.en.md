# SPK ephemerides and SP3 orbits

**English** | [Русский](spk-ephemerides.ru.md) · [All passports](index.en.md)

## What it is

The library's own reader and writer of NAIF SPK files — the planetary ephemerides DE440/DE441
of JPL and EPM of IAA RAS, small-body files of JPL Horizons, spacecraft trajectories. A segment
gives the state of a target relative to a centre over an interval in a frame named by a NAIF
code; the reader joins two bodies through a chain of covering segments and sums the branches
with compensated summation, and with a graph of frames it takes segments in any frames. Files
are memory-mapped and read on demand. The writer lays out segments of types 9, 13, 2 and 3
from states or a provider, and copies segments of any type word for word. SP3-c/d precise
orbits (IGS, ILRS) are read and written in the same header.

## Source

- NAIF "SPK Required Reading" and "DAF Required Reading" (JPL): the file layout, the segment
  types, the priority rule (later file, later segment). No SPICE code is used; CSPICE N0067
  (through spiceypy) served only as an external oracle that wrote the test files and their
  reference states.
- Type 15: two-body motion from the periapsis, then the node and the periapsis turned in
  proportion to the true anomaly passed, as in SPICE (the J₂ precession of the type).
- SP3-c and SP3-d: the format descriptions of S. Hilla (2010, 2016).
- DE440: R. S. Park et al. (2021); the time ephemeris of `de440t.bsp` is used by
  [time scales](time-scales.en.md).

## Interface

Header [`myb/kin.h`](../../include/myb/kin.h):

- reading: `myb_spk`, `myb_spk_create`, `myb_spk_destroy`, `myb_spk_load_file`,
  `myb_spk_load_memory`, `myb_spk_segment_info`, `myb_spk_segment_count`, `myb_spk_segment`,
  `myb_spk_segment_size`, `myb_spk_segment_data`, `myb_spk_state`, `myb_spk_state_in`;
- as a provider: `myb_state_function`, `myb_spk_body` (its member `graph` for SPK across
  frames), `myb_spk_body_state`; the time ephemeris `myb_spk_time`,
  `myb_spk_time_difference`;
- writing: `myb_spk_writer`, `myb_spk_writer_create`, `myb_spk_writer_destroy`,
  `myb_spk_segment_header`, `myb_spk_writer_add_states` (types 9, 13),
  `myb_spk_writer_add_chebyshev` (types 2, 3), `myb_spk_writer_add_raw`,
  `myb_spk_writer_size`, `myb_spk_writer_write`, `myb_spk_writer_write_file`;
- SP3: `myb_sp3`, `myb_sp3_load`, `myb_sp3_load_file`, `myb_sp3_destroy`,
  `myb_sp3_satellite_count`, `myb_sp3_satellite`, `myb_sp3_epoch_count`,
  `myb_sp3_time_system`, `myb_sp3_frame`, `myb_sp3_epoch`, `myb_sp3_state`, `myb_sp3_clock`,
  `myb_sp3_sigma`, `myb_sp3_accuracy`, `myb_sp3_flags` (`MYB_SP3_CLOCK_EVENT`,
  `MYB_SP3_CLOCK_PREDICTED`, `MYB_SP3_MANEUVER`, `MYB_SP3_ORBIT_PREDICTED`), `myb_sp3_data`,
  `myb_sp3_write`.

The graph of frames (`myb_frame_graph`) is in [frames and geodesy](frames-geodesy.en.md);
SPK type 10 (two-line elements) is evaluated in `ext/sgp4` — see [SGP4 and TLE](sgp4-tle.en.md).

## Inputs and outputs

- Time is the ephemeris' argument as a `myb_time` on its scale (TDB for the JPL and IAA
  ephemerides): its two parts keep t − t_mid exact in the Chebyshev and interpolation
  arguments, which one double of seconds (0.1 µs resolution today) does not.
- States in m and m/s (the files hold km and km/s). `myb_spk_state` reports the frame code
  shared by all the segments used; `target == center` gives zeros.
- `myb_spk_state_in`: every segment's state is turned into the frame asked for at t —
  velocities with the rotation of the frames — before the sum; t on the scale the graph was
  created with.
- The identifiers are the file's: the library attaches no meaning to them.
- Writer: states in m and m/s, epochs in seconds past J2000 TDB; one division by 1000 on the
  way in and one multiplication on the way out, so a state read back agrees to an ulp. Types 9
  (Lagrange, degree + 1 states) and 13 (Hermite of odd degree, (degree + 1)/2 states with
  velocities), degree 1–27; types 2 (positions) and 3 (positions and velocities): Chebyshev
  polynomials through a provider at the degree + 1 Chebyshev nodes of equal parts, with an
  optional report of the largest deviation half-way between the nodes (degree + 2 more
  provider calls a part). Raw segments are written as given, unchecked. The output is one
  little-endian DAF file, into a buffer or a file.
- SP3: positions in m, velocities in m/s, clocks in s (the file keeps µs) and s/s; epochs are
  calendar labels in the file's time system, which the caller converts. A position 0, 0, 0 is
  reported missing. The writer keeps six decimals of km and dm/s (1 mm, 0.1 µm/s) and of µs
  (1 ps, |clock| < 1 s). Standard deviations are powers of the header's bases (1.25, 1.025): a
  value read is within 12 % (1.2 % for clocks) of the producer's.

## Domain of validity

- Each file covers its own interval (DE440: 1550–2650); a query outside every covering chain
  gives `MYB_E_RANGE`.
- `myb_spk_state`: `MYB_E_FORMAT` for segments in different frames (without a graph), a type
  that is loaded but not evaluated, or invalid data; `MYB_E_DOMAIN` when the two-body
  propagation of types 5 and 15 fails.
- Loading is atomic: `MYB_E_IO`, `MYB_E_FORMAT` (not a DAF/SPK file, a corrupt FTP string, an
  inconsistent segment — lengths and monotonic epochs are checked at load), `MYB_E_NOMEM`,
  `MYB_E_OVERFLOW`. Little- and big-endian files are read on any host.
- Writer: `MYB_E_DOMAIN` (the header, the type, the degree, too few states, epochs not
  increasing, values not finite, stop ≤ start); `MYB_E_OVERFLOW` beyond 2³¹ words.
- SP3: `MYB_E_FORMAT` (diag gives the line); the writer refuses an invalid description with
  `MYB_E_DOMAIN` and a small buffer with `MYB_E_OVERFLOW`; SP3-c holds at most 85 satellites.

## Accuracy

| What | Reference | Result |
|---|---|---|
| ~7000 states of types 1, 2, 3, 5, 8, 9, 12, 13, 14, 17–21 and chains | CSPICE N0067 on files it wrote | positions ≤ 6.3·10⁻¹⁶, velocities ≤ 10⁻¹⁴, Hermite velocities 1.1·10⁻¹⁴ (relative) |
| Hermite velocities (types 12, 13, 19) | the interpolant in 300 bits | ≤ 1.6·10⁻¹⁶ (CSPICE: up to 2.3·10⁻¹⁴, it interpolates on absolute epochs) |
| Type 15, 440 states (circle to hyperbola) | its rule in 60 digits | position 9·10⁻¹⁶, velocity 1.5·10⁻¹⁵; arcs of 3400 revolutions 1.3·10⁻¹⁴ (CSPICE 1.6·10⁻¹¹) |
| DE440, 11 354 positions and velocities of 1550–2650 | JPL `testpo.440` | 7.1·10⁻¹⁵ au (≈1 mm), 2·10⁻¹⁷ au/d; JPL tolerance 10⁻¹³ |
| 433 Eros, type 21 file from Horizons | Horizons vectors | 1.8·10⁻¹² (0.3 m) — CSPICE gives the same: Horizons prints its own integration |
| SPK across frames, 336 states (ECLIPJ2000, B1950, GALACTIC, FK4, J2000, a PCK frame) | CSPICE `spkgeo_c` | 1.3·10⁻¹⁵ (velocities 3·10⁻¹⁴ — Hermite segments, CSPICE's difference) |
| Written file of 34 segments | CSPICE reading it | positions 3.9·10⁻¹⁶, velocities 1.1·10⁻¹³ |
| Types 9, 13 at their epochs | the states written | within 1.8 ulp of the vector's length |
| Type 2, degree 12, 1000 s | the two-body orbit | 6·10⁻⁸ m |
| 44 raw copies of 14 types | the originals | bit for bit, in the library and in CSPICE (308 states) |
| ILRS SP3 of LAGEOS-1 read and written | the file | 258 720 characters of epoch and state lines — character for character |
| SP3 with clocks, accuracy codes, flags | the file | all 37 lines the same |

Tests: [`test_spk`](../../tests/unit/test_spk.c),
[`test_spk_writer`](../../tests/unit/test_spk_writer.c),
[`test_testpo`](../../tests/unit/test_testpo.c) (skipped without the downloads),
[`test_frame_graph`](../../tests/unit/test_frame_graph.c),
[`test_sp3`](../../tests/unit/test_sp3.c).

## Cost

Not measured.

## Limitations

- Types read: 1, 2, 3, 5, 8, 9, 12, 13, 14, 15, 17, 18, 19, 20, 21. Others load, and a query
  that needs one of their segments fails with `MYB_E_FORMAT`; type 10 is evaluated outside the
  core (`ext/sgp4`) from the raw words.
- Written natively: types 9, 13, 2, 3. Types 1 and 21 (the difference arrays of JPL's
  integrator, which only it produces), 10, 5, 8, 12, 14, 17–20 only as raw copies, whose
  words are neither converted nor checked.
- Quirks of SPICE found and followed: its writer stores the type 19 selection flag as −1, not
  0 as documented; the type 21 difference line has MAXDIM words, a record 4·MAXDIM + 11; the
  type 17 velocity is the full time derivative of the position, derived analytically.
- The nutations of DE440 (target 14) are not in its SPK conversion.
- SP3 epochs are labels: converting them (e.g. `myb_timesys_from_utc`) is the caller's.
- `myb_spk_load_file` must not run while another thread queries the set; queries themselves
  are thread-safe.

## Data

Downloaded (`cmake/FetchData.cmake`, SHA-256 checked): `de440.bsp` (JPL DE440, NAIF SPK,
1550–2650), `testpo.440` and `header.440` (JPL SSD); with `-DLARGE=ON` `de440t.bsp` (152 MB)
and `sb441-n373s.bsp` (the 373 small bodies that perturb DE440, 982 MB). JPL ephemerides are
free and not redistributed. In the repository: the test files `tests/data/spk/` (written with
CSPICE as an oracle; `eros.bsp` from JPL Horizons) and the ILRS LAGEOS-1 SP3 excerpt
`tests/data/sp3/` (freely available with acknowledgement of the ILRS).
