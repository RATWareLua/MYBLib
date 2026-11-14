# Binary PCK, CK and spacecraft clocks

**English** | [Русский](pck-ck-sclk.ru.md) · [All passports](index.en.md)

## What it is

Orientation read from NAIF binary files. A binary PCK gives a body-fixed frame relative to an
inertial base frame as three Euler angles — the high-precision Earth (ITRF93) and the lunar
principal-axis frame of the JPL ephemerides; the library reads and writes it. A CK gives the
rotation of a spacecraft or instrument frame (the C-matrix) as quaternions at instants of the
spacecraft clock, optionally with angular velocities; segment types 1–6 are read and written,
and which pointing a query gets follows SPICE's `ckgpav_c` type by type. An SCLK of type 1
ties the clock's ticks to a time scale. All three enter the frame graph and an environment as
orientation providers.

## Source

- NAIF "PCK Required Reading", "CK Required Reading", "SCLK Required Reading", "DAF Required
  Reading" (JPL). No SPICE code is used: CSPICE N0067 (through spiceypy) was the external
  oracle that wrote test files (`ckw05_c` among them) and their references.
- The DAF layout module is shared with the SPK writer ([SPK ephemerides](spk-ephemerides.en.md)).
- The lunar frame: `moon_pa_de440_200625.bpc` of DE440; the Earth: NAIF's
  `earth_latest_high_prec.bpc` (ITRF93 from current EOP).

## Interface

Header [`myb/frames.h`](../../include/myb/frames.h), binary PCK part:

- reading: `myb_pck`, `myb_pck_create`, `myb_pck_destroy`, `myb_pck_load_file`,
  `myb_pck_load_memory`, `myb_pck_euler`, `myb_pck_rotation`;
- writing: `myb_pck_writer`, `myb_pck_writer_create`, `myb_pck_writer_destroy`,
  `myb_pck_writer_add_chebyshev` (type 2), `myb_pck_writer_add_raw`, `myb_pck_writer_size`,
  `myb_pck_writer_write`, `myb_pck_writer_write_file`;
- provider: `myb_pck_orientation`, `myb_orientation_pck`.

Header [`myb/ck.h`](../../include/myb/ck.h):

- clocks: `myb_sclk`, `MYB_SCLK_MAX_FIELDS`, `myb_sclk_create`, `myb_sclk_destroy`,
  `myb_sclk_time_system`, `myb_sclk_partitions`, `myb_sclk_to_time`, `myb_sclk_from_time`,
  `myb_sclk_encode`, `myb_sclk_decode`;
- CK reading: `myb_ck`, `myb_ck_create`, `myb_ck_destroy`, `myb_ck_load_file`,
  `myb_ck_load_memory`, `myb_ck_segment_info`, `myb_ck_segment_count`, `myb_ck_segment`,
  `myb_ck_pointing`, `myb_ck_pointing_at`;
- provider: `myb_ck_orientation`, `myb_orientation_ck`;
- CK writing: `myb_ck_writer`, `myb_ck_writer_create`, `myb_ck_writer_destroy`,
  `myb_ck_segment_header`, `myb_ck_writer_add_pointings` (types 1, 3),
  `myb_ck_writer_add_intervals` (2), `myb_ck_writer_add_chebyshev` (4), `myb_ck_packets`,
  `MYB_CK_MAX_WINDOW`, `myb_ck_writer_add_packets` (5), `myb_ck_writer_add_minisegments` (6),
  `myb_ck_writer_size`, `myb_ck_writer_write`, `myb_ck_writer_write_file`.

## Inputs and outputs

- PCK: Euler angles (φ, δ, w) in rad and rates in rad/s with M = R3(w) R1(δ) R3(φ), the
  rotation base frame → body frame; t in seconds past J2000 on the file's scale (TDB for NAIF
  files). Segment types 2, 3 and 20 are read. `myb_pck_orientation.base` is the rotation from
  the inertial axes in use to the segments' base frame (e.g. the TK frame ECLIPJ2000, code 17,
  for ITRF93); NULL when they coincide.
- PCK writer: Chebyshev polynomials of degree 1–27 on equal parts of [start, stop] (TDB), fitted
  to the Euler angles of a provider's rotation at the Chebyshev nodes; the angles are followed
  from node to node; `max_error` (optional, degree + 2 more calls a part) is the largest angle
  between the provider and the polynomials half-way between the nodes.
- CK: quaternions are SPICE's, q = (cos(a/2), sin(a/2) u) for the C-matrix (base → instrument);
  angular velocities of the instrument frame relative to the base, in the base frame, rad/s.
  Times are encoded ticks of the spacecraft clock.
- `myb_ck_pointing_at`: the segments of the instrument in order of precedence (later file,
  later segment); the first with a pointing within `tolerance` ticks wins — inside an interval
  of types 2 and 3 at the instant asked, elsewhere at the nearest instant of the segment. At a
  tie half-way across a gap: types 1–3 take the later record, 4 and 5 the earlier, 6 the earlier
  unless the gap is within the tolerance of the next mini-segment.
- SCLK type 1: fields with moduli and offsets, partitions, a table of (ticks, time, rate);
  ticks ↔ time on the clock's scale (1 = TDB, 2 = TT); strings ↔ ticks.
- `myb_orientation_ck`: the argument is the instant on the clock's scale (register the provider
  with that scale in a graph of frames).
- CK writer: types 1 and 3 from quaternions (made unit) and optional rates; type 2 intervals of
  constant angular velocity; type 4 Chebyshev series (1–19 coefficients a component) on
  intervals; types 5 and 6 packets of four subtypes (Hermite or Lagrange, with or without
  angular velocity) interpolated in sliding windows (even, 2–12 Hermite, 2–24 Lagrange).

## Domain of validity

- PCK: `MYB_E_RANGE` where no segment covers t, `MYB_E_FORMAT` for another type or invalid
  data. The writer refuses (`MYB_E_DOMAIN`) a frame whose pole comes to the z axis of the base,
  where the angles are undefined, or a turn of half a revolution or more between two nodes.
- SCLK: only type 1 (`MYB_E_FORMAT` otherwise); at most 10 fields; `MYB_E_RANGE` before the
  first table row or after the last partition; `myb_sclk_encode` gives `MYB_E_FORMAT` for a
  malformed string and `MYB_E_DOMAIN` for no such partition or a count outside it.
- CK: `MYB_E_RANGE` when no segment gives a pointing, `MYB_E_DOMAIN` for a negative or
  non-finite tolerance, `MYB_E_FORMAT` for invalid data; segment types other than 1–6 load and
  are passed over. The provider gives `MYB_E_FORMAT` for a pointing in another frame than
  `base`.
- Writers: `MYB_E_DOMAIN` (the header, names of 1–40 characters, ticks not increasing, zero
  quaternions, windows, counts, boundaries), `MYB_E_NOMEM`, `MYB_E_OVERFLOW`.

## Accuracy

| What | Reference | Result |
|---|---|---|
| 930 lunar libration values of DE440 from the binary PCK | JPL `testpo.440` | 1 ulp (the angle ψ accumulates to 3.5·10⁴ rad) |
| Site vectors through the ITRF93 binary PCK, 2 sites × 10 dates | JPL Horizons | 0.077 m |
| PCK writer: a uniform rotation | the provider | 7·10⁻¹⁵ |
| PCK writer: a librating frame, degree 5, 9, 13 | the provider | 8·10⁻⁷, 3·10⁻¹¹, 1.5·10⁻¹⁵ rad |
| The written PCK read by CSPICE (`pxform` through a class 2 frame) | the library | 7·10⁻¹⁶ |
| SCLK, two clocks, 519 conversions | CSPICE `sct2e_c`, `sce2c_c`, `scdecd_c`, `scencd_c` | time 6·10⁻⁸ s (the resolution of a double of seconds from J2000); strings and ticks exact |
| CK types 1, 2, 3, 1302 requests | CSPICE `ckgpav_c`, `ckgp_c` | the same 966 found and 336 refused; matrix 7·10⁻¹⁶, angular velocity 9·10⁻¹⁹ rad/s, the instant exact |
| An instrument frame in the graph | CSPICE `sxform_c` | 1.3·10⁻¹⁵ |
| CK types 1, 2, 3 written by the library, 648 requests | CSPICE reading it | the same answers, matrix 3·10⁻¹⁶ |
| CK type 5 written by `ckw05_c` (four subtypes, six segments) | CSPICE | the same 1043 found and 287 refused; matrix 4.4·10⁻¹⁶, angular velocity 6.5·10⁻¹⁷ rad/s |
| CK types 4, 5, 6 written by the library | CSPICE reading it | the same 1953 found and 559 refused; matrix 4.4·10⁻¹⁶ |

Tests: [`test_ck`](../../tests/unit/test_ck.c),
[`test_pck_writer`](../../tests/unit/test_pck_writer.c),
[`test_testpo`](../../tests/unit/test_testpo.c) (librations; skipped without the downloads),
[`test_earth_frames`](../../tests/unit/test_earth_frames.c) (ITRF93; skipped without the
file), [`test_frame_graph`](../../tests/unit/test_frame_graph.c) (a PCK frame in the graph).

## Cost

Not measured.

## Limitations

- Binary PCK types 2, 3, 20 are read; the writer builds type 2 and copies others raw.
- The Euler angles of a written PCK are followed from node to node: the base must keep the
  pole of the frame away from its z axis (for the Earth the ecliptic frame 17 is the usual
  base).
- The ITRF93 PCK is defined relative to ECLIPJ2000, not J2000: the base rotation comes from
  `myb_tkframe` and `data/frames/naif-inertial.tk` (see [frames and geodesy](frames-geodesy.en.md)).
- SCLK of type 1 only. A table whose times are not continuous has no single inverse: the row
  is chosen by the time, as SPICE does.
- Two errors of CSPICE found on the way, where the library deliberately differs: its type 5
  reader keeps the interval of the previous request and can answer from a wrong boundary (the
  library keeps no state); its type 6 reader, from a gap to the next mini-segment, rounds the
  clock rate to a whole number (at a rate of 2.5 the angular velocity is divided by 3).
- The CK provider answers on the clock's scale: a graph or environment on another scale
  converts the instant through the time system.

## Data

Downloaded (`cmake/FetchData.cmake`): `moon_pa_de440_200625.bpc` (NAIF, fixed release,
SHA-256 checked) and `earth_latest_high_prec.bpc` (NAIF, ITRF93 from current EOP, a living
file updated daily). In the repository: `data/frames/naif-inertial.tk` (the SPICE inertial
frames, ECLIPJ2000 among them) and the test files `tests/data/ck/` (`test.bc`, `t5.bc` written
through CSPICE's CK writers, the clock kernel `test.tsc`) with their CSPICE references. NASA/JPL
data products, used unchanged.
