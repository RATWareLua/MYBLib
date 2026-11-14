# Numerics and text formats

**English** | [Русский](numerics-formats.ru.md) · [All passports](index.en.md)

## What it is

The numerical floor under every model and the text formats data come in and results go out.
Elementary functions are correctly rounded (the same bits on every platform), `fma` is the
library's own; double-word ("double-double") arithmetic, error-free transformations and
compensated sums Sum2/Dot2 carry precision where a double is not enough; Carlson's elliptic
integrals, Legendre's forms and Jacobi's functions; a small dense linear algebra with a fixed
order of operations. Numbers are read by a correctly rounding parser and printed by an exact
printer, both independent of the C runtime and the locale. Data files are NAIF text kernels;
tables of numbers are CSV and its relatives; orbits are exchanged as the CCSDS messages OEM,
OPM, OMM and OCM in both notations, KVN and XML.

## Source

- Correctly rounded functions: CORE-MATH (P. Zimmermann, A. Sibidanov et al.; MIT), revision
  `660067dfd3e751c8f07acc31c899ff410aa7206e`, 19 functions imported by
  `third_party/core-math/import.cmake` with three mechanical edits (the "impossible hard case"
  path of atan2 became a fail-fast check); its `fma` is routed to `myb_fma`, the library's own
  exact implementation in 128-bit integer arithmetic.
- Double-word arithmetic and compensated sums: Dekker (1971), Knuth TAOCP vol. 2 (1997), Joldes,
  Muller, Popescu (ACM TOMS 44, 2017), Muller & Rideau (ACM TOMS 48, 2022), Ogita, Rump, Oishi
  (SIAM J. Sci. Comput. 26, 2005).
- Elliptic: Carlson's duplication (Numer. Algorithms 10, 13, 1995; R_J in the form of Carlson
  1988 with R_C), DLMF 19 and 22, A&S 17.7.15 for Π at n < 0.
- Linear algebra: Parlett–Reinsch balancing, Hessenberg reduction and the implicit Francis
  double-shift QR after Golub & Van Loan (Alg. 7.5.2), inverse iteration.
- erf and erfc (internal, used by the charging currents): derived from `s_erf.c` of FreeBSD's
  libm (Sun's fdlibm, the SunPro notice kept; see `THIRD_PARTY_NOTICES.md`).
- Parsing: Clinger's fast path (PLDI 1990), a double-word candidate and an exact big-integer
  comparison with the midpoints; the 768-digit bound of Lemire (2021). Printing: digit
  generation after Steele and White (PLDI 1990) on big integers.
- Text kernels: NAIF "Kernel Required Reading". CCSDS 502.0-B-2/B-3 (Orbit Data Messages),
  505.0-B (XML Specification for Navigation Data Messages).

## Interface

- [`myb/cmath.h`](../../include/myb/cmath.h): `myb_fma`, `myb_sin`, `myb_cos`, `myb_tan`,
  `myb_asin`, `myb_acos`, `myb_atan`, `myb_atan2`, `myb_exp`, `myb_expm1`, `myb_log`,
  `myb_log1p`, `myb_sinh`, `myb_cosh`, `myb_tanh`, `myb_asinh`, `myb_acosh`, `myb_atanh`,
  `myb_cbrt`, `myb_hypot`.
- [`myb/num.h`](../../include/myb/num.h): `myb_two_sum`, `myb_two_prod`; `myb_dd`,
  `myb_dd_from_double`, `myb_dd_to_double`, `myb_dd_add`, `myb_dd_sub`, `myb_dd_add_d`,
  `myb_dd_mul`, `myb_dd_mul_d`, `myb_dd_div`, `myb_dd_div_d`, `myb_dd_sqrt`, `myb_dd_neg`;
  `myb_sum_acc`, `myb_sum_acc_add`, `myb_sum_acc_value`, `myb_sum2`, `myb_dot2`;
  `myb_carlson_rf`, `myb_carlson_rd`, `myb_carlson_rj`, `myb_carlson_rc`, `myb_elliptic_k`,
  `myb_elliptic_e`, `myb_elliptic_pi`, `myb_elliptic_f_inc`, `myb_elliptic_e_inc`,
  `myb_elliptic_pi_inc`, `myb_jacobi`; `myb_lu_solve`, `myb_cholesky`, `myb_cholesky_solve`,
  `myb_least_squares`, `myb_eigenvalues`, `myb_eigenvector_real`.
- [`myb/io.h`](../../include/myb/io.h): `myb_parse_double`; `MYB_DOUBLE_CHARS`,
  `myb_format_double`, `myb_format_fixed`, `myb_format_scientific`; `myb_table`,
  `myb_table_load`, `myb_table_load_file`, `myb_table_destroy`, `myb_table_rows`,
  `myb_table_columns`, `myb_table_name`, `myb_table_row`, `myb_table_write`; `myb_kpool`,
  `myb_kvar_type`, `myb_kpool_create`, `myb_kpool_destroy`, `myb_kpool_load_text`,
  `myb_kpool_load_file`, `myb_kpool_describe`, `myb_kpool_numbers`, `myb_kpool_string`.
- [`myb/ccsds.h`](../../include/myb/ccsds.h): `myb_odm_notation` (`MYB_ODM_KVN`,
  `MYB_ODM_XML`), `MYB_ODM_TEXT`, `MYB_OPM_MAX_MANEUVERS`; OEM — `myb_oem`, `myb_oem_load`,
  `myb_oem_load_file`, `myb_oem_destroy`, `myb_oem_header`, `myb_oem_segment_count`,
  `myb_oem_meta`, `myb_oem_state_count`, `myb_oem_state`, `myb_oem_covariance_count`,
  `myb_oem_covariance`, `myb_oem_segment_data`, `myb_oem_data`, `myb_oem_write`,
  `myb_oem_write_as`; OPM — `myb_opm`, `myb_opm_maneuver`, `myb_opm_parse`, `myb_opm_write`,
  `myb_opm_write_as`; OMM — `myb_omm`, `myb_omm_parse`, `myb_omm_write`, `myb_omm_write_as`;
  OCM — `myb_ocm`, `myb_ocm_kind`, `myb_ocm_load`, `myb_ocm_load_file`, `myb_ocm_destroy`,
  `myb_ocm_header`, `myb_ocm_meta`, `myb_ocm_block_count`, `myb_ocm_value`,
  `myb_ocm_line_count`, `myb_ocm_line_text`, `myb_ocm_line`, `myb_ocm_state`, `myb_ocm_pair`,
  `myb_ocm_block_data`, `myb_ocm_data`, `myb_ocm_write`, `myb_ocm_format_line`.

TLE and the OMM ⇄ TLE conversion are in `ext/sgp4` ([SGP4 and TLE](sgp4-tle.en.md)); SPK and
SP3 files in [SPK ephemerides](spk-ephemerides.en.md).

## Inputs and outputs

- `cmath`: every function returns RN(f(x)); special values follow IEEE 754 / C Annex F;
  floating-point flags and errno are not part of the contract. The build enforces binary64
  with `FLT_EVAL_METHOD == 0`, `-ffp-contract=off` and no fast-math; `test_no_libm_imports`
  checks that the library imports no inexact libm function.
- Double-word numbers (hi, lo) with hi = RN(hi + lo), about 106 bits; proven relative bounds
  (u = 2⁻⁵³): add/sub 3u², add_d 2u², mul 7u², mul_d 1.5u² + 4u³, div 15u² + 56u³, div_d
  3.5u²; sqrt (one Karp–Markstein correction) has no published bound and is required < 4u² by
  the tests. FMA-free variants on purpose (the x86-64-v2 baseline has no hardware FMA).
  `myb_two_prod` is exact if no overflow and e_a + e_b ≥ −970.
- Sum2: |result − S| ≤ u|S| + γ²ₙ₋₁ Σ|xᵢ|; Dot2: ≤ u|x·y| + γ²ₙ |x|·|y|; the order is the
  caller's, so the result is deterministic.
- Elliptic: x, y, z ≥ 0 with at most one zero, p > 0; m = k² ≤ 1, n < 1; incomplete forms for
  any real φ by quasi-periodicity; Jacobi for 0 ≤ m ≤ 1.
- Linear algebra: row-major matrices of order n ≤ 64, a fixed order of operations (the same
  bits everywhere); the functions keep their own stack below 1 KB.
- `myb_parse_double`: decimal (with `eEdD` exponents), hexadecimal, `inf`, `infinity`, `nan`;
  no leading white space skipped; the decimal separator is always '.'; round to nearest, ties
  to even, overflow to infinity, gradual underflow; 800 significant digits kept, beyond them a
  sticky flag.
- Printing: the shortest decimal that reads back (Python's `repr` form, positional for
  10⁻⁴ ≤ |x| < 10¹⁶, at most 24 characters); `decimals` 0–60 as `%.Nf`; `digits` 1–60 as
  `%.(digits−1)e`; ties to even; the full length is always reported.
- Tables: commas, semicolons, tabs or blanks; `#` comments; a first row of non-numbers names
  the columns (31 characters kept); an empty field is a NaN; written with the shortest
  decimals, at most 1024 columns.
- Kernels: `NAME = value`, `NAME = ( … )`, `NAME += ( … )` inside `\begindata … \begintext`;
  Fortran `D` exponents; strings in single quotes; @-dates kept as strings; loading is atomic;
  paths are UTF-8; JPL's own kernels (`gm_de440.tpc`, `pck00011.tpc`) load unchanged.
- CCSDS: in the file km, km/s, km/s², degrees, km³/s²; in the API SI (OMM mean elements stay in
  the units of the message: they mean something only with their theory). Epochs are calendar
  labels on the message's `TIME_SYSTEM`, converted by the caller. Covariances are 6 × 6 by rows
  (the file holds the lower triangle). The reader tells the notation by the first character
  ('<' for XML). A number is written as the shortest decimal in the file's unit that the
  reader turns back into the same SI value; seconds of epochs get as many decimals as they
  need, at least three. OEM is written as version 2.0, OCM as version 3.0. The
  writers do not allocate: text goes into the caller's buffer.

## Domain of validity

- Double-word operations need finite arguments and results that neither overflow nor
  underflow. Elliptic arguments outside their domains give NaN.
- Linear algebra: `MYB_E_DOMAIN` for a singular, non-positive-definite, rank-deficient or
  non-finite matrix; `myb_eigenvalues` gives `MYB_E_NOCONV` after 30 n sweeps.
- Parsing: `MYB_E_FORMAT` when no number starts at the text (outputs untouched).
- Printing: `MYB_E_DOMAIN` for decimals or digits out of range, `MYB_E_OVERFLOW` for a short
  buffer (the text is cut but terminated).
- Tables and kernels: `MYB_E_FORMAT` (diag gives the line), `MYB_E_NOMEM`, `MYB_E_OVERFLOW`,
  `MYB_E_IO`; kernel queries `MYB_E_NOT_FOUND`, `MYB_E_FORMAT` (wrong type), `MYB_E_RANGE`
  (string index).
- CCSDS readers: `MYB_E_FORMAT` (diag gives the line: not such a message, a missing mandatory
  keyword, a bad number or date, more than 16 OPM manoeuvres); writers `MYB_E_DOMAIN` (no
  segments, missing text, values not finite, invalid dates) and `MYB_E_OVERFLOW` (call again
  with length + 1 bytes); `myb_ocm_state` gives `MYB_E_DOMAIN` for a `TRAJ_TYPE` other than
  CARTP, CARTPV, CARTPVA.

## Accuracy

| What | Reference | Result |
|---|---|---|
| The 19 CORE-MATH functions: random inputs over the whole domain, subnormals, arguments near kπ/2, samples of the CORE-MATH worst cases | GNU MPFR | every result correctly rounded, bit for bit |
| `myb_fma`: fixed cases, deep cancellation, operands binades apart | hardware FMA (FMA3 CPUs) | bit for bit (mingw-w64's `fma()` is wrong on about 1 input in 6) |
| `myb_two_sum`, `myb_two_prod`, 200 000 inputs each | MPFR | exact |
| Double-word operations | MPFR | within the proven bounds above (sqrt < 4u²) |
| Sum2, Dot2, 200 ill-conditioned trials each | MPFR | no bound violated |
| Carlson, Legendre, Jacobi, 1650 values | mpmath, 60 digits | ≤ 5.6 ulp; the amplitude ≤ 7.4 ulp (accounting for \|u\|) |
| erf, erfc, 23 values | mpmath | within 2 ulp (`test_charging`) |
| LU, Cholesky, least squares; eigenvalues of companion matrices, rotations, a spectrum 1500 … 1/1500, 50 random matrices, a Jordan block | exact solutions, known roots, trace and determinant | solutions to 10⁻¹⁴ |
| `myb_parse_double`: random short and long decimals, hex, exact midpoints and strings next to them, longer than 800 digits | `mpfr_strtofr` | every one correctly rounded |
| Printing in the three forms, 6124 values | CPython 3.13 | not one difference |
| Shortest form, 200 000 random doubles | reading back | the same bits |
| Table of 3500 random numbers | writing and reading back | bit for bit |
| OEM written and read, 2000 random states (12 000 numbers) | the states | exact but 183 numbers (1.5 %), 1 ulp off (no number of km times 1000 reaches them); epochs bit for bit |
| Covariances of OEM, OPM, OMM in both notations | written and read | to 1 ulp of the km² ↔ m² conversion |
| OCM, 8 blocks of the 7 kinds | read, written in both notations, read again | the same message |

Tests: [`test_mpfr_cmath`](../../tests/mpfr/test_mpfr_cmath.c),
[`test_fma`](../../tests/unit/test_fma.c), [`test_mpfr_num`](../../tests/mpfr/test_mpfr_num.c),
[`test_num`](../../tests/unit/test_num.c), [`test_elliptic`](../../tests/unit/test_elliptic.c),
[`test_linalg`](../../tests/unit/test_linalg.c),
[`test_charging`](../../tests/unit/test_charging.c) (erf),
[`test_mpfr_parse`](../../tests/mpfr/test_mpfr_parse.c),
[`test_parse`](../../tests/unit/test_parse.c), [`test_format`](../../tests/unit/test_format.c),
[`test_table`](../../tests/unit/test_table.c), [`test_kpool`](../../tests/unit/test_kpool.c)
(syntax, errors with lines, atomic loading, a Cyrillic path),
[`test_ccsds`](../../tests/unit/test_ccsds.c) (563 checks).

## Cost

Not measured.

## Limitations

- No `pow` among the correctly rounded functions (deferred in phase 2); erf and erfc are
  internal, not part of the public header.
- The double-word square root has no published error bound; only the tests' sample bounds it.
- Linear algebra is for the small systems (n ≤ 64) of the library's own solvers.
- XML of CCSDS: the elements and attributes of the messages are read, a value may be one CDATA
  section, and a document type with an internal subset is passed over — but the entities it
  declares are not known (`&name;` stays as written); a CDATA section mixed with other text,
  two sections in a row and an unterminated one are refused.
- Comments and the user-defined keywords of OEM, OPM, OMM are passed over. OCM is given as it
  is — keywords by name, lines as text and numbers in the message's units; only Cartesian
  trajectories are also given in SI.
- Converting the epochs of a message to a time scale is the caller's
  ([time scales](time-scales.en.md)).

## Data

None at run time. Test references in the repository: samples of the CORE-MATH worst cases
`tests/data/core-math/*.wc` (MIT), `tests/data/elliptic/reference.txt` (computed with mpmath
1.4.1 at 60 digits), `tests/data/io/format.ref` (CPython 3.13). GNU MPFR and GMP are linked
into the test executables only.
