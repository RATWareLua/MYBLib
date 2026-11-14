# Gravity fields in spherical harmonics

**English** | [Русский](spherical-harmonics.ru.md) · [All passports](index.en.md)

## What it is

The potential, the acceleration and the gradient of the acceleration of an extended body,
expanded in fully normalised spherical harmonics in its body-fixed frame: up to degree 2700
for the potential and the acceleration, 2698 for the gradient, with the poles as ordinary
points. The coefficients come from ICGEM files, with their time-variable terms, or from
arrays; corrections (tides, the pole tide) are added on top until the time changes again. A
ready-made force term evaluates the field in the inertial axes through the orientation of the
body, keeps the coefficients at the instant, gives its analytic Jacobian and, for the central
body, the indirect term of its figure.

## Source

- The evaluation: the forward column recursion of S. A. Holmes, W. E. Featherstone,
  *J. Geodesy* 76, 279 (2002), with the sectorial seeds P̄ₘₘ scaled by 2⁹³⁰, so that degree
  2700 stays within range at every latitude. The derivatives have no division by cos φ: the
  columns P̄ₙₘ / cos φ are recursed from their own seeds, and
  dP̄ₙₘ/dφ = eₙₘ P̄ₙ,ₘ₊₁ − m sin φ P̄ₙₘ / cos φ.
- The gradient: the second derivatives of a solid harmonic of degree n are solid harmonics of
  degree n + 2 (Cunningham 1970; Montenbruck & Gill, eq. 3.33); the five independent
  components are series of the same form, summed by the same column recursion.
- The files: the ICGEM formats 1.0 and 2.0 (`gfc`; `gfct`, `trnd`, `acos`, `asin`, epochs
  `yyyymmdd[.hhmm]`).
- Implemented from the papers; no third-party code.

## Interface

Header [`myb/grav.h`](../../include/myb/grav.h):

- the field: `myb_gravity_field_create`, `myb_gravity_field_load_icgem`,
  `myb_gravity_field_load_icgem_file`, `myb_gravity_field_clone`,
  `myb_gravity_field_destroy`, `myb_gravity_field_degree`, `myb_gravity_field_gm`,
  `myb_gravity_field_radius`, `myb_gravity_field_tide_system`, `myb_gravity_field_set_time`,
  `myb_gravity_field_add_correction`, `myb_gravity_field_coefficient`;
- the evaluation: `myb_gravity_work_create`, `myb_gravity_work_destroy`,
  `myb_gravity_field_eval`, `myb_gravity_field_gradient`.

Header [`myb/env.h`](../../include/myb/env.h): the structure `myb_field_term` with
`myb_term_field`, `myb_term_field_jacobian`, `myb_term_field_indirect`.

The tidal changes of the coefficients: [IERS 2010 tides](tides-iers.en.md); a field as a level
of a system of bodies (`myb_system_set_field`, `myb_system_add_field`):
[system of bodies](system-of-bodies.en.md).

## Inputs and outputs

- V(r, φ, λ) = GM/r Σₙ (R/r)ⁿ Σₘ P̄ₙₘ(sin φ) (C̄ₙₘ cos mλ + S̄ₙₘ sin mλ): geocentric latitude,
  fully normalised associated Legendre functions in the geodesy convention, without the
  Condon–Shortley phase. Coefficient arrays are indexed n (n + 1)/2 + m.
- Positions in m in the body-fixed frame. The potential in m²/s², **positive** (GM/r for a
  point mass); the acceleration in m/s², its gradient; the gradient in s⁻²,
  `gradient[3 i + j]` = ∂aᵢ/∂rⱼ, symmetric and, outside the body, trace-free.
- A range of degrees n_min … n_max (n_max = −1: all; n_min = 2 leaves out the central and the
  degree-1 terms).
- ICGEM keys: `earth_gravity_constant` (or `gravity_constant`), `radius`, `max_degree`, `norm`
  (`fully_normalized` only), `tide_system`, `format`; truncation by degree on loading; absent
  coefficients are zero. The declared tide system (`tide_free`, `zero_tide`, `mean_tide` or
  empty) is reported, not converted.
- Effective coefficients = static + time-variable terms at t (TT, as the models use it; trends
  per year, periods in years) + the corrections added since the last `set_time`.
- The term: x relative to the origin of the environment in its inertial axes; the field turned
  by the orientation of the body. For a body that is not at the origin (with a body at the
  origin) the same field at the origin is subtracted — the indirect part. With `update` the
  coefficients follow time; with `tides` the tidal changes raised by `raisers` are applied once
  per instant, at the instant on the scale of the body's orientation (TT for the Earth).
- `myb_term_field_indirect`: a = Σⱼ (GMⱼ / GM) g(rⱼ), g the field without its degrees 0 and 1
  at the raiser j — the acceleration of the origin by the figure of the central body. It does
  not depend on the object.

## Domain of validity

- Degrees 0 … 2700 (`myb_gravity_field_create` refuses 2701); the gradient up to n_max = 2698;
  the work must be created for at least n_max.
- The series converge outside the Brillouin sphere of the body. Inside it they may diverge:
  for a 4 × 2.4 × 1.6 box the exact spherical series goes from 3.8·10⁻² to 0.26 of error from
  degree 4 to 16 (`test_ellipsoidal`) — for such bodies see
  [irregular bodies](irregular-bodies-gravity.en.md).
- `MYB_E_DOMAIN`: r = 0 or not finite, degrees outside the field or the work, gm ≤ 0,
  radius ≤ 0, coefficients or corrections not finite, a correction outside the field; the
  indirect term of a body that is not at the origin.
- `MYB_E_FORMAT` (the diagnostic gives the line): no GM, a normalisation other than
  `fully_normalized`, no `end_of_head`, a short line, m > n, `gfct` without its epoch, a
  periodic term without its period, a trend without `gfct`. Also `MYB_E_NOMEM`,
  `MYB_E_OVERFLOW`.
- `set_time`, corrections and tides modify the field: not concurrently with evaluations. One
  field (`myb_gravity_field_clone`) and one work per thread.

## Accuracy

Tests: [`test_mpfr_grav`](../../tests/mpfr/test_mpfr_grav.c),
[`test_grav`](../../tests/unit/test_grav.c), [`test_gradient`](../../tests/unit/test_gradient.c),
[`test_env`](../../tests/unit/test_env.c), [`test_orbit`](../../tests/unit/test_orbit.c),
[`test_nongrav`](../../tests/unit/test_nongrav.c), [`test_lageos`](../../tests/unit/test_lageos.c).

| What | Reference | Result |
|---|---|---|
| Potential / acceleration, degree 120, random coefficients, at a pole, next to it, at the equator, from the surface to 3R | MPFR 256 bits, the acceleration as differences of that potential (`test_mpfr_grav`) | 4.7·10⁻¹⁶ / 7.5·10⁻¹⁶ |
| The same at degree 2190, latitude 62° (unscaled seeds would underflow) | MPFR 256 bits (`test_mpfr_grav`) | 4.8·10⁻¹⁶ / 2.5·10⁻¹⁶ |
| Point mass and J₂, 200 points | closed forms (`test_grav`) | below 7·10⁻¹⁶ and 2·10⁻¹⁵ (the bounds of the test) |
| Degrees 0..1 plus 2..2 | all degrees at once (`test_grav`) | 4·10⁻¹⁶ (bound) |
| Time-variable terms (`gfct`, `trnd`, `acos`, `asin`) | the sums by hand (`test_grav`) | 10⁻²⁶ absolute (bound) |
| Gradient, degree 8, 32 cases at all latitudes and at the poles | second derivatives of the series in mpmath, 60 digits (`test_gradient`) | 2.1·10⁻¹⁵ of the largest component, 6·10⁻¹⁶ at the poles |
| Gradient, degree 60 | central differences of the acceleration (`test_gradient`) | 2·10⁻¹⁰ |
| The term: with and without tides, a part of the degrees, a body away from the origin | the same functions called by hand (`test_env`) | bit for bit |
| Jacobian of the term with tides | central differences (`test_env`) | 1·10⁻¹⁰ |
| State transition matrix, rotating field of degree 8, two revolutions | symplecticity; differences of whole propagations (`test_orbit`) | 6·10⁻¹⁷; 5·10⁻¹⁰ (10⁻¹² with a Jacobian by differences) |
| Indirect term: the Earth's J₂ answering the Moon | J₂ in closed form (`test_nongrav`) | −1.4834·10⁻¹¹ m/s², to 10⁻¹² |
| LAGEOS-1 over a week, GOCO06s to degree 30 with trends and annual terms inside the full model | the ILRS orbit (`test_lageos`) | RMS 1.57 cm after the fit |

## Cost

Not measured.

## Limitations

- Fully normalised coefficients only. The tide system of a file is reported, not converted;
  the permanent tide of a zero-tide field is handled by the tides
  ([IERS 2010 tides](tides-iers.en.md)).
- Time-variable terms of the ICGEM 2.0 kinds only: a static part with its epoch, a trend, and
  cosine and sine terms with their periods.
- No special treatment inside the Brillouin sphere: the series is used as it is.
- A field object holds state (the time, the corrections): one object per thread.
- The indirect term has no Jacobian (`myb_term_no_jacobian`): it does not depend on the object.

## Data

- None is required: coefficients can be given as arrays.
- In the repository, for the tests: `tests/data/grav/EGM2008_n20.gfc` (EGM2008 of
  N. K. Pavlis et al., NGA, truncated to degree 20, from ICGEM; coefficients unchanged),
  `synthetic2.gfc` (a format 2.0 file with time-variable terms), `gradient.ref` (the mpmath
  reference of the gradient).
- Downloaded by `cmake/FetchData.cmake`: `GOCO06s.gfc` (degree 300 with trends and annual
  terms; A. Kvas et al. 2019, ICGEM, CC BY 4.0), used by `test_lageos`.
