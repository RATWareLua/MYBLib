# Tides of a gravity field (IERS 2010)

**English** | [Русский](tides-iers.ru.md) · [All passports](index.en.md)

## What it is

The tidal changes of the coefficients of a gravity field after the IERS Conventions 2010,
chapter 6: the solid tide in its three steps (the Love numbers, their frequency dependence, the
permanent tide), the solid and the ocean pole tides about the secular pole, the ocean tides
from the Stokes coefficients of a model such as FES2004, the secondary waves of the ocean tide
by a linear admittance, and waves added on top (the mean S1, the atmospheric S1 and S2). The
formulas are general and every number comes from a kernel — for the Earth
`data/grav/iers2010-tides.tk`. The changes are added to the field as corrections until its next
change of time.

## Source

- IERS Conventions (2010), IERS Technical Note 36, chapter 6 (version of 2018-02-01): eqs. 6.6,
  6.7 and Table 6.3 (step 1: complex Love numbers, k⁽⁺⁾ for degree 4); eq. 6.8 and
  Tables 6.5a–c (step 2: 71 rows); eqs. 6.13, 6.14 (step 3); sections 6.4, 6.5 and eq. 6.24
  (pole tides); eq. 6.15 (ocean tides); eq. 6.16 and Table 6.7 (secondary waves). Chapter 7,
  eq. 21 (the secular pole).
- The kernel holds only numbers, with the SHA-256 of the source PDFs. The 71 rows of
  Tables 6.5 were converted by a script that checks, in every row, the Doodson number against
  its six multipliers and those against the five Delaunay multipliers; the amplitudes are
  copied unchanged. Two pivots printed in Table 6.7 are not waves of the table and are
  corrected in the kernel with a note: 135.455 for 135.655 (Q1, the pivot of 117.655) and
  237.755 for 235.755 (2N2, the pivot of 245.555 and 245.645).
- Ocean tides: files of the form of the IERS `fes2004_Cnm-Snm.dat` (the FES2004 model of
  F. Lyard et al., converted to Stokes coefficients for the Conventions).
- Implemented from the Conventions; no IERS code.

## Interface

Header [`myb/grav.h`](../../include/myb/grav.h):

- `myb_tides_create`, `myb_tides_destroy`, `myb_tides_load_ocean`,
  `myb_tides_load_ocean_file`, `myb_tides_add_waves`, `myb_tides_add_waves_file`,
  `myb_tides_ocean_waves`, `myb_tides_secondary_waves`, `myb_tides_apply`;
- structures `myb_tide_raiser`, `myb_tides_input`; the enum `myb_permanent_tide`
  (`MYB_PERMANENT_TIDE_AUTO`, `_KEEP`, `_REMOVE`); flags `MYB_TIDES_SOLID`,
  `MYB_TIDES_FREQUENCY`, `MYB_TIDES_OCEAN`, `MYB_TIDES_POLE_SOLID`, `MYB_TIDES_POLE_OCEAN`,
  `MYB_TIDES_ADMITTANCE`, `MYB_TIDES_ALL`.

In a force model the members `tides`, `tide_flags`, `permanent` and `raisers` of
`myb_field_term` apply them once per instant ([spherical harmonics](spherical-harmonics.en.md)).
The tides in the displacement of a station: [site displacement](site-displacement.en.md).

## Inputs and outputs

- A raiser: its GM (m³/s²) and its position in the body-fixed frame of the field (m); GM and R
  of eq. 6.6 are those of the field.
- Step 1: ΔC̄ₙₘ − iΔS̄ₙₘ = kₙₘ/(2n + 1) Σⱼ (GMⱼ/GM)(R/rⱼ)ⁿ⁺¹ P̄ₙₘ(sin φⱼ) e^(−imλⱼ), n = 2, 3;
  degree 4 from the degree-2 tide with k₂ₘ⁽⁺⁾.
- Step 2: C̄₂₀, C̄₂₁, S̄₂₁, C̄₂₂, S̄₂₂ as Poisson series in χ = GMST + π and the Delaunay
  arguments. The input time `t` is UT1 (the kernel counts all its arguments in Julian centuries
  of UT1); it is also the date of the secular pole.
- Step 3: for a zero-tide field the permanent part A₀H₀ Re k₂₀ is removed from C̄₂₀. `AUTO`
  takes it from the field: `zero_tide` → remove, `tide_free` → keep; a `mean_tide` field needs
  an explicit choice.
- Pole tides: ΔC̄₂₁ = c₁(m₁ + a₁m₂), ΔS̄₂₁ = c₂(m₂ + a₂m₁), m₁ = xp − xs, m₂ = −(yp − ys) about
  the secular pole of chapter 7; xp, yp in rad.
- Ocean tides (eq. 6.15): ΔC̄ₙₘ − iΔS̄ₙₘ = Σ_f (C⁺ − iS⁺) e^(iθ_f) + (C⁻ + iS⁻) e^(−iθ_f), from
  rows (Doodson number, name, n, m, C⁺, S⁺, C⁻, S⁻ in units of 10⁻¹¹); degree 1 is skipped,
  S̄ₙ₀ is not changed, coefficients above the field's degree are skipped.
- Secondary waves (`MYB_TIDES_OCEAN | MYB_TIDES_ADMITTANCE`, eq. 6.16): a wave the file does
  not give takes its coefficients from two pivot waves it gives, the admittance linear in the
  frequency, H the astronomical amplitudes. No new coefficients: the phasors of the secondary
  waves are added to the phasor of their pivot.
- `myb_tides_add_waves`: waves of the same form added to those loaded; a wave already loaded
  counts twice (its coefficients add).
- Kernel variables: `<prefix>_ARGUMENT_*`, the series `_C20` … `_S22`, `_LOVE` (n, m, Re k,
  Im k, k⁺), `_PERMANENT` (A₀, H₀), `_POLE_SOLID`, `_POLE_OCEAN` (c₁, a₁, c₂, a₂),
  `_SECULAR_POLE` (mas), `_ADMITTANCE` (Doodson numbers and H of the wave and its two pivots).

## Domain of validity

- Ocean tides of degree 2 … 360 (`MYB_E_DOMAIN` otherwise); secondary waves only where both
  pivots are loaded and the wave itself is not.
- `MYB_E_DOMAIN`: a raiser at the origin, input that is not finite, `AUTO` with a field that is
  neither zero-tide nor tide-free, a prefix longer than 20 characters.
- `MYB_E_FORMAT` (the diagnostic names the variable or the line), `MYB_E_NOT_FOUND`,
  `MYB_E_IO`, `MYB_E_NOMEM`; a refused text leaves the waves loaded before.

## Accuracy

Tests: [`test_tides`](../../tests/unit/test_tides.c), [`test_env`](../../tests/unit/test_env.c),
[`test_lageos`](../../tests/unit/test_lageos.c).

| What | Reference | Result |
|---|---|---|
| Step 1, degrees 2 and 3, the Moon and the Sun | the addition theorem k GMⱼ R²ⁿ⁺¹/(rⁿ⁺¹ rⱼⁿ⁺¹) Pₙ(cos ψ) (`test_tides`) | 2·10⁻¹⁶ |
| Step 2: K1, M2 and the 18.6-year term | the K1 example of the IERS text with ERFA's GMST (`test_tides`) | 1.3·10⁻¹⁶ |
| Step 3 | A₀H₀ Re k₂₀ (`test_tides`) | 10⁻²² absolute (bound) |
| Pole tides at 2025.5 | the formulas with the secular pole (`test_tides`) | 10⁻²² absolute (bound) |
| Ocean tides, an excerpt of FES2004 | eq. 6.15 with ERFA's fundamental arguments (`test_tides`) | 3·10⁻¹⁶ (the rest is GMST of 1982 against 2006, 0.056″) |
| Secondary wave L2 carried by M2 and K2 | the same wave given in the file with the coefficients of eq. 6.16 (`test_tides`) | 2·10⁻¹⁰ |
| FES2004 with Table 6.7 | — (`test_tides`) | 63 secondary waves on the 18 given |
| Waves added with `myb_tides_add_waves` | one file with both (`test_tides`) | bit for bit |
| The field term with tides | `myb_tides_apply` called by hand (`test_env`) | bit for bit |
| LAGEOS-1 over a week, the full model (FES2004 to degree 20, all flags) | the ILRS orbit (`test_lageos`) | RMS 1.57 cm after the fit; 2.05 cm without the secondary waves |

What the tides weigh on LAGEOS, measured with the model of 2.0 cm before the secondary waves
(RMS of the fit without the part): without the solid tides 21 cm, without the ocean tides
6.3 cm.

## Cost

Not measured.

## Limitations

- χ in the kernel is GMST of 1982, as in the IERS routines (0.056″ from ERFA's GMST of 2006 in
  2025); all arguments of step 2 and of the ocean tides are taken on UT1.
- Step 1 covers degrees 2 and 3 (and 4 through k⁽⁺⁾); step 2 the degree-2 coefficients of
  Tables 6.5; the pole tides C̄₂₁, S̄₂₁ only.
- One ocean model at a time, of the FES2004 form. The ILRS analysis centres use EOT11a and the
  atmospheric tides; on LAGEOS the once-per-revolution empirical terms absorb the difference.
- The atmospheric S1 and S2 enter only as geopotential coefficients given by the caller through
  `myb_tides_add_waves`; no model of them is built in.
- The kernel holds the Earth's numbers; another body needs its own kernel.
- Tides modify the field: one field per thread.

## Data

- In the repository: `data/grav/iers2010-tides.tk` (prefix `IERS10_TIDE`) — the Love numbers
  (Table 6.3), Tables 6.5a–c, the permanent and pole tide constants of chapter 6, the secular
  pole of chapter 7 and the 63 rows of Table 6.7, from the IERS Conventions (2010, version of
  2018), numbers only, with the SHA-256 of the PDFs; `tests/data/grav/fes2004_excerpt.dat`, an
  unchanged excerpt of FES2004.
- Downloaded by `cmake/FetchData.cmake`: `fes2004_Cnm-Snm.dat` (FES2004 as Stokes coefficients,
  IERS Conventions Centre; 18 waves).
- `test_lageos` also needs `GOCO06s.gfc`, `de440.bsp`, `header.440`, `eopc04.1962-now` and the
  ILRS excerpt `tests/data/sp3/ilrsa.lageos1.250607.excerpt.sp3`.
