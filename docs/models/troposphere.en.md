# Tropospheric delay

**English** | [Русский](troposphere.ru.md) · [All passports](index.en.md)

## What it is

The delay of radio and optical signals in the Earth's neutral atmosphere, in an optional
module: zenith hydrostatic and wet delays, mapping functions to a slant ray, horizontal
gradients and the empirical models of the surface meteorology. The slant delay is
m_h(e) ZHD + m_w(e) ZWD + the gradient term; it lengthens the light time of
[observables](observables-light-time.en.md) by delay / c. The models are those of the IERS
Conventions 2010, chapter 9 (Saastamoinen, Mendes–Pavlis, FCULa, VMF1, GMF, GPT, Chen–Herring)
and of Landskron and Böhm 2018 (VMF3, GPT3 with its gradients, which replace APG), with the wet
delay of Askne and Nordius.

## Source

- J. Saastamoinen (1972) in the form of J. L. Davis et al., *Radio Sci.* 20, 1593 (1985);
  V. B. Mendes, E. C. Pavlis, *Geophys. Res. Lett.* 31, L14602 (2004); V. B. Mendes et al.,
  *GRL* 29, 1414 (2002) (FCULa); J. Böhm, B. Werl, H. Schuh, *J. Geophys. Res.* 111, B02406
  (2006) (VMF1); J. Böhm et al., *GRL* 33, L07304 (2006) (GMF); J. Böhm, R. Heinkelmann,
  H. Schuh, *J. Geod.* 81, 679 (2007) (GPT); A. E. Niell, *JGR* 101, 3227 (1996) (height
  correction); G. Chen, T. A. Herring, *JGR* 102, 20489 (1997) (gradients); D. Landskron,
  J. Böhm, *J. Geod.* 92, 349 (2018) (VMF3, GPT3); J. Askne, H. Nordius, *Radio Sci.* 22,
  379 (1987).
- `troposphere.c` is written from the published equations; it is not a translation of the
  IERS Conventions software, whose licence is not BSD-compatible
  ([`ext/troposphere/NOTICE`](../../ext/troposphere/NOTICE)). The coefficient tables of GMF
  and GPT are numbers from the DATA statements of `GMF.F` and `GPT.F` (SHA-256 in
  `troposphere_tables.inc`); the coefficients of b and c of VMF3 are numbers from
  `vmf3_ht.f90` of TU Wien (SHA-256 in `vmf3_tables.inc`). The routines of TU Wien carry no
  licence terms: VMF3/GPT3 are written from the paper, the tables and the grid are data, the
  routines serve only for comparison ([PLAN §9](../PLAN.en.md#9-licenses)).

## Interface

Header [`myb_troposphere.h`](../../ext/troposphere/myb_troposphere.h) (target
`myb::troposphere`, not part of libmyb):

- zenith delays `myb_tropo_zhd_saastamoinen`, `myb_tropo_zwd_saastamoinen`,
  `myb_tropo_zenith_optical`, `myb_tropo_zwd_askne`;
- mapping functions `myb_tropo_fcula`, `myb_tropo_vmf1`, `myb_tropo_gmf`, `myb_tropo_vmf3`;
- gradients `myb_tropo_gradient_delay`, `myb_tropo_gradient_delay_parts`;
- meteorology `myb_tropo_gpt`; GPT3 `myb_gpt3_load`, `myb_gpt3_destroy`, `myb_gpt3_step`,
  `myb_gpt3_at` with the structure `myb_gpt3_values`.

## Inputs and outputs

- Angles in rad, heights in m, pressures in hPa, temperatures in K (GPT and GPT3 give theirs
  in °C), time as an MJD (the seasonal terms need no better than a day); delays in m.
- ZHD = 0.0022768 P / (1 − 0.00266 cos 2φ − 0.28·10⁻⁶ h) (φ geodetic latitude, h
  ellipsoidal height); Saastamoinen's ZWD = 0.002277 (1255 / T + 0.05) e.
- Mendes–Pavlis: latitude, height, pressure, water-vapour pressure, wavelength (µm,
  0.3–1.2), CO₂ content (ppm; the Conventions use 375) → hydrostatic and non-hydrostatic
  zenith delays. FCULa: elevation, latitude, height above mean sea level, temperature.
- VMF1 and VMF3: the coefficients a_h, a_w of the site at the epoch (of the Vienna service,
  or of GPT3 for VMF3), MJD, latitude (and longitude for VMF3), ellipsoidal height (0: no
  height correction) and zenith distance → m_h, m_w. GMF: MJD, latitude, longitude, height,
  zenith distance. GPT: pressure, temperature (°C), geoid undulation.
- Gradients: north and east gradients in m, azimuth from north towards east;
  m_g = 1 / (sin e tan e + 0.0032); the parts version gives hydrostatic and wet gradients
  their own constants, 0.0031 and 0.0007 (the gradients of GPT3).
- Askne–Nordius: water-vapour pressure, mean temperature of the vapour, its decrease factor
  λ > −1: 10⁻⁶ (k₂′ + k₃ / T_m) R_d / ((λ + 1) g_m) e.
- GPT3 (`gpt3_5.grd` or `gpt3_1.grd`): pressure, temperature (°C) and its lapse rate (K/km),
  the mean temperature and the pressure of the water vapour, λ, a_h and a_w, the geoid
  undulation and the four gradients — the four nodes around the point, each brought to the
  height of the point, linear between them (within half a step of a pole: the node itself);
  `constant ≠ 0` gives the means without the seasons. A longitude a turn away is the same
  point.

## Domain of validity

- Refused with `MYB_E_DOMAIN`: input that is not finite, a zenith distance beyond 90° (GMF,
  VMF1, VMF3), an elevation not in (0, 90°] (FCULa, the gradients), a negative pressure, a
  temperature that is not positive, λ ≤ −1, a GPT height outside −1…40 km, a GPT3 latitude
  outside ±90°. Mendes–Pavlis is meant for 0.3–1.2 µm (the header) and refuses a wavelength
  outside 0.25–2 µm.
- `myb_gpt3_load`: `MYB_E_FORMAT` (a torn line, a count that is no grid, nodes out of order;
  the line in `diag`), `MYB_E_NOMEM`.
- Saastamoinen's ZWD is good to a few centimetres (the header); estimate the wet delay from
  the data when better is needed.

## Accuracy

| What | Reference | Result |
|---|---|---|
| FCULa, VMF1 with and without the height, GMF, GPT | the test cases in the headers of `FCUL_A.F`, `VMF1.F`, `VMF1_HT.F`, `GMF.F`, `GPT.F` (19 digits) | 10⁻¹⁶…10⁻¹⁵ (bounds 2·10⁻¹⁵; GPT temperature and undulation 2·10⁻¹⁴) |
| Mendes–Pavlis ZHD and ZWD | the case of `FCUL_ZD_HPA.F` | the published case is 1.96·10⁻⁶ smaller in both (4 µm): see Limitations |
| GPT3, 20 sites and dates | `gpt3_5.f90` of TU Wien (gfortran 15.3 on Gentoo) | **6·10⁻¹³** of the values |
| VMF3 | `vmf3_ht.f90` | 1.3·10⁻¹⁵ |
| Askne–Nordius wet delay | `asknewet.f` | 9·10⁻¹⁶ |
| The GPT3 gradients as a delay | the formula with 0.0031 and 0.0007 | 10⁻¹⁷ m |
| Mapping at the zenith; the closed forms of Saastamoinen | 1; the formulas | 10⁻¹⁵ |

The 20 GPT3 cases lie between nodes, at a node, across the first meridian and the date line,
within half a step of the poles, at sea level and on Everest, in a leap year after February,
and as means without the seasons; the sums of the seasons cancel in the small quantities, so
the last digits of cos and sin of the two libraries show there. The GPT3, VMF3 and
Askne–Nordius comparison runs only with `gpt3_5.grd` downloaded (skipped otherwise).

Test: `test_troposphere` ([tests/ext/test_troposphere.c](../../tests/ext/test_troposphere.c)).

## Cost

Not measured.

## Limitations

- The published case of `FCUL_ZD_HPA.F` does not match its own inputs: it is 1.96361·10⁻⁶
  (ZHD) and 2.00404·10⁻⁶ (ZWD) smaller; executing the routine's own statements with its
  single-precision literals `0.00266` and `3.759` makes both ratios 1 + 1.96361·10⁻⁶ — a
  common factor f(φ, H), as if the case had been run 7 m lower.
- No reader of the products of the Vienna service: a_h, a_w of VMF1 and VMF3 are the
  caller's numbers. Only the 5° grid of GPT3 is tested; the 1° grid is read by the same
  reader.
- No ray tracing; the gradients are a mapping function of horizontal gradients estimated
  elsewhere (or given by GPT3).
- The delay is not added to the measurements by itself: the measurements of
  [estimation](estimation.en.md) take media delays as corrections of the observed value.

## Data

- Built in: the GMF and GPT tables (`ext/troposphere/troposphere_tables.inc`) and the VMF3
  coefficients (`ext/troposphere/vmf3_tables.inc`).
- `gpt3_5.grd` — the GPT3 grid of TU Wien, 5°, 1.1 MB, in `cmake/FetchData.cmake` with its
  SHA-256; the tests skip GPT3 and VMF3 without it.
- `tests/data/tropo/gpt3_vmf3.ref` — outputs of the routines of TU Wien at the 20 cases, for
  comparison only.
