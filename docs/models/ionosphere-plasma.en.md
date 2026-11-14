# Ionosphere and plasma delay

**English** | [Русский](ionosphere-plasma.ru.md) · [All passports](index.en.md)

## What it is

The dispersive delay of radio signals by free electrons. In an optional module for the Earth:
global IONEX maps of the vertical total electron content (two- and three-dimensional) and
their interpolation; the single-layer model (pierce point, mapping function) that turns a map
into the slant content of a ray; the higher orders of the delay after IERS 2010 §9.4 — the
integrals along a ray, the coefficients s₁, s₂, s₃, the delays of the carrier phase and of the
code — over a Chapman layer on the two-dimensional maps or a grid of densities, and the
bending of Jakowski. In the core, for any body: the electron column of a power-law plasma
profile along a segment (the solar corona and wind) and the first-order group delay K N / f².
The delays are added to the light time of [observables](observables-light-time.en.md); the
neutral atmosphere is in [troposphere](troposphere.en.md).

## Source

- S. Schaer, W. Gurtner, J. Feltens, "IONEX: The IONosphere Map EXchange Format Version 1"
  (1998) and IONEX 1.1 (2015): the format, bilinear interpolation in space (eq. 4), linear in
  time between maps rotated with the Sun (eq. 3); 3-D maps — the values are the densities
  times DHGT (IONEX 1.0, table 2).
- S. Schaer, dissertation, University of Bern (1999): the single-layer model and the modified
  model of CODE.
- IERS Conventions 2010, §9.4, eq. 9.19–9.23, 9.37, 9.40, 9.46, 9.47; M. Fritsche et al.,
  *GRL* 32, L23311 (2005); N. Jakowski et al. (1994) for the bending. K₁ = 40.308,
  K₂ = 1.1283·10¹², K₃ = 812.38, K₄ = 1.5792·10²² follow from e, mₑ, ε₀ of CODATA 2022
  (`data/constants/codata-2022.tk`).
- Plasma profiles of the kind of D. O. Muhleman, J. D. Anderson, *ApJ* 247, 1093 (1981).
- Written by the authors of myb; no code of others
  ([`ext/ionosphere/NOTICE`](../../ext/ionosphere/NOTICE)). `myb_plasma_column` uses the
  15-point Gauss–Kronrod nodes and weights of QUADPACK (public domain). The magnetic field is
  the IGRF of [`myb/force.h`](../../include/myb/force.h)
  ([magnetic field](magnetic-field-charging.en.md)).

## Interface

Header [`myb_ionosphere.h`](../../ext/ionosphere/myb_ionosphere.h) (target
`myb::ionosphere`, not part of libmyb):

- maps `myb_ionex_load`, `myb_ionex_destroy`, `myb_ionex_shell`, `myb_ionex_dimension`,
  `myb_ionex_vtec`, `myb_ionex_density_grid`;
- single layer `myb_iono_pierce_point`;
- densities in space (`myb_electron_density` callbacks) `myb_iono_layer_density`,
  `myb_iono_grid_density`;
- orders `myb_iono_ray_integrals`, `myb_iono_orders_of`, `myb_iono_phase_delay`,
  `myb_iono_code_delay`, `myb_iono_bending`;
- structures `myb_iono_layer`, `myb_iono_grid`, `myb_iono_ray`, `myb_iono_orders`.

Header [`myb/observe.h`](../../include/myb/observe.h) (in libmyb): `myb_plasma_column`,
`myb_plasma_group_delay`, the structure `myb_plasma_profile` (up to `MYB_PLASMA_TERMS` = 8
terms).

## Inputs and outputs

- Angles in rad, heights and radii in m, time an MJD in UT; content in electrons/m²
  (1 TECU = 10¹⁶ m⁻²; about 16 cm of delay per TECU at L1).
- Single layer: sin z′ = R / (R + H) sin(αz), M = 1 / cos z′, z = 90° − elevation, azimuth
  from north towards east. α = 1, R = 6371 km, H = 450 km is the standard model of the IGS
  maps (a file's own values: `myb_ionex_shell`); CODE's modified model takes H = 506.7 km,
  α = 0.9782.
- Maps: RMS and height maps are skipped, 9999 is a missing value. A 3-D map keeps the
  densities (value / DHGT); its vertical content is the sum of density × DHGT. The grid of
  densities is in geocentric latitude, longitude and height above a sphere, at given epochs
  (`[map][height][latitude][longitude]`, m⁻³): linear between eight nodes and two epochs,
  each rotated with the Sun; nothing below the lowest and above the highest node; a grid that
  spans a turn goes round.
- Chapman layer over the 2-D maps: N(h) = N_max exp((1 − z − e^−z) / 2), z = (h − peak) /
  scale, N_max = VTEC / (scale √(2πe)); peak 300–450 km, scale 60–100 km (the header).
- Ray integrals along the straight line receiver–transmitter (Earth-fixed, m) between the radii
  `bottom` and `top`, four-point Gauss panels (1…4096): ∫N dl (m⁻²), ∫N B cos θ dl (T m⁻²),
  ∫N² dl (m⁻⁵), ∫N B² (1 + cos² θ) dl (T² m⁻²); θ is the angle of the field to the direction
  transmitter → receiver.
- s₁ (m Hz²), s₂ (m Hz³), s₃ (m Hz⁴); phase delay −s₁/f² − s₂/f³ − s₃/f⁴, code delay
  s₁/f² + 2s₂/f³ + 3s₃/f⁴ (m); the sign of s₂ is that of a right-hand circularly polarised
  signal (GNSS). The bending adds 2.495·10⁸ ((1 − 0.8592 cos² E)^−½ − 1) S² mm MHz⁴ (S in
  TECU) to s₃, returned in m Hz⁴.
- Plasma: n_e(r) = Σ n_k (R / r)^p_k (n_k ≥ 0, p_k > 0), positions relative to the centre of
  the body (m); the column in m⁻² along the straight segment, adaptive Gauss–Kronrod in the
  angle seen from the centre to a relative 10⁻¹³ (radial segments in closed form). Group delay
  +K N / f², phase delay −K N / f² (m), K = e² / (8π² ε₀ mₑ) = 40.308 m³ s⁻².

## Domain of validity

- `myb_ionex_load`: `MYB_E_FORMAT` (the line in `diag`; a dimension other than 2 and 3, a
  height off the grid), `MYB_E_NOMEM`. `myb_ionex_vtec`: `MYB_E_RANGE` (a time outside the
  maps, a latitude outside the grid, a missing value).
- `myb_iono_grid_density`: `MYB_E_RANGE` (date or latitude outside, a missing node around the
  place), `MYB_E_DOMAIN` (a grid that is none); `myb_ionex_density_grid`: `MYB_E_DOMAIN` for
  2-D maps. The layer: `MYB_E_DOMAIN` for a scale that is not positive.
- `myb_iono_pierce_point` refuses a negative elevation; `myb_iono_ray_integrals` refuses ends
  at one place, radii or panels that are none, input that is not finite (`MYB_E_DOMAIN`).
  Phase and code delays are NaN for a frequency that is not positive, the bending for invalid
  input.
- `myb_plasma_column`: `MYB_E_DOMAIN` (an invalid profile, a segment through the centre,
  input that is not finite), `MYB_E_NOCONV` (the quadrature budget). The header of
  `myb/observe.h` puts the orders beyond the first below 1 mm above 1 GHz on Earth paths;
  the example below gives +6.1 mm for the second order at 1575.42 MHz — use the module's
  higher orders for GNSS.

## Accuracy

| What | Reference | Result |
|---|---|---|
| IONEX 2-D interpolation, 60 points | an independent Python implementation of the IONEX description | 2·10⁻¹¹ — the resolution of an MJD in a double, 0.6 µs |
| 3-D maps, at the nodes and between them | the made field (changes of the exponent) | **2.2·10⁻¹⁵** |
| 3-D maps, between epochs and across the date line | the same | 1.6·10⁻¹² (resolution of the MJD) |
| Pierce point and mapping | the ray–shell intersection in Cartesian coordinates | 9·10⁻¹⁶ |
| Vertical ray through the Chapman layer | the content of the map; √(e/2π) = 0.6577 for ∫N² dl / (N_max S) | 10⁻⁷ |
| The four integrals, dipole field (IGRF 2020, degree 1) | direct sums on 200 000 steps | **1.3·10⁻⁹** |
| 16, 32, 64 panels on 3950 km | 256 panels | below 3·10⁻³, 10⁻⁵, 10⁻⁷ |
| K₁…K₄ from CODATA 2022 | the constants printed in the Conventions | to their printed digits |
| A grid of 25 km, 10° × 15° sampled from the layer | the layer | within 3 % |
| Plasma column: corona + wind, 1 → 1.5 au past 5 R☉; wind r^−2.3 past 50 R☉ | mpmath, 40 digits | 6·10⁻¹⁶ |
| Column for r⁻², r⁻⁴, radial r⁻², r⁻¹ | closed forms | 10⁻¹⁵ |

The header adds: panels as wide as the scale height give 4·10⁻⁸ of the content, twice as wide
3·10⁻⁶. At 30° of elevation the slant content differs from the single layer by 2.3 % and
∫N B cos θ dl by 4.7 % — the difference the integrals are computed for. An example on the
synthetic map: 33 TECU, 30°, 1575.42 MHz — the first order 9.12 m, the second +6.1 mm, the
third +0.17 mm; the bending at 100 TECU and 10° — 0.59 mm.

Tests: `test_ionosphere` ([tests/ext/test_ionosphere.c](../../tests/ext/test_ionosphere.c)),
`test_plasma` ([tests/unit/test_plasma.c](../../tests/unit/test_plasma.c)).

## Cost

Not measured.

## Limitations

- All IONEX checks are on synthetic files written by the authors (`synthetic.inx`,
  `synthetic3d.inx`), not on products of an analysis centre; there are no public 3-D
  products to check against.
- The profile over 2-D maps is a Chapman layer with the caller's peak and scale; the
  Conventions name no conventional profile, nor a conventional bending model (Jakowski's is a
  fit to ray tracing at the frequencies of GPS).
- The ray is a straight line; no ray tracing.
- The plasma column is first order only and its profile spherically symmetric; the delays are
  not added to the measurements by themselves — [estimation](estimation.en.md) takes them as
  corrections of the observed value.

## Data

- No downloads. IONEX files are the user's (IGS analysis centres: CODE, JPL, ESA, UPC, …).
- In the repository: `tests/data/ionosphere/synthetic.inx`, `synthetic3d.inx` and
  `expected.txt` (made by the authors of myb); `data/constants/codata-2022.tk` for e, mₑ, ε₀.
- The field for the higher orders: `igrf14coeffs.txt` (FetchData) or any
  `myb_magnetic_field`.
