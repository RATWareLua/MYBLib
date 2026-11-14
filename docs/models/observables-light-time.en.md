# Observables and light time

**English** | [Русский](observables-light-time.ru.md) · [All passports](index.en.md)

## What it is

The geometry of a signal between two moving points: the light-time equation for reception and
transmission with the relativistic (Shapiro) delay of any number of bodies, range and range
rate, round trips, aberration by the Lorentz transformation, light deflection in the PPN form,
right ascension and declination. On top of it, the measurements of a ground or space station
with the light time (range, range rate, RA/Dec, azimuth and elevation, one-way or two-way),
which carry exact derivatives in the state for [estimation](estimation.en.md). The dispersive
delays of the media are separate: [ionosphere and plasma](ionosphere-plasma.en.md),
[troposphere](troposphere.en.md); the moved place of a station is
[site displacement](site-displacement.en.md).

## Source

- I. I. Shapiro, *Phys. Rev. Lett.* 13, 789 (1964); T. D. Moyer, *Formulation for Observed and
  Computed Values of Deep Space Network Data Types*, DESCANSO Monograph 2 (2003), eq. 8-38,
  8-55 and §8.3.1 (the body taken at its position when the ray passes it).
- The aberration is the exact Lorentz transformation; the deflection is the PPN form of ERFA's
  `eraLd`.
- Written from the formulas; CSPICE and ERFA serve as references only.

## Interface

Header [`myb/observe.h`](../../include/myb/observe.h):

- functions `myb_light_time`, `myb_shapiro_delay`, `myb_round_trip`, `myb_aberration`,
  `myb_light_deflection`, `myb_direction_radec`;
- structures `myb_light_body`, `myb_light_model`, `myb_light_path`; the enumeration
  `myb_light_direction` (`MYB_LIGHT_RECEPTION`, `MYB_LIGHT_TRANSMISSION`);
- the participants are providers, `myb_state_function` of
  [`myb/kin.h`](../../include/myb/kin.h).

Header [`myb/estim.h`](../../include/myb/estim.h): `myb_measure_light_range`,
`myb_measure_light_range_rate`, `myb_measure_light_radec`, `myb_measure_light_azel` with the
context `myb_light_measure` (a `myb_station`).

## Inputs and outputs

- One inertial frame (the BCRS for the solar system) with its time coordinate (TDB); states in
  m and m/s; c in m/s; PPN γ (1 in general relativity).
- Light time: c (t − t_e) = |x_o(t) − x_T(t_e)| + c Σ Δt_B, solved by fixed-point iteration
  (contraction |v_T| / c) to the resolution of the time.
  Δt_B = (1 + γ) gm_B / c³ ln[(r₁ + r₂ + r₁₂ + (1 + γ) gm_B / c²) /
  (r₁ + r₂ − r₁₂ + (1 + γ) gm_B / c²)] — the gm terms inside the logarithm keep the
  significance near the limb. `myb_light_path` gives the target's time and state, the light
  time, the range, the delay included, the range rate at the observer (Newtonian light time)
  and the iterations.
- Round trip: the station receives at t; `total` = t − t_transmit, with the bounce time.
- Aberration: p′ = (p / g + (1 + p·b / (1 + 1/g)) b) / |…|, b = v / c, 1/g = √(1 − b²), exact
  in v / c. Deflection: p + (1 + γ) gm / (c² e_m) / max(q·(q + e), d_lim) (p × (e × q)),
  not normalised. RA in [0, 2π).
- Light measurements: what the station receives at the model time t left the object τ
  earlier, c τ = |x(t − τ) − s(t)|; two-way, the station sent it τ_up before,
  c τ_up = |x(t − τ) − s(t − τ − τ_up)|. Where the object was follows from its state at t by
  its motion about a point mass gm over τ, to the third order (the remainder is 10⁻¹⁴ m on a
  low orbit, 10⁻⁸ m at the distance of the Moon); gm = 0 moves it in a straight line. Range:
  one-way c τ, two-way c (τ + τ_up) / 2 (m); range rate in the time of reception (m/s);
  RA/Dec and azimuth (from north towards east) / elevation, with `aberration` as the moving
  station sees them. The derivatives in the state are exact, those of the implicit τ included
  (dual numbers through the fixed point).

## Domain of validity

- `myb_light_time`: `MYB_E_DOMAIN` (invalid model, e.g. c ≤ 0), `MYB_E_NOCONV` (no
  convergence in 32 iterations), the statuses of the providers.
- Light measurements: `MYB_E_DOMAIN` (c not positive, gm negative, the object at the station or
  at the origin, angles at the pole of their axes), `MYB_E_NOCONV` (no light time: an end
  faster than light), the statuses of the environment.
- The measurements' motion over τ is a point-mass Kepler arc: good for orbits about the origin
  of the environment, to the remainders above.

## Accuracy

| What | Reference | Result |
|---|---|---|
| Light time, Moon, Sun, Mars, Jupiter from the Earth, 1974–2026, 208 rows (CN, XCN) | CSPICE N0067 `spkezr` on DE440 | light time 4·10⁻¹² s; position 1.7 mm at 10¹² m; velocity 5·10⁻¹⁰ m/s |
| Lorentz aberration | SPICE's first-order stellar aberration | 2.5·10⁻⁹ rad (the (v/c)² SPICE leaves out) |
| Aberration | ERFA `eraAb` | bit for bit |
| Deflection | ERFA `eraLd` | 9·10⁻⁹ of the deflection (ERFA's SRS constant) |
| Delay, closed form | quadrature of the 1PN null line | 9·10⁻¹⁵ |
| Range rate | differences of light-time solutions | 3·10⁻¹¹ (relative) |
| Two-way Earth–Moon | 2 \|x_M − x_E\| / c at the bounce + 2 delays | 2.5·10⁻⁸ s (order (v/c)²) |
| Light measurements: orbits of 800, 36 000 and 380 000 km, one-way and two-way — range | `myb_light_time`, `myb_round_trip` on the exact Kepler motion | **≤ 1.4·10⁻⁸ m** |
| — range rate | differences of the ranges of the equation | 4·10⁻⁷ m/s (the limit of the difference) |
| — direction (with the station's aberration: against `myb_aberration`) | the same equation | 9·10⁻¹⁶ rad |
| — derivatives | central differences | within what a difference allows |
| An orbit from 90 two-way ranges | the orbit that made them | 2·10⁻⁹ m; read as geometric ranges, 122 m off |

The one-way delay of a ray grazing the Sun from 1 to 9.5 au is 141.6 µs, as
2 GM / c³ ln(4 r₁ r₂ / b²) gives. Tests: `test_observe`
([tests/unit/test_observe.c](../../tests/unit/test_observe.c); the CSPICE rows need
`de440.bsp`), `test_light_measure`
([tests/unit/test_light_measure.c](../../tests/unit/test_light_measure.c)).

## Cost

Not measured.

## Limitations

- The delaying body is taken at t − |x_B − x_o| / c, Moyer's approximation of its position
  when the ray passes; the deflection is first order in gm with a limiter near the limb.
- The range rate of `myb_light_path` uses the Newtonian light time.
- The light measurements carry no relativistic delay and no media delays — those are added as
  corrections of the observed value — and no time-scale conversion: everything is on the model
  time of the environment.
- The object's motion over the light time is a point-mass arc about the origin of the
  environment: other forces over τ are left out (the remainders above are for that motion).

## Data

None for the formulas. The CSPICE comparison needs `de440.bsp` (FetchData) and
`tests/data/observe/spice_lighttime.txt` (CSPICE N0067 outputs through spiceypy, hex floats).
