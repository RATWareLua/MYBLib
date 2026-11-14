# Lagged tides: tidal dissipation

**English** | [Русский](tides-lagged.ru.md) · [All passports](index.en.md)

## What it is

Tidal dissipation as a delay of the tidal bulge. Two models. Mignard's constant time lag: a
point mass raises a tide on a rotating body, the bulge is where the raiser was a time Δt ago,
and the raiser feels a force across the radius that changes its orbit — with the analytic
Jacobian and the torque on the body (what the orbit gains, the rotation loses). And, in a
system of bodies, the tide by band of the JPL ephemerides: the long-period, diurnal and
semi-diurnal bands each with its Love number, its orbital lag and its rotational lag, acting on
a third body or on the raiser itself, with the reaction on the body.

## Source

- F. Mignard, *Moon and Planets* 20, 301 (1979).
- W. M. Folkner, J. G. Williams, D. H. Boggs, R. S. Park, P. Kuchynka, IPN Progress Report
  42-196 (2014), eqs. 31–32: the tide by band with orbital and rotational lags, the model of
  DE440.
- Implemented from the papers; no third-party code.

## Interface

- Header [`myb/force.h`](../../include/myb/force.h): `myb_force_tidal_lag`,
  `myb_torque_tidal_lag`.
- Header [`myb/env.h`](../../include/myb/env.h): the structure `myb_tidal_lag_term` with
  `myb_term_tidal_lag`, `myb_term_tidal_lag_jacobian`.
- Header [`myb/system.h`](../../include/myb/system.h): the structure `myb_system_tide`,
  `myb_system_set_tide`, `myb_system_add_tide` — see [system of bodies](system-of-bodies.en.md).

The tide that distorts a body whose rotation is integrated (`myb_system_add_distortion`) and the
lunar librations: [figure–figure and rotation](figure-figure-rotation.en.md) and
[system of bodies](system-of-bodies.en.md). The tides without lag of a field:
[IERS 2010 tides](tides-iers.en.md).

## Inputs and outputs

- Mignard: the body (GM, Love number k₂, radius R, angular velocity ω in rad/s in the axes of
  r), the raiser's GM, r and v of the raiser relative to the body (m, m/s), the lag Δt (s):
  a −= 3 k₂ gm_raiser (1 + gm_raiser/gm_body) R⁵/r⁸ { [r] + Δt [2 (r·v) r/r² + r × ω + v] };
  the term [r], the tide without delay that the tides of a field already give, only with
  `with_static`. For a circular orbit in the equator da/dt = 6 k₂ Δt (m/M)(R/a)⁵ n a (ω − n).
  The tide the body raises on the raiser is the same call with the roles exchanged.
- `myb_torque_tidal_lag`: the torque of the same tide on the body, per unit of its mass (the
  unit of `myb/attitude.h`), in the axes of r; written, not added.
- The term: the motion of a moon about its planet; the rotation of the body from its
  orientation, `gm_object` the GM of the object, `k2`, `lag` (s), `with_static`.
- The tide by band: a body T at r_T raises, at r,
  U = G m_T R⁵/(2 |r|⁵) Σⱼ k₂[j] Wⱼ(r, r*ⱼ) / |r*ⱼ|⁵ with
  W₀ = (3z² − |r|²)(3z*² − |r*|²)/2, W₁ = 6 z z* (ρ·ρ*), W₂ = 3 (2 (ρ·ρ*)² − ρ² ρ*²)/2
  (z along the pole — the third axis of the body's frame, ρ across it); for equal bands this is
  k₂ G m_T R⁵ P₂(cos) / (|r|³ |r*|³). The raiser is taken where it was:
  r*ⱼ = R_pole(rotation_rate · rotation_lag[j]) r_T(t − orbital_lag[j]), the earlier position
  from a Taylor series to the third order of the motion at t (the jerk of the two bodies alone;
  no history is kept). The point is accelerated by the gradient and the body takes the
  reaction; the raiser and the point may be one.
- `myb_system_tide`: `radius` (m, > 0), `k2[3]`, `orbital_lag[3]` (s, may be negative),
  `rotation_lag[3]` (s; [0] unused — the long-period band does not turn), `rotation_rate`
  (rad/s).

## Domain of validity

- Mignard: first order in the lag; `MYB_E_DOMAIN` for r = 0 and invalid input.
- The tide by band: the body needs a frame given by an orientation by the time the
  accelerations are asked for. `myb_system_set_tide` refuses (`MYB_E_DOMAIN`) a body whose
  rotation is integrated, a radius that is not positive, values that are not finite. Up to
  `MYB_SYSTEM_MAX_TIDES` = 64 triples (body, raiser, point).
- The earlier position from the series is measured for a lag of 900 s on an orbit of 27 days:
  the lags must be short against the orbital period.

## Accuracy

Tests: [`test_nongrav`](../../tests/unit/test_nongrav.c),
[`test_system`](../../tests/unit/test_system.c),
[`test_de440_system`](../../tests/unit/test_de440_system.c).

| What | Reference | Result |
|---|---|---|
| The Moon on the Earth, k₂ = 0.3, 600 s | da/dt of the formula (`test_nongrav`) | 3.79 cm a year |
| The acceleration across the radius | da/dt = 2 T/n of the formula (`test_nongrav`) | 10⁻¹⁴ (bound) |
| A synchronous orbit; the static part | zero; 3 k₂ gm (1 + m/M) R⁵/a⁷ along the radius (`test_nongrav`) | exactly zero; 10⁻¹⁴ (bound) |
| Ten revolutions with a lag of 6·10⁵ s | the drift of the formula (`test_nongrav`) | 28.348358 m against 28.348364 m; a retrograde spin pulls the orbit in, as the formula says (10⁻⁵ bound) |
| The torque | the angular momentum the orbit gains (`test_nongrav`) | 10⁻¹⁵ (bound); the Moon brakes the Earth |
| Jacobian of the term, with and without the static part | central differences (`test_nongrav`) | 2·10⁻⁹ |
| Tide by band, equal bands, no lag | the gradient of k₂ G m R⁵ P₂/(r³ r*³) by differences (`test_system`) | 7·10⁻¹¹ |
| Three different bands | differences of the potential band by band (`test_system`) | 10⁻⁸ (bound) |
| The raiser 900 s earlier | Kepler motion of the pair (`test_system`) | 3·10⁻¹¹ |
| The lagging part, lag of 60 s | Mignard's first-order formula (`test_system`) | 5·10⁻³ at a lag angle of 4·10⁻³ |
| DE440 again: the tides of the Earth from the Moon and the Sun on the Moon with the five lags of DE440 | DE440 (`test_de440_system`) | the Moon 0.011 m after 30 days, 0.088 m after a year, 0.19 m after 5 years; 54.3 m after 5 years without the tides of the Earth |

## Cost

Not measured.

## Limitations

- Mignard: one Love number of degree 2, a constant time lag, the first order in the lag, a
  point-mass raiser.
- The tide by band: degree 2 only, three bands about the pole of the body, a constant rotation
  rate; earlier positions from the series at t, the jerk of the two bodies alone.
- The term gives only the tide that the object raises on the body itself (a moon about its
  planet).

## Data

None: the Love numbers and the lags are the caller's (`test_de440_system` reads `K2E0`…`K2E2`,
`TAUE0`…`TAUE2`, `TAUER1`, `TAUER2` of `header.440`).
