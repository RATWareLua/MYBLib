# Micrometeoroid impacts

**English** | [Русский](micrometeoroids.ru.md) · [All passports](index.en.md)

## What it is

Random impacts of interplanetary dust on an object, taken as a Poisson process of impulses.
The rate follows from the interplanetary flux model of Grün et al. (1985) and the mean
cross-section of the object; every impact draws a mass from the distribution of the flux, a
speed and a direction, and gives the object the momentum β m v. The impacts are applied to a
propagated orbit as jumps of the velocity.

## Source

- The flux: E. Grün, H. A. Zook, H. Fechtig, R. H. Giese, *Icarus* 62, 244 (1985), eq. A3 —
  the cumulative flux F(m) of particles of mass m or more at 1 au, on one side of a plate
  turned every way. The coefficients c4…c10, g4…g10 are printed in the header.
- The scaling with the distance from the Sun, (r / 1 au)^−p: the exponent is the caller's;
  the header cites r^−1.3 for the density of the zodiacal cloud (C. Leinert, I. Richter,
  E. Pitz, B. Planck, *A&A* 103, 177, 1981).
- Implemented from the paper; no third-party code, no data files.

## Interface

Header [`myb/meteoroid.h`](../../include/myb/meteoroid.h):

- functions `myb_meteoroid_flux`, `myb_meteoroid_rate`, `myb_meteoroid_mass`,
  `myb_meteoroid_next`, `myb_meteoroid_propagate`;
- structures `myb_meteoroid_options`, `myb_meteoroid_impact`;
- callbacks `myb_meteoroid_speed_function`, `myb_meteoroid_direction_function`;
- the generator `myb_rng` of [`myb/uncert.h`](../../include/myb/uncert.h) and the kick
  `myb_orbit_kick` of [`myb/orbit.h`](../../include/myb/orbit.h).

## Inputs and outputs

- Masses in kg (the paper's formula is in grams; the library converts), the flux in impacts
  per m² per s, the rate λ in impacts per s, the cross-section in m², speeds in m/s,
  impulses in N s, the change of velocity in m/s.
- The cross-section is the **mean** projected area (π R² for a sphere, a quarter of the
  surface of any convex body): the object takes 4 A F impacts a second.
- Rate between `mass_min` and `mass_max`: λ = 4 A (F(m_min) − F(m_max)) (r / 1 au)^−p.
- The direction is the way the particle moves, in the inertial axes of the model; isotropic
  by default. The speed is one value for all impacts (20 km/s by default, the mean impact speed
  Grün et al. take at 1 au) or drawn by the caller's function.
- An impact gives the momentum β m v along the direction of the particle (β ≥ 1 counts the
  ejecta; 1 by default) and the change of velocity β m v / M.
- Everything is drawn from xoshiro256** in a fixed order (wait, mass, speed, direction): one
  seed gives the same impacts bit for bit on every platform.

## Domain of validity

- Masses from about 10⁻¹⁸ g to 100 g (10⁻²¹ to 0.1 kg) — the range of the flux model.
- `myb_meteoroid_flux` returns NaN for a mass that is not positive and finite.
- Refused with `MYB_E_DOMAIN`: invalid options (non-positive cross-section or object mass,
  `mass_max ≤ mass_min`), u outside [0, 1) in `myb_meteoroid_mass`, a time that is not
  finite, a speed or direction drawn by the caller that is not valid, `t_end` behind the time
  reached.
- `MYB_E_NOCONV`: the inversion of the mass distribution does not converge, or
  `myb_meteoroid_propagate` meets far more impacts than λT allows
  (1000 + λT + 20 √(λT)).

## Accuracy

All in `test_meteoroid` ([tests/unit/test_meteoroid.c](../../tests/unit/test_meteoroid.c)):

| What | Reference | Result |
|---|---|---|
| F(m), 10⁻²¹ … 0.1 kg | eq. A3 evaluated with the C library's `pow` | 6·10⁻¹⁵ |
| Number of impacts, 1 m², 10⁻¹⁵ … 10⁻⁹ kg | λT = 50 000 | 50 085 (0.4σ) |
| 2000 runs with a mean of 20 | Poisson statistics | mean 20.06 (0.6σ), variance 19.26 (−1.2σ) |
| Inverse of the mass distribution | F itself | 3·10⁻¹⁵ of F |
| 50 000 masses drawn | Kolmogorov–Smirnov | √n·D = 1.17 (bound at 0.001: 1.95) |
| Total impulse (β = 2, 20 km/s, 1000 runs of ~100 impacts) | λT β v E[m], E[m] by quadrature of F | 0.3σ; a component 0.9σ from zero; its variance −0.5σ |
| 310 impacts applied to a free motion with `myb_orbit_kick` | free motion with the same kicks | 7·10⁻¹⁶ |

The rate scales with the distance exactly as (r / 1 au)^−p (2^−1.3 checked to 10⁻¹⁵).

## Cost

Not measured.

## Limitations

- Only the interplanetary background of Grün et al. 1985 is built in: no orbital debris, no
  meteor streams, no apex enhancement — a stream or an apex enters only through the caller's
  direction and speed functions.
- The scaling with the distance is a single power law; the dependence of what the object meets
  on its own speed through the cloud is not modelled (the header says so).
- An impact changes the velocity of the centre of mass only: no torque, no change of the
  attitude, no erosion or damage.
- A later call of `myb_meteoroid_propagate` draws from its own start: it has the same
  statistics, not the same impacts, as one call over both intervals.

## Data

None: the coefficients of the flux are constants of the code.
