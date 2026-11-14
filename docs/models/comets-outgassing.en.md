# Comets and outgassing

**English** | [Русский](comets-outgassing.ru.md) · [All passports](index.en.md)

## What it is

The non-gravitational acceleration of comets and asteroids in the form the orbit fits of the
catalogues use: three constants along the radius from the Sun, across it in the orbital plane
towards the motion, and along the orbital angular momentum, each times one function of the
distance,
a = g(r) (A1 r̂ + A2 t̂ + A3 n̂), g(r) = α (r/r0)⁻ᵐ (1 + (r/r0)ⁿ)⁻ᵏ.
Two forms of g are ready: the sublimation of water ice of Marsden, Sekanina and Yeomans and a
power law (with A2 alone and exponent 2 it is the Yarkovsky drift as fitted to asteroids). The
activity may lag the heating: g is then taken where the body was a given time earlier on its
orbit.

## Source

- B. G. Marsden, Z. Sekanina, D. K. Yeomans, *AJ* 78, 211 (1973): r0 = 2.808 au, m = 2.15,
  n = 5.093, k = 4.6142, α such that g(1 au) = 1 (α = 0.1112620426, the ten digits of the JPL
  catalogue).
- Power law g = (r0/r)^exponent; with A2 alone and the exponent 2 the drift is
  da/dt = 2 A2 (1 − e²) (r0/p)² / n, p the semi-latus rectum (D. Farnocchia et al., *Icarus* 224,
  1, 2013).
- Delayed outgassing: D. K. Yeomans, P. W. Chodas, *AJ* 98, 1083 (1989).
- Implemented from the papers; no third-party code, no data files.

## Interface

Header [`myb/force.h`](../../include/myb/force.h):

- structure `myb_nongrav` (A1, A2, A3, r0, α, m, n, k);
- functions `myb_nongrav_marsden`, `myb_nongrav_power`, `myb_nongrav_g`, `myb_force_nongrav`,
  `myb_force_nongrav_delayed`.

Header [`myb/env.h`](../../include/myb/env.h): `myb_term_nongrav`, `myb_term_nongrav_jacobian`
(context `myb_nongrav_term`: the source body, the model, the delay). The delay propagates the
two-body orbit with `myb_kepler_propagate`
([Kepler and regularisation](kepler-regularization.en.md)).
The Yarkovsky effect from the thermal properties of a body is in
[Yarkovsky and YORP](yarkovsky-yorp.en.md).

## Inputs and outputs

- A1, A2, A3 in m/s² where g = 1 (at 1 au for Marsden's g, at r0 for the power law); r0 in m.
- r and v of the body relative to the Sun (m, m/s); `myb_force_nongrav` adds to `a[3]`.
- Directions: r̂ from the Sun; t̂ across r̂ in the orbital plane, towards the motion; n̂ along
  r × v.
- `myb_nongrav_g(model, r)`: g at the distance r; NaN at r = 0.
- Delay (s): g is taken at the distance the body had `delay` seconds earlier on the two-body
  orbit of GM = gm through r, v; the directions are those of the moment. A positive delay puts
  the peak of the outgassing after the perihelion. In the term the GM is that of the source body.

## Domain of validity

- `MYB_E_DOMAIN`: r = 0, r parallel to v (no orbital plane), values not finite; with a delay also
  a GM that is not positive or a delay that is not finite, and the errors of
  `myb_kepler_propagate`.
- Marsden's g describes water ice: about r0 = 2.808 au it turns down, and far beyond it falls as
  r^−(m + n k).

## Accuracy

All in [test_nongrav](../../tests/unit/test_nongrav.c):

| What | Reference | Result |
|---|---|---|
| Marsden's g at 1 au | 1 | 1 + 2·10⁻⁹ (the ten-digit α) |
| Marsden's g at r0 | α 2⁻ᵏ | 10⁻¹⁴ |
| Fall between 20 and 40 au | r^−(m + n k) | 10⁻³ in the exponent |
| Power law of exponent 2 at 2 r0 | 1/4 | 10⁻¹⁶ |
| Components of the force | projections on r̂, t̂, n̂ | 10⁻¹² |
| A2 = 10⁻¹⁰ m/s², a = 1.3 au, e = 0.3, 10 revolutions | 2 A2 (1 − e²) (r0/p)² / n = 9.681647·10⁻⁴ m/s | 9.681642·10⁻⁴ m/s (5·10⁻⁷) |
| A comet with A1 = 3·10⁻⁹, A2 = 4·10⁻¹⁰ m/s² and Marsden's g, 10 revolutions | Gauss's equation averaged over the orbit by a quadrature | 5·10⁻⁶ |
| Delay 2·10⁶ s, the body 2·10⁶ s after the perihelion (a = 1.3 au, e = 0.5) | undelayed at the perihelion, 2.579·10⁻⁸ m/s² | 10⁻⁹ (the test's bound) |
| Delay 0 | `myb_force_nongrav` | bit for bit |
| `myb_term_nongrav_jacobian` | central differences | 1.4·10⁻¹⁰ and 5·10⁻¹¹ |

[test_spacecraft](../../tests/unit/test_spacecraft.c) fits the radial law A1 (1 au/r)² of JPL's
solution #11 to the Tesla Roadster over 7 years against Horizons (RMS 306 m,
A1 = 5.95·10⁻⁸ m/s²), but with a term written in the test, not with `myb_term_nongrav`.

## Cost

Not measured.

## Limitations

- The model is the empirical one of the orbit fits: one function of the distance for all three
  components, constant coefficients; no jets, no rotation of the nucleus, no asymmetry about the
  perihelion other than the delay.
- Only g lags: the delay uses the two-body orbit about the Sun, and the directions of the force
  stay those of the moment.
- With a delay the Jacobian of the term is central differences of the term.
- The ready-made g is that of water ice; other volatiles need the caller's r0, m, n, k, α.

## Data

None: A1, A2, A3 and the form of g are the caller's (e.g. from a catalogue orbit solution).
