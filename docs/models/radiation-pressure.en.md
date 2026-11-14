# Radiation pressure and shadows

**English** | [Русский](radiation-pressure.ru.md) · [All passports](index.en.md)

## What it is

The push of the light of a source on an object treated as a sphere (cannonball) or as flat
plates, scaled by the part of the source's disk that the object still sees past the bodies in
front of it. The shadow is computed on the sky sphere with the angular sizes of both disks,
quadratic limb darkening, several occulters at once and, optionally, the refraction and
extinction in the atmosphere of an occulter. For a moving object the terms of the first order in
v/c (aberration and Doppler shift) give the Poynting–Robertson drag. Sunlight reflected and heat
emitted by a planet follow the model of Knocke, Ries and Tapley. Pure functions in `force.h`,
ready-made terms of a force model with analytic Jacobians in `env.h`.

## Source

- Sphere and plates: D. Vallado, *Fundamentals of Astrodynamics and Applications* (2013),
  §9.6.3; O. Montenbruck, E. Gill, *Satellite Orbits* (2000), eq. 3.75. Thermal re-emission of
  the absorbed light as a Lambertian emitter on the lit side.
- Moving object: J. A. Burns, P. L. Lamy, S. Soter, *Icarus* 40, 1 (1979), eq. 2. The decay of
  a and e is checked against the averaged equations of Wyatt and Whipple (1950).
- Shadow: exact geometry of the disks on the sky sphere (the occulter's disk is a spherical cap)
  and the brightness I(μ) = 1 − u1 (1 − μ) − u2 (1 − μ)², integrated over annuli of the source's
  disk with exact occulted arcs and Gauss–Legendre rules between the kinks of the integrand.
  Several occulters: every annulus loses the union of the arcs they cover.
- Refraction: the theory of the illumination of the eclipsed Moon (F. Link, *Eclipse Phenomena
  in Astronomy*, 1969) as it enters the radiation pressure on satellites (D. Vokrouhlický,
  P. Farinella, F. Mignard, *A&A* 280, 295, 1993); bending (n − 1) √(2π (R + h) / H), the first
  order in the refractivity.
- Planet: G. Knocke, J. Ries, B. Tapley, AIAA 88-4292 (1988) — zonal albedo and emissivity
  a0 + a1 P1(sin φ) + a2 P2(sin φ), a1 seasonal.
- Implemented from the papers; no third-party code, no data files.

## Interface

Header [`myb/force.h`](../../include/myb/force.h):

- functions `myb_radiation_flux`, `myb_force_radiation_sphere`, `myb_force_radiation_moving`,
  `myb_force_radiation_plates`, `myb_shadow_fraction`, `myb_shadow_margin`,
  `myb_shadow_fraction_many`, `myb_refraction_bending`, `myb_shadow_fraction_refracted`,
  `myb_force_planet_radiation`, `myb_force_planet_radiation_gradient`;
- structures `myb_plate`, `myb_refraction`, `myb_planet_radiation`.

Header [`myb/env.h`](../../include/myb/env.h):

- terms `myb_term_radiation`, `myb_term_radiation_jacobian` (context `myb_radiation_term`),
  `myb_term_planet_radiation`, `myb_term_planet_radiation_jacobian` (context
  `myb_planet_radiation_term`); the object `myb_object` (mass, area, cr, a mass function).

The edges of the shadows as breaks of the integrator's step (`myb_term_radiation_breaks`,
`myb_eclipse_event` of [`myb/events.h`](../../include/myb/events.h)) are in
[events and geometry](events-geometry.en.md); plates turned by an attitude law, meshes and the
recoil of emitted radiation in [radiation on surfaces](radiation-surfaces.en.md).

## Inputs and outputs

- SI throughout: m, m/s, W, W m⁻²; the pressure is flux / c (N m⁻²); c is a parameter of the
  call or of the term. The functions of `force.h` add to `a[3]`; the terms write it.
- u is the unit direction of propagation, from the source to the object. Sphere:
  a += P cr (A/m) u with cr = 1 + reflectivity.
- Plate: area, outward unit normal, fractions reflected specularly and diffusely (the rest is
  absorbed); `thermal` ≠ 0 re-emits the absorbed part at once from the lit side. Only plates lit
  from the front (n·u < 0) count; plates, u and a share one frame (usually the body frame).
- Moving object: v is the velocity relative to the source;
  a += P cr (A/m) [(1 − u·v/c) u − v/c].
- Shadow: positions of the source and of the occulter relative to the object, their radii, limb
  darkening (u1, u2) or NULL for a uniform disk. Result 0…1; 0 inside the occulter, 1 when the
  occulter is not in front (beyond the source or behind the object).
- `myb_shadow_margin` (rad): c − (a + b) with `inner` = 0 (zero at the edge of the penumbra),
  c − |b − a| otherwise (zero at the edge of the umbra or of the annular phase); a, b the angular
  radii, c the angle between the centres; positive towards the light, 1 where the occulter is not
  in front, −1 inside it.
- Atmosphere of an occulter: refractivity n − 1 and vertical optical depth at the surface, scale
  height H of an exponential atmosphere; the source is taken far beyond the occulter.
- Planet: object and source in the planet's body-fixed frame (z the rotation axis), the solar
  flux at the planet, the seasonal phase ω (t − t0), cr A/m of the object as a sphere; `albedo` =
  (a0, c0, c1, c2, a2), `emissivity` likewise. The visible cap is a central disk and `rings`
  rings of 6k elements of equal projected weight. The gradient is written as dadr[3 i + j].
- Term of the source: a = ν L / (4π d² c) cr (A/m) u. ν is the product of the separate
  fractions; where two or more airless occulters cover parts of the disk at once, the union of
  what they cover (one occulter on the disk: the product, to the bit). An occulter with an
  atmosphere (`atmospheres`, scale height > 0) enters as a factor. `moving` ≠ 0 adds the v/c
  terms. The planet term turns the state by the planet's orientation and takes the flux from the
  source's luminosity.

## Domain of validity

- Plates: `MYB_E_DOMAIN` for a mass ≤ 0, a negative area, a normal not unit to 10⁻¹², fractions
  outside [0, 1] or summing to more than 1.
- Planet radiation: `MYB_E_DOMAIN` for an object inside the planet, `rings` outside 1…50,
  invalid values.
- Terms: `MYB_E_DOMAIN` for the object at the centre of a body, invalid values of the structure
  (a mass or c that is not positive), an occulter that is the source, occulters with a source
  without a radius; the statuses of the providers.
- `myb_shadow_fraction_many` takes at most 8 occulters: more is a broken precondition and the
  library stops on its check (no status is returned).
- Refraction is of the first order in the refractivity: within a few per cent at the surface of
  the Earth, where a grazing ray is bent by 1.1°. Rays whose tangent point lies under the surface,
  or under the level where the air bends light as fast as the planet curves, do not come out.
- The moving term is of the first order in v/c.

## Accuracy

Tests [test_force](../../tests/unit/test_force.c), [test_shadow](../../tests/unit/test_shadow.c),
[test_env](../../tests/unit/test_env.c), [test_geometry](../../tests/unit/test_geometry.c),
[test_nongrav](../../tests/unit/test_nongrav.c), [test_lageos](../../tests/unit/test_lageos.c):

| What | Reference | Result |
|---|---|---|
| Plate: absorber, mirror, diffuse, thermal (`test_force`) | P A/m, twice that, (1 + 2/3) of it | to rounding |
| Sphere of 320 000 plates, diffuse or specular 0.6 (`test_force`) | cr = 1 + 4ρ/9 and 1 | within 2·10⁻⁵ (the test's bound) |
| Shadow fraction, uniform and limb-darkened: penumbra at 420 km, at LAGEOS, from the distance of the Moon, annular eclipse (`test_force`) | exact area of two spherical caps; 40-digit quadrature (mpmath) | 1.1·10⁻¹⁴ |
| Zero of `myb_shadow_margin` 40 000 km behind the Earth (`test_geometry`) | where the fraction leaves 1 / reaches 0 | within 1 m of the edge |
| Two occulters overlapping on the Sun's disk (`test_shadow`) | count of points of the disk | 1.7·10⁻⁶ (the product of the fractions: 7·10⁻²) |
| Two occulters apart; one occulter given as many (`test_shadow`) | f1 + f2 − 1; `myb_shadow_fraction` | 10⁻¹¹ (bound); bit for bit |
| Grazing ray at the Earth, n − 1 = 2.77·10⁻⁴, H = 8430 m (`test_shadow`) | two horizontal refractions of 34.5′ | 65.7′ |
| Next to no air (`test_shadow`) | airless fraction | 3·10⁻¹⁰ |
| Clear air without extinction, plane 40 000 km behind the Earth (`test_shadow`) | flux of the rays that meet the ground | 1.000559 against 1.000554 cross-sections |
| 6000 km from the axis, inside the airless umbra (`test_shadow`) | — | clear air leaves 2.2 % of the light |
| Uniform Lambertian planet (`test_force`) | M (R/r)², radial | 0.4/rings² |
| `myb_force_planet_radiation_gradient` (`test_shadow`) | central differences of the acceleration | 4.5·10⁻⁴, the accuracy of the differences over the rings (bound 5·10⁻³) |
| Every term (`test_env`) | the functions of `force.h` called by hand with the providers' values | bit for bit |
| Term with two bodies on the disk; with an atmosphere (`test_shadow`) | `myb_shadow_fraction_many`; `myb_shadow_fraction_refracted` | 10⁻¹⁴ |
| `myb_term_radiation_jacobian` in full light | central differences | 6·10⁻¹¹ |
| Poynting–Robertson, β = 0.2, e = 0.3, 20 revolutions: a falls by 7.272665·10⁻³ au, e by 2.190050·10⁻³ (`test_nongrav`) | averaged equations of Wyatt and Whipple, integrated | the same seven digits |
| Jacobian of the moving term (`test_nongrav`) | central differences | 2·10⁻¹⁰ (position), 10⁻⁹ (velocity) |
| LAGEOS-1 over a week, sunlight with the penumbrae of the Earth and the Moon and Knocke's albedo/IR among the forces (`test_lageos`) | ILRS combined orbit | 1.6 cm RMS for the whole model |

In LAGEOS the shadows are worth 1 mm of the RMS (2.1 cm without them against 2.0 cm, measured
before the secondary waves of the ocean tides went in). With the edges of the shadows declared
as breaks the residual of the first guess, 25.8887 m, is the same at ε = 10⁻¹⁰, 10⁻¹¹, 10⁻¹².

## Cost

Not measured for the functions. Ending the integrator's steps at the edges of the shadows
(withdrawing the step) costs +26 % of the time of the LAGEOS fit (docs/PLAN.en.md, "The user's
decision of 2026-11-11").

## Limitations

- The source and the occulters are spheres: no oblateness of the Earth in its shadow.
- Plates cast no shadows on each other (a mesh with self-shadowing: radiation on surfaces); the
  thermal re-emission is instantaneous; the optical fractions are constants.
- The atmosphere of an occulter is exponential, with one scale height and one optical depth, and
  only to the first order in the refractivity; the union of arcs is for airless occulters.
- In a penumbra, where the fraction changes over tens of kilometres, the Jacobian of the term is
  central differences of the term itself (steps of 6·10⁻⁶ of the distance from the origin).
- The planet's light: Lambertian reflection and thermal emission spread uniformly over the
  sphere, the zonal coefficients of Knocke's model, the object as a sphere. The rotation of the
  gradient into the inertial axes (`myb_term_planet_radiation_jacobian`) is not checked by a
  test of its own; test_lageos leaves that gradient out.

## Data

None: luminosities, radii, limb darkening and the coefficients of a planet are the caller's.
