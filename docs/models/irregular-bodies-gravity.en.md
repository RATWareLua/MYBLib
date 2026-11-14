# Gravity of irregular bodies: polyhedra, mascons, ellipsoidal harmonics

**English** | [Русский](irregular-bodies-gravity.ru.md) · [All passports](index.en.md)

## What it is

Three forms of the field of a body far from a sphere — an asteroid, a comet nucleus, a small
moon. The exact field of a homogeneous closed polyhedron, valid inside, outside and on the
surface; a sum of point masses (mascons); and the exterior expansion in ellipsoidal harmonics
about a triaxial reference ellipsoid, which converges down to the ellipsoid itself — inside the
Brillouin sphere of an elongated body, where the spherical-harmonic series diverges.

## Source

- The polyhedron: R. A. Werner, D. J. Scheeres, *Celest. Mech. Dyn. Astron.* 65, 313 (1997);
  the Laplacian from the solid angles.
- Ellipsoidal harmonics: R. Garmier, J.-P. Barriot, *Celest. Mech. Dyn. Astron.* 79, 235
  (2001). The Lamé functions of the four classes K, L, M, N come from three-term recurrences
  derived by substitution into Lamé's equation (the eigenvalues as roots of a continuant, then
  the roots of the polynomial Q in (0, k²) by Stieltjes). The exterior term
  F(λ) E(μ) E(ν) is evaluated as R(λ) H(x, y, z), with H the Cartesian solid-harmonic
  polynomial of Niven and R(λ) an integral with a bounded integrand: the gradient has no
  singular set of the ellipsoidal coordinates.
- Mascons: a compensated sum.
- Implemented from the papers; no third-party code.

## Interface

Header [`myb/grav.h`](../../include/myb/grav.h):

- the polyhedron: `myb_polyhedron_create`, `myb_polyhedron_destroy`, `myb_polyhedron_volume`,
  `myb_polyhedron_eval`;
- mascons: `myb_mascons_eval`;
- ellipsoidal harmonics: `myb_ellipsoidal_field_create`, `myb_ellipsoidal_field_destroy`,
  `myb_ellipsoidal_field_degree`, `myb_ellipsoidal_field_set_coefficients`,
  `myb_ellipsoidal_field_eval`, `myb_ellipsoidal_field_basis`, `myb_ellipsoidal_lame`.

There is no ready-made force term for these fields in `myb/env.h`: the caller wraps them in a
term function. The same mesh serves the geometry of [events and geometry](events-geometry.en.md)
and the ground of [contact and landing](contact-landing.en.md).

## Inputs and outputs

- The polyhedron: vertices in m, faces as three vertex indices, counter-clockwise seen from
  outside; the density as G ρ (`g_rho`, s⁻²); r in m in the frame of the mesh. The potential is
  positive (m²/s²), the acceleration in m/s², the Laplacian −4π G ρ inside, 0 outside and, on
  the surface, that of the local solid angle; any output may be NULL. `myb_polyhedron_volume`
  gives the volume (m³) and the centre of mass (m).
- Mascons: `gm[k]` (m³/s²) at `positions[3k .. 3k+2]` (m, body frame); potential
  Σ gm / |p − r| and acceleration, summed with compensation in array order.
- Ellipsoidal harmonics: semi-axes a > b > c > 0 (m), h² = a² − b², k² = a² − c²;
  U = (GM/a) Σₙ Σₚ αₙₚ [Fₙₚ(λ) / Fₙₚ(a)] Eₙₚ(μ) Eₙₚ(ν) / a²ⁿ, λ ≥ a ≥ μ ≥ … the ellipsoidal
  coordinates; (degree + 1)² coefficients αₙₚ indexed n² + p (NULL: zeros); within a degree the
  2n + 1 functions are ordered by class K, L, M, N. Eₙₚ has a monic polynomial factor (mⁿ).
  Positions in the frame of the ellipsoid's axes. `myb_ellipsoidal_field_basis` gives the
  potential of each unit coefficient divided by GM/a, for estimating α from data.

## Domain of validity

- The polyhedron: a closed, consistently oriented mesh — every edge shared by exactly two faces
  in opposite directions. `MYB_E_DOMAIN`: an index out of range, a degenerate face, an open
  mesh, a face or the whole mesh turned inwards, `g_rho` negative, a vertex or r not finite.
  Valid inside, outside and on the surface (edge and face terms take their limits). The terms
  grow like r² while the potential falls like 1/r: the relative rounding is about
  10⁻¹⁶ (r / size)², so far from the body use its spherical-harmonic expansion.
- Mascons: `MYB_E_DOMAIN` at a mascon or for r not finite.
- Ellipsoidal harmonics: on and outside the reference ellipsoid (`MYB_E_DOMAIN` inside);
  degree ≤ 20 (21 is refused); axes strictly ordered (a = b is refused); `MYB_E_NOCONV` when the
  Lamé functions are not resolved.

## Accuracy

Tests: [`test_grav`](../../tests/unit/test_grav.c),
[`test_ellipsoidal`](../../tests/unit/test_ellipsoidal.c).

| What | Reference | Result |
|---|---|---|
| Polyhedron, 6 points inside and outside a prism | closed form of the rectangular prism (D. Nagy, G. Papp, J. Benedek 2000) in 40 digits (`test_grav`) | potential 3·10⁻¹⁶, acceleration 7.6·10⁻¹⁶, times max(1, (r/size)²) |
| Laplacian | −4π G ρ inside, 0 outside (`test_grav`) | 5·10⁻¹³ (bound) |
| On an edge, a vertex, a face diagonal, a face | the field 10⁻⁹ m away (`test_grav`) | continuous: 10⁻⁸ in the potential, 10⁻⁶ in the acceleration (bounds) |
| 12 mascons at 3R | their own spherical-harmonic series to degree 40, with an independent Legendre recursion (`test_grav`) | 1.2·10⁻¹⁵ |
| Every ellipsoidal basis term to degree 16 | harmonicity, a fourth-order finite-difference Laplacian (`test_ellipsoidal`) | below 10⁻⁶ of the second derivatives (bound) |
| Niven's form on the reference ellipsoid | E(μ) E(ν) / a²ⁿ with the coordinates from the confocal cubic (`test_ellipsoidal`) | 3·10⁻¹² |
| Ellipsoidal acceleration | differences of the potential (`test_ellipsoidal`) | 10⁻⁸ (bound) |
| A 4 × 2.4 × 1.6 box inside its Brillouin sphere, outside the ellipsoid, degree 4 → 16 | the exact field of the box as a polyhedron (`test_ellipsoidal`) | ellipsoidal 8.9·10⁻³ → 6.5·10⁻⁵; the exact spherical series of the same box diverges, 3.8·10⁻² → 0.26 |

## Cost

Not measured.

## Limitations

- The polyhedron has one density; an inhomogeneous body needs several meshes or added mascons.
  It gives the potential, the acceleration and the Laplacian, not the gradient tensor.
- Far from the body the polyhedron loses digits as (r / size)².
- Ellipsoidal harmonics are exterior only. Harmonicity is checked to degree 16 while degree 20
  is allowed: degrees 17 to 20 are not measured. On the box the corners, close to the reference
  ellipsoid, set the slow rate of convergence.
- No ready-made term and no analytic Jacobian of these fields in `myb/env.h`; no estimator of
  α beyond the basis functions.

## Data

None: the mesh, the masses and the coefficients are the caller's.
