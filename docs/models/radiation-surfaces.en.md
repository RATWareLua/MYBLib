# Radiation on surfaces

**English** | [Русский](radiation-surfaces.ru.md) · [All passports](index.en.md)

## What it is

Light on an object that is more than a sphere: flat plates turned by an attitude law, or a closed
triangular mesh whose parts cast shadows on each other, with the force and the torque. The same
header carries the recoil of radiation that the object sends out itself (the heat of a
radioisotope generator, the beam of an antenna), a mesh whose faces re-emit the absorbed heat
late, as a surface with thermal inertia does on a turning body, and the light on a body of a
system of bodies with the torque of the moment. The terms of a force model have analytic
Jacobians, the attitude law included.

## Source

- Plates: the formula of [radiation pressure](radiation-pressure.en.md) (Vallado 2013, §9.6.3;
  Montenbruck & Gill 2000, eq. 3.75), applied in the body frame given by the attitude law.
- Mesh: every face is such a plate; with self-shadowing a face is cut into k² triangles and a
  part is lit only if the ray from its centre to the source meets nothing of the mesh
  (`myb_polyhedron_ray`, a tree of bounding boxes: faces halved along the long side of their box
  down to four in a leaf). With a purely thermal surface this is the zero-conductivity model of
  D. P. Rubincam, *Icarus* 148, 2 (2000) — its YORP use is in
  [Yarkovsky and YORP](yarkovsky-yorp.en.md).
- Thermal inertia on a turning mesh: the answer of a half-space to a periodic heating,
  1 / (1 + θ √k (1 + i) / 2) for the harmonic k of the rotation (the transfer function of
  `myb_force_yorp`; D. Čapek, D. Vokrouhlický, *Icarus* 172, 526, 2004).
- Recoil: a = −η P / (m c) d, η = 1 for a beam, 2/3 for a flat Lambertian radiator.
- Implemented by the library; no third-party code, no data files.

## Interface

Header [`myb/surface.h`](../../include/myb/surface.h):

- functions `myb_force_radiation_mesh`, `myb_force_radiation_mesh_lagged`;
- attitude laws `myb_attitude_fixed`, `myb_attitude_orbital`, `myb_attitude_sun_nadir`,
  `myb_attitude_spin` with `myb_orbital_attitude`, `myb_sun_attitude`, `myb_spin_attitude`; the
  callback `myb_attitude_function`;
- terms `myb_term_surface_radiation`, `myb_term_surface_radiation_jacobian`
  (`myb_surface_radiation_term`), `myb_term_recoil`, `myb_term_recoil_jacobian`
  (`myb_recoil_term`).

Header [`myb/force.h`](../../include/myb/force.h): `myb_force_emission`,
`myb_force_radiation_plates`, `myb_plate`. Header [`myb/env.h`](../../include/myb/env.h): the
mass function of `myb_object` (a mass that changes). Header
[`myb/system.h`](../../include/myb/system.h): `myb_system_light`, `myb_system_light_action`,
`myb_system_add_action` — see [the system of bodies](system-of-bodies.en.md). The mesh is a
`myb_polyhedron` of [`myb/grav.h`](../../include/myb/grav.h).

## Inputs and outputs

- SI; pressure = flux / c (N m⁻²); force in N, torque in N m about the origin of the mesh; the
  mesh functions write `force` and `torque` (either may be NULL).
- An attitude law gives the matrix M (rows) with body = M inertial, at the model time t and the
  state (x, v) relative to the origin of the environment. Orbital: x from the body to the object,
  z along the orbital angular momentum. Sun–nadir (yaw steering): z towards the body, y across
  the plane of the body, the object and the source, x in that plane on the side of the source.
  Spin: body = R3(phase + rate t) M0 inertial, z the axis.
- Mesh: closed (several closed parts may make one mesh); every face takes the light with the
  `specular`, `diffuse`, `thermal` of `surface` (its area and normal are not used); u is the
  direction of propagation in the frame of the mesh. `shadows` = 0: every face turned to the
  source is lit; `shadows` = k ≥ 1: k² parts a face, one ray each.
- Lagged mesh: `axis` (unit, in the frame of the mesh, in the sense of the rotation), θ the
  thermal parameter at the rotation frequency, `phases` steps of a turn from the moment on;
  the absorbed power of every face (every part) is split into the harmonics of the rotation
  below that of Nyquist, and what it emits now is their sum at the phase of the moment. The
  momentum of the light absorbed and reflected is that of the moment. θ = 0 is the plain mesh.
- Surface term: plates, or a mesh with `surface` and `shadows`; occulters and `atmospheres` as
  in `myb_radiation_term`; the mass from the object (its mass function at t if given).
- Recoil term: power (W), efficiency η, `direction` (unit, body frame: where the radiation
  goes), `decay` > 0 makes the power fall as exp(−(t − epoch) / decay); a = −η P / (m c) d
  turned by the attitude law. `myb_force_emission` adds the same in a given frame.
- System light: the source body, the mesh of the lit body about its centre of mass, the flux at
  a distance, c, the mass, `shadows`; θ > 0 with `spin_axis` and `phases` (0: 64) for the late
  re-emission. The action gives the acceleration and the torque of the moment, the momentum of
  the light absorbed included.

## Domain of validity

- Mesh: `MYB_E_DOMAIN` for u not unit to 10⁻¹², fractions outside [0, 1], values not finite.
- Lagged mesh: also an axis not unit to 10⁻¹², θ negative or not finite, `phases` outside
  8…256.
- Sun–nadir law: `MYB_E_DOMAIN` when the body, the object and the source are in line.
- Recoil term: `MYB_E_DOMAIN` for a direction that is not unit.
- System light: `MYB_E_DOMAIN` for a lit body without a frame, a source that is the body itself,
  invalid values; the errors of the mesh.
- Plates: as in [radiation pressure](radiation-pressure.en.md).

## Accuracy

Tests [test_surface](../../tests/unit/test_surface.c), [test_shadow](../../tests/unit/test_shadow.c),
[test_nongrav](../../tests/unit/test_nongrav.c), [test_system](../../tests/unit/test_system.c):

| What | Reference | Result |
|---|---|---|
| A cube as a mesh (`test_surface`) | six plates of `myb_force_radiation_plates` | 10⁻¹⁶; no torque |
| A cube before a slab of 3×3 m, 24² parts a face (`test_surface`) | the silhouette, 9 m² | 9.000000000000 m² |
| A sphere of 1280 facets, Lambertian (`test_surface`) | (1 + 4/9) π R² p | 0.995 |
| Attitude laws (`test_surface`) | orthonormal, right-handed, axes as defined | orthonormal to 10⁻¹⁵ |
| A black plate facing the Sun in the term (`test_surface`) | `myb_term_radiation`, cr = 1 | 4·10⁻¹⁶ |
| A black cube as the mesh of the term (`test_surface`) | the sphere of its cross-section | 10⁻¹⁴ (bound) |
| Plates behind a planet with an atmosphere (`test_shadow`) | `myb_shadow_fraction_refracted` | 10⁻¹⁴ |
| Jacobian, light on plates, the four laws and one of the caller's (`test_surface`) | Richardson differences of the term | 1·10⁻¹⁰…7·10⁻¹⁰ |
| Jacobian, light on a mesh, with and without rays | the same | 1·10⁻¹⁰…5·10⁻¹⁰ |
| Jacobian, recoil | the same | 6·10⁻¹²…3·10⁻¹¹ |
| 63 W on 241 kg along a beam (`test_nongrav`) | the Pioneer anomaly, 8.74·10⁻¹⁰ m/s² | 8.72·10⁻¹⁰ m/s² |
| Recoil term on a spinning body, a quarter turn, decay of one e-fold (`test_nongrav`) | the direction turned, the power × e⁻¹ | to rounding |
| Light on what is left after a burn (`test_nongrav`) | ratio of accelerations = ratio of masses | 10⁻¹⁴ |
| Lagged mesh at θ = 0 (`test_surface`) | `myb_force_radiation_mesh` | 1.5·10⁻¹⁵ of the pressure |
| Lagged mesh, mean torque over a turn, θ = 0, 0.7, 3 (`test_surface`) | `myb_force_yorp` | 3.5·10⁻¹⁶ |
| Light action in a system, torque through the equations of rotation (`test_system`) | the mesh itself | 8·10⁻¹⁵ |
| Re-emitted torque about the axis over a turn (`test_system`) | `myb_force_yorp` | 10⁻¹⁴ (bound) |

Why the torque of the moment matters, measured: over a turn with the Sun at one latitude the
re-emission gives −6.572·10⁻⁴ p R³ about the axis and the light absorbed −1.026·10⁻³ p R³; the
second vanishes only after the average over the orbit.

## Cost

A ray through the tree of bounding boxes: 1.0 µs against 640 µs over every face (16 128 faces,
3000 rays, [test_geometry](../../tests/unit/test_geometry.c)). Self-shadowing casts one ray per
part (k² a lit face). The Jacobian with an attitude law of the caller's takes 12 calls of the
law (central differences over 6·10⁻⁶ of |x| and |v|), not of the term.

## Limitations

- Plates have no places: none hides another; only a mesh casts shadows on itself.
- A part of a face is lit or dark as its centre is: a coarse division is off at the edge of a
  shadow (with 5² parts a face the test finds the lit area off the 9 m² silhouette, by less than
  0.5 m²).
- Without the lagged mesh the re-emission is instantaneous; the lagged mesh takes the rotation
  as steady and the source as still over a turn.
- The terms move the centre of mass only; torques come from the mesh functions and from the
  light action of a system. The recoil is fixed in the body; its power is constant or decays
  exponentially.
- Jacobians: where the force is not smooth (a grazed plate, the parts of a mesh lit now) the
  derivative is that of the side the object is on; in a penumbra it is central differences of
  the term; a caller's attitude law is differentiated by central differences of the law.

## Data

None: plates, meshes, optical properties and powers are the caller's.
