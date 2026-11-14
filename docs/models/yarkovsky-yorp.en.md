# Yarkovsky and YORP

**English** | [Русский](yarkovsky-yorp.ru.md) · [All passports](index.en.md)

## What it is

The recoil of the heat that a small body re-emits. The Yarkovsky effect pushes its centre of
mass: the library gives the drift of the semi-major axis from the thermal properties of the body
(the linear theory of heat conduction in a sphere), the answer of a conducting sphere of any
size, the diurnal acceleration at any instant on any orbit with its gradient, and the seasonal
acceleration on an elliptic orbit of eccentricity up to 0.95 as a Fourier series in the mean
anomaly. The YORP effect turns the body: the torque of a triangular mesh that re-emits as a
Lambertian surface, averaged over a rotation, without conduction (Rubincam's model) and with
it.

## Source

- D. Vokrouhlický, *A&A* 335, 1093 (1998) and 344, 362 (1999): the linear theory of the diurnal
  and seasonal drifts of a body much larger than the depth of the thermal wave (κ₁ = κ₂ = κ₃ = ½).
- The answer of a sphere to a periodic heating, 1 / (1 + λ ζ(z)), ζ = z j₁′(z) / j₁(z): the
  solution of the linear heat conduction (Vokrouhlický 1998). For a radius below about 0.28 of
  the depth of the thermal wave a series to z¹² replaces the closed form, which loses eight
  digits there; the two meet at 10⁻¹².
- Seasonal effect on eccentric orbits: the linear theory of D. Vokrouhlický, P. Farinella,
  *AJ* 118, 3049 (1999), each harmonic k of the heating along the spin axis answered at the
  frequency k n.
- YORP: D. P. Rubincam, *Icarus* 148, 2 (2000) — the zero-conductivity model; with conduction
  the axial torque does not depend on the thermal parameter (D. Čapek, D. Vokrouhlický,
  *Icarus* 172, 526, 2004).
- Implemented from the papers; no third-party code, no data files.

## Interface

Header [`myb/force.h`](../../include/myb/force.h):

- structure `myb_yarkovsky`; functions `myb_yarkovsky_drift`, `myb_yarkovsky_response`,
  `myb_force_yarkovsky`, `myb_force_yarkovsky_gradient`;
- the seasonal series `myb_yarkovsky_season` (opaque): `myb_yarkovsky_season_create`,
  `myb_yarkovsky_season_destroy`, `myb_force_yarkovsky_season`,
  `myb_force_yarkovsky_season_gradient`, `myb_yarkovsky_season_drift`,
  `myb_yarkovsky_season_harmonics`.

Header [`myb/surface.h`](../../include/myb/surface.h): `myb_force_yorp`; the instantaneous
force and torque of the same mesh, `myb_force_radiation_mesh` and
`myb_force_radiation_mesh_lagged`, are in [radiation on surfaces](radiation-surfaces.en.md).

Header [`myb/env.h`](../../include/myb/env.h): `myb_term_yarkovsky`,
`myb_term_yarkovsky_jacobian` (context `myb_yarkovsky_term`). The empirical A2 law that orbit
fits use for asteroids is in [comets and outgassing](comets-outgassing.en.md).

## Inputs and outputs

- `myb_yarkovsky`: radius (m), density (kg m⁻³), absorptivity (1 − Bond albedo), emissivity,
  thermal inertia √(K ρ C) (J m⁻² s⁻½ K⁻¹), rotation rate (rad/s), obliquity (rad, from the
  orbital angular momentum to the spin), conductivity (W m⁻¹ K⁻¹; 0: the body is taken much
  larger than the thermal wave, otherwise its depth K / (Γ √ν) is compared with the radius).
- `myb_yarkovsky_drift`: the solar flux on the orbit (W m⁻²), the mean motion (rad/s), c; gives
  da/dt (m/s) of the diurnal effect, outwards for a prograde rotation, and of the seasonal one,
  always inwards (either output may be NULL).
- `myb_yarkovsky_response(ratio, θ)`: the complex amplitude relative to a body without thermal
  inertia; the imaginary part is negative (the answer lags); ratio ≤ 0 gives the large-body
  limit 1 / (1 + θ (1 + i) / 2).
- `myb_force_yarkovsky`: a += −(4α/9) Φ [Re(R) s − Im(R) (spin × s)], s the part of the unit
  direction to the Sun across the spin axis, Φ the light's force per unit mass at the flux, R the
  response at the rotation frequency; `obliquity` is not used. The gradient takes the vector d
  from the body to the Sun (m), the flux falling as 1/|d|², and writes dadr[3 i + j].
- Season: the state of the body relative to the Sun (m, m/s) fixes the orbit, gm of the Sun, the
  spin axis (unit), the flux at a distance. The acceleration is along the spin axis; the place
  on the orbit is read from the direction of r in the plane of that orbit. The series is cut
  where its terms fall below 10⁻¹⁴ of the first.
- `myb_force_yorp`: the mesh turns about its z axis; the Sun at a colatitude from that axis; θ
  the thermal parameter at the rotation frequency (0: Rubincam's model); `phases` ≥ 8 equal steps
  of a turn. torque (N m): [0] across the axis in the plane of the axis and the Sun, on the Sun's
  side; [1] completing; [2] along the axis.
- Term: the spin axis fixed in the inertial axes; the diurnal part and, with `season`, the
  seasonal one; position and velocity relative to the source body.

## Domain of validity

- With `conductivity` = 0 the body must be much larger than the thermal wave: metres for the
  day, tens of metres for the year — asteroids of a hundred metres and more. With a conductivity
  any size is allowed.
- Season: elliptic orbits with e ≤ 0.95. `myb_yarkovsky_season_create` refuses with
  `MYB_E_DOMAIN` an invalid body, a spin that is not a unit vector, an orbit that is not an
  ellipse of e ≤ 0.95 (tested at e = 0.954 and on a hyperbola); `MYB_E_NOMEM`.
- `MYB_E_DOMAIN` from the drift, the force and the gradient for invalid values (a negative
  radius, a spin axis that is not unit, d = 0).
- `myb_force_yorp`: `MYB_E_DOMAIN` for θ negative or `phases` below 8.
- The theory is linear: the temperature is linearised about its mean.

## Accuracy

Tests [test_nongrav](../../tests/unit/test_nongrav.c) and
[test_surface](../../tests/unit/test_surface.c):

| What | Reference | Result |
|---|---|---|
| A kilometre of rock at 2.5 au, a rotation of 6 h, obliquity 0 | of the order of 10⁻⁴ au/Myr | 2.2·10⁻⁴ au/Myr outwards; seasonal 0 |
| The same | the formula of the linear theory written out | 10⁻¹³ |
| Retrograde; on its side | minus the prograde drift; the seasonal alone, inwards | 10⁻¹⁵ |
| Twice the radius; the thermal parameter √2 | half the drift; the maximum of the formula | 10⁻¹⁵; 10⁻¹³ |
| Response of a conducting sphere, radius 0.02…10 depths of the wave | heat equation solved in mpmath (30 digits) | 5·10⁻¹³ |
| 5 m of rock with K = 0.5; a grain of 2 mm | the large-body limit | within 2 % (bound); the diurnal drift falls below half of it |
| Diurnal term, 5 m of rock, circular orbit, 3 revolutions | `myb_yarkovsky_drift` | 10⁻⁴ (the test's bound) |
| Seasonal accelerations and mean drift, e = 0…0.85, 5 bodies | mpmath, 25 digits, other quadrature nodes, the response from the Bessel function | 8·10⁻¹⁴; 4·10⁻¹³ |
| Seasonal mean drift on a circle | the closed form of `myb_yarkovsky_drift` | 5·10⁻¹⁵ |
| `myb_term_yarkovsky_jacobian` | central differences | 1.3·10⁻¹⁰, with the season 5·10⁻¹⁰ |
| YORP without conduction, a symmetric body (sphere of facets) | zero | 6·10⁻¹⁷ |
| A body with "vanes" and its mirror image | opposite torques | −6.2190·10⁻³ and +6.2190·10⁻³ p R³, equal to 10⁻¹² |
| `myb_force_yorp` at θ = 0 | the mesh averaged over a turn, re-emission alone | 10⁻¹⁴ (bound) |
| Axial torque at θ = 0, 1, 30 | independent of θ (Čapek and Vokrouhlický) | −2.1924·10⁻³ p R³, the same bits |
| Transverse torque at θ = 1 | that of θ = 0 times 1/(1.5 + 0.5 i): 1/√2.5 of it, turned by atan(0.5/1.5) | 10⁻¹² |

## Cost

Not measured as time. The seasonal series keeps 2 harmonics on a circle, 38 at e = 0.3, 111 at
0.6, 551 at 0.85 (docs/PLAN.en.md); the header gives 70 at e = 0.5 and 3000 at 0.95.

## Limitations

- The Yarkovsky force is that of a sphere in the linear theory; the spin axis is fixed in space
  and its own evolution (by YORP) is not coupled to the orbit.
- `myb_force_yarkovsky` is the diurnal part only. The seasonal part belongs to one orbit: the
  series is made for the orbit given at its creation, and a body whose orbit drifts away needs a
  new one. `myb_yarkovsky_drift` is the mean on a circular orbit.
- YORP is the mean over a rotation of the re-emitted light alone; the momentum of the light
  absorbed is left out (over a rotation and an orbit it averages away). The torque of the moment,
  absorbed light included, is that of the mesh functions and of the light action of a system
  ([radiation on surfaces](radiation-surfaces.en.md)). The faces re-emit as Lambertian
  radiators.

## Data

None: the thermal properties, the mesh and the spin are the caller's.
