# Aerodynamics: drag, gas–surface interaction, entry

**English** | [Русский](aerodynamics.ru.md) · [All passports](index.en.md)

## What it is

The forces of an atmosphere on an object, from the thin gas of an orbit to the dense air of an
entry. The drag of a sphere in an atmosphere that turns with its body, with winds; a layered
exponential atmosphere; the free-molecular flow on plates (Sentman), with drag and lift; the
transition to the continuum by the Knudsen number; the modified Newtonian continuum on plates and
on meshes that hide their parts from the flow; the convective heating of the stagnation point
(Sutton–Graves) and the ablation of the mass it causes, the mass carried in the state.

## Source

- Drag: a = −½ ρ (c_d A/m) |w| w, w the velocity relative to the atmosphere. Orbital decay is
  checked against D. King-Hele, *Satellite Orbits in an Atmosphere* (1987), eq. 4.9.
- Free-molecular flow with diffuse re-emission: L. H. Sentman, Lockheed LMSC-448514 (1961), in
  the form of E. Doornbos, *Thermospheric Density and Wind Determination from Satellite Dynamics*
  (2012), eqs. 3.50–3.53. The agreement of the plate formula, integrated over a sphere, with the
  closed form of the sphere was checked to 10⁻³⁰ in mpmath before the code was written.
- Transition: R. G. Wilmoth, R. A. Mitcheltree, J. N. Moss, *J. Spacecraft Rockets* 36, 436
  (1999): C = C_cont + (C_free − C_cont) sin²(π (3 + log₁₀ Kn) / 8), between Kn = 10⁻³ and 10;
  mean free path M / (√2 π d² ρ N_A).
- Modified Newtonian flow: L. Lees (1955); J. D. Anderson, *Hypersonic and High-Temperature Gas
  Dynamics* (2006), ch. 3. C_p,max from Rayleigh's pitot formula (NACA Report 1135, 1953); below
  M = 1 the isentropic stagnation pressure as a continuation.
- Heating: K. Sutton, R. A. Graves, NASA TR R-376 (1971), q = k √(ρ / R_n) V³,
  k = 1.7415·10⁻⁴ kg^½ m⁻¹ for the air of the Earth. Ablation: dm/dt = −q A / Q*.
- Implemented from the papers; no third-party code, no data files.

## Interface

Header [`myb/force.h`](../../include/myb/force.h): `myb_force_drag`, the callbacks
`myb_density_function` and `myb_wind_function`, `myb_exponential_atmosphere`
(`MYB_ATMOSPHERE_LAYERS`), `myb_density_exponential`.

Header [`myb/env.h`](../../include/myb/env.h): `myb_term_drag`, `myb_term_drag_jacobian`
(context `myb_drag_term`).

Header [`myb/surface.h`](../../include/myb/surface.h):

- free molecules and transition: `myb_gas`, `myb_force_gas_plates`, `myb_gas_sphere_cd`,
  `myb_gas_knudsen`, `myb_gas_bridge`;
- continuum: `myb_air`, `myb_air_function`, `myb_air_heat_ratio`, `myb_newtonian_cp_max`,
  `myb_force_newtonian_plates`, `myb_force_newtonian_mesh`, `myb_continuum`;
- the term `myb_term_surface_drag`, `myb_term_surface_drag_jacobian` (`myb_surface_drag_term`);
- entry: `myb_heating_sutton_graves`, `myb_ablation`, `myb_ablation_rate`, `myb_ablation_rhs`,
  `myb_ablation_mass`.

Density models from the literature are in [thermosphere models](thermosphere-models.en.md) and
[JB2008](jb2008.en.md); the attitude laws of the surface term in
[radiation on surfaces](radiation-surfaces.en.md).

## Inputs and outputs

- SI: ρ in kg m⁻³, speeds in m/s, temperatures in K, molar masses in kg/mol, lengths in m,
  heating in W m⁻², mass rate in kg/s. `myb_force_*` add to `a[3]`; the mesh function writes
  force (N) and torque (N m about the origin of the mesh); the terms write.
- Density and wind functions take the body-fixed position and the instant; the wind is the
  velocity of the gas relative to the rotating body, in its body-fixed frame. The drag term
  reads the density at the instant on the scale of the body's orientation and uses the velocity
  in the body-fixed frame (v − ω × r for a uniform rotation) less the wind.
- Exponential atmosphere: geodetic heights `base[k]` over the ellipsoid, increasing;
  ρ = density[k] exp(−(h − base[k]) / scale[k]) up to `base[k + 1]`; below `base[0]` the first
  layer continues; up to 32 layers; t is not used.
- Gas: its temperature and molar mass, the wall temperature, the accommodation of energy
  (0…1). Plates: area and outward normal; every plate counts, the faces turned away included;
  a += ½ ρ V² / m Σ A_i (C_D,i d + C_L,i l_i), v relative to the gas in the frame of the plates.
  `myb_gas_sphere_cd` is referred to the cross-section.
- Air of the continuum: temperature, molar mass, γ; M = V / √(γ R T / M_mol). An element facing
  the flow takes ½ ρ V² C_p,max sin²θ along −n, sin θ = n·v/V; one turned away takes nothing.
  C_p,max is 1.8394 as M → ∞, 1.8088 at M = 5, 1.6573 at M = 2 for γ = 1.4.
- Surface drag term: the free molecules on its plates; with `length` > 0 the bridge
  a = W a_free + (1 − W) a_cont, W = `myb_gas_bridge(Kn)`. The continuum is either the drag
  −½ ρ |w| w c_d A / m of the single number `continuum_cd` (with `continuum_area`) or, with
  `continuum`, the Newtonian flow on the plates or on a mesh; the air of `myb_continuum.air`,
  when given, also serves the free molecules and the Knudsen number.
- Ablation: the mass is a seventh component of the state y = (x, v, m); `myb_ablation_rhs` is a
  right-hand side for an integrator of first-order systems, and `myb_ablation_mass` as the mass
  function of the object makes every term act on the mass of the state being evaluated. The
  radius of the nose and the heated area stay as given.

## Domain of validity

- Drag term: `MYB_E_DOMAIN` for a density that is negative or not finite, invalid values of the
  object; the statuses of the density and wind functions pass through.
- Exponential atmosphere: `MYB_E_DOMAIN` for an invalid table (no layers or more than 32, bases
  not increasing, a scale that is not positive), an invalid ellipsoid, a non-finite position.
- Free molecules: `MYB_E_DOMAIN` for an invalid gas (accommodation outside [0, 1]), mass,
  density or plate (normal not unit), v = 0. The Knudsen number is infinite at zero density; the
  bridge weight is 0 below Kn = 10⁻³ and 1 above 10.
- `myb_newtonian_cp_max`: `MYB_E_DOMAIN` for γ ≤ 1 or a Mach number negative or not finite.
  The theory is one of hypersonic speeds; below M = 1 the branch only gives a descent through the
  speed of sound a force, it does not make the theory right there. `myb_air_heat_ratio` returns 0
  for amounts that are not valid.
- Surface term with a continuum: `MYB_E_DOMAIN` for an air of zero temperature, γ ≤ 1, a
  molecule diameter of 0 with a length.
- Sutton–Graves: `MYB_E_DOMAIN` for negative or non-finite values, a nose radius that is not
  positive. Ablation: `MYB_E_DOMAIN` for invalid values (a heat of ablation of 0) and for a mass
  that is not positive (the object is consumed).

## Accuracy

Tests [test_env](../../tests/unit/test_env.c), [test_orbit](../../tests/unit/test_orbit.c),
[test_surface](../../tests/unit/test_surface.c), [test_entry](../../tests/unit/test_entry.c):

| What | Reference | Result |
|---|---|---|
| Drag term (`test_env`) | `myb_force_drag` by hand with the providers' values; −½ ρ c_d A/m \|w\| w with w = v − ω × r | bit for bit; 1.4·10⁻¹⁶ |
| Exponential atmosphere at six heights, below the base and in every layer (`test_env`) | its formula | bit for bit |
| `myb_term_drag_jacobian` | central differences | 4·10⁻⁸ (position), 2·10⁻¹¹ (velocity); with the gradient of a wind 2.4·10⁻⁹ |
| Decay of a circular orbit at 400 km, 10 revolutions (`test_orbit`) | first order, 2π c_d A/m ρ a² a revolution | −25.408 m against −25.403 m |
| A wind to the east in proportion to the distance from the axis (`test_surface`) | a faster rotation of the atmosphere | 10⁻¹⁴ |
| Sentman's plates over a sphere, speed ratios 0.3…31, accommodation 0, 0.8, 1 (`test_surface`) | the closed form of the sphere | drag 5·10⁻¹⁵, lift 2·10⁻¹⁵ |
| C_D of a sphere at 7.6 km/s in atomic oxygen at 900 K, accommodation 0.8, wall at 300 K | — | 2.41 |
| A plate across a cold fast flow | 2 + √π v_out / v_in | 10⁻⁴ (bound) |
| Air at 100 km, 5.6·10⁻⁷ kg m⁻³, 28.4 g/mol, d = 3.7·10⁻¹⁰ m | — | free path 0.14 m; a body of 2 m: Kn = 0.069, free-molecular weight 0.44 |
| Surface term at Kn = 10⁷; at a free path of a millimetre | the free molecules; the continuum drag given | bit for bit; 10⁻¹⁴ |
| `myb_term_surface_drag_jacobian`, five attitude laws, a wind, the transition | Richardson differences of the term | 2.4·10⁻⁹ |
| C_p,max, M = 1…100 (`test_entry`) | Rayleigh's pitot formula written out | 10⁻¹⁴ (bound) |
| p₀₂/p∞ at M = 2, 3, 5, 10; p/p₀ at M = 0.5 | NACA 1135, tables II and I | to their printed digits |
| C_p,max as M → ∞ (γ = 1.4 and 5/3); at M = 1 from both branches | closed forms | 10⁻¹⁴ |
| A flat plate at 5…90° | C_p,max sin²α (normal), sin³α (drag), sin²α cos α (lift) | 10⁻¹⁵ |
| Plates over a sphere, a quadrature exact for cos³ | C_D = C_p,max / 2 | 4.4·10⁻¹⁶ |
| A sphere of 1280 facets | C_p,max / 2 | within 1 % (bound) |
| A pyramid of 64 faces, 20° | C_p,max sin²δ_f of its faces | 2.2·10⁻¹⁶ |
| Two cubes in a line, the rays cast | the front face alone takes the flow | 10⁻¹³ (bound) |
| The Newtonian sphere in the bridge, density 10⁻¹²…10⁻² | `continuum_cd` = C_p,max / 2 | 10⁻¹⁴ |
| The air of a function; a mesh as the continuum | the parts computed apart; `myb_force_newtonian_mesh` | 10⁻¹⁴ |
| Jacobian of the term with a continuum, 90 km | Richardson differences | 5.2·10⁻⁶ |
| Sutton–Graves at 7.5 km/s, 3·10⁻⁴ kg m⁻³, a nose of 1 m | 1.27253·10⁶ W m⁻² by hand | to its six digits |
| 10 kg at 12 km/s in air of 10⁻⁴ kg m⁻³, 200 s, the mass in the state (DOP853) | m = m₀ exp(−(B/2C)(v₀² − v²)) | 1.9·10⁻¹⁴ |

## Cost

Not measured as time. With a continuum the Jacobian of the surface term is central differences
of the term over 6·10⁻⁶ of |x| and |v|: 12 calls of the term. A part of the bridge whose weight
is 0 is not computed. The flow shadows of a mesh cast one ray per part, as in
[radiation on surfaces](radiation-surfaces.en.md).

## Limitations

- The sphere's drag has one c_d A/m; the atmosphere turns rigidly with its body, plus the wind
  of a function.
- Sentman's flow: diffuse re-emission with an accommodation of the energy. Plates have no
  places, so none hides another from the flow; only a mesh in the Newtonian continuum casts flow
  shadows. A term whose continuum is a mesh and that has no plates feels nothing in a
  free-molecular flow.
- The Newtonian theory: a perfect gas of one γ, the base pressure taken as p∞, nothing on the
  faces turned away; hypersonic only.
- The bridge between the regimes is the empirical form of Wilmoth et al.
- Heating: the convective heating of the stagnation point only; the heat radiated by the shock
  layer, which takes over above some 10 km/s for a large nose, is not in it. The nose radius and
  the heated area do not change as the mass burns away; the mass in the state needs an
  integrator of first-order systems (a term cannot keep it).
- Jacobians: the gradients of the density and of the wind are central differences of their
  functions over 10⁻⁷ of the distance from the body's centre; with a continuum the Jacobian is
  differences of the term, and a step across a change of the Mach number or of the parts the
  flow reaches gives the slope of that step.

## Data

None: the layers of the exponential atmosphere, the gas and the air are the caller's.
