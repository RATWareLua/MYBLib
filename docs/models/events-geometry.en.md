# Events and surface geometry

**English** | [Русский](events-geometry.ru.md) · [All passports](index.en.md)

## What it is

Three tools around a propagated orbit. Ready-made event functions g(t, y) whose sign changes
mark what is looked for — apsides and close approaches, nodes, a distance, an altitude, rise and
set over a station, the edges of a shadow, the impact on a shape. Breaks of the step: where a
force is continuous but not smooth (the edge of a shadow), a step that straddles the instant
loses its order, so the force model declares the instant and the integrator's step is withdrawn
and retaken exactly up to the root. And the geometry of a body's surface in its body-fixed
frame: a ray against a triaxial ellipsoid or a triangular mesh, the nearest point and the
altitude, normals, the limb and the terminator, the solid angle of a mesh, the view corrected
for light time, and the angles of illumination.

## Source

- Ray and ellipsoid in coordinates divided by the semi-axes, through the point of the line
  nearest the centre; the nearest point as the largest root of Eberly's equation by a monotone
  Newton iteration from a lower bound.
- Ray and mesh: Möller–Trumbore with faces widened by a part in 10¹², in a tree of bounding
  boxes; the inside/outside of the nearest point by the angle-weighted pseudo-normal of
  J. A. Baerentzen, H. Aanaes (*IEEE Trans. Vis. Comput. Graph.* 11, 243, 2005).
- The edges of a shadow: `myb_shadow_margin` of the radiation model (see
  [radiation pressure](radiation-pressure.en.md)); the view with light time uses the exact
  aberration of `myb/observe.h` ([observables and light time](observables-light-time.en.md)).
- Written from these descriptions; CSPICE N0067 (`gfoclt_c`, `surfpt_c`, `nearpt_c`,
  `sincpt_c`, `subpnt_c`, `ilumin_c`, `edlimb_c`, `edterm_c`, `spkpos_c`) is only the test
  oracle.

## Interface

- [`myb/events.h`](../../include/myb/events.h): `MYB_EVENT_ORIGIN`; `myb_range_event`,
  `myb_event_range`, `myb_event_range_rate`; `myb_plane_event`, `myb_event_plane`;
  `myb_altitude_event`, `myb_event_altitude`; `myb_elevation_event`, `myb_event_elevation`;
  `myb_eclipse_event`, `myb_event_eclipse`, `myb_term_radiation_breaks`; `myb_shape_event`,
  `myb_event_shape`. They are used with `myb_orbit_event` of
  [`myb/orbit.h`](../../include/myb/orbit.h).
- [`myb/model.h`](../../include/myb/model.h): `myb_break_function`, `MYB_MODEL_MAX_BREAKS`,
  `myb_force_model_add_break`, `myb_force_model_break_count`, `myb_force_model_break`; the
  scheme of the breaks `myb_break_scheme` (`MYB_BREAKS_RETAKE`, `MYB_BREAKS_PREDICT`),
  `myb_force_model_set_break_scheme`, `myb_force_model_break_scheme`.
- [`myb/geometry.h`](../../include/myb/geometry.h): `myb_ray_ellipsoid`,
  `myb_ellipsoid_nearest`, `myb_ellipsoid_normal`, `myb_ellipsoid_limb`,
  `myb_ellipsoid_terminator`; `myb_polyhedron_face_count`, `myb_polyhedron_face_normal`,
  `myb_polyhedron_ray`, `myb_polyhedron_solid_angle`, `myb_polyhedron_nearest`;
  `myb_view_body`, `myb_view_intercept`; `myb_illumination_angles`.

## Inputs and outputs

- Events: g(t, y) of the model time t and the state y = (x, v) of the object; bodies come from
  an environment ([`myb/env.h`](../../include/myb/env.h)) by index, `MYB_EVENT_ORIGIN` is its
  origin. `myb_event_range` is |x − x_body| − range (m); `myb_event_range_rate` the rate of the
  distance (m/s; zero and growing at the closest approach); `myb_event_plane` normal ·
  (x − x_body) (ascending node growing); `myb_event_altitude` the geodetic height over the
  body's ellipsoid minus a threshold (the body needs an orientation); `myb_event_elevation`
  the elevation over a `myb_station` minus a threshold (rad); `myb_event_eclipse` the
  `myb_shadow_margin` of a source and an occulter seen from the object (rad): `inner == 0` —
  the penumbra, otherwise the umbra or the annular phase; `myb_event_shape` 2π minus the solid
  angle of a mesh (positive outside).
- Breaks: functions of the same form declared to the force model, so `myb_orbit`, the
  estimator and the filter all see them. After a step the breaks are examined on its dense
  output; a step that contains a sign change is withdrawn and retaken from its beginning up to
  the root, with the polynomial of the withdrawn step as the prediction; the integrator is not
  restarted (the state, Φ and S continue). That is `MYB_BREAKS_RETAKE`, the default. With
  `MYB_BREAKS_PREDICT` (IAS15; from the next `myb_orbit_start`) the break functions are examined
  before a step on the polynomial of the last step carried over the step to come, a predicted
  sign change becomes the end of the step, and once taken the step is reshaped (warm, in a sweep
  or two) until the root of its own polynomial is at its end to 2⁻³⁰ of the step — no step is
  taken across a break and thrown away; the step after a break is planned at most a quarter of
  the one before it, a break the prediction misses is met as by withdrawing. DOP853 and RKN12
  restart at a break in either case. `myb_term_radiation_breaks` declares two breaks per
  occulter of a radiation term.
- Geometry: everything in the body-fixed frame, centre at the origin, in m; the caller turns
  inertial vectors with the body's orientation. Ellipsoids by their semi-axes along x, y, z;
  meshes are the `myb_polyhedron` of `myb/grav.h`. `myb_ellipsoid_nearest` gives the point and
  the height along the normal (negative inside); `myb_polyhedron_nearest` the point, the signed
  distance, the outward normal and the face. The limb is an ellipse (centre, major, minor
  semi-axes); the terminator `count` points in half-planes about the line to the source, the
  edge of the full shadow (`penumbral == 0`) or of the full light.
- View: the body taken where it was, and turned as it was, a light time τ earlier (iterated on
  the environment's states); with `aberration != 0` the direction turned by the observer's
  motion with the exact transformation. Illumination: phase, incidence, emission in [0, π].

## Domain of validity

- An event function has no status: when a provider of the environment fails it gives NaN —
  no sign change, and the event is not found.
- At most `MYB_MODEL_MAX_BREAKS` = 32 breaks per model (`MYB_E_OVERFLOW`, also from
  `myb_term_radiation_breaks`).
- `myb_ray_ellipsoid`, `myb_polyhedron_ray` return 0 for a miss, a zero direction or
  non-finite input (radii must be positive); from inside they give the point where the ray
  leaves.
- `myb_ellipsoid_nearest`: `MYB_E_DOMAIN` (non-finite input), `MYB_E_NOCONV`; deep inside, in
  a plane of symmetry, of two equally near points (or a circle) the one on the positive side of
  the shortest axis is given. `myb_ellipsoid_limb`: `MYB_E_DOMAIN` for an observer on or inside
  the ellipsoid; `myb_ellipsoid_terminator`: for a source not outside or a negative radius.
- `myb_polyhedron_solid_angle` is 4π inside and 0 outside; on the surface itself anything
  between. `myb_polyhedron_nearest`, `myb_illumination_angles`, `myb_view_*`: `MYB_E_DOMAIN`
  for non-finite input or degenerate vectors, plus the statuses of the environment.

## Accuracy

| What | Reference | Result |
|---|---|---|
| Final position over ε = 10⁻⁹…10⁻¹², 12 000 km orbit, 0.2 m²/kg, two days, 54 edges | the same at ε = 10⁻¹³ | 6.4 m without breaks, 8·10⁻⁷ m with them |
| The same, the breaks by prediction (`MYB_BREAKS_PREDICT`) | the same; withdrawing | within 8.4·10⁻⁷ m (withdrawing 8.2·10⁻⁷ m), the two 4·10⁻⁷ m apart at 10⁻¹³; without the quarter step after an edge 1.2·10⁻⁴ m at 10⁻⁹ (the predicted edges right to 3·10⁻⁷ s) |
| LAGEOS week, both schemes, ε = 10⁻⁹ … 10⁻¹² | each other | the first-guess residual within 3·10⁻⁶ m of withdrawing at 10⁻¹², the fit 1.5685 cm either way |
| Φ across eight edges; ∂x/∂c_r | differences of whole propagations | 7·10⁻⁸ (largest element 25); 2·10⁻⁷ m of 49 m |
| LAGEOS week, residual of the first guess at ε = 10⁻¹⁰, 10⁻¹¹, 10⁻¹² | itself | 25.9587 m at all three in phase 16 (without breaks 25.9520 and 25.9582 m); 25.8887 m since the secondary ocean-tide waves of phase 20 |
| Apsides; nodes and a sphere, ellipse e = 0.3 | Kepler's equation | 2·10⁻¹¹ s; 10⁻¹² s |
| Altitude event; elevation event | the threshold; `myb_measure_azel` | 10⁻⁹ m; 2·10⁻¹⁴ rad |
| Impact on an octahedron | its face | 2·10⁻⁹ m |
| 28 penumbra and umbra contacts in a day | CSPICE `gfoclt_c` (convergence 10⁻⁹ s) | within 4.4·10⁻¹⁰ s |
| Ray and ellipsoid, five shapes, 300 rays | `surfpt_c` | 5·10⁻¹⁶ of the distance |
| Nearest point, 300 points including interior ones | `nearpt_c` | 4·10⁻¹⁶ of the size, altitude 3·10⁻¹⁶ |
| Ray and mesh, convex and dented icosahedron, 8000 rays | planes of the faces and areas | no disagreement, distances 2·10⁻¹⁵ |
| Tree of boxes, 16 128 faces, 3000 rays | a pass over every face | no difference in distance or face |
| Nearest point of a mesh: a cube, 4000 points; a star-shaped mesh, 3000 points | closed form; a search over all faces | < 10⁻¹⁵ and < 10⁻¹⁴ (tolerances), no point on the wrong side |
| Surface through kernels, 20 epochs | `sincpt_c`, `subpnt_c`, `ilumin_c` | points 10⁻⁸ m (the SPK interpolation), altitude 3·10⁻⁹ m, angles 7·10⁻¹⁶ rad |
| Limb; terminator | `edlimb_c`; the horizon condition | 2·10⁻⁷ m; 5·10⁻⁵ m (the points of `edterm_c`: 1.4·10⁻⁴ m) |
| A ray with light time; the Sun from the object | `sincpt_c`, `spkpos_c` with CN | 7·10⁻⁹ m and 5·10⁻¹² s; 9·10⁻⁵ m |

Tests: [`test_events`](../../tests/unit/test_events.c) (breaks with IAS15, DOP853 and RKN12,
events, CSPICE), [`test_geometry`](../../tests/unit/test_geometry.c),
[`test_contact`](../../tests/unit/test_contact.c) (`myb_polyhedron_nearest`),
[`test_lageos`](../../tests/unit/test_lageos.c) (176 edges of the shadows of the Earth and the
Moon in the week).

## Cost

- Breaks: withdrawing the step costs +26 % of the time of the LAGEOS fit (the first version,
  a cold restart at every edge, cost +70 %); LAGEOS takes 2290 steps a week instead of 2177.
- Prediction against withdrawing: on the orbit of two days 557, 742, 993, 1362 steps against
  517, 701, 953, 1308 at ε = 10⁻⁹ … 10⁻¹² (a step more an edge, the quarter step after it); on
  the LAGEOS fit 160 570, 188 045, 240 809, 314 696 evaluations of the model against 169 778,
  196 804, 244 341, 318 308 at ε = 10⁻⁹, 10⁻¹⁰, 10⁻¹¹, 10⁻¹² (5.4, 4.5, 1.4, 1.1 % fewer), and no
  gain in time (the least of three rounds on a shared machine, withdrawing / predicting: 28.6 /
  26.0, 31.9 / 35.7, 40.5 / 44.3, 53.3 / 50.9 s; one case spread over a factor of two). Most of
  the 26 % is the shorter steps around the edges, which every scheme pays.
- A ray on a mesh of 16 128 faces through the tree of boxes: 1.0 µs, against 640 µs over every
  face.

## Limitations

- The event geometry is instantaneous: no light time; the elevation has no refraction.
- `myb_event_shape` jumps at the surface: the crossing is bracketed and halved, some fifty
  evaluations, each a pass over the faces.
- By default each edge of a shadow costs one withdrawn step. Ending the step short of the edge
  predicted from the polynomial (phase 16) lost: 7·10⁻⁴ m at ε = 10⁻⁹ and 5·10⁻⁶ m at 10⁻¹²
  against 8·10⁻⁷ m for withdrawing. The scheme that reshapes the step to the edge its own
  polynomial finds (`MYB_BREAKS_PREDICT`) is as accurate, but saves only 1–5 % of the
  evaluations and no time that shows, so it is not the default.
- With aberration `myb_view_body` differs from CSPICE's `CN+S` by 13 m at 1 au: CSPICE takes
  the classical formula, the library the exact Lorentz transformation (they differ by (v/c)²,
  25 m). Without the light-time correction a point on the Earth is 8 m off.
- A ray through an edge shared by two faces meets the first of them in the mesh; of faces
  equally near a point, `myb_polyhedron_nearest` gives the first.
- In phase 22 `myb_ellipsoid_nearest` stopped refusing points just off a sphere (its Newton
  steps fell below the resolution); every value it returned before is the same.

## Data

None: the shapes, the environment and the stations are the caller's. The CSPICE references
`tests/data/events/cspice.ref` and `tests/data/geometry/ellipsoid_cspice.ref` were made with
CSPICE N0067 as an external oracle.
