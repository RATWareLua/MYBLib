# Contact and landing

**English** | [Русский](contact-landing.ru.md) · [All passports](index.en.md)

## What it is

The contact of a body with the surface of a shape: the landing, the rebound, sliding and coming
to rest. The body touches the surface with contact spheres fixed in its frame (one about the
centre for a ball, the pads of the legs of a lander); a sphere that presses in is pushed back by
the normal force of Hunt and Crossley, whose restitution has a closed form, and rubbed by a
Coulomb friction made smooth, or by a friction that sticks — a spring with a slip limit whose
deflection is a state of its own. The surface is a plane, a triaxial ellipsoid or a closed
triangular mesh. Three ways in: the force on a body in a given state, a term of a force model
with an event for the touch, and a level of a system of bodies with the torques on integrated
rotations. It is integrated with DOP853, or with IAS15 in `myb_orbit` given the timescale of the
jerk and the breaks at the touches.

## Source

- Normal force: K. H. Hunt, F. R. E. Crossley, *J. Appl. Mech.* 42, 440 (1975),
  f_n = k dⁿ (1 + χ ḋ) — zero where the surfaces touch, never pulling.
- Restitution: D. W. Marhefka, D. E. Orin, *IEEE Trans. Syst. Man Cybern. A* 29, 566 (1999),
  χ v_in − ln(1 + χ v_in) = −χ v_out − ln(1 − χ v_out), whatever k and n.
- Friction: f_t = −μ f_n u / √(|u|² + v_s²). A switch between sticking and sliding needs a state
  of its own (P. R. Dahl, 1968; C. Canudas de Wit et al., *IEEE Trans. Autom. Control* 40, 419,
  1995), which a term of a force model, a function of the state alone, cannot carry.
- Friction that sticks: the elasto-plastic element of Prandtl, the "spring with a slip limit" of
  the discrete element method (P. A. Cundall, O. D. L. Strack, *Géotechnique* 29, 47, 1979):
  f_t = −P_L(k_t z + c_t u), L = μ f_n, P_L(F) = F min(1, L/|F|), and the deflection follows the
  point, z′ = u − φ ẑ max(0, ẑ·u) − ẑ max(0, k_t|z| − 2L)/c_t, φ rising from 0 at 0.95 L to 1 at
  L (the band, the drawing back beyond twice the limit and the energy statement are the
  library's own, written in the header).
- IAS15 through a contact: the timescale τ₁² = 2|a′|² / (|a″|² + |a′||a‴|), one derivative up
  from that of D. Pham, H. Rein, D. S. Spiegel (2024) — see [IAS15](ias15.en.md).
- Nearest point of a mesh: `myb_polyhedron_nearest`; of an ellipsoid: `myb_ellipsoid_nearest`
  (both of [`myb/geometry.h`](../../include/myb/geometry.h)).
- Implemented from the papers; no third-party code, no data files.

## Interface

Header [`myb/contact.h`](../../include/myb/contact.h):

- the surface: `myb_contact_shape` (`MYB_CONTACT_PLANE`, `MYB_CONTACT_ELLIPSOID`,
  `MYB_CONTACT_MESH`), `myb_contact_surface`, `myb_contact_surface_nearest`;
- the law: `myb_contact_law` (with `kind`, `tangential_stiffness`, `tangential_damping`),
  `myb_contact_friction` (`MYB_CONTACT_FRICTION_SMOOTH`, `MYB_CONTACT_FRICTION_STICK`),
  `myb_contact_point`, `myb_contact_result`, `myb_contact_force`, `myb_contact_force_stick`,
  `myb_contact_restitution`, `myb_contact_damping`;
- the term: `myb_contact_term`, `myb_term_contact`, the event and break `myb_event_contact`.

Header [`myb/system.h`](../../include/myb/system.h): `myb_system_contact`,
`myb_system_add_contact`, `myb_system_pack_contact`, `myb_system_unpack_contact` — see
[encounters in a system](system-encounters.en.md). Contact problems are integrated with DOP853
(`myb_orbit` with `MYB_ORBIT_DOP853`, or `myb_dop853` of
[`myb/integ.h`](../../include/myb/integ.h); [classical integrators](classical-integrators.en.md)),
or with IAS15 in `myb_orbit`: `myb_orbit_options.timescale = MYB_IAS15_TIMESCALE_JERK`
(`myb_ias15_timescale` of `myb/integ.h`), `myb_event_contact` declared a break
(`myb_force_model_add_break`; either scheme of `myb_force_model_set_break_scheme`,
[`myb/model.h`](../../include/myb/model.h)).

## Inputs and outputs

- Every force is per unit mass of the body: k in m^(1−n) s⁻² (for the legs of a lander the
  stiffness of one leg in N m⁻ⁿ over the mass), n ≥ 1 (1 a linear spring, 1.5 the sphere of
  Hertz), χ ≥ 0 in s/m, μ ≥ 0, v_s > 0 in m/s where μ > 0.
- A sphere presses in by d = r − h, h the height of its centre over the surface (the distance to
  the nearest point, negative inside), and is pushed along the outward normal there; on a mesh
  the nearest point may be on a face, an edge or a corner, and the normal is the line from it to
  the centre of the sphere.
- `myb_contact_force`: x and v of the centre of mass in the frame of the surface (m, m/s), the
  attitude (rows: the body axes; NULL: the same axes), ω of the body relative to the surface in
  the body axes (NULL: none). The result: the force per unit mass in the axes of the surface
  (m s⁻²), the torque per unit mass about the centre of mass in the body axes (m² s⁻²), the
  energy of the springs Σ k d^(n+1) / (n + 1), the gap (the least height of a sphere less its
  radius, negative while one presses in), the number of points touching.
- Restitution e = v_out / v_in for a contact met at v_in along the normal with no other force
  during it: 1 at low speeds (1 − e = 2 χ v_in / 3 to the first order), 1/(χ v_in) at high ones;
  `myb_contact_damping` gives χ for e at a speed (0 for e = 1).
- Term: the object of `myb_orbit` as a sphere of `radius` about its centre, on the surface of a
  body of an environment given in that body's frame; the object is taken into that frame with the
  body's state, orientation and rate of rotation, so a turning body carries its surface along.
  With `env` = NULL the surface is fixed at the origin of the model. The body is not pushed back;
  no torque.
- `myb_event_contact`: the gap; its fall through zero is the touch, its rise the parting; an
  event function and a break of a force model with the term as context; NaN where the
  environment fails.
- System: the ground is pushed back by the ratio of the masses at the same points, so momentum
  and angular momentum are kept (a ground of gm 0 is not pushed); the torques turn the bodies
  whose rotation is integrated.
- Friction that sticks (`kind` = `MYB_CONTACT_FRICTION_STICK`): k_t > 0 in s⁻² and c_t ≥ 0 in
  1/s per unit mass where μ > 0 (v_s is not used). The deflection z of each point — the point
  less its anchor, m, in the axes of the surface — is carried with the state: in a system 3
  coordinates a point after every other kind (`myb_system_pack_contact`; zero — the springs at
  rest); with `myb_contact_force_stick` the caller's integrator gives `deflection` and gets
  `deflection_rate`. Below 0.95 of the limit the anchor holds, and a body on a slope below the
  angle of friction stays where it is, held by springs deflected by its weight (g sin α / k_t of
  a single point); at the limit the anchor is dragged and the force is μ f_n along the spring —
  Coulomb's law exactly; a spring left beyond twice the limit is drawn back at the rate
  k_t / c_t; a point off the surface keeps its deflection and is released against the limit when
  it touches again. The result's energy adds Σ k_t |z|²/2.
- Energy: dE/dt = −Σ (k χ dⁿ ḋ² + μ f_n |u|² / √(|u|² + v_s²)) ≤ 0; a contact never adds
  energy. With the friction that sticks, the energy of its springs counted, dE/dt ≤ 0 in every
  case (exact on a plane: on a curved surface the tangent plane turns under a point that
  slides).

## Domain of validity

- `MYB_E_DOMAIN`: an invalid surface (a plane normal of zero length, an ellipsoid radius that is
  not positive, no mesh), an invalid law (k ≤ 0, n < 1, χ < 0, μ > 0 with no v_s, values not
  finite) or point, no points, a state that is not finite; `MYB_E_NOCONV` from the ellipsoid.
- `myb_contact_restitution`: `MYB_E_DOMAIN` for negative or non-finite values (χ = 0 or v = 0
  give 1). `myb_contact_damping`: `MYB_E_DOMAIN` for e outside (0, 1], v not positive.
- The friction that sticks needs its state: `myb_contact_force` and `myb_term_contact` refuse it
  (`MYB_E_DOMAIN`); `myb_contact_force_stick` refuses a deflection that is missing or not
  finite. Keep c_t > 0: it damps the springs where the anchor holds and draws back a spring
  overstretched.
- IAS15 goes only in `myb_orbit` and only given three things. The timescale
  `MYB_IAS15_TIMESCALE_JERK`: with the timescale of the acceleration, which passes through zero
  twice in every bounce, a ball dropped on a plane stops IAS15 with `MYB_E_STEP_MIN` within its
  first contact. The breaks at the touch and the parting (`myb_event_contact`), so that no step
  holds the kink of a touch. Some other force on the body during the contact (gravity): where
  nothing else acts the acceleration leaves zero as a power of the time at both ends of a
  contact, and no timescale is defined there. A system of bodies with contacts goes with DOP853:
  it does not end its steps at the touches.

## Accuracy

All in [test_contact](../../tests/unit/test_contact.c), DOP853 at 10⁻¹³ unless said otherwise:

| What | Reference | Result |
|---|---|---|
| Restitution → damping → restitution, e = 0.02…0.999999 | the relation of Marhefka and Orin | 10⁻¹⁴ (bound); 1 − e → 2χv/3 and e → 1/(χv) at the limits |
| A ball (Hertz, k = 10⁹, χ = 0.5 s/m) thrown at a plane without gravity | the closed form | 3.5·10⁻¹³ |
| The same dropped under 1.62 and 0.162 m/s² at the same speed | the closed form | −2.2·10⁻⁴ and −2.2·10⁻⁵ (ratio 10.003): the impulse of the gravity during the contact |
| A ball dropped 0.2 m on the Moon | rest at the depth (g/k)^(1/n) | 24 bounces; at rest to 2·10⁻¹⁶ m; the energy never up by more than 10⁻¹⁶ m²/s² from one millisecond to the next |
| A block sliding at 2 m/s, μ = 0.4 | the closed form of the smooth law; Coulomb's block | 2·10⁻¹⁵ s; half the speed 6·10⁻⁸ s after Coulomb's; 0.509685 m slid against 0.509684 m |
| A cube as a mesh: rebound from a face | from a plane | 5·10⁻¹³ |
| From an edge, from a corner | back along the line to the ball; the closed form | 1.2·10⁻¹⁵; 10⁻¹² |
| A ball on the equator of a turning asteroid (a term with an environment) | the speed of the ground | 2·10⁻¹¹ m/s |
| A lander on four legs, a slope of 15°, μ = 0.6 (a system, its rotation integrated) | the creep speed v_s t / √(1 − t²), t = tan α / μ | 6·10⁻¹³; the legs carry its weight to 10⁻⁹; the energy never grows |
| The same on 40° | g (sin α − μ cos α) | 9·10⁻⁸ |
| A spinning ball off a turning ellipsoid of its system | momentum; angular momentum of orbits and spins | 2·10⁻¹⁶; 1.2·10⁻¹³ |
| Nearest point of a cube; of a star-shaped mesh of 528 faces | its closed form; a search of every face | exact; 4·10⁻¹⁶, inside and outside as the solid angle says at 3000 points |

Friction that sticks (k_t = 10⁴ s⁻², c_t = 100 s⁻¹, μ = 0.6, g = 1.62 m/s²; a system under
DOP853):

| What | Reference | Result |
|---|---|---|
| A block let go on 25° | rest at the deflection of its weight, 6.8464·10⁻⁵ m | to 10⁻⁹ of it (its overshoot, 0.90 of the limit, stays below the band); 7.5·10⁻¹⁸ m moved in the last 30 s of a minute; the same place to 1.6·10⁻¹⁸ m at rtol 10⁻⁹ … 10⁻¹² |
| The slope turned slowly to a part in a thousand below the angle of friction; a part above it | it holds; g (sin α − μ cos α) | 1.11 slip deflections down, 5.7·10⁻¹⁵ m in its last 40 s; slides to 2.5·10⁻¹³ |
| The lander on 15° (k_t = 1000 s⁻², c_t = 16 s⁻¹ a leg) | at rest | 9.3·10⁻¹⁴ m/s after 20 s (the smooth law creeps at 5.0·10⁻⁴ m/s), 1.6·10⁻¹⁴ m in its last 10 s; the legs bear its weight to 10⁻¹¹; the same place to 1.1·10⁻⁹ m at rtol 10⁻⁹ … 10⁻¹¹ |
| The same on 40° | g (sin α − μ cos α) | 8.3·10⁻⁸; the energy never grows by more than 4.4·10⁻¹⁶ m²/s² |

Under IAS15 (ε = 10⁻⁹, 10⁻¹⁰, 10⁻¹², the timescale of a′, the touches predicted) against DOP853
at rtol 10⁻¹²:

| What | Reference | Result |
|---|---|---|
| The ball dropped 0.2 m on the Moon | the 24 touches of DOP853 | 24 touches, at 7.0·10⁻⁸, 1.4·10⁻⁸, 7.6·10⁻⁹ s of them (the touches withdrawn: 5.0·10⁻⁸, 2.6·10⁻⁸, 3.3·10⁻⁹ s); at rest at its height to 1.5·10⁻¹⁶ m; the energy never up by more than 1.1·10⁻¹⁶ m²/s² |
| The sliding block | its three speeds by DOP853 (and the closed form) | at 7.9·10⁻¹¹, 1.0·10⁻¹¹, 9.1·10⁻¹⁴ s |
| The same ball with the timescale of the acceleration alone | — | `MYB_E_STEP_MIN` at the first contact |

## Cost

Steps counted in `test_contact` (no timings): the lander on 15° takes 3402 steps of DOP853 with
the friction that sticks against 7653 with the smooth law, on 40° 80 854 (the direction of a
spring of 0.2 mm follows a velocity of metres a second); the dropped ball under IAS15 13 789,
18 778 and 34 922 steps at ε = 10⁻⁹, 10⁻¹⁰, 10⁻¹² (the touches withdrawn: 14 039, 19 024,
35 180; DOP853 6729); the sliding block about 2300 (DOP853 489).

## Limitations

- The smooth friction does not stick: on a slope below the angle of friction a body creeps at
  v_s t / √(1 − t²). For |u| ≫ v_s it falls short of μ f_n by (v_s/u)²/2. Near rest it is
  stiff, μ f_n / v_s per unit of speed, and that sets the step.
- The friction that sticks has a state: no term of a force model can carry it, only a system
  or the caller's integrator. Its rate has corners (the band, the projection at the limit):
  continuous, not smooth, so DOP853 goes through measuring its error; a slide on 40° costs it
  tens of thousands of steps.
- No adhesion: a force that would be negative is held at zero.
- The closed-form restitution leaves out the other forces during the contact (the impulse of
  gravity, of the order g τ / v).
- A point meets the surface once, at its nearest point: a sphere caught in a concave corner of a
  mesh feels the nearer wall alone. The surface is at rest in its own frame.
- The term moves the object's centre only and does not push the body back; a body that turns
  needs the system level.
- IAS15 only in `myb_orbit`, with the timescale of the jerk, the breaks at the touches and some
  other force during the contact; a system with contacts goes with DOP853 (the system does not
  end its steps at the touches). Contact under RKN12(10) is not tested.

## Data

None: the surface, the law and the contact points are the caller's.
