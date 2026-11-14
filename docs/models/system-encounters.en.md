# Encounters, collisions and mergers in a system

**English** | [Русский](system-encounters.ru.md) · [All passports](index.en.md)

## What it is

Events between the bodies of a [system of bodies](system-of-bodies.en.md). Bodies given a radius
are spheres that may touch; `myb_system_advance` steps the IAS15 integrator of the system
towards a final time and stops at the first event — the contact of two spheres, found on the
polynomial of the step and made the end of the step (IAS15 withdraws the step and retakes it to
that instant), or the closest approach of two bodies that pass within a watched distance. What
happens then is the caller's: merge the two (`myb_system_merge`, rotating bodies and PN pairs
with spins included), let them rebound (`myb_system_bounce`), or go on. Only the pairs whose
paths over the step can meet are examined: boxes from the step's own polynomial, swept along one
axis. A body resting on or sliding over the surface of another — a lander on its legs — is a
separate level, `myb_system_add_contact`.

## Source

- The event search, the merger and the rebound are implemented in the library from the
  conservation of mass, momentum and angular momentum; no third-party code.
- The contact level uses the law of [`myb/contact.h`](../../include/myb/contact.h): the normal
  force of K. H. Hunt and F. R. E. Crossley (*J. Appl. Mech.* 42, 440, 1975), whose restitution
  follows D. W. Marhefka, D. E. Orin (*IEEE Trans. Syst. Man Cybern. A* 29, 566, 1999), and a
  smooth Coulomb friction — see [contact and landing](contact-landing.en.md).

## Interface

Header [`myb/system.h`](../../include/myb/system.h):

- functions `myb_system_set_radius`, `myb_system_radius`, `myb_system_present`,
  `myb_system_advance`, `myb_system_merge`, `myb_system_bounce`, `myb_system_add_contact`,
  `myb_system_pack_contact`, `myb_system_unpack_contact`;
- structures `myb_system_event` (kinds `MYB_SYSTEM_ARRIVED`, `MYB_SYSTEM_CONTACT`,
  `MYB_SYSTEM_APPROACH` of `myb_system_event_kind`), `myb_system_merger`,
  `myb_system_contact`;
- the integrator from `myb_system_integrator`, started with `myb_system_acceleration`
  ([IAS15](ias15.en.md)); the surface, law and points of the contact level
  (`myb_contact_surface`, `myb_contact_law`, `myb_contact_point`) from `myb/contact.h`.

## Inputs and outputs

- Radii and the watched distance in m. A radius of 0 (the default) makes a body a point: it
  takes part only in closest approaches, when `watch` > 0.
- `myb_system_event`: the kind, the pair a < b, the time t of the integrator (left there) and the
  distance of the centres in m. A contact is the instant at which the distance falls to the sum
  of the radii while the bodies approach.
- `myb_system_merger`, with G as a factor like the invariants of `myb/dyn.h`: the energy of the
  relative motion G m_a m_b (|v|²/(2(m_a + m_b)) − 1/|r|) in m⁵ s⁻⁴, its angular momentum
  G m_a m_b/(m_a + m_b) r × v in m⁵ s⁻³ (r = x_b − x_a, v = v_b − v_a), and the spin of the
  remnant — that angular momentum plus the spins brought: G m (I ω) of integrated rotations in
  inertial axes and χ G² m²/c of a pair beyond 1PN.
- A merger: a takes the mass of both, the position and velocity of their centre of mass and the
  radius of a sphere of both volumes, ∛(r_a³ + r_b³); it keeps its own model; b is absent from
  then on (its mass is 0, its part of the state means nothing) and every level that names it is
  dropped. A rotating a keeps its attitude and its shape scaled to the new radius (moments per
  unit mass as the square, the radius of its distortion as the first power of the ratio of the
  radii), turns as a rigid whole, the core with the mantle, and takes the whole spin of the
  merger. A PN pair that names a or b ends: its terms stop and its spins are set to zero in the
  state. The integrator is restarted.
- A rebound reverses the component of the relative velocity along the line of centres and
  multiplies it by the restitution (1 elastic, 0 they stay together); the momentum is kept; the
  integrator is restarted.
- The contact level: the surface is at rest in the frame of the ground, about its centre; the
  points are spheres in the frame of the body, about its centre; the law is per unit mass of the
  body. The friction acts on the velocity of each point relative to the ground under it (a
  turning ground carries its surface along); the ground is pushed back by the ratio of the masses
  at the same points (a ground of gm 0 is not pushed); the torques turn the bodies whose rotation
  is integrated. A body without a frame keeps the axes of the system and does not turn. A
  friction that sticks (`MYB_CONTACT_FRICTION_STICK`) adds 3 coordinates a point to the state,
  after every other kind: their rates are the deflections of the springs in the axes of the
  surface (`myb_system_pack_contact`, `myb_system_unpack_contact`; zero — the springs at rest),
  their accelerations the rates of the deflections; a merger that drops the contact leaves them
  in the state, meaning nothing.

## Domain of validity

- Refused with `MYB_E_DOMAIN`: a radius negative or not finite; `watch` negative or not finite;
  in `myb_system_merge` wrong indices, an absent body, b with an integrated rotation and a
  without (merge the other way); in `myb_system_bounce` a restitution outside 0…1 or coincident
  centres; in `myb_system_add_contact` the same body twice, an invalid surface, law or point, no
  points. More than 16 contacts — `MYB_E_OVERFLOW`.
- `myb_system_advance` needs an integrator started with `myb_system_acceleration` and this
  system; other errors are those of `myb_ias15_step`.
- A pair already in contact when the call starts (an event left alone) is not reported again
  until it has separated; a closest approach is not reported twice.
- A closest approach is found where the range rate turns from negative to non-negative inside a
  step: a pass with its minimum inside a step is not missed.
- A system with contacts is integrated with DOP853, `myb_system_acceleration` taken as a
  first-order system (x, v)′ = (v, a): IAS15 needs its steps ended at the touches, which the
  system does not do (a term of a force model under `myb_orbit` goes with IAS15:
  [contact and landing](contact-landing.en.md)); the events of `myb_system_advance` are found
  with IAS15 only. `myb_system_partials` takes a system with contacts by differences of whole
  propagations ([system of bodies](system-of-bodies.en.md)).

## Accuracy

In `test_system` ([tests/unit/test_system.c](../../tests/unit/test_system.c)); figures from the
plan (phases 18 and 22) where it gives them, otherwise the test's bound:

| What | Reference | Result |
|---|---|---|
| Two spheres (60 m and 40 m) let go at rest 5000 m apart | radial Kepler fall in closed form | contact at 554.689657182 s, 10⁻¹² relative; centres 100 m − 9·10⁻¹² m apart |
| Elastic rebound, then the same time | the start | 5000.000000000 m apart, at rest |
| Merger after the second fall | mass, momentum; the energy of the pair at rest, −G m_a m_b/d | energy to 9 digits; angular momentum 0; radius ∛(60³ + 40³) to 10⁻¹² m |
| Fly-by, e = 1.6, pericentre 7000 km | the pericentre of the hyperbola | 10⁻¹¹ in time, 10⁻¹³ in distance |
| The same with spheres reaching 7500 km | the hyperbola at r = 7500 km | test bound: 10⁻¹¹ in time, 10⁻⁶ m in distance |
| Spinning body hit off-centre | total angular momentum, orbit and spins | 1.4·10⁻¹⁶; the remnant's spin is the reported one to 10⁻¹⁴ (bound) |
| PN pair with spins merging | orbital angular momentum + χ G² m²/c of each | exact (plan); test bound 10⁻¹⁴ |
| 300 spheres of 2 m in straight lines (the neighbour search) | the billiard worked out pair by pair in closed form | the same first five contacts and bounces, times to 5.6·10⁻¹³ s |

The contact level, in `test_contact` ([tests/unit/test_contact.c](../../tests/unit/test_contact.c),
integrated with DOP853):

| What | Reference | Result |
|---|---|---|
| Lander on four legs, slope 15°, μ = 0.6 | creep speed of the smooth law in closed form | 6·10⁻¹³; the legs carry its weight to 10⁻⁹ |
| The same on 40° | g (sin α − μ cos α) | 9·10⁻⁸ |
| Spinning ball off a turning ellipsoid | conservation | momentum 2·10⁻¹⁶, angular momentum 1.2·10⁻¹³ |
| The lander with the friction that sticks (k_t = 1000 s⁻², c_t = 16 s⁻¹ a leg), 15° | at rest | 9.3·10⁻¹⁴ m/s after 20 s (the smooth law creeps at 5.0·10⁻⁴ m/s); the legs bear its weight to 10⁻¹¹ |
| The same on 40° | g (sin α − μ cos α) | 8.3·10⁻⁸; the energy never grows by more than 4.4·10⁻¹⁶ m²/s² |
| A block on 25° with the friction that sticks | rest at the deflection of its weight | to 10⁻⁹ of it |

## Cost

With 1500 spheres the neighbour search gives the same events 3.5 times faster than the loop over
all pairs: 7.3 s against 22.8 s (plan, phase 22); what remains is the direct O(N²) sum of the
gravity.

## Limitations

- For events the bodies are spheres; shapes are not used (a collision with a shape is in
  [events and geometry](events-geometry.en.md)).
- A merger is perfectly inelastic: the energy it takes goes into heat, which the system does not
  follow. The spin of a body whose rotation is not integrated is not known and counts as none;
  the distortions of the inertia are left out of the spins brought.
- A rebound is that of smooth hard spheres: no friction between them, no change of their spins.
- A merged body keeps its place in the state: the number of coordinates and the cost do not
  fall.
- The neighbour search cuts the pairs examined, not the gravity: that stays the direct sum
  (O(N²)) unless the approximate FMM level of [N bodies](nbody.en.md) is chosen.
- The contact level is integrated with DOP853, the event search of `myb_system_advance` with
  IAS15: the two are not combined, and contacts of a system under IAS15 are not done (the system
  does not end its steps at the touches). With the smooth friction a body below the angle of
  friction creeps; with the friction that sticks it stands.

## Data

None.
