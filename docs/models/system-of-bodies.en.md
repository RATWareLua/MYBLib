# System of bodies with models

**English** | [Русский](system-of-bodies.ru.md) · [All passports](index.en.md)

## What it is

A system of massive bodies that are more than points: the point masses of `myb_nbody` (their
graph, the 1PN terms, the threads) with further levels of the interaction graph asked for by
name — the field of a body acting on another taken as a point, the quadrupoles of two bodies,
the tide by band with lags, the distortion of a body by a tide and by its own spin, the terms
beyond 1PN of a pair, actions given as functions (light on a mesh among them), the dragging of
the frame by a spin, a mass that changes, bodies whose motion an ephemeris gives. The rotation of
chosen bodies (a mantle that may hold a fluid core) is integrated together with the motion, by
one IAS15 integrator. The criterion of the model is DE440 integrated again from its initial
state, the Moon with its librations. Encounters, collisions, mergers and contact are in
[encounters in the system](system-encounters.en.md).

## Source

- The rotation of a mantle with a fluid core, its distortion and the tide by band with orbital
  and rotational lags: the model of the lunar ephemerides — W. M. Folkner, J. G. Williams,
  D. H. Boggs, R. S. Park, P. Kuchynka, IPN Progress Report 42-196 (2014), sections III and IV,
  eqs. 31–32 (the tide), eq. 36 (the moments of the Moon in the test); R. S. Park et al.,
  *Astron. J.* 161, 105 (2021).
- The dragging of the frame: J. Lense, H. Thirring (1918), through `myb_rel_lense_thirring`
  (see [1PN relativity](relativity-1pn.en.md)).
- The small bodies of DE440 in the test: the file of perturbers sb441-n373s.bsp (D. Farnocchia,
  JPL IOM 392R-21-005, 2021) and the ring of the Kuiper belt of Park et al. (2021, section 2.2).
- Implemented from the papers; no third-party code. The levels reuse the library's own models:
  [spherical harmonics](spherical-harmonics.en.md),
  [figure–figure](figure-figure-rotation.en.md), [PN binaries](pn-binary.en.md),
  [radiation on surfaces](radiation-surfaces.en.md), [IAS15](ias15.en.md).

## Interface

Header [`myb/system.h`](../../include/myb/system.h):

- the system: `myb_system_create`, `myb_system_destroy`, `myb_system_bodies`,
  `myb_system_nbody` (the `myb_nbody` inside, see [N bodies](nbody.en.md));
- the model of a body: `myb_system_set_orientation`, `myb_system_set_field`,
  `myb_system_set_rotation`, `myb_system_set_tide`, `myb_system_set_mass_rate`;
- the levels: `myb_system_add_field`, `myb_system_add_figures`, `myb_system_add_tide`,
  `myb_system_add_distortion`, `myb_system_add_pn_pair`, `myb_system_add_action`,
  `myb_system_light_action`, `myb_system_add_frame_dragging`, `myb_system_add_guided`,
  `myb_system_guided`;
- the state and the integration: `myb_system_coordinates`, `myb_system_pack_rotation`,
  `myb_system_unpack_rotation`, `myb_system_pack_spins`, `myb_system_unpack_spins`,
  `myb_system_integrator`, `myb_system_acceleration`, `myb_system_inertia`;
- the derivatives: `myb_system_partials`;
- structures `myb_system_rotation`, `myb_system_tide`, `myb_system_pn_pair`,
  `myb_system_view`, `myb_system_light`; the callback `myb_system_action`; the limits
  `MYB_SYSTEM_MAX_*`.

## Inputs and outputs

- SI units: gravitational parameters in m³/s², positions in m, velocities in m/s, times in
  seconds of the integration (providers read them as `myb_time_from_seconds`). The frame and the
  origin of the system are the caller's (the DE440 test: the barycentre and the axes of DE440).
- `myb_system_acceleration` has the signature of `myb_accel_function`: n =
  `myb_system_coordinates` coordinates in blocks of three, a state holds the coordinates and then
  their rates (2n values) — 3 per body (m); 9 per integrated rotation (the attitude quaternion of
  [`myb/attitude.h`](../../include/myb/attitude.h), 2 unused values, 3 values whose rate is the
  angular velocity of the core in the mantle frame); 6 per PN pair with spins (values whose
  rates are the two dimensionless spins). The quaternion is a second-order coordinate:
  q″ = (q′ω + qω′)/2, ω = 2 q̄ q′ the angular velocity of the mantle in its own frame.
- An orientation from outside gives the rotation from the inertial frame of the system to the
  body. An action writes the acceleration (m/s², inertial axes) and the torque per unit mass
  (m²/s², body axes); the view gives the frame as the rotation inertial → body, 9 values by rows.
- Moments of inertia are per unit mass (m²): A ≤ B ≤ C of the undistorted body; the mantle is
  diag(A − A_c, B − A_c, C − C_c) plus the tidal part −k₂R⁵ Σ (m_T/m)(r rᵀ − |r|²/3)/|r|⁵ (r at
  t − lag in the mantle frame of t − lag) and the spin part k₂R⁵/(3Gm)(ω ωᵀ − (|ω|² − n²)/3 −
  n² z zᵀ); the core is coupled by N = k_v(ω_c − ω) + (C_c − A_c)(z·ω_c)(z × ω_c). The degree 2 of
  the field of such a body comes from this tensor, the degrees from 3 from its field.
- The tide by band (j = 0 long-period, 1 diurnal, 2 semi-diurnal, about the third axis of the
  body frame): Love numbers k₂ⱼ, orbital lags in s (may be negative), rotational lags in s (the
  long-period one is not used), the rotation rate in rad/s. The raiser is taken where it was,
  from a Taylor series of the motion at t to the third order.
- The mass rate is relative, in 1/s: gm(t) = gm (1 + rate (t − epoch)). The spin of frame
  dragging is the angular momentum per unit mass (m²/s) in the axes of the system, copied.
- `myb_system_partials`: phi is 2n × 2n, row-major, ∂y₁ᵢ/∂y₀ⱼ; the sensitivity 2n × count,
  ∂y₁ᵢ/∂gm_k; steps of 10⁻⁶ of the largest coordinate of each kind (or rate), or of the gm.

## Domain of validity

- Limits per system: 16 rotations, 256 fields, 32 pairs of figures, 64 tides, 32 distortions,
  16 PN pairs, 64 actions, 64 draggings; beyond them `MYB_E_OVERFLOW`.
- `MYB_E_DOMAIN` when setting up: wrong indices, the same body twice, a body without mass where
  one is needed, a field degree that is negative, 1 or above the field's, a rotation for a body
  with an orientation from outside, non-positive moments, k₂ without a radius, a distortion of a
  body whose rotation is not integrated, unknown PN orders, no speed of light, `MYB_PN_3` or `MYB_PN_4` without r0, a
  negative κ, frame dragging before `myb_nbody_set_relativity`, a negative guided gm, values that
  are not finite. A tide is accepted only for a body whose frame is given from outside (the code
  refuses it for a body whose rotation is integrated: that body takes the tide through its
  distortion).
- `MYB_E_DOMAIN` from `myb_system_acceleration`: a level that needs a frame the body does not
  have, coincident bodies, an inertia tensor that cannot be inverted, a mass that would become
  negative, a guided place that is not finite or lies on a body; the statuses of providers and
  fields are passed on.
- The tails of `MYB_PN_4` and `MYB_PN_4RR` need the pair's osculating orbit to stay bound: the
  integration fails otherwise.
- `myb_system_partials`: valid while the motion does not branch (no encounter between t0 and
  t1); t1 = t0, an index or a gm of 0, non-finite input — `MYB_E_DOMAIN`.
- A system is used by one thread at a time.

## Accuracy

Levels, in `test_system` ([tests/unit/test_system.c](../../tests/unit/test_system.c)); the
figures are those of the plan (phases 18 and 22), the test bounds are looser:

| What | Reference | Result |
|---|---|---|
| Field of a body on two points | the same sum by hand | the same bits; momentum 2·10⁻¹⁶ |
| Tide by band, no lag | gradient of the potential by differences | 7·10⁻¹¹ |
| Earlier position of the raiser (900 s, a 27-day orbit) | Kepler motion | 3·10⁻¹¹ |
| Lagging part of the tide (lag angle 4·10⁻³) | Mignard's first order (`myb_force_tidal_lag`) | 5·10⁻³ (the order of the lag angle) |
| Free top over 30 days (some 25 turns) | exact torque-free motion | 3·10⁻¹³; \|q\| − 1 = 5·10⁻¹⁵ |
| Figure turned by a point mass, 40 days | orbit + spin angular momentum | exchange 4·10⁻⁵, sum kept to 4·10⁻¹⁶ |
| PN pair among three bodies, with spins | `myb_pn_binary_evaluate` | 9·10⁻¹⁷ |
| Light on a mesh: the torque through the equations of rotation | the mesh itself | 8·10⁻¹⁵ |
| Light with thermal inertia (`test_surface`) | θ = 0: the mesh of the moment; mean over a turn: `myb_force_yorp` | 1.5·10⁻¹⁵; 3.5·10⁻¹⁶ (θ = 0, 0.7, 3) |
| Derivatives, two bodies over 4000 s | closed transition matrix (`myb_kepler_stm`) | 1.8·10⁻¹¹ of the largest element; ∂/∂gm 6·10⁻¹¹ |

Bounds checked by the test without a figure in the plan: a massless body about a guided centre
returns after a period to 10⁻¹¹ of the radius; a guided planet pulls as the same planet
integrated to 10⁻¹⁸ m/s² (of 6·10⁻³ m/s²); the mass rate scales the pull by 1 − 9·10⁻⁸ after a
million years to 10⁻¹⁵; frame dragging against `myb_rel_lense_thirring`, the momentum kept, to
10⁻¹⁷ m/s².

DE440 integrated again (`test_de440_system`,
[tests/unit/test_de440_system.c](../../tests/unit/test_de440_system.c)): only the state at
JD 2440400.5 is taken (positions from de440.bsp, the angles and angular velocities of the Moon
and its core from header.440); the model is EIH, J₂ of the Sun on every body, J₂…J₅ of the Earth
and the field of the Moon to degree 6, the quadrupoles of the Earth and the Moon, the tides of
the Earth from the Moon and the Sun on the Moon with the five lags of DE440, the rotation of the
Moon with its core. Planets heliocentric, the Moon geocentric, the Sun barycentric; differences
after 5 years (plan, phase 20; the small bodies are 343 asteroids and 30 Kuiper-belt bodies as
guided bodies, and a ring of 36 equal masses at 44 au):

| Body | without the small bodies | with the small bodies, the loss of mass of the Sun and Lense–Thirring |
|---|---|---|
| Sun (from the barycentre) | 1160 m | 5.7 cm |
| Mercury | 22 m | 2.1 cm |
| Venus | 82 m | 1.0 cm |
| Earth | 43 m | 21 cm |
| Moon (from the Earth) | 15 cm | 4.9 cm |
| Mars | 4022 m | 1.3 cm |
| Jupiter, Saturn | 2440 m | 0.46 and 0.64 m |
| Uranus, Neptune, Pluto | 1626–2789 m | 1.2; 6.1; 4.2 m |

Before the small bodies, the loss of mass and Lense–Thirring (phase 18): the Moon 1.1 cm after
30 days, 8.8 cm after a year, 19 cm after 5 years (54 m without the tides of the Earth); Mercury
0.77 m and Venus 0.56 m after a year; the attitude of the Moon against moon_pa_de440_200625.bpc
0.003″, 0.041″ and 0.16″ (1.4 m at the surface). The rate of the rotational deformation moves
the attitude of the Moon by 0.0003″ and its position by 1 cm in 5 years. The test bounds the Moon
by 5 cm, 1 m and 10 m, its attitude by 0.06″ (1 year) and 0.5″ (5 years); with the small bodies
Mercury, Venus and Mars by 2 cm after a year and 10 cm after 5 years.

## Cost

- The DE440 model over five years: 8 s for both runs, with and without the tides of the Earth
  (1828 steps, 77 498 evaluations; plan, phase 18).
- A guided body: its provider is called once per instant of a step, 8 times a step in place of
  20–40 (plan, phase 20).
- `myb_system_partials`: two propagations per column (header).

## Limitations

- The tide of a body whose rotation is integrated is not a level: such a body answers only
  through the distortion of its inertia tensor, of degree 2, with one k₂ and one lag (not by
  band); its degrees from 3 are rigid.
- Earlier positions for the lags come from a Taylor series to the third order (the jerk of the
  two bodies alone); no history is kept.
- Frame dragging keeps the spin it was given (a pole that stays put); the mass rate is linear in
  time and what leaves does not push the body.
- A guided body feels nothing: the momentum of the system is not kept, by as much as the bodies
  would have moved it.
- The light of the moment with thermal inertia takes the rotation as steady over a turn.
- In the DE440 run: the pressure of light on the Earth–Moon system and Vondrák's precession are
  not modelled (the Earth 21 cm); the phase of the points of the Kuiper-belt ring is not
  published (the outer planets 0.46–6.1 m); the pole of the Earth is the CIP of IAU 2006/2000A,
  where DE440 has a simpler model.
- The derivatives are central differences of whole propagations: a part in 10⁹…10¹¹ of a
  derivative, two propagations per column, wrong across an encounter.
- Phase 21 adds: the variational equations of the whole system (in place of differences of
  propagations).

## Data

None for the library. The DE440 test needs, from
[`cmake/FetchData.cmake`](../../cmake/FetchData.cmake) (NASA/JPL products, checked by SHA-256,
not part of the repository): de440.bsp, header.440, moon_pa_de440_200625.bpc; the run with the
small bodies also needs sb441-n373s.bsp (982 MB, only with `-DLARGE=ON`; skipped without it).
The pole of the Earth uses the in-repo kernel `data/frames/iau2006-cip.tk`. Without the
downloads the test is skipped.
