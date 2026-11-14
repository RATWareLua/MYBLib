# Low thrust: burn arcs, Edelbaum and the Q-law

**English** | [Русский](low-thrust-qlaw.ru.md) · [All passports](index.en.md)

## What it is

Continuous thrust as a term of a force model. Burn arcs of constant thrust with a linear mass
flow and a steering law (along the velocity, inertial, fixed in the radial–transverse–normal
frame, Edelbaum's yaw law); the Gauss equations in modified equinoctial elements; Edelbaum's
circle-to-circle transfer in closed form; and the Q-law, a feedback law that steers towards any
target orbit without an optimisation or a first guess — as a term of a force model, or flown
with the mass as a seventh component of the state so that the coasts it decides on are counted.
The mass left after the burns can feed the surface forces of the object.

## Source

- Burn arcs: the rocket equation (dv = c ln(m0/m1), exact for the acceleration integrated along
  an arc); c = Isp g0 with g0 = 9.80665 m/s².
- Edelbaum's transfer and yaw law: J. Edelbaum, *ARS J.* 31, 1079 (1961); the yaw law as in
  J. A. Kechichian, *J. Guid. Control Dyn.* 20, 803 (1997).
- Gauss equations in modified equinoctial elements: M. J. H. Walker, B. Ireland, J. Owens,
  *Celest. Mech.* 36, 409 (1985).
- The Q-law and its effectivity: A. E. Petropoulos, AIAA 2004-5089.
- Implemented from the publications; no third-party code, no data files. The optimal transfers
  are in [optimal control](optimal-control.en.md).

## Interface

Header [`myb/thrust.h`](../../include/myb/thrust.h):

- burn arcs: `myb_thrust_init`, `myb_thrust_mass`, `myb_thrust_mass_at`,
  `myb_thrust_acceleration` (a `myb_term_function`), `myb_thrust_rtn`; structures
  `myb_thrust_arc`, `myb_thrust`; enum `myb_steering`; `MYB_THRUST_MAX_ARCS` (16);
- `myb_equinoctial_rates`; `myb_edelbaum`, `myb_edelbaum_state`;
- the Q-law: `myb_qlaw_direction`, `myb_qlaw_arrived`, `myb_term_qlaw`, `myb_qlaw_flight_rhs`,
  `myb_qlaw_fly`; structures `myb_qlaw`, `myb_qlaw_term`, `myb_qlaw_flight`, `myb_qlaw_arrival`.

## Inputs and outputs

- Thrust in N, exhaust velocity in m/s, masses in kg, times in s on the time of the force
  model; x and v relative to the central body in its inertial axes.
- A burn arc [start, stop]: the mass falls at T/c, the acceleration T/m(t) grows; the mass is a
  function of time, so the state of the integrator stays (x, v). Arcs may not overlap.
- Steering: `MYB_STEER_VELOCITY`, `MYB_STEER_ANTI_VELOCITY`, `MYB_STEER_INERTIAL` (a fixed
  direction), `MYB_STEER_RTN` (fixed radial, transverse, normal components), `MYB_STEER_EDELBAUM`
  (thrust at β out of the plane, tan β = v0 sin β0 / (v0 cos β0 − Δv), Δv = c ln(m(start)/m(t)),
  the out-of-plane sign following cos u, the inclination changing in the sense given).
- `myb_equinoctial_rates`: the rates of p, f, g, h, k, L for an acceleration in R, T, N
  (Walker et al.).
- `myb_edelbaum`: ΔV = √(v0² − 2 v0 v1 cos(π Δi/2) + v1²) and β0 for circular speeds v0, v1 and
  the plane change Δi; `myb_edelbaum_state`: the speed and Δi after an increment Δv.
- The Q-law: the target a (m), e, i, Ω, ω (rad) with weights (0 leaves an element free), an
  optional lowest periapsis radius with P = exp(100 (1 − r_p/rp_min)). Q is the distance in each
  element over the largest rate at which thrust of unit acceleration can change it anywhere on
  the orbit — about the square of the time still needed (s²); the thrust points where Q falls
  fastest. `myb_qlaw_direction` returns the unit R, T, N direction, Q, its rate (s, never
  positive) and the effectivity — that rate over the best at any of 72 places of the orbit, in
  [0, 1].
- `myb_term_qlaw`: an engine of constant thrust that points as the law says and is off where the
  effectivity is below the threshold (0.2 … 0.5 is usual; 0 — always on) and once every element
  is within its tolerance; the mass falls from t0 as if it never paused.
- `myb_qlaw_flight_rhs`, `myb_qlaw_fly`: y = (x, v, m), the mass falling only while the engine
  burns; the field without the engine is the point mass of the law or any force model; DOP853
  (rtol 10⁻¹⁰, atol 10⁻⁴ in m, m/s, kg), the arrival looked for every `step` seconds.
- `myb_thrust_mass_at`: the mass as a `myb_object.mass_function`
  ([`myb/env.h`](../../include/myb/env.h)), so that the surface forces act on the mass that is
  left.

## Domain of validity

- Edelbaum's relations are orbit-averaged between circular orbits: their relative accuracy is of
  the order of the ratio of the thrust to the gravitational acceleration.
- The Q-law needs an ellipse; it is singular where a targeted angle is not defined (ω on a
  circle, Ω in the equator).
- The engine of the Q-law switches sharply: an integrator with steps of its own choosing
  (`MYB_ORBIT_DOP853` of `myb/orbit.h`) is the one to use.
- Refused with `MYB_E_DOMAIN`: a mass ≤ 0, more than 16 arcs, an invalid, overlapping or
  unordered arc, a zero direction, more propellant than mass; speeds not positive (Edelbaum);
  an invalid law, an orbit that is not an ellipse, non-finite input; the mass exhausted, a
  threshold outside [0, 1); t_limit ≤ t0, a step ≤ 0, a mass not positive (`myb_qlaw_fly`).

## Accuracy

`test_thrust` ([tests/unit/test_thrust.c](../../tests/unit/test_thrust.c)); the Q-law with 1 N on
1000 kg at 30 km/s:

| What | Reference | Result |
|---|---|---|
| A burn in free space | c ln(m0/m1) | 1.7·10⁻¹⁴ |
| Constant radial acceleration 0.9 μ/(8 r0²) | the closed-form apocentre (μ − √(μ² − 8aμr0²))/(4ar0) | 3·10⁻¹³ |
| 10 days of RTN thrust, Gauss equations in MEE | Cowell (IAS15) | 0.44 mm at 8617 km |
| Edelbaum's transfer 7000 → 7200 km with 1° of plane, flown | 7200 km, i = 27.5° | a = 7200.001 km, i = 27.5025° (thrust/gravity = 2.4·10⁻⁴) |
| Q-law, 7000 → 10 000 km | Edelbaum 1232.6 m/s | 1234.3 m/s (+0.14 %) in 14.00 days |
| Q-law, the same with the plane turned by 0.1 rad | Edelbaum's optimum 1641 m/s | 1746 m/s (6 % above) |
| Q-law, ellipse e = 0.725 → circle of 42 164 km, rp_min 6600 km, threshold 0.3 | — | 26.6 days; the periapsis never below 6738 km |
| The same with the mass in the state | the term without memory, 2389 m/s | 26.5 days, the engine on for 22.6; 2020 m/s spent |
| An engine without pauses, mass in the state | the rate of the engine | to the step of the look (300 s) |

`test_nongrav` ([tests/unit/test_nongrav.c](../../tests/unit/test_nongrav.c)): with
`myb_thrust_mass_at` the light pushes harder on what is left after a burn — the ratio of the
accelerations equals the ratio of the masses to 10⁻¹⁴.

## Cost

Without the effectivity threshold an engine that always burns dithers near the target: 100 km
from it, half a million integrator steps a day (the finding behind the threshold).

## Limitations

- Burn arcs have constant thrust (no throttling) and at most 16 arcs; the mass is a function of
  time, not a state.
- The Q-law is a feedback law, not an optimum: 0.14 % above Edelbaum coplanar, 6 % above with a
  plane change. It targets the five elements a, e, i, Ω, ω — not the place on the orbit.
- `myb_term_qlaw` has no memory of when it burned: with a threshold it overestimates the mass
  spent (2389 against 2020 m/s above); `myb_qlaw_fly` carries the mass in the state.
- Edelbaum's relations hold between circles in the orbit-averaged approximation only.

## Data

None.
