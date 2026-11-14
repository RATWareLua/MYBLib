# Figure–figure coupling and rotational dynamics

**English** | [Русский](figure-figure-rotation.ru.md) · [All passports](index.en.md)

## What it is

Two parts of the dynamics of extended bodies. The mutual potential of two figures beyond the
terms "the field of one at the centre of the other": its first term that couples the two
figures, U₂₂ between their quadrupoles, with the accelerations of both centres of mass and the
torques. And the rotation of rigid bodies: the conventions of attitude (quaternions, the twelve
Euler sequences), principal axes, Euler's equations for the body angular momentum, the
gravitational torques of a point mass, of a whole field and of two quadrupoles, the exact
torque-free top, a symplectic splitting integrator, and the equations of a mantle with a fluid
core — the model of the Moon of the JPL ephemerides, checked on its physical librations against
DE440.

## Source

- U₂₂: the STF multipole expansion of the mutual potential, T. Damour, M. Soffel, C. Xu,
  *Phys. Rev. D* 43, 3273 (1991).
- The torque-free top: L. D. Landau, E. M. Lifshitz, *Mechanics*, §37; E. Celledoni, F. Fassò,
  N. Säfström, A. Zanna, *SIAM J. Sci. Comput.* 30, 2084 (2008) — Jacobi elliptic functions and
  the incomplete integral Π of the library.
- The splitting: J. Touma, J. Wisdom, *Astron. J.* 107, 1189 (1994), with Yoshida's weights.
- The fluid core: H. Poincaré 1910; W. M. Folkner, J. G. Williams, D. H. Boggs, R. S. Park,
  P. Kuchynka, IPN Progress Report 42-196 (2014), §IV.
- The quaternion of a matrix: Shepperd's method.
- Implemented from the papers; no third-party code.

## Interface

Header [`myb/grav.h`](../../include/myb/grav.h): `myb_quadrupole_from_coefficients`,
`myb_quadrupole_rotate`, `myb_figure_figure`.

Header [`myb/attitude.h`](../../include/myb/attitude.h):

- quaternions: `myb_quat_multiply`, `myb_quat_conjugate`, `myb_quat_normalize`,
  `myb_quat_to_matrix`, `myb_quat_from_matrix`, `myb_quat_rotate`,
  `myb_quat_from_rotation_vector`, `myb_quat_to_rotation_vector`, `myb_quat_slerp`,
  `myb_quat_derivative`;
- Euler angles: `myb_euler_to_matrix`, `myb_euler_from_matrix`; `myb_principal_axes`;
- equations of motion: `myb_rigid_body_rates`, `myb_core_mantle_rates`;
- torques: `myb_torque_point_mass`, `myb_torque_field`, `myb_torque_figure_figure`;
- the exact top: the structure `myb_free_rotation`, `myb_free_rotation_init`,
  `myb_free_rotation_state`;
- splitting: the callback `myb_torque_function`, `myb_rotation_split`.

A rotation integrated together with the motion, figures as a level of the graph
(`myb_system_set_rotation`, `myb_system_add_figures`, `myb_system_add_distortion`):
[system of bodies](system-of-bodies.en.md). The family of splitting integrators:
[symplectic integrators](symplectic-integrators.en.md). The torque of a lagging tide
(`myb_torque_tidal_lag`): [lagged tides](tides-lagged.en.md).

## Inputs and outputs

- An attitude is the frame rotation M from the inertial frame to the body frame
  (v_body = M v_inertial; the rows of M are the body axes), as `myb_rotation.m`. The quaternion
  q = (w, x, y, z) is scalar first, Hamilton product, M(q) v = q̄ v q; a frame turned by a about
  the unit axis n has q = (cos a/2, n sin a/2), M = Rₙ(a). Turns about the body's own axes
  compose as q = q₁ q₂, M(q₁ q₂) = M(q₂) M(q₁). The body angular velocity ω gives
  dM/dt = −[ω×] M and dq/dt = q (0, ω)/2. `myb_quat_from_matrix` returns w ≥ 0.
- Euler angles: axes (i, j, k) of the twelve sequences (proper, e.g. 3-1-3, and Tait–Bryan,
  e.g. 3-2-1); M = R_k(c) R_j(b) R_i(a). From a matrix: b in [0, π] (proper) or [−π/2, π/2]
  (Tait–Bryan), a and c in (−π, π]; at gimbal lock a = 0.
- Principal axes: moments ascending, `axes` the rotation from the given frame to the principal
  frame (rows, det +1).
- Dynamics in principal axes with the body angular momentum L = I ω as the state:
  dL/dt = L × ω + torque. Inertia and torque absolute (kg m², N m) or per unit mass of the body
  (m², m²/s²); the gravitational torques come per unit mass.
- Torques in the body frame; r is the attracting body relative to the centre of mass of the
  torqued body (m). Point mass: 3 gm/r⁵ (r × I r), the quadrupole term. Whole field:
  −(m_other/m_body) r × a, a the acceleration of the body's own field at r (e.g. by
  `myb_gravity_field_eval`). Two quadrupoles: −2 axial(∂U/∂Q_A Q_A − Q_A ∂U/∂Q_A) / gm_a.
- Quadrupoles: trace-free Q = ∫ ρ (3 x xᵀ − |x|² 1) dV = tr(I) 1 − 3 I, passed as G Q
  (m⁵/s²), row-major; from the fully normalised C̄₂₀, C̄₂₁, C̄₂₂, S̄₂₁, S̄₂₂ with GM and the
  reference radius; into inertial axes as gq = mᵀ gq_body m.
- U₂₂ = G/(36 r⁵) [105 (n·Q_A·n)(n·Q_B·n) − 60 n·Q_A Q_B·n + 6 tr(Q_A Q_B)], r = x_B − x_A;
  `myb_figure_figure` adds the accelerations of both centres of mass (in the frame of r).
- The core: a mantle of inertia I_m (body frame; its rate dI_m or NULL) and a core of moments
  (A_c, A_c, C_c) fixed in the mantle, coupled by
  N_cmb = k_v (ω_c − ω_m) + (C_c − A_c)(z·ω_c)(z × ω_c); the rates of ω_m and ω_c in the
  mantle frame.

## Domain of validity

- U₂₂ is of relative size (size_A size_B / r²)² of the point-mass term; the next terms are of
  relative order (size/r)⁶ for bodies without octupoles, 5 otherwise.
- `MYB_E_DOMAIN`: `myb_figure_figure` and `myb_torque_figure_figure` with r = 0, input not
  finite or gm ≤ 0; an invalid Euler sequence; an inertia tensor that is not symmetric or not
  finite; a singular mantle tensor or core moments that are not positive.
- The exact top: moments positive and finite (l = 0 is allowed); the separatrix
  (|L|² = 2 T I₂ with the motion not steady) is refused with `MYB_E_DOMAIN`, and the errors grow
  near it like those of the Jacobi functions.
- `myb_rotation_split`: orders 2, 4, 6, 8, steps ≥ 0, torques that depend on the time and the
  attitude but not on L; `MYB_E_DOMAIN` otherwise.

## Accuracy

Tests: [`test_grav`](../../tests/unit/test_grav.c),
[`test_attitude`](../../tests/unit/test_attitude.c),
[`test_moon_libration`](../../tests/unit/test_moon_libration.c),
[`test_de440_rel`](../../tests/unit/test_de440_rel.c).

| What | Reference | Result |
|---|---|---|
| U₂₂ on two clouds of point masses without odd moments, 20 trials | the exact double sum (`test_grav`) | the error of the mutual force falls 29 times at 12 sizes and 95 times at 24; the rest falls as r⁻⁸ |
| Q from C̄₂ₘ, S̄₂ₘ, and its rotation | degree 2 of the field (`test_grav`) | 10⁻¹⁴ (bound) |
| Quaternions: matrix round trip, composition, exp(log) | identities (`test_attitude`) | 5·10⁻¹⁶, 10⁻¹⁵, 5·10⁻¹⁶ (bounds) |
| 12 Euler sequences: matrix round trip, angles, gimbal lock | identities (`test_attitude`) | 10⁻¹⁵, 10⁻¹³, 10⁻¹⁵ (bounds) |
| Principal axes | diagonalisation (`test_attitude`) | 2·10⁻¹⁵ (bound) |
| Torque of two quadrupoles | the derivative with respect to a rotation; the balance of the angular momentum of the bodies and the orbit (`test_attitude`) | 10⁻⁸ (bound); 6.6·10⁻¹⁶ |
| Exact top: circulation about the largest and the smallest axis, symmetric bodies, steady rotation | DOP853 at tolerance 10⁻¹⁴ (`test_attitude`) | 2·10⁻¹⁴ rad; the inertial angular momentum constant to 4·10⁻¹⁶ |
| Symmetric top | closed form (`test_attitude`) | 7·10⁻¹⁵ rad |
| Splitting | exact free motion; a gravity-gradient torque (`test_attitude`) | exact without torque (10⁻¹² bound); observed orders 2.00 / 3.99 / 5.98 / 7.95 |
| Lunar librations over two years from the DE440 initial conditions | `moon_pa_de440_200625.bpc` (`test_moon_libration`) | 0.040″ the first year, 0.071″ over two (0.6 m at the surface); angular velocity 6.6·10⁻⁸ rad/day |
| DE440 with figures: EIH, J₂…J₅ of the Earth, the Moon to degree 4, J₂ of the Sun, Earth–Moon U₂₂; one year | DE440 (`test_de440_rel`) | Moon 0.65 m (23 km without figures), Mercury 14 m (106 m without the solar J₂), Earth 124 m (270 m) |

Every part of the libration model is needed: without the tidal deformation 0.64″, without its
rate 0.97″, without the friction of the core 0.11″, without U₂₂ 0.085″, without the rotational
deformation 0.077″, and with the tide taken in the mantle frame of t instead of t − τ — 29″.
After ten years of the DE440 model with figures the Moon is off by 239 m — the secular tidal
acceleration of DE440 (25.8″/cy² gives 240 m), which that model has no tides for.

## Cost

Not measured.

## Limitations

- Only the l = l′ = 2 coupling of two figures is given; the higher couplings (octupole with
  quadrupole and beyond) are not. The terms l′ = 0 come from a field of any degree.
- `myb_torque_point_mass` is the quadrupole term from the inertia tensor; `myb_torque_field` is
  exact for the whole field.
- The splitting integrator takes torques from a potential only (no dependence on L) and a fixed
  step.
- The core is axisymmetric (A_c, A_c, C_c) and fixed in the mantle frame.
- `test_moon_libration` integrates the attitude alone with DOP853, the positions from DE440; the
  rotation integrated together with the orbits is measured in
  [system of bodies](system-of-bodies.en.md).

## Data

None for the functions. The tests read `header.440`, `de440.bsp` and
`moon_pa_de440_200625.bpc` (JPL and NAIF), downloaded by `cmake/FetchData.cmake`.
