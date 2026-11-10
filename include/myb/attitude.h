/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/attitude.h — orientation and rotational dynamics of rigid bodies.
 *
 * Conventions. An attitude is the frame rotation M from an inertial frame to the body frame
 * (v_body = M v_inertial; the rows of M are the body axes in inertial coordinates), as
 * myb_rotation.m elsewhere in the library. Its quaternion q = (w, x, y, z) is scalar first,
 * Hamilton product, with M(q) v = conj(q) * v * q: a frame rotated by the angle a about the
 * unit axis n has q = (cos(a/2), n sin(a/2)) and M = R_n(a) (the elementary frame rotations
 * R1, R2, R3 of myb/frames.h). Successive rotations about the body's own axes compose as
 * q = q1 * q2 (first q1, then q2): M(q1 * q2) = M(q2) M(q1). The body angular velocity omega
 * (body coordinates) gives dM/dt = -[omega x] M and dq/dt = q * (0, omega) / 2.
 *
 * Rotational dynamics in principal axes (inertia diag(I1, I2, I3)) with the body angular
 * momentum L = I omega as the state: dL/dt = L x omega + torque. Inertia and torque may be
 * absolute (kg m^2, N m) or per unit mass of the body (m^2, m^2/s^2): the motion depends only
 * on their ratio, and the gravitational torques below come per unit mass. */
#ifndef MYB_ATTITUDE_H_
#define MYB_ATTITUDE_H_

#include <stdint.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Quaternions ------------------------------------------------------------------------- */

MYB_API void myb_quat_multiply(const double a[4], const double b[4], double out[4]);
MYB_API void myb_quat_conjugate(const double q[4], double out[4]);
/* Scales q to unit norm and returns the former norm (q unchanged if it is 0). */
MYB_API double myb_quat_normalize(double q[4]);
/* M(q) of a nonzero q (divided by |q|^2, so exact for any scale). */
MYB_API void myb_quat_to_matrix(const double q[4], double m[9]);
/* The unit quaternion of a rotation matrix, w >= 0 (Shepperd's method: the largest of
 * w, x, y, z from the diagonal, the others from the off-diagonal sums). */
MYB_API void myb_quat_from_matrix(const double m[9], double q[4]);
/* out = M(q) v (out may alias v). */
MYB_API void myb_quat_rotate(const double q[4], const double v[3], double out[3]);
/* Rotation vector phi = angle * axis <-> quaternion (exponential and logarithm, accurate for
 * small angles; the logarithm returns the angle in [0, pi]). */
MYB_API void myb_quat_from_rotation_vector(const double phi[3], double q[4]);
MYB_API void myb_quat_to_rotation_vector(const double q[4], double phi[3]);
/* Spherical linear interpolation along the shorter arc, t in [0, 1]. */
MYB_API void myb_quat_slerp(const double a[4], const double b[4], double t, double out[4]);
/* dq/dt = q * (0, omega) / 2 for the body angular velocity omega. */
MYB_API void myb_quat_derivative(const double q[4], const double omega[3], double dq[4]);

/* ---- Euler angles ---------------------------------------------------------------------------
 * axes = (i, j, k), each 1, 2 or 3, i != j, j != k: the twelve sequences, proper (i == k, e.g.
 * 3-1-3) and Tait-Bryan (all different, e.g. 3-2-1). M = R_k(c) R_j(b) R_i(a) for the angles
 * (a, b, c): first a about axis i, then b about the new j, then c about the new k. From a
 * matrix: b in [0, pi] (proper) or [-pi/2, pi/2] (Tait-Bryan), a and c in (-pi, pi]; at the
 * singular b (gimbal lock) a = 0. Errors: MYB_E_DOMAIN (an invalid sequence). */
MYB_API MYB_NODISCARD myb_status myb_euler_to_matrix(const int axes[3], const double angles[3],
                                                     double m[9]);
MYB_API MYB_NODISCARD myb_status myb_euler_from_matrix(const int axes[3], const double m[9],
                                                       double angles[3]);

/* ---- Inertia ------------------------------------------------------------------------------- */

/* Principal moments (ascending) and axes of a symmetric inertia tensor: axes holds the rotation
 * from the given frame to the principal frame (rows = principal axes, det +1), so that
 * axes inertia axes^T = diag(moments). Cyclic Jacobi rotations to below roundoff.
 * Errors: MYB_E_DOMAIN (not symmetric, not finite). */
MYB_API MYB_NODISCARD myb_status myb_principal_axes(const double inertia[9], double moments[3],
                                                    double axes[9]);

/* ---- Equations of motion --------------------------------------------------------------------- */

/* dq/dt and dL/dt for principal moments, attitude q (unit), body angular momentum l and body
 * torque (may be NULL: free). */
MYB_API void myb_rigid_body_rates(const double moments[3], const double q[4], const double l[3],
                                  const double torque[3], double dq[4], double dl[3]);

/* Rotation of a body with a fluid core (H. Poincare 1910; the model of the lunar ephemerides,
 * W. M. Folkner et al., IPN Progress Report 42-196, 2014, §IV): a mantle of inertia tensor I_m
 * (body frame, may vary with tides and spin: dI_m its rate, NULL = constant) and a core of
 * principal moments (A_c, A_c, C_c) fixed in the mantle frame, rotating with omega_c (given in
 * the mantle frame), coupled by the torque on the mantle
 *   N_cmb = k_v (omega_c - omega_m) + (C_c - A_c) (z . omega_c) (z x omega_c):
 *   I_m omega_m' = N + N_cmb - dI_m omega_m - omega_m x I_m omega_m,
 *   I_c omega_c' = -omega_m x I_c omega_c - N_cmb   (components in the rotating mantle frame).
 * N is the external torque on the mantle. Units as above (per unit mass of the whole body, or
 * absolute). Errors: MYB_E_DOMAIN (I_m singular, core moments not positive). */
MYB_API MYB_NODISCARD myb_status myb_core_mantle_rates(const double inertia_mantle[9],
                                                       const double inertia_rate[9],
                                                       const double core[2], double kv,
                                                       const double omega_mantle[3],
                                                       const double omega_core[3],
                                                       const double torque[3],
                                                       double rate_mantle[3],
                                                       double rate_core[3]);

/* ---- Gravitational torques (body frame) -------------------------------------------------------
 * r: position of the attracting body relative to the centre of mass of the torqued body, in its
 * body frame (m). */

/* Point mass gm on a body of inertia diag(moments): 3 gm / r^5 (r x I r) (the quadrupole term;
 * same units as the inertia times s^-2). */
MYB_API void myb_torque_point_mass(double gm, const double r[3], const double moments[3],
                                   double torque[3]);

/* Point mass of mass ratio m_other / m_body on a body whose own gravity field gives the
 * acceleration a (m/s^2) at r: torque per unit mass of the body, -(m_other / m_body) r x a.
 * Exact for the whole field (e.g. myb_gravity_field_eval in the body frame). */
MYB_API void myb_torque_field(double mass_ratio, const double r[3], const double a[3],
                              double torque[3]);

/* Torque per unit mass of A from the quadrupole-quadrupole potential U22 of two bodies
 * (myb/grav.h, myb_figure_figure; G Q row-major, both in the frame of r = x_B - x_A):
 * -2 axial(dU/dQ_A Q_A - Q_A dU/dQ_A) / gm_a, so that A's spin and the orbit exchange angular
 * momentum exactly. Errors: MYB_E_DOMAIN (r = 0, gm_a <= 0, non-finite input). */
MYB_API MYB_NODISCARD myb_status myb_torque_figure_figure(double gm_a, const double gq_a[9],
                                                          const double gq_b[9],
                                                          const double r[3], double torque[3]);

/* ---- Torque-free rotation, exact ---------------------------------------------------------------
 * The Euler-Poinsot motion in closed form (L. D. Landau, E. M. Lifshitz, Mechanics, §37;
 * E. Celledoni, F. Fasso, N. Safstrom, A. Zanna, SIAM J. Sci. Comput. 30, 2084, 2008): with the
 * axes ordered so that the body circulates around axis 3 (I2 the middle moment), the body
 * angular momentum is (a1 cn u, s a2 sn u, s3 a3 dn u), u = u0 + nu t, Jacobi functions of
 * parameter k^2, and the attitude is R3(psi) R1(theta) R3(phi) Z with cos theta = L3 / |L|,
 * tan psi = L1 / L2 and the precession about the fixed momentum
 *   phi = |L| t / I3 + |L| (I3 - I1) / (I1 I3 nu) [Pi(n; am u, k^2) - Pi(n; am u0, k^2)],
 *   n = -I3 (I2 - I1) / (I1 (I3 - I2)) <= 0.
 * Symmetric bodies (k = 0), steady rotations and spheres are included; the separatrix
 * (|L|^2 = 2 T I2 with the motion not steady) is not. Errors grow like the Jacobi functions'
 * near it. */

typedef struct myb_free_rotation {
  int kind;          /* internal: 0 general, 1 steady */
  double perm[9];    /* principal -> ordered axes (signed permutation, det +1) */
  double inertia[3]; /* ordered moments */
  double g;          /* |L| */
  double z[9];       /* inertial -> frame of L (L along +z), in ordered axes */
  double a[3];       /* amplitudes */
  double sign2, sign3;
  double nu, k2, n, u0, pi0, rate0, rate_pi;
  double omega[3];   /* steady: constant body angular velocity (ordered axes) */
  double m0[9];      /* steady: initial attitude (ordered axes) */
} myb_free_rotation;

/* From principal moments (> 0), the attitude q and the body angular momentum l (principal axes)
 * at t = 0. Errors: MYB_E_DOMAIN (moments not positive or not finite, l = 0 is allowed, the
 * separatrix). */
MYB_API MYB_NODISCARD myb_status myb_free_rotation_init(const double moments[3], const double q[4],
                                                        const double l[3],
                                                        myb_free_rotation* rotation);

/* Attitude and body angular momentum at t (either may be NULL). */
MYB_API void myb_free_rotation_state(const myb_free_rotation* rotation, double t, double q[4],
                                     double l[3]);

/* ---- Splitting integrator ---------------------------------------------------------------------
 * For torques from a potential (depending on the time and the attitude, not on L): the
 * symmetric composition of kick(h/2) - exact free rotation(h) - kick(h/2) (J. Touma,
 * J. Wisdom, Astron. J. 107, 1189, 1994; Celledoni et al. 2008), with Yoshida's weights for
 * orders 2, 4, 6, 8 (myb/integ.h). Symplectic and time-reversible: |L| and the energy stay
 * bounded over any number of steps; with zero torque the result is exact. */

typedef myb_status (*myb_torque_function)(void* context, double t, const double q[4],
                                          double torque[3]);

/* Advances (q, l) from t0 by `steps` steps of h. Errors: MYB_E_DOMAIN (order not 2, 4, 6, 8,
 * steps < 0, invalid moments, a free rotation on the separatrix), statuses of the torque. */
MYB_API MYB_NODISCARD myb_status myb_rotation_split(const double moments[3],
                                                    myb_torque_function torque, void* context,
                                                    int order, double t0, double h,
                                                    int64_t steps, double q[4], double l[3]);

#ifdef __cplusplus
}
#endif

#endif /* MYB_ATTITUDE_H_ */
