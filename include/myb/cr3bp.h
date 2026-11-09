/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/cr3bp.h — the circular restricted three-body problem.
 *
 * Units: the distance of the primaries, their total mass and their mean motion are 1; the
 * rotating frame has its origin at the barycentre, the larger primary at (-mu, 0, 0), the
 * smaller at (1 - mu, 0, 0), mu = m2 / (m1 + m2). Equations of motion
 *   x'' - 2 y' = dOmega/dx,  y'' + 2 x' = dOmega/dy,  z'' = dOmega/dz,
 *   Omega = (x^2 + y^2) / 2 + (1 - mu) / r1 + mu / r2,
 * Jacobi constant C = 2 Omega - |v|^2 (V. Szebehely, Theory of Orbits, 1967).
 *
 * Periodic orbits symmetric with respect to the xz-plane (halo, Lyapunov, distant retrograde)
 * or to the x-axis (vertical, axial) are corrected by single shooting on the half period
 * (K. C. Howell, Celest. Mech. 32, 53, 1984): the state transition matrix comes from the
 * variational equations (myb/uncert.h) with the analytic Jacobian. Third-order halo orbits of
 * D. L. Richardson (Celest. Mech. 22, 241, 1980) give first guesses. Stability from the
 * eigenvalues of the monodromy matrix; invariant manifolds from its real eigenvectors (W. S.
 * Koon, M. W. Lo, J. E. Marsden, S. D. Ross, Dynamical Systems, the Three-Body Problem and
 * Space Mission Design, 2011). */
#ifndef MYB_CR3BP_H_
#define MYB_CR3BP_H_

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Acceleration in the rotating frame; Jacobian d a / d x and d a / d v (row-major). */
MYB_API void myb_cr3bp_acceleration(double mu, const double x[3], const double v[3], double a[3]);
MYB_API void myb_cr3bp_jacobian(double mu, const double x[3], double dadx[9], double dadv[9]);
MYB_API double myb_cr3bp_jacobi(double mu, const double y[6]);

/* Terms for a force model (myb/model.h), context = a pointer to mu. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_term(void* mu, double t, const double x[3],
                                                const double v[3], double a[3]);
MYB_API MYB_NODISCARD myb_status myb_cr3bp_term_jacobian(void* mu, double t, const double x[3],
                                                         const double v[3], double dadx[9],
                                                         double dadv[9]);

/* The five Lagrange points (L1 between the primaries, L2 beyond the smaller, L3 beyond the
 * larger, L4 and L5 at +-60 degrees). MYB_E_DOMAIN: mu outside (0, 0.5]. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_lagrange(double mu, double points[5][3]);

typedef enum myb_cr3bp_symmetry {
  MYB_CR3BP_SYMMETRY_XZ = 1, /* y0 = (x, 0, z, 0, vy, 0); at T/2: y = 0, vx = vz = 0 */
  MYB_CR3BP_SYMMETRY_X = 2   /* y0 = (x, 0, 0, 0, vy, vz); at T/2: y = 0, z = 0, vx = 0 */
} myb_cr3bp_symmetry;

typedef struct myb_cr3bp_correction {
  myb_cr3bp_symmetry symmetry;
  int fixed;          /* the component of y0 held fixed: XZ: 0 (x), 2 (z) or 4 (vy);
                         X: 0 (x), 4 (vy) or 5 (vz) */
  int crossing;       /* the half period ends at this crossing of y = 0 (1 = the first) */
  double tolerance;   /* on the constraints; 0 = 1e-12 */
  int max_iterations; /* 0 = 30 */
} myb_cr3bp_correction;

/* Corrects y0 (in place) to a periodic orbit of the given symmetry and returns its period.
 * Errors: MYB_E_DOMAIN (y0 not of the symmetric form, invalid settings), MYB_E_NOCONV (no
 * convergence, or no crossing within 50 time units). */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_correct(double mu,
                                                   const myb_cr3bp_correction* correction,
                                                   double y0[6], double* period);

/* Monodromy matrix (state transition matrix over one period) and its eigenvalues; index =
 * (|lambda_max| + 1 / |lambda_max|) / 2 (1 for a linearly stable orbit), the stability index
 * of the JPL catalogue. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_monodromy(double mu, const double y0[6], double period,
                                                     double monodromy[36]);
MYB_API MYB_NODISCARD myb_status myb_cr3bp_stability(const double monodromy[36], double re[6],
                                                     double im[6], double* index);

/* A state on the unstable (unstable != 0) or stable manifold of the periodic orbit at the
 * fraction phase (0 .. 1) of the period: the orbit's state displaced by epsilon (in length
 * units, along the position part of the eigenvector of the monodromy carried to that phase)
 * on the side sign (+1 or -1). MYB_E_DOMAIN: the orbit is not unstable. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_manifold(double mu, const double y0[6], double period,
                                                    double phase, int unstable, int sign,
                                                    double epsilon, double y[6]);

/* Third-order halo orbit of Richardson about L1 (point = 1) or L2 (2) with out-of-plane
 * amplitude az (in units of the distance of the primaries), branch +1 (northern, z > 0 at
 * the start) or -1: y0 at the xz-plane crossing with vx = vz = 0 (approximately), and its
 * period. MYB_E_DOMAIN: invalid input or an amplitude without a halo orbit. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_richardson(double mu, int point, double az, int branch,
                                                      double y0[6], double* period);

/* ---- Frames ------------------------------------------------------------------------------------
 * The rotating frame of two real bodies: origin at their barycentre, x along body 1 -> body 2,
 * z along their orbital angular momentum, lengths in units of their instantaneous distance
 * and time in units of 1 / n, n = sqrt((gm1 + gm2) / distance^3) (the "pulsating" frame of
 * ephemeris-based libration-point studies; the rotation of the orbital plane is neglected in
 * the velocities). States in m, m/s relative to any common origin. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_from_ephemeris(const double body1[6],
                                                          const double body2[6], double gm1,
                                                          double gm2, const double state[6],
                                                          double rotating[6]);
MYB_API MYB_NODISCARD myb_status myb_cr3bp_to_ephemeris(const double body1[6],
                                                        const double body2[6], double gm1,
                                                        double gm2, const double rotating[6],
                                                        double state[6]);

#ifdef __cplusplus
}
#endif

#endif /* MYB_CR3BP_H_ */
