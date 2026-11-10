/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/rel.h — relativity at the first post-Newtonian order: PPN parameters, the relativistic
 * accelerations of a satellite in the local frame of its central body, the transformation
 * between the global (barycentric) and a local (body-centred) reference system, and the rate
 * of proper time.
 *
 * The N-body equations of motion of Einstein, Infeld and Hoffmann with the PPN parameters are
 * part of the N-body system (myb/dyn.h, MYB_INTERACTION_EIH).
 *
 * Sources: IERS Conventions 2010, IERS Technical Note 36, ch. 10 (eq. 10.6 - 10.12);
 * IAU 2000 Resolution B1.3 and M. Soffel et al., AJ 126, 2687 (2003); C. M. Will, Theory and
 * Experiment in Gravitational Physics, 2nd ed. (2018). The library knows no body: the central
 * body, its spin and the external body come from the caller. */
#ifndef MYB_REL_H_
#define MYB_REL_H_

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Speed of light (m/s; e.g. MYB_CONST_C of data/constants/codata-2022.tk) and the PPN
 * parameters beta and gamma (1 and 1 in general relativity). */
typedef struct myb_ppn {
  double c;
  double beta;
  double gamma;
} myb_ppn;

/* MYB_OK for c > 0 and finite c, beta, gamma; MYB_E_DOMAIN otherwise. */
MYB_API MYB_NODISCARD myb_status myb_ppn_validate(const myb_ppn* ppn);

/* ---- Satellite in the local frame of its central body (IERS 2010 eq. 10.12) ----------------
 * r, v: position (m) and velocity (m/s) of the satellite relative to the central body of
 * gravitational parameter gm (m^3 s^-2). Each function returns MYB_E_DOMAIN for r = 0 or
 * non-finite input; a receives the acceleration (m s^-2). */

/* Schwarzschild field of the central body:
 *   gm / (c^2 r^3) { [2 (beta + gamma) gm / r - gamma v.v] r + 2 (1 + gamma) (r.v) v }. */
MYB_API MYB_NODISCARD myb_status myb_rel_schwarzschild(const myb_ppn* ppn, double gm,
                                                       const double r[3], const double v[3],
                                                       double a[3]);

/* Lense-Thirring (frame dragging) by the spin of the central body, j = its angular momentum
 * per unit mass (m^2/s) in the axes of r:
 *   (1 + gamma) gm / (c^2 r^3) [ (3 / r^2) (r x v) (r.j) + v x j ]. */
MYB_API MYB_NODISCARD myb_status myb_rel_lense_thirring(const myb_ppn* ppn, double gm,
                                                        const double j[3], const double r[3],
                                                        const double v[3], double a[3]);

/* de Sitter (geodetic) precession of the local frame in the field of an external body of
 * gravitational parameter gm_external; R, V: position and velocity of the central body
 * relative to the external one:
 *   (1 + 2 gamma) [ V x (-gm_external R / (c^2 R^3)) ] x v. */
MYB_API MYB_NODISCARD myb_status myb_rel_de_sitter(const myb_ppn* ppn, double gm_external,
                                                   const double R[3], const double V[3],
                                                   const double v[3], double a[3]);

/* ---- Global (BCRS) and local (GCRS-like) coordinates (IAU 2000 Resolution B1.3) ------------
 * The local, kinematically non-rotating system of a body E whose barycentric position,
 * velocity and acceleration at the coordinate time t (TCB) are xe, ve, ae, and w_ext the
 * Newtonian potential of the other bodies at xe (sum of GM/r, m^2 s^-2). For an event (t, x),
 * with r = x - xe(t),
 *   X = r + (1/c^2) [ 1/2 ve (ve.r) + w_ext r + r (ae.r) - 1/2 ae r^2 ] + O(c^-4),
 *   T = t - (1/c^2) [ A(t) + ve.r ] + O(c^-4),   dA/dt = ve^2 / 2 + w_ext(xe),
 * T being TCG for the Earth; A(t) is TCB - TCG at the geocentre (myb/time.h). Coordinates
 * are in the TCB- and TCG-compatible units of the Resolution (rescale for TDB, TT). */
typedef struct myb_rel_local_frame {
  double c;     /* m/s */
  double xe[3]; /* m */
  double ve[3]; /* m/s */
  double ae[3]; /* m/s^2 */
  double w_ext; /* m^2/s^2 */
} myb_rel_local_frame;

/* X from x at the same t; out may alias x. */
MYB_API void myb_rel_global_to_local(const myb_rel_local_frame* frame, const double x[3],
                                     double out[3]);

/* x from X (inverse to O(c^-4) by one fixed-point step); out may alias X. */
MYB_API void myb_rel_local_to_global(const myb_rel_local_frame* frame, const double X[3],
                                     double out[3]);

/* The position-dependent part of T - t: -(ve.r) / c^2 (s), with r = x - xe. */
MYB_API double myb_rel_local_time_shift(const myb_rel_local_frame* frame, const double x[3]);

/* ---- Proper time ---------------------------------------------------------------------------
 * dtau/dt = 1 - (v^2 / 2 + U) / c^2 for a clock moving with speed v (m/s) in a coordinate
 * system where the Newtonian potential at the clock is U (m^2 s^-2, positive; IERS 2010
 * eq. 10.6 with U including the tidal terms of the external bodies, eq. 10.7 without). */
MYB_API double myb_rel_proper_time_rate(double c, double v2, double potential);

/* ---- Compact binaries beyond the first order -----------------------------------------------
 * Relative acceleration of two bodies in harmonic coordinates in the centre-of-mass frame
 * (L. Blanchet, Living Rev. Relativ. 17, 2, 2014, eq. 203-204; L. E. Kidder, Phys. Rev. D 52,
 * 821, 1995, eq. 2.2): with m = m1 + m2, nu = m1 m2 / m^2, x = x1 - x2, r = |x|, n = x / r,
 * v = v1 - v2, rdot = n.v and G m written gm,
 *   a = -(gm / r^2) [(1 + A) n + B v],
 *   A1 = -3/2 rdot^2 nu + v^2 + 3 nu v^2 - (gm/r)(4 + 2 nu),          B1 = -4 rdot + 2 rdot nu,
 *   A2 = 15/8 rdot^4 nu - 45/8 rdot^4 nu^2 - 9/2 rdot^2 nu v^2 + 6 rdot^2 nu^2 v^2 + 3 nu v^4
 *        - 4 nu^2 v^4 + (gm/r)(-2 rdot^2 - 25 rdot^2 nu - 2 rdot^2 nu^2 - 13/2 nu v^2
 *        + 2 nu^2 v^2) + (gm/r)^2 (9 + 87/4 nu),
 *   B2 = 9/2 rdot^3 nu + 3 rdot^3 nu^2 - 15/2 rdot nu v^2 - 2 rdot nu^2 v^2
 *        + (gm/r)(2 rdot + 41/2 rdot nu + 4 rdot nu^2),
 *   A5/2 = -24/5 rdot nu v^2 (gm/r) - 136/15 rdot nu (gm/r)^2,        (radiation reaction)
 *   B5/2 = 8/5 nu v^2 (gm/r) + 24/5 nu (gm/r)^2,
 * the terms of order k divided by c^(2k). In general relativity only (beta = gamma = 1). The
 * 1PN part is the two-body limit of the EIH equations (myb/dyn.h). Spins are not included. */

enum {
  MYB_PN_1 = 1,  /* 1PN, conservative */
  MYB_PN_2 = 2,  /* 2PN, conservative */
  MYB_PN_25 = 4  /* 2.5PN, radiation reaction */
};

/* Relative acceleration (m/s^2) of body 1 with respect to body 2 for gravitational parameters
 * gm1, gm2, c (m/s) and the relative position x (m) and velocity v (m/s); `orders` selects the
 * post-Newtonian terms added to the Newtonian one (0: Newton). Errors: MYB_E_DOMAIN (r = 0,
 * gm1 or gm2 <= 0, c <= 0, non-finite input). */
MYB_API MYB_NODISCARD myb_status myb_pn_binary_acceleration(double gm1, double gm2, double c,
                                                            unsigned orders, const double x[3],
                                                            const double v[3], double a[3]);

/* Peters-Mathews (Phys. Rev. 131, 435, 1963) secular change of the period of a binary by the
 * quadrupole gravitational radiation (leading order, s/s):
 *   dP/dt = -(192 pi / 5) (2 pi G m_c / (P c^3))^(5/3) f(e),
 *   f(e) = (1 + 73/24 e^2 + 37/96 e^4) / (1 - e^2)^(7/2),  G m_c = (gm1 gm2)^(3/5) / (gm1 + gm2)^(1/5).
 * NaN for invalid input. */
MYB_API double myb_peters_period_derivative(double gm1, double gm2, double c, double period,
                                            double e);

#ifdef __cplusplus
}
#endif

#endif /* MYB_REL_H_ */
