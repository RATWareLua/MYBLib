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
 * (L. Blanchet, "Post-Newtonian theory for gravitational waves", Living Rev. Relativ., revision
 * of 2024, eq. calAcoeff and calBcoeff = eq. 203-204 of the 2014 edition; to 2.5PN also
 * L. E. Kidder, Phys. Rev. D 52, 821, 1995, eq. 2.2): with m = m1 + m2, nu = m1 m2 / m^2,
 * x = x1 - x2, r = |x|, n = x / r, v = v1 - v2, rdot = n.v, u = G m / r,
 *   a = -(G m / r^2) [(1 + A) n + B v],
 *   A1 = -3/2 rdot^2 nu + v^2 + 3 nu v^2 - u (4 + 2 nu),          B1 = -4 rdot + 2 rdot nu,
 *   A2 = 15/8 rdot^4 nu - 45/8 rdot^4 nu^2 - 9/2 rdot^2 nu v^2 + 6 rdot^2 nu^2 v^2 + 3 nu v^4
 *        - 4 nu^2 v^4 + u (-2 rdot^2 - 25 rdot^2 nu - 2 rdot^2 nu^2 - 13/2 nu v^2
 *        + 2 nu^2 v^2) + u^2 (9 + 87/4 nu),
 *   B2 = 9/2 rdot^3 nu + 3 rdot^3 nu^2 - 15/2 rdot nu v^2 - 2 rdot nu^2 v^2
 *        + u (2 rdot + 41/2 rdot nu + 4 rdot nu^2),
 *   A5/2 = -24/5 rdot nu v^2 u - 136/15 rdot nu u^2,                 (radiation reaction)
 *   B5/2 = 8/5 nu v^2 u + 24/5 nu u^2,
 *   A3, B3: the 3PN terms of the reference (cubic in nu, with pi^2 and ln(r / r0'), where r0' is
 *        a gauge length of the harmonic coordinates: it changes the coordinates, not the
 *        physics — gauge-invariant quantities such as Omega_phi / Omega_r do not depend on it),
 *   A7/2, B7/2: radiation reaction at the first relative order (Iyer & Will 1993-1995 in the
 *        harmonic gauge),
 * the terms of order k divided by c^(2k). Not included: the 4PN terms (with the hereditary tail
 * integral) and beyond. In general relativity only (beta = gamma = 1). The 1PN part is the
 * two-body limit of the EIH equations (myb/dyn.h).
 *
 * Spins (Kidder 1995, eq. 2.2c, 2.2e, 2.4; covariant spin supplementary condition) enter as
 * dimensionless vectors chi_A = c S_A / (G m_A^2) (|chi| <= 1 for black holes):
 *   spin-orbit (1.5PN for |chi| ~ 1), with S = S1 + S2, Delta = m (S2 / m2 - S1 / m1),
 *   dm = m1 - m2:  a_SO = (G / (c^2 r^3)) { 6 n [(n x v).(2 S + (dm/m) Delta)]
 *                        - v x (7 S + 3 (dm/m) Delta) + 3 rdot n x (3 S + (dm/m) Delta) },
 *   spin-spin (2PN):  a_SS = -(3 G / (mu c^2 r^4)) { n (S1.S2) + S1 (n.S2) + S2 (n.S1)
 *                        - 5 n (n.S1)(n.S2) },
 *   spin-induced quadrupoles (E. Poisson, Phys. Rev. D 57, 5287, 1998; Q_A = -kappa_A S_A^2 /
 *   (m_A c^2), kappa = 1 for a Kerr black hole, 2-10 for neutron stars):
 *     a_QM = -sum_A (3 kappa_A G m / (2 m_A^2 c^2 r^4)) [2 (n.S_A) S_A + (S_A^2 - 5 (n.S_A)^2) n],
 * and the spins precess, dS_A/dt = Omega_A x S_A, with
 *   Omega_1 = (G / (c^2 r^3)) { (2 + 3 m2 / (2 m1)) mu x x v - S2 + 3 (n.S2) n
 *             + 3 kappa_1 (m2 / m1) (n.S1) n }   (and 1 <-> 2),
 * so |chi_A| is constant. The first term is the geodetic precession (Barker & O'Connell 1975):
 * averaged over an orbit, (2 + 3 m2 / (2 m1)) G mu n_orbit / (c^2 a (1 - e^2)). Not included:
 * the spin-orbit terms of higher order (2.5PN and beyond), radiation reaction by the spins. */

enum {
  MYB_PN_1 = 1,    /* 1PN, conservative */
  MYB_PN_2 = 2,    /* 2PN, conservative */
  MYB_PN_25 = 4,   /* 2.5PN, radiation reaction */
  MYB_PN_3 = 8,    /* 3PN, conservative (needs the gauge length r0) */
  MYB_PN_35 = 16,  /* 3.5PN, radiation reaction */
  MYB_PN_SO = 32,  /* spin-orbit coupling and geodetic precession of the spins */
  MYB_PN_SS = 64,  /* spin-spin coupling and spin-induced quadrupoles, with their precession */
  MYB_PN_ALL = 127 /* every term above */
};

typedef struct myb_pn_binary {
  double gm1;      /* G m1 (m^3/s^2) */
  double gm2;      /* G m2 */
  double c;        /* speed of light (m/s) */
  unsigned orders; /* MYB_PN_* terms added to the Newtonian acceleration */
  double r0;       /* gauge length r0' of the 3PN logarithms (m, > 0; needed with MYB_PN_3) */
  double kappa1;   /* spin-induced quadrupole constants (>= 0; 0 omits the quadrupoles) */
  double kappa2;
} myb_pn_binary;

/* Relative acceleration a (m/s^2) of body 1 with respect to body 2 for the relative position x
 * (m) and velocity v (m/s), and the rates (1/s) of the dimensionless spins chi1, chi2. chi1,
 * chi2 may be NULL (no spin); dchi1, dchi2 may be NULL (not wanted). Errors: MYB_E_DOMAIN
 * (r = 0, gm1 or gm2 <= 0, c <= 0, non-finite input, unknown bits in orders, MYB_PN_3 without
 * r0 > 0, MYB_PN_SS with a negative kappa). */
MYB_API MYB_NODISCARD myb_status myb_pn_binary_evaluate(const myb_pn_binary* binary,
                                                        const double x[3], const double v[3],
                                                        const double chi1[3],
                                                        const double chi2[3], double a[3],
                                                        double dchi1[3], double dchi2[3]);

/* The same without spins and without the 3PN terms (which need r0): myb_pn_binary_evaluate
 * with r0 unset. Errors: those of myb_pn_binary_evaluate. */
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
