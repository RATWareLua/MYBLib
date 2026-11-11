/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/kepler.h — the two-body problem: Kepler's equation, anomalies, orbital elements and
 * propagation of a state by universal variables. Pure mathematics of one central body with
 * gravitational parameter mu (m^3 s^-2); vectors are in any inertial frame (m, m/s).
 *
 * Sources: R. H. Battin, An Introduction to the Mathematics and Methods of Astrodynamics
 * (1999); D. A. Vallado, Fundamentals of Astrodynamics and Applications, 5th ed. (2022);
 * O. Montenbruck, E. Gill, Satellite Orbits (2000); M. J. H. Walker, B. Ireland, J. Owens,
 * "A set of modified equinoctial orbit elements", Celest. Mech. 36 (1985); R. A. Broucke,
 * P. J. Cefola, "On the equinoctial orbit elements", Celest. Mech. 5 (1972). */
#ifndef MYB_KEPLER_H_
#define MYB_KEPLER_H_

#include <stddef.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Position (m) and velocity (m/s) relative to the central body. */
typedef struct myb_pv {
  double r[3];
  double v[3];
} myb_pv;

/* ---- Kepler's equation ---------------------------------------------------------------------
 * The residuals are written without cancellation, (1 - e) E + e (E - sin E) - M with the
 * Stumpff series for E - sin E, so that near-parabolic orbits keep full accuracy. Solutions
 * use Newton's method kept inside a guaranteed bracket; MYB_E_NOCONV is never expected. */

/* E - e sin E = M for 0 <= e < 1 and any finite M; E is in the same revolution as M.
 * MYB_E_DOMAIN for e outside [0, 1) or non-finite input. */
MYB_API MYB_NODISCARD myb_status myb_kepler_elliptic(double mean_anomaly, double e,
                                                     double* eccentric_anomaly);

/* e sinh H - H = M for e > 1 and any finite M. MYB_E_DOMAIN otherwise. */
MYB_API MYB_NODISCARD myb_status myb_kepler_hyperbolic(double mean_anomaly, double e,
                                                       double* hyperbolic_anomaly);

/* Barker's equation D + D^3 / 3 = M, D = tan(nu / 2), M = 2 sqrt(mu / p^3) (t - T_periapsis):
 * closed form D = 2 sinh(asinh(3M / 2) / 3). */
MYB_API double myb_kepler_parabolic(double mean_anomaly);

/* ---- Anomalies -----------------------------------------------------------------------------
 * Angles in radians; true anomalies are returned in (-pi, pi]. MYB_E_DOMAIN for an
 * eccentricity outside the conic's range, or a true anomaly beyond the asymptotes of a
 * hyperbola. */
MYB_API MYB_NODISCARD myb_status myb_anomaly_true_from_eccentric(double eccentric_anomaly, double e,
                                                                 double* nu);
MYB_API MYB_NODISCARD myb_status myb_anomaly_eccentric_from_true(double nu, double e,
                                                                 double* eccentric_anomaly);
MYB_API MYB_NODISCARD myb_status myb_anomaly_mean_from_eccentric(double eccentric_anomaly, double e,
                                                                 double* mean_anomaly);
MYB_API MYB_NODISCARD myb_status myb_anomaly_true_from_hyperbolic(double hyperbolic_anomaly,
                                                                  double e, double* nu);
MYB_API MYB_NODISCARD myb_status myb_anomaly_hyperbolic_from_true(double nu, double e,
                                                                  double* hyperbolic_anomaly);
MYB_API MYB_NODISCARD myb_status myb_anomaly_mean_from_hyperbolic(double hyperbolic_anomaly,
                                                                  double e, double* mean_anomaly);

/* ---- Orbital elements ----------------------------------------------------------------------- */

/* Classical elements. The semi-latus rectum p is used instead of the semi-major axis so that
 * the parabola (e = 1) is included.
 * Conventions for the undefined angles (exactly circular or exactly equatorial orbits):
 *   e == 0: argp = 0 and nu is the argument of latitude (from the node);
 *   h parallel to z (i = 0 or pi): raan = 0 and argp is measured from the x axis;
 *   both: nu is the true longitude.
 * Near these cases the classical angles are ill-conditioned (only their sums are well
 * defined); use the modified equinoctial elements there. */
typedef struct myb_elements {
  double p;    /* semi-latus rectum h^2 / mu, m */
  double e;    /* eccentricity */
  double i;    /* inclination, [0, pi] */
  double raan; /* right ascension of the ascending node, [0, 2 pi) */
  double argp; /* argument of periapsis, [0, 2 pi) */
  double nu;   /* true anomaly, (-pi, pi] */
} myb_elements;

/* Modified equinoctial elements (Walker, Ireland, Owens 1985), non-singular for every orbit
 * except i = pi exactly:
 *   f = e cos(argp + raan), g = e sin(argp + raan),
 *   h = tan(i/2) cos(raan), k = tan(i/2) sin(raan), L = raan + argp + nu (true longitude). */
typedef struct myb_equinoctial {
  double p;
  double f;
  double g;
  double h;
  double k;
  double l; /* (-pi, pi] */
} myb_equinoctial;

/* Errors: MYB_E_DOMAIN for mu <= 0, non-finite input, r = 0, or zero angular momentum
 * (rectilinear motion has no orbital plane). */
MYB_API MYB_NODISCARD myb_status myb_elements_from_pv(double mu, const myb_pv* pv,
                                                      myb_elements* elements);
/* MYB_E_DOMAIN for p <= 0, e < 0, or nu beyond the asymptotes (1 + e cos nu <= 0). */
MYB_API MYB_NODISCARD myb_status myb_pv_from_elements(double mu, const myb_elements* elements,
                                                      myb_pv* pv);

/* As above; additionally MYB_E_DOMAIN for i = pi (retrograde equatorial). */
MYB_API MYB_NODISCARD myb_status myb_equinoctial_from_pv(double mu, const myb_pv* pv,
                                                         myb_equinoctial* elements);
MYB_API MYB_NODISCARD myb_status myb_pv_from_equinoctial(double mu, const myb_equinoctial* elements,
                                                         myb_pv* pv);

/* Mean elements by averaging: the elements of a perturbed orbit oscillate about values that
 * change slowly; the means of the semi-major axis and of the equinoctial f, g, h, k over
 * samples that span whole revolutions are free of the oscillation to the first order, and
 * the mean longitude M + argp + raan is fitted by a line — its value at the mean of the times
 * and its rate, the mean motion with the secular drifts in it. Nothing is assumed about the
 * forces. (The mean elements of an analytic theory are another thing: for those of SGP4 see
 * myb_tle_fit of ext/sgp4.) */
typedef struct myb_mean_elements {
  double a;      /* m */
  double f;      /* e cos(argp + raan) */
  double g;      /* e sin(argp + raan) */
  double h;      /* tan(i/2) cos(raan) */
  double k;      /* tan(i/2) sin(raan) */
  double lambda; /* mean longitude at `epoch`, [0, 2 pi) */
  double rate;   /* of the mean longitude (rad/s) */
  double epoch;  /* the mean of the times (s) */
} myb_mean_elements;

/* From count >= 3 states (6 values each: m, m/s) at increasing times (s), neighbours less
 * than a revolution apart, the whole spanning whole revolutions as closely as the caller can
 * arrange (a part of a revolution left over biases the means by that part of the
 * oscillation). Errors: MYB_E_DOMAIN (mu <= 0, too few states, times not increasing, a state
 * that is no ellipse), those of myb_equinoctial_from_pv. */
MYB_API MYB_NODISCARD myb_status myb_mean_elements_average(double mu, const double* times,
                                                           const double* states, size_t count,
                                                           myb_mean_elements* mean);

/* ---- Propagation ------------------------------------------------------------------------- */

/* State after dt seconds (either sign) on the two-body orbit through pv0, for every conic
 * (ellipse, parabola, hyperbola, and rectilinear motion while r > 0), by the universal
 * Kepler equation in the variable chi with the universal functions U0..U3:
 *   sqrt(mu) dt = r0 U1 + sigma0 U2 + U3,   sigma0 = r0.v0 / sqrt(mu),
 *   f = 1 - U2/r0,  g = (r0 U1 + sigma0 U2) / sqrt(mu)   (no cancellation in g),
 *   fdot = -sqrt(mu) U1 / (r r0),  gdot = 1 - U2 / r.
 * For ellipses dt is first reduced modulo the period; alpha = 2/r0 - v0^2/mu and the period
 * are computed in double-word arithmetic, so that no phase error accumulates with the periods
 * removed. A hyperbola started far from its periapsis (|H0| > 1) is restarted from the
 * periapsis, where the equation has no cancellation (measured: relative error <= 5e-15
 * against a 400-bit MPFR reference for ellipses, also over 1 to 30 periods, hyperbolas,
 * near-parabolic and near-circular orbits).
 * Errors: MYB_E_DOMAIN (mu <= 0, r0 = 0, non-finite input, or a state that overflows),
 * MYB_E_NOCONV (iteration budget exhausted; not expected). */
MYB_API MYB_NODISCARD myb_status myb_kepler_propagate(double mu, const myb_pv* pv0, double dt,
                                                      myb_pv* pv);

/* State after dt (as myb_kepler_propagate, same routes) and the state transition matrix
 * stm[6 i + j] = d y_i / d y0_j, y = (r, v), in closed form: the Lagrange coefficients
 * differentiated through the universal Kepler equation (cf. S. W. Shepperd, Celest. Mech. 35,
 * 129, 1985; R. H. Battin 1999, §9.7). The whole periods removed for ellipses contribute their
 * dependence on the energy explicitly (the secular term); a hyperbola far from its periapsis
 * is composed of two legs from the periapsis, Phi = Phi(p -> t) Phi(p -> t0)^-1, the inverse
 * of a symplectic matrix taken exactly. Errors as myb_kepler_propagate. */
MYB_API MYB_NODISCARD myb_status myb_kepler_stm(double mu, const myb_pv* pv0, double dt,
                                                myb_pv* pv, double stm[36]);

#ifdef __cplusplus
}
#endif

#endif /* MYB_KEPLER_H_ */
