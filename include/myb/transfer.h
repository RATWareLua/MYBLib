/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/transfer.h — transfers between bodies: hyperbolic excess and the B-plane of an approach,
 * patched-conic relations, and a general shooting (differential-correction) solver.
 * The Lambert problem is in myb/lambert.h.
 *
 * B-plane (W. Kizner, JPL EPD-457, 1959; R. Sergeyevsky, G. Snyder, R. Cunniff, JPL Publ.
 * 82-43, 1983; D. Vallado, Fundamentals of Astrodynamics, 2013, §12.6): S along the incoming
 * asymptote, T = S x k / |S x k| for a reference pole k (the ecliptic or the equator of the
 * target), R = S x T; B from the target to the asymptote, perpendicular to S. */
#ifndef MYB_TRANSFER_H_
#define MYB_TRANSFER_H_

#include <stddef.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct myb_bplane {
  double v_infinity;   /* m/s */
  double c3;           /* m^2/s^2, v_infinity^2 */
  double s[3];         /* unit vectors S, T, R */
  double t[3];
  double r[3];
  double b[3];         /* the B vector (m) */
  double b_t;          /* B.T, B.R (m) */
  double b_r;
  double b_magnitude;  /* |B| (m) */
  double theta;        /* angle from T to B, atan2(B.R, B.T) (rad) */
  double periapsis;    /* radius (m) */
  double eccentricity;
  double time_to_periapsis; /* two-body, from the given state (s; negative after it) */
} myb_bplane;

/* The B-plane of the two-body hyperbola through (r, v) relative to the target of gravitational
 * parameter gm. MYB_E_DOMAIN: not hyperbolic, or S parallel to the pole. */
MYB_API MYB_NODISCARD myb_status myb_bplane_compute(double gm, const double r[3],
                                                    const double v[3], const double pole[3],
                                                    myb_bplane* bplane);

/* Patched conics. Sphere of influence a (m/M)^(2/5) (Laplace); Hill radius a (1 - e)
 * (m / (3M))^(1/3); turn angle of a flyby 2 asin(1 / (1 + r_p v_inf^2 / gm)); the impulse
 * from a circular orbit of radius r to a hyperbola of excess v_inf,
 * sqrt(v_inf^2 + 2 gm / r) - sqrt(gm / r). */
MYB_API double myb_sphere_of_influence(double a, double mass_ratio);
MYB_API double myb_hill_radius(double a, double e, double mass_ratio);
MYB_API double myb_flyby_turn(double gm, double v_infinity, double periapsis);
MYB_API double myb_departure_impulse(double gm, double radius, double v_infinity);

/* ---- Shooting ------------------------------------------------------------------------------------
 * Newton's method for f(u) = 0 with n unknowns and m equations: the minimum-norm step for
 * m < n, least squares for m > n; the step is halved (at most 10 times) while |f| does not
 * decrease. The residual function fills f (m values) and, when jacobian is not NULL, the
 * m x n Jacobian (row-major) — from the state transition matrix of myb/uncert.h for
 * trajectory targeting. */

typedef myb_status (*myb_shoot_function)(void* context, const double* u, double* f,
                                         double* jacobian);

typedef struct myb_shoot_options {
  double tolerance;   /* on max |f_i|; > 0 */
  int max_iterations; /* 0 = 20 */
} myb_shoot_options;

/* Doubles of work memory myb_shoot needs for n unknowns and m equations. */
MYB_API size_t myb_shoot_work_size(size_t n, size_t m);

/* u (n values) is updated in place; work: myb_shoot_work_size(n, m) doubles from the caller
 * (no allocation, little stack). *iterations and *residual (max |f_i|) may be NULL.
 * Errors: MYB_E_DOMAIN (sizes outside 1 .. 64, invalid options), MYB_E_NOCONV, statuses of
 * the function. */
MYB_API MYB_NODISCARD myb_status myb_shoot(myb_shoot_function function, void* context, size_t n,
                                           size_t m, double* u,
                                           const myb_shoot_options* options, double* work,
                                           int* iterations, double* residual);

#ifdef __cplusplus
}
#endif

#endif /* MYB_TRANSFER_H_ */
