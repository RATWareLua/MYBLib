/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/regular.h — regularised equations of the perturbed two-body problem.
 *
 * Close approaches and highly eccentric orbits make the Cartesian equations stiff near the
 * pericentre (the step must follow 1/r^(3/2)). Two remedies:
 *
 * Kustaanheimo-Stiefel (P. Kustaanheimo, E. Stiefel, J. Reine Angew. Math. 218, 204, 1965;
 * E. Stiefel, G. Scheifele, Linear and Regular Celestial Mechanics, Springer 1971): x = L(u) u
 * with u in R^4 and the KS matrix
 *   L(u) = | u1 -u2 -u3  u4 |
 *          | u2  u1 -u4 -u3 |
 *          | u3  u4  u1  u2 |
 *          | u4 -u3  u2 -u1 |,
 * the fictitious time s with dt = r ds (r = |u|^2), u' = du/ds = L(u)^T v / 2, and
 *   u'' + (h / 2) u = (r / 2) L(u)^T P,   h' = -2 u'.L(u)^T P,   t' = |u|^2,
 * h = mu / r - |v|^2 / 2 the binding energy and P the perturbing acceleration (the central
 * term -mu x / r^3 excluded). The state is y = (u, u', h, t), 10 values. Unperturbed, u is a
 * harmonic oscillator of frequency sqrt(h / 2): no singularity at r = 0 and a step size
 * independent of the eccentricity.
 *
 * The time element (Stiefel & Scheifele, sec. 18): t' = |u|^2 swings between the pericentre and
 * the apocentre, and the error of its integration becomes an error along the orbit. On a bound
 * orbit (h > 0) the time is replaced by
 *   tau = t + (u . u') / h,
 *   tau' = mu / (2 h) + [(r / 2) (u . L^T P) + (2 / h) (u . u') (u' . L^T P)] / h,
 * which is constant in Kepler motion: the time is then integrated exactly, and under
 * perturbations with an error of their size. The state is y = (u, u', h, tau).
 *
 * Sundman (K. F. Sundman, Acta Math. 36, 105, 1913): Cartesian coordinates in the time s with
 * dt = r^alpha ds (alpha = 1: steps uniform in the eccentric anomaly; 2: in the true anomaly,
 * r^2 d(theta) being proportional to dt; 3/2: between the two, the intermediate anomaly of
 * P. Nacozy, Celest. Mech. 16, 309, 1977; 0: the time itself); the state is y = (x, v, t),
 * 7 values: dx/ds = r^alpha v, dv/ds = r^alpha a, dt/ds = r^alpha.
 *
 * Both right-hand sides have the signature of myb_ode_function (myb/integ.h) with a
 * myb_regular_problem as the context, for use with DOP853; myb_regular_propagate does that up to
 * a given physical time. */
#ifndef MYB_REGULAR_H_
#define MYB_REGULAR_H_

#include <stdint.h>

#include "myb/core.h"
#include "myb/integ.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Perturbing acceleration P(t, x, v) (m/s^2, frame of x), without the central term. */
typedef myb_status (*myb_perturbation_function)(void* context, double t, const double x[3],
                                                const double v[3], double a[3]);

typedef struct myb_regular_problem {
  double mu;                              /* gravitational parameter of the centre, > 0 */
  myb_perturbation_function perturbation; /* NULL: pure Kepler motion */
  void* context;
  double alpha;                           /* Sundman exponent, >= 0 (the Sundman form only) */
} myb_regular_problem;

/* (t, x, v) -> y = (u, u', h, t). Of the circle of u giving x, the one with u4 = 0 for
 * x1 >= 0, with u3 = 0 otherwise (no cancellation in either); x = 0 gives u = 0. The bilinear
 * relation u4 u1' - u3 u2' + u2 u3' - u1 u4' = 0 holds. */
MYB_API void myb_ks_from_cartesian(double mu, double t, const double x[3], const double v[3],
                                   double y[10]);
/* y -> (t, x, v); t may be NULL. */
MYB_API void myb_ks_to_cartesian(const double y[10], double* t, double x[3], double v[3]);

/* d y / ds of the KS equations; context: const myb_regular_problem*. */
MYB_API MYB_NODISCARD myb_status myb_ks_rhs(void* problem, double s, const double* y, double* dy);

/* The same with the time element in place of the time: y = (u, u', h, tau). Errors:
 * MYB_E_DOMAIN (h <= 0: the orbit is not bound). */
MYB_API MYB_NODISCARD myb_status myb_ks_element_from_cartesian(double mu, double t,
                                                               const double x[3],
                                                               const double v[3], double y[10]);
MYB_API MYB_NODISCARD myb_status myb_ks_element_to_cartesian(const double y[10], double* t,
                                                             double x[3], double v[3]);
MYB_API MYB_NODISCARD myb_status myb_ks_element_rhs(void* problem, double s, const double* y,
                                                    double* dy);

/* d y / ds of the Sundman-transformed Cartesian equations, y = (x, v, t). Any exponent
 * alpha >= 0: 1, 3/2 and 2 by multiplications and a square root, the others as
 * exp(alpha log r). Errors: MYB_E_DOMAIN (r = 0, alpha negative or not finite). */
MYB_API MYB_NODISCARD myb_status myb_sundman_rhs(void* problem, double s, const double* y,
                                                 double* dy);

typedef enum myb_regular_kind {
  MYB_REGULAR_KS = 0,
  MYB_REGULAR_SUNDMAN = 1,
  MYB_REGULAR_KS_ELEMENT = 2 /* KS with the time element: bound orbits only */
} myb_regular_kind;

/* Propagates (x0, v0) from t0 to t1 with DOP853 in the regularised time and returns the state at
 * t1 exactly (the event t(s) = t1 located on the dense output). *evaluations (may be NULL)
 * receives the count of right-hand sides. Errors: MYB_E_DOMAIN (mu <= 0, invalid input, an
 * orbit that is not bound with the time element), statuses of DOP853 and of the perturbation.
 * With the time element an orbit that a perturbation unbinds on the way ends with
 * MYB_E_DOMAIN or, the derivative of the element growing without bound as h comes to zero,
 * with MYB_E_STEP_MIN: use MYB_REGULAR_KS there. */
MYB_API MYB_NODISCARD myb_status myb_regular_propagate(myb_regular_kind kind,
                                                       const myb_regular_problem* problem,
                                                       const myb_dop853_options* options,
                                                       double t0, const double x0[3],
                                                       const double v0[3], double t1,
                                                       double x1[3], double v1[3],
                                                       int64_t* evaluations);

#ifdef __cplusplus
}
#endif

#endif /* MYB_REGULAR_H_ */
