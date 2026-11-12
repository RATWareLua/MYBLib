/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/integ.h — numerical integration of ordinary differential equations.
 *
 * The integrators know nothing about physics: they advance y' = f(t, y) for a state vector of
 * any dimension n given by a callback. Time t is a double relative to an epoch chosen by the
 * caller (e.g. a myb_time); keep |t| moderate for full resolution.
 *
 * Reproducibility: every operation has a fixed order; the state and the time are updated with
 * compensated summation (the rounding error of each update is carried to the next one), so the
 * accumulated round-off grows like a random walk and not linearly with the number of steps.
 * Memory is allocated once, when the integrator is created; stepping never allocates. */
#ifndef MYB_INTEG_H_
#define MYB_INTEG_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Right-hand side: dydt = f(t, y) for n components. Any status other than MYB_OK stops the
 * integration and is returned to the caller unchanged. */
typedef myb_status (*myb_ode_function)(void* context, double t, const double* y, double* dydt);

/* Scalar event function g(t, y); events are its sign changes. */
typedef double (*myb_event_function)(void* context, double t, const double* y);

/* ---- DOP853 -------------------------------------------------------------------------------
 * Explicit Runge-Kutta 8(5,3) of Dormand and Prince with Hairer's combined error estimate and a
 * dense output of order 7 (E. Hairer, S. P. Norsett, G. Wanner, Solving ODE I, 2nd ed., 1993,
 * §II.5, §II.6). A general-purpose high-accuracy method for smooth first-order systems
 * (e.g. orbits with their variational equations). Step-size control: Hairer's controller with
 * factor 0.9 err^(-1/8) bounded to [0.2, 10]; error norm
 *   err = sqrt(mean_i ((e_i / s_i)^2)),   s_i = atol + rtol max(|y_i|, |y_new_i|). */

typedef struct myb_dop853 myb_dop853;

typedef struct myb_dop853_options {
  double rtol;      /* relative tolerance, > 0 */
  double atol;      /* absolute tolerance, > 0 */
  double h_initial; /* first step size; 0 = automatic (Hairer's starting algorithm) */
  double h_max;     /* largest step magnitude; 0 = unlimited */
  double h_min;     /* smallest step magnitude before MYB_E_STEP_MIN; 0 = none */
  int64_t max_steps;   /* budget of accepted steps per call of myb_dop853_integrate; 0 = 10^9 */
} myb_dop853_options;

typedef struct myb_dop853_stats {
  int64_t accepted;
  int64_t rejected;
  int64_t evaluations; /* calls of the right-hand side */
} myb_dop853_stats;

/* Creates an integrator for n components. Errors: MYB_E_DOMAIN (invalid options, n == 0),
 * MYB_E_NOMEM, MYB_E_OVERFLOW. allocator may be NULL. */
MYB_API MYB_NODISCARD myb_status myb_dop853_create(size_t n, const myb_dop853_options* options,
                                                   const myb_allocator* allocator,
                                                   myb_dop853** integrator);
MYB_API void myb_dop853_destroy(myb_dop853* integrator); /* NULL is allowed */

/* Sets the problem and the initial state (copied). Evaluates f(t0, y0) once. */
MYB_API MYB_NODISCARD myb_status myb_dop853_start(myb_dop853* integrator, myb_ode_function f,
                                                  void* context, double t0, const double* y0);

/* One accepted step towards t_limit (never beyond it; the direction of integration is the sign
 * of t_limit - t). Rejected attempts are retried inside. Errors: MYB_E_STEP_MIN, MYB_E_NOCONV
 * (budget), a status from f. */
MYB_API MYB_NODISCARD myb_status myb_dop853_step(myb_dop853* integrator, double t_limit);

/* Steps until t_end exactly. */
MYB_API MYB_NODISCARD myb_status myb_dop853_integrate(myb_dop853* integrator, double t_end);

MYB_API double myb_dop853_time(const myb_dop853* integrator);
MYB_API const double* myb_dop853_state(const myb_dop853* integrator); /* n values */
MYB_API myb_dop853_stats myb_dop853_statistics(const myb_dop853* integrator);

/* Dense output: the state at any t within the last accepted step (order 7). The first call
 * after a step evaluates f three more times. MYB_E_RANGE if t is outside the last step or no
 * step was taken. */
MYB_API MYB_NODISCARD myb_status myb_dop853_dense(myb_dop853* integrator, double t, double* y);

/* Root of g(t, y(t)) inside the last accepted step, if g changes sign between its ends, located
 * on the dense output by the Illinois variant of regula falsi to the resolution of t.
 * *found = 0 without a sign change. y_event may be NULL. */
MYB_API MYB_NODISCARD myb_status myb_dop853_locate(myb_dop853* integrator, myb_event_function g,
                                                   void* context, int* found, double* t_event,
                                                   double* y_event);

/* ---- IAS15 --------------------------------------------------------------------------------
 * Implicit integrator of order 15 with adaptive steps for second-order systems
 * x'' = a(t, x, v), v = x' (H. Rein, D. S. Spiegel, MNRAS 446, 1424, 2015): Gauss-Radau
 * quadrature on 8 nodes in Everhart's formulation, the implicit equations solved by
 * predictor-corrector iteration until the correction reaches round-off. Step size from the
 * timescale of D. Pham, H. Rein, D. S. Spiegel (Open J. Astrophys. 7, 2024, eq. 15-17):
 *   tau^2 = 2 |a|^2 / (|a'|^2 + |a| |a''|),   dt = (7! epsilon)^(1/7) min tau,
 * a, a', a'' at the end of the step from the step's own polynomial, each measured by the norm
 * of eq. 7: the maximum over blocks of components (e.g. 3: the Cartesian vector of one body) of
 * the Euclidean norm of the block, so a body whose acceleration passes through zero (the middle
 * body of a symmetric collinear configuration) does not stall the step; the next step is at
 * most 4 times the last one, and a step whose successor would be below a quarter of it is
 * rejected and redone.
 * With the default epsilon = 1e-9 the truncation error stays below round-off, so the energy
 * error of a conservative system grows like a random walk (Brouwer's law).
 *
 * The state is y = (x[0..n-1], v[0..n-1]): 2n values. The acceleration may depend on v.
 *
 * Limitation of the published criterion: tau is relative to the norm of a, so it shrinks to
 * zero where the acceleration of the whole state (of a whole group, when the options name
 * groups) vanishes — every block at once. A system of
 * bodies never gets there, but a one-dimensional oscillator, whose acceleration changes sign,
 * stalls with MYB_E_STEP_MIN: use DOP853 there. */

typedef myb_status (*myb_accel_function)(void* context, double t, const double* x,
                                         const double* v, double* a);

typedef struct myb_ias15 myb_ias15;

typedef struct myb_ias15_options {
  double epsilon;   /* precision parameter, > 0; 0 = 1e-9 */
  double h_initial; /* first step; 0 = automatic (see myb_ias15_start) */
  double h_max;     /* largest step magnitude; 0 = unlimited */
  double h_min;     /* smallest step magnitude before MYB_E_STEP_MIN; 0 = none */
  int64_t max_steps;   /* budget of accepted steps per call of myb_ias15_integrate; 0 = 10^9 */
  size_t block;     /* components per block of the timescale; 0 = one block of all n;
                       must divide n */
  /* Groups of components of different kinds (positions in metres next to the quaternion of
   * a rotation, say): group_count sizes, each a multiple of `block`, their sum n; NULL and 0:
   * one group. The norms that set the step and judge the convergence of the iteration are
   * maxima over the blocks of a group, so that one group measures the others in its own
   * units and the others disappear beside it; with groups each is judged by itself, the
   * shortest timescale and the worst convergence counting. The sizes are copied. */
  const size_t* groups;
  size_t group_count;
} myb_ias15_options;

typedef struct myb_ias15_stats {
  int64_t accepted;
  int64_t rejected;
  int64_t evaluations; /* calls of the acceleration */
  int64_t iterations;  /* predictor-corrector sweeps (7 evaluations each) */
  int64_t unconverged; /* steps whose iteration stopped at the limit of 12 sweeps */
} myb_ias15_stats;

/* Creates an integrator for n coordinates (state of 2n values). Errors: MYB_E_DOMAIN (invalid
 * options, n == 0, block not dividing n), MYB_E_NOMEM, MYB_E_OVERFLOW. allocator may be NULL. */
MYB_API MYB_NODISCARD myb_status myb_ias15_create(size_t n, const myb_ias15_options* options,
                                                  const myb_allocator* allocator,
                                                  myb_ias15** integrator);
MYB_API void myb_ias15_destroy(myb_ias15* integrator); /* NULL is allowed */

/* Sets the problem and the initial state y0 = (x0, v0) (copied); evaluates a(t0) once.
 * Automatic first step: (7! epsilon)^(1/7) times the smallest of sqrt(|x|/|a|) and |v|/|a| over
 * the blocks (zero terms ignored); a first step that is too long is rejected and redone. */
MYB_API MYB_NODISCARD myb_status myb_ias15_start(myb_ias15* integrator, myb_accel_function a,
                                                 void* context, double t0, const double* y0);

/* One accepted step towards t_limit (never beyond it; the direction is the sign of
 * t_limit - t). Errors: MYB_E_STEP_MIN, MYB_E_NOCONV (100 rejections in a row), a status
 * from the acceleration. */
MYB_API MYB_NODISCARD myb_status myb_ias15_step(myb_ias15* integrator, double t_limit);

/* Steps until t_end exactly. */
MYB_API MYB_NODISCARD myb_status myb_ias15_integrate(myb_ias15* integrator, double t_end);

MYB_API double myb_ias15_time(const myb_ias15* integrator);
MYB_API const double* myb_ias15_state(const myb_ias15* integrator); /* 2n values: x, then v */
MYB_API myb_ias15_stats myb_ias15_statistics(const myb_ias15* integrator);

/* Dense output: y = (x, v) at any t within the last accepted step, from the step's polynomial
 * (no evaluations). MYB_E_RANGE if t is outside the last step or no step was taken. */
MYB_API MYB_NODISCARD myb_status myb_ias15_dense(myb_ias15* integrator, double t, double* y);

/* As myb_dop853_locate; g receives the 2n-value state (x, v). */
MYB_API MYB_NODISCARD myb_status myb_ias15_locate(myb_ias15* integrator, myb_event_function g,
                                                  void* context, int* found, double* t_event,
                                                  double* y_event);

/* ---- Symplectic compositions ----------------------------------------------------------------
 * Fixed-step symmetric compositions of the kick-drift-kick leapfrog for x'' = a(t, x) (forces
 * derived from a potential, possibly time-dependent; no velocity dependence): orders 2, 4, 6
 * and 8 by Yoshida's triple jump (H. Yoshida, Phys. Lett. A 150, 262, 1990), with 1, 3, 9, 27
 * evaluations per step; and orders 6 and 8 with 7 and 15 evaluations and smaller errors
 * (MYB_SYMPLECTIC_6_S7: Yoshida's solution A; MYB_SYMPLECTIC_8_S15: the method p8 s15 of
 * Hairer, Lubich and Wanner, Geometric Numerical Integration, V.3). On a Kepler orbit of
 * e = 0.3 their error for equal work is 1.75 times smaller than that of the triple jump of
 * order 6 and 200 times smaller than of order 8 (the work for one accuracy: 10 % and 50 %
 * less). The map is symplectic and time-reversible: the energy error stays
 * bounded over any number of steps instead of drifting, and the phase error grows linearly;
 * the price is a fixed step (no error control). Positions, velocities and time are updated with
 * compensated summation. The state is y = (x[0..n-1], v[0..n-1]). */

typedef myb_status (*myb_position_accel_function)(void* context, double t, const double* x,
                                                  double* a);

typedef struct myb_symplectic myb_symplectic;

typedef struct myb_symplectic_stats {
  int64_t steps;
  int64_t evaluations;
} myb_symplectic_stats;

enum { MYB_SYMPLECTIC_6_S7 = 607, MYB_SYMPLECTIC_8_S15 = 815 };

/* order: 2, 4, 6 or 8 (the triple jump), MYB_SYMPLECTIC_6_S7 or MYB_SYMPLECTIC_8_S15. Errors:
 * MYB_E_DOMAIN (other order, n == 0), MYB_E_NOMEM, MYB_E_OVERFLOW. */
MYB_API MYB_NODISCARD myb_status myb_symplectic_create(size_t n, int order,
                                                       const myb_allocator* allocator,
                                                       myb_symplectic** integrator);
MYB_API void myb_symplectic_destroy(myb_symplectic* integrator); /* NULL is allowed */

/* Sets the problem and y0 = (x0, v0) (copied); evaluates a(t0) once. */
MYB_API MYB_NODISCARD myb_status myb_symplectic_start(myb_symplectic* integrator,
                                                      myb_position_accel_function a,
                                                      void* context, double t0, const double* y0);

/* One step of length h (either sign, not 0). A status from the acceleration leaves the state
 * in the middle of the step: start again. */
MYB_API MYB_NODISCARD myb_status myb_symplectic_step(myb_symplectic* integrator, double h);

/* Equal steps of at most h_max from the current time to exactly t_end. */
MYB_API MYB_NODISCARD myb_status myb_symplectic_integrate(myb_symplectic* integrator,
                                                          double t_end, double h_max);

MYB_API double myb_symplectic_time(const myb_symplectic* integrator);
MYB_API const double* myb_symplectic_state(const myb_symplectic* integrator); /* 2n values */
MYB_API myb_symplectic_stats myb_symplectic_statistics(const myb_symplectic* integrator);

/* ---- RKN12(10) ---------------------------------------------------------------------------------
 * Explicit Runge-Kutta-Nystrom pair of orders 12 and 10 for x'' = a(t, x) (forces without
 * velocity dependence) by J. R. Dormand, M. E. A. El-Mikkawy and P. J. Prince (IMA J. Numer.
 * Anal. 7, 423, 1987): 17 stages, 17 evaluations per accepted step, error estimate from
 * positions and velocities, controller 0.9 err^(-1/11) bounded to [0.2, 5];
 *   err = sqrt(mean_i ((e_i / s_i)^2)),   s_i = atol + rtol max(|y_i|, |y_new_i|)
 * over the 2n components of y = (x, v). Efficient for very tight tolerances. The pair has no
 * interpolant of its own: myb_rkn12_dense takes one more step of the pair from the start of the
 * last step (16 evaluations a request), as accurate as the steps — or, with
 * options.interpolant, builds one on the first request inside a step: the acceleration as a
 * polynomial of degree 12 on the Chebyshev-Lobatto nodes of the step, integrated twice, its
 * inner values refined by Picard's iteration from the quintic of Hermite of the ends until the
 * places settle to a thousandth of the tolerance (11 evaluations a round, usually 4), ending
 * exactly at the end of the step; every further request in the step costs no evaluation (a
 * step too long for the degree, not settled in 12 rounds, falls back to the step of the
 * pair). Inside the steps of a two-body orbit it is as accurate as their ends (1.9e-11 s
 * against 2.2e-11 s at rtol 1e-14); building it costs 44 to 49 evaluations a step, three
 * requests of the other way. */

typedef struct myb_rkn12 myb_rkn12;

typedef struct myb_rkn12_options {
  double rtol;       /* relative tolerance, > 0 */
  double atol;       /* absolute tolerance, > 0 */
  double h_initial;  /* first step; 0 = automatic (Hairer's starting algorithm) */
  double h_max;      /* largest step magnitude; 0 = unlimited */
  double h_min;      /* smallest step magnitude before MYB_E_STEP_MIN; 0 = none */
  int64_t max_steps; /* budget of accepted steps per call of myb_rkn12_integrate; 0 = 10^9 */
  int interpolant;   /* != 0: the dense output from the interpolant of the step (above) */
} myb_rkn12_options;

typedef struct myb_rkn12_stats {
  int64_t accepted;
  int64_t rejected;
  int64_t evaluations; /* calls of the acceleration, dense output included */
} myb_rkn12_stats;

/* Creates an integrator for n coordinates (state of 2n values). Errors: MYB_E_DOMAIN
 * (invalid options, n == 0), MYB_E_NOMEM, MYB_E_OVERFLOW. allocator may be NULL. */
MYB_API MYB_NODISCARD myb_status myb_rkn12_create(size_t n, const myb_rkn12_options* options,
                                                  const myb_allocator* allocator,
                                                  myb_rkn12** integrator);
MYB_API void myb_rkn12_destroy(myb_rkn12* integrator); /* NULL is allowed */

/* Sets the problem and y0 = (x0, v0) (copied); evaluates a(t0) once. */
MYB_API MYB_NODISCARD myb_status myb_rkn12_start(myb_rkn12* integrator,
                                                 myb_position_accel_function a, void* context,
                                                 double t0, const double* y0);

/* One accepted step towards t_limit (never beyond it). Errors: MYB_E_STEP_MIN, MYB_E_NOCONV
 * (100 rejections in a row), a status from the acceleration. */
MYB_API MYB_NODISCARD myb_status myb_rkn12_step(myb_rkn12* integrator, double t_limit);

/* Steps until t_end exactly. */
MYB_API MYB_NODISCARD myb_status myb_rkn12_integrate(myb_rkn12* integrator, double t_end);

MYB_API double myb_rkn12_time(const myb_rkn12* integrator);
MYB_API const double* myb_rkn12_state(const myb_rkn12* integrator); /* 2n values: x, then v */
MYB_API myb_rkn12_stats myb_rkn12_statistics(const myb_rkn12* integrator);

/* State at any t within the last accepted step. MYB_E_RANGE outside it. */
MYB_API MYB_NODISCARD myb_status myb_rkn12_dense(myb_rkn12* integrator, double t, double* y);

/* As myb_dop853_locate; every iteration costs one dense output. */
MYB_API MYB_NODISCARD myb_status myb_rkn12_locate(myb_rkn12* integrator, myb_event_function g,
                                                  void* context, int* found, double* t_event,
                                                  double* y_event);

/* ---- Gauss-Jackson ----------------------------------------------------------------------------
 * Fixed-step multistep predictor-corrector of order 8 for x'' = a(t, x, v) in summed ordinate
 * form (M. M. Berry, L. M. Healy, J. Astronaut. Sci. 52, 331, 2004): positions and velocities
 * come from compensated second and first sums of the accelerations, so round-off grows slowly
 * over long arcs; a step costs one evaluation for the predictor and one per corrector
 * iteration (usually 1-3). The start finds the nodes t0 + j h, j = -4..4, by iterating the
 * method's own formulas from a Taylor guess, so the acceleration must be defined up to 4 steps
 * before t0 as well. Nodes lie on the grid t0 + k h; between them the dense output integrates
 * the window's interpolant of the accelerations exactly (degree 8). No error control: the
 * statistics report the largest predictor-corrector difference of the positions, an estimate
 * of the local error, to choose h. The state is y = (x[0..n-1], v[0..n-1]). */

typedef struct myb_gauss_jackson myb_gauss_jackson;

typedef struct myb_gauss_jackson_options {
  double h;            /* step, nonzero; its sign is the direction of integration */
  double tolerance;    /* corrector convergence (relative change of a); 0 = 1e-15 */
  int max_corrections; /* corrector iterations per step; 0 = 10 */
} myb_gauss_jackson_options;

typedef struct myb_gauss_jackson_stats {
  int64_t steps;
  int64_t evaluations;
  int64_t unconverged;   /* steps whose corrector stopped at max_corrections */
  double max_difference; /* largest |x_corrected - x_predicted| / |x| of a step */
} myb_gauss_jackson_stats;

/* Errors: MYB_E_DOMAIN (h = 0 or not finite, n == 0, negative options), MYB_E_NOMEM,
 * MYB_E_OVERFLOW. allocator may be NULL. */
MYB_API MYB_NODISCARD myb_status myb_gauss_jackson_create(
    size_t n, const myb_gauss_jackson_options* options, const myb_allocator* allocator,
    myb_gauss_jackson** integrator);
MYB_API void myb_gauss_jackson_destroy(myb_gauss_jackson* integrator); /* NULL is allowed */

/* Sets the problem and y0 = (x0, v0) at t0 (copied) and runs the start. Errors: a status from
 * the acceleration, MYB_E_NOCONV (the start iteration does not converge: h is too long). */
MYB_API MYB_NODISCARD myb_status myb_gauss_jackson_start(myb_gauss_jackson* integrator,
                                                         myb_accel_function a, void* context,
                                                         double t0, const double* y0);

/* Advances the state to the next node (one step of h). */
MYB_API MYB_NODISCARD myb_status myb_gauss_jackson_step(myb_gauss_jackson* integrator);

/* Steps until the node reaches t_end or just passes it; the state at t_end itself then comes
 * from myb_gauss_jackson_dense. MYB_E_DOMAIN if t_end lies behind the current node. */
MYB_API MYB_NODISCARD myb_status myb_gauss_jackson_integrate(myb_gauss_jackson* integrator,
                                                             double t_end);

MYB_API double myb_gauss_jackson_time(const myb_gauss_jackson* integrator);
MYB_API const double* myb_gauss_jackson_state(const myb_gauss_jackson* integrator); /* 2n */
MYB_API myb_gauss_jackson_stats myb_gauss_jackson_statistics(
    const myb_gauss_jackson* integrator);

/* State at any t within the window: up to 4 steps around its centre, which covers the last
 * 4 steps (MYB_E_RANGE otherwise). No evaluations. */
MYB_API MYB_NODISCARD myb_status myb_gauss_jackson_dense(myb_gauss_jackson* integrator, double t,
                                                         double* y);

/* As myb_dop853_locate, within the last step. */
MYB_API MYB_NODISCARD myb_status myb_gauss_jackson_locate(myb_gauss_jackson* integrator,
                                                          myb_event_function g, void* context,
                                                          int* found, double* t_event,
                                                          double* y_event);

#ifdef __cplusplus
}
#endif

#endif /* MYB_INTEG_H_ */
