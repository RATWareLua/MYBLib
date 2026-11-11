/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/optimal.h — optimal transfers with low thrust: the indirect method (the conditions of
 * Pontryagin's principle solved by shooting), the direct method of Sims and Flanagan, and the
 * small solver of nonlinear programs behind the latter.
 *
 * The engine gives a thrust up to T (N) with the effective exhaust velocity c (m/s):
 *   x'' = g(t, x, x') + (T u / m) d,   m' = -T u / c,   0 <= u <= 1,  |d| = 1.
 * Both methods join two given states — a rendezvous. What comes out is a local optimum: the
 * one the first guess leads to. */
#ifndef MYB_OPTIMAL_H_
#define MYB_OPTIMAL_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/model.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- The indirect method -----------------------------------------------------------------------
 * With the costates l_x, l_v, l_m of the position, the velocity and the mass, the thrust
 * that makes the Hamiltonian least points against l_v, and for the least propellant
 *   H = l_x . v + l_v . g - (T u / m) |l_v| - l_m T u / c + l_0 (T / c) [u - eps u (1 - u)]
 * is least in u at
 *   u = 0 for S > eps,  1 for S < -eps,  (eps - S) / (2 eps) between,
 *   S = 1 - (c |l_v| / m + l_m) / l_0   (the switching function).
 * eps = 1 is the transfer of the least "energy" (the integral of u^2), smooth and easy to
 * find; eps -> 0 the transfer of the least propellant, with the engine either full on or off;
 * the solver goes from the first to the second, each solution the first guess of the next
 * (R. Bertrand, R. Epenoy, Optim. Control Appl. Meth. 23, 171, 2002). The costates obey
 *   l_x' = -(dg/dx)^T l_v,  l_v' = -l_x - (dg/dv)^T l_v,  l_m' = -T u |l_v| / m^2,
 * and the unknowns are their values at the start together with l_0, the multiplier of the
 * cost, all eight made a unit vector: then a first guess taken at random on the sphere has a
 * fair chance (F. Jiang, H. Baoyin, J. Li, J. Guid. Control Dyn. 35, 245, 2012). The
 * conditions are the state wanted at the arrival and l_m = 0 there (the mass left is free).
 * For the least time the engine is always on and the time of the arrival is one more
 * unknown, with the condition H = l_x . v_T + l_v . a_T at the arrival for a target that
 * moves with the velocity v_T and the acceleration a_T.
 *
 * Inside, lengths are in units of the first distance from the centre, times in units of
 * sqrt(length^3 / mu), masses in units of the first mass: the costates, the Hamiltonian and
 * the tolerance are in those units. The equations are integrated by DOP853 and the Jacobian
 * of the conditions is taken by central differences. */

typedef enum myb_indirect_goal {
  MYB_INDIRECT_FUEL = 0, /* the least propellant in a given time */
  MYB_INDIRECT_TIME = 1  /* the least time */
} myb_indirect_goal;

typedef struct myb_indirect_problem {
  double mu;              /* of the centre (m^3 s^-2): the field when `field` is NULL, and the
                             units of the solver */
  myb_force_model* field; /* everything but the engine, the centre included (the costates use
                             its Jacobian); NULL: the point mass mu */
  double thrust;          /* N: the most the engine gives */
  double exhaust;         /* m/s */
  double t0;              /* s */
  double start[6];        /* x, v at t0 (m, m/s) */
  double mass0;           /* kg */
  double t1;              /* FUEL: the time of the arrival; TIME: a first guess of it */
  double target[6];       /* FUEL: x, v wanted at t1; TIME: the state at t0 of a target that
                             moves on its conic about mu */
  myb_indirect_goal goal;
  double smoothing;       /* FUEL: the eps to end at, in (0, 1]; 0 = 1e-4 */
  double tolerance;       /* on the conditions at the arrival; 0 = 1e-9 */
  int tries;              /* first guesses to try when none is given; 0 = 100 */
} myb_indirect_problem;

typedef struct myb_indirect_solution {
  double costate[8]; /* at t0: l_0, l_x, l_v, l_m; a unit vector, l_0 > 0 */
  double t1;         /* the time of the arrival (s) */
  double mass1;      /* kg there */
  double smoothing;  /* FUEL: the eps of this solution */
  double residual;   /* the largest of the conditions left */
  int tries;         /* first guesses used */
  int stages;        /* FUEL: steps taken from eps = 1 */
} myb_indirect_solution;

/* Solves the problem. guess (may be NULL): the eight costates of a solution nearby — of the
 * same problem with other numbers, at the eps asked — from which the conditions are solved
 * at once; without one, first guesses are drawn from a fixed sequence (the same result on
 * every run and platform) until one leads to a solution, for FUEL at eps = 1, which is then
 * followed to the eps asked. Errors: MYB_E_DOMAIN (invalid numbers), MYB_E_NOMEM,
 * MYB_E_NOCONV (no first guess led anywhere — the transfer may be impossible in the time
 * given; or the way to the eps asked broke off: *solution then holds the last one found, its
 * `smoothing` saying where). allocator may be NULL. */
MYB_API MYB_NODISCARD myb_status myb_indirect_solve(const myb_indirect_problem* problem,
                                                    const double* guess,
                                                    const myb_allocator* allocator,
                                                    myb_indirect_solution* solution);

typedef struct myb_indirect_point {
  double y[7];         /* x, v (m, m/s) and the mass (kg) */
  double direction[3]; /* of the thrust, unit */
  double throttle;     /* u, 0 .. 1 */
  double switching;    /* S: the engine burns where it is negative (FUEL) */
  double hamiltonian;  /* constant along a solution in a field that does not depend on time */
} myb_indirect_point;

/* The solution at `count` times in increasing order between t0 and the arrival. Errors:
 * MYB_E_DOMAIN (times out of order or outside the transfer), MYB_E_NOMEM, those of the
 * field. */
MYB_API MYB_NODISCARD myb_status myb_indirect_sample(const myb_indirect_problem* problem,
                                                     const myb_indirect_solution* solution,
                                                     const double* times, size_t count,
                                                     const myb_allocator* allocator,
                                                     myb_indirect_point* points);

/* ---- A solver of nonlinear programs ------------------------------------------------------------
 *   minimise f(x)  with  c_i(x) = 0 for the first `equalities` constraints,
 *                        c_i(x) <= 0 for the `inequalities` after them,  lower <= x <= upper,
 * for up to a few hundred unknowns, with dense matrices: what the direct method below needs,
 * and no replacement for the large sparse solvers. The method of multipliers (M. R. Hestenes,
 * M. J. D. Powell, 1969; R. T. Rockafellar, 1973, for the inequalities): the augmented
 * Lagrangian
 *   L = f + (rho / 2) [sum_eq (c_i + l_i / rho)^2 + sum_ineq max(0, c_i + l_i / rho)^2]
 * is minimised in the box with the multipliers l fixed, then l moves by rho c, and rho grows
 * tenfold whenever the constraints have not come twice closer (the scheme of E. G. Birgin,
 * J. M. Martinez, Practical Augmented Lagrangian Methods, 2014). Inside, a quasi-Newton
 * iteration with the damped BFGS matrix of M. J. D. Powell (1978): Newton's step for the
 * variables that are free, those that the gradient pushes against a bound stay on it, and
 * the step is halved along its projection into the box until the function falls enough
 * (D. P. Bertsekas, SIAM J. Control Optim. 20, 221, 1982). The unknowns and the functions
 * should be of the order of one. */

enum { MYB_NLP_MAX_UNKNOWNS = 256, MYB_NLP_MAX_CONSTRAINTS = 256 };

/* The objective *f and the constraints c (equalities first) at x; when gradient is not NULL
 * also the gradient of f (n values) and the Jacobian of c (constraints x n, row-major). Any
 * status but MYB_OK at the first point ends the solution with it; later it makes the solver
 * shorten its step. */
typedef myb_status (*myb_nlp_function)(void* context, const double* x, double* f, double* c,
                                       double* gradient, double* jacobian);

typedef struct myb_nlp_problem {
  size_t n;            /* unknowns, 1 .. MYB_NLP_MAX_UNKNOWNS */
  size_t equalities;   /* together with the inequalities at most MYB_NLP_MAX_CONSTRAINTS */
  size_t inequalities;
  const double* lower; /* n bounds, -HUGE_VAL for none; NULL: none at all */
  const double* upper; /* n bounds, HUGE_VAL for none; NULL: none at all */
  myb_nlp_function function;
  void* context;
} myb_nlp_problem;

typedef struct myb_nlp_options {
  double feasibility; /* on the constraints left (and c_i l_i of the inequalities); 0 = 1e-8 */
  double optimality;  /* on the largest component of x - P(x - grad L), P the projection into
                         the box, L the Lagrangian; 0 = 1e-6 */
  int max_outer;      /* updates of the multipliers; 0 = 50 */
  int max_inner;      /* steps between two updates; 0 = 1000 */
} myb_nlp_options;

typedef struct myb_nlp_result {
  double objective;
  double violation;    /* of the constraints, as `feasibility` measures it */
  double stationarity; /* as `optimality` measures it */
  int outer;
  int inner;           /* steps in all */
  int evaluations;     /* calls of the function */
} myb_nlp_result;

/* x: the first guess on entry (moved into the box), the solution on return; multipliers (may
 * be NULL): those of the constraints. options may be NULL. Errors: MYB_E_DOMAIN (sizes,
 * bounds that cross, negative options, a function that is not finite), MYB_E_NOMEM,
 * MYB_E_NOCONV (x and *result hold what was reached), the statuses of the function at the
 * first point. allocator may be NULL. */
MYB_API MYB_NODISCARD myb_status myb_nlp_solve(const myb_nlp_problem* problem,
                                               const myb_nlp_options* options, double* x,
                                               double* multipliers,
                                               const myb_allocator* allocator,
                                               myb_nlp_result* result);

/* ---- The direct method of Sims and Flanagan ----------------------------------------------------
 * The transfer is cut into equal segments in time, and the thrust of each is replaced by one
 * impulse in its middle, of at most what the engine gives in a segment; between the impulses
 * the motion is a conic (J. A. Sims, S. N. Flanagan, AAS 99-338, 1999). The first half of
 * the segments is flown forward from the start, the rest backward from the arrival, and the
 * two must meet: six conditions on the throttles, tau_i with |tau_i| <= 1, one vector for a
 * segment,
 *   dv_i = tau_i T dt / m_i,   m_i = m_0 - (T dt / c) (|tau_1| + ... + |tau_i| / 2),
 * the mass in the middle of the segment, which falls by |tau_i| T dt / c over it. Among the
 * throttles that meet, the ones that spend the least are looked for with myb_nlp_solve: the
 * unknowns there are the size of each throttle between 0 and 1 and two angles of its
 * direction, so that the engine being off is a bound and not a corner. The field is the point
 * mass mu. As the segments grow in number the result tends to the continuous optimum (that
 * of the indirect method above); with 30 a per cent of the propellant or so is typical. */

enum { MYB_SIMS_FLANAGAN_MAX_SEGMENTS = 64 };

typedef struct myb_sims_flanagan {
  double mu;       /* m^3 s^-2 */
  double thrust;   /* N */
  double exhaust;  /* m/s */
  double t0;       /* s */
  double start[6]; /* x, v at t0 */
  double mass0;    /* kg */
  double t1;       /* s, > t0 */
  double target[6]; /* x, v at t1 */
  size_t segments; /* 2 .. MYB_SIMS_FLANAGAN_MAX_SEGMENTS */
} myb_sims_flanagan;

/* The distance at the meeting point between the two halves — forward minus backward:
 * position (m) and velocity (m/s) — for the throttles given (3 per segment), the mass at the
 * arrival *mass1 (may be NULL), and the Jacobian of the mismatch by the throttles (6 x 3
 * segments, row-major; may be NULL; where a throttle is zero the derivative of its size is
 * taken as zero). Errors: MYB_E_DOMAIN (invalid numbers, a throttle longer than 1 + 1e-9,
 * more propellant than mass), those of the Kepler propagation. */
MYB_API MYB_NODISCARD myb_status myb_sims_flanagan_mismatch(const myb_sims_flanagan* leg,
                                                            const double* throttles,
                                                            double mismatch[6], double* mass1,
                                                            double* jacobian);

typedef struct myb_sims_flanagan_result {
  double mass1;        /* kg at the arrival */
  double mismatch[6];  /* left at the meeting point (m, m/s) */
  myb_nlp_result nlp;
} myb_sims_flanagan_result;

/* The throttles that join the two states with the least propellant. throttles (3 per
 * segment): the first guess on entry — all zeros ask for the one built in, half thrust along
 * the velocity of the two coasts — the solution on return. options (may be NULL: 1e-9 and
 * 1e-6) are in the units of the solver: the first distance from the centre, the time
 * sqrt(length^3 / mu), the first mass. Errors: MYB_E_DOMAIN, MYB_E_NOMEM, MYB_E_NOCONV (the
 * throttles and *result hold what was reached). allocator may be NULL. */
MYB_API MYB_NODISCARD myb_status myb_sims_flanagan_solve(const myb_sims_flanagan* leg,
                                                         const myb_nlp_options* options,
                                                         double* throttles,
                                                         const myb_allocator* allocator,
                                                         myb_sims_flanagan_result* result);

#ifdef __cplusplus
}
#endif

#endif /* MYB_OPTIMAL_H_ */
