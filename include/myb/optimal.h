/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/optimal.h — optimal transfers with low thrust: the indirect method (the conditions of
 * Pontryagin's principle solved by shooting), the direct method of Sims and Flanagan, chains
 * of its legs joined by flybys, and the small solver of nonlinear programs behind the latter.
 *
 * The engine gives a thrust up to T (N) with the effective exhaust velocity c (m/s):
 *   x'' = g(t, x, x') + (T u / m) d,   m' = -T u / c,   0 <= u <= 1,  |d| = 1.
 * The methods join two given states — a rendezvous — or, the indirect one, a state and an
 * orbit; the chains join bodies. What comes out is a local optimum: the one the first guess
 * leads to. */
#ifndef MYB_OPTIMAL_H_
#define MYB_OPTIMAL_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/kin.h"
#include "myb/model.h"
#include "myb/time.h"

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
 * of the conditions is taken by central differences.
 *
 * Insertion into an orbit with the phase free (MYB_INDIRECT_ORBIT): the arrival may be at any
 * point of the conic through `target` about mu, y_T(s), the state of that conic a time s after
 * `target`. The phase s is one more unknown, and one more condition takes its place, the
 * transversality condition of Pontryagin's principle: the costates at the arrival are
 * orthogonal to the tangent of the orbit,
 *   l_x . v_T + l_v . a_T = 0,   (v_T, a_T) = dy_T / ds = (velocity, -mu r_T / |r_T|^3),
 * moving the point of arrival along the orbit gains nothing to the first order. Five elements
 * of the arrival are fixed and the sixth, where on the orbit, is chosen by the method. For the
 * least time the orbit does not move with the time, and the condition on the Hamiltonian at
 * the arrival becomes H = 0. The orbit is an ellipse; whatever `field` holds, it is the conic
 * of mu that is the target. Measured: from 7000 km to the circle of 8000 km turned by
 * 0.06 rad, in 4500 s at 0.3 m/s^2 (eps = 1e-4), the method arrives 4.3546 rad along the circle
 * and leaves 783.920 kg, against 780.076 for the rendezvous at 4.4 rad and 744.5 .. 777.8 at
 * four other points within 0.3 rad; the rendezvous at the point chosen gives the same mass to
 * 2e-7 kg; the transversality condition from the costates of the samples 3e-12 of its terms,
 * the Hamiltonian constant to 2e-11. The least time into the circle of 8000 km is 3011.759 s,
 * against 3633.451 s to a target that moves on it, with H = 0 to 2e-15 all the way. */

typedef enum myb_indirect_goal {
  MYB_INDIRECT_FUEL = 0, /* the least propellant in a given time */
  MYB_INDIRECT_TIME = 1  /* the least time */
} myb_indirect_goal;

typedef enum myb_indirect_arrival {
  MYB_INDIRECT_RENDEZVOUS = 0, /* the state `target` (FUEL) or the target that moves (TIME) */
  MYB_INDIRECT_ORBIT = 1       /* anywhere on the orbit through `target`: the phase free */
} myb_indirect_arrival;

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
  myb_indirect_arrival arrival; /* RENDEZVOUS (the default) or ORBIT: `target` is then any
                                   state of the orbit, an ellipse about mu */
  double phase;           /* ORBIT, with a guess of the costates: a first guess of the phase
                             (s; as myb_indirect_solution.phase) */
} myb_indirect_problem;

typedef struct myb_indirect_solution {
  double costate[8]; /* at t0: l_0, l_x, l_v, l_m; a unit vector, l_0 > 0 */
  double t1;         /* the time of the arrival (s) */
  double mass1;      /* kg there */
  double smoothing;  /* FUEL: the eps of this solution */
  double residual;   /* the largest of the conditions left */
  int tries;         /* first guesses used */
  int stages;        /* FUEL: steps taken from eps = 1 */
  double phase;      /* ORBIT: where the transfer ends — the time along the target orbit from
                        `target` (s), in [0, its period) */
} myb_indirect_solution;

/* Solves the problem. guess (may be NULL): the eight costates of a solution nearby — of the
 * same problem with other numbers, at the eps asked — from which the conditions are solved
 * at once (for ORBIT with problem->phase the first guess of the phase); without one, first
 * guesses are drawn from a fixed sequence (the same result on every run and platform) until
 * one leads to a solution, for FUEL at eps = 1, which is then followed to the eps asked; for
 * ORBIT the phase is then first guessed as the point of the orbit in the direction halfway
 * between where the start would coast in the time given and where the orbit would carry the
 * point of it nearest to the start. Errors: MYB_E_DOMAIN (invalid numbers, a target orbit
 * that is no ellipse), MYB_E_NOMEM, MYB_E_NOCONV (no first guess led anywhere — the transfer
 * may be impossible in the time given; or the way to the eps asked broke off: *solution then
 * holds the last one found, its `smoothing` saying where). allocator may be NULL. */
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
  double costate[7];   /* l_x, l_v, l_m in the units of the solver (-l_v is the primer
                          vector) */
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
  double weight;      /* the first rho, at most 1e12; 0 = 10 max(1, |f|) / max(1, the half sum
                         of the squares of the violations) at the first point, in 1e-8 .. 1e8.
                         From a point that is nearly a solution, or with an objective that
                         would pull the first steps far from the constraints, a larger one
                         keeps them close (myb_chain_solve takes 1e4) */
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
 * bounds that cross, negative options, a weight above 1e12, a function that is not finite),
 * MYB_E_NOMEM,
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
 * direction, so that the engine being off is a bound and not a corner. As the segments grow
 * in number the result tends to the continuous optimum (that of the indirect method above);
 * with 30 a per cent of the propellant or so is typical.
 *
 * The field between the impulses is the point mass mu — the conics — or a force model: then
 * each half segment is integrated by DOP853 (tolerance 1e-13 in the units of the solver)
 * together with its transition matrix, from the variational equations Phi' = [0 I; dg/dx
 * dg/dv] Phi with the Jacobian of the model, and the derivatives of the mismatch are those
 * transition matrices multiplied from the meeting point outward. The model is evaluated at
 * its own time, t0 + the time along the leg. Measured on the transfer of the indirect method
 * with 30 segments: the point mass as a model gives the mismatch of the conics to 8e-16 of
 * the orbit and the same least propellant, 220.176 kg (9e-7 kg apart); with the J2 of the
 * Earth added random throttles meet 41 km elsewhere, the least propellant is 219.892 kg, the
 * halves meet to 4e-3 m and 7e-6 m/s, and the Jacobian of the mismatch agrees with central
 * differences to 1.5e-7 of its size (9e-8 for the conics). */

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
  myb_force_model* field; /* everything but the engine, the centre included (the transition
                             matrices use its Jacobian); NULL: the point mass mu, the conics.
                             mu still sets the units of the solver */
} myb_sims_flanagan;

/* The distance at the meeting point between the two halves — forward minus backward:
 * position (m) and velocity (m/s) — for the throttles given (3 per segment), the mass at the
 * arrival *mass1 (may be NULL), and the Jacobian of the mismatch by the throttles (6 x 3
 * segments, row-major; may be NULL; where a throttle is zero the derivative of its size is
 * taken as zero). With a field the integrator and its work are allocated by malloc. Errors:
 * MYB_E_DOMAIN (invalid numbers, a throttle longer than 1 + 1e-9, more propellant than mass),
 * those of the Kepler propagation, with a field MYB_E_NOMEM and those of the integrator and
 * of the model. */
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

/* ---- Chains of legs with flybys -----------------------------------------------------------------
 * Legs of Sims and Flanagan from body to body — the departure, the bodies flown by, the
 * arrival — joined by flybys, all of it one nonlinear program for myb_nlp_solve (the flybys
 * are part of the method from its start: J. A. Sims, S. N. Flanagan, AAS 99-338, 1999). A body
 * is a provider of its state relative to the centre; a leg k starts at body k at
 * the time t_k with the body's velocity plus an excess velocity w_k^+ and ends at body k + 1 at
 * t_k + T_k with the body's velocity plus w_{k+1}^-. A flyby is a turn of the excess velocity
 * in the field of the body alone, quicker than anything the engine does: unpowered, it keeps
 * the excess speed,
 *   |w^-| - |w^+| = 0,
 * and it turns the excess velocity by at most the angle of a hyperbola that passes the body at
 * the least distance r_p allowed, sin(delta_max / 2) = 1 / e, e = 1 + r_p |w^-| |w^+| / gm:
 *   1 - 2 / e^2 - w^- . w^+ / (|w^-| |w^+|) <= 0   (cos delta >= cos delta_max).
 * The engine is the same on every leg and its mass carries from leg to leg; the launcher gives
 * an excess speed of at most `launch` at the departure, and at most `arrival` is left at the
 * arrival (0: a rendezvous). The speeds enter as themselves and not as their squares: the
 * gradient of a square is as small as the speed (0.08 in the units of the solver for 2.5 km/s
 * about the Sun), and the solver had to raise its weight a thousandfold before it held such a
 * constraint to 1e-9.
 *
 * The unknowns, in this order (myb_chain_size gives their number):
 *   the time of the departure t_0 (s), the time of flight of each leg T_k (s),
 *   the excess velocity at the departure (m/s), then at each body flown by the excess velocity
 *   coming in and going out, then the excess velocity at the arrival — 2 legs vectors of 3,
 *   the throttles of the segments of each leg in turn, 3 per segment (as Sims and Flanagan).
 * The times lie between bounds (equal bounds hold one fixed). The constraints: the mismatch of
 * each leg (6: m, m/s), the excess speed of each flyby (m/s), then the inequalities: the turn
 * of each flyby, the excess speed at the departure and at the arrival less its limit (m/s).
 * What is minimised is the propellant (kg). Inside, the units are those of Sims and Flanagan:
 * the distance of the first body from the centre at the earliest departure, the time
 * sqrt(length^3 / mu), the first mass. The derivatives are analytic but for the acceleration
 * of a body, taken from its provider's velocity by central differences (a ten-thousandth of
 * the unit of time), and they go through the time of every impulse: a segment that grows
 * longer moves its impulses, its end and the bodies, the impulses grow with it, and in a field
 * that changes with time the whole leg moved in time is another leg.
 *
 * The solver takes the throttles in two passes. The first sees each as a vector tau, the
 * impulse q tau, with a size u of its own for the masses and |tau|^2 <= u^2: the direction of
 * an engine that is off is not lost. As a size and two angles (the form of
 * myb_sims_flanagan_solve) a throttle that reaches zero loses the gradient of its direction:
 * from the built-in guess the first steps turned the engine off over most of the first leg
 * pointing the wrong way, and no weight brought it back (the departure stayed 278 m/s above
 * the launcher). But the cone is degenerate where the engine is off, and there the solver
 * crawls; so the first pass stops at 1e-6 after 10 updates of 300 steps, and the second takes
 * a size and two angles from where it left, with the directions right. Both start with the
 * weight 1e4 (myb_nlp_options.weight): with the weight chosen from the first point (10) the
 * propellant pulled the first steps far from the constraints.
 *
 * Measured (test_chain): the Earth, Venus and Mars on circles, Mars in a plane 1.85 degrees
 * off, 1 N on 1000 kg at 30 km/s, 2.5 km/s from the launcher, a rendezvous with Mars, 10
 * segments a leg, the times free within 20 days: the derivatives against central differences
 * 1.4e-8 of their sizes (1.5e-8 in a field with a companion of a tenth of the Sun's mass at
 * 3 au, which moves with time), the propellant 4.7e-11; the Sun as a model meets as the
 * conics to 3.5e-14. By Venus: 834.470 kg left (5429 m/s), the legs meet to 30 m and 6e-6 m/s,
 * the excess speeds at Venus 3e-6 m/s apart, turned by 1.2112 rad of the 1.2863 allowed (7201
 * km from the centre), the launcher's 2500 m/s used to 7e-8 m/s; 10 + 5 updates of the
 * multipliers, 14 000 evaluations. The same chain kept 10^6 km from Venus (a turn of 0.018 rad
 * at most) leaves 731.350 kg (9386 m/s) from the guess and from the chain by Venus alike: the
 * flyby saves 103 kg. Found on the way: the direct transfer in these windows — Mars 59 degrees
 * ahead after 390 days — is hard for the solver: coplanar it gave 712.8 kg at the middle of
 * the windows and 763.7 kg with the dates free (2.7 km/s more than by Venus), with Mars off
 * the plane it ended in a local optimum that spends four fifths of the mass. */

enum { MYB_CHAIN_MAX_LEGS = 8 };

typedef struct myb_chain_body {
  myb_state_function state; /* x, v relative to the centre (m, m/s) in the axes of the legs,
                               asked at epoch + t (myb/kin.h) */
  void* context;
  double gm;        /* m^3 s^-2: of a body flown by */
  double periapsis; /* m: the least distance of a flyby from the centre of the body */
} myb_chain_body;

typedef struct myb_chain {
  double mu;              /* of the centre (m^3 s^-2): the field when `field` is NULL, and the
                             units of the solver */
  myb_force_model* field; /* of the legs, everything but the engine; NULL: the point mass mu.
                             Its time is the time of the chain */
  myb_time epoch;         /* the instant of the time 0 of the chain */
  double thrust;          /* N */
  double exhaust;         /* m/s */
  double mass0;           /* kg at the departure */
  size_t legs;            /* 1 .. MYB_CHAIN_MAX_LEGS */
  size_t segments[MYB_CHAIN_MAX_LEGS]; /* of each leg, 2 .. MYB_SIMS_FLANAGAN_MAX_SEGMENTS */
  myb_chain_body bodies[MYB_CHAIN_MAX_LEGS + 1]; /* the departure, those flown by, the
                                                     arrival: legs + 1 of them */
  double departure[2];    /* the earliest and the latest time of the departure (s) */
  double flight[MYB_CHAIN_MAX_LEGS][2]; /* the shortest and the longest flight of each leg (s) */
  double launch;          /* the largest excess speed at the departure (m/s), >= 0 */
  double arrival;         /* the largest excess speed at the arrival (m/s); 0: a rendezvous */
} myb_chain;

/* The number of the unknowns (7 legs + 1 + 3 segments in all) and of the constraints (6 legs +
 * legs - 1 equalities, legs + 1 inequalities). Errors: MYB_E_DOMAIN (sizes, bounds, values that
 * are not finite or negative, a body without a provider, a body flown by without gm or least
 * distance, more than MYB_NLP_MAX_UNKNOWNS unknowns or MYB_NLP_MAX_CONSTRAINTS constraints in
 * the first pass of the solver: 7 legs + 1 + 4 segments and 8 legs + segments). Any output may
 * be NULL. */
MYB_API MYB_NODISCARD myb_status myb_chain_size(const myb_chain* chain, size_t* unknowns,
                                                size_t* equalities, size_t* inequalities);

/* The first guess built in: the times in the middle of their bounds, the excess velocities
 * those of the conics of Lambert (zero revolutions, in the sense of the motion of the first
 * body) between the bodies at those times, cut down to `launch` and `arrival` at the ends, the
 * engine off. x: myb_chain_size unknowns. Errors: those of myb_chain_size, of the providers
 * and of the problem of Lambert. */
MYB_API MYB_NODISCARD myb_status myb_chain_guess(const myb_chain* chain, double* x);

/* The propellant spent (kg; may be NULL) and the constraints (may be NULL) at x, with, when
 * not NULL, the gradient of the propellant (unknowns values) and the Jacobian of the
 * constraints (constraints x unknowns, row-major), by the unknowns in their units. Where a
 * throttle is zero the derivative of its size is taken as zero. allocator may be NULL.
 * Errors: those of myb_chain_size, MYB_E_DOMAIN (a throttle longer than 1 + 1e-9, times out of
 * the order of the legs, more propellant than mass, a flyby without excess speed), MYB_E_NOMEM,
 * those of the providers, the propagation and the model. */
MYB_API MYB_NODISCARD myb_status myb_chain_evaluate(const myb_chain* chain, const double* x,
                                                    double* propellant, double* constraints,
                                                    double* gradient, double* jacobian,
                                                    const myb_allocator* allocator);

typedef struct myb_chain_result {
  double mass[MYB_CHAIN_MAX_LEGS + 1];      /* kg at each body, mass0 at the departure */
  double excess[MYB_CHAIN_MAX_LEGS + 1][2]; /* m/s: the excess speed at each body coming in and
                                               going out (zero where there is none) */
  double turn[MYB_CHAIN_MAX_LEGS + 1];      /* rad: of the excess velocity at each flyby */
  double periapsis[MYB_CHAIN_MAX_LEGS + 1]; /* m: the distance of the hyperbola of that turn
                                               (zero at the ends, HUGE_VAL for no turn) */
  double mismatch[MYB_CHAIN_MAX_LEGS][6];   /* left at the meeting point of each leg (m, m/s) */
  myb_nlp_result first; /* the first pass, the throttles as vectors */
  myb_nlp_result nlp;   /* the second, a size and two angles each */
} myb_chain_result;

/* The chain with the least propellant. x: the first guess on entry — all zeros ask for that of
 * myb_chain_guess; with no throttle given the solver starts from half thrust along the motion
 * of the coast from the nearer end of each leg — the solution on return. options (may be NULL:
 * 1e-9 and 1e-6, the weight 1e4) are those of the second pass, in the units of the solver; a
 * weight given there is that of both passes. Errors: those of myb_chain_evaluate,
 * MYB_E_NOCONV (x and *result hold what was reached). allocator may be NULL. */
MYB_API MYB_NODISCARD myb_status myb_chain_solve(const myb_chain* chain,
                                                 const myb_nlp_options* options, double* x,
                                                 const myb_allocator* allocator,
                                                 myb_chain_result* result);

#ifdef __cplusplus
}
#endif

#endif /* MYB_OPTIMAL_H_ */
