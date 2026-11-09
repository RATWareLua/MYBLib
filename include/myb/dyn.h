/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/dyn.h — systems of bodies: the Newtonian N-body problem with an interaction graph, its
 * invariants, and the Wisdom-Holman family of splitting integrators (WHFast, SABA, SBAB).
 *
 * Bodies are given by GM (m^3 s^-2) only. States are arrays of 3 values per body (x, y, z of
 * body 0, then body 1, ...; m, m/s, m/s^2) in an inertial frame. Bodies 0 .. massive-1 attract;
 * the following test bodies feel the massive ones and attract nothing.
 *
 * Accelerations are direct sums over all pairs, without approximation: for receiver i
 *   a_i = sum_j GM_j (x_j - x_i) / |x_j - x_i|^3,
 * the terms added in increasing j with compensated summation (Sum2: Ogita, Rump, Oishi, SIAM
 * J. Sci. Comput. 26, 1955, 2005), so the sum is accurate to about one ulp whatever the number
 * and spread of the terms. The order of operations is fixed per receiver, so the result is the
 * same, bit for bit, for the symmetric single-threaded evaluation (each pair computed once,
 * used for both bodies) and for the per-receiver one (four receivers per vector block, without
 * FMA, each receiver summed by one thread), for any number of threads and any ISA level.
 *
 * Sources: J. Wisdom, M. Holman, AJ 102, 1528 (1991); J. Wisdom, M. Holman, J. Touma, Fields
 * Inst. Commun. 10, 217 (1996) (symplectic correctors); H. Rein, D. Tamayo, MNRAS 452, 376
 * (2015) (WHFast); J. Laskar, P. Robutel, Celest. Mech. Dyn. Astron. 80, 39 (2001) (SABA,
 * SBAB). */
#ifndef MYB_DYN_H_
#define MYB_DYN_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"
#include "myb/num.h"
#include "myb/pool.h"
#include "myb/rel.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- The N-body system --------------------------------------------------------------------- */

typedef struct myb_nbody myb_nbody;

/* How a source body acts on a receiver. Later model levels (gravity fields, tides) extend this
 * list. The default for every pair is MYB_INTERACTION_NEWTON, or MYB_INTERACTION_EIH once
 * myb_nbody_set_relativity has been called; pairs set explicitly keep their level. */
typedef enum myb_interaction {
  MYB_INTERACTION_NONE = 0,   /* the source does not act on the receiver */
  MYB_INTERACTION_NEWTON = 1, /* point masses */
  MYB_INTERACTION_EIH = 2     /* point masses with the 1PN terms (below) */
} myb_interaction;

/* massive >= 1 bodies with gm[0 .. massive-1] (finite, >= 0) and test >= 0 test bodies.
 * Errors: MYB_E_DOMAIN (massive == 0, a GM negative or not finite), MYB_E_OVERFLOW,
 * MYB_E_NOMEM. A system is used by one thread at a time. */
MYB_API MYB_NODISCARD myb_status myb_nbody_create(size_t massive, size_t test, const double* gm,
                                                  const myb_allocator* allocator,
                                                  myb_nbody** nbody);

/* Accepts NULL. */
MYB_API void myb_nbody_destroy(myb_nbody* nbody);

MYB_API size_t myb_nbody_count(const myb_nbody* nbody);   /* massive + test */
MYB_API size_t myb_nbody_massive(const myb_nbody* nbody); /* massive */
MYB_API double myb_nbody_gm(const myb_nbody* nbody, size_t body); /* 0 for test bodies */

/* Threads for the accelerations (NULL: the calling thread only). The pool must outlive its use
 * by the system; the results do not depend on it. */
MYB_API void myb_nbody_set_pool(myb_nbody* nbody, myb_pool* pool);

/* Level of the action of source on receiver (a directed edge: excluding only one direction of
 * a pair breaks the conservation of momentum, which is then the caller's choice). The source
 * must be massive and differ from the receiver. Errors: MYB_E_DOMAIN (indices, level,
 * MYB_INTERACTION_EIH without myb_nbody_set_relativity), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_nbody_set_interaction(myb_nbody* nbody, size_t receiver,
                                                           size_t source, myb_interaction level);

/* First post-Newtonian N-body equations (Einstein-Infeld-Hoffmann with the PPN parameters;
 * T. D. Moyer, DESCANSO Monograph 2 (2003) eq. 4-61, the equations of the JPL ephemerides,
 * W. M. Folkner et al., IPN PR 42-196 (2014) eq. 1). For a receiver i and a source j with
 * r_ij = |x_j - x_i|, the Newtonian term GM_j (x_j - x_i) / r_ij^3 is joined by
 *   GM_j (x_j - x_i) / (c^2 r_ij^3) { -2 (beta + gamma) U_i - (2 beta - 1) U_j + gamma v_i^2
 *       + (1 + gamma) v_j^2 - 2 (1 + gamma) v_i.v_j - 3/2 [(x_i - x_j).v_j / r_ij]^2
 *       + 1/2 (x_j - x_i).a_j }
 *   + GM_j / (c^2 r_ij^3) [(x_i - x_j).((2 + 2 gamma) v_i - (1 + 2 gamma) v_j)] (v_i - v_j)
 *   + (3 + 4 gamma) GM_j a_j / (2 c^2 r_ij),
 * with U_k the Newtonian potential and a_k the Newtonian acceleration of body k from its own
 * sources (the graph), v_k its velocity. Checked against the Euler-Lagrange equations of the
 * global PPN Lagrangian (C. Bernus et al., arXiv:2301.01186, app. B; tests/data/rel).
 *
 * myb_nbody_set_relativity stores the parameters and makes EIH the default level of every
 * pair (NULL: back to Newtonian). The accelerations then depend on the velocities: only
 * myb_nbody_acceleration_xv computes them (myb_nbody_acceleration fails with MYB_E_DOMAIN),
 * and the Wisdom-Holman integrators reject the system. Errors: MYB_E_DOMAIN (invalid
 * parameters). */
MYB_API MYB_NODISCARD myb_status myb_nbody_set_relativity(myb_nbody* nbody, const myb_ppn* ppn);

/* a = accelerations at positions x (3 per body) of the system nbody (a myb_nbody*: the
 * signature is myb_position_accel_function, so the system can be integrated directly by
 * RKN12, Gauss-Jackson and the symplectic compositions; t is not used by Newtonian gravity).
 * Separations must lie within 1e-150 .. 1e150 m. Errors: MYB_E_DOMAIN (two interacting bodies
 * at the same place, or a non-finite position or result); a is then unspecified. No memory
 * is allocated. */
MYB_API MYB_NODISCARD myb_status myb_nbody_acceleration(void* nbody, double t, const double* x,
                                                        double* a);

/* The same with velocities (used by the 1PN terms): the signature of myb_accel_function, for
 * IAS15 and Gauss-Jackson. With 1PN terms no allocation either; the Newtonian part and the 1PN
 * part are summed per receiver in a fixed order (Sum2), so the bits do not depend on the
 * number of threads. */
MYB_API MYB_NODISCARD myb_status myb_nbody_acceleration_xv(void* nbody, double t, const double* x,
                                                           const double* v, double* a);

/* Conserved quantities of the massive bodies, in double-word arithmetic (energy to ~1e-30
 * relative, so that integration errors of 1e-16 are resolved). With GM in place of masses
 * every quantity carries one factor G:
 *   energy           G E = sum_i GM_i |v_i|^2 / 2 - sum_{i<j} GM_i GM_j / |x_i - x_j|  (m^5 s^-4)
 *   momentum         G P = sum_i GM_i v_i                                             (m^4 s^-3)
 *   angular_momentum G L = sum_i GM_i x_i x v_i                                       (m^5 s^-3)
 *   center, center_velocity: the centre of mass and its velocity (sum GM_i > 0 needed).
 * Every pair is counted (the interaction graph is not consulted). Errors: MYB_E_DOMAIN
 * (coincident massive bodies, sum of GM zero, non-finite input). */
typedef struct myb_nbody_invariants {
  myb_dd energy;
  double momentum[3];
  double angular_momentum[3];
  double center[3];
  double center_velocity[3];
} myb_nbody_invariants;

MYB_API MYB_NODISCARD myb_status myb_nbody_invariants_of(const myb_nbody* nbody, const double* x,
                                                         const double* v,
                                                         myb_nbody_invariants* invariants);

/* First post-Newtonian conserved quantities of the massive bodies for the parameters of
 * myb_nbody_set_relativity (C. Bernus et al., arXiv:2301.01186, app. B), with GM in place of
 * masses (every quantity carries one factor G), U_A = sum_{B != A} GM_B / r_AB:
 *   energy   h = sum GM v^2/2 - 1/2 sum sum' GM_A GM_B / r_AB + 1/c^2 [ sum 3 GM v^4 / 8
 *              + (2 beta - 1)/2 sum GM_A U_A^2 + 1/2 sum sum' GM_A GM_B / r_AB
 *                ((2 gamma + 1) v_A^2 - (4 gamma + 3)/2 v_A.v_B - 1/2 (v_A.n)(v_B.n)) ]
 *   momentum P = sum GM_A v_A [1 + (v_A^2 - U_A) / (2 c^2)]
 *              - 1/(2 c^2) sum sum' GM_A GM_B / r_AB (n_AB.v_A) n_AB
 *   centre  (the relativistic barycentre, origin of the BCRS of the JPL ephemerides) and its
 *            velocity: sum GM*_A x_A / sum GM*_A and P / sum GM*_A,
 *            GM*_A = GM_A [1 + (v_A^2 - U_A) / (2 c^2)].
 * Errors: MYB_E_DOMAIN (no relativity parameters, coincident bodies, non-finite input). */
typedef struct myb_nbody_pn_invariants {
  myb_dd energy;
  double momentum[3];
  double center[3];
  double center_velocity[3];
} myb_nbody_pn_invariants;

MYB_API MYB_NODISCARD myb_status myb_nbody_pn_invariants_of(const myb_nbody* nbody,
                                                            const double* x, const double* v,
                                                            myb_nbody_pn_invariants* invariants);

/* ---- Wisdom-Holman splitting integrators ----------------------------------------------------
 * H = H_Kepler + H_interaction in canonical Jacobi coordinates (Wisdom & Holman 1991): body 0
 * is the central body (GM > 0), the massive bodies follow in order (best accuracy from the
 * inside out), test bodies are referred to the centre of mass of all massive bodies. With
 * eta_i = GM_0 + ... + GM_i, Jacobi body i >= 1 follows the Kepler problem with GM = eta_i
 * (drift: the universal-variable propagator, all conics), and the kick is the full N-body
 * acceleration less those Kepler terms (whatever the interaction graph says, so that the
 * splitting stays consistent). The Jacobi state is accumulated with compensated summation,
 * so that rounding stays below one ulp per step.
 *
 * Schemes (A = drift, B = kick, step h, every coefficient a fraction of h):
 *   MYB_WH_DKD   A(1/2) B(1) A(1/2)            (Wisdom-Holman leapfrog, = SABA1)
 *   MYB_WH_SABAn A(c1) B(d1) ... B(dn) A(c1)   B at the n Gauss-Legendre nodes of the step
 *   MYB_WH_SBABn B(d1) A(c2) ... A(c2) B(d1)   B at the n+1 Gauss-Lobatto nodes (SBAB1 = KDK)
 * With e the ratio of the perturbation to the Kepler part, the error is O(e h^(2n) + e^2 h^2)
 * for SABAn and SBABn (Laskar & Robutel 2001). For MYB_WH_DKD a symplectic corrector of odd
 * order k (3 .. 17) makes it O(e h^(k+1) + e^2 h^2) at the cost of 2 (k - 1) extra kicks per
 * reported state only: the integrator advances the uncorrected (kernel) state and the
 * corrector is applied to the states it reports (Wisdom, Holman & Touma 1996; coefficients derived in
 * dyn/wh_coefficients.h). Close encounters between bodies break the splitting: use IAS15. */

typedef enum myb_wh_scheme {
  MYB_WH_DKD = 0,
  MYB_WH_SABA2 = 1,
  MYB_WH_SABA3 = 2,
  MYB_WH_SABA4 = 3,
  MYB_WH_SBAB1 = 4,
  MYB_WH_SBAB2 = 5,
  MYB_WH_SBAB3 = 6,
  MYB_WH_SBAB4 = 7
} myb_wh_scheme;

typedef struct myb_wh_options {
  double h;             /* step (s), finite and nonzero; its sign is the forward direction */
  myb_wh_scheme scheme;
  int corrector;        /* 0, or an odd order 3 .. 17 (MYB_WH_DKD only) */
} myb_wh_options;

typedef struct myb_wh_stats {
  int64_t steps;       /* kernel steps (either direction) */
  int64_t evaluations; /* N-body acceleration evaluations, corrector and dense output included */
  int64_t drifts;      /* Kepler drifts of the whole system */
} myb_wh_stats;

typedef struct myb_wh myb_wh;

/* The system must outlive the integrator. Errors: MYB_E_DOMAIN (options, GM_0 = 0, 1PN
 * terms in the system), MYB_E_OVERFLOW, MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_wh_create(myb_nbody* system, const myb_wh_options* options,
                                               const myb_allocator* allocator, myb_wh** wh);

/* Accepts NULL. */
MYB_API void myb_wh_destroy(myb_wh* wh);

/* Starts at time t0 from inertial positions x and velocities v (3 per body of the system).
 * Errors: MYB_E_DOMAIN (non-finite input, coincident bodies, 1PN terms in the system), as for
 * the drift and the kick. */
MYB_API MYB_NODISCARD myb_status myb_wh_start(myb_wh* wh, double t0, const double* x,
                                              const double* v);

/* One step of h. The time of the k-th node is t0 + k h (computed from k, not accumulated). */
MYB_API MYB_NODISCARD myb_status myb_wh_step(myb_wh* wh);

/* Steps (forward or backward) to the node nearest to t_end; myb_wh_dense then reaches t_end
 * itself. Errors: as myb_wh_step; MYB_E_DOMAIN for a non-finite t_end or more than 2^53
 * steps. */
MYB_API MYB_NODISCARD myb_status myb_wh_integrate(myb_wh* wh, double t_end);

/* Time of the current node. */
MYB_API double myb_wh_time(const myb_wh* wh);

/* Corrected inertial state at the current node. */
MYB_API MYB_NODISCARD myb_status myb_wh_state(myb_wh* wh, double* x, double* v);

/* Inertial state at t with |t - node| <= |h|: one step of length t - node of the same scheme
 * (and corrector) from the corrected node state, on scratch memory; the integration itself is
 * not changed. Errors: MYB_E_RANGE (t too far from the node), as myb_wh_step. */
MYB_API MYB_NODISCARD myb_status myb_wh_dense(myb_wh* wh, double t, double* x, double* v);

MYB_API myb_wh_stats myb_wh_statistics(const myb_wh* wh);

#ifdef __cplusplus
}
#endif

#endif /* MYB_DYN_H_ */
