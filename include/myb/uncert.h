/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/uncert.h — uncertainty of a trajectory: variational equations (state transition matrix
 * and sensitivities to parameters), propagation of covariances, Monte Carlo with a
 * deterministic random number generator.
 *
 * Variational equations: for x'' = a(t, x, v, p) the columns of Phi = d(x, v)/d(x0, v0) and of
 * S = d(x, v)/dp obey dx'' = A_x dx + A_v dv (+ A_p for S), integrated with the state by IAS15
 * (one system of 3 (7 + n_p) coordinates); A from the force model (myb/model.h: analytic
 * Jacobians where the terms give them, central differences otherwise). */
#ifndef MYB_UNCERT_H_
#define MYB_UNCERT_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"
#include "myb/model.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Variational equations ------------------------------------------------------------------ */

typedef struct myb_variational myb_variational;

typedef struct myb_variational_options {
  double epsilon; /* IAS15 precision parameter; 0 = 1e-10 */
  double h_max;   /* 0 = unlimited */
  int64_t max_steps;
} myb_variational_options;

/* The model must outlive the object; its parameters at creation time are those differentiated.
 * Errors: MYB_E_DOMAIN, MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_variational_create(myb_force_model* model,
                                                        const myb_variational_options* options,
                                                        const myb_allocator* allocator,
                                                        myb_variational** variational);
MYB_API void myb_variational_destroy(myb_variational* variational); /* accepts NULL */

/* Starts at t0 from y0 = (x0, v0) with Phi = I and S = 0. */
MYB_API MYB_NODISCARD myb_status myb_variational_start(myb_variational* variational, double t0,
                                                       const double y0[6]);
/* One step towards t_limit, or up to t_end; errors as IAS15 and the model. */
MYB_API MYB_NODISCARD myb_status myb_variational_step(myb_variational* variational,
                                                      double t_limit);
MYB_API MYB_NODISCARD myb_status myb_variational_integrate(myb_variational* variational,
                                                           double t_end);
MYB_API double myb_variational_time(const myb_variational* variational);

/* The current state y (6), Phi (6 x 6, row-major: d y_i / d y0_j) and S (6 x n_p); any output
 * may be NULL. */
MYB_API void myb_variational_get(const myb_variational* variational, double y[6], double phi[36],
                                 double* sensitivity);
/* The same at any t within the last step (dense output). MYB_E_RANGE outside it. */
MYB_API MYB_NODISCARD myb_status myb_variational_dense(myb_variational* variational, double t,
                                                       double y[6], double phi[36],
                                                       double* sensitivity);

/* ---- Covariances ------------------------------------------------------------------------------ */

/* out = M P M^T for n x n matrices (out may not alias). */
MYB_API void myb_covariance_transform(size_t n, const double* m, const double* p, double* out);

/* A 6 x 6 state covariance in the radial, transverse, normal axes of (r, v) (rotation of both
 * the position and the velocity blocks; the rotation rate of the axes is not included).
 * MYB_E_DOMAIN: r and v parallel. */
MYB_API MYB_NODISCARD myb_status myb_covariance_rtn(const double r[3], const double v[3],
                                                    const double p[36], double out[36]);

/* ---- Deterministic random numbers --------------------------------------------------------------
 * xoshiro256** (D. Blackman, S. Vigna, ACM TOMS 47, 36, 2021), seeded by SplitMix64; normal
 * deviates by the Box-Muller transform with myb's correctly rounded log, sin, cos: the same
 * sequence on every platform. */

typedef struct myb_rng {
  uint64_t s[4];
  double spare; /* the second deviate of Box-Muller */
  int has_spare;
} myb_rng;

MYB_API void myb_rng_seed(myb_rng* rng, uint64_t seed);
MYB_API uint64_t myb_rng_next(myb_rng* rng);
MYB_API double myb_rng_uniform(myb_rng* rng); /* [0, 1), 53 random bits */
MYB_API double myb_rng_normal(myb_rng* rng);  /* standard normal */

/* A sample of N(mean, L L^T): out = mean + L z, with L lower triangular (myb_cholesky). */
MYB_API void myb_rng_gaussian(myb_rng* rng, size_t n, const double* mean, const double* l,
                              double* out);

/* Mean and covariance (divisor count - 1) of count samples of dimension n stored one after the
 * other; compensated sums. MYB_E_DOMAIN: count < 2. */
MYB_API MYB_NODISCARD myb_status myb_sample_statistics(size_t count, size_t n,
                                                       const double* samples, double* mean,
                                                       double* covariance);

#ifdef __cplusplus
}
#endif

#endif /* MYB_UNCERT_H_ */
