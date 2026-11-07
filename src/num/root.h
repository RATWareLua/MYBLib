/* SPDX-License-Identifier: BSD-3-Clause */
/* num/root.h — Newton's method kept inside a bracket, for increasing functions.
 *
 * Contract: f is increasing on [lo, hi], f(lo) <= 0 <= f(hi) (NaN counts as "beyond the
 * root"), df > 0 is its derivative. Each Newton step that leaves the bracket is replaced by
 * bisection, so the iteration cannot diverge. Termination:
 *   - f(x) == 0, or the Newton step is below the resolution of x (next == x);
 *   - once a step is below 2^-50 |x| (quadratic convergence reached), one more Newton
 *     correction is taken if it stays in the bracket, and returned;
 *   - the bracket shrinks to the resolution of the doubles.
 * The callback must not call back into the solver (no recursion). */
#ifndef MYB_SRC_NUM_ROOT_H_
#define MYB_SRC_NUM_ROOT_H_

#include "myb/core.h"

typedef struct myb_internal_root_function {
  void (*eval)(const void* context, double x, double* f, double* df);
  const void* context;
} myb_internal_root_function;

/* MYB_E_NOCONV if max_iterations are exhausted (start x is clamped into the bracket). */
myb_status myb_internal_newton_bracketed(const myb_internal_root_function* function, double lo,
                                         double hi, double x, int max_iterations, double* root);

#endif /* MYB_SRC_NUM_ROOT_H_ */
