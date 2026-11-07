/* SPDX-License-Identifier: BSD-3-Clause */
/* Safeguarded Newton iteration. The safeguard is the classical one of "rtsafe" (W. H. Press
 * et al., Numerical Recipes, §9.4; algorithm re-implemented here): a Newton step is taken only
 * if it stays in the bracket and is at most half of the step before the previous one;
 * otherwise the bracket is bisected. The steps therefore shrink at least geometrically,
 * whatever the shape of the function (inflections included). */
#include "num/root.h"

#include <math.h> /* fabs */

#include "core/base.h"

myb_status myb_internal_newton_bracketed(const myb_internal_root_function* function, double lo,
                                         double hi, double x, int max_iterations, double* root) {
  MYB_CHECK(function != NULL && function->eval != NULL && root != NULL);
  MYB_CHECK(lo <= hi);
  if (!(x >= lo && x <= hi)) {
    x = 0.5 * (lo + hi);
  }
  double step_before_last = hi - lo;
  double last_step = step_before_last;
  int converged = 0;
  /* bound: max_iterations; the step halves at least every second iteration */
  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    double f = 0.0;
    double df = 0.0;
    function->eval(function->context, x, &f, &df);
    if (f == 0.0) {
      *root = x;
      return MYB_OK;
    }
    if (f < 0.0) {
      lo = x;
    } else {
      hi = x;
    }
    const double next = x - f / df;
    if (next == x) {
      *root = x;
      return MYB_OK;
    }
    const int inside = next >= lo && next <= hi;
    const int fast = fabs(2.0 * f) <= fabs(step_before_last * df);
    if (inside && fast) {
      if (converged) {
        *root = next; /* one Newton correction after reaching rounding level */
        return MYB_OK;
      }
      step_before_last = last_step;
      last_step = next - x;
      converged = fabs(last_step) <= 0x1p-50 * fabs(next);
      x = next;
    } else {
      step_before_last = last_step;
      last_step = 0.5 * (hi - lo);
      converged = 0;
      x = lo + last_step;
    }
    if (hi - lo <= 0x1p-53 * (fabs(hi) + fabs(lo))) {
      *root = x;
      return MYB_OK;
    }
  }
  return MYB_E_NOCONV;
}
