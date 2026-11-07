/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/num.h — error-free transformations, double-word ("double-double") arithmetic and
 * compensated summation.
 *
 * Notation: u = 2^-53 is the unit roundoff of binary64, RN() is rounding to nearest even.
 * A double-word number x = (hi, lo) represents hi + lo exactly, with hi = RN(hi + lo);
 * its precision is about 106 bits.
 *
 * Sources:
 *   [D71]  T. J. Dekker, "A floating-point technique for extending the available
 *          precision", Numer. Math. 18 (1971) 224-242.
 *   [K97]  D. E. Knuth, TAOCP vol. 2, 3rd ed. (1997), §4.2.2 (TwoSum).
 *   [JMP17] M. Joldes, J.-M. Muller, V. Popescu, "Tight and rigorous error bounds for basic
 *          building blocks of double-word arithmetic", ACM TOMS 44(2), 2017.
 *   [MR22] J.-M. Muller, L. Rideau, "Formalization of double-word arithmetic, and comments
 *          on [JMP17]", ACM TOMS 48(1), 2022.
 *   [ORO05] T. Ogita, S. M. Rump, S. Oishi, "Accurate sum and dot product",
 *          SIAM J. Sci. Comput. 26(6), 2005.
 *
 * The error bounds quoted below are the proven ones from the papers; tests/mpfr measures the
 * actual maximum error against MPFR and fails if a bound is exceeded. */
#ifndef MYB_NUM_H_
#define MYB_NUM_H_

#include <stddef.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Error-free transformations ------------------------------------------------------ */

/* s = RN(a + b) and *err = (a + b) - s exactly (TwoSum, [K97]); no condition on a, b
 * except that no overflow occurs. */
MYB_API double myb_two_sum(double a, double b, double* err);

/* p = RN(a * b) and *err = a*b - p exactly (Dekker's product with Veltkamp splitting,
 * [D71]). Exact provided no overflow occurs and e_a + e_b >= -970 (no underflow of the
 * error term), where e_x is the binary exponent of x. Results are identical to the FMA
 * formulation err = fma(a, b, -p) whenever both are exact. */
MYB_API double myb_two_prod(double a, double b, double* err);

/* ---- Double-word arithmetic ---------------------------------------------------------- */

typedef struct myb_dd {
  double hi;
  double lo;
} myb_dd;

/* Exact conversions. */
MYB_API myb_dd myb_dd_from_double(double x);
MYB_API double myb_dd_to_double(myb_dd x); /* RN(hi + lo) == hi */

/* Algorithms (names as in [JMP17]) and their proven relative error bounds, u = 2^-53:
 *   add, sub   AccurateDWPlusDW   3u^2 + O(u^3)   ([JMP17], proof revisited in [MR22])
 *   add_d      DWPlusFP           2u^2
 *   mul        DWTimesDW1         7u^2
 *   mul_d      DWTimesFP1         1.5u^2 + 4u^3
 *   div        DWDivDW2           15u^2 + 56u^3
 *   div_d      DWDivFP1           3.5u^2
 *   sqrt       one Karp-Markstein correction with an exact residual: no published bound for
 *              this variant; tests/mpfr requires < 4u^2 over its sample.
 * The FMA-free variants are used on purpose: the x86-64-v2 baseline has no hardware FMA.
 * Arguments must be finite and results must not overflow or underflow. */
MYB_API myb_dd myb_dd_add(myb_dd x, myb_dd y);
MYB_API myb_dd myb_dd_sub(myb_dd x, myb_dd y);
MYB_API myb_dd myb_dd_add_d(myb_dd x, double y);
MYB_API myb_dd myb_dd_mul(myb_dd x, myb_dd y);
MYB_API myb_dd myb_dd_mul_d(myb_dd x, double y);
MYB_API myb_dd myb_dd_div(myb_dd x, myb_dd y);
MYB_API myb_dd myb_dd_div_d(myb_dd x, double y);
MYB_API myb_dd myb_dd_sqrt(myb_dd x); /* x >= 0 */
MYB_API myb_dd myb_dd_neg(myb_dd x);

/* ---- Compensated summation ----------------------------------------------------------- */

/* Running sum in the style of Sum2 [ORO05]: the result is as accurate as if the
 * sum were computed in twice the working precision and then rounded:
 *   |result - S| <= u|S| + gamma_{n-1}^2 * sum|x_i|,   gamma_k = k*u / (1 - k*u).
 * Order of additions is the caller's order: the result is deterministic. */
typedef struct myb_sum_acc {
  double sum;
  double comp;
} myb_sum_acc;

/* Zero-initialised myb_sum_acc ({0}) is an empty sum. */
MYB_API void myb_sum_acc_add(myb_sum_acc* acc, double x);
MYB_API double myb_sum_acc_value(const myb_sum_acc* acc);

/* Sum2 / Dot2 over arrays [ORO05]; n == 0 gives 0.
 * Dot2: |result - x.y| <= u|x.y| + gamma_n^2 * |x|.|y|. */
MYB_API double myb_sum2(const double* x, size_t n);
MYB_API double myb_dot2(const double* x, const double* y, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* MYB_NUM_H_ */
