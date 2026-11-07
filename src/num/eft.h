/* SPDX-License-Identifier: BSD-3-Clause */
/* num/eft.h — error-free transformations for internal use (inline, library flags only).
 *
 * All of them are exact under the stated conditions and use only +, -, * of binary64, so
 * they give the same bits on every IEEE-754 platform with FLT_EVAL_METHOD == 0 and without
 * contraction (both enforced by core/fpguard.h and cmake/MybFlags.cmake).
 * Sources: Knuth TAOCP vol. 2 §4.2.2 (TwoSum); Dekker 1971 (Fast2Sum, Mul12, splitting);
 * Muller et al., Handbook of Floating-Point Arithmetic, 2nd ed. (2018), ch. 4. */
#ifndef MYB_SRC_NUM_EFT_H_
#define MYB_SRC_NUM_EFT_H_

/* s = RN(a + b), *err = a + b - s exactly. Branch-free, no condition on a and b. */
static inline double myb_eft_two_sum(double a, double b, double* err) {
  const double s = a + b;
  const double a_rounded = s - b;
  const double b_rounded = s - a_rounded;
  const double a_error = a - a_rounded;
  const double b_error = b - b_rounded;
  *err = a_error + b_error;
  return s;
}

/* s = RN(a + b), *err = a + b - s exactly, provided exponent(a) >= exponent(b) (true when
 * |a| >= |b|) or a == 0. */
static inline double myb_eft_fast_two_sum(double a, double b, double* err) {
  const double s = a + b;
  const double b_rounded = s - a;
  *err = b - b_rounded;
  return s;
}

/* Veltkamp splitting: a = *hi + *lo exactly, each part fitting in 26 bits, so that products
 * of parts are exact. For |a| > 2^996 the splitter product a * (2^27 + 1) would overflow, so
 * such a is scaled by 2^-28 first (exact) and the parts are scaled back (exact). */
static inline void myb_eft_split(double a, double* hi, double* lo) {
  constexpr double splitter = 0x1p27 + 1.0;
  constexpr double big = 0x1p996;
  double scale_back = 1.0;
  if (a > big || a < -big) {
    a *= 0x1p-28;
    scale_back = 0x1p28;
  }
  const double c = splitter * a;
  const double h = c - (c - a);
  *hi = h * scale_back;
  *lo = (a - h) * scale_back;
}

/* p = RN(a * b), *err = a * b - p exactly (Dekker's Mul12), provided a * b does not overflow
 * and e_a + e_b >= -970 (no underflow of the error term). */
static inline double myb_eft_two_prod(double a, double b, double* err) {
  const double p = a * b;
  double a_hi;
  double a_lo;
  double b_hi;
  double b_lo;
  myb_eft_split(a, &a_hi, &a_lo);
  myb_eft_split(b, &b_hi, &b_lo);
  *err = (((a_hi * b_hi - p) + a_hi * b_lo) + a_lo * b_hi) + a_lo * b_lo;
  return p;
}

#endif /* MYB_SRC_NUM_EFT_H_ */
