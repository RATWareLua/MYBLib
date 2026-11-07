/* SPDX-License-Identifier: BSD-3-Clause */
/* Double-word arithmetic. Algorithm names and error bounds: M. Joldes, J.-M. Muller,
 * V. Popescu, ACM TOMS 44(2), 2017 [JMP17]; see include/myb/num.h for the bounds. */
#include <math.h> /* sqrt only: correctly rounded by IEEE-754, a single instruction here */

#include "core/base.h"
#include "myb/num.h"
#include "num/eft.h"

double myb_two_sum(double a, double b, double* err) {
  MYB_CHECK(err != NULL);
  return myb_eft_two_sum(a, b, err);
}

double myb_two_prod(double a, double b, double* err) {
  MYB_CHECK(err != NULL);
  return myb_eft_two_prod(a, b, err);
}

myb_dd myb_dd_from_double(double x) {
  return (myb_dd){.hi = x, .lo = 0.0};
}

double myb_dd_to_double(myb_dd x) {
  return x.hi + x.lo;
}

myb_dd myb_dd_neg(myb_dd x) {
  return (myb_dd){.hi = -x.hi, .lo = -x.lo};
}

/* AccurateDWPlusDW [JMP17]. */
myb_dd myb_dd_add(myb_dd x, myb_dd y) {
  double s_lo;
  const double s_hi = myb_eft_two_sum(x.hi, y.hi, &s_lo);
  double t_lo;
  const double t_hi = myb_eft_two_sum(x.lo, y.lo, &t_lo);
  const double c = s_lo + t_hi;
  double v_lo;
  const double v_hi = myb_eft_fast_two_sum(s_hi, c, &v_lo);
  const double w = t_lo + v_lo;
  double z_lo;
  const double z_hi = myb_eft_fast_two_sum(v_hi, w, &z_lo);
  return (myb_dd){.hi = z_hi, .lo = z_lo};
}

myb_dd myb_dd_sub(myb_dd x, myb_dd y) {
  return myb_dd_add(x, myb_dd_neg(y));
}

/* DWPlusFP [JMP17]. */
myb_dd myb_dd_add_d(myb_dd x, double y) {
  double s_lo;
  const double s_hi = myb_eft_two_sum(x.hi, y, &s_lo);
  const double v = x.lo + s_lo;
  double z_lo;
  const double z_hi = myb_eft_fast_two_sum(s_hi, v, &z_lo);
  return (myb_dd){.hi = z_hi, .lo = z_lo};
}

/* DWTimesDW1 [JMP17]: the FMA-free variant (x86-64-v2 has no hardware FMA). */
myb_dd myb_dd_mul(myb_dd x, myb_dd y) {
  double c_lo1;
  const double c_hi = myb_eft_two_prod(x.hi, y.hi, &c_lo1);
  const double t_lo1 = x.hi * y.lo;
  const double t_lo2 = x.lo * y.hi;
  const double c_lo2 = t_lo1 + t_lo2;
  const double c_lo3 = c_lo1 + c_lo2;
  double z_lo;
  const double z_hi = myb_eft_fast_two_sum(c_hi, c_lo3, &z_lo);
  return (myb_dd){.hi = z_hi, .lo = z_lo};
}

/* DWTimesFP1 [JMP17]. */
static myb_dd dd_times_fp1(myb_dd x, double y) {
  double c_lo1;
  const double c_hi = myb_eft_two_prod(x.hi, y, &c_lo1);
  const double c_lo2 = x.lo * y;
  double t_lo1;
  const double t_hi = myb_eft_fast_two_sum(c_hi, c_lo2, &t_lo1);
  const double t_lo2 = t_lo1 + c_lo1;
  double z_lo;
  const double z_hi = myb_eft_fast_two_sum(t_hi, t_lo2, &z_lo);
  return (myb_dd){.hi = z_hi, .lo = z_lo};
}

myb_dd myb_dd_mul_d(myb_dd x, double y) {
  return dd_times_fp1(x, y);
}

/* DWDivFP1 [JMP17]. */
myb_dd myb_dd_div_d(myb_dd x, double y) {
  const double t_hi = x.hi / y;
  double p_lo;
  const double p_hi = myb_eft_two_prod(t_hi, y, &p_lo);
  const double d_hi = x.hi - p_hi;
  const double d_t = d_hi - p_lo;
  const double d = d_t + x.lo;
  const double t_lo = d / y;
  double z_lo;
  const double z_hi = myb_eft_fast_two_sum(t_hi, t_lo, &z_lo);
  return (myb_dd){.hi = z_hi, .lo = z_lo};
}

/* DWDivDW2 [JMP17]. x.hi - r.hi is exact (Sterbenz lemma). */
myb_dd myb_dd_div(myb_dd x, myb_dd y) {
  const double t_hi = x.hi / y.hi;
  const myb_dd r = dd_times_fp1(y, t_hi);
  const double p_hi = x.hi - r.hi;
  const double d_lo = x.lo - r.lo;
  const double d = p_hi + d_lo;
  const double t_lo = d / y.hi;
  double z_lo;
  const double z_hi = myb_eft_fast_two_sum(t_hi, t_lo, &z_lo);
  return (myb_dd){.hi = z_hi, .lo = z_lo};
}

/* One Karp-Markstein correction of s = RN(sqrt(x.hi)): the residual x - s^2 is formed with
 * the exact product s*s = p_hi + p_lo; x.hi - p_hi is exact by the Sterbenz lemma. */
myb_dd myb_dd_sqrt(myb_dd x) {
  if (x.hi == 0.0) {
    return (myb_dd){.hi = x.hi, .lo = 0.0};
  }
  const double s = sqrt(x.hi);
  double p_lo;
  const double p_hi = myb_eft_two_prod(s, s, &p_lo);
  const double d = x.hi - p_hi;
  const double r = (d - p_lo) + x.lo;
  const double e = r / (2.0 * s);
  double z_lo;
  const double z_hi = myb_eft_fast_two_sum(s, e, &z_lo);
  return (myb_dd){.hi = z_hi, .lo = z_lo};
}
