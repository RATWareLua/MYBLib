/* SPDX-License-Identifier: BSD-3-Clause */
/* The numerical environment the library relies on, checked at run time. */
#include <fenv.h>

#include "core/fpguard.h"
#include "mtest.h"

double myb_test_probe_mul_add(double a, double b, double c);

/* a*b = 1 + 2^-29 + 2^-60 exactly; RN(a*b) drops 2^-60, so the unfused result is 0 and the
 * fused one is 2^-60. */
static void test_no_contraction(void) {
  if (!mt_cpu_has_fma()) {
    printf("no FMA on this CPU: the contraction probe is vacuous here\n");
  }
  const double a = 1.0 + 0x1p-30;
  const double c = -(1.0 + 0x1p-29);
  MT_CHECK_SAME(myb_test_probe_mul_add(a, a, c), 0.0);
}

static void test_rounding_mode(void) {
  MT_CHECK(fegetround() == FE_TONEAREST);
}

/* Subnormals must not be flushed to zero (no FTZ/DAZ). */
static void test_gradual_underflow(void) {
  volatile double tiny = 0x1p-1022;
  const double half = tiny * 0.5;
  MT_CHECK_SAME(half, 0x1p-1023);
  volatile double smallest = 0x1p-1074;
  MT_CHECK(smallest * 1.0 != 0.0);
}

int main(void) {
  test_no_contraction();
  test_rounding_mode();
  test_gradual_underflow();
  return mt_finish("test_fpenv");
}
