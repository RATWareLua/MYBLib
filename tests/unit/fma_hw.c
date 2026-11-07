/* SPDX-License-Identifier: BSD-3-Clause */
/* Hardware FMA (compiled with -mfma): the reference for myb_fma. Called only after
 * mt_cpu_has_fma() confirmed that the instruction is available. */
#include "core/fpguard.h"

double myb_test_hw_fma(double a, double b, double c);

double myb_test_hw_fma(double a, double b, double c) {
#if defined(__FMA__)
  return __builtin_fma(a, b, c);
#else
  (void)a, (void)b, (void)c;
  return 0.0; /* never called: the test skips the comparison without FMA */
#endif
}
