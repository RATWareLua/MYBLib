/* SPDX-License-Identifier: BSD-3-Clause */
#include "myb/cmath.h"

#include "cmath/fma_soft.h"
#include "core/base.h"
#include "myb_coremath.h"

double myb_fma(double a, double b, double c) {
  return myb_internal_fma(a, b, c);
}

double myb_sin(double x) {
  return myb_cm_sin(x);
}

double myb_cos(double x) {
  return myb_cm_cos(x);
}
