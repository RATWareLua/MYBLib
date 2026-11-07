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

double myb_tan(double x) {
  return myb_cm_tan(x);
}

double myb_asin(double x) {
  return myb_cm_asin(x);
}

double myb_acos(double x) {
  return myb_cm_acos(x);
}

double myb_atan(double x) {
  return myb_cm_atan(x);
}

double myb_exp(double x) {
  return myb_cm_exp(x);
}

double myb_expm1(double x) {
  return myb_cm_expm1(x);
}

double myb_log(double x) {
  return myb_cm_log(x);
}

double myb_log1p(double x) {
  return myb_cm_log1p(x);
}

double myb_sinh(double x) {
  return myb_cm_sinh(x);
}

double myb_cosh(double x) {
  return myb_cm_cosh(x);
}

double myb_tanh(double x) {
  return myb_cm_tanh(x);
}

double myb_asinh(double x) {
  return myb_cm_asinh(x);
}

double myb_acosh(double x) {
  return myb_cm_acosh(x);
}

double myb_atanh(double x) {
  return myb_cm_atanh(x);
}

double myb_cbrt(double x) {
  return myb_cm_cbrt(x);
}

double myb_atan2(double y, double x) {
  return myb_cm_atan2(y, x);
}

double myb_hypot(double x, double y) {
  return myb_cm_hypot(x, y);
}
