/* SPDX-License-Identifier: BSD-3-Clause */
/* num/dual.h — a value together with its derivatives with respect to six variables (the
 * position and the velocity of an object): the chain rule carried along with the arithmetic,
 * so that the Jacobian of a formula is the formula itself written on these numbers. Internal;
 * the order of the operations is fixed, so the results are reproducible. */
#ifndef MYB_SRC_NUM_DUAL_H_
#define MYB_SRC_NUM_DUAL_H_

#include <math.h> /* sqrt */

#include "myb/cmath.h"

enum { MYB_DUAL_VARIABLES = 6 };

typedef struct myb_dual {
  double v;
  double d[MYB_DUAL_VARIABLES];
} myb_dual;

static inline myb_dual myb_dual_constant(double v) {
  myb_dual r = {.v = v};
  return r;
}

/* The variable number index (0 .. 5) with the value v. */
static inline myb_dual myb_dual_variable(double v, int index) {
  myb_dual r = {.v = v};
  r.d[index] = 1.0;
  return r;
}

static inline myb_dual myb_dual_add(myb_dual a, myb_dual b) {
  myb_dual r;
  r.v = a.v + b.v;
  for (int i = 0; i < MYB_DUAL_VARIABLES; ++i) {
    r.d[i] = a.d[i] + b.d[i];
  }
  return r;
}

static inline myb_dual myb_dual_sub(myb_dual a, myb_dual b) {
  myb_dual r;
  r.v = a.v - b.v;
  for (int i = 0; i < MYB_DUAL_VARIABLES; ++i) {
    r.d[i] = a.d[i] - b.d[i];
  }
  return r;
}

static inline myb_dual myb_dual_scale(myb_dual a, double s) {
  myb_dual r;
  r.v = a.v * s;
  for (int i = 0; i < MYB_DUAL_VARIABLES; ++i) {
    r.d[i] = a.d[i] * s;
  }
  return r;
}

static inline myb_dual myb_dual_mul(myb_dual a, myb_dual b) {
  myb_dual r;
  r.v = a.v * b.v;
  for (int i = 0; i < MYB_DUAL_VARIABLES; ++i) {
    r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  }
  return r;
}

static inline myb_dual myb_dual_div(myb_dual a, myb_dual b) {
  myb_dual r;
  r.v = a.v / b.v;
  for (int i = 0; i < MYB_DUAL_VARIABLES; ++i) {
    r.d[i] = (a.d[i] - r.v * b.d[i]) / b.v;
  }
  return r;
}

static inline myb_dual myb_dual_sqrt(myb_dual a) {
  myb_dual r;
  r.v = sqrt(a.v);
  for (int i = 0; i < MYB_DUAL_VARIABLES; ++i) {
    r.d[i] = 0.5 * a.d[i] / r.v;
  }
  return r;
}

static inline myb_dual myb_dual_exp(myb_dual a) {
  myb_dual r;
  r.v = myb_exp(a.v);
  for (int i = 0; i < MYB_DUAL_VARIABLES; ++i) {
    r.d[i] = r.v * a.d[i];
  }
  return r;
}

/* a^p for a > 0 and a constant p. */
static inline myb_dual myb_dual_power(myb_dual a, double p) {
  myb_dual r;
  r.v = myb_exp(p * myb_log(a.v));
  const double slope = p * r.v / a.v;
  for (int i = 0; i < MYB_DUAL_VARIABLES; ++i) {
    r.d[i] = slope * a.d[i];
  }
  return r;
}

/* ---- Vectors of three --------------------------------------------------------------------- */

static inline myb_dual myb_dual_dot(const myb_dual a[3], const myb_dual b[3]) {
  return myb_dual_add(myb_dual_add(myb_dual_mul(a[0], b[0]), myb_dual_mul(a[1], b[1])),
                      myb_dual_mul(a[2], b[2]));
}

static inline void myb_dual_cross(const myb_dual a[3], const myb_dual b[3], myb_dual out[3]) {
  const myb_dual x = myb_dual_sub(myb_dual_mul(a[1], b[2]), myb_dual_mul(a[2], b[1]));
  const myb_dual y = myb_dual_sub(myb_dual_mul(a[2], b[0]), myb_dual_mul(a[0], b[2]));
  const myb_dual z = myb_dual_sub(myb_dual_mul(a[0], b[1]), myb_dual_mul(a[1], b[0]));
  out[0] = x;
  out[1] = y;
  out[2] = z;
}

static inline myb_dual myb_dual_norm(const myb_dual a[3]) {
  return myb_dual_sqrt(myb_dual_dot(a, a));
}

/* The position x - s and the velocity v - s' as variables 0 .. 2 and 3 .. 5. */
static inline void myb_dual_state(const double x[3], const double v[3], const double* s,
                                  myb_dual r[3], myb_dual w[3]) {
  for (int i = 0; i < 3; ++i) {
    r[i] = myb_dual_variable(x[i] - (s != NULL ? s[i] : 0.0), i);
    w[i] = myb_dual_variable(v[i] - (s != NULL ? s[3 + i] : 0.0), 3 + i);
  }
}

/* dadx[3 i + j] and dadv[3 i + j] from the three components of an acceleration. */
static inline void myb_dual_jacobian(const myb_dual a[3], double dadx[9], double dadv[9]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      dadx[3 * i + j] = a[i].d[j];
      dadv[3 * i + j] = a[i].d[3 + j];
    }
  }
}

#endif /* MYB_SRC_NUM_DUAL_H_ */
