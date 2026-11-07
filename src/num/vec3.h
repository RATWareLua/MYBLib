/* SPDX-License-Identifier: BSD-3-Clause */
/* num/vec3.h — 3-vector helpers for internal use (arrays of 3 doubles). The operation order
 * is fixed, so results are reproducible. Magnitudes up to 1e150 do not overflow in norms. */
#ifndef MYB_SRC_NUM_VEC3_H_
#define MYB_SRC_NUM_VEC3_H_

#include <math.h> /* sqrt: correctly rounded by IEEE-754 */

static inline double myb_vec3_dot(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline double myb_vec3_norm(const double a[3]) {
  return sqrt(myb_vec3_dot(a, a));
}

static inline void myb_vec3_cross(const double a[3], const double b[3], double out[3]) {
  const double x = a[1] * b[2] - a[2] * b[1];
  const double y = a[2] * b[0] - a[0] * b[2];
  const double z = a[0] * b[1] - a[1] * b[0];
  out[0] = x;
  out[1] = y;
  out[2] = z;
}

/* out = a * s + b * t */
static inline void myb_vec3_combine(const double a[3], double s, const double b[3], double t,
                                    double out[3]) {
  for (int i = 0; i < 3; ++i) {
    out[i] = a[i] * s + b[i] * t;
  }
}

static inline int myb_vec3_finite(const double a[3]) {
  return (a[0] - a[0] == 0.0) && (a[1] - a[1] == 0.0) && (a[2] - a[2] == 0.0);
}

#endif /* MYB_SRC_NUM_VEC3_H_ */
