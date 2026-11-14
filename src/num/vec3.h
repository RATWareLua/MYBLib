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

/* |a| for any finite a: the plain norm (the same bits as myb_vec3_norm) where its squares
 * neither overflow nor lose digits, otherwise scaled by the largest component — infinite only
 * when |a| itself is beyond the largest double; NaN for a NaN component. */
static inline double myb_vec3_norm_scaled(const double a[3]) {
  const double plain = myb_vec3_norm(a);
  if (plain >= 0x1p-511 && plain <= 0x1p+511) {
    return plain;
  }
  const double m = fmax(fmax(fabs(a[0]), fabs(a[1])), fabs(a[2]));
  if (isnan(a[0]) || isnan(a[1]) || isnan(a[2])) {
    return (double)NAN;
  }
  if (m == 0.0 || !isfinite(m)) {
    return m;
  }
  const double x = a[0] / m;
  const double y = a[1] / m;
  const double z = a[2] / m;
  return m * sqrt(x * x + y * y + z * z);
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
