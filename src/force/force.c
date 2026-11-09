/* SPDX-License-Identifier: BSD-3-Clause */
/* Point-mass perturbations, drag, empirical accelerations. */
#include <math.h>

#include "core/base.h"
#include "myb/force.h"
#include "num/vec3.h"

void myb_force_third_body(double gm, const double s[3], const double r[3], double a[3]) {
  MYB_CHECK(s != NULL && r != NULL && a != NULL);
  /* (s - r)/|d|^3 - s/|s|^3 = -(r + f(q) s)/|d|^3, q = r.(r - 2 s)/s.s,
   * f(q) = (1 + q)^(3/2) - 1 = q (3 + 3q + q^2) / (1 + (1 + q)^(3/2)) */
  const double ss = myb_vec3_dot(s, s);
  const double q = (r[0] * (r[0] - 2.0 * s[0]) + r[1] * (r[1] - 2.0 * s[1]) +
                    r[2] * (r[2] - 2.0 * s[2])) / ss;
  const double p = (1.0 + q) * sqrt(1.0 + q);
  const double f = q * (3.0 + q * (3.0 + q)) / (1.0 + p);
  const double d[3] = {s[0] - r[0], s[1] - r[1], s[2] - r[2]};
  const double dd = myb_vec3_norm(d);
  const double k = -gm / (dd * dd * dd);
  for (int i = 0; i < 3; ++i) {
    a[i] += k * (r[i] + f * s[i]);
  }
}

void myb_force_drag(double density, double cd_area_over_mass, const double v[3], double a[3]) {
  MYB_CHECK(v != NULL && a != NULL);
  const double k = -0.5 * density * cd_area_over_mass * myb_vec3_norm(v);
  for (int i = 0; i < 3; ++i) {
    a[i] += k * v[i];
  }
}

myb_status myb_force_empirical_rtn(const double r[3], const double v[3], const double a_rtn[3],
                                   double a[3]) {
  MYB_CHECK(r != NULL && v != NULL && a_rtn != NULL && a != NULL);
  double n[3];
  myb_vec3_cross(r, v, n);
  const double rr = myb_vec3_norm(r);
  const double nn = myb_vec3_norm(n);
  if (!(rr > 0.0) || !(nn > 0.0) || !isfinite(rr) || !isfinite(nn)) {
    return MYB_E_DOMAIN;
  }
  double e_r[3];
  double e_n[3];
  for (int i = 0; i < 3; ++i) {
    e_r[i] = r[i] / rr;
    e_n[i] = n[i] / nn;
  }
  double e_t[3];
  myb_vec3_cross(e_n, e_r, e_t);
  for (int i = 0; i < 3; ++i) {
    a[i] += a_rtn[0] * e_r[i] + a_rtn[1] * e_t[i] + a_rtn[2] * e_n[i];
  }
  return MYB_OK;
}
