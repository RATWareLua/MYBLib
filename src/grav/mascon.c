/* SPDX-License-Identifier: BSD-3-Clause */
/* Mass concentrations: a field as a sum of point masses (lunar and asteroid gravity models,
 * corrections to a spherical-harmonic field). Sums are compensated (Sum2) and taken in the
 * order of the array, so the result does not depend on how the caller splits the work. */
#include <math.h>

#include "core/base.h"
#include "myb/grav.h"
#include "myb/num.h"

myb_status myb_mascons_eval(const double* gm, const double* positions, size_t count,
                            const double r[3], double* potential, double a[3]) {
  MYB_CHECK(r != NULL && (count == 0 || (gm != NULL && positions != NULL)));
  if (!(r[0] - r[0] == 0.0) || !(r[1] - r[1] == 0.0) || !(r[2] - r[2] == 0.0)) {
    return MYB_E_DOMAIN;
  }
  myb_sum_acc v = {0};
  myb_sum_acc g[3] = {{0}};
  for (size_t k = 0; k < count; ++k) {
    const double* p = positions + 3 * k;
    const double d[3] = {p[0] - r[0], p[1] - r[1], p[2] - r[2]};
    const double dd = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (!(dd > 0.0)) {
      return MYB_E_DOMAIN; /* at a mascon */
    }
    const double k1 = gm[k] / dd;
    myb_sum_acc_add(&v, k1);
    const double k3 = k1 / (dd * dd);
    for (int i = 0; i < 3; ++i) {
      myb_sum_acc_add(&g[i], k3 * d[i]);
    }
  }
  if (potential != NULL) {
    *potential = myb_sum_acc_value(&v);
  }
  if (a != NULL) {
    for (int i = 0; i < 3; ++i) {
      a[i] = myb_sum_acc_value(&g[i]);
    }
  }
  return MYB_OK;
}
