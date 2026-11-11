/* SPDX-License-Identifier: BSD-3-Clause */
/* Mean elements by averaging (myb/kepler.h). */
#include <math.h> /* sqrt, floor, ceil */
#include <stddef.h>

#include "core/check.h"
#include "myb/cmath.h"
#include "myb/kepler.h"

static const double two_pi = 6.283185307179586476925286766559;

myb_status myb_mean_elements_average(double mu, const double* times, const double* states,
                                     size_t count, myb_mean_elements* mean) {
  MYB_CHECK(times != NULL && states != NULL && mean != NULL);
  if (!(mu > 0.0) || !isfinite(mu) || count < 3) {
    return MYB_E_DOMAIN;
  }
  double sum[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double t_sum = 0.0;
  for (size_t k = 0; k < count; ++k) {
    if (!isfinite(times[k]) || (k > 0 && !(times[k] > times[k - 1]))) {
      return MYB_E_DOMAIN;
    }
    t_sum += times[k];
  }
  const double t_mean = t_sum / (double)count;
  /* the mean longitude of each sample, unwrapped along the series, and the line through it:
   * sums for the least squares in (t - t_mean) */
  double previous = 0.0;
  double turns = 0.0;
  double s_l = 0.0;
  double s_tl = 0.0;
  double s_tt = 0.0;
  for (size_t k = 0; k < count; ++k) {
    myb_pv pv;
    for (int i = 0; i < 3; ++i) {
      pv.r[i] = states[6 * k + (size_t)i];
      pv.v[i] = states[6 * k + 3 + (size_t)i];
    }
    myb_equinoctial el;
    myb_status st = myb_equinoctial_from_pv(mu, &pv, &el);
    if (st != MYB_OK) {
      return st;
    }
    const double e2 = el.f * el.f + el.g * el.g;
    if (!(e2 < 1.0)) {
      return MYB_E_DOMAIN;
    }
    const double e = sqrt(e2);
    const double along = e > 0.0 ? myb_atan2(el.g, el.f) : 0.0; /* raan + argp */
    double eccentric = 0.0;
    double anomaly = 0.0;
    st = myb_anomaly_eccentric_from_true(el.l - along, e, &eccentric);
    if (st == MYB_OK) {
      st = myb_anomaly_mean_from_eccentric(eccentric, e, &anomaly);
    }
    if (st != MYB_OK) {
      return st;
    }
    double lambda = anomaly + along;
    lambda -= two_pi * floor(lambda / two_pi);
    if (k > 0) {
      /* the longitude only grows: whole turns are added until it is past the sample before
       * (samples less than a revolution apart) */
      const double behind = previous - (lambda + turns);
      if (behind > 0.0) {
        turns += two_pi * ceil(behind / two_pi);
      }
    }
    lambda += turns;
    previous = lambda;
    sum[0] += el.p / (1.0 - e2);
    sum[1] += el.f;
    sum[2] += el.g;
    sum[3] += el.h;
    sum[4] += el.k;
    const double dt = times[k] - t_mean;
    s_l += lambda;
    s_tl += dt * lambda;
    s_tt += dt * dt;
  }
  if (!(s_tt > 0.0)) {
    return MYB_E_DOMAIN;
  }
  mean->a = sum[0] / (double)count;
  mean->f = sum[1] / (double)count;
  mean->g = sum[2] / (double)count;
  mean->h = sum[3] / (double)count;
  mean->k = sum[4] / (double)count;
  mean->epoch = t_mean;
  mean->rate = s_tl / s_tt;
  double lambda0 = s_l / (double)count;
  lambda0 -= two_pi * floor(lambda0 / two_pi);
  mean->lambda = lambda0;
  return MYB_OK;
}
