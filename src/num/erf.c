/* SPDX-License-Identifier: BSD-3-Clause AND SunPro */
/* The error function and its complement: the fdlibm algorithm (FreeBSD lib/msun/src/s_erf.c,
 * revision of the main branch d02883b3d8d0943fdb98cbc365074263f918908b, SHA-256
 * 5a2f06bd1d0452456cda86fa1879b956215a7396c29fd8306491f4a9cd9a4c21), with exp replaced by the
 * correctly rounded myb_exp and the word access by memcpy. Error about 1 ulp; the same bits on
 * every platform. The original notice:
 *
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
#include "num/erf.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "myb/cmath.h"

static const double erx = 8.45062911510467529297e-01;
static const double efx = 1.28379167095512586316e-01;
static const double efx8 = 1.02703333676410069053e+00;
static const double pp0 = 1.28379167095512558561e-01;
static const double pp1 = -3.25042107247001499370e-01;
static const double pp2 = -2.84817495755985104766e-02;
static const double pp3 = -5.77027029648944159157e-03;
static const double pp4 = -2.37630166566501626084e-05;
static const double qq1 = 3.97917223959155352819e-01;
static const double qq2 = 6.50222499887672944485e-02;
static const double qq3 = 5.08130628187576562776e-03;
static const double qq4 = 1.32494738004321644526e-04;
static const double qq5 = -3.96022827877536812320e-06;
static const double pa0 = -2.36211856075265944077e-03;
static const double pa1 = 4.14856118683748331666e-01;
static const double pa2 = -3.72207876035701323847e-01;
static const double pa3 = 3.18346619901161753674e-01;
static const double pa4 = -1.10894694282396677476e-01;
static const double pa5 = 3.54783043256182359371e-02;
static const double pa6 = -2.16637559486879084300e-03;
static const double qa1 = 1.06420880400844228286e-01;
static const double qa2 = 5.40397917702171048937e-01;
static const double qa3 = 7.18286544141962662868e-02;
static const double qa4 = 1.26171219808761642112e-01;
static const double qa5 = 1.36370839120290507362e-02;
static const double qa6 = 1.19844998467991074170e-02;
static const double ra0 = -9.86494403484714822705e-03;
static const double ra1 = -6.93858572707181764372e-01;
static const double ra2 = -1.05586262253232909814e+01;
static const double ra3 = -6.23753324503260060396e+01;
static const double ra4 = -1.62396669462573470355e+02;
static const double ra5 = -1.84605092906711035994e+02;
static const double ra6 = -8.12874355063065934246e+01;
static const double ra7 = -9.81432934416914548592e+00;
static const double sa1 = 1.96512716674392571292e+01;
static const double sa2 = 1.37657754143519042600e+02;
static const double sa3 = 4.34565877475229228821e+02;
static const double sa4 = 6.45387271733267880336e+02;
static const double sa5 = 4.29008140027567833386e+02;
static const double sa6 = 1.08635005541779435134e+02;
static const double sa7 = 6.57024977031928170135e+00;
static const double sa8 = -6.04244152148580987438e-02;
static const double rb0 = -9.86494292470009928597e-03;
static const double rb1 = -7.99283237680523006574e-01;
static const double rb2 = -1.77579549177547519889e+01;
static const double rb3 = -1.60636384855821916062e+02;
static const double rb4 = -6.37566443368389627722e+02;
static const double rb5 = -1.02509513161107724954e+03;
static const double rb6 = -4.83519191608651397019e+02;
static const double sb1 = 3.03380607434824582924e+01;
static const double sb2 = 3.25792512996573918826e+02;
static const double sb3 = 1.53672958608443695994e+03;
static const double sb4 = 3.19985821950859553908e+03;
static const double sb5 = 2.55305040643316442583e+03;
static const double sb6 = 4.74528541206955367215e+02;
static const double sb7 = -2.24409524465858183362e+01;

static int32_t high_word(double x) {
  uint64_t u = 0;
  memcpy(&u, &x, sizeof u);
  return (int32_t)(u >> 32);
}

static double clear_low_word(double x) {
  uint64_t u = 0;
  memcpy(&u, &x, sizeof u);
  u &= 0xffffffff00000000ULL;
  double z = 0.0;
  memcpy(&z, &u, sizeof z);
  return z;
}

/* erfc(x) for 1.25 <= x < 28 (x positive): exp(-x^2 - 0.5625 + R/S) / x */
static double tail(double x, int32_t ix) {
  const double s = 1.0 / (x * x);
  double r = 0.0;
  double q = 0.0;
  if (ix < 0x4006DB6E) { /* |x| < 1/0.35 */
    r = ra0 + s * (ra1 + s * (ra2 + s * (ra3 + s * (ra4 + s * (ra5 + s * (ra6 + s * ra7))))));
    q = 1.0 + s * (sa1 + s * (sa2 + s * (sa3 + s * (sa4 + s * (sa5 + s * (sa6 + s * (sa7 +
                                                                              s * sa8)))))));
  } else {
    r = rb0 + s * (rb1 + s * (rb2 + s * (rb3 + s * (rb4 + s * (rb5 + s * rb6)))));
    q = 1.0 + s * (sb1 + s * (sb2 + s * (sb3 + s * (sb4 + s * (sb5 + s * (sb6 + s * sb7))))));
  }
  const double z = clear_low_word(x);
  return myb_exp(-z * z - 0.5625) * myb_exp((z - x) * (z + x) + r / q) / x;
}

double myb_internal_erf(double x) {
  const int32_t hx = high_word(x);
  const int32_t ix = hx & 0x7fffffff;
  if (ix >= 0x7ff00000) {
    return isnan(x) ? x : (x > 0.0 ? 1.0 : -1.0);
  }
  if (ix < 0x3feb0000) { /* |x| < 0.84375 */
    if (ix < 0x3e300000) { /* |x| < 2^-28 */
      if (ix < 0x00800000) {
        return (8.0 * x + efx8 * x) / 8.0;
      }
      return x + efx * x;
    }
    const double z = x * x;
    const double r = pp0 + z * (pp1 + z * (pp2 + z * (pp3 + z * pp4)));
    const double s = 1.0 + z * (qq1 + z * (qq2 + z * (qq3 + z * (qq4 + z * qq5))));
    return x + x * (r / s);
  }
  if (ix < 0x3ff40000) { /* 0.84375 <= |x| < 1.25 */
    const double s = fabs(x) - 1.0;
    const double p = pa0 + s * (pa1 + s * (pa2 + s * (pa3 + s * (pa4 + s * (pa5 + s * pa6)))));
    const double q = 1.0 + s * (qa1 + s * (qa2 + s * (qa3 + s * (qa4 + s * (qa5 + s * qa6)))));
    return hx >= 0 ? erx + p / q : -erx - p / q;
  }
  if (ix >= 0x40180000) { /* |x| >= 6 */
    return hx >= 0 ? 1.0 : -1.0;
  }
  const double r = tail(fabs(x), ix);
  return hx >= 0 ? 1.0 - r : r - 1.0;
}

double myb_internal_erfc(double x) {
  const int32_t hx = high_word(x);
  const int32_t ix = hx & 0x7fffffff;
  if (ix >= 0x7ff00000) {
    return isnan(x) ? x : (x > 0.0 ? 0.0 : 2.0);
  }
  if (ix < 0x3feb0000) { /* |x| < 0.84375 */
    if (ix < 0x3c700000) { /* |x| < 2^-56 */
      return 1.0 - x;
    }
    const double z = x * x;
    const double r = pp0 + z * (pp1 + z * (pp2 + z * (pp3 + z * pp4)));
    const double s = 1.0 + z * (qq1 + z * (qq2 + z * (qq3 + z * (qq4 + z * qq5))));
    const double y = r / s;
    if (hx < 0x3fd00000) { /* x < 1/4 */
      return 1.0 - (x + x * y);
    }
    double t = x * y;
    t += (x - 0.5);
    return 0.5 - t;
  }
  if (ix < 0x3ff40000) { /* 0.84375 <= |x| < 1.25 */
    const double s = fabs(x) - 1.0;
    const double p = pa0 + s * (pa1 + s * (pa2 + s * (pa3 + s * (pa4 + s * (pa5 + s * pa6)))));
    const double q = 1.0 + s * (qa1 + s * (qa2 + s * (qa3 + s * (qa4 + s * (qa5 + s * qa6)))));
    if (hx >= 0) {
      const double z = 1.0 - erx;
      return z - p / q;
    }
    const double z = erx + p / q;
    return 1.0 + z;
  }
  if (ix < 0x403c0000) { /* |x| < 28 */
    if (hx < 0 && ix >= 0x40180000) { /* x < -6 */
      return 2.0;
    }
    /* the original switches the rational approximation at 0x4006DB6D here (one ulp-sized
     * difference from erf; kept) */
    const double ax = fabs(x);
    const double s = 1.0 / (ax * ax);
    double r = 0.0;
    double q = 0.0;
    if (ix < 0x4006DB6D) {
      r = ra0 + s * (ra1 + s * (ra2 + s * (ra3 + s * (ra4 + s * (ra5 + s * (ra6 + s * ra7))))));
      q = 1.0 + s * (sa1 + s * (sa2 + s * (sa3 + s * (sa4 + s * (sa5 + s * (sa6 + s * (sa7 +
                                                                                s * sa8)))))));
    } else {
      r = rb0 + s * (rb1 + s * (rb2 + s * (rb3 + s * (rb4 + s * (rb5 + s * rb6)))));
      q = 1.0 + s * (sb1 + s * (sb2 + s * (sb3 + s * (sb4 + s * (sb5 + s * (sb6 + s * sb7))))));
    }
    const double z = clear_low_word(ax);
    const double e = myb_exp(-z * z - 0.5625) * myb_exp((z - ax) * (z + ax) + r / q);
    return hx > 0 ? e / ax : 2.0 - e / ax;
  }
  return hx > 0 ? 0.0 : 2.0;
}
