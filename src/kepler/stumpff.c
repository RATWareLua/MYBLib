/* SPDX-License-Identifier: BSD-3-Clause */
#include "kepler/stumpff.h"

#include <math.h> /* sqrt only */

#include "core/base.h"
#include "myb/cmath.h"
#include "myb/num.h"

/* 1/(2k+3)! for k = 0..13: the c3 series sum (-z)^k / (2k+3)! converges for all z; for
 * |z| < 1 the 14th term is below 2^-90 of the first. Folded exactly by the compiler as
 * IEEE divisions; the accumulated rounding of the constants is below 15 ulp of each
 * coefficient, i.e. far below the truncation of the last terms. */
static constexpr double inverse_factorial_odd[14] = {
    1.0 / 6.0,
    1.0 / 6.0 / 4.0 / 5.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0 / 10.0 / 11.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0 / 10.0 / 11.0 / 12.0 / 13.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0 / 10.0 / 11.0 / 12.0 / 13.0 / 14.0 / 15.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0 / 10.0 / 11.0 / 12.0 / 13.0 / 14.0 / 15.0 / 16.0 /
        17.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0 / 10.0 / 11.0 / 12.0 / 13.0 / 14.0 / 15.0 / 16.0 /
        17.0 / 18.0 / 19.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0 / 10.0 / 11.0 / 12.0 / 13.0 / 14.0 / 15.0 / 16.0 /
        17.0 / 18.0 / 19.0 / 20.0 / 21.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0 / 10.0 / 11.0 / 12.0 / 13.0 / 14.0 / 15.0 / 16.0 /
        17.0 / 18.0 / 19.0 / 20.0 / 21.0 / 22.0 / 23.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0 / 10.0 / 11.0 / 12.0 / 13.0 / 14.0 / 15.0 / 16.0 /
        17.0 / 18.0 / 19.0 / 20.0 / 21.0 / 22.0 / 23.0 / 24.0 / 25.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0 / 10.0 / 11.0 / 12.0 / 13.0 / 14.0 / 15.0 / 16.0 /
        17.0 / 18.0 / 19.0 / 20.0 / 21.0 / 22.0 / 23.0 / 24.0 / 25.0 / 26.0 / 27.0,
    1.0 / 6.0 / 4.0 / 5.0 / 6.0 / 7.0 / 8.0 / 9.0 / 10.0 / 11.0 / 12.0 / 13.0 / 14.0 / 15.0 / 16.0 /
        17.0 / 18.0 / 19.0 / 20.0 / 21.0 / 22.0 / 23.0 / 24.0 / 25.0 / 26.0 / 27.0 / 28.0 / 29.0,
};

double myb_internal_stumpff_c3(double z) {
  if (z > -1.0 && z < 1.0) {
    /* Horner in -z from the smallest term. */
    double sum = 0.0;
    for (int k = 13; k >= 0; --k) {
      sum = inverse_factorial_odd[k] - z * sum;
    }
    return sum;
  }
  if (z > 0.0) {
    const double s = sqrt(z);
    return (s - myb_sin(s)) / (z * s);
  }
  const double s = sqrt(-z);
  return (myb_sinh(s) - s) / (-z * s);
}

double myb_internal_stumpff_c2(double z) {
  if (z == 0.0) {
    return 0.5;
  }
  if (z > 0.0) {
    const double s = sqrt(z);
    const double q = myb_sin(0.5 * s) / s;
    return 2.0 * q * q;
  }
  const double s = sqrt(-z);
  const double q = myb_sinh(0.5 * s) / s;
  return 2.0 * q * q;
}

/* sum_k (-z)^k / (2k + first)! / (1 / first!) in nested form; for |z| < 4 the 16th factor is
 * below 2^-60. */
static double stumpff_series(double z, int first) {
  double inner = 1.0;
  for (int k = 16; k >= 1; --k) {
    const double m = (double)(2 * k + first - 1) * (double)(2 * k + first);
    inner = 1.0 - z * inner / m;
  }
  return inner;
}

double myb_internal_stumpff_c4(double z) {
  if (z > -4.0 && z < 4.0) {
    return stumpff_series(z, 4) / 24.0;
  }
  return (0.5 - myb_internal_stumpff_c2(z)) / z; /* cancellation below a factor 3.5 */
}

double myb_internal_stumpff_c5(double z) {
  if (z > -4.0 && z < 4.0) {
    return stumpff_series(z, 5) / 120.0;
  }
  return (1.0 / 6.0 - myb_internal_stumpff_c3(z)) / z;
}

double myb_internal_stumpff_c0(double z) {
  if (z >= 0.0) {
    return myb_cos(sqrt(z));
  }
  return myb_cosh(sqrt(-z));
}

double myb_internal_stumpff_c1(double z) {
  if (z == 0.0) {
    return 1.0;
  }
  if (z > 0.0) {
    const double s = sqrt(z);
    return myb_sin(s) / s;
  }
  const double s = sqrt(-z);
  return myb_sinh(s) / s;
}

/* 2*pi = two_pi_hi + two_pi_lo to 106 bits. */
static const double two_pi_hi = 0x1.921fb54442d18p+2;
static const double two_pi_lo = 0x1.1a62633145c07p-52;
static const double pi_hi = 0x1.921fb54442d18p+1;

double myb_internal_add_turns(double x, double turns) {
  double err;
  const double p = myb_two_prod(turns, two_pi_hi, &err);
  myb_dd sum = myb_dd_add_d(myb_dd_from_double(p), x);
  sum = myb_dd_add_d(sum, err + turns * two_pi_lo);
  return myb_dd_to_double(sum);
}

double myb_internal_reduce_angle(double x, double* turns) {
  MYB_CHECK(turns != NULL);
  double k = __builtin_roundeven(x / two_pi_hi);
  double r = myb_internal_add_turns(x, -k);
  /* The quotient may round across a half turn: one correction step at most. */
  if (r > pi_hi) {
    k += 1.0;
    r = myb_internal_add_turns(x, -k);
  } else if (r < -pi_hi) {
    k -= 1.0;
    r = myb_internal_add_turns(x, -k);
  }
  *turns = k;
  return r;
}
