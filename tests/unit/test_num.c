/* SPDX-License-Identifier: BSD-3-Clause */
/* Error-free transformations, double-word arithmetic and compensated sums: properties that
 * hold exactly and can be checked without a multiple-precision library. Error bounds against
 * MPFR are in tests/mpfr/test_mpfr_num.c. */
#include <float.h>
#include <math.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/cmath.h"
#include "myb/num.h"

static const double u = 0x1p-53;

static void test_two_sum(void) {
  mt_rng rng = {.state = 11};
  for (int i = 0; i < 200000; ++i) {
    const double a = mt_double(&rng, -200, 200);
    const double b = mt_double(&rng, -200, 200);
    double err;
    const double s = myb_two_sum(a, b, &err);
    MT_CHECK_SAME(s, a + b);
    /* Independent algorithm: Fast2Sum with the operands ordered by magnitude. */
    const double big = fabs(a) >= fabs(b) ? a : b;
    const double small = fabs(a) >= fabs(b) ? b : a;
    const double s2 = big + small;
    const double err2 = small - (s2 - big);
    MT_CHECK_SAME(err, err2);
  }
  double err;
  MT_CHECK_SAME(myb_two_sum(1.0, 0x1p-60, &err), 1.0);
  MT_CHECK_SAME(err, 0x1p-60);
  /* ulp(2^1023) = 2^971: adding half an ulp is a tie, rounded to the even 2^1023; the error
   * term carries the whole addend. */
  MT_CHECK_SAME(myb_two_sum(0x1p1023, 0x1p970, &err), 0x1p1023);
  MT_CHECK_SAME(err, 0x1p970);
  MT_CHECK_SAME(myb_two_sum(DBL_MAX, -0x1p1023, &err), 0x1.ffffffffffffep1022);
  MT_CHECK_SAME(err, 0.0);
}

static void test_two_prod(void) {
  mt_rng rng = {.state = 12};
  for (int i = 0; i < 200000; ++i) {
    /* e_a + e_b >= -970 and no overflow: the transformation is exact, and the error term
     * equals RN(a*b - p) = fma(a, b, -p) (myb_fma is checked against hardware FMA). */
    const int range = (i % 3 == 0) ? 1000 : 300;
    const double a = mt_double(&rng, -range / 2, range / 2);
    const double b = mt_double(&rng, -range / 2, range / 2);
    double err;
    const double p = myb_two_prod(a, b, &err);
    MT_CHECK_SAME(p, a * b);
    MT_CHECK_SAME(err, myb_fma(a, b, -p));
  }
  /* Operands above 2^996 take the scaled splitting path. */
  for (int i = 0; i < 20000; ++i) {
    const double a = mt_double(&rng, 997, 1010);
    const double b = mt_double(&rng, -40, -10);
    double err;
    const double p = myb_two_prod(a, b, &err);
    MT_CHECK_SAME(p, a * b);
    MT_CHECK_SAME(err, myb_fma(a, b, -p));
  }
}

/* |dd(x) - 1| <= k*u^2, computed exactly enough with one more double-word subtraction. */
static int dd_near_one(myb_dd x, double k) {
  const myb_dd d = myb_dd_add_d(x, -1.0);
  return fabs(d.hi) <= k * u * u;
}

static void test_dd_identities(void) {
  const myb_dd one = myb_dd_from_double(1.0);
  const myb_dd third = myb_dd_div_d(one, 3.0);
  MT_CHECK(third.hi == 1.0 / 3.0);
  MT_CHECK(third.lo != 0.0);
  MT_CHECK(dd_near_one(myb_dd_mul_d(third, 3.0), 6.0));

  mt_rng rng = {.state = 13};
  for (int i = 0; i < 100000; ++i) {
    const double hi = mt_double(&rng, -100, 100);
    const double lo = hi * 0x1p-53 * (double)(int64_t)(mt_next(&rng) % 2001 - 1000) / 1000.0;
    double e;
    const double h = myb_two_sum(hi, lo, &e);
    const myb_dd x = {.hi = h, .lo = e};
    const myb_dd y = {.hi = mt_double(&rng, -100, 100), .lo = 0.0};

    const myb_dd zero = myb_dd_sub(x, x);
    MT_CHECK(zero.hi == 0.0 && zero.lo == 0.0);
    const myb_dd same = myb_dd_mul(x, one);
    MT_CHECK(same.hi == x.hi && same.lo == x.lo);
    MT_CHECK(dd_near_one(myb_dd_div(x, x), 16.0));
    /* (x*y)/y = x within the sum of the bounds of mul and div. */
    const myb_dd back = myb_dd_div(myb_dd_mul(x, y), y);
    MT_CHECK(dd_near_one(myb_dd_div(back, x), 40.0));
    const myb_dd abs_x = x.hi < 0.0 ? myb_dd_neg(x) : x;
    const myb_dd r = myb_dd_sqrt(abs_x);
    MT_CHECK(dd_near_one(myb_dd_div(myb_dd_mul(r, r), abs_x), 40.0));
  }
  const myb_dd negzero = myb_dd_sqrt(myb_dd_from_double(-0.0));
  MT_CHECK_SAME(negzero.hi, -0.0);
}

static void test_sums(void) {
  const double cancel[] = {1e16, 1.0, -1e16};
  MT_CHECK_SAME(myb_sum2(cancel, 3), 1.0);
  MT_CHECK_SAME(myb_sum2(cancel, 0), 0.0);
  MT_CHECK_SAME(myb_dot2(cancel, cancel, 0), 0.0);

  /* sum = t exactly (pairs x, -x around it); Sum2 bound: |res - t| <= u|t| + g^2 sum|x|. */
  enum { N = 2001 };

  static double x[N];
  static double y[N];
  mt_rng rng = {.state = 14};
  for (int trial = 0; trial < 50; ++trial) {
    double abs_sum = 0.0;
    for (int i = 0; i < N - 1; i += 2) {
      x[i] = mt_double(&rng, -10, 60);
      x[i + 1] = -x[i];
      abs_sum += 2.0 * fabs(x[i]);
    }
    const double t = mt_double(&rng, -5, 5);
    x[N - 1] = t;
    abs_sum += fabs(t);
    /* Deterministic shuffle so that the cancelling pairs are far apart. */
    for (int i = N - 1; i > 0; --i) {
      const int j = (int)(mt_next(&rng) % (uint64_t)(i + 1));
      const double tmp = x[i];
      x[i] = x[j];
      x[j] = tmp;
    }
    const double gamma = (double)(N - 1) * u / (1.0 - (double)(N - 1) * u);
    const double bound = u * fabs(t) + gamma * gamma * abs_sum;
    const double sum = myb_sum2(x, N);
    MT_CHECK(fabs(sum - t) <= bound);

    myb_sum_acc acc = {};
    for (int i = 0; i < N; ++i) {
      myb_sum_acc_add(&acc, x[i]);
    }
    MT_CHECK_SAME(myb_sum_acc_value(&acc), sum);

    /* Dot product of the same construction: x . y with y = 1 except the pair partners. */
    for (int i = 0; i < N; ++i) {
      y[i] = 1.0;
    }
    MT_CHECK_SAME(myb_dot2(x, y, N), sum);
  }
}

int main(void) {
  test_two_sum();
  test_two_prod();
  test_dd_identities();
  test_sums();
  return mt_finish("test_num");
}
