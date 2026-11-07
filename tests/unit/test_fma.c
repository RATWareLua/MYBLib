/* SPDX-License-Identifier: BSD-3-Clause */
/* myb_fma: fixed cases with known exact answers (every platform) and a bit-for-bit
 * comparison with hardware FMA on several input distributions (CPUs with FMA3). */
#include <float.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/cmath.h"

double myb_test_hw_fma(double a, double b, double c);

static void test_known_values(void) {
  /* Single rounding: a*b = 1 + 2^-29 + 2^-60 exactly. */
  MT_CHECK_SAME(myb_fma(1.0 + 0x1p-30, 1.0 + 0x1p-30, -(1.0 + 0x1p-29)), 0x1p-60);

  /* Inputs on which mingw-w64's fma() is wrong; answers from hardware FMA. */
  MT_CHECK_SAME(myb_fma(-0x1.01dce4e7bfb79p-17, -0x1.a792e1af470eap-21, 0x1.474ffb8e8ab15p-29),
                0x1.48254f874f47bp-29);
  MT_CHECK_SAME(myb_fma(0x1.a4ef3496f4112p-16, 0x1.17f104eb30026p+6, -0x1.0273b10a6d4afp-26),
                0x1.cc4bfd23c8686p-10);
  MT_CHECK_SAME(myb_fma(-0x1.591e85a840907p-11, 0x1.7d0a0e9b446d4p-7, 0x1.eb56cc3fecd4dp-9),
                0x1.ea55f4420605fp-9);

  /* Signed zeros. */
  MT_CHECK_SAME(myb_fma(0.0, 1.0, -0.0), 0.0);
  MT_CHECK_SAME(myb_fma(-0.0, 1.0, -0.0), -0.0);
  MT_CHECK_SAME(myb_fma(-0.0, 1.0, 0.0), 0.0);
  MT_CHECK_SAME(myb_fma(1.0, 1.0, -1.0), 0.0); /* exact cancellation: +0 */
  MT_CHECK_SAME(myb_fma(-1.0, 1.0, 1.0), 0.0);
  MT_CHECK_SAME(myb_fma(0x1p-600, -0x1p-600, 0.0), -0.0); /* underflow keeps the sign */
  MT_CHECK_SAME(myb_fma(0x1p-600, -0x1p-600, -0.0), -0.0);

  /* Overflow: the product 2^1024 is not representable, the fused result is. */
  MT_CHECK_SAME(myb_fma(0x1p1023, 2.0, -0x1p1023), 0x1p1023);
  MT_CHECK_SAME(myb_fma(0x1p1023, 2.0, 0.0), (double)INFINITY);
  MT_CHECK_SAME(myb_fma(DBL_MAX, 1.0, DBL_MAX), (double)INFINITY);
  MT_CHECK_SAME(myb_fma(-DBL_MAX, 2.0, DBL_MAX), -DBL_MAX);
  /* ulp(DBL_MAX) = 2^971. Half an ulp is a tie; DBL_MAX has an odd significand, so the tie
   * rounds up, the carry leaves the exponent range and gives infinity. */
  MT_CHECK_SAME(myb_fma(DBL_MAX, 1.0, 0x1p970), (double)INFINITY);
  MT_CHECK_SAME(myb_fma(DBL_MAX, 1.0, 0x1p969), DBL_MAX); /* a quarter ulp: rounds down */

  /* Gradual underflow, ties to even at the bottom of the range. */
  MT_CHECK_SAME(myb_fma(0x1p-537, 0x1p-537, 0.0), 0x1p-1074);
  MT_CHECK_SAME(myb_fma(0x1p-538, 0x1p-537, 0.0), 0.0);         /* 2^-1075: tie, even */
  MT_CHECK_SAME(myb_fma(0x1.8p-538, 0x1p-537, 0.0), 0x1p-1074); /* 0.75 ulp: up */
  MT_CHECK_SAME(myb_fma(0x1p-1074, 1.5, 0.0), 0x1p-1073);       /* 1.5 ulp: tie, even */
  MT_CHECK_SAME(myb_fma(0x1p-1074, 2.5, 0.0), 0x1p-1073);       /* 2.5 ulp: tie, even */
  MT_CHECK_SAME(myb_fma(0x1p-1022, 0x1p-52, 0.0), 0x1p-1074);
  /* Largest subnormal + smallest subnormal = smallest normal (carry into the exponent). */
  MT_CHECK_SAME(myb_fma(0x0.fffffffffffffp-1022, 1.0, 0x1p-1074), 0x1p-1022);
  MT_CHECK_SAME(myb_fma(-0x1p-1074, 0.5, -0x1p-1074), -0x1p-1073); /* 1.5 ulp: tie, even */

  /* Infinities and NaN. */
  MT_CHECK(myb_fma((double)INFINITY, 0.0, 1.0) != myb_fma((double)INFINITY, 0.0, 1.0));
  MT_CHECK(myb_fma((double)INFINITY, 1.0, -(double)INFINITY) !=
           myb_fma((double)INFINITY, 1.0, -(double)INFINITY));
  MT_CHECK_SAME(myb_fma(1.0, 1.0, (double)INFINITY), (double)INFINITY);
  MT_CHECK_SAME(myb_fma(DBL_MAX, 2.0, -(double)INFINITY), -(double)INFINITY);
  MT_CHECK(myb_fma((double)NAN, 1.0, 1.0) != myb_fma((double)NAN, 1.0, 1.0));
  MT_CHECK(myb_fma(1.0, 1.0, (double)NAN) != myb_fma(1.0, 1.0, (double)NAN));
}

typedef void (*fma_generator)(mt_rng* rng, double* a, double* b, double* c);

static void gen_any(mt_rng* rng, double* a, double* b, double* c) {
  *a = mt_any_double(rng);
  *b = mt_any_double(rng);
  *c = mt_any_double(rng);
}

static void gen_moderate(mt_rng* rng, double* a, double* b, double* c) {
  *a = mt_double(rng, -60, 60);
  *b = mt_double(rng, -60, 60);
  *c = mt_double(rng, -120, 120);
}

/* c close to -a*b: deep cancellation, the hardest case for alignment and sticky bits. */
static void gen_cancellation(mt_rng* rng, double* a, double* b, double* c) {
  *a = mt_double(rng, -30, 30);
  *b = mt_double(rng, -30, 30);
  const double p = *a * *b;
  const uint64_t delta = mt_next(rng) % 64;
  const uint64_t bits = mt_bits(p);
  *c = -mt_from_bits((mt_next(rng) & 1) != 0 ? bits + delta : bits - delta);
}

static void gen_subnormal(mt_rng* rng, double* a, double* b, double* c) {
  *a = mt_double(rng, -560, -500);
  *b = mt_double(rng, -560, -500);
  *c = (mt_next(rng) & 1) != 0 ? mt_from_bits(mt_next(rng) & UINT64_C(0x800FFFFFFFFFFFFF))
                               : mt_double(rng, -1022, -1000);
}

static void gen_overflow(mt_rng* rng, double* a, double* b, double* c) {
  *a = mt_double(rng, 500, 520);
  *b = mt_double(rng, 500, 520);
  *c = mt_double(rng, 1000, 1023);
}

/* Operands many binades apart: sticky-bit path. */
static void gen_far(mt_rng* rng, double* a, double* b, double* c) {
  *a = mt_double(rng, -200, 200);
  *b = mt_double(rng, -200, 200);
  *c = mt_double(rng, -1022, 1023);
}

static void compare_with_hardware(const char* name, fma_generator generate, long count,
                                  uint64_t seed) {
  mt_rng rng = {.state = seed};
  long mismatches = 0;
  for (long i = 0; i < count; ++i) {
    double a;
    double b;
    double c;
    generate(&rng, &a, &b, &c);
    const double expected = myb_test_hw_fma(a, b, c);
    const double actual = myb_fma(a, b, c);
    mtest.checks++;
    if (!mt_same(actual, expected)) {
      mtest.failures++;
      if (++mismatches <= 5) {
        fprintf(stderr, "%s: fma(%a, %a, %a) = %a, hardware %a\n", name, a, b, c, actual, expected);
      }
    }
  }
  printf("%-12s %ld inputs, %ld mismatches with hardware FMA\n", name, count, mismatches);
}

int main(void) {
  test_known_values();
  if (!mt_cpu_has_fma()) {
    printf("no FMA on this CPU: hardware comparison skipped\n");
    const int status = mt_finish("test_fma");
    return status != 0 ? status : MTEST_SKIP;
  }

  enum { COUNT = 2000000 };

  compare_with_hardware("any-bits", gen_any, COUNT, 1);
  compare_with_hardware("moderate", gen_moderate, COUNT, 2);
  compare_with_hardware("cancel", gen_cancellation, COUNT, 3);
  compare_with_hardware("subnormal", gen_subnormal, COUNT, 4);
  compare_with_hardware("overflow", gen_overflow, COUNT, 5);
  compare_with_hardware("far-apart", gen_far, COUNT, 6);
  return mt_finish("test_fma");
}
