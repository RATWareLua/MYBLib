/* SPDX-License-Identifier: BSD-3-Clause */
/* Every function of myb/cmath.h is correctly rounded: compared bit for bit with GNU MPFR on
 * random inputs over the whole domain (normal and subnormal), special values, hard argument
 * reductions near multiples of pi/2, and samples of the CORE-MATH worst cases for correct
 * rounding (tests/data/core-math). */
#include <float.h>
#include <mpfr.h>
#include <stdlib.h>
#include <string.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/cmath.h"

typedef int (*mpfr_unary)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
typedef int (*mpfr_binary)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);

typedef enum sign_mode { BOTH_SIGNS, POSITIVE_ONLY } sign_mode;

typedef struct unary {
  const char* name;
  double (*f)(double);
  mpfr_unary reference;
  int emin; /* binary exponent range of random inputs */
  int emax;
  int emax_negative; /* for negative inputs (log1p: |x| < 1) */
  sign_mode signs;
} unary;

typedef struct binary {
  const char* name;
  double (*f)(double, double);
  mpfr_binary reference;
} binary;

static long shown = 0;

/* RN(f(x)) in binary64 with gradual underflow (emin/emax of binary64 are set in main). */
static double reference_unary(mpfr_unary f, double x) {
  mpfr_t y;
  mpfr_init2(y, 53);
  mpfr_set_d(y, x, MPFR_RNDN);
  const int inexact = f(y, y, MPFR_RNDN);
  mpfr_subnormalize(y, inexact, MPFR_RNDN);
  const double r = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear(y);
  return r;
}

static double reference_binary(mpfr_binary f, double a, double b) {
  mpfr_t x;
  mpfr_t y;
  mpfr_t r;
  mpfr_inits2(53, x, y, r, (mpfr_ptr)0);
  mpfr_set_d(x, a, MPFR_RNDN);
  mpfr_set_d(y, b, MPFR_RNDN);
  const int inexact = f(r, x, y, MPFR_RNDN);
  mpfr_subnormalize(r, inexact, MPFR_RNDN);
  const double result = mpfr_get_d(r, MPFR_RNDN);
  mpfr_clears(x, y, r, (mpfr_ptr)0);
  return result;
}

static long check_unary(const unary* s, double x) {
  const double expected = reference_unary(s->reference, x);
  const double actual = s->f(x);
  mtest.checks++;
  if (mt_same(actual, expected)) {
    return 0;
  }
  mtest.failures++;
  if (++shown <= 20) {
    fprintf(stderr, "%s(%a) = %a, MPFR %a\n", s->name, x, actual, expected);
  }
  return 1;
}

static long check_binary(const binary* s, double a, double b) {
  const double expected = reference_binary(s->reference, a, b);
  const double actual = s->f(a, b);
  mtest.checks++;
  if (mt_same(actual, expected)) {
    return 0;
  }
  mtest.failures++;
  if (++shown <= 20) {
    fprintf(stderr, "%s(%a, %a) = %a, MPFR %a\n", s->name, a, b, actual, expected);
  }
  return 1;
}

static const double specials[] = {0.0,
                                  -0.0,
                                  1.0,
                                  -1.0,
                                  0.5,
                                  -0.5,
                                  2.0,
                                  -2.0,
                                  DBL_MIN,
                                  -DBL_MIN,
                                  0x1p-1074,
                                  -0x1p-1074,
                                  DBL_MAX,
                                  -DBL_MAX,
                                  0x1p-26,
                                  -0x1p-27,
                                  709.78,
                                  -745.2,
                                  710.5,
                                  (double)INFINITY,
                                  -(double)INFINITY,
                                  (double)NAN};

static FILE* open_worst(const char* name) {
  char path[1024];
  const int n = snprintf(path, sizeof path, "%s/core-math/%s.wc", MYB_TEST_DATA_DIR, name);
  MT_CHECK(n > 0 && (size_t)n < sizeof path);
  FILE* f = fopen(path, "r");
  MT_CHECK(f != NULL);
  return f;
}

static void test_unary(const unary* s, uint64_t seed) {
  long wrong = 0;
  long count = 0;
  for (size_t i = 0; i < sizeof specials / sizeof specials[0]; ++i) {
    wrong += check_unary(s, specials[i]);
    ++count;
  }
  mt_rng rng = {.state = seed};
  for (int i = 0; i < 60000; ++i) {
    const int negative = s->signs == BOTH_SIGNS && (mt_next(&rng) & 1) != 0;
    const int emax = negative ? s->emax_negative : s->emax;
    double x = mt_double(&rng, s->emin, emax);
    x = negative ? -__builtin_fabs(x) : __builtin_fabs(x);
    wrong += check_unary(s, x);
    ++count;
  }
  if (s->emin <= -1022) {
    for (int i = 0; i < 5000; ++i) {
      const uint64_t sign = s->signs == BOTH_SIGNS ? (mt_next(&rng) & (UINT64_C(1) << 63)) : 0;
      wrong += check_unary(s, mt_from_bits(sign | (mt_next(&rng) & UINT64_C(0x000FFFFFFFFFFFFF))));
      ++count;
    }
  }
  FILE* f = open_worst(s->name);
  long worst = 0;
  char line[256];
  /* bound: one iteration per line of a finite file */
  while (f != NULL && fgets(line, sizeof line, f) != NULL) {
    if (line[0] == '#') {
      continue;
    }
    const double x = strtod(line, NULL);
    wrong += check_unary(s, x);
    worst += 1;
    if (s->signs == BOTH_SIGNS) {
      wrong += check_unary(s, -x);
      worst += 1;
    }
  }
  if (f != NULL) {
    fclose(f);
  }
  MT_CHECK(worst > 0);
  printf("%-6s %7ld inputs (%4ld worst cases), %ld not correctly rounded\n", s->name, count + worst,
         worst, wrong);
}

static void test_binary(const binary* s, uint64_t seed) {
  long wrong = 0;
  long count = 0;
  const size_t n = sizeof specials / sizeof specials[0];
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      wrong += check_binary(s, specials[i], specials[j]);
      ++count;
    }
  }
  mt_rng rng = {.state = seed};
  for (int i = 0; i < 60000; ++i) {
    const double a = mt_double(&rng, -1022, 1023);
    /* Second argument at a random distance in binades, including nearly equal ones. */
    const int shift = mt_int(&rng, -80, 80);
    const double b = __builtin_ldexp(mt_double(&rng, -1, 0) * a, shift);
    wrong += check_binary(s, a, (mt_next(&rng) & 1) != 0 ? b : -b);
    ++count;
  }
  FILE* f = open_worst(s->name);
  long worst = 0;
  char line[256];
  /* bound: one iteration per line of a finite file */
  while (f != NULL && fgets(line, sizeof line, f) != NULL) {
    if (line[0] == '#') {
      continue;
    }
    char* comma = NULL;
    const double a = strtod(line, &comma);
    const double b = strtod(comma + 1, NULL);
    wrong += check_binary(s, a, b);
    wrong += check_binary(s, -a, b);
    wrong += check_binary(s, a, -b);
    worst += 3;
  }
  if (f != NULL) {
    fclose(f);
  }
  MT_CHECK(worst > 0);
  printf("%-6s %7ld inputs (%4ld worst cases), %ld not correctly rounded\n", s->name, count + worst,
         worst, wrong);
}

/* Doubles nearest to k*pi/2 and their neighbours: the hardest argument reductions. */
static void test_reduction(const unary* s) {
  long wrong = 0;
  long count = 0;
  mpfr_t t;
  mpfr_init2(t, 2000);
  for (long k = 1; k <= 4000; ++k) {
    mpfr_const_pi(t, MPFR_RNDN);
    mpfr_mul_si(t, t, k, MPFR_RNDN);
    mpfr_div_2ui(t, t, 1, MPFR_RNDN);
    const double centre = mpfr_get_d(t, MPFR_RNDN);
    for (int delta = -3; delta <= 3; ++delta) {
      const double x = mt_from_bits(mt_bits(centre) + (uint64_t)(int64_t)delta);
      wrong += check_unary(s, x) + check_unary(s, -x);
      count += 2;
    }
  }
  mpfr_clear(t);
  printf("%-6s %7ld inputs near k*pi/2, %ld not correctly rounded\n", s->name, count, wrong);
}

int main(void) {
  mpfr_set_emin(-1073);
  mpfr_set_emax(1024);
  const unary unaries[] = {
      {"sin", myb_sin, mpfr_sin, -1022, 1023, 1023, BOTH_SIGNS},
      {"cos", myb_cos, mpfr_cos, -1022, 1023, 1023, BOTH_SIGNS},
      {"tan", myb_tan, mpfr_tan, -1022, 1023, 1023, BOTH_SIGNS},
      {"asin", myb_asin, mpfr_asin, -1022, -1, -1, BOTH_SIGNS},
      {"acos", myb_acos, mpfr_acos, -1022, -1, -1, BOTH_SIGNS},
      {"atan", myb_atan, mpfr_atan, -1022, 1023, 1023, BOTH_SIGNS},
      {"exp", myb_exp, mpfr_exp, -1022, 10, 10, BOTH_SIGNS},
      {"expm1", myb_expm1, mpfr_expm1, -1022, 10, 10, BOTH_SIGNS},
      {"log", myb_log, mpfr_log, -1022, 1023, 1023, POSITIVE_ONLY},
      {"log1p", myb_log1p, mpfr_log1p, -1022, 1023, -1, BOTH_SIGNS},
      {"sinh", myb_sinh, mpfr_sinh, -1022, 10, 10, BOTH_SIGNS},
      {"cosh", myb_cosh, mpfr_cosh, -1022, 10, 10, BOTH_SIGNS},
      {"tanh", myb_tanh, mpfr_tanh, -1022, 6, 6, BOTH_SIGNS},
      {"asinh", myb_asinh, mpfr_asinh, -1022, 1023, 1023, BOTH_SIGNS},
      {"acosh", myb_acosh, mpfr_acosh, 0, 1023, 1023, POSITIVE_ONLY},
      {"atanh", myb_atanh, mpfr_atanh, -1022, -1, -1, BOTH_SIGNS},
      {"cbrt", myb_cbrt, mpfr_cbrt, -1022, 1023, 1023, BOTH_SIGNS},
  };
  for (size_t i = 0; i < sizeof unaries / sizeof unaries[0]; ++i) {
    test_unary(&unaries[i], 1000 + i);
  }
  for (size_t i = 0; i < 3; ++i) {
    test_reduction(&unaries[i]);
  }
  const binary binaries[] = {
      {"atan2", myb_atan2, mpfr_atan2},
      {"hypot", myb_hypot, mpfr_hypot},
  };
  for (size_t i = 0; i < sizeof binaries / sizeof binaries[0]; ++i) {
    test_binary(&binaries[i], 2000 + i);
  }
  return mt_finish("test_mpfr_cmath");
}
