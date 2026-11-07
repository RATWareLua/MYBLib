/* SPDX-License-Identifier: BSD-3-Clause */
/* myb_parse_double: syntax, consumed length, special values, and values checked against
 * Clang's own literals (which are correctly rounded). Exhaustive accuracy checks against MPFR
 * are in tests/mpfr/test_mpfr_parse.c. */
#include <float.h>
#include <locale.h>
#include <string.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/io.h"

static double parse(const char* text, size_t* consumed) {
  double value = -12345.0;
  size_t used = 0;
  const myb_status status = myb_parse_double(text, strlen(text), &value, &used);
  MT_CHECK(status == MYB_OK);
  if (consumed != NULL) {
    *consumed = used;
  }
  return value;
}

static void check_full(const char* text, double expected) {
  size_t used = 0;
  const double value = parse(text, &used);
  mtest.checks++;
  if (!mt_same(value, expected) || used != strlen(text)) {
    mtest.failures++;
    fprintf(stderr, "parse(\"%s\") = %a (consumed %zu), expected %a\n", text, value, used,
            expected);
  }
}

static void test_values(void) {
  check_full("0", 0.0);
  check_full("-0", -0.0);
  check_full("0.0e10", 0.0);
  check_full("-0.000e-999999999999", -0.0);
  check_full("1", 1.0);
  check_full("+1.5", 1.5);
  check_full(".5", 0.5);
  check_full("5.", 5.0);
  check_full("0.1", 0.1);
  check_full("3.141592653589793238462643383279502884197",
             3.141592653589793238462643383279502884197);
  check_full("1e23", 1e23); /* first power of ten that is not exact */
  check_full("8.98846567431158e307", 8.98846567431158e307);
  check_full("1.7976931348623157e308", DBL_MAX);
  check_full("1.7976931348623158e308", DBL_MAX); /* below DBL_MAX + ulp/2 */
  check_full("1.7976931348623159e308", (double)INFINITY);
  check_full("1e309", (double)INFINITY);
  check_full("-1e99999", -(double)INFINITY);
  check_full("2.2250738585072014e-308", DBL_MIN);
  check_full("2.2250738585072011e-308", 0x0.fffffffffffffp-1022); /* largest subnormal */
  check_full("4.9406564584124654e-324", 0x1p-1074);
  check_full("2.4703282292062328e-324", 0x1p-1074); /* above half of the smallest subnormal */
  check_full("2.4703282292062327e-324", 0.0);       /* below half */
  check_full("1e-400", 0.0);
  check_full("123456789012345678901234567890", 123456789012345678901234567890.0);
  check_full("0.000000000000000000000000000000000000000000001234",
             0.000000000000000000000000000000000000000000001234);
  /* Fortran 'D' exponents (NAIF kernels, IERS files). */
  check_full("1.0D0", 1.0);
  check_full("-6.55D-5", -6.55e-5);
  check_full("3.98600435507d+05", 3.98600435507e+05);
  /* Hexadecimal. */
  check_full("0x1.8p1", 3.0);
  check_full("-0x1p-1074", -0x1p-1074);
  check_full("0x1.fffffffffffff8p0", 2.0); /* tie, rounds to even */
  check_full("0x1.fffffffffffff7ffffp0", 0x1.fffffffffffffp0);
  check_full("0x.8", 0.5);
  /* Words. */
  check_full("inf", (double)INFINITY);
  check_full("-Infinity", -(double)INFINITY);
  size_t used = 0;
  const double nan = parse("NaN", &used);
  MT_CHECK(nan != nan && used == 3);
}

static void test_partial(void) {
  size_t used = 0;
  MT_CHECK_SAME(parse("1e", &used), 1.0); /* exponent without digits is not consumed */
  MT_CHECK(used == 1);
  MT_CHECK_SAME(parse("2.5e+x", &used), 2.5);
  MT_CHECK(used == 3);
  MT_CHECK_SAME(parse("7,5", &used), 7.0); /* ',' is never a decimal separator */
  MT_CHECK(used == 1);
  MT_CHECK_SAME(parse("0x", &used), 0.0); /* "0" then 'x' */
  MT_CHECK(used == 1);
  MT_CHECK_SAME(parse("12)", &used), 12.0);
  MT_CHECK(used == 2);
  MT_CHECK_SAME(parse("infinite", &used), (double)INFINITY);
  MT_CHECK(used == 3);
}

static void test_errors(void) {
  const char* bad[] = {"", "-", "+", ".", "e5", ".e5", "x", " 1", "--1", "+-1"};
  for (size_t i = 0; i < sizeof bad / sizeof bad[0]; ++i) {
    double value = 42.0;
    size_t used = 99;
    MT_CHECK(myb_parse_double(bad[i], strlen(bad[i]), &value, &used) == MYB_E_FORMAT);
    MT_CHECK(value == 42.0 && used == 99);
  }
  /* The length is respected: no reading past it. */
  double value = 0.0;
  size_t used = 0;
  MT_CHECK(myb_parse_double("123456", 3, &value, &used) == MYB_OK);
  MT_CHECK(value == 123.0 && used == 3);
}

/* The decimal separator does not follow the C locale. */
static void test_locale(void) {
  const char* names[] = {"de_DE.UTF-8", "de_DE", "German_Germany.1252", "ru_RU.UTF-8", "Russian"};
  const char* set = NULL;
  for (size_t i = 0; i < sizeof names / sizeof names[0] && set == NULL; ++i) {
    set = setlocale(LC_NUMERIC, names[i]);
  }
  if (set == NULL) {
    printf("no locale with a decimal comma available: locale check skipped\n");
    return;
  }
  check_full("1.25", 1.25);
  size_t used = 0;
  MT_CHECK_SAME(parse("1,25", &used), 1.0);
  MT_CHECK(used == 1);
  setlocale(LC_NUMERIC, "C");
}

int main(void) {
  test_values();
  test_partial();
  test_errors();
  test_locale();
  return mt_finish("test_parse");
}
