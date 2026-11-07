/* SPDX-License-Identifier: BSD-3-Clause */
/* myb_parse_double is correctly rounded: compared with mpfr_strtofr (exact for any length) on
 * random short and long decimals, on hexadecimal strings, and on the hardest inputs — exact
 * midpoints between neighbouring doubles and strings just above or below them, including
 * strings longer than the parser's 800-digit window. */
#include <mpfr.h>
#include <stdlib.h>
#include <string.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/io.h"

enum { BUFFER = 2048 };

/* Correctly rounded binary64 value of text: MPFR with the exponent range of binary64 and
 * subnormalisation (the method of CORE-MATH's checks). The range is restored afterwards,
 * because the midpoint construction needs the full MPFR range. */
static double reference(const char* text) {
  const mpfr_exp_t saved_emin = mpfr_get_emin();
  const mpfr_exp_t saved_emax = mpfr_get_emax();
  mpfr_set_emin(-1073);
  mpfr_set_emax(1024);
  mpfr_t x;
  mpfr_init2(x, 53);
  const int inexact = mpfr_strtofr(x, text, NULL, 0, MPFR_RNDN);
  mpfr_subnormalize(x, inexact, MPFR_RNDN);
  const double r = mpfr_get_d(x, MPFR_RNDN);
  mpfr_clear(x);
  mpfr_set_emin(saved_emin);
  mpfr_set_emax(saved_emax);
  return r;
}

static long wrong_total = 0;

static void check(const char* label, const char* text) {
  double value = 0.0;
  size_t used = 0;
  const myb_status status = myb_parse_double(text, strlen(text), &value, &used);
  const double expected = reference(text);
  mtest.checks++;
  if (status != MYB_OK || used != strlen(text) || !mt_same(value, expected)) {
    mtest.failures++;
    if (++wrong_total <= 10) {
      fprintf(stderr, "%s: \"%.80s%s\" -> %a, MPFR %a\n", label, text,
              strlen(text) > 80 ? "..." : "", value, expected);
    }
  }
}

/* Random decimal with the given number of significant digits and decimal exponent range. */
static void random_decimal(mt_rng* rng, int digits, int emin, int emax, char* out) {
  size_t n = 0;
  if ((mt_next(rng) & 1) != 0) {
    out[n++] = '-';
  }
  out[n++] = (char)('1' + (int)(mt_next(rng) % 9));
  out[n++] = '.';
  for (int i = 1; i < digits; ++i) {
    out[n++] = (char)('0' + (int)(mt_next(rng) % 10));
  }
  snprintf(out + n, (size_t)BUFFER - n, "e%d", mt_int(rng, emin, emax));
}

static void test_random_decimals(void) {
  mt_rng rng = {.state = 77};
  static char text[BUFFER];
  long count = 0;
  for (int i = 0; i < 300000; ++i) {
    random_decimal(&rng, mt_int(&rng, 1, 20), -345, 310, text);
    check("short", text);
    ++count;
  }
  for (int i = 0; i < 3000; ++i) {
    random_decimal(&rng, mt_int(&rng, 20, 1000), -345, 310, text);
    check("long", text);
    ++count;
  }
  printf("random decimals: %ld strings\n", count);
}

static void test_random_hex(void) {
  mt_rng rng = {.state = 78};
  static const char digits[] = "0123456789abcdef";
  static char text[BUFFER];
  for (int i = 0; i < 100000; ++i) {
    size_t n = 0;
    text[n++] = '0';
    text[n++] = 'x';
    text[n++] = digits[1 + mt_next(&rng) % 15];
    text[n++] = '.';
    const int count = mt_int(&rng, 1, 30);
    for (int j = 0; j < count; ++j) {
      text[n++] = digits[mt_next(&rng) % 16];
    }
    snprintf(text + n, sizeof text - n, "p%d", mt_int(&rng, -1100, 1030));
    check("hex", text);
  }
  printf("random hex: 100000 strings\n");
}

/* digits: decimal digit string (no point), as from mpfr_get_str. Adds or subtracts one unit in
 * the last place, in place; the string keeps its length (no carry beyond the first digit is
 * needed for the values used here). */
static void step_last_digit(char* digits, int up) {
  /* bound: the string length */
  for (size_t i = strlen(digits); i > 0; --i) {
    char* c = &digits[i - 1];
    if (up) {
      if (*c != '9') {
        *c += 1;
        return;
      }
      *c = '0';
    } else {
      if (*c != '0') {
        *c -= 1;
        return;
      }
      *c = '9';
    }
  }
}

/* "0.<digits>e<exp>" for MPFR's convention value = 0.digits * 10^exp. */
static void format_scientific(const char* digits, long exponent, const char* suffix, char* out) {
  snprintf(out, BUFFER, "0.%s%se%ld", digits, suffix, exponent);
}

static void test_midpoints(void) {
  mt_rng rng = {.state = 79};
  mpfr_t x;
  mpfr_t next;
  mpfr_t mid;
  mpfr_inits2(200, x, next, mid, (mpfr_ptr)0);
  static char text[BUFFER];
  long count = 0;
  for (int i = 0; i < 4000; ++i) {
    double a;
    const int kind = i % 4;
    if (kind == 0) {
      a = mt_from_bits(mt_next(&rng) & UINT64_C(0x000FFFFFFFFFFFFF)); /* subnormal */
    } else if (kind == 1) {
      a = mt_double(&rng, 1000, 1023); /* near overflow */
    } else {
      a = mt_double(&rng, -1022, 1022);
    }
    a = a < 0.0 ? -a : a;
    double b = mt_from_bits(mt_bits(a) + 1);
    if (b == (double)INFINITY) {
      b = a;
      a = mt_from_bits(mt_bits(a) - 1);
    }
    mpfr_set_d(x, a, MPFR_RNDN);
    mpfr_set_d(next, b, MPFR_RNDN);
    mpfr_add(mid, x, next, MPFR_RNDN); /* exact at 200 bits */
    mpfr_div_2ui(mid, mid, 1, MPFR_RNDN);
    mpfr_exp_t exponent = 0;
    char* digits = mpfr_get_str(NULL, &exponent, 10, 790, mid, MPFR_RNDN); /* exact */
    char variant[BUFFER];

    format_scientific(digits, (long)exponent, "", text);
    check("midpoint", text);
    format_scientific(digits, (long)exponent, "0000000001", text); /* beyond 800 digits */
    check("midpoint+tail", text);

    strcpy(variant, digits);
    step_last_digit(variant, 1);
    format_scientific(variant, (long)exponent, "", text);
    check("midpoint+1", text);
    strcpy(variant, digits);
    step_last_digit(variant, 0);
    format_scientific(variant, (long)exponent, "", text);
    check("midpoint-1", text);
    mpfr_free_str(digits);
    count += 4;
  }
  mpfr_clears(x, next, mid, (mpfr_ptr)0);
  printf("midpoint cases: %ld strings\n", count);
}

int main(void) {
  test_random_decimals();
  test_random_hex();
  test_midpoints();
  return mt_finish("test_mpfr_parse");
}
