/* SPDX-License-Identifier: BSD-3-Clause */
/* myb_sin, myb_cos are correctly rounded: compared bit for bit with GNU MPFR on random inputs
 * of every magnitude, on hard argument reductions near multiples of pi/2, and on a sample of
 * the CORE-MATH worst cases for correct rounding (tests/data/core-math). */
#include <float.h>
#include <mpfr.h>
#include <stdlib.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/cmath.h"

typedef double (*function)(double);
typedef int (*mpfr_function)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);

typedef struct subject {
  const char* name;
  function under_test;
  mpfr_function reference;
} subject;

/* RN(f(x)) for binary64, with subnormal results rounded like binary64 (emin/emax of double
 * are set in main, as in CORE-MATH's own checks). */
static double reference_value(mpfr_function f, double x) {
  mpfr_t y;
  mpfr_init2(y, 53);
  mpfr_set_d(y, x, MPFR_RNDN);
  const int inexact = f(y, y, MPFR_RNDN);
  mpfr_subnormalize(y, inexact, MPFR_RNDN);
  const double result = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear(y);
  return result;
}

static long check_one(const subject* s, double x, long* shown) {
  const double expected = reference_value(s->reference, x);
  const double actual = s->under_test(x);
  mtest.checks++;
  if (mt_same(actual, expected)) {
    return 0;
  }
  mtest.failures++;
  if (++*shown <= 5) {
    fprintf(stderr, "%s(%a) = %a, MPFR %a\n", s->name, x, actual, expected);
  }
  return 1;
}

static void check_range(const subject* s, const char* label, int emin, int emax, long count,
                        uint64_t seed) {
  mt_rng rng = {.state = seed};
  long wrong = 0;
  long shown = 0;
  for (long i = 0; i < count; ++i) {
    wrong += check_one(s, mt_double(&rng, emin, emax), &shown);
  }
  printf("%s %-10s %7ld inputs, %ld not correctly rounded\n", s->name, label, count, wrong);
}

static void check_subnormals(const subject* s, long count, uint64_t seed) {
  mt_rng rng = {.state = seed};
  long wrong = 0;
  long shown = 0;
  for (long i = 0; i < count; ++i) {
    wrong += check_one(s, mt_from_bits(mt_next(&rng) & UINT64_C(0x800FFFFFFFFFFFFF)), &shown);
  }
  printf("%s %-10s %7ld inputs, %ld not correctly rounded\n", s->name, "subnormal", count, wrong);
}

/* Doubles nearest to k*pi/2 and their neighbours: the hardest argument reductions. */
static void check_near_half_pi_multiples(const subject* s) {
  long wrong = 0;
  long shown = 0;
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
      wrong += check_one(s, x, &shown);
      wrong += check_one(s, -x, &shown);
      count += 2;
    }
  }
  mpfr_clear(t);
  printf("%s %-10s %7ld inputs, %ld not correctly rounded\n", s->name, "k*pi/2", count, wrong);
}

static void check_specials(const subject* s) {
  const double values[] = {
      0.0,        -0.0,     DBL_MIN, -DBL_MIN, 0x1p-1074,        -0x1p-1074,
      DBL_MAX,    -DBL_MAX, 0x1p-26, 0x1p-27,  0x1p31,           0x1p31 * 3.0,
      1.0,        -1.0,     0.5,     2.0,      (double)INFINITY, -(double)INFINITY,
      (double)NAN};
  long wrong = 0;
  long shown = 0;
  for (size_t i = 0; i < sizeof values / sizeof values[0]; ++i) {
    wrong += check_one(s, values[i], &shown);
  }
  printf("%s %-10s %7zu inputs, %ld not correctly rounded\n", s->name, "special",
         sizeof values / sizeof values[0], wrong);
}

/* One hexadecimal double per line, '#' starts a comment line. Both x and -x are checked. */
static void check_worst_cases(const subject* s, const char* file_name) {
  char path[1024];
  const int n = snprintf(path, sizeof path, "%s/core-math/%s", MYB_TEST_DATA_DIR, file_name);
  MT_CHECK(n > 0 && (size_t)n < sizeof path);
  FILE* f = fopen(path, "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  long wrong = 0;
  long shown = 0;
  long count = 0;
  char line[256];
  /* bound: one iteration per line of a finite file */
  while (fgets(line, sizeof line, f) != NULL) {
    if (line[0] == '#' || line[0] == '\n') {
      continue;
    }
    char* end = NULL;
    const double x = strtod(line, &end);
    MT_CHECK(end != line);
    wrong += check_one(s, x, &shown);
    wrong += check_one(s, -x, &shown);
    count += 2;
  }
  fclose(f);
  MT_CHECK(count > 0);
  printf("%s %-10s %7ld inputs, %ld not correctly rounded (%s)\n", s->name, "worst", count, wrong,
         file_name);
}

static void check_subject(const subject* s, const char* worst_file, uint64_t seed) {
  check_specials(s);
  check_subnormals(s, 20000, seed + 1);
  check_range(s, "tiny", -1022, -27, 100000, seed + 2);
  check_range(s, "small", -26, 0, 200000, seed + 3);
  check_range(s, "moderate", 1, 30, 200000, seed + 4);
  check_range(s, "large", 31, 1023, 50000, seed + 5);
  check_near_half_pi_multiples(s);
  check_worst_cases(s, worst_file);
}

int main(void) {
  mpfr_set_emin(-1073);
  mpfr_set_emax(1024);
  const subject sin_subject = {"sin", myb_sin, mpfr_sin};
  const subject cos_subject = {"cos", myb_cos, mpfr_cos};
  check_subject(&sin_subject, "sin.wc", 100);
  check_subject(&cos_subject, "cos.wc", 200);
  return mt_finish("test_mpfr_cmath");
}
