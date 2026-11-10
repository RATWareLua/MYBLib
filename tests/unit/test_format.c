/* SPDX-License-Identifier: BSD-3-Clause */
/* Printing of numbers (myb_format_double, myb_format_fixed, myb_format_scientific):
 *  - against CPython's formatting (David Gay's dtoa: correctly rounded, the shortest repr) on
 *    6 100 doubles — random bit patterns, powers of two and their neighbours, subnormals,
 *    decimal fractions, exact ties (tests/data/io/format.ref);
 *  - every value printed in the shortest form reads back to the same bits;
 *  - special values, buffer sizes, errors. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/io.h"

static void test_reference(void) {
  FILE* file = fopen(MYB_TEST_DATA_DIR "/io/format.ref", "rb");
  MT_CHECK(file != NULL);
  if (file == NULL) {
    return;
  }
  static char line[2048];
  long cases = 0;
  long wrong_short = 0;
  long wrong_fixed = 0;
  long wrong_scientific = 0;
  /* bound: the lines of the file */
  while (fgets(line, sizeof line, file) != NULL) {
    if (line[0] == '#') {
      continue;
    }
    char* field[6];
    int count = 0;
    char* p = line;
    /* bound: six fields */
    while (count < 6 && *p != '\0' && *p != '\n') {
      field[count++] = p;
      while (*p != ' ' && *p != '\n' && *p != '\0') {
        ++p;
      }
      if (*p != '\0') {
        *p++ = '\0';
      }
    }
    MT_CHECK(count == 6);
    if (count != 6) {
      continue;
    }
    const double value = mt_from_bits(strtoull(field[0], NULL, 16));
    ++cases;
    char text[MYB_DOUBLE_CHARS];
    const size_t n = myb_format_double(value, text);
    if (n != strlen(field[1]) || strcmp(text, field[1]) != 0) {
      if (wrong_short++ < 5) {
        fprintf(stderr, "shortest: %s, expected %s\n", text, field[1]);
      }
    }
    double back = 0.0;
    size_t used = 0;
    MT_CHECK(myb_parse_double(text, n, &back, &used) == MYB_OK && used == n);
    MT_CHECK(mt_bits(back) == mt_bits(value));
    static char wide[512];
    size_t length = 0;
    if (strcmp(field[3], "-") != 0) {
      MT_CHECK(myb_format_fixed(value, atoi(field[2]), wide, sizeof wide, &length) == MYB_OK);
      if (length != strlen(field[3]) || strcmp(wide, field[3]) != 0) {
        if (wrong_fixed++ < 5) {
          fprintf(stderr, "fixed %s: %s, expected %s\n", field[2], wide, field[3]);
        }
      }
    }
    MT_CHECK(myb_format_scientific(value, atoi(field[4]), wide, sizeof wide, &length) == MYB_OK);
    if (length != strlen(field[5]) || strcmp(wide, field[5]) != 0) {
      if (wrong_scientific++ < 5) {
        fprintf(stderr, "scientific %s: %s, expected %s\n", field[4], wide, field[5]);
      }
    }
  }
  fclose(file);
  printf("printing against CPython: %ld values; differences: shortest %ld, fixed %ld, "
         "scientific %ld\n",
         cases, wrong_short, wrong_fixed, wrong_scientific);
  MT_CHECK(cases > 6000);
  MT_CHECK(wrong_short == 0);
  MT_CHECK(wrong_fixed == 0);
  MT_CHECK(wrong_scientific == 0);
}

static void test_round_trip(void) {
  mt_rng rng = {.state = 2026};
  long worst_length = 0;
  for (int k = 0; k < 200000; ++k) {
    const double value = mt_any_double(&rng);
    if (value != value || value - value != 0.0) {
      continue;
    }
    char text[MYB_DOUBLE_CHARS];
    const size_t n = myb_format_double(value, text);
    double back = 0.0;
    size_t used = 0;
    MT_CHECK(myb_parse_double(text, n, &back, &used) == MYB_OK && used == n);
    MT_CHECK(mt_bits(back) == mt_bits(value));
    worst_length = (long)n > worst_length ? (long)n : worst_length;
  }
  printf("shortest form: 200 000 random doubles read back to the same bits; longest text %ld "
         "characters\n",
         worst_length);
  MT_CHECK(worst_length <= 24);
}

static void test_special(void) {
  char text[MYB_DOUBLE_CHARS];
  MT_CHECK(myb_format_double(0.0, text) == 3 && strcmp(text, "0.0") == 0);
  MT_CHECK(myb_format_double(-0.0, text) == 4 && strcmp(text, "-0.0") == 0);
  MT_CHECK(myb_format_double((double)NAN, text) == 3 && strcmp(text, "nan") == 0);
  MT_CHECK(myb_format_double((double)INFINITY, text) == 3 && strcmp(text, "inf") == 0);
  MT_CHECK(myb_format_double(-(double)INFINITY, text) == 4 && strcmp(text, "-inf") == 0);
  MT_CHECK(myb_format_double(0.1, text) == 3 && strcmp(text, "0.1") == 0);
  MT_CHECK(myb_format_double(1e16, text) > 0 && strcmp(text, "1e+16") == 0);
  MT_CHECK(myb_format_double(123456789012345.0, text) > 0 &&
           strcmp(text, "123456789012345.0") == 0);
  MT_CHECK(myb_format_double(5e-324, text) > 0 && strcmp(text, "5e-324") == 0);
  MT_CHECK(myb_format_double(-1.5e-7, text) > 0 && strcmp(text, "-1.5e-07") == 0);

  char out[64];
  size_t length = 0;
  MT_CHECK(myb_format_fixed(2.5, 0, out, sizeof out, &length) == MYB_OK &&
           strcmp(out, "2") == 0); /* a tie goes to the even digit */
  MT_CHECK(myb_format_fixed(3.5, 0, out, sizeof out, &length) == MYB_OK && strcmp(out, "4") == 0);
  MT_CHECK(myb_format_fixed(0.125, 2, out, sizeof out, &length) == MYB_OK &&
           strcmp(out, "0.12") == 0);
  MT_CHECK(myb_format_fixed(1.005, 2, out, sizeof out, &length) == MYB_OK &&
           strcmp(out, "1.00") == 0); /* 1.005 is below 1.005: 1.00499999999999989... */
  MT_CHECK(myb_format_fixed(-0.001, 2, out, sizeof out, &length) == MYB_OK &&
           strcmp(out, "-0.00") == 0);
  MT_CHECK(myb_format_fixed(999.9996, 3, out, sizeof out, &length) == MYB_OK &&
           strcmp(out, "1000.000") == 0 && length == 8);
  MT_CHECK(myb_format_fixed((double)NAN, 3, out, sizeof out, &length) == MYB_OK &&
           strcmp(out, "nan") == 0);
  MT_CHECK(myb_format_scientific(0.0, 4, out, sizeof out, &length) == MYB_OK &&
           strcmp(out, "0.000e+00") == 0);
  MT_CHECK(myb_format_scientific(9.9996, 4, out, sizeof out, &length) == MYB_OK &&
           strcmp(out, "1.000e+01") == 0);
  MT_CHECK(myb_format_scientific(-299792458.0, 6, out, sizeof out, &length) == MYB_OK &&
           strcmp(out, "-2.99792e+08") == 0);
  MT_CHECK(myb_format_scientific(1e-310, 3, out, sizeof out, &length) == MYB_OK &&
           strcmp(out, "1.00e-310") == 0);
  /* sizes: the full length is reported, the text is cut and terminated */
  char small[5];
  MT_CHECK(myb_format_fixed(123.456, 3, small, sizeof small, &length) == MYB_E_OVERFLOW);
  MT_CHECK(length == 7 && strcmp(small, "123.") == 0);
  MT_CHECK(myb_format_fixed(123.456, 3, NULL, 0, &length) == MYB_E_OVERFLOW && length == 7);
  MT_CHECK(myb_format_fixed(1.0, -1, out, sizeof out, &length) == MYB_E_DOMAIN);
  MT_CHECK(myb_format_fixed(1.0, 61, out, sizeof out, &length) == MYB_E_DOMAIN);
  MT_CHECK(myb_format_scientific(1.0, 0, out, sizeof out, &length) == MYB_E_DOMAIN);
  MT_CHECK(myb_format_fixed(1.7976931348623157e308, 60, out, sizeof out, &length) ==
           MYB_E_OVERFLOW);
  MT_CHECK(length == 309 + 1 + 60);
}

int main(void) {
  test_reference();
  test_round_trip();
  test_special();
  return mt_finish("test_format");
}
