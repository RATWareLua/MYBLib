/* SPDX-License-Identifier: BSD-3-Clause */
/* Tables of numbers as text (myb_table_*): separators, names, comments, empty fields; a table
 * written and read back gives every number to the bit; errors with their lines. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/io.h"

static myb_status load(const char* text, myb_table** table, myb_diag* diag) {
  return myb_table_load(text, strlen(text), NULL, table, diag);
}

static void test_reading(void) {
  static const char text[] = "# a table\n"
                             "t, x ,a_name_that_is_longer_than_the_31_characters_kept\n"
                             "0,1.5,-2e3\n"
                             "1,,3\n"
                             "\n"
                             "2;4;5 # the rest is a comment\n"
                             "3\t6\t7\r\n"
                             "  4 8   9  \n"
                             "5,10,\n"
                             ",,\n"
                             "6, 0x1p-3 ,-1e400"; /* what myb_parse_double reads; no last newline */
  myb_table* table = NULL;
  myb_diag diag = {0};
  MT_CHECK(load(text, &table, &diag) == MYB_OK);
  if (table == NULL) {
    printf("line %ld: %s\n", diag.line, diag.message);
    return;
  }
  MT_CHECK(myb_table_rows(table) == 8 && myb_table_columns(table) == 3);
  MT_CHECK(myb_table_row(table, 7)[0] == 6.0 && myb_table_row(table, 7)[1] == 0.125 &&
           myb_table_row(table, 7)[2] == -(double)INFINITY);
  MT_CHECK(strcmp(myb_table_name(table, 0), "t") == 0);
  MT_CHECK(strcmp(myb_table_name(table, 1), "x") == 0);
  MT_CHECK(strcmp(myb_table_name(table, 2), "a_name_that_is_longer_than_the_") == 0);
  static const double expected[7][3] = {{0, 1.5, -2e3}, {1, 0, 3}, {2, 4, 5}, {3, 6, 7},
                                        {4, 8, 9},      {5, 10, 0}, {0, 0, 0}};
  static const int missing[7][3] = {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}, {0, 0, 0},
                                    {0, 0, 0}, {0, 0, 1}, {1, 1, 1}};
  for (size_t r = 0; r < 7; ++r) {
    const double* row = myb_table_row(table, r);
    for (size_t c = 0; c < 3; ++c) {
      if (missing[r][c]) {
        MT_CHECK(row[c] != row[c]);
      } else {
        MT_CHECK(row[c] == expected[r][c]);
      }
    }
  }
  myb_table_destroy(table);
  myb_table_destroy(NULL);

  /* no names */
  table = NULL;
  MT_CHECK(load("1 2\n3 4\n", &table, &diag) == MYB_OK);
  if (table != NULL) {
    MT_CHECK(myb_table_rows(table) == 2 && myb_table_columns(table) == 2);
    MT_CHECK(myb_table_name(table, 1)[0] == '\0' && myb_table_row(table, 1)[0] == 3.0);
    myb_table_destroy(table);
  }
  /* names alone are a table without rows */
  table = NULL;
  MT_CHECK(load("a,b,c\n", &table, &diag) == MYB_OK);
  if (table != NULL) {
    MT_CHECK(myb_table_rows(table) == 0 && myb_table_columns(table) == 3);
    MT_CHECK(strcmp(myb_table_name(table, 2), "c") == 0);
    myb_table_destroy(table);
  }
  /* errors */
  table = NULL;
  MT_CHECK(load("", &table, &diag) == MYB_E_FORMAT && diag.message[0] != '\0');
  MT_CHECK(load("# nothing\n\n", &table, &diag) == MYB_E_FORMAT);
  MT_CHECK(load("1,2\n3\n", &table, &diag) == MYB_E_FORMAT && diag.line == 2);
  MT_CHECK(load("a,b\n1,2\nc,d\n", &table, &diag) == MYB_E_FORMAT && diag.line == 3);
  MT_CHECK(load("1,2\n3,x\n", &table, &diag) == MYB_E_FORMAT && diag.line == 2);
  MT_CHECK(load("a,b\n1,2,3\n", &table, &diag) == MYB_E_FORMAT && diag.line == 2);
  MT_CHECK(table == NULL);
  MT_CHECK(myb_table_load_file("no-such-file.csv", NULL, &table, &diag) == MYB_E_IO);
}

static uint64_t next_random(uint64_t* state) {
  *state = *state * 6364136223846793005ull + 1442695040888963407ull;
  uint64_t z = *state;
  z ^= z >> 33;
  z *= 0xff51afd7ed558ccdull;
  z ^= z >> 33;
  return z;
}

static void test_round_trip(void) {
  enum { ROWS = 500, COLUMNS = 7 };
  static double values[ROWS * COLUMNS];
  uint64_t state = 20261110;
  for (size_t k = 0; k < ROWS * COLUMNS; ++k) {
    double v = 0.0;
    /* bound: finite doubles are 2046/2048 of all bit patterns */
    do {
      const uint64_t bits = next_random(&state);
      memcpy(&v, &bits, sizeof v);
    } while (!isfinite(v));
    values[k] = v;
  }
  values[3] = (double)NAN;                 /* inside a row */
  values[COLUMNS - 1] = (double)NAN;       /* at the end of a row */
  values[2 * COLUMNS] = (double)NAN;       /* at the start */
  values[5 * COLUMNS + 1] = 0.0;
  values[5 * COLUMNS + 2] = -0.0;
  values[5 * COLUMNS + 3] = 5e-324;
  values[5 * COLUMNS + 4] = 1.7976931348623157e308;
  static const char* const names[COLUMNS] = {"t", "x", "y", "z", "vx", "vy", "vz"};
  size_t length = 0;
  MT_CHECK(myb_table_write(names, COLUMNS, values, ROWS, NULL, 0, &length) == MYB_E_OVERFLOW);
  char* text = malloc(length + 1);
  MT_CHECK(text != NULL);
  if (text == NULL) {
    return;
  }
  MT_CHECK(myb_table_write(names, COLUMNS, values, ROWS, text, length + 1, NULL) == MYB_OK);
  MT_CHECK(strlen(text) == length && memcmp(text, "t,x,y,z,vx,vy,vz\n", 17) == 0);
  myb_table* table = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_table_load(text, length, NULL, &table, &diag) == MYB_OK);
  if (table != NULL) {
    MT_CHECK(myb_table_rows(table) == ROWS && myb_table_columns(table) == COLUMNS);
    MT_CHECK(strcmp(myb_table_name(table, 4), "vx") == 0);
    size_t differences = 0;
    for (size_t r = 0; r < ROWS; ++r) {
      const double* row = myb_table_row(table, r);
      for (size_t c = 0; c < COLUMNS; ++c) {
        const double v = values[r * COLUMNS + c];
        if (v != v ? row[c] == row[c] : memcmp(&row[c], &v, sizeof v) != 0) {
          ++differences;
        }
      }
    }
    printf("table of %d x %d random doubles written and read back: %zu differences\n", ROWS,
           COLUMNS, differences);
    MT_CHECK(differences == 0);
    myb_table_destroy(table);
  } else {
    printf("line %ld: %s\n", diag.line, diag.message);
  }
  free(text);
  /* without names; invalid descriptions */
  char small[64];
  static const double two[4] = {1.0, 2.5, -3.0, 1e-7};
  MT_CHECK(myb_table_write(NULL, 2, two, 2, small, sizeof small, &length) == MYB_OK);
  MT_CHECK(strcmp(small, "1.0,2.5\n-3.0,1e-07\n") == 0 && length == 19);
  MT_CHECK(myb_table_write(NULL, 0, two, 2, small, sizeof small, &length) == MYB_E_DOMAIN);
  static const char* const hole[2] = {"a", NULL};
  MT_CHECK(myb_table_write(hole, 2, two, 2, small, sizeof small, &length) == MYB_E_DOMAIN);
  MT_CHECK(myb_table_write(NULL, 2, two, 2, small, 10, &length) == MYB_E_OVERFLOW && length == 19);
}

int main(void) {
  test_reading();
  test_round_trip();
  return mt_finish("test_table");
}
