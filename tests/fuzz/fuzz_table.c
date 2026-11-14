/* SPDX-License-Identifier: BSD-3-Clause */
/* Tables of numbers (myb_table_load): the table read, written back (myb_table_write, with and
 * without its names) and read again must give the same numbers, bit for bit, and the same
 * names. */
#include <stdlib.h>

#include "fuzz.h"
#include "myb/io.h"

MYB_FUZZ_NO_INITIALIZE

/* A name that the writer refuses as part of a row of names that are all numbers or empty. */
static int number_or_empty(const char* name) {
  const size_t n = strlen(name);
  double value = 0.0;
  size_t used = 0;
  return n == 0 || (myb_parse_double(name, n, &value, &used) == MYB_OK && used == n);
}

static void check_round_trip(const myb_table* table, int with_names) {
  const size_t rows = myb_table_rows(table);
  const size_t columns = myb_table_columns(table);
  if (columns > 1024 || rows > 4096 || (rows == 0 && !with_names)) {
    return;
  }
  const char** names = NULL;
  int numeric = 1;
  if (with_names) {
    names = malloc(columns * sizeof *names);
    MYB_CHECK(names != NULL);
    for (size_t c = 0; c < columns; ++c) {
      names[c] = myb_table_name(table, c);
      numeric = numeric && number_or_empty(names[c]);
    }
  }
  double* values = malloc((rows > 0 ? rows : 1) * columns * sizeof(double));
  MYB_CHECK(values != NULL);
  for (size_t r = 0; r < rows; ++r) {
    memcpy(values + r * columns, myb_table_row(table, r), columns * sizeof(double));
  }
  size_t length = 0;
  const myb_status sized = myb_table_write(names, columns, values, rows, NULL, 0, &length);
  if (with_names && numeric) {
    /* names cut to 31 characters that became numbers: the writer refuses them */
    MYB_CHECK(sized == MYB_E_DOMAIN);
    free(values);
    free((void*)names);
    return;
  }
  MYB_CHECK(sized == MYB_E_OVERFLOW);
  char* text = malloc(length + 1);
  MYB_CHECK(text != NULL);
  MYB_CHECK(myb_table_write(names, columns, values, rows, text, length + 1, &length) == MYB_OK);
  myb_table* back = NULL;
  myb_diag diag = {0};
  MYB_CHECK(myb_table_load(text, length, NULL, &back, &diag) == MYB_OK);
  MYB_CHECK(myb_table_rows(back) == rows && myb_table_columns(back) == columns);
  for (size_t r = 0; r < rows; ++r) {
    const double* row = myb_table_row(back, r);
    for (size_t c = 0; c < columns; ++c) {
      const double v = values[r * columns + c];
      MYB_CHECK(memcmp(&row[c], &v, sizeof v) == 0 || (v != v && row[c] != row[c]));
    }
  }
  for (size_t c = 0; with_names && c < columns; ++c) {
    MYB_CHECK(strcmp(myb_table_name(back, c), names[c]) == 0);
  }
  myb_table_destroy(back);
  free(text);
  free(values);
  free((void*)names);
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_table* table = NULL;
  myb_diag diag = {0};
  if (myb_table_load((const char*)data, size, NULL, &table, &diag) != MYB_OK) {
    MYB_CHECK(table == NULL);
    return 0;
  }
  const size_t columns = myb_table_columns(table);
  MYB_CHECK(columns > 0);
  int named = 0;
  for (size_t c = 0; c < columns; ++c) {
    MYB_CHECK(strlen(myb_table_name(table, c)) <= 31);
    named = named || myb_table_name(table, c)[0] != '\0';
  }
  check_round_trip(table, 0);
  if (named) {
    check_round_trip(table, 1);
  }
  myb_table_destroy(table);
  return 0;
}
