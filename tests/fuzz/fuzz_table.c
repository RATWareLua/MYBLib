/* SPDX-License-Identifier: BSD-3-Clause */
/* Tables of numbers (myb_table_load): the table read, written back (myb_table_write) and read
 * again must give the same numbers, bit for bit. */
#include <stdlib.h>

#include "fuzz.h"
#include "myb/io.h"

MYB_FUZZ_NO_INITIALIZE

static void check_round_trip(const myb_table* table) {
  const size_t rows = myb_table_rows(table);
  const size_t columns = myb_table_columns(table);
  /* a row of one NaN is written as an empty line, which a reader passes over */
  if (rows == 0 || columns < 2 || columns > 1024 || rows > 4096) {
    return;
  }
  double* values = malloc(rows * columns * sizeof(double));
  MYB_CHECK(values != NULL);
  for (size_t r = 0; r < rows; ++r) {
    memcpy(values + r * columns, myb_table_row(table, r), columns * sizeof(double));
  }
  size_t length = 0;
  MYB_CHECK(myb_table_write(NULL, columns, values, rows, NULL, 0, &length) == MYB_E_OVERFLOW);
  char* text = malloc(length + 1);
  MYB_CHECK(text != NULL);
  MYB_CHECK(myb_table_write(NULL, columns, values, rows, text, length + 1, &length) == MYB_OK);
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
  myb_table_destroy(back);
  free(text);
  free(values);
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
  for (size_t c = 0; c < columns && c < 64; ++c) {
    MYB_CHECK(strlen(myb_table_name(table, c)) <= 31);
  }
  check_round_trip(table);
  myb_table_destroy(table);
  return 0;
}
