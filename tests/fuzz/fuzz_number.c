/* SPDX-License-Identifier: BSD-3-Clause */
/* The number parser (myb_parse_double) at every place of the input where a number may start,
 * and the printers on what it reads: the shortest decimal and 17 significant digits must read
 * back to the same bits, the fixed form must fit its stated length. */
#include <math.h>

#include "fuzz.h"
#include "myb/io.h"

MYB_FUZZ_NO_INITIALIZE

static int same_bits(double a, double b) {
  uint64_t x;
  uint64_t y;
  memcpy(&x, &a, sizeof x);
  memcpy(&y, &b, sizeof y);
  return x == y || (a != a && b != b);
}

static void check_printers(double value, uint64_t hash) {
  char shortest[MYB_DOUBLE_CHARS];
  const size_t n = myb_format_double(value, shortest);
  MYB_CHECK(n > 0 && n < MYB_DOUBLE_CHARS && shortest[n] == '\0');
  double back = 0.0;
  size_t used = 0;
  MYB_CHECK(myb_parse_double(shortest, n, &back, &used) == MYB_OK && used == n);
  MYB_CHECK(same_bits(back, value));
  char text[512];
  size_t length = 0;
  if (myb_format_scientific(value, 17, text, sizeof text, &length) == MYB_OK && isfinite(value)) {
    MYB_CHECK(length == strlen(text));
    MYB_CHECK(myb_parse_double(text, length, &back, &used) == MYB_OK && used == length);
    MYB_CHECK(same_bits(back, value));
  }
  const int decimals = (int)(hash % 61u);
  const myb_status fixed = myb_format_fixed(value, decimals, text, sizeof text, &length);
  MYB_CHECK(fixed == MYB_OK || fixed == MYB_E_OVERFLOW);
  MYB_CHECK(fixed == MYB_E_OVERFLOW ? length + 1 > sizeof text : length == strlen(text));
  const int digits = 1 + (int)((hash >> 8) % 60u);
  const myb_status scientific = myb_format_scientific(value, digits, text, sizeof text, &length);
  MYB_CHECK(scientific == MYB_OK && length == strlen(text));
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  const char* text = (const char*)data;
  const uint64_t hash = fuzz_hash(data, size);
  size_t at = 0;
  int numbers = 0;
  /* bound: every iteration moves on by at least one character; at most 512 numbers */
  while (at < size && numbers < 512) {
    double value = 0.0;
    size_t used = 0;
    if (myb_parse_double(text + at, size - at, &value, &used) == MYB_OK) {
      MYB_CHECK(used >= 1 && used <= size - at);
      check_printers(value, hash ^ (uint64_t)numbers);
      at += used;
      ++numbers;
    } else {
      at += 1;
    }
  }
  return 0;
}
