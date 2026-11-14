/* SPDX-License-Identifier: BSD-3-Clause */
/* NAIF text kernels (myb_kpool_load_text): the text loaded twice into one pool (appends act
 * on what the first load left, a failed load leaves the pool as it was), then every
 * variable named in it asked for its type, its numbers and its strings. */
#include "fuzz.h"
#include "myb/io.h"

MYB_FUZZ_NO_INITIALIZE

static int is_name_char(uint8_t c) {
  return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
         c == '_' || c == '/' || c == '-' || c == '.';
}

/* Asks the pool for each word of the text that could be a variable name (at most 128). */
static void ask_names(const myb_kpool* pool, const uint8_t* data, size_t size) {
  size_t at = 0;
  int asked = 0;
  /* bound: every iteration moves on by at least one character */
  while (at < size && asked < 128) {
    if (!is_name_char(data[at])) {
      ++at;
      continue;
    }
    size_t end = at;
    /* bound: the rest of the text */
    while (end < size && is_name_char(data[end])) {
      ++end;
    }
    char name[40];
    if (end - at < sizeof name) {
      memcpy(name, data + at, end - at);
      name[end - at] = '\0';
      myb_kvar_type type = 0;
      size_t count = 0;
      if (myb_kpool_describe(pool, name, &type, &count) == MYB_OK) {
        MYB_CHECK(type == MYB_KVAR_NUMBERS || type == MYB_KVAR_STRINGS);
        MYB_CHECK(count > 0);
        const double* values = NULL;
        size_t n = 0;
        const myb_status numbers = myb_kpool_numbers(pool, name, &values, &n);
        MYB_CHECK((numbers == MYB_OK) == (type == MYB_KVAR_NUMBERS));
        double sum = 0.0;
        for (size_t i = 0; numbers == MYB_OK && i < n; ++i) {
          sum += values[i];
        }
        (void)sum;
        for (size_t i = 0; i <= count && i < 64; ++i) {
          const char* value = NULL;
          size_t length = 0;
          const myb_status string = myb_kpool_string(pool, name, i, &value, &length);
          if (string == MYB_OK) {
            MYB_CHECK(type == MYB_KVAR_STRINGS && i < count);
            for (size_t j = 0; j < length; ++j) {
              sum += (double)value[j];
            }
          }
        }
      }
      ++asked;
    }
    at = end;
  }
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_kpool* pool = NULL;
  MYB_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  myb_diag diag = {0};
  const myb_status first = myb_kpool_load_text(pool, (const char*)data, size, &diag);
  if (first == MYB_OK) {
    ask_names(pool, data, size);
    /* the second load appends to what the first one left, or fails and changes nothing */
    (void)myb_kpool_load_text(pool, (const char*)data, size, &diag);
  }
  ask_names(pool, data, size);
  myb_kpool_destroy(pool);
  return 0;
}
