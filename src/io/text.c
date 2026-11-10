/* SPDX-License-Identifier: BSD-3-Clause */
#include <string.h>

#include "core/base.h"
#include "io/text.h"
#include "myb/io.h"

enum { NUMBER_CHARS = 448 }; /* the longest fixed form: 309 digits, the point, 60 decimals */

void myb_internal_text_char(myb_internal_text* t, char c) {
  if (t->length + 1 < t->size) {
    t->out[t->length] = c;
  }
  t->length++;
}

void myb_internal_text_span(myb_internal_text* t, const char* s, size_t length) {
  for (size_t i = 0; i < length; ++i) {
    myb_internal_text_char(t, s[i]);
  }
}

void myb_internal_text_put(myb_internal_text* t, const char* s) {
  myb_internal_text_span(t, s, strlen(s));
}

void myb_internal_text_fill(myb_internal_text* t, char c, size_t count) {
  for (size_t i = 0; i < count; ++i) {
    myb_internal_text_char(t, c);
  }
}

static void aligned(myb_internal_text* t, const char* s, size_t length, int width) {
  if (width > 0 && length < (size_t)width) {
    myb_internal_text_fill(t, ' ', (size_t)width - length);
  }
  myb_internal_text_span(t, s, length);
}

void myb_internal_text_int(myb_internal_text* t, int64_t value, int width, char pad) {
  char digits[24];
  int count = 0;
  const int negative = value < 0;
  uint64_t magnitude = negative ? (uint64_t)0 - (uint64_t)value : (uint64_t)value;
  /* bound: 20 decimal digits */
  do {
    digits[count++] = (char)('0' + magnitude % 10);
    magnitude /= 10;
  } while (magnitude != 0 && count < 20);
  const int length = count + negative;
  if (negative && pad == '0') {
    myb_internal_text_char(t, '-');
  }
  for (int i = length; i < width; ++i) {
    myb_internal_text_char(t, pad);
  }
  if (negative && pad != '0') {
    myb_internal_text_char(t, '-');
  }
  for (int i = count - 1; i >= 0; --i) {
    myb_internal_text_char(t, digits[i]);
  }
}

void myb_internal_text_double(myb_internal_text* t, double value) {
  char number[MYB_DOUBLE_CHARS];
  const size_t length = myb_format_double(value, number);
  myb_internal_text_span(t, number, length);
}

void myb_internal_text_fixed(myb_internal_text* t, double value, int decimals, int width) {
  char number[NUMBER_CHARS];
  size_t length = 0;
  const myb_status status = myb_format_fixed(value, decimals, number, sizeof number, &length);
  MYB_CHECK(status == MYB_OK);
  aligned(t, number, length, width);
}

void myb_internal_text_scientific(myb_internal_text* t, double value, int digits, int width) {
  char number[NUMBER_CHARS];
  size_t length = 0;
  const myb_status status = myb_format_scientific(value, digits, number, sizeof number, &length);
  MYB_CHECK(status == MYB_OK);
  aligned(t, number, length, width);
}

myb_status myb_internal_text_finish(myb_internal_text* t, size_t* length) {
  if (length != NULL) {
    *length = t->length;
  }
  if (t->length + 1 > t->size) {
    if (t->size > 0) {
      t->out[t->size - 1] = '\0';
    }
    return MYB_E_OVERFLOW;
  }
  t->out[t->length] = '\0';
  return MYB_OK;
}
