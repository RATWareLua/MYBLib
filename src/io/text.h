/* SPDX-License-Identifier: BSD-3-Clause */
/* io/text.h — building text in a caller's buffer for the writers of file formats: characters
 * beyond the buffer are counted, not stored, so one pass tells the size a file needs and the
 * next writes it. Numbers go through myb_format_* (exact, the same on every platform). */
#ifndef MYB_SRC_IO_TEXT_H_
#define MYB_SRC_IO_TEXT_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"

typedef struct myb_internal_text {
  char* out;     /* may be NULL when size is 0 */
  size_t size;   /* bytes of out, the terminating NUL included */
  size_t length; /* characters produced so far */
} myb_internal_text;

void myb_internal_text_char(myb_internal_text* t, char c);
void myb_internal_text_put(myb_internal_text* t, const char* s);
void myb_internal_text_span(myb_internal_text* t, const char* s, size_t length);

/* `count` copies of c. */
void myb_internal_text_fill(myb_internal_text* t, char c, size_t count);

/* A decimal integer right-aligned in at least `width` characters, padded with `pad`. */
void myb_internal_text_int(myb_internal_text* t, int64_t value, int width, char pad);

/* The shortest decimal that reads back (myb_format_double). */
void myb_internal_text_double(myb_internal_text* t, double value);

/* `decimals` digits after the point, right-aligned in at least `width` characters. */
void myb_internal_text_fixed(myb_internal_text* t, double value, int decimals, int width);

/* `digits` significant digits with an exponent, right-aligned in at least `width`. */
void myb_internal_text_scientific(myb_internal_text* t, double value, int digits, int width);

/* Terminates the text. *length (may be NULL): the characters of the whole text.
 * MYB_E_OVERFLOW if the buffer could not hold it (what fits is kept, terminated). */
[[nodiscard]] myb_status myb_internal_text_finish(myb_internal_text* t, size_t* length);

#endif /* MYB_SRC_IO_TEXT_H_ */
