/* SPDX-License-Identifier: BSD-3-Clause */
/* NAIF text kernel reader: syntax, replacement and appending, errors with line numbers,
 * atomic loading, and a UTF-8 (Cyrillic) file path. */
#include <string.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/io.h"

static myb_kpool* new_pool(void) {
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  return pool;
}

static myb_status load(myb_kpool* pool, const char* text, myb_diag* diag) {
  return myb_kpool_load_text(pool, text, strlen(text), diag);
}

static int string_is(const myb_kpool* pool, const char* name, size_t index, const char* expected) {
  const char* value = NULL;
  size_t length = 0;
  if (myb_kpool_string(pool, name, index, &value, &length) != MYB_OK) {
    return 0;
  }
  return length == strlen(expected) && memcmp(value, expected, length) == 0;
}

static void test_syntax(void) {
  static const char text[] =
      "KPL/PCK\n"
      "Commentary may contain anything = ( 'not data' ).\n"
      "\\begindata\n"
      "   SCALAR   = 42\n"
      "   LIST     = ( 1, 2.5D0  -3e2\n"
      "                0x10 )\n"
      "   NAMES    = ( 'TAI' 'it''s' '' )\n"
      "   DATES    = @1972-JAN-1\n"
      "   LIST    += ( 7 )\n"
      "\\begintext\n"
      "More commentary.\n"
      "  \\begindata  \n"
      "SCALAR = ( 43 )\n";
  myb_kpool* pool = new_pool();
  myb_diag diag;
  MT_CHECK(load(pool, text, &diag) == MYB_OK);
  MT_CHECK(diag.message[0] == '\0');

  const double* values = NULL;
  size_t count = 0;
  MT_CHECK(myb_kpool_numbers(pool, "SCALAR", &values, &count) == MYB_OK);
  MT_CHECK(count == 1 && values[0] == 43.0); /* the later '=' replaces */
  MT_CHECK(myb_kpool_numbers(pool, "LIST", &values, &count) == MYB_OK);
  MT_CHECK(count == 5);
  if (count == 5) {
    MT_CHECK(values[0] == 1.0 && values[1] == 2.5 && values[2] == -300.0 && values[3] == 16.0 &&
             values[4] == 7.0);
  }
  myb_kvar_type type = 0;
  MT_CHECK(myb_kpool_describe(pool, "NAMES", &type, &count) == MYB_OK);
  MT_CHECK(type == MYB_KVAR_STRINGS && count == 3);
  MT_CHECK(string_is(pool, "NAMES", 0, "TAI"));
  MT_CHECK(string_is(pool, "NAMES", 1, "it's"));
  MT_CHECK(string_is(pool, "NAMES", 2, ""));
  MT_CHECK(string_is(pool, "DATES", 0, "@1972-JAN-1"));

  const char* s = NULL;
  size_t length = 0;
  MT_CHECK(myb_kpool_string(pool, "NAMES", 3, &s, &length) == MYB_E_RANGE);
  MT_CHECK(myb_kpool_string(pool, "LIST", 0, &s, &length) == MYB_E_FORMAT);
  MT_CHECK(myb_kpool_numbers(pool, "NAMES", &values, &count) == MYB_E_FORMAT);
  MT_CHECK(myb_kpool_numbers(pool, "MISSING", &values, &count) == MYB_E_NOT_FOUND);
  MT_CHECK(myb_kpool_describe(pool, "MISSING", NULL, NULL) == MYB_E_NOT_FOUND);
  /* Text outside data blocks is never data. */
  MT_CHECK(myb_kpool_describe(pool, "Commentary", NULL, NULL) == MYB_E_NOT_FOUND);

  /* Appending across loads moves the old values in front of the new ones. */
  MT_CHECK(load(pool, "\\begindata\nOTHER = 1\nLIST += ( 8 9 )\n", NULL) == MYB_OK);
  MT_CHECK(myb_kpool_numbers(pool, "LIST", &values, &count) == MYB_OK);
  MT_CHECK(count == 7 && values[0] == 1.0 && values[4] == 7.0 && values[5] == 8.0 &&
           values[6] == 9.0);
  /* '+=' on a new variable defines it. */
  MT_CHECK(load(pool, "\\begindata\nFRESH += 'a'\n", NULL) == MYB_OK);
  MT_CHECK(string_is(pool, "FRESH", 0, "a"));
  myb_kpool_destroy(pool);
}

typedef struct bad_case {
  const char* text;
  long line;
  const char* fragment; /* expected part of the diagnostic */
} bad_case;

static void test_errors(void) {
  static const bad_case cases[] = {
      {"\\begindata\nA = ( 1 2\n", 2, "unterminated '(' list"}, /* line of the '(' */
      {"\\begindata\nA = 'open\n", 2, "unterminated string"},
      {"\\begindata\nA = ( 1 'x' )\n", 2, "mixed"},
      {"\\begindata\n\nA = ( 1.0.0 )\n", 3, "'1.0.0' is not a number"},
      {"\\begindata\nA 1\n", 2, "expected '=' or '+='"},
      {"\\begindata\nA = ( )\n", 2, "empty value list"},
      {"\\begindata\nA_NAME_THAT_IS_LONGER_THAN_32_CHARS = 1\n", 2, "longer than 32"},
      {"\\begindata\n= 1\n", 2, "expected a variable name"},
      {"\\begindata\nA = 1\nA += 'x'\n", 3, "mixed"},
      {"\\begindata\nA = ( 1\n\\begintext\n", 2, "unterminated '(' list"},
  };
  for (size_t i = 0; i < sizeof cases / sizeof cases[0]; ++i) {
    myb_kpool* pool = new_pool();
    MT_CHECK(load(pool, "\\begindata\nKEEP = 5\n", NULL) == MYB_OK);
    myb_diag diag;
    const myb_status status = load(pool, cases[i].text, &diag);
    mtest.checks++;
    if (status != MYB_E_FORMAT || diag.line != cases[i].line ||
        strstr(diag.message, cases[i].fragment) == NULL) {
      mtest.failures++;
      fprintf(stderr, "case %zu: status %s, line %ld, message '%s'\n", i, myb_status_name(status),
              diag.line, diag.message);
    }
    /* Atomic: the pool still holds exactly what it held before. */
    const double* values = NULL;
    size_t count = 0;
    MT_CHECK(myb_kpool_numbers(pool, "KEEP", &values, &count) == MYB_OK && count == 1 &&
             values[0] == 5.0);
    MT_CHECK(myb_kpool_describe(pool, "A", NULL, NULL) == MYB_E_NOT_FOUND);
    myb_kpool_destroy(pool);
  }
}

static void test_files(void) {
  myb_kpool* pool = new_pool();
  myb_diag diag;
  MT_CHECK(myb_kpool_load_file(pool,
                               MYB_TEST_DATA_DIR "/kernels/\xD1\x8F\xD0\xB4\xD1\x80\xD0\xBE.tk",
                               &diag) == MYB_OK);
  const double* values = NULL;
  size_t count = 0;
  MT_CHECK(myb_kpool_numbers(pool, "TEST_PATH_CHECK", &values, &count) == MYB_OK && count == 2 &&
           values[0] == 1.5 && values[1] == 2.5);
  MT_CHECK(string_is(pool, "TEST_LABEL", 0, "\xD1\x8F\xD0\xB4\xD1\x80\xD0\xBE"));
  MT_CHECK(myb_kpool_load_file(pool, MYB_TEST_DATA_DIR "/kernels/missing.tk", &diag) == MYB_E_IO);
  MT_CHECK(strstr(diag.message, "cannot read") != NULL);
  myb_kpool_destroy(pool);
  myb_kpool_destroy(NULL);
}

int main(void) {
  test_syntax();
  test_errors();
  test_files();
  return mt_finish("test_kpool");
}
