/* SPDX-License-Identifier: BSD-3-Clause */
/* Tables of numbers as text: rows of values separated by commas, semicolons, tabs or blanks
 * (CSV and its relatives). */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "io/text.h"
#include "myb/io.h"
#include "platform/platform.h"

enum { MAX_COLUMNS = 1024, NAME_CHARS = 32 };

static const size_t max_file = (size_t)1 << 31;

struct myb_table {
  myb_allocator allocator;
  size_t rows;
  size_t columns;
  double* values;
  char (*names)[NAME_CHARS];
};

void myb_table_destroy(myb_table* table) {
  if (table == NULL) {
    return;
  }
  const myb_allocator a = table->allocator;
  myb_internal_release(&a, table->values);
  myb_internal_release(&a, table->names);
  myb_internal_release(&a, table);
}

static int separator(char c) {
  return c == ',' || c == ';' || c == '\t';
}

static int space(char c) {
  return c == ' ' || c == '\r';
}

typedef struct fields {
  const char* s;
  size_t n;
  size_t at;
  int pending; /* a separator was just passed: a field follows even if it is empty */
} fields;

/* The next field of a line: its bounds, blanks trimmed; returns 0 at the end of the line.
 * Fields are separated by one comma, semicolon or tab, or by blanks. */
static int next_field(fields* f, size_t* begin, size_t* length) {
  size_t i = f->at;
  while (i < f->n && space(f->s[i])) {
    ++i;
  }
  if (i >= f->n && !f->pending) {
    return 0;
  }
  size_t end = i;
  while (end < f->n && !separator(f->s[end]) && !space(f->s[end])) {
    ++end;
  }
  *begin = i;
  *length = end - i;
  size_t next = end;
  while (next < f->n && space(f->s[next])) {
    ++next;
  }
  f->pending = next < f->n && separator(f->s[next]);
  f->at = f->pending ? next + 1 : next;
  return 1;
}

/* A field as a number: an empty one is a NaN. Returns 0 if it is not a number. */
static int field_value(const char* s, size_t length, double* value) {
  *value = (double)NAN;
  if (length == 0) {
    return 1;
  }
  size_t used = 0;
  return myb_parse_double(s, length, value, &used) == MYB_OK && used == length;
}

/* One pass: counts rows and columns, and with a table to fill stores names and values. */
static myb_status table_pass(const char* text, size_t size, myb_table* table, size_t* rows,
                             size_t* columns, myb_diag* diag) {
  size_t at = 0;
  long number = 0;
  size_t row = 0;
  size_t width = 0;
  int header_seen = 0;
  /* bound: the size of the text */
  while (at < size) {
    size_t end = at;
    while (end < size && text[end] != '\n') {
      ++end;
    }
    const char* s = text + at;
    size_t n = end - at;
    at = end + 1;
    ++number;
    for (size_t i = 0; i < n; ++i) {
      if (s[i] == '#') {
        n = i;
        break;
      }
    }
    /* the fields of the line: how many, and whether all are numbers */
    fields f = {.s = s, .n = n};
    size_t begin = 0;
    size_t length = 0;
    size_t count = 0;
    int numeric = 1;
    /* bound: the fields of the line */
    while (next_field(&f, &begin, &length)) {
      double v = 0.0;
      if (!field_value(s + begin, length, &v)) {
        numeric = 0;
      }
      ++count;
    }
    if (count == 0) {
      continue;
    }
    if (count > MAX_COLUMNS) {
      myb_internal_diag(diag, number, "table: more than %d columns", MAX_COLUMNS);
      return MYB_E_FORMAT;
    }
    if (width != 0 && count != width) {
      myb_internal_diag(diag, number, "table: %zu fields where %zu are expected", count, width);
      return MYB_E_FORMAT;
    }
    width = count;
    if (!numeric && (row > 0 || header_seen)) {
      myb_internal_diag(diag, number, "table: a field is not a number");
      return MYB_E_FORMAT;
    }
    if (table != NULL) {
      f = (fields){.s = s, .n = n};
      for (size_t k = 0; k < count && next_field(&f, &begin, &length); ++k) {
        if (numeric) {
          (void)field_value(s + begin, length, &table->values[row * width + k]);
        } else {
          const size_t keep = length < NAME_CHARS - 1 ? length : NAME_CHARS - 1;
          memcpy(table->names[k], s + begin, keep);
          table->names[k][keep] = '\0';
        }
      }
    }
    if (numeric) {
      ++row;
    } else {
      header_seen = 1;
    }
  }
  *rows = row;
  *columns = width;
  return MYB_OK;
}

myb_status myb_table_load(const char* text, size_t size, const myb_allocator* allocator,
                          myb_table** table, myb_diag* diag) {
  MYB_CHECK(text != NULL && table != NULL);
  myb_internal_diag_clear(diag);
  size_t rows = 0;
  size_t columns = 0;
  myb_status status = table_pass(text, size, NULL, &rows, &columns, diag);
  if (status != MYB_OK) {
    return status;
  }
  if (columns == 0) {
    myb_internal_diag(diag, 0, "table: no values");
    return MYB_E_FORMAT;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  status = myb_internal_alloc_array(&a, 1, sizeof(myb_table), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_table* t = block;
  memset(t, 0, sizeof *t);
  t->allocator = a;
  t->rows = rows;
  t->columns = columns;
  if (rows > SIZE_MAX / sizeof(double) / columns) {
    myb_table_destroy(t);
    return MYB_E_OVERFLOW;
  }
  status = myb_internal_alloc_array(&a, rows * columns, sizeof(double), &block);
  if (status == MYB_OK) {
    t->values = block;
    status = myb_internal_alloc_array(&a, columns, NAME_CHARS, &block);
  }
  if (status == MYB_OK) {
    t->names = block;
    memset(t->names, 0, columns * NAME_CHARS);
    status = table_pass(text, size, t, &rows, &columns, diag);
  }
  if (status != MYB_OK) {
    myb_table_destroy(t);
    return status;
  }
  *table = t;
  return MYB_OK;
}

myb_status myb_table_load_file(const char* path, const myb_allocator* allocator,
                               myb_table** table, myb_diag* diag) {
  MYB_CHECK(path != NULL && table != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  char* data = NULL;
  size_t size = 0;
  const myb_status status = myb_platform_read_file(path, &a, max_file, &data, &size);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot read '%s'", path);
    return status;
  }
  const myb_status loaded = myb_table_load(data, size, allocator, table, diag);
  myb_internal_release(&a, data);
  return loaded;
}

size_t myb_table_rows(const myb_table* table) {
  MYB_CHECK(table != NULL);
  return table->rows;
}

size_t myb_table_columns(const myb_table* table) {
  MYB_CHECK(table != NULL);
  return table->columns;
}

const char* myb_table_name(const myb_table* table, size_t column) {
  MYB_CHECK(table != NULL && column < table->columns);
  return table->names[column];
}

const double* myb_table_row(const myb_table* table, size_t row) {
  MYB_CHECK(table != NULL && row < table->rows);
  return table->values + row * table->columns;
}

myb_status myb_table_write(const char* const* names, size_t columns, const double* values,
                           size_t rows, char* buffer, size_t size, size_t* length) {
  MYB_CHECK((values != NULL || rows == 0) && (buffer != NULL || size == 0));
  if (columns == 0 || columns > MAX_COLUMNS) {
    return MYB_E_DOMAIN;
  }
  myb_internal_text t = {.out = buffer, .size = size};
  if (names != NULL) {
    for (size_t c = 0; c < columns; ++c) {
      if (names[c] == NULL) {
        return MYB_E_DOMAIN;
      }
      if (c > 0) {
        myb_internal_text_char(&t, ',');
      }
      myb_internal_text_put(&t, names[c]);
    }
    myb_internal_text_char(&t, '\n');
  }
  for (size_t r = 0; r < rows; ++r) {
    for (size_t c = 0; c < columns; ++c) {
      if (c > 0) {
        myb_internal_text_char(&t, ',');
      }
      const double v = values[r * columns + c];
      if (v == v) { /* a NaN is an empty field */
        myb_internal_text_double(&t, v);
      }
    }
    myb_internal_text_char(&t, '\n');
  }
  return myb_internal_text_finish(&t, length);
}
