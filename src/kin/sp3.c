/* SPDX-License-Identifier: BSD-3-Clause */
/* SP3-c / SP3-d precise orbit files (IGS, ILRS): satellite positions and velocities in an
 * Earth-fixed frame at tabulated epochs. Format: S. Hilla, "The Extended Standard Product 3
 * Orbit Format (SP3-d)", IGS, 2016; SP3-c (2010) is the same for the records read here. */
#include <math.h>
#include <stdckdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/io.h"
#include "myb/kin.h"
#include "platform/platform.h"

enum { MAX_SATELLITES = 999, ID_LENGTH = 4 };

struct myb_sp3 {
  myb_allocator allocator;
  size_t satellites;
  size_t epochs;
  char (*ids)[ID_LENGTH];
  myb_calendar* labels;
  double* position; /* [epoch][satellite][3], m */
  double* velocity; /* m/s */
  uint8_t* present; /* bit 0: position, bit 1: velocity */
  char time_system[4];
  char frame[6];
};

typedef struct cursor {
  const char* text;
  size_t size;
  size_t at;
  long line;
} cursor;

static int next_line(cursor* c, const char** s, size_t* n) {
  if (c->at >= c->size) {
    return 0;
  }
  const size_t start = c->at;
  size_t end = start;
  /* bound: the size of the text */
  while (end < c->size && c->text[end] != '\n') {
    ++end;
  }
  c->at = end + 1;
  c->line++;
  size_t length = end - start;
  if (length > 0 && c->text[start + length - 1] == '\r') {
    --length;
  }
  *s = c->text + start;
  *n = length;
  return 1;
}

/* The field [from, from + width) of a line, blanks trimmed. */
static void field(const char* s, size_t n, size_t from, size_t width, const char** f, size_t* len) {
  if (from >= n) {
    *f = s + n;
    *len = 0;
    return;
  }
  size_t end = from + width < n ? from + width : n;
  /* bound: the field */
  while (from < end && s[from] == ' ') {
    ++from;
  }
  /* bound: the field */
  while (end > from && s[end - 1] == ' ') {
    --end;
  }
  *f = s + from;
  *len = end - from;
}

static int number(const char* s, size_t n, size_t from, size_t width, double* value) {
  const char* f = NULL;
  size_t len = 0;
  field(s, n, from, width, &f, &len);
  size_t used = 0;
  return len > 0 && myb_parse_double(f, len, value, &used) == MYB_OK && used == len &&
         isfinite(*value);
}

static int integer(const char* s, size_t n, size_t from, size_t width, long lo, long hi,
                   long* value) {
  double v = 0.0;
  if (!number(s, n, from, width, &v) || v != floor(v) || v < (double)lo || v > (double)hi) {
    return 0;
  }
  *value = (long)v;
  return 1;
}

/* Epoch fields of the first line and of epoch lines ("*  yyyy mm dd hh mm ss.ssssssss"). */
static int epoch(const char* s, size_t n, myb_calendar* c) {
  long y = 0;
  long mo = 0;
  long d = 0;
  long h = 0;
  long mi = 0;
  double sec = 0.0;
  if (!integer(s, n, 3, 4, 1900, 2999, &y) || !integer(s, n, 8, 2, 1, 12, &mo) ||
      !integer(s, n, 11, 2, 1, 31, &d) || !integer(s, n, 14, 2, 0, 23, &h) ||
      !integer(s, n, 17, 2, 0, 59, &mi) || !number(s, n, 20, 11, &sec) || sec < 0.0 ||
      sec >= 61.0) {
    return 0;
  }
  *c = (myb_calendar){y, (int)mo, (int)d, (int)h, (int)mi, sec};
  return 1;
}

static myb_status fail(myb_diag* diag, long line, const char* what) {
  myb_internal_diag(diag, line, "SP3: %s", what);
  return MYB_E_FORMAT;
}

void myb_sp3_destroy(myb_sp3* sp3) {
  if (sp3 == NULL) {
    return;
  }
  const myb_allocator a = sp3->allocator;
  myb_internal_release(&a, sp3->ids);
  myb_internal_release(&a, sp3->labels);
  myb_internal_release(&a, sp3->position);
  myb_internal_release(&a, sp3->velocity);
  myb_internal_release(&a, sp3->present);
  myb_internal_release(&a, sp3);
}

myb_status myb_sp3_load(const char* text, size_t size, const myb_allocator* allocator,
                        myb_sp3** sp3, myb_diag* diag) {
  MYB_CHECK(text != NULL && sp3 != NULL);
  myb_internal_diag_clear(diag);
  cursor c = {text, size, 0, 0};
  const char* s = NULL;
  size_t n = 0;
  if (!next_line(&c, &s, &n) || n < 60 || s[0] != '#' || (s[1] != 'c' && s[1] != 'd')) {
    return fail(diag, 1, "not an SP3-c/d file (first line #c or #d)");
  }
  myb_calendar first;
  long declared = 0;
  if (!epoch(s, n, &first) || !integer(s, n, 32, 7, 1, 10000000, &declared)) {
    return fail(diag, 1, "invalid first line");
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_sp3), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_sp3* f = block;
  memset(f, 0, sizeof *f);
  f->allocator = a;
  {
    const char* fr = NULL;
    size_t len = 0;
    field(s, n, 46, 5, &fr, &len);
    memcpy(f->frame, fr, len < 5 ? len : 5);
  }
  /* header: satellites, time system */
  long count = -1;
  size_t listed = 0;
  int percent_c = 0;
  int header = 1;
  cursor records = c;
  /* bound: one line per iteration */
  while (header && status == MYB_OK) {
    records = c; /* the start of the line about to be read */
    if (!next_line(&c, &s, &n)) {
      break;
    }
    if (n >= 1 && s[0] == '*') {
      header = 0;
      break;
    }
    if (n >= 2 && s[0] == '+' && s[1] == ' ') {
      if (count < 0) {
        if (!integer(s, n, 3, 3, 1, MAX_SATELLITES, &count)) {
          status = fail(diag, c.line, "invalid number of satellites");
          break;
        }
        status = myb_internal_alloc_array(&a, (size_t)count, ID_LENGTH, &block);
        if (status != MYB_OK) {
          break;
        }
        f->ids = block;
      }
      for (size_t k = 0; k < 17 && listed < (size_t)count; ++k) {
        const size_t col = 9 + 3 * k;
        if (col + 3 > n) {
          break;
        }
        memcpy(f->ids[listed], s + col, 3);
        f->ids[listed][3] = '\0';
        for (int i = 0; i < 3; ++i) {
          if (f->ids[listed][i] == ' ') {
            f->ids[listed][i] = '0'; /* "G 1" of old files is G01 */
          }
        }
        ++listed;
      }
      continue;
    }
    if (n >= 12 && s[0] == '%' && s[1] == 'c' && percent_c++ == 0) {
      memcpy(f->time_system, s + 9, 3);
      f->time_system[3] = '\0';
    }
  }
  if (status == MYB_OK && (count < 1 || listed != (size_t)count)) {
    status = fail(diag, c.line, "satellite list incomplete");
  }
  if (status == MYB_OK && header) {
    status = fail(diag, c.line, "no epoch");
  }
  if (status != MYB_OK) {
    myb_sp3_destroy(f);
    return status;
  }
  f->satellites = (size_t)count;
  /* records: count the epochs, then fill */
  size_t epochs = 0;
  {
    cursor k = records;
    /* bound: one line per iteration */
    while (next_line(&k, &s, &n)) {
      if (n >= 1 && s[0] == '*') {
        ++epochs;
      }
    }
  }
  size_t cells = 0;
  if (ckd_mul(&cells, epochs, f->satellites)) {
    myb_sp3_destroy(f);
    return MYB_E_OVERFLOW;
  }
  status = myb_internal_alloc_array(&a, epochs, sizeof(myb_calendar), &block);
  if (status == MYB_OK) {
    f->labels = block;
    status = myb_internal_alloc_array(&a, 3 * cells, sizeof(double), &block);
  }
  if (status == MYB_OK) {
    f->position = block;
    status = myb_internal_alloc_array(&a, 3 * cells, sizeof(double), &block);
  }
  if (status == MYB_OK) {
    f->velocity = block;
    status = myb_internal_alloc_array(&a, cells, 1, &block);
  }
  if (status != MYB_OK) {
    myb_sp3_destroy(f);
    return status;
  }
  f->present = block;
  memset(f->present, 0, cells);
  f->epochs = epochs;
  c = records;
  size_t e = SIZE_MAX;
  /* bound: one line per iteration */
  while (status == MYB_OK && next_line(&c, &s, &n)) {
    if (n == 0 || (n >= 2 && s[0] == '/' && s[1] == '*')) {
      continue;
    }
    if (n >= 3 && memcmp(s, "EOF", 3) == 0) {
      break;
    }
    if (s[0] == '*') {
      ++e;
      if (!epoch(s, n, &f->labels[e])) {
        status = fail(diag, c.line, "invalid epoch line");
      }
      continue;
    }
    if (s[0] != 'P' && s[0] != 'V') {
      continue; /* EP/EV correlation records and others */
    }
    if (e == SIZE_MAX || n < 46) {
      status = fail(diag, c.line, "record before the first epoch or too short");
      break;
    }
    char id[ID_LENGTH];
    memcpy(id, s + 1, 3);
    id[3] = '\0';
    for (int i = 0; i < 3; ++i) {
      if (id[i] == ' ') {
        id[i] = '0';
      }
    }
    size_t sat = SIZE_MAX;
    for (size_t k = 0; k < f->satellites; ++k) {
      if (memcmp(f->ids[k], id, 3) == 0) {
        sat = k;
        break;
      }
    }
    if (sat == SIZE_MAX) {
      status = fail(diag, c.line, "satellite not in the header");
      break;
    }
    double v[3];
    if (!number(s, n, 4, 14, &v[0]) || !number(s, n, 18, 14, &v[1]) ||
        !number(s, n, 32, 14, &v[2])) {
      status = fail(diag, c.line, "invalid coordinates");
      break;
    }
    if (v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0) {
      continue; /* missing (bad or absent) value */
    }
    const size_t cell = e * f->satellites + sat;
    if (s[0] == 'P') {
      for (int i = 0; i < 3; ++i) {
        f->position[3 * cell + (size_t)i] = v[i] * 1e3; /* km */
      }
      f->present[cell] |= 1u;
    } else {
      for (int i = 0; i < 3; ++i) {
        f->velocity[3 * cell + (size_t)i] = v[i] * 1e-1; /* dm/s */
      }
      f->present[cell] |= 2u;
    }
  }
  if (status != MYB_OK) {
    myb_sp3_destroy(f);
    return status;
  }
  *sp3 = f;
  return MYB_OK;
}

myb_status myb_sp3_load_file(const char* path, const myb_allocator* allocator, myb_sp3** sp3,
                             myb_diag* diag) {
  MYB_CHECK(path != NULL && sp3 != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &a, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded = myb_sp3_load((const char*)map.data, map.size, allocator, sp3, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

size_t myb_sp3_satellite_count(const myb_sp3* sp3) {
  MYB_CHECK(sp3 != NULL);
  return sp3->satellites;
}

const char* myb_sp3_satellite(const myb_sp3* sp3, size_t index) {
  MYB_CHECK(sp3 != NULL && index < sp3->satellites);
  return sp3->ids[index];
}

size_t myb_sp3_epoch_count(const myb_sp3* sp3) {
  MYB_CHECK(sp3 != NULL);
  return sp3->epochs;
}

const char* myb_sp3_time_system(const myb_sp3* sp3) {
  MYB_CHECK(sp3 != NULL);
  return sp3->time_system;
}

const char* myb_sp3_frame(const myb_sp3* sp3) {
  MYB_CHECK(sp3 != NULL);
  return sp3->frame;
}

void myb_sp3_epoch(const myb_sp3* sp3, size_t epoch_index, myb_calendar* label) {
  MYB_CHECK(sp3 != NULL && label != NULL && epoch_index < sp3->epochs);
  *label = sp3->labels[epoch_index];
}

int myb_sp3_state(const myb_sp3* sp3, size_t satellite, size_t epoch_index, double position[3],
                  double velocity[3]) {
  MYB_CHECK(sp3 != NULL && satellite < sp3->satellites && epoch_index < sp3->epochs);
  const size_t cell = epoch_index * sp3->satellites + satellite;
  const int flags = sp3->present[cell];
  if (position != NULL) {
    memcpy(position, sp3->position + 3 * cell, 3 * sizeof(double));
  }
  if (velocity != NULL) {
    memcpy(velocity, sp3->velocity + 3 * cell, 3 * sizeof(double));
  }
  return flags;
}
