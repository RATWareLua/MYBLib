/* SPDX-License-Identifier: BSD-3-Clause */
/* SP3-c / SP3-d precise orbit files (IGS, ILRS): satellite positions, velocities and clocks
 * in an Earth-fixed frame at tabulated epochs, with their accuracy codes. Format: S. Hilla, "The Extended Standard Product 3
 * Orbit Format (SP3-d)", IGS, 2016; SP3-c (2010) is the same for the records read here. */
#include <math.h>
#include <stdckdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "io/text.h"
#include "myb/cmath.h"
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
  uint8_t* present; /* bit 0: position, bit 1: velocity, bit 2: clock, bit 3: clock rate */
  double* clock;    /* [epoch][satellite][2]: s, s/s */
  double* sigma;    /* [epoch][satellite][8], as myb_sp3_sigma */
  uint8_t* flags;   /* MYB_SP3_* */
  double* accuracy; /* [satellite], m */
  double base[2];   /* of the exponents: positions and velocities, clocks */
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

/* base^n units: what an accuracy code stands for; 0 for a blank field or a base not given. */
static double from_code(const char* s, size_t n, size_t from, size_t width, double base,
                        double unit) {
  long code = 0;
  if (!(base > 1.0) || !integer(s, n, from, width, 0, 999, &code)) {
    return 0.0;
  }
  return unit * myb_exp((double)code * myb_log(base));
}

/* The code nearest to a value: -1 if there is none to write. */
static int to_code(double value, double unit, double base, int lowest, int highest) {
  if (!(value > 0.0) || !isfinite(value)) {
    return -1;
  }
  const double n = floor(myb_log(value / unit) / myb_log(base) + 0.5);
  return n < (double)lowest ? lowest : n > (double)highest ? highest : (int)n;
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
  myb_internal_release(&a, sp3->clock);
  myb_internal_release(&a, sp3->sigma);
  myb_internal_release(&a, sp3->flags);
  myb_internal_release(&a, sp3->accuracy);
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
  int percent_f = 0;
  size_t rated = 0;
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
        status = myb_internal_alloc_array(&a, (size_t)count, sizeof(double), &block);
        if (status != MYB_OK) {
          break;
        }
        f->accuracy = block;
        memset(f->accuracy, 0, (size_t)count * sizeof(double));
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
    if (n >= 2 && s[0] == '+' && s[1] == '+' && count > 0) {
      /* 2^code mm for each satellite, in the order of the list; 0: unknown */
      for (size_t k = 0; k < 17 && rated < (size_t)count; ++k, ++rated) {
        long code = 0;
        if (integer(s, n, 9 + 3 * k, 3, 1, 999, &code)) {
          f->accuracy[rated] = 1e-3 * myb_exp((double)code * myb_log(2.0));
        }
      }
      continue;
    }
    if (n >= 12 && s[0] == '%' && s[1] == 'c' && percent_c++ == 0) {
      memcpy(f->time_system, s + 9, 3);
      f->time_system[3] = '\0';
    }
    if (n >= 26 && s[0] == '%' && s[1] == 'f' && percent_f++ == 0) {
      double base = 0.0;
      if (number(s, n, 3, 10, &base)) {
        f->base[0] = base;
      }
      if (number(s, n, 14, 12, &base)) {
        f->base[1] = base;
      }
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
    status = myb_internal_alloc_array(&a, 2 * cells, sizeof(double), &block);
  }
  if (status == MYB_OK) {
    f->clock = block;
    status = myb_internal_alloc_array(&a, 8 * cells, sizeof(double), &block);
  }
  if (status == MYB_OK) {
    f->sigma = block;
    status = myb_internal_alloc_array(&a, cells, 1, &block);
  }
  if (status == MYB_OK) {
    f->flags = block;
    status = myb_internal_alloc_array(&a, cells, 1, &block);
  }
  if (status != MYB_OK) {
    myb_sp3_destroy(f);
    return status;
  }
  f->present = block;
  memset(f->present, 0, cells);
  memset(f->flags, 0, cells);
  /* bound: the cells */
  for (size_t k = 0; k < 3 * cells; ++k) {
    f->position[k] = 0.0;
    f->velocity[k] = 0.0;
  }
  for (size_t k = 0; k < 2 * cells; ++k) {
    f->clock[k] = 0.0;
  }
  for (size_t k = 0; k < 8 * cells; ++k) {
    f->sigma[k] = 0.0;
  }
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
    const size_t cell = e * f->satellites + sat;
    const int is_position = s[0] == 'P';
    /* the clock (microseconds) or its rate (1e-4 microseconds a second): 999999.999999
     * stands for none */
    double clock = 0.0;
    if (number(s, n, 46, 14, &clock) && fabs(clock) < 999999.0) {
      f->clock[2 * cell + (is_position ? 0u : 1u)] = clock * (is_position ? 1e-6 : 1e-10);
      f->present[cell] |= is_position ? 4u : 8u;
    }
    double* sigma = f->sigma + 8 * cell + (is_position ? 0u : 4u);
    const double unit = is_position ? 1e-3 : 1e-7;         /* mm, 1e-4 mm/s */
    const double clock_unit = is_position ? 1e-12 : 1e-16; /* ps, 1e-4 ps/s */
    for (size_t i = 0; i < 3; ++i) {
      sigma[i] = from_code(s, n, 61 + 3 * i, 2, f->base[0], unit);
    }
    sigma[3] = from_code(s, n, 70, 3, f->base[1], clock_unit);
    if (is_position) {
      unsigned flags = 0;
      flags |= n > 74 && s[74] == 'E' ? (unsigned)MYB_SP3_CLOCK_EVENT : 0u;
      flags |= n > 75 && s[75] == 'P' ? (unsigned)MYB_SP3_CLOCK_PREDICTED : 0u;
      flags |= n > 78 && s[78] == 'M' ? (unsigned)MYB_SP3_MANEUVER : 0u;
      flags |= n > 79 && s[79] == 'P' ? (unsigned)MYB_SP3_ORBIT_PREDICTED : 0u;
      f->flags[cell] = (uint8_t)flags;
    }
    if (v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0) {
      continue; /* missing (bad or absent) value */
    }
    if (is_position) {
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
  const int flags = sp3->present[cell] & 3;
  if (position != NULL) {
    memcpy(position, sp3->position + 3 * cell, 3 * sizeof(double));
  }
  if (velocity != NULL) {
    memcpy(velocity, sp3->velocity + 3 * cell, 3 * sizeof(double));
  }
  return flags;
}

int myb_sp3_clock(const myb_sp3* sp3, size_t satellite, size_t epoch_index, double* clock,
                  double* rate) {
  MYB_CHECK(sp3 != NULL && satellite < sp3->satellites && epoch_index < sp3->epochs);
  const size_t cell = epoch_index * sp3->satellites + satellite;
  if (clock != NULL) {
    *clock = sp3->clock[2 * cell];
  }
  if (rate != NULL) {
    *rate = sp3->clock[2 * cell + 1];
  }
  return (sp3->present[cell] >> 2) & 3;
}

void myb_sp3_sigma(const myb_sp3* sp3, size_t satellite, size_t epoch_index, double sigma[8]) {
  MYB_CHECK(sp3 != NULL && sigma != NULL && satellite < sp3->satellites &&
            epoch_index < sp3->epochs);
  memcpy(sigma, sp3->sigma + 8 * (epoch_index * sp3->satellites + satellite),
         8 * sizeof(double));
}

double myb_sp3_accuracy(const myb_sp3* sp3, size_t satellite) {
  MYB_CHECK(sp3 != NULL && satellite < sp3->satellites);
  return sp3->accuracy[satellite];
}

unsigned myb_sp3_flags(const myb_sp3* sp3, size_t satellite, size_t epoch_index) {
  MYB_CHECK(sp3 != NULL && satellite < sp3->satellites && epoch_index < sp3->epochs);
  return sp3->flags[epoch_index * sp3->satellites + satellite];
}

/* ---- Writing ---------------------------------------------------------------------------------- */

/* The columns after the coordinates: the clock, the four accuracy codes and, for a position,
 * the flags. Without codes and flags the line ends at the clock, as files without them do. */
static void record_tail(myb_internal_text* t, int have_clock, double clock, const double* sigma,
                        double unit, double clock_unit, int with_flags, unsigned flags) {
  if (have_clock) {
    myb_internal_text_fixed(t, clock, 6, 14);
  } else {
    myb_internal_text_put(t, " 999999.999999");
  }
  int codes[4] = {-1, -1, -1, -1};
  if (sigma != NULL) {
    for (int i = 0; i < 3; ++i) {
      codes[i] = to_code(sigma[i], unit, 1.25, 0, 99);
    }
    codes[3] = to_code(sigma[3], clock_unit, 1.025, 0, 999);
  }
  if (codes[0] < 0 && codes[1] < 0 && codes[2] < 0 && codes[3] < 0 && flags == 0) {
    myb_internal_text_char(t, '\n');
    return;
  }
  for (int i = 0; i < 4; ++i) {
    myb_internal_text_char(t, ' ');
    if (codes[i] < 0) {
      myb_internal_text_put(t, i < 3 ? "  " : "   ");
    } else {
      myb_internal_text_int(t, codes[i], i < 3 ? 2 : 3, ' ');
    }
  }
  if (with_flags && flags != 0) {
    char tail[8] = {' ',
                    flags & (unsigned)MYB_SP3_CLOCK_EVENT ? 'E' : ' ',
                    flags & (unsigned)MYB_SP3_CLOCK_PREDICTED ? 'P' : ' ',
                    ' ',
                    ' ',
                    flags & (unsigned)MYB_SP3_MANEUVER ? 'M' : ' ',
                    flags & (unsigned)MYB_SP3_ORBIT_PREDICTED ? 'P' : ' ',
                    '\0'};
    size_t length = 7;
    /* bound: the seven columns */
    while (tail[length - 1] == ' ') {
      --length;
    }
    tail[length] = '\0';
    myb_internal_text_put(t, tail);
  }
  myb_internal_text_char(t, '\n');
}

static int sp3_text(const char* s, size_t max) {
  if (s == NULL) {
    return 0;
  }
  const size_t n = strlen(s);
  if (n == 0 || n > max) {
    return 0;
  }
  for (size_t i = 0; i < n; ++i) {
    if (s[i] < 32 || s[i] > 126) {
      return 0;
    }
  }
  return 1;
}

/* s left-justified in `width` characters. */
static void padded(myb_internal_text* t, const char* s, size_t width) {
  const size_t n = strlen(s);
  myb_internal_text_span(t, s, n < width ? n : width);
  if (n < width) {
    myb_internal_text_fill(t, ' ', width - n);
  }
}

static void sp3_epoch(myb_internal_text* t, const myb_calendar* c) {
  myb_internal_text_int(t, c->year, 4, ' ');
  myb_internal_text_char(t, ' ');
  myb_internal_text_int(t, c->month, 2, ' ');
  myb_internal_text_char(t, ' ');
  myb_internal_text_int(t, c->day, 2, ' ');
  myb_internal_text_char(t, ' ');
  myb_internal_text_int(t, c->hour, 2, ' ');
  myb_internal_text_char(t, ' ');
  myb_internal_text_int(t, c->minute, 2, ' ');
  myb_internal_text_char(t, ' ');
  myb_internal_text_fixed(t, c->second, 8, 11);
}

myb_status myb_sp3_write(const myb_sp3_data* data, char* buffer, size_t size, size_t* length) {
  MYB_CHECK(data != NULL && (buffer != NULL || size == 0));
  const myb_sp3_data* d = data;
  if ((d->version != 'c' && d->version != 'd') || !sp3_text(d->time_system, 3) ||
      !sp3_text(d->frame, 5) || !sp3_text(d->agency, 4) || !sp3_text(d->orbit_type, 3) ||
      !sp3_text(d->data_used, 5) || d->file_type < 'A' || d->file_type > 'Z' ||
      d->satellite_count == 0 || d->satellite_count > MAX_SATELLITES ||
      (d->version == 'c' && d->satellite_count > 85) || d->epoch_count == 0 ||
      d->epoch_count > 9999999 || d->satellites == NULL || d->epochs == NULL ||
      d->positions == NULL || !(d->interval > 0.0) || !isfinite(d->interval) ||
      (d->clock_rates != NULL && d->velocities == NULL)) {
    return MYB_E_DOMAIN;
  }
  /* a clock of a second or more does not fit the field */
  for (size_t k = 0; d->clocks != NULL && k < d->epoch_count * d->satellite_count; ++k) {
    if (isfinite(d->clocks[k]) && !(fabs(d->clocks[k]) < 0.999999)) {
      return MYB_E_DOMAIN;
    }
  }
  for (size_t k = 0; d->clock_rates != NULL && k < d->epoch_count * d->satellite_count; ++k) {
    if (isfinite(d->clock_rates[k]) && !(fabs(d->clock_rates[k]) < 0.999999e-4)) {
      return MYB_E_DOMAIN;
    }
  }
  for (size_t s = 0; s < d->satellite_count; ++s) {
    if (!sp3_text(d->satellites[s], 3) || strlen(d->satellites[s]) != 3) {
      return MYB_E_DOMAIN;
    }
  }
  /* the first epoch as week and seconds of the GPS count and as a Julian date: labels of the
   * header, computed from the calendar date whatever the time system */
  myb_time first;
  myb_status status = myb_calendar_to_time(&d->epochs[0], &first);
  if (status != MYB_OK) {
    return status;
  }
  const int64_t shifted = first.s + 43200; /* seconds since 2000-01-01 0h */
  const int64_t days = shifted >= 0 ? shifted / 86400 : -((-shifted + 86399) / 86400);
  const int64_t mjd = 51544 + days;
  const double seconds_of_day = (double)(shifted - days * 86400) + first.f;
  const int64_t since_gps = mjd - 44244;
  const int64_t week = since_gps >= 0 ? since_gps / 7 : -((-since_gps + 6) / 7);
  const double seconds_of_week = (double)((since_gps - 7 * week) * 86400) + seconds_of_day;
  myb_internal_text t = {.out = buffer, .size = size};
  /* line 1 */
  myb_internal_text_char(&t, '#');
  myb_internal_text_char(&t, d->version);
  myb_internal_text_char(&t, d->velocities != NULL ? 'V' : 'P');
  sp3_epoch(&t, &d->epochs[0]);
  myb_internal_text_char(&t, ' ');
  myb_internal_text_int(&t, (int64_t)d->epoch_count, 7, ' ');
  myb_internal_text_char(&t, ' ');
  padded(&t, d->data_used, 5);
  myb_internal_text_char(&t, ' ');
  padded(&t, d->frame, 5);
  myb_internal_text_char(&t, ' ');
  padded(&t, d->orbit_type, 3);
  myb_internal_text_char(&t, ' ');
  padded(&t, d->agency, 4);
  myb_internal_text_char(&t, '\n');
  /* line 2 */
  myb_internal_text_put(&t, "## ");
  myb_internal_text_int(&t, week, 4, ' ');
  myb_internal_text_char(&t, ' ');
  myb_internal_text_fixed(&t, seconds_of_week, 8, 15);
  myb_internal_text_char(&t, ' ');
  myb_internal_text_fixed(&t, d->interval, 8, 14);
  myb_internal_text_char(&t, ' ');
  myb_internal_text_int(&t, mjd, 5, ' ');
  myb_internal_text_char(&t, ' ');
  myb_internal_text_fixed(&t, seconds_of_day / 86400.0, 13, 15);
  myb_internal_text_char(&t, '\n');
  /* satellites, 17 a line, at least five lines; then as many lines of accuracy codes */
  const size_t lines = d->satellite_count <= 85 ? 5 : (d->satellite_count + 16) / 17;
  for (size_t l = 0; l < lines; ++l) {
    myb_internal_text_put(&t, "+  ");
    if (l == 0) {
      myb_internal_text_int(&t, (int64_t)d->satellite_count, 3, ' ');
    } else {
      myb_internal_text_put(&t, "   ");
    }
    myb_internal_text_put(&t, "   ");
    for (size_t k = 0; k < 17; ++k) {
      const size_t s = 17 * l + k;
      myb_internal_text_put(&t, s < d->satellite_count ? d->satellites[s] : "  0");
    }
    myb_internal_text_char(&t, '\n');
  }
  for (size_t l = 0; l < lines; ++l) {
    myb_internal_text_put(&t, "++       ");
    for (size_t k = 0; k < 17; ++k) {
      const size_t s = 17 * l + k;
      /* 2^code mm; 0 stands for unknown, so the least known accuracy is 2 mm */
      const int code = d->accuracy != NULL && s < d->satellite_count
                           ? to_code(d->accuracy[s], 1e-3, 2.0, 1, 999)
                           : -1;
      myb_internal_text_int(&t, code < 0 ? 0 : code, 3, ' ');
    }
    myb_internal_text_char(&t, '\n');
  }
  myb_internal_text_put(&t, "%c ");
  myb_internal_text_char(&t, d->file_type);
  myb_internal_text_put(&t, "  cc ");
  padded(&t, d->time_system, 3);
  myb_internal_text_put(&t, " ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n");
  myb_internal_text_put(&t, "%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n");
  /* the bases of the accuracy codes of the records */
  myb_internal_text_put(&t, d->sigmas != NULL
                                ? "%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n"
                                : "%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n");
  myb_internal_text_put(&t, "%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n");
  for (int k = 0; k < 2; ++k) {
    myb_internal_text_put(&t, "%i    0    0    0    0      0      0      0      0         0\n");
  }
  myb_internal_text_put(&t, "/* written by myb\n/*\n/*\n/*\n");
  /* epochs: positions in km, velocities in dm/s, clocks in microseconds, six decimals */
  for (size_t e = 0; e < d->epoch_count; ++e) {
    myb_internal_text_put(&t, "*  ");
    sp3_epoch(&t, &d->epochs[e]);
    myb_internal_text_char(&t, '\n');
    for (size_t s = 0; s < d->satellite_count; ++s) {
      const double* r = d->positions + 3 * (e * d->satellite_count + s);
      const int have = isfinite(r[0]) && isfinite(r[1]) && isfinite(r[2]);
      myb_internal_text_char(&t, 'P');
      myb_internal_text_put(&t, d->satellites[s]);
      for (int i = 0; i < 3; ++i) {
        myb_internal_text_fixed(&t, have ? r[i] / 1000.0 : 0.0, 6, 14);
      }
      const size_t cell = e * d->satellite_count + s;
      const int timed = d->clocks != NULL && isfinite(d->clocks[cell]);
      record_tail(&t, timed, timed ? d->clocks[cell] * 1e6 : 0.0,
                  d->sigmas != NULL ? d->sigmas + 8 * cell : NULL, 1e-3, 1e-12, 1,
                  d->flags != NULL ? d->flags[cell] : 0u);
      if (d->velocities != NULL) {
        const double* v = d->velocities + 3 * (e * d->satellite_count + s);
        const int moving = isfinite(v[0]) && isfinite(v[1]) && isfinite(v[2]);
        myb_internal_text_char(&t, 'V');
        myb_internal_text_put(&t, d->satellites[s]);
        for (int i = 0; i < 3; ++i) {
          myb_internal_text_fixed(&t, moving ? v[i] * 10.0 : 0.0, 6, 14);
        }
        const int drifting = d->clock_rates != NULL && isfinite(d->clock_rates[cell]);
        record_tail(&t, drifting, drifting ? d->clock_rates[cell] * 1e10 : 0.0,
                    d->sigmas != NULL ? d->sigmas + 8 * cell + 4 : NULL, 1e-7, 1e-16, 0, 0u);
      }
    }
  }
  myb_internal_text_put(&t, "EOF\n");
  return myb_internal_text_finish(&t, length);
}
