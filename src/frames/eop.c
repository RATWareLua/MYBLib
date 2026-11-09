/* SPDX-License-Identifier: BSD-3-Clause */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/frames.h"
#include "num/poisson.h"
#include "platform/platform.h"

static const double arcsec = 4.848136811095359935899141e-6; /* rad */
static const double mjd_j2000 = 51544.5;
static const size_t max_file = (size_t)1 << 30;

enum { FIELDS = 7 }; /* mjd, xp, yp, dut1, lod, dx, dy */

struct myb_eop {
  myb_allocator allocator;
  size_t count;
  double* data; /* FIELDS arrays of count values */
  int tides;
  myb_internal_poisson_group tidal; /* XP, YP, UT1 */
};

static double* column(const myb_eop* eop, int field) {
  return eop->data + (size_t)field * eop->count;
}

myb_status myb_eop_create(const myb_allocator* allocator, myb_eop** eop) {
  MYB_CHECK(eop != NULL);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_eop), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_eop* e = block;
  memset(e, 0, sizeof *e);
  e->allocator = a;
  *eop = e;
  return MYB_OK;
}

void myb_eop_destroy(myb_eop* eop) {
  if (eop == NULL) {
    return;
  }
  const myb_allocator a = eop->allocator;
  myb_internal_release(&a, eop->data);
  myb_internal_poisson_free(&eop->tidal, &a);
  myb_internal_release(&a, eop);
}

/* ---- Parsing ----------------------------------------------------------------------------- */

typedef struct row {
  double v[FIELDS];
} row;

/* A growing array of rows. */
typedef struct rows {
  const myb_allocator* allocator;
  row* items;
  size_t count;
  size_t capacity;
} rows;

static myb_status push(rows* r, const row* item) {
  if (r->count == r->capacity) {
    const size_t grown = r->capacity == 0 ? 1024 : 2 * r->capacity;
    void* block = NULL;
    const myb_status status = myb_internal_alloc_array(r->allocator, grown, sizeof(row), &block);
    if (status != MYB_OK) {
      return status;
    }
    if (r->count > 0) {
      memcpy(block, r->items, r->count * sizeof(row));
    }
    myb_internal_release(r->allocator, r->items);
    r->items = block;
    r->capacity = grown;
  }
  r->items[r->count++] = *item;
  return MYB_OK;
}

/* Moves rows into the table (columns), replacing it. */
static myb_status install(myb_eop* eop, rows* r) {
  void* block = NULL;
  const myb_status status =
      myb_internal_alloc_array(&eop->allocator, r->count * FIELDS, sizeof(double), &block);
  if (status != MYB_OK) {
    return status;
  }
  double* data = block;
  for (size_t i = 0; i < r->count; ++i) {
    for (int f = 0; f < FIELDS; ++f) {
      data[(size_t)f * r->count + i] = r->items[i].v[f];
    }
  }
  myb_internal_release(&eop->allocator, eop->data);
  eop->data = data;
  eop->count = r->count;
  return MYB_OK;
}

/* Parses a number at text[0 .. length) after blanks; NaN for a blank field. Returns 0 on a
 * malformed field. */
static int field(const char* text, size_t length, double* value) {
  size_t i = 0;
  /* bound: at most length characters */
  while (i < length && text[i] == ' ') {
    ++i;
  }
  if (i == length) {
    *value = (double)NAN;
    return 1;
  }
  size_t used = 0;
  if (myb_parse_double(text + i, length - i, value, &used) != MYB_OK) {
    return 0;
  }
  i += used;
  /* bound: at most length characters */
  while (i < length && text[i] == ' ') {
    ++i;
  }
  return i == length && isfinite(*value);
}

/* Fixed columns [first, last] (1-based, inclusive) of a line; missing columns are blank. */
static int columns(const char* line, size_t length, size_t first, size_t last, double* value) {
  if (first > length) {
    *value = (double)NAN;
    return 1;
  }
  const size_t end = last < length ? last : length;
  return field(line + first - 1, end - first + 1, value);
}

static double mjd_of_date(int64_t year, int month, int day, int* ok) {
  const myb_calendar c = {year, month, day, 0, 0, 0.0};
  myb_time t;
  *ok = myb_calendar_to_time(&c, &t) == MYB_OK;
  return mjd_j2000 + (double)t.s / 86400.0;
}

typedef enum file_format { FORMAT_C04, FORMAT_FINALS } file_format;

/* One data line into a row; 0 = not a data row (comment), 1 = row, -1 = error, 2 = end of
 * usable data (finals: no more polar motion or UT1). */
static int parse_line(file_format format, const char* line, size_t length, row* out) {
  if (format == FORMAT_C04) {
    if (length == 0 || line[0] == '#') {
      return 0;
    }
    double v[13];
    size_t i = 0;
    for (int k = 0; k < 13; ++k) {
      /* bound: at most the line length */
      while (i < length && (line[i] == ' ' || line[i] == '\t' || line[i] == '\r')) {
        ++i;
      }
      size_t used = 0;
      if (i == length || myb_parse_double(line + i, length - i, &v[k], &used) != MYB_OK ||
          !isfinite(v[k])) {
        return -1;
      }
      i += used;
    }
    int ok = 0;
    const double mjd = mjd_of_date((int64_t)v[0], (int)v[1], (int)v[2], &ok);
    if (!ok || v[3] != 0.0 || mjd != v[4]) {
      return -1;
    }
    out->v[0] = v[4];
    out->v[1] = v[5] * arcsec;
    out->v[2] = v[6] * arcsec;
    out->v[3] = v[7];
    out->v[5] = v[8] * arcsec;
    out->v[6] = v[9] * arcsec;
    out->v[4] = v[12];
    return 1;
  }
  /* finals2000A: fixed columns of the Bulletin A values */
  if (length < 15) {
    return length == 0 ? 0 : -1;
  }
  double mjd;
  double xp;
  double yp;
  double dut1;
  double lod;
  double dx;
  double dy;
  if (!columns(line, length, 8, 15, &mjd) || !isfinite(mjd) ||
      !columns(line, length, 19, 27, &xp) || !columns(line, length, 38, 46, &yp) ||
      !columns(line, length, 59, 68, &dut1) || !columns(line, length, 80, 86, &lod) ||
      !columns(line, length, 98, 106, &dx) || !columns(line, length, 117, 125, &dy)) {
    return -1;
  }
  if (!isfinite(xp) || !isfinite(yp) || !isfinite(dut1)) {
    return 2;
  }
  out->v[0] = mjd;
  out->v[1] = xp * arcsec;
  out->v[2] = yp * arcsec;
  out->v[3] = dut1;
  out->v[4] = lod * 1e-3;          /* ms -> s; NaN stays NaN */
  out->v[5] = dx * 1e-3 * arcsec;  /* mas -> rad */
  out->v[6] = dy * 1e-3 * arcsec;
  return 1;
}

static myb_status load(myb_eop* eop, file_format format, const char* text, size_t size,
                       myb_diag* diag) {
  MYB_CHECK(eop != NULL && (text != NULL || size == 0));
  myb_internal_diag_clear(diag);
  rows r = {&eop->allocator, NULL, 0, 0};
  myb_status status = MYB_OK;
  size_t start = 0;
  long number = 0;
  /* bound: each iteration consumes one line, at most size + 1 lines */
  while (start < size && status == MYB_OK) {
    size_t end = start;
    /* bound: at most size characters */
    while (end < size && text[end] != '\n') {
      ++end;
    }
    ++number;
    size_t length = end - start;
    if (length > 0 && text[start + length - 1] == '\r') {
      --length;
    }
    row item;
    const int kind = parse_line(format, text + start, length, &item);
    if (kind == -1) {
      myb_internal_diag(diag, number, "malformed EOP line");
      status = MYB_E_FORMAT;
    } else if (kind == 2) {
      break;
    } else if (kind == 1) {
      if (r.count > 0 && !(item.v[0] > r.items[r.count - 1].v[0])) {
        myb_internal_diag(diag, number, "MJD not increasing");
        status = MYB_E_FORMAT;
      } else {
        status = push(&r, &item);
      }
    }
    start = end + 1;
  }
  if (status == MYB_OK && r.count < 4) {
    myb_internal_diag(diag, 0, "fewer than 4 EOP rows");
    status = MYB_E_FORMAT;
  }
  if (status == MYB_OK) {
    status = install(eop, &r);
  }
  myb_internal_release(&eop->allocator, r.items);
  return status;
}

static myb_status load_file(myb_eop* eop, file_format format, const char* path, myb_diag* diag) {
  MYB_CHECK(eop != NULL && path != NULL);
  myb_internal_diag_clear(diag);
  char* text = NULL;
  size_t size = 0;
  myb_status status = myb_platform_read_file(path, &eop->allocator, max_file, &text, &size);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot read '%s'", path);
    return status;
  }
  status = load(eop, format, text, size, diag);
  myb_internal_release(&eop->allocator, text);
  return status;
}

myb_status myb_eop_load_c04(myb_eop* eop, const char* text, size_t size, myb_diag* diag) {
  return load(eop, FORMAT_C04, text, size, diag);
}

myb_status myb_eop_load_c04_file(myb_eop* eop, const char* path, myb_diag* diag) {
  return load_file(eop, FORMAT_C04, path, diag);
}

myb_status myb_eop_load_finals(myb_eop* eop, const char* text, size_t size, myb_diag* diag) {
  return load(eop, FORMAT_FINALS, text, size, diag);
}

myb_status myb_eop_load_finals_file(myb_eop* eop, const char* path, myb_diag* diag) {
  return load_file(eop, FORMAT_FINALS, path, diag);
}

myb_status myb_eop_set_tides(myb_eop* eop, const myb_kpool* pool, const char* prefix,
                             myb_diag* diag) {
  MYB_CHECK(eop != NULL && pool != NULL && prefix != NULL);
  myb_internal_diag_clear(diag);
  char names[3][128];
  const char* list[3] = {names[0], names[1], names[2]};
  const char* suffixes[3] = {"XP", "YP", "UT1"};
  for (int i = 0; i < 3; ++i) {
    snprintf(names[i], sizeof names[i], "%s_%s", prefix, suffixes[i]);
  }
  myb_internal_poisson_group group;
  const myb_status status =
      myb_internal_poisson_load(pool, prefix, list, 3, &eop->allocator, &group, diag);
  if (status != MYB_OK) {
    return status;
  }
  myb_internal_poisson_free(&eop->tidal, &eop->allocator);
  eop->tidal = group;
  eop->tides = 1;
  return MYB_OK;
}

size_t myb_eop_rows(const myb_eop* eop, double* first, double* last) {
  MYB_CHECK(eop != NULL);
  if (first != NULL) {
    *first = eop->count > 0 ? column(eop, 0)[0] : 0.0;
  }
  if (last != NULL) {
    *last = eop->count > 0 ? column(eop, 0)[eop->count - 1] : 0.0;
  }
  return eop->count;
}

/* ---- Interpolation ----------------------------------------------------------------------- */

/* The reading of a UTC label as if UTC were continuous (a leap second 23:59:60.x reads as
 * 23:59:59 + 1.x). */
static myb_status label_reading(const myb_calendar* utc, myb_time* t) {
  if (utc->second < 60.0) {
    return myb_calendar_to_time(utc, t);
  }
  myb_calendar c = *utc;
  c.second = 59.0;
  const myb_status status = myb_calendar_to_time(&c, t);
  if (status == MYB_OK) {
    *t = myb_time_add(*t, utc->second - 59.0);
  }
  return status;
}

myb_status myb_eop_at(const myb_eop* eop, const myb_timesys* system, myb_time t, myb_scale scale,
                      myb_eop_values* values, myb_time* ut1) {
  MYB_CHECK(eop != NULL && system != NULL && values != NULL);
  myb_calendar utc;
  myb_status status = myb_timesys_to_utc(system, t, scale, &utc);
  if (status != MYB_OK) {
    return status;
  }
  myb_time label;
  status = label_reading(&utc, &label);
  if (status != MYB_OK) {
    return status;
  }
  const double mjd = mjd_j2000 + myb_time_seconds(label) / 86400.0;
  const double* m = column(eop, 0);
  /* j: the last row with m[j] <= mjd; nodes j-1 .. j+2 */
  size_t lo = 0;
  size_t hi = eop->count;
  /* bound: the interval halves */
  while (lo < hi) {
    const size_t mid = lo + (hi - lo) / 2;
    if (m[mid] <= mjd) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  if (lo < 2 || lo + 2 > eop->count) {
    return MYB_E_RANGE;
  }
  const size_t first = lo - 2;
  /* The nodes as instants of the caller's scale (abscissas in days from t, continuous across
   * leap seconds, unlike UTC labels), and UT1 minus that scale at each of them. */
  double x[4];
  double d[4];
  for (int i = 0; i < 4; ++i) {
    const size_t r = first + (size_t)i;
    const double dut1 = column(eop, 3)[r];
    if (!isfinite(dut1)) {
      return MYB_E_RANGE;
    }
    const myb_time node_label = myb_time_from_seconds((m[r] - mjd_j2000) * 86400.0);
    myb_calendar node;
    myb_time_to_calendar(node_label, &node);
    myb_time node_t;
    status = myb_timesys_from_utc(system, &node, scale, &node_t);
    if (status != MYB_OK) {
      return status;
    }
    x[i] = myb_time_diff(node_t, t).hi / 86400.0;
    d[i] = myb_time_diff(node_label, node_t).hi + dut1;
  }
  double w[4];
  for (int i = 0; i < 4; ++i) {
    double p = 1.0;
    for (int k = 0; k < 4; ++k) {
      if (k != i) {
        p *= x[k] / (x[k] - x[i]);
      }
    }
    w[i] = p;
  }
  double offset = 0.0;
  for (int i = 0; i < 4; ++i) {
    offset += w[i] * d[i];
  }
  double v[FIELDS] = {0};
  for (int f = 1; f < FIELDS; ++f) {
    if (f == 3) {
      continue;
    }
    double sum = 0.0;
    int missing = 0;
    for (int i = 0; i < 4; ++i) {
      const double value = column(eop, f)[first + (size_t)i];
      missing |= !isfinite(value);
      sum += w[i] * value;
    }
    if (missing && (f == 1 || f == 2)) {
      return MYB_E_RANGE;
    }
    v[f] = missing ? 0.0 : sum;
  }
  if (eop->tides) {
    double tide[3];
    double rate[3];
    /* argument time: UT1 (continuous; chi = GMST(UT1) + pi), where the IERS routines take
     * the UTC MJD, which differs by |UT1 - UTC| < 0.9 s and repeats during a leap second */
    const double ut1_days = myb_time_seconds(myb_time_add(t, offset)) / 86400.0;
    myb_internal_poisson_eval(&eop->tidal, ut1_days / 36525.0, tide, rate);
    v[1] += tide[0];
    v[2] += tide[1];
    offset += tide[2];
    v[4] -= rate[2] / 36525.0; /* LOD = -dUT1/d(day) */
  }
  const myb_time u = myb_time_add(t, offset);
  values->xp = v[1];
  values->yp = v[2];
  values->dut1 = myb_time_diff(u, label).hi;
  values->lod = v[4];
  values->dx = v[5];
  values->dy = v[6];
  if (ut1 != NULL) {
    *ut1 = u;
  }
  return MYB_OK;
}
