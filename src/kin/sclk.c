/* SPDX-License-Identifier: BSD-3-Clause */
/* Spacecraft clocks of type 1 (myb/ck.h), after "SCLK Required Reading" (NAIF, JPL): the
 * conventions — partition bounds rounded to whole ticks, a count that ends a partition
 * belonging to the next one, the row of the table chosen by the time on the way back — were
 * checked against CSPICE N0067 (tests/data/ck). */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/ck.h"

enum { MAX_FIELDS = MYB_SCLK_MAX_FIELDS, MAX_DIGITS = 18 };

struct myb_sclk {
  myb_allocator allocator;
  int system; /* 1 TDB, 2 TT */
  int fields;
  double moduli[MAX_FIELDS];
  double offsets[MAX_FIELDS];
  double ticks_per[MAX_FIELDS]; /* ticks in one count of the field */
  int width[MAX_FIELDS];        /* digits of the field's largest value */
  char delimiter;
  size_t partitions;
  double* part_start; /* three arrays of `partitions`: the bounds, and the encoded ticks */
  double* part_end;   /* at the end of each partition */
  double* part_total;
  size_t rows;
  double* table; /* three a row: encoded ticks, seconds past J2000, seconds per count */
};

static double nearest(double x) {
  return floor(x + 0.5);
}

/* The numbers of the variable <prefix>_<id>; MYB_E_NOT_FOUND with the name in diag. */
static myb_status numbers(const myb_kpool* pool, const char* prefix, int64_t id,
                          const double** values, size_t* count, int required, myb_diag* diag) {
  char name[64];
  snprintf(name, sizeof name, "%s_%" PRId64, prefix, id);
  const myb_status status = myb_kpool_numbers(pool, name, values, count);
  if (status == MYB_E_NOT_FOUND && !required) {
    *values = NULL;
    *count = 0;
    return MYB_OK;
  }
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "%s: %s", name,
                      status == MYB_E_NOT_FOUND ? "not defined" : "not numeric");
    return status;
  }
  for (size_t i = 0; i < *count; ++i) {
    if (!isfinite((*values)[i])) {
      myb_internal_diag(diag, 0, "%s: non-finite value", name);
      return MYB_E_FORMAT;
    }
  }
  return MYB_OK;
}

static myb_status invalid(myb_diag* diag, int64_t id, const char* what) {
  myb_internal_diag(diag, 0, "clock %" PRId64 ": %s", id, what);
  return MYB_E_FORMAT;
}

myb_status myb_sclk_create(const myb_kpool* pool, int64_t id, const myb_allocator* allocator,
                           myb_sclk** sclk, myb_diag* diag) {
  MYB_CHECK(pool != NULL && sclk != NULL);
  myb_internal_diag_clear(diag);
  const double* type = NULL;
  const double* system = NULL;
  const double* fields = NULL;
  const double* moduli = NULL;
  const double* offsets = NULL;
  const double* delimiter = NULL;
  const double* starts = NULL;
  const double* ends = NULL;
  const double* table = NULL;
  size_t n = 0;
  size_t n_system = 0;
  size_t n_moduli = 0;
  size_t n_offsets = 0;
  size_t n_starts = 0;
  size_t n_ends = 0;
  size_t n_table = 0;
  myb_status status = numbers(pool, "SCLK_DATA_TYPE", id, &type, &n, 1, diag);
  if (status != MYB_OK) {
    return status;
  }
  if (n != 1 || type[0] != 1.0) {
    return invalid(diag, id, "a clock of another type than 1");
  }
  status = numbers(pool, "SCLK01_TIME_SYSTEM", id, &system, &n_system, 0, diag);
  if (status == MYB_OK) {
    status = numbers(pool, "SCLK01_N_FIELDS", id, &fields, &n, 1, diag);
  }
  if (status == MYB_OK && n != 1) {
    status = invalid(diag, id, "SCLK01_N_FIELDS is not one number");
  }
  if (status == MYB_OK) {
    status = numbers(pool, "SCLK01_MODULI", id, &moduli, &n_moduli, 1, diag);
  }
  if (status == MYB_OK) {
    status = numbers(pool, "SCLK01_OFFSETS", id, &offsets, &n_offsets, 1, diag);
  }
  if (status == MYB_OK) {
    status = numbers(pool, "SCLK01_OUTPUT_DELIM", id, &delimiter, &n, 1, diag);
  }
  if (status == MYB_OK && (n != 1 || !(delimiter[0] >= 1.0 && delimiter[0] <= 5.0) ||
                           delimiter[0] != floor(delimiter[0]))) {
    status = invalid(diag, id, "SCLK01_OUTPUT_DELIM is not one of 1 .. 5");
  }
  if (status == MYB_OK) {
    status = numbers(pool, "SCLK_PARTITION_START", id, &starts, &n_starts, 1, diag);
  }
  if (status == MYB_OK) {
    status = numbers(pool, "SCLK_PARTITION_END", id, &ends, &n_ends, 1, diag);
  }
  if (status == MYB_OK) {
    status = numbers(pool, "SCLK01_COEFFICIENTS", id, &table, &n_table, 1, diag);
  }
  if (status != MYB_OK) {
    return status;
  }
  const double count = fields[0];
  if (!(count >= 1.0 && count <= (double)MAX_FIELDS) || count != floor(count) ||
      n_moduli != (size_t)count || n_offsets != (size_t)count) {
    return invalid(diag, id, "the fields, their moduli and offsets do not agree");
  }
  if (n_system > 1 || (n_system == 1 && system[0] != 1.0 && system[0] != 2.0)) {
    return invalid(diag, id, "SCLK01_TIME_SYSTEM is neither 1 (TDB) nor 2 (TT)");
  }
  if (n_starts == 0 || n_starts != n_ends) {
    return invalid(diag, id, "the partitions do not have a start and an end each");
  }
  if (n_table == 0 || n_table % 3 != 0) {
    return invalid(diag, id, "SCLK01_COEFFICIENTS is not a list of triples");
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  status = myb_internal_alloc_array(&a, 1, sizeof(myb_sclk), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_sclk* s = block;
  memset(s, 0, sizeof *s);
  s->allocator = a;
  s->system = n_system == 1 ? (int)system[0] : 1;
  s->fields = (int)count;
  s->delimiter = ".:-, "[(int)delimiter[0] - 1];
  s->partitions = n_starts;
  s->rows = n_table / 3;
  status = myb_internal_alloc_array(&a, 3 * n_starts + n_table, sizeof(double), &block);
  if (status != MYB_OK) {
    myb_internal_release(&a, s);
    return status;
  }
  s->part_start = block;
  s->part_end = s->part_start + n_starts;
  s->part_total = s->part_end + n_starts;
  s->table = s->part_total + n_starts;
  memcpy(s->table, table, n_table * sizeof *table);
  const char* fault = NULL;
  double per = 1.0;
  for (int i = s->fields; i-- > 0;) {
    s->moduli[i] = moduli[i];
    s->offsets[i] = offsets[i];
    s->ticks_per[i] = per;
    if (!(moduli[i] >= 1.0) || moduli[i] != floor(moduli[i]) || offsets[i] != floor(offsets[i])) {
      fault = "a modulus or an offset that is not a whole number";
      break;
    }
    /* the digits of the largest value of the field */
    const double largest = moduli[i] - 1.0 + offsets[i];
    int width = 1;
    double power = 10.0;
    /* bound: MAX_DIGITS digits */
    for (; width < MAX_DIGITS && largest >= power; ++width) {
      power *= 10.0;
    }
    s->width[i] = width;
    per *= moduli[i];
    if (!(per < 0x1p62)) {
      fault = "moduli whose product does not fit the ticks";
      break;
    }
  }
  double total = 0.0;
  for (size_t p = 0; p < n_starts && fault == NULL; ++p) {
    s->part_start[p] = nearest(starts[p]);
    s->part_end[p] = nearest(ends[p]);
    if (!(s->part_end[p] >= s->part_start[p]) || !(s->part_start[p] >= 0.0)) {
      fault = "a partition that ends before it starts";
    }
    total += nearest(s->part_end[p] - s->part_start[p]);
    s->part_total[p] = total;
  }
  for (size_t r = 0; r < s->rows && fault == NULL; ++r) {
    if (r > 0 && !(s->table[3 * r] > s->table[3 * (r - 1)])) {
      fault = "a table whose ticks do not increase";
    }
    if (!(s->table[3 * r + 2] > 0.0)) {
      fault = "a rate that is not positive";
    }
  }
  if (fault != NULL) {
    myb_sclk_destroy(s);
    return invalid(diag, id, fault);
  }
  *sclk = s;
  return MYB_OK;
}

void myb_sclk_destroy(myb_sclk* sclk) {
  if (sclk == NULL) {
    return;
  }
  const myb_allocator a = sclk->allocator;
  myb_internal_release(&a, sclk->part_start);
  myb_internal_release(&a, sclk);
}

int myb_sclk_time_system(const myb_sclk* sclk) {
  MYB_CHECK(sclk != NULL);
  return sclk->system;
}

size_t myb_sclk_partitions(const myb_sclk* sclk) {
  MYB_CHECK(sclk != NULL);
  return sclk->partitions;
}

/* Ticks in a count of the most significant field. */
static double per_count(const myb_sclk* s) {
  return s->ticks_per[0];
}

myb_status myb_sclk_to_time(const myb_sclk* sclk, double ticks, myb_time* t) {
  MYB_CHECK(sclk != NULL && t != NULL);
  const myb_sclk* s = sclk;
  if (!(ticks >= s->table[0]) || !(ticks <= s->part_total[s->partitions - 1])) {
    return MYB_E_RANGE;
  }
  size_t lo = 0;
  size_t hi = s->rows;
  /* bound: the interval halves; the last row at or before the ticks */
  while (hi - lo > 1) {
    const size_t mid = lo + (hi - lo) / 2;
    if (s->table[3 * mid] <= ticks) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  const double* row = s->table + 3 * lo;
  *t = myb_time_add(myb_time_from_seconds(row[1]), row[2] * (ticks - row[0]) / per_count(s));
  return MYB_OK;
}

myb_status myb_sclk_from_time(const myb_sclk* sclk, myb_time t, double* ticks) {
  MYB_CHECK(sclk != NULL && ticks != NULL);
  const myb_sclk* s = sclk;
  if (myb_time_compare(t, myb_time_from_seconds(s->table[1])) < 0) {
    return MYB_E_RANGE;
  }
  size_t lo = 0;
  size_t hi = s->rows;
  /* bound: the interval halves; the last row whose time is at or before t */
  while (hi - lo > 1) {
    const size_t mid = lo + (hi - lo) / 2;
    if (myb_time_compare(myb_time_from_seconds(s->table[3 * mid + 1]), t) <= 0) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  const double* row = s->table + 3 * lo;
  const double since = myb_dd_to_double(myb_time_diff(t, myb_time_from_seconds(row[1])));
  const double value = row[0] + since / row[2] * per_count(s);
  if (!(value <= s->part_total[s->partitions - 1])) {
    return MYB_E_RANGE;
  }
  *ticks = value;
  return MYB_OK;
}

static int is_delimiter(char c) {
  return c == '.' || c == ':' || c == '-' || c == ',' || c == ' ';
}

/* A whole number of at most MAX_DIGITS digits at text[*at]; 0 if there is none. */
static int whole_number(const char* text, size_t length, size_t* at, double* value) {
  double v = 0.0;
  int digits = 0;
  /* bound: the characters of the text */
  for (; *at < length && text[*at] >= '0' && text[*at] <= '9'; ++*at) {
    v = 10.0 * v + (double)(text[*at] - '0');
    if (++digits > MAX_DIGITS) {
      return 0;
    }
  }
  *value = v;
  return digits > 0;
}

myb_status myb_sclk_encode(const myb_sclk* sclk, const char* text, size_t length,
                           double* ticks) {
  MYB_CHECK(sclk != NULL && text != NULL && ticks != NULL);
  const myb_sclk* s = sclk;
  size_t at = 0;
  size_t end = length;
  /* bound: the characters of the text */
  while (at < end && text[at] == ' ') {
    ++at;
  }
  while (end > at && text[end - 1] == ' ') {
    --end;
  }
  /* the partition, if the string names one */
  size_t partition = 0;
  for (size_t i = at; i < end; ++i) {
    if (text[i] == '/') {
      size_t p = at;
      double number = 0.0;
      if (!whole_number(text, i, &p, &number) || p != i) {
        return MYB_E_FORMAT;
      }
      if (!(number >= 1.0) || !(number <= (double)s->partitions)) {
        return MYB_E_DOMAIN;
      }
      partition = (size_t)number;
      at = i + 1;
      break;
    }
  }
  while (at < end && text[at] == ' ') {
    ++at;
  }
  double raw = 0.0;
  int field = 0;
  /* bound: the fields of the clock */
  for (; field < s->fields && at < end; ++field) {
    double value = 0.0;
    if (!whole_number(text, end, &at, &value)) {
      return MYB_E_FORMAT;
    }
    raw += (value - s->offsets[field]) * s->ticks_per[field];
    if (at < end) {
      if (!is_delimiter(text[at])) {
        return MYB_E_FORMAT;
      }
      ++at;
      /* blanks may follow a delimiter */
      while (at < end && text[at] == ' ') {
        ++at;
      }
      if (at == end) {
        return MYB_E_FORMAT; /* a delimiter with nothing after it */
      }
    }
  }
  if (field == 0 || at != end) {
    return MYB_E_FORMAT;
  }
  if (partition == 0) {
    /* bound: the partitions; the first that holds the count */
    for (size_t p = 0; p < s->partitions && partition == 0; ++p) {
      if (raw >= s->part_start[p] && raw <= s->part_end[p]) {
        partition = p + 1;
      }
    }
    if (partition == 0) {
      return MYB_E_DOMAIN;
    }
  }
  const size_t p = partition - 1;
  if (!(raw >= s->part_start[p]) || !(raw <= s->part_end[p])) {
    return MYB_E_DOMAIN;
  }
  *ticks = raw - s->part_start[p] + (p > 0 ? s->part_total[p - 1] : 0.0);
  return MYB_OK;
}

myb_status myb_sclk_decode(const myb_sclk* sclk, double ticks, char* text, size_t size,
                           size_t* length) {
  MYB_CHECK(sclk != NULL && text != NULL);
  const myb_sclk* s = sclk;
  if (!(ticks >= 0.0) || !(ticks <= s->part_total[s->partitions - 1])) {
    return MYB_E_RANGE;
  }
  /* a count at the end of a partition is the start of the next one */
  size_t p = 0;
  /* bound: the partitions */
  while (p + 1 < s->partitions && s->part_total[p] <= ticks) {
    ++p;
  }
  double raw = nearest(ticks - (p > 0 ? s->part_total[p - 1] : 0.0) + s->part_start[p]);
  char out[32 + MAX_FIELDS * (MAX_DIGITS + 1)];
  size_t n = (size_t)snprintf(out, sizeof out, "%u/", (unsigned)(p + 1));
  for (int i = 0; i < s->fields; ++i) {
    const double value = floor(raw / s->ticks_per[i]);
    raw -= value * s->ticks_per[i];
    n += (size_t)snprintf(out + n, sizeof out - n, "%0*.0f", s->width[i], value + s->offsets[i]);
    if (i + 1 < s->fields) {
      out[n++] = s->delimiter;
    }
  }
  out[n] = '\0';
  if (length != NULL) {
    *length = n;
  }
  if (n + 1 > size) {
    return MYB_E_OVERFLOW;
  }
  memcpy(text, out, n + 1);
  return MYB_OK;
}
