/* SPDX-License-Identifier: BSD-3-Clause */
/* Time scale systems: a tree of scales defined by kernel variables, conversions along the
 * tree, and UTC labels (leap seconds; pre-1972 rate offsets and fractional steps).
 *
 * UTC model (the one of IAU SOFA iauDtf2d/iauUtctai, so that results can be compared):
 * for the UTC day with Modified Julian Date D,
 *   dat0  = TAI - UTC at 0h of D       (row in force on D)
 *   dlod  = drift rate of that row      (s/day; 0 since 1972)
 *   dleap = step of TAI - UTC at the end of D (row of D+1 minus row of D at 0h of D+1)
 * the day has 86400 + dleap label seconds (the last minute 60 + dleap), and the label u
 * seconds after 0h is the instant  TAI = 0h(D) + dat0 + u * (1 + dlod / 86400). */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/time.h"
#include "num/series.h"
#include "time/calendar.h"

enum {
  MAX_SCALES = 32,
  SCALE_NAME_MAX = 16,
  SERIES_NAME_MAX = 15,
  VAR_NAME_SIZE = 64,
  SERIES_ITERATIONS = 1, /* |dS/dX| < 1e-9, |S| < 2e-3 s: the first value is within 2e-12 s,
                            one fixed-point step within 2e-21 s (resolution of myb_time 1e-16 s) */
  MJD_J2000_DAY = 51544, /* J2000.0 = MJD 51544.5 */
  MAX_NODES = 1 << 26    /* of a tabulated relation */
};

typedef enum relation { REL_ROOT, REL_OFFSET, REL_RATE, REL_RATE_INVERSE, REL_SERIES } relation;

typedef struct scale_def {
  char name[SCALE_NAME_MAX + 1];
  int parent; /* -1 for the root */
  int depth;
  relation kind;
  double offset;
  double rate;
  myb_dd one_minus_rate;
  myb_time t0;
  double series_unit;
  myb_internal_series series;
  double* owned_terms;
  /* a SERIES relation given by a function in place of the series, and its values at equal
   * steps (myb_timesys_set_difference, myb_timesys_tabulate): nodes -1 .. count of the
   * table, table[0] the node before the start */
  myb_scale_difference function;
  void* function_context;
  double* table;
  size_t intervals;
  myb_time table_start;
  double table_step;
} scale_def;

typedef struct utc_row {
  int64_t start_mjd;
  double offset;
  double drift_mjd;
  double drift_rate;
} utc_row;

struct myb_timesys {
  myb_allocator allocator;
  int count;
  scale_def scale[MAX_SCALES];
  utc_row* utc;
  size_t utc_count;
  int utc_base;
  int64_t utc_valid_until;
};

/* ---- Reading kernel variables ---------------------------------------------------------- */

typedef struct reader {
  const myb_kpool* pool;
  myb_diag* diag;
} reader;

static myb_status missing_or_bad(reader* r, const char* name, myb_status status, const char* what) {
  myb_internal_diag(r->diag, 0, "%s: %s", name,
                    status == MYB_E_NOT_FOUND ? "variable not found" : what);
  return status == MYB_E_NOT_FOUND ? MYB_E_NOT_FOUND : MYB_E_FORMAT;
}

static myb_status read_numbers(reader* r, const char* name, size_t expected, const double** values,
                               size_t* count) {
  const myb_status status = myb_kpool_numbers(r->pool, name, values, count);
  if (status != MYB_OK) {
    return missing_or_bad(r, name, status, "expected numbers");
  }
  if (expected != 0 && *count != expected) {
    myb_internal_diag(r->diag, 0, "%s: expected %zu value(s), found %zu", name, expected, *count);
    return MYB_E_FORMAT;
  }
  return MYB_OK;
}

static myb_status read_number(reader* r, const char* name, double* value) {
  const double* values = NULL;
  size_t count = 0;
  const myb_status status = read_numbers(r, name, 1, &values, &count);
  if (status == MYB_OK) {
    *value = values[0];
  }
  return status;
}

/* String value copied into out (NUL-terminated); MYB_E_FORMAT if longer than max. */
static myb_status read_string(reader* r, const char* name, size_t index, char* out, size_t max) {
  const char* value = NULL;
  size_t length = 0;
  const myb_status status = myb_kpool_string(r->pool, name, index, &value, &length);
  if (status != MYB_OK) {
    return missing_or_bad(r, name, status, "expected a string");
  }
  if (length > max) {
    myb_internal_diag(r->diag, 0, "%s: '%.*s' is longer than %zu characters", name, (int)length,
                      value, max);
    return MYB_E_FORMAT;
  }
  memcpy(out, value, length);
  out[length] = '\0';
  return MYB_OK;
}

static int valid_identifier(const char* s) {
  if (s[0] == '\0') {
    return 0;
  }
  for (const char* c = s; *c != '\0'; ++c) {
    if (!((*c >= 'A' && *c <= 'Z') || (*c >= '0' && *c <= '9') || *c == '_')) {
      return 0;
    }
  }
  return 1;
}

static void var_name(char* out, const char* prefix, const char* name, const char* suffix) {
  const int n = snprintf(out, VAR_NAME_SIZE, "%s%s%s", prefix, name, suffix);
  MYB_CHECK(n > 0 && n < VAR_NAME_SIZE);
}

static int find_scale(const myb_timesys* s, const char* name, int limit) {
  for (int i = 0; i < limit; ++i) {
    if (strcmp(s->scale[i].name, name) == 0) {
      return i;
    }
  }
  return -1;
}

/* ---- Building ---------------------------------------------------------------------------- */

static myb_status read_rate(reader* r, scale_def* d) {
  char name[VAR_NAME_SIZE];
  var_name(name, "MYB_TIME_", d->name, "_RATE");
  myb_status status = read_number(r, name, &d->rate);
  if (status != MYB_OK) {
    return status;
  }
  if (!(d->rate > -1e-3 && d->rate < 1e-3)) {
    myb_internal_diag(r->diag, 0, "%s: rate outside (-1e-3, 1e-3)", name);
    return MYB_E_FORMAT;
  }
  double err;
  const double hi = myb_two_sum(1.0, -d->rate, &err);
  d->one_minus_rate = (myb_dd){.hi = hi, .lo = err};
  const double* jd = NULL;
  size_t count = 0;
  var_name(name, "MYB_TIME_", d->name, "_T0_JD");
  status = read_numbers(r, name, 2, &jd, &count);
  if (status != MYB_OK) {
    return status;
  }
  d->t0 = myb_time_from_jd(jd[0], jd[1]);
  var_name(name, "MYB_TIME_", d->name, "_OFFSET");
  return read_number(r, name, &d->offset);
}

static myb_status read_series(reader* r, const myb_allocator* a, scale_def* d) {
  char name[VAR_NAME_SIZE];
  char series[SERIES_NAME_MAX + 1];
  var_name(name, "MYB_TIME_", d->name, "_SERIES");
  myb_status status = read_string(r, name, 0, series, SERIES_NAME_MAX);
  if (status != MYB_OK) {
    return status;
  }
  if (!valid_identifier(series)) {
    myb_internal_diag(r->diag, 0, "%s: '%s' is not a valid series name", name, series);
    return MYB_E_FORMAT;
  }
  var_name(name, "MYB_SERIES_", series, "_UNIT");
  status = read_number(r, name, &d->series_unit);
  if (status != MYB_OK) {
    return status;
  }
  if (!(d->series_unit > 0.0)) {
    myb_internal_diag(r->diag, 0, "%s: the unit must be positive", name);
    return MYB_E_FORMAT;
  }
  const double* terms = NULL;
  size_t count = 0;
  var_name(name, "MYB_SERIES_", series, "_TERMS");
  status = read_numbers(r, name, 0, &terms, &count);
  if (status != MYB_OK) {
    return status;
  }
  if (count % 4 != 0 || !myb_internal_series_valid(terms, count / 4)) {
    myb_internal_diag(r->diag, 0,
                      "%s: expected groups of 4 numbers (k A w phi) with integer k in 0..%d", name,
                      MYB_SERIES_MAX_POWER);
    return MYB_E_FORMAT;
  }
  void* block = NULL;
  status = myb_internal_alloc_array(a, count, sizeof(double), &block);
  if (status != MYB_OK) {
    return status;
  }
  d->owned_terms = block;
  memcpy(d->owned_terms, terms, count * sizeof(double));
  d->series = (myb_internal_series){.count = count / 4, .terms = d->owned_terms};
  return MYB_OK;
}

static myb_status read_scale(reader* r, myb_timesys* s, int index) {
  scale_def* d = &s->scale[index];
  char name[VAR_NAME_SIZE];
  char text[SCALE_NAME_MAX + 1];
  var_name(name, "MYB_TIME_", d->name, "_FROM");
  myb_status status = read_string(r, name, 0, text, SCALE_NAME_MAX);
  if (status != MYB_OK) {
    return status;
  }
  d->parent = find_scale(s, text, index);
  if (d->parent < 0) {
    myb_internal_diag(r->diag, 0, "%s: '%s' is not a scale listed before %s", name, text, d->name);
    return MYB_E_FORMAT;
  }
  d->depth = s->scale[d->parent].depth + 1;
  var_name(name, "MYB_TIME_", d->name, "_KIND");
  status = read_string(r, name, 0, text, SCALE_NAME_MAX);
  if (status != MYB_OK) {
    return status;
  }
  if (strcmp(text, "OFFSET") == 0) {
    d->kind = REL_OFFSET;
    var_name(name, "MYB_TIME_", d->name, "_OFFSET");
    return read_number(r, name, &d->offset);
  }
  if (strcmp(text, "RATE") == 0 || strcmp(text, "RATE_INVERSE") == 0) {
    d->kind = text[4] == '\0' ? REL_RATE : REL_RATE_INVERSE;
    return read_rate(r, d);
  }
  if (strcmp(text, "SERIES") == 0) {
    d->kind = REL_SERIES;
    return read_series(r, &s->allocator, d);
  }
  myb_internal_diag(r->diag, 0, "%s: unknown kind '%s'", name, text);
  return MYB_E_FORMAT;
}

static myb_status read_utc(reader* r, myb_timesys* s) {
  if (myb_kpool_describe(r->pool, "MYB_UTC_START_MJD", NULL, NULL) == MYB_E_NOT_FOUND) {
    s->utc_base = -1;
    return MYB_OK;
  }
  char base[SCALE_NAME_MAX + 1];
  myb_status status = read_string(r, "MYB_UTC_BASE", 0, base, SCALE_NAME_MAX);
  if (status != MYB_OK) {
    return status;
  }
  s->utc_base = find_scale(s, base, s->count);
  if (s->utc_base < 0) {
    myb_internal_diag(r->diag, 0, "MYB_UTC_BASE: '%s' is not a scale of the system", base);
    return MYB_E_FORMAT;
  }
  const double* start = NULL;
  size_t n = 0;
  status = read_numbers(r, "MYB_UTC_START_MJD", 0, &start, &n);
  const char* columns[] = {"MYB_UTC_OFFSET", "MYB_UTC_DRIFT_MJD", "MYB_UTC_DRIFT_RATE"};
  const double* column[3] = {NULL, NULL, NULL};
  for (int c = 0; c < 3 && status == MYB_OK; ++c) {
    size_t count = 0;
    status = read_numbers(r, columns[c], n, &column[c], &count);
  }
  double valid_until = 0.0;
  if (status == MYB_OK) {
    status = read_number(r, "MYB_UTC_VALID_UNTIL_MJD", &valid_until);
  }
  if (status != MYB_OK) {
    return status;
  }
  for (size_t i = 0; i < n; ++i) {
    const int integral = start[i] == (double)(int64_t)start[i];
    if (!integral || (i > 0 && !(start[i] > start[i - 1]))) {
      myb_internal_diag(r->diag, 0, "MYB_UTC_START_MJD: value %zu is not an increasing integer",
                        i + 1);
      return MYB_E_FORMAT;
    }
  }
  if (valid_until != (double)(int64_t)valid_until || !(valid_until > start[n - 1])) {
    myb_internal_diag(r->diag, 0, "MYB_UTC_VALID_UNTIL_MJD: must be an integer after the last row");
    return MYB_E_FORMAT;
  }
  void* block = NULL;
  status = myb_internal_alloc_array(&s->allocator, n, sizeof(utc_row), &block);
  if (status != MYB_OK) {
    return status;
  }
  s->utc = block;
  s->utc_count = n;
  for (size_t i = 0; i < n; ++i) {
    s->utc[i] = (utc_row){.start_mjd = (int64_t)start[i],
                          .offset = column[0][i],
                          .drift_mjd = column[1][i],
                          .drift_rate = column[2][i]};
  }
  s->utc_valid_until = (int64_t)valid_until;
  return MYB_OK;
}

static myb_status build(reader* r, myb_timesys* s) {
  myb_kvar_type type = 0;
  size_t n = 0;
  myb_status status = myb_kpool_describe(r->pool, "MYB_TIME_SCALES", &type, &n);
  if (status != MYB_OK || type != MYB_KVAR_STRINGS) {
    return missing_or_bad(r, "MYB_TIME_SCALES", status, "expected a list of scale names");
  }
  if (n > MAX_SCALES) {
    myb_internal_diag(r->diag, 0, "MYB_TIME_SCALES: more than %d scales", MAX_SCALES);
    return MYB_E_FORMAT;
  }
  for (size_t i = 0; i < n; ++i) {
    scale_def* d = &s->scale[i];
    status = read_string(r, "MYB_TIME_SCALES", i, d->name, SCALE_NAME_MAX);
    if (status != MYB_OK) {
      return status;
    }
    if (!valid_identifier(d->name) || find_scale(s, d->name, (int)i) >= 0) {
      myb_internal_diag(r->diag, 0, "MYB_TIME_SCALES: '%s' is not a new name of [A-Z0-9_]",
                        d->name);
      return MYB_E_FORMAT;
    }
    s->count = (int)i + 1;
    if (i == 0) {
      d->kind = REL_ROOT;
      d->parent = -1;
      d->depth = 0;
      continue;
    }
    status = read_scale(r, s, (int)i);
    if (status != MYB_OK) {
      return status;
    }
  }
  return read_utc(r, s);
}

myb_status myb_timesys_create(const myb_kpool* pool, const myb_allocator* allocator,
                              myb_timesys** system, myb_diag* diag) {
  MYB_CHECK(pool != NULL && system != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_timesys), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_timesys* s = block;
  memset(s, 0, sizeof *s); /* no compound literal: the structure is 3.5 KB */
  s->allocator = a;
  s->utc_base = -1;
  reader r = {.pool = pool, .diag = diag};
  status = build(&r, s);
  if (status != MYB_OK) {
    myb_timesys_destroy(s);
    return status;
  }
  *system = s;
  return MYB_OK;
}

void myb_timesys_destroy(myb_timesys* system) {
  if (system == NULL) {
    return;
  }
  const myb_allocator a = system->allocator;
  for (int i = 0; i < MAX_SCALES; ++i) {
    myb_internal_release(&a, system->scale[i].owned_terms);
    myb_internal_release(&a, system->scale[i].table);
  }
  myb_internal_release(&a, system->utc);
  myb_internal_release(&a, system);
}

myb_status myb_timesys_scale(const myb_timesys* system, const char* name, myb_scale* scale) {
  MYB_CHECK(system != NULL && name != NULL && scale != NULL);
  const int index = find_scale(system, name, system->count);
  if (index < 0) {
    return MYB_E_NOT_FOUND;
  }
  *scale = (myb_scale){.index = index};
  return MYB_OK;
}

static void check_scale(const myb_timesys* system, myb_scale scale) {
  MYB_CHECK(system != NULL && scale.index >= 0 && scale.index < system->count);
}

const char* myb_timesys_scale_name(const myb_timesys* system, myb_scale scale) {
  check_scale(system, scale);
  return system->scale[scale.index].name;
}

/* ---- Relations ----------------------------------------------------------------------------- */

/* X = F - rate (F - T0) + offset */
static myb_time apply_rate(const scale_def* d, myb_time f) {
  const myb_dd delta = myb_time_diff(f, d->t0);
  return myb_time_add_dd(f, myb_dd_add_d(myb_dd_mul_d(delta, -d->rate), d->offset));
}

/* F from X = F - rate (F - T0) + offset:  F = T0 + (X - offset - T0) / (1 - rate) */
static myb_time solve_rate(const scale_def* d, myb_time x) {
  const myb_dd delta = myb_dd_add_d(myb_time_diff(x, d->t0), -d->offset);
  return myb_time_add_dd(d->t0, myb_dd_div(delta, d->one_minus_rate));
}

/* The difference scale - parent of a SERIES relation at the instant x of the scale, from
 * what defines it: the function if one was set, the series otherwise. */
static myb_status defined_difference(const scale_def* d, myb_time x, double* value) {
  if (d->function != NULL) {
    const myb_status status = d->function(d->function_context, x, value);
    if (status == MYB_OK && !isfinite(*value)) {
      return MYB_E_DOMAIN;
    }
    return status;
  }
  *value = myb_internal_series_eval(&d->series, myb_time_seconds(x) / d->series_unit);
  return MYB_OK;
}

/* The same through the table where it covers x: the cubic through the four nodes around. */
static myb_status difference_at(const scale_def* d, myb_time x, double* value) {
  if (d->table != NULL) {
    const double h = d->table_step;
    const double since = myb_dd_to_double(myb_time_diff(x, d->table_start));
    if (since >= 0.0 && since <= h * (double)d->intervals) {
      double k = floor(since / h);
      if (k > (double)(d->intervals - 1)) {
        k = (double)(d->intervals - 1);
      }
      const double* f = d->table + (size_t)k; /* the nodes k - 1 .. k + 2 */
      const double s = since / h - k;
      *value = -s * (s - 1.0) * (s - 2.0) / 6.0 * f[0] +
               (s + 1.0) * (s - 1.0) * (s - 2.0) / 2.0 * f[1] -
               (s + 1.0) * s * (s - 2.0) / 2.0 * f[2] + (s + 1.0) * s * (s - 1.0) / 6.0 * f[3];
      return MYB_OK;
    }
  }
  return defined_difference(d, x, value);
}

/* Instant on scale d from the instant on its parent. */
static myb_status down(const scale_def* d, myb_time parent_time, myb_time* out) {
  switch (d->kind) {
    case REL_OFFSET:
      *out = myb_time_add(parent_time, d->offset);
      return MYB_OK;
    case REL_RATE:
      *out = apply_rate(d, parent_time);
      return MYB_OK;
    case REL_RATE_INVERSE:
      *out = solve_rate(d, parent_time);
      return MYB_OK;
    case REL_SERIES: {
      double difference = 0.0;
      myb_status status = difference_at(d, parent_time, &difference);
      myb_time x = myb_time_add(parent_time, difference);
      for (int i = 0; i < SERIES_ITERATIONS && status == MYB_OK; ++i) {
        status = difference_at(d, x, &difference);
        x = myb_time_add(parent_time, difference);
      }
      *out = x;
      return status;
    }
    case REL_ROOT:
      break;
  }
  MYB_UNREACHABLE();
}

/* Instant on the parent of d from the instant on d. */
static myb_status up(const scale_def* d, myb_time t, myb_time* out) {
  switch (d->kind) {
    case REL_OFFSET:
      *out = myb_time_add(t, -d->offset);
      return MYB_OK;
    case REL_RATE:
      *out = solve_rate(d, t);
      return MYB_OK;
    case REL_RATE_INVERSE:
      *out = apply_rate(d, t);
      return MYB_OK;
    case REL_SERIES: {
      double difference = 0.0;
      const myb_status status = difference_at(d, t, &difference);
      *out = myb_time_add(t, -difference);
      return status;
    }
    case REL_ROOT:
      break;
  }
  MYB_UNREACHABLE();
}

myb_status myb_timesys_set_difference(myb_timesys* system, myb_scale scale,
                                      myb_scale_difference function, void* context) {
  check_scale(system, scale);
  scale_def* d = &system->scale[scale.index];
  if (d->kind != REL_SERIES) {
    return MYB_E_DOMAIN;
  }
  /* a table holds the values of what defined the relation when it was made */
  myb_internal_release(&system->allocator, d->table);
  d->table = NULL;
  d->intervals = 0;
  d->function = function;
  d->function_context = context;
  return MYB_OK;
}

myb_status myb_timesys_tabulate(myb_timesys* system, myb_scale scale, myb_time t0, myb_time t1,
                                double step) {
  check_scale(system, scale);
  scale_def* d = &system->scale[scale.index];
  const double h = step > 0.0 ? step : 3600.0;
  const double span = myb_dd_to_double(myb_time_diff(t1, t0));
  if (d->kind != REL_SERIES || !(step >= 0.0) || !isfinite(step) || !(span > 0.0) ||
      !(span / h < (double)MAX_NODES - 4.0)) {
    return MYB_E_DOMAIN;
  }
  const size_t intervals = (size_t)ceil(span / h);
  void* block = NULL;
  myb_status status =
      myb_internal_alloc_array(&system->allocator, intervals + 3, sizeof(double), &block);
  if (status != MYB_OK) {
    return status;
  }
  double* table = block;
  for (size_t k = 0; k < intervals + 3 && status == MYB_OK; ++k) {
    status = defined_difference(d, myb_time_add(t0, h * ((double)k - 1.0)), &table[k]);
  }
  if (status != MYB_OK) {
    myb_internal_release(&system->allocator, table);
    return status;
  }
  myb_internal_release(&system->allocator, d->table);
  d->table = table;
  d->intervals = intervals;
  d->table_start = t0;
  d->table_step = h;
  return MYB_OK;
}

myb_status myb_timesys_convert(const myb_timesys* system, myb_time t, myb_scale from, myb_scale to,
                               myb_time* out) {
  check_scale(system, from);
  check_scale(system, to);
  MYB_CHECK(out != NULL);
  const scale_def* sc = system->scale;
  int a = from.index;
  int b = to.index;
  int path[MAX_SCALES];
  int path_length = 0;
  myb_time x = t;
  /* bound: every step moves one level towards the root; depth < MAX_SCALES */
  for (int step = 0; step < 2 * MAX_SCALES && a != b; ++step) {
    if (sc[a].depth >= sc[b].depth) {
      const myb_status status = up(&sc[a], x, &x);
      if (status != MYB_OK) {
        return status;
      }
      a = sc[a].parent;
    } else {
      path[path_length++] = b;
      b = sc[b].parent;
    }
  }
  MYB_CHECK(a == b);
  for (int i = path_length - 1; i >= 0; --i) {
    const myb_status status = down(&sc[path[i]], x, &x);
    if (status != MYB_OK) {
      return status;
    }
  }
  *out = x;
  return MYB_OK;
}

/* ---- UTC ------------------------------------------------------------------------------------ */

typedef struct utc_day {
  double dat0;  /* TAI - UTC at 0h */
  double dlod;  /* s/day */
  double dleap; /* step at the end of the day */
} utc_day;

static int utc_row_index(const myb_timesys* s, int64_t mjd) {
  int found = -1;
  for (size_t i = 0; i < s->utc_count && s->utc[i].start_mjd <= mjd; ++i) {
    found = (int)i;
  }
  return found;
}

static double utc_dat(const utc_row* row, double mjd) {
  return row->offset + (mjd - row->drift_mjd) * row->drift_rate;
}

static myb_status utc_day_of(const myb_timesys* s, int64_t mjd, utc_day* day) {
  const int row = utc_row_index(s, mjd);
  if (row < 0 || mjd >= s->utc_valid_until) {
    return MYB_E_RANGE;
  }
  const utc_row* r = &s->utc[row];
  day->dat0 = utc_dat(r, (double)mjd);
  day->dlod = r->drift_rate;
  day->dleap = 0.0;
  if (mjd + 1 < s->utc_valid_until) {
    const int next = utc_row_index(s, mjd + 1);
    if (next != row) {
      day->dleap = utc_dat(&s->utc[next], (double)(mjd + 1)) - utc_dat(r, (double)(mjd + 1));
    }
  }
  return MYB_OK;
}

/* Base-scale instant of 0h UTC of day mjd plus u label seconds. */
static myb_time utc_instant(int64_t mjd, const utc_day* day, myb_dd u) {
  const int64_t midnight = (mjd - MJD_J2000_DAY) * 86400 - 43200;
  myb_time t = myb_time_add((myb_time){.s = midnight, .f = 0.0}, day->dat0);
  t = myb_time_add_dd(t, u);
  if (day->dlod != 0.0) {
    t = myb_time_add(t, u.hi * (day->dlod / 86400.0));
  }
  return t;
}

myb_status myb_timesys_from_utc(const myb_timesys* system, const myb_calendar* utc, myb_scale to,
                                myb_time* out) {
  MYB_CHECK(system != NULL && utc != NULL && out != NULL);
  check_scale(system, to);
  if (system->utc_base < 0) {
    return MYB_E_NOT_FOUND;
  }
  const myb_calendar c = *utc;
  if (!myb_internal_valid_date(c.year, c.month, c.day) || c.hour < 0 || c.hour > 23 ||
      c.minute < 0 || c.minute > 59) {
    return MYB_E_DOMAIN;
  }
  const int64_t mjd = myb_internal_days_from_civil(c.year, c.month, c.day) + MYB_MJD_OF_UNIX_EPOCH;
  utc_day day;
  const myb_status status = utc_day_of(system, mjd, &day);
  if (status != MYB_OK) {
    return status;
  }
  const double limit = (c.hour == 23 && c.minute == 59) ? 60.0 + day.dleap : 60.0;
  if (!(c.second >= 0.0 && c.second < limit)) {
    return MYB_E_DOMAIN;
  }
  double err;
  const double u_hi = myb_two_sum((double)(c.hour * 3600 + c.minute * 60), c.second, &err);
  const myb_time base_time = utc_instant(mjd, &day, (myb_dd){.hi = u_hi, .lo = err});
  return myb_timesys_convert(system, base_time, (myb_scale){.index = system->utc_base}, to, out);
}

myb_status myb_timesys_to_utc(const myb_timesys* system, myb_time t, myb_scale from,
                              myb_calendar* utc) {
  MYB_CHECK(system != NULL && utc != NULL);
  check_scale(system, from);
  if (system->utc_base < 0) {
    return MYB_E_NOT_FOUND;
  }
  myb_time base_time;
  myb_status status =
      myb_timesys_convert(system, t, from, (myb_scale){.index = system->utc_base}, &base_time);
  if (status != MYB_OK) {
    return status;
  }
  /* TAI - UTC is positive and below a day: the UTC day is the base day or the one before. */
  const int64_t base_day = myb_internal_floor_div(base_time.s + 43200, 86400) + MJD_J2000_DAY;
  int64_t mjd = 0;
  utc_day day = {};
  myb_time midnight = {};
  status = MYB_E_RANGE;
  for (int64_t candidate = base_day; candidate >= base_day - 1 && status != MYB_OK; --candidate) {
    utc_day d;
    if (utc_day_of(system, candidate, &d) != MYB_OK) {
      continue;
    }
    const myb_time start = utc_instant(candidate, &d, (myb_dd){.hi = 0.0, .lo = 0.0});
    if (myb_time_compare(base_time, start) >= 0) {
      mjd = candidate;
      day = d;
      midnight = start;
      status = MYB_OK;
    }
  }
  if (status != MYB_OK) {
    return MYB_E_RANGE;
  }
  myb_dd u = myb_time_diff(base_time, midnight);
  if (day.dlod != 0.0) {
    u = myb_dd_div(u, myb_dd_add_d(myb_dd_from_double(1.0), day.dlod / 86400.0));
  }
  /* Whole label seconds and fraction. */
  double whole = __builtin_floor(u.hi);
  myb_dd fraction = myb_dd_add_d(u, -whole);
  if (fraction.hi < 0.0) {
    whole -= 1.0;
    fraction = myb_dd_add_d(fraction, 1.0);
  }
  const int64_t seconds = (int64_t)whole;
  myb_calendar c = {};
  myb_internal_civil_from_days(mjd - MYB_MJD_OF_UNIX_EPOCH, &c.year, &c.month, &c.day);
  double limit = 60.0;
  if (seconds < 86340) {
    c.hour = (int)(seconds / 3600);
    c.minute = (int)(seconds % 3600 / 60);
    c.second = (double)(seconds % 60) + myb_dd_to_double(fraction);
  } else {
    c.hour = 23;
    c.minute = 59;
    c.second = (double)(seconds - 86340) + myb_dd_to_double(fraction);
    limit = 60.0 + day.dleap;
  }
  if (c.second >= limit) {
    /* Rounding of the fraction, or the 1e-9 s discontinuity of the pre-1972 model on days
     * with both a rate and a step: keep the label valid with the largest double below the
     * limit (limit > 0, so the encoding minus one is its predecessor). */
    uint64_t bits;
    memcpy(&bits, &limit, sizeof bits);
    bits -= 1;
    memcpy(&c.second, &bits, sizeof bits);
  }
  *utc = c;
  return MYB_OK;
}
