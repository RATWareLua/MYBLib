/* SPDX-License-Identifier: BSD-3-Clause */
/* Space weather (myb/weather.h): a table of days read by the names of its columns. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/io.h"
#include "myb/time.h"
#include "myb/weather.h"
#include "platform/platform.h"

enum {
  MAX_COLUMNS = 96,
  MAX_FIELD = 31,
  /* what a day keeps, as indices of its values */
  KP = 0,            /* 8 */
  AP = 8,            /* 8 */
  AP_DAILY = 16,
  F107_OBSERVED = 17,
  F107_ADJUSTED = 18,
  OBSERVED_CENTRED = 19,
  OBSERVED_LAST = 20,
  ADJUSTED_CENTRED = 21,
  ADJUSTED_LAST = 22,
  VALUES = 23,
  /* the columns that are not numbers */
  COLUMN_DATE = VALUES,
  COLUMN_KIND = VALUES + 1,
  KNOWN = VALUES + 2,
  SPANS = 8 /* of three hours in a day */
};

static const char* const column_names[KNOWN] = {
    "KP1", "KP2", "KP3", "KP4", "KP5", "KP6", "KP7", "KP8", "AP1", "AP2", "AP3", "AP4",
    "AP5", "AP6", "AP7", "AP8", "AP_AVG", "F10.7_OBS", "F10.7_ADJ", "F10.7_OBS_CENTER81",
    "F10.7_OBS_LAST81", "F10.7_ADJ_CENTER81", "F10.7_ADJ_LAST81", "DATE", "F10.7_DATA_TYPE"};

typedef struct day_row {
  int64_t mjd;
  int kind;
  double value[VALUES];
} day_row;

struct myb_space_weather {
  myb_allocator allocator;
  day_row* days;
  size_t count;
};

/* The field of a line that starts at *at: its bounds without the blanks around it; moves
 * *at behind the comma. 0 at the end of the line. */
static int next_field(const char* line, size_t length, size_t* at, const char** field,
                      size_t* size) {
  if (*at > length) {
    return 0;
  }
  size_t begin = *at;
  size_t end = begin;
  /* bound: the characters of the line */
  while (end < length && line[end] != ',') {
    ++end;
  }
  *at = end + 1;
  while (begin < end && (line[begin] == ' ' || line[begin] == '\t')) {
    ++begin;
  }
  while (end > begin && (line[end - 1] == ' ' || line[end - 1] == '\t')) {
    --end;
  }
  *field = line + begin;
  *size = end - begin;
  return 1;
}

static int digits(const char* s, size_t n, int64_t* value) {
  int64_t v = 0;
  for (size_t i = 0; i < n; ++i) {
    if (s[i] < '0' || s[i] > '9') {
      return 0;
    }
    v = 10 * v + (s[i] - '0');
  }
  *value = v;
  return n > 0;
}

/* YYYY-MM-DD -> the Modified Julian Date of the day. */
static int parse_date(const char* s, size_t n, int64_t* mjd) {
  int64_t year = 0;
  int64_t month = 0;
  int64_t day = 0;
  if (n != 10 || s[4] != '-' || s[7] != '-' || !digits(s, 4, &year) ||
      !digits(s + 5, 2, &month) || !digits(s + 8, 2, &day)) {
    return 0;
  }
  const myb_calendar calendar = {.year = year, .month = (int)month, .day = (int)day};
  myb_time t;
  if (myb_calendar_to_time(&calendar, &t) != MYB_OK) {
    return 0;
  }
  /* 0h of the day is half a day before a noon: J2000.0 is MJD 51544.5 */
  *mjd = 51544 + (t.s + 43200) / 86400;
  return (t.s + 43200) % 86400 == 0;
}

static myb_status fail(myb_diag* diag, long line, const char* what) {
  myb_internal_diag(diag, line, "line %ld: %s", line, what);
  return MYB_E_FORMAT;
}

myb_status myb_space_weather_load(const char* text, size_t size, const myb_allocator* allocator,
                                  myb_space_weather** weather, myb_diag* diag) {
  MYB_CHECK((text != NULL || size == 0) && weather != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  /* a day a line at most */
  size_t lines = 1;
  for (size_t i = 0; i < size; ++i) {
    lines += text[i] == '\n';
  }
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, lines, sizeof(day_row), &block);
  if (status != MYB_OK) {
    return status;
  }
  day_row* days = block;
  size_t count = 0;
  int column_of[MAX_COLUMNS]; /* what the column of the file is: an index of column_names */
  int columns = -1;           /* -1: the names have not been read yet */
  long line_number = 0;
  size_t at = 0;
  /* bound: the lines of the text */
  while (at < size && status == MYB_OK) {
    size_t end = at;
    while (end < size && text[end] != '\n') {
      ++end;
    }
    const char* line = text + at;
    size_t length = end - at;
    at = end + 1;
    ++line_number;
    while (length > 0 && (line[length - 1] == '\r' || line[length - 1] == ' ')) {
      --length;
    }
    if (length == 0 || line[0] == '#') {
      continue;
    }
    size_t cursor = 0;
    const char* field = NULL;
    size_t n = 0;
    if (columns < 0) {
      columns = 0;
      int has_date = 0;
      /* bound: the fields of the line */
      while (next_field(line, length, &cursor, &field, &n)) {
        if (columns == MAX_COLUMNS) {
          status = fail(diag, line_number, "too many columns");
          break;
        }
        column_of[columns] = -1;
        for (int k = 0; k < KNOWN; ++k) {
          if (n == strlen(column_names[k]) && memcmp(field, column_names[k], n) == 0) {
            column_of[columns] = k;
            has_date = has_date || k == COLUMN_DATE;
          }
        }
        ++columns;
      }
      if (status == MYB_OK && !has_date) {
        status = fail(diag, line_number, "no column DATE in the first line");
      }
      continue;
    }
    day_row* d = &days[count];
    d->kind = MYB_SPACE_WEATHER_OBSERVED;
    for (int k = 0; k < VALUES; ++k) {
      d->value[k] = (double)NAN;
    }
    int dated = 0;
    int column = 0;
    /* bound: the fields of the line */
    for (; status == MYB_OK && next_field(line, length, &cursor, &field, &n); ++column) {
      const int what = column < columns ? column_of[column] : -1;
      if (what < 0 || n == 0) {
        continue;
      }
      if (what == COLUMN_DATE) {
        dated = parse_date(field, n, &d->mjd);
        if (!dated) {
          status = fail(diag, line_number, "a date that is not YYYY-MM-DD");
        }
      } else if (what == COLUMN_KIND) {
        if (n == 3 && memcmp(field, "INT", 3) == 0) {
          d->kind = MYB_SPACE_WEATHER_INTERPOLATED;
        } else if (n == 3 && field[0] == 'P') {
          d->kind = MYB_SPACE_WEATHER_PREDICTED; /* PRD, PRM */
        }
      } else {
        double value = 0.0;
        size_t used = 0;
        if (myb_parse_double(field, n, &value, &used) != MYB_OK || used != n ||
            !isfinite(value)) {
          status = fail(diag, line_number, "a field that is not a number");
        }
        d->value[what] = what < KP + SPANS ? value / 10.0 : value; /* Kp comes in tenths */
      }
    }
    if (status == MYB_OK && !dated) {
      status = fail(diag, line_number, "a row without a date");
    }
    if (status == MYB_OK && count > 0 && d->mjd <= days[count - 1].mjd) {
      status = fail(diag, line_number, "days out of order");
    }
    ++count;
  }
  if (status == MYB_OK && count == 0) {
    status = fail(diag, line_number, "no days");
  }
  void* object = NULL;
  if (status == MYB_OK) {
    status = myb_internal_alloc_array(&a, 1, sizeof(myb_space_weather), &object);
  }
  if (status != MYB_OK) {
    myb_internal_release(&a, days);
    return status;
  }
  myb_space_weather* w = object;
  w->allocator = a;
  w->days = days;
  w->count = count;
  *weather = w;
  return MYB_OK;
}

myb_status myb_space_weather_load_file(const char* path, const myb_allocator* allocator,
                                       myb_space_weather** weather, myb_diag* diag) {
  MYB_CHECK(path != NULL && weather != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &a, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded =
      myb_space_weather_load((const char*)map.data, map.size, allocator, weather, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

void myb_space_weather_destroy(myb_space_weather* weather) {
  if (weather == NULL) {
    return;
  }
  const myb_allocator a = weather->allocator;
  myb_internal_release(&a, weather->days);
  myb_internal_release(&a, weather);
}

size_t myb_space_weather_days(const myb_space_weather* weather, int64_t* first, int64_t* last) {
  MYB_CHECK(weather != NULL);
  if (first != NULL) {
    *first = weather->days[0].mjd;
  }
  if (last != NULL) {
    *last = weather->days[weather->count - 1].mjd;
  }
  return weather->count;
}

/* The row of the day, or NULL. */
static const day_row* find(const myb_space_weather* w, int64_t mjd) {
  size_t lo = 0;
  size_t hi = w->count;
  /* bound: the interval halves */
  while (lo < hi) {
    const size_t mid = lo + (hi - lo) / 2;
    if (w->days[mid].mjd < mjd) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo < w->count && w->days[lo].mjd == mjd ? &w->days[lo] : NULL;
}

myb_status myb_space_weather_day_at(const myb_space_weather* weather, int64_t mjd,
                                    myb_space_weather_day* day) {
  MYB_CHECK(weather != NULL && day != NULL);
  const day_row* d = find(weather, mjd);
  if (d == NULL) {
    return MYB_E_RANGE;
  }
  day->mjd = d->mjd;
  memcpy(day->kp, d->value + KP, sizeof day->kp);
  memcpy(day->ap, d->value + AP, sizeof day->ap);
  day->ap_daily = d->value[AP_DAILY];
  day->f107_observed = d->value[F107_OBSERVED];
  day->f107_adjusted = d->value[F107_ADJUSTED];
  day->f107_observed_81[0] = d->value[OBSERVED_CENTRED];
  day->f107_observed_81[1] = d->value[OBSERVED_LAST];
  day->f107_adjusted_81[0] = d->value[ADJUSTED_CENTRED];
  day->f107_adjusted_81[1] = d->value[ADJUSTED_LAST];
  day->kind = d->kind;
  return MYB_OK;
}

/* ap of the three-hour span number `span` counted from the start of the day mjd (negative:
 * of the days before); a NaN where the table has none. */
static double ap_of_span(const myb_space_weather* w, int64_t mjd, int span) {
  int64_t day = mjd;
  int s = span;
  /* bound: at most eight days back */
  for (int k = 0; k < 8 && s < 0; ++k) {
    s += SPANS;
    --day;
  }
  const day_row* d = s >= 0 && s < SPANS ? find(w, day) : NULL;
  return d != NULL ? d->value[AP + s] : (double)NAN;
}

myb_status myb_space_weather_inputs_at(const myb_space_weather* weather, int64_t mjd,
                                       double seconds, int adjusted,
                                       myb_space_weather_inputs* inputs) {
  MYB_CHECK(weather != NULL && inputs != NULL);
  if (!(seconds >= 0.0) || !(seconds <= 86401.0)) {
    return MYB_E_DOMAIN;
  }
  const day_row* today = find(weather, mjd);
  const day_row* before = find(weather, mjd - 1);
  if (today == NULL || before == NULL) {
    return MYB_E_RANGE;
  }
  inputs->f107 = before->value[adjusted ? F107_ADJUSTED : F107_OBSERVED];
  inputs->f107_average = today->value[adjusted ? ADJUSTED_CENTRED : OBSERVED_CENTRED];
  inputs->ap = today->value[AP_DAILY];
  if (!isfinite(inputs->f107) || !isfinite(inputs->f107_average) || !isfinite(inputs->ap)) {
    return MYB_E_RANGE;
  }
  int span = (int)(seconds / 10800.0);
  span = span > SPANS - 1 ? SPANS - 1 : span; /* the leap second belongs to the last span */
  inputs->kp = today->value[KP + span];
  double history[7] = {inputs->ap, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int k = 0; k < 4; ++k) {
    history[1 + k] = ap_of_span(weather, mjd, span - k);
  }
  for (int k = 4; k < 12; ++k) {
    history[5] += ap_of_span(weather, mjd, span - k) / 8.0;
    history[6] += ap_of_span(weather, mjd, span - k - 8) / 8.0;
  }
  inputs->has_history = 1;
  for (int k = 0; k < 7; ++k) {
    inputs->has_history = inputs->has_history && isfinite(history[k]);
  }
  for (int k = 0; k < 7; ++k) {
    inputs->ap_history[k] = inputs->has_history ? history[k] : inputs->ap;
  }
  return MYB_OK;
}

/* ---- What a model of the upper atmosphere is asked with ------------------------------------- */

/* The Modified Julian Date of 0h of a calendar day. */
static myb_status day_number(int64_t year, int month, int day, int64_t* mjd) {
  const myb_calendar calendar = {.year = year, .month = month, .day = day};
  myb_time t;
  const myb_status status = myb_calendar_to_time(&calendar, &t);
  if (status == MYB_OK) {
    *mjd = 51544 + (t.s + 43200) / 86400;
  }
  return status;
}

myb_status myb_thermosphere_point_at(const myb_thermosphere_setting* setting, myb_time t,
                                     const double r[3], myb_thermosphere_point* point) {
  MYB_CHECK(setting != NULL && setting->timesys != NULL && r != NULL && point != NULL);
  myb_calendar utc;
  myb_status status = myb_timesys_to_utc(setting->timesys, t, setting->scale, &utc);
  int64_t today = 0;
  int64_t new_year = 0;
  if (status == MYB_OK) {
    status = day_number(utc.year, utc.month, utc.day, &today);
  }
  if (status == MYB_OK) {
    status = day_number(utc.year, 1, 1, &new_year);
  }
  double longitude = 0.0;
  double latitude = 0.0;
  double height = 0.0;
  if (status == MYB_OK) {
    status = myb_geodetic_from_cartesian(&setting->ellipsoid, r, &longitude, &latitude, &height);
  }
  if (status != MYB_OK) {
    return status;
  }
  const double seconds = 3600.0 * utc.hour + 60.0 * utc.minute + utc.second;
  myb_space_weather_inputs indices = {.f107 = setting->f107,
                                      .f107_average = setting->f107_average,
                                      .ap = setting->ap};
  for (int i = 0; i < 7; ++i) {
    indices.ap_history[i] = setting->ap;
  }
  if (setting->weather != NULL) {
    status = myb_space_weather_inputs_at(setting->weather, today, seconds, setting->adjusted,
                                         &indices);
    if (status != MYB_OK) {
      return status;
    }
  }
  point->mjd = today;
  point->day_of_year = (int)(today - new_year) + 1;
  point->seconds = seconds < 86400.0 ? seconds : 86399.999;
  point->height = height;
  point->latitude = latitude;
  point->longitude = longitude;
  point->local_solar_time = seconds / 3600.0 + longitude * 57.29577951308232 / 15.0;
  point->indices = indices;
  return MYB_OK;
}
