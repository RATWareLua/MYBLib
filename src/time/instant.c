/* SPDX-License-Identifier: BSD-3-Clause */
/* Instants (whole seconds + fraction), Julian dates and the proleptic Gregorian calendar.
 * Calendar arithmetic: H. Hinnant, "chrono-Compatible Low-Level Date Algorithms" (2013,
 * public domain): integer-only, valid for every year, no loops. */
#include <math.h>
#include <stdckdint.h>

#include "core/base.h"
#include "myb/time.h"
#include "num/eft.h"
#include "time/calendar.h"

static const double limit = 0x1p62;

/* s + hi + lo, normalised. hi and lo are any finite doubles with |hi|, |lo| < 2^62. */
static myb_time normalize(int64_t s, double hi, double lo) {
  MYB_CHECK(hi > -limit && hi < limit && lo > -limit && lo < limit);
  const double hi_whole = floor(hi);
  const double lo_whole = floor(lo);
  int64_t total = s;
  MYB_CHECK(!ckd_add(&total, total, (int64_t)hi_whole));
  MYB_CHECK(!ckd_add(&total, total, (int64_t)lo_whole));
  const double h = hi - hi_whole; /* exact, in [0, 1) */
  const double l = lo - lo_whole; /* exact, in [0, 1) */
  double f = h + l;               /* in [0, 2] after rounding */
  if (f >= 1.0) {
    MYB_CHECK(!ckd_add(&total, total, (int64_t)1));
    const double g = (h - 1.0) + l; /* h - 1 is exact */
    /* g < 0 only when h + l is just below 1 and rounded up to 1: the nearest instant with a
     * fraction in [0, 1) is the next whole second. */
    f = g >= 0.0 ? g : 0.0;
  }
  return (myb_time){.s = total, .f = f};
}

myb_time myb_time_from_seconds(double seconds) {
  return normalize(0, seconds, 0.0);
}

myb_time myb_time_from_dd(myb_dd seconds) {
  return normalize(0, seconds.hi, seconds.lo);
}

myb_time myb_time_add(myb_time t, double seconds) {
  double err;
  const double sum = myb_eft_two_sum(t.f, seconds, &err);
  return normalize(t.s, sum, err);
}

myb_time myb_time_add_dd(myb_time t, myb_dd seconds) {
  const myb_dd sum = myb_dd_add_d(seconds, t.f);
  return normalize(t.s, sum.hi, sum.lo);
}

myb_dd myb_time_diff(myb_time a, myb_time b) {
  int64_t whole;
  MYB_CHECK(!ckd_sub(&whole, a.s, b.s));
  MYB_CHECK(whole > -(INT64_C(1) << 53) && whole < (INT64_C(1) << 53));
  double err;
  const double fraction = myb_eft_two_sum(a.f, -b.f, &err); /* exact pair */
  return myb_dd_add_d(myb_dd_add_d(myb_dd_from_double((double)whole), fraction), err);
}

double myb_time_seconds(myb_time t) {
  MYB_CHECK(t.s > -(INT64_C(1) << 53) && t.s < (INT64_C(1) << 53));
  return (double)t.s + t.f;
}

int myb_time_compare(myb_time a, myb_time b) {
  if (a.s != b.s) {
    return a.s < b.s ? -1 : 1;
  }
  if (a.f != b.f) {
    return a.f < b.f ? -1 : 1;
  }
  return 0;
}

/* ---- Julian dates -------------------------------------------------------------------------- */

static const double jd_j2000 = 2451545.0;
static const double seconds_per_day = 86400.0;

myb_time myb_time_from_jd(double jd1, double jd2) {
  double err;
  const double days = myb_eft_two_sum(jd1, -jd_j2000, &err);
  myb_dd total = myb_dd_add_d((myb_dd){.hi = days, .lo = err}, jd2);
  total = myb_dd_mul_d(total, seconds_per_day);
  return myb_time_from_dd(total);
}

void myb_time_to_jd(myb_time t, double* jd1, double* jd2) {
  MYB_CHECK(jd1 != NULL && jd2 != NULL);
  const int64_t days = myb_internal_floor_div(t.s, 86400);
  const int64_t rest = t.s - days * 86400; /* 0 .. 86399 */
  *jd1 = jd_j2000 + (double)days;
  *jd2 = ((double)rest + t.f) / seconds_per_day;
}

/* ---- Calendar -------------------------------------------------------------------------------- */

int64_t myb_internal_floor_div(int64_t a, int64_t b) {
  const int64_t q = a / b;
  return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

int64_t myb_internal_days_from_civil(int64_t year, int month, int day) {
  const int64_t y = year - (month <= 2 ? 1 : 0);
  const int64_t era = myb_internal_floor_div(y, 400);
  const int64_t yoe = y - era * 400;                         /* [0, 399] */
  const int64_t mp = month > 2 ? month - 3 : month + 9;      /* March = 0 */
  const int64_t doy = (153 * mp + 2) / 5 + day - 1;          /* [0, 365] */
  const int64_t doe = yoe * 365 + yoe / 4 - yoe / 100 + doy; /* [0, 146096] */
  return era * 146097 + doe - 719468;                        /* days since 1970-01-01 */
}

void myb_internal_civil_from_days(int64_t days, int64_t* year, int* month, int* day) {
  const int64_t z = days + 719468;
  const int64_t era = myb_internal_floor_div(z, 146097);
  const int64_t doe = z - era * 146097;                                      /* [0, 146096] */
  const int64_t yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365; /* [0, 399] */
  const int64_t doy = doe - (365 * yoe + yoe / 4 - yoe / 100);               /* [0, 365] */
  const int64_t mp = (5 * doy + 2) / 153;                                    /* [0, 11] */
  *day = (int)(doy - (153 * mp + 2) / 5 + 1);
  *month = (int)(mp < 10 ? mp + 3 : mp - 9);
  *year = yoe + era * 400 + (*month <= 2 ? 1 : 0);
}

int myb_internal_days_in_month(int64_t year, int month) {
  static const int days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  const int leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
  return month == 2 && leap ? 29 : days[month - 1];
}

int myb_internal_valid_date(int64_t year, int month, int day) {
  return year >= -1000000000 && year <= 1000000000 && month >= 1 && month <= 12 && day >= 1 &&
         day <= myb_internal_days_in_month(year, month);
}

/* Days from 1970-01-01 to 2000-01-01: J2000.0 is that day at 12:00. */
static const int64_t days_to_j2000_day = 10957;

myb_status myb_calendar_to_time(const myb_calendar* calendar, myb_time* t) {
  MYB_CHECK(calendar != NULL && t != NULL);
  const myb_calendar c = *calendar;
  if (!myb_internal_valid_date(c.year, c.month, c.day) || c.hour < 0 || c.hour > 23 ||
      c.minute < 0 || c.minute > 59 || !(c.second >= 0.0 && c.second < 60.0)) {
    return MYB_E_DOMAIN;
  }
  const int64_t days = myb_internal_days_from_civil(c.year, c.month, c.day) - days_to_j2000_day;
  const int64_t whole = days * 86400 + c.hour * 3600 + c.minute * 60 - 43200;
  *t = myb_time_add((myb_time){.s = whole, .f = 0.0}, c.second);
  return MYB_OK;
}

void myb_time_to_calendar(myb_time t, myb_calendar* calendar) {
  MYB_CHECK(calendar != NULL);
  const int64_t from_midnight = t.s + 43200;
  const int64_t days = myb_internal_floor_div(from_midnight, 86400);
  const int64_t second_of_day = from_midnight - days * 86400; /* 0 .. 86399 */
  myb_calendar c = {};
  myb_internal_civil_from_days(days + days_to_j2000_day, &c.year, &c.month, &c.day);
  c.hour = (int)(second_of_day / 3600);
  c.minute = (int)(second_of_day % 3600 / 60);
  c.second = (double)(second_of_day % 60) + t.f;
  if (c.second >= 60.0) {
    /* 59 + f rounded up to 60: keep the label valid (error below 2^-47 s). */
    c.second = 0x1.dffffffffffffp5;
  }
  *calendar = c;
}
