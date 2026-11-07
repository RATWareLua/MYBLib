/* SPDX-License-Identifier: BSD-3-Clause */
/* Time scales against IAU SOFA, through its BSD-licensed copy ERFA (tests/third_party/erfa,
 * used only here). ERFA works with two-part Julian dates, whose resolution limits the
 * agreement to ~10^-11 s; the requirement of the plan is 1 ns. Maximum differences are
 * printed. */
#include <math.h>

#include "core/fpguard.h"
#include "erfa.h"
#include "mtest.h"
#include "myb/time.h"

static myb_timesys* system_ = NULL;

static myb_scale scale(const char* name) {
  myb_scale s = {.index = 0};
  MT_CHECK(myb_timesys_scale(system_, name, &s) == MYB_OK);
  return s;
}

static double seconds(myb_time a, myb_time b) {
  return myb_dd_to_double(myb_time_diff(a, b));
}

static void report(const char* what, double worst, double tolerance, long count) {
  printf("%-28s %7ld cases, max |difference| %.3g s (tolerance %.0e s)\n", what, count, worst,
         tolerance);
  MT_CHECK(worst <= tolerance);
}

static void test_calendar(void) {
  mt_rng rng = {.state = 601};
  double worst = 0.0;
  for (int i = 0; i < 200000; ++i) {
    const int year = mt_int(&rng, -4700, 3000);
    const int month = mt_int(&rng, 1, 12);
    const int day = mt_int(&rng, 1, 28);
    double d1 = 0.0;
    double d2 = 0.0;
    MT_CHECK(eraCal2jd(year, month, day, &d1, &d2) == 0);
    const myb_calendar c = {.year = year, .month = month, .day = day};
    myb_time t = {};
    MT_CHECK(myb_calendar_to_time(&c, &t) == MYB_OK);
    const double e = fabs(seconds(t, myb_time_from_jd(d1, d2)));
    worst = e > worst ? e : worst;
  }
  report("calendar (eraCal2jd)", worst, 0.0, 200000);
}

/* TDB - TT itself: our series S(TDB) against eraDtdb(TDB) without topocentric terms. */
static void test_tdb_tt(void) {
  mt_rng rng = {.state = 602};
  double worst = 0.0;
  const long count = 20000;
  for (long i = 0; i < count; ++i) {
    const double jd2 = (mt_int(&rng, 0, 1) != 0 ? 1.0 : -1.0) * 365.25 * 450.0 *
                       ((double)(mt_next(&rng) >> 11) * 0x1p-53); /* 1550 .. 2450 */
    const myb_time tdb = myb_time_from_jd(2451545.0, jd2);
    myb_time tt = {};
    MT_CHECK(myb_timesys_convert(system_, tdb, scale("TDB"), scale("TT"), &tt) == MYB_OK);
    const double ours = seconds(tdb, tt);
    const double theirs = eraDtdb(2451545.0, jd2, 0.0, 0.0, 0.0, 0.0);
    const double e = fabs(ours - theirs);
    worst = e > worst ? e : worst;
  }
  /* ERFA calls the platform's sin (not correctly rounded, differs between C runtimes) and sums
   * in another order: 1e-14 s leaves room for that and is still 1e5 times below 1 ns. */
  report("TDB-TT series (eraDtdb)", worst, 1e-14, count);
}

typedef int (*erfa_conversion)(double, double, double*, double*);

/* ERFA keeps the input split in its output: with (2451545.0, +-1e5 days) the result's second
 * part would be rounded to ulp(1e5 days) ~ 1e-6 s. The "date & time" split of the SOFA
 * documentation (whole days + fraction of a day) gives ERFA its full ~1e-11 s resolution. */
static void test_rate(const char* what, const char* from, const char* to, erfa_conversion f,
                      uint64_t seed) {
  mt_rng rng = {.state = seed};
  double worst = 0.0;
  const long count = 50000;
  for (long i = 0; i < count; ++i) {
    const double jd1 = 2451545.0 + (double)mt_int(&rng, -110000, 110000);
    const double jd2 = (double)(mt_next(&rng) >> 11) * 0x1p-53;
    const myb_time t = myb_time_from_jd(jd1, jd2);
    myb_time ours = {};
    MT_CHECK(myb_timesys_convert(system_, t, scale(from), scale(to), &ours) == MYB_OK);
    double r1 = 0.0;
    double r2 = 0.0;
    MT_CHECK(f(jd1, jd2, &r1, &r2) == 0);
    const double e = fabs(seconds(ours, myb_time_from_jd(r1, r2)));
    worst = e > worst ? e : worst;
  }
  report(what, worst, 1e-10, count);
}

static void test_utc(void) {
  mt_rng rng = {.state = 605};
  double worst = 0.0;
  long count = 0;
  for (int i = 0; i < 100000; ++i) {
    const int year = mt_int(&rng, 1960, 2026);
    const int month = mt_int(&rng, 1, 12);
    const int day = mt_int(&rng, 1, 28);
    const int hour = mt_int(&rng, 0, 23);
    const int minute = mt_int(&rng, 0, 59);
    const double second = 60.0 * ((double)(mt_next(&rng) >> 11) * 0x1p-53);
    double u1 = 0.0;
    double u2 = 0.0;
    double a1 = 0.0;
    double a2 = 0.0;
    MT_CHECK(eraDtf2d("UTC", year, month, day, hour, minute, second, &u1, &u2) >= 0);
    MT_CHECK(eraUtctai(u1, u2, &a1, &a2) >= 0);
    const myb_calendar c = {
        .year = year, .month = month, .day = day, .hour = hour, .minute = minute, .second = second};
    myb_time tai = {};
    MT_CHECK(myb_timesys_from_utc(system_, &c, scale("TAI"), &tai) == MYB_OK);
    const double e = fabs(seconds(tai, myb_time_from_jd(a1, a2)));
    worst = e > worst ? e : worst;
    ++count;
  }
  /* Inside every positive leap second and on the step days before 1972. */
  static const int leap_days[][3] = {{1961, 7, 31}, {1961, 12, 31}, {1963, 10, 31}, {1964, 3, 31},
                                     {1965, 6, 30}, {1972, 6, 30},  {1972, 12, 31}, {1979, 12, 31},
                                     {1985, 6, 30}, {1998, 12, 31}, {2005, 12, 31}, {2008, 12, 31},
                                     {2012, 6, 30}, {2015, 6, 30},  {2016, 12, 31}};
  for (size_t i = 0; i < sizeof leap_days / sizeof leap_days[0]; ++i) {
    for (int k = 0; k < 100; ++k) {
      const double second = 59.0 + 0.0094 * k; /* up to 59.93: valid on every listed day */
      double u1 = 0.0;
      double u2 = 0.0;
      double a1 = 0.0;
      double a2 = 0.0;
      const int* d = leap_days[i];
      MT_CHECK(eraDtf2d("UTC", d[0], d[1], d[2], 23, 59, second, &u1, &u2) >= 0);
      MT_CHECK(eraUtctai(u1, u2, &a1, &a2) >= 0);
      const myb_calendar c = {
          .year = d[0], .month = d[1], .day = d[2], .hour = 23, .minute = 59, .second = second};
      myb_time tai = {};
      MT_CHECK(myb_timesys_from_utc(system_, &c, scale("TAI"), &tai) == MYB_OK);
      const double e = fabs(seconds(tai, myb_time_from_jd(a1, a2)));
      worst = e > worst ? e : worst;
      ++count;
    }
  }
  report("UTC -> TAI (eraUtctai)", worst, 1e-9, count);
}

static myb_timesys* load(void) {
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/time-scales.tk", NULL) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/tdb-tt-fb1990.tk", NULL) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/utc-tai.tk", NULL) == MYB_OK);
  myb_timesys* s = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &s, NULL) == MYB_OK);
  myb_kpool_destroy(pool);
  return s;
}

int main(void) {
  system_ = load();
  if (system_ == NULL) {
    return mt_finish("test_erfa_time");
  }
  test_calendar();
  test_tdb_tt();
  test_rate("TT -> TCG (eraTttcg)", "TT", "TCG", eraTttcg, 603);
  test_rate("TCG -> TT (eraTcgtt)", "TCG", "TT", eraTcgtt, 604);
  test_rate("TDB -> TCB (eraTdbtcb)", "TDB", "TCB", eraTdbtcb, 606);
  test_rate("TCB -> TDB (eraTcbtdb)", "TCB", "TDB", eraTcbtdb, 607);
  test_rate("TAI -> TT (eraTaitt)", "TAI", "TT", eraTaitt, 608);
  test_utc();
  myb_timesys_destroy(system_);
  return mt_finish("test_erfa_time");
}
