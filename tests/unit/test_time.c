/* SPDX-License-Identifier: BSD-3-Clause */
/* Instants, calendar, time scale systems and UTC: exact properties and known values.
 * Comparison with IAU SOFA (through ERFA) is in tests/erfa/test_erfa_time.c. */
#include <math.h>
#include <string.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/time.h"

static double seconds_between(myb_time a, myb_time b) {
  return myb_dd_to_double(myb_time_diff(a, b));
}

/* ---- Instants ------------------------------------------------------------------------------ */

static void test_instants(void) {
  const myb_time t = myb_time_from_seconds(-1.25);
  MT_CHECK(t.s == -2 && t.f == 0.75);
  const myb_time u = myb_time_add(t, 1e9 + 0.5);
  MT_CHECK(u.s == 999999999 && u.f == 0.25);
  MT_CHECK(myb_time_compare(t, u) < 0 && myb_time_compare(u, t) > 0 && myb_time_compare(u, u) == 0);
  MT_CHECK_SAME(seconds_between(u, t), 1e9 + 0.5);

  /* Resolution far beyond a double: 10^10 s plus 10^-15 s survives. */
  const myb_time big = myb_time_add(myb_time_from_seconds(1e10), 1e-15);
  MT_CHECK(big.s == 10000000000 && big.f == 1e-15);
  MT_CHECK_SAME(seconds_between(big, myb_time_from_seconds(1e10)), 1e-15);

  /* A fraction that rounds to 1 becomes the next whole second. */
  const myb_time edge = myb_time_add((myb_time){.s = 5, .f = 0.0}, 1.0 - 0x1p-60);
  MT_CHECK(edge.s == 6 && edge.f == 0.0);

  /* Double-word input keeps the low part. */
  const myb_time dd = myb_time_from_dd((myb_dd){.hi = 3e9, .lo = 1e-7});
  MT_CHECK(dd.s == 3000000000 && fabs(dd.f - 1e-7) < 1e-22);

  mt_rng rng = {.state = 501};
  for (int i = 0; i < 100000; ++i) {
    const myb_time a = myb_time_from_seconds(mt_double(&rng, -20, 33));
    const double step = mt_double(&rng, -30, 30);
    const myb_time b = myb_time_add(a, step);
    MT_CHECK(b.f >= 0.0 && b.f < 1.0);
    MT_CHECK(fabs(seconds_between(b, a) - step) <= 0x1p-52 * (1.0 + fabs(step)));
  }
}

/* ---- Julian dates and calendar ------------------------------------------------------------- */

static myb_calendar cal(int64_t y, int mo, int d, int h, int mi, double s) {
  return (myb_calendar){.year = y, .month = mo, .day = d, .hour = h, .minute = mi, .second = s};
}

static myb_time from_cal(myb_calendar c) {
  myb_time t = {};
  MT_CHECK(myb_calendar_to_time(&c, &t) == MYB_OK);
  return t;
}

static void test_calendar(void) {
  MT_CHECK(myb_time_compare(from_cal(cal(2000, 1, 1, 12, 0, 0.0)), (myb_time){}) == 0);
  const myb_time unix_epoch = from_cal(cal(1970, 1, 1, 0, 0, 0.0));
  MT_CHECK(unix_epoch.s == -946728000 && unix_epoch.f == 0.0);
  double jd1 = 0.0;
  double jd2 = 0.0;
  myb_time_to_jd(unix_epoch, &jd1, &jd2);
  MT_CHECK(jd1 + jd2 == 2440587.5);
  MT_CHECK(myb_time_compare(myb_time_from_jd(2440587.5, 0.0), unix_epoch) == 0);
  MT_CHECK(myb_time_compare(myb_time_from_jd(2400000.5, 51544.5), (myb_time){}) == 0);
  /* Gregorian reform date and year 0 (1 BC, a leap year). */
  MT_CHECK(myb_time_compare(from_cal(cal(1582, 10, 15, 12, 0, 0.0)),
                            myb_time_from_jd(2299161.0, 0.0)) == 0);
  myb_time t = {};
  myb_calendar c = cal(0, 2, 29, 0, 0, 0.0);
  MT_CHECK(myb_calendar_to_time(&c, &t) == MYB_OK);
  c = cal(1900, 2, 29, 0, 0, 0.0);
  MT_CHECK(myb_calendar_to_time(&c, &t) == MYB_E_DOMAIN);
  c = cal(2023, 4, 31, 0, 0, 0.0);
  MT_CHECK(myb_calendar_to_time(&c, &t) == MYB_E_DOMAIN);
  c = cal(2023, 1, 1, 0, 0, 60.0);
  MT_CHECK(myb_calendar_to_time(&c, &t) == MYB_E_DOMAIN);
  c = cal(2023, 13, 1, 0, 0, 0.0);
  MT_CHECK(myb_calendar_to_time(&c, &t) == MYB_E_DOMAIN);
  c = cal(2023, 1, 1, 24, 0, 0.0);
  MT_CHECK(myb_calendar_to_time(&c, &t) == MYB_E_DOMAIN);

  /* Round trips over ±10^6 years. */
  mt_rng rng = {.state = 502};
  for (int i = 0; i < 200000; ++i) {
    const myb_time a = {
        .s = (int64_t)(mt_next(&rng) % UINT64_C(63113904000000)) - INT64_C(31556952000000),
        .f = (double)(mt_next(&rng) >> 11) * 0x1p-53};
    myb_calendar k;
    myb_time_to_calendar(a, &k);
    myb_time b = {};
    MT_CHECK(myb_calendar_to_time(&k, &b) == MYB_OK);
    MT_CHECK(fabs(seconds_between(a, b)) <= 0x1p-47);
    myb_time_to_jd(a, &jd1, &jd2);
    MT_CHECK(jd2 >= 0.0 && jd2 < 1.0);
    MT_CHECK(fabs(seconds_between(myb_time_from_jd(jd1, jd2), a)) <= 2e-11);
  }
}

/* ---- Time scale systems ---------------------------------------------------------------------- */

static myb_timesys* load_system(void) {
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const char* files[] = {MYB_DATA_DIR "/time/time-scales.tk", MYB_DATA_DIR "/time/tdb-tt-fb1990.tk",
                         MYB_DATA_DIR "/time/utc-tai.tk"};
  for (size_t i = 0; i < sizeof files / sizeof files[0]; ++i) {
    myb_diag diag;
    const myb_status status = myb_kpool_load_file(pool, files[i], &diag);
    MT_CHECK(status == MYB_OK);
    if (status != MYB_OK) {
      fprintf(stderr, "%s: %s\n", files[i], diag.message);
    }
  }
  myb_timesys* system = NULL;
  myb_diag diag;
  const myb_status status = myb_timesys_create(pool, NULL, &system, &diag);
  MT_CHECK(status == MYB_OK);
  if (status != MYB_OK) {
    fprintf(stderr, "timesys: %s\n", diag.message);
  }
  myb_kpool_destroy(pool);
  return system;
}

static myb_scale scale(const myb_timesys* s, const char* name) {
  myb_scale sc = {.index = 0};
  MT_CHECK(myb_timesys_scale(s, name, &sc) == MYB_OK);
  return sc;
}

static myb_time convert(const myb_timesys* s, myb_time t, const char* from, const char* to) {
  myb_time out = {};
  MT_CHECK(myb_timesys_convert(s, t, scale(s, from), scale(s, to), &out) == MYB_OK);
  return out;
}

static void test_scales(const myb_timesys* s) {
  myb_scale sc;
  MT_CHECK(myb_timesys_scale(s, "UTC", &sc) == MYB_E_NOT_FOUND);
  MT_CHECK(strcmp(myb_timesys_scale_name(s, scale(s, "TCB")), "TCB") == 0);

  const myb_time t = myb_time_from_seconds(8.5e8 + 0.125);
  /* 32.184 is not a binary number: the stored instant differs by its rounding only. */
  MT_CHECK(fabs(seconds_between(convert(s, t, "TAI", "TT"), t) - 32.184) < 1e-15);
  MT_CHECK_SAME(seconds_between(convert(s, t, "TAI", "GPS"), t), -19.0);
  MT_CHECK_SAME(seconds_between(convert(s, t, "GPS", "BDT"), t), -14.0);

  /* TCG = TT at T0 = 1977-01-01T00:00:32.184; TCB - TDB = -TDB0 at T0 (TCB). */
  const myb_time t0 = myb_time_from_jd(2443144.5, 0.0003725);
  MT_CHECK(fabs(seconds_between(convert(s, t0, "TT", "TCG"), t0)) < 1e-15);
  MT_CHECK(fabs(seconds_between(convert(s, t0, "TCB", "TDB"), t0) - (-6.55e-5)) < 1e-15);
  /* Rate: TCG - TT = L_G / (1 - L_G) * (TT - T0). */
  const myb_time later = myb_time_add(t0, 1e9);
  const double lg = 6.969290134e-10;
  MT_CHECK(fabs(seconds_between(convert(s, later, "TT", "TCG"), later) - lg / (1.0 - lg) * 1e9) <
           1e-12);

  /* Every pair converts there and back to within 10^-15 s; TT-TDB stays below 2 ms. */
  const char* names[] = {"TAI", "TT", "TCG", "TDB", "TCB", "GPS", "GST", "BDT"};
  mt_rng rng = {.state = 503};
  double worst = 0.0;
  for (int i = 0; i < 3000; ++i) {
    const myb_time a = myb_time_from_seconds(mt_double(&rng, 0, 34) * ((i & 1) ? 1.0 : -1.0));
    for (size_t p = 0; p < 8; ++p) {
      for (size_t q = 0; q < 8; ++q) {
        const myb_time b = convert(s, a, names[p], names[q]);
        const double e = fabs(seconds_between(convert(s, b, names[q], names[p]), a));
        worst = e > worst ? e : worst;
      }
    }
    MT_CHECK(fabs(seconds_between(convert(s, a, "TT", "TDB"), a)) < 2e-3);
  }
  MT_CHECK(worst < 1e-15);
  printf("scale round trips: worst %.3g s\n", worst);
}

static myb_time from_utc(const myb_timesys* s, myb_calendar c, const char* to,
                         myb_status expected) {
  myb_time t = {};
  const myb_status status = myb_timesys_from_utc(s, &c, scale(s, to), &t);
  mtest.checks++;
  if (status != expected) {
    mtest.failures++;
    fprintf(stderr, "from_utc %lld-%02d-%02dT%02d:%02d:%g: %s, expected %s\n", (long long)c.year,
            c.month, c.day, c.hour, c.minute, c.second, myb_status_name(status),
            myb_status_name(expected));
  }
  return t;
}

static void test_utc(const myb_timesys* s) {
  /* The leap second at the end of 2016: TAI - UTC goes from 36 s to 37 s. */
  const myb_time leap = from_utc(s, cal(2016, 12, 31, 23, 59, 60.5), "TAI", MYB_OK);
  MT_CHECK(myb_time_compare(leap, from_cal(cal(2017, 1, 1, 0, 0, 36.5))) == 0);
  const myb_time after = from_utc(s, cal(2017, 1, 1, 0, 0, 0.0), "TAI", MYB_OK);
  MT_CHECK(myb_time_compare(after, from_cal(cal(2017, 1, 1, 0, 0, 37.0))) == 0);
  myb_calendar back;
  MT_CHECK(myb_timesys_to_utc(s, leap, scale(s, "TAI"), &back) == MYB_OK);
  MT_CHECK(back.year == 2016 && back.month == 12 && back.day == 31 && back.hour == 23 &&
           back.minute == 59 && back.second == 60.5);

  /* Invalid labels and coverage. */
  from_utc(s, cal(2016, 12, 30, 23, 59, 60.0), "TAI", MYB_E_DOMAIN); /* no leap that day */
  from_utc(s, cal(2016, 12, 31, 23, 59, 61.0), "TAI", MYB_E_DOMAIN);
  from_utc(s, cal(2016, 12, 31, 23, 58, 60.0), "TAI", MYB_E_DOMAIN);
  from_utc(s, cal(1959, 12, 31, 0, 0, 0.0), "TAI", MYB_E_RANGE);
  from_utc(s, cal(2027, 6, 27, 23, 59, 59.0), "TAI", MYB_OK);
  from_utc(s, cal(2027, 6, 28, 0, 0, 0.0), "TAI", MYB_E_RANGE); /* expiry of the table */
  /* 1961-07-31: TAI - UTC stepped by -0.05 s at its end, the last minute has 59.95 s. */
  from_utc(s, cal(1961, 7, 31, 23, 59, 59.94), "TAI", MYB_OK);
  from_utc(s, cal(1961, 7, 31, 23, 59, 59.96), "TAI", MYB_E_DOMAIN);

  /* Pre-1972: 1965-06-01 0h UTC = MJD 38912 (38761 + 151), row of 1965-03-01:
   * TAI - UTC = 3.6401300 + (38912 - 38761) * 0.0012960. */
  const myb_time old = from_utc(s, cal(1965, 6, 1, 0, 0, 0.0), "TAI", MYB_OK);
  const double expected = 3.6401300 + (38912.0 - 38761.0) * 0.0012960;
  MT_CHECK(fabs(seconds_between(old, from_cal(cal(1965, 6, 1, 0, 0, 0.0))) - expected) < 1e-12);

  /* Round trips TAI -> UTC -> TAI over 1960-2027. */
  mt_rng rng = {.state = 504};
  const myb_time first = from_cal(cal(1960, 1, 1, 0, 0, 2.0));
  const double span = seconds_between(from_cal(cal(2027, 6, 27, 0, 0, 0.0)), first);
  double worst = 0.0;
  for (int i = 0; i < 200000; ++i) {
    const myb_time tai = myb_time_add(first, span * ((double)(mt_next(&rng) >> 11) * 0x1p-53));
    myb_calendar label;
    MT_CHECK(myb_timesys_to_utc(s, tai, scale(s, "TAI"), &label) == MYB_OK);
    myb_time again = {};
    MT_CHECK(myb_timesys_from_utc(s, &label, scale(s, "TAI"), &again) == MYB_OK);
    const double e = fabs(seconds_between(again, tai));
    worst = e > worst ? e : worst;
  }
  MT_CHECK(worst < 1e-12);
  printf("UTC round trips: worst %.3g s\n", worst);
  /* Every leap second label round-trips exactly. */
  myb_calendar k;
  MT_CHECK(myb_timesys_to_utc(s, from_utc(s, cal(1972, 6, 30, 23, 59, 60.25), "TT", MYB_OK),
                              scale(s, "TT"), &k) == MYB_OK);
  MT_CHECK(k.minute == 59 && k.second == 60.25);
}

static void test_bad_definitions(void) {
  struct {
    const char* text;
    myb_status status;
    const char* fragment;
  } cases[] = {
      {"\\begindata\nMYB_TIME_SCALES = ( 'A' 'B' )\n", MYB_E_NOT_FOUND, "MYB_TIME_B_FROM"},
      {"\\begindata\nMYB_TIME_SCALES = ( 'A' 'B' )\nMYB_TIME_B_FROM = 'C'\n", MYB_E_FORMAT,
       "not a scale listed before"},
      {"\\begindata\nMYB_TIME_SCALES = ( 'A' 'B' )\nMYB_TIME_B_FROM = 'A'\n"
       "MYB_TIME_B_KIND = 'MAGIC'\n",
       MYB_E_FORMAT, "unknown kind"},
      {"\\begindata\nMYB_TIME_SCALES = ( 'A' 'A' )\n", MYB_E_FORMAT, "not a new name"},
      {"\\begindata\nMYB_TIME_SCALES = ( 'A' 'b' )\n", MYB_E_FORMAT, "not a new name"},
  };

  for (size_t i = 0; i < sizeof cases / sizeof cases[0]; ++i) {
    myb_kpool* pool = NULL;
    MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
    MT_CHECK(myb_kpool_load_text(pool, cases[i].text, strlen(cases[i].text), NULL) == MYB_OK);
    myb_timesys* system = NULL;
    myb_diag diag;
    const myb_status status = myb_timesys_create(pool, NULL, &system, &diag);
    mtest.checks++;
    if (status != cases[i].status || strstr(diag.message, cases[i].fragment) == NULL) {
      mtest.failures++;
      fprintf(stderr, "bad definition %zu: %s '%s'\n", i, myb_status_name(status), diag.message);
    }
    MT_CHECK(system == NULL);
    myb_kpool_destroy(pool);
  }
}

int main(void) {
  test_instants();
  test_calendar();
  myb_timesys* system = load_system();
  if (system != NULL) {
    test_scales(system);
    test_utc(system);
  }
  myb_timesys_destroy(system);
  test_bad_definitions();
  return mt_finish("test_time");
}
