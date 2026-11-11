/* SPDX-License-Identifier: BSD-3-Clause */
/* Space weather: a table of days in the layout of the CelesTrak files (tests/data/weather:
 * rows made for the test, with observed, interpolated and predicted days and the monthly
 * predictions that give the flux alone), what a day holds, what the models of the MSIS
 * family take at an instant, and the files that are refused. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/io.h"
#include "myb/weather.h"

#define WEATHER MYB_TEST_DATA_DIR "/weather/sw-excerpt.csv"

enum { MARCH_1 = 61100 }; /* the Modified Julian Date of 2026-03-01 */

static void test_table(void) {
  myb_space_weather* w = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_space_weather_load_file(WEATHER, NULL, &w, &diag) == MYB_OK);
  int64_t first = 0;
  int64_t last = 0;
  MT_CHECK(myb_space_weather_days(w, &first, &last) == 15);
  MT_CHECK(first == MARCH_1 && last == MARCH_1 + 61);
  myb_space_weather_day day;
  MT_CHECK(myb_space_weather_day_at(w, MARCH_1 + 3, &day) == MYB_OK);
  /* 2026-03-04: Kp comes in tenths */
  MT_CHECK(day.mjd == MARCH_1 + 3 && day.kp[0] == 1.3 && day.kp[4] == 4.7 && day.ap[4] == 39.0);
  MT_CHECK(day.ap_daily == 19.0 && day.f107_observed == 161.4 && day.f107_adjusted == 158.9);
  MT_CHECK(day.f107_observed_81[0] == 149.6 && day.f107_observed_81[1] == 151.9 &&
           day.f107_adjusted_81[0] == 147.1 && day.f107_adjusted_81[1] == 147.9);
  MT_CHECK(day.kind == MYB_SPACE_WEATHER_OBSERVED);
  MT_CHECK(myb_space_weather_day_at(w, MARCH_1 + 10, &day) == MYB_OK &&
           day.kind == MYB_SPACE_WEATHER_INTERPOLATED);
  MT_CHECK(myb_space_weather_day_at(w, MARCH_1 + 11, &day) == MYB_OK &&
           day.kind == MYB_SPACE_WEATHER_PREDICTED);
  /* a monthly prediction: the flux alone */
  MT_CHECK(myb_space_weather_day_at(w, MARCH_1 + 31, &day) == MYB_OK);
  MT_CHECK(day.kind == MYB_SPACE_WEATHER_PREDICTED && day.f107_observed == 144.0 &&
           isnan(day.ap_daily) && isnan(day.kp[3]));
  /* a day in a gap of the table, a day outside it */
  MT_CHECK(myb_space_weather_day_at(w, MARCH_1 + 20, &day) == MYB_E_RANGE);
  MT_CHECK(myb_space_weather_day_at(w, MARCH_1 - 1, &day) == MYB_E_RANGE);

  /* 2026-03-05 at 13:00 UTC: the fifth span of three hours */
  myb_space_weather_inputs in;
  MT_CHECK(myb_space_weather_inputs_at(w, MARCH_1 + 4, 13.0 * 3600.0, 0, &in) == MYB_OK);
  MT_CHECK(in.f107 == 161.4);         /* of the day before */
  MT_CHECK(in.f107_average == 150.0); /* centred on the day */
  MT_CHECK(in.ap == 6.0 && in.kp == 1.0 && in.has_history);
  /* ap of the span and of the three before it; then the eight from 12 to 33 hours before
   * (the first span of the day back to the second of the day before), and the eight before
   * those */
  MT_CHECK(in.ap_history[0] == 6.0 && in.ap_history[1] == 4.0 && in.ap_history[2] == 5.0 &&
           in.ap_history[3] == 6.0 && in.ap_history[4] == 7.0);
  const double late = (9.0 + 12.0 + 18.0 + 32.0 + 39.0 + 27.0 + 12.0 + 6.0) / 8.0;
  const double early = (5.0 + 4.0 + 5.0 + 4.0 + 3.0 + 2.0 + 2.0 + 3.0) / 8.0;
  MT_CHECK(fabs(in.ap_history[5] - late) < 1e-15 && fabs(in.ap_history[6] - early) < 1e-15);
  /* the flux adjusted to 1 au */
  MT_CHECK(myb_space_weather_inputs_at(w, MARCH_1 + 4, 13.0 * 3600.0, 1, &in) == MYB_OK);
  MT_CHECK(in.f107 == 158.9 && in.f107_average == 147.6);
  /* the last instant of a day, and a leap second after it: the last span */
  MT_CHECK(myb_space_weather_inputs_at(w, MARCH_1 + 4, 86400.5, 0, &in) == MYB_OK);
  MT_CHECK(in.kp == 2.0 && in.ap_history[1] == 7.0);
  /* early in the table the 57 hours are not there: the daily value seven times */
  MT_CHECK(myb_space_weather_inputs_at(w, MARCH_1 + 1, 3600.0, 0, &in) == MYB_OK);
  MT_CHECK(!in.has_history && in.ap == 8.0 && in.ap_history[1] == 8.0 &&
           in.ap_history[6] == 8.0 && in.kp == 2.7);
  /* no day before; a day without Ap; a gap; seconds that are none */
  MT_CHECK(myb_space_weather_inputs_at(w, MARCH_1, 0.0, 0, &in) == MYB_E_RANGE);
  MT_CHECK(myb_space_weather_inputs_at(w, MARCH_1 + 31, 0.0, 0, &in) == MYB_E_RANGE);
  MT_CHECK(myb_space_weather_inputs_at(w, MARCH_1 + 13, 0.0, 0, &in) == MYB_E_RANGE);
  MT_CHECK(myb_space_weather_inputs_at(w, MARCH_1 + 4, -1.0, 0, &in) == MYB_E_DOMAIN);
  MT_CHECK(myb_space_weather_inputs_at(w, MARCH_1 + 4, 90000.0, 0, &in) == MYB_E_DOMAIN);
  printf("space weather: %d days of a table by the names of its columns; at 13:00 of "
         "2026-03-05 F10.7 = %.1f (81 days: %.1f), Ap = %.0f, the history of ap %.0f %.0f "
         "%.0f %.0f %.3f %.3f\n",
         15, 161.4, 150.0, 6.0, 4.0, 5.0, 6.0, 7.0, late, early);
  myb_space_weather_destroy(w);
  myb_space_weather_destroy(NULL);
}

static myb_status load(const char* text, myb_space_weather** w, myb_diag* diag) {
  return myb_space_weather_load(text, strlen(text), NULL, w, diag);
}

static void test_layouts(void) {
  myb_space_weather* w = NULL;
  myb_diag diag = {0};
  /* the columns in any order, others among them, comments, blank lines, CR LF */
  const char* mixed = "# a table\r\nF10.7_OBS, other ,DATE,AP_AVG\r\n\r\n"
                      "70.5,x,1999-12-31,4\r\n71.5,,2000-01-01,\r\n";
  MT_CHECK(load(mixed, &w, &diag) == MYB_OK);
  int64_t first = 0;
  MT_CHECK(myb_space_weather_days(w, &first, NULL) == 2 && first == 51543);
  myb_space_weather_day day;
  MT_CHECK(myb_space_weather_day_at(w, 51544, &day) == MYB_OK);
  MT_CHECK(day.f107_observed == 71.5 && isnan(day.ap_daily) && isnan(day.f107_adjusted));
  myb_space_weather_destroy(w);
  /* refused, with the line */
  MT_CHECK(load("F10.7_OBS,AP_AVG\n70,4\n", &w, &diag) == MYB_E_FORMAT && diag.line == 1);
  MT_CHECK(load("DATE,F10.7_OBS\n2000-01-01,70\n2000-1-2,71\n", &w, &diag) == MYB_E_FORMAT &&
           diag.line == 3);
  MT_CHECK(load("DATE,F10.7_OBS\n2000-02-30,70\n", &w, &diag) == MYB_E_FORMAT);
  MT_CHECK(load("DATE,F10.7_OBS\n2000-01-02,70\n2000-01-01,71\n", &w, &diag) == MYB_E_FORMAT &&
           diag.line == 3);
  MT_CHECK(load("DATE,F10.7_OBS\n2000-01-01,seventy\n", &w, &diag) == MYB_E_FORMAT);
  MT_CHECK(load("DATE,F10.7_OBS\n,70\n", &w, &diag) == MYB_E_FORMAT);
  MT_CHECK(load("DATE,F10.7_OBS\n", &w, &diag) == MYB_E_FORMAT);
  MT_CHECK(load("", &w, &diag) == MYB_E_FORMAT);
  MT_CHECK(myb_space_weather_load_file(MYB_TEST_DATA_DIR "/weather/none.csv", NULL, &w, &diag) ==
           MYB_E_IO);
}

/* What a model of the upper atmosphere is asked with: the day, the time, the place and the
 * indices of an instant on TT at a place over 30 E. */
static void test_point(void) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/time-scales.tk", &diag) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/tdb-tt-fb1990.tk", &diag) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/utc-tai.tk", &diag) == MYB_OK);
  myb_timesys* sys = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_scale tt;
  MT_CHECK(myb_timesys_scale(sys, "TT", &tt) == MYB_OK);
  myb_space_weather* w = NULL;
  MT_CHECK(myb_space_weather_load_file(WEATHER, NULL, &w, &diag) == MYB_OK);
  const myb_ellipsoid wgs84 = {6378137.0, 1.0 / 298.257223563};
  const myb_calendar utc = {.year = 2026, .month = 3, .day = 5, .hour = 13};
  myb_time t;
  MT_CHECK(myb_timesys_from_utc(sys, &utc, tt, &t) == MYB_OK);
  const double degree = 3.141592653589793 / 180.0;
  double r[3];
  MT_CHECK(myb_geodetic_to_cartesian(&wgs84, 30.0 * degree, 51.6 * degree, 420.0e3, r) ==
           MYB_OK);
  myb_thermosphere_setting setting = {.weather = w, .timesys = sys, .scale = tt,
                                      .ellipsoid = wgs84};
  myb_thermosphere_point p;
  MT_CHECK(myb_thermosphere_point_at(&setting, t, r, &p) == MYB_OK);
  MT_CHECK(p.mjd == MARCH_1 + 4 && p.day_of_year == 64 && fabs(p.seconds - 46800.0) < 1e-6);
  MT_CHECK(fabs(p.height - 420.0e3) < 1e-6 && fabs(p.latitude - 51.6 * degree) < 1e-12 &&
           fabs(p.longitude - 30.0 * degree) < 1e-12);
  MT_CHECK(fabs(p.local_solar_time - 15.0) < 1e-9);
  MT_CHECK(p.indices.f107 == 161.4 && p.indices.f107_average == 150.0 && p.indices.ap == 6.0 &&
           p.indices.has_history && p.indices.ap_history[1] == 4.0);
  /* the flux adjusted to 1 au; no table: the three numbers, the daily Ap for the history */
  setting.adjusted = 1;
  MT_CHECK(myb_thermosphere_point_at(&setting, t, r, &p) == MYB_OK);
  MT_CHECK(p.indices.f107 == 158.9 && p.indices.f107_average == 147.6);
  setting = (myb_thermosphere_setting){.f107 = 120.0, .f107_average = 110.0, .ap = 9.0,
                                       .timesys = sys, .scale = tt, .ellipsoid = wgs84};
  MT_CHECK(myb_thermosphere_point_at(&setting, myb_time_add(t, 400.0 * 86400.0), r, &p) ==
           MYB_OK);
  MT_CHECK(p.indices.f107 == 120.0 && p.indices.ap == 9.0 && !p.indices.has_history &&
           p.indices.ap_history[6] == 9.0 && p.day_of_year == 99);
  /* outside the table */
  setting.weather = w;
  MT_CHECK(myb_thermosphere_point_at(&setting, myb_time_add(t, 400.0 * 86400.0), r, &p) ==
           MYB_E_RANGE);
  myb_space_weather_destroy(w);
  myb_timesys_destroy(sys);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_table();
  test_layouts();
  test_point();
  return mt_finish("test_weather");
}
