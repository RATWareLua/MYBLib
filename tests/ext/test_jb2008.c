/* SPDX-License-Identifier: BSD-3-Clause */
/* ext/jb2008: the clean implementation of JB2008.
 *  - the static model against the tables of Jacchia (SAO Special Report 313, 1970, Table 5,
 *    exospheric temperature 1000 K): temperatures, number densities and densities;
 *  - the readers of the files of indices (files made here in their published formats) and the
 *    lags of the indices;
 *  - properties: nothing below 90 km, the temperature continuous through the joins of the
 *    corrections, the storm term added to the exospheric temperature, the density falling
 *    with the height, the atmosphere of a force model against the model called by hand;
 *  - if the file data/download/jb2008_set_outputs.txt is there (made locally by running the
 *    program of SET as a black box; never distributed), the outputs against it. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "jb2008_internal.h"
#include "mtest.h"
#include "myb/cmath.h"
#include "myb/io.h"
#include "myb_jb2008.h"

static const double pi = 3.141592653589793;

/* Jacchia 1970, Table 5, T_inf = 1000 K: z, T, log10 n (cm^-3) of N2, O2, O, Ar, He, H (0: not
 * printed), log10 rho (g cm^-3). */
static const double table5[][9] = {
    {90, 183.0, 0, 0, 0, 0, 0, 0, -8.461},
    {100, 195.9, 0, 0, 0, 0, 0, 0, -9.262},
    {110, 252.9, 12.1973, 11.4800, 11.3812, 10.1487, 7.6948, 0, -10.015},
    {120, 358.1, 11.5831, 10.8000, 10.9656, 9.3373, 7.5350, 0, -10.624},
    {125, 421.8, 11.3344, 10.5260, 10.7931, 9.0129, 7.4655, 0, -10.868},
    {150, 690.1, 10.5033, 9.6070, 10.2267, 7.9188, 7.2448, 0, -11.667},
    {200, 889.2, 9.5513, 8.5352, 9.6358, 6.6081, 7.0562, 0, -12.524},
    {300, 978.9, 8.1057, 6.8899, 8.7922, 4.5643, 6.8298, 0, -13.641},
    {400, 994.3, 6.7947, 5.3934, 8.0406, 2.6977, 6.6392, 0, -14.488},
    {500, 997.9, 5.5371, 3.9571, 7.3217, 0.9050, 6.4588, 4.4359, -15.227},
    {600, 999.0, 4.3192, 2.5660, 6.6259, 0, 6.2845, 4.3916, -15.899},
    {800, 999.7, 1.9884, 0, 5.2946, 0, 5.9513, 4.3074, -16.950},
    {1000, 999.9, 0, 0, 4.0362, 0, 5.6365, 4.2281, -17.496},
    {1500, 1000.0, 0, 0, 1.1712, 0, 4.9198, 4.0475, -18.243},
    {2500, 1000.0, 0, 0, 0, 0, 3.7291, 3.7477, -19.347}};

static void test_static(void) {
  double worst_t = 0.0;
  double worst_n = 0.0; /* N2, O2, O, Ar, H */
  double worst_he = 0.0;
  double worst_rho = 0.0;
  for (size_t k = 0; k < sizeof table5 / sizeof table5[0]; ++k) {
    const double* row = table5[k];
    double temperature = 0.0;
    double n[6];
    const double rho = myb_internal_jb2008_static(1000.0, row[0], &temperature, n);
    worst_t = fmax(worst_t, fabs(temperature - row[1]));
    for (int s = 0; s < 6; ++s) {
      if (row[2 + s] == 0.0) {
        continue;
      }
      const double d = fabs(log10(n[s]) - row[2 + s]);
      if (s == 4) {
        worst_he = fmax(worst_he, d);
      } else {
        worst_n = fmax(worst_n, d);
      }
    }
    worst_rho = fmax(worst_rho, fabs(log10(rho) - row[8]));
  }
  printf("the static model against Jacchia 1970, Table 5 (1000 K, 90 to 2500 km): temperature "
         "%.3f K, log n of N2, O2, O, Ar, H %.1e, of He %.1e, log rho %.1e\n",
         worst_t, worst_n, worst_he, worst_rho);
  /* the table's own rounding: 0.05 K, 5e-5 in log n, 5e-4 in log rho; helium of the table is
   * 0.0012 above its own equations with the printed abundance at sea level */
  MT_CHECK(worst_t < 0.06 && worst_n < 3e-4 && worst_he < 1.5e-3 && worst_rho < 1.2e-3);
}

static const char solar_text[] =
    "# YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c  Ssrc\n"
    "  2004 365   2453370.0 101.0 111.0 102.0 112.0 103.0 113.0 104.0 114.0  1B11\n"
    "  2004 366   2453371.0 201.0 211.0 202.0 212.0 203.0 213.0 204.0 214.0  1B11\n"
    "  2005   1   2453372.0 301.0 311.0 302.0 312.0 303.0 313.0 304.0 314.0  1B11\n"
    "  2005   2   2453373.0 401.0 411.0 402.0 412.0 403.0 413.0 404.0 414.0  1B11\n"
    "  2005   3   2453374.0 501.0 511.0 502.0 512.0 503.0 513.0 504.0 514.0  1B11\n"
    "  2005   4   2453375.0 601.0 611.0 602.0 612.0 603.0 613.0 604.0 614.0  1B11\n"
    "  2005   5   2453376.0 701.0 711.0 702.0 712.0 703.0 713.0 704.0 714.0  1B11\n";

static void test_indices(void) {
  myb_jb2008_solar* solar = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_jb2008_solar_load(solar_text, sizeof solar_text - 1, NULL, &solar, &diag) ==
           MYB_OK);
  if (solar == NULL) {
    return;
  }
  /* 2005-01-05 (MJD 53375) at 07:12: F10 and S10 of 01-04, M10 of 01-03, Y10 of 2004-12-31 */
  myb_jb2008_indices x = {.dtc = 7.0};
  MT_CHECK(myb_jb2008_solar_at(solar, 53375.3, &x) == MYB_OK);
  MT_CHECK(x.f10 == 601.0 && x.f10_average == 611.0 && x.s10 == 602.0 && x.s10_average == 612.0 &&
           x.m10 == 503.0 && x.m10_average == 513.0 && x.y10 == 204.0 && x.y10_average == 214.0 &&
           x.dtc == 7.0);
  MT_CHECK(myb_jb2008_solar_at(solar, 53373.9, &x) == MYB_E_RANGE);
  MT_CHECK(myb_jb2008_solar_at(solar, 53377.0, &x) == MYB_E_RANGE);
  myb_jb2008_solar_destroy(solar);
  const char gap[] = "  2005   1   2453372.0 1 2 3 4 5 6 7 8 1B11\n"
                     "  2005   3   2453374.0 1 2 3 4 5 6 7 8 1B11\n";
  solar = NULL;
  MT_CHECK(myb_jb2008_solar_load(gap, sizeof gap - 1, NULL, &solar, &diag) == MYB_E_FORMAT &&
           diag.line == 2 && solar == NULL);
  /* DTCFILE: the value of the hour */
  char storm_text[512];
  int at = snprintf(storm_text, sizeof storm_text, "DTC 2005   1");
  for (int h = 0; h < 24; ++h) {
    at += snprintf(storm_text + at, sizeof storm_text - (size_t)at, " %3d", 10 + h);
  }
  at += snprintf(storm_text + at, sizeof storm_text - (size_t)at, "\nDTC 2005   2");
  for (int h = 0; h < 24; ++h) {
    at += snprintf(storm_text + at, sizeof storm_text - (size_t)at, " %3d", 100 + h);
  }
  at += snprintf(storm_text + at, sizeof storm_text - (size_t)at, "\n");
  myb_jb2008_storm* storm = NULL;
  MT_CHECK(myb_jb2008_storm_load(storm_text, (size_t)at, NULL, &storm, &diag) == MYB_OK);
  if (storm != NULL) {
    double dtc = 0.0;
    MT_CHECK(myb_jb2008_storm_at(storm, 53371.0 + 13.7 / 24.0, &dtc) == MYB_OK && dtc == 23.0);
    MT_CHECK(myb_jb2008_storm_at(storm, 53372.0, &dtc) == MYB_OK && dtc == 100.0);
    MT_CHECK(myb_jb2008_storm_at(storm, 53373.0, &dtc) == MYB_E_RANGE);
    myb_jb2008_storm_destroy(storm);
  }
  const char bad[] = "DTX 2005   1 1 2 3\n";
  MT_CHECK(myb_jb2008_storm_load(bad, sizeof bad - 1, NULL, &storm, &diag) == MYB_E_FORMAT);
}

static void test_properties(void) {
  const myb_jb2008_indices x = {.f10 = 140.0, .f10_average = 130.0, .s10 = 135.0,
                                .s10_average = 128.0, .m10 = 132.0, .m10_average = 129.0,
                                .y10 = 126.0, .y10_average = 127.0, .dtc = 30.0};
  myb_jb2008_output out;
  MT_CHECK(myb_jb2008(54000.5, 1.0, 0.2, 1.5, 0.3, 89.0, &x, &out) == MYB_E_DOMAIN);
  /* the temperature through the joins of the corrections */
  static const double joins[6] = {120.0, 200.0, 240.0, 300.0, 600.0, 800.0};
  double worst_jump = 0.0;
  for (int k = 0; k < 6; ++k) {
    myb_jb2008_output below;
    myb_jb2008_output above;
    MT_CHECK(myb_jb2008(54000.5, 1.0, 0.2, 1.5, 0.3, joins[k] - 1e-7, &x, &below) == MYB_OK);
    MT_CHECK(myb_jb2008(54000.5, 1.0, 0.2, 1.5, 0.3, joins[k] + 1e-7, &x, &above) == MYB_OK);
    worst_jump = fmax(worst_jump, fabs(above.temperature - below.temperature));
  }
  MT_CHECK(worst_jump < 1e-4);
  /* the storm term adds to the exospheric temperature; the density falls with height */
  myb_jb2008_indices quiet = x;
  quiet.dtc = 0.0;
  myb_jb2008_output calm;
  MT_CHECK(myb_jb2008(54000.5, 1.0, 0.2, 1.5, 0.3, 400.0, &x, &out) == MYB_OK);
  MT_CHECK(myb_jb2008(54000.5, 1.0, 0.2, 1.5, 0.3, 400.0, &quiet, &calm) == MYB_OK);
  MT_CHECK(fabs(out.exospheric_temperature - calm.exospheric_temperature - 30.0) < 1e-9 &&
           out.density > calm.density);
  double previous = 1.0;
  int falling = 1;
  for (double z = 90.0; z <= 4000.0; z += 7.3) {
    MT_CHECK(myb_jb2008(54000.5, 1.0, 0.2, 1.5, 0.3, z, &x, &out) == MYB_OK);
    falling = falling && out.density < previous;
    previous = out.density;
  }
  MT_CHECK(falling);
}

/* The atmosphere of a force model against the model called with the place, the Sun and the
 * sidereal time worked out here as the header describes them. */
static void test_atmosphere(void) {
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
  myb_jb2008_solar* solar = NULL;
  MT_CHECK(myb_jb2008_solar_load(solar_text, sizeof solar_text - 1, NULL, &solar, &diag) ==
           MYB_OK);
  const myb_ellipsoid wgs84 = {6378137.0, 1.0 / 298.257223563};
  /* 2005-01-05 07:12 UTC, 420 km above 51.6 N, 30 E */
  const myb_calendar utc = {.year = 2005, .month = 1, .day = 5, .hour = 7, .minute = 12};
  myb_time t;
  MT_CHECK(myb_timesys_from_utc(sys, &utc, tt, &t) == MYB_OK);
  const double degree = pi / 180.0;
  double r[3];
  MT_CHECK(myb_geodetic_to_cartesian(&wgs84, 30.0 * degree, 51.6 * degree, 420.0e3, r) ==
           MYB_OK);
  myb_jb2008_atmosphere air = {.solar = solar, .fixed = {.dtc = 12.0}, .timesys = sys,
                                     .scale = tt, .ellipsoid = wgs84};
  double density = 0.0;
  MT_CHECK(myb_jb2008_density(&air, t, r, &density) == MYB_OK);
  /* by hand */
  const double mjd = 53375.0 + 7.2 / 24.0;
  const double n = mjd - 51544.5;
  const double l = (280.460 + 0.9856474 * n) * degree;
  const double g = (357.528 + 0.9856003 * n) * degree;
  const double lambda = l + (1.915 * sin(g) + 0.020 * sin(2.0 * g)) * degree;
  const double eps = (23.439 - 0.0000004 * n) * degree;
  const double ra = atan2(cos(eps) * sin(lambda), cos(lambda));
  const double dec = asin(sin(eps) * sin(lambda));
  const double tc = n / 36525.0;
  double gmst = 280.46061837 + 360.98564736629 * n + 0.000387933 * tc * tc - tc * tc * tc / 38710000.0;
  gmst -= 360.0 * floor(gmst / 360.0);
  double lon = 0.0;
  double lat = 0.0;
  double h = 0.0;
  MT_CHECK(myb_geodetic_from_cartesian(&wgs84, r, &lon, &lat, &h) == MYB_OK);
  const myb_jb2008_indices x = {.f10 = 601.0, .f10_average = 611.0, .s10 = 602.0,
                                .s10_average = 612.0, .m10 = 503.0, .m10_average = 513.0,
                                .y10 = 104.0, .y10_average = 114.0, .dtc = 12.0};
  myb_jb2008_output out;
  MT_CHECK(myb_jb2008(mjd, ra, dec, gmst * degree + lon, atan2(r[2], hypot(r[0], r[1])), h / 1000.0,
                      &x, &out) == MYB_OK);
  printf("the model as an atmosphere, 420 km above 51.6 N on 2005-01-05 07:12 UTC: %.6e kg/m^3 "
         "(called by hand: %.6e)\n",
         density, out.density);
  MT_CHECK(fabs(density / out.density - 1.0) < 1e-6);
  myb_jb2008_solar_destroy(solar);
  myb_timesys_destroy(sys);
}

/* Against the outputs of the program of SET, if they were made locally. */
static void test_black_box(void) {
  FILE* f = fopen(MYB_DATA_DIR "/download/jb2008_set_outputs.txt", "r");
  if (f == NULL) {
    printf("data/download/jb2008_set_outputs.txt is not there (outputs of the program of SET, "
           "made locally): not compared\n");
    return;
  }
  char line[1024];
  int count = 0;
  double worst_tinf = 0.0;
  double worst_t = 0.0;
  double worst_rho = 0.0;
  double sum_rho = 0.0;
  for (int guard = 0; guard < 100000 && fgets(line, sizeof line, f) != NULL; ++guard) {
    if (line[0] == '#') {
      continue;
    }
    double v[18];
    char* p = line;
    int k = 0;
    for (; k < 18; ++k) {
      char* end = NULL;
      v[k] = strtod(p, &end);
      if (end == p) {
        break;
      }
      p = end;
    }
    if (k != 18) {
      continue;
    }
    const myb_jb2008_indices x = {.f10 = v[6], .f10_average = v[7], .s10 = v[8],
                                  .s10_average = v[9], .m10 = v[10], .m10_average = v[11],
                                  .y10 = v[12], .y10_average = v[13], .dtc = v[14]};
    myb_jb2008_output out;
    MT_CHECK(myb_jb2008(v[0], v[1], v[2], v[3], v[4], v[5], &x, &out) == MYB_OK);
    worst_tinf = fmax(worst_tinf, fabs(out.exospheric_temperature - v[15]));
    worst_t = fmax(worst_t, fabs(out.temperature - v[16]));
    const double d = fabs(log10(out.density / v[17]));
    worst_rho = fmax(worst_rho, d);
    sum_rho += d * d;
    ++count;
  }
  fclose(f);
  printf("against the program of SET (%d inputs, 90 to 6000 km): exospheric temperature %.1e K, "
         "temperature %.1e K, log10 rho %.1e at worst, %.1e rms\n",
         count, worst_tinf, worst_t, worst_rho, sqrt(sum_rho / (count > 0 ? count : 1)));
  MT_CHECK(count > 500 && worst_tinf < 1e-5 && worst_t < 1e-3 && worst_rho < 6e-3);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_static();
  test_indices();
  test_properties();
  test_atmosphere();
  test_black_box();
  return mt_finish("test_jb2008");
}
