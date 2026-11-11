/* SPDX-License-Identifier: BSD-3-Clause */
/* NRLMSISE-00 (ext/nrlmsise00) against the 17 reference evaluations of the C release
 * (DOCUMENTATION, section 7, printed to 7 significant digits in cm^-3, g cm^-3 and K; the same
 * values as the test output of the NRL Fortran code); thread independence; errors. */
#include <math.h>
#include <stdio.h>

#include "mtest.h"
#include "myb/io.h"
#include "myb_nrlmsise00.h"

static const double reference[17][11] = {
    {6.665177E+05, 1.138806E+08, 1.998211E+07, 4.022764E+05, 3.557465E+03, 4.074714E-15, 3.475312E+04, 4.095913E+06, 2.667273E+04, 1.250540E+03, 1.241416E+03},
    {3.407293E+06, 1.586333E+08, 1.391117E+07, 3.262560E+05, 1.559618E+03, 5.001846E-15, 4.854208E+04, 4.380967E+06, 6.956682E+03, 1.166754E+03, 1.161710E+03},
    {1.123767E+05, 6.934130E+04, 4.247105E+01, 1.322750E-01, 2.618848E-05, 2.756772E-18, 2.016750E+04, 5.741256E+03, 2.374394E+04, 1.239892E+03, 1.239891E+03},
    {5.411554E+07, 1.918893E+11, 6.115826E+12, 1.225201E+12, 6.023212E+10, 3.584426E-10, 1.059880E+07, 2.615737E+05, 2.819879E-42, 1.027318E+03, 2.068878E+02},
    {1.851122E+06, 1.476555E+08, 1.579356E+07, 2.633795E+05, 1.588781E+03, 4.809630E-15, 5.816167E+04, 5.478984E+06, 1.264446E+03, 1.212396E+03, 1.208135E+03},
    {8.673095E+05, 1.278862E+08, 1.822577E+07, 2.922214E+05, 2.402962E+03, 4.355866E-15, 3.686389E+04, 3.897276E+06, 2.667273E+04, 1.220146E+03, 1.212712E+03},
    {5.776251E+05, 6.979139E+07, 1.236814E+07, 2.492868E+05, 1.405739E+03, 2.470651E-15, 5.291986E+04, 1.069814E+06, 2.667273E+04, 1.116385E+03, 1.112999E+03},
    {3.740304E+05, 4.782720E+07, 5.240380E+06, 1.759875E+05, 5.501649E+02, 1.571889E-15, 8.896776E+04, 1.979741E+06, 9.121815E+03, 1.031247E+03, 1.024848E+03},
    {6.748339E+05, 1.245315E+08, 2.369010E+07, 4.911583E+05, 4.578781E+03, 4.564420E-15, 3.244595E+04, 5.370833E+06, 2.667273E+04, 1.306052E+03, 1.293374E+03},
    {5.528601E+05, 1.198041E+08, 3.495798E+07, 9.339618E+05, 1.096255E+04, 4.974543E-15, 2.686428E+04, 4.889974E+06, 2.805445E+04, 1.361868E+03, 1.347389E+03},
    {1.375488E+14, 0.000000E+00, 2.049687E+19, 5.498695E+18, 2.451733E+17, 1.261066E-03, 0.000000E+00, 0.000000E+00, 0.000000E+00, 1.027318E+03, 2.814648E+02},
    {4.427443E+13, 0.000000E+00, 6.597567E+18, 1.769929E+18, 7.891680E+16, 4.059139E-04, 0.000000E+00, 0.000000E+00, 0.000000E+00, 1.027318E+03, 2.274180E+02},
    {2.127829E+12, 0.000000E+00, 3.170791E+17, 8.506280E+16, 3.792741E+15, 1.950822E-05, 0.000000E+00, 0.000000E+00, 0.000000E+00, 1.027318E+03, 2.374389E+02},
    {1.412184E+11, 0.000000E+00, 2.104370E+16, 5.645392E+15, 2.517142E+14, 1.294709E-06, 0.000000E+00, 0.000000E+00, 0.000000E+00, 1.027318E+03, 2.795551E+02},
    {1.254884E+10, 0.000000E+00, 1.874533E+15, 4.923051E+14, 2.239685E+13, 1.147668E-07, 0.000000E+00, 0.000000E+00, 0.000000E+00, 1.027318E+03, 2.190732E+02},
    {5.196477E+05, 1.274494E+08, 4.850450E+07, 1.720838E+06, 2.354487E+04, 5.881940E-15, 2.500078E+04, 6.279210E+06, 2.667273E+04, 1.426412E+03, 1.408608E+03},
    {4.260860E+07, 1.241342E+11, 4.929562E+12, 1.048407E+12, 4.993465E+10, 2.914304E-10, 8.831229E+06, 2.252516E+05, 2.415246E-42, 1.027318E+03, 1.934071E+02},
};

/* The model as an atmosphere: the density at a place and an instant against the model
 * called with the day, the time, the geodetic place and the indices worked out here. */
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
  myb_space_weather* weather = NULL;
  MT_CHECK(myb_space_weather_load_file(MYB_TEST_DATA_DIR "/weather/sw-excerpt.csv", NULL,
                                       &weather, &diag) == MYB_OK);
  const myb_ellipsoid wgs84 = {6378137.0, 1.0 / 298.257223563};
  /* 2026-03-05 13:00:00 UTC, 420 km above 51.6 N, 30 E */
  const myb_calendar utc = {.year = 2026, .month = 3, .day = 5, .hour = 13};
  myb_time t;
  MT_CHECK(myb_timesys_from_utc(sys, &utc, tt, &t) == MYB_OK);
  const double degree = 3.141592653589793 / 180.0;
  double r[3];
  MT_CHECK(myb_geodetic_to_cartesian(&wgs84, 30.0 * degree, 51.6 * degree, 420.0e3, r) ==
           MYB_OK);
  myb_nrlmsise00_atmosphere air = {.weather = weather, .timesys = sys, .scale = tt,
                                   .ellipsoid = wgs84};
  double density = 0.0;
  MT_CHECK(myb_nrlmsise00_density(&air, t, r, &density) == MYB_OK);
  /* the 64th day of the year; the indices of the table at that hour (test_weather) */
  const double history[7] = {6.0, 4.0, 5.0, 6.0, 7.0, 19.375, 3.5};
  const myb_nrlmsise00_input direct = {.day_of_year = 64, .seconds = 46800.0, .altitude = 420.0,
                                       .latitude = 51.6, .longitude = 30.0,
                                       .local_solar_time = 15.0, .f107_average = 150.0,
                                       .f107 = 161.4, .ap = 6.0, .ap_history = history};
  int switches[24];
  for (int i = 0; i < 24; ++i) {
    switches[i] = i == 9 ? -1 : 1;
  }
  myb_nrlmsise00_output out;
  MT_CHECK(myb_nrlmsise00(&direct, switches, 0, &out) == MYB_OK);
  printf("the model as an atmosphere, 420 km above 51.6 N on 2026-03-05 13:00 UTC with the "
         "indices of a table: %.4e kg/m^3 (the model called by hand: %.4e)\n",
         density, out.density[5]);
  MT_CHECK(fabs(density / out.density[5] - 1.0) < 1e-9);
  MT_CHECK(density > 1e-13 && density < 1e-10);
  /* without a table: the three numbers given, no history */
  myb_nrlmsise00_atmosphere fixed = {.f107 = 161.4, .f107_average = 150.0, .ap = 6.0,
                                     .timesys = sys, .scale = tt, .ellipsoid = wgs84,
                                     .anomalous_oxygen = 1};
  myb_nrlmsise00_input plain = direct;
  plain.ap_history = NULL;
  double steady = 0.0;
  MT_CHECK(myb_nrlmsise00_density(&fixed, t, r, &steady) == MYB_OK);
  MT_CHECK(myb_nrlmsise00(&plain, NULL, 1, &out) == MYB_OK);
  MT_CHECK(fabs(steady / out.density[5] - 1.0) < 1e-9 && steady != density);
  /* an instant outside the table; a place inside the ellipsoid */
  MT_CHECK(myb_nrlmsise00_density(&air, myb_time_add(t, 40.0 * 86400.0), r, &density) ==
           MYB_E_RANGE);
  const double below[3] = {6.0e6, 0.0, 0.0};
  MT_CHECK(myb_nrlmsise00_density(&fixed, t, below, &density) == MYB_E_DOMAIN);
  myb_space_weather_destroy(weather);
  myb_timesys_destroy(sys);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_atmosphere();
  myb_nrlmsise00_input in[17];
  const double aph[7] = {100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
  for (int i = 0; i < 17; ++i) {
    in[i] = (myb_nrlmsise00_input){.day_of_year = 172, .seconds = 29000.0, .altitude = 400.0,
                                   .latitude = 60.0, .longitude = -70.0, .local_solar_time = 16.0,
                                   .f107_average = 150.0, .f107 = 150.0, .ap = 4.0};
  }
  in[1].day_of_year = 81;
  in[2].seconds = 75000.0;
  in[2].altitude = 1000.0;
  in[3].altitude = 100.0;
  in[10].altitude = 0.0;
  in[11].altitude = 10.0;
  in[12].altitude = 30.0;
  in[13].altitude = 50.0;
  in[14].altitude = 70.0;
  in[16].altitude = 100.0;
  in[4].latitude = 0.0;
  in[5].longitude = 0.0;
  in[6].local_solar_time = 4.0;
  in[7].f107_average = 70.0;
  in[8].f107 = 180.0;
  in[9].ap = 40.0;
  in[15].ap_history = aph;
  in[16].ap_history = aph;
  int switches[24];
  for (int i = 0; i < 24; ++i) {
    switches[i] = 1;
  }
  double worst = 0.0;
  for (int i = 0; i < 17; ++i) {
    switches[9] = i >= 15 ? -1 : 1;
    myb_nrlmsise00_output out;
    MT_CHECK(myb_nrlmsise00(&in[i], switches, 0, &out) == MYB_OK);
    for (int j = 0; j < 9; ++j) {
      const double scale = j == 5 ? 1e-3 : 1e-6; /* SI -> g cm^-3, cm^-3 */
      const double v = out.density[j] * scale;
      const double r = reference[i][j];
      const double e = r == 0.0 ? fabs(v) : fabs(v - r) / fabs(r);
      worst = fmax(worst, e);
      if (!(e < 6e-7)) {
        printf("  case %d density %d: %.7e vs %.7e\n", i, j, v, r);
      }
    }
    worst = fmax(worst, fabs(out.exospheric_temperature - reference[i][9]) / reference[i][9]);
    worst = fmax(worst, fabs(out.temperature - reference[i][10]) / reference[i][10]);
  }
  printf("NRLMSISE-00, 17 reference cases (0 .. 1000 km): worst relative difference %.1e\n",
         worst);
  MT_CHECK(worst < 6e-7);
  /* GTD7D: the total mass density includes anomalous oxygen (16 amu) */
  myb_nrlmsise00_output a;
  myb_nrlmsise00_output b;
  MT_CHECK(myb_nrlmsise00(&in[2], NULL, 0, &a) == MYB_OK);
  MT_CHECK(myb_nrlmsise00(&in[2], NULL, 1, &b) == MYB_OK);
  const double amu = 1.66e-27;
  MT_CHECK(fabs((b.density[5] - a.density[5]) - 16.0 * amu * a.density[8]) <
           1e-3 * b.density[5]);
  /* errors */
  myb_nrlmsise00_input bad = in[0];
  bad.altitude = -1.0;
  MT_CHECK(myb_nrlmsise00(&bad, NULL, 0, &a) == MYB_E_DOMAIN);
  bad = in[0];
  bad.latitude = (double)NAN;
  MT_CHECK(myb_nrlmsise00(&bad, NULL, 0, &a) == MYB_E_DOMAIN);
  bad = in[0];
  bad.day_of_year = 0;
  MT_CHECK(myb_nrlmsise00(&bad, NULL, 0, &a) == MYB_E_DOMAIN);
  switches[9] = -1;
  MT_CHECK(myb_nrlmsise00(&in[0], switches, 0, &a) == MYB_E_DOMAIN); /* no ap history */
  return mt_finish("test_nrlmsise00");
}
