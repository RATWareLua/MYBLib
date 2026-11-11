/* SPDX-License-Identifier: BSD-3-Clause */
/* The motion of a point of the crust (myb/site.h) against the test cases that come with the
 * IERS Conventions 2010, chapter 7: the tide of the solid Earth (the header of
 * DEHANTTIDEINEL.F), the load of the ocean tides (the header of HARDISP.F), the load of the
 * ocean pole tide (opoleloadcmcor.test, with nine nodes of the grid of coefficients); the
 * pole tide, the permanent part and the atmospheric tides against the equations themselves. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/env.h"
#include "myb/estim.h"
#include "myb/frames.h"
#include "myb/site.h"

#define SITE_KERNEL MYB_DATA_DIR "/frames/iers2010-site.tk"
#define BLQ MYB_TEST_DATA_DIR "/site/loading.blq"
#define GRID MYB_TEST_DATA_DIR "/site/opole-excerpt.txt"

static const double pi = 3.141592653589793;
static const double arcsecond = 3.141592653589793 / 648000.0;

static myb_site_tides* model(const char* extra) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const myb_status loaded = myb_kpool_load_file(pool, SITE_KERNEL, &diag);
  if (loaded != MYB_OK) {
    printf("the kernel: status %d, line %ld, %s\n", (int)loaded, diag.line, diag.message);
  }
  MT_CHECK(loaded == MYB_OK);
  if (extra != NULL) {
    MT_CHECK(myb_kpool_load_text(pool, extra, strlen(extra), &diag) == MYB_OK);
  }
  myb_site_tides* tides = NULL;
  const myb_status status = myb_site_tides_create(pool, "IERS10_SITE", NULL, &tides, &diag);
  if (status != MYB_OK) {
    printf("the model: status %d, %s\n", (int)status, diag.message);
  }
  MT_CHECK(status == MYB_OK);
  myb_kpool_destroy(pool);
  return tides;
}

static myb_time at(int64_t year, int month, int day, double seconds) {
  const myb_calendar calendar = {.year = year, .month = month, .day = day};
  myb_time t;
  MT_CHECK(myb_calendar_to_time(&calendar, &t) == MYB_OK);
  return myb_time_add(t, seconds);
}

/* The tide of the solid Earth: the four cases of the routine. Its masses are ratios to that
 * of the Earth; its date is UTC taken as the angle of rotation, its arguments are on TT.
 *
 * `header` is what the header of the routine prints. The routine adds the precession since
 * J2000 to the mean longitude of the Moon after it has taken the lunar time from it (S = S +
 * PR in STEP2DIU): every diurnal wave is then turned by its multiple of 1.4 degrees a
 * century — K1, whose argument is the sidereal angle chi itself, by 0.13 degrees in 2009,
 * 0.03 mm of its 12 mm. The Conventions define the arguments by chi and the Delaunay
 * arguments (Table 7.3a prints both sets of multipliers), and so does the library;
 * `arithmetic` is the arithmetic of the routine without that addition, done apart in double
 * precision. For the fourth case the header repeats the numbers of the third: its Sun is at
 * 0.1 au and the tide it raises is 18 m. */
typedef struct solid_case {
  double site[3];
  double sun[3];
  double moon[3];
  int year;
  int month;
  int day;
  double tt_utc;
  double header[3];
  double arithmetic[3];
} solid_case;

static const solid_case solid_cases[4] = {
    {{4075578.385, 931852.890, 4801570.154},
     {137859926952.015, 54228127881.4350, 23509422341.6960},
     {-179996231.920342, -312468450.131567, -169288918.592160},
     2009, 4, 13, 66.184,
     {0.7700420357108125891e-01, 0.6304056321824967613e-01, 0.5516568152597246810e-01},
     {7.699135937343193e-02, 6.303680995958384e-02, 5.514868535039304e-02}},
    {{1112189.660, -4842955.026, 3985352.284},
     {-54537460436.2357, 130244288385.279, 56463429031.5996},
     {300396716.912, 243238281.451, 120548075.939},
     2012, 7, 13, 67.184,
     {-0.2036831479592075833e-01, 0.5658254776225972449e-01, -0.7597679676871742227e-01},
     {-2.037475680049865e-02, 5.660634098762001e-02, -7.599714781455122e-02}},
    {{1112200.5696, -4842957.8511, 3985345.9122},
     {100210282451.6279, 103055630398.3160, 56855096480.4475},
     {369817604.4348, 1897917.5258, 120804980.8284},
     2015, 7, 15, 68.184,
     {0.00509570869172363845, 0.0828663025983528700, -0.0636634925404189617},
     {5.088396231902818e-03, 8.289436298634660e-02, -6.368680397587144e-02}},
    {{1112152.8166, -4842857.5435, 3985496.1783},
     {8382471154.1312895, 10512408445.356153, -5360583240.3763866},
     {380934092.93550891, 2871428.1904491195, 79015680.553570181},
     2017, 1, 15, 69.184,
     {0.0050957086917236384, 0.082866302598352870, -0.063663492540418962},
     {-1.821734956506847e+01, -2.350537969200471e+01, 1.209763704210483e+01}},
};

static void test_solid(const myb_site_tides* tides) {
  double worst = 0.0;
  double worst_header = 0.0;
  for (int k = 0; k < 4; ++k) {
    const solid_case* c = &solid_cases[k];
    myb_tide_raiser raisers[2] = {{.gm = 332946.0482}, {.gm = 0.0123000371}};
    memcpy(raisers[0].r, c->sun, sizeof c->sun);
    memcpy(raisers[1].r, c->moon, sizeof c->moon);
    const myb_site_input input = {.gm = 1.0, .raisers = raisers, .count = 2,
                                  .t = at(c->year, c->month, c->day, 0.0),
                                  .dynamical = c->tt_utc,
                                  .flags = MYB_SITE_SOLID | MYB_SITE_FREQUENCY};
    double d[3];
    MT_CHECK(myb_site_displacement(tides, c->site, NULL, NULL, &input, d) == MYB_OK);
    double e = 0.0;
    double printed = 0.0;
    for (int i = 0; i < 3; ++i) {
      e = fmax(e, fabs(d[i] - c->arithmetic[i]));
      printed = fmax(printed, fabs(d[i] - c->header[i]));
    }
    printf("  the tide of the solid Earth, case %d: %+.6f %+.6f %+.6f m; from the arithmetic "
           "of the routine %.1e m, from its header %.1e m\n",
           k + 1, d[0], d[1], d[2], e, printed);
    worst = fmax(worst, e);
    if (k < 3) {
      worst_header = fmax(worst_header, printed);
    }
  }
  printf("the tide of the solid Earth (steps 1 and 2), four cases: %.1e m from the arithmetic "
         "of the routine; %.1e m from the values of its header, which turn the diurnal waves "
         "by the precession\n",
         worst, worst_header);
  MT_CHECK(worst < 2e-7);
  MT_CHECK(worst_header < 4e-5);

  /* what the parts are: a point at 45 N on the meridian of the Moon, the Moon in the plane
   * of the equator */
  const double a = 6378136.6;
  const double site[3] = {a * sqrt(0.5), 0.0, a * sqrt(0.5)};
  const myb_tide_raiser moon = {.gm = 0.0123000371, .r = {3.844e8, 0.0, 0.0}};
  myb_site_input input = {.gm = 1.0, .raisers = &moon, .count = 1, .flags = MYB_SITE_SOLID};
  double d[3];
  MT_CHECK(myb_site_displacement(tides, site, NULL, NULL, &input, d) == MYB_OK);
  /* eq. 5 with c = cos 45: along the radius h2 (3 c^2 - 1) / 2 = h2 / 4, across it towards
   * the Moon 3 l2 c sin 45 = 1.5 l2; degree 3 and the small terms are within a few percent */
  const double f2 = 0.0123000371 * a * pow(a / 3.844e8, 3.0);
  const double radial = (d[0] + d[2]) * sqrt(0.5);
  const double south = (d[0] - d[2]) * sqrt(0.5);
  printf("a point at 45 N under a Moon on the equator: %.4f m up, %.4f m to the south (eq. 5 "
         "alone: %.4f, %.4f)\n",
         radial, south, f2 * 0.6078 / 4.0, f2 * 1.5 * 0.0847);
  MT_CHECK(fabs(radial / (f2 * 0.6078 / 4.0) - 1.0) < 0.05);
  MT_CHECK(fabs(south / (f2 * 1.5 * 0.0847) - 1.0) < 0.05);
  MT_CHECK(fabs(d[1]) < 1e-2 * fabs(radial));

  /* the permanent part (eq. 14): no raiser, a mean-tide place */
  input.count = 0;
  input.flags = MYB_SITE_SOLID | MYB_SITE_MEAN_TIDE;
  double worst_permanent = 0.0;
  for (int k = 0; k <= 6; ++k) {
    const double edge = k == 0 ? 0.5 : (k == 6 ? -0.5 : 0.0);
    const double phi = (-90.0 + 30.0 * k + edge) * pi / 180.0;
    const double p[3] = {a * cos(phi), 0.0, a * sin(phi)};
    MT_CHECK(myb_site_displacement(tides, p, NULL, NULL, &input, d) == MYB_OK);
    const double p2 = 1.5 * sin(phi) * sin(phi) - 0.5;
    const double up = d[0] * cos(phi) + d[2] * sin(phi);
    const double north = -d[0] * sin(phi) + d[2] * cos(phi);
    /* taken out: the opposite of what eq. 14 adds */
    worst_permanent = fmax(worst_permanent, fabs(up + (-0.1206 + 0.0001 * p2) * p2));
    worst_permanent =
        fmax(worst_permanent, fabs(north + (-0.0252 - 0.0001 * p2) * sin(2.0 * phi)));
  }
  printf("the permanent part against eq. 14 (printed to 0.1 mm): %.1e m\n", worst_permanent);
  MT_CHECK(worst_permanent < 1e-4);

  /* refused */
  const double pole[3] = {0.0, 0.0, 6356752.0};
  const double centre[3] = {0.0, 0.0, 0.0};
  input.flags = MYB_SITE_SOLID;
  MT_CHECK(myb_site_displacement(tides, pole, NULL, NULL, &input, d) == MYB_E_DOMAIN);
  MT_CHECK(myb_site_displacement(tides, centre, NULL, NULL, &input, d) == MYB_E_DOMAIN);
  const myb_tide_raiser nowhere = {.gm = 1.0, .r = {0.0, 0.0, 0.0}};
  input.raisers = &nowhere;
  input.count = 1;
  MT_CHECK(myb_site_displacement(tides, site, NULL, NULL, &input, d) == MYB_E_DOMAIN);
  input.count = 0;
  input.gm = 0.0;
  MT_CHECK(myb_site_displacement(tides, site, NULL, NULL, &input, d) == MYB_E_DOMAIN);
}

/* The pole tide: eq. 26 to 28 at a point, with the secular pole of eq. 21. */
static void test_pole(const myb_site_tides* tides) {
  const double lambda = 0.7;
  const double phi = -0.4; /* geocentric */
  const double a = 6371000.0;
  const double site[3] = {a * cos(phi) * cos(lambda), a * cos(phi) * sin(lambda), a * sin(phi)};
  const myb_time t = at(2015, 6, 1, 0.0);
  const double xp = 0.323 * arcsecond;
  const double yp = 0.156 * arcsecond;
  const myb_site_input input = {.gm = 1.0, .t = t, .xp = xp, .yp = yp, .flags = MYB_SITE_POLE};
  double d[3];
  MT_CHECK(myb_site_displacement(tides, site, NULL, NULL, &input, d) == MYB_OK);
  const double years = myb_time_seconds(t) / 31557600.0;
  const double m1 = 0.323 - (55.0 + 1.677 * years) * 1e-3;
  const double m2 = -(0.156 - (320.5 + 3.460 * years) * 1e-3);
  const double theta = pi / 2.0 - phi;
  const double s_theta = -9e-3 * cos(2.0 * theta) * (m1 * cos(lambda) + m2 * sin(lambda));
  const double s_lambda = 9e-3 * cos(theta) * (m1 * sin(lambda) - m2 * cos(lambda));
  const double s_r = -33e-3 * sin(2.0 * theta) * (m1 * cos(lambda) + m2 * sin(lambda));
  /* eq. 27: the transpose of R on (S_theta, S_lambda, S_r) */
  const double expected[3] = {
      cos(theta) * cos(lambda) * s_theta - sin(lambda) * s_lambda + sin(theta) * cos(lambda) * s_r,
      cos(theta) * sin(lambda) * s_theta + cos(lambda) * s_lambda + sin(theta) * sin(lambda) * s_r,
      -sin(theta) * s_theta + cos(theta) * s_r};
  double e = 0.0;
  for (int i = 0; i < 3; ++i) {
    e = fmax(e, fabs(d[i] - expected[i]));
  }
  printf("the pole tide against eq. 26-28 (m1 = %.3f\", m2 = %.3f\"): %+.5f %+.5f %+.5f m, "
         "%.1e m\n",
         m1, m2, d[0], d[1], d[2], e);
  MT_CHECK(e < 1e-17 && fabs(s_r) > 1e-3);
}

/* The load of the ocean tides: HARDISP at Onsala and Reykjavik from 2009-06-25 01:10:45
 * UTC, 24 hours; up, to the south, to the west as the program prints them (six decimals;
 * the program works in single precision).
 *
 * The program takes the lunar time as 360 UT - D, that is, the angle of rotation as UT plus
 * the mean longitude of the Sun, and that longitude on TT. The chi = GMST + pi of the
 * Conventions has in its place the mean Sun of the sidereal time, 20.5" behind (the constant
 * of aberration), taken on UT1: in June 2009 the program's angle is 23.6" ahead, 1.57 s of
 * rotation. Both are compared: the rotation taken that much later, which is the program's
 * arithmetic, and the argument as the Conventions define it. */
static const double program_ahead = 1.5688; /* s */
static const double hardisp[48][3] = {
    {0.003094, -0.001538, -0.000895},  {0.001812, -0.000950, -0.000193},
    {0.000218, -0.000248, 0.000421},   {-0.001104, 0.000404, 0.000741},
    {-0.001668, 0.000863, 0.000646},   {-0.001209, 0.001042, 0.000137},
    {0.000235, 0.000926, -0.000667},   {0.002337, 0.000580, -0.001555},
    {0.004554, 0.000125, -0.002278},   {0.006271, -0.000291, -0.002615},
    {0.006955, -0.000537, -0.002430},  {0.006299, -0.000526, -0.001706},
    {0.004305, -0.000244, -0.000559},  {0.001294, 0.000245, 0.000793},
    {-0.002163, 0.000819, 0.002075},   {-0.005375, 0.001326, 0.003024},
    {-0.007695, 0.001622, 0.003448},   {-0.008669, 0.001610, 0.003272},
    {-0.008143, 0.001262, 0.002557},   {-0.006290, 0.000633, 0.001477},
    {-0.003566, -0.000155, 0.000282},  {-0.000593, -0.000941, -0.000766},
    {0.001992, -0.001561, -0.001457},  {0.003689, -0.001889, -0.001680},
    {-0.005940, -0.001245, -0.000278}, {0.013516, -0.001086, 0.003212},
    {0.029599, -0.000353, 0.005483},   {0.038468, 0.000699, 0.005997},
    {0.038098, 0.001721, 0.004690},    {0.028780, 0.002363, 0.001974},
    {0.013016, 0.002371, -0.001369},   {-0.005124, 0.001653, -0.004390},
    {-0.021047, 0.000310, -0.006225},  {-0.030799, -0.001383, -0.006313},
    {-0.032056, -0.003048, -0.004549}, {-0.024698, -0.004288, -0.001314},
    {-0.010814, -0.004794, 0.002623},  {0.005849, -0.004416, 0.006291},
    {0.020857, -0.003208, 0.008766},   {0.030226, -0.001413, 0.009402},
    {0.031437, 0.000594, 0.007996},    {0.024079, 0.002389, 0.004844},
    {0.009945, 0.003606, 0.000663},    {-0.007426, 0.004022, -0.003581},
    {-0.023652, 0.003601, -0.006911},  {-0.034618, 0.002505, -0.008585},
    {-0.037515, 0.001044, -0.008270},  {-0.031544, -0.000402, -0.006125},
};

static void test_ocean(const myb_site_tides* tides) {
  MT_CHECK(myb_site_tides_waves(tides) == 342 && myb_site_tides_given_waves(tides) == 11);
  static const char* const names[2] = {"onsala", "Reykjavik"};
  myb_diag diag = {0};
  for (int k = 0; k < 2; ++k) {
    myb_ocean_loading* loading = NULL;
    MT_CHECK(myb_ocean_loading_load_file(tides, BLQ, names[k], NULL, &loading, &diag) == MYB_OK);
    double worst = 0.0;
    double defined = 0.0;
    double largest = 0.0;
    for (int hour = 0; hour < 24; ++hour) {
      const myb_time t = at(2009, 6, 25, 3600.0 + 10.0 * 60.0 + 45.0 + 3600.0 * hour);
      const double* e = hardisp[24 * k + hour];
      double uws[3];
      myb_ocean_loading_at(tides, loading, myb_time_add(t, program_ahead),
                           66.184 - program_ahead, uws);
      worst = fmax(worst, fabs(uws[0] - e[0]));
      worst = fmax(worst, fabs(uws[2] - e[1]));
      worst = fmax(worst, fabs(uws[1] - e[2]));
      myb_ocean_loading_at(tides, loading, t, 66.184, uws);
      defined = fmax(defined, fabs(uws[0] - e[0]));
      defined = fmax(defined, fabs(uws[2] - e[1]));
      defined = fmax(defined, fabs(uws[1] - e[2]));
      largest = fmax(largest, fabs(e[0]));
    }
    printf("the load of the ocean tides at %s, 24 hours (up to %.1f mm): from the program's "
           "values %.1e m with its lunar time, %.1e m with the argument of the Conventions\n",
           names[k], 1e3 * largest, worst, defined);
    MT_CHECK(worst < 1.5e-6);
    MT_CHECK(defined < 1e-5);
    /* as a part of the displacement: along the axes of the point */
    const double lon = (k == 0 ? 11.9264 : -21.9555) * pi / 180.0;
    const double lat = (k == 0 ? 57.3958 : 64.1388) * pi / 180.0;
    const myb_ellipsoid grs80 = {6378137.0, 1.0 / 298.257222101};
    double site[3];
    MT_CHECK(myb_geodetic_to_cartesian(&grs80, lon, lat, 0.0, site) == MYB_OK);
    myb_rotation enu;
    myb_topocentric_enu(lon, lat, &enu);
    const myb_site_loads loads = {.ocean = loading};
    const myb_time t = at(2009, 6, 25, 4245.0);
    const myb_site_input input = {.gm = 1.0, .t = t, .dynamical = 66.184,
                                  .flags = MYB_SITE_OCEAN};
    double d[3];
    MT_CHECK(myb_site_displacement(tides, site, enu.m, &loads, &input, d) == MYB_OK);
    double uws[3];
    myb_ocean_loading_at(tides, loading, t, 66.184, uws);
    double local[3];
    myb_rotation_apply_position(&enu, d, local);
    MT_CHECK(fabs(local[0] + uws[1]) < 1e-17 && fabs(local[1] + uws[2]) < 1e-17 &&
             fabs(local[2] - uws[0]) < 1e-17);
    /* the loads need the axes */
    MT_CHECK(myb_site_displacement(tides, site, NULL, &loads, &input, d) == MYB_E_DOMAIN);
    myb_ocean_loading_destroy(loading);
  }
  /* M2 alone, the other given waves zero: the admittance falls to nothing at N2 and S2 on
   * either side, and the band swings by the amplitude of M2 and a little of the waves next
   * to it (those of the Moon's node) */
  double amplitudes[33] = {0};
  double phases[33] = {0};
  amplitudes[0] = 0.01; /* M2, up */
  phases[0] = 30.0;
  myb_ocean_loading* one = NULL;
  MT_CHECK(myb_ocean_loading_create(tides, amplitudes, phases, NULL, &one) == MYB_OK);
  double low = 1.0;
  double high = -1.0;
  for (int k = 0; k < 4000; ++k) {
    double uws[3];
    myb_ocean_loading_at(tides, one, at(2020, 1, 1, 1223.0 * k), 69.184, uws);
    low = fmin(low, uws[0]);
    high = fmax(high, uws[0]);
    MT_CHECK(uws[1] == 0.0 && uws[2] == 0.0);
  }
  printf("M2 alone (10 mm up, the other given waves zero): the band swings from %.2f to "
         "%.2f mm\n",
         1e3 * low, 1e3 * high);
  MT_CHECK(high > 0.0100 && high < 0.0105 && low < -0.0100 && low > -0.0105);
  myb_ocean_loading_destroy(one);
  amplitudes[5] = -1.0;
  MT_CHECK(myb_ocean_loading_create(tides, amplitudes, phases, NULL, &one) == MYB_E_DOMAIN);
  /* files that are refused */
  myb_ocean_loading* none = NULL;
  MT_CHECK(myb_ocean_loading_load_file(tides, BLQ, "nowhere", NULL, &none, &diag) ==
           MYB_E_NOT_FOUND);
  MT_CHECK(myb_ocean_loading_load_file(tides, BLQ, "broken", NULL, &none, &diag) ==
               MYB_E_FORMAT &&
           diag.line == 31);
  MT_CHECK(myb_ocean_loading_load_file(tides, MYB_TEST_DATA_DIR "/site/none.blq", "onsala", NULL,
                                       &none, &diag) == MYB_E_IO);
  myb_ocean_loading_destroy(NULL);
}

/* The load of the ocean pole tide: rows of the test file of the Conventions (the modified
 * Julian date, xp, yp in arcseconds; radial, north, east in metres) at 43.75 S, 232.25 E.
 * The file was made with G = 6.673e-11 and the mean pole of the Conventions of 2003, whose
 * years it counts from the midnight before J2000.0: the dates are given that half a day
 * later here (the mean pole is all the date serves; half a day of it is 3e-8 m). */
static const double ocean_pole[6][6] = {
    {52640, -8.844100e-02, 1.879860e-01, 3.690801e-04, 4.801230e-05, 3.133391e-05},
    {52760, -6.682700e-02, 5.080160e-01, -1.175011e-03, -8.890547e-05, -2.727074e-05},
    {52890, 2.649760e-01, 3.742030e-01, 7.803496e-04, 4.528833e-05, 2.518984e-06},
    {53020, -5.402000e-03, 1.627060e-01, 8.327240e-04, 7.944129e-05, 3.795510e-05},
    {53160, -7.757900e-02, 4.723840e-01, -1.015743e-03, -7.374887e-05, -2.005379e-05},
    {53370, 1.501520e-01, 2.390430e-01, 1.053583e-03, 8.177465e-05, 2.678392e-05},
};

static void test_ocean_pole(void) {
  myb_site_tides* tides =
      model("\\begindata\n"
            "IERS10_SITE_OCEAN_POLE = ( 6.673D-11 6378136.6 1025 9.7803278 7.292115D-5\n"
            "                           3.986004418D14 0.6870 0.0036 )\n"
            "IERS10_SITE_SECULAR_POLE = ( 54.0 0.83 357.0 3.95 )\n"
            "\\begintext\n");
  myb_site_grid* grid = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_site_grid_load_file(GRID, NULL, &grid, &diag) == MYB_OK);
  MT_CHECK(myb_site_grid_values(grid) == 6);
  const double lon = 232.25 * pi / 180.0;
  const double lat = -43.75 * pi / 180.0;
  myb_site_loads loads = {0};
  MT_CHECK(myb_site_grid_at(grid, lon, lat, loads.ocean_pole) == MYB_OK);
  MT_CHECK(fabs(loads.ocean_pole[0] - 0.216192) < 1e-15 &&
           fabs(loads.ocean_pole[5] - 0.010389) < 1e-15);
  const double a = 6371000.0;
  const double site[3] = {a * cos(lat) * cos(lon), a * cos(lat) * sin(lon), a * sin(lat)};
  myb_rotation enu;
  myb_topocentric_enu(lon, lat, &enu);
  double worst = 0.0;
  for (int k = 0; k < 6; ++k) {
    const double* row = ocean_pole[k];
    const myb_site_input input = {.gm = 1.0,
                                  .t = myb_time_from_seconds((row[0] - 51544.0) * 86400.0),
                                  .xp = row[1] * arcsecond,
                                  .yp = row[2] * arcsecond,
                                  .flags = MYB_SITE_OCEAN_POLE};
    double d[3];
    MT_CHECK(myb_site_displacement(tides, site, enu.m, &loads, &input, d) == MYB_OK);
    double local[3]; /* east, north, up */
    myb_rotation_apply_position(&enu, d, local);
    worst = fmax(worst, fabs(local[2] - row[3]) / fabs(row[3]));
    worst = fmax(worst, fabs(local[1] - row[4]) / fabs(row[4]));
    worst = fmax(worst, fabs(local[0] - row[5]) / fabs(row[5]));
  }
  printf("the load of the ocean pole tide, six dates of the test file (up to 1.2 mm; seven "
         "digits printed): %.1e of the value\n",
         worst);
  MT_CHECK(worst < 2e-6);
  /* the grid between its nodes and outside them */
  double v[MYB_SITE_GRID_MAX_VALUES];
  MT_CHECK(myb_site_grid_at(grid, 232.0 * pi / 180.0, -44.0 * pi / 180.0, v) == MYB_OK);
  MT_CHECK(fabs(v[0] - 0.25 * (0.218082 + 0.216194 + 0.218095 + 0.216192)) < 1e-13);
  MT_CHECK(myb_site_grid_at(grid, 232.75 * pi / 180.0, -43.25 * pi / 180.0, v) == MYB_OK);
  MT_CHECK(fabs(v[0] - 0.214219) < 1e-13 && fabs(v[5] - 0.009972) < 1e-13);
  MT_CHECK(myb_site_grid_at(grid, 233.0 * pi / 180.0, -44.0 * pi / 180.0, v) == MYB_E_RANGE);
  MT_CHECK(myb_site_grid_at(grid, 232.0 * pi / 180.0, -43.0 * pi / 180.0, v) == MYB_E_RANGE);
  MT_CHECK(myb_site_grid_at(grid, (double)NAN, 0.0, v) == MYB_E_DOMAIN);
  myb_site_grid_destroy(grid);
  myb_site_grid_destroy(NULL);
  /* a grid round the body: across its first meridian */
  const char* round = "Number_longitude_Grid_Points = 4\nFirst_longitude_degrees = 45\n"
                      "Longitude_step_degrees = 90\nNumber_latitude_grid_points = 2\n"
                      "First_latitude_degrees = -10\nLatitude_step_degrees = 20\n-----\n"
                      "45 -10 1\n135 -10 2\n225 -10 3\n315 -10 4\n"
                      "45 10 5\n135 10 6\n225 10 7\n315 10 8\n";
  MT_CHECK(myb_site_grid_load(round, strlen(round), NULL, &grid, &diag) == MYB_OK);
  MT_CHECK(myb_site_grid_at(grid, 0.0, 0.0, v) == MYB_OK && fabs(v[0] - 4.5) < 1e-13);
  MT_CHECK(myb_site_grid_at(grid, -pi / 2.0, -10.0 * pi / 180.0, v) == MYB_OK &&
           fabs(v[0] - 3.5) < 1e-13);
  myb_site_grid_destroy(grid);
  /* grids that are refused */
  const char* brief = "Number_longitude_Grid_Points = 4\n-----\n45 -10 1\n";
  MT_CHECK(myb_site_grid_load(brief, strlen(brief), NULL, &grid, &diag) == MYB_E_FORMAT);
  const char* torn = "Number_longitude_Grid_Points = 2\nFirst_longitude_degrees = 0\n"
                     "Longitude_step_degrees = 1\nNumber_latitude_grid_points = 2\n"
                     "First_latitude_degrees = 0\nLatitude_step_degrees = 1\n-----\n"
                     "0 0 1\n1 0 2\n0 1 3\n";
  MT_CHECK(myb_site_grid_load(torn, strlen(torn), NULL, &grid, &diag) == MYB_E_FORMAT);
  const char* stray = "Number_longitude_Grid_Points = 2\nFirst_longitude_degrees = 0\n"
                      "Longitude_step_degrees = 1\nNumber_latitude_grid_points = 2\n"
                      "First_latitude_degrees = 0\nLatitude_step_degrees = 1\n-----\n"
                      "0 0 1\n1 0 2\n0 1 3\n2 1 4\n";
  MT_CHECK(myb_site_grid_load(stray, strlen(stray), NULL, &grid, &diag) == MYB_E_FORMAT &&
           diag.line == 11);
  MT_CHECK(myb_site_grid_load_file(MYB_TEST_DATA_DIR "/site/none.txt", NULL, &grid, &diag) ==
           MYB_E_IO);
  myb_site_tides_destroy(tides);
}

/* The atmospheric tides S1 and S2: eq. 19 at a time of the day. */
static void test_atmosphere(const myb_site_tides* tides) {
  const double lon = 0.3;
  const double lat = 0.9;
  const myb_ellipsoid grs80 = {6378137.0, 1.0 / 298.257222101};
  double site[3];
  MT_CHECK(myb_geodetic_to_cartesian(&grs80, lon, lat, 100.0, site) == MYB_OK);
  myb_rotation enu;
  myb_topocentric_enu(lon, lat, &enu);
  const myb_site_loads loads = {.atmosphere = {1.1e-3, -0.4e-3, 0.6e-3, 0.2e-3, 1e-4, 2e-4, -3e-4,
                                               4e-4, -5e-5, 6e-5, 7e-5, -8e-5}};
  /* 07:30 of UT1 */
  const double day = 7.5 / 24.0;
  const myb_site_input input = {.gm = 1.0, .t = at(2021, 3, 4, 7.5 * 3600.0),
                                .flags = MYB_SITE_ATMOSPHERE};
  double d[3];
  MT_CHECK(myb_site_displacement(tides, site, enu.m, &loads, &input, d) == MYB_OK);
  double local[3];
  myb_rotation_apply_position(&enu, d, local);
  const double c1 = cos(2.0 * pi * day);
  const double s1 = sin(2.0 * pi * day);
  const double c2 = cos(4.0 * pi * day);
  const double s2 = sin(4.0 * pi * day);
  const double* a = loads.atmosphere;
  MT_CHECK(fabs(local[2] - (a[0] * c1 + a[1] * s1 + a[2] * c2 + a[3] * s2)) < 1e-17);
  MT_CHECK(fabs(local[0] - (a[4] * c1 + a[5] * s1 + a[6] * c2 + a[7] * s2)) < 1e-17);
  MT_CHECK(fabs(local[1] - (a[8] * c1 + a[9] * s1 + a[10] * c2 + a[11] * s2)) < 1e-17);
  /* no loads: nothing */
  MT_CHECK(myb_site_displacement(tides, site, NULL, NULL, &input, d) == MYB_OK);
  MT_CHECK(d[0] == 0.0 && d[1] == 0.0 && d[2] == 0.0);
}

/* A station that moves: a planet that turns, with a rotation time 69 s behind the time of
 * the model and a pole off its secular place, and two bodies going round it. */
typedef struct turning {
  double rate;
} turning;

static myb_status turning_orientation(void* context, myb_time t, myb_orientation* o) {
  const turning* s = context;
  const double angles[3] = {0.1, 0.05, 2.0 + s->rate * myb_time_seconds(t)};
  const double rates[3] = {0.0, 0.0, s->rate};
  myb_rotation_euler313(angles, rates, &o->rotation);
  o->ut = myb_time_add(t, -69.184);
  o->xp = 0.2 * arcsecond;
  o->yp = 0.4 * arcsecond;
  return MYB_OK;
}

typedef struct circling {
  double radius;
  double rate;
  double tilt;
} circling;

static myb_status circling_state(void* context, myb_time t, double state[6]) {
  const circling* c = context;
  const double angle = 0.7 + c->rate * myb_time_seconds(t);
  const double x = c->radius * myb_cos(angle);
  const double y = c->radius * myb_sin(angle);
  const double ct = myb_cos(c->tilt);
  const double st = myb_sin(c->tilt);
  state[0] = x;
  state[1] = y * ct;
  state[2] = y * st;
  state[3] = -c->rate * y;
  state[4] = c->rate * x * ct;
  state[5] = c->rate * x * st;
  return MYB_OK;
}

static void test_station(const myb_site_tides* tides) {
  turning spin = {7.292115e-5};
  circling moon = {3.844e8, 2.6617e-6, 0.4};
  circling sun = {1.496e11, 1.991e-7, 0.41};
  const myb_env_options options = {.epoch = at(2018, 3, 1, 0.0)};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  size_t planet = 0;
  size_t raisers[2] = {0, 0};
  const myb_env_body earth = {.gm = 3.986004418e14, .orientation = turning_orientation,
                              .orientation_context = &spin};
  const myb_env_body second = {.gm = 4.9028e12, .state = circling_state, .state_context = &moon};
  const myb_env_body third = {.gm = 1.32712440018e20, .state = circling_state,
                              .state_context = &sun};
  MT_CHECK(myb_env_add_body(env, &earth, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &second, &raisers[0]) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &third, &raisers[1]) == MYB_OK);
  const myb_ellipsoid grs80 = {6378137.0, 1.0 / 298.257222101};
  myb_station fixed;
  myb_station moving;
  MT_CHECK(myb_station_init(&fixed, env, planet, &grs80, 0.4, 0.8, 250.0) == MYB_OK);
  MT_CHECK(myb_station_init(&moving, env, planet, &grs80, 0.4, 0.8, 250.0) == MYB_OK);
  MT_CHECK(fixed.motion == NULL);
  const myb_site_loads loads = {.atmosphere = {1e-3, 2e-4, -5e-4, 3e-4}};
  const myb_station_motion motion = {.tides = tides, .loads = &loads,
                                     .raisers = {raisers[0], raisers[1]}, .raiser_count = 2,
                                     .flags = MYB_SITE_ALL, .rate = 1};
  moving.motion = &motion;
  const double t = 12345.0;
  double a[6];
  double b[6];
  MT_CHECK(myb_station_state(&fixed, t, a) == MYB_OK);
  MT_CHECK(myb_station_state(&moving, t, b) == MYB_OK);
  /* the same displacement by hand: the raisers in the body-fixed frame, UT1 and the pole of
   * the orientation, TT - UT1 */
  myb_orientation o;
  MT_CHECK(turning_orientation(&spin, myb_time_add(options.epoch, t), &o) == MYB_OK);
  myb_tide_raiser by_hand[2];
  double s[6];
  MT_CHECK(circling_state(&moon, myb_time_add(options.epoch, t), s) == MYB_OK);
  myb_rotation_apply_position(&o.rotation, s, by_hand[0].r);
  by_hand[0].gm = second.gm;
  MT_CHECK(circling_state(&sun, myb_time_add(options.epoch, t), s) == MYB_OK);
  myb_rotation_apply_position(&o.rotation, s, by_hand[1].r);
  by_hand[1].gm = third.gm;
  const myb_site_input input = {.gm = earth.gm, .raisers = by_hand, .count = 2, .t = o.ut,
                                .dynamical = 69.184, .xp = o.xp, .yp = o.yp,
                                .flags = MYB_SITE_ALL};
  double d[3];
  MT_CHECK(myb_site_displacement(tides, fixed.position, fixed.enu, &loads, &input, d) == MYB_OK);
  const double moved[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
  double in_body[3];
  myb_rotation_apply_position(&o.rotation, moved, in_body);
  double e = 0.0;
  for (int i = 0; i < 3; ++i) {
    e = fmax(e, fabs(in_body[i] - d[i]));
  }
  const double size = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  printf("a station that moves: %.4f m from its place; against the displacement by hand "
         "%.1e m\n",
         size, e);
  MT_CHECK(e < 5e-9 && size > 0.01 && size < 0.6);
  /* its velocity: what the displacement adds to the velocity of the fixed station against
   * the difference of what it adds to the place a second before and after (whole seconds:
   * the angle of the planet above is exact at them) */
  double before[6];
  double after[6];
  double fixed_before[6];
  double fixed_after[6];
  MT_CHECK(myb_station_state(&moving, t - 1.0, before) == MYB_OK);
  MT_CHECK(myb_station_state(&moving, t + 1.0, after) == MYB_OK);
  MT_CHECK(myb_station_state(&fixed, t - 1.0, fixed_before) == MYB_OK);
  MT_CHECK(myb_station_state(&fixed, t + 1.0, fixed_after) == MYB_OK);
  double worst_rate = 0.0;
  double rate = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double v = ((after[i] - fixed_after[i]) - (before[i] - fixed_before[i])) / 2.0;
    worst_rate = fmax(worst_rate, fabs((b[3 + i] - a[3 + i]) - v));
    rate = fmax(rate, fabs(b[3 + i] - a[3 + i]));
  }
  printf("  the velocity the displacement adds, %.1e m/s, against the difference of the "
         "places: %.1e m/s\n",
         rate, worst_rate);
  MT_CHECK(worst_rate < 2e-9 && rate > 1e-6);
  /* without the rate the velocity is that of the place alone */
  myb_station_motion still = motion;
  still.rate = 0;
  moving.motion = &still;
  double c[6];
  MT_CHECK(myb_station_state(&moving, t, c) == MYB_OK);
  MT_CHECK(c[0] == b[0] && c[1] == b[1] && c[2] == b[2] && c[3] != b[3]);
  /* the measurements take the moved place: the range to a point above the station */
  const double y[6] = {a[0] * 1.1, a[1] * 1.1, a[2] * 1.1, 0.0, 0.0, 0.0};
  double near = 0.0;
  double far = 0.0;
  MT_CHECK(myb_measure_range(&fixed, t, y, &near, NULL) == MYB_OK);
  MT_CHECK(myb_measure_range(&moving, t, y, &far, NULL) == MYB_OK);
  const double up = (moved[0] * a[0] + moved[1] * a[1] + moved[2] * a[2]) /
                    sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  printf("  the range to a point above it changes by %.4f m (the displacement along the "
         "radius: %.4f m)\n",
         near - far, up);
  MT_CHECK(fabs((near - far) - up) < 1e-6);
  /* a raiser that is the body itself: refused */
  myb_station_motion wrong = motion;
  wrong.raisers[0] = planet;
  moving.motion = &wrong;
  MT_CHECK(myb_station_state(&moving, t, c) == MYB_E_DOMAIN);
  myb_env_destroy(env);
}

static void test_kernels(void) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  myb_site_tides* tides = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_site_tides_create(pool, "IERS10_SITE", NULL, &tides, &diag) == MYB_E_NOT_FOUND);
  MT_CHECK(myb_kpool_load_file(pool, SITE_KERNEL, &diag) == MYB_OK);
  MT_CHECK(myb_site_tides_create(pool, "A_PREFIX_THAT_IS_TOO_LONG", NULL, &tides, &diag) ==
           MYB_E_DOMAIN);
  const char* odd = "\\begindata\nIERS10_SITE_LOVE = ( 0.6 0.08 )\n\\begintext\n";
  MT_CHECK(myb_kpool_load_text(pool, odd, strlen(odd), &diag) == MYB_OK);
  MT_CHECK(myb_site_tides_create(pool, "IERS10_SITE", NULL, &tides, &diag) == MYB_E_FORMAT);
  MT_CHECK(strstr(diag.message, "IERS10_SITE_LOVE") != NULL);
  myb_kpool_destroy(pool);
  myb_site_tides_destroy(NULL);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  myb_site_tides* tides = model(NULL);
  test_solid(tides);
  test_pole(tides);
  test_ocean(tides);
  test_ocean_pole();
  test_atmosphere(tides);
  test_station(tides);
  test_kernels();
  myb_site_tides_destroy(tides);
  return mt_finish("test_site");
}
