/* SPDX-License-Identifier: BSD-3-Clause */
/* EOP series: parsing of EOP 20 C04 and finals2000A excerpts (tests/data/eop), interpolation
 * with the diurnal/semidiurnal corrections against an independent evaluation the IERS way
 * (interp.f LAGINT, PMUT1_OCEANS, PMSDNUT2, UTLIBR from the original Fortran tables,
 * tests/data/eop/reference.txt), continuity of UT1 across the leap second of 2016-12-31,
 * errors. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/frames.h"
#include "myb/io.h"
#include "myb/time.h"

#define EOP_DIR MYB_TEST_DATA_DIR "/eop/"

static const double arcsec = 4.848136811095359935899141e-6;

static myb_timesys* load_system(void) {
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/time-scales.tk", NULL) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/tdb-tt-fb1990.tk", NULL) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/utc-tai.tk", NULL) == MYB_OK);
  myb_timesys* system = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &system, NULL) == MYB_OK);
  myb_kpool_destroy(pool);
  return system;
}

static myb_eop* load_eop(int finals, int tides) {
  myb_eop* eop = NULL;
  myb_diag diag;
  MT_CHECK(myb_eop_create(NULL, &eop) == MYB_OK);
  const myb_status status =
      finals ? myb_eop_load_finals_file(eop, EOP_DIR "finals2000A-excerpt.txt", &diag)
             : myb_eop_load_c04_file(eop, EOP_DIR "eopc04-excerpt.txt", &diag);
  if (status != MYB_OK) {
    fprintf(stderr, "EOP: %s line %ld: %s\n", myb_status_name(status), diag.line, diag.message);
  }
  MT_CHECK(status == MYB_OK);
  if (tides) {
    myb_kpool* pool = NULL;
    MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
    MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/frames/eop-tides.tk", &diag) == MYB_OK);
    MT_CHECK(myb_eop_set_tides(eop, pool, "IERS10_EOPT", &diag) == MYB_OK);
    myb_kpool_destroy(pool);
  }
  return eop;
}

/* TT instant of a UTC MJD. */
static myb_time tt_of_mjd(const myb_timesys* sys, myb_scale tt, double mjd) {
  myb_calendar c;
  myb_time_to_calendar(myb_time_from_seconds((mjd - 51544.5) * 86400.0), &c);
  myb_time t = {0};
  MT_CHECK(myb_timesys_from_utc(sys, &c, tt, &t) == MYB_OK);
  return t;
}

static void test_reference(const myb_timesys* sys, myb_scale tt) {
  myb_eop* eop = load_eop(0, 1);
  double first = 0.0;
  double last = 0.0;
  MT_CHECK(myb_eop_rows(eop, &first, &last) == 407 && first == 57740.0 && last == 59220.0);
  FILE* f = fopen(EOP_DIR "reference.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    myb_eop_destroy(eop);
    return;
  }
  char line[512];
  double worst_pm = 0.0;
  double worst_ut1 = 0.0;
  double worst_lod = 0.0;
  double worst_dxy = 0.0;
  int rows = 0;
  /* bound: one iteration per line of the file */
  while (fgets(line, sizeof line, f) != NULL) {
    if (line[0] == '#') {
      continue;
    }
    double r[7];
    if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &r[0], &r[1], &r[2], &r[3], &r[4], &r[5],
               &r[6]) != 7) {
      continue;
    }
    myb_eop_values v;
    myb_time ut1;
    const myb_time t = tt_of_mjd(sys, tt, r[0]);
    MT_CHECK(myb_eop_at(eop, sys, t, tt, &v, &ut1) == MYB_OK);
    worst_pm = fmax(worst_pm, fmax(fabs(v.xp / arcsec - r[1]), fabs(v.yp / arcsec - r[2])));
    worst_ut1 = fmax(worst_ut1, fabs(v.dut1 - r[3]));
    worst_lod = fmax(worst_lod, fabs(v.lod - r[4]));
    worst_dxy = fmax(worst_dxy, fmax(fabs(v.dx / arcsec - r[5]), fabs(v.dy / arcsec - r[6])));
    /* UT1 = UTC + (UT1 - UTC) */
    const double utc_seconds = (r[0] - 51544.5) * 86400.0;
    MT_CHECK(fabs(myb_time_seconds(ut1) - utc_seconds - v.dut1) < 1e-6);
    ++rows;
  }
  fclose(f);
  printf("vs IERS-way evaluation (%d dates): polar motion %.2e\", UT1-UTC %.2e s, LOD %.2e s, "
         "dX,dY %.2e\"\n",
         rows, worst_pm, worst_ut1, worst_lod, worst_dxy);
  MT_CHECK(rows == 10);
  /* The same formulas up to the order of operations and the reference's libm. UT1 is
   * interpolated as UT1 - TT (~ -69 s, one ulp 1.4e-14 s); the reference's LOD comes from
   * finite differences (1e-4 day) of its UT1 corrections, good to ~1e-10 s. */
  MT_CHECK(worst_pm < 5e-14 && worst_ut1 < 5e-14 && worst_lod < 3e-10 && worst_dxy < 1e-16);
  myb_eop_destroy(eop);
}

static void test_leap_second(const myb_timesys* sys, myb_scale tt) {
  myb_eop* eop = load_eop(0, 1);
  /* 2016-12-31T23:59:58 .. 2017-01-01T00:00:02 UTC: 5 SI seconds including 23:59:60 */
  const myb_calendar before = {2016, 12, 31, 23, 59, 58.0};
  myb_time t0;
  MT_CHECK(myb_timesys_from_utc(sys, &before, tt, &t0) == MYB_OK);
  myb_time previous_ut1 = {0};
  double worst = 0.0;
  for (int k = 0; k <= 10; ++k) {
    const myb_time t = myb_time_add(t0, 0.5 * k);
    myb_eop_values v;
    myb_time ut1;
    MT_CHECK(myb_eop_at(eop, sys, t, tt, &v, &ut1) == MYB_OK);
    if (k > 0) {
      /* UT1 advances by 0.5 s (1 - LOD/86400) whatever UTC does (the tidal terms change UT1
       * by < 1e-10 s in 0.5 s) */
      const double du = myb_time_diff(ut1, previous_ut1).hi;
      worst = fmax(worst, fabs(du - 0.5 * (1.0 - v.lod / 86400.0)));
    }
    previous_ut1 = ut1;
  }
  printf("UT1 steps across the leap second: max deviation %.2e s\n", worst);
  MT_CHECK(worst < 1e-9);
  myb_eop_values a;
  myb_eop_values b;
  const myb_calendar c1 = {2016, 12, 31, 23, 59, 59.0};
  const myb_calendar c2 = {2017, 1, 1, 0, 0, 1.0};
  myb_time t1;
  myb_time t2;
  MT_CHECK(myb_timesys_from_utc(sys, &c1, tt, &t1) == MYB_OK);
  MT_CHECK(myb_timesys_from_utc(sys, &c2, tt, &t2) == MYB_OK);
  MT_CHECK(myb_eop_at(eop, sys, t1, tt, &a, NULL) == MYB_OK);
  MT_CHECK(myb_eop_at(eop, sys, t2, tt, &b, NULL) == MYB_OK);
  printf("UT1-UTC: %.7f s at 23:59:59, %.7f s at 00:00:01 (leap second between)\n", a.dut1,
         b.dut1);
  MT_CHECK(fabs((b.dut1 - a.dut1) - 1.0) < 1e-4);
  myb_eop_destroy(eop);
}

static void test_finals(const myb_timesys* sys, myb_scale tt) {
  myb_eop* eop = load_eop(1, 0);
  MT_CHECK(myb_eop_rows(eop, NULL, NULL) == 31);
  /* at 0h of 2016-12-19 (MJD 57741) the row's own values come back exactly */
  myb_eop_values v;
  MT_CHECK(myb_eop_at(eop, sys, tt_of_mjd(sys, tt, 57741.0), tt, &v, NULL) == MYB_OK);
  MT_CHECK(v.xp == 0.108374 * arcsec && v.yp == 0.265923 * arcsec);
  MT_CHECK(fabs(v.dut1 - (-0.3952822)) < 1e-12);
  MT_CHECK(v.lod == 1.2693 * 1e-3);
  MT_CHECK(v.dx == -0.051 * 1e-3 * arcsec && v.dy == -0.287 * 1e-3 * arcsec);
  myb_eop_destroy(eop);
}

static void test_errors(const myb_timesys* sys, myb_scale tt) {
  myb_eop* eop = NULL;
  myb_diag diag;
  MT_CHECK(myb_eop_create(NULL, &eop) == MYB_OK);
  MT_CHECK(myb_eop_load_c04_file(eop, EOP_DIR "missing.txt", &diag) == MYB_E_IO);
  static const char bad[] =
      "# header\n"
      "2016  12  18   0  57740.00    0.109755    0.266098  -0.3940414   -0.000001   -0.000336"
      "   -0.001164   -0.000310   0.0012548\n"
      "2016  12  19   0  57745.00    0.108415    0.265923  -0.3952916   -0.000007   -0.000282"
      "   -0.001490   -0.000019   0.0012547\n";
  MT_CHECK(myb_eop_load_c04(eop, bad, strlen(bad), &diag) == MYB_E_FORMAT && diag.line == 3);
  MT_CHECK(myb_eop_rows(eop, NULL, NULL) == 0);
  myb_eop_values v;
  MT_CHECK(myb_eop_at(eop, sys, tt_of_mjd(sys, tt, 57750.0), tt, &v, NULL) == MYB_E_RANGE);
  myb_eop_destroy(eop);
  eop = load_eop(0, 0);
  /* two rows are needed on either side */
  MT_CHECK(myb_eop_at(eop, sys, tt_of_mjd(sys, tt, 57740.5), tt, &v, NULL) == MYB_E_RANGE);
  MT_CHECK(myb_eop_at(eop, sys, tt_of_mjd(sys, tt, 57741.5), tt, &v, NULL) == MYB_OK);
  MT_CHECK(myb_eop_at(eop, sys, tt_of_mjd(sys, tt, 59219.5), tt, &v, NULL) == MYB_E_RANGE);
  MT_CHECK(myb_eop_at(eop, sys, tt_of_mjd(sys, tt, 57760.5), tt, &v, NULL) == MYB_OK);
  /* in the gap of the excerpt between 2017-01-17 and 2019-12-27 the nodes are far apart but
   * still four: allowed (the caller supplies the table) */
  myb_eop_destroy(eop);
  myb_eop_destroy(NULL);
}

int main(void) {
  myb_timesys* sys = load_system();
  myb_scale tt = {0};
  MT_CHECK(myb_timesys_scale(sys, "TT", &tt) == MYB_OK);
  test_reference(sys, tt);
  test_leap_second(sys, tt);
  test_finals(sys, tt);
  test_errors(sys, tt);
  myb_timesys_destroy(sys);
  return mt_finish("test_eop");
}
