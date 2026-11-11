/* SPDX-License-Identifier: BSD-3-Clause */
/* TDB - TT from the time ephemeris of DE440 (data/download/de440t.bsp; skipped without it):
 *  - what the file holds is told from the motion itself: the rate of TT - TDB at the
 *    geocentre must be L_C - (v^2 / 2 + U) / c^2, v the velocity of the Earth and U the
 *    potential of the other bodies there, taken from the same ephemeris;
 *  - against the analytic series of Fairhead and Bretagnon that the time system uses by
 *    default;
 *  - as the relation of TDB in a time system, and tabulated. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/io.h"
#include "myb/kin.h"
#include "myb/time.h"

#define DOWNLOAD MYB_DATA_DIR "/download/"

static double seconds_between(myb_time a, myb_time b) {
  return myb_dd_to_double(myb_time_diff(a, b));
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  myb_spk* spk = NULL;
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  if (myb_spk_load_file(spk, DOWNLOAD "de440t.bsp", NULL) != MYB_OK) {
    printf("de440t.bsp not downloaded: skipped\n");
    myb_spk_destroy(spk);
    return MTEST_SKIP;
  }
  myb_spk_time ephemeris = {.spk = spk};

  /* the rate of TT - TDB against the motion of the Earth. GM of DE440 (km^3/s^2): the
   * Sun, the planets' systems 1 .. 9 with the Earth and the Moon apart */
  static const struct {
    int32_t id;
    double gm;
  } bodies[] = {{10, 132712440041.279419}, {1, 22031.868551}, {2, 324858.592000},
                {301, 4902.800118},        {4, 42828.375816}, {5, 126712764.100000},
                {6, 37940584.841800},      {7, 5794556.400000}, {8, 6836527.100580},
                {9, 975.500000}};
  const double c = 299792.458;          /* km/s */
  const double lc = 1.48082686741e-8;   /* the mean rate of TCB - TCG (IAU 2000 B1.9) */
  double worst_rate = 0.0;
  double largest_rate = 0.0;
  for (int k = 0; k < 400; ++k) {
    const myb_time t = myb_time_from_seconds(-1.5e9 + 7.6e6 * k + 12345.0);
    double earth[6];
    double clock[6];
    MT_CHECK(myb_spk_state(spk, 399, 0, t, earth, NULL) == MYB_OK);
    MT_CHECK(myb_spk_state(spk, 1000000001, 1000000000, t, clock, NULL) == MYB_OK);
    double potential = 0.0;
    for (size_t b = 0; b < sizeof bodies / sizeof bodies[0]; ++b) {
      double other[6];
      MT_CHECK(myb_spk_state(spk, bodies[b].id, 399, t, other, NULL) == MYB_OK);
      potential += bodies[b].gm /
                   (sqrt(other[0] * other[0] + other[1] * other[1] + other[2] * other[2]) / 1e3);
    }
    const double v2 = (earth[3] * earth[3] + earth[4] * earth[4] + earth[5] * earth[5]) / 1e6;
    const double expected = lc - (0.5 * v2 + potential) / (c * c);
    const double rate = clock[3] / 1e3; /* s/s where a velocity in km/s would be */
    worst_rate = fmax(worst_rate, fabs(rate - expected));
    largest_rate = fmax(largest_rate, fabs(rate));
  }
  printf("DE440t: the rate of TT - TDB at the geocentre against L_C - (v^2 / 2 + U) / c^2 of "
         "the same ephemeris, 400 dates of 1952 - 2048: %.1e (the rate itself up to %.1e)\n",
         worst_rate, largest_rate);
  MT_CHECK(largest_rate > 3.0e-10 && largest_rate < 3.6e-10);
  MT_CHECK(worst_rate < 3e-16);

  /* against the series */
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/time-scales.tk", &diag) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/tdb-tt-fb1990.tk", &diag) == MYB_OK);
  myb_timesys* sys = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_scale tt;
  myb_scale tdb;
  MT_CHECK(myb_timesys_scale(sys, "TT", &tt) == MYB_OK);
  MT_CHECK(myb_timesys_scale(sys, "TDB", &tdb) == MYB_OK);
  double worst_series = 0.0;
  double mean = 0.0;
  for (int k = 0; k < 400; ++k) {
    const myb_time t = myb_time_from_seconds(-1.5e9 + 7.6e6 * k + 12345.0);
    myb_time on_tt;
    double from_file = 0.0;
    MT_CHECK(myb_timesys_convert(sys, t, tdb, tt, &on_tt) == MYB_OK);
    MT_CHECK(myb_spk_time_difference(&ephemeris, t, &from_file) == MYB_OK);
    const double apart = from_file - seconds_between(t, on_tt);
    worst_series = fmax(worst_series, fabs(apart));
    mean += apart / 400.0;
  }
  printf("  TDB - TT of the file against the series of Fairhead and Bretagnon: within %.1f ns "
         "(%.1f ns on the average)\n",
         1e9 * worst_series, 1e9 * mean);
  MT_CHECK(worst_series < 1e-7);

  /* as the relation of TDB: there and back; outside the file the status of the reader */
  MT_CHECK(myb_timesys_set_difference(sys, tdb, myb_spk_time_difference, &ephemeris) == MYB_OK);
  const myb_time start = myb_time_from_seconds(8.0e8 + 0.5);
  myb_time there;
  myb_time back;
  double value = 0.0;
  MT_CHECK(myb_timesys_convert(sys, start, tt, tdb, &there) == MYB_OK);
  MT_CHECK(myb_timesys_convert(sys, there, tdb, tt, &back) == MYB_OK);
  MT_CHECK(myb_spk_time_difference(&ephemeris, there, &value) == MYB_OK);
  MT_CHECK(fabs(seconds_between(there, start) - value) < 1e-15);
  MT_CHECK(fabs(seconds_between(back, start)) < 1e-15);
  MT_CHECK(myb_timesys_convert(sys, myb_time_from_seconds(-2.0e10), tt, tdb, &there) ==
           MYB_E_RANGE);
  /* tabulated */
  myb_time direct[50];
  for (int k = 0; k < 50; ++k) {
    MT_CHECK(myb_timesys_convert(sys, myb_time_add(start, 86400.0 * k + 777.0), tt, tdb,
                                 &direct[k]) == MYB_OK);
  }
  MT_CHECK(myb_timesys_tabulate(sys, tdb, myb_time_add(start, -10.0),
                                myb_time_add(start, 51.0 * 86400.0), 0.0) == MYB_OK);
  double worst_table = 0.0;
  for (int k = 0; k < 50; ++k) {
    MT_CHECK(myb_timesys_convert(sys, myb_time_add(start, 86400.0 * k + 777.0), tt, tdb,
                                 &there) == MYB_OK);
    worst_table = fmax(worst_table, fabs(seconds_between(there, direct[k])));
  }
  printf("  as the relation of TDB in the time system; at steps of an hour: %.1e s\n",
         worst_table);
  MT_CHECK(worst_table < 2e-15);
  myb_timesys_destroy(sys);
  myb_spk_destroy(spk);
  return mt_finish("test_de440_time");
}
