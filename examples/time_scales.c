/* SPDX-License-Identifier: BSD-3-Clause */
/* One instant on several time scales.
 *
 * The scales are not built into the library: they come from the data kernels installed with it
 * (the CMake package gives their directory as myb_DATA_DIR). The program reads them, takes a
 * UTC date and prints the same instant on TAI, TT, TDB and GPS. */
#include <math.h>
#include <stdio.h>

#include "myb/myb.h"

#ifndef MYB_EXAMPLE_DATA_DIR
#error "define MYB_EXAMPLE_DATA_DIR: the directory of the data kernels"
#endif

int main(void) {
  static const char* const kernels[3] = {MYB_EXAMPLE_DATA_DIR "/time/time-scales.tk",
                                         MYB_EXAMPLE_DATA_DIR "/time/tdb-tt-fb1990.tk",
                                         MYB_EXAMPLE_DATA_DIR "/time/utc-tai.tk"};
  myb_kpool* pool = NULL;
  myb_diag diag;
  myb_status status = myb_kpool_create(NULL, &pool);
  for (int k = 0; k < 3 && status == MYB_OK; ++k) {
    status = myb_kpool_load_file(pool, kernels[k], &diag);
    if (status != MYB_OK) {
      fprintf(stderr, "%s: line %ld: %s\n", kernels[k], diag.line, diag.message);
    }
  }
  myb_timesys* system = NULL;
  if (status == MYB_OK) {
    status = myb_timesys_create(pool, NULL, &system, &diag);
  }
  myb_kpool_destroy(pool); /* the system keeps what it needs */
  if (status != MYB_OK) {
    fprintf(stderr, "time scales: %s\n", myb_status_name(status));
    return 1;
  }

  const myb_calendar utc = {2025, 6, 1, 0, 0, 0.0};
  static const char* const names[4] = {"TAI", "TT", "TDB", "GPS"};
  myb_time instant[4];
  printf("UTC %04lld-%02d-%02d %02d:%02d:%09.6f\n", (long long)utc.year, utc.month, utc.day,
         utc.hour, utc.minute, utc.second);
  for (int k = 0; k < 4; ++k) {
    myb_scale scale;
    status = myb_timesys_scale(system, names[k], &scale);
    if (status == MYB_OK) {
      status = myb_timesys_from_utc(system, &utc, scale, &instant[k]);
    }
    if (status != MYB_OK) {
      fprintf(stderr, "%s: %s\n", names[k], myb_status_name(status));
      return 1;
    }
    myb_calendar label;
    myb_time_to_calendar(instant[k], &label);
    printf("%-3s %04lld-%02d-%02d %02d:%02d:%09.6f\n", names[k], (long long)label.year,
           label.month, label.day, label.hour, label.minute, label.second);
  }
  const double tt_minus_tai = myb_dd_to_double(myb_time_diff(instant[1], instant[0]));
  const double tdb_minus_tt = myb_dd_to_double(myb_time_diff(instant[2], instant[1]));
  printf("TT - TAI = %.3f s (by definition 32.184), TDB - TT = %.6f ms\n", tt_minus_tai,
         1e3 * tdb_minus_tt);
  myb_timesys_destroy(system);
  return fabs(tt_minus_tai - 32.184) < 1e-12 && fabs(tdb_minus_tt) < 2e-3 ? 0 : 1;
}
