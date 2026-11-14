/* SPDX-License-Identifier: BSD-3-Clause */
/* The time scales of a kernel (myb_timesys_create): the input is loaded over the definitions
 * of data/time/time-scales.tk and a series of TDB - TT (a UTC table like data/time/utc-tai.tk,
 * a series like data/time/tdb-tt-fb1990.tk, or definitions that replace the shipped ones), and
 * the system built from the pool converts an instant between its scales and to and from UTC
 * labels. */
#include "fuzz.h"
#include "myb/io.h"
#include "myb/time.h"

MYB_FUZZ_NO_INITIALIZE

static const char time_scales[] = {
#embed <time/time-scales.tk>
};

static void use_system(const myb_timesys* system, uint64_t hash) {
  static const char* const names[] = {"TAI", "TT", "TDB", "TCG", "TCB", "UTC", "GPS", "UT1"};
  enum { NAMES = sizeof names / sizeof names[0] };
  myb_scale scales[NAMES];
  int found[NAMES];
  for (int i = 0; i < NAMES; ++i) {
    found[i] = myb_timesys_scale(system, names[i], &scales[i]) == MYB_OK;
    if (found[i]) {
      MYB_CHECK(myb_timesys_scale_name(system, scales[i]) != NULL);
    }
  }
  const double seconds = (double)(int64_t)(hash % 4000000000u) - 2.0e9;
  const myb_time t = myb_time_from_seconds(seconds);
  for (int i = 0; i < NAMES; ++i) {
    for (int j = 0; j < NAMES; ++j) {
      if (found[i] && found[j]) {
        myb_time out;
        (void)myb_timesys_convert(system, t, scales[i], scales[j], &out);
      }
    }
    if (found[i]) {
      myb_calendar utc;
      if (myb_timesys_to_utc(system, t, scales[i], &utc) == MYB_OK) {
        myb_time back;
        (void)myb_timesys_from_utc(system, &utc, scales[i], &back);
      }
      const myb_calendar label = {1972 + (int64_t)(hash % 80u), 1 + (int)((hash >> 8) % 12u),
                                  1 + (int)((hash >> 16) % 28u), 23, 59, 60.5};
      myb_time out;
      (void)myb_timesys_from_utc(system, &label, scales[i], &out);
    }
  }
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_kpool* pool = NULL;
  MYB_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  myb_diag diag = {0};
  MYB_CHECK(myb_kpool_load_text(pool, time_scales, sizeof time_scales, &diag) == MYB_OK);
  MYB_CHECK(myb_kpool_load_text(pool, FUZZ_TDB_SERIES, sizeof FUZZ_TDB_SERIES - 1, &diag) ==
            MYB_OK);
  (void)myb_kpool_load_text(pool, (const char*)data, size, &diag);
  myb_timesys* system = NULL;
  if (myb_timesys_create(pool, NULL, &system, &diag) == MYB_OK) {
    use_system(system, fuzz_hash(data, size));
    myb_timesys_destroy(system);
  }
  myb_kpool_destroy(pool);
  return 0;
}
