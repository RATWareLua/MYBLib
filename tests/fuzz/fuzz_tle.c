/* SPDX-License-Identifier: BSD-3-Clause */
/* Two-line element sets (ext/sgp4: myb_tle_parse) on every pair of consecutive lines of the
 * input, with and without the checksums: an element set read is formatted (myb_tle_format),
 * which must read back to the same numbers, turned into a mean-elements message and
 * propagated by SGP4 in both modes. */
#include <math.h>

#include "fuzz.h"
#include "myb/ccsds.h"
#include "myb_sgp4.h"

MYB_FUZZ_NO_INITIALIZE

enum { MAX_LINES = 512 };

static int same(double a, double b) {
  return memcmp(&a, &b, sizeof a) == 0;
}

static void use_tle(const myb_tle* tle) {
  /* the epoch a date of the calendar (a day 366 of a common year was 32 December) */
  static const int days[12] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  myb_calendar epoch;
  myb_tle_epoch(tle, &epoch);
  const int leap = (epoch.year % 4 == 0 && epoch.year % 100 != 0) || epoch.year % 400 == 0;
  MYB_CHECK(epoch.month >= 1 && epoch.month <= 12 && epoch.day >= 1 &&
            epoch.day <= days[epoch.month - 1] - (epoch.month == 2 && !leap ? 1 : 0) &&
            epoch.hour >= 0 && epoch.hour <= 23 && epoch.minute >= 0 && epoch.minute <= 59 &&
            epoch.second > -1e-6 && epoch.second < 60.0 + 1e-6);
  char line1[70];
  char line2[70];
  if (myb_tle_format(tle, line1, line2) == MYB_OK) {
    MYB_CHECK(strlen(line1) == 69 && strlen(line2) == 69);
    myb_tle back;
    myb_diag diag = {0};
    MYB_CHECK(myb_tle_parse(line1, line2, 0, &back, &diag) == MYB_OK);
    MYB_CHECK(same(back.epoch_day, tle->epoch_day) && same(back.eccentricity, tle->eccentricity) &&
              same(back.inclination, tle->inclination) && same(back.node, tle->node) &&
              same(back.mean_motion, tle->mean_motion) && same(back.bstar, tle->bstar) &&
              back.epoch_year == tle->epoch_year && back.revolution == tle->revolution);
  }
  /* the message of the set and the set of the message: the same epoch */
  myb_omm omm;
  if (myb_omm_from_tle(tle, &omm) == MYB_OK) {
    myb_tle again;
    if (myb_tle_from_omm(&omm, &again) == MYB_OK) {
      MYB_CHECK(again.epoch_year == tle->epoch_year &&
                fabs(again.epoch_day - tle->epoch_day) < 1e-8);
    }
  }
  for (int mode = 0; mode < 2; ++mode) {
    myb_sgp4 sgp4;
    if (myb_sgp4_init(tle, MYB_SGP4_WGS72, (myb_sgp4_mode)mode, &sgp4) != MYB_OK) {
      continue;
    }
    const double minutes[3] = {0.0, 1440.0, -4320.0};
    for (int i = 0; i < 3; ++i) {
      double r[3];
      double v[3];
      myb_sgp4_error error = MYB_SGP4_OK;
      if (myb_sgp4_propagate(&sgp4, minutes[i], r, v, &error) == MYB_OK) {
        for (int k = 0; k < 3; ++k) {
          MYB_CHECK(isfinite(r[k]) && isfinite(v[k])); /* a state of numbers, or an error */
        }
      }
    }
  }
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  char lines[2][96];
  if (!fuzz_line(data, size, 0, lines[0], sizeof lines[0])) {
    return 0;
  }
  /* bound: at most MAX_LINES lines */
  for (size_t k = 1; k < MAX_LINES && fuzz_line(data, size, k, lines[k & 1], sizeof lines[0]);
       ++k) {
    const char* first = lines[(k - 1) & 1];
    const char* second = lines[k & 1];
    for (unsigned flags = 0; flags < 2; ++flags) {
      myb_tle tle;
      myb_diag diag = {0};
      if (myb_tle_parse(first, second, flags ? MYB_TLE_NO_CHECKSUM : 0, &tle, &diag) == MYB_OK) {
        use_tle(&tle);
      }
    }
  }
  return 0;
}
