/* SPDX-License-Identifier: BSD-3-Clause */
/* Two-line element sets (ext/sgp4: myb_tle_parse) on every pair of consecutive lines of the
 * input, with and without the checksums: an element set read is formatted (myb_tle_format),
 * which must read back to the same numbers, turned into a mean-elements message and
 * propagated by SGP4 in both modes. */
#include "fuzz.h"
#include "myb/ccsds.h"
#include "myb_sgp4.h"

MYB_FUZZ_NO_INITIALIZE

enum { MAX_LINES = 512 };

static int same(double a, double b) {
  return memcmp(&a, &b, sizeof a) == 0;
}

static void use_tle(const myb_tle* tle) {
  myb_calendar epoch;
  myb_tle_epoch(tle, &epoch);
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
  myb_omm omm;
  (void)myb_omm_from_tle(tle, &omm);
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
      (void)myb_sgp4_propagate(&sgp4, minutes[i], r, v, &error);
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
