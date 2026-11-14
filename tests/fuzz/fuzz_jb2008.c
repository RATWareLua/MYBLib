/* SPDX-License-Identifier: BSD-3-Clause */
/* The index files of JB2008 (ext/jb2008): SOLFSMY (myb_jb2008_solar_load) and DTCFILE
 * (myb_jb2008_storm_load), both readers on every input, each table asked for the days it
 * holds and the days around them. */
#include "fuzz.h"
#include "myb_jb2008.h"

MYB_FUZZ_NO_INITIALIZE

/* The first day of a record "YYYY DDD" (after "DTC " for the storms), as an MJD. */
static double first_day(const uint8_t* data, size_t size, int storm) {
  char line[96];
  /* bound: at most 64 lines */
  for (size_t k = 0; k < 64 && fuzz_line(data, size, k, line, sizeof line); ++k) {
    const char* p = line;
    if (storm) {
      if (strncmp(p, "DTC", 3) != 0) {
        continue;
      }
      p += 3;
    } else if (line[0] == '#') {
      continue;
    }
    long numbers[2] = {0, 0};
    int ok = 1;
    for (int i = 0; i < 2 && ok; ++i) {
      /* bound: the line */
      while (*p == ' ') {
        ++p;
      }
      int digits = 0;
      /* bound: at most 6 digits */
      while (*p >= '0' && *p <= '9' && digits < 6) {
        numbers[i] = 10 * numbers[i] + (*p - '0');
        ++p;
        ++digits;
      }
      ok = digits > 0;
    }
    if (ok && numbers[0] > 1900 && numbers[0] < 2200) {
      const long y = numbers[0] - 1;
      /* the MJD of January 0 of the year, Gregorian */
      const long jan0 = 365 * y + y / 4 - y / 100 + y / 400 - 678576;
      return (double)(jan0 + numbers[1]);
    }
  }
  return 60000.0;
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_diag diag = {0};
  myb_jb2008_solar* solar = NULL;
  if (myb_jb2008_solar_load((const char*)data, size, NULL, &solar, &diag) == MYB_OK) {
    const double mjd = first_day(data, size, 0);
    for (int d = -6; d < 12; ++d) {
      myb_jb2008_indices indices;
      memset(&indices, 0, sizeof indices);
      (void)myb_jb2008_solar_at(solar, mjd + d + 0.5, &indices);
    }
    myb_jb2008_solar_destroy(solar);
  }
  myb_jb2008_storm* storm = NULL;
  if (myb_jb2008_storm_load((const char*)data, size, NULL, &storm, &diag) == MYB_OK) {
    const double mjd = first_day(data, size, 1);
    for (int h = -30; h < 80; ++h) {
      double dtc = 0.0;
      (void)myb_jb2008_storm_at(storm, mjd + h / 24.0, &dtc);
    }
    myb_jb2008_storm_destroy(storm);
  }
  return 0;
}
