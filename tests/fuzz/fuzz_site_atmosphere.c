/* SPDX-License-Identifier: BSD-3-Clause */
/* The atmospheric tides S1 and S2 of a point (myb_site_atmosphere_load): the files of the
 * IERS loading calculator and the station files of TU Wien, for the names the input holds. */
#include <math.h>

#include "fuzz.h"
#include "myb/site.h"
#include "site_names.h"

MYB_FUZZ_NO_INITIALIZE

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  char names[SITE_NAMES][SITE_NAME_CHARS];
  const int count = site_names(data, size, names);
  for (int i = 0; i < count; ++i) {
    double coefficients[12];
    myb_diag diag = {0};
    if (myb_site_atmosphere_load((const char*)data, size, names[i], coefficients, &diag) ==
        MYB_OK) {
      for (int k = 0; k < 12; ++k) {
        MYB_CHECK(isfinite(coefficients[k]));
      }
    }
  }
  return 0;
}
