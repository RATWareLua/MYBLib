/* SPDX-License-Identifier: BSD-3-Clause */
/* GPT3 grids (ext/troposphere: myb_gpt3_load): the values at places over the globe and at the
 * poles, with and without the seasons, and the VMF3 mapping functions from them. */
#include "fuzz.h"
#include "myb_troposphere.h"

MYB_FUZZ_NO_INITIALIZE

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_gpt3* grid = NULL;
  myb_diag diag = {0};
  if (myb_gpt3_load((const char*)data, size, NULL, &grid, &diag) != MYB_OK) {
    MYB_CHECK(grid == NULL);
    return 0;
  }
  (void)myb_gpt3_step(grid);
  const double latitudes[5] = {-1.5707963267948966, -0.5, 0.0, 0.9, 1.5707963267948966};
  const double longitudes[5] = {-3.141592653589793, -1.0, 0.0, 2.5, 6.283185307179586};
  for (int i = 0; i < 5; ++i) {
    for (int constant = 0; constant < 2; ++constant) {
      myb_gpt3_values values;
      if (myb_gpt3_at(grid, 58849.25 + i, latitudes[i], longitudes[(i + constant) % 5],
                      120.0 * i, constant, &values) == MYB_OK) {
        double mh = 0.0;
        double mw = 0.0;
        (void)myb_tropo_vmf3(values.ah, values.aw, 58849.25, latitudes[i], longitudes[i],
                             120.0 * i, 1.2, &mh, &mw);
      }
    }
  }
  myb_gpt3_destroy(grid);
  return 0;
}
