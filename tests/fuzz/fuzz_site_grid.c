/* SPDX-License-Identifier: BSD-3-Clause */
/* Grids of coefficients on the globe: the ocean pole tide load (myb_site_grid_load) and the
 * atmospheric tides of TU Wien (myb_site_atmosphere_grid_load), both readers on every input,
 * each grid read asked for its values at places inside, at its edges and outside. */
#include "fuzz.h"
#include "myb/site.h"

MYB_FUZZ_NO_INITIALIZE

static void use_grid(const myb_site_grid* grid) {
  const size_t values = myb_site_grid_values(grid);
  MYB_CHECK(values >= 1 && values <= MYB_SITE_GRID_MAX_VALUES);
  const double longitudes[6] = {-3.141592653589793, -0.001, 0.0, 1.0, 3.141592653589793, 6.3};
  const double latitudes[6] = {-1.5707963267948966, -1.56, 0.0, 0.5, 1.56, 1.5707963267948966};
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double out[MYB_SITE_GRID_MAX_VALUES];
      (void)myb_site_grid_at(grid, longitudes[i], latitudes[j], out);
    }
  }
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_diag diag = {0};
  myb_site_grid* grid = NULL;
  if (myb_site_grid_load((const char*)data, size, NULL, &grid, &diag) == MYB_OK) {
    use_grid(grid);
    myb_site_grid_destroy(grid);
  }
  grid = NULL;
  if (myb_site_atmosphere_grid_load((const char*)data, size, NULL, &grid, &diag) == MYB_OK) {
    use_grid(grid);
    myb_site_grid_destroy(grid);
  }
  return 0;
}
