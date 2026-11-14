/* SPDX-License-Identifier: BSD-3-Clause */
/* IONEX maps of the total electron content (ext/ionosphere: myb_ionex_load): the vertical
 * content at places and times around the epoch of the first map the header states, and for
 * three-dimensional maps the grid of densities (myb_ionex_density_grid) and a density from
 * it and from the layer over the maps. */
#include "fuzz.h"
#include "myb/time.h"
#include "myb_ionosphere.h"

MYB_FUZZ_NO_INITIALIZE

/* The MJD of the "EPOCH OF FIRST MAP" line (its first three numbers), or 60310 (2024). */
static double first_epoch(const uint8_t* data, size_t size) {
  char line[96];
  /* bound: at most 256 lines */
  for (size_t k = 0; k < 256 && fuzz_line(data, size, k, line, sizeof line); ++k) {
    if (strstr(line, "EPOCH OF FIRST MAP") == NULL) {
      continue;
    }
    long numbers[3] = {0, 0, 0};
    const char* p = line;
    for (int i = 0; i < 3; ++i) {
      /* bound: the line */
      while (*p == ' ') {
        ++p;
      }
      long v = 0;
      int digits = 0;
      /* bound: at most 9 digits */
      while (*p >= '0' && *p <= '9' && digits < 9) {
        v = 10 * v + (*p - '0');
        ++p;
        ++digits;
      }
      numbers[i] = v;
    }
    const myb_calendar c = {numbers[0], (int)numbers[1], (int)numbers[2], 0, 0, 0.0};
    myb_time t;
    if (myb_calendar_to_time(&c, &t) == MYB_OK) {
      return 51544.5 + myb_time_seconds(t) / 86400.0;
    }
  }
  return 60310.0;
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_ionex* maps = NULL;
  myb_diag diag = {0};
  if (myb_ionex_load((const char*)data, size, NULL, &maps, &diag) != MYB_OK) {
    MYB_CHECK(maps == NULL);
    return 0;
  }
  double radius = 0.0;
  double height = 0.0;
  myb_ionex_shell(maps, &radius, &height);
  const int dimension = myb_ionex_dimension(maps);
  MYB_CHECK(dimension == 2 || dimension == 3);
  const double mjd = first_epoch(data, size);
  const double latitudes[4] = {-1.5, -0.3, 0.7, 1.57};
  const double longitudes[4] = {-3.1, 0.0, 1.0, 3.14};
  for (int t = -1; t <= 9; ++t) {
    for (int i = 0; i < 4; ++i) {
      double vtec = 0.0;
      (void)myb_ionex_vtec(maps, mjd + t / 8.0, latitudes[i], longitudes[(i + t + 4) & 3],
                           &vtec);
    }
  }
  const double r[3] = {radius + height, 1.0e5, 2.0e5};
  myb_iono_layer layer = {.maps = maps, .peak = 350.0e3, .scale = 70.0e3};
  double density = 0.0;
  (void)myb_iono_layer_density(&layer, mjd + 0.1, r, &density);
  if (dimension == 3) {
    myb_iono_grid grid;
    if (myb_ionex_density_grid(maps, &grid) == MYB_OK) {
      (void)myb_iono_grid_density(&grid, mjd + 0.1, r, &density);
    }
  }
  myb_ionex_destroy(maps);
  return 0;
}
