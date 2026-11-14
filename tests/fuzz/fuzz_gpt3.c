/* SPDX-License-Identifier: BSD-3-Clause */
/* GPT3 grids (ext/troposphere: myb_gpt3_load): the values at places over the globe and at the
 * poles, with and without the seasons, and the VMF3 mapping functions from them; longitudes of
 * many turns and heights far from the nodes, which must give numbers or an error. */
#include <math.h>

#include "fuzz.h"
#include "myb_troposphere.h"

MYB_FUZZ_NO_INITIALIZE

static void values_are_numbers(const myb_gpt3_values* v) {
  const double all[13] = {v->pressure, v->temperature, v->lapse_rate, v->mean_temperature,
                          v->vapour_pressure, v->ah, v->aw, v->lambda, v->undulation,
                          v->gn_h, v->ge_h, v->gn_w, v->ge_w};
  for (int k = 0; k < 13; ++k) {
    MYB_CHECK(isfinite(all[k]));
  }
}

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
        values_are_numbers(&values);
        double mh = 0.0;
        double mw = 0.0;
        (void)myb_tropo_vmf3(values.ah, values.aw, 58849.25, latitudes[i], longitudes[i],
                             120.0 * i, 1.2, &mh, &mw);
      }
    }
  }
  const double far_longitudes[4] = {-1.0e18, 1.0e18, 7.0e15, 1.0e307};
  const double far_heights[4] = {-1.0e7, 1.0e7, -5.0e4, 9.0e3};
  for (int i = 0; i < 4; ++i) {
    myb_gpt3_values values;
    if (myb_gpt3_at(grid, 58849.25, 0.3, far_longitudes[i], 100.0, 0, &values) == MYB_OK) {
      values_are_numbers(&values);
    }
    if (myb_gpt3_at(grid, 58849.25, -0.7, 1.0, far_heights[i], 0, &values) == MYB_OK) {
      values_are_numbers(&values);
    }
  }
  myb_gpt3_destroy(grid);
  return 0;
}
