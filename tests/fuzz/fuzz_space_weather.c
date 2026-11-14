/* SPDX-License-Identifier: BSD-3-Clause */
/* Space weather tables of the CelesTrak form (myb_space_weather_load): every day read (at
 * most 512) and the indices of the thermosphere models at three times of it, with the
 * observed and with the adjusted flux. */
#include "fuzz.h"
#include "myb/weather.h"

MYB_FUZZ_NO_INITIALIZE

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_space_weather* weather = NULL;
  myb_diag diag = {0};
  if (myb_space_weather_load((const char*)data, size, NULL, &weather, &diag) != MYB_OK) {
    MYB_CHECK(weather == NULL);
    return 0;
  }
  int64_t first = 0;
  int64_t last = 0;
  const size_t days = myb_space_weather_days(weather, &first, &last);
  MYB_CHECK(days > 0 && first <= last);
  const int64_t span = last - first;
  const int64_t step = span / 512 + 1;
  const double seconds[4] = {0.0, 43200.0, 86399.5, 86401.0};
  /* bound: at most 512 + 2 days */
  for (int64_t mjd = first - 1; mjd <= last + 1; mjd += step) {
    myb_space_weather_day day;
    (void)myb_space_weather_day_at(weather, mjd, &day);
    for (int i = 0; i < 4; ++i) {
      myb_space_weather_inputs inputs;
      (void)myb_space_weather_inputs_at(weather, mjd, seconds[i], i & 1, &inputs);
    }
  }
  myb_space_weather_destroy(weather);
  return 0;
}
