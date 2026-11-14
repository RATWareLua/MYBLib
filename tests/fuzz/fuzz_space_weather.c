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
    if (myb_space_weather_day_at(weather, mjd, &day) == MYB_OK) {
      /* a value is missing (NaN) or inside its scale */
      for (int k = 0; k < 8; ++k) {
        MYB_CHECK(day.kp[k] != day.kp[k] || (day.kp[k] >= 0.0 && day.kp[k] <= 9.0));
        MYB_CHECK(day.ap[k] != day.ap[k] || (day.ap[k] >= 0.0 && day.ap[k] <= 400.0));
      }
      MYB_CHECK(day.f107_observed != day.f107_observed || day.f107_observed >= 0.0);
    }
    for (int i = 0; i < 4; ++i) {
      myb_space_weather_inputs inputs;
      if (myb_space_weather_inputs_at(weather, mjd, seconds[i], i & 1, &inputs) == MYB_OK) {
        const double given[3] = {inputs.f107, inputs.f107_average, inputs.ap};
        MYB_CHECK(fuzz_finite(given, 3) && fuzz_finite(inputs.ap_history, 7));
        MYB_CHECK(inputs.f107 >= 0.0 && inputs.ap >= 0.0);
      }
    }
  }
  myb_space_weather_destroy(weather);
  return 0;
}
