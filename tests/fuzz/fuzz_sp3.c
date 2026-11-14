/* SPDX-License-Identifier: BSD-3-Clause */
/* SP3-c and SP3-d precise orbits (myb_sp3_load): the satellites, epochs, frame and time
 * system, and of every record the state, the clock, the standard deviations and the flags. */
#include "fuzz.h"
#include "myb/kin.h"

MYB_FUZZ_NO_INITIALIZE

enum { MAX_RECORDS = 1 << 14 };

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_sp3* sp3 = NULL;
  myb_diag diag = {0};
  if (myb_sp3_load((const char*)data, size, NULL, &sp3, &diag) != MYB_OK) {
    MYB_CHECK(sp3 == NULL);
    return 0;
  }
  const size_t satellites = myb_sp3_satellite_count(sp3);
  const size_t epochs = myb_sp3_epoch_count(sp3);
  MYB_CHECK(myb_sp3_time_system(sp3) != NULL && myb_sp3_frame(sp3) != NULL);
  size_t records = 0;
  for (size_t s = 0; s < satellites; ++s) {
    MYB_CHECK(myb_sp3_satellite(sp3, s) != NULL);
    (void)myb_sp3_accuracy(sp3, s);
  }
  for (size_t e = 0; e < epochs && records < MAX_RECORDS; ++e) {
    myb_calendar label;
    myb_sp3_epoch(sp3, e, &label);
    for (size_t s = 0; s < satellites && records < MAX_RECORDS; ++s, ++records) {
      double position[3];
      double velocity[3];
      const int given = myb_sp3_state(sp3, s, e, position, velocity);
      MYB_CHECK(given >= 0 && given <= 3);
      double clock = 0.0;
      double rate = 0.0;
      const int clocks = myb_sp3_clock(sp3, s, e, &clock, &rate);
      MYB_CHECK(clocks >= 0 && clocks <= 3);
      double sigma[8];
      myb_sp3_sigma(sp3, s, e, sigma);
      (void)myb_sp3_flags(sp3, s, e);
    }
  }
  myb_sp3_destroy(sp3);
  return 0;
}
