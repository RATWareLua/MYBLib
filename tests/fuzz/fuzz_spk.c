/* SPDX-License-Identifier: BSD-3-Clause */
/* SPK ephemerides (myb_spk_load_memory): the image loaded, every segment described and its
 * words read, the state of each target relative to its centre asked at the ends and the
 * middle of its coverage and at an instant drawn from the input; segments of type 10
 * (two-line elements) also through ext/sgp4 (myb_sgp4_spk10_state), with the time scales
 * and the equinox model of data/. */
#include <math.h>

#include "fuzz.h"
#include "myb/frames.h"
#include "myb/io.h"
#include "myb/kin.h"
#include "myb/time.h"
#include "myb_sgp4.h"

MYB_FUZZ_NO_INITIALIZE

static const char kernels[] = {
#embed <time/time-scales.tk>
    , '\n',
#embed <time/utc-tai.tk>
    , '\n',
#embed <frames/iau1980.tk>
};

enum { MAX_SEGMENTS = 64 };

/* The type-10 segments through the theory of ext/sgp4. */
static void use_type10(const myb_spk* spk, size_t index, const myb_spk_segment_info* info) {
  myb_kpool* pool = NULL;
  MYB_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  myb_diag diag = {0};
  MYB_CHECK(myb_kpool_load_text(pool, kernels, sizeof kernels, &diag) == MYB_OK);
  MYB_CHECK(myb_kpool_load_text(pool, FUZZ_TDB_SERIES, sizeof FUZZ_TDB_SERIES - 1, &diag) ==
            MYB_OK);
  myb_timesys* system = NULL;
  myb_equinox* equinox = NULL;
  MYB_CHECK(myb_timesys_create(pool, NULL, &system, &diag) == MYB_OK);
  MYB_CHECK(myb_equinox_create(pool, "IAU1980", NULL, &equinox, &diag) == MYB_OK);
  const double probes[3] = {info->start, 0.5 * (info->start + info->stop), info->stop};
  for (int i = 0; i < 3; ++i) {
    double state[6];
    myb_time t;
    if (fuzz_instant(probes[i], &t)) {
      (void)myb_sgp4_spk10_state(spk, index, system, equinox, t, state);
    }
  }
  myb_equinox_destroy(equinox);
  myb_timesys_destroy(system);
  myb_kpool_destroy(pool);
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_spk* spk = NULL;
  MYB_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  myb_diag diag = {0};
  if (myb_spk_load_memory(spk, data, size, &diag) == MYB_OK) {
    const uint64_t hash = fuzz_hash(data, size);
    const size_t count = myb_spk_segment_count(spk);
    int type10 = 0;
    for (size_t i = 0; i < count && i < MAX_SEGMENTS; ++i) {
      myb_spk_segment_info info;
      myb_spk_segment(spk, i, &info);
      MYB_CHECK(strlen(info.name) <= 40);
      const size_t words = myb_spk_segment_size(spk, i);
      double values[16];
      const size_t n = words < 16 ? words : 16;
      MYB_CHECK(myb_spk_segment_data(spk, i, 0, n, values) == MYB_OK);
      if (words > 0) {
        MYB_CHECK(myb_spk_segment_data(spk, i, words - 1, 1, values) == MYB_OK);
      }
      MYB_CHECK(myb_spk_segment_data(spk, i, words, 1, values) == MYB_E_RANGE);
      const double fraction = (double)(hash % 1024u) / 1023.0;
      const double probes[5] = {info.start, info.stop, 0.5 * (info.start + info.stop),
                                info.start + fraction * (info.stop - info.start),
                                nextafter(info.stop, (double)INFINITY)};
      for (int k = 0; k < 5; ++k) {
        double state[6];
        int32_t frame = 0;
        myb_time t;
        if (fuzz_instant(probes[k], &t)) {
          (void)myb_spk_state(spk, info.target, info.center, t, state, &frame);
          (void)myb_spk_state(spk, info.center, info.target, t, state, NULL);
        }
      }
      if (info.type == 10 && type10 < 2) {
        use_type10(spk, i, &info);
        ++type10;
      }
    }
    /* the same image a second time: precedence over the first */
    if (count < 8 && myb_spk_load_memory(spk, data, size, &diag) == MYB_OK && count > 0) {
      myb_spk_segment_info info;
      myb_spk_segment(spk, 0, &info);
      double state[6];
      myb_time t;
      if (fuzz_instant(info.start, &t)) {
        (void)myb_spk_state(spk, info.target, 0, t, state, NULL);
      }
    }
  }
  myb_spk_destroy(spk);
  return 0;
}
