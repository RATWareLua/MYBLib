/* SPDX-License-Identifier: BSD-3-Clause */
/* Spacecraft clock kernels (myb_sclk_create over a text kernel): the clock the input defines
 * converts ticks to times and back at the ends of its partitions, and clock strings — fixed
 * ones and each line of the input — to ticks and back. */
#include <math.h>

#include "fuzz.h"
#include "myb/ck.h"
#include "myb/io.h"

MYB_FUZZ_NO_INITIALIZE

static void use_ticks(const myb_sclk* sclk, double ticks) {
  myb_time t;
  if (myb_sclk_to_time(sclk, ticks, &t) == MYB_OK) {
    double back = 0.0;
    (void)myb_sclk_from_time(sclk, t, &back);
  }
  char text[128];
  size_t length = 0;
  if (myb_sclk_decode(sclk, ticks, text, sizeof text, &length) == MYB_OK) {
    MYB_CHECK(length == strlen(text));
    double back = 0.0;
    (void)myb_sclk_encode(sclk, text, length, &back);
  }
  char tiny[4];
  (void)myb_sclk_decode(sclk, ticks, tiny, sizeof tiny, &length);
}

static void use_clock(const myb_sclk* sclk, const uint8_t* data, size_t size) {
  (void)myb_sclk_time_system(sclk);
  (void)myb_sclk_partitions(sclk);
  static const char* const strings[] = {"1/0", "0", "1/1:2:3", "2/99999999.255.7", "1/-1",
                                        "3/1", "", "  1 / 12 ", "1/1.2.3.4.5.6.7.8.9.10.11",
                                        "1/18446744073709551616", "/", "1/:"};
  for (size_t i = 0; i < sizeof strings / sizeof strings[0]; ++i) {
    double ticks = 0.0;
    if (myb_sclk_encode(sclk, strings[i], strlen(strings[i]), &ticks) == MYB_OK) {
      use_ticks(sclk, ticks);
    }
  }
  char line[96];
  /* bound: at most 64 lines */
  for (size_t k = 0; k < 64 && fuzz_line(data, size, k, line, sizeof line); ++k) {
    double ticks = 0.0;
    if (myb_sclk_encode(sclk, line, strlen(line), &ticks) == MYB_OK) {
      use_ticks(sclk, ticks);
    }
  }
  const double probes[] = {-1.0, 0.0, 0.5, 1.0e6, 4.0e9, 1.0e15, 1.0e300, (double)INFINITY, (double)NAN};
  for (size_t i = 0; i < sizeof probes / sizeof probes[0]; ++i) {
    use_ticks(sclk, probes[i]);
    myb_time t;
    double ticks = 0.0;
    if (fuzz_instant(probes[i] * 1.0e-3, &t)) {
      (void)myb_sclk_from_time(sclk, t, &ticks);
    }
  }
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_kpool* pool = NULL;
  MYB_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  myb_diag diag = {0};
  if (myb_kpool_load_text(pool, (const char*)data, size, &diag) == MYB_OK) {
    const int64_t ids[3] = {fuzz_id_after(data, size, "SCLK_DATA_TYPE_", 82),
                            fuzz_id_after(data, size, "SCLK01_COEFFICIENTS_", 82), 9099};
    for (int i = 0; i < 3; ++i) {
      myb_sclk* sclk = NULL;
      if (myb_sclk_create(pool, ids[i], NULL, &sclk, &diag) == MYB_OK) {
        use_clock(sclk, data, size);
        myb_sclk_destroy(sclk);
      }
    }
  }
  myb_kpool_destroy(pool);
  return 0;
}
