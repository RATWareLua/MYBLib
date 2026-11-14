/* SPDX-License-Identifier: BSD-3-Clause */
/* Daily EOP series (myb_eop_load_c04 and myb_eop_load_finals on every input): the rows read
 * interpolated (myb_eop_at) at their ends, at their middle and at an instant drawn from the
 * input, on TT through the time scales and the UTC table of data/time. */
#include "fuzz.h"
#include "myb/frames.h"
#include "myb/io.h"
#include "myb/time.h"

MYB_FUZZ_NO_INITIALIZE

static const char time_kernels[] = {
#embed <time/time-scales.tk>
    , '\n',
#embed <time/utc-tai.tk>
};

static void use_table(const myb_eop* eop, const myb_timesys* system, myb_scale tt,
                      uint64_t hash) {
  double first = 0.0;
  double last = 0.0;
  if (myb_eop_rows(eop, &first, &last) == 0) {
    return;
  }
  const double fraction = (double)(hash % 1024u) / 1023.0;
  const double days[6] = {first, last, 0.5 * (first + last), first + fraction * (last - first),
                          first - 1.0, last + 0.75};
  for (int i = 0; i < 6; ++i) {
    myb_time t;
    if (!fuzz_instant((days[i] - 51544.5) * 86400.0, &t)) {
      continue;
    }
    myb_eop_values values;
    myb_time ut1;
    if (myb_eop_at(eop, system, t, tt, &values, &ut1) == MYB_OK) {
      /* what the rows can hold: angles below a radian, LOD below a day, UT1 within a day of
       * TT (and so of UTC, within a day and a minute) */
      const double all[6] = {values.xp, values.yp, values.dut1, values.lod, values.dx, values.dy};
      MYB_CHECK(fuzz_finite(all, 6));
      MYB_CHECK(values.xp > -1.0 && values.xp < 1.0 && values.yp > -1.0 && values.yp < 1.0);
      MYB_CHECK(values.dut1 > -2.0 * 86400.0 && values.dut1 < 2.0 * 86400.0);
    }
  }
}

/* The time system, built only for an input that a reader accepts. */
static myb_timesys* time_system(void) {
  myb_kpool* pool = NULL;
  MYB_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  myb_diag diag = {0};
  MYB_CHECK(myb_kpool_load_text(pool, time_kernels, sizeof time_kernels, &diag) == MYB_OK);
  MYB_CHECK(myb_kpool_load_text(pool, FUZZ_TDB_SERIES, sizeof FUZZ_TDB_SERIES - 1, &diag) ==
            MYB_OK);
  myb_timesys* system = NULL;
  MYB_CHECK(myb_timesys_create(pool, NULL, &system, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  return system;
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  const uint64_t hash = fuzz_hash(data, size);
  myb_diag diag = {0};
  myb_eop* tables[2] = {NULL, NULL};
  MYB_CHECK(myb_eop_create(NULL, &tables[0]) == MYB_OK);
  MYB_CHECK(myb_eop_create(NULL, &tables[1]) == MYB_OK);
  const int read[2] = {myb_eop_load_c04(tables[0], (const char*)data, size, &diag) == MYB_OK,
                       myb_eop_load_finals(tables[1], (const char*)data, size, &diag) == MYB_OK};
  if (read[0] || read[1]) {
    myb_timesys* system = time_system();
    myb_scale tt;
    MYB_CHECK(myb_timesys_scale(system, "TT", &tt) == MYB_OK);
    for (int k = 0; k < 2; ++k) {
      if (read[k]) {
        use_table(tables[k], system, tt, hash);
      }
    }
    myb_timesys_destroy(system);
  }
  myb_eop_destroy(tables[1]);
  myb_eop_destroy(tables[0]);
  return 0;
}
