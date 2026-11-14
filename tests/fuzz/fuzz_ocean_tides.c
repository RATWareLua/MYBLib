/* SPDX-License-Identifier: BSD-3-Clause */
/* Ocean tide coefficients of the geopotential (myb_tides_load_ocean, the FES form) over the
 * tides of data/grav/iers2010-tides.tk: the input loaded to a degree drawn from it, its waves
 * added once more (myb_tides_add_waves), then the tides applied to a small field. */
#include "fuzz.h"
#include "myb/grav.h"
#include "myb/io.h"

MYB_FUZZ_NO_INITIALIZE

static const char tide_kernel[] = {
#embed <grav/iers2010-tides.tk>
};

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_kpool* pool = NULL;
  MYB_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  myb_diag diag = {0};
  MYB_CHECK(myb_kpool_load_text(pool, tide_kernel, sizeof tide_kernel, &diag) == MYB_OK);
  myb_tides* tides = NULL;
  MYB_CHECK(myb_tides_create(pool, "IERS10_TIDE", NULL, &tides, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  const uint64_t hash = fuzz_hash(data, size);
  const int max_degree = 2 + (int)(hash % 12u);
  if (myb_tides_load_ocean(tides, (const char*)data, size, max_degree, &diag) == MYB_OK) {
    (void)myb_tides_ocean_waves(tides);
    (void)myb_tides_add_waves(tides, (const char*)data, size, max_degree, &diag);
    (void)myb_tides_secondary_waves(tides);
    enum { DEGREE = 8, COUNT = (DEGREE + 1) * (DEGREE + 2) / 2 };
    double c[COUNT] = {1.0};
    double s[COUNT] = {0.0};
    myb_gravity_field* field = NULL;
    MYB_CHECK(myb_gravity_field_create(3.986004418e14, 6378136.3, DEGREE, c, s, NULL, &field) ==
              MYB_OK);
    const myb_tide_raiser raisers[2] = {{.gm = 1.32712440041e20, .r = {1.4e11, 2.0e10, 1.0e10}},
                                        {.gm = 4.9028e12, .r = {3.0e8, 2.0e8, 1.0e8}}};
    const myb_tides_input input = {.raisers = raisers, .count = 2,
                                   .t = myb_time_from_seconds((double)(int64_t)(hash >> 35)),
                                   .flags = MYB_TIDES_OCEAN | MYB_TIDES_ADMITTANCE,
                                   .permanent = MYB_PERMANENT_TIDE_KEEP};
    (void)myb_tides_apply(tides, field, &input);
    myb_gravity_field_destroy(field);
  }
  myb_tides_destroy(tides);
  return 0;
}
