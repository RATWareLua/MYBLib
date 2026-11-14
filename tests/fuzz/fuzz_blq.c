/* SPDX-License-Identifier: BSD-3-Clause */
/* Ocean loading tables of the BLQ form (myb_ocean_loading_load) over the site tides of
 * data/frames/iers2010-site.tk: the records of the names the input holds, and the
 * displacement each gives at an instant drawn from the input. */
#include "fuzz.h"
#include "myb/io.h"
#include "myb/site.h"
#include "site_names.h"

MYB_FUZZ_NO_INITIALIZE

static const char site_kernel[] = {
#embed <frames/iers2010-site.tk>
};

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_kpool* pool = NULL;
  MYB_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  myb_diag diag = {0};
  MYB_CHECK(myb_kpool_load_text(pool, site_kernel, sizeof site_kernel, &diag) == MYB_OK);
  myb_site_tides* tides = NULL;
  MYB_CHECK(myb_site_tides_create(pool, "IERS10_SITE", NULL, &tides, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  char names[SITE_NAMES][SITE_NAME_CHARS];
  const int count = site_names(data, size, names);
  const myb_time t = myb_time_from_seconds((double)(int64_t)(fuzz_hash(data, size) % 1000000000u));
  for (int i = 0; i < count; ++i) {
    myb_ocean_loading* loading = NULL;
    if (myb_ocean_loading_load(tides, (const char*)data, size, names[i], NULL, &loading, &diag) ==
        MYB_OK) {
      double displacement[3];
      myb_ocean_loading_at(tides, loading, t, 69.184, displacement);
      myb_ocean_loading_destroy(loading);
    }
  }
  myb_site_tides_destroy(tides);
  return 0;
}
