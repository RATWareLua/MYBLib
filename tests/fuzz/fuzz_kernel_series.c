/* SPDX-License-Identifier: BSD-3-Clause */
/* The models given as Poisson series in kernels — the CIP of IAU 2006/2000A (myb_cip_create),
 * the equinox-based systems (myb_equinox_create), the diurnal terms of the EOP
 * (myb_eop_set_tides), the tides of the geopotential (myb_tides_create) and of a site
 * (myb_site_tides_create) — built from the input under the prefix it uses (the letters before
 * its first _ARGUMENT_COUNT) and under the shipped prefixes, then evaluated at one instant. */
#include "fuzz.h"
#include "myb/frames.h"
#include "myb/grav.h"
#include "myb/io.h"
#include "myb/site.h"

MYB_FUZZ_NO_INITIALIZE

/* The prefix of the first "<prefix>_ARGUMENT_COUNT" of the text (at most 20 characters). */
static void input_prefix(const uint8_t* data, size_t size, char prefix[24]) {
  static const char key[] = "_ARGUMENT_COUNT";
  const size_t n = sizeof key - 1;
  prefix[0] = '\0';
  for (size_t i = 1; i + n <= size; ++i) {
    if (memcmp(data + i, key, n) != 0) {
      continue;
    }
    size_t start = i;
    /* bound: at most 20 characters back */
    while (start > 0 && i - start < 20 &&
           ((data[start - 1] >= 'A' && data[start - 1] <= 'Z') ||
            (data[start - 1] >= '0' && data[start - 1] <= '9') || data[start - 1] == '_')) {
      --start;
    }
    memcpy(prefix, data + start, i - start);
    prefix[i - start] = '\0';
    return;
  }
}

static void use_models(const myb_kpool* pool, const char* prefix, myb_time t) {
  myb_diag diag = {0};
  myb_cip* cip = NULL;
  if (myb_cip_create(pool, prefix, NULL, &cip, &diag) == MYB_OK) {
    double xys[3];
    double rates[3];
    myb_cip_xys(cip, t, xys, rates);
    /* the series of a kernel that was read give numbers (bounded coefficients and scales) */
    MYB_CHECK(fuzz_finite(xys, 3) && fuzz_finite(rates, 3));
    (void)myb_cip_era(cip, t, NULL);
    (void)myb_cip_sprime(cip, t);
    myb_rotation rotation;
    myb_cip_terrestrial(cip, t, t, NULL, &rotation);
    myb_cip_destroy(cip);
  }
  myb_equinox* equinox = NULL;
  if (myb_equinox_create(pool, prefix, NULL, &equinox, &diag) == MYB_OK) {
    myb_rotation rotation;
    myb_equinox_rotation(equinox, MYB_EQUINOX_TEME, t, &rotation);
    myb_equinox_teme_pef(equinox, t, &rotation);
    myb_equinox_destroy(equinox);
  }
  myb_eop* eop = NULL;
  MYB_CHECK(myb_eop_create(NULL, &eop) == MYB_OK);
  (void)myb_eop_set_tides(eop, pool, prefix, &diag);
  myb_eop_destroy(eop);
  myb_tides* tides = NULL;
  if (myb_tides_create(pool, prefix, NULL, &tides, &diag) == MYB_OK) {
    (void)myb_tides_secondary_waves(tides);
    myb_tides_destroy(tides);
  }
  myb_site_tides* site = NULL;
  if (myb_site_tides_create(pool, prefix, NULL, &site, &diag) == MYB_OK) {
    (void)myb_site_tides_waves(site);
    (void)myb_site_tides_given_waves(site);
    const myb_tide_raiser raisers[2] = {{.gm = 332946.0482, .r = {1.4e11, 2.0e10, 1.0e10}},
                                        {.gm = 0.0123000371, .r = {3.0e8, 2.0e8, 1.0e8}}};
    const myb_site_input input = {.gm = 1.0, .raisers = raisers, .count = 2, .t = t,
                                  .dynamical = 69.184, .xp = 1e-6, .yp = 2e-6,
                                  .flags = MYB_SITE_SOLID | MYB_SITE_FREQUENCY | MYB_SITE_POLE};
    const double position[3] = {4075578.385, 931852.890, 4801570.154};
    double displacement[3];
    (void)myb_site_displacement(site, position, NULL, NULL, &input, displacement);
    myb_site_tides_destroy(site);
  }
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_kpool* pool = NULL;
  MYB_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  myb_diag diag = {0};
  if (myb_kpool_load_text(pool, (const char*)data, size, &diag) == MYB_OK) {
    const myb_time t = myb_time_from_seconds((double)(int64_t)(fuzz_hash(data, size) % 2000000000u));
    char prefix[24];
    input_prefix(data, size, prefix);
    if (prefix[0] != '\0') {
      use_models(pool, prefix, t);
    }
    static const char* const shipped[] = {"IAU2006A", "IAU1980", "IERS10_EOPT", "IERS10_TIDE",
                                          "IERS10_SITE"};
    for (size_t i = 0; i < sizeof shipped / sizeof shipped[0]; ++i) {
      if (strcmp(prefix, shipped[i]) != 0) {
        use_models(pool, shipped[i], t);
      }
    }
  }
  myb_kpool_destroy(pool);
  return 0;
}
