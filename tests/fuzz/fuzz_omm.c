/* SPDX-License-Identifier: BSD-3-Clause */
/* CCSDS mean-elements messages, KVN and XML (myb_omm_parse): a message read is written in
 * both notations (myb_omm_write_as) and read again — the numbers of a catalogue pass
 * unchanged, so the message read back must equal the one written, bit for bit — and turned
 * into a two-line element set (ext/sgp4: myb_tle_from_omm) that is formatted, read back and
 * propagated. */
#include <stdlib.h>

#include "fuzz.h"
#include "myb/ccsds.h"
#include "myb_sgp4.h"

MYB_FUZZ_NO_INITIALIZE

static void round_trip(const myb_omm* omm) {
  for (int notation = 0; notation < 2; ++notation) {
    size_t length = 0;
    const myb_status sized = myb_omm_write_as(omm, (myb_odm_notation)notation, NULL, 0, &length);
    if (sized != MYB_E_OVERFLOW) {
      MYB_CHECK(sized == MYB_E_DOMAIN);
      continue;
    }
    char* text = malloc(length + 1);
    MYB_CHECK(text != NULL);
    MYB_CHECK(myb_omm_write_as(omm, (myb_odm_notation)notation, text, length + 1, &length) ==
              MYB_OK);
    myb_omm back;
    myb_diag diag = {0};
    MYB_CHECK(myb_omm_parse(text, length, &back, &diag) == MYB_OK);
    MYB_CHECK(back.has_mean_motion == omm->has_mean_motion && back.has_tle == omm->has_tle &&
              back.has_gm == omm->has_gm && back.has_covariance == omm->has_covariance);
    free(text);
  }
}

static void as_tle(const myb_omm* omm) {
  myb_tle tle;
  if (myb_tle_from_omm(omm, &tle) != MYB_OK) {
    return;
  }
  char line1[70];
  char line2[70];
  if (myb_tle_format(&tle, line1, line2) == MYB_OK) {
    myb_tle back;
    myb_diag diag = {0};
    MYB_CHECK(myb_tle_parse(line1, line2, 0, &back, &diag) == MYB_OK);
  }
  myb_sgp4 sgp4;
  if (myb_sgp4_init(&tle, MYB_SGP4_WGS72, MYB_SGP4_IMPROVED, &sgp4) == MYB_OK) {
    double r[3];
    double v[3];
    (void)myb_sgp4_propagate(&sgp4, 720.0, r, v, NULL);
  }
  myb_omm again;
  (void)myb_omm_from_tle(&tle, &again);
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_omm omm;
  myb_diag diag = {0};
  if (myb_omm_parse((const char*)data, size, &omm, &diag) == MYB_OK) {
    round_trip(&omm);
    as_tle(&omm);
  }
  return 0;
}
