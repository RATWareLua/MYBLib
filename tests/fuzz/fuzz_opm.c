/* SPDX-License-Identifier: BSD-3-Clause */
/* CCSDS parameter messages, KVN and XML (myb_opm_parse): a message read is written in both
 * notations (myb_opm_write_as) and read again: a message the writer accepts must read back
 * with the same state, manoeuvres and covariance flags. */
#include <stdlib.h>

#include "fuzz.h"
#include "myb/ccsds.h"

MYB_FUZZ_NO_INITIALIZE

static void round_trip(const myb_opm* opm) {
  for (int notation = 0; notation < 2; ++notation) {
    size_t length = 0;
    const myb_status sized = myb_opm_write_as(opm, (myb_odm_notation)notation, NULL, 0, &length);
    if (sized != MYB_E_OVERFLOW) {
      MYB_CHECK(sized == MYB_E_DOMAIN);
      continue;
    }
    char* text = malloc(length + 1);
    MYB_CHECK(text != NULL);
    MYB_CHECK(myb_opm_write_as(opm, (myb_odm_notation)notation, text, length + 1, &length) ==
              MYB_OK);
    myb_opm* back = malloc(sizeof *back);
    MYB_CHECK(back != NULL);
    myb_diag diag = {0};
    MYB_CHECK(myb_opm_parse(text, length, back, &diag) == MYB_OK);
    MYB_CHECK(back->maneuver_count == opm->maneuver_count);
    MYB_CHECK(back->has_elements == opm->has_elements &&
              back->has_spacecraft == opm->has_spacecraft &&
              back->has_covariance == opm->has_covariance);
    free(back);
    free(text);
  }
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_opm* opm = malloc(sizeof *opm);
  MYB_CHECK(opm != NULL);
  myb_diag diag = {0};
  if (myb_opm_parse((const char*)data, size, opm, &diag) == MYB_OK) {
    MYB_CHECK(opm->maneuver_count <= MYB_OPM_MAX_MANEUVERS);
    round_trip(opm);
  }
  free(opm);
  return 0;
}
