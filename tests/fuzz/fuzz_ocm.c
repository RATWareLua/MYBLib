/* SPDX-License-Identifier: BSD-3-Clause */
/* CCSDS comprehensive messages, KVN and XML (myb_ocm_load): the keywords of the header, the
 * metadata and the blocks of every kind, each line of data as text, as numbers and, for
 * trajectories, as a state in SI. */
#include "fuzz.h"
#include "myb/ccsds.h"

MYB_FUZZ_NO_INITIALIZE

enum { MAX_BLOCKS = 16, MAX_LINES = 256 };

static const char* const keys[] = {"CCSDS_OCM_VERS", "CREATION_DATE", "ORIGINATOR", "OBJECT_NAME",
                                   "TIME_SYSTEM", "EPOCH_TZERO", "TRAJ_TYPE", "TRAJ_UNITS",
                                   "CENTER_NAME", "COV_TYPE", "MAN_COMPOSITION", "TRAJ_REF_FRAME",
                                   "USER_DEFINED_CONSOLE_POC", "COMMENT"};
enum { KEYS = sizeof keys / sizeof keys[0] };

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_ocm* ocm = NULL;
  myb_diag diag = {0};
  if (myb_ocm_load((const char*)data, size, NULL, &ocm, &diag) != MYB_OK) {
    MYB_CHECK(ocm == NULL);
    return 0;
  }
  for (int i = 0; i < KEYS; ++i) {
    (void)myb_ocm_header(ocm, keys[i]);
    (void)myb_ocm_meta(ocm, keys[i]);
  }
  for (int kind = MYB_OCM_TRAJ; kind <= MYB_OCM_USER; ++kind) {
    const size_t blocks = myb_ocm_block_count(ocm, (myb_ocm_kind)kind);
    for (size_t b = 0; b < blocks && b < MAX_BLOCKS; ++b) {
      for (int i = 0; i < KEYS; ++i) {
        (void)myb_ocm_value(ocm, (myb_ocm_kind)kind, b, keys[i]);
      }
      const size_t lines = myb_ocm_line_count(ocm, (myb_ocm_kind)kind, b);
      for (size_t l = 0; l < lines && l < MAX_LINES; ++l) {
        MYB_CHECK(myb_ocm_line_text(ocm, (myb_ocm_kind)kind, b, l) != NULL);
        int absolute = 0;
        myb_calendar epoch;
        double offset = 0.0;
        double values[16];
        size_t count = 0;
        (void)myb_ocm_line(ocm, (myb_ocm_kind)kind, b, l, &absolute, &epoch, &offset, values, 16,
                           &count);
        (void)myb_ocm_line(ocm, (myb_ocm_kind)kind, b, l, &absolute, &epoch, &offset, NULL, 0,
                           &count);
        if (kind == MYB_OCM_TRAJ) {
          double state[9];
          int n = 0;
          if (myb_ocm_state(ocm, b, l, &absolute, &epoch, &offset, state, &n) == MYB_OK) {
            MYB_CHECK(n == 3 || n == 6 || n == 9);
          }
        }
      }
    }
  }
  myb_ocm_destroy(ocm);
  return 0;
}
