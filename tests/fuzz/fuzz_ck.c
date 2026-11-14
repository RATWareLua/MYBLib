/* SPDX-License-Identifier: BSD-3-Clause */
/* CK orientation files (myb_ck_load_memory): every segment described, the pointing of its
 * instrument asked at the ends and the middle of its coverage and at an instant drawn from
 * the input, with and without a tolerance and angular velocities; what is given with MYB_OK
 * must be finite. */
#include <math.h>

#include "fuzz.h"
#include "myb/ck.h"

MYB_FUZZ_NO_INITIALIZE

enum { MAX_SEGMENTS = 64 };

/* A pointing given with MYB_OK is a rotation of finite numbers, whatever the file held. */
static void check_pointing(const myb_ck_pointing* p) {
  MYB_CHECK(isfinite(p->ticks));
  for (int i = 0; i < 4; ++i) {
    MYB_CHECK(isfinite(p->q[i]));
  }
  MYB_CHECK(fabs(p->q[0] * p->q[0] + p->q[1] * p->q[1] + p->q[2] * p->q[2] + p->q[3] * p->q[3] -
                 1.0) < 1e-12);
  for (int i = 0; i < 3; ++i) {
    MYB_CHECK(isfinite(p->rate[i]));
  }
  for (int i = 0; i < 9; ++i) {
    MYB_CHECK(isfinite(p->rotation.m[i]) && isfinite(p->rotation.dm[i]));
  }
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_ck* ck = NULL;
  MYB_CHECK(myb_ck_create(NULL, &ck) == MYB_OK);
  myb_diag diag = {0};
  if (myb_ck_load_memory(ck, data, size, &diag) == MYB_OK) {
    const double fraction = (double)(fuzz_hash(data, size) % 1024u) / 1023.0;
    const size_t count = myb_ck_segment_count(ck);
    for (size_t i = 0; i < count && i < MAX_SEGMENTS; ++i) {
      myb_ck_segment_info info;
      myb_ck_segment(ck, i, &info);
      MYB_CHECK(strlen(info.name) <= 40);
      const double span = info.stop - info.start;
      const double probes[5] = {info.start, info.stop, 0.5 * (info.start + info.stop),
                                info.start + fraction * span, nextafter(info.stop, (double)INFINITY)};
      const double tolerances[3] = {0.0, 1.0, isfinite(span) && span > 0.0 ? span : 1.0e9};
      for (int k = 0; k < 5; ++k) {
        for (int j = 0; j < 3; ++j) {
          myb_ck_pointing pointing;
          if (myb_ck_pointing_at(ck, info.instrument, probes[k], tolerances[j], j == 1,
                                 &pointing) == MYB_OK) {
            check_pointing(&pointing);
          }
        }
      }
    }
  }
  myb_ck_destroy(ck);
  return 0;
}
