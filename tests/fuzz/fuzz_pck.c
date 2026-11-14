/* SPDX-License-Identifier: BSD-3-Clause */
/* Binary PCK files (myb_pck_load_memory): the Euler angles and the rotation of every frame
 * the summaries of the image name, at the ends and the middle of each segment's coverage. */
#include <math.h>

#include "daf_probe.h"
#include "myb/frames.h"

MYB_FUZZ_NO_INITIALIZE

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_pck* pck = NULL;
  MYB_CHECK(myb_pck_create(NULL, &pck) == MYB_OK);
  myb_diag diag = {0};
  if (myb_pck_load_memory(pck, data, size, &diag) == MYB_OK) {
    daf_probe probe;
    daf_probe_read(data, size, &probe);
    const double fraction = (double)(fuzz_hash(data, size) % 1024u) / 1023.0;
    for (int k = 0; k < probe.count; ++k) {
      const double start = probe.dc[k][0];
      const double stop = probe.dc[k][1];
      const double probes[5] = {start, stop, 0.5 * (start + stop),
                                start + fraction * (stop - start), nextafter(start, -(double)INFINITY)};
      for (int i = 0; i < 5; ++i) {
        double angles[3];
        double rates[3];
        int32_t base = 0;
        myb_time t;
        if (!fuzz_instant(probes[i], &t)) {
          continue;
        }
        (void)myb_pck_euler(pck, probe.ic[k][0], t, angles, rates, &base);
        myb_rotation rotation;
        (void)myb_pck_rotation(pck, probe.ic[k][0], t, &rotation, NULL);
      }
    }
  }
  myb_pck_destroy(pck);
  return 0;
}
