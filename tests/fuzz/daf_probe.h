/* SPDX-License-Identifier: BSD-3-Clause */
/* daf_probe.h — the summaries of a DAF image (SPK, binary PCK, CK) for the harnesses that
 * query what a reader loaded: the codes and the coverage of at most MAX_PROBES arrays, read
 * with the library's own DAF layer (io/daf.h). */
#ifndef MYB_TESTS_FUZZ_DAF_PROBE_H_
#define MYB_TESTS_FUZZ_DAF_PROBE_H_

#include <stdlib.h>

#include "fuzz.h"
#include "io/daf.h"

enum { MAX_PROBES = 32 };

typedef struct daf_probe {
  int count;
  int32_t ic[MAX_PROBES][6]; /* the first integers of each summary */
  double dc[MAX_PROBES][2];  /* the first doubles: the coverage */
} daf_probe;

/* Fills probe from the summaries of data[0 .. size); count 0 if it is not a DAF. */
static inline void daf_probe_read(const uint8_t* data, size_t size, daf_probe* probe) {
  probe->count = 0;
  myb_internal_daf daf;
  if (myb_internal_daf_open(data, size, &daf, NULL) != MYB_OK) {
    return;
  }
  myb_internal_daf_summary* summary = malloc(sizeof *summary);
  MYB_CHECK(summary != NULL);
  myb_internal_daf_cursor cursor = {0};
  /* bound: at most MAX_PROBES summaries */
  while (probe->count < MAX_PROBES) {
    int found = 0;
    if (myb_internal_daf_next(&daf, &cursor, summary, &found, NULL) != MYB_OK || !found) {
      break;
    }
    const int k = probe->count++;
    for (int i = 0; i < 6; ++i) {
      probe->ic[k][i] = i < daf.ni ? summary->ic[i] : 0;
    }
    for (int i = 0; i < 2; ++i) {
      probe->dc[k][i] = i < daf.nd ? summary->dc[i] : 0.0;
    }
  }
  free(summary);
}

#endif /* MYB_TESTS_FUZZ_DAF_PROBE_H_ */
