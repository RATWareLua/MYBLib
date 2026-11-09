/* SPDX-License-Identifier: BSD-3-Clause */
/* kin/segset.h — a set of DAF segments loaded from SPK or binary PCK files, indexed by the
 * identifier they describe and ordered by priority (later files first, later segments of a
 * file first: the SPICE rule).
 *
 * SPK summaries: start, stop; target, center, frame, type, begin, end (ND = 2, NI = 6).
 * Binary PCK summaries: start, stop; frame class ID, base frame, type, begin, end (ND = 2,
 * NI = 5); the PCK types 2, 3 and 20 have the layouts of the SPK types of the same numbers
 * (three Euler angles in place of the position, an angle scale in place of DSCALE). */
#ifndef MYB_SRC_KIN_SEGSET_H_
#define MYB_SRC_KIN_SEGSET_H_

#include <stddef.h>
#include <stdint.h>

#include "kin/spk_segment.h"
#include "myb/core.h"
#include "myb/time.h"
#include "platform/platform.h"

typedef enum myb_internal_segset_kind {
  MYB_SEGSET_SPK = 0,
  MYB_SEGSET_PCK = 1
} myb_internal_segset_kind;

typedef struct myb_internal_segset {
  myb_allocator allocator;
  myb_internal_segset_kind kind;
  myb_platform_mapping* files; /* zero entries for memory supplied by the caller */
  size_t file_count;
  myb_internal_spk_segment* segments; /* load order */
  size_t count;
  size_t* order; /* by identifier, then by decreasing priority */
} myb_internal_segset;

void myb_internal_segset_init(myb_internal_segset* set, const myb_allocator* allocator,
                              myb_internal_segset_kind kind);
void myb_internal_segset_free(myb_internal_segset* set);

/* Atomic: on error the set is unchanged. Errors: MYB_E_IO, MYB_E_FORMAT, MYB_E_NOMEM,
 * MYB_E_OVERFLOW. */
[[nodiscard]] myb_status myb_internal_segset_load_memory(myb_internal_segset* set,
                                                         const void* data, size_t size,
                                                         myb_diag* diag);
[[nodiscard]] myb_status myb_internal_segset_load_file(myb_internal_segset* set, const char* path,
                                                       myb_diag* diag);

/* The highest-priority segment for `id` (SPK target, PCK frame class) covering t, or NULL. */
const myb_internal_spk_segment* myb_internal_segset_find(const myb_internal_segset* set,
                                                         int32_t id, myb_time t);

#endif /* MYB_SRC_KIN_SEGSET_H_ */
