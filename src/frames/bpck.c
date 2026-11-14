/* SPDX-License-Identifier: BSD-3-Clause */
#include <math.h>

#include "core/alloc.h"
#include "core/base.h"
#include "kin/segset.h"
#include "kin/spk_segment.h"
#include "myb/frames.h"

struct myb_pck {
  myb_internal_segset set;
};

myb_status myb_pck_create(const myb_allocator* allocator, myb_pck** pck) {
  MYB_CHECK(pck != NULL);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_pck), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_pck* p = block;
  myb_internal_segset_init(&p->set, &a, MYB_SEGSET_PCK);
  *pck = p;
  return MYB_OK;
}

void myb_pck_destroy(myb_pck* pck) {
  if (pck == NULL) {
    return;
  }
  const myb_allocator a = pck->set.allocator;
  myb_internal_segset_free(&pck->set);
  myb_internal_release(&a, pck);
}

myb_status myb_pck_load_file(myb_pck* pck, const char* path, myb_diag* diag) {
  MYB_CHECK(pck != NULL);
  return myb_internal_segset_load_file(&pck->set, path, diag);
}

myb_status myb_pck_load_memory(myb_pck* pck, const void* data, size_t size, myb_diag* diag) {
  MYB_CHECK(pck != NULL);
  return myb_internal_segset_load_memory(&pck->set, data, size, diag);
}

myb_status myb_pck_euler(const myb_pck* pck, int32_t id, myb_time t, double angles[3],
                         double rates[3], int32_t* base) {
  MYB_CHECK(pck != NULL && angles != NULL && rates != NULL);
  const myb_internal_spk_segment* s = myb_internal_segset_find(&pck->set, id, t);
  if (s == NULL) {
    return MYB_E_RANGE;
  }
  double state[6];
  const myb_status status = myb_internal_spk_segment_state(s, t, state);
  if (status != MYB_OK) {
    return status;
  }
  for (int i = 0; i < 3; ++i) {
    angles[i] = state[i];
    rates[i] = state[i + 3];
  }
  if (base != NULL) {
    *base = s->frame;
  }
  return MYB_OK;
}

myb_status myb_pck_rotation(const myb_pck* pck, int32_t id, myb_time t, myb_rotation* rotation,
                            int32_t* base) {
  MYB_CHECK(rotation != NULL);
  double angles[3];
  double rates[3];
  const myb_status status = myb_pck_euler(pck, id, t, angles, rates, base);
  if (status != MYB_OK) {
    return status;
  }
  myb_rotation_euler313(angles, rates, rotation);
  /* finite angles and rates, but the rates of a damaged file can be large enough for the sum
   * of their terms in the derivative of the matrix to leave the range of doubles */
  for (int i = 0; i < 9; ++i) {
    if (!isfinite(rotation->m[i]) || !isfinite(rotation->dm[i])) {
      return MYB_E_FORMAT;
    }
  }
  return MYB_OK;
}

myb_status myb_pck_rotation_many(const myb_pck* pck, int32_t id, const myb_time* t, size_t count,
                                 myb_rotation* rotations, int32_t* bases, size_t* failed) {
  MYB_CHECK(pck != NULL && (count == 0 || (t != NULL && rotations != NULL)));
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_pck_rotation(pck, id, t[k], &rotations[k], bases != NULL ? &bases[k] : NULL);
    if (status != MYB_OK) {
      break;
    }
  }
  if (failed != NULL) {
    *failed = k;
  }
  return status;
}
