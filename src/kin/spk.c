/* SPDX-License-Identifier: BSD-3-Clause */
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "kin/segset.h"
#include "kin/spk_segment.h"
#include "myb/kin.h"
#include "myb/num.h"

enum { MAX_CHAIN = 64 }; /* bodies between a body and the root of its tree (SPICE: 20) */

struct myb_spk {
  myb_internal_segset set;
};

myb_status myb_spk_create(const myb_allocator* allocator, myb_spk** spk) {
  MYB_CHECK(spk != NULL);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_spk), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_spk* s = block;
  myb_internal_segset_init(&s->set, &a, MYB_SEGSET_SPK);
  *spk = s;
  return MYB_OK;
}

void myb_spk_destroy(myb_spk* spk) {
  if (spk == NULL) {
    return;
  }
  const myb_allocator a = spk->set.allocator;
  myb_internal_segset_free(&spk->set);
  myb_internal_release(&a, spk);
}

myb_status myb_spk_load_memory(myb_spk* spk, const void* data, size_t size, myb_diag* diag) {
  MYB_CHECK(spk != NULL);
  return myb_internal_segset_load_memory(&spk->set, data, size, diag);
}

myb_status myb_spk_load_file(myb_spk* spk, const char* path, myb_diag* diag) {
  MYB_CHECK(spk != NULL);
  return myb_internal_segset_load_file(&spk->set, path, diag);
}

size_t myb_spk_segment_count(const myb_spk* spk) {
  MYB_CHECK(spk != NULL);
  return spk->set.count;
}

void myb_spk_segment(const myb_spk* spk, size_t index, myb_spk_segment_info* info) {
  MYB_CHECK(spk != NULL && info != NULL && index < spk->set.count);
  const myb_internal_spk_segment* s = &spk->set.segments[index];
  info->target = s->target;
  info->center = s->center;
  info->frame = s->frame;
  info->type = s->type;
  info->start = s->start;
  info->stop = s->stop;
  info->file = s->file;
  memcpy(info->name, s->name, sizeof info->name);
}

/* Bodies from `body` towards the root at t: chain[0] = body, segment k links chain[k] to
 * chain[k + 1]. Returns the number of bodies, or 0 if the chain is too long (a cycle). */
static int build_chain(const myb_spk* spk, int32_t body, myb_time t, int32_t* chain,
                       const myb_internal_spk_segment** links) {
  chain[0] = body;
  for (int k = 0; k < MAX_CHAIN - 1; ++k) {
    const myb_internal_spk_segment* s = myb_internal_segset_find(&spk->set, chain[k], t);
    if (s == NULL) {
      return k + 1;
    }
    links[k] = s;
    chain[k + 1] = s->center;
  }
  return 0;
}

myb_status myb_spk_state(const myb_spk* spk, int32_t target, int32_t center, myb_time t,
                         double state[6], int32_t* frame) {
  MYB_CHECK(spk != NULL && state != NULL);
  memset(state, 0, 6 * sizeof *state);
  if (frame != NULL) {
    *frame = 0;
  }
  if (target == center) {
    return MYB_OK;
  }
  int32_t up[MAX_CHAIN];
  int32_t down[MAX_CHAIN];
  const myb_internal_spk_segment* up_links[MAX_CHAIN];
  const myb_internal_spk_segment* down_links[MAX_CHAIN];
  const int nu = build_chain(spk, target, t, up, up_links);
  const int nd = build_chain(spk, center, t, down, down_links);
  if (nu == 0 || nd == 0) {
    return MYB_E_FORMAT;
  }
  /* The first body of the center's chain that is also in the target's chain. */
  int iu = -1;
  int id = -1;
  for (int j = 0; j < nd && iu < 0; ++j) {
    for (int i = 0; i < nu; ++i) {
      if (up[i] == down[j]) {
        iu = i;
        id = j;
        break;
      }
    }
  }
  if (iu < 0) {
    return MYB_E_RANGE;
  }
  myb_sum_acc sum[6] = {{0}};
  int32_t common_frame = 0;
  int have_frame = 0;
  for (int pass = 0; pass < 2; ++pass) {
    const int n = pass == 0 ? iu : id;
    const myb_internal_spk_segment* const* links = pass == 0 ? up_links : down_links;
    const double sign = pass == 0 ? 1.0 : -1.0;
    for (int k = 0; k < n; ++k) {
      if (have_frame && links[k]->frame != common_frame) {
        return MYB_E_FORMAT;
      }
      common_frame = links[k]->frame;
      have_frame = 1;
      double part[6];
      const myb_status status = myb_internal_spk_segment_state(links[k], t, part);
      if (status != MYB_OK) {
        return status;
      }
      for (int c = 0; c < 6; ++c) {
        myb_sum_acc_add(&sum[c], sign * part[c]);
      }
    }
  }
  for (int c = 0; c < 6; ++c) {
    state[c] = 1000.0 * myb_sum_acc_value(&sum[c]); /* km -> m */
  }
  if (frame != NULL) {
    *frame = common_frame;
  }
  return MYB_OK;
}
