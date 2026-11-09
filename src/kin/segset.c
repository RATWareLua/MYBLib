/* SPDX-License-Identifier: BSD-3-Clause */
#include "kin/segset.h"

#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "io/daf.h"

void myb_internal_segset_init(myb_internal_segset* set, const myb_allocator* allocator,
                              myb_internal_segset_kind kind) {
  MYB_CHECK(set != NULL);
  memset(set, 0, sizeof *set);
  set->allocator = myb_internal_allocator(allocator);
  set->kind = kind;
}

void myb_internal_segset_free(myb_internal_segset* set) {
  if (set == NULL) {
    return;
  }
  for (size_t i = 0; i < set->file_count; ++i) {
    myb_platform_unmap_file(&set->files[i]);
  }
  myb_internal_release(&set->allocator, set->files);
  myb_internal_release(&set->allocator, set->segments);
  myb_internal_release(&set->allocator, set->order);
  set->files = NULL;
  set->segments = NULL;
  set->order = NULL;
  set->file_count = 0;
  set->count = 0;
}

/* order[] sorted by (id ascending, index descending): a stable bottom-up merge sort. */
static int before(const myb_internal_spk_segment* segments, size_t a, size_t b) {
  if (segments[a].target != segments[b].target) {
    return segments[a].target < segments[b].target;
  }
  return a > b;
}

static void sort_order(const myb_internal_spk_segment* segments, size_t* order, size_t* scratch,
                       size_t n) {
  /* bound: the run width doubles, at most 64 passes */
  for (size_t width = 1; width < n; width *= 2) {
    for (size_t lo = 0; lo < n; lo += 2 * width) {
      const size_t mid = lo + width < n ? lo + width : n;
      const size_t hi = lo + 2 * width < n ? lo + 2 * width : n;
      size_t i = lo;
      size_t j = mid;
      for (size_t k = lo; k < hi; ++k) {
        if (i < mid && (j >= hi || !before(segments, order[j], order[i]))) {
          scratch[k] = order[i++];
        } else {
          scratch[k] = order[j++];
        }
      }
    }
    memcpy(order, scratch, n * sizeof *order);
  }
}

/* Fills a segment from a summary; MYB_E_FORMAT for an invalid descriptor or layout. */
static myb_status describe(const myb_internal_segset* set, const myb_internal_daf* daf,
                           const myb_internal_daf_summary* summary, size_t file,
                           myb_internal_spk_segment* s, myb_diag* diag) {
  memset(s, 0, sizeof *s);
  s->daf = *daf;
  s->start = summary->dc[0];
  s->stop = summary->dc[1];
  s->file = file;
  int32_t begin;
  int32_t end;
  if (set->kind == MYB_SEGSET_SPK) {
    s->target = summary->ic[0];
    s->center = summary->ic[1];
    s->frame = summary->ic[2];
    s->type = summary->ic[3];
    begin = summary->ic[4];
    end = summary->ic[5];
  } else {
    s->target = summary->ic[0];
    s->frame = summary->ic[1];
    s->type = summary->ic[2];
    begin = summary->ic[3];
    end = summary->ic[4];
  }
  const size_t name_length = strlen(summary->name);
  memcpy(s->name, summary->name, name_length < 40 ? name_length : 40);
  if (!(s->start >= -0x1p62 && s->start <= s->stop && s->stop <= 0x1p62) || begin < 1 ||
      end < begin || (uint64_t)end > daf->words) {
    myb_internal_diag(diag, 0, "segment '%s': invalid descriptor", s->name);
    return MYB_E_FORMAT;
  }
  s->begin = (uint64_t)begin;
  s->end = (uint64_t)end;
  if (set->kind == MYB_SEGSET_PCK && s->type != 2 && s->type != 3 && s->type != 20) {
    s->supported = 0;
    return MYB_OK;
  }
  return myb_internal_spk_segment_parse(s, diag);
}

/* Reads the summaries of a DAF image into newly allocated segments. */
static myb_status read_segments(const myb_internal_segset* set, const unsigned char* data,
                                size_t size, myb_internal_spk_segment** out, size_t* out_count,
                                myb_diag* diag) {
  myb_internal_daf daf;
  myb_status status = myb_internal_daf_open(data, size, &daf, diag);
  if (status != MYB_OK) {
    return status;
  }
  const char* kind = set->kind == MYB_SEGSET_SPK ? "SPK" : "PCK";
  const int ni = set->kind == MYB_SEGSET_SPK ? 6 : 5;
  if ((daf.kind[0] != '\0' && strcmp(daf.kind, kind) != 0) || daf.nd != 2 || daf.ni != ni) {
    myb_internal_diag(diag, 0, "not a binary %s file (DAF/%s, ND = %d, NI = %d)", kind,
                      daf.kind, daf.nd, daf.ni);
    return MYB_E_FORMAT;
  }
  const myb_allocator* a = &set->allocator;
  myb_internal_spk_segment* segments = NULL;
  size_t count = 0;
  size_t capacity = 0;
  myb_internal_daf_cursor cursor = {0};
  myb_internal_daf_summary summary;
  /* bound: a summary occupies at least 5 words of the file, so at most words / 5 summaries */
  for (uint64_t guard = 0; guard <= daf.words / 5; ++guard) {
    int found = 0;
    status = myb_internal_daf_next(&daf, &cursor, &summary, &found, diag);
    if (status != MYB_OK || !found) {
      break;
    }
    if (count == capacity) {
      const size_t grown = capacity == 0 ? 16 : 2 * capacity;
      void* block = NULL;
      status = myb_internal_alloc_array(a, grown, sizeof *segments, &block);
      if (status != MYB_OK) {
        break;
      }
      if (count > 0) {
        memcpy(block, segments, count * sizeof *segments);
      }
      myb_internal_release(a, segments);
      segments = block;
      capacity = grown;
    }
    status = describe(set, &daf, &summary, set->file_count, &segments[count], diag);
    if (status != MYB_OK) {
      break;
    }
    ++count;
  }
  if (status != MYB_OK) {
    myb_internal_release(a, segments);
    return status;
  }
  *out = segments;
  *out_count = count;
  return MYB_OK;
}

static myb_status load(myb_internal_segset* set, const unsigned char* data, size_t size,
                       myb_platform_mapping mapping, myb_diag* diag) {
  myb_internal_diag_clear(diag);
  const myb_allocator* a = &set->allocator;
  myb_internal_spk_segment* added = NULL;
  size_t added_count = 0;
  myb_internal_spk_segment* segments = NULL;
  size_t* order = NULL;
  size_t* scratch = NULL;
  myb_platform_mapping* files = NULL;
  myb_status status = read_segments(set, data, size, &added, &added_count, diag);
  if (status != MYB_OK) {
    goto fail;
  }
  const size_t total = set->count + added_count;
  void* block = NULL;
  status = myb_internal_alloc_array(a, total, sizeof *segments, &block);
  if (status != MYB_OK) {
    goto fail;
  }
  segments = block;
  status = myb_internal_alloc_array(a, total, sizeof *order, &block);
  if (status != MYB_OK) {
    goto fail;
  }
  order = block;
  status = myb_internal_alloc_array(a, total, sizeof *scratch, &block);
  if (status != MYB_OK) {
    goto fail;
  }
  scratch = block;
  status = myb_internal_alloc_array(a, set->file_count + 1, sizeof *files, &block);
  if (status != MYB_OK) {
    goto fail;
  }
  files = block;
  if (set->count > 0) {
    memcpy(segments, set->segments, set->count * sizeof *segments);
  }
  if (added_count > 0) {
    memcpy(segments + set->count, added, added_count * sizeof *segments);
  }
  for (size_t i = 0; i < total; ++i) {
    order[i] = i;
  }
  sort_order(segments, order, scratch, total);
  if (set->file_count > 0) {
    memcpy(files, set->files, set->file_count * sizeof *files);
  }
  files[set->file_count] = mapping;
  myb_internal_release(a, set->segments);
  myb_internal_release(a, set->order);
  myb_internal_release(a, set->files);
  set->segments = segments;
  set->order = order;
  set->files = files;
  set->count = total;
  set->file_count += 1;
  myb_internal_release(a, scratch);
  myb_internal_release(a, added);
  return MYB_OK;
fail:
  myb_internal_release(a, files);
  myb_internal_release(a, scratch);
  myb_internal_release(a, order);
  myb_internal_release(a, segments);
  myb_internal_release(a, added);
  return status;
}

myb_status myb_internal_segset_load_memory(myb_internal_segset* set, const void* data,
                                           size_t size, myb_diag* diag) {
  MYB_CHECK(set != NULL && (data != NULL || size == 0));
  const myb_platform_mapping none = {0};
  return load(set, data, size, none, diag);
}

myb_status myb_internal_segset_load_file(myb_internal_segset* set, const char* path,
                                         myb_diag* diag) {
  MYB_CHECK(set != NULL && path != NULL);
  myb_internal_diag_clear(diag);
  myb_platform_mapping mapping = {0};
  myb_status status = myb_platform_map_file(path, &set->allocator, &mapping);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  status = load(set, mapping.data, mapping.size, mapping, diag);
  if (status != MYB_OK) {
    myb_platform_unmap_file(&mapping);
  }
  return status;
}

const myb_internal_spk_segment* myb_internal_segset_find(const myb_internal_segset* set,
                                                         int32_t id, myb_time t) {
  size_t lo = 0;
  size_t hi = set->count;
  /* bound: the interval halves */
  while (lo < hi) {
    const size_t mid = lo + (hi - lo) / 2;
    if (set->segments[set->order[mid]].target < id) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  for (size_t k = lo; k < set->count && set->segments[set->order[k]].target == id; ++k) {
    const myb_internal_spk_segment* s = &set->segments[set->order[k]];
    if (myb_time_compare(t, myb_time_from_seconds(s->start)) >= 0 &&
        myb_time_compare(t, myb_time_from_seconds(s->stop)) <= 0) {
      return s;
    }
  }
  return NULL;
}
