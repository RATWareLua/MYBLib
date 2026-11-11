/* SPDX-License-Identifier: BSD-3-Clause */
/* The graph of frames (myb/frames.h): every frame is given in another by a provider of its
 * orientation; the rotation between any two is composed through their first common
 * ancestor. */
#include <inttypes.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/frames.h"

enum { MAX_DEPTH = MYB_FRAME_GRAPH_MAX_DEPTH, FIRST_CAPACITY = 16 };

typedef struct frame_entry {
  int32_t id;
  int32_t base;
  myb_orientation_function orientation; /* NULL: the constant rotation below */
  void* context;
  int has_scale; /* the provider's argument is on `scale`, not on the graph's */
  myb_scale scale;
  myb_rotation rotation;
} frame_entry;

struct myb_frame_graph {
  myb_allocator allocator;
  const myb_timesys* timesys;
  myb_scale own;
  frame_entry* frames;
  size_t count;
  size_t capacity;
};

myb_status myb_frame_graph_create(const myb_frame_graph_options* options,
                                  const myb_allocator* allocator, myb_frame_graph** graph) {
  MYB_CHECK(graph != NULL);
  const myb_allocator a = myb_internal_allocator(allocator);
  myb_scale own = {0};
  const myb_timesys* timesys = options != NULL ? options->timesys : NULL;
  if (timesys != NULL) {
    if (options->scale == NULL) {
      return MYB_E_DOMAIN;
    }
    const myb_status status = myb_timesys_scale(timesys, options->scale, &own);
    if (status != MYB_OK) {
      return status;
    }
  }
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_frame_graph), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_frame_graph* g = block;
  memset(g, 0, sizeof *g);
  g->allocator = a;
  g->timesys = timesys;
  g->own = own;
  *graph = g;
  return MYB_OK;
}

void myb_frame_graph_destroy(myb_frame_graph* graph) {
  if (graph == NULL) {
    return;
  }
  const myb_allocator a = graph->allocator;
  myb_internal_release(&a, graph->frames);
  myb_internal_release(&a, graph);
}

static const frame_entry* find(const myb_frame_graph* graph, int32_t id) {
  for (size_t i = 0; i < graph->count; ++i) {
    if (graph->frames[i].id == id) {
      return &graph->frames[i];
    }
  }
  return NULL;
}

/* A new entry for the frame, or NULL with the status of what went wrong. */
static myb_status new_entry(myb_frame_graph* graph, int32_t id, int32_t base,
                            frame_entry** entry) {
  if (id == 0 || base == 0 || id == base || find(graph, id) != NULL) {
    return MYB_E_DOMAIN;
  }
  if (graph->count == graph->capacity) {
    const size_t capacity = graph->capacity == 0 ? FIRST_CAPACITY : 2 * graph->capacity;
    void* block = NULL;
    const myb_status status =
        myb_internal_alloc_array(&graph->allocator, capacity, sizeof(frame_entry), &block);
    if (status != MYB_OK) {
      return status;
    }
    if (graph->count > 0) {
      memcpy(block, graph->frames, graph->count * sizeof(frame_entry));
    }
    myb_internal_release(&graph->allocator, graph->frames);
    graph->frames = block;
    graph->capacity = capacity;
  }
  frame_entry* e = &graph->frames[graph->count];
  memset(e, 0, sizeof *e);
  e->id = id;
  e->base = base;
  *entry = e;
  return MYB_OK;
}

myb_status myb_frame_graph_add(myb_frame_graph* graph, int32_t id, int32_t base,
                               myb_orientation_function orientation, void* context,
                               const char* scale) {
  MYB_CHECK(graph != NULL && orientation != NULL);
  myb_scale provider = {0};
  if (scale != NULL) {
    if (graph->timesys == NULL) {
      return MYB_E_DOMAIN;
    }
    const myb_status status = myb_timesys_scale(graph->timesys, scale, &provider);
    if (status != MYB_OK) {
      return status;
    }
  }
  frame_entry* e = NULL;
  const myb_status status = new_entry(graph, id, base, &e);
  if (status != MYB_OK) {
    return status;
  }
  e->orientation = orientation;
  e->context = context;
  e->has_scale = scale != NULL;
  e->scale = provider;
  ++graph->count;
  return MYB_OK;
}

myb_status myb_frame_graph_add_fixed(myb_frame_graph* graph, int32_t id, int32_t base,
                                     const myb_rotation* rotation) {
  MYB_CHECK(graph != NULL && rotation != NULL);
  frame_entry* e = NULL;
  const myb_status status = new_entry(graph, id, base, &e);
  if (status != MYB_OK) {
    return status;
  }
  e->rotation = *rotation;
  ++graph->count;
  return MYB_OK;
}

myb_status myb_frame_graph_add_tkframe(myb_frame_graph* graph, const myb_kpool* pool, int32_t id,
                                       myb_diag* diag) {
  MYB_CHECK(graph != NULL && pool != NULL);
  myb_internal_diag_clear(diag);
  myb_rotation rotation;
  char relative[33];
  myb_status status = myb_tkframe(pool, id, &rotation, relative, diag);
  if (status != MYB_OK) {
    return status;
  }
  char name[64] = "FRAME_";
  memcpy(name + 6, relative, strlen(relative) + 1);
  const double* code = NULL;
  size_t count = 0;
  status = myb_kpool_numbers(pool, name, &code, &count);
  if (status != MYB_OK || count != 1 || !(code[0] >= -2147483647.0 && code[0] <= 2147483647.0) ||
      code[0] != (double)(int32_t)code[0]) {
    myb_internal_diag(diag, 0, "%s: the code of the frame TKFRAME_%" PRId32 " is given in",
                      name, id);
    return status != MYB_OK ? status : MYB_E_FORMAT;
  }
  return myb_frame_graph_add_fixed(graph, id, (int32_t)code[0], &rotation);
}

/* The frames from `id` to the root of its tree: chain[0] = id. Returns how many, or 0 for a
 * chain longer than the limit (a cycle). */
static int chain_of(const myb_frame_graph* graph, int32_t id, int32_t* chain,
                    const frame_entry** entries) {
  chain[0] = id;
  for (int k = 0; k < MAX_DEPTH; ++k) {
    const frame_entry* e = find(graph, chain[k]);
    entries[k] = e;
    if (e == NULL) {
      return k + 1;
    }
    if (k + 1 < MAX_DEPTH) {
      chain[k + 1] = e->base;
    }
  }
  return 0;
}

/* The rotation base -> frame of an entry at t (on the scale of the graph). */
static myb_status entry_rotation(const myb_frame_graph* graph, const frame_entry* e, myb_time t,
                                 myb_rotation* rotation) {
  if (e->orientation == NULL) {
    *rotation = e->rotation;
    return MYB_OK;
  }
  myb_time at = t;
  if (e->has_scale) {
    const myb_status status = myb_timesys_convert(graph->timesys, t, graph->own, e->scale, &at);
    if (status != MYB_OK) {
      return status;
    }
  }
  myb_orientation o;
  const myb_status status = e->orientation(e->context, at, &o);
  if (status != MYB_OK) {
    return status;
  }
  *rotation = o.rotation;
  return MYB_OK;
}

myb_status myb_frame_graph_rotation(const myb_frame_graph* graph, int32_t from, int32_t to,
                                    myb_time t, myb_rotation* rotation) {
  MYB_CHECK(graph != NULL && rotation != NULL);
  memset(rotation, 0, sizeof *rotation);
  rotation->m[0] = rotation->m[4] = rotation->m[8] = 1.0;
  if (from == to) {
    return MYB_OK;
  }
  int32_t up[MAX_DEPTH];
  int32_t down[MAX_DEPTH];
  const frame_entry* up_entries[MAX_DEPTH];
  const frame_entry* down_entries[MAX_DEPTH];
  const int nu = chain_of(graph, from, up, up_entries);
  const int nd = chain_of(graph, to, down, down_entries);
  if (nu == 0 || nd == 0) {
    return MYB_E_FORMAT;
  }
  /* the first frame of the chain of `to` that the chain of `from` has too */
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
    return MYB_E_NOT_FOUND;
  }
  /* from -> common: the inverses of the frames of its chain, the nearest first */
  for (int k = 0; k < iu; ++k) {
    myb_rotation step;
    const myb_status status = entry_rotation(graph, up_entries[k], t, &step);
    if (status != MYB_OK) {
      return status;
    }
    myb_rotation_invert(&step, &step);
    myb_rotation_compose(&step, rotation, rotation);
  }
  /* common -> to: the frames of its chain, the farthest first */
  for (int k = id; k-- > 0;) {
    myb_rotation step;
    const myb_status status = entry_rotation(graph, down_entries[k], t, &step);
    if (status != MYB_OK) {
      return status;
    }
    myb_rotation_compose(&step, rotation, rotation);
  }
  return MYB_OK;
}

myb_status myb_orientation_graph(void* context, myb_time t, myb_orientation* orientation) {
  const myb_graph_orientation* c = context;
  MYB_CHECK(c != NULL && c->graph != NULL && orientation != NULL);
  const myb_status status =
      myb_frame_graph_rotation(c->graph, c->from, c->to, t, &orientation->rotation);
  orientation->ut = t;
  orientation->xp = 0.0;
  orientation->yp = 0.0;
  return status;
}
