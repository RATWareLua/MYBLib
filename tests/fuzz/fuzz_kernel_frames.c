/* SPDX-License-Identifier: BSD-3-Clause */
/* Frames and body constants of text kernels: rotation models of the WGCCRE
 * (myb_wgccre_create, pck00011.tpc), fixed frames of frames kernels (myb_tkframe, a graph of
 * frames through myb_frame_graph_add_tkframe) and parameters (myb_param_from_kpool), for the
 * bodies and frames the input names and for a few fixed codes. */
#include "fuzz.h"
#include "myb/body.h"
#include "myb/frames.h"
#include "myb/io.h"

MYB_FUZZ_NO_INITIALIZE

static void use_body(const myb_kpool* pool, int64_t body, uint64_t hash) {
  myb_wgccre* model = NULL;
  myb_diag diag = {0};
  if (myb_wgccre_create(pool, body, NULL, &model, &diag) != MYB_OK) {
    return;
  }
  (void)myb_wgccre_frame(model);
  const double seconds[3] = {0.0, (double)(int64_t)(hash % 2000000000u) - 1.0e9, 3.0e9};
  for (int i = 0; i < 3; ++i) {
    double angles[3];
    double rates[3];
    myb_wgccre_euler(model, myb_time_from_seconds(seconds[i]), angles, rates);
    myb_rotation rotation;
    myb_wgccre_rotation(model, myb_time_from_seconds(seconds[i]), &rotation);
  }
  myb_wgccre_destroy(model);
}

static void use_frame(const myb_kpool* pool, int64_t id) {
  myb_rotation rotation;
  char relative[33];
  myb_diag diag = {0};
  if (myb_tkframe(pool, id, &rotation, relative, &diag) == MYB_OK) {
    MYB_CHECK(strlen(relative) <= 32);
  }
  if (id == 0 || id < INT32_MIN || id > INT32_MAX) {
    return;
  }
  myb_frame_graph* graph = NULL;
  MYB_CHECK(myb_frame_graph_create(NULL, NULL, &graph) == MYB_OK);
  if (myb_frame_graph_add_tkframe(graph, pool, (int32_t)id, &diag) == MYB_OK) {
    myb_rotation between;
    (void)myb_frame_graph_rotation(graph, (int32_t)id, 1, myb_time_from_seconds(0.0), &between);
  }
  myb_frame_graph_destroy(graph);
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_kpool* pool = NULL;
  MYB_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  myb_diag diag = {0};
  if (myb_kpool_load_text(pool, (const char*)data, size, &diag) == MYB_OK) {
    const uint64_t hash = fuzz_hash(data, size);
    const int64_t bodies[4] = {fuzz_id_after(data, size, "BODY", 399), 399, 301, 10};
    for (int i = 0; i < 4; ++i) {
      use_body(pool, bodies[i], hash);
    }
    const int64_t frames[4] = {fuzz_id_after(data, size, "TKFRAME_", 1400000), 1400000,
                               fuzz_id_after(data, size, "FRAME_", 17), 17};
    for (int i = 0; i < 4; ++i) {
      use_frame(pool, frames[i]);
    }
    myb_param param;
    (void)myb_param_from_kpool(pool, "BODY399_GM", NULL, 1.0e9, &param);
    (void)myb_param_from_kpool(pool, "MYB_CONST_G", "MYB_CONST_G_SIGMA", 1.0, &param);
  }
  myb_kpool_destroy(pool);
  return 0;
}
