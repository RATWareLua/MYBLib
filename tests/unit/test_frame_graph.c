/* SPDX-License-Identifier: BSD-3-Clause */
/* The graph of frames and SPK chains through different frames, against CSPICE N0067
 * (tests/data/frames/graph_cspice.ref: spkgeo_c and sxform_c on an SPK whose segments are in
 * ECLIPJ2000, B1950, GALACTIC, FK4, J2000 and a turning frame of a binary PCK); the time
 * scales of the providers; the errors. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/frames.h"
#include "myb/io.h"
#include "myb/kin.h"

#define FRAMES_DIR MYB_TEST_DATA_DIR "/frames/"

enum { TURNING = 1400099, TURNING_CLASS = 31099 };

static int next_number(const char** p, double* value) {
  /* bound: skips at most the rest of the line */
  while (**p == ' ' || **p == '\t') {
    ++*p;
  }
  size_t used = 0;
  if (myb_parse_double(*p, strlen(*p), value, &used) != MYB_OK) {
    return 0;
  }
  *p += used;
  return 1;
}

static double norm3(const double* v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/* The inertial frames of data/frames/naif-inertial.tk and the turning frame of the PCK. */
static myb_frame_graph* build(const myb_frame_graph_options* options, myb_pck_orientation* turning) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/frames/naif-inertial.tk", &diag) == MYB_OK);
  myb_frame_graph* graph = NULL;
  MT_CHECK(myb_frame_graph_create(options, NULL, &graph) == MYB_OK);
  const int32_t fixed[5] = {2, 3, 13, 17, 18};
  for (int i = 0; i < 5; ++i) {
    MT_CHECK(myb_frame_graph_add_tkframe(graph, pool, fixed[i], &diag) == MYB_OK);
  }
  /* defined already; not in the kernel */
  MT_CHECK(myb_frame_graph_add_tkframe(graph, pool, 17, &diag) == MYB_E_DOMAIN);
  MT_CHECK(myb_frame_graph_add_tkframe(graph, pool, 99, &diag) == MYB_E_NOT_FOUND);
  myb_kpool_destroy(pool);
  MT_CHECK(myb_frame_graph_add(graph, TURNING, 1, myb_orientation_pck, turning, NULL) == MYB_OK);
  return graph;
}

static void test_reference(void) {
  myb_pck* pck = NULL;
  myb_spk* spk = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_pck_create(NULL, &pck) == MYB_OK);
  MT_CHECK(myb_pck_load_file(pck, FRAMES_DIR "graph.bpc", &diag) == MYB_OK);
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  MT_CHECK(myb_spk_load_file(spk, FRAMES_DIR "graph.bsp", &diag) == MYB_OK);
  myb_pck_orientation turning = {pck, TURNING_CLASS, NULL};
  myb_frame_graph* graph = build(NULL, &turning);
  FILE* f = fopen(FRAMES_DIR "graph_cspice.ref", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  char line[2048];
  int states = 0;
  int rotations = 0;
  double worst_r = 0.0;
  double worst_v = 0.0;
  double worst_m = 0.0;
  double worst_dm = 0.0;
  /* bound: one iteration per line of the file */
  while (fgets(line, sizeof line, f) != NULL) {
    const char* p = line + 1;
    double v[24];
    if (line[0] == 'S') {
      int ok = 1;
      for (int i = 0; i < 10; ++i) {
        ok = ok && next_number(&p, &v[i]);
      }
      MT_CHECK(ok);
      double state[6];
      MT_CHECK(myb_spk_state_in(spk, graph, (int32_t)v[0], (int32_t)v[1],
                                myb_time_from_seconds(v[3]), (int32_t)v[2], state) == MYB_OK);
      double dr[3];
      double dv[3];
      for (int i = 0; i < 3; ++i) {
        dr[i] = state[i] / 1000.0 - v[4 + i];
        dv[i] = state[3 + i] / 1000.0 - v[7 + i];
      }
      worst_r = fmax(worst_r, norm3(dr) / norm3(v + 4));
      worst_v = fmax(worst_v, norm3(dv) / norm3(v + 7));
      ++states;
    } else if (line[0] == 'R') {
      int ok = 1;
      for (int i = 0; i < 21; ++i) {
        ok = ok && next_number(&p, &v[i]);
      }
      MT_CHECK(ok);
      myb_rotation rotation;
      MT_CHECK(myb_frame_graph_rotation(graph, (int32_t)v[0], (int32_t)v[1],
                                        myb_time_from_seconds(v[2]), &rotation) == MYB_OK);
      for (int i = 0; i < 9; ++i) {
        worst_m = fmax(worst_m, fabs(rotation.m[i] - v[3 + i]));
        worst_dm = fmax(worst_dm, fabs(rotation.dm[i] - v[12 + i]) / 7.2921e-5);
      }
      ++rotations;
    }
  }
  fclose(f);
  printf("frame graph against CSPICE: %d states through segments in six frames, r %.1e, v "
         "%.1e; %d rotations between seven frames, matrix %.1e, its derivative %.1e of the "
         "turning rate\n",
         states, worst_r, worst_v, rotations, worst_m, worst_dm);
  MT_CHECK(states == 336 && rotations == 126);
  MT_CHECK(worst_r < 5e-15 && worst_v < 5e-14);
  MT_CHECK(worst_m < 1e-15 && worst_dm < 1e-14);

  /* the set alone refuses segments in different frames; as a provider with the graph it
   * gives the same state */
  const myb_time t = myb_time_from_seconds(6.0e8 + 1000.0);
  double direct[6];
  double provided[6];
  int32_t frame = 0;
  MT_CHECK(myb_spk_state(spk, 1602, 399, t, direct, &frame) == MYB_E_FORMAT);
  MT_CHECK(myb_spk_state(spk, 1601, 399, t, direct, &frame) == MYB_OK && frame == 17);
  MT_CHECK(myb_spk_state_in(spk, graph, 1602, 399, t, 13, direct) == MYB_OK);
  myb_spk_body body = {.spk = spk, .target = 1602, .center = 399, .frame = 13, .graph = graph};
  MT_CHECK(myb_spk_body_state(&body, t, provided) == MYB_OK);
  MT_CHECK(memcmp(direct, provided, sizeof direct) == 0);
  /* a frame the graph does not know */
  MT_CHECK(myb_spk_state_in(spk, graph, 1602, 399, t, 777, direct) == MYB_E_NOT_FOUND);
  /* the graph as a provider of an orientation */
  myb_graph_orientation pair = {graph, 17, TURNING};
  myb_orientation o;
  myb_rotation r;
  MT_CHECK(myb_orientation_graph(&pair, t, &o) == MYB_OK);
  MT_CHECK(myb_frame_graph_rotation(graph, 17, TURNING, t, &r) == MYB_OK);
  MT_CHECK(memcmp(&o.rotation, &r, sizeof r) == 0);
  /* there and back is the unit matrix, derivative included */
  myb_rotation back;
  MT_CHECK(myb_frame_graph_rotation(graph, TURNING, 17, t, &back) == MYB_OK);
  myb_rotation_compose(&back, &r, &r);
  double off = 0.0;
  for (int i = 0; i < 9; ++i) {
    off = fmax(off, fabs(r.m[i] - (i % 4 == 0 ? 1.0 : 0.0)));
    off = fmax(off, fabs(r.dm[i]) / 7.2921e-5);
  }
  MT_CHECK(off < 1e-15);
  myb_frame_graph_destroy(graph);
  myb_spk_destroy(spk);
  myb_pck_destroy(pck);
}

/* A provider that tells the instant it was asked for. */
typedef struct witness {
  myb_time asked;
  int calls;
} witness;

static myb_status witness_orientation(void* context, myb_time t, myb_orientation* orientation) {
  witness* w = context;
  w->asked = t;
  ++w->calls;
  memset(orientation, 0, sizeof *orientation);
  orientation->rotation.m[0] = orientation->rotation.m[4] = orientation->rotation.m[8] = 1.0;
  orientation->ut = t;
  return MYB_OK;
}

static myb_status failing_orientation(void* context, myb_time t, myb_orientation* orientation) {
  (void)context;
  (void)t;
  (void)orientation;
  return MYB_E_RANGE;
}

static void test_scales_and_errors(void) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/time-scales.tk", &diag) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/tdb-tt-fb1990.tk", &diag) == MYB_OK);
  myb_timesys* sys = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_scale tt;
  myb_scale tdb;
  MT_CHECK(myb_timesys_scale(sys, "TT", &tt) == MYB_OK);
  MT_CHECK(myb_timesys_scale(sys, "TDB", &tdb) == MYB_OK);
  /* a graph on TDB with one provider on TT and one on the scale of the graph */
  const myb_frame_graph_options options = {.timesys = sys, .scale = "TDB"};
  myb_frame_graph* graph = NULL;
  MT_CHECK(myb_frame_graph_create(&options, NULL, &graph) == MYB_OK);
  witness on_tt = {0};
  witness on_tdb = {0};
  MT_CHECK(myb_frame_graph_add(graph, 10, 1, witness_orientation, &on_tt, "TT") == MYB_OK);
  MT_CHECK(myb_frame_graph_add(graph, 11, 10, witness_orientation, &on_tdb, NULL) == MYB_OK);
  const myb_time t = myb_time_from_seconds(7.0e8 + 0.25);
  myb_rotation r;
  MT_CHECK(myb_frame_graph_rotation(graph, 11, 1, t, &r) == MYB_OK);
  myb_time expected;
  MT_CHECK(myb_timesys_convert(sys, t, tdb, tt, &expected) == MYB_OK);
  const double shift = myb_dd_to_double(myb_time_diff(on_tt.asked, t));
  printf("a provider on TT in a graph on TDB is asked %.9f s earlier (TDB - TT there)\n", -shift);
  MT_CHECK(on_tt.calls == 1 && on_tdb.calls == 1);
  MT_CHECK(myb_time_compare(on_tt.asked, expected) == 0 && fabs(shift) > 1e-4 &&
           fabs(shift) < 2e-3);
  MT_CHECK(myb_time_compare(on_tdb.asked, t) == 0);
  /* the same frame on both sides: nothing is asked */
  MT_CHECK(myb_frame_graph_rotation(graph, 11, 11, t, &r) == MYB_OK && on_tt.calls == 1);
  /* refused: a code of 0, a frame on itself, a second definition, a scale that is not known */
  MT_CHECK(myb_frame_graph_add(graph, 0, 1, witness_orientation, &on_tt, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_frame_graph_add(graph, 12, 12, witness_orientation, &on_tt, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_frame_graph_add(graph, 10, 1, witness_orientation, &on_tt, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_frame_graph_add(graph, 12, 1, witness_orientation, &on_tt, "XYZ") ==
           MYB_E_NOT_FOUND);
  /* two trees: no way from one to the other */
  myb_rotation unit = {.m = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}};
  MT_CHECK(myb_frame_graph_add_fixed(graph, 21, 20, &unit) == MYB_OK);
  MT_CHECK(myb_frame_graph_rotation(graph, 21, 11, t, &r) == MYB_E_NOT_FOUND);
  /* a cycle is found by its length */
  MT_CHECK(myb_frame_graph_add_fixed(graph, 31, 32, &unit) == MYB_OK);
  MT_CHECK(myb_frame_graph_add_fixed(graph, 32, 31, &unit) == MYB_OK);
  MT_CHECK(myb_frame_graph_rotation(graph, 31, 1, t, &r) == MYB_E_FORMAT);
  /* the status of a provider comes back */
  MT_CHECK(myb_frame_graph_add(graph, 41, 1, failing_orientation, NULL, NULL) == MYB_OK);
  MT_CHECK(myb_frame_graph_rotation(graph, 41, 10, t, &r) == MYB_E_RANGE);
  myb_frame_graph_destroy(graph);
  /* a time system needs the name of the scale; without one no provider may name a scale */
  const myb_frame_graph_options nameless = {.timesys = sys};
  MT_CHECK(myb_frame_graph_create(&nameless, NULL, &graph) == MYB_E_DOMAIN);
  MT_CHECK(myb_frame_graph_create(NULL, NULL, &graph) == MYB_OK);
  MT_CHECK(myb_frame_graph_add(graph, 10, 1, witness_orientation, &on_tt, "TT") == MYB_E_DOMAIN);
  /* a chain of the greatest length allowed */
  for (int32_t k = 100; k < 100 + MYB_FRAME_GRAPH_MAX_DEPTH - 1; ++k) {
    MT_CHECK(myb_frame_graph_add_fixed(graph, k + 1, k, &unit) == MYB_OK);
  }
  MT_CHECK(myb_frame_graph_rotation(graph, 100 + MYB_FRAME_GRAPH_MAX_DEPTH - 1, 100, t, &r) ==
           MYB_OK);
  MT_CHECK(myb_frame_graph_add_fixed(graph, 100 + MYB_FRAME_GRAPH_MAX_DEPTH, 100 +
                                     MYB_FRAME_GRAPH_MAX_DEPTH - 1, &unit) == MYB_OK);
  MT_CHECK(myb_frame_graph_rotation(graph, 100 + MYB_FRAME_GRAPH_MAX_DEPTH, 100, t, &r) ==
           MYB_E_FORMAT);
  myb_frame_graph_destroy(graph);
  myb_frame_graph_destroy(NULL);
  myb_timesys_destroy(sys);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_reference();
  test_scales_and_errors();
  return mt_finish("test_frame_graph");
}
