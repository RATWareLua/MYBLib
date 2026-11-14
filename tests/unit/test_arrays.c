/* SPDX-License-Identifier: BSD-3-Clause */
/* The array forms (the functions named *_many) against the single calls they repeat: the
 * outputs of one call over a series are the bytes of the single call in a loop (memcmp), for
 * the time scales and the calendar, SPK states, the rotations and orientation providers of
 * the frames, Kepler's propagation and the conversions of elements. count 0 does nothing and
 * accepts NULL arrays; a failing element stops the call with its status and its index, the
 * outputs before it filled as the loop fills them. (NULL arrays with count > 0 break the
 * contract: tests/unit/array_death.c.) */
#include <math.h>
#include <string.h>

#include "mtest.h"
#include "myb/myb.h"

enum { N = 12, TURNING = 1400099, TURNING_CLASS = 31099 };

static const size_t untouched = (size_t)-7;

static myb_kpool* pool_of(const char* const* files, size_t count) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  for (size_t i = 0; i < count; ++i) {
    MT_CHECK(myb_kpool_load_file(pool, files[i], &diag) == MYB_OK);
  }
  return pool;
}

/* ---- Time --------------------------------------------------------------------------------- */

static void test_time(void) {
  const char* files[] = {MYB_DATA_DIR "/time/time-scales.tk", MYB_DATA_DIR "/time/tdb-tt-fb1990.tk",
                         MYB_DATA_DIR "/time/utc-tai.tk"};
  myb_kpool* pool = pool_of(files, 3);
  myb_timesys* sys = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_scale tt = {0};
  myb_scale tdb = {0};
  myb_scale tai = {0};
  myb_scale tcb = {0};
  MT_CHECK(myb_timesys_scale(sys, "TT", &tt) == MYB_OK);
  MT_CHECK(myb_timesys_scale(sys, "TDB", &tdb) == MYB_OK);
  MT_CHECK(myb_timesys_scale(sys, "TAI", &tai) == MYB_OK);
  MT_CHECK(myb_timesys_scale(sys, "TCB", &tcb) == MYB_OK);

  /* seconds -> instants and back */
  double seconds[N];
  for (int k = 0; k < N; ++k) {
    seconds[k] = 8.0e8 + 3.7e6 * k + 0.123456789 * k;
  }
  myb_time many[N];
  myb_time loop[N];
  memset(many, 0, sizeof many);
  memset(loop, 0, sizeof loop);
  myb_time_from_seconds_many(seconds, N, many);
  for (int k = 0; k < N; ++k) {
    loop[k] = myb_time_from_seconds(seconds[k]);
  }
  MT_CHECK(memcmp(many, loop, sizeof many) == 0);
  double back_many[N];
  double back_loop[N];
  myb_time_seconds_many(many, N, back_many);
  for (int k = 0; k < N; ++k) {
    back_loop[k] = myb_time_seconds(loop[k]);
  }
  MT_CHECK(memcmp(back_many, back_loop, sizeof back_many) == 0);

  /* calendars */
  myb_calendar labels_many[N];
  myb_calendar labels_loop[N];
  memset(labels_many, 0, sizeof labels_many);
  memset(labels_loop, 0, sizeof labels_loop);
  myb_time_to_calendar_many(many, N, labels_many);
  for (int k = 0; k < N; ++k) {
    myb_time_to_calendar(many[k], &labels_loop[k]);
  }
  MT_CHECK(memcmp(labels_many, labels_loop, sizeof labels_many) == 0);
  myb_time from_labels_many[N];
  myb_time from_labels_loop[N];
  size_t failed = untouched;
  MT_CHECK(myb_calendar_to_time_many(labels_many, N, from_labels_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_calendar_to_time(&labels_many[k], &from_labels_loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(from_labels_many, from_labels_loop, sizeof from_labels_many) == 0);
  myb_calendar wrong[N];
  memcpy(wrong, labels_many, sizeof wrong);
  wrong[5].month = 13;
  failed = untouched;
  MT_CHECK(myb_calendar_to_time_many(wrong, N, from_labels_many, &failed) == MYB_E_DOMAIN);
  MT_CHECK(failed == 5);
  MT_CHECK(memcmp(from_labels_many, from_labels_loop, 5 * sizeof(myb_time)) == 0);

  /* conversions between scales, also in place */
  const myb_scale pairs[4][2] = {{tt, tdb}, {tdb, tt}, {tt, tai}, {tdb, tcb}};
  for (int p = 0; p < 4; ++p) {
    myb_time out_many[N];
    myb_time out_loop[N];
    memset(out_many, 0, sizeof out_many);
    memset(out_loop, 0, sizeof out_loop);
    failed = untouched;
    MT_CHECK(myb_timesys_convert_many(sys, many, N, pairs[p][0], pairs[p][1], out_many, &failed) ==
             MYB_OK);
    MT_CHECK(failed == N);
    for (int k = 0; k < N; ++k) {
      MT_CHECK(myb_timesys_convert(sys, many[k], pairs[p][0], pairs[p][1], &out_loop[k]) == MYB_OK);
    }
    MT_CHECK(memcmp(out_many, out_loop, sizeof out_many) == 0);
    myb_time in_place[N];
    memcpy(in_place, many, sizeof in_place);
    MT_CHECK(myb_timesys_convert_many(sys, in_place, N, pairs[p][0], pairs[p][1], in_place, NULL) ==
             MYB_OK);
    MT_CHECK(memcmp(in_place, out_loop, sizeof in_place) == 0);
  }

  /* UTC labels both ways */
  myb_calendar utc_many[N];
  myb_calendar utc_loop[N];
  memset(utc_many, 0, sizeof utc_many);
  memset(utc_loop, 0, sizeof utc_loop);
  failed = untouched;
  MT_CHECK(myb_timesys_to_utc_many(sys, many, N, tt, utc_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_timesys_to_utc(sys, many[k], tt, &utc_loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(utc_many, utc_loop, sizeof utc_many) == 0);
  myb_time tdb_many[N];
  myb_time tdb_loop[N];
  memset(tdb_many, 0, sizeof tdb_many);
  memset(tdb_loop, 0, sizeof tdb_loop);
  failed = untouched;
  MT_CHECK(myb_timesys_from_utc_many(sys, utc_many, N, tdb, tdb_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_timesys_from_utc(sys, &utc_many[k], tdb, &tdb_loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(tdb_many, tdb_loop, sizeof tdb_many) == 0);
  /* 23:59:60 on a day without a leap second, then a label before the table */
  myb_calendar bad[N];
  memcpy(bad, utc_many, sizeof bad);
  bad[7].hour = 23;
  bad[7].minute = 59;
  bad[7].second = 60.5;
  failed = untouched;
  MT_CHECK(myb_timesys_from_utc_many(sys, bad, N, tdb, tdb_many, &failed) == MYB_E_DOMAIN);
  MT_CHECK(failed == 7);
  MT_CHECK(memcmp(tdb_many, tdb_loop, 7 * sizeof(myb_time)) == 0);
  bad[0].year = 1900;
  MT_CHECK(myb_timesys_from_utc_many(sys, bad, N, tdb, tdb_many, &failed) == MYB_E_RANGE);
  MT_CHECK(failed == 0);

  /* nothing to do: NULL arrays, the index is 0 */
  failed = untouched;
  MT_CHECK(myb_timesys_convert_many(sys, NULL, 0, tt, tdb, NULL, &failed) == MYB_OK);
  MT_CHECK(failed == 0);
  failed = untouched;
  MT_CHECK(myb_timesys_from_utc_many(sys, NULL, 0, tt, NULL, &failed) == MYB_OK);
  MT_CHECK(failed == 0);
  MT_CHECK(myb_timesys_to_utc_many(sys, NULL, 0, tt, NULL, NULL) == MYB_OK);
  MT_CHECK(myb_calendar_to_time_many(NULL, 0, NULL, NULL) == MYB_OK);
  myb_time_from_seconds_many(NULL, 0, NULL);
  myb_time_seconds_many(NULL, 0, NULL);
  myb_time_to_calendar_many(NULL, 0, NULL);
  myb_timesys_destroy(sys);
}

/* ---- SPK ---------------------------------------------------------------------------------- */

static void test_spk(void) {
  myb_spk* spk = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  MT_CHECK(myb_spk_load_file(spk, MYB_TEST_DATA_DIR "/spk/t02_03.bsp", &diag) == MYB_OK);
  myb_spk_segment_info info;
  myb_spk_segment(spk, 0, &info);
  myb_time t[N];
  for (int k = 0; k < N; ++k) {
    t[k] = myb_time_from_seconds(info.start + (info.stop - info.start) * (k + 0.5) / (double)N);
  }
  double many[6 * N];
  double loop[6 * N];
  int32_t frames_many[N];
  int32_t frames_loop[N];
  memset(many, 0, sizeof many);
  memset(loop, 0, sizeof loop);
  size_t failed = untouched;
  MT_CHECK(myb_spk_state_many(spk, info.target, info.center, t, N, many, frames_many, &failed) ==
           MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_spk_state(spk, info.target, info.center, t[k], loop + 6 * k, &frames_loop[k]) ==
             MYB_OK);
  }
  MT_CHECK(memcmp(many, loop, sizeof many) == 0);
  MT_CHECK(memcmp(frames_many, frames_loop, sizeof frames_many) == 0);
  MT_CHECK(myb_spk_state_many(spk, info.target, info.center, t, N, many, NULL, NULL) == MYB_OK);
  MT_CHECK(memcmp(many, loop, sizeof many) == 0);

  /* through a graph of frames (here the segments' own frame: the identity) */
  myb_frame_graph* graph = NULL;
  MT_CHECK(myb_frame_graph_create(NULL, NULL, &graph) == MYB_OK);
  double in_many[6 * N];
  double in_loop[6 * N];
  memset(in_many, 0, sizeof in_many);
  memset(in_loop, 0, sizeof in_loop);
  failed = untouched;
  MT_CHECK(myb_spk_state_in_many(spk, graph, info.target, info.center, t, N, info.frame, in_many,
                                 &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_spk_state_in(spk, graph, info.target, info.center, t[k], info.frame,
                              in_loop + 6 * k) == MYB_OK);
  }
  MT_CHECK(memcmp(in_many, in_loop, sizeof in_many) == 0);

  /* an epoch beyond the coverage */
  myb_time beyond[N];
  memcpy(beyond, t, sizeof beyond);
  beyond[9] = myb_time_from_seconds(info.stop + 1.0e9);
  failed = untouched;
  MT_CHECK(myb_spk_state_many(spk, info.target, info.center, beyond, N, many, NULL, &failed) ==
           MYB_E_RANGE);
  MT_CHECK(failed == 9);
  MT_CHECK(memcmp(many, loop, 9 * 6 * sizeof(double)) == 0);
  failed = untouched;
  MT_CHECK(myb_spk_state_in_many(spk, graph, info.target, info.center, beyond, N, info.frame,
                                 in_many, &failed) == MYB_E_RANGE);
  MT_CHECK(failed == 9);
  failed = untouched;
  MT_CHECK(myb_spk_state_many(spk, info.target, info.center, NULL, 0, NULL, NULL, &failed) ==
           MYB_OK);
  MT_CHECK(failed == 0);
  MT_CHECK(myb_spk_state_in_many(spk, graph, info.target, info.center, NULL, 0, 1, NULL, NULL) ==
           MYB_OK);
  myb_frame_graph_destroy(graph);
  myb_spk_destroy(spk);
}

/* ---- Frames ------------------------------------------------------------------------------- */

static void same_rotations(const myb_rotation* a, const myb_rotation* b, size_t count) {
  MT_CHECK(memcmp(a, b, count * sizeof(myb_rotation)) == 0);
}

static void test_earth(void) {
  const char* files[] = {MYB_DATA_DIR "/time/time-scales.tk", MYB_DATA_DIR "/time/tdb-tt-fb1990.tk",
                         MYB_DATA_DIR "/time/utc-tai.tk",     MYB_DATA_DIR "/frames/iau2006-cip.tk",
                         MYB_DATA_DIR "/frames/iau1980.tk",   MYB_DATA_DIR "/frames/eop-tides.tk"};
  myb_kpool* pool = pool_of(files, 6);
  myb_diag diag = {0};
  myb_timesys* sys = NULL;
  myb_cip* cip = NULL;
  myb_equinox* equinox = NULL;
  myb_eop* eop = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, &diag) == MYB_OK);
  MT_CHECK(myb_cip_create(pool, "IAU2006A", NULL, &cip, &diag) == MYB_OK);
  MT_CHECK(myb_equinox_create(pool, "IAU1980", NULL, &equinox, &diag) == MYB_OK);
  MT_CHECK(myb_eop_create(NULL, &eop) == MYB_OK);
  MT_CHECK(myb_eop_load_c04_file(eop, MYB_TEST_DATA_DIR "/eop/eopc04-excerpt.txt", &diag) ==
           MYB_OK);
  MT_CHECK(myb_eop_set_tides(eop, pool, "IERS10_EOPT", &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_scale tt = {0};
  MT_CHECK(myb_timesys_scale(sys, "TT", &tt) == MYB_OK);

  /* TT instants over 2020 (inside the excerpt of the C04 series) */
  myb_time t[N];
  for (int k = 0; k < N; ++k) {
    t[k] = myb_time_from_seconds(6.32e8 + 2.5e6 * k + 0.375 * k);
  }

  /* the pole and the origin */
  double xys_many[3 * N];
  double xys_loop[3 * N];
  double rates_many[3 * N];
  double rates_loop[3 * N];
  myb_cip_xys_many(cip, t, N, xys_many, rates_many);
  for (int k = 0; k < N; ++k) {
    myb_cip_xys(cip, t[k], xys_loop + 3 * k, rates_loop + 3 * k);
  }
  MT_CHECK(memcmp(xys_many, xys_loop, sizeof xys_many) == 0);
  MT_CHECK(memcmp(rates_many, rates_loop, sizeof rates_many) == 0);
  myb_cip_xys_many(cip, t, N, xys_many, NULL);
  MT_CHECK(memcmp(xys_many, xys_loop, sizeof xys_many) == 0);

  /* the EOP and UT1 */
  myb_eop_values eop_many[N];
  myb_eop_values eop_loop[N];
  myb_time ut1_many[N];
  myb_time ut1_loop[N];
  memset(ut1_many, 0, sizeof ut1_many);
  memset(ut1_loop, 0, sizeof ut1_loop);
  size_t failed = untouched;
  MT_CHECK(myb_eop_at_many(eop, sys, t, N, tt, eop_many, ut1_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_eop_at(eop, sys, t[k], tt, &eop_loop[k], &ut1_loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(eop_many, eop_loop, sizeof eop_many) == 0);
  MT_CHECK(memcmp(ut1_many, ut1_loop, sizeof ut1_many) == 0);
  MT_CHECK(myb_eop_at_many(eop, sys, t, N, tt, eop_many, NULL, NULL) == MYB_OK);
  MT_CHECK(memcmp(eop_many, eop_loop, sizeof eop_many) == 0);
  myb_time outside[N];
  memcpy(outside, t, sizeof outside);
  outside[4] = myb_time_from_seconds(3.0e8); /* 2009: before the excerpt */
  failed = untouched;
  MT_CHECK(myb_eop_at_many(eop, sys, outside, N, tt, eop_many, ut1_many, &failed) == MYB_E_RANGE);
  MT_CHECK(failed == 4);
  MT_CHECK(memcmp(eop_many, eop_loop, 4 * sizeof(myb_eop_values)) == 0);

  /* celestial -> terrestrial, with the EOP and without */
  myb_rotation many[N];
  myb_rotation loop[N];
  myb_cip_terrestrial_many(cip, t, ut1_loop, eop_loop, N, many);
  for (int k = 0; k < N; ++k) {
    myb_cip_terrestrial(cip, t[k], ut1_loop[k], &eop_loop[k], &loop[k]);
  }
  same_rotations(many, loop, N);
  myb_cip_terrestrial_many(cip, t, ut1_loop, NULL, N, many);
  for (int k = 0; k < N; ++k) {
    myb_cip_terrestrial(cip, t[k], ut1_loop[k], NULL, &loop[k]);
  }
  same_rotations(many, loop, N);

  /* states carried by the rotations, also in place */
  double states[6 * N];
  for (int k = 0; k < 6 * N; ++k) {
    states[k] = (k % 6 < 3 ? 6.9e6 : 7.4e3) * (1.0 + 0.01 * k) * (k % 2 == 0 ? 1.0 : -1.0);
  }
  double carried_many[6 * N];
  double carried_loop[6 * N];
  myb_rotation_apply_state_many(many, states, N, carried_many);
  for (int k = 0; k < N; ++k) {
    myb_rotation_apply_state(&many[k], states + 6 * k, carried_loop + 6 * k);
  }
  MT_CHECK(memcmp(carried_many, carried_loop, sizeof carried_many) == 0);
  double in_place[6 * N];
  memcpy(in_place, states, sizeof in_place);
  myb_rotation_apply_state_many(many, in_place, N, in_place);
  MT_CHECK(memcmp(in_place, carried_loop, sizeof in_place) == 0);
  myb_rotation_apply_position_many(many, states, N, carried_many);
  for (int k = 0; k < N; ++k) {
    myb_rotation_apply_position(&many[k], states + 3 * k, carried_loop + 3 * k);
  }
  MT_CHECK(memcmp(carried_many, carried_loop, 3 * N * sizeof(double)) == 0);

  /* the equinox-based systems */
  const myb_equinox_frame kinds[3] = {MYB_EQUINOX_MOD, MYB_EQUINOX_TOD, MYB_EQUINOX_TEME};
  for (int i = 0; i < 3; ++i) {
    myb_equinox_rotation_many(equinox, kinds[i], t, N, many);
    for (int k = 0; k < N; ++k) {
      myb_equinox_rotation(equinox, kinds[i], t[k], &loop[k]);
    }
    same_rotations(many, loop, N);
  }
  myb_equinox_teme_pef_many(equinox, ut1_loop, N, many);
  for (int k = 0; k < N; ++k) {
    myb_equinox_teme_pef(equinox, ut1_loop[k], &loop[k]);
  }
  same_rotations(many, loop, N);

  /* an orientation provider sampled: the model with observed EOP */
  myb_cip_orientation earth = {cip, eop, sys, tt};
  myb_orientation o_many[N];
  myb_orientation o_loop[N];
  memset(o_many, 0, sizeof o_many);
  memset(o_loop, 0, sizeof o_loop);
  failed = untouched;
  MT_CHECK(myb_orientation_many(myb_orientation_cip, &earth, t, N, o_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_orientation_cip(&earth, t[k], &o_loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(o_many, o_loop, sizeof o_many) == 0);
  failed = untouched;
  MT_CHECK(myb_orientation_many(myb_orientation_cip, &earth, outside, N, o_many, &failed) ==
           MYB_E_RANGE);
  MT_CHECK(failed == 4);
  MT_CHECK(memcmp(o_many, o_loop, 4 * sizeof(myb_orientation)) == 0);

  /* nothing to do */
  failed = untouched;
  MT_CHECK(myb_eop_at_many(eop, sys, NULL, 0, tt, NULL, NULL, &failed) == MYB_OK);
  MT_CHECK(failed == 0);
  failed = untouched;
  MT_CHECK(myb_orientation_many(myb_orientation_cip, &earth, NULL, 0, NULL, &failed) == MYB_OK);
  MT_CHECK(failed == 0);
  myb_cip_xys_many(cip, NULL, 0, NULL, NULL);
  myb_cip_terrestrial_many(cip, NULL, NULL, NULL, 0, NULL);
  myb_equinox_rotation_many(equinox, MYB_EQUINOX_TOD, NULL, 0, NULL);
  myb_equinox_teme_pef_many(equinox, NULL, 0, NULL);
  myb_rotation_apply_state_many(NULL, NULL, 0, NULL);
  myb_rotation_apply_position_many(NULL, NULL, 0, NULL);
  myb_eop_destroy(eop);
  myb_equinox_destroy(equinox);
  myb_cip_destroy(cip);
  myb_timesys_destroy(sys);
}

static void test_bodies(void) {
  /* a rotation model of the WGCCRE: Mars */
  const char* wgccre_files[] = {MYB_TEST_DATA_DIR "/wgccre/pck00011.tpc"};
  myb_kpool* pool = pool_of(wgccre_files, 1);
  myb_diag diag = {0};
  myb_wgccre* mars = NULL;
  MT_CHECK(myb_wgccre_create(pool, 499, NULL, &mars, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_time t[N];
  for (int k = 0; k < N; ++k) {
    t[k] = myb_time_from_seconds(-3.0e8 + 5.1e7 * k + 0.0625 * k);
  }
  myb_rotation many[N];
  myb_rotation loop[N];
  myb_wgccre_rotation_many(mars, t, N, many);
  for (int k = 0; k < N; ++k) {
    myb_wgccre_rotation(mars, t[k], &loop[k]);
  }
  same_rotations(many, loop, N);
  myb_wgccre_orientation mars_orientation = {mars, NULL};
  myb_orientation o_many[N];
  myb_orientation o_loop[N];
  memset(o_many, 0, sizeof o_many);
  memset(o_loop, 0, sizeof o_loop);
  MT_CHECK(myb_orientation_many(myb_orientation_wgccre, &mars_orientation, t, N, o_many, NULL) ==
           MYB_OK);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_orientation_wgccre(&mars_orientation, t[k], &o_loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(o_many, o_loop, sizeof o_many) == 0);
  myb_wgccre_rotation_many(mars, NULL, 0, NULL);
  myb_wgccre_destroy(mars);

  /* a frame of a binary PCK, alone and in a graph with the inertial frames of NAIF */
  myb_pck* pck = NULL;
  MT_CHECK(myb_pck_create(NULL, &pck) == MYB_OK);
  MT_CHECK(myb_pck_load_file(pck, MYB_TEST_DATA_DIR "/frames/graph.bpc", &diag) == MYB_OK);
  for (int k = 0; k < N; ++k) {
    t[k] = myb_time_from_seconds(600002748.0 + 160.0 * k + 0.25 * k);
  }
  int32_t bases_many[N];
  int32_t bases_loop[N];
  size_t failed = untouched;
  MT_CHECK(myb_pck_rotation_many(pck, TURNING_CLASS, t, N, many, bases_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_pck_rotation(pck, TURNING_CLASS, t[k], &loop[k], &bases_loop[k]) == MYB_OK);
  }
  same_rotations(many, loop, N);
  MT_CHECK(memcmp(bases_many, bases_loop, sizeof bases_many) == 0);
  MT_CHECK(myb_pck_rotation_many(pck, TURNING_CLASS, t, N, many, NULL, NULL) == MYB_OK);
  same_rotations(many, loop, N);
  myb_time beyond[N];
  memcpy(beyond, t, sizeof beyond);
  beyond[6] = myb_time_from_seconds(600002748.0 + 1.0e10);
  failed = untouched;
  MT_CHECK(myb_pck_rotation_many(pck, TURNING_CLASS, beyond, N, many, NULL, &failed) ==
           MYB_E_RANGE);
  MT_CHECK(failed == 6);
  same_rotations(many, loop, 6);

  const char* frame_files[] = {MYB_DATA_DIR "/frames/naif-inertial.tk"};
  pool = pool_of(frame_files, 1);
  myb_frame_graph* graph = NULL;
  MT_CHECK(myb_frame_graph_create(NULL, NULL, &graph) == MYB_OK);
  MT_CHECK(myb_frame_graph_add_tkframe(graph, pool, 17, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_pck_orientation turning = {pck, TURNING_CLASS, NULL};
  MT_CHECK(myb_frame_graph_add(graph, TURNING, 1, myb_orientation_pck, &turning, NULL) == MYB_OK);
  failed = untouched;
  MT_CHECK(myb_frame_graph_rotation_many(graph, 17, TURNING, t, N, many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_frame_graph_rotation(graph, 17, TURNING, t[k], &loop[k]) == MYB_OK);
  }
  same_rotations(many, loop, N);
  failed = untouched;
  MT_CHECK(myb_frame_graph_rotation_many(graph, 17, TURNING, beyond, N, many, &failed) ==
           MYB_E_RANGE);
  MT_CHECK(failed == 6);
  same_rotations(many, loop, 6);
  failed = untouched;
  MT_CHECK(myb_frame_graph_rotation_many(graph, 17, 99, t, N, many, &failed) == MYB_E_NOT_FOUND);
  MT_CHECK(failed == 0);
  failed = untouched;
  MT_CHECK(myb_frame_graph_rotation_many(graph, 17, TURNING, NULL, 0, NULL, &failed) == MYB_OK);
  MT_CHECK(failed == 0);
  MT_CHECK(myb_pck_rotation_many(pck, TURNING_CLASS, NULL, 0, NULL, NULL, NULL) == MYB_OK);
  myb_frame_graph_destroy(graph);
  myb_pck_destroy(pck);
}

/* ---- Kepler ------------------------------------------------------------------------------- */

static void test_kepler(void) {
  const double mu = 3.986004418e14;
  const myb_pv pv0 = {{6.9e6, 1.2e5, -3.4e5}, {-150.0, 7.3e3, 1.1e3}};
  double dt[N];
  for (int k = 0; k < N; ++k) {
    dt[k] = -2.0e4 + 6.1e3 * k + 0.001 * k;
  }
  dt[3] = 0.0;
  myb_pv many[N];
  myb_pv loop[N];
  size_t failed = untouched;
  MT_CHECK(myb_kepler_propagate_many(mu, &pv0, dt, N, many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_kepler_propagate(mu, &pv0, dt[k], &loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(many, loop, sizeof many) == 0);
  double stm_many[36 * N];
  double stm_loop[36 * N];
  myb_pv with_stm[N];
  failed = untouched;
  MT_CHECK(myb_kepler_stm_many(mu, &pv0, dt, N, with_stm, stm_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  myb_pv with_stm_loop[N];
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_kepler_stm(mu, &pv0, dt[k], &with_stm_loop[k], stm_loop + 36 * k) == MYB_OK);
  }
  MT_CHECK(memcmp(with_stm, with_stm_loop, sizeof with_stm) == 0);
  MT_CHECK(memcmp(stm_many, stm_loop, sizeof stm_many) == 0);
  double bad_dt[N];
  memcpy(bad_dt, dt, sizeof bad_dt);
  bad_dt[8] = (double)NAN;
  failed = untouched;
  MT_CHECK(myb_kepler_propagate_many(mu, &pv0, bad_dt, N, many, &failed) == MYB_E_DOMAIN);
  MT_CHECK(failed == 8);
  MT_CHECK(memcmp(many, loop, 8 * sizeof(myb_pv)) == 0);
  failed = untouched;
  MT_CHECK(myb_kepler_stm_many(mu, &pv0, bad_dt, N, with_stm, stm_many, &failed) == MYB_E_DOMAIN);
  MT_CHECK(failed == 8);
  failed = untouched;
  MT_CHECK(myb_kepler_propagate_many(-mu, &pv0, dt, N, many, &failed) == MYB_E_DOMAIN);
  MT_CHECK(failed == 0);

  /* elements of the states and states of the elements */
  myb_elements el_many[N];
  myb_elements el_loop[N];
  failed = untouched;
  MT_CHECK(myb_elements_from_pv_many(mu, loop, N, el_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_elements_from_pv(mu, &loop[k], &el_loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(el_many, el_loop, sizeof el_many) == 0);
  myb_pv back_many[N];
  myb_pv back_loop[N];
  failed = untouched;
  MT_CHECK(myb_pv_from_elements_many(mu, el_loop, N, back_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_pv_from_elements(mu, &el_loop[k], &back_loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(back_many, back_loop, sizeof back_many) == 0);
  myb_equinoctial eq_many[N];
  myb_equinoctial eq_loop[N];
  failed = untouched;
  MT_CHECK(myb_equinoctial_from_pv_many(mu, loop, N, eq_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_equinoctial_from_pv(mu, &loop[k], &eq_loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(eq_many, eq_loop, sizeof eq_many) == 0);
  failed = untouched;
  MT_CHECK(myb_pv_from_equinoctial_many(mu, eq_loop, N, back_many, &failed) == MYB_OK);
  MT_CHECK(failed == N);
  for (int k = 0; k < N; ++k) {
    MT_CHECK(myb_pv_from_equinoctial(mu, &eq_loop[k], &back_loop[k]) == MYB_OK);
  }
  MT_CHECK(memcmp(back_many, back_loop, sizeof back_many) == 0);

  /* failures: a negative semi-latus rectum, a state at the centre, a negative p again */
  myb_elements bad_el[N];
  memcpy(bad_el, el_loop, sizeof bad_el);
  bad_el[2].p = -1.0;
  failed = untouched;
  MT_CHECK(myb_pv_from_elements_many(mu, bad_el, N, back_many, &failed) == MYB_E_DOMAIN);
  MT_CHECK(failed == 2);
  myb_pv bad_pv[N];
  memcpy(bad_pv, loop, sizeof bad_pv);
  memset(&bad_pv[10], 0, sizeof bad_pv[10]);
  failed = untouched;
  MT_CHECK(myb_elements_from_pv_many(mu, bad_pv, N, el_many, &failed) == MYB_E_DOMAIN);
  MT_CHECK(failed == 10);
  MT_CHECK(memcmp(el_many, el_loop, 10 * sizeof(myb_elements)) == 0);
  failed = untouched;
  MT_CHECK(myb_equinoctial_from_pv_many(mu, bad_pv, N, eq_many, &failed) == MYB_E_DOMAIN);
  MT_CHECK(failed == 10);
  myb_equinoctial bad_eq[N];
  memcpy(bad_eq, eq_loop, sizeof bad_eq);
  bad_eq[1].p = -1.0;
  failed = untouched;
  MT_CHECK(myb_pv_from_equinoctial_many(mu, bad_eq, N, back_many, &failed) == MYB_E_DOMAIN);
  MT_CHECK(failed == 1);

  /* nothing to do */
  failed = untouched;
  MT_CHECK(myb_kepler_propagate_many(mu, &pv0, NULL, 0, NULL, &failed) == MYB_OK);
  MT_CHECK(failed == 0);
  MT_CHECK(myb_kepler_stm_many(mu, &pv0, NULL, 0, NULL, NULL, NULL) == MYB_OK);
  MT_CHECK(myb_elements_from_pv_many(mu, NULL, 0, NULL, NULL) == MYB_OK);
  MT_CHECK(myb_pv_from_elements_many(mu, NULL, 0, NULL, NULL) == MYB_OK);
  MT_CHECK(myb_equinoctial_from_pv_many(mu, NULL, 0, NULL, NULL) == MYB_OK);
  MT_CHECK(myb_pv_from_equinoctial_many(mu, NULL, 0, NULL, NULL) == MYB_OK);
}

int main(void) {
  test_time();
  test_spk();
  test_earth();
  test_bodies();
  test_kepler();
  return mt_finish("test_arrays");
}
