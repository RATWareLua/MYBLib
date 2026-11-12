/* SPDX-License-Identifier: BSD-3-Clause */
/* Spacecraft clocks and CK files against CSPICE N0067 (tests/data/ck, made by the oracle
 * script cspice_ck.py):
 *  - two clocks: ticks to time and back, strings to ticks and back;
 *  - a CK with segments of types 1, 2 and 3: what is found within a tolerance, at which
 *    instant, the matrix and the angular velocity (ckgpav_c, ckgp_c);
 *  - an instrument frame in a graph of frames (sxform_c);
 *  - a CK written here and read by CSPICE (ck_writer_cspice.ref, made from the file this
 *    test writes with MYB_CK_WRITER_DUMP=<path>);
 *  - a CK of type 5 written by CSPICE (t5.bc, t5_cspice.ref: cspice_ck456.py), and a CK of
 *    types 4, 5 and 6 written here and read by CSPICE (ck456_cspice.ref, from the file this
 *    test writes with MYB_CK456_DUMP=<path>);
 *  - what is refused. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/ck.h"
#include "myb/cmath.h"
#include "myb/frames.h"
#include "myb/io.h"
#include "myb/uncert.h"

#define CK_DIR MYB_TEST_DATA_DIR "/ck/"

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

/* What the comparison of pointings found. */
typedef struct pointing_stats {
  int found;
  int missing;
  double matrix;
  double rate;
  double ticks;
} pointing_stats;

/* One line "P instrument ticks tolerance need_rate found [ticks_out matrix rate]". */
static void check_pointing(const myb_ck* ck, const char* line, pointing_stats* stats) {
  const char* p = line + 1;
  double v[5];
  int ok = 1;
  for (int i = 0; i < 5; ++i) {
    ok = ok && next_number(&p, &v[i]);
  }
  MT_CHECK(ok);
  myb_ck_pointing pointing;
  const myb_status status =
      myb_ck_pointing_at(ck, (int32_t)v[0], v[1], v[2], (int)v[3], &pointing);
  if (v[4] == 0.0) {
    if (status != MYB_E_RANGE) {
      fprintf(stderr, "expected nothing: %s", line);
    }
    MT_CHECK(status == MYB_E_RANGE);
    ++stats->missing;
    return;
  }
  double e[13];
  for (int i = 0; i < 13; ++i) {
    ok = ok && next_number(&p, &e[i]);
  }
  MT_CHECK(ok);
  if (status != MYB_OK) {
    fprintf(stderr, "expected a pointing: %s", line);
  }
  MT_CHECK(status == MYB_OK);
  if (status != MYB_OK) {
    return;
  }
  ++stats->found;
  stats->ticks = fmax(stats->ticks, fabs(pointing.ticks - e[0]));
  for (int i = 0; i < 9; ++i) {
    stats->matrix = fmax(stats->matrix, fabs(pointing.rotation.m[i] - e[1 + i]));
  }
  if (v[3] != 0.0) {
    for (int i = 0; i < 3; ++i) {
      stats->rate = fmax(stats->rate, fabs(pointing.rate[i] - e[10 + i]));
    }
    MT_CHECK(pointing.has_rate);
  }
  /* the quaternion and the matrix say the same; the derivative is that of the rate */
  const double* q = pointing.q;
  MT_CHECK(fabs(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3] - 1.0) < 1e-15);
  MT_CHECK(fabs(pointing.rotation.m[0] - (1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]))) < 1e-15);
}

static void test_reference(void) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, CK_DIR "test.tsc", &diag) == MYB_OK);
  myb_sclk* tdb_clock = NULL;
  myb_sclk* tt_clock = NULL;
  MT_CHECK(myb_sclk_create(pool, 9099, NULL, &tdb_clock, &diag) == MYB_OK);
  MT_CHECK(myb_sclk_create(pool, 9098, NULL, &tt_clock, &diag) == MYB_OK);
  /* a clock the pool does not describe: refused, the output left as it was */
  myb_sclk* none = tt_clock;
  MT_CHECK(myb_sclk_create(pool, 9097, NULL, &none, &diag) == MYB_E_NOT_FOUND);
  MT_CHECK(none == tt_clock);
  myb_kpool_destroy(pool);
  MT_CHECK(myb_sclk_time_system(tdb_clock) == 1 && myb_sclk_time_system(tt_clock) == 2);
  MT_CHECK(myb_sclk_partitions(tdb_clock) == 2 && myb_sclk_partitions(tt_clock) == 1);
  myb_ck* ck = NULL;
  MT_CHECK(myb_ck_create(NULL, &ck) == MYB_OK);
  MT_CHECK(myb_ck_load_file(ck, CK_DIR "test.bc", &diag) == MYB_OK);
  MT_CHECK(myb_ck_segment_count(ck) == 7);
  myb_ck_segment_info info;
  myb_ck_segment(ck, 4, &info);
  MT_CHECK(info.instrument == -9099005 && info.type == 2 && info.frame == 13 && info.has_rates &&
           strcmp(info.name, "type 2") == 0 && info.start == 6.0e6);
  /* the instrument frame in a graph: the clock is on TDB, as the graph */
  myb_frame_graph* graph = NULL;
  MT_CHECK(myb_frame_graph_create(NULL, NULL, &graph) == MYB_OK);
  myb_ck_orientation instrument = {.ck = ck, .sclk = tdb_clock, .instrument = -9099001, .base = 1};
  MT_CHECK(myb_frame_graph_add(graph, -9099001, 1, myb_orientation_ck, &instrument, NULL) ==
           MYB_OK);

  FILE* f = fopen(CK_DIR "ck_cspice.ref", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  char line[2048];
  int clocks = 0;
  int frames = 0;
  double worst_time[2] = {0.0, 0.0};
  double worst_ticks[2] = {0.0, 0.0};
  double worst_frame = 0.0;
  double worst_frame_rate = 0.0;
  pointing_stats stats = {0};
  /* bound: one iteration per line of the file */
  while (fgets(line, sizeof line, f) != NULL) {
    const char kind = line[0];
    const char* p = line + 1;
    double v[20];
    if (kind == 'T' || kind == 'E' || kind == 'D' || kind == 'N') {
      MT_CHECK(next_number(&p, &v[0]) && next_number(&p, &v[1]));
      const myb_sclk* clock = v[0] == 9099.0 ? tdb_clock : tt_clock;
      const int which = v[0] == 9099.0 ? 0 : 1;
      ++clocks;
      if (kind == 'T') {
        MT_CHECK(next_number(&p, &v[2]));
        myb_time t;
        MT_CHECK(myb_sclk_to_time(clock, v[1], &t) == MYB_OK);
        worst_time[which] =
            fmax(worst_time[which],
                 fabs(myb_dd_to_double(myb_time_diff(t, myb_time_from_seconds(v[2])))));
      } else if (kind == 'E') {
        MT_CHECK(next_number(&p, &v[2]));
        double ticks = 0.0;
        MT_CHECK(myb_sclk_from_time(clock, myb_time_from_seconds(v[1]), &ticks) == MYB_OK);
        worst_ticks[which] = fmax(worst_ticks[which], fabs(ticks - v[2]));
      } else {
        while (*p == ' ') {
          ++p;
        }
        size_t length = strlen(p);
        while (length > 0 && (p[length - 1] == '\n' || p[length - 1] == '\r')) {
          --length;
        }
        if (kind == 'D') {
          char text[64];
          size_t written = 0;
          MT_CHECK(myb_sclk_decode(clock, v[1], text, sizeof text, &written) == MYB_OK);
          if (written != length || memcmp(text, p, length) != 0) {
            fprintf(stderr, "decode %.17g: '%s' against %s", v[1], text, p);
          }
          MT_CHECK(written == length && memcmp(text, p, length) == 0);
        } else {
          double ticks = -1.0;
          const myb_status status = myb_sclk_encode(clock, p, length, &ticks);
          if (status != MYB_OK || ticks != v[1]) {
            fprintf(stderr, "encode %s: %.17g (%s) against %.17g\n", p, ticks,
                    myb_status_name(status), v[1]);
          }
          MT_CHECK(status == MYB_OK && ticks == v[1]);
        }
      }
    } else if (kind == 'P') {
      check_pointing(ck, line, &stats);
    } else if (kind == 'F') {
      int ok = 1;
      for (int i = 0; i < 19; ++i) {
        ok = ok && next_number(&p, &v[i]);
      }
      MT_CHECK(ok);
      myb_rotation r;
      MT_CHECK(myb_frame_graph_rotation(graph, 1, -9099001, myb_time_from_seconds(v[0]), &r) ==
               MYB_OK);
      for (int i = 0; i < 9; ++i) {
        worst_frame = fmax(worst_frame, fabs(r.m[i] - v[1 + i]));
        worst_frame_rate = fmax(worst_frame_rate, fabs(r.dm[i] - v[10 + i]));
      }
      ++frames;
    }
  }
  fclose(f);
  printf("spacecraft clocks against CSPICE: %d conversions; time %.1e s (TDB clock), %.1e s "
         "(TT clock, through the TDB of SPICE); ticks back %.1e and %.1e\n",
         clocks, worst_time[0], worst_time[1], worst_ticks[0], worst_ticks[1]);
  printf("CK types 1, 2, 3 against CSPICE: %d pointings found and %d refused as there; matrix "
         "%.1e, angular velocity %.1e rad/s, instant %.1e ticks; the instrument frame in a "
         "graph: %d rotations, matrix %.1e, derivative %.1e 1/s\n",
         stats.found, stats.missing, stats.matrix, stats.rate, stats.ticks, frames, worst_frame,
         worst_frame_rate);
  MT_CHECK(clocks > 400 && stats.found == 966 && stats.missing == 336 && frames == 40);
  /* a double of seconds past J2000 resolves 1.2e-7 s now; the TT clock goes through the
   * TDB - TT of SPICE and back */
  MT_CHECK(worst_time[0] < 2e-7 && worst_time[1] < 1e-6);
  MT_CHECK(worst_ticks[0] < 1e-3 && worst_ticks[1] < 0.05);
  MT_CHECK(stats.matrix < 2e-15 && stats.rate < 1e-18 && stats.ticks == 0.0);
  MT_CHECK(worst_frame < 1e-12 && worst_frame_rate < 1e-15);

  /* refused */
  myb_ck_pointing pointing;
  MT_CHECK(myb_ck_pointing_at(ck, -9099001, 6.1e6, -1.0, 0, &pointing) == MYB_E_DOMAIN);
  MT_CHECK(myb_ck_pointing_at(ck, -1, 6.1e6, 0.0, 0, &pointing) == MYB_E_RANGE);
  /* rates asked of a segment without them */
  MT_CHECK(myb_ck_pointing_at(ck, -9099002, 6.2e6, 0.0, 1, &pointing) == MYB_E_RANGE);
  MT_CHECK(myb_ck_pointing_at(ck, -9099002, 6.2e6, 0.0, 0, &pointing) == MYB_OK &&
           !pointing.has_rate && pointing.frame == 17);
  /* a base frame other than the one promised to the graph */
  myb_ck_orientation wrong = {.ck = ck, .sclk = tdb_clock, .instrument = -9099002, .base = 1};
  myb_orientation o;
  myb_time inside;
  MT_CHECK(myb_sclk_to_time(tdb_clock, 6.2e6, &inside) == MYB_OK);
  MT_CHECK(myb_orientation_ck(&wrong, inside, &o) == MYB_E_FORMAT);
  wrong.base = 0;
  MT_CHECK(myb_orientation_ck(&wrong, inside, &o) == MYB_OK);
  /* outside the clock, outside the file */
  MT_CHECK(myb_orientation_ck(&wrong, myb_time_from_seconds(1.0e8), &o) == MYB_E_RANGE);
  MT_CHECK(myb_orientation_ck(&wrong, myb_time_from_seconds(6.0e8 + 5.0e4), &o) == MYB_E_RANGE);
  MT_CHECK(myb_ck_load_file(ck, CK_DIR "test.tsc", &diag) == MYB_E_FORMAT);
  MT_CHECK(myb_ck_load_file(ck, CK_DIR "none.bc", &diag) == MYB_E_IO);
  MT_CHECK(myb_ck_segment_count(ck) == 7);
  /* clock strings */
  double ticks = 0.0;
  myb_time t;
  char text[64];
  MT_CHECK(myb_sclk_encode(tdb_clock, "3/100:5:3", 9, &ticks) == MYB_E_DOMAIN);
  MT_CHECK(myb_sclk_encode(tdb_clock, "1/2000:0:0", 10, &ticks) == MYB_E_DOMAIN); /* past its end */
  MT_CHECK(myb_sclk_encode(tdb_clock, "100:5:3:1", 9, &ticks) == MYB_E_FORMAT);
  MT_CHECK(myb_sclk_encode(tdb_clock, "100:x", 5, &ticks) == MYB_E_FORMAT);
  MT_CHECK(myb_sclk_encode(tdb_clock, "100:", 4, &ticks) == MYB_E_FORMAT);
  MT_CHECK(myb_sclk_encode(tdb_clock, "", 0, &ticks) == MYB_E_FORMAT);
  MT_CHECK(myb_sclk_decode(tdb_clock, 100.0, text, 8, NULL) == MYB_E_OVERFLOW);
  MT_CHECK(myb_sclk_decode(tdb_clock, -1.0, text, sizeof text, NULL) == MYB_E_RANGE);
  MT_CHECK(myb_sclk_to_time(tdb_clock, -1.0, &t) == MYB_E_RANGE);
  MT_CHECK(myb_sclk_to_time(tdb_clock, 3.0e9, &t) == MYB_E_RANGE);
  MT_CHECK(myb_sclk_from_time(tdb_clock, myb_time_from_seconds(5.9e8), &ticks) == MYB_E_RANGE);
  myb_frame_graph_destroy(graph);
  myb_ck_destroy(ck);
  myb_sclk_destroy(tdb_clock);
  myb_sclk_destroy(tt_clock);
  myb_sclk_destroy(NULL);
  myb_ck_destroy(NULL);
}

/* A smooth orientation for the writer: three turns whose angles go with tau. */
static void sample(double tau, double phase, double q[4], double rate[3]) {
  const double a = 0.5 * (phase + 0.9 * tau);
  const double b = 0.5 * (0.4 + 0.5 * myb_sin(2.1 * tau + phase));
  const double qa[4] = {myb_cos(a), 0.0, 0.0, myb_sin(a)};
  const double qb[4] = {myb_cos(b), myb_sin(b), 0.0, 0.0};
  q[0] = 3.0 * (qa[0] * qb[0] - qa[3] * qb[3] * 0.0);
  q[1] = 3.0 * (qa[0] * qb[1]);
  q[2] = 3.0 * (qa[3] * qb[1]);
  q[3] = 3.0 * (qa[3] * qb[0]);
  rate[0] = 2.0e-3 * myb_cos(1.7 * tau + phase);
  rate[1] = 1.0e-3 * myb_sin(0.9 * tau);
  rate[2] = 3.0e-3 + 1.0e-3 * tau;
}

static void test_writer(void) {
  enum { N = 250, INTERVALS = 8 };
  static double ticks[N];
  static double q[4 * N];
  static double rates[3 * N];
  for (int k = 0; k < N; ++k) {
    ticks[k] = 5.0e6 + 9000.0 * k + 1000.0 * myb_sin(0.7 * k);
    sample((double)k / (double)N, 0.3, q + 4 * k, rates + 3 * k);
  }
  const double starts[3] = {ticks[0], ticks[101], ticks[220]};
  myb_ck_writer* writer = NULL;
  MT_CHECK(myb_ck_writer_create("myb test: CK writer", NULL, &writer) == MYB_OK);
  MT_CHECK(myb_ck_writer_size(writer) == 0);
  const myb_ck_segment_header h3 = {-77001, 1, "type 3, three intervals"};
  const myb_ck_segment_header h3b = {-77002, 17, "type 3 without rates"};
  const myb_ck_segment_header h1 = {-77003, 1, "type 1"};
  const myb_ck_segment_header h1b = {-77004, 2, "type 1 without rates"};
  const myb_ck_segment_header h2 = {-77005, 13, "type 2"};
  MT_CHECK(myb_ck_writer_add_pointings(writer, &h3, 3, ticks, q, rates, N, starts, 3) == MYB_OK);
  MT_CHECK(myb_ck_writer_add_pointings(writer, &h3b, 3, ticks, q, NULL, 40, NULL, 0) == MYB_OK);
  MT_CHECK(myb_ck_writer_add_pointings(writer, &h1, 1, ticks, q, rates, N, NULL, 0) == MYB_OK);
  MT_CHECK(myb_ck_writer_add_pointings(writer, &h1b, 1, ticks, q, NULL, 7, NULL, 0) == MYB_OK);
  double begin[INTERVALS];
  double end[INTERVALS];
  double per_tick[INTERVALS];
  for (int k = 0; k < INTERVALS; ++k) {
    begin[k] = 5.0e6 + 90000.0 * k;
    end[k] = begin[k] + 60000.0;
    per_tick[k] = 1.0 / 2048.0 * (1.0 + 1.0e-5 * k);
  }
  MT_CHECK(myb_ck_writer_add_intervals(writer, &h2, begin, end, q, rates, per_tick, INTERVALS) ==
           MYB_OK);
  /* refused: nothing is added */
  const size_t before = myb_ck_writer_size(writer);
  const myb_ck_segment_header nameless = {-77009, 1, ""};
  MT_CHECK(myb_ck_writer_add_pointings(writer, &nameless, 3, ticks, q, rates, N, NULL, 0) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_ck_writer_add_pointings(writer, &h3, 2, ticks, q, rates, N, NULL, 0) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_ck_writer_add_pointings(writer, &h3, 3, ticks, q, rates, 0, NULL, 0) ==
           MYB_E_DOMAIN);
  const double off_the_ticks[2] = {ticks[0], ticks[5] + 1.0};
  MT_CHECK(myb_ck_writer_add_pointings(writer, &h3, 3, ticks, q, rates, N, off_the_ticks, 2) ==
           MYB_E_DOMAIN);
  const double late_start[1] = {ticks[1]};
  MT_CHECK(myb_ck_writer_add_pointings(writer, &h3, 3, ticks, q, rates, N, late_start, 1) ==
           MYB_E_DOMAIN);
  const double backwards[3] = {5.0, 4.0, 6.0};
  MT_CHECK(myb_ck_writer_add_pointings(writer, &h3, 1, backwards, q, rates, 3, NULL, 0) ==
           MYB_E_DOMAIN);
  const double zero[4] = {0.0, 0.0, 0.0, 0.0};
  MT_CHECK(myb_ck_writer_add_pointings(writer, &h3, 1, ticks, zero, NULL, 1, NULL, 0) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_ck_writer_add_intervals(writer, &h2, end, begin, q, rates, per_tick, INTERVALS) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_ck_writer_size(writer) == before);

  const size_t bytes = myb_ck_writer_size(writer);
  MT_CHECK(bytes > 0 && bytes % 1024 == 0);
  unsigned char* image = malloc(bytes);
  MT_CHECK(image != NULL);
  MT_CHECK(myb_ck_writer_write(writer, image, bytes - 1) == MYB_E_OVERFLOW);
  MT_CHECK(myb_ck_writer_write(writer, image, bytes) == MYB_OK);
  MT_CHECK(memcmp(image, "DAF/CK  ", 8) == 0);
  const char* dump = getenv("MYB_CK_WRITER_DUMP");
  if (dump != NULL) {
    MT_CHECK(myb_ck_writer_write_file(writer, dump) == MYB_OK);
    printf("dump: %s\n", dump);
  }
  myb_ck* ck = NULL;
  myb_diag diag;
  MT_CHECK(myb_ck_create(NULL, &ck) == MYB_OK);
  MT_CHECK(myb_ck_load_memory(ck, image, bytes, &diag) == MYB_OK);
  MT_CHECK(myb_ck_segment_count(ck) == 5);
  myb_ck_segment_info info;
  myb_ck_segment(ck, 0, &info);
  MT_CHECK(info.instrument == -77001 && info.type == 3 && info.has_rates && info.frame == 1 &&
           info.start == ticks[0] && info.stop == ticks[N - 1]);
  myb_ck_segment(ck, 3, &info);
  MT_CHECK(info.instrument == -77004 && info.type == 1 && !info.has_rates && info.frame == 2);
  /* at its own instants a pointing comes back as it went in, made unit */
  double worst_own = 0.0;
  for (int k = 0; k < N; ++k) {
    myb_ck_pointing pointing;
    MT_CHECK(myb_ck_pointing_at(ck, -77001, ticks[k], 0.0, 1, &pointing) == MYB_OK);
    for (int i = 0; i < 4; ++i) {
      worst_own = fmax(worst_own, fabs(pointing.q[i] - q[4 * k + i] / 3.0));
    }
    MT_CHECK(memcmp(pointing.rate, rates + 3 * k, 3 * sizeof(double)) == 0);
  }
  MT_CHECK(worst_own < 3e-16);
  /* CSPICE reads the file */
  FILE* f = fopen(CK_DIR "ck_writer_cspice.ref", "r");
  MT_CHECK(f != NULL);
  pointing_stats stats = {0};
  if (f != NULL) {
    char line[2048];
    /* bound: one iteration per line of the file */
    while (fgets(line, sizeof line, f) != NULL) {
      if (line[0] == 'P') {
        check_pointing(ck, line, &stats);
      }
    }
    fclose(f);
  }
  printf("a CK written here and read by CSPICE: %d pointings found and %d refused as there; "
         "matrix %.1e, angular velocity %.1e rad/s\n",
         stats.found, stats.missing, stats.matrix, stats.rate);
  MT_CHECK(stats.found > 300 && stats.missing > 20);
  MT_CHECK(stats.matrix < 2e-15 && stats.rate < 1e-18 && stats.ticks == 0.0);
  myb_ck_destroy(ck);
  free(image);
  myb_ck_writer_destroy(writer);
  myb_ck_writer_destroy(NULL);
  MT_CHECK(myb_ck_writer_create("", NULL, &writer) == MYB_E_DOMAIN);
}

/* Every line "P" of a reference file against a CK. */
static pointing_stats compare_file(const myb_ck* ck, const char* path) {
  pointing_stats stats = {0};
  FILE* f = fopen(path, "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return stats;
  }
  char line[2048];
  /* bound: one iteration per line of the file */
  while (fgets(line, sizeof line, f) != NULL) {
    if (line[0] == 'P') {
      check_pointing(ck, line, &stats);
    }
  }
  fclose(f);
  return stats;
}

static void test_type5_by_cspice(void) {
  myb_ck* ck = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_ck_create(NULL, &ck) == MYB_OK);
  MT_CHECK(myb_ck_load_file(ck, CK_DIR "t5.bc", &diag) == MYB_OK);
  MT_CHECK(myb_ck_segment_count(ck) == 6);
  const pointing_stats stats = compare_file(ck, CK_DIR "t5_cspice.ref");
  printf("CK type 5 written by CSPICE, four subtypes: %d pointings found and %d refused as "
         "there; matrix %.1e, angular velocity %.1e rad/s, instant %.1e ticks\n",
         stats.found, stats.missing, stats.matrix, stats.rate, stats.ticks);
  MT_CHECK(stats.found > 900 && stats.missing > 200);
  MT_CHECK(stats.matrix < 2e-15 && stats.rate < 5e-16 && stats.ticks == 0.0);
  myb_ck_destroy(ck);
}

/* A smooth turning for the packets of types 5 and 6: the quaternion (of length 3), its
 * derivative per second, the angular velocity and its derivative. */
static void turning(double seconds, double phase, double out[14]) {
  const double h = 0.5;
  double q[3][4];
  double r[3][3];
  for (int k = 0; k < 3; ++k) {
    sample((seconds + (k - 1) * h) / 100000.0, phase, q[k], r[k]);
  }
  for (int i = 0; i < 4; ++i) {
    out[i] = q[1][i];
    out[4 + i] = (q[2][i] - q[0][i]) / (2.0 * h);
  }
  for (int i = 0; i < 3; ++i) {
    out[8 + i] = r[1][i];
    out[11 + i] = (r[2][i] - r[0][i]) / (2.0 * h);
  }
}

/* Packets of a subtype at the ticks, `seconds_per_tick` apart in time. */
static void make_packets(int subtype, const double* ticks, int count, double seconds_per_tick,
                         double phase, double* packets) {
  static const int size[4] = {8, 4, 14, 7};
  for (int k = 0; k < count; ++k) {
    double all[14];
    turning((ticks[k] - ticks[0]) * seconds_per_tick, phase, all);
    double* p = packets + size[subtype] * k;
    if (subtype == 0) {
      memcpy(p, all, 8 * sizeof(double));
    } else if (subtype == 1) {
      memcpy(p, all, 4 * sizeof(double));
    } else if (subtype == 2) {
      memcpy(p, all, 14 * sizeof(double));
    } else {
      memcpy(p, all, 4 * sizeof(double));
      memcpy(p + 4, all + 8, 3 * sizeof(double));
    }
  }
}

/* A line "R" of ck456_cspice.ref: a request in a gap answered from the next mini-segment of
 * a type 6 segment whose rate is not whole. The pointing is the same; CSPICE divides the
 * angular velocity by the rate rounded, the library by the rate. */
static void check_rounded(const myb_ck* ck, const char* line, int* lines, double* worst) {
  const char* p = line + 1;
  double v[20];
  int count = 0;
  /* bound: the numbers of the line */
  while (count < 20 && next_number(&p, &v[count])) {
    ++count;
  }
  MT_CHECK(count == 20 && v[4] == 1.0);
  if (count != 20) {
    return;
  }
  myb_ck_pointing pointing;
  MT_CHECK(myb_ck_pointing_at(ck, (int32_t)v[0], v[1], v[2], (int)v[3], &pointing) == MYB_OK);
  double e = fabs(pointing.ticks - v[5]);
  for (int i = 0; i < 9; ++i) {
    e = fmax(e, fabs(pointing.rotation.m[i] - v[6 + i]));
  }
  if (v[3] != 0.0) {
    const double rate = v[18];
    const double rounded = v[19];
    double differs = 0.0;
    for (int i = 0; i < 3; ++i) {
      e = fmax(e, fabs(pointing.rate[i] - v[15 + i] * rounded / rate));
      differs = fmax(differs, fabs(pointing.rate[i] - v[15 + i]));
    }
    MT_CHECK(differs > 1e-6);
  }
  *worst = fmax(*worst, e);
  ++*lines;
}

static void test_types456_writer(void) {
  myb_rng rng;
  myb_rng_seed(&rng, 20261112);
  myb_ck_writer* writer = NULL;
  MT_CHECK(myb_ck_writer_create("myb test: CK types 4, 5 and 6", NULL, &writer) == MYB_OK);
  /* type 4: series of varying lengths on intervals with gaps, with rates and without; and
   * 150 intervals for the directory of references */
  enum { MAX_COEFFICIENTS = 150 * 7 * 19 };
  static double coefficients[MAX_COEFFICIENTS];
  static double starts[150];
  static double stops[150];
  static int counts[150 * 7];
  const struct {
    int32_t instrument;
    int intervals;
    int has_rates;
    double first;
  } plan4[3] = {{-9099041, 6, 1, 2.0e7}, {-9099042, 3, 0, 2.1e7}, {-9099046, 150, 1, 2.2e7}};
  for (int s = 0; s < 3; ++s) {
    size_t used = 0;
    for (int k = 0; k < plan4[s].intervals; ++k) {
      starts[k] = plan4[s].first + 50000.0 * k;
      stops[k] = starts[k] + (plan4[s].intervals > 10 ? 50000.0 : 30000.0 + 3000.0 * k);
      for (int i = 0; i < 7; ++i) {
        const int n = i < 4 ? 3 + (k + 2 * i) % 10 : (plan4[s].has_rates ? 2 + (k + i) % 6 : 0);
        counts[7 * k + i] = n;
        for (int j = 0; j < n; ++j) {
          static const double centre[4] = {0.8, 0.3, -0.4, 0.2};
          const double scale = i < 4 ? 0.1 / (j + 1) : 1.0e-3 / (j + 1);
          double c = (2.0 * myb_rng_uniform(&rng) - 1.0) * scale;
          if (j == 0 && i < 4) {
            c += centre[i];
          }
          coefficients[used++] = c;
        }
      }
    }
    char name[40];
    snprintf(name, sizeof name, "type 4, %d intervals", plan4[s].intervals);
    const myb_ck_segment_header h = {plan4[s].instrument, 1, name};
    MT_CHECK(myb_ck_writer_add_chebyshev(writer, &h, starts, stops, counts, coefficients,
                                         (size_t)plan4[s].intervals,
                                         plan4[s].has_rates) == MYB_OK);
  }
  /* type 5, subtype 0, with a single-pointing interval */
  enum { N5 = 130 };
  static double ticks5[N5];
  static double packets5[8 * N5];
  for (int k = 0; k < N5; ++k) {
    ticks5[k] = 2.5e7 + 1000.0 * k + 200.0 * (k % 4);
  }
  make_packets(0, ticks5, N5, 0.1, 0.2, packets5);
  const double starts5[4] = {ticks5[0], ticks5[50], ticks5[51], ticks5[100]};
  const myb_ck_packets set5 = {0, 6, 0.1, ticks5, packets5, N5};
  const myb_ck_segment_header h5 = {-9099043, 1, "type 5 written here"};
  MT_CHECK(myb_ck_writer_add_packets(writer, &h5, &set5, 1, ticks5[0], ticks5[N5 - 1], starts5,
                                     4) == MYB_OK);
  /* type 6: four mini-segments of the four subtypes, with padding and a gap, the later one
   * taken on a boundary; the third has a rate that is not whole */
  static double ticks6[4][40];
  static double packets6[4][14 * 40];
  const struct {
    int subtype;
    int window;
    double rate;
    double first;
    double step;
    int count;
  } plan6[4] = {{0, 4, 1.0, 0.0, 1000.0, 40},
                {3, 8, 2.0, 34000.0, 1500.0, 18},
                {1, 6, 2.5, 70000.0, 800.0, 40},
                {2, 4, 0.5, 99000.0, 1200.0, 30}};
  myb_ck_packets minis[4];
  for (int m = 0; m < 4; ++m) {
    for (int k = 0; k < plan6[m].count; ++k) {
      ticks6[m][k] = 3.0e7 + plan6[m].first + plan6[m].step * k;
    }
    make_packets(plan6[m].subtype, ticks6[m], plan6[m].count, plan6[m].rate, 0.1 * m,
                 packets6[m]);
    minis[m] = (myb_ck_packets){plan6[m].subtype, plan6[m].window, plan6[m].rate, ticks6[m],
                                packets6[m], (size_t)plan6[m].count};
  }
  const double bounds6[5] = {3.0e7, 3.0e7 + 35000.0, 3.0e7 + 70000.0, 3.0e7 + 100000.0,
                             3.0e7 + 130000.0};
  const myb_ck_segment_header h6 = {-9099044, 1, "type 6, four subtypes"};
  MT_CHECK(myb_ck_writer_add_minisegments(writer, &h6, minis, bounds6, 4, 1, 1, 3.0e7 + 1000.0,
                                          3.0e7 + 129000.0) == MYB_OK);
  /* type 6 with 120 short mini-segments, the earlier one taken on a boundary */
  enum { MINIS = 120 };
  static double ticks_many[MINIS][3];
  static double packets_many[MINIS][7 * 3];
  static myb_ck_packets many[MINIS];
  static double bounds_many[MINIS + 1];
  for (int m = 0; m < MINIS; ++m) {
    bounds_many[m] = 4.0e7 + 3000.0 * m;
    for (int k = 0; k < 3; ++k) {
      ticks_many[m][k] = bounds_many[m] + 1000.0 * k + (k == 2 ? 500.0 : 0.0);
    }
    make_packets(3, ticks_many[m], 3, 0.25, 0.01 * m, packets_many[m]);
    many[m] = (myb_ck_packets){3, 2, 0.25, ticks_many[m], packets_many[m], 3};
  }
  bounds_many[MINIS] = 4.0e7 + 3000.0 * (double)MINIS;
  const myb_ck_segment_header h6b = {-9099045, 1, "type 6, 120 mini-segments"};
  MT_CHECK(myb_ck_writer_add_minisegments(writer, &h6b, many, bounds_many, MINIS, 1, 0,
                                          bounds_many[0], bounds_many[MINIS]) == MYB_OK);

  /* refused, nothing added */
  const size_t before = myb_ck_writer_size(writer);
  const myb_ck_segment_header h = {-9099049, 1, "refused"};
  int bad_counts[7] = {3, 3, 3, 3, 0, 0, 0};
  const double one_start = 1.0;
  const double one_stop = 2.0;
  MT_CHECK(myb_ck_writer_add_chebyshev(writer, &h, &one_start, &one_stop, bad_counts,
                                       coefficients, 1, 1) == MYB_E_DOMAIN);
  bad_counts[0] = 20;
  MT_CHECK(myb_ck_writer_add_chebyshev(writer, &h, &one_start, &one_stop, bad_counts,
                                       coefficients, 1, 0) == MYB_E_DOMAIN);
  bad_counts[0] = 3;
  MT_CHECK(myb_ck_writer_add_chebyshev(writer, &h, &one_stop, &one_start, bad_counts,
                                       coefficients, 1, 0) == MYB_E_DOMAIN);
  myb_ck_packets odd = set5;
  odd.window = 5;
  MT_CHECK(myb_ck_writer_add_packets(writer, &h, &odd, 1, ticks5[0], ticks5[N5 - 1], NULL, 0) ==
           MYB_E_DOMAIN);
  odd = set5;
  odd.window = 14;
  MT_CHECK(myb_ck_writer_add_packets(writer, &h, &odd, 1, ticks5[0], ticks5[N5 - 1], NULL, 0) ==
           MYB_E_DOMAIN);
  odd = set5;
  odd.subtype = 4;
  MT_CHECK(myb_ck_writer_add_packets(writer, &h, &odd, 1, ticks5[0], ticks5[N5 - 1], NULL, 0) ==
           MYB_E_DOMAIN);
  odd = set5;
  odd.seconds_per_tick = 0.0;
  MT_CHECK(myb_ck_writer_add_packets(writer, &h, &odd, 1, ticks5[0], ticks5[N5 - 1], NULL, 0) ==
           MYB_E_DOMAIN);
  const double off[2] = {ticks5[0], ticks5[3] + 1.0};
  MT_CHECK(myb_ck_writer_add_packets(writer, &h, &set5, 1, ticks5[0], ticks5[N5 - 1], off, 2) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_ck_writer_add_packets(writer, &h, &set5, 1, ticks5[5], ticks5[2], NULL, 0) ==
           MYB_E_DOMAIN);
  double shifted[5];
  memcpy(shifted, bounds6, sizeof shifted);
  shifted[1] = 3.0e7 + 33000.0; /* before the first tick of the second mini-segment */
  MT_CHECK(myb_ck_writer_add_minisegments(writer, &h, minis, shifted, 4, 1, 1, 3.0e7 + 1000.0,
                                          3.0e7 + 129000.0) == MYB_E_DOMAIN);
  MT_CHECK(myb_ck_writer_add_minisegments(writer, &h, minis, bounds6, 4, 1, 1, 3.0e7 - 1.0,
                                          3.0e7 + 129000.0) == MYB_E_DOMAIN);
  myb_ck_packets single[1] = {minis[0]};
  single[0].count = 1;
  MT_CHECK(myb_ck_writer_add_minisegments(writer, &h, single, bounds6, 1, 1, 1, 3.0e7,
                                          3.0e7 + 1.0) == MYB_E_DOMAIN);
  MT_CHECK(myb_ck_writer_size(writer) == before);

  const size_t bytes = myb_ck_writer_size(writer);
  unsigned char* image = malloc(bytes);
  MT_CHECK(image != NULL && myb_ck_writer_write(writer, image, bytes) == MYB_OK);
  const char* dump = getenv("MYB_CK456_DUMP");
  if (dump != NULL) {
    MT_CHECK(myb_ck_writer_write_file(writer, dump) == MYB_OK);
    printf("dump: %s\n", dump);
  }
  myb_ck* ck = NULL;
  myb_diag diag;
  MT_CHECK(myb_ck_create(NULL, &ck) == MYB_OK);
  MT_CHECK(myb_ck_load_memory(ck, image, bytes, &diag) == MYB_OK);
  MT_CHECK(myb_ck_segment_count(ck) == 6);
  /* type 4 at the middle of an interval (s = 0, where T_j = 1, 0, -1, 0, 1, ...): the series
   * of the last segment written, made unit */
  myb_ck_pointing pointing;
  MT_CHECK(myb_ck_pointing_at(ck, -9099046, 2.2e7 + 25000.0, 0.0, 1, &pointing) == MYB_OK);
  {
    double q[4];
    double rate[3];
    size_t at = 0;
    for (int i = 0; i < 7; ++i) {
      double sum = 0.0;
      for (int j = 0; j < counts[i]; ++j) {
        static const double t_at_zero[4] = {1.0, 0.0, -1.0, 0.0};
        sum += coefficients[at + (size_t)j] * t_at_zero[j % 4];
      }
      at += (size_t)counts[i];
      if (i < 4) {
        q[i] = sum;
      } else {
        rate[i - 4] = sum;
      }
    }
    const double norm = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    double e = 0.0;
    for (int i = 0; i < 4; ++i) {
      e = fmax(e, fabs(pointing.q[i] - q[i] / norm));
    }
    for (int i = 0; i < 3; ++i) {
      e = fmax(e, fabs(pointing.rate[i] - rate[i]));
    }
    MT_CHECK(e < 1e-15 && pointing.ticks == 2.2e7 + 25000.0);
  }
  const pointing_stats stats = compare_file(ck, CK_DIR "ck456_cspice.ref");
  int rounded = 0;
  double worst_rounded = 0.0;
  FILE* f = fopen(CK_DIR "ck456_cspice.ref", "r");
  if (f != NULL) {
    char line[2048];
    /* bound: one iteration per line of the file */
    while (fgets(line, sizeof line, f) != NULL) {
      if (line[0] == 'R') {
        check_rounded(ck, line, &rounded, &worst_rounded);
      }
    }
    fclose(f);
  }
  printf("CK types 4, 5 and 6 written here and read by CSPICE: %d pointings found and %d "
         "refused as there; matrix %.1e, angular velocity %.1e rad/s, instant %.1e ticks; %d "
         "requests in a gap answered from a mini-segment of rate 2.5, where CSPICE divides by "
         "3: %.1e\n",
         stats.found, stats.missing, stats.matrix, stats.rate, stats.ticks, rounded,
         worst_rounded);
  MT_CHECK(stats.found > 500 && stats.missing > 50 && rounded > 0);
  MT_CHECK(stats.matrix < 2e-15 && stats.rate < 1e-16 && stats.ticks == 0.0 &&
           worst_rounded < 2e-16);
  myb_ck_destroy(ck);
  free(image);
  myb_ck_writer_destroy(writer);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_reference();
  test_writer();
  test_type5_by_cspice();
  test_types456_writer();
  return mt_finish("test_ck");
}
