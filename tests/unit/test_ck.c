/* SPDX-License-Identifier: BSD-3-Clause */
/* Spacecraft clocks and CK files against CSPICE N0067 (tests/data/ck, made by the oracle
 * script cspice_ck.py):
 *  - two clocks: ticks to time and back, strings to ticks and back;
 *  - a CK with segments of types 1, 2 and 3: what is found within a tolerance, at which
 *    instant, the matrix and the angular velocity (ckgpav_c, ckgp_c);
 *  - an instrument frame in a graph of frames (sxform_c);
 *  - a CK written here and read by CSPICE (ck_writer_cspice.ref, made from the file this
 *    test writes with MYB_CK_WRITER_DUMP=<path>);
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

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_reference();
  test_writer();
  return mt_finish("test_ck");
}
