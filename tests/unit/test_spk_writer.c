/* SPDX-License-Identifier: BSD-3-Clause */
/* Writing SPK files (myb_spk_writer): a file of 34 segments (two summary records) — discrete
 * states of types 9 and 13, Chebyshev polynomials of types 2 and 3 fitted to a two-body orbit —
 * is written, read back by myb's own reader and compared
 *  - with what was written: the states at the epochs of types 9 and 13 to rounding, the orbit
 *    between them and under the polynomials to the accuracy the writer reports;
 *  - with CSPICE N0067 reading the very same file (spkgeo_c through spiceypy on the bytes this
 *    test writes; tests/data/spk/writer_cspice.ref).
 * MYB_SPK_WRITER_DUMP=<path> writes the file and the list of check epochs for the oracle. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/io.h"
#include "myb/kepler.h"
#include "myb/kin.h"

enum { STATES = 250, SMALL = 30, CHECKS = 40 };

static const double mu = 3.986004418e14;

typedef struct orbit {
  myb_pv pv0;
} orbit;

static myb_status orbit_state(void* context, myb_time t, double y[6]) {
  const orbit* o = context;
  myb_pv pv;
  const myb_status status = myb_kepler_propagate(mu, &o->pv0, myb_time_seconds(t), &pv);
  memcpy(y, pv.r, 3 * sizeof(double));
  memcpy(y + 3, pv.v, 3 * sizeof(double));
  return status;
}

static void make_orbit(orbit* o, double a, double e, double nu) {
  const myb_elements el = {.p = a * (1.0 - e * e), .e = e, .i = 0.9, .raan = 0.4, .argp = 1.1,
                           .nu = nu};
  MT_CHECK(myb_pv_from_elements(mu, &el, &o->pv0) == MYB_OK);
}

static double epoch_of(int k) {
  return 100.0 * k + 7.0 * (double)((k * 37) % 11) / 11.0;
}

static double distance(const double a[6], const double b[6], double* dv) {
  *dv = hypot(hypot(a[3] - b[3], a[4] - b[4]), a[5] - b[5]);
  return hypot(hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2]);
}

/* The same relative to the lengths of the position and of the velocity of b. */
static double relative(const double a[6], const double b[6], double* dv) {
  const double dx = distance(a, b, dv);
  *dv /= hypot(hypot(b[3], b[4]), b[5]);
  return dx / hypot(hypot(b[0], b[1]), b[2]);
}

int main(void) {
  orbit o9;
  orbit o13;
  orbit o2;
  orbit o3;
  make_orbit(&o9, 7.2e6, 0.05, 0.3);
  make_orbit(&o13, 2.66e7, 0.7, -2.0);
  make_orbit(&o2, 4.2164e7, 0.001, 1.0);
  make_orbit(&o3, 1.2e7, 0.3, 2.5);

  myb_spk_writer* writer = NULL;
  MT_CHECK(myb_spk_writer_create("myb test_spk_writer", NULL, &writer) == MYB_OK);
  MT_CHECK(myb_spk_writer_size(writer) == 0);
  static double epochs[STATES];
  static double states9[6 * STATES];
  static double states13[6 * STATES];
  for (int k = 0; k < STATES; ++k) {
    epochs[k] = epoch_of(k);
    MT_CHECK(orbit_state(&o9, myb_time_from_seconds(epochs[k]), states9 + 6 * k) == MYB_OK);
    MT_CHECK(orbit_state(&o13, myb_time_from_seconds(epochs[k]), states13 + 6 * k) == MYB_OK);
  }
  const myb_spk_segment_header h9 = {-100, 399, 1, "LAGRANGE 7"};
  const myb_spk_segment_header h13 = {-101, 399, 1, "HERMITE 7"};
  const myb_spk_segment_header h2 = {-102, 399, 1, "CHEBYSHEV POSITION 12"};
  const myb_spk_segment_header h3 = {-103, 399, 1, "CHEBYSHEV STATE 10"};
  MT_CHECK(myb_spk_writer_add_states(writer, &h9, 9, 7, epochs, states9, STATES) == MYB_OK);
  MT_CHECK(myb_spk_writer_add_states(writer, &h13, 13, 7, epochs, states13, STATES) == MYB_OK);
  double error2 = -1.0;
  double error3 = -1.0;
  MT_CHECK(myb_spk_writer_add_chebyshev(writer, &h2, 2, orbit_state, &o2, 0.0, 20000.0, 20, 12,
                                        &error2) == MYB_OK);
  MT_CHECK(myb_spk_writer_add_chebyshev(writer, &h3, 3, orbit_state, &o3, -3000.0, 13000.0, 16,
                                        10, &error3) == MYB_OK);
  /* more segments than one summary record holds */
  for (int k = 0; k < SMALL; ++k) {
    char name[41];
    snprintf(name, sizeof name, "SMALL %d", k);
    const myb_spk_segment_header h = {-200 - k, 399, 1, name};
    MT_CHECK(myb_spk_writer_add_states(writer, &h, 9, 1, epochs + k, states9 + 6 * k, 2) ==
             MYB_OK);
  }
  const size_t size = myb_spk_writer_size(writer);
  MT_CHECK(size > 0 && size % 1024 == 0);
  unsigned char* file = malloc(size);
  MT_CHECK(file != NULL);
  if (file == NULL) {
    return 1;
  }
  MT_CHECK(myb_spk_writer_write(writer, file, size - 1) == MYB_E_OVERFLOW);
  MT_CHECK(myb_spk_writer_write(writer, file, size) == MYB_OK);
  MT_CHECK(memcmp(file, "DAF/SPK ", 8) == 0 && memcmp(file + 88, "LTL-IEEE", 8) == 0);
  /* the same bytes through the file system */
  const char* path = "test_spk_writer.bsp";
  MT_CHECK(myb_spk_writer_write_file(writer, path) == MYB_OK);
  FILE* stored = fopen(path, "rb");
  MT_CHECK(stored != NULL);
  if (stored != NULL) {
    unsigned char* copy = malloc(size + 1);
    MT_CHECK(copy != NULL && fread(copy, 1, size + 1, stored) == size);
    MT_CHECK(copy != NULL && memcmp(copy, file, size) == 0);
    free(copy);
    fclose(stored);
  }

  /* read back */
  myb_spk* spk = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  MT_CHECK(myb_spk_load_memory(spk, file, size, &diag) == MYB_OK);
  MT_CHECK(myb_spk_segment_count(spk) == 4 + SMALL);
  myb_spk_segment_info info;
  myb_spk_segment(spk, 0, &info);
  MT_CHECK(info.target == -100 && info.center == 399 && info.frame == 1 && info.type == 9 &&
           info.start == epochs[0] && info.stop == epochs[STATES - 1] &&
           strcmp(info.name, "LAGRANGE 7") == 0);
  myb_spk_segment(spk, 3, &info);
  MT_CHECK(info.target == -103 && info.type == 3 && info.start == -3000.0 &&
           info.stop == 13000.0 && strcmp(info.name, "CHEBYSHEV STATE 10") == 0);
  myb_spk_segment(spk, 4 + SMALL - 1, &info);
  MT_CHECK(info.target == -200 - (SMALL - 1) && info.type == 9 &&
           strcmp(info.name, "SMALL 29") == 0);

  /* types 9 and 13: at the epochs the states come back (one division and one multiplication
   * by 1000: an ulp); between them the orbit */
  double worst_ulp = 0.0;
  double between9 = 0.0;
  double between13 = 0.0;
  for (int k = 0; k < STATES; ++k) {
    double y[6];
    double dv = 0.0;
    MT_CHECK(myb_spk_state(spk, -100, 399, myb_time_from_seconds(epochs[k]), y, NULL) == MYB_OK);
    worst_ulp = fmax(worst_ulp, relative(y, states9 + 6 * k, &dv) / 2.220446049250313e-16);
    worst_ulp = fmax(worst_ulp, dv / 2.220446049250313e-16);
    MT_CHECK(myb_spk_state(spk, -101, 399, myb_time_from_seconds(epochs[k]), y, NULL) == MYB_OK);
    worst_ulp = fmax(worst_ulp, relative(y, states13 + 6 * k, &dv) / 2.220446049250313e-16);
    worst_ulp = fmax(worst_ulp, dv / 2.220446049250313e-16);
    if (k + 1 < STATES) {
      const double t = 0.5 * (epochs[k] + epochs[k + 1]);
      double truth[6];
      MT_CHECK(myb_spk_state(spk, -100, 399, myb_time_from_seconds(t), y, NULL) == MYB_OK);
      MT_CHECK(orbit_state(&o9, myb_time_from_seconds(t), truth) == MYB_OK);
      between9 = fmax(between9, distance(y, truth, &dv));
      MT_CHECK(myb_spk_state(spk, -101, 399, myb_time_from_seconds(t), y, NULL) == MYB_OK);
      MT_CHECK(orbit_state(&o13, myb_time_from_seconds(t), truth) == MYB_OK);
      between13 = fmax(between13, distance(y, truth, &dv));
    }
  }
  printf("types 9 and 13: states at the epochs within %.2f ulp of their length of those "
         "written; between epochs 100 s apart %.1e m (Lagrange 7) and %.1e m (Hermite 7) from "
         "the orbit\n",
         worst_ulp, between9, between13);
  MT_CHECK(worst_ulp <= 4.0);
  MT_CHECK(between9 < 0.1 && between13 < 1e-2);

  /* types 2 and 3: the polynomials against the provider, within what the writer reported */
  double worst2 = 0.0;
  double worst3 = 0.0;
  double worst3v = 0.0;
  double worst2v = 0.0;
  for (int k = 0; k <= 400; ++k) {
    double y[6];
    double truth[6];
    double dv = 0.0;
    const double t2 = 20000.0 * k / 400.0;
    MT_CHECK(myb_spk_state(spk, -102, 399, myb_time_from_seconds(t2), y, NULL) == MYB_OK);
    MT_CHECK(orbit_state(&o2, myb_time_from_seconds(t2), truth) == MYB_OK);
    worst2 = fmax(worst2, distance(y, truth, &dv));
    worst2v = fmax(worst2v, dv);
    const double t3 = -3000.0 + 16000.0 * k / 400.0;
    MT_CHECK(myb_spk_state(spk, -103, 399, myb_time_from_seconds(t3), y, NULL) == MYB_OK);
    MT_CHECK(orbit_state(&o3, myb_time_from_seconds(t3), truth) == MYB_OK);
    worst3 = fmax(worst3, distance(y, truth, &dv));
    worst3v = fmax(worst3v, dv);
  }
  printf("type 2 (degree 12, 1000 s): %.1e m, %.1e m/s from the orbit (reported %.1e m); type 3 "
         "(degree 10, 1000 s): %.1e m, %.1e m/s (reported %.1e m)\n",
         worst2, worst2v, error2, worst3, worst3v, error3);
  MT_CHECK(error2 >= 0.0 && error3 >= 0.0);
  MT_CHECK(worst2 <= 2.0 * error2 + 1e-6 && worst3 <= 2.0 * error3 + 1e-6);
  MT_CHECK(worst2 < 1e-3 && worst3 < 1.0);
  /* a higher degree fits better */
  myb_spk_writer* finer = NULL;
  double coarse_error = 0.0;
  double fine_error = 0.0;
  MT_CHECK(myb_spk_writer_create(NULL, NULL, &finer) == MYB_OK);
  MT_CHECK(myb_spk_writer_add_chebyshev(finer, &h3, 3, orbit_state, &o3, -3000.0, 13000.0, 16, 6,
                                        &coarse_error) == MYB_OK);
  MT_CHECK(myb_spk_writer_add_chebyshev(finer, &h3, 3, orbit_state, &o3, -3000.0, 13000.0, 16, 14,
                                        &fine_error) == MYB_OK);
  MT_CHECK(coarse_error > 100.0 * error3 && fine_error < error3);
  myb_spk_writer_destroy(finer);

  /* check epochs: for the oracle, and against its answers */
  const char* dump = getenv("MYB_SPK_WRITER_DUMP");
  FILE* list = dump != NULL ? fopen(dump, "w") : NULL;
  FILE* reference = fopen(MYB_TEST_DATA_DIR "/spk/writer_cspice.ref", "rb");
  const int targets[5] = {-100, -101, -102, -103, -207};
  const double first[5] = {0.0, 0.0, 0.0, -3000.0, epochs[7]};
  const double last[5] = {24900.0, 24900.0, 20000.0, 13000.0, epochs[8]};
  double worst_position = 0.0;
  double worst_velocity = 0.0;
  int compared = 0;
  for (int s = 0; s < 5; ++s) {
    for (int k = 0; k < CHECKS; ++k) {
      /* inside the coverage; the first and the last at its ends */
      const double t = first[s] + (last[s] - first[s]) * ((double)k / (CHECKS - 1));
      double y[6];
      MT_CHECK(myb_spk_state(spk, targets[s], 399, myb_time_from_seconds(t), y, NULL) == MYB_OK);
      if (list != NULL) {
        char text[MYB_DOUBLE_CHARS];
        (void)myb_format_double(t, text);
        fprintf(list, "%d %s\n", targets[s], text);
      }
      if (reference != NULL) {
        char line[512];
        double v[8];
        int n = 0;
        if (fgets(line, sizeof line, reference) != NULL) {
          const char* p = line;
          /* bound: eight numbers */
          while (n < 8) {
            while (*p == ' ') {
              ++p;
            }
            size_t used = 0;
            if (myb_parse_double(p, strlen(p), &v[n], &used) != MYB_OK) {
              break;
            }
            p += used;
            ++n;
          }
        }
        MT_CHECK(n == 8 && (int)v[0] == targets[s] && v[1] == t);
        if (n == 8) {
          const double r = hypot(hypot(y[0], y[1]), y[2]);
          const double speed = hypot(hypot(y[3], y[4]), y[5]);
          const double theirs[6] = {1e3 * v[2], 1e3 * v[3], 1e3 * v[4],
                                    1e3 * v[5], 1e3 * v[6], 1e3 * v[7]};
          double dv = 0.0;
          worst_position = fmax(worst_position, distance(y, theirs, &dv) / r);
          worst_velocity = fmax(worst_velocity, dv / speed);
          ++compared;
        }
      }
    }
  }
  if (list != NULL) {
    fclose(list);
    printf("check epochs written to %s; the file is %s\n", dump, path);
  }
  if (reference != NULL) {
    fclose(reference);
    printf("CSPICE N0067 on the same file (%d states): positions %.1e, velocities %.1e\n",
           compared, worst_position, worst_velocity);
    MT_CHECK(compared == 5 * CHECKS);
    MT_CHECK(worst_position < 1e-14);
    MT_CHECK(worst_velocity < 1e-12); /* the two differentiate the polynomials differently */
  } else {
    MT_CHECK(dump != NULL); /* the reference belongs to the repository */
  }

  /* errors */
  const myb_spk_segment_header bad_name = {-1, 399, 1, ""};
  const myb_spk_segment_header long_name = {-1, 399, 1,
                                            "A NAME OF MORE THAN FORTY CHARACTERS IN ALL"};
  MT_CHECK(myb_spk_writer_add_states(writer, &bad_name, 9, 1, epochs, states9, 2) == MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_add_states(writer, &long_name, 9, 1, epochs, states9, 2) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_add_states(writer, &h9, 8, 1, epochs, states9, 2) == MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_add_states(writer, &h9, 13, 4, epochs, states9, 20) == MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_add_states(writer, &h9, 9, 28, epochs, states9, STATES) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_add_states(writer, &h9, 9, 7, epochs, states9, 5) == MYB_E_DOMAIN);
  const double backwards[2] = {10.0, 5.0};
  MT_CHECK(myb_spk_writer_add_states(writer, &h9, 9, 1, backwards, states9, 2) == MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_add_chebyshev(writer, &h2, 5, orbit_state, &o2, 0.0, 1.0, 1, 3, NULL) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_add_chebyshev(writer, &h2, 2, orbit_state, &o2, 1.0, 1.0, 1, 3, NULL) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_add_chebyshev(writer, &h2, 2, orbit_state, &o2, 0.0, 1.0, 0, 3, NULL) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_spk_segment_count(spk) == 4 + SMALL); /* nothing was added by the failures */
  myb_spk_writer* empty = NULL;
  MT_CHECK(myb_spk_writer_create("\tname", NULL, &empty) == MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_create(NULL, NULL, &empty) == MYB_OK);
  MT_CHECK(myb_spk_writer_write(empty, file, size) == MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_write_file(empty, path) == MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_write_file(writer, "no-such-directory/x.bsp") == MYB_E_IO);
  myb_spk_writer_destroy(empty);
  myb_spk_writer_destroy(NULL);

  myb_spk_destroy(spk);
  free(file);
  myb_spk_writer_destroy(writer);
  if (dump == NULL) {
    remove(path); /* the oracle reads it otherwise */
  }
  return mt_finish("test_spk_writer");
}
