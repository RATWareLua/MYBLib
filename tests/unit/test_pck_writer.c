/* SPDX-License-Identifier: BSD-3-Clause */
/* Segments written word for word and binary PCK files (myb_spk_writer_add_raw, myb/kin.h;
 * myb_pck_writer, myb/frames.h):
 *  - the test kernels of SPK types 1, 2, 3, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20 copied segment
 *    by segment: the copy is read as the original, bit for bit;
 *  - a frame that turns uniformly about a tilted pole is given back by the reader to
 *    rounding, whatever the degree; a frame that librates, to the error the writer reports,
 *    which falls with the degree;
 *  - the file is read by CSPICE (pxform_c through a frame of class 2;
 *    tests/data/frames/pck_writer_cspice.ref, made by the oracle from the file this test
 *    writes with MYB_PCK_WRITER_DUMP=<path>);
 *  - what is refused. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/frames.h"
#include "myb/kin.h"

static const double pi = 3.141592653589793;

/* ---- Raw copies ------------------------------------------------------------------------------ */

static void test_raw_copy(void) {
  static const char* const files[] = {"t01.bsp", "t02_03.bsp", "t05.bsp", "t08_09_12_13.bsp",
                                      "t14.bsp", "t17.bsp",    "t18_19.bsp", "t20.bsp",
                                      "chain_a.bsp", "eros.bsp"};
  enum { FILES = sizeof files / sizeof files[0], MOST_WORDS = 1 << 20 };
  static double words[MOST_WORDS];
  static unsigned char image[16 << 20];
  unsigned long types = 0;
  int segments = 0;
  int states = 0;
  int different = 0;
  for (int file = 0; file < (int)FILES; ++file) {
    char path[512];
    snprintf(path, sizeof path, MYB_TEST_DATA_DIR "/spk/%s", files[file]);
    myb_spk* spk = NULL;
    myb_diag diag;
    MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
    MT_CHECK(myb_spk_load_file(spk, path, &diag) == MYB_OK);
    myb_spk_writer* writer = NULL;
    MT_CHECK(myb_spk_writer_create(files[file], NULL, &writer) == MYB_OK);
    const size_t count = myb_spk_segment_count(spk);
    for (size_t index = 0; index < count; ++index) {
      myb_spk_segment_info info;
      myb_spk_segment(spk, index, &info);
      const size_t size = myb_spk_segment_size(spk, index);
      MT_CHECK(size <= MOST_WORDS);
      MT_CHECK(myb_spk_segment_data(spk, index, 0, size, words) == MYB_OK);
      const myb_spk_segment_header header = {
          .target = info.target, .center = info.center, .frame = info.frame, .name = info.name};
      MT_CHECK(myb_spk_writer_add_raw(writer, &header, info.type, info.start, info.stop, words,
                                      size) == MYB_OK);
      types |= 1ul << info.type;
      ++segments;
    }
    const size_t bytes = myb_spk_writer_size(writer);
    MT_CHECK(bytes > 0 && bytes <= sizeof image);
    MT_CHECK(myb_spk_writer_write(writer, image, sizeof image) == MYB_OK);
    myb_spk* copy = NULL;
    MT_CHECK(myb_spk_create(NULL, &copy) == MYB_OK);
    MT_CHECK(myb_spk_load_memory(copy, image, bytes, &diag) == MYB_OK);
    MT_CHECK(myb_spk_segment_count(copy) == count);
    for (size_t index = 0; index < count; ++index) {
      myb_spk_segment_info info;
      myb_spk_segment(spk, index, &info);
      for (int k = 0; k < 7; ++k) {
        const double et = info.start + (info.stop - info.start) * ((double)k + 0.5) / 7.0;
        double a[6] = {0};
        double b[6] = {0};
        const myb_status sa =
            myb_spk_state(spk, info.target, info.center, myb_time_from_seconds(et), a, NULL);
        const myb_status sb =
            myb_spk_state(copy, info.target, info.center, myb_time_from_seconds(et), b, NULL);
        different += sa != sb || memcmp(a, b, sizeof a) != 0;
        states += sa == MYB_OK;
      }
    }
    const char* dump = getenv("MYB_PCK_WRITER_DUMP");
    if (dump != NULL) {
      char out[600];
      snprintf(out, sizeof out, "%s.copy.%s", dump, files[file]);
      MT_CHECK(myb_spk_writer_write_file(writer, out) == MYB_OK);
    }
    myb_spk_destroy(copy);
    myb_spk_writer_destroy(writer);
    myb_spk_destroy(spk);
  }
  int kinds = 0;
  for (int type = 1; type < 32; ++type) {
    kinds += (types >> type) & 1ul ? 1 : 0;
  }
  printf("raw copies: %d segments of %d types in %d files, %d states of the copies, %d differ "
         "from the files\n",
         segments, kinds, (int)FILES, states, different);
  MT_CHECK(kinds >= 12 && states > 100 && different == 0);
  MT_CHECK((types >> 1) & 1ul && (types >> 5) & 1ul && (types >> 20) & 1ul);

  /* what is refused */
  myb_spk_writer* writer = NULL;
  MT_CHECK(myb_spk_writer_create(NULL, NULL, &writer) == MYB_OK);
  const myb_spk_segment_header header = {.target = -1, .center = 399, .frame = 1, .name = "x"};
  const double word = 1.0;
  MT_CHECK(myb_spk_writer_add_raw(writer, &header, 0, 0.0, 1.0, &word, 1) == MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_add_raw(writer, &header, 1, 1.0, 0.0, &word, 1) == MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_add_raw(writer, &header, 1, 0.0, 1.0, &word, 0) == MYB_E_DOMAIN);
  MT_CHECK(myb_spk_writer_size(writer) == 0);
  myb_spk_writer_destroy(writer);
}

/* ---- Binary PCK ------------------------------------------------------------------------------ */

typedef struct turning {
  double angles[3];    /* at t = 0 */
  double rates[3];     /* rad/s */
  double amplitude[3]; /* of the libration */
  double frequency[3]; /* rad/s */
} turning;

static void turning_angles(const turning* f, double t, double angles[3], double rates[3]) {
  for (int i = 0; i < 3; ++i) {
    angles[i] = f->angles[i] + f->rates[i] * t + f->amplitude[i] * sin(f->frequency[i] * t);
    rates[i] = f->rates[i] + f->amplitude[i] * f->frequency[i] * cos(f->frequency[i] * t);
  }
}

static myb_status turning_orientation(void* context, myb_time t, myb_orientation* orientation) {
  double angles[3];
  double rates[3];
  turning_angles(context, myb_time_seconds(t), angles, rates);
  memset(orientation, 0, sizeof *orientation);
  myb_rotation_euler313(angles, rates, &orientation->rotation);
  orientation->ut = t;
  return MYB_OK;
}

static myb_status failing_orientation(void* context, myb_time t, myb_orientation* orientation) {
  (void)context;
  (void)t;
  (void)orientation;
  return MYB_E_RANGE;
}

/* The largest element of the difference of the rotations (and of their rates, times a day)
 * of the reader and of the provider over the span. */
static double compare(const myb_pck* pck, int32_t frame, int32_t base, turning* f, double start,
                      double stop, double* worst_rate) {
  double worst = 0.0;
  *worst_rate = 0.0;
  for (int k = 0; k < 400; ++k) {
    const double t = start + (stop - start) * ((double)k + 0.31) / 400.0;
    myb_rotation read;
    int32_t from = 0;
    MT_CHECK(myb_pck_rotation(pck, frame, myb_time_from_seconds(t), &read, &from) == MYB_OK);
    MT_CHECK(from == base);
    myb_orientation o;
    MT_CHECK(turning_orientation(f, myb_time_from_seconds(t), &o) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      worst = fmax(worst, fabs(read.m[i] - o.rotation.m[i]));
      *worst_rate = fmax(*worst_rate, 86400.0 * fabs(read.dm[i] - o.rotation.dm[i]));
    }
  }
  return worst;
}

static void test_pck(void) {
  const double start = -43200.0;
  const double stop = 3.0 * 86400.0;
  /* a planet: the pole fixed, a turn in 0.9 day */
  turning uniform = {.angles = {0.3, 0.4, 1.0}, .rates = {0.0, 0.0, 2.0 * pi / (0.9 * 86400.0)}};
  /* a moon: a turn in 27.3 days, the pole and the meridian librate */
  turning librating = {.angles = {0.05, 0.41, 0.7},
                       .rates = {-1.0e-8, 0.0, 2.0 * pi / (27.3 * 86400.0)},
                       .amplitude = {0.02, 0.01, 0.003},
                       .frequency = {2.0e-5, 4.1e-5, 6.3e-5}};
  myb_pck_writer* writer = NULL;
  MT_CHECK(myb_pck_writer_create("myb test", NULL, &writer) == MYB_OK);
  double error = -1.0;
  MT_CHECK(myb_pck_writer_add_chebyshev(writer, 31099, 1, "a planet", turning_orientation,
                                        &uniform, start, stop, 14, 3, &error) == MYB_OK);
  printf("PCK of a uniform rotation, degree 3: the writer reports %.1e rad\n", error);
  MT_CHECK(error >= 0.0 && error < 1e-14);
  double errors[3];
  static const int degrees[3] = {5, 9, 13};
  for (int k = 0; k < 3; ++k) {
    myb_pck_writer* trial = NULL;
    MT_CHECK(myb_pck_writer_create(NULL, NULL, &trial) == MYB_OK);
    MT_CHECK(myb_pck_writer_add_chebyshev(trial, 31098, 17, "a moon", turning_orientation,
                                          &librating, start, stop, 7, degrees[k],
                                          &errors[k]) == MYB_OK);
    myb_pck_writer_destroy(trial);
  }
  printf("PCK of a librating frame, half a day a part: degree 5 %.1e rad, 9 %.1e, 13 %.1e\n",
         errors[0], errors[1], errors[2]);
  MT_CHECK(errors[0] > 1e-7 && errors[1] < 1e-3 * errors[0] && errors[2] < 1e-12);
  MT_CHECK(myb_pck_writer_add_chebyshev(writer, 31098, 17, "a moon", turning_orientation,
                                        &librating, start, stop, 7, 13, &error) == MYB_OK);
  const size_t bytes = myb_pck_writer_size(writer);
  MT_CHECK(bytes > 0 && bytes % 1024 == 0);
  unsigned char* image = malloc(bytes);
  MT_CHECK(image != NULL);
  if (image == NULL) {
    myb_pck_writer_destroy(writer);
    return;
  }
  MT_CHECK(myb_pck_writer_write(writer, image, bytes - 1) == MYB_E_OVERFLOW);
  MT_CHECK(myb_pck_writer_write(writer, image, bytes) == MYB_OK);
  MT_CHECK(memcmp(image, "DAF/PCK ", 8) == 0);
  const char* dump = getenv("MYB_PCK_WRITER_DUMP");
  if (dump != NULL) {
    MT_CHECK(myb_pck_writer_write_file(writer, dump) == MYB_OK);
    printf("dump: %s\n", dump);
  }
  myb_pck* pck = NULL;
  myb_diag diag;
  MT_CHECK(myb_pck_create(NULL, &pck) == MYB_OK);
  MT_CHECK(myb_pck_load_memory(pck, image, bytes, &diag) == MYB_OK);
  double rate = 0.0;
  const double planet = compare(pck, 31099, 1, &uniform, start, stop, &rate);
  printf("read back: the uniform rotation within %.1e (rates, a day: %.1e)", planet, rate);
  MT_CHECK(planet < 1e-13 && rate < 1e-11);
  const double moon = compare(pck, 31098, 17, &librating, start, stop, &rate);
  printf(", the librating frame within %.1e (%.1e); the writer reported %.1e rad\n", moon, rate,
         error);
  MT_CHECK(moon < 1e-12 && rate < 1e-9 && error < 1e-12);
  /* outside the span, another frame */
  myb_rotation r;
  MT_CHECK(myb_pck_rotation(pck, 31099, myb_time_from_seconds(stop + 10.0), &r, NULL) ==
           MYB_E_RANGE);
  MT_CHECK(myb_pck_rotation(pck, 31000, myb_time_from_seconds(0.0), &r, NULL) == MYB_E_RANGE);

  /* against CSPICE: pxform_c from the base frame at the listed epochs */
  FILE* reference = fopen(MYB_TEST_DATA_DIR "/frames/pck_writer_cspice.ref", "r");
  MT_CHECK(reference != NULL);
  if (reference != NULL) {
    char line[1024];
    double worst = 0.0;
    int count = 0;
    /* bound: the lines of the file */
    while (fgets(line, sizeof line, reference) != NULL) {
      int frame = 0;
      double et = 0.0;
      double m[9];
      if (sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &frame, &et, &m[0], &m[1],
                 &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], &m[8]) != 11) {
        continue;
      }
      MT_CHECK(myb_pck_rotation(pck, frame, myb_time_from_seconds(et), &r, NULL) == MYB_OK);
      for (int i = 0; i < 9; ++i) {
        worst = fmax(worst, fabs(r.m[i] - m[i]));
      }
      ++count;
    }
    fclose(reference);
    printf("the file read by CSPICE (pxform, %d epochs): rotations within %.1e\n", count, worst);
    MT_CHECK(count == 40 && worst < 1e-13);
  }

  /* a segment copied word for word */
  myb_pck_writer* raw = NULL;
  MT_CHECK(myb_pck_writer_create(NULL, NULL, &raw) == MYB_OK);
  const double words[12] = {0.0,  43200.0, 0.3, 0.0, 0.4, 0.0,
                            1.0,  3.0,     -43200.0, 86400.0, 8.0, 1.0};
  MT_CHECK(myb_pck_writer_add_raw(raw, 31097, 1, "constant", 2, -43200.0, 43200.0, words, 12) ==
           MYB_OK);
  MT_CHECK(myb_pck_writer_add_raw(raw, 31097, 1, "constant", 0, -43200.0, 43200.0, words, 12) ==
           MYB_E_DOMAIN);
  const size_t raw_bytes = myb_pck_writer_size(raw);
  unsigned char* raw_image = malloc(raw_bytes);
  MT_CHECK(raw_image != NULL);
  if (raw_image != NULL) {
    MT_CHECK(myb_pck_writer_write(raw, raw_image, raw_bytes) == MYB_OK);
    myb_pck* small = NULL;
    MT_CHECK(myb_pck_create(NULL, &small) == MYB_OK);
    MT_CHECK(myb_pck_load_memory(small, raw_image, raw_bytes, &diag) == MYB_OK);
    double angles[3];
    double rates[3];
    MT_CHECK(myb_pck_euler(small, 31097, myb_time_from_seconds(21600.0), angles, rates, NULL) ==
             MYB_OK);
    /* 1 + 3 x at x = 0.5 */
    MT_CHECK(fabs(angles[0] - 0.3) < 1e-15 && fabs(angles[1] - 0.4) < 1e-15 &&
             fabs(angles[2] - 2.5) < 1e-15 && fabs(rates[2] - 3.0 / 43200.0) < 1e-18);
    myb_pck_destroy(small);
    free(raw_image);
  }
  myb_pck_writer_destroy(raw);

  /* what is refused: the pole on the z axis of the base frame, the arguments, the provider */
  turning upright = {.angles = {0.3, 0.0, 1.0}, .rates = {0.0, 0.0, 7.0e-5}};
  MT_CHECK(myb_pck_writer_add_chebyshev(writer, 31096, 1, "upright", turning_orientation,
                                        &upright, start, stop, 4, 5, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_pck_writer_add_chebyshev(writer, 31096, 1, "x", turning_orientation, &uniform,
                                        start, stop, 4, 0, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_pck_writer_add_chebyshev(writer, 31096, 1, "x", turning_orientation, &uniform,
                                        start, stop, 4, 28, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_pck_writer_add_chebyshev(writer, 31096, 1, "x", turning_orientation, &uniform,
                                        stop, start, 4, 5, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_pck_writer_add_chebyshev(writer, 31096, 1, "x", turning_orientation, &uniform,
                                        start, stop, 0, 5, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_pck_writer_add_chebyshev(writer, 31096, 1, NULL, turning_orientation, &uniform,
                                        start, stop, 4, 5, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_pck_writer_add_chebyshev(writer, 31096, 1, "x", failing_orientation, NULL, start,
                                        stop, 4, 5, NULL) == MYB_E_RANGE);
  MT_CHECK(myb_pck_writer_size(writer) == bytes); /* nothing was added */
  myb_pck_writer* none = NULL;
  MT_CHECK(myb_pck_writer_create("\x01", NULL, &none) == MYB_E_DOMAIN);
  myb_pck_writer_destroy(NULL);
  myb_pck_destroy(pck);
  free(image);
  myb_pck_writer_destroy(writer);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_raw_copy();
  test_pck();
  return mt_finish("test_pck_writer");
}
