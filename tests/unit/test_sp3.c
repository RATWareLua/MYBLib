/* SPDX-License-Identifier: BSD-3-Clause */
/* The SP3 writer (myb_sp3_write) against a file of the ILRS and against the reader:
 *  - the excerpt of the ILRS combined orbit of LAGEOS-1 read and written from what was read:
 *    every epoch, position and velocity line comes back character for character, and so do
 *    the two first lines of the header and the lines of satellites, accuracies and types;
 *  - SP3-d with positions only, a missing position, more than 85 satellites, an epoch that is
 *    not a midnight (week, seconds of the week, the fraction of the day in the header);
 *  - clocks, their rates, the accuracy codes of the header and of the records and the flags:
 *    a file in the layout of the format description read value by value, written from what
 *    was read and compared line by line;
 *  - invalid descriptions and a buffer too small. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/kin.h"

static char* read_all(const char* path, size_t* size) {
  FILE* f = fopen(path, "rb");
  if (f == NULL) {
    return NULL;
  }
  fseek(f, 0, SEEK_END);
  const long n = ftell(f);
  fseek(f, 0, SEEK_SET);
  char* text = malloc((size_t)n + 1);
  if (text != NULL && fread(text, 1, (size_t)n, f) != (size_t)n) {
    free(text);
    text = NULL;
  }
  if (text != NULL) {
    text[n] = '\0';
    *size = (size_t)n;
  }
  fclose(f);
  return text;
}

/* The start of line k (from 0) of a text. */
static const char* line_start(const char* text, int k) {
  const char* p = text;
  for (int i = 0; i < k && p != NULL; ++i) {
    p = strchr(p, '\n');
    if (p != NULL) {
      ++p;
    }
  }
  return p;
}

/* The lines of epochs and states: from the first epoch line to the end, without "EOF". */
static const char* data_of(const char* text, size_t* length) {
  const char* begin = strstr(text, "\n*  ");
  if (begin == NULL) {
    *length = 0;
    return NULL;
  }
  ++begin;
  size_t n = strlen(begin);
  if (n >= 4 && memcmp(begin + n - 4, "EOF\n", 4) == 0) {
    n -= 4;
  }
  *length = n;
  return begin;
}

static void test_excerpt(void) {
  const char* path = MYB_TEST_DATA_DIR "/sp3/ilrsa.lageos1.250607.excerpt.sp3";
  size_t size = 0;
  char* original = read_all(path, &size);
  MT_CHECK(original != NULL);
  if (original == NULL) {
    return;
  }
  myb_sp3* sp3 = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_sp3_load(original, size, NULL, &sp3, &diag) == MYB_OK);
  if (sp3 == NULL) {
    free(original);
    return;
  }
  const size_t epochs = myb_sp3_epoch_count(sp3);
  MT_CHECK(epochs == 1680 && myb_sp3_satellite_count(sp3) == 1);
  myb_calendar* labels = malloc(epochs * sizeof *labels);
  double* positions = malloc(3 * epochs * sizeof(double));
  double* velocities = malloc(3 * epochs * sizeof(double));
  MT_CHECK(labels != NULL && positions != NULL && velocities != NULL);
  if (labels == NULL || positions == NULL || velocities == NULL) {
    return;
  }
  for (size_t e = 0; e < epochs; ++e) {
    myb_sp3_epoch(sp3, e, &labels[e]);
    MT_CHECK(myb_sp3_state(sp3, 0, e, positions + 3 * e, velocities + 3 * e) == 3);
  }
  const char* const ids[1] = {myb_sp3_satellite(sp3, 0)};
  const myb_sp3_data data = {.version = 'c',
                             .file_type = 'L',
                             .time_system = myb_sp3_time_system(sp3),
                             .frame = myb_sp3_frame(sp3),
                             .orbit_type = "FIT",
                             .agency = "COMB",
                             .data_used = "  SLR",
                             .satellite_count = 1,
                             .satellites = ids,
                             .epoch_count = epochs,
                             .epochs = labels,
                             .interval = 120.0,
                             .positions = positions,
                             .velocities = velocities};
  size_t length = 0;
  MT_CHECK(myb_sp3_write(&data, NULL, 0, &length) == MYB_E_OVERFLOW && length > 100000);
  char* text = malloc(length + 1);
  MT_CHECK(text != NULL);
  if (text == NULL) {
    return;
  }
  size_t written = 0;
  MT_CHECK(myb_sp3_write(&data, text, length, &written) == MYB_E_OVERFLOW);
  MT_CHECK(myb_sp3_write(&data, text, length + 1, &written) == MYB_OK);
  MT_CHECK(written == length && strlen(text) == length);
  /* the header: the first line but for the number of epochs (the excerpt keeps the header of
   * the whole file), the second whole, the lines of satellites, accuracies and types */
  static const char first[] = "#cV2025  6  1  0  0  0.00000000    1680   SLR SLR20 FIT COMB\n";
  MT_CHECK(memcmp(text, first, sizeof first - 1) == 0);
  MT_CHECK(memcmp(original, first, 32) == 0 && memcmp(original + 39, first + 39, 22) == 0);
  const char* mine2 = line_start(text, 1);
  const char* mine18 = line_start(text, 18);
  const char* theirs2 = line_start(original, 1);
  const char* theirs18 = line_start(original, 18);
  MT_CHECK(mine2 != NULL && mine18 != NULL && theirs2 != NULL && theirs18 != NULL);
  if (mine2 != NULL && mine18 != NULL && theirs2 != NULL && theirs18 != NULL) {
    MT_CHECK(mine18 - mine2 == theirs18 - theirs2);
    MT_CHECK(memcmp(mine2, theirs2, (size_t)(theirs18 - theirs2)) == 0);
  }
  /* every epoch and state line */
  size_t mine_length = 0;
  size_t theirs_length = 0;
  const char* mine = data_of(text, &mine_length);
  const char* theirs = data_of(original, &theirs_length);
  MT_CHECK(mine != NULL && theirs != NULL && mine_length == theirs_length);
  if (mine != NULL && theirs != NULL && mine_length == theirs_length) {
    MT_CHECK(memcmp(mine, theirs, mine_length) == 0);
  }
  printf("SP3 excerpt: %zu characters of epochs and states written as in the file\n",
         mine_length);
  /* and the reader gives the same numbers from the written file */
  myb_sp3* back = NULL;
  MT_CHECK(myb_sp3_load(text, length, NULL, &back, &diag) == MYB_OK);
  if (back != NULL) {
    MT_CHECK(myb_sp3_epoch_count(back) == epochs && myb_sp3_satellite_count(back) == 1);
    MT_CHECK(strcmp(myb_sp3_time_system(back), "UTC") == 0);
    MT_CHECK(strcmp(myb_sp3_frame(back), "SLR20") == 0);
    MT_CHECK(strcmp(myb_sp3_satellite(back, 0), "L51") == 0);
    int same = 1;
    for (size_t e = 0; e < epochs; ++e) {
      myb_calendar label;
      double r[3];
      double v[3];
      myb_sp3_epoch(back, e, &label);
      same = same && myb_sp3_state(back, 0, e, r, v) == 3 &&
             memcmp(&label, &labels[e], sizeof label) == 0 &&
             memcmp(r, positions + 3 * e, sizeof r) == 0 &&
             memcmp(v, velocities + 3 * e, sizeof v) == 0;
    }
    MT_CHECK(same);
    myb_sp3_destroy(back);
  }
  free(text);
  free(labels);
  free(positions);
  free(velocities);
  myb_sp3_destroy(sp3);
  free(original);
}

static void test_small(void) {
  enum { SATELLITES = 100, EPOCHS = 3 };
  static char names[SATELLITES][4];
  static const char* ids[SATELLITES];
  static double positions[3 * EPOCHS * SATELLITES];
  for (int s = 0; s < SATELLITES; ++s) {
    snprintf(names[s], sizeof names[s], "%c%02d", s < 50 ? 'G' : 'E', s % 50 + 1);
    ids[s] = names[s];
  }
  for (int e = 0; e < EPOCHS; ++e) {
    for (int s = 0; s < SATELLITES; ++s) {
      double* r = positions + 3 * (e * SATELLITES + s);
      r[0] = 2.0e7 + 1234.5678912 * s + 0.123456789 * e;
      r[1] = -1.5e7 + 987.654321 * e - 3.3 * s;
      r[2] = 1.0e6 * (s % 7 - 3) + 0.000449 * (e + 1);
    }
  }
  positions[3 * (1 * SATELLITES + 4) + 1] = (double)NAN; /* epoch 1, satellite 4: missing */
  const myb_calendar labels[EPOCHS] = {{2025, 6, 1, 12, 0, 0.0},
                                       {2025, 6, 1, 12, 15, 0.0},
                                       {2025, 6, 1, 12, 30, 0.0}};
  myb_sp3_data data = {.version = 'd',
                       .file_type = 'M',
                       .time_system = "GPS",
                       .frame = "IGS20",
                       .orbit_type = "FIT",
                       .agency = "MYB",
                       .data_used = "ORBIT",
                       .satellite_count = SATELLITES,
                       .satellites = ids,
                       .epoch_count = EPOCHS,
                       .epochs = labels,
                       .interval = 900.0,
                       .positions = positions};
  static char text[40000];
  size_t length = 0;
  MT_CHECK(myb_sp3_write(&data, text, sizeof text, &length) == MYB_OK && strlen(text) == length);
  const size_t whole = length;
  MT_CHECK(memcmp(text, "#dP2025  6  1 12  0  0.00000000       3 ORBIT IGS20 FIT MYB \n", 61) == 0);
  /* 2025-06-01 is MJD 60827, the first day of GPS week 2369 */
  MT_CHECK(memcmp(line_start(text, 1),
                  "## 2369  43200.00000000   900.00000000 60827 0.5000000000000\n", 61) == 0);
  MT_CHECK(memcmp(line_start(text, 2), "+  100   G01G02G03", 18) == 0);
  /* six lines of 17 satellites, six of accuracies */
  MT_CHECK(memcmp(line_start(text, 7), "+        E36E37", 15) == 0);
  MT_CHECK(memcmp(line_start(text, 8), "++       ", 9) == 0);
  MT_CHECK(memcmp(line_start(text, 14), "%c M  cc GPS ccc", 16) == 0);
  myb_sp3* sp3 = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_sp3_load(text, length, NULL, &sp3, &diag) == MYB_OK);
  if (sp3 == NULL) {
    printf("line %ld: %s\n", diag.line, diag.message);
    return;
  }
  MT_CHECK(myb_sp3_satellite_count(sp3) == SATELLITES && myb_sp3_epoch_count(sp3) == EPOCHS);
  MT_CHECK(strcmp(myb_sp3_satellite(sp3, 99), "E50") == 0);
  MT_CHECK(strcmp(myb_sp3_time_system(sp3), "GPS") == 0);
  double worst = 0.0;
  for (size_t e = 0; e < EPOCHS; ++e) {
    myb_calendar label;
    myb_sp3_epoch(sp3, e, &label);
    MT_CHECK(memcmp(&label, &labels[e], sizeof label) == 0);
    for (size_t s = 0; s < SATELLITES; ++s) {
      double r[3];
      double v[3];
      const int flags = myb_sp3_state(sp3, s, e, r, v);
      if (e == 1 && s == 4) {
        MT_CHECK(flags == 0);
        continue;
      }
      MT_CHECK(flags == 1);
      for (int i = 0; i < 3; ++i) {
        worst = fmax(worst, fabs(r[i] - positions[3 * (e * SATELLITES + s) + (size_t)i]));
      }
    }
  }
  /* six decimals of a kilometre: half a millimetre */
  printf("SP3-d written and read back: positions within %.3e m\n", worst);
  MT_CHECK(worst <= 5.0e-4 * (1.0 + 1.0e-9));
  myb_sp3_destroy(sp3);

  /* invalid descriptions; a buffer too small */
  MT_CHECK(myb_sp3_write(&data, text, 100, &length) == MYB_E_OVERFLOW && length == whole);
  MT_CHECK(strlen(text) == 99);
  myb_sp3_data bad = data;
  bad.version = 'c'; /* at most 85 satellites */
  MT_CHECK(myb_sp3_write(&bad, text, sizeof text, &length) == MYB_E_DOMAIN);
  bad = data;
  bad.version = 'a';
  MT_CHECK(myb_sp3_write(&bad, text, sizeof text, &length) == MYB_E_DOMAIN);
  bad = data;
  bad.interval = 0.0;
  MT_CHECK(myb_sp3_write(&bad, text, sizeof text, &length) == MYB_E_DOMAIN);
  bad = data;
  bad.time_system = "GPST";
  MT_CHECK(myb_sp3_write(&bad, text, sizeof text, &length) == MYB_E_DOMAIN);
  bad = data;
  bad.epoch_count = 0;
  MT_CHECK(myb_sp3_write(&bad, text, sizeof text, &length) == MYB_E_DOMAIN);
  bad = data;
  bad.positions = NULL;
  MT_CHECK(myb_sp3_write(&bad, text, sizeof text, &length) == MYB_E_DOMAIN);
  static const char* const short_id[1] = {"G1"};
  bad = data;
  bad.satellite_count = 1;
  bad.satellites = short_id;
  MT_CHECK(myb_sp3_write(&bad, text, sizeof text, &length) == MYB_E_DOMAIN);
}

/* The records of the example of the format description (S. Hilla, SP3-c, 2010, and SP3-d,
 * 2016), cut to three satellites, with velocity records added in the same layout. */
static const char clocks_text[] =
    "#cV2001  8  8  0  0  0.00000000       2 ORBIT IGS97 HLM IGS \n"
    "## 1126 259200.00000000   900.00000000 52129 0.0000000000000\n"
    "+    3   G01G02G03  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "++         7  8  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "%c G  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n"
    "%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n"
    "%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n"
    "%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n"
    "%i    0    0    0    0      0      0      0      0         0\n"
    "%i    0    0    0    0      0      0      0      0         0\n"
    "/* written by myb\n"
    "/*\n"
    "/*\n"
    "/*\n"
    "*  2001  8  8  0  0  0.00000000\n"
    "PG01 -11044.805800 -10475.672350  21929.418200    189.163300 18 18 18 219\n"
    "VG01  20298.880364 -18462.044804   1381.387685     -4.534317 14 14 14 191\n"
    "PG02 -12593.593500  10170.327650 -20354.534400    -55.976000 18 18 18 219 EP  MP\n"
    "VG02  -9481.923808 -25832.652567  -7277.160056      8.801258 14 14 14 191\n"
    "PG03   9335.606450 -21952.990750 -11817.338450 999999.999999\n"
    "VG03  12497.392894  -8482.260298  25230.348435 999999.999999\n"
    "*  2001  8  8  0 15  0.00000000\n"
    "PG01 -11044.805800 -10475.672350  21929.418200    189.163300 18 18 18 219      P\n"
    "VG01  20298.880364 -18462.044804   1381.387685     -4.534317 14 14 14 191\n"
    "PG02 -12593.593500  10170.327650 -20354.534400    -55.976000 17    19     E\n"
    "VG02  -9481.923808 -25832.652567  -7277.160056      8.801258\n"
    "PG03      0.000000      0.000000      0.000000     54.862300          201\n"
    "VG03      0.000000      0.000000      0.000000 999999.999999\n"
    "EOF\n";

static void test_clocks(void) {
  enum { SATELLITES = 3, EPOCHS = 2, CELLS = SATELLITES * EPOCHS };
  myb_sp3* sp3 = NULL;
  myb_diag diag;
  MT_CHECK(myb_sp3_load(clocks_text, sizeof clocks_text - 1, NULL, &sp3, &diag) == MYB_OK);
  if (sp3 == NULL) {
    return;
  }
  MT_CHECK(myb_sp3_satellite_count(sp3) == SATELLITES && myb_sp3_epoch_count(sp3) == EPOCHS);
  /* the header: 2^7 and 2^8 mm, unknown */
  MT_CHECK(fabs(myb_sp3_accuracy(sp3, 0) - 0.128) < 1e-16);
  MT_CHECK(fabs(myb_sp3_accuracy(sp3, 1) - 0.256) < 1e-16);
  MT_CHECK(myb_sp3_accuracy(sp3, 2) == 0.0);
  /* G01 at the first epoch: 189.1633 us, -4.534317e-4 us/s; 1.25^18 mm, 1.025^219 ps,
   * 1.25^14 of 1e-4 mm/s, 1.025^191 of 1e-4 ps/s */
  double clock = 0.0;
  double rate = 0.0;
  MT_CHECK(myb_sp3_clock(sp3, 0, 0, &clock, &rate) == 3);
  MT_CHECK(fabs(clock - 189.1633e-6) < 1e-19 && fabs(rate + 4.534317e-10) < 1e-24);
  double sigma[8];
  myb_sp3_sigma(sp3, 0, 0, sigma);
  const double position_sigma = 1e-3 * pow(1.25, 18.0);
  const double clock_sigma = 1e-12 * pow(1.025, 219.0);
  const double velocity_sigma = 1e-7 * pow(1.25, 14.0);
  const double rate_sigma = 1e-16 * pow(1.025, 191.0);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(sigma[i] / position_sigma - 1.0) < 1e-14);
    MT_CHECK(fabs(sigma[4 + i] / velocity_sigma - 1.0) < 1e-14);
  }
  MT_CHECK(fabs(sigma[3] / clock_sigma - 1.0) < 1e-13 && fabs(sigma[7] / rate_sigma - 1.0) < 1e-13);
  printf("SP3 accuracy codes: position %.4f m, clock %.3e s, velocity %.3e m/s, rate %.3e s/s\n",
         sigma[0], sigma[3], sigma[4], sigma[7]);
  MT_CHECK(myb_sp3_flags(sp3, 0, 0) == 0);
  MT_CHECK(myb_sp3_flags(sp3, 1, 0) == (MYB_SP3_CLOCK_EVENT | MYB_SP3_CLOCK_PREDICTED |
                                        MYB_SP3_MANEUVER | MYB_SP3_ORBIT_PREDICTED));
  MT_CHECK(myb_sp3_flags(sp3, 0, 1) == MYB_SP3_ORBIT_PREDICTED);
  MT_CHECK(myb_sp3_flags(sp3, 1, 1) == MYB_SP3_CLOCK_EVENT);
  /* no clock; codes left blank; a clock without a position */
  MT_CHECK(myb_sp3_clock(sp3, 2, 0, &clock, &rate) == 0 && clock == 0.0 && rate == 0.0);
  MT_CHECK(myb_sp3_state(sp3, 2, 0, NULL, NULL) == 3);
  myb_sp3_sigma(sp3, 2, 0, sigma);
  for (int i = 0; i < 8; ++i) {
    MT_CHECK(sigma[i] == 0.0);
  }
  myb_sp3_sigma(sp3, 1, 1, sigma);
  MT_CHECK(sigma[0] > 0.0 && sigma[1] == 0.0 && sigma[2] > sigma[0] && sigma[3] == 0.0 &&
           sigma[4] == 0.0 && sigma[7] == 0.0);
  MT_CHECK(myb_sp3_state(sp3, 2, 1, NULL, NULL) == 0);
  MT_CHECK(myb_sp3_clock(sp3, 2, 1, &clock, NULL) == 1 && fabs(clock - 54.8623e-6) < 1e-19);
  myb_sp3_sigma(sp3, 2, 1, sigma);
  MT_CHECK(sigma[0] == 0.0 && fabs(sigma[3] / (1e-12 * pow(1.025, 201.0)) - 1.0) < 1e-13);

  /* written from what was read: the same lines */
  static const char* ids[SATELLITES] = {"G01", "G02", "G03"};
  myb_calendar labels[EPOCHS];
  double positions[3 * CELLS];
  double velocities[3 * CELLS];
  double clocks[CELLS];
  double rates[CELLS];
  double sigmas[8 * CELLS];
  double accuracy[SATELLITES];
  unsigned char flags[CELLS];
  for (size_t e = 0; e < EPOCHS; ++e) {
    myb_sp3_epoch(sp3, e, &labels[e]);
    for (size_t s = 0; s < SATELLITES; ++s) {
      const size_t cell = e * SATELLITES + s;
      const int have = myb_sp3_state(sp3, s, e, positions + 3 * cell, velocities + 3 * cell);
      if (!(have & 1)) {
        positions[3 * cell] = (double)NAN;
      }
      if (!(have & 2)) {
        velocities[3 * cell] = (double)NAN;
      }
      const int timed = myb_sp3_clock(sp3, s, e, &clocks[cell], &rates[cell]);
      if (!(timed & 1)) {
        clocks[cell] = (double)NAN;
      }
      if (!(timed & 2)) {
        rates[cell] = (double)NAN;
      }
      myb_sp3_sigma(sp3, s, e, sigmas + 8 * cell);
      flags[cell] = (unsigned char)myb_sp3_flags(sp3, s, e);
    }
  }
  for (size_t s = 0; s < SATELLITES; ++s) {
    accuracy[s] = myb_sp3_accuracy(sp3, s);
  }
  myb_sp3_data data = {.version = 'c',
                       .file_type = 'G',
                       .time_system = "GPS",
                       .frame = "IGS97",
                       .orbit_type = "HLM",
                       .agency = "IGS",
                       .data_used = "ORBIT",
                       .satellite_count = SATELLITES,
                       .satellites = ids,
                       .epoch_count = EPOCHS,
                       .epochs = labels,
                       .interval = 900.0,
                       .positions = positions,
                       .velocities = velocities,
                       .clocks = clocks,
                       .clock_rates = rates,
                       .accuracy = accuracy,
                       .sigmas = sigmas,
                       .flags = flags};
  static char text[8192];
  size_t length = 0;
  MT_CHECK(myb_sp3_write(&data, text, sizeof text, &length) == MYB_OK);
  int different = 0;
  int lines = 0;
  const char* a = clocks_text;
  const char* b = text;
  /* bound: the lines of the text */
  while (*a != '\0' && *b != '\0') {
    const char* ea = strchr(a, '\n');
    const char* eb = strchr(b, '\n');
    if (ea == NULL || eb == NULL) {
      break;
    }
    if (ea - a != eb - b || memcmp(a, b, (size_t)(ea - a)) != 0) {
      if (different++ < 5) {
        printf("line %d:\n  file    %.*s\n  written %.*s\n", lines + 1, (int)(ea - a), a,
               (int)(eb - b), b);
      }
    }
    ++lines;
    a = ea + 1;
    b = eb + 1;
  }
  printf("SP3 with clocks, accuracy codes and flags: %d lines written, %d differ from the file\n",
         lines, different);
  MT_CHECK(lines == 37 && different == 0 && *a == '\0' && *b == '\0');

  /* a sigma between two codes goes to the nearer; one too large for the field to its end */
  sigmas[0] = 1e-3 * pow(1.25, 18.4);
  sigmas[1] = 1e-3 * pow(1.25, 18.6);
  sigmas[2] = 1.0e9;
  sigmas[3] = 0.4e-12;
  accuracy[2] = 0.0005; /* less than the least the header can state: 2 mm */
  MT_CHECK(myb_sp3_write(&data, text, sizeof text, &length) == MYB_OK);
  myb_sp3* again = NULL;
  MT_CHECK(myb_sp3_load(text, length, NULL, &again, &diag) == MYB_OK);
  if (again != NULL) {
    myb_sp3_sigma(again, 0, 0, sigma);
    MT_CHECK(fabs(sigma[0] / (1e-3 * pow(1.25, 18.0)) - 1.0) < 1e-13);
    MT_CHECK(fabs(sigma[1] / (1e-3 * pow(1.25, 19.0)) - 1.0) < 1e-13);
    MT_CHECK(fabs(sigma[2] / (1e-3 * pow(1.25, 99.0)) - 1.0) < 1e-13);
    MT_CHECK(fabs(sigma[3] / 1e-12 - 1.0) < 1e-13);
    MT_CHECK(fabs(myb_sp3_accuracy(again, 2) - 0.002) < 1e-17);
    myb_sp3_destroy(again);
  }
  /* what the fields cannot hold */
  clocks[0] = 1.0;
  MT_CHECK(myb_sp3_write(&data, text, sizeof text, &length) == MYB_E_DOMAIN);
  clocks[0] = 1e-4;
  rates[0] = 1e-4;
  MT_CHECK(myb_sp3_write(&data, text, sizeof text, &length) == MYB_E_DOMAIN);
  rates[0] = 0.0;
  data.velocities = NULL; /* rates need the velocity records */
  MT_CHECK(myb_sp3_write(&data, text, sizeof text, &length) == MYB_E_DOMAIN);
  myb_sp3_destroy(sp3);
}

int main(void) {
  test_excerpt();
  test_small();
  test_clocks();
  return mt_finish("test_sp3");
}
