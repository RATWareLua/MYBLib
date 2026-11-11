/* SPDX-License-Identifier: BSD-3-Clause */
/* SGP4 (ext/sgp4) against the verification output of Vallado et al. (2006): the 33 element sets
 * of SGP4-VER.TLE propagated over their time spans as the verification program does (WGS-72,
 * improved mode), compared line by line with tcppver.out (positions printed to 1e-8 km,
 * velocities to 1e-9 km/s), and the sequence of error codes [1, 1, 6, 6, 4, 3, 6]. Then the TLE
 * parser (Alpha-5, checksums, epochs), the modes and the error statuses; the element sets
 * written back (every set of the verification file reads back to the same values) and the
 * mean-elements message of the CCSDS turned into an element set and back. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb_sgp4.h"

#define MAX_LINES 800

typedef struct out_line {
  int header;    /* "NNNNN xx" */
  long satnum;
  double f[7];   /* tsince, r (km), v (km/s) */
} out_line;

static char* read_all(const char* path) {
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
  }
  fclose(f);
  return text;
}

/* Splits text into lines in place (LF or CRLF); returns the count. */
static int split_lines(char* text, char** lines, int max) {
  int n = 0;
  char* p = text;
  while (*p != '\0' && n < max) {
    lines[n++] = p;
    char* e = strchr(p, '\n');
    if (e == NULL) {
      break;
    }
    if (e > p && e[-1] == '\r') {
      e[-1] = '\0';
    }
    *e = '\0';
    p = e + 1;
  }
  return n;
}

static out_line out[MAX_LINES];
static int out_count = 0;
static int errors[16];
static int error_count = 0;
static double previous[7];

static void emit(double tsince, const double r[3], const double v[3]) {
  if (out_count < MAX_LINES) {
    out[out_count] = (out_line){.header = 0, .f = {tsince, r[0], r[1], r[2], v[0], v[1], v[2]}};
    ++out_count;
  }
  memcpy(previous, out[out_count - 1].f, sizeof previous);
}

static int invoke(const myb_sgp4* s, double t, double r[3], double v[3]) {
  myb_sgp4_error e = MYB_SGP4_OK;
  const myb_status st = myb_sgp4_propagate_km(s, t, r, v, &e);
  MT_CHECK((st == MYB_OK) == (e == MYB_SGP4_OK || e == MYB_SGP4_DECAYED));
  return (int)e;
}

static void run_verification(void) {
  char* tle_text = read_all(MYB_TEST_DATA_DIR "/sgp4/SGP4-VER.TLE");
  char* ref_text = read_all(MYB_TEST_DATA_DIR "/sgp4/tcppver.out");
  MT_CHECK(tle_text != NULL && ref_text != NULL);
  if (tle_text == NULL || ref_text == NULL) {
    free(tle_text);
    free(ref_text);
    return;
  }
  static char* tle_lines[400];
  static char* ref_lines[MAX_LINES];
  const int nt = split_lines(tle_text, tle_lines, 400);
  const int nr = split_lines(ref_text, ref_lines, MAX_LINES);

  int satellites = 0;
  for (int i = 0; i + 1 < nt; ++i) {
    if (tle_lines[i][0] != '1') {
      continue;
    }
    const char* line1 = tle_lines[i];
    const char* line2 = tle_lines[++i];
    myb_tle tle;
    myb_diag diag = {0};
    MT_CHECK(myb_tle_parse(line1, line2, MYB_TLE_NO_CHECKSUM, &tle, &diag) == MYB_OK);
    if (diag.message[0] != '\0') {
      printf("%s\n", diag.message);
    }
    myb_sgp4 s;
    MT_CHECK(myb_sgp4_init(&tle, MYB_SGP4_WGS72, MYB_SGP4_IMPROVED, &s) == MYB_OK);
    ++satellites;
    out[out_count++] = (out_line){.header = 1, .satnum = tle.catalog_number};
    double r[3];
    double v[3];
    int e = invoke(&s, 0.0, r, v);
    if (e != 0 && e != 6) {
      errors[error_count++] = e;
      double r0[3] = {previous[1], previous[2], previous[3]};
      double v0[3] = {previous[4], previous[5], previous[6]};
      emit(0.0, r0, v0); /* "(Use previous data line)" */
      continue;
    }
    emit(0.0, r, v);
    double tstart = 0.0;
    double tend = 0.0;
    double tstep = 0.0;
    MT_CHECK(sscanf(line2 + 69, "%lf %lf %lf", &tstart, &tend, &tstep) == 3);
    double tsince = tstart;
    int stopped = 0;
    for (int k = 0; k < 10000 && tsince <= tend; ++k) {
      if (tsince == tstart && tstart == 0.0) {
        tsince += tstep;
        continue;
      }
      e = invoke(&s, tsince, r, v);
      if (e != 0) {
        errors[error_count++] = e;
        stopped = 1;
        break;
      }
      emit(tsince, r, v);
      tsince += tstep;
    }
    if (!stopped && tsince - tend < tstep - 1e-6) {
      e = invoke(&s, tend, r, v);
      if (e != 0) {
        errors[error_count++] = e;
      } else {
        emit(tend, r, v);
      }
    }
  }
  MT_CHECK(satellites == 33);

  /* the reference */
  MT_CHECK(nr == 700 && out_count == 700);
  double max_r = 0.0;
  double max_v = 0.0;
  double max_t = 0.0;
  double max_long = 0.0;
  int mismatched = 0;
  for (int i = 0; i < nr && i < out_count; ++i) {
    const char* line = ref_lines[i];
    if (strstr(line, "xx") != NULL) {
      long satnum = -1;
      MT_CHECK(sscanf(line, "%ld", &satnum) == 1);
      if (!out[i].header || out[i].satnum != satnum) {
        ++mismatched;
      }
      continue;
    }
    double f[7];
    if (out[i].header || sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &f[0], &f[1], &f[2], &f[3],
                                &f[4], &f[5], &f[6]) != 7) {
      ++mismatched;
      continue;
    }
    max_t = fmax(max_t, fabs(out[i].f[0] - f[0]));
    if (getenv("MYB_SGP4_DUMP") != NULL) {
      double dr = 0.0;
      for (int k = 1; k <= 3; ++k) {
        dr = fmax(dr, fabs(out[i].f[k] - f[k]));
      }
      if (dr > 1e-8 || out[i].f[0] != f[0]) {
        printf("line %d: t %.8f / %.8f dr %.3g\n", i + 1, out[i].f[0], f[0], dr);
      }
    }
    const int long_span = fabs(f[0]) > 1e5;
    for (int k = 1; k <= 3; ++k) {
      if (long_span) {
        max_long = fmax(max_long, fabs(out[i].f[k] - f[k]));
      } else {
        max_r = fmax(max_r, fabs(out[i].f[k] - f[k]));
      }
      max_v = fmax(max_v, fabs(out[i].f[k + 3] - f[k + 3]));
    }
  }
  printf("tcppver.out: 700 lines, max |dr| = %.3g km (%.3g km beyond 1e5 min), max |dv| = %.3g "
         "km/s\n",
         max_r, max_long, max_v);
  MT_CHECK(mismatched == 0);
  MT_CHECK(max_t <= 5e-9); /* times accumulate the steps; printed to 1e-8 min */
  /* the printed precision: 5e-9 km and 5e-10 km/s (half a unit of the last digit), plus the
   * rounding of the reference's own libm on 3.5 years of the deep-space SL-12 R/B (e = 0.79,
   * relative 5e-12; python-sgp4 differs there by the same 1.2e-7 km) */
  MT_CHECK(max_r <= 1.0e-8);
  MT_CHECK(max_long <= 2.0e-7);
  MT_CHECK(max_v <= 1.0e-9);
  static const int expected_errors[7] = {1, 1, 6, 6, 4, 3, 6};
  MT_CHECK(error_count == 7);
  for (int i = 0; i < 7 && i < error_count; ++i) {
    MT_CHECK(errors[i] == expected_errors[i]);
  }
  free(tle_text);
  free(ref_text);
}

static void test_parser(void) {
  /* Vanguard 1, with checksums */
  const char* l1 = "1 00005U 58002B   20287.20333880 -.00000016  00000-0 -22483-4 0  9998";
  const char* l2 = "2 00005  34.2443 225.5254 1845686 162.2516 205.2356 10.84869164218149";
  myb_tle t;
  myb_diag diag = {0};
  MT_CHECK(myb_tle_parse(l1, l2, 0, &t, &diag) == MYB_OK);
  MT_CHECK(t.catalog_number == 5 && strcmp(t.satnum, "00005") == 0);
  MT_CHECK(t.classification == 'U' && strcmp(t.designator, "58002B") == 0);
  MT_CHECK(t.epoch_year == 2020 && t.epoch_day == 287.20333880);
  MT_CHECK(t.ndot == -0.00000016 && t.nddot == 0.0);
  MT_CHECK(t.bstar == -0.22483 * 1e-4);
  MT_CHECK(t.element_number == 999 && t.revolution == 21814);
  MT_CHECK(t.inclination == 34.2443 && t.node == 225.5254 && t.eccentricity == 0.1845686);
  MT_CHECK(t.argument_of_perigee == 162.2516 && t.mean_anomaly == 205.2356);
  MT_CHECK(t.mean_motion == 10.84869164);
  myb_calendar c;
  myb_tle_epoch(&t, &c);
  /* day 287 of 2020 (leap year) = 13 October; 0.20333880 d = 4 h 52 min 48.472 s */
  MT_CHECK(c.year == 2020 && c.month == 10 && c.day == 13 && c.hour == 4 && c.minute == 52);
  MT_CHECK(fabs(c.second - 48.4723200) < 1e-6);

  /* a wrong checksum, a wrong line, differing numbers */
  char bad1[80];
  strcpy(bad1, l1);
  bad1[68] = '7';
  MT_CHECK(myb_tle_parse(bad1, l2, 0, &t, &diag) == MYB_E_FORMAT);
  MT_CHECK(diag.message[0] != '\0');
  MT_CHECK(myb_tle_parse(l2, l1, 0, &t, &diag) == MYB_E_FORMAT);
  char bad2[80];
  strcpy(bad2, l2);
  bad2[6] = '6';
  bad2[68] = '\0'; /* no checksum */
  MT_CHECK(myb_tle_parse(l1, bad2, 0, &t, &diag) == MYB_E_FORMAT);

  /* Alpha-5 (A0001 = 100001, Z9999 = 339999) */
  char a1[80];
  char a2[80];
  strcpy(a1, l1);
  strcpy(a2, l2);
  a1[68] = a2[68] = '\0';
  memcpy(a1 + 2, "A0001", 5);
  memcpy(a2 + 2, "A0001", 5);
  MT_CHECK(myb_tle_parse(a1, a2, 0, &t, &diag) == MYB_OK && t.catalog_number == 100001);
  memcpy(a1 + 2, "Z9999", 5);
  memcpy(a2 + 2, "Z9999", 5);
  MT_CHECK(myb_tle_parse(a1, a2, 0, &t, &diag) == MYB_OK && t.catalog_number == 339999);
  memcpy(a1 + 2, "J1234", 5);
  memcpy(a2 + 2, "J1234", 5);
  MT_CHECK(myb_tle_parse(a1, a2, 0, &t, &diag) == MYB_OK && t.catalog_number == 181234);
}

static void test_modes(void) {
  const char* l1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
  const char* l2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
  myb_tle t;
  MT_CHECK(myb_tle_parse(l1, l2, 0, &t, NULL) == MYB_OK);
  myb_sgp4 a;
  myb_sgp4 b;
  myb_sgp4 c;
  MT_CHECK(myb_sgp4_init(&t, MYB_SGP4_WGS72, MYB_SGP4_IMPROVED, &a) == MYB_OK);
  MT_CHECK(myb_sgp4_init(&t, MYB_SGP4_WGS72, MYB_SGP4_AFSPC, &b) == MYB_OK);
  MT_CHECK(myb_sgp4_init(&t, MYB_SGP4_WGS84, MYB_SGP4_IMPROVED, &c) == MYB_OK);
  double ra[3], va[3], rb[3], vb[3], rc[3], vc[3], rk[3], vk[3];
  MT_CHECK(myb_sgp4_propagate(&a, 4320.0, ra, va, NULL) == MYB_OK);
  MT_CHECK(myb_sgp4_propagate(&b, 4320.0, rb, vb, NULL) == MYB_OK);
  MT_CHECK(myb_sgp4_propagate(&c, 4320.0, rc, vc, NULL) == MYB_OK);
  MT_CHECK(myb_sgp4_propagate_km(&a, 4320.0, rk, vk, NULL) == MYB_OK);
  double dab = 0.0;
  double dac = 0.0;
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(ra[i] == 1000.0 * rk[i] && va[i] == 1000.0 * vk[i]);
    dab = fmax(dab, fabs(ra[i] - rb[i]));
    dac = fmax(dac, fabs(ra[i] - rc[i]));
  }
  /* near-Earth: the modes differ only in the sidereal time (unused); the gravity models by
   * metres to kilometres */
  MT_CHECK(dab == 0.0);
  MT_CHECK(dac > 1.0 && dac < 1e5);

  /* purity: the deep-space integrator gives the same bits whatever was computed before */
  const char* d1 = "1 09880U 77021A   06176.56157475  .00000421  00000-0  10000-3 0  9814";
  const char* d2 = "2 09880  64.5968 349.3786 7069051 270.0229  16.3320  2.00813614112380";
  MT_CHECK(myb_tle_parse(d1, d2, 0, &t, NULL) == MYB_OK);
  MT_CHECK(myb_sgp4_init(&t, MYB_SGP4_WGS72, MYB_SGP4_IMPROVED, &a) == MYB_OK);
  MT_CHECK(a.deep == 1 && a.irez == 2);
  MT_CHECK(myb_sgp4_propagate(&a, 2880.0, ra, va, NULL) == MYB_OK);
  MT_CHECK(myb_sgp4_propagate(&a, -1440.0, rb, vb, NULL) == MYB_OK);
  MT_CHECK(myb_sgp4_propagate(&a, 2880.0, rc, vc, NULL) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK_SAME(rc[i], ra[i]);
    MT_CHECK_SAME(vc[i], va[i]);
  }

  /* errors */
  myb_sgp4_error e = MYB_SGP4_OK;
  MT_CHECK(myb_sgp4_propagate(&a, (double)NAN, ra, va, &e) == MYB_E_DOMAIN);
  myb_sgp4_constants wgs72;
  myb_sgp4_constants_of(MYB_SGP4_WGS72, &wgs72);
  MT_CHECK(myb_sgp4_init_elements(&wgs72, MYB_SGP4_IMPROVED, 20000.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                                  0.05, 0.0, &a) == MYB_E_DOMAIN);
  MT_CHECK(myb_sgp4_init_elements(&wgs72, MYB_SGP4_IMPROVED, 20000.0, 0.0, 0.1, 0.0, 0.0, 0.0,
                                  -0.05, 0.0, &a) == MYB_E_DOMAIN);
  MT_CHECK(fabs(wgs72.xke - 0.07436691613317342) < 1e-17);

  /* GMST of 1982 at J2000.0 UT1: 280.46061837504 deg */
  MT_CHECK(fabs(myb_sgp4_gmst(2451545.0, 0.0) - 280.46061837504 * 3.14159265358979323846 / 180.0) <
           1e-12);
}

/* ---- SPK type 10 against CSPICE N0067 ---------------------------------------------------- */

static int parse_hex_line(const char* line, char* kind, long* body, int* set, double* values,
                          int max) {
  int used = 0;
  if (sscanf(line, "%15s %ld %d%n", kind, body, set, &used) != 3) {
    return -1;
  }
  const char* p = line + used;
  int n = 0;
  while (n < max) {
    char* end = NULL;
    const double v = strtod(p, &end);
    if (end == p) {
      break;
    }
    values[n++] = v;
    p = end;
  }
  return n;
}

/* The epoch of a set as SPICE gives it to SGP4: TDB -> UTC Julian date in one double. */
static double spice_epoch(const myb_timesys* sys, myb_scale tdb, double et) {
  myb_calendar utc;
  MT_CHECK(myb_timesys_to_utc(sys, myb_time_from_seconds(et), tdb, &utc) == MYB_OK);
  myb_time t;
  myb_calendar midnight = utc;
  midnight.hour = 12;
  midnight.minute = 0;
  midnight.second = 0.0;
  MT_CHECK(myb_calendar_to_time(&midnight, &t) == MYB_OK);
  const double days = (double)(t.s / 86400); /* whole days from 2000 January 1.5 */
  const double seconds = days * 86400.0 - 43200.0 + (utc.hour * 3600.0 + utc.minute * 60.0 +
                                                     utc.second);
  return 2451545.0 + seconds / 86400.0 - 2433281.5;
}

static void test_spk10(void) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const char* kernels[] = {MYB_DATA_DIR "/time/time-scales.tk", MYB_DATA_DIR "/time/tdb-tt-fb1990.tk",
                           MYB_DATA_DIR "/time/utc-tai.tk", MYB_DATA_DIR "/frames/iau1980.tk"};
  for (int i = 0; i < 4; ++i) {
    MT_CHECK(myb_kpool_load_file(pool, kernels[i], &diag) == MYB_OK);
  }
  myb_timesys* sys = NULL;
  myb_equinox* eq = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, &diag) == MYB_OK);
  MT_CHECK(myb_equinox_create(pool, "IAU1980", NULL, &eq, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_spk* spk = NULL;
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  MT_CHECK(myb_spk_load_file(spk, MYB_TEST_DATA_DIR "/spk/t10.bsp", &diag) == MYB_OK);
  char* text = read_all(MYB_TEST_DATA_DIR "/spk/t10_reference.txt");
  MT_CHECK(text != NULL && sys != NULL && eq != NULL && spk != NULL);
  if (text == NULL || sys == NULL || eq == NULL || spk == NULL) {
    free(text);
    myb_spk_destroy(spk);
    myb_equinox_destroy(eq);
    myb_timesys_destroy(sys);
    return;
  }
  MT_CHECK(myb_spk_segment_count(spk) == 2);
  myb_scale tdb;
  MT_CHECK(myb_timesys_scale(sys, "TDB", &tdb) == MYB_OK);

  static char* lines[200];
  const int n = split_lines(text, lines, 200);
  double geophs[8] = {0};
  double elems[2][3][10] = {{{0}}};
  double worst_teme[2] = {0.0, 0.0}; /* near-Earth, deep space (km) */
  double worst_spk[2] = {0.0, 0.0};
  double worst_spk_v[2] = {0.0, 0.0};
  int evsgp4 = 0;
  int states = 0;
  for (int i = 0; i < n; ++i) {
    if (strncmp(lines[i], "# GEOPHS", 8) == 0) {
      const char* p = lines[i] + 8;
      for (int k = 0; k < 8; ++k) {
        char* end = NULL;
        geophs[k] = strtod(p, &end);
        p = end;
      }
      continue;
    }
    if (lines[i][0] == '#' || lines[i][0] == '\0') {
      continue;
    }
    char kind[16];
    long body = 0;
    int set = 0;
    double v[12];
    const int count = parse_hex_line(lines[i], kind, &body, &set, v, 12);
    const int which = body == -900005 ? 0 : 1;
    if (strcmp(kind, "ELEMS") == 0 && count == 11 && set >= 0 && set < 3) {
      memcpy(elems[which][set], v + 1, sizeof elems[which][set]);
    } else if (strcmp(kind, "EVSGP4") == 0 && count == 7) {
      const double* e = elems[which][set];
      const myb_sgp4_constants c = {
          .radius = geophs[6], .xke = geophs[3], .j2 = geophs[0], .j3 = geophs[1], .j4 = geophs[2]};
      myb_sgp4 s;
      MT_CHECK(myb_sgp4_init_elements(&c, MYB_SGP4_AFSPC, spice_epoch(sys, tdb, e[9]), e[2], e[5],
                                      e[6], e[3], e[7], e[8], e[4], &s) == MYB_OK);
      double r[3];
      double vel[3];
      MT_CHECK(myb_sgp4_propagate_km(&s, (v[0] - e[9]) / 60.0, r, vel, NULL) == MYB_OK);
      for (int k = 0; k < 3; ++k) {
        worst_teme[which] = fmax(worst_teme[which], fabs(r[k] - v[1 + k]));
      }
      ++evsgp4;
    } else if (strcmp(kind, "SPK") == 0 && count == 7) {
      double st[6];
      MT_CHECK(myb_sgp4_spk10_state(spk, (size_t)which, sys, eq, myb_time_from_seconds(v[0]),
                                    st) == MYB_OK);
      for (int k = 0; k < 3; ++k) {
        worst_spk[which] = fmax(worst_spk[which], fabs(st[k] / 1000.0 - v[1 + k]));
        worst_spk_v[which] = fmax(worst_spk_v[which], fabs(st[3 + k] / 1000.0 - v[4 + k]));
      }
      ++states;
    }
  }
  printf("SPK type 10 vs CSPICE N0067: EVSGP4 (%d) near-Earth %.3g km, deep space %.3g km; "
         "states (%d) near-Earth %.3g km %.3g km/s, deep space %.3g km %.3g km/s\n",
         evsgp4, worst_teme[0], worst_teme[1], states, worst_spk[0], worst_spk_v[0], worst_spk[1],
         worst_spk_v[1]);
  MT_CHECK(evsgp4 == 20 && states == 79);
  /* measured: EVSGP4 4e-11 / 4.4e-9 km, states 2.7e-12 km 9e-14 km/s / 1.9e-9 km 1e-12 km/s.
   * Deep space feels the epoch (lunar and solar terms): SPICE's TDB -> UTC uses its one-term
   * TDB - TT, myb the Fairhead-Bretagnon series; both round the Julian date to 40 us. */
  MT_CHECK(worst_teme[0] < 1e-10);
  MT_CHECK(worst_teme[1] < 1e-8);
  MT_CHECK(worst_spk[0] < 1e-10 && worst_spk_v[0] < 1e-12);
  MT_CHECK(worst_spk[1] < 1e-8 && worst_spk_v[1] < 1e-11);

  /* the library's own SPK reader does not evaluate type 10; other types are refused here */
  double st[6];
  MT_CHECK(myb_spk_state(spk, -900005, 399, myb_time_from_seconds(elems[0][0][9]), st, NULL) ==
           MYB_E_FORMAT);
  myb_spk* other = NULL;
  MT_CHECK(myb_spk_create(NULL, &other) == MYB_OK);
  MT_CHECK(myb_spk_load_file(other, MYB_TEST_DATA_DIR "/spk/t01.bsp", &diag) == MYB_OK);
  MT_CHECK(myb_sgp4_spk10_state(other, 0, sys, eq, myb_time_from_seconds(0.0), st) ==
           MYB_E_FORMAT);
  myb_spk_destroy(other);
  free(text);
  myb_spk_destroy(spk);
  myb_equinox_destroy(eq);
  myb_timesys_destroy(sys);
}

static int same_tle(const myb_tle* a, const myb_tle* b) {
  return strcmp(a->satnum, b->satnum) == 0 && a->catalog_number == b->catalog_number &&
         a->classification == b->classification && strcmp(a->designator, b->designator) == 0 &&
         a->epoch_year == b->epoch_year && a->epoch_day == b->epoch_day && a->ndot == b->ndot &&
         a->nddot == b->nddot && a->bstar == b->bstar && a->ephemeris_type == b->ephemeris_type &&
         a->element_number == b->element_number && a->inclination == b->inclination &&
         a->node == b->node && a->eccentricity == b->eccentricity &&
         a->argument_of_perigee == b->argument_of_perigee && a->mean_anomaly == b->mean_anomaly &&
         a->mean_motion == b->mean_motion && a->revolution == b->revolution;
}

static void test_format(void) {
  /* a set of the catalogue comes back character for character */
  const char* l1 = "1 00005U 58002B   20287.20333880 -.00000016  00000-0 -22483-4 0  9998";
  const char* l2 = "2 00005  34.2443 225.5254 1845686 162.2516 205.2356 10.84869164218149";
  myb_tle t;
  myb_tle back;
  myb_diag diag = {0};
  char f1[70];
  char f2[70];
  MT_CHECK(myb_tle_parse(l1, l2, 0, &t, &diag) == MYB_OK);
  MT_CHECK(myb_tle_format(&t, f1, f2) == MYB_OK);
  MT_CHECK(strcmp(f1, l1) == 0 && strcmp(f2, l2) == 0);

  /* every set of the verification file: written from what was read, read again with its
   * checksum, the same values to the bit */
  char* text = read_all(MYB_TEST_DATA_DIR "/sgp4/SGP4-VER.TLE");
  MT_CHECK(text != NULL);
  if (text == NULL) {
    return;
  }
  static char* lines[400];
  const int n = split_lines(text, lines, 400);
  int sets = 0;
  int identical = 0;
  for (int i = 0; i + 1 < n; ++i) {
    if (lines[i][0] != '1') {
      continue;
    }
    const char* line1 = lines[i];
    const char* line2 = lines[++i];
    MT_CHECK(myb_tle_parse(line1, line2, MYB_TLE_NO_CHECKSUM, &t, &diag) == MYB_OK);
    const myb_status status = myb_tle_format(&t, f1, f2);
    MT_CHECK(status == MYB_OK);
    if (status != MYB_OK) {
      printf("not written: %s\n", line1);
      continue;
    }
    MT_CHECK(strlen(f1) == 69 && strlen(f2) == 69);
    MT_CHECK(myb_tle_parse(f1, f2, 0, &back, &diag) == MYB_OK);
    MT_CHECK(same_tle(&t, &back));
    ++sets;
    if (memcmp(f1, line1, 68) == 0 && memcmp(f2, line2, 68) == 0) {
      ++identical;
    } else if (getenv("MYB_TLE_SHOW") != NULL) {
      printf("%.68s\n%.68s\n%.68s\n%.68s\n\n", line1, f1, line2, f2);
    }
  }
  free(text);
  printf("element sets written back: %d, character for character: %d\n", sets, identical);
  /* two sets of the file are not as the catalogue writes them: one has no ephemeris type in
   * its column, one a zero drag term as 00000-0 */
  MT_CHECK(sets == 33 && identical == 31);

  /* Alpha-5 from the number; what does not fit a field */
  MT_CHECK(myb_tle_parse(l1, l2, 0, &t, &diag) == MYB_OK);
  t.satnum[0] = '\0';
  t.catalog_number = 181234;
  MT_CHECK(myb_tle_format(&t, f1, f2) == MYB_OK);
  MT_CHECK(memcmp(f1 + 2, "J1234", 5) == 0 && memcmp(f2 + 2, "J1234", 5) == 0);
  MT_CHECK(myb_tle_parse(f1, f2, 0, &back, &diag) == MYB_OK && back.catalog_number == 181234);
  t.catalog_number = 340000;
  MT_CHECK(myb_tle_format(&t, f1, f2) == MYB_E_DOMAIN);
  MT_CHECK(myb_tle_parse(l1, l2, 0, &t, &diag) == MYB_OK);
  myb_tle bad = t;
  bad.eccentricity = 1.0;
  MT_CHECK(myb_tle_format(&bad, f1, f2) == MYB_E_DOMAIN);
  bad = t;
  bad.eccentricity = 0.99999996; /* rounds to 1 */
  MT_CHECK(myb_tle_format(&bad, f1, f2) == MYB_E_DOMAIN);
  bad = t;
  bad.bstar = 1.0e12;
  MT_CHECK(myb_tle_format(&bad, f1, f2) == MYB_E_DOMAIN);
  bad = t;
  bad.ndot = 1.5;
  MT_CHECK(myb_tle_format(&bad, f1, f2) == MYB_E_DOMAIN);
  bad = t;
  bad.ndot = 0.999999999;
  MT_CHECK(myb_tle_format(&bad, f1, f2) == MYB_E_DOMAIN);
  bad = t;
  bad.inclination = 1000.0;
  MT_CHECK(myb_tle_format(&bad, f1, f2) == MYB_E_DOMAIN);
  bad = t;
  bad.mean_motion = (double)NAN;
  MT_CHECK(myb_tle_format(&bad, f1, f2) == MYB_E_DOMAIN);
  bad = t;
  bad.revolution = 100000;
  MT_CHECK(myb_tle_format(&bad, f1, f2) == MYB_E_DOMAIN);

  /* a mean-elements message of the CCSDS: the same numbers as an element set */
  static const char omm_text[] =
      "CCSDS_OMM_VERS = 2.0\nCREATION_DATE = 2020-06-01T10:00:00\nORIGINATOR = 18 SPCS\n"
      "OBJECT_NAME = ISS (ZARYA)\nOBJECT_ID = 1998-067A\nCENTER_NAME = EARTH\n"
      "REF_FRAME = TEME\nTIME_SYSTEM = UTC\nMEAN_ELEMENT_THEORY = SGP4\n"
      "EPOCH = 2020-06-01T03:22:38.887392\nMEAN_MOTION = 15.49438371\n"
      "ECCENTRICITY = 0.0002335\nINCLINATION = 51.6443\nRA_OF_ASC_NODE = 93.0000\n"
      "ARG_OF_PERICENTER = 45.4201\nMEAN_ANOMALY = 314.6988\nEPHEMERIS_TYPE = 0\n"
      "CLASSIFICATION_TYPE = U\nNORAD_CAT_ID = 25544\nELEMENT_SET_NO = 999\n"
      "REV_AT_EPOCH = 22954\nBSTAR = 0.000024949\nMEAN_MOTION_DOT = 0.00000907\n"
      "MEAN_MOTION_DDOT = 0\n";
  myb_omm omm;
  MT_CHECK(myb_omm_parse(omm_text, strlen(omm_text), &omm, &diag) == MYB_OK);
  MT_CHECK(myb_tle_from_omm(&omm, &t) == MYB_OK);
  MT_CHECK(t.catalog_number == 25544 && strcmp(t.designator, "98067A") == 0);
  MT_CHECK(t.mean_motion == 15.49438371 && t.bstar == 0.000024949 && t.ndot == 0.00000907);
  MT_CHECK(myb_tle_format(&t, f1, f2) == MYB_OK);
  MT_CHECK(memcmp(f1, "1 25544U 98067A   20153.14072786  .00000907  00000-0  24949-4 0  999",
                  68) == 0);
  MT_CHECK(memcmp(f2, "2 25544  51.6443  93.0000 0002335  45.4201 314.6988 15.4943837122954",
                  68) == 0);
  MT_CHECK(myb_tle_parse(f1, f2, 0, &back, &diag) == MYB_OK);
  /* the theory gets the same elements from the message as from the lines */
  myb_sgp4 from_message;
  myb_sgp4 from_lines;
  MT_CHECK(myb_sgp4_init(&t, MYB_SGP4_WGS72, MYB_SGP4_IMPROVED, &from_message) == MYB_OK);
  MT_CHECK(myb_sgp4_init(&back, MYB_SGP4_WGS72, MYB_SGP4_IMPROVED, &from_lines) == MYB_OK);
  double r1[3];
  double v1[3];
  double r2[3];
  double v2[3];
  MT_CHECK(invoke(&from_message, 720.0, r1, v1) == 0 && invoke(&from_lines, 720.0, r2, v2) == 0);
  double apart = 0.0;
  for (int i = 0; i < 3; ++i) {
    apart = fmax(apart, fabs(r1[i] - r2[i]));
  }
  /* the same motion from the epoch; the epoch of the lines keeps 1e-8 day */
  MT_CHECK(apart == 0.0);
  MT_CHECK(fabs(t.epoch_day - back.epoch_day) <= 0.5e-8 * (1.0 + 1.0e-6));
  MT_CHECK(fabs(t.epoch_day - (153.0 + 12158.887392 / 86400.0)) < 1.0e-12);
  /* and the element set as a message */
  myb_omm made;
  MT_CHECK(myb_omm_from_tle(&back, &made) == MYB_OK);
  MT_CHECK(strcmp(made.object_id, "1998-067A") == 0 && strcmp(made.ref_frame, "TEME") == 0 &&
           strcmp(made.mean_element_theory, "SGP4") == 0 && strcmp(made.time_system, "UTC") == 0);
  MT_CHECK(made.has_mean_motion && made.has_tle && made.norad_cat_id == 25544 &&
           made.element_set_no == 999 && made.rev_at_epoch == 22954);
  MT_CHECK(made.mean_motion == omm.mean_motion && made.eccentricity == omm.eccentricity &&
           made.inclination == omm.inclination && made.bstar == omm.bstar &&
           made.mean_motion_dot == omm.mean_motion_dot);
  MT_CHECK(made.epoch.year == 2020 && made.epoch.month == 6 && made.epoch.day == 1 &&
           made.epoch.hour == 3 && made.epoch.minute == 22);
  MT_CHECK(fabs(made.epoch.second - 38.887392) < 1.0e-3);
  char text2[2048];
  size_t length = 0;
  MT_CHECK(myb_omm_write(&made, text2, sizeof text2, &length) == MYB_OK);
  MT_CHECK(myb_omm_parse(text2, length, &omm, &diag) == MYB_OK);
  MT_CHECK(myb_tle_from_omm(&omm, &t) == MYB_OK);
  char g1[70];
  char g2[70];
  MT_CHECK(myb_tle_format(&t, g1, g2) == MYB_OK);
  MT_CHECK(strcmp(g1, f1) == 0 && strcmp(g2, f2) == 0);
  /* a message that is not an element set */
  omm.has_tle = 0;
  MT_CHECK(myb_tle_from_omm(&omm, &t) == MYB_E_DOMAIN);
  omm.has_tle = 1;
  omm.has_mean_motion = 0;
  MT_CHECK(myb_tle_from_omm(&omm, &t) == MYB_E_DOMAIN);
}

/* Writing type 10: the words of the segments that CSPICE spkw10_c wrote (t10.bsp: three and
 * two element sets; t10_many.bsp: 230, with a directory of the epochs) are rebuilt from their
 * constants and packets; the file written from them is read as the original; the packet of
 * a two-line set against getelm_c. */
static void test_spk10_writer(void) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const char* kernels[] = {MYB_DATA_DIR "/time/time-scales.tk", MYB_DATA_DIR "/time/tdb-tt-fb1990.tk",
                           MYB_DATA_DIR "/time/utc-tai.tk", MYB_DATA_DIR "/frames/iau1980.tk"};
  for (int i = 0; i < 4; ++i) {
    MT_CHECK(myb_kpool_load_file(pool, kernels[i], &diag) == MYB_OK);
  }
  myb_timesys* sys = NULL;
  myb_equinox* eq = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, &diag) == MYB_OK);
  MT_CHECK(myb_equinox_create(pool, "IAU1980", NULL, &eq, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  if (sys == NULL || eq == NULL) {
    myb_equinox_destroy(eq);
    myb_timesys_destroy(sys);
    return;
  }
  enum { MOST_WORDS = 4096, MOST_SETS = 256 };
  static double words[MOST_WORDS];
  static double rebuilt[MOST_WORDS];
  static double packets[14 * MOST_SETS];
  static unsigned char image[65536];
  static const char* const files[2] = {MYB_TEST_DATA_DIR "/spk/t10.bsp",
                                       MYB_TEST_DATA_DIR "/spk/t10_many.bsp"};
  int segments = 0;
  int sets = 0;
  int different = 0;
  int states = 0;
  double worst_nutation = 0.0;
  double first_packets[2][14];
  for (int file = 0; file < 2; ++file) {
    myb_spk* spk = NULL;
    MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
    MT_CHECK(myb_spk_load_file(spk, files[file], &diag) == MYB_OK);
    myb_spk_writer* writer = NULL;
    MT_CHECK(myb_spk_writer_create("copy", NULL, &writer) == MYB_OK);
    const size_t count = myb_spk_segment_count(spk);
    for (size_t index = 0; index < count; ++index) {
      myb_spk_segment_info info;
      myb_spk_segment(spk, index, &info);
      const size_t size = myb_spk_segment_size(spk, index);
      MT_CHECK(info.type == 10 && size <= MOST_WORDS);
      MT_CHECK(myb_spk_segment_data(spk, index, 0, size, words) == MYB_OK);
      /* the sets: the meta data say how many; each packet follows one word */
      const size_t n = (size_t)words[size - 17 + 6];
      MT_CHECK(n >= 1 && n <= MOST_SETS && myb_sgp4_spk10_size(n) == size);
      for (size_t k = 0; k < n; ++k) {
        memcpy(packets + 14 * k, words + 8 + 15 * k + 1, 14 * sizeof(double));
        /* the nutation SPKW10 stored against that of the library at the epoch */
        double nutation[2];
        double rates[2];
        myb_equinox_nutation(eq, myb_time_from_seconds(packets[14 * k + 9]), nutation, rates);
        worst_nutation = fmax(worst_nutation, fabs(nutation[1] - packets[14 * k + 10]));
        worst_nutation = fmax(worst_nutation, fabs(nutation[0] - packets[14 * k + 11]));
        worst_nutation = fmax(worst_nutation, 86400.0 * fabs(rates[1] - packets[14 * k + 12]));
        worst_nutation = fmax(worst_nutation, 86400.0 * fabs(rates[0] - packets[14 * k + 13]));
      }
      if (file == 0 && index < 2) {
        memcpy(first_packets[index], packets, sizeof first_packets[index]);
      }
      MT_CHECK(myb_sgp4_spk10_words(words, packets, n, rebuilt, MOST_WORDS) == MYB_OK);
      different += memcmp(words, rebuilt, size * sizeof(double)) != 0;
      const myb_spk_segment_header header = {
          .target = info.target, .center = info.center, .frame = info.frame, .name = info.name};
      MT_CHECK(myb_spk_writer_add_raw(writer, &header, 10, info.start, info.stop, rebuilt, size) ==
               MYB_OK);
      ++segments;
      sets += (int)n;
      /* too small a buffer, epochs out of order */
      MT_CHECK(myb_sgp4_spk10_words(words, packets, n, rebuilt, size - 1) == MYB_E_OVERFLOW);
      if (n > 1) {
        const double keep = packets[14 + 9];
        packets[14 + 9] = packets[9];
        MT_CHECK(myb_sgp4_spk10_words(words, packets, n, rebuilt, MOST_WORDS) == MYB_E_DOMAIN);
        packets[14 + 9] = keep;
      }
    }
    /* the file written from the words: the same states, bit for bit */
    const size_t bytes = myb_spk_writer_size(writer);
    MT_CHECK(bytes > 0 && bytes <= sizeof image);
    MT_CHECK(myb_spk_writer_write(writer, image, sizeof image) == MYB_OK);
    myb_spk* copy = NULL;
    MT_CHECK(myb_spk_create(NULL, &copy) == MYB_OK);
    MT_CHECK(myb_spk_load_memory(copy, image, bytes, &diag) == MYB_OK);
    MT_CHECK(myb_spk_segment_count(copy) == count);
    for (size_t index = 0; index < count; ++index) {
      myb_spk_segment_info info;
      myb_spk_segment_info same;
      myb_spk_segment(spk, index, &info);
      myb_spk_segment(copy, index, &same);
      MT_CHECK(info.target == same.target && info.center == same.center &&
               info.frame == same.frame && same.type == 10 && info.start == same.start &&
               info.stop == same.stop && strcmp(info.name, same.name) == 0);
      for (int k = 0; k <= 12; ++k) {
        const double et = info.start + (info.stop - info.start) * (double)k / 12.0;
        double a[6];
        double b[6];
        MT_CHECK(myb_sgp4_spk10_state(spk, index, sys, eq, myb_time_from_seconds(et), a) ==
                 MYB_OK);
        MT_CHECK(myb_sgp4_spk10_state(copy, index, sys, eq, myb_time_from_seconds(et), b) ==
                 MYB_OK);
        MT_CHECK(memcmp(a, b, sizeof a) == 0);
        ++states;
      }
    }
    myb_spk_destroy(copy);
    myb_spk_writer_destroy(writer);
    myb_spk_destroy(spk);
  }
  printf("SPK type 10 written: %d segments of %d element sets rebuilt from their packets, %d "
         "differ from the words of SPKW10; %d states of the copies as of the files; nutation of "
         "the packets within %.1e rad (rates: a day)\n",
         segments, sets, different, states, worst_nutation);
  MT_CHECK(segments == 3 && sets == 235 && different == 0 && states == 39);
  MT_CHECK(worst_nutation < 1e-9);
  MT_CHECK(myb_sgp4_spk10_size(0) == 0 && myb_sgp4_spk10_size(1) == 8 + 15 + 1 + 17 &&
           myb_sgp4_spk10_size(101) == 8 + 101 * 16 + 1 + 17);

  /* the packet of a two-line set against getelm_c: the first sets of the two series are
   * Vanguard 1 and Molniya 2-14 of the verification file */
  char* tle_text = read_all(MYB_TEST_DATA_DIR "/sgp4/SGP4-VER.TLE");
  MT_CHECK(tle_text != NULL);
  if (tle_text != NULL) {
    static char* lines[400];
    const int n = split_lines(tle_text, lines, 400);
    static const char* const numbers[2] = {"00005", "08195"};
    double worst = 0.0;
    double worst_epoch = 0.0;
    int found = 0;
    for (int which = 0; which < 2; ++which) {
      for (int i = 0; i + 1 < n; ++i) {
        if (lines[i][0] != '1' || strncmp(lines[i] + 2, numbers[which], 5) != 0) {
          continue;
        }
        myb_tle tle;
        MT_CHECK(myb_tle_parse(lines[i], lines[i + 1], MYB_TLE_NO_CHECKSUM, &tle, &diag) ==
                 MYB_OK);
        double packet[14];
        MT_CHECK(myb_sgp4_spk10_packet(&tle, sys, eq, packet) == MYB_OK);
        const double* expected = first_packets[which];
        for (int k = 0; k < 9; ++k) {
          const double scale = fabs(expected[k]) > 0.0 ? fabs(expected[k]) : 1.0;
          worst = fmax(worst, fabs(packet[k] - expected[k]) / scale);
        }
        worst_epoch = fmax(worst_epoch, fabs(packet[9] - expected[9]));
        for (int k = 10; k < 14; ++k) {
          MT_CHECK(fabs(packet[k] - expected[k]) < 1e-9);
        }
        /* without the nutation model: zeros */
        MT_CHECK(myb_sgp4_spk10_packet(&tle, sys, NULL, packet) == MYB_OK);
        MT_CHECK(packet[10] == 0.0 && packet[13] == 0.0 && packet[9] != 0.0);
        ++found;
        break;
      }
    }
    printf("packets of two-line sets against getelm_c: elements within %.1e, epochs within "
           "%.1e s (SPICE's TDB - TT has one term)\n",
           worst, worst_epoch);
    MT_CHECK(found == 2 && worst < 1e-14 && worst_epoch < 1e-4);
    free(tle_text);
  }
  double constants[8];
  myb_sgp4_constants model;
  myb_sgp4_constants_of(MYB_SGP4_WGS72, &model);
  myb_sgp4_spk10_constants(&model, constants);
  MT_CHECK(constants[0] == model.j2 && constants[3] == model.xke && constants[4] == 120.0 &&
           constants[5] == 78.0 && constants[6] == model.radius && constants[7] == 1.0);
  myb_equinox_destroy(eq);
  myb_timesys_destroy(sys);
}

/* The largest relative difference of the six elements of two sets (angles by their sines
 * and cosines, so that 359.99 and 0.01 degrees are close). */
static double elements_apart(const myb_tle* a, const myb_tle* b) {
  const double rad = 3.141592653589793 / 180.0;
  double worst = fabs(a->mean_motion / b->mean_motion - 1.0);
  worst = fmax(worst, fabs(a->eccentricity - b->eccentricity));
  worst = fmax(worst, fabs(a->inclination - b->inclination) * rad);
  /* node, and the sum that stays defined on a circle */
  worst = fmax(worst, fabs(sin((a->node - b->node) * rad)));
  const double la = a->node + a->argument_of_perigee + a->mean_anomaly;
  const double lb = b->node + b->argument_of_perigee + b->mean_anomaly;
  worst = fmax(worst, fabs(sin((la - lb) * rad)));
  return worst;
}

static void test_fit(void) {
  /* positions of the theory itself: the set that made them comes back */
  struct {
    const char* name;
    const char* l1;
    const char* l2;
    double span; /* minutes */
    int bstar;
  } cases[3] = {
      {"low orbit, 2 days, with B*",
       "1 25544U 98067A   20153.14072786  .00000907  00000-0  24949-4 0  9994",
       "2 25544  51.6443  93.0000 0002335  45.4201 314.6988 15.49438371229542", 2880.0, 1},
      {"eccentric, 3 days",
       "1 00005U 58002B   20287.20333880 -.00000016  00000-0 -22483-4 0  9998",
       "2 00005  34.2443 225.5254 1845686 162.2516 205.2356 10.84869164218149", 4320.0, 0},
      {"12-hour orbit (deep space), 4 days", NULL, NULL, 5760.0, 0}};
  /* the third from the verification file: the first set of the deep-space branch of the
   * theory (a period above 225 minutes) */
  char* text = read_all(MYB_TEST_DATA_DIR "/sgp4/SGP4-VER.TLE");
  MT_CHECK(text != NULL);
  if (text == NULL) {
    return;
  }
  static char* lines[400];
  const int line_count = split_lines(text, lines, 400);
  for (int i = 0; i + 1 < line_count && cases[2].l1 == NULL; ++i) {
    myb_tle candidate;
    myb_diag unused = {0};
    if (lines[i][0] == '1' && lines[i + 1][0] == '2' &&
        myb_tle_parse(lines[i], lines[i + 1], 0, &candidate, &unused) == MYB_OK &&
        candidate.mean_motion < 3.0 && candidate.eccentricity < 0.75) {
      cases[2].l1 = lines[i];
      cases[2].l2 = lines[i + 1];
    }
  }
  MT_CHECK(cases[2].l1 != NULL);
  for (int k = 0; k < 3 && cases[2].l1 != NULL; ++k) {
    myb_tle truth;
    myb_diag diag = {0};
    MT_CHECK(myb_tle_parse(cases[k].l1, cases[k].l2, 0, &truth, &diag) == MYB_OK);
    myb_sgp4 theory;
    MT_CHECK(myb_sgp4_init(&truth, MYB_SGP4_WGS72, MYB_SGP4_IMPROVED, &theory) == MYB_OK);
    enum { POINTS = 97 };
    static double minutes[POINTS];
    static double positions[3 * POINTS];
    double v[3];
    double v0[3];
    double r0[3];
    for (int i = 0; i < POINTS; ++i) {
      minutes[i] = cases[k].span * i / (POINTS - 1);
      MT_CHECK(myb_sgp4_propagate(&theory, minutes[i], positions + 3 * i, v, NULL) == MYB_OK);
      if (i == 0) {
        memcpy(r0, positions, sizeof r0);
        memcpy(v0, v, sizeof v0);
      }
    }
    /* the guess: the osculating elements at the epoch, the identification of the set, a
     * drag term off by a half */
    myb_tle fitted = truth;
    fitted.bstar = cases[k].bstar ? 1.5 * truth.bstar : truth.bstar;
    MT_CHECK(myb_tle_from_state(r0, v0, MYB_SGP4_WGS72, &fitted) == MYB_OK);
    const double guess = elements_apart(&fitted, &truth);
    const myb_tle_fit_options options = {.bstar = cases[k].bstar, .gravity = MYB_SGP4_WGS72};
    double rms = 0.0;
    int iterations = 0;
    MT_CHECK(myb_tle_fit(minutes, positions, POINTS, &options, &fitted, &rms, &iterations) ==
             MYB_OK);
    printf("fit of an element set, %-34s: from a guess off by %.1e to %.1e in %d steps; "
           "residuals %.1e m",
           cases[k].name, guess, elements_apart(&fitted, &truth), iterations, rms);
    if (cases[k].bstar) {
      printf("; B* %.5e (%.5e)", fitted.bstar, truth.bstar);
      MT_CHECK(fabs(fitted.bstar / truth.bstar - 1.0) < 1e-5);
    }
    printf("\n");
    MT_CHECK(guess > 1e-5);
    MT_CHECK(elements_apart(&fitted, &truth) < 1e-9);
    MT_CHECK(rms < 1e-3);
    MT_CHECK(strcmp(fitted.designator, truth.designator) == 0 &&
             fitted.epoch_day == truth.epoch_day);
  }
  /* refused */
  myb_tle tle;
  myb_diag diag = {0};
  MT_CHECK(myb_tle_parse(cases[0].l1, cases[0].l2, 0, &tle, &diag) == MYB_OK);
  free(text);
  const double two[2] = {0.0, 1.0};
  const double points[6] = {7.0e6, 0.0, 0.0, 0.0, 7.0e6, 0.0};
  MT_CHECK(myb_tle_fit(two, points, 2, NULL, &tle, NULL, NULL) == MYB_E_DOMAIN);
  const double fast[3] = {0.0, 12000.0, 0.0};
  const double here[3] = {7.0e6, 0.0, 0.0};
  MT_CHECK(myb_tle_from_state(here, fast, MYB_SGP4_WGS72, &tle) == MYB_E_DOMAIN);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  run_verification();
  test_parser();
  test_format();
  test_modes();
  test_spk10();
  test_spk10_writer();
  test_fit();
  return mt_finish("test_sgp4");
}
