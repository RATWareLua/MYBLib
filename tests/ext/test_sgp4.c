/* SPDX-License-Identifier: BSD-3-Clause */
/* SGP4 (ext/sgp4) against the verification output of Vallado et al. (2006): the 33 element sets
 * of SGP4-VER.TLE propagated over their time spans as the verification program does (WGS-72,
 * improved mode), compared line by line with tcppver.out (positions printed to 1e-8 km,
 * velocities to 1e-9 km/s), and the sequence of error codes [1, 1, 6, 6, 4, 3, 6]. Then the TLE
 * parser (Alpha-5, checksums, epochs), the modes and the error statuses. */
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

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  run_verification();
  test_parser();
  test_modes();
  test_spk10();
  return mt_finish("test_sgp4");
}
