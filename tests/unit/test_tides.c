/* SPDX-License-Identifier: BSD-3-Clause */
/* Tides of a gravity field (IERS Conventions 2010, ch. 6):
 *  - step 1 against the addition theorem: with equal real Love numbers the tidal potential of
 *    degree n must be k_n GM_j R^(2n+1) / (r^(n+1) r_j^(n+1)) P_n(cos psi), whatever the
 *    normalisation and the split into orders;
 *  - step 2 with the example of the IERS text (K1: dC21 = 470.9e-12 sin chi - 30.2e-12 cos chi),
 *    M2 and the 18.6-year term, chi = GMST + pi from ERFA (eraGmst06) and Omega from eraFaom03;
 *  - the permanent tide, the pole tides (formulas of sections 6.4, 6.5 and the secular pole);
 *  - ocean tides (eq. 6.15) from an excerpt of FES2004, the Doodson arguments built from ERFA's
 *    fundamental arguments (tau = chi - s, s = F + Omega);
 *  - errors. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "erfa.h"
#include "mtest.h"
#include "myb/grav.h"
#include "myb/io.h"
#include "myb/time.h"

static const char arguments[] =
    "\\begindata\n"
    "T_ARGUMENT_COUNT = 6\n"
    "T_ARGUMENT_SCALE = ( 4.848136811095359935899141D-6 4.848136811095359935899141D-6\n"
    "  4.848136811095359935899141D-6 4.848136811095359935899141D-6\n"
    "  4.848136811095359935899141D-6 4.848136811095359935899141D-6 )\n"
    "T_ARGUMENT_TURN = ( 1296000 1296000 1296000 1296000 1296000 1296000 )\n"
    "T_ARGUMENT_POLYNOMIAL = (\n"
    "  1657658.22615 47466002772.19299 1.39656 -9.3D-5 0\n"
    "  485868.249036 1717915923.2178 31.8792 0.051635 -0.00024470\n"
    "  1287104.793048 129596581.0481 -0.5532 0.000136 -0.00001149\n"
    "  335779.526232 1739527262.8478 -12.7512 -0.001037 0.00000417\n"
    "  1072260.703692 1602961601.2090 -6.3706 0.006593 -0.00003169\n"
    "  450160.398036 -6962890.5431 7.4722 0.007702 -0.00005939 )\n";

/* Step 2 reduced to three rows of Tables 6.5: 55,565 (C20), K1 (C21, S21), M2 (C22, S22). */
static const char reduced[] =
    "T_C20_SCALE = 1D-12\nT_C20_POLYNOMIAL = ( 0 )\n"
    "T_C20_TERMS = ( 0 6.7 16.6 0 0 0 0 0 -1 )\n"
    "T_C21_SCALE = 1D-12\nT_C21_POLYNOMIAL = ( 0 )\n"
    "T_C21_TERMS = ( 0 470.9 -30.2 1 0 0 0 0 0 )\n"
    "T_S21_SCALE = 1D-12\nT_S21_POLYNOMIAL = ( 0 )\n"
    "T_S21_TERMS = ( 0 30.2 470.9 1 0 0 0 0 0 )\n"
    "T_C22_SCALE = 1D-12\nT_C22_POLYNOMIAL = ( 0 )\n"
    "T_C22_TERMS = ( 0 0 -1.2 2 0 0 -2 0 -2 )\n"
    "T_S22_SCALE = 1D-12\nT_S22_POLYNOMIAL = ( 0 )\n"
    "T_S22_TERMS = ( 0 1.2 0 2 0 0 -2 0 -2 )\n";

static myb_kpool* pool_of(const char* love) {
  char text[4096];
  snprintf(text, sizeof text,
           "%s%s%s\nT_PERMANENT = ( 4.4228D-8 -0.31460 )\n"
           "T_POLE_SOLID = ( -1.333D-9 0.0115 -1.333D-9 -0.0115 )\n"
           "T_POLE_OCEAN = ( -2.1778D-10 -0.01724 -1.7232D-10 -0.03365 )\n"
           "T_SECULAR_POLE = ( 55.0 1.677 320.5 3.460 )\n",
           arguments, reduced, love);
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_text(pool, text, strlen(text), &diag) == MYB_OK);
  return pool;
}

static myb_gravity_field* empty_field(int degree, const char* tide) {
  myb_gravity_field* f = NULL;
  /* through the ICGEM reader to get a tide system */
  char text[512];
  snprintf(text, sizeof text,
           "earth_gravity_constant 3.986004415e14\nradius 6378136.3\nmax_degree %d\n"
           "norm fully_normalized\ntide_system %s\nend_of_head\ngfc 0 0 1.0 0 0 0\n",
           degree, tide);
  myb_diag diag = {0};
  MT_CHECK(myb_gravity_field_load_icgem(text, strlen(text), 0, NULL, &f, &diag) == MYB_OK);
  return f;
}

static double legendre(int n, double x) {
  return n == 2 ? 0.5 * (3.0 * x * x - 1.0) : 0.5 * (5.0 * x * x * x - 3.0 * x);
}

static void test_addition_theorem(void) {
  myb_kpool* pool = pool_of("T_LOVE = ( 2 0 0.3 0 0  2 1 0.3 0 0  2 2 0.3 0 0\n"
                            "  3 0 0.09 0 0  3 1 0.09 0 0  3 2 0.09 0 0  3 3 0.09 0 0 )\n");
  myb_tides* tides = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_tides_create(pool, "T", NULL, &tides, &diag) == MYB_OK);
  if (tides == NULL) {
    printf("%s\n", diag.message);
    myb_kpool_destroy(pool);
    return;
  }
  myb_gravity_field* f = empty_field(4, "tide_free");
  myb_gravity_work* w = NULL;
  MT_CHECK(myb_gravity_work_create(4, NULL, &w) == MYB_OK);

  const double R = 6378136.3;
  const myb_tide_raiser raisers[2] = {{4.902800e12, {2.1e8, -3.0e8, 1.1e8}},
                                      {1.32712440018e20, {-1.1e11, 0.7e11, 0.4e11}}};
  const myb_tides_input in = {.raisers = raisers, .count = 2, .flags = MYB_TIDES_SOLID};
  MT_CHECK(myb_tides_apply(tides, f, &in) == MYB_OK);
  const double sat[3][3] = {{7.0e6, 1.0e5, -2.0e5}, {1.0e3, -2.0e3, 1.2e7}, {-5e6, 6e6, 3e6}};
  double worst = 0.0;
  for (int p = 0; p < 3; ++p) {
    const double r = sqrt(sat[p][0] * sat[p][0] + sat[p][1] * sat[p][1] + sat[p][2] * sat[p][2]);
    for (int n = 2; n <= 3; ++n) {
      double v = 0.0;
      MT_CHECK(myb_gravity_field_eval(f, w, sat[p], n, n, &v, NULL) == MYB_OK);
      double expected = 0.0;
      double scale = 0.0;
      for (int j = 0; j < 2; ++j) {
        const double* b = raisers[j].r;
        const double rj = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
        const double cpsi = (b[0] * sat[p][0] + b[1] * sat[p][1] + b[2] * sat[p][2]) / (r * rj);
        const double k = n == 2 ? 0.3 : 0.09;
        const double term = k * raisers[j].gm * pow(R, 2 * n + 1) / pow(r * rj, n + 1);
        expected += term * legendre(n, cpsi);
        scale += fabs(term);
      }
      worst = fmax(worst, fabs(v - expected) / scale);
    }
  }
  printf("step 1 vs the addition theorem (degrees 2, 3; Moon and Sun): %.1e of the scale\n",
         worst);
  MT_CHECK(worst < 1e-14);
  /* the field's own GM does not enter the result: dC scales with GM_j / GM */
  myb_gravity_work_destroy(w);
  myb_gravity_field_destroy(f);
  myb_tides_destroy(tides);
  myb_kpool_destroy(pool);
}

static myb_time utc_like(int y, int mo, int d, int h, int mi, double s) {
  const myb_calendar c = {y, mo, d, h, mi, s};
  myb_time t;
  MT_CHECK(myb_calendar_to_time(&c, &t) == MYB_OK);
  return t;
}

static void test_step2_and_poles(void) {
  myb_kpool* pool = pool_of("T_LOVE = ( 2 0 0.30190 0 -0.00089 )\n");
  myb_tides* tides = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_tides_create(pool, "T", NULL, &tides, &diag) == MYB_OK);
  if (tides == NULL) {
    printf("%s\n", diag.message);
    myb_kpool_destroy(pool);
    return;
  }
  double worst = 0.0;
  for (int k = 0; k < 40; ++k) {
    const myb_time t = myb_time_add(utc_like(2025, 6, 1, 0, 0, 0.0), 3600.0 * 37.3 * k);
    myb_gravity_field* f = empty_field(2, "tide_free");
    const myb_tides_input in = {.t = t, .flags = MYB_TIDES_FREQUENCY};
    MT_CHECK(myb_tides_apply(tides, f, &in) == MYB_OK);
    double jd1 = 0.0;
    double jd2 = 0.0;
    myb_time_to_jd(t, &jd1, &jd2);
    /* t is taken as UT1; TT 69 s later (irrelevant at this level) */
    const double chi = eraGmst06(jd1, jd2, jd1, jd2 + 69.184 / 86400.0) + 3.141592653589793;
    const double tc = ((jd1 - 2451545.0) + jd2) / 36525.0;
    const double om = eraFaom03(tc);
    const double f_ = eraFaf03(tc);
    double c = 0.0;
    double s = 0.0;
    MT_CHECK(myb_gravity_field_coefficient(f, 2, 1, &c, &s) == MYB_OK);
    worst = fmax(worst, fabs(c - (470.9e-12 * sin(chi) - 30.2e-12 * cos(chi))));
    worst = fmax(worst, fabs(s - (470.9e-12 * cos(chi) + 30.2e-12 * sin(chi))));
    const double m2 = 2.0 * chi - 2.0 * f_ - 2.0 * om; /* 2 tau */
    MT_CHECK(myb_gravity_field_coefficient(f, 2, 2, &c, &s) == MYB_OK);
    worst = fmax(worst, fabs(c - (-1.2e-12 * cos(m2))));
    worst = fmax(worst, fabs(s - (1.2e-12 * sin(m2))));
    MT_CHECK(myb_gravity_field_coefficient(f, 2, 0, &c, &s) == MYB_OK);
    /* 55,565: theta = -Omega (N' = -Omega), ip 16.6, op -6.7 */
    worst = fmax(worst, fabs(c - (16.6e-12 * cos(-om) + 6.7e-12 * sin(-om))));
    myb_gravity_field_destroy(f);
  }
  printf("step 2 (K1, M2, 18.6 yr) vs the IERS example with ERFA's GMST: %.1e\n", worst);
  MT_CHECK(worst < 2e-15);
  /* permanent tide: removed for a zero-tide field, kept for a tide-free one */
  const myb_tide_raiser none[1] = {{0.0, {1.0, 0.0, 0.0}}};
  myb_gravity_field* zt = empty_field(2, "zero_tide");
  myb_gravity_field* tf = empty_field(2, "tide_free");
  myb_gravity_field* mt = empty_field(2, "mean_tide");
  myb_tides_input in = {.raisers = none, .count = 1, .flags = MYB_TIDES_SOLID};
  MT_CHECK(myb_tides_apply(tides, zt, &in) == MYB_OK);
  MT_CHECK(myb_tides_apply(tides, tf, &in) == MYB_OK);
  MT_CHECK(myb_tides_apply(tides, mt, &in) == MYB_E_DOMAIN);
  in.permanent = MYB_PERMANENT_TIDE_KEEP;
  MT_CHECK(myb_tides_apply(tides, mt, &in) == MYB_OK);
  double c0 = 0.0;
  double c1 = 0.0;
  double s = 0.0;
  MT_CHECK(myb_gravity_field_coefficient(zt, 2, 0, &c0, &s) == MYB_OK);
  MT_CHECK(myb_gravity_field_coefficient(tf, 2, 0, &c1, &s) == MYB_OK);
  printf("permanent tide removed from a zero-tide field: %.5e\n", c0 - c1);
  MT_CHECK(fabs((c0 - c1) - 4.4228e-8 * 0.31460 * 0.30190) < 1e-22);
  myb_gravity_field_destroy(zt);
  myb_gravity_field_destroy(tf);
  myb_gravity_field_destroy(mt);
  /* pole tides at 2025.5 with xp = 0.15", yp = 0.40" */
  const double as = 4.84813681109536e-6;
  const myb_time t = utc_like(2025, 7, 2, 12, 0, 0.0);
  const double years = myb_time_seconds(t) / 31557600.0;
  const double m1 = 0.15 - (55.0 + 1.677 * years) * 1e-3;
  const double m2 = -(0.40 - (320.5 + 3.460 * years) * 1e-3);
  myb_gravity_field* f = empty_field(2, "tide_free");
  in = (myb_tides_input){.t = t, .xp = 0.15 * as, .yp = 0.40 * as,
                         .flags = MYB_TIDES_POLE_SOLID | MYB_TIDES_POLE_OCEAN};
  MT_CHECK(myb_tides_apply(tides, f, &in) == MYB_OK);
  double c = 0.0;
  MT_CHECK(myb_gravity_field_coefficient(f, 2, 1, &c, &s) == MYB_OK);
  const double ec = -1.333e-9 * (m1 + 0.0115 * m2) - 2.1778e-10 * (m1 - 0.01724 * m2);
  const double es = -1.333e-9 * (m2 - 0.0115 * m1) - 1.7232e-10 * (m2 - 0.03365 * m1);
  printf("pole tides: dC21 %.4e, dS21 %.4e\n", c, s);
  MT_CHECK(fabs(c - ec) < 1e-22 && fabs(s - es) < 1e-22);
  myb_gravity_field_destroy(f);
  /* errors */
  f = empty_field(2, "tide_free");
  const myb_tide_raiser origin[1] = {{1e12, {0.0, 0.0, 0.0}}};
  in = (myb_tides_input){.raisers = origin, .count = 1, .flags = MYB_TIDES_SOLID};
  MT_CHECK(myb_tides_apply(tides, f, &in) == MYB_E_DOMAIN);
  in = (myb_tides_input){.xp = (double)NAN, .flags = MYB_TIDES_POLE_SOLID};
  MT_CHECK(myb_tides_apply(tides, f, &in) == MYB_E_DOMAIN);
  myb_gravity_field_destroy(f);
  myb_tides_destroy(tides);
  myb_kpool_destroy(pool);
  /* incomplete kernels */
  myb_kpool* broken = pool_of("T_LOVE = ( 4 0 0.3 0 0 )\n");
  MT_CHECK(myb_tides_create(broken, "T", NULL, &tides, &diag) == MYB_E_FORMAT);
  myb_kpool_destroy(broken);
  broken = pool_of("");
  MT_CHECK(myb_tides_create(broken, "T", NULL, &tides, &diag) == MYB_E_NOT_FOUND);
  myb_kpool_destroy(broken);
  myb_tides_destroy(NULL);
}

static void test_ocean(void) {
  myb_kpool* pool = pool_of("T_LOVE = ( 2 0 0.30190 0 -0.00089 )\n");
  myb_tides* tides = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_tides_create(pool, "T", NULL, &tides, &diag) == MYB_OK);
  if (tides == NULL) {
    myb_kpool_destroy(pool);
    return;
  }
  MT_CHECK(myb_tides_load_ocean_file(tides, MYB_TEST_DATA_DIR "/grav/fes2004_excerpt.dat", 4,
                                     &diag) == MYB_OK);
  MT_CHECK(myb_tides_ocean_waves(tides) == 4);
  /* expected from eq. 6.15 with the Doodson arguments from ERFA */
  static const struct {
    double tau, s, h, p, np, ps; /* multipliers */
    int n, m;
    double cp, sp, cm, sm;
  } rows[] = {
      {0, 0, 0, 0, 1, 0, 2, 0, -6.58128, 0.0, 0.0, 0.0},             /* Om1 */
      {0, 2, 0, 0, 0, 0, 3, 1, 0.30789, 0.41601, 0.52036, -0.33269},  /* Mf */
      {1, 1, 0, 0, 0, 0, 2, 1, 20.22034, 18.62745, -9.55722, -10.28875}, /* K1 */
      {2, 0, 0, 0, 0, 0, 4, 4, -26.73398, 17.10630, -12.59975, -11.57187}, /* M2 */
      {2, 0, 0, 0, 0, 0, 2, 2, -39.36214, 46.75729, 9.57270, 5.24459},     /* M2 */
  };
  const myb_time t = utc_like(2025, 6, 3, 7, 30, 0.0);
  double jd1 = 0.0;
  double jd2 = 0.0;
  myb_time_to_jd(t, &jd1, &jd2);
  const double tc = ((jd1 - 2451545.0) + jd2) / 36525.0;
  const double chi = eraGmst06(jd1, jd2, jd1, jd2 + 69.184 / 86400.0) + 3.141592653589793;
  const double F = eraFaf03(tc);
  const double D = eraFad03(tc);
  const double l = eraFal03(tc);
  const double lp = eraFalp03(tc);
  const double om = eraFaom03(tc);
  const double s_ = F + om;
  const double beta[6] = {chi - s_, s_, s_ - D, s_ - l, -om, s_ - D - lp};
  myb_gravity_field* f = empty_field(4, "tide_free");
  const myb_tides_input in = {.t = t, .flags = MYB_TIDES_OCEAN};
  MT_CHECK(myb_tides_apply(tides, f, &in) == MYB_OK);
  /* (n, m) = (2, 0) gets Om1, Mf, K1, M2: compare the sum of all four from the file */
  double worst = 0.0;
  for (size_t k = 0; k < sizeof rows / sizeof rows[0]; ++k) {
    if (rows[k].n == 2 && rows[k].m == 0) {
      continue;
    }
    const double mult[6] = {rows[k].tau, rows[k].s, rows[k].h, rows[k].p, rows[k].np, rows[k].ps};
    double theta = 0.0;
    for (int i = 0; i < 6; ++i) {
      theta += mult[i] * beta[i];
    }
    const double ec = 1e-11 * ((rows[k].cp + rows[k].cm) * cos(theta) +
                               (rows[k].sp + rows[k].sm) * sin(theta));
    const double es = 1e-11 * ((rows[k].sp - rows[k].sm) * cos(theta) -
                               (rows[k].cp - rows[k].cm) * sin(theta));
    double c = 0.0;
    double s = 0.0;
    MT_CHECK(myb_gravity_field_coefficient(f, rows[k].n, rows[k].m, &c, &s) == MYB_OK);
    if (rows[k].n == 3 && rows[k].m == 1) {
      /* (3, 1) also gets K1 and M2: add them */
      static const double other[2][10] = {
          {1, 1, 0, 0, 0, 0, 7.35715, -5.00724, 11.03003, 11.53259},
          {2, 0, 0, 0, 0, 0, -2.97668, 12.33708, 18.47672, 24.89254}};
      double ec2 = ec;
      double es2 = es;
      for (int o = 0; o < 2; ++o) {
        double th = 0.0;
        for (int i = 0; i < 6; ++i) {
          th += other[o][i] * beta[i];
        }
        ec2 += 1e-11 * ((other[o][6] + other[o][8]) * cos(th) + (other[o][7] + other[o][9]) * sin(th));
        es2 += 1e-11 * ((other[o][7] - other[o][9]) * cos(th) - (other[o][6] - other[o][8]) * sin(th));
      }
      worst = fmax(worst, fmax(fabs(c - ec2), fabs(s - es2)));
      continue;
    }
    if (rows[k].n == 2 && rows[k].m == 1) {
      /* (2, 1) also gets Mf and M2 */
      static const double other[2][10] = {
          {0, 2, 0, 0, 0, 0, -0.07331, -0.41232, 0.18849, 0.19707},
          {2, 0, 0, 0, 0, 0, -12.07164, -4.38919, -3.09008, 1.50139}};
      double ec2 = ec;
      double es2 = es;
      for (int o = 0; o < 2; ++o) {
        double th = 0.0;
        for (int i = 0; i < 6; ++i) {
          th += other[o][i] * beta[i];
        }
        ec2 += 1e-11 * ((other[o][6] + other[o][8]) * cos(th) + (other[o][7] + other[o][9]) * sin(th));
        es2 += 1e-11 * ((other[o][7] - other[o][9]) * cos(th) - (other[o][6] - other[o][8]) * sin(th));
      }
      worst = fmax(worst, fmax(fabs(c - ec2), fabs(s - es2)));
      continue;
    }
    if (rows[k].n == 2 && rows[k].m == 2) {
      /* (2, 2) also gets Mf and K1 */
      static const double other[2][10] = {
          {0, 2, 0, 0, 0, 0, -0.68073, -0.10677, -0.04377, 0.09242},
          {1, 1, 0, 0, 0, 0, -1.73444, -18.68962, 0.42290, 3.10616}};
      double ec2 = ec;
      double es2 = es;
      for (int o = 0; o < 2; ++o) {
        double th = 0.0;
        for (int i = 0; i < 6; ++i) {
          th += other[o][i] * beta[i];
        }
        ec2 += 1e-11 * ((other[o][6] + other[o][8]) * cos(th) + (other[o][7] + other[o][9]) * sin(th));
        es2 += 1e-11 * ((other[o][7] - other[o][9]) * cos(th) - (other[o][6] - other[o][8]) * sin(th));
      }
      worst = fmax(worst, fmax(fabs(c - ec2), fabs(s - es2)));
      continue;
    }
    if (rows[k].n == 4 && rows[k].m == 4) {
      /* (4, 4) also gets Mf and K1 */
      static const double other[2][10] = {
          {0, 2, 0, 0, 0, 0, 0.19325, -0.57540, 0.11220, 0.58875},
          {1, 1, 0, 0, 0, 0, -7.87675, -6.82474, 15.02278, 0.65801}};
      double ec2 = ec;
      double es2 = es;
      for (int o = 0; o < 2; ++o) {
        double th = 0.0;
        for (int i = 0; i < 6; ++i) {
          th += other[o][i] * beta[i];
        }
        ec2 += 1e-11 * ((other[o][6] + other[o][8]) * cos(th) + (other[o][7] + other[o][9]) * sin(th));
        es2 += 1e-11 * ((other[o][7] - other[o][9]) * cos(th) - (other[o][6] - other[o][8]) * sin(th));
      }
      worst = fmax(worst, fmax(fabs(c - ec2), fabs(s - es2)));
    }
  }
  /* S_n0 unchanged, degree 1 skipped */
  double c = 0.0;
  double s = 0.0;
  MT_CHECK(myb_gravity_field_coefficient(f, 2, 0, &c, &s) == MYB_OK);
  MT_CHECK(s == 0.0);
  MT_CHECK(myb_gravity_field_coefficient(f, 1, 0, &c, &s) == MYB_OK);
  MT_CHECK(c == 0.0);
  printf("ocean tides (FES2004 excerpt) vs eq. 6.15 with ERFA's arguments: %.1e\n", worst);
  /* the kernel's chi uses GMST of 1982 (as the IERS routines), ERFA's GMST of 2006: they differ
   * by 0.056" in 2025, i.e. 5e-7 rad in 2 chi, times amplitudes up to 6e-10 */
  MT_CHECK(worst < 6e-16);
  myb_gravity_field_destroy(f);
  /* truncation and errors */
  MT_CHECK(myb_tides_load_ocean_file(tides, MYB_TEST_DATA_DIR "/grav/fes2004_excerpt.dat", 2,
                                     &diag) == MYB_OK);
  MT_CHECK(myb_tides_load_ocean(tides, "255.555 M2 2 2 1.0 2.0\n", 23, 4, &diag) ==
           MYB_E_FORMAT);
  MT_CHECK(diag.line == 1);
  MT_CHECK(myb_tides_load_ocean(tides, "255.555 M2 2 3 1.0 2.0 3.0 4.0\n", 31, 4, &diag) ==
           MYB_E_FORMAT);
  MT_CHECK(myb_tides_load_ocean(tides, "header only\n", 12, 4, &diag) == MYB_E_FORMAT);
  MT_CHECK(myb_tides_load_ocean(tides, "255.555 M2 2 2 1 2 3 4\n", 23, 1, &diag) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_tides_ocean_waves(tides) == 4); /* a failed load keeps the previous waves */
  /* the full FES2004 file, if downloaded */
  myb_tides* full = NULL;
  MT_CHECK(myb_tides_create(pool, "T", NULL, &full, &diag) == MYB_OK);
  const myb_status st = myb_tides_load_ocean_file(
      full, MYB_DATA_DIR "/download/fes2004_Cnm-Snm.dat", 100, &diag);
  if (st == MYB_OK) {
    printf("FES2004: %zu waves\n", myb_tides_ocean_waves(full));
    MT_CHECK(myb_tides_ocean_waves(full) == 18);
  }
  myb_tides_destroy(full);
  myb_tides_destroy(tides);
  myb_kpool_destroy(pool);
}

static void test_iers_kernel(void) {
  myb_kpool* pool = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/grav/iers2010-tides.tk", &diag) == MYB_OK);
  myb_tides* tides = NULL;
  MT_CHECK(myb_tides_create(pool, "IERS10_TIDE", NULL, &tides, &diag) == MYB_OK);
  if (tides == NULL) {
    printf("%s\n", diag.message);
    myb_kpool_destroy(pool);
    return;
  }
  /* step 2 over a month: dC21 is dominated by K1 (470.9e-12) and P1, psi1 */
  double peak = 0.0;
  for (int k = 0; k < 24 * 30; ++k) {
    const myb_time t = myb_time_add(utc_like(2025, 6, 1, 0, 0, 0.0), 3600.0 * k);
    myb_gravity_field* f = empty_field(2, "tide_free");
    const myb_tides_input in = {.t = t, .flags = MYB_TIDES_FREQUENCY};
    MT_CHECK(myb_tides_apply(tides, f, &in) == MYB_OK);
    double c = 0.0;
    double s = 0.0;
    MT_CHECK(myb_gravity_field_coefficient(f, 2, 1, &c, &s) == MYB_OK);
    peak = fmax(peak, hypot(c, s));
    myb_gravity_field_destroy(f);
  }
  printf("IERS 2010 step 2 over a month: |dC21 - i dS21| up to %.3e\n", peak);
  MT_CHECK(peak > 4.0e-10 && peak < 6.5e-10);
  myb_tides_destroy(tides);
  myb_kpool_destroy(pool);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_addition_theorem();
  test_step2_and_poles();
  test_ocean();
  test_iers_kernel();
  return mt_finish("test_tides");
}
