/* SPDX-License-Identifier: BSD-3-Clause */
/* Gravity fields: the ICGEM reader (format 1.0 excerpt of EGM2008, a synthetic format 2.0 file
 * with time-variable terms), closed forms (point mass, J2), the split by degrees, the poles,
 * errors; the polyhedron against the closed form of a rectangular prism (D. Nagy, G. Papp,
 * J. Benedek, J. Geodesy 74, 552, 2000) evaluated with 40 digits (mpmath), the solid-angle
 * Laplacian, mesh validation. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/grav.h"
#include "myb/time.h"

static double rel(const double a[3], const double b[3]) {
  const double d = hypot(hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2]);
  return d / hypot(hypot(b[0], b[1]), b[2]);
}

/* ---- ICGEM ----------------------------------------------------------------------------------- */

static void test_icgem_egm2008(void) {
  myb_gravity_field* f = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_gravity_field_load_icgem_file(MYB_TEST_DATA_DIR "/grav/EGM2008_n20.gfc", 0, NULL,
                                             &f, &diag) == MYB_OK);
  if (f == NULL) {
    printf("%s\n", diag.message);
    return;
  }
  MT_CHECK(myb_gravity_field_degree(f) == 20);
  MT_CHECK_SAME(myb_gravity_field_gm(f), 0.3986004415E+15);
  MT_CHECK_SAME(myb_gravity_field_radius(f), 0.63781363E+07);
  MT_CHECK(strcmp(myb_gravity_field_tide_system(f), "tide_free") == 0);
  double c = 0.0;
  double s = 0.0;
  MT_CHECK(myb_gravity_field_coefficient(f, 2, 0, &c, &s) == MYB_OK);
  MT_CHECK_SAME(c, -0.484165143790815e-03);
  MT_CHECK(myb_gravity_field_coefficient(f, 2, 1, &c, &s) == MYB_OK);
  MT_CHECK_SAME(c, -0.206615509074176e-09);
  MT_CHECK_SAME(s, 0.138441389137979e-08);
  MT_CHECK(myb_gravity_field_coefficient(f, 0, 0, &c, &s) == MYB_OK);
  MT_CHECK_SAME(c, 1.0); /* "1.0d0": Fortran exponent letter */
  MT_CHECK(myb_gravity_field_coefficient(f, 1, 0, &c, &s) == MYB_OK);
  MT_CHECK(c == 0.0 && s == 0.0); /* absent lines are zero */
  MT_CHECK(myb_gravity_field_coefficient(f, 21, 0, &c, &s) == MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_field_coefficient(f, 3, 4, &c, &s) == MYB_E_DOMAIN);
  /* surface gravity: about 9.798 at the equator, 9.832 at the poles (with the rotation of the
   * Earth: 9.780 / 9.832) */
  myb_gravity_work* w = NULL;
  MT_CHECK(myb_gravity_work_create(20, NULL, &w) == MYB_OK);
  double a[3];
  double v = 0.0;
  MT_CHECK(myb_gravity_field_eval(f, w, (const double[3]){6378137.0, 0.0, 0.0}, 0, -1, &v, a) ==
           MYB_OK);
  const double g_eq = -a[0];
  MT_CHECK(myb_gravity_field_eval(f, w, (const double[3]){0.0, 0.0, 6356752.3}, 0, -1, &v, a) ==
           MYB_OK);
  const double g_pole = -a[2];
  printf("EGM2008 n<=20: |g| equator %.5f (with rotation %.5f), pole %.5f m/s^2\n", g_eq,
         g_eq - 7.292115e-5 * 7.292115e-5 * 6378137.0, g_pole);
  MT_CHECK(fabs(g_eq - 7.292115e-5 * 7.292115e-5 * 6378137.0 - 9.7803) < 2e-3);
  MT_CHECK(fabs(g_pole - 9.8322) < 2e-3);
  myb_gravity_work_destroy(w);
  myb_gravity_field_destroy(f);
  /* truncation */
  MT_CHECK(myb_gravity_field_load_icgem_file(MYB_TEST_DATA_DIR "/grav/EGM2008_n20.gfc", 8, NULL,
                                             &f, &diag) == MYB_OK);
  MT_CHECK(myb_gravity_field_degree(f) == 8);
  myb_gravity_field_destroy(f);
}

static const char head[] = "earth_gravity_constant 3.986004415e+14\n"
                           "radius 6378136.3\n"
                           "max_degree 2\n"
                           "norm fully_normalized\n"
                           "end_of_head ====\n";

static myb_status load_text(const char* text, myb_gravity_field** f, myb_diag* diag) {
  *f = NULL;
  const myb_status st = myb_gravity_field_load_icgem(text, strlen(text), 0, NULL, f, diag);
  myb_gravity_field_destroy(*f);
  *f = NULL;
  return st;
}

static void test_icgem_errors(void) {
  myb_gravity_field* f = NULL;
  myb_diag diag = {0};
  char buf[1024];
  snprintf(buf, sizeof buf, "%sgfc 2 0 -4.8e-4 0 0 0\n", head);
  MT_CHECK(load_text(buf, &f, &diag) == MYB_OK);
  MT_CHECK(load_text("radius 1\nmax_degree 2\nnorm fully_normalized\nend_of_head\n", &f,
                     &diag) == MYB_E_FORMAT); /* no GM */
  MT_CHECK(load_text("earth_gravity_constant 1\nradius 1\nmax_degree 2\nnorm unnormalized\n"
                     "end_of_head\n",
                     &f, &diag) == MYB_E_FORMAT);
  MT_CHECK(load_text("earth_gravity_constant 1\nradius 1\nmax_degree 2\nnorm fully_normalized\n",
                     &f, &diag) == MYB_E_FORMAT); /* no end_of_head */
  snprintf(buf, sizeof buf, "%sgfc 2 0 -4.8e-4\n", head);
  MT_CHECK(load_text(buf, &f, &diag) == MYB_E_FORMAT);
  MT_CHECK(diag.line == 6);
  snprintf(buf, sizeof buf, "%sgfc 2 3 1e-6 0 0 0\n", head);
  MT_CHECK(load_text(buf, &f, &diag) == MYB_E_FORMAT); /* m > n */
  snprintf(buf, sizeof buf, "%sgfc 2 0 x 0 0 0\n", head);
  MT_CHECK(load_text(buf, &f, &diag) == MYB_E_FORMAT);
  snprintf(buf, sizeof buf, "%sgfct 2 0 -4.8e-4 0 0 0\n", head);
  MT_CHECK(load_text(buf, &f, &diag) == MYB_E_FORMAT); /* no epoch */
  snprintf(buf, sizeof buf, "%sgfct 2 0 -4.8e-4 0 0 0 20100101\nacos 2 0 1e-11 0 0 0\n", head);
  MT_CHECK(load_text(buf, &f, &diag) == MYB_E_FORMAT); /* no period */
  snprintf(buf, sizeof buf, "%sgfc 2 0 -4.8e-4 0 0 0\ntrnd 2 0 1e-11 0 0 0\n", head);
  MT_CHECK(load_text(buf, &f, &diag) == MYB_E_FORMAT); /* trend without gfct */
  /* coefficients above max_degree are ignored */
  snprintf(buf, sizeof buf, "%sgfc 3 0 1e-6 0 0 0\n", head);
  MT_CHECK(load_text(buf, &f, &diag) == MYB_OK);
  MT_CHECK(myb_gravity_field_load_icgem_file(MYB_TEST_DATA_DIR "/grav/missing.gfc", 0, NULL, &f,
                                             &diag) != MYB_OK);
}

static myb_time epoch(int y, int mo, int d, int h, int mi) {
  const myb_calendar c = {y, mo, d, h, mi, 0.0};
  myb_time t;
  MT_CHECK(myb_calendar_to_time(&c, &t) == MYB_OK);
  return t;
}

static void test_icgem_time_variable(void) {
  myb_gravity_field* f = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_gravity_field_load_icgem_file(MYB_TEST_DATA_DIR "/grav/synthetic2.gfc", 0, NULL,
                                             &f, &diag) == MYB_OK);
  if (f == NULL) {
    printf("%s\n", diag.message);
    return;
  }
  MT_CHECK(strcmp(myb_gravity_field_tide_system(f), "zero_tide") == 0);
  const double gfct20 = -4.841694947139e-04;
  const double trnd20 = -2.154662425317e-11;
  const double acos20 = 8.331958055607e-11;
  const double asin20 = 7.523339889891e-11;
  double c = 0.0;
  double s = 0.0;
  /* at t0: static + acos */
  myb_gravity_field_set_time(f, epoch(2010, 1, 1, 0, 0));
  MT_CHECK(myb_gravity_field_coefficient(f, 2, 0, &c, &s) == MYB_OK);
  MT_CHECK(fabs(c - (gfct20 + acos20)) < 1e-26);
  /* a quarter of a year later: + trend / 4 + asin */
  myb_time t = epoch(2010, 1, 1, 0, 0);
  t = myb_time_add(t, 0.25 * 365.25 * 86400.0);
  myb_gravity_field_set_time(f, t);
  MT_CHECK(myb_gravity_field_coefficient(f, 2, 0, &c, &s) == MYB_OK);
  MT_CHECK(fabs(c - (gfct20 + 0.25 * trnd20 + asin20)) < 1e-26);
  /* ten years later: the trend dominates */
  t = myb_time_add(epoch(2010, 1, 1, 0, 0), 10.0 * 365.25 * 86400.0);
  myb_gravity_field_set_time(f, t);
  MT_CHECK(myb_gravity_field_coefficient(f, 2, 0, &c, &s) == MYB_OK);
  MT_CHECK(fabs(c - (gfct20 + 10.0 * trnd20 + acos20)) < 1e-25);
  /* C21 with epoch 2010-01-01 12:00, half-year period: at 12:00 + 1/8 year the phase is pi/2 */
  t = myb_time_add(epoch(2010, 1, 1, 12, 0), 0.125 * 365.25 * 86400.0);
  myb_gravity_field_set_time(f, t);
  MT_CHECK(myb_gravity_field_coefficient(f, 2, 1, &c, &s) == MYB_OK);
  MT_CHECK(fabs(c - (-3.383273346550e-10 + 0.125 * 1.0e-12)) < 1e-24);
  MT_CHECK(fabs(s - (1.453697346412e-09 - 0.125 * 2.0e-12)) < 1e-24);
  /* corrections last until the next set_time */
  myb_gravity_field_set_time(f, epoch(2010, 1, 1, 0, 0));
  MT_CHECK(myb_gravity_field_add_correction(f, 3, 3, 1e-9, -2e-9) == MYB_OK);
  MT_CHECK(myb_gravity_field_add_correction(f, 2, 0, 1e-10, 0.0) == MYB_OK);
  MT_CHECK(myb_gravity_field_coefficient(f, 3, 3, &c, &s) == MYB_OK);
  MT_CHECK_SAME(c, 7.212144036300e-07 + 1e-9);
  MT_CHECK_SAME(s, 1.414349261929e-06 - 2e-9);
  MT_CHECK(myb_gravity_field_add_correction(f, 4, 0, 1e-9, 0.0) == MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_field_add_correction(f, 2, 0, (double)NAN, 0.0) == MYB_E_DOMAIN);
  myb_gravity_field_set_time(f, epoch(2010, 1, 1, 0, 0));
  MT_CHECK(myb_gravity_field_coefficient(f, 3, 3, &c, &s) == MYB_OK);
  MT_CHECK_SAME(c, 7.212144036300e-07);
  MT_CHECK(myb_gravity_field_coefficient(f, 2, 0, &c, &s) == MYB_OK);
  MT_CHECK(fabs(c - (gfct20 + acos20)) < 1e-26);
  myb_gravity_field_destroy(f);
}

/* ---- Closed forms ---------------------------------------------------------------------------- */

static void test_closed_forms(void) {
  const double gm = 3.986004418e14;
  const double R = 6378137.0;
  const double J2 = 1.0826359e-3;
  double c[6] = {1.0, 0.0, 0.0, -J2 / sqrt(5.0), 0.0, 0.0};
  double s[6] = {0.0};
  myb_gravity_field* point = NULL;
  myb_gravity_field* j2 = NULL;
  myb_gravity_work* w = NULL;
  MT_CHECK(myb_gravity_field_create(gm, R, 0, c, s, NULL, &point) == MYB_OK);
  MT_CHECK(myb_gravity_field_create(gm, R, 2, c, s, NULL, &j2) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(2, NULL, &w) == MYB_OK);
  mt_rng rng = {7};
  double worst_point = 0.0;
  double worst_j2 = 0.0;
  for (int k = 0; k < 200; ++k) {
    double r[3];
    for (int i = 0; i < 3; ++i) {
      r[i] = R * (double)mt_int(&rng, -3000000, 3000000) / 1e6;
    }
    const double rr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    if (rr < R) {
      continue;
    }
    double v = 0.0;
    double a[3];
    MT_CHECK(myb_gravity_field_eval(point, w, r, 0, -1, &v, a) == MYB_OK);
    const double ap[3] = {-gm * r[0] / (rr * rr * rr), -gm * r[1] / (rr * rr * rr),
                          -gm * r[2] / (rr * rr * rr)};
    worst_point = fmax(worst_point, fmax(rel(a, ap), fabs(v - gm / rr) / (gm / rr)));
    /* J2 alone (n_min = 2) against the closed form */
    MT_CHECK(myb_gravity_field_eval(j2, w, r, 2, -1, &v, a) == MYB_OK);
    const double z2 = r[2] * r[2] / (rr * rr);
    const double kj = 1.5 * J2 * gm * R * R / (rr * rr * rr * rr * rr);
    const double aj[3] = {kj * r[0] * (5.0 * z2 - 1.0), kj * r[1] * (5.0 * z2 - 1.0),
                          kj * r[2] * (5.0 * z2 - 3.0)};
    const double vj = -0.5 * J2 * gm * R * R / (rr * rr * rr) * (3.0 * z2 - 1.0);
    /* the potential relative to its scale J2 GM R^2 / r^3 (3 z^2 - 1 crosses zero) */
    worst_j2 = fmax(worst_j2,
                    fmax(rel(a, aj), fabs(v - vj) / (J2 * gm * R * R / (rr * rr * rr))));
    /* n_min / n_max split: (0..1) + (2..2) = all */
    double v01 = 0.0;
    double a01[3];
    double v2 = 0.0;
    double a2[3];
    double va = 0.0;
    double aa[3];
    MT_CHECK(myb_gravity_field_eval(j2, w, r, 0, 1, &v01, a01) == MYB_OK);
    MT_CHECK(myb_gravity_field_eval(j2, w, r, 2, 2, &v2, a2) == MYB_OK);
    MT_CHECK(myb_gravity_field_eval(j2, w, r, 0, -1, &va, aa) == MYB_OK);
    const double sum[3] = {a01[0] + a2[0], a01[1] + a2[1], a01[2] + a2[2]};
    MT_CHECK(rel(sum, aa) < 4e-16);
  }
  printf("closed forms: point mass %.2e, J2 %.2e relative\n", worst_point, worst_j2);
  MT_CHECK(worst_point < 7e-16 && worst_j2 < 2e-15);
  /* poles: finite, continuous */
  double a_pole[3];
  double a_near[3];
  double v = 0.0;
  MT_CHECK(myb_gravity_field_eval(j2, w, (const double[3]){0.0, 0.0, -R}, 0, -1, &v, a_pole) ==
           MYB_OK);
  MT_CHECK(myb_gravity_field_eval(j2, w, (const double[3]){1e-3, 0.0, -R}, 0, -1, &v, a_near) ==
           MYB_OK);
  MT_CHECK(a_pole[0] == 0.0 && a_pole[1] == 0.0 && rel(a_near, a_pole) < 1e-9);
  /* errors */
  MT_CHECK(myb_gravity_field_eval(j2, w, (const double[3]){0.0, 0.0, 0.0}, 0, -1, &v, a_pole) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_field_eval(j2, w, (const double[3]){(double)INFINITY, 0.0, 0.0}, 0, -1, &v,
                                  a_pole) == MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_field_eval(j2, w, (const double[3]){R, 0.0, 0.0}, 0, 3, &v, a_pole) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_field_eval(j2, w, (const double[3]){R, 0.0, 0.0}, 4, -1, &v, a_pole) ==
           MYB_E_DOMAIN);
  myb_gravity_work* small = NULL;
  MT_CHECK(myb_gravity_work_create(1, NULL, &small) == MYB_OK);
  MT_CHECK(myb_gravity_field_eval(j2, small, (const double[3]){R, 0.0, 0.0}, 0, -1, &v, a_pole) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_field_eval(j2, small, (const double[3]){R, 0.0, 0.0}, 0, 1, NULL, NULL) ==
           MYB_OK);
  myb_gravity_work_destroy(small);
  MT_CHECK(myb_gravity_field_create(-1.0, R, 2, c, s, NULL, &point) == MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_field_create(gm, R, 2701, c, s, NULL, &point) == MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_work_create(-1, NULL, &small) == MYB_E_DOMAIN);
  c[4] = (double)NAN;
  MT_CHECK(myb_gravity_field_create(gm, R, 2, c, s, NULL, &point) == MYB_E_DOMAIN);
  myb_gravity_field_destroy(j2);
  myb_gravity_field_destroy(point);
  myb_gravity_work_destroy(w);
  myb_gravity_field_destroy(NULL);
  myb_gravity_work_destroy(NULL);
}

/* The gradient of the potential by central differences (an independent check of the
 * derivative formulas at full degree with the real field). */
static void test_gradient_egm2008(void) {
  myb_gravity_field* f = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_gravity_field_load_icgem_file(MYB_TEST_DATA_DIR "/grav/EGM2008_n20.gfc", 0, NULL,
                                             &f, &diag) == MYB_OK);
  myb_gravity_work* w = NULL;
  MT_CHECK(myb_gravity_work_create(20, NULL, &w) == MYB_OK);
  const double pts[4][3] = {
      {4.1e6, -3.3e6, 3.9e6}, {0.0, 1.0, 6.6e6}, {-7.0e6, 1.0e5, -2.0e3}, {1.5e6, 2.5e6, -6.0e6}};
  double worst = 0.0;
  for (int k = 0; k < 4; ++k) {
    double a[3];
    MT_CHECK(myb_gravity_field_eval(f, w, pts[k], 2, -1, NULL, a) == MYB_OK);
    double fd[3];
    for (int i = 0; i < 3; ++i) {
      double p[3] = {pts[k][0], pts[k][1], pts[k][2]};
      const double h = 2.0;
      double vp = 0.0;
      double vm = 0.0;
      p[i] += h;
      MT_CHECK(myb_gravity_field_eval(f, w, p, 2, -1, &vp, NULL) == MYB_OK);
      p[i] -= 2.0 * h;
      MT_CHECK(myb_gravity_field_eval(f, w, p, 2, -1, &vm, NULL) == MYB_OK);
      fd[i] = (vp - vm) / (2.0 * h);
    }
    worst = fmax(worst, rel(fd, a));
  }
  printf("EGM2008 n<=20, degrees >= 2: gradient vs differences of the potential %.1e\n", worst);
  MT_CHECK(worst < 1e-8);
  myb_gravity_work_destroy(w);
  myb_gravity_field_destroy(f);
}

/* ---- Polyhedron ------------------------------------------------------------------------------ */

static const double box[3][2] = {{-1.0, 2.0}, {-0.5, 1.5}, {-2.0, 0.7}};
static const size_t box_faces[12][3] = {{0, 2, 3}, {0, 3, 1}, {4, 5, 7}, {4, 7, 6},
                                        {0, 1, 5}, {0, 5, 4}, {2, 6, 7}, {2, 7, 3},
                                        {0, 4, 6}, {0, 6, 2}, {1, 3, 7}, {1, 7, 5}};

static void box_vertices(double v[8][3]) {
  for (int k = 0; k < 8; ++k) {
    v[k][0] = box[0][k & 1];
    v[k][1] = box[1][(k >> 1) & 1];
    v[k][2] = box[2][(k >> 2) & 1];
  }
}

static void test_polyhedron_prism(void) {
  double v[8][3];
  box_vertices(v);
  myb_polyhedron* p = NULL;
  MT_CHECK(myb_polyhedron_create(&v[0][0], 8, &box_faces[0][0], 12, 1.0, NULL, &p) == MYB_OK);
  double center[3];
  MT_CHECK(fabs(myb_polyhedron_volume(p, center) - 3.0 * 2.0 * 2.7) < 1e-14);
  MT_CHECK(fabs(center[0] - 0.5) < 1e-15 && fabs(center[1] - 0.5) < 1e-15 &&
           fabs(center[2] + 0.65) < 1e-15);
  /* closed form, 40 digits (G rho = 1) */
  static const struct {
    double r[3];
    double v;
    double a[3];
    int inside;
  } ref[] = {
      {{3.1, -2.2, 1.7}, 3.662261170497838, {-0.47248005806037906, 0.52466333295984846,
                                             -0.43596638331423144}, 0},
      {{0.3, 0.2, -0.4}, 14.588222107797312, {0.58867896241474771, 1.7440570738968916,
                                              -0.89374899198448437}, 1},
      {{1.9, 1.4, 0.69}, 8.0321881816689897, {-2.6408628026377004, -2.2967207941881204,
                                              -3.3259175753470448}, 1},
      {{-1.25, 0.1, -0.3}, 9.1518124606982928, {4.7980173375269149, 1.0552568436867686,
                                                -0.63742288946247598}, 0},
      {{40.0, -25.0, 31.0}, 0.28584059181545919, {-0.0035146864623785232, 0.0022698578583676437,
                                                  -0.00281657250214161}, 0},
      {{0.5, 0.5, 0.71}, 11.113588621497331, {0.0, 0.0, -6.3947819458214915}, 0},
  };
  double worst_v = 0.0;
  double worst_a = 0.0;
  for (size_t k = 0; k < sizeof ref / sizeof ref[0]; ++k) {
    double u = 0.0;
    double a[3];
    double lap = 0.0;
    MT_CHECK(myb_polyhedron_eval(p, ref[k].r, &u, a, &lap) == MYB_OK);
    /* rounding grows like (r / size)^2 (see grav.h) */
    const double rr = sqrt(ref[k].r[0] * ref[k].r[0] + ref[k].r[1] * ref[k].r[1] +
                           ref[k].r[2] * ref[k].r[2]);
    const double scale = fmax(1.0, rr * rr / 9.0);
    const double ev = fabs(u - ref[k].v) / ref[k].v / scale;
    const double ea = rel(a, ref[k].a) / scale;
    worst_v = fmax(worst_v, ev);
    worst_a = fmax(worst_a, ea);
    const double expected_lap = ref[k].inside ? -4.0 * 3.141592653589793 : 0.0;
    MT_CHECK(fabs(lap - expected_lap) < 5e-13);
    if (!(fabs(lap - expected_lap) < 5e-13) || ev > 2e-15 || ea > 2e-15) {
      printf("  point %zu: potential %.1e, acceleration %.1e, Laplacian %.3e\n", k, ev, ea, lap);
    }
  }
  printf("polyhedron (prism) vs closed form: potential %.1e, acceleration %.1e relative / "
         "max(1, (r / 3)^2)\n",
         worst_v, worst_a);
  MT_CHECK(worst_v < 2e-15 && worst_a < 2e-15);
  /* far away: a point mass of G M = g_rho V, up to the degree-2 term ~ (size / r)^2 */
  const double far[3] = {300.0, -100.0, 200.0};
  double a[3];
  double u = 0.0;
  MT_CHECK(myb_polyhedron_eval(p, far, &u, a, NULL) == MYB_OK);
  const double d[3] = {far[0] - 0.5, far[1] - 0.5, far[2] + 0.65};
  const double dd = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  MT_CHECK(fabs(u - 16.2 / dd) / (16.2 / dd) < (3.0 / dd) * (3.0 / dd));
  /* on the surface: edges, a vertex, a face diagonal, a face; the potential and the
   * acceleration are continuous there */
  const double surface[5][3] = {
      {0.5, -0.5, -2.0}, {2.0, 1.5, 0.7}, {0.5, 0.5, 0.7}, {0.25, 0.75, 0.7}, {2.0, 0.1, -1.0}};
  for (int k = 0; k < 5; ++k) {
    double us = 0.0;
    double as[3];
    double un = 0.0;
    double an[3];
    double near[3] = {surface[k][0], surface[k][1], surface[k][2]};
    near[0] += 1e-9;
    near[1] += 1e-9;
    near[2] += 1e-9;
    MT_CHECK(myb_polyhedron_eval(p, surface[k], &us, as, NULL) == MYB_OK);
    MT_CHECK(myb_polyhedron_eval(p, near, &un, an, NULL) == MYB_OK);
    MT_CHECK(fabs(us - un) < 1e-8 && rel(as, an) < 1e-6);
  }
  MT_CHECK(myb_polyhedron_eval(p, (const double[3]){(double)NAN, 0.0, 0.0}, &u, a, NULL) ==
           MYB_E_DOMAIN);
  myb_polyhedron_destroy(p);
  myb_polyhedron_destroy(NULL);
}

static void test_polyhedron_validation(void) {
  double v[8][3];
  box_vertices(v);
  size_t faces[12][3];
  memcpy(faces, box_faces, sizeof faces);
  myb_polyhedron* p = NULL;
  /* open mesh */
  MT_CHECK(myb_polyhedron_create(&v[0][0], 8, &faces[0][0], 11, 1.0, NULL, &p) == MYB_E_DOMAIN);
  /* one face reversed */
  faces[3][1] = 7;
  faces[3][2] = 4;
  faces[3][0] = 6;
  MT_CHECK(myb_polyhedron_create(&v[0][0], 8, &faces[0][0], 12, 1.0, NULL, &p) == MYB_E_DOMAIN);
  /* all reversed: inward orientation */
  for (int f = 0; f < 12; ++f) {
    faces[f][0] = box_faces[f][0];
    faces[f][1] = box_faces[f][2];
    faces[f][2] = box_faces[f][1];
  }
  MT_CHECK(myb_polyhedron_create(&v[0][0], 8, &faces[0][0], 12, 1.0, NULL, &p) == MYB_E_DOMAIN);
  /* index out of range, degenerate face */
  memcpy(faces, box_faces, sizeof faces);
  faces[5][2] = 8;
  MT_CHECK(myb_polyhedron_create(&v[0][0], 8, &faces[0][0], 12, 1.0, NULL, &p) == MYB_E_DOMAIN);
  memcpy(faces, box_faces, sizeof faces);
  faces[5][2] = faces[5][1];
  MT_CHECK(myb_polyhedron_create(&v[0][0], 8, &faces[0][0], 12, 1.0, NULL, &p) == MYB_E_DOMAIN);
  memcpy(faces, box_faces, sizeof faces);
  MT_CHECK(myb_polyhedron_create(&v[0][0], 8, &faces[0][0], 12, -1.0, NULL, &p) == MYB_E_DOMAIN);
  v[3][1] = (double)INFINITY;
  MT_CHECK(myb_polyhedron_create(&v[0][0], 8, &faces[0][0], 12, 1.0, NULL, &p) == MYB_E_DOMAIN);
}

/* ---- Two extended bodies ------------------------------------------------------------------- */

enum { CLUSTER = 5 };

/* A body of 2 CLUSTER point masses at +-p_k (no odd moments), G = 1. */
typedef struct cluster {
  double m[2 * CLUSTER];
  double x[2 * CLUSTER][3];
  double mass;
  double q[9]; /* sum m (3 x x^T - |x|^2 1) */
} cluster;

static void make_cluster(cluster* c, mt_rng* rng, const double offset[3]) {
  c->mass = 0.0;
  memset(c->q, 0, sizeof c->q);
  for (int k = 0; k < CLUSTER; ++k) {
    const double mk = 0.5 + (double)mt_int(rng, 0, 1000) / 1000.0;
    double p[3];
    for (int i = 0; i < 3; ++i) {
      p[i] = (double)mt_int(rng, -1000, 1000) / 1000.0;
    }
    for (int s = 0; s < 2; ++s) {
      const int j = 2 * k + s;
      c->m[j] = mk;
      for (int i = 0; i < 3; ++i) {
        c->x[j][i] = offset[i] + (s == 0 ? p[i] : -p[i]);
      }
      c->mass += mk;
      const double pp = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
      for (int i = 0; i < 3; ++i) {
        for (int l = 0; l < 3; ++l) {
          c->q[3 * i + l] += mk * (3.0 * p[i] * p[l] - (i == l ? pp : 0.0));
        }
      }
    }
  }
}

/* Force on B from the points of A at the points (or the centre) of B. */
static void force_points(const cluster* a, const double (*y)[3], const double* my, int ny,
                         double f[3]) {
  f[0] = f[1] = f[2] = 0.0;
  for (int i = 0; i < 2 * CLUSTER; ++i) {
    for (int j = 0; j < ny; ++j) {
      const double d[3] = {a->x[i][0] - y[j][0], a->x[i][1] - y[j][1], a->x[i][2] - y[j][2]};
      const double r = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      const double k = a->m[i] * my[j] / (r * r * r);
      for (int l = 0; l < 3; ++l) {
        f[l] += k * d[l];
      }
    }
  }
}

static void test_figure_figure(void) {
  mt_rng rng = {91};
  double worst_with[2] = {0.0, 0.0};
  double worst_without[2] = {0.0, 0.0};
  for (int trial = 0; trial < 20; ++trial) {
    for (int far = 0; far < 2; ++far) {
      const double distance = far ? 24.0 : 12.0;
      cluster a;
      cluster b;
      const double origin[3] = {0.0, 0.0, 0.0};
      make_cluster(&a, &rng, origin);
      double dir[3];
      for (int i = 0; i < 3; ++i) {
        dir[i] = (double)mt_int(&rng, -1000, 1000) / 1000.0;
      }
      const double dn = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
      const double r[3] = {distance * dir[0] / dn, distance * dir[1] / dn, distance * dir[2] / dn};
      make_cluster(&b, &rng, r);
      /* exact */
      double exact[3];
      force_points(&a, (const double(*)[3])b.x, b.m, 2 * CLUSTER, exact);
      /* field of A at B's centre (all degrees), of B at A's centre, minus the point term
       * counted twice */
      double fa[3];
      const double centre_b[1][3] = {{r[0], r[1], r[2]}};
      const double mass_b[1] = {b.mass};
      force_points(&a, centre_b, mass_b, 1, fa);
      double fb[3];
      const double centre_a[1][3] = {{0.0, 0.0, 0.0}};
      const double mass_a[1] = {a.mass};
      force_points(&b, centre_a, mass_a, 1, fb); /* force on A from B */
      const double rr = distance;
      double expansion[3];
      for (int i = 0; i < 3; ++i) {
        const double point = -a.mass * b.mass * r[i] / (rr * rr * rr);
        expansion[i] = fa[i] - fb[i] - point;
      }
      const double pn = a.mass * b.mass / (rr * rr);
      double e0 = 0.0;
      for (int i = 0; i < 3; ++i) {
        e0 += (expansion[i] - exact[i]) * (expansion[i] - exact[i]);
      }
      /* + U22: accelerations times masses (G = 1) */
      double acc_a[3] = {0.0, 0.0, 0.0};
      double acc_b[3] = {0.0, 0.0, 0.0};
      MT_CHECK(myb_figure_figure(a.mass, a.q, b.mass, b.q, r, acc_a, acc_b) == MYB_OK);
      double e1 = 0.0;
      for (int i = 0; i < 3; ++i) {
        MT_CHECK(fabs(acc_a[i] * a.mass + acc_b[i] * b.mass) < 1e-15 * pn); /* action-reaction */
        const double d = expansion[i] + acc_b[i] * b.mass - exact[i];
        e1 += d * d;
      }
      worst_without[far] = fmax(worst_without[far], sqrt(e0) / pn);
      worst_with[far] = fmax(worst_with[far], sqrt(e1) / pn);
    }
  }
  printf("figure-figure, point clusters at 12 and 24 sizes: without U22 %.1e -> %.1e (ratio "
         "%.1f), with U22 %.1e -> %.1e (ratio %.1f)\n",
         worst_without[0], worst_without[1], worst_without[0] / worst_without[1], worst_with[0],
         worst_with[1], worst_with[0] / worst_with[1]);
  /* the remaining terms fall like r^-8 against r^-6 for U22 (force, relative to the point) */
  MT_CHECK(worst_with[0] < 0.15 * worst_without[0] && worst_with[1] < 0.04 * worst_without[1]);
  MT_CHECK(worst_with[0] / worst_with[1] > 30.0);
  static const double a_zero_q[9] = {0.0};
  MT_CHECK(myb_figure_figure(1.0, a_zero_q, 1.0, a_zero_q, (const double[3]){0, 0, 0}, NULL,
                             NULL) == MYB_E_DOMAIN);
  /* quadrupole from coefficients: the degree-2 potential of a field equals G n.Q.n / (2 r^3) */
  const double c2[3] = {-4.84165e-4, 2.4e-10, 2.43938e-6};
  const double s2[2] = {1.4e-9, -1.40027e-6};
  double coef_c[6] = {1.0, 0.0, 0.0, c2[0], c2[1], c2[2]};
  double coef_s[6] = {0.0, 0.0, 0.0, 0.0, s2[0], s2[1]};
  myb_gravity_field* f = NULL;
  myb_gravity_work* w = NULL;
  MT_CHECK(myb_gravity_field_create(3.986e14, 6.378e6, 2, coef_c, coef_s, NULL, &f) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(2, NULL, &w) == MYB_OK);
  double gq[9];
  myb_quadrupole_from_coefficients(3.986e14, 6.378e6, c2, s2, gq);
  MT_CHECK(fabs(gq[0] + gq[4] + gq[8]) < 3e-16 * (fabs(gq[0]) + fabs(gq[4]) + fabs(gq[8])));
  const double x[3] = {3.1e6, -5.2e6, 4.4e6};
  const double xr = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  double v2 = 0.0;
  MT_CHECK(myb_gravity_field_eval(f, w, x, 2, 2, &v2, NULL) == MYB_OK);
  double qx[3];
  for (int i = 0; i < 3; ++i) {
    qx[i] = gq[3 * i] * x[0] + gq[3 * i + 1] * x[1] + gq[3 * i + 2] * x[2];
  }
  const double vq = (x[0] * qx[0] + x[1] * qx[1] + x[2] * qx[2]) / (2.0 * xr * xr * xr * xr * xr);
  printf("quadrupole from coefficients vs the field's degree 2: %.1e relative\n",
         fabs(vq - v2) / fabs(v2));
  MT_CHECK(fabs(vq - v2) < 1e-14 * fabs(v2));
  /* rotation: Q of rotated points */
  const double ang = 0.7;
  const double m[9] = {cos(ang), sin(ang), 0.0, -sin(ang), cos(ang), 0.0, 0.0, 0.0, 1.0};
  double gq_in[9];
  myb_quadrupole_rotate(m, gq, gq_in);
  /* n.Q_in.n = (M n).Q_body.(M n) */
  double mn[3];
  for (int i = 0; i < 3; ++i) {
    mn[i] = m[3 * i] * x[0] + m[3 * i + 1] * x[1] + m[3 * i + 2] * x[2];
  }
  double a1 = 0.0;
  double a2 = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      a1 += x[i] * gq_in[3 * i + j] * x[j];
      a2 += mn[i] * gq[3 * i + j] * mn[j];
    }
  }
  MT_CHECK(fabs(a1 - a2) < 1e-14 * fabs(a2));
  myb_gravity_work_destroy(w);
  myb_gravity_field_destroy(f);
}

/* Mascons against their own spherical-harmonic expansion (exact coefficients of point masses,
 * C_nm - i S_nm = 1 / (M (2n + 1)) sum m (r/R)^n Pbar_nm(sin phi) e^(-i m lambda)), converging
 * outside the sphere that holds them; the mascon sum is the reference of the expansion. */
static void test_mascons(void) {
  enum { K = 12, N = 40 };
  double gm[K];
  double pos[3 * K];
  mt_rng rng = {5};
  double total = 0.0;
  for (int k = 0; k < K; ++k) {
    gm[k] = 1e9 * (1.0 + (double)mt_int(&rng, 0, 1000) / 1000.0);
    total += gm[k];
    for (int i = 0; i < 3; ++i) {
      pos[3 * k + i] = (double)mt_int(&rng, -1000, 1000); /* within 1.8 km */
    }
  }
  /* coefficients by the recursion of the field itself: a field of one mascon at a time */
  static double c[(N + 1) * (N + 2) / 2];
  static double s[(N + 1) * (N + 2) / 2];
  const double R = 2000.0;
  for (int k = 0; k < K; ++k) {
    const double* p = pos + 3 * k;
    const double r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    const double t = p[2] / r;
    const double u = sqrt(p[0] * p[0] + p[1] * p[1]) / r;
    const double lam = atan2(p[1], p[0]);
    /* Pbar_nm by the standard column recursion (independent of myb's) */
    double pmm = 1.0;
    for (int m = 0; m <= N; ++m) {
      if (m == 1) {
        pmm = sqrt(3.0) * u;
      } else if (m > 1) {
        pmm *= sqrt((2.0 * m + 1.0) / (2.0 * m)) * u;
      }
      double p0 = pmm;
      double p1 = m < N ? sqrt(2.0 * m + 3.0) * t * pmm : 0.0;
      for (int n = m; n <= N; ++n) {
        double pn;
        if (n == m) {
          pn = p0;
        } else if (n == m + 1) {
          pn = p1;
        } else {
          const double a = sqrt((2.0 * n - 1.0) * (2.0 * n + 1.0) / ((double)(n - m) * (n + m)));
          const double b = sqrt((2.0 * n + 1.0) * (n + m - 1.0) * (n - m - 1.0) /
                                ((double)(n - m) * (n + m) * (2.0 * n - 3.0)));
          pn = a * t * p1 - b * p0;
          p0 = p1;
          p1 = pn;
        }
        const double w = gm[k] / total * pow(r / R, n) / (2.0 * n + 1.0) * pn;
        c[n * (n + 1) / 2 + m] += w * cos(m * lam);
        s[n * (n + 1) / 2 + m] += w * sin(m * lam);
      }
    }
  }
  myb_gravity_field* f = NULL;
  myb_gravity_work* w = NULL;
  MT_CHECK(myb_gravity_field_create(total, R, N, c, s, NULL, &f) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(N, NULL, &w) == MYB_OK);
  double worst = 0.0;
  for (int k = 0; k < 20; ++k) {
    double x[3];
    for (int i = 0; i < 3; ++i) {
      x[i] = (double)mt_int(&rng, -1000, 1000);
    }
    const double xn = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    for (int i = 0; i < 3; ++i) {
      x[i] *= 3.0 * R / xn; /* at 3 R: the series to 40 converges to ~1e-19 */
    }
    double vm = 0.0;
    double am[3];
    double vs = 0.0;
    double as[3];
    MT_CHECK(myb_mascons_eval(gm, pos, K, x, &vm, am) == MYB_OK);
    MT_CHECK(myb_gravity_field_eval(f, w, x, 0, -1, &vs, as) == MYB_OK);
    worst = fmax(worst, fmax(fabs(vm - vs) / vm, rel(as, am)));
  }
  printf("mascons vs their spherical-harmonic expansion to degree %d at 3 R: %.1e\n", N, worst);
  MT_CHECK(worst < 1e-14);
  MT_CHECK(myb_mascons_eval(gm, pos, K, pos + 3, NULL, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_mascons_eval(NULL, NULL, 0, pos, NULL, NULL) == MYB_OK);
  myb_gravity_work_destroy(w);
  myb_gravity_field_destroy(f);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_figure_figure();
  test_mascons();
  test_icgem_egm2008();
  test_icgem_errors();
  test_icgem_time_variable();
  test_closed_forms();
  test_gradient_egm2008();
  test_polyhedron_prism();
  test_polyhedron_validation();
  return mt_finish("test_grav");
}
