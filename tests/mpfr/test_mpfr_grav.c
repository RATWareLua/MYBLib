/* SPDX-License-Identifier: BSD-3-Clause */
/* Spherical-harmonic fields against MPFR (256 bits): the potential summed with the fully
 * normalised Legendre functions in multiple precision (their recursion checked against closed
 * forms), the acceleration as central differences of that potential (independent of myb's
 * derivative formulas). Degree 120 with random coefficients at the surface and above, at the
 * equator, at a pole and next to it; degree 2190 at latitude 62 deg, where the sectorial seeds
 * u^m underflow without the 2^930 scaling while the columns still matter. */
#include <math.h>
#include <mpfr.h>
#include <stdio.h>
#include <stdlib.h>

#include "mtest.h"
#include "myb/grav.h"

enum { P = 256 };

typedef struct model {
  int degree;
  double gm;
  double radius;
  double* c;
  double* s;
} model;

static size_t idx(int n, int m) {
  return (size_t)n * (size_t)(n + 1) / 2 + (size_t)m;
}

/* V at (x, y, z) in MPFR; out = V. */
static void potential(const model* md, mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t out) {
  mpfr_t r, t, u, rho, ratio, pw, cl, sl, cm, sm, tmp, tmp2, sum, colsum;
  mpfr_inits2(P, r, t, u, rho, ratio, pw, cl, sl, cm, sm, tmp, tmp2, sum, colsum, (mpfr_ptr)0);
  mpfr_hypot(rho, x, y, MPFR_RNDN);
  mpfr_hypot(r, rho, z, MPFR_RNDN);
  mpfr_div(t, z, r, MPFR_RNDN);
  mpfr_div(u, rho, r, MPFR_RNDN);
  if (mpfr_zero_p(rho)) {
    mpfr_set_ui(cl, 1, MPFR_RNDN);
    mpfr_set_ui(sl, 0, MPFR_RNDN);
  } else {
    mpfr_div(cl, x, rho, MPFR_RNDN);
    mpfr_div(sl, y, rho, MPFR_RNDN);
  }
  mpfr_set_d(ratio, md->radius, MPFR_RNDN);
  mpfr_div(ratio, ratio, r, MPFR_RNDN);
  const int N = md->degree;
  mpfr_t* col = malloc((size_t)(N + 2) * sizeof(mpfr_t));
  for (int n = 0; n <= N + 1; ++n) {
    mpfr_init2(col[n], P);
  }
  mpfr_t seed;
  mpfr_init2(seed, P);
  mpfr_set_ui(seed, 1, MPFR_RNDN);
  mpfr_set_ui(cm, 1, MPFR_RNDN);
  mpfr_set_ui(sm, 0, MPFR_RNDN);
  mpfr_set_ui(sum, 0, MPFR_RNDN);
  for (int m = 0; m <= N; ++m) {
    if (m == 1) {
      mpfr_sqrt_ui(tmp, 3, MPFR_RNDN);
      mpfr_mul(seed, seed, tmp, MPFR_RNDN);
      mpfr_mul(seed, seed, u, MPFR_RNDN);
    } else if (m >= 2) {
      mpfr_set_ui(tmp, 2 * (unsigned long)m + 1, MPFR_RNDN);
      mpfr_div_ui(tmp, tmp, 2 * (unsigned long)m, MPFR_RNDN);
      mpfr_sqrt(tmp, tmp, MPFR_RNDN);
      mpfr_mul(seed, seed, tmp, MPFR_RNDN);
      mpfr_mul(seed, seed, u, MPFR_RNDN);
    }
    if (m >= 1) { /* cos(m lambda), sin(m lambda) */
      mpfr_mul(tmp, cm, cl, MPFR_RNDN);
      mpfr_mul(tmp2, sm, sl, MPFR_RNDN);
      mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
      mpfr_mul(tmp2, sm, cl, MPFR_RNDN);
      mpfr_fma(sm, cm, sl, tmp2, MPFR_RNDN);
      mpfr_set(cm, tmp, MPFR_RNDN);
    }
    mpfr_set(col[m], seed, MPFR_RNDN);
    if (m + 1 <= N) {
      mpfr_sqrt_ui(tmp, 2 * (unsigned long)m + 3, MPFR_RNDN);
      mpfr_mul(tmp, tmp, t, MPFR_RNDN);
      mpfr_mul(col[m + 1], tmp, seed, MPFR_RNDN);
    }
    for (int n = m + 2; n <= N; ++n) {
      const double dn = n;
      const double dm = m;
      /* a = sqrt((2n-1)(2n+1)/((n-m)(n+m))), b = sqrt((2n+1)(n+m-1)(n-m-1)/((n-m)(n+m)(2n-3))) */
      mpfr_set_d(tmp, (2 * dn - 1) * (2 * dn + 1), MPFR_RNDN);
      mpfr_div_d(tmp, tmp, (dn - dm) * (dn + dm), MPFR_RNDN);
      mpfr_sqrt(tmp, tmp, MPFR_RNDN);
      mpfr_mul(tmp, tmp, t, MPFR_RNDN);
      mpfr_mul(tmp, tmp, col[n - 1], MPFR_RNDN);
      mpfr_set_d(tmp2, (2 * dn + 1) * (dn + dm - 1) * (dn - dm - 1), MPFR_RNDN);
      mpfr_div_d(tmp2, tmp2, (dn - dm) * (dn + dm) * (2 * dn - 3), MPFR_RNDN);
      mpfr_sqrt(tmp2, tmp2, MPFR_RNDN);
      mpfr_mul(tmp2, tmp2, col[n - 2], MPFR_RNDN);
      mpfr_sub(col[n], tmp, tmp2, MPFR_RNDN);
    }
    /* sum_n ratio^n Pbar_nm (C cos + S sin) */
    mpfr_set_ui(colsum, 0, MPFR_RNDN);
    mpfr_pow_ui(pw, ratio, (unsigned long)m, MPFR_RNDN);
    for (int n = m; n <= N; ++n) {
      const size_t i = idx(n, m);
      mpfr_mul_d(tmp, cm, md->c[i], MPFR_RNDN);
      mpfr_mul_d(tmp2, sm, md->s[i], MPFR_RNDN);
      mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
      mpfr_mul(tmp, tmp, col[n], MPFR_RNDN);
      mpfr_mul(tmp, tmp, pw, MPFR_RNDN);
      mpfr_add(colsum, colsum, tmp, MPFR_RNDN);
      mpfr_mul(pw, pw, ratio, MPFR_RNDN);
    }
    mpfr_add(sum, sum, colsum, MPFR_RNDN);
  }
  mpfr_set_d(tmp, md->gm, MPFR_RNDN);
  mpfr_div(tmp, tmp, r, MPFR_RNDN);
  mpfr_mul(out, sum, tmp, MPFR_RNDN);
  for (int n = 0; n <= N + 1; ++n) {
    mpfr_clear(col[n]);
  }
  free(col);
  mpfr_clears(seed, r, t, u, rho, ratio, pw, cl, sl, cm, sm, tmp, tmp2, sum, colsum, (mpfr_ptr)0);
}

/* Reference potential and acceleration (central differences, h = 1e-25 r). */
static void reference(const model* md, const double xyz[3], double* v, double a[3]) {
  mpfr_t x[3], h, vp, vm, out;
  for (int k = 0; k < 3; ++k) {
    mpfr_init2(x[k], P);
    mpfr_set_d(x[k], xyz[k], MPFR_RNDN);
  }
  mpfr_inits2(P, h, vp, vm, out, (mpfr_ptr)0);
  potential(md, x[0], x[1], x[2], out);
  *v = mpfr_get_d(out, MPFR_RNDN);
  const double r = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
  for (int k = 0; k < 3; ++k) {
    mpfr_set_d(h, r, MPFR_RNDN);
    mpfr_mul_d(h, h, 1e-25, MPFR_RNDN);
    mpfr_add(x[k], x[k], h, MPFR_RNDN);
    potential(md, x[0], x[1], x[2], vp);
    mpfr_sub(x[k], x[k], h, MPFR_RNDN);
    mpfr_sub(x[k], x[k], h, MPFR_RNDN);
    potential(md, x[0], x[1], x[2], vm);
    mpfr_add(x[k], x[k], h, MPFR_RNDN);
    mpfr_sub(out, vp, vm, MPFR_RNDN);
    mpfr_div(out, out, h, MPFR_RNDN);
    mpfr_div_ui(out, out, 2, MPFR_RNDN);
    a[k] = mpfr_get_d(out, MPFR_RNDN);
  }
  for (int k = 0; k < 3; ++k) {
    mpfr_clear(x[k]);
  }
  mpfr_clears(h, vp, vm, out, (mpfr_ptr)0);
}

static void random_model(model* md, int degree, mt_rng* rng) {
  md->degree = degree;
  md->gm = 3.986004415e14;
  md->radius = 6378136.3;
  const size_t count = idx(degree + 1, 0);
  md->c = calloc(count, sizeof(double));
  md->s = calloc(count, sizeof(double));
  md->c[0] = 1.0;
  for (int n = 2; n <= degree; ++n) {
    for (int m = 0; m <= n; ++m) {
      const double scale = 1e-5 / ((double)n * n);
      md->c[idx(n, m)] = scale * ((double)mt_int(rng, -1000000, 1000000) / 1e6);
      md->s[idx(n, m)] = m == 0 ? 0.0 : scale * ((double)mt_int(rng, -1000000, 1000000) / 1e6);
    }
  }
  md->c[idx(2, 0)] = -4.841651437908150e-04;
}

static void check_points(const model* md, const double (*pts)[3], int count, double* worst_v,
                         double* worst_a) {
  myb_gravity_field* field = NULL;
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_field_create(md->gm, md->radius, md->degree, md->c, md->s, NULL, &field) ==
           MYB_OK);
  MT_CHECK(myb_gravity_work_create(md->degree, NULL, &work) == MYB_OK);
  for (int k = 0; k < count; ++k) {
    double v = 0.0;
    double a[3];
    double rv = 0.0;
    double ra[3];
    MT_CHECK(myb_gravity_field_eval(field, work, pts[k], 0, -1, &v, a) == MYB_OK);
    reference(md, pts[k], &rv, ra);
    const double na = sqrt(ra[0] * ra[0] + ra[1] * ra[1] + ra[2] * ra[2]);
    const double ea = sqrt((a[0] - ra[0]) * (a[0] - ra[0]) + (a[1] - ra[1]) * (a[1] - ra[1]) +
                           (a[2] - ra[2]) * (a[2] - ra[2])) / na;
    *worst_v = fmax(*worst_v, fabs(v - rv) / fabs(rv));
    *worst_a = fmax(*worst_a, ea);
  }
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
}

static void test_legendre_closed_forms(void) {
  /* the reference recursion: Pbar_20, Pbar_21, Pbar_22, Pbar_31, Pbar_33 through the potential of
   * single-coefficient models at a point */
  const double xyz[3] = {3.1e6, -4.2e6, 5.3e6};
  const double r = sqrt(3.1e6 * 3.1e6 + 4.2e6 * 4.2e6 + 5.3e6 * 5.3e6);
  const double t = 5.3e6 / r;
  const double u = sqrt(3.1e6 * 3.1e6 + 4.2e6 * 4.2e6) / r;
  const double lam = atan2(-4.2e6, 3.1e6);
  const struct {
    int n, m;
    double p;
  } cases[] = {{2, 0, sqrt(5.0) * (3 * t * t - 1) / 2}, {2, 1, sqrt(15.0) * t * u},
               {2, 2, sqrt(15.0) / 2 * u * u}, {3, 1, sqrt(42.0) / 4 * u * (5 * t * t - 1)},
               {3, 3, sqrt(70.0) / 4 * u * u * u}};
  for (size_t k = 0; k < sizeof cases / sizeof cases[0]; ++k) {
    model md = {.degree = 3, .gm = 1.0, .radius = r};
    md.c = calloc(idx(4, 0), sizeof(double));
    md.s = calloc(idx(4, 0), sizeof(double));
    md.c[idx(cases[k].n, cases[k].m)] = 1.0;
    double v = 0.0;
    double a[3];
    reference(&md, xyz, &v, a);
    const double expected = cases[k].p * cos(cases[k].m * lam) / r;
    MT_CHECK(fabs(v - expected) < 1e-15 * fabs(expected) + 1e-30);
    free(md.c);
    free(md.s);
  }
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_legendre_closed_forms();
  mt_rng rng = {120};
  model md;
  random_model(&md, 120, &rng);
  const double R = md.radius;
  double pts[16][3] = {
      {R, 0, 0},           {0, 0, R * 1.001},     {0, 0, -2.0 * R},  {1e-10 * R, 0, R},
      {R * 0.7, R * 0.7, 0}, {-R, 2 * R, 1.5 * R}, {1.2 * R, -0.3 * R, 0.4 * R},
      {0, -3 * R, 0.2 * R}, {2e3, -1e3, R * 1.0003}, {R * 0.3, R * 0.2, -R * 0.93}};
  double worst_v = 0.0;
  double worst_a = 0.0;
  check_points(&md, (const double(*)[3])pts, 10, &worst_v, &worst_a);
  printf("degree 120, 10 points (pole, next to it, equator, surface .. 3 R): potential %.2e, "
         "acceleration %.2e relative\n",
         worst_v, worst_a);
  MT_CHECK(worst_v < 5e-16 && worst_a < 3e-15);
  free(md.c);
  free(md.s);
  model big;
  random_model(&big, 2190, &rng);
  const double lat = 62.0 * 3.141592653589793 / 180.0;
  double far[1][3] = {{R * 1.0002 * cos(lat) * cos(0.7), R * 1.0002 * cos(lat) * sin(0.7),
                       R * 1.0002 * sin(lat)}};
  double wv = 0.0;
  double wa = 0.0;
  check_points(&big, (const double(*)[3])far, 1, &wv, &wa);
  printf("degree 2190 at latitude 62 deg: potential %.2e, acceleration %.2e relative\n", wv, wa);
  MT_CHECK(wv < 5e-16 && wa < 1e-14);
  free(big.c);
  free(big.s);
  return mt_finish("test_mpfr_grav");
}
