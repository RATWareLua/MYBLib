/* SPDX-License-Identifier: BSD-3-Clause */
/* Regularisation: the KS transformation (round trip, bilinear relation, both branches), Kepler
 * orbits of eccentricity 0.9 to 0.999 over 20 revolutions against the closed form for Cowell,
 * KS and Sundman with DOP853 at the same tolerance (accuracy and work), a J2-perturbed orbit
 * against IAS15, and the passage of a rectilinear orbit through the centre. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/integ.h"
#include "myb/kepler.h"
#include "myb/regular.h"

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static void test_transformation(void) {
  mt_rng rng = {9};
  double worst_x = 0.0;
  double worst_v = 0.0;
  double worst_bil = 0.0;
  double worst_h = 0.0;
  for (int k = 0; k < 2000; ++k) {
    double x[3];
    double v[3];
    for (int i = 0; i < 3; ++i) {
      x[i] = ((double)(mt_next(&rng) >> 11) * 0x1p-53 - 0.5) * 2e7;
      v[i] = ((double)(mt_next(&rng) >> 11) * 0x1p-53 - 0.5) * 2e4;
    }
    double y[10];
    myb_ks_from_cartesian(3.986e14, 12.5, x, v, y);
    double x2[3];
    double v2[3];
    double t;
    myb_ks_to_cartesian(y, &t, x2, v2);
    MT_CHECK(t == 12.5);
    for (int i = 0; i < 3; ++i) {
      worst_x = fmax(worst_x, fabs(x2[i] - x[i]) / norm3(x));
      worst_v = fmax(worst_v, fabs(v2[i] - v[i]) / norm3(v));
    }
    const double* u = y;
    const double* up = y + 4;
    const double bil = u[3] * up[0] - u[2] * up[1] + u[1] * up[2] - u[0] * up[3];
    worst_bil = fmax(worst_bil, fabs(bil) / (norm3(x) * norm3(v)));
    const double r = norm3(x);
    worst_h = fmax(worst_h, fabs(y[8] - (3.986e14 / r - 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))));
  }
  printf("KS transformation: x %.2e, v %.2e, bilinear %.2e (relative)\n", worst_x, worst_v,
         worst_bil);
  MT_CHECK(worst_x < 5e-16 && worst_v < 5e-16 && worst_bil < 2e-16);
  (void)worst_h;
}

typedef struct cowell_context {
  double mu;
  double j2r2; /* J2 R^2 */
} cowell_context;

static myb_status j2(void* context, double t, const double x[3], const double v[3], double a[3]) {
  (void)t;
  (void)v;
  const cowell_context* c = context;
  const double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  const double r = sqrt(r2);
  const double z2 = x[2] * x[2] / r2;
  const double f = -1.5 * c->mu * c->j2r2 / (r2 * r2 * r);
  a[0] = f * x[0] * (1.0 - 5.0 * z2);
  a[1] = f * x[1] * (1.0 - 5.0 * z2);
  a[2] = f * x[2] * (3.0 - 5.0 * z2);
  return MYB_OK;
}

static myb_status cowell_rhs(void* context, double t, const double* y, double* dy) {
  const cowell_context* c = context;
  const double r = norm3(y);
  const double k = c->mu / (r * r * r);
  double a[3] = {0, 0, 0};
  if (c->j2r2 != 0.0) {
    (void)j2(context, t, y, y + 3, a);
  }
  for (int i = 0; i < 3; ++i) {
    dy[i] = y[3 + i];
    dy[3 + i] = a[i] - k * y[i];
  }
  return MYB_OK;
}

static myb_status ias_accel(void* context, double t, const double* x, const double* v, double* a) {
  const cowell_context* c = context;
  const double r = norm3(x);
  const double k = c->mu / (r * r * r);
  (void)j2(context, t, x, v, a);
  for (int i = 0; i < 3; ++i) {
    a[i] -= k * x[i];
  }
  return MYB_OK;
}

static myb_status cowell(const cowell_context* c, const myb_dop853_options* o, const double x0[3],
                         const double v0[3], double t1, double x1[3], double v1[3],
                         int64_t* evaluations) {
  myb_dop853* d = NULL;
  myb_status st = myb_dop853_create(6, o, NULL, &d);
  if (st != MYB_OK) {
    return st;
  }
  double y0[6] = {x0[0], x0[1], x0[2], v0[0], v0[1], v0[2]};
  cowell_context cc = *c;
  st = myb_dop853_start(d, cowell_rhs, &cc, 0.0, y0);
  if (st == MYB_OK) {
    st = myb_dop853_integrate(d, t1);
  }
  const double* y = myb_dop853_state(d);
  memcpy(x1, y, 3 * sizeof *x1);
  memcpy(v1, y + 3, 3 * sizeof *v1);
  *evaluations = myb_dop853_statistics(d).evaluations;
  myb_dop853_destroy(d);
  return st;
}

static void test_kepler(void) {
  const double mu = 1.0;
  const double eccentricities[3] = {0.9, 0.99, 0.999};
  const myb_dop853_options o = {.rtol = 1e-12, .atol = 1e-15};
  for (int k = 0; k < 3; ++k) {
    const double e = eccentricities[k];
    /* start at apocentre, a = 1: 20 periods of 2 pi */
    const double ra = 1.0 + e;
    const double x0[3] = {ra * 0.6, ra * 0.8, 0.0};
    const double va = sqrt(mu * (1.0 - e) / ra);
    const double v0[3] = {-0.8 * va * 0.9, 0.6 * va * 0.9, 0.4359 * va};
    const myb_pv pv0 = {{x0[0], x0[1], x0[2]}, {v0[0], v0[1], v0[2]}};
    const double t1 = 20.0 * 2.0 * 3.141592653589793;
    myb_pv exact;
    MT_CHECK(myb_kepler_propagate(mu, &pv0, t1, &exact) == MYB_OK);
    const double scale = norm3(exact.r);
    double x[3];
    double v[3];
    int64_t nc = 0;
    int64_t nk = 0;
    int64_t ns = 0;
    const cowell_context cc = {.mu = mu};
    MT_CHECK(cowell(&cc, &o, x0, v0, t1, x, v, &nc) == MYB_OK);
    double ec = 0.0;
    for (int i = 0; i < 3; ++i) {
      ec = fmax(ec, fabs(x[i] - exact.r[i]) / scale);
    }
    const myb_regular_problem p = {.mu = mu, .alpha = 1.0};
    MT_CHECK(myb_regular_propagate(MYB_REGULAR_KS, &p, &o, 0.0, x0, v0, t1, x, v, &nk) == MYB_OK);
    double ek = 0.0;
    for (int i = 0; i < 3; ++i) {
      ek = fmax(ek, fabs(x[i] - exact.r[i]) / scale);
    }
    MT_CHECK(myb_regular_propagate(MYB_REGULAR_SUNDMAN, &p, &o, 0.0, x0, v0, t1, x, v, &ns) ==
             MYB_OK);
    double es = 0.0;
    for (int i = 0; i < 3; ++i) {
      es = fmax(es, fabs(x[i] - exact.r[i]) / scale);
    }
    printf("e = %.3f, 20 periods, DOP853 1e-12: Cowell %.2e (%lld f), KS %.2e (%lld f), Sundman "
           "%.2e (%lld f)\n",
           e, ec, (long long)nc, ek, (long long)nk, es, (long long)ns);
    /* KS: the same work and accuracy whatever the eccentricity */
    MT_CHECK(ek < 3e-9);
    MT_CHECK(nk < 30000);
    if (e >= 0.99) {
      MT_CHECK(nk * 3 < nc); /* at least three times less work than Cowell */
    }
  }
}

static void test_perturbed(void) {
  /* J2 orbit of the Earth, e = 0.7, 3 days: KS and Sundman against IAS15 at epsilon 1e-12 */
  const double mu = 3.986004418e14;
  cowell_context cc = {.mu = mu, .j2r2 = 1.08263e-3 * 6378137.0 * 6378137.0};
  const double rp = 7000e3;
  const double e = 0.7;
  const double a = rp / (1.0 - e);
  const double vp = sqrt(mu * (1.0 + e) / rp);
  const double x0[3] = {rp, 0.0, 0.0};
  const double v0[3] = {0.0, vp * cos(0.9), vp * sin(0.9)};
  const double t1 = 3.0 * 86400.0;
  (void)a;
  const myb_dop853_options o = {.rtol = 1e-14, .atol = 1e-6};
  const myb_regular_problem p = {.mu = mu, .perturbation = j2, .context = &cc, .alpha = 1.0};
  double xk[3];
  double vk[3];
  double xs[3];
  double vs[3];
  int64_t nk = 0;
  int64_t ns = 0;
  MT_CHECK(myb_regular_propagate(MYB_REGULAR_KS, &p, &o, 0.0, x0, v0, t1, xk, vk, &nk) == MYB_OK);
  MT_CHECK(myb_regular_propagate(MYB_REGULAR_SUNDMAN, &p, &o, 0.0, x0, v0, t1, xs, vs, &ns) ==
           MYB_OK);
  double xc[3];
  double vc[3];
  int64_t nc = 0;
  /* the reference: IAS15 at its round-off limit */
  myb_ias15* ias = NULL;
  const myb_ias15_options io = {.epsilon = 1e-12};
  MT_CHECK(myb_ias15_create(3, &io, NULL, &ias) == MYB_OK);
  double y0[6] = {x0[0], x0[1], x0[2], v0[0], v0[1], v0[2]};
  MT_CHECK(myb_ias15_start(ias, ias_accel, &cc, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(ias, t1) == MYB_OK);
  memcpy(xc, myb_ias15_state(ias), sizeof xc);
  memcpy(vc, myb_ias15_state(ias) + 3, sizeof vc);
  nc = myb_ias15_statistics(ias).evaluations;
  myb_ias15_destroy(ias);
  double dk = 0.0;
  double ds = 0.0;
  for (int i = 0; i < 3; ++i) {
    dk = fmax(dk, fabs(xk[i] - xc[i]));
    ds = fmax(ds, fabs(xs[i] - xc[i]));
  }
  printf("J2, e = 0.7, 3 days, DOP853 1e-14: KS vs IAS15 %.2e m (%lld f), Sundman %.2e m "
         "(%lld f); IAS15 %lld f\n",
         dk, (long long)nk, ds, (long long)ns, (long long)nc);
  /* three independent integrations agree to 1 mm in 23000 km after 7 revolutions (4e-11) */
  MT_CHECK(dk < 3e-3);
  MT_CHECK(ds < 3e-3);
}

static void test_collision(void) {
  /* rectilinear fall from rest at x = 1 (mu = 1): reaches the centre at t = pi / (2 sqrt 2) and
   * bounces back to x = 1 at t = pi / sqrt 2 (the regularised continuation) */
  const myb_regular_problem p = {.mu = 1.0, .alpha = 1.0};
  const myb_dop853_options o = {.rtol = 1e-13, .atol = 1e-16};
  const double x0[3] = {1.0, 0.0, 0.0};
  const double v0[3] = {0.0, 0.0, 0.0};
  double x[3];
  double v[3];
  const double period = 3.141592653589793 / sqrt(2.0);
  MT_CHECK(myb_regular_propagate(MYB_REGULAR_KS, &p, &o, 0.0, x0, v0, period, x, v, NULL) ==
           MYB_OK);
  printf("rectilinear orbit through the centre after one bounce: x = %.15f (1), |v| = %.2e\n",
         x[0], norm3(v));
  MT_CHECK(fabs(x[0] - 1.0) < 1e-11 && fabs(x[1]) < 1e-12 && norm3(v) < 1e-6);
  /* half way: at the centre */
  MT_CHECK(myb_regular_propagate(MYB_REGULAR_KS, &p, &o, 0.0, x0, v0, 0.5 * period, x, v, NULL) ==
           MYB_OK);
  /* x ~ (9 mu / 2)^(1/3) dt^(2/3) near the centre: 1e-14 of time gives 1e-9 */
  printf("at the collision time: |x| = %.2e\n", norm3(x));
  MT_CHECK(norm3(x) < 1e-8);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_transformation();
  test_kepler();
  test_perturbed();
  test_collision();
  return mt_finish("test_regular");
}
