/* SPDX-License-Identifier: BSD-3-Clause */
/* Force models, variational equations, covariances, random numbers:
 *  - the state transition matrix of the variational equations against the closed form of the
 *    two-body problem (myb_kepler_stm), with the analytic Jacobian of a point mass and with
 *    central differences;
 *  - sensitivities to a parameter against differences of whole propagations;
 *  - a covariance propagated linearly against 4000 Monte Carlo samples;
 *  - xoshiro256** against its reference values, normal deviates, the model's bookkeeping. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/integ.h"
#include "myb/kepler.h"
#include "myb/model.h"
#include "myb/num.h"
#include "myb/uncert.h"

static const double gm_earth = 3.986004418e14;

typedef struct central {
  double gm;
} central;

static myb_status central_accel(void* context, double t, const double x[3], const double v[3],
                                double a[3]) {
  (void)t;
  (void)v;
  const central* c = context;
  const double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  const double k = -c->gm / (r * r * r);
  for (int i = 0; i < 3; ++i) {
    a[i] = k * x[i];
  }
  return MYB_OK;
}

static myb_status central_jacobian(void* context, double t, const double x[3], const double v[3],
                                   double dadx[9], double dadv[9]) {
  (void)t;
  (void)v;
  const central* c = context;
  myb_point_mass_jacobian(c->gm, x, dadx);
  memset(dadv, 0, 9 * sizeof(double));
  return MYB_OK;
}

/* A constant acceleration along a fixed direction scaled by a parameter (like a pressure). */
typedef struct push {
  double k;
} push;

static myb_status push_accel(void* context, double t, const double x[3], const double v[3],
                             double a[3]) {
  (void)t;
  (void)x;
  const push* p = context;
  /* and a velocity-dependent part (drag-like) to exercise d a / d v */
  const double s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  a[0] = p->k * 0.6 - 1e-9 * s * v[0];
  a[1] = p->k * 0.8 - 1e-9 * s * v[1];
  a[2] = -1e-9 * s * v[2];
  return MYB_OK;
}

static double max_rel_stm(const double* a, const double* b) {
  /* blockwise relative error (the blocks differ in units) */
  double worst = 0.0;
  for (int bi = 0; bi < 2; ++bi) {
    for (int bj = 0; bj < 2; ++bj) {
      double nb = 0.0;
      double nd = 0.0;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          const int k = (3 * bi + i) * 6 + 3 * bj + j;
          nb = fmax(nb, fabs(b[k]));
          nd = fmax(nd, fabs(a[k] - b[k]));
        }
      }
      worst = fmax(worst, nd / nb);
    }
  }
  return worst;
}

static void test_two_body_stm(void) {
  central c = {gm_earth};
  const myb_pv pv0 = {{7.0e6, 1.0e5, -3.0e5}, {-200.0, 7.3e3, 1.9e3}};
  const double y0[6] = {pv0.r[0], pv0.r[1], pv0.r[2], pv0.v[0], pv0.v[1], pv0.v[2]};
  for (int analytic = 0; analytic < 2; ++analytic) {
    myb_force_model* m = NULL;
    MT_CHECK(myb_force_model_create(NULL, &m) == MYB_OK);
    MT_CHECK(myb_force_model_add(m, "earth", central_accel, analytic ? central_jacobian : NULL,
                                 &c, NULL) == MYB_OK);
    myb_variational* var = NULL;
    MT_CHECK(myb_variational_create(m, NULL, NULL, &var) == MYB_OK);
    MT_CHECK(myb_variational_start(var, 0.0, y0) == MYB_OK);
    double worst = 0.0;
    for (int k = 1; k <= 5; ++k) {
      const double t = 4000.0 * k; /* about 4 orbits at the end */
      MT_CHECK(myb_variational_integrate(var, t) == MYB_OK);
      double y[6];
      double phi[36];
      myb_variational_get(var, y, phi, NULL);
      myb_pv pv;
      double ref[36];
      MT_CHECK(myb_kepler_stm(gm_earth, &pv0, t, &pv, ref) == MYB_OK);
      worst = fmax(worst, max_rel_stm(phi, ref));
    }
    printf("two-body STM, variational equations vs closed form over 4 orbits (%s Jacobian): "
           "%.1e\n",
           analytic ? "analytic" : "central-difference", worst);
    MT_CHECK(worst < (analytic ? 1e-10 : 1e-8));
    myb_variational_destroy(var);
    myb_force_model_destroy(m);
  }
}

static void propagate(myb_force_model* m, const double y0[6], double t, double y[6]) {
  myb_ias15* in = NULL;
  const myb_ias15_options o = {.epsilon = 1e-10};
  MT_CHECK(myb_ias15_create(3, &o, NULL, &in) == MYB_OK);
  MT_CHECK(myb_ias15_start(in, myb_force_model_accel, m, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(in, t) == MYB_OK);
  memcpy(y, myb_ias15_state(in), 6 * sizeof(double));
  myb_ias15_destroy(in);
}

static void test_sensitivity(void) {
  central c = {gm_earth};
  push p = {2e-7};
  myb_force_model* m = NULL;
  MT_CHECK(myb_force_model_create(NULL, &m) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "earth", central_accel, central_jacobian, &c, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "push", push_accel, NULL, &p, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(m, "k", &p.k, 1e-9, NULL) == MYB_OK);
  const double y0[6] = {7.0e6, 1.0e5, -3.0e5, -200.0, 7.3e3, 1.9e3};
  myb_variational* var = NULL;
  MT_CHECK(myb_variational_create(m, NULL, NULL, &var) == MYB_OK);
  MT_CHECK(myb_variational_start(var, 0.0, y0) == MYB_OK);
  const double t = 20000.0;
  MT_CHECK(myb_variational_integrate(var, t) == MYB_OK);
  double y[6];
  double phi[36];
  double s[6];
  myb_variational_get(var, y, phi, s);
  /* d y / d k by central differences of two propagations */
  double yp[6];
  double ym[6];
  const double h = 2e-8;
  p.k = 2e-7 + h;
  propagate(m, y0, t, yp);
  p.k = 2e-7 - h;
  propagate(m, y0, t, ym);
  p.k = 2e-7;
  double worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double fd = (yp[i] - ym[i]) / (2.0 * h);
    worst = fmax(worst, fabs(s[i] - fd) / fabs(fd));
  }
  /* and one column of Phi (d y / d vx0) */
  double y1[6];
  double y2[6];
  double q[6];
  memcpy(q, y0, sizeof q);
  q[3] += 1e-3;
  propagate(m, q, t, y1);
  q[3] -= 2e-3;
  propagate(m, q, t, y2);
  double worst_phi = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double fd = (y1[i] - y2[i]) / 2e-3;
    worst_phi = fmax(worst_phi, fabs(phi[i * 6 + 3] - fd) / fabs(fd));
  }
  printf("sensitivity to a pressure-like parameter vs differences of propagations: %.1e; a "
         "column of Phi with drag: %.1e\n",
         worst, worst_phi);
  MT_CHECK(worst < 1e-6 && worst_phi < 1e-6);
  myb_variational_destroy(var);
  myb_force_model_destroy(m);
}

static void test_covariance(void) {
  /* linear propagation of a covariance (closed-form Phi) against Monte Carlo */
  const myb_pv pv0 = {{7.0e6, 1.0e5, -3.0e5}, {-200.0, 7.3e3, 1.9e3}};
  const double sigma[6] = {10.0, 10.0, 10.0, 0.01, 0.01, 0.01};
  double p0[36];
  memset(p0, 0, sizeof p0);
  for (int i = 0; i < 6; ++i) {
    p0[i * 6 + i] = sigma[i] * sigma[i];
  }
  p0[0 * 6 + 4] = p0[4 * 6 + 0] = 0.5 * sigma[0] * sigma[4]; /* a correlation */
  const double t = 6000.0;
  myb_pv pv;
  double phi[36];
  MT_CHECK(myb_kepler_stm(gm_earth, &pv0, t, &pv, phi) == MYB_OK);
  double p[36];
  myb_covariance_transform(6, phi, p0, p);
  double l[36];
  memcpy(l, p0, sizeof l);
  MT_CHECK(myb_cholesky(6, l) == MYB_OK);
  enum { COUNT = 4000 };
  static double samples[COUNT * 6];
  myb_rng rng;
  myb_rng_seed(&rng, 2026);
  const double mean0[6] = {pv0.r[0], pv0.r[1], pv0.r[2], pv0.v[0], pv0.v[1], pv0.v[2]};
  for (int k = 0; k < COUNT; ++k) {
    double y[6];
    myb_rng_gaussian(&rng, 6, mean0, l, y);
    const myb_pv s0 = {{y[0], y[1], y[2]}, {y[3], y[4], y[5]}};
    myb_pv s1;
    MT_CHECK(myb_kepler_propagate(gm_earth, &s0, t, &s1) == MYB_OK);
    memcpy(samples + 6 * k, s1.r, 3 * sizeof(double));
    memcpy(samples + 6 * k + 3, s1.v, 3 * sizeof(double));
  }
  double mean[6];
  double cov[36];
  MT_CHECK(myb_sample_statistics(COUNT, 6, samples, mean, cov) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 6; ++i) {
    worst = fmax(worst, fabs(sqrt(cov[i * 6 + i]) / sqrt(p[i * 6 + i]) - 1.0));
  }
  printf("covariance after %.0f s: Monte Carlo (%d samples) vs linear, sigma ratios within %.1f%%\n",
         t, COUNT, 100.0 * worst);
  MT_CHECK(worst < 0.05); /* sampling error of a sigma: 1/sqrt(2 N) = 1.1 % */
  /* RTN: trace preserved by a rotation */
  double rtn[36];
  MT_CHECK(myb_covariance_rtn(pv.r, pv.v, p, rtn) == MYB_OK);
  double tr0 = 0.0;
  double tr1 = 0.0;
  for (int i = 0; i < 3; ++i) {
    tr0 += p[i * 6 + i];
    tr1 += rtn[i * 6 + i];
  }
  MT_CHECK(fabs(tr0 - tr1) < 1e-12 * tr0);
  /* the along-track position variance dominates after an orbit with a velocity uncertainty */
  MT_CHECK(rtn[1 * 6 + 1] > rtn[0 * 6 + 0] && rtn[1 * 6 + 1] > rtn[2 * 6 + 2]);
}

static void test_rng(void) {
  /* reference: xoshiro256** from state {1, 2, 3, 4} (the algorithm of Blackman & Vigna) */
  myb_rng r = {.s = {1, 2, 3, 4}};
  MT_CHECK(myb_rng_next(&r) == UINT64_C(11520));
  MT_CHECK(myb_rng_next(&r) == UINT64_C(0));
  MT_CHECK(myb_rng_next(&r) == UINT64_C(1509978240));
  MT_CHECK(myb_rng_next(&r) == UINT64_C(1215971899390074240));
  /* normal deviates: moments of 1e6 samples */
  myb_rng_seed(&r, 7);
  double s1 = 0.0;
  double s2 = 0.0;
  double s4 = 0.0;
  const int n = 1000000;
  for (int k = 0; k < n; ++k) {
    const double z = myb_rng_normal(&r);
    s1 += z;
    s2 += z * z;
    s4 += z * z * z * z;
  }
  const double m1 = s1 / n;
  const double m2 = s2 / n;
  const double m4 = s4 / n;
  printf("normal deviates: mean %.4f, variance %.4f, fourth moment %.3f\n", m1, m2, m4);
  MT_CHECK(fabs(m1) < 0.005 && fabs(m2 - 1.0) < 0.005 && fabs(m4 - 3.0) < 0.03);
  /* the same seed gives the same sequence */
  myb_rng a;
  myb_rng b;
  myb_rng_seed(&a, 42);
  myb_rng_seed(&b, 42);
  for (int k = 0; k < 100; ++k) {
    MT_CHECK(mt_bits(myb_rng_normal(&a)) == mt_bits(myb_rng_normal(&b)));
  }
  double u = 0.0;
  for (int k = 0; k < 1000; ++k) {
    u = myb_rng_uniform(&a);
    MT_CHECK(u >= 0.0 && u < 1.0);
  }
}

static void test_model_bookkeeping(void) {
  central c = {gm_earth};
  push p = {1e-7};
  myb_force_model* m = NULL;
  MT_CHECK(myb_force_model_create(NULL, &m) == MYB_OK);
  size_t i0 = 99;
  size_t i1 = 99;
  MT_CHECK(myb_force_model_add(m, "earth", central_accel, NULL, &c, &i0) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "push", push_accel, NULL, &p, &i1) == MYB_OK);
  MT_CHECK(i0 == 0 && i1 == 1 && myb_force_model_count(m) == 2);
  MT_CHECK(myb_force_model_add(m, "push", push_accel, NULL, &p, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_model_add(m, "", push_accel, NULL, &p, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_model_add(m, "a name that is far too long for it", push_accel, NULL, &p,
                               NULL) == MYB_E_DOMAIN);
  size_t found = 0;
  MT_CHECK(myb_force_model_find(m, "push", &found) == MYB_OK && found == 1);
  MT_CHECK(myb_force_model_find(m, "none", &found) == MYB_E_NOT_FOUND);
  const double x[3] = {7e6, 0.0, 0.0};
  const double v[3] = {0.0, 7.5e3, 0.0};
  double a[3];
  double ae[3];
  double ap[3];
  MT_CHECK(myb_force_model_acceleration(m, 0.0, x, v, a) == MYB_OK);
  MT_CHECK(myb_force_model_term(m, 0, 0.0, x, v, ae) == MYB_OK);
  MT_CHECK(myb_force_model_term(m, 1, 0.0, x, v, ap) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(a[i] == ae[i] + ap[i] || fabs(a[i] - (ae[i] + ap[i])) <= 1e-16 * fabs(ae[i]));
  }
  myb_force_model_set_enabled(m, 1, 0);
  MT_CHECK(!myb_force_model_enabled(m, 1));
  MT_CHECK(myb_force_model_acceleration(m, 0.0, x, v, a) == MYB_OK);
  MT_CHECK(a[0] == ae[0] && a[1] == ae[1]);
  MT_CHECK(myb_force_model_add_parameter(m, "k", &p.k, 0.0, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_model_add_parameter(m, "k", &p.k, 1e-9, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(m, "k", &p.k, 1e-9, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_model_parameter_count(m) == 1 && myb_force_model_parameter(m, 0) == &p.k);
  MT_CHECK(myb_force_model_set_scales(m, 0.0, 1.0) == MYB_E_DOMAIN);
  /* numeric Jacobian of the point mass vs analytic */
  myb_force_model_set_enabled(m, 1, 1);
  myb_force_model_set_enabled(m, 1, 0);
  double jx[9];
  double jv[9];
  MT_CHECK(myb_force_model_jacobian(m, 0.0, x, v, jx, jv, NULL) == MYB_OK);
  double ref[9];
  myb_point_mass_jacobian(gm_earth, x, ref);
  double worst = 0.0;
  for (int i = 0; i < 9; ++i) {
    worst = fmax(worst, fabs(jx[i] - ref[i]) / fabs(ref[0]));
  }
  MT_CHECK(worst < 1e-9);
  myb_force_model_destroy(m);
  myb_force_model_destroy(NULL);
  myb_variational_destroy(NULL);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_model_bookkeeping();
  test_rng();
  test_two_body_stm();
  test_sensitivity();
  test_covariance();
  return mt_finish("test_uncert");
}
