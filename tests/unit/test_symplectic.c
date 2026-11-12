/* SPDX-License-Identifier: BSD-3-Clause */
/* Symplectic compositions: order conditions of the weights, observed orders 2, 4, 6, 8 (the
 * triple jump and the schemes of 7 and 15 stages), the error of the new schemes against the
 * triple jump for equal work, bounded energy error over 10^4 orbits (no drift, unlike a
 * non-symplectic method), time reversibility, a time-dependent force, error statuses. */
#include <math.h>

#include "core/fpguard.h"
#include "integ/symplectic_coefficients.h"
#include "mtest.h"
#include "myb/integ.h"
#include "myb/kepler.h"

static const double pi = 3.141592653589793;
static const double mu = 1.0; /* units: a = 1, period 2 pi */

static void test_weights(void) {
  /* sum c = 1; a symmetric composition of order p >= 4 needs sum c^3 = 0, order 6 also
   * sum c^5 = 0, order 8 also sum c^7 = 0 (with further conditions met by construction). */
  const struct {
    const double* w;
    int count;
    int order;
  } sets[] = {{myb_symplectic_order4, 3, 4},      {myb_symplectic_order6, 9, 6},
              {myb_symplectic_order8, 27, 8},     {myb_symplectic_order6_s7, 7, 6},
              {myb_symplectic_order8_s15, 15, 8}};
  for (int k = 0; k < 5; ++k) {
    double s1 = 0.0;
    double s3 = 0.0;
    double s5 = 0.0;
    double s7 = 0.0;
    double size = 0.0;
    for (int i = 0; i < sets[k].count; ++i) {
      const double c = sets[k].w[i];
      MT_CHECK(c == sets[k].w[sets[k].count - 1 - i]); /* symmetric */
      s1 += c;
      s3 += c * c * c;
      s5 += c * c * c * c * c;
      s7 += c * c * c * c * c * c * c;
      size = fmax(size, fabs(c));
    }
    MT_CHECK(fabs(s1 - 1.0) < 1e-14);
    MT_CHECK(fabs(s3) < 1e-14 * size * size * size);
    if (sets[k].order >= 6) {
      MT_CHECK(fabs(s5) < 1e-14 * pow(size, 5));
    }
    if (sets[k].order >= 8) {
      MT_CHECK(fabs(s7) < 1e-14 * pow(size, 7));
    }
  }
}

static myb_status kepler(void* context, double t, const double* x, double* a) {
  (void)context;
  (void)t;
  const double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  const double k = -mu / (r * r * r);
  for (int i = 0; i < 3; ++i) {
    a[i] = k * x[i];
  }
  return MYB_OK;
}

static myb_pv orbit(double e) {
  const myb_elements el = {.p = 1.0 - e * e, .e = e, .i = 0.6, .raan = 1.1, .argp = 0.4, .nu = 0.3};
  myb_pv pv = {};
  MT_CHECK(myb_pv_from_elements(mu, &el, &pv) == MYB_OK);
  return pv;
}

static double energy(const double* y) {
  return 0.5 * (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]) -
         mu / sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
}

static myb_symplectic* make(int order) {
  myb_symplectic* s = NULL;
  MT_CHECK(myb_symplectic_create(3, order, NULL, &s) == MYB_OK);
  return s;
}

static void test_orders(void) {
  /* One orbit of e = 0.3 with N steps against the analytic solution. */
  const myb_pv start = orbit(0.3);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  myb_pv exact = {};
  MT_CHECK(myb_kepler_propagate(mu, &start, 2.0 * pi, &exact) == MYB_OK);
  const int orders[6] = {2, 4, 6, 8, MYB_SYMPLECTIC_6_S7, MYB_SYMPLECTIC_8_S15};
  const int expected[6] = {2, 4, 6, 8, 6, 8};
  const int base[6] = {400, 200, 100, 50, 100, 50}; /* above round-off at the finest step */
  double finest[6] = {0.0};
  for (int q = 0; q < 6; ++q) {
    double previous = 0.0;
    double order_seen = 0.0;
    for (int k = 0; k < 3; ++k) {
      const int steps = base[q] << k;
      myb_symplectic* s = make(orders[q]);
      MT_CHECK(myb_symplectic_start(s, kepler, NULL, 0.0, y0) == MYB_OK);
      MT_CHECK(myb_symplectic_integrate(s, 2.0 * pi, 2.0 * pi / steps * 1.0000001) == MYB_OK);
      MT_CHECK(myb_symplectic_statistics(s).steps == steps);
      const double* y = myb_symplectic_state(s);
      const double err = hypot(hypot(y[0] - exact.r[0], y[1] - exact.r[1]), y[2] - exact.r[2]);
      if (k > 0) {
        order_seen = log2(previous / err);
      }
      previous = err;
      myb_symplectic_destroy(s);
    }
    printf("order %d%s: observed %.2f (error %.2e at %d steps per orbit)\n", expected[q],
           q < 4 ? "" : (q == 4 ? ", 7 stages" : ", 15 stages"), order_seen, previous,
           base[q] << 2);
    MT_CHECK(fabs(order_seen - expected[q]) < 0.3);
    finest[q] = previous;
  }
  /* for equal work: the error of the new scheme at the same steps, scaled to as many
   * evaluations as the triple jump spends, (stages ratio)^p */
  const double gain6 = finest[2] / (finest[4] * pow(9.0 / 7.0, 6.0));
  const double gain8 = finest[3] / (finest[5] * pow(27.0 / 15.0, 8.0));
  printf("for equal work the error of the triple jump is %.0f times (order 6) and %.0f times "
         "(order 8) that of the 7- and 15-stage schemes\n",
         gain6, gain8);
  MT_CHECK(gain6 > 1.4 && gain8 > 100.0);
}

static void test_energy(void) {
  /* 5000 orbits of e = 0.5 at order 6 with 200 steps per orbit: the energy error oscillates
   * with the orbit and stays bounded; it does not grow from the first 100 orbits to the last. */
  const myb_pv start = orbit(0.5);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  const double e0 = energy(y0);
  myb_symplectic* s = make(6);
  MT_CHECK(myb_symplectic_start(s, kepler, NULL, 0.0, y0) == MYB_OK);
  const double h = 2.0 * pi / 200.0;
  double early = 0.0;
  double late = 0.0;
  for (long i = 1; i <= 1000000; ++i) {
    if (myb_symplectic_step(s, h) != MYB_OK) {
      MT_CHECK(0);
      break;
    }
    if (i <= 20000 || i > 980000) {
      const double d = fabs(energy(myb_symplectic_state(s)) / e0 - 1.0);
      if (i <= 20000) {
        early = fmax(early, d);
      } else {
        late = fmax(late, d);
      }
    }
  }
  printf("order 6, 200 steps per orbit: max relative energy error %.2e (orbits 1-100), %.2e "
         "(orbits 4901-5000)\n",
         early, late);
  MT_CHECK(late < 1.1 * early); /* bounded: no growth over 5000 orbits */
  MT_CHECK(late < 1e-6);        /* the triple jump has large error constants */
  myb_symplectic_destroy(s);
}

static void test_reversibility(void) {
  const myb_pv start = orbit(0.7);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  myb_symplectic* s = make(8);
  MT_CHECK(myb_symplectic_start(s, kepler, NULL, 0.0, y0) == MYB_OK);
  const double h = 2.0 * pi / 300.0;
  for (int i = 0; i < 3000; ++i) {
    MT_CHECK(myb_symplectic_step(s, h) == MYB_OK);
  }
  for (int i = 0; i < 3000; ++i) {
    MT_CHECK(myb_symplectic_step(s, -h) == MYB_OK);
  }
  const double* y = myb_symplectic_state(s);
  double worst = 0.0;
  for (int k = 0; k < 6; ++k) {
    worst = fmax(worst, fabs(y[k] - y0[k]));
  }
  printf("10 orbits forward and back (order 8): largest difference %.2e\n", worst);
  MT_CHECK(worst < 1e-11); /* round-off of 81000 stages, sheared by the orbit */
  MT_CHECK(myb_symplectic_time(s) == 0.0 || fabs(myb_symplectic_time(s)) < 1e-15);
  myb_symplectic_destroy(s);
}

/* x'' = -x + cos(2 t): x = cos t + (cos t - cos 2t) / 3 for x(0) = 1, v(0) = 0. */
static myb_status driven(void* context, double t, const double* x, double* a) {
  (void)context;
  a[0] = -x[0] + cos(2.0 * t);
  return MYB_OK;
}

static void test_time_dependent(void) {
  myb_symplectic* s = NULL;
  MT_CHECK(myb_symplectic_create(1, 6, NULL, &s) == MYB_OK);
  const double y0[2] = {1.0, 0.0};
  MT_CHECK(myb_symplectic_start(s, driven, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_symplectic_integrate(s, 10.0, 0.01) == MYB_OK);
  MT_CHECK(myb_symplectic_time(s) == 10.0);
  const double exact = cos(10.0) + (cos(10.0) - cos(20.0)) / 3.0;
  printf("driven oscillator (order 6, h = 0.01): error %.2e\n",
         fabs(myb_symplectic_state(s)[0] - exact));
  MT_CHECK(fabs(myb_symplectic_state(s)[0] - exact) < 1e-12);
  myb_symplectic_destroy(s);
}

static myb_status failing(void* context, double t, const double* x, double* a) {
  (void)context;
  (void)x;
  a[0] = 0.0;
  return t > 0.5 ? MYB_E_DOMAIN : MYB_OK;
}

static void test_errors(void) {
  myb_symplectic* s = NULL;
  MT_CHECK(myb_symplectic_create(3, 5, NULL, &s) == MYB_E_DOMAIN);
  MT_CHECK(myb_symplectic_create(0, 4, NULL, &s) == MYB_E_DOMAIN);
  MT_CHECK(myb_symplectic_create(1, 4, NULL, &s) == MYB_OK);
  const double y0[2] = {0.0, 1.0};
  MT_CHECK(myb_symplectic_start(s, failing, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_symplectic_step(s, 0.0) == MYB_E_DOMAIN);
  MT_CHECK(myb_symplectic_integrate(s, 1.0, -0.1) == MYB_E_DOMAIN);
  MT_CHECK(myb_symplectic_integrate(s, 1.0, 0.1) == MYB_E_DOMAIN); /* from the callback */
  myb_symplectic_destroy(s);
  myb_symplectic_destroy(NULL);
}

int main(void) {
  test_weights();
  test_orders();
  test_energy();
  test_reversibility();
  test_time_dependent();
  test_errors();
  return mt_finish("test_symplectic");
}
