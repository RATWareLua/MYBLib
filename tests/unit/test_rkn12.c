/* SPDX-License-Identifier: BSD-3-Clause */
/* RKN12(10): order conditions of the coefficients (in double), observed order 12, adaptive
 * accuracy against the analytic two-body solution, dense output and events, backward
 * integration, intermediate stops, error statuses. */
#include <math.h>

#include "core/fpguard.h"
#include "integ/rkn12_coefficients.h"
#include "mtest.h"
#include "myb/integ.h"
#include "myb/kepler.h"

static const double pi = 3.141592653589793;
static const double mu = 3.986004418e14;

static void test_coefficients(void) {
  const double* c = myb_rkn12_c;
  for (int i = 0; i < MYB_RKN12_STAGES; ++i) {
    double row = 0.0;
    double size = 0.0;
    for (int j = 0; j < i; ++j) {
      row += myb_rkn12_a[i][j];
      size += fabs(myb_rkn12_a[i][j]);
    }
    MT_CHECK(fabs(row - 0.5 * c[i] * c[i]) <= 4e-16 * size); /* rounding of the sum */
    MT_CHECK(fabs(myb_rkn12_b[i] - myb_rkn12_bp[i] * (1.0 - c[i])) < 1e-15);
  }
  for (int k = 0; k <= 11; ++k) {
    double qb = 0.0;
    double qbp = 0.0;
    double qbt = 0.0;
    for (int i = 0; i < MYB_RKN12_STAGES; ++i) {
      const double ck = pow(c[i], k);
      qb += myb_rkn12_b[i] * ck;
      qbp += myb_rkn12_bp[i] * ck;
      qbt += myb_rkn12_bt[i] * ck;
    }
    if (k <= 10) {
      MT_CHECK(fabs(qb - 1.0 / ((k + 1.0) * (k + 2.0))) < 1e-15);
    }
    MT_CHECK(fabs(qbp - 1.0 / (k + 1.0)) < 1e-15);
    if (k <= 8) {
      MT_CHECK(fabs(qbt) < 1e-15); /* the order-10 solution shares these quadratures */
    }
  }
}

static myb_status two_body(void* context, double t, const double* x, double* a) {
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
  const myb_elements el = {.p = 1.2e7 * (1.0 - e * e), .e = e, .i = 0.6, .raan = 1.1, .argp = 0.4,
                           .nu = 0.3};
  myb_pv pv = {};
  MT_CHECK(myb_pv_from_elements(mu, &el, &pv) == MYB_OK);
  return pv;
}

static double period(void) {
  return 2.0 * pi * sqrt(1.2e7 * 1.2e7 * 1.2e7 / mu);
}

static double norm3(const double* v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static double time_offset(const double* y, const myb_pv* exact) {
  const double d[3] = {y[0] - exact->r[0], y[1] - exact->r[1], y[2] - exact->r[2]};
  return norm3(d) / norm3(exact->v);
}

static myb_rkn12* make(double rtol, double atol, double h_fixed) {
  const myb_rkn12_options options = {
      .rtol = rtol, .atol = atol, .h_initial = h_fixed, .h_max = h_fixed};
  myb_rkn12* s = NULL;
  MT_CHECK(myb_rkn12_create(3, &options, NULL, &s) == MYB_OK);
  return s;
}

static double run(double e, double rtol, double h_fixed, double orbits, myb_rkn12_stats* stats) {
  const myb_pv start = orbit(e);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  myb_rkn12* s = make(rtol, h_fixed > 0.0 ? 1e10 : rtol, h_fixed);
  MT_CHECK(myb_rkn12_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  const double t_end = orbits * period();
  MT_CHECK(myb_rkn12_integrate(s, t_end) == MYB_OK);
  myb_pv exact = {};
  MT_CHECK(myb_kepler_propagate(mu, &start, t_end, &exact) == MYB_OK);
  if (stats != NULL) {
    *stats = myb_rkn12_statistics(s);
  }
  const double err = time_offset(myb_rkn12_state(s), &exact);
  myb_rkn12_destroy(s);
  return err;
}

static void test_order(void) {
  /* Fixed steps (tolerances that accept every step), one orbit of e = 0.5. */
  double previous = 0.0;
  double order = 0.0;
  for (int k = 0; k < 3; ++k) {
    const int steps = 10 << k; /* the error stays above round-off */
    const double err = run(0.5, 1e10, period() / steps, 1.0, NULL);
    if (k > 0) {
      order = log2(previous / err);
      printf("%3d steps per orbit: time offset %.3e s, observed order %.2f\n", steps, err, order);
    }
    previous = err;
  }
  MT_CHECK(order > 11.0 && order < 13.5);
}

static void test_adaptive(void) {
  double previous = 1e300;
  for (int k = 8; k <= 14; k += 2) {
    myb_rkn12_stats st = {};
    const double err = run(0.9, pow(10.0, -k), 0.0, 10.0, &st);
    printf("rtol 1e-%d: 10 orbits e = 0.9, time offset %.3e s (%.1e of the period), %lld steps "
           "(%lld rejected)\n",
           k, err, err / period(), (long long)st.accepted, (long long)st.rejected);
    MT_CHECK(err < previous);
    previous = err;
  }
  MT_CHECK(previous / period() < 1e-12);
}

static double radial_velocity(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[0] * y[3] + y[1] * y[4] + y[2] * y[5];
}

static void test_dense_and_events(void) {
  const myb_pv start = orbit(0.6);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  myb_rkn12* s = make(1e-13, 1e-13, 0.0);
  MT_CHECK(myb_rkn12_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  const double t_end = 10.0 * period();
  double worst_dense = 0.0;
  double worst_event = 0.0;
  int events = 0;
  double t_previous = 0.0;
  /* bound: the integration ends at t_end; at most 10^6 steps */
  for (int i = 0; i < 1000000 && myb_rkn12_time(s) < t_end; ++i) {
    MT_CHECK(myb_rkn12_step(s, t_end) == MYB_OK);
    const double t_now = myb_rkn12_time(s);
    const double tm = t_previous + 0.37 * (t_now - t_previous);
    double y[6];
    MT_CHECK(myb_rkn12_dense(s, tm, y) == MYB_OK);
    myb_pv exact = {};
    MT_CHECK(myb_kepler_propagate(mu, &start, tm, &exact) == MYB_OK);
    worst_dense = fmax(worst_dense, time_offset(y, &exact));
    int found = 0;
    double t_event = 0.0;
    MT_CHECK(myb_rkn12_locate(s, radial_velocity, NULL, &found, &t_event, NULL) == MYB_OK);
    if (found) {
      MT_CHECK(myb_kepler_propagate(mu, &start, t_event, &exact) == MYB_OK);
      const double rv = exact.r[0] * exact.v[0] + exact.r[1] * exact.v[1] + exact.r[2] * exact.v[2];
      const double v2 = exact.v[0] * exact.v[0] + exact.v[1] * exact.v[1] + exact.v[2] * exact.v[2];
      worst_event = fmax(worst_event, fabs(rv / (v2 - mu / norm3(exact.r))));
      ++events;
    }
    t_previous = t_now;
  }
  printf("dense output: time offset %.3e s; %d apsis events, worst timing error %.3e s\n",
         worst_dense, events, worst_event);
  MT_CHECK(events == 20);
  MT_CHECK(worst_dense < 1e-8);
  MT_CHECK(worst_event < 1e-7);
  double y[6];
  MT_CHECK(myb_rkn12_dense(s, t_end + 1.0, y) == MYB_E_RANGE);
  myb_rkn12_destroy(s);
}

static void test_backward_and_stops(void) {
  const myb_pv start = orbit(0.4);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  myb_rkn12* s = make(1e-14, 1e-14, 0.0);
  MT_CHECK(myb_rkn12_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_rkn12_integrate(s, 3.0 * period()) == MYB_OK);
  const long long forward = (long long)myb_rkn12_statistics(s).accepted;
  MT_CHECK(myb_rkn12_integrate(s, 0.0) == MYB_OK);
  MT_CHECK(myb_rkn12_time(s) == 0.0);
  const double* y = myb_rkn12_state(s);
  const double d[3] = {y[0] - y0[0], y[1] - y0[1], y[2] - y0[2]};
  printf("3 orbits forward and back: relative position difference %.3e (%lld steps forward)\n",
         norm3(d) / norm3(y0), forward);
  MT_CHECK(norm3(d) / norm3(y0) < 1e-12);
  myb_rkn12_destroy(s);
  s = make(1e-14, 1e-14, 0.0);
  MT_CHECK(myb_rkn12_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  for (int k = 1; k <= 30; ++k) {
    MT_CHECK(myb_rkn12_integrate(s, 3.0 * period() * k / 30.0) == MYB_OK);
  }
  const long long with_stops = (long long)myb_rkn12_statistics(s).accepted;
  printf("with 30 stops: %lld steps\n", with_stops);
  MT_CHECK(with_stops <= forward + 30);
  myb_rkn12_destroy(s);
}

static myb_status failing(void* context, double t, const double* x, double* a) {
  (void)context;
  (void)x;
  a[0] = 1.0;
  return t > 0.5 ? MYB_E_DOMAIN : MYB_OK;
}

static myb_status singular(void* context, double t, const double* x, double* a) {
  (void)context;
  (void)x;
  a[0] = 2.0 / ((1.0 - t) * (1.0 - t) * (1.0 - t)); /* x = 1/(1-t) */
  return MYB_OK;
}

static void test_errors(void) {
  myb_rkn12* s = NULL;
  const myb_rkn12_options bad = {.rtol = 0.0, .atol = 1e-9};
  MT_CHECK(myb_rkn12_create(1, &bad, NULL, &s) == MYB_E_DOMAIN);
  const myb_rkn12_options good = {.rtol = 1e-10, .atol = 1e-10, .h_min = 1e-12};
  MT_CHECK(myb_rkn12_create(1, &good, NULL, &s) == MYB_OK);
  const double y0[2] = {1.0, 1.0};
  MT_CHECK(myb_rkn12_start(s, failing, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_rkn12_integrate(s, 1.0) == MYB_E_DOMAIN);
  MT_CHECK(myb_rkn12_start(s, singular, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_rkn12_integrate(s, 2.0) == MYB_E_STEP_MIN);
  MT_CHECK(myb_rkn12_time(s) < 1.0);
  myb_rkn12_destroy(s);
  myb_rkn12_destroy(NULL);
}

int main(void) {
  test_coefficients();
  test_order();
  test_adaptive();
  test_dense_and_events();
  test_backward_and_stops();
  test_errors();
  return mt_finish("test_rkn12");
}
