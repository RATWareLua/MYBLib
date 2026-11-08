/* SPDX-License-Identifier: BSD-3-Clause */
/* Gauss-Jackson: structure of the coefficients, observed order 8, accuracy against the
 * analytic two-body solution (100 orbits), dense output and events, a velocity-dependent force,
 * backward integration, failure of the start for a step that is too long, error statuses. */
#include <math.h>

#include "core/fpguard.h"
#include "integ/gj_coefficients.h"
#include "mtest.h"
#include "myb/integ.h"
#include "myb/kepler.h"

static const double pi = 3.141592653589793;
static const double mu = 3.986004418e14;

static void test_coefficients(void) {
  /* On a constant acceleration P1(D) = 1/2 and P2(D) = 1/12 at every window position; the
   * Lagrange basis sums to 1 (only the constant term survives). */
  for (int j = 0; j < 10; ++j) {
    double sv = 0.0;
    double sx = 0.0;
    for (int k = 0; k < 9; ++k) {
      sv += myb_gj_velocity[j][k];
      sx += myb_gj_position[j][k];
    }
    MT_CHECK(fabs(sv - 0.5) < 1e-13);
    MT_CHECK(fabs(sx - 1.0 / 12.0) < 1e-13);
  }
  for (int m = 0; m < 9; ++m) {
    double sum = 0.0;
    for (int k = 0; k < 9; ++k) {
      sum += myb_gj_lagrange[k][m];
    }
    MT_CHECK(fabs(sum - (m == 0 ? 1.0 : 0.0)) < 1e-15);
  }
}

static myb_status two_body(void* context, double t, const double* x, const double* v, double* a) {
  (void)context;
  (void)t;
  (void)v;
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

static myb_gauss_jackson* make(size_t n, double h) {
  const myb_gauss_jackson_options options = {.h = h};
  myb_gauss_jackson* s = NULL;
  MT_CHECK(myb_gauss_jackson_create(n, &options, NULL, &s) == MYB_OK);
  return s;
}

static double run(double e, double h, double orbits, myb_gauss_jackson_stats* stats) {
  const myb_pv start = orbit(e);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  myb_gauss_jackson* s = make(3, h);
  MT_CHECK(myb_gauss_jackson_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  const double t_end = orbits * period();
  MT_CHECK(myb_gauss_jackson_integrate(s, t_end) == MYB_OK);
  double y[6];
  MT_CHECK(myb_gauss_jackson_dense(s, t_end, y) == MYB_OK);
  myb_pv exact = {};
  MT_CHECK(myb_kepler_propagate(mu, &start, t_end, &exact) == MYB_OK);
  if (stats != NULL) {
    *stats = myb_gauss_jackson_statistics(s);
  }
  myb_gauss_jackson_destroy(s);
  return time_offset(y, &exact);
}

static void test_order(void) {
  /* 10 orbits of e = 0.1. The window interpolates the acceleration to degree 8, so positions
   * are exact to degree 10: the observed order of the global error approaches 10. */
  double previous = 0.0;
  double order = 0.0;
  for (int k = 0; k < 3; ++k) {
    const int steps = 40 << k;
    const double err = run(0.1, period() / steps, 10.0, NULL);
    if (k > 0) {
      order = log2(previous / err);
      printf("%4d steps per orbit: time offset %.3e s, observed order %.2f\n", steps, err, order);
    }
    previous = err;
  }
  MT_CHECK(order > 9.5 && order < 10.8);
}

static void test_accuracy(void) {
  myb_gauss_jackson_stats st = {};
  const double err = run(0.1, period() / 200.0, 100.0, &st);
  printf("100 orbits e = 0.1, 200 steps per orbit: time offset %.3e s (%.1e of the period); "
         "%.2f evaluations per step, predictor-corrector difference %.1e, unconverged %lld\n",
         err, err / period(), (double)st.evaluations / (double)st.steps, st.max_difference,
         (long long)st.unconverged);
  MT_CHECK(err / period() < 1e-12);
  MT_CHECK(st.unconverged == 0);
  MT_CHECK(st.max_difference > 0.0 && st.max_difference < 1e-10);
}

static double radial_velocity(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[0] * y[3] + y[1] * y[4] + y[2] * y[5];
}

static void test_dense_and_events(void) {
  const myb_pv start = orbit(0.3);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  const double h = period() / 300.0;
  myb_gauss_jackson* s = make(3, h);
  MT_CHECK(myb_gauss_jackson_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  double worst_dense = 0.0;
  double worst_node = 0.0;
  double worst_event = 0.0;
  int events = 0;
  for (int i = 0; i < 3000; ++i) { /* 10 orbits */
    MT_CHECK(myb_gauss_jackson_step(s) == MYB_OK);
    const double t = myb_gauss_jackson_time(s);
    double y[6];
    /* Between nodes, against the analytic solution. */
    const double tm = t - 0.37 * h;
    MT_CHECK(myb_gauss_jackson_dense(s, tm, y) == MYB_OK);
    myb_pv exact = {};
    MT_CHECK(myb_kepler_propagate(mu, &start, tm, &exact) == MYB_OK);
    worst_dense = fmax(worst_dense, time_offset(y, &exact));
    /* At the node, against the state (same polynomial, different rounding). */
    MT_CHECK(myb_gauss_jackson_dense(s, t, y) == MYB_OK);
    const double* state = myb_gauss_jackson_state(s);
    const double d[3] = {y[0] - state[0], y[1] - state[1], y[2] - state[2]};
    worst_node = fmax(worst_node, norm3(d) / norm3(state));
    int found = 0;
    double t_event = 0.0;
    MT_CHECK(myb_gauss_jackson_locate(s, radial_velocity, NULL, &found, &t_event, NULL) ==
             MYB_OK);
    if (found) {
      MT_CHECK(myb_kepler_propagate(mu, &start, t_event, &exact) == MYB_OK);
      const double rv = exact.r[0] * exact.v[0] + exact.r[1] * exact.v[1] + exact.r[2] * exact.v[2];
      const double v2 = exact.v[0] * exact.v[0] + exact.v[1] * exact.v[1] + exact.v[2] * exact.v[2];
      worst_event = fmax(worst_event, fabs(rv / (v2 - mu / norm3(exact.r))));
      ++events;
    }
  }
  printf("dense output: time offset %.3e s between nodes, %.1e relative at nodes; %d apsis "
         "events, worst timing error %.3e s\n",
         worst_dense, worst_node, events, worst_event);
  MT_CHECK(events == 20);
  MT_CHECK(worst_dense < 1e-8);
  MT_CHECK(worst_node < 5e-14);
  MT_CHECK(worst_event < 1e-7);
  double y[6];
  MT_CHECK(myb_gauss_jackson_dense(s, myb_gauss_jackson_time(s) + h, y) == MYB_E_RANGE);
  myb_gauss_jackson_destroy(s);
}

/* Charged particle in uniform fields (see test_ias15.c): a velocity-dependent force. */
static const double gyro = 2.0;
static const double field[3] = {0.3, -0.1, 0.05};

static myb_status lorentz(void* context, double t, const double* x, const double* v, double* a) {
  (void)context;
  (void)t;
  (void)x;
  a[0] = gyro * v[1] + field[0];
  a[1] = -gyro * v[0] + field[1];
  a[2] = field[2];
  return MYB_OK;
}

static void test_lorentz(void) {
  const double y0[6] = {1.0, 0.0, 0.0, 0.2, 0.5, 0.1};
  myb_gauss_jackson* s = make(3, 0.01);
  MT_CHECK(myb_gauss_jackson_start(s, lorentz, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_gauss_jackson_integrate(s, 20.0) == MYB_OK);
  double y[6];
  MT_CHECK(myb_gauss_jackson_dense(s, 20.0, y) == MYB_OK);
  const double t = 20.0;
  const double vd[2] = {field[1] / gyro, -field[0] / gyro};
  const double u0[2] = {y0[3] - vd[0], y0[4] - vd[1]};
  const double c = cos(gyro * t);
  const double sn = sin(gyro * t);
  const double exact[3] = {y0[0] + vd[0] * t + (u0[0] * sn + u0[1] * (1.0 - c)) / gyro,
                           y0[1] + vd[1] * t + (-u0[0] * (1.0 - c) + u0[1] * sn) / gyro,
                           y0[2] + y0[5] * t + 0.5 * field[2] * t * t};
  const double d[3] = {y[0] - exact[0], y[1] - exact[1], y[2] - exact[2]};
  printf("charged particle (v-dependent force), h = 0.01: relative error %.3e\n",
         norm3(d) / norm3(exact));
  MT_CHECK(norm3(d) / norm3(exact) < 1e-12);
  myb_gauss_jackson_destroy(s);
}

static void test_backward(void) {
  const myb_pv start = orbit(0.2);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  const double h = -period() / 200.0;
  myb_gauss_jackson* s = make(3, h);
  MT_CHECK(myb_gauss_jackson_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  const double t_end = -3.0 * period();
  MT_CHECK(myb_gauss_jackson_integrate(s, t_end) == MYB_OK);
  MT_CHECK(myb_gauss_jackson_time(s) <= t_end);
  double y[6];
  MT_CHECK(myb_gauss_jackson_dense(s, t_end, y) == MYB_OK);
  myb_pv exact = {};
  MT_CHECK(myb_kepler_propagate(mu, &start, t_end, &exact) == MYB_OK);
  printf("3 orbits backward: time offset %.3e s\n", time_offset(y, &exact));
  MT_CHECK(time_offset(y, &exact) < 1e-8);
  MT_CHECK(myb_gauss_jackson_integrate(s, 0.0) == MYB_E_DOMAIN); /* the grid runs one way */
  myb_gauss_jackson_destroy(s);
}

static void test_errors(void) {
  myb_gauss_jackson* s = NULL;
  const myb_gauss_jackson_options zero = {.h = 0.0};
  MT_CHECK(myb_gauss_jackson_create(3, &zero, NULL, &s) == MYB_E_DOMAIN);
  /* A step of a third of the period: either the start iteration diverges, or the steps report
   * a predictor-corrector difference of the size of the orbit (the method has no error
   * control: this indicator is the signal to shorten h). */
  s = make(3, period() / 3.0);
  const myb_pv start = orbit(0.5);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  const myb_status status = myb_gauss_jackson_start(s, two_body, NULL, 0.0, y0);
  MT_CHECK(status == MYB_OK || status == MYB_E_NOCONV);
  if (status == MYB_OK) {
    myb_status stepped = MYB_OK;
    for (int i = 0; i < 6 && stepped == MYB_OK; ++i) {
      stepped = myb_gauss_jackson_step(s);
    }
    const double difference = myb_gauss_jackson_statistics(s).max_difference;
    printf("step of P/3: predictor-corrector difference %.2e%c", difference, 10);
    MT_CHECK(stepped != MYB_OK || difference > 1e-3);
  }
  myb_gauss_jackson_destroy(s);
  myb_gauss_jackson_destroy(NULL);
}

int main(void) {
  test_coefficients();
  test_order();
  test_accuracy();
  test_dense_and_events();
  test_lorentz();
  test_backward();
  test_errors();
  return mt_finish("test_gauss_jackson");
}
