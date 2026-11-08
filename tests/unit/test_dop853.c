/* SPDX-License-Identifier: BSD-3-Clause */
/* DOP853: Butcher conditions of the coefficients, observed order of convergence, adaptive
 * accuracy against the analytic two-body solution (including 10^4 revolutions), dense output,
 * event location, the state transition matrix from variational equations, error statuses. */
#include <math.h>
#include <string.h>

#include "core/fpguard.h"
#include "integ/dop853_coefficients.h"
#include "mtest.h"
#include "myb/integ.h"
#include "myb/kepler.h"

static const double mu = 3.986004418e14;

/* ---- Coefficients -------------------------------------------------------------------------- */

static double b_dot(const double* v) {
  double s = 0.0;
  for (int i = 0; i < MYB_DOP853_STAGES; ++i) {
    s += myb_dop853_a[12][i] * v[i];
  }
  return s;
}

/* (A v)_i over the method stages. */
static void a_times(const double* v, double* out) {
  for (int i = 0; i < MYB_DOP853_STAGES; ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) {
      s += myb_dop853_a[i][j] * v[j];
    }
    out[i] = s;
  }
}

static void test_coefficients(void) {
  /* Row sums: sum_j a_ij = c_i for every stage, extra stages included. */
  double worst = 0.0;
  for (int i = 1; i < MYB_DOP853_STAGES_EXTENDED; ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) {
      s += myb_dop853_a[i][j];
    }
    worst = fmax(worst, fabs(s - myb_dop853_c[i]));
  }
  MT_CHECK(worst < 1e-14);
  /* Quadrature conditions sum b c^(k-1) = 1/k, k = 1..8 (order 8). */
  double c[MYB_DOP853_STAGES];
  double ck[MYB_DOP853_STAGES];
  for (int i = 0; i < MYB_DOP853_STAGES; ++i) {
    c[i] = myb_dop853_c[i];
    ck[i] = 1.0;
  }
  for (int k = 1; k <= 8; ++k) {
    MT_CHECK(fabs(b_dot(ck) - 1.0 / k) < 1e-14);
    for (int i = 0; i < MYB_DOP853_STAGES; ++i) {
      ck[i] *= c[i];
    }
  }
  /* All trees of orders 3, 4 and 5 (they involve every entry of A). */
  double ac[12], ac2[12], aac[12], ac3[12], aac2[12], aaac[12], cac[12], aca[12];
  double c2[12], c3[12], v[12];
  for (int i = 0; i < 12; ++i) {
    c2[i] = c[i] * c[i];
    c3[i] = c2[i] * c[i];
  }
  a_times(c, ac);
  a_times(c2, ac2);
  a_times(c3, ac3);
  a_times(ac, aac);
  a_times(ac2, aac2);
  a_times(aac, aaac);
  for (int i = 0; i < 12; ++i) {
    cac[i] = c[i] * ac[i];
  }
  a_times(cac, aca);
  struct {
    double value;
    double expected;
  } trees[] = {
      {b_dot(ac), 1.0 / 6.0},
      {0, 1.0 / 8.0}, {0, 1.0 / 12.0}, {b_dot(aac), 1.0 / 24.0},
      {0, 1.0 / 10.0}, {0, 1.0 / 15.0}, {0, 1.0 / 30.0}, {0, 1.0 / 20.0},
      {b_dot(ac3), 1.0 / 20.0}, {b_dot(aca), 1.0 / 40.0}, {b_dot(aac2), 1.0 / 60.0},
      {b_dot(aaac), 1.0 / 120.0},
  };
  for (int i = 0; i < 12; ++i) v[i] = c[i] * ac[i];
  trees[1].value = b_dot(v); /* b c (A c) */
  trees[2].value = b_dot(ac2);
  for (int i = 0; i < 12; ++i) v[i] = c2[i] * ac[i];
  trees[4].value = b_dot(v); /* b c^2 (A c) */
  for (int i = 0; i < 12; ++i) v[i] = c[i] * ac2[i];
  trees[5].value = b_dot(v); /* b c (A c^2) */
  for (int i = 0; i < 12; ++i) v[i] = c[i] * aac[i];
  trees[6].value = b_dot(v); /* b c (A A c) */
  for (int i = 0; i < 12; ++i) v[i] = ac[i] * ac[i];
  trees[7].value = b_dot(v); /* b (A c)^2 */
  for (size_t i = 0; i < sizeof trees / sizeof trees[0]; ++i) {
    MT_CHECK(fabs(trees[i].value - trees[i].expected) < 1e-14);
  }
  /* The error estimators compare two consistent methods: their weights sum to 0. */
  double e5 = 0.0;
  for (int i = 0; i < 13; ++i) {
    e5 += myb_dop853_e5[i];
  }
  MT_CHECK(fabs(e5) < 1e-14);
}

/* ---- Two-body problem ----------------------------------------------------------------------- */

static myb_status two_body(void* context, double t, const double* y, double* dydt) {
  (void)context;
  (void)t;
  const double r = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
  const double k = -mu / (r * r * r);
  for (int i = 0; i < 3; ++i) {
    dydt[i] = y[3 + i];
    dydt[3 + i] = k * y[i];
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

static double period(double e) {
  const double a = 1.2e7;
  (void)e;
  return 2.0 * 3.141592653589793 * sqrt(a * a * a / mu);
}

static double position_error(const double* y, const myb_pv* exact) {
  const double d[3] = {y[0] - exact->r[0], y[1] - exact->r[1], y[2] - exact->r[2]};
  return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) /
         sqrt(exact->r[0] * exact->r[0] + exact->r[1] * exact->r[1] + exact->r[2] * exact->r[2]);
}

static myb_dop853* make(size_t n, double rtol, double atol, double h_fixed) {
  const myb_dop853_options options = {
      .rtol = rtol, .atol = atol, .h_initial = h_fixed, .h_max = h_fixed};
  myb_dop853* s = NULL;
  MT_CHECK(myb_dop853_create(n, &options, NULL, &s) == MYB_OK);
  return s;
}

static double run_two_body(double e, double rtol, double h_fixed, double revolutions,
                           myb_dop853_stats* stats) {
  const myb_pv start = orbit(e);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  myb_dop853* s = make(6, rtol, h_fixed > 0.0 ? 1e10 : rtol, h_fixed);
  MT_CHECK(myb_dop853_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  const double t_end = revolutions * period(e);
  MT_CHECK(myb_dop853_integrate(s, t_end) == MYB_OK);
  myb_pv exact = {};
  MT_CHECK(myb_kepler_propagate(mu, &start, t_end, &exact) == MYB_OK);
  const double err = position_error(myb_dop853_state(s), &exact);
  if (stats != NULL) {
    *stats = myb_dop853_statistics(s);
  }
  myb_dop853_destroy(s);
  return err;
}

static void test_order(void) {
  /* Fixed steps (tolerances that accept every step): error ~ h^8. */
  /* The observed order approaches 8 from below as the step enters the asymptotic regime; the
   * step counts stay above the round-off floor (~1e-13 near 400 steps per orbit). */
  double previous = 0.0;
  double previous_order = 0.0;
  for (int k = 0; k < 4; ++k) {
    const int steps = 25 << k;
    const double err = run_two_body(0.5, 1e10, period(0.5) / steps, 1.0, NULL);
    if (k > 0) {
      const double order = log2(previous / err);
      printf("fixed step: %4d steps per orbit, error %.3e, observed order %.2f\n", steps, err, order);
      MT_CHECK(order > 7.2 && order < 8.6 && order > previous_order);
      previous_order = order;
    }
    previous = err;
  }
  MT_CHECK(previous_order > 7.7);
}

static void test_adaptive(void) {
  double previous = 1.0;
  for (int k = 6; k <= 13; k += 1) {
    const double tol = pow(10.0, -k);
    myb_dop853_stats st;
    const double err = run_two_body(0.7, tol, 0.0, 10.0, &st);
    printf("adaptive rtol 1e-%02d: 10 orbits e = 0.7, error %.3e, %lld steps (%lld rejected)\n", k,
           err, (long long)st.accepted, (long long)st.rejected);
    MT_CHECK(err < previous);
    previous = err;
  }
  myb_dop853_stats st;
  const double err = run_two_body(0.3, 1e-13, 0.0, 10000.0, &st);
  /* A non-symplectic method lets the energy drift linearly, so the phase error grows
   * quadratically with time; long-term work belongs to IAS15 or symplectic methods. */
  printf("10^4 orbits e = 0.3 at rtol 1e-13: relative position error %.3e, %lld steps\n", err,
         (long long)st.accepted);
  MT_CHECK(err < 1e-4);
}

static double radial_velocity(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[0] * y[3] + y[1] * y[4] + y[2] * y[5];
}

static void test_dense_and_events(void) {
  const myb_pv start = orbit(0.6);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  myb_dop853* s = make(6, 1e-12, 1e-12, 0.0);
  MT_CHECK(myb_dop853_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  const double p = period(0.6);
  double worst_dense = 0.0;
  double worst_event = 0.0;
  int events = 0;
  const double t_end = 20.0 * p;
  double t_previous = 0.0;
  myb_pv step_start = start;
  /* bound: the integration ends at t_end; at most 10^6 steps */
  for (int i = 0; i < 1000000 && myb_dop853_time(s) < t_end; ++i) {
    MT_CHECK(myb_dop853_step(s, t_end) == MYB_OK);
    const double t_now = myb_dop853_time(s);
    /* Interpolation error alone: compared with the analytic propagation from the state at the
     * start of the step (the global error of the integration is not part of it). */
    const double tm = t_previous + 0.37 * (t_now - t_previous);
    double y[6];
    MT_CHECK(myb_dop853_dense(s, tm, y) == MYB_OK);
    myb_pv exact = {};
    MT_CHECK(myb_kepler_propagate(mu, &step_start, tm - t_previous, &exact) == MYB_OK);
    worst_dense = fmax(worst_dense, position_error(y, &exact));
    /* Events: radial velocity changes sign at periapsis and apoapsis. On the local analytic
     * trajectory the located time must give r.v = 0; the timing error is (r.v) / d(r.v)/dt,
     * d(r.v)/dt = v^2 - mu / r. */
    int found = 0;
    double t_event = 0.0;
    MT_CHECK(myb_dop853_locate(s, radial_velocity, NULL, &found, &t_event, y) == MYB_OK);
    if (found) {
      myb_pv at = {};
      MT_CHECK(myb_kepler_propagate(mu, &step_start, t_event - t_previous, &at) == MYB_OK);
      const double rv = at.r[0] * at.v[0] + at.r[1] * at.v[1] + at.r[2] * at.v[2];
      const double v2 = at.v[0] * at.v[0] + at.v[1] * at.v[1] + at.v[2] * at.v[2];
      const double r = sqrt(at.r[0] * at.r[0] + at.r[1] * at.r[1] + at.r[2] * at.r[2]);
      worst_event = fmax(worst_event, fabs(rv / (v2 - mu / r)));
      ++events;
    }
    const double* y_now = myb_dop853_state(s);
    step_start = (myb_pv){.r = {y_now[0], y_now[1], y_now[2]}, .v = {y_now[3], y_now[4], y_now[5]}};
    t_previous = t_now;
  }
  printf("dense output: worst relative position error %.3e; %d apsis events, worst timing "
         "error %.3e s (period %.0f s)\n",
         worst_dense, events, worst_event, p);
  MT_CHECK(events == 40);
  MT_CHECK(worst_dense < 1e-11);
  MT_CHECK(worst_event < 1e-6);
  myb_dop853_destroy(s);
}

/* ---- Variational equations: state transition matrix of the two-body problem ------------------ */

/* y = (r, v, Phi[36] row-major); Phi' = [[0, I], [G, 0]] Phi, G = mu / r^5 (3 r r^T - r^2 I). */
static myb_status two_body_stm(void* context, double t, const double* y, double* dydt) {
  (void)context;
  const myb_status status = two_body(NULL, t, y, dydt);
  const double r2 = y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
  const double r = sqrt(r2);
  const double k = mu / (r2 * r2 * r);
  double g[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      g[i][j] = k * (3.0 * y[i] * y[j] - (i == j ? r2 : 0.0));
    }
  }
  const double* phi = y + 6;
  double* d = dydt + 6;
  for (int c = 0; c < 6; ++c) {
    for (int i = 0; i < 3; ++i) {
      d[i * 6 + c] = phi[(3 + i) * 6 + c];
      double sum = 0.0;
      for (int j = 0; j < 3; ++j) {
        sum += g[i][j] * phi[j * 6 + c];
      }
      d[(3 + i) * 6 + c] = sum;
    }
  }
  return status;
}

static void test_stm(void) {
  const myb_pv start = orbit(0.4);
  double y0[42] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  for (int i = 0; i < 6; ++i) {
    y0[6 + i * 6 + i] = 1.0;
  }
  myb_dop853* s = make(42, 1e-13, 1e-13, 0.0);
  MT_CHECK(myb_dop853_start(s, two_body_stm, NULL, 0.0, y0) == MYB_OK);
  const double t = 0.37 * period(0.4);
  MT_CHECK(myb_dop853_integrate(s, t) == MYB_OK);
  const double* phi = myb_dop853_state(s) + 6;
  /* Against the closed-form matrix (verified against MPFR in test_mpfr_kepler), in the scales
   * of the state: positions in |r0|, velocities in |v0|. */
  myb_pv end = {};
  double exact[36];
  MT_CHECK(myb_kepler_stm(mu, &start, t, &end, exact) == MYB_OK);
  const double rs = sqrt(start.r[0] * start.r[0] + start.r[1] * start.r[1] + start.r[2] * start.r[2]);
  const double vs = sqrt(start.v[0] * start.v[0] + start.v[1] * start.v[1] + start.v[2] * start.v[2]);
  double worst = 0.0;
  for (int i = 0; i < 6; ++i) {
    for (int c = 0; c < 6; ++c) {
      const double scale = (c < 3 ? rs : vs) / (i < 3 ? rs : vs);
      worst = fmax(worst, fabs(phi[i * 6 + c] - exact[i * 6 + c]) * scale);
    }
  }
  printf("state transition matrix (variational equations, rtol 1e-13) vs closed form: worst "
         "scaled difference %.3e\n",
         worst);
  MT_CHECK(worst < 1e-10);
  myb_dop853_destroy(s);
}

/* ---- Backward integration and intermediate stops ------------------------------------------------- */

static void test_backward_and_stops(void) {
  const double p = period(0.4);
  const myb_pv start = orbit(0.4);
  const double y0[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
  myb_dop853* s = make(6, 1e-13, 1e-13, 0.0);
  MT_CHECK(myb_dop853_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_dop853_integrate(s, 3.0 * p) == MYB_OK);
  const long long forward = myb_dop853_statistics(s).accepted;
  MT_CHECK(myb_dop853_integrate(s, 0.0) == MYB_OK);
  MT_CHECK(myb_dop853_time(s) == 0.0);
  const double back = position_error(myb_dop853_state(s), &start);
  printf("3 orbits forward and back at rtol 1e-13: relative position difference %.3e\n", back);
  MT_CHECK(back < 1e-10);
  myb_dop853_destroy(s);
  /* A step cut at a stop must not shrink the following ones. */
  s = make(6, 1e-13, 1e-13, 0.0);
  MT_CHECK(myb_dop853_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  for (int k = 1; k <= 30; ++k) {
    MT_CHECK(myb_dop853_integrate(s, 3.0 * p * k / 30.0) == MYB_OK);
  }
  const long long with_stops = myb_dop853_statistics(s).accepted;
  printf("with 30 stops: %lld steps (%lld without)\n", with_stops, forward);
  MT_CHECK(with_stops <= forward + 30);
  myb_dop853_destroy(s);
}

/* ---- Errors ------------------------------------------------------------------------------------ */

static myb_status failing(void* context, double t, const double* y, double* dydt) {
  (void)context;
  (void)y;
  dydt[0] = 1.0;
  return t > 0.5 ? MYB_E_DOMAIN : MYB_OK;
}

static myb_status singular(void* context, double t, const double* y, double* dydt) {
  (void)context;
  (void)y;
  dydt[0] = 1.0 / ((1.0 - t) * (1.0 - t)); /* y = 1/(1-t) blows up at t = 1 */
  return MYB_OK;
}

static void test_errors(void) {
  myb_dop853* s = NULL;
  const myb_dop853_options bad = {.rtol = 0.0, .atol = 1e-9};
  MT_CHECK(myb_dop853_create(1, &bad, NULL, &s) == MYB_E_DOMAIN);
  const myb_dop853_options good = {.rtol = 1e-9, .atol = 1e-9, .h_min = 1e-12};
  MT_CHECK(myb_dop853_create(1, &good, NULL, &s) == MYB_OK);
  const double y0[1] = {1.0};
  MT_CHECK(myb_dop853_start(s, failing, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_dop853_integrate(s, 1.0) == MYB_E_DOMAIN); /* the callback's status */
  MT_CHECK(myb_dop853_start(s, singular, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_dop853_integrate(s, 2.0) == MYB_E_STEP_MIN);
  MT_CHECK(myb_dop853_time(s) < 1.0);
  double y[1];
  MT_CHECK(myb_dop853_dense(s, 5.0, y) == MYB_E_RANGE);
  myb_dop853_destroy(s);
}

int main(void) {
  test_coefficients();
  test_order();
  test_adaptive();
  test_dense_and_events();
  test_stm();
  test_backward_and_stops();
  test_errors();
  return mt_finish("test_dop853");
}
