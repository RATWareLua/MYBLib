/* SPDX-License-Identifier: BSD-3-Clause */
/* IAS15: consistency of the Gauss-Radau constants, observed order of convergence, accuracy
 * against the analytic two-body solution (ellipses up to e = 0.99, hyperbola), Brouwer's law of
 * the energy error, dense output, event location, a velocity-dependent force, blocks of the
 * timescale, backward integration, intermediate stops, error statuses. */
#include <math.h>
#include <string.h>

#include "core/fpguard.h"
#include "integ/ias15_coefficients.h"
#include "mtest.h"
#include "myb/integ.h"
#include "myb/kepler.h"

static const double pi = 3.141592653589793;

/* ---- Constants -------------------------------------------------------------------------------- */

static double legendre_sum(double h) { /* P_7(x) + P_8(x), x = 2h - 1 */
  const double x = 2.0 * h - 1.0;
  double p0 = 1.0;
  double p1 = x;
  double p7 = 0.0;
  for (int k = 1; k < 8; ++k) {
    const double p2 = ((2.0 * k + 1.0) * x * p1 - k * p0) / (k + 1.0);
    p0 = p1;
    p1 = p2;
    if (k == 6) {
      p7 = p1;
    }
  }
  return p7 + p1;
}

static void test_constants(void) {
  MT_CHECK(myb_ias15_h[0] == 0.0);
  for (int i = 1; i < 8; ++i) {
    MT_CHECK(myb_ias15_h[i] > myb_ias15_h[i - 1] && myb_ias15_h[i] < 1.0);
    MT_CHECK(fabs(legendre_sum(myb_ias15_h[i])) < 1e-13);
    for (int k = 0; k < i; ++k) {
      MT_CHECK(fabs(myb_ias15_r[i][k] * (myb_ias15_h[i] - myb_ias15_h[k]) - 1.0) < 1e-15);
    }
  }
  /* C D = I */
  double worst = 0.0;
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      double sum = 0.0;
      for (int k = 0; k < 7; ++k) {
        sum += myb_ias15_c[i][k] * myb_ias15_d[k][j];
      }
      worst = fmax(worst, fabs(sum - (i == j ? 1.0 : 0.0)));
    }
  }
  MT_CHECK(worst < 1e-14);
  /* Column j of C holds h (h - h_1) ... (h - h_j): it vanishes at the nodes h_1..h_j. */
  for (int j = 0; j < 7; ++j) {
    for (int node = 1; node <= j; ++node) {
      const double h = myb_ias15_h[node];
      double value = 0.0;
      double power = h;
      double scale = 0.0;
      for (int k = 0; k <= j; ++k) {
        value += myb_ias15_c[k][j] * power;
        scale += fabs(myb_ias15_c[k][j] * power);
        power *= h;
      }
      MT_CHECK(fabs(value) <= 1e-15 * scale);
    }
  }
}

/* ---- Harmonic oscillator: order of convergence ------------------------------------------------ */

static myb_status oscillator(void* context, double t, const double* x, const double* v,
                             double* a) {
  (void)context;
  (void)t;
  (void)v;
  a[0] = -x[0];
  return MYB_OK;
}

static myb_ias15* make(size_t n, double epsilon, double h_fixed, size_t block) {
  const myb_ias15_options options = {
      .epsilon = epsilon, .h_initial = h_fixed, .h_max = h_fixed, .block = block};
  myb_ias15* s = NULL;
  MT_CHECK(myb_ias15_create(n, &options, NULL, &s) == MYB_OK);
  return s;
}

static void test_order(void) {
  /* One step from x = 1, v = 0: the local error is O(dt^16). Steps of 1.5 to 2.4 radians keep
   * the error above round-off (at the larger ones the iteration is not fully converged). */
  const double steps[] = {2.4, 1.92, 1.536};
  double previous = 0.0;
  for (int k = 0; k < 3; ++k) {
    myb_ias15* s = make(1, 1e30, steps[k], 0);
    const double y0[2] = {1.0, 0.0};
    MT_CHECK(myb_ias15_start(s, oscillator, NULL, 0.0, y0) == MYB_OK);
    MT_CHECK(myb_ias15_step(s, 100.0) == MYB_OK);
    MT_CHECK(myb_ias15_time(s) == steps[k]);
    const double* y = myb_ias15_state(s);
    const double err = hypot(y[0] - cos(steps[k]), y[1] + sin(steps[k]));
    if (k > 0) {
      const double order = log(previous / err) / log(steps[k - 1] / steps[k]);
      printf("one step of %.3f: local error %.3e, observed order %.2f\n", steps[k], err, order);
      MT_CHECK(order > 15.0 && order < 17.0);
    }
    previous = err;
    myb_ias15_destroy(s);
  }
}

/* ---- Two-body problem -------------------------------------------------------------------------- */

static const double mu = 3.986004418e14;

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

static myb_pv orbit(double e, double nu) {
  const double p = e < 1.0 ? 1.2e7 * (1.0 - e * e) : 1.2e7 * (e * e - 1.0);
  const myb_elements el = {.p = p, .e = e, .i = 0.6, .raan = 1.1, .argp = 0.4, .nu = nu};
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

/* Distance from the analytic state, as a time offset along the orbit: |dr| / |v|. It removes
 * the amplification of a phase error by the speed at periapsis. */
static double time_offset(const double* y, const myb_pv* exact) {
  const double d[3] = {y[0] - exact->r[0], y[1] - exact->r[1], y[2] - exact->r[2]};
  return norm3(d) / norm3(exact->v);
}

static void load(const myb_pv* pv, double* y) {
  for (int i = 0; i < 3; ++i) {
    y[i] = pv->r[i];
    y[3 + i] = pv->v[i];
  }
}

static void test_two_body(void) {
  const double p = period();
  const double eccentricities[] = {0.0, 0.5, 0.9, 0.99};
  for (int q = 0; q < 4; ++q) {
    const myb_pv start = orbit(eccentricities[q], 0.3);
    double y0[6];
    load(&start, y0);
    myb_ias15* s = make(3, 0.0, 0.0, 0);
    MT_CHECK(myb_ias15_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
    MT_CHECK(myb_ias15_integrate(s, 10.0 * p) == MYB_OK);
    myb_pv exact = {};
    MT_CHECK(myb_kepler_propagate(mu, &start, 10.0 * p, &exact) == MYB_OK);
    const double dt = time_offset(myb_ias15_state(s), &exact);
    const myb_ias15_stats st = myb_ias15_statistics(s);
    printf("e = %.2f, 10 orbits: time offset %.3e s (%.1e of the period), %lld steps, %lld "
           "rejected, %.2f sweeps per attempt\n",
           eccentricities[q], dt, dt / p, (long long)st.accepted, (long long)st.rejected,
           (double)st.iterations / (double)(st.accepted + st.rejected));
    MT_CHECK(dt / p < 1e-12);
    MT_CHECK(st.unconverged == 0);
    myb_ias15_destroy(s);
  }
  /* 10^4 orbits of e = 0.3 (criterion of phase 3): the error stays at the level of the
   * accumulated round-off (Brouwer), far below that of the non-symplectic DOP853. */
  {
    const myb_pv start = orbit(0.3, 0.3);
    double y0[6];
    load(&start, y0);
    myb_ias15* s = make(3, 0.0, 0.0, 0);
    MT_CHECK(myb_ias15_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
    MT_CHECK(myb_ias15_integrate(s, 1e4 * p) == MYB_OK);
    myb_pv exact = {};
    MT_CHECK(myb_kepler_propagate(mu, &start, 1e4 * p, &exact) == MYB_OK);
    const double dt = time_offset(myb_ias15_state(s), &exact);
    printf("e = 0.30, 10^4 orbits: time offset %.3e s (%.1e of the period), %lld steps\n", dt,
           dt / p, (long long)myb_ias15_statistics(s).accepted);
    MT_CHECK(dt / p < 1e-10);
    myb_ias15_destroy(s);
  }
  /* Hyperbolic flyby, e = 2: from 3 h before to 3 h after periapsis. */
  const myb_pv before = orbit(2.0, -2.0);
  double y0[6];
  load(&before, y0);
  myb_ias15* s = make(3, 0.0, 0.0, 0);
  MT_CHECK(myb_ias15_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(s, 6.0 * 3600.0) == MYB_OK);
  myb_pv exact = {};
  MT_CHECK(myb_kepler_propagate(mu, &before, 6.0 * 3600.0, &exact) == MYB_OK);
  const double dt = time_offset(myb_ias15_state(s), &exact);
  printf("hyperbola e = 2, 6 h: time offset %.3e s, %lld steps\n", dt,
         (long long)myb_ias15_statistics(s).accepted);
  MT_CHECK(dt < 1e-9);
  myb_ias15_destroy(s);
}

/* ---- Brouwer's law ------------------------------------------------------------------------------ */

static double energy(const double* y) {
  return 0.5 * (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]) - mu / norm3(y);
}

static void test_brouwer(void) {
  /* Ensemble of orbits (e = 0.5, different starting anomalies): with the truncation error
   * below round-off, the r.m.s. relative energy error grows like t^(1/2) (a random walk), not
   * like t. The exactly circular orbit with constant steps is excluded on purpose: there the
   * rounding errors of consecutive steps are identical and do not average out. */
  enum { RUNS = 12, POINTS = 3 };
  const double orbits[POINTS] = {10.0, 100.0, 1000.0};
  double sum2[POINTS] = {};
  const double p = period();
  for (int run = 0; run < RUNS; ++run) {
    const myb_pv start = orbit(0.5, 0.3 + 0.5 * run);
    double y0[6];
    load(&start, y0);
    const double e0 = energy(y0);
    myb_ias15* s = make(3, 0.0, 0.0, 0);
    MT_CHECK(myb_ias15_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
    for (int k = 0; k < POINTS; ++k) {
      MT_CHECK(myb_ias15_integrate(s, orbits[k] * p) == MYB_OK);
      const double d = (energy(myb_ias15_state(s)) - e0) / e0;
      sum2[k] += d * d;
    }
    myb_ias15_destroy(s);
  }
  const double rms_first = sqrt(sum2[0] / (double)RUNS);
  const double rms_last = sqrt(sum2[POINTS - 1] / (double)RUNS);
  const double slope = log(rms_last / rms_first) / log(orbits[POINTS - 1] / orbits[0]);
  printf("Brouwer: r.m.s. relative energy error %.3e (10 orbits), %.3e (100), %.3e (1000); "
         "slope %.2f\n",
         rms_first, sqrt(sum2[1] / (double)RUNS), rms_last, slope);
  MT_CHECK(slope > 0.25 && slope < 0.75);
  MT_CHECK(rms_last < 3e-14);
}

/* ---- Dense output and events -------------------------------------------------------------------- */

static double radial_velocity(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[0] * y[3] + y[1] * y[4] + y[2] * y[5];
}

static void test_dense_and_events(void) {
  const double p = period();
  const myb_pv start = orbit(0.6, 0.3);
  double y0[6];
  load(&start, y0);
  myb_ias15* s = make(3, 0.0, 0.0, 0);
  MT_CHECK(myb_ias15_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  double y[6];
  MT_CHECK(myb_ias15_dense(s, 0.0, y) == MYB_E_RANGE); /* no step yet */
  double worst_dense = 0.0;
  double worst_event = 0.0;
  double worst_end = 0.0;
  int events = 0;
  const double t_end = 20.0 * p;
  double t_previous = 0.0;
  myb_pv step_start = start;
  /* bound: the integration ends at t_end; at most 10^6 steps */
  for (int i = 0; i < 1000000 && myb_ias15_time(s) < t_end; ++i) {
    MT_CHECK(myb_ias15_step(s, t_end) == MYB_OK);
    const double t_now = myb_ias15_time(s);
    /* Interpolation alone: against the analytic propagation from the start of the step. */
    for (int k = 1; k < 4; ++k) {
      const double tm = t_previous + 0.25 * k * (t_now - t_previous);
      MT_CHECK(myb_ias15_dense(s, tm, y) == MYB_OK);
      myb_pv exact = {};
      MT_CHECK(myb_kepler_propagate(mu, &step_start, tm - t_previous, &exact) == MYB_OK);
      worst_dense = fmax(worst_dense, time_offset(y, &exact));
    }
    /* At the end of the step the dense output meets the state, up to the resolution of the
     * double t_now (the integrator keeps the time compensated): offset in units of ulp(t). */
    MT_CHECK(myb_ias15_dense(s, t_now, y) == MYB_OK);
    const double* state = myb_ias15_state(s);
    const double d_end[3] = {y[0] - state[0], y[1] - state[1], y[2] - state[2]};
    const double ulp_t = nextafter(t_now, (double)INFINITY) - t_now;
    worst_end = fmax(worst_end, norm3(d_end) / norm3(state + 3) / ulp_t);
    int found = 0;
    double t_event = 0.0;
    MT_CHECK(myb_ias15_locate(s, radial_velocity, NULL, &found, &t_event, y) == MYB_OK);
    if (found) {
      myb_pv at = {};
      MT_CHECK(myb_kepler_propagate(mu, &step_start, t_event - t_previous, &at) == MYB_OK);
      const double rv = at.r[0] * at.v[0] + at.r[1] * at.v[1] + at.r[2] * at.v[2];
      const double v2 = at.v[0] * at.v[0] + at.v[1] * at.v[1] + at.v[2] * at.v[2];
      worst_event = fmax(worst_event, fabs(rv / (v2 - mu / norm3(at.r))));
      ++events;
    }
    step_start = (myb_pv){.r = {state[0], state[1], state[2]}, .v = {state[3], state[4], state[5]}};
    t_previous = t_now;
  }
  printf("dense output: worst time offset %.3e s, at the step end %.2f ulp(t); %d apsis "
         "events, worst timing error %.3e s\n",
         worst_dense, worst_end, events, worst_event);
  MT_CHECK(events == 40);
  MT_CHECK(worst_dense < 1e-9);
  MT_CHECK(worst_end <= 2.0);
  MT_CHECK(worst_event < 1e-8);
  MT_CHECK(myb_ias15_dense(s, t_end + 1.0, y) == MYB_E_RANGE);
  myb_ias15_destroy(s);
}

/* ---- Velocity-dependent force: charged particle in uniform fields ------------------------------ */

/* a = W (v x z) + g: gyration with frequency W about B = B z and a uniform acceleration g (an
 * electric field). With the drift v_d = (g_y, -g_x, 0) / W, u = v_perp - v_d rotates by -W t and
 * z moves with constant acceleration g_z; |a|^2 = W^2 |u|^2 + g_z^2 never vanishes. */
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

static void lorentz_exact(const double* y0, double t, double* y) {
  const double vd[2] = {field[1] / gyro, -field[0] / gyro};
  const double u0[2] = {y0[3] - vd[0], y0[4] - vd[1]};
  const double c = cos(gyro * t);
  const double s = sin(gyro * t);
  y[3] = vd[0] + u0[0] * c + u0[1] * s;
  y[4] = vd[1] - u0[0] * s + u0[1] * c;
  y[5] = y0[5] + field[2] * t;
  y[0] = y0[0] + vd[0] * t + (u0[0] * s + u0[1] * (1.0 - c)) / gyro;
  y[1] = y0[1] + vd[1] * t + (-u0[0] * (1.0 - c) + u0[1] * s) / gyro;
  y[2] = y0[2] + y0[5] * t + 0.5 * field[2] * t * t;
}

static void test_lorentz(void) {
  const double y0[6] = {1.0, 0.0, 0.0, 0.2, 0.5, 0.1};
  myb_ias15* s = make(3, 0.0, 0.0, 0);
  MT_CHECK(myb_ias15_start(s, lorentz, NULL, 0.0, y0) == MYB_OK);
  double worst = 0.0;
  for (int k = 1; k <= 20; ++k) { /* about 6 gyrations */
    const double t = 1.0 * k;
    MT_CHECK(myb_ias15_integrate(s, t) == MYB_OK);
    double exact[6];
    lorentz_exact(y0, t, exact);
    const double* y = myb_ias15_state(s);
    for (int i = 0; i < 6; ++i) {
      const double scale = i < 3 ? norm3(exact) : norm3(exact + 3);
      worst = fmax(worst, fabs(y[i] - exact[i]) / scale);
    }
  }
  printf("charged particle (v-dependent force): worst relative error %.3e over 20 time units, "
         "%lld steps\n",
         worst, (long long)myb_ias15_statistics(s).accepted);
  MT_CHECK(worst < 1e-13);
  myb_ias15_destroy(s);
}

/* ---- Blocks: two independent orbits with periods 1 : 64 in one state ----------------------------- */

static myb_status two_orbits(void* context, double t, const double* x, const double* v,
                             double* a) {
  (void)context;
  for (int body = 0; body < 2; ++body) {
    const myb_status status = two_body(NULL, t, x + 3 * body, v + 3 * body, a + 3 * body);
    if (status != MYB_OK) {
      return status;
    }
  }
  return MYB_OK;
}

static void test_blocks(void) {
  const myb_pv inner = orbit(0.2, 0.1);
  myb_pv outer = orbit(0.2, 2.0);
  for (int i = 0; i < 3; ++i) { /* 16 times the distance: 64 times the period */
    outer.r[i] *= 16.0;
    outer.v[i] *= 0.25;
  }
  /* state: x_inner, x_outer, v_inner, v_outer */
  const double y0[12] = {inner.r[0], inner.r[1], inner.r[2], outer.r[0], outer.r[1], outer.r[2],
                         inner.v[0], inner.v[1], inner.v[2], outer.v[0], outer.v[1], outer.v[2]};
  const double t_end = 3.0 * period();
  myb_pv exact_inner = {};
  myb_pv exact_outer = {};
  MT_CHECK(myb_kepler_propagate(mu, &inner, t_end, &exact_inner) == MYB_OK);
  MT_CHECK(myb_kepler_propagate(mu, &outer, t_end, &exact_outer) == MYB_OK);
  const size_t blocks[2] = {3, 0};
  long long steps[2] = {};
  for (int q = 0; q < 2; ++q) {
    myb_ias15* s = make(6, 0.0, 0.0, blocks[q]);
    MT_CHECK(myb_ias15_start(s, two_orbits, NULL, 0.0, y0) == MYB_OK);
    MT_CHECK(myb_ias15_integrate(s, t_end) == MYB_OK);
    const double* y = myb_ias15_state(s);
    const double yi[6] = {y[0], y[1], y[2], y[6], y[7], y[8]};
    const double yo[6] = {y[3], y[4], y[5], y[9], y[10], y[11]};
    const double ei = time_offset(yi, &exact_inner);
    const double eo = time_offset(yo, &exact_outer);
    steps[q] = myb_ias15_statistics(s).accepted;
    printf("blocks of %zu: inner offset %.2e s, outer %.2e s, %lld steps\n", blocks[q], ei, eo,
           steps[q]);
    MT_CHECK(ei < 1e-10 && eo < 1e-10);
    myb_ias15_destroy(s);
  }
  myb_ias15* s = NULL;
  const myb_ias15_options bad = {.block = 4};
  MT_CHECK(myb_ias15_create(6, &bad, NULL, &s) == MYB_E_DOMAIN);
}

/* ---- Backward integration and intermediate stops ------------------------------------------------- */

static void test_backward_and_stops(void) {
  const double p = period();
  const myb_pv start = orbit(0.4, 0.3);
  double y0[6];
  load(&start, y0);
  myb_ias15* s = make(3, 0.0, 0.0, 0);
  MT_CHECK(myb_ias15_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(s, 3.0 * p) == MYB_OK);
  const long long forward = myb_ias15_statistics(s).accepted;
  MT_CHECK(myb_ias15_integrate(s, 0.0) == MYB_OK);
  const double* y = myb_ias15_state(s);
  double back = 0.0;
  for (int k = 0; k < 6; ++k) {
    back = fmax(back, fabs(y[k] - y0[k]) / (k < 3 ? norm3(y0) : norm3(y0 + 3)));
  }
  printf("3 orbits forward and back: largest relative difference %.3e (%lld steps each way)\n",
         back, forward);
  MT_CHECK(back < 1e-13);
  myb_ias15_destroy(s);
  /* 10 output points per orbit through integrate(): a step cut at a stop must not shrink the
   * following ones, so a stop costs at most one extra step. */
  s = make(3, 0.0, 0.0, 0);
  MT_CHECK(myb_ias15_start(s, two_body, NULL, 0.0, y0) == MYB_OK);
  for (int k = 1; k <= 30; ++k) {
    MT_CHECK(myb_ias15_integrate(s, 3.0 * p * k / 30.0) == MYB_OK);
  }
  myb_pv exact = {};
  MT_CHECK(myb_kepler_propagate(mu, &start, 3.0 * p, &exact) == MYB_OK);
  const long long with_stops = myb_ias15_statistics(s).accepted;
  printf("with 30 stops: %lld steps, time offset %.3e s\n", with_stops,
         time_offset(myb_ias15_state(s), &exact));
  MT_CHECK(with_stops <= forward + 30);
  MT_CHECK(time_offset(myb_ias15_state(s), &exact) < 1e-9);
  myb_ias15_destroy(s);
}

/* ---- Errors ------------------------------------------------------------------------------------- */

static myb_status failing(void* context, double t, const double* x, const double* v, double* a) {
  (void)context;
  (void)x;
  (void)v;
  a[0] = 1.0;
  return t > 0.5 ? MYB_E_DOMAIN : MYB_OK;
}

static myb_status singular(void* context, double t, const double* x, const double* v, double* a) {
  (void)context;
  (void)x;
  (void)v;
  a[0] = 2.0 / ((1.0 - t) * (1.0 - t) * (1.0 - t)); /* x = 1/(1-t) */
  return MYB_OK;
}

static void test_errors(void) {
  myb_ias15* s = NULL;
  const myb_ias15_options negative = {.epsilon = -1.0};
  MT_CHECK(myb_ias15_create(1, &negative, NULL, &s) == MYB_E_DOMAIN);
  const myb_ias15_options zero = {};
  MT_CHECK(myb_ias15_create(0, &zero, NULL, &s) == MYB_E_DOMAIN);
  const myb_ias15_options good = {.h_min = 1e-12};
  MT_CHECK(myb_ias15_create(1, &good, NULL, &s) == MYB_OK);
  const double y0[2] = {1.0, 1.0};
  MT_CHECK(myb_ias15_start(s, failing, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(s, 1.0) == MYB_E_DOMAIN);
  MT_CHECK(myb_ias15_start(s, singular, NULL, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(s, 2.0) == MYB_E_STEP_MIN);
  MT_CHECK(myb_ias15_time(s) < 1.0);
  const double nan_state[2] = {(double)NAN, 0.0};
  MT_CHECK(myb_ias15_start(s, singular, NULL, 0.0, nan_state) == MYB_E_DOMAIN);
  myb_ias15_destroy(s);
  myb_ias15_destroy(NULL);
}

int main(void) {
  test_constants();
  test_order();
  test_two_body();
  test_brouwer();
  test_dense_and_events();
  test_lorentz();
  test_blocks();
  test_backward_and_stops();
  test_errors();
  return mt_finish("test_ias15");
}
