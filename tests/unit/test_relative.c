/* SPDX-License-Identifier: BSD-3-Clause */
/* Relative motion (myb/relative.h):
 *  - the relative state from two inertial states and back;
 *  - the linear map against the exact motion of two orbits: the error falls as the square of
 *    the separation, for a chief on a circle, an ellipse of e = 0.7 and a hyperbola; on a
 *    circle it is the closed form of Clohessy and Wiltshire;
 *  - a two-impulse rendezvous flown with the exact motion;
 *  - relative orbital elements: constant in two-body motion but for the drift of the second,
 *    and the relative position they give;
 *  - mean elements by averaging (myb/kepler.h): an orbit under J2 — the mean semi-major axis,
 *    eccentricity and inclination stay while the osculating ones swing, and the node and the
 *    mean longitude drift at the rates of the first-order theory. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/integ.h"
#include "myb/relative.h"

static const double MU = 3.986004418e14;
static const double PI = 3.141592653589793;

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static myb_pv orbit(double a, double e, double i, double raan, double argp, double nu) {
  const myb_elements el = {.p = a * (1.0 - e * e), .e = e, .i = i, .raan = raan, .argp = argp,
                           .nu = nu};
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(MU, &el, &pv) == MYB_OK);
  return pv;
}

/* The relative state after dt by the exact motion of both, and by the linear map: the
 * difference in position (m). */
static double map_error(const myb_pv* chief, const double relative[6], double dt) {
  myb_pv deputy;
  MT_CHECK(myb_relative_to_inertial(chief, relative, &deputy) == MYB_OK);
  myb_pv chief1;
  myb_pv deputy1;
  MT_CHECK(myb_kepler_propagate(MU, chief, dt, &chief1) == MYB_OK);
  MT_CHECK(myb_kepler_propagate(MU, &deputy, dt, &deputy1) == MYB_OK);
  double exact[6];
  MT_CHECK(myb_relative_from_inertial(&chief1, &deputy1, exact) == MYB_OK);
  double phi[36];
  MT_CHECK(myb_relative_stm(MU, chief, dt, phi) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    double linear = 0.0;
    for (int j = 0; j < 6; ++j) {
      linear += phi[6 * i + j] * relative[j];
    }
    worst = fmax(worst, fabs(linear - exact[i]));
  }
  return worst;
}

static void test_frame(void) {
  const myb_pv chief = orbit(8.0e6, 0.3, 0.9, 1.0, 2.0, 0.7);
  const double relative[6] = {120.0, -3400.0, 560.0, 0.3, -0.2, 1.1};
  myb_pv deputy;
  MT_CHECK(myb_relative_to_inertial(&chief, relative, &deputy) == MYB_OK);
  double back[6];
  MT_CHECK(myb_relative_from_inertial(&chief, &deputy, back) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 6; ++i) {
    worst = fmax(worst, fabs(back[i] - relative[i]) / fabs(relative[i]));
  }
  /* the separation keeps its length, and the chief itself is at rest at the origin */
  double d[3];
  for (int i = 0; i < 3; ++i) {
    d[i] = deputy.r[i] - chief.r[i];
  }
  double self[6];
  MT_CHECK(myb_relative_from_inertial(&chief, &chief, self) == MYB_OK);
  printf("relative state: there and back %.2e; |separation| %.2e\n", worst,
         fabs(norm3(d) / norm3(relative) - 1.0));
  MT_CHECK(worst < 1e-9 && fabs(norm3(d) / norm3(relative) - 1.0) < 1e-12);
  for (int i = 0; i < 6; ++i) {
    MT_CHECK(self[i] == 0.0);
  }
  /* a deputy 1 km ahead on the same circle is at rest in the turning frame */
  const myb_pv circle = orbit(7.0e6, 0.0, 0.5, 0.2, 0.0, 0.0);
  const myb_pv ahead = orbit(7.0e6, 0.0, 0.5, 0.2, 0.0, 1000.0 / 7.0e6);
  MT_CHECK(myb_relative_from_inertial(&circle, &ahead, back) == MYB_OK);
  MT_CHECK(fabs(back[1] - 1000.0) < 1e-4 && fabs(back[3]) < 1e-9 && fabs(back[4]) < 1e-9);
  const myb_pv centre = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
  MT_CHECK(myb_relative_from_inertial(&centre, &chief, back) == MYB_E_DOMAIN);
}

static void test_linear_map(void) {
  const struct {
    const char* name;
    myb_pv chief;
    double dt;
  } cases[3] = {{"circle, 0.8 of a period", orbit(7.0e6, 0.0, 0.9, 1.0, 0.0, 0.3), 4660.0},
                {"ellipse e = 0.7 through the pericentre", orbit(2.4e7, 0.7, 1.1, 0.4, 2.0, -1.8),
                 9000.0},
                {"hyperbola e = 1.5", orbit(-1.5e7, 1.5, 0.6, 2.0, 1.0, -1.0), 5000.0}};
  for (int k = 0; k < 3; ++k) {
    const double shape[6] = {0.6, -1.0, 0.4, 3.0e-4, 5.0e-4, -8.0e-4};
    double errors[3];
    for (int s = 0; s < 3; ++s) {
      const double size = 10.0 * pow(10.0, s); /* 10 m, 100 m, 1 km */
      double relative[6];
      for (int i = 0; i < 6; ++i) {
        relative[i] = size * shape[i];
      }
      errors[s] = map_error(&cases[k].chief, relative, cases[k].dt);
    }
    printf("linear map, %-40s: separations 10 m, 100 m, 1 km leave %.2e, %.2e, %.2e m "
           "(ratios %.0f, %.0f)\n",
           cases[k].name, errors[0], errors[1], errors[2], errors[1] / errors[0],
           errors[2] / errors[1]);
    MT_CHECK(errors[2] < 100.0 && errors[0] < 0.01);
    MT_CHECK(fabs(errors[2] / errors[1] - 100.0) < 15.0);
  }
  /* on a circle: Clohessy-Wiltshire */
  const myb_pv circle = orbit(7.0e6, 0.0, 0.9, 1.0, 0.0, 0.3);
  const double n = sqrt(MU / (7.0e6 * 7.0e6 * 7.0e6));
  double phi[36];
  double cw[36];
  MT_CHECK(myb_relative_stm(MU, &circle, 4660.0, phi) == MYB_OK);
  myb_cw_stm(n, 4660.0, cw);
  double worst = 0.0;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      /* the blocks have the units 1, s, 1/s, 1 */
      const double unit = (i < 3 ? 1.0 : 1.0 / n) * (j < 3 ? 1.0 : n);
      worst = fmax(worst, fabs(phi[6 * i + j] - cw[6 * i + j]) * unit);
    }
  }
  printf("  on a circle against Clohessy-Wiltshire: %.2e\n", worst);
  MT_CHECK(worst < 1e-10);
}

static void test_rendezvous(void) {
  /* a deputy 5 km behind and 300 m below a chief on an ellipse of e = 0.2; to arrive 100 m
   * behind it, at rest, after 40 minutes */
  const myb_pv chief = orbit(7.5e6, 0.2, 0.9, 1.0, 0.5, 1.0);
  const double relative[6] = {-300.0, -5000.0, 80.0, 0.0, 0.5, 0.0};
  const double target[3] = {0.0, -100.0, 0.0};
  const double dt = 2400.0;
  double dv0[3];
  double dv1[3];
  MT_CHECK(myb_relative_rendezvous(MU, &chief, dt, relative, target, NULL, dv0, dv1) == MYB_OK);
  /* flown exactly */
  double after[6];
  memcpy(after, relative, sizeof after);
  for (int i = 0; i < 3; ++i) {
    after[3 + i] += dv0[i];
  }
  myb_pv deputy;
  myb_pv chief1;
  myb_pv deputy1;
  MT_CHECK(myb_relative_to_inertial(&chief, after, &deputy) == MYB_OK);
  MT_CHECK(myb_kepler_propagate(MU, &chief, dt, &chief1) == MYB_OK);
  MT_CHECK(myb_kepler_propagate(MU, &deputy, dt, &deputy1) == MYB_OK);
  double arrived[6];
  MT_CHECK(myb_relative_from_inertial(&chief1, &deputy1, arrived) == MYB_OK);
  double miss = 0.0;
  double speed = 0.0;
  for (int i = 0; i < 3; ++i) {
    miss = fmax(miss, fabs(arrived[i] - target[i]));
    speed = fmax(speed, fabs(arrived[3 + i] + dv1[i]));
  }
  printf("rendezvous from 5 km in 40 minutes: impulses %.3f and %.3f m/s; flown exactly it "
         "misses by %.2f m and is left with %.1e m/s\n",
         norm3(dv0), norm3(dv1), miss, speed);
  /* the linear map is good to (5 km / 7500 km) of 5 km, carried over 40 minutes: metres */
  MT_CHECK(miss < 30.0 && speed < 0.05);
  /* after a whole revolution of a circle the position does not depend on every velocity */
  const myb_pv circle = orbit(7.0e6, 0.0, 0.9, 1.0, 0.0, 0.3);
  const double period = 2.0 * PI * sqrt(7.0e6 * 7.0e6 * 7.0e6 / MU);
  MT_CHECK(myb_relative_rendezvous(MU, &circle, period, relative, target, NULL, dv0, dv1) ==
           MYB_E_DOMAIN);
}

static void test_roe(void) {
  /* a deputy with slightly different elements: its relative elements stay as they are while
   * both orbits go round, but for the drift of dlambda by -3/2 da per radian */
  const double a = 7.0e6;
  const myb_pv chief = orbit(a, 0.001, 0.9, 1.0, 0.7, 0.2);
  const myb_pv deputy = orbit(a + 30.0, 0.0012, 0.9 + 2.0e-5, 1.0 - 3.0e-5, 0.7 + 0.01, 0.2 - 0.0095);
  double roe0[6];
  MT_CHECK(myb_roe_from_states(MU, &chief, &deputy, roe0) == MYB_OK);
  const double n = sqrt(MU / (a * a * a));
  double worst = 0.0;
  double worst_position = 0.0;
  for (int k = 1; k <= 12; ++k) {
    const double dt = 500.0 * k;
    myb_pv c1;
    myb_pv d1;
    MT_CHECK(myb_kepler_propagate(MU, &chief, dt, &c1) == MYB_OK);
    MT_CHECK(myb_kepler_propagate(MU, &deputy, dt, &d1) == MYB_OK);
    double roe[6];
    MT_CHECK(myb_roe_from_states(MU, &c1, &d1, roe) == MYB_OK);
    for (int i = 0; i < 6; ++i) {
      const double expected = roe0[i] + (i == 1 ? -1.5 * roe0[0] * n * dt : 0.0);
      worst = fmax(worst, fabs(roe[i] - expected));
    }
    /* the position they give against the exact one */
    myb_elements el;
    MT_CHECK(myb_elements_from_pv(MU, &c1, &el) == MYB_OK);
    double eccentric = 0.0;
    double mean = 0.0;
    MT_CHECK(myb_anomaly_eccentric_from_true(el.nu, el.e, &eccentric) == MYB_OK);
    MT_CHECK(myb_anomaly_mean_from_eccentric(eccentric, el.e, &mean) == MYB_OK);
    double position[3];
    myb_roe_to_position(a, el.argp + mean, roe, position);
    double exact[6];
    MT_CHECK(myb_relative_from_inertial(&c1, &d1, exact) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      worst_position = fmax(worst_position, fabs(position[i] - exact[i]));
    }
  }
  printf("relative elements (da %.2e, dlambda %.2e, de %.2e, di %.2e): constant to %.1e over "
         "a revolution; the position they give within %.2f m of the exact one (separation "
         "%.0f m, e = 0.001)\n",
         roe0[0], roe0[1], hypot(roe0[2], roe0[3]), hypot(roe0[4], roe0[5]), worst,
         worst_position, a * hypot(roe0[1], hypot(roe0[2], roe0[3])));
  MT_CHECK(worst < 1e-9);
  MT_CHECK(worst_position < 10.0);
  const myb_pv open = orbit(-1.5e7, 1.5, 0.6, 2.0, 1.0, -1.0);
  MT_CHECK(myb_roe_from_states(MU, &open, &deputy, roe0) == MYB_E_DOMAIN);
}

static const double J2 = 1.08262668e-3;
static const double RE = 6378137.0;

static myb_status oblate(void* context, double t, const double* x, const double* v, double* a) {
  (void)context;
  (void)t;
  (void)v;
  const double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  const double r = sqrt(r2);
  const double z2 = x[2] * x[2] / r2;
  const double k = MU / (r2 * r);
  const double f = 1.5 * J2 * RE * RE / r2;
  a[0] = -k * x[0] * (1.0 + f * (1.0 - 5.0 * z2));
  a[1] = -k * x[1] * (1.0 + f * (1.0 - 5.0 * z2));
  a[2] = -k * x[2] * (1.0 + f * (3.0 - 5.0 * z2));
  return MYB_OK;
}

static void test_mean_elements(void) {
  const double a0 = 7.2e6;
  const double inclination = 0.9;
  const myb_pv start = orbit(a0, 0.01, inclination, 1.0, 0.5, 0.0);
  const double period = 2.0 * PI * sqrt(a0 * a0 * a0 / MU);
  enum { PER_REV = 180, REVS = 6, WINDOWS = 3 };
  static double times[PER_REV * REVS];
  static double states[6 * PER_REV * REVS];
  myb_ias15* in = NULL;
  const myb_ias15_options o = {.epsilon = 1e-11};
  MT_CHECK(myb_ias15_create(3, &o, NULL, &in) == MYB_OK);
  double y0[6];
  memcpy(y0, start.r, 3 * sizeof(double));
  memcpy(y0 + 3, start.v, 3 * sizeof(double));
  MT_CHECK(myb_ias15_start(in, oblate, NULL, 0.0, y0) == MYB_OK);
  myb_mean_elements mean[WINDOWS];
  double swing = 0.0;
  for (int w = 0; w < WINDOWS; ++w) {
    double lo = 1e300;
    double hi = 0.0;
    for (int k = 0; k < PER_REV * REVS; ++k) {
      times[k] = period * ((double)(w * REVS) + (k + 0.5) / (double)PER_REV);
      MT_CHECK(myb_ias15_integrate(in, times[k]) == MYB_OK);
      memcpy(states + 6 * k, myb_ias15_state(in), 6 * sizeof(double));
      myb_pv pv;
      memcpy(pv.r, states + 6 * k, 3 * sizeof(double));
      memcpy(pv.v, states + 6 * k + 3, 3 * sizeof(double));
      myb_elements el;
      MT_CHECK(myb_elements_from_pv(MU, &pv, &el) == MYB_OK);
      const double a = el.p / (1.0 - el.e * el.e);
      lo = fmin(lo, a);
      hi = fmax(hi, a);
    }
    swing = fmax(swing, hi - lo);
    MT_CHECK(myb_mean_elements_average(MU, times, states, PER_REV * REVS, &mean[w]) == MYB_OK);
  }
  myb_ias15_destroy(in);
  /* the means stay from one window of six revolutions to the next */
  const double da = fmax(fabs(mean[1].a - mean[0].a), fabs(mean[2].a - mean[1].a));
  const double e0 = hypot(mean[0].f, mean[0].g);
  const double e2 = hypot(mean[2].f, mean[2].g);
  const double i0 = 2.0 * atan(hypot(mean[0].h, mean[0].k));
  const double i2 = 2.0 * atan(hypot(mean[2].h, mean[2].k));
  /* the node falls back at -3/2 J2 (R / p)^2 n cos i, and the mean longitude runs at
   * n (1 + 3/2 J2 (R / p)^2 [sqrt(1 - e^2) (1 - 3/2 sin^2 i) + 2 - 5/2 sin^2 i - cos i]) */
  const double n = sqrt(MU / (mean[0].a * mean[0].a * mean[0].a));
  const double p = mean[0].a * (1.0 - e0 * e0);
  const double k2 = 1.5 * J2 * RE * RE / (p * p);
  const double s2 = sin(i0) * sin(i0);
  const double node_rate = -k2 * n * cos(i0);
  const double lambda_rate =
      n * (1.0 + k2 * (sqrt(1.0 - e0 * e0) * (1.0 - 1.5 * s2) + 2.0 - 2.5 * s2 - cos(i0)));
  const double node0 = atan2(mean[0].k, mean[0].h);
  const double node2 = atan2(mean[2].k, mean[2].h);
  const double measured = (node2 - node0) / (mean[2].epoch - mean[0].epoch);
  printf("mean elements under J2 (the osculating a swings by %.0f m): the mean a keeps to "
         "%.2f m, e to %.1e, i to %.1e rad; the node drifts at %.6e rad/s (first order "
         "%.6e), the mean longitude at %.9e (%.9e)\n",
         swing, da, fabs(e2 - e0), fabs(i2 - i0), measured, node_rate, mean[1].rate,
         lambda_rate);
  MT_CHECK(swing > 5.0e3);
  MT_CHECK(da < 1.0 && fabs(e2 - e0) < 3e-6 && fabs(i2 - i0) < 1e-7);
  MT_CHECK(fabs(measured / node_rate - 1.0) < 3e-3);
  MT_CHECK(fabs(mean[1].rate / lambda_rate - 1.0) < 3e-6);
  MT_CHECK(fabs(mean[1].epoch - period * 9.0) < 1e-6);
  /* refused */
  MT_CHECK(myb_mean_elements_average(MU, times, states, 2, &mean[0]) == MYB_E_DOMAIN);
  times[1] = times[0];
  MT_CHECK(myb_mean_elements_average(MU, times, states, 10, &mean[0]) == MYB_E_DOMAIN);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_frame();
  test_linear_map();
  test_rendezvous();
  test_roe();
  test_mean_elements();
  return mt_finish("test_relative");
}
