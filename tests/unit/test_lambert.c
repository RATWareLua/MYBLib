/* SPDX-License-Identifier: BSD-3-Clause */
/* Lambert's problem: every solution is verified independently by propagating (r1, v1) over
 * the transfer time with the universal-variable propagator, which must arrive at r2 with v2;
 * the direction of motion and the number of revolutions are checked too. */
#include <math.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/kepler.h"
#include "myb/lambert.h"

static double norm3(const double a[3]) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

static double rel_diff3(const double a[3], const double b[3]) {
  const double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  return norm3(d) / norm3(b);
}

static void test_curtis(void) {
  /* H. D. Curtis, Orbital Mechanics for Engineering Students, Example 5.2 (km, km/s, 1 h). */
  const double mu = 398600e9;
  const double r1[3] = {5000e3, 10000e3, 2100e3};
  const double r2[3] = {-14600e3, 2500e3, 7000e3};
  myb_lambert_solution solutions[1];
  size_t count = 0;
  MT_CHECK(myb_lambert(mu, r1, r2, 3600.0, 0, 0, solutions, 1, &count) == MYB_OK);
  MT_CHECK(count == 1);
  const double v1[3] = {-5.9925e3, 1.9254e3, 3.2456e3};
  const double v2[3] = {-3.3125e3, -4.1966e3, -0.38529e3};
  for (int k = 0; k < 3; ++k) {
    MT_CHECK(fabs(solutions[0].v1[k] - v1[k]) < 0.06);
    MT_CHECK(fabs(solutions[0].v2[k] - v2[k]) < 0.06);
  }
  printf("Curtis 5.2: v1 = (%.4f, %.4f, %.4f) km/s, v2 = (%.4f, %.4f, %.5f) km/s\n",
         solutions[0].v1[0] / 1e3, solutions[0].v1[1] / 1e3, solutions[0].v1[2] / 1e3,
         solutions[0].v2[0] / 1e3, solutions[0].v2[1] / 1e3, solutions[0].v2[2] / 1e3);
}

static void random_vector(mt_rng* rng, double magnitude, double out[3]) {
  double v[3];
  double n = 0.0;
  /* bound: rejection of the (rare) tiny vectors, at most 100 tries */
  for (int i = 0; i < 100 && !(n > 0.1); ++i) {
    for (int k = 0; k < 3; ++k) {
      v[k] = 2.0 * (double)(mt_next(rng) >> 11) * 0x1p-53 - 1.0;
    }
    n = norm3(v);
  }
  for (int k = 0; k < 3; ++k) {
    out[k] = v[k] / n * magnitude;
  }
}

static void test_random(void) {
  const double mu = 3.986004418e14;
  mt_rng rng = {.state = 901};
  double worst_r = 0.0;
  double worst_v = 0.0;
  long solutions_checked = 0;
  long multi = 0;
  for (int i = 0; i < 20000; ++i) {
    double r1[3];
    double r2[3];
    random_vector(&rng, 6.6e6 * (1.0 + 5.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53), r1);
    random_vector(&rng, 6.6e6 * (1.0 + 5.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53), r2);
    double cross[3] = {r1[1] * r2[2] - r1[2] * r2[1], r1[2] * r2[0] - r1[0] * r2[2],
                       r1[0] * r2[1] - r1[1] * r2[0]};
    if (norm3(cross) < 0.02 * norm3(r1) * norm3(r2)) {
      continue; /* nearly collinear: the transfer plane is ill-defined */
    }
    const double tof = 300.0 * exp((double)(mt_next(&rng) >> 11) * 0x1p-53 * 7.0);
    const int retrograde = (int)(mt_next(&rng) & 1);
    myb_lambert_solution solutions[7];
    size_t count = 0;
    MT_CHECK(myb_lambert(mu, r1, r2, tof, retrograde, 3, solutions, 7, &count) == MYB_OK);
    MT_CHECK(count >= 1 && count % 2 == 1);
    for (size_t j = 0; j < count; ++j) {
      const myb_lambert_solution* sol = &solutions[j];
      myb_pv start = {.r = {r1[0], r1[1], r1[2]}, .v = {sol->v1[0], sol->v1[1], sol->v1[2]}};
      myb_pv end = {};
      MT_CHECK(myb_kepler_propagate(mu, &start, tof, &end) == MYB_OK);
      const double er = rel_diff3(end.r, r2);
      const double ev = rel_diff3(end.v, sol->v2);
      worst_r = er > worst_r ? er : worst_r;
      worst_v = ev > worst_v ? ev : worst_v;
      /* Direction of motion: h_z > 0 for prograde. */
      const double hz = r1[0] * sol->v1[1] - r1[1] * sol->v1[0];
      MT_CHECK(retrograde ? hz < 0.0 : hz > 0.0);
      /* Revolutions: an elliptic solution with M revolutions has M P < tof < (M + 1) P. */
      const double energy =
          0.5 * (sol->v1[0] * sol->v1[0] + sol->v1[1] * sol->v1[1] + sol->v1[2] * sol->v1[2]) -
          mu / norm3(r1);
      if (energy < 0.0) {
        const double a = -mu / (2.0 * energy);
        const double period = 2.0 * 3.141592653589793 * sqrt(a * a * a / mu);
        MT_CHECK((int)(tof / period) == sol->revolutions);
      } else {
        MT_CHECK(sol->revolutions == 0);
      }
      multi += sol->revolutions > 0;
      ++solutions_checked;
    }
  }
  printf(
      "Lambert: %ld solutions (%ld multi-revolution), propagated to r2: worst relative error "
      "%.3g in r, %.3g in v\n",
      solutions_checked, multi, worst_r, worst_v);
  MT_CHECK(multi > 0);
  MT_CHECK(worst_r < 1e-10 && worst_v < 1e-10);
}

static void test_errors(void) {
  const double mu = 3.986004418e14;
  const double r1[3] = {7e6, 0.0, 0.0};
  const double r2[3] = {14e6, 0.0, 0.0};
  const double r3[3] = {0.0, 7e6, 0.0};
  myb_lambert_solution solutions[3];
  size_t count = 9;
  MT_CHECK(myb_lambert(mu, r1, r2, 1000.0, 0, 1, solutions, 3, &count) == MYB_E_DOMAIN);
  MT_CHECK(count == 0);
  MT_CHECK(myb_lambert(mu, r1, r3, 0.0, 0, 1, solutions, 3, &count) == MYB_E_DOMAIN);
  MT_CHECK(myb_lambert(mu, r1, r3, -5.0, 0, 1, solutions, 3, &count) == MYB_E_DOMAIN);
  MT_CHECK(myb_lambert(-mu, r1, r3, 1000.0, 0, 1, solutions, 3, &count) == MYB_E_DOMAIN);
}

int main(void) {
  test_curtis();
  test_random();
  test_errors();
  return mt_finish("test_lambert");
}
