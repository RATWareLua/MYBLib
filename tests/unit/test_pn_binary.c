/* SPDX-License-Identifier: BSD-3-Clause */
/* Compact binaries beyond 1PN (myb_pn_binary_acceleration):
 *   - the gauge-invariant ratio K = Omega_phi / Omega_r of a circular orbit, K = 1 + 3x +
 *     (27/2 - 7 nu) x^2 + O(x^3) with x = (G m Omega_phi / c^3)^(2/3) (Damour, Jaranowski,
 *     Schafer 2000; the test-mass limit is (1 - 6x)^(-1/2) of Schwarzschild), measured on
 *     near-circular orbits integrated with the 1PN + 2PN equations;
 *   - the circular inspiral by radiation reaction against Peters' dOmega/dt;
 *   - PSR B1913+16 (Hulse-Taylor): the period derivative from the times of periastron over
 *     200 orbits of the 2.5PN equations against Peters-Mathews (-2.4026e-12 for the masses of
 *     Weisberg & Huang 2016). */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "myb/num.h"
#include "myb/rel.h"

typedef struct binary {
  double gm1;
  double gm2;
  double c;
  unsigned orders;
} binary;

static myb_status rhs(void* context, double t, const double* y, double* dy) {
  (void)t;
  const binary* b = context;
  dy[0] = y[3];
  dy[1] = y[4];
  dy[2] = y[5];
  return myb_pn_binary_acceleration(b->gm1, b->gm2, b->c, b->orders, y, y + 3, dy + 3);
}

static double radial_velocity(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[0] * y[3] + y[1] * y[4] + y[2] * y[5];
}

static myb_status accel(void* context, double t, const double* x, const double* v, double* a) {
  (void)t;
  const binary* b = context;
  return myb_pn_binary_acceleration(b->gm1, b->gm2, b->c, b->orders, x, v, a);
}

/* Integrates with IAS15 (no secular energy drift: Brouwer's law) until `count` periastra
 * (r.v from - to +) and returns their times and the unwrapped longitudes there. */
static double extreme_radius[2]; /* the first periastron and apastron radii found */

static int periastra(binary* b, const double y0[6], int count, double* times,
                     double* longitudes) {
  extreme_radius[0] = extreme_radius[1] = 0.0;
  const myb_ias15_options o = {.epsilon = 1e-9};
  myb_ias15* d = NULL;
  MT_CHECK(myb_ias15_create(3, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_ias15_start(d, accel, b, 0.0, y0) == MYB_OK);
  double phi = myb_atan2(y0[1], y0[0]);
  double last = phi;
  int found_count = 0;
  double prev = radial_velocity(NULL, 0.0, y0);
  for (long k = 0; k < 50000000L && found_count < count; ++k) {
    if (myb_ias15_step(d, 1e300) != MYB_OK) {
      break;
    }
    const double* y = myb_ias15_state(d);
    const double now = radial_velocity(NULL, 0.0, y);
    if (prev > 0.0 && now <= 0.0 && extreme_radius[1] == 0.0) {
      int found = 0;
      double te = 0.0;
      double ye[6];
      MT_CHECK(myb_ias15_locate(d, radial_velocity, NULL, &found, &te, ye) == MYB_OK);
      if (found) {
        extreme_radius[1] = sqrt(ye[0] * ye[0] + ye[1] * ye[1] + ye[2] * ye[2]);
      }
    }
    if (prev < 0.0 && now >= 0.0) {
      int found = 0;
      double te = 0.0;
      double ye[6];
      MT_CHECK(myb_ias15_locate(d, radial_velocity, NULL, &found, &te, ye) == MYB_OK);
      if (found && extreme_radius[0] == 0.0) {
        extreme_radius[0] = sqrt(ye[0] * ye[0] + ye[1] * ye[1] + ye[2] * ye[2]);
      }
      if (found) {
        double a = myb_atan2(ye[1], ye[0]);
        double da = a - last;
        da -= 6.283185307179586 * floor(da / 6.283185307179586 + 0.5);
        times[found_count] = te;
        longitudes[found_count] = phi + da;
        ++found_count;
      }
    }
    /* unwrap the longitude along the steps */
    double a = myb_atan2(y[1], y[0]);
    double da = a - last;
    da -= 6.283185307179586 * floor(da / 6.283185307179586 + 0.5);
    phi += da;
    last = a;
    prev = now;
  }
  myb_ias15_destroy(d);
  return found_count;
}

/* Speed of the circular orbit of radius r at 2PN (rdot = 0): v^2 = (gm/r)(1 + A(v^2)). */
static double circular_speed(double gm, double nu, double c, double r) {
  const double u = gm / r;
  double v2 = u;
  for (int k = 0; k < 50; ++k) {
    const double a1 = (v2 + 3.0 * nu * v2 - u * (4.0 + 2.0 * nu)) / (c * c);
    const double a2 = (3.0 * nu * v2 * v2 - 4.0 * nu * nu * v2 * v2 +
                       u * (-6.5 * nu * v2 + 2.0 * nu * nu * v2) + u * u * (9.0 + 87.0 / 4.0 * nu)) /
                      (c * c * c * c);
    v2 = u * (1.0 + a1 + a2);
  }
  return sqrt(v2);
}

static void test_periastron_advance(void) {
  /* units G m = c = 1 */
  const double nus[2] = {0.25, 0.01};
  const double xs[2] = {1e-3, 3e-4};
  for (int a = 0; a < 2; ++a) {
    const double nu = nus[a];
    /* gm1 gm2 / (gm1 + gm2)^2 = nu, gm1 + gm2 = 1 */
    const double d = sqrt(1.0 - 4.0 * nu);
    binary b = {.gm1 = 0.5 * (1.0 + d), .gm2 = 0.5 * (1.0 - d), .c = 1.0,
                .orders = MYB_PN_1 | MYB_PN_2};
    for (int k = 0; k < 2; ++k) {
      const double r0 = 1.0 / xs[k];
      const double v = circular_speed(1.0, nu, 1.0, r0) * (1.0 + 1e-4);
      const double y0[6] = {r0, 0.0, 0.0, 0.0, v, 0.0};
      double times[3];
      double lon[3];
      MT_CHECK(periastra(&b, y0, 3, times, lon) == 3);
      const double pr = times[2] - times[1];
      const double dphi = lon[2] - lon[1];
      const double kk = dphi / 6.283185307179586;
      const double omega = dphi / pr;
      const double x = myb_cbrt(omega * omega);
      const double second = (kk - 1.0 - 3.0 * x) / (x * x);
      const double expected = 13.5 - 7.0 * nu;
      printf("nu %.2f, x %.1e: K = %.12f, (K - 1 - 3x) / x^2 = %.5f (27/2 - 7 nu = %.5f), "
             "difference %.2e\n",
             nu, x, kk, second, expected, fabs(second - expected));
      /* the remainder is O(x): truncation of the 2PN equations and the residual eccentricity */
      MT_CHECK(fabs(second - expected) < 120.0 * x);
    }
  }
}

static void test_circular_inspiral(void) {
  const double nu = 0.25;
  binary b = {.gm1 = 0.5, .gm2 = 0.5, .c = 1.0, .orders = MYB_PN_1 | MYB_PN_2 | MYB_PN_25};
  const double r0 = 2e4; /* x = 5e-5 */
  const double v = circular_speed(1.0, nu, 1.0, r0);
  const double y0[6] = {r0, 0.0, 0.0, 0.0, v, 0.0};
  /* the mean angular velocity over successive turns from the integration */
  const myb_dop853_options o = {.rtol = 1e-15, .atol = 1e-18, .h_max = 1000.0};
  myb_dop853* d = NULL;
  MT_CHECK(myb_dop853_create(6, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_dop853_start(d, rhs, &b, 0.0, y0) == MYB_OK);
  const double period = 6.283185307179586 * r0 * sqrt(r0);
  const int turns = 20;
  double omega[2];
  double tmid[2];
  for (int s = 0; s < 2; ++s) {
    /* angular velocity at two instants, half a turn averaged: |x x v| / r^2 */
    const double t = s == 0 ? 2.0 * period : (2.0 + turns) * period;
    MT_CHECK(myb_dop853_integrate(d, t) == MYB_OK);
    const double* y = myb_dop853_state(d);
    const double hz = y[0] * y[4] - y[1] * y[3];
    const double rr = y[0] * y[0] + y[1] * y[1];
    omega[s] = hz / rr;
    tmid[s] = t;
  }
  myb_dop853_destroy(d);
  const double rate = (omega[1] - omega[0]) / (tmid[1] - tmid[0]);
  const double om = 0.5 * (omega[0] + omega[1]);
  const double peters = 96.0 / 5.0 * nu * myb_exp(11.0 / 3.0 * myb_log(om));
  printf("circular inspiral, x = 5e-5: dOmega/dt %.6e, Peters %.6e, ratio %.6f\n", rate, peters,
         rate / peters);
  MT_CHECK(fabs(rate / peters - 1.0) < 2e-3);
}

static void test_hulse_taylor(void) {
  const double gm_sun = 1.32712440041e20;
  const double c = 299792458.0;
  const double m1 = 1.438;
  const double m2 = 1.390;
  const double pb = 0.322997448918 * 86400.0;
  const double e = 0.6171340;
  binary b = {.gm1 = m1 * gm_sun, .gm2 = m2 * gm_sun, .c = c,
              .orders = MYB_PN_1 | MYB_PN_2 | MYB_PN_25};
  const double gm = b.gm1 + b.gm2;
  const double n = 6.283185307179586 / pb;
  const double a = myb_cbrt(gm / (n * n));
  const double rp = a * (1.0 - e);
  const double vp = sqrt(gm / a * (1.0 + e) / (1.0 - e));
  const double y0[6] = {rp, 0.0, 0.0, 0.0, vp, 0.0};
  enum { ORBITS = 200 };
  static double times[ORBITS];
  static double lon[ORBITS];
  MT_CHECK(periastra(&b, y0, ORBITS, times, lon) == ORBITS);
  /* least squares t_k = t0 + P k + q k^2 */
  double ata[9] = {0};
  double atb[3] = {0};
  for (int k = 0; k < ORBITS; ++k) {
    const double row[3] = {1.0, (double)k, (double)k * (double)k};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        ata[3 * i + j] += row[i] * row[j];
      }
      atb[i] += row[i] * (times[k] - times[0] - pb * k);
    }
  }
  MT_CHECK(myb_lu_solve(3, ata, atb, 1) == MYB_OK);
  const double p0 = pb + atb[1];
  const double pdot = 2.0 * atb[2] / p0;
  /* the eccentricity of the integrated (1PN) orbit from its radii, not the Newtonian one of the
   * initial state: at periastron (v/c)^2 v^2 / |E| ~ 1e-4 changes it by ~4e-5, which f(e) of
   * Peters-Mathews (d ln f / de = 8.8) would turn into 3.5e-4 */
  const double e_orbit = (extreme_radius[1] - extreme_radius[0]) /
                         (extreme_radius[1] + extreme_radius[0]);
  const double peters = myb_peters_period_derivative(b.gm1, b.gm2, c, p0, e_orbit);
  const double peters_nominal = myb_peters_period_derivative(b.gm1, b.gm2, c, pb, e);
  /* the periastron advance per year */
  const double advance = (lon[ORBITS - 1] - lon[0]) / (times[ORBITS - 1] - times[0]) -
                         6.283185307179586 / p0;
  const double deg_per_year = advance * 365.25 * 86400.0 * 180.0 / 3.141592653589793;
  printf("PSR B1913+16, 200 orbits: dP/dt %.6e, Peters-Mathews %.6e for the orbit's e = %.7f "
         "(ratio %.6f); nominal %.5e; periastron advance %.4f deg/yr (observed 4.2266)\n",
         pdot, peters, e_orbit, pdot / peters, peters_nominal, deg_per_year);
  MT_CHECK(fabs(pdot / peters - 1.0) < 5e-5);
  /* the GR prediction -2.40263e-12 of Weisberg & Huang (2016) uses unrounded masses; with the
   * rounded 1.438 and 1.390 it is 0.05 % smaller */
  MT_CHECK(fabs(peters_nominal + 2.4026e-12) < 0.002e-12);
  MT_CHECK(fabs(deg_per_year - 4.2266) < 0.003);
}

static void test_domain(void) {
  const double x[3] = {1.0, 0.0, 0.0};
  const double v[3] = {0.0, 1.0, 0.0};
  const double z[3] = {0.0, 0.0, 0.0};
  double a[3];
  MT_CHECK(myb_pn_binary_acceleration(1.0, 1.0, 1.0, 0, z, v, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_pn_binary_acceleration(0.0, 1.0, 1.0, 0, x, v, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_pn_binary_acceleration(1.0, 1.0, 0.0, 0, x, v, a) == MYB_E_DOMAIN);
  /* Newton */
  MT_CHECK(myb_pn_binary_acceleration(0.3, 0.7, 1.0, 0, x, v, a) == MYB_OK);
  MT_CHECK(a[0] == -1.0 && a[1] == 0.0 && a[2] == 0.0);
  MT_CHECK(isnan(myb_peters_period_derivative(1.0, 1.0, 1.0, 1.0, 1.0)));
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_periastron_advance();
  test_circular_inspiral();
  test_hulse_taylor();
  test_domain();
  return mt_finish("test_pn_binary");
}
