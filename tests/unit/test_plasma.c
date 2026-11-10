/* SPDX-License-Identifier: BSD-3-Clause */
/* Electron columns of power-law plasma profiles (myb_plasma_column) against closed forms and
 * mpmath quadrature (40 digits, in the path length and in the angle, which agree to 25 digits),
 * and the first-order dispersive delay (myb_plasma_group_delay). */
#include <math.h>
#include <stdio.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/observe.h"

#define AU 1.495978707e11
#define RSUN 6.957e8

static double rel(double a, double b) {
  return fabs(a - b) / fabs(b);
}

static double angle(const double* p1, const double* p2, double* b_out, double* t1, double* t2) {
  const double d[3] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
  const double len = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  const double u[3] = {d[0] / len, d[1] / len, d[2] / len};
  const double c[3] = {p1[1] * u[2] - p1[2] * u[1], p1[2] * u[0] - p1[0] * u[2],
                       p1[0] * u[1] - p1[1] * u[0]};
  const double b = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
  const double s0 = -(p1[0] * u[0] + p1[1] * u[1] + p1[2] * u[2]);
  *b_out = b;
  *t1 = myb_atan2(-s0, b);
  *t2 = myb_atan2(len - s0, b);
  return len;
}

static void test_reference(void) {
  /* a corona of Allen-Baumbach type plus a wind ~ r^-2.3, Earth at 1 AU to a spacecraft at
   * 1.5 AU behind the Sun, closest approach 5 solar radii */
  myb_plasma_profile prof = {.radius = RSUN, .count = 3,
                             .density = {2.99e14, 1.55e14, 3.44e11},
                             .power = {16.0, 6.0, 2.3}};
  const double p1[3] = {AU, 0.0, 0.0};
  const double p2[3] = {-1.5 * AU, 5.0 * RSUN, 0.3 * RSUN};
  double n = 0.0;
  MT_CHECK(myb_plasma_column(&prof, p1, p2, &n) == MYB_OK);
  printf("corona + wind, 1 AU -> 1.5 AU past 5 solar radii: %.15e m^-2 (mpmath "
         "4.194194573219541213889e21), relative %.1e\n",
         n, rel(n, 4.194194573219541213889438e21));
  MT_CHECK(rel(n, 4.194194573219541213889438e21) < 1e-12);
  myb_plasma_profile wind = {.radius = RSUN, .count = 1, .density = {3.44e11}, .power = {2.3}};
  const double p3[3] = {-1.2 * AU, 50.0 * RSUN, 0.0};
  MT_CHECK(myb_plasma_column(&wind, p1, p3, &n) == MYB_OK);
  printf("wind r^-2.3 past 50 solar radii: %.15e (mpmath 1.062798097046474398483612e19), "
         "relative %.1e\n",
         n, rel(n, 1.062798097046474398483612e19));
  MT_CHECK(rel(n, 1.062798097046474398483612e19) < 1e-12);
}

static void test_closed_forms(void) {
  const double p1[3] = {AU, 0.2 * AU, -0.1 * AU};
  const double p2[3] = {-0.7 * AU, -0.1 * AU, 0.05 * AU};
  double b = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  (void)angle(p1, p2, &b, &t1, &t2);
  /* n (R/r)^2: n R^2 (t2 - t1) / b */
  myb_plasma_profile prof = {.radius = RSUN, .count = 1, .density = {1e12}, .power = {2.0}};
  double n = 0.0;
  MT_CHECK(myb_plasma_column(&prof, p1, p2, &n) == MYB_OK);
  const double e2 = 1e12 * RSUN * RSUN * (t2 - t1) / b;
  MT_CHECK(rel(n, e2) < 1e-14);
  /* n (R/r)^4: n R^4 / b^3 [t/2 + sin 2t / 4] */
  prof.power[0] = 4.0;
  MT_CHECK(myb_plasma_column(&prof, p1, p2, &n) == MYB_OK);
  const double r4 = RSUN * RSUN * RSUN * RSUN;
  const double e4 = 1e12 * r4 / (b * b * b) *
                    ((t2 - t1) / 2.0 + (myb_sin(2.0 * t2) - myb_sin(2.0 * t1)) / 4.0);
  MT_CHECK(rel(n, e4) < 1e-13);
  printf("closed forms r^-2, r^-4: %.1e, %.1e\n", rel(n, e4), rel(n, e4));
  /* radial paths: n R^2 (1 / r1 - 1 / r2) and n R log(r2 / r1) */
  const double q1[3] = {2.0 * RSUN, 0.0, 0.0};
  const double q2[3] = {300.0 * RSUN, 0.0, 0.0};
  prof.power[0] = 2.0;
  MT_CHECK(myb_plasma_column(&prof, q1, q2, &n) == MYB_OK);
  MT_CHECK(rel(n, 1e12 * RSUN * RSUN * (1.0 / (2.0 * RSUN) - 1.0 / (300.0 * RSUN))) < 1e-14);
  prof.power[0] = 1.0;
  MT_CHECK(myb_plasma_column(&prof, q2, q1, &n) == MYB_OK);
  MT_CHECK(rel(n, 1e12 * RSUN * myb_log(150.0)) < 1e-14);
  /* through the centre: no answer */
  const double q3[3] = {-3.0 * RSUN, 0.0, 0.0};
  MT_CHECK(myb_plasma_column(&prof, q1, q3, &n) == MYB_E_DOMAIN);
  prof.power[0] = -1.0;
  MT_CHECK(myb_plasma_column(&prof, p1, p2, &n) == MYB_E_DOMAIN);
}

static void test_delay(void) {
  /* CODATA 2022: K = e^2 / (8 pi^2 eps0 m_e) = 40.308 193 m^3 s^-2 */
  const double e = 1.602176634e-19;
  const double me = 9.1093837139e-31;
  const double eps0 = 8.8541878188e-12;
  const double k = myb_plasma_group_delay(1.0, 1.0, e, me, eps0);
  printf("K = %.9f m^3/s^2; 10 TECU at 8.4 GHz: %.6f m\n", k,
         myb_plasma_group_delay(1e17, 8.4e9, e, me, eps0));
  MT_CHECK(fabs(k - 40.308193) < 1e-6);
  MT_CHECK(fabs(myb_plasma_group_delay(1e17, 8.4e9, e, me, eps0) - k * 1e17 / 7.056e19) < 1e-15);
  MT_CHECK(isnan(myb_plasma_group_delay(1e17, 0.0, e, me, eps0)));
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_reference();
  test_closed_forms();
  test_delay();
  return mt_finish("test_plasma");
}
