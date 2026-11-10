/* SPDX-License-Identifier: BSD-3-Clause */
/* ext/troposphere against the test cases published in the headers of the IERS Conventions
 * (2010) chapter 9 routines (FCUL_ZD_HPA.F, FCUL_A.F, VMF1.F, VMF1_HT.F, GMF.F, GPT.F; the
 * authors' code, printed to 19 digits), and properties of the models. */
#include <math.h>
#include <stdio.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb_troposphere.h"

#define DEG (3.14159265358979323846 / 180.0)

static double rel(double a, double b) {
  return fabs(a - b) / fabs(b);
}

static void test_iers_cases(void) {
  /* FCUL_ZD_HPA: McDonald Observatory, 14 August 2009, 532 nm. The published case is 2.0e-6
   * smaller in both components. Evaluating the routine's own statements (with its literals
   * 0.00266 and 3.759 in single precision, as Fortran reads them) makes the two ratios equal,
   * 1 + 1.96361e-6: a common factor, i.e. the term f(phi, H) of the denominator, as if the case
   * had been run with the height 7 m lower — the case does not match its own inputs. 2e-6 of
   * the delay is 4 micrometres. */
  double zhd = 0.0;
  double zwd = 0.0;
  MT_CHECK(myb_tropo_zenith_optical(30.67166667 * DEG, 2010.344, 798.4188, 14.322, 0.532, 375.0,
                                    &zhd, &zwd) == MYB_OK);
  printf("Mendes-Pavlis: ZHD %.16f (IERS 1.932992176591644462), ZWD %.16e (2.233748255158703871e-3)"
         ", ratios - 1: %.6e %.6e\n",
         zhd, zwd, zhd / 1.932992176591644462 - 1.0, zwd / 0.2233748255158703871e-2 - 1.0);
  MT_CHECK(fabs(zhd / 1.932992176591644462 - 1.0 - 1.96361e-6) < 1e-10);
  MT_CHECK(fabs(zwd / 0.2233748255158703871e-2 - 1.0 - 2.00404e-6) < 1e-10);
  MT_CHECK(rel(zhd + zwd, 1.935225924846803114) < 3e-6);

  double m = 0.0;
  MT_CHECK(myb_tropo_fcula(15.0 * DEG, 30.67166667 * DEG, 2075.0, 300.15, &m) == MYB_OK);
  printf("FCULa: %.16f (IERS 3.800243667312344087)\n", m);
  MT_CHECK(rel(m, 3.800243667312344087) < 2e-15);

  double mh = 0.0;
  double mw = 0.0;
  MT_CHECK(myb_tropo_vmf1(0.00127683, 0.00060955, 55055.0, 0.6708665767, 0.0, 1.278564131, &mh,
                          &mw) == MYB_OK);
  printf("VMF1: %.16f %.16f (IERS 3.424342122738070593 3.448299714692572238)\n", mh, mw);
  MT_CHECK(rel(mh, 3.424342122738070593) < 2e-15);
  MT_CHECK(rel(mw, 3.448299714692572238) < 2e-15);
  MT_CHECK(myb_tropo_vmf1(0.00127683, 0.00060955, 55055.0, 0.6708665767, 824.17, 1.278564131, &mh,
                          &mw) == MYB_OK);
  printf("VMF1 with height: %.16f (IERS 3.425088087972572470)\n", mh);
  MT_CHECK(rel(mh, 3.425088087972572470) < 2e-15);
  MT_CHECK(rel(mw, 3.448299714692572238) < 2e-15);

  MT_CHECK(myb_tropo_gmf(55055.0, 0.6708665767, -1.393397187, 844.715, 1.278564131, &mh, &mw) ==
           MYB_OK);
  printf("GMF: %.16f %.16f (IERS 3.425245519339138678 3.449589116182419257)\n", mh, mw);
  MT_CHECK(rel(mh, 3.425245519339138678) < 2e-15);
  MT_CHECK(rel(mw, 3.449589116182419257) < 2e-15);

  double p = 0.0;
  double t = 0.0;
  double u = 0.0;
  MT_CHECK(myb_tropo_gpt(55055.0, 0.6708665767, -1.393397187, 812.546, &p, &t, &u) == MYB_OK);
  printf("GPT: %.13f hPa, %.14f C, %.14f m (IERS 918.0710638757363995 19.31914181012882992 "
         "-42.19185643717770517)\n",
         p, t, u);
  MT_CHECK(rel(p, 918.0710638757363995) < 2e-15);
  MT_CHECK(rel(t, 19.31914181012882992) < 2e-14);
  MT_CHECK(rel(u, -42.19185643717770517) < 2e-14);
}

static void test_properties(void) {
  /* the mapping functions are 1 at the zenith (the height correction vanishes there) */
  double mh = 0.0;
  double mw = 0.0;
  MT_CHECK(myb_tropo_gmf(58000.0, -0.6, 2.0, 1500.0, 0.0, &mh, &mw) == MYB_OK);
  MT_CHECK(fabs(mh - 1.0) < 1e-15 && fabs(mw - 1.0) < 1e-15);
  MT_CHECK(myb_tropo_vmf1(0.0012, 0.0006, 58000.0, -0.6, 1500.0, 0.0, &mh, &mw) == MYB_OK);
  MT_CHECK(fabs(mh - 1.0) < 1e-15 && fabs(mw - 1.0) < 1e-15);
  double m = 0.0;
  MT_CHECK(myb_tropo_fcula(1.5707963267948966, 0.3, 100.0, 290.0, &m) == MYB_OK);
  MT_CHECK(fabs(m - 1.0) < 1e-15);
  /* and grow like 1 / sin e towards the horizon: at 5 degrees between 10 and 11.5 */
  MT_CHECK(myb_tropo_gmf(58000.0, 0.8, 0.2, 100.0, 85.0 * DEG, &mh, &mw) == MYB_OK);
  MT_CHECK(mh > 10.0 && mh < 11.5 && mw > 10.0 && mw < 11.5);
  /* Saastamoinen: 2.3070 m at sea level, 45 degrees, standard pressure */
  double zhd = 0.0;
  MT_CHECK(myb_tropo_zhd_saastamoinen(1013.25, 45.0 * DEG, 0.0, &zhd) == MYB_OK);
  MT_CHECK(fabs(zhd - 0.0022768 * 1013.25) < 1e-15);
  double zwd = 0.0;
  MT_CHECK(myb_tropo_zwd_saastamoinen(288.15, 10.0, &zwd) == MYB_OK);
  MT_CHECK(fabs(zwd - 0.002277 * (1255.0 / 288.15 + 0.05) * 10.0) < 1e-15);
  /* gradients: zero at the zenith, cos/sin of the azimuth, 1 / (sin e tan e + 0.0032) */
  double d = 0.0;
  MT_CHECK(myb_tropo_gradient_delay(0.001, 0.0005, 1.5707963267948966, 0.3, &d) == MYB_OK);
  MT_CHECK(fabs(d) < 1e-18);
  MT_CHECK(myb_tropo_gradient_delay(0.001, 0.0, 10.0 * DEG, 0.0, &d) == MYB_OK);
  const double e = 10.0 * DEG;
  MT_CHECK(fabs(d - 0.001 / (myb_sin(e) * myb_tan(e) + 0.0032)) < 1e-17);
  /* the GPT pressure at sea level is about 1013 hPa */
  double p = 0.0;
  double t = 0.0;
  double u = 0.0;
  MT_CHECK(myb_tropo_gpt(58000.0, 0.0, 0.0, 0.0, &p, &t, &u) == MYB_OK);
  MT_CHECK(p > 1000.0 && p < 1025.0 && t > 15.0 && t < 35.0);
  /* domain */
  MT_CHECK(myb_tropo_gmf(58000.0, 0.8, 0.2, 100.0, 1.6, &mh, &mw) == MYB_E_DOMAIN);
  MT_CHECK(myb_tropo_fcula(0.0, 0.3, 100.0, 290.0, &m) == MYB_E_DOMAIN);
  MT_CHECK(myb_tropo_zenith_optical(0.3, 100.0, 1000.0, 10.0, 5.0, 375.0, &zhd, &zwd) ==
           MYB_E_DOMAIN);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_iers_cases();
  test_properties();
  return mt_finish("test_troposphere");
}
