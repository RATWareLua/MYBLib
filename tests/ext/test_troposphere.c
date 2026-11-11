/* SPDX-License-Identifier: BSD-3-Clause */
/* ext/troposphere against the test cases published in the headers of the IERS Conventions
 * (2010) chapter 9 routines (FCUL_ZD_HPA.F, FCUL_A.F, VMF1.F, VMF1_HT.F, GMF.F, GPT.F; the
 * authors' code, printed to 19 digits), and properties of the models. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

/* VMF3 and GPT3 against the routines of TU Wien (gpt3_5.f90, vmf3_ht.f90, asknewet.f run
 * apart: tests/data/tropo/gpt3_vmf3.ref) on the grid gpt3_5.grd, if downloaded: 20 sites and
 * dates — between nodes, at a node, across the first meridian and the date line, within half
 * a step of the poles, at sea level and on Everest, a leap year after February, the means
 * without the seasons. */
static void test_gpt3(void) {
  FILE* file = fopen(MYB_DATA_DIR "/download/gpt3_5.grd", "rb");
  if (file == NULL) {
    printf("GPT3: data/download/gpt3_5.grd is not there (cmake -P cmake/FetchData.cmake): "
           "skipped\n");
    return;
  }
  enum { CAPACITY = 2 * 1024 * 1024 };
  char* text = malloc(CAPACITY);
  MT_CHECK(text != NULL);
  const size_t size = text != NULL ? fread(text, 1, CAPACITY, file) : 0;
  fclose(file);
  myb_gpt3* grid = NULL;
  myb_diag diag = {0};
  MT_CHECK(size > 0 && size < CAPACITY);
  const myb_status loaded = myb_gpt3_load(text, size, NULL, &grid, &diag);
  if (loaded != MYB_OK) {
    printf("GPT3: status %d, line %ld: %s%c", (int)loaded, diag.line, diag.message, 10);
    free(text);
    MT_CHECK(loaded == MYB_OK);
    return;
  }
  MT_CHECK(fabs(myb_gpt3_step(grid) - 5.0 * DEG) < 1e-15);
  FILE* reference = fopen(MYB_TEST_DATA_DIR "/tropo/gpt3_vmf3.ref", "r");
  MT_CHECK(reference != NULL);
  char line[1024];
  int cases = 0;
  double worst[3] = {0.0, 0.0, 0.0}; /* GPT3, VMF3, the wet delay */
  while (reference != NULL && fgets(line, sizeof line, reference) != NULL) {
    if (line[0] == '#') {
      continue;
    }
    double in[6];
    double r[16];
    int read = 0;
    const char* at = line;
    for (int k = 0; k < 22; ++k) {
      char* end = NULL;
      const double value = strtod(at, &end);
      if (end == at) {
        break;
      }
      if (k < 6) {
        in[k] = value;
      } else {
        r[k - 6] = value;
      }
      at = end;
      ++read;
    }
    MT_CHECK(read == 22);
    if (read != 22) {
      continue;
    }
    myb_gpt3_values v;
    MT_CHECK(myb_gpt3_at(grid, in[0], in[1] * DEG, in[2] * DEG, in[3], (int)in[4], &v) ==
             MYB_OK);
    const double mine[13] = {v.pressure, v.temperature, v.lapse_rate, v.mean_temperature,
                             v.vapour_pressure, v.ah, v.aw, v.lambda, v.undulation, v.gn_h,
                             v.ge_h, v.gn_w, v.ge_w};
    for (int k = 0; k < 13; ++k) {
      /* by what the quantity is worth: the temperature in degrees Celsius passes zero */
      const double scale = k == 1 ? 273.15 : fmax(fabs(r[k]), 1e-6);
      worst[0] = fmax(worst[0], fabs(mine[k] - r[k]) / scale);
    }
    double mh = 0.0;
    double mw = 0.0;
    double zwd = 0.0;
    MT_CHECK(myb_tropo_vmf3(r[5], r[6], in[0], in[1] * DEG, in[2] * DEG, in[3], in[5], &mh,
                            &mw) == MYB_OK);
    MT_CHECK(myb_tropo_zwd_askne(r[4], r[3], r[7], &zwd) == MYB_OK);
    worst[1] = fmax(worst[1], fmax(rel(mh, r[13]), rel(mw, r[14])));
    worst[2] = fmax(worst[2], rel(zwd, r[15]));
    ++cases;
  }
  if (reference != NULL) {
    fclose(reference);
  }
  printf("GPT3 on its 5 degree grid, %d sites and dates against the routines of TU Wien: "
         "%.1e of the values; VMF3 %.1e; the wet delay of Askne and Nordius %.1e\n",
         cases, worst[0], worst[1], worst[2]);
  MT_CHECK(cases == 20);
  /* the sums of the seasons cancel in the small quantities: the last digits of cos and sin
   * of the two libraries show there */
  MT_CHECK(worst[0] < 3e-12 && worst[1] < 2e-14 && worst[2] < 1e-14);

  /* the gradients of GPT3 as a delay: each part with its own mapping function */
  myb_gpt3_values v;
  MT_CHECK(myb_gpt3_at(grid, 59000.0, 48.2 * DEG, 16.37 * DEG, 156.0, 0, &v) == MYB_OK);
  double delay = 0.0;
  const double e = 7.0 * DEG;
  const double az = 0.6;
  MT_CHECK(myb_tropo_gradient_delay_parts(v.gn_h, v.ge_h, v.gn_w, v.ge_w, e, az, &delay) ==
           MYB_OK);
  const double tangent = myb_sin(e) * myb_tan(e);
  const double expected =
      (v.gn_h * myb_cos(az) + v.ge_h * myb_sin(az)) / (tangent + 0.0031) +
      (v.gn_w * myb_cos(az) + v.ge_w * myb_sin(az)) / (tangent + 0.0007);
  printf("  Vienna: %.1f hPa, %.1f C, e %.1f hPa; the gradients at 7 degrees of elevation: "
         "%.1f mm\n",
         v.pressure, v.temperature, v.vapour_pressure, 1e3 * delay);
  MT_CHECK(fabs(delay - expected) < 1e-17 && fabs(delay) > 1e-4 && fabs(delay) < 0.05);
  /* a longitude a turn away is the same point */
  myb_gpt3_values turned;
  MT_CHECK(myb_gpt3_at(grid, 59000.0, 48.2 * DEG, (16.37 - 360.0) * DEG, 156.0, 0, &turned) ==
           MYB_OK);
  MT_CHECK(fabs(turned.pressure - v.pressure) < 1e-9);
  /* refused */
  MT_CHECK(myb_gpt3_at(grid, 59000.0, 1.6, 0.0, 0.0, 0, &v) == MYB_E_DOMAIN);
  MT_CHECK(myb_gpt3_at(grid, (double)NAN, 0.5, 0.0, 0.0, 0, &v) == MYB_E_DOMAIN);
  double mh = 0.0;
  double mw = 0.0;
  MT_CHECK(myb_tropo_vmf3(0.0012, 0.0006, 59000.0, 0.8, 0.2, 100.0, 1.6, &mh, &mw) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_tropo_zwd_askne(10.0, 270.0, -1.0, &delay) == MYB_E_DOMAIN);
  MT_CHECK(myb_tropo_gradient_delay_parts(1e-4, 1e-4, 1e-4, 1e-4, 0.0, 0.0, &delay) ==
           MYB_E_DOMAIN);
  myb_gpt3_destroy(grid);
  /* grids that are refused: a torn line, a count that is no grid, nodes out of order */
  myb_gpt3* none = NULL;
  size_t second = 0;
  for (size_t i = 0, lines = 0; i < size && lines < 2; ++i) {
    if (text[i] == '\n') {
      ++lines;
      second = i + 1;
    }
  }
  MT_CHECK(myb_gpt3_load(text, second, NULL, &none, &diag) == MYB_E_FORMAT);
  text[second + 3] = '9'; /* the latitude of the second node */
  MT_CHECK(myb_gpt3_load(text, size, NULL, &none, &diag) == MYB_E_FORMAT && diag.line == 3);
  MT_CHECK(myb_gpt3_load("% nothing\n", 10, NULL, &none, &diag) == MYB_E_FORMAT);
  myb_gpt3_destroy(NULL);
  free(text);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_iers_cases();
  test_properties();
  test_gpt3();
  return mt_finish("test_troposphere");
}
