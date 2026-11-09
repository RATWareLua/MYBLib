/* SPDX-License-Identifier: BSD-3-Clause */
/* The Earth rotation angle of myb_cip_era against MPFR (256 bits) evaluating the same
 * expression, 2 pi (c0 + f + c1 (day + f)) mod 2 pi with the constants of
 * data/frames/iau2006-cip.tk, over +-300 years: the error must stay at the rounding of the
 * final angle (SOFA's eraEra00 sums ~200 turns in one double and is off by up to 3e-13 rad). */
#include <math.h>
#include <mpfr.h>
#include <stdio.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/frames.h"
#include "myb/io.h"

enum { P = 256 };

int main(void) {
  myb_kpool* pool = NULL;
  myb_cip* cip = NULL;
  myb_diag diag;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/frames/iau2006-cip.tk", &diag) == MYB_OK);
  MT_CHECK(myb_cip_create(pool, "IAU2006A", NULL, &cip, &diag) == MYB_OK);
  const double* era = NULL;
  size_t n = 0;
  MT_CHECK(myb_kpool_numbers(pool, "IAU2006A_ERA", &era, &n) == MYB_OK && n == 2);
  if (cip == NULL || era == NULL) {
    return mt_finish("test_mpfr_frames");
  }
  mpfr_t c0;
  mpfr_t c1;
  mpfr_t x;
  mpfr_t y;
  mpfr_t pi2;
  mpfr_inits2(P, c0, c1, x, y, pi2, (mpfr_ptr)0);
  mpfr_set_d(c0, era[0], MPFR_RNDN);
  mpfr_set_d(c1, era[1], MPFR_RNDN);
  mpfr_const_pi(pi2, MPFR_RNDN);
  mpfr_mul_ui(pi2, pi2, 2, MPFR_RNDN);
  mt_rng rng = {42};
  double worst = 0.0;
  for (int k = 0; k < 20000; ++k) {
    const int64_t day = mt_int(&rng, -110000, 110000);
    const double seconds = (double)mt_int(&rng, 0, 86399) + (double)mt_int(&rng, 0, 1 << 20) /
                                                                (double)(1 << 20);
    const myb_time t = myb_time_add((myb_time){day * 86400, 0.0}, seconds);
    const double angle = myb_cip_era(cip, t, NULL);
    /* exact fraction of the day of t, then turns = c0 + f + c1 (day + f) */
    mpfr_set_si(x, (long)(t.s - day * 86400), MPFR_RNDN);
    mpfr_add_d(x, x, t.f, MPFR_RNDN);
    mpfr_div_ui(x, x, 86400, MPFR_RNDN); /* f */
    mpfr_add_si(y, x, (long)day, MPFR_RNDN);
    mpfr_mul(y, y, c1, MPFR_RNDN);
    mpfr_add(y, y, x, MPFR_RNDN);
    mpfr_add(y, y, c0, MPFR_RNDN);
    mpfr_frac(y, y, MPFR_RNDN);
    if (mpfr_sgn(y) < 0) {
      mpfr_add_ui(y, y, 1, MPFR_RNDN);
    }
    mpfr_mul(y, y, pi2, MPFR_RNDN);
    mpfr_sub_d(y, y, angle, MPFR_RNDN);
    double error = fabs(mpfr_get_d(y, MPFR_RNDN));
    error = fmin(error, fabs(error - 6.283185307179586));
    worst = fmax(worst, error);
  }
  printf("ERA: max |error| %.2e rad over 20000 instants, years 1700-2300\n", worst);
  MT_CHECK(worst < 2e-15);

  /* Prime meridian W = w0 + w1 d of a fast rotator (Vesta's constants of pck00011) */
  static const char vesta[] =
      "\\begindata\n"
      "BODY2000004_POLE_RA = ( 309.031 0 0 )\n"
      "BODY2000004_POLE_DEC = ( 42.235 0 0 )\n"
      "BODY2000004_PM = ( 285.39 1617.3329428 0 )\n";
  myb_kpool* text = NULL;
  myb_wgccre* model = NULL;
  MT_CHECK(myb_kpool_create(NULL, &text) == MYB_OK);
  MT_CHECK(myb_kpool_load_text(text, vesta, sizeof vesta - 1, &diag) == MYB_OK);
  MT_CHECK(myb_wgccre_create(text, 2000004, NULL, &model, &diag) == MYB_OK);
  mpfr_t w1;
  mpfr_init2(w1, P);
  mpfr_set_d(w1, 1617.3329428, MPFR_RNDN);
  double worst_w = 0.0;
  for (int k = 0; k < 5000 && model != NULL; ++k) {
    const myb_time t = {(int64_t)mt_int(&rng, -1000000, 1000000) * 3000, 0.375};
    double angles[3];
    double rates[3];
    myb_wgccre_euler(model, t, angles, rates);
    mpfr_set_d(x, (double)t.s, MPFR_RNDN); /* exact; long is 32 bits on Windows */
    mpfr_add_d(x, x, t.f, MPFR_RNDN);
    mpfr_div_ui(x, x, 86400, MPFR_RNDN); /* d */
    mpfr_mul(x, x, w1, MPFR_RNDN);
    mpfr_add_d(x, x, 285.39, MPFR_RNDN);
    mpfr_div_ui(x, x, 360, MPFR_RNDN);
    mpfr_frac(x, x, MPFR_RNDN);
    if (mpfr_sgn(x) < 0) {
      mpfr_add_ui(x, x, 1, MPFR_RNDN);
    }
    mpfr_mul(x, x, pi2, MPFR_RNDN);
    mpfr_sub_d(x, x, angles[2], MPFR_RNDN);
    double error = fabs(mpfr_get_d(x, MPFR_RNDN));
    error = fmin(error, fabs(error - 6.283185307179586));
    worst_w = fmax(worst_w, error);
  }
  printf("W of a fast rotator: max |error| %.2e rad over 5000 instants, +-95 years\n", worst_w);
  MT_CHECK(worst_w < 3e-15);
  mpfr_clear(w1);
  myb_wgccre_destroy(model);
  myb_kpool_destroy(text);
  mpfr_clears(c0, c1, x, y, pi2, (mpfr_ptr)0);
  myb_cip_destroy(cip);
  myb_kpool_destroy(pool);
  return mt_finish("test_mpfr_frames");
}
