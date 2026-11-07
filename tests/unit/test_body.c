/* SPDX-License-Identifier: BSD-3-Clause */
/* Parameters with uncertainty, GM from mass, parameters from kernels, body validation. */
#include <math.h>
#include <string.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/body.h"

static void test_product(void) {
  myb_param p = {};
  MT_CHECK(myb_param_product((myb_param){3.0, 0.3}, (myb_param){4.0, 0.4}, &p) == MYB_OK);
  MT_CHECK(p.value == 12.0);
  MT_CHECK(fabs(p.sigma - sqrt(1.2 * 1.2 + 1.2 * 1.2)) < 1e-15);
  /* Huge and tiny magnitudes: no overflow in the uncertainty. */
  MT_CHECK(myb_param_product((myb_param){1e200, 1e199}, (myb_param){1e-100, 0.0}, &p) == MYB_OK);
  MT_CHECK(fabs(p.sigma / 1e99 - 1.0) < 1e-15);
  MT_CHECK(myb_param_product((myb_param){1.0, -1.0}, (myb_param){1.0, 0.0}, &p) == MYB_E_DOMAIN);
  MT_CHECK(myb_param_product((myb_param){(double)NAN, 0.0}, (myb_param){1.0, 0.0}, &p) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_param_product((myb_param){1e300, 0.0}, (myb_param){1e300, 0.0}, &p) == MYB_E_DOMAIN);
}

static void test_gm(void) {
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/constants/codata-2022.tk", NULL) == MYB_OK);
  myb_param g = {};
  MT_CHECK(myb_param_from_kpool(pool, "MYB_CONST_G", "MYB_CONST_G_SIGMA", 1.0, &g) == MYB_OK);
  MT_CHECK(g.value == 6.67430e-11 && g.sigma == 0.00015e-11);
  myb_param c = {};
  MT_CHECK(myb_param_from_kpool(pool, "MYB_CONST_C", NULL, 1.0, &c) == MYB_OK);
  MT_CHECK(c.value == 299792458.0 && c.sigma == 0.0);

  /* A 1000 kg spacecraft: GM inherits the 2.2e-5 relative uncertainty of G. */
  myb_param gm = {};
  MT_CHECK(myb_gm_from_mass((myb_param){1000.0, 0.0}, g, &gm) == MYB_OK);
  MT_CHECK(fabs(gm.value - 6.67430e-8) < 1e-22);
  MT_CHECK(fabs(gm.sigma / gm.value - 0.00015 / 6.67430) < 1e-15);
  MT_CHECK(myb_gm_from_mass((myb_param){-1.0, 0.0}, g, &gm) == MYB_E_DOMAIN);

  /* Unit conversion while reading: km^3 s^-2 -> m^3 s^-2. */
  static const char text[] =
      "\\begindata\nBODY999_GM = ( 3.98600435507D+05 )\nBODY999_GM_SIGMA = 0.0003\n"
      "TWO = ( 1 2 )\nWORD = 'x'\n";
  MT_CHECK(myb_kpool_load_text(pool, text, strlen(text), NULL) == MYB_OK);
  MT_CHECK(myb_param_from_kpool(pool, "BODY999_GM", "BODY999_GM_SIGMA", 1e9, &gm) == MYB_OK);
  MT_CHECK(gm.value == 3.98600435507e5 * 1e9 && gm.sigma == 0.0003 * 1e9);
  MT_CHECK(myb_param_from_kpool(pool, "MISSING", NULL, 1.0, &gm) == MYB_E_NOT_FOUND);
  MT_CHECK(myb_param_from_kpool(pool, "TWO", NULL, 1.0, &gm) == MYB_E_FORMAT);
  MT_CHECK(myb_param_from_kpool(pool, "WORD", NULL, 1.0, &gm) == MYB_E_FORMAT);
  myb_kpool_destroy(pool);
}

static void test_validate(void) {
  myb_body body = {.id = 399, .gm = {3.986004418e14, 8e5}, .motion = MYB_MOTION_DYNAMIC};
  MT_CHECK(myb_body_validate(&body) == MYB_OK);
  body.gm.value = 0.0; /* a test particle */
  MT_CHECK(myb_body_validate(&body) == MYB_OK);
  body.gm.value = -1.0;
  MT_CHECK(myb_body_validate(&body) == MYB_E_DOMAIN);
  body.gm = (myb_param){1.0, -0.1};
  MT_CHECK(myb_body_validate(&body) == MYB_E_DOMAIN);
  body.gm = (myb_param){(double)INFINITY, 0.0};
  MT_CHECK(myb_body_validate(&body) == MYB_E_DOMAIN);
  body.gm = (myb_param){1.0, 0.0};
  body.motion = (myb_motion)7;
  MT_CHECK(myb_body_validate(&body) == MYB_E_DOMAIN);
}

int main(void) {
  test_product();
  test_gm();
  test_validate();
  return mt_finish("test_body");
}
