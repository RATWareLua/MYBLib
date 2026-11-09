/* SPDX-License-Identifier: BSD-3-Clause */
/* WGCCRE rotation models read from NAIF's pck00011.tpc against CSPICE N0067 (sxform from J2000
 * to the body-fixed frames of 18 bodies: the Sun, planets, satellites with nutation-precession
 * terms including Mars' quadratic phase angles, an asteroid and a comet with its own epoch;
 * tests/data/wgccre/reference.txt), and the rotation helpers. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/frames.h"
#include "myb/io.h"

#define DIR MYB_TEST_DATA_DIR "/wgccre/"

static void test_reference(void) {
  myb_kpool* pool = NULL;
  myb_diag diag;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, DIR "pck00011.tpc", &diag) == MYB_OK);
  FILE* f = fopen(DIR "reference.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    myb_kpool_destroy(pool);
    return;
  }
  char line[2048];
  double worst_m = 0.0;
  double worst_dm = 0.0;
  long worst_body = 0;
  int rows = 0;
  /* bound: one iteration per line of the file */
  while (fgets(line, sizeof line, f) != NULL) {
    if (line[0] == '#') {
      continue;
    }
    long body = 0;
    double et = 0.0;
    double r[18];
    int used = 0;
    if (sscanf(line, "%ld %lf%n", &body, &et, &used) != 2) {
      continue;
    }
    const char* p = line + used;
    for (int i = 0; i < 18; ++i) {
      int n = 0;
      MT_CHECK(sscanf(p, "%lf%n", &r[i], &n) == 1);
      p += n;
    }
    myb_wgccre* model = NULL;
    const myb_status status = myb_wgccre_create(pool, body, NULL, &model, &diag);
    if (status != MYB_OK) {
      fprintf(stderr, "body %ld: %s (%s)\n", body, myb_status_name(status), diag.message);
    }
    MT_CHECK(status == MYB_OK && myb_wgccre_frame(model) == 1);
    if (model == NULL) {
      continue;
    }
    myb_rotation m;
    myb_wgccre_rotation(model, myb_time_from_seconds(et), &m);
    double e = 0.0;
    double ed = 0.0;
    double scale = 0.0;
    for (int i = 0; i < 9; ++i) {
      e = fmax(e, fabs(m.m[i] - r[i]));
      ed = fmax(ed, fabs(m.dm[i] - r[9 + i]));
      scale = fmax(scale, fabs(r[9 + i]));
    }
    ed /= scale; /* relative to the rotation rate */
    if (e > worst_m) {
      worst_m = e;
      worst_body = body;
    }
    worst_dm = fmax(worst_dm, ed);
    myb_wgccre_destroy(model);
    ++rows;
  }
  fclose(f);
  myb_kpool_destroy(pool);
  printf("vs CSPICE (%d rotations): matrix %.2e (worst body %ld), derivative %.2e relative\n",
         rows, worst_m, worst_body, worst_dm);
  MT_CHECK(rows == 126);
  /* The difference grows with the rotation rate: 1e-13 for Venus, 4e-13 for the Moon, 2.5e-11
   * for the Earth, 1.1e-10 for Vesta (w1 d ~ 6e7 degrees at +-100 years). CSPICE evaluates
   * W = w0 + w1 d in one double (one ulp of 6e7 degrees = 1.3e-10 rad); myb keeps w1 d in
   * double-word arithmetic and is exact to ~1e-15 rad (test_mpfr_frames). */
  MT_CHECK(worst_m < 2e-10 && worst_dm < 1e-9);
}

static void test_errors(void) {
  myb_kpool* pool = NULL;
  myb_diag diag;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  static const char text[] =
      "\\begindata\n"
      "BODY901_POLE_RA = ( 10 1 )\n"
      "BODY901_POLE_DEC = ( 20 )\n"
      "BODY901_PM = ( 30 360 )\n"
      "BODY901_NUT_PREC_RA = ( 1 2 )\n"
      "BODY9_NUT_PREC_ANGLES = ( 0 1 )\n"
      "BODY902_POLE_RA = ( 10 1 0 0 )\n"
      "BODY902_POLE_DEC = ( 20 )\n"
      "BODY902_PM = ( 30 360 )\n";
  MT_CHECK(myb_kpool_load_text(pool, text, strlen(text), &diag) == MYB_OK);
  myb_wgccre* model = NULL;
  MT_CHECK(myb_wgccre_create(pool, 903, NULL, &model, &diag) == MYB_E_NOT_FOUND);
  MT_CHECK(myb_wgccre_create(pool, 901, NULL, &model, &diag) == MYB_E_FORMAT);
  MT_CHECK(strstr(diag.message, "more nutation terms") != NULL);
  MT_CHECK(myb_wgccre_create(pool, 902, NULL, &model, &diag) == MYB_E_FORMAT);
  myb_kpool_destroy(pool);
  myb_wgccre_destroy(NULL);
}

static void test_rotations(void) {
  /* euler313 against its definition, compose/invert identities */
  const double angles[3] = {0.3, 1.1, -2.5};
  const double rates[3] = {1e-3, -2e-4, 7e-5};
  myb_rotation a;
  myb_rotation_euler313(angles, rates, &a);
  myb_rotation inv;
  myb_rotation_invert(&a, &inv);
  myb_rotation id;
  myb_rotation_compose(&inv, &a, &id);
  double e = 0.0;
  double ed = 0.0;
  for (int i = 0; i < 9; ++i) {
    e = fmax(e, fabs(id.m[i] - ((i % 4) == 0 ? 1.0 : 0.0)));
    ed = fmax(ed, fabs(id.dm[i]));
  }
  MT_CHECK(e < 1e-15 && ed < 1e-18);
  /* derivative against central differences of the angles */
  const double h = 1e-3;
  double ap[3];
  double am[3];
  for (int i = 0; i < 3; ++i) {
    ap[i] = angles[i] + h * rates[i];
    am[i] = angles[i] - h * rates[i];
  }
  myb_rotation p;
  myb_rotation m;
  myb_rotation_euler313(ap, rates, &p);
  myb_rotation_euler313(am, rates, &m);
  double worst = 0.0;
  for (int i = 0; i < 9; ++i) {
    worst = fmax(worst, fabs(a.dm[i] - (p.m[i] - m.m[i]) / (2.0 * h)));
  }
  MT_CHECK(worst < 1e-12);
  /* a state through a rotation and back */
  const double state[6] = {7000e3, -1200e3, 300e3, 1.5e3, 7.4e3, -0.2e3};
  double out[6];
  double back[6];
  myb_rotation_apply_state(&a, state, out);
  myb_rotation_apply_state(&inv, out, back);
  for (int i = 0; i < 6; ++i) {
    MT_CHECK(fabs(back[i] - state[i]) < 1e-8 * (i < 3 ? 1.0 : 1e-3) * 1e3);
  }
}

int main(void) {
  test_reference();
  test_errors();
  test_rotations();
  return mt_finish("test_wgccre");
}
