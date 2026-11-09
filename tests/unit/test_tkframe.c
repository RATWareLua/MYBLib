/* SPDX-License-Identifier: BSD-3-Clause */
/* NAIF TK frames: the inertial frames of data/frames/naif-inertial.tk, chained to J2000,
 * against CSPICE N0067 (pxform from J2000), and the MATRIX and QUATERNION forms. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/frames.h"
#include "myb/io.h"

/* CSPICE pxform('J2000', name, 0.0), rows */
static const struct {
  const char* name;
  int id;
  double m[9];
} expected[] = {
    {"B1950", 2, {0.9999257079523629, 0.01117893812642769, 0.0048590038414544285,
                  -0.011178938137770135, 0.9999375133499887, -2.7157926258510777e-05,
                  -0.00485900381535927, -2.716259471424704e-05, 0.9999881946023742}},
    {"FK4", 3, {0.9999256794956877, 0.011181483239171792, 0.004859003772314385,
                -0.01118148322046629, 0.9999374848933135, -2.7170293744002025e-05,
                -0.00485900381535927, -2.716259471424704e-05, 0.9999881946023742}},
    {"GALACTIC", 13, {-0.054875539395742516, -0.8734371047275961, -0.4838349917700252,
                      0.49410945362774383, -0.44482959429757496, 0.7469822486998919,
                      -0.8676661356833737, -0.19807638961301985, 0.4559837945214199}},
    {"ECLIPJ2000", 17, {1.0, 0.0, 0.0, 0.0, 0.9174820620691818, 0.3977771559319137, 0.0,
                        -0.3977771559319137, 0.9174820620691818}},
    {"ECLIPB1950", 18, {0.9999257079523629, 0.01117893812642769, 0.0048590038414544285,
                        -0.012189277138214924, 0.9173688178789828, 0.3978515722052201,
                        -9.940500920351154e-06, -0.3978812427417045, 0.9174369278459982}},
};

/* J2000 -> frame `id` by following the RELATIVE links (at most 4 here). */
static int chain(const myb_kpool* pool, int64_t id, myb_rotation* out) {
  myb_rotation total;
  myb_diag diag;
  char relative[33];
  MT_CHECK(myb_tkframe(pool, id, &total, relative, &diag) == MYB_OK);
  for (int depth = 0; depth < 4 && strcmp(relative, "J2000") != 0; ++depth) {
    char name[64];
    snprintf(name, sizeof name, "FRAME_%s", relative);
    const double* v = NULL;
    size_t n = 0;
    MT_CHECK(myb_kpool_numbers(pool, name, &v, &n) == MYB_OK && n == 1);
    if (v == NULL) {
      return 0;
    }
    myb_rotation step;
    MT_CHECK(myb_tkframe(pool, (int64_t)v[0], &step, relative, &diag) == MYB_OK);
    myb_rotation_compose(&total, &step, &total);
  }
  *out = total;
  return strcmp(relative, "J2000") == 0;
}

static void test_inertial(void) {
  myb_kpool* pool = NULL;
  myb_diag diag;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/frames/naif-inertial.tk", &diag) == MYB_OK);
  double worst = 0.0;
  for (size_t k = 0; k < sizeof expected / sizeof expected[0]; ++k) {
    myb_rotation m;
    MT_CHECK(chain(pool, expected[k].id, &m));
    for (int i = 0; i < 9; ++i) {
      worst = fmax(worst, fabs(m.m[i] - expected[k].m[i]));
    }
  }
  printf("built-in inertial frames vs CSPICE: max difference %.2e\n", worst);
  MT_CHECK(worst < 1e-15);
  myb_kpool_destroy(pool);
}

static void test_forms(void) {
  /* the same rotation as ANGLES, MATRIX (by columns, TK -> relative) and QUATERNION */
  static const char text[] =
      "\\begindata\n"
      "TKFRAME_-1_RELATIVE = 'J2000'\n"
      "TKFRAME_-1_SPEC = 'ANGLES'\n"
      "TKFRAME_-1_ANGLES = ( 30 0 0 )\n"
      "TKFRAME_-1_AXES = ( 3 1 3 )\n"
      "TKFRAME_-1_UNITS = 'DEGREES'\n"
      "TKFRAME_-2_RELATIVE = 'J2000'\n"
      "TKFRAME_-2_SPEC = 'MATRIX'\n"
      "TKFRAME_-2_MATRIX = ( 0.8660254037844387 -0.5 0  0.5 0.8660254037844387 0  0 0 1 )\n"
      "TKFRAME_-3_RELATIVE = 'J2000'\n"
      "TKFRAME_-3_SPEC = 'QUATERNION'\n"
      "TKFRAME_-3_Q = ( 0.9659258262890683 0 0 -0.25881904510252074 )\n"
      "TKFRAME_-4_RELATIVE = 'J2000'\n"
      "TKFRAME_-4_SPEC = 'MATRIX'\n"
      "TKFRAME_-4_MATRIX = ( 1 0 0 0 1 0 0 0 -1 )\n"
      "TKFRAME_-5_RELATIVE = 'J2000'\n"
      "TKFRAME_-5_SPEC = 'ANGLES'\n"
      "TKFRAME_-5_ANGLES = ( 1 2 3 )\n"
      "TKFRAME_-5_AXES = ( 3 1 3 )\n"
      "TKFRAME_-5_UNITS = 'FURLONGS'\n";
  myb_kpool* pool = NULL;
  myb_diag diag;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_text(pool, text, strlen(text), &diag) == MYB_OK);
  myb_rotation a;
  myb_rotation b;
  myb_rotation c;
  char relative[33];
  MT_CHECK(myb_tkframe(pool, -1, &a, relative, &diag) == MYB_OK && strcmp(relative, "J2000") == 0);
  MT_CHECK(myb_tkframe(pool, -2, &b, relative, &diag) == MYB_OK);
  MT_CHECK(myb_tkframe(pool, -3, &c, relative, &diag) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 9; ++i) {
    worst = fmax(worst, fmax(fabs(a.m[i] - b.m[i]), fabs(a.m[i] - c.m[i])));
  }
  MT_CHECK(worst < 1e-15);
  /* M = [30 deg]_3 maps TK to J2000, so J2000's x axis has TK coordinates (cos 30, sin 30, 0) */
  const double x[3] = {1.0, 0.0, 0.0};
  double in_tk[3];
  myb_rotation_apply_position(&a, x, in_tk);
  MT_CHECK(fabs(in_tk[0] - 0.8660254037844387) < 1e-15 && fabs(in_tk[1] - 0.5) < 1e-15);
  MT_CHECK(myb_tkframe(pool, -4, &a, relative, &diag) == MYB_E_FORMAT); /* reflection */
  MT_CHECK(myb_tkframe(pool, -5, &a, relative, &diag) == MYB_E_FORMAT); /* unit */
  MT_CHECK(myb_tkframe(pool, -6, &a, relative, &diag) == MYB_E_NOT_FOUND);
  myb_kpool_destroy(pool);
}

int main(void) {
  test_inertial();
  test_forms();
  return mt_finish("test_tkframe");
}
