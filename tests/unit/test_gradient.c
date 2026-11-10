/* SPDX-License-Identifier: BSD-3-Clause */
/* The gradient of the acceleration of a field of spherical harmonics
 * (myb_gravity_field_gradient):
 *  - against the second derivatives of the series as defined, differentiated by mpmath with 60
 *    digits (tests/data/grav/gradient.ref: a field of degree 8, points at all latitudes, the
 *    poles themselves, several ranges of degrees);
 *  - against central differences of the acceleration for a field of degree 60, near a pole
 *    too;
 *  - the ranges of degrees add up; symmetry and vanishing trace; errors. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/grav.h"
#include "myb/io.h"

enum { MAX_DEGREE = 60, MAX_COUNT = (MAX_DEGREE + 1) * (MAX_DEGREE + 2) / 2 };

static double largest(const double g[9]) {
  double m = 0.0;
  for (int i = 0; i < 9; ++i) {
    m = fmax(m, fabs(g[i]));
  }
  return m;
}

/* Numbers of one line after its first word. */
static int numbers(const char* line, double* values, int max) {
  int count = 0;
  const char* p = line;
  /* bound: the characters of the line */
  while (*p != '\0' && count < max) {
    while (*p == ' ') {
      ++p;
    }
    size_t used = 0;
    if (myb_parse_double(p, strlen(p), &values[count], &used) != MYB_OK) {
      break;
    }
    ++count;
    p += used;
  }
  return count;
}

static void test_reference(void) {
  FILE* file = fopen(MYB_TEST_DATA_DIR "/grav/gradient.ref", "rb");
  MT_CHECK(file != NULL);
  if (file == NULL) {
    return;
  }
  static double c[MAX_COUNT];
  static double s[MAX_COUNT];
  double gm = 0.0;
  double radius = 0.0;
  int degree = -1;
  myb_gravity_field* field = NULL;
  myb_gravity_work* work = NULL;
  char line[512];
  double worst = 0.0;
  double worst_pole = 0.0;
  int cases = 0;
  /* bound: the lines of the file */
  while (fgets(line, sizeof line, file) != NULL) {
    double v[16];
    if (strncmp(line, "field ", 6) == 0) {
      MT_CHECK(numbers(line + 6, v, 3) == 3);
      gm = v[0];
      radius = v[1];
      degree = (int)v[2];
      MT_CHECK(degree >= 0 && degree <= MAX_DEGREE);
    } else if (strncmp(line, "c ", 2) == 0) {
      MT_CHECK(numbers(line + 2, v, 4) == 4);
      const int n = (int)v[0];
      const int m = (int)v[1];
      MT_CHECK(n <= degree && m <= n);
      c[n * (n + 1) / 2 + m] = v[2];
      s[n * (n + 1) / 2 + m] = v[3];
    } else if (strncmp(line, "g ", 2) == 0) {
      if (field == NULL) {
        MT_CHECK(myb_gravity_field_create(gm, radius, degree, c, s, NULL, &field) == MYB_OK);
        MT_CHECK(myb_gravity_work_create(degree, NULL, &work) == MYB_OK);
      }
      MT_CHECK(numbers(line + 2, v, 11) == 11);
      const double r[3] = {v[2], v[3], v[4]};
      const double expected[9] = {v[5], v[6], v[7], v[6], v[8], v[9], v[7], v[9], v[10]};
      double g[9];
      MT_CHECK(myb_gravity_field_gradient(field, work, r, (int)v[0], (int)v[1], g) == MYB_OK);
      const double scale = largest(expected);
      double error = 0.0;
      for (int i = 0; i < 9; ++i) {
        error = fmax(error, fabs(g[i] - expected[i]) / scale);
      }
      if (r[0] == 0.0 && r[1] == 0.0) {
        worst_pole = fmax(worst_pole, error);
      } else {
        worst = fmax(worst, error);
      }
      MT_CHECK_SAME(g[1], g[3]);
      MT_CHECK_SAME(g[2], g[6]);
      MT_CHECK_SAME(g[5], g[7]);
      MT_CHECK(fabs(g[0] + g[4] + g[8]) <= 4e-16 * scale);
      ++cases;
    }
  }
  fclose(file);
  printf("gradient against mpmath (%d cases, degree %d): %.1e of the largest component, at the "
         "poles %.1e\n",
         cases, degree, worst, worst_pole);
  MT_CHECK(cases == 32);
  MT_CHECK(worst < 5e-15);
  MT_CHECK(worst_pole < 2e-15);

  /* the ranges of degrees add up */
  const double r[3] = {5.1e6, -4.2e6, 3.3e6};
  double all[9];
  double low[9];
  double high[9];
  MT_CHECK(myb_gravity_field_gradient(field, work, r, 0, -1, all) == MYB_OK);
  MT_CHECK(myb_gravity_field_gradient(field, work, r, 0, 2, low) == MYB_OK);
  MT_CHECK(myb_gravity_field_gradient(field, work, r, 3, -1, high) == MYB_OK);
  for (int i = 0; i < 9; ++i) {
    MT_CHECK(fabs(all[i] - (low[i] + high[i])) <= 4e-16 * largest(all));
  }
  /* the point mass alone: GM (3 r r^T - r^2) / r^5 */
  double point[9];
  MT_CHECK(myb_gravity_field_gradient(field, work, r, 0, 0, point) == MYB_OK);
  const double rr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double exact =
          gm * (3.0 * r[i] * r[j] - (i == j ? rr * rr : 0.0)) / (rr * rr * rr * rr * rr);
      MT_CHECK(fabs(point[3 * i + j] - exact) <= 2e-15 * largest(point));
    }
  }
  /* errors */
  const double origin[3] = {0.0, 0.0, 0.0};
  double g[9];
  MT_CHECK(myb_gravity_field_gradient(field, work, origin, 0, -1, g) == MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_field_gradient(field, work, r, 0, degree + 1, g) == MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_field_gradient(field, work, r, -1, 2, g) == MYB_E_DOMAIN);
  const double bad[3] = {(double)NAN, 0.0, 1.0};
  MT_CHECK(myb_gravity_field_gradient(field, work, bad, 0, -1, g) == MYB_E_DOMAIN);
  myb_gravity_work* small = NULL;
  MT_CHECK(myb_gravity_work_create(3, NULL, &small) == MYB_OK);
  MT_CHECK(myb_gravity_field_gradient(field, small, r, 0, -1, g) == MYB_E_DOMAIN);
  MT_CHECK(myb_gravity_field_gradient(field, small, r, 0, 3, g) == MYB_OK);
  myb_gravity_work_destroy(small);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
}

/* A field of degree 60 with coefficients of the size of the Earth's (Kaula's rule): the
 * gradient against central differences of the acceleration. */
static void test_differences(void) {
  static double c[MAX_COUNT];
  static double s[MAX_COUNT];
  mt_rng rng = {.state = 77};
  c[0] = 1.0;
  for (int n = 2; n <= MAX_DEGREE; ++n) {
    for (int m = 0; m <= n; ++m) {
      const double size = 1e-5 / ((double)n * (double)n);
      c[n * (n + 1) / 2 + m] = size * ((double)(mt_next(&rng) >> 11) / 4503599627370496.0 - 1.0);
      s[n * (n + 1) / 2 + m] =
          m == 0 ? 0.0 : size * ((double)(mt_next(&rng) >> 11) / 4503599627370496.0 - 1.0);
    }
  }
  myb_gravity_field* field = NULL;
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_field_create(3.986004415e14, 6378136.3, MAX_DEGREE, c, s, NULL, &field) ==
           MYB_OK);
  MT_CHECK(myb_gravity_work_create(MAX_DEGREE, NULL, &work) == MYB_OK);
  const double points[5][3] = {{6.6e6, 1.1e6, -2.0e6},
                               {-2.0e6, -6.0e6, 2.9e6},
                               {3.0e3, 4.0e3, 6.7e6}, /* 5 km from the axis */
                               {6.9e6, 0.0, 1.0},
                               {1.0e7, 2.0e7, -3.0e7}};
  double worst = 0.0;
  for (int k = 0; k < 5; ++k) {
    const double* r = points[k];
    double g[9];
    MT_CHECK(myb_gravity_field_gradient(field, work, r, 0, -1, g) == MYB_OK);
    const double h = 20.0; /* m: truncation (h / r)^2 ~ 1e-11 */
    for (int j = 0; j < 3; ++j) {
      double plus[3] = {r[0], r[1], r[2]};
      double minus[3] = {r[0], r[1], r[2]};
      plus[j] += h;
      minus[j] -= h;
      double ap[3];
      double am[3];
      MT_CHECK(myb_gravity_field_eval(field, work, plus, 0, -1, NULL, ap) == MYB_OK);
      MT_CHECK(myb_gravity_field_eval(field, work, minus, 0, -1, NULL, am) == MYB_OK);
      for (int i = 0; i < 3; ++i) {
        const double numeric = (ap[i] - am[i]) / (2.0 * h);
        worst = fmax(worst, fabs(g[3 * i + j] - numeric) / largest(g));
      }
    }
    /* the part above the point mass and the oblateness alone, whose differences are not
     * hidden by the main term */
    double high[9];
    MT_CHECK(myb_gravity_field_gradient(field, work, r, 3, -1, high) == MYB_OK);
    double worst_high = 0.0;
    for (int j = 0; j < 3; ++j) {
      double plus[3] = {r[0], r[1], r[2]};
      double minus[3] = {r[0], r[1], r[2]};
      plus[j] += h;
      minus[j] -= h;
      double ap[3];
      double am[3];
      MT_CHECK(myb_gravity_field_eval(field, work, plus, 3, -1, NULL, ap) == MYB_OK);
      MT_CHECK(myb_gravity_field_eval(field, work, minus, 3, -1, NULL, am) == MYB_OK);
      for (int i = 0; i < 3; ++i) {
        const double numeric = (ap[i] - am[i]) / (2.0 * h);
        worst_high = fmax(worst_high, fabs(high[3 * i + j] - numeric) / largest(high));
      }
    }
    MT_CHECK(worst_high < 1e-6);
  }
  printf("gradient against central differences of the acceleration (degree %d): %.1e\n",
         MAX_DEGREE, worst);
  MT_CHECK(worst < 1e-9);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
}

int main(void) {
  test_reference();
  test_differences();
  return mt_finish("test_gradient");
}
