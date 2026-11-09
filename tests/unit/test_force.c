/* SPDX-License-Identifier: BSD-3-Clause */
/* Forces of myb/force.h:
 *  - third body against the direct formula in 50-digit arithmetic (mpmath), down to |r|/|s| of
 *    1e-11 where the direct double formula loses every digit;
 *  - plates: the classical limits (absorber, mirror, diffuse 1 + 2/3) and a sphere tessellated
 *    into plates against cr = 1 + 4/9 rho_d (diffuse) and 1 (specular);
 *  - shadow fractions against the exact area of two overlapping spherical caps (uniform disk)
 *    and a 40-digit quadrature (quadratic limb darkening), mpmath: penumbra near the Earth,
 *    at LAGEOS, from the Moon's distance (annular and total eclipses);
 *  - radiation of a planet: a uniformly emitting Lambertian sphere gives M (R/r)^2 radially;
 *    convergence in the number of rings;
 *  - drag, empirical accelerations, errors. */
#include <math.h>
#include <stdio.h>

#include "mtest.h"
#include "myb/force.h"

static double rel(const double a[3], const double b[3]) {
  const double d = hypot(hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2]);
  return d / hypot(hypot(b[0], b[1]), b[2]);
}

static void test_third_body(void) {
  static const struct {
    double gm;
    double s[3];
    double r[3];
    double a[3];
  } ref[] = {
      {4902800000000.0, {380000000.0, -120000000.0, 50000000.0}, {7000000.0, 100000.0, -200000.0},
       {9.0015440530060294e-7, -4.6798215563166731e-7, 2.0758253371113761e-7}},
      {1.32712440018e+20, {140000000000.0, 50000000000.0, -20000000000.0},
       {-4000000.0, 5000000.0, 1000000.0},
       {-8.4949414408605173e-8, -2.8311355374245176e-7, -4.716001181395922e-9}},
      {1.32712440018e+20, {140000000000.0, 50000000000.0, -20000000000.0}, {1.1, -0.7, 0.3},
       {3.9689211691041806e-14, 5.7148270467562399e-14, -2.3645752276032365e-14}},
      {300000000000000.0, {10000000.0, 0, 0}, {6000000.0, 3000000.0, 0}, {6.6, -7.2, 0.0}},
  };
  double worst = 0.0;
  for (size_t k = 0; k < sizeof ref / sizeof ref[0]; ++k) {
    double a[3] = {0.0, 0.0, 0.0};
    myb_force_third_body(ref[k].gm, ref[k].s, ref[k].r, a);
    worst = fmax(worst, rel(a, ref[k].a));
  }
  printf("third body vs 50 digits (|r|/|s| from 0.6 to 1e-11): %.1e relative\n", worst);
  MT_CHECK(worst < 2e-15);
}

static void test_plates(void) {
  const double u[3] = {0.0, 0.0, -1.0}; /* light travelling down */
  const double up[3] = {0.0, 0.0, 1.0}; /* plate facing the source */
  myb_plate p = {.area = 2.0, .normal = {up[0], up[1], up[2]}};
  double a[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_radiation_plates(1e-5, 4.0, &p, 1, u, a) == MYB_OK);
  MT_CHECK(fabs(a[2] + 1e-5 * 2.0 / 4.0) < 1e-22 && a[0] == 0.0); /* absorber: P A / m along u */
  p.specular = 1.0;
  a[2] = 0.0;
  MT_CHECK(myb_force_radiation_plates(1e-5, 4.0, &p, 1, u, a) == MYB_OK);
  MT_CHECK(fabs(a[2] + 2.0 * 1e-5 * 2.0 / 4.0) < 1e-22); /* mirror: twice */
  p.specular = 0.0;
  p.diffuse = 1.0;
  a[2] = 0.0;
  MT_CHECK(myb_force_radiation_plates(1e-5, 4.0, &p, 1, u, a) == MYB_OK);
  MT_CHECK(fabs(a[2] + (1.0 + 2.0 / 3.0) * 1e-5 * 2.0 / 4.0) < 1e-21);
  p.diffuse = 0.0;
  p.thermal = 1;
  a[2] = 0.0;
  MT_CHECK(myb_force_radiation_plates(1e-5, 4.0, &p, 1, u, a) == MYB_OK);
  MT_CHECK(fabs(a[2] + (1.0 + 2.0 / 3.0) * 1e-5 * 2.0 / 4.0) < 1e-21);
  /* lit from behind: nothing */
  p.normal[2] = -1.0;
  a[2] = 0.0;
  MT_CHECK(myb_force_radiation_plates(1e-5, 4.0, &p, 1, u, a) == MYB_OK);
  MT_CHECK(a[2] == 0.0);
  /* a sphere of plates: cr = 1 + 4/9 rho_d for diffuse reflection, 1 for specular */
#define NLAT 400
#define NLON 800
  static myb_plate sphere[NLAT * NLON];
  const double pi = 3.141592653589793;
  for (int mode = 0; mode < 2; ++mode) {
    double cross_section = 0.0;
    for (int i = 0; i < NLAT; ++i) {
      const double t0 = pi * i / NLAT;
      const double t1 = pi * (i + 1) / NLAT;
      const double tm = 0.5 * (t0 + t1);
      for (int j = 0; j < NLON; ++j) {
        const double ph = 2.0 * pi * (j + 0.5) / NLON;
        myb_plate* q = &sphere[i * NLON + j];
        q->area = (cos(t0) - cos(t1)) * 2.0 * pi / NLON; /* unit sphere */
        q->normal[0] = sin(tm) * cos(ph);
        q->normal[1] = sin(tm) * sin(ph);
        q->normal[2] = cos(tm);
        const double nn = sqrt(q->normal[0] * q->normal[0] + q->normal[1] * q->normal[1] +
                               q->normal[2] * q->normal[2]);
        for (int k = 0; k < 3; ++k) {
          q->normal[k] /= nn;
        }
        q->specular = mode == 1 ? 0.6 : 0.0;
        q->diffuse = mode == 0 ? 0.6 : 0.0;
        q->thermal = 0;
        if (q->normal[2] > 0.0) {
          cross_section += q->area * q->normal[2];
        }
      }
    }
    double s[3] = {0.0, 0.0, 0.0};
    MT_CHECK(myb_force_radiation_plates(1.0, 1.0, sphere, (size_t)NLAT * NLON, u, s) == MYB_OK);
    const double cr = -s[2] / cross_section;
    const double expected = mode == 0 ? 1.0 + 4.0 / 9.0 * 0.6 : 1.0;
    printf("sphere of %d plates, %s 0.6: cr %.6f (expected %.6f)\n", NLAT * NLON,
           mode == 0 ? "diffuse" : "specular", cr, expected);
    MT_CHECK(fabs(cr - expected) < 2e-5 && fabs(s[0]) < 1e-12 && fabs(s[1]) < 1e-12);
  }
  /* sphere model */
  double b[3] = {0.0, 0.0, 0.0};
  myb_force_radiation_sphere(4.56e-6, 1.13, 6.95e-4, u, b);
  MT_CHECK(fabs(b[2] + 4.56e-6 * 1.13 * 6.95e-4) < 1e-24);
  /* errors */
  p = (myb_plate){.area = 1.0, .normal = {0.0, 0.0, 1.0}, .specular = 0.7, .diffuse = 0.4};
  MT_CHECK(myb_force_radiation_plates(1e-5, 1.0, &p, 1, u, a) == MYB_E_DOMAIN);
  p = (myb_plate){.area = 1.0, .normal = {0.0, 0.1, 1.0}};
  MT_CHECK(myb_force_radiation_plates(1e-5, 1.0, &p, 1, u, a) == MYB_E_DOMAIN);
  p = (myb_plate){.area = -1.0, .normal = {0.0, 0.0, 1.0}};
  MT_CHECK(myb_force_radiation_plates(1e-5, 1.0, &p, 1, u, a) == MYB_E_DOMAIN);
  p = (myb_plate){.area = 1.0, .normal = {0.0, 0.0, 1.0}};
  MT_CHECK(myb_force_radiation_plates(1e-5, 0.0, &p, 1, u, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_radiation_plates(1e-5, 1.0, NULL, 0, u, a) == MYB_OK);
}

static void test_shadow(void) {
  const double au = 149597870700.0;
  const double rs = 695700000.0;
  const double source[3] = {au, 0.0, 0.0};
  const double limb[2] = {0.93, -0.23};
  static const struct {
    double o[3];
    double radius;
    double uniform;
    double darkened;
  } ref[] = {
      {{2342994.220787208, 6383602.280950603, 0.0}, 6378137.0, 0.80461692005220465,
       0.83236517056547392},
      {{2372655.6182774897, 6372637.234070073, 0.0}, 6378137.0, 0.19562003699712403,
       0.16783288473447793},
      {{10455210.333411586, 6421952.731383469, 0.0}, 6378137.0, 0.98137369348986361,
       0.98876743943900355},
      {{10510073.314036904, 6331765.862187998, 0.0}, 6378137.0, 0.0066848621333324279,
       0.0036324437841829946},
      {{383952792.9852993, 6021026.389903262, 0.0}, 6378137.0, 0.40275901455244469,
       0.38442845695419421},
      {{383931158.1226081, 7270888.709976264, 0.0}, 6378137.0, 0.82199884793492595,
       0.84689164330288182},
      {{383999999.52, 19199.999992, 0.0}, 1737400.0, 0.053449307640285286,
       0.027971356076620286},
      {{383999808.000016, 383999.9360000032, 0.0}, 6378137.0, 0.0, 0.0},
  };
  double worst = 0.0;
  for (size_t k = 0; k < sizeof ref / sizeof ref[0]; ++k) {
    const double fu = myb_shadow_fraction(source, rs, ref[k].o, ref[k].radius, NULL);
    const double fl = myb_shadow_fraction(source, rs, ref[k].o, ref[k].radius, limb);
    worst = fmax(worst, fmax(fabs(fu - ref[k].uniform), fabs(fl - ref[k].darkened)));
  }
  printf("shadow fractions (penumbra at 420 km, at LAGEOS, from the Moon, annular): %.1e\n",
         worst);
  MT_CHECK(worst < 2e-13);
  /* full light, umbra, inside the occulter, occulter behind */
  MT_CHECK(myb_shadow_fraction(source, rs, (const double[3]){-7e6, 0.0, 0.0}, 6378137.0, NULL) ==
           1.0);
  MT_CHECK(myb_shadow_fraction(source, rs, (const double[3]){7e6, 0.0, 0.0}, 6378137.0, NULL) ==
           0.0);
  MT_CHECK(myb_shadow_fraction(source, rs, (const double[3]){5e6, 0.0, 0.0}, 6378137.0, NULL) ==
           0.0);
  MT_CHECK(myb_shadow_fraction(source, rs, (const double[3]){2.0 * au, 0.0, 0.0}, 6378137.0,
                               NULL) == 1.0);
  /* continuity across the penumbra: monotonic in the separation */
  double last = 0.0;
  int monotonic = 1;
  for (int k = 0; k <= 200; ++k) {
    const double y = 6.33e6 + 1e3 * k;
    const double o[3] = {sqrt(1.227e7 * 1.227e7 - y * y), y, 0.0};
    const double f = myb_shadow_fraction(source, rs, o, 6378137.0, limb);
    if (f < last) {
      monotonic = 0;
    }
    last = f;
  }
  MT_CHECK(monotonic && last == 1.0);
}

static void test_planet_radiation(void) {
  const myb_planet_radiation ir_only = {
      .radius = 6378137.0, .albedo = {0.0}, .emissivity = {1.0, 0.0, 0.0, 0.0, 0.0}};
  const double r[3] = {3.0e6, -4.0e6, 1.1e7}; /* about LAGEOS */
  const double sun[3] = {1.4e11, 5.0e10, 1.0e10};
  const double flux = 1361.0;
  double prev = 0.0;
  double err = 0.0;
  for (int rings = 10; rings <= 40; rings += 30) {
    double a[3] = {0.0, 0.0, 0.0};
    MT_CHECK(myb_force_planet_radiation(&ir_only, r, sun, flux, 0.0, 1.0, rings, a) == MYB_OK);
    const double rr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    const double m = 0.25 * flux * (6378137.0 / rr) * (6378137.0 / rr) / 299792458.0;
    const double e[3] = {m * r[0] / rr, m * r[1] / rr, m * r[2] / rr};
    err = rel(a, e);
    printf("uniform Lambertian emitter, %d rings: %.1e relative\n", rings, err);
    if (rings == 40) {
      MT_CHECK(err < 0.5 * prev);
    }
    prev = err;
  }
  MT_CHECK(err < 3e-4); /* midpoint rule in the rings: about 0.4 / rings^2 */
  /* Knocke's Earth: albedo and IR at LAGEOS, about 1e-10 m/s^2 for cr A/m = 1.13 * 6.95e-4 */
  const myb_planet_radiation knocke = {.radius = 6378137.0,
                                       .albedo = {0.34, 0.0, 0.10, 0.0, 0.29},
                                       .emissivity = {0.68, 0.0, -0.07, 0.0, -0.18},
                                       .omega = 2.0 * 3.141592653589793 / (365.25 * 86400.0)};
  double a[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_planet_radiation(&knocke, r, sun, flux, 1.0, 1.13 * 6.95e-4, 20, a) ==
           MYB_OK);
  const double na = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  printf("Knocke albedo + IR at LAGEOS: %.3e m/s^2\n", na);
  MT_CHECK(na > 3e-11 && na < 3e-10);
  MT_CHECK(myb_force_planet_radiation(&knocke, (const double[3]){1e6, 0.0, 0.0}, sun, flux, 0.0,
                                      1.0, 20, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_planet_radiation(&knocke, r, sun, flux, 0.0, 1.0, 0, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_planet_radiation(&knocke, r, sun, flux, 0.0, 1.0, 51, a) == MYB_E_DOMAIN);
}

static void test_drag_empirical(void) {
  double a[3] = {0.0, 0.0, 0.0};
  myb_force_drag(1e-12, 2.2 * 0.01, (const double[3]){7500.0, 0.0, 0.0}, a);
  MT_CHECK(fabs(a[0] / (-0.5 * 1e-12 * 0.022 * 7500.0 * 7500.0) - 1.0) < 1e-15 && a[1] == 0.0);
  const double r[3] = {7e6, 0.0, 0.0};
  const double v[3] = {0.0, 7.5e3, 1.0e3};
  double e[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_empirical_rtn(r, v, (const double[3]){1.0, 2.0, 3.0}, e) == MYB_OK);
  /* R = x, N = (r x v)/|r x v| = (0, -1e3, 7.5e3)/n, T = N x R */
  const double n = hypot(1e3, 7.5e3);
  const double expected[3] = {1.0, 2.0 * 7.5e3 / n - 3.0 * 1e3 / n, 2.0 * 1e3 / n + 3.0 * 7.5e3 / n};
  MT_CHECK(rel(e, expected) < 1e-15);
  MT_CHECK(myb_force_empirical_rtn(r, (const double[3]){1.0, 0.0, 0.0}, (const double[3]){1, 1, 1},
                                   e) == MYB_E_DOMAIN);
  MT_CHECK(fabs(myb_radiation_flux(3.828e26, 149597870700.0) - 1361.2) < 0.1);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_third_body();
  test_plates();
  test_shadow();
  test_planet_radiation();
  test_drag_empirical();
  return mt_finish("test_force");
}
