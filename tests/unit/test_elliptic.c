/* SPDX-License-Identifier: BSD-3-Clause */
/* Elliptic integrals and Jacobi functions against mpmath (60 digits, tests/data/elliptic):
 * relative errors in units of 2^-53, plus identities (Legendre's relation, sn^2 + cn^2 = 1,
 * quasi-periodicity) and the domain. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/num.h"

static const double ulp = 0x1p-53;

typedef struct worst {
  const char* name;
  double rel;
  int count;
} worst;

static void note(worst* w, double value, double expected, double scale) {
  /* relative, with an absolute floor `scale` for values crossing zero */
  const double d = fabs(value - expected) / fmax(fabs(expected), scale) / ulp;
  if (!(d <= w->rel)) {
    w->rel = isnan(d) ? (double)INFINITY : fmax(w->rel, d);
  }
  ++w->count;
}

static void test_reference(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/elliptic/reference.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  worst rf = {"R_F", 0, 0}, rd = {"R_D", 0, 0}, rj = {"R_J", 0, 0}, rc = {"R_C", 0, 0};
  worst k = {"K", 0, 0}, e = {"E", 0, 0}, pi = {"Pi", 0, 0};
  worst finc = {"F(phi)", 0, 0}, einc = {"E(phi)", 0, 0}, piinc = {"Pi(phi)", 0, 0};
  worst am = {"am", 0, 0}, sn = {"sn", 0, 0}, cn = {"cn", 0, 0}, dn = {"dn", 0, 0};
  static char line[512];
  while (fgets(line, sizeof line, f) != NULL) {
    if (line[0] == '#') {
      continue;
    }
    char name[16];
    double v[6];
    int used = 0;
    if (sscanf(line, "%15s%n", name, &used) != 1) {
      continue;
    }
    int n = 0;
    const char* p = line + used;
    while (n < 6) {
      char* end = NULL;
      const double x = strtod(p, &end);
      if (end == p) {
        break;
      }
      v[n++] = x;
      p = end;
    }
    if (strcmp(name, "RF") == 0 && n == 4) {
      note(&rf, myb_carlson_rf(v[0], v[1], v[2]), v[3], 0.0);
    } else if (strcmp(name, "RD") == 0 && n == 4) {
      note(&rd, myb_carlson_rd(v[0], v[1], v[2]), v[3], 0.0);
    } else if (strcmp(name, "RJ") == 0 && n == 5) {
      note(&rj, myb_carlson_rj(v[0], v[1], v[2], v[3]), v[4], 0.0);
    } else if (strcmp(name, "RC") == 0 && n == 3) {
      note(&rc, myb_carlson_rc(v[0], v[1]), v[2], 0.0);
    } else if (strcmp(name, "K") == 0 && n == 2) {
      note(&k, myb_elliptic_k(v[0]), v[1], 0.0);
    } else if (strcmp(name, "E") == 0 && n == 2) {
      note(&e, myb_elliptic_e(v[0]), v[1], 0.0);
    } else if (strcmp(name, "PI") == 0 && n == 3) {
      note(&pi, myb_elliptic_pi(v[0], v[1]), v[2], 0.0);
    } else if (strcmp(name, "FINC") == 0 && n == 3) {
      note(&finc, myb_elliptic_f_inc(v[0], v[1]), v[2], 1.0);
    } else if (strcmp(name, "EINC") == 0 && n == 3) {
      note(&einc, myb_elliptic_e_inc(v[0], v[1]), v[2], 1.0);
    } else if (strcmp(name, "PIINC") == 0 && n == 4) {
      note(&piinc, myb_elliptic_pi_inc(v[0], v[1], v[2]), v[3], 1.0);
    } else if (strcmp(name, "JACOBI") == 0 && n == 6) {
      double a, s, c, d;
      myb_jacobi(v[0], v[1], &a, &s, &c, &d);
      /* the amplitude and the functions carry the reduction of u: absolute errors in units of
       * |u| (conditioning) */
      const double scale = fmax(1.0, fabs(v[0]));
      note(&am, a, v[2], scale);
      note(&sn, s, v[3], scale);
      note(&cn, c, v[4], scale);
      note(&dn, d, v[5], scale);
    }
  }
  fclose(f);
  const worst* all[] = {&rf, &rd, &rj, &rc, &k, &e, &pi, &finc, &einc, &piinc, &am, &sn, &cn, &dn};
  /* bounds in units of roundoff: measured 4.4 4.0 4.5 3.0 4.6 4.2 5.0 5.6 4.3 4.4 7.4 2.8 4.9 2.6 */
  static const double bound[] = {6, 6, 6, 6, 6, 6, 7, 7, 6, 6, 9, 6, 7, 6};
  for (size_t i = 0; i < sizeof all / sizeof all[0]; ++i) {
    printf("%-8s %4d values, max error %.2f ulp\n", all[i]->name, all[i]->count, all[i]->rel);
    MT_CHECK(all[i]->count >= 150);
    MT_CHECK(all[i]->rel <= bound[i]);
  }
}

static void test_identities(void) {
  mt_rng rng = {11};
  double legendre = 0.0;
  double pythagoras = 0.0;
  double period = 0.0;
  for (int i = 0; i < 1000; ++i) {
    const double m = (double)mt_int(&rng, 1, 999999) / 1e6;
    /* Legendre: E K' + E' K - K K' = pi / 2 */
    const double k1 = myb_elliptic_k(m);
    const double k2 = myb_elliptic_k(1.0 - m);
    const double e1 = myb_elliptic_e(m);
    const double e2 = myb_elliptic_e(1.0 - m);
    legendre = fmax(legendre, fabs(e1 * k2 + e2 * k1 - k1 * k2 - 1.5707963267948966) / 1.5707963267948966);
    const double u = (double)mt_int(&rng, -100000, 100000) / 1000.0;
    double s, c, d, a;
    myb_jacobi(u, m, &a, &s, &c, &d);
    pythagoras = fmax(pythagoras, fabs(s * s + c * c - 1.0));
    pythagoras = fmax(pythagoras, fabs(d * d + m * s * s - 1.0));
    /* am(u + 2K) = am(u) + pi, and F(am(u)) = u */
    double a2;
    myb_jacobi(u + 2.0 * k1, m, &a2, NULL, NULL, NULL);
    period = fmax(period, fabs(a2 - a - 3.141592653589793) / fmax(1.0, fabs(a)));
    period = fmax(period, fabs(myb_elliptic_f_inc(a, m) - u) / fmax(1.0, fabs(u)));
  }
  printf("identities: Legendre %.2e, sn^2+cn^2 and dn^2+m sn^2 %.2e, periodicity %.2e\n",
         legendre, pythagoras, period);
  MT_CHECK(legendre < 2e-15);
  MT_CHECK(pythagoras < 5e-16);
  MT_CHECK(period < 2e-14);

  /* special and invalid arguments */
  MT_CHECK(myb_elliptic_k(0.0) == 0x1.921fb54442d18p+0);
  MT_CHECK(myb_elliptic_e(0.0) == 0x1.921fb54442d18p+0);
  MT_CHECK(myb_elliptic_e(1.0) == 1.0);
  MT_CHECK(isinf(myb_elliptic_k(1.0)));
  MT_CHECK(isnan(myb_elliptic_k(1.5)));
  MT_CHECK(isnan(myb_carlson_rf(0.0, 0.0, 1.0)));
  MT_CHECK(isnan(myb_carlson_rf(-1.0, 1.0, 1.0)));
  MT_CHECK(isnan(myb_carlson_rj(1.0, 1.0, 1.0, 0.0)));
  MT_CHECK(isnan(myb_elliptic_pi(1.0, 0.5)));
  MT_CHECK(isnan(myb_elliptic_pi_inc(2.0, 1.0, 0.5))); /* 1 - n sin^2 < 0 */
  double a, s, c, d;
  myb_jacobi(0.7, 1.0, &a, &s, &c, &d); /* m = 1: hyperbolic */
  /* against the correctly rounded functions (the C library's tanh and cosh may be off by an ulp) */
  MT_CHECK(s == myb_tanh(0.7) && c == 1.0 / myb_cosh(0.7) && c == d);
  myb_jacobi(0.7, 0.0, &a, &s, &c, &d); /* m = 0: circular */
  MT_CHECK(a == 0.7 && d == 1.0);
  myb_jacobi(0.7, 1.5, &a, &s, &c, &d);
  MT_CHECK(isnan(a) && isnan(s) && isnan(c) && isnan(d));
  /* Carlson's test values (1995, Table): R_F(1, 2, 0) = 1.3110287771461, R_C(0, 1/4) = pi */
  MT_CHECK(fabs(myb_carlson_rf(1.0, 2.0, 0.0) - 1.3110287771460599052) < 4e-16);
  MT_CHECK(fabs(myb_carlson_rc(0.0, 0.25) - 3.141592653589793) < 1e-15);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_reference();
  test_identities();
  return mt_finish("test_elliptic");
}
