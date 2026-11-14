/* SPDX-License-Identifier: BSD-3-Clause */
/* Compact binaries: the spin-spin terms at next-to-leading order (myb_pn_binary_evaluate,
 * MYB_PN_SS3; A. Bohe, G. Faye, S. Marsat, E. K. Porter, CQG 32, 195010, 2015), units
 * G m = c = 1 unless said otherwise:
 *   - the acceleration and the precession against a machine transcription of the LaTeX source
 *     of the paper (tests/data/pn/ss_terms.txt), and the leading order of rel.c (Kidder's S1 S2
 *     term and Poisson's spin-induced quadrupoles) against the same transcription;
 *   - circular orbits with aligned spins: the S^2 coefficients of the gauge-invariant energy E(x)
 *     at x^2 and x^3 (eq. Ecirc) measured from the orbits of the library's acceleration and the
 *     paper's energy;
 *   - integrated eccentric orbits with oblique spins: the part of the change of the energy that
 *     is quadratic in the spins falls as the remainder O(S^2 / c^8) with the NLO terms (by 32 when
 *     the orbit is doubled), as O(S^2 / c^6) without them (16).
 * The paper prints no angular momentum; the transcription admits one (an ansatz fitted through
 * 1/c^6 leaves a residual of 5e-16 of the terms, 1e-2 without the NLO precession), which checks
 * the precession offline; here the precession is checked against the transcription. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "myb/rel.h"
#include "pn4_test_terms.h"
#include "pn_ss_test_terms.h"
#include "rel/rel_internal.h"

#define TWO_PI 6.283185307179586
#define PI2 9.8696044010893586188

/* G m = c = 1 with the given nu */
static myb_pn_binary binary_of(double nu, unsigned orders, double r0) {
  const double d = sqrt(1.0 - 4.0 * nu);
  const myb_pn_binary b = {
      .gm1 = 0.5 * (1.0 + d), .gm2 = 0.5 * (1.0 - d), .c = 1.0, .orders = orders, .r0 = r0};
  return b;
}

static double norm3(const double* a) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

static void cross3(const double* a, const double* b, double* out) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

static double dot3(const double* a, const double* b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* The next data row of f (at most max numbers): their count, 0 at the end of the file. */
static int next_row(FILE* f, double* w, int max) {
  static char line[8192];
  for (int guard = 0; guard < 100000; ++guard) {
    if (fgets(line, sizeof line, f) == NULL) {
      return 0;
    }
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') {
      continue;
    }
    char* p = line;
    int k = 0;
    for (; k < max; ++k) {
      char* end = NULL;
      w[k] = strtod(p, &end);
      if (end == p) {
        break;
      }
      p = end;
    }
    return k;
  }
  return 0;
}

/* ---- the terms against the machine transcription ------------------------------------------- */

static void test_terms(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/pn/ss_terms.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  double worst_lo = 0.0;
  double worst_a = 0.0;
  double worst_w = 0.0;
  double worst_e = 0.0;
  int rows = 0;
  double w[40];
  for (int guard = 0; guard < 1000; ++guard) {
    const int k = next_row(f, w, 40);
    if (k == 0) {
      break;
    }
    MT_CHECK(k == 29);
    const double x1 = w[0];
    const double x2 = 1.0 - x1;
    const double* x = w + 3;
    const double* v = w + 6;
    const double* c1 = w + 9;
    const double* c2 = w + 12;
    myb_pn_binary b = {
        .gm1 = x1, .gm2 = x2, .c = 1.0, .orders = MYB_PN_SS, .kappa1 = w[1], .kappa2 = w[2]};
    double a[3];
    double d1[3];
    double d2[3];
    MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, c1, c2, a, d1, d2) == MYB_OK);
    double diff[3] = {a[0] - w[15], a[1] - w[16], a[2] - w[17]};
    worst_lo = fmax(worst_lo, norm3(diff) / norm3(w + 15));
    b.orders = MYB_PN_SS3;
    MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, c1, c2, a, d1, d2) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      diff[i] = a[i] - w[18 + i];
    }
    worst_a = fmax(worst_a, norm3(diff) / norm3(w + 18));
    /* dchi_A = Omega_A x chi_A */
    double e1[3];
    double e2[3];
    cross3(w + 21, c1, e1);
    cross3(w + 24, c2, e2);
    double f1[3] = {d1[0] - e1[0], d1[1] - e1[1], d1[2] - e1[2]};
    double f2[3] = {d2[0] - e2[0], d2[1] - e2[1], d2[2] - e2[2]};
    worst_w = fmax(worst_w, fmax(norm3(f1) / norm3(e1), norm3(f2) / norm3(e2)));
    /* the energy of the test header against the same transcription */
    const double r = norm3(x);
    const double n[3] = {x[0] / r, x[1] / r, x[2] / r};
    double s[3];
    double z[3];
    for (int i = 0; i < 3; ++i) {
      s[i] = x1 * x1 * c1[i] + x2 * x2 * c2[i];
      z[i] = x2 * c2[i] - x1 * c1[i];
    }
    const double nu = x1 * x2;
    const double dm = x1 - x2;
    const double kp = w[1] + w[2];
    const double km = w[1] - w[2];
    const double args_e4 =
        pn_ss_energy4(nu, dm, kp, km, 1.0 / r, dot3(n, v), dot3(v, v), dot3(n, s), dot3(n, z),
                      dot3(v, s), dot3(v, z), dot3(s, s), dot3(z, z), dot3(s, z));
    const double args_e6 =
        pn_ss_energy6(nu, dm, kp, km, 1.0 / r, dot3(n, v), dot3(v, v), dot3(n, s), dot3(n, z),
                      dot3(v, s), dot3(v, z), dot3(s, s), dot3(z, z), dot3(s, z));
    worst_e = fmax(worst_e,
                   fmax(fabs(args_e4 - w[27]) / fabs(w[27]), fabs(args_e6 - w[28]) / fabs(w[28])));
    ++rows;
  }
  fclose(f);
  printf(
      "spin-spin vs the machine transcription of arXiv:1501.01529 (%d states): leading order of "
      "rel.c %.1e, NLO acceleration %.1e, NLO precession %.1e; energy of the test header %.1e\n",
      rows, worst_lo, worst_a, worst_w, worst_e);
  MT_CHECK(rows == 30);
  MT_CHECK(worst_lo < 1e-13);
  MT_CHECK(worst_a < 1e-13);
  MT_CHECK(worst_w < 1e-13);
  MT_CHECK(worst_e < 1e-12);
}

/* ---- circular orbits ------------------------------------------------------------------------ */

/* v^2 of the circular orbit of radius r (x = (r, 0, 0), v = (0, v, 0)) of b with the spins */
static double circular_v2(const myb_pn_binary* b, double r, const double* chi1,
                          const double* chi2) {
  double v2 = 1.0 / r;
  for (int k = 0; k < 80; ++k) {
    const double x[3] = {r, 0.0, 0.0};
    const double v[3] = {0.0, sqrt(v2), 0.0};
    double a[3];
    MT_CHECK(myb_pn_binary_evaluate(b, x, v, chi1, chi2, a, NULL, NULL) == MYB_OK);
    v2 = -a[0] * r;
  }
  return v2;
}

enum { FIT = 5 };

/* the polynomial of degree FIT - 1 through the points (x_k, q_k): its coefficients c[k] */
static void poly_fit(const double* x, const double* q, double* c) {
  double d[FIT];
  for (int i = 0; i < FIT; ++i) {
    d[i] = q[i];
  }
  for (int j = 1; j < FIT; ++j) {
    for (int i = FIT - 1; i >= j; --i) {
      d[i] = (d[i] - d[i - 1]) / (x[i] - x[i - j]);
    }
  }
  double poly[FIT] = {0.0};
  poly[0] = d[FIT - 1];
  for (int i = FIT - 2; i >= 0; --i) {
    double next[FIT] = {0.0};
    for (int k = 0; k < FIT - 1; ++k) {
      next[k + 1] += poly[k];
      next[k] -= x[i] * poly[k];
    }
    next[0] += d[i];
    for (int k = 0; k < FIT; ++k) {
      poly[k] = next[k];
    }
  }
  for (int k = 0; k < FIT; ++k) {
    c[k] = poly[k];
  }
}

/* E / mu at the state with the spins: spin-free through 4PN (instantaneous), spin-orbit through
 * 3.5PN, spin-spin through 1/c^6 (or 1/c^4 without `nlo`) */
static double energy_of(const myb_pn_binary* b, const double* x, const double* v, const double* c1,
                        const double* c2, int nlo) {
  const double x1 = b->gm1;
  const double x2 = b->gm2;
  const double nu = x1 * x2;
  const double dm = x1 - x2;
  const double r = norm3(x);
  const double n[3] = {x[0] / r, x[1] / r, x[2] / r};
  double s[3];
  double z[3];
  for (int i = 0; i < 3; ++i) {
    s[i] = x1 * x1 * c1[i] + x2 * x2 * c2[i];
    z[i] = x2 * c2[i] - x1 * c1[i];
  }
  double sxv[3];
  double zxv[3];
  cross3(s, v, sxv);
  cross3(z, v, zxv);
  const double nsv = dot3(n, sxv);
  const double nzv = dot3(n, zxv);
  const double rd = dot3(n, v);
  const double v2 = dot3(v, v);
  double e = pn4_energy_4pn(nu, 1.0 / r, rd, v2, myb_log(r / b->r0));
  e += pn4_spin_energy(3, nu, dm, 1.0 / r, rd, v2, nsv, nzv);
  e += pn4_spin_energy(5, nu, dm, 1.0 / r, rd, v2, nsv, nzv);
  e += pn4_spin_energy(7, nu, dm, 1.0 / r, rd, v2, nsv, nzv);
  const double kp = b->kappa1 + b->kappa2;
  const double km = b->kappa1 - b->kappa2;
  e += pn_ss_energy4(nu, dm, kp, km, 1.0 / r, rd, v2, dot3(n, s), dot3(n, z), dot3(v, s),
                     dot3(v, z), dot3(s, s), dot3(z, z), dot3(s, z));
  if (nlo) {
    e += pn_ss_energy6(nu, dm, kp, km, 1.0 / r, rd, v2, dot3(n, s), dot3(n, z), dot3(v, s),
                       dot3(v, z), dot3(s, s), dot3(z, z), dot3(s, z));
  }
  return e;
}

static void test_energy_circular(void) {
  static const unsigned orbit = MYB_PN_1 | MYB_PN_2 | MYB_PN_3 | MYB_PN_SO | MYB_PN_SO25 |
                                MYB_PN_SO35 | MYB_PN_SS | MYB_PN_SS3;
  static const double kappas[2][2] = {{1.0, 1.0}, {2.5, 7.0}};
  for (int kc = 0; kc < 2; ++kc) {
    const double nu = kc == 0 ? 0.2 : 0.16;
    const double dm = sqrt(1.0 - 4.0 * nu);
    myb_pn_binary b = binary_of(nu, orbit, 1.0);
    b.kappa1 = kappas[kc][0];
    b.kappa2 = kappas[kc][1];
    const double x1 = b.gm1;
    const double x2 = b.gm2;
    const double chi1 = 0.6;
    const double chi2 = -0.45;
    double xs[FIT];
    double qs[FIT];
    for (int k = 0; k < FIT; ++k) {
      const double r = 60.0 * (double)(1 << k);
      double rr[3];
      double xx[3];
      for (int s = 0; s < 3; ++s) { /* +chi, -chi, no spin */
        const double sign = s == 0 ? 1.0 : (s == 1 ? -1.0 : 0.0);
        const double c1[3] = {0.0, 0.0, sign * chi1};
        const double c2[3] = {0.0, 0.0, sign * chi2};
        const double v2 = circular_v2(&b, r, c1, c2);
        const double x[3] = {r, 0.0, 0.0};
        const double v[3] = {0.0, sqrt(v2), 0.0};
        const double xv = myb_cbrt(v2 / (r * r));
        const double e = energy_of(&b, x, v, c1, c2, 1);
        /* the known spin-free terms through 3PN and spin-linear terms through 3.5PN (Eofx) */
        const double e1 = -0.75 - nu / 12.0;
        const double e2 = -27.0 / 8.0 + 19.0 / 8.0 * nu - nu * nu / 24.0;
        const double e3 = -675.0 / 64.0 + (34445.0 / 576.0 - 205.0 / 96.0 * PI2) * nu -
                          155.0 / 96.0 * nu * nu - 35.0 / 5184.0 * nu * nu * nu;
        const double sl = x1 * x1 * c1[2] + x2 * x2 * c2[2];
        const double zl = x2 * c2[2] - x1 * c1[2];
        const double s15 = 14.0 / 3.0 * sl + 2.0 * dm * zl;
        const double s25 = (11.0 - 61.0 / 9.0 * nu) * sl + dm * (3.0 - 10.0 / 3.0 * nu) * zl;
        const double s35 = (135.0 / 4.0 - 367.0 / 4.0 * nu + 29.0 / 12.0 * nu * nu) * sl +
                           dm * (27.0 / 4.0 - 39.0 * nu + 5.0 / 4.0 * nu * nu) * zl;
        const double h = sqrt(xv);
        rr[s] = e / (-0.5 * xv) - 1.0 - e1 * xv - e2 * xv * xv - e3 * xv * xv * xv -
                xv * h * (s15 + xv * (s25 + xv * s35));
        xx[s] = xv;
      }
      xs[k] = xx[2];
      qs[k] = (0.5 * (rr[0] + rr[1]) - rr[2]) / (xs[k] * xs[k]);
    }
    double c[FIT];
    poly_fit(xs, qs, c);
    const double sl = x1 * x1 * chi1 + x2 * x2 * chi2;
    const double zl = x2 * chi2 - x1 * chi1;
    const double kp = b.kappa1 + b.kappa2;
    const double km = b.kappa1 - b.kappa2;
    const double ref2 = pn_ss_ecirc2(nu, dm, kp, km, sl, zl);
    const double ref3 = pn_ss_ecirc3(nu, dm, kp, km, sl, zl);
    printf(
        "E(x) with aligned spins, part quadratic in the spins, nu %.2f, kappa %.1f %.1f: x^2 %.9f "
        "(Ecirc %.9f), x^3 %.6f (Ecirc %.6f)\n",
        nu, b.kappa1, b.kappa2, c[0], ref2, c[1], ref3);
    MT_CHECK(fabs(c[0] / ref2 - 1.0) < 1e-6);
    MT_CHECK(fabs(c[1] / ref3 - 1.0) < 2e-3);
  }
}

/* ---- conservation along integrated orbits ---------------------------------------------------- */

typedef struct spin_run {
  myb_pn_binary b;
} spin_run;

/* 9 "positions" (x, p1, p2) whose velocities are (v, chi1, chi2); the instantaneous 4PN terms
 * are added (without the tail) so that the spin-free energy is conserved beyond the order looked
 * at */
static myb_status spin_accel(void* context, double t, const double* x, const double* v, double* a) {
  (void)t;
  const spin_run* q = context;
  const myb_status st = myb_pn_binary_evaluate(&q->b, x, v, v + 3, v + 6, a, a + 3, a + 6);
  if (st == MYB_OK) {
    const double r = norm3(x);
    double a4 = 0.0;
    double b4 = 0.0;
    myb_internal_pn4_coefficients(q->b.gm1 * q->b.gm2, 1.0 / r, dot3(x, v) / r, dot3(v, v),
                                  myb_log(r / q->b.r0), &a4, &b4);
    for (int i = 0; i < 3; ++i) {
      a[i] -= (a4 * x[i] / r + b4 * v[i]) / (r * r);
    }
  }
  return st;
}

enum { GRID = 300 };

/* the energy (with or without the NLO spin-spin energy) at GRID + 1 times over 2 orbits */
static void energy_history(spin_run* q, const double y0[18], double period, int nlo, double* e) {
  const myb_ias15_options o = {.epsilon = 1e-10, .block = 3};
  myb_ias15* d = NULL;
  MT_CHECK(myb_ias15_create(9, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_ias15_start(d, spin_accel, q, 0.0, y0) == MYB_OK);
  int g = 0;
  for (long k = 0; k < 10000000L && g <= GRID; ++k) {
    double y[18];
    const double tg = 2.0 * period * (double)g / (double)GRID;
    if (g == 0) {
      memcpy(y, y0, sizeof y);
    } else if (tg <= myb_ias15_time(d)) {
      MT_CHECK(myb_ias15_dense(d, tg, y) == MYB_OK);
    } else {
      MT_CHECK(myb_ias15_step(d, 1e300) == MYB_OK);
      continue;
    }
    e[g] = energy_of(&q->b, y, y + 9, y + 12, y + 15, nlo);
    ++g;
  }
  MT_CHECK(g == GRID + 1);
  myb_ias15_destroy(d);
}

static void test_conservation(void) {
  /* G m = c = 1, nu = 0.2, kappa 1 and 4; an orbit of rp = 12, ra = 24 (and twice as large) in
   * the xy plane, oblique spins of about 0.6: the part of the energy change quadratic in the
   * spins, (run(+chi) + run(-chi)) / 2 - run(0), is the remainder of the spin-spin sector,
   * O(S^2 / c^6) with the leading terms only (E of order 1/r changes by r^-4: 16 when the orbit is
   * doubled) and O(S^2 / c^8) with the NLO terms (r^-5: 32). The spin-free and spin-linear
   * remainders cancel in the combination up to products of second order. */
  static double ep[GRID + 1];
  static double em[GRID + 1];
  static double e0[GRID + 1];
  const double nu = 0.2;
  const double chi1[3] = {0.3, -0.35, 0.4};
  const double chi2[3] = {-0.25, 0.4, 0.35};
  double dnlo[2];
  double dlo[2];
  for (int size = 0; size < 2; ++size) {
    const double rp = 12.0 * (double)(1 << size);
    const double ra = 2.0 * rp;
    const double sma = 0.5 * (rp + ra);
    const double vp = sqrt(ra / (sma * rp));
    const double period = TWO_PI * sma * sqrt(sma);
    for (int nlo = 1; nlo >= 0; --nlo) {
      const unsigned orders = MYB_PN_1 | MYB_PN_2 | MYB_PN_3 | MYB_PN_SO | MYB_PN_SO25 |
                              MYB_PN_SO35 | MYB_PN_SS | (nlo ? MYB_PN_SS3 : 0u);
      spin_run q = {.b = binary_of(nu, orders, 1.0)};
      q.b.kappa1 = 1.0;
      q.b.kappa2 = 4.0;
      double y0[18] = {0};
      y0[0] = rp;
      y0[10] = vp;
      for (int s = 0; s < 3; ++s) {
        const double sign = s == 0 ? 1.0 : (s == 1 ? -1.0 : 0.0);
        for (int i = 0; i < 3; ++i) {
          y0[12 + i] = sign * chi1[i];
          y0[15 + i] = sign * chi2[i];
        }
        energy_history(&q, y0, period, nlo, s == 0 ? ep : (s == 1 ? em : e0));
      }
      double worst = 0.0;
      for (int g = 1; g <= GRID; ++g) {
        const double d2 = 0.5 * ((ep[g] - ep[0]) + (em[g] - em[0])) - (e0[g] - e0[0]);
        worst = fmax(worst, fabs(d2));
      }
      if (nlo) {
        dnlo[size] = worst;
      } else {
        dlo[size] = worst;
      }
      printf(
          "spin-spin, rp %.0f, ra %.0f, %s: part of the energy change quadratic in the spins %.2e "
          "of %.3e\n",
          rp, ra, nlo ? "with the NLO terms" : "leading order only", worst, fabs(ep[0]));
    }
  }
  const double rn = dnlo[0] / dnlo[1];
  const double rl = dlo[0] / dlo[1];
  printf(
      "spin-spin, orbit x 2: the quadratic part of the energy change falls by %.1f with the NLO "
      "terms (expected 32), by %.1f with the leading order only (16)\n",
      rn, rl);
  MT_CHECK(dnlo[1] < 0.2 * dlo[1]);
  MT_CHECK(rn > 24.0 && rn < 45.0);
  MT_CHECK(rl > 12.0 && rl < 21.0);
}

/* ---- domain -----------------------------------------------------------------------------------
 */

static void test_domain(void) {
  const double x[3] = {10.0, 0.0, 0.0};
  const double v[3] = {0.0, 0.3, 0.0};
  double a[3];
  double d1[3];
  double d2[3];
  myb_pn_binary b = binary_of(0.25, MYB_PN_SS3, 1.0);
  b.kappa1 = -1.0;
  MT_CHECK(myb_pn_binary_evaluate(&b, x, v, NULL, NULL, a, NULL, NULL) == MYB_E_DOMAIN);
  b.kappa1 = 1.0;
  MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, NULL, NULL, a, d1, d2) == MYB_OK);
  MT_CHECK(a[0] == 0.0 && a[1] == 0.0 && a[2] == 0.0);
  MT_CHECK(d1[0] == 0.0 && d1[1] == 0.0 && d1[2] == 0.0);
  /* aligned spins: planar motion, no precession */
  b.orders = MYB_PN_ALL | MYB_PN_SS3;
  b.kappa2 = 3.0;
  const double xs[3] = {12.0, 5.0, 0.0};
  const double vs[3] = {-0.1, 0.22, 0.0};
  const double c1[3] = {0.0, 0.0, 0.9};
  const double c2[3] = {0.0, 0.0, -0.5};
  MT_CHECK(myb_pn_binary_evaluate(&b, xs, vs, c1, c2, a, d1, d2) == MYB_OK);
  MT_CHECK(a[2] == 0.0);
  MT_CHECK(d1[0] == 0.0 && d1[1] == 0.0 && d1[2] == 0.0);
  MT_CHECK(d2[0] == 0.0 && d2[1] == 0.0 && d2[2] == 0.0);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  /* the helpers of pn4_test_terms.h that this test does not call */
  (void)&pn4_energy_circular;
  (void)&pn4_angular_4pn;
  (void)&pn4_spin_angular;
  test_terms();
  test_energy_circular();
  test_conservation();
  test_domain();
  return mt_finish("test_pn_ss");
}
