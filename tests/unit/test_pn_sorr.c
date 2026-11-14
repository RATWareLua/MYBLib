/* SPDX-License-Identifier: BSD-3-Clause */
/* Compact binaries: the spin-orbit radiation reaction at next-to-leading order (MYB_PN_SORR5) and
 * the spin-orbit part of the dissipative tail (MYB_PN_SOTAIL), units G m = c = 1:
 *   - MYB_PN_SORR5 against the exact rational solution of its balance construction
 *     (tests/data/pn/sorr5_terms.txt);
 *   - circular orbits with aligned spins: the part of the energy loss linear in the spins at
 *     x^(3/2) and x^(5/2) relative to the quadrupole flux, against the published flux of
 *     gr-qc/0605140 in the conserved-norm spins (eq. FxK: -4 S_l - 5/4 dm Sigma_l and
 *     (-9/2 + 272/9 nu) S_l + (-13/16 + 43/4 nu) dm Sigma_l), and the angular momentum, dJ/dt =
 *     -F / Omega;
 *   - an integrated eccentric orbit with aligned spins: the energy lost over a radial period,
 *     its part linear in the spins, against the integral of the published general-orbit flux of
 *     the same paper (eqs. flux15PN, flux25PN, its spins converted to the conserved-norm ones);
 *   - the tail: on spin-orbit circular orbits the spin-linear part of the dissipative tail of
 *     MYB_PN_4RR (on the osculating Kepler orbit, -1189/24 pi S_l - 675/24 pi dm Sigma_l at x^3, as
 *     derived), and with MYB_PN_SOTAIL the published tail-induced spin-orbit flux, -16 pi S_l -
 *     31/6 pi dm Sigma_l (arXiv:1104.5659 eq. 3PNspintail). */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "myb/rel.h"
#include "pn4_test_terms.h"
#include "pn_flux_test_terms.h"
#include "rel/rel_internal.h"

#define PI 3.14159265358979323846
#define TWO_PI 6.283185307179586

static myb_pn_binary binary_of(double nu, unsigned orders) {
  const double d = sqrt(1.0 - 4.0 * nu);
  const myb_pn_binary b = {
      .gm1 = 0.5 * (1.0 + d), .gm2 = 0.5 * (1.0 - d), .c = 1.0, .orders = orders, .r0 = 1.0};
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

static double circular_v2(const myb_pn_binary* b, double r, const double* c1, const double* c2) {
  double v2 = 1.0 / r;
  for (int k = 0; k < 80; ++k) {
    const double x[3] = {r, 0.0, 0.0};
    const double v[3] = {0.0, sqrt(v2), 0.0};
    double a[3];
    MT_CHECK(myb_pn_binary_evaluate(b, x, v, c1, c2, a, NULL, NULL) == MYB_OK);
    v2 = -a[0] * r;
  }
  return v2;
}

enum { FIT = 5 };

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

/* E / mu: spin-free instantaneous terms through 4PN (r0 = r), spin-orbit through 1/c^5 */
static double energy_of(const myb_pn_binary* b, const double* x, const double* v, const double* c1,
                        const double* c2) {
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
  const double rd = dot3(n, v);
  const double v2 = dot3(v, v);
  double e = pn4_energy_4pn(nu, 1.0 / r, rd, v2, 0.0);
  e += pn4_spin_energy(3, nu, dm, 1.0 / r, rd, v2, dot3(n, sxv), dot3(n, zxv));
  e += pn4_spin_energy(5, nu, dm, 1.0 / r, rd, v2, dot3(n, sxv), dot3(n, zxv));
  return e;
}

/* dE/dv: v plus the central differences of E - E_N (step 1e-6 of |v|) */
static void denergy_dv(const myb_pn_binary* b, const double* x, const double* v, const double* c1,
                       const double* c2, double* out) {
  const double r = norm3(x);
  const double h = 1e-6 * norm3(v);
  for (int i = 0; i < 3; ++i) {
    double vp[3] = {v[0], v[1], v[2]};
    double vm[3] = {v[0], v[1], v[2]};
    vp[i] += h;
    vm[i] -= h;
    const double ep = energy_of(b, x, vp, c1, c2) - (0.5 * dot3(vp, vp) - 1.0 / r);
    const double em = energy_of(b, x, vm, c1, c2) - (0.5 * dot3(vm, vm) - 1.0 / r);
    out[i] = v[i] + (ep - em) / (2.0 * h);
  }
}

/* J_z / mu: (x x v)_z J/J_N (4PN, instantaneous) + l_S (1/c^3, 1/c^5); S / nu is constant here */
static double angular_z(const myb_pn_binary* b, const double* x, const double* v, const double* c1,
                        const double* c2) {
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
  const double rd = dot3(n, v);
  const double v2 = dot3(v, v);
  double h[3];
  cross3(x, v, h);
  double out = pn4_angular_4pn(nu, 1.0 / r, rd, v2, 0.0) * h[2];
  for (int order = 3; order <= 5; order += 2) {
    double cn = 0.0;
    double cv = 0.0;
    double cs = 0.0;
    double cz = 0.0;
    pn4_spin_angular(order, nu, dm, 1.0 / r, rd, v2, dot3(n, s), dot3(n, z), dot3(s, v), dot3(z, v),
                     &cn, &cv, &cs, &cz);
    out += cn * n[2] + cv * v[2] + cs * s[2] + cz * z[2];
  }
  return out;
}

/* ---- the code against the balance solution ------------------------------------------------- */

static void test_terms(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/pn/sorr5_terms.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  double worst = 0.0;
  int rows = 0;
  double w[24];
  for (int guard = 0; guard < 100; ++guard) {
    const int k = next_row(f, w, 24);
    if (k == 0) {
      break;
    }
    MT_CHECK(k == 16);
    const myb_pn_binary b = {.gm1 = w[0], .gm2 = 1.0 - w[0], .c = 1.0, .orders = MYB_PN_SORR5};
    double a[3];
    MT_CHECK(myb_internal_pn_binary_terms(&b, w + 1, w + 4, w + 7, w + 10, a, NULL, NULL) ==
             MYB_OK);
    const double d[3] = {a[0] - w[13], a[1] - w[14], a[2] - w[15]};
    worst = fmax(worst, norm3(d) / norm3(w + 13));
    ++rows;
  }
  fclose(f);
  printf(
      "spin-orbit radiation reaction at NLO vs the rational balance solution (%d states): "
      "%.1e\n",
      rows, worst);
  MT_CHECK(rows == 20);
  MT_CHECK(worst < 1e-13);
}

/* ---- circular orbits ---------------------------------------------------------------------------
 */

/* (P / (-(32/5) nu x^5) - 1 - f1 x) / x^(3/2) on the circular orbit of radius r with aligned
 * spins (sign times chi), P = dE/dv . a_RR, and the same for -dJ_z/dt Omega; the orbit of N +
 * 1PN + SO (LO, NLO), a_RR of `reaction` */
static void loss_factors(double nu, double r, double sign, unsigned reaction, const double* chi1,
                         const double* chi2, double* fe, double* fj, double* xout) {
  myb_pn_binary b = binary_of(nu, MYB_PN_1 | MYB_PN_SO | MYB_PN_SO25);
  const double c1[3] = {0.0, 0.0, sign * chi1[2]};
  const double c2[3] = {0.0, 0.0, sign * chi2[2]};
  const double v2 = circular_v2(&b, r, c1, c2);
  const double x[3] = {r, 0.0, 0.0};
  const double v[3] = {0.0, sqrt(v2), 0.0};
  b.orders = reaction;
  double a[3];
  MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, c1, c2, a, NULL, NULL) == MYB_OK);
  double g[3];
  b.orders = MYB_PN_1 | MYB_PN_SO | MYB_PN_SO25;
  denergy_dv(&b, x, v, c1, c2, g);
  const double p = dot3(g, a);
  const double om = sqrt(v2) / r;
  const double xv = myb_cbrt(om * om);
  const double x5 = xv * xv * xv * xv * xv;
  const double f1 = -1247.0 / 336.0 - 35.0 / 12.0 * nu;
  *fe = (p / (-32.0 / 5.0 * nu * x5) - 1.0 - f1 * xv) / (xv * sqrt(xv));
  /* dJ_z/dt = dJ_z/dv . a; on a circular orbit dJ/dt = dE/dt / Omega */
  const double hh = 1e-6 * sqrt(v2) / norm3(a);
  double vp[3];
  double vm[3];
  for (int i = 0; i < 3; ++i) {
    vp[i] = v[i] + hh * a[i];
    vm[i] = v[i] - hh * a[i];
  }
  const double dj = (angular_z(&b, x, vp, c1, c2) - angular_z(&b, x, vm, c1, c2)) / (2.0 * hh);
  *fj = (dj * om / (-32.0 / 5.0 * nu * x5) - 1.0 - f1 * xv) / (xv * sqrt(xv));
  *xout = xv;
}

static void test_circular(void) {
  static const unsigned reaction = MYB_PN_25 | MYB_PN_35 | MYB_PN_SORR | MYB_PN_SORR5;
  const double nus[2] = {0.2, 0.1};
  for (int i = 0; i < 2; ++i) {
    const double nu = nus[i];
    const double dm = sqrt(1.0 - 4.0 * nu);
    const myb_pn_binary b0 = binary_of(nu, 0);
    const double chi1[3] = {0.0, 0.0, 0.3};
    const double chi2[3] = {0.0, 0.0, 0.5};
    const double sl = b0.gm1 * b0.gm1 * chi1[2] + b0.gm2 * b0.gm2 * chi2[2];
    const double zl = b0.gm2 * chi2[2] - b0.gm1 * chi1[2];
    for (int with = 1; with >= 0; --with) {
      const unsigned rr = with ? reaction : (reaction & ~(unsigned)MYB_PN_SORR5);
      double xs[FIT];
      double qe[FIT];
      double qj[FIT];
      for (int k = 0; k < FIT; ++k) {
        const double r = 40.0 * (double)(1 << k);
        double fep, fjp, xp, fem, fjm, xm;
        loss_factors(nu, r, 1.0, rr, chi1, chi2, &fep, &fjp, &xp);
        loss_factors(nu, r, -1.0, rr, chi1, chi2, &fem, &fjm, &xm);
        xs[k] = 0.5 * (xp + xm);
        qe[k] = 0.5 * (fep - fem);
        qj[k] = 0.5 * (fjp - fjm);
      }
      double ce[FIT];
      double cj[FIT];
      poly_fit(xs, qe, ce);
      poly_fit(xs, qj, cj);
      const double lo = -4.0 * sl - 1.25 * dm * zl;
      const double nlo =
          (-4.5 + 272.0 / 9.0 * nu) * sl + (-13.0 / 16.0 + 43.0 / 4.0 * nu) * dm * zl;
      printf(
          "spin-orbit radiation reaction%s, circular orbits, nu %.1f: spin-linear energy loss "
          "x^(3/2) %.9f (FxK %.9f), x^(5/2) %.5f (FxK %.5f); angular momentum loss times "
          "Omega x^(3/2) %.9f, x^(5/2) %.5f\n",
          with ? " through NLO" : " at LO only", nu, ce[0], lo, ce[1], nlo, cj[0], cj[1]);
      /* the fit determines the x^(3/2) coefficient to about 1e-6 and the next to about 1e-3 of the
       * part of it that MYB_PN_SORR5 contributes (about 12) */
      MT_CHECK(fabs(ce[0] - lo) < 1e-5 * fabs(lo));
      MT_CHECK(fabs(cj[0] - lo) < 1e-5 * fabs(lo));
      if (with) {
        MT_CHECK(fabs(ce[1] - nlo) < 1e-2);
        MT_CHECK(fabs(cj[1] - nlo) < 1e-2);
      } else {
        MT_CHECK(fabs(ce[1] - nlo) > 5.0);
      }
    }
  }
}

/* ---- an integrated eccentric orbit
 * --------------------------------------------------------------- */

typedef struct spin_run {
  myb_pn_binary b;
} spin_run;

static myb_status spin_accel(void* context, double t, const double* x, const double* v, double* a) {
  (void)t;
  const spin_run* q = context;
  return myb_pn_binary_evaluate(&q->b, x, v, v + 3, v + 6, a, a + 3, a + 6);
}

static double rd_event(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return dot3(y, y + 9);
}

/* the published flux F (G m = c = 1) at the state, nlo: with its 2.5PN spin-orbit term */
static double flux_of(const myb_pn_binary* b, const double* y, int nlo) {
  const double x1 = b->gm1;
  const double x2 = b->gm2;
  const double nu = x1 * x2;
  const double dm = x1 - x2;
  const double* x = y;
  const double* v = y + 9;
  const double r = norm3(x);
  const double ri = 1.0 / r;
  const double n[3] = {x[0] * ri, x[1] * ri, x[2] * ri};
  const double rd = dot3(n, v);
  const double v2 = dot3(v, v);
  double s1[3];
  double s2[3];
  for (int i = 0; i < 3; ++i) {
    s1[i] = x1 * x1 * y[12 + i];
    s2[i] = x2 * x2 * y[15 + i];
  }
  /* (n, S, v), (n, Sigma, v) of S_A and of dS_A = -(1/2)(v_A.S_A) v_A + (m_B / r) S_A */
  double sc[2][2];
  for (int pass = 0; pass < 2; ++pass) {
    double a1[3];
    double a2[3];
    for (int i = 0; i < 3; ++i) {
      if (pass == 0) {
        a1[i] = s1[i];
        a2[i] = s2[i];
      } else {
        a1[i] = -0.5 * x2 * x2 * dot3(v, s1) * v[i] + x2 * ri * s1[i];
        a2[i] = -0.5 * x1 * x1 * dot3(v, s2) * v[i] + x1 * ri * s2[i];
      }
    }
    double s[3];
    double z[3];
    for (int i = 0; i < 3; ++i) {
      s[i] = a1[i] + a2[i];
      z[i] = a2[i] / x2 - a1[i] / x1;
    }
    double sxv[3];
    double zxv[3];
    cross3(s, v, sxv);
    cross3(z, v, zxv);
    sc[pass][0] = dot3(n, sxv);
    sc[pass][1] = dot3(n, zxv);
  }
  double f = pn_flux_fN(rd, v2) + pn_flux_f1PN(nu, ri, rd, v2) +
             pn_flux_fS15(dm, ri, rd, v2, sc[0][0], sc[0][1]);
  if (nlo) {
    f += pn_flux_fS25(nu, dm, ri, rd, v2, sc[0][0], sc[0][1]) -
         pn_flux_fS15(dm, ri, rd, v2, sc[1][0], sc[1][1]);
  }
  return 8.0 / 15.0 * nu * nu * ri * ri * ri * ri * f;
}

enum { NGRID = 8000 };

/* one radial period from the periapsis: the energy lost and the integrals of the published flux
 * (with and without its NLO spin term) */
static void one_period(const myb_pn_binary* b, const double y0[18], double period, double* lost,
                       double* f_nlo, double* f_lo) {
  spin_run q = {.b = *b};
  const myb_ias15_options o = {.epsilon = 1e-11, .block = 3};
  myb_ias15* d = NULL;
  MT_CHECK(myb_ias15_create(9, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_ias15_start(d, spin_accel, &q, 0.0, y0) == MYB_OK);
  const double e0 = energy_of(b, y0, y0 + 9, y0 + 12, y0 + 15);
  const double dt = 2.5 * period / (double)NGRID; /* the PN orbit is longer than Kepler's */
  double fn_prev = flux_of(b, y0, 1);
  double fl_prev = flux_of(b, y0, 0);
  double in = 0.0;
  double il = 0.0;
  double prev_rd = 0.0;
  int passed_apo = 0;
  for (int k = 1; k <= NGRID; ++k) {
    MT_CHECK(myb_ias15_integrate(d, dt * (double)k) == MYB_OK);
    const double* y = myb_ias15_state(d);
    const double now_rd = rd_event(NULL, 0.0, y);
    if (prev_rd > 0.0 && now_rd <= 0.0) {
      passed_apo = 1;
    }
    if (passed_apo && prev_rd < 0.0 && now_rd >= 0.0) {
      int found = 0;
      double te = 0.0;
      double ye[18];
      MT_CHECK(myb_ias15_locate(d, rd_event, NULL, &found, &te, ye) == MYB_OK);
      MT_CHECK(found);
      /* the last partial interval, trapezoidal */
      const double h = te - dt * (double)(k - 1);
      in += 0.5 * h * (fn_prev + flux_of(b, ye, 1));
      il += 0.5 * h * (fl_prev + flux_of(b, ye, 0));
      *lost = (e0 - energy_of(b, ye, ye + 9, ye + 12, ye + 15)) * b->gm1 * b->gm2;
      *f_nlo = in;
      *f_lo = il;
      myb_ias15_destroy(d);
      return;
    }
    const double fn = flux_of(b, y, 1);
    const double fl = flux_of(b, y, 0);
    in += 0.5 * dt * (fn + fn_prev);
    il += 0.5 * dt * (fl + fl_prev);
    fn_prev = fn;
    fl_prev = fl;
    prev_rd = now_rd;
  }
  MT_CHECK(0 && "no periapsis");
  myb_ias15_destroy(d);
}

static void test_eccentric(void) {
  /* G m = c = 1, nu = 0.2, rp = 30 and 60, ra = 2 rp, spins along the orbital normal (so that the
   * conservative remainders, functions of the state invariant about the spin axis, return to
   * their values at the next periapsis): the energy lost from one periapsis to the next, its
   * part linear in the spins ((+chi) - (-chi)) / 2, against the integral of the published flux
   * through its NLO spin-orbit term. The comparison is limited by the next order (NNLO in the
   * spins, x^2 relative), which neither side has completely: with MYB_PN_SORR5 the relative
   * difference must fall as x^2 (by 4 when the orbit is doubled), without it as x (by 2). */
  const double nu = 0.2;
  const double chi1 = 0.6;
  const double chi2 = 0.8;
  double rel[2][2];
  for (int size = 0; size < 2; ++size) {
    const double rp = 30.0 * (double)(1 << size);
    const double ra = 2.0 * rp;
    const double sma = 0.5 * (rp + ra);
    const double period = TWO_PI * sma * sqrt(sma);
    for (int with = 1; with >= 0; --with) {
      const unsigned orders = MYB_PN_1 | MYB_PN_2 | MYB_PN_SO | MYB_PN_SO25 | MYB_PN_25 |
                              MYB_PN_35 | MYB_PN_SORR | (with ? MYB_PN_SORR5 : 0u);
      const myb_pn_binary b = binary_of(nu, orders);
      double lost[2];
      double fn[2];
      double fl[2];
      for (int s = 0; s < 2; ++s) {
        const double sign = s == 0 ? 1.0 : -1.0;
        double y0[18] = {0};
        y0[0] = rp;
        y0[10] = sqrt(ra / (sma * rp));
        y0[14] = sign * chi1;
        y0[17] = sign * chi2;
        one_period(&b, y0, period, &lost[s], &fn[s], &fl[s]);
      }
      const double ll = 0.5 * (lost[0] - lost[1]);
      const double ln = 0.5 * (fn[0] - fn[1]);
      const double lo = 0.5 * (fl[0] - fl[1]);
      rel[size][with] = ll / ln - 1.0;
      printf(
          "spin-orbit radiation reaction%s, eccentric orbit rp %.0f, ra %.0f, one radial "
          "period: energy lost %.10e / %.10e (+chi / -chi), spin-linear part %.6e; published "
          "flux integral, spin-linear part %.6e (ratio %.5f), its leading spin-orbit term alone "
          "%.6e\n",
          with ? " through NLO" : " at LO only", rp, ra, lost[0], lost[1], ll, ln, ll / ln, lo);
    }
  }
  printf(
      "spin-orbit radiation reaction, eccentric orbits, orbit x 2: the relative difference "
      "from the published flux falls by %.2f with MYB_PN_SORR5 (expected 4, the NNLO remainder), "
      "by %.2f without it (2, the missing NLO term)\n",
      rel[0][1] / rel[1][1], rel[0][0] / rel[1][0]);
  MT_CHECK(fabs(rel[0][1]) < 0.3 * fabs(rel[0][0]));
  MT_CHECK(rel[0][1] / rel[1][1] > 3.0 && rel[0][1] / rel[1][1] < 6.5);
  MT_CHECK(rel[0][0] / rel[1][0] > 1.6 && rel[0][0] / rel[1][0] < 2.5);
}

/* ---- the spin-orbit tail
 * ------------------------------------------------------------------------- */

/* (P / (-(32/5) nu x^5) - 4 pi x^(3/2)) / x^3 on the spin-orbit circular orbit (N + SO), P the
 * power of the dissipative terms `reaction` with E = E_N + E_SO */
static double tail_factor(double nu, double r, double sign, unsigned reaction, const double* chi1,
                          const double* chi2, double* xout) {
  myb_pn_binary b = binary_of(nu, MYB_PN_SO);
  const double c1[3] = {0.0, 0.0, sign * chi1[2]};
  const double c2[3] = {0.0, 0.0, sign * chi2[2]};
  const double v2 = circular_v2(&b, r, c1, c2);
  const double x[3] = {r, 0.0, 0.0};
  const double v[3] = {0.0, sqrt(v2), 0.0};
  b.orders = reaction;
  double a[3];
  MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, c1, c2, a, NULL, NULL) == MYB_OK);
  const double x1 = b.gm1;
  const double x2 = b.gm2;
  const double dm = x1 - x2;
  double sz[3];
  for (int i = 0; i < 3; ++i) {
    sz[i] = x1 * x1 * c1[i] + x2 * x2 * c2[i] + dm * (x2 * c2[i] - x1 * c1[i]);
  }
  const double n[3] = {1.0, 0.0, 0.0};
  double nsz[3];
  cross3(n, sz, nsz);
  /* d/dv of e_S = -(1/r^2) [(n, S, v) + dm (n, Sigma, v)] is -(1/r^2) n x (S + dm Sigma) */
  const double p = dot3(v, a) - dot3(nsz, a) / (r * r);
  const double xv = myb_cbrt(v2 / (r * r));
  const double x5 = xv * xv * xv * xv * xv;
  const double tail = (reaction & MYB_PN_4RR) ? 4.0 * PI * xv * sqrt(xv) : 0.0;
  *xout = xv;
  return (p / (-32.0 / 5.0 * nu * x5) - tail) / (xv * xv * xv);
}

static void test_tail(void) {
  const double nu = 0.2;
  const double dm = sqrt(1.0 - 4.0 * nu);
  const myb_pn_binary b0 = binary_of(nu, 0);
  const double chi1[3] = {0.0, 0.0, 0.4};
  const double chi2[3] = {0.0, 0.0, 0.7};
  const double sl = b0.gm1 * b0.gm1 * chi1[2] + b0.gm2 * b0.gm2 * chi2[2];
  const double zl = b0.gm2 * chi2[2] - b0.gm1 * chi1[2];
  static const unsigned sets[3] = {MYB_PN_4RR, MYB_PN_4RR | MYB_PN_SOTAIL, MYB_PN_SOTAIL};
  const double expected[3] = {-1189.0 / 24.0 * PI * sl - 675.0 / 24.0 * PI * dm * zl,
                              -16.0 * PI * sl - 31.0 / 6.0 * PI * dm * zl,
                              -16.0 * PI * sl - 31.0 / 6.0 * PI * dm * zl};
  const char* names[3] = {"MYB_PN_4RR (osculating Kepler orbit)", "MYB_PN_4RR + MYB_PN_SOTAIL",
                          "MYB_PN_SOTAIL alone"};
  for (int k = 0; k < 3; ++k) {
    double xs[FIT];
    double qs[FIT];
    for (int j = 0; j < FIT; ++j) {
      const double r = 200.0 * (double)(1 << j);
      double xp = 0.0;
      double xm = 0.0;
      const double fp = tail_factor(nu, r, 1.0, sets[k], chi1, chi2, &xp);
      const double fm = tail_factor(nu, r, -1.0, sets[k], chi1, chi2, &xm);
      xs[j] = 0.5 * (xp + xm);
      qs[j] = 0.5 * (fp - fm);
    }
    double c[FIT];
    poly_fit(xs, qs, c);
    printf(
        "spin-orbit tail, circular orbits, nu 0.2, %s: spin-linear energy loss at x^3 %.6f "
        "(expected %.6f)\n",
        names[k], c[0], expected[k]);
    MT_CHECK(fabs(c[0] / expected[k] - 1.0) < 1e-4);
  }
}

static void test_domain(void) {
  const double x[3] = {10.0, 0.0, 0.0};
  const double v[3] = {0.0, 0.3, 0.0};
  double a[3];
  myb_pn_binary b = binary_of(0.25, MYB_PN_SORR5 | MYB_PN_SOTAIL);
  MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, NULL, NULL, a, NULL, NULL) == MYB_OK);
  MT_CHECK(a[0] == 0.0 && a[1] == 0.0 && a[2] == 0.0);
  /* radial motion: no orbital frequency, the tail term vanishes */
  const double vr[3] = {0.2, 0.0, 0.0};
  const double c1[3] = {0.1, 0.2, 0.3};
  b.orders = MYB_PN_SOTAIL;
  MT_CHECK(myb_internal_pn_binary_terms(&b, x, vr, c1, c1, a, NULL, NULL) == MYB_OK);
  MT_CHECK(a[0] == 0.0 && a[1] == 0.0 && a[2] == 0.0);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  (void)&pn4_energy_circular;
  test_terms();
  test_circular();
  test_eccentric();
  test_tail();
  test_domain();
  return mt_finish("test_pn_sorr");
}
