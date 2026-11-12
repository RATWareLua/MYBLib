/* SPDX-License-Identifier: BSD-3-Clause */
/* Compact binaries at 4PN with the tail, the higher spin-orbit terms and the radiation
 * reaction of the tail and of the spins (myb_pn_binary_evaluate, MYB_PN_4 ... MYB_PN_SORR),
 * units G m = c = 1 unless said otherwise:
 *   - the coefficients against machine transcriptions of the LaTeX sources
 *     (tests/data/pn/pn4_coefficients.txt, so_terms.txt) and the tail against an mpmath
 *     evaluation of its Fourier series (tail_terms.txt); the tail on a circular orbit against
 *     its closed forms (arXiv:1711.00283 eqs. omtail and eomcirc);
 *   - 4PN physics: the generalised Kepler law through 4PN with the tail logarithm ln(16 gamma)
 *     and the gauge logarithms (eq. keplerlaw); the gauge-invariant energy of circular orbits
 *     E(x) through 4PN with -3969/128 and (448/15) nu ln(16 x) (eq. Ecirc), from the
 *     instantaneous energy of the same paper and the tail energy; the periastron advance
 *     K = Omega_phi / Omega_r through 4PN with ln x, ln 2, ln 3 (eq. Kcirc), measured on
 *     integrated near-circular orbits, and its independence of the gauge length r0;
 *   - 4PN radiation reaction: the factor 1 + 4 pi x^(3/2) of the damping on a circular orbit,
 *     and on an orbit of e = 0.5 the energy lost in a period against the tail flux
 *     (4 pi / 5) G^2 M n^7 sum_p p^7 |I_p|^2 / c^8 from an independent Fourier analysis;
 *   - spin-orbit: the spin-linear coefficients of E(x) at 1.5, 2.5 and 3.5PN (arXiv:1212.5520
 *     eq. Eofx) and of the precession frequency of circular orbits (eq. precessionCO); the
 *     conservation of the energy and of the total angular momentum (its component along the
 *     orbital normal) with the spins through NNLO along integrated eccentric precessing orbits,
 *     to the expected remainder O(S / c^9);
 *   - spin-orbit radiation reaction: on circular orbits with aligned spins the energy loss has
 *     Kidder's spin-orbit factor 1 + x^(3/2) (-4 S_l - 5/4 dm Sigma_l) exactly; the energy lost
 *     by an integrated quasi-circular orbit against the flux formula. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "myb/kepler.h"
#include "myb/rel.h"
#include "pn4_test_terms.h"
#include "rel/rel_internal.h"

#define PI 3.14159265358979323846
#define TWO_PI 6.283185307179586
#define PI2 9.8696044010893586188
#define EULER 0.57721566490153286061
#define LN2 0.69314718055994530942
#define LN3 1.09861228866810969140

static const unsigned conservative4 = MYB_PN_1 | MYB_PN_2 | MYB_PN_3 | MYB_PN_4;
static const unsigned spins_nnlo =
    MYB_PN_1 | MYB_PN_2 | MYB_PN_3 | MYB_PN_SO | MYB_PN_SO25 | MYB_PN_SO35;

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

/* v^2 of the circular orbit of radius r (x = (r, 0, 0), v = (0, v, 0)) of b with the spins
 * (NULL: none): Omega^2 r = -a_x by fixed-point iteration. */
static double circular_v2(const myb_pn_binary* b, double r, const double* chi1,
                          const double* chi2) {
  double v2 = 1.0 / r;
  for (int k = 0; k < 60; ++k) {
    const double x[3] = {r, 0.0, 0.0};
    const double v[3] = {0.0, sqrt(v2), 0.0};
    double a[3];
    MT_CHECK(myb_pn_binary_evaluate(b, x, v, chi1, chi2, a, NULL, NULL) == MYB_OK);
    v2 = -a[0] * r;
  }
  return v2;
}

/* ---- the coefficients against the machine transcriptions ------------------------------------ */

static void test_coefficients(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/pn/pn4_coefficients.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  static const unsigned flags[3] = {MYB_PN_1, MYB_PN_2, MYB_PN_3};
  double worst[4] = {0.0};
  int rows = 0;
  double w[16];
  for (int guard = 0; guard < 1000; ++guard) {
    const int k = next_row(f, w, 16);
    if (k == 0) {
      break;
    }
    MT_CHECK(k == 13);
    const double nu = w[0];
    const double r = w[1];
    const double rd = w[3];
    const double vt = w[4];
    const double x[3] = {r, 0.0, 0.0};
    const double v[3] = {rd, vt, 0.0};
    for (int o = 0; o < 3; ++o) {
      myb_pn_binary b = binary_of(nu, flags[o], w[2]);
      double a[3];
      MT_CHECK(myb_pn_binary_evaluate(&b, x, v, NULL, NULL, a, NULL, NULL) == MYB_OK);
      const double bb = -a[1] * r * r / vt;
      const double aa = -a[0] * r * r - 1.0 - bb * rd;
      worst[o] = fmax(worst[o], fmax(fabs(aa - w[5 + 2 * o]) / (1.0 + fabs(w[5 + 2 * o])),
                                     fabs(bb - w[6 + 2 * o]) / (1.0 + fabs(w[6 + 2 * o]))));
    }
    double a4 = 0.0;
    double b4 = 0.0;
    myb_internal_pn4_coefficients(nu, 1.0 / r, rd, rd * rd + vt * vt, myb_log(r / w[2]), &a4, &b4);
    worst[3] = fmax(worst[3], fmax(fabs(a4 - w[11]) / (1.0 + fabs(w[11])),
                                   fabs(b4 - w[12]) / (1.0 + fabs(w[12]))));
    ++rows;
  }
  fclose(f);
  printf(
      "A, B vs the machine transcription of arXiv:1711.00283 (%d states): 1PN %.1e, 2PN "
      "%.1e, 3PN %.1e (the code of rel.c, transcribed from Blanchet's review), 4PN %.1e\n",
      rows, worst[0], worst[1], worst[2], worst[3]);
  MT_CHECK(rows == 40);
  for (int o = 0; o < 4; ++o) {
    MT_CHECK(worst[o] < 1e-14);
  }
}

static void test_spin_terms(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/pn/so_terms.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  static const unsigned flags[4] = {MYB_PN_SO, MYB_PN_SO25, MYB_PN_SO35, MYB_PN_SORR};
  static const int column[4] = {13, 16, 19, 28};
  double worst_a[4] = {0.0};
  double worst_w[3] = {0.0};
  int rows = 0;
  double w[40];
  for (int guard = 0; guard < 1000; ++guard) {
    const int k = next_row(f, w, 40);
    if (k == 0) {
      break;
    }
    MT_CHECK(k == 31);
    const double x1 = w[0];
    const double* x = w + 1;
    const double* v = w + 4;
    const double* c1 = w + 7;
    const double* c2 = w + 10;
    const double r = norm3(x);
    double nxv[3];
    cross3(x, v, nxv);
    for (int i = 0; i < 3; ++i) {
      nxv[i] /= r;
    }
    for (int o = 0; o < 4; ++o) {
      const myb_pn_binary b = {.gm1 = x1, .gm2 = 1.0 - x1, .c = 1.0, .orders = flags[o]};
      double a[3];
      double d1[3];
      double d2[3];
      MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, c1, c2, a, d1, d2) == MYB_OK);
      const double* ref = w + column[o];
      double diff[3] = {a[0] - ref[0], a[1] - ref[1], a[2] - ref[2]};
      worst_a[o] = fmax(worst_a[o], norm3(diff) / norm3(ref));
      if (o < 3) {
        /* dchi_A = alpha_A (n x v) x chi_A, alpha at the order of the flag */
        double e1[3];
        double e2[3];
        double o1[3];
        double o2[3];
        for (int i = 0; i < 3; ++i) {
          o1[i] = w[22 + o] * nxv[i];
          o2[i] = w[25 + o] * nxv[i];
        }
        cross3(o1, c1, e1);
        cross3(o2, c2, e2);
        double f1[3] = {d1[0] - e1[0], d1[1] - e1[1], d1[2] - e1[2]};
        double f2[3] = {d2[0] - e2[0], d2[1] - e2[1], d2[2] - e2[2]};
        worst_w[o] = fmax(worst_w[o], fmax(norm3(f1) / norm3(e1), norm3(f2) / norm3(e2)));
      }
    }
    ++rows;
  }
  fclose(f);
  printf(
      "spin-orbit vs the machine transcription of arXiv:1212.5520 (%d states): acceleration "
      "LO %.1e (Kidder's code of rel.c), NLO %.1e, NNLO %.1e; precession LO %.1e, 2PN %.1e, "
      "3PN %.1e; Will's radiation reaction %.1e\n",
      rows, worst_a[0], worst_a[1], worst_a[2], worst_w[0], worst_w[1], worst_w[2], worst_a[3]);
  MT_CHECK(rows == 30);
  for (int o = 0; o < 4; ++o) {
    MT_CHECK(worst_a[o] < 1e-13);
  }
  for (int o = 0; o < 3; ++o) {
    MT_CHECK(worst_w[o] < 1e-13);
  }
}

/* the tail part of MYB_PN_4: the terms minus the instantaneous A_4PN, B_4PN */
static void tail_part(const myb_pn_binary* b, const double* x, const double* v, double* out) {
  double a[3];
  MT_CHECK(myb_internal_pn_binary_terms(b, x, v, NULL, NULL, a, NULL, NULL) == MYB_OK);
  const double r = norm3(x);
  const double rd = dot3(x, v) / r;
  double a4 = 0.0;
  double b4 = 0.0;
  myb_internal_pn4_coefficients(b->gm1 * b->gm2, 1.0 / r, rd, dot3(v, v), myb_log(r / b->r0), &a4,
                                &b4);
  for (int i = 0; i < 3; ++i) {
    out[i] = a[i] + (a4 * x[i] / r + b4 * v[i]) / (r * r);
  }
}

static void test_tail_reference(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/pn/tail_terms.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  int rows = 0;
  double w[20];
  for (int guard = 0; guard < 100; ++guard) {
    const int k = next_row(f, w, 20);
    if (k == 0) {
      break;
    }
    MT_CHECK(k == 14);
    const double nu = w[0];
    const double* x = w + 1;
    const double* v = w + 4;
    myb_pn_binary b = binary_of(nu, MYB_PN_4, 1.0);
    double t[3];
    tail_part(&b, x, v, t);
    double dc[3] = {t[0] - w[7], t[1] - w[8], t[2] - w[9]};
    b.orders = MYB_PN_4RR;
    double ad[3];
    MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, NULL, NULL, ad, NULL, NULL) == MYB_OK);
    double dd[3] = {ad[0] - w[10], ad[1] - w[11], ad[2] - w[12]};
    const double ec = norm3(dc) / norm3(w + 7);
    const double ed = norm3(dd) / norm3(w + 10);
    printf(
        "tail vs its Fourier series with mpmath, e %.2f, r %.1f: conservative %.1e, "
        "dissipative %.1e\n",
        w[13], norm3(x), ec, ed);
    MT_CHECK(ec < 1e-9 && ed < 1e-9);
    ++rows;
  }
  fclose(f);
  MT_CHECK(rows == 7);
}

static void test_tail_circular(void) {
  const double nus[2] = {0.25, 0.1};
  const double radii[3] = {20.0, 300.0, 5000.0};
  double worst_c = 0.0;
  double worst_d = 0.0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double nu = nus[i];
      const double r = radii[j];
      const double g = 1.0 / r;
      /* a Newtonian circular orbit, inclined */
      const double x[3] = {0.6 * r, 0.0, 0.8 * r};
      const double vt = sqrt(g);
      const double v[3] = {0.0, vt, 0.0};
      myb_pn_binary b = binary_of(nu, MYB_PN_4, 3.0);
      double t[3];
      tail_part(&b, x, v, t);
      /* eq. (omtail): a = -omega^2 x, omega^2 = (128/5) nu gamma^4 [ln(16 gamma) + 2 gE + 1/2] /
       * r^3 */
      const double om2 =
          128.0 / 5.0 * nu * g * g * g * g * (myb_log(16.0 * g) + 2.0 * EULER + 0.5) * g * g * g;
      b.orders = MYB_PN_4RR;
      double ad[3];
      MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, NULL, NULL, ad, NULL, NULL) == MYB_OK);
      /* eq. (eomcirc): -(32/5) nu gamma^4 4 pi gamma^(3/2) v */
      const double kd = -32.0 / 5.0 * nu * g * g * g * g * 4.0 * PI * g * sqrt(g);
      for (int k = 0; k < 3; ++k) {
        worst_c = fmax(worst_c, fabs(t[k] + om2 * x[k]) / (om2 * r));
        worst_d = fmax(worst_d, fabs(ad[k] - kd * v[k]) / fabs(kd * vt));
      }
    }
  }
  printf(
      "tail on circular orbits vs the closed forms: conservative (omtail) %.1e, dissipative "
      "4 pi x^(3/2) (eomcirc) %.1e\n",
      worst_c, worst_d);
  MT_CHECK(worst_c < 1e-12);
  MT_CHECK(worst_d < 1e-12);
}

/* ---- 4PN circular orbits ----------------------------------------------------------------------
 */

static void test_kepler_law(void) {
  const double nus[2] = {0.25, 0.1};
  const double logs[2] = {0.0, 2.0};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      const double nu = nus[i];
      const double lr = logs[j];
      double d[3];
      double expected0 = 0.0;
      for (int k = 0; k < 3; ++k) {
        const double g = 5e-4 * (double)(1 << k);
        const double r = 1.0 / g;
        myb_pn_binary b = binary_of(nu, conservative4, r * myb_exp(-lr));
        const double om2r3 = circular_v2(&b, r, NULL, NULL) * r;
        const double c1 = -3.0 + nu;
        const double c2 = 6.0 + 41.0 / 4.0 * nu + nu * nu;
        const double c3 = -10.0 + (-75707.0 / 840.0 + 41.0 / 64.0 * PI2 + 22.0 * lr) * nu +
                          9.5 * nu * nu + nu * nu * nu;
        /* eq. (keplerlaw) with r'_0 = r''_0 = r0 */
        const double c4 =
            15.0 +
            nu * (19644217.0 / 33600.0 + 163.0 / 1024.0 * PI2 + 256.0 / 5.0 * EULER +
                  128.0 / 5.0 * myb_log(16.0 * g) + (82.0 - 372.0) * lr) +
            nu * nu * (44329.0 / 336.0 - 1907.0 / 64.0 * PI2 + (-992.0 / 3.0 + 720.0) * lr) +
            51.0 / 4.0 * nu * nu * nu + nu * nu * nu * nu;
        const double res = (om2r3 - 1.0 - c1 * g - c2 * g * g - c3 * g * g * g) / (g * g * g * g);
        d[k] = res - c4;
        if (k == 0) {
          expected0 = c4;
        }
      }
      const double d0 = (8.0 * d[0] - 6.0 * d[1] + d[2]) / 3.0;
      printf(
          "Kepler law through 4PN, nu %.2f, ln(r/r0) %.0f: 4PN coefficient at gamma = 5e-4 "
          "%.4f (eq. keplerlaw %.4f), the difference extrapolated to gamma = 0 %.1e\n",
          nu, lr, expected0 + d[0], expected0, d0);
      MT_CHECK(fabs(d0) < 0.05); /* of a coefficient of 12 to 128, its tail part ~ -100 nu */
      MT_CHECK(fabs(d[0]) < 0.2);
    }
  }
}

static void test_energy_circular(void) {
  const double nus[2] = {0.25, 0.1};
  for (int i = 0; i < 2; ++i) {
    const double nu = nus[i];
    double e4g[2];
    for (int gauge = 0; gauge < 2; ++gauge) { /* r0 = r and r0 = r / 1000 */
      const double lr = gauge == 0 ? 0.0 : myb_log(1000.0);
      double d[3];
      double dn[3];
      double expected0 = 0.0;
      for (int k = 0; k < 3; ++k) {
        const double g = 5e-4 * (double)(1 << k);
        const double r = 1.0 / g;
        myb_pn_binary b = binary_of(nu, conservative4, r * myb_exp(-lr));
        const double v2 = circular_v2(&b, r, NULL, NULL);
        const double x = myb_cbrt(v2 / (r * r)); /* (v / r)^(2/3) */
        /* the tail energy of a circular orbit, E0tail + EDC of eqs. (E0tail), (EDC) (EAC = 0):
         * -(I3)^2 [ln(16 gamma) + 2 gE + 2/5] with (I3)^2 = 32 nu^2 gamma^5 */
        const double e_tail =
            -32.0 * nu * g * g * g * g * g * (myb_log(16.0 * g) + 2.0 * EULER + 0.4);
        const double e_inst = pn4_energy_circular(nu, g, v2, lr);
        const double e1 = -0.75 - nu / 12.0;
        const double e2 = -27.0 / 8.0 + 19.0 / 8.0 * nu - nu * nu / 24.0;
        const double e3 = -675.0 / 64.0 + (34445.0 / 576.0 - 205.0 / 96.0 * PI2) * nu -
                          155.0 / 96.0 * nu * nu - 35.0 / 5184.0 * nu * nu * nu;
        /* eq. (Ecirc) */
        const double e4 = -3969.0 / 128.0 +
                          (-123671.0 / 5760.0 + 9037.0 / 1536.0 * PI2 + 896.0 / 15.0 * EULER +
                           448.0 / 15.0 * myb_log(16.0 * x)) *
                              nu +
                          (-498449.0 / 3456.0 + 3157.0 / 576.0 * PI2) * nu * nu +
                          301.0 / 1728.0 * nu * nu * nu + 77.0 / 31104.0 * nu * nu * nu * nu;
        const double x4 = x * x * x * x;
        const double res =
            ((e_inst + e_tail) / (-0.5 * x) - 1.0 - e1 * x - e2 * x * x - e3 * x * x * x) / x4;
        const double res_n =
            (e_inst / (-0.5 * x) - 1.0 - e1 * x - e2 * x * x - e3 * x * x * x) / x4;
        d[k] = res - e4;
        dn[k] = res_n - e4;
        if (k == 0) {
          expected0 = e4;
        }
      }
      const double d0 = (8.0 * d[0] - 6.0 * d[1] + d[2]) / 3.0;
      e4g[gauge] = expected0 + d0;
      printf(
          "E(x) through 4PN, nu %.2f, r0 = r%s: 4PN coefficient at x ~ 5e-4 %.5f (eq. Ecirc "
          "%.5f), difference extrapolated to x = 0 %.1e; without the tail energy %.3f\n",
          nu, gauge == 0 ? "" : " / 1000", expected0 + d[0], expected0, d0, expected0 + dn[0]);
      MT_CHECK(fabs(d0) < 3e-3 * fabs(expected0));
      MT_CHECK(fabs(dn[0] - d[0]) > 0.1 * fabs(expected0)); /* the tail energy matters */
    }
    /* r0 changes E(gamma) by 93 nu ln 1000 gamma^5 (eq. Ecircgam) but not E(x) */
    printf(
        "E(x), nu %.2f: extrapolated 4PN coefficients with r0 = r and r0 = r / 1000 differ by "
        "%.1e\n",
        nu, e4g[0] - e4g[1]);
    MT_CHECK(fabs(e4g[0] - e4g[1]) < 1e-3 * fabs(e4g[0]));
  }
}

/* ---- the periastron advance at 4PN ------------------------------------------------------------
 */

typedef struct run {
  myb_pn_binary b;
} run;

static myb_status orbit_accel(void* context, double t, const double* x, const double* v,
                              double* a) {
  (void)t;
  const run* q = context;
  return myb_pn_binary_evaluate(&q->b, x, v, NULL, NULL, a, NULL, NULL);
}

static double rdot_event(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return dot3(y, y + 3);
}

/* K over `periods` radial periods of a near-circular orbit; x from Omega_phi */
static void measure_k(const myb_pn_binary* b, double r, int periods, double* k_out, double* x_out) {
  run q = {.b = *b};
  const double vt = sqrt(circular_v2(b, r, NULL, NULL)) * (1.0 + 1.5e-5);
  const double y0[6] = {r, 0.0, 0.0, 0.0, vt, 0.0};
  const myb_ias15_options o = {.epsilon = 1e-9};
  myb_ias15* d = NULL;
  MT_CHECK(myb_ias15_create(3, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_ias15_start(d, orbit_accel, &q, 0.0, y0) == MYB_OK);
  double phi = 0.0;
  double last = 0.0;
  double prev = rdot_event(NULL, 0.0, y0);
  double t_first = 0.0;
  double phi_first = 0.0;
  double t_last = 0.0;
  double phi_last = 0.0;
  int found_count = 0;
  for (long k = 0; k < 10000000L && found_count <= periods; ++k) {
    MT_CHECK(myb_ias15_step(d, 1e300) == MYB_OK);
    const double* y = myb_ias15_state(d);
    const double now = rdot_event(NULL, 0.0, y);
    if (prev < 0.0 && now >= 0.0) {
      int found = 0;
      double te = 0.0;
      double ye[6];
      MT_CHECK(myb_ias15_locate(d, rdot_event, NULL, &found, &te, ye) == MYB_OK);
      if (found) {
        const double a = myb_atan2(ye[1], ye[0]);
        double da = a - last;
        da -= TWO_PI * floor(da / TWO_PI + 0.5);
        if (found_count == 0) {
          t_first = te;
          phi_first = phi + da;
        }
        t_last = te;
        phi_last = phi + da;
        ++found_count;
      }
    }
    const double a = myb_atan2(y[1], y[0]);
    double da = a - last;
    da -= TWO_PI * floor(da / TWO_PI + 0.5);
    phi += da;
    last = a;
    prev = now;
  }
  myb_ias15_destroy(d);
  const double dphi = phi_last - phi_first;
  *k_out = dphi / (TWO_PI * periods);
  *x_out = myb_cbrt((dphi / (t_last - t_first)) * (dphi / (t_last - t_first)));
}

static void test_periastron_advance_4pn(void) {
  const double nus[2] = {0.25, 0.1};
  for (int i = 0; i < 2; ++i) {
    const double nu = nus[i];
    const double k3 = 67.5 + (-649.0 / 4.0 + 123.0 / 32.0 * PI2) * nu + 7.0 * nu * nu;
    double k4g[2];
    for (int gauge = 0; gauge < 2; ++gauge) {
      double dlt[3];
      double expected0 = 0.0;
      double tail0 = 0.0;
      for (int s = 0; s < 3; ++s) {
        const double r = 400.0 / (double)(1 << s);
        myb_pn_binary b = binary_of(nu, conservative4, gauge == 0 ? r : r / 1000.0);
        double kk = 0.0;
        double x = 0.0;
        measure_k(&b, r, 8, &kk, &x);
        const double res =
            (kk - 1.0 - 3.0 * x - (13.5 - 7.0 * nu) * x * x - k3 * x * x * x) / (x * x * x * x);
        /* eq. (Kcirc) */
        const double lx = myb_log(x);
        const double k4 = 2835.0 / 8.0 +
                          (-275941.0 / 360.0 + 48007.0 / 3072.0 * PI2 - 1256.0 / 15.0 * lx -
                           592.0 / 15.0 * LN2 - 1458.0 / 5.0 * LN3 - 2512.0 / 15.0 * EULER) *
                              nu +
                          (5861.0 / 12.0 - 451.0 / 32.0 * PI2) * nu * nu -
                          98.0 / 27.0 * nu * nu * nu;
        dlt[s] = res - k4;
        if (s == 0) {
          expected0 = k4;
          /* eq. (Ktailcirc) */
          tail0 = (352.0 / 5.0 - 1256.0 / 15.0 * lx - 592.0 / 15.0 * LN2 - 1458.0 / 5.0 * LN3 -
                   2512.0 / 15.0 * EULER) *
                  nu;
        }
      }
      /* dlt(x) = O(x) + O(x^2): Richardson on x, 2x, 4x */
      const double d0 = (8.0 * dlt[0] - 6.0 * dlt[1] + dlt[2]) / 3.0;
      k4g[gauge] = expected0 + d0;
      printf(
          "K through 4PN, nu %.2f, r0 = r%s: (K - K_3PN) / x^4 - k4 = %.3f (x 0.0025), %.3f "
          "(x 0.005), %.3f (x 0.01); extrapolated %.4f of k4 %.3f (its tail part %.3f)\n",
          nu, gauge == 0 ? "" : " / 1000", dlt[0], dlt[1], dlt[2], d0, expected0, tail0);
      MT_CHECK(fabs(d0) < 0.01 * fabs(expected0));
      MT_CHECK(fabs(d0) < 0.1 * fabs(tail0));
    }
    MT_CHECK(fabs(k4g[0] - k4g[1]) < 0.01 * fabs(k4g[0]));
  }
}

/* ---- the radiation reaction of the tail -------------------------------------------------------
 */

static void test_tail_reaction_circular(void) {
  const double nus[2] = {0.25, 0.1};
  for (int i = 0; i < 2; ++i) {
    const double nu = nus[i];
    const double g = 1e-5;
    const double r = 1.0 / g;
    myb_pn_binary b = binary_of(nu, MYB_PN_1 | MYB_PN_2 | MYB_PN_3, r);
    const double vt = sqrt(circular_v2(&b, r, NULL, NULL));
    const double x[3] = {r, 0.0, 0.0};
    const double v[3] = {0.0, vt, 0.0};
    b.orders = MYB_PN_4RR;
    double a[3];
    MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, NULL, NULL, a, NULL, NULL) == MYB_OK);
    /* a = -(32/5) nu g^4 [... + 4 pi g^(3/2)] v (eq. eomcirc) */
    const double coefficient = a[1] / vt / (-32.0 / 5.0 * nu * g * g * g * g * g * sqrt(g));
    printf(
        "4PN radiation reaction on a circular orbit, nu %.2f, gamma 1e-5: coefficient of "
        "x^(3/2) %.6f (4 pi = %.6f)\n",
        nu, coefficient, 4.0 * PI);
    MT_CHECK(fabs(coefficient / (4.0 * PI) - 1.0) < 3e-4);
    MT_CHECK(fabs(a[0]) < 1e-3 * fabs(a[1]));
  }
}

static double radial_event(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return dot3(y, y + 3);
}

/* sum_{p>=1} p |J_p|^2 and sum_{p>=1} |J_p|^2 of the Fourier coefficients J_p (in the mean
 * anomaly) of I3_ij = (2 nu / r^2) [-4 n^<i v^j> + 3 rdot n^<i n^j>] along the Kepler orbit
 * (a, e) with G m = 1 (arXiv:1711.00283, footnote to eq. acctailCM): samples uniform in the
 * mean anomaly, Kepler's equation solved by myb_kepler_elliptic. */
static void tail_flux_sums(double nu, double a, double e, double* sum_p, double* sum_1) {
  enum { M = 4096, P = 600 };

  static double samples[M][4]; /* xx, xy, yy, zz */
  const double n = 1.0 / (a * sqrt(a));
  const double b = a * sqrt(1.0 - e * e);
  for (int k = 0; k < M; ++k) {
    const double l = TWO_PI * (double)k / (double)M;
    double ea = 0.0;
    MT_CHECK(myb_kepler_elliptic(l, e, &ea) == MYB_OK);
    const double ce = myb_cos(ea);
    const double se = myb_sin(ea);
    const double w = 1.0 - e * ce;
    const double pos[3] = {a * (ce - e), b * se, 0.0};
    const double vel[3] = {-a * n * se / w, b * n * ce / w, 0.0};
    const double r = norm3(pos);
    const double nn[3] = {pos[0] / r, pos[1] / r, 0.0};
    const double rd = dot3(nn, vel);
    const double nv = rd; /* n.v */
    double t[3][3];
    for (int p = 0; p < 3; ++p) {
      for (int q = 0; q < 3; ++q) {
        const double delta = p == q ? 1.0 : 0.0;
        const double sym_nv = 0.5 * (nn[p] * vel[q] + nn[q] * vel[p]) - delta * nv / 3.0;
        const double sym_nn = nn[p] * nn[q] - delta / 3.0;
        t[p][q] = 2.0 * nu / (r * r) * (-4.0 * sym_nv + 3.0 * rd * sym_nn);
      }
    }
    samples[k][0] = t[0][0];
    samples[k][1] = t[0][1];
    samples[k][2] = t[1][1];
    samples[k][3] = t[2][2];
  }
  double s_p = 0.0;
  double s_1 = 0.0;
  for (int p = 1; p <= P; ++p) {
    double re[4] = {0.0};
    double im[4] = {0.0};
    for (int k = 0; k < M; ++k) {
      const double ang = TWO_PI * (double)((long)p * k % M) / (double)M;
      const double c = myb_cos(ang);
      const double s = myb_sin(ang);
      for (int j = 0; j < 4; ++j) {
        re[j] += samples[k][j] * c;
        im[j] -= samples[k][j] * s;
      }
    }
    /* |J_p|^2 = sum_ij |J_p,ij|^2: xx, yy, zz once, xy twice */
    double m2 = 0.0;
    for (int j = 0; j < 4; ++j) {
      const double q = (re[j] * re[j] + im[j] * im[j]) / ((double)M * (double)M);
      m2 += j == 1 ? 2.0 * q : q;
    }
    s_p += (double)p * m2;
    s_1 += m2;
  }
  *sum_p = s_p;
  *sum_1 = s_1;
}

static void test_tail_reaction_eccentric(void) {
  const double nu = 0.25;
  const double e = 0.5;
  const double sma = 200.0;
  const double rp = sma * (1.0 - e);
  const double vp = sqrt((1.0 + e) / rp);
  const double n = 1.0 / (sma * sqrt(sma));
  double sum_p = 0.0;
  double sum_1 = 0.0;
  tail_flux_sums(nu, sma, e, &sum_p, &sum_1);
  /* Parseval against Peters-Mathews: (2/5) sum |J_p|^2 = <I3:I3> / 5 = (32/5) nu^2 f(e) / a^5 */
  const double e2 = e * e;
  const double fe = (1.0 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e2 * e2) / pow(1.0 - e2, 3.5);
  const double flux_n = 0.4 * sum_1;
  const double peters = 32.0 / 5.0 * nu * nu * fe / pow(sma, 5.0);
  /* the tail flux, from eq. (FHtaildiss) averaged: (4 pi / 5) n sum_p p |J_p|^2 */
  const double flux_tail = 4.0 * PI / 5.0 * n * sum_p;
  run q = {.b = binary_of(nu, MYB_PN_4RR, 1.0)};
  const double y0[6] = {rp, 0.0, 0.0, 0.0, vp, 0.0};
  const myb_ias15_options o = {.epsilon = 1e-10};
  myb_ias15* d = NULL;
  MT_CHECK(myb_ias15_create(3, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_ias15_start(d, orbit_accel, &q, 0.0, y0) == MYB_OK);
  /* to the next periapsis: rdot from - to + after the apoapsis */
  double prev = radial_event(NULL, 0.0, y0);
  int passed_apo = 0;
  double te = 0.0;
  double ye[6] = {0.0};
  int done = 0;
  for (long k = 0; k < 1000000L && !done; ++k) {
    MT_CHECK(myb_ias15_step(d, 1e300) == MYB_OK);
    const double now = radial_event(NULL, 0.0, myb_ias15_state(d));
    if (prev > 0.0 && now <= 0.0) {
      passed_apo = 1;
    }
    if (passed_apo && prev < 0.0 && now >= 0.0) {
      int found = 0;
      MT_CHECK(myb_ias15_locate(d, radial_event, NULL, &found, &te, ye) == MYB_OK);
      done = found;
    }
    prev = now;
  }
  myb_ias15_destroy(d);
  MT_CHECK(done);
  const double e_start = 0.5 * vp * vp - 1.0 / rp;
  const double e_end = 0.5 * dot3(ye + 3, ye + 3) - 1.0 / norm3(ye);
  const double lost = -(e_end - e_start) * nu; /* E = mu E_N, mu = nu */
  const double expected = flux_tail * te;
  printf(
      "4PN tail radiation reaction, e = 0.5, a = 200: Newtonian flux by Fourier %.10e vs "
      "Peters-Mathews %.10e; energy lost in a period %.8e vs the tail flux %.8e "
      "(ratio %.7f; tail flux over 4 pi n F_N, 1 on a circle: %.4f)\n",
      flux_n, peters, lost, expected, lost / expected,
      flux_tail / (4.0 * PI * pow(n, 1.0) * flux_n));
  MT_CHECK(fabs(flux_n / peters - 1.0) < 1e-12);
  MT_CHECK(fabs(lost / expected - 1.0) < 1e-4);
}

/* ---- spin-orbit -------------------------------------------------------------------------------
 */

/* E / mu with the spins through 3.5PN (spin-orbit) and 4PN (no spin, instantaneous part) at
 * the state (x, v) */
static double spin_energy(const myb_pn_binary* b, const double* x, const double* v,
                          const double* c1, const double* c2, int with7) {
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
  if (b->orders & MYB_PN_SO) {
    e += pn4_spin_energy(3, nu, dm, 1.0 / r, rd, v2, nsv, nzv);
  }
  if (b->orders & MYB_PN_SO25) {
    e += pn4_spin_energy(5, nu, dm, 1.0 / r, rd, v2, nsv, nzv);
  }
  if ((b->orders & MYB_PN_SO35) && with7) {
    e += pn4_spin_energy(7, nu, dm, 1.0 / r, rd, v2, nsv, nzv);
  }
  return e;
}

/* J / mu with the spins: (x x v) J/J_N(4PN instantaneous) + l_S(1.5, 2.5, 3.5PN) + S / nu */
static void spin_angular(const myb_pn_binary* b, const double* x, const double* v, const double* c1,
                         const double* c2, int with7, double* out) {
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
  const double jj = pn4_angular_4pn(nu, 1.0 / r, rd, v2, myb_log(r / b->r0));
  for (int i = 0; i < 3; ++i) {
    out[i] = jj * h[i] + s[i] / nu;
  }
  static const int orders[3] = {3, 5, 7};
  static const unsigned flags[3] = {MYB_PN_SO, MYB_PN_SO25, MYB_PN_SO35};
  for (int k = 0; k < 3; ++k) {
    if (!(b->orders & flags[k]) || (k == 2 && !with7)) {
      continue;
    }
    double cn = 0.0;
    double cv = 0.0;
    double cs = 0.0;
    double cz = 0.0;
    pn4_spin_angular(orders[k], nu, dm, 1.0 / r, rd, v2, dot3(n, s), dot3(n, z), dot3(s, v),
                     dot3(z, v), &cn, &cv, &cs, &cz);
    for (int i = 0; i < 3; ++i) {
      out[i] += cn * n[i] + cv * v[i] + cs * s[i] + cz * z[i];
    }
  }
}

enum { FIT = 5 };

/* the polynomial of degree FIT - 1 through the points (x_k, q_k): its coefficients c[k] */
static void poly_fit(const double* x, const double* q, double* c) {
  /* Newton's divided differences, then the nested form expanded */
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
  for (int i = FIT - 2; i >= 0; --i) { /* poly = poly * (t - x_i) + d_i */
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

static void test_spin_energy_circular(void) {
  const double nu = 0.2;
  const double dm = sqrt(1.0 - 4.0 * nu);
  myb_pn_binary b = binary_of(nu, spins_nnlo, 1.0);
  const double x1 = b.gm1;
  const double x2 = b.gm2;
  const double chi1[3] = {0.0, 0.0, 0.5};
  const double chi2[3] = {0.0, 0.0, 0.8};
  const double sl = x1 * x1 * chi1[2] + x2 * x2 * chi2[2];
  const double zl = x2 * chi2[2] - x1 * chi1[2];
  double xs[FIT];
  double qs[FIT];
  for (int k = 0; k < FIT; ++k) {
    const double r = 100.0 * (double)(1 << k);
    double rr[2];
    double xx[2];
    for (int s = 0; s < 2; ++s) {
      const double sign = s == 0 ? 1.0 : -1.0;
      const double c1[3] = {0.0, 0.0, sign * chi1[2]};
      const double c2[3] = {0.0, 0.0, sign * chi2[2]};
      const double v2 = circular_v2(&b, r, c1, c2);
      const double x[3] = {r, 0.0, 0.0};
      const double v[3] = {0.0, sqrt(v2), 0.0};
      const double xv = myb_cbrt(v2 / (r * r));
      const double e = spin_energy(&b, x, v, c1, c2, 1);
      const double e1 = -0.75 - nu / 12.0;
      const double e2 = -27.0 / 8.0 + 19.0 / 8.0 * nu - nu * nu / 24.0;
      rr[s] = e / (-0.5 * xv) - 1.0 - e1 * xv - e2 * xv * xv;
      xx[s] = xv;
    }
    xs[k] = 0.5 * (xx[0] + xx[1]);
    qs[k] = 0.5 * (rr[0] - rr[1]) / (xs[k] * sqrt(xs[k]));
  }
  double c[FIT];
  poly_fit(xs, qs, c);
  /* eq. (Eofx) in G m^2 = 1 */
  const double s15 = 14.0 / 3.0 * sl + 2.0 * dm * zl;
  const double s25 = (11.0 - 61.0 / 9.0 * nu) * sl + dm * (3.0 - 10.0 / 3.0 * nu) * zl;
  const double s35 = (135.0 / 4.0 - 367.0 / 4.0 * nu + 29.0 / 12.0 * nu * nu) * sl +
                     dm * (27.0 / 4.0 - 39.0 * nu + 5.0 / 4.0 * nu * nu) * zl;
  printf(
      "E(x) with aligned spins, spin-linear part: x^(3/2) %.8f (Eofx %.8f), x^(5/2) %.6f "
      "(%.6f), x^(7/2) %.4f (%.4f)\n",
      c[0], s15, c[1], s25, c[2], s35);
  MT_CHECK(fabs(c[0] / s15 - 1.0) < 1e-6);
  MT_CHECK(fabs(c[1] / s25 - 1.0) < 1e-4);
  MT_CHECK(fabs(c[2] / s35 - 1.0) < 1e-2);
}

static void test_precession_circular(void) {
  const double nus[2] = {0.2, 0.1};
  for (int i = 0; i < 2; ++i) {
    const double nu = nus[i];
    myb_pn_binary b = binary_of(nu, spins_nnlo, 1.0);
    for (int body = 0; body < 2; ++body) {
      const double dm = (body == 0 ? 1.0 : -1.0) * sqrt(1.0 - 4.0 * nu);
      double xs[FIT];
      double qs[FIT];
      for (int k = 0; k < FIT; ++k) {
        const double r = 100.0 * (double)(1 << k);
        const double tiny[3] = {1e-7, 0.0, 0.0};
        const double zero[3] = {0.0, 0.0, 0.0};
        const double* c1 = body == 0 ? tiny : zero;
        const double* c2 = body == 0 ? zero : tiny;
        const double v2 = circular_v2(&b, r, c1, c2);
        const double x[3] = {r, 0.0, 0.0};
        const double v[3] = {0.0, sqrt(v2), 0.0};
        double a[3];
        double d1[3];
        double d2[3];
        MT_CHECK(myb_pn_binary_evaluate(&b, x, v, c1, c2, a, d1, d2) == MYB_OK);
        const double om = norm3(body == 0 ? d1 : d2) / 1e-7;
        const double xv = myb_cbrt(v2 / (r * r));
        xs[k] = xv;
        qs[k] = om / (xv * xv * sqrt(xv));
      }
      double c[FIT];
      poly_fit(xs, qs, c);
      /* eq. (precessionCO) */
      const double c0 = 0.75 + 0.5 * nu - 0.75 * dm;
      const double c1 =
          9.0 / 16.0 + 1.25 * nu - nu * nu / 24.0 + dm * (-9.0 / 16.0 + 5.0 / 8.0 * nu);
      const double c2 = 27.0 / 32.0 + 3.0 / 16.0 * nu - 105.0 / 32.0 * nu * nu -
                        nu * nu * nu / 48.0 +
                        dm * (-27.0 / 32.0 + 39.0 / 8.0 * nu - 5.0 / 32.0 * nu * nu);
      printf(
          "precession of spin %d on circular orbits, nu %.2f: x^(5/2) %.10f (precessionCO "
          "%.10f), x^(7/2) %.7f (%.7f), x^(9/2) %.4f (%.4f)\n",
          body + 1, nu, c[0], c0, c[1], c1, c[2], c2);
      MT_CHECK(fabs(c[0] - c0) < 1e-8 * fabs(c0));
      MT_CHECK(fabs(c[1] - c1) < 1e-5 * fabs(c1) + 1e-7);
      MT_CHECK(fabs(c[2] - c2) < 2e-2 * fabs(c2) + 1e-3);
    }
  }
}

typedef struct spin_run {
  myb_pn_binary b;
  int inst4; /* add the instantaneous 4PN terms (without the tail) */
} spin_run;

/* 9 "positions" (x, p1, p2) whose velocities are (v, chi1, chi2) */
static myb_status spin_accel(void* context, double t, const double* x, const double* v, double* a) {
  (void)t;
  const spin_run* q = context;
  const myb_status st = myb_pn_binary_evaluate(&q->b, x, v, v + 3, v + 6, a, a + 3, a + 6);
  if (st == MYB_OK && q->inst4) {
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

enum { GRID = 400 };

/* the energy and angular momentum (with and without the 3.5PN terms of E and L) at GRID + 1
 * times over `orbits` orbits of period `period` */
static void spin_history(spin_run* q, const double y0[18], double period, int orbits, double* e7,
                         double* e5, double (*j7)[3], double (*j5)[3]) {
  const myb_ias15_options o = {.epsilon = 1e-10, .block = 3};
  myb_ias15* d = NULL;
  MT_CHECK(myb_ias15_create(9, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_ias15_start(d, spin_accel, q, 0.0, y0) == MYB_OK);
  /* free steps, the states at the grid times from the dense output of each step */
  int g = 0;
  for (long k = 0; k < 10000000L && g <= GRID; ++k) {
    double y[18];
    const double tg = period * orbits * (double)g / (double)GRID;
    if (g == 0) {
      memcpy(y, y0, sizeof y);
    } else if (tg <= myb_ias15_time(d)) {
      MT_CHECK(myb_ias15_dense(d, tg, y) == MYB_OK);
    } else {
      MT_CHECK(myb_ias15_step(d, 1e300) == MYB_OK);
      continue;
    }
    e7[g] = spin_energy(&q->b, y, y + 9, y + 12, y + 15, 1);
    e5[g] = spin_energy(&q->b, y, y + 9, y + 12, y + 15, 0);
    spin_angular(&q->b, y, y + 9, y + 12, y + 15, 1, j7[g]);
    spin_angular(&q->b, y, y + 9, y + 12, y + 15, 0, j5[g]);
    ++g;
  }
  MT_CHECK(g == GRID + 1);
  myb_ias15_destroy(d);
}

static void test_spin_conservation(void) {
  /* G m = c = 1, nu = 0.2; an orbit of rp = 15, ra = 30 (and twice as large) in the xy plane,
   * oblique spins (|chi| ~ 0.04, so that the terms cubic in the spins stay far below):
   * the spin-linear part of the changes of E and J, (run(+chi) - run(-chi)) / 2, is the remainder
   * of the spin-orbit terms, O(S / c^9), against O(S / c^7) of the 3.5PN terms. Its expected
   * scaling with the orbit: E (of order 1/r) changes by r^-5.5 (2^5.5 = 45 when the orbit is
   * doubled; without the 3.5PN energy r^-4.5, 22.6); J (of order r^1/2) by r^-4 (16; without the
   * 3.5PN angular momentum r^-3, 8).
   *
   * The spin-free terms must be conserved to an order beyond that: the orbit of +chi and that of
   * -chi differ at the order of the spins, so a spin-free change of E or J of relative order eps
   * makes a spin-linear one of order eps times the spin effect. Through 3PN the spin-free
   * remainder is 4PN, and it tilts with the orbital plane (whose normal differs by S / (c L)
   * between the two runs, a "0.5PN" spin effect): a spin-linear change of J of the order of the
   * one looked for. The test therefore adds the instantaneous 4PN terms to the acceleration
   * (spin_accel, without the tail, whose energy is not local) and takes E and J through 4PN; at
   * rp = 15 the spin-free 5PN remainder is still not small (PN coefficients grow as fast as G m /
   * (r c^2) falls), and its tilt still dominates the components of J in the orbital plane:
   * the component of J along the normal of the initial plane, which the tilt does not reach at
   * first order, is the one compared. */
  static double e7p[GRID + 1];
  static double e5p[GRID + 1];
  static double e7m[GRID + 1];
  static double e5m[GRID + 1];
  static double j7p[GRID + 1][3];
  static double j5p[GRID + 1][3];
  static double j7m[GRID + 1][3];
  static double j5m[GRID + 1][3];
  const double nu = 0.2;
  const double chi1[3] = {0.015, -0.025, 0.03};
  const double chi2[3] = {0.035, 0.01, -0.02};
  double de7[2];
  double de5[2];
  double dj7[2];
  double dj5[2];
  for (int size = 0; size < 2; ++size) {
    const double rp = 15.0 * (double)(1 << size);
    const double ra = 2.0 * rp;
    const double sma = 0.5 * (rp + ra);
    const double vp = sqrt(ra / (sma * rp));
    const double period = TWO_PI * sma * sqrt(sma);
    spin_run q = {.b = binary_of(nu, spins_nnlo, 1.0), .inst4 = 1};
    double y0[18] = {0};
    y0[0] = rp;
    y0[10] = vp;
    memcpy(y0 + 12, chi1, sizeof chi1);
    memcpy(y0 + 15, chi2, sizeof chi2);
    spin_history(&q, y0, period, 2, e7p, e5p, j7p, j5p);
    for (int i = 12; i < 18; ++i) {
      y0[i] = -y0[i];
    }
    spin_history(&q, y0, period, 2, e7m, e5m, j7m, j5m);
    de7[size] = de5[size] = dj7[size] = dj5[size] = 0.0;
    for (int g = 1; g <= GRID; ++g) {
      de7[size] = fmax(de7[size], fabs((e7p[g] - e7p[0]) - (e7m[g] - e7m[0])) / 2.0);
      de5[size] = fmax(de5[size], fabs((e5p[g] - e5p[0]) - (e5m[g] - e5m[0])) / 2.0);
      dj7[size] = fmax(dj7[size], fabs((j7p[g][2] - j7p[0][2]) - (j7m[g][2] - j7m[0][2])) / 2.0);
      dj5[size] = fmax(dj5[size], fabs((j5p[g][2] - j5p[0][2]) - (j5m[g][2] - j5m[0][2])) / 2.0);
    }
    printf(
        "spin-orbit through NNLO, rp %.0f, ra %.0f: spin-linear change of E %.2e (without "
        "its 3.5PN terms %.2e) of %.3e; of J along the normal %.2e (without l_S 3.5PN %.2e) "
        "of %.3f\n",
        rp, ra, de7[size], de5[size], fabs(e7p[0]), dj7[size], dj5[size], j7p[0][2]);
  }
  /* the 3.5PN terms are needed, and the remainder falls faster than they do, as expected */
  MT_CHECK(de7[0] < 0.2 * de5[0]);
  MT_CHECK(dj7[0] < 0.2 * dj5[0]);
  const double re7 = de7[0] / de7[1];
  const double re5 = de5[0] / de5[1];
  const double rj7 = dj7[0] / dj7[1];
  const double rj5 = dj5[0] / dj5[1];
  printf(
      "spin-orbit through NNLO, orbit x 2: the spin-linear change of E falls by %.1f (%.1f "
      "without the 3.5PN energy; expected 45 and 22.6), of J by %.1f (%.1f; 16 and 8)\n",
      re7, re5, rj7, rj5);
  MT_CHECK(re7 > 35.0 && re5 < 30.0);
  MT_CHECK(rj7 > 12.0 && rj5 < 11.0);
}

/* ---- spin-orbit radiation reaction ------------------------------------------------------------
 */

/* (P / (-(32/5) nu x^5) - 1) / x^(3/2) on the circular orbit of radius r with aligned spins, P =
 * dE/dt = (dE/dv).a_RR, E / mu = v^2 / 2 - 1 / r + e_S(1.5PN), a_RR of `reaction`, the orbit of
 * the conservative orders `orbit`. */
static double spin_flux_factor(double nu, double r, double sign, unsigned orbit, unsigned reaction,
                               const double* chi1, const double* chi2) {
  myb_pn_binary b = binary_of(nu, orbit, 1.0);
  const double c1[3] = {sign * chi1[0], sign * chi1[1], sign * chi1[2]};
  const double c2[3] = {sign * chi2[0], sign * chi2[1], sign * chi2[2]};
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
  return (p / (-32.0 / 5.0 * nu * x5) - 1.0) / (xv * sqrt(xv));
}

static void test_spin_flux_circular(void) {
  /* Newton + the leading spin-orbit terms + 2.5PN + Will's spin-orbit reaction: in this model
   * the only scale is r, so the spin-linear part of dE/dt at fixed x is exactly Kidder's factor
   * times x^(3/2) of the Newtonian flux (no higher-order spin-linear term); what remains is cubic
   * in the spins, (chi x^(3/2))^2 relative, 1e-6 at r = 30 and 1e-9 at r = 300 */
  const double nu = 0.2;
  const double dm = sqrt(1.0 - 4.0 * nu);
  const myb_pn_binary b0 = binary_of(nu, 0, 1.0);
  const double chi1[3] = {0.0, 0.0, 0.1};
  const double chi2[3] = {0.0, 0.0, 0.2};
  const double sl = b0.gm1 * b0.gm1 * chi1[2] + b0.gm2 * b0.gm2 * chi2[2];
  const double zl = b0.gm2 * chi2[2] - b0.gm1 * chi1[2];
  const double kidder = -4.0 * sl - 1.25 * dm * zl;
  const double radii[3] = {30.0, 300.0, 3000.0};
  for (int k = 0; k < 3; ++k) {
    const double kp =
        spin_flux_factor(nu, radii[k], 1.0, MYB_PN_SO, MYB_PN_25 | MYB_PN_SORR, chi1, chi2);
    const double km =
        spin_flux_factor(nu, radii[k], -1.0, MYB_PN_SO, MYB_PN_25 | MYB_PN_SORR, chi1, chi2);
    const double lin = 0.5 * (kp - km);
    const double bp = spin_flux_factor(nu, radii[k], 1.0, MYB_PN_SO, MYB_PN_25, chi1, chi2);
    const double bm = spin_flux_factor(nu, radii[k], -1.0, MYB_PN_SO, MYB_PN_25, chi1, chi2);
    printf(
        "spin-orbit radiation reaction, circular, r %.0f: spin-linear factor of dE/dt %.12f "
        "(Kidder %.12f; without a_SO-RR %.6f)\n",
        radii[k], lin, kidder, 0.5 * (bp - bm));
    MT_CHECK(fabs(lin / kidder - 1.0) < (k == 0 ? 1e-5 : 1e-8));
    MT_CHECK(fabs(0.5 * (bp - bm) / kidder - 1.0) > 0.1);
  }
}

static void test_spin_flux_integrated(void) {
  /* N + SO + 2.5PN + SO radiation reaction, aligned spins, a quasi-circular orbit at r = 60:
   * the energy E / mu = v^2/2 - 1/r + e_S(1.5PN) lost over 3 orbits against the integral of the
   * flux (32/5) nu x^5 [1 + x^(3/2) (-4 S_l - 5/4 dm Sigma_l)] along the orbit, x from the
   * instantaneous angular velocity */
  enum { STEPS = 600 };

  const double nu = 0.2;
  const double dm = sqrt(1.0 - 4.0 * nu);
  const double r = 60.0;
  const double chi1[3] = {0.0, 0.0, 0.5};
  const double chi2[3] = {0.0, 0.0, 0.8};
  double lost[2][2];
  double flux[2][2];
  for (int with = 0; with < 2; ++with) {
    for (int s = 0; s < 2; ++s) {
      const double sign = s == 0 ? 1.0 : -1.0;
      spin_run q = {.b = binary_of(nu, MYB_PN_SO, 1.0)};
      const double c1[3] = {0.0, 0.0, sign * chi1[2]};
      const double c2[3] = {0.0, 0.0, sign * chi2[2]};
      const double v2 = circular_v2(&q.b, r, c1, c2);
      q.b.orders = MYB_PN_SO | MYB_PN_25 | (with ? MYB_PN_SORR : 0u);
      const double g = 1.0 / r;
      const double sl = q.b.gm1 * q.b.gm1 * c1[2] + q.b.gm2 * q.b.gm2 * c2[2];
      const double zl = q.b.gm2 * c2[2] - q.b.gm1 * c1[2];
      const double kfac = -4.0 * sl - 1.25 * dm * zl;
      double y0[18] = {0};
      y0[0] = r;
      y0[9] = -64.0 / 5.0 * nu * g * g * g;
      y0[10] = sqrt(v2);
      memcpy(y0 + 12, c1, sizeof c1);
      memcpy(y0 + 15, c2, sizeof c2);
      const double period = TWO_PI / sqrt(v2 / (r * r));
      const myb_ias15_options o = {.epsilon = 1e-10, .block = 3};
      myb_ias15* d = NULL;
      MT_CHECK(myb_ias15_create(9, &o, NULL, &d) == MYB_OK);
      MT_CHECK(myb_ias15_start(d, spin_accel, &q, 0.0, y0) == MYB_OK);
      double e0 = 0.0;
      double integral = 0.0;
      double f_prev = 0.0;
      double e_now = 0.0;
      for (int k = 0; k <= STEPS; ++k) {
        if (k > 0) {
          MT_CHECK(myb_ias15_integrate(d, 3.0 * period * (double)k / (double)STEPS) == MYB_OK);
        }
        const double* y = myb_ias15_state(d);
        const double* xx = y;
        const double* vv = y + 9;
        const double rr = norm3(xx);
        double h[3];
        cross3(xx, vv, h);
        const double om = norm3(h) / (rr * rr);
        const double xv = myb_cbrt(om * om);
        const double x5 = xv * xv * xv * xv * xv;
        const double f = 32.0 / 5.0 * nu * x5 * (1.0 + kfac * xv * sqrt(xv));
        const double v2n = dot3(vv, vv);
        const double n[3] = {xx[0] / rr, xx[1] / rr, xx[2] / rr};
        double sv3[3];
        double z[3];
        for (int i = 0; i < 3; ++i) {
          sv3[i] = q.b.gm1 * q.b.gm1 * y[12 + i] + q.b.gm2 * q.b.gm2 * y[15 + i];
          z[i] = q.b.gm2 * y[15 + i] - q.b.gm1 * y[12 + i];
        }
        double sxv[3];
        double zxv[3];
        cross3(sv3, vv, sxv);
        cross3(z, vv, zxv);
        e_now = 0.5 * v2n - 1.0 / rr +
                pn4_spin_energy(3, nu, dm, 1.0 / rr, dot3(n, vv), v2n, dot3(n, sxv), dot3(n, zxv));
        if (k == 0) {
          e0 = e_now;
        } else {
          integral += 0.5 * (f + f_prev) * 3.0 * period / (double)STEPS;
        }
        f_prev = f;
      }
      myb_ias15_destroy(d);
      lost[with][s] = e0 - e_now;
      flux[with][s] = integral;
    }
  }
  for (int with = 0; with < 2; ++with) {
    const double lin_lost = 0.5 * (lost[with][0] - lost[with][1]);
    const double lin_flux = 0.5 * (flux[with][0] - flux[with][1]);
    printf(
        "spin-orbit radiation reaction%s, integrated quasi-circular orbit r 60, 3 orbits: "
        "energy lost %.9e / %.9e (+chi / -chi), flux integral %.9e / %.9e; spin-linear part "
        "%.6e vs %.6e (ratio %.5f)\n",
        with ? "" : " left out", lost[with][0], lost[with][1], flux[with][0], flux[with][1],
        lin_lost, lin_flux, lin_lost / lin_flux);
    if (with) {
      MT_CHECK(fabs(lost[with][0] / flux[with][0] - 1.0) < 1e-4);
      MT_CHECK(fabs(lost[with][1] / flux[with][1] - 1.0) < 1e-4);
      MT_CHECK(fabs(lin_lost / lin_flux - 1.0) < 1e-2);
    } else {
      MT_CHECK(fabs(lin_lost / lin_flux - 1.0) > 0.1);
    }
  }
}

/* ---- domain and cost ----------------------------------------------------------------------- */

static void test_domain(void) {
  const double x[3] = {10.0, 0.0, 0.0};
  const double vb[3] = {0.0, 0.3, 0.0};  /* bound: v^2 < 2 / r */
  const double vu[3] = {0.0, 0.46, 0.0}; /* unbound */
  double a[3];
  myb_pn_binary b = binary_of(0.25, MYB_PN_4, 0.0);
  MT_CHECK(myb_pn_binary_evaluate(&b, x, vb, NULL, NULL, a, NULL, NULL) == MYB_E_DOMAIN);
  b.r0 = 1.0;
  MT_CHECK(myb_pn_binary_evaluate(&b, x, vb, NULL, NULL, a, NULL, NULL) == MYB_OK);
  MT_CHECK(myb_pn_binary_evaluate(&b, x, vu, NULL, NULL, a, NULL, NULL) == MYB_E_DOMAIN);
  b.orders = MYB_PN_4RR;
  MT_CHECK(myb_pn_binary_evaluate(&b, x, vu, NULL, NULL, a, NULL, NULL) == MYB_E_DOMAIN);
  /* e = 0.8 at periapsis: allowed; e = 0.9: more than 1024 harmonics */
  const double v8[3] = {0.0, sqrt(1.8 / 10.0), 0.0};
  const double v9[3] = {0.0, sqrt(1.9 / 10.0), 0.0};
  MT_CHECK(myb_pn_binary_evaluate(&b, x, v8, NULL, NULL, a, NULL, NULL) == MYB_OK);
  MT_CHECK(myb_pn_binary_evaluate(&b, x, v9, NULL, NULL, a, NULL, NULL) == MYB_E_DOMAIN);
  b.orders = 4096;
  MT_CHECK(myb_pn_binary_evaluate(&b, x, vb, NULL, NULL, a, NULL, NULL) == MYB_E_DOMAIN);
  /* the spin terms without spins vanish; the old entry point refuses 4PN (no r0) */
  b.orders = MYB_PN_SO25 | MYB_PN_SO35 | MYB_PN_SORR;
  double d1[3];
  double d2[3];
  MT_CHECK(myb_internal_pn_binary_terms(&b, x, vb, NULL, NULL, a, d1, d2) == MYB_OK);
  MT_CHECK(a[0] == 0.0 && a[1] == 0.0 && a[2] == 0.0);
  MT_CHECK(myb_pn_binary_acceleration(0.5, 0.5, 1.0, MYB_PN_4, x, vb, a) == MYB_E_DOMAIN);
  /* aligned spins: no precession, planar motion with every term */
  b = binary_of(0.2, MYB_PN_ALL | MYB_PN_4 | MYB_PN_4RR | MYB_PN_SO25 | MYB_PN_SO35 | MYB_PN_SORR,
                1.0);
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
  test_coefficients();
  test_spin_terms();
  test_tail_reference();
  test_tail_circular();
  test_kepler_law();
  test_energy_circular();
  test_periastron_advance_4pn();
  test_tail_reaction_circular();
  test_tail_reaction_eccentric();
  test_spin_energy_circular();
  test_precession_circular();
  test_spin_conservation();
  test_spin_flux_circular();
  test_spin_flux_integrated();
  test_domain();
  return mt_finish("test_pn4");
}
