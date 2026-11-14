/* SPDX-License-Identifier: BSD-3-Clause */
/* Compact binaries: the 4PN tail on any osculating conic (MYB_PN_TAILX: the hereditary integrals
 * of arXiv:1711.00283 eqs. acctailCM and acctaildiss in the time domain), units G m = c = 1:
 *   - on ellipses the time domain against the Fourier series of MYB_PN_4 (e = 0 to 0.8): the
 *     same tail, and the switch between them at e ~ 0.83 continuous;
 *   - against an independent mpmath evaluation (tests/data/pn/tail_open.txt: hyperbolic or
 *     eccentric anomaly, tanh-sinh quadrature) on very eccentric ellipses, near-parabolic conics
 *     and hyperbolas;
 *   - continuity through e = 1: on conics of the same periapsis distance and true anomaly with
 *     e = 1 -+ delta the tail differs by O(delta);
 *   - hyperbolic encounters along the Newtonian orbit: the energy radiated by the Burke-Thorne
 *     force against the published closed form (Turner 1977, eq. EGW of arXiv:1706.06877, a check
 *     of the quadrature over the orbit), the energy and angular momentum radiated by the
 *     dissipative tail against the frequency-domain values from the published Fourier transforms
 *     of the quadrupole (Bini and Damour, arXiv:1706.06877, eq. rel_x2y2xx_final), and the
 *     scattering angle of the conservative tail against (1/nu) dW/dj of the same paper with the
 *     time-dependent Hadamard scale 2 r of the force, and against its leading order in 1/e
 *     (eq. WLO1) at large eccentricity. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/rel.h"
#include "rel/rel_internal.h"

#define PI 3.14159265358979323846
#define TWO_PI 6.283185307179586

static double norm3(const double* a) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

[[maybe_unused]] static void cross3(const double* a, const double* b, double* out) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

[[maybe_unused]] static double dot3(const double* a, const double* b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double rel3(const double* a, const double* b) {
  const double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  return norm3(d) / norm3(b);
}

/* The state of true anomaly f on the conic of periapsis distance q and eccentricity e (G m = 1),
 * rotated out of the xy plane (inclination 0.7, node 0.4, argument of periapsis 1.1). */
static void conic_state(double q, double e, double f, double x[3], double v[3]) {
  const double p = q * (1.0 + e);
  const double r = p / (1.0 + e * myb_cos(f));
  const double pos[3] = {r * myb_cos(f), r * myb_sin(f), 0.0};
  const double k = 1.0 / sqrt(p);
  const double vel[3] = {-k * myb_sin(f), k * (e + myb_cos(f)), 0.0};
  const double ci = myb_cos(0.7);
  const double si = myb_sin(0.7);
  const double cn = myb_cos(0.4);
  const double sn = myb_sin(0.4);
  const double cw = myb_cos(1.1);
  const double sw = myb_sin(1.1);
  const double* in[2] = {pos, vel};
  double* out[2] = {x, v};
  for (int s = 0; s < 2; ++s) {
    const double* u = in[s];
    const double x1 = cw * u[0] - sw * u[1];
    const double y1 = sw * u[0] + cw * u[1];
    const double y2 = ci * y1 - si * u[2];
    const double z2 = si * y1 + ci * u[2];
    out[s][0] = cn * x1 - sn * y2;
    out[s][1] = sn * x1 + cn * y2;
    out[s][2] = z2;
  }
}

/* ---- ellipses: the time domain against the Fourier series ----------------------------------- */

static void test_ellipses(void) {
  /* the difference relative to the tail at the periapsis of the same orbit: near the apoapsis of
   * an eccentric orbit the tail is 1e3-1e5 times smaller, and the rounding of the periapsis part of
   * either method dominates there */
  static const double es[6] = {0.0, 0.05, 0.3, 0.5, 0.7, 0.8};
  static const double fs[3] = {0.0, 1.9, 3.0};
  double worst_c = 0.0;
  double worst_d = 0.0;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double q = 12.0 + 3.0 * (double)j;
      double x[3];
      double v[3];
      conic_state(q, es[i], fs[j], x, v);
      double xp[3];
      double vp[3];
      conic_state(q, es[i], 0.0, xp, vp);
      double pc[3];
      double pd[3];
      MT_CHECK(myb_internal_pn_tail(0.2, xp, vp, 0, pc, pd) == MYB_OK);
      double fc[3];
      double fd[3];
      double tc[3];
      double td[3];
      MT_CHECK(myb_internal_pn_tail(0.2, x, v, 0, fc, fd) == MYB_OK);
      MT_CHECK(myb_internal_pn_tail(0.2, x, v, 1, tc, td) == MYB_OK);
      const double dc[3] = {tc[0] - fc[0], tc[1] - fc[1], tc[2] - fc[2]};
      const double dd[3] = {td[0] - fd[0], td[1] - fd[1], td[2] - fd[2]};
      const double ec = norm3(dc) / norm3(pc);
      const double ed = norm3(dd) / norm3(pd);
      worst_c = fmax(worst_c, ec);
      worst_d = fmax(worst_d, ed);
      printf(
          "tail on an ellipse, e %.2f, q %.0f, true anomaly %.1f: time domain vs Fourier "
          "series, relative to the tail at the periapsis: conservative %.1e, dissipative "
          "%.1e (of themselves %.1e, %.1e)\n",
          es[i], q, fs[j], ec, ed, rel3(tc, fc), rel3(td, fd));
    }
  }
  printf(
      "tail on ellipses, e 0 to 0.8, 18 states: time domain vs Fourier series, worst "
      "conservative %.1e, dissipative %.1e\n",
      worst_c, worst_d);
  MT_CHECK(worst_c < 1e-12);
  MT_CHECK(worst_d < 1e-12);
}

/* through the public entry point: below about e = 0.83 the Fourier series, above the time
 * domain with MYB_PN_TAILX (MYB_E_DOMAIN without it) */
static void test_switch(void) {
  const myb_pn_binary b0 = {.gm1 = 0.6, .gm2 = 0.4, .c = 1.0, .orders = MYB_PN_4RR};
  myb_pn_binary bx = b0;
  bx.orders = MYB_PN_4RR | MYB_PN_TAILX;
  double prev[3] = {0.0, 0.0, 0.0};
  double worst = 0.0;
  int fourier = 0;
  int time_domain = 0;
  for (int k = 0; k <= 20; ++k) {
    const double e = 0.82 + 0.001 * (double)k;
    double x[3];
    double v[3];
    conic_state(20.0, e, 0.5, x, v);
    double a0[3];
    double ax[3];
    const myb_status s0 = myb_internal_pn_binary_terms(&b0, x, v, NULL, NULL, a0, NULL, NULL);
    MT_CHECK(myb_internal_pn_binary_terms(&bx, x, v, NULL, NULL, ax, NULL, NULL) == MYB_OK);
    if (s0 == MYB_OK) {
      ++fourier;
      for (int i = 0; i < 3; ++i) {
        MT_CHECK(a0[i] == ax[i]); /* below the switch MYB_PN_TAILX changes nothing */
      }
    } else {
      MT_CHECK(s0 == MYB_E_DOMAIN);
      ++time_domain;
    }
    if (k > 0) {
      /* the change from one step to the next is smooth: compare with the neighbours */
      worst = fmax(worst, rel3(ax, prev));
    }
    memcpy(prev, ax, sizeof prev);
  }
  printf(
      "public entry, e 0.820 to 0.840 by 0.001: %d states by the Fourier series, %d in the "
      "time domain; largest relative change between neighbours %.2e\n",
      fourier, time_domain, worst);
  MT_CHECK(fourier > 0 && time_domain > 0);
  MT_CHECK(worst < 0.02);
}

/* ---- through e = 1 ----------------------------------------------------------------------------
 */

static void test_parabolic(void) {
  static const double deltas[4] = {1e-3, 1e-5, 1e-7, 1e-9};
  double prev_c = 0.0;
  double prev_d = 0.0;
  for (int k = 0; k < 4; ++k) {
    double c[2][3];
    double d[2][3];
    for (int s = 0; s < 2; ++s) {
      const double e = 1.0 + (s == 0 ? -deltas[k] : deltas[k]);
      double x[3];
      double v[3];
      conic_state(20.0, e, 0.6, x, v);
      MT_CHECK(myb_internal_pn_tail(0.25, x, v, 1, c[s], d[s]) == MYB_OK);
    }
    const double jc = rel3(c[1], c[0]);
    const double jd = rel3(d[1], d[0]);
    printf(
        "tail across e = 1, q 20, true anomaly 0.6, e = 1 -+ %.0e: relative jump "
        "conservative %.2e, dissipative %.2e\n",
        deltas[k], jc, jd);
    if (k > 0) {
      /* O(delta), down to the rounding of the state */
      MT_CHECK(jc < 0.02 * prev_c + 1e-11);
      MT_CHECK(jd < 0.02 * prev_d + 1e-11);
    }
    prev_c = jc;
    prev_d = jd;
  }
}

/* ---- against the mpmath reference ---------------------------------------------------------------
 */

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

static void test_reference(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/pn/tail_open.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  int rows = 0;
  double worst_c = 0.0;
  double worst_d = 0.0;
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
    double tc[3];
    double td[3];
    MT_CHECK(myb_internal_pn_tail(nu, x, v, 1, tc, td) == MYB_OK);
    /* and through the public entry point: the dissipative part alone */
    const double d = sqrt(1.0 - 4.0 * nu);
    const myb_pn_binary b = {.gm1 = 0.5 * (1.0 + d),
                             .gm2 = 0.5 * (1.0 - d),
                             .c = 1.0,
                             .orders = MYB_PN_4RR | MYB_PN_TAILX};
    double ad[3];
    MT_CHECK(myb_internal_pn_binary_terms(&b, x, v, NULL, NULL, ad, NULL, NULL) == MYB_OK);
    /* relative to the tail at the periapsis of the same conic (near the apoapsis of an ellipse
     * the tail is far smaller than the L1 norm of its integrand) */
    double h[3];
    cross3(x, v, h);
    double ev[3];
    cross3(v, h, ev);
    const double r = norm3(x);
    for (int i = 0; i < 3; ++i) {
      ev[i] -= x[i] / r;
    }
    const double ecc = norm3(ev);
    const double q = dot3(h, h) / (1.0 + ecc);
    double qhat[3];
    cross3(h, ev, qhat);
    const double vq = sqrt((1.0 + ecc) / q) / norm3(qhat);
    double xp[3];
    double vp[3];
    for (int i = 0; i < 3; ++i) {
      xp[i] = q * ev[i] / ecc;
      vp[i] = vq * qhat[i];
    }
    double pc[3];
    double pd[3];
    MT_CHECK(myb_internal_pn_tail(nu, xp, vp, 1, pc, pd) == MYB_OK);
    const double dc[3] = {tc[0] - w[7], tc[1] - w[8], tc[2] - w[9]};
    const double dd[3] = {td[0] - w[10], td[1] - w[11], td[2] - w[12]};
    const double dda[3] = {ad[0] - w[10], ad[1] - w[11], ad[2] - w[12]};
    const double ec = norm3(dc) / fmax(norm3(w + 7), norm3(pc));
    const double ed = fmax(norm3(dd), norm3(dda)) / fmax(norm3(w + 10), norm3(pd));
    printf(
        "tail vs mpmath, e %.7f, r %.1f: conservative %.1e, dissipative %.1e (of themselves "
        "%.1e, %.1e)\n",
        w[13], r, ec, ed, rel3(tc, w + 7), rel3(td, w + 10));
    worst_c = fmax(worst_c, ec);
    worst_d = fmax(worst_d, ed);
    ++rows;
  }
  fclose(f);
  printf(
      "tail on very eccentric, near-parabolic and hyperbolic conics vs mpmath (%d states): "
      "worst conservative %.1e, dissipative %.1e\n",
      rows, worst_c, worst_d);
  MT_CHECK(rows == 15);
  MT_CHECK(worst_c < 1e-12);
  MT_CHECK(worst_d < 1e-12);
}

/* ---- hyperbolic encounters
 * ------------------------------------------------------------------------ */

enum { GL = 16 };

/* Gauss-Legendre nodes and weights on [-1, 1] by Newton on P_16 */
static void gauss_legendre(double* xg, double* wg) {
  for (int i = 0; i < GL; ++i) {
    double z = myb_cos(PI * ((double)i + 0.75) / ((double)GL + 0.5));
    double dp = 0.0;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0;
      double p1 = z;
      for (int k = 2; k <= GL; ++k) {
        const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / (double)k;
        p0 = p1;
        p1 = p2;
      }
      dp = (double)GL * (z * p1 - p0) / (z * z - 1.0);
      const double dz = p1 / dp;
      z -= dz;
      if (fabs(dz) < 1e-16) {
        break;
      }
    }
    xg[i] = z;
    wg[i] = 2.0 / ((1.0 - z * z) * dp * dp);
  }
}

/* the encounter of eccentricity e and speed at infinity vinf (G m = c = 1) in the xy plane,
 * periapsis on +x, counterclockwise: the state at the hyperbolic anomaly H and dt/dH */
static void hyperbola(double e, double vinf, double h_anom, double x[3], double v[3],
                      double* dtdh) {
  const double a = 1.0 / (vinf * vinf);
  const double se = sqrt(e * e - 1.0);
  const double ch = myb_cosh(h_anom);
  const double sh = myb_sinh(h_anom);
  const double r = a * (e * ch - 1.0);
  x[0] = a * (e - ch);
  x[1] = a * se * sh;
  x[2] = 0.0;
  const double k = sqrt(a) / r;
  v[0] = -k * sh;
  v[1] = k * se * ch;
  v[2] = 0.0;
  *dtdh = r / (a * vinf * vinf * vinf);
}

typedef struct encounter_result {
  double de;      /* energy radiated by the dissipative tail (positive) */
  double dj;      /* change of the orbital angular momentum along z by it */
  double dchi;    /* change of the scattering angle by the conservative tail */
  double de_cons; /* net energy change by the conservative tail (should vanish) */
} encounter_result;

/* The integrands at the hyperbolic anomaly H: v.a_diss, (x x a_diss)_z, v.a_cons and the rate
 * of the Runge-Lenz vector a_cons x h + v x (x x a_cons) (components 3, 4), times dt/dH. */
static void encounter_rates(double nu, double e, double vinf, double hh, double out[6]) {
  double x[3];
  double v[3];
  double dtdh = 0.0;
  hyperbola(e, vinf, hh, x, v, &dtdh);
  double fc[3];
  double fd[3];
  MT_CHECK(myb_internal_pn_tail(nu, x, v, 1, fc, fd) == MYB_OK);
  double xf[3];
  cross3(x, fd, xf);
  double h[3];
  cross3(x, v, h);
  double t1[3];
  double t2[3];
  double t3[3];
  cross3(fc, h, t1);
  cross3(x, fc, t2);
  cross3(v, t2, t3);
  out[0] = -nu * dot3(v, fd) * dtdh;
  out[1] = nu * xf[2] * dtdh;
  out[2] = dot3(v, fc) * dtdh;
  out[3] = (t1[0] + t3[0]) * dtdh;
  out[4] = (t1[1] + t3[1]) * dtdh;
  out[5] = 0.0;
}

/* Integrals along the Newtonian hyperbola (16 Gauss-Legendre nodes on panels of 1/8 in H) between
 * -H_k and H_k for the times T_k = 40, 80, 160 times the passage time j / v_inf^2 (and 1 / n if
 * longer), each completed by the leading asymptotic tail of the tail forces, (24/5) DI2 v / |t|^3
 * (times sgn t for the dissipative one; DI2 = 2 nu (v_+^<i v_+^j> - v_-^<i v_-^j>), the jump of
 * the second derivative of the quadrupole), which gives the integrals of x x a and of the
 * Runge-Lenz rate, converging as 1 / T, a remainder (a + b ln T) / T^2: the three cutoffs
 * determine it and the limit. Far beyond the far field cannot be evaluated: there the tail is
 * t^-3 relative to the L1 norm of its integrand. */
static encounter_result encounter(double nu, double e, double vinf) {
  double xg[GL];
  double wg[GL];
  gauss_legendre(xg, wg);
  const double j = sqrt(e * e - 1.0) / vinf;
  const double n = vinf * vinf * vinf;
  const double tau = fmax(j / (vinf * vinf), 1.0 / n);
  double hk[3];
  for (int k = 0; k < 3; ++k) {
    const double target = 40.0 * (double)(1 << k) * tau * n; /* e sinh H - H */
    double hh = myb_asinh(target / e);
    for (int it = 0; it < 60; ++it) {
      hh -= (e * myb_sinh(hh) - hh - target) / (e * myb_cosh(hh) - 1.0);
    }
    hk[k] = hh;
  }
  /* the asymptotic velocities and the jump of I2 */
  double xp[3];
  double vp[3];
  double xm[3];
  double vm[3];
  double dummy = 0.0;
  hyperbola(e, vinf, 60.0, xp, vp, &dummy);
  hyperbola(e, vinf, -60.0, xm, vm, &dummy);
  double di2[3][3];
  const double vp2 = dot3(vp, vp);
  const double vm2 = dot3(vm, vm);
  for (int i = 0; i < 3; ++i) {
    for (int jj = 0; jj < 3; ++jj) {
      const double delta = i == jj ? 1.0 : 0.0;
      di2[i][jj] =
          2.0 * nu * ((vp[i] * vp[jj] - delta * vp2 / 3.0) - (vm[i] * vm[jj] - delta * vm2 / 3.0));
    }
  }
  double tails[6] = {0.0}; /* times 1 / T */
  for (int s = 0; s < 2; ++s) {
    const double* vv = s == 0 ? vp : vm;
    double g[3];
    for (int i = 0; i < 3; ++i) {
      g[i] = 24.0 / 5.0 * (di2[i][0] * vv[0] + di2[i][1] * vv[1] + di2[i][2] * vv[2]);
    }
    double vg[3];
    cross3(vv, g, vg);
    tails[1] += nu * vg[2];
    double vvg[3];
    cross3(vv, vg, vvg);
    const double sg = s == 0 ? 1.0 : -1.0;
    tails[3] += sg * vvg[0];
    tails[4] += sg * vvg[1];
  }
  double part[3][6] = {{0.0}};
  double prev_h = 0.0;
  double acc[6] = {0.0};
  for (int k = 0; k < 3; ++k) {
    /* from prev_h to hk[k] on both sides */
    const int panels = (int)(8.0 * (hk[k] - prev_h)) + 1;
    const double width = (hk[k] - prev_h) / (double)panels;
    for (int pnl = 0; pnl < panels; ++pnl) {
      const double c = prev_h + width * ((double)pnl + 0.5);
      for (int gq = 0; gq < GL; ++gq) {
        for (int sd = 0; sd < 2; ++sd) {
          const double hh = (sd == 0 ? 1.0 : -1.0) * (c + 0.5 * width * xg[gq]);
          double rt[6];
          encounter_rates(nu, e, vinf, hh, rt);
          for (int q = 0; q < 6; ++q) {
            acc[q] += 0.5 * width * wg[gq] * rt[q];
          }
        }
      }
    }
    prev_h = hk[k];
    const double tk = (e * myb_sinh(hk[k]) - hk[k]) / n;
    for (int q = 0; q < 6; ++q) {
      part[k][q] = acc[q] + tails[q] / tk;
    }
  }
  /* I(T) = I + (a + b ln T) / T^2 at T, 2T, 4T */
  double lim[6];
  for (int q = 0; q < 6; ++q) {
    double tk[3];
    for (int k = 0; k < 3; ++k) {
      tk[k] = (e * myb_sinh(hk[k]) - hk[k]) / n;
    }
    /* eliminate: rows [1, 1/T^2, ln T / T^2] */
    double m[3][4];
    for (int k = 0; k < 3; ++k) {
      m[k][0] = 1.0;
      m[k][1] = 1.0 / (tk[k] * tk[k]);
      m[k][2] = myb_log(tk[k]) / (tk[k] * tk[k]);
      m[k][3] = part[k][q];
    }
    for (int c = 0; c < 3; ++c) { /* Gauss-Jordan on a 3 x 4 system */
      const double pv = m[c][c];
      for (int cc = c; cc < 4; ++cc) {
        m[c][cc] /= pv;
      }
      for (int r = 0; r < 3; ++r) {
        if (r != c) {
          const double f = m[r][c];
          for (int cc = c; cc < 4; ++cc) {
            m[r][cc] -= f * m[c][cc];
          }
        }
      }
    }
    lim[q] = m[0][3];
  }
  encounter_result res = {.de = lim[0], .dj = lim[1], .de_cons = lim[2]};
  /* the periapsis direction is +x, the motion counterclockwise: d chi = dA_y / e (the change of
   * e, dA_x, moves the outgoing asymptote by dA_x / (e sqrt(e^2 - 1))) */
  res.dchi = lim[4] / e + lim[3] / (e * sqrt(e * e - 1.0));
  return res;
}

static void test_encounters(void) {
  /* reference values from the published Fourier transforms (Bini and Damour 2017, eq.
   * rel_x2y2xx_final, with the factor 1/e^2 of (y^2)_w removed: as printed, (y^2)_w has e^-2 of
   * the norm of y^2 in the time domain, and the corrected transforms reproduce Turner's
   * Delta E_GW, eq. EGW, to 2e-21), and the scattering angle (1/nu) dW_r/dj of their eq.
   * W_fourier with the scale 2 r(t) (scratchpad enc_ref.py, mpmath) */
  static const double cases[3][6] = {
      /* e, vinf, Delta E_tail, Delta J_tail, chi_tail, Delta E_GW (Turner) */
      {1.5, 0.2, 4.89364320055063277e-5, -1.79735901023480376e-3, -1.15223461772237e-3, 0.0},
      {3.0, 0.1, 6.49921949729459704e-11, -1.21243674297879726e-7, 3.75195987788463e-9, 0.0},
      {10.0, 0.1, 8.6000042064855418e-14, -9.122378827144159e-10, 1.46127065285951e-11, 0.0}};
  const double nu = 0.25;
  for (int k = 0; k < 3; ++k) {
    const double e = cases[k][0];
    const double vinf = cases[k][1];
    const encounter_result r = encounter(nu, e, vinf);
    const double ee = r.de / cases[k][2] - 1.0;
    const double ej = fabs(r.dj) / fabs(cases[k][3]) - 1.0;
    const double ec = r.dchi / cases[k][4] - 1.0;
    printf(
        "hyperbolic encounter e %.1f, v_inf %.1f, nu 0.25: energy radiated by the dissipative "
        "tail %.10e (frequency domain %.10e, %.1e), angular momentum %.10e (%.10e, %.1e); "
        "scattering angle of the conservative tail %.8e (%.8e, %.1e); net energy change by "
        "the conservative tail %.1e of %.1e\n",
        e, vinf, r.de, cases[k][2], ee, r.dj, cases[k][3], ej, r.dchi, cases[k][4], ec, r.de_cons,
        r.de);
    MT_CHECK(fabs(ee) < 1e-5);
    MT_CHECK(fabs(ej) < 1e-4);
    MT_CHECK(r.dj < 0.0);
    MT_CHECK(fabs(ec) < 1e-4);
  }
  /* large eccentricity: the leading order in 1/e of Bini and Damour (eq. WLO1) with the scale
   * 2 r(t), chi = -(2 pi nu vinf^4 / (5 j^4)) [315/4 + 37 ln(vinf / 2)] (the log(j) terms cancel
   * between W and the change of scale); the difference falls as 1/e */
  double prev = 0.0;
  for (int k = 0; k < 2; ++k) {
    const double e = k == 0 ? 10.0 : 40.0;
    const double vinf = 0.1;
    const double j = sqrt(e * e - 1.0) / vinf;
    const double lo = -(2.0 * PI * nu * vinf * vinf * vinf * vinf / (5.0 * j * j * j * j)) *
                      (315.0 / 4.0 + 37.0 * myb_log(0.5 * vinf));
    const encounter_result r = encounter(nu, e, vinf);
    const double diff = r.dchi / lo - 1.0;
    printf(
        "hyperbolic encounter e %.0f, v_inf 0.1: scattering angle of the conservative tail "
        "%.6e, leading order in 1/e %.6e: ratio - 1 = %.3f (times e: %.2f)\n",
        e, r.dchi, lo, diff, diff * e);
    if (k == 1) {
      MT_CHECK(fabs(diff) < 0.4 * fabs(prev));
    }
    prev = diff;
  }
}

/* the cost of one evaluation (both parts) in the time domain, printed only */
static void test_timing(void) {
  /* e, true anomaly, time domain (1) or Fourier series (0) */
  static const double cases[9][3] = {{0.5, 0.0, 1.0},        {0.9, 0.0, 1.0}, {0.995, 1.0, 1.0},
                                     {1.0 - 1e-6, 0.6, 1.0}, {1.5, 0.0, 1.0}, {3.0, 1.2, 1.0},
                                     {0.5, 0.0, 0.0},        {0.7, 0.0, 0.0}, {0.82, 0.0, 0.0}};
  for (int k = 0; k < 9; ++k) {
    double x[3];
    double v[3];
    conic_state(20.0, cases[k][0], cases[k][1], x, v);
    const int td = cases[k][2] > 0.5;
    double c[3];
    double d[3];
    const clock_t t0 = clock();
    int reps = 0;
    for (; reps < 200; ++reps) {
      MT_CHECK(myb_internal_pn_tail(0.25, x, v, td, c, d) == MYB_OK);
      if ((double)(clock() - t0) > 0.3 * (double)CLOCKS_PER_SEC) {
        ++reps;
        break;
      }
    }
    const double ms = 1e3 * (double)(clock() - t0) / (double)CLOCKS_PER_SEC / (double)reps;
    printf("%s, e %.6f, q 20, true anomaly %.1f: %.2f ms per evaluation\n",
           td ? "time domain" : "Fourier series", cases[k][0], cases[k][1], ms);
  }
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_timing();
  test_ellipses();
  test_switch();
  test_parabolic();
  test_reference();
  test_encounters();
  return mt_finish("test_pn_tail_open");
}
