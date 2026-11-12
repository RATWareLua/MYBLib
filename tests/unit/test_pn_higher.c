/* SPDX-License-Identifier: BSD-3-Clause */
/* Compact binaries at 3PN and 3.5PN and with spins (myb_pn_binary_evaluate):
 *   - the coefficients A, B of every order against a machine transcription of the LaTeX source
 *     of Blanchet's review (tests/data/pn/blanchet_coefficients.txt, sympy + mpmath);
 *   - circular orbits: the generalised Kepler law Omega^2(r) through 3PN with its logarithm of
 *     the gauge length r0 (Blanchet, eq. keplerlaw), and the radiation-reaction force on them,
 *     whose 1PN relative correction is -743/336 - 11/4 nu (eq. aieom, without the 4PN tail);
 *   - the inspiral rate dr/dt = -(64/5) nu (G m)^3 / (c^5 r^3) [1 + (-1751/336 - 7/4 nu) G m /
 *     (r c^2)] (eq. rOmdot), from the time an integrated orbit takes to shrink;
 *   - the gauge-invariant K = Omega_phi / Omega_r = 1 + 3x + (27/2 - 7 nu) x^2 + (135/2 +
 *     (-649/4 + 123/32 pi^2) nu + 7 nu^2) x^3 (Damour, Jaranowski, Schafer 2000), measured on
 *     near-circular orbits, and its independence of r0;
 *   - spins: the test-mass limit of the spin-orbit term is the Lense-Thirring acceleration
 *     (myb_rel_lense_thirring); the spin-spin and spin-quadrupole dynamics conserve the energy
 *     and the total angular momentum exactly, the spin-orbit dynamics conserve Kidder's
 *     E_N + E_SO and L_N + L_SO + S1 + S2 to second order in the spins; aligned spins do not
 *     precess; the geodetic precession of pulsar B of PSR J0737-3039A/B against Barker &
 *     O'Connell and the predicted 5.0734 deg/yr (Breton et al. 2008). */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "myb/rel.h"

#define TWO_PI 6.283185307179586
#define PI2 9.8696044010893586188

static const unsigned conservative = MYB_PN_1 | MYB_PN_2 | MYB_PN_3;
static const unsigned all_orbital = MYB_PN_1 | MYB_PN_2 | MYB_PN_25 | MYB_PN_3 | MYB_PN_35;

/* G m = c = 1 with the given nu */
static myb_pn_binary binary_of(double nu, unsigned orders, double r0) {
  const double d = sqrt(1.0 - 4.0 * nu);
  const myb_pn_binary b = {.gm1 = 0.5 * (1.0 + d), .gm2 = 0.5 * (1.0 - d), .c = 1.0,
                           .orders = orders, .r0 = r0};
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

/* ---- the coefficients against the machine transcription ------------------------------------ */

static void test_coefficients(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/pn/blanchet_coefficients.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  static const unsigned flags[5] = {MYB_PN_1, MYB_PN_2, MYB_PN_25, MYB_PN_3, MYB_PN_35};
  char line[1024];
  int rows = 0;
  double worst[5] = {0.0};
  for (int guard = 0; guard < 1000 && fgets(line, sizeof line, f) != NULL; ++guard) {
    if (line[0] == '#') {
      continue;
    }
    double w[15];
    char* p = line;
    int k = 0;
    for (; k < 15; ++k) {
      char* end = NULL;
      w[k] = strtod(p, &end);
      if (end == p) {
        break;
      }
      p = end;
    }
    MT_CHECK(k == 15);
    if (k != 15) {
      break;
    }
    const double r = w[1];
    const double rd = w[3];
    const double vt = w[4];
    const double x[3] = {r, 0.0, 0.0};
    const double v[3] = {rd, vt, 0.0};
    for (int o = 0; o < 5; ++o) {
      myb_pn_binary b = binary_of(w[0], flags[o], w[2]);
      double a[3];
      MT_CHECK(myb_pn_binary_evaluate(&b, x, v, NULL, NULL, a, NULL, NULL) == MYB_OK);
      /* a = -(1 / r^2) [(1 + A) n + B v], n = (1, 0, 0) */
      const double bb = -a[1] * r * r / vt;
      const double aa = -a[0] * r * r - 1.0 - bb * rd;
      const double ea = fabs(aa - w[5 + 2 * o]) / (1.0 + fabs(w[5 + 2 * o]));
      const double eb = fabs(bb - w[6 + 2 * o]) / (1.0 + fabs(w[6 + 2 * o]));
      worst[o] = fmax(worst[o], fmax(ea, eb));
      MT_CHECK(a[2] == 0.0);
    }
    ++rows;
  }
  fclose(f);
  printf("coefficients vs the machine transcription of Blanchet (%d states): 1PN %.1e, 2PN %.1e, "
         "2.5PN %.1e, 3PN %.1e, 3.5PN %.1e\n",
         rows, worst[0], worst[1], worst[2], worst[3], worst[4]);
  MT_CHECK(rows == 40);
  for (int o = 0; o < 5; ++o) {
    MT_CHECK(worst[o] < 1e-14);
  }
}

/* ---- circular orbits ------------------------------------------------------------------------ */

/* v^2 of the circular orbit of radius r of the conservative equations (rdot = 0):
 * Omega^2 r = -a_x, v^2 = Omega^2 r^2 (fixed point; dA/dv^2 = O(1/c^2)). */
static double circular_v2(const myb_pn_binary* b, double r) {
  double v2 = 1.0 / r;
  for (int k = 0; k < 60; ++k) {
    const double x[3] = {r, 0.0, 0.0};
    const double v[3] = {0.0, sqrt(v2), 0.0};
    double a[3];
    MT_CHECK(myb_pn_binary_evaluate(b, x, v, NULL, NULL, a, NULL, NULL) == MYB_OK);
    v2 = -a[0] * r;
  }
  return v2;
}

static void test_kepler_law(void) {
  const double nus[2] = {0.25, 0.1};
  const double logs[2] = {0.0, 2.0};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      const double nu = nus[i];
      const double lr = logs[j];
      double d[3];
      const double g0 = 5e-4;
      for (int k = 0; k < 3; ++k) {
        const double g = g0 * (double)(1 << k);
        const double r = 1.0 / g;
        myb_pn_binary b = binary_of(nu, conservative, r * myb_exp(-lr));
        const double om2r3 = circular_v2(&b, r) * r; /* Omega^2 r^3 / (G m) = v^2 r */
        const double c1 = -3.0 + nu;
        const double c2 = 6.0 + 41.0 / 4.0 * nu + nu * nu;
        d[k] = (om2r3 - 1.0 - c1 * g - c2 * g * g) / (g * g * g);
      }
      /* d = c3 + O(g) + O(g^2): Richardson on g, 2g, 4g */
      const double c3 = (8.0 * d[0] - 6.0 * d[1] + d[2]) / 3.0;
      const double expected = -10.0 + (-75707.0 / 840.0 + 41.0 / 64.0 * PI2 + 22.0 * lr) * nu +
                              9.5 * nu * nu + nu * nu * nu;
      printf("Kepler law, nu %.2f, ln(r/r0) %.0f: 3PN coefficient %.6f (Blanchet %.6f)\n", nu, lr,
             c3, expected);
      MT_CHECK(fabs(c3 - expected) < 1e-3);
    }
  }
}

static void test_circular_radiation_reaction(void) {
  const double nus[2] = {0.25, 0.1};
  for (int i = 0; i < 2; ++i) {
    const double nu = nus[i];
    const double g = 1e-5;
    const double r = 1.0 / g;
    myb_pn_binary b = binary_of(nu, conservative, r);
    const double vt = sqrt(circular_v2(&b, r));
    const double rdot =
        -64.0 / 5.0 * nu * g * g * g * (1.0 + g * (-1751.0 / 336.0 - 1.75 * nu));
    b.orders = all_orbital;
    const double x[3] = {r, 0.0, 0.0};
    const double v[3] = {rdot, vt, 0.0};
    double a[3];
    MT_CHECK(myb_pn_binary_evaluate(&b, x, v, NULL, NULL, a, NULL, NULL) == MYB_OK);
    /* a = -Omega^2 x - (32/5) nu g^4 [1 + g (-743/336 - 11/4 nu)] v */
    const double beta = a[1] / vt;
    const double beta_n = -32.0 / 5.0 * nu * g * g * g * g;
    const double coefficient = (beta / beta_n - 1.0) / g;
    const double expected = -743.0 / 336.0 - 2.75 * nu;
    printf("radiation reaction on a circular orbit, nu %.2f: 1PN coefficient %.6f (Blanchet "
           "%.6f)\n",
           nu, coefficient, expected);
    MT_CHECK(fabs(coefficient - expected) < 1e-3);
  }
}

/* ---- integrations ----------------------------------------------------------------------------- */

typedef struct run {
  myb_pn_binary b;
  double r_event;
} run;

static myb_status orbit_accel(void* context, double t, const double* x, const double* v,
                              double* a) {
  (void)t;
  const run* q = context;
  return myb_pn_binary_evaluate(&q->b, x, v, NULL, NULL, a, NULL, NULL);
}

static double radius_event(void* context, double t, const double* y) {
  (void)t;
  const run* q = context;
  return norm3(y) - q->r_event;
}

static double rdot_event(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return dot3(y, y + 3);
}

static double t_shrink(double c1, double nu, double ra, double rb) {
  /* integral of r^4 / (r + c1) from rb to ra, times 5 / (64 nu) */
  const double fa = ra * ra * ra * ra / 4.0 - c1 * ra * ra * ra / 3.0 + c1 * c1 * ra * ra / 2.0 -
                    c1 * c1 * c1 * ra + c1 * c1 * c1 * c1 * myb_log(ra + c1);
  const double fb = rb * rb * rb * rb / 4.0 - c1 * rb * rb * rb / 3.0 + c1 * c1 * rb * rb / 2.0 -
                    c1 * c1 * c1 * rb + c1 * c1 * c1 * c1 * myb_log(rb + c1);
  return 5.0 / (64.0 * nu) * (fa - fb);
}

static void test_inspiral_rate(void) {
  const double nu = 0.25;
  const double r_start = 200.0;
  const double r_end = 190.0;
  run q = {.b = binary_of(nu, conservative, r_start), .r_event = r_end};
  const double g = 1.0 / r_start;
  const double c1 = -1751.0 / 336.0 - 1.75 * nu;
  const double vt = sqrt(circular_v2(&q.b, r_start));
  const double rdot = -64.0 / 5.0 * nu * g * g * g * (1.0 + c1 * g);
  q.b.orders = all_orbital;
  const double y0[6] = {r_start, 0.0, 0.0, rdot, vt, 0.0};
  const myb_ias15_options o = {.epsilon = 1e-9};
  myb_ias15* d = NULL;
  MT_CHECK(myb_ias15_create(3, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_ias15_start(d, orbit_accel, &q, 0.0, y0) == MYB_OK);
  double t_end = 0.0;
  for (long k = 0; k < 10000000L; ++k) {
    MT_CHECK(myb_ias15_step(d, 1e300) == MYB_OK);
    if (norm3(myb_ias15_state(d)) < r_end) {
      int found = 0;
      double ye[6];
      MT_CHECK(myb_ias15_locate(d, radius_event, &q, &found, &t_end, ye) == MYB_OK);
      MT_CHECK(found);
      break;
    }
  }
  myb_ias15_destroy(d);
  /* the measured time against the formula with and without its 1PN term */
  const double t0 = t_shrink(0.0, nu, r_start, r_end);
  const double t1 = t_shrink(c1, nu, r_start, r_end);
  const double fraction = (t_end - t0) / (t1 - t0);
  printf("inspiral 200 -> 190 G m / c^2, nu 0.25: %.6e (leading order %.6e, with the 1PN term "
         "%.6e): fraction of the 1PN term %.4f\n",
         t_end, t0, t1, fraction);
  /* the remainder is the 2PN correction of dr/dt, O(G m / (r c^2)) of the 1PN term */
  MT_CHECK(fabs(fraction - 1.0) < 0.05);
}

/* K over `periods` radial periods of a near-circular orbit; x from Omega_phi */
static void measure_k(myb_pn_binary* b, double r, int periods, double* k_out, double* x_out) {
  run q = {.b = *b};
  const double vt = sqrt(circular_v2(b, r)) * (1.0 + 1.5e-5);
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
        double a = myb_atan2(ye[1], ye[0]);
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

static void test_periastron_advance_3pn(void) {
  const double nus[2] = {0.25, 0.1};
  for (int i = 0; i < 2; ++i) {
    const double nu = nus[i];
    const double expected = 67.5 + (-649.0 / 4.0 + 123.0 / 32.0 * PI2) * nu + 7.0 * nu * nu;
    double k3[2];
    for (int g = 0; g < 2; ++g) { /* r0 = r and r0 = r / 1000 */
      double res[2];
      for (int s = 0; s < 2; ++s) {
        const double r = s == 0 ? 1000.0 : 2000.0;
        myb_pn_binary b = binary_of(nu, conservative, g == 0 ? r : r / 1000.0);
        double kk = 0.0;
        double x = 0.0;
        measure_k(&b, r, 8, &kk, &x);
        res[s] = (kk - 1.0 - 3.0 * x - (13.5 - 7.0 * nu) * x * x) / (x * x * x);
      }
      /* res = k3 + O(x): extrapolate from x = 1e-3, 5e-4 */
      k3[g] = 2.0 * res[1] - res[0];
      printf("K, nu %.2f, r0 = r%s: (K - 1 - 3x - k2 x^2) / x^3 = %.4f (x 1e-3), %.4f (x 5e-4), "
             "extrapolated %.4f (DJS %.4f)\n",
             nu, g == 0 ? "" : " / 1000", res[0], res[1], k3[g], expected);
      MT_CHECK(fabs(k3[g] - expected) < 0.01 * expected);
    }
    /* r0 changes the coordinates (ln 1000 shifts the Kepler law by 22 nu ln 1000 x^3 ~ 38 x^3)
     * but not K */
    MT_CHECK(fabs(k3[0] - k3[1]) < 0.002 * expected);
  }
}

/* ---- spins ------------------------------------------------------------------------------------ */

typedef struct spin_run {
  myb_pn_binary b;
} spin_run;

/* 9 "positions" (x, p1, p2) whose velocities are (v, chi1, chi2) */
static myb_status spin_accel(void* context, double t, const double* x, const double* v,
                             double* a) {
  (void)t;
  const spin_run* q = context;
  return myb_pn_binary_evaluate(&q->b, x, v, v + 3, v + 6, a, a + 3, a + 6);
}

static void test_lense_thirring_limit(void) {
  const double gm1 = 1.0;
  const double gm2 = 1e-12;
  const myb_pn_binary b = {.gm1 = gm1, .gm2 = gm2, .c = 1.0, .orders = MYB_PN_SO};
  const double chi[3] = {0.3, 0.5, 0.8};
  const double x[3] = {7.0, -3.0, 2.0};
  const double v[3] = {0.1, 0.25, -0.12};
  double a[3];
  double a0[3];
  MT_CHECK(myb_pn_binary_evaluate(&b, x, v, chi, NULL, a, NULL, NULL) == MYB_OK);
  myb_pn_binary bn = b;
  bn.orders = 0;
  MT_CHECK(myb_pn_binary_evaluate(&bn, x, v, chi, NULL, a0, NULL, NULL) == MYB_OK);
  /* angular momentum per unit mass of the central body: S1 / m1 = chi gm1 / c */
  const double j[3] = {chi[0] * gm1, chi[1] * gm1, chi[2] * gm1};
  const myb_ppn ppn = {.c = 1.0, .beta = 1.0, .gamma = 1.0};
  double lt[3];
  MT_CHECK(myb_rel_lense_thirring(&ppn, gm1 + gm2, j, x, v, lt) == MYB_OK);
  double err = 0.0;
  for (int i = 0; i < 3; ++i) {
    err = fmax(err, fabs((a[i] - a0[i]) - lt[i]));
  }
  printf("spin-orbit, test-mass limit vs Lense-Thirring: %.2e of %.2e\n", err, norm3(lt));
  MT_CHECK(err < 1e-10 * norm3(lt));
}

static void test_aligned_spins(void) {
  const double nu = 0.2;
  myb_pn_binary b = binary_of(nu, MYB_PN_ALL, 1.0);
  b.kappa1 = 1.0;
  b.kappa2 = 5.0;
  const double x[3] = {12.0, 5.0, 0.0};
  const double v[3] = {-0.1, 0.22, 0.0};
  const double chi1[3] = {0.0, 0.0, 0.9};
  const double chi2[3] = {0.0, 0.0, -0.5};
  double a[3];
  double d1[3];
  double d2[3];
  MT_CHECK(myb_pn_binary_evaluate(&b, x, v, chi1, chi2, a, d1, d2) == MYB_OK);
  MT_CHECK(a[2] == 0.0);
  MT_CHECK(d1[0] == 0.0 && d1[1] == 0.0 && d1[2] == 0.0);
  MT_CHECK(d2[0] == 0.0 && d2[1] == 0.0 && d2[2] == 0.0);
}

/* Kidder's E_SO, L_SO (G = c = 1, m = 1) and the spin-spin and spin-quadrupole energies */
typedef struct spin_terms {
  double e_n;
  double e_so;
  double e_ss;
  double e_qm;
  double l_n[3];
  double l_so[3];
  double s[3];
} spin_terms;

static spin_terms spin_terms_of(const myb_pn_binary* b, const double* y) {
  const double m1 = b->gm1;
  const double m2 = b->gm2;
  const double mu = m1 * m2;
  const double* x = y;
  const double* v = y + 9;
  const double* c1 = y + 12;
  const double* c2 = y + 15;
  const double r = norm3(x);
  const double n[3] = {x[0] / r, x[1] / r, x[2] / r};
  double s1[3];
  double s2[3];
  double s[3];
  double sd[3]; /* (dm / m) Delta */
  for (int i = 0; i < 3; ++i) {
    s1[i] = c1[i] * m1 * m1;
    s2[i] = c2[i] * m2 * m2;
    s[i] = s1[i] + s2[i];
    sd[i] = (m1 - m2) * (s2[i] / m2 - s1[i] / m1);
  }
  spin_terms t;
  memset(&t, 0, sizeof t);
  t.e_n = mu * (0.5 * dot3(v, v) - 1.0 / r);
  double xv[3];
  cross3(x, v, xv);
  for (int i = 0; i < 3; ++i) {
    t.l_n[i] = mu * xv[i];
    t.s[i] = s[i];
  }
  double ssd[3];
  for (int i = 0; i < 3; ++i) {
    ssd[i] = s[i] + sd[i];
  }
  const double r3 = r * r * r;
  t.e_so = dot3(t.l_n, ssd) / r3;
  t.e_ss = (3.0 * dot3(n, s1) * dot3(n, s2) - dot3(s1, s2)) / r3;
  t.e_qm = 1.5 * m2 * b->kappa1 / (m1 * r3) * (dot3(n, s1) * dot3(n, s1) - dot3(s1, s1) / 3.0) +
           1.5 * m1 * b->kappa2 / (m2 * r3) * (dot3(n, s2) * dot3(n, s2) - dot3(s2, s2) / 3.0);
  /* L_SO = (mu / m) { (m / r) n x [n x (3 S + sd)] - 1/2 v x [v x (S + sd)] } */
  double u3[3];
  double inner[3];
  double outer[3];
  for (int i = 0; i < 3; ++i) {
    u3[i] = 3.0 * s[i] + sd[i];
  }
  cross3(n, u3, inner);
  cross3(n, inner, outer);
  double vin[3];
  double vout[3];
  cross3(v, ssd, vin);
  cross3(v, vin, vout);
  for (int i = 0; i < 3; ++i) {
    t.l_so[i] = mu * (outer[i] / r - 0.5 * vout[i]);
  }
  return t;
}

/* Integrates `orbits` orbits of period `period` and returns the largest change of the energy
 * and of the total angular momentum with and without the spin terms. */
static void spin_conservation(spin_run* q, const double y0[18], double period, int orbits,
                              double* de_with, double* de_without, double* dj_with,
                              double* dj_without) {
  const myb_ias15_options o = {.epsilon = 1e-9, .block = 3};
  myb_ias15* d = NULL;
  MT_CHECK(myb_ias15_create(9, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_ias15_start(d, spin_accel, q, 0.0, y0) == MYB_OK);
  const spin_terms t0 = spin_terms_of(&q->b, y0);
  const int with_so = (q->b.orders & MYB_PN_SO) != 0;
  const int with_ss = (q->b.orders & MYB_PN_SS) != 0;
  const double e_spin0 = (with_so ? t0.e_so : 0.0) + (with_ss ? t0.e_ss + t0.e_qm : 0.0);
  *de_with = *de_without = *dj_with = *dj_without = 0.0;
  for (long k = 0; k < 10000000L && myb_ias15_time(d) < orbits * period; ++k) {
    MT_CHECK(myb_ias15_step(d, orbits * period) == MYB_OK);
    const spin_terms t = spin_terms_of(&q->b, myb_ias15_state(d));
    const double e_spin = (with_so ? t.e_so : 0.0) + (with_ss ? t.e_ss + t.e_qm : 0.0);
    *de_with = fmax(*de_with, fabs(t.e_n + e_spin - t0.e_n - e_spin0));
    *de_without = fmax(*de_without, fabs(t.e_n - t0.e_n));
    double jw[3];
    double jo[3];
    for (int i = 0; i < 3; ++i) {
      const double so = with_so ? t.l_so[i] - t0.l_so[i] : 0.0;
      jw[i] = t.l_n[i] - t0.l_n[i] + so + t.s[i] - t0.s[i];
      jo[i] = t.l_n[i] - t0.l_n[i] + t.s[i] - t0.s[i];
    }
    *dj_with = fmax(*dj_with, norm3(jw));
    *dj_without = fmax(*dj_without, norm3(jo));
  }
  myb_ias15_destroy(d);
}

static void test_spin_conservation(void) {
  /* G m = c = 1, nu = 0.2; an orbit of rp = 30, ra = 70 in the xy plane, spins oblique */
  const double nu = 0.2;
  const double rp = 30.0;
  const double ra = 70.0;
  const double sma = 0.5 * (rp + ra);
  const double vp = sqrt(ra / (sma * rp));
  const double period = TWO_PI * sma * sqrt(sma);
  double y0[18] = {0};
  y0[0] = rp;
  y0[10] = vp;
  const double chi1[3] = {0.3, -0.5, 0.6};
  const double chi2[3] = {0.7, 0.2, -0.4};
  memcpy(y0 + 12, chi1, sizeof chi1);
  memcpy(y0 + 15, chi2, sizeof chi2);
  {
    /* Newton + spin-spin + spin quadrupoles: a Hamiltonian system, E_N + E_SS + E_QM and
     * L_N + S1 + S2 conserved exactly */
    spin_run q = {.b = binary_of(nu, MYB_PN_SS, 1.0)};
    q.b.kappa1 = 1.0;
    q.b.kappa2 = 6.0;
    double y1[18];
    memcpy(y1, y0, sizeof y1);
    y1[0] = 5.0; /* a tight orbit for large spin-spin terms: rp = 5, ra = 11.67 */
    y1[10] = sqrt(1.4 / 5.0);
    double dew;
    double deo;
    double djw;
    double djo;
    spin_conservation(&q, y1, TWO_PI * 8.333 * sqrt(8.333), 3, &dew, &deo, &djw, &djo);
    const spin_terms t = spin_terms_of(&q.b, y1);
    printf("spin-spin + quadrupoles: energy change %.2e (without E_SS, E_QM %.2e) of %.3e; "
           "angular momentum %.2e of %.3e\n",
           dew, deo, fabs(t.e_n), djw, norm3(t.l_n));
    MT_CHECK(dew < 1e-13 * fabs(t.e_n));
    MT_CHECK(deo > 1e-7 * fabs(t.e_n));
    MT_CHECK(djw < 1e-13 * norm3(t.l_n));
  }
  {
    /* Newton + spin-orbit. Kidder's E_N + E_SO is conserved up to terms of second order in the
     * spins (spin-orbit squared, of the order of the spin-spin terms): with the spins scaled by
     * 1/10 its change must fall by 100, while without E_SO by 10. L_N + L_SO + S is conserved at
     * first order in the spins up to the precession of the spins inside L_SO, a term
     * O(G mu / (c^2 r)) smaller than L_SO (beyond the order of the equations): the ratio of
     * the change to the change without L_SO falls as 1/r. */
    spin_run q = {.b = binary_of(nu, MYB_PN_SO, 1.0)};
    double dew[3];
    double deo[3];
    double djw[3];
    double djo[3];
    for (int k = 0; k < 3; ++k) {
      double y1[18];
      memcpy(y1, y0, sizeof y1);
      double p = period;
      if (k == 1) {
        for (int i = 12; i < 18; ++i) {
          y1[i] *= 0.1;
        }
      } else if (k == 2) { /* the same orbit 10 times larger */
        y1[0] *= 10.0;
        y1[10] /= sqrt(10.0);
        p *= 10.0 * sqrt(10.0);
      }
      spin_conservation(&q, y1, p, 2, &dew[k], &deo[k], &djw[k], &djo[k]);
      printf("spin-orbit, %s: energy change %.2e (without E_SO %.2e); angular momentum "
             "change %.2e (without L_SO %.2e)\n",
             k == 0 ? "rp 30, ra 70" : k == 1 ? "spins / 10" : "orbit x 10", dew[k], deo[k],
             djw[k], djo[k]);
    }
    MT_CHECK(dew[0] / dew[1] > 80.0 && dew[0] / dew[1] < 120.0);
    MT_CHECK(deo[0] / deo[1] > 8.0 && deo[0] / deo[1] < 12.0);
    MT_CHECK(dew[1] < 0.01 * deo[1]);
    const double ratio_near = djw[0] / djo[0];
    const double ratio_far = djw[2] / djo[2];
    MT_CHECK(ratio_near < 0.1);
    MT_CHECK(ratio_far / ratio_near > 0.08 && ratio_far / ratio_near < 0.12);
  }
}

static double spin_rdot_event(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return dot3(y, y + 9);
}

static void test_double_pulsar(void) {
  /* PSR J0737-3039A/B (Kramer et al. 2006): m_A = 1.3381, m_B = 1.2489, Pb = 0.10225156248 d,
   * e = 0.0877775; body 1 = pulsar B, whose spin precesses */
  const double gm_sun = 1.32712440041e20;
  const double c = 299792458.0;
  const double ma = 1.3381;
  const double mb = 1.2489;
  const double pb = 0.10225156248 * 86400.0;
  const double e = 0.0877775;
  spin_run q = {.b = {.gm1 = mb * gm_sun, .gm2 = ma * gm_sun, .c = c,
                      .orders = MYB_PN_1 | MYB_PN_2 | MYB_PN_SO}};
  const double gm = q.b.gm1 + q.b.gm2;
  const double n = TWO_PI / pb;
  const double a = myb_cbrt(gm / (n * n));
  double y0[18] = {0};
  y0[0] = a * (1.0 - e);
  y0[10] = sqrt(gm / a * (1.0 + e) / (1.0 - e));
  y0[12] = 1e-4; /* a small spin of B in the orbital plane */
  const myb_ias15_options o = {.epsilon = 1e-9, .block = 3};
  myb_ias15* d = NULL;
  MT_CHECK(myb_ias15_create(9, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_ias15_start(d, spin_accel, &q, 0.0, y0) == MYB_OK);
  enum { ORBITS = 1000 };
  int found_count = 0;
  double prev = spin_rdot_event(NULL, 0.0, y0);
  double t_first = 0.0;
  double t_last = 0.0;
  double phi_first = 0.0;
  double phi_last = 0.0;
  double r_min = 1e300;
  double r_max = 0.0;
  for (long k = 0; k < 10000000L && found_count < ORBITS; ++k) {
    MT_CHECK(myb_ias15_step(d, 1e300) == MYB_OK);
    const double* y = myb_ias15_state(d);
    const double now = spin_rdot_event(NULL, 0.0, y);
    if ((prev < 0.0) != (now < 0.0)) {
      int found = 0;
      double te = 0.0;
      double ye[18];
      MT_CHECK(myb_ias15_locate(d, spin_rdot_event, NULL, &found, &te, ye) == MYB_OK);
      if (found) {
        const double rr = norm3(ye);
        r_min = fmin(r_min, rr);
        r_max = fmax(r_max, rr);
        if (prev < 0.0) { /* periastron: the same orbital phase each time */
          const double phi = myb_atan2(ye[13], ye[12]);
          if (found_count == 0) {
            t_first = te;
            phi_first = phi;
          }
          t_last = te;
          phi_last = phi;
          ++found_count;
        }
      }
    }
    prev = now;
  }
  myb_ias15_destroy(d);
  MT_CHECK(found_count == ORBITS);
  const double year = 365.25 * 86400.0;
  const double deg = 180.0 / 3.141592653589793;
  const double rate = (phi_last - phi_first) / (t_last - t_first) * year * deg;
  /* Barker & O'Connell for the integrated orbit (its radii, the mean motion of the periastra) */
  const double e_orbit = (r_max - r_min) / (r_max + r_min);
  const double a_orbit = 0.5 * (r_max + r_min);
  const double n_orbit = TWO_PI * (found_count - 1) / (t_last - t_first);
  const double gmu = q.b.gm1 * q.b.gm2 / gm;
  const double bo = (2.0 + 1.5 * ma / mb) * gmu / (c * c * a_orbit) * n_orbit /
                    (1.0 - e_orbit * e_orbit) * year * deg;
  printf("PSR J0737-3039B, %d orbits: geodetic precession %.5f deg/yr, Barker-O'Connell %.5f "
         "(ratio %.7f), predicted 5.0734 (Breton et al. 2008)\n",
         ORBITS, rate, bo, rate / bo);
  MT_CHECK(fabs(rate / bo - 1.0) < 1e-4);
  MT_CHECK(fabs(rate - 5.0734) < 0.0010);
}

static void test_domain(void) {
  const double x[3] = {10.0, 0.0, 0.0};
  const double v[3] = {0.0, 0.3, 0.0};
  double a[3];
  myb_pn_binary b = binary_of(0.25, MYB_PN_3, 0.0);
  MT_CHECK(myb_pn_binary_evaluate(&b, x, v, NULL, NULL, a, NULL, NULL) == MYB_E_DOMAIN);
  b.r0 = 1.0;
  MT_CHECK(myb_pn_binary_evaluate(&b, x, v, NULL, NULL, a, NULL, NULL) == MYB_OK);
  b.orders = 1u << 20; /* an unknown bit (128 and above up to 2048 are the higher terms) */
  MT_CHECK(myb_pn_binary_evaluate(&b, x, v, NULL, NULL, a, NULL, NULL) == MYB_E_DOMAIN);
  b.orders = MYB_PN_SS;
  b.kappa1 = -1.0;
  MT_CHECK(myb_pn_binary_evaluate(&b, x, v, NULL, NULL, a, NULL, NULL) == MYB_E_DOMAIN);
  b.kappa1 = 1.0;
  const double bad[3] = {(double)NAN, 0.0, 0.0};
  MT_CHECK(myb_pn_binary_evaluate(&b, x, v, bad, NULL, a, NULL, NULL) == MYB_E_DOMAIN);
  /* the old entry point refuses 3PN (no gauge length) */
  MT_CHECK(myb_pn_binary_acceleration(0.5, 0.5, 1.0, MYB_PN_3, x, v, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_pn_binary_acceleration(0.5, 0.5, 1.0, MYB_PN_1 | MYB_PN_2 | MYB_PN_25, x, v, a) ==
           MYB_OK);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_coefficients();
  test_kepler_law();
  test_circular_radiation_reaction();
  test_inspiral_rate();
  test_periastron_advance_3pn();
  test_lense_thirring_limit();
  test_aligned_spins();
  test_spin_conservation();
  test_double_pulsar();
  test_domain();
  return mt_finish("test_pn_higher");
}
