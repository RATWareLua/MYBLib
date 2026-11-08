/* SPDX-License-Identifier: BSD-3-Clause */
/* Two-body problem: Kepler's equation, anomalies, elements, propagation. Exact solutions
 * against MPFR are in tests/mpfr/test_mpfr_kepler.c. */
#include <math.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/cmath.h"
#include "myb/kepler.h"

static const double u = 0x1p-53;
static const double mu_earth = 3.986004418e14;

static double eccentricity(mt_rng* rng, int i) {
  switch (i % 4) {
    case 0:
      return (double)(mt_next(rng) >> 11) * 0x1p-53; /* [0, 1) */
    case 1:
      return 1.0 - __builtin_ldexp(1.0, -mt_int(rng, 1, 52)); /* up to 1 - 2^-52 */
    case 2:
      return __builtin_ldexp(1.0, -mt_int(rng, 1, 60)); /* nearly circular */
    default:
      return 0.9 + 0.1 * (double)(mt_next(rng) >> 11) * 0x1p-53;
  }
}

static void test_kepler_equation(void) {
  mt_rng rng = {.state = 701};
  double worst = 0.0;
  for (int i = 0; i < 200000; ++i) {
    const double e = eccentricity(&rng, i);
    const double m = (i % 3 == 0) ? mt_double(&rng, -40, 2) : mt_double(&rng, -2, 12);
    double ecc = 0.0;
    MT_CHECK(myb_kepler_elliptic(m, e, &ecc) == MYB_OK);
    double back = 0.0;
    MT_CHECK(myb_anomaly_mean_from_eccentric(ecc, e, &back) == MYB_OK);
    const double err = fabs(back - m) / (fabs(m) + 1e-300);
    worst = err > worst ? err : worst;
    MT_CHECK(fabs(back - m) <= 8.0 * u * fabs(m) + 1e-300);
  }
  printf("elliptic Kepler: worst relative residual %.3g\n", worst);

  worst = 0.0;
  for (int i = 0; i < 200000; ++i) {
    const double e = (i % 2 == 0) ? 1.0 + __builtin_ldexp(1.0, -mt_int(&rng, 1, 52))
                                  : 1.0 + fabs(mt_double(&rng, -10, 10));
    const double m = mt_double(&rng, -40, 40);
    double h = 0.0;
    MT_CHECK(myb_kepler_hyperbolic(m, e, &h) == MYB_OK);
    double back = 0.0;
    MT_CHECK(myb_anomaly_mean_from_hyperbolic(h, e, &back) == MYB_OK);
    /* Rounding H to a double already moves M by up to (e cosh H - 1) |H| u: the residual is
     * judged against that conditioning (the root itself is checked in ulps against MPFR). */
    const double conditioning = (e * myb_cosh(h) - 1.0) * fabs(h) * u;
    const double err = fabs(back - m) / (8.0 * u * fabs(m) + 2.0 * conditioning);
    worst = err > worst ? err : worst;
    MT_CHECK(err <= 1.0);
  }
  printf("hyperbolic Kepler: worst residual %.3g of its conditioning bound\n", worst);

  for (int i = 0; i < 100000; ++i) {
    const double m = mt_double(&rng, -300, 300);
    const double d = myb_kepler_parabolic(m);
    /* Rounding D to a double moves M by up to (1 + D^2) |D| u; evaluating adds ~4u |M|. */
    MT_CHECK(fabs((d + d * d * d / 3.0) - m) <= (1.0 + d * d) * fabs(d) * u + 4.0 * u * fabs(m));
  }

  double out = 0.0;
  MT_CHECK(myb_kepler_elliptic(1.0, 1.0, &out) == MYB_E_DOMAIN);
  MT_CHECK(myb_kepler_elliptic(1.0, -0.1, &out) == MYB_E_DOMAIN);
  MT_CHECK(myb_kepler_elliptic((double)INFINITY, 0.1, &out) == MYB_E_DOMAIN);
  MT_CHECK(myb_kepler_hyperbolic(1.0, 1.0, &out) == MYB_E_DOMAIN);
  MT_CHECK(myb_kepler_elliptic(0.0, 0.5, &out) == MYB_OK && out == 0.0);
  MT_CHECK(myb_kepler_elliptic(-0.0, 0.5, &out) == MYB_OK);
}

static double angle_diff(double a, double b) {
  const double d = myb_atan2(myb_sin(a - b), myb_cos(a - b));
  return fabs(d);
}

static void test_anomalies(void) {
  mt_rng rng = {.state = 702};
  for (int i = 0; i < 100000; ++i) {
    const double e = eccentricity(&rng, i);
    const double nu = (double)(mt_next(&rng) >> 11) * 0x1p-53 * 6.0 - 3.0;
    double ecc = 0.0;
    double back = 0.0;
    MT_CHECK(myb_anomaly_eccentric_from_true(nu, e, &ecc) == MYB_OK);
    MT_CHECK(myb_anomaly_true_from_eccentric(ecc, e, &back) == MYB_OK);
    MT_CHECK(back > -3.2 && back <= 3.2);
    /* Near e = 1 the map squeezes the anomaly: compare in E, where it is well-conditioned. */
    double again = 0.0;
    MT_CHECK(myb_anomaly_eccentric_from_true(back, e, &again) == MYB_OK);
    MT_CHECK(angle_diff(again, ecc) <= 1e-14);
  }
  for (int i = 0; i < 100000; ++i) {
    const double e = 1.0 + fabs(mt_double(&rng, -20, 5));
    const double limit = myb_acos(-1.0 / e);
    const double nu = limit * ((double)(mt_next(&rng) >> 11) * 0x1p-53 * 1.98 - 0.99);
    double h = 0.0;
    double back = 0.0;
    MT_CHECK(myb_anomaly_hyperbolic_from_true(nu, e, &h) == MYB_OK);
    MT_CHECK(myb_anomaly_true_from_hyperbolic(h, e, &back) == MYB_OK);
    MT_CHECK(fabs(back - nu) <= 1e-12 * (1.0 + fabs(nu)) || fabs(h) > 30.0);
  }
  double h = 0.0;
  MT_CHECK(myb_anomaly_hyperbolic_from_true(3.0, 1.5, &h) == MYB_E_DOMAIN); /* beyond */
}

static double rel_diff3(const double a[3], const double b[3]) {
  const double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) /
         sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
}

static myb_pv random_state(mt_rng* rng, int kind) {
  /* kind 0: ellipse, 1: hyperbola, 2: near-parabola, 3: nearly circular equatorial-ish */
  myb_elements el = {};
  el.p = 6.6e6 * (1.0 + 10.0 * (double)(mt_next(rng) >> 11) * 0x1p-53);
  if (kind == 0) {
    el.e = 0.95 * (double)(mt_next(rng) >> 11) * 0x1p-53;
  } else if (kind == 1) {
    el.e = 1.0 + fabs(mt_double(rng, -8, 4));
  } else if (kind == 2) {
    el.e = 1.0 + mt_double(rng, -40, -20) * ((mt_next(rng) & 1) ? 1.0 : -1.0);
  } else {
    el.e = fabs(mt_double(rng, -40, -10));
  }
  el.i = (kind == 3) ? mt_double(rng, -40, -10) : 3.14 * (double)(mt_next(rng) >> 11) * 0x1p-53;
  el.raan = 6.28 * (double)(mt_next(rng) >> 11) * 0x1p-53;
  el.argp = 6.28 * (double)(mt_next(rng) >> 11) * 0x1p-53;
  /* Near-parabolic ellipses far from the periapsis (1 + e cos nu ~ 1 - e) make (p, e, nu)
   * ill-conditioned: an ulp of nu moves r by r (1 - e)^-1 u. They are kept within 2.5 rad. */
  const double limit = el.e > 1.0 ? 0.9 * myb_acos(-1.0 / el.e) : (kind == 2 ? 2.5 : 3.14);
  el.nu = limit * (2.0 * (double)(mt_next(rng) >> 11) * 0x1p-53 - 1.0);
  myb_pv pv = {};
  MT_CHECK(myb_pv_from_elements(mu_earth, &el, &pv) == MYB_OK);
  return pv;
}

static void test_elements(void) {
  mt_rng rng = {.state = 703};
  double worst = 0.0;
  for (int i = 0; i < 100000; ++i) {
    const myb_pv pv = random_state(&rng, i % 4);
    myb_elements el = {};
    MT_CHECK(myb_elements_from_pv(mu_earth, &pv, &el) == MYB_OK);
    myb_pv back = {};
    MT_CHECK(myb_pv_from_elements(mu_earth, &el, &back) == MYB_OK);
    myb_equinoctial eq = {};
    MT_CHECK(myb_equinoctial_from_pv(mu_earth, &pv, &eq) == MYB_OK);
    myb_pv back_eq = {};
    MT_CHECK(myb_pv_from_equinoctial(mu_earth, &eq, &back_eq) == MYB_OK);
    const double e = fmax(fmax(rel_diff3(back.r, pv.r), rel_diff3(back.v, pv.v)),
                          fmax(rel_diff3(back_eq.r, pv.r), rel_diff3(back_eq.v, pv.v)));
    worst = e > worst ? e : worst;
    MT_CHECK(el.i >= 0.0 && el.i <= 3.1415926535897932 && el.raan >= 0.0 &&
             el.raan < 6.2831853071795865 && el.argp >= 0.0 && el.argp < 6.2831853071795865);
  }
  MT_CHECK(worst < 1e-12);
  printf("element round trips: worst relative error %.3g\n", worst);

  /* Exactly equatorial orbit: raan = 0, argp measured from the x axis. */
  const myb_elements eq_in = {.p = 7e6, .e = 0.1, .i = 0.0, .raan = 0.0, .argp = 0.5, .nu = 0.3};
  myb_pv pv = {};
  MT_CHECK(myb_pv_from_elements(mu_earth, &eq_in, &pv) == MYB_OK);
  pv.r[2] = 0.0;
  pv.v[2] = 0.0;
  myb_elements el = {};
  MT_CHECK(myb_elements_from_pv(mu_earth, &pv, &el) == MYB_OK);
  MT_CHECK(el.i == 0.0 && el.raan == 0.0 && fabs(el.argp - 0.5) < 1e-14 &&
           fabs(el.nu - 0.3) < 1e-14);

  /* Retrograde equatorial: equinoctial elements are singular. */
  pv.v[0] = -pv.v[0];
  pv.v[1] = -pv.v[1];
  myb_equinoctial eq = {};
  MT_CHECK(myb_equinoctial_from_pv(mu_earth, &pv, &eq) == MYB_E_DOMAIN);
  /* Rectilinear motion has no plane. */
  const myb_pv radial = {.r = {7e6, 0.0, 0.0}, .v = {100.0, 0.0, 0.0}};
  MT_CHECK(myb_elements_from_pv(mu_earth, &radial, &el) == MYB_E_DOMAIN);
}

static double energy(const myb_pv* s) {
  const double v2 = s->v[0] * s->v[0] + s->v[1] * s->v[1] + s->v[2] * s->v[2];
  return 0.5 * v2 - mu_earth / sqrt(s->r[0] * s->r[0] + s->r[1] * s->r[1] + s->r[2] * s->r[2]);
}

/* The scale at which the energy is computed: its terms v^2/2 and mu/r. */
static double energy_scale(const myb_pv* s) {
  const double v2 = s->v[0] * s->v[0] + s->v[1] * s->v[1] + s->v[2] * s->v[2];
  return 0.5 * v2 + mu_earth / sqrt(s->r[0] * s->r[0] + s->r[1] * s->r[1] + s->r[2] * s->r[2]);
}

static void momentum(const myb_pv* s, double h[3]) {
  h[0] = s->r[1] * s->v[2] - s->r[2] * s->v[1];
  h[1] = s->r[2] * s->v[0] - s->r[0] * s->v[2];
  h[2] = s->r[0] * s->v[1] - s->r[1] * s->v[0];
}

static void test_propagation(void) {
  /* Vallado, Fundamentals of Astrodynamics, Example 2-4 (KEPLER), km and km/s. */
  const myb_pv v0 = {.r = {1131.340e3, -2282.343e3, 6672.423e3},
                     .v = {-5.64305e3, 4.30333e3, 2.42879e3}};
  myb_pv v1 = {};
  MT_CHECK(myb_kepler_propagate(398600.4418e9, &v0, 40.0 * 60.0, &v1) == MYB_OK);
  const double r_expected[3] = {-4219.7527e3, 4363.0292e3, -3958.7666e3};
  const double v_expected[3] = {3.689866e3, -1.916735e3, -6.112511e3};
  for (int k = 0; k < 3; ++k) {
    MT_CHECK(fabs(v1.r[k] - r_expected[k]) < 0.6);  /* half of the last printed digit, m */
    MT_CHECK(fabs(v1.v[k] - v_expected[k]) < 6e-4); /* m/s */
  }
  printf("Vallado 2-4: r = (%.4f, %.4f, %.4f) km, v = (%.6f, %.6f, %.6f) km/s\n", v1.r[0] / 1e3,
         v1.r[1] / 1e3, v1.r[2] / 1e3, v1.v[0] / 1e3, v1.v[1] / 1e3, v1.v[2] / 1e3);

  mt_rng rng = {.state = 704};
  double worst_invariant = 0.0;
  double worst_reverse = 0.0;
  double worst_compose = 0.0;
  for (int i = 0; i < 30000; ++i) {
    const int kind = i % 4;
    const myb_pv s0 = random_state(&rng, kind);
    const double t1 = mt_double(&rng, 0, 15) * ((mt_next(&rng) & 1) ? 1.0 : -1.0);
    const double t2 = mt_double(&rng, 0, 15) * ((mt_next(&rng) & 1) ? 1.0 : -1.0);
    myb_pv a = {};
    myb_pv b = {};
    myb_pv c = {};
    myb_pv back = {};
    if (myb_kepler_propagate(mu_earth, &s0, t1, &a) != MYB_OK) {
      MT_CHECK(kind == 1 || kind == 2); /* hyperbolas may leave the double range */
      continue;
    }
    MT_CHECK(myb_kepler_propagate(mu_earth, &a, t2, &b) == MYB_OK);
    MT_CHECK(myb_kepler_propagate(mu_earth, &s0, t1 + t2, &c) == MYB_OK);
    MT_CHECK(myb_kepler_propagate(mu_earth, &a, -t1, &back) == MYB_OK);
    double h0[3];
    double h1[3];
    momentum(&s0, h0);
    momentum(&a, h1);
    /* h = r x v is measured on the far state too, where r and v are nearly parallel: its
     * rounding scale is |r||v|, not |h|. */
    const double dh[3] = {h1[0] - h0[0], h1[1] - h0[1], h1[2] - h0[2]};
    const double h_scale = fmax(sqrt(a.r[0] * a.r[0] + a.r[1] * a.r[1] + a.r[2] * a.r[2]) *
                                    sqrt(a.v[0] * a.v[0] + a.v[1] * a.v[1] + a.v[2] * a.v[2]),
                                sqrt(h0[0] * h0[0] + h0[1] * h0[1] + h0[2] * h0[2]));
    const double inv =
        fmax(fabs(energy(&a) - energy(&s0)) / fmax(energy_scale(&a), energy_scale(&s0)),
             sqrt(dh[0] * dh[0] + dh[1] * dh[1] + dh[2] * dh[2]) / h_scale);
    worst_invariant = inv > worst_invariant ? inv : worst_invariant;
    const double rev = rel_diff3(back.r, s0.r);
    worst_reverse = rev > worst_reverse ? rev : worst_reverse;
    const double comp = rel_diff3(b.r, c.r);
    worst_compose = comp > worst_compose ? comp : worst_compose;
  }
  printf("propagation: invariants %.3g, reversibility %.3g, composition %.3g (relative)\n",
         worst_invariant, worst_reverse, worst_compose);
  MT_CHECK(worst_invariant < 1e-12);
  MT_CHECK(worst_reverse < 1e-9);
  MT_CHECK(worst_compose < 1e-9);

  /* Independent route for ellipses: elements -> mean anomaly + n dt -> Kepler -> state. */
  double worst_route = 0.0;
  for (int i = 0; i < 20000; ++i) {
    const myb_pv s0 = random_state(&rng, 0);
    myb_elements el = {};
    MT_CHECK(myb_elements_from_pv(mu_earth, &s0, &el) == MYB_OK);
    if (el.e < 0.01) {
      continue;
    }
    const double a_axis = el.p / (1.0 - el.e * el.e);
    const double n = sqrt(mu_earth / (a_axis * a_axis * a_axis));
    const double dt = mt_double(&rng, 0, 12);
    double ecc0 = 0.0;
    double m0 = 0.0;
    double ecc = 0.0;
    MT_CHECK(myb_anomaly_eccentric_from_true(el.nu, el.e, &ecc0) == MYB_OK);
    MT_CHECK(myb_anomaly_mean_from_eccentric(ecc0, el.e, &m0) == MYB_OK);
    MT_CHECK(myb_kepler_elliptic(m0 + n * dt, el.e, &ecc) == MYB_OK);
    MT_CHECK(myb_anomaly_true_from_eccentric(ecc, el.e, &el.nu) == MYB_OK);
    myb_pv via_elements = {};
    MT_CHECK(myb_pv_from_elements(mu_earth, &el, &via_elements) == MYB_OK);
    myb_pv via_universal = {};
    MT_CHECK(myb_kepler_propagate(mu_earth, &s0, dt, &via_universal) == MYB_OK);
    const double d = rel_diff3(via_universal.r, via_elements.r);
    worst_route = d > worst_route ? d : worst_route;
  }
  printf("universal vs Kepler-equation route: worst relative difference %.3g\n", worst_route);
  MT_CHECK(worst_route < 1e-11);

  myb_pv out = {};
  const myb_pv zero = {};
  MT_CHECK(myb_kepler_propagate(mu_earth, &zero, 10.0, &out) == MYB_E_DOMAIN);
  MT_CHECK(myb_kepler_propagate(-1.0, &v0, 10.0, &out) == MYB_E_DOMAIN);
  MT_CHECK(myb_kepler_propagate(mu_earth, &v0, 0.0, &out) == MYB_OK && out.r[0] == v0.r[0]);
}

/* ---- State transition matrix ------------------------------------------------------------------ */

/* Scaled copy: positions in |r0|, velocities in |v0| (the matrix becomes dimensionless). */
static void scale_stm(const double stm[36], const myb_pv* s0, double out[36]) {
  const double rs = sqrt(s0->r[0] * s0->r[0] + s0->r[1] * s0->r[1] + s0->r[2] * s0->r[2]);
  const double vs = sqrt(s0->v[0] * s0->v[0] + s0->v[1] * s0->v[1] + s0->v[2] * s0->v[2]);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      out[6 * i + j] = stm[6 * i + j] * (j < 3 ? rs : vs) / (i < 3 ? rs : vs);
    }
  }
}

/* max |Phi^T J Phi - J| / max |Phi|^2 (J pairs r_k with v_k; the scaled matrix of scale_stm is
 * symplectic with the same J). */
static double symplectic_defect(const double phi[36]) {
  double defect = 0.0;
  double size = 0.0;
  for (int i = 0; i < 36; ++i) {
    size = fmax(size, fabs(phi[i]));
  }
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double s = 0.0;
      for (int k = 0; k < 3; ++k) {
        s += phi[6 * k + i] * phi[6 * (3 + k) + j] - phi[6 * (3 + k) + i] * phi[6 * k + j];
      }
      const double target = (j == i + 3) ? 1.0 : (i == j + 3 ? -1.0 : 0.0);
      defect = fmax(defect, fabs(s - target));
    }
  }
  return defect / (size * size);
}

static void test_stm(void) {
  mt_rng rng = {.state = 1985};
  double worst_symplectic = 0.0;
  double worst_composition = 0.0;
  for (int i = 0; i < 3000; ++i) {
    myb_elements el = {};
    el.p = 7e6 * (1.0 + 5.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53);
    const int kind = i % 3;
    el.e = kind == 0 ? 0.95 * (double)(mt_next(&rng) >> 11) * 0x1p-53
                     : (kind == 1 ? 1.0 + fabs(mt_double(&rng, -6, 2))
                                  : 1.0 + mt_double(&rng, -40, -12));
    el.i = 3.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    el.raan = 6.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    el.argp = 6.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    const double limit = el.e > 1.0 ? 0.9 * acos(-1.0 / el.e) : 3.0;
    el.nu = limit * (2.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53 - 1.0);
    myb_pv s0 = {};
    MT_CHECK(myb_pv_from_elements(mu_earth, &el, &s0) == MYB_OK);
    const double t1 = mt_double(&rng, 0, 13);
    const double t2 = mt_double(&rng, 0, 13);
    myb_pv s1 = {};
    myb_pv s2 = {};
    myb_pv direct = {};
    double phi1[36];
    double phi2[36];
    double phi_direct[36];
    if (myb_kepler_stm(mu_earth, &s0, t1, &s1, phi1) != MYB_OK ||
        myb_kepler_stm(mu_earth, &s1, t2, &s2, phi2) != MYB_OK ||
        myb_kepler_stm(mu_earth, &s0, t1 + t2, &direct, phi_direct) != MYB_OK) {
      continue;
    }
    double scaled[36];
    scale_stm(phi_direct, &s0, scaled);
    worst_symplectic = fmax(worst_symplectic, symplectic_defect(scaled));
    /* Phi(t1 + t2) = Phi(t2; s1) Phi(t1; s0) */
    double product[36];
    for (int r = 0; r < 6; ++r) {
      for (int c = 0; c < 6; ++c) {
        double sum = 0.0;
        for (int k = 0; k < 6; ++k) {
          sum += phi2[6 * r + k] * phi1[6 * k + c];
        }
        product[6 * r + c] = sum;
      }
    }
    double a[36];
    double b[36];
    scale_stm(product, &s0, a);
    scale_stm(phi_direct, &s0, b);
    double diff = 0.0;
    double size = 0.0;
    for (int k = 0; k < 36; ++k) {
      diff = fmax(diff, fabs(a[k] - b[k]));
      size = fmax(size, fabs(b[k]));
    }
    worst_composition = fmax(worst_composition, diff / size);
  }
  printf("state transition matrix: symplectic defect %.2e, composition %.2e (relative)\n",
         worst_symplectic, worst_composition);
  MT_CHECK(worst_symplectic < 1e-13);
  MT_CHECK(worst_composition < 1e-9); /* the composition carries the conditioning of s1 */
  /* dt = 0: identity and the same state; errors as the propagator. */
  const myb_elements el = {.p = 8e6, .e = 0.1, .i = 0.5, .raan = 0.2, .argp = 0.3, .nu = 0.4};
  myb_pv s0 = {};
  MT_CHECK(myb_pv_from_elements(mu_earth, &el, &s0) == MYB_OK);
  myb_pv out = {};
  double phi[36];
  MT_CHECK(myb_kepler_stm(mu_earth, &s0, 0.0, &out, phi) == MYB_OK);
  for (int k = 0; k < 36; ++k) {
    MT_CHECK(phi[k] == (k % 7 == 0 ? 1.0 : 0.0));
  }
  MT_CHECK(out.r[0] == s0.r[0] && out.v[2] == s0.v[2]);
  MT_CHECK(myb_kepler_stm(-1.0, &s0, 10.0, &out, phi) == MYB_E_DOMAIN);
  const myb_pv origin = {};
  MT_CHECK(myb_kepler_stm(mu_earth, &origin, 10.0, &out, phi) == MYB_E_DOMAIN);
}

int main(void) {
  test_kepler_equation();
  test_anomalies();
  test_elements();
  test_propagation();
  test_stm();
  return mt_finish("test_kepler");
}
