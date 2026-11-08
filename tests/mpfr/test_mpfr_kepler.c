/* SPDX-License-Identifier: BSD-3-Clause */
/* Kepler's equation and two-body propagation against high-precision references in MPFR.
 *  - Roots E, H, D of the elliptic, hyperbolic and parabolic equations are compared in ulps with
 *    the exact root of the same (exact double) input, found by Newton's method at 300 bits.
 *  - Propagation is compared with an independent universal-variable propagation at 400 bits
 *    (bracketing + bisection/Newton in MPFR) of the same double input state.
 *  - The state transition matrix is compared with central differences of that propagation. */
#include <math.h>
#include <mpfr.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/kepler.h"

enum { ROOT_PREC = 300, PROP_PREC = 400 };

/* |x - exact| in ulps of RN(exact). */
static double ulps(double x, mpfr_srcptr exact) {
  mpfr_t d;
  mpfr_init2(d, ROOT_PREC);
  mpfr_set_d(d, x, MPFR_RNDN);
  mpfr_sub(d, d, exact, MPFR_RNDN);
  const double r = mpfr_get_d(exact, MPFR_RNDN);
  const double ulp = fabs(nextafter(fabs(r), (double)INFINITY) - fabs(r));
  const double e = fabs(mpfr_get_d(d, MPFR_RNDN)) / ulp;
  mpfr_clear(d);
  return e;
}

typedef enum conic { ELLIPSE, HYPERBOLA, PARABOLA } conic;

/* Newton at ROOT_PREC for f(x) = x - e sin x - m, e sinh x - x - m, x + x^3/3 - m. */
static void exact_root(conic kind, double m, double e, double start, mpfr_ptr x) {
  mpfr_t f;
  mpfr_t df;
  mpfr_t t;
  mpfr_inits2(ROOT_PREC, f, df, t, (mpfr_ptr)0);
  mpfr_set_d(x, start, MPFR_RNDN);
  for (int i = 0; i < 40; ++i) {
    if (kind == ELLIPSE) {
      mpfr_sin(t, x, MPFR_RNDN);
      mpfr_mul_d(t, t, e, MPFR_RNDN);
      mpfr_sub(f, x, t, MPFR_RNDN);
      mpfr_sub_d(f, f, m, MPFR_RNDN);
      mpfr_cos(t, x, MPFR_RNDN);
      mpfr_mul_d(t, t, e, MPFR_RNDN);
      mpfr_d_sub(df, 1.0, t, MPFR_RNDN);
    } else if (kind == HYPERBOLA) {
      mpfr_sinh(t, x, MPFR_RNDN);
      mpfr_mul_d(t, t, e, MPFR_RNDN);
      mpfr_sub(f, t, x, MPFR_RNDN);
      mpfr_sub_d(f, f, m, MPFR_RNDN);
      mpfr_cosh(t, x, MPFR_RNDN);
      mpfr_mul_d(t, t, e, MPFR_RNDN);
      mpfr_sub_ui(df, t, 1, MPFR_RNDN);
    } else {
      mpfr_pow_ui(t, x, 3, MPFR_RNDN);
      mpfr_div_ui(t, t, 3, MPFR_RNDN);
      mpfr_add(f, x, t, MPFR_RNDN);
      mpfr_sub_d(f, f, m, MPFR_RNDN);
      mpfr_sqr(t, x, MPFR_RNDN);
      mpfr_add_ui(df, t, 1, MPFR_RNDN);
    }
    mpfr_div(t, f, df, MPFR_RNDN);
    mpfr_sub(x, x, t, MPFR_RNDN);
  }
  mpfr_clears(f, df, t, (mpfr_ptr)0);
}

static void test_roots(void) {
  mt_rng rng = {.state = 801};
  mpfr_t x;
  mpfr_init2(x, ROOT_PREC);
  double worst[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 30000; ++i) {
    double e = 0.0;
    switch (i % 3) {
      case 0:
        e = (double)(mt_next(&rng) >> 11) * 0x1p-53;
        break;
      case 1:
        e = 1.0 - ldexp(1.0, -mt_int(&rng, 1, 52));
        break;
      default:
        e = ldexp(1.0, -mt_int(&rng, 1, 60));
        break;
    }
    const double m = (i % 4 == 0) ? mt_double(&rng, -60, -2)
                                  : 3.14159 * (2.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53 - 1.0);
    double ecc = 0.0;
    MT_CHECK(myb_kepler_elliptic(m, e, &ecc) == MYB_OK);
    exact_root(ELLIPSE, m, e, ecc, x);
    const double err = ulps(ecc, x);
    worst[0] = err > worst[0] ? err : worst[0];

    const double eh = (i % 2 == 0) ? 1.0 + ldexp(1.0, -mt_int(&rng, 1, 52))
                                   : 1.0 + fabs(mt_double(&rng, -10, 10));
    const double mh = mt_double(&rng, -60, 40);
    double h = 0.0;
    MT_CHECK(myb_kepler_hyperbolic(mh, eh, &h) == MYB_OK);
    exact_root(HYPERBOLA, mh, eh, h, x);
    const double errh = ulps(h, x);
    worst[1] = errh > worst[1] ? errh : worst[1];

    const double mp = mt_double(&rng, -300, 300);
    const double d = myb_kepler_parabolic(mp);
    exact_root(PARABOLA, mp, 0.0, d, x);
    const double errp = ulps(d, x);
    worst[2] = errp > worst[2] ? errp : worst[2];
  }
  mpfr_clear(x);
  printf("roots vs exact (ulps): elliptic %.3g, hyperbolic %.3g, parabolic %.3g\n", worst[0],
         worst[1], worst[2]);
  MT_CHECK(worst[0] <= 2.0);
  MT_CHECK(worst[1] <= 2.0);
  MT_CHECK(worst[2] <= 2.0);
}

/* ---- Propagation reference ------------------------------------------------------------------ */

typedef struct ref {
  mpfr_t r0;
  mpfr_t sigma0;
  mpfr_t alpha;
  mpfr_t target;
  mpfr_t t1;
  mpfr_t t2;
  mpfr_t z;
  mpfr_t s;
  mpfr_t u0;
  mpfr_t u1;
  mpfr_t u2;
  mpfr_t u3;
} ref;

static void ref_init(ref* k) {
  mpfr_inits2(PROP_PREC, k->r0, k->sigma0, k->alpha, k->target, k->t1, k->t2, k->z, k->s, k->u0,
              k->u1, k->u2, k->u3, (mpfr_ptr)0);
}

static void ref_clear(ref* k) {
  mpfr_clears(k->r0, k->sigma0, k->alpha, k->target, k->t1, k->t2, k->z, k->s, k->u0, k->u1, k->u2,
              k->u3, (mpfr_ptr)0);
}

/* U0..U3 at chi from the closed forms (cancellation is harmless at 400 bits). */
static void ref_u(ref* k, mpfr_srcptr chi) {
  mpfr_sqr(k->z, chi, MPFR_RNDN);
  mpfr_mul(k->z, k->z, k->alpha, MPFR_RNDN); /* z = alpha chi^2 */
  if (mpfr_zero_p(k->alpha)) {
    mpfr_set_ui(k->u0, 1, MPFR_RNDN);
    mpfr_set(k->u1, chi, MPFR_RNDN);
    mpfr_sqr(k->u2, chi, MPFR_RNDN);
    mpfr_div_ui(k->u2, k->u2, 2, MPFR_RNDN);
    mpfr_pow_ui(k->u3, chi, 3, MPFR_RNDN);
    mpfr_div_ui(k->u3, k->u3, 6, MPFR_RNDN);
    return;
  }
  /* With q = sqrt(|alpha|), s = q chi: U0 = cos s, U1 = sin s / q, U2 = (1 - cos s) / alpha,
   * U3 = (s - sin s) / q^3 (hyperbolic functions and signs for alpha < 0). */
  const int ellipse = mpfr_sgn(k->alpha) > 0;
  mpfr_abs(k->t1, k->alpha, MPFR_RNDN);
  mpfr_sqrt(k->t1, k->t1, MPFR_RNDN); /* q */
  mpfr_mul(k->s, k->t1, chi, MPFR_RNDN);
  if (ellipse) {
    mpfr_cos(k->u0, k->s, MPFR_RNDN);
    mpfr_sin(k->t2, k->s, MPFR_RNDN);
  } else {
    mpfr_cosh(k->u0, k->s, MPFR_RNDN);
    mpfr_sinh(k->t2, k->s, MPFR_RNDN);
  }
  mpfr_div(k->u1, k->t2, k->t1, MPFR_RNDN);
  mpfr_ui_sub(k->u2, 1, k->u0, MPFR_RNDN);
  mpfr_div(k->u2, k->u2, k->alpha, MPFR_RNDN);
  if (ellipse) {
    mpfr_sub(k->u3, k->s, k->t2, MPFR_RNDN);
  } else {
    mpfr_sub(k->u3, k->t2, k->s, MPFR_RNDN);
  }
  mpfr_pow_ui(k->t1, k->t1, 3, MPFR_RNDN);
  mpfr_div(k->u3, k->u3, k->t1, MPFR_RNDN);
}

/* F(chi) into out (uses U values just computed). */
static void ref_f(ref* k, mpfr_ptr out) {
  mpfr_mul(out, k->r0, k->u1, MPFR_RNDN);
  mpfr_mul(k->t1, k->sigma0, k->u2, MPFR_RNDN);
  mpfr_add(out, out, k->t1, MPFR_RNDN);
  mpfr_add(out, out, k->u3, MPFR_RNDN);
  mpfr_sub(out, out, k->target, MPFR_RNDN);
}

static void ref_radius(ref* k, mpfr_ptr out) {
  mpfr_mul(out, k->r0, k->u0, MPFR_RNDN);
  mpfr_mul(k->t1, k->sigma0, k->u1, MPFR_RNDN);
  mpfr_add(out, out, k->t1, MPFR_RNDN);
  mpfr_add(out, out, k->u2, MPFR_RNDN);
}

static void dot_m(mpfr_ptr out, mpfr_t a[3], mpfr_t b[3], mpfr_ptr tmp) {
  mpfr_mul(out, a[0], b[0], MPFR_RNDN);
  for (int i = 1; i < 3; ++i) {
    mpfr_mul(tmp, a[i], b[i], MPFR_RNDN);
    mpfr_add(out, out, tmp, MPFR_RNDN);
  }
}

/* Reference state y = (r, v) after dt from the MPFR state (r0, v0); returns 0 if the bracket
 * search fails (overflow). */
static int reference_core(double mu, mpfr_t r0v[3], mpfr_t v0v[3], double dt, mpfr_t y[6]) {
  ref k;
  ref_init(&k);
  mpfr_t mu_m;
  mpfr_t sqrt_mu;
  mpfr_t lo;
  mpfr_t hi;
  mpfr_t chi;
  mpfr_t f;
  mpfr_t r;
  mpfr_t tmp;
  mpfr_inits2(PROP_PREC, mu_m, sqrt_mu, lo, hi, chi, f, r, tmp, (mpfr_ptr)0);
  mpfr_set_d(mu_m, mu, MPFR_RNDN);
  mpfr_sqrt(sqrt_mu, mu_m, MPFR_RNDN);
  dot_m(k.r0, r0v, r0v, tmp);
  mpfr_sqrt(k.r0, k.r0, MPFR_RNDN);
  dot_m(k.sigma0, r0v, v0v, tmp);
  mpfr_div(k.sigma0, k.sigma0, sqrt_mu, MPFR_RNDN);
  dot_m(k.alpha, v0v, v0v, tmp);
  mpfr_div(k.alpha, k.alpha, mu_m, MPFR_RNDN);
  mpfr_ui_div(tmp, 2, k.r0, MPFR_RNDN);
  mpfr_sub(k.alpha, tmp, k.alpha, MPFR_RNDN);
  mpfr_mul_d(k.target, sqrt_mu, dt, MPFR_RNDN);

  /* Bracket by doubling, then bisection-safeguarded Newton at high precision. */
  const int sign = dt > 0.0 ? 1 : -1;
  mpfr_set_zero(lo, 1);
  mpfr_abs(hi, k.target, MPFR_RNDN);
  mpfr_div(hi, hi, k.r0, MPFR_RNDN);
  if (sign < 0) {
    mpfr_neg(hi, hi, MPFR_RNDN);
  }
  int found = 0;
  for (int i = 0; i < 3000 && !found; ++i) {
    ref_u(&k, hi);
    ref_f(&k, f);
    if (!mpfr_number_p(f)) {
      break;
    }
    if (mpfr_sgn(f) * sign > 0) {
      found = 1;
    } else {
      mpfr_set(lo, hi, MPFR_RNDN);
      mpfr_mul_ui(hi, hi, 2, MPFR_RNDN);
    }
  }
  if (!found) {
    mpfr_clears(mu_m, sqrt_mu, lo, hi, chi, f, r, tmp, (mpfr_ptr)0);
    ref_clear(&k);
    return 0;
  }
  if (sign < 0) {
    mpfr_swap(lo, hi);
  }
  mpfr_add(chi, lo, hi, MPFR_RNDN);
  mpfr_div_ui(chi, chi, 2, MPFR_RNDN);
  for (int i = 0; i < 2000; ++i) {
    ref_u(&k, chi);
    ref_f(&k, f);
    if (mpfr_sgn(f) < 0) {
      mpfr_set(lo, chi, MPFR_RNDN);
    } else {
      mpfr_set(hi, chi, MPFR_RNDN);
    }
    ref_radius(&k, r);
    mpfr_div(tmp, f, r, MPFR_RNDN);
    mpfr_sub(tmp, chi, tmp, MPFR_RNDN);
    if (mpfr_cmp(tmp, lo) <= 0 || mpfr_cmp(tmp, hi) >= 0) {
      mpfr_add(tmp, lo, hi, MPFR_RNDN);
      mpfr_div_ui(tmp, tmp, 2, MPFR_RNDN);
    }
    mpfr_sub(f, tmp, chi, MPFR_RNDN);
    mpfr_set(chi, tmp, MPFR_RNDN);
    if (mpfr_zero_p(f) || mpfr_get_exp(f) < mpfr_get_exp(chi) - (PROP_PREC - 40)) {
      break;
    }
  }
  ref_u(&k, chi);
  ref_radius(&k, r);
  /* f = 1 - U2/r0, g = (r0 U1 + sigma0 U2)/sqrt(mu), fdot = -sqrt(mu) U1/(r r0), gdot = 1 - U2/r */
  mpfr_t fv;
  mpfr_t gv;
  mpfr_t fd;
  mpfr_t gd;
  mpfr_inits2(PROP_PREC, fv, gv, fd, gd, (mpfr_ptr)0);
  mpfr_div(fv, k.u2, k.r0, MPFR_RNDN);
  mpfr_ui_sub(fv, 1, fv, MPFR_RNDN);
  mpfr_mul(gv, k.r0, k.u1, MPFR_RNDN);
  mpfr_mul(tmp, k.sigma0, k.u2, MPFR_RNDN);
  mpfr_add(gv, gv, tmp, MPFR_RNDN);
  mpfr_div(gv, gv, sqrt_mu, MPFR_RNDN);
  mpfr_mul(fd, sqrt_mu, k.u1, MPFR_RNDN);
  mpfr_div(fd, fd, r, MPFR_RNDN);
  mpfr_div(fd, fd, k.r0, MPFR_RNDN);
  mpfr_neg(fd, fd, MPFR_RNDN);
  mpfr_div(gd, k.u2, r, MPFR_RNDN);
  mpfr_ui_sub(gd, 1, gd, MPFR_RNDN);
  for (int i = 0; i < 3; ++i) {
    mpfr_mul(y[i], fv, r0v[i], MPFR_RNDN);
    mpfr_mul(f, gv, v0v[i], MPFR_RNDN);
    mpfr_add(y[i], y[i], f, MPFR_RNDN);
    mpfr_mul(y[3 + i], fd, r0v[i], MPFR_RNDN);
    mpfr_mul(f, gd, v0v[i], MPFR_RNDN);
    mpfr_add(y[3 + i], y[3 + i], f, MPFR_RNDN);
  }
  mpfr_clears(fv, gv, fd, gd, (mpfr_ptr)0);
  mpfr_clears(mu_m, sqrt_mu, lo, hi, chi, f, r, tmp, (mpfr_ptr)0);
  ref_clear(&k);
  return 1;
}

/* Reference state after dt from a double state, rounded to doubles. */
static int reference_propagate(double mu, const myb_pv* s0, double dt, double r_out[3],
                               double v_out[3]) {
  mpfr_t r0v[3];
  mpfr_t v0v[3];
  mpfr_t y[6];
  for (int i = 0; i < 3; ++i) {
    mpfr_init2(r0v[i], PROP_PREC);
    mpfr_init2(v0v[i], PROP_PREC);
    mpfr_set_d(r0v[i], s0->r[i], MPFR_RNDN);
    mpfr_set_d(v0v[i], s0->v[i], MPFR_RNDN);
  }
  for (int i = 0; i < 6; ++i) {
    mpfr_init2(y[i], PROP_PREC);
  }
  const int ok = reference_core(mu, r0v, v0v, dt, y);
  for (int i = 0; ok && i < 3; ++i) {
    r_out[i] = mpfr_get_d(y[i], MPFR_RNDN);
    v_out[i] = mpfr_get_d(y[3 + i], MPFR_RNDN);
  }
  for (int i = 0; i < 3; ++i) {
    mpfr_clear(r0v[i]);
    mpfr_clear(v0v[i]);
  }
  for (int i = 0; i < 6; ++i) {
    mpfr_clear(y[i]);
  }
  return ok;
}

/* Reference state transition matrix: central differences of the 400-bit propagation with
 * steps of 2^-100 of the position and velocity scales (truncation ~2^-200, cancellation
 * leaves ~300 bits). */
static int reference_stm(double mu, const myb_pv* s0, double dt, double stm[36]) {
  mpfr_t r0v[3];
  mpfr_t v0v[3];
  mpfr_t plus[6];
  mpfr_t minus[6];
  mpfr_t h;
  mpfr_init2(h, PROP_PREC);
  for (int i = 0; i < 3; ++i) {
    mpfr_init2(r0v[i], PROP_PREC);
    mpfr_init2(v0v[i], PROP_PREC);
  }
  for (int i = 0; i < 6; ++i) {
    mpfr_init2(plus[i], PROP_PREC);
    mpfr_init2(minus[i], PROP_PREC);
  }
  const double r_scale = sqrt(s0->r[0] * s0->r[0] + s0->r[1] * s0->r[1] + s0->r[2] * s0->r[2]);
  const double v_scale = sqrt(s0->v[0] * s0->v[0] + s0->v[1] * s0->v[1] + s0->v[2] * s0->v[2]);
  int ok = 1;
  for (int j = 0; ok && j < 6; ++j) {
    mpfr_set_d(h, (j < 3 ? r_scale : v_scale) * 0x1p-100, MPFR_RNDN);
    for (int side = 0; ok && side < 2; ++side) {
      for (int i = 0; i < 3; ++i) {
        mpfr_set_d(r0v[i], s0->r[i], MPFR_RNDN);
        mpfr_set_d(v0v[i], s0->v[i], MPFR_RNDN);
      }
      mpfr_ptr component = j < 3 ? r0v[j] : v0v[j - 3];
      if (side == 0) {
        mpfr_add(component, component, h, MPFR_RNDN);
      } else {
        mpfr_sub(component, component, h, MPFR_RNDN);
      }
      ok = reference_core(mu, r0v, v0v, dt, side == 0 ? plus : minus);
    }
    for (int i = 0; ok && i < 6; ++i) {
      mpfr_sub(plus[i], plus[i], minus[i], MPFR_RNDN);
      mpfr_div(plus[i], plus[i], h, MPFR_RNDN);
      mpfr_div_ui(plus[i], plus[i], 2, MPFR_RNDN);
      stm[6 * i + j] = mpfr_get_d(plus[i], MPFR_RNDN);
    }
  }
  mpfr_clear(h);
  for (int i = 0; i < 3; ++i) {
    mpfr_clear(r0v[i]);
    mpfr_clear(v0v[i]);
  }
  for (int i = 0; i < 6; ++i) {
    mpfr_clear(plus[i]);
    mpfr_clear(minus[i]);
  }
  return ok;
}

static double rel3(const double a[3], const double b[3]) {
  const double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) /
         sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
}

static void test_propagation(void) {
  static const char* names[4] = {"ellipse", "hyperbola", "near-parabola", "near-circle"};
  const double mu = 3.986004418e14;
  mt_rng rng = {.state = 802};
  double worst[4] = {0.0, 0.0, 0.0, 0.0};
  long count[4] = {0, 0, 0, 0};
  for (int i = 0; i < 2000; ++i) {
    const int kind = i % 4;
    myb_elements el = {};
    el.p = 6.6e6 * (1.0 + 10.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53);
    if (kind == 0) {
      el.e = 0.95 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    } else if (kind == 1) {
      el.e = 1.0 + fabs(mt_double(&rng, -8, 3));
    } else if (kind == 2) {
      el.e = 1.0 + mt_double(&rng, -40, -20);
    } else {
      el.e = fabs(mt_double(&rng, -40, -10));
    }
    el.i = 3.14 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    el.raan = 6.28 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    el.argp = 6.28 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    const double limit = el.e > 1.0 ? 0.9 * acos(-1.0 / el.e) : (kind == 2 ? 2.5 : 3.14);
    el.nu = limit * (2.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53 - 1.0);
    myb_pv s0 = {};
    if (myb_pv_from_elements(mu, &el, &s0) != MYB_OK) {
      continue;
    }
    /* Up to ~1 period for ellipses: the reference does not reduce periods. */
    const double dt = mt_double(&rng, 0, 14);
    myb_pv ours = {};
    double r_ref[3];
    double v_ref[3];
    if (myb_kepler_propagate(mu, &s0, dt, &ours) != MYB_OK ||
        !reference_propagate(mu, &s0, dt, r_ref, v_ref)) {
      continue;
    }
    const double e = fmax(rel3(ours.r, r_ref), rel3(ours.v, v_ref));
    worst[kind] = e > worst[kind] ? e : worst[kind];
    count[kind]++;
  }
  for (int k = 0; k < 4; ++k) {
    printf("propagation %-13s %5ld cases, max relative error %.3g (MPFR, 400 bits)\n", names[k],
           count[k], worst[k]);
    MT_CHECK(count[k] > 0);
  }
  MT_CHECK(worst[0] < 1e-14 && worst[1] < 1e-14 && worst[2] < 1e-14 && worst[3] < 1e-14);
  /* Ellipses over 1 to 30 periods: the removed periods must not add a phase error beyond the
   * conditioning of the input (the period is computed from alpha in double-word arithmetic). */
  mt_rng turns_rng = {.state = 1930};
  double worst_turns = 0.0;
  long turns_count = 0;
  for (int i = 0; i < 200; ++i) {
    myb_elements el = {};
    el.p = 6.6e6 * (1.0 + 10.0 * (double)(mt_next(&turns_rng) >> 11) * 0x1p-53);
    el.e = 0.95 * (double)(mt_next(&turns_rng) >> 11) * 0x1p-53;
    el.i = 3.14 * (double)(mt_next(&turns_rng) >> 11) * 0x1p-53;
    el.raan = 6.28 * (double)(mt_next(&turns_rng) >> 11) * 0x1p-53;
    el.argp = 6.28 * (double)(mt_next(&turns_rng) >> 11) * 0x1p-53;
    el.nu = 3.14 * (2.0 * (double)(mt_next(&turns_rng) >> 11) * 0x1p-53 - 1.0);
    myb_pv s0 = {};
    if (myb_pv_from_elements(mu, &el, &s0) != MYB_OK) {
      continue;
    }
    const double a = el.p / (1.0 - el.e * el.e);
    const double turns = 1.0 + 29.0 * (double)(mt_next(&turns_rng) >> 11) * 0x1p-53;
    const double dt = turns * 2.0 * acos(-1.0) * sqrt(a * a * a / mu);
    myb_pv ours = {};
    double r_ref[3];
    double v_ref[3];
    if (myb_kepler_propagate(mu, &s0, dt, &ours) != MYB_OK ||
        !reference_propagate(mu, &s0, dt, r_ref, v_ref)) {
      continue;
    }
    const double e = fmax(rel3(ours.r, r_ref), rel3(ours.v, v_ref));
    worst_turns = e > worst_turns ? e : worst_turns;
    turns_count++;
  }
  printf("propagation 1-30 turns    %5ld cases, max relative error %.3g (MPFR, 400 bits)\n",
         turns_count, worst_turns);
  MT_CHECK(turns_count > 100);
  MT_CHECK(worst_turns < 1e-14);
}

/* Largest difference of two state transition matrices in the scales of the initial state
 * (positions in |r0|, velocities in |v0|), relative to the largest scaled entry of ref. */
static double stm_error(const double ours[36], const double ref[36], const myb_pv* s0) {
  const double rs = sqrt(s0->r[0] * s0->r[0] + s0->r[1] * s0->r[1] + s0->r[2] * s0->r[2]);
  const double vs = sqrt(s0->v[0] * s0->v[0] + s0->v[1] * s0->v[1] + s0->v[2] * s0->v[2]);
  double diff = 0.0;
  double size = 0.0;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      const double scale = (j < 3 ? rs : vs) / (i < 3 ? rs : vs);
      diff = fmax(diff, fabs(ours[6 * i + j] - ref[6 * i + j]) * scale);
      size = fmax(size, fabs(ref[6 * i + j]) * scale);
    }
  }
  return diff / size;
}

static void test_stm(void) {
  static const char* names[5] = {"ellipse", "hyperbola", "near-parabola", "near-circle",
                                 "ellipse, 1-30 turns"};
  const double mu = 3.986004418e14;
  mt_rng rng = {.state = 1985};
  double worst[5] = {};
  long count[5] = {};
  for (int i = 0; i < 60; ++i) {
    const int kind = i % 5;
    myb_elements el = {};
    el.p = 6.6e6 * (1.0 + 10.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53);
    if (kind == 0 || kind == 4) {
      el.e = 0.95 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    } else if (kind == 1) {
      el.e = 1.0 + fabs(mt_double(&rng, -8, 3));
    } else if (kind == 2) {
      el.e = 1.0 + mt_double(&rng, -40, -20);
    } else {
      el.e = fabs(mt_double(&rng, -40, -10));
    }
    el.i = 3.14 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    el.raan = 6.28 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    el.argp = 6.28 * (double)(mt_next(&rng) >> 11) * 0x1p-53;
    const double limit = el.e > 1.0 ? 0.9 * acos(-1.0 / el.e) : (kind == 2 ? 2.5 : 3.14);
    el.nu = limit * (2.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53 - 1.0);
    myb_pv s0 = {};
    if (myb_pv_from_elements(mu, &el, &s0) != MYB_OK) {
      continue;
    }
    double dt = mt_double(&rng, 0, 14);
    if (kind == 4) {
      const double a = el.p / (1.0 - el.e * el.e);
      const double period = 2.0 * acos(-1.0) * sqrt(a * a * a / mu);
      dt = (1.0 + 29.0 * (double)(mt_next(&rng) >> 11) * 0x1p-53) * period;
    }
    myb_pv ours = {};
    myb_pv propagated = {};
    double stm[36];
    double ref[36];
    if (myb_kepler_stm(mu, &s0, dt, &ours, stm) != MYB_OK || !reference_stm(mu, &s0, dt, ref)) {
      continue;
    }
    /* The state is the one of myb_kepler_propagate, bit for bit. */
    MT_CHECK(myb_kepler_propagate(mu, &s0, dt, &propagated) == MYB_OK);
    for (int k = 0; k < 3; ++k) {
      MT_CHECK_SAME(ours.r[k], propagated.r[k]);
      MT_CHECK_SAME(ours.v[k], propagated.v[k]);
    }
    const double e = stm_error(stm, ref, &s0);
    worst[kind] = e > worst[kind] ? e : worst[kind];
    count[kind]++;
  }
  for (int k = 0; k < 5; ++k) {
    printf("state transition matrix %-20s %3ld cases, max relative error %.3g (MPFR differences)\n",
           names[k], count[k], worst[k]);
    MT_CHECK(count[k] > 0);
    MT_CHECK(worst[k] < 1e-14);
  }
}

int main(void) {
  test_roots();
  test_propagation();
  test_stm();
  return mt_finish("test_mpfr_kepler");
}
