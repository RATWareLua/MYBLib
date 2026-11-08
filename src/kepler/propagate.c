/* SPDX-License-Identifier: BSD-3-Clause */
/* Two-body propagation in universal variables (Battin 1999, §4.5; Goodyear 1965; Shepperd
 * 1985). With alpha = 2/r0 - v0^2/mu (reciprocal semi-major axis) and z = alpha chi^2:
 *   U0 = c0(z), U1 = chi c1(z), U2 = chi^2 c2(z), U3 = chi^3 c3(z),
 *   F(chi) = r0 U1 + sigma0 U2 + U3 - sqrt(mu) dt,   F'(chi) = r = r0 U0 + sigma0 U1 + U2 > 0.
 * F is increasing, so the root is bracketed (the bracket is found by doubling) and Newton's
 * method is safeguarded by bisection. The Stumpff forms of kepler/stumpff.c keep every U
 * accurate for all conics, near-parabolic ones included. */
#include <math.h> /* sqrt, fabs */

#include "core/base.h"
#include "kepler/stumpff.h"
#include "myb/cmath.h"
#include "myb/kepler.h"
#include "myb/num.h"
#include "num/root.h"
#include "num/vec3.h"

enum { MAX_DOUBLINGS = 2100, MAX_ITERATIONS = 200 };

static const double two_pi_hi = 0x1.921fb54442d18p+2;

static int finite(double x) {
  return x - x == 0.0;
}

typedef struct universal {
  double r0;
  double sigma0;
  double alpha;
  double target; /* sqrt(mu) dt */
} universal;

typedef struct u_values {
  double u0;
  double u1;
  double u2;
  double u3;
} u_values;

static u_values u_at(const universal* k, double chi) {
  const double z = k->alpha * chi * chi;
  const double chi2 = chi * chi;
  return (u_values){.u0 = myb_internal_stumpff_c0(z),
                    .u1 = chi * myb_internal_stumpff_c1(z),
                    .u2 = chi2 * myb_internal_stumpff_c2(z),
                    .u3 = chi2 * chi * myb_internal_stumpff_c3(z)};
}

static double residual(const universal* k, const u_values* u) {
  return k->r0 * u->u1 + k->sigma0 * u->u2 + u->u3 - k->target;
}

static double radius(const universal* k, const u_values* u) {
  return k->r0 * u->u0 + k->sigma0 * u->u1 + u->u2;
}

/* Largest |chi| at which every term of F stays finite. For hyperbolas the terms grow like
 * scale * e^s / 2 with s = sqrt(-alpha) |chi| and scale = max(r0 |a|^(1/2), |sigma0| |a|,
 * |a|^(3/2)), a = 1/alpha; the cap keeps them below ~1e300. A root beyond it would give a
 * state that a double cannot hold. */
static double chi_cap(const universal* k) {
  if (k->alpha >= 0.0) {
    return 0x1p300; /* bounded (ellipse) or polynomial (parabola) terms */
  }
  const double root_a = 1.0 / sqrt(-k->alpha);
  const double a = root_a * root_a;
  double scale = k->r0 * root_a;
  scale = fabs(k->sigma0) * a > scale ? fabs(k->sigma0) * a : scale;
  scale = a * root_a > scale ? a * root_a : scale;
  scale = scale > 1.0 ? scale : 1.0;
  const double s_max = 690.0 - myb_log(scale);
  return s_max > 0.0 ? s_max * root_a : 0.0;
}

/* Bracket [lo, hi] with F(lo) <= 0 <= F(hi), searched by doubling up to chi_cap. */
static myb_status bracket(const universal* k, double* lo, double* hi) {
  const double direction = k->target > 0.0 ? 1.0 : -1.0;
  const double cap = chi_cap(k);
  if (!(cap > 0.0)) {
    return MYB_E_DOMAIN;
  }
  double near = 0.0;
  double far = fabs(k->target) / k->r0;
  far = far < cap ? far : cap;
  /* bound: MAX_DOUBLINGS; |far| doubles each step until the root or the cap is reached */
  for (int i = 0; i < MAX_DOUBLINGS; ++i) {
    const u_values u = u_at(k, direction * far);
    const double f = residual(k, &u) * direction;
    if (f > 0.0) {
      *lo = direction > 0.0 ? near : -far;
      *hi = direction > 0.0 ? far : -near;
      return MYB_OK;
    }
    if (!(f <= 0.0) || far >= cap) {
      return MYB_E_DOMAIN;
    }
    near = far;
    far = 2.0 * far < cap ? 2.0 * far : cap;
  }
  return MYB_E_DOMAIN;
}

static void universal_eval(const void* context, double chi, double* f, double* df) {
  const universal* k = context;
  const u_values u = u_at(k, chi);
  *f = residual(k, &u);
  *df = radius(k, &u);
}

static myb_status solve(const universal* k, double* chi_out) {
  double lo = 0.0;
  double hi = 0.0;
  const myb_status status = bracket(k, &lo, &hi);
  if (status != MYB_OK) {
    return status;
  }
  /* Start: the elliptic guess chi = sqrt(mu) dt alpha (Vallado) when it lies in the bracket. */
  double chi = k->alpha > 0.0 ? k->target * k->alpha : 0.5 * (lo + hi);
  if (!(chi > lo && chi < hi)) {
    chi = 0.5 * (lo + hi);
  }
  const myb_internal_root_function function = {.eval = universal_eval, .context = k};
  return myb_internal_newton_bracketed(&function, lo, hi, chi, MAX_ITERATIONS, chi_out);
}

/* Universal propagation of pv0 (|r0| = r0, alpha = 2/r0 - v0^2/mu) by t, without period
 * reduction. */
static myb_status universal_step(double mu, const myb_pv* pv0, double r0, double alpha, double t,
                                 myb_pv* pv) {
  const double sqrt_mu = sqrt(mu);
  const universal k = {.r0 = r0,
                       .sigma0 = myb_vec3_dot(pv0->r, pv0->v) / sqrt_mu,
                       .alpha = alpha,
                       .target = sqrt_mu * t};
  double chi = 0.0;
  const myb_status status = solve(&k, &chi);
  if (status != MYB_OK) {
    return status;
  }
  const u_values u = u_at(&k, chi);
  const double r = radius(&k, &u);
  const double f = 1.0 - u.u2 / r0;
  const double g = (r0 * u.u1 + k.sigma0 * u.u2) / sqrt_mu;
  const double fdot = -sqrt_mu * u.u1 / (r * r0);
  const double gdot = 1.0 - u.u2 / r;
  myb_pv out;
  myb_vec3_combine(pv0->r, f, pv0->v, g, out.r);
  myb_vec3_combine(pv0->r, fdot, pv0->v, gdot, out.v);
  if (!(r > 0.0) || !myb_vec3_finite(out.r) || !myb_vec3_finite(out.v)) {
    return MYB_E_DOMAIN;
  }
  *pv = out;
  return MYB_OK;
}

/* Hyperbola far from its periapsis (|H0| > 1): periapsis state and time since periapsis.
 * Propagating from there avoids the cancellation of the universal equation when it starts
 * far out (its terms grow like e^s while sqrt(mu) dt stays moderate): at the periapsis
 * sigma0 = 0 and all terms have the sign of chi. Time since periapsis:
 * sinh H0 = r0.v0 / (e sqrt(mu |a|)), M0 = e sinh H0 - H0 (cancellation-free form),
 * tau0 = M0 / n. |H0| > 1 implies r0 / |a| = e cosh H0 - 1 > 0.5, so alpha (from the energy)
 * is well-conditioned there; near-parabolic orbits (|a| huge, H0 small) keep the direct
 * universal route, which does not need alpha accurately. Returns 0 when not applicable. */
static int hyperbola_periapsis(double mu, const myb_pv* pv0, double r0, double alpha,
                               myb_pv* periapsis, double* tau0) {
  double h[3];
  myb_vec3_cross(pv0->r, pv0->v, h);
  const double h_norm = myb_vec3_norm(h);
  if (!(h_norm > 0.0)) {
    return 0;
  }
  double vxh[3];
  double e_vector[3];
  myb_vec3_cross(pv0->v, h, vxh);
  myb_vec3_combine(vxh, 1.0 / mu, pv0->r, -1.0 / r0, e_vector);
  const double e = myb_vec3_norm(e_vector);
  if (!(e > 1.0)) {
    return 0;
  }
  const double rp = h_norm * h_norm / (mu * (1.0 + e));
  double e_hat[3];
  double q_hat[3];
  myb_vec3_combine(e_vector, 1.0 / e, e_vector, 0.0, e_hat);
  myb_vec3_cross(h, e_hat, q_hat);
  myb_vec3_combine(q_hat, 1.0 / h_norm, q_hat, 0.0, q_hat);
  myb_vec3_combine(e_hat, rp, e_hat, 0.0, periapsis->r);
  myb_vec3_combine(q_hat, h_norm / rp, q_hat, 0.0, periapsis->v);
  const double a_abs = -1.0 / alpha;
  const double h0 = myb_asinh(myb_vec3_dot(pv0->r, pv0->v) / (e * sqrt(mu * a_abs)));
  if (!(fabs(h0) > 1.0)) {
    return 0;
  }
  double m0 = 0.0;
  if (myb_anomaly_mean_from_hyperbolic(h0, e, &m0) != MYB_OK) {
    return 0;
  }
  *tau0 = m0 / sqrt(mu / (a_abs * a_abs * a_abs));
  return 1;
}

static const double two_pi_lo = 0x1.1a62633145c07p-52; /* 2 pi = two_pi_hi + two_pi_lo */

/* |a|^2 as a double-word sum of exact products. */
static myb_dd norm2_dd(const double a[3]) {
  myb_dd sum = myb_dd_from_double(0.0);
  for (int i = 0; i < 3; ++i) {
    double err = 0.0;
    const double p = myb_two_prod(a[i], a[i], &err);
    sum = myb_dd_add(sum, (myb_dd){.hi = p, .lo = err});
  }
  return sum;
}

/* alpha = 2/|r0| - |v0|^2/mu in double-word arithmetic. The terms cancel by up to a factor
 * 2a/r0 (39 at the periapsis of e = 0.95); in double that leaves several ulps of error in the
 * period, i.e. a phase error growing with every period removed (measured before: 2.8e-11
 * relative after 30 periods; now at the conditioning of the input). */
static myb_dd alpha_dd(double mu, const myb_pv* pv0) {
  const myb_dd r0 = myb_dd_sqrt(norm2_dd(pv0->r));
  return myb_dd_sub(myb_dd_div(myb_dd_from_double(2.0), r0),
                    myb_dd_div_d(norm2_dd(pv0->v), mu));
}

/* Ellipse: dt minus the nearest whole number of periods P = 2 pi / (sqrt(mu) alpha^(3/2)),
 * in double-word arithmetic. *period receives P rounded to double. */
static double reduce_periods(double mu, myb_dd alpha, double dt, double* turns,
                             double* period) {
  const myb_dd two_pi = {.hi = two_pi_hi, .lo = two_pi_lo};
  const myb_dd denominator = myb_dd_mul(
      myb_dd_mul(myb_dd_sqrt(myb_dd_from_double(mu)), alpha), myb_dd_sqrt(alpha));
  const myb_dd p = myb_dd_div(two_pi, denominator);
  *period = p.hi;
  *turns = __builtin_roundeven(dt / p.hi);
  if (*turns == 0.0) {
    return dt;
  }
  return myb_dd_to_double(myb_dd_sub(myb_dd_from_double(dt), myb_dd_mul_d(p, *turns)));
}

myb_status myb_kepler_propagate(double mu, const myb_pv* pv0, double dt, myb_pv* pv) {
  MYB_CHECK(pv0 != NULL && pv != NULL);
  if (!(mu > 0.0) || !finite(mu) || !finite(dt) || !myb_vec3_finite(pv0->r) ||
      !myb_vec3_finite(pv0->v)) {
    return MYB_E_DOMAIN;
  }
  const double r0 = myb_vec3_norm(pv0->r);
  if (!(r0 > 0.0)) {
    return MYB_E_DOMAIN;
  }
  if (dt == 0.0) {
    *pv = *pv0;
    return MYB_OK;
  }
  const myb_dd alpha_w = alpha_dd(mu, pv0);
  const double alpha = alpha_w.hi;
  if (alpha < 0.0) {
    myb_pv periapsis;
    double tau0 = 0.0;
    if (hyperbola_periapsis(mu, pv0, r0, alpha, &periapsis, &tau0)) {
      const double rp = myb_vec3_norm(periapsis.r);
      const double alpha_p = 2.0 / rp - myb_vec3_dot(periapsis.v, periapsis.v) / mu;
      const double tau = tau0 + dt;
      if (tau == 0.0) {
        *pv = periapsis;
        return MYB_OK;
      }
      return universal_step(mu, &periapsis, rp, alpha_p, tau, pv);
    }
    return universal_step(mu, pv0, r0, alpha, dt, pv);
  }
  double t = dt;
  if (alpha > 0.0) {
    double turns = 0.0;
    double period = 0.0;
    t = reduce_periods(mu, alpha_w, dt, &turns, &period);
    if (t == 0.0) {
      *pv = *pv0;
      return MYB_OK;
    }
  }
  return universal_step(mu, pv0, r0, alpha, t, pv);
}

/* ---- State transition matrix ----------------------------------------------------------------
 * The state after t is r = F r0 + G v0, v = Ft r0 + Gt v0 with
 *   F = 1 - U2/r0, G = (r0 U1 + sigma0 U2)/sqrt(mu), Ft = -sqrt(mu) U1/(r r0), Gt = 1 - U2/r,
 * functions of the scalars p = (r0, sigma0, alpha) of the initial state, directly and through
 * the root chi of the universal Kepler equation K = r0 U1 + sigma0 U2 + U3 - sqrt(mu) t = 0
 * (dK/dchi = r, dK/dr0 = U1, dK/dsigma0 = U2, dK/dalpha = r0 A1 + sigma0 A2 + A3, with
 * A_n = dU_n/dalpha = -(chi U_{n+1} - n U_{n+2}) / 2 at fixed chi). Hence
 *   Phi = [F I, G I; Ft I, Gt I] + [r0; 0] grad F + ... ,
 * every gradient being a combination of grad r0 = (r0/|r0|, 0), grad sigma0 = (v0, r0)/sqrt(mu)
 * and grad alpha = (-2 r0/|r0|^3, -2 v0/mu). */

typedef struct u_values5 {
  double u[6]; /* U0..U5 */
} u_values5;

static u_values5 u5_at(double alpha, double chi) {
  const double z = alpha * chi * chi;
  const double chi2 = chi * chi;
  u_values5 out;
  out.u[0] = myb_internal_stumpff_c0(z);
  out.u[1] = chi * myb_internal_stumpff_c1(z);
  out.u[2] = chi2 * myb_internal_stumpff_c2(z);
  out.u[3] = chi2 * chi * myb_internal_stumpff_c3(z);
  out.u[4] = chi2 * chi2 * myb_internal_stumpff_c4(z);
  out.u[5] = chi2 * chi2 * chi * myb_internal_stumpff_c5(z);
  return out;
}

/* d/dp of a quantity for p = r0, sigma0, alpha. */
typedef struct partials {
  double r0;
  double sigma0;
  double alpha;
} partials;

static void identity6(double stm[36]) {
  for (int i = 0; i < 36; ++i) {
    stm[i] = (i % 7 == 0) ? 1.0 : 0.0;
  }
}

/* Row vector over the 6 initial components: d/dp r0 * grad r0 + ... */
static void gradient(const partials* d, const double grad_r0[6], const double grad_s[6],
                     const double grad_a[6], double out[6]) {
  for (int j = 0; j < 6; ++j) {
    out[j] = d->r0 * grad_r0[j] + d->sigma0 * grad_s[j] + d->alpha * grad_a[j];
  }
}

static myb_status universal_stm(double mu, const myb_pv* pv0, double r0, double alpha, double t,
                                myb_pv* pv, double stm[36]) {
  if (t == 0.0) {
    *pv = *pv0;
    identity6(stm);
    return MYB_OK;
  }
  const double sqrt_mu = sqrt(mu);
  const universal k = {.r0 = r0,
                       .sigma0 = myb_vec3_dot(pv0->r, pv0->v) / sqrt_mu,
                       .alpha = alpha,
                       .target = sqrt_mu * t};
  double chi = 0.0;
  const myb_status status = solve(&k, &chi);
  if (status != MYB_OK) {
    return status;
  }
  const u_values5 uv = u5_at(alpha, chi);
  const double* u = uv.u;
  const double s0 = k.sigma0;
  const double r = r0 * u[0] + s0 * u[1] + u[2];
  /* A_n = dU_n/dalpha at fixed chi */
  const double a0 = -0.5 * chi * u[1];
  const double a1 = -0.5 * (chi * u[2] - u[3]);
  const double a2 = -0.5 * (chi * u[3] - 2.0 * u[4]);
  const double a3 = -0.5 * (chi * u[4] - 3.0 * u[5]);
  /* chi_p = -K_p / r */
  const partials chi_p = {.r0 = -u[1] / r,
                          .sigma0 = -u[2] / r,
                          .alpha = -(r0 * a1 + s0 * a2 + a3) / r};
  /* r = r0 U0 + sigma0 U1 + U2; dr/dchi = -alpha r0 U1 + sigma0 U0 + U1 */
  const double r_chi = -alpha * r0 * u[1] + s0 * u[0] + u[1];
  const partials r_p = {.r0 = u[0] + r_chi * chi_p.r0,
                        .sigma0 = u[1] + r_chi * chi_p.sigma0,
                        .alpha = r0 * a0 + s0 * a1 + a2 + r_chi * chi_p.alpha};
  const double f = 1.0 - u[2] / r0;
  const double g = (r0 * u[1] + s0 * u[2]) / sqrt_mu;
  const double ft = -sqrt_mu * u[1] / (r * r0);
  const double gt = 1.0 - u[2] / r;
  /* dU1 = U0 dchi + A1 dalpha, dU2 = U1 dchi + A2 dalpha */
  const partials u1_p = {.r0 = u[0] * chi_p.r0,
                         .sigma0 = u[0] * chi_p.sigma0,
                         .alpha = u[0] * chi_p.alpha + a1};
  const partials u2_p = {.r0 = u[1] * chi_p.r0,
                         .sigma0 = u[1] * chi_p.sigma0,
                         .alpha = u[1] * chi_p.alpha + a2};
  const partials f_p = {.r0 = -u2_p.r0 / r0 + u[2] / (r0 * r0),
                        .sigma0 = -u2_p.sigma0 / r0,
                        .alpha = -u2_p.alpha / r0};
  const partials g_p = {.r0 = (u[1] + r0 * u1_p.r0 + s0 * u2_p.r0) / sqrt_mu,
                        .sigma0 = (r0 * u1_p.sigma0 + u[2] + s0 * u2_p.sigma0) / sqrt_mu,
                        .alpha = (r0 * u1_p.alpha + s0 * u2_p.alpha) / sqrt_mu};
  /* Ft = -sqrt(mu) U1 / (r r0): dFt = Ft (dU1/U1 - dr/r - dr0/r0), written without 1/U1 */
  const double c = -sqrt_mu / (r * r0);
  const partials ft_p = {.r0 = c * u1_p.r0 - ft * (r_p.r0 / r + 1.0 / r0),
                         .sigma0 = c * u1_p.sigma0 - ft * r_p.sigma0 / r,
                         .alpha = c * u1_p.alpha - ft * r_p.alpha / r};
  /* Gt = 1 - U2 / r */
  const partials gt_p = {.r0 = -u2_p.r0 / r + u[2] * r_p.r0 / (r * r),
                         .sigma0 = -u2_p.sigma0 / r + u[2] * r_p.sigma0 / (r * r),
                         .alpha = -u2_p.alpha / r + u[2] * r_p.alpha / (r * r)};
  double grad_r0[6];
  double grad_s[6];
  double grad_a[6];
  for (int j = 0; j < 3; ++j) {
    grad_r0[j] = pv0->r[j] / r0;
    grad_r0[3 + j] = 0.0;
    grad_s[j] = pv0->v[j] / sqrt_mu;
    grad_s[3 + j] = pv0->r[j] / sqrt_mu;
    grad_a[j] = -2.0 * pv0->r[j] / (r0 * r0 * r0);
    grad_a[3 + j] = -2.0 * pv0->v[j] / mu;
  }
  double df[6];
  double dg[6];
  double dft[6];
  double dgt[6];
  gradient(&f_p, grad_r0, grad_s, grad_a, df);
  gradient(&g_p, grad_r0, grad_s, grad_a, dg);
  gradient(&ft_p, grad_r0, grad_s, grad_a, dft);
  gradient(&gt_p, grad_r0, grad_s, grad_a, dgt);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      const double diag_r = (j == i) ? f : (j == 3 + i) ? g : 0.0;
      const double diag_v = (j == i) ? ft : (j == 3 + i) ? gt : 0.0;
      stm[6 * i + j] = diag_r + pv0->r[i] * df[j] + pv0->v[i] * dg[j];
      stm[6 * (3 + i) + j] = diag_v + pv0->r[i] * dft[j] + pv0->v[i] * dgt[j];
    }
  }
  myb_pv out;
  myb_vec3_combine(pv0->r, f, pv0->v, g, out.r);
  myb_vec3_combine(pv0->r, ft, pv0->v, gt, out.v);
  if (!(r > 0.0) || !myb_vec3_finite(out.r) || !myb_vec3_finite(out.v)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 36; ++i) {
    if (!finite(stm[i])) {
      return MYB_E_DOMAIN;
    }
  }
  *pv = out;
  return MYB_OK;
}

/* a = b c (6 x 6, row-major); a may not alias b or c. */
static void multiply6(const double b[36], const double c[36], double a[36]) {
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double sum = 0.0;
      for (int k = 0; k < 6; ++k) {
        sum += b[6 * i + k] * c[6 * k + j];
      }
      a[6 * i + j] = sum;
    }
  }
}

/* Inverse of a symplectic matrix [A B; C D]: [D^T -B^T; -C^T A^T] (exact, no division). */
static void symplectic_inverse(const double m[36], double out[36]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      out[6 * i + j] = m[6 * (3 + j) + (3 + i)];            /*  D^T */
      out[6 * i + (3 + j)] = -m[6 * j + (3 + i)];           /* -B^T */
      out[6 * (3 + i) + j] = -m[6 * (3 + j) + i];           /* -C^T */
      out[6 * (3 + i) + (3 + j)] = m[6 * j + i];            /*  A^T */
    }
  }
}

myb_status myb_kepler_stm(double mu, const myb_pv* pv0, double dt, myb_pv* pv, double stm[36]) {
  MYB_CHECK(pv0 != NULL && pv != NULL && stm != NULL);
  if (!(mu > 0.0) || !finite(mu) || !finite(dt) || !myb_vec3_finite(pv0->r) ||
      !myb_vec3_finite(pv0->v)) {
    return MYB_E_DOMAIN;
  }
  const double r0 = myb_vec3_norm(pv0->r);
  if (!(r0 > 0.0)) {
    return MYB_E_DOMAIN;
  }
  if (dt == 0.0) {
    *pv = *pv0;
    identity6(stm);
    return MYB_OK;
  }
  const myb_dd alpha_w = alpha_dd(mu, pv0);
  const double alpha = alpha_w.hi;
  if (alpha < 0.0) {
    myb_pv periapsis;
    double tau0 = 0.0;
    if (hyperbola_periapsis(mu, pv0, r0, alpha, &periapsis, &tau0)) {
      /* Phi(t0 -> t) = Phi(p -> t) Phi(p -> t0)^-1: both legs start at the periapsis. */
      const double rp = myb_vec3_norm(periapsis.r);
      const double alpha_p = 2.0 / rp - myb_vec3_dot(periapsis.v, periapsis.v) / mu;
      double to_end[36];
      double to_start[36];
      double back[36];
      myb_pv start;
      myb_status status = universal_stm(mu, &periapsis, rp, alpha_p, tau0 + dt, pv, to_end);
      if (status != MYB_OK) {
        return status;
      }
      status = universal_stm(mu, &periapsis, rp, alpha_p, tau0, &start, to_start);
      if (status != MYB_OK) {
        return status;
      }
      symplectic_inverse(to_start, back);
      multiply6(to_end, back, stm);
      return MYB_OK;
    }
    return universal_stm(mu, pv0, r0, alpha, dt, pv, stm);
  }
  if (alpha > 0.0) {
    double turns = 0.0;
    double period = 0.0;
    const double t = reduce_periods(mu, alpha_w, dt, &turns, &period);
    if (turns != 0.0) {
      const myb_status status = universal_stm(mu, pv0, r0, alpha, t, pv, stm);
      if (status != MYB_OK) {
        return status;
      }
      /* The reduced time t = dt - turns P(alpha) depends on the initial state through the
       * period: dt/dx0 = turns 3 P / (2 alpha) grad alpha, carried by the end velocity and
       * acceleration. */
      const double r = myb_vec3_norm(pv->r);
      const double w = turns * 1.5 * period / alpha;
      for (int j = 0; j < 6; ++j) {
        const double grad_a = j < 3 ? -2.0 * pv0->r[j] / (r0 * r0 * r0) : -2.0 * pv0->v[j - 3] / mu;
        for (int i = 0; i < 3; ++i) {
          stm[6 * i + j] += pv->v[i] * w * grad_a;
          stm[6 * (3 + i) + j] += -mu * pv->r[i] / (r * r * r) * w * grad_a;
        }
      }
      return MYB_OK;
    }
  }
  return universal_stm(mu, pv0, r0, alpha, dt, pv, stm);
}
