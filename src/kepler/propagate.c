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
  const double alpha = 2.0 / r0 - myb_vec3_dot(pv0->v, pv0->v) / mu;
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
    /* Whole periods are removed exactly in double-word arithmetic. */
    const double period = two_pi_hi / (sqrt(mu) * alpha * sqrt(alpha));
    const double turns = __builtin_roundeven(dt / period);
    if (turns != 0.0) {
      double err;
      const double p = myb_two_prod(turns, period, &err);
      t = myb_dd_to_double(myb_dd_add_d(myb_dd_add_d(myb_dd_from_double(dt), -p), -err));
    }
    if (t == 0.0) {
      *pv = *pv0;
      return MYB_OK;
    }
  }
  return universal_step(mu, pv0, r0, alpha, t, pv);
}
