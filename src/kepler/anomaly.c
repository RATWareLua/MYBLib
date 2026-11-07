/* SPDX-License-Identifier: BSD-3-Clause */
/* Kepler's equation (elliptic, hyperbolic, parabolic) and conversions between anomalies.
 *
 * Elliptic residual   f(E) = (1 - e) E + e (E - sin E) - M,   f'(E) = (1 - e) + e (1 - cos E)
 * hyperbolic residual f(H) = (e - 1) H + e (sinh H - H) - M,  f'(H) = (e - 1) + e (cosh H - 1)
 * with E - sin E = E^3 c3(E^2), 1 - cos E = E^2 c2(E^2) (and the hyperbolic counterparts):
 * every term is computed with full relative accuracy, also for e -> 1 and small anomalies,
 * where the textbook form E - e sin E - M loses all digits.
 * Both residuals increase monotonically; Newton steps are kept inside a bracket that is
 * known to contain the root (falling back to bisection), so convergence is guaranteed. */
#include <math.h> /* sqrt, fabs, cbrt-free: exact operations only */

#include "core/base.h"
#include "kepler/stumpff.h"
#include "myb/cmath.h"
#include "myb/kepler.h"
#include "num/root.h"

enum { MAX_ITERATIONS = 100 };

static const double pi_hi = 0x1.921fb54442d18p+1;

static int finite(double x) {
  return x - x == 0.0;
}

/* Real root of t^3 + p t - q = 0 for p >= 0, q >= 0, in the stable hyperbolic form
 * t = 2 sqrt(p/3) sinh(asinh((3q / 2p) sqrt(3/p)) / 3). */
static double cubic_root(double p, double q) {
  if (p == 0.0) {
    return myb_cbrt(q);
  }
  const double a = sqrt(p / 3.0);
  return 2.0 * a * myb_sinh(myb_asinh(1.5 * q / (p * a)) / 3.0);
}

typedef struct solver {
  double e;
  double m;  /* target, >= 0 */
  int conic; /* 1 ellipse, -1 hyperbola */
} solver;

static double residual(const solver* s, double x) {
  const double z = s->conic * x * x;
  if (x < 1.0) {
    const double third = x * x * x * myb_internal_stumpff_c3(z);
    return (s->conic > 0 ? (1.0 - s->e) * x + s->e * third : (s->e - 1.0) * x + s->e * third) -
           s->m;
  }
  return s->conic > 0 ? x - s->e * myb_sin(x) - s->m : s->e * myb_sinh(x) - x - s->m;
}

static double slope(const solver* s, double x) {
  const double z = s->conic * x * x;
  const double second = x * x * myb_internal_stumpff_c2(z);
  return (s->conic > 0 ? 1.0 - s->e : s->e - 1.0) + s->e * second;
}

static void kepler_eval(const void* context, double x, double* f, double* df) {
  const solver* s = context;
  *f = residual(s, x);
  *df = slope(s, x);
}

/* Root of the residual in [lo, hi] starting from x. */
static myb_status solve(const solver* s, double lo, double hi, double x, double* root) {
  const myb_internal_root_function function = {.eval = kepler_eval, .context = s};
  return myb_internal_newton_bracketed(&function, lo, hi, x, MAX_ITERATIONS, root);
}

myb_status myb_kepler_elliptic(double mean_anomaly, double e, double* eccentric_anomaly) {
  MYB_CHECK(eccentric_anomaly != NULL);
  if (!finite(mean_anomaly) || !(e >= 0.0 && e < 1.0)) {
    return MYB_E_DOMAIN;
  }
  double turns = 0.0;
  const double r = myb_internal_reduce_angle(mean_anomaly, &turns);
  const double m = fabs(r);
  double root = m;
  if (e != 0.0 && m != 0.0) {
    const solver s = {.e = e, .m = m, .conic = 1};
    /* E - M = e sin E lies in [0, e] for M in [0, pi]. */
    const double start =
        (e > 0.8 && m < 1.0) ? cubic_root(6.0 * (1.0 - e), 6.0 * m) : m + e * myb_sin(m);
    const myb_status status = solve(&s, m, m + e, start, &root);
    if (status != MYB_OK) {
      return status;
    }
  }
  *eccentric_anomaly = myb_internal_add_turns(r < 0.0 ? -root : root, turns);
  return MYB_OK;
}

myb_status myb_kepler_hyperbolic(double mean_anomaly, double e, double* hyperbolic_anomaly) {
  MYB_CHECK(hyperbolic_anomaly != NULL);
  if (!finite(mean_anomaly) || !(e > 1.0) || !finite(e)) {
    return MYB_E_DOMAIN;
  }
  const double m = fabs(mean_anomaly);
  double root = 0.0;
  if (m != 0.0) {
    const solver s = {.e = e, .m = m, .conic = -1};
    /* e sinh H >= M (lower) and (e - 1) sinh H <= M, H^3 / 6 <= M (upper). */
    const double lo = myb_asinh(m / e);
    const double upper_a = myb_asinh(m / (e - 1.0));
    const double upper_b = myb_cbrt(6.0 * m);
    const double hi = upper_a < upper_b ? upper_a : upper_b;
    const double start =
        (e < 1.2 && m < 1.0) ? cubic_root(6.0 * (e - 1.0), 6.0 * m) : myb_log(2.0 * m / e + 1.8);
    const myb_status status = solve(&s, lo, hi, start, &root);
    if (status != MYB_OK) {
      return status;
    }
  }
  *hyperbolic_anomaly = mean_anomaly < 0.0 ? -root : root;
  return MYB_OK;
}

double myb_kepler_parabolic(double mean_anomaly) {
  const double d = 2.0 * myb_sinh(myb_asinh(1.5 * mean_anomaly) / 3.0);
  /* The closed form loses about (asinh(3M/2)/3) ulp for large |M| through the conditioning
   * of sinh; one Newton step on D + D^3/3 - M (derivative 1 + D^2 > 0) restores ~1 ulp. */
  const double residual = d + d * d * d / 3.0 - mean_anomaly;
  return d - residual / (1.0 + d * d);
}

/* ---- Anomalies ------------------------------------------------------------------------------ */

static double principal(double angle) { /* [-pi, pi] -> (-pi, pi] */
  return angle <= -pi_hi ? pi_hi : angle;
}

myb_status myb_anomaly_true_from_eccentric(double eccentric_anomaly, double e, double* nu) {
  MYB_CHECK(nu != NULL);
  if (!finite(eccentric_anomaly) || !(e >= 0.0 && e < 1.0)) {
    return MYB_E_DOMAIN;
  }
  double turns = 0.0;
  const double r = myb_internal_reduce_angle(eccentric_anomaly, &turns);
  const double half = 0.5 * r;
  *nu = principal(2.0 * myb_atan2(sqrt(1.0 + e) * myb_sin(half), sqrt(1.0 - e) * myb_cos(half)));
  return MYB_OK;
}

myb_status myb_anomaly_eccentric_from_true(double nu, double e, double* eccentric_anomaly) {
  MYB_CHECK(eccentric_anomaly != NULL);
  if (!finite(nu) || !(e >= 0.0 && e < 1.0)) {
    return MYB_E_DOMAIN;
  }
  double turns = 0.0;
  const double r = myb_internal_reduce_angle(nu, &turns);
  const double half = 0.5 * r;
  *eccentric_anomaly =
      principal(2.0 * myb_atan2(sqrt(1.0 - e) * myb_sin(half), sqrt(1.0 + e) * myb_cos(half)));
  return MYB_OK;
}

myb_status myb_anomaly_mean_from_eccentric(double eccentric_anomaly, double e,
                                           double* mean_anomaly) {
  MYB_CHECK(mean_anomaly != NULL);
  if (!finite(eccentric_anomaly) || !(e >= 0.0 && e < 1.0)) {
    return MYB_E_DOMAIN;
  }
  double turns = 0.0;
  const double r = myb_internal_reduce_angle(eccentric_anomaly, &turns);
  const solver s = {.e = e, .m = 0.0, .conic = 1};
  const double m = residual(&s, fabs(r));
  *mean_anomaly = myb_internal_add_turns(r < 0.0 ? -m : m, turns);
  return MYB_OK;
}

myb_status myb_anomaly_true_from_hyperbolic(double hyperbolic_anomaly, double e, double* nu) {
  MYB_CHECK(nu != NULL);
  if (!finite(hyperbolic_anomaly) || !(e > 1.0) || !finite(e)) {
    return MYB_E_DOMAIN;
  }
  *nu = 2.0 * myb_atan2(sqrt(e + 1.0) * myb_tanh(0.5 * hyperbolic_anomaly), sqrt(e - 1.0));
  return MYB_OK;
}

myb_status myb_anomaly_hyperbolic_from_true(double nu, double e, double* hyperbolic_anomaly) {
  MYB_CHECK(hyperbolic_anomaly != NULL);
  if (!finite(nu) || !(e > 1.0) || !finite(e)) {
    return MYB_E_DOMAIN;
  }
  double turns = 0.0;
  const double r = myb_internal_reduce_angle(nu, &turns);
  const double half = 0.5 * r;
  const double t = sqrt(e - 1.0) * myb_sin(half) / (sqrt(e + 1.0) * myb_cos(half));
  if (!(t > -1.0 && t < 1.0)) {
    return MYB_E_DOMAIN; /* on or beyond the asymptotes */
  }
  *hyperbolic_anomaly = 2.0 * myb_atanh(t);
  return MYB_OK;
}

myb_status myb_anomaly_mean_from_hyperbolic(double hyperbolic_anomaly, double e,
                                            double* mean_anomaly) {
  MYB_CHECK(mean_anomaly != NULL);
  if (!finite(hyperbolic_anomaly) || !(e > 1.0) || !finite(e)) {
    return MYB_E_DOMAIN;
  }
  const solver s = {.e = e, .m = 0.0, .conic = -1};
  const double m = residual(&s, fabs(hyperbolic_anomaly));
  *mean_anomaly = hyperbolic_anomaly < 0.0 ? -m : m;
  return MYB_OK;
}
