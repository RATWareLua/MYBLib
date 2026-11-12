/* SPDX-License-Identifier: BSD-3-Clause */
/* Initial orbit determination (myb/iod.h): Gibbs, Herrick-Gibbs; from three sights Gauss,
 * Laplace, the double-r iteration and Gooding, each with the same refinement. */
#include <math.h> /* sqrt, fabs */
#include <string.h>

#include "core/check.h"
#include "myb/cmath.h"
#include "myb/iod.h"
#include "myb/kepler.h"
#include "myb/lambert.h"
#include "myb/num.h"
#include "num/vec3.h"

static const double TWO_PI = 6.283185307179586;

enum {
  REFINEMENTS = 50,  /* Newton steps of the refinement of a candidate */
  SCAN = 400,        /* intervals in which the roots of the polynomial of Gauss are looked for */
  BISECTIONS = 200,  /* more than the bits of a double */
  OUTER = 60,        /* Newton steps of the double-r iteration and of Gooding's method */
  HALVINGS = 30      /* halvings of a step that leaves the lines of sight */
};

static int usable(const double v[3]) {
  return myb_vec3_finite(v) && myb_vec3_dot(v, v) > 0.0;
}

myb_status myb_iod_gibbs(double mu, const double r1[3], const double r2[3], const double r3[3],
                         double v2[3], double* out_of_plane) {
  MYB_CHECK(r1 != NULL && r2 != NULL && r3 != NULL && v2 != NULL);
  if (!(mu > 0.0) || !isfinite(mu) || !usable(r1) || !usable(r2) || !usable(r3)) {
    return MYB_E_DOMAIN;
  }
  const double n1 = myb_vec3_norm(r1);
  const double n2 = myb_vec3_norm(r2);
  const double n3 = myb_vec3_norm(r3);
  double c12[3];
  double c23[3];
  double c31[3];
  myb_vec3_cross(r1, r2, c12);
  myb_vec3_cross(r2, r3, c23);
  myb_vec3_cross(r3, r1, c31);
  double n[3];
  double d[3];
  double s[3];
  for (int i = 0; i < 3; ++i) {
    n[i] = n1 * c23[i] + n2 * c31[i] + n3 * c12[i];
    d[i] = c12[i] + c23[i] + c31[i];
    s[i] = (n2 - n3) * r1[i] + (n3 - n1) * r2[i] + (n1 - n2) * r3[i];
  }
  const double nn = myb_vec3_norm(n);
  const double dn = myb_vec3_norm(d);
  const double plane = myb_vec3_norm(c23);
  if (!(nn > 0.0) || !(dn > 0.0) || !(plane > 0.0) || !isfinite(nn + dn)) {
    return MYB_E_DOMAIN; /* on one line */
  }
  if (out_of_plane != NULL) {
    /* the angle between r1 and the plane of r2, r3: the complement of that to its normal */
    const double sine = myb_vec3_dot(r1, c23) / (n1 * plane);
    *out_of_plane = myb_asin(sine > 1.0 ? 1.0 : (sine < -1.0 ? -1.0 : sine));
  }
  double b[3];
  myb_vec3_cross(d, r2, b);
  const double l = sqrt(mu / (nn * dn));
  for (int i = 0; i < 3; ++i) {
    v2[i] = l / n2 * b[i] + l * s[i];
  }
  return MYB_OK;
}

myb_status myb_iod_herrick_gibbs(double mu, const double r1[3], double t1, const double r2[3],
                                 double t2, const double r3[3], double t3, double v2[3]) {
  MYB_CHECK(r1 != NULL && r2 != NULL && r3 != NULL && v2 != NULL);
  if (!(mu > 0.0) || !isfinite(mu) || !usable(r1) || !usable(r2) || !usable(r3) ||
      !(t1 < t2) || !(t2 < t3) || !isfinite(t1) || !isfinite(t3)) {
    return MYB_E_DOMAIN;
  }
  const double d21 = t2 - t1;
  const double d32 = t3 - t2;
  const double d31 = t3 - t1;
  const double n1 = myb_vec3_norm(r1);
  const double n2 = myb_vec3_norm(r2);
  const double n3 = myb_vec3_norm(r3);
  const double k1 = -d32 * (1.0 / (d21 * d31) + mu / (12.0 * n1 * n1 * n1));
  const double k2 = (d32 - d21) * (1.0 / (d21 * d32) + mu / (12.0 * n2 * n2 * n2));
  const double k3 = d21 * (1.0 / (d32 * d31) + mu / (12.0 * n3 * n3 * n3));
  for (int i = 0; i < 3; ++i) {
    v2[i] = k1 * r1[i] + k2 * r2[i] + k3 * r3[i];
  }
  return MYB_OK;
}

/* ---- Gauss ------------------------------------------------------------------------------------ */

typedef struct gauss {
  double mu;
  double tau1; /* t1 - t2 */
  double tau3; /* t3 - t2 */
  double d0;
  double d[3][3]; /* d[i][j] = observer_i . p_j */
  const double (*observer)[3];
  const double (*sight)[3];
} gauss;

static double polynomial(double a, double b, double c, double x) {
  const double x3 = x * x * x;
  return (x3 * x3 * x * x + a * x3 * x3) + (b * x3 + c);
}

/* The three positions from the ranges. */
static void positions(const gauss* g, const double range[3], double r[3][3]) {
  for (int k = 0; k < 3; ++k) {
    for (int i = 0; i < 3; ++i) {
      r[k][i] = g->observer[k][i] + range[k] * g->sight[k][i];
    }
  }
}

/* The velocity at the middle time from the positions and the Lagrange coefficients of the two
 * others. 0 if the coefficients are degenerate. */
static int velocity(const double r[3][3], double f1, double g1, double f3, double g3,
                    double v[3]) {
  const double den = f1 * g3 - f3 * g1;
  if (!(den != 0.0) || !isfinite(den)) {
    return 0;
  }
  for (int i = 0; i < 3; ++i) {
    v[i] = (-f3 * r[0][i] + f1 * r[2][i]) / den;
  }
  return 1;
}

/* Two unit vectors across the sight s, with it a right-handed triad. */
static void across(const double s[3], double a[3], double b[3]) {
  /* the axis least along the sight */
  int axis = 0;
  for (int i = 1; i < 3; ++i) {
    if (fabs(s[i]) < fabs(s[axis])) {
      axis = i;
    }
  }
  double e[3] = {0.0, 0.0, 0.0};
  e[axis] = 1.0;
  myb_vec3_cross(s, e, a);
  const double n = myb_vec3_norm(a);
  for (int i = 0; i < 3; ++i) {
    a[i] /= n;
  }
  myb_vec3_cross(s, a, b);
}

static int refine(const gauss* g, myb_pv middle, myb_iod_solution* out);

/* One candidate of Gauss: from the distance x at the middle time, through the series of Gauss,
 * a state at the middle time, then refined. 1 if it settled with positive ranges. */
static int candidate(const gauss* g, double a_coefficient, double b_coefficient, double x,
                     myb_iod_solution* out) {
  const double mu = g->mu;
  const double tau1 = g->tau1;
  const double tau3 = g->tau3;
  const double tau = tau3 - tau1;
  const double x3 = x * x * x;
  double range[3];
  range[1] = a_coefficient + mu * b_coefficient / x3;
  range[0] = ((6.0 * (g->d[2][0] * tau1 / tau3 + g->d[1][0] * tau / tau3) * x3 +
               mu * g->d[2][0] * (tau * tau - tau1 * tau1) * tau1 / tau3) /
                  (6.0 * x3 + mu * (tau * tau - tau3 * tau3)) -
              g->d[0][0]) / g->d0;
  range[2] = ((6.0 * (g->d[0][2] * tau3 / tau1 - g->d[1][2] * tau / tau1) * x3 +
               mu * g->d[0][2] * (tau * tau - tau3 * tau3) * tau3 / tau1) /
                  (6.0 * x3 + mu * (tau * tau - tau1 * tau1)) -
              g->d[2][2]) / g->d0;
  double r[3][3];
  positions(g, range, r);
  /* the series of Gauss for the Lagrange coefficients */
  const double f1 = 1.0 - 0.5 * mu * tau1 * tau1 / x3;
  const double f3 = 1.0 - 0.5 * mu * tau3 * tau3 / x3;
  const double g1 = tau1 - mu * tau1 * tau1 * tau1 / (6.0 * x3);
  const double g3 = tau3 - mu * tau3 * tau3 * tau3 / (6.0 * x3);
  myb_pv middle;
  memcpy(middle.r, r[1], sizeof middle.r);
  if (!velocity(r, f1, g1, f3, g3, middle.v) || !myb_vec3_finite(middle.r) ||
      !myb_vec3_finite(middle.v)) {
    return 0;
  }
  return refine(g, middle, out);
}

/* The orbit through the three sights from a state at the middle time: Newton's method on the
 * six conditions (the object lies on each sight: two components across it vanish) for the six
 * numbers (r, v) at the middle time, with the transition matrix of the two-body motion. 1 if
 * it settled with positive ranges. */
static int refine(const gauss* g, myb_pv middle, myb_iod_solution* out) {
  const double mu = g->mu;
  double range[3];
  const double taus[3] = {g->tau1, 0.0, g->tau3};
  double basis[3][2][3];
  for (int k = 0; k < 3; ++k) {
    across(g->sight[k], basis[k][0], basis[k][1]);
  }
  int settled = 0;
  int steps = 0;
  /* bound: REFINEMENTS */
  for (steps = 0; steps < REFINEMENTS && !settled; ++steps) {
    double jacobian[36];
    double residual[6];
    double worst_angle = 0.0;
    for (int k = 0; k < 3; ++k) {
      myb_pv at = middle;
      double stm[36] = {0.0};
      if (k == 1) {
        for (int i = 0; i < 6; ++i) {
          stm[7 * i] = 1.0;
        }
      } else if (myb_kepler_stm(mu, &middle, taus[k], &at, stm) != MYB_OK) {
        return 0;
      }
      double d[3];
      for (int i = 0; i < 3; ++i) {
        d[i] = at.r[i] - g->observer[k][i];
      }
      range[k] = myb_vec3_dot(d, g->sight[k]);
      for (int side = 0; side < 2; ++side) {
        const double* e = basis[k][side];
        const int row = 2 * k + side;
        residual[row] = -myb_vec3_dot(e, d);
        worst_angle = fmax(worst_angle, fabs(residual[row]) / fabs(range[k]));
        for (int j = 0; j < 6; ++j) {
          jacobian[6 * row + j] = e[0] * stm[j] + e[1] * stm[6 + j] + e[2] * stm[12 + j];
        }
      }
    }
    if (!isfinite(worst_angle)) {
      return 0;
    }
    if (myb_lu_solve(6, jacobian, residual, 1) != MYB_OK) {
      return 0;
    }
    double change = 0.0;
    const double size_r = myb_vec3_norm(middle.r);
    const double size_v = myb_vec3_norm(middle.v);
    for (int i = 0; i < 3; ++i) {
      middle.r[i] += residual[i];
      middle.v[i] += residual[3 + i];
      change = fmax(change, fabs(residual[i]) / size_r);
      change = fmax(change, fabs(residual[3 + i]) / size_v);
    }
    if (!isfinite(change)) {
      return 0;
    }
    /* settled: the step is at the rounding of the state, the sights are met */
    settled = change < 1e-13 && worst_angle < 1e-9;
  }
  if (!settled) {
    return 0;
  }
  /* the ranges of the orbit settled on */
  for (int k = 0; k < 3; ++k) {
    myb_pv at = middle;
    if (k != 1 && myb_kepler_propagate(mu, &middle, taus[k], &at) != MYB_OK) {
      return 0;
    }
    double d[3];
    for (int i = 0; i < 3; ++i) {
      d[i] = at.r[i] - g->observer[k][i];
    }
    range[k] = myb_vec3_dot(d, g->sight[k]);
    if (!(range[k] > 0.0)) {
      return 0;
    }
  }
  memcpy(out->r, middle.r, sizeof out->r);
  memcpy(out->v, middle.v, sizeof out->v);
  memcpy(out->range, range, sizeof out->range);
  out->iterations = steps;
  return 1;
}

/* The positive roots of x^8 + a x^6 + b x^3 + c (at most three: three changes of sign):
 * sign changes over a scan in the logarithm of the distance, from a thousandth to a thousand
 * times the scale of the problem, each closed by bisection. */
static size_t positive_roots(double a, double b, double c, double scale,
                             double roots[MYB_IOD_MAX_SOLUTIONS]) {
  size_t found = 0;
  double before = scale * 1e-3;
  double f_before = polynomial(a, b, c, before);
  const double ratio = myb_exp(myb_log(1e6) / (double)SCAN);
  for (int k = 0; k < SCAN && found < MYB_IOD_MAX_SOLUTIONS; ++k) {
    const double after = before * ratio;
    const double f_after = polynomial(a, b, c, after);
    if ((f_before < 0.0) != (f_after < 0.0)) {
      double lo = before;
      double hi = after;
      /* bound: BISECTIONS halvings */
      for (int step = 0; step < BISECTIONS; ++step) {
        const double mid = lo + 0.5 * (hi - lo);
        if (mid == lo || mid == hi) {
          break;
        }
        if ((polynomial(a, b, c, mid) < 0.0) == (f_before < 0.0)) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      roots[found++] = hi;
    }
    before = after;
    f_before = f_after;
  }
  return found;
}

/* Two roots may settle on one orbit: solutions[*count] is kept only if it is a new one. */
static void keep(myb_iod_solution* solutions, size_t* count) {
  int fresh = 1;
  const double* r = solutions[*count].r;
  for (size_t s = 0; s < *count; ++s) {
    double d[3];
    for (int i = 0; i < 3; ++i) {
      d[i] = r[i] - solutions[s].r[i];
    }
    fresh = fresh && myb_vec3_norm(d) > 1e-9 * myb_vec3_norm(r);
  }
  *count += fresh ? 1u : 0u;
}

/* The checks common to the methods from three sights; fills the problem. */
static myb_status sights_problem(double mu, const double t[3], const double observer[3][3],
                                 const double sight[3][3], gauss* g) {
  if (!(mu > 0.0) || !isfinite(mu) || !(t[0] < t[1]) || !(t[1] < t[2]) || !isfinite(t[0]) ||
      !isfinite(t[2])) {
    return MYB_E_DOMAIN;
  }
  for (int k = 0; k < 3; ++k) {
    if (!myb_vec3_finite(observer[k]) || !myb_vec3_finite(sight[k]) ||
        !(fabs(myb_vec3_dot(sight[k], sight[k]) - 1.0) < 1e-9)) {
      return MYB_E_DOMAIN;
    }
  }
  memset(g, 0, sizeof *g);
  g->mu = mu;
  g->tau1 = t[0] - t[1];
  g->tau3 = t[2] - t[1];
  g->observer = observer;
  g->sight = sight;
  return MYB_OK;
}

myb_status myb_iod_gauss(double mu, const double t[3], const double observer[3][3],
                         const double sight[3][3],
                         myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS], size_t* count) {
  MYB_CHECK(t != NULL && observer != NULL && sight != NULL && solutions != NULL && count != NULL);
  *count = 0;
  gauss g;
  const myb_status checked = sights_problem(mu, t, observer, sight, &g);
  if (checked != MYB_OK) {
    return checked;
  }
  double p[3][3];
  myb_vec3_cross(sight[1], sight[2], p[0]);
  myb_vec3_cross(sight[0], sight[2], p[1]);
  myb_vec3_cross(sight[0], sight[1], p[2]);
  g.d0 = myb_vec3_dot(sight[0], p[0]);
  if (!(fabs(g.d0) > 1e-14)) {
    return MYB_E_DOMAIN; /* the three sights lie in a plane */
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      g.d[i][j] = myb_vec3_dot(observer[i], p[j]);
    }
  }
  const double tau = g.tau3 - g.tau1;
  const double a_coefficient =
      (-g.d[0][1] * g.tau3 / tau + g.d[1][1] + g.d[2][1] * g.tau1 / tau) / g.d0;
  const double b_coefficient =
      (g.d[0][1] * (g.tau3 * g.tau3 - tau * tau) * g.tau3 / tau +
       g.d[2][1] * (tau * tau - g.tau1 * g.tau1) * g.tau1 / tau) / (6.0 * g.d0);
  const double e = myb_vec3_dot(observer[1], sight[1]);
  const double r22 = myb_vec3_dot(observer[1], observer[1]);
  const double a = -(a_coefficient * a_coefficient + 2.0 * a_coefficient * e + r22);
  const double b = -2.0 * mu * b_coefficient * (a_coefficient + e);
  const double c = -mu * mu * b_coefficient * b_coefficient;
  if (!isfinite(a + b + c)) {
    return MYB_E_DOMAIN;
  }
  double roots[MYB_IOD_MAX_SOLUTIONS];
  const size_t found = positive_roots(a, b, c, fmax(sqrt(r22), myb_cbrt(mu * tau * tau)), roots);
  for (size_t k = 0; k < found; ++k) {
    if (candidate(&g, a_coefficient, b_coefficient, roots[k], &solutions[*count])) {
      keep(solutions, count);
    }
  }
  return *count > 0 ? MYB_OK : MYB_E_NOCONV;
}

/* ---- Laplace -------------------------------------------------------------------------------- */

/* The first and second derivatives at the middle time of the quadratic through three values
 * at tau1, 0, tau3. */
static void derivatives(double tau1, double tau3, const double y[3][3], double first[3],
                        double second[3]) {
  const double w1[3] = {-tau3 / (tau1 * (tau1 - tau3)), -(tau1 + tau3) / (tau1 * tau3),
                        -tau1 / ((tau3 - tau1) * tau3)};
  const double w2[3] = {2.0 / (tau1 * (tau1 - tau3)), 2.0 / (tau1 * tau3),
                        2.0 / ((tau3 - tau1) * tau3)};
  for (int i = 0; i < 3; ++i) {
    first[i] = w1[0] * y[0][i] + w1[1] * y[1][i] + w1[2] * y[2][i];
    second[i] = w2[0] * y[0][i] + w2[1] * y[1][i] + w2[2] * y[2][i];
  }
}

static double triple(const double a[3], const double b[3], const double c[3]) {
  double cross[3];
  myb_vec3_cross(b, c, cross);
  return myb_vec3_dot(a, cross);
}

myb_status myb_iod_laplace(double mu, const double t[3], const double observer[3][3],
                           const double sight[3][3],
                           myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS], size_t* count) {
  MYB_CHECK(t != NULL && observer != NULL && sight != NULL && solutions != NULL && count != NULL);
  *count = 0;
  gauss g;
  const myb_status checked = sights_problem(mu, t, observer, sight, &g);
  if (checked != MYB_OK) {
    return checked;
  }
  /* the sight and the observer at the middle time with their derivatives */
  const double* l = sight[1];
  const double* place = observer[1];
  double l1[3];
  double l2[3];
  double place1[3];
  double place2[3];
  derivatives(g.tau1, g.tau3, sight, l1, l2);
  derivatives(g.tau1, g.tau3, observer, place1, place2);
  /* r = R + rho L, r'' = -mu r / r^3: across L and L' the range, across L and L'' its rate */
  const double d = 2.0 * triple(l, l1, l2);
  if (!(fabs(d) > 0.0) || !isfinite(d)) {
    return MYB_E_DOMAIN; /* the sights do not turn: no curvature to measure */
  }
  const double a_range = -2.0 * triple(l, l1, place2) / d;
  const double b_range = -2.0 * mu * triple(l, l1, place) / d;
  const double c = myb_vec3_dot(l, place);
  const double r22 = myb_vec3_dot(place, place);
  /* r^2 = rho^2 + 2 rho c + R^2 with rho = a + b / r^3: a polynomial of the eighth degree */
  const double pa = -(a_range * a_range + 2.0 * a_range * c + r22);
  const double pb = -2.0 * b_range * (a_range + c);
  const double pc = -b_range * b_range;
  if (!isfinite(pa + pb + pc)) {
    return MYB_E_DOMAIN;
  }
  const double tau = g.tau3 - g.tau1;
  double roots[MYB_IOD_MAX_SOLUTIONS];
  const size_t found =
      positive_roots(pa, pb, pc, fmax(sqrt(r22), myb_cbrt(mu * tau * tau)), roots);
  for (size_t k = 0; k < found; ++k) {
    const double r3 = roots[k] * roots[k] * roots[k];
    const double rho = a_range + b_range / r3;
    const double rate = 2.0 * (triple(l, l2, place2) + mu / r3 * triple(l, l2, place)) / d;
    if (!(rho > 0.0)) {
      continue;
    }
    myb_pv middle;
    for (int i = 0; i < 3; ++i) {
      middle.r[i] = place[i] + rho * l[i];
      middle.v[i] = place1[i] + rate * l[i] + rho * l1[i];
    }
    if (myb_vec3_finite(middle.v) && refine(&g, middle, &solutions[*count])) {
      keep(solutions, count);
    }
  }
  return *count > 0 ? MYB_OK : MYB_E_NOCONV;
}

/* ---- The double-r iteration --------------------------------------------------------------- */

/* The point of sight k at the distance r from the centre (the far one of the two). 0 if the
 * line does not reach it. */
static int at_distance(const gauss* g, int k, double r, double point[3]) {
  const double c = myb_vec3_dot(g->sight[k], g->observer[k]);
  const double disc = c * c - myb_vec3_dot(g->observer[k], g->observer[k]) + r * r;
  if (!(disc >= 0.0) || !(r > 0.0)) {
    return 0;
  }
  const double rho = -c + sqrt(disc);
  if (!(rho > 0.0)) {
    return 0;
  }
  for (int i = 0; i < 3; ++i) {
    point[i] = g->observer[k][i] + rho * g->sight[k][i];
  }
  return 1;
}

/* The mean anomaly (elliptic or hyperbolic) of the true anomaly nu. */
static int mean_of(double nu, double e, double* mean) {
  double anomaly = 0.0;
  if (e < 1.0) {
    return myb_anomaly_eccentric_from_true(nu, e, &anomaly) == MYB_OK &&
           myb_anomaly_mean_from_eccentric(anomaly, e, mean) == MYB_OK;
  }
  return myb_anomaly_hyperbolic_from_true(nu, e, &anomaly) == MYB_OK &&
         myb_anomaly_mean_from_hyperbolic(anomaly, e, mean) == MYB_OK;
}

/* The times of flight from the second point to the first and the third along the conic
 * through the three points (the orbit of Gibbs) for the distances r1, r2 at the first two
 * times; the middle state. 0 where there is none (a line that misses the distance, the points
 * on one line, a nearly parabolic conic). */
static int double_r_times(const gauss* g, double r1, double r2, double times[2], myb_pv* middle) {
  double p[3][3];
  if (!at_distance(g, 0, r1, p[0]) || !at_distance(g, 1, r2, p[1])) {
    return 0;
  }
  /* the third point: where its sight meets the plane of the first two */
  double w[3];
  myb_vec3_cross(p[0], p[1], w);
  const double across3 = myb_vec3_dot(g->sight[2], w);
  if (!(fabs(across3) > 0.0)) {
    return 0;
  }
  const double rho3 = -myb_vec3_dot(g->observer[2], w) / across3;
  if (!(rho3 > 0.0)) {
    return 0;
  }
  for (int i = 0; i < 3; ++i) {
    p[2][i] = g->observer[2][i] + rho3 * g->sight[2][i];
  }
  memcpy(middle->r, p[1], sizeof middle->r);
  if (myb_iod_gibbs(g->mu, p[0], p[1], p[2], middle->v, NULL) != MYB_OK) {
    return 0;
  }
  myb_elements el;
  if (myb_elements_from_pv(g->mu, middle, &el) != MYB_OK || !(fabs(el.e - 1.0) > 1e-9)) {
    return 0;
  }
  double h[3];
  myb_vec3_cross(middle->r, middle->v, h);
  double m2 = 0.0;
  if (!mean_of(el.nu, el.e, &m2)) {
    return 0;
  }
  const double a = el.p / (1.0 - el.e * el.e);
  const double n = sqrt(g->mu / fabs(a * a * a));
  for (int s = 0; s < 2; ++s) {
    const double* q = p[s == 0 ? 0 : 2];
    double cross[3];
    myb_vec3_cross(middle->r, q, cross);
    /* the angle from the second point, in the sense of the motion: back to the first, on to
     * the third */
    double angle = myb_atan2(myb_vec3_dot(cross, h) / myb_vec3_norm(h), myb_vec3_dot(middle->r, q));
    if (s == 0 && angle > 0.0) {
      angle -= TWO_PI;
    } else if (s == 1 && angle < 0.0) {
      angle += TWO_PI;
    }
    double m = 0.0;
    if (!mean_of(el.nu + angle, el.e, &m)) {
      return 0;
    }
    double dm = m - m2;
    if (el.e < 1.0) {
      /* the mean anomalies are given within a turn: the difference in the sense of the angle */
      dm -= TWO_PI * floor(dm / TWO_PI);
      if (s == 0) {
        dm -= TWO_PI;
      }
    }
    times[s] = dm / n;
  }
  return isfinite(times[0]) && isfinite(times[1]);
}

myb_status myb_iod_double_r(double mu, const double t[3], const double observer[3][3],
                            const double sight[3][3], const double distance[2],
                            myb_iod_solution* solution) {
  MYB_CHECK(t != NULL && observer != NULL && sight != NULL && distance != NULL &&
            solution != NULL);
  gauss g;
  const myb_status checked = sights_problem(mu, t, observer, sight, &g);
  if (checked != MYB_OK) {
    return checked;
  }
  if (!(distance[0] > 0.0) || !(distance[1] > 0.0) || !isfinite(distance[0] + distance[1])) {
    return MYB_E_DOMAIN;
  }
  double r[2] = {distance[0], distance[1]};
  const double span = g.tau3 - g.tau1;
  myb_pv middle;
  double times[2];
  if (!double_r_times(&g, r[0], r[1], times, &middle)) {
    return MYB_E_NOCONV;
  }
  int settled = 0;
  /* Newton's method on the two times of flight, the derivatives by central differences */
  for (int step = 0; step < OUTER && !settled; ++step) { /* bound: OUTER */
    const double f[2] = {times[0] - g.tau1, times[1] - g.tau3};
    if (fabs(f[0]) + fabs(f[1]) < 1e-11 * span) {
      settled = 1;
      break;
    }
    double jacobian[4];
    for (int j = 0; j < 2; ++j) {
      const double h = 1e-6 * r[j];
      double up[2] = {r[0], r[1]};
      double down[2] = {r[0], r[1]};
      up[j] += h;
      down[j] -= h;
      double tu[2];
      double td[2];
      myb_pv ignored;
      if (!double_r_times(&g, up[0], up[1], tu, &ignored) ||
          !double_r_times(&g, down[0], down[1], td, &ignored)) {
        return MYB_E_NOCONV;
      }
      jacobian[j] = (tu[0] - td[0]) / (2.0 * h);
      jacobian[2 + j] = (tu[1] - td[1]) / (2.0 * h);
    }
    const double det = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
    if (!(fabs(det) > 0.0) || !isfinite(det)) {
      return MYB_E_NOCONV;
    }
    const double dr[2] = {-(jacobian[3] * f[0] - jacobian[1] * f[1]) / det,
                          -(-jacobian[2] * f[0] + jacobian[0] * f[1]) / det};
    /* damped: a step is halved until the lines reach the distances and the error of the times
     * falls */
    const double before = fabs(f[0]) + fabs(f[1]);
    double scale = 1.0;
    int moved = 0;
    for (int k = 0; k < HALVINGS && !moved; ++k) { /* bound: HALVINGS */
      const double next[2] = {r[0] + scale * dr[0], r[1] + scale * dr[1]};
      double trial[2];
      myb_pv state;
      if (double_r_times(&g, next[0], next[1], trial, &state) &&
          fabs(trial[0] - g.tau1) + fabs(trial[1] - g.tau3) < (1.0 - 0.25 * scale) * before) {
        r[0] = next[0];
        r[1] = next[1];
        times[0] = trial[0];
        times[1] = trial[1];
        middle = state;
        moved = 1;
      }
      scale *= 0.5;
    }
    if (!moved) {
      return MYB_E_NOCONV;
    }
  }
  if (!settled && !(fabs(times[0] - g.tau1) + fabs(times[1] - g.tau3) < 1e-9 * span)) {
    return MYB_E_NOCONV;
  }
  return refine(&g, middle, solution) ? MYB_OK : MYB_E_NOCONV;
}

/* ---- Gooding ------------------------------------------------------------------------------ */

/* The transfer from the first sight at the range rho1 to the third at rho3 (Lambert's problem
 * with the revolutions and the branch asked, in the sense `retrograde`), the state at the
 * middle time on it and the two components of its miss across the second sight with their
 * derivatives in the two ranges. 0 where there is none. */
static int gooding_miss(const gauss* g, double rho1, double rho3, int revolutions, int branch,
                        int retrograde, const double basis[2][3], double miss[2],
                        double jacobian[4], myb_pv* middle) {
  myb_pv first;
  double third[3];
  for (int i = 0; i < 3; ++i) {
    first.r[i] = g->observer[0][i] + rho1 * g->sight[0][i];
    third[i] = g->observer[2][i] + rho3 * g->sight[2][i];
  }
  myb_lambert_solution found[1 + 2 * 8];
  size_t count = 0;
  if (myb_lambert(g->mu, first.r, third, g->tau3 - g->tau1, retrograde, revolutions, found,
                  1 + 2 * 8, &count) != MYB_OK) {
    return 0;
  }
  int chosen = -1;
  for (size_t k = 0; k < count; ++k) {
    if (found[k].revolutions == revolutions && (revolutions == 0 || found[k].branch == branch)) {
      chosen = (int)k;
    }
  }
  if (chosen < 0) {
    return 0;
  }
  memcpy(first.v, found[chosen].v1, sizeof first.v);
  /* the transfer and the arc to the middle time with their transition matrices */
  myb_pv end;
  double whole[36];
  double part[36];
  if (myb_kepler_stm(g->mu, &first, g->tau3 - g->tau1, &end, whole) != MYB_OK ||
      myb_kepler_stm(g->mu, &first, -g->tau1, middle, part) != MYB_OK) {
    return 0;
  }
  /* v1 as the ends move: d v1 / d r3 = B^-1, d v1 / d r1 = -B^-1 A (A, B the position rows of
   * the whole transfer) */
  double b[9];
  double b_inverse[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      b[3 * i + j] = whole[6 * i + 3 + j];
    }
  }
  if (myb_lu_solve(3, b, b_inverse, 3) != MYB_OK) {
    return 0;
  }
  double dv_dr1[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double sum = 0.0;
      for (int k = 0; k < 3; ++k) {
        sum += b_inverse[3 * i + k] * whole[6 * k + j];
      }
      dv_dr1[3 * i + j] = -sum;
    }
  }
  /* d r2 / d rho1 = (C + D dv1/dr1) L1, d r2 / d rho3 = D B^-1 L3 (C, D of the part) */
  double by1[3];
  double by3[3];
  for (int i = 0; i < 3; ++i) {
    by1[i] = 0.0;
    by3[i] = 0.0;
    for (int j = 0; j < 3; ++j) {
      double through_v1 = 0.0;
      double through_v3 = 0.0;
      for (int k = 0; k < 3; ++k) {
        through_v1 += part[6 * i + 3 + k] * dv_dr1[3 * k + j];
        through_v3 += part[6 * i + 3 + k] * b_inverse[3 * k + j];
      }
      by1[i] += (part[6 * i + j] + through_v1) * g->sight[0][j];
      by3[i] += through_v3 * g->sight[2][j];
    }
  }
  double d[3];
  for (int i = 0; i < 3; ++i) {
    d[i] = middle->r[i] - g->observer[1][i];
  }
  if (!(myb_vec3_dot(d, g->sight[1]) > 0.0)) {
    return 0; /* behind the second observer */
  }
  for (int s = 0; s < 2; ++s) {
    miss[s] = myb_vec3_dot(basis[s], d);
    jacobian[2 * s] = myb_vec3_dot(basis[s], by1);
    jacobian[2 * s + 1] = myb_vec3_dot(basis[s], by3);
  }
  return isfinite(miss[0] + miss[1] + jacobian[0] + jacobian[1] + jacobian[2] + jacobian[3]);
}

myb_status myb_iod_gooding(double mu, const double t[3], const double observer[3][3],
                           const double sight[3][3], const double range[2], int revolutions,
                           int branch, myb_iod_solution* solution) {
  MYB_CHECK(t != NULL && observer != NULL && sight != NULL && range != NULL && solution != NULL);
  gauss g;
  const myb_status checked = sights_problem(mu, t, observer, sight, &g);
  if (checked != MYB_OK) {
    return checked;
  }
  if (!(range[0] > 0.0) || !(range[1] > 0.0) || !isfinite(range[0] + range[1]) ||
      revolutions < 0 || revolutions > 8 || (revolutions > 0 && branch != -1 && branch != 1)) {
    return MYB_E_DOMAIN;
  }
  double basis[2][3];
  across(sight[1], basis[0], basis[1]);
  double rho[2] = {range[0], range[1]};
  double miss[2];
  double jacobian[4];
  myb_pv middle;
  /* the sense of the motion: the one that misses the second sight less at the guess */
  int retrograde = -1;
  double least = (double)INFINITY;
  for (int sense = 0; sense < 2; ++sense) {
    if (gooding_miss(&g, rho[0], rho[1], revolutions, branch, sense, basis, miss, jacobian,
                     &middle) &&
        myb_hypot(miss[0], miss[1]) < least) {
      least = myb_hypot(miss[0], miss[1]);
      retrograde = sense;
    }
  }
  if (retrograde < 0) {
    return MYB_E_NOCONV;
  }
  if (!gooding_miss(&g, rho[0], rho[1], revolutions, branch, retrograde, basis, miss, jacobian,
                    &middle)) {
    return MYB_E_NOCONV;
  }
  int settled = 0;
  for (int step = 0; step < OUTER && !settled; ++step) { /* bound: OUTER */
    const double distance = myb_vec3_norm(middle.r);
    if (myb_hypot(miss[0], miss[1]) < 1e-12 * distance) {
      settled = 1;
      break;
    }
    const double det = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
    if (!(fabs(det) > 0.0) || !isfinite(det)) {
      return MYB_E_NOCONV;
    }
    const double dr[2] = {-(jacobian[3] * miss[0] - jacobian[1] * miss[1]) / det,
                          -(-jacobian[2] * miss[0] + jacobian[0] * miss[1]) / det};
    /* damped: a step is halved until the ranges are positive, the transfer has the solution
     * asked and the miss falls */
    const double before = myb_hypot(miss[0], miss[1]);
    double scale = 1.0;
    int moved = 0;
    for (int k = 0; k < HALVINGS && !moved; ++k) { /* bound: HALVINGS */
      const double next[2] = {rho[0] + scale * dr[0], rho[1] + scale * dr[1]};
      double trial_miss[2];
      double trial_jacobian[4];
      myb_pv state;
      if (next[0] > 0.0 && next[1] > 0.0 &&
          gooding_miss(&g, next[0], next[1], revolutions, branch, retrograde, basis, trial_miss,
                       trial_jacobian, &state) &&
          myb_hypot(trial_miss[0], trial_miss[1]) < (1.0 - 0.25 * scale) * before) {
        rho[0] = next[0];
        rho[1] = next[1];
        memcpy(miss, trial_miss, sizeof miss);
        memcpy(jacobian, trial_jacobian, sizeof jacobian);
        middle = state;
        moved = 1;
      }
      scale *= 0.5;
    }
    if (!moved) {
      return MYB_E_NOCONV;
    }
  }
  if (!settled && !(myb_hypot(miss[0], miss[1]) < 1e-9 * myb_vec3_norm(middle.r))) {
    return MYB_E_NOCONV;
  }
  return refine(&g, middle, solution) ? MYB_OK : MYB_E_NOCONV;
}
