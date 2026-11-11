/* SPDX-License-Identifier: BSD-3-Clause */
/* Initial orbit determination (myb/iod.h): Gibbs, Herrick-Gibbs, Gauss with refinement. */
#include <math.h> /* sqrt, fabs */
#include <string.h>

#include "core/check.h"
#include "myb/cmath.h"
#include "myb/iod.h"
#include "myb/kepler.h"
#include "myb/num.h"
#include "num/vec3.h"

enum {
  REFINEMENTS = 50,  /* Newton steps of the refinement of a candidate */
  SCAN = 400,        /* intervals in which the roots of the polynomial of Gauss are looked for */
  BISECTIONS = 200   /* more than the bits of a double */
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

/* One candidate: from the distance x at the middle time, through the series of Gauss, to the
 * orbit that passes through the three sights — Newton's method on the six conditions (the
 * object lies on each sight: two components across it vanish) for the six numbers (r, v) at
 * the middle time, with the transition matrix of the two-body motion. 1 if it settled with
 * positive ranges. */
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
  const double taus[3] = {tau1, 0.0, tau3};
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

myb_status myb_iod_gauss(double mu, const double t[3], const double observer[3][3],
                         const double sight[3][3],
                         myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS], size_t* count) {
  MYB_CHECK(t != NULL && observer != NULL && sight != NULL && solutions != NULL && count != NULL);
  *count = 0;
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
  gauss g = {.mu = mu, .tau1 = t[0] - t[1], .tau3 = t[2] - t[1], .observer = observer,
             .sight = sight};
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
  /* the positive roots: sign changes over a scan in the logarithm of the distance, from a
   * thousandth to a thousand times the scale of the problem */
  const double scale = fmax(sqrt(r22), myb_cbrt(mu * tau * tau));
  double before = scale * 1e-3;
  double f_before = polynomial(a, b, c, before);
  const double ratio = myb_exp(myb_log(1e6) / (double)SCAN);
  for (int k = 0; k < SCAN && *count < MYB_IOD_MAX_SOLUTIONS; ++k) {
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
      if (candidate(&g, a_coefficient, b_coefficient, hi, &solutions[*count])) {
        /* two roots may settle on one orbit: it is given once */
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
    }
    before = after;
    f_before = f_after;
  }
  return *count > 0 ? MYB_OK : MYB_E_NOCONV;
}
