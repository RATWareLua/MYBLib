/* SPDX-License-Identifier: BSD-3-Clause */
/* B-plane, patched conics, shooting. */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "myb/cmath.h"
#include "myb/num.h"
#include "myb/transfer.h"
#include "num/vec3.h"

myb_status myb_bplane_compute(double gm, const double r[3], const double v[3], const double pole[3],
                              myb_bplane* bplane) {
  MYB_CHECK(r != NULL && v != NULL && pole != NULL && bplane != NULL);
  const double rr = myb_vec3_norm(r);
  const double vv = myb_vec3_dot(v, v);
  if (!(gm > 0.0) || !(rr > 0.0) || !isfinite(rr) || !isfinite(vv)) {
    return MYB_E_DOMAIN;
  }
  const double vinf2 = vv - 2.0 * gm / rr;
  if (!(vinf2 > 0.0)) {
    return MYB_E_DOMAIN; /* not hyperbolic */
  }
  double h[3];
  myb_vec3_cross(r, v, h);
  const double hh = myb_vec3_norm(h);
  const double rv = myb_vec3_dot(r, v);
  double e[3];
  for (int i = 0; i < 3; ++i) {
    e[i] = ((vv - gm / rr) * r[i] - rv * v[i]) / gm;
  }
  const double ee = myb_vec3_norm(e);
  if (!(hh > 0.0) || !(ee > 1.0)) {
    return MYB_E_DOMAIN;
  }
  myb_bplane* b = bplane;
  memset(b, 0, sizeof *b);
  b->v_infinity = sqrt(vinf2);
  b->c3 = vinf2;
  b->eccentricity = ee;
  double eh[3];
  double hu[3];
  for (int i = 0; i < 3; ++i) {
    eh[i] = e[i] / ee;
    hu[i] = h[i] / hh;
  }
  double he[3];
  myb_vec3_cross(hu, eh, he);
  const double ca = 1.0 / ee;
  const double sa = sqrt((ee - 1.0) * (ee + 1.0)) / ee;
  for (int i = 0; i < 3; ++i) {
    b->s[i] = ca * eh[i] + sa * he[i];
  }
  double t[3];
  myb_vec3_cross(b->s, pole, t);
  const double tn = myb_vec3_norm(t);
  if (!(tn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    b->t[i] = t[i] / tn;
  }
  myb_vec3_cross(b->s, b->t, b->r);
  /* |B| = h / v_inf, direction S x h */
  b->b_magnitude = hh / b->v_infinity;
  double sxh[3];
  myb_vec3_cross(b->s, hu, sxh);
  for (int i = 0; i < 3; ++i) {
    b->b[i] = b->b_magnitude * sxh[i];
  }
  b->b_t = myb_vec3_dot(b->b, b->t);
  b->b_r = myb_vec3_dot(b->b, b->r);
  b->theta = myb_atan2(b->b_r, b->b_t);
  const double a = -gm / vinf2;
  b->periapsis = -a * (ee - 1.0);
  /* time from periapsis by the hyperbolic anomaly: cosh F = (1 - r / a) / e */
  const double ch = (1.0 - rr / a) / ee;
  const double fh = myb_acosh(ch > 1.0 ? ch : 1.0);
  const double mean = ee * myb_sinh(fh) - fh;
  const double n = sqrt(gm / (-a * a * a));
  const double since = (rv >= 0.0 ? 1.0 : -1.0) * mean / n;
  b->time_to_periapsis = -since;
  return MYB_OK;
}

double myb_sphere_of_influence(double a, double mass_ratio) {
  return a * myb_exp(0.4 * myb_log(mass_ratio));
}

double myb_hill_radius(double a, double e, double mass_ratio) {
  return a * (1.0 - e) * myb_cbrt(mass_ratio / 3.0);
}

double myb_flyby_turn(double gm, double v_infinity, double periapsis) {
  return 2.0 * myb_asin(1.0 / (1.0 + periapsis * v_infinity * v_infinity / gm));
}

double myb_departure_impulse(double gm, double radius, double v_infinity) {
  return sqrt(v_infinity * v_infinity + 2.0 * gm / radius) - sqrt(gm / radius);
}

/* ---- Shooting --------------------------------------------------------------------------------- */

enum { MAX_SHOOT = 64 };

static double max_abs(const double* f, size_t m) {
  double r = 0.0;
  for (size_t i = 0; i < m; ++i) {
    r = fmax(r, fabs(f[i]));
  }
  return r;
}

static size_t larger(size_t a, size_t b) {
  return a > b ? a : b;
}

size_t myb_shoot_work_size(size_t n, size_t m) {
  return m * n + 2 * m + 2 * n + larger(larger(n * n, m * m), m * n) + larger(m, n);
}

/* The Newton step du for J du = -f (minimum norm, square or least squares); scratch of
 * max(n^2, m^2, m n) + max(m, n) doubles. */
static myb_status newton_step(size_t n, size_t m, const double* jac, const double* f, double* du,
                              double* scratch) {
  double* a = scratch;
  double* y = scratch + larger(larger(n * n, m * m), m * n);
  if (m == n) {
    memcpy(a, jac, n * n * sizeof(double));
    for (size_t i = 0; i < n; ++i) {
      du[i] = -f[i];
    }
    return myb_lu_solve(n, a, du, 1);
  }
  if (m < n) {
    /* du = -J^T (J J^T)^-1 f */
    for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < m; ++j) {
        double s = 0.0;
        for (size_t k = 0; k < n; ++k) {
          s += jac[i * n + k] * jac[j * n + k];
        }
        a[i * m + j] = s;
      }
    }
    for (size_t i = 0; i < m; ++i) {
      y[i] = -f[i];
    }
    const myb_status st = myb_lu_solve(m, a, y, 1);
    if (st != MYB_OK) {
      return st;
    }
    for (size_t k = 0; k < n; ++k) {
      double s = 0.0;
      for (size_t i = 0; i < m; ++i) {
        s += jac[i * n + k] * y[i];
      }
      du[k] = s;
    }
    return MYB_OK;
  }
  memcpy(a, jac, m * n * sizeof(double));
  for (size_t i = 0; i < m; ++i) {
    y[i] = -f[i];
  }
  const myb_status st = myb_least_squares(m, n, a, y, NULL);
  if (st == MYB_OK) {
    memcpy(du, y, n * sizeof(double));
  }
  return st;
}

myb_status myb_shoot(myb_shoot_function function, void* context, size_t n, size_t m, double* u,
                     const myb_shoot_options* options, double* work, int* iterations,
                     double* residual) {
  MYB_CHECK(function != NULL && u != NULL && options != NULL && work != NULL);
  if (n < 1 || m < 1 || n > MAX_SHOOT || m > MAX_SHOOT || !(options->tolerance > 0.0) ||
      options->max_iterations < 0) {
    return MYB_E_DOMAIN;
  }
  const int max_it = options->max_iterations > 0 ? options->max_iterations : 20;
  double* jac = work;
  double* f = jac + m * n;
  double* ft = f + m;
  double* du = ft + m;
  double* trial = du + n;
  double* scratch = trial + n;
  myb_status st = function(context, u, f, jac);
  if (st != MYB_OK) {
    return st;
  }
  double norm = max_abs(f, m);
  int it = 0;
  for (; it < max_it && !(norm <= options->tolerance); ++it) {
    st = newton_step(n, m, jac, f, du, scratch);
    if (st != MYB_OK) {
      return MYB_E_NOCONV;
    }
    /* backtracking: the full step, then halves while |f| does not decrease */
    double scale = 1.0;
    int accepted = 0;
    for (int half = 0; half <= 10 && !accepted; ++half) {
      for (size_t k = 0; k < n; ++k) {
        trial[k] = u[k] + scale * du[k];
      }
      st = function(context, trial, ft, NULL);
      if (st == MYB_OK && max_abs(ft, m) < norm) {
        accepted = 1;
      } else {
        scale *= 0.5;
      }
    }
    if (!accepted) {
      return MYB_E_NOCONV;
    }
    memcpy(u, trial, n * sizeof(double));
    st = function(context, u, f, jac);
    if (st != MYB_OK) {
      return st;
    }
    norm = max_abs(f, m);
  }
  if (iterations != NULL) {
    *iterations = it;
  }
  if (residual != NULL) {
    *residual = norm;
  }
  return norm <= options->tolerance ? MYB_OK : MYB_E_NOCONV;
}
