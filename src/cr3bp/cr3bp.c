/* SPDX-License-Identifier: BSD-3-Clause */
/* The circular restricted three-body problem: dynamics, Lagrange points, periodic orbits,
 * stability, manifolds, Richardson's halo orbits, ephemeris frames. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "myb/cr3bp.h"
#include "myb/model.h"
#include "myb/num.h"
#include "myb/transfer.h"
#include "myb/uncert.h"

static const double pi = 3.141592653589793;

static double ipow(double x, int n) {
  double r = 1.0;
  for (int i = 0; i < n; ++i) {
    r *= x;
  }
  return r;
}

void myb_cr3bp_acceleration(double mu, const double x[3], const double v[3], double a[3]) {
  MYB_CHECK(x != NULL && v != NULL && a != NULL);
  const double d1[3] = {x[0] + mu, x[1], x[2]};
  const double d2[3] = {x[0] - 1.0 + mu, x[1], x[2]};
  const double r1 = sqrt(d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]);
  const double r2 = sqrt(d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]);
  const double k1 = (1.0 - mu) / (r1 * r1 * r1);
  const double k2 = mu / (r2 * r2 * r2);
  a[0] = 2.0 * v[1] + x[0] - k1 * d1[0] - k2 * d2[0];
  a[1] = -2.0 * v[0] + x[1] - k1 * d1[1] - k2 * d2[1];
  a[2] = -k1 * d1[2] - k2 * d2[2];
}

void myb_cr3bp_jacobian(double mu, const double x[3], double dadx[9], double dadv[9]) {
  MYB_CHECK(x != NULL && dadx != NULL && dadv != NULL);
  const double d1[3] = {x[0] + mu, x[1], x[2]};
  const double d2[3] = {x[0] - 1.0 + mu, x[1], x[2]};
  const double r1s = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2];
  const double r2s = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2];
  const double r1 = sqrt(r1s);
  const double r2 = sqrt(r2s);
  const double m1 = (1.0 - mu) / (r1s * r1);
  const double m2 = mu / (r2s * r2);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double delta = i == j ? 1.0 : 0.0;
      dadx[3 * i + j] = (i == j && i < 2 ? 1.0 : 0.0) + m1 * (3.0 * d1[i] * d1[j] / r1s - delta) +
                        m2 * (3.0 * d2[i] * d2[j] / r2s - delta);
    }
  }
  memset(dadv, 0, 9 * sizeof(double));
  dadv[1] = 2.0;
  dadv[3] = -2.0;
}

double myb_cr3bp_jacobi(double mu, const double y[6]) {
  MYB_CHECK(y != NULL);
  const double d1[3] = {y[0] + mu, y[1], y[2]};
  const double d2[3] = {y[0] - 1.0 + mu, y[1], y[2]};
  const double r1 = sqrt(d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]);
  const double r2 = sqrt(d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]);
  const double omega = 0.5 * (y[0] * y[0] + y[1] * y[1]) + (1.0 - mu) / r1 + mu / r2;
  return 2.0 * omega - (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]);
}

myb_status myb_cr3bp_term(void* mu, double t, const double x[3], const double v[3], double a[3]) {
  (void)t;
  myb_cr3bp_acceleration(*(const double*)mu, x, v, a);
  return MYB_OK;
}

myb_status myb_cr3bp_term_jacobian(void* mu, double t, const double x[3], const double v[3],
                                   double dadx[9], double dadv[9]) {
  (void)t;
  (void)v;
  myb_cr3bp_jacobian(*(const double*)mu, x, dadx, dadv);
  return MYB_OK;
}

/* ---- Lagrange points --------------------------------------------------------------------------- */

static double collinear_force(double mu, double x) {
  const double d1 = x + mu;
  const double d2 = x - 1.0 + mu;
  return x - (1.0 - mu) * d1 / (fabs(d1) * d1 * d1) - mu * d2 / (fabs(d2) * d2 * d2);
}

/* The root of the increasing function on (lo, hi) by bisection to the resolution of doubles. */
static double collinear_root(double mu, double lo, double hi) {
  /* bound: at most 200 halvings */
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (mid <= lo || mid >= hi) {
      break;
    }
    if (collinear_force(mu, mid) < 0.0) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

myb_status myb_cr3bp_lagrange(double mu, double points[5][3]) {
  MYB_CHECK(points != NULL);
  if (!(mu > 0.0) || !(mu <= 0.5)) {
    return MYB_E_DOMAIN;
  }
  const double gap = 1e-9 * fmin(1.0, myb_cbrt(mu));
  memset(points, 0, 5 * sizeof points[0]);
  points[0][0] = collinear_root(mu, -mu + gap, 1.0 - mu - gap);
  points[1][0] = collinear_root(mu, 1.0 - mu + gap, 3.0);
  points[2][0] = collinear_root(mu, -3.0, -mu - gap);
  points[3][0] = 0.5 - mu;
  points[3][1] = 0.8660254037844386;
  points[4][0] = 0.5 - mu;
  points[4][1] = -0.8660254037844386;
  return MYB_OK;
}

/* ---- Propagation with the state transition matrix ---------------------------------------------- */

typedef struct propagator {
  double mu;
  myb_force_model* model;
  myb_variational* var;
} propagator;

static myb_status propagator_open(propagator* p, double mu) {
  memset(p, 0, sizeof *p);
  p->mu = mu;
  myb_status st = myb_force_model_create(NULL, &p->model);
  if (st == MYB_OK) {
    st = myb_force_model_add(p->model, "cr3bp", myb_cr3bp_term, myb_cr3bp_term_jacobian, &p->mu,
                             NULL);
  }
  const myb_variational_options o = {.epsilon = 1e-12, .h_max = 0.02};
  if (st == MYB_OK) {
    st = myb_variational_create(p->model, &o, NULL, &p->var);
  }
  return st;
}

static void propagator_close(propagator* p) {
  myb_variational_destroy(p->var);
  myb_force_model_destroy(p->model);
}

/* Integrates from y0 to the n-th crossing of y = 0 (t > 0); state and Phi there. */
static myb_status to_crossing(propagator* p, const double y0[6], int n, double* t_cross,
                              double y[6], double phi[36]) {
  myb_status st = myb_variational_start(p->var, 0.0, y0);
  if (st != MYB_OK) {
    return st;
  }
  int count = 0;
  double t_prev = 0.0;
  double y_prev = 0.0;
  int first = 1;
  /* bound: steps of at most 0.05 until t = 50 */
  while (st == MYB_OK) {
    st = myb_variational_step(p->var, 50.0);
    if (st != MYB_OK) {
      break;
    }
    const double t = myb_variational_time(p->var);
    double s[6];
    myb_variational_get(p->var, s, NULL, NULL);
    if (!first && ((y_prev < 0.0) != (s[1] < 0.0)) && s[1] != 0.0) {
      ++count;
      if (count == n) {
        /* bisection on the dense output, then the state there */
        double lo = t_prev;
        double hi = t;
        const double sign_lo = y_prev;
        for (int it = 0; it < 100; ++it) {
          const double mid = 0.5 * (lo + hi);
          if (mid <= lo || mid >= hi) {
            break;
          }
          double m[6];
          st = myb_variational_dense(p->var, mid, m, NULL, NULL);
          if (st != MYB_OK) {
            return st;
          }
          if ((m[1] < 0.0) == (sign_lo < 0.0) && m[1] != 0.0) {
            lo = mid;
          } else {
            hi = mid;
          }
        }
        *t_cross = 0.5 * (lo + hi);
        return myb_variational_dense(p->var, *t_cross, y, phi, NULL);
      }
    }
    first = 0;
    t_prev = t;
    y_prev = s[1];
    if (t >= 50.0) {
      return MYB_E_NOCONV;
    }
  }
  return st;
}

myb_status myb_cr3bp_correct(double mu, const myb_cr3bp_correction* correction, double y0[6],
                             double* period) {
  MYB_CHECK(correction != NULL && y0 != NULL && period != NULL);
  const myb_cr3bp_correction* c = correction;
  int free_var[2];
  int constraint[2];
  if (c->symmetry == MYB_CR3BP_SYMMETRY_XZ) {
    if (y0[1] != 0.0 || y0[3] != 0.0 || y0[5] != 0.0 ||
        !(c->fixed == 0 || c->fixed == 2 || c->fixed == 4)) {
      return MYB_E_DOMAIN;
    }
    const int all[3] = {0, 2, 4};
    int k = 0;
    for (int i = 0; i < 3; ++i) {
      if (all[i] != c->fixed) {
        free_var[k++] = all[i];
      }
    }
    constraint[0] = 3; /* vx */
    constraint[1] = 5; /* vz */
  } else if (c->symmetry == MYB_CR3BP_SYMMETRY_X) {
    if (y0[1] != 0.0 || y0[2] != 0.0 || y0[3] != 0.0 ||
        !(c->fixed == 0 || c->fixed == 4 || c->fixed == 5)) {
      return MYB_E_DOMAIN;
    }
    const int all[3] = {0, 4, 5};
    int k = 0;
    for (int i = 0; i < 3; ++i) {
      if (all[i] != c->fixed) {
        free_var[k++] = all[i];
      }
    }
    constraint[0] = 2; /* z */
    constraint[1] = 3; /* vx */
  } else {
    return MYB_E_DOMAIN;
  }
  const int crossing = c->crossing > 0 ? c->crossing : 1;
  const double tolerance = c->tolerance > 0.0 ? c->tolerance : 1e-12;
  const int max_iterations = c->max_iterations > 0 ? c->max_iterations : 30;
  if (!(mu > 0.0 && mu <= 0.5)) {
    return MYB_E_DOMAIN;
  }
  propagator p;
  myb_status st = propagator_open(&p, mu);
  if (st != MYB_OK) {
    propagator_close(&p);
    return st;
  }
  st = MYB_E_NOCONV;
  for (int it = 0; it < max_iterations; ++it) {
    double t = 0.0;
    double y[6];
    double phi[36];
    const myb_status s = to_crossing(&p, y0, crossing, &t, y, phi);
    if (s != MYB_OK) {
      st = s;
      break;
    }
    const double f[2] = {y[constraint[0]], y[constraint[1]]};
    if (fabs(f[0]) <= tolerance && fabs(f[1]) <= tolerance) {
      *period = 2.0 * t;
      st = MYB_OK;
      break;
    }
    /* derivatives of the state at the crossing */
    double a[3];
    myb_cr3bp_acceleration(mu, y, y + 3, a);
    const double ydot[6] = {y[3], y[4], y[5], a[0], a[1], a[2]};
    double m[4];
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        m[2 * i + j] = phi[constraint[i] * 6 + free_var[j]] -
                       ydot[constraint[i]] / ydot[1] * phi[1 * 6 + free_var[j]];
      }
    }
    double delta[2] = {-f[0], -f[1]};
    const myb_status ls = myb_lu_solve(2, m, delta, 1);
    if (ls != MYB_OK) {
      st = MYB_E_NOCONV;
      break;
    }
    y0[free_var[0]] += delta[0];
    y0[free_var[1]] += delta[1];
  }
  propagator_close(&p);
  return st;
}

myb_status myb_cr3bp_monodromy(double mu, const double y0[6], double period, double monodromy[36]) {
  MYB_CHECK(y0 != NULL && monodromy != NULL);
  if (!(mu > 0.0 && mu <= 0.5) || !(period > 0.0) || !isfinite(period)) {
    return MYB_E_DOMAIN;
  }
  propagator p;
  myb_status st = propagator_open(&p, mu);
  if (st == MYB_OK) {
    st = myb_variational_start(p.var, 0.0, y0);
  }
  if (st == MYB_OK) {
    st = myb_variational_integrate(p.var, period);
  }
  if (st == MYB_OK) {
    myb_variational_get(p.var, NULL, monodromy, NULL);
  }
  propagator_close(&p);
  return st;
}

myb_status myb_cr3bp_stability(const double monodromy[36], double re[6], double im[6],
                               double* index) {
  MYB_CHECK(monodromy != NULL && re != NULL && im != NULL);
  double m[36];
  memcpy(m, monodromy, sizeof m);
  const myb_status st = myb_eigenvalues(6, m, re, im);
  if (st != MYB_OK) {
    return st;
  }
  if (index != NULL) {
    double big = 0.0;
    for (int k = 0; k < 6; ++k) {
      big = fmax(big, myb_hypot(re[k], im[k]));
    }
    *index = 0.5 * (big + 1.0 / big);
  }
  return MYB_OK;
}

myb_status myb_cr3bp_manifold(double mu, const double y0[6], double period, double phase,
                              int unstable, int sign, double epsilon, double y[6]) {
  MYB_CHECK(y0 != NULL && y != NULL);
  if (!(phase >= 0.0 && phase <= 1.0) || !(epsilon > 0.0) || (sign != 1 && sign != -1)) {
    return MYB_E_DOMAIN;
  }
  double mono[36];
  myb_status st = myb_cr3bp_monodromy(mu, y0, period, mono);
  if (st != MYB_OK) {
    return st;
  }
  double re[6];
  double im[6];
  double index = 0.0;
  st = myb_cr3bp_stability(mono, re, im, &index);
  if (st != MYB_OK) {
    return st;
  }
  /* the real eigenvalue of largest (unstable) or smallest (stable) modulus */
  int pick = -1;
  for (int k = 0; k < 6; ++k) {
    if (im[k] != 0.0 || !(fabs(re[k]) > 0.0)) {
      continue;
    }
    if (pick < 0 || (unstable ? fabs(re[k]) > fabs(re[pick]) : fabs(re[k]) < fabs(re[pick]))) {
      pick = k;
    }
  }
  if (pick < 0 || (unstable ? !(fabs(re[pick]) > 1.0 + 1e-9) : !(fabs(re[pick]) < 1.0 - 1e-9))) {
    return MYB_E_DOMAIN;
  }
  double v[6];
  double work[36];
  st = myb_eigenvector_real(6, mono, re[pick], v, work);
  if (st != MYB_OK) {
    return st;
  }
  /* carry the eigenvector to the phase: Phi(t) v */
  propagator p;
  st = propagator_open(&p, mu);
  double state[6];
  double phi[36];
  if (st == MYB_OK) {
    st = myb_variational_start(p.var, 0.0, y0);
  }
  if (st == MYB_OK && phase > 0.0) {
    st = myb_variational_integrate(p.var, phase * period);
  }
  if (st == MYB_OK) {
    myb_variational_get(p.var, state, phi, NULL);
  }
  propagator_close(&p);
  if (st != MYB_OK) {
    return st;
  }
  double w[6];
  for (int i = 0; i < 6; ++i) {
    double s = 0.0;
    for (int j = 0; j < 6; ++j) {
      s += phi[i * 6 + j] * v[j];
    }
    w[i] = s;
  }
  const double norm = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
  if (!(norm > 0.0)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 6; ++i) {
    y[i] = state[i] + sign * epsilon * w[i] / norm;
  }
  return MYB_OK;
}

/* ---- Richardson's third-order halo orbits -------------------------------------------------------- */

myb_status myb_cr3bp_richardson(double mu, int point, double az, int branch, double y0[6],
                                double* period) {
  MYB_CHECK(y0 != NULL && period != NULL);
  if (!(mu > 0.0 && mu <= 0.5) || (point != 1 && point != 2) || !(az > 0.0) ||
      (branch != 1 && branch != -1)) {
    return MYB_E_DOMAIN;
  }
  double lp[5][3];
  myb_status st = myb_cr3bp_lagrange(mu, lp);
  if (st != MYB_OK) {
    return st;
  }
  const double xl = lp[point - 1][0];
  const double g = fabs(xl - (1.0 - mu)); /* distance to the smaller primary */
  /* c_n = (1 / g^3) ((+-1)^n mu + (-1)^n (1 - mu) g^(n+1) / (1 -+ g)^(n+1)) */
  double c[5];
  for (int n = 2; n <= 4; ++n) {
    const double sn = (n % 2 == 0) ? 1.0 : -1.0;
    const double gn = ipow(g, n + 1);
    if (point == 1) {
      c[n] = (mu + sn * (1.0 - mu) * gn / ipow(1.0 - g, n + 1)) / (g * g * g);
    } else {
      c[n] = (sn * mu + sn * (1.0 - mu) * gn / ipow(1.0 + g, n + 1)) / (g * g * g);
    }
  }
  const double c2 = c[2];
  const double c3 = c[3];
  const double c4 = c[4];
  /* in-plane frequency: lambda^4 + (c2 - 2) lambda^2 - (c2 - 1)(1 + 2 c2) = 0 */
  const double bq = c2 - 2.0;
  const double cq = -(c2 - 1.0) * (1.0 + 2.0 * c2);
  const double l2s = 0.5 * (-bq + sqrt(bq * bq - 4.0 * cq));
  const double lam = sqrt(l2s);
  const double k = (lam * lam + 1.0 + 2.0 * c2) / (2.0 * lam);
  const double delta = lam * lam - c2;
  const double d1 = 3.0 * lam * lam / k * (k * (6.0 * lam * lam - 1.0) - 2.0 * lam);
  const double d2 = 8.0 * lam * lam / k * (k * (11.0 * lam * lam - 1.0) - 2.0 * lam);
  const double a21 = 3.0 * c3 * (k * k - 2.0) / (4.0 * (1.0 + 2.0 * c2));
  const double a22 = 3.0 * c3 / (4.0 * (1.0 + 2.0 * c2));
  const double a23 = -3.0 * c3 * lam / (4.0 * k * d1) * (3.0 * k * k * k * lam - 6.0 * k * (k - lam) + 4.0);
  const double a24 = -3.0 * c3 * lam / (4.0 * k * d1) * (2.0 + 3.0 * k * lam);
  const double b21 = -3.0 * c3 * lam / (2.0 * d1) * (3.0 * k * lam - 4.0);
  const double b22 = 3.0 * c3 * lam / d1;
  const double d21 = -c3 / (2.0 * lam * lam);
  const double a31 = -9.0 * lam / (4.0 * d2) * (4.0 * c3 * (k * a23 - b21) + k * c4 * (4.0 + k * k)) +
                     (9.0 * lam * lam + 1.0 - c2) / (2.0 * d2) *
                         (3.0 * c3 * (2.0 * a23 - k * b21) + c4 * (2.0 + 3.0 * k * k));
  const double a32 = -1.0 / d2 *
                     (9.0 * lam / 4.0 * (4.0 * c3 * (k * a24 - b22) + k * c4) +
                      1.5 * (9.0 * lam * lam + 1.0 - c2) * (c3 * (k * b22 + d21 - 2.0 * a24) - c4));
  const double b31 = 3.0 / (8.0 * d2) *
                     (8.0 * lam * (3.0 * c3 * (k * b21 - 2.0 * a23) - c4 * (2.0 + 3.0 * k * k)) +
                      (9.0 * lam * lam + 1.0 + 2.0 * c2) *
                          (4.0 * c3 * (k * a23 - b21) + k * c4 * (4.0 + k * k)));
  const double b32 = 1.0 / d2 *
                     (9.0 * lam * (c3 * (k * b22 + d21 - 2.0 * a24) - c4) +
                      3.0 / 8.0 * (9.0 * lam * lam + 1.0 + 2.0 * c2) *
                          (4.0 * c3 * (k * a24 - b22) + k * c4));
  const double d31 = 3.0 / (64.0 * lam * lam) * (4.0 * c3 * a24 + c4);
  const double d32 = 3.0 / (64.0 * lam * lam) * (4.0 * c3 * (a23 - d21) + c4 * (4.0 + k * k));
  const double sden = 2.0 * lam * (lam * (1.0 + k * k) - 2.0 * k);
  const double s1 = (1.5 * c3 * (2.0 * a21 * (k * k - 2.0) - a23 * (k * k + 2.0) - 2.0 * k * b21) -
                     3.0 / 8.0 * c4 * (3.0 * k * k * k * k - 8.0 * k * k + 8.0)) /
                    sden;
  const double s2 = (1.5 * c3 * (2.0 * a22 * (k * k - 2.0) + a24 * (k * k + 2.0) + 2.0 * k * b22 + 5.0 * d21) +
                     3.0 / 8.0 * c4 * (12.0 - k * k)) /
                    sden;
  const double a1 = -1.5 * c3 * (2.0 * a21 + a23 + 5.0 * d21) - 3.0 / 8.0 * c4 * (12.0 - k * k);
  const double a2 = 1.5 * c3 * (a24 - 2.0 * a22) + 9.0 / 8.0 * c4;
  const double l1 = a1 + 2.0 * lam * lam * s1;
  const double l2 = a2 + 2.0 * lam * lam * s2;
  const double az_g = az / g;
  const double ax2 = (-l2 * az_g * az_g - delta) / l1;
  if (!(ax2 > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double ax = sqrt(ax2);
  const double azz = az_g;
  const double dm = branch == 1 ? 1.0 : -1.0;
  const double omega = 1.0 + s1 * ax * ax + s2 * azz * azz;
  /* tau1 = 0: the xz-plane crossing */
  const double x = a21 * ax * ax + a22 * azz * azz - ax + (a23 * ax * ax - a24 * azz * azz) +
                   (a31 * ax * ax * ax - a32 * ax * azz * azz);
  const double z = dm * azz + dm * d21 * ax * azz * (1.0 - 3.0) + dm * (d32 * azz * ax * ax - d31 * azz * azz * azz);
  const double ydot = lam * omega *
                      (k * ax + 2.0 * (b21 * ax * ax - b22 * azz * azz) +
                       3.0 * (b31 * ax * ax * ax - b32 * ax * azz * azz));
  /* Richardson's frame: origin at the point, unit g, x away from the larger primary for L2 and
   * towards it for L1 */
  const double side = 1.0;
  y0[0] = xl + side * g * x;
  y0[1] = 0.0;
  y0[2] = g * z;
  y0[3] = 0.0;
  y0[4] = side * g * ydot;
  y0[5] = 0.0;
  *period = 2.0 * pi / (lam * omega);
  return MYB_OK;
}

/* ---- Ephemeris frames ------------------------------------------------------------------------------ */

typedef struct pulsating {
  double origin[6]; /* barycentre */
  double e[9];      /* rows: x, y, z axes */
  double length;
  double rate;      /* d length / dt */
  double omega;     /* rotation rate about e_z */
  double n;         /* time unit 1 / n */
} pulsating;

static myb_status pulsating_frame(const double b1[6], const double b2[6], double gm1, double gm2,
                                  pulsating* f) {
  if (!(gm1 > 0.0) || !(gm2 > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double m = gm1 + gm2;
  for (int i = 0; i < 6; ++i) {
    f->origin[i] = (gm1 * b1[i] + gm2 * b2[i]) / m;
  }
  const double r[3] = {b2[0] - b1[0], b2[1] - b1[1], b2[2] - b1[2]};
  const double v[3] = {b2[3] - b1[3], b2[4] - b1[4], b2[5] - b1[5]};
  const double rr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  const double h[3] = {r[1] * v[2] - r[2] * v[1], r[2] * v[0] - r[0] * v[2], r[0] * v[1] - r[1] * v[0]};
  const double hh = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
  if (!(rr > 0.0) || !(hh > 0.0) || !isfinite(rr) || !isfinite(hh)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    f->e[i] = r[i] / rr;
    f->e[6 + i] = h[i] / hh;
  }
  f->e[3] = f->e[7] * f->e[2] - f->e[8] * f->e[1];
  f->e[4] = f->e[8] * f->e[0] - f->e[6] * f->e[2];
  f->e[5] = f->e[6] * f->e[1] - f->e[7] * f->e[0];
  f->length = rr;
  f->rate = (r[0] * v[0] + r[1] * v[1] + r[2] * v[2]) / rr;
  f->omega = hh / (rr * rr);
  f->n = sqrt(m / (rr * rr * rr));
  return MYB_OK;
}

myb_status myb_cr3bp_from_ephemeris(const double body1[6], const double body2[6], double gm1,
                                    double gm2, const double state[6], double rotating[6]) {
  MYB_CHECK(body1 != NULL && body2 != NULL && state != NULL && rotating != NULL);
  pulsating f;
  const myb_status st = pulsating_frame(body1, body2, gm1, gm2, &f);
  if (st != MYB_OK) {
    return st;
  }
  double d[3];
  double dv[3];
  for (int i = 0; i < 3; ++i) {
    d[i] = state[i] - f.origin[i];
    dv[i] = state[3 + i] - f.origin[3 + i];
  }
  double p[3];
  double q[3];
  for (int i = 0; i < 3; ++i) {
    p[i] = f.e[3 * i] * d[0] + f.e[3 * i + 1] * d[1] + f.e[3 * i + 2] * d[2];
    q[i] = f.e[3 * i] * dv[0] + f.e[3 * i + 1] * dv[1] + f.e[3 * i + 2] * dv[2];
  }
  /* velocity in the rotating axes: q - omega e_z x p; then the pulsation and the time unit */
  const double w[3] = {q[0] + f.omega * p[1], q[1] - f.omega * p[0], q[2]};
  for (int i = 0; i < 3; ++i) {
    rotating[i] = p[i] / f.length;
    rotating[3 + i] = (w[i] / f.length - rotating[i] * f.rate / f.length) / f.n;
  }
  return MYB_OK;
}

myb_status myb_cr3bp_to_ephemeris(const double body1[6], const double body2[6], double gm1,
                                  double gm2, const double rotating[6], double state[6]) {
  MYB_CHECK(body1 != NULL && body2 != NULL && state != NULL && rotating != NULL);
  pulsating f;
  const myb_status st = pulsating_frame(body1, body2, gm1, gm2, &f);
  if (st != MYB_OK) {
    return st;
  }
  double p[3];
  double w[3];
  for (int i = 0; i < 3; ++i) {
    p[i] = rotating[i] * f.length;
    w[i] = (rotating[3 + i] * f.n + rotating[i] * f.rate / f.length) * f.length;
  }
  const double q[3] = {w[0] - f.omega * p[1], w[1] + f.omega * p[0], w[2]};
  for (int i = 0; i < 3; ++i) {
    state[i] = f.origin[i] + f.e[i] * p[0] + f.e[3 + i] * p[1] + f.e[6 + i] * p[2];
    state[3 + i] = f.origin[3 + i] + f.e[i] * q[0] + f.e[3 + i] * q[1] + f.e[6 + i] * q[2];
  }
  return MYB_OK;
}

/* ---- Multiple shooting and continuation -------------------------------------------------------- */

enum { MAX_UNKNOWNS = 6 * MYB_CR3BP_MAX_PATCHES + 1 };

typedef struct shooting {
  propagator p;
  size_t patches;
  int dropped;           /* the component of the closing condition left to the integral */
  double reference[6];   /* the first state the phase is counted from */
  double flow[6];        /* the flow there */
  int arclength;         /* the last condition: (u - predicted) . tangent = 0 */
  const double* tangent;
  const double* predicted;
} shooting;

static void flow_at(double mu, const double y[6], double f[6]) {
  f[0] = y[3];
  f[1] = y[4];
  f[2] = y[5];
  myb_cr3bp_acceleration(mu, y, y + 3, f + 3);
}

/* The conditions and their Jacobian for u = (states, period): 6 patches equations, one more
 * with the arclength. */
static myb_status shooting_residual(void* context, const double* u, double* f, double* jacobian) {
  shooting* s = context;
  const size_t patches = s->patches;
  const size_t n = 6 * patches + 1;
  const size_t m = 6 * patches + (s->arclength ? 1u : 0u);
  const double period = u[n - 1];
  if (!(period > 0.0) || !isfinite(period)) {
    return MYB_E_DOMAIN;
  }
  const double part = period / (double)patches;
  if (jacobian != NULL) {
    memset(jacobian, 0, m * n * sizeof *jacobian);
  }
  size_t row = 0;
  for (size_t k = 0; k < patches; ++k) {
    const double* start = u + 6 * k;
    const size_t next = k + 1 < patches ? k + 1 : 0;
    double end[6];
    double phi[36];
    myb_status st = myb_variational_start(s->p.var, 0.0, start);
    if (st == MYB_OK) {
      st = myb_variational_integrate(s->p.var, part);
    }
    if (st != MYB_OK) {
      return st;
    }
    myb_variational_get(s->p.var, end, phi, NULL);
    double rate[6];
    flow_at(s->p.mu, end, rate);
    for (int i = 0; i < 6; ++i) {
      if (next == 0 && i == s->dropped) {
        continue;
      }
      f[row] = end[i] - u[6 * next + (size_t)i];
      if (jacobian != NULL) {
        double* line = jacobian + row * n;
        for (int j = 0; j < 6; ++j) {
          line[6 * k + (size_t)j] += phi[6 * i + j];
        }
        line[6 * next + (size_t)i] -= 1.0;
        line[n - 1] = rate[i] / (double)patches;
      }
      ++row;
    }
  }
  /* the phase */
  f[row] = 0.0;
  for (int i = 0; i < 6; ++i) {
    f[row] += (u[i] - s->reference[i]) * s->flow[i];
    if (jacobian != NULL) {
      jacobian[row * n + (size_t)i] = s->flow[i];
    }
  }
  ++row;
  if (s->arclength) {
    f[row] = 0.0;
    for (size_t j = 0; j < n; ++j) {
      f[row] += (u[j] - s->predicted[j]) * s->tangent[j];
      if (jacobian != NULL) {
        jacobian[row * n + j] = s->tangent[j];
      }
    }
  }
  return MYB_OK;
}

static int orbit_valid(double mu, const myb_cr3bp_orbit* orbit) {
  if (!(mu > 0.0 && mu <= 0.5) || orbit->patches < 1 || orbit->patches > MYB_CR3BP_MAX_PATCHES ||
      !(orbit->period > 0.0) || !isfinite(orbit->period)) {
    return 0;
  }
  for (size_t i = 0; i < 6 * orbit->patches; ++i) {
    if (!isfinite(orbit->states[i])) {
      return 0;
    }
  }
  return 1;
}

static void shooting_setup(shooting* s, double mu, const myb_cr3bp_orbit* orbit) {
  s->patches = orbit->patches;
  memcpy(s->reference, orbit->states, sizeof s->reference);
  flow_at(mu, orbit->states, s->flow);
  /* the largest velocity component at the first cut is left to the Jacobi integral */
  s->dropped = 3;
  for (int i = 4; i < 6; ++i) {
    if (fabs(orbit->states[i]) > fabs(orbit->states[s->dropped])) {
      s->dropped = i;
    }
  }
  s->arclength = 0;
  s->tangent = NULL;
  s->predicted = NULL;
}

/* The tangent of the family at the solution u: J t = 0 with `previous` . t = 1, made unit.
 * `work` holds (n + 1) n doubles. 0 if no row closes the system. */
static int family_tangent(shooting* s, const double* u, const double* previous, double* work,
                          double* tangent) {
  const size_t n = 6 * s->patches + 1;
  const size_t m = n - 1;
  double f[MAX_UNKNOWNS];
  /* rows to close the system with: the previous tangent, or the period, then coordinates */
  const size_t candidates[4] = {n - 1, 0, 2, 4};
  const int tries = previous != NULL ? 1 : 4;
  for (int attempt = 0; attempt < tries; ++attempt) {
    s->arclength = 0;
    if (shooting_residual(s, u, f, work) != MYB_OK) {
      return 0;
    }
    double* last = work + m * n;
    double rhs[MAX_UNKNOWNS];
    memset(rhs, 0, sizeof rhs);
    rhs[n - 1] = 1.0;
    if (previous != NULL) {
      memcpy(last, previous, n * sizeof *last);
    } else {
      memset(last, 0, n * sizeof *last);
      last[candidates[attempt]] = 1.0;
    }
    if (myb_lu_solve(n, work, rhs, 1) != MYB_OK) {
      continue;
    }
    double norm = 0.0;
    for (size_t j = 0; j < n; ++j) {
      norm += rhs[j] * rhs[j];
    }
    norm = sqrt(norm);
    if (!(norm > 0.0) || !isfinite(norm) || norm > 1e12) {
      continue; /* the row chosen is across the tangent */
    }
    for (size_t j = 0; j < n; ++j) {
      tangent[j] = rhs[j] / norm;
    }
    return 1;
  }
  return 0;
}

myb_status myb_cr3bp_orbit_from_state(double mu, const double y0[6], double period,
                                      size_t patches, myb_cr3bp_orbit* orbit) {
  MYB_CHECK(y0 != NULL && orbit != NULL);
  if (!(mu > 0.0 && mu <= 0.5) || patches < 1 || patches > MYB_CR3BP_MAX_PATCHES ||
      !(period > 0.0) || !isfinite(period)) {
    return MYB_E_DOMAIN;
  }
  memset(orbit, 0, sizeof *orbit);
  orbit->patches = patches;
  orbit->period = period;
  memcpy(orbit->states, y0, 6 * sizeof *y0);
  propagator p;
  myb_status st = propagator_open(&p, mu);
  if (st == MYB_OK) {
    st = myb_variational_start(p.var, 0.0, y0);
  }
  for (size_t k = 1; k < patches && st == MYB_OK; ++k) {
    st = myb_variational_integrate(p.var, period * (double)k / (double)patches);
    if (st == MYB_OK) {
      myb_variational_get(p.var, orbit->states + 6 * k, NULL, NULL);
    }
  }
  propagator_close(&p);
  return st;
}

/* Newton on the conditions of s for the unknowns u (n of them, m conditions), then the
 * tangent of the family. */
static myb_status shooting_solve(shooting* s, double* u, const double* previous_tangent,
                                 double tolerance, int max_iterations, double* tangent,
                                 int* has_tangent) {
  const size_t n = 6 * s->patches + 1;
  const size_t m = 6 * s->patches + (s->arclength ? 1u : 0u);
  const myb_allocator a = myb_internal_allocator(NULL);
  const size_t shoot_size = myb_shoot_work_size(n, m);
  const size_t tangent_size = (n + 1) * n;
  void* block = NULL;
  myb_status st = myb_internal_alloc_array(&a, shoot_size + tangent_size, sizeof(double), &block);
  if (st != MYB_OK) {
    return st;
  }
  double* work = block;
  const myb_shoot_options options = {.tolerance = tolerance > 0.0 ? tolerance : 1e-11,
                                     .max_iterations = max_iterations > 0 ? max_iterations : 30};
  st = myb_shoot(shooting_residual, s, n, m, u, &options, work, NULL, NULL);
  if (st == MYB_OK) {
    *has_tangent = family_tangent(s, u, previous_tangent, work + shoot_size, tangent);
  }
  myb_internal_release(&a, block);
  return st;
}

myb_status myb_cr3bp_orbit_correct(double mu, myb_cr3bp_orbit* orbit, double tolerance,
                                   int max_iterations) {
  MYB_CHECK(orbit != NULL);
  if (!orbit_valid(mu, orbit) || !(tolerance >= 0.0) || max_iterations < 0) {
    return MYB_E_DOMAIN;
  }
  shooting s;
  myb_status st = propagator_open(&s.p, mu);
  if (st != MYB_OK) {
    propagator_close(&s.p);
    return st;
  }
  shooting_setup(&s, mu, orbit);
  const size_t n = 6 * orbit->patches + 1;
  double u[MAX_UNKNOWNS];
  memcpy(u, orbit->states, (n - 1) * sizeof *u);
  u[n - 1] = orbit->period;
  double tangent[MAX_UNKNOWNS];
  int has_tangent = 0;
  st = shooting_solve(&s, u, orbit->has_tangent ? orbit->tangent : NULL, tolerance,
                      max_iterations, tangent, &has_tangent);
  propagator_close(&s.p);
  if (st != MYB_OK) {
    return st;
  }
  memcpy(orbit->states, u, (n - 1) * sizeof *u);
  orbit->period = u[n - 1];
  if (has_tangent) {
    /* the sense: as before, or towards a longer period, or a larger first coordinate */
    double along = 0.0;
    if (orbit->has_tangent) {
      for (size_t j = 0; j < n; ++j) {
        along += tangent[j] * orbit->tangent[j];
      }
    } else {
      along = tangent[n - 1] != 0.0 ? tangent[n - 1] : tangent[0];
    }
    for (size_t j = 0; j < n; ++j) {
      orbit->tangent[j] = along < 0.0 ? -tangent[j] : tangent[j];
    }
  }
  orbit->has_tangent = has_tangent;
  return MYB_OK;
}

myb_status myb_cr3bp_orbit_continue(double mu, myb_cr3bp_orbit* orbit, double step,
                                    double tolerance, int max_iterations) {
  MYB_CHECK(orbit != NULL);
  if (!orbit_valid(mu, orbit) || !orbit->has_tangent || !isfinite(step) || step == 0.0 ||
      !(tolerance >= 0.0) || max_iterations < 0) {
    return MYB_E_DOMAIN;
  }
  shooting s;
  myb_status st = propagator_open(&s.p, mu);
  if (st != MYB_OK) {
    propagator_close(&s.p);
    return st;
  }
  shooting_setup(&s, mu, orbit);
  const size_t n = 6 * orbit->patches + 1;
  double predicted[MAX_UNKNOWNS];
  double u[MAX_UNKNOWNS];
  for (size_t j = 0; j < n - 1; ++j) {
    predicted[j] = orbit->states[j] + step * orbit->tangent[j];
  }
  predicted[n - 1] = orbit->period + step * orbit->tangent[n - 1];
  memcpy(u, predicted, n * sizeof *u);
  s.arclength = 1;
  s.tangent = orbit->tangent;
  s.predicted = predicted;
  double tangent[MAX_UNKNOWNS];
  int has_tangent = 0;
  st = shooting_solve(&s, u, orbit->tangent, tolerance, max_iterations, tangent, &has_tangent);
  propagator_close(&s.p);
  if (st != MYB_OK) {
    return st;
  }
  if (!has_tangent) {
    return MYB_E_NOCONV;
  }
  double along = 0.0;
  for (size_t j = 0; j < n; ++j) {
    along += tangent[j] * orbit->tangent[j];
  }
  memcpy(orbit->states, u, (n - 1) * sizeof *u);
  orbit->period = u[n - 1];
  for (size_t j = 0; j < n; ++j) {
    orbit->tangent[j] = along < 0.0 ? -tangent[j] : tangent[j];
  }
  return MYB_OK;
}
