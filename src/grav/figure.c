/* SPDX-License-Identifier: BSD-3-Clause */
/* Mutual gravitation of two extended bodies beyond "field of one at the centre of the other":
 * the quadrupole-quadrupole term of the mutual potential.
 *
 * With the trace-free quadrupoles Q = integral rho (3 x x^T - |x|^2 1) dV (= tr(I) 1 - 3 I for
 * the inertia tensor I) and the STF moments M^ij = Q^ij / 3, the mutual potential of two bodies
 * (positive, U = G M_A M_B / r + ...; r = x_B - x_A, n = r / r) is
 *   U = G sum_{l,l'} (-1)^l' / (l! l'!) M_A^L M_B^L' d_LL'(1/r)
 * (T. Damour, M. Soffel, C. Xu, Phys. Rev. D 43, 3273, 1991; W. Hou, D. J. Scheeres, X. Xin,
 * Celest. Mech. Dyn. Astron. 127, 369, 2017). The terms with l' = 0 or l = 0 are the field of one
 * body at the centre of mass of the other, which myb evaluates with its spherical harmonics to
 * any degree; the first term coupling the two figures is l = l' = 2:
 *   U22 = G / (36 r^5) [105 (n.QA.n)(n.QB.n) - 60 n.QA QB.n + 6 tr(QA QB)],
 * from d_ijkl(1/r) contracted with two trace-free tensors. The force on B is grad_r U22. */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "myb/grav.h"
#include "num/vec3.h"

void myb_quadrupole_from_coefficients(double gm, double radius, const double c[3],
                                      const double s[2], double gq[9]) {
  MYB_CHECK(c != NULL && s != NULL && gq != NULL);
  /* unnormalised coefficients (geodesy normalisation, no Condon-Shortley phase) */
  const double c20 = 2.23606797749979 * c[0];   /* sqrt(5) */
  const double c21 = 1.2909944487358056 * c[1]; /* sqrt(5/3) */
  const double c22 = 0.6454972243679028 * c[2]; /* sqrt(5/12) */
  const double s21 = 1.2909944487358056 * s[0];
  const double s22 = 0.6454972243679028 * s[1];
  /* V2 = GM R^2 / r^3 n.A.n = G n.Q.n / (2 r^3) */
  const double k = 2.0 * gm * radius * radius;
  gq[0] = k * (-0.5 * c20 + 3.0 * c22);
  gq[4] = k * (-0.5 * c20 - 3.0 * c22);
  gq[8] = k * c20;
  gq[1] = gq[3] = k * 3.0 * s22;
  gq[2] = gq[6] = k * 1.5 * c21;
  gq[5] = gq[7] = k * 1.5 * s21;
}

static void mat_vec(const double m[9], const double v[3], double out[3]) {
  for (int i = 0; i < 3; ++i) {
    out[i] = m[3 * i] * v[0] + m[3 * i + 1] * v[1] + m[3 * i + 2] * v[2];
  }
}

myb_status myb_figure_figure(double gm_a, const double gq_a[9], double gm_b, const double gq_b[9],
                             const double r[3], double a_a[3], double a_b[3]) {
  MYB_CHECK(gq_a != NULL && gq_b != NULL && r != NULL);
  const double rr = myb_vec3_norm(r);
  if (!(rr > 0.0) || !isfinite(rr) || !(gm_a > 0.0) || !(gm_b > 0.0) || !isfinite(gm_a) ||
      !isfinite(gm_b)) {
    return MYB_E_DOMAIN;
  }
  const double n[3] = {r[0] / rr, r[1] / rr, r[2] / rr};
  double an[3];
  double bn[3];
  mat_vec(gq_a, n, an);
  mat_vec(gq_b, n, bn);
  const double nan_ = myb_vec3_dot(n, an);
  const double nbn = myb_vec3_dot(n, bn);
  /* (A B + B A) n and n.A B.n */
  double abn[3];
  double ban[3];
  mat_vec(gq_a, bn, abn);
  mat_vec(gq_b, an, ban);
  const double nabn = myb_vec3_dot(n, abn);
  double tr = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      tr += gq_a[3 * i + j] * gq_b[3 * j + i];
    }
  }
  const double f = 105.0 * nan_ * nbn - 60.0 * nabn + 6.0 * tr;
  double g[3];
  for (int i = 0; i < 3; ++i) {
    g[i] = 210.0 * (an[i] * nbn + nan_ * bn[i]) - 60.0 * (abn[i] + ban[i]);
  }
  const double ng = myb_vec3_dot(n, g);
  /* G^2 QA QB / 36 r^6 [g - (n.g + 5 f) n]: per unit mass of B divide by G M_B, of A by G M_A */
  const double r6 = rr * rr * rr * rr * rr * rr;
  double force[3];
  for (int i = 0; i < 3; ++i) {
    force[i] = (g[i] - (ng + 5.0 * f) * n[i]) / (36.0 * r6);
  }
  if (a_b != NULL) {
    for (int i = 0; i < 3; ++i) {
      a_b[i] += force[i] / gm_b;
    }
  }
  if (a_a != NULL) {
    for (int i = 0; i < 3; ++i) {
      a_a[i] -= force[i] / gm_a;
    }
  }
  return MYB_OK;
}

void myb_quadrupole_rotate(const double m[9], const double gq_body[9], double gq[9]) {
  MYB_CHECK(m != NULL && gq_body != NULL && gq != NULL);
  /* gq = M^T Q_body M for M: inertial -> body (rows) */
  double t[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      t[3 * i + j] = gq_body[3 * i] * m[j] + gq_body[3 * i + 1] * m[3 + j] +
                     gq_body[3 * i + 2] * m[6 + j];
    }
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      gq[3 * i + j] = m[i] * t[j] + m[3 + i] * t[3 + j] + m[6 + i] * t[6 + j];
    }
  }
}
