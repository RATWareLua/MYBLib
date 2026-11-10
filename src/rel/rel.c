/* SPDX-License-Identifier: BSD-3-Clause */
#include <math.h>

#include "core/base.h"
#include "myb/cmath.h"
#include "myb/rel.h"
#include "num/vec3.h"

myb_status myb_ppn_validate(const myb_ppn* ppn) {
  MYB_CHECK(ppn != NULL);
  if (!(ppn->c > 0.0 && ppn->c < 1e300) || !isfinite(ppn->beta) || !isfinite(ppn->gamma)) {
    return MYB_E_DOMAIN;
  }
  return MYB_OK;
}

static int usable(const myb_ppn* ppn, const double* a, const double* b) {
  return myb_ppn_validate(ppn) == MYB_OK && myb_vec3_finite(a) && myb_vec3_finite(b);
}

myb_status myb_rel_schwarzschild(const myb_ppn* ppn, double gm, const double r[3],
                                 const double v[3], double a[3]) {
  MYB_CHECK(ppn != NULL && r != NULL && v != NULL && a != NULL);
  const double rn = myb_vec3_norm(r);
  if (!usable(ppn, r, v) || !isfinite(gm) || !(rn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double c2 = ppn->c * ppn->c;
  const double k = gm / (c2 * rn * rn * rn);
  const double radial = 2.0 * (ppn->beta + ppn->gamma) * gm / rn - ppn->gamma * myb_vec3_dot(v, v);
  const double along = 2.0 * (1.0 + ppn->gamma) * myb_vec3_dot(r, v);
  for (int i = 0; i < 3; ++i) {
    a[i] = k * (radial * r[i] + along * v[i]);
  }
  return MYB_OK;
}

myb_status myb_rel_lense_thirring(const myb_ppn* ppn, double gm, const double j[3],
                                  const double r[3], const double v[3], double a[3]) {
  MYB_CHECK(ppn != NULL && j != NULL && r != NULL && v != NULL && a != NULL);
  const double rn = myb_vec3_norm(r);
  if (!usable(ppn, r, v) || !myb_vec3_finite(j) || !isfinite(gm) || !(rn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double c2 = ppn->c * ppn->c;
  const double k = (1.0 + ppn->gamma) * gm / (c2 * rn * rn * rn);
  double rxv[3];
  double vxj[3];
  myb_vec3_cross(r, v, rxv);
  myb_vec3_cross(v, j, vxj);
  const double s = 3.0 * myb_vec3_dot(r, j) / (rn * rn);
  for (int i = 0; i < 3; ++i) {
    a[i] = k * (s * rxv[i] + vxj[i]);
  }
  return MYB_OK;
}

myb_status myb_rel_de_sitter(const myb_ppn* ppn, double gm_external, const double R[3],
                             const double V[3], const double v[3], double a[3]) {
  MYB_CHECK(ppn != NULL && R != NULL && V != NULL && v != NULL && a != NULL);
  const double rn = myb_vec3_norm(R);
  if (!usable(ppn, R, V) || !myb_vec3_finite(v) || !isfinite(gm_external) || !(rn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double c2 = ppn->c * ppn->c;
  const double k = -gm_external / (c2 * rn * rn * rn);
  double field[3];
  double inner[3];
  for (int i = 0; i < 3; ++i) {
    field[i] = k * R[i];
  }
  myb_vec3_cross(V, field, inner);
  myb_vec3_cross(inner, v, a);
  for (int i = 0; i < 3; ++i) {
    a[i] *= 1.0 + 2.0 * ppn->gamma;
  }
  return MYB_OK;
}

/* [a x]: the matrix of the cross product, ([a x] b) = a x b. */
static void cross_matrix(const double a[3], double m[9]) {
  m[0] = 0.0;
  m[1] = -a[2];
  m[2] = a[1];
  m[3] = a[2];
  m[4] = 0.0;
  m[5] = -a[0];
  m[6] = -a[1];
  m[7] = a[0];
  m[8] = 0.0;
}

myb_status myb_rel_schwarzschild_jacobian(const myb_ppn* ppn, double gm, const double r[3],
                                          const double v[3], double dadr[9], double dadv[9]) {
  MYB_CHECK(ppn != NULL && r != NULL && v != NULL);
  const double rn = myb_vec3_norm(r);
  if (!usable(ppn, r, v) || !isfinite(gm) || !(rn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  /* a = k (A r + B v), k = gm / (c^2 r^3), A = 2 (beta + gamma) gm / r - gamma v.v,
   * B = 2 (1 + gamma) r.v */
  const double c2 = ppn->c * ppn->c;
  const double k = gm / (c2 * rn * rn * rn);
  const double radial = 2.0 * (ppn->beta + ppn->gamma) * gm / rn - ppn->gamma * myb_vec3_dot(v, v);
  const double along = 2.0 * (1.0 + ppn->gamma) * myb_vec3_dot(r, v);
  const double d_radial = -2.0 * (ppn->beta + ppn->gamma) * gm / (rn * rn * rn);
  for (int i = 0; i < 3; ++i) {
    const double ai = radial * r[i] + along * v[i];
    for (int j = 0; j < 3; ++j) {
      if (dadr != NULL) {
        dadr[3 * i + j] = -3.0 * k * ai * r[j] / (rn * rn) +
                          k * ((i == j ? radial : 0.0) + r[i] * d_radial * r[j] +
                               2.0 * (1.0 + ppn->gamma) * v[i] * v[j]);
      }
      if (dadv != NULL) {
        dadv[3 * i + j] = k * (-2.0 * ppn->gamma * r[i] * v[j] + (i == j ? along : 0.0) +
                               2.0 * (1.0 + ppn->gamma) * v[i] * r[j]);
      }
    }
  }
  return MYB_OK;
}

myb_status myb_rel_lense_thirring_jacobian(const myb_ppn* ppn, double gm, const double j[3],
                                           const double r[3], const double v[3], double dadr[9],
                                           double dadv[9]) {
  MYB_CHECK(ppn != NULL && j != NULL && r != NULL && v != NULL);
  const double rn = myb_vec3_norm(r);
  if (!usable(ppn, r, v) || !myb_vec3_finite(j) || !isfinite(gm) || !(rn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  /* a = K (3 h p / r^5 + q / r^3), K = (1 + gamma) gm / c^2, h = r x v, p = r.j, q = v x j */
  const double c2 = ppn->c * ppn->c;
  const double big_k = (1.0 + ppn->gamma) * gm / c2;
  const double r2 = rn * rn;
  const double r3 = r2 * rn;
  const double r5 = r3 * r2;
  double h[3];
  double q[3];
  myb_vec3_cross(r, v, h);
  myb_vec3_cross(v, j, q);
  const double p = myb_vec3_dot(r, j);
  double vx[9];
  double rx[9];
  double jx[9];
  cross_matrix(v, vx);
  cross_matrix(r, rx);
  cross_matrix(j, jx);
  for (int i = 0; i < 3; ++i) {
    for (int k = 0; k < 3; ++k) {
      if (dadr != NULL) {
        /* d h / d r = -[v x], d p / d r = j */
        dadr[3 * i + k] = big_k * (3.0 * (h[i] * j[k] - p * vx[3 * i + k]) / r5 -
                                   15.0 * h[i] * p * r[k] / (r5 * r2) - 3.0 * q[i] * r[k] / r5);
      }
      if (dadv != NULL) {
        /* d h / d v = [r x], d q / d v = -[j x] */
        dadv[3 * i + k] = big_k * (3.0 * p * rx[3 * i + k] / r5 - jx[3 * i + k] / r3);
      }
    }
  }
  return MYB_OK;
}

myb_status myb_rel_de_sitter_jacobian(const myb_ppn* ppn, double gm_external, const double R[3],
                                      const double V[3], double dadv[9]) {
  MYB_CHECK(ppn != NULL && R != NULL && V != NULL && dadv != NULL);
  const double rn = myb_vec3_norm(R);
  if (!usable(ppn, R, V) || !isfinite(gm_external) || !(rn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  /* a = Omega x v, Omega = (1 + 2 gamma) V x (-gm R / (c^2 R^3)) */
  const double c2 = ppn->c * ppn->c;
  const double k = -gm_external / (c2 * rn * rn * rn);
  double field[3];
  double omega[3];
  for (int i = 0; i < 3; ++i) {
    field[i] = k * R[i];
  }
  myb_vec3_cross(V, field, omega);
  for (int i = 0; i < 3; ++i) {
    omega[i] *= 1.0 + 2.0 * ppn->gamma;
  }
  cross_matrix(omega, dadv);
  return MYB_OK;
}

/* (1/c^2) [ 1/2 ve (ve.r) + w_ext r + r (ae.r) - 1/2 ae r^2 ] */
static void local_correction(const myb_rel_local_frame* f, const double r[3], double out[3]) {
  const double c2 = f->c * f->c;
  const double vr = myb_vec3_dot(f->ve, r);
  const double ar = myb_vec3_dot(f->ae, r);
  const double r2 = myb_vec3_dot(r, r);
  for (int i = 0; i < 3; ++i) {
    out[i] = (0.5 * f->ve[i] * vr + f->w_ext * r[i] + r[i] * ar - 0.5 * f->ae[i] * r2) / c2;
  }
}

void myb_rel_global_to_local(const myb_rel_local_frame* frame, const double x[3],
                             double out[3]) {
  MYB_CHECK(frame != NULL && x != NULL && out != NULL && frame->c > 0.0);
  double r[3];
  for (int i = 0; i < 3; ++i) {
    r[i] = x[i] - frame->xe[i];
  }
  double d[3];
  local_correction(frame, r, d);
  for (int i = 0; i < 3; ++i) {
    out[i] = r[i] + d[i];
  }
}

void myb_rel_local_to_global(const myb_rel_local_frame* frame, const double X[3], double out[3]) {
  MYB_CHECK(frame != NULL && X != NULL && out != NULL && frame->c > 0.0);
  const double local[3] = {X[0], X[1], X[2]};
  double d[3];
  local_correction(frame, local, d);
  for (int i = 0; i < 3; ++i) {
    out[i] = frame->xe[i] + (local[i] - d[i]);
  }
}

double myb_rel_local_time_shift(const myb_rel_local_frame* frame, const double x[3]) {
  MYB_CHECK(frame != NULL && x != NULL && frame->c > 0.0);
  double r[3];
  for (int i = 0; i < 3; ++i) {
    r[i] = x[i] - frame->xe[i];
  }
  return -myb_vec3_dot(frame->ve, r) / (frame->c * frame->c);
}

double myb_rel_proper_time_rate(double c, double v2, double potential) {
  MYB_CHECK(c > 0.0);
  return 1.0 - (0.5 * v2 + potential) / (c * c);
}

/* ---- Compact binaries ---------------------------------------------------------------------- */

/* pi^2 */
#define MYB_PI2 9.8696044010893586188

/* Coefficients A and B of a = -(gm / r^2) [(1 + A) n + B v] (myb/rel.h), each order divided by
 * its power of c. u = gm / r, rd = n.v, v2 = v.v, log_r = ln(r / r0). */
static void pn_coefficients(unsigned orders, double nu, double u, double rd, double v2, double c,
                            double log_r, double* aa_out, double* bb_out) {
  const double nu2 = nu * nu;
  const double nu3 = nu2 * nu;
  const double rd2 = rd * rd;
  const double rd4 = rd2 * rd2;
  const double v4 = v2 * v2;
  const double c2 = c * c;
  const double c4 = c2 * c2;
  double aa = 0.0;
  double bb = 0.0;
  if (orders & MYB_PN_1) {
    aa += (-1.5 * rd2 * nu + v2 + 3.0 * nu * v2 - u * (4.0 + 2.0 * nu)) / c2;
    bb += (-4.0 * rd + 2.0 * rd * nu) / c2;
  }
  if (orders & MYB_PN_2) {
    const double a2 = 15.0 / 8.0 * rd4 * nu - 45.0 / 8.0 * rd4 * nu2 - 4.5 * rd2 * nu * v2 +
                      6.0 * rd2 * nu2 * v2 + 3.0 * nu * v4 - 4.0 * nu2 * v4 +
                      u * (-2.0 * rd2 - 25.0 * rd2 * nu - 2.0 * rd2 * nu2 - 6.5 * nu * v2 +
                           2.0 * nu2 * v2) +
                      u * u * (9.0 + 87.0 / 4.0 * nu);
    const double b2 = 4.5 * rd2 * rd * nu + 3.0 * rd2 * rd * nu2 - 7.5 * rd * nu * v2 -
                      2.0 * rd * nu2 * v2 + u * (2.0 * rd + 20.5 * rd * nu + 4.0 * rd * nu2);
    aa += a2 / c4;
    bb += b2 / c4;
  }
  if (orders & MYB_PN_25) {
    const double c5 = c4 * c;
    aa += (-24.0 / 5.0 * rd * nu * v2 * u - 136.0 / 15.0 * rd * nu * u * u) / c5;
    bb += (8.0 / 5.0 * nu * v2 * u + 24.0 / 5.0 * nu * u * u) / c5;
  }
  if (orders & MYB_PN_3) {
    const double rd6 = rd4 * rd2;
    const double v6 = v4 * v2;
    const double u2 = u * u;
    const double a0 = rd6 * (-35.0 / 16.0 * nu + 175.0 / 16.0 * nu2 - 175.0 / 16.0 * nu3) +
                      rd4 * v2 * (7.5 * nu - 135.0 / 4.0 * nu2 + 255.0 / 8.0 * nu3) +
                      rd2 * v4 * (-7.5 * nu + 237.0 / 8.0 * nu2 - 22.5 * nu3) +
                      v6 * (11.0 / 4.0 * nu - 49.0 / 4.0 * nu2 + 13.0 * nu3);
    const double a1 = rd4 * (79.0 * nu - 34.5 * nu2 - 30.0 * nu3) +
                      rd2 * v2 * (-121.0 * nu + 16.0 * nu2 + 20.0 * nu3) +
                      v4 * (75.0 / 4.0 * nu + 8.0 * nu2 - 10.0 * nu3);
    const double a2 = rd2 * (1.0 + 32573.0 / 168.0 * nu + 11.0 / 8.0 * nu2 - 7.0 * nu3 +
                             615.0 / 64.0 * MYB_PI2 * nu - 110.0 * nu * log_r) +
                      v2 * (-26987.0 / 840.0 * nu + nu3 - 123.0 / 64.0 * MYB_PI2 * nu +
                            22.0 * nu * log_r);
    const double a3 = -16.0 - 437.0 / 4.0 * nu - 35.5 * nu2 + 41.0 / 16.0 * MYB_PI2 * nu;
    const double b0 = rd4 * rd * (-45.0 / 8.0 * nu + 15.0 * nu2 + 3.75 * nu3) +
                      rd2 * rd * v2 * (12.0 * nu - 111.0 / 4.0 * nu2 - 12.0 * nu3) +
                      rd * v4 * (-65.0 / 8.0 * nu + 19.0 * nu2 + 6.0 * nu3);
    const double b1 = rd2 * rd * (329.0 / 6.0 * nu + 29.5 * nu2 + 18.0 * nu3) +
                      rd * v2 * (-15.0 * nu - 27.0 * nu2 - 10.0 * nu3);
    const double b2 = rd * (-4.0 - 18169.0 / 840.0 * nu + 25.0 * nu2 + 8.0 * nu3 -
                            123.0 / 32.0 * MYB_PI2 * nu + 44.0 * nu * log_r);
    const double c6 = c4 * c2;
    aa += (a0 + u * a1 + u2 * a2 + u2 * u * a3) / c6;
    bb += (b0 + u * b1 + u2 * b2) / c6;
  }
  if (orders & MYB_PN_35) {
    const double u2 = u * u;
    const double a35 = u * nu * rd *
                       (u2 * (3956.0 / 35.0 + 184.0 / 5.0 * nu) +
                        u * v2 * (692.0 / 35.0 - 724.0 / 15.0 * nu) + v4 * (366.0 / 35.0 + 12.0 * nu) +
                        u * rd2 * (294.0 / 5.0 + 376.0 / 5.0 * nu) - v2 * rd2 * (114.0 + 12.0 * nu) +
                        112.0 * rd4);
    const double b35 = u * nu *
                       (u2 * (-1060.0 / 21.0 - 104.0 / 5.0 * nu) +
                        u * v2 * (164.0 / 21.0 + 148.0 / 5.0 * nu) +
                        v4 * (-626.0 / 35.0 - 12.0 / 5.0 * nu) +
                        u * rd2 * (-82.0 / 3.0 - 848.0 / 15.0 * nu) +
                        v2 * rd2 * (678.0 / 5.0 + 12.0 / 5.0 * nu) - 120.0 * rd4);
    const double c7 = c4 * c2 * c;
    aa += a35 / c7;
    bb += b35 / c7;
  }
  *aa_out = aa;
  *bb_out = bb;
}

/* out += s * (Omega x chi) */
static void precess(const double omega[3], const double chi[3], double out[3]) {
  double t[3];
  myb_vec3_cross(omega, chi, t);
  for (int i = 0; i < 3; ++i) {
    out[i] += t[i];
  }
}

myb_status myb_pn_binary_evaluate(const myb_pn_binary* b, const double x[3], const double v[3],
                                  const double chi1[3], const double chi2[3], double a[3],
                                  double dchi1[3], double dchi2[3]) {
  MYB_CHECK(b != NULL && x != NULL && v != NULL && a != NULL);
  static const double zero[3] = {0.0, 0.0, 0.0};
  const double* s1 = chi1 != NULL ? chi1 : zero;
  const double* s2 = chi2 != NULL ? chi2 : zero;
  const double gm1 = b->gm1;
  const double gm2 = b->gm2;
  const double c = b->c;
  const unsigned orders = b->orders;
  const double r2 = myb_vec3_dot(x, x);
  const double r = sqrt(r2);
  if (!(gm1 > 0.0) || !(gm2 > 0.0) || !(c > 0.0) || !(r > 0.0) || !isfinite(r) ||
      !isfinite(gm1 + gm2 + c) || !myb_vec3_finite(v) || !myb_vec3_finite(s1) ||
      !myb_vec3_finite(s2) || (orders & ~(unsigned)MYB_PN_ALL) != 0u) {
    return MYB_E_DOMAIN;
  }
  if ((orders & MYB_PN_3) && !(b->r0 > 0.0 && b->r0 < 1e300)) {
    return MYB_E_DOMAIN;
  }
  if ((orders & MYB_PN_SS) && !(b->kappa1 >= 0.0 && b->kappa1 < 1e300 && b->kappa2 >= 0.0 &&
                                b->kappa2 < 1e300)) {
    return MYB_E_DOMAIN;
  }
  const double gm = gm1 + gm2;
  const double nu = gm1 * gm2 / (gm * gm);
  const double n[3] = {x[0] / r, x[1] / r, x[2] / r};
  const double rd = myb_vec3_dot(n, v);
  const double v2 = myb_vec3_dot(v, v);
  const double u = gm / r;
  const double log_r = (orders & MYB_PN_3) ? myb_log(r / b->r0) : 0.0;
  double aa = 0.0;
  double bb = 0.0;
  pn_coefficients(orders, nu, u, rd, v2, c, log_r, &aa, &bb);
  const double f = gm / r2;
  for (int i = 0; i < 3; ++i) {
    a[i] = -f * ((1.0 + aa) * n[i] + bb * v[i]);
  }
  double w1[3] = {0.0, 0.0, 0.0}; /* precession angular velocities of the spins */
  double w2[3] = {0.0, 0.0, 0.0};
  const double c3 = c * c * c;
  const double r3 = r2 * r;
  if (orders & MYB_PN_SO) {
    /* G S / c^2 (m^3/s) and (dm / m) G Delta / c^2 */
    double sig[3];
    double sigd[3];
    for (int i = 0; i < 3; ++i) {
      sig[i] = (gm1 * gm1 * s1[i] + gm2 * gm2 * s2[i]) / c3;
      sigd[i] = (gm1 - gm2) * (gm2 * s2[i] - gm1 * s1[i]) / c3;
    }
    double nxv[3];
    myb_vec3_cross(n, v, nxv);
    double p[3];
    double q[3];
    double t[3];
    for (int i = 0; i < 3; ++i) {
      p[i] = 2.0 * sig[i] + sigd[i];
      q[i] = 7.0 * sig[i] + 3.0 * sigd[i];
      t[i] = 3.0 * sig[i] + sigd[i];
    }
    double vxq[3];
    double nxt[3];
    myb_vec3_cross(v, q, vxq);
    myb_vec3_cross(n, t, nxt);
    const double k = 6.0 * myb_vec3_dot(nxv, p);
    for (int i = 0; i < 3; ++i) {
      a[i] += (k * n[i] - vxq[i] + 3.0 * rd * nxt[i]) / r3;
    }
    /* geodetic precession: (G mu / (c^2 r^3)) (2 + 3 m_B / (2 m_A)) x x v */
    double l[3];
    myb_vec3_cross(x, v, l);
    const double g = gm1 * gm2 / (gm * c * c * r3);
    for (int i = 0; i < 3; ++i) {
      w1[i] += g * (2.0 + 1.5 * gm2 / gm1) * l[i];
      w2[i] += g * (2.0 + 1.5 * gm1 / gm2) * l[i];
    }
  }
  if (orders & MYB_PN_SS) {
    const double n1 = myb_vec3_dot(n, s1);
    const double n2 = myb_vec3_dot(n, s2);
    const double s12 = myb_vec3_dot(s1, s2);
    const double c4 = c3 * c;
    const double r4 = r2 * r2;
    const double kss = -3.0 * gm * gm1 * gm2 / (c4 * r4);
    const double k1 = -1.5 * b->kappa1 * gm * gm1 * gm1 / (c4 * r4);
    const double k2 = -1.5 * b->kappa2 * gm * gm2 * gm2 / (c4 * r4);
    const double q1 = myb_vec3_dot(s1, s1) - 5.0 * n1 * n1;
    const double q2 = myb_vec3_dot(s2, s2) - 5.0 * n2 * n2;
    for (int i = 0; i < 3; ++i) {
      a[i] += kss * (n[i] * (s12 - 5.0 * n1 * n2) + s1[i] * n2 + s2[i] * n1) +
              k1 * (2.0 * n1 * s1[i] + q1 * n[i]) + k2 * (2.0 * n2 * s2[i] + q2 * n[i]);
    }
    const double h = 1.0 / (c3 * r3);
    const double qm1 = 3.0 * b->kappa1 * gm1 * gm2 * n1 * h;
    const double qm2 = 3.0 * b->kappa2 * gm1 * gm2 * n2 * h;
    for (int i = 0; i < 3; ++i) {
      w1[i] += gm2 * gm2 * h * (3.0 * n2 * n[i] - s2[i]) + qm1 * n[i];
      w2[i] += gm1 * gm1 * h * (3.0 * n1 * n[i] - s1[i]) + qm2 * n[i];
    }
  }
  if (dchi1 != NULL) {
    dchi1[0] = dchi1[1] = dchi1[2] = 0.0;
    precess(w1, s1, dchi1);
  }
  if (dchi2 != NULL) {
    dchi2[0] = dchi2[1] = dchi2[2] = 0.0;
    precess(w2, s2, dchi2);
  }
  return MYB_OK;
}

myb_status myb_pn_binary_acceleration(double gm1, double gm2, double c, unsigned orders,
                                      const double x[3], const double v[3], double a[3]) {
  MYB_CHECK(x != NULL && v != NULL && a != NULL);
  const myb_pn_binary b = {.gm1 = gm1, .gm2 = gm2, .c = c, .orders = orders, .r0 = (double)NAN};
  return myb_pn_binary_evaluate(&b, x, v, NULL, NULL, a, NULL, NULL);
}

double myb_peters_period_derivative(double gm1, double gm2, double c, double period, double e) {
  if (!(gm1 > 0.0) || !(gm2 > 0.0) || !(c > 0.0) || !(period > 0.0) || !(e >= 0.0 && e < 1.0)) {
    return (double)NAN;
  }
  const double two_pi = 6.283185307179586;
  const double gmc = myb_exp(0.6 * myb_log(gm1 * gm2) - 0.2 * myb_log(gm1 + gm2));
  const double e2 = e * e;
  const double one = 1.0 - e2;
  const double fe = (1.0 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e2 * e2) / (one * one * one * sqrt(one));
  const double y = two_pi * gmc / (period * c * c * c);
  return -192.0 * 3.141592653589793 / 5.0 * myb_exp(5.0 / 3.0 * myb_log(y)) * fe;
}
