/* SPDX-License-Identifier: BSD-3-Clause */
/* Quaternions, Euler angles, principal axes, rigid-body equations and torques. */
#include <math.h>
#include <string.h>

#include "core/check.h"
#include "frames/rotation.h"
#include "myb/attitude.h"
#include "myb/cmath.h"
#include "myb/num.h"
#include "num/vec3.h"

/* ---- Quaternions ---------------------------------------------------------------------------- */

void myb_quat_multiply(const double a[4], const double b[4], double out[4]) {
  MYB_CHECK(a != NULL && b != NULL && out != NULL);
  const double w = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
  const double x = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2];
  const double y = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1];
  const double z = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0];
  out[0] = w;
  out[1] = x;
  out[2] = y;
  out[3] = z;
}

void myb_quat_conjugate(const double q[4], double out[4]) {
  MYB_CHECK(q != NULL && out != NULL);
  out[0] = q[0];
  out[1] = -q[1];
  out[2] = -q[2];
  out[3] = -q[3];
}

double myb_quat_normalize(double q[4]) {
  MYB_CHECK(q != NULL);
  const double n = myb_hypot(myb_hypot(q[0], q[1]), myb_hypot(q[2], q[3]));
  if (n > 0.0) {
    for (int i = 0; i < 4; ++i) {
      q[i] /= n;
    }
  }
  return n;
}

void myb_quat_to_matrix(const double q[4], double m[9]) {
  MYB_CHECK(q != NULL && m != NULL);
  const double w = q[0];
  const double x = q[1];
  const double y = q[2];
  const double z = q[3];
  const double n2 = w * w + x * x + y * y + z * z;
  MYB_CHECK(n2 > 0.0);
  const double s = 1.0 / n2;
  m[0] = (w * w + x * x - y * y - z * z) * s;
  m[1] = 2.0 * (x * y + w * z) * s;
  m[2] = 2.0 * (x * z - w * y) * s;
  m[3] = 2.0 * (x * y - w * z) * s;
  m[4] = (w * w - x * x + y * y - z * z) * s;
  m[5] = 2.0 * (y * z + w * x) * s;
  m[6] = 2.0 * (x * z + w * y) * s;
  m[7] = 2.0 * (y * z - w * x) * s;
  m[8] = (w * w - x * x - y * y + z * z) * s;
}

void myb_quat_from_matrix(const double m[9], double q[4]) {
  MYB_CHECK(m != NULL && q != NULL);
  const double tr = m[0] + m[4] + m[8];
  double w;
  double x;
  double y;
  double z;
  /* 4w^2 = 1 + tr, 4x^2 = 1 + m00 - m11 - m22, ...; m12 - m21 = 4wx, m20 - m02 = 4wy,
   * m01 - m10 = 4wz, m01 + m10 = 4xy, m02 + m20 = 4xz, m12 + m21 = 4yz */
  if (tr >= m[0] && tr >= m[4] && tr >= m[8]) {
    const double f = 2.0 * sqrt(1.0 + tr);
    w = 0.25 * f;
    x = (m[5] - m[7]) / f;
    y = (m[6] - m[2]) / f;
    z = (m[1] - m[3]) / f;
  } else if (m[0] >= m[4] && m[0] >= m[8]) {
    const double f = 2.0 * sqrt(1.0 + m[0] - m[4] - m[8]);
    x = 0.25 * f;
    w = (m[5] - m[7]) / f;
    y = (m[1] + m[3]) / f;
    z = (m[2] + m[6]) / f;
  } else if (m[4] >= m[8]) {
    const double f = 2.0 * sqrt(1.0 - m[0] + m[4] - m[8]);
    y = 0.25 * f;
    w = (m[6] - m[2]) / f;
    x = (m[1] + m[3]) / f;
    z = (m[5] + m[7]) / f;
  } else {
    const double f = 2.0 * sqrt(1.0 - m[0] - m[4] + m[8]);
    z = 0.25 * f;
    w = (m[1] - m[3]) / f;
    x = (m[2] + m[6]) / f;
    y = (m[5] + m[7]) / f;
  }
  const double sign = w < 0.0 ? -1.0 : 1.0;
  q[0] = sign * w;
  q[1] = sign * x;
  q[2] = sign * y;
  q[3] = sign * z;
  myb_quat_normalize(q);
}

void myb_quat_rotate(const double q[4], const double v[3], double out[3]) {
  double m[9];
  myb_quat_to_matrix(q, m);
  myb_internal_mat3_vec(m, v, out);
}

void myb_quat_from_rotation_vector(const double phi[3], double q[4]) {
  MYB_CHECK(phi != NULL && q != NULL);
  const double a = myb_vec3_norm(phi);
  double c;
  double f; /* sin(a/2) / a */
  if (a < 1e-3) {
    const double a2 = a * a;
    c = 1.0 - a2 / 8.0 + a2 * a2 / 384.0 - a2 * a2 * a2 / 46080.0;
    f = 0.5 - a2 / 48.0 + a2 * a2 / 3840.0 - a2 * a2 * a2 / 645120.0;
  } else {
    c = myb_cos(0.5 * a);
    f = myb_sin(0.5 * a) / a;
  }
  q[0] = c;
  q[1] = f * phi[0];
  q[2] = f * phi[1];
  q[3] = f * phi[2];
}

void myb_quat_to_rotation_vector(const double q[4], double phi[3]) {
  MYB_CHECK(q != NULL && phi != NULL);
  const double sign = q[0] < 0.0 ? -1.0 : 1.0;
  const double w = sign * q[0];
  const double v[3] = {sign * q[1], sign * q[2], sign * q[3]};
  const double s = myb_vec3_norm(v);
  if (s == 0.0) {
    phi[0] = phi[1] = phi[2] = 0.0;
    return;
  }
  const double f = 2.0 * myb_atan2(s, w) / s;
  for (int i = 0; i < 3; ++i) {
    phi[i] = f * v[i];
  }
}

void myb_quat_slerp(const double a[4], const double b[4], double t, double out[4]) {
  MYB_CHECK(a != NULL && b != NULL && out != NULL);
  double bb[4] = {b[0], b[1], b[2], b[3]};
  double dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
  if (dot < 0.0) {
    for (int i = 0; i < 4; ++i) {
      bb[i] = -bb[i];
    }
  }
  double dm = 0.0;
  double dp = 0.0;
  for (int i = 0; i < 4; ++i) {
    dm += (a[i] - bb[i]) * (a[i] - bb[i]);
    dp += (a[i] + bb[i]) * (a[i] + bb[i]);
  }
  /* angle between the unit 4-vectors, without acos */
  const double theta = 2.0 * myb_atan2(sqrt(dm), sqrt(dp));
  double wa;
  double wb;
  if (theta < 1e-6) {
    wa = 1.0 - t;
    wb = t;
  } else {
    const double s = myb_sin(theta);
    wa = myb_sin((1.0 - t) * theta) / s;
    wb = myb_sin(t * theta) / s;
  }
  for (int i = 0; i < 4; ++i) {
    out[i] = wa * a[i] + wb * bb[i];
  }
  myb_quat_normalize(out);
}

void myb_quat_derivative(const double q[4], const double omega[3], double dq[4]) {
  MYB_CHECK(q != NULL && omega != NULL && dq != NULL);
  const double w[4] = {0.0, 0.5 * omega[0], 0.5 * omega[1], 0.5 * omega[2]};
  myb_quat_multiply(q, w, dq);
}

/* ---- Euler angles --------------------------------------------------------------------------- */

static int valid_sequence(const int axes[3]) {
  for (int i = 0; i < 3; ++i) {
    if (axes[i] < 1 || axes[i] > 3) {
      return 0;
    }
  }
  return axes[0] != axes[1] && axes[1] != axes[2];
}

myb_status myb_euler_to_matrix(const int axes[3], const double angles[3], double m[9]) {
  MYB_CHECK(axes != NULL && angles != NULL && m != NULL);
  if (!valid_sequence(axes)) {
    return MYB_E_DOMAIN;
  }
  double r[9];
  double dr[9];
  double w[9];
  myb_internal_rot_axis(axes[0] - 1, angles[0], m, dr);
  for (int s = 1; s < 3; ++s) {
    myb_internal_rot_axis(axes[s] - 1, angles[s], r, dr);
    myb_internal_mat3_mul(r, m, w);
    memcpy(m, w, sizeof w);
  }
  return MYB_OK;
}

/* +1 for a cyclic order of three distinct axes (0-based), -1 otherwise. */
static double parity(int p, int q) {
  return (q - p + 3) % 3 == 1 ? 1.0 : -1.0;
}

myb_status myb_euler_from_matrix(const int axes[3], const double m[9], double angles[3]) {
  MYB_CHECK(axes != NULL && m != NULL && angles != NULL);
  if (!valid_sequence(axes)) {
    return MYB_E_DOMAIN;
  }
  const int i = axes[0] - 1;
  const int j = axes[1] - 1;
  const int k = axes[2] - 1;
#define M(r, c) m[3 * (r) + (c)]
  double a;
  double b;
  double c;
  double off; /* sin b (proper) or cos b (Tait-Bryan): 0 at the gimbal lock */
  if (i == k) {
    const int l = 3 - i - j;
    const double e = parity(i, j);
    off = myb_hypot(M(i, j), M(i, l));
    b = myb_atan2(off, M(i, i));
    a = myb_atan2(M(i, j), -e * M(i, l));
    c = myb_atan2(M(j, i), e * M(l, i));
  } else {
    const double e = parity(i, j);
    off = myb_hypot(M(k, j), M(k, k));
    b = myb_atan2(e * M(k, i), off);
    a = myb_atan2(-e * M(k, j), M(k, k));
    c = myb_atan2(-e * M(j, i), M(i, i));
  }
#undef M
  if (off < 1e-14) {
    /* gimbal lock: only a combination of a and c is defined; a = 0, c from M R_j(b)^T */
    a = 0.0;
    double rb[9];
    double drb[9];
    double rbt[9];
    double rc[9];
    myb_internal_rot_axis(j, b, rb, drb);
    myb_internal_mat3_transpose(rb, rbt);
    myb_internal_mat3_mul(m, rbt, rc);
    const int k1 = (k + 1) % 3;
    const int k2 = (k + 2) % 3;
    c = myb_atan2(rc[3 * k1 + k2], rc[3 * k1 + k1]);
  }
  angles[0] = a;
  angles[1] = b;
  angles[2] = c;
  return MYB_OK;
}

/* ---- Principal axes ------------------------------------------------------------------------- */

myb_status myb_principal_axes(const double inertia[9], double moments[3], double axes[9]) {
  MYB_CHECK(inertia != NULL && moments != NULL && axes != NULL);
  double a[9];
  double scale = 0.0;
  for (int i = 0; i < 9; ++i) {
    if (!isfinite(inertia[i])) {
      return MYB_E_DOMAIN;
    }
    a[i] = inertia[i];
    scale = fmax(scale, fabs(a[i]));
  }
  if (a[1] != a[3] || a[2] != a[6] || a[5] != a[7]) {
    return MYB_E_DOMAIN;
  }
  double v[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; /* columns: eigenvectors */
  static const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
  for (int sweep = 0; sweep < 50; ++sweep) { /* bound: quadratic convergence, ~6 sweeps */
    const double off = fabs(a[1]) + fabs(a[2]) + fabs(a[5]);
    if (off <= 0x1p-60 * scale) {
      break;
    }
    for (int r = 0; r < 3; ++r) {
      const int p = pairs[r][0];
      const int q = pairs[r][1];
      const double apq = a[3 * p + q];
      if (apq == 0.0) {
        continue;
      }
      const double tau = (a[3 * q + q] - a[3 * p + p]) / (2.0 * apq);
      const double t = (tau >= 0.0 ? 1.0 : -1.0) / (fabs(tau) + sqrt(1.0 + tau * tau));
      const double cs = 1.0 / sqrt(1.0 + t * t);
      const double sn = t * cs;
      /* A <- J^T A J with J the rotation in the (p, q) plane */
      for (int x = 0; x < 3; ++x) {
        const double ap = a[3 * x + p];
        const double aq = a[3 * x + q];
        a[3 * x + p] = cs * ap - sn * aq;
        a[3 * x + q] = sn * ap + cs * aq;
      }
      for (int x = 0; x < 3; ++x) {
        const double ap = a[3 * p + x];
        const double aq = a[3 * q + x];
        a[3 * p + x] = cs * ap - sn * aq;
        a[3 * q + x] = sn * ap + cs * aq;
      }
      a[3 * p + q] = a[3 * q + p] = 0.0;
      for (int x = 0; x < 3; ++x) {
        const double vp = v[3 * x + p];
        const double vq = v[3 * x + q];
        v[3 * x + p] = cs * vp - sn * vq;
        v[3 * x + q] = sn * vp + cs * vq;
      }
    }
  }
  /* ascending order */
  int order[3] = {0, 1, 2};
  for (int x = 0; x < 2; ++x) {
    for (int y = x + 1; y < 3; ++y) {
      if (a[4 * order[y]] < a[4 * order[x]]) {
        const int s = order[x];
        order[x] = order[y];
        order[y] = s;
      }
    }
  }
  for (int r = 0; r < 3; ++r) {
    moments[r] = a[4 * order[r]];
    for (int x = 0; x < 3; ++x) {
      axes[3 * r + x] = v[3 * x + order[r]];
    }
  }
  /* a right-handed frame */
  double z[3];
  myb_vec3_cross(axes, axes + 3, z);
  if (myb_vec3_dot(z, axes + 6) < 0.0) {
    for (int x = 0; x < 3; ++x) {
      axes[6 + x] = -axes[6 + x];
    }
  }
  return MYB_OK;
}

/* ---- Equations of motion and torques ------------------------------------------------------- */

void myb_rigid_body_rates(const double moments[3], const double q[4], const double l[3],
                          const double torque[3], double dq[4], double dl[3]) {
  MYB_CHECK(moments != NULL && q != NULL && l != NULL && dq != NULL && dl != NULL);
  const double omega[3] = {l[0] / moments[0], l[1] / moments[1], l[2] / moments[2]};
  myb_quat_derivative(q, omega, dq);
  myb_vec3_cross(l, omega, dl);
  if (torque != NULL) {
    for (int i = 0; i < 3; ++i) {
      dl[i] += torque[i];
    }
  }
}

void myb_torque_point_mass(double gm, const double r[3], const double moments[3],
                           double torque[3]) {
  MYB_CHECK(r != NULL && moments != NULL && torque != NULL);
  const double r2 = myb_vec3_dot(r, r);
  const double ir[3] = {moments[0] * r[0], moments[1] * r[1], moments[2] * r[2]};
  double c[3];
  myb_vec3_cross(r, ir, c);
  const double f = 3.0 * gm / (r2 * r2 * sqrt(r2));
  for (int i = 0; i < 3; ++i) {
    torque[i] = f * c[i];
  }
}

void myb_torque_field(double mass_ratio, const double r[3], const double a[3], double torque[3]) {
  MYB_CHECK(r != NULL && a != NULL && torque != NULL);
  double c[3];
  myb_vec3_cross(r, a, c);
  for (int i = 0; i < 3; ++i) {
    torque[i] = -mass_ratio * c[i];
  }
}

myb_status myb_torque_figure_figure(double gm_a, const double gq_a[9], const double gq_b[9],
                                    const double r[3], double torque[3]) {
  MYB_CHECK(gq_a != NULL && gq_b != NULL && r != NULL && torque != NULL);
  const double rr = myb_vec3_norm(r);
  if (!(rr > 0.0) || !isfinite(rr) || !(gm_a > 0.0)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 9; ++i) {
    if (!isfinite(gq_a[i]) || !isfinite(gq_b[i])) {
      return MYB_E_DOMAIN;
    }
  }
  const double n[3] = {r[0] / rr, r[1] / rr, r[2] / rr};
  double b[3];
  myb_internal_mat3_vec(gq_b, n, b);
  const double nbn = myb_vec3_dot(n, b);
  const double r2 = rr * rr;
  const double f = 1.0 / (36.0 * r2 * r2 * rr);
  /* G = dU/dQ_A (symmetric) of U = f [105 (n.Qa n)(n.Qb n) - 60 n.Qa Qb n + 6 tr(Qa Qb)] */
  double g[9];
  for (int x = 0; x < 3; ++x) {
    for (int y = 0; y < 3; ++y) {
      g[3 * x + y] = f * (105.0 * nbn * n[x] * n[y] - 30.0 * (n[x] * b[y] + b[x] * n[y]) +
                          6.0 * gq_b[3 * x + y]);
    }
  }
  /* dU = tr((Qa G - G Qa) W) for dQa = W Qa - Qa W, W = [dtheta x]: the torque is the axial
   * part of A = Qa G - G Qa, (A12 - A21, A20 - A02, A01 - A10) */
  double qg[9];
  double gq[9];
  myb_internal_mat3_mul(gq_a, g, qg);
  myb_internal_mat3_mul(g, gq_a, gq);
  double am[9];
  for (int x = 0; x < 9; ++x) {
    am[x] = qg[x] - gq[x];
  }
  torque[0] = (am[5] - am[7]) / gm_a;
  torque[1] = (am[6] - am[2]) / gm_a;
  torque[2] = (am[1] - am[3]) / gm_a;
  return MYB_OK;
}

myb_status myb_core_mantle_rates(const double inertia_mantle[9], const double inertia_rate[9],
                                 const double core[2], double kv, const double omega_mantle[3],
                                 const double omega_core[3], const double torque[3],
                                 double rate_mantle[3], double rate_core[3]) {
  MYB_CHECK(inertia_mantle != NULL && core != NULL && omega_mantle != NULL &&
            omega_core != NULL && torque != NULL && rate_mantle != NULL && rate_core != NULL);
  const double ac = core[0];
  const double cc = core[1];
  if (!(ac > 0.0) || !(cc > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double* wm = omega_mantle;
  const double* wc = omega_core;
  /* N_cmb */
  double ncmb[3];
  const double zc = wc[2];
  const double zxw[3] = {-wc[1], wc[0], 0.0}; /* z x omega_c */
  for (int i = 0; i < 3; ++i) {
    ncmb[i] = kv * (wc[i] - wm[i]) + (cc - ac) * zc * zxw[i];
  }
  /* mantle */
  double l[3];
  myb_internal_mat3_vec(inertia_mantle, wm, l);
  double gyro[3];
  myb_vec3_cross(wm, l, gyro);
  double rhs[3];
  double di[3] = {0.0, 0.0, 0.0};
  if (inertia_rate != NULL) {
    myb_internal_mat3_vec(inertia_rate, wm, di);
  }
  for (int i = 0; i < 3; ++i) {
    rhs[i] = torque[i] + ncmb[i] - di[i] - gyro[i];
  }
  double a[9];
  memcpy(a, inertia_mantle, sizeof a);
  if (myb_lu_solve(3, a, rhs, 1) != MYB_OK) {
    return MYB_E_DOMAIN;
  }
  memcpy(rate_mantle, rhs, 3 * sizeof *rhs);
  /* core: I_c diagonal (A_c, A_c, C_c) */
  const double lc[3] = {ac * wc[0], ac * wc[1], cc * wc[2]};
  double gc[3];
  myb_vec3_cross(wm, lc, gc);
  rate_core[0] = (-gc[0] - ncmb[0]) / ac;
  rate_core[1] = (-gc[1] - ncmb[1]) / ac;
  rate_core[2] = (-gc[2] - ncmb[2]) / cc;
  return MYB_OK;
}
