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

myb_status myb_pn_binary_acceleration(double gm1, double gm2, double c, unsigned orders,
                                      const double x[3], const double v[3], double a[3]) {
  MYB_CHECK(x != NULL && v != NULL && a != NULL);
  const double r2 = myb_vec3_dot(x, x);
  const double r = sqrt(r2);
  if (!(gm1 > 0.0) || !(gm2 > 0.0) || !(c > 0.0) || !(r > 0.0) || !isfinite(r) ||
      !isfinite(gm1 + gm2 + c) || !isfinite(myb_vec3_dot(v, v))) {
    return MYB_E_DOMAIN;
  }
  const double gm = gm1 + gm2;
  const double nu = gm1 * gm2 / (gm * gm);
  const double nu2 = nu * nu;
  const double n[3] = {x[0] / r, x[1] / r, x[2] / r};
  const double rd = myb_vec3_dot(n, v);
  const double v2 = myb_vec3_dot(v, v);
  const double u = gm / r; /* gm / r */
  const double c2 = c * c;
  double aa = 0.0;
  double bb = 0.0;
  if (orders & MYB_PN_1) {
    aa += (-1.5 * rd * rd * nu + v2 + 3.0 * nu * v2 - u * (4.0 + 2.0 * nu)) / c2;
    bb += (-4.0 * rd + 2.0 * rd * nu) / c2;
  }
  if (orders & MYB_PN_2) {
    const double rd2 = rd * rd;
    const double a2 = 15.0 / 8.0 * rd2 * rd2 * nu - 45.0 / 8.0 * rd2 * rd2 * nu2 -
                      4.5 * rd2 * nu * v2 + 6.0 * rd2 * nu2 * v2 + 3.0 * nu * v2 * v2 -
                      4.0 * nu2 * v2 * v2 +
                      u * (-2.0 * rd2 - 25.0 * rd2 * nu - 2.0 * rd2 * nu2 - 6.5 * nu * v2 +
                           2.0 * nu2 * v2) +
                      u * u * (9.0 + 87.0 / 4.0 * nu);
    const double b2 = 4.5 * rd2 * rd * nu + 3.0 * rd2 * rd * nu2 - 7.5 * rd * nu * v2 -
                      2.0 * rd * nu2 * v2 + u * (2.0 * rd + 20.5 * rd * nu + 4.0 * rd * nu2);
    aa += a2 / (c2 * c2);
    bb += b2 / (c2 * c2);
  }
  if (orders & MYB_PN_25) {
    const double c5 = c2 * c2 * c;
    aa += (-24.0 / 5.0 * rd * nu * v2 * u - 136.0 / 15.0 * rd * nu * u * u) / c5;
    bb += (8.0 / 5.0 * nu * v2 * u + 24.0 / 5.0 * nu * u * u) / c5;
  }
  const double f = gm / r2;
  for (int i = 0; i < 3; ++i) {
    a[i] = -f * ((1.0 + aa) * n[i] + bb * v[i]);
  }
  return MYB_OK;
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
