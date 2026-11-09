/* SPDX-License-Identifier: BSD-3-Clause */
#include <math.h>

#include "core/base.h"
#include "myb/frames.h"
#include "num/vec3.h"

/* Unit vector u = a / |a| and its derivative from da/dt: (da - u (u . da)) / |a|. */
static int unit(const double a[3], const double da[3], double u[3], double du[3]) {
  const double n = myb_vec3_norm(a);
  if (!(n > 0.0) || !isfinite(n)) {
    return 0;
  }
  for (int i = 0; i < 3; ++i) {
    u[i] = a[i] / n;
  }
  const double along = myb_vec3_dot(u, da);
  for (int i = 0; i < 3; ++i) {
    du[i] = (da[i] - u[i] * along) / n;
  }
  return 1;
}

static void cross_rate(const double a[3], const double da[3], const double b[3],
                       const double db[3], double c[3], double dc[3]) {
  double t1[3];
  double t2[3];
  myb_vec3_cross(a, b, c);
  myb_vec3_cross(da, b, t1);
  myb_vec3_cross(a, db, t2);
  for (int i = 0; i < 3; ++i) {
    dc[i] = t1[i] + t2[i];
  }
}

static void set_rows(myb_rotation* rotation, const double* rows[3], const double* rates[3],
                     const double signs[3]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      rotation->m[3 * i + j] = signs[i] * rows[i][j];
      rotation->dm[3 * i + j] = signs[i] * rates[i][j];
    }
  }
}

myb_status myb_orbital_rotation(myb_orbital_frame kind, const double state[6],
                                const double acceleration[3], myb_rotation* rotation) {
  MYB_CHECK(state != NULL && acceleration != NULL && rotation != NULL);
  MYB_CHECK(kind == MYB_ORBITAL_RTN || kind == MYB_ORBITAL_VNC || kind == MYB_ORBITAL_LVLH);
  const double* r = state;
  const double* v = state + 3;
  if (!myb_vec3_finite(r) || !myb_vec3_finite(v) || !myb_vec3_finite(acceleration)) {
    return MYB_E_DOMAIN;
  }
  /* angular momentum h = r x v, dh/dt = r x a */
  double h[3];
  double dh[3];
  myb_vec3_cross(r, v, h);
  myb_vec3_cross(r, acceleration, dh);
  double n[3];
  double dn[3];
  if (!unit(h, dh, n, dn)) {
    return MYB_E_DOMAIN;
  }
  if (kind == MYB_ORBITAL_VNC) {
    double vv[3];
    double dvv[3];
    if (!unit(v, acceleration, vv, dvv)) {
      return MYB_E_DOMAIN;
    }
    double c[3];
    double dc[3];
    cross_rate(vv, dvv, n, dn, c, dc);
    const double* rows[3] = {vv, n, c};
    const double* rates[3] = {dvv, dn, dc};
    const double signs[3] = {1.0, 1.0, 1.0};
    set_rows(rotation, rows, rates, signs);
    return MYB_OK;
  }
  double rr[3];
  double drr[3];
  if (!unit(r, v, rr, drr)) {
    return MYB_E_DOMAIN;
  }
  double t[3];
  double dt[3];
  cross_rate(n, dn, rr, drr, t, dt);
  if (kind == MYB_ORBITAL_RTN) {
    const double* rows[3] = {rr, t, n};
    const double* rates[3] = {drr, dt, dn};
    const double signs[3] = {1.0, 1.0, 1.0};
    set_rows(rotation, rows, rates, signs);
  } else {
    /* LVLH: x = T, y = -N, z = -R */
    const double* rows[3] = {t, n, rr};
    const double* rates[3] = {dt, dn, drr};
    const double signs[3] = {1.0, -1.0, -1.0};
    set_rows(rotation, rows, rates, signs);
  }
  return MYB_OK;
}
