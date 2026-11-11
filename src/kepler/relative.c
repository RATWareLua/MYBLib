/* SPDX-License-Identifier: BSD-3-Clause */
/* Relative motion (myb/relative.h): the turning frame of the chief, the linear map of the
 * relative state, rendezvous, relative orbital elements. */
#include <math.h> /* sqrt, fabs */
#include <string.h>

#include "core/check.h"
#include "myb/cmath.h"
#include "myb/num.h"
#include "myb/relative.h"
#include "num/vec3.h"

static const double two_pi = 6.283185307179586476925286766559;

/* The rows R, T, N of the turning frame of the chief and its rate of turn about N. */
static int frame_of(const myb_pv* chief, double c[9], double* rate) {
  const double* r = chief->r;
  const double* v = chief->v;
  double h[3];
  myb_vec3_cross(r, v, h);
  const double rn = myb_vec3_norm(r);
  const double hn = myb_vec3_norm(h);
  if (!(rn > 0.0) || !(hn > 0.0) || !isfinite(rn) || !isfinite(hn)) {
    return 0;
  }
  for (int i = 0; i < 3; ++i) {
    c[i] = r[i] / rn;
    c[6 + i] = h[i] / hn;
  }
  myb_vec3_cross(c + 6, c, c + 3);
  *rate = hn / (rn * rn);
  return 1;
}

static void rows_times(const double c[9], const double v[3], double out[3]) {
  for (int i = 0; i < 3; ++i) {
    out[i] = c[3 * i] * v[0] + c[3 * i + 1] * v[1] + c[3 * i + 2] * v[2];
  }
}

static void columns_times(const double c[9], const double v[3], double out[3]) {
  for (int i = 0; i < 3; ++i) {
    out[i] = c[i] * v[0] + c[3 + i] * v[1] + c[6 + i] * v[2];
  }
}

myb_status myb_relative_from_inertial(const myb_pv* chief, const myb_pv* deputy,
                                      double relative[6]) {
  MYB_CHECK(chief != NULL && deputy != NULL && relative != NULL);
  double c[9];
  double rate = 0.0;
  if (!frame_of(chief, c, &rate) || !myb_vec3_finite(deputy->r) || !myb_vec3_finite(deputy->v)) {
    return MYB_E_DOMAIN;
  }
  double d[3];
  double dv[3];
  for (int i = 0; i < 3; ++i) {
    d[i] = deputy->r[i] - chief->r[i];
    dv[i] = deputy->v[i] - chief->v[i];
  }
  rows_times(c, d, relative);
  rows_times(c, dv, relative + 3);
  /* less omega x rho, omega = rate N */
  relative[3] += rate * relative[1];
  relative[4] -= rate * relative[0];
  return MYB_OK;
}

myb_status myb_relative_to_inertial(const myb_pv* chief, const double relative[6],
                                    myb_pv* deputy) {
  MYB_CHECK(chief != NULL && deputy != NULL && relative != NULL);
  double c[9];
  double rate = 0.0;
  if (!frame_of(chief, c, &rate) || !myb_vec3_finite(relative) ||
      !myb_vec3_finite(relative + 3)) {
    return MYB_E_DOMAIN;
  }
  const double turning[3] = {relative[3] - rate * relative[1], relative[4] + rate * relative[0],
                             relative[5]};
  double d[3];
  double dv[3];
  columns_times(c, relative, d);
  columns_times(c, turning, dv);
  for (int i = 0; i < 3; ++i) {
    deputy->r[i] = chief->r[i] + d[i];
    deputy->v[i] = chief->v[i] + dv[i];
  }
  return MYB_OK;
}

/* The 6 x 6 map of an inertial variation into the turning frame (inverse == 0), or back. */
static void frame_map(const double c[9], double rate, int inverse, double t[36]) {
  memset(t, 0, 36 * sizeof *t);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double e = inverse ? c[3 * j + i] : c[3 * i + j];
      t[6 * i + j] = e;
      t[6 * (3 + i) + 3 + j] = e;
    }
  }
  /* -[omega x] C below the diagonal (to the frame), C^T [omega x] (back);
   * [omega x] = rate [[0, -1, 0], [1, 0, 0], [0, 0, 0]] */
  for (int j = 0; j < 3; ++j) {
    if (!inverse) {
      t[6 * 3 + j] = rate * c[3 + j];  /* row R of the velocity: +rate T */
      t[6 * 4 + j] = -rate * c[j];     /* row T: -rate R */
    } else {
      t[6 * (3 + j) + 0] = rate * c[3 + j]; /* C^T [omega x]: column R gets rate T */
      t[6 * (3 + j) + 1] = -rate * c[j];    /* column T gets -rate R */
    }
  }
}

static void multiply6(const double a[36], const double b[36], double out[36]) {
  double t[36];
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double s = 0.0;
      for (int k = 0; k < 6; ++k) {
        s += a[6 * i + k] * b[6 * k + j];
      }
      t[6 * i + j] = s;
    }
  }
  memcpy(out, t, sizeof t);
}

myb_status myb_relative_stm(double mu, const myb_pv* chief, double dt, double phi[36]) {
  MYB_CHECK(chief != NULL && phi != NULL);
  double c0[9];
  double c1[9];
  double rate0 = 0.0;
  double rate1 = 0.0;
  myb_pv end;
  double kepler[36];
  if (!frame_of(chief, c0, &rate0)) {
    return MYB_E_DOMAIN;
  }
  const myb_status st = myb_kepler_stm(mu, chief, dt, &end, kepler);
  if (st != MYB_OK) {
    return st;
  }
  if (!frame_of(&end, c1, &rate1)) {
    return MYB_E_DOMAIN;
  }
  double to[36];
  double back[36];
  frame_map(c1, rate1, 0, to);
  frame_map(c0, rate0, 1, back);
  multiply6(kepler, back, phi);
  multiply6(to, phi, phi);
  return MYB_OK;
}

void myb_cw_stm(double n, double dt, double phi[36]) {
  MYB_CHECK(phi != NULL && n != 0.0);
  const double nt = n * dt;
  const double s = myb_sin(nt);
  const double c = myb_cos(nt);
  memset(phi, 0, 36 * sizeof *phi);
  phi[0] = 4.0 - 3.0 * c;
  phi[3] = s / n;
  phi[4] = 2.0 * (1.0 - c) / n;
  phi[6] = 6.0 * (s - nt);
  phi[7] = 1.0;
  phi[9] = -2.0 * (1.0 - c) / n;
  phi[10] = (4.0 * s - 3.0 * nt) / n;
  phi[14] = c;
  phi[17] = s / n;
  phi[18] = 3.0 * n * s;
  phi[21] = c;
  phi[22] = 2.0 * s;
  phi[24] = -6.0 * n * (1.0 - c);
  phi[27] = -2.0 * s;
  phi[28] = 4.0 * c - 3.0;
  phi[32] = -n * s;
  phi[35] = c;
}

myb_status myb_relative_rendezvous(double mu, const myb_pv* chief, double dt,
                                   const double relative[6], const double target[3],
                                   const double target_velocity[3], double dv0[3],
                                   double dv1[3]) {
  MYB_CHECK(chief != NULL && relative != NULL && target != NULL && dv0 != NULL && dv1 != NULL);
  double phi[36];
  const myb_status st = myb_relative_stm(mu, chief, dt, phi);
  if (st != MYB_OK) {
    return st;
  }
  /* target = phi_rr r0 + phi_rv v: the velocity needed at t0 */
  double matrix[9];
  double needed[3];
  double scale = 0.0;
  for (int i = 0; i < 3; ++i) {
    needed[i] = target[i];
    for (int j = 0; j < 3; ++j) {
      matrix[3 * i + j] = phi[6 * i + 3 + j];
      scale = fmax(scale, fabs(matrix[3 * i + j]));
      needed[i] -= phi[6 * i + j] * relative[j];
    }
  }
  /* the position must depend on the velocity in every direction: the determinant against
   * the size of the block */
  const double det = matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7]) -
                     matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6]) +
                     matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
  if (!(fabs(det) > 1e-9 * scale * scale * scale) || !myb_vec3_finite(needed)) {
    return MYB_E_DOMAIN;
  }
  if (myb_lu_solve(3, matrix, needed, 1) != MYB_OK) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    dv0[i] = needed[i] - relative[3 + i];
    double arriving = 0.0;
    for (int j = 0; j < 3; ++j) {
      arriving += phi[6 * (3 + i) + j] * relative[j] + phi[6 * (3 + i) + 3 + j] * needed[j];
    }
    dv1[i] = (target_velocity != NULL ? target_velocity[i] : 0.0) - arriving;
  }
  return MYB_OK;
}

/* ---- Relative orbital elements ---------------------------------------------------------------- */

static double wrap(double angle) {
  /* into (-pi, pi]: differences of angles of nearby orbits */
  const double turns = floor((angle + 0.5 * two_pi) / two_pi);
  return angle - turns * two_pi;
}

/* a, e, i, raan, argp and the mean argument of latitude of an ellipse. */
static int mean_elements(double mu, const myb_pv* pv, double out[6]) {
  myb_elements el;
  if (myb_elements_from_pv(mu, pv, &el) != MYB_OK || !(el.e < 1.0)) {
    return 0;
  }
  double eccentric = 0.0;
  double mean = 0.0;
  if (myb_anomaly_eccentric_from_true(el.nu, el.e, &eccentric) != MYB_OK ||
      myb_anomaly_mean_from_eccentric(eccentric, el.e, &mean) != MYB_OK) {
    return 0;
  }
  out[0] = el.p / (1.0 - el.e * el.e);
  out[1] = el.e;
  out[2] = el.i;
  out[3] = el.raan;
  out[4] = el.argp;
  out[5] = el.argp + mean;
  return 1;
}

myb_status myb_roe_from_states(double mu, const myb_pv* chief, const myb_pv* deputy,
                               double roe[6]) {
  MYB_CHECK(chief != NULL && deputy != NULL && roe != NULL);
  double c[6];
  double d[6];
  if (!mean_elements(mu, chief, c) || !mean_elements(mu, deputy, d)) {
    return MYB_E_DOMAIN;
  }
  const double dnode = wrap(d[3] - c[3]);
  roe[0] = (d[0] - c[0]) / c[0];
  roe[1] = wrap(d[5] - c[5]) + dnode * myb_cos(c[2]);
  roe[2] = d[1] * myb_cos(d[4]) - c[1] * myb_cos(c[4]);
  roe[3] = d[1] * myb_sin(d[4]) - c[1] * myb_sin(c[4]);
  roe[4] = d[2] - c[2];
  roe[5] = dnode * myb_sin(c[2]);
  return MYB_OK;
}

void myb_roe_to_position(double a, double u, const double roe[6], double position[3]) {
  MYB_CHECK(roe != NULL && position != NULL);
  const double s = myb_sin(u);
  const double c = myb_cos(u);
  position[0] = a * (roe[0] - roe[2] * c - roe[3] * s);
  position[1] = a * (roe[1] + 2.0 * roe[2] * s - 2.0 * roe[3] * c);
  position[2] = a * (roe[4] * s - roe[5] * c);
}
