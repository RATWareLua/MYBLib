/* SPDX-License-Identifier: BSD-3-Clause */
/* Torque-free rotation in closed form and the splitting integrator built on it. */
#include <math.h>
#include <string.h>

#include "core/check.h"
#include "frames/rotation.h"
#include "integ/symplectic_coefficients.h"
#include "myb/attitude.h"
#include "myb/cmath.h"
#include "myb/num.h"
#include "num/vec3.h"

enum { KIND_GENERAL = 0, KIND_STEADY = 1 };

static void rot3(int axis, double c, double s, double m[9]) {
  /* frame rotation about `axis` by the angle with cosine c and sine s */
  memset(m, 0, 9 * sizeof *m);
  const int i = (axis + 1) % 3;
  const int j = (axis + 2) % 3;
  m[4 * axis] = 1.0;
  m[4 * i] = c;
  m[3 * i + j] = s;
  m[3 * j + i] = -s;
  m[4 * j] = c;
}

/* M(t) in ordered axes for the general motion: R3(psi) R1(theta) R3(phi) Z. */
static void general_attitude(const myb_free_rotation* f, const double mo[3], double phi,
                             double m[9]) {
  const double rho = myb_hypot(mo[0], mo[1]);
  const double st = rho / f->g;
  const double ct = mo[2] / f->g;
  double sp = 0.0;
  double cp = 1.0;
  if (rho > 0.0) {
    sp = mo[0] / rho;
    cp = mo[1] / rho;
  }
  double r_psi[9];
  double r_theta[9];
  double r_phi[9];
  double w[9];
  rot3(2, cp, sp, r_psi);
  rot3(0, ct, st, r_theta);
  rot3(2, myb_cos(phi), myb_sin(phi), r_phi);
  myb_internal_mat3_mul(r_phi, f->z, w);
  myb_internal_mat3_mul(r_theta, w, m);
  myb_internal_mat3_mul(r_psi, m, w);
  memcpy(m, w, sizeof w);
}

/* Body momentum (ordered axes) and the precession angle at t. */
static void general_motion(const myb_free_rotation* f, double t, double mo[3], double* phi) {
  const double u = f->u0 + f->nu * t;
  double am;
  double sn;
  double cn;
  double dn;
  myb_jacobi(u, f->k2, &am, &sn, &cn, &dn);
  mo[0] = f->a[0] * cn;
  mo[1] = f->sign2 * f->a[1] * sn;
  mo[2] = f->sign3 * f->a[2] * dn;
  *phi = f->rate0 * t + f->rate_pi * (myb_elliptic_pi_inc(f->n, am, f->k2) - f->pi0);
}

myb_status myb_free_rotation_init(const double moments[3], const double q[4], const double l[3],
                                  myb_free_rotation* rotation) {
  MYB_CHECK(moments != NULL && q != NULL && l != NULL && rotation != NULL);
  for (int i = 0; i < 3; ++i) {
    if (!(moments[i] > 0.0) || !isfinite(moments[i]) || !isfinite(l[i])) {
      return MYB_E_DOMAIN;
    }
  }
  myb_free_rotation* f = rotation;
  memset(f, 0, sizeof *f);
  double m0[9];
  myb_quat_to_matrix(q, m0);

  /* steady rotation: omega parallel to L (a principal axis, an equal-moment plane, a sphere) */
  const double om[3] = {l[0] / moments[0], l[1] / moments[1], l[2] / moments[2]};
  double cross[3];
  myb_vec3_cross(om, l, cross);
  if (myb_vec3_norm(cross) == 0.0) {
    f->kind = KIND_STEADY;
    for (int i = 0; i < 9; ++i) {
      f->perm[i] = (i % 4 == 0) ? 1.0 : 0.0;
      f->m0[i] = m0[i];
    }
    memcpy(f->omega, om, sizeof om);
    memcpy(f->inertia, moments, 3 * sizeof *moments);
    f->g = myb_vec3_norm(l);
    return MYB_OK;
  }

  /* order the axes: 2 = the middle moment, 3 = the axis the momentum circulates around */
  int lo = 0;
  int mid = 1;
  int hi = 2;
  if (moments[lo] > moments[mid]) {
    const int s = lo;
    lo = mid;
    mid = s;
  }
  if (moments[mid] > moments[hi]) {
    const int s = mid;
    mid = hi;
    hi = s;
  }
  if (moments[lo] > moments[mid]) {
    const int s = lo;
    lo = mid;
    mid = s;
  }
  int axis[3];
  if (moments[lo] == moments[mid]) {
    axis[0] = lo; /* symmetric: the distinct axis is 3 */
    axis[1] = mid;
    axis[2] = hi;
  } else if (moments[mid] == moments[hi]) {
    axis[0] = mid;
    axis[1] = hi;
    axis[2] = lo;
  } else {
    /* G^2 - 2 T I_mid without cancellation */
    const double i1 = moments[lo];
    const double i2 = moments[mid];
    const double i3 = moments[hi];
    const double d = l[lo] * l[lo] * (i1 - i2) / i1 + l[hi] * l[hi] * (i3 - i2) / i3;
    if (d > 0.0) {
      axis[0] = lo;
      axis[1] = mid;
      axis[2] = hi;
    } else if (d < 0.0) {
      axis[0] = hi;
      axis[1] = mid;
      axis[2] = lo;
    } else {
      return MYB_E_DOMAIN; /* separatrix */
    }
  }
  memset(f->perm, 0, sizeof f->perm);
  for (int r = 0; r < 3; ++r) {
    f->perm[3 * r + axis[r]] = 1.0;
  }
  /* right-handed: det of a permutation matrix is its sign */
  double c01[3];
  myb_vec3_cross(f->perm, f->perm + 3, c01);
  if (myb_vec3_dot(c01, f->perm + 6) < 0.0) {
    for (int x = 0; x < 3; ++x) {
      f->perm[3 + x] = -f->perm[3 + x];
    }
  }
  double mo[3];
  myb_internal_mat3_vec(f->perm, l, mo);
  for (int r = 0; r < 3; ++r) {
    f->inertia[r] = moments[axis[r]];
  }
  const double i1 = f->inertia[0];
  const double i2 = f->inertia[1];
  const double i3 = f->inertia[2];
  f->g = myb_vec3_norm(mo);

  /* amplitudes, frequency and parameter without cancellation (ratios of differences keep their
   * signs in either order of the moments) */
  const double a1 = sqrt(mo[0] * mo[0] + mo[1] * mo[1] * i1 * (i3 - i2) / (i2 * (i3 - i1)));
  const double a2 = sqrt(mo[0] * mo[0] * i2 * (i3 - i1) / (i1 * (i3 - i2)) + mo[1] * mo[1]);
  const double a3 = sqrt(mo[2] * mo[2] + mo[1] * mo[1] * i3 * (i2 - i1) / (i2 * (i3 - i1)));
  const double g2_2ti1 = mo[1] * mo[1] * (i2 - i1) / i2 + mo[2] * mo[2] * (i3 - i1) / i3;
  const double ti3_g2 = mo[0] * mo[0] * (i3 - i1) / i1 + mo[1] * mo[1] * (i3 - i2) / i2;
  f->nu = sqrt((i3 - i2) * g2_2ti1 / (i1 * i2 * i3));
  f->k2 = (i2 - i1) * ti3_g2 / ((i3 - i2) * g2_2ti1);
  f->n = -i3 * (i2 - i1) / (i1 * (i3 - i2));
  if (!(f->k2 >= 0.0 && f->k2 < 1.0) || !(f->nu > 0.0) || !isfinite(f->n)) {
    return MYB_E_DOMAIN;
  }
  f->a[0] = a1;
  f->a[1] = a2;
  f->a[2] = a3;
  f->sign3 = mo[2] < 0.0 ? -1.0 : 1.0;
  f->sign2 = (i3 > i2 ? 1.0 : -1.0) * f->sign3;
  /* the phase: sn u0 = L2 / (s a2), cn u0 = L1 / a1 */
  double am0 = 0.0;
  if (a1 > 0.0 && a2 > 0.0) {
    am0 = myb_atan2(f->sign2 * mo[1] / a2, mo[0] / a1);
  }
  f->u0 = myb_elliptic_f_inc(am0, f->k2);
  f->pi0 = myb_elliptic_pi_inc(f->n, am0, f->k2);
  f->rate0 = f->g / i3;
  f->rate_pi = f->g * (i3 - i1) / (i1 * i3 * f->nu);

  /* Z = (R3(psi0) R1(theta0))^T M0 in ordered axes, so that M(0) = M0 */
  double mo0[9];
  myb_internal_mat3_mul(f->perm, m0, mo0);
  static const double identity[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  memcpy(f->z, identity, sizeof identity);
  double rpt[9];
  general_attitude(f, mo, 0.0, rpt); /* R3(psi0) R1(theta0) */
  double rpt_t[9];
  myb_internal_mat3_transpose(rpt, rpt_t);
  myb_internal_mat3_mul(rpt_t, mo0, f->z);
  return MYB_OK;
}

void myb_free_rotation_state(const myb_free_rotation* rotation, double t, double q[4],
                             double l[3]) {
  MYB_CHECK(rotation != NULL);
  const myb_free_rotation* f = rotation;
  double m[9];
  double lb[3];
  if (f->kind == KIND_STEADY) {
    /* M(t) = R_omega(|omega| t) M0, about the body axis omega */
    const double w = myb_vec3_norm(f->omega);
    double r[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    if (w > 0.0) {
      const double phi[3] = {f->omega[0] * t, f->omega[1] * t, f->omega[2] * t};
      double qr[4];
      myb_quat_from_rotation_vector(phi, qr);
      myb_quat_to_matrix(qr, r);
    }
    myb_internal_mat3_mul(r, f->m0, m);
    for (int i = 0; i < 3; ++i) {
      lb[i] = f->omega[i] * f->inertia[i];
    }
  } else {
    double mo[3];
    double phi;
    general_motion(f, t, mo, &phi);
    double mt[9];
    general_attitude(f, mo, phi, mt);
    /* back to principal axes: M = P^T M_ordered, L = P^T L_ordered */
    double pt[9];
    myb_internal_mat3_transpose(f->perm, pt);
    myb_internal_mat3_mul(pt, mt, m);
    myb_internal_mat3_vec(pt, mo, lb);
  }
  if (q != NULL) {
    myb_quat_from_matrix(m, q);
  }
  if (l != NULL) {
    memcpy(l, lb, sizeof lb);
  }
}

/* ---- Splitting ------------------------------------------------------------------------------- */

myb_status myb_rotation_split(const double moments[3], myb_torque_function torque, void* context,
                              int order, double t0, double h, int64_t steps, double q[4],
                              double l[3]) {
  MYB_CHECK(moments != NULL && torque != NULL && q != NULL && l != NULL);
  const double* weights = NULL;
  int stages = 0;
  switch (order) {
  case 2:
    weights = myb_symplectic_order2;
    stages = 1;
    break;
  case 4:
    weights = myb_symplectic_order4;
    stages = 3;
    break;
  case 6:
    weights = myb_symplectic_order6;
    stages = 9;
    break;
  case 8:
    weights = myb_symplectic_order8;
    stages = 27;
    break;
  default:
    return MYB_E_DOMAIN;
  }
  if (steps < 0 || !isfinite(h) || !isfinite(t0)) {
    return MYB_E_DOMAIN;
  }
  double t = t0;
  double tc = 0.0; /* compensation of t */
  for (int64_t step = 0; step < steps; ++step) {
    for (int s = 0; s < stages; ++s) {
      const double hs = weights[s] * h;
      double tau[3];
      myb_status st = torque(context, t + tc, q, tau);
      if (st != MYB_OK) {
        return st;
      }
      for (int i = 0; i < 3; ++i) {
        l[i] += 0.5 * hs * tau[i];
      }
      myb_free_rotation f;
      st = myb_free_rotation_init(moments, q, l, &f);
      if (st != MYB_OK) {
        return st;
      }
      myb_free_rotation_state(&f, hs, q, l);
      double err;
      t = myb_two_sum(t, hs, &err);
      tc += err;
      st = torque(context, t + tc, q, tau);
      if (st != MYB_OK) {
        return st;
      }
      for (int i = 0; i < 3; ++i) {
        l[i] += 0.5 * hs * tau[i];
      }
    }
  }
  return MYB_OK;
}
