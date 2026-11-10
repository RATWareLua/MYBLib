/* SPDX-License-Identifier: BSD-3-Clause */
/* Attitude: quaternion conventions and identities, the twelve Euler sequences, principal axes,
 * torques (point mass, field, figure-figure with the balance of angular momentum), the exact
 * torque-free rotation against DOP853 at tolerance 1e-14 and against the symmetric top in
 * closed form, and the splitting integrator (exact without torque, orders 2-8 with a
 * gravity-gradient torque). */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/attitude.h"
#include "myb/frames.h"
#include "myb/grav.h"
#include "myb/integ.h"

static const double PI = 3.141592653589793;

/* uniform in [0, 1) */
static double uniform(mt_rng* rng) {
  return (double)(mt_next(rng) >> 11) * 0x1p-53;
}

static void random_quat(mt_rng* rng, double q[4]) {
  for (int i = 0; i < 4; ++i) {
    q[i] = 2.0 * uniform(rng) - 1.0;
  }
  myb_quat_normalize(q);
}

static double mat_diff(const double a[9], const double b[9]) {
  double d = 0.0;
  for (int i = 0; i < 9; ++i) {
    d = fmax(d, fabs(a[i] - b[i]));
  }
  return d;
}

static void mat_mul(const double a[9], const double b[9], double out[9]) {
  double w[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      w[3 * i + j] = a[3 * i] * b[j] + a[3 * i + 1] * b[3 + j] + a[3 * i + 2] * b[6 + j];
    }
  }
  memcpy(out, w, sizeof w);
}

/* rotation angle between two attitudes */
static double attitude_angle(const double q1[4], const double q2[4]) {
  double c[4];
  double p[4];
  myb_quat_conjugate(q1, c);
  myb_quat_multiply(c, q2, p);
  double phi[3];
  myb_quat_to_rotation_vector(p, phi);
  return sqrt(phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]);
}

static void test_quaternions(void) {
  mt_rng rng = {1};
  double worst_rt = 0.0;
  double worst_comp = 0.0;
  double worst_ortho = 0.0;
  double worst_log = 0.0;
  for (int k = 0; k < 2000; ++k) {
    double a[4];
    double b[4];
    random_quat(&rng, a);
    random_quat(&rng, b);
    double ma[9];
    double mb[9];
    myb_quat_to_matrix(a, ma);
    myb_quat_to_matrix(b, mb);
    double back[4];
    myb_quat_from_matrix(ma, back);
    const double s = a[0] < 0.0 ? -1.0 : 1.0;
    for (int i = 0; i < 4; ++i) {
      worst_rt = fmax(worst_rt, fabs(back[i] - s * a[i]));
    }
    double ab[4];
    myb_quat_multiply(a, b, ab);
    double mab[9];
    double mbma[9];
    myb_quat_to_matrix(ab, mab);
    mat_mul(mb, ma, mbma);
    worst_comp = fmax(worst_comp, mat_diff(mab, mbma)); /* M(a * b) = M(b) M(a) */
    double mmt[9];
    double mt[9];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        mt[3 * i + j] = ma[3 * j + i];
      }
    }
    mat_mul(ma, mt, mmt);
    const double id[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    worst_ortho = fmax(worst_ortho, mat_diff(mmt, id));
    double phi[3];
    myb_quat_to_rotation_vector(a, phi);
    double e[4];
    myb_quat_from_rotation_vector(phi, e);
    for (int i = 0; i < 4; ++i) {
      worst_log = fmax(worst_log, fabs(e[i] - s * a[i]));
    }
  }
  printf("quaternions: matrix round trip %.2e, composition %.2e, orthonormal %.2e, exp(log) %.2e\n",
         worst_rt, worst_comp, worst_ortho, worst_log);
  MT_CHECK(worst_rt < 5e-16 && worst_comp < 1e-15 && worst_ortho < 1e-15 && worst_log < 5e-16);

  /* the convention: q = (cos a/2, 0, 0, sin a/2) is the frame rotation R3(a) */
  const double a = 0.7;
  const double q[4] = {cos(0.35), 0.0, 0.0, sin(0.35)};
  double m[9];
  myb_quat_to_matrix(q, m);
  double r3[9];
  const int seq[3] = {3, 1, 3};
  const double ang[3] = {a, 0.0, 0.0};
  MT_CHECK(myb_euler_to_matrix(seq, ang, r3) == MYB_OK);
  MT_CHECK(mat_diff(m, r3) < 2e-16);
  MT_CHECK(fabs(m[1] - sin(a)) < 2e-16 && fabs(m[3] + sin(a)) < 2e-16);

  /* small rotation vectors keep their relative precision */
  const double tiny[3] = {1e-20, -2e-20, 3e-20};
  double qt[4];
  double back[3];
  myb_quat_from_rotation_vector(tiny, qt);
  myb_quat_to_rotation_vector(qt, back);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(back[i] - tiny[i]) <= 1e-36);
  }

  /* slerp: endpoints, and the midpoint of a rotation about z is half the angle */
  double qa[4] = {1, 0, 0, 0};
  double qb[4] = {cos(1.0), 0, 0, sin(1.0)};
  double out[4];
  myb_quat_slerp(qa, qb, 0.0, out);
  MT_CHECK(attitude_angle(out, qa) < 1e-15);
  myb_quat_slerp(qa, qb, 1.0, out);
  MT_CHECK(attitude_angle(out, qb) < 1e-15);
  myb_quat_slerp(qa, qb, 0.5, out);
  MT_CHECK(fabs(out[0] - cos(0.5)) < 2e-16 && fabs(out[3] - sin(0.5)) < 2e-16);

  /* kinematics: q(t) = q0 * exp(omega t / 2) has dq/dt = q * (0, omega) / 2 and M' = -[w x] M */
  double q0[4];
  random_quat(&rng, q0);
  const double omega[3] = {0.3, -0.2, 0.5};
  const double h = 1e-5;
  double qp[4];
  double qm[4];
  const double pp[3] = {omega[0] * h, omega[1] * h, omega[2] * h};
  const double pm[3] = {-omega[0] * h, -omega[1] * h, -omega[2] * h};
  double ep[4];
  double em[4];
  myb_quat_from_rotation_vector(pp, ep);
  myb_quat_from_rotation_vector(pm, em);
  myb_quat_multiply(q0, ep, qp);
  myb_quat_multiply(q0, em, qm);
  double dq[4];
  myb_quat_derivative(q0, omega, dq);
  double worst_dq = 0.0;
  for (int i = 0; i < 4; ++i) {
    worst_dq = fmax(worst_dq, fabs(dq[i] - (qp[i] - qm[i]) / (2.0 * h)));
  }
  MT_CHECK(worst_dq < 1e-10);
  double m0[9];
  double mp[9];
  double mm[9];
  myb_quat_to_matrix(q0, m0);
  myb_quat_to_matrix(qp, mp);
  myb_quat_to_matrix(qm, mm);
  /* -[w x] M0 */
  const double wx[9] = {0, -omega[2], omega[1], omega[2], 0, -omega[0], -omega[1], omega[0], 0};
  double wm[9];
  mat_mul(wx, m0, wm);
  double worst_dm = 0.0;
  for (int i = 0; i < 9; ++i) {
    worst_dm = fmax(worst_dm, fabs(-wm[i] - (mp[i] - mm[i]) / (2.0 * h)));
  }
  MT_CHECK(worst_dm < 1e-10);
}

static void test_euler(void) {
  mt_rng rng = {2};
  double worst_matrix = 0.0;
  double worst_angles = 0.0;
  double worst_lock = 0.0;
  int sequences = 0;
  for (int i = 1; i <= 3; ++i) {
    for (int j = 1; j <= 3; ++j) {
      for (int k = 1; k <= 3; ++k) {
        const int axes[3] = {i, j, k};
        if (i == j || j == k) {
          double m[9];
          const double z[3] = {0, 0, 0};
          MT_CHECK(myb_euler_to_matrix(axes, z, m) == MYB_E_DOMAIN);
          continue;
        }
        ++sequences;
        const int proper = i == k;
        for (int s = 0; s < 300; ++s) {
          double ang[3];
          ang[0] = (uniform(&rng) * 2.0 - 1.0) * 3.0;
          ang[2] = (uniform(&rng) * 2.0 - 1.0) * 3.0;
          ang[1] = proper ? 0.05 + 3.0 * uniform(&rng) : 1.5 * (2.0 * uniform(&rng) - 1.0);
          double m[9];
          double back[3];
          double m2[9];
          MT_CHECK(myb_euler_to_matrix(axes, ang, m) == MYB_OK);
          MT_CHECK(myb_euler_from_matrix(axes, m, back) == MYB_OK);
          MT_CHECK(myb_euler_to_matrix(axes, back, m2) == MYB_OK);
          worst_matrix = fmax(worst_matrix, mat_diff(m, m2));
          for (int x = 0; x < 3; ++x) {
            worst_angles = fmax(worst_angles, fabs(back[x] - ang[x]));
          }
        }
        /* gimbal lock: b = 0 or pi (proper), +-pi/2 (Tait-Bryan) */
        const double locks[2] = {proper ? 0.0 : PI / 2, proper ? PI : -PI / 2};
        for (int x = 0; x < 2; ++x) {
          const double ang[3] = {0.4, locks[x], -1.1};
          double m[9];
          double back[3];
          double m2[9];
          MT_CHECK(myb_euler_to_matrix(axes, ang, m) == MYB_OK);
          MT_CHECK(myb_euler_from_matrix(axes, m, back) == MYB_OK);
          MT_CHECK(myb_euler_to_matrix(axes, back, m2) == MYB_OK);
          worst_lock = fmax(worst_lock, mat_diff(m, m2));
        }
      }
    }
  }
  printf("Euler: %d sequences, matrix round trip %.2e, angles %.2e, at gimbal lock %.2e\n",
         sequences, worst_matrix, worst_angles, worst_lock);
  MT_CHECK(sequences == 12);
  MT_CHECK(worst_matrix < 1e-15 && worst_angles < 1e-13 && worst_lock < 1e-15);

  /* 3-1-3 agrees with the WGCCRE helper of myb/frames.h */
  const double ang[3] = {1.1, 0.4, -2.3};
  myb_rotation r;
  myb_rotation_euler313(ang, ang, &r);
  double m[9];
  const int axes[3] = {3, 1, 3};
  MT_CHECK(myb_euler_to_matrix(axes, ang, m) == MYB_OK);
  MT_CHECK(mat_diff(m, r.m) < 1e-16);
}

static void test_principal_axes(void) {
  mt_rng rng = {3};
  double worst = 0.0;
  double worst_m = 0.0;
  for (int k = 0; k < 500; ++k) {
    double q[4];
    random_quat(&rng, q);
    double r[9];
    myb_quat_to_matrix(q, r);
    double d[3] = {1.0 + uniform(&rng), 2.0 + uniform(&rng),
                   3.0 + uniform(&rng)};
    if (k % 50 == 0) {
      d[1] = d[0]; /* a symmetric body */
    }
    /* I = R^T D R */
    double in[9];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        in[3 * i + j] = r[i] * d[0] * r[j] + r[3 + i] * d[1] * r[3 + j] + r[6 + i] * d[2] * r[6 + j];
      }
    }
    in[3] = in[1];
    in[6] = in[2];
    in[7] = in[5];
    double moments[3];
    double axes[9];
    MT_CHECK(myb_principal_axes(in, moments, axes) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      worst_m = fmax(worst_m, fabs(moments[i] - d[i]) / d[i]);
    }
    /* axes I axes^T = diag(moments) */
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        double s = 0.0;
        for (int x = 0; x < 3; ++x) {
          for (int y = 0; y < 3; ++y) {
            s += axes[3 * i + x] * in[3 * x + y] * axes[3 * j + y];
          }
        }
        worst = fmax(worst, fabs(s - (i == j ? moments[i] : 0.0)) / 4.0);
      }
    }
    double c[3] = {axes[1] * axes[5] - axes[2] * axes[4], axes[2] * axes[3] - axes[0] * axes[5],
                   axes[0] * axes[4] - axes[1] * axes[3]};
    MT_CHECK(c[0] * axes[6] + c[1] * axes[7] + c[2] * axes[8] > 0.999);
  }
  printf("principal axes: diagonalisation %.2e, moments %.2e\n", worst, worst_m);
  MT_CHECK(worst < 2e-15 && worst_m < 2e-15);
  const double bad[9] = {1, 0.1, 0, 0.2, 1, 0, 0, 0, 1};
  double mo[3];
  double ax[9];
  MT_CHECK(myb_principal_axes(bad, mo, ax) == MYB_E_DOMAIN);
}

/* ---- Torques -------------------------------------------------------------------------------- */

/* U22 as in myb/grav.h with G Q, without the factor 1/G */
static double u22(const double ga[9], const double gb[9], const double r[3]) {
  const double rr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  const double n[3] = {r[0] / rr, r[1] / rr, r[2] / rr};
  double an[3] = {0, 0, 0};
  double bn[3] = {0, 0, 0};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      an[i] += ga[3 * i + j] * n[j];
      bn[i] += gb[3 * i + j] * n[j];
    }
  }
  double nan = 0.0, nbn = 0.0, nabn = 0.0, tr = 0.0;
  for (int i = 0; i < 3; ++i) {
    nan += n[i] * an[i];
    nbn += n[i] * bn[i];
    nabn += an[i] * bn[i];
    for (int j = 0; j < 3; ++j) {
      tr += ga[3 * i + j] * gb[3 * j + i];
    }
  }
  return (105.0 * nan * nbn - 60.0 * nabn + 6.0 * tr) / (36.0 * pow(rr, 5));
}

static void rotate_tensor(const double q[4], const double t[9], double out[9]) {
  double m[9];
  myb_quat_to_matrix(q, m);
  /* active rotation of the body: R = M^T; R t R^T */
  double w[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double s = 0.0;
      for (int x = 0; x < 3; ++x) {
        for (int y = 0; y < 3; ++y) {
          s += m[3 * x + i] * t[3 * x + y] * m[3 * y + j];
        }
      }
      w[3 * i + j] = s;
    }
  }
  memcpy(out, w, sizeof w);
}

static void test_torques(void) {
  /* point mass: 3 gm / r^5 (r x I r) = -(m_p / m_b) r x grad U2 with the MacCullagh potential
   * U2 = GM (tr i - 3 r.i.r / r^2) / (2 r^3) of the body (inertia per unit mass i) */
  const double gm_body = 4.9e12;
  const double gm_p = 3.986e14;
  const double inertia[3] = {1.0e12, 1.05e12, 1.2e12};
  const double r[3] = {2.1e8, -1.3e8, 0.9e8};
  double tau_pm[3];
  myb_torque_point_mass(gm_p, r, inertia, tau_pm);
  double grad[3];
  for (int c = 0; c < 3; ++c) {
    double rp[3] = {r[0], r[1], r[2]};
    double rm[3] = {r[0], r[1], r[2]};
    const double h = 1e3;
    rp[c] += h;
    rm[c] -= h;
    double up = 0.0;
    double um = 0.0;
    for (int s = 0; s < 2; ++s) {
      const double* x = s == 0 ? rp : rm;
      const double x2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
      const double xix = inertia[0] * x[0] * x[0] + inertia[1] * x[1] * x[1] + inertia[2] * x[2] * x[2];
      const double u = gm_body * ((inertia[0] + inertia[1] + inertia[2]) - 3.0 * xix / x2) /
                       (2.0 * x2 * sqrt(x2));
      if (s == 0) {
        up = u;
      } else {
        um = u;
      }
    }
    grad[c] = (up - um) / (2e3);
  }
  double tau_f[3];
  myb_torque_field(gm_p / gm_body, r, grad, tau_f);
  double worst = 0.0;
  const double scale = sqrt(tau_pm[0] * tau_pm[0] + tau_pm[1] * tau_pm[1] + tau_pm[2] * tau_pm[2]);
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(tau_f[i] - tau_pm[i]) / scale);
  }
  printf("torques: point mass vs field gradient %.2e\n", worst);
  MT_CHECK(worst < 1e-7); /* the central difference */

  /* figure-figure: against the derivative of U22 under a rotation of A, and the balance of
   * angular momentum gm_a tau_a + gm_b tau_b + gm_b r x a_b = 0 */
  const double ga[9] = {3e15, 2e14, -1e14, 2e14, -1e15, 5e13, -1e14, 5e13, -2e15};
  const double gb[9] = {-4e19, 1e18, 3e18, 1e18, 2e19, -2e18, 3e18, -2e18, 2e19};
  const double gm_a = 4.9e12;
  const double gm_b = 3.986e14;
  const double rab[3] = {3.1e8, 2.2e8, -0.4e8};
  double ta[3];
  double tb[3];
  MT_CHECK(myb_torque_figure_figure(gm_a, ga, gb, rab, ta) == MYB_OK);
  const double rba[3] = {-rab[0], -rab[1], -rab[2]};
  MT_CHECK(myb_torque_figure_figure(gm_b, gb, ga, rba, tb) == MYB_OK);
  double worst_fd = 0.0;
  const double tnorm = sqrt(ta[0] * ta[0] + ta[1] * ta[1] + ta[2] * ta[2]);
  for (int c = 0; c < 3; ++c) {
    const double h = 1e-6;
    double phi[3] = {0, 0, 0};
    phi[c] = h;
    double qp[4];
    myb_quat_from_rotation_vector(phi, qp);
    phi[c] = -h;
    double qm[4];
    myb_quat_from_rotation_vector(phi, qm);
    double gap[9];
    double gam[9];
    rotate_tensor(qp, ga, gap);
    rotate_tensor(qm, ga, gam);
    const double d = (u22(gap, gb, rab) - u22(gam, gb, rab)) / (2.0 * h) / gm_a;
    worst_fd = fmax(worst_fd, fabs(d - ta[c]) / tnorm);
  }
  double aa[3] = {0, 0, 0};
  double ab[3] = {0, 0, 0};
  MT_CHECK(myb_figure_figure(gm_a, ga, gm_b, gb, rab, aa, ab) == MYB_OK);
  const double rxa[3] = {rab[1] * ab[2] - rab[2] * ab[1], rab[2] * ab[0] - rab[0] * ab[2],
                         rab[0] * ab[1] - rab[1] * ab[0]};
  double worst_l = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double total = gm_a * ta[i] + gm_b * tb[i] + gm_b * rxa[i];
    worst_l = fmax(worst_l, fabs(total) / (gm_a * tnorm));
  }
  printf("figure-figure torque: vs rotation derivative %.2e, angular momentum balance %.2e\n",
         worst_fd, worst_l);
  MT_CHECK(worst_fd < 1e-8);
  MT_CHECK(worst_l < 1e-14);
}

/* ---- Torque-free rotation ------------------------------------------------------------------ */

typedef struct body_context {
  double moments[3];
  double gm;
  double r_inertial[3];
} body_context;

static myb_status rigid_rhs(void* context, double t, const double* y, double* dy) {
  (void)t;
  const body_context* b = context;
  double tau[3] = {0, 0, 0};
  if (b->gm != 0.0) {
    double rb[3];
    myb_quat_rotate(y, b->r_inertial, rb);
    myb_torque_point_mass(b->gm, rb, b->moments, tau);
  }
  myb_rigid_body_rates(b->moments, y, y + 4, tau, dy, dy + 4);
  return MYB_OK;
}

static myb_status torque_fn(void* context, double t, const double q[4], double torque[3]) {
  (void)t;
  const body_context* b = context;
  double rb[3];
  myb_quat_rotate(q, b->r_inertial, rb);
  myb_torque_point_mass(b->gm, rb, b->moments, torque);
  return MYB_OK;
}

static myb_status zero_torque(void* context, double t, const double q[4], double torque[3]) {
  (void)context;
  (void)t;
  (void)q;
  torque[0] = torque[1] = torque[2] = 0.0;
  return MYB_OK;
}

static void reference(body_context* b, const double q0[4], const double l0[3], double t,
                      double q[4], double l[3]) {
  const myb_dop853_options o = {.rtol = 1e-14, .atol = 1e-16};
  myb_dop853* d = NULL;
  MT_CHECK(myb_dop853_create(7, &o, NULL, &d) == MYB_OK);
  double y0[7] = {q0[0], q0[1], q0[2], q0[3], l0[0], l0[1], l0[2]};
  MT_CHECK(myb_dop853_start(d, rigid_rhs, b, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_dop853_integrate(d, t) == MYB_OK);
  const double* y = myb_dop853_state(d);
  memcpy(q, y, 4 * sizeof *q);
  memcpy(l, y + 4, 3 * sizeof *l);
  myb_quat_normalize(q);
  myb_dop853_destroy(d);
}

static void test_free_rotation(void) {
  mt_rng rng = {4};
  /* cases: around the largest axis, around the smallest, symmetric (both orders), steady */
  const double moments[6][3] = {{1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}, {3.0, 1.0, 2.0},
                                {2.0, 2.0, 3.0}, {3.0, 1.5, 3.0}, {1.0, 2.0, 3.0}};
  const double ls[6][3] = {{0.3, 0.2, 1.0}, {1.0, 0.3, 0.2}, {0.4, -0.3, 0.9},
                           {0.6, -0.4, 0.7}, {-0.5, 0.8, 0.3}, {0.0, 0.0, 1.0}};
  double worst_att = 0.0;
  double worst_l = 0.0;
  double worst_inv = 0.0;
  for (int c = 0; c < 6; ++c) {
    body_context b = {.moments = {moments[c][0], moments[c][1], moments[c][2]}, .gm = 0.0};
    double q0[4];
    random_quat(&rng, q0);
    myb_free_rotation f;
    MT_CHECK(myb_free_rotation_init(b.moments, q0, ls[c], &f) == MYB_OK);
    for (int k = 1; k <= 4; ++k) {
      const double t = 7.3 * k; /* several periods */
      double q[4];
      double l[3];
      myb_free_rotation_state(&f, t, q, l);
      double qr[4];
      double lr[3];
      reference(&b, q0, ls[c], t, qr, lr);
      worst_att = fmax(worst_att, attitude_angle(q, qr));
      for (int i = 0; i < 3; ++i) {
        worst_l = fmax(worst_l, fabs(l[i] - lr[i]));
      }
      /* the inertial momentum M^T L is constant */
      double m[9];
      double m0[9];
      myb_quat_to_matrix(q, m);
      myb_quat_to_matrix(q0, m0);
      for (int i = 0; i < 3; ++i) {
        const double li = m[i] * l[0] + m[3 + i] * l[1] + m[6 + i] * l[2];
        const double l0i = m0[i] * ls[c][0] + m0[3 + i] * ls[c][1] + m0[6 + i] * ls[c][2];
        worst_inv = fmax(worst_inv, fabs(li - l0i));
      }
    }
    double q[4];
    double l[3];
    myb_free_rotation_state(&f, 0.0, q, l);
    MT_CHECK(attitude_angle(q, q0) < 1e-15);
  }
  printf("free rotation vs DOP853 (1e-14): attitude %.2e rad, momentum %.2e; inertial momentum "
         "%.2e\n",
         worst_att, worst_l, worst_inv);
  MT_CHECK(worst_att < 1e-13); /* the reference's own tolerance */
  MT_CHECK(worst_l < 1e-13);
  MT_CHECK(worst_inv < 2e-14);

  /* symmetric top in closed form: M(t) = R3(lambda t) M0 R_l(G t / A), lambda = (1/C - 1/A) L3 */
  const double ms[3] = {2.0, 2.0, 5.0};
  const double l0[3] = {0.7, -0.2, 1.3};
  double q0[4];
  random_quat(&rng, q0);
  myb_free_rotation f;
  MT_CHECK(myb_free_rotation_init(ms, q0, l0, &f) == MYB_OK);
  double m0[9];
  myb_quat_to_matrix(q0, m0);
  double ln[3];
  for (int i = 0; i < 3; ++i) {
    ln[i] = m0[i] * l0[0] + m0[3 + i] * l0[1] + m0[6 + i] * l0[2];
  }
  const double g = sqrt(l0[0] * l0[0] + l0[1] * l0[1] + l0[2] * l0[2]);
  double worst_sym = 0.0;
  for (int k = 1; k <= 20; ++k) {
    const double t = 3.7 * k;
    const double lambda = (1.0 / ms[2] - 1.0 / ms[0]) * l0[2];
    const double pb[3] = {0.0, 0.0, lambda * t};
    const double pl[3] = {ln[0] / g * g / ms[0] * t, ln[1] / g * g / ms[0] * t, ln[2] / g * g / ms[0] * t};
    double qb[4];
    double ql[4];
    myb_quat_from_rotation_vector(pb, qb);
    myb_quat_from_rotation_vector(pl, ql);
    /* M = M(qb) M0 M(ql) = M(ql * q0 * qb) */
    double w[4];
    double qx[4];
    myb_quat_multiply(ql, q0, w);
    myb_quat_multiply(w, qb, qx);
    double q[4];
    myb_free_rotation_state(&f, t, q, NULL);
    worst_sym = fmax(worst_sym, attitude_angle(q, qx));
  }
  printf("symmetric top vs closed form: %.2e rad\n", worst_sym);
  MT_CHECK(worst_sym < 1e-13);

  /* the separatrix |L|^2 = 2 T I2, exactly: moments (1, 3, 6), L = (1, 0.7, 2) */
  const double sep_m[3] = {1.0, 3.0, 6.0};
  const double sep_l[3] = {1.0, 0.7, 2.0};
  MT_CHECK(myb_free_rotation_init(sep_m, q0, sep_l, &f) == MYB_E_DOMAIN);
  const double bad[3] = {1.0, -2.0, 3.0};
  MT_CHECK(myb_free_rotation_init(bad, q0, l0, &f) == MYB_E_DOMAIN);
}

static void test_splitting(void) {
  mt_rng rng = {5};
  body_context b = {.moments = {1.0, 1.6, 2.1}, .gm = 0.0, .r_inertial = {2.0, 1.0, -0.5}};
  double q0[4];
  random_quat(&rng, q0);
  const double l0[3] = {0.3, 0.5, 1.1};
  /* no torque: exact */
  double q[4];
  double l[3];
  memcpy(q, q0, sizeof q);
  memcpy(l, l0, sizeof l);
  MT_CHECK(myb_rotation_split(b.moments, zero_torque, &b, 4, 0.0, 0.37, 50, q, l) == MYB_OK);
  myb_free_rotation f;
  MT_CHECK(myb_free_rotation_init(b.moments, q0, l0, &f) == MYB_OK);
  double qe[4];
  double le[3];
  myb_free_rotation_state(&f, 0.37 * 50, qe, le);
  const double free_err = attitude_angle(q, qe);
  printf("splitting without torque vs exact: %.2e rad\n", free_err);
  MT_CHECK(free_err < 1e-12);

  /* gravity gradient of a point mass fixed in space: orders */
  b.gm = 0.8;
  const double t_end = 6.0;
  double qr[4];
  double lr[3];
  reference(&b, q0, l0, t_end, qr, lr);
  const int orders[4] = {2, 4, 6, 8};
  for (int o = 0; o < 4; ++o) {
    double errors[2];
    for (int s = 0; s < 2; ++s) {
      const int64_t steps = (s == 0 ? 40 : 80) * (o == 0 ? 4 : 1);
      memcpy(q, q0, sizeof q);
      memcpy(l, l0, sizeof l);
      MT_CHECK(myb_rotation_split(b.moments, torque_fn, &b, orders[o], 0.0, t_end / (double)steps,
                                  steps, q, l) == MYB_OK);
      errors[s] = attitude_angle(q, qr);
    }
    const double observed = log2(errors[0] / errors[1]);
    printf("splitting order %d: error %.2e -> %.2e, observed order %.2f\n", orders[o], errors[0],
           errors[1], observed);
    if (errors[1] > 1e-11) {
      MT_CHECK(observed > orders[o] - 0.4);
    }
  }
  MT_CHECK(myb_rotation_split(b.moments, torque_fn, &b, 3, 0.0, 0.1, 1, q, l) == MYB_E_DOMAIN);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_quaternions();
  test_euler();
  test_principal_axes();
  test_torques();
  test_free_rotation();
  test_splitting();
  return mt_finish("test_attitude");
}
