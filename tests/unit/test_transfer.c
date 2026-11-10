/* SPDX-License-Identifier: BSD-3-Clause */
/* Transfers and effect budgets:
 *  - the B-plane of a hyperbola built from its periapsis: |B| = r_p sqrt(1 + 2 gm / (r_p v^2)),
 *    B perpendicular to S, the time to periapsis, B.T and B.R of a known geometry;
 *  - patched-conic numbers (sphere of influence and Hill radius of the Earth, departure
 *    impulses, turn angles);
 *  - shooting: the Lambert problem solved by Newton with the closed-form state transition
 *    matrix against myb_lambert; a minimum-norm step; targeting the B-plane of an Earth flyby
 *    perturbed by the Moon, with the state transition matrix of the variational equations;
 *  - the budget of a low orbit: J2 and a missing drag flagged, a negligible term flagged,
 *    effects linear in a small force. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/budget.h"
#include "myb/kepler.h"
#include "myb/lambert.h"
#include "myb/model.h"
#include "myb/transfer.h"
#include "myb/uncert.h"

static const double gm_earth = 3.986004418e14;

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static void test_bplane(void) {
  const double rp = 7.0e6;
  const double vinf = 5000.0;
  const double vp = sqrt(vinf * vinf + 2.0 * gm_earth / rp);
  /* periapsis on the x axis, motion along y in a plane tilted about x by 30 degrees */
  const double inc = 30.0 * 3.141592653589793 / 180.0;
  const myb_pv peri = {{rp, 0.0, 0.0}, {0.0, vp * cos(inc), vp * sin(inc)}};
  myb_pv before;
  MT_CHECK(myb_kepler_propagate(gm_earth, &peri, -2.0 * 86400.0, &before) == MYB_OK);
  const double k[3] = {0.0, 0.0, 1.0};
  myb_bplane b;
  MT_CHECK(myb_bplane_compute(gm_earth, before.r, before.v, k, &b) == MYB_OK);
  const double bmag = rp * sqrt(1.0 + 2.0 * gm_earth / (rp * vinf * vinf));
  const double e = 1.0 + rp * vinf * vinf / gm_earth;
  MT_CHECK(fabs(b.b_magnitude - bmag) < 1e-9 * bmag);
  MT_CHECK(fabs(b.periapsis - rp) < 1e-8 * rp && fabs(b.eccentricity - e) < 1e-12 * e);
  MT_CHECK(fabs(b.v_infinity - vinf) < 1e-12 * vinf && fabs(b.c3 / (vinf * vinf) - 1.0) < 1e-12);
  MT_CHECK(fabs(b.time_to_periapsis - 2.0 * 86400.0) < 1e-6);
  MT_CHECK(fabs(b.b[0] * b.s[0] + b.b[1] * b.s[1] + b.b[2] * b.s[2]) < 1e-9 * bmag);
  MT_CHECK(fabs(hypot(b.b_t, b.b_r) - bmag) < 1e-9 * bmag);
  /* geometry, independently: far in the past the motion is along the asymptote; S is the
   * direction of the velocity there and B the part of the position perpendicular to S
   * (both converge as 1/r: 4e-7 at 1e5 days) */
  myb_pv far;
  MT_CHECK(myb_kepler_propagate(gm_earth, &peri, -1e5 * 86400.0, &far) == MYB_OK);
  const double vn = norm3(far.v);
  double s_far[3];
  for (int i = 0; i < 3; ++i) {
    s_far[i] = far.v[i] / vn;
  }
  const double rs = far.r[0] * s_far[0] + far.r[1] * s_far[1] + far.r[2] * s_far[2];
  double b_far[3];
  for (int i = 0; i < 3; ++i) {
    b_far[i] = far.r[i] - rs * s_far[i];
    MT_CHECK(fabs(s_far[i] - b.s[i]) < 1e-6);
  }
  const double bt_far = b_far[0] * b.t[0] + b_far[1] * b.t[1] + b_far[2] * b.t[2];
  const double br_far = b_far[0] * b.r[0] + b_far[1] * b.r[1] + b_far[2] * b.r[2];
  MT_CHECK(fabs(bt_far - b.b_t) < 1e-5 * bmag && fabs(br_far - b.b_r) < 1e-5 * bmag);
  /* after periapsis the time is negative; the same B-plane (same asymptote) */
  myb_pv after;
  MT_CHECK(myb_kepler_propagate(gm_earth, &peri, 3600.0, &after) == MYB_OK);
  myb_bplane b2;
  MT_CHECK(myb_bplane_compute(gm_earth, after.r, after.v, k, &b2) == MYB_OK);
  MT_CHECK(fabs(b2.time_to_periapsis + 3600.0) < 1e-6);
  MT_CHECK(fabs(b2.b_t - b.b_t) < 1e-6 && fabs(b2.b_r - b.b_r) < 1e-6);
  /* an ellipse is rejected */
  const double r[3] = {7e6, 0.0, 0.0};
  const double v[3] = {0.0, 7.5e3, 0.0};
  MT_CHECK(myb_bplane_compute(gm_earth, r, v, k, &b) == MYB_E_DOMAIN);
  printf("B-plane: |B| %.6f km for r_p 7000 km, v_inf 5 km/s; theta %.6f deg\n", bmag / 1e3,
         b2.theta * 180.0 / 3.141592653589793);
}

static void test_patched(void) {
  const double au = 1.495978707e11;
  const double earth = 3.003489596e-6; /* (Earth + Moon) / Sun */
  const double soi = myb_sphere_of_influence(au, earth);
  const double hill = myb_hill_radius(au, 0.0167, earth);
  printf("Earth: sphere of influence %.0f km, Hill radius %.0f km\n", soi / 1e3, hill / 1e3);
  MT_CHECK(fabs(soi / 1e3 - 925000.0) < 2000.0);
  MT_CHECK(fabs(hill / 1e3 - 1471000.0) < 5000.0);
  /* departure from a 300 km circular orbit with C3 = 9 km^2/s^2 */
  const double dv = myb_departure_impulse(gm_earth, 6678.137e3, 3000.0);
  MT_CHECK(fabs(dv - (sqrt(9e6 + 2.0 * gm_earth / 6678.137e3) - sqrt(gm_earth / 6678.137e3))) <
           1e-9);
  MT_CHECK(dv > 3500.0 && dv < 3700.0);
  /* turn angle = 2 asin(1/e) */
  const double turn = myb_flyby_turn(gm_earth, 5000.0, 7e6);
  const double e = 1.0 + 7e6 * 25e6 / gm_earth;
  MT_CHECK(fabs(turn - 2.0 * asin(1.0 / e)) < 1e-15);
}

/* ---- Shooting ------------------------------------------------------------------------------ */

typedef struct lambert_shot {
  double r1[3];
  double r2[3];
  double dt;
} lambert_shot;

static myb_status lambert_residual(void* context, const double* u, double* f, double* jac) {
  const lambert_shot* s = context;
  const myb_pv pv0 = {{s->r1[0], s->r1[1], s->r1[2]}, {u[0], u[1], u[2]}};
  myb_pv pv;
  double phi[36];
  const myb_status st = myb_kepler_stm(gm_earth, &pv0, s->dt, &pv, phi);
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 3; ++i) {
    f[i] = pv.r[i] - s->r2[i];
  }
  if (jac != NULL) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        jac[3 * i + j] = phi[i * 6 + 3 + j];
      }
    }
  }
  return MYB_OK;
}

typedef struct min_norm_problem {
  double target;
} min_norm_problem;

static myb_status plane_residual(void* context, const double* u, double* f, double* jac) {
  const min_norm_problem* p = context;
  f[0] = u[0] + 2.0 * u[1] - 3.0 * u[2] - p->target;
  if (jac != NULL) {
    jac[0] = 1.0;
    jac[1] = 2.0;
    jac[2] = -3.0;
  }
  return MYB_OK;
}

/* B-plane targeting: Earth plus the Moon (fixed at its position for the few days). */
typedef struct earth_term {
  double gm;
} earth_term;

static myb_status earth_accel(void* context, double t, const double x[3], const double v[3],
                              double a[3]) {
  (void)t;
  (void)v;
  const double r = norm3(x);
  const double k = -((const earth_term*)context)->gm / (r * r * r);
  for (int i = 0; i < 3; ++i) {
    a[i] = k * x[i];
  }
  return MYB_OK;
}

static myb_status earth_jac(void* context, double t, const double x[3], const double v[3],
                            double dadx[9], double dadv[9]) {
  (void)t;
  (void)v;
  myb_point_mass_jacobian(((const earth_term*)context)->gm, x, dadx);
  memset(dadv, 0, 9 * sizeof(double));
  return MYB_OK;
}

static const double moon_pos[3] = {3.0e8, -2.4e8, 0.5e8};

static myb_status moon_accel(void* context, double t, const double x[3], const double v[3],
                             double a[3]) {
  (void)context;
  (void)t;
  (void)v;
  const double gm = 4.9028e12;
  const double d[3] = {moon_pos[0] - x[0], moon_pos[1] - x[1], moon_pos[2] - x[2]};
  const double dd = norm3(d);
  const double s = norm3(moon_pos);
  for (int i = 0; i < 3; ++i) {
    a[i] = gm * (d[i] / (dd * dd * dd) - moon_pos[i] / (s * s * s));
  }
  return MYB_OK;
}

typedef struct bplane_shot {
  myb_force_model* model;
  double y0[6];
  double duration;
  double target[2]; /* B.T, B.R */
} bplane_shot;

static myb_status bplane_of(const double y[6], double* bt, double* br) {
  const double k[3] = {0.0, 0.0, 1.0};
  myb_bplane b;
  const myb_status st = myb_bplane_compute(gm_earth, y, y + 3, k, &b);
  *bt = b.b_t;
  *br = b.b_r;
  return st;
}

static myb_status bplane_residual(void* context, const double* u, double* f, double* jac) {
  bplane_shot* s = context;
  double y0[6];
  memcpy(y0, s->y0, sizeof y0);
  for (int i = 0; i < 3; ++i) {
    y0[3 + i] += u[i];
  }
  myb_variational* var = NULL;
  myb_status st = myb_variational_create(s->model, NULL, NULL, &var);
  if (st == MYB_OK) {
    st = myb_variational_start(var, 0.0, y0);
  }
  if (st == MYB_OK) {
    st = myb_variational_integrate(var, s->duration);
  }
  double y[6];
  double phi[36];
  if (st == MYB_OK) {
    myb_variational_get(var, y, phi, NULL);
  }
  myb_variational_destroy(var);
  if (st != MYB_OK) {
    return st;
  }
  double bt = 0.0;
  double br = 0.0;
  st = bplane_of(y, &bt, &br);
  if (st != MYB_OK) {
    return st;
  }
  f[0] = bt - s->target[0];
  f[1] = br - s->target[1];
  if (jac != NULL) {
    /* d(B.T, B.R)/d y by central differences of the two-body B-plane, times Phi[:, v0] */
    double db[2][6];
    for (int j = 0; j < 6; ++j) {
      double yp[6];
      double ym[6];
      memcpy(yp, y, sizeof yp);
      memcpy(ym, y, sizeof ym);
      const double h = j < 3 ? 1.0 : 1e-3;
      yp[j] += h;
      ym[j] -= h;
      double tp = 0.0;
      double rp = 0.0;
      double tm = 0.0;
      double rm = 0.0;
      if (bplane_of(yp, &tp, &rp) != MYB_OK || bplane_of(ym, &tm, &rm) != MYB_OK) {
        return MYB_E_DOMAIN;
      }
      db[0][j] = (tp - tm) / (2.0 * h);
      db[1][j] = (rp - rm) / (2.0 * h);
    }
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 3; ++j) {
        double sum = 0.0;
        for (int k = 0; k < 6; ++k) {
          sum += db[i][k] * phi[k * 6 + 3 + j];
        }
        jac[3 * i + j] = sum;
      }
    }
  }
  return MYB_OK;
}

static void test_shooting(void) {
  double work[64];
  MT_CHECK(myb_shoot_work_size(3, 3) <= 64 && myb_shoot_work_size(3, 2) <= 64);
  /* Lambert by shooting */
  lambert_shot ls = {{7.0e6, 1.0e6, -2.0e5}, {-3.0e6, 9.0e6, 2.0e6}, 5000.0};
  myb_lambert_solution sol[4];
  size_t count = 0;
  MT_CHECK(myb_lambert(gm_earth, ls.r1, ls.r2, ls.dt, 0, 0, sol, 4, &count) == MYB_OK ||
           count > 0);
  double u[3];
  const double r1n = norm3(ls.r1);
  for (int i = 0; i < 3; ++i) {
    u[i] = ls.r2[i] - ls.r1[i];
  }
  const double un = norm3(u);
  for (int i = 0; i < 3; ++i) {
    u[i] *= sqrt(gm_earth / r1n) / un; /* a crude first guess */
  }
  const myb_shoot_options o = {.tolerance = 1e-6};
  int its = 0;
  double res = 0.0;
  MT_CHECK(myb_shoot(lambert_residual, &ls, 3, 3, u, &o, work, &its, &res) == MYB_OK);
  double dv = 0.0;
  if (count > 0) {
    for (int i = 0; i < 3; ++i) {
      dv = fmax(dv, fabs(u[i] - sol[0].v1[i]));
    }
  }
  printf("Lambert by shooting with the closed-form STM: %d iterations, |v - lambert| %.1e m/s\n",
         its, dv);
  MT_CHECK(count > 0 && dv < 1e-8);
  /* minimum-norm step: u = J^T (J J^T)^-1 target from 0 */
  min_norm_problem mp = {14.0};
  double w[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_shoot(plane_residual, &mp, 3, 1, w, &o, work, &its, &res) == MYB_OK);
  MT_CHECK(fabs(w[0] - 1.0) < 1e-14 && fabs(w[1] - 2.0) < 1e-14 && fabs(w[2] + 3.0) < 1e-14);
  /* B-plane targeting of an Earth flyby with the Moon */
  earth_term et = {gm_earth};
  myb_force_model* m = NULL;
  MT_CHECK(myb_force_model_create(NULL, &m) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "earth", earth_accel, earth_jac, &et, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "moon", moon_accel, NULL, NULL, NULL) == MYB_OK);
  bplane_shot bs = {.model = m,
                    .y0 = {-6.0e8, 5.0e7, 1.0e7, 4500.0, 300.0, -50.0},
                    .duration = 1.3e5,
                    .target = {8000e3, -3000e3}};
  double dvel[3] = {0.0, 0.0, 0.0};
  const myb_shoot_options ob = {.tolerance = 1.0};
  MT_CHECK(myb_shoot(bplane_residual, &bs, 3, 2, dvel, &ob, work, &its, &res) == MYB_OK);
  printf("B-plane targeting (Earth + Moon, variational STM): %d iterations, miss %.2e m, "
         "correction %.4f m/s\n",
         its, res, norm3(dvel));
  MT_CHECK(res <= 1.0 && its <= 8);
  myb_force_model_destroy(m);
  /* errors */
  MT_CHECK(myb_shoot(plane_residual, &mp, 0, 1, w, &o, work, &its, &res) == MYB_E_DOMAIN);
  const myb_shoot_options bad = {.tolerance = 0.0};
  MT_CHECK(myb_shoot(plane_residual, &mp, 3, 1, w, &bad, work, &its, &res) == MYB_E_DOMAIN);
}

/* ---- Budget --------------------------------------------------------------------------------- */

typedef struct j2_term {
  double gm;
  double j2;
  double r;
} j2_term;

static myb_status j2_accel(void* context, double t, const double x[3], const double v[3],
                           double a[3]) {
  (void)t;
  (void)v;
  const j2_term* j = context;
  const double r = norm3(x);
  const double z2 = x[2] * x[2] / (r * r);
  const double k = 1.5 * j->j2 * j->gm * j->r * j->r / (r * r * r * r * r);
  a[0] = k * x[0] * (5.0 * z2 - 1.0);
  a[1] = k * x[1] * (5.0 * z2 - 1.0);
  a[2] = k * x[2] * (5.0 * z2 - 3.0);
  return MYB_OK;
}

static myb_status drag_accel(void* context, double t, const double x[3], const double v[3],
                             double a[3]) {
  (void)context;
  (void)t;
  const double h = norm3(x) - 6378137.0;
  const double rho = 3.7e-12 * exp(-(h - 400e3) / 60e3); /* exponential, 400 km */
  const double w = 7.292115e-5;
  const double vr[3] = {v[0] + w * x[1], v[1] - w * x[0], v[2]};
  const double s = norm3(vr);
  const double k = -0.5 * rho * 2.2 * 0.01 * s;
  for (int i = 0; i < 3; ++i) {
    a[i] = k * vr[i];
  }
  return MYB_OK;
}

typedef struct constant_term {
  double a[3];
} constant_term;

static myb_status constant_accel(void* context, double t, const double x[3], const double v[3],
                                 double a[3]) {
  (void)t;
  (void)x;
  (void)v;
  memcpy(a, ((const constant_term*)context)->a, 3 * sizeof(double));
  return MYB_OK;
}

static void test_budget(void) {
  earth_term et = {gm_earth};
  j2_term j2 = {gm_earth, 1.0826e-3, 6378137.0};
  constant_term srp = {{4e-8, -3e-8, 1e-8}};
  constant_term tiny = {{1e-14, 0.0, 0.0}};
  myb_force_model* m = NULL;
  MT_CHECK(myb_force_model_create(NULL, &m) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "earth", earth_accel, earth_jac, &et, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "J2", j2_accel, NULL, &j2, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "drag", drag_accel, NULL, NULL, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "srp", constant_accel, NULL, &srp, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "tiny", constant_accel, NULL, &tiny, NULL) == MYB_OK);
  myb_force_model_set_enabled(m, 2, 0); /* drag left out */
  const double r = 6378137.0 + 400e3;
  const double vc = sqrt(gm_earth / r);
  const double y0[6] = {r, 0.0, 0.0, 0.0, vc * cos(0.9), vc * sin(0.9)};
  myb_budget_entry e[5];
  const myb_budget_options o = {.tolerance = 1.0, .samples = 48};
  MT_CHECK(myb_budget(m, 0.0, y0, 86400.0, &o, e) == MYB_OK);
  printf("budget of a 400 km orbit over a day (tolerance 1 m):\n");
  for (int k = 0; k < 5; ++k) {
    printf("  %-6s %-3s  max |a| %.3e m/s^2  effect %.3e m  %s\n", e[k].name,
           e[k].enabled ? "on" : "off", e[k].max_acceleration, e[k].effect,
           e[k].flag == MYB_BUDGET_MISSING ? "MISSING"
                                           : (e[k].flag == MYB_BUDGET_NEGLIGIBLE ? "negligible" : ""));
  }
  MT_CHECK(e[1].effect > 1e5 && e[1].flag == MYB_BUDGET_OK);     /* J2: hundreds of km */
  MT_CHECK(e[2].flag == MYB_BUDGET_MISSING && e[2].effect > 100.0); /* drag left out */
  MT_CHECK(e[4].flag == MYB_BUDGET_NEGLIGIBLE && e[4].effect < 0.01);
  MT_CHECK(myb_force_model_enabled(m, 2) == 0); /* configuration restored */
  /* the effect of a small force is linear in it */
  const double first = e[3].effect;
  for (int i = 0; i < 3; ++i) {
    srp.a[i] *= 2.0;
  }
  MT_CHECK(myb_budget(m, 0.0, y0, 86400.0, &o, e) == MYB_OK);
  printf("  doubled pressure: effect x %.4f\n", e[3].effect / first);
  MT_CHECK(fabs(e[3].effect / first - 2.0) < 0.01);
  const myb_budget_options bad = {.tolerance = 0.0};
  MT_CHECK(myb_budget(m, 0.0, y0, 86400.0, &bad, e) == MYB_E_DOMAIN);
  myb_force_model_destroy(m);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_bplane();
  test_patched();
  test_shooting();
  test_budget();
  return mt_finish("test_transfer");
}
