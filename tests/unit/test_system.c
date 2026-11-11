/* SPDX-License-Identifier: BSD-3-Clause */
/* The system of bodies with models (myb/system.h), level by level:
 *  - a field on a point against the same sum made by hand; momentum kept;
 *  - the tide: equal bands without lag against the gradient of k2 G m R^5 P2 / (r^3 r*^3) by
 *    differences; separate bands against differences of their own potential; with a small lag
 *    against Mignard's formula (myb_force_tidal_lag); the earlier position of the raiser
 *    against Kepler motion;
 *  - an integrated rotation: without torque against the exact free top; with the torque of a
 *    point mass, the angular momentum of orbit and spin together;
 *  - the terms beyond 1PN of a pair against myb_pn_binary_evaluate, spins included;
 *  - encounters: the instant of contact of two falling spheres against the closed form, a
 *    merger (mass, momentum, what it takes out), an elastic rebound (energy), the closest
 *    approach of a fly-by against the pericentre of its hyperbola;
 *  - actions: light on the mesh of a turning body — the push and the torque of the moment
 *    against the mesh itself, through the equations of rotation. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/attitude.h"
#include "myb/force.h"
#include "myb/kepler.h"
#include "myb/rel.h"
#include "myb/system.h"

static const double PI = 3.141592653589793;

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static double dot3(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross3(const double a[3], const double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

/* A frame turned by a fixed quaternion. */
typedef struct fixed_frame {
  double m[9];
  int calls;
} fixed_frame;

static myb_status fixed_orientation(void* context, myb_time t, myb_orientation* o) {
  (void)t;
  fixed_frame* f = context;
  ++f->calls;
  memset(o, 0, sizeof *o);
  memcpy(o->rotation.m, f->m, sizeof f->m);
  return MYB_OK;
}

/* ---- Fields ----------------------------------------------------------------------------------- */

static void test_field(void) {
  const double gm[3] = {3.986004418e14, 4.9028e12, 1.0e9};
  /* J2 and a sectorial term, normalised */
  const double c[6] = {1.0, 0.0, 0.0, -4.84165e-4, 0.0, 2.4394e-6};
  const double sn[6] = {0.0, 0.0, 0.0, 0.0, 0.0, -1.4003e-6};
  myb_gravity_field* field = NULL;
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_field_create(gm[0], 6378136.3, 2, c, sn, NULL, &field) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(2, NULL, &work) == MYB_OK);
  fixed_frame frame = {.calls = 0};
  const double q[4] = {0.9, 0.1, -0.3, 0.3};
  myb_quat_to_matrix(q, frame.m);
  myb_system* s = NULL;
  MT_CHECK(myb_system_create(3, gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_system_bodies(s) == 3 && myb_system_coordinates(s) == 9);
  /* a level that needs a frame the body does not have */
  const double x[9] = {1.0e5, -2.0e5, 3.0e4, 2.1e8, 3.0e8, -1.0e8, -4.0e7, 1.0e7, 2.0e7};
  const double v[9] = {0.0, 0.0, 0.0, -700.0, 600.0, 100.0, 2000.0, 1500.0, -300.0};
  double a[9];
  MT_CHECK(myb_system_add_field(s, 0, 1) == MYB_OK);
  MT_CHECK(myb_system_acceleration(s, 0.0, x, v, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_system_set_orientation(s, 0, fixed_orientation, &frame) == MYB_OK);
  MT_CHECK(myb_system_acceleration(s, 0.0, x, v, a) == MYB_E_DOMAIN); /* no field yet */
  MT_CHECK(myb_system_set_field(s, 0, field, 0) == MYB_OK);
  MT_CHECK(myb_system_add_field(s, 0, 2) == MYB_OK);
  MT_CHECK(myb_system_acceleration(s, 0.0, x, v, a) == MYB_OK);
  /* by hand */
  double expected[9];
  MT_CHECK(myb_nbody_acceleration_xv(myb_system_nbody(s), 0.0, x, v, expected) == MYB_OK);
  for (int point = 1; point <= 2; ++point) {
    double d[3];
    double db[3];
    double ab[3];
    for (int i = 0; i < 3; ++i) {
      d[i] = x[3 * point + i] - x[i];
    }
    for (int i = 0; i < 3; ++i) {
      db[i] = frame.m[3 * i] * d[0] + frame.m[3 * i + 1] * d[1] + frame.m[3 * i + 2] * d[2];
    }
    MT_CHECK(myb_gravity_field_eval(field, work, db, 2, -1, NULL, ab) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      const double ai = frame.m[i] * ab[0] + frame.m[3 + i] * ab[1] + frame.m[6 + i] * ab[2];
      expected[3 * point + i] += ai;
      expected[i] -= gm[point] / gm[0] * ai;
    }
  }
  double worst = 0.0;
  double momentum = 0.0;
  for (int i = 0; i < 3; ++i) {
    double sum = 0.0;
    double scale = 0.0;
    for (int b = 0; b < 3; ++b) {
      worst = fmax(worst, fabs(a[3 * b + i] - expected[3 * b + i]) / norm3(expected + 3 * b));
      sum += gm[b] * a[3 * b + i];
      scale += fabs(gm[b] * a[3 * b + i]);
    }
    momentum = fmax(momentum, fabs(sum) / scale);
  }
  printf("field of a body on two points: against the sum by hand %.2e, momentum %.2e\n", worst,
         momentum);
  MT_CHECK(worst < 4e-16 && momentum < 4e-16);
  /* what is refused */
  MT_CHECK(myb_system_add_field(s, 0, 0) == MYB_E_DOMAIN);
  MT_CHECK(myb_system_add_field(s, 3, 0) == MYB_E_DOMAIN);
  MT_CHECK(myb_system_set_field(s, 0, field, 3) == MYB_E_DOMAIN);
  MT_CHECK(myb_system_set_field(s, 0, field, 1) == MYB_E_DOMAIN);
  myb_system_destroy(s);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
}

/* ---- Tides ------------------------------------------------------------------------------------ */

/* The potential of the tide by band, as the header writes it. */
static double tide_potential(const myb_system_tide* tide, double gm_raiser, const double pole[3],
                             const double r[3], const double raiser[3][3]) {
  const double r2 = dot3(r, r);
  const double z = dot3(r, pole);
  double u = 0.0;
  for (int j = 0; j < 3; ++j) {
    const double* rs = raiser[j];
    const double s2 = dot3(rs, rs);
    const double zs = dot3(rs, pole);
    double rho[3];
    double rhos[3];
    for (int i = 0; i < 3; ++i) {
      rho[i] = r[i] - z * pole[i];
      rhos[i] = rs[i] - zs * pole[i];
    }
    const double cross = dot3(rho, rhos);
    double w = 0.0;
    if (j == 0) {
      w = 0.5 * (3.0 * z * z - r2) * (3.0 * zs * zs - s2);
    } else if (j == 1) {
      w = 6.0 * z * zs * cross;
    } else {
      w = 1.5 * (2.0 * cross * cross - dot3(rho, rho) * dot3(rhos, rhos));
    }
    u += tide->k2[j] * w / (s2 * s2 * sqrt(s2));
  }
  const double radius5 = tide->radius * tide->radius * tide->radius * tide->radius * tide->radius;
  return 0.5 * gm_raiser * radius5 * u / (r2 * r2 * sqrt(r2));
}

/* The acceleration of the tide alone on every body: the system less its point masses. */
static void tide_only(myb_system* s, const double* x, const double* v, double* a) {
  double point[9];
  MT_CHECK(myb_system_acceleration(s, 0.0, x, v, a) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration_xv(myb_system_nbody(s), 0.0, x, v, point) == MYB_OK);
  for (int i = 0; i < 9; ++i) {
    a[i] -= point[i];
  }
}

static void test_tide(void) {
  /* a planet, its moon and a far star */
  const double gm[3] = {3.986004418e14, 4.9028e12, 1.32712440018e20};
  fixed_frame frame = {.calls = 0};
  const double q[4] = {0.97, 0.12, -0.2, 0.05};
  myb_quat_to_matrix(q, frame.m);
  const double pole[3] = {frame.m[6], frame.m[7], frame.m[8]};
  const double x[9] = {0.0, 0.0, 0.0, 3.1e8, 2.0e8, 0.9e8, -1.2e11, 0.8e11, 0.3e11};
  const double v[9] = {0.0, 0.0, 0.0, -550.0, 800.0, 120.0, -16000.0, -24000.0, 300.0};
  /* a radius of 200 000 km (the formula does not mind) lifts the tide of the star, 1e-14
   * m/s^2 for the real one, above the rounding of the point masses it is taken out of */
  myb_system_tide tide = {.radius = 2.0e8, .k2 = {0.3, 0.3, 0.3}};
  myb_system* s = NULL;
  MT_CHECK(myb_system_create(3, gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_system_set_orientation(s, 0, fixed_orientation, &frame) == MYB_OK);
  MT_CHECK(myb_system_set_tide(s, 0, &tide) == MYB_OK);
  MT_CHECK(myb_system_add_tide(s, 0, 2, 1) == MYB_OK); /* raised by the star, on the moon */
  double a[9];
  tide_only(s, x, v, a);
  /* equal bands, no lag: the gradient of k2 G m R^5 P2(cos) / (r^3 r*^3) at the moon */
  const double* rs = x + 6;
  double worst = 0.0;
  double grad[3];
  for (int i = 0; i < 3; ++i) {
    double u[2];
    for (int side = 0; side < 2; ++side) {
      double r[3] = {x[3], x[4], x[5]};
      const double h = 2.0e3;
      r[i] += side == 0 ? h : -h;
      const double rr = norm3(r);
      const double ss = norm3(rs);
      const double cosine = dot3(r, rs) / (rr * ss);
      u[side] = 0.3 * gm[2] * pow(tide.radius, 5.0) * 0.5 * (3.0 * cosine * cosine - 1.0) /
                (rr * rr * rr * ss * ss * ss);
    }
    grad[i] = (u[0] - u[1]) / 4.0e3;
  }
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(a[3 + i] - grad[i]) / norm3(grad));
  }
  double momentum = 0.0;
  for (int i = 0; i < 3; ++i) {
    momentum = fmax(momentum, fabs(gm[0] * a[i] + gm[1] * a[3 + i]) / fabs(gm[1] * a[3 + i]));
    MT_CHECK(a[6 + i] == 0.0); /* the raiser itself feels nothing of it */
  }
  printf("tide, equal bands: against the gradient of k2 G m R^5 P2 / (r^3 r*^3) %.2e, "
         "momentum %.2e\n",
         worst, momentum);
  MT_CHECK(worst < 1e-8 && momentum < 1e-8);

  /* separate bands: against differences of the potential written band by band */
  tide.k2[0] = 0.335;
  tide.k2[1] = 0.277;
  tide.k2[2] = 0.283;
  MT_CHECK(myb_system_set_tide(s, 0, &tide) == MYB_OK);
  tide_only(s, x, v, a);
  const double same[3][3] = {{rs[0], rs[1], rs[2]}, {rs[0], rs[1], rs[2]}, {rs[0], rs[1], rs[2]}};
  worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    double r[3] = {x[3], x[4], x[5]};
    r[i] += 2.0e3;
    const double up = tide_potential(&tide, gm[2], pole, r, same);
    r[i] -= 4.0e3;
    const double down = tide_potential(&tide, gm[2], pole, r, same);
    grad[i] = (up - down) / 4.0e3;
  }
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(a[3 + i] - grad[i]) / norm3(grad));
  }
  printf("tide, three bands: against differences of the potential %.2e\n", worst);
  MT_CHECK(worst < 1e-8);
  myb_system_destroy(s);

  /* the tide a moon raises acting on that moon, with a small lag: Mignard's formula. The
   * planet turns about its pole; the bands answer alike, the orbital and the rotational lag
   * are the same time */
  const double gm2[2] = {gm[0], gm[1]};
  const double rate = 7.292115e-5;
  const double lag = 60.0;
  myb_system_tide lagged = {.radius = 2.0e8, .k2 = {0.3, 0.3, 0.3},
                            .orbital_lag = {lag, lag, lag}, .rotation_lag = {0.0, lag, lag},
                            .rotation_rate = rate};
  MT_CHECK(myb_system_create(2, gm2, NULL, &s) == MYB_OK);
  MT_CHECK(myb_system_set_orientation(s, 0, fixed_orientation, &frame) == MYB_OK);
  MT_CHECK(myb_system_set_tide(s, 0, &lagged) == MYB_OK);
  MT_CHECK(myb_system_add_tide(s, 0, 1, 1) == MYB_OK);
  double a2[6];
  double point[6];
  MT_CHECK(myb_system_acceleration(s, 0.0, x, v, a2) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration_xv(myb_system_nbody(s), 0.0, x, v, point) == MYB_OK);
  double relative_tide[3];
  for (int i = 0; i < 3; ++i) {
    relative_tide[i] = (a2[3 + i] - point[3 + i]) - (a2[i] - point[i]);
  }
  const double omega[3] = {rate * pole[0], rate * pole[1], rate * pole[2]};
  double mignard[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_tidal_lag(gm2[0], gm2[1], 0.3, lagged.radius, lag, omega, x + 3, v + 3, 1,
                               mignard) == MYB_OK);
  double static_part[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_tidal_lag(gm2[0], gm2[1], 0.3, lagged.radius, 0.0, omega, x + 3, v + 3, 1,
                               static_part) == MYB_OK);
  double lag_part[3];
  double lag_ours[3];
  worst = 0.0;
  /* the lagging parts compared: the series of Mignard stops at the first order of the lag,
   * whose angle is 60 s * omega = 4e-3 */
  MT_CHECK(myb_system_set_tide(s, 0, &(myb_system_tide){.radius = lagged.radius,
                                                        .k2 = {0.3, 0.3, 0.3}}) == MYB_OK);
  double a0[6];
  MT_CHECK(myb_system_acceleration(s, 0.0, x, v, a0) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    lag_part[i] = mignard[i] - static_part[i];
    lag_ours[i] = relative_tide[i] - ((a0[3 + i] - point[3 + i]) - (a0[i] - point[i]));
  }
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(lag_ours[i] - lag_part[i]) / norm3(lag_part));
  }
  printf("tide with a lag of 60 s: the lagging part against Mignard's first order %.2e (the "
         "angle of the lag is %.1e)\n",
         worst, rate * lag);
  MT_CHECK(worst < 3.0 * rate * lag);
  /* the lag turns the orbit the way the planet turns: a pull along the motion when the
   * planet turns faster than the moon goes round */
  double h[3];
  cross3(x + 3, v + 3, h);
  MT_CHECK(dot3(h, pole) > 0.0);
  double along[3];
  cross3(pole, x + 3, along);
  MT_CHECK(dot3(lag_ours, along) > 0.0);
  myb_system_destroy(s);

  /* where the raiser was: the series of the system against Kepler motion of the pair, a lag
   * of a quarter of an hour on an orbit of 27 days */
  myb_system_tide early = {.radius = 2.0e8, .k2 = {0.3, 0.0, 0.0},
                           .orbital_lag = {900.0, 0.0, 0.0}};
  MT_CHECK(myb_system_create(2, gm2, NULL, &s) == MYB_OK);
  MT_CHECK(myb_system_set_orientation(s, 0, fixed_orientation, &frame) == MYB_OK);
  MT_CHECK(myb_system_set_tide(s, 0, &early) == MYB_OK);
  MT_CHECK(myb_system_add_tide(s, 0, 1, 1) == MYB_OK);
  MT_CHECK(myb_system_acceleration(s, 0.0, x, v, a2) == MYB_OK);
  const myb_pv now = {{x[3], x[4], x[5]}, {v[3], v[4], v[5]}};
  myb_pv before;
  MT_CHECK(myb_kepler_propagate(gm2[0] + gm2[1], &now, -900.0, &before) == MYB_OK);
  const double there[3][3] = {{before.r[0], before.r[1], before.r[2]},
                              {before.r[0], before.r[1], before.r[2]},
                              {before.r[0], before.r[1], before.r[2]}};
  worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    double r[3] = {x[3], x[4], x[5]};
    r[i] += 2.0e3;
    const double up = tide_potential(&early, gm2[1], pole, r, there);
    r[i] -= 4.0e3;
    const double down = tide_potential(&early, gm2[1], pole, r, there);
    grad[i] = (up - down) / 4.0e3;
  }
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs((a2[3 + i] - point[3 + i]) - grad[i]) / norm3(grad));
  }
  printf("tide of a raiser 900 s ago: the series against Kepler motion %.2e\n", worst);
  MT_CHECK(worst < 1e-8);
  /* refused: a tide for a body whose rotation is integrated, and the other way round */
  const myb_system_rotation top = {.moments = {1.0, 2.0, 3.0}};
  MT_CHECK(myb_system_set_rotation(s, 0, &top) == MYB_E_DOMAIN);
  MT_CHECK(myb_system_set_rotation(s, 1, &top) == MYB_OK);
  MT_CHECK(myb_system_set_tide(s, 1, &early) == MYB_E_DOMAIN);
  MT_CHECK(myb_system_add_tide(s, 1, 0, 0) == MYB_E_DOMAIN);
  myb_system_destroy(s);
}

/* ---- Rotation --------------------------------------------------------------------------------- */

static void test_rotation(void) {
  /* a top far from everything: the exact free rotation */
  const double gm[2] = {4.9028e12, 3.986004418e14};
  const myb_system_rotation top = {.moments = {0.9e12, 1.0e12, 1.3e12}};
  myb_system* s = NULL;
  MT_CHECK(myb_system_create(2, gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_system_set_rotation(s, 0, &top) == MYB_OK);
  const size_t n = myb_system_coordinates(s);
  MT_CHECK(n == 6 + 9);
  double y[2 * 15];
  memset(y, 0, sizeof y);
  const double x0[6] = {0.0, 0.0, 0.0, 3.0e8, 2.0e8, 1.0e8};
  const double v0[6] = {0.0, 0.0, 0.0, -600.0, 750.0, 200.0};
  memcpy(y, x0, sizeof x0);
  memcpy(y + n, v0, sizeof v0);
  double q0[4] = {0.8, 0.2, -0.4, 0.3};
  (void)myb_quat_normalize(q0);
  const double w0[3] = {1.0e-5, -2.0e-5, 6.0e-5};
  MT_CHECK(myb_system_pack_rotation(s, 0, q0, w0, NULL, y) == MYB_OK);
  double q[4];
  double w[3];
  double wc[3];
  MT_CHECK(myb_system_unpack_rotation(s, 0, y, q, w, wc) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(w[i] - w0[i]) < 1e-20 && wc[i] == 0.0);
  }
  const double l0[3] = {top.moments[0] * w0[0], top.moments[1] * w0[1], top.moments[2] * w0[2]};
  myb_free_rotation exact;
  MT_CHECK(myb_free_rotation_init(top.moments, q0, l0, &exact) == MYB_OK);
  const myb_ias15_options o = {.epsilon = 1e-10};
  myb_ias15* in = NULL;
  MT_CHECK(myb_system_integrator(s, &o, NULL, &in) == MYB_OK);
  MT_CHECK(myb_ias15_start(in, myb_system_acceleration, s, 0.0, y) == MYB_OK);
  const double t1 = 30.0 * 86400.0; /* some 25 turns */
  MT_CHECK(myb_ias15_integrate(in, t1) == MYB_OK);
  double qe[4];
  double le[3];
  myb_free_rotation_state(&exact, t1, qe, le);
  MT_CHECK(myb_system_unpack_rotation(s, 0, myb_ias15_state(in), q, w, NULL) == MYB_OK);
  double worst = 0.0;
  const double sign = q[0] * qe[0] + q[1] * qe[1] + q[2] * qe[2] + q[3] * qe[3] < 0.0 ? -1.0 : 1.0;
  for (int i = 0; i < 4; ++i) {
    worst = fmax(worst, fabs(q[i] - sign * qe[i]));
  }
  double worst_l = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst_l = fmax(worst_l, fabs(top.moments[i] * w[i] - le[i]) / norm3(le));
  }
  const double* state = myb_ias15_state(in);
  const double drift = fabs(sqrt(state[6] * state[6] + state[7] * state[7] + state[8] * state[8] +
                                 state[9] * state[9]) - 1.0);
  printf("free top over 30 days (%lld steps): attitude against the exact motion %.2e, angular "
         "momentum %.2e, |q| - 1 = %.1e\n",
         (long long)myb_ias15_statistics(in).accepted, worst, worst_l, drift);
  MT_CHECK(worst < 2e-11 && worst_l < 2e-11 && drift < 1e-13);
  myb_ias15_destroy(in);
  myb_system_destroy(s);

  /* the same top turned by a point mass through its figure: the orbit and the spin exchange
   * angular momentum and keep their sum; the pair keeps its momentum */
  MT_CHECK(myb_system_create(2, gm, NULL, &s) == MYB_OK);
  /* moments of a body of 1700 km: 0.39 R^2 and a triaxial figure */
  const double r2 = 1.7e6 * 1.7e6;
  const myb_system_rotation moon = {.moments = {0.390 * r2, 0.3902 * r2, 0.3905 * r2}};
  MT_CHECK(myb_system_set_rotation(s, 0, &moon) == MYB_OK);
  MT_CHECK(myb_system_add_field(s, 0, 1) == MYB_OK);
  memset(y, 0, sizeof y);
  /* the pair about its centre of mass, 60 000 km apart */
  const double total = gm[0] + gm[1];
  const double sep[3] = {6.0e7, 0.0, 0.0};
  const double speed = sqrt(total / 6.0e7);
  const double vel[3] = {0.0, speed * 0.9, speed * 0.3};
  for (int i = 0; i < 3; ++i) {
    y[i] = -gm[1] / total * sep[i];
    y[3 + i] = gm[0] / total * sep[i];
    y[n + (size_t)i] = -gm[1] / total * vel[i];
    y[n + 3 + (size_t)i] = gm[0] / total * vel[i];
  }
  const double spin[3] = {2.0e-6, -1.0e-6, 2.6e-5};
  MT_CHECK(myb_system_pack_rotation(s, 0, q0, spin, NULL, y) == MYB_OK);
  double first[3] = {0.0, 0.0, 0.0};
  double worst_sum = 0.0;
  double exchanged = 0.0;
  double spin_first[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_system_integrator(s, &o, NULL, &in) == MYB_OK);
  MT_CHECK(myb_ias15_start(in, myb_system_acceleration, s, 0.0, y) == MYB_OK);
  for (int k = 0; k <= 40; ++k) {
    if (k > 0) {
      MT_CHECK(myb_ias15_integrate(in, k * 86400.0) == MYB_OK);
    }
    const double* z = myb_ias15_state(in);
    double sum[3] = {0.0, 0.0, 0.0};
    for (int b = 0; b < 2; ++b) {
      double l[3];
      cross3(z + 3 * b, z + n + 3 * b, l);
      for (int i = 0; i < 3; ++i) {
        sum[i] += gm[b] * l[i];
      }
    }
    /* the spin: G m M^T (I w), the inertia tensor being that of the state */
    double inertia[9];
    MT_CHECK(myb_system_inertia(s, 0, k * 86400.0, z, inertia) == MYB_OK);
    MT_CHECK(myb_system_unpack_rotation(s, 0, z, q, w, NULL) == MYB_OK);
    double m[9];
    myb_quat_to_matrix(q, m);
    double lb[3];
    double li[3];
    for (int i = 0; i < 3; ++i) {
      lb[i] = inertia[3 * i] * w[0] + inertia[3 * i + 1] * w[1] + inertia[3 * i + 2] * w[2];
    }
    for (int i = 0; i < 3; ++i) {
      li[i] = gm[0] * (m[i] * lb[0] + m[3 + i] * lb[1] + m[6 + i] * lb[2]);
      sum[i] += li[i];
    }
    if (k == 0) {
      memcpy(first, sum, sizeof first);
      memcpy(spin_first, li, sizeof li);
    }
    for (int i = 0; i < 3; ++i) {
      worst_sum = fmax(worst_sum, fabs(sum[i] - first[i]) / norm3(first));
      exchanged = fmax(exchanged, fabs(li[i] - spin_first[i]) / norm3(first));
    }
  }
  const double* z = myb_ias15_state(in);
  double p = 0.0;
  for (int i = 0; i < 3; ++i) {
    p = fmax(p, fabs(gm[0] * z[n + (size_t)i] + gm[1] * z[n + 3 + (size_t)i]) / (gm[0] * speed));
  }
  printf("a figure turned by a point mass, 40 days: orbit and spin exchange %.2e of the angular "
         "momentum, their sum keeps to %.2e; momentum %.2e\n",
         exchanged, worst_sum, p);
  MT_CHECK(exchanged > 1e-10);
  MT_CHECK(worst_sum < 1e-4 * exchanged && worst_sum < 1e-13);
  MT_CHECK(p < 1e-14);
  myb_ias15_destroy(in);
  /* refused */
  MT_CHECK(myb_system_pack_rotation(s, 1, q0, spin, NULL, y) == MYB_E_DOMAIN);
  const double zero[4] = {0.0, 0.0, 0.0, 0.0};
  MT_CHECK(myb_system_pack_rotation(s, 0, zero, spin, NULL, y) == MYB_E_DOMAIN);
  const myb_system_rotation bad = {.moments = {1.0, 0.0, 3.0}};
  MT_CHECK(myb_system_set_rotation(s, 1, &bad) == MYB_E_DOMAIN);
  const myb_system_rotation soft = {.moments = {1.0, 2.0, 3.0}, .k2 = 0.02};
  MT_CHECK(myb_system_set_rotation(s, 1, &soft) == MYB_E_DOMAIN); /* k2 without a radius */
  MT_CHECK(myb_system_add_distortion(s, 1, 0) == MYB_E_DOMAIN);   /* its rotation is not integrated */
  myb_system_destroy(s);
}

/* ---- Pairs ------------------------------------------------------------------------------------ */

static void test_pn_pair(void) {
  const double c = 299792458.0;
  const double gm[3] = {1.4 * 1.32712440018e20, 1.3 * 1.32712440018e20, 1.0e15};
  const double x[9] = {4.0e8, 1.0e8, -2.0e8, -5.0e8, -1.0e8, 3.0e8, 1.0e11, 2.0e11, 0.0};
  const double v[9] = {-1.0e5, 2.5e5, 0.5e5, 1.2e5, -2.6e5, -0.6e5, 0.0, 1.0e3, 0.0};
  myb_system* s = NULL;
  MT_CHECK(myb_system_create(3, gm, NULL, &s) == MYB_OK);
  /* no speed of light yet */
  myb_system_pn_pair model = {.orders = MYB_PN_1 | MYB_PN_2 | MYB_PN_25};
  MT_CHECK(myb_system_add_pn_pair(s, 0, 1, &model) == MYB_E_DOMAIN);
  model.c = c;
  MT_CHECK(myb_system_add_pn_pair(s, 0, 1, &model) == MYB_OK);
  MT_CHECK(myb_system_coordinates(s) == 9);
  double a[9];
  double point[9];
  MT_CHECK(myb_system_acceleration(s, 0.0, x, v, a) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration_xv(myb_system_nbody(s), 0.0, x, v, point) == MYB_OK);
  double xr[3];
  double vr[3];
  for (int i = 0; i < 3; ++i) {
    xr[i] = x[i] - x[3 + i];
    vr[i] = v[i] - v[3 + i];
  }
  double whole[3];
  double newton[3];
  MT_CHECK(myb_pn_binary_acceleration(gm[0], gm[1], c, model.orders, xr, vr, whole) == MYB_OK);
  MT_CHECK(myb_pn_binary_acceleration(gm[0], gm[1], c, 0u, xr, vr, newton) == MYB_OK);
  double worst = 0.0;
  double momentum = 0.0;
  double size = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double ours = (a[i] - point[i]) - (a[3 + i] - point[3 + i]);
    worst = fmax(worst, fabs(ours - (whole[i] - newton[i])) / norm3(newton));
    size = fmax(size, fabs(whole[i] - newton[i]) / norm3(newton));
    momentum = fmax(momentum, fabs(gm[0] * (a[i] - point[i]) + gm[1] * (a[3 + i] - point[3 + i])) /
                                  (gm[0] * norm3(newton)));
    MT_CHECK(a[6 + i] == point[6 + i]);
  }
  printf("terms beyond Newton of a pair among three bodies (%.1e of the attraction): against "
         "the binary %.2e, momentum %.2e\n",
         size, worst, momentum);
  MT_CHECK(worst < 4e-16 && momentum < 1e-15 && size > 1e-7);
  myb_system_destroy(s);

  /* with spins: their rates are those of the binary, and they are part of the state */
  MT_CHECK(myb_system_create(3, gm, NULL, &s) == MYB_OK);
  myb_system_pn_pair spinning = {.orders = MYB_PN_ALL, .c = c, .r0 = 1.0e7, .kappa1 = 4.0,
                                 .kappa2 = 5.0};
  MT_CHECK(myb_system_add_pn_pair(s, 1, 0, &spinning) == MYB_OK);
  const size_t n = myb_system_coordinates(s);
  MT_CHECK(n == 15);
  double y[30];
  memset(y, 0, sizeof y);
  memcpy(y, x, sizeof x);
  memcpy(y + n, v, sizeof v);
  const double chi1[3] = {0.1, -0.2, 0.3};
  const double chi2[3] = {-0.4, 0.05, 0.2};
  MT_CHECK(myb_system_pack_spins(s, 0, chi1, chi2, y) == MYB_OK);
  double back1[3];
  double back2[3];
  MT_CHECK(myb_system_unpack_spins(s, 0, y, back1, back2) == MYB_OK);
  MT_CHECK(memcmp(back1, chi1, sizeof chi1) == 0 && memcmp(back2, chi2, sizeof chi2) == 0);
  double full[15];
  MT_CHECK(myb_system_acceleration(s, 0.0, y, y + n, full) == MYB_OK);
  const myb_pn_binary binary = {.gm1 = gm[1], .gm2 = gm[0], .c = c, .orders = MYB_PN_ALL,
                                .r0 = 1.0e7, .kappa1 = 4.0, .kappa2 = 5.0};
  double d1[3];
  double d2[3];
  for (int i = 0; i < 3; ++i) {
    xr[i] = -xr[i];
    vr[i] = -vr[i];
  }
  MT_CHECK(myb_pn_binary_evaluate(&binary, xr, vr, chi1, chi2, whole, d1, d2) == MYB_OK);
  worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(full[9 + i] - d1[i]) / norm3(d1));
    worst = fmax(worst, fabs(full[12 + i] - d2[i]) / norm3(d2));
    const double ours = (full[3 + i] - point[3 + i]) - (full[i] - point[i]);
    worst = fmax(worst, fabs(ours - (whole[i] + newton[i])) / norm3(newton));
  }
  printf("  with spins: rates of the spins and the acceleration against the binary %.2e\n", worst);
  MT_CHECK(worst < 4e-16);
  MT_CHECK(myb_system_pack_spins(s, 1, chi1, chi2, y) == MYB_E_DOMAIN);
  const myb_system_pn_pair nothing = {.orders = 0u, .c = c};
  MT_CHECK(myb_system_add_pn_pair(s, 0, 2, &nothing) == MYB_E_DOMAIN);
  const myb_system_pn_pair gauge = {.orders = MYB_PN_3, .c = c};
  MT_CHECK(myb_system_add_pn_pair(s, 0, 2, &gauge) == MYB_E_DOMAIN);
  myb_system_destroy(s);
}

/* ---- Encounters ------------------------------------------------------------------------------- */

static void test_encounters(void) {
  /* two spheres let go at rest: the instant of contact is that of the radial Kepler fall,
   * t = sqrt(d^3 / (8 mu)) (acos(sqrt(w)) + sqrt(w (1 - w))) 2, w = contact / d */
  const double gm[2] = {4.0e5, 1.0e5};
  myb_system* s = NULL;
  MT_CHECK(myb_system_create(2, gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_system_set_radius(s, 0, 60.0) == MYB_OK);
  MT_CHECK(myb_system_set_radius(s, 1, 40.0) == MYB_OK);
  const double d = 5000.0;
  const double total = gm[0] + gm[1];
  double y[12] = {-gm[1] / total * d, 0.0, 0.0, gm[0] / total * d, 0.0, 0.0,
                  0.0,                0.0, 0.0, 0.0,               0.0, 0.0};
  const myb_ias15_options o = {.epsilon = 1e-10};
  myb_ias15* in = NULL;
  MT_CHECK(myb_system_integrator(s, &o, NULL, &in) == MYB_OK);
  MT_CHECK(myb_ias15_start(in, myb_system_acceleration, s, 0.0, y) == MYB_OK);
  myb_system_event event;
  MT_CHECK(myb_system_advance(s, in, 1.0e6, 0.0, &event) == MYB_OK);
  const double w = 100.0 / d;
  const double fall = sqrt(d * d * d / (2.0 * total)) * (acos(sqrt(w)) + sqrt(w * (1.0 - w)));
  printf("two spheres falling together: contact at %.9f s (closed form %.9f), centres %.12f m "
         "apart\n",
         event.t, fall, event.distance);
  MT_CHECK(event.kind == MYB_SYSTEM_CONTACT && event.a == 0 && event.b == 1);
  MT_CHECK(fabs(event.t / fall - 1.0) < 1e-12);
  MT_CHECK(fabs(event.distance - 100.0) < 1e-9 && event.distance <= 100.0);
  MT_CHECK(myb_ias15_time(in) == event.t);
  /* an elastic rebound: back to where they started, at rest, after the same time */
  const double* z = myb_ias15_state(in);
  const double speed_in = fabs(z[9] - z[6]);
  MT_CHECK(myb_system_bounce(s, in, 0, 1, 1.0) == MYB_OK);
  z = myb_ias15_state(in);
  MT_CHECK(fabs((z[9] - z[6]) - speed_in) < 1e-12 * speed_in);
  MT_CHECK(fabs(gm[0] * z[6] + gm[1] * z[9]) < 1e-12 * gm[0] * speed_in);
  MT_CHECK(myb_system_advance(s, in, 2.0 * fall, 0.0, &event) == MYB_OK);
  MT_CHECK(event.kind == MYB_SYSTEM_ARRIVED);
  z = myb_ias15_state(in);
  printf("  after an elastic rebound and the same time: %.9f m apart (5000), speed %.1e m/s\n",
         z[3] - z[0], fabs(z[9] - z[6]));
  MT_CHECK(fabs((z[3] - z[0]) - d) < 1e-7 && fabs(z[9] - z[6]) < 1e-9);
  /* they fall again and merge */
  MT_CHECK(myb_system_advance(s, in, 1.0e6, 0.0, &event) == MYB_OK);
  MT_CHECK(event.kind == MYB_SYSTEM_CONTACT);
  MT_CHECK(fabs(event.t / (3.0 * fall) - 1.0) < 1e-10);
  z = myb_ias15_state(in);
  const double centre = (gm[0] * z[0] + gm[1] * z[3]) / total;
  const double speed = fabs(z[9] - z[6]);
  myb_system_merger merger;
  MT_CHECK(myb_system_merge(s, in, 0, 1, &merger) == MYB_OK);
  z = myb_ias15_state(in);
  MT_CHECK(myb_system_present(s, 0) && !myb_system_present(s, 1));
  MT_CHECK(myb_nbody_gm(myb_system_nbody(s), 0) == total && myb_nbody_gm(myb_system_nbody(s), 1) == 0.0);
  MT_CHECK(fabs(z[0] - centre) < 1e-12 && fabs(z[6]) < 1e-12 * speed);
  MT_CHECK(fabs(myb_system_radius(s, 0) - cbrt(60.0 * 60.0 * 60.0 + 40.0 * 40.0 * 40.0)) < 1e-12);
  /* what the merger took: the energy the pair had at rest, 5000 m apart */
  printf("  merger: energy taken %.9e (the pair at rest: %.9e), angular momentum %.1e\n",
         merger.energy, -gm[0] * gm[1] / d, norm3(merger.angular_momentum));
  MT_CHECK(fabs(merger.energy / (-gm[0] * gm[1] / d) - 1.0) < 1e-9);
  MT_CHECK(norm3(merger.angular_momentum) == 0.0);
  /* the merged body goes on alone */
  MT_CHECK(myb_system_advance(s, in, event.t + 100.0, 0.0, &event) == MYB_OK);
  MT_CHECK(event.kind == MYB_SYSTEM_ARRIVED);
  MT_CHECK(myb_system_merge(s, in, 0, 1, NULL) == MYB_E_DOMAIN); /* absent */
  myb_ias15_destroy(in);
  myb_system_destroy(s);

  /* a fly-by: the closest approach against the pericentre of the hyperbola; with radii that
   * reach it, the contact on the way in */
  const double gm2[2] = {3.986004418e14, 0.0};
  MT_CHECK(myb_system_create(2, gm2, NULL, &s) == MYB_OK);
  const double rp = 7.0e6;
  const double e = 1.6;
  const double p = rp * (1.0 + e);
  const double nu0 = -2.0; /* rad: on the way in */
  const double r0 = p / (1.0 + e * cos(nu0));
  const double k = sqrt(gm2[0] / p);
  double y2[12] = {0.0, 0.0, 0.0, r0 * cos(nu0), r0 * sin(nu0), 0.0,
                   0.0, 0.0, 0.0, -k * sin(nu0), k * (e + cos(nu0)), 0.0};
  MT_CHECK(myb_system_integrator(s, &o, NULL, &in) == MYB_OK);
  MT_CHECK(myb_ias15_start(in, myb_system_acceleration, s, 0.0, y2) == MYB_OK);
  MT_CHECK(myb_system_advance(s, in, 1.0e5, 1.0e7, &event) == MYB_OK);
  /* the time from nu0 to the pericentre: Kepler's equation of the hyperbola */
  const double f0 = 2.0 * atanh(sqrt((e - 1.0) / (e + 1.0)) * tan(0.5 * nu0));
  const double a_h = rp / (e - 1.0);
  const double t_peri = -sqrt(a_h * a_h * a_h / gm2[0]) * (e * sinh(f0) - f0);
  printf("fly-by: closest approach at %.9f s (hyperbola %.9f), %.6f m (%.0f)\n", event.t, t_peri,
         event.distance, rp);
  MT_CHECK(event.kind == MYB_SYSTEM_APPROACH && event.a == 0 && event.b == 1);
  MT_CHECK(fabs(event.t / t_peri - 1.0) < 1e-11);
  MT_CHECK(fabs(event.distance / rp - 1.0) < 1e-13);
  /* going on: it is not found twice, and nothing else happens */
  MT_CHECK(myb_system_advance(s, in, 2.0 * t_peri, 1.0e7, &event) == MYB_OK);
  MT_CHECK(event.kind == MYB_SYSTEM_ARRIVED);
  /* beyond the watched distance: not reported */
  MT_CHECK(myb_ias15_start(in, myb_system_acceleration, s, 0.0, y2) == MYB_OK);
  MT_CHECK(myb_system_advance(s, in, 2.0 * t_peri, 6.9e6, &event) == MYB_OK);
  MT_CHECK(event.kind == MYB_SYSTEM_ARRIVED);
  /* spheres of 7500 km together: the contact comes before the pericentre, where
   * r = 7500 km on the way in */
  MT_CHECK(myb_system_set_radius(s, 0, 6.4e6) == MYB_OK);
  MT_CHECK(myb_system_set_radius(s, 1, 1.1e6) == MYB_OK);
  MT_CHECK(myb_ias15_start(in, myb_system_acceleration, s, 0.0, y2) == MYB_OK);
  MT_CHECK(myb_system_advance(s, in, 2.0 * t_peri, 1.0e7, &event) == MYB_OK);
  const double nu_c = -acos((p / 7.5e6 - 1.0) / e);
  const double f_c = 2.0 * atanh(sqrt((e - 1.0) / (e + 1.0)) * tan(0.5 * nu_c));
  const double t_c = t_peri + sqrt(a_h * a_h * a_h / gm2[0]) * (e * sinh(f_c) - f_c);
  printf("  with spheres reaching 7500 km: contact at %.9f s (hyperbola %.9f)\n", event.t, t_c);
  MT_CHECK(event.kind == MYB_SYSTEM_CONTACT);
  MT_CHECK(fabs(event.t / t_c - 1.0) < 1e-11);
  MT_CHECK(fabs(event.distance - 7.5e6) < 1e-6);
  /* left alone, the pair is not reported again while it overlaps; then it parts */
  MT_CHECK(myb_system_advance(s, in, 2.0 * t_peri, 0.0, &event) == MYB_OK);
  MT_CHECK(event.kind == MYB_SYSTEM_ARRIVED);
  MT_CHECK(myb_system_advance(s, in, 0.0, -1.0, &event) == MYB_E_DOMAIN);
  MT_CHECK(myb_system_bounce(s, in, 0, 1, 1.5) == MYB_E_DOMAIN);
  myb_ias15_destroy(in);
  myb_system_destroy(s);
}

/* ---- Actions ---------------------------------------------------------------------------------- */

static myb_status push(void* context, double t, const myb_system_view* view, double acceleration[3],
                       double torque[3]) {
  (void)t;
  const double* amount = context;
  /* towards body 0, and a torque about the third axis of the body */
  for (int i = 0; i < 3; ++i) {
    acceleration[i] = *amount * (view->x[i] - view->x[3 * view->body + (size_t)i]);
  }
  torque[2] = *amount;
  return view->frame != NULL ? MYB_OK : MYB_E_DOMAIN;
}

static void test_actions(void) {
  const double c = 299792458.0;
  const double au = 1.495978707e11;
  /* a tetrahedron, faces counter-clockwise seen from outside */
  const double vertices[12] = {1.0, 0.0, -0.3, -0.5, 0.9, -0.3, -0.6, -0.8, -0.3, 0.1, 0.2, 1.0};
  const size_t faces[12] = {0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3};
  myb_polyhedron* shape = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, 4, faces, 4, 1.0, NULL, &shape) == MYB_OK);
  const double gm[2] = {1.32712440018e20, 3.3e-8}; /* the Sun and half a tonne */
  const myb_system_rotation top = {.moments = {0.30, 0.35, 0.42}};
  myb_system* s = NULL;
  MT_CHECK(myb_system_create(2, gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_system_set_rotation(s, 1, &top) == MYB_OK);
  myb_system_light light = {.source = 0, .shape = shape,
                            .surface = {.specular = 0.1, .diffuse = 0.2, .thermal = 1},
                            .flux = 1361.0, .distance = au, .c = c, .mass = 500.0};
  MT_CHECK(myb_system_add_action(s, 1, myb_system_light_action, &light) == MYB_OK);
  const size_t n = myb_system_coordinates(s);
  double y[2 * 15];
  memset(y, 0, sizeof y);
  const double place[3] = {0.9 * au, -0.7 * au, 0.2 * au};
  const double speed[3] = {1.6e4, 2.0e4, -3.0e3};
  memcpy(y + 3, place, sizeof place);
  memcpy(y + n + 3, speed, sizeof speed);
  double q[4] = {0.7, -0.3, 0.5, 0.4};
  (void)myb_quat_normalize(q);
  const double w[3] = {3.0e-4, -1.0e-4, 8.0e-4};
  MT_CHECK(myb_system_pack_rotation(s, 1, q, w, NULL, y) == MYB_OK);
  double a[15];
  double point[6];
  MT_CHECK(myb_system_acceleration(s, 0.0, y, y + n, a) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration_xv(myb_system_nbody(s), 0.0, y, y + n, point) == MYB_OK);
  /* by hand: the direction of the light in the frame of the body, the mesh, and back */
  double m[9];
  myb_quat_to_matrix(q, m);
  const double r = norm3(place);
  double u[3];
  for (int i = 0; i < 3; ++i) {
    u[i] = (m[3 * i] * place[0] + m[3 * i + 1] * place[1] + m[3 * i + 2] * place[2]) / r;
  }
  const double pressure = 1361.0 * (au / r) * (au / r) / c;
  double force[3];
  double moment[3];
  MT_CHECK(myb_force_radiation_mesh(shape, &light.surface, pressure, u, 0, force, moment) ==
           MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double expected = (m[i] * force[0] + m[3 + i] * force[1] + m[6 + i] * force[2]) / 500.0;
    worst = fmax(worst, fabs((a[3 + i] - point[3 + i]) - expected) / (norm3(force) / 500.0));
  }
  /* the torque through the equations of rotation: (0, w') = conj(q) (2 q'' - q' (0, w)), and
   * I w' + w x I w is the torque per unit mass */
  const double* dq = y + n + 6;
  const double w4[4] = {0.0, w[0], w[1], w[2]};
  double product[4];
  double rest[4];
  double conj[4];
  double rate4[4];
  myb_quat_multiply(dq, w4, product);
  for (int i = 0; i < 4; ++i) {
    rest[i] = 2.0 * a[6 + i] - product[i];
  }
  myb_quat_conjugate(q, conj);
  myb_quat_multiply(conj, rest, rate4);
  const double lw[3] = {top.moments[0] * w[0], top.moments[1] * w[1], top.moments[2] * w[2]};
  double gyro[3];
  cross3(w, lw, gyro);
  double worst_torque = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double felt = top.moments[i] * rate4[1 + i] + gyro[i];
    worst_torque = fmax(worst_torque, fabs(felt - moment[i] / 500.0) / (norm3(moment) / 500.0));
  }
  printf("light on the mesh of a turning body (%.2e m/s^2, torque %.2e m^2/s^2 per unit "
         "mass): the push against the mesh %.2e, the torque through the equations of rotation "
         "%.2e\n",
         norm3(force) / 500.0, norm3(moment) / 500.0, worst, worst_torque);
  MT_CHECK(worst < 1e-9); /* a difference of accelerations 1e7 times larger */
  MT_CHECK(worst_torque < 1e-9);
  MT_CHECK(fabs(rate4[0]) < 1e-20);
  /* what the light absorbed adds about the axis, over a turn with the Sun at one colatitude:
   * the mean that myb_force_yorp gives is that of the re-emission alone */
  const myb_plate warm = {.thermal = 1};
  const myb_plate black = {0};
  const int phases = 96;
  const double colatitude = 1.1;
  double absorbed = 0.0;
  double emitted = 0.0;
  for (int k = 0; k < phases; ++k) {
    const double phi = 2.0 * PI * (k + 0.5) / phases;
    const double dir[3] = {-sin(colatitude) * cos(phi), sin(colatitude) * sin(phi),
                           -cos(colatitude)};
    double all[3];
    double dark[3];
    MT_CHECK(myb_force_radiation_mesh(shape, &warm, 1.0, dir, 0, NULL, all) == MYB_OK);
    MT_CHECK(myb_force_radiation_mesh(shape, &black, 1.0, dir, 0, NULL, dark) == MYB_OK);
    absorbed += dark[2] / phases;
    emitted += (all[2] - dark[2]) / phases;
  }
  double yorp[3];
  MT_CHECK(myb_force_yorp(shape, 1.0, colatitude, 0.0, 0, phases, yorp) == MYB_OK);
  printf("  about the axis over a turn (p R^3): re-emitted %.6e (myb_force_yorp %.6e), "
         "absorbed %.6e\n",
         emitted, yorp[2], absorbed);
  MT_CHECK(fabs(emitted - yorp[2]) < 1e-14);
  /* any function: a pull and a torque, added to what is there */
  double amount = 1.0e-18;
  MT_CHECK(myb_system_add_action(s, 1, push, &amount) == MYB_OK);
  double b[15];
  MT_CHECK(myb_system_acceleration(s, 0.0, y, y + n, b) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs((b[3 + i] - a[3 + i]) + amount * place[i]) < 1e-12 * amount * r);
  }
  /* refused: a body without a frame for the light, a source that is the body itself */
  MT_CHECK(myb_system_add_action(s, 2, push, &amount) == MYB_E_DOMAIN);
  myb_system* bare = NULL;
  MT_CHECK(myb_system_create(2, gm, NULL, &bare) == MYB_OK);
  MT_CHECK(myb_system_add_action(bare, 1, myb_system_light_action, &light) == MYB_OK);
  MT_CHECK(myb_system_acceleration(bare, 0.0, y, y + n, b) == MYB_E_DOMAIN);
  myb_system_destroy(bare);
  light.source = 1;
  MT_CHECK(myb_system_acceleration(s, 0.0, y, y + n, b) == MYB_E_DOMAIN);
  myb_system_destroy(s);
  myb_polyhedron_destroy(shape);
}

/* A body on a circle about the origin, as a provider: the context is its radius, its angular
 * rate and how many times it was asked. */
typedef struct circling {
  double radius;
  double rate;
  int calls;
  int fail;
} circling;

static myb_status circling_state(void* context, myb_time t, double state[6]) {
  circling* c = context;
  ++c->calls;
  if (c->fail) {
    return MYB_E_RANGE;
  }
  const double angle = c->rate * myb_time_seconds(t);
  state[0] = c->radius * cos(angle);
  state[1] = c->radius * sin(angle);
  state[2] = 0.0;
  state[3] = -c->radius * c->rate * sin(angle);
  state[4] = c->radius * c->rate * cos(angle);
  state[5] = 0.0;
  return MYB_OK;
}

/* Bodies with a given motion, a mass that changes, the dragging of the frame by a spin. */
static void test_guided_and_drags(void) {
  const double gm_sun = 1.32712440018e20;
  const double au = 1.495978707e11;
  /* a body of no mass about a guided centre at rest: Kepler's circle */
  {
    const double gm[1] = {0.0};
    myb_system* s = NULL;
    MT_CHECK(myb_system_create(1, gm, NULL, &s) == MYB_OK);
    circling centre = {0.0, 0.0, 0, 0};
    size_t index = 99;
    MT_CHECK(myb_system_add_guided(s, gm_sun, circling_state, &centre, &index) == MYB_OK);
    MT_CHECK(index == 0 && myb_system_guided(s) == 1);
    const double speed = sqrt(gm_sun / au);
    const double period = 2.0 * PI * au / speed;
    double y[6] = {au, 0.0, 0.0, 0.0, speed, 0.0};
    myb_ias15* in = NULL;
    MT_CHECK(myb_system_integrator(s, NULL, NULL, &in) == MYB_OK);
    MT_CHECK(myb_ias15_start(in, myb_system_acceleration, s, 0.0, y) == MYB_OK);
    MT_CHECK(myb_ias15_integrate(in, period) == MYB_OK);
    const double* z = myb_ias15_state(in);
    const long long steps = (long long)myb_ias15_statistics(in).accepted;
    const long long evaluations = (long long)myb_ias15_statistics(in).evaluations;
    printf("a body about a centre with a given motion: back after a period to %.1e of the "
           "radius; the provider asked %d times for %lld evaluations of %lld steps\n",
           sqrt((z[0] - au) * (z[0] - au) + z[1] * z[1]) / au, centre.calls, evaluations, steps);
    MT_CHECK(sqrt((z[0] - au) * (z[0] - au) + z[1] * z[1]) / au < 1e-11);
    /* the places are kept for the instants of a step: far fewer calls than evaluations */
    MT_CHECK(centre.calls < evaluations / 2 && centre.calls >= 7 * steps);
    /* the status of the provider comes back; a place on a body is refused */
    centre.fail = 1;
    double a[3];
    MT_CHECK(myb_system_acceleration(s, 12345.0, y, y + 3, a) == MYB_E_RANGE);
    centre.fail = 0;
    const double at_centre[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    MT_CHECK(myb_system_acceleration(s, 23456.0, at_centre, at_centre + 3, a) == MYB_E_DOMAIN);
    MT_CHECK(myb_system_add_guided(s, -1.0, circling_state, &centre, NULL) == MYB_E_DOMAIN);
    myb_ias15_destroy(in);
    myb_system_destroy(s);
  }
  /* a guided body pulls as an integrated one would: a planet among the bodies against the
   * same planet given its circle */
  {
    const double gm_planet = 1.26686534e17;
    const double far = 5.2 * au;
    const double rate = sqrt(gm_sun / (far * far * far));
    const double gm[3] = {gm_sun, 0.0, gm_planet};
    const double gm_alone[2] = {gm_sun, 0.0};
    myb_system* whole = NULL;
    myb_system* led = NULL;
    MT_CHECK(myb_system_create(3, gm, NULL, &whole) == MYB_OK);
    MT_CHECK(myb_system_create(2, gm_alone, NULL, &led) == MYB_OK);
    circling planet = {far, rate, 0, 0};
    MT_CHECK(myb_system_add_guided(led, gm_planet, circling_state, &planet, NULL) == MYB_OK);
    const double y[18] = {0.0, 0.0, 0.0, 0.7 * au, 0.3 * au, 0.1 * au, far, 0.0, 0.0,
                          0.0, 0.0, 0.0, -1.0e4, 2.0e4, 0.0,      0.0, far * rate, 0.0};
    const double y_led[12] = {0.0, 0.0, 0.0, 0.7 * au, 0.3 * au, 0.1 * au,
                              0.0, 0.0, 0.0, -1.0e4, 2.0e4,      0.0};
    double a[9];
    double b[6];
    MT_CHECK(myb_system_acceleration(whole, 0.0, y, y + 9, a) == MYB_OK);
    MT_CHECK(myb_system_acceleration(led, 0.0, y_led, y_led + 6, b) == MYB_OK);
    double worst = 0.0;
    for (int i = 0; i < 6; ++i) {
      worst = fmax(worst, fabs(a[i] - b[i]));
    }
    MT_CHECK(worst < 1e-18); /* of accelerations of 6e-3 and 2e-7 m/s^2 */
    myb_system_destroy(whole);
    myb_system_destroy(led);
  }
  /* a mass that changes, and the spin of the Sun */
  {
    const double gm[2] = {gm_sun, 2.2032e13};
    myb_system* s = NULL;
    MT_CHECK(myb_system_create(2, gm, NULL, &s) == MYB_OK);
    const double y[12] = {0.0, 0.0, 0.0, 0.31 * au, 0.2 * au, 0.05 * au,
                          0.0, 0.0, 0.0, -3.0e4,    4.0e4,    5.0e3};
    double before[6];
    double after[6];
    MT_CHECK(myb_system_acceleration(s, 0.0, y, y + 6, before) == MYB_OK);
    const double year = 365.25 * 86400.0;
    MT_CHECK(myb_system_set_mass_rate(s, 0, -9.0e-14 / year, 1000.0) == MYB_OK);
    MT_CHECK(myb_system_acceleration(s, 1000.0, y, y + 6, after) == MYB_OK);
    MT_CHECK(memcmp(before, after, sizeof before) == 0); /* at the epoch: the mass given */
    MT_CHECK(myb_system_acceleration(s, 1000.0 + 1.0e6 * year, y, y + 6, after) == MYB_OK);
    /* the planet is pulled by 9e-8 less after a million years; the Sun as before */
    MT_CHECK(fabs(after[3] / before[3] - (1.0 - 9.0e-8)) < 1e-15);
    MT_CHECK(fabs(after[0] / before[0] - 1.0) < 1e-15);
    MT_CHECK(myb_system_acceleration(s, 1000.0 + 2.0e13 * year, y, y + 6, after) == MYB_E_DOMAIN);
    MT_CHECK(myb_system_set_mass_rate(s, 0, 0.0, 0.0) == MYB_OK);
    MT_CHECK(myb_system_set_mass_rate(s, 5, 1.0, 0.0) == MYB_E_DOMAIN);
    /* the dragging of the frame needs the relativity of the point masses */
    const double spin[3] = {2.9e10, -8.2e10, 8.6e10};
    MT_CHECK(myb_system_add_frame_dragging(s, 0, 1, spin) == MYB_E_DOMAIN);
    const myb_ppn ppn = {299792458.0, 1.0, 1.0};
    MT_CHECK(myb_nbody_set_relativity(myb_system_nbody(s), &ppn) == MYB_OK);
    MT_CHECK(myb_system_acceleration(s, 0.0, y, y + 6, before) == MYB_OK);
    MT_CHECK(myb_system_add_frame_dragging(s, 0, 0, spin) == MYB_E_DOMAIN);
    MT_CHECK(myb_system_add_frame_dragging(s, 0, 1, spin) == MYB_OK);
    MT_CHECK(myb_system_acceleration(s, 0.0, y, y + 6, after) == MYB_OK);
    double expected[3];
    MT_CHECK(myb_rel_lense_thirring(&ppn, gm[0], spin, y + 3, y + 9, expected) == MYB_OK);
    double worst = 0.0;
    double size = 0.0;
    for (int i = 0; i < 3; ++i) {
      worst = fmax(worst, fabs(after[3 + i] - before[3 + i] - expected[i]));
      /* the momentum is kept: gm_sun a_sun + gm_planet a_planet does not change */
      worst = fmax(worst, fabs((after[i] - before[i]) * gm[0] / gm[1] + expected[i]));
      size = fmax(size, fabs(expected[i]));
    }
    printf("the spin of the Sun on a body at 0.37 au: %.2e m/s^2, %.1e of its attraction; in "
           "the system to %.1e m/s^2\n",
           size, size / (gm_sun / (0.1386 * au * au)), worst);
    /* the difference of two accelerations of 6e-3 m/s^2 resolves 1e-18 */
    MT_CHECK(size > 1e-14 && size < 1e-12 && worst < 1e-17);
    myb_system_destroy(s);
  }
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_field();
  test_tide();
  test_rotation();
  test_pn_pair();
  test_encounters();
  test_actions();
  test_guided_and_drags();
  return mt_finish("test_system");
}
