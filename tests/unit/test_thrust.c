/* SPDX-License-Identifier: BSD-3-Clause */
/* Continuous thrust: the rocket equation in free space; the escape radius of a constant radial
 * acceleration (closed form); Gauss equations in modified equinoctial elements against Cowell;
 * Edelbaum's analytic relations and a numerical LEO transfer along his yaw law. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/integ.h"
#include "myb/env.h"
#include "myb/kepler.h"
#include "myb/model.h"
#include "myb/orbit.h"
#include "myb/thrust.h"

static const double MU = 3.986004418e14;
static const double PI = 3.141592653589793;

typedef struct cowell_ctx {
  double mu;
  myb_thrust* thrust;
} cowell_ctx;

static myb_status accel(void* context, double t, const double* x, const double* v, double* a) {
  const cowell_ctx* c = context;
  myb_status st = myb_thrust_acceleration(c->thrust, t, x, v, a);
  if (st != MYB_OK) {
    return st;
  }
  if (c->mu > 0.0) {
    const double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    const double k = c->mu / (r * r * r);
    for (int i = 0; i < 3; ++i) {
      a[i] -= k * x[i];
    }
  }
  return MYB_OK;
}

static myb_status cowell_rhs(void* context, double t, const double* y, double* dy) {
  double a[3];
  const myb_status st = accel(context, t, y, y + 3, a);
  for (int i = 0; i < 3; ++i) {
    dy[i] = y[3 + i];
    dy[3 + i] = a[i];
  }
  return st;
}

static void test_rocket_equation(void) {
  const myb_thrust_arc arc = {.start = 100.0, .stop = 5100.0, .thrust = 2.0, .exhaust = 30000.0,
                              .steering = MYB_STEER_INERTIAL, .direction = {0.0, 3.0, 4.0}};
  myb_thrust th;
  MT_CHECK(myb_thrust_init(&th, 50.0, &arc, 1) == MYB_OK);
  MT_CHECK(myb_thrust_mass(&th, 0.0) == 50.0);
  MT_CHECK(fabs(myb_thrust_mass(&th, 6000.0) - (50.0 - 2.0 / 30000.0 * 5000.0)) < 1e-14);
  cowell_ctx c = {.mu = 0.0, .thrust = &th};
  myb_dop853* d = NULL;
  const myb_dop853_options o = {.rtol = 1e-14, .atol = 1e-12, .h_max = 1000.0};
  MT_CHECK(myb_dop853_create(6, &o, NULL, &d) == MYB_OK);
  const double y0[6] = {0, 0, 0, 0, 0, 0};
  MT_CHECK(myb_dop853_start(d, cowell_rhs, &c, 0.0, y0) == MYB_OK);
  /* arc boundaries are discontinuities: integrate piecewise */
  MT_CHECK(myb_dop853_integrate(d, 100.0) == MYB_OK);
  MT_CHECK(myb_dop853_start(d, cowell_rhs, &c, 100.0, myb_dop853_state(d)) == MYB_OK);
  MT_CHECK(myb_dop853_integrate(d, 5100.0) == MYB_OK);
  const double* y = myb_dop853_state(d);
  const double dv = sqrt(y[3] * y[3] + y[4] * y[4] + y[5] * y[5]);
  const double expected = 30000.0 * log(50.0 / myb_thrust_mass(&th, 5100.0));
  printf("rocket equation: dv %.15g m/s, c ln(m0/m1) %.15g, relative %.2e; direction %.2e\n", dv,
         expected, fabs(dv - expected) / expected, fabs(y[3]));
  MT_CHECK(fabs(dv - expected) / expected < 1e-13);
  MT_CHECK(fabs(y[3]) < 1e-15 && fabs(y[4] / y[5] - 0.75) < 1e-14);
  myb_dop853_destroy(d);

  /* invalid arcs */
  myb_thrust_arc bad = arc;
  bad.direction[0] = bad.direction[1] = bad.direction[2] = 0.0;
  MT_CHECK(myb_thrust_init(&th, 50.0, &bad, 1) == MYB_E_DOMAIN);
  bad = arc;
  bad.thrust = 1e3; /* more propellant than mass */
  MT_CHECK(myb_thrust_init(&th, 50.0, &bad, 1) == MYB_E_DOMAIN);
  myb_thrust_arc two[2] = {arc, arc};
  MT_CHECK(myb_thrust_init(&th, 50.0, two, 2) == MYB_E_DOMAIN); /* overlapping */
}

static double radial_velocity(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[0] * y[3] + y[1] * y[4] + y[2] * y[5];
}

static void test_radial_thrust(void) {
  /* constant radial acceleration a from a circular orbit r0: bounded iff a <= mu / (8 r0^2),
   * the apocentre r = (mu - sqrt(mu^2 - 8 a mu r0^2)) / (4 a r0) (Tsien 1953; Prussing &
   * Coverstone-Carroll 1998); constant acceleration: a huge exhaust velocity */
  const double r0 = 7000e3;
  const double a = 0.9 * MU / (8.0 * r0 * r0);
  const double mass = 1000.0;
  const myb_thrust_arc arc = {.start = 0.0, .stop = 1e9, .thrust = a * mass, .exhaust = 1e30,
                              .steering = MYB_STEER_RTN, .direction = {1.0, 0.0, 0.0}};
  myb_thrust th;
  MT_CHECK(myb_thrust_init(&th, mass, &arc, 1) == MYB_OK);
  cowell_ctx c = {.mu = MU, .thrust = &th};
  const double vc = sqrt(MU / r0);
  const double y0[6] = {r0, 0, 0, 0, vc * 0.6, vc * 0.8};
  myb_dop853* d = NULL;
  const myb_dop853_options o = {.rtol = 1e-13, .atol = 1e-8};
  MT_CHECK(myb_dop853_create(6, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_dop853_start(d, cowell_rhs, &c, 0.0, y0) == MYB_OK);
  double y[6];
  int found = 0;
  double te = 0.0;
  for (int k = 0; k < 100000 && !found; ++k) {
    MT_CHECK(myb_dop853_step(d, 1e9) == MYB_OK);
    if (radial_velocity(NULL, 0.0, myb_dop853_state(d)) < 0.0) {
      MT_CHECK(myb_dop853_locate(d, radial_velocity, NULL, &found, &te, y) == MYB_OK);
    }
  }
  MT_CHECK(found);
  const double rmax = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
  const double expected = (MU - sqrt(MU * MU - 8.0 * a * MU * r0 * r0)) / (4.0 * a * r0);
  printf("constant radial acceleration 0.9 mu/(8 r0^2): apocentre %.6f km, closed form %.6f km, "
         "relative %.2e\n",
         rmax / 1e3, expected / 1e3, fabs(rmax - expected) / expected);
  MT_CHECK(fabs(rmax - expected) / expected < 1e-10);
  myb_dop853_destroy(d);
}

typedef struct mee_ctx {
  myb_thrust* thrust;
} mee_ctx;

static myb_status mee_rhs(void* context, double t, const double* y, double* dy) {
  const mee_ctx* c = context;
  const myb_equinoctial e = {y[0], y[1], y[2], y[3], y[4], y[5]};
  myb_pv pv;
  myb_status st = myb_pv_from_equinoctial(MU, &e, &pv);
  if (st != MYB_OK) {
    return st;
  }
  double artn[3];
  st = myb_thrust_rtn(c->thrust, t, pv.r, pv.v, artn);
  if (st != MYB_OK) {
    return st;
  }
  myb_equinoctial_rates(MU, &e, artn, dy);
  return MYB_OK;
}

static void test_mee_vs_cowell(void) {
  const myb_thrust_arc arc = {.start = 0.0, .stop = 10.0 * 86400.0, .thrust = 0.5,
                              .exhaust = 25000.0, .steering = MYB_STEER_RTN,
                              .direction = {0.2, 1.0, 0.3}};
  myb_thrust th;
  MT_CHECK(myb_thrust_init(&th, 500.0, &arc, 1) == MYB_OK);
  const myb_pv pv0 = {{6900e3, 300e3, 100e3}, {-200.0, 7000.0, 2500.0}};
  const double t1 = 10.0 * 86400.0;
  /* Cowell, IAS15 */
  cowell_ctx c = {.mu = MU, .thrust = &th};
  myb_ias15* ias = NULL;
  const myb_ias15_options io = {.epsilon = 1e-12};
  MT_CHECK(myb_ias15_create(3, &io, NULL, &ias) == MYB_OK);
  const double y0[6] = {pv0.r[0], pv0.r[1], pv0.r[2], pv0.v[0], pv0.v[1], pv0.v[2]};
  MT_CHECK(myb_ias15_start(ias, accel, &c, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(ias, t1) == MYB_OK);
  double xc[3];
  memcpy(xc, myb_ias15_state(ias), sizeof xc);
  myb_ias15_destroy(ias);
  /* Gauss equations in MEE, DOP853 */
  myb_equinoctial e0;
  MT_CHECK(myb_equinoctial_from_pv(MU, &pv0, &e0) == MYB_OK);
  mee_ctx m = {.thrust = &th};
  myb_dop853* d = NULL;
  const myb_dop853_options o = {.rtol = 1e-14, .atol = 1e-12};
  MT_CHECK(myb_dop853_create(6, &o, NULL, &d) == MYB_OK);
  const double ye[6] = {e0.p, e0.f, e0.g, e0.h, e0.k, e0.l};
  MT_CHECK(myb_dop853_start(d, mee_rhs, &m, 0.0, ye) == MYB_OK);
  MT_CHECK(myb_dop853_integrate(d, t1) == MYB_OK);
  const double* y = myb_dop853_state(d);
  const myb_equinoctial e1 = {y[0], y[1], y[2], y[3], y[4], y[5]};
  myb_pv pv1;
  MT_CHECK(myb_pv_from_equinoctial(MU, &e1, &pv1) == MYB_OK);
  const myb_dop853_stats ms = myb_dop853_statistics(d);
  myb_dop853_destroy(d);
  double diff = 0.0;
  for (int i = 0; i < 3; ++i) {
    diff = fmax(diff, fabs(pv1.r[i] - xc[i]));
  }
  const double rr = sqrt(xc[0] * xc[0] + xc[1] * xc[1] + xc[2] * xc[2]);
  printf("10 days of RTN thrust: Gauss/MEE (%lld f) vs Cowell IAS15: %.2e m of %.0f km\n",
         (long long)ms.evaluations, diff, rr / 1e3);
  MT_CHECK(diff / rr < 1e-9);
}

static double inclination(const double x[3], const double v[3]) {
  const double h[3] = {x[1] * v[2] - x[2] * v[1], x[2] * v[0] - x[0] * v[2],
                       x[0] * v[1] - x[1] * v[0]};
  return acos(h[2] / sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]));
}

static void test_edelbaum(void) {
  /* the analytic relations */
  double dv;
  double beta0;
  const double v0 = sqrt(MU / 7000e3);
  const double v1 = sqrt(MU / 7200e3);
  const double di = -1.0 * PI / 180.0;
  MT_CHECK(myb_edelbaum(v0, v1, fabs(di), &dv, &beta0) == MYB_OK);
  double v;
  double dii;
  myb_edelbaum_state(v0, beta0, dv, &v, &dii);
  printf("Edelbaum: 7000 -> 7200 km, 1 deg: dV = %.6f m/s, beta0 = %.6f deg; end state v %.2e, "
         "di %.2e\n",
         dv, beta0 * 180.0 / PI, fabs(v - v1) / v1, fabs(dii - fabs(di)));
  MT_CHECK(fabs(v - v1) / v1 < 1e-14 && fabs(dii - fabs(di)) < 1e-14);
  MT_CHECK(myb_edelbaum(0.0, v1, di, &dv, &beta0) == MYB_E_DOMAIN);
  /* no plane change: the Hohmann-less spiral dV = v0 - v1 */
  double dv0;
  double b0;
  MT_CHECK(myb_edelbaum(v0, v1, 0.0, &dv0, &b0) == MYB_OK);
  MT_CHECK(fabs(dv0 - (v0 - v1)) < 1e-11 && b0 == 0.0);

  /* numerical transfer with the yaw law: 0.2 N on 100 kg (2e-3 m/s^2, 2.4e-4 of gravity) */
  const double c = 30000.0;
  const double thrust = 0.2;
  const double m0 = 100.0;
  const double m1 = m0 * exp(-dv / c);
  const double duration = (m0 - m1) * c / thrust;
  const myb_thrust_arc arc = {.start = 0.0, .stop = duration, .thrust = thrust, .exhaust = c,
                              .steering = MYB_STEER_EDELBAUM, .edelbaum_v0 = v0,
                              .edelbaum_beta0 = beta0, .edelbaum_sense = -1.0};
  myb_thrust th;
  MT_CHECK(myb_thrust_init(&th, m0, &arc, 1) == MYB_OK);
  const double i0 = 28.5 * PI / 180.0;
  const myb_pv pv0 = {{7000e3, 0.0, 0.0}, {0.0, v0 * cos(i0), v0 * sin(i0)}};
  myb_equinoctial e0;
  MT_CHECK(myb_equinoctial_from_pv(MU, &pv0, &e0) == MYB_OK);
  mee_ctx m = {.thrust = &th};
  myb_dop853* d = NULL;
  const myb_dop853_options o = {.rtol = 1e-12, .atol = 1e-10};
  MT_CHECK(myb_dop853_create(6, &o, NULL, &d) == MYB_OK);
  const double ye[6] = {e0.p, e0.f, e0.g, e0.h, e0.k, e0.l};
  MT_CHECK(myb_dop853_start(d, mee_rhs, &m, 0.0, ye) == MYB_OK);
  MT_CHECK(myb_dop853_integrate(d, duration) == MYB_OK);
  const double* y = myb_dop853_state(d);
  const myb_equinoctial e1 = {y[0], y[1], y[2], y[3], y[4], y[5]};
  myb_pv pv1;
  MT_CHECK(myb_pv_from_equinoctial(MU, &e1, &pv1) == MYB_OK);
  myb_dop853_destroy(d);
  const double ecc = sqrt(e1.f * e1.f + e1.g * e1.g);
  const double a1 = e1.p / (1.0 - ecc * ecc);
  const double inc = inclination(pv1.r, pv1.v);
  printf("Edelbaum transfer (%.2f days): a %.3f km (7200), e %.2e, i %.5f deg (27.5)\n",
         duration / 86400.0, a1 / 1e3, ecc, inc * 180.0 / PI);
  /* the averaged theory is exact to the order of the thrust-to-gravity ratio (2.4e-4) */
  MT_CHECK(fabs(a1 - 7200e3) / 7200e3 < 2e-6); /* measured 1.4e-7 */
  MT_CHECK(fabs(inc * 180.0 / PI - 27.5) < 0.005); /* measured 0.0025 deg */
  MT_CHECK(ecc < 1e-3);
}

/* A transfer flown with the Q-law: the engine on until every targeted element is within its
 * tolerance. Returns the time taken; *dv the velocity increment spent, *lowest the lowest
 * periapsis radius met, *end the elements at the arrival. */
static double fly_qlaw(const myb_qlaw* law, const double tolerance[5], double effectivity,
                       const myb_elements* start, double budget, double* dv, double* lowest,
                       myb_elements* end) {
  myb_point_mass_term centre = {MU};
  myb_qlaw_term engine = {.law = law, .thrust = 1.0, .exhaust = 30000.0, .mass0 = 1000.0,
                          .effectivity = effectivity};
  memcpy(engine.tolerance, tolerance, sizeof engine.tolerance);
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "centre", myb_term_point_mass, myb_term_point_mass_jacobian,
                               &centre, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "engine", myb_term_qlaw, NULL, &engine, NULL) == MYB_OK);
  const myb_orbit_options options = {.integrator = MYB_ORBIT_DOP853, .rtol = 1e-10, .atol = 1e-4};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, &options, NULL, &orbit) == MYB_OK);
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(MU, start, &pv) == MYB_OK);
  double y[6];
  memcpy(y, pv.r, 3 * sizeof(double));
  memcpy(y + 3, pv.v, 3 * sizeof(double));
  MT_CHECK(myb_orbit_start(orbit, 0.0, y, budget) == MYB_OK);
  double t = 0.0;
  *lowest = 1e300;
  int arrived = 0;
  /* bound: the budget of time in steps of five minutes */
  while (!arrived && t < budget) {
    t = fmin(t + 300.0, budget);
    if (myb_orbit_state(orbit, t, y) != MYB_OK) {
      break;
    }
    memcpy(pv.r, y, 3 * sizeof(double));
    memcpy(pv.v, y + 3, 3 * sizeof(double));
    MT_CHECK(myb_elements_from_pv(MU, &pv, end) == MYB_OK);
    *lowest = fmin(*lowest, end->p / (1.0 + end->e));
    arrived = myb_qlaw_arrived(law, end, tolerance);
  }
  MT_CHECK(arrived);
  *dv = engine.exhaust * log(engine.mass0 / (engine.mass0 - engine.thrust / engine.exhaust * t));
  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);
  return t;
}

static void test_qlaw(void) {
  const double day = 86400.0;
  /* the direction of the law for the semi-major axis alone on a circle: along the motion */
  const myb_qlaw raise = {.mu = MU, .target = {1.0e7, 0.0, 0.0, 0.0, 0.0},
                          .weight = {1.0, 0.0, 0.0, 0.0, 0.0}};
  const myb_elements low = {.p = 7.0e6 * (1.0 - 1e-6), .e = 1e-3, .i = 0.3, .raan = 1.0,
                            .argp = 0.5, .nu = 2.0};
  double direction[3];
  double q = 0.0;
  double rate = 0.0;
  double effectivity = 0.0;
  MT_CHECK(myb_qlaw_direction(&raise, &low, direction, &q, &rate, &effectivity) == MYB_OK);
  /* mostly along the motion; a part along the radius, because a larger eccentricity would
   * let the semi-major axis be changed faster: the law sees that too */
  MT_CHECK(direction[1] > 0.98 && fabs(direction[2]) < 1e-12 && q > 0.0 && rate < 0.0);
  MT_CHECK(effectivity > 0.5 && effectivity <= 1.0);
  /* Q is of the order of the square of the time still needed: for this transfer at unit
   * acceleration, the difference of the circular speeds */
  const double v0 = sqrt(MU / 7.0e6);
  const double v1 = sqrt(MU / 1.0e7);
  printf("Q-law, 7000 -> 10 000 km: sqrt(Q) = %.1f s per m/s^2 against the difference of the "
         "circular speeds %.1f m/s\n",
         sqrt(q), v0 - v1);
  MT_CHECK(fabs(sqrt(q) / (v0 - v1) - 1.0) < 0.4);
  /* flown: against Edelbaum's velocity increment */
  const double loose[5] = {2.0e3, 1e-2, 2e-4, 1e-3, 1e-3};
  double dv = 0.0;
  double lowest = 0.0;
  myb_elements end;
  const double time = fly_qlaw(&raise, loose, 0.0, &low, 40.0 * day, &dv, &lowest, &end);
  printf("  flown in %.2f days: %.1f m/s (Edelbaum %.1f); arrives with e = %.4f\n", time / day,
         dv, v0 - v1, end.e);
  MT_CHECK(fabs(dv / (v0 - v1) - 1.0) < 0.02);
  MT_CHECK(end.e < 0.01);
  /* with a change of the plane of 0.1 rad: within a few per cent of Edelbaum's optimum, and
   * not below it by more than the law's bang at the end allows */
  const myb_qlaw tilt = {.mu = MU, .target = {1.0e7, 0.0, 0.2, 0.0, 0.0},
                         .weight = {1.0, 0.0, 1.0, 0.0, 0.0}};
  double optimum = 0.0;
  double beta0 = 0.0;
  MT_CHECK(myb_edelbaum(v0, v1, 0.1, &optimum, &beta0) == MYB_OK);
  const double time_tilt = fly_qlaw(&tilt, loose, 0.0, &low, 60.0 * day, &dv, &lowest, &end);
  printf("  with the plane turned by 0.1 rad: %.2f days, %.1f m/s (Edelbaum %.1f, ratio "
         "%.3f); i = %.5f\n",
         time_tilt / day, dv, optimum, dv / optimum, end.i);
  MT_CHECK(dv / optimum > 0.97 && dv / optimum < 1.15);
  MT_CHECK(fabs(end.i - 0.2) < 3e-4);
  /* from a transfer ellipse to a high circle, the periapsis held above 6600 km */
  const myb_qlaw circularise = {.mu = MU, .target = {4.2164e7, 0.0, 0.05, 0.0, 0.0},
                                .weight = {1.0, 1.0, 1.0, 0.0, 0.0}, .rp_min = 6.6e6,
                                .rp_weight = 1.0};
  const myb_elements ellipse = {.p = 2.45e7 * (1.0 - 0.725 * 0.725), .e = 0.725, .i = 0.12,
                                .raan = 0.3, .argp = 3.0, .nu = 0.0};
  const double tight[5] = {5.0e4, 5e-3, 5e-4, 1e-3, 1e-3};
  const double time_geo =
      fly_qlaw(&circularise, tight, 0.3, &ellipse, 150.0 * day, &dv, &lowest, &end);
  /* with a threshold of 0.3 on the effectivity: without it the engine dithers 100 km from
   * the target and advances no more. The velocity increment is that of an engine that
   * never paused: more than was spent */
  printf("  transfer ellipse (e = 0.725) to a circle at 42 164 km: %.1f days, no more than "
         "%.0f m/s; the periapsis never below %.0f km; arrives with a = %.0f km, e = %.4f, "
         "i = %.4f\n",
         time_geo / day, dv, lowest / 1e3, end.p / (1.0 - end.e * end.e) / 1e3, end.e, end.i);
  MT_CHECK(lowest > 6.6e6);
  MT_CHECK(time_geo < 40.0 * day && dv < 3500.0);
  MT_CHECK(fabs(end.p / (1.0 - end.e * end.e) - 4.2164e7) <= 5.0e4 && end.e <= 5e-3 &&
           fabs(end.i - 0.05) <= 5e-4);
  /* the same flight with the mass in the state: the pauses cost nothing */
  const myb_qlaw_term engine = {.law = &circularise, .thrust = 1.0, .exhaust = 30000.0,
                                .tolerance = {5.0e4, 5e-3, 5e-4, 1e-3, 1e-3},
                                .effectivity = 0.3};
  const myb_qlaw_flight flight = {.engine = &engine};
  myb_pv from;
  MT_CHECK(myb_pv_from_elements(MU, &ellipse, &from) == MYB_OK);
  double y0[7];
  memcpy(y0, from.r, 3 * sizeof(double));
  memcpy(y0 + 3, from.v, 3 * sizeof(double));
  y0[6] = 1000.0;
  myb_qlaw_arrival arrival;
  MT_CHECK(myb_qlaw_fly(&flight, 0.0, y0, 150.0 * day, 300.0, NULL, &arrival) == MYB_OK);
  const double spent = 30000.0 * log(1000.0 / arrival.y[6]);
  const double burning = (1000.0 - arrival.y[6]) * 30000.0 / 1.0;
  printf("  the same with the mass in the state: %.1f days, the engine on for %.1f of them, "
         "%.0f m/s spent (%.0f if it had never paused), %.2f kg left of 1000\n",
         arrival.t / day, burning / day, spent, dv, arrival.y[6]);
  MT_CHECK(arrival.arrived);
  MT_CHECK(fabs(arrival.t - time_geo) < 0.5 * day);
  MT_CHECK(spent < 0.95 * dv && spent > 1500.0);
  MT_CHECK(burning < 0.95 * arrival.t);
  /* an engine that never pauses spends thrust / exhaust all the way: the rocket equation */
  const myb_qlaw_term steady = {.law = &raise, .thrust = 1.0, .exhaust = 30000.0,
                                .tolerance = {2.0e3, 1e-2, 2e-4, 1e-3, 1e-3}};
  myb_point_mass_term centre = {MU};
  myb_force_model* field = NULL;
  MT_CHECK(myb_force_model_create(NULL, &field) == MYB_OK);
  MT_CHECK(myb_force_model_add(field, "centre", myb_term_point_mass, myb_term_point_mass_jacobian,
                               &centre, NULL) == MYB_OK);
  const myb_qlaw_flight plain = {.engine = &steady};
  const myb_qlaw_flight modelled = {.engine = &steady, .field = field};
  MT_CHECK(myb_pv_from_elements(MU, &low, &from) == MYB_OK);
  memcpy(y0, from.r, 3 * sizeof(double));
  memcpy(y0 + 3, from.v, 3 * sizeof(double));
  myb_qlaw_arrival a1;
  myb_qlaw_arrival a2;
  MT_CHECK(myb_qlaw_fly(&plain, 0.0, y0, 40.0 * day, 300.0, NULL, &a1) == MYB_OK);
  MT_CHECK(myb_qlaw_fly(&modelled, 0.0, y0, 40.0 * day, 300.0, NULL, &a2) == MYB_OK);
  printf("  the rise to 10 000 km with the mass in the state: %.2f days, %.6f kg left "
         "(%.6f had it burned to the last look); with the field as a model: %.1e m apart\n",
         a1.t / day, a1.y[6], 1000.0 - a1.t / 30000.0,
         sqrt((a1.y[0] - a2.y[0]) * (a1.y[0] - a2.y[0]) + (a1.y[1] - a2.y[1]) * (a1.y[1] - a2.y[1]) +
              (a1.y[2] - a2.y[2]) * (a1.y[2] - a2.y[2])));
  MT_CHECK(a1.arrived && a2.arrived && a1.t == a2.t);
  MT_CHECK(fabs(a1.t - time) <= 300.0);
  /* the engine stops at the arrival, which is seen at the next look: up to 300 s later */
  MT_CHECK(a1.y[6] - (1000.0 - a1.t / 30000.0) >= -1e-6 &&
           a1.y[6] - (1000.0 - a1.t / 30000.0) <= 300.0 / 30000.0 + 1e-6);
  MT_CHECK(fabs(a1.y[0] - a2.y[0]) + fabs(a1.y[1] - a2.y[1]) + fabs(a1.y[2] - a2.y[2]) < 1.0);
  myb_force_model_destroy(field);
  /* a limit of time before the arrival; refusals */
  MT_CHECK(myb_qlaw_fly(&plain, 0.0, y0, 1.0 * day, 300.0, NULL, &a1) == MYB_OK);
  MT_CHECK(!a1.arrived && a1.t == 1.0 * day);
  MT_CHECK(myb_qlaw_fly(&plain, 0.0, y0, 0.0, 300.0, NULL, &a1) == MYB_E_DOMAIN);
  MT_CHECK(myb_qlaw_fly(&plain, 0.0, y0, 1.0 * day, 0.0, NULL, &a1) == MYB_E_DOMAIN);
  y0[6] = 0.0;
  MT_CHECK(myb_qlaw_fly(&plain, 0.0, y0, 1.0 * day, 300.0, NULL, &a1) == MYB_E_DOMAIN);
  /* at the target nothing is asked of the engine */
  const myb_elements there = {.p = 1.0e7, .e = 0.0, .i = 0.3, .raan = 1.0, .argp = 0.0, .nu = 1.0};
  MT_CHECK(myb_qlaw_arrived(&raise, &there, loose));
  /* refused */
  const myb_qlaw nothing = {.mu = MU, .target = {1.0e7, 0.0, 0.0, 0.0, 0.0}};
  MT_CHECK(myb_qlaw_direction(&nothing, &low, direction, NULL, NULL, NULL) == MYB_E_DOMAIN);
  const myb_elements open = {.p = 2.0e7, .e = 1.2, .i = 0.3, .raan = 1.0, .argp = 0.5, .nu = 0.1};
  MT_CHECK(myb_qlaw_direction(&raise, &open, direction, NULL, NULL, NULL) == MYB_E_DOMAIN);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_rocket_equation();
  test_radial_thrust();
  test_mee_vs_cowell();
  test_edelbaum();
  test_qlaw();
  return mt_finish("test_thrust");
}
