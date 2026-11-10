/* SPDX-License-Identifier: BSD-3-Clause */
/* Continuous thrust: the rocket equation in free space; the escape radius of a constant radial
 * acceleration (closed form); Gauss equations in modified equinoctial elements against Cowell;
 * Edelbaum's analytic relations and a numerical LEO transfer along his yaw law. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/integ.h"
#include "myb/kepler.h"
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

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_rocket_equation();
  test_radial_thrust();
  test_mee_vs_cowell();
  test_edelbaum();
  return mt_finish("test_thrust");
}
