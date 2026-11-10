/* SPDX-License-Identifier: BSD-3-Clause */
/* Small forces of §4 (myb/force.h, myb/env.h, myb/orbit.h) against what the papers derive
 * from them:
 *  - light on a moving grain: the Poynting-Robertson decay of a and e (Wyatt & Whipple 1950;
 *    Burns, Lamy & Soter 1979, eq. 47) over twenty revolutions; the Jacobian of the term;
 *  - the accelerations of comets and asteroids: g(r) of Marsden, Sekanina & Yeomans; the
 *    Yarkovsky drift da/dt = 2 A2 (1 - e^2) (r0 / p)^2 / n (Farnocchia et al. 2013); the drift
 *    under g(r) of a comet against Gauss's equation averaged over the orbit;
 *  - the Yarkovsky drift of the linear theory: its limits, signs and scaling;
 *  - the lagging tide: da/dt = 6 k2 lag (m / M) (R / a)^5 n a (omega - n) (Mignard 1979), and
 *    3.8 cm a year for the Moon;
 *  - the recoil of emitted power; a mass that changes; impulses and the transition matrix
 *    across them (a Hohmann transfer). */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/myb.h"

static const double pi = 3.141592653589793;
static const double two_pi = 6.283185307179586;
static const double au = 1.495978707e11;
static const double gm_sun = 1.32712440018e20;
static const double c_light = 299792458.0;

static double norm3(const double v[3]) {
  return hypot(hypot(v[0], v[1]), v[2]);
}

typedef struct spinner {
  double rate;
} spinner;

static myb_status spinner_orientation(void* context, myb_time t, myb_orientation* o) {
  const spinner* s = context;
  const double w = s->rate * myb_time_seconds(t);
  const double c = myb_cos(w);
  const double sn = myb_sin(w);
  const double m[9] = {c, sn, 0.0, -sn, c, 0.0, 0.0, 0.0, 1.0};
  const double dm[9] = {-sn * s->rate, c * s->rate, 0.0, -c * s->rate, -sn * s->rate, 0.0,
                        0.0,           0.0,         0.0};
  memcpy(o->rotation.m, m, sizeof m);
  memcpy(o->rotation.dm, dm, sizeof dm);
  o->ut = t;
  o->xp = 0.0;
  o->yp = 0.0;
  return MYB_OK;
}

static myb_status fixed_state(void* context, myb_time t, double y[6]) {
  (void)t;
  memcpy(y, context, 6 * sizeof(double));
  return MYB_OK;
}

/* a and e of the state y about mu. */
static void elements(double mu, const double y[6], double* a, double* e) {
  myb_pv pv;
  memcpy(pv.r, y, 3 * sizeof(double));
  memcpy(pv.v, y + 3, 3 * sizeof(double));
  myb_elements el;
  MT_CHECK(myb_elements_from_pv(mu, &pv, &el) == MYB_OK);
  *e = el.e;
  *a = el.p / (1.0 - el.e * el.e);
}

static void state_of(double mu, double a, double e, double inclination, double y[6]) {
  const myb_elements el = {.p = a * (1.0 - e * e), .e = e, .i = inclination, .raan = 0.3,
                           .argp = 0.7, .nu = 0.0};
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(mu, &el, &pv) == MYB_OK);
  memcpy(y, pv.r, 3 * sizeof(double));
  memcpy(y + 3, pv.v, 3 * sizeof(double));
}

/* The state at the n-th pericentre after the start (the start is one). */
static double pericentre(myb_force_model* model, const double y0[6], double span, int n,
                         double y[6]) {
  const myb_orbit_options options = {.epsilon = 1e-11};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, &options, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, span) == MYB_OK);
  myb_range_event centre = {.body = MYB_EVENT_ORIGIN};
  /* away from the pericentre of the start before looking for the next one */
  MT_CHECK(myb_orbit_state(orbit, span / (4.4 * n), y) == MYB_OK);
  double t = 0.0;
  for (int k = 0; k < n; ++k) {
    int found = 0;
    MT_CHECK(myb_orbit_event(orbit, myb_event_range_rate, &centre, span, +1, &found, &t, y) ==
             MYB_OK);
    MT_CHECK(found);
  }
  myb_orbit_destroy(orbit);
  return t;
}

/* ---- Poynting-Robertson ----------------------------------------------------------------------- */

static void test_poynting_robertson(void) {
  /* the function: at rest the force of a sphere; moving away, less by twice v / c */
  const double u[3] = {0.6, 0.0, 0.8};
  const double rest[3] = {0.0, 0.0, 0.0};
  double a[3] = {0.0, 0.0, 0.0};
  double ref[3] = {0.0, 0.0, 0.0};
  myb_force_radiation_moving(4.5e-6, 1.3, 0.02, u, rest, c_light, a);
  myb_force_radiation_sphere(4.5e-6, 1.3, 0.02, u, ref);
  MT_CHECK(a[0] == ref[0] && a[1] == ref[1] && a[2] == ref[2]);
  const double away[3] = {0.6 * 3.0e4, 0.0, 0.8 * 3.0e4};
  double moving[3] = {0.0, 0.0, 0.0};
  myb_force_radiation_moving(4.5e-6, 1.3, 0.02, u, away, c_light, moving);
  MT_CHECK(fabs(norm3(moving) / norm3(ref) - (1.0 - 2.0 * 3.0e4 / c_light)) < 1e-15);
  const double across[3] = {0.8 * 3.0e4, 0.0, -0.6 * 3.0e4};
  double drag[3] = {0.0, 0.0, 0.0};
  myb_force_radiation_moving(4.5e-6, 1.3, 0.02, u, across, c_light, drag);
  const double along = (drag[0] * across[0] + drag[2] * across[2]) / 3.0e4;
  MT_CHECK(fabs(along / norm3(ref) + 3.0e4 / c_light) < 1e-15); /* the drag: -v / c */

  /* a grain about the Sun with beta = 0.2 */
  const double beta = 0.2;
  const double luminosity = 3.828e26;
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body sun = {.id = 10, .gm = gm_sun, .luminosity = luminosity};
  size_t i_sun = 0;
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  /* beta = L cr (A / m) / (4 pi c mu) */
  myb_object grain = {.mass = 1.0, .cr = 1.0};
  grain.area = beta * 4.0 * pi * c_light * gm_sun / luminosity;
  myb_point_mass_term central = {.gm = gm_sun};
  myb_radiation_term light = {.env = env, .source = i_sun, .object = &grain, .c = c_light,
                              .moving = 1};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "sun", myb_term_point_mass, myb_term_point_mass_jacobian,
                               &central, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "light", myb_term_radiation, myb_term_radiation_jacobian,
                               &light, NULL) == MYB_OK);
  const double mu = gm_sun * (1.0 - beta); /* what the grain feels: the elements are about it */
  const double a0 = au;
  const double e0 = 0.3;
  double y0[6];
  state_of(mu, a0, e0, 0.4, y0);
  const int revolutions = 20;
  const double period = two_pi * sqrt(a0 * a0 * a0 / mu);
  double y[6];
  const double t = pericentre(model, y0, 1.1 * revolutions * period, revolutions, y);
  double a1 = 0.0;
  double e1 = 0.0;
  elements(mu, y, &a1, &e1);
  /* the averaged equations, integrated over the same time (Runge-Kutta 4, 2000 steps) */
  const double alpha = beta * gm_sun / c_light;
  double as = a0;
  double es = e0;
  const int steps = 2000;
  const double h = t / steps;
  for (int k = 0; k < steps; ++k) {
    double ka[4];
    double ke[4];
    double aa = as;
    double ee = es;
    for (int stage = 0; stage < 4; ++stage) {
      const double q = 1.0 - ee * ee;
      ka[stage] = -alpha * (2.0 + 3.0 * ee * ee) / (aa * q * sqrt(q));
      ke[stage] = -2.5 * alpha * ee / (aa * aa * sqrt(q));
      const double w = stage == 2 ? 1.0 : 0.5;
      aa = as + w * h * ka[stage];
      ee = es + w * h * ke[stage];
    }
    as += h * (ka[0] + 2.0 * ka[1] + 2.0 * ka[2] + ka[3]) / 6.0;
    es += h * (ke[0] + 2.0 * ke[1] + 2.0 * ke[2] + ke[3]) / 6.0;
  }
  printf("Poynting-Robertson, beta 0.2, %d revolutions: a falls by %.6e au (averaged theory "
         "%.6e), e by %.6e (%.6e)\n",
         revolutions, (a0 - a1) / au, (a0 - as) / au, e0 - e1, e0 - es);
  MT_CHECK(fabs((a0 - a1) / (a0 - as) - 1.0) < 2e-4);
  MT_CHECK(fabs((e0 - e1) / (e0 - es) - 1.0) < 2e-4);

  /* the Jacobian of the moving term against differences */
  double dadx[9];
  double dadv[9];
  MT_CHECK(myb_term_radiation_jacobian(&light, 0.0, y0, y0 + 3, dadx, dadv) == MYB_OK);
  double worst_x = 0.0;
  double worst_v = 0.0;
  double size_x = 0.0;
  double size_v = 0.0;
  for (int j = 0; j < 3; ++j) {
    double plus[3];
    double minus[3];
    double q[6];
    memcpy(q, y0, sizeof q);
    const double hx = 1e-5 * au;
    q[j] = y0[j] + hx;
    MT_CHECK(myb_term_radiation(&light, 0.0, q, q + 3, plus) == MYB_OK);
    q[j] = y0[j] - hx;
    MT_CHECK(myb_term_radiation(&light, 0.0, q, q + 3, minus) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      worst_x = fmax(worst_x, fabs((plus[i] - minus[i]) / (2.0 * hx) - dadx[3 * i + j]));
      size_x = fmax(size_x, fabs(dadx[3 * i + j]));
    }
    memcpy(q, y0, sizeof q);
    const double hv = 10.0;
    q[3 + j] = y0[3 + j] + hv;
    MT_CHECK(myb_term_radiation(&light, 0.0, q, q + 3, plus) == MYB_OK);
    q[3 + j] = y0[3 + j] - hv;
    MT_CHECK(myb_term_radiation(&light, 0.0, q, q + 3, minus) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      worst_v = fmax(worst_v, fabs((plus[i] - minus[i]) / (2.0 * hv) - dadv[3 * i + j]));
      size_v = fmax(size_v, fabs(dadv[3 * i + j]));
    }
  }
  printf("  Jacobian of the moving term: position %.1e, velocity %.1e of the largest element\n",
         worst_x / size_x, worst_v / size_v);
  MT_CHECK(worst_x < 1e-8 * size_x && worst_v < 1e-8 * size_v);
  myb_force_model_destroy(model);
  myb_env_destroy(env);
}

/* ---- Comets and asteroids --------------------------------------------------------------------- */

/* da/dt of Gauss's equation averaged over the orbit for radial and transverse accelerations
 * A1 g(r) and A2 g(r). */
static double gauss_drift(const myb_nongrav* model, double mu, double a, double e) {
  const double n = sqrt(mu / (a * a * a));
  const double p = a * (1.0 - e * e);
  const int points = 4000;
  double sum = 0.0;
  for (int k = 0; k < points; ++k) {
    const double ecc = two_pi * (k + 0.5) / points;
    const double r = a * (1.0 - e * cos(ecc));
    const double sin_f = sqrt(1.0 - e * e) * sin(ecc) * a / r;
    const double g = myb_nongrav_g(model, r);
    /* dM = (r / a) dE */
    sum += (e * sin_f * model->a[0] * g + (p / r) * model->a[1] * g) * (r / a);
  }
  return 2.0 / (n * sqrt(1.0 - e * e)) * sum / points;
}

static void test_small_bodies(void) {
  myb_nongrav comet;
  myb_nongrav_marsden(1.0e-8, 2.0e-9, -3.0e-10, &comet);
  printf("g of Marsden, Sekanina and Yeomans: g(1 au) - 1 = %.1e\n",
         myb_nongrav_g(&comet, au) - 1.0);
  MT_CHECK(fabs(myb_nongrav_g(&comet, au) - 1.0) < 1e-8); /* alpha has ten digits */
  /* at r0: alpha 2^-k; far beyond it the fall steepens to r^-(m + n k) */
  MT_CHECK(fabs(myb_nongrav_g(&comet, 2.808 * au) / (0.1112620426 * pow(2.0, -4.6142)) - 1.0) <
           1e-14);
  const double far = myb_nongrav_g(&comet, 40.0 * au) / myb_nongrav_g(&comet, 20.0 * au);
  MT_CHECK(fabs(log2(far) + 2.15 + 5.093 * 4.6142) < 1e-3);
  MT_CHECK(myb_nongrav_g(&comet, 0.0) != myb_nongrav_g(&comet, 0.0)); /* NaN */
  myb_nongrav inverse_square;
  myb_nongrav_power(0.0, 1.0e-10, 0.0, au, 2.0, &inverse_square);
  MT_CHECK(fabs(myb_nongrav_g(&inverse_square, 2.0 * au) - 0.25) < 1e-16);
  /* the directions */
  const double r[3] = {1.2 * au, 0.3 * au, -0.2 * au};
  const double v[3] = {-4.0e3, 2.6e4, 1.0e3};
  double a[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_nongrav(&comet, r, v, a) == MYB_OK);
  const double rn = norm3(r);
  const double h[3] = {r[1] * v[2] - r[2] * v[1], r[2] * v[0] - r[0] * v[2],
                       r[0] * v[1] - r[1] * v[0]};
  const double hn = norm3(h);
  const double g = myb_nongrav_g(&comet, rn);
  const double radial = (a[0] * r[0] + a[1] * r[1] + a[2] * r[2]) / rn;
  const double normal = (a[0] * h[0] + a[1] * h[1] + a[2] * h[2]) / hn;
  const double transverse = sqrt(fmax(0.0, norm3(a) * norm3(a) - radial * radial - normal * normal));
  MT_CHECK(fabs(radial / (1.0e-8 * g) - 1.0) < 1e-13);
  MT_CHECK(fabs(normal / (-3.0e-10 * g) - 1.0) < 1e-12);
  MT_CHECK(fabs(transverse / (2.0e-9 * g) - 1.0) < 1e-12);
  /* towards the motion */
  MT_CHECK(a[0] * v[0] + a[1] * v[1] + a[2] * v[2] - radial * (r[0] * v[0] + r[1] * v[1] +
                                                                  r[2] * v[2]) / rn > 0.0);
  const double zero[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_nongrav(&comet, zero, v, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_nongrav(&comet, r, r, a) == MYB_E_DOMAIN); /* no orbital plane */

  /* the drifts, through the term of a model: the Sun at the origin */
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body sun = {.id = 10, .gm = gm_sun};
  size_t i_sun = 0;
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  myb_point_mass_term central = {.gm = gm_sun};
  myb_nongrav_term term = {.env = env, .source = i_sun, .model = inverse_square};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "sun", myb_term_point_mass, myb_term_point_mass_jacobian,
                               &central, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "nongrav", myb_term_nongrav, NULL, &term, NULL) == MYB_OK);
  const double a0 = 1.3 * au;
  const double e0 = 0.3;
  double y0[6];
  state_of(gm_sun, a0, e0, 0.2, y0);
  const double period = two_pi * sqrt(a0 * a0 * a0 / gm_sun);
  const int revolutions = 10;
  double y[6];
  double t = pericentre(model, y0, 1.1 * revolutions * period, revolutions, y);
  double a1 = 0.0;
  double e1 = 0.0;
  elements(gm_sun, y, &a1, &e1);
  const double n = two_pi / period;
  const double p = a0 * (1.0 - e0 * e0);
  const double yarkovsky = 2.0 * 1.0e-10 * (1.0 - e0 * e0) * (au / p) * (au / p) / n;
  printf("Yarkovsky, A2 = 1e-10 m/s^2: da/dt %.9e m/s over %d revolutions, the formula %.9e\n",
         (a1 - a0) / t, revolutions, yarkovsky);
  MT_CHECK(fabs((a1 - a0) / t / yarkovsky - 1.0) < 1e-6);
  MT_CHECK(fabs(gauss_drift(&inverse_square, gm_sun, a0, e0) / yarkovsky - 1.0) < 1e-12);
  /* a comet: A1 and A2 with the g of water ice, against the averaged Gauss equation */
  myb_nongrav_marsden(3.0e-9, 4.0e-10, 0.0, &term.model);
  t = pericentre(model, y0, 1.1 * revolutions * period, revolutions, y);
  elements(gm_sun, y, &a1, &e1);
  const double expected = gauss_drift(&term.model, gm_sun, a0, e0);
  printf("comet, A1 = 3e-9, A2 = 4e-10 m/s^2: da/dt %.9e m/s, Gauss's equation averaged "
         "%.9e\n",
         (a1 - a0) / t, expected);
  MT_CHECK(fabs((a1 - a0) / t / expected - 1.0) < 1e-5);
  myb_force_model_destroy(model);
  myb_env_destroy(env);
}

static void test_yarkovsky_theory(void) {
  /* a kilometre of rock at 2.5 au, six hours of rotation */
  myb_yarkovsky body = {.radius = 500.0, .density = 2500.0, .absorptivity = 0.9,
                        .emissivity = 0.9, .thermal_inertia = 200.0,
                        .rotation_rate = two_pi / 21600.0, .obliquity = 0.0};
  const double a = 2.5 * au;
  const double flux = 1361.0 / 6.25;
  const double n = sqrt(gm_sun / (a * a * a));
  double diurnal = 0.0;
  double seasonal = 0.0;
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &diurnal, &seasonal) == MYB_OK);
  const double per_myr = 3.15576e13 / au; /* m/s -> au per million years */
  printf("Yarkovsky, a kilometre of rock at 2.5 au, obliquity 0: %.2e au/Myr (diurnal), %.1e "
         "(seasonal)\n",
         diurnal * per_myr, seasonal * per_myr);
  /* outwards for a prograde rotation, of the order of 1e-4 au per million years */
  MT_CHECK(diurnal > 0.0 && diurnal * per_myr > 2e-5 && diurnal * per_myr < 5e-4);
  MT_CHECK(seasonal == 0.0);
  /* against the formula written out */
  const double sigma = 5.670374419e-8;
  const double phi = 3.0 * flux / (4.0 * 500.0 * 2500.0 * c_light);
  const double t_star = pow(0.9 * flux / (0.9 * sigma), 0.25);
  const double theta = 200.0 * sqrt(body.rotation_rate) / (0.9 * sigma * t_star * t_star * t_star);
  const double f = -0.5 * theta / (1.0 + theta + 0.5 * theta * theta);
  MT_CHECK(fabs(diurnal / (-(8.0 / 9.0) * 0.9 * phi / n * f) - 1.0) < 1e-13);
  /* retrograde: inwards, the same size; on its side: the seasonal effect alone, inwards */
  double retro = 0.0;
  body.obliquity = pi;
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &retro, &seasonal) == MYB_OK);
  MT_CHECK(fabs(retro + diurnal) < 1e-15 * diurnal && fabs(seasonal) < 1e-30);
  body.obliquity = pi / 2.0;
  double side = 0.0;
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &side, &seasonal) == MYB_OK);
  MT_CHECK(fabs(side) < 1e-15 * diurnal && seasonal < 0.0);
  /* twice the radius, half the drift; no thermal inertia or no rotation, none */
  body.obliquity = 0.0;
  body.radius = 1000.0;
  double larger = 0.0;
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &larger, NULL) == MYB_OK);
  MT_CHECK(fabs(larger / diurnal - 0.5) < 1e-15);
  body.thermal_inertia = 0.0;
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &larger, NULL) == MYB_OK);
  MT_CHECK(larger == 0.0);
  /* the largest effect: thermal parameter sqrt(2) */
  body.thermal_inertia = sqrt(2.0) * 0.9 * sigma * t_star * t_star * t_star /
                         sqrt(body.rotation_rate);
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &larger, NULL) == MYB_OK);
  const double best = (8.0 / 9.0) * 0.9 * (phi / 2.0) / n * (sqrt(2.0) - 1.0) / 2.0;
  MT_CHECK(fabs(larger / best - 1.0) < 1e-13);
  body.radius = -1.0;
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &larger, NULL) == MYB_E_DOMAIN);
}

/* ---- The lagging tide ------------------------------------------------------------------------- */

static void test_tidal_lag(void) {
  const double gm_earth = 3.986004418e14;
  const double gm_moon = 4.9028e12;
  const double radius = 6378136.6;
  const double a0 = 3.844e8;
  const double mu = gm_earth + gm_moon;
  const double n = sqrt(mu / (a0 * a0 * a0));
  const double spin = 7.292115e-5;
  /* the Moon: k2 0.3, ten minutes of delay */
  const double k2 = 0.3;
  const double rate = 6.0 * k2 * 600.0 * (gm_moon / gm_earth) * pow(radius / a0, 5.0) * n * a0 *
                      (spin - n);
  printf("the lagging tide of the Moon on the Earth (k2 0.3, 600 s): da/dt = %.2f cm a year\n",
         rate * 3.15576e7 * 100.0);
  MT_CHECK(rate * 3.15576e9 > 3.7 && rate * 3.15576e9 < 3.9);
  /* the function: across the radius, forwards, for a spin faster than the orbit */
  const double omega[3] = {0.0, 0.0, spin};
  const double r[3] = {a0, 0.0, 0.0};
  const double v[3] = {0.0, n * a0, 0.0};
  double a[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_tidal_lag(gm_earth, gm_moon, k2, radius, 600.0, omega, r, v, 0, a) ==
           MYB_OK);
  MT_CHECK(a[0] == 0.0 && a[2] == 0.0 && a[1] > 0.0);
  MT_CHECK(fabs(2.0 * a[1] / n / rate - 1.0) < 1e-14); /* da/dt = 2 T / n */
  /* a synchronous orbit: nothing; with the tide that does not lag: along the radius */
  const double sync[3] = {0.0, spin * a0, 0.0};
  double none[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_tidal_lag(gm_earth, gm_moon, k2, radius, 600.0, omega, r, sync, 0, none) ==
           MYB_OK);
  MT_CHECK(none[0] == 0.0 && none[1] == 0.0 && none[2] == 0.0);
  double whole[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_tidal_lag(gm_earth, gm_moon, k2, radius, 600.0, omega, r, sync, 1, whole) ==
           MYB_OK);
  const double tide = 3.0 * k2 * gm_moon * (1.0 + gm_moon / gm_earth) * pow(radius, 5.0) /
                      pow(a0, 7.0);
  MT_CHECK(fabs(whole[0] / -tide - 1.0) < 1e-14 && whole[1] == 0.0);
  MT_CHECK(myb_force_tidal_lag(gm_earth, gm_moon, k2, radius, 600.0, omega, none, v, 0, a) ==
           MYB_E_DOMAIN);

  /* the orbit: a delay a thousand times longer, so that ten revolutions show the drift */
  spinner earth_spin = {.rate = spin};
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body earth = {.id = 399, .gm = gm_earth, .radius = radius,
                              .orientation = spinner_orientation,
                              .orientation_context = &earth_spin};
  size_t i_earth = 0;
  MT_CHECK(myb_env_add_body(env, &earth, &i_earth) == MYB_OK);
  myb_point_mass_term central = {.gm = mu};
  myb_tidal_lag_term term = {.env = env, .body = i_earth, .gm_object = gm_moon, .k2 = k2,
                             .lag = 6.0e5};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "earth", myb_term_point_mass,
                               myb_term_point_mass_jacobian, &central, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "tide", myb_term_tidal_lag, NULL, &term, NULL) == MYB_OK);
  const double y0[6] = {a0, 0.0, 0.0, 0.0, n * a0, 0.0};
  const myb_orbit_options orbit_options = {.epsilon = 1e-11};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, &orbit_options, NULL, &orbit) == MYB_OK);
  const double span = 10.0 * two_pi / n;
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, span) == MYB_OK);
  double y[6];
  MT_CHECK(myb_orbit_state(orbit, span, y) == MYB_OK);
  double a1 = 0.0;
  double e1 = 0.0;
  elements(mu, y, &a1, &e1);
  printf("  ten revolutions with a delay of 6e5 s: a grows by %.6f m, the formula %.6f m; e "
         "%.1e\n",
         a1 - a0, 1000.0 * rate * span, e1);
  MT_CHECK(fabs((a1 - a0) / (1000.0 * rate * span) - 1.0) < 1e-5);
  /* a retrograde spin pulls the orbit in */
  earth_spin.rate = -spin;
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, span) == MYB_OK);
  MT_CHECK(myb_orbit_state(orbit, span, y) == MYB_OK);
  elements(mu, y, &a1, &e1);
  MT_CHECK(fabs((a1 - a0) / (1000.0 * rate * span * (-spin - n) / (spin - n)) - 1.0) < 1e-5);
  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);
  myb_env_destroy(env);
}

/* ---- Recoil, a changing mass, impulses -------------------------------------------------------- */

static void test_recoil_and_mass(void) {
  /* Pioneer 10: 63 W of its heat leaving along the antenna would give the acceleration seen */
  const double along[3] = {0.0, 0.0, 1.0};
  double a[3] = {0.0, 0.0, 0.0};
  myb_force_emission(63.0, 1.0, 241.0, c_light, along, a);
  printf("recoil of 63 W on 241 kg: %.2e m/s^2 (the Pioneer anomaly: 8.74e-10)\n", -a[2]);
  MT_CHECK(fabs(a[2] + 63.0 / (241.0 * c_light)) < 1e-24 && -a[2] > 8.5e-10 && -a[2] < 9.0e-10);
  /* the term: a body that spins about z carries the direction round; the power decays */
  myb_object craft = {.mass = 241.0};
  myb_spin_attitude spin = {.m0 = {1, 0, 0, 0, 1, 0, 0, 0, 1}, .rate = 0.5, .phase = 0.0};
  myb_recoil_term recoil = {.object = &craft, .c = c_light, .power = 2500.0, .efficiency = 0.1,
                            .direction = {1.0, 0.0, 0.0}, .attitude = myb_attitude_spin,
                            .attitude_context = &spin};
  const double x[3] = {7.0e6, 0.0, 0.0};
  const double v[3] = {0.0, 7.5e3, 0.0};
  const double size = 0.1 * 2500.0 / (241.0 * c_light);
  MT_CHECK(myb_term_recoil(&recoil, 0.0, x, v, a) == MYB_OK);
  MT_CHECK(fabs(a[0] + size) < 1e-24 && a[1] == 0.0 && a[2] == 0.0);
  MT_CHECK(myb_term_recoil(&recoil, pi, x, v, a) == MYB_OK); /* a quarter of a turn */
  MT_CHECK(fabs(a[0]) < 1e-24 && fabs(a[1] + size) < 1e-24);
  recoil.decay = 100.0;
  recoil.epoch = pi - 100.0;
  MT_CHECK(myb_term_recoil(&recoil, pi, x, v, a) == MYB_OK);
  MT_CHECK(fabs(a[1] + size * exp(-1.0)) < 1e-24);
  recoil.direction[0] = 2.0;
  MT_CHECK(myb_term_recoil(&recoil, pi, x, v, a) == MYB_E_DOMAIN);

  /* a mass that burns away: the light pushes what is left harder */
  const myb_thrust_arc burn = {.start = 100.0, .stop = 200.0, .thrust = 400.0, .exhaust = 3000.0,
                               .steering = MYB_STEER_INERTIAL, .direction = {1.0, 0.0, 0.0}};
  myb_thrust thrust;
  MT_CHECK(myb_thrust_init(&thrust, 1000.0, &burn, 1) == MYB_OK);
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body sun = {.id = 10, .gm = gm_sun, .luminosity = 3.828e26};
  size_t i_sun = 0;
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  myb_object tank = {.mass = 1000.0, .area = 10.0, .cr = 1.2, .mass_function = myb_thrust_mass_at,
                     .mass_context = &thrust};
  myb_radiation_term light = {.env = env, .source = i_sun, .object = &tank, .c = c_light};
  const double at[3] = {au, 0.0, 0.0};
  double before[3];
  double after[3];
  MT_CHECK(myb_term_radiation(&light, 50.0, at, v, before) == MYB_OK);
  MT_CHECK(myb_term_radiation(&light, 500.0, at, v, after) == MYB_OK);
  const double left = 1000.0 - 400.0 / 3000.0 * 100.0;
  MT_CHECK(myb_thrust_mass_at(&thrust, 500.0) == myb_thrust_mass(&thrust, 500.0));
  MT_CHECK(fabs(myb_thrust_mass(&thrust, 500.0) - left) < 1e-9);
  MT_CHECK(fabs(after[0] / before[0] - 1000.0 / left) < 1e-14);
  tank.mass_function = NULL;
  double constant[3];
  MT_CHECK(myb_term_radiation(&light, 500.0, at, v, constant) == MYB_OK);
  MT_CHECK(constant[0] == before[0]);
  myb_env_destroy(env);
}

static void test_impulses(void) {
  const double mu = 3.986004418e14;
  myb_point_mass_term central = {.gm = mu};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "earth", myb_term_point_mass,
                               myb_term_point_mass_jacobian, &central, NULL) == MYB_OK);
  /* Hohmann: 7000 km to 42 164 km */
  const double r1 = 7.0e6;
  const double r2 = 4.2164e7;
  const double at = 0.5 * (r1 + r2);
  const double v1 = sqrt(mu / r1);
  const double dv1 = sqrt(mu * (2.0 / r1 - 1.0 / at)) - v1;
  const double dv2 = sqrt(mu / r2) - sqrt(mu * (2.0 / r2 - 1.0 / at));
  const double half = pi * sqrt(at * at * at / mu);
  const double y0[6] = {r1, 0.0, 0.0, 0.0, v1, 0.0};
  const double kick1[3] = {0.0, dv1, 0.0};
  const double kick2[3] = {0.0, -dv2, 0.0}; /* at the apocentre the motion is along -y */
  const double span = half + 30000.0;
  double final[6] = {0.0};
  double phi[36];
  for (int variational = 0; variational < 2; ++variational) {
    const myb_orbit_options options = {.epsilon = 1e-11, .variational = variational};
    myb_orbit* orbit = NULL;
    MT_CHECK(myb_orbit_create(model, &options, NULL, &orbit) == MYB_OK);
    MT_CHECK(myb_orbit_start(orbit, 0.0, y0, span) == MYB_OK);
    MT_CHECK(myb_orbit_kick(orbit, kick1) == MYB_OK);
    double y[6];
    MT_CHECK(myb_orbit_state(orbit, 0.0, y) == MYB_OK); /* the state after the impulse */
    MT_CHECK(y[4] == v1 + dv1 && y[0] == r1);
    MT_CHECK(myb_orbit_state(orbit, half, y) == MYB_OK);
    MT_CHECK(fabs(y[0] + r2) < 1e-4 && fabs(y[1]) < 1e-3);
    MT_CHECK(myb_orbit_kick(orbit, kick2) == MYB_OK);
    if (variational) {
      MT_CHECK(myb_orbit_partials(orbit, span, y, phi, NULL) == MYB_OK);
      for (int i = 0; i < 6; ++i) {
        MT_CHECK(fabs(y[i] - final[i]) < 1e-5 * (i < 3 ? 1.0 : 1e-3));
      }
    } else {
      MT_CHECK(myb_orbit_state(orbit, span, final) == MYB_OK);
      double a = 0.0;
      double e = 0.0;
      elements(mu, final, &a, &e);
      printf("Hohmann transfer by two impulses: the final orbit %.2e m from 42 164 km, "
             "e = %.1e\n",
             a - r2, e);
      MT_CHECK(fabs(a - r2) < 1e-3 && e < 1e-11);
      const double bad[3] = {(double)NAN, 0.0, 0.0};
      MT_CHECK(myb_orbit_kick(orbit, bad) == MYB_E_DOMAIN);
    }
    myb_orbit_destroy(orbit);
  }
  /* the transition matrix across the impulses against differences of the whole manoeuvre */
  double worst = 0.0;
  double largest = 0.0;
  for (int j = 0; j < 6; ++j) {
    const double h = j < 3 ? 10.0 : 1.0e-2;
    double ends[2][6];
    for (int side = 0; side < 2; ++side) {
      double q[6];
      memcpy(q, y0, sizeof q);
      q[j] += side == 0 ? h : -h;
      const myb_orbit_options options = {.epsilon = 1e-12};
      myb_orbit* orbit = NULL;
      MT_CHECK(myb_orbit_create(model, &options, NULL, &orbit) == MYB_OK);
      MT_CHECK(myb_orbit_start(orbit, 0.0, q, span) == MYB_OK);
      MT_CHECK(myb_orbit_kick(orbit, kick1) == MYB_OK);
      double y[6];
      MT_CHECK(myb_orbit_state(orbit, half, y) == MYB_OK);
      MT_CHECK(myb_orbit_kick(orbit, kick2) == MYB_OK);
      MT_CHECK(myb_orbit_state(orbit, span, ends[side]) == MYB_OK);
      myb_orbit_destroy(orbit);
    }
    for (int i = 0; i < 6; ++i) {
      const double scale = (i < 3 ? 1.0 : 1.0e-3) * (j < 3 ? 1.0 : 1.0e3);
      worst = fmax(worst, fabs((ends[0][i] - ends[1][i]) / (2.0 * h) - phi[6 * i + j]) / scale);
      largest = fmax(largest, fabs(phi[6 * i + j]) / scale);
    }
  }
  printf("  the transition matrix across the impulses: within %.1e of differences (largest "
         "element %.1f)\n",
         worst, largest);
  MT_CHECK(worst < 1e-6 * largest);
  myb_force_model_destroy(model);
}

/* The Jacobian of a term against central differences: the largest error over the largest
 * element, for the position and for the velocity. */
static void jacobian_check(myb_term_function function, myb_term_jacobian jacobian, void* term,
                           const double x[3], const double v[3], double* error_x,
                           double* error_v) {
  double dadx[9];
  double dadv[9];
  MT_CHECK(jacobian(term, 100.0, x, v, dadx, dadv) == MYB_OK);
  double worst[2] = {0.0, 0.0};
  double size[2] = {0.0, 0.0};
  for (int j = 0; j < 6; ++j) {
    double q[6] = {x[0], x[1], x[2], v[0], v[1], v[2]};
    const double h = 1e-5 * (j < 3 ? norm3(x) : norm3(v));
    double plus[3];
    double minus[3];
    q[j] += h;
    MT_CHECK(function(term, 100.0, q, q + 3, plus) == MYB_OK);
    q[j] -= 2.0 * h;
    MT_CHECK(function(term, 100.0, q, q + 3, minus) == MYB_OK);
    const double* analytic = j < 3 ? dadx : dadv;
    for (int i = 0; i < 3; ++i) {
      const double value = analytic[3 * i + j % 3];
      worst[j / 3] = fmax(worst[j / 3], fabs((plus[i] - minus[i]) / (2.0 * h) - value));
      size[j / 3] = fmax(size[j / 3], fabs(value));
    }
  }
  *error_x = size[0] > 0.0 ? worst[0] / size[0] : worst[0];
  *error_v = size[1] > 0.0 ? worst[1] / size[1] : worst[1];
}

static void test_jacobians(void) {
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  spinner turning = {.rate = 7.292115e-5};
  const myb_env_body sun = {.id = 10, .gm = gm_sun, .radius = 6.957e8,
                            .orientation = spinner_orientation, .orientation_context = &turning};
  size_t i_sun = 0;
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  const double x[3] = {1.2 * au, 0.3 * au, -0.2 * au};
  const double v[3] = {-4.0e3, 2.6e4, 1.0e3};
  double ex = 0.0;
  double ev = 0.0;
  /* comets and asteroids */
  myb_nongrav_term comet = {.env = env, .source = i_sun};
  myb_nongrav_marsden(1.0e-8, 2.0e-9, -3.0e-10, &comet.model);
  jacobian_check(myb_term_nongrav, myb_term_nongrav_jacobian, &comet, x, v, &ex, &ev);
  printf("Jacobians against differences: comet %.1e, %.1e", ex, ev);
  MT_CHECK(ex < 1e-8 && ev < 1e-8);
  myb_nongrav_power(0.0, 1.0e-10, 0.0, au, 2.0, &comet.model);
  jacobian_check(myb_term_nongrav, myb_term_nongrav_jacobian, &comet, x, v, &ex, &ev);
  MT_CHECK(ex < 1e-8 && ev < 1e-8);
  /* with a delay: differences of the term inside, the same answer to their accuracy */
  myb_nongrav_marsden(1.0e-8, 2.0e-9, -3.0e-10, &comet.model);
  comet.delay = 3.0e6;
  jacobian_check(myb_term_nongrav, myb_term_nongrav_jacobian, &comet, x, v, &ex, &ev);
  MT_CHECK(ex < 1e-4 && ev < 1e-4);
  /* the lagging tide, with and without the tide that does not lag */
  myb_tidal_lag_term tide = {.env = env, .body = i_sun, .gm_object = 1.0e17, .k2 = 0.3,
                             .lag = 600.0};
  const double near[3] = {3.0e9, -1.0e9, 4.0e8};
  const double fast[3] = {5.0e4, 1.9e5, -2.0e4};
  jacobian_check(myb_term_tidal_lag, myb_term_tidal_lag_jacobian, &tide, near, fast, &ex, &ev);
  printf("; lagging tide %.1e, %.1e", ex, ev);
  MT_CHECK(ex < 1e-8 && ev < 1e-8);
  tide.with_static = 1;
  jacobian_check(myb_term_tidal_lag, myb_term_tidal_lag_jacobian, &tide, near, fast, &ex, &ev);
  MT_CHECK(ex < 1e-8 && ev < 1e-8);
  /* empirical accelerations */
  myb_empirical_term empirical = {.constant = {1e-9, -2e-9, 3e-10},
                                  .cosine = {4e-10, 5e-9, -6e-9},
                                  .sine = {-7e-10, 8e-10, 9e-9}};
  const double leo[3] = {-2.0e6, 6.5e6, 1.5e6};
  const double leo_v[3] = {-7.0e3, -2.0e3, 1.5e3};
  jacobian_check(myb_term_empirical, myb_term_empirical_jacobian, &empirical, leo, leo_v, &ex,
                 &ev);
  printf("; empirical %.1e, %.1e\n", ex, ev);
  MT_CHECK(ex < 1e-8 && ev < 1e-8);
  /* an orbit in the plane the node is counted on: the node is the x axis there */
  const double flat[3] = {7.0e6, 1.0e6, 0.0};
  const double flat_v[3] = {-1.0e3, 7.4e3, 0.0};
  double dadx[9];
  double dadv[9];
  MT_CHECK(myb_term_empirical_jacobian(&empirical, 0.0, flat, flat_v, dadx, dadv) == MYB_OK);
  MT_CHECK(myb_term_empirical_jacobian(&empirical, 0.0, flat, flat, dadx, dadv) == MYB_E_DOMAIN);
  myb_env_destroy(env);
}

static void test_delay_and_yarkovsky(void) {
  /* the outgassing behind the heating: a time `delay` after the perihelion the comet is
   * pushed as hard as an undelayed one at the perihelion */
  myb_nongrav comet;
  myb_nongrav_marsden(1.0e-8, 0.0, 0.0, &comet);
  double y0[6];
  state_of(gm_sun, 1.3 * au, 0.5, 0.2, y0); /* at the perihelion */
  const double delay = 2.0e6;
  myb_pv start;
  myb_pv later;
  memcpy(start.r, y0, sizeof start.r);
  memcpy(start.v, y0 + 3, sizeof start.v);
  MT_CHECK(myb_kepler_propagate(gm_sun, &start, delay, &later) == MYB_OK);
  double a[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_nongrav_delayed(&comet, gm_sun, delay, later.r, later.v, a) == MYB_OK);
  const double at_perihelion = 1.0e-8 * myb_nongrav_g(&comet, norm3(y0));
  double now[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_nongrav(&comet, later.r, later.v, now) == MYB_OK);
  printf("outgassing 2e6 s behind the heating: %.3e m/s^2 that long after the perihelion "
         "(at the perihelion without delay %.3e, here without delay %.3e)\n",
         norm3(a), at_perihelion, norm3(now));
  MT_CHECK(fabs(norm3(a) / at_perihelion - 1.0) < 1e-9);
  MT_CHECK(norm3(now) < 0.99 * norm3(a));
  double same[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_nongrav_delayed(&comet, gm_sun, 0.0, later.r, later.v, same) == MYB_OK);
  MT_CHECK(same[0] == now[0] && same[1] == now[1] && same[2] == now[2]);
  MT_CHECK(myb_force_nongrav_delayed(&comet, -1.0, delay, later.r, later.v, same) ==
           MYB_E_DOMAIN);

  /* the answer of a conducting sphere: values from the solution of the heat equation in
   * mpmath (1 / (1 + lambda z j1'(z) / j1(z)), 30 digits) */
  static const double reference[8][4] = {
      {0.3, 0.5, 0.458974052846, -0.00446909135886},
      {1.0, 1.0, 0.580451335842, -0.0479120817778},
      {2.5, 3.0, 0.395242313955, -0.197762483771},
      {10.0, 0.7, 0.716457682328, -0.190092729588},
      /* small bodies: the series, and both sides of where it hands over to the closed form */
      {0.02, 0.5, 0.0535398697259013, -4.05386816020604e-6},
      {0.03, 1.5, 0.0275062753018137, -4.81494239033833e-6},
      {0.28, 1.5, 0.2088143943414327, -0.00259029118415612},
      {0.29, 1.5, 0.2146655551158022, -0.002835311889725602}};
  double worst = 0.0;
  double worst_small = 0.0;
  for (int k = 0; k < 8; ++k) {
    double re = 0.0;
    double im = 0.0;
    myb_yarkovsky_response(reference[k][0], reference[k][1], &re, &im);
    if (k < 4) {
      worst = fmax(worst, fmax(fabs(re - reference[k][2]), fabs(im - reference[k][3])));
    } else {
      worst_small = fmax(worst_small, fmax(fabs(re / reference[k][2] - 1.0),
                                           fabs(im / reference[k][3] - 1.0)));
    }
  }
  double re = 0.0;
  double im = 0.0;
  /* a very large body: the half-space; a very small one: the same temperature all over */
  myb_yarkovsky_response(0.0, 2.0, &re, &im);
  MT_CHECK(fabs(re - 0.4) < 1e-15 && fabs(im + 0.2) < 1e-15);
  double re_big = 0.0;
  double im_big = 0.0;
  myb_yarkovsky_response(400.0, 2.0, &re_big, &im_big);
  MT_CHECK(fabs(re_big - 0.4) < 2e-3 && fabs(im_big + 0.2) < 2e-3);
  myb_yarkovsky_response(1.0e-3, 2.0, &re, &im);
  MT_CHECK(fabs(im) < 1e-6 && re > 0.0 && re < 1.0);
  printf("the answer of a conducting sphere: within %.1e of the solution of the heat "
         "equation; small bodies, each part within %.1e of itself\n",
         worst, worst_small);
  MT_CHECK(worst < 1e-11);
  MT_CHECK(worst_small < 1e-9);

  /* a small body drifts less than the large-body limit says; the term on a circular orbit
   * gives the drift of the formula */
  myb_yarkovsky body = {.radius = 5.0, .density = 2500.0, .absorptivity = 0.9,
                        .emissivity = 0.9, .thermal_inertia = 200.0,
                        .rotation_rate = two_pi / 600.0, .obliquity = 0.6};
  const double radius_orbit = 1.2 * au;
  const double flux = 1361.0 / (1.2 * 1.2);
  const double n = sqrt(gm_sun / (radius_orbit * radius_orbit * radius_orbit));
  double large = 0.0;
  double finite = 0.0;
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &large, NULL) == MYB_OK);
  body.conductivity = 0.5; /* the daily wave reaches 2.4 cm: a body of 5 m is large; ... */
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &finite, NULL) == MYB_OK);
  MT_CHECK(fabs(finite / large - 1.0) < 0.02);
  body.radius = 0.002; /* ... one of 2 mm is not */
  double small = 0.0;
  double small_large = 0.0;
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &small, NULL) == MYB_OK);
  body.conductivity = 0.0;
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &small_large, NULL) == MYB_OK);
  printf("a grain of 2 mm: the diurnal drift %.2f of what the large-body limit gives\n",
         small / small_large);
  MT_CHECK(small > 0.0 && small < 0.5 * small_large);
  body.radius = 5.0;
  body.conductivity = 0.5;

  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body sun = {.id = 10, .gm = gm_sun, .luminosity = 1361.0 * 4.0 * pi * au * au};
  size_t i_sun = 0;
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  myb_point_mass_term central = {.gm = gm_sun};
  /* the orbit in the xy plane; the spin 0.6 rad from its normal */
  myb_yarkovsky_term term = {.env = env, .source = i_sun, .body = body,
                             .spin = {sin(0.6), 0.0, cos(0.6)}, .c = c_light};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "sun", myb_term_point_mass, myb_term_point_mass_jacobian,
                               &central, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "yarkovsky", myb_term_yarkovsky, NULL, &term, NULL) ==
           MYB_OK);
  const double circular[6] = {radius_orbit, 0.0, 0.0, 0.0, n * radius_orbit, 0.0};
  const myb_orbit_options orbit_options = {.epsilon = 1e-11};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, &orbit_options, NULL, &orbit) == MYB_OK);
  const double span = 3.0 * two_pi / n;
  MT_CHECK(myb_orbit_start(orbit, 0.0, circular, span) == MYB_OK);
  double y[6];
  MT_CHECK(myb_orbit_state(orbit, span, y) == MYB_OK);
  double a1 = 0.0;
  double e1 = 0.0;
  elements(gm_sun, y, &a1, &e1);
  MT_CHECK(myb_yarkovsky_drift(&body, flux, n, c_light, &finite, NULL) == MYB_OK);
  printf("the diurnal Yarkovsky term on a circular orbit, 5 m of rock: da/dt %.6e m/s, the "
         "formula %.6e\n",
         (a1 - radius_orbit) / span, finite);
  MT_CHECK(fabs((a1 - radius_orbit) / span / finite - 1.0) < 1e-4);
  /* the force: across the spin axis only; no thermal inertia: straight away from the Sun */
  const double to_sun[3] = {-1.0, 0.0, 0.0};
  double f[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_yarkovsky(&body, term.spin, flux, c_light, to_sun, f) == MYB_OK);
  MT_CHECK(fabs(f[0] * term.spin[0] + f[2] * term.spin[2]) < 1e-20 * (1.0 + norm3(f) * 1e20));
  MT_CHECK(f[1] != 0.0); /* the afternoon side */
  body.thermal_inertia = 0.0;
  double cold[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_yarkovsky(&body, term.spin, flux, c_light, to_sun, cold) == MYB_OK);
  MT_CHECK(cold[1] == 0.0 && cold[0] > 0.0);
  const double skew[3] = {1.0, 1.0, 0.0};
  MT_CHECK(myb_force_yarkovsky(&body, skew, flux, c_light, to_sun, cold) == MYB_E_DOMAIN);
  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);
  myb_env_destroy(env);
}

static void test_tide_and_figure(void) {
  /* what the orbit gains in angular momentum, the body loses */
  const double gm_earth = 3.986004418e14;
  const double gm_moon = 4.9028e12;
  const double omega[3] = {1.0e-6, -2.0e-6, 7.292115e-5};
  const double r[3] = {3.6e8, 1.2e8, -3.0e7};
  const double v[3] = {-300.0, 950.0, 80.0};
  double a[3] = {0.0, 0.0, 0.0};
  double torque[3];
  MT_CHECK(myb_force_tidal_lag(gm_earth, gm_moon, 0.3, 6378136.6, 600.0, omega, r, v, 0, a) ==
           MYB_OK);
  MT_CHECK(myb_torque_tidal_lag(gm_earth, gm_moon, 0.3, 6378136.6, 600.0, omega, r, v,
                                torque) == MYB_OK);
  /* per unit mass of the body: the orbit's rate is (m / (M + m)) r x a */
  const double share = gm_moon / (gm_earth + gm_moon);
  const double orbit[3] = {share * (r[1] * a[2] - r[2] * a[1]), share * (r[2] * a[0] - r[0] * a[2]),
                           share * (r[0] * a[1] - r[1] * a[0])};
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(orbit[i] + torque[i]) <= 1e-15 * norm3(orbit));
  }
  MT_CHECK(torque[2] < 0.0); /* the Moon brakes the Earth */

  /* the figure of the Earth answering the Moon: with J2 alone and the Moon in the equator,
   * (m / M) (-3/2) J2 mu R^2 / r^4 along the radius */
  spinner still = {.rate = 0.0};
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body earth = {.id = 399, .gm = gm_earth, .radius = 6378136.6,
                              .orientation = spinner_orientation, .orientation_context = &still};
  size_t i_earth = 0;
  MT_CHECK(myb_env_add_body(env, &earth, &i_earth) == MYB_OK);
  /* the Moon at rest on the x axis: a fixed state */
  static double moon_state[6] = {3.844e8, 0.0, 0.0, 0.0, 0.0, 0.0};
  const myb_env_body moon = {.id = 301, .gm = gm_moon, .state = fixed_state,
                             .state_context = moon_state};
  size_t i_moon = 0;
  MT_CHECK(myb_env_add_body(env, &moon, &i_moon) == MYB_OK);
  double cbar[6] = {1.0, 0.0, 0.0, -4.84165e-4, 0.0, 0.0};
  double sbar[6] = {0.0};
  myb_gravity_field* field = NULL;
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_field_create(gm_earth, 6378136.3, 2, cbar, sbar, NULL, &field) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(2, NULL, &work) == MYB_OK);
  myb_field_term figure = {.env = env, .body = i_earth, .field = field, .work = work,
                           .n_max = -1, .raisers = {i_moon}, .raiser_count = 1};
  const double x[3] = {7.0e6, 0.0, 0.0};
  MT_CHECK(myb_term_field_indirect(&figure, 0.0, x, v, a) == MYB_OK);
  const double j2 = 4.84165e-4 * sqrt(5.0);
  const double expected = -(gm_moon / gm_earth) * 1.5 * j2 * gm_earth * 6378136.3 * 6378136.3 /
                          pow(3.844e8, 4.0);
  printf("the Earth's oblateness answering the Moon: %.4e m/s^2 on the motion about the "
         "Earth (J2 alone: %.4e)\n",
         a[0], expected);
  MT_CHECK(fabs(a[0] / expected - 1.0) < 1e-12 && fabs(a[1]) < 1e-25 && fabs(a[2]) < 1e-25);
  /* it does not depend on the object */
  const double elsewhere[3] = {-1.0e7, 4.0e6, 2.0e6};
  double b[3];
  MT_CHECK(myb_term_field_indirect(&figure, 0.0, elsewhere, v, b) == MYB_OK);
  MT_CHECK(b[0] == a[0] && b[1] == a[1] && b[2] == a[2]);
  /* a field of a body that is not at the origin has no such term */
  myb_field_term other = figure;
  other.body = i_moon;
  MT_CHECK(myb_term_field_indirect(&other, 0.0, x, v, b) == MYB_E_DOMAIN);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
  myb_env_destroy(env);
}

/* The seasonal Yarkovsky effect on eccentric orbits against mpmath (25 digits; the harmonics
 * of the heating by sums of another number of steps, the answer of the sphere from the
 * spherical Bessel function; yark_season.py): the acceleration along the spin axis at three
 * places of the orbit and the mean drift. */
static void test_yarkovsky_season(void) {
  typedef struct reference {
    double e;
    double radius;
    double conductivity;
    double along[3]; /* at the mean anomalies 0.3, 2.0, 4.5 */
    double drift;
    int harmonics;
  } reference;
  static const reference cases[5] = {
      /* a circle, a large body */
      {0.0, 50.0, 0.0, {6.6647062495195162e-13, -1.5977238553063296e-12, 1.0020253707919178e-12}, -7.8307881793096769e-7, 2},
      /* e = 0.3, a large body */
      {0.3, 50.0, 0.0, {4.8837553081927934e-13, -1.1797327396349474e-12, 1.7713584717537545e-13}, -8.8854035815320799e-7, 46},
      /* e = 0.3, half a metre of rock that conducts poorly */
      {0.3, 0.5, 0.1, {3.1929300034211549e-11, -1.0552607108912313e-10, 1.8775606349864179e-11}, -8.6239283937240573e-6, 46},
      /* e = 0.6, two metres of rock */
      {0.6, 2.0, 1.0, {-6.0675717336783525e-11, -1.3757378620906335e-11, -3.5516955987858973e-12}, -1.5393699193811222e-6, 136},
      /* e = 0.85 */
      {0.85, 10.0, 2.5, {-2.4870895112111151e-11, -2.0796952259320378e-12, -1.2879502995995185e-12}, -5.6987126380019621e-6, 682}};
  static const double anomalies[3] = {0.3, 2.0, 4.5};
  const double gm = 1.32712440018e20;
  const double a = 2.2 * au;
  /* the axes of the orbit and the spin in them */
  const double p[3] = {0.36, 0.48, 0.8};
  const double q[3] = {-0.8, 0.6, 0.0};
  const double h[3] = {-0.48, -0.64, 0.6};
  double spin[3];
  for (int i = 0; i < 3; ++i) {
    spin[i] = 0.48 * p[i] - 0.6 * q[i] + 0.64 * h[i];
  }
  double worst_along = 0.0;
  double worst_drift = 0.0;
  for (int k = 0; k < 5; ++k) {
    const reference* ref = &cases[k];
    const myb_yarkovsky body = {.radius = ref->radius, .density = 2500.0, .absorptivity = 0.9,
                                .emissivity = 0.9, .thermal_inertia = 300.0,
                                .rotation_rate = 1.0e-3, .conductivity = ref->conductivity};
    /* the state at the true anomaly 1.1 (at the pericentre for the circle, which has none) */
    const double e = ref->e;
    const double f = e > 0.0 ? 1.1 : 0.0;
    const double semi = a * (1.0 - e * e);
    const double r = semi / (1.0 + e * cos(f));
    const double speed = sqrt(gm / semi);
    double state[6];
    for (int i = 0; i < 3; ++i) {
      state[i] = r * (cos(f) * p[i] + sin(f) * q[i]);
      state[3 + i] = speed * (-sin(f) * p[i] + (e + cos(f)) * q[i]);
    }
    myb_yarkovsky_season* season = NULL;
    MT_CHECK(myb_yarkovsky_season_create(&body, spin, gm, state, 1361.0, au, c_light, NULL,
                                         &season) == MYB_OK);
    if (season == NULL) {
      continue;
    }
    double scale = 0.0;
    for (int i = 0; i < 3; ++i) {
      scale = fmax(scale, fabs(ref->along[i]));
    }
    for (int i = 0; i < 3; ++i) {
      /* the place of the mean anomaly: Kepler's equation by Newton's steps */
      double ea = anomalies[i];
      for (int it = 0; it < 60; ++it) {
        ea -= (ea - e * sin(ea) - anomalies[i]) / (1.0 - e * cos(ea));
      }
      double at[3];
      for (int j = 0; j < 3; ++j) {
        at[j] = a * ((cos(ea) - e) * p[j] + sqrt(1.0 - e * e) * sin(ea) * q[j]);
      }
      double acc[3] = {0.0, 0.0, 0.0};
      myb_force_yarkovsky_season(season, at, acc);
      const double along = acc[0] * spin[0] + acc[1] * spin[1] + acc[2] * spin[2];
      worst_along = fmax(worst_along, fabs(along - ref->along[i]) / scale);
      /* along the axis and nothing else */
      for (int j = 0; j < 3; ++j) {
        MT_CHECK(fabs(acc[j] - along * spin[j]) <= 1e-15 * fabs(along));
      }
    }
    const double drift = myb_yarkovsky_season_drift(season);
    worst_drift = fmax(worst_drift, fabs(drift / ref->drift - 1.0));
    MT_CHECK(drift < 0.0);
    printf("  e = %.2f, R = %g m: %zu harmonics (mpmath took %d), da/dt = %.6e m/s = %.3f km "
           "in a million years\n",
           e, ref->radius, myb_yarkovsky_season_harmonics(season), ref->harmonics, drift,
           drift * 3.15576e13 / 1000.0);
    if (k == 0) {
      /* a circle: the closed form of the drift */
      const double n = sqrt(gm / a) / a;
      myb_yarkovsky tilted = body;
      tilted.obliquity = acos(0.64);
      double seasonal = 0.0;
      MT_CHECK(myb_yarkovsky_drift(&tilted, 1361.0 / (2.2 * 2.2), n, c_light, NULL, &seasonal) ==
               MYB_OK);
      printf("  on a circle: the mean drift within %.1e of the closed form\n",
             fabs(drift / seasonal - 1.0));
      MT_CHECK(fabs(drift / seasonal - 1.0) < 1e-13);
    }
    myb_yarkovsky_season_destroy(season);
  }
  printf("the seasonal Yarkovsky effect on eccentric orbits against mpmath: accelerations "
         "within %.1e of their size, drifts within %.1e\n",
         worst_along, worst_drift);
  MT_CHECK(worst_along < 1e-11 && worst_drift < 1e-11);

  /* in a model: the term adds it to the diurnal part */
  const myb_yarkovsky body = {.radius = 2.0, .density = 2500.0, .absorptivity = 0.9,
                              .emissivity = 0.9, .thermal_inertia = 300.0,
                              .rotation_rate = 1.0e-3, .conductivity = 1.0};
  const double e = 0.6;
  const double state[6] = {a * (1.0 - e) * p[0], a * (1.0 - e) * p[1], a * (1.0 - e) * p[2],
                           sqrt(gm * (1.0 + e) / (a * (1.0 - e))) * q[0],
                           sqrt(gm * (1.0 + e) / (a * (1.0 - e))) * q[1],
                           sqrt(gm * (1.0 + e) / (a * (1.0 - e))) * q[2]};
  myb_yarkovsky_season* season = NULL;
  MT_CHECK(myb_yarkovsky_season_create(&body, spin, gm, state, 1361.0, au, c_light, NULL,
                                       &season) == MYB_OK);
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  /* the luminosity that gives 1361 W/m^2 at 1 au */
  const myb_env_body sun = {.id = 10, .gm = gm, .luminosity = 1361.0 * 4.0 * pi * au * au};
  size_t i_sun = 0;
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  myb_yarkovsky_term term = {.env = env, .source = i_sun, .body = body,
                             .spin = {spin[0], spin[1], spin[2]}, .c = c_light};
  const double x[3] = {0.3 * a * p[0] + 1.1 * a * q[0], 0.3 * a * p[1] + 1.1 * a * q[1],
                       0.3 * a * p[2] + 1.1 * a * q[2]};
  double diurnal[3];
  double both[3];
  double seasonal[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_term_yarkovsky(&term, 0.0, x, state + 3, diurnal) == MYB_OK);
  term.season = season;
  MT_CHECK(myb_term_yarkovsky(&term, 0.0, x, state + 3, both) == MYB_OK);
  myb_force_yarkovsky_season(season, x, seasonal);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(both[i] - diurnal[i] - seasonal[i]) <= 1e-15 * fabs(both[i]) + 1e-30);
  }
  MT_CHECK(fabs(seasonal[0]) + fabs(seasonal[1]) + fabs(seasonal[2]) > 1e-13);
  myb_env_destroy(env);
  myb_yarkovsky_season_destroy(season);

  /* what is refused: an orbit too eccentric or open, a spin that is not a unit vector */
  double fast[6];
  memcpy(fast, state, sizeof fast);
  for (int i = 3; i < 6; ++i) {
    fast[i] *= 1.105; /* e = 0.954 */
  }
  season = NULL;
  MT_CHECK(myb_yarkovsky_season_create(&body, spin, gm, fast, 1361.0, au, c_light, NULL,
                                       &season) == MYB_E_DOMAIN);
  for (int i = 3; i < 6; ++i) {
    fast[i] *= 1.2; /* a hyperbola */
  }
  MT_CHECK(myb_yarkovsky_season_create(&body, spin, gm, fast, 1361.0, au, c_light, NULL,
                                       &season) == MYB_E_DOMAIN);
  const double long_spin[3] = {0.0, 0.0, 1.1};
  MT_CHECK(myb_yarkovsky_season_create(&body, long_spin, gm, state, 1361.0, au, c_light, NULL,
                                       &season) == MYB_E_DOMAIN);
  MT_CHECK(season == NULL);
  myb_yarkovsky_season_destroy(NULL);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_jacobians();
  test_delay_and_yarkovsky();
  test_yarkovsky_season();
  test_tide_and_figure();
  test_poynting_robertson();
  test_small_bodies();
  test_yarkovsky_theory();
  test_tidal_lag();
  test_recoil_and_mass();
  test_impulses();
  return mt_finish("test_nongrav");
}
