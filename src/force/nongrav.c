/* SPDX-License-Identifier: BSD-3-Clause */
/* Small non-gravitational forces (myb/force.h): radiation on a moving object
 * (Poynting-Robertson), the accelerations of comets and of the Yarkovsky effect, the drift of
 * the linear Yarkovsky theory, the recoil of what an object emits, the lagging tide. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "myb/force.h"
#include "myb/kepler.h"
#include "num/vec3.h"

static const double au = 1.495978707e11;
static const double stefan_boltzmann = 5.670374419e-8; /* W m^-2 K^-4, exact since 2019 */

/* J. A. Burns, P. L. Lamy, S. Soter, Icarus 40, 1 (1979), eq. 2. */
void myb_force_radiation_moving(double pressure, double cr, double area_over_mass,
                                const double u[3], const double v[3], double c, double a[3]) {
  MYB_CHECK(u != NULL && v != NULL && a != NULL && c > 0.0);
  const double k = pressure * cr * area_over_mass;
  const double radial = myb_vec3_dot(u, v) / c;
  for (int i = 0; i < 3; ++i) {
    a[i] += k * ((1.0 - radial) * u[i] - v[i] / c);
  }
}

/* x^y for x > 0. */
static double power(double x, double y) {
  return myb_exp(y * myb_log(x));
}

void myb_nongrav_marsden(double a1, double a2, double a3, myb_nongrav* model) {
  MYB_CHECK(model != NULL);
  /* B. G. Marsden, Z. Sekanina, D. K. Yeomans, Astron. J. 78, 211 (1973): the sublimation of
   * water ice; alpha makes g(1 au) = 1 */
  *model = (myb_nongrav){.a = {a1, a2, a3},
                         .r0 = 2.808 * au,
                         .alpha = 0.1112620426,
                         .m = 2.15,
                         .n = 5.093,
                         .k = 4.6142};
}

void myb_nongrav_power(double a1, double a2, double a3, double r0, double exponent,
                       myb_nongrav* model) {
  MYB_CHECK(model != NULL);
  *model = (myb_nongrav){.a = {a1, a2, a3}, .r0 = r0, .alpha = 1.0, .m = exponent};
}

double myb_nongrav_g(const myb_nongrav* model, double r) {
  MYB_CHECK(model != NULL);
  if (!(r > 0.0) || !(model->r0 > 0.0)) {
    return (double)NAN;
  }
  const double x = r / model->r0;
  double g = model->alpha * power(x, -model->m);
  if (model->k != 0.0) {
    g *= power(1.0 + power(x, model->n), -model->k);
  }
  return g;
}

myb_status myb_force_nongrav(const myb_nongrav* model, const double r[3], const double v[3],
                             double a[3]) {
  MYB_CHECK(model != NULL && r != NULL && v != NULL && a != NULL);
  const double rn = myb_vec3_norm(r);
  double h[3];
  myb_vec3_cross(r, v, h);
  const double hn = myb_vec3_norm(h);
  const double g = myb_nongrav_g(model, rn);
  if (!(rn > 0.0) || !(hn > 0.0) || !isfinite(g) || !isfinite(hn)) {
    return MYB_E_DOMAIN;
  }
  const double radial[3] = {r[0] / rn, r[1] / rn, r[2] / rn};
  const double normal[3] = {h[0] / hn, h[1] / hn, h[2] / hn};
  double transverse[3];
  myb_vec3_cross(normal, radial, transverse);
  for (int i = 0; i < 3; ++i) {
    a[i] += g * (model->a[0] * radial[i] + model->a[1] * transverse[i] +
                 model->a[2] * normal[i]);
  }
  return MYB_OK;
}

myb_status myb_force_nongrav_delayed(const myb_nongrav* model, double gm, double delay,
                                     const double r[3], const double v[3], double a[3]) {
  MYB_CHECK(model != NULL && r != NULL && v != NULL && a != NULL);
  if (delay == 0.0) {
    return myb_force_nongrav(model, r, v, a);
  }
  if (!isfinite(delay) || !(gm > 0.0) || !isfinite(gm)) {
    return MYB_E_DOMAIN;
  }
  myb_pv now;
  myb_pv past;
  memcpy(now.r, r, sizeof now.r);
  memcpy(now.v, v, sizeof now.v);
  const myb_status status = myb_kepler_propagate(gm, &now, -delay, &past);
  if (status != MYB_OK) {
    return status;
  }
  const double g_now = myb_nongrav_g(model, myb_vec3_norm(r));
  const double g_past = myb_nongrav_g(model, myb_vec3_norm(past.r));
  if (!isfinite(g_now) || !isfinite(g_past) || !(g_now > 0.0)) {
    return MYB_E_DOMAIN;
  }
  double present[3] = {0.0, 0.0, 0.0};
  const myb_status force = myb_force_nongrav(model, r, v, present);
  if (force != MYB_OK) {
    return force;
  }
  for (int i = 0; i < 3; ++i) {
    a[i] += present[i] * (g_past / g_now);
  }
  return MYB_OK;
}

/* The answer of a conducting sphere to a periodic heating. With T = c j1(k r) e^(i nu t)
 * inside, k^2 = -i nu / kappa, and the surface condition K dT/dr + 4 eps sigma T_av^3 T =
 * (heating), linearised about the mean temperature T_av = T_star / sqrt(2), the surface
 * amplitude is that of a body without inertia divided by 1 + lambda zeta(z). In the real
 * functions of x = sqrt(2) ratio (D. Vokrouhlicky, Astron. Astrophys. 335, 1093, 1998):
 *   E e^(i delta) = (A + i B) / (C + i D),  answer = E e^(i delta) / (1 + lambda),
 * every term divided here by e^x so that large bodies do not overflow. */
void myb_yarkovsky_response(double ratio, double theta, double* real, double* imaginary) {
  MYB_CHECK(real != NULL && imaginary != NULL);
  if (!(theta > 0.0)) {
    *real = 1.0;
    *imaginary = 0.0;
    return;
  }
  double re = 0.0; /* 1 + lambda zeta */
  double im = 0.0;
  if (!(ratio > 0.0) || ratio > 500.0) {
    re = 1.0 + 0.5 * theta;
    im = 0.5 * theta;
  } else {
    const double x = 1.4142135623730951 * ratio;
    const double lambda = theta / x;
    if (x < 0.4) {
      /* zeta = 1 - z^2/5 - z^4/175 - 2 z^6/7875 - 37 z^8/3031875 - 118 z^10/197071875
       *        - 5506 z^12/186232921875, z^2 = -i q, q = x^2 / 2: the closed form below loses
       * to cancellation what this gains, and they meet at 1e-12 */
      const double q = 0.5 * x * x;
      const double q2 = q * q;
      const double zeta_re =
          1.0 + q2 * (1.0 / 175.0 + q2 * (-37.0 / 3031875.0 + q2 * (5506.0 / 186232921875.0)));
      const double zeta_im =
          q * (1.0 / 5.0 + q2 * (-2.0 / 7875.0 + q2 * (118.0 / 197071875.0)));
      re = 1.0 + lambda * zeta_re;
      im = lambda * zeta_im;
    } else {
      const double e = myb_exp(-x);
      const double cx = myb_cos(x);
      const double sx = myb_sin(x);
      const double w = lambda / (1.0 + lambda);
      const double a = -(x + 2.0) * e - ((x - 2.0) * cx - x * sx);
      const double b = -x * e - (x * cx + (x - 2.0) * sx);
      const double cc = a + w * (3.0 * (x + 2.0) * e + (3.0 * (x - 2.0) * cx + x * (x - 3.0) * sx));
      const double d = b + w * (x * (x + 3.0) * e - (x * (x - 3.0) * cx - 3.0 * (x - 2.0) * sx));
      /* (A + i B) / ((C + i D) (1 + lambda)) */
      const double norm = (cc * cc + d * d) * (1.0 + lambda);
      *real = (a * cc + b * d) / norm;
      *imaginary = (b * cc - a * d) / norm;
      return;
    }
  }
  const double norm = re * re + im * im;
  *real = re / norm;
  *imaginary = -im / norm;
}

/* The thermal parameter of the body at a frequency and the radius in depths of the wave. */
static void thermal_state(const myb_yarkovsky* b, double flux, double frequency, double* theta,
                          double* ratio) {
  const double t4 = b->absorptivity * flux / (b->emissivity * stefan_boltzmann);
  const double t3 = power(t4, 0.75);
  const double root = sqrt(frequency);
  *theta = b->thermal_inertia * root / (b->emissivity * stefan_boltzmann * t3);
  *ratio = b->conductivity > 0.0 && b->thermal_inertia > 0.0 && root > 0.0
               ? b->radius * b->thermal_inertia * root / b->conductivity
               : 0.0;
}

static int valid_body(const myb_yarkovsky* b) {
  return b->radius > 0.0 && b->density > 0.0 && b->absorptivity > 0.0 &&
         b->absorptivity <= 1.0 && b->emissivity > 0.0 && b->emissivity <= 1.0 &&
         b->thermal_inertia >= 0.0 && b->rotation_rate >= 0.0 && b->conductivity >= 0.0 &&
         isfinite(b->radius) && isfinite(b->density) && isfinite(b->thermal_inertia) &&
         isfinite(b->rotation_rate) && isfinite(b->obliquity) && isfinite(b->conductivity);
}

myb_status myb_force_yarkovsky(const myb_yarkovsky* body, const double spin[3], double flux,
                               double c, const double to_sun[3], double a[3]) {
  MYB_CHECK(body != NULL && spin != NULL && to_sun != NULL && a != NULL);
  const double s2 = myb_vec3_dot(spin, spin);
  const double u2 = myb_vec3_dot(to_sun, to_sun);
  if (!valid_body(body) || !(flux > 0.0) || !isfinite(flux) || !(c > 0.0) ||
      !(fabs(s2 - 1.0) <= 1e-12) || !(fabs(u2 - 1.0) <= 1e-12)) {
    return MYB_E_DOMAIN;
  }
  double theta = 0.0;
  double ratio = 0.0;
  thermal_state(body, flux, body->rotation_rate, &theta, &ratio);
  double re = 1.0;
  double im = 0.0;
  myb_yarkovsky_response(ratio, theta, &re, &im);
  const double along = myb_vec3_dot(to_sun, spin);
  double across[3];
  double turned[3];
  for (int i = 0; i < 3; ++i) {
    across[i] = to_sun[i] - along * spin[i];
  }
  myb_vec3_cross(spin, across, turned);
  const double scale =
      (4.0 / 9.0) * body->absorptivity * 3.0 * flux / (4.0 * body->radius * body->density * c);
  for (int i = 0; i < 3; ++i) {
    a[i] -= scale * (re * across[i] - im * turned[i]);
  }
  return MYB_OK;
}

/* D. Vokrouhlicky, Astron. Astrophys. 335, 1093 (1998) and 344, 362 (1999), in the form of
 * W. F. Bottke et al., Annu. Rev. Earth Planet. Sci. 34, 157 (2006), eqs. 1-4, for a body much
 * larger than the depth the thermal wave reaches (kappa_1 = kappa_2 = kappa_3 = 1/2), or,
 * with the conductivity given, for any size (myb_yarkovsky_response). */
myb_status myb_yarkovsky_drift(const myb_yarkovsky* body, double flux, double mean_motion,
                               double c, double* diurnal, double* seasonal) {
  MYB_CHECK(body != NULL);
  const myb_yarkovsky* b = body;
  if (!(b->radius > 0.0) || !(b->density > 0.0) || !(b->absorptivity > 0.0) ||
      !(b->absorptivity <= 1.0) || !(b->emissivity > 0.0) || !(b->emissivity <= 1.0) ||
      !(b->thermal_inertia >= 0.0) || !(b->rotation_rate >= 0.0) || !(flux > 0.0) ||
      !(mean_motion > 0.0) || !(c > 0.0) || !isfinite(b->radius) || !isfinite(b->density) ||
      !isfinite(b->thermal_inertia) || !isfinite(b->rotation_rate) || !isfinite(b->obliquity) ||
      !isfinite(flux) || !isfinite(mean_motion) || !(b->conductivity >= 0.0) ||
      !isfinite(b->conductivity)) {
    return MYB_E_DOMAIN;
  }
  /* the force of the light per unit mass, the subsolar temperature, the thermal parameters */
  const double phi = 3.0 * flux / (4.0 * b->radius * b->density * c);
  const double t4 = b->absorptivity * flux / (b->emissivity * stefan_boltzmann);
  const double t3 = power(t4, 0.75);
  const double unit = b->thermal_inertia / (b->emissivity * stefan_boltzmann * t3);
  const double theta[2] = {unit * sqrt(b->rotation_rate), unit * sqrt(mean_motion)};
  double f[2];
  for (int k = 0; k < 2; ++k) {
    f[k] = -0.5 * theta[k] / (1.0 + theta[k] + 0.5 * theta[k] * theta[k]);
  }
  if (b->conductivity > 0.0 && b->thermal_inertia > 0.0) {
    const double frequency[2] = {b->rotation_rate, mean_motion};
    for (int k = 0; k < 2; ++k) {
      double re = 1.0;
      double ratio = 0.0;
      double unused = 0.0;
      thermal_state(b, flux, frequency[k], &unused, &ratio);
      myb_yarkovsky_response(ratio, theta[k], &re, &f[k]);
    }
  }
  const double cg = myb_cos(b->obliquity);
  const double sg = myb_sin(b->obliquity);
  const double scale = b->absorptivity * phi / mean_motion;
  if (diurnal != NULL) {
    *diurnal = -(8.0 / 9.0) * scale * f[0] * cg;
  }
  if (seasonal != NULL) {
    *seasonal = (4.0 / 9.0) * scale * f[1] * sg * sg;
  }
  return MYB_OK;
}

void myb_force_emission(double power_emitted, double efficiency, double mass, double c,
                        const double direction[3], double a[3]) {
  MYB_CHECK(direction != NULL && a != NULL && mass > 0.0 && c > 0.0);
  const double k = efficiency * power_emitted / (mass * c);
  for (int i = 0; i < 3; ++i) {
    a[i] -= k * direction[i];
  }
}

/* F. Mignard, Moon and Planets 20, 301 (1979), eq. 5, and 23, 185 (1980): the tide raised on
 * a body by a point mass, delayed by a constant time; the force on the raiser, turned into the
 * acceleration of its motion relative to the body by the factor 1 + m / M. */
myb_status myb_force_tidal_lag(double gm_body, double gm_raiser, double k2, double radius,
                               double lag, const double omega[3], const double r[3],
                               const double v[3], int with_static, double a[3]) {
  MYB_CHECK(omega != NULL && r != NULL && v != NULL && a != NULL);
  const double r2 = myb_vec3_dot(r, r);
  if (!(gm_body > 0.0) || !(gm_raiser >= 0.0) || !(radius > 0.0) || !(r2 > 0.0) ||
      !isfinite(gm_body) || !isfinite(gm_raiser) || !isfinite(k2) || !isfinite(radius) ||
      !isfinite(lag) || !isfinite(r2) || !myb_vec3_finite(omega) || !myb_vec3_finite(v)) {
    return MYB_E_DOMAIN;
  }
  const double radius2 = radius * radius;
  const double r4 = r2 * r2;
  /* 3 k2 G m (1 + m / M) R^5 / r^8 */
  const double scale =
      3.0 * k2 * gm_raiser * (1.0 + gm_raiser / gm_body) * (radius2 * radius2 * radius) / (r4 * r4);
  double spin[3];
  myb_vec3_cross(r, omega, spin);
  const double rv = myb_vec3_dot(r, v);
  for (int i = 0; i < 3; ++i) {
    const double lagging = lag * (2.0 * rv * r[i] / r2 + spin[i] + v[i]);
    a[i] -= scale * ((with_static ? r[i] : 0.0) + lagging);
  }
  return MYB_OK;
}

myb_status myb_torque_tidal_lag(double gm_body, double gm_raiser, double k2, double radius,
                                double lag, const double omega[3], const double r[3],
                                const double v[3], double torque[3]) {
  MYB_CHECK(torque != NULL);
  double a[3] = {0.0, 0.0, 0.0};
  const myb_status status =
      myb_force_tidal_lag(gm_body, gm_raiser, k2, radius, lag, omega, r, v, 0, a);
  if (status != MYB_OK) {
    return status;
  }
  /* the force on the raiser is its mass times a / (1 + m / M); the torque on the body,
   * -r x F, per unit of the body's mass */
  double lever[3];
  myb_vec3_cross(r, a, lever);
  const double scale = -gm_raiser / (gm_body + gm_raiser);
  for (int i = 0; i < 3; ++i) {
    torque[i] = scale * lever[i];
  }
  return MYB_OK;
}

/* ---- The seasonal effect on an eccentric orbit ---------------------------------------------- */

enum { SEASON_MOST_HARMONICS = 4096 };

struct myb_yarkovsky_season {
  myb_allocator allocator;
  double e;
  double p[3]; /* to the perihelion */
  double q[3]; /* across, in the plane, in the sense of the motion */
  double spin[3];
  double scale; /* (4 alpha / 9) of the light's force per unit mass at the semi-major axis */
  double drift;
  size_t harmonics;
  double* series; /* 2 a harmonic: the factors of cos(k M) and of sin(k M) */
};

void myb_yarkovsky_season_destroy(myb_yarkovsky_season* season) {
  if (season == NULL) {
    return;
  }
  const myb_allocator a = season->allocator;
  myb_internal_release(&a, season->series);
  myb_internal_release(&a, season);
}

/* The sum of the series at the mean anomaly M: cos(k M) and sin(k M) by turning the pair of
 * the harmonic before (the error grows as k, not as k^2). */
static double season_sum(const myb_yarkovsky_season* s, double mean_anomaly) {
  const double c1 = myb_cos(mean_anomaly);
  const double s1 = myb_sin(mean_anomaly);
  double ck = 1.0;
  double sk = 0.0;
  double sum = 0.0;
  for (size_t k = 0; k < s->harmonics; ++k) {
    const double next = ck * c1 - sk * s1;
    sk = sk * c1 + ck * s1;
    ck = next;
    sum += s->series[2 * k] * ck + s->series[2 * k + 1] * sk;
  }
  return sum;
}

myb_status myb_yarkovsky_season_create(const myb_yarkovsky* body, const double spin[3], double gm,
                                       const double state[6], double flux, double distance,
                                       double c, const myb_allocator* allocator,
                                       myb_yarkovsky_season** season) {
  MYB_CHECK(body != NULL && spin != NULL && state != NULL && season != NULL);
  static const double two_pi = 6.283185307179586;
  const double s2 = myb_vec3_dot(spin, spin);
  if (!valid_body(body) || !(fabs(s2 - 1.0) <= 1e-12) || !(gm > 0.0) || !isfinite(gm) ||
      !(flux > 0.0) || !isfinite(flux) || !(distance > 0.0) || !isfinite(distance) ||
      !(c > 0.0) || !myb_vec3_finite(state) || !myb_vec3_finite(state + 3)) {
    return MYB_E_DOMAIN;
  }
  /* the orbit: the semi-major axis, the eccentricity vector, the axes of the plane */
  const double* r = state;
  const double* v = state + 3;
  const double rn = myb_vec3_norm(r);
  const double v2 = myb_vec3_dot(v, v);
  double h[3];
  myb_vec3_cross(r, v, h);
  const double hn = myb_vec3_norm(h);
  if (!(rn > 0.0) || !(hn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double inverse_a = 2.0 / rn - v2 / gm;
  if (!(inverse_a > 0.0) || !isfinite(inverse_a)) {
    return MYB_E_DOMAIN;
  }
  const double a = 1.0 / inverse_a;
  double vh[3];
  double ev[3];
  myb_vec3_cross(v, h, vh);
  for (int i = 0; i < 3; ++i) {
    ev[i] = vh[i] / gm - r[i] / rn;
  }
  const double e = myb_vec3_norm(ev);
  if (!(e <= 0.95)) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator al = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&al, 1, sizeof(myb_yarkovsky_season), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_yarkovsky_season* s = block;
  memset(s, 0, sizeof *s);
  s->allocator = al;
  s->e = e;
  for (int i = 0; i < 3; ++i) {
    s->p[i] = e > 1e-12 ? ev[i] / e : r[i] / rn; /* a circle: from where the body is */
    s->spin[i] = spin[i];
  }
  const double unit_h[3] = {h[0] / hn, h[1] / hn, h[2] / hn};
  myb_vec3_cross(unit_h, s->p, s->q);
  const double sp = myb_vec3_dot(spin, s->p);
  const double sq = myb_vec3_dot(spin, s->q);
  const double root = sqrt((1.0 - e) * (1.0 + e));
  const double mean_motion = sqrt(gm / a) / a;
  const double flux_a = flux * (distance / a) * (distance / a);
  s->scale = (4.0 / 9.0) * body->absorptivity * 3.0 * flux_a /
             (4.0 * body->radius * body->density * c);
  /* the harmonics of (a/r)^2 cos f and (a/r)^2 sin f fall as exp(-k rho(e)) */
  size_t harmonics = 2;
  if (e > 1e-8) {
    const double rho = myb_log((1.0 + root) / e) - root;
    const double needed = 32.3 / rho + 3.0; /* 1e-14 */
    harmonics = needed < (double)SEASON_MOST_HARMONICS ? (size_t)needed : SEASON_MOST_HARMONICS;
  }
  s->harmonics = harmonics;
  /* the series, then room for the harmonics of the heating */
  status = myb_internal_alloc_array(&al, 4 * harmonics, sizeof(double), &block);
  if (status != MYB_OK) {
    myb_yarkovsky_season_destroy(s);
    return status;
  }
  s->series = block;
  double* heating = s->series + 2 * harmonics; /* a_k, b_k */
  for (size_t k = 0; k < 2 * harmonics; ++k) {
    heating[k] = 0.0;
  }
  /* a_k = (1/pi) int (a/r)^2 cos f cos(k M) dM, b_k likewise with the sines; in the
   * eccentric anomaly, dM = (r/a) dE: a sum over equal steps of E, exact to rounding for a
   * periodic integrand once the steps resolve it */
  const size_t nodes = 4 * harmonics + 64;
  for (size_t j = 0; j < nodes; ++j) {
    const double ea = two_pi * (double)j / (double)nodes;
    const double ce = myb_cos(ea);
    const double se = myb_sin(ea);
    const double w = 1.0 / ((1.0 - e * ce) * (1.0 - e * ce));
    const double ga = (ce - e) * w;
    const double gb = root * se * w;
    const double m = ea - e * se;
    const double c1 = myb_cos(m);
    const double s1 = myb_sin(m);
    double ck = 1.0;
    double sk = 0.0;
    for (size_t k = 0; k < harmonics; ++k) {
      const double next = ck * c1 - sk * s1;
      sk = sk * c1 + ck * s1;
      ck = next;
      heating[2 * k] += ga * ck;
      heating[2 * k + 1] += gb * sk;
    }
  }
  /* each harmonic answered at its frequency, about the mean temperature of the orbit; the
   * heating is -Phi (s . r/r): cos(k M) is answered with re cos - im sin, sin(k M) with
   * re sin + im cos */
  const double mean_flux = flux_a / root;
  for (size_t k = 0; k < harmonics; ++k) {
    const double ak = heating[2 * k] * 2.0 / (double)nodes;
    const double bk = heating[2 * k + 1] * 2.0 / (double)nodes;
    double theta = 0.0;
    double ratio = 0.0;
    thermal_state(body, mean_flux, (double)(k + 1) * mean_motion, &theta, &ratio);
    double re = 1.0;
    double im = 0.0;
    myb_yarkovsky_response(ratio, theta, &re, &im);
    s->series[2 * k] = sp * ak * re + sq * bk * im;
    s->series[2 * k + 1] = -sp * ak * im + sq * bk * re;
  }
  /* the mean drift: da/dt = 2 / (n sqrt(1 - e^2)) (e sin f a_R + (1 + e cos f) a_T) */
  double sum = 0.0;
  for (size_t j = 0; j < nodes; ++j) {
    const double ea = two_pi * (double)j / (double)nodes;
    const double ce = myb_cos(ea);
    const double se = myb_sin(ea);
    const double cf = (ce - e) / (1.0 - e * ce);
    const double sf = root * se / (1.0 - e * ce);
    const double along = s->scale * season_sum(s, ea - e * se);
    const double radial = along * (sp * cf + sq * sf);
    const double transverse = along * (-sp * sf + sq * cf);
    sum += (e * sf * radial + (1.0 + e * cf) * transverse) * (1.0 - e * ce);
  }
  s->drift = 2.0 / (mean_motion * root) * sum / (double)nodes;
  *season = s;
  return MYB_OK;
}

void myb_force_yarkovsky_season(const myb_yarkovsky_season* season, const double r[3],
                                double a[3]) {
  MYB_CHECK(season != NULL && r != NULL && a != NULL);
  const myb_yarkovsky_season* s = season;
  /* the true anomaly of the direction of r in the plane, then the mean anomaly */
  const double x = myb_vec3_dot(r, s->p);
  const double y = myb_vec3_dot(r, s->q);
  const double rho = myb_hypot(x, y);
  if (!(rho > 0.0) || !isfinite(rho)) {
    return;
  }
  const double cf = x / rho;
  const double sf = y / rho;
  const double root = sqrt((1.0 - s->e) * (1.0 + s->e));
  const double ea = myb_atan2(root * sf, s->e + cf);
  const double along = s->scale * season_sum(s, ea - s->e * myb_sin(ea));
  for (int i = 0; i < 3; ++i) {
    a[i] += along * s->spin[i];
  }
}

double myb_yarkovsky_season_drift(const myb_yarkovsky_season* season) {
  MYB_CHECK(season != NULL);
  return season->drift;
}

size_t myb_yarkovsky_season_harmonics(const myb_yarkovsky_season* season) {
  MYB_CHECK(season != NULL);
  return season->harmonics;
}
