/* SPDX-License-Identifier: BSD-3-Clause */
/* Cartesian state <-> classical and modified equinoctial elements.
 *
 * Angles come from atan2 of vector components (never acos), which keeps them accurate over
 * the whole range. Classical: Vallado (2022) alg. 9/10, Montenbruck & Gill (2000) §2.2;
 * equinoctial frame: Broucke & Cefola (1972); modified equinoctial elements: Walker, Ireland
 * & Owens (1985). */
#include <math.h> /* sqrt, fabs */

#include "core/base.h"
#include "kepler/stumpff.h"
#include "myb/cmath.h"
#include "myb/kepler.h"
#include "num/vec3.h"

static const double pi_hi = 0x1.921fb54442d18p+1;

static int finite(double x) {
  return x - x == 0.0;
}

static double positive_angle(double a) { /* (-pi, pi] -> [0, 2 pi) */
  return a < 0.0 ? myb_internal_add_turns(a, 1.0) : a;
}

static double principal(double a) { /* any -> (-pi, pi] */
  double turns = 0.0;
  const double r = myb_internal_reduce_angle(a, &turns);
  return r <= -pi_hi ? pi_hi : r;
}

/* Angular momentum h, eccentricity vector, |r|; MYB_E_DOMAIN for invalid input. */
typedef struct orbit_vectors {
  double h[3];
  double h_norm;
  double e[3];
  double r_norm;
} orbit_vectors;

static myb_status orbit_vectors_from(double mu, const myb_pv* pv, orbit_vectors* o) {
  if (!(mu > 0.0) || !finite(mu) || !myb_vec3_finite(pv->r) || !myb_vec3_finite(pv->v)) {
    return MYB_E_DOMAIN;
  }
  o->r_norm = myb_vec3_norm(pv->r);
  myb_vec3_cross(pv->r, pv->v, o->h);
  o->h_norm = myb_vec3_norm(o->h);
  if (!(o->r_norm > 0.0) || !(o->h_norm > 0.0)) {
    return MYB_E_DOMAIN;
  }
  /* e = (v x h) / mu - r / |r| */
  double vxh[3];
  myb_vec3_cross(pv->v, o->h, vxh);
  myb_vec3_combine(vxh, 1.0 / mu, pv->r, -1.0 / o->r_norm, o->e);
  return MYB_OK;
}

myb_status myb_elements_from_pv(double mu, const myb_pv* pv, myb_elements* elements) {
  MYB_CHECK(pv != NULL && elements != NULL);
  orbit_vectors o;
  const myb_status status = orbit_vectors_from(mu, pv, &o);
  if (status != MYB_OK) {
    return status;
  }
  myb_elements el = {};
  el.p = o.h_norm * o.h_norm / mu;
  el.e = myb_vec3_norm(o.e);
  const double h_xy = myb_hypot(o.h[0], o.h[1]);
  el.i = myb_atan2(h_xy, o.h[2]);

  /* Node vector n = z x h = (-h_y, h_x, 0); for an equatorial orbit the x axis replaces it. */
  const int equatorial = h_xy == 0.0;
  const double n[3] = {equatorial ? 1.0 : -o.h[1], equatorial ? 0.0 : o.h[0], 0.0};
  el.raan = equatorial ? 0.0 : positive_angle(myb_atan2(o.h[0], -o.h[1]));

  /* Angles in the orbital plane, measured in the direction of motion:
   * angle(a -> b) = atan2(h . (a x b) / |h|, a . b). */
  double cross[3];
  if (el.e == 0.0) {
    el.argp = 0.0;
    myb_vec3_cross(n, pv->r, cross);
    el.nu = principal(myb_atan2(myb_vec3_dot(o.h, cross) / o.h_norm, myb_vec3_dot(n, pv->r)));
  } else {
    myb_vec3_cross(n, o.e, cross);
    el.argp = positive_angle(myb_atan2(myb_vec3_dot(o.h, cross) / o.h_norm, myb_vec3_dot(n, o.e)));
    myb_vec3_cross(o.e, pv->r, cross);
    el.nu = principal(myb_atan2(myb_vec3_dot(o.h, cross) / o.h_norm, myb_vec3_dot(o.e, pv->r)));
  }
  *elements = el;
  return MYB_OK;
}

myb_status myb_pv_from_elements(double mu, const myb_elements* elements, myb_pv* pv) {
  MYB_CHECK(elements != NULL && pv != NULL);
  const myb_elements el = *elements;
  if (!(mu > 0.0) || !finite(mu) || !(el.p > 0.0) || !finite(el.p) || !(el.e >= 0.0) ||
      !finite(el.e) || !finite(el.i) || !finite(el.raan) || !finite(el.argp) || !finite(el.nu)) {
    return MYB_E_DOMAIN;
  }
  const double cos_nu = myb_cos(el.nu);
  const double sin_nu = myb_sin(el.nu);
  const double w = 1.0 + el.e * cos_nu;
  if (!(w > 0.0)) {
    return MYB_E_DOMAIN; /* beyond the asymptotes of a hyperbola */
  }
  const double r = el.p / w;
  const double co = myb_cos(el.raan);
  const double so = myb_sin(el.raan);
  const double cw = myb_cos(el.argp);
  const double sw = myb_sin(el.argp);
  const double ci = myb_cos(el.i);
  const double si = myb_sin(el.i);
  /* Perifocal unit vectors P (to periapsis) and Q (90 degrees ahead). */
  const double p_hat[3] = {co * cw - so * sw * ci, so * cw + co * sw * ci, sw * si};
  const double q_hat[3] = {-co * sw - so * cw * ci, -so * sw + co * cw * ci, cw * si};
  const double speed = sqrt(mu / el.p);
  myb_vec3_combine(p_hat, r * cos_nu, q_hat, r * sin_nu, pv->r);
  myb_vec3_combine(p_hat, -speed * sin_nu, q_hat, speed * (el.e + cos_nu), pv->v);
  return MYB_OK;
}

/* Equinoctial frame (Broucke & Cefola): f_hat towards the reference direction, g_hat 90
 * degrees ahead in the orbital plane, both from h = tan(i/2) cos(raan), k = tan(i/2) sin(raan). */
static void equinoctial_frame(double h, double k, double f_hat[3], double g_hat[3]) {
  const double s2 = 1.0 + h * h + k * k;
  f_hat[0] = (1.0 - k * k + h * h) / s2;
  f_hat[1] = 2.0 * h * k / s2;
  f_hat[2] = -2.0 * k / s2;
  g_hat[0] = 2.0 * h * k / s2;
  g_hat[1] = (1.0 + k * k - h * h) / s2;
  g_hat[2] = 2.0 * h / s2;
}

myb_status myb_equinoctial_from_pv(double mu, const myb_pv* pv, myb_equinoctial* elements) {
  MYB_CHECK(pv != NULL && elements != NULL);
  orbit_vectors o;
  const myb_status status = orbit_vectors_from(mu, pv, &o);
  if (status != MYB_OK) {
    return status;
  }
  const double unit[3] = {o.h[0] / o.h_norm, o.h[1] / o.h_norm, o.h[2] / o.h_norm};
  const double denominator = 1.0 + unit[2];
  if (!(denominator > 0.0)) {
    return MYB_E_DOMAIN; /* i = pi: h and k are infinite */
  }
  myb_equinoctial el = {};
  el.p = o.h_norm * o.h_norm / mu;
  el.h = -unit[1] / denominator;
  el.k = unit[0] / denominator;
  double f_hat[3];
  double g_hat[3];
  equinoctial_frame(el.h, el.k, f_hat, g_hat);
  el.f = myb_vec3_dot(o.e, f_hat);
  el.g = myb_vec3_dot(o.e, g_hat);
  el.l = principal(myb_atan2(myb_vec3_dot(pv->r, g_hat), myb_vec3_dot(pv->r, f_hat)));
  *elements = el;
  return MYB_OK;
}

myb_status myb_pv_from_equinoctial(double mu, const myb_equinoctial* elements, myb_pv* pv) {
  MYB_CHECK(elements != NULL && pv != NULL);
  const myb_equinoctial el = *elements;
  if (!(mu > 0.0) || !finite(mu) || !(el.p > 0.0) || !finite(el.p) || !finite(el.f) ||
      !finite(el.g) || !finite(el.h) || !finite(el.k) || !finite(el.l)) {
    return MYB_E_DOMAIN;
  }
  const double cl = myb_cos(el.l);
  const double sl = myb_sin(el.l);
  const double w = 1.0 + el.f * cl + el.g * sl;
  if (!(w > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double r = el.p / w;
  const double speed = sqrt(mu / el.p);
  double f_hat[3];
  double g_hat[3];
  equinoctial_frame(el.h, el.k, f_hat, g_hat);
  myb_vec3_combine(f_hat, r * cl, g_hat, r * sl, pv->r);
  myb_vec3_combine(f_hat, -speed * (el.g + sl), g_hat, speed * (el.f + cl), pv->v);
  return MYB_OK;
}
