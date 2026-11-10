/* SPDX-License-Identifier: BSD-3-Clause */
/* Continuous thrust: burn arcs, steering, Gauss equations, Edelbaum. */
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "core/check.h"
#include "myb/cmath.h"
#include "myb/thrust.h"
#include "num/vec3.h"

static const double PI = 3.141592653589793;

myb_status myb_thrust_init(myb_thrust* thrust, double mass0, const myb_thrust_arc* arcs,
                           size_t count) {
  MYB_CHECK(thrust != NULL && (arcs != NULL || count == 0));
  if (!(mass0 > 0.0) || !isfinite(mass0) || count > MYB_THRUST_MAX_ARCS) {
    return MYB_E_DOMAIN;
  }
  memset(thrust, 0, sizeof *thrust);
  thrust->mass0 = mass0;
  double used = 0.0;
  for (size_t i = 0; i < count; ++i) {
    myb_thrust_arc a = arcs[i];
    if (!(a.stop > a.start) || !isfinite(a.start) || !isfinite(a.stop) || !(a.thrust > 0.0) ||
        !(a.exhaust > 0.0) || !isfinite(a.thrust) || !isfinite(a.exhaust) ||
        (i > 0 && a.start < arcs[i - 1].stop)) {
      return MYB_E_DOMAIN;
    }
    if (a.steering == MYB_STEER_INERTIAL || a.steering == MYB_STEER_RTN) {
      const double n = myb_vec3_norm(a.direction);
      if (!(n > 0.0) || !isfinite(n)) {
        return MYB_E_DOMAIN;
      }
      for (int k = 0; k < 3; ++k) {
        a.direction[k] /= n;
      }
    } else if (a.steering == MYB_STEER_EDELBAUM) {
      if (!(a.edelbaum_v0 > 0.0) || !isfinite(a.edelbaum_beta0) ||
          (a.edelbaum_sense != 1.0 && a.edelbaum_sense != -1.0)) {
        return MYB_E_DOMAIN;
      }
    } else if (a.steering != MYB_STEER_VELOCITY && a.steering != MYB_STEER_ANTI_VELOCITY) {
      return MYB_E_DOMAIN;
    }
    used += a.thrust / a.exhaust * (a.stop - a.start);
    thrust->arcs[i] = a;
  }
  if (!(used < mass0)) {
    return MYB_E_DOMAIN;
  }
  thrust->count = count;
  return MYB_OK;
}

double myb_thrust_mass(const myb_thrust* thrust, double t) {
  MYB_CHECK(thrust != NULL);
  double m = thrust->mass0;
  for (size_t i = 0; i < thrust->count; ++i) {
    const myb_thrust_arc* a = &thrust->arcs[i];
    if (t <= a->start) {
      break;
    }
    const double burnt = (t < a->stop ? t : a->stop) - a->start;
    m -= a->thrust / a->exhaust * burnt;
  }
  return m;
}

double myb_thrust_mass_at(void* thrust, double t) {
  return myb_thrust_mass(thrust, t);
}

/* The unit thrust direction in the frame of x for arc a at t (mass m). */
static myb_status direction(const myb_thrust_arc* a, double t, double m, double mass_start,
                            const double x[3], const double v[3], double d[3]) {
  (void)t;
  double rhat[3];
  double nhat[3];
  double that[3];
  double h[3];
  const double r = myb_vec3_norm(x);
  myb_vec3_cross(x, v, h);
  const double hn = myb_vec3_norm(h);
  switch (a->steering) {
  case MYB_STEER_VELOCITY:
  case MYB_STEER_ANTI_VELOCITY: {
    const double vn = myb_vec3_norm(v);
    if (!(vn > 0.0)) {
      return MYB_E_DOMAIN;
    }
    const double s = a->steering == MYB_STEER_VELOCITY ? 1.0 : -1.0;
    for (int k = 0; k < 3; ++k) {
      d[k] = s * v[k] / vn;
    }
    return MYB_OK;
  }
  case MYB_STEER_INERTIAL:
    memcpy(d, a->direction, 3 * sizeof *d);
    return MYB_OK;
  case MYB_STEER_RTN:
  case MYB_STEER_EDELBAUM:
    break;
  }
  if (!(r > 0.0) || !(hn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  for (int k = 0; k < 3; ++k) {
    rhat[k] = x[k] / r;
    nhat[k] = h[k] / hn;
  }
  myb_vec3_cross(nhat, rhat, that);
  if (a->steering == MYB_STEER_RTN) {
    for (int k = 0; k < 3; ++k) {
      d[k] = a->direction[0] * rhat[k] + a->direction[1] * that[k] + a->direction[2] * nhat[k];
    }
    return MYB_OK;
  }
  /* Edelbaum: yaw from the velocity increment of the arc so far */
  const double dv = a->exhaust * myb_log(mass_start / m);
  const double sb0 = myb_sin(a->edelbaum_beta0);
  const double cb0 = myb_cos(a->edelbaum_beta0);
  const double beta = myb_atan2(a->edelbaum_v0 * sb0, a->edelbaum_v0 * cb0 - dv);
  /* cos u: x against the ascending node (z x h; the x axis for an equatorial orbit) */
  double node[3] = {-h[1], h[0], 0.0};
  double nn = myb_vec3_norm(node);
  if (!(nn > 0.0)) {
    node[0] = 1.0;
    node[1] = 0.0;
    nn = 1.0;
  }
  const double cosu = myb_vec3_dot(x, node) / (r * nn);
  const double side = cosu >= 0.0 ? 1.0 : -1.0;
  const double vn = myb_vec3_norm(v);
  const double cb = myb_cos(beta);
  const double sb = myb_sin(beta);
  for (int k = 0; k < 3; ++k) {
    d[k] = cb * v[k] / vn + sb * a->edelbaum_sense * side * nhat[k];
  }
  return MYB_OK;
}

myb_status myb_thrust_acceleration(void* thrust, double t, const double x[3], const double v[3],
                                   double a[3]) {
  MYB_CHECK(thrust != NULL && x != NULL && v != NULL && a != NULL);
  const myb_thrust* th = thrust;
  a[0] = a[1] = a[2] = 0.0;
  double mass_start = th->mass0;
  for (size_t i = 0; i < th->count; ++i) {
    const myb_thrust_arc* arc = &th->arcs[i];
    if (t < arc->start) {
      break;
    }
    if (t < arc->stop) {
      const double m = mass_start - arc->thrust / arc->exhaust * (t - arc->start);
      double d[3];
      const myb_status st = direction(arc, t, m, mass_start, x, v, d);
      if (st != MYB_OK) {
        return st;
      }
      const double f = arc->thrust / m;
      for (int k = 0; k < 3; ++k) {
        a[k] = f * d[k];
      }
      return MYB_OK;
    }
    mass_start -= arc->thrust / arc->exhaust * (arc->stop - arc->start);
  }
  return MYB_OK;
}

myb_status myb_thrust_rtn(const myb_thrust* thrust, double t, const double x[3],
                          const double v[3], double a_rtn[3]) {
  double a[3];
  const myb_status st = myb_thrust_acceleration((void*)(uintptr_t)thrust, t, x, v, a);
  if (st != MYB_OK) {
    return st;
  }
  double h[3];
  myb_vec3_cross(x, v, h);
  const double r = myb_vec3_norm(x);
  const double hn = myb_vec3_norm(h);
  if (!(r > 0.0) || !(hn > 0.0)) {
    return MYB_E_DOMAIN;
  }
  double rhat[3];
  double nhat[3];
  double that[3];
  for (int k = 0; k < 3; ++k) {
    rhat[k] = x[k] / r;
    nhat[k] = h[k] / hn;
  }
  myb_vec3_cross(nhat, rhat, that);
  a_rtn[0] = myb_vec3_dot(a, rhat);
  a_rtn[1] = myb_vec3_dot(a, that);
  a_rtn[2] = myb_vec3_dot(a, nhat);
  return MYB_OK;
}

void myb_equinoctial_rates(double mu, const myb_equinoctial* e, const double a_rtn[3],
                           double rates[6]) {
  MYB_CHECK(e != NULL && a_rtn != NULL && rates != NULL);
  const double cl = myb_cos(e->l);
  const double sl = myb_sin(e->l);
  const double w = 1.0 + e->f * cl + e->g * sl;
  const double s2 = 1.0 + e->h * e->h + e->k * e->k;
  const double q = sqrt(e->p / mu);
  const double z = e->h * sl - e->k * cl;
  const double ar = a_rtn[0];
  const double at = a_rtn[1];
  const double an = a_rtn[2];
  rates[0] = 2.0 * e->p * q * at / w;
  rates[1] = q * (ar * sl + ((w + 1.0) * cl + e->f) * at / w - z * e->g * an / w);
  rates[2] = q * (-ar * cl + ((w + 1.0) * sl + e->g) * at / w + z * e->f * an / w);
  rates[3] = q * s2 * cl * an / (2.0 * w);
  rates[4] = q * s2 * sl * an / (2.0 * w);
  rates[5] = sqrt(mu * e->p) * (w / e->p) * (w / e->p) + q * z * an / w;
}

myb_status myb_edelbaum(double v0, double v1, double di, double* dv, double* beta0) {
  MYB_CHECK(dv != NULL && beta0 != NULL);
  if (!(v0 > 0.0) || !(v1 > 0.0) || !isfinite(v0) || !isfinite(v1) || !isfinite(di)) {
    return MYB_E_DOMAIN;
  }
  const double c = myb_cos(0.5 * PI * di);
  const double s = myb_sin(0.5 * PI * di);
  /* v0^2 - 2 v0 v1 c + v1^2 = (v0 - v1)^2 + 2 v0 v1 (1 - c), 1 - c = 2 sin^2(pi di / 4) */
  const double h = myb_sin(0.25 * PI * di);
  *dv = sqrt((v0 - v1) * (v0 - v1) + 4.0 * v0 * v1 * h * h);
  *beta0 = myb_atan2(s, v0 / v1 - c);
  return MYB_OK;
}

void myb_edelbaum_state(double v0, double beta0, double dv, double* v, double* di) {
  const double cb = myb_cos(beta0);
  const double sb = myb_sin(beta0);
  if (v != NULL) {
    *v = sqrt(v0 * v0 - 2.0 * v0 * dv * cb + dv * dv);
  }
  if (di != NULL) {
    *di = 2.0 / PI * (myb_atan((dv - v0 * cb) / (v0 * sb)) + 0.5 * PI - beta0);
  }
}
