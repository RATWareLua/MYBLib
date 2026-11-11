/* SPDX-License-Identifier: BSD-3-Clause */
/* Measurements with the light time (myb/estim.h): what a station receives at t left the
 * object a light time earlier, and, two-way, left the station two light times earlier.
 *
 * The estimator gives a measurement the state y of the object at the time of the
 * observation. Where the object was when the light left it follows from y by its motion
 * about a point mass at the origin over the light time tau, to the third order:
 *   p = x - v tau + a tau^2 / 2 - j tau^3 / 6,   v_e = v - a tau + j tau^2 / 2,
 * a = -gm x / r^3 and j its rate along the motion. What is left out is tau^4 of the fourth
 * derivative: 1e-14 m for a low orbit of the Earth, 1e-8 m at the distance of the Moon.
 *
 * Every quantity is carried with its derivatives in y (num/dual.h); the light-time equation
 * c tau = |p(tau) - s| is solved by fixed point on such numbers, so the derivatives of tau
 * are those of the implicit function. */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "myb/cmath.h"
#include "myb/estim.h"
#include "num/dual.h"

enum { ITERATIONS = 24 };

static const double two_pi = 6.283185307179586;

typedef struct light_path {
  myb_dual down;        /* the light time of the leg to the station, s */
  myb_dual up;          /* of the leg from the station (two-way) */
  myb_dual down_rate;   /* their rates in the time of reception */
  myb_dual up_rate;
  myb_dual d[3];        /* from the station at t to the object when the light left it */
  double station[6];    /* at t */
  myb_orientation axes; /* of the station's body at t */
} light_path;

static myb_dual dual_atan2(myb_dual y, myb_dual x) {
  myb_dual r;
  r.v = myb_atan2(y.v, x.v);
  const double h2 = x.v * x.v + y.v * y.v;
  for (int i = 0; i < MYB_DUAL_VARIABLES; ++i) {
    r.d[i] = (x.v * y.d[i] - y.v * x.d[i]) / h2;
  }
  return r;
}

/* The fixed point tau = |p(tau) - e(tau)| / c from `tau`: p the place of the object, e that
 * of the other end, both written by `ends`. */
typedef struct leg {
  const myb_dual* x;    /* the object at t: place, velocity, acceleration, its rate */
  const myb_dual* v;
  const myb_dual* a;
  const myb_dual* j;
  const double* other;  /* the station: its state at the instant `shift` before t */
  myb_dual lead;        /* how long before `other` the station is taken, apart from tau: the
                           down leg has the station at t (0), the up leg tau_down earlier */
  double shift;
  int up;               /* the station is taken tau itself earlier still */
} leg;

/* p = x - v T + a T^2 / 2 - j T^3 / 6 at T before t, and the place of the station's end. */
static void ends(const leg* g, myb_dual back, myb_dual tau, myb_dual p[3], myb_dual e[3]) {
  const myb_dual t2 = myb_dual_scale(myb_dual_mul(back, back), 0.5);
  const myb_dual t3 = myb_dual_scale(myb_dual_mul(t2, back), 1.0 / 3.0);
  for (int i = 0; i < 3; ++i) {
    p[i] = myb_dual_sub(myb_dual_add(myb_dual_sub(g->x[i], myb_dual_mul(g->v[i], back)),
                                     myb_dual_mul(g->a[i], t2)),
                        myb_dual_mul(g->j[i], t3));
    /* the station moves in a straight line over what is left of the shift: microseconds */
    myb_dual ahead = myb_dual_sub(myb_dual_constant(g->shift), g->lead);
    if (g->up) {
      ahead = myb_dual_sub(ahead, tau);
    }
    e[i] = myb_dual_add(myb_dual_constant(g->other[i]), myb_dual_scale(ahead, g->other[3 + i]));
  }
}

/* Solves a leg: `back` is how long before t the object is taken (tau itself for the down
 * leg, the given tau_down for the up leg). */
static myb_status solve_leg(const leg* g, const myb_dual* fixed_back, double c, myb_dual* tau,
                            myb_dual p[3], myb_dual e[3]) {
  int settled = 0;
  /* bound: ITERATIONS; the error falls by the speed of the ends over that of light a turn */
  for (int k = 0; k < ITERATIONS; ++k) {
    ends(g, fixed_back != NULL ? *fixed_back : *tau, *tau, p, e);
    const myb_dual between[3] = {myb_dual_sub(p[0], e[0]), myb_dual_sub(p[1], e[1]),
                                 myb_dual_sub(p[2], e[2])};
    const myb_dual length = myb_dual_norm(between);
    if (!(length.v > 0.0) || !isfinite(length.v)) {
      return MYB_E_DOMAIN;
    }
    const myb_dual next = myb_dual_scale(length, 1.0 / c);
    const double moved = fabs(next.v - tau->v);
    *tau = next;
    /* two more turns after the value has settled bring its derivatives to the same */
    if (moved <= 4e-16 * next.v) {
      if (++settled == 3) {
        ends(g, fixed_back != NULL ? *fixed_back : *tau, *tau, p, e);
        return MYB_OK;
      }
    } else {
      settled = 0;
    }
  }
  return MYB_E_NOCONV;
}

static myb_status solve(const myb_light_measure* m, double t, const double y[6],
                        light_path* path) {
  MYB_CHECK(m != NULL && m->station != NULL && m->station->env != NULL && y != NULL);
  if (!(m->c > 0.0) || !isfinite(m->c) || !(m->gm >= 0.0) || !isfinite(m->gm)) {
    return MYB_E_DOMAIN;
  }
  myb_status status = myb_station_state(m->station, t, path->station);
  if (status == MYB_OK) {
    status = myb_env_orientation(m->station->env, t, m->station->body, &path->axes);
  }
  if (status != MYB_OK) {
    return status;
  }
  myb_dual x[3];
  myb_dual v[3];
  myb_dual a[3];
  myb_dual j[3];
  myb_dual_state(y, y + 3, NULL, x, v);
  for (int i = 0; i < 3; ++i) {
    a[i] = myb_dual_constant(0.0);
    j[i] = myb_dual_constant(0.0);
  }
  if (m->gm > 0.0) {
    /* a = -gm x / r^3, j = -gm (v / r^3 - 3 x (x.v) / r^5) */
    const myb_dual r2 = myb_dual_dot(x, x);
    const myb_dual r = myb_dual_sqrt(r2);
    if (!(r.v > 0.0) || !isfinite(r.v)) {
      return MYB_E_DOMAIN;
    }
    const myb_dual k3 = myb_dual_div(myb_dual_constant(-m->gm), myb_dual_mul(r2, r));
    const myb_dual along = myb_dual_scale(myb_dual_div(myb_dual_dot(x, v), r2), 3.0);
    for (int i = 0; i < 3; ++i) {
      a[i] = myb_dual_mul(x[i], k3);
      j[i] = myb_dual_mul(myb_dual_sub(v[i], myb_dual_mul(x[i], along)), k3);
    }
  }
  const double from[3] = {y[0] - path->station[0], y[1] - path->station[1],
                          y[2] - path->station[2]};
  const double first = sqrt(from[0] * from[0] + from[1] * from[1] + from[2] * from[2]) / m->c;
  if (!(first > 0.0) || !isfinite(first)) {
    return MYB_E_DOMAIN;
  }
  /* the leg to the station */
  const leg down = {.x = x, .v = v, .a = a, .j = j, .other = path->station,
                    .lead = myb_dual_constant(0.0)};
  myb_dual p[3];
  myb_dual e[3];
  path->down = myb_dual_constant(first);
  status = solve_leg(&down, NULL, m->c, &path->down, p, e);
  if (status != MYB_OK) {
    return status;
  }
  /* the velocity of the object then, and the direction station -> object */
  myb_dual ve[3];
  const myb_dual half = myb_dual_scale(myb_dual_mul(path->down, path->down), 0.5);
  for (int i = 0; i < 3; ++i) {
    ve[i] = myb_dual_add(myb_dual_sub(v[i], myb_dual_mul(a[i], path->down)),
                         myb_dual_mul(j[i], half));
    path->d[i] = myb_dual_sub(p[i], e[i]);
  }
  const myb_dual range = myb_dual_norm(path->d);
  myb_dual u[3];
  myb_dual closing[3];
  for (int i = 0; i < 3; ++i) {
    u[i] = myb_dual_div(path->d[i], range);
    closing[i] = myb_dual_sub(ve[i], myb_dual_constant(path->station[3 + i]));
  }
  /* d tau / dt = u.(v_e - v_s) / (c + u.v_e) */
  path->down_rate = myb_dual_div(myb_dual_dot(u, closing),
                                 myb_dual_add(myb_dual_constant(m->c), myb_dual_dot(u, ve)));
  path->up = myb_dual_constant(0.0);
  path->up_rate = myb_dual_constant(0.0);
  if (!m->two_way) {
    return MYB_OK;
  }
  /* the leg from the station: the station two light times before t, then where it was when
   * the light that met the object left it */
  const double shift = 2.0 * path->down.v;
  double earlier[6];
  status = myb_station_state(m->station, t - shift, earlier);
  if (status != MYB_OK) {
    return status;
  }
  const leg up = {.x = x, .v = v, .a = a, .j = j, .other = earlier, .lead = path->down,
                  .shift = shift, .up = 1};
  path->up = path->down;
  status = solve_leg(&up, &path->down, m->c, &path->up, p, e);
  if (status != MYB_OK) {
    return status;
  }
  myb_dual w[3];
  for (int i = 0; i < 3; ++i) {
    w[i] = myb_dual_sub(p[i], e[i]);
  }
  const myb_dual reach = myb_dual_norm(w);
  /* the velocity of the station when it sent the light: that of `earlier` carried over what
   * is left of the shift by the mean acceleration of the station between there and t — the
   * instant depends on the state through the light times, and the velocity with it */
  const myb_dual ahead =
      myb_dual_sub(myb_dual_sub(myb_dual_constant(shift), path->down), path->up);
  myb_dual sent[3];
  myb_dual about[3];
  for (int i = 0; i < 3; ++i) {
    const double acceleration = (path->station[3 + i] - earlier[3 + i]) / shift;
    w[i] = myb_dual_div(w[i], reach);
    sent[i] = myb_dual_add(myb_dual_constant(earlier[3 + i]), myb_dual_scale(ahead, acceleration));
    about[i] = myb_dual_sub(ve[i], sent[i]);
  }
  const myb_dual station_part = myb_dual_dot(w, sent);
  /* d tau_up / dt = (1 - d tau_down / dt) w.(v_e - v_s') / (c - w.v_s') */
  path->up_rate = myb_dual_div(
      myb_dual_mul(myb_dual_sub(myb_dual_constant(1.0), path->down_rate),
                   myb_dual_dot(w, about)),
      myb_dual_sub(myb_dual_constant(m->c), station_part));
  return MYB_OK;
}

static void put(myb_dual value, double* out, double* partials) {
  *out = value.v;
  if (partials != NULL) {
    memcpy(partials, value.d, sizeof value.d);
  }
}

myb_status myb_measure_light_range(void* light, double t, const double y[6], double* values,
                                   double* partials) {
  const myb_light_measure* m = light;
  MYB_CHECK(values != NULL);
  light_path path;
  const myb_status status = solve(m, t, y, &path);
  if (status != MYB_OK) {
    return status;
  }
  const myb_dual range = m->two_way
                             ? myb_dual_scale(myb_dual_add(path.down, path.up), 0.5 * m->c)
                             : myb_dual_scale(path.down, m->c);
  put(range, values, partials);
  return MYB_OK;
}

myb_status myb_measure_light_range_rate(void* light, double t, const double y[6],
                                        double* values, double* partials) {
  const myb_light_measure* m = light;
  MYB_CHECK(values != NULL);
  light_path path;
  const myb_status status = solve(m, t, y, &path);
  if (status != MYB_OK) {
    return status;
  }
  const myb_dual rate =
      m->two_way ? myb_dual_scale(myb_dual_add(path.down_rate, path.up_rate), 0.5 * m->c)
                 : myb_dual_scale(path.down_rate, m->c);
  put(rate, values, partials);
  return MYB_OK;
}

/* The direction the light comes from: station -> object when the light left it, and, with
 * the aberration, as the moving station sees it (the Lorentz transformation of
 * myb_aberration, not normalised). */
static void direction(const myb_light_measure* m, const light_path* path, myb_dual out[3]) {
  if (!m->aberration) {
    memcpy(out, path->d, 3 * sizeof(myb_dual));
    return;
  }
  const myb_dual length = myb_dual_norm(path->d);
  const double b[3] = {path->station[3] / m->c, path->station[4] / m->c,
                       path->station[5] / m->c};
  const double root = sqrt(1.0 - (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]));
  myb_dual unit[3];
  myb_dual along = myb_dual_constant(0.0);
  for (int i = 0; i < 3; ++i) {
    unit[i] = myb_dual_div(path->d[i], length);
    along = myb_dual_add(along, myb_dual_scale(unit[i], b[i]));
  }
  const myb_dual factor =
      myb_dual_add(myb_dual_constant(1.0), myb_dual_scale(along, 1.0 / (1.0 + root)));
  for (int i = 0; i < 3; ++i) {
    out[i] = myb_dual_add(myb_dual_scale(unit[i], root), myb_dual_scale(factor, b[i]));
  }
}

/* Longitude in [0, 2 pi) and latitude of a direction. MYB_E_DOMAIN on the pole. */
static myb_status angles(const myb_dual p[3], double* values, double* partials) {
  const myb_dual h =
      myb_dual_sqrt(myb_dual_add(myb_dual_mul(p[0], p[0]), myb_dual_mul(p[1], p[1])));
  if (!(h.v > 0.0) || !isfinite(h.v)) {
    return MYB_E_DOMAIN;
  }
  myb_dual longitude = dual_atan2(p[1], p[0]);
  if (longitude.v < 0.0) {
    longitude.v += two_pi;
  }
  const myb_dual latitude = dual_atan2(p[2], h);
  put(longitude, &values[0], partials);
  put(latitude, &values[1], partials != NULL ? partials + 6 : NULL);
  return MYB_OK;
}

myb_status myb_measure_light_radec(void* light, double t, const double y[6], double* values,
                                   double* partials) {
  const myb_light_measure* m = light;
  MYB_CHECK(values != NULL);
  light_path path;
  const myb_status status = solve(m, t, y, &path);
  if (status != MYB_OK) {
    return status;
  }
  myb_dual p[3];
  direction(m, &path, p);
  return angles(p, values, partials);
}

myb_status myb_measure_light_azel(void* light, double t, const double y[6], double* values,
                                  double* partials) {
  const myb_light_measure* m = light;
  MYB_CHECK(values != NULL);
  light_path path;
  const myb_status status = solve(m, t, y, &path);
  if (status != MYB_OK) {
    return status;
  }
  myb_dual p[3];
  direction(m, &path, p);
  /* into the body-fixed frame at t, then along north, east, up of the station */
  const double* r = path.axes.rotation.m;
  const double* e = m->station->enu;
  myb_dual fixed[3];
  for (int i = 0; i < 3; ++i) {
    fixed[i] = myb_dual_add(myb_dual_add(myb_dual_scale(p[0], r[3 * i]),
                                         myb_dual_scale(p[1], r[3 * i + 1])),
                            myb_dual_scale(p[2], r[3 * i + 2]));
  }
  myb_dual local[3]; /* east, north, up */
  for (int i = 0; i < 3; ++i) {
    local[i] = myb_dual_add(myb_dual_add(myb_dual_scale(fixed[0], e[3 * i]),
                                         myb_dual_scale(fixed[1], e[3 * i + 1])),
                            myb_dual_scale(fixed[2], e[3 * i + 2]));
  }
  const myb_dual north_east_up[3] = {local[1], local[0], local[2]};
  return angles(north_east_up, values, partials);
}
