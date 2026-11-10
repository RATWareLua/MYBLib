/* SPDX-License-Identifier: BSD-3-Clause */
/* Light time, relativistic delay, round trips, aberration, deflection. */
#include <math.h>
#include <string.h>

#include "core/check.h"
#include "myb/cmath.h"
#include "myb/num.h"
#include "myb/observe.h"
#include "num/vec3.h"

enum { MAX_ITERATIONS = 32 };

double myb_shapiro_delay(double gm, double c, double gamma, const double x1[3],
                         const double x2[3], const double xb[3]) {
  MYB_CHECK(x1 != NULL && x2 != NULL && xb != NULL);
  double d1[3];
  double d2[3];
  double d12[3];
  for (int i = 0; i < 3; ++i) {
    d1[i] = x1[i] - xb[i];
    d2[i] = x2[i] - xb[i];
    d12[i] = x2[i] - x1[i];
  }
  const double r1 = myb_vec3_norm(d1);
  const double r2 = myb_vec3_norm(d2);
  const double r12 = myb_vec3_norm(d12);
  const double k = (1.0 + gamma) * gm / (c * c); /* m */
  /* r1 + r2 - r12 loses digits for a grazing ray: (r1 + r2)^2 - r12^2 = 2 (r1 r2 + d1.d2) */
  const double s = r1 + r2;
  const double minus = 2.0 * (r1 * r2 + myb_vec3_dot(d1, d2)) / (s + r12);
  return k / c * myb_log((s + r12 + k) / (minus + k));
}

/* The delays of all bodies for a ray from xt (target end) to xo (observer end at t). */
static myb_status delays(const myb_light_model* m, myb_time t, const double xo[3],
                         const double xt[3], double* total) {
  double sum = 0.0;
  for (size_t b = 0; b < m->count; ++b) {
    const myb_light_body* body = &m->bodies[b];
    double sb[6];
    myb_status st = body->state(body->context, t, sb);
    if (st != MYB_OK) {
      return st;
    }
    /* the body when the ray passes it: t - |x_B - x_o| / c */
    double d[3] = {sb[0] - xo[0], sb[1] - xo[1], sb[2] - xo[2]};
    st = body->state(body->context, myb_time_add(t, -myb_vec3_norm(d) / m->c), sb);
    if (st != MYB_OK) {
      return st;
    }
    sum += myb_shapiro_delay(body->gm, m->c, m->gamma, xt, xo, sb);
  }
  *total = sum;
  return MYB_OK;
}

myb_status myb_light_time(const myb_light_model* model, myb_state_function target, void* context,
                          myb_time t, const double observer[6], myb_light_direction direction,
                          myb_light_path* path) {
  MYB_CHECK(model != NULL && target != NULL && observer != NULL && path != NULL &&
            (model->bodies != NULL || model->count == 0));
  if (!(model->c > 0.0) || !isfinite(model->c) || !isfinite(model->gamma) ||
      (direction != MYB_LIGHT_RECEPTION && direction != MYB_LIGHT_TRANSMISSION)) {
    return MYB_E_DOMAIN;
  }
  const double sign = direction == MYB_LIGHT_RECEPTION ? -1.0 : 1.0;
  const double c = model->c;
  double tau = 0.0;
  double state[6];
  myb_status st = target(context, t, state);
  if (st != MYB_OK) {
    return st;
  }
  double delay = 0.0;
  int it = 0;
  int converged = 0;
  for (; it < MAX_ITERATIONS && !converged; ++it) {
    double d[3] = {state[0] - observer[0], state[1] - observer[1], state[2] - observer[2]};
    const double range = myb_vec3_norm(d);
    if (model->count > 0) {
      st = delays(model, t, observer, state, &delay);
      if (st != MYB_OK) {
        return st;
      }
    }
    const double next = range / c + delay;
    /* converged when the change is a few units of roundoff */
    converged = fabs(next - tau) <= 0x1p-50 * fmax(next, 1e-9);
    tau = next;
    st = target(context, myb_time_add(t, sign * tau), state);
    if (st != MYB_OK) {
      return st;
    }
  }
  if (!converged) {
    return MYB_E_NOCONV;
  }
  double d[3] = {state[0] - observer[0], state[1] - observer[1], state[2] - observer[2]};
  const double range = myb_vec3_norm(d);
  double n[3] = {d[0] / range, d[1] / range, d[2] / range}; /* observer -> target */
  const double nvt = myb_vec3_dot(n, state + 3);
  double dv[3] = {state[3] - observer[3], state[4] - observer[4], state[5] - observer[5]};
  /* d range / dt with dt_target / dt = 1 + sign range_rate / c */
  const double rate = myb_vec3_dot(n, dv) / (1.0 - sign * nvt / c);
  path->t_target = myb_time_add(t, sign * tau);
  memcpy(path->target, state, sizeof state);
  path->light_time = tau;
  path->range = range;
  path->delay = delay;
  path->range_rate = rate;
  path->iterations = it;
  return MYB_OK;
}

typedef struct bounce_target {
  myb_state_function station;
  void* context;
} bounce_target;

myb_status myb_round_trip(const myb_light_model* model, myb_state_function station,
                          void* station_context, myb_state_function target,
                          void* target_context, myb_time t, double* total, myb_time* t_bounce,
                          myb_time* t_transmit) {
  MYB_CHECK(model != NULL && station != NULL && target != NULL && total != NULL);
  double receiver[6];
  myb_status st = station(station_context, t, receiver);
  if (st != MYB_OK) {
    return st;
  }
  /* down-link: the target at the bounce */
  myb_light_path down;
  st = myb_light_time(model, target, target_context, t, receiver, MYB_LIGHT_RECEPTION, &down);
  if (st != MYB_OK) {
    return st;
  }
  /* up-link: the station seen from the target at the bounce (received there) */
  myb_light_path up;
  st = myb_light_time(model, station, station_context, down.t_target, down.target,
                      MYB_LIGHT_RECEPTION, &up);
  if (st != MYB_OK) {
    return st;
  }
  *total = myb_dd_to_double(myb_time_diff(t, up.t_target));
  if (t_bounce != NULL) {
    *t_bounce = down.t_target;
  }
  if (t_transmit != NULL) {
    *t_transmit = up.t_target;
  }
  return MYB_OK;
}

void myb_aberration(const double p[3], const double v[3], double c, double out[3]) {
  MYB_CHECK(p != NULL && v != NULL && out != NULL && c > 0.0);
  const double b[3] = {v[0] / c, v[1] / c, v[2] / c};
  const double bm1 = sqrt(1.0 - myb_vec3_dot(b, b));
  const double pdv = myb_vec3_dot(p, b);
  const double w1 = 1.0 + pdv / (1.0 + bm1);
  double w[3];
  for (int i = 0; i < 3; ++i) {
    w[i] = p[i] * bm1 + w1 * b[i];
  }
  const double r = myb_vec3_norm(w);
  for (int i = 0; i < 3; ++i) {
    out[i] = w[i] / r;
  }
}

void myb_light_deflection(double gm, double c, double gamma, const double p[3], const double q[3],
                          const double e[3], double em, double dlim, double out[3]) {
  MYB_CHECK(p != NULL && q != NULL && e != NULL && out != NULL && c > 0.0 && em > 0.0);
  double qpe[3];
  for (int i = 0; i < 3; ++i) {
    qpe[i] = q[i] + e[i];
  }
  const double qdqpe = myb_vec3_dot(q, qpe);
  const double w = (1.0 + gamma) * gm / (c * c) / em / fmax(qdqpe, dlim);
  double eq[3];
  double peq[3];
  myb_vec3_cross(e, q, eq);
  myb_vec3_cross(p, eq, peq);
  for (int i = 0; i < 3; ++i) {
    out[i] = p[i] + w * peq[i];
  }
}

void myb_direction_radec(const double p[3], double* ra, double* dec) {
  MYB_CHECK(p != NULL && ra != NULL && dec != NULL);
  const double rxy = myb_hypot(p[0], p[1]);
  double a = rxy > 0.0 ? myb_atan2(p[1], p[0]) : 0.0;
  if (a < 0.0) {
    a += 6.283185307179586;
  }
  *ra = a;
  *dec = (rxy > 0.0 || p[2] != 0.0) ? myb_atan2(p[2], rxy) : 0.0;
}
