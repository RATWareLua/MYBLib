/* SPDX-License-Identifier: BSD-3-Clause */
/* The Q-law (myb/thrust.h): a feedback law for low thrust from a measure of the distance to
 * the target orbit. */
#include <math.h> /* sqrt, fabs, floor */
#include <stdint.h>
#include <string.h>

#include "core/check.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "myb/thrust.h"
#include "num/vec3.h"

enum {
  ELEMENTS = 5, /* a, e, i, raan, argp */
  ANOMALIES = 72 /* places on the orbit where the best rate of Q is looked for */
};

static const double two_pi = 6.283185307179586476925286766559;

static const double penalty_steepness = 100.0;

/* The shortest way round between two angles. */
static double angle_apart(double a, double b) {
  const double c = myb_cos(a - b);
  return myb_acos(c > 1.0 ? 1.0 : (c < -1.0 ? -1.0 : c));
}

/* Q for the elements oe = (a, e, i, raan, argp) and a unit thrust acceleration: the sum over
 * the targeted elements of the squared distance over the largest rate at which the element
 * can be changed anywhere on the orbit (Petropoulos 2004, eqs. 4-13). */
static double q_of(const myb_qlaw* law, const double oe[ELEMENTS]) {
  const double mu = law->mu;
  const double a = oe[0];
  const double e = oe[1] > 1e-8 ? oe[1] : 1e-8;
  const double inc = oe[2];
  const double argp = oe[4];
  const double p = a * (1.0 - e * e);
  const double h = sqrt(mu * p);
  const double sin_i = fabs(myb_sin(inc)) > 1e-8 ? fabs(myb_sin(inc)) : 1e-8;
  const double sw = myb_sin(argp);
  const double cw = myb_cos(argp);
  double q = 0.0;
  if (law->weight[0] > 0.0) {
    const double rate = 2.0 * a * sqrt(a * (1.0 + e) / (mu * (1.0 - e)));
    const double d = a - law->target[0];
    /* the scaling that keeps the semi-major axis from running away (m = 3, n = 4, r = 2) */
    const double x = d / (3.0 * law->target[0]);
    const double scaling = sqrt(1.0 + x * x * x * x);
    q += law->weight[0] * scaling * (d / rate) * (d / rate);
  }
  if (law->weight[1] > 0.0) {
    const double rate = 2.0 * p / h;
    const double d = oe[1] - law->target[1];
    q += law->weight[1] * (d / rate) * (d / rate);
  }
  if (law->weight[2] > 0.0) {
    const double rate = p / (h * (sqrt(1.0 - e * e * sw * sw) - e * fabs(cw)));
    const double d = inc - law->target[2];
    q += law->weight[2] * (d / rate) * (d / rate);
  }
  const double node_rate = p / (h * sin_i * (sqrt(1.0 - e * e * cw * cw) - e * fabs(sw)));
  if (law->weight[3] > 0.0) {
    const double d = angle_apart(oe[3], law->target[3]);
    q += law->weight[3] * (d / node_rate) * (d / node_rate);
  }
  if (law->weight[4] > 0.0) {
    /* in the plane: the true anomaly at which the rate is largest; out of it: through the
     * node; the two blended as Petropoulos does (b = 0.01) */
    const double alpha = (1.0 - e * e) / (2.0 * e * e * e);
    const double beta = sqrt(alpha * alpha + 1.0 / 27.0);
    double cos_nu = myb_cbrt(alpha + beta) - myb_cbrt(beta - alpha) - 1.0 / e;
    cos_nu = cos_nu > 1.0 ? 1.0 : (cos_nu < -1.0 ? -1.0 : cos_nu);
    const double r = p / (1.0 + e * cos_nu);
    const double in_plane =
        sqrt(p * p * cos_nu * cos_nu + (p + r) * (p + r) * (1.0 - cos_nu * cos_nu)) / (e * h);
    const double out_of_plane = node_rate * fabs(myb_cos(inc));
    const double rate = (in_plane + 0.01 * out_of_plane) / 1.01;
    const double d = angle_apart(argp, law->target[4]);
    q += law->weight[4] * (d / rate) * (d / rate);
  }
  if (law->rp_min > 0.0) {
    const double weight = law->rp_weight > 0.0 ? law->rp_weight : 1.0;
    q *= 1.0 + weight * myb_exp(penalty_steepness * (1.0 - a * (1.0 - e) / law->rp_min));
  }
  return q;
}

static int law_valid(const myb_qlaw* law) {
  double sum = law->mu + law->rp_min + law->rp_weight;
  int any = 0;
  for (int k = 0; k < ELEMENTS; ++k) {
    sum += law->target[k] + law->weight[k];
    if (!(law->weight[k] >= 0.0)) {
      return 0;
    }
    any = any || law->weight[k] > 0.0;
  }
  return any && isfinite(sum) && law->mu > 0.0 && law->target[0] > 0.0 && law->rp_min >= 0.0 &&
         law->rp_weight >= 0.0;
}

/* The rates of the elements per unit of radial, transverse and normal acceleration at the
 * true anomaly nu (the Gauss equations), carried through the gradient of Q: d[j] = sum_k
 * gradient[k] d(oe_k)/dt per unit acceleration along j. */
static void rate_of_q(double mu, const myb_elements* orbit, double nu,
                      const double gradient[ELEMENTS], double d[3]) {
  const double e = orbit->e > 1e-8 ? orbit->e : 1e-8;
  const double p = orbit->p;
  const double a = p / (1.0 - e * e);
  const double h = sqrt(mu * p);
  const double sn = myb_sin(nu);
  const double cn = myb_cos(nu);
  const double r = p / (1.0 + orbit->e * cn);
  const double su = myb_sin(nu + orbit->argp);
  const double cu = myb_cos(nu + orbit->argp);
  const double sin_i = myb_sin(orbit->i);
  const double safe_sin = fabs(sin_i) > 1e-8 ? sin_i : (sin_i < 0.0 ? -1e-8 : 1e-8);
  const double rates[ELEMENTS][3] = {
      {2.0 * a * a / h * orbit->e * sn, 2.0 * a * a / h * p / r, 0.0},
      {p * sn / h, ((p + r) * cn + r * orbit->e) / h, 0.0},
      {0.0, 0.0, r * cu / h},
      {0.0, 0.0, r * su / (h * safe_sin)},
      {-p * cn / (e * h), (p + r) * sn / (e * h), -r * su * myb_cos(orbit->i) / (h * safe_sin)}};
  d[0] = d[1] = d[2] = 0.0;
  for (int k = 0; k < ELEMENTS; ++k) {
    for (int j = 0; j < 3; ++j) {
      d[j] += gradient[k] * rates[k][j];
    }
  }
}

myb_status myb_qlaw_direction(const myb_qlaw* law, const myb_elements* now, double direction[3],
                              double* q, double* q_rate, double* effectivity) {
  MYB_CHECK(law != NULL && now != NULL && direction != NULL);
  if (!law_valid(law) || !(now->e < 1.0) || !(now->e >= 0.0) || !(now->p > 0.0) ||
      !isfinite(now->p + now->e + now->i + now->raan + now->argp + now->nu)) {
    return MYB_E_DOMAIN;
  }
  const double mu = law->mu;
  const double e = now->e > 1e-8 ? now->e : 1e-8;
  const double a = now->p / (1.0 - e * e);
  const double oe[ELEMENTS] = {a, now->e, now->i, now->raan, now->argp};
  /* the gradient of Q by central differences in the elements */
  double gradient[ELEMENTS];
  const double steps[ELEMENTS] = {1e-6 * a, 1e-7, 1e-7, 1e-7, 1e-7};
  for (int k = 0; k < ELEMENTS; ++k) {
    double up[ELEMENTS];
    double down[ELEMENTS];
    memcpy(up, oe, sizeof up);
    memcpy(down, oe, sizeof down);
    up[k] += steps[k];
    down[k] -= steps[k];
    if (k == 1 && down[k] < 0.0) { /* at a circle: one-sided */
      down[k] = oe[k];
    }
    gradient[k] = (q_of(law, up) - q_of(law, down)) / (up[k] - down[k]);
  }
  double d[3];
  rate_of_q(mu, now, now->nu, gradient, d);
  const double size = myb_vec3_norm(d);
  if (!isfinite(size)) {
    return MYB_E_DOMAIN;
  }
  for (int j = 0; j < 3; ++j) {
    direction[j] = size > 0.0 ? -d[j] / size : 0.0;
  }
  if (q != NULL) {
    *q = q_of(law, oe);
  }
  if (q_rate != NULL) {
    *q_rate = -size;
  }
  if (effectivity != NULL) {
    /* against the best place of the orbit for it: the gradient does not depend on the place */
    double best = size;
    for (int k = 0; k < ANOMALIES; ++k) {
      double there[3];
      rate_of_q(mu, now, two_pi * (double)k / (double)ANOMALIES, gradient, there);
      const double at = myb_vec3_norm(there);
      best = at > best ? at : best;
    }
    *effectivity = best > 0.0 ? size / best : 0.0;
  }
  return MYB_OK;
}

int myb_qlaw_arrived(const myb_qlaw* law, const myb_elements* now, const double tolerance[5]) {
  MYB_CHECK(law != NULL && now != NULL && tolerance != NULL);
  const double a = now->p / (1.0 - now->e * now->e);
  const double apart[ELEMENTS] = {fabs(a - law->target[0]), fabs(now->e - law->target[1]),
                                  fabs(now->i - law->target[2]),
                                  angle_apart(now->raan, law->target[3]),
                                  angle_apart(now->argp, law->target[4])};
  for (int k = 0; k < ELEMENTS; ++k) {
    if (law->weight[k] > 0.0 && !(apart[k] <= tolerance[k])) {
      return 0;
    }
  }
  return 1;
}

static int engine_valid(const myb_qlaw_term* term) {
  return term->law != NULL && term->thrust >= 0.0 && term->exhaust > 0.0 &&
         term->effectivity >= 0.0 && term->effectivity < 1.0 &&
         isfinite(term->thrust + term->exhaust);
}

/* The acceleration of the engine at (x, v) for the mass given; *burning (may be NULL) says
 * whether it is on. */
static myb_status engine_push(const myb_qlaw_term* term, double mass, const double x[3],
                              const double v[3], double a[3], int* burning) {
  a[0] = a[1] = a[2] = 0.0;
  if (burning != NULL) {
    *burning = 0;
  }
  myb_pv pv;
  memcpy(pv.r, x, sizeof pv.r);
  memcpy(pv.v, v, sizeof pv.v);
  myb_elements now;
  myb_status st = myb_elements_from_pv(term->law->mu, &pv, &now);
  if (st != MYB_OK) {
    return st;
  }
  if (!(now.e < 1.0)) {
    return MYB_E_DOMAIN;
  }
  if (myb_qlaw_arrived(term->law, &now, term->tolerance)) {
    return MYB_OK; /* the engine is off */
  }
  if (!(mass > 0.0)) {
    return MYB_E_DOMAIN; /* nothing left to push */
  }
  double direction[3];
  double effectivity = 1.0;
  st = myb_qlaw_direction(term->law, &now, direction, NULL, NULL,
                          term->effectivity > 0.0 ? &effectivity : NULL);
  if (st != MYB_OK) {
    return st;
  }
  if (effectivity < term->effectivity) {
    return MYB_OK; /* a place where thrust gains little: the engine waits */
  }
  /* radial, transverse, normal -> the axes of x */
  double radial[3];
  double normal[3];
  double transverse[3];
  const double rn = myb_vec3_norm(x);
  myb_vec3_cross(x, v, normal);
  const double hn = myb_vec3_norm(normal);
  for (int i = 0; i < 3; ++i) {
    radial[i] = x[i] / rn;
    normal[i] /= hn;
  }
  myb_vec3_cross(normal, radial, transverse);
  const double f = term->thrust / mass;
  for (int i = 0; i < 3; ++i) {
    a[i] = f * (direction[0] * radial[i] + direction[1] * transverse[i] +
                direction[2] * normal[i]);
  }
  if (burning != NULL) {
    *burning = 1;
  }
  return MYB_OK;
}

myb_status myb_term_qlaw(void* context, double t, const double x[3], const double v[3],
                         double a[3]) {
  MYB_CHECK(context != NULL && x != NULL && v != NULL && a != NULL);
  const myb_qlaw_term* term = context;
  a[0] = a[1] = a[2] = 0.0;
  if (!engine_valid(term) || !(term->mass0 > 0.0) || !isfinite(term->mass0 + term->t0)) {
    return MYB_E_DOMAIN;
  }
  const double burned = t > term->t0 ? term->thrust / term->exhaust * (t - term->t0) : 0.0;
  return engine_push(term, term->mass0 - burned, x, v, a, NULL);
}

/* ---- The mass as a part of the state --------------------------------------------------------- */

myb_status myb_qlaw_flight_rhs(void* context, double t, const double* y, double* dydt) {
  MYB_CHECK(context != NULL && y != NULL && dydt != NULL);
  const myb_qlaw_flight* flight = context;
  if (flight->engine == NULL || !engine_valid(flight->engine)) {
    return MYB_E_DOMAIN;
  }
  double field[3];
  if (flight->field != NULL) {
    const myb_status st = myb_force_model_acceleration(flight->field, t, y, y + 3, field);
    if (st != MYB_OK) {
      return st;
    }
  } else {
    const double r = myb_vec3_norm(y);
    if (!(r > 0.0)) {
      return MYB_E_DOMAIN;
    }
    const double k = -flight->engine->law->mu / (r * r * r);
    for (int i = 0; i < 3; ++i) {
      field[i] = k * y[i];
    }
  }
  double push[3];
  int burning = 0;
  const myb_status st = engine_push(flight->engine, y[6], y, y + 3, push, &burning);
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 3; ++i) {
    dydt[i] = y[3 + i];
    dydt[3 + i] = field[i] + push[i];
  }
  dydt[6] = burning ? -flight->engine->thrust / flight->engine->exhaust : 0.0;
  return MYB_OK;
}

myb_status myb_qlaw_fly(const myb_qlaw_flight* flight, double t0, const double y0[7],
                        double t_limit, double step, const myb_allocator* allocator,
                        myb_qlaw_arrival* arrival) {
  MYB_CHECK(flight != NULL && y0 != NULL && arrival != NULL);
  memset(arrival, 0, sizeof *arrival);
  if (flight->engine == NULL || !engine_valid(flight->engine) || !(t_limit > t0) ||
      !(step > 0.0) || !isfinite(t0) || !isfinite(t_limit) || !isfinite(step) ||
      !((t_limit - t0) / step < 1e12) || !(y0[6] > 0.0) || !(flight->rtol >= 0.0) ||
      !(flight->atol >= 0.0)) {
    return MYB_E_DOMAIN;
  }
  myb_qlaw_flight own = *flight; /* the integrator keeps a pointer that is not to a constant */
  const myb_dop853_options options = {.rtol = flight->rtol > 0.0 ? flight->rtol : 1e-10,
                                      .atol = flight->atol > 0.0 ? flight->atol : 1e-4};
  myb_dop853* integrator = NULL;
  myb_status st = myb_dop853_create(7, &options, allocator, &integrator);
  if (st != MYB_OK) {
    return st;
  }
  st = myb_dop853_start(integrator, myb_qlaw_flight_rhs, &own, t0, y0);
  double t = t0;
  memcpy(arrival->y, y0, sizeof arrival->y);
  const int64_t looks = (int64_t)floor((t_limit - t0) / step) + 1;
  /* bound: the looks of the interval, at most 1e12 */
  for (int64_t k = 0; k <= looks && st == MYB_OK; ++k) {
    myb_pv pv;
    memcpy(pv.r, arrival->y, sizeof pv.r);
    memcpy(pv.v, arrival->y + 3, sizeof pv.v);
    myb_elements now;
    st = myb_elements_from_pv(flight->engine->law->mu, &pv, &now);
    if (st != MYB_OK) {
      break;
    }
    if (now.e < 1.0 && myb_qlaw_arrived(flight->engine->law, &now, flight->engine->tolerance)) {
      arrival->arrived = 1;
      break;
    }
    if (!(t < t_limit)) {
      break;
    }
    const double next = t0 + step * (double)(k + 1);
    st = myb_dop853_integrate(integrator, next < t_limit ? next : t_limit);
    /* on a failure the integrator stands where its last step ended: the state and the time
     * returned are those */
    t = myb_dop853_time(integrator);
    memcpy(arrival->y, myb_dop853_state(integrator), sizeof arrival->y);
  }
  arrival->t = t;
  myb_dop853_destroy(integrator);
  return st;
}
