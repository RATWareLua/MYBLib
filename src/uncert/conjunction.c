/* SPDX-License-Identifier: BSD-3-Clause */
/* Closest approaches of two objects and the probability of their collision
 * (myb/conjunction.h). */
#include <math.h> /* sqrt, fabs */
#include <stdint.h>
#include <string.h>

#include "core/check.h"
#include "myb/cmath.h"
#include "myb/conjunction.h"
#include "myb/time.h"
#include "num/erf.h"
#include "num/vec3.h"

enum {
  BISECTIONS = 200, /* more than the bits of a double */
  DOUBLINGS = 11,   /* of the panels of the quadrature: 4 .. 4096 */
  CHAN_TERMS = 10000,
  GOLDEN_STEPS = 120
};

static const double pi = 3.14159265358979323846264338327950288;

/* ---- Closest approaches ----------------------------------------------------------------------- */

typedef struct pair_of_providers {
  myb_state_function first;
  void* first_context;
  myb_state_function second;
  void* second_context;
} pair_of_providers;

/* The state of the second from the first at t. */
static myb_status relative_at(const pair_of_providers* p, double t, double rel[6]) {
  double a[6];
  double b[6];
  myb_status st = p->first(p->first_context, myb_time_from_seconds(t), a);
  if (st == MYB_OK) {
    st = p->second(p->second_context, myb_time_from_seconds(t), b);
  }
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 6; ++i) {
    rel[i] = b[i] - a[i];
  }
  return MYB_OK;
}

static void fill(myb_approach* out, double t, const double rel[6]) {
  out->t = t;
  out->distance = myb_vec3_norm(rel);
  memcpy(out->position, rel, 3 * sizeof *rel);
  memcpy(out->velocity, rel + 3, 3 * sizeof *rel);
}

myb_status myb_approach_find(myb_state_function first, void* first_context,
                             myb_state_function second, void* second_context, double t0,
                             double t1, double step, myb_approach* approaches, size_t capacity,
                             size_t* count) {
  MYB_CHECK(first != NULL && second != NULL && count != NULL &&
            (approaches != NULL || capacity == 0));
  *count = 0;
  if (!(t1 > t0) || !(step > 0.0) || !isfinite(t0) || !isfinite(t1) || !isfinite(step) ||
      !((t1 - t0) / step < 1e12)) {
    return MYB_E_DOMAIN;
  }
  const pair_of_providers p = {first, first_context, second, second_context};
  double rel[6];
  myb_status st = relative_at(&p, t0, rel);
  if (st != MYB_OK) {
    return st;
  }
  double before = t0;
  double rate_before = myb_vec3_dot(rel, rel + 3);
  const int64_t steps = (int64_t)floor((t1 - t0) / step) + 1;
  /* bound: the steps of the interval, at most 1e12 */
  for (int64_t k = 1; k <= steps && before < t1; ++k) {
    const double after = k == steps || t0 + step * (double)k > t1 ? t1 : t0 + step * (double)k;
    st = relative_at(&p, after, rel);
    if (st != MYB_OK) {
      return st;
    }
    const double rate_after = myb_vec3_dot(rel, rel + 3);
    if (rate_before < 0.0 && rate_after >= 0.0) {
      double lo = before;
      double hi = after;
      /* bound: BISECTIONS halvings */
      for (int i = 0; i < BISECTIONS; ++i) {
        const double mid = lo + 0.5 * (hi - lo);
        if (mid == lo || mid == hi) {
          break;
        }
        double at[6];
        st = relative_at(&p, mid, at);
        if (st != MYB_OK) {
          return st;
        }
        if (myb_vec3_dot(at, at + 3) < 0.0) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      if (*count < capacity) {
        double at[6];
        st = relative_at(&p, hi, at);
        if (st != MYB_OK) {
          return st;
        }
        fill(&approaches[*count], hi, at);
      }
      ++*count;
    }
    before = after;
    rate_before = rate_after;
  }
  return MYB_OK;
}

/* The cubic through the relative states at the two ends of a step h, at the fraction s of
 * it: position and velocity. */
static void hermite(const double r0[6], const double r1[6], double h, double s, double out[6]) {
  const double s2 = s * s;
  const double s3 = s2 * s;
  const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
  const double h10 = s3 - 2.0 * s2 + s;
  const double h01 = -2.0 * s3 + 3.0 * s2;
  const double h11 = s3 - s2;
  const double d00 = 6.0 * s2 - 6.0 * s;
  const double d10 = 3.0 * s2 - 4.0 * s + 1.0;
  const double d01 = -d00;
  const double d11 = 3.0 * s2 - 2.0 * s;
  for (int i = 0; i < 3; ++i) {
    out[i] = h00 * r0[i] + h10 * h * r0[3 + i] + h01 * r1[i] + h11 * h * r1[3 + i];
    out[3 + i] = (d00 * r0[i] + d01 * r1[i]) / h + d10 * r0[3 + i] + d11 * r1[3 + i];
  }
}

/* The two orbits of a search, with the covariances of their states when there are any. */
typedef struct pair_of_orbits {
  myb_orbit* first;
  myb_orbit* second;
  const double* first_covariance; /* NULL: the approaches alone */
  const double* second_covariance;
  double radius;
} pair_of_orbits;

/* The two transition matrices at one end of a step. */
typedef struct transitions {
  double first[36];
  double second[36];
} transitions;

/* The state of the second from the first at t, and the transition matrices when wanted. */
static myb_status pair_at(const pair_of_orbits* p, double t, double rel[6], transitions* phi) {
  double a[6];
  double b[6];
  myb_status st = phi != NULL ? myb_orbit_partials(p->first, t, a, phi->first, NULL)
                              : myb_orbit_state(p->first, t, a);
  if (st == MYB_OK) {
    st = phi != NULL ? myb_orbit_partials(p->second, t, b, phi->second, NULL)
                     : myb_orbit_state(p->second, t, b);
  }
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 6; ++i) {
    rel[i] = b[i] - a[i];
  }
  return MYB_OK;
}

/* out += R P R^T for the rows of the positions R of the transition matrix at the fraction s
 * of a step h, on the cubic through its two ends: the rows of the velocities are the rates
 * of those of the positions. */
static void add_position_covariance(const double phi0[36], const double phi1[36], double h,
                                    double s, const double p[36], double out[9]) {
  const double s2 = s * s;
  const double s3 = s2 * s;
  const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
  const double h10 = s3 - 2.0 * s2 + s;
  const double h01 = -2.0 * s3 + 3.0 * s2;
  const double h11 = s3 - s2;
  double rows[18];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      rows[6 * i + j] = h00 * phi0[6 * i + j] + h10 * h * phi0[6 * (3 + i) + j] +
                        h01 * phi1[6 * i + j] + h11 * h * phi1[6 * (3 + i) + j];
    }
  }
  double rp[18];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      double sum = 0.0;
      for (int k = 0; k < 6; ++k) {
        sum += rows[6 * i + k] * p[6 * k + j];
      }
      rp[6 * i + j] = sum;
    }
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double sum = 0.0;
      for (int k = 0; k < 6; ++k) {
        sum += rp[6 * i + k] * rows[6 * j + k];
      }
      out[3 * i + j] += sum;
    }
  }
}

/* The approaches of the two orbits up to t1; with covariances every one is made a
 * conjunction. */
static myb_status scan_orbits(const pair_of_orbits* p, double t1, double step,
                              myb_approach* approaches, myb_conjunction* conjunctions,
                              size_t capacity, size_t* count) {
  *count = 0;
  const int carried = p->first_covariance != NULL;
  const double t0 = myb_orbit_time(p->first);
  if (myb_orbit_time(p->second) != t0 || !(t1 > t0) || !(step > 0.0) || !isfinite(t1) ||
      !isfinite(step) || !((t1 - t0) / step < 1e12)) {
    return MYB_E_DOMAIN;
  }
  double before[6];
  transitions phi_before;
  transitions phi_after;
  myb_status st = pair_at(p, t0, before, carried ? &phi_before : NULL);
  if (st != MYB_OK) {
    return st;
  }
  double t_before = t0;
  const int64_t steps = (int64_t)floor((t1 - t0) / step) + 1;
  /* bound: the steps of the interval, at most 1e12 */
  for (int64_t k = 1; k <= steps && t_before < t1; ++k) {
    const double t_after = k == steps || t0 + step * (double)k > t1 ? t1 : t0 + step * (double)k;
    double after[6];
    st = pair_at(p, t_after, after, carried ? &phi_after : NULL);
    if (st != MYB_OK) {
      return st;
    }
    if (myb_vec3_dot(before, before + 3) < 0.0 && myb_vec3_dot(after, after + 3) >= 0.0) {
      const double h = t_after - t_before;
      double lo = 0.0;
      double hi = 1.0;
      /* bound: BISECTIONS halvings */
      for (int i = 0; i < BISECTIONS; ++i) {
        const double mid = lo + 0.5 * (hi - lo);
        if (mid == lo || mid == hi) {
          break;
        }
        double at[6];
        hermite(before, after, h, mid, at);
        if (myb_vec3_dot(at, at + 3) < 0.0) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      if (*count < capacity) {
        double at[6];
        hermite(before, after, h, hi, at);
        if (carried) {
          myb_conjunction* c = &conjunctions[*count];
          memset(c, 0, sizeof *c);
          fill(&c->approach, t_before + hi * h, at);
          add_position_covariance(phi_before.first, phi_after.first, h, hi, p->first_covariance,
                                  c->covariance);
          add_position_covariance(phi_before.second, phi_after.second, h, hi,
                                  p->second_covariance, c->covariance);
          st = myb_encounter_plane(c->approach.position, c->approach.velocity, c->covariance,
                                   p->radius, &c->encounter);
          if (st == MYB_OK) {
            st = myb_collision_probability(&c->encounter, &c->probability);
          }
          if (st != MYB_OK) {
            return st;
          }
        } else {
          fill(&approaches[*count], t_before + hi * h, at);
        }
      }
      ++*count;
    }
    memcpy(before, after, sizeof before);
    if (carried) {
      phi_before = phi_after;
    }
    t_before = t_after;
  }
  return MYB_OK;
}

myb_status myb_approach_orbits(myb_orbit* first, myb_orbit* second, double t1, double step,
                               myb_approach* approaches, size_t capacity, size_t* count) {
  MYB_CHECK(first != NULL && second != NULL && count != NULL &&
            (approaches != NULL || capacity == 0));
  const pair_of_orbits p = {.first = first, .second = second};
  return scan_orbits(&p, t1, step, approaches, NULL, capacity, count);
}

myb_status myb_conjunction_find(myb_orbit* first, const double first_covariance[36],
                                myb_orbit* second, const double second_covariance[36],
                                double radius, double t1, double step,
                                myb_conjunction* conjunctions, size_t capacity, size_t* count) {
  MYB_CHECK(first != NULL && second != NULL && first_covariance != NULL &&
            second_covariance != NULL && count != NULL &&
            (conjunctions != NULL || capacity == 0));
  const pair_of_orbits p = {first, second, first_covariance, second_covariance, radius};
  return scan_orbits(&p, t1, step, NULL, conjunctions, capacity, count);
}

/* ---- The probability of collision ------------------------------------------------------------- */

static int valid(const myb_encounter* e) {
  return e->sigma[0] > 0.0 && e->sigma[1] > 0.0 && e->radius > 0.0 &&
         isfinite(e->miss[0] + e->miss[1] + e->sigma[0] + e->sigma[1] + e->radius);
}

myb_status myb_encounter_plane(const double position[3], const double velocity[3],
                               const double covariance[9], double radius,
                               myb_encounter* encounter) {
  MYB_CHECK(position != NULL && velocity != NULL && covariance != NULL && encounter != NULL);
  const double speed = myb_vec3_norm(velocity);
  double scale = 0.0;
  for (int i = 0; i < 9; ++i) {
    scale = fmax(scale, fabs(covariance[i]));
  }
  if (!(speed > 0.0) || !isfinite(speed) || !myb_vec3_finite(position) || !(radius > 0.0) ||
      !isfinite(radius) || !isfinite(scale) || !(scale > 0.0)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < i; ++j) {
      if (!(fabs(covariance[3 * i + j] - covariance[3 * j + i]) <= 1e-12 * scale)) {
        return MYB_E_DOMAIN;
      }
    }
  }
  /* two unit vectors across the velocity */
  double w[3];
  for (int i = 0; i < 3; ++i) {
    w[i] = velocity[i] / speed;
  }
  int axis = 0;
  for (int i = 1; i < 3; ++i) {
    if (fabs(w[i]) < fabs(w[axis])) {
      axis = i;
    }
  }
  double e[3] = {0.0, 0.0, 0.0};
  e[axis] = 1.0;
  double x[3];
  double y[3];
  myb_vec3_cross(w, e, x);
  const double xn = myb_vec3_norm(x);
  for (int i = 0; i < 3; ++i) {
    x[i] /= xn;
  }
  myb_vec3_cross(w, x, y);
  /* the covariance and the miss in that plane */
  double cx[3];
  double cy[3];
  for (int i = 0; i < 3; ++i) {
    cx[i] = covariance[3 * i] * x[0] + covariance[3 * i + 1] * x[1] + covariance[3 * i + 2] * x[2];
    cy[i] = covariance[3 * i] * y[0] + covariance[3 * i + 1] * y[1] + covariance[3 * i + 2] * y[2];
  }
  const double pxx = myb_vec3_dot(x, cx);
  const double pxy = myb_vec3_dot(x, cy);
  const double pyy = myb_vec3_dot(y, cy);
  const double mx = myb_vec3_dot(position, x);
  const double my = myb_vec3_dot(position, y);
  /* its principal axes: eigenvalues of [[pxx, pxy], [pxy, pyy]] */
  const double half = 0.5 * (pxx - pyy);
  const double root = sqrt(half * half + pxy * pxy);
  const double mean = 0.5 * (pxx + pyy);
  const double large = mean + root;
  const double small = mean - root;
  if (!(small > 0.0) || !isfinite(large)) {
    return MYB_E_DOMAIN;
  }
  /* the axis of the larger: (pxy, large - pxx) or (large - pyy, pxy), whichever is longer */
  double ax = large - pyy;
  double ay = pxy;
  if (fabs(large - pxx) > fabs(ax)) {
    ax = pxy;
    ay = large - pxx;
  }
  const double an = sqrt(ax * ax + ay * ay);
  if (an > 0.0) {
    ax /= an;
    ay /= an;
  } else { /* a circle: any axes */
    ax = 1.0;
    ay = 0.0;
  }
  encounter->miss[0] = ax * mx + ay * my;
  encounter->miss[1] = -ay * mx + ax * my;
  encounter->sigma[0] = sqrt(large);
  encounter->sigma[1] = sqrt(small);
  encounter->radius = radius;
  return MYB_OK;
}

/* The Gaussian of the encounter over the chord of the disk at x = R sin(theta), times the
 * width of the strip: the integrand in theta. */
static double strip(const myb_encounter* e, double theta) {
  const double radius = e->radius;
  const double x = radius * myb_sin(theta);
  const double w = radius * myb_cos(theta);
  const double dx = (x - e->miss[0]) / e->sigma[0];
  const double along = myb_exp(-0.5 * dx * dx);
  const double my = fabs(e->miss[1]);
  const double s = 1.4142135623730951 * e->sigma[1];
  double across = 0.0;
  if (my > w) { /* the chord lies to one side of the centre of the Gaussian: tails */
    across = 0.5 * (myb_internal_erfc((my - w) / s) - myb_internal_erfc((my + w) / s));
  } else {
    across = 0.5 * (myb_internal_erf((w - my) / s) + myb_internal_erf((w + my) / s));
  }
  return along * across * w;
}

myb_status myb_collision_probability(const myb_encounter* encounter, double* probability) {
  MYB_CHECK(encounter != NULL && probability != NULL);
  if (!valid(encounter)) {
    return MYB_E_DOMAIN;
  }
  /* Gauss-Legendre, 8 nodes a panel */
  static const double node[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267,
                                 0.9602898564975363};
  static const double weight[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745,
                                   0.1012285362903763};
  const double factor = 1.0 / (sqrt(2.0 * pi) * encounter->sigma[0]);
  double last = -1.0;
  int panels = 4;
  for (int doubling = 0; doubling < DOUBLINGS; ++doubling) {
    const double width = pi / panels;
    double sum = 0.0;
    for (int p = 0; p < panels; ++p) {
      const double centre = -0.5 * pi + (p + 0.5) * width;
      for (int k = 0; k < 4; ++k) {
        const double offset = 0.5 * width * node[k];
        sum += weight[k] * (strip(encounter, centre - offset) + strip(encounter, centre + offset));
      }
    }
    const double value = factor * 0.5 * width * sum;
    if (!isfinite(value)) {
      return MYB_E_DOMAIN;
    }
    if (fabs(value - last) <= 1e-13 * value || (value < 1e-300 && last < 1e-300 && last >= 0.0)) {
      *probability = value > 1.0 ? 1.0 : value;
      return MYB_OK;
    }
    last = value;
    panels *= 2;
  }
  return MYB_E_NOCONV;
}

double myb_collision_probability_chan(const myb_encounter* encounter) {
  MYB_CHECK(encounter != NULL);
  if (!valid(encounter)) {
    return (double)NAN;
  }
  const double sx = encounter->sigma[0];
  const double sy = encounter->sigma[1];
  const double u = 0.5 * encounter->radius * encounter->radius / (sx * sy);
  const double v = 0.5 * (encounter->miss[0] * encounter->miss[0] / (sx * sx) +
                          encounter->miss[1] * encounter->miss[1] / (sy * sy));
  if (!(u < 700.0) || !(v < 700.0)) {
    return (double)NAN; /* e^-u or e^-v is beyond the doubles */
  }
  /* e^-v sum_m v^m / m! (1 - e^-u sum_{k <= m} u^k / k!) with the sums exchanged, so that
   * every term is positive: sum_{k >= 1} p_u(k) F_v(k - 1), p the terms and F the partial
   * sums of the two Poisson series */
  double p_u = myb_exp(-u); /* p_u(0) */
  double p_v = myb_exp(-v); /* p_v(0) */
  double f_v = 0.0;         /* F_v(k - 1) */
  double sum = 0.0;
  for (int k = 1; k < CHAN_TERMS; ++k) {
    f_v += p_v;         /* now F_v(k - 1) */
    p_v *= v / k;       /* p_v(k) */
    p_u *= u / k;       /* p_u(k) */
    const double part = p_u * f_v;
    sum += part;
    if (k > u && part <= 1e-17 * sum) {
      break;
    }
  }
  return sum > 1.0 ? 1.0 : sum;
}

myb_status myb_collision_probability_max(const myb_encounter* encounter, double* probability,
                                         double* factor) {
  MYB_CHECK(encounter != NULL && probability != NULL && factor != NULL);
  if (!valid(encounter)) {
    return MYB_E_DOMAIN;
  }
  const double miss = sqrt(encounter->miss[0] * encounter->miss[0] +
                           encounter->miss[1] * encounter->miss[1]);
  if (miss <= encounter->radius) { /* the centres pass inside the disk: certain as it shrinks */
    *probability = 1.0;
    *factor = 0.0;
    return MYB_OK;
  }
  /* the factor that puts the miss at sqrt(2) sigma is near the maximum; searched over three
   * decades each way in the logarithm of the factor by the golden section */
  const double mx = encounter->miss[0] / encounter->sigma[0];
  const double my = encounter->miss[1] / encounter->sigma[1];
  const double guess = sqrt(0.5 * (mx * mx + my * my));
  const double ratio = 0.6180339887498949;
  double lo = myb_log(guess) - myb_log(1e3);
  double hi = myb_log(guess) + myb_log(1e3);
  double x1 = hi - ratio * (hi - lo);
  double x2 = lo + ratio * (hi - lo);
  double f1 = 0.0;
  double f2 = 0.0;
  myb_encounter scaled = *encounter;
  myb_status st = MYB_OK;
  for (int pass = 0; pass < 2 && st == MYB_OK; ++pass) {
    const double f = myb_exp(pass == 0 ? x1 : x2);
    scaled.sigma[0] = f * encounter->sigma[0];
    scaled.sigma[1] = f * encounter->sigma[1];
    st = myb_collision_probability(&scaled, pass == 0 ? &f1 : &f2);
  }
  /* bound: GOLDEN_STEPS */
  for (int k = 0; k < GOLDEN_STEPS && st == MYB_OK && hi - lo > 1e-12; ++k) {
    double x = 0.0;
    double* f = NULL;
    if (f1 < f2) {
      lo = x1;
      x1 = x2;
      f1 = f2;
      x2 = lo + ratio * (hi - lo);
      x = x2;
      f = &f2;
    } else {
      hi = x2;
      x2 = x1;
      f2 = f1;
      x1 = hi - ratio * (hi - lo);
      x = x1;
      f = &f1;
    }
    const double scale = myb_exp(x);
    scaled.sigma[0] = scale * encounter->sigma[0];
    scaled.sigma[1] = scale * encounter->sigma[1];
    st = myb_collision_probability(&scaled, f);
  }
  if (st != MYB_OK) {
    return st;
  }
  *probability = f1 > f2 ? f1 : f2;
  *factor = myb_exp(f1 > f2 ? x1 : x2);
  return MYB_OK;
}
