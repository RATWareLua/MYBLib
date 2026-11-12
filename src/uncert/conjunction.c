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

/* ---- Long encounters ---------------------------------------------------------------------- */

enum {
  SPHERE_LEVELS = 7,   /* doublings of the grid on the sphere */
  PIECES = 4000,       /* of the adaptive integral in time */
  PIECE_STACK = 64     /* halvings of one panel at most */
};

/* Gauss-Legendre, 8 nodes on [-1, 1] (the halves) */
static const double gl8_node[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267,
                                   0.9602898564975363};
static const double gl8_weight[4] = {0.3626837833783620, 0.3137066458778873,
                                     0.2223810344533745, 0.1012285362903763};

/* Gauss-Kronrod 15 (the nodes of Gauss 7 are the odd ones) */
static const double gk_node[8] = {0.991455371120812639206854697526329,
                                  0.949107912342758524526189684047851,
                                  0.864864423359769072789712788640926,
                                  0.741531185599394439863864773280788,
                                  0.586087235467691130294144845693013,
                                  0.405845151377397166906606412076961,
                                  0.207784955007898467600689403773245,
                                  0.0};
static const double gk_weight[8] = {0.022935322010529224963732008058970,
                                    0.063092092629978553290700663189204,
                                    0.104790010322250183839876322541518,
                                    0.140653259715525918745189590510238,
                                    0.169004726639267902826583426598550,
                                    0.190350578064785409913256402421014,
                                    0.204432940075298892414161999234649,
                                    0.209482141084727828012999174891714};
static const double g7_weight[4] = {0.129484966168869693270611432679082,
                                    0.279705391489276667901467771423780,
                                    0.381830050505118944950369775488975,
                                    0.417959183673469387755102040816327};

/* The Cholesky factor of a 3 x 3 matrix (lower, row-major); 0 if not positive definite. */
static int cholesky3(const double a[9], double l[9]) {
  memset(l, 0, 9 * sizeof *l);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = a[3 * i + j];
      for (int k = 0; k < j; ++k) {
        sum -= l[3 * i + k] * l[3 * j + k];
      }
      if (i == j) {
        if (!(sum > 0.0)) {
          return 0;
        }
        l[3 * i + i] = sqrt(sum);
      } else {
        l[3 * i + j] = sum / l[3 * j + j];
      }
    }
  }
  return 1;
}

/* The mean inward flux of a Gaussian of mean m and deviation s: E[(-u)^+]. */
static double inward(double m, double s) {
  if (!(s > 0.0)) {
    return m < 0.0 ? -m : 0.0;
  }
  const double x = m / s;
  /* s phi(x) - m Phi(-x), Phi(-x) = erfc(x / sqrt 2) / 2 */
  return s * 0.3989422804014327 * myb_exp(-0.5 * x * x) -
         m * 0.5 * myb_internal_erfc(x * 0.7071067811865476);
}

myb_status myb_collision_rate(const double mean[6], const double covariance[36], double radius,
                              double* rate) {
  MYB_CHECK(mean != NULL && covariance != NULL && rate != NULL);
  if (!(radius > 0.0) || !isfinite(radius)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 36; ++i) {
    if (!isfinite(covariance[i]) || (i < 6 && !isfinite(mean[i]))) {
      return MYB_E_DOMAIN;
    }
  }
  /* A (position), B (position x velocity), D (velocity) */
  double a[9];
  double b[9];
  double d[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      a[3 * i + j] = 0.5 * (covariance[6 * i + j] + covariance[6 * j + i]);
      b[3 * i + j] = 0.5 * (covariance[6 * i + 3 + j] + covariance[6 * (3 + j) + i]);
      d[3 * i + j] = 0.5 * (covariance[6 * (3 + i) + 3 + j] + covariance[6 * (3 + j) + 3 + i]);
    }
  }
  double l[9];
  if (!cholesky3(a, l)) {
    return MYB_E_DOMAIN;
  }
  /* A^-1 by the factor; K = B^T A^-1 (the velocity given the place), C = D - B^T A^-1 B */
  double inverse[9];
  for (int col = 0; col < 3; ++col) {
    double y[3];
    for (int i = 0; i < 3; ++i) {
      double sum = i == col ? 1.0 : 0.0;
      for (int k = 0; k < i; ++k) {
        sum -= l[3 * i + k] * y[k];
      }
      y[i] = sum / l[3 * i + i];
    }
    for (int i = 2; i >= 0; --i) {
      double sum = y[i];
      for (int k = i + 1; k < 3; ++k) {
        sum -= l[3 * k + i] * inverse[3 * k + col];
      }
      inverse[3 * i + col] = sum / l[3 * i + i];
    }
  }
  double k_gain[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double sum = 0.0;
      for (int m = 0; m < 3; ++m) {
        sum += b[3 * m + i] * inverse[3 * m + j];
      }
      k_gain[3 * i + j] = sum;
    }
  }
  double c[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double sum = d[3 * i + j];
      for (int m = 0; m < 3; ++m) {
        sum -= k_gain[3 * i + m] * b[3 * m + j];
      }
      c[3 * i + j] = sum;
    }
  }
  /* the grid on the sphere has its pole along the mean relative velocity at the centre (given
   * the place there): the flux turns on and off across the great circle at right angles to
   * it — sharply where the velocity is nearly certain — and that circle is then the equator,
   * a boundary of the panels */
  double axis[3];
  for (int i = 0; i < 3; ++i) {
    double ki = 0.0;
    for (int m = 0; m < 3; ++m) {
      ki -= k_gain[3 * i + m] * mean[m];
    }
    axis[i] = mean[3 + i] + ki;
  }
  double frame[9]; /* rows: two directions across the axis, then the axis */
  const double size = myb_vec3_norm(axis);
  if (size > 0.0 && isfinite(size)) {
    for (int i = 0; i < 3; ++i) {
      frame[6 + i] = axis[i] / size;
    }
  } else {
    frame[6] = 0.0;
    frame[7] = 0.0;
    frame[8] = 1.0;
  }
  {
    int least = 0;
    for (int i = 1; i < 3; ++i) {
      if (fabs(frame[6 + i]) < fabs(frame[6 + least])) {
        least = i;
      }
    }
    double e[3] = {0.0, 0.0, 0.0};
    e[least] = 1.0;
    myb_vec3_cross(frame + 6, e, frame);
    const double n0 = myb_vec3_norm(frame);
    for (int i = 0; i < 3; ++i) {
      frame[i] /= n0;
    }
    myb_vec3_cross(frame + 6, frame, frame + 3);
  }
  const double det = l[0] * l[4] * l[8];
  const double norm = 1.0 / (15.749609945722419 * det); /* (2 pi)^(3/2) sqrt(det A) */
  double last = -1.0;
  int panels = 2;
  int around = 16;
  for (int level = 0; level < SPHERE_LEVELS; ++level) { /* bound: SPHERE_LEVELS */
    const double width = 2.0 / (double)panels;
    double sum = 0.0;
    for (int p = 0; p < panels; ++p) {
      const double centre = -1.0 + ((double)p + 0.5) * width;
      for (int q = 0; q < 8; ++q) {
        const double z = centre + (q < 4 ? -1.0 : 1.0) * 0.5 * width * gl8_node[q % 4];
        const double wz = 0.5 * width * gl8_weight[q % 4];
        const double ring = sqrt(fmax(0.0, 1.0 - z * z));
        double row = 0.0;
        for (int j = 0; j < around; ++j) {
          const double phi = 2.0 * pi * ((double)j + 0.5) / (double)around;
          const double cx = ring * myb_cos(phi);
          const double cy = ring * myb_sin(phi);
          const double n[3] = {cx * frame[0] + cy * frame[3] + z * frame[6],
                               cx * frame[1] + cy * frame[4] + z * frame[7],
                               cx * frame[2] + cy * frame[5] + z * frame[8]};
          double dx[3];
          for (int i = 0; i < 3; ++i) {
            dx[i] = radius * n[i] - mean[i];
          }
          double quad = 0.0;
          double mu = 0.0;
          for (int i = 0; i < 3; ++i) {
            double ai = 0.0;
            double ki = 0.0;
            for (int m = 0; m < 3; ++m) {
              ai += inverse[3 * i + m] * dx[m];
              ki += k_gain[3 * i + m] * dx[m];
            }
            quad += dx[i] * ai;
            mu += n[i] * (mean[3 + i] + ki);
          }
          double s2 = 0.0;
          for (int i = 0; i < 3; ++i) {
            for (int m = 0; m < 3; ++m) {
              s2 += n[i] * c[3 * i + m] * n[m];
            }
          }
          row += myb_exp(-0.5 * quad) * inward(mu, sqrt(fmax(s2, 0.0)));
        }
        sum += wz * row * (2.0 * pi / (double)around);
      }
    }
    const double value = radius * radius * norm * sum;
    if (!isfinite(value)) {
      return MYB_E_DOMAIN;
    }
    if (level > 0 && (fabs(value - last) <= 1e-12 * value || (value < 1e-300 && last < 1e-300))) {
      *rate = value;
      return MYB_OK;
    }
    last = value;
    panels *= 2;
    around *= 2;
  }
  return MYB_E_NOCONV;
}

/* The two rules of Gauss and Kronrod over [a, b] from the rates at their 15 nodes, in the
 * order of the time. */
static void kronrod_sums(const double rates[15], double a, double b, double* kronrod,
                         double* gauss) {
  const double half = 0.5 * (b - a);
  double k = gk_weight[7] * rates[7];
  double g = g7_weight[3] * rates[7];
  for (int j = 0; j < 7; ++j) {
    const double pair = rates[j] + rates[14 - j];
    k += gk_weight[j] * pair;
    if (j % 2 == 1) {
      g += g7_weight[j / 2] * pair;
    }
  }
  *kronrod = half * k;
  *gauss = half * g;
}

/* The time of node j (0 .. 14, increasing) of the rule over [a, b]. */
static double kronrod_time(double a, double b, int j) {
  const double centre = 0.5 * (a + b);
  const double half = 0.5 * (b - a);
  return j < 7 ? centre - half * gk_node[j] : (j == 7 ? centre : centre + half * gk_node[14 - j]);
}

static myb_status rate_of(myb_relative_gaussian_function f, void* context, double t,
                          double radius, double* rate) {
  double mean[6];
  double covariance[36];
  const myb_status status = f(context, t, mean, covariance);
  if (status != MYB_OK) {
    return status;
  }
  return myb_collision_rate(mean, covariance, radius, rate);
}

myb_status myb_collision_probability_long(myb_relative_gaussian_function f, void* context,
                                          double t0, double t1, double step, double radius,
                                          double* probability, double* error) {
  MYB_CHECK(f != NULL && probability != NULL);
  if (!(t1 > t0) || !(step > 0.0) || !isfinite(t0) || !isfinite(t1) || !isfinite(step)) {
    return MYB_E_DOMAIN;
  }
  const double panels_wanted = ceil((t1 - t0) / step);
  if (!(panels_wanted < (double)PIECES)) {
    return MYB_E_DOMAIN;
  }
  const int panels = (int)panels_wanted;
  /* first pass: the panels, to know the size of the total */
  double total = 0.0;
  double spread = 0.0;
  int pieces = 0;
  double stack_a[PIECE_STACK];
  double stack_b[PIECE_STACK];
  for (int p = 0; p < panels; ++p) {
    const double a0 = t0 + (t1 - t0) * ((double)p / (double)panels);
    const double b0 = p + 1 == panels ? t1 : t0 + (t1 - t0) * ((double)(p + 1) / (double)panels);
    int depth = 0;
    stack_a[depth] = a0;
    stack_b[depth] = b0;
    ++depth;
    while (depth > 0) { /* bound: PIECES pieces in all */
      --depth;
      const double a = stack_a[depth];
      const double b = stack_b[depth];
      if (++pieces > PIECES) {
        return MYB_E_NOCONV;
      }
      double rates[15];
      for (int j = 0; j < 15; ++j) {
        const myb_status status = rate_of(f, context, kronrod_time(a, b, j), radius, &rates[j]);
        if (status != MYB_OK) {
          return status;
        }
      }
      double k = 0.0;
      double g = 0.0;
      kronrod_sums(rates, a, b, &k, &g);
      const double gap = fabs(k - g);
      /* halved while the two disagree against the total so far and this piece's own size */
      if (gap > 1e-10 * fmax(fabs(total) + fabs(k), 1e-300) && depth + 2 <= PIECE_STACK &&
          b - a > 1e-9 * (t1 - t0)) {
        const double middle = 0.5 * (a + b);
        stack_a[depth] = middle;
        stack_b[depth] = b;
        ++depth;
        stack_a[depth] = a;
        stack_b[depth] = middle;
        ++depth;
        continue;
      }
      total += k;
      spread += gap;
    }
  }
  *probability = total;
  if (error != NULL) {
    *error = spread;
  }
  return MYB_OK;
}

/* The relative Gaussian of two orbits at t: the difference of the states and the sum of
 * Phi P Phi^T. */
static myb_status orbits_gaussian(const pair_of_orbits* p, double t, double mean[6],
                                  double covariance[36]) {
  transitions phi;
  myb_status status = pair_at(p, t, mean, &phi);
  if (status != MYB_OK) {
    return status;
  }
  memset(covariance, 0, 36 * sizeof *covariance);
  for (int body = 0; body < 2; ++body) {
    const double* m = body == 0 ? phi.first : phi.second;
    const double* cov = body == 0 ? p->first_covariance : p->second_covariance;
    double mp[36];
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        double sum = 0.0;
        for (int k = 0; k < 6; ++k) {
          sum += m[6 * i + k] * cov[6 * k + j];
        }
        mp[6 * i + j] = sum;
      }
    }
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        double sum = 0.0;
        for (int k = 0; k < 6; ++k) {
          sum += mp[6 * i + k] * m[6 * j + k];
        }
        covariance[6 * i + j] += sum;
      }
    }
  }
  return MYB_OK;
}

myb_status myb_conjunction_long(myb_orbit* first, const double first_covariance[36],
                                myb_orbit* second, const double second_covariance[36],
                                double radius, double t1, double step, double* probability,
                                double* error) {
  MYB_CHECK(first != NULL && second != NULL && first_covariance != NULL &&
            second_covariance != NULL && probability != NULL);
  const double t0 = myb_orbit_time(first);
  if (myb_orbit_time(second) != t0 || !(step > 0.0) || !isfinite(t1) || !isfinite(step) ||
      t1 == t0) {
    return MYB_E_DOMAIN;
  }
  const double span = fabs(t1 - t0);
  const double panels_wanted = ceil(span / step);
  if (!(panels_wanted < (double)PIECES)) {
    return MYB_E_DOMAIN;
  }
  const int panels = (int)panels_wanted;
  const pair_of_orbits pair = {.first = first,
                               .second = second,
                               .first_covariance = first_covariance,
                               .second_covariance = second_covariance,
                               .radius = radius};
  double total = 0.0;
  double spread = 0.0;
  for (int p = 0; p < panels; ++p) {
    const double a = t0 + (t1 - t0) * ((double)p / (double)panels);
    const double b = p + 1 == panels ? t1 : t0 + (t1 - t0) * ((double)(p + 1) / (double)panels);
    double rates[15];
    for (int j = 0; j < 15; ++j) {
      /* the nodes in the order of the propagation, backward as well */
      const int node = t1 > t0 ? j : 14 - j;
      double mean[6];
      double covariance[36];
      myb_status status = orbits_gaussian(&pair, kronrod_time(a, b, node), mean, covariance);
      if (status == MYB_OK) {
        status = myb_collision_rate(mean, covariance, radius, &rates[node]);
      }
      if (status != MYB_OK) {
        return status;
      }
    }
    double k = 0.0;
    double g = 0.0;
    kronrod_sums(rates, a, b, &k, &g);
    total += fabs(k);
    spread += fabs(k - g);
  }
  /* the end of the propagation */
  double end[6];
  const myb_status status = pair_at(&pair, t1, end, NULL);
  if (status != MYB_OK) {
    return status;
  }
  *probability = total;
  if (error != NULL) {
    *error = spread;
  }
  return MYB_OK;
}
