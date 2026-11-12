/* SPDX-License-Identifier: BSD-3-Clause */
/* Contact of a body with a surface (myb/contact.h): the nearest point of a plane, an ellipsoid
 * or a mesh, the force of Hunt and Crossley with Coulomb's friction made smooth, and the
 * restitution of the law in closed form. */
#include "myb/contact.h"

#include <math.h> /* sqrt, fabs */
#include <string.h>

#include "core/base.h"
#include "force/contact_internal.h"
#include "myb/cmath.h"
#include "myb/geometry.h"
#include "num/root.h"
#include "num/vec3.h"

enum {
  SERIES_TERMS = 32,     /* of the series below, for arguments under SERIES_LIMIT: 0.25^30 */
  SOLVE_ITERATIONS = 200 /* of the bracketed Newton iterations: more than the bits of a double */
};

static const double SERIES_LIMIT = 0.25;

/* ---- Validity ------------------------------------------------------------------------------ */

static int finite_scalar(double x) {
  return x - x == 0.0;
}

static int valid_surface(const myb_contact_surface* s) {
  switch (s->shape) {
    case MYB_CONTACT_PLANE:
      return myb_vec3_finite(s->point) && myb_vec3_finite(s->normal) &&
             myb_vec3_norm(s->normal) > 0.0;
    case MYB_CONTACT_ELLIPSOID:
      return s->radii[0] > 0.0 && s->radii[1] > 0.0 && s->radii[2] > 0.0 &&
             myb_vec3_finite(s->radii);
    case MYB_CONTACT_MESH:
      return s->mesh != NULL;
    default:
      return 0;
  }
}

static int valid_law(const myb_contact_law* law) {
  return law->stiffness > 0.0 && finite_scalar(law->stiffness) && law->exponent >= 1.0 &&
         finite_scalar(law->exponent) && law->damping >= 0.0 && finite_scalar(law->damping) &&
         law->friction >= 0.0 && finite_scalar(law->friction) &&
         (law->friction == 0.0 || (law->slip > 0.0 && finite_scalar(law->slip)));
}

static int valid_point(const myb_contact_point* point) {
  return myb_vec3_finite(point->offset) && point->radius >= 0.0 && finite_scalar(point->radius);
}

int myb_internal_contact_valid(const myb_contact_surface* surface, const myb_contact_law* law,
                               const myb_contact_point* points, size_t count) {
  if (surface == NULL || law == NULL || points == NULL || count == 0 || !valid_surface(surface) ||
      !valid_law(law)) {
    return 0;
  }
  for (size_t k = 0; k < count; ++k) {
    if (!valid_point(&points[k])) {
      return 0;
    }
  }
  return 1;
}

/* ---- The nearest point --------------------------------------------------------------------- */

/* For a valid surface and a finite position. */
static myb_status nearest(const myb_contact_surface* s, const double position[3], double point[3],
                          double normal[3], double* height) {
  switch (s->shape) {
    case MYB_CONTACT_PLANE: {
      const double length = myb_vec3_norm(s->normal);
      for (int i = 0; i < 3; ++i) {
        normal[i] = s->normal[i] / length;
      }
      const double d[3] = {position[0] - s->point[0], position[1] - s->point[1],
                           position[2] - s->point[2]};
      const double h = myb_vec3_dot(normal, d);
      for (int i = 0; i < 3; ++i) {
        point[i] = position[i] - h * normal[i];
      }
      *height = h;
      return MYB_OK;
    }
    case MYB_CONTACT_ELLIPSOID: {
      const myb_status st = myb_ellipsoid_nearest(s->radii, position, point, height);
      if (st != MYB_OK) {
        return st;
      }
      myb_ellipsoid_normal(s->radii, point, normal);
      return MYB_OK;
    }
    case MYB_CONTACT_MESH:
      return myb_polyhedron_nearest(s->mesh, position, point, normal, height, NULL);
    default:
      MYB_UNREACHABLE();
  }
}

myb_status myb_contact_surface_nearest(const myb_contact_surface* surface, const double position[3],
                                       double point[3], double normal[3], double* height) {
  MYB_CHECK(surface != NULL && position != NULL);
  if (!valid_surface(surface) || !myb_vec3_finite(position)) {
    return MYB_E_DOMAIN;
  }
  double p[3];
  double n[3];
  double h = 0.0;
  const myb_status st = nearest(surface, position, p, n, &h);
  if (st != MYB_OK) {
    return st;
  }
  if (point != NULL) {
    memcpy(point, p, sizeof p);
  }
  if (normal != NULL) {
    memcpy(normal, n, sizeof n);
  }
  if (height != NULL) {
    *height = h;
  }
  return MYB_OK;
}

/* ---- The force ------------------------------------------------------------------------------- */

/* d^n for d > 0: exact for the usual exponents. */
static double depth_power(double d, double n) {
  if (n == 1.0) {
    return d;
  }
  if (n == 1.5) {
    return d * sqrt(d);
  }
  if (n == 2.0) {
    return d * d;
  }
  return myb_exp(n * myb_log(d));
}

myb_status myb_contact_force(const myb_contact_surface* surface, const myb_contact_law* law,
                             const myb_contact_point* points, size_t count, const double x[3],
                             const double v[3], const double attitude[9], const double omega[3],
                             myb_contact_result* result) {
  MYB_CHECK(surface != NULL && law != NULL && (points != NULL || count == 0) && x != NULL &&
            v != NULL && result != NULL);
  if (!myb_internal_contact_valid(surface, law, points, count) || !myb_vec3_finite(x) ||
      !myb_vec3_finite(v)) {
    return MYB_E_DOMAIN;
  }
  /* m: surface axes -> body axes; w: the angular velocity in the axes of the surface */
  double m[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  double w[3] = {0.0, 0.0, 0.0};
  if (attitude != NULL) {
    for (int i = 0; i < 9; ++i) {
      if (!finite_scalar(attitude[i])) {
        return MYB_E_DOMAIN;
      }
    }
    memcpy(m, attitude, sizeof m);
  }
  if (omega != NULL) {
    if (!myb_vec3_finite(omega)) {
      return MYB_E_DOMAIN;
    }
    for (int i = 0; i < 3; ++i) {
      w[i] = m[i] * omega[0] + m[3 + i] * omega[1] + m[6 + i] * omega[2];
    }
  }
  double force[3] = {0.0, 0.0, 0.0};
  double moment[3] = {0.0, 0.0, 0.0}; /* about the centre of mass, axes of the surface */
  double energy = 0.0;
  double gap = (double)INFINITY;
  size_t touching = 0;
  const double n = law->exponent;
  for (size_t k = 0; k < count; ++k) {
    const myb_contact_point* cp = &points[k];
    double centre[3];
    for (int i = 0; i < 3; ++i) {
      centre[i] =
          x[i] + (m[i] * cp->offset[0] + m[3 + i] * cp->offset[1] + m[6 + i] * cp->offset[2]);
    }
    double s[3];
    double normal[3];
    double h = 0.0;
    const myb_status st = nearest(surface, centre, s, normal, &h);
    if (st != MYB_OK) {
      return st;
    }
    const double clearance = h - cp->radius;
    gap = clearance < gap ? clearance : gap;
    if (!(clearance < 0.0)) {
      continue;
    }
    const double depth = -clearance;
    double lever[3];
    double turn[3];
    double u[3];
    for (int i = 0; i < 3; ++i) {
      lever[i] = s[i] - x[i];
    }
    myb_vec3_cross(w, lever, turn);
    for (int i = 0; i < 3; ++i) {
      u[i] = v[i] + turn[i]; /* the velocity of the point of the body at s */
    }
    const double along = myb_vec3_dot(u, normal);
    const double spring = law->stiffness * depth_power(depth, n);
    double fn = spring * (1.0 + law->damping * -along); /* the depth grows at -along */
    if (!(fn > 0.0)) {
      fn = 0.0; /* thrown off faster than 1 / chi: no adhesion */
    }
    double f[3];
    for (int i = 0; i < 3; ++i) {
      f[i] = fn * normal[i];
    }
    if (law->friction > 0.0) {
      double slide[3];
      for (int i = 0; i < 3; ++i) {
        slide[i] = u[i] - along * normal[i];
      }
      const double speed2 = myb_vec3_dot(slide, slide);
      const double scale = law->friction * fn / sqrt(speed2 + law->slip * law->slip);
      for (int i = 0; i < 3; ++i) {
        f[i] -= scale * slide[i];
      }
    }
    double arm[3];
    myb_vec3_cross(lever, f, arm);
    for (int i = 0; i < 3; ++i) {
      force[i] += f[i];
      moment[i] += arm[i];
    }
    energy += spring * depth / (n + 1.0);
    ++touching;
  }
  memcpy(result->force, force, sizeof force);
  for (int i = 0; i < 3; ++i) {
    result->torque[i] = m[3 * i] * moment[0] + m[3 * i + 1] * moment[1] + m[3 * i + 2] * moment[2];
  }
  result->energy = energy;
  result->gap = gap;
  result->touching = touching;
  return MYB_OK;
}

/* ---- The restitution of the law ---------------------------------------------------------------
 * With y = chi v_in and w = chi v_out the relation is c(y) = phi(w):
 *   c(y)   = y - ln(1 + y)  = sum_{k>=2} (-1)^k y^k / k,
 *   phi(w) = -w - ln(1 - w) = sum_{k>=2} w^k / k,       0 <= w < 1,
 * both near y^2 / 2 for small arguments, where the closed forms lose the digits that the
 * series keep. */

/* sum_{k>=2} s^k y^k / k for s = -1 (c) or +1 (phi), y < SERIES_LIMIT: Horner from the end. */
static double log_series(double y, double s) {
  double p = 0.0;
  for (int k = SERIES_TERMS + 1; k >= 2; --k) {
    const double sign = (k % 2 == 0) ? 1.0 : s;
    p = sign / (double)k + y * p;
  }
  return y * y * p;
}

static double excess_up(double y) {
  return y < SERIES_LIMIT ? log_series(y, -1.0) : y - myb_log1p(y);
}

static double excess_down(double w) {
  return w < SERIES_LIMIT ? log_series(w, 1.0) : -w - myb_log1p(-w);
}

typedef struct outgoing {
  double c; /* c(y) */
} outgoing;

/* phi(w) - c, increasing in w on [0, 1). */
static void outgoing_eval(const void* context, double w, double* f, double* df) {
  const outgoing* o = context;
  *f = excess_down(w) - o->c;
  *df = w / (1.0 - w);
}

myb_status myb_contact_restitution(double damping, double speed, double* restitution) {
  MYB_CHECK(restitution != NULL);
  const double y = damping * speed;
  if (!(damping >= 0.0) || !(speed >= 0.0) || !finite_scalar(y)) {
    return MYB_E_DOMAIN;
  }
  if (y == 0.0) {
    *restitution = 1.0;
    return MYB_OK;
  }
  const outgoing o = {excess_up(y)};
  const myb_internal_root_function function = {outgoing_eval, &o};
  /* the root is below y (phi(y) >= c(y)) and below 1 - exp(-(1 + c)) (phi there is
   * c + exp(-(1 + c))): Newton's iteration on the convex phi from there does not overshoot */
  const double hi = y < 1.0 ? y : 1.0;
  const double start = -myb_expm1(-(1.0 + o.c));
  double w = 0.0;
  const myb_status st = myb_internal_newton_bracketed(&function, 0.0, hi, start < hi ? start : hi,
                                                      SOLVE_ITERATIONS, &w);
  if (st != MYB_OK) {
    return st;
  }
  *restitution = w / y;
  return MYB_OK;
}

typedef struct incoming {
  double e;
  double log_e; /* ln e */
} incoming;

/* G(y) = phi(e y) - c(y) and its derivative: increasing from y* = (1 - e) / e, where it is
 * negative, to +infinity at 1 / e. For a small y the two series are taken as one, whose
 * coefficients e^k - (-1)^k keep 1 - e^k from cancellation. */
static void incoming_eval(const void* context, double y, double* f, double* df) {
  const incoming* in = context;
  const double e = in->e;
  if (y < SERIES_LIMIT) {
    double p = 0.0;
    for (int k = SERIES_TERMS + 1; k >= 2; --k) {
      const double b =
          (k % 2 == 0) ? myb_expm1((double)k * in->log_e) : myb_exp((double)k * in->log_e) + 1.0;
      p = b / (double)k + y * p;
    }
    *f = y * y * p;
  } else {
    *f = excess_down(e * y) - excess_up(y);
  }
  *df = e * e * y / (1.0 - e * y) - y / (1.0 + y);
}

myb_status myb_contact_damping(double restitution, double speed, double* damping) {
  MYB_CHECK(damping != NULL);
  const double e = restitution;
  if (!(e > 0.0 && e <= 1.0) || !(speed > 0.0) || !finite_scalar(speed)) {
    return MYB_E_DOMAIN;
  }
  if (e == 1.0) {
    *damping = 0.0;
    return MYB_OK;
  }
  const incoming in = {e, myb_log1p(e - 1.0)};
  const myb_internal_root_function function = {incoming_eval, &in};
  const double lo = (1.0 - e) / e;
  const double hi = 1.0 / e;
  double y = 0.0;
  /* near e = 1 the root is 3 (1 - e) / 2 */
  const myb_status st =
      myb_internal_newton_bracketed(&function, lo, hi, 1.5 * lo, SOLVE_ITERATIONS, &y);
  if (st != MYB_OK) {
    return st;
  }
  *damping = y / speed;
  return MYB_OK;
}
