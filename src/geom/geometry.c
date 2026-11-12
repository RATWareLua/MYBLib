/* SPDX-License-Identifier: BSD-3-Clause */
/* The geometry of a triaxial ellipsoid and the angles of illumination (myb/geometry.h). The
 * functions of the triangular mesh are with the mesh, in grav/polyhedron.c. */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "myb/cmath.h"
#include "myb/geometry.h"
#include "num/vec3.h"

enum { NEAREST_ITERATIONS = 100 };

static int valid_radii(const double radii[3]) {
  return radii[0] > 0.0 && radii[1] > 0.0 && radii[2] > 0.0 && isfinite(radii[0]) &&
         isfinite(radii[1]) && isfinite(radii[2]);
}

int myb_ray_ellipsoid(const double radii[3], const double origin[3], const double direction[3],
                      double* distance, double point[3]) {
  MYB_CHECK(radii != NULL && origin != NULL && direction != NULL && valid_radii(radii));
  if (!myb_vec3_finite(origin) || !myb_vec3_finite(direction)) {
    return 0;
  }
  const double length = myb_vec3_norm(direction);
  if (!(length > 0.0)) {
    return 0;
  }
  /* in coordinates divided by the semi-axes the ellipsoid is the unit sphere; d keeps the
   * parameter of the ray in metres */
  double o[3];
  double d[3];
  for (int i = 0; i < 3; ++i) {
    o[i] = origin[i] / radii[i];
    d[i] = direction[i] / length / radii[i];
  }
  const double a = myb_vec3_dot(d, d);
  const double middle = -myb_vec3_dot(o, d) / a; /* the point of the line nearest the centre */
  const double p[3] = {o[0] + middle * d[0], o[1] + middle * d[1], o[2] + middle * d[2]};
  const double h2 = myb_vec3_dot(p, p);
  if (!(h2 <= 1.0)) {
    return 0;
  }
  const double half = sqrt((1.0 - h2) / a);
  double s = middle - half;
  if (s < 0.0) {
    s = middle + half; /* from inside: where the ray leaves */
    if (s < 0.0) {
      return 0;
    }
  }
  if (distance != NULL) {
    *distance = s;
  }
  if (point != NULL) {
    for (int i = 0; i < 3; ++i) {
      point[i] = origin[i] + s * (direction[i] / length);
    }
  }
  return 1;
}

void myb_ellipsoid_normal(const double radii[3], const double point[3], double normal[3]) {
  MYB_CHECK(radii != NULL && point != NULL && normal != NULL && valid_radii(radii));
  /* the gradient of sum (x_i / a_i)^2, scaled first so that the squares cannot overflow */
  double g[3];
  for (int i = 0; i < 3; ++i) {
    g[i] = point[i] / radii[i] / radii[i];
  }
  const double scale = fmax(fabs(g[0]), fmax(fabs(g[1]), fabs(g[2])));
  if (!(scale > 0.0) || !isfinite(scale)) {
    normal[0] = normal[1] = normal[2] = (double)NAN;
    return;
  }
  for (int i = 0; i < 3; ++i) {
    g[i] /= scale;
  }
  const double norm = myb_vec3_norm(g);
  for (int i = 0; i < 3; ++i) {
    normal[i] = g[i] / norm;
  }
}

/* The nearest point x of the ellipsoid to y has y - x along the normal at x:
 *   x_i = a_i^2 y_i / (t + a_i^2),   F(t) = sum (a_i y_i / (t + a_i^2))^2 - 1 = 0,
 * the largest root t (D. Eberly, "Distance from a point to an ellipse, an ellipsoid, or a
 * hyperellipsoid", 2013). On t > -c^2, c the shortest semi-axis, F falls and is convex, so
 * Newton's iteration from a point where F >= 0 climbs to the root without passing it. Such a
 * point: each term alone gives F(t) >= (a_i y_i / (t + a_i^2))^2 - 1, zero at
 * t_i = a_i |y_i| - a_i^2; outside the ellipsoid also F(0) >= 0.
 *
 * A position with no component along the shortest axes has no pole of F at -c^2. If F is
 * still >= 0 there, the root lies above as before. If not, the position is inside the evolute:
 * the nearest points leave the plane, t = -c^2, and what is left of the unit sum goes to the
 * shortest axes — two points, or a circle of them; the one on the positive side of the first
 * of those axes is given. */
myb_status myb_ellipsoid_nearest(const double radii[3], const double position[3],
                                 double point[3], double* height) {
  MYB_CHECK(radii != NULL && position != NULL && valid_radii(radii));
  if (!myb_vec3_finite(position)) {
    return MYB_E_DOMAIN;
  }
  double a2[3];
  double q = 0.0;
  double shortest = radii[0] * radii[0];
  for (int i = 0; i < 3; ++i) {
    a2[i] = radii[i] * radii[i];
    const double z = position[i] / radii[i];
    q += z * z;
    shortest = fmin(shortest, a2[i]);
  }
  double t = -shortest;
  int pole = 0; /* a component along a shortest axis: F grows without bound at -c^2 */
  for (int i = 0; i < 3; ++i) {
    if (position[i] != 0.0) {
      t = fmax(t, radii[i] * fabs(position[i]) - a2[i]);
      pole = pole || a2[i] == shortest;
    }
  }
  if (q >= 1.0) {
    t = fmax(t, 0.0);
  }
  double x[3] = {0.0, 0.0, 0.0};
  if (!pole && t == -shortest) {
    double rest = 1.0; /* 1 - sum over the longer axes of (x_i / a_i)^2 at t = -c^2 */
    for (int i = 0; i < 3; ++i) {
      if (a2[i] != shortest) {
        x[i] = a2[i] * position[i] / (a2[i] - shortest);
        rest -= (x[i] / radii[i]) * (x[i] / radii[i]);
      }
    }
    if (rest > 0.0) {
      for (int i = 0; i < 3; ++i) {
        if (a2[i] == shortest) {
          x[i] = radii[i] * sqrt(rest);
          break;
        }
      }
      if (point != NULL) {
        memcpy(point, x, sizeof x);
      }
      if (height != NULL) {
        const double d[3] = {position[0] - x[0], position[1] - x[1], position[2] - x[2]};
        *height = -myb_vec3_norm(d);
      }
      return MYB_OK;
    }
  }
  int converged = 0;
  double last_step = 0.0;
  /* bound: NEAREST_ITERATIONS; the iteration is monotone and ends when it stops moving. Near a
   * sphere F may round to a positive value far below the resolution of 1, whose Newton step
   * is then below the resolution of t + a_i^2, the only way t enters F: F comes back nearly
   * the same and t creeps by such steps until the budget is spent. A t that still moves by
   * less than the resolution of the smallest t + a_i^2 is the root to the resolution of what
   * it enters, and is kept */
  for (int iteration = 0; iteration < NEAREST_ITERATIONS && !converged; ++iteration) {
    double f = -1.0;
    double slope = 0.0;
    for (int i = 0; i < 3; ++i) {
      if (position[i] != 0.0) {
        const double w = radii[i] * position[i] / (t + a2[i]);
        f += w * w;
        slope -= 2.0 * w * w / (t + a2[i]);
      }
    }
    const double next = slope < 0.0 ? t - f / slope : t;
    converged = !(next > t);
    if (!converged) {
      last_step = next - t;
      t = next;
    }
  }
  if (!converged && !(last_step <= 0x1p-51 * (t + shortest))) {
    return MYB_E_NOCONV;
  }
  double offset2 = 0.0;
  for (int i = 0; i < 3; ++i) {
    if (position[i] != 0.0) {
      const double w = position[i] / (t + a2[i]);
      offset2 += w * w;
      x[i] = a2[i] * w;
    }
  }
  if (point != NULL) {
    memcpy(point, x, sizeof x);
  }
  if (height != NULL) {
    *height = t * sqrt(offset2); /* y - x = t y_i / (t + a_i^2) */
  }
  return MYB_OK;
}

/* In coordinates divided by the semi-axes the ellipsoid is the unit sphere and the limb seen
 * from o the circle of the plane p.o = 1: centre o / |o|^2, radius sqrt(1 - 1 / |o|^2).
 * Stretched back, two perpendicular radii of that circle are conjugate half-diameters of
 * the limb; its own axes are those turned by half the angle atan2(2 u.v, u.u - v.v). */
myb_status myb_ellipsoid_limb(const double radii[3], const double observer[3], double centre[3],
                              double major[3], double minor[3]) {
  MYB_CHECK(radii != NULL && observer != NULL && centre != NULL && major != NULL &&
            minor != NULL && valid_radii(radii));
  if (!myb_vec3_finite(observer)) {
    return MYB_E_DOMAIN;
  }
  const double o[3] = {observer[0] / radii[0], observer[1] / radii[1], observer[2] / radii[2]};
  const double d2 = myb_vec3_dot(o, o);
  if (!(d2 > 1.0) || !isfinite(d2)) {
    return MYB_E_DOMAIN;
  }
  const double d = sqrt(d2);
  const double z[3] = {o[0] / d, o[1] / d, o[2] / d};
  const double helper[3] = {fabs(z[0]) < 0.6 ? 1.0 : 0.0, fabs(z[0]) < 0.6 ? 0.0 : 1.0, 0.0};
  double e1[3];
  double e2[3];
  myb_vec3_cross(helper, z, e1);
  const double n1 = myb_vec3_norm(e1);
  for (int i = 0; i < 3; ++i) {
    e1[i] /= n1;
  }
  myb_vec3_cross(z, e1, e2);
  const double radius = sqrt((d - 1.0) * (d + 1.0)) / d;
  double u[3];
  double v[3];
  for (int i = 0; i < 3; ++i) {
    centre[i] = radii[i] * o[i] / d2;
    u[i] = radii[i] * radius * e1[i];
    v[i] = radii[i] * radius * e2[i];
  }
  const double turn =
      0.5 * myb_atan2(2.0 * myb_vec3_dot(u, v), myb_vec3_dot(u, u) - myb_vec3_dot(v, v));
  const double ct = myb_cos(turn);
  const double st = myb_sin(turn);
  for (int i = 0; i < 3; ++i) {
    major[i] = ct * u[i] + st * v[i];
    minor[i] = -st * u[i] + ct * v[i];
  }
  if (myb_vec3_norm(minor) > myb_vec3_norm(major)) {
    for (int i = 0; i < 3; ++i) {
      const double swap = major[i];
      major[i] = minor[i];
      minor[i] = -swap;
    }
  }
  return MYB_OK;
}

/* The point of the ellipsoid in the direction w (unit) and how far the centre of the source
 * stands over the plane tangent there, less what the edge wanted asks. */
static double horizon_excess(const double radii[3], const double w[3], const double source[3],
                             double level, double point[3]) {
  double q = 0.0;
  for (int i = 0; i < 3; ++i) {
    q += (w[i] / radii[i]) * (w[i] / radii[i]);
  }
  const double reach = 1.0 / sqrt(q);
  double normal[3];
  for (int i = 0; i < 3; ++i) {
    point[i] = reach * w[i];
  }
  myb_ellipsoid_normal(radii, point, normal);
  const double to_source[3] = {source[0] - point[0], source[1] - point[1], source[2] - point[2]};
  return myb_vec3_dot(normal, to_source) - level;
}

myb_status myb_ellipsoid_terminator(const double radii[3], const double source[3],
                                    double source_radius, int penumbral, size_t count,
                                    double* points) {
  MYB_CHECK(radii != NULL && source != NULL && (points != NULL || count == 0) &&
            valid_radii(radii));
  double q = 0.0;
  for (int i = 0; i < 3; ++i) {
    q += (source[i] / radii[i]) * (source[i] / radii[i]);
  }
  if (!myb_vec3_finite(source) || !(q > 1.0) || !(source_radius >= 0.0) ||
      !isfinite(source_radius)) {
    return MYB_E_DOMAIN;
  }
  const double sn = myb_vec3_norm(source);
  const double axis[3] = {source[0] / sn, source[1] / sn, source[2] / sn};
  /* across the axis: towards +z, or +x for a source on the z axis */
  const double helper[3] = {fabs(axis[2]) < 1.0 - 1e-12 ? 0.0 : 1.0, 0.0,
                            fabs(axis[2]) < 1.0 - 1e-12 ? 1.0 : 0.0};
  const double along = myb_vec3_dot(helper, axis);
  double e1[3];
  double e2[3];
  for (int i = 0; i < 3; ++i) {
    e1[i] = helper[i] - along * axis[i];
  }
  const double n1 = myb_vec3_norm(e1);
  for (int i = 0; i < 3; ++i) {
    e1[i] /= n1;
  }
  myb_vec3_cross(axis, e1, e2);
  /* the centre of the source over the horizon by its radius: all of it is seen; under it by
   * the radius: none */
  const double level = penumbral ? source_radius : -source_radius;
  for (size_t k = 0; k < count; ++k) {
    const double azimuth = 6.283185307179586 * (double)k / (double)count;
    const double ca = myb_cos(azimuth);
    const double sa = myb_sin(azimuth);
    double across[3];
    for (int i = 0; i < 3; ++i) {
      across[i] = ca * e1[i] + sa * e2[i];
    }
    /* from under the source (the excess is positive) to the far side (negative): halving;
     * bound: 80 halvings */
    double lo = 0.0;
    double hi = 3.141592653589793;
    double point[3] = {0.0, 0.0, 0.0};
    for (int it = 0; it < 80; ++it) {
      const double mid = 0.5 * (lo + hi);
      const double cm = myb_cos(mid);
      const double sm = myb_sin(mid);
      const double w[3] = {cm * axis[0] + sm * across[0], cm * axis[1] + sm * across[1],
                           cm * axis[2] + sm * across[2]};
      if (horizon_excess(radii, w, source, level, point) > 0.0) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    memcpy(points + 3 * k, point, sizeof point);
  }
  return MYB_OK;
}

static double angle_between(const double u[3], const double v[3]) {
  double c[3];
  myb_vec3_cross(u, v, c);
  return myb_atan2(myb_vec3_norm(c), myb_vec3_dot(u, v));
}

myb_status myb_illumination_angles(const double point[3], const double normal[3],
                                   const double source[3], const double observer[3],
                                   double* phase, double* incidence, double* emission) {
  MYB_CHECK(point != NULL && normal != NULL && source != NULL && observer != NULL);
  double s[3];
  double o[3];
  for (int i = 0; i < 3; ++i) {
    s[i] = source[i] - point[i];
    o[i] = observer[i] - point[i];
  }
  if (!myb_vec3_finite(s) || !myb_vec3_finite(o) || !myb_vec3_finite(normal) ||
      !(myb_vec3_norm(s) > 0.0) || !(myb_vec3_norm(o) > 0.0) || !(myb_vec3_norm(normal) > 0.0)) {
    return MYB_E_DOMAIN;
  }
  if (phase != NULL) {
    *phase = angle_between(s, o);
  }
  if (incidence != NULL) {
    *incidence = angle_between(normal, s);
  }
  if (emission != NULL) {
    *emission = angle_between(normal, o);
  }
  return MYB_OK;
}
