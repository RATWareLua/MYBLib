/* SPDX-License-Identifier: BSD-3-Clause */
/* Radiation pressure: spheres, plates, shadows, the radiation of a planet. */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "myb/cmath.h"
#include "myb/force.h"
#include "num/gauss24.h"
#include "num/vec3.h"

static const double pi = 3.141592653589793;
static const double speed_of_light = 299792458.0; /* m s^-1, exact (SI) */

double myb_radiation_flux(double luminosity, double distance) {
  return luminosity / (4.0 * pi * distance * distance);
}

void myb_force_radiation_sphere(double pressure, double cr, double area_over_mass,
                                const double u[3], double a[3]) {
  MYB_CHECK(u != NULL && a != NULL);
  const double k = pressure * cr * area_over_mass;
  a[0] += k * u[0];
  a[1] += k * u[1];
  a[2] += k * u[2];
}

myb_status myb_force_radiation_plates(double pressure, double mass, const myb_plate* plates,
                                      size_t count, const double u[3], double a[3]) {
  MYB_CHECK(u != NULL && a != NULL && (plates != NULL || count == 0));
  if (!(mass > 0.0) || !isfinite(mass) || !isfinite(pressure)) {
    return MYB_E_DOMAIN;
  }
  for (size_t k = 0; k < count; ++k) {
    const myb_plate* p = &plates[k];
    const double nn = myb_vec3_dot(p->normal, p->normal);
    if (!(p->area >= 0.0) || !isfinite(p->area) || !(fabs(nn - 1.0) <= 2e-12) ||
        !(p->specular >= 0.0 && p->specular <= 1.0) ||
        !(p->diffuse >= 0.0 && p->diffuse <= 1.0) || p->specular + p->diffuse > 1.0) {
      return MYB_E_DOMAIN;
    }
  }
  double f[3] = {0.0, 0.0, 0.0};
  for (size_t k = 0; k < count; ++k) {
    const myb_plate* p = &plates[k];
    const double c = -myb_vec3_dot(p->normal, u);
    if (!(c > 0.0)) {
      continue; /* lit from behind or edge-on */
    }
    const double absorbed = 1.0 - p->specular - p->diffuse;
    double normal = 2.0 * (p->specular * c + p->diffuse / 3.0);
    if (p->thermal) {
      normal += 2.0 / 3.0 * absorbed;
    }
    const double k_u = pressure * p->area * c * (1.0 - p->specular);
    const double k_n = -pressure * p->area * c * normal;
    for (int i = 0; i < 3; ++i) {
      f[i] += k_u * u[i] + k_n * p->normal[i];
    }
  }
  for (int i = 0; i < 3; ++i) {
    a[i] += f[i] / mass;
  }
  return MYB_OK;
}

/* ---- Shadows --------------------------------------------------------------------------------- */

typedef struct shadow_geometry {
  double a;  /* angular radius of the source */
  double b;  /* angular radius of the occulter */
  double c;  /* angular separation of their centres */
  double sa; /* sin(a) */
  double sc; /* sin(c) */
  double u1;
  double u2;
} shadow_geometry;

/* Brightness of the source's disk at angular distance rho from its centre. */
static double brightness(const shadow_geometry* g, double rho) {
  if (g->u1 == 0.0 && g->u2 == 0.0) {
    return 1.0;
  }
  const double x = myb_sin(rho) / g->sa;
  const double mu = sqrt(fmax(0.0, (1.0 - x) * (1.0 + x)));
  const double v = 1.0 - mu;
  return 1.0 - g->u1 * v - g->u2 * v * v;
}

/* Half of the angle of the circle of radius rho around the source's centre that lies inside the
 * occulter's cap: points at azimuth phi are inside when cos(phi) > q, with
 *   (1 - q) / 2 = sin((b + c - rho) / 2) sin((b - c + rho) / 2) / (sin(c) sin(rho)),
 * written without the cancellation of cos(b) - cos(c) cos(rho). */
static double occulted_half_angle(const shadow_geometry* g, double rho) {
  const double sr = myb_sin(rho);
  if (g->sc * sr == 0.0) {
    return rho + g->c < g->b ? pi : 0.0;
  }
  const double h =
      myb_sin(0.5 * (g->b + g->c - rho)) * myb_sin(0.5 * (g->b - g->c + rho)) / (g->sc * sr);
  if (h <= 0.0) {
    return 0.0;
  }
  if (h >= 1.0) {
    return pi;
  }
  return 2.0 * myb_asin(sqrt(h));
}

/* Integrals over [lo, hi] of I(rho) sin(rho) and I(rho) sin(rho) 2 phi_in(rho), with
 * rho = lo + (hi - lo)(1 - cos t)/2 so that square-root behaviour at both ends (the kinks of
 * phi_in, the limb of the source) becomes smooth in t. */
static void integrate(const shadow_geometry* g, double lo, double hi, double* total,
                      double* blocked) {
  const double half = 0.5 * (hi - lo);
  for (int k = 0; k < 24; ++k) {
    const double x = k < 12 ? myb_gauss24_x[k] : -myb_gauss24_x[k - 12];
    const double w = myb_gauss24_w[k < 12 ? k : k - 12];
    const double t = 0.5 * pi * (x + 1.0);
    const double rho = lo + half * (1.0 - myb_cos(t));
    const double jac = half * myb_sin(t) * 0.5 * pi;
    const double base = w * jac * brightness(g, rho) * myb_sin(rho);
    *total += 2.0 * pi * base;
    *blocked += 2.0 * occulted_half_angle(g, rho) * base;
  }
}

double myb_shadow_fraction(const double source[3], double source_radius,
                           const double occulter[3], double occulter_radius,
                           const double limb[2]) {
  MYB_CHECK(source != NULL && occulter != NULL);
  const double ds = myb_vec3_norm(source);
  const double d_o = myb_vec3_norm(occulter);
  if (d_o <= occulter_radius) {
    return 0.0;
  }
  if (!(ds > source_radius) || d_o > ds) {
    return 1.0;
  }
  double cross[3];
  myb_vec3_cross(source, occulter, cross);
  shadow_geometry g;
  g.a = myb_asin(source_radius / ds);
  g.b = myb_asin(occulter_radius / d_o);
  g.c = myb_atan2(myb_vec3_norm(cross), myb_vec3_dot(source, occulter));
  g.sa = myb_sin(g.a);
  g.sc = myb_sin(g.c);
  g.u1 = limb != NULL ? limb[0] : 0.0;
  g.u2 = limb != NULL ? limb[1] : 0.0;
  if (g.c >= g.a + g.b) {
    return 1.0;
  }
  if (g.c + g.a <= g.b) {
    return 0.0;
  }
  /* kinks of phi_in inside (0, a): |c - b| <= c + b */
  const double cuts[2] = {fabs(g.c - g.b), g.c + g.b};
  int n = 0;
  double points[4];
  points[n++] = 0.0;
  for (int k = 0; k < 2; ++k) {
    if (cuts[k] > points[n - 1] && cuts[k] < g.a) {
      points[n++] = cuts[k];
    }
  }
  points[n++] = g.a;
  double total = 0.0;
  double blocked = 0.0;
  for (int k = 0; k + 1 < n; ++k) {
    integrate(&g, points[k], points[k + 1], &total, &blocked);
  }
  const double f = 1.0 - blocked / total;
  return f < 0.0 ? 0.0 : (f > 1.0 ? 1.0 : f);
}

/* ---- Radiation of a planet ------------------------------------------------------------------- */

myb_status myb_force_planet_radiation(const myb_planet_radiation* planet, const double r[3],
                                      const double sun[3], double solar_flux, double phase,
                                      double cr_area_over_mass, int rings, double a[3]) {
  MYB_CHECK(planet != NULL && r != NULL && sun != NULL && a != NULL);
  const double R = planet->radius;
  const double rr = myb_vec3_norm(r);
  const double ss = myb_vec3_norm(sun);
  if (rings < 1 || rings > 50 || !(R > 0.0) || !(rr > R) || !isfinite(rr) || !(ss > 0.0) ||
      !isfinite(ss) || !isfinite(solar_flux) || !isfinite(phase) ||
      !isfinite(cr_area_over_mass)) {
    return MYB_E_DOMAIN;
  }
  const double* al = planet->albedo;
  const double* em = planet->emissivity;
  const double cp = myb_cos(phase);
  const double sp = myb_sin(phase);
  const double a1 = al[1] + al[2] * cp + al[3] * sp;
  const double e1 = em[1] + em[2] * cp + em[3] * sp;
  /* basis: z along r, x and y across */
  double ez[3] = {r[0] / rr, r[1] / rr, r[2] / rr};
  double ex[3];
  const double helper[3] = {fabs(ez[0]) < 0.9 ? 1.0 : 0.0, fabs(ez[0]) < 0.9 ? 0.0 : 1.0, 0.0};
  myb_vec3_cross(helper, ez, ex);
  const double nx = myb_vec3_norm(ex);
  for (int i = 0; i < 3; ++i) {
    ex[i] /= nx;
  }
  double ey[3];
  myb_vec3_cross(ez, ex, ey);
  const double s_hat[3] = {sun[0] / ss, sun[1] / ss, sun[2] / ss};
  const double cap = myb_acos(R / rr); /* central angle of the visible cap */
  const double step = cap / (2.0 * rings + 1.0);
  double sum[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k <= rings; ++k) {
    const double lo = k == 0 ? 0.0 : (2.0 * k - 1.0) * step;
    const double hi = (2.0 * k + 1.0) * step;
    const double mid = k == 0 ? 0.0 : 2.0 * k * step;
    const int segments = k == 0 ? 1 : 6 * k;
    const double area = R * R * (myb_cos(lo) - myb_cos(hi)) * (2.0 * pi / segments);
    const double cg = myb_cos(mid);
    const double sg = myb_sin(mid);
    /* azimuths (j + 1/2) 2 pi / segments by rotation of the phasor: 300 products at most,
     * an error of a few 1e-14 rad, far below the discretisation */
    const double step_angle = 2.0 * pi / segments;
    const double cd = myb_cos(step_angle);
    const double sd = myb_sin(step_angle);
    double c = myb_cos(0.5 * step_angle);
    double s = myb_sin(0.5 * step_angle);
    for (int j = 0; j < segments; ++j) {
      if (j > 0) {
        const double cn = c * cd - s * sd;
        s = s * cd + c * sd;
        c = cn;
      }
      double n[3];
      for (int i = 0; i < 3; ++i) {
        n[i] = cg * ez[i] + sg * (c * ex[i] + s * ey[i]);
      }
      double d[3] = {r[0] - R * n[0], r[1] - R * n[1], r[2] - R * n[2]};
      const double dist = myb_vec3_norm(d);
      for (int i = 0; i < 3; ++i) {
        d[i] /= dist;
      }
      const double cos_view = myb_vec3_dot(n, d);
      if (!(cos_view > 0.0)) {
        continue;
      }
      const double sin_lat = n[2];
      const double p2 = 0.5 * (3.0 * sin_lat * sin_lat - 1.0);
      const double albedo = al[0] + a1 * sin_lat + al[4] * p2;
      const double emissivity = em[0] + e1 * sin_lat + em[4] * p2;
      const double cos_sun = myb_vec3_dot(n, s_hat);
      const double exitance =
          solar_flux * (albedo * (cos_sun > 0.0 ? cos_sun : 0.0) + 0.25 * emissivity);
      const double irradiance = exitance * cos_view * area / (pi * dist * dist);
      for (int i = 0; i < 3; ++i) {
        sum[i] += irradiance * d[i];
      }
    }
  }
  const double k = cr_area_over_mass / speed_of_light;
  for (int i = 0; i < 3; ++i) {
    a[i] += k * sum[i];
  }
  return MYB_OK;
}

double myb_shadow_margin(const double source[3], double source_radius,
                         const double occulter[3], double occulter_radius, int inner) {
  MYB_CHECK(source != NULL && occulter != NULL);
  const double ds = myb_vec3_norm(source);
  const double d_o = myb_vec3_norm(occulter);
  if (d_o <= occulter_radius) {
    return -1.0;
  }
  if (!(ds > source_radius) || d_o > ds) {
    return 1.0;
  }
  double cross[3];
  myb_vec3_cross(source, occulter, cross);
  const double a = myb_asin(source_radius / ds);
  const double b = myb_asin(occulter_radius / d_o);
  const double c = myb_atan2(myb_vec3_norm(cross), myb_vec3_dot(source, occulter));
  return inner ? c - fabs(b - a) : c - (a + b);
}

/* ---- The gradient of the radiation of a planet ------------------------------------------------- */

myb_status myb_force_planet_radiation_gradient(const myb_planet_radiation* planet,
                                               const double r[3], const double sun[3],
                                               double solar_flux, double phase,
                                               double cr_area_over_mass, int rings,
                                               double dadr[9]) {
  MYB_CHECK(planet != NULL && r != NULL && sun != NULL && dadr != NULL);
  const double R = planet->radius;
  const double rr = myb_vec3_norm(r);
  const double ss = myb_vec3_norm(sun);
  if (rings < 1 || rings > 50 || !(R > 0.0) || !(rr > R) || !isfinite(rr) || !(ss > 0.0) ||
      !isfinite(ss) || !isfinite(solar_flux) || !isfinite(phase) ||
      !isfinite(cr_area_over_mass)) {
    return MYB_E_DOMAIN;
  }
  const double* al = planet->albedo;
  const double* em = planet->emissivity;
  const double cp = myb_cos(phase);
  const double sp = myb_sin(phase);
  const double a1 = al[1] + al[2] * cp + al[3] * sp;
  const double e1 = em[1] + em[2] * cp + em[3] * sp;
  double ez[3] = {r[0] / rr, r[1] / rr, r[2] / rr};
  double ex[3];
  const double helper[3] = {fabs(ez[0]) < 0.9 ? 1.0 : 0.0, fabs(ez[0]) < 0.9 ? 0.0 : 1.0, 0.0};
  myb_vec3_cross(helper, ez, ex);
  const double nx = myb_vec3_norm(ex);
  for (int i = 0; i < 3; ++i) {
    ex[i] /= nx;
  }
  double ey[3];
  myb_vec3_cross(ez, ex, ey);
  const double s_hat[3] = {sun[0] / ss, sun[1] / ss, sun[2] / ss};
  const double cap = myb_acos(R / rr);
  const double step = cap / (2.0 * rings + 1.0);
  double sum[9] = {0.0};
  for (int k = 0; k <= rings; ++k) {
    const double lo = k == 0 ? 0.0 : (2.0 * k - 1.0) * step;
    const double hi = (2.0 * k + 1.0) * step;
    const double mid = k == 0 ? 0.0 : 2.0 * k * step;
    const int segments = k == 0 ? 1 : 6 * k;
    const double area = R * R * (myb_cos(lo) - myb_cos(hi)) * (2.0 * pi / segments);
    const double cg = myb_cos(mid);
    const double sg = myb_sin(mid);
    const double step_angle = 2.0 * pi / segments;
    const double cd = myb_cos(step_angle);
    const double sd = myb_sin(step_angle);
    double c = myb_cos(0.5 * step_angle);
    double s = myb_sin(0.5 * step_angle);
    for (int j = 0; j < segments; ++j) {
      if (j > 0) {
        const double cn = c * cd - s * sd;
        s = s * cd + c * sd;
        c = cn;
      }
      double n[3];
      for (int i = 0; i < 3; ++i) {
        n[i] = cg * ez[i] + sg * (c * ex[i] + s * ey[i]);
      }
      double d[3] = {r[0] - R * n[0], r[1] - R * n[1], r[2] - R * n[2]};
      const double dist = myb_vec3_norm(d);
      for (int i = 0; i < 3; ++i) {
        d[i] /= dist;
      }
      const double cos_view = myb_vec3_dot(n, d);
      if (!(cos_view > 0.0)) {
        continue;
      }
      const double sin_lat = n[2];
      const double p2 = 0.5 * (3.0 * sin_lat * sin_lat - 1.0);
      const double albedo = al[0] + a1 * sin_lat + al[4] * p2;
      const double emissivity = em[0] + e1 * sin_lat + em[4] * p2;
      const double cos_sun = myb_vec3_dot(n, s_hat);
      const double exitance =
          solar_flux * (albedo * (cos_sun > 0.0 ? cos_sun : 0.0) + 0.25 * emissivity);
      /* of E (n.D) D / (pi |D|^4), D = r - R n:
       * E / (pi |D|^3) [ d n^T + cos (1 - 4 d d^T) ] */
      const double weight = exitance * area / (pi * dist * dist * dist);
      for (int i = 0; i < 3; ++i) {
        for (int q = 0; q < 3; ++q) {
          sum[3 * i + q] +=
              weight * (d[i] * n[q] + cos_view * ((i == q ? 1.0 : 0.0) - 4.0 * d[i] * d[q]));
        }
      }
    }
  }
  const double k = cr_area_over_mass / speed_of_light;
  for (int i = 0; i < 9; ++i) {
    dadr[i] = k * sum[i];
  }
  return MYB_OK;
}

/* ---- Several occulters at once ----------------------------------------------------------------- */

enum {
  MAX_OCCULTERS = 8,
  MAX_CUTS = 2 + 2 * MAX_OCCULTERS + MAX_OCCULTERS * (MAX_OCCULTERS - 1)
};

typedef struct occulter_view {
  shadow_geometry g; /* a, sa, u1, u2 of the source; b, c, sc of this occulter */
  double azimuth;    /* of its centre about the centre of the source */
  double unit[3];
} occulter_view;

/* The length of the union of arcs on a circle, each given by the azimuth of its middle and
 * its half-length (0 .. pi). */
static double arcs_union(const double* middle, const double* half, size_t count) {
  double lo[2 * MAX_OCCULTERS];
  double hi[2 * MAX_OCCULTERS];
  size_t n = 0;
  for (size_t k = 0; k < count; ++k) {
    if (!(half[k] > 0.0)) {
      continue;
    }
    if (half[k] >= pi) {
      return 2.0 * pi;
    }
    const double m = middle[k] - 2.0 * pi * floor(middle[k] / (2.0 * pi));
    const double a = m - half[k];
    const double b = m + half[k];
    if (a < 0.0) {
      lo[n] = a + 2.0 * pi;
      hi[n++] = 2.0 * pi;
      lo[n] = 0.0;
      hi[n++] = b;
    } else if (b > 2.0 * pi) {
      lo[n] = a;
      hi[n++] = 2.0 * pi;
      lo[n] = 0.0;
      hi[n++] = b - 2.0 * pi;
    } else {
      lo[n] = a;
      hi[n++] = b;
    }
  }
  /* by their beginnings, then merged; bound: 16 intervals */
  for (size_t i = 1; i < n; ++i) {
    const double l = lo[i];
    const double h = hi[i];
    size_t j = i;
    while (j > 0 && lo[j - 1] > l) {
      lo[j] = lo[j - 1];
      hi[j] = hi[j - 1];
      --j;
    }
    lo[j] = l;
    hi[j] = h;
  }
  double total = 0.0;
  double end = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double from = lo[i] > end ? lo[i] : end;
    if (hi[i] > from) {
      total += hi[i] - from;
      end = hi[i];
    }
  }
  return total < 2.0 * pi ? total : 2.0 * pi;
}

static void integrate_many(const occulter_view* views, size_t count, double lo, double hi,
                           double* total, double* blocked) {
  const shadow_geometry* g = &views[0].g;
  const double half = 0.5 * (hi - lo);
  for (int k = 0; k < 24; ++k) {
    const double x = k < 12 ? myb_gauss24_x[k] : -myb_gauss24_x[k - 12];
    const double w = myb_gauss24_w[k < 12 ? k : k - 12];
    const double t = 0.5 * pi * (x + 1.0);
    const double rho = lo + half * (1.0 - myb_cos(t));
    const double jac = half * myb_sin(t) * 0.5 * pi;
    const double base = w * jac * brightness(g, rho) * myb_sin(rho);
    double middle[MAX_OCCULTERS];
    double arc[MAX_OCCULTERS];
    for (size_t q = 0; q < count; ++q) {
      middle[q] = views[q].azimuth;
      arc[q] = occulted_half_angle(&views[q].g, rho);
    }
    *total += 2.0 * pi * base;
    *blocked += arcs_union(middle, arc, count) * base;
  }
}

double myb_shadow_fraction_many(const double source[3], double source_radius,
                                const double* occulters, const double* radii, size_t count,
                                const double limb[2]) {
  MYB_CHECK(source != NULL && (count == 0 || (occulters != NULL && radii != NULL)) &&
            count <= MAX_OCCULTERS);
  const double ds = myb_vec3_norm(source);
  if (!(ds > source_radius)) {
    /* as myb_shadow_fraction: only an object inside an occulter is in the dark */
    for (size_t k = 0; k < count; ++k) {
      if (myb_vec3_norm(occulters + 3 * k) <= radii[k]) {
        return 0.0;
      }
    }
    return 1.0;
  }
  const double a = myb_asin(source_radius / ds);
  const double s_hat[3] = {source[0] / ds, source[1] / ds, source[2] / ds};
  /* two directions across the direction to the source */
  const double helper[3] = {fabs(s_hat[0]) < 0.9 ? 1.0 : 0.0, fabs(s_hat[0]) < 0.9 ? 0.0 : 1.0,
                            0.0};
  double e1[3];
  double e2[3];
  myb_vec3_cross(helper, s_hat, e1);
  const double e1n = myb_vec3_norm(e1);
  for (int i = 0; i < 3; ++i) {
    e1[i] /= e1n;
  }
  myb_vec3_cross(s_hat, e1, e2);
  occulter_view views[MAX_OCCULTERS];
  size_t active = 0;
  size_t only = 0;
  for (size_t k = 0; k < count; ++k) {
    const double* o = occulters + 3 * k;
    const double d_o = myb_vec3_norm(o);
    if (d_o <= radii[k]) {
      return 0.0;
    }
    if (d_o > ds) {
      continue;
    }
    double cross[3];
    myb_vec3_cross(source, o, cross);
    const double b = myb_asin(radii[k] / d_o);
    const double c = myb_atan2(myb_vec3_norm(cross), myb_vec3_dot(source, o));
    if (c >= a + b) {
      continue;
    }
    if (c + a <= b) {
      return 0.0;
    }
    occulter_view* v = &views[active];
    v->g = (shadow_geometry){.a = a,
                             .b = b,
                             .c = c,
                             .sa = myb_sin(a),
                             .sc = myb_sin(c),
                             .u1 = limb != NULL ? limb[0] : 0.0,
                             .u2 = limb != NULL ? limb[1] : 0.0};
    for (int i = 0; i < 3; ++i) {
      v->unit[i] = o[i] / d_o;
    }
    v->azimuth = myb_atan2(myb_vec3_dot(v->unit, e2), myb_vec3_dot(v->unit, e1));
    only = k;
    ++active;
  }
  if (active == 0) {
    return 1.0;
  }
  if (active == 1) {
    return myb_shadow_fraction(source, source_radius, occulters + 3 * only, radii[only], limb);
  }
  /* where the arcs change their kind: the limbs of the occulters, and the points where two
   * of them cross */
  double cuts[MAX_CUTS];
  size_t n = 0;
  cuts[n++] = 0.0;
  cuts[n++] = a;
  for (size_t k = 0; k < active; ++k) {
    cuts[n++] = fabs(views[k].g.c - views[k].g.b);
    cuts[n++] = views[k].g.c + views[k].g.b;
    for (size_t q = k + 1; q < active; ++q) {
      const double d = myb_vec3_dot(views[k].unit, views[q].unit);
      const double span = (1.0 - d) * (1.0 + d);
      if (!(span > 1e-30)) {
        continue;
      }
      const double cbk = myb_cos(views[k].g.b);
      const double cbq = myb_cos(views[q].g.b);
      const double alpha = (cbk - d * cbq) / span;
      const double beta = (cbq - d * cbk) / span;
      const double g2 = (1.0 - alpha * alpha - beta * beta - 2.0 * alpha * beta * d) / span;
      if (!(g2 >= 0.0)) {
        continue;
      }
      double normal[3];
      myb_vec3_cross(views[k].unit, views[q].unit, normal);
      const double gamma = sqrt(g2);
      for (int sign = -1; sign <= 1; sign += 2) {
        double p[3];
        for (int i = 0; i < 3; ++i) {
          p[i] = alpha * views[k].unit[i] + beta * views[q].unit[i] + sign * gamma * normal[i];
        }
        double across[3];
        myb_vec3_cross(p, s_hat, across);
        cuts[n++] = myb_atan2(myb_vec3_norm(across), myb_vec3_dot(p, s_hat));
      }
    }
  }
  /* in order; bound: MAX_CUTS */
  for (size_t i = 1; i < n; ++i) {
    const double value = cuts[i];
    size_t j = i;
    while (j > 0 && cuts[j - 1] > value) {
      cuts[j] = cuts[j - 1];
      --j;
    }
    cuts[j] = value;
  }
  double total = 0.0;
  double blocked = 0.0;
  double from = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double to = cuts[i] < a ? cuts[i] : a;
    if (to > from) {
      integrate_many(views, active, from, to, &total, &blocked);
      from = to;
    }
  }
  const double f = 1.0 - blocked / total;
  return f < 0.0 ? 0.0 : (f > 1.0 ? 1.0 : f);
}

/* ---- Refraction in the atmosphere of the occulter ---------------------------------------------- */

typedef struct refraction_view {
  double r;      /* from the object to the centre of the occulter */
  double radius; /* of the occulter */
  double nu0;
  double scale;
  double tau0;
  double h_min;  /* the lowest tangent height from which light comes out */
} refraction_view;

double myb_refraction_bending(const myb_refraction* atmosphere, double radius, double height) {
  MYB_CHECK(atmosphere != NULL);
  if (!(atmosphere->scale_height > 0.0) || !(radius + height > 0.0)) {
    return 0.0;
  }
  return atmosphere->refractivity * myb_exp(-height / atmosphere->scale_height) *
         sqrt(2.0 * pi * (radius + height) / atmosphere->scale_height);
}

/* The impact parameter of the ray whose tangent point is at the height h (Bouguer:
 * n rho sin i is constant along the ray) and its derivative. */
static double impact_of(const refraction_view* v, double h, double* slope) {
  const double nu = v->nu0 * myb_exp(-h / v->scale);
  *slope = 1.0 + nu - (v->radius + h) * nu / v->scale;
  return (v->radius + h) * (1.0 + nu);
}

/* The tangent height of the ray of impact parameter p >= that of h_min: the impact parameter
 * grows with the height and is convex above h_min, so Newton's iteration from above comes
 * down to the root without passing it. */
static double tangent_height(const refraction_view* v, double p) {
  double h = p - v->radius;
  /* bound: 60 iterations; ends when the height stops falling */
  for (int k = 0; k < 60; ++k) {
    double slope = 0.0;
    const double excess = impact_of(v, h, &slope) - p;
    if (!(slope > 0.0)) {
      break;
    }
    const double next = h - excess / slope;
    if (!(next < h)) {
      break;
    }
    h = next > v->h_min ? next : v->h_min;
    if (h == v->h_min) {
      break;
    }
  }
  return h;
}

/* The true angle from the centre of the occulter of the light seen at the apparent angle
 * theta, its derivative with respect to theta and the transmission along the ray. */
static double true_angle(const refraction_view* v, double theta, double* slope,
                         double* transmission) {
  const double p = v->r * myb_sin(theta);
  const double h = tangent_height(v, p);
  const double shell = sqrt(2.0 * pi * (v->radius + h) / v->scale);
  const double thin = myb_exp(-h / v->scale);
  const double bend = v->nu0 * thin * shell;
  double dp_dh = 0.0;
  (void)impact_of(v, h, &dp_dh);
  const double dbend_dh = bend * (0.5 / (v->radius + h) - 1.0 / v->scale);
  *slope = dp_dh > 0.0 ? 1.0 - dbend_dh * v->r * myb_cos(theta) / dp_dh : (double)INFINITY;
  *transmission = v->tau0 > 0.0 ? myb_exp(-v->tau0 * thin * shell) : 1.0;
  return theta - bend;
}

/* Brightness of the source at the angle whose cosine from its centre is given. */
static double brightness_at(const shadow_geometry* g, double cosine) {
  if (g->u1 == 0.0 && g->u2 == 0.0) {
    return 1.0;
  }
  const double sine2 = (1.0 - cosine) * (1.0 + cosine);
  const double mu = sqrt(fmax(0.0, 1.0 - sine2 / (g->sa * g->sa)));
  const double w = 1.0 - mu;
  return 1.0 - g->u1 * w - g->u2 * w * w;
}

/* The light of the circle of radius psi about the centre of the occulter that lies on the
 * source: twice the integral of the brightness over the half-arc. */
static double source_arc(const shadow_geometry* g, double psi) {
  /* the half-angle of that circle inside the cap of the source: the roles of
   * occulted_half_angle exchanged */
  const shadow_geometry swapped = {.b = g->a, .c = g->c, .sc = g->sc};
  const double half = occulted_half_angle(&swapped, psi);
  if (!(half > 0.0)) {
    return 0.0;
  }
  if (g->u1 == 0.0 && g->u2 == 0.0) {
    return 2.0 * half;
  }
  const double cc = myb_cos(g->c);
  const double cpsi = myb_cos(psi);
  const double spsi = myb_sin(psi);
  double sum = 0.0;
  /* lambda = half (1 - u^2): the square root of the brightness at the limb becomes smooth */
  for (int k = 0; k < 24; ++k) {
    const double x = k < 12 ? myb_gauss24_x[k] : -myb_gauss24_x[k - 12];
    const double w = myb_gauss24_w[k < 12 ? k : k - 12];
    const double u = 0.5 * (x + 1.0);
    const double lambda = half * (1.0 - u * u);
    const double cosine = cc * cpsi + g->sc * spsi * myb_cos(lambda);
    sum += 0.5 * w * 2.0 * half * u * brightness_at(g, cosine);
  }
  return 2.0 * sum;
}

double myb_shadow_fraction_refracted(const double source[3], double source_radius,
                                     const double occulter[3], double occulter_radius,
                                     const myb_refraction* atmosphere, const double limb[2]) {
  MYB_CHECK(source != NULL && occulter != NULL);
  if (atmosphere == NULL || !(atmosphere->scale_height > 0.0) ||
      (atmosphere->refractivity == 0.0 && atmosphere->optical_depth == 0.0)) {
    return myb_shadow_fraction(source, source_radius, occulter, occulter_radius, limb);
  }
  const double ds = myb_vec3_norm(source);
  const double d_o = myb_vec3_norm(occulter);
  if (d_o <= occulter_radius) {
    return 0.0;
  }
  if (!(ds > source_radius) || d_o > ds || !(atmosphere->refractivity >= 0.0) ||
      !(atmosphere->optical_depth >= 0.0)) {
    return 1.0;
  }
  double cross[3];
  myb_vec3_cross(source, occulter, cross);
  shadow_geometry g;
  g.a = myb_asin(source_radius / ds);
  g.b = myb_asin(occulter_radius / d_o);
  g.c = myb_atan2(myb_vec3_norm(cross), myb_vec3_dot(source, occulter));
  g.sa = myb_sin(g.a);
  g.sc = myb_sin(g.c);
  g.u1 = limb != NULL ? limb[0] : 0.0;
  g.u2 = limb != NULL ? limb[1] : 0.0;
  refraction_view v = {.r = d_o,
                       .radius = occulter_radius,
                       .nu0 = atmosphere->refractivity,
                       .scale = atmosphere->scale_height,
                       .tau0 = atmosphere->optical_depth};
  /* light that passes more than 40 scale heights over the surface is not touched */
  const double near = g.c - g.a;
  if (!(near < 0.5 * pi) ||
      (near > 0.0 && d_o * myb_sin(near) - occulter_radius > 40.0 * v.scale)) {
    return 1.0;
  }
  /* the level under which the atmosphere bends light faster than the body curves: rays
   * below it spiral in. d p / d h falls to zero there; bound: 80 halvings */
  double slope = 0.0;
  (void)impact_of(&v, 0.0, &slope);
  if (!(slope > 0.0)) {
    double lo = 0.0;
    double hi = 60.0 * v.scale;
    for (int k = 0; k < 80; ++k) {
      const double mid = 0.5 * (lo + hi);
      (void)impact_of(&v, mid, &slope);
      if (slope > 0.0) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    v.h_min = hi;
  }
  const double p_min = impact_of(&v, v.h_min, &slope);
  if (!(d_o > p_min)) {
    return 0.0;
  }
  const double theta_limb = myb_asin(p_min / d_o);
  double transmission = 1.0;
  const double psi_limb = true_angle(&v, theta_limb, &slope, &transmission);
  const double lo = near > psi_limb ? near : psi_limb;
  const double hi = g.c + g.a;
  if (!(hi > lo)) {
    return 0.0;
  }
  double total = 0.0;
  double unused = 0.0;
  integrate(&g, 0.0, g.a, &total, &unused);
  /* over the true angle psi from the centre of the occulter, in eight parts, each with the
   * substitution that smooths the ends of the disk */
  double seen = 0.0;
  for (int part = 0; part < 8; ++part) {
    const double from = lo + (hi - lo) * part / 8.0;
    const double to = lo + (hi - lo) * (part + 1) / 8.0;
    const double half = 0.5 * (to - from);
    for (int k = 0; k < 24; ++k) {
      const double x = k < 12 ? myb_gauss24_x[k] : -myb_gauss24_x[k - 12];
      const double w = myb_gauss24_w[k < 12 ? k : k - 12];
      const double t = 0.5 * pi * (x + 1.0);
      const double psi = from + half * (1.0 - myb_cos(t));
      const double jac = half * myb_sin(t) * 0.5 * pi;
      /* the apparent angle of this light: the true angle grows with the apparent one and is
       * concave, so Newton's iteration from below climbs to it; bound: 60 iterations */
      double theta = psi > theta_limb ? psi : theta_limb;
      double weight = 0.0;
      for (int it = 0; it < 60; ++it) {
        const double value = true_angle(&v, theta, &slope, &transmission);
        weight = transmission * myb_sin(theta) / slope;
        if (!(slope > 0.0) || !isfinite(slope)) {
          weight = 0.0;
          break;
        }
        const double next = theta + (psi - value) / slope;
        if (!(next > theta) || next >= 0.5 * pi) {
          break;
        }
        theta = next;
      }
      seen += w * jac * weight * source_arc(&g, psi);
    }
  }
  const double f = seen / total;
  return f < 0.0 ? 0.0 : (f > 1.0 ? 1.0 : f);
}
