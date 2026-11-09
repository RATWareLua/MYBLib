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
