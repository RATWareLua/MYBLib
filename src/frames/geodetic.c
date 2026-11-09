/* SPDX-License-Identifier: BSD-3-Clause */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "myb/cmath.h"
#include "myb/frames.h"
#include "num/vec3.h"

static const double half_pi = 0x1.921fb54442d18p+0;

static int valid(const myb_ellipsoid* e) {
  return e != NULL && e->a > 0.0 && e->a < 1e300 && e->f >= 0.0 && e->f < 1.0;
}

myb_status myb_geodetic_to_cartesian(const myb_ellipsoid* ellipsoid, double longitude,
                                     double latitude, double height, double xyz[3]) {
  MYB_CHECK(xyz != NULL);
  if (!valid(ellipsoid) || !isfinite(longitude) || !(fabs(latitude) <= half_pi + 1e-15) ||
      !isfinite(height)) {
    return MYB_E_DOMAIN;
  }
  const double a = ellipsoid->a;
  const double f = ellipsoid->f;
  const double e2 = f * (2.0 - f);
  const double sp = myb_sin(latitude);
  const double cp = myb_cos(latitude);
  const double n = a / sqrt(1.0 - e2 * sp * sp); /* prime vertical radius */
  const double r = (n + height) * cp;
  xyz[0] = r * myb_cos(longitude);
  xyz[1] = r * myb_sin(longitude);
  xyz[2] = (n * (1.0 - f) * (1.0 - f) + height) * sp;
  return MYB_OK;
}

myb_status myb_geodetic_from_cartesian(const myb_ellipsoid* ellipsoid, const double xyz[3],
                                       double* longitude, double* latitude, double* height) {
  MYB_CHECK(xyz != NULL && longitude != NULL && latitude != NULL && height != NULL);
  if (!valid(ellipsoid) || !myb_vec3_finite(xyz)) {
    return MYB_E_DOMAIN;
  }
  const double a = ellipsoid->a;
  const double f = ellipsoid->f;
  const double e2 = f * (2.0 - f);
  const double e4 = e2 * e2;
  const double x = xyz[0];
  const double y = xyz[1];
  const double z = xyz[2];
  const double rho = myb_hypot(x, y);
  *longitude = (x == 0.0 && y == 0.0) ? 0.0 : myb_atan2(y, x);
  if (rho == 0.0 && z == 0.0) {
    *latitude = half_pi;
    *height = -a * (1.0 - f);
    return MYB_OK;
  }
  /* Vermeille (2011), section 4: p, q, r, the evolute test 8 r^3 + e^4 p q, then u. */
  const double p = (rho / a) * (rho / a);
  const double q = (1.0 - e2) * (z / a) * (z / a);
  const double r = (p + q - e4) / 6.0;
  const double epq = e4 * p * q;
  const double evolute = 8.0 * r * r * r + epq;
  if (z == 0.0 && p < e4) {
    /* On the equatorial plane inside the evolute the normal from (rho, 0) meets the
     * ellipsoid at +-phi with rho = a e^2 cos(phi) / sqrt(1 - e^2 sin^2 phi); the northern
     * solution: */
    const double c = rho * sqrt(1.0 - e2) / sqrt(e2 * (a * a * e2 - rho * rho));
    const double phi = myb_acos(c > 1.0 ? 1.0 : c);
    const double n = a / sqrt(1.0 - e2 * myb_sin(phi) * myb_sin(phi));
    *latitude = phi;
    *height = (c > 0.0 ? rho / c : 0.0) - n;
    if (c == 0.0) {
      *height = -a * (1.0 - f);
    }
    return MYB_OK;
  }
  double u;
  if (evolute > 0.0) {
    const double t = myb_cbrt(sqrt(evolute) + sqrt(epq));
    const double t1 = t * t;
    u = r + 0.5 * t1 + 2.0 * r * r / t1;
  } else {
    const double angle = myb_atan2(sqrt(epq), sqrt(-evolute) + sqrt(-8.0 * r * r * r));
    u = 2.0 * r * myb_cos(angle * (2.0 / 3.0));
  }
  const double v = sqrt(u * u + e4 * q);
  const double w = e2 * (u + v - q) / (2.0 * v);
  const double k = (u + v) / (sqrt(w * w + u + v) + w);
  const double d = k * rho / (k + e2);
  const double s = myb_hypot(d, z);
  *latitude = 2.0 * myb_atan2(z, d + s);
  *height = (k + e2 - 1.0) / k * s;
  return MYB_OK;
}

void myb_topocentric_enu(double longitude, double latitude, myb_rotation* rotation) {
  MYB_CHECK(rotation != NULL);
  const double sl = myb_sin(longitude);
  const double cl = myb_cos(longitude);
  const double sp = myb_sin(latitude);
  const double cp = myb_cos(latitude);
  const double m[9] = {-sl, cl, 0.0, -sp * cl, -sp * sl, cp, cp * cl, cp * sl, sp};
  memcpy(rotation->m, m, sizeof m);
  memset(rotation->dm, 0, sizeof rotation->dm);
}
