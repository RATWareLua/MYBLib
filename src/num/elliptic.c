/* SPDX-License-Identifier: BSD-3-Clause */
/* Elliptic integrals (Carlson's symmetric forms and Legendre's) and Jacobi's functions.
 *
 * Carlson's algorithms (B. C. Carlson, Numer. Algorithms 10, 13, 1995; DLMF 19.36): duplication
 * until the scaled spread 4^-m Q of the arguments falls below their mean A_m, with
 * Q = (3r)^(-1/6) max |A_0 - x_i| (R_F), (r/4)^(-1/6) ... (R_D, R_J), (3r)^(-1/8) (R_C) and
 * r = 2^-64, then the Taylor series of DLMF 19.36.1-2 to fifth order, whose truncation is then
 * below r. The duplication count is log4 of the spread of the arguments: at most ~1100 over the
 * range of binary64, checked.
 *
 * Jacobi's amplitude by the arithmetic-geometric mean (DLMF 22.20(ii)) after the reduction of
 * u by multiples of 2K (am(u + 2jK) = am(u) + j pi); sn, cn from the amplitude, dn as
 * sqrt((1 - m) + m cn^2) (no cancellation). */
#include <math.h>

#include "core/check.h"
#include "myb/cmath.h"
#include "myb/num.h"

static const double pi_hi = 0x1.921fb54442d18p+1;
static const double pi_lo = 0x1.1a62633145c07p-53;

enum { MAX_DUPLICATIONS = 1200 };

static double max4(double a, double b, double c, double d) {
  return fmax(fmax(a, b), fmax(c, d));
}

double myb_carlson_rf(double x, double y, double z) {
  if (!(x >= 0.0 && y >= 0.0 && z >= 0.0) || (x == 0.0) + (y == 0.0) + (z == 0.0) > 1 ||
      !isfinite(x + y + z)) {
    return (double)NAN;
  }
  const double a0 = (x + y + z) / 3.0;
  /* (3 r)^(-1/6), r = 2^-64 */
  const double q = 1.0e3 * max4(fabs(a0 - x), fabs(a0 - y), fabs(a0 - z), 0.0);
  double a = a0;
  double scale = 1.0; /* 4^-m */
  double xm = x;
  double ym = y;
  double zm = z;
  int done = 0;
  for (int m = 0; m < MAX_DUPLICATIONS; ++m) {
    if (q * scale < fabs(a)) {
      done = 1;
      break;
    }
    const double sx = sqrt(xm);
    const double sy = sqrt(ym);
    const double sz = sqrt(zm);
    const double lambda = sx * sy + sx * sz + sy * sz;
    xm = (xm + lambda) * 0.25;
    ym = (ym + lambda) * 0.25;
    zm = (zm + lambda) * 0.25;
    a = (a + lambda) * 0.25;
    scale *= 0.25;
  }
  MYB_CHECK(done);
  const double X = (a0 - x) * scale / a;
  const double Y = (a0 - y) * scale / a;
  const double Z = -(X + Y);
  const double e2 = X * Y - Z * Z;
  const double e3 = X * Y * Z;
  /* DLMF 19.36.1 */
  const double s = 1.0 - e2 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0 - 3.0 * e2 * e3 / 44.0 -
                   5.0 * e2 * e2 * e2 / 208.0 + 3.0 * e3 * e3 / 104.0 + e2 * e2 * e3 / 16.0;
  return s / sqrt(a);
}

double myb_carlson_rc(double x, double y) {
  if (!(x >= 0.0 && y > 0.0) || !isfinite(x + y)) {
    return (double)NAN;
  }
  const double a0 = (x + 2.0 * y) / 3.0;
  /* (3 r)^(-1/8), r = 2^-64 */
  const double q = 1.8e2 * fabs(a0 - x);
  double a = a0;
  double scale = 1.0;
  double xm = x;
  double ym = y;
  int done = 0;
  for (int m = 0; m < MAX_DUPLICATIONS; ++m) {
    if (q * scale < fabs(a)) {
      done = 1;
      break;
    }
    const double lambda = 2.0 * sqrt(xm) * sqrt(ym) + ym;
    xm = (xm + lambda) * 0.25;
    ym = (ym + lambda) * 0.25;
    a = (a + lambda) * 0.25;
    scale *= 0.25;
  }
  MYB_CHECK(done);
  const double s = (y - a0) * scale / a;
  const double t = 1.0 + s * s * (3.0 / 10.0 + s * (1.0 / 7.0 + s * (3.0 / 8.0 + s * (9.0 / 22.0 +
                    s * (159.0 / 208.0 + s * (9.0 / 8.0))))));
  return t / sqrt(a);
}

/* R_J with p > 0; R_D(x, y, z) = R_J(x, y, z, z). */
double myb_carlson_rj(double x, double y, double z, double p) {
  if (!(x >= 0.0 && y >= 0.0 && z >= 0.0 && p > 0.0) ||
      (x == 0.0) + (y == 0.0) + (z == 0.0) > 1 || !isfinite(x + y + z + p)) {
    return (double)NAN;
  }
  const double a0 = (x + y + z + 2.0 * p) / 5.0;
  /* (r / 4)^(-1/6), r = 2^-64 */
  const double q = 2.1e3 * max4(fabs(a0 - x), fabs(a0 - y), fabs(a0 - z), fabs(a0 - p));
  double a = a0;
  double scale = 1.0;
  double xm = x;
  double ym = y;
  double zm = z;
  double pm = p;
  double sum = 0.0;
  double comp = 0.0;
  int done = 0;
  for (int m = 0; m < MAX_DUPLICATIONS; ++m) {
    if (q * scale < fabs(a)) {
      done = 1;
      break;
    }
    const double sx = sqrt(xm);
    const double sy = sqrt(ym);
    const double sz = sqrt(zm);
    const double lambda = sx * sy + sx * sz + sy * sz;
    /* duplication: R_J(x, y, z, p) = 2 R_J(x + l, y + l, z + l, p + l) + 3 R_C(alpha, beta)
     * (Carlson 1988), both arguments positive, so no cancellation in 1 + e of the 1995 form */
    const double t = pm * (sx + sy + sz) + sx * sy * sz;
    const double u = pm + lambda;
    const double term = scale * myb_carlson_rc(t * t, pm * u * u);
    double err;
    sum = myb_two_sum(sum, term, &err);
    comp += err;
    xm = (xm + lambda) * 0.25;
    ym = (ym + lambda) * 0.25;
    zm = (zm + lambda) * 0.25;
    pm = (pm + lambda) * 0.25;
    a = (a + lambda) * 0.25;
    scale *= 0.25;
  }
  MYB_CHECK(done);
  const double X = (a0 - x) * scale / a;
  const double Y = (a0 - y) * scale / a;
  const double Z = (a0 - z) * scale / a;
  const double P = -(X + Y + Z) / 2.0;
  const double e2 = X * Y + X * Z + Y * Z - 3.0 * P * P;
  const double e3 = X * Y * Z + 2.0 * e2 * P + 4.0 * P * P * P;
  const double e4 = (2.0 * X * Y * Z + e2 * P + 3.0 * P * P * P) * P;
  const double e5 = X * Y * Z * P * P;
  /* DLMF 19.36.2 with the terms of Carlson (1995) */
  const double s = 1.0 - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2 * e2 / 88.0 - 3.0 * e4 / 22.0 -
                   9.0 * e2 * e3 / 52.0 + 3.0 * e5 / 26.0;
  return scale * s / (a * sqrt(a)) + 3.0 * (sum + comp);
}

double myb_carlson_rd(double x, double y, double z) {
  if (!(z > 0.0)) {
    return (double)NAN;
  }
  return myb_carlson_rj(x, y, z, z);
}

/* ---- Legendre's forms --------------------------------------------------------------------- */

double myb_elliptic_k(double m) {
  if (!(m < 1.0) || isnan(m)) {
    return m == 1.0 ? (double)INFINITY : (double)NAN;
  }
  return myb_carlson_rf(0.0, 1.0 - m, 1.0);
}

/* E(m) = 2 R_G(0, 1, 1 - m) = y R_F(0, 1, y) + m y R_D(0, 1, y) / 3, y = 1 - m: both terms
 * positive for 0 <= m < 1 (R_G with the smallest argument last; DLMF 19.21.10); for m < 0 the
 * usual R_F - m R_D / 3 has positive terms. */
double myb_elliptic_e(double m) {
  if (!(m <= 1.0)) {
    return (double)NAN;
  }
  if (m == 1.0) {
    return 1.0;
  }
  const double y = 1.0 - m;
  if (m < 0.0) {
    return myb_carlson_rf(0.0, y, 1.0) - m / 3.0 * myb_carlson_rd(0.0, y, 1.0);
  }
  return y * myb_carlson_rf(0.0, 1.0, y) + m * y / 3.0 * myb_carlson_rd(0.0, 1.0, y);
}

/* Pi for 0 <= s, phi in [0, pi / 2] given by s = sin phi, c = cos phi, and n < 1 or
 * n s^2 < 1. Delta^2 = 1 - m s^2 = c^2 + (1 - m) s^2 and 1 - n s^2 = c^2 + (1 - n) s^2 without
 * cancellation. For n < 0 <= m the terms R_F and n R_J cancel (Pi -> 0 as n -> -inf): then
 * A&S 17.7.15 maps to N = (m - n) / (1 - n) in (m, 1), all terms positive. */
static double pi_part(double n, double s, double c, double m) {
  const double s2 = s * s;
  const double d2 = c * c + (1.0 - m) * s2;
  if (n < 0.0 && m >= 0.0) {
    const double big = (m - n) / (1.0 - n);
    const double big_m1 = (1.0 - m) / (1.0 - n); /* 1 - N */
    const double pn = s * myb_carlson_rf(c * c, d2, 1.0) +
                      big / 3.0 * s * s2 * myb_carlson_rj(c * c, d2, 1.0, c * c + big_m1 * s2);
    const double f = s * myb_carlson_rf(c * c, d2, 1.0);
    const double p2 = sqrt(-n * big);
    const double w = sqrt(-n / ((m - n) * (1.0 - n)));
    return pn * ((-n) * (1.0 - m) / ((1.0 - n) * (m - n))) + f * (m / (m - n)) +
           w * myb_atan(p2 * s * c / sqrt(d2));
  }
  const double p = n <= 1.0 ? c * c + (1.0 - n) * s2 : 1.0 - n * s2;
  return s * myb_carlson_rf(c * c, d2, 1.0) + n / 3.0 * s * s2 * myb_carlson_rj(c * c, d2, 1.0, p);
}

double myb_elliptic_pi(double n, double m) {
  if (!(m < 1.0 && n < 1.0)) {
    return (double)NAN;
  }
  return pi_part(n, 1.0, 0.0, m);
}

/* phi = j pi + r, |r| <= pi / 2 (r with one rounding). */
static double reduce_phi(double phi, double* j) {
  const double q = floor(phi / pi_hi + 0.5);
  const double r = myb_fma(-q, pi_lo, myb_fma(-q, pi_hi, phi));
  *j = q;
  return r;
}

double myb_elliptic_f_inc(double phi, double m) {
  if (!(m <= 1.0) || !isfinite(phi)) {
    return (double)NAN;
  }
  double j;
  const double r = reduce_phi(phi, &j);
  const double s = myb_sin(r);
  const double c = myb_cos(r);
  const double part = s * myb_carlson_rf(c * c, c * c + (1.0 - m) * s * s, 1.0);
  if (j == 0.0) {
    return part;
  }
  return m < 1.0 ? 2.0 * j * myb_elliptic_k(m) + part : copysign((double)INFINITY, phi);
}

/* E(phi, m) = (1 - m) s R_F(c^2, d^2, 1) + m (1 - m) s^3 R_D(c^2, 1, d^2) / 3 + m s c / d
 * (DLMF 19.25.10 after homogeneity), positive terms for 0 <= m <= 1. */
double myb_elliptic_e_inc(double phi, double m) {
  if (!(m <= 1.0) || !isfinite(phi)) {
    return (double)NAN;
  }
  double j;
  const double r = reduce_phi(phi, &j);
  const double s = myb_sin(r);
  const double c = myb_cos(r);
  const double s2 = s * s;
  const double d2 = c * c + (1.0 - m) * s2;
  double part;
  if (m < 0.0) {
    part = s * myb_carlson_rf(c * c, d2, 1.0) - m / 3.0 * s * s2 * myb_carlson_rd(c * c, d2, 1.0);
  } else if (d2 == 0.0) {
    part = s; /* m = 1, phi = +-pi/2 */
  } else {
    const double y = 1.0 - m;
    part = y * s * myb_carlson_rf(c * c, d2, 1.0) +
           m * y / 3.0 * s * s2 * myb_carlson_rd(c * c, 1.0, d2) + m * s * fabs(c) / sqrt(d2);
  }
  return j == 0.0 ? part : 2.0 * j * myb_elliptic_e(m) + part;
}

double myb_elliptic_pi_inc(double n, double phi, double m) {
  if (!(m <= 1.0) || !isfinite(phi) || isnan(n)) {
    return (double)NAN;
  }
  double j;
  const double r = reduce_phi(phi, &j);
  const double s = myb_sin(r);
  const double c = myb_cos(r);
  if (!(n * s * s < 1.0) || (j != 0.0 && !(n < 1.0 && m < 1.0))) {
    return (double)NAN;
  }
  const double part = copysign(pi_part(n, fabs(s), c, m), s);
  return j == 0.0 ? part : 2.0 * j * myb_elliptic_pi(n, m) + part;
}

/* ---- Jacobi's functions ------------------------------------------------------------------- */

enum { MAX_AGM = 64 };

/* am(u) for |u| <= K by the arithmetic-geometric mean. */
static double amplitude(double u, double m) {
  if (m == 0.0) {
    return u;
  }
  double a[MAX_AGM];
  double c[MAX_AGM];
  a[0] = 1.0;
  double b = sqrt(1.0 - m);
  c[0] = sqrt(m);
  int n = 0;
  while (n + 1 < MAX_AGM && fabs(c[n]) > 0x1p-51 * a[n]) { /* bound: MAX_AGM; quadratic */
    a[n + 1] = (a[n] + b) * 0.5;
    c[n + 1] = (a[n] - b) * 0.5;
    b = sqrt(a[n] * b);
    ++n;
  }
  MYB_CHECK(n + 1 < MAX_AGM);
  double phi = ldexp(a[n] * u, n);
  for (int k = n; k >= 1; --k) {
    phi = 0.5 * (phi + myb_asin(c[k] / a[k] * myb_sin(phi)));
  }
  return phi;
}

void myb_jacobi(double u, double m, double* am, double* sn, double* cn, double* dn) {
  if (!(m >= 0.0 && m <= 1.0) || !isfinite(u)) {
    const double nan = (double)NAN;
    if (am != NULL) {
      *am = nan;
    }
    if (sn != NULL) {
      *sn = nan;
    }
    if (cn != NULL) {
      *cn = nan;
    }
    if (dn != NULL) {
      *dn = nan;
    }
    return;
  }
  double phi;
  double s;
  double c;
  if (m == 1.0) {
    /* am = gd(u): sn = tanh u, cn = dn = sech u */
    phi = myb_atan(myb_sinh(u));
    s = myb_tanh(u);
    c = 1.0 / myb_cosh(u);
  } else {
    const double k = myb_elliptic_k(m);
    const double j = floor(u / (2.0 * k) + 0.5);
    const double r = myb_fma(-2.0 * j, k, u);
    const double phi_r = amplitude(r, m);
    phi = myb_fma(j, pi_hi, phi_r) + j * pi_lo;
    const double half = 0.5 * j;
    const double sign = half == floor(half) ? 1.0 : -1.0;
    s = sign * myb_sin(phi_r);
    c = sign * myb_cos(phi_r);
  }
  if (am != NULL) {
    *am = phi;
  }
  if (sn != NULL) {
    *sn = s;
  }
  if (cn != NULL) {
    *cn = c;
  }
  if (dn != NULL) {
    *dn = sqrt((1.0 - m) + m * c * c);
  }
}
