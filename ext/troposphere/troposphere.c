/* SPDX-License-Identifier: BSD-3-Clause */
/* Tropospheric delay models of the IERS Conventions (2010), chapter 9 (myb_troposphere.h),
 * written from the published equations; the coefficient tables of GMF and GPT are data
 * (troposphere_tables.inc). */
#include <math.h>

#include "myb/cmath.h"
#include "myb_troposphere.h"

#include "troposphere_tables.inc"

#define TWO_PI 6.283185307179586476925287

/* Day of year counted from 28 January 1980 (MJD 44239 = 1 January 1980), as A. E. Niell (1996)
 * and the models built on his seasonal phase. */
static double niell_doy(double mjd) {
  return mjd - 44239.0 + 1.0 - 28.0;
}

/* Marini's continued fraction normalised to 1 at the zenith. */
static double marini(double sine, double a, double b, double c) {
  const double top = 1.0 + a / (1.0 + b / (1.0 + c));
  return top / (sine + a / (sine + b / (sine + c)));
}

/* Niell's height correction of the hydrostatic mapping function (height in m). */
static double niell_height(double sine, double height) {
  const double a_ht = 2.53e-5;
  const double b_ht = 5.49e-3;
  const double c_ht = 1.14e-3;
  return (1.0 / sine - marini(sine, a_ht, b_ht, c_ht)) * (height / 1000.0);
}

/* The c coefficient of the hydrostatic VMF1 and GMF (seasonal, by hemisphere). */
static double hydrostatic_c(double mjd, double latitude) {
  const double c0 = 0.062;
  const int south = latitude < 0.0;
  const double phase = south ? 3.14159265358979323846 : 0.0;
  const double c11 = south ? 0.007 : 0.005;
  const double c10 = south ? 0.002 : 0.001;
  return c0 + ((myb_cos(niell_doy(mjd) / 365.25 * TWO_PI + phase) + 1.0) * c11 / 2.0 + c10) *
                  (1.0 - myb_cos(latitude));
}

static int finite_all(double a, double b, double c, double d) {
  return isfinite(a) && isfinite(b) && isfinite(c) && isfinite(d);
}

myb_status myb_tropo_zhd_saastamoinen(double pressure, double latitude, double height,
                                      double* zhd) {
  if (zhd == NULL || !finite_all(pressure, latitude, height, 0.0) || !(pressure >= 0.0)) {
    return MYB_E_DOMAIN;
  }
  *zhd = 0.0022768 * pressure / (1.0 - 0.00266 * myb_cos(2.0 * latitude) - 0.28e-6 * height);
  return MYB_OK;
}

myb_status myb_tropo_zwd_saastamoinen(double temperature, double vapour_pressure, double* zwd) {
  if (zwd == NULL || !(temperature > 0.0) || !isfinite(temperature) ||
      !(vapour_pressure >= 0.0) || !isfinite(vapour_pressure)) {
    return MYB_E_DOMAIN;
  }
  *zwd = 0.002277 * (1255.0 / temperature + 0.05) * vapour_pressure;
  return MYB_OK;
}

myb_status myb_tropo_zenith_optical(double latitude, double height, double pressure,
                                    double vapour_pressure, double wavelength, double co2,
                                    double* zhd, double* zwd) {
  if (zhd == NULL || zwd == NULL || !finite_all(latitude, height, pressure, vapour_pressure) ||
      !(wavelength >= 0.25 && wavelength <= 2.0) || !(co2 >= 0.0 && co2 < 1e6) ||
      !(pressure >= 0.0) || !(vapour_pressure >= 0.0)) {
    return MYB_E_DOMAIN;
  }
  /* Mendes & Pavlis (2004), eq. 20 (dispersion of the hydrostatic refractivity, Ciddor 1996),
   * 24, 26, 32, 38; sigma is the wave number in 1/micrometre */
  const double k0 = 238.0185;
  const double k1 = 19990.975;
  const double k2 = 57.362;
  const double k3 = 579.55174;
  const double w0 = 295.235;
  const double w1 = 2.6422;
  const double w2 = -0.032380;
  const double w3 = 0.004028;
  const double sigma = 1.0 / wavelength;
  const double s2 = sigma * sigma;
  const double f = 1.0 - 0.00266 * myb_cos(2.0 * latitude) - 0.00028e-3 * height;
  const double corr = 1.0 + 0.534e-6 * (co2 - 450.0);
  const double fh = 0.01 * corr *
                    ((k1 * (k0 + s2)) / ((k0 - s2) * (k0 - s2)) +
                     k3 * (k2 + s2) / ((k2 - s2) * (k2 - s2)));
  const double fnh = 0.003101 * (w0 + 3.0 * w1 * s2 + 5.0 * w2 * s2 * s2 + 7.0 * w3 * s2 * s2 * s2);
  *zhd = 2.416579e-3 * fh * pressure / f;
  *zwd = 1e-4 * (5.316 * fnh - 3.759 * fh) * vapour_pressure / f;
  return MYB_OK;
}

myb_status myb_tropo_fcula(double elevation, double latitude, double height, double temperature,
                           double* mapping) {
  if (mapping == NULL || !finite_all(elevation, latitude, height, temperature) ||
      !(elevation > 0.0 && elevation <= 1.5707963267948966) || !(temperature > 0.0)) {
    return MYB_E_DOMAIN;
  }
  /* Mendes et al. (2002), IERS Conventions (2010) table 9.1: a_i = a_i0 + a_i1 t + a_i2 cos phi
   * + a_i3 H with t in degrees Celsius and H in metres */
  const double t = temperature - 273.15;
  const double cphi = myb_cos(latitude);
  const double a1 = 0.121008e-2 + 0.17295e-5 * t + 0.3191e-4 * cphi - 0.18478e-7 * height;
  const double a2 = 0.304965e-2 + 0.2346e-5 * t - 0.1035e-3 * cphi - 0.1856e-7 * height;
  const double a3 = 0.68777e-1 + 0.1972e-4 * t - 0.3458e-2 * cphi + 0.1060e-6 * height;
  *mapping = marini(myb_sin(elevation), a1, a2, a3);
  return MYB_OK;
}

static int zenith_ok(double z) {
  return z >= 0.0 && z < 1.5707963267948966;
}

myb_status myb_tropo_vmf1(double ah, double aw, double mjd, double latitude, double height,
                          double zenith_distance, double* mh, double* mw) {
  if (mh == NULL || mw == NULL || !finite_all(ah, aw, mjd, latitude) || !isfinite(height) ||
      !zenith_ok(zenith_distance)) {
    return MYB_E_DOMAIN;
  }
  const double sine = myb_cos(zenith_distance);
  *mh = marini(sine, ah, 0.0029, hydrostatic_c(mjd, latitude)) + niell_height(sine, height);
  *mw = marini(sine, aw, 0.00146, 0.04391);
  return MYB_OK;
}

/* Unnormalised Legendre functions times cos(m lon), sin(m lon) to degree and order 9 at the
 * point of the unit sphere (x, y, z), in the order n = 0..9, m = 0..n (recurrences in the
 * Cartesian form, without the Condon-Shortley phase). */
static void harmonics9(double latitude, double longitude, double v[55], double w[55]) {
  enum { N = 9 };
  double vv[N + 1][N + 1] = {{0.0}};
  double ww[N + 1][N + 1] = {{0.0}};
  const double cl = myb_cos(latitude);
  const double x = cl * myb_cos(longitude);
  const double y = cl * myb_sin(longitude);
  const double z = myb_sin(latitude);
  vv[0][0] = 1.0;
  vv[1][0] = z;
  for (int n = 2; n <= N; ++n) {
    vv[n][0] = ((2 * n - 1) * z * vv[n - 1][0] - (n - 1) * vv[n - 2][0]) / n;
  }
  for (int m = 1; m <= N; ++m) {
    vv[m][m] = (2 * m - 1) * (x * vv[m - 1][m - 1] - y * ww[m - 1][m - 1]);
    ww[m][m] = (2 * m - 1) * (x * ww[m - 1][m - 1] + y * vv[m - 1][m - 1]);
    if (m < N) {
      vv[m + 1][m] = (2 * m + 1) * z * vv[m][m];
      ww[m + 1][m] = (2 * m + 1) * z * ww[m][m];
    }
    for (int n = m + 2; n <= N; ++n) {
      vv[n][m] = ((2 * n - 1) * z * vv[n - 1][m] - (n + m - 1) * vv[n - 2][m]) / (n - m);
      ww[n][m] = ((2 * n - 1) * z * ww[n - 1][m] - (n + m - 1) * ww[n - 2][m]) / (n - m);
    }
  }
  int i = 0;
  for (int n = 0; n <= N; ++n) {
    for (int m = 0; m <= n; ++m) {
      v[i] = vv[n][m];
      w[i] = ww[n][m];
      ++i;
    }
  }
}

static double expand(const double* a, const double* b, const double* v, const double* w) {
  double s = 0.0;
  for (int i = 0; i < 55; ++i) {
    s += a[i] * v[i] + b[i] * w[i];
  }
  return s;
}

myb_status myb_tropo_gmf(double mjd, double latitude, double longitude, double height,
                         double zenith_distance, double* mh, double* mw) {
  if (mh == NULL || mw == NULL || !finite_all(mjd, latitude, longitude, height) ||
      !zenith_ok(zenith_distance)) {
    return MYB_E_DOMAIN;
  }
  double v[55];
  double w[55];
  harmonics9(latitude, longitude, v, w);
  const double season = myb_cos(niell_doy(mjd) / 365.25 * TWO_PI);
  const double ah = (expand(gmf_ah_mean, gmf_bh_mean, v, w) +
                     expand(gmf_ah_amp, gmf_bh_amp, v, w) * season) *
                    1e-5;
  const double aw = (expand(gmf_aw_mean, gmf_bw_mean, v, w) +
                     expand(gmf_aw_amp, gmf_bw_amp, v, w) * season) *
                    1e-5;
  const double sine = myb_cos(zenith_distance);
  *mh = marini(sine, ah, 0.0029, hydrostatic_c(mjd, latitude)) + niell_height(sine, height);
  *mw = marini(sine, aw, 0.00146, 0.04391);
  return MYB_OK;
}

myb_status myb_tropo_gpt(double mjd, double latitude, double longitude, double height,
                         double* pressure, double* temperature, double* undulation) {
  if (pressure == NULL || temperature == NULL || undulation == NULL ||
      !finite_all(mjd, latitude, longitude, height) || !(height > -1000.0 && height < 40000.0)) {
    return MYB_E_DOMAIN;
  }
  double v[55];
  double w[55];
  harmonics9(latitude, longitude, v, w);
  const double season = myb_cos(niell_doy(mjd) / 365.25 * TWO_PI);
  const double und = expand(gpt_a_geoid, gpt_b_geoid, v, w);
  const double orthometric = height - und;
  const double p0 = expand(gpt_ap_mean, gpt_bp_mean, v, w) +
                    expand(gpt_ap_amp, gpt_bp_amp, v, w) * season;
  const double t0 = expand(gpt_at_mean, gpt_bt_mean, v, w) +
                    expand(gpt_at_amp, gpt_bt_amp, v, w) * season;
  /* height corrections of the model: barometric (1 - 2.26e-5 H)^5.225 and -6.5 K/km */
  *pressure = p0 * myb_exp(5.225 * myb_log(1.0 - 0.0000226 * orthometric));
  *temperature = t0 - 0.0065 * orthometric;
  *undulation = und;
  return MYB_OK;
}

myb_status myb_tropo_gradient_delay(double gn, double ge, double elevation, double azimuth,
                                    double* delay) {
  if (delay == NULL || !finite_all(gn, ge, elevation, azimuth) ||
      !(elevation > 0.0 && elevation <= 1.5707963267948966)) {
    return MYB_E_DOMAIN;
  }
  const double se = myb_sin(elevation);
  const double mg = 1.0 / (se * se / myb_cos(elevation) + 0.0032);
  *delay = mg * (gn * myb_cos(azimuth) + ge * myb_sin(azimuth));
  return MYB_OK;
}
