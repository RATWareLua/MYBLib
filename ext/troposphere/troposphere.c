/* SPDX-License-Identifier: BSD-3-Clause */
/* Tropospheric delay models of the IERS Conventions (2010), chapter 9, and VMF3 with GPT3
 * (myb_troposphere.h), written from the published equations; the coefficient tables of GMF,
 * GPT and VMF3 are data (troposphere_tables.inc, vmf3_tables.inc). */
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/check.h"
#include "myb/cmath.h"
#include "myb/io.h"
#include "myb/time.h"
#include "myb_troposphere.h"

#include "troposphere_tables.inc"
#include "vmf3_tables.inc"

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

myb_status myb_tropo_gradient_delay_parts(double gn_h, double ge_h, double gn_w, double ge_w,
                                          double elevation, double azimuth, double* delay) {
  if (delay == NULL || !finite_all(gn_h, ge_h, gn_w, ge_w) || !isfinite(elevation) ||
      !isfinite(azimuth) || !(elevation > 0.0 && elevation <= 1.5707963267948966)) {
    return MYB_E_DOMAIN;
  }
  const double se = myb_sin(elevation);
  const double tangent = se * se / myb_cos(elevation); /* sin e tan e */
  const double ca = myb_cos(azimuth);
  const double sa = myb_sin(azimuth);
  *delay = (gn_h * ca + ge_h * sa) / (tangent + 0.0031) +
           (gn_w * ca + ge_w * sa) / (tangent + 0.0007);
  return MYB_OK;
}

myb_status myb_tropo_zwd_askne(double vapour_pressure, double mean_temperature, double lambda,
                               double* zwd) {
  if (zwd == NULL || !isfinite(vapour_pressure) || !(mean_temperature > 0.0) ||
      !isfinite(mean_temperature) || !(lambda > -1.0) || !isfinite(lambda)) {
    return MYB_E_DOMAIN;
  }
  const double k1 = 77.604; /* K/hPa */
  const double k2 = 64.79;
  const double k3 = 377600.0; /* K^2/hPa */
  const double k2p = k2 - k1 * 18.0152 / 28.9644;
  const double rd = 8.3143 / 28.965e-3; /* of dry air, J/(kg K) */
  *zwd = 1.0e-6 * (k2p + k3 / mean_temperature) * rd / (lambda + 1.0) / 9.80665 * vapour_pressure;
  return MYB_OK;
}

/* ---- VMF3 and GPT3 --------------------------------------------------------------------------- */

/* The day of the year, 1 at the midnight that begins it, with the fraction of the day: the
 * argument of the seasons of VMF3 and GPT3. */
static myb_status year_day(double mjd, double* day) {
  if (!(fabs(mjd) < 1e7)) {
    return MYB_E_DOMAIN;
  }
  const double whole = floor(mjd);
  const myb_time midnight = {.s = ((int64_t)whole - 51544) * 86400 - 43200, .f = 0.0};
  myb_calendar date;
  myb_time_to_calendar(midnight, &date);
  const myb_calendar first = {.year = date.year, .month = 1, .day = 1};
  myb_time start;
  const myb_status status = myb_calendar_to_time(&first, &start);
  if (status != MYB_OK) {
    return status;
  }
  *day = (double)((midnight.s - start.s) / 86400) + 1.0 + (mjd - whole);
  return MYB_OK;
}

/* cos and sin of one and two turns a year at the day: what the mean and the amplitudes A1,
 * B1, A2, B2 of a seasonal quantity are multiplied by. */
static void seasons(double day, int constant, double out[5]) {
  out[0] = 1.0;
  for (int i = 1; i < 5; ++i) {
    out[i] = 0.0;
  }
  if (!constant) {
    const double year = day / 365.25 * TWO_PI;
    out[1] = myb_cos(year);
    out[2] = myb_sin(year);
    out[3] = myb_cos(2.0 * year);
    out[4] = myb_sin(2.0 * year);
  }
}

static double seasonal(const double a[5], const double s[5]) {
  return a[0] + a[1] * s[1] + a[2] * s[2] + a[3] * s[3] + a[4] * s[4];
}

enum { VMF3_DEGREE = 12, VMF3_TERMS = 91 };

/* As harmonics9, to degree and order 12, the functions n, m at n (n + 1) / 2 + m. */
static void harmonics12(double latitude, double longitude, double v[VMF3_TERMS],
                        double w[VMF3_TERMS]) {
  const double cl = myb_cos(latitude);
  const double x = cl * myb_cos(longitude);
  const double y = cl * myb_sin(longitude);
  const double z = myb_sin(latitude);
#define AT(n, m) ((n) * ((n) + 1) / 2 + (m))
  v[AT(0, 0)] = 1.0;
  w[AT(0, 0)] = 0.0;
  v[AT(1, 0)] = z;
  w[AT(1, 0)] = 0.0;
  for (int n = 2; n <= VMF3_DEGREE; ++n) {
    v[AT(n, 0)] = ((2 * n - 1) * z * v[AT(n - 1, 0)] - (n - 1) * v[AT(n - 2, 0)]) / n;
    w[AT(n, 0)] = 0.0;
  }
  for (int m = 1; m <= VMF3_DEGREE; ++m) {
    v[AT(m, m)] = (2 * m - 1) * (x * v[AT(m - 1, m - 1)] - y * w[AT(m - 1, m - 1)]);
    w[AT(m, m)] = (2 * m - 1) * (x * w[AT(m - 1, m - 1)] + y * v[AT(m - 1, m - 1)]);
    if (m < VMF3_DEGREE) {
      v[AT(m + 1, m)] = (2 * m + 1) * z * v[AT(m, m)];
      w[AT(m + 1, m)] = (2 * m + 1) * z * w[AT(m, m)];
    }
    for (int n = m + 2; n <= VMF3_DEGREE; ++n) {
      v[AT(n, m)] =
          ((2 * n - 1) * z * v[AT(n - 1, m)] - (n + m - 1) * v[AT(n - 2, m)]) / (n - m);
      w[AT(n, m)] =
          ((2 * n - 1) * z * w[AT(n - 1, m)] - (n + m - 1) * w[AT(n - 2, m)]) / (n - m);
    }
  }
#undef AT
}

/* A coefficient b or c at the place and the season: the five expansions, then the season. */
static double vmf3_coefficient(const double cosine[VMF3_TERMS][5],
                               const double sine[VMF3_TERMS][5], const double* v,
                               const double* w, const double s[5]) {
  double part[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < VMF3_TERMS; ++i) {
    for (int k = 0; k < 5; ++k) {
      part[k] += cosine[i][k] * v[i] + sine[i][k] * w[i];
    }
  }
  return seasonal(part, s);
}

myb_status myb_tropo_vmf3(double ah, double aw, double mjd, double latitude, double longitude,
                          double height, double zenith_distance, double* mh, double* mw) {
  if (mh == NULL || mw == NULL || !finite_all(ah, aw, mjd, latitude) || !isfinite(longitude) ||
      !isfinite(height) || !zenith_ok(zenith_distance)) {
    return MYB_E_DOMAIN;
  }
  double day = 0.0;
  const myb_status status = year_day(mjd, &day);
  if (status != MYB_OK) {
    return status;
  }
  double s[5];
  seasons(day, 0, s);
  double v[VMF3_TERMS];
  double w[VMF3_TERMS];
  harmonics12(latitude, longitude, v, w);
  const double bh = vmf3_coefficient(vmf3_cos_bh, vmf3_sin_bh, v, w, s);
  const double bw = vmf3_coefficient(vmf3_cos_bw, vmf3_sin_bw, v, w, s);
  const double ch = vmf3_coefficient(vmf3_cos_ch, vmf3_sin_ch, v, w, s);
  const double cw = vmf3_coefficient(vmf3_cos_cw, vmf3_sin_cw, v, w, s);
  const double sine = myb_cos(zenith_distance);
  *mh = marini(sine, ah, bh, ch) + niell_height(sine, height);
  *mw = marini(sine, aw, bw, cw);
  return MYB_OK;
}

/* A node of the grid of GPT3: the seasonal quantities as the mean and the four amplitudes,
 * in the units of the model (the file has some of them a thousand or a hundred thousand
 * times larger). */
typedef struct gpt3_node {
  double pressure[5];    /* Pa */
  double temperature[5]; /* K */
  double humidity[5];    /* specific, kg/kg */
  double lapse[5];       /* K/m */
  double undulation;     /* m */
  double height;         /* of the node above the geoid, m */
  double ah[5];
  double aw[5];
  double lambda[5];
  double mean[5];        /* the mean temperature of the water vapour, K */
  double gradient[4][5]; /* north and east, hydrostatic; north and east, wet: m */
} gpt3_node;

struct myb_gpt3 {
  myb_allocator allocator;
  double step;    /* degrees */
  size_t rows;    /* of latitude, from the north */
  size_t columns; /* of longitude, from the first meridian eastwards */
  gpt3_node* nodes;
};

enum { GPT3_NUMBERS = 64 };

/* The numbers of a line; -1 for a field that is not a number or a wrong count. */
static int gpt3_line(const char* line, size_t length, double out[GPT3_NUMBERS]) {
  int n = 0;
  size_t at = 0;
  while (at < length) { /* bound: the characters of the line */
    if (line[at] == ' ' || line[at] == '\t' || line[at] == '\r') {
      ++at;
      continue;
    }
    size_t end = at;
    while (end < length && line[end] != ' ' && line[end] != '\t' && line[end] != '\r') {
      ++end;
    }
    size_t used = 0;
    double value = 0.0;
    if (n == GPT3_NUMBERS || myb_parse_double(line + at, end - at, &value, &used) != MYB_OK ||
        used != end - at || !isfinite(value)) {
      return -1;
    }
    out[n++] = value;
    at = end;
  }
  return n;
}

void myb_gpt3_destroy(myb_gpt3* grid) {
  if (grid == NULL) {
    return;
  }
  const myb_allocator a = grid->allocator;
  myb_internal_release(&a, grid->nodes);
  myb_internal_release(&a, grid);
}

static myb_status gpt3_fail(myb_gpt3* grid, myb_diag* diag, long line, const char* what) {
  myb_internal_diag(diag, line, "%s", what);
  myb_gpt3_destroy(grid);
  return MYB_E_FORMAT;
}

myb_status myb_gpt3_load(const char* text, size_t size, const myb_allocator* allocator,
                         myb_gpt3** grid, myb_diag* diag) {
  MYB_CHECK((text != NULL || size == 0) && grid != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_gpt3), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_gpt3* g = block;
  memset(g, 0, sizeof *g);
  g->allocator = a;
  /* two passes: count the nodes, then read them */
  size_t nodes = 0;
  for (int pass = 0; pass < 2; ++pass) {
    size_t at = 0;
    size_t node = 0;
    long number = 0;
    while (at < size) { /* bound: the lines of the text */
      size_t end = at;
      while (end < size && text[end] != '\n') {
        ++end;
      }
      const char* line = text + at;
      size_t length = end - at;
      at = end + 1;
      ++number;
      while (length > 0 && (line[0] == ' ' || line[0] == '\t')) {
        ++line;
        --length;
      }
      if (length == 0 || line[0] == '%' || line[0] == '\r') {
        continue; /* the heading and empty lines */
      }
      if (pass == 0) {
        ++node;
        continue;
      }
      double row[GPT3_NUMBERS];
      if (gpt3_line(line, length, row) != GPT3_NUMBERS) {
        return gpt3_fail(g, diag, number, "a node of GPT3 has 64 numbers");
      }
      const double latitude = 90.0 - g->step / 2.0 - (double)(node / g->columns) * g->step;
      const double longitude = g->step / 2.0 + (double)(node % g->columns) * g->step;
      /* the file writes the longitudes beyond 180 degrees as negative ones */
      const double east = row[1] < 0.0 ? row[1] + 360.0 : row[1];
      if (fabs(row[0] - latitude) > 1e-6 || fabs(east - longitude) > 1e-6) {
        return gpt3_fail(g, diag, number, "a node that is not the next of the grid");
      }
      gpt3_node* n = &g->nodes[node];
      for (int k = 0; k < 5; ++k) {
        n->pressure[k] = row[2 + k];
        n->temperature[k] = row[7 + k];
        n->humidity[k] = row[12 + k] / 1000.0;
        n->lapse[k] = row[17 + k] / 1000.0;
        n->ah[k] = row[24 + k] / 1000.0;
        n->aw[k] = row[29 + k] / 1000.0;
        n->lambda[k] = row[34 + k];
        n->mean[k] = row[39 + k];
        for (int q = 0; q < 4; ++q) {
          n->gradient[q][k] = row[44 + 5 * q + k] / 100000.0;
        }
      }
      n->undulation = row[22];
      n->height = row[23];
      ++node;
    }
    if (pass == 0) {
      /* a grid of equal steps over the whole body: 2 s^2 nodes for 180 / s rows */
      nodes = node;
      double rows = 0.0;
      if (nodes >= 2) {
        rows = floor(sqrt((double)nodes / 2.0) + 0.5);
      }
      if (rows < 2.0 || rows > 2000.0 || (size_t)(2.0 * rows * rows) != nodes) {
        return gpt3_fail(g, diag, number, "the nodes are not a grid of equal steps");
      }
      g->rows = (size_t)rows;
      g->columns = 2 * g->rows;
      g->step = 180.0 / rows;
      status = myb_internal_alloc_array(&a, nodes, sizeof(gpt3_node), &block);
      if (status != MYB_OK) {
        myb_gpt3_destroy(g);
        return status;
      }
      g->nodes = block;
    }
  }
  *grid = g;
  return MYB_OK;
}

double myb_gpt3_step(const myb_gpt3* grid) {
  MYB_CHECK(grid != NULL);
  return grid->step * (TWO_PI / 360.0);
}

/* What a node gives at the height of the point and the season. */
typedef struct gpt3_point {
  double v[13]; /* in the order of myb_gpt3_values */
} gpt3_point;

static void gpt3_at_node(const gpt3_node* n, double height, const double s[5], gpt3_point* out) {
  const double gm = 9.80665;
  const double molar = 28.965e-3;
  const double gas = 8.3143;
  const double t0 = seasonal(n->temperature, s);
  const double p0 = seasonal(n->pressure, s);
  const double q = seasonal(n->humidity, s);
  const double lapse = seasonal(n->lapse, s);
  const double lambda = seasonal(n->lambda, s);
  const double above = height - n->undulation - n->height;
  const double virtual_temperature = t0 * (1.0 + 0.6077 * q);
  const double pressure = p0 * myb_exp(-gm * molar / (gas * virtual_temperature) * above) / 100.0;
  const double e0 = q * p0 / (0.622 + 0.378 * q) / 100.0;
  out->v[0] = pressure;
  out->v[1] = t0 + lapse * above - 273.15;
  out->v[2] = lapse * 1000.0;
  out->v[3] = seasonal(n->mean, s);
  /* the vapour pressure at the height of the point (Askne and Nordius 1987, eq. 14) */
  out->v[4] = e0 * myb_exp((lambda + 1.0) * myb_log(100.0 * pressure / p0));
  out->v[5] = seasonal(n->ah, s);
  out->v[6] = seasonal(n->aw, s);
  out->v[7] = lambda;
  out->v[8] = n->undulation;
  for (int k = 0; k < 4; ++k) {
    out->v[9 + k] = seasonal(n->gradient[k], s);
  }
}

myb_status myb_gpt3_at(const myb_gpt3* grid, double mjd, double latitude, double longitude,
                       double height, int constant, myb_gpt3_values* values) {
  MYB_CHECK(grid != NULL && values != NULL);
  if (!isfinite(mjd) || !isfinite(longitude) || !isfinite(height) ||
      !(fabs(latitude) <= 1.5707963267948966)) {
    return MYB_E_DOMAIN;
  }
  double day = 0.0;
  const myb_status status = year_day(mjd, &day);
  if (status != MYB_OK) {
    return status;
  }
  double s[5];
  seasons(day, constant, s);
  const double step = grid->step;
  const double degrees = 360.0 / TWO_PI;
  double east = longitude * degrees;
  east -= 360.0 * floor(east / 360.0);
  const double polar = 90.0 - latitude * degrees;
  /* the node of the cell the point is in, counted from 1, and where the point is in it */
  double row = floor((polar + step) / step);
  double column = floor((east + step) / step);
  const double across = (polar - (row * step - step / 2.0)) / step;
  const double along = (east - (column * step - step / 2.0)) / step;
  const double rows = (double)grid->rows;
  const double columns = (double)grid->columns;
  if (row > rows) {
    row = rows;
  }
  if (column > columns) {
    column = 1.0;
  }
  const size_t r0 = (size_t)row - 1;
  const size_t c0 = (size_t)column - 1;
  gpt3_point p[4];
  gpt3_at_node(&grid->nodes[r0 * grid->columns + c0], height, s, &p[0]);
  double out[13];
  if (polar > step / 2.0 && polar < 180.0 - step / 2.0) {
    /* between four nodes: the next towards the point in latitude and in longitude */
    const size_t r1 = across < 0.0 ? r0 - 1 : r0 + 1;
    const size_t c1 = along < 0.0 ? (c0 == 0 ? grid->columns - 1 : c0 - 1)
                                  : (c0 + 1 == grid->columns ? 0 : c0 + 1);
    gpt3_at_node(&grid->nodes[r1 * grid->columns + c0], height, s, &p[1]);
    gpt3_at_node(&grid->nodes[r0 * grid->columns + c1], height, s, &p[2]);
    gpt3_at_node(&grid->nodes[r1 * grid->columns + c1], height, s, &p[3]);
    const double far_row = fabs(across);
    const double near_row = 1.0 - far_row;
    const double far_column = fabs(along);
    const double near_column = 1.0 - far_column;
    for (int k = 0; k < 13; ++k) {
      const double first = near_row * p[0].v[k] + far_row * p[1].v[k];
      const double second = near_row * p[2].v[k] + far_row * p[3].v[k];
      out[k] = near_column * first + far_column * second;
    }
  } else {
    memcpy(out, p[0].v, sizeof out); /* by the pole: the node itself */
  }
  values->pressure = out[0];
  values->temperature = out[1];
  values->lapse_rate = out[2];
  values->mean_temperature = out[3];
  values->vapour_pressure = out[4];
  values->ah = out[5];
  values->aw = out[6];
  values->lambda = out[7];
  values->undulation = out[8];
  values->gn_h = out[9];
  values->ge_h = out[10];
  values->gn_w = out[11];
  values->ge_w = out[12];
  return MYB_OK;
}
