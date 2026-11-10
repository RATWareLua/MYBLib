/* SPDX-License-Identifier: BSD-3-Clause */
/* IONEX maps of the vertical total electron content and the single-layer model
 * (myb_ionosphere.h). Written from the IONEX description (Schaer, Gurtner, Feltens 1998). */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/check.h"
#include "myb/cmath.h"
#include "myb/io.h"
#include "myb_ionosphere.h"

#define DEG (3.14159265358979323846 / 180.0)
#define TECU 1e16

struct myb_ionex {
  myb_allocator allocator;
  double radius; /* m */
  double height; /* m */
  double lat1;   /* degrees */
  double dlat;
  double lon1;
  double dlon;
  int nlat;
  int nlon;
  int nmaps;
  int global; /* the longitudes cover 360 degrees */
  double* epochs; /* MJD of each map */
  double* tec;    /* electrons/m^2, [map][lat][lon]; NaN where missing */
};

myb_status myb_iono_pierce_point(double latitude, double longitude, double azimuth,
                                 double elevation, double radius, double height, double alpha,
                                 double* pierce_latitude, double* pierce_longitude,
                                 double* mapping) {
  if (pierce_latitude == NULL || pierce_longitude == NULL || mapping == NULL ||
      !isfinite(latitude + longitude + azimuth + elevation + radius + height + alpha) ||
      !(elevation > 0.0 && elevation <= 1.5707963267948966) || !(radius > 0.0) ||
      !(height > 0.0) || !(alpha > 0.0 && alpha <= 1.0) || !(fabs(latitude) <= 1.5707963267948966)) {
    return MYB_E_DOMAIN;
  }
  const double z = 1.5707963267948966 - elevation;
  const double k = radius / (radius + height);
  const double sz = k * myb_sin(alpha * z);
  *mapping = 1.0 / sqrt(1.0 - sz * sz);
  const double psi = z - myb_asin(k * myb_sin(z));
  const double sl = myb_sin(latitude);
  const double cl = myb_cos(latitude);
  const double sp = myb_sin(psi);
  const double cp = myb_cos(psi);
  const double s = sl * cp + cl * sp * myb_cos(azimuth);
  const double lat = myb_asin(fmax(-1.0, fmin(1.0, s)));
  *pierce_latitude = lat;
  *pierce_longitude = longitude + myb_atan2(myb_sin(azimuth) * sp * cl, cp - sl * s);
  return MYB_OK;
}

/* ---- IONEX ----------------------------------------------------------------------------------- */

/* days since 1858-11-17 of a Gregorian date (H. Hinnant's days_from_civil) plus the time */
static double mjd_of(long y, long m, long d, long hh, long mm, double ss) {
  y -= m <= 2;
  const long era = (y >= 0 ? y : y - 399) / 400;
  const long yoe = y - era * 400;
  const long doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
  const long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  const long days = era * 146097 + doe - 719468;
  return (double)(days + 40587) + ((double)hh * 3600.0 + (double)mm * 60.0 + ss) / 86400.0;
}

/* the fixed-width field [start, start + width) of a line (blanks allowed around the number) */
static int field(const char* line, size_t len, size_t start, size_t width, double* value) {
  if (start >= len) {
    return 0;
  }
  size_t a = start;
  size_t b = start + width < len ? start + width : len;
  while (a < b && line[a] == ' ') { /* bound: width */
    ++a;
  }
  while (b > a && (line[b - 1] == ' ' || line[b - 1] == '\r')) { /* bound: width */
    --b;
  }
  if (a == b) {
    return 0;
  }
  size_t used = 0;
  return myb_parse_double(line + a, b - a, value, &used) == MYB_OK && used == b - a;
}

static int label_is(const char* line, size_t len, const char* label) {
  const size_t n = strlen(label);
  if (len < 60 + n) {
    return 0;
  }
  return memcmp(line + 60, label, n) == 0;
}

static int read_epoch(const char* line, size_t len, double* mjd) {
  double f[6];
  for (int i = 0; i < 6; ++i) {
    if (!field(line, len, (size_t)(6 * i), 6, &f[i])) {
      return 0;
    }
  }
  *mjd = mjd_of((long)f[0], (long)f[1], (long)f[2], (long)f[3], (long)f[4], f[5]);
  return 1;
}

static long count_of(double first, double last, double step) {
  const double n = (last - first) / step;
  const double r = floor(n + 0.5);
  if (!(step != 0.0) || !(r >= 0.0) || fabs(n - r) > 1e-6 || r > 100000.0) {
    return -1;
  }
  return (long)r + 1;
}

myb_status myb_ionex_load(const char* text, size_t size, const myb_allocator* allocator,
                          myb_ionex** maps, myb_diag* diag) {
  MYB_CHECK(text != NULL && maps != NULL);
  myb_internal_diag_clear(diag);
  *maps = NULL;
  const myb_allocator a = myb_internal_allocator(allocator);
  myb_ionex* x = NULL;
  void* block = NULL;
  myb_status st = myb_internal_alloc_array(&a, 1, sizeof(myb_ionex), &block);
  if (st != MYB_OK) {
    return st;
  }
  x = block;
  memset(x, 0, sizeof *x);
  x->allocator = a;
  double hgt[3] = {0.0, 0.0, 0.0};
  double lat[3] = {0.0, 0.0, 0.0};
  double lon[3] = {0.0, 0.0, 0.0};
  int have = 0; /* bits: 1 maps, 2 radius, 4 hgt, 8 lat, 16 lon */
  double dimension = 2.0;
  double exponent = -1.0;
  int in_header = 1;
  int in_tec = 0;    /* inside a TEC map */
  int in_other = 0;  /* inside an RMS or height map */
  long map = -1;
  long row = -1;
  long col = 0;
  long line_number = 0;
  size_t pos = 0;
  for (size_t guard = 0; guard <= size && pos < size; ++guard) { /* bound: the lines */
    const char* line = text + pos;
    const char* nl = memchr(line, '\n', size - pos);
    const size_t len = nl != NULL ? (size_t)(nl - line) : size - pos;
    pos += len + (nl != NULL ? 1 : 0);
    ++line_number;
    if (in_header) {
      if (label_is(line, len, "END OF HEADER")) {
        if ((have & 31) != 31 || dimension != 2.0) {
          myb_internal_diag(diag, line_number,
                            dimension != 2.0 ? "three-dimensional maps are not supported"
                                             : "incomplete header");
          st = dimension != 2.0 ? MYB_E_DOMAIN : MYB_E_FORMAT;
          goto fail;
        }
        const long nlat = count_of(lat[0], lat[1], lat[2]);
        const long nlon = count_of(lon[0], lon[1], lon[2]);
        if (nlat < 2 || nlon < 2 || x->nmaps < 1) {
          myb_internal_diag(diag, line_number, "invalid grid");
          st = MYB_E_FORMAT;
          goto fail;
        }
        x->nlat = (int)nlat;
        x->nlon = (int)nlon;
        x->lat1 = lat[0];
        x->dlat = lat[2];
        x->lon1 = lon[0];
        x->dlon = lon[2];
        x->global = fabs(fabs(lon[1] - lon[0]) - 360.0) < 1e-9;
        x->height = hgt[0] * 1000.0;
        size_t cells = 0;
        if (__builtin_mul_overflow((size_t)x->nmaps, (size_t)nlat * (size_t)nlon, &cells)) {
          st = MYB_E_OVERFLOW;
          goto fail;
        }
        void* e = NULL;
        st = myb_internal_alloc_array(&a, (size_t)x->nmaps, sizeof(double), &e);
        if (st != MYB_OK) {
          goto fail;
        }
        x->epochs = e;
        void* t = NULL;
        st = myb_internal_alloc_array(&a, cells, sizeof(double), &t);
        if (st != MYB_OK) {
          goto fail;
        }
        x->tec = t;
        for (size_t i = 0; i < cells; ++i) {
          x->tec[i] = (double)NAN;
        }
        for (int i = 0; i < x->nmaps; ++i) {
          x->epochs[i] = (double)NAN;
        }
        in_header = 0;
      } else if (label_is(line, len, "# OF MAPS IN FILE")) {
        double n = 0.0;
        if (!field(line, len, 0, 6, &n) || !(n >= 1.0 && n <= 10000.0)) {
          myb_internal_diag(diag, line_number, "invalid number of maps");
          st = MYB_E_FORMAT;
          goto fail;
        }
        x->nmaps = (int)n;
        have |= 1;
      } else if (label_is(line, len, "BASE RADIUS")) {
        double r = 0.0;
        if (!field(line, len, 0, 8, &r) || !(r > 0.0)) {
          myb_internal_diag(diag, line_number, "invalid base radius");
          st = MYB_E_FORMAT;
          goto fail;
        }
        x->radius = r * 1000.0;
        have |= 2;
      } else if (label_is(line, len, "MAP DIMENSION")) {
        if (!field(line, len, 0, 6, &dimension)) {
          st = MYB_E_FORMAT;
          goto fail;
        }
      } else if (label_is(line, len, "EXPONENT")) {
        if (!field(line, len, 0, 6, &exponent)) {
          st = MYB_E_FORMAT;
          goto fail;
        }
      } else if (label_is(line, len, "HGT1 / HGT2 / DHGT") ||
                 label_is(line, len, "LAT1 / LAT2 / DLAT") ||
                 label_is(line, len, "LON1 / LON2 / DLON")) {
        double* f = label_is(line, len, "HGT1") ? hgt : label_is(line, len, "LAT1") ? lat : lon;
        for (int i = 0; i < 3; ++i) {
          if (!field(line, len, (size_t)(2 + 6 * i), 6, &f[i])) {
            myb_internal_diag(diag, line_number, "invalid grid definition");
            st = MYB_E_FORMAT;
            goto fail;
          }
        }
        have |= f == hgt ? 4 : f == lat ? 8 : 16;
      }
      continue;
    }
    if (in_other) {
      if (label_is(line, len, "END OF RMS MAP") || label_is(line, len, "END OF HEIGHT MAP")) {
        in_other = 0;
      }
      continue;
    }
    if (!in_tec) {
      if (label_is(line, len, "START OF TEC MAP")) {
        double n = 0.0;
        if (!field(line, len, 0, 6, &n) || !(n >= 1.0 && n <= (double)x->nmaps)) {
          myb_internal_diag(diag, line_number, "invalid map number");
          st = MYB_E_FORMAT;
          goto fail;
        }
        map = (long)n - 1;
        in_tec = 1;
        row = -1;
      } else if (label_is(line, len, "START OF RMS MAP") ||
                 label_is(line, len, "START OF HEIGHT MAP")) {
        in_other = 1;
      } else if (label_is(line, len, "EXPONENT")) {
        if (!field(line, len, 0, 6, &exponent)) {
          st = MYB_E_FORMAT;
          goto fail;
        }
      }
      continue;
    }
    /* inside a TEC map */
    if (label_is(line, len, "END OF TEC MAP")) {
      in_tec = 0;
      continue;
    }
    if (label_is(line, len, "EPOCH OF CURRENT MAP")) {
      if (!read_epoch(line, len, &x->epochs[map])) {
        myb_internal_diag(diag, line_number, "invalid epoch");
        st = MYB_E_FORMAT;
        goto fail;
      }
      continue;
    }
    if (label_is(line, len, "EXPONENT")) {
      if (!field(line, len, 0, 6, &exponent)) {
        st = MYB_E_FORMAT;
        goto fail;
      }
      continue;
    }
    if (label_is(line, len, "LAT/LON1/LON2/DLON/H")) {
      double f[5];
      for (int i = 0; i < 5; ++i) {
        if (!field(line, len, (size_t)(2 + 6 * i), 6, &f[i])) {
          myb_internal_diag(diag, line_number, "invalid latitude record");
          st = MYB_E_FORMAT;
          goto fail;
        }
      }
      const double q = (f[0] - x->lat1) / x->dlat;
      row = (long)floor(q + 0.5);
      if (row < 0 || row >= x->nlat || fabs(q - (double)row) > 1e-6 || f[1] != lon[0] ||
          f[2] != lon[1] || f[3] != lon[2]) {
        myb_internal_diag(diag, line_number, "latitude record off the grid");
        st = MYB_E_FORMAT;
        goto fail;
      }
      col = 0;
      continue;
    }
    if (row < 0) {
      myb_internal_diag(diag, line_number, "values before a latitude record");
      st = MYB_E_FORMAT;
      goto fail;
    }
    /* 16I5 values */
    /* values in units of 10^exponent TECU: an exact power of ten, one division or product */
    double power = 1.0;
    for (int k = 0; k < (int)fabs(exponent) && k < 22; ++k) {
      power *= 10.0;
    }
    for (size_t k = 0; k < 16 && col < x->nlon; ++k) {
      double v = 0.0;
      if (!field(line, len, 5 * k, 5, &v)) {
        break;
      }
      const size_t idx = ((size_t)map * (size_t)x->nlat + (size_t)row) * (size_t)x->nlon +
                         (size_t)col;
      const double tecu = exponent < 0.0 ? v / power : v * power;
      x->tec[idx] = v == 9999.0 ? (double)NAN : tecu * TECU;
      ++col;
    }
  }
  if (in_header) {
    myb_internal_diag(diag, line_number, "no END OF HEADER");
    st = MYB_E_FORMAT;
    goto fail;
  }
  for (int i = 0; i < x->nmaps; ++i) {
    if (!isfinite(x->epochs[i]) || (i > 0 && !(x->epochs[i] > x->epochs[i - 1]))) {
      myb_internal_diag(diag, 0, "map %d missing or out of order", i + 1);
      st = MYB_E_FORMAT;
      goto fail;
    }
  }
  *maps = x;
  return MYB_OK;
fail:
  myb_ionex_destroy(x);
  return st;
}

void myb_ionex_destroy(myb_ionex* maps) {
  if (maps == NULL) {
    return;
  }
  const myb_allocator a = maps->allocator;
  myb_internal_release(&a, maps->epochs);
  myb_internal_release(&a, maps->tec);
  myb_internal_release(&a, maps);
}

void myb_ionex_shell(const myb_ionex* maps, double* radius, double* height) {
  MYB_CHECK(maps != NULL && radius != NULL && height != NULL);
  *radius = maps->radius;
  *height = maps->height;
}

/* bilinear value of map m at latitude, longitude in degrees; NaN outside or missing */
static double map_value(const myb_ionex* x, int m, double lat_deg, double lon_deg) {
  double lon = lon_deg;
  if (x->global) {
    const double lo = fmin(x->lon1, x->lon1 + (x->nlon - 1) * x->dlon);
    lon = lo + (lon - lo) - 360.0 * floor((lon - lo) / 360.0);
  }
  const double p = (lon - x->lon1) / x->dlon;
  const double q = (lat_deg - x->lat1) / x->dlat;
  if (!(p >= -1e-12 && p <= x->nlon - 1 + 1e-12) || !(q >= -1e-12 && q <= x->nlat - 1 + 1e-12)) {
    return (double)NAN;
  }
  int i = (int)floor(p);
  int j = (int)floor(q);
  i = i < 0 ? 0 : i > x->nlon - 2 ? x->nlon - 2 : i;
  j = j < 0 ? 0 : j > x->nlat - 2 ? x->nlat - 2 : j;
  const double fp = p - i;
  const double fq = q - j;
  const double* t = x->tec + (size_t)m * (size_t)x->nlat * (size_t)x->nlon;
  const size_t n = (size_t)x->nlon;
  const size_t i0 = (size_t)i;
  const size_t j0 = (size_t)j;
  const double e00 = t[j0 * n + i0];
  const double e10 = t[j0 * n + i0 + 1];
  const double e01 = t[(j0 + 1) * n + i0];
  const double e11 = t[(j0 + 1) * n + i0 + 1];
  return (1.0 - fp) * (1.0 - fq) * e00 + fp * (1.0 - fq) * e10 + fq * (1.0 - fp) * e01 +
         fp * fq * e11;
}

myb_status myb_ionex_vtec(const myb_ionex* maps, double mjd, double latitude, double longitude,
                          double* vtec) {
  MYB_CHECK(maps != NULL && vtec != NULL);
  if (!isfinite(mjd + latitude + longitude)) {
    return MYB_E_DOMAIN;
  }
  const myb_ionex* x = maps;
  if (!(mjd >= x->epochs[0] && mjd <= x->epochs[x->nmaps - 1]) || x->nmaps < 1) {
    return MYB_E_RANGE;
  }
  int i = 0;
  for (; i + 2 < x->nmaps && mjd > x->epochs[i + 1]; ++i) { /* bound: nmaps */
  }
  const double lat = latitude / DEG;
  const double lon = longitude / DEG;
  double value = 0.0;
  if (x->nmaps == 1) {
    value = map_value(x, 0, lat, lon);
  } else {
    const double t0 = x->epochs[i];
    const double t1 = x->epochs[i + 1];
    const double e0 = map_value(x, i, lat, lon + (mjd - t0) * 360.0);
    const double e1 = map_value(x, i + 1, lat, lon + (mjd - t1) * 360.0);
    value = ((t1 - mjd) * e0 + (mjd - t0) * e1) / (t1 - t0);
  }
  if (!isfinite(value)) {
    return MYB_E_RANGE;
  }
  *vtec = value;
  return MYB_OK;
}
