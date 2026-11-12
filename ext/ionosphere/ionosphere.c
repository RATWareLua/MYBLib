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
#define TWO_PI 6.283185307179586476925287

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
  int dimension;  /* 2 or 3 */
  int nhgt;       /* the heights of three-dimensional maps; 1 */
  double dhgt;    /* their step, m */
  double* epochs; /* MJD of each map */
  double* tec;    /* electrons/m^2, [map][lat][lon]; of three-dimensional maps the densities,
                     electrons/m^3, [map][height][lat][lon]; NaN where missing */
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
  long level = 0;
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
        if ((have & 31) != 31 || (dimension != 2.0 && dimension != 3.0)) {
          myb_internal_diag(diag, line_number,
                            (have & 31) != 31 ? "incomplete header" : "a map dimension not 2 or 3");
          st = MYB_E_FORMAT;
          goto fail;
        }
        const long nlat = count_of(lat[0], lat[1], lat[2]);
        const long nlon = count_of(lon[0], lon[1], lon[2]);
        const long nhgt = dimension == 3.0 ? count_of(hgt[0], hgt[1], hgt[2]) : 1;
        if (nlat < 2 || nlon < 2 || x->nmaps < 1 ||
            (dimension == 3.0 && (nhgt < 2 || !(hgt[2] > 0.0)))) {
          myb_internal_diag(diag, line_number, "invalid grid");
          st = MYB_E_FORMAT;
          goto fail;
        }
        x->dimension = (int)dimension;
        x->nhgt = (int)nhgt;
        x->dhgt = hgt[2] * 1000.0;
        x->nlat = (int)nlat;
        x->nlon = (int)nlon;
        x->lat1 = lat[0];
        x->dlat = lat[2];
        x->lon1 = lon[0];
        x->dlon = lon[2];
        x->global = fabs(fabs(lon[1] - lon[0]) - 360.0) < 1e-9;
        x->height = hgt[0] * 1000.0;
        size_t cells = 0;
        if (__builtin_mul_overflow((size_t)x->nmaps * (size_t)nhgt, (size_t)nlat * (size_t)nlon,
                                   &cells)) {
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
        level = 0;
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
      if (x->dimension == 3) {
        const double z = (f[4] - hgt[0]) / hgt[2];
        level = (long)floor(z + 0.5);
        if (level < 0 || level >= x->nhgt || fabs(z - (double)level) > 1e-6) {
          myb_internal_diag(diag, line_number, "a height off the grid");
          st = MYB_E_FORMAT;
          goto fail;
        }
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
      const size_t idx =
          (((size_t)map * (size_t)x->nhgt + (size_t)level) * (size_t)x->nlat + (size_t)row) *
              (size_t)x->nlon +
          (size_t)col;
      const double tecu = exponent < 0.0 ? v / power : v * power;
      /* of three-dimensional maps the density at the node times DHGT */
      const double content = tecu * TECU;
      x->tec[idx] = v == 9999.0 ? (double)NAN : x->dimension == 3 ? content / x->dhgt : content;
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

int myb_ionex_dimension(const myb_ionex* maps) {
  MYB_CHECK(maps != NULL);
  return maps->dimension;
}

/* bilinear value of level k of map m at latitude, longitude in degrees; NaN outside or
 * missing */
static double level_value(const myb_ionex* x, int m, int k, double lat_deg, double lon_deg) {
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
  const double* t =
      x->tec + ((size_t)m * (size_t)x->nhgt + (size_t)k) * (size_t)x->nlat * (size_t)x->nlon;
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

/* the vertical content of map m: of three-dimensional maps the sum of the densities times
 * DHGT over the heights */
static double map_value(const myb_ionex* x, int m, double lat_deg, double lon_deg) {
  if (x->dimension != 3) {
    return level_value(x, m, 0, lat_deg, lon_deg);
  }
  double sum = 0.0;
  for (int k = 0; k < x->nhgt; ++k) {
    sum += level_value(x, m, k, lat_deg, lon_deg);
  }
  return sum * x->dhgt;
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
    /* at the epoch of a map the other one has no part: a value missing there is not needed */
    if (!isfinite(value)) {
      if (mjd == t0 && isfinite(e0)) {
        value = e0;
      } else if (mjd == t1 && isfinite(e1)) {
        value = e1;
      }
    }
  }
  if (!isfinite(value)) {
    return MYB_E_RANGE;
  }
  *vtec = value;
  return MYB_OK;
}

/* ---- The density of electrons in space -------------------------------------------------------- */

/* The geocentric latitude and longitude (rad) and the distance of a place. */
static int spherical(const double r[3], double* latitude, double* longitude, double* distance) {
  const double h = sqrt(r[0] * r[0] + r[1] * r[1]);
  const double d = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  if (!(d > 0.0) || !isfinite(d)) {
    return 0;
  }
  *latitude = myb_atan2(r[2], h);
  *longitude = h > 0.0 ? myb_atan2(r[1], r[0]) : 0.0;
  *distance = d;
  return 1;
}

myb_status myb_iono_layer_density(void* layer, double mjd, const double r[3], double* density) {
  const myb_iono_layer* l = layer;
  MYB_CHECK(l != NULL && l->maps != NULL && r != NULL && density != NULL);
  double latitude = 0.0;
  double longitude = 0.0;
  double distance = 0.0;
  if (!(l->scale > 0.0) || !isfinite(l->scale) || !isfinite(l->peak) ||
      !spherical(r, &latitude, &longitude, &distance)) {
    return MYB_E_DOMAIN;
  }
  double vtec = 0.0;
  const myb_status status = myb_ionex_vtec(l->maps, mjd, latitude, longitude, &vtec);
  if (status != MYB_OK) {
    return status;
  }
  const double z = (distance - l->maps->radius - l->peak) / l->scale;
  /* sqrt(2 pi e): the integral of the layer over z */
  *density = vtec / (l->scale * 4.132731354122493) * myb_exp(0.5 * (1.0 - z - myb_exp(-z)));
  return MYB_OK;
}

/* The value of the map m of a grid at a height and a place, in steps of the grid. */
static myb_status grid_value(const myb_iono_grid* g, int m, double h, double q, double p,
                             double* value) {
  if (!(q >= -1e-9 && q <= g->latitudes - 1 + 1e-9)) {
    return MYB_E_RANGE;
  }
  /* round: the columns span a turn, or a turn and a step with the last the first again */
  const double tiny = 1e-9 * g->longitude_step;
  const int closed = fabs((double)(g->longitudes - 1) * g->longitude_step - TWO_PI) < tiny;
  const int round = closed || fabs((double)g->longitudes * g->longitude_step - TWO_PI) < tiny;
  const int period = closed ? g->longitudes - 1 : g->longitudes;
  const double turn = (double)period;
  if (round) {
    p -= turn * floor(p / turn);
  } else if (!(p >= -1e-9 && p <= turn - 1.0 + 1e-9)) {
    return MYB_E_RANGE;
  }
  int k = (int)floor(h);
  int j = (int)floor(q);
  int i = (int)floor(p);
  k = k < 0 ? 0 : k > g->heights - 2 ? g->heights - 2 : k;
  j = j < 0 ? 0 : j > g->latitudes - 2 ? g->latitudes - 2 : j;
  i = i < 0 ? 0 : i > period - (round ? 1 : 2) ? period - (round ? 1 : 2) : i;
  const int i1 = i + 1 == period ? 0 : i + 1;
  const double fh = h - k;
  const double fq = q - j;
  const double fp = p - i;
  const size_t row = (size_t)g->longitudes;
  const size_t level = row * (size_t)g->latitudes;
  const double* d = g->density + (size_t)m * level * (size_t)g->heights;
  double sum = 0.0;
  for (int a = 0; a < 2; ++a) {
    for (int b = 0; b < 2; ++b) {
      const double* line = d + (size_t)(k + a) * level + (size_t)(j + b) * row;
      const double along = (1.0 - fp) * line[i] + fp * line[i1];
      sum += (a ? fh : 1.0 - fh) * (b ? fq : 1.0 - fq) * along;
    }
  }
  if (!isfinite(sum)) {
    return MYB_E_RANGE; /* a missing value around */
  }
  *value = sum;
  return MYB_OK;
}

myb_status myb_ionex_density_grid(const myb_ionex* maps, myb_iono_grid* grid) {
  MYB_CHECK(maps != NULL && grid != NULL);
  if (maps->dimension != 3 || !(maps->dlon > 0.0)) {
    return MYB_E_DOMAIN;
  }
  grid->radius = maps->radius;
  grid->latitude = maps->lat1 * DEG;
  grid->latitude_step = maps->dlat * DEG;
  grid->latitudes = maps->nlat;
  grid->longitude = maps->lon1 * DEG;
  grid->longitude_step = maps->dlon * DEG;
  grid->longitudes = maps->nlon;
  grid->height = maps->height;
  grid->height_step = maps->dhgt;
  grid->heights = maps->nhgt;
  grid->epochs = maps->epochs;
  grid->maps = maps->nmaps;
  grid->density = maps->tec;
  return MYB_OK;
}

myb_status myb_iono_grid_density(void* grid, double mjd, const double r[3], double* density) {
  const myb_iono_grid* g = grid;
  MYB_CHECK(g != NULL && r != NULL && density != NULL);
  double latitude = 0.0;
  double longitude = 0.0;
  double distance = 0.0;
  if (g->epochs == NULL || g->density == NULL || g->maps < 1 || g->latitudes < 2 ||
      g->longitudes < 2 || g->heights < 2 || !(g->longitude_step > 0.0) ||
      !(g->height_step > 0.0) || g->latitude_step == 0.0 || !isfinite(mjd) ||
      !spherical(r, &latitude, &longitude, &distance)) {
    return MYB_E_DOMAIN;
  }
  if (!(mjd >= g->epochs[0] && mjd <= g->epochs[g->maps - 1])) {
    return MYB_E_RANGE;
  }
  const double h = (distance - g->radius - g->height) / g->height_step;
  if (h < 0.0 || h > (double)(g->heights - 1)) {
    *density = 0.0;
    return MYB_OK;
  }
  const double q = (latitude - g->latitude) / g->latitude_step;
  const double p = (longitude - g->longitude) / g->longitude_step;
  int i = 0;
  for (; i + 2 < g->maps && mjd > g->epochs[i + 1]; ++i) { /* bound: the maps */
  }
  if (g->maps == 1) {
    return grid_value(g, 0, h, q, p, density);
  }
  /* each map rotated with the Sun to the date, as the maps of the vertical content */
  const double t0 = g->epochs[i];
  const double t1 = g->epochs[i + 1];
  const double turn = TWO_PI / g->longitude_step; /* steps of longitude in a day */
  double first = 0.0;
  double second = 0.0;
  const myb_status s0 = grid_value(g, i, h, q, p + (mjd - t0) * turn, &first);
  const myb_status s1 = grid_value(g, i + 1, h, q, p + (mjd - t1) * turn, &second);
  if (s0 == MYB_OK && s1 == MYB_OK) {
    *density = ((t1 - mjd) * first + (mjd - t0) * second) / (t1 - t0);
    return MYB_OK;
  }
  /* at the epoch of a map the other one has no part: a value missing there is not needed */
  if (mjd == t0 && s0 == MYB_OK) {
    *density = first;
    return MYB_OK;
  }
  if (mjd == t1 && s1 == MYB_OK) {
    *density = second;
    return MYB_OK;
  }
  return s0 != MYB_OK ? s0 : s1;
}

/* ---- The integrals along a ray --------------------------------------------------------------- */

myb_status myb_iono_ray_integrals(myb_electron_density density, void* context,
                                  const myb_magnetic_field* field, myb_gravity_work* work,
                                  double mjd, const double receiver[3],
                                  const double transmitter[3], double bottom, double top,
                                  int panels, myb_iono_ray* ray) {
  MYB_CHECK(density != NULL && receiver != NULL && transmitter != NULL && ray != NULL &&
            (field == NULL || work != NULL));
  double u[3];
  double length2 = 0.0;
  double along = 0.0; /* receiver . u */
  double reach2 = 0.0;
  for (int i = 0; i < 3; ++i) {
    u[i] = transmitter[i] - receiver[i];
    length2 += u[i] * u[i];
    reach2 += receiver[i] * receiver[i];
  }
  const double length = sqrt(length2);
  if (!(length > 0.0) || !isfinite(length) || !isfinite(reach2) || !isfinite(mjd) ||
      !(bottom > 0.0) || !(top > bottom) || !isfinite(top) || panels < 1 || panels > 4096) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    u[i] /= length;
    along += receiver[i] * u[i];
  }
  /* |receiver + s u|^2 = reach2 + 2 along s + s^2: inside the sphere of radius rho for s
   * between -along -+ sqrt(along^2 - reach2 + rho^2) */
  const double outer = along * along - reach2 + top * top;
  *ray = (myb_iono_ray){0.0, 0.0, 0.0, 0.0};
  if (!(outer > 0.0)) {
    return MYB_OK; /* the line passes above the ionosphere */
  }
  const double far = sqrt(outer);
  double part[2][2] = {{fmax(0.0, -along - far), fmin(length, -along + far)}, {0.0, 0.0}};
  int parts = part[0][1] > part[0][0] ? 1 : 0;
  const double inner = along * along - reach2 + bottom * bottom;
  if (parts == 1 && inner > 0.0) {
    /* what is below the ionosphere is taken out: one part may become two */
    const double near = sqrt(inner);
    const double in = -along - near;
    const double out = -along + near;
    const double begin = part[0][0];
    const double end = part[0][1];
    parts = 0;
    if (in > begin) {
      part[parts][0] = begin;
      part[parts][1] = fmin(end, in);
      ++parts;
    }
    if (out < end) {
      part[parts][0] = fmax(begin, out);
      part[parts][1] = end;
      ++parts;
    }
  }
  /* the nodes and the weights of the four-point Gauss rule on (-1, 1) */
  static const double node[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
                                 0.8611363115940526};
  static const double weight[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
                                   0.3478548451374538};
  double sum[4] = {0.0, 0.0, 0.0, 0.0};
  for (int k = 0; k < parts; ++k) {
    const double half = 0.5 * (part[k][1] - part[k][0]) / (double)panels;
    for (int panel = 0; panel < panels; ++panel) {
      const double middle = part[k][0] + (2.0 * (double)panel + 1.0) * half;
      for (int q = 0; q < 4; ++q) {
        const double s = middle + half * node[q];
        const double r[3] = {receiver[0] + s * u[0], receiver[1] + s * u[1],
                             receiver[2] + s * u[2]};
        double n = 0.0;
        myb_status status = density(context, mjd, r, &n);
        if (status != MYB_OK) {
          return status;
        }
        if (!isfinite(n)) {
          return MYB_E_DOMAIN;
        }
        const double w = weight[q] * half;
        sum[0] += w * n;
        sum[2] += w * n * n;
        if (field != NULL && n != 0.0) {
          double b[3];
          status = myb_magnetic_field_eval(field, work, r, b);
          if (status != MYB_OK) {
            return status;
          }
          /* the direction of propagation is -u: transmitter -> receiver */
          const double cosine = -(b[0] * u[0] + b[1] * u[1] + b[2] * u[2]);
          const double b2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
          sum[1] += w * n * cosine;          /* B cos(theta) = B . k */
          sum[3] += w * n * (b2 + cosine * cosine);
        }
      }
    }
  }
  ray->content = sum[0];
  ray->field = sum[1];
  ray->square = sum[2];
  ray->field_square = sum[3];
  return MYB_OK;
}

void myb_iono_orders_of(const myb_iono_ray* ray, double e, double m_e, double eps0,
                        myb_iono_orders* orders) {
  MYB_CHECK(ray != NULL && orders != NULL);
  const double pi = 0.5 * TWO_PI;
  const double pi2 = pi * pi;
  const double e2 = e * e;
  orders->s1 = e2 / (8.0 * pi2 * m_e * eps0) * ray->content;
  orders->s2 = e2 * e / (16.0 * pi2 * pi * m_e * m_e * eps0) * ray->field;
  orders->s3 = e2 * e2 / (128.0 * pi2 * pi2 * m_e * m_e * eps0 * eps0) * ray->square +
               e2 * e2 / (64.0 * pi2 * pi2 * m_e * m_e * m_e * eps0) * ray->field_square;
}

double myb_iono_phase_delay(const myb_iono_orders* orders, double frequency) {
  MYB_CHECK(orders != NULL);
  if (!(frequency > 0.0) || !isfinite(frequency)) {
    return (double)NAN;
  }
  const double f2 = frequency * frequency;
  return -orders->s1 / f2 - orders->s2 / (f2 * frequency) - orders->s3 / (f2 * f2);
}

double myb_iono_code_delay(const myb_iono_orders* orders, double frequency) {
  MYB_CHECK(orders != NULL);
  if (!(frequency > 0.0) || !isfinite(frequency)) {
    return (double)NAN;
  }
  const double f2 = frequency * frequency;
  return orders->s1 / f2 + 2.0 * orders->s2 / (f2 * frequency) + 3.0 * orders->s3 / (f2 * f2);
}

double myb_iono_bending(double content, double elevation) {
  if (!(content >= 0.0) || !isfinite(content) ||
      !(elevation >= 0.0 && elevation <= 1.5707963267948966)) {
    return (double)NAN;
  }
  const double c = myb_cos(elevation);
  /* mm MHz^4 with the content in TECU: 1e-3 m 1e24 Hz^4 (1e-16)^2 */
  return 2.495e8 * 1e-11 * (1.0 / sqrt(1.0 - 0.8592 * c * c) - 1.0) * content * content;
}
