/* SPDX-License-Identifier: BSD-3-Clause */
/* The motion of a point of the crust (myb/site.h): the tide of the solid body, the pole
 * tide and the loads, after the IERS Conventions 2010, chapter 7. Written from the equations
 * of the chapter; the numbers are those of the kernel. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "myb/io.h"
#include "myb/site.h"
#include "num/poisson.h"
#include "platform/platform.h"

enum { ARGUMENTS = 6, BANDS = 3, DIURNAL_SERIES = 4, ZONAL_SERIES = 2, ROW_MAX = 16 };

static const double seconds_per_century = 3155760000.0;
static const double seconds_per_year = 31557600.0;
static const double arcseconds_per_radian = 206264.80624709636;
static const double two_pi = 6.283185307179586;
static const double degree = 0.017453292519943295;
static const double zonal_norm = 0.63078313050504; /* sqrt(5 / (4 pi)) */

typedef struct wave {
  int8_t multiplier[ARGUMENTS]; /* chi, l, l', F, D, Omega */
  int8_t band;                  /* 0: long-period, 1: diurnal, 2: semidiurnal */
  double amplitude;             /* of the wave in the potential */
  double frequency;             /* cycles per day */
} wave;

struct myb_site_tides {
  myb_allocator allocator;
  myb_internal_poisson_group diurnal;
  int has_diurnal;
  myb_internal_poisson_group zonal;
  int has_zonal;
  double radius;
  double love[6];      /* h(0), h(2), l(0), l(2), h3, l3 */
  double l1[2];        /* diurnal, semidiurnal */
  double imaginary[4]; /* hI, lI diurnal; hI, lI semidiurnal */
  double permanent;    /* H0 */
  double pole[3];      /* radial, south, east: mm per arcsecond */
  double secular_pole[4];
  double ocean_pole;   /* K, m per radian */
  double gamma[2];
  double shift[BANDS]; /* rad */
  wave* waves;
  size_t wave_count;
  uint16_t given[MYB_SITE_MAX_GIVEN_WAVES];
  size_t given_count;
};

struct myb_ocean_loading {
  myb_allocator allocator;
  size_t count;   /* the waves of the model it was made with */
  double* cosine; /* 3 x count: up, west, south */
  double* sine;
};

struct myb_site_grid {
  myb_allocator allocator;
  size_t longitudes;
  size_t latitudes;
  size_t values;
  double first_longitude; /* degrees */
  double longitude_step;
  double first_latitude;
  double latitude_step;
  int round;              /* the longitudes go round the body */
  double* data;           /* [latitude][longitude][value] */
};

static int finite(double x) {
  return x - x == 0.0;
}

static myb_status bad(myb_diag* diag, const char* prefix, const char* key, const char* what) {
  myb_internal_diag(diag, 0, "%s_%s: %s", prefix, key, what);
  return MYB_E_FORMAT;
}

static myb_status numbers(const myb_kpool* pool, const char* prefix, const char* key,
                          const double** v, size_t* n, myb_diag* diag) {
  char name[64];
  snprintf(name, sizeof name, "%s_%s", prefix, key);
  const myb_status status = myb_kpool_numbers(pool, name, v, n);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "%s: %s", name,
                      status == MYB_E_NOT_FOUND ? "not defined" : "not numeric");
    return status;
  }
  for (size_t i = 0; i < *n; ++i) {
    if (!finite((*v)[i])) {
      return bad(diag, prefix, key, "non-finite value");
    }
  }
  return MYB_OK;
}

static myb_status fixed(const myb_kpool* pool, const char* prefix, const char* key, size_t count,
                        double* out, myb_diag* diag) {
  const double* v = NULL;
  size_t n = 0;
  const myb_status status = numbers(pool, prefix, key, &v, &n, diag);
  if (status != MYB_OK) {
    return status;
  }
  if (n != count) {
    return bad(diag, prefix, key, "wrong number of values");
  }
  memcpy(out, v, count * sizeof(double));
  return MYB_OK;
}

/* Doodson's multipliers (tau, s, h, p, N', ps) as those of (chi, l, l', F, D, Omega):
 * tau = chi - s, s = F + Omega, h = s - D, p = s - l, N' = -Omega, ps = s - D - l'. Returns 0
 * for numbers that are not small integers. */
static int doodson(const double* row, int8_t out[ARGUMENTS]) {
  int d[ARGUMENTS];
  for (int i = 0; i < ARGUMENTS; ++i) {
    if (!(fabs(row[i]) <= 20.0) || row[i] != floor(row[i])) {
      return 0;
    }
    d[i] = (int)row[i];
  }
  const int k = d[1] - d[0] + d[2] + d[3] + d[5];
  const int m[ARGUMENTS] = {d[0], -d[3], -d[5], k, -(d[2] + d[5]), k - d[4]};
  for (int i = 0; i < ARGUMENTS; ++i) {
    out[i] = (int8_t)m[i]; /* |m| <= 100 */
  }
  return 1;
}

static myb_status load_waves(myb_site_tides* t, const myb_kpool* pool, const char* prefix,
                             myb_diag* diag) {
  const double* v = NULL;
  size_t n = 0;
  myb_status status = numbers(pool, prefix, "WAVES", &v, &n, diag);
  if (status != MYB_OK) {
    return status;
  }
  if (n == 0 || n % 7 != 0 || n / 7 > 60000) {
    return bad(diag, prefix, "WAVES", "expected rows of 7 numbers (tau s h p N' ps amplitude)");
  }
  void* block = NULL;
  status = myb_internal_alloc_array(&t->allocator, n / 7, sizeof(wave), &block);
  if (status != MYB_OK) {
    return status;
  }
  t->waves = block;
  t->wave_count = n / 7;
  /* the rates of the arguments, cycles per day */
  double rate[ARGUMENTS];
  for (int k = 0; k < ARGUMENTS; ++k) {
    rate[k] = t->diurnal.argument_polynomial[k][1] * t->diurnal.argument_scale[k] / two_pi /
              36525.0;
  }
  for (size_t i = 0; i < t->wave_count; ++i) {
    const double* row = v + 7 * i;
    wave* w = &t->waves[i];
    if (!doodson(row, w->multiplier) || row[0] < 0.0 || row[0] >= (double)BANDS ||
        row[6] == 0.0) {
      return bad(diag, prefix, "WAVES", "a row that is not a wave of the three bands");
    }
    w->band = (int8_t)row[0];
    w->amplitude = row[6];
    w->frequency = 0.0;
    for (int k = 0; k < ARGUMENTS; ++k) {
      w->frequency += (double)w->multiplier[k] * rate[k];
    }
  }
  status = numbers(pool, prefix, "BLQ", &v, &n, diag);
  if (status != MYB_OK) {
    return status;
  }
  if (n == 0 || n % 6 != 0 || n / 6 > (size_t)MYB_SITE_MAX_GIVEN_WAVES) {
    return bad(diag, prefix, "BLQ", "expected 1 .. 16 rows of 6 multipliers");
  }
  t->given_count = n / 6;
  for (size_t g = 0; g < t->given_count; ++g) {
    int8_t m[ARGUMENTS];
    size_t found = t->wave_count;
    if (doodson(v + 6 * g, m)) {
      for (size_t i = 0; i < t->wave_count && found == t->wave_count; ++i) {
        if (memcmp(m, t->waves[i].multiplier, sizeof m) == 0) {
          found = i;
        }
      }
    }
    if (found == t->wave_count) {
      return bad(diag, prefix, "BLQ", "a wave that is not among WAVES");
    }
    for (size_t e = 0; e < g; ++e) {
      const wave* other = &t->waves[t->given[e]];
      if (other->band == t->waves[found].band &&
          other->frequency == t->waves[found].frequency) {
        return bad(diag, prefix, "BLQ", "two waves of one frequency");
      }
    }
    t->given[g] = (uint16_t)found;
  }
  return MYB_OK;
}

void myb_site_tides_destroy(myb_site_tides* tides) {
  if (tides == NULL) {
    return;
  }
  const myb_allocator a = tides->allocator;
  if (tides->has_diurnal) {
    myb_internal_poisson_free(&tides->diurnal, &a);
  }
  if (tides->has_zonal) {
    myb_internal_poisson_free(&tides->zonal, &a);
  }
  myb_internal_release(&a, tides->waves);
  myb_internal_release(&a, tides);
}

static myb_status load_series(const myb_kpool* pool, const char* prefix,
                              const char* const* names, int count, const myb_allocator* a,
                              myb_internal_poisson_group* group, myb_diag* diag) {
  char full[DIURNAL_SERIES][64];
  const char* list[DIURNAL_SERIES];
  for (int k = 0; k < count; ++k) {
    snprintf(full[k], sizeof full[k], "%s_%s", prefix, names[k]);
    list[k] = full[k];
  }
  myb_status status = myb_internal_poisson_load(pool, prefix, list, count, a, group, diag);
  if (status == MYB_OK && group->arguments != ARGUMENTS) {
    myb_internal_poisson_free(group, a);
    status = bad(diag, prefix, "ARGUMENT_COUNT", "expected 6 (chi, l, l', F, D, Omega)");
  }
  return status;
}

myb_status myb_site_tides_create(const myb_kpool* pool, const char* prefix,
                                 const myb_allocator* allocator, myb_site_tides** tides,
                                 myb_diag* diag) {
  MYB_CHECK(pool != NULL && prefix != NULL && tides != NULL);
  myb_internal_diag_clear(diag);
  if (strlen(prefix) > 20) {
    myb_internal_diag(diag, 0, "prefix longer than 20 characters");
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_site_tides), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_site_tides* t = block;
  memset(t, 0, sizeof *t);
  t->allocator = a;
  static const char* const diurnal[DIURNAL_SERIES] = {"BAND1_RS", "BAND1_RC", "BAND1_TS",
                                                      "BAND1_TC"};
  static const char* const zonal[ZONAL_SERIES] = {"BAND0_R", "BAND0_T"};
  status = load_series(pool, prefix, diurnal, DIURNAL_SERIES, &a, &t->diurnal, diag);
  if (status == MYB_OK) {
    t->has_diurnal = 1;
    status = load_series(pool, prefix, zonal, ZONAL_SERIES, &a, &t->zonal, diag);
  }
  if (status == MYB_OK) {
    t->has_zonal = 1;
    status = fixed(pool, prefix, "RADIUS", 1, &t->radius, diag);
  }
  if (status == MYB_OK && !(t->radius > 0.0)) {
    status = bad(diag, prefix, "RADIUS", "not positive");
  }
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "LOVE", 6, t->love, diag);
  }
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "L1", 2, t->l1, diag);
  }
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "IMAGINARY", 4, t->imaginary, diag);
  }
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "PERMANENT", 1, &t->permanent, diag);
  }
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "POLE", 3, t->pole, diag);
  }
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "SECULAR_POLE", 4, t->secular_pole, diag);
  }
  double c[8] = {0};
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "OCEAN_POLE", 8, c, diag);
  }
  if (status == MYB_OK) {
    /* G, aE, rho_w, ge, Omega, GM, gamma: K = 4 pi G aE rho_w Hp / (3 ge),
     * Hp = sqrt(8 pi / 15) Omega^2 aE^4 / GM */
    if (!(c[3] > 0.0) || !(c[5] > 0.0)) {
      status = bad(diag, prefix, "OCEAN_POLE", "ge and GM must be positive");
    } else {
      const double hp = sqrt(4.0 * two_pi / 15.0) * c[4] * c[4] * c[1] * c[1] * c[1] * c[1] /
                        c[5];
      t->ocean_pole = 2.0 * two_pi * c[0] * c[1] * c[2] * hp / (3.0 * c[3]);
      t->gamma[0] = c[6];
      t->gamma[1] = c[7];
    }
  }
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "WAVE_SHIFT", BANDS, t->shift, diag);
  }
  if (status == MYB_OK) {
    for (int b = 0; b < BANDS; ++b) {
      t->shift[b] *= degree;
    }
    status = load_waves(t, pool, prefix, diag);
  }
  if (status != MYB_OK) {
    myb_site_tides_destroy(t);
    return status;
  }
  *tides = t;
  return MYB_OK;
}

size_t myb_site_tides_waves(const myb_site_tides* tides) {
  MYB_CHECK(tides != NULL);
  return tides->wave_count;
}

size_t myb_site_tides_given_waves(const myb_site_tides* tides) {
  MYB_CHECK(tides != NULL);
  return tides->given_count;
}

/* ---- The load of the ocean tides ------------------------------------------------------------- */

/* The second derivatives s of the cubic through (x, u), n points with x increasing, whose
 * end slopes are those of the parabolas through the three end points; zeros (straight lines)
 * for three points or fewer. */
static void spline(const double* x, const double* u, int n, double* s) {
  for (int i = 0; i < n; ++i) {
    s[i] = 0.0;
  }
  if (n <= 3) {
    return;
  }
  /* the slope at the first of three points of the parabola through them */
  const double a1 = x[1] - x[0];
  const double a2 = x[2] - x[0];
  const double first =
      ((u[1] - u[0]) / (a1 * a1) - (u[2] - u[0]) / (a2 * a2)) / (1.0 / a1 - 1.0 / a2);
  const double b1 = x[n - 2] - x[n - 1];
  const double b2 = x[n - 3] - x[n - 1];
  const double last = ((u[n - 2] - u[n - 1]) / (b1 * b1) - (u[n - 3] - u[n - 1]) / (b2 * b2)) /
                      (1.0 / b1 - 1.0 / b2);
  /* h[i-1] s[i-1] + 2 (h[i-1] + h[i]) s[i] + h[i] s[i+1] = 6 (slope[i] - slope[i-1]), with the
   * end slopes in the first and the last rows: three diagonals, eliminated downwards */
  double diagonal[ROW_MAX];
  double right[ROW_MAX];
  diagonal[0] = 2.0 * (x[1] - x[0]);
  right[0] = 6.0 * ((u[1] - u[0]) / (x[1] - x[0]) - first);
  for (int i = 1; i < n; ++i) {
    const double below = x[i] - x[i - 1];
    const double factor = below / diagonal[i - 1];
    if (i < n - 1) {
      const double above = x[i + 1] - x[i];
      diagonal[i] = 2.0 * (below + above) - factor * below;
      right[i] = 6.0 * ((u[i + 1] - u[i]) / above - (u[i] - u[i - 1]) / below) -
                 factor * right[i - 1];
    } else {
      diagonal[i] = 2.0 * below - factor * below;
      right[i] = 6.0 * (last - (u[i] - u[i - 1]) / below) - factor * right[i - 1];
    }
  }
  s[n - 1] = right[n - 1] / diagonal[n - 1];
  for (int i = n - 2; i >= 0; --i) {
    s[i] = (right[i] - (x[i + 1] - x[i]) * s[i + 1]) / diagonal[i];
  }
}

/* The value at y: the end values outside the points. */
static double spline_at(const double* x, const double* u, const double* s, int n, double y) {
  if (y <= x[0]) {
    return u[0];
  }
  if (y >= x[n - 1]) {
    return u[n - 1];
  }
  int k = 1;
  while (k < n - 1 && x[k] < y) { /* bound: the points */
    ++k;
  }
  const double h = x[k] - x[k - 1];
  const double right = x[k] - y;
  const double left = y - x[k - 1];
  return (s[k - 1] * right * right * right + s[k] * left * left * left) / (6.0 * h) +
         left * (u[k] / h - s[k] * h / 6.0) + right * (u[k - 1] / h - s[k - 1] * h / 6.0);
}

/* One component: the amplitudes and lags of the given waves to the cosine and sine parts of
 * every wave. */
static void component(const myb_site_tides* t, const double* amplitude, const double* phase,
                      double* cosine, double* sine) {
  for (int band = 0; band < BANDS; ++band) {
    double x[ROW_MAX];
    double re[ROW_MAX];
    double im[ROW_MAX];
    int n = 0;
    for (size_t g = 0; g < t->given_count; ++g) {
      const wave* w = &t->waves[t->given[g]];
      if (w->band != band) {
        continue;
      }
      /* in its place among those before it, by frequency */
      int at = n;
      while (at > 0 && x[at - 1] > w->frequency) { /* bound: the given waves */
        x[at] = x[at - 1];
        re[at] = re[at - 1];
        im[at] = im[at - 1];
        --at;
      }
      const double scale = amplitude[g] / fabs(w->amplitude);
      x[at] = w->frequency;
      re[at] = scale * myb_cos(phase[g] * degree);
      im[at] = -scale * myb_sin(phase[g] * degree);
      ++n;
    }
    double sr[ROW_MAX];
    double si[ROW_MAX];
    spline(x, re, n, sr);
    spline(x, im, n, si);
    const double cs = myb_cos(t->shift[band]);
    const double sn = myb_sin(t->shift[band]);
    for (size_t i = 0; i < t->wave_count; ++i) {
      const wave* w = &t->waves[i];
      if (w->band != band) {
        continue;
      }
      if (n == 0) {
        cosine[i] = 0.0;
        sine[i] = 0.0;
        continue;
      }
      const double zr = spline_at(x, re, sr, n, w->frequency);
      const double zi = spline_at(x, im, si, n, w->frequency);
      /* H Re(Z exp(i (theta + shift))) = H (Wr cos(theta) - Wi sin(theta)) */
      cosine[i] = w->amplitude * (zr * cs - zi * sn);
      sine[i] = -w->amplitude * (zr * sn + zi * cs);
    }
  }
}

myb_status myb_ocean_loading_create(const myb_site_tides* tides, const double* amplitudes,
                                    const double* phases, const myb_allocator* allocator,
                                    myb_ocean_loading** loading) {
  MYB_CHECK(tides != NULL && amplitudes != NULL && phases != NULL && loading != NULL);
  const size_t given = tides->given_count;
  for (size_t i = 0; i < 3 * given; ++i) {
    if (!(amplitudes[i] >= 0.0) || !finite(amplitudes[i]) || !finite(phases[i])) {
      return MYB_E_DOMAIN;
    }
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_ocean_loading), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_ocean_loading* l = block;
  memset(l, 0, sizeof *l);
  l->allocator = a;
  l->count = tides->wave_count;
  status = myb_internal_alloc_array(&a, 6 * l->count, sizeof(double), &block);
  if (status != MYB_OK) {
    myb_internal_release(&a, l);
    return status;
  }
  l->cosine = block;
  l->sine = l->cosine + 3 * l->count;
  for (size_t c = 0; c < 3; ++c) {
    component(tides, amplitudes + c * given, phases + c * given, l->cosine + c * l->count,
              l->sine + c * l->count);
  }
  *loading = l;
  return MYB_OK;
}

void myb_ocean_loading_destroy(myb_ocean_loading* loading) {
  if (loading == NULL) {
    return;
  }
  const myb_allocator a = loading->allocator;
  myb_internal_release(&a, loading->cosine);
  myb_internal_release(&a, loading);
}

/* What the first argument, the angle of the rotation, gains in `dynamical` seconds: the
 * arguments are taken that much later than the rotation. */
static double rotation_lag(const myb_site_tides* tides, double dynamical) {
  return tides->diurnal.argument_polynomial[0][1] * tides->diurnal.argument_scale[0] /
         seconds_per_century * dynamical;
}

void myb_ocean_loading_at(const myb_site_tides* tides, const myb_ocean_loading* loading,
                          myb_time t, double dynamical, double displacement[3]) {
  MYB_CHECK(tides != NULL && loading != NULL && displacement != NULL &&
            loading->count == tides->wave_count);
  double angle[MYB_POISSON_MAX_ARGUMENTS];
  myb_internal_poisson_arguments(&tides->diurnal,
                                 (myb_time_seconds(t) + dynamical) / seconds_per_century, angle);
  angle[0] -= rotation_lag(tides, dynamical);
  double sum[3] = {0.0, 0.0, 0.0};
  const size_t n = loading->count;
  for (size_t i = 0; i < n; ++i) {
    const wave* w = &tides->waves[i];
    double theta = 0.0;
    for (int k = 0; k < ARGUMENTS; ++k) {
      theta += (double)w->multiplier[k] * angle[k];
    }
    const double c = myb_cos(theta);
    const double s = myb_sin(theta);
    for (size_t j = 0; j < 3; ++j) {
      sum[j] += loading->cosine[j * n + i] * c + loading->sine[j * n + i] * s;
    }
  }
  for (int j = 0; j < 3; ++j) {
    displacement[j] = sum[j];
  }
}

/* ---- Text ------------------------------------------------------------------------------------ */

/* The next line of the text: its span without the end of line. Returns 0 at the end. */
static int next_line(const char* text, size_t size, size_t* at, const char** line,
                     size_t* length) {
  if (*at >= size) {
    return 0;
  }
  size_t end = *at;
  while (end < size && text[end] != '\n') { /* bound: the characters of the text */
    ++end;
  }
  *line = text + *at;
  *length = end - *at;
  *at = end + 1;
  while (*length > 0 && ((*line)[*length - 1] == '\r' || (*line)[*length - 1] == ' ' ||
                         (*line)[*length - 1] == '\t')) {
    --*length;
  }
  while (*length > 0 && ((*line)[0] == ' ' || (*line)[0] == '\t')) {
    ++*line;
    --*length;
  }
  return 1;
}

/* The numbers of a line separated by blanks: their count, or -1 for a field that is not a
 * number or more than `capacity` of them. */
static int line_numbers(const char* line, size_t length, double* out, int capacity) {
  int n = 0;
  size_t at = 0;
  while (at < length) { /* bound: the characters of the line */
    if (line[at] == ' ' || line[at] == '\t') {
      ++at;
      continue;
    }
    size_t end = at;
    while (end < length && line[end] != ' ' && line[end] != '\t') {
      ++end;
    }
    size_t used = 0;
    double value = 0.0;
    if (n == capacity || myb_parse_double(line + at, end - at, &value, &used) != MYB_OK ||
        used != end - at || !finite(value)) {
      return -1;
    }
    out[n++] = value;
    at = end;
  }
  return n;
}

static char lower(char c) {
  return c >= 'A' && c <= 'Z' ? (char)(c - 'A' + 'a') : c;
}

/* Whether the span equals the word, without regard to case. */
static int same_word(const char* span, size_t length, const char* word) {
  if (strlen(word) != length) {
    return 0;
  }
  for (size_t i = 0; i < length; ++i) {
    if (lower(span[i]) != lower(word[i])) {
      return 0;
    }
  }
  return 1;
}

myb_status myb_ocean_loading_load(const myb_site_tides* tides, const char* text, size_t size,
                                  const char* site, const myb_allocator* allocator,
                                  myb_ocean_loading** loading, myb_diag* diag) {
  MYB_CHECK(tides != NULL && (text != NULL || size == 0) && site != NULL && loading != NULL);
  myb_internal_diag_clear(diag);
  const int given = (int)tides->given_count;
  double table[6 * MYB_SITE_MAX_GIVEN_WAVES];
  int rows = -1; /* -1: looking for the name; then the rows read */
  size_t at = 0;
  long number = 0;
  const char* line = NULL;
  size_t length = 0;
  while (rows < 6 && next_line(text, size, &at, &line, &length)) { /* bound: the lines */
    ++number;
    if (length == 0 || (length >= 2 && line[0] == '$' && line[1] == '$')) {
      continue;
    }
    if (rows < 0) {
      if (same_word(line, length, site)) {
        rows = 0;
      }
      continue;
    }
    if (line_numbers(line, length, table + rows * given, given) != given) {
      myb_internal_diag(diag, number, "expected %d numbers of the point '%s'", given, site);
      return MYB_E_FORMAT;
    }
    ++rows;
  }
  if (rows < 0) {
    myb_internal_diag(diag, 0, "no point '%s'", site);
    return MYB_E_NOT_FOUND;
  }
  if (rows < 6) {
    myb_internal_diag(diag, number, "the record of '%s' ends after %d of 6 rows", site, rows);
    return MYB_E_FORMAT;
  }
  const myb_status status =
      myb_ocean_loading_create(tides, table, table + 3 * given, allocator, loading);
  if (status == MYB_E_DOMAIN) {
    myb_internal_diag(diag, number, "a negative amplitude of the point '%s'", site);
    return MYB_E_FORMAT;
  }
  return status;
}

myb_status myb_ocean_loading_load_file(const myb_site_tides* tides, const char* path,
                                       const char* site, const myb_allocator* allocator,
                                       myb_ocean_loading** loading, myb_diag* diag) {
  MYB_CHECK(path != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &a, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded = myb_ocean_loading_load(tides, (const char*)map.data, map.size, site,
                                                   allocator, loading, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

/* ---- The load of the atmospheric tides ------------------------------------------------------ */

myb_status myb_site_atmosphere_load(const char* text, size_t size, const char* site,
                                    double coefficients[12], myb_diag* diag) {
  MYB_CHECK((text != NULL || size == 0) && site != NULL && coefficients != NULL);
  myb_internal_diag_clear(diag);
  double rows[3][4];
  int row = -1; /* -1: looking for the name; then the rows read */
  size_t at = 0;
  long number = 0;
  const char* line = NULL;
  size_t length = 0;
  while (row < 3 && next_line(text, size, &at, &line, &length)) { /* bound: the lines */
    ++number;
    if (length == 0 || (length >= 2 && line[0] == '$' && line[1] == '$')) {
      continue;
    }
    if (row < 0) {
      size_t word = 0;
      while (word < length && line[word] != ' ' && line[word] != '\t') { /* bound: the line */
        ++word;
      }
      if (!same_word(line, word, site)) {
        continue;
      }
      if (word == length) {
        row = 0; /* the name on its line: three rows follow */
        continue;
      }
      /* the name and 18 numbers: radial, east, north of S1, S2, S3 */
      double numbers[18];
      if (line_numbers(line + word, length - word, numbers, 18) != 18) {
        myb_internal_diag(diag, number, "expected the name and 18 numbers of the point '%s'",
                          site);
        return MYB_E_FORMAT;
      }
      for (int c = 0; c < 3; ++c) {
        for (int k = 0; k < 4; ++k) {
          coefficients[4 * c + k] = 1e-3 * numbers[6 * c + k];
        }
      }
      return MYB_OK;
    }
    if (line_numbers(line, length, rows[row], 4) != 4) {
      myb_internal_diag(diag, number, "expected 4 numbers of the point '%s'", site);
      return MYB_E_FORMAT;
    }
    ++row;
  }
  if (row < 0) {
    myb_internal_diag(diag, 0, "no point '%s'", site);
    return MYB_E_NOT_FOUND;
  }
  if (row < 3) {
    myb_internal_diag(diag, number, "the record of '%s' ends after %d of 3 rows", site, row);
    return MYB_E_FORMAT;
  }
  /* radial, north-south, east-west -> up, east, north; millimetres -> metres */
  static const int from[3] = {0, 2, 1};
  for (int c = 0; c < 3; ++c) {
    for (int k = 0; k < 4; ++k) {
      coefficients[4 * c + k] = 1e-3 * rows[from[c]][k];
    }
  }
  return MYB_OK;
}

myb_status myb_site_atmosphere_load_file(const char* path, const char* site,
                                         double coefficients[12], myb_diag* diag) {
  MYB_CHECK(path != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(NULL);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &a, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded =
      myb_site_atmosphere_load((const char*)map.data, map.size, site, coefficients, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

/* ---- A grid of coefficients ------------------------------------------------------------------ */

void myb_site_grid_destroy(myb_site_grid* grid) {
  if (grid == NULL) {
    return;
  }
  const myb_allocator a = grid->allocator;
  myb_internal_release(&a, grid->data);
  myb_internal_release(&a, grid);
}

static myb_status grid_fail(myb_site_grid* grid, myb_diag* diag, long line, const char* what) {
  myb_internal_diag(diag, line, "%s", what);
  myb_site_grid_destroy(grid);
  return MYB_E_FORMAT;
}

myb_status myb_site_grid_load(const char* text, size_t size, const myb_allocator* allocator,
                              myb_site_grid** grid, myb_diag* diag) {
  MYB_CHECK((text != NULL || size == 0) && grid != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_site_grid), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_site_grid* g = block;
  memset(g, 0, sizeof *g);
  g->allocator = a;
  /* the header: name = value */
  static const char* const keys[6] = {"Number_longitude_Grid_Points", "First_longitude_degrees",
                                      "Longitude_step_degrees",       "Number_latitude_grid_points",
                                      "First_latitude_degrees",       "Latitude_step_degrees"};
  double header[6] = {0};
  unsigned seen = 0;
  size_t at = 0;
  long number = 0;
  const char* line = NULL;
  size_t length = 0;
  int dashes = 0;
  while (!dashes && next_line(text, size, &at, &line, &length)) { /* bound: the lines */
    ++number;
    if (length >= 3 && line[0] == '-' && line[1] == '-' && line[2] == '-') {
      dashes = 1;
      continue;
    }
    size_t equals = 0;
    while (equals < length && line[equals] != '=') {
      ++equals;
    }
    if (equals == length) {
      continue;
    }
    size_t name = equals;
    while (name > 0 && (line[name - 1] == ' ' || line[name - 1] == '\t')) {
      --name;
    }
    for (int k = 0; k < 6; ++k) {
      if (same_word(line, name, keys[k])) {
        double value = 0.0;
        if (line_numbers(line + equals + 1, length - equals - 1, &value, 1) != 1) {
          return grid_fail(g, diag, number, "the value of a header line is not a number");
        }
        header[k] = value;
        seen |= 1u << k;
      }
    }
  }
  if (!dashes || seen != 63u) {
    return grid_fail(g, diag, number, "the header does not give the grid");
  }
  if (!(header[0] >= 2.0) || !(header[3] >= 2.0) || header[0] > 100000.0 ||
      header[3] > 100000.0 || header[0] != floor(header[0]) || header[3] != floor(header[3]) ||
      !(header[2] > 0.0) || !(header[5] > 0.0)) {
    return grid_fail(g, diag, number, "the counts and steps of the grid");
  }
  g->longitudes = (size_t)header[0];
  g->first_longitude = header[1];
  g->longitude_step = header[2];
  g->latitudes = (size_t)header[3];
  g->first_latitude = header[4];
  g->latitude_step = header[5];
  g->round = fabs(header[0] * header[2] - 360.0) < 1e-6 * header[2];
  const size_t nodes = g->longitudes * g->latitudes;
  double row[2 + MYB_SITE_GRID_MAX_VALUES];
  for (size_t node = 0; node < nodes; ++node) {
    int n = 0;
    do { /* bound: the lines */
      if (!next_line(text, size, &at, &line, &length)) {
        return grid_fail(g, diag, number, "fewer rows than nodes");
      }
      ++number;
    } while (length == 0);
    n = line_numbers(line, length, row, 2 + MYB_SITE_GRID_MAX_VALUES);
    if (n < 3 || (node > 0 && (size_t)(n - 2) != g->values)) {
      return grid_fail(g, diag, number, "a row of the grid");
    }
    if (node == 0) {
      g->values = (size_t)(n - 2);
      status = myb_internal_alloc_array(&a, nodes * g->values, sizeof(double), &block);
      if (status != MYB_OK) {
        myb_site_grid_destroy(g);
        return status;
      }
      g->data = block;
    }
    const double longitude =
        g->first_longitude + (double)(node % g->longitudes) * g->longitude_step;
    const double latitude = g->first_latitude + (double)(node / g->longitudes) * g->latitude_step;
    if (fabs(row[0] - longitude) > 1e-6 || fabs(row[1] - latitude) > 1e-6) {
      return grid_fail(g, diag, number, "a row that is not the next node of the grid");
    }
    memcpy(g->data + node * g->values, row + 2, g->values * sizeof(double));
  }
  *grid = g;
  return MYB_OK;
}

myb_status myb_site_grid_load_file(const char* path, const myb_allocator* allocator,
                                   myb_site_grid** grid, myb_diag* diag) {
  MYB_CHECK(path != NULL && grid != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &a, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded =
      myb_site_grid_load((const char*)map.data, map.size, allocator, grid, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

/* ---- The grid of the atmospheric tides -------------------------------------------------------- */

/* latitude, longitude, then (radial, east, north) x (S1, S2, S3) x (cos, sin) */
enum { TIDAL_COLUMNS = 20 };

static myb_status tidal_rows(const char* text, size_t size, myb_site_grid* g, myb_diag* diag) {
  size_t at = 0;
  long number = 0;
  const char* line = NULL;
  size_t length = 0;
  double row[TIDAL_COLUMNS];
  size_t node = 0;
  const size_t nodes = g->longitudes * g->latitudes;
  while (node < nodes && next_line(text, size, &at, &line, &length)) { /* bound: the lines */
    ++number;
    if (length == 0) {
      continue;
    }
    if (line_numbers(line, length, row, TIDAL_COLUMNS) != TIDAL_COLUMNS) {
      myb_internal_diag(diag, number, "expected %d numbers", (int)TIDAL_COLUMNS);
      return MYB_E_FORMAT;
    }
    const double latitude = g->first_latitude + (double)(node / g->longitudes) * g->latitude_step;
    const double longitude =
        g->first_longitude + (double)(node % g->longitudes) * g->longitude_step;
    if (fabs(row[0] - latitude) > 1e-6 || fabs(row[1] - longitude) > 1e-6) {
      myb_internal_diag(diag, number, "a row that is not the next node of the grid");
      return MYB_E_FORMAT;
    }
    double* values = g->data + node * g->values;
    for (int c = 0; c < 3; ++c) { /* millimetres -> metres; S3 left */
      for (int k = 0; k < 4; ++k) {
        values[4 * c + k] = 1e-3 * row[2 + 6 * c + k];
      }
    }
    ++node;
  }
  return MYB_OK;
}

myb_status myb_site_atmosphere_grid_load(const char* text, size_t size,
                                         const myb_allocator* allocator, myb_site_grid** grid,
                                         myb_diag* diag) {
  MYB_CHECK((text != NULL || size == 0) && grid != NULL);
  myb_internal_diag_clear(diag);
  /* the shape: the first two rows give the longitude step, the first row of another
   * latitude the count of longitudes and the latitude step */
  size_t at = 0;
  long number = 0;
  const char* line = NULL;
  size_t length = 0;
  double row[TIDAL_COLUMNS];
  double first[2] = {0.0, 0.0};
  double second_longitude = 0.0;
  double second_latitude = 0.0;
  size_t rows = 0;
  size_t longitudes = 0;
  while (next_line(text, size, &at, &line, &length)) { /* bound: the lines */
    ++number;
    if (length == 0) {
      continue;
    }
    if (line_numbers(line, length, row, TIDAL_COLUMNS) != TIDAL_COLUMNS) {
      myb_internal_diag(diag, number, "expected %d numbers", (int)TIDAL_COLUMNS);
      return MYB_E_FORMAT;
    }
    if (rows == 0) {
      first[0] = row[0];
      first[1] = row[1];
    } else if (rows == 1) {
      second_longitude = row[1];
    }
    if (longitudes == 0 && rows > 0 && row[0] != first[0]) {
      longitudes = rows;
      second_latitude = row[0];
    }
    ++rows;
  }
  if (longitudes < 2 || rows % longitudes != 0 || rows / longitudes > 100000 ||
      longitudes > 100000 || !(second_longitude > first[1]) ||
      !(second_latitude != first[0])) {
    myb_internal_diag(diag, number, "the rows do not make a grid");
    return MYB_E_FORMAT;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_site_grid), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_site_grid* g = block;
  memset(g, 0, sizeof *g);
  g->allocator = a;
  g->longitudes = longitudes;
  g->latitudes = rows / longitudes;
  g->values = 12;
  g->first_latitude = first[0];
  g->first_longitude = first[1];
  g->longitude_step = second_longitude - first[1];
  g->latitude_step = second_latitude - first[0];
  g->round = fabs((double)longitudes * g->longitude_step - 360.0) < 1e-6 * g->longitude_step;
  status = myb_internal_alloc_array(&a, rows * g->values, sizeof(double), &block);
  if (status != MYB_OK) {
    myb_site_grid_destroy(g);
    return status;
  }
  g->data = block;
  status = tidal_rows(text, size, g, diag);
  if (status != MYB_OK) {
    myb_site_grid_destroy(g);
    return status;
  }
  *grid = g;
  return MYB_OK;
}

myb_status myb_site_atmosphere_grid_load_file(const char* path, const myb_allocator* allocator,
                                              myb_site_grid** grid, myb_diag* diag) {
  MYB_CHECK(path != NULL && grid != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &a, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded =
      myb_site_atmosphere_grid_load((const char*)map.data, map.size, allocator, grid, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

size_t myb_site_grid_values(const myb_site_grid* grid) {
  MYB_CHECK(grid != NULL);
  return grid->values;
}

myb_status myb_site_grid_at(const myb_site_grid* grid, double longitude, double latitude,
                            double* values) {
  MYB_CHECK(grid != NULL && values != NULL);
  if (!finite(longitude) || !finite(latitude)) {
    return MYB_E_DOMAIN;
  }
  const double tiny = 1e-9; /* of a step: a point on the edge is inside */
  double x = (longitude / degree - grid->first_longitude) / grid->longitude_step;
  double y = (latitude / degree - grid->first_latitude) / grid->latitude_step;
  const double columns = (double)grid->longitudes;
  if (grid->round) {
    x -= columns * floor(x / columns);
    if (x >= columns) {
      x = 0.0;
    }
  } else {
    if (x < -tiny || x > columns - 1.0 + tiny) {
      return MYB_E_RANGE;
    }
    x = fmin(fmax(x, 0.0), columns - 1.0);
  }
  const double rows = (double)grid->latitudes;
  if (y < -tiny || y > rows - 1.0 + tiny) {
    return MYB_E_RANGE;
  }
  y = fmin(fmax(y, 0.0), rows - 1.0);
  size_t i0 = (size_t)x;
  size_t j0 = (size_t)y;
  if (!grid->round && i0 > grid->longitudes - 2) {
    i0 = grid->longitudes - 2;
  }
  if (j0 > grid->latitudes - 2) {
    j0 = grid->latitudes - 2;
  }
  const size_t i1 = (i0 + 1) % grid->longitudes;
  const double fx = x - (double)i0;
  const double fy = y - (double)j0;
  const size_t n = grid->values;
  const double* low = grid->data + j0 * grid->longitudes * n;
  const double* high = low + grid->longitudes * n;
  for (size_t k = 0; k < n; ++k) {
    const double south = low[i0 * n + k] + fx * (low[i1 * n + k] - low[i0 * n + k]);
    const double north = high[i0 * n + k] + fx * (high[i1 * n + k] - high[i0 * n + k]);
    values[k] = south + fy * (north - south);
  }
  return MYB_OK;
}

/* ---- The displacement ------------------------------------------------------------------------ */

myb_status myb_site_displacement(const myb_site_tides* tides, const double position[3],
                                 const double* enu, const myb_site_loads* loads,
                                 const myb_site_input* input, double displacement[3]) {
  MYB_CHECK(tides != NULL && position != NULL && input != NULL && displacement != NULL &&
            (input->count == 0 || input->raisers != NULL));
  const int flags = input->flags;
  const double x = position[0];
  const double y = position[1];
  const double z = position[2];
  const double rho = sqrt(x * x + y * y);
  const double r = sqrt(x * x + y * y + z * z);
  if (!finite(r) || !(rho > 0.0) || !(input->gm > 0.0) || !finite(input->gm) ||
      !finite(input->xp) || !finite(input->yp) || !finite(input->t.f) ||
      !finite(input->dynamical)) {
    return MYB_E_DOMAIN;
  }
  const double sp = z / r; /* of the geocentric latitude and of the longitude */
  const double cp = rho / r;
  const double cl = x / rho;
  const double sl = y / rho;
  const double s2p = 2.0 * sp * cp;
  const double c2p = cp * cp - sp * sp;
  const double s2l = 2.0 * sl * cl;
  const double c2l = cl * cl - sl * sl;
  const double p2 = 1.5 * sp * sp - 0.5;
  const double up[3] = {cp * cl, cp * sl, sp};
  const double east[3] = {-sl, cl, 0.0};
  const double north[3] = {-sp * cl, -sp * sl, cp};
  const double h2 = tides->love[0] + tides->love[1] * p2;
  const double l2 = tides->love[2] + tides->love[3] * p2;
  double d[3] = {0.0, 0.0, 0.0}; /* the part along the raisers and the loads */
  double dr = 0.0;               /* along the radius, to the north, to the east */
  double dn = 0.0;
  double de = 0.0;
  if (flags & MYB_SITE_SOLID) {
    const double h3 = tides->love[4];
    const double l3 = tides->love[5];
    for (size_t j = 0; j < input->count; ++j) {
      const double* b = input->raisers[j].r;
      const double gm = input->raisers[j].gm;
      const double b2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
      const double dist = sqrt(b2);
      if (!(dist > 0.0) || !finite(dist) || !(gm >= 0.0) || !finite(gm)) {
        return MYB_E_DOMAIN;
      }
      const double ratio = tides->radius / dist;
      const double f2 = gm / input->gm * tides->radius * ratio * ratio * ratio;
      const double f3 = f2 * ratio;
      const double c = (b[0] * up[0] + b[1] * up[1] + b[2] * up[2]) / dist;
      /* eq. 5 and 6: along the radius and along the direction of the raiser */
      const double towards = f2 * 3.0 * l2 * c + f3 * l3 * (7.5 * c * c - 1.5);
      dr += f2 * (h2 * (1.5 * c * c - 0.5) - 3.0 * l2 * c * c) +
            f3 * (h3 * (2.5 * c * c * c - 1.5 * c) - l3 * (7.5 * c * c - 1.5) * c);
      for (int i = 0; i < 3; ++i) {
        d[i] += towards * b[i] / dist;
      }
      /* sin(2 Phi) sin and cos of (lambda - lambda_j); cos^2(Phi) sin and cos of twice it */
      const double sd = 2.0 * b[2] * (sl * b[0] - cl * b[1]) / b2;
      const double cd = 2.0 * b[2] * (cl * b[0] + sl * b[1]) / b2;
      const double ss = (s2l * (b[0] * b[0] - b[1] * b[1]) - c2l * 2.0 * b[0] * b[1]) / b2;
      const double cs = (c2l * (b[0] * b[0] - b[1] * b[1]) + s2l * 2.0 * b[0] * b[1]) / b2;
      /* the imaginary parts, eq. 10 and 11 */
      dr += f2 * (-0.75 * tides->imaginary[0] * sd * s2p -
                  0.75 * tides->imaginary[2] * ss * cp * cp);
      dn += f2 * (-1.5 * tides->imaginary[1] * sd * c2p + 0.75 * tides->imaginary[3] * ss * s2p);
      de += f2 * (-1.5 * tides->imaginary[1] * cd * sp - 1.5 * tides->imaginary[3] * cs * cp);
      /* l(1), eq. 8 and 9 (P21 = 1.5 sin(2 Phi), P22 = 3 cos^2(Phi)) */
      dn += f2 * (-1.5 * tides->l1[0] * sp * sp * cd - 1.5 * tides->l1[1] * sp * cp * cs);
      de += f2 * (1.5 * tides->l1[0] * sp * c2p * sd - 1.5 * tides->l1[1] * sp * sp * cp * ss);
    }
    if (flags & MYB_SITE_MEAN_TIDE) {
      dr -= zonal_norm * tides->permanent * h2 * p2;
      dn -= zonal_norm * tides->permanent * l2 * 1.5 * s2p;
    }
  }
  const double seconds = myb_time_seconds(input->t);
  if (flags & MYB_SITE_FREQUENCY) {
    /* the series on TT; the rotation, which every diurnal row has once, was that much
     * earlier: the same as a point that much to the west */
    const double centuries = (seconds + input->dynamical) / seconds_per_century;
    const double lag = rotation_lag(tides, input->dynamical);
    const double cg = myb_cos(lag);
    const double sg = myb_sin(lag);
    const double ca = cl * cg + sl * sg;
    const double sa = sl * cg - cl * sg;
    double v[MYB_POISSON_MAX_SERIES];
    double rates[MYB_POISSON_MAX_SERIES];
    myb_internal_poisson_eval(&tides->diurnal, centuries, v, rates);
    dr += s2p * (v[0] * ca + v[1] * sa);
    dn += c2p * (v[2] * ca + v[3] * sa);
    de += sp * (v[3] * ca - v[2] * sa);
    myb_internal_poisson_eval(&tides->zonal, centuries, v, rates);
    dr += p2 * v[0];
    dn += s2p * v[1];
  }
  /* the wobble about the secular pole, arcseconds */
  const double years = seconds / seconds_per_year;
  const double m1 = input->xp * arcseconds_per_radian -
                    (tides->secular_pole[0] + tides->secular_pole[1] * years) * 1e-3;
  const double m2 = -(input->yp * arcseconds_per_radian -
                      (tides->secular_pole[2] + tides->secular_pole[3] * years) * 1e-3);
  if (flags & MYB_SITE_POLE) {
    /* eq. 26 with the colatitude: sin(2 theta) = sin(2 phi), cos(2 theta) = -cos(2 phi) */
    const double along = m1 * cl + m2 * sl;
    dr += 1e-3 * tides->pole[0] * s2p * along;
    dn += 1e-3 * tides->pole[1] * c2p * along; /* south: -cos(2 phi), north: the opposite */
    de += 1e-3 * tides->pole[2] * sp * (m1 * sl - m2 * cl);
  }
  if (loads != NULL && (flags & (MYB_SITE_OCEAN | MYB_SITE_ATMOSPHERE | MYB_SITE_OCEAN_POLE))) {
    double e = 0.0; /* east, north, up */
    double n = 0.0;
    double u = 0.0;
    if ((flags & MYB_SITE_OCEAN) && loads->ocean != NULL) {
      double uws[3];
      myb_ocean_loading_at(tides, loads->ocean, input->t, input->dynamical, uws);
      u += uws[0];
      e -= uws[1];
      n -= uws[2];
    }
    if (flags & MYB_SITE_ATMOSPHERE) {
      /* the day of UT1 from its midnight */
      const int64_t whole = (input->t.s + 43200) % 86400;
      const double day = ((double)(whole < 0 ? whole + 86400 : whole) + input->t.f) / 86400.0;
      const double c1 = myb_cos(two_pi * day);
      const double s1 = myb_sin(two_pi * day);
      const double c2 = c1 * c1 - s1 * s1;
      const double s2 = 2.0 * s1 * c1;
      const double* a = loads->atmosphere;
      u += a[0] * c1 + a[1] * s1 + a[2] * c2 + a[3] * s2;
      e += a[4] * c1 + a[5] * s1 + a[6] * c2 + a[7] * s2;
      n += a[8] * c1 + a[9] * s1 + a[10] * c2 + a[11] * s2;
    }
    if (flags & MYB_SITE_OCEAN_POLE) {
      const double k = tides->ocean_pole / arcseconds_per_radian;
      const double real = k * (m1 * tides->gamma[0] + m2 * tides->gamma[1]);
      const double imag = k * (m2 * tides->gamma[0] - m1 * tides->gamma[1]);
      const double* c = loads->ocean_pole;
      u += real * c[0] + imag * c[1];
      n += real * c[2] + imag * c[3];
      e += real * c[4] + imag * c[5];
    }
    if (e != 0.0 || n != 0.0 || u != 0.0) {
      if (enu == NULL) {
        return MYB_E_DOMAIN;
      }
      for (int i = 0; i < 3; ++i) {
        d[i] += e * enu[i] + n * enu[3 + i] + u * enu[6 + i];
      }
    }
  }
  for (int i = 0; i < 3; ++i) {
    displacement[i] = d[i] + dr * up[i] + dn * north[i] + de * east[i];
    if (!finite(displacement[i])) {
      return MYB_E_DOMAIN;
    }
  }
  return MYB_OK;
}
