/* SPDX-License-Identifier: BSD-3-Clause */
/* Tidal changes of a gravity field after IERS Conventions 2010, chapter 6, with every number
 * from data: solid tides (Love numbers, step 1; frequency-dependent corrections as Poisson
 * series, step 2; the permanent tide, step 3), ocean tides from harmonic files of the form of
 * the IERS fes2004_Cnm-Snm.dat (eq. 6.15), solid and ocean pole tides (sections 6.4, 6.5). */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "grav/field_internal.h"
#include "myb/cmath.h"
#include "myb/grav.h"
#include "myb/io.h"
#include "num/poisson.h"
#include "platform/platform.h"

enum { STEP2_SERIES = 5, ARGUMENTS = 6, OCEAN_MAX_DEGREE = 360 };

static const double seconds_per_century = 3155760000.0;
static const double seconds_per_year = 31557600.0;
static const double arcseconds_per_radian = 206264.80624709636;

typedef struct ocean_wave {
  int8_t multiplier[ARGUMENTS]; /* chi, l, l', F, D, Omega */
  double doodson;               /* as written in the file, for identification */
} ocean_wave;

typedef struct ocean_entry {
  uint32_t wave;
  int16_t n;
  int16_t m;
  double cp; /* C+_f,nm, S+, C-, S- */
  double sp;
  double cm;
  double sm;
} ocean_entry;

/* A secondary wave of the table of admittance: its arguments, those of its two pivot waves
 * and the parts of their coefficients it takes (eq. 6.16). */
typedef struct admittance_row {
  int8_t wave[ARGUMENTS];
  int8_t pivot[2][ARGUMENTS];
  double part[2];
} admittance_row;

/* A secondary wave attached to a pivot wave that is loaded: where the pivot has exp(i theta),
 * the sum has factor exp(i theta_f) more. */
typedef struct ocean_link {
  int8_t multiplier[ARGUMENTS];
  double factor;
} ocean_link;

struct myb_tides {
  myb_allocator allocator;
  myb_internal_poisson_group step2;
  int has_step2;
  double k_re[4][4]; /* [n][m], n = 2, 3 */
  double k_im[4][4];
  double k_plus[3];
  int has_love[4][4];
  double permanent; /* A0 H0, multiplied by Re k20 when used */
  double pole_solid[4];
  double pole_ocean[4];
  double secular_pole[4];
  ocean_wave* waves;
  size_t wave_count;
  ocean_entry* entries;
  size_t entry_count;
  admittance_row* table;
  size_t table_count;
  ocean_link* links;     /* by pivot wave */
  uint32_t* link_start;  /* wave_count + 1; NULL: none */
  size_t secondary;      /* waves attached */
};

static void from_digits(const int d[6], int8_t out[ARGUMENTS]);
static myb_status load_admittance(myb_tides* t, const myb_kpool* pool, const char* prefix,
                                  myb_diag* diag);

static myb_status numbers(const myb_kpool* pool, const char* prefix, const char* key,
                          const double** v, size_t* n, myb_diag* diag) {
  char name[64];
  snprintf(name, sizeof name, "%s_%s", prefix, key);
  const myb_status status = myb_kpool_numbers(pool, name, v, n);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "%s: %s", name,
                      status == MYB_E_NOT_FOUND ? "not defined" : "not numeric");
  }
  return status;
}

static myb_status bad(myb_diag* diag, const char* prefix, const char* key, const char* what) {
  myb_internal_diag(diag, 0, "%s_%s: %s", prefix, key, what);
  return MYB_E_FORMAT;
}

static int finite(double x) {
  return x - x == 0.0;
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
  for (size_t i = 0; i < n; ++i) {
    if (!finite(v[i])) {
      return bad(diag, prefix, key, "non-finite value");
    }
    out[i] = v[i];
  }
  return MYB_OK;
}

void myb_tides_destroy(myb_tides* tides) {
  if (tides == NULL) {
    return;
  }
  const myb_allocator a = tides->allocator;
  if (tides->has_step2) {
    myb_internal_poisson_free(&tides->step2, &a);
  }
  myb_internal_release(&a, tides->waves);
  myb_internal_release(&a, tides->entries);
  myb_internal_release(&a, tides->table);
  myb_internal_release(&a, tides->links);
  myb_internal_release(&a, tides->link_start);
  myb_internal_release(&a, tides);
}

myb_status myb_tides_create(const myb_kpool* pool, const char* prefix,
                            const myb_allocator* allocator, myb_tides** tides, myb_diag* diag) {
  MYB_CHECK(pool != NULL && prefix != NULL && tides != NULL);
  myb_internal_diag_clear(diag);
  if (strlen(prefix) > 20) {
    myb_internal_diag(diag, 0, "prefix longer than 20 characters");
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_tides), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_tides* t = block;
  memset(t, 0, sizeof *t);
  t->allocator = a;
  static const char* const names[STEP2_SERIES] = {"C20", "C21", "S21", "C22", "S22"};
  char full[STEP2_SERIES][64];
  const char* list[STEP2_SERIES];
  for (int k = 0; k < STEP2_SERIES; ++k) {
    snprintf(full[k], sizeof full[k], "%s_%s", prefix, names[k]);
    list[k] = full[k];
  }
  status = myb_internal_poisson_load(pool, prefix, list, STEP2_SERIES, &a, &t->step2, diag);
  if (status == MYB_OK) {
    t->has_step2 = 1;
    if (t->step2.arguments != ARGUMENTS) {
      status = bad(diag, prefix, "ARGUMENT_COUNT", "expected 6 (chi, l, l', F, D, Omega)");
    }
  }
  /* Love numbers */
  const double* v = NULL;
  size_t n = 0;
  if (status == MYB_OK) {
    status = numbers(pool, prefix, "LOVE", &v, &n, diag);
  }
  if (status == MYB_OK && (n == 0 || n % 5 != 0)) {
    status = bad(diag, prefix, "LOVE", "expected rows of 5 numbers (n m Re Im k+)");
  }
  for (size_t r = 0; status == MYB_OK && r < n / 5; ++r) {
    const double* row = v + 5 * r;
    if (!(row[0] == 2.0 || row[0] == 3.0) || !(row[1] >= 0.0 && row[1] <= row[0]) ||
        row[1] != floor(row[1]) || !finite(row[2]) || !finite(row[3]) || !finite(row[4]) ||
        (row[0] == 3.0 && row[4] != 0.0)) {
      status = bad(diag, prefix, "LOVE", "invalid row (n = 2 or 3, 0 <= m <= n; k+ only for n = 2)");
      break;
    }
    const int nn = (int)row[0];
    const int mm = (int)row[1];
    t->k_re[nn][mm] = row[2];
    t->k_im[nn][mm] = row[3];
    t->has_love[nn][mm] = 1;
    if (nn == 2) {
      t->k_plus[mm] = row[4];
    }
  }
  double permanent[2] = {0.0, 0.0};
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "PERMANENT", 2, permanent, diag);
  }
  t->permanent = permanent[0] * permanent[1];
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "POLE_SOLID", 4, t->pole_solid, diag);
  }
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "POLE_OCEAN", 4, t->pole_ocean, diag);
  }
  if (status == MYB_OK) {
    status = fixed(pool, prefix, "SECULAR_POLE", 4, t->secular_pole, diag);
  }
  if (status == MYB_OK) {
    status = load_admittance(t, pool, prefix, diag);
  }
  if (status != MYB_OK) {
    myb_tides_destroy(t);
    return status;
  }
  *tides = t;
  return MYB_OK;
}

/* ---- Ocean tide files ------------------------------------------------------------------------ */

/* A Doodson number written as d1d2d3.d4d5d6 (leading zero possibly omitted): the multipliers of
 * (tau, s, h, p, N', ps) are (d1, d2 - 5, ..., d6 - 5); in the arguments (chi, l, l', F, D,
 * Omega), with tau = chi - s, s = F + Omega, h = s - D, p = s - l, N' = -Omega, ps = s - D - l'. */
static int doodson(const char* s, size_t n, int8_t out[ARGUMENTS]) {
  int digits[6];
  int count = 0;
  int before = 0;
  int dot = 0;
  for (size_t i = 0; i < n; ++i) {
    if (s[i] == '.') {
      if (dot) {
        return 0;
      }
      dot = 1;
      continue;
    }
    if (s[i] < '0' || s[i] > '9' || count == 6) {
      return 0;
    }
    if (!dot) {
      ++before;
    }
    digits[count++] = s[i] - '0';
  }
  if (!dot || count != before + 3 || before < 2 || before > 3) {
    return 0;
  }
  int d[6] = {0, 0, 0, 0, 0, 0};
  const int shift = 3 - before; /* "55.565" stands for 055.565 */
  for (int i = 0; i < count; ++i) {
    d[i + shift] = digits[i];
  }
  from_digits(d, out);
  return 1;
}

/* The same from the number itself. */
static int doodson_number(double code, int8_t out[ARGUMENTS]) {
  const double scaled = floor(code * 1000.0 + 0.5);
  if (!(scaled > 0.0) || !(scaled < 1e6) || fabs(code * 1000.0 - scaled) > 1e-6) {
    return 0;
  }
  int rest = (int)scaled;
  int d[6];
  for (int i = 5; i >= 0; --i) {
    d[i] = rest % 10;
    rest /= 10;
  }
  from_digits(d, out);
  return 1;
}

static void from_digits(const int d[6], int8_t out[ARGUMENTS]) {
  const int n1 = d[0];
  const int n2 = d[1] - 5;
  const int n3 = d[2] - 5;
  const int n4 = d[3] - 5;
  const int n5 = d[4] - 5;
  const int n6 = d[5] - 5;
  out[0] = (int8_t)n1;
  out[1] = (int8_t)(-n4);
  out[2] = (int8_t)(-n6);
  out[3] = (int8_t)((n2 - n1) + n3 + n4 + n6);
  out[4] = (int8_t)(-n3 - n6);
  out[5] = (int8_t)((n2 - n1) + n3 + n4 - n5 + n6);
}

/* The table of admittance, if the kernel has one: <prefix>_ADMITTANCE, rows of (Doodson
 * number, amplitude) of a secondary wave and of its two pivot waves. */
static myb_status load_admittance(myb_tides* t, const myb_kpool* pool, const char* prefix,
                                  myb_diag* diag) {
  char name[64];
  snprintf(name, sizeof name, "%s_ADMITTANCE", prefix);
  const double* v = NULL;
  size_t n = 0;
  const myb_status found = myb_kpool_numbers(pool, name, &v, &n);
  if (found == MYB_E_NOT_FOUND) {
    return MYB_OK;
  }
  if (found != MYB_OK) {
    myb_internal_diag(diag, 0, "%s: not numeric", name);
    return found;
  }
  if (n == 0 || n % 6 != 0) {
    return bad(diag, prefix, "ADMITTANCE", "expected rows of 6 numbers (wave H, pivot H twice)");
  }
  void* block = NULL;
  const myb_status status =
      myb_internal_alloc_array(&t->allocator, n / 6, sizeof(admittance_row), &block);
  if (status != MYB_OK) {
    return status;
  }
  t->table = block;
  t->table_count = n / 6;
  for (size_t r = 0; r < t->table_count; ++r) {
    const double* row = v + 6 * r;
    admittance_row* a = &t->table[r];
    if (!doodson_number(row[0], a->wave) || !doodson_number(row[2], a->pivot[0]) ||
        !doodson_number(row[4], a->pivot[1]) || !finite(row[1]) || !finite(row[3]) ||
        !finite(row[5]) || row[3] == 0.0 || row[5] == 0.0) {
      return bad(diag, prefix, "ADMITTANCE", "a row that is not three waves with amplitudes");
    }
    /* the rates of the three arguments, in any one unit */
    double rate[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < ARGUMENTS; ++k) {
      const double unit = t->step2.argument_polynomial[k][1] * t->step2.argument_scale[k];
      rate[0] += (double)a->wave[k] * unit;
      rate[1] += (double)a->pivot[0][k] * unit;
      rate[2] += (double)a->pivot[1][k] * unit;
    }
    if (rate[2] == rate[1]) {
      return bad(diag, prefix, "ADMITTANCE", "two pivot waves of one frequency");
    }
    a->part[0] = (rate[2] - rate[0]) / (rate[2] - rate[1]) * row[1] / row[3];
    a->part[1] = (rate[0] - rate[1]) / (rate[2] - rate[1]) * row[1] / row[5];
  }
  return MYB_OK;
}

/* The loaded wave with these arguments: its index, or wave_count. */
static size_t find_wave(const myb_tides* t, const int8_t m[ARGUMENTS]) {
  for (size_t i = 0; i < t->wave_count; ++i) {
    if (memcmp(t->waves[i].multiplier, m, ARGUMENTS) == 0) {
      return i;
    }
  }
  return t->wave_count;
}

/* Attaches to the loaded waves the secondary waves of the table whose two pivots are loaded
 * and which are not loaded themselves. */
static myb_status relink(myb_tides* t) {
  myb_internal_release(&t->allocator, t->links);
  myb_internal_release(&t->allocator, t->link_start);
  t->links = NULL;
  t->link_start = NULL;
  t->secondary = 0;
  if (t->table_count == 0 || t->wave_count == 0) {
    return MYB_OK;
  }
  void* block = NULL;
  myb_status status =
      myb_internal_alloc_array(&t->allocator, t->wave_count + 1, sizeof(uint32_t), &block);
  if (status != MYB_OK) {
    return status;
  }
  uint32_t* start = block;
  memset(start, 0, (t->wave_count + 1) * sizeof(uint32_t));
  /* two passes: count the links of each pivot, then place them */
  for (int pass = 0; pass < 2; ++pass) {
    for (size_t r = 0; r < t->table_count; ++r) {
      const admittance_row* a = &t->table[r];
      const size_t first = find_wave(t, a->pivot[0]);
      const size_t second = find_wave(t, a->pivot[1]);
      if (first == t->wave_count || second == t->wave_count ||
          find_wave(t, a->wave) != t->wave_count) {
        continue;
      }
      const size_t pivot[2] = {first, second};
      for (int k = 0; k < 2; ++k) {
        if (pass == 1) {
          ocean_link* link = &t->links[start[pivot[k]]];
          memcpy(link->multiplier, a->wave, ARGUMENTS);
          link->factor = a->part[k];
        }
        ++start[pivot[k]];
      }
      if (pass == 0) {
        ++t->secondary;
      }
    }
    if (pass == 0) {
      /* counts -> the offset of each pivot's links */
      uint32_t total = 0;
      for (size_t i = 0; i <= t->wave_count; ++i) {
        const uint32_t count = start[i];
        start[i] = total;
        total += count;
      }
      if (total == 0) {
        myb_internal_release(&t->allocator, start);
        return MYB_OK;
      }
      status = myb_internal_alloc_array(&t->allocator, total, sizeof(ocean_link), &block);
      if (status != MYB_OK) {
        myb_internal_release(&t->allocator, start);
        t->secondary = 0;
        return status;
      }
      t->links = block;
    }
  }
  /* placing moved every offset to the end of its pivot's links: back by one pivot */
  for (size_t i = t->wave_count; i > 0; --i) {
    start[i] = start[i - 1];
  }
  start[0] = 0;
  t->link_start = start;
  return MYB_OK;
}

static int split(const char* s, size_t n, const char** tok, size_t* len, int max) {
  int count = 0;
  size_t i = 0;
  /* bound: every iteration consumes at least one character */
  while (i < n && count < max) {
    /* bound: blanks */
    while (i < n && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r')) {
      ++i;
    }
    if (i == n) {
      break;
    }
    const size_t start = i;
    /* bound: a token */
    while (i < n && s[i] != ' ' && s[i] != '\t' && s[i] != '\r') {
      ++i;
    }
    tok[count] = s + start;
    len[count] = i - start;
    ++count;
  }
  return count;
}

static int number(const char* s, size_t n, double* value) {
  size_t used = 0;
  return myb_parse_double(s, n, value, &used) == MYB_OK && used == n && finite(*value);
}

/* The waves and coefficients of a text, in arrays of the allocator of `tides`. */
static myb_status read_waves(const myb_tides* tides, const char* text, size_t size,
                             int max_degree, myb_diag* diag, ocean_wave** out_waves,
                             size_t* out_wave_count, ocean_entry** out_entries,
                             size_t* out_entry_count) {
  if (max_degree < 2 || max_degree > OCEAN_MAX_DEGREE) {
    myb_internal_diag(diag, 0, "maximum degree outside 2 .. 360");
    return MYB_E_DOMAIN;
  }
  /* two passes: count, then fill */
  size_t waves = 0;
  size_t entries = 0;
  ocean_wave* w = NULL;
  ocean_entry* e = NULL;
  myb_status status = MYB_OK;
  for (int pass = 0; pass < 2 && status == MYB_OK; ++pass) {
    size_t at = 0;
    long line_number = 0;
    size_t nw = 0;
    size_t ne = 0;
    double last = -1.0;
    /* bound: one line per iteration, at most the size of the text */
    while (at < size && status == MYB_OK) {
      const size_t start = at;
      /* bound: one line */
      while (at < size && text[at] != '\n') {
        ++at;
      }
      const size_t length = at - start;
      ++at;
      ++line_number;
      const char* tok[10];
      size_t len[10];
      const int k = split(text + start, length, tok, len, 10);
      int8_t mult[ARGUMENTS];
      if (k == 0 || !doodson(tok[0], len[0], mult)) {
        continue; /* header and comment lines */
      }
      double values[6];
      int ok = k == 8;
      for (int i = 0; ok && i < 6; ++i) {
        ok = number(tok[2 + i], len[2 + i], &values[i]);
      }
      if (!ok || values[0] != floor(values[0]) || values[1] != floor(values[1]) ||
          values[0] < 0.0 || values[1] < 0.0 || values[1] > values[0] ||
          values[0] > (double)OCEAN_MAX_DEGREE) {
        myb_internal_diag(diag, line_number, "ocean tides: malformed line");
        status = MYB_E_FORMAT;
        break;
      }
      double code = 0.0;
      (void)number(tok[0], len[0], &code);
      if (code != last) { /* a new wave */
        if (pass == 1) {
          memcpy(w[nw].multiplier, mult, sizeof mult);
          w[nw].doodson = code;
        }
        ++nw;
        last = code;
      }
      const int n = (int)values[0];
      const int m = (int)values[1];
      if (n < 2 || n > max_degree) {
        continue; /* degree 1 (geocentre motion) is not a change of the field about the CM */
      }
      if (pass == 1) {
        e[ne] = (ocean_entry){.wave = (uint32_t)(nw - 1),
                              .n = (int16_t)n,
                              .m = (int16_t)m,
                              .cp = values[2] * 1e-11,
                              .sp = values[3] * 1e-11,
                              .cm = values[4] * 1e-11,
                              .sm = values[5] * 1e-11};
      }
      ++ne;
    }
    if (status != MYB_OK) {
      break;
    }
    if (pass == 0) {
      if (nw == 0 || ne == 0) {
        myb_internal_diag(diag, 0, "ocean tides: no coefficients");
        status = MYB_E_FORMAT;
        break;
      }
      waves = nw;
      entries = ne;
      void* block = NULL;
      status = myb_internal_alloc_array(&tides->allocator, waves, sizeof(ocean_wave), &block);
      if (status == MYB_OK) {
        w = block;
        status = myb_internal_alloc_array(&tides->allocator, entries, sizeof(ocean_entry), &block);
        if (status == MYB_OK) {
          e = block;
        }
      }
    }
  }
  if (status != MYB_OK) {
    myb_internal_release(&tides->allocator, w);
    myb_internal_release(&tides->allocator, e);
    return status;
  }
  *out_waves = w;
  *out_wave_count = waves;
  *out_entries = e;
  *out_entry_count = entries;
  return MYB_OK;
}

myb_status myb_tides_load_ocean(myb_tides* tides, const char* text, size_t size, int max_degree,
                                myb_diag* diag) {
  MYB_CHECK(tides != NULL && text != NULL);
  myb_internal_diag_clear(diag);
  ocean_wave* w = NULL;
  ocean_entry* e = NULL;
  size_t waves = 0;
  size_t entries = 0;
  const myb_status status =
      read_waves(tides, text, size, max_degree, diag, &w, &waves, &e, &entries);
  if (status != MYB_OK) {
    return status;
  }
  myb_internal_release(&tides->allocator, tides->waves);
  myb_internal_release(&tides->allocator, tides->entries);
  tides->waves = w;
  tides->wave_count = waves;
  tides->entries = e;
  tides->entry_count = entries;
  return relink(tides);
}

myb_status myb_tides_add_waves(myb_tides* tides, const char* text, size_t size, int max_degree,
                               myb_diag* diag) {
  MYB_CHECK(tides != NULL && text != NULL);
  myb_internal_diag_clear(diag);
  ocean_wave* w = NULL;
  ocean_entry* e = NULL;
  size_t waves = 0;
  size_t entries = 0;
  myb_status status = read_waves(tides, text, size, max_degree, diag, &w, &waves, &e, &entries);
  if (status != MYB_OK) {
    return status;
  }
  /* those loaded, then the new ones */
  void* block = NULL;
  ocean_wave* all_waves = NULL;
  ocean_entry* all_entries = NULL;
  status = myb_internal_alloc_array(&tides->allocator, tides->wave_count + waves,
                                    sizeof(ocean_wave), &block);
  if (status == MYB_OK) {
    all_waves = block;
    status = myb_internal_alloc_array(&tides->allocator, tides->entry_count + entries,
                                      sizeof(ocean_entry), &block);
  }
  if (status == MYB_OK) {
    all_entries = block;
    if (tides->wave_count > 0) {
      memcpy(all_waves, tides->waves, tides->wave_count * sizeof(ocean_wave));
    }
    if (tides->entry_count > 0) {
      memcpy(all_entries, tides->entries, tides->entry_count * sizeof(ocean_entry));
    }
    memcpy(all_waves + tides->wave_count, w, waves * sizeof(ocean_wave));
    for (size_t k = 0; k < entries; ++k) {
      all_entries[tides->entry_count + k] = e[k];
      all_entries[tides->entry_count + k].wave += (uint32_t)tides->wave_count;
    }
    myb_internal_release(&tides->allocator, tides->waves);
    myb_internal_release(&tides->allocator, tides->entries);
    tides->waves = all_waves;
    tides->wave_count += waves;
    tides->entries = all_entries;
    tides->entry_count += entries;
    status = relink(tides);
  } else {
    myb_internal_release(&tides->allocator, all_waves);
  }
  myb_internal_release(&tides->allocator, w);
  myb_internal_release(&tides->allocator, e);
  return status;
}

myb_status myb_tides_add_waves_file(myb_tides* tides, const char* path, int max_degree,
                                    myb_diag* diag) {
  MYB_CHECK(tides != NULL && path != NULL);
  myb_internal_diag_clear(diag);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &tides->allocator, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded =
      myb_tides_add_waves(tides, (const char*)map.data, map.size, max_degree, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

myb_status myb_tides_load_ocean_file(myb_tides* tides, const char* path, int max_degree,
                                     myb_diag* diag) {
  MYB_CHECK(tides != NULL && path != NULL);
  myb_internal_diag_clear(diag);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &tides->allocator, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded =
      myb_tides_load_ocean(tides, (const char*)map.data, map.size, max_degree, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

/* ---- Application ----------------------------------------------------------------------------- */

static void add(myb_gravity_field* f, int degree, int n, int m, double dc, double ds) {
  if (n <= degree) {
    const myb_status status = myb_gravity_field_add_correction(f, n, m, dc, ds);
    MYB_CHECK(status == MYB_OK);
  }
}

myb_status myb_tides_apply(const myb_tides* tides, myb_gravity_field* field,
                           const myb_tides_input* input) {
  MYB_CHECK(tides != NULL && field != NULL && input != NULL);
  MYB_CHECK(input->count == 0 || input->raisers != NULL);
  const int degree = myb_gravity_field_degree(field);
  const double gm = myb_gravity_field_gm(field);
  const double radius = myb_gravity_field_radius(field);
  const int flags = input->flags;
  int remove_permanent = 0;
  if (flags & MYB_TIDES_SOLID) {
    const char* system = myb_gravity_field_tide_system(field);
    if (input->permanent == MYB_PERMANENT_TIDE_AUTO) {
      if (strcmp(system, "zero_tide") == 0) {
        remove_permanent = 1;
      } else if (strcmp(system, "tide_free") != 0) {
        return MYB_E_DOMAIN; /* mean tide or undeclared: the caller must decide */
      }
    } else {
      remove_permanent = input->permanent == MYB_PERMANENT_TIDE_REMOVE;
    }
  }
  if (!finite(input->xp) || !finite(input->yp)) {
    return MYB_E_DOMAIN;
  }
  for (size_t j = 0; j < input->count; ++j) {
    const myb_tide_raiser* b = &input->raisers[j];
    const double r2 = b->r[0] * b->r[0] + b->r[1] * b->r[1] + b->r[2] * b->r[2];
    if (!finite(b->gm) || !(r2 > 0.0) || !finite(r2)) {
      return MYB_E_DOMAIN;
    }
  }
  /* step 1: X_nm - i Y_nm = 1/(2n + 1) sum_j (GM_j / GM)(R / r_j)^(n+1) Pbar_nm e^(-i m lambda) */
  if (flags & MYB_TIDES_SOLID) {
    double x[4][4] = {{0.0}};
    double y[4][4] = {{0.0}};
    for (size_t j = 0; j < input->count; ++j) {
      const myb_tide_raiser* b = &input->raisers[j];
      const double r = sqrt(b->r[0] * b->r[0] + b->r[1] * b->r[1] + b->r[2] * b->r[2]);
      const double t = b->r[2] / r;     /* sin(latitude) */
      const double c1 = b->r[0] / r;    /* u cos(lambda) */
      const double s1 = b->r[1] / r;    /* u sin(lambda) */
      const double c2 = c1 * c1 - s1 * s1; /* u^2 cos(2 lambda) */
      const double s2 = 2.0 * c1 * s1;
      const double c3 = c2 * c1 - s2 * s1;
      const double s3 = s2 * c1 + c2 * s1;
      const double q = radius / r;
      const double w2 = b->gm / gm * q * q * q / 5.0;
      const double w3 = b->gm / gm * q * q * q * q / 7.0;
      /* Pbar_nm / u^m */
      const double p20 = 1.118033988749895 * (3.0 * t * t - 1.0); /* sqrt(5) / 2 */
      const double p21 = 3.872983346207417 * t;                    /* sqrt(15) */
      const double p22 = 1.9364916731037085;                       /* sqrt(15) / 2 */
      const double p30 = 1.3228756555322954 * (5.0 * t * t - 3.0) * t; /* sqrt(7) / 2 */
      const double p31 = 1.6201851746019651 * (5.0 * t * t - 1.0);     /* sqrt(42) / 4 */
      const double p32 = 5.123475382979799 * t;                        /* sqrt(105) / 2 */
      const double p33 = 2.091650066335189;                            /* sqrt(70) / 4 */
      x[2][0] += w2 * p20;
      x[2][1] += w2 * p21 * c1;
      y[2][1] += w2 * p21 * s1;
      x[2][2] += w2 * p22 * c2;
      y[2][2] += w2 * p22 * s2;
      x[3][0] += w3 * p30;
      x[3][1] += w3 * p31 * c1;
      y[3][1] += w3 * p31 * s1;
      x[3][2] += w3 * p32 * c2;
      y[3][2] += w3 * p32 * s2;
      x[3][3] += w3 * p33 * c3;
      y[3][3] += w3 * p33 * s3;
    }
    for (int n = 2; n <= 3; ++n) {
      for (int m = 0; m <= n; ++m) {
        if (!tides->has_love[n][m]) {
          continue;
        }
        /* dC - i dS = (kR + i kI)(X - i Y) */
        const double kr = tides->k_re[n][m];
        const double ki = tides->k_im[n][m];
        add(field, degree, n, m, kr * x[n][m] + ki * y[n][m],
            m == 0 ? 0.0 : kr * y[n][m] - ki * x[n][m]);
        if (n == 2 && tides->k_plus[m] != 0.0) {
          add(field, degree, 4, m, tides->k_plus[m] * x[2][m],
              m == 0 ? 0.0 : tides->k_plus[m] * y[2][m]);
        }
      }
    }
    if (remove_permanent) {
      add(field, degree, 2, 0, -tides->permanent * tides->k_re[2][0], 0.0);
    }
  }
  const double centuries = myb_time_seconds(input->t) / seconds_per_century;
  if (flags & MYB_TIDES_FREQUENCY) {
    double v[MYB_POISSON_MAX_SERIES];
    double rates[MYB_POISSON_MAX_SERIES];
    myb_internal_poisson_eval(&tides->step2, centuries, v, rates);
    add(field, degree, 2, 0, v[0], 0.0);
    add(field, degree, 2, 1, v[1], v[2]);
    add(field, degree, 2, 2, v[3], v[4]);
  }
  if (flags & (MYB_TIDES_POLE_SOLID | MYB_TIDES_POLE_OCEAN)) {
    const double years = myb_time_seconds(input->t) / seconds_per_year;
    const double xs = (tides->secular_pole[0] + tides->secular_pole[1] * years) * 1e-3;
    const double ys = (tides->secular_pole[2] + tides->secular_pole[3] * years) * 1e-3;
    const double m1 = input->xp * arcseconds_per_radian - xs;
    const double m2 = -(input->yp * arcseconds_per_radian - ys);
    double dc = 0.0;
    double ds = 0.0;
    if (flags & MYB_TIDES_POLE_SOLID) {
      const double* c = tides->pole_solid;
      dc += c[0] * (m1 + c[1] * m2);
      ds += c[2] * (m2 + c[3] * m1);
    }
    if (flags & MYB_TIDES_POLE_OCEAN) {
      const double* c = tides->pole_ocean;
      dc += c[0] * (m1 + c[1] * m2);
      ds += c[2] * (m2 + c[3] * m1);
    }
    add(field, degree, 2, 1, dc, ds);
  }
  if ((flags & MYB_TIDES_OCEAN) && tides->entry_count > 0) {
    double angle[MYB_POISSON_MAX_ARGUMENTS];
    myb_internal_poisson_arguments(&tides->step2, centuries, angle);
    size_t current = SIZE_MAX;
    double ct = 1.0;
    double st = 0.0;
    double* fc = myb_internal_field_c(field);
    double* fs = myb_internal_field_s(field);
    int touched = -1;
    for (size_t k = 0; k < tides->entry_count; ++k) {
      const ocean_entry* e = &tides->entries[k];
      if (e->n > degree) {
        continue;
      }
      if (e->wave != current) {
        current = e->wave;
        double theta = 0.0;
        for (int i = 0; i < ARGUMENTS; ++i) {
          theta += tides->waves[current].multiplier[i] * angle[i];
        }
        ct = myb_cos(theta);
        st = myb_sin(theta);
        if ((flags & MYB_TIDES_ADMITTANCE) && tides->link_start != NULL) {
          /* the secondary waves that take a part of this wave's coefficients (eq. 6.16) */
          for (uint32_t q = tides->link_start[current]; q < tides->link_start[current + 1];
               ++q) {
            const ocean_link* link = &tides->links[q];
            double phase = 0.0;
            for (int i = 0; i < ARGUMENTS; ++i) {
              phase += link->multiplier[i] * angle[i];
            }
            ct += link->factor * myb_cos(phase);
            st += link->factor * myb_sin(phase);
          }
        }
      }
      /* dC - i dS = sum (C+ - i S+) e^(i theta) + (C- + i S-) e^(-i theta) (eq. 6.15) */
      const double dc = (e->cp + e->cm) * ct + (e->sp + e->sm) * st;
      const double ds = e->m == 0 ? 0.0 : (e->sp - e->sm) * ct - (e->cp - e->cm) * st;
      const size_t i = (size_t)e->n * (size_t)(e->n + 1) / 2 + (size_t)e->m;
      fc[i] += dc;
      fs[i] += ds;
      touched = e->n > touched ? e->n : touched;
    }
    if (touched >= 0) {
      myb_internal_field_touched(field, touched);
    }
  }
  return MYB_OK;
}

size_t myb_tides_ocean_waves(const myb_tides* tides) {
  MYB_CHECK(tides != NULL);
  return tides->wave_count;
}

size_t myb_tides_secondary_waves(const myb_tides* tides) {
  MYB_CHECK(tides != NULL);
  return tides->secondary;
}
