/* SPDX-License-Identifier: BSD-3-Clause */
#include "num/poisson.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "num/eft.h"

static myb_status numbers(const myb_kpool* pool, const char* name, const double** values,
                          size_t* count, myb_diag* diag) {
  const myb_status status = myb_kpool_numbers(pool, name, values, count);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "%s: %s", name,
                      status == MYB_E_NOT_FOUND ? "not defined" : "not numeric");
  }
  return status;
}

static myb_status bad(myb_diag* diag, const char* name, const char* what) {
  myb_internal_diag(diag, 0, "%s: %s", name, what);
  return MYB_E_FORMAT;
}

/* x as an integer in [lo, hi]. */
static int as_int(double x, int lo, int hi, int* out) {
  if (!(x >= lo && x <= hi) || x != (double)(int)x) {
    return 0;
  }
  *out = (int)x;
  return 1;
}

static myb_status load_arguments(const myb_kpool* pool, const char* prefix,
                                 myb_internal_poisson_group* g, myb_diag* diag) {
  char name[128];
  const double* v = NULL;
  size_t n = 0;
  snprintf(name, sizeof name, "%s_ARGUMENT_COUNT", prefix);
  myb_status status = numbers(pool, name, &v, &n, diag);
  if (status != MYB_OK) {
    return status;
  }
  int count = 0;
  if (n != 1 || !as_int(v[0], 1, MYB_POISSON_MAX_ARGUMENTS, &count)) {
    return bad(diag, name, "expected one integer in 1 .. 16");
  }
  g->arguments = count;
  snprintf(name, sizeof name, "%s_ARGUMENT_POLYNOMIAL", prefix);
  status = numbers(pool, name, &v, &n, diag);
  if (status != MYB_OK) {
    return status;
  }
  if (n != (size_t)count * (MYB_POISSON_ARGUMENT_DEGREE + 1)) {
    return bad(diag, name, "expected 5 coefficients per argument");
  }
  for (int k = 0; k < count; ++k) {
    for (int i = 0; i <= MYB_POISSON_ARGUMENT_DEGREE; ++i) {
      const double c = v[(size_t)k * (MYB_POISSON_ARGUMENT_DEGREE + 1) + (size_t)i];
      if (!isfinite(c)) {
        return bad(diag, name, "non-finite coefficient");
      }
      g->argument_polynomial[k][i] = c;
    }
  }
  const char* suffixes[2] = {"SCALE", "TURN"};
  for (int which = 0; which < 2; ++which) {
    snprintf(name, sizeof name, "%s_ARGUMENT_%s", prefix, suffixes[which]);
    status = numbers(pool, name, &v, &n, diag);
    if (status != MYB_OK) {
      return status;
    }
    if (n != (size_t)count) {
      return bad(diag, name, "expected one value per argument");
    }
    for (int k = 0; k < count; ++k) {
      if (!(v[k] > 0.0 && v[k] < 1e300)) {
        return bad(diag, name, "values must be positive and finite");
      }
      if (which == 0) {
        g->argument_scale[k] = v[k];
      } else {
        g->argument_turn[k] = v[k];
      }
    }
  }
  return MYB_OK;
}

/* Scale and polynomial of series i; the terms are validated and returned (kernel memory). */
static myb_status load_series(const myb_kpool* pool, const char* prefix, int i,
                              myb_internal_poisson_group* g, const double** terms,
                              size_t* count, myb_diag* diag) {
  char name[128];
  const double* v = NULL;
  size_t n = 0;
  snprintf(name, sizeof name, "%s_SCALE", prefix);
  myb_status status = numbers(pool, name, &v, &n, diag);
  if (status != MYB_OK) {
    return status;
  }
  if (n != 1 || !(v[0] > 0.0 && v[0] < 1e300)) {
    return bad(diag, name, "expected one positive number");
  }
  g->scale[i] = v[0];
  snprintf(name, sizeof name, "%s_POLYNOMIAL", prefix);
  status = numbers(pool, name, &v, &n, diag);
  if (status != MYB_OK) {
    return status;
  }
  if (n < 1 || n > MYB_POISSON_MAX_POWER + 1) {
    return bad(diag, name, "expected 1 .. 9 coefficients");
  }
  for (size_t k = 0; k < n; ++k) {
    if (!isfinite(v[k])) {
      return bad(diag, name, "non-finite coefficient");
    }
    g->polynomial[i][k] = v[k];
  }
  g->degree[i] = (int)n - 1;
  snprintf(name, sizeof name, "%s_TERMS", prefix);
  status = numbers(pool, name, &v, &n, diag);
  if (status != MYB_OK) {
    return status;
  }
  const size_t width = 3 + (size_t)g->arguments;
  if (n % width != 0 || n / width > UINT32_MAX / 4) {
    return bad(diag, name, "expected 3 + ARGUMENT_COUNT numbers per term");
  }
  for (size_t r = 0; r < n / width; ++r) {
    const double* row = v + r * width;
    int value = 0;
    if (!as_int(row[0], 0, MYB_POISSON_MAX_POWER, &value) || !isfinite(row[1]) ||
        !isfinite(row[2])) {
      return bad(diag, name, "invalid power or amplitude");
    }
    for (int k = 0; k < g->arguments; ++k) {
      if (!as_int(row[3 + k], -MYB_POISSON_MAX_MULTIPLE, MYB_POISSON_MAX_MULTIPLE, &value)) {
        return bad(diag, name, "argument multiplier not an integer in -64 .. 64");
      }
      const int magnitude = value < 0 ? -value : value;
      if (magnitude > g->max_multiple[k]) {
        g->max_multiple[k] = magnitude;
      }
    }
  }
  *terms = v;
  *count = n / width;
  return MYB_OK;
}

typedef struct term_ref {
  const double* row;
  int series;
} term_ref;

/* Lexicographic order of the argument vectors (multipliers n_1 .. n_count). */
static int vector_before(const term_ref* a, const term_ref* b, int count) {
  for (int k = 0; k < count; ++k) {
    if (a->row[3 + k] != b->row[3 + k]) {
      return a->row[3 + k] < b->row[3 + k];
    }
  }
  return 0;
}

/* Stable bottom-up merge sort of the terms by argument vector. */
static void sort_terms(term_ref* terms, term_ref* scratch, size_t n, int count) {
  /* bound: the run width doubles, at most 64 passes */
  for (size_t width = 1; width < n; width *= 2) {
    for (size_t lo = 0; lo < n; lo += 2 * width) {
      const size_t mid = lo + width < n ? lo + width : n;
      const size_t hi = lo + 2 * width < n ? lo + 2 * width : n;
      size_t i = lo;
      size_t j = mid;
      for (size_t k = lo; k < hi; ++k) {
        if (i < mid && (j >= hi || !vector_before(&terms[j], &terms[i], count))) {
          scratch[k] = terms[i++];
        } else {
          scratch[k] = terms[j++];
        }
      }
    }
    memcpy(terms, scratch, n * sizeof *terms);
  }
}

void myb_internal_poisson_free(myb_internal_poisson_group* group,
                               const myb_allocator* allocator) {
  if (group == NULL) {
    return;
  }
  myb_internal_release(allocator, group->pair_start);
  myb_internal_release(allocator, group->pair_argument);
  myb_internal_release(allocator, group->pair_multiple);
  myb_internal_release(allocator, group->entry_start);
  myb_internal_release(allocator, group->entry_series);
  myb_internal_release(allocator, group->entry_power);
  myb_internal_release(allocator, group->entry_amplitude);
  memset(group, 0, sizeof *group);
}

myb_status myb_internal_poisson_load(const myb_kpool* pool, const char* prefix,
                                     const char* const* names, int count,
                                     const myb_allocator* allocator,
                                     myb_internal_poisson_group* group, myb_diag* diag) {
  MYB_CHECK(pool != NULL && prefix != NULL && names != NULL && count >= 1 &&
            count <= MYB_POISSON_MAX_SERIES && allocator != NULL && group != NULL);
  memset(group, 0, sizeof *group);
  myb_internal_poisson_group* g = group;
  myb_status status = load_arguments(pool, prefix, g, diag);
  if (status != MYB_OK) {
    return status;
  }
  g->series = count;
  const double* rows[MYB_POISSON_MAX_SERIES] = {NULL};
  size_t counts[MYB_POISSON_MAX_SERIES] = {0};
  size_t total = 0;
  for (int i = 0; i < count; ++i) {
    status = load_series(pool, names[i], i, g, &rows[i], &counts[i], diag);
    if (status != MYB_OK) {
      return status;
    }
    total += counts[i];
  }
  term_ref* terms = NULL;
  term_ref* scratch = NULL;
  void* block = NULL;
  status = myb_internal_alloc_array(allocator, total, 2 * sizeof(term_ref), &block);
  if (status != MYB_OK) {
    return status;
  }
  terms = block;
  scratch = terms + total;
  const size_t width = 3 + (size_t)g->arguments;
  size_t next = 0;
  for (int i = 0; i < count; ++i) {
    for (size_t r = 0; r < counts[i]; ++r) {
      terms[next++] = (term_ref){rows[i] + r * width, i};
    }
  }
  sort_terms(terms, scratch, total, g->arguments);
  size_t unique = 0;
  size_t pairs = 0;
  for (size_t r = 0; r < total; ++r) {
    if (r == 0 || vector_before(&terms[r - 1], &terms[r], g->arguments)) {
      ++unique;
      for (int k = 0; k < g->arguments; ++k) {
        pairs += terms[r].row[3 + k] != 0.0;
      }
    }
  }
  void* blocks[7] = {NULL};
  const size_t sizes[7][2] = {{unique + 1, sizeof(uint32_t)}, {pairs, 1}, {pairs, 1},
                              {unique + 1, sizeof(uint32_t)}, {total, 1}, {total, 1},
                              {total, 2 * sizeof(double)}};
  for (int b = 0; b < 7 && status == MYB_OK; ++b) {
    status = myb_internal_alloc_array(allocator, sizes[b][0], sizes[b][1], &blocks[b]);
  }
  if (status != MYB_OK) {
    for (int b = 0; b < 7; ++b) {
      myb_internal_release(allocator, blocks[b]);
    }
    myb_internal_release(allocator, terms);
    return status;
  }
  g->pair_start = blocks[0];
  g->pair_argument = blocks[1];
  g->pair_multiple = blocks[2];
  g->entry_start = blocks[3];
  g->entry_series = blocks[4];
  g->entry_power = blocks[5];
  g->entry_amplitude = blocks[6];
  g->unique = unique;
  size_t u = 0;
  size_t p = 0;
  for (size_t r = 0; r < total; ++r) {
    if (r == 0 || vector_before(&terms[r - 1], &terms[r], g->arguments)) {
      g->pair_start[u] = (uint32_t)p;
      g->entry_start[u] = (uint32_t)r;
      for (int k = 0; k < g->arguments; ++k) {
        const double n = terms[r].row[3 + k];
        if (n != 0.0) {
          g->pair_argument[p] = (uint8_t)k;
          g->pair_multiple[p] = (int8_t)n;
          ++p;
        }
      }
      ++u;
    }
    g->entry_series[r] = (uint8_t)terms[r].series;
    g->entry_power[r] = (uint8_t)terms[r].row[0];
    g->entry_amplitude[2 * r] = terms[r].row[1];
    g->entry_amplitude[2 * r + 1] = terms[r].row[2];
  }
  g->pair_start[unique] = (uint32_t)pairs;
  g->entry_start[unique] = (uint32_t)total;
  myb_internal_release(allocator, terms);
  return MYB_OK;
}

/* x reduced to [0, turn) to within one rounding: x - q turn with q turn split exactly. */
static double reduce(double x, double turn) {
  const double q = floor(x / turn);
  double low;
  const double high = myb_eft_two_prod(q, turn, &low);
  double r = (x - high) - low;
  if (r < 0.0) {
    r += turn;
  } else if (r >= turn) {
    r -= turn;
  }
  return r;
}

void myb_internal_poisson_arguments(const myb_internal_poisson_group* g, double t,
                                    double* angles) {
  MYB_CHECK(g != NULL && angles != NULL);
  for (int k = 0; k < g->arguments; ++k) {
    const double* c = g->argument_polynomial[k];
    const double value = c[0] + t * (c[1] + t * (c[2] + t * (c[3] + t * c[4])));
    angles[k] = reduce(value, g->argument_turn[k]) * g->argument_scale[k];
  }
}

void myb_internal_poisson_eval(const myb_internal_poisson_group* g, double t, double* values,
                               double* rates) {
  MYB_CHECK(g != NULL && values != NULL && rates != NULL);
  /* phasors exp(i m F_k) and the rates dF_k/dt */
  double re[MYB_POISSON_MAX_ARGUMENTS][MYB_POISSON_MAX_MULTIPLE + 1];
  double im[MYB_POISSON_MAX_ARGUMENTS][MYB_POISSON_MAX_MULTIPLE + 1];
  double argument_rate[MYB_POISSON_MAX_ARGUMENTS];
  for (int k = 0; k < g->arguments; ++k) {
    const double* c = g->argument_polynomial[k];
    const double value = c[0] + t * (c[1] + t * (c[2] + t * (c[3] + t * c[4])));
    const double rate = c[1] + t * (2.0 * c[2] + t * (3.0 * c[3] + t * 4.0 * c[4]));
    const double angle = reduce(value, g->argument_turn[k]) * g->argument_scale[k];
    argument_rate[k] = rate * g->argument_scale[k];
    re[k][0] = 1.0;
    im[k][0] = 0.0;
    for (int m = 1; m <= g->max_multiple[k]; ++m) {
      const double multiple = (double)m * angle;
      re[k][m] = myb_cos(multiple);
      im[k][m] = myb_sin(multiple);
    }
  }
  /* sums per series and power: compensated for the values, plain for the rates */
  double sum[MYB_POISSON_MAX_SERIES][MYB_POISSON_MAX_POWER + 1] = {{0}};
  double comp[MYB_POISSON_MAX_SERIES][MYB_POISSON_MAX_POWER + 1] = {{0}};
  double dsum[MYB_POISSON_MAX_SERIES][MYB_POISSON_MAX_POWER + 1] = {{0}};
  for (size_t u = 0; u < g->unique; ++u) {
    double pr = 1.0;
    double pi = 0.0;
    double theta_dot = 0.0;
    for (uint32_t q = g->pair_start[u]; q < g->pair_start[u + 1]; ++q) {
      const int k = g->pair_argument[q];
      const int m = g->pair_multiple[q];
      const int magnitude = m < 0 ? -m : m;
      const double zr = re[k][magnitude];
      const double zi = m < 0 ? -im[k][magnitude] : im[k][magnitude];
      const double r = pr * zr - pi * zi;
      pi = pr * zi + pi * zr;
      pr = r;
      theta_dot += (double)m * argument_rate[k];
    }
    for (uint32_t e = g->entry_start[u]; e < g->entry_start[u + 1]; ++e) {
      const int s = g->entry_series[e];
      const int j = g->entry_power[e];
      const double as = g->entry_amplitude[2 * e];
      const double ac = g->entry_amplitude[2 * e + 1];
      double err;
      sum[s][j] = myb_eft_two_sum(sum[s][j], as * pi + ac * pr, &err);
      comp[s][j] += err;
      dsum[s][j] += theta_dot * (as * pr - ac * pi); /* rates: plain sums suffice */
    }
  }
  /* Q = sum_j t^j p_j and dQ/dt = sum_j (j t^(j-1) p_j + t^j dp_j), p_j = polynomial part plus
   * the periodic terms of power j; by Horner from the highest power. */
  for (int s = 0; s < g->series; ++s) {
    double v = 0.0;
    double d = 0.0;
    double periodic_rate = 0.0;
    for (int j = MYB_POISSON_MAX_POWER; j >= 0; --j) {
      const double p = (j <= g->degree[s] ? g->polynomial[s][j] : 0.0) + (sum[s][j] + comp[s][j]);
      v = v * t + p;
      if (j >= 1) {
        d = d * t + (double)j * p;
      }
      periodic_rate = periodic_rate * t + dsum[s][j];
    }
    values[s] = v * g->scale[s];
    rates[s] = (d + periodic_rate) * g->scale[s];
  }
}
