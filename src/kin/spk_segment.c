/* SPDX-License-Identifier: BSD-3-Clause */
#include "kin/spk_segment.h"

#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "myb/kepler.h"
#include "myb/num.h"
#include "num/vec3.h"

static const double pi = 0x1.921fb54442d18p+1;
static const double two_pi = 0x1.921fb54442d18p+2;

/* ---- Reading ----------------------------------------------------------------------------- */

static double word(const myb_internal_spk_segment* s, uint64_t address) {
  MYB_DCHECK(address >= s->begin && address <= s->end);
  return myb_internal_daf_double(&s->daf, address);
}

/* Epochs are seconds past J2000; the bound keeps myb_time arithmetic exact. */
static int epoch_ok(double x) {
  return x >= -0x1p62 && x <= 0x1p62;
}

static int positive(double x) {
  return x > 0.0 && x <= 0x1p1000;
}

/* x as a count in [0, max]. */
static int as_count(double x, uint64_t max, uint64_t* out) {
  if (!(x >= 0.0 && x <= (double)max) || x != (double)(uint64_t)x) {
    return 0;
  }
  *out = (uint64_t)x;
  return 1;
}

/* RN(t - x) for an epoch x: the difference is exact in double-word form, then rounded once. */
static double since(myb_time t, double x) {
  return myb_time_diff(t, myb_time_from_seconds(x)).hi;
}

static int compare_epoch(myb_time t, double x) {
  return myb_time_compare(t, myb_time_from_seconds(x));
}

static myb_status format_error(myb_diag* diag, const myb_internal_spk_segment* s,
                               const char* what) {
  myb_internal_diag(diag, 0, "SPK segment '%s' (type %ld, target %ld): %s", s->name, (long)s->type,
                    (long)s->target, what);
  return MYB_E_FORMAT;
}

/* ---- Epoch lists ------------------------------------------------------------------------- */

static double epoch_at(const myb_internal_spk_segment* s, uint64_t epochs, uint64_t i) {
  if (epochs == 0) {
    return s->init + (double)i * s->step;
  }
  return word(s, epochs + i);
}

/* Number of epochs <= t among the n epochs (binary search). */
static uint64_t count_le(const myb_internal_spk_segment* s, uint64_t epochs, uint64_t n,
                         myb_time t) {
  uint64_t lo = 0;
  uint64_t hi = n;
  /* bound: the interval halves, at most 64 iterations */
  while (lo < hi) {
    const uint64_t mid = lo + (hi - lo) / 2;
    if (compare_epoch(t, epoch_at(s, epochs, mid)) >= 0) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

/* Number of epochs < t. */
static uint64_t count_lt(const myb_internal_spk_segment* s, uint64_t epochs, uint64_t n,
                         myb_time t) {
  uint64_t lo = 0;
  uint64_t hi = n;
  /* bound: the interval halves, at most 64 iterations */
  while (lo < hi) {
    const uint64_t mid = lo + (hi - lo) / 2;
    if (compare_epoch(t, epoch_at(s, epochs, mid)) > 0) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

/* The n epochs at `epochs` are finite and strictly increasing. */
static int increasing(const myb_internal_spk_segment* s, uint64_t epochs, uint64_t n) {
  double previous = -HUGE_VAL;
  for (uint64_t i = 0; i < n; ++i) {
    const double e = word(s, epochs + i);
    if (!epoch_ok(e) || !(e > previous)) {
      return 0;
    }
    previous = e;
  }
  return 1;
}

/* ---- Chebyshev series -------------------------------------------------------------------- */

/* sum_{k<n} c_k T_k(x) (Clenshaw) for the coefficients at `address`. */
static double chebyshev(const myb_internal_spk_segment* s, uint64_t address, int n, double x) {
  double b1 = 0.0;
  double b2 = 0.0;
  for (int k = n - 1; k >= 1; --k) {
    const double b0 = 2.0 * x * b1 - b2 + word(s, address + (uint64_t)k);
    b2 = b1;
    b1 = b0;
  }
  return x * b1 - b2 + word(s, address);
}

/* The same and its derivative with respect to x. */
static void chebyshev_derivative(const myb_internal_spk_segment* s, uint64_t address, int n,
                                 double x, double* value, double* derivative) {
  double b1 = 0.0;
  double b2 = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;
  for (int k = n - 1; k >= 1; --k) {
    const double b0 = 2.0 * x * b1 - b2 + word(s, address + (uint64_t)k);
    const double d0 = 2.0 * b1 + 2.0 * x * d1 - d2;
    b2 = b1;
    b1 = b0;
    d2 = d1;
    d1 = d0;
  }
  *value = x * b1 - b2 + word(s, address);
  *derivative = b1 + x * d1 - d2;
}

/* ---- Types 2 and 3: Chebyshev, equal intervals ------------------------------------------- */

static myb_status parse_chebyshev(myb_internal_spk_segment* s, int components, myb_diag* diag) {
  const uint64_t length = s->end - s->begin + 1;
  uint64_t rsize = 0;
  uint64_t n = 0;
  if (length < 4 + 2 + (uint64_t)components) {
    return format_error(diag, s, "too short");
  }
  s->init = word(s, s->end - 3);
  s->step = word(s, s->end - 2);
  if (!epoch_ok(s->init) || !positive(s->step) || !as_count(word(s, s->end - 1), length, &rsize) ||
      !as_count(word(s, s->end), length, &n) || n == 0 || rsize < 2 + (uint64_t)components ||
      (rsize - 2) % (uint64_t)components != 0 || n > (length - 4) / rsize ||
      n * rsize + 4 != length) {
    return format_error(diag, s, "inconsistent directory");
  }
  const uint64_t coefficients = (rsize - 2) / (uint64_t)components;
  if (coefficients > MYB_SPK_MAX_COEFFICIENTS) {
    return format_error(diag, s, "polynomial degree too high");
  }
  s->count = n;
  s->size = rsize;
  s->window = (int)coefficients;
  return MYB_OK;
}

/* Index of the equal-length record containing t (the later one at a common boundary). */
static uint64_t record_index(const myb_internal_spk_segment* s, double offset) {
  const double q = offset / s->step;
  if (!(q >= 0.0)) {
    return 0;
  }
  if (q >= (double)(s->count - 1)) {
    return s->count - 1;
  }
  return (uint64_t)q;
}

static myb_status state_chebyshev(const myb_internal_spk_segment* s, myb_time t,
                                  double state[6]) {
  const uint64_t record = s->begin + record_index(s, since(t, s->init)) * s->size;
  const double mid = word(s, record);
  const double radius = word(s, record + 1);
  if (!epoch_ok(mid) || !positive(radius)) {
    return MYB_E_FORMAT;
  }
  const double x = since(t, mid) / radius;
  const int n = s->window;
  for (int j = 0; j < 3; ++j) {
    const uint64_t address = record + 2 + (uint64_t)(j * n);
    if (s->type == 2) {
      double value;
      double derivative;
      chebyshev_derivative(s, address, n, x, &value, &derivative);
      state[j] = value;
      state[j + 3] = derivative / radius;
    } else {
      state[j] = chebyshev(s, address, n, x);
      state[j + 3] = chebyshev(s, address + (uint64_t)(3 * n), n, x);
    }
  }
  return MYB_OK;
}

/* ---- Type 20: Chebyshev velocity, position by integration -------------------------------- */

static myb_status parse_type20(myb_internal_spk_segment* s, myb_diag* diag) {
  const uint64_t length = s->end - s->begin + 1;
  uint64_t rsize = 0;
  uint64_t n = 0;
  if (length < 7 + 6) {
    return format_error(diag, s, "too short");
  }
  s->dscale = word(s, s->end - 6);
  s->tscale = word(s, s->end - 5);
  const double jd = word(s, s->end - 4);
  const double fraction = word(s, s->end - 3);
  const double days = word(s, s->end - 2);
  if (!positive(s->dscale) || !positive(s->tscale) || !(fabs(jd) < 1e9) ||
      !(fabs(fraction) < 1.0) || !(days > 0.0 && days < 1e9) ||
      !as_count(word(s, s->end - 1), length, &rsize) || !as_count(word(s, s->end), length, &n) ||
      n == 0 || rsize < 6 || rsize % 3 != 0 || n > (length - 7) / rsize ||
      n * rsize + 7 != length) {
    return format_error(diag, s, "inconsistent directory");
  }
  if (rsize / 3 > MYB_SPK_MAX_COEFFICIENTS) {
    return format_error(diag, s, "polynomial degree too high");
  }
  s->init_time = myb_time_from_jd(jd, fraction);
  s->step = days * 86400.0;
  s->count = n;
  s->size = rsize;
  s->window = (int)(rsize / 3) - 1; /* velocity coefficients per component */
  return MYB_OK;
}

static myb_status state_type20(const myb_internal_spk_segment* s, myb_time t, double state[6]) {
  const myb_dd offset = myb_time_diff(t, s->init_time);
  const uint64_t index = record_index(s, offset.hi);
  /* t - mid, mid = start + (index + 1/2) step, in double-word arithmetic */
  double low;
  const double high = myb_two_prod((double)index + 0.5, s->step, &low);
  const myb_dd relative = myb_dd_sub(offset, (myb_dd){high, low});
  const double radius = 0.5 * s->step;
  const double x = relative.hi / radius;
  const int n = s->window;
  const uint64_t record = s->begin + index * s->size;
  for (int j = 0; j < 3; ++j) {
    const uint64_t address = record + (uint64_t)(j * (n + 1));
    double c[MYB_SPK_MAX_COEFFICIENTS + 2] = {0};
    for (int k = 0; k < n; ++k) {
      c[k] = word(s, address + (uint64_t)k);
    }
    /* Antiderivative sum C_k T_k with C_1 = c_0 - c_2 / 2, C_k = (c_{k-1} - c_{k+1}) / 2k, and
     * C_0 chosen so that it vanishes at x = 0 (the midpoint). */
    double integral[MYB_SPK_MAX_COEFFICIENTS + 2] = {0};
    integral[1] = c[0] - 0.5 * c[2];
    for (int k = 2; k <= n; ++k) {
      integral[k] = (c[k - 1] - c[k + 1]) / (2.0 * k);
    }
    double at_zero = 0.0;
    for (int k = 2; k <= n; k += 2) {
      at_zero += (k % 4 == 0) ? integral[k] : -integral[k];
    }
    integral[0] = -at_zero;
    double b1 = 0.0;
    double b2 = 0.0;
    double v1 = 0.0;
    double v2 = 0.0;
    for (int k = n; k >= 1; --k) {
      const double b0 = 2.0 * x * b1 - b2 + integral[k];
      b2 = b1;
      b1 = b0;
      if (k < n) {
        const double v0 = 2.0 * x * v1 - v2 + c[k];
        v2 = v1;
        v1 = v0;
      }
    }
    const double primitive = x * b1 - b2 + integral[0];
    const double velocity = x * v1 - v2 + c[0];
    const double midpoint = word(s, address + (uint64_t)n);
    state[j] = (midpoint + primitive * (radius / s->tscale)) * s->dscale;
    state[j + 3] = velocity * (s->dscale / s->tscale);
  }
  return MYB_OK;
}

/* ---- Types 1 and 21: modified difference arrays ------------------------------------------ */

static myb_status parse_mda(myb_internal_spk_segment* s, myb_diag* diag) {
  const uint64_t length = s->end - s->begin + 1;
  const uint64_t trailer = s->type == 1 ? 1 : 2;
  uint64_t n = 0;
  uint64_t dimension = 15;
  if (length < trailer + 2 || !as_count(word(s, s->end), length, &n) || n == 0) {
    return format_error(diag, s, "inconsistent directory");
  }
  if (s->type == 21 && (!as_count(word(s, s->end - 1), MYB_SPK_MAX_MDA, &dimension) ||
                        dimension == 0)) {
    return format_error(diag, s, "difference line size outside 1 .. 25");
  }
  const uint64_t rsize = 4 * dimension + 11;
  if (n > length / (rsize + 1) || n * (rsize + 1) + n / 100 + trailer != length) {
    return format_error(diag, s, "inconsistent directory");
  }
  s->count = n;
  s->size = rsize;
  s->window = (int)dimension;
  s->epochs = s->begin + n * rsize;
  if (!increasing(s, s->epochs, n)) {
    return format_error(diag, s, "epochs not increasing");
  }
  return MYB_OK;
}

static myb_status state_mda(const myb_internal_spk_segment* s, myb_time t, double state[6]) {
  enum { MAXTRM = MYB_SPK_MAX_MDA };

  uint64_t index = count_lt(s, s->epochs, s->count, t); /* first final epoch >= t */
  if (index >= s->count) {
    index = s->count - 1;
  }
  const uint64_t record = s->begin + index * s->size;
  const uint64_t dim = (uint64_t)s->window;
  /* Arrays indexed from 1, as in the description of the algorithm. */
  double g[MAXTRM + 1];
  double dt[3][MAXTRM + 1];
  double fc[MAXTRM + 2];
  double wc[MAXTRM + 1];
  double w[MAXTRM + 3];
  double reference_position[3];
  double reference_velocity[3];
  int kq[3];
  const double tl = word(s, record);
  for (uint64_t j = 1; j <= dim; ++j) {
    g[j] = word(s, record + j);
  }
  for (int i = 0; i < 3; ++i) {
    reference_position[i] = word(s, record + dim + 1 + (uint64_t)(2 * i));
    reference_velocity[i] = word(s, record + dim + 2 + (uint64_t)(2 * i));
    for (uint64_t j = 1; j <= dim; ++j) {
      dt[i][j] = word(s, record + dim + 6 + (uint64_t)i * dim + j);
    }
  }
  uint64_t kqmax1 = 0;
  if (!epoch_ok(tl) || !as_count(word(s, record + 4 * dim + 7), dim + 1, &kqmax1) ||
      kqmax1 < 2) {
    return MYB_E_FORMAT;
  }
  for (int i = 0; i < 3; ++i) {
    uint64_t value = 0;
    if (!as_count(word(s, record + 4 * dim + 8 + (uint64_t)i), kqmax1 - 1, &value)) {
      return MYB_E_FORMAT;
    }
    kq[i] = (int)value;
  }
  const int kmax = (int)kqmax1;
  const double delta = since(t, tl);
  double tp = delta;
  for (int j = 1; j <= kmax - 2; ++j) {
    if (g[j] == 0.0) {
      return MYB_E_FORMAT;
    }
    fc[j + 1] = tp / g[j];
    wc[j] = delta / g[j];
    tp = delta + g[j];
  }
  for (int j = 1; j <= kmax; ++j) {
    w[j] = 1.0 / (double)j;
  }
  int ks = kmax - 1;
  int ks1 = ks - 1;
  int jx = 0;
  /* bound: ks decreases by one per iteration, from kmax - 1 down to 1 */
  while (ks >= 2) {
    ++jx;
    for (int j = 1; j <= jx; ++j) {
      w[j + ks] = fc[j + 1] * w[j + ks1] - wc[j] * w[j + ks];
    }
    ks = ks1;
    --ks1;
  }
  for (int i = 0; i < 3; ++i) {
    double sum = 0.0;
    for (int j = kq[i]; j >= 1; --j) {
      sum += dt[i][j] * w[j + ks];
    }
    state[i] = reference_position[i] + delta * (reference_velocity[i] + delta * sum);
  }
  for (int j = 1; j <= jx; ++j) {
    w[j + ks] = fc[j + 1] * w[j + ks1] - wc[j] * w[j + ks];
  }
  --ks;
  for (int i = 0; i < 3; ++i) {
    double sum = 0.0;
    for (int j = kq[i]; j >= 1; --j) {
      sum += dt[i][j] * w[j + ks];
    }
    state[i + 3] = reference_velocity[i] + delta * sum;
  }
  return MYB_OK;
}

/* ---- Interpolation ------------------------------------------------------------------------- */

/* Lagrange interpolation at 0 of the m points (d_i, y_i) (Neville's algorithm). */
static double lagrange(int m, const double* d, const double* y) {
  double p[MYB_SPK_MAX_WINDOW];
  for (int i = 0; i < m; ++i) {
    p[i] = y[i];
  }
  for (int k = 1; k < m; ++k) {
    for (int i = 0; i + k < m; ++i) {
      p[i] = (d[i + k] * p[i] - d[i] * p[i + 1]) / (d[i + k] - d[i]);
    }
  }
  return p[0];
}

void myb_internal_hermite(int m, const double* d, const double* y, const double* yd,
                          double* value, double* derivative) {
  double z[2 * MYB_SPK_MAX_WINDOW];
  double c[2 * MYB_SPK_MAX_WINDOW];
  const int n = 2 * m;
  for (int i = 0; i < m; ++i) {
    z[2 * i] = d[i];
    z[2 * i + 1] = d[i];
    c[2 * i] = y[i];
    c[2 * i + 1] = y[i];
  }
  for (int k = 1; k < n; ++k) {
    for (int j = n - 1; j >= k; --j) {
      if (k == 1 && (j % 2) == 1) {
        c[j] = yd[j / 2];
      } else {
        c[j] = (c[j] - c[j - 1]) / (z[j] - z[j - k]);
      }
    }
  }
  double p = c[n - 1];
  double dp = 0.0;
  for (int k = n - 2; k >= 0; --k) {
    dp = p - z[k] * dp;
    p = c[k] - z[k] * p;
  }
  *value = p;
  *derivative = dp;
}

/* ---- Types 8, 9, 12, 13: Lagrange and Hermite, window shifted inside the segment ---------- */

static myb_status parse_discrete(myb_internal_spk_segment* s, myb_diag* diag) {
  const uint64_t length = s->end - s->begin + 1;
  const int equal = s->type == 8 || s->type == 12;
  const uint64_t trailer = equal ? 4 : 2;
  uint64_t n = 0;
  uint64_t value = 0;
  if (length < trailer + 6 || !as_count(word(s, s->end), length, &n) || n == 0 ||
      !as_count(word(s, s->end - 1), MYB_SPK_MAX_WINDOW, &value)) {
    return format_error(diag, s, "inconsistent directory");
  }
  /* type 8, 9: polynomial degree; type 12, 13: window size - 1 */
  const uint64_t window = value + 1;
  if (window > MYB_SPK_MAX_WINDOW || window > n) {
    return format_error(diag, s, "window larger than the data or than supported");
  }
  if (equal) {
    s->init = word(s, s->end - 3);
    s->step = word(s, s->end - 2);
    if (!epoch_ok(s->init) || !positive(s->step) || n > (length - 4) / 6 ||
        n * 6 + 4 != length) {
      return format_error(diag, s, "inconsistent directory");
    }
    s->epochs = 0;
  } else {
    if (n > (length - 2) / 7 || n * 7 + (n - 1) / 100 + 2 != length) {
      return format_error(diag, s, "inconsistent directory");
    }
    s->epochs = s->begin + 6 * n;
    if (!increasing(s, s->epochs, n)) {
      return format_error(diag, s, "epochs not increasing");
    }
  }
  s->count = n;
  s->size = 6;
  s->window = (int)window;
  return MYB_OK;
}

/* First state of the window of `size` states centred on t: an even window has t between its
 * two middle epochs, an odd one is centred on the epoch nearest to t; the window is shifted to
 * stay inside the segment. */
static uint64_t centred_window(const myb_internal_spk_segment* s, myb_time t, uint64_t size) {
  const uint64_t n = s->count;
  const uint64_t le = count_le(s, s->epochs, n, t);
  int64_t first;
  if (size % 2 == 0) {
    first = (int64_t)le - (int64_t)(size / 2);
  } else {
    uint64_t near;
    if (le == 0) {
      near = 0;
    } else if (le == n) {
      near = n - 1;
    } else {
      const double below = since(t, epoch_at(s, s->epochs, le - 1));
      const double above = -since(t, epoch_at(s, s->epochs, le));
      near = below < above ? le - 1 : le;
    }
    first = (int64_t)near - (int64_t)(size / 2);
  }
  if (first < 0) {
    first = 0;
  }
  if ((uint64_t)first > n - size) {
    first = (int64_t)(n - size);
  }
  return (uint64_t)first;
}

static myb_status state_discrete(const myb_internal_spk_segment* s, myb_time t,
                                 double state[6]) {
  const int m = s->window;
  const uint64_t first = centred_window(s, t, (uint64_t)m);
  double d[MYB_SPK_MAX_WINDOW];
  for (int i = 0; i < m; ++i) {
    d[i] = -since(t, epoch_at(s, s->epochs, first + (uint64_t)i));
  }
  const uint64_t base = s->begin + first * 6;
  for (int j = 0; j < 3; ++j) {
    double y[MYB_SPK_MAX_WINDOW];
    double yd[MYB_SPK_MAX_WINDOW];
    for (int i = 0; i < m; ++i) {
      y[i] = word(s, base + (uint64_t)(6 * i + j));
      yd[i] = word(s, base + (uint64_t)(6 * i + j + 3));
    }
    if (s->type == 8 || s->type == 9) {
      state[j] = lagrange(m, d, y);
      state[j + 3] = lagrange(m, d, yd);
    } else {
      myb_internal_hermite(m, d, y, yd, &state[j], &state[j + 3]);
    }
  }
  return MYB_OK;
}

/* ---- Types 18 and 19: ESOC/DDID packets, window truncated at the ends --------------------- */

typedef struct mini_segment {
  uint64_t begin;  /* first packet */
  uint64_t epochs; /* first epoch */
  uint64_t count;
  int subtype;
  int window;
  int packet;      /* words per packet */
} mini_segment;

/* Layout of a type 18 segment or a type 19 mini-segment in words [first, last]; with
 * `validate`, also checks that its epochs increase (once, at load time). */
static int parse_packets(const myb_internal_spk_segment* s, uint64_t first, uint64_t last,
                         int with_subtype2, int validate, mini_segment* mini) {
  if (last < first || last - first + 1 < 3 + 7) {
    return 0;
  }
  const uint64_t length = last - first + 1;
  uint64_t subtype = 0;
  uint64_t window = 0;
  uint64_t n = 0;
  if (!as_count(word(s, last - 2), with_subtype2 ? 2 : 1, &subtype) ||
      !as_count(word(s, last - 1), MYB_SPK_MAX_WINDOW, &window) ||
      !as_count(word(s, last), length, &n) || n == 0 || window < 2 || window % 2 != 0) {
    return 0;
  }
  const uint64_t packet = subtype == 0 ? 12 : 6;
  if (n > (length - 3) / (packet + 1) || n * (packet + 1) + (n - 1) / 100 + 3 != length) {
    return 0;
  }
  mini->begin = first;
  mini->epochs = first + n * packet;
  mini->count = n;
  mini->subtype = (int)subtype;
  mini->window = (int)window;
  mini->packet = (int)packet;
  return !validate || increasing(s, mini->epochs, n);
}

static void state_packets(const myb_internal_spk_segment* s, const mini_segment* mini,
                          myb_time t, double state[6]) {
  const uint64_t n = mini->count;
  const int64_t low = (int64_t)count_le(s, mini->epochs, n, t) - 1;
  int64_t first = low - mini->window / 2 + 1;
  int64_t last = low + mini->window / 2;
  if (first < 0) {
    first = 0;
  }
  if (last > (int64_t)n - 1) {
    last = (int64_t)n - 1;
  }
  const int m = (int)(last - first + 1);
  double d[MYB_SPK_MAX_WINDOW];
  for (int i = 0; i < m; ++i) {
    d[i] = -since(t, word(s, mini->epochs + (uint64_t)first + (uint64_t)i));
  }
  const uint64_t base = mini->begin + (uint64_t)first * (uint64_t)mini->packet;
  const uint64_t stride = (uint64_t)mini->packet;
  for (int j = 0; j < 3; ++j) {
    double a[MYB_SPK_MAX_WINDOW];
    double b[MYB_SPK_MAX_WINDOW];
    for (int i = 0; i < m; ++i) {
      a[i] = word(s, base + (uint64_t)i * stride + (uint64_t)j);
      b[i] = word(s, base + (uint64_t)i * stride + (uint64_t)j + 3);
    }
    if (mini->subtype == 0) {
      double c[MYB_SPK_MAX_WINDOW];
      double e[MYB_SPK_MAX_WINDOW];
      for (int i = 0; i < m; ++i) {
        c[i] = word(s, base + (uint64_t)i * stride + (uint64_t)j + 6);
        e[i] = word(s, base + (uint64_t)i * stride + (uint64_t)j + 9);
      }
      double unused;
      myb_internal_hermite(m, d, a, b, &state[j], &unused);
      myb_internal_hermite(m, d, c, e, &state[j + 3], &unused);
    } else if (mini->subtype == 1) {
      state[j] = lagrange(m, d, a);
      state[j + 3] = lagrange(m, d, b);
    } else {
      myb_internal_hermite(m, d, a, b, &state[j], &state[j + 3]);
    }
  }
}

static myb_status parse_type18(myb_internal_spk_segment* s, myb_diag* diag) {
  mini_segment mini;
  if (!parse_packets(s, s->begin, s->end, 0, 1, &mini)) {
    return format_error(diag, s, "inconsistent packets");
  }
  return MYB_OK;
}

/* Type 19: intervals, each with its mini-segment. */
static myb_status parse_type19(myb_internal_spk_segment* s, myb_diag* diag) {
  const uint64_t length = s->end - s->begin + 1;
  uint64_t n = 0;
  /* Boundary selection: 1 selects the later interval; the documented 0 and the -1 written by
   * SPICE's own writer select the earlier one. */
  const double flag = length >= 6 ? word(s, s->end - 1) : 0.0;
  if (length < 6 || !as_count(word(s, s->end), length, &n) || n == 0 ||
      !(flag == 1.0 || flag == 0.0 || flag == -1.0) || n > (length - 4) / 2) {
    return format_error(diag, s, "inconsistent directory");
  }
  const uint64_t pointers = s->end - 1 - (n + 1);
  const double stop = word(s, pointers + n);
  uint64_t first_pointer = 0;
  uint64_t stop_pointer = 0;
  if (!as_count(word(s, pointers), length, &first_pointer) || first_pointer != 1 ||
      !as_count(stop, length, &stop_pointer)) {
    return format_error(diag, s, "inconsistent pointers");
  }
  const uint64_t bounds = s->begin + stop_pointer - 1;
  if (stop_pointer < 2 || bounds + (n + 1) + n / 100 != pointers) {
    return format_error(diag, s, "inconsistent directory");
  }
  if (!increasing(s, bounds, n + 1)) {
    return format_error(diag, s, "interval bounds not increasing");
  }
  uint64_t previous = 1;
  for (uint64_t i = 1; i <= n; ++i) {
    uint64_t next = 0;
    mini_segment mini;
    if (!as_count(word(s, pointers + i), length, &next) || next <= previous ||
        !parse_packets(s, s->begin + previous - 1, s->begin + next - 2, 1, 1, &mini)) {
      return format_error(diag, s, "inconsistent mini-segment");
    }
    previous = next;
  }
  s->count = n;
  s->epochs = bounds;
  s->subtype = flag == 1.0; /* later interval on a common boundary */
  return MYB_OK;
}

static myb_status state_type19(const myb_internal_spk_segment* s, myb_time t, double state[6]) {
  const uint64_t n = s->count;
  const uint64_t le = count_le(s, s->epochs, n + 1, t);
  uint64_t interval = le == 0 ? 0 : le - 1;
  if (interval >= n) {
    interval = n - 1;
  } else if (!s->subtype && interval > 0 &&
             compare_epoch(t, word(s, s->epochs + interval)) == 0) {
    --interval; /* on a common boundary the earlier interval was selected */
  }
  const uint64_t pointers = s->end - 1 - (n + 1);
  uint64_t first = 0;
  uint64_t next = 0;
  mini_segment mini;
  MYB_CHECK(as_count(word(s, pointers + interval), UINT64_MAX >> 12, &first) &&
            as_count(word(s, pointers + interval + 1), UINT64_MAX >> 12, &next) &&
            parse_packets(s, s->begin + first - 1, s->begin + next - 2, 1, 0, &mini));
  state_packets(s, &mini, t, state);
  return MYB_OK;
}

/* ---- Type 5: discrete states with two-body propagation ------------------------------------ */

static myb_status parse_type5(myb_internal_spk_segment* s, myb_diag* diag) {
  const uint64_t length = s->end - s->begin + 1;
  uint64_t n = 0;
  if (length < 2 + 7 || !as_count(word(s, s->end), length, &n) || n == 0 ||
      n > (length - 2) / 7 || n * 7 + n / 100 + 2 != length) {
    return format_error(diag, s, "inconsistent directory");
  }
  s->gm = word(s, s->end - 1);
  if (!positive(s->gm)) {
    return format_error(diag, s, "GM not positive");
  }
  s->count = n;
  s->epochs = s->begin + 6 * n;
  if (!increasing(s, s->epochs, n)) {
    return format_error(diag, s, "epochs not increasing");
  }
  return MYB_OK;
}

static myb_status propagate(const myb_internal_spk_segment* s, uint64_t i, myb_time t,
                            myb_pv* pv) {
  myb_pv start;
  for (int j = 0; j < 3; ++j) {
    start.r[j] = word(s, s->begin + 6 * i + (uint64_t)j);
    start.v[j] = word(s, s->begin + 6 * i + (uint64_t)j + 3);
  }
  return myb_kepler_propagate(s->gm, &start, since(t, word(s, s->epochs + i)), pv);
}

static myb_status state_type5(const myb_internal_spk_segment* s, myb_time t, double state[6]) {
  const uint64_t n = s->count;
  const uint64_t le = count_le(s, s->epochs, n, t);
  myb_pv p1;
  myb_status status;
  if (le == 0 || le == n) {
    status = propagate(s, le == 0 ? 0 : n - 1, t, &p1);
    if (status == MYB_OK) {
      memcpy(state, p1.r, sizeof p1.r);
      memcpy(state + 3, p1.v, sizeof p1.v);
    }
    return status;
  }
  myb_pv p2;
  status = propagate(s, le - 1, t, &p1);
  if (status == MYB_OK) {
    status = propagate(s, le, t, &p2);
  }
  if (status != MYB_OK) {
    return status;
  }
  const double t1 = word(s, s->epochs + le - 1);
  const double t2 = word(s, s->epochs + le);
  const double span = t2 - t1;
  const double argument = since(t, t1) * pi / span;
  const double w = 0.5 + 0.5 * myb_cos(argument);
  const double dw = -0.5 * myb_sin(argument) * (pi / span);
  for (int j = 0; j < 3; ++j) {
    state[j] = w * p1.r[j] + (1.0 - w) * p2.r[j];
    state[j + 3] = w * p1.v[j] + (1.0 - w) * p2.v[j] + dw * (p1.r[j] - p2.r[j]);
  }
  return MYB_OK;
}

/* ---- Type 14: Chebyshev packets of unequal length (generic segment) ----------------------- */

static myb_status parse_type14(myb_internal_spk_segment* s, myb_diag* diag) {
  enum { META = 17 };

  const uint64_t length = s->end - s->begin + 1;
  uint64_t nmeta = 0;
  uint64_t m[META + 1] = {0};
  if (!as_count(word(s, s->end), length, &nmeta) || nmeta < META) {
    return format_error(diag, s, "invalid generic segment meta data");
  }
  for (int i = 1; i <= META; ++i) {
    if (!as_count(word(s, s->end - nmeta + (uint64_t)i), length, &m[i])) {
      return format_error(diag, s, "invalid generic segment meta data");
    }
  }
  const uint64_t refbas = m[6];
  const uint64_t nref = m[7];
  const uint64_t npdr = m[9];
  const uint64_t pktbas = m[11];
  const uint64_t npkt = m[12];
  const uint64_t pktsz = m[15];
  const uint64_t pktoff = m[16];
  if (npkt == 0 || nref != npkt || npdr != 0 || pktsz < 8 || (pktsz - 2) % 6 != 0 ||
      (pktsz - 2) / 6 > MYB_SPK_MAX_COEFFICIENTS || npkt > length / (pktsz + pktoff) ||
      pktbas + npkt * (pktsz + pktoff) > length || refbas + nref > length) {
    return format_error(diag, s, "inconsistent generic segment");
  }
  s->count = npkt;
  s->size = pktsz + pktoff;
  s->offset = (int64_t)pktoff;
  s->data = s->begin + pktbas;
  s->window = (int)((pktsz - 2) / 6);
  s->epochs = s->begin + refbas;
  if (!increasing(s, s->epochs, nref)) {
    return format_error(diag, s, "reference epochs not increasing");
  }
  return MYB_OK;
}

static myb_status state_type14(const myb_internal_spk_segment* s, myb_time t, double state[6]) {
  const uint64_t le = count_le(s, s->epochs, s->count, t);
  const uint64_t index = le == 0 ? 0 : le - 1;
  const uint64_t packet =
      s->data + index * s->size + (uint64_t)s->offset;
  const double mid = word(s, packet);
  const double radius = word(s, packet + 1);
  if (!epoch_ok(mid) || !positive(radius)) {
    return MYB_E_FORMAT;
  }
  const double x = since(t, mid) / radius;
  const int n = s->window;
  for (int j = 0; j < 6; ++j) {
    state[j] = chebyshev(s, packet + 2 + (uint64_t)(j * n), n, x);
  }
  return MYB_OK;
}

/* ---- Type 17: equinoctial elements with precessing node and periapsis --------------------- */

static myb_status parse_type17(myb_internal_spk_segment* s, myb_diag* diag) {
  if (s->end - s->begin + 1 != 12) {
    return format_error(diag, s, "length is not 12");
  }
  for (uint64_t i = 0; i < 12; ++i) {
    if (!isfinite(word(s, s->begin + i))) {
      return format_error(diag, s, "non-finite element");
    }
  }
  const double h = word(s, s->begin + 2);
  const double k = word(s, s->begin + 3);
  if (!epoch_ok(word(s, s->begin)) || !positive(word(s, s->begin + 1)) ||
      !(myb_hypot(h, k) <= 0.9)) {
    return format_error(diag, s, "elements outside the model (a > 0, e <= 0.9)");
  }
  return MYB_OK;
}

static myb_status state_type17(const myb_internal_spk_segment* s, myb_time t, double state[6]) {
  double e[12];
  for (int i = 0; i < 12; ++i) {
    e[i] = word(s, s->begin + (uint64_t)i);
  }
  const double dt = since(t, e[0]);
  const double a = e[1];
  const double rate_periapsis = e[7];
  const double rate_longitude = e[8];
  const double rate_node = e[9];
  /* (h, k) turn with the longitude of periapsis, (p, q) with the node */
  const double cp = myb_cos(rate_periapsis * dt);
  const double sp = myb_sin(rate_periapsis * dt);
  const double h = e[2] * cp + e[3] * sp;
  const double k = e[3] * cp - e[2] * sp;
  const double cn = myb_cos(rate_node * dt);
  const double sn = myb_sin(rate_node * dt);
  const double p = e[5] * cn + e[6] * sn;
  const double q = e[6] * cn - e[5] * sn;
  const double longitude = e[4] + rate_longitude * dt;
  const double hdot = rate_periapsis * k;
  const double kdot = -rate_periapsis * h;
  const double pdot = rate_node * q;
  const double qdot = -rate_node * p;
  /* Kepler's equation in equinoctial form, longitude = F + h cos F - k sin F, through the
   * classical form with the longitude of periapsis. */
  const double eccentricity = myb_hypot(h, k);
  const double varpi = eccentricity > 0.0 ? myb_atan2(h, k) : 0.0;
  double anomaly;
  const myb_status status = myb_kepler_elliptic(longitude - varpi, eccentricity, &anomaly);
  if (status != MYB_OK) {
    return MYB_E_FORMAT;
  }
  const double f_angle = anomaly + varpi;
  const double sf = myb_sin(f_angle);
  const double cf = myb_cos(f_angle);
  const double beta = 1.0 / (1.0 + sqrt((1.0 - eccentricity) * (1.0 + eccentricity)));
  const double fdot = (rate_longitude - hdot * cf + kdot * sf) / (1.0 - h * sf - k * cf);
  const double hk = h * k * beta;
  const double hkdot = (hdot * k + h * kdot) * beta;
  const double x = a * ((1.0 - h * h * beta) * cf + hk * sf - k);
  const double y = a * ((1.0 - k * k * beta) * sf + hk * cf - h);
  const double xdot = a * (-2.0 * h * hdot * beta * cf - (1.0 - h * h * beta) * sf * fdot +
                           hkdot * sf + hk * cf * fdot - kdot);
  const double ydot = a * (-2.0 * k * kdot * beta * sf + (1.0 - k * k * beta) * cf * fdot +
                           hkdot * cf - hk * sf * fdot - hdot);
  const double d = 1.0 / (1.0 + p * p + q * q);
  const double f[3] = {d * (1.0 - p * p + q * q), d * 2.0 * p * q, -d * 2.0 * p};
  const double g[3] = {d * 2.0 * p * q, d * (1.0 + p * p - q * q), d * 2.0 * q};
  const double fd[3] = {d * (-2.0 * p * pdot + 2.0 * q * qdot), d * 2.0 * (pdot * q + p * qdot),
                        -d * 2.0 * pdot};
  const double gd[3] = {d * 2.0 * (pdot * q + p * qdot), d * (2.0 * p * pdot - 2.0 * q * qdot),
                        d * 2.0 * qdot};
  double r[3];
  double v[3];
  for (int i = 0; i < 3; ++i) {
    r[i] = x * f[i] + y * g[i];
    v[i] = xdot * f[i] + x * fd[i] + ydot * g[i] + y * gd[i];
  }
  /* Equatorial frame of the pole (ra, dec): X = Z_inertial x pole, Z = pole. */
  const double cr = myb_cos(e[10]);
  const double sr = myb_sin(e[10]);
  const double cd = myb_cos(e[11]);
  const double sd = myb_sin(e[11]);
  const double ex[3] = {-sr, cr, 0.0};
  const double ey[3] = {-sd * cr, -sd * sr, cd};
  const double ez[3] = {cd * cr, cd * sr, sd};
  for (int i = 0; i < 3; ++i) {
    state[i] = r[0] * ex[i] + r[1] * ey[i] + r[2] * ez[i];
    state[i + 3] = v[0] * ex[i] + v[1] * ey[i] + v[2] * ez[i];
  }
  return MYB_OK;
}

/* ---- Type 15: a conic with precessing node and periapsis ----------------------------------- */

enum { TYPE15_WORDS = 16 };

static myb_status parse_type15(myb_internal_spk_segment* s, myb_diag* diag) {
  if (s->end - s->begin + 1 != TYPE15_WORDS) {
    return format_error(diag, s, "length is not 16");
  }
  double w[TYPE15_WORDS];
  for (uint64_t i = 0; i < TYPE15_WORDS; ++i) {
    w[i] = word(s, s->begin + i);
    if (!isfinite(w[i])) {
      return format_error(diag, s, "non-finite element");
    }
  }
  if (!epoch_ok(w[0]) || !positive(w[7]) || !(w[8] >= 0.0) || !positive(w[13]) ||
      !positive(myb_vec3_norm(w + 1)) || !positive(myb_vec3_norm(w + 4)) ||
      !positive(myb_vec3_norm(w + 10))) {
    return format_error(diag, s, "elements outside the model (p > 0, e >= 0, GM > 0, poles)");
  }
  return MYB_OK;
}

/* v turned about the unit axis by the angle. */
static void turn_about(const double axis[3], double angle, double v[3]) {
  const double c = myb_cos(angle);
  const double sn = myb_sin(angle);
  const double along = myb_vec3_dot(axis, v);
  double across[3];
  myb_vec3_cross(axis, v, across);
  for (int i = 0; i < 3; ++i) {
    v[i] = c * v[i] + sn * across[i] + (1.0 - c) * along * axis[i];
  }
}

/* The conic through the periapsis at the epoch, then the first-order secular effect of the
 * oblateness of the centre: the periapsis turned about the pole of the orbit and the node
 * about the pole of the centre, by angles that grow with the true anomaly gone through (as
 * SPICE does it; with the time they would not follow the motion round an eccentric orbit).
 * The velocity is turned with the position: the rates of the two angles are left out, as in
 * SPICE.
 *
 * On an ellipse the mean anomaly is formed and reduced to one revolution in double-word
 * arithmetic: after thousands of revolutions the last bit of the mean motion is otherwise
 * what the position resolves (CSPICE is 2e-11 off the rule after 3400 revolutions, 1e-13
 * after 50; this evaluation stays at the rounding of one). */
static myb_status state_type15(const myb_internal_spk_segment* s, myb_time t, double state[6]) {
  double w[TYPE15_WORDS];
  for (int i = 0; i < TYPE15_WORDS; ++i) {
    w[i] = word(s, s->begin + (uint64_t)i);
  }
  double orbit_pole[3];
  double periapsis[3];
  double centre_pole[3];
  const double n1 = myb_vec3_norm(w + 1);
  const double n4 = myb_vec3_norm(w + 4);
  const double n10 = myb_vec3_norm(w + 10);
  for (int i = 0; i < 3; ++i) {
    orbit_pole[i] = w[1 + i] / n1;
    periapsis[i] = w[4 + i] / n4;
    centre_pole[i] = w[10 + i] / n10;
  }
  const double p = w[7];
  const double ecc = w[8];
  const double flag = w[9];
  const double gm = w[13];
  const double j2 = w[14];
  const double radius = w[15];
  const double dt = since(t, w[0]);
  double along[3];
  myb_vec3_cross(orbit_pole, periapsis, along);
  myb_pv now;
  double anomaly = 0.0; /* the true anomaly gone through since the periapsis */
  if (ecc < 1.0) {
    /* n = sqrt(gm (1 - e^2)^3 / p^3), M = n dt, both to twice the precision */
    const myb_dd one_me2 = myb_dd_mul(myb_dd_add_d(myb_dd_from_double(-ecc), 1.0),
                                      myb_dd_add_d(myb_dd_from_double(ecc), 1.0));
    const myb_dd ratio = myb_dd_div_d(one_me2, p);
    const myb_dd motion = myb_dd_mul(ratio, myb_dd_sqrt(myb_dd_mul_d(ratio, gm)));
    const myb_dd mean = myb_dd_mul_d(motion, dt);
    const double turns = floor(myb_dd_to_double(mean) / two_pi + 0.5);
    const myb_dd circle = {two_pi, 0x1.1a62633145c07p-52};
    const double part = myb_dd_to_double(myb_dd_sub(mean, myb_dd_mul_d(circle, turns)));
    double e_anomaly = 0.0;
    if (myb_kepler_elliptic(part, ecc, &e_anomaly) != MYB_OK) {
      return MYB_E_FORMAT;
    }
    const double a = p / myb_dd_to_double(one_me2);
    const double b = p / sqrt(myb_dd_to_double(one_me2));
    const double n = myb_dd_to_double(motion);
    const double ce = myb_cos(e_anomaly);
    const double se = myb_sin(e_anomaly);
    const double x = a * (ce - ecc);
    const double y = b * se;
    const double rate = n / (1.0 - ecc * ce);
    const double vx = -a * rate * se;
    const double vy = b * rate * ce;
    for (int i = 0; i < 3; ++i) {
      now.r[i] = x * periapsis[i] + y * along[i];
      now.v[i] = vx * periapsis[i] + vy * along[i];
    }
    anomaly = myb_atan2(y, x) + two_pi * turns;
  } else {
    const double rp = p / (1.0 + ecc);
    const double speed = sqrt(gm / rp) * sqrt(1.0 + ecc);
    myb_pv start;
    for (int i = 0; i < 3; ++i) {
      start.r[i] = rp * periapsis[i];
      start.v[i] = speed * along[i];
    }
    const myb_status status = myb_kepler_propagate(gm, &start, dt, &now);
    if (status != MYB_OK) {
      return status;
    }
  }
  if (flag != 3.0 && j2 != 0.0 && ecc < 1.0) {
    const double cos_inc = myb_vec3_dot(centre_pole, orbit_pole);
    const double z = anomaly * 1.5 * j2 * (radius / p) * (radius / p);
    if (flag != 1.0) {
      turn_about(orbit_pole, z * (2.5 * cos_inc * cos_inc - 0.5), now.r);
      turn_about(orbit_pole, z * (2.5 * cos_inc * cos_inc - 0.5), now.v);
    }
    if (flag != 2.0) {
      turn_about(centre_pole, -z * cos_inc, now.r);
      turn_about(centre_pole, -z * cos_inc, now.v);
    }
  }
  memcpy(state, now.r, sizeof now.r);
  memcpy(state + 3, now.v, sizeof now.v);
  return MYB_OK;
}

/* ---- Dispatch ------------------------------------------------------------------------------ */

myb_status myb_internal_spk_segment_parse(myb_internal_spk_segment* s, myb_diag* diag) {
  MYB_CHECK(s != NULL && s->begin >= 1 && s->begin <= s->end && s->end <= s->daf.words);
  s->supported = 1;
  s->epochs = 0;
  switch (s->type) {
    case 1:
    case 21:
      return parse_mda(s, diag);
    case 2:
      return parse_chebyshev(s, 3, diag);
    case 3:
      return parse_chebyshev(s, 6, diag);
    case 5:
      return parse_type5(s, diag);
    case 8:
    case 9:
    case 12:
    case 13:
      return parse_discrete(s, diag);
    case 14:
      return parse_type14(s, diag);
    case 15:
      return parse_type15(s, diag);
    case 17:
      return parse_type17(s, diag);
    case 18:
      return parse_type18(s, diag);
    case 19:
      return parse_type19(s, diag);
    case 20:
      return parse_type20(s, diag);
    default:
      s->supported = 0;
      return MYB_OK;
  }
}

myb_status myb_internal_spk_segment_state(const myb_internal_spk_segment* s, myb_time t,
                                          double state[6]) {
  MYB_CHECK(s != NULL && state != NULL);
  myb_status status = MYB_E_FORMAT;
  switch (s->supported ? s->type : 0) {
    case 1:
    case 21:
      status = state_mda(s, t, state);
      break;
    case 2:
    case 3:
      status = state_chebyshev(s, t, state);
      break;
    case 5:
      status = state_type5(s, t, state);
      break;
    case 8:
    case 9:
    case 12:
    case 13:
      status = state_discrete(s, t, state);
      break;
    case 14:
      status = state_type14(s, t, state);
      break;
    case 15:
      status = state_type15(s, t, state);
      break;
    case 17:
      status = state_type17(s, t, state);
      break;
    case 18: {
      mini_segment mini;
      MYB_CHECK(parse_packets(s, s->begin, s->end, 0, 0, &mini));
      state_packets(s, &mini, t, state);
      status = MYB_OK;
      break;
    }
    case 19:
      status = state_type19(s, t, state);
      break;
    case 20:
      status = state_type20(s, t, state);
      break;
    default:
      return MYB_E_FORMAT;
  }
  if (status == MYB_OK) {
    for (int i = 0; i < 6; ++i) {
      if (!isfinite(state[i])) {
        return MYB_E_FORMAT;
      }
    }
  }
  return status;
}
