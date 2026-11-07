/* SPDX-License-Identifier: BSD-3-Clause */
/* Correctly rounded conversion of decimal and hexadecimal text to binary64, independent of
 * the C locale and of the C runtime's strtod.
 *
 * Decimal method:
 *  1. The significant digits form an integer D (at most 800 digits; any further non-zero
 *     digit sets a sticky flag), the value is V = D * 10^q.
 *  2. Clinger's fast path: D < 10^15 and |q| <= 22 — D and 10^|q| are exact doubles, so a
 *     single IEEE multiplication or division is the correctly rounded result.
 *  3. Otherwise a candidate is computed in double-word arithmetic with a separate binary
 *     exponent (relative error below 10^-17, far below an ulp), placed on the binary64 grid
 *     (normal or subnormal) and decided exactly: V is compared with the midpoints to the
 *     neighbouring doubles using big integers (D * 10^q against K * 2^f, each power moved to
 *     the side where its exponent is non-negative). Ties round to even.
 * 800 digits suffice: a midpoint between two doubles has at most 767 significant digits, so
 * digits beyond the 800th can only break an exact tie upwards (the sticky flag).
 * References: W. D. Clinger, "How to read floating point numbers accurately", PLDI 1990;
 * D. Lemire, "Number parsing at a gigabyte per second", Softw. Pract. Exp. 51(8), 2021
 * (the 768-digit bound); Muller et al., Handbook of Floating-Point Arithmetic, 2nd ed. */
#include <stdint.h>
#include <string.h>

#include "cmath/fma_soft.h"
#include "core/base.h"
#include "io/bigint.h"
#include "myb/io.h"
#include "myb/num.h"

enum {
  MAX_DIGITS = 800,
  EXPONENT_LIMIT = 100000 /* |written exponent| saturates here; the result is 0 or inf */
};

static const uint64_t two_52 = UINT64_C(1) << 52;
static const uint64_t two_53 = UINT64_C(1) << 53;

/* V = D * 10^q, D = digit[0] digit[1] ... digit[count-1] without leading/trailing zeros. */
typedef struct decimal {
  uint8_t digit[MAX_DIGITS];
  int count;
  int64_t q;
  int truncated;
} decimal;

/* ---- Scanning ------------------------------------------------------------------------------- */

static int lower(char c) {
  return (c >= 'A' && c <= 'Z') ? c - 'A' + 'a' : c;
}

static int is_digit(char c) {
  return c >= '0' && c <= '9';
}

static int hex_value(char c) {
  if (is_digit(c)) {
    return c - '0';
  }
  const int l = lower(c);
  return (l >= 'a' && l <= 'f') ? l - 'a' + 10 : -1;
}

/* Case-insensitive match of word at the start of s. */
static int match_word(const char* s, size_t len, const char* word) {
  const size_t n = strlen(word);
  if (len < n) {
    return 0;
  }
  for (size_t i = 0; i < n; ++i) {
    if (lower(s[i]) != word[i]) {
      return 0;
    }
  }
  return 1;
}

static int is_marker(char c, const char* markers) {
  for (const char* m = markers; *m != '\0'; ++m) {
    if (*m == c) {
      return 1;
    }
  }
  return 0;
}

/* Optional "[marker][+-]digits" at s[*pos]: adds its value to *exponent, saturated at
 * EXPONENT_LIMIT. Consumes nothing unless a digit follows the marker (strtod rule). */
static void scan_exponent(const char* s, size_t len, size_t* pos, const char* markers,
                          int64_t* exponent) {
  size_t i = *pos;
  if (i >= len || !is_marker(s[i], markers)) {
    return;
  }
  ++i;
  int negative = 0;
  if (i < len && (s[i] == '+' || s[i] == '-')) {
    negative = s[i] == '-';
    ++i;
  }
  if (i >= len || !is_digit(s[i])) {
    return;
  }
  int64_t value = 0;
  for (; i < len && is_digit(s[i]); ++i) {
    value = value * 10 + (s[i] - '0');
    if (value > EXPONENT_LIMIT) {
      value = EXPONENT_LIMIT;
    }
  }
  *exponent += negative ? -value : value;
  *pos = i;
}

/* Mantissa "digits[.digits]" or ".digits" and an optional exponent with marker e, E, d or D
 * (Fortran style, used by NAIF kernels). Returns characters consumed, 0 without a digit. */
static size_t scan_decimal(const char* s, size_t len, decimal* d) {
  size_t i = 0;
  int seen_digit = 0;
  int seen_point = 0;
  int64_t fraction_digits = 0;
  int64_t dropped_digits = 0;
  d->count = 0;
  d->truncated = 0;
  for (; i < len; ++i) {
    const char c = s[i];
    if (c == '.' && !seen_point) {
      seen_point = 1;
      continue;
    }
    if (!is_digit(c)) {
      break;
    }
    seen_digit = 1;
    if (seen_point) {
      fraction_digits += 1;
    }
    if (d->count == 0 && c == '0') {
      continue; /* leading zero: does not change D */
    }
    if (d->count < MAX_DIGITS) {
      d->digit[d->count++] = (uint8_t)(c - '0');
    } else {
      dropped_digits += 1;
      d->truncated |= c != '0';
    }
  }
  if (!seen_digit) {
    return 0;
  }
  int64_t exponent = 0;
  scan_exponent(s, len, &i, "eEdD", &exponent);
  d->q = exponent - fraction_digits + dropped_digits;
  /* bound: count decreases each iteration */
  while (d->count > 0 && d->digit[d->count - 1] == 0) {
    d->count -= 1;
    d->q += 1;
  }
  return i;
}

/* ---- Double-word value with a separate binary exponent ------------------------------------- */

typedef struct ddx {
  myb_dd m; /* 1 <= m.hi < 2 */
  int64_t e;
} ddx;

static double pow2(int exponent) { /* exact for -1022 <= exponent <= 1023 */
  return myb_fma_from_bits((uint64_t)(exponent + 1023) << 52);
}

static ddx ddx_normalize(myb_dd m, int64_t e) {
  const int k = (int)((myb_fma_bits(m.hi) >> 52) & 0x7FF) - 1023;
  const double scale = pow2(-k);
  return (ddx){.m = {.hi = m.hi * scale, .lo = m.lo * scale}, .e = e + k};
}

static ddx ddx_mul(ddx a, ddx b) {
  return ddx_normalize(myb_dd_mul(a.m, b.m), a.e + b.e);
}

static ddx ddx_div(ddx a, ddx b) {
  return ddx_normalize(myb_dd_div(a.m, b.m), a.e - b.e);
}

static ddx ddx_pow10(int64_t k) {
  ddx result = {.m = {.hi = 1.0, .lo = 0.0}, .e = 0};
  ddx base = {.m = {.hi = 1.25, .lo = 0.0}, .e = 3}; /* 10 = 1.25 * 2^3 */
  /* bound: one iteration per bit of k (k < 2^11 here) */
  for (int64_t left = k; left > 0; left >>= 1) {
    if ((left & 1) != 0) {
      result = ddx_mul(result, base);
    }
    base = ddx_mul(base, base);
  }
  return result;
}

/* ---- Exact decision on the binary64 grid --------------------------------------------------- */

/* A non-negative binary64 value m * 2^e2 in canonical form: normal m in [2^52, 2^53) with
 * e2 >= -1074, subnormal m < 2^52 with e2 = -1074; e2 > 971 stands for infinity. */
typedef struct grid {
  uint64_t m;
  int e2;
} grid;

static grid grid_next_up(grid x) {
  x.m += 1;
  if (x.m == two_53) {
    x.m = two_52;
    x.e2 += 1;
  }
  return x;
}

static grid grid_next_down(grid x) {
  if (x.m == two_52 && x.e2 > -1074) {
    return (grid){.m = two_53 - 1, .e2 = x.e2 - 1};
  }
  x.m -= 1;
  return x;
}

static double grid_to_double(int negative, grid x) {
  const uint64_t sign = negative ? UINT64_C(1) << 63 : 0;
  if (x.e2 > 971) {
    return myb_fma_from_bits(sign | UINT64_C(0x7FF0000000000000));
  }
  if (x.m < two_52) {
    return myb_fma_from_bits(sign | x.m);
  }
  /* biased exponent e2 + 1075, hidden bit removed: ((e2 + 1075) << 52) + m - 2^52 */
  return myb_fma_from_bits(sign | (((uint64_t)(x.e2 + 1074) << 52) + x.m));
}

static grid candidate_from(ddx v) {
  if (v.e > 1023) {
    return (grid){.m = two_53 - 1, .e2 = 971}; /* DBL_MAX; the decision may move to inf */
  }
  if (v.e >= -1022) {
    return (grid){.m = (uint64_t)(v.m.hi * 0x1p52), .e2 = (int)v.e - 52};
  }
  const int64_t shift = v.e + 1074; /* < 52 */
  if (shift < -1) {
    return (grid){.m = 0, .e2 = -1074};
  }
  const double t = v.m.hi * pow2((int)shift); /* < 2^52, exact scaling */
  return (grid){.m = (uint64_t)(t + 0.5), .e2 = -1074};
}

/* sign(D * 10^q + sticky - K * 2^f). */
static int compare_with(const myb_bigint* digits, int64_t q, int truncated, uint64_t k, int64_t f) {
  myb_bigint left = *digits;
  myb_bigint right;
  myb_bigint_set(&right, k);
  if (q >= 0) {
    myb_bigint_mul_pow10(&left, q);
  } else {
    myb_bigint_mul_pow10(&right, -q);
  }
  if (f >= 0) {
    myb_bigint_shift_left(&right, f);
  } else {
    myb_bigint_shift_left(&left, -f);
  }
  const int c = myb_bigint_compare(&left, &right);
  return (c == 0 && truncated) ? 1 : c;
}

static double decide(int negative, const decimal* d) {
  myb_bigint digits;
  myb_bigint_set(&digits, 0);
  uint64_t top = 0; /* the first (up to) 19 digits */
  for (int i = 0; i < d->count; ++i) {
    myb_bigint_mul_add(&digits, 10, d->digit[i]);
    if (i < 19) {
      top = top * 10 + d->digit[i];
    }
  }
  const int used = d->count < 19 ? d->count : 19;
  const int64_t e10 = d->q + (d->count - used); /* V ~ top * 10^e10 */

  const myb_dd top_dd = myb_dd_add_d(myb_dd_from_double((double)(top >> 11) * 2048.0),
                                     (double)(top & 0x7FF)); /* exact */
  ddx approx = ddx_normalize(top_dd, 0);
  approx = e10 >= 0 ? ddx_mul(approx, ddx_pow10(e10)) : ddx_div(approx, ddx_pow10(-e10));

  grid x = candidate_from(approx);
  /* bound: the candidate is at most one unit away; more than 8 steps cannot happen */
  for (int step = 0; step < 8; ++step) {
    if (x.e2 > 971) {
      return grid_to_double(negative, x);
    }
    const int up = compare_with(&digits, d->q, d->truncated, 2 * x.m + 1, (int64_t)x.e2 - 1);
    if (up > 0) {
      x = grid_next_up(x);
      continue;
    }
    if (up == 0) {
      return grid_to_double(negative, (x.m & 1) == 0 ? x : grid_next_up(x));
    }
    if (x.m == 0) {
      return grid_to_double(negative, x);
    }
    const int boundary = x.m == two_52 && x.e2 > -1074;
    const uint64_t k = boundary ? 4 * x.m - 1 : 2 * x.m - 1;
    const int64_t f = boundary ? (int64_t)x.e2 - 2 : (int64_t)x.e2 - 1;
    const int down = compare_with(&digits, d->q, d->truncated, k, f);
    if (down < 0) {
      x = grid_next_down(x);
      continue;
    }
    if (down == 0) {
      return grid_to_double(negative, (x.m & 1) == 0 ? x : grid_next_down(x));
    }
    return grid_to_double(negative, x);
  }
  MYB_UNREACHABLE();
}

static double decimal_to_double(int negative, const decimal* d) {
  const double signed_zero = negative ? -0.0 : 0.0;
  if (d->count == 0) {
    return signed_zero;
  }
  /* V lies in [10^(q + count - 1), 10^(q + count)). */
  const int64_t magnitude = d->q + d->count;
  if (magnitude - 1 >= 309) { /* V >= 1e309 > DBL_MAX + ulp/2 */
    return negative ? -(double)INFINITY : (double)INFINITY;
  }
  if (magnitude <= -324) { /* V < 1e-324 < 2^-1075: rounds to zero */
    return signed_zero;
  }
  static const double exact_pow10[23] = {1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,
                                         1e8,  1e9,  1e10, 1e11, 1e12, 1e13, 1e14, 1e15,
                                         1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
  if (!d->truncated && d->count <= 15 && d->q >= -22 && d->q <= 22) {
    uint64_t value = 0;
    for (int i = 0; i < d->count; ++i) {
      value = value * 10 + d->digit[i];
    }
    const double v = (double)value; /* exact: value < 10^15 < 2^53 */
    const double r = d->q >= 0 ? v * exact_pow10[d->q] : v / exact_pow10[-d->q];
    return negative ? -r : r;
  }
  return decide(negative, d);
}

/* ---- Hexadecimal ---------------------------------------------------------------------------- */

/* "0x" hexdigits [. hexdigits] [p[+-]digits]; s starts at "0x". Returns characters consumed. */
static size_t scan_hex(const char* s, size_t len, int negative, double* out) {
  size_t i = 2;
  uint64_t mantissa = 0;
  int kept = 0; /* significant hex digits in mantissa, at most 16 */
  int sticky = 0;
  int64_t exponent = 0;
  int seen_digit = 0;
  int seen_point = 0;
  for (; i < len; ++i) {
    if (s[i] == '.' && !seen_point) {
      seen_point = 1;
      continue;
    }
    const int v = hex_value(s[i]);
    if (v < 0) {
      break;
    }
    seen_digit = 1;
    if (kept == 0 && v == 0) {
      exponent -= seen_point ? 4 : 0;
      continue;
    }
    if (kept < 16) {
      mantissa = (mantissa << 4) | (uint64_t)v;
      kept += 1;
      exponent -= seen_point ? 4 : 0;
    } else {
      sticky |= v != 0;
      exponent += seen_point ? 0 : 4;
    }
  }
  if (!seen_digit) {
    return 0;
  }
  scan_exponent(s, len, &i, "pP", &exponent);
  if (mantissa == 0) {
    *out = negative ? -0.0 : 0.0;
    return i;
  }
  if (exponent > EXPONENT_LIMIT) {
    exponent = EXPONENT_LIMIT;
  } else if (exponent < -EXPONENT_LIMIT) {
    exponent = -EXPONENT_LIMIT;
  }
  /* 8 spare bits below the mantissa keep the sticky bit under the rounding position. */
  const myb_u128 r = ((myb_u128)mantissa << 8) | (myb_u128)(sticky ? 1 : 0);
  *out = myb_fma_round(negative, r, (int)exponent - 8);
  return i;
}

/* ---- Entry point ----------------------------------------------------------------------------- */

static void finish(size_t used, size_t* consumed) {
  if (consumed != NULL) {
    *consumed = used;
  }
}

myb_status myb_parse_double(const char* text, size_t length, double* value, size_t* consumed) {
  MYB_CHECK(text != NULL || length == 0);
  MYB_CHECK(value != NULL);
  size_t i = 0;
  int negative = 0;
  if (i < length && (text[i] == '+' || text[i] == '-')) {
    negative = text[i] == '-';
    ++i;
  }
  const char* s = text + i;
  const size_t rest = length - i;
  if (match_word(s, rest, "infinity") || match_word(s, rest, "inf")) {
    *value = negative ? -(double)INFINITY : (double)INFINITY;
    finish(i + (match_word(s, rest, "infinity") ? 8 : 3), consumed);
    return MYB_OK;
  }
  if (match_word(s, rest, "nan")) {
    *value = negative ? -(double)NAN : (double)NAN;
    finish(i + 3, consumed);
    return MYB_OK;
  }
  if (rest >= 3 && s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {
    double result = 0.0;
    const size_t used = scan_hex(s, rest, negative, &result);
    if (used != 0) {
      *value = result;
      finish(i + used, consumed);
      return MYB_OK;
    }
  }
  decimal d;
  const size_t used = scan_decimal(s, rest, &d);
  if (used == 0) {
    return MYB_E_FORMAT;
  }
  *value = decimal_to_double(negative, &d);
  finish(i + used, consumed);
  return MYB_OK;
}
