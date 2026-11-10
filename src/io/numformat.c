/* SPDX-License-Identifier: BSD-3-Clause */
/* Printing of binary64 numbers: the shortest decimal that reads back to the same value, and
 * correctly rounded fixed and scientific forms.
 *
 * Everything rests on one exact operation in integer arithmetic: for value = M 2^E and a
 * decimal shift p, the integer D = RN(|value| 10^p) (ties to even), computed as a quotient of
 * two big integers R / S by generating its decimal digits most significant first (G. L.
 * Steele, J. L. White, "How to print floating-point numbers accurately", PLDI 1990). The
 * shortest form tries 1, 2, ... significant digits and takes the first whose correctly rounded
 * value, or the neighbour on the side of the true value, reads back exactly (myb_parse_double):
 * of all decimals of that length inside the rounding interval it is the closest. */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "io/bigint.h"
#include "myb/io.h"

enum {
  MAX_DIGITS = 400, /* 309 digits before the point and 60 after it, with room to spare */
  MAX_PLACES = 60,
  SHORTEST_DIGITS = 17
};

/* Number of significant bits. */
static int bit_length(const myb_bigint* b) {
  if (b->used == 0) {
    return 0;
  }
  uint64_t top = b->limb[b->used - 1];
  int bits = 0;
  /* bound: 64 bits */
  while (top != 0) {
    top >>= 1;
    ++bits;
  }
  return 64 * (b->used - 1) + bits;
}

/* floor(x * log10(2)) for |x| < 2^20, in integers: 1292913986 / 2^32 = 0.30102999566... */
static int floor_log10_pow2(int x) {
  return (int)(((int64_t)x * 1292913986) >> 32);
}

/* The decimal digits of D = RN(|value| 10^p) for a finite value other than zero, most
 * significant first without leading zeros ("0" for D = 0). Returns their number, or -1 if
 * there are more than MAX_DIGITS. *up: 1 if D lies above |value| 10^p, -1 if below, 0 if
 * equal. */
static int scaled_round(double value, int p, char digits[MAX_DIGITS + 1], int* up) {
  uint64_t bits = 0;
  memcpy(&bits, &value, sizeof bits);
  const uint64_t fraction = bits & ((UINT64_C(1) << 52) - 1);
  const int exponent = (int)((bits >> 52) & 0x7FF);
  const uint64_t m = exponent == 0 ? fraction : (fraction | (UINT64_C(1) << 52));
  const int e = exponent == 0 ? -1074 : exponent - 1075;
  myb_bigint r;
  myb_bigint s;
  myb_bigint_set(&r, m);
  myb_bigint_set(&s, 1);
  if (e > 0) {
    myb_bigint_shift_left(&r, e);
  } else {
    myb_bigint_shift_left(&s, -e);
  }
  if (p > 0) {
    myb_bigint_mul_pow10(&r, p);
  } else {
    myb_bigint_mul_pow10(&s, -p);
  }
  /* nd: the number of digits of floor(R / S); t = S 10^nd > R */
  int nd = 0;
  const int b = bit_length(&r) - bit_length(&s); /* R / S in [2^(b-1), 2^(b+1)) */
  if (b > 1) {
    nd = floor_log10_pow2(b - 1); /* 10^nd <= R / S */
  }
  myb_bigint t = s;
  myb_bigint_mul_pow10(&t, nd);
  /* bound: the estimate is at most two digits short */
  for (int tries = 0; tries < 4 && myb_bigint_compare(&r, &t) >= 0; ++tries) {
    myb_bigint_mul_add(&t, 10, 0);
    ++nd;
  }
  MYB_CHECK(myb_bigint_compare(&r, &t) < 0);
  *up = 0;
  if (nd == 0) {
    /* D is 0 or 1: 2 R against S */
    myb_bigint twice = r;
    myb_bigint_mul_add(&twice, 2, 0);
    const int cmp = myb_bigint_compare(&twice, &s);
    digits[0] = cmp > 0 ? '1' : '0';
    digits[1] = '\0';
    *up = cmp > 0 ? 1 : -1;
    return 1;
  }
  if (nd > MAX_DIGITS) {
    return -1;
  }
  /* t = S 10^(nd - 1): the weight of the first digit */
  t = s;
  myb_bigint_mul_pow10(&t, nd - 1);
  for (int i = 0; i < nd; ++i) {
    int d = 0;
    /* bound: a digit is at most 9 */
    while (d < 9 && myb_bigint_compare(&r, &t) >= 0) {
      myb_bigint_sub(&r, &t);
      ++d;
    }
    MYB_CHECK(myb_bigint_compare(&r, &t) < 0);
    digits[i] = (char)('0' + d);
    myb_bigint_mul_add(&r, 10, 0);
  }
  digits[nd] = '\0';
  /* the rest is R / (10 t): against one half */
  myb_bigint half = t;
  myb_bigint_mul_add(&half, 5, 0);
  const int cmp = myb_bigint_compare(&r, &half);
  const int odd = (digits[nd - 1] - '0') & 1;
  if (cmp > 0 || (cmp == 0 && odd)) {
    *up = 1;
    int i = nd - 1;
    /* bound: the digits */
    while (i >= 0 && digits[i] == '9') {
      digits[i] = '0';
      --i;
    }
    if (i >= 0) {
      digits[i] = (char)(digits[i] + 1);
    } else {
      if (nd + 1 > MAX_DIGITS) {
        return -1;
      }
      memmove(digits + 1, digits, (size_t)nd + 1);
      digits[0] = '1';
      ++nd;
    }
  } else if (r.used != 0) {
    *up = -1;
  }
  return nd;
}

/* floor(log2 |value|) of a finite value other than zero. */
static int binary_exponent(double value) {
  uint64_t bits = 0;
  memcpy(&bits, &value, sizeof bits);
  const uint64_t fraction = bits & ((UINT64_C(1) << 52) - 1);
  const int exponent = (int)((bits >> 52) & 0x7FF);
  if (exponent != 0) {
    return exponent - 1023;
  }
  int top = -1;
  /* bound: 52 bits */
  for (int bit = 51; bit >= 0 && top < 0; --bit) {
    if ((fraction >> bit) & 1u) {
      top = bit;
    }
  }
  return top - 1074;
}

/* `count` significant digits of |value| (finite, not zero), correctly rounded, with the decimal
 * exponent k of the first: |value| ~ d1.d2 ... 10^k. */
static void significant(double value, int count, char digits[MAX_DIGITS + 1], int* k, int* up) {
  int guess = floor_log10_pow2(binary_exponent(value));
  /* bound: the guess is within one of the exponent */
  for (int tries = 0; tries < 4; ++tries) {
    const int n = scaled_round(value, count - 1 - guess, digits, up);
    MYB_CHECK(n > 0);
    if (n == count && digits[0] != '0') {
      *k = guess;
      return;
    }
    if (n == count + 1) {
      int power = digits[0] == '1';
      for (int i = 1; i <= count && power; ++i) {
        power = digits[i] == '0';
      }
      if (power) { /* 10^count: the value itself, or rounded up to it */
        digits[count] = '\0';
        *k = guess + 1;
        return;
      }
    }
    if (n > count) {
      ++guess;
    } else {
      --guess;
    }
  }
  MYB_UNREACHABLE();
}

typedef struct writer {
  char* out;
  size_t size;
  size_t length;
} writer;

static void put(writer* w, char c) {
  if (w->length + 1 < w->size) {
    w->out[w->length] = c;
  }
  w->length++;
}

static void put_text(writer* w, const char* text) {
  for (size_t i = 0; text[i] != '\0'; ++i) {
    put(w, text[i]);
  }
}

static void put_exponent(writer* w, int k) {
  put(w, 'e');
  put(w, k < 0 ? '-' : '+');
  const int a = k < 0 ? -k : k;
  if (a >= 100) {
    put(w, (char)('0' + a / 100));
  }
  put(w, (char)('0' + (a / 10) % 10));
  put(w, (char)('0' + a % 10));
}

static myb_status finish(writer* w, size_t* length) {
  if (length != NULL) {
    *length = w->length;
  }
  if (w->length + 1 > w->size) {
    if (w->size > 0) {
      w->out[w->size - 1] = '\0';
    }
    return MYB_E_OVERFLOW;
  }
  w->out[w->length] = '\0';
  return MYB_OK;
}

/* "nan", "inf", "-inf"; returns 0 for finite values. */
static int put_special(writer* w, double value) {
  if (value != value) {
    put_text(w, "nan");
    return 1;
  }
  if (value - value != 0.0) {
    put_text(w, value < 0.0 ? "-inf" : "inf");
    return 1;
  }
  return 0;
}

static int negative(double value) {
  uint64_t bits = 0;
  memcpy(&bits, &value, sizeof bits);
  return (int)(bits >> 63);
}

/* "d[.ddd]e+XX" from digits and the exponent of the first. */
static void put_scientific(writer* w, const char* digits, int count, int k) {
  put(w, digits[0]);
  if (count > 1) {
    put(w, '.');
    for (int i = 1; i < count; ++i) {
      put(w, digits[i]);
    }
  }
  put_exponent(w, k);
}

size_t myb_format_double(double value, char out[MYB_DOUBLE_CHARS]) {
  MYB_CHECK(out != NULL);
  writer w = {.out = out, .size = MYB_DOUBLE_CHARS};
  if (put_special(&w, value)) {
    (void)finish(&w, NULL);
    return w.length;
  }
  if (negative(value)) {
    put(&w, '-');
  }
  if (value == 0.0) {
    put_text(&w, "0.0");
    (void)finish(&w, NULL);
    return w.length;
  }
  char digits[MAX_DIGITS + 1];
  int count = 0;
  int k = 0;
  /* bound: 17 digits always read back */
  for (int n = 1; n <= SHORTEST_DIGITS && count == 0; ++n) {
    int up = 0;
    significant(value, n, digits, &k, &up);
    /* the correctly rounded digits, then the neighbour on the side of the value */
    for (int attempt = 0; attempt < 2 && count == 0; ++attempt) {
      char text[SHORTEST_DIGITS + 16];
      writer trial = {.out = text, .size = sizeof text};
      put_scientific(&trial, digits, n, k);
      (void)finish(&trial, NULL);
      double back = 0.0;
      size_t used = 0;
      if (myb_parse_double(text, trial.length, &back, &used) == MYB_OK &&
          back == (value < 0.0 ? -value : value)) {
        count = n;
        break;
      }
      if (attempt == 1 || up == 0) {
        break;
      }
      /* one unit of the last digit towards the value, if that keeps the number of digits */
      int i = n - 1;
      if (up > 0) {
        /* bound: the digits */
        while (i >= 0 && digits[i] == '0') {
          digits[i] = '9';
          --i;
        }
        if (i < 0 || (i == 0 && digits[0] == '1' && n > 1)) {
          break;
        }
        digits[i] = (char)(digits[i] - 1);
        if (digits[0] == '0') {
          break;
        }
      } else {
        /* bound: the digits */
        while (i >= 0 && digits[i] == '9') {
          digits[i] = '0';
          --i;
        }
        if (i < 0) {
          break;
        }
        digits[i] = (char)(digits[i] + 1);
      }
    }
  }
  MYB_CHECK(count > 0);
  /* bound: the digits */
  while (count > 1 && digits[count - 1] == '0') {
    --count;
  }
  if (k < -4 || k > 15) {
    put_scientific(&w, digits, count, k);
  } else if (k < 0) {
    put_text(&w, "0.");
    for (int i = 0; i < -k - 1; ++i) {
      put(&w, '0');
    }
    for (int i = 0; i < count; ++i) {
      put(&w, digits[i]);
    }
  } else {
    for (int i = 0; i <= k; ++i) {
      put(&w, i < count ? digits[i] : '0');
    }
    put(&w, '.');
    if (count > k + 1) {
      for (int i = k + 1; i < count; ++i) {
        put(&w, digits[i]);
      }
    } else {
      put(&w, '0');
    }
  }
  (void)finish(&w, NULL);
  return w.length;
}

myb_status myb_format_fixed(double value, int decimals, char* out, size_t size, size_t* length) {
  MYB_CHECK(out != NULL || size == 0);
  if (decimals < 0 || decimals > MAX_PLACES) {
    return MYB_E_DOMAIN;
  }
  writer w = {.out = out, .size = size};
  if (put_special(&w, value)) {
    return finish(&w, length);
  }
  if (negative(value)) {
    put(&w, '-');
  }
  char digits[MAX_DIGITS + 1] = "0";
  int count = 1;
  if (value != 0.0) {
    int up = 0;
    count = scaled_round(value, decimals, digits, &up);
    if (count < 0) {
      return MYB_E_DOMAIN;
    }
  }
  /* the digits are those of an integer: the last `decimals` of them follow the point */
  const int whole = count - decimals;
  if (whole <= 0) {
    put(&w, '0');
  } else {
    for (int i = 0; i < whole; ++i) {
      put(&w, digits[i]);
    }
  }
  if (decimals > 0) {
    put(&w, '.');
    for (int i = whole; i < count; ++i) {
      put(&w, i < 0 ? '0' : digits[i]);
    }
  }
  return finish(&w, length);
}

myb_status myb_format_scientific(double value, int digits_wanted, char* out, size_t size,
                                 size_t* length) {
  MYB_CHECK(out != NULL || size == 0);
  if (digits_wanted < 1 || digits_wanted > MAX_PLACES) {
    return MYB_E_DOMAIN;
  }
  writer w = {.out = out, .size = size};
  if (put_special(&w, value)) {
    return finish(&w, length);
  }
  if (negative(value)) {
    put(&w, '-');
  }
  char digits[MAX_DIGITS + 1];
  int k = 0;
  if (value == 0.0) {
    for (int i = 0; i < digits_wanted; ++i) {
      digits[i] = '0';
    }
    digits[digits_wanted] = '\0';
  } else {
    int up = 0;
    significant(value, digits_wanted, digits, &k, &up);
  }
  put_scientific(&w, digits, digits_wanted, k);
  return finish(&w, length);
}
