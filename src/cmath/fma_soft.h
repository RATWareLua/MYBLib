/* SPDX-License-Identifier: BSD-3-Clause */
/* cmath/fma_soft.h — correctly rounded fused multiply-add in integer arithmetic.
 *
 * myb_internal_fma(a, b, c) = RN(a*b + c), round to nearest, ties to even, for every input
 * including subnormals, signed zeros, infinities and NaN, with IEEE-754 overflow and
 * gradual underflow. It does not depend on hardware FMA (absent from the x86-64-v2
 * baseline) nor on the C runtime (mingw-w64's fma() is not correctly rounded).
 *
 * Method: the exact product of the two 53-bit significands (106 bits) and the significand
 * of c are both normalised to bit 125 of a 128-bit integer, aligned with a sticky bit and
 * added. The low 20 (product) and 73 (c) bits of the normalised values are zero, so an
 * alignment shift that could lose bits implies a distance of more than 20 binades between
 * the operands; then no cancellation can bring the sticky bit near the rounding position,
 * which stays at least 70 bits above it. The rounded integer is written straight into the
 * binary64 encoding, where a carry out of the significand correctly increments the exponent
 * (and produces infinity at the top of the range). Verified bit-for-bit against hardware FMA
 * by tests/unit/test_fma.c. */
#ifndef MYB_SRC_CMATH_FMA_SOFT_H_
#define MYB_SRC_CMATH_FMA_SOFT_H_

#include <stdint.h>
#include <string.h>

typedef unsigned _BitInt(128) myb_u128;

static inline uint64_t myb_fma_bits(double x) {
  uint64_t u;
  memcpy(&u, &x, sizeof u);
  return u;
}

static inline double myb_fma_from_bits(uint64_t u) {
  double x;
  memcpy(&x, &u, sizeof x);
  return x;
}

/* |x| = *significand * 2^*exponent for a finite non-zero encoding (sign ignored). */
static inline void myb_fma_unpack(uint64_t bits, uint64_t* significand, int* exponent) {
  const uint64_t fraction = bits & ((UINT64_C(1) << 52) - 1);
  const int biased = (int)((bits >> 52) & 0x7FF);
  if (biased == 0) {
    *significand = fraction;
    *exponent = -1074;
  } else {
    *significand = fraction | (UINT64_C(1) << 52);
    *exponent = biased - 1075;
  }
}

/* Index of the most significant set bit; x != 0. */
static inline int myb_fma_msb(myb_u128 x) {
  const uint64_t high = (uint64_t)(x >> 64);
  if (high != 0) {
    return 127 - __builtin_clzll(high);
  }
  return 63 - __builtin_clzll((uint64_t)x);
}

/* x >> distance, with every shifted-out non-zero bit ORed into bit 0 (sticky). */
static inline myb_u128 myb_fma_shift_sticky(myb_u128 x, int distance) {
  if (distance == 0) {
    return x;
  }
  if (distance >= 128) {
    return x != 0 ? 1 : 0;
  }
  const myb_u128 lost = x & ((((myb_u128)1) << distance) - 1);
  return (x >> distance) | (lost != 0 ? 1 : 0);
}

/* Rounds (-1)^negative * r * 2^e to binary64; r != 0. */
static inline double myb_fma_round(int negative, myb_u128 r, int e) {
  const uint64_t sign = negative != 0 ? UINT64_C(1) << 63 : 0;
  const int n = myb_fma_msb(r);
  const int result_exponent = n + e; /* r * 2^e lies in [2^result_exponent, 2 * that) */
  if (result_exponent > 1023) {
    return myb_fma_from_bits(sign | UINT64_C(0x7FF0000000000000));
  }
  int shift;     /* bits of r below the last bit of the result */
  uint64_t base; /* encoding offset: (biased exponent - 1) << 52 for normal results */
  if (result_exponent >= -1022) {
    shift = n - 52;
    base = (uint64_t)(result_exponent + 1022) << 52;
  } else {
    shift = -1074 - e;
    base = 0;
  }
  uint64_t q;
  if (shift <= 0) {
    q = (uint64_t)(r << -shift);
  } else if (shift > n + 1) {
    q = 0; /* below half of the smallest subnormal */
  } else if (shift == n + 1) {
    q = r > (((myb_u128)1) << n) ? 1 : 0; /* tie at exactly half rounds to even zero */
  } else {
    q = (uint64_t)(r >> shift);
    const myb_u128 rest = r & ((((myb_u128)1) << shift) - 1);
    const myb_u128 half = ((myb_u128)1) << (shift - 1);
    if (rest > half || (rest == half && (q & 1) != 0)) {
      q += 1;
    }
  }
  return myb_fma_from_bits(sign | (base + q));
}

static inline double myb_internal_fma(double a, double b, double c) {
  const uint64_t a_bits = myb_fma_bits(a);
  const uint64_t b_bits = myb_fma_bits(b);
  const uint64_t c_bits = myb_fma_bits(c);
  const uint64_t exponent_mask = UINT64_C(0x7FF0000000000000);

  /* Infinite or NaN factor: the product is exactly inf or NaN, one rounding suffices. */
  if ((a_bits & exponent_mask) == exponent_mask || (b_bits & exponent_mask) == exponent_mask) {
    return a * b + c;
  }
  /* Finite factors and infinite or NaN addend: the result is the addend (quietened). */
  if ((c_bits & exponent_mask) == exponent_mask) {
    return c + c;
  }
  /* Zero factor: the product is an exact signed zero, one rounding suffices. */
  if (a == 0.0 || b == 0.0) {
    return a * b + c;
  }

  uint64_t a_significand;
  uint64_t b_significand;
  int a_exponent;
  int b_exponent;
  myb_fma_unpack(a_bits, &a_significand, &a_exponent);
  myb_fma_unpack(b_bits, &b_significand, &b_exponent);
  myb_u128 product = (myb_u128)a_significand * (myb_u128)b_significand;
  const int product_shift = 125 - myb_fma_msb(product);
  product <<= product_shift;
  const int product_exponent = a_exponent + b_exponent - product_shift;
  const int product_negative = (int)((a_bits ^ b_bits) >> 63);

  if (c == 0.0) {
    return myb_fma_round(product_negative, product, product_exponent);
  }

  uint64_t c_significand;
  int c_exponent;
  myb_fma_unpack(c_bits, &c_significand, &c_exponent);
  myb_u128 addend = (myb_u128)c_significand;
  const int addend_shift = 125 - myb_fma_msb(addend);
  addend <<= addend_shift;
  const int addend_exponent = c_exponent - addend_shift;
  const int addend_negative = (int)(c_bits >> 63);

  myb_u128 big;
  myb_u128 small;
  int exponent;
  int big_negative;
  int small_negative;
  if (product_exponent >= addend_exponent) {
    big = product;
    big_negative = product_negative;
    small = myb_fma_shift_sticky(addend, product_exponent - addend_exponent);
    small_negative = addend_negative;
    exponent = product_exponent;
  } else {
    big = addend;
    big_negative = addend_negative;
    small = myb_fma_shift_sticky(product, addend_exponent - product_exponent);
    small_negative = product_negative;
    exponent = addend_exponent;
  }

  myb_u128 sum;
  int negative;
  if (big_negative == small_negative) {
    sum = big + small;
    negative = big_negative;
  } else if (big >= small) {
    sum = big - small;
    negative = big_negative;
  } else {
    sum = small - big;
    negative = small_negative;
  }
  if (sum == 0) {
    return 0.0; /* exact cancellation: +0 in round to nearest */
  }
  return myb_fma_round(negative, sum, exponent);
}

#endif /* MYB_SRC_CMATH_FMA_SOFT_H_ */
