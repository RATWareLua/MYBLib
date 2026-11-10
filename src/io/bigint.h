/* SPDX-License-Identifier: BSD-3-Clause */
/* io/bigint.h — non-negative integers of fixed capacity for exact decimal conversion.
 *
 * Capacity: 68 limbs = 4352 bits. The largest operand of the number parser is
 * 10^1125 * (2^54 - 1) < 2^3793 (see io/numparse.c), so the capacity is never exceeded for
 * valid use; exceeding it is a bug and terminates through MYB_CHECK. */
#ifndef MYB_SRC_IO_BIGINT_H_
#define MYB_SRC_IO_BIGINT_H_

#include <stdint.h>

#include "core/base.h"

enum { MYB_BIGINT_LIMBS = 68 };

typedef unsigned _BitInt(128) myb_bigint_wide;

typedef struct myb_bigint {
  uint64_t limb[MYB_BIGINT_LIMBS]; /* little-endian */
  int used;                        /* limbs in use; the top used limb is non-zero */
} myb_bigint;

static inline void myb_bigint_set(myb_bigint* b, uint64_t value) {
  b->limb[0] = value;
  b->used = value != 0 ? 1 : 0;
}

/* b = b * factor + addend. */
static inline void myb_bigint_mul_add(myb_bigint* b, uint64_t factor, uint64_t addend) {
  uint64_t carry = addend;
  for (int i = 0; i < b->used; ++i) {
    const myb_bigint_wide product = (myb_bigint_wide)b->limb[i] * factor + carry;
    b->limb[i] = (uint64_t)product;
    carry = (uint64_t)(product >> 64);
  }
  if (carry != 0) {
    MYB_CHECK(b->used < MYB_BIGINT_LIMBS);
    b->limb[b->used++] = carry;
  }
}

/* b = b * 10^exponent, exponent >= 0. */
static inline void myb_bigint_mul_pow10(myb_bigint* b, int64_t exponent) {
  static const uint64_t pow10[20] = {1u,
                                     10u,
                                     100u,
                                     1000u,
                                     10000u,
                                     100000u,
                                     1000000u,
                                     10000000u,
                                     100000000u,
                                     1000000000u,
                                     10000000000u,
                                     100000000000u,
                                     1000000000000u,
                                     10000000000000u,
                                     100000000000000u,
                                     1000000000000000u,
                                     10000000000000000u,
                                     100000000000000000u,
                                     1000000000000000000u,
                                     10000000000000000000u};
  MYB_CHECK(exponent >= 0 && exponent <= 4 * 64 * MYB_BIGINT_LIMBS);
  /* bound: exponent / 19 iterations, at most 58 for the checked range */
  for (int64_t left = exponent; left > 0; left -= 19) {
    myb_bigint_mul_add(b, pow10[left >= 19 ? 19 : left], 0);
  }
}

/* b = b * 2^bits, bits >= 0. */
static inline void myb_bigint_shift_left(myb_bigint* b, int64_t bits) {
  MYB_CHECK(bits >= 0);
  if (b->used == 0 || bits == 0) {
    return;
  }
  const int limbs = (int)(bits / 64);
  const int shift = (int)(bits % 64);
  MYB_CHECK(bits / 64 < MYB_BIGINT_LIMBS && b->used + limbs + 1 <= MYB_BIGINT_LIMBS);
  int used = b->used + limbs;
  b->limb[used] = 0;
  for (int i = b->used - 1; i >= 0; --i) {
    const uint64_t value = b->limb[i];
    if (shift != 0) {
      b->limb[i + limbs + 1] |= value >> (64 - shift);
      b->limb[i + limbs] = value << shift;
    } else {
      b->limb[i + limbs] = value;
    }
  }
  for (int i = 0; i < limbs; ++i) {
    b->limb[i] = 0;
  }
  if (b->limb[used] != 0) {
    ++used;
  }
  b->used = used;
}

/* a = a - b for a >= b. */
static inline void myb_bigint_sub(myb_bigint* a, const myb_bigint* b) {
  MYB_CHECK(a->used >= b->used);
  uint64_t borrow = 0;
  for (int i = 0; i < a->used; ++i) {
    const uint64_t subtrahend = i < b->used ? b->limb[i] : 0;
    const uint64_t first = a->limb[i] - subtrahend;
    const uint64_t borrow_first = a->limb[i] < subtrahend;
    const uint64_t second = first - borrow;
    const uint64_t borrow_second = first < borrow;
    a->limb[i] = second;
    borrow = borrow_first | borrow_second;
  }
  MYB_CHECK(borrow == 0);
  /* bound: the limbs */
  while (a->used > 0 && a->limb[a->used - 1] == 0) {
    a->used--;
  }
}

/* -1, 0, 1 as a < b, a == b, a > b. */
static inline int myb_bigint_compare(const myb_bigint* a, const myb_bigint* b) {
  if (a->used != b->used) {
    return a->used < b->used ? -1 : 1;
  }
  for (int i = a->used - 1; i >= 0; --i) {
    if (a->limb[i] != b->limb[i]) {
      return a->limb[i] < b->limb[i] ? -1 : 1;
    }
  }
  return 0;
}

#endif /* MYB_SRC_IO_BIGINT_H_ */
