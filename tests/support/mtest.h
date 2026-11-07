/* SPDX-License-Identifier: BSD-3-Clause */
/* mtest.h — minimal test support: checks with bit-exact comparison, a deterministic
 * random generator and CPU feature detection. Each test is one executable: exit code 0 =
 * pass, 1 = fail, 77 = skipped (CTest SKIP_RETURN_CODE). */
#ifndef MYB_TESTS_MTEST_H_
#define MYB_TESTS_MTEST_H_

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

enum { MTEST_SKIP = 77 };

typedef struct mtest_state {
  long checks;
  long failures;
} mtest_state;

static mtest_state mtest = {};

static inline uint64_t mt_bits(double x) {
  uint64_t u;
  memcpy(&u, &x, sizeof u);
  return u;
}

static inline double mt_from_bits(uint64_t u) {
  double x;
  memcpy(&x, &u, sizeof x);
  return x;
}

/* Same bit pattern, or both NaN (payloads are not part of any contract). */
static inline int mt_same(double a, double b) {
  return mt_bits(a) == mt_bits(b) || (a != a && b != b);
}

static inline void mt_report_failure(const char* file, int line, const char* text) {
  mtest.failures++;
  fprintf(stderr, "%s:%d: check failed: %s\n", file, line, text);
}

#define MT_CHECK(cond)                              \
  do {                                              \
    mtest.checks++;                                 \
    if (!(cond)) {                                  \
      mt_report_failure(__FILE__, __LINE__, #cond); \
    }                                               \
  } while (0)

#define MT_CHECK_SAME(actual, expected)                                                      \
  do {                                                                                       \
    const double mt_a_ = (actual);                                                           \
    const double mt_e_ = (expected);                                                         \
    mtest.checks++;                                                                          \
    if (!mt_same(mt_a_, mt_e_)) {                                                            \
      mtest.failures++;                                                                      \
      fprintf(stderr, "%s:%d: %s = %a (0x%016" PRIx64 "), expected %a (0x%016" PRIx64 ")\n", \
              __FILE__, __LINE__, #actual, mt_a_, mt_bits(mt_a_), mt_e_, mt_bits(mt_e_));    \
    }                                                                                        \
  } while (0)

static inline int mt_finish(const char* name) {
  printf("%s: %ld checks, %ld failures\n", name, mtest.checks, mtest.failures);
  return mtest.failures == 0 ? 0 : 1;
}

/* ---- Deterministic random numbers: SplitMix64 (Steele, Lea, Flood, OOPSLA 2014) ---------- */
typedef struct mt_rng {
  uint64_t state;
} mt_rng;

static inline uint64_t mt_next(mt_rng* rng) {
  uint64_t z = (rng->state += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}

/* Integer uniformly distributed in [lo, hi]. */
static inline int mt_int(mt_rng* rng, int lo, int hi) {
  const uint64_t span = (uint64_t)((int64_t)hi - (int64_t)lo + 1);
  return (int)((int64_t)lo + (int64_t)(mt_next(rng) % span));
}

/* Random sign, random 52-bit fraction, binary exponent uniformly in [emin, emax]
 * (normal numbers; emin >= -1022, emax <= 1023). */
static inline double mt_double(mt_rng* rng, int emin, int emax) {
  const uint64_t fraction = mt_next(rng) >> 12;
  const uint64_t sign = mt_next(rng) & (UINT64_C(1) << 63);
  const uint64_t biased = (uint64_t)(mt_int(rng, emin, emax) + 1023);
  return mt_from_bits(sign | (biased << 52) | fraction);
}

/* Any bit pattern: normals, subnormals, zeros, infinities and NaNs. */
static inline double mt_any_double(mt_rng* rng) {
  return mt_from_bits(mt_next(rng));
}

/* ---- CPU features (x86-64) ------------------------------------------------------------- */
#if defined(__x86_64__)
#include <cpuid.h>

/* FMA3 usable: CPUID.1:ECX.FMA[12] and OSXSAVE[27], and the OS saves YMM state (XCR0). */
static inline int mt_cpu_has_fma(void) {
  unsigned eax = 0;
  unsigned ebx = 0;
  unsigned ecx = 0;
  unsigned edx = 0;
  if (__get_cpuid(1, &eax, &ebx, &ecx, &edx) == 0) {
    return 0;
  }
  if ((ecx & (1u << 12)) == 0 || (ecx & (1u << 27)) == 0) {
    return 0;
  }
  unsigned xcr0_lo = 0;
  unsigned xcr0_hi = 0;
  __asm__("xgetbv" : "=a"(xcr0_lo), "=d"(xcr0_hi) : "c"(0));
  (void)xcr0_hi;
  return (xcr0_lo & 6u) == 6u;
}
#else
static inline int mt_cpu_has_fma(void) {
  return 0;
}
#endif

#endif /* MYB_TESTS_MTEST_H_ */
