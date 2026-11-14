/* SPDX-License-Identifier: BSD-3-Clause */
/* fuzz.h — what the fuzz harnesses of tests/fuzz share.
 *
 * A harness is one function, LLVMFuzzerTestOneInput, that hands the bytes of one input to a
 * reader of the library and then asks the object read for what it holds, so that the values
 * of a damaged file travel as far as a caller would carry them. A reader fed garbage must
 * return a myb_status: a crash, a read outside the input, a hang or a failed MYB_CHECK is a
 * bug of the reader. What a harness checks beyond that it checks with MYB_CHECK, the
 * library's own fail-fast path, so that libFuzzer and the replay main (replay.c) both see a
 * failure as the end of the process.
 *
 * The same harness links libFuzzer (MYB_FUZZ with Clang on Linux: -fsanitize=fuzzer) or the
 * replay main (every other build: the seeds and a few variants of them run as a test). */
#ifndef MYB_TESTS_FUZZ_FUZZ_H_
#define MYB_TESTS_FUZZ_FUZZ_H_

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "core/check.h"
#include "myb/core.h"
#include "myb/time.h"

/* The entry points libFuzzer calls (and replay.c, which calls them the same way). */
int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size);
int LLVMFuzzerInitialize(int* argc, char*** argv);

/* TDB - TT as the largest term of data/time/tdb-tt-fb1990.tk alone: what the time scales of
 * data/time/time-scales.tk need, for a harness whose input is not that series (the series of
 * 787 terms is 50 KB that every input would parse again). */
#define FUZZ_TDB_SERIES                                                                \
  "\n\\begindata\nMYB_SERIES_TDB_TT_UNIT = 31557600000.0\n"                            \
  "MYB_SERIES_TDB_TT_TERMS = ( 0 1656.674564e-6 6283.075849991 6.240054195 )\n"      \
  "\\begintext\n"

/* A harness that needs no initialisation defines the hook with this. */
#define MYB_FUZZ_NO_INITIALIZE                         \
  int LLVMFuzzerInitialize(int* argc, char*** argv) {  \
    (void)argc;                                        \
    (void)argv;                                        \
    return 0;                                          \
  }

/* An instant of `seconds` past J2000 for a query, or 0 when the number cannot be one: a
 * myb_time takes finite seconds below 2^62 in magnitude (beyond that the call is a contract
 * violation of the harness, not of the reader), and a harness keeps well inside. */
static inline int fuzz_instant(double seconds, myb_time* t) {
  if (!(seconds > -0x1p60 && seconds < 0x1p60)) {
    return 0;
  }
  *t = myb_time_from_seconds(seconds);
  return 1;
}

/* 1 when the n values are all finite: a reader that returns MYB_OK returns numbers, never an
 * infinity or a NaN that a damaged file overflowed into. */
static inline int fuzz_finite(const double* values, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    if (!(values[i] - values[i] == 0.0)) {
      return 0;
    }
  }
  return 1;
}

/* FNV-1a of the input: parameters of a harness (a degree, an instant) are drawn from it, so
 * that the seeds stay the files themselves and every input still reaches several settings. */
static inline uint64_t fuzz_hash(const uint8_t* data, size_t size) {
  uint64_t h = UINT64_C(14695981039346656037);
  for (size_t i = 0; i < size; ++i) {
    h = (h ^ data[i]) * UINT64_C(1099511628211);
  }
  return h;
}

/* The k-th line of the text (0-based; the line ends before '\n', a '\r' before it dropped),
 * copied into out (NUL-terminated, cut to cap - 1 characters); 0 when the text has fewer
 * lines. */
static inline int fuzz_line(const uint8_t* data, size_t size, size_t k, char* out, size_t cap) {
  size_t at = 0;
  size_t line = 0;
  /* bound: one line per iteration, the text is finite */
  while (at < size && line < k) {
    const void* newline = memchr(data + at, '\n', size - at);
    if (newline == NULL) {
      return 0;
    }
    at = (size_t)((const uint8_t*)newline - data) + 1;
    ++line;
  }
  if (at >= size && !(at == 0 && k == 0 && size > 0)) {
    return 0;
  }
  size_t end = at;
  /* bound: the rest of the text */
  while (end < size && data[end] != '\n') {
    ++end;
  }
  if (end > at && data[end - 1] == '\r') {
    --end;
  }
  size_t n = end - at;
  if (n > cap - 1) {
    n = cap - 1;
  }
  memcpy(out, data + at, n);
  out[n] = '\0';
  return 1;
}

/* The integer written right after the first occurrence of `key` in the text (digits, an
 * optional '-'), or `fallback` when there is none: how a harness finds the body or the clock
 * a kernel defines (BODY399_POLE_RA, SCLK01_N_FIELDS_82). */
static inline int64_t fuzz_id_after(const uint8_t* data, size_t size, const char* key,
                                    int64_t fallback) {
  const size_t n = strlen(key);
  if (size < n) {
    return fallback;
  }
  /* bound: every position of the text */
  for (size_t i = 0; i + n <= size; ++i) {
    if (memcmp(data + i, key, n) != 0) {
      continue;
    }
    size_t j = i + n;
    int negative = 0;
    if (j < size && data[j] == '-') {
      negative = 1;
      ++j;
    }
    int64_t value = 0;
    int digits = 0;
    /* bound: at most 15 digits */
    while (j < size && digits < 15 && data[j] >= '0' && data[j] <= '9') {
      value = value * 10 + (data[j] - '0');
      ++digits;
      ++j;
    }
    if (digits > 0) {
      return negative ? -value : value;
    }
  }
  return fallback;
}

#endif /* MYB_TESTS_FUZZ_FUZZ_H_ */
