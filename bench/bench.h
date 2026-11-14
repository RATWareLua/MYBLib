/* SPDX-License-Identifier: BSD-3-Clause */
/* bench.h — the timing harness of bench/: a case is a function that does its work a given
 * number of times; the harness finds how many times make one repeat last at least the
 * minimum time, takes the repeats, and prints one row of a table with the median time of one
 * operation, the fastest and the slowest repeat and their spread around the median.
 *
 * Time comes from the monotonic clock of the platform layer (myb_platform_monotonic_ns).
 * Results vary with the load of the machine: the spread column says by how much, and the
 * median is what to compare. */
#ifndef MYB_BENCH_BENCH_H_
#define MYB_BENCH_BENCH_H_

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "platform/platform.h"

enum { BENCH_MAX_REPEATS = 31 };

typedef struct bench_options {
  int repeats;        /* 1 .. BENCH_MAX_REPEATS */
  double min_seconds; /* the shortest repeat */
  int smoke;          /* the short run of the test: small sizes, few repeats */
  const char* filter; /* only the cases whose group or name contains it; NULL: all */
  int threads;        /* the most threads the parallel cases use */
  const char* data;   /* the directory data/ of the source tree */
} bench_options;

/* The work of a case done `iterations` times; returns 0 on success (a failed case is
 * reported, the others go on). */
typedef int (*bench_body)(void* context, long iterations);

/* "12.3 us" for a time in nanoseconds, into out. */
static inline void bench_format_time(double ns, char out[16]) {
  if (ns < 1e3) {
    snprintf(out, 16, "%.1f ns", ns);
  } else if (ns < 1e6) {
    snprintf(out, 16, "%.2f us", ns / 1e3);
  } else if (ns < 1e9) {
    snprintf(out, 16, "%.2f ms", ns / 1e6);
  } else {
    snprintf(out, 16, "%.3f s", ns / 1e9);
  }
}

static inline int bench_selected(const bench_options* o, const char* group, const char* name) {
  return o->filter == NULL || strstr(group, o->filter) != NULL || strstr(name, o->filter) != NULL;
}

static inline void bench_header(void) {
  printf("| group | case | median | fastest | slowest | spread | n x repeats | note |\n");
  printf("|---|---|---:|---:|---:|---:|---:|---|\n");
}

/* Times body and prints its row, with note (may be NULL) in the last column. Returns the
 * median time of one iteration in nanoseconds, or a negative number if the case failed or is
 * not selected. */
static inline double bench_measure(const bench_options* o, const char* group, const char* name,
                                   bench_body body, void* context, const char* note) {
  if (!bench_selected(o, group, name)) {
    return -1.0;
  }
  /* calibration: double the iterations until one repeat lasts the minimum time */
  long iterations = 1;
  double elapsed = 0.0;
  /* bound: at most 40 doublings */
  for (int k = 0; k < 40; ++k) {
    const uint64_t start = myb_platform_monotonic_ns();
    if (body(context, iterations) != 0) {
      printf("| %s | %s | failed | | | | | |\n", group, name);
      fflush(stdout);
      return -1.0;
    }
    elapsed = (double)(myb_platform_monotonic_ns() - start) * 1e-9;
    if (elapsed >= o->min_seconds || iterations > (1L << 28)) { /* long is 32 bits on Windows */
      break;
    }
    iterations *= elapsed > 0.0 && o->min_seconds / elapsed < 2.0 ? 2 : 4;
  }
  double per[BENCH_MAX_REPEATS];
  const int repeats = o->repeats < 1 ? 1 : o->repeats > BENCH_MAX_REPEATS ? BENCH_MAX_REPEATS
                                                                           : o->repeats;
  for (int r = 0; r < repeats; ++r) {
    const uint64_t start = myb_platform_monotonic_ns();
    if (body(context, iterations) != 0) {
      printf("| %s | %s | failed | | | | | |\n", group, name);
      fflush(stdout);
      return -1.0;
    }
    per[r] = (double)(myb_platform_monotonic_ns() - start) / (double)iterations;
  }
  /* insertion sort of at most BENCH_MAX_REPEATS values */
  for (int i = 1; i < repeats; ++i) {
    const double v = per[i];
    int j = i;
    /* bound: i values */
    while (j > 0 && per[j - 1] > v) {
      per[j] = per[j - 1];
      --j;
    }
    per[j] = v;
  }
  const double median =
      repeats % 2 == 1 ? per[repeats / 2] : 0.5 * (per[repeats / 2 - 1] + per[repeats / 2]);
  char t_median[16];
  char t_fast[16];
  char t_slow[16];
  bench_format_time(median, t_median);
  bench_format_time(per[0], t_fast);
  bench_format_time(per[repeats - 1], t_slow);
  const double spread = median > 0.0 ? 100.0 * (per[repeats - 1] - per[0]) / median : 0.0;
  printf("| %s | %s | %s | %s | %s | %.0f %% | %ld x %d | %s |\n", group, name, t_median, t_fast,
         t_slow, spread, iterations, repeats, note != NULL ? note : "");
  fflush(stdout);
  return median;
}

#endif /* MYB_BENCH_BENCH_H_ */
