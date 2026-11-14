/* SPDX-License-Identifier: BSD-3-Clause */
/* The fast multipole method of myb_nbody_set_approximation on the threads of a pool
 * (myb_nbody_set_pool):
 *   - the bits: the accelerations with 1, 2, 3 and 4 threads are those of the evaluation without
 *     a pool, bit for bit, and so are the statistics of the walk — for a uniform cube with test
 *     bodies, a Plummer sphere, a cube with a dense clump in a corner, and a tree forced deep
 *     (small leaves, expansions wherever the criterion allows them), at several orders and
 *     opening angles; below the size from which the threads are used the evaluation is the
 *     serial one anyway;
 *   - the cost: microseconds per body for 10^4 and 10^5 bodies at p 8, theta 0.5 with 1, 2 and 4
 *     threads, against the serial FMM and the direct sum on the same threads (printed; the
 *     machine may be shared, so each time is the shortest of a few). */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dyn/fmm.h"
#include "dyn/nbody_internal.h"
#include "mtest.h"
#include "myb/dyn.h"
#include "myb/pool.h"

static const double pi = 3.141592653589793;

enum { MAX_TEAM = 4 };

static double uniform(mt_rng* rng) {
  return (double)(mt_next(rng) >> 11) * 0x1p-53;
}

/* Wall-clock seconds: clock() is the wall clock on Windows (whose C runtime of MinGW lacks
 * timespec_get), but the processor time of all threads on POSIX. */
static double now(void) {
#if defined(_WIN32)
  return (double)clock() / CLOCKS_PER_SEC;
#else
  struct timespec ts;
  MT_CHECK(timespec_get(&ts, TIME_UTC) == TIME_UTC);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
#endif
}

typedef struct cloud {
  size_t massive;
  size_t test;
  double* gm;
  double* x; /* 3 (massive + test) */
} cloud;

static cloud cloud_alloc(size_t massive, size_t test) {
  cloud c = {.massive = massive, .test = test};
  c.gm = malloc(massive * sizeof(double));
  c.x = calloc(3 * (massive + test), sizeof(double));
  MT_CHECK(c.gm != NULL && c.x != NULL);
  return c;
}

static void cloud_free(cloud* c) {
  free(c->gm);
  free(c->x);
}

/* A uniform cube of side 2e11 m, GM 0.5e16 .. 1.5e16 m^3 s^-2; with `clump` > 0 that fraction
 * of the bodies packed into a corner 1000 times smaller. */
static cloud make_cube(size_t massive, size_t test, double clump, uint64_t seed) {
  mt_rng rng = {seed};
  cloud c = cloud_alloc(massive, test);
  for (size_t i = 0; i < massive; ++i) {
    c.gm[i] = 1e16 * (0.5 + uniform(&rng));
  }
  const size_t n = massive + test;
  const size_t packed = (size_t)(clump * (double)n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < 3; ++k) {
      const double u = 2.0 * uniform(&rng) - 1.0;
      c.x[3 * i + k] = i < packed ? 9e10 + 1e8 * u : 1e11 * u;
    }
  }
  return c;
}

/* A Plummer sphere in N-body units (as in test_fmm.c, positions only). */
static cloud make_plummer(size_t n, uint64_t seed) {
  mt_rng rng = {seed};
  cloud c = cloud_alloc(n, 0);
  const double length = 3.0 * pi / 16.0;
  for (size_t i = 0; i < n; ++i) {
    c.gm[i] = 1.0 / (double)n;
    double m = 0.0;
    do {
      m = uniform(&rng);
    } while (m < 1e-10 || m > 0.999);
    const double r = 1.0 / sqrt(pow(m, -2.0 / 3.0) - 1.0);
    const double cz = 2.0 * uniform(&rng) - 1.0;
    const double ph = 2.0 * pi * uniform(&rng);
    const double sz = sqrt(1.0 - cz * cz);
    c.x[3 * i] = length * r * sz * cos(ph);
    c.x[3 * i + 1] = length * r * sz * sin(ph);
    c.x[3 * i + 2] = length * r * cz;
  }
  return c;
}

/* ---- The bits ------------------------------------------------------------------------------ */

typedef struct bits_case {
  const char* name;
  int order;
  double theta;
  size_t leaf;   /* 0: as chosen by the order */
  size_t direct; /* 0: as chosen by the order */
} bits_case;

static void same_bits(myb_pool* pools[MAX_TEAM + 1], const cloud* c, const bits_case* k) {
  const size_t n = c->massive + c->test;
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(c->massive, c->test, c->gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_nbody_set_approximation(s, k->order, k->theta) == MYB_OK);
  myb_internal_fmm_tune(myb_internal_nbody_fmm(s), k->leaf, k->direct);
  double* serial = malloc(3 * n * sizeof(double));
  double* a = malloc(3 * n * sizeof(double));
  MT_CHECK(serial != NULL && a != NULL);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c->x, serial) == MYB_OK);
  const myb_internal_fmm_stats ref = myb_internal_fmm_statistics(myb_internal_nbody_fmm(s));
  size_t differ = 0;
  for (int team = 1; team <= MAX_TEAM; ++team) {
    myb_nbody_set_pool(s, pools[team]);
    for (int repeat = 0; repeat < 2; ++repeat) { /* the second with the lists already grown */
      memset(a, 0, 3 * n * sizeof(double));
      MT_CHECK(myb_nbody_acceleration(s, 0.0, c->x, a) == MYB_OK);
      for (size_t i = 0; i < 3 * n; ++i) {
        differ += mt_bits(a[i]) != mt_bits(serial[i]);
      }
    }
    const myb_internal_fmm_stats st = myb_internal_fmm_statistics(myb_internal_nbody_fmm(s));
    MT_CHECK(st.cells == ref.cells && st.leaves == ref.leaves && st.depth == ref.depth &&
             st.m2l == ref.m2l && st.p2p == ref.p2p);
  }
  MT_CHECK(differ == 0);
  printf(
      "  %-28s N = %6zu (%5zu test), p %2d, theta %.1f: %5zu cells, %7zu M2L, %10zu pairs; "
      "1 .. %d threads: %zu components differ from the serial bits\n",
      k->name, c->massive, c->test, k->order, k->theta, ref.cells, ref.m2l, ref.p2p, MAX_TEAM,
      differ);
  myb_nbody_set_pool(s, NULL);
  free(a);
  free(serial);
  myb_nbody_destroy(s);
}

static void test_bits(myb_pool* pools[MAX_TEAM + 1]) {
  printf("the bits of the threads against the evaluation without a pool:\n");
  cloud cube = make_cube(20000, 3000, 0.0, 5);
  const bits_case c1 = {"uniform cube, test bodies", 8, 0.5, 0, 0};
  same_bits(pools, &cube, &c1);
  const bits_case c2 = {"uniform cube, deep tree", 4, 0.4, 8, 1};
  same_bits(pools, &cube, &c2);
  cloud_free(&cube);
  cloud plummer = make_plummer(20000, 6);
  const bits_case c3 = {"Plummer sphere", 6, 0.7, 0, 0};
  same_bits(pools, &plummer, &c3);
  const bits_case c4 = {"Plummer sphere, p 12", 12, 0.3, 0, 0};
  same_bits(pools, &plummer, &c4);
  cloud_free(&plummer);
  cloud clumped = make_cube(12000, 500, 0.4, 7);
  const bits_case c5 = {"cube with a dense corner", 8, 0.5, 16, 0};
  same_bits(pools, &clumped, &c5);
  cloud_free(&clumped);
  cloud small = make_cube(1500, 100, 0.0, 8); /* below the size of the threads */
  const bits_case c6 = {"small cube (serial anyway)", 8, 0.5, 8, 32};
  same_bits(pools, &small, &c6);
  cloud_free(&small);
}

/* ---- The cost ------------------------------------------------------------------------------ */

/* The shortest of `repeat` evaluations, seconds. */
static double timed(myb_nbody* s, const double* x, double* a, int repeat,
                    myb_internal_nbody_path path) {
  double best = 1e300;
  for (int r = 0; r < repeat; ++r) {
    const double t0 = now();
    MT_CHECK(myb_internal_nbody_acceleration_path(s, x, a, path) == MYB_OK);
    const double t = now() - t0;
    best = t < best ? t : best;
  }
  return best;
}

static void test_cost(myb_pool* pools[MAX_TEAM + 1]) {
  static const int teams[3] = {1, 2, 4};
  printf(
      "cost, p 8, theta 0.5: microseconds per body (wall clock, the shortest of a few); the "
      "direct sum measured for 10^4 bodies and scaled as N^2 for 10^5\n");
  printf(
      "                     FMM serial   FMM 1 thread  2 threads     4 threads   |  direct 1"
      "      2         4 threads\n");
  for (int kind = 0; kind < 2; ++kind) {
    double direct_per_pair[3] = {0.0, 0.0, 0.0};
    for (size_t size = 10000; size <= 100000; size *= 10) {
      cloud c = kind == 0 ? make_cube(size, 0, 0.0, 31) : make_plummer(size, 41);
      const size_t n = c.massive;
      double* a = malloc(3 * n * sizeof(double));
      MT_CHECK(a != NULL);
      myb_nbody* s = NULL;
      MT_CHECK(myb_nbody_create(n, 0, c.gm, NULL, &s) == MYB_OK);
      double direct[3];
      for (int t = 0; t < 3; ++t) {
        if (size == 10000) {
          myb_nbody_set_pool(s, pools[teams[t]]);
          direct[t] = timed(s, c.x, a, 2, MYB_NBODY_PATH_RECEIVER);
          direct_per_pair[t] = direct[t] / ((double)n * (double)n);
        } else {
          direct[t] = direct_per_pair[t] * (double)n * (double)n;
        }
      }
      MT_CHECK(myb_nbody_set_approximation(s, 8, 0.5) == MYB_OK);
      const int repeat = size == 10000 ? 5 : 2;
      myb_nbody_set_pool(s, NULL);
      const double serial = timed(s, c.x, a, repeat, MYB_NBODY_PATH_AUTO);
      double fmm[3];
      for (int t = 0; t < 3; ++t) {
        myb_nbody_set_pool(s, pools[teams[t]]);
        fmm[t] = timed(s, c.x, a, repeat, MYB_NBODY_PATH_AUTO);
      }
      const double per = 1e6 / (double)n;
      printf(
          "  %-7s N = %6zu: %6.2f      %6.2f      %6.2f (%.2fx) %6.2f (%.2fx) | %8.1f  %8.1f  "
          "%8.1f%s\n",
          kind == 0 ? "cube" : "Plummer", n, serial * per, fmm[0] * per, fmm[1] * per,
          serial / fmm[1], fmm[2] * per, serial / fmm[2], direct[0] * per, direct[1] * per,
          direct[2] * per, size == 10000 ? "" : " (scaled)");
      myb_nbody_set_pool(s, NULL);
      myb_nbody_destroy(s);
      free(a);
      cloud_free(&c);
    }
  }
}

int main(void) {
  myb_pool* pools[MAX_TEAM + 1] = {NULL};
  for (int t = 1; t <= MAX_TEAM; ++t) {
    MT_CHECK(myb_pool_create(t, NULL, &pools[t]) == MYB_OK);
  }
  test_bits(pools);
  test_cost(pools);
  for (int t = 1; t <= MAX_TEAM; ++t) {
    myb_pool_destroy(pools[t]);
  }
  return mt_finish("test_fmm_threads");
}
