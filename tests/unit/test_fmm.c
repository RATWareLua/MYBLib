/* SPDX-License-Identifier: BSD-3-Clause */
/* The fast multipole method of myb_nbody_set_approximation against the exact direct sums:
 *   - refusals and switching: invalid orders and opening angles are refused, order 0 gives
 *     back the exact sums bit for bit, the same input gives the same bits;
 *   - the relative force error against the order p and the opening angle theta for a uniform
 *     cube and a Plummer sphere of N = 100 .. 20000 bodies (the tables are printed): it falls
 *     with p down to the rounding of the direct sums;
 *   - momentum: sum GM_i a_i is zero to rounding (every interaction of two cells is mutual);
 *   - the cost against N up to 10^5 (timings printed; the work per body — pairs summed
 *     directly and interactions of cells — checked to stay bounded, i.e. O(N));
 *   - the graph: a receiver with an excluded source keeps its exact sum, bit for bit; test
 *     bodies; the 1PN terms and a system of bodies on top of the approximation;
 *   - the energy of a Plummer sphere over two crossing times, with the FMM and with the direct
 *     sums. */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dyn/fmm.h"
#include "dyn/nbody_internal.h"
#include "mtest.h"
#include "myb/dyn.h"
#include "myb/integ.h"
#include "myb/pool.h"
#include "myb/rel.h"
#include "myb/system.h"

static const double pi = 3.141592653589793;

static double uniform(mt_rng* rng) {
  return (double)(mt_next(rng) >> 11) * 0x1p-53;
}

static double seconds_since(clock_t t0) {
  return (double)(clock() - t0) / CLOCKS_PER_SEC;
}

/* ---- Clouds of bodies -------------------------------------------------------------------- */

typedef struct cloud {
  size_t massive;
  size_t test;
  double* gm;
  double* x; /* 3 (massive + test) */
  double* v;
} cloud;

static cloud cloud_alloc(size_t massive, size_t test) {
  cloud c = {.massive = massive, .test = test};
  c.gm = malloc(massive * sizeof(double));
  c.x = calloc(3 * (massive + test), sizeof(double));
  c.v = calloc(3 * (massive + test), sizeof(double));
  MT_CHECK(c.gm != NULL && c.x != NULL && c.v != NULL);
  return c;
}

static void cloud_free(cloud* c) {
  free(c->gm);
  free(c->x);
  free(c->v);
}

/* A uniform cube of side 2e11 m, GM 0.5e16 .. 1.5e16 m^3 s^-2, speeds up to 3e4 m/s. */
static cloud make_cube(size_t massive, size_t test, uint64_t seed) {
  mt_rng rng = {seed};
  cloud c = cloud_alloc(massive, test);
  for (size_t i = 0; i < massive; ++i) {
    c.gm[i] = 1e16 * (0.5 + uniform(&rng));
  }
  for (size_t i = 0; i < 3 * (massive + test); ++i) {
    c.x[i] = 1e11 * (2.0 * uniform(&rng) - 1.0);
    c.v[i] = 3e4 * (2.0 * uniform(&rng) - 1.0);
  }
  return c;
}

/* A Plummer sphere in N-body units (G M = 1, virial radius 1, E = -1/4): S. J. Aarseth,
 * M. Henon, R. Wielen, Astron. Astrophys. 37, 183 (1974); the radius is cut where 99.9 % of
 * the mass lies within; the centre of mass at rest at the origin. */
static cloud make_plummer(size_t n, uint64_t seed) {
  mt_rng rng = {seed};
  cloud c = cloud_alloc(n, 0);
  const double length = 3.0 * pi / 16.0;
  const double speed = sqrt(16.0 / (3.0 * pi));
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
    /* the speed by von Neumann's rejection on g(q) = q^2 (1 - q^2)^(7/2) */
    double q = 0.0;
    double g = 0.0;
    do {
      q = uniform(&rng);
      g = 0.1 * uniform(&rng);
    } while (g > q * q * pow(1.0 - q * q, 3.5));
    const double ve = q * sqrt(2.0) * pow(1.0 + r * r, -0.25);
    const double vz = 2.0 * uniform(&rng) - 1.0;
    const double vp = 2.0 * pi * uniform(&rng);
    const double vs = sqrt(1.0 - vz * vz);
    c.v[3 * i] = speed * ve * vs * cos(vp);
    c.v[3 * i + 1] = speed * ve * vs * sin(vp);
    c.v[3 * i + 2] = speed * ve * vz;
  }
  for (size_t k = 0; k < 3; ++k) {
    double sx = 0.0;
    double sv = 0.0;
    for (size_t i = 0; i < n; ++i) {
      sx += c.x[3 * i + k];
      sv += c.v[3 * i + k];
    }
    for (size_t i = 0; i < n; ++i) {
      c.x[3 * i + k] -= sx / (double)n;
      c.v[3 * i + k] -= sv / (double)n;
    }
  }
  return c;
}

/* ---- Measures ---------------------------------------------------------------------------- */

typedef struct error_stats {
  double rms;   /* sqrt(mean (|da| / |a|)^2) */
  double worst; /* max |da| / |a| */
} error_stats;

/* Bodies [begin, end). */
static error_stats compare(size_t begin, size_t end, const double* a, const double* ref) {
  double sum = 0.0;
  double worst = 0.0;
  for (size_t i = begin; i < end; ++i) {
    double d2 = 0.0;
    double r2 = 0.0;
    for (size_t k = 0; k < 3; ++k) {
      const double d = a[3 * i + k] - ref[3 * i + k];
      d2 += d * d;
      r2 += ref[3 * i + k] * ref[3 * i + k];
    }
    const double e = sqrt(d2 / r2);
    sum += e * e;
    worst = e > worst ? e : worst;
  }
  return (error_stats){.rms = sqrt(sum / (double)(end - begin)), .worst = worst};
}

/* |sum GM_i a_i| / sum GM_i |a_i| over the massive bodies. */
static double momentum_balance(size_t n, const double* gm, const double* a) {
  double s[3] = {0.0, 0.0, 0.0};
  double scale = 0.0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < 3; ++k) {
      s[k] += gm[i] * a[3 * i + k];
    }
    scale += gm[i] *
             sqrt(a[3 * i] * a[3 * i] + a[3 * i + 1] * a[3 * i + 1] + a[3 * i + 2] * a[3 * i + 2]);
  }
  return sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]) / scale;
}

/* The exact reference: per receiver, Sum2, on the threads of pool. */
static void direct_reference(myb_nbody* s, myb_pool* pool, const double* x, double* a) {
  myb_nbody_set_pool(s, pool);
  MT_CHECK(myb_internal_nbody_acceleration_path(s, x, a, MYB_NBODY_PATH_RECEIVER) == MYB_OK);
  myb_nbody_set_pool(s, NULL);
}

/* The shortest of `repeat` evaluations (the machine may be shared). */
static double timed(myb_nbody* s, const double* x, double* a, int repeat) {
  double best = 1e300;
  for (int r = 0; r < repeat; ++r) {
    const clock_t t0 = clock();
    MT_CHECK(myb_nbody_acceleration(s, 0.0, x, a) == MYB_OK);
    const double t = seconds_since(t0);
    best = t < best ? t : best;
  }
  return best;
}

/* ---- Refusals, switching, determinism ----------------------------------------------------- */

static void test_switching(void) {
  cloud c = make_cube(300, 20, 3);
  const size_t n = c.massive + c.test;
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(c.massive, c.test, c.gm, NULL, &s) == MYB_OK);
  double theta = 1.0;
  MT_CHECK(myb_nbody_approximation(s, &theta) == 0 && theta == 0.0);
  const int bad_orders[] = {-1, MYB_NBODY_FMM_MAX_ORDER + 1, 1000};
  for (size_t k = 0; k < sizeof bad_orders / sizeof bad_orders[0]; ++k) {
    MT_CHECK(myb_nbody_set_approximation(s, bad_orders[k], 0.5) == MYB_E_DOMAIN);
  }
  const double bad_thetas[] = {0.0, -0.5, 1.0, 1.5, (double)NAN, (double)INFINITY};
  for (size_t k = 0; k < sizeof bad_thetas / sizeof bad_thetas[0]; ++k) {
    MT_CHECK(myb_nbody_set_approximation(s, 4, bad_thetas[k]) == MYB_E_DOMAIN);
  }
  MT_CHECK(myb_nbody_approximation(s, NULL) == 0);
  double* exact = malloc(3 * n * sizeof(double));
  double* a = malloc(3 * n * sizeof(double));
  double* b = malloc(3 * n * sizeof(double));
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, exact) == MYB_OK);

  MT_CHECK(myb_nbody_set_approximation(s, 6, 0.4) == MYB_OK);
  MT_CHECK(myb_nbody_approximation(s, &theta) == 6 && theta == 0.4);
  /* a few levels for 320 bodies, and expansions wherever the criterion allows them */
  myb_internal_fmm_tune(myb_internal_nbody_fmm(s), 8, 1);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, a) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, b) == MYB_OK);
  size_t differ = 0;
  for (size_t i = 0; i < 3 * n; ++i) {
    MT_CHECK_SAME(b[i], a[i]); /* the same input, the same bits */
    differ += mt_bits(a[i]) != mt_bits(exact[i]);
  }
  const myb_internal_fmm_stats st = myb_internal_fmm_statistics(myb_internal_nbody_fmm(s));
  MT_CHECK(st.m2l > 0 && st.depth >= 2 && differ > 0);
  const error_stats e = compare(0, n, a, exact);
  MT_CHECK(e.rms < 1e-4);
  /* a refused change keeps the setting */
  MT_CHECK(myb_nbody_set_approximation(s, 30, 0.5) == MYB_E_DOMAIN);
  MT_CHECK(myb_nbody_approximation(s, NULL) == 6);
  /* the Wisdom-Holman splitting does not take approximate sums */
  {
    const myb_wh_options wo = {.h = 1e-3, .scheme = MYB_WH_DKD};
    myb_wh* wh = NULL;
    MT_CHECK(myb_wh_create(s, &wo, NULL, &wh) == MYB_OK);
    MT_CHECK(myb_wh_start(wh, 0.0, c.x, c.x) == MYB_E_DOMAIN);
    myb_wh_destroy(wh);
  }
  /* order 0: back to the exact sums, bit for bit (theta not read) */
  MT_CHECK(myb_nbody_set_approximation(s, 0, (double)NAN) == MYB_OK);
  MT_CHECK(myb_nbody_approximation(s, &theta) == 0 && theta == 0.0);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, b) == MYB_OK);
  for (size_t i = 0; i < 3 * n; ++i) {
    MT_CHECK_SAME(b[i], exact[i]);
  }
  printf(
      "switching: p 6 theta 0.4, 320 bodies, %zu cells in %zu levels, %zu M2L: rms error "
      "%.1e; order 0 gives back the exact bits\n",
      st.cells, st.depth + 1, st.m2l, e.rms);

  /* a system smaller than a leaf: all pairs summed directly */
  MT_CHECK(myb_nbody_set_approximation(s, 4, 0.5) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, a) == MYB_OK);
  MT_CHECK(myb_nbody_set_approximation(s, 0, 0.0) == MYB_OK);
  myb_nbody* small = NULL;
  MT_CHECK(myb_nbody_create(20, 0, c.gm, NULL, &small) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration(small, 0.0, c.x, exact) == MYB_OK);
  MT_CHECK(myb_nbody_set_approximation(small, 12, 0.5) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration(small, 0.0, c.x, a) == MYB_OK);
  const error_stats e20 = compare(0, 20, a, exact);
  MT_CHECK(e20.worst < 1e-14);
  MT_CHECK(myb_internal_fmm_statistics(myb_internal_nbody_fmm(small)).m2l == 0);
  myb_nbody_destroy(small);

  /* a non-finite position */
  MT_CHECK(myb_nbody_set_approximation(s, 4, 0.5) == MYB_OK);
  c.x[7] = (double)NAN;
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, a) == MYB_E_DOMAIN);
  free(exact);
  free(a);
  free(b);
  myb_nbody_destroy(s);
  cloud_free(&c);
}

/* ---- The error against p and theta -------------------------------------------------------- */

enum { ORDERS = 10, THETAS = 3 };

static const int orders[ORDERS] = {1, 2, 3, 4, 6, 8, 10, 12, 14, 16};
static const double thetas[THETAS] = {0.3, 0.5, 0.7};

static void error_table(myb_pool* pool, const char* name, const cloud* c,
                        double rms[ORDERS][THETAS]) {
  const size_t n = c->massive;
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(n, 0, c->gm, NULL, &s) == MYB_OK);
  double* ref = malloc(3 * n * sizeof(double));
  double* a = malloc(3 * n * sizeof(double));
  direct_reference(s, pool, c->x, ref);
  const double direct = timed(s, c->x, a, 1);
  printf(
      "%s, N = %zu: relative force error, rms and largest (ms per evaluation; the direct "
      "sum %.0f ms)\n   p",
      name, n, 1e3 * direct);
  for (size_t t = 0; t < THETAS; ++t) {
    printf("         theta %.1f        ", thetas[t]);
  }
  printf("\n");
  double momentum = 0.0;
  for (size_t o = 0; o < ORDERS; ++o) {
    printf("  %2d", orders[o]);
    for (size_t t = 0; t < THETAS; ++t) {
      MT_CHECK(myb_nbody_set_approximation(s, orders[o], thetas[t]) == MYB_OK);
      const double time = timed(s, c->x, a, 2);
      const error_stats e = compare(0, n, a, ref);
      rms[o][t] = e.rms;
      const double m = momentum_balance(n, c->gm, a);
      momentum = m > momentum ? m : momentum;
      printf("   %.1e %.1e (%5.0f)", e.rms, e.worst, 1e3 * time);
    }
    printf("\n");
  }
  printf("  momentum: |sum GM a| / sum GM |a| <= %.1e over the table\n", momentum);
  MT_CHECK(momentum < 1e-13);
  free(a);
  free(ref);
  myb_nbody_destroy(s);
}

/* The error falls with p until the rounding of the sums. */
static void check_convergence(double rms[ORDERS][THETAS]) {
  for (size_t t = 0; t < THETAS; ++t) {
    for (size_t o = 1; o < ORDERS; ++o) {
      if (rms[o - 1][t] > 1e-13) {
        MT_CHECK(rms[o][t] < rms[o - 1][t]);
      }
    }
  }
}

static void test_errors(myb_pool* pool) {
  double rms[ORDERS][THETAS];
  cloud cube = make_cube(10000, 0, 9);
  error_table(pool, "uniform cube", &cube, rms);
  check_convergence(rms);
  MT_CHECK(rms[3][1] < 2e-3);  /* p 4, theta 0.5 */
  MT_CHECK(rms[5][1] < 4e-5);  /* p 8 */
  MT_CHECK(rms[7][1] < 2e-6);  /* p 12 */
  MT_CHECK(rms[7][0] < 1e-10); /* p 12, theta 0.3 */
  cloud_free(&cube);
  cloud plummer = make_plummer(10000, 7);
  error_table(pool, "Plummer sphere", &plummer, rms);
  check_convergence(rms);
  MT_CHECK(rms[3][1] < 1e-2);
  MT_CHECK(rms[5][1] < 2e-4);
  MT_CHECK(rms[7][1] < 5e-6);
  MT_CHECK(rms[7][0] < 1e-10);
  cloud_free(&plummer);

  /* the error hardly depends on N */
  const size_t sizes[] = {100, 1000, 20000};
  const int few[] = {4, 8, 12};
  printf("the error against N (theta 0.5), rms for p = 4, 8, 12:\n");
  for (size_t k = 0; k < 3; ++k) {
    for (int kind = 0; kind < 2; ++kind) {
      cloud c = kind == 0 ? make_cube(sizes[k], 0, 11 + k) : make_plummer(sizes[k], 21 + k);
      const size_t n = c.massive;
      myb_nbody* s = NULL;
      MT_CHECK(myb_nbody_create(n, 0, c.gm, NULL, &s) == MYB_OK);
      double* ref = malloc(3 * n * sizeof(double));
      double* a = malloc(3 * n * sizeof(double));
      direct_reference(s, pool, c.x, ref);
      printf("  %-7s N = %5zu:", kind == 0 ? "cube" : "Plummer", n);
      for (size_t o = 0; o < 3; ++o) {
        MT_CHECK(myb_nbody_set_approximation(s, few[o], 0.5) == MYB_OK);
        if (n <= 1000) {
          myb_internal_fmm_tune(myb_internal_nbody_fmm(s), 8, 32); /* small: force the levels */
        }
        MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, a) == MYB_OK);
        const error_stats e = compare(0, n, a, ref);
        printf("  %.1e", e.rms);
        MT_CHECK(e.rms < (few[o] == 4 ? 1e-2 : few[o] == 8 ? 3e-4 : 1e-5));
        MT_CHECK(momentum_balance(n, c.gm, a) < 1e-13);
      }
      printf("\n");
      free(a);
      free(ref);
      myb_nbody_destroy(s);
      cloud_free(&c);
    }
  }
}

/* ---- The cost against N ------------------------------------------------------------------- */

static void test_scaling(void) {
  const size_t sizes[] = {1000, 2000, 5000, 10000, 20000, 50000, 100000};

  enum { SIZES = sizeof sizes / sizeof sizes[0] };

  printf(
      "cost against N, p 8, theta 0.5, one thread: microseconds per body for the FMM and "
      "for the direct sum (measured to 10000, N^2 beyond), the work per body\n");
  for (int kind = 0; kind < 2; ++kind) {
    double direct_per_pair = 0.0;
    double pairs_per_body[SIZES];
    double m2l_per_body[SIZES];
    double fmm_per_body[SIZES];
    for (size_t k = 0; k < SIZES; ++k) {
      const size_t n = sizes[k];
      cloud c = kind == 0 ? make_cube(n, 0, 31 + k) : make_plummer(n, 41 + k);
      myb_nbody* s = NULL;
      MT_CHECK(myb_nbody_create(n, 0, c.gm, NULL, &s) == MYB_OK);
      double* a = malloc(3 * n * sizeof(double));
      double direct = 0.0;
      if (n <= 10000) {
        direct = timed(s, c.x, a, 1);
        direct_per_pair = direct / ((double)n * (double)n);
      } else {
        direct = direct_per_pair * (double)n * (double)n;
      }
      MT_CHECK(myb_nbody_set_approximation(s, 8, 0.5) == MYB_OK);
      const double fmm = timed(s, c.x, a, n <= 20000 ? 3 : 1);
      const myb_internal_fmm_stats st = myb_internal_fmm_statistics(myb_internal_nbody_fmm(s));
      pairs_per_body[k] = (double)st.p2p / (double)n;
      m2l_per_body[k] = (double)st.m2l / (double)n;
      fmm_per_body[k] = fmm / (double)n;
      printf(
          "  %-7s N = %6zu: FMM %6.1f us, direct %7.1f us (%5.1fx); %5.0f pairs and %5.2f "
          "M2L per body, %zu cells, %zu levels\n",
          kind == 0 ? "cube" : "Plummer", n, 1e6 * fmm / (double)n, 1e6 * direct / (double)n,
          direct / fmm, pairs_per_body[k], m2l_per_body[k], st.cells, st.depth + 1);
      free(a);
      myb_nbody_destroy(s);
      cloud_free(&c);
    }
    /* O(N): from 10^4 to 10^5 bodies the work per body stays within a band (it oscillates
     * with the filling of the leaves at each level of the octree) while the direct sum's grows
     * tenfold */
    double lo = pairs_per_body[3];
    double hi = pairs_per_body[3];
    for (size_t k = 3; k < SIZES; ++k) {
      lo = pairs_per_body[k] < lo ? pairs_per_body[k] : lo;
      hi = pairs_per_body[k] > hi ? pairs_per_body[k] : hi;
      MT_CHECK(m2l_per_body[k] < 10.0);
    }
    printf("    10^4 -> 10^5: time per body x%.2f, pairs per body within a factor %.2f\n",
           fmm_per_body[SIZES - 1] / fmm_per_body[3], hi / lo);
    MT_CHECK(hi < 3.0 * lo);
  }
}

/* ---- The graph, test bodies, 1PN, a system ------------------------------------------------- */

static void test_graph(myb_pool* pool) {
  cloud c = make_cube(3000, 500, 51);
  const size_t n = c.massive + c.test;
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(c.massive, c.test, c.gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_nbody_set_interaction(s, 5, 17, MYB_INTERACTION_NONE) == MYB_OK);
  MT_CHECK(myb_nbody_set_interaction(s, 3100, 3, MYB_INTERACTION_NONE) == MYB_OK);
  /* a test body on top of another: they do not interact */
  for (size_t k = 0; k < 3; ++k) {
    c.x[3 * 3200 + k] = c.x[3 * 3201 + k];
  }
  double* ref = malloc(3 * n * sizeof(double));
  double* a = malloc(3 * n * sizeof(double));
  direct_reference(s, pool, c.x, ref);
  MT_CHECK(myb_nbody_set_approximation(s, 8, 0.5) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, a) == MYB_OK);
  for (size_t k = 0; k < 3; ++k) {
    MT_CHECK_SAME(a[3 * 5 + k], ref[3 * 5 + k]);
    MT_CHECK_SAME(a[3 * 3100 + k], ref[3 * 3100 + k]);
  }
  const error_stats em = compare(0, c.massive, a, ref);
  const error_stats et = compare(c.massive, n, a, ref);
  printf(
      "graph: 3000 bodies and 500 test bodies, p 8, theta 0.5: rms error %.1e (massive), "
      "%.1e (test); the two receivers with an excluded source exact to the bit\n",
      em.rms, et.rms);
  MT_CHECK(em.rms < 4e-5 && et.rms < 4e-5);

  /* 1PN on top: the potentials and the 1PN terms stay direct sums */
  const myb_ppn gr = {.c = 299792458.0, .beta = 1.0, .gamma = 1.0};
  MT_CHECK(myb_nbody_set_relativity(s, &gr) == MYB_OK);
  double* newton = malloc(3 * n * sizeof(double));
  myb_nbody_set_pool(s, pool);
  MT_CHECK(myb_internal_nbody_acceleration_xv_path(s, c.x, c.v, ref, MYB_NBODY_PATH_RECEIVER) ==
           MYB_OK);
  myb_nbody_set_pool(s, NULL);
  MT_CHECK(myb_nbody_acceleration_xv(s, 0.0, c.x, c.v, a) == MYB_OK);
  MT_CHECK(myb_nbody_set_relativity(s, NULL) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, newton) == MYB_OK);
  double pn = 0.0;
  double pn_error = 0.0;
  for (size_t i = 0; i < n; ++i) {
    if (i == 5 || i == 3100) {
      continue;
    }
    for (size_t k = 0; k < 3; ++k) {
      const double d = a[3 * i + k] - newton[3 * i + k]; /* the 1PN part with the FMM */
      pn = fabs(d) > pn ? fabs(d) : pn;
    }
  }
  const error_stats epn = compare(0, n, a, ref);
  pn_error = epn.rms;
  printf(
      "  with the 1PN terms: rms error %.1e of the whole acceleration (1PN part up to "
      "%.1e m/s^2)\n",
      pn_error, pn);
  MT_CHECK(pn_error < 4e-5 && pn > 0.0);
  free(newton);
  free(a);
  free(ref);
  myb_nbody_destroy(s);
  cloud_free(&c);

  /* a system of bodies takes the approximation of its point masses */
  cloud d = make_cube(2000, 0, 61);
  myb_system* sys = NULL;
  MT_CHECK(myb_system_create(d.massive, d.gm, NULL, &sys) == MYB_OK);
  const size_t coords = myb_system_coordinates(sys);
  MT_CHECK(coords == 3 * d.massive);
  double* exact = malloc(coords * sizeof(double));
  double* approx = malloc(coords * sizeof(double));
  MT_CHECK(myb_system_acceleration(sys, 0.0, d.x, d.v, exact) == MYB_OK);
  MT_CHECK(myb_nbody_set_approximation(myb_system_nbody(sys), 8, 0.5) == MYB_OK);
  MT_CHECK(myb_system_acceleration(sys, 0.0, d.x, d.v, approx) == MYB_OK);
  const error_stats es = compare(0, d.massive, approx, exact);
  printf("  a system of 2000 bodies: rms error %.1e\n", es.rms);
  MT_CHECK(es.rms < 4e-5 && es.rms > 0.0);
  free(exact);
  free(approx);
  myb_system_destroy(sys);
  cloud_free(&d);
}

/* ---- Energy over a few crossing times -----------------------------------------------------
 * Without softening the energy error of a fixed-step integration of a cluster is that of its
 * close encounters, whatever evaluates the forces. The energy that the approximation itself
 * injects is therefore also measured apart, along the trajectory of the direct integration:
 *   dE = integral of sum_i GM_i v_i . (a_FMM,i - a_i) dt   (G E, as in myb_nbody_invariants),
 * by the trapezoidal rule over every eighth step. In N-body units the dynamical time
 * sqrt(R^3 / GM) is 1 and the crossing time 2 sqrt(2): the run lasts 6 dynamical times. A close
 * pass of two bodies at 1e-4 (an orbit of 1e-4 against the step of 4e-3) spoils the energy of
 * the direct integration by tens of per cent within the first half unit already. */

enum { PROBES = 3 };

typedef struct drift {
  double energy_end; /* |E(T) - E(0)| / |E(0)| */
  double energy_max; /* the largest over the run */
  double momentum;   /* |P(T) - P(0)| / sum GM |v| */
  double seconds;
  double injected[PROBES]; /* |dE| / |E(0)| of each probe */
} drift;

/* sum GM_i v_i . (a_i - b_i) */
static double power(size_t n, const double* gm, const double* v, const double* a, const double* b) {
  double p = 0.0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < 3; ++k) {
      p += gm[i] * v[3 * i + k] * (a[3 * i + k] - b[3 * i + k]);
    }
  }
  return p;
}

/* Yoshida's fourth order with the step h to steps h; probes (NULL or PROBES approximate
 * systems) have their injected energy measured along the way. */
static drift integrate(myb_nbody* s, const cloud* c, double h, int steps, myb_nbody** probes) {
  const size_t n = c->massive;
  myb_symplectic* sym = NULL;
  MT_CHECK(myb_symplectic_create(3 * n, 4, NULL, &sym) == MYB_OK);
  double* y = malloc(6 * n * sizeof(double));
  double* a = malloc(3 * n * sizeof(double));
  double* b = malloc(3 * n * sizeof(double));
  memcpy(y, c->x, 3 * n * sizeof(double));
  memcpy(y + 3 * n, c->v, 3 * n * sizeof(double));
  myb_nbody_invariants i0 = {};
  MT_CHECK(myb_nbody_invariants_of(s, c->x, c->v, &i0) == MYB_OK);
  const double e0 = fabs(myb_dd_to_double(i0.energy));
  double vscale = 0.0;
  for (size_t i = 0; i < n; ++i) {
    vscale += c->gm[i] * sqrt(c->v[3 * i] * c->v[3 * i] + c->v[3 * i + 1] * c->v[3 * i + 1] +
                              c->v[3 * i + 2] * c->v[3 * i + 2]);
  }
  drift out = {};
  double last[PROBES] = {0.0, 0.0, 0.0};
  const int every = 8;
  const clock_t t0 = clock();
  MT_CHECK(myb_symplectic_start(sym, myb_nbody_acceleration, s, 0.0, y) == MYB_OK);
  for (int step = 0; step <= steps; ++step) {
    if (step > 0) {
      MT_CHECK(myb_symplectic_step(sym, h) == MYB_OK);
    }
    const double* state = myb_symplectic_state(sym);
    if (probes != NULL && step % every == 0) {
      MT_CHECK(myb_nbody_acceleration(s, 0.0, state, a) == MYB_OK);
      for (size_t k = 0; k < PROBES; ++k) {
        MT_CHECK(myb_nbody_acceleration(probes[k], 0.0, state, b) == MYB_OK);
        const double p = power(n, c->gm, state + 3 * n, b, a);
        if (step > 0) {
          out.injected[k] += 0.5 * (last[k] + p) * every * h;
        }
        last[k] = p;
      }
    }
    if (step > 0 && step % (steps / 12) == 0) {
      myb_nbody_invariants inv = {};
      MT_CHECK(myb_nbody_invariants_of(s, state, state + 3 * n, &inv) == MYB_OK);
      const double de = fabs(myb_dd_to_double(myb_dd_sub(inv.energy, i0.energy))) / e0;
      out.energy_max = de > out.energy_max ? de : out.energy_max;
      out.energy_end = de;
      double dp = 0.0;
      for (size_t j = 0; j < 3; ++j) {
        const double d = inv.momentum[j] - i0.momentum[j];
        dp += d * d;
      }
      out.momentum = sqrt(dp) / vscale;
    }
  }
  out.seconds = seconds_since(t0);
  for (size_t k = 0; k < PROBES; ++k) {
    out.injected[k] = fabs(out.injected[k]) / e0;
  }
  free(b);
  free(a);
  free(y);
  myb_symplectic_destroy(sym);
  return out;
}

static void test_energy(myb_pool* pool) {
  /* N-body units: the crossing time is 2 sqrt(2) */
  cloud c = make_plummer(512, 71);
  const size_t n = c.massive;
  const double h = 0x1p-8;
  const int steps = 6 * 256;
  const int p[PROBES] = {4, 8, 12};
  myb_nbody* s = NULL;
  myb_nbody* probes[PROBES] = {NULL, NULL, NULL};
  MT_CHECK(myb_nbody_create(n, 0, c.gm, NULL, &s) == MYB_OK);
  for (size_t k = 0; k < PROBES; ++k) {
    MT_CHECK(myb_nbody_create(n, 0, c.gm, NULL, &probes[k]) == MYB_OK);
    MT_CHECK(myb_nbody_set_approximation(probes[k], p[k], 0.5) == MYB_OK);
    /* levels for 512 bodies, and the expansions wherever the criterion allows them */
    myb_internal_fmm_tune(myb_internal_nbody_fmm(probes[k]), 16, 1);
  }
  myb_nbody_set_pool(s, pool);
  const drift exact = integrate(s, &c, h, steps, probes);
  myb_nbody_set_pool(s, NULL);
  printf(
      "Plummer sphere of 512 bodies, Yoshida's 4th order, h = 1/256, over 6 dynamical "
      "(%.1f crossing) times: relative energy error at the end (largest), momentum change\n"
      "  direct sums          %.1e (%.1e), %.1e  [%.1f s, all threads]\n",
      steps * h / (2.0 * sqrt(2.0)), exact.energy_end, exact.energy_max, exact.momentum,
      exact.seconds);
  /* the approximation as it comes (most pairs of 512 bodies are cheaper summed directly) */
  myb_nbody* plain = NULL;
  MT_CHECK(myb_nbody_create(n, 0, c.gm, NULL, &plain) == MYB_OK);
  MT_CHECK(myb_nbody_set_approximation(plain, 8, 0.5) == MYB_OK);
  myb_internal_fmm_tune(myb_internal_nbody_fmm(plain), 16, 0);
  const drift d = integrate(plain, &c, h, steps, NULL);
  printf("  FMM p 8, theta 0.5   %.1e (%.1e), %.1e  [%.1f s, one thread]\n", d.energy_end,
         d.energy_max, d.momentum, d.seconds);
  MT_CHECK(d.momentum < 1e-13);
  MT_CHECK(d.energy_max < 10.0 * exact.energy_max + 1e-6);
  myb_nbody_destroy(plain);
  printf(
      "  energy injected along the direct trajectory by the expansions alone (theta 0.5, "
      "every pair the criterion admits), relative: p 4 %.1e, p 8 %.1e, p 12 %.1e\n",
      exact.injected[0], exact.injected[1], exact.injected[2]);
  MT_CHECK(exact.injected[1] < exact.injected[0] && exact.injected[2] < exact.injected[1]);
  MT_CHECK(exact.injected[1] < 1e-5 && exact.injected[2] < 3e-7);
  for (size_t k = 0; k < PROBES; ++k) {
    myb_nbody_destroy(probes[k]);
  }
  myb_nbody_destroy(s);
  cloud_free(&c);
}

int main(void) {
  myb_pool* pool = NULL;
  MT_CHECK(myb_pool_create(0, NULL, &pool) == MYB_OK);
  test_switching();
  test_errors(pool);
  test_graph(pool);
  test_energy(pool);
  myb_pool_destroy(pool); /* idle workers stay out of the timings below */
  test_scaling();
  return mt_finish("test_fmm");
}
