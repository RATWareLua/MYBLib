/* SPDX-License-Identifier: BSD-3-Clause */
/* N-body system: accelerations against a double-double reference, bit-identical results of the
 * symmetric and per-receiver paths for any number of threads, the interaction graph, the
 * invariants, and classical solutions:
 *   - Lagrange's equilateral solution (rigid rotation, exact),
 *   - the figure-eight choreography (Chenciner & Montgomery 2000; C. Moore 1993), its initial
 *     condition refined here by Gauss-Newton shooting from the 6-digit values of Suvakov &
 *     Dmitrasinovic (PRL 110, 114301, 2013) and compared with Simo's,
 *   - Burrau's Pythagorean problem (Szebehely & Peters, AJ 72, 876, 1967): a binary of the
 *     masses 4 and 5, the mass 3 escapes. */
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "core/fpguard.h"
#include "dyn/nbody_internal.h"
#include "mtest.h"
#include "myb/dyn.h"
#include "myb/integ.h"
#include "myb/num.h"
#include "myb/pool.h"
#include "num/eft.h"

static const double pi = 3.141592653589793;

static double uniform(mt_rng* rng) {
  return (double)(mt_next(rng) >> 11) * 0x1p-53;
}

/* ---- Accelerations --------------------------------------------------------------------- */

typedef struct cluster {
  size_t massive;
  size_t test;
  double* gm;
  double* x;
} cluster;

static cluster make_cluster(size_t massive, size_t test, uint64_t seed) {
  mt_rng rng = {seed};
  cluster c = {.massive = massive, .test = test};
  c.gm = malloc(massive * sizeof(double));
  c.x = malloc(3 * (massive + test) * sizeof(double));
  MT_CHECK(c.gm != NULL && c.x != NULL);
  for (size_t i = 0; i < massive; ++i) {
    c.gm[i] = pow(10.0, 5.0 + 15.0 * uniform(&rng)); /* 1e5 .. 1e20 m^3 s^-2 */
  }
  for (size_t i = 0; i < 3 * (massive + test); ++i) {
    c.x[i] = 1e11 * (2.0 * uniform(&rng) - 1.0);
  }
  return c;
}

static void free_cluster(cluster* c) {
  free(c->gm);
  free(c->x);
}

/* Reference: every term and the sum in double-double. Also the sum of |terms| per component
 * (the scale of the rounding errors that any summation must make in the terms). */
static void reference(const cluster* c, size_t i, double out[3], double scale[3]) {
  myb_dd sum[3] = {myb_dd_from_double(0.0), myb_dd_from_double(0.0), myb_dd_from_double(0.0)};
  for (int k = 0; k < 3; ++k) {
    scale[k] = 0.0;
  }
  const double* xi = c->x + 3 * i;
  for (size_t j = 0; j < c->massive; ++j) {
    if (j == i) {
      continue;
    }
    const double* xj = c->x + 3 * j;
    myb_dd d[3];
    for (int k = 0; k < 3; ++k) {
      double err = 0.0;
      d[k].hi = myb_two_sum(xj[k], -xi[k], &err);
      d[k].lo = err;
    }
    const myb_dd r2 =
        myb_dd_add(myb_dd_add(myb_dd_mul(d[0], d[0]), myb_dd_mul(d[1], d[1])), myb_dd_mul(d[2], d[2]));
    const myb_dd r3 = myb_dd_mul(r2, myb_dd_sqrt(r2));
    for (int k = 0; k < 3; ++k) {
      const myb_dd term = myb_dd_div(myb_dd_mul_d(d[k], c->gm[j]), r3);
      sum[k] = myb_dd_add(sum[k], term);
      scale[k] += fabs(term.hi);
    }
  }
  for (int k = 0; k < 3; ++k) {
    out[k] = myb_dd_to_double(sum[k]);
  }
}

static void test_accuracy(void) {
  cluster c = make_cluster(400, 50, 11);
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(c.massive, c.test, c.gm, NULL, &s) == MYB_OK);
  const size_t n = c.massive + c.test;
  double* a = malloc(3 * n * sizeof(double));
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, a) == MYB_OK);
  double worst = 0.0;
  double worst_plain = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double ref[3];
    double scale[3];
    reference(&c, i, ref, scale);
    /* plain summation of the same terms, for comparison */
    double plain[3] = {0.0, 0.0, 0.0};
    for (size_t j = 0; j < c.massive; ++j) {
      if (j == i) {
        continue;
      }
      const double dx = c.x[3 * j] - c.x[3 * i];
      const double dy = c.x[3 * j + 1] - c.x[3 * i + 1];
      const double dz = c.x[3 * j + 2] - c.x[3 * i + 2];
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double sj = c.gm[j] / (r2 * sqrt(r2));
      plain[0] += sj * dx;
      plain[1] += sj * dy;
      plain[2] += sj * dz;
    }
    for (int k = 0; k < 3; ++k) {
      worst = fmax(worst, fabs(a[3 * i + (size_t)k] - ref[k]) / scale[k]);
      worst_plain = fmax(worst_plain, fabs(plain[k] - ref[k]) / scale[k]);
    }
  }
  printf("accelerations (400 + 50 bodies, GM 1e5..1e20): error / sum|terms| %.2e "
         "(plain summation %.2e)\n",
         worst, worst_plain);
  MT_CHECK(worst < 4e-16);
  free(a);
  myb_nbody_destroy(s);
  free_cluster(&c);
}

static int same_arrays(const double* a, const double* b, size_t count) {
  for (size_t i = 0; i < count; ++i) {
    if (!mt_same(a[i], b[i])) {
      return 0;
    }
  }
  return 1;
}

static void test_paths_and_threads(void) {
  cluster c = make_cluster(700, 33, 12);
  const size_t n = c.massive + c.test;
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(c.massive, c.test, c.gm, NULL, &s) == MYB_OK);
  double* sym = malloc(3 * n * sizeof(double));
  double* other = malloc(3 * n * sizeof(double));
  MT_CHECK(myb_internal_nbody_acceleration_path(s, c.x, sym, MYB_NBODY_PATH_SYMMETRIC) == MYB_OK);
  MT_CHECK(myb_internal_nbody_acceleration_path(s, c.x, other, MYB_NBODY_PATH_RECEIVER) ==
           MYB_OK);
  MT_CHECK(same_arrays(sym, other, 3 * n));
  const int threads[] = {1, 2, 3, 4, 8, 13};
  int identical = 1;
  for (size_t t = 0; t < sizeof threads / sizeof threads[0]; ++t) {
    myb_pool* pool = NULL;
    MT_CHECK(myb_pool_create(threads[t], NULL, &pool) == MYB_OK);
    myb_nbody_set_pool(s, pool);
    MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, other) == MYB_OK);
    identical = identical && same_arrays(sym, other, 3 * n);
    myb_nbody_set_pool(s, NULL);
    myb_pool_destroy(pool);
  }
  MT_CHECK(identical);
  /* Throughput (printed only). */
  myb_pool* pool = NULL;
  MT_CHECK(myb_pool_create(0, NULL, &pool) == MYB_OK);
  const clock_t t0 = clock();
  for (int r = 0; r < 20; ++r) {
    MT_CHECK(myb_internal_nbody_acceleration_path(s, c.x, sym, MYB_NBODY_PATH_SYMMETRIC) ==
             MYB_OK);
  }
  const clock_t t1 = clock();
  myb_nbody_set_pool(s, pool);
  for (int r = 0; r < 20; ++r) {
    MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, other) == MYB_OK);
  }
  const clock_t t2 = clock();
  const double pairs = (double)c.massive * (double)(c.massive - 1) / 2.0;
  printf("700 + 33 bodies: symmetric %.1f ns per pair; %d threads %.2f ms per evaluation "
         "(symmetric %.2f ms); identical for 1..13 threads: %s\n",
         1e9 * (double)(t1 - t0) / CLOCKS_PER_SEC / 20.0 / pairs, myb_pool_threads(pool),
         1e3 * (double)(t2 - t1) / CLOCKS_PER_SEC / 20.0,
         1e3 * (double)(t1 - t0) / CLOCKS_PER_SEC / 20.0, identical ? "yes" : "NO");
  myb_nbody_set_pool(s, NULL);
  myb_pool_destroy(pool);
  free(sym);
  free(other);
  myb_nbody_destroy(s);
  free_cluster(&c);
}

/* Vector blocks of four receivers against the scalar receiver path and the symmetric path,
 * for every size up to 17 massive and 6 test bodies (all tails, blocks straddling the end of
 * the massive bodies), GM = 0 sources, and parts of the parallel loop that start anywhere. */
static void test_vector_paths(void) {
  int identical = 1;
  myb_pool* pools[3] = {NULL, NULL, NULL};
  MT_CHECK(myb_pool_create(3, NULL, &pools[1]) == MYB_OK);
  MT_CHECK(myb_pool_create(5, NULL, &pools[2]) == MYB_OK);
  for (size_t m = 1; m <= 17; ++m) {
    for (size_t t = 0; t <= 6; ++t) {
      cluster c = make_cluster(m, t, 1000 * m + t);
      if (m > 2) {
        c.gm[m / 2] = 0.0;
      }
      const size_t n = m + t;
      myb_nbody* s = NULL;
      MT_CHECK(myb_nbody_create(m, t, c.gm, NULL, &s) == MYB_OK);
      double* ref = malloc(3 * n * sizeof(double));
      double* got = malloc(3 * n * sizeof(double));
      MT_CHECK(myb_internal_nbody_acceleration_path(s, c.x, ref, MYB_NBODY_PATH_RECEIVER_SCALAR) ==
               MYB_OK);
      MT_CHECK(myb_internal_nbody_acceleration_path(s, c.x, got, MYB_NBODY_PATH_SYMMETRIC) ==
               MYB_OK);
      identical = identical && same_arrays(ref, got, 3 * n);
      for (int p = 0; p < 3; ++p) {
        myb_nbody_set_pool(s, pools[p]);
        MT_CHECK(myb_internal_nbody_acceleration_path(s, c.x, got, MYB_NBODY_PATH_RECEIVER) ==
                 MYB_OK);
        identical = identical && same_arrays(ref, got, 3 * n);
      }
      myb_nbody_set_pool(s, NULL);
      free(ref);
      free(got);
      myb_nbody_destroy(s);
      free_cluster(&c);
    }
  }
  MT_CHECK(identical);
  /* Throughput of the serial paths (printed only). */
  cluster c = make_cluster(2000, 0, 77);
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(c.massive, 0, c.gm, NULL, &s) == MYB_OK);
  double* a = malloc(3 * c.massive * sizeof(double));
  const myb_internal_nbody_path paths[3] = {MYB_NBODY_PATH_SYMMETRIC, MYB_NBODY_PATH_RECEIVER,
                                            MYB_NBODY_PATH_RECEIVER_SCALAR};
  const char* names[3] = {"symmetric", "receiver, vector", "receiver, scalar"};
  const double pairs = (double)c.massive * (double)(c.massive - 1) / 2.0;
  printf("2000 bodies, one thread, time per interacting pair:");
  for (int k = 0; k < 3; ++k) {
    const clock_t t0 = clock();
    for (int r = 0; r < 3; ++r) {
      MT_CHECK(myb_internal_nbody_acceleration_path(s, c.x, a, paths[k]) == MYB_OK);
    }
    printf(" %s %.1f ns;", names[k], 1e9 * (double)(clock() - t0) / CLOCKS_PER_SEC / 3.0 / pairs);
  }
  printf("\n");
  free(a);
  myb_nbody_destroy(s);
  free_cluster(&c);
  myb_pool_destroy(pools[1]);
  myb_pool_destroy(pools[2]);
}

/* The library's arithmetic, written out: Sum2 of the terms in increasing j. */
static void manual_receiver(const cluster* c, size_t i, const unsigned char* excluded,
                            double out[3]) {
  double sum[3] = {0.0, 0.0, 0.0};
  double comp[3] = {0.0, 0.0, 0.0};
  for (size_t j = 0; j < c->massive; ++j) {
    if (j == i || excluded[j]) {
      continue;
    }
    const double dx = c->x[3 * j] - c->x[3 * i];
    const double dy = c->x[3 * j + 1] - c->x[3 * i + 1];
    const double dz = c->x[3 * j + 2] - c->x[3 * i + 2];
    const double r2 = dx * dx + dy * dy + dz * dz;
    const double sj = c->gm[j] * (1.0 / (r2 * sqrt(r2)));
    const double t[3] = {sj * dx, sj * dy, sj * dz};
    for (int k = 0; k < 3; ++k) {
      double e = 0.0;
      sum[k] = myb_eft_two_sum(sum[k], t[k], &e);
      comp[k] += e;
    }
  }
  for (int k = 0; k < 3; ++k) {
    out[k] = sum[k] + comp[k];
  }
}

static void test_graph(void) {
  cluster c = make_cluster(60, 10, 13);
  const size_t n = c.massive + c.test;
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(c.massive, c.test, c.gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_internal_nbody_reciprocal(s));
  /* Errors. */
  MT_CHECK(myb_nbody_set_interaction(s, 3, 3, MYB_INTERACTION_NONE) == MYB_E_DOMAIN);
  MT_CHECK(myb_nbody_set_interaction(s, 3, 60, MYB_INTERACTION_NONE) == MYB_E_DOMAIN);
  MT_CHECK(myb_nbody_set_interaction(s, 70, 3, MYB_INTERACTION_NONE) == MYB_E_DOMAIN);
  MT_CHECK(myb_nbody_set_interaction(s, 2, 3, (myb_interaction)7) == MYB_E_DOMAIN);
  /* A random directed graph: about one pair in five excluded, test bodies too. */
  unsigned char* excluded = calloc(n * c.massive, 1);
  mt_rng rng = {99};
  for (size_t r = 0; r < n; ++r) {
    for (size_t j = 0; j < c.massive; ++j) {
      if (j != r && mt_next(&rng) % 5 == 0) {
        excluded[r * c.massive + j] = 1;
        MT_CHECK(myb_nbody_set_interaction(s, r, j, MYB_INTERACTION_NONE) == MYB_OK);
      }
    }
  }
  /* Re-including one edge, excluding one twice. */
  MT_CHECK(myb_nbody_set_interaction(s, 5, 7, MYB_INTERACTION_NONE) == MYB_OK);
  MT_CHECK(myb_nbody_set_interaction(s, 5, 7, MYB_INTERACTION_NONE) == MYB_OK);
  MT_CHECK(myb_nbody_set_interaction(s, 5, 7, MYB_INTERACTION_NEWTON) == MYB_OK);
  excluded[5 * c.massive + 7] = 0;
  MT_CHECK(!myb_internal_nbody_reciprocal(s));
  double* a = malloc(3 * n * sizeof(double));
  double* b = malloc(3 * n * sizeof(double));
  MT_CHECK(myb_nbody_acceleration(s, 0.0, c.x, a) == MYB_OK);
  int exact = 1;
  for (size_t i = 0; i < n; ++i) {
    double ref[3];
    manual_receiver(&c, i, excluded + i * c.massive, ref);
    for (int k = 0; k < 3; ++k) {
      if (!mt_same(a[3 * i + (size_t)k], ref[k]) && exact) {
        printf("graph: receiver %zu component %d: %a, by hand %a\n", i, k, a[3 * i + (size_t)k],
               ref[k]);
      }
      exact = exact && mt_same(a[3 * i + (size_t)k], ref[k]);
    }
  }
  MT_CHECK(exact);
  /* The symmetric path cannot skip single directions: it falls back to the receiver path. */
  MT_CHECK(myb_internal_nbody_acceleration_path(s, c.x, b, MYB_NBODY_PATH_SYMMETRIC) == MYB_OK);
  MT_CHECK(same_arrays(a, b, 3 * n));
  myb_pool* pool = NULL;
  MT_CHECK(myb_pool_create(5, NULL, &pool) == MYB_OK);
  myb_nbody_set_pool(s, pool);
  MT_CHECK(myb_internal_nbody_acceleration_path(s, c.x, b, MYB_NBODY_PATH_RECEIVER) == MYB_OK);
  MT_CHECK(same_arrays(a, b, 3 * n));
  myb_nbody_set_pool(s, NULL);
  myb_pool_destroy(pool);
  free(a);
  free(b);
  free(excluded);
  myb_nbody_destroy(s);
  free_cluster(&c);
  /* Reciprocity: mirrored massive exclusions keep it, test bodies do not matter. */
  const double gm[3] = {1.0, 2.0, 3.0};
  MT_CHECK(myb_nbody_create(3, 1, gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_nbody_set_interaction(s, 1, 2, MYB_INTERACTION_NONE) == MYB_OK);
  MT_CHECK(!myb_internal_nbody_reciprocal(s));
  MT_CHECK(myb_nbody_set_interaction(s, 2, 1, MYB_INTERACTION_NONE) == MYB_OK);
  MT_CHECK(myb_nbody_set_interaction(s, 3, 0, MYB_INTERACTION_NONE) == MYB_OK);
  MT_CHECK(myb_internal_nbody_reciprocal(s));
  /* Coincident bodies: an error unless they do not interact. */
  const double x[12] = {0, 0, 0, 1, 0, 0, 1, 0, 0, 5, 5, 5};
  double acc[12];
  MT_CHECK(myb_nbody_acceleration(s, 0.0, x, acc) == MYB_OK);
  MT_CHECK(myb_nbody_set_interaction(s, 1, 2, MYB_INTERACTION_NEWTON) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, x, acc) == MYB_E_DOMAIN);
  myb_nbody_destroy(s);
}

static void test_errors_and_invariants(void) {
  myb_nbody* s = NULL;
  const double bad[2] = {1.0, -1.0};
  MT_CHECK(myb_nbody_create(0, 0, bad, NULL, &s) == MYB_E_DOMAIN);
  MT_CHECK(myb_nbody_create(2, 0, bad, NULL, &s) == MYB_E_DOMAIN);
  /* Two bodies on a circular orbit: E, L and P in closed form. */
  const double gm[2] = {3.0, 1.0};
  MT_CHECK(myb_nbody_create(2, 0, gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_nbody_count(s) == 2 && myb_nbody_massive(s) == 2 && myb_nbody_gm(s, 1) == 1.0);
  const double d = 2.0;
  const double omega = sqrt(4.0 / (d * d * d));
  const double x[6] = {-0.5, 0.0, 0.0, 1.5, 0.0, 0.0}; /* centre of mass at the origin */
  const double v[6] = {0.0, -0.5 * omega, 0.0, 0.0, 1.5 * omega, 0.0};
  myb_nbody_invariants inv = {};
  MT_CHECK(myb_nbody_invariants_of(s, x, v, &inv) == MYB_OK);
  const double kinetic = 0.5 * (3.0 * 0.25 + 1.0 * 2.25) * omega * omega;
  const double energy = kinetic - 3.0 / d;
  printf("circular pair: energy %.17g (closed form %.17g)\n", myb_dd_to_double(inv.energy),
         energy);
  MT_CHECK(fabs(myb_dd_to_double(inv.energy) - energy) < 1e-15);
  MT_CHECK(fabs(inv.momentum[1]) < 1e-16 && inv.momentum[0] == 0.0);
  MT_CHECK(fabs(inv.angular_momentum[2] - (3.0 * 0.25 + 2.25) * omega) < 1e-15);
  MT_CHECK(fabs(inv.center[0]) < 1e-16 && fabs(inv.center_velocity[1]) < 1e-16);
  double a[6];
  MT_CHECK(myb_nbody_acceleration(s, 0.0, x, a) == MYB_OK);
  MT_CHECK_SAME(a[0], 1.0 * (1.0 / (4.0 * 2.0)) * 2.0);
  MT_CHECK_SAME(a[3], -(3.0 * (1.0 / (4.0 * 2.0)) * 2.0));
  const double same[6] = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  MT_CHECK(myb_nbody_invariants_of(s, same, v, &inv) == MYB_E_DOMAIN);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, same, a) == MYB_E_DOMAIN);
  myb_nbody_destroy(s);
  myb_nbody_destroy(NULL);
}

/* ---- Classical solutions ------------------------------------------------------------------- */

static myb_ias15* ias15_for(myb_nbody* s, double t0, const double* x, const double* v,
                            double epsilon) {
  const size_t n = 3 * myb_nbody_count(s);
  const myb_ias15_options o = {.epsilon = epsilon, .block = 3};
  myb_ias15* in = NULL;
  MT_CHECK(myb_ias15_create(n, &o, NULL, &in) == MYB_OK);
  double* y = malloc(2 * n * sizeof(double));
  for (size_t i = 0; i < n; ++i) {
    y[i] = x[i];
    y[n + i] = v[i];
  }
  MT_CHECK(myb_ias15_start(in, myb_nbody_acceleration_xv, s, t0, y) == MYB_OK);
  free(y);
  return in;
}

static void test_lagrange(void) {
  const double gm[3] = {1.0, 1e-3, 1e-6}; /* Routh: 27 sum m_i m_j < M^2, linearly stable */
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(3, 0, gm, NULL, &s) == MYB_OK);
  const double total = gm[0] + gm[1] + gm[2];
  double x[9];
  double v[9];
  double c[2] = {0.0, 0.0};
  for (int k = 0; k < 3; ++k) {
    const double angle = pi / 2.0 + 2.0 * pi * k / 3.0;
    x[3 * k] = cos(angle) / sqrt(3.0); /* side 1 */
    x[3 * k + 1] = sin(angle) / sqrt(3.0);
    x[3 * k + 2] = 0.0;
    c[0] += gm[k] * x[3 * k] / total;
    c[1] += gm[k] * x[3 * k + 1] / total;
  }
  const double omega = sqrt(total);
  for (int k = 0; k < 3; ++k) {
    x[3 * k] -= c[0];
    x[3 * k + 1] -= c[1];
    v[3 * k] = -omega * x[3 * k + 1];
    v[3 * k + 1] = omega * x[3 * k];
    v[3 * k + 2] = 0.0;
  }
  myb_ias15* in = ias15_for(s, 0.0, x, v, 0.0);
  const double t_end = 20.0 * 2.0 * pi / omega;
  MT_CHECK(myb_ias15_integrate(in, t_end) == MYB_OK);
  const double* y = myb_ias15_state(in);
  double worst = 0.0;
  const double ca = cos(omega * t_end);
  const double sa = sin(omega * t_end);
  for (int k = 0; k < 3; ++k) {
    const double ex = ca * x[3 * k] - sa * x[3 * k + 1];
    const double ey = sa * x[3 * k] + ca * x[3 * k + 1];
    worst = fmax(worst, hypot(y[3 * k] - ex, y[3 * k + 1] - ey)); /* side 1 */
  }
  printf("Lagrange equilateral (GM 1, 1e-3, 1e-6), 20 periods by IAS15: deviation from rigid "
         "rotation %.2e of the side\n",
         worst);
  MT_CHECK(worst < 3e-13);
  myb_ias15_destroy(in);
  myb_nbody_destroy(s);
}

/* Figure eight in the parametrization of Suvakov & Dmitrasinovic: x1 = (-1, 0), x2 = (1, 0),
 * x3 = 0, v1 = v2 = (p1, p2), v3 = -2 (p1, p2); GM = 1. */
static void eight_start(double p1, double p2, double x[9], double v[9]) {
  const double xs[9] = {-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double vs[9] = {p1, p2, 0.0, p1, p2, 0.0, -2.0 * p1, -2.0 * p2, 0.0};
  for (int i = 0; i < 9; ++i) {
    x[i] = xs[i];
    v[i] = vs[i];
  }
}

/* F = y(T) - y(0) (18 values); field = (v, a) at y(T). */
static void eight_residual(myb_nbody* s, double p1, double p2, double period, double f[18],
                           double field[18]) {
  double x[9];
  double v[9];
  eight_start(p1, p2, x, v);
  myb_ias15* in = ias15_for(s, 0.0, x, v, 0.0);
  MT_CHECK(myb_ias15_integrate(in, period) == MYB_OK);
  const double* y = myb_ias15_state(in);
  for (int i = 0; i < 9; ++i) {
    f[i] = y[i] - x[i];
    f[9 + i] = y[9 + i] - v[i];
  }
  if (field != NULL) {
    double a[9];
    MT_CHECK(myb_nbody_acceleration(s, period, y, a) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      field[i] = y[9 + i];
      field[9 + i] = a[i];
    }
  }
  myb_ias15_destroy(in);
}

static double norm18(const double f[18]) {
  double sum = 0.0;
  for (int i = 0; i < 18; ++i) {
    sum += f[i] * f[i];
  }
  return sqrt(sum);
}

/* Solves the 3 x 3 system m z = r by Gaussian elimination with pivoting. */
static void solve3(double m[3][3], double r[3], double z[3]) {
  for (int p = 0; p < 3; ++p) {
    int best = p;
    for (int i = p + 1; i < 3; ++i) {
      if (fabs(m[i][p]) > fabs(m[best][p])) {
        best = i;
      }
    }
    for (int j = 0; j < 3; ++j) {
      const double t = m[p][j];
      m[p][j] = m[best][j];
      m[best][j] = t;
    }
    const double t = r[p];
    r[p] = r[best];
    r[best] = t;
    for (int i = p + 1; i < 3; ++i) {
      const double f = m[i][p] / m[p][p];
      for (int j = p; j < 3; ++j) {
        m[i][j] -= f * m[p][j];
      }
      r[i] -= f * r[p];
    }
  }
  for (int i = 2; i >= 0; --i) {
    double sum = r[i];
    for (int j = i + 1; j < 3; ++j) {
      sum -= m[i][j] * z[j];
    }
    z[i] = sum / m[i][i];
  }
}

static void test_figure_eight(void) {
  const double gm[3] = {1.0, 1.0, 1.0};
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(3, 0, gm, NULL, &s) == MYB_OK);
  double u[3] = {0.347111, 0.532728, 6.324449}; /* Suvakov & Dmitrasinovic 2013, 6 digits */
  double f[18];
  double field[18];
  eight_residual(s, u[0], u[1], u[2], f, NULL);
  printf("figure eight: residual of the published 6-digit start %.2e\n", norm18(f));
  double residual = 0.0;
  /* bound: Gauss-Newton converges in 3-4 iterations; at most 8 */
  for (int iteration = 0; iteration < 8; ++iteration) {
    eight_residual(s, u[0], u[1], u[2], f, field);
    residual = norm18(f);
    double jac[18][3];
    const double delta = 1e-6;
    for (int c = 0; c < 2; ++c) {
      double up[3] = {u[0], u[1], u[2]};
      double down[3] = {u[0], u[1], u[2]};
      up[c] += delta;
      down[c] -= delta;
      double fu[18];
      double fd[18];
      eight_residual(s, up[0], up[1], up[2], fu, NULL);
      eight_residual(s, down[0], down[1], down[2], fd, NULL);
      for (int i = 0; i < 18; ++i) {
        jac[i][c] = (fu[i] - fd[i]) / (2.0 * delta);
      }
    }
    for (int i = 0; i < 18; ++i) {
      jac[i][2] = field[i]; /* dF/dT */
    }
    double normal[3][3] = {};
    double rhs[3] = {};
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        for (int i = 0; i < 18; ++i) {
          normal[a][b] += jac[i][a] * jac[i][b];
        }
      }
      for (int i = 0; i < 18; ++i) {
        rhs[a] -= jac[i][a] * f[i];
      }
    }
    double step[3];
    solve3(normal, rhs, step);
    for (int k = 0; k < 3; ++k) {
      u[k] += step[k];
    }
    if (fabs(step[0]) + fabs(step[1]) + fabs(step[2]) < 1e-14) {
      break;
    }
  }
  eight_residual(s, u[0], u[1], u[2], f, NULL);
  residual = norm18(f);
  printf("figure eight refined: p1 = %.16f, p2 = %.16f, T = %.16f, residual after one period "
         "%.2e\n",
         u[0], u[1], u[2], residual);
  MT_CHECK(residual < 1e-12);
  /* Simo's initial condition (x1 = -x2 = (0.97000436, -0.24308753), v3 = -2 v1 =
   * (-0.93240737, -0.86473146), T = 6.32591398), rotated and scaled to this parametrization:
   * positions / L, velocities * sqrt(L), periods / L^1.5, L = |x1|. */
  const double sx = 0.97000436;
  const double sy = -0.24308753;
  const double length = hypot(sx, sy);
  const double angle = atan2(-sy, -sx); /* direction of Simo's x2 */
  const double svx = 0.93240737 / 2.0;
  const double svy = 0.86473146 / 2.0;
  const double p1 = (cos(angle) * svx + sin(angle) * svy) * sqrt(length);
  const double p2 = (-sin(angle) * svx + cos(angle) * svy) * sqrt(length);
  const double period = 6.32591398 / pow(length, 1.5);
  /* Simo's orientation is the time-reversed orbit here: v -> -v. */
  printf("Simo's start in this parametrization: p1 = %.9f, p2 = %.9f, T = %.9f\n", -p1, -p2,
         period);
  MT_CHECK(fabs(-p1 - u[0]) < 1e-8 && fabs(-p2 - u[1]) < 1e-8 && fabs(period - u[2]) < 5e-8);
  /* Ten periods of the refined orbit: a choreography, stable. */
  double x[9];
  double v[9];
  eight_start(u[0], u[1], x, v);
  myb_ias15* in = ias15_for(s, 0.0, x, v, 0.0);
  double worst = 0.0;
  myb_nbody_invariants inv0 = {};
  MT_CHECK(myb_nbody_invariants_of(s, x, v, &inv0) == MYB_OK);
  for (int k = 1; k <= 10; ++k) {
    MT_CHECK(myb_ias15_integrate(in, k * u[2]) == MYB_OK);
    const double* y = myb_ias15_state(in);
    for (int i = 0; i < 9; ++i) {
      worst = fmax(worst, fabs(y[i] - x[i]));
    }
  }
  const double* y = myb_ias15_state(in);
  myb_nbody_invariants inv = {};
  MT_CHECK(myb_nbody_invariants_of(s, y, y + 9, &inv) == MYB_OK);
  const double de = myb_dd_to_double(myb_dd_sub(inv.energy, inv0.energy)) /
                    fabs(myb_dd_to_double(inv0.energy));
  printf("figure eight, 10 periods: return to the start within %.2e, relative energy error "
         "%.2e, angular momentum %.1e\n",
         worst, de, inv.angular_momentum[2]);
  MT_CHECK(worst < 1e-10);
  MT_CHECK(fabs(de) < 1e-14);
  myb_ias15_destroy(in);
  myb_nbody_destroy(s);
}

/* Pairwise outcome at the end: which pair is bound, and whether the third body escapes. */
static void pythagorean_outcome(const double gm[3], const double* y, int* bound_a, int* bound_b,
                                int* escapes) {
  *bound_a = -1;
  *bound_b = -1;
  for (int i = 0; i < 3; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      double dx[3];
      double dv[3];
      for (int c = 0; c < 3; ++c) {
        dx[c] = y[3 * j + c] - y[3 * i + c];
        dv[c] = y[9 + 3 * j + c] - y[9 + 3 * i + c];
      }
      const double r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
      const double e = 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]) - (gm[i] + gm[j]) / r;
      if (e < 0.0) {
        *bound_a = i;
        *bound_b = j;
      }
    }
  }
  *escapes = 0;
  if (*bound_a >= 0) {
    const int third = 3 - *bound_a - *bound_b;
    const double mb = gm[*bound_a] + gm[*bound_b];
    double rc[3];
    double vc[3];
    for (int c = 0; c < 3; ++c) {
      rc[c] = y[3 * third + c] -
              (gm[*bound_a] * y[3 * *bound_a + c] + gm[*bound_b] * y[3 * *bound_b + c]) / mb;
      vc[c] = y[9 + 3 * third + c] - (gm[*bound_a] * y[9 + 3 * *bound_a + c] +
                                      gm[*bound_b] * y[9 + 3 * *bound_b + c]) /
                                         mb;
    }
    const double r = sqrt(rc[0] * rc[0] + rc[1] * rc[1] + rc[2] * rc[2]);
    const double e = 0.5 * (vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]) - (mb + gm[third]) / r;
    const double radial = rc[0] * vc[0] + rc[1] * vc[1] + rc[2] * vc[2];
    *escapes = e > 0.0 && radial > 0.0;
  }
}

static void test_pythagorean(void) {
  const double gm[3] = {3.0, 4.0, 5.0};
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(3, 0, gm, NULL, &s) == MYB_OK);
  const double x[9] = {1.0, 3.0, 0.0, -2.0, -1.0, 0.0, 1.0, -1.0, 0.0};
  const double v[9] = {};
  myb_nbody_invariants inv0 = {};
  MT_CHECK(myb_nbody_invariants_of(s, x, v, &inv0) == MYB_OK);
  MT_CHECK(fabs(myb_dd_to_double(inv0.energy) + 769.0 / 60.0) < 1e-14);
  const double t_early = 10.0; /* before the closest encounters (t = 15.8) */
  const double t_end = 70.0;
  double early[2][18];
  double finals[2][18];
  double energy_error[2];
  double closest = 1e300;
  double closest_t = 0.0;
  /* IAS15 at its default precision, and RKN12(10) at a tight tolerance: two independent
   * methods through the close encounters of a chaotic problem. */
  for (int method = 0; method < 2; ++method) {
    const double* y = NULL;
    myb_ias15* in = NULL;
    myb_rkn12* rk = NULL;
    if (method == 0) {
      in = ias15_for(s, 0.0, x, v, 0.0);
      for (int phase = 0; phase < 2; ++phase) {
        const double t_stop = phase == 0 ? t_early : t_end;
        /* bound: at most 10^6 steps */
        for (int i = 0; i < 1000000 && myb_ias15_time(in) < t_stop; ++i) {
          MT_CHECK(myb_ias15_step(in, t_stop) == MYB_OK);
          y = myb_ias15_state(in);
          for (int a = 0; a < 3; ++a) {
            for (int b = a + 1; b < 3; ++b) {
              const double r = hypot(y[3 * b] - y[3 * a], y[3 * b + 1] - y[3 * a + 1]);
              if (r < closest) {
                closest = r;
                closest_t = myb_ias15_time(in);
              }
            }
          }
        }
        if (phase == 0) {
          memcpy(early[method], myb_ias15_state(in), sizeof early[method]);
        }
      }
      y = myb_ias15_state(in);
      const myb_ias15_stats st = myb_ias15_statistics(in);
      printf("Pythagorean problem by IAS15: %lld steps; closest approach %.3e at t = %.3f\n",
             (long long)st.accepted, closest, closest_t);
    } else {
      const myb_rkn12_options o = {.rtol = 1e-14, .atol = 1e-16};
      MT_CHECK(myb_rkn12_create(9, &o, NULL, &rk) == MYB_OK);
      double y0[18];
      for (int i = 0; i < 9; ++i) {
        y0[i] = x[i];
        y0[9 + i] = v[i];
      }
      MT_CHECK(myb_rkn12_start(rk, myb_nbody_acceleration, s, 0.0, y0) == MYB_OK);
      MT_CHECK(myb_rkn12_integrate(rk, t_early) == MYB_OK);
      memcpy(early[method], myb_rkn12_state(rk), sizeof early[method]);
      MT_CHECK(myb_rkn12_integrate(rk, t_end) == MYB_OK);
      y = myb_rkn12_state(rk);
      printf("Pythagorean problem by RKN12(10): %lld steps\n",
             (long long)myb_rkn12_statistics(rk).accepted);
    }
    for (int i = 0; i < 18; ++i) {
      finals[method][i] = y[i];
    }
    myb_nbody_invariants inv = {};
    MT_CHECK(myb_nbody_invariants_of(s, y, y + 9, &inv) == MYB_OK);
    energy_error[method] = myb_dd_to_double(myb_dd_sub(inv.energy, inv0.energy)) /
                           fabs(myb_dd_to_double(inv0.energy));
    myb_ias15_destroy(in);
    myb_rkn12_destroy(rk);
  }
  double difference = 0.0;
  double difference_early = 0.0;
  for (int i = 0; i < 9; ++i) {
    difference = fmax(difference, fabs(finals[0][i] - finals[1][i]));
    difference_early = fmax(difference_early, fabs(early[0][i] - early[1][i]));
  }
  int a = 0;
  int b = 0;
  int escapes = 0;
  pythagorean_outcome(gm, finals[0], &a, &b, &escapes);
  printf("positions of IAS15 and RKN12 differ by %.2e at t = 10 and by %.2e at t = 70 "
         "(chaos amplifies the round-off of the encounters)\n",
         difference_early, difference);
  printf("at t = 70: bound pair (GM %.0f, %.0f), GM %.0f escapes: %s, velocity (%.4f, %.4f); "
         "relative energy errors %.1e, %.1e\n",
         a >= 0 ? gm[a] : 0.0, b >= 0 ? gm[b] : 0.0, a >= 0 ? gm[3 - a - b] : 0.0,
         escapes ? "yes" : "no", finals[0][9], finals[0][10], energy_error[0], energy_error[1]);
  MT_CHECK(difference_early < 1e-11);
  MT_CHECK(a == 1 && b == 2 && escapes);
  int a2 = 0;
  int b2 = 0;
  int escapes2 = 0;
  pythagorean_outcome(gm, finals[1], &a2, &b2, &escapes2);
  MT_CHECK(a2 == 1 && b2 == 2 && escapes2);
  MT_CHECK(fabs(energy_error[0]) < 1e-10 && fabs(energy_error[1]) < 1e-10);
  MT_CHECK(difference < 0.2); /* same branch of the outcome, not the same digits */
  myb_nbody_destroy(s);
}

/* A whole integration with and without threads: the same bits. */
static void test_integration_threads(void) {
  cluster c = make_cluster(300, 0, 14);
  const size_t n = c.massive;
  for (size_t i = 0; i < n; ++i) {
    c.gm[i] = 1e-3 * (1.0 + uniform(&(mt_rng){i}));
    for (int k = 0; k < 3; ++k) {
      c.x[3 * i + (size_t)k] *= 1e-11; /* unit cube */
    }
  }
  double* v = calloc(3 * n, sizeof(double));
  double* finals[2] = {malloc(6 * n * sizeof(double)), malloc(6 * n * sizeof(double))};
  for (int run = 0; run < 2; ++run) {
    myb_nbody* s = NULL;
    MT_CHECK(myb_nbody_create(n, 0, c.gm, NULL, &s) == MYB_OK);
    myb_pool* pool = NULL;
    if (run == 1) {
      MT_CHECK(myb_pool_create(6, NULL, &pool) == MYB_OK);
      myb_nbody_set_pool(s, pool);
    }
    myb_ias15* in = ias15_for(s, 0.0, c.x, v, 0.0);
    for (int i = 0; i < 40; ++i) {
      MT_CHECK(myb_ias15_step(in, 1e9) == MYB_OK);
    }
    memcpy(finals[run], myb_ias15_state(in), 6 * n * sizeof(double));
    myb_ias15_destroy(in);
    myb_nbody_destroy(s);
    myb_pool_destroy(pool);
  }
  MT_CHECK(same_arrays(finals[0], finals[1], 6 * n));
  free(finals[0]);
  free(finals[1]);
  free(v);
  free_cluster(&c);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0); /* keep the output of a fatal check */
  test_errors_and_invariants();
  test_accuracy();
  test_paths_and_threads();
  test_vector_paths();
  test_graph();
  test_integration_threads();
  test_lagrange();
  test_figure_eight();
  test_pythagorean();
  return mt_finish("test_nbody");
}
