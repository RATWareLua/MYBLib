/* SPDX-License-Identifier: BSD-3-Clause */
/* Newtonian N-body accelerations by direct summation over all pairs.
 *
 * Every term of receiver i is GM_j q (x_j - x_i) with q = 1 / (r2 sqrt(r2)) and
 * r2 = ((dx dx + dy dy) + dz dz), d = x_j - x_i, and the terms enter a Sum2 accumulator
 * (TwoSum per term, errors summed apart) in increasing j. The symmetric path computes the pair
 * (i, j), i < j, once: the term of j is -(GM_i q dx), which equals GM_i q (x_i - x_j) bit for
 * bit because x_i - x_j = -(x_j - x_i) exactly and r2 is the same. Receiver j sees the
 * sources i < j at the outer steps 0 .. j-1 and the sources k > j at its own outer step, i.e.
 * in increasing order, as on the receiver path. Both paths therefore give the same bits, and
 * the receiver path gives the same bits whatever thread sums which receiver. The receiver path
 * sums blocks of four receivers in vector lanes (sum_receivers4, the same operations lane by
 * lane), the remainder one by one (sum_receiver, the scalar reference).
 *
 * Exclusions are kept sorted by (receiver, source); a receiver sums the gaps between its
 * excluded sources (and itself), so the inner loops have no branch. */
#include <math.h> /* sqrt */
#include <stdckdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "core/pool.h"
#include "dyn/nbody_internal.h"
#include "myb/dyn.h"
#include "myb/num.h"
#include "num/eft.h"

enum {
  PARALLEL_MIN_PAIRS = 1 << 14, /* below this many pairs the pool does not pay */
  PART_PAIRS = 1 << 14,         /* pairs per part of the parallel loop */
  VECTOR_MIN_BODIES = 16,       /* from here the vector receiver path is the faster one */
};

typedef struct exclusion {
  size_t receiver;
  size_t source;
} exclusion;

typedef struct accumulator {
  double s[3];
  double c[3];
} accumulator;

struct myb_nbody {
  myb_allocator allocator;
  size_t massive;
  size_t count;
  double* gm; /* massive */
  double* px; /* count each: positions of the current evaluation, by component */
  double* py;
  double* pz;
  accumulator* acc; /* count: the symmetric path */
  myb_pool* pool;
  exclusion* excluded; /* sorted by (receiver, source) */
  size_t excluded_count;
  size_t excluded_capacity;
  size_t* first; /* count + 1: excluded[first[i] .. first[i + 1]) have receiver i */
  int graph_valid;
  int reciprocal;
};

static int finite(double x) {
  return x - x == 0.0;
}

myb_status myb_nbody_create(size_t massive, size_t test, const double* gm,
                            const myb_allocator* allocator, myb_nbody** nbody) {
  MYB_CHECK(nbody != NULL && gm != NULL);
  if (massive == 0) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < massive; ++i) {
    if (!finite(gm[i]) || gm[i] < 0.0) {
      return MYB_E_DOMAIN;
    }
  }
  size_t count = 0;
  if (ckd_add(&count, massive, test)) {
    return MYB_E_OVERFLOW;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_nbody), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_nbody* s = block;
  *s = (myb_nbody){.allocator = a, .massive = massive, .count = count};
  size_t doubles = 0;
  size_t firsts = 0;
  if (ckd_mul(&doubles, count, (size_t)3) || ckd_add(&doubles, doubles, massive) ||
      ckd_add(&firsts, count, (size_t)1)) {
    status = MYB_E_OVERFLOW;
    goto fail;
  }
  status = myb_internal_alloc_array(&a, doubles, sizeof(double), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  s->gm = block;
  s->px = s->gm + massive;
  s->py = s->px + count;
  s->pz = s->py + count;
  memcpy(s->gm, gm, massive * sizeof(double));
  status = myb_internal_alloc_array(&a, count, sizeof(accumulator), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  s->acc = block;
  status = myb_internal_alloc_array(&a, firsts, sizeof(size_t), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  s->first = block;
  *nbody = s;
  return MYB_OK;

fail:
  myb_internal_release(&a, s->first);
  myb_internal_release(&a, s->acc);
  myb_internal_release(&a, s->gm);
  myb_internal_release(&a, s);
  return status;
}

void myb_nbody_destroy(myb_nbody* nbody) {
  if (nbody == NULL) {
    return;
  }
  const myb_allocator a = nbody->allocator;
  myb_internal_release(&a, nbody->excluded);
  myb_internal_release(&a, nbody->first);
  myb_internal_release(&a, nbody->acc);
  myb_internal_release(&a, nbody->gm);
  myb_internal_release(&a, nbody);
}

size_t myb_nbody_count(const myb_nbody* nbody) {
  MYB_CHECK(nbody != NULL);
  return nbody->count;
}

size_t myb_nbody_massive(const myb_nbody* nbody) {
  MYB_CHECK(nbody != NULL);
  return nbody->massive;
}

double myb_nbody_gm(const myb_nbody* nbody, size_t body) {
  MYB_CHECK(nbody != NULL && body < nbody->count);
  return body < nbody->massive ? nbody->gm[body] : 0.0;
}

void myb_nbody_set_pool(myb_nbody* nbody, myb_pool* pool) {
  MYB_CHECK(nbody != NULL);
  nbody->pool = pool;
}

/* ---- Interaction graph ------------------------------------------------------------------- */

static int exclusion_less(const exclusion* a, size_t receiver, size_t source) {
  return a->receiver < receiver || (a->receiver == receiver && a->source < source);
}

/* First position whose key is not less than (receiver, source). */
static size_t exclusion_lower_bound(const myb_nbody* s, size_t receiver, size_t source) {
  size_t lo = 0;
  size_t hi = s->excluded_count;
  /* bound: the interval halves, at most 64 iterations */
  while (lo < hi) {
    const size_t mid = lo + (hi - lo) / 2;
    if (exclusion_less(&s->excluded[mid], receiver, source)) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

static int exclusion_found(const myb_nbody* s, size_t at, size_t receiver, size_t source) {
  return at < s->excluded_count && s->excluded[at].receiver == receiver &&
         s->excluded[at].source == source;
}

myb_status myb_nbody_set_interaction(myb_nbody* nbody, size_t receiver, size_t source,
                                     myb_interaction level) {
  MYB_CHECK(nbody != NULL);
  myb_nbody* s = nbody;
  if (receiver >= s->count || source >= s->massive || receiver == source ||
      (level != MYB_INTERACTION_NONE && level != MYB_INTERACTION_NEWTON)) {
    return MYB_E_DOMAIN;
  }
  const size_t at = exclusion_lower_bound(s, receiver, source);
  const int present = exclusion_found(s, at, receiver, source);
  if (level == MYB_INTERACTION_NEWTON) {
    if (present) {
      memmove(&s->excluded[at], &s->excluded[at + 1],
              (s->excluded_count - at - 1) * sizeof(exclusion));
      s->excluded_count--;
      s->graph_valid = 0;
    }
    return MYB_OK;
  }
  if (present) {
    return MYB_OK;
  }
  if (s->excluded_count == s->excluded_capacity) {
    size_t capacity = s->excluded_capacity == 0 ? 16 : s->excluded_capacity;
    if (s->excluded_capacity != 0 && ckd_mul(&capacity, capacity, (size_t)2)) {
      return MYB_E_OVERFLOW;
    }
    void* block = NULL;
    const myb_status status =
        myb_internal_alloc_array(&s->allocator, capacity, sizeof(exclusion), &block);
    if (status != MYB_OK) {
      return status;
    }
    if (s->excluded_count != 0) {
      memcpy(block, s->excluded, s->excluded_count * sizeof(exclusion));
    }
    myb_internal_release(&s->allocator, s->excluded);
    s->excluded = block;
    s->excluded_capacity = capacity;
  }
  memmove(&s->excluded[at + 1], &s->excluded[at], (s->excluded_count - at) * sizeof(exclusion));
  s->excluded[at] = (exclusion){.receiver = receiver, .source = source};
  s->excluded_count++;
  s->graph_valid = 0;
  return MYB_OK;
}

/* Index of the exclusions by receiver, and whether the massive part is reciprocal. */
static void graph_update(myb_nbody* s) {
  if (s->graph_valid) {
    return;
  }
  size_t e = 0;
  for (size_t i = 0; i <= s->count; ++i) {
    /* bound: e only grows, to excluded_count */
    while (e < s->excluded_count && s->excluded[e].receiver < i) {
      ++e;
    }
    s->first[i] = e;
  }
  s->reciprocal = 1;
  for (size_t k = 0; k < s->excluded_count; ++k) {
    const exclusion* x = &s->excluded[k];
    if (x->receiver < s->massive &&
        !exclusion_found(s, exclusion_lower_bound(s, x->source, x->receiver), x->source,
                         x->receiver)) {
      s->reciprocal = 0;
    }
  }
  s->graph_valid = 1;
}

int myb_internal_nbody_reciprocal(myb_nbody* nbody) {
  MYB_CHECK(nbody != NULL);
  graph_update(nbody);
  return nbody->reciprocal;
}

/* ---- Accelerations ---------------------------------------------------------------------- */

static inline void accumulate(accumulator* a, double x, double y, double z) {
  double e = 0.0;
  a->s[0] = myb_eft_two_sum(a->s[0], x, &e);
  a->c[0] += e;
  a->s[1] = myb_eft_two_sum(a->s[1], y, &e);
  a->c[1] += e;
  a->s[2] = myb_eft_two_sum(a->s[2], z, &e);
  a->c[2] += e;
}

static inline double inverse_cube(double dx, double dy, double dz) {
  const double r2 = dx * dx + dy * dy + dz * dz;
  return 1.0 / (r2 * sqrt(r2));
}

/* Sources [begin, end) acting on a receiver at (xi, yi, zi). */
static void sum_sources(const myb_nbody* s, double xi, double yi, double zi, size_t begin,
                        size_t end, accumulator* acc) {
  const double* gm = s->gm;
  const double* px = s->px;
  const double* py = s->py;
  const double* pz = s->pz;
  for (size_t j = begin; j < end; ++j) {
    const double dx = px[j] - xi;
    const double dy = py[j] - yi;
    const double dz = pz[j] - zi;
    const double sj = gm[j] * inverse_cube(dx, dy, dz);
    accumulate(acc, sj * dx, sj * dy, sj * dz);
  }
}

/* Receiver path: all sources of receiver i in increasing order, skipping i itself and its
 * excluded sources. */
static void sum_receiver(const myb_nbody* s, size_t i, double* a) {
  const double xi = s->px[i];
  const double yi = s->py[i];
  const double zi = s->pz[i];
  accumulator acc = {};
  size_t e = s->excluded_count == 0 ? 0 : s->first[i];
  const size_t e_end = s->excluded_count == 0 ? 0 : s->first[i + 1];
  size_t j = 0;
  const size_t passes = e_end - e + 2;
  /* bound: every pass moves j past one skipped index; at most (exclusions of i) + 2 passes */
  for (size_t pass = 0; pass < passes && j < s->massive; ++pass) {
    size_t skip = s->massive;
    if (i >= j && i < skip) {
      skip = i;
    }
    if (e < e_end && s->excluded[e].source < skip) {
      skip = s->excluded[e].source;
    }
    sum_sources(s, xi, yi, zi, j, skip, &acc);
    if (e < e_end && s->excluded[e].source == skip) {
      ++e;
    }
    j = skip + 1;
  }
  MYB_CHECK(j >= s->massive && e == e_end); /* every source seen, every exclusion consumed */
  for (int k = 0; k < 3; ++k) {
    a[3 * i + (size_t)k] = acc.s[k] + acc.c[k];
  }
}

/* ---- Four receivers at once ------------------------------------------------------------------
 * Clang's portable vectors (no intrinsics: two SSE registers on x86-64-v2, one AVX register on
 * v3, NEON pairs on AArch64). Every lane performs exactly the scalar operations of
 * sum_sources in the same order (sub, mul, add, sqrt, mul, div, mul, mul, TwoSum), and IEEE-754
 * defines each of them lane by lane, so a lane gives the bits of sum_receiver. A lane's own
 * index is masked: its term becomes +0 (bitwise and), and adding +0 to a Sum2 accumulator
 * changes no bit because neither s nor c can ever hold -0 (s starts at +0 and x + y is -0 in
 * round-to-nearest only for -0 + -0; the error a - RN(..) of TwoSum is +0 when zero). */
typedef double v4d __attribute__((vector_size(32)));
typedef long long v4i __attribute__((vector_size(32)));

enum { LANES = 4 };

typedef struct accumulator4 {
  v4d s[3];
  v4d c[3];
} accumulator4;

static inline v4d splat(double x) {
  return (v4d){x, x, x, x};
}

static inline void accumulate4(accumulator4* a, int k, v4d x) {
  const v4d s = a->s[k] + x;
  const v4d a_rounded = s - x;
  const v4d b_rounded = s - a_rounded;
  const v4d a_error = a->s[k] - a_rounded;
  const v4d b_error = x - b_rounded;
  a->s[k] = s;
  a->c[k] += a_error + b_error;
}

/* Sources [begin, end) on the four receivers at (xi, yi, zi); lanes whose receiver index
 * equals the source are masked when mask is set (receivers i0 .. i0 + 3). */
static void sum_sources4(const myb_nbody* s, v4d xi, v4d yi, v4d zi, size_t begin, size_t end,
                         int mask, size_t i0, accumulator4* acc) {
  const double* gm = s->gm;
  const double* px = s->px;
  const double* py = s->py;
  const double* pz = s->pz;
  const v4i lane_index = {(long long)i0, (long long)i0 + 1, (long long)i0 + 2,
                          (long long)i0 + 3};
  for (size_t j = begin; j < end; ++j) {
    const v4d dx = splat(px[j]) - xi;
    const v4d dy = splat(py[j]) - yi;
    const v4d dz = splat(pz[j]) - zi;
    const v4d r2 = dx * dx + dy * dy + dz * dz;
    const v4d q = splat(1.0) / (r2 * __builtin_elementwise_sqrt(r2));
    const v4d sj = splat(gm[j]) * q;
    v4d t[3] = {sj * dx, sj * dy, sj * dz};
    if (mask) {
      const v4i keep = lane_index != (v4i){(long long)j, (long long)j, (long long)j, (long long)j};
      for (int k = 0; k < 3; ++k) {
        t[k] = (v4d)((v4i)t[k] & keep);
      }
    }
    for (int k = 0; k < 3; ++k) {
      accumulate4(acc, k, t[k]);
    }
  }
}

/* Receivers i0 .. i0 + 3, none of which has an exclusion. */
static void sum_receivers4(const myb_nbody* s, size_t i0, double* a) {
  const v4d xi = {s->px[i0], s->px[i0 + 1], s->px[i0 + 2], s->px[i0 + 3]};
  const v4d yi = {s->py[i0], s->py[i0 + 1], s->py[i0 + 2], s->py[i0 + 3]};
  const v4d zi = {s->pz[i0], s->pz[i0 + 1], s->pz[i0 + 2], s->pz[i0 + 3]};
  accumulator4 acc = {};
  const size_t m = s->massive;
  const size_t lo = i0 < m ? i0 : m;          /* sources below every lane: no mask */
  const size_t hi = i0 + LANES < m ? i0 + LANES : m; /* [lo, hi): one lane is the source */
  sum_sources4(s, xi, yi, zi, 0, lo, 0, i0, &acc);
  sum_sources4(s, xi, yi, zi, lo, hi, 1, i0, &acc);
  sum_sources4(s, xi, yi, zi, hi, m, 0, i0, &acc);
  for (int k = 0; k < 3; ++k) {
    const v4d sum = acc.s[k] + acc.c[k];
    for (size_t lane = 0; lane < LANES; ++lane) {
      a[3 * (i0 + lane) + (size_t)k] = sum[lane];
    }
  }
}

typedef struct receiver_job {
  const myb_nbody* s;
  double* a;
  int vector;
} receiver_job;

/* Receivers [begin, end): blocks of four without exclusions on the vector path, the rest one
 * by one; the bits do not depend on the blocking. */
static void receiver_range(void* context, size_t begin, size_t end) {
  const receiver_job* job = context;
  const myb_nbody* s = job->s;
  size_t i = begin;
  /* bound: i advances by 1 or LANES up to end */
  while (i < end) {
    const int free_block =
        job->vector && end - i >= LANES &&
        (s->excluded_count == 0 || s->first[i] == s->first[i + LANES]);
    if (free_block) {
      sum_receivers4(s, i, job->a);
      i += LANES;
    } else {
      sum_receiver(s, i, job->a);
      i += 1;
    }
  }
}

/* Symmetric path: each massive pair once; then the test bodies. Needs no exclusion. */
static void sum_symmetric(myb_nbody* s, double* a) {
  accumulator* acc = s->acc;
  const size_t m = s->massive;
  memset(acc, 0, m * sizeof(accumulator));
  const double* gm = s->gm;
  const double* px = s->px;
  const double* py = s->py;
  const double* pz = s->pz;
  for (size_t i = 0; i < m; ++i) {
    const double xi = px[i];
    const double yi = py[i];
    const double zi = pz[i];
    const double gi = gm[i];
    accumulator ai = acc[i];
    for (size_t j = i + 1; j < m; ++j) {
      const double dx = px[j] - xi;
      const double dy = py[j] - yi;
      const double dz = pz[j] - zi;
      const double q = inverse_cube(dx, dy, dz);
      const double si = gi * q;
      const double sj = gm[j] * q;
      accumulate(&ai, sj * dx, sj * dy, sj * dz);
      accumulate(&acc[j], -(si * dx), -(si * dy), -(si * dz));
    }
    for (int k = 0; k < 3; ++k) {
      a[3 * i + (size_t)k] = ai.s[k] + ai.c[k];
    }
  }
  for (size_t i = m; i < s->count; ++i) {
    sum_receiver(s, i, a);
  }
}

myb_status myb_internal_nbody_acceleration_path(myb_nbody* nbody, const double* x, double* a,
                                                myb_internal_nbody_path path) {
  MYB_CHECK(nbody != NULL && x != NULL && a != NULL);
  myb_nbody* s = nbody;
  graph_update(s);
  for (size_t i = 0; i < s->count; ++i) {
    s->px[i] = x[3 * i];
    s->py[i] = x[3 * i + 1];
    s->pz[i] = x[3 * i + 2];
  }
  if (path == MYB_NBODY_PATH_AUTO) {
    size_t pairs = 0;
    const int huge = ckd_mul(&pairs, s->count, s->massive);
    const int pooled = s->pool != NULL && myb_pool_threads(s->pool) > 1;
    const int parallel = pooled && (huge || pairs >= PARALLEL_MIN_PAIRS);
    /* The vector blocks beat the symmetric path from a few dozen bodies on, although they
     * compute every pair twice (measured, x86-64-v2, 2000 bodies: 12.2 against 15.8 ns per
     * pair); the symmetric path stays for small systems. The bits are the same either way. */
    path = (parallel || s->excluded_count != 0 || s->massive >= VECTOR_MIN_BODIES)
               ? MYB_NBODY_PATH_RECEIVER
               : MYB_NBODY_PATH_SYMMETRIC;
  }
  if (path == MYB_NBODY_PATH_SYMMETRIC && s->excluded_count == 0) {
    sum_symmetric(s, a);
  } else {
    receiver_job job = {.s = s, .a = a, .vector = path != MYB_NBODY_PATH_RECEIVER_SCALAR};
    const size_t grain = s->massive >= PART_PAIRS ? 1 : PART_PAIRS / s->massive;
    myb_internal_pool_for(s->pool, s->count, grain, receiver_range, &job);
  }
  for (size_t i = 0; i < 3 * s->count; ++i) {
    if (!finite(a[i])) {
      return MYB_E_DOMAIN;
    }
  }
  return MYB_OK;
}

myb_status myb_nbody_acceleration(void* nbody, double t, const double* x, double* a) {
  (void)t;
  return myb_internal_nbody_acceleration_path(nbody, x, a, MYB_NBODY_PATH_AUTO);
}

myb_status myb_nbody_acceleration_xv(void* nbody, double t, const double* x, const double* v,
                                     double* a) {
  (void)t;
  (void)v;
  return myb_internal_nbody_acceleration_path(nbody, x, a, MYB_NBODY_PATH_AUTO);
}

/* ---- Invariants ---------------------------------------------------------------------------- */

static myb_dd dd_product(double a, double b) {
  double err = 0.0;
  const double p = myb_eft_two_prod(a, b, &err);
  return (myb_dd){.hi = p, .lo = err};
}

static myb_dd dd_difference(double a, double b) {
  double err = 0.0;
  const double d = myb_eft_two_sum(a, -b, &err);
  return (myb_dd){.hi = d, .lo = err};
}

/* a.b for double-word 3-vectors. */
static myb_dd dd_dot(const myb_dd a[3], const myb_dd b[3]) {
  myb_dd sum = myb_dd_mul(a[0], b[0]);
  sum = myb_dd_add(sum, myb_dd_mul(a[1], b[1]));
  return myb_dd_add(sum, myb_dd_mul(a[2], b[2]));
}

myb_status myb_nbody_invariants_of(const myb_nbody* nbody, const double* x, const double* v,
                                   myb_nbody_invariants* invariants) {
  MYB_CHECK(nbody != NULL && x != NULL && v != NULL && invariants != NULL);
  const myb_nbody* s = nbody;
  const size_t m = s->massive;
  for (size_t i = 0; i < 3 * m; ++i) {
    if (!finite(x[i]) || !finite(v[i])) {
      return MYB_E_DOMAIN;
    }
  }
  myb_dd kinetic = myb_dd_from_double(0.0);
  myb_dd potential = myb_dd_from_double(0.0);
  myb_dd gm_sum = myb_dd_from_double(0.0);
  myb_dd p[3];
  myb_dd l[3];
  myb_dd c[3];
  for (int k = 0; k < 3; ++k) {
    p[k] = myb_dd_from_double(0.0);
    l[k] = myb_dd_from_double(0.0);
    c[k] = myb_dd_from_double(0.0);
  }
  for (size_t i = 0; i < m; ++i) {
    const double g = s->gm[i];
    const double* xi = x + 3 * i;
    const double* vi = v + 3 * i;
    gm_sum = myb_dd_add_d(gm_sum, g);
    const myb_dd v2 = myb_dd_add(myb_dd_add(dd_product(vi[0], vi[0]), dd_product(vi[1], vi[1])),
                                 dd_product(vi[2], vi[2]));
    kinetic = myb_dd_add(kinetic, myb_dd_mul_d(v2, 0.5 * g));
    for (int k = 0; k < 3; ++k) {
      p[k] = myb_dd_add(p[k], dd_product(g, vi[k]));
      c[k] = myb_dd_add(c[k], dd_product(g, xi[k]));
      const int k1 = (k + 1) % 3;
      const int k2 = (k + 2) % 3;
      const myb_dd cross = myb_dd_sub(dd_product(xi[k1], vi[k2]), dd_product(xi[k2], vi[k1]));
      l[k] = myb_dd_add(l[k], myb_dd_mul_d(cross, g));
    }
    for (size_t j = i + 1; j < m; ++j) {
      const double* xj = x + 3 * j;
      const myb_dd d[3] = {dd_difference(xj[0], xi[0]), dd_difference(xj[1], xi[1]),
                           dd_difference(xj[2], xi[2])};
      const myb_dd r2 = dd_dot(d, d);
      if (!(r2.hi > 0.0)) {
        return MYB_E_DOMAIN;
      }
      potential = myb_dd_add(potential, myb_dd_div(dd_product(g, s->gm[j]), myb_dd_sqrt(r2)));
    }
  }
  if (!(gm_sum.hi > 0.0)) {
    return MYB_E_DOMAIN;
  }
  myb_nbody_invariants out = {.energy = myb_dd_sub(kinetic, potential)};
  for (int k = 0; k < 3; ++k) {
    out.momentum[k] = myb_dd_to_double(p[k]);
    out.angular_momentum[k] = myb_dd_to_double(l[k]);
    out.center[k] = myb_dd_to_double(myb_dd_div(c[k], gm_sum));
    out.center_velocity[k] = myb_dd_to_double(myb_dd_div(p[k], gm_sum));
  }
  if (!finite(out.energy.hi)) {
    return MYB_E_DOMAIN;
  }
  *invariants = out;
  return MYB_OK;
}
