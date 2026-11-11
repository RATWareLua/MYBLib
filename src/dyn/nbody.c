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
 * excluded sources (and itself), so the inner loops have no branch.
 *
 * With the 1PN terms (myb_nbody_set_relativity) two more passes follow the Newtonian one, both
 * per receiver in increasing source order: the potentials U_i, then the 1PN terms, which need
 * the Newtonian accelerations and potentials of the sources. The 1PN terms enter their own
 * Sum2 accumulator and are added to the Newtonian acceleration once. */
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

/* A pair whose level was set explicitly to NEWTON or EIH (NONE pairs are exclusions). */
typedef struct pn_edge {
  size_t receiver;
  size_t source;
  int level;
} pn_edge;

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
  /* 1PN */
  int relativity; /* parameters set: EIH is the default level */
  myb_ppn ppn;
  pn_edge* explicit_pn; /* sorted by (receiver, source) */
  size_t explicit_count;
  size_t explicit_capacity;
  size_t* pn_first;  /* count + 1, as first[] for explicit_pn */
  double* potential; /* count: U of each body (massive sources) */
  double* newton;    /* 3 count: Newtonian accelerations of the current evaluation */
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
  if (ckd_mul(&doubles, count, (size_t)7) || ckd_add(&doubles, doubles, massive) ||
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
  s->potential = s->pz + count;
  s->newton = s->potential + count;
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
  status = myb_internal_alloc_array(&a, firsts, sizeof(size_t), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  s->pn_first = block;
  *nbody = s;
  return MYB_OK;

fail:
  myb_internal_release(&a, s->pn_first);
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
  myb_internal_release(&a, nbody->explicit_pn);
  myb_internal_release(&a, nbody->pn_first);
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

void myb_internal_nbody_set_gm(myb_nbody* nbody, size_t body, double gm) {
  MYB_CHECK(nbody != NULL && body < nbody->massive && gm >= 0.0 && gm <= 0x1p1023);
  nbody->gm[body] = gm;
}

int myb_internal_nbody_relativity(const myb_nbody* nbody, myb_ppn* ppn) {
  MYB_CHECK(nbody != NULL && ppn != NULL);
  if (!nbody->relativity) {
    return 0;
  }
  *ppn = nbody->ppn;
  return 1;
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

/* Removes the exclusion (receiver, source) if present. */
static void exclusion_remove(myb_nbody* s, size_t receiver, size_t source) {
  const size_t at = exclusion_lower_bound(s, receiver, source);
  if (exclusion_found(s, at, receiver, source)) {
    memmove(&s->excluded[at], &s->excluded[at + 1],
            (s->excluded_count - at - 1) * sizeof(exclusion));
    s->excluded_count--;
    s->graph_valid = 0;
  }
}

static myb_status exclusion_add(myb_nbody* s, size_t receiver, size_t source) {
  const size_t at = exclusion_lower_bound(s, receiver, source);
  if (exclusion_found(s, at, receiver, source)) {
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

/* First explicit 1PN-level edge whose key is not less than (receiver, source). */
static size_t pn_lower_bound(const myb_nbody* s, size_t receiver, size_t source) {
  size_t lo = 0;
  size_t hi = s->explicit_count;
  /* bound: the interval halves, at most 64 iterations */
  while (lo < hi) {
    const size_t mid = lo + (hi - lo) / 2;
    const pn_edge* e = &s->explicit_pn[mid];
    if (e->receiver < receiver || (e->receiver == receiver && e->source < source)) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

/* Records (or, for level < 0, forgets) the explicit level of a pair. */
static myb_status pn_set(myb_nbody* s, size_t receiver, size_t source, int level) {
  const size_t at = pn_lower_bound(s, receiver, source);
  const int present = at < s->explicit_count && s->explicit_pn[at].receiver == receiver &&
                      s->explicit_pn[at].source == source;
  if (level < 0) {
    if (present) {
      memmove(&s->explicit_pn[at], &s->explicit_pn[at + 1],
              (s->explicit_count - at - 1) * sizeof(pn_edge));
      s->explicit_count--;
      s->graph_valid = 0;
    }
    return MYB_OK;
  }
  if (present) {
    s->explicit_pn[at].level = level;
    return MYB_OK;
  }
  if (s->explicit_count == s->explicit_capacity) {
    size_t capacity = s->explicit_capacity == 0 ? 16 : s->explicit_capacity;
    if (s->explicit_capacity != 0 && ckd_mul(&capacity, capacity, (size_t)2)) {
      return MYB_E_OVERFLOW;
    }
    void* block = NULL;
    const myb_status status =
        myb_internal_alloc_array(&s->allocator, capacity, sizeof(pn_edge), &block);
    if (status != MYB_OK) {
      return status;
    }
    if (s->explicit_count != 0) {
      memcpy(block, s->explicit_pn, s->explicit_count * sizeof(pn_edge));
    }
    myb_internal_release(&s->allocator, s->explicit_pn);
    s->explicit_pn = block;
    s->explicit_capacity = capacity;
  }
  memmove(&s->explicit_pn[at + 1], &s->explicit_pn[at],
          (s->explicit_count - at) * sizeof(pn_edge));
  s->explicit_pn[at] = (pn_edge){.receiver = receiver, .source = source, .level = level};
  s->explicit_count++;
  s->graph_valid = 0;
  return MYB_OK;
}

myb_status myb_nbody_set_interaction(myb_nbody* nbody, size_t receiver, size_t source,
                                     myb_interaction level) {
  MYB_CHECK(nbody != NULL);
  myb_nbody* s = nbody;
  if (receiver >= s->count || source >= s->massive || receiver == source ||
      (level != MYB_INTERACTION_NONE && level != MYB_INTERACTION_NEWTON &&
       level != MYB_INTERACTION_EIH) ||
      (level == MYB_INTERACTION_EIH && !s->relativity)) {
    return MYB_E_DOMAIN;
  }
  if (level == MYB_INTERACTION_NONE) {
    const myb_status status = exclusion_add(s, receiver, source);
    if (status != MYB_OK) {
      return status;
    }
    return pn_set(s, receiver, source, -1);
  }
  /* reserve the explicit edge first, so that a failure leaves the graph unchanged */
  const myb_status status = pn_set(s, receiver, source, (int)level);
  if (status != MYB_OK) {
    return status;
  }
  exclusion_remove(s, receiver, source);
  return MYB_OK;
}

myb_status myb_nbody_set_relativity(myb_nbody* nbody, const myb_ppn* ppn) {
  MYB_CHECK(nbody != NULL);
  myb_nbody* s = nbody;
  if (ppn == NULL) {
    /* back to Newtonian: pairs set to EIH become Newtonian */
    for (size_t k = 0; k < s->explicit_count; ++k) {
      if (s->explicit_pn[k].level == MYB_INTERACTION_EIH) {
        s->explicit_pn[k].level = MYB_INTERACTION_NEWTON;
      }
    }
    s->relativity = 0;
    return MYB_OK;
  }
  if (myb_ppn_validate(ppn) != MYB_OK) {
    return MYB_E_DOMAIN;
  }
  s->ppn = *ppn;
  s->relativity = 1;
  return MYB_OK;
}

int myb_internal_nbody_velocity_dependent(const myb_nbody* nbody) {
  MYB_CHECK(nbody != NULL);
  return nbody->relativity;
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
  size_t p = 0;
  for (size_t i = 0; i <= s->count; ++i) {
    /* bound: p only grows, to explicit_count */
    while (p < s->explicit_count && s->explicit_pn[p].receiver < i) {
      ++p;
    }
    s->pn_first[i] = p;
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

/* Vectors never cross a function boundary by value: on x86-64 SysV a 32-byte vector argument
 * would take a different ABI with and without AVX (-Wpsabi). */
#define MYB_SPLAT(x) ((v4d){(x), (x), (x), (x)})

static inline void accumulate4(accumulator4* a, int k, const v4d* x) {
  const v4d s = a->s[k] + *x;
  const v4d a_rounded = s - *x;
  const v4d b_rounded = s - a_rounded;
  const v4d a_error = a->s[k] - a_rounded;
  const v4d b_error = *x - b_rounded;
  a->s[k] = s;
  a->c[k] += a_error + b_error;
}

/* Sources [begin, end) on the four receivers at (xi, yi, zi); lanes whose receiver index
 * equals the source are masked when mask is set (receivers i0 .. i0 + 3). */
static void sum_sources4(const myb_nbody* s, const v4d receiver[3], size_t begin, size_t end,
                         int mask, size_t i0, accumulator4* acc) {
  const v4d xi = receiver[0];
  const v4d yi = receiver[1];
  const v4d zi = receiver[2];
  const double* gm = s->gm;
  const double* px = s->px;
  const double* py = s->py;
  const double* pz = s->pz;
  const v4i lane_index = {(long long)i0, (long long)i0 + 1, (long long)i0 + 2,
                          (long long)i0 + 3};
  for (size_t j = begin; j < end; ++j) {
    const v4d dx = MYB_SPLAT(px[j]) - xi;
    const v4d dy = MYB_SPLAT(py[j]) - yi;
    const v4d dz = MYB_SPLAT(pz[j]) - zi;
    const v4d r2 = dx * dx + dy * dy + dz * dz;
    const v4d q = MYB_SPLAT(1.0) / (r2 * __builtin_elementwise_sqrt(r2));
    const v4d sj = MYB_SPLAT(gm[j]) * q;
    v4d t[3] = {sj * dx, sj * dy, sj * dz};
    if (mask) {
      const v4i keep = lane_index != (v4i){(long long)j, (long long)j, (long long)j, (long long)j};
      for (int k = 0; k < 3; ++k) {
        t[k] = (v4d)((v4i)t[k] & keep);
      }
    }
    for (int k = 0; k < 3; ++k) {
      accumulate4(acc, k, &t[k]);
    }
  }
}

/* Receivers i0 .. i0 + 3, none of which has an exclusion. */
static void sum_receivers4(const myb_nbody* s, size_t i0, double* a) {
  const v4d receiver[3] = {{s->px[i0], s->px[i0 + 1], s->px[i0 + 2], s->px[i0 + 3]},
                           {s->py[i0], s->py[i0 + 1], s->py[i0 + 2], s->py[i0 + 3]},
                           {s->pz[i0], s->pz[i0 + 1], s->pz[i0 + 2], s->pz[i0 + 3]}};
  accumulator4 acc = {};
  const size_t m = s->massive;
  const size_t lo = i0 < m ? i0 : m;          /* sources below every lane: no mask */
  const size_t hi = i0 + LANES < m ? i0 + LANES : m; /* [lo, hi): one lane is the source */
  sum_sources4(s, receiver, 0, lo, 0, i0, &acc);
  sum_sources4(s, receiver, lo, hi, 1, i0, &acc);
  sum_sources4(s, receiver, hi, m, 0, i0, &acc);
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

/* ---- 1PN terms ------------------------------------------------------------------------------ */

typedef struct pn_job {
  const myb_nbody* s;
  const double* v;
  double* a;
} pn_job;

/* U_i = sum of GM_j / r_ij over the massive sources of i (Sum2, increasing j). */
static void potential_range(void* context, size_t begin, size_t end) {
  const pn_job* job = context;
  const myb_nbody* s = job->s;
  for (size_t i = begin; i < end; ++i) {
    const double xi = s->px[i];
    const double yi = s->py[i];
    const double zi = s->pz[i];
    size_t e = s->excluded_count == 0 ? 0 : s->first[i];
    const size_t e_end = s->excluded_count == 0 ? 0 : s->first[i + 1];
    double sum = 0.0;
    double comp = 0.0;
    for (size_t j = 0; j < s->massive; ++j) {
      /* bound: e only grows, to e_end */
      while (e < e_end && s->excluded[e].source < j) {
        ++e;
      }
      if (j == i || (e < e_end && s->excluded[e].source == j)) {
        continue;
      }
      const double dx = s->px[j] - xi;
      const double dy = s->py[j] - yi;
      const double dz = s->pz[j] - zi;
      double err = 0.0;
      sum = myb_eft_two_sum(sum, s->gm[j] / sqrt(dx * dx + dy * dy + dz * dz), &err);
      comp += err;
    }
    s->potential[i] = sum + comp;
  }
}

/* a_i = Newtonian a_i + the 1PN terms of its EIH sources (Sum2, increasing j). */
static void pn_range(void* context, size_t begin, size_t end) {
  const pn_job* job = context;
  const myb_nbody* s = job->s;
  const double beta = s->ppn.beta;
  const double gamma = s->ppn.gamma;
  const double inv_c2 = 1.0 / (s->ppn.c * s->ppn.c);
  const int default_level = s->relativity ? MYB_INTERACTION_EIH : MYB_INTERACTION_NEWTON;
  for (size_t i = begin; i < end; ++i) {
    const double xi = s->px[i];
    const double yi = s->py[i];
    const double zi = s->pz[i];
    const double* vi = job->v + 3 * i;
    const double ui = s->potential[i];
    const double vi2 = vi[0] * vi[0] + vi[1] * vi[1] + vi[2] * vi[2];
    size_t e = s->excluded_count == 0 ? 0 : s->first[i];
    const size_t e_end = s->excluded_count == 0 ? 0 : s->first[i + 1];
    size_t p = s->explicit_count == 0 ? 0 : s->pn_first[i];
    const size_t p_end = s->explicit_count == 0 ? 0 : s->pn_first[i + 1];
    accumulator acc = {};
    for (size_t j = 0; j < s->massive; ++j) {
      /* bound: e and p only grow, to e_end and p_end */
      while (e < e_end && s->excluded[e].source < j) {
        ++e;
      }
      while (p < p_end && s->explicit_pn[p].source < j) {
        ++p;
      }
      if (j == i || (e < e_end && s->excluded[e].source == j)) {
        continue;
      }
      const int level =
          (p < p_end && s->explicit_pn[p].source == j) ? s->explicit_pn[p].level : default_level;
      if (level != MYB_INTERACTION_EIH) {
        continue;
      }
      const double dx = s->px[j] - xi;
      const double dy = s->py[j] - yi;
      const double dz = s->pz[j] - zi;
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double r = sqrt(r2);
      const double* vj = job->v + 3 * j;
      const double* aj = s->newton + 3 * j;
      const double vj2 = vj[0] * vj[0] + vj[1] * vj[1] + vj[2] * vj[2];
      const double vivj = vi[0] * vj[0] + vi[1] * vj[1] + vi[2] * vj[2];
      const double dvj = dx * vj[0] + dy * vj[1] + dz * vj[2];
      const double daj = dx * aj[0] + dy * aj[1] + dz * aj[2];
      const double f = -2.0 * (beta + gamma) * ui - (2.0 * beta - 1.0) * s->potential[j] +
                       gamma * vi2 + (1.0 + gamma) * vj2 - 2.0 * (1.0 + gamma) * vivj -
                       1.5 * dvj * dvj / r2 + 0.5 * daj;
      const double w0 = (2.0 + 2.0 * gamma) * vi[0] - (1.0 + 2.0 * gamma) * vj[0];
      const double w1 = (2.0 + 2.0 * gamma) * vi[1] - (1.0 + 2.0 * gamma) * vj[1];
      const double w2 = (2.0 + 2.0 * gamma) * vi[2] - (1.0 + 2.0 * gamma) * vj[2];
      const double g = -(dx * w0 + dy * w1 + dz * w2);
      const double s1 = s->gm[j] / (r2 * r);
      const double s3 = (1.5 + 2.0 * gamma) * s->gm[j] / r;
      const double tx = s1 * (f * dx + g * (vi[0] - vj[0])) + s3 * aj[0];
      const double ty = s1 * (f * dy + g * (vi[1] - vj[1])) + s3 * aj[1];
      const double tz = s1 * (f * dz + g * (vi[2] - vj[2])) + s3 * aj[2];
      accumulate(&acc, tx * inv_c2, ty * inv_c2, tz * inv_c2);
    }
    for (int k = 0; k < 3; ++k) {
      job->a[3 * i + (size_t)k] = s->newton[3 * i + (size_t)k] + (acc.s[k] + acc.c[k]);
    }
  }
}

myb_status myb_internal_nbody_acceleration_xv_path(myb_nbody* nbody, const double* x,
                                                   const double* v, double* a,
                                                   myb_internal_nbody_path path) {
  MYB_CHECK(nbody != NULL && x != NULL && a != NULL);
  myb_nbody* s = nbody;
  if (s->relativity && v == NULL) {
    return MYB_E_DOMAIN;
  }
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
  double* newtonian = s->relativity ? s->newton : a;
  if (path == MYB_NBODY_PATH_SYMMETRIC && s->excluded_count == 0) {
    sum_symmetric(s, newtonian);
  } else {
    receiver_job job = {.s = s, .a = newtonian, .vector = path != MYB_NBODY_PATH_RECEIVER_SCALAR};
    const size_t grain = s->massive >= PART_PAIRS ? 1 : PART_PAIRS / s->massive;
    myb_internal_pool_for(s->pool, s->count, grain, receiver_range, &job);
  }
  if (s->relativity) {
    for (size_t i = 0; i < 3 * s->count; ++i) {
      if (!finite(v[i])) {
        return MYB_E_DOMAIN;
      }
    }
    pn_job job = {.s = s, .v = v, .a = a};
    const size_t grain = s->massive >= PART_PAIRS ? 1 : PART_PAIRS / s->massive;
    myb_internal_pool_for(s->pool, s->count, grain, potential_range, &job);
    myb_internal_pool_for(s->pool, s->count, grain, pn_range, &job);
  }
  for (size_t i = 0; i < 3 * s->count; ++i) {
    if (!finite(a[i])) {
      return MYB_E_DOMAIN;
    }
  }
  return MYB_OK;
}

myb_status myb_internal_nbody_acceleration_path(myb_nbody* nbody, const double* x, double* a,
                                                myb_internal_nbody_path path) {
  return myb_internal_nbody_acceleration_xv_path(nbody, x, NULL, a, path);
}

myb_status myb_nbody_acceleration(void* nbody, double t, const double* x, double* a) {
  (void)t;
  return myb_internal_nbody_acceleration_xv_path(nbody, x, NULL, a, MYB_NBODY_PATH_AUTO);
}

myb_status myb_nbody_acceleration_xv(void* nbody, double t, const double* x, const double* v,
                                     double* a) {
  (void)t;
  MYB_CHECK(v != NULL);
  return myb_internal_nbody_acceleration_xv_path(nbody, x, v, a, MYB_NBODY_PATH_AUTO);
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

myb_status myb_nbody_pn_invariants_of(const myb_nbody* nbody, const double* x, const double* v,
                                      myb_nbody_pn_invariants* invariants) {
  MYB_CHECK(nbody != NULL && x != NULL && v != NULL && invariants != NULL);
  const myb_nbody* s = nbody;
  if (!s->relativity) {
    return MYB_E_DOMAIN;
  }
  myb_nbody_invariants newtonian;
  const myb_status status = myb_nbody_invariants_of(nbody, x, v, &newtonian);
  if (status != MYB_OK) {
    return status;
  }
  const size_t m = s->massive;
  const double inv_c2 = 1.0 / (s->ppn.c * s->ppn.c);
  const double beta = s->ppn.beta;
  const double gamma = s->ppn.gamma;
  myb_sum_acc energy = {0};
  myb_sum_acc momentum[3] = {{0}};
  myb_sum_acc center[3] = {{0}};
  myb_sum_acc weight = {0};
  for (size_t a = 0; a < m; ++a) {
    const double* xa = x + 3 * a;
    const double* va = v + 3 * a;
    const double ga = s->gm[a];
    const double va2 = va[0] * va[0] + va[1] * va[1] + va[2] * va[2];
    /* U_a and the pair terms of a */
    myb_sum_acc potential = {0};
    myb_sum_acc pair_energy = {0};
    myb_sum_acc pair_momentum[3] = {{0}};
    for (size_t b = 0; b < m; ++b) {
      if (b == a) {
        continue;
      }
      const double* xb = x + 3 * b;
      const double* vb = v + 3 * b;
      const double d[3] = {xa[0] - xb[0], xa[1] - xb[1], xa[2] - xb[2]};
      const double r = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      const double n[3] = {d[0] / r, d[1] / r, d[2] / r};
      const double k = s->gm[b] / r; /* GM_b / r_ab */
      const double nva = n[0] * va[0] + n[1] * va[1] + n[2] * va[2];
      const double nvb = n[0] * vb[0] + n[1] * vb[1] + n[2] * vb[2];
      const double vavb = va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2];
      myb_sum_acc_add(&potential, k);
      myb_sum_acc_add(&pair_energy, 0.5 * ga * k *
                                        ((2.0 * gamma + 1.0) * va2 -
                                         (2.0 * gamma + 1.5) * vavb - 0.5 * nva * nvb));
      for (int c = 0; c < 3; ++c) {
        myb_sum_acc_add(&pair_momentum[c], -0.5 * ga * k * nva * n[c]);
      }
    }
    const double ua = myb_sum_acc_value(&potential);
    const double factor = 1.0 + 0.5 * (va2 - ua) * inv_c2; /* GM*_a / GM_a */
    myb_sum_acc_add(&energy, 0.375 * ga * va2 * va2 * inv_c2);
    myb_sum_acc_add(&energy, 0.5 * (2.0 * beta - 1.0) * ga * ua * ua * inv_c2);
    myb_sum_acc_add(&energy, myb_sum_acc_value(&pair_energy) * inv_c2);
    myb_sum_acc_add(&weight, ga * factor);
    for (int c = 0; c < 3; ++c) {
      myb_sum_acc_add(&momentum[c], ga * factor * va[c]);
      myb_sum_acc_add(&momentum[c], myb_sum_acc_value(&pair_momentum[c]) * inv_c2);
      myb_sum_acc_add(&center[c], ga * factor * xa[c]);
    }
  }
  const double total = myb_sum_acc_value(&weight);
  if (!(total > 0.0)) {
    return MYB_E_DOMAIN;
  }
  myb_nbody_pn_invariants out = {
      .energy = myb_dd_add_d(newtonian.energy, myb_sum_acc_value(&energy))};
  for (int c = 0; c < 3; ++c) {
    out.momentum[c] = myb_sum_acc_value(&momentum[c]);
    out.center[c] = myb_sum_acc_value(&center[c]) / total;
    out.center_velocity[c] = out.momentum[c] / total;
  }
  if (!finite(out.energy.hi)) {
    return MYB_E_DOMAIN;
  }
  *invariants = out;
  return MYB_OK;
}
