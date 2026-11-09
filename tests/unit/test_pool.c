/* SPDX-License-Identifier: BSD-3-Clause */
/* myb_pool: every index of a parallel loop is visited exactly once, the number of parts follows
 * the grain and the thread count, and thousands of consecutive jobs on parked threads keep
 * working (start/finish handshake). */
#include <stdatomic.h>
#include <stdlib.h>

#include "core/fpguard.h"
#include "core/pool.h"
#include "mtest.h"
#include "myb/pool.h"

typedef struct marks {
  int* hits;
  atomic_int calls;
  atomic_long sum;
} marks;

static void mark(void* context, size_t begin, size_t end) {
  marks* m = context;
  atomic_fetch_add(&m->calls, 1);
  long sum = 0;
  for (size_t i = begin; i < end; ++i) {
    m->hits[i]++;
    sum += (long)i;
  }
  atomic_fetch_add(&m->sum, sum);
}

static void test_create(void) {
  myb_pool* pool = NULL;
  MT_CHECK(myb_pool_create(-1, NULL, &pool) == MYB_E_DOMAIN);
  MT_CHECK(myb_pool_create(1025, NULL, &pool) == MYB_E_DOMAIN);
  MT_CHECK(myb_pool_create(0, NULL, &pool) == MYB_OK);
  printf("pool of the machine: %d threads\n", myb_pool_threads(pool));
  MT_CHECK(myb_pool_threads(pool) >= 1);
  myb_pool_destroy(pool);
  myb_pool_destroy(NULL);
}

static void test_cover(void) {
  const int thread_counts[] = {1, 2, 3, 4, 7, 16};
  const size_t sizes[] = {0, 1, 2, 5, 100, 1001};
  const size_t grains[] = {1, 7, 1000};
  int* hits = calloc(1001, sizeof(int));
  MT_CHECK(hits != NULL);
  for (size_t t = 0; t < sizeof thread_counts / sizeof thread_counts[0]; ++t) {
    myb_pool* pool = NULL;
    MT_CHECK(myb_pool_create(thread_counts[t], NULL, &pool) == MYB_OK);
    MT_CHECK(myb_pool_threads(pool) == thread_counts[t]);
    for (size_t s = 0; s < sizeof sizes / sizeof sizes[0]; ++s) {
      for (size_t g = 0; g < sizeof grains / sizeof grains[0]; ++g) {
        const size_t count = sizes[s];
        for (size_t i = 0; i < count; ++i) {
          hits[i] = 0;
        }
        marks m = {.hits = hits};
        atomic_init(&m.calls, 0);
        atomic_init(&m.sum, 0);
        myb_internal_pool_for(pool, count, grains[g], mark, &m);
        int once = 1;
        for (size_t i = 0; i < count; ++i) {
          once = once && hits[i] == 1;
        }
        MT_CHECK(once);
        const size_t wanted = (count + grains[g] - 1) / grains[g];
        const size_t parts = wanted < (size_t)thread_counts[t] ? wanted : (size_t)thread_counts[t];
        const int expected = count == 0 ? 0 : (parts <= 1 ? 1 : (int)parts);
        MT_CHECK(atomic_load(&m.calls) == expected);
        MT_CHECK(atomic_load(&m.sum) == (long)(count * (count - (count > 0)) / 2) ||
                 count == 0);
      }
    }
    myb_pool_destroy(pool);
  }
  /* No pool: one call on the calling thread. */
  marks m = {.hits = hits};
  atomic_init(&m.calls, 0);
  atomic_init(&m.sum, 0);
  for (size_t i = 0; i < 50; ++i) {
    hits[i] = 0;
  }
  myb_internal_pool_for(NULL, 50, 1, mark, &m);
  MT_CHECK(atomic_load(&m.calls) == 1 && hits[49] == 1);
  free(hits);
}

static void test_many_jobs(void) {
  myb_pool* pool = NULL;
  MT_CHECK(myb_pool_create(4, NULL, &pool) == MYB_OK);
  int* hits = calloc(64, sizeof(int));
  MT_CHECK(hits != NULL);
  marks m = {.hits = hits};
  atomic_init(&m.calls, 0);
  atomic_init(&m.sum, 0);
  for (int job = 0; job < 5000; ++job) {
    myb_internal_pool_for(pool, 64, 1, mark, &m);
  }
  int all = 1;
  for (int i = 0; i < 64; ++i) {
    all = all && hits[i] == 5000;
  }
  MT_CHECK(all);
  MT_CHECK(atomic_load(&m.calls) == 4 * 5000);
  MT_CHECK(atomic_load(&m.sum) == 5000L * (63 * 64 / 2));
  free(hits);
  myb_pool_destroy(pool);
}

int main(void) {
  test_create();
  test_cover();
  test_many_jobs();
  return mt_finish("test_pool");
}
