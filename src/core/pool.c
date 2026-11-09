/* SPDX-License-Identifier: BSD-3-Clause */
#include "core/pool.h"

#include "core/alloc.h"
#include "core/base.h"
#include "platform/platform.h"

enum { MAX_THREADS = 1024 };

struct myb_pool {
  myb_allocator allocator;
  int threads;
  myb_platform_team* team;
};

myb_status myb_pool_create(int threads, const myb_allocator* allocator, myb_pool** pool) {
  MYB_CHECK(pool != NULL);
  if (threads == 0) {
    threads = myb_platform_processor_count();
    threads = threads > MAX_THREADS ? MAX_THREADS : threads;
  }
  if (threads < 1 || threads > MAX_THREADS) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_pool), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_pool* p = block;
  *p = (myb_pool){.allocator = a, .threads = threads};
  status = myb_platform_team_create(threads, &p->allocator, &p->team);
  if (status != MYB_OK) {
    goto fail;
  }
  *pool = p;
  return MYB_OK;

fail:
  myb_internal_release(&a, p);
  return status;
}

void myb_pool_destroy(myb_pool* pool) {
  if (pool == NULL) {
    return;
  }
  const myb_allocator a = pool->allocator;
  myb_platform_team_destroy(pool->team);
  myb_internal_release(&a, pool);
}

int myb_pool_threads(const myb_pool* pool) {
  MYB_CHECK(pool != NULL);
  return pool->threads;
}

typedef struct range_job {
  myb_internal_range_function* body;
  void* context;
  size_t count;
  size_t parts;
} range_job;

/* Part k of [0, count) cut into parts nearly equal pieces (the first count % parts are one
 * longer); no product that could overflow. */
static void run_part(void* context, int member) {
  const range_job* job = context;
  const size_t k = (size_t)member;
  if (k >= job->parts) {
    return;
  }
  const size_t base = job->count / job->parts;
  const size_t extra = job->count % job->parts;
  const size_t begin = k * base + (k < extra ? k : extra);
  const size_t end = begin + base + (k < extra ? 1 : 0);
  job->body(job->context, begin, end);
}

void myb_internal_pool_for(myb_pool* pool, size_t count, size_t grain,
                           myb_internal_range_function* body, void* context) {
  MYB_CHECK(body != NULL && grain > 0);
  if (count == 0) {
    return;
  }
  const size_t wanted = count / grain + (count % grain != 0 ? 1 : 0);
  const size_t threads = pool == NULL ? 1 : (size_t)pool->threads;
  const size_t parts = wanted < threads ? wanted : threads;
  if (parts <= 1) {
    body(context, 0, count);
    return;
  }
  range_job job = {.body = body, .context = context, .count = count, .parts = parts};
  myb_platform_team_run(pool->team, run_part, &job);
}
