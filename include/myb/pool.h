/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/pool.h — a fixed set of threads for the data-parallel loops of myb (N-body forces).
 *
 * Results never depend on the pool: work is split so that every output element is computed by
 * one thread with the same operations in the same order as without a pool, so the bits are
 * identical for any number of threads (checked by tests/repro). The threads are started once
 * and parked between jobs. A pool runs one job at a time: objects that share a pool must not
 * be used from several threads concurrently (a violation terminates the process). */
#ifndef MYB_POOL_H_
#define MYB_POOL_H_

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct myb_pool myb_pool;

/* threads: threads working on a job, the calling thread included (1 = no extra thread,
 * 0 = one per logical processor), at most 1024. Errors: MYB_E_DOMAIN (threads out of range),
 * MYB_E_NOMEM, MYB_E_IO (the system refused to start a thread). */
MYB_API MYB_NODISCARD myb_status myb_pool_create(int threads, const myb_allocator* allocator,
                                                 myb_pool** pool);

/* Stops the threads; accepts NULL. */
MYB_API void myb_pool_destroy(myb_pool* pool);

/* Threads working on a job, the calling thread included. */
MYB_API int myb_pool_threads(const myb_pool* pool);

#ifdef __cplusplus
}
#endif

#endif /* MYB_POOL_H_ */
