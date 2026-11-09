/* SPDX-License-Identifier: BSD-3-Clause */
/* core/pool.h — the data-parallel loop of myb over a myb_pool. */
#ifndef MYB_SRC_CORE_POOL_H_
#define MYB_SRC_CORE_POOL_H_

#include <stddef.h>

#include "myb/pool.h"

typedef void myb_internal_range_function(void* context, size_t begin, size_t end);

/* Calls body(context, begin, end) on consecutive parts that cover [0, count) exactly once,
 * in parallel on the threads of pool, and returns when all parts are done. The number of
 * parts is min(threads, ceil(count / grain)); pool == NULL, one thread or one part means a
 * single call on the calling thread. Parts must be independent: the caller guarantees that
 * the result does not depend on how [0, count) is cut. */
void myb_internal_pool_for(myb_pool* pool, size_t count, size_t grain,
                           myb_internal_range_function* body, void* context);

#endif /* MYB_SRC_CORE_POOL_H_ */
