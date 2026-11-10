/* SPDX-License-Identifier: BSD-3-Clause */
/* core/alloc.h — allocation through a myb_allocator (NULL = malloc/free) with overflow-checked
 * array sizes, and recording of diagnostics. */
#ifndef MYB_SRC_CORE_ALLOC_H_
#define MYB_SRC_CORE_ALLOC_H_

#include <stddef.h>

#include "myb/core.h"

/* The functions of this header are exported from the shared library (MYB_API) for the
 * libraries of ext/, which allocate and report the same way; not part of the public API. */

/* The allocator to use: *requested, or the C runtime's malloc/free when requested is NULL. */
MYB_API myb_allocator myb_internal_allocator(const myb_allocator* requested);

/* count * size bytes (count == 0 gives a valid 1-byte block). MYB_E_OVERFLOW if the size does
 * not fit in size_t, MYB_E_NOMEM if the allocator fails. *out is untouched on error. */
[[nodiscard]] MYB_API myb_status myb_internal_alloc_array(const myb_allocator* allocator,
                                                          size_t count, size_t size, void** out);

MYB_API void myb_internal_release(const myb_allocator* allocator, void* pointer);

/* Records a diagnostic (printf-style) if diag is not NULL; the message is truncated to fit. */
[[gnu::format(printf, 3, 4)]] MYB_API void myb_internal_diag(myb_diag* diag, long line,
                                                             const char* format, ...);

/* Clears a diagnostic (empty message, line 0) if diag is not NULL. */
MYB_API void myb_internal_diag_clear(myb_diag* diag);

#endif /* MYB_SRC_CORE_ALLOC_H_ */
