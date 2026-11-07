/* SPDX-License-Identifier: BSD-3-Clause */
/* platform/platform.h — the only interface between myb and the operating system.
 *
 * Implementations: platform_win32.c (Windows), platform_posix.c (Linux and other POSIX).
 * No other file of the library includes an OS header. */
#ifndef MYB_SRC_PLATFORM_PLATFORM_H_
#define MYB_SRC_PLATFORM_PLATFORM_H_

#include <stddef.h>

#include "myb/core.h"

/* Reads a whole file into memory obtained from allocator (already resolved, not NULL).
 * path is UTF-8 on every platform (Windows converts it to UTF-16, so non-ASCII paths work).
 * The buffer gets a terminating NUL that is not counted in *size.
 * Errors: MYB_E_IO (cannot open or read), MYB_E_OVERFLOW (larger than max_size),
 * MYB_E_NOMEM. Outputs are untouched on error. */
[[nodiscard]] myb_status myb_platform_read_file(const char* path, const myb_allocator* allocator,
                                                size_t max_size, char** data, size_t* size);

/* Writes bytes to the standard error descriptor with one unbuffered call; best effort, used
 * only on the fatal-error path. Must not allocate and must not take stdio locks. */
void myb_platform_write_stderr(const char* data, size_t size);

/* Terminates the process immediately, bypassing exception handlers and atexit hooks. */
[[noreturn]] void myb_platform_fatal_exit(void);

#endif /* MYB_SRC_PLATFORM_PLATFORM_H_ */
