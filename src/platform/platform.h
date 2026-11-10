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

/* Creates the file at path (UTF-8), or replaces it, with the size bytes of data. Errors:
 * MYB_E_IO (cannot create or write: what was written of a failed attempt is left behind),
 * MYB_E_NOMEM (the path, on Windows). */
[[nodiscard]] myb_status myb_platform_write_file(const char* path, const myb_allocator* allocator,
                                                 const void* data, size_t size);

/* A read-only view of a whole file, mapped into the address space: large data files (planetary
 * ephemerides reach gigabytes) are read on demand by the virtual memory system, without a
 * system call per access. The view stays valid until myb_platform_unmap_file; the file must
 * not be truncated meanwhile (an access beyond the new end is a fatal fault, by design). */
typedef struct myb_platform_mapping {
  const unsigned char* data;
  size_t size;
} myb_platform_mapping;

/* path is UTF-8 on every platform. Errors: MYB_E_IO (cannot open or map, empty file),
 * MYB_E_OVERFLOW (does not fit the address space). Output untouched on error. */
[[nodiscard]] myb_status myb_platform_map_file(const char* path, const myb_allocator* allocator,
                                               myb_platform_mapping* mapping);

/* Releases a view from myb_platform_map_file; accepts a zeroed mapping. */
void myb_platform_unmap_file(myb_platform_mapping* mapping);

/* Writes bytes to the standard error descriptor with one unbuffered call; best effort, used
 * only on the fatal-error path. Must not allocate and must not take stdio locks. */
void myb_platform_write_stderr(const char* data, size_t size);

/* Terminates the process immediately, bypassing exception handlers and atexit hooks. */
[[noreturn]] void myb_platform_fatal_exit(void);

/* ---- Threads ----------------------------------------------------------------------------------
 * A team of members threads: member 0 is the thread that calls myb_platform_team_run, members
 * 1 .. members-1 are worker threads started once by create and parked on a condition variable
 * between jobs (no thread is created or destroyed per job). */
typedef struct myb_platform_team myb_platform_team;
typedef void myb_platform_job(void* context, int member);

/* members >= 1 (1 starts no thread). allocator is resolved (not NULL) and must outlive the
 * team. Errors: MYB_E_NOMEM, MYB_E_IO (the system refused to start a thread). */
[[nodiscard]] myb_status myb_platform_team_create(int members, const myb_allocator* allocator,
                                                  myb_platform_team** team);

/* Stops and joins the workers; accepts NULL. Must not be called while a job runs. */
void myb_platform_team_destroy(myb_platform_team* team);

/* Runs job(context, k) for every member k exactly once, k = 0 on the calling thread, and
 * returns when all have returned. Everything the jobs wrote happens before the return. One job
 * at a time per team (a second concurrent call is a contract violation). */
void myb_platform_team_run(myb_platform_team* team, myb_platform_job* job, void* context);

/* Logical processors available to the process (>= 1). */
int myb_platform_processor_count(void);

#endif /* MYB_SRC_PLATFORM_PLATFORM_H_ */
