/* SPDX-License-Identifier: BSD-3-Clause */
/* core/check.h — fail-fast checks of contracts and invariants.
 *
 *   MYB_CHECK(cond)    — checked in every build, including release. For preconditions of the
 *                        public API and invariants whose violation corrupts results.
 *   MYB_DCHECK(cond)   — checked only without NDEBUG; the expression is always compiled.
 *   MYB_UNREACHABLE()  — a branch that cannot be taken: terminates, never undefined behaviour.
 *
 * On failure one line "file:line: function: CHECK failed: condition" is written to stderr
 * with a single low-level write (no stdio lock, no allocation) and the process terminates
 * immediately (__fastfail on Windows, abort() elsewhere). Expected errors are never reported
 * this way: they are returned as myb_status. */
#ifndef MYB_SRC_CORE_CHECK_H_
#define MYB_SRC_CORE_CHECK_H_

#include "myb/core.h"

/* Exported from the shared library (MYB_API) for the libraries of ext/, which are built with
 * the same checks; not part of the public API. */
[[noreturn]] [[gnu::cold]] MYB_API void myb_internal_check_failed(const char* file, int line,
                                                                  const char* func,
                                                                  const char* what);

#define MYB_CHECK(cond)                                               \
  do {                                                                \
    if (!__builtin_expect(!!(cond), 1)) {                             \
      myb_internal_check_failed(__FILE__, __LINE__, __func__, #cond); \
    }                                                                 \
  } while (0)

#if defined(NDEBUG)
#define MYB_DCHECK_ENABLED 0
#else
#define MYB_DCHECK_ENABLED 1
#endif

#define MYB_DCHECK(cond)                                              \
  do {                                                                \
    if (MYB_DCHECK_ENABLED && !(cond)) {                              \
      myb_internal_check_failed(__FILE__, __LINE__, __func__, #cond); \
    }                                                                 \
  } while (0)

#define MYB_UNREACHABLE() myb_internal_check_failed(__FILE__, __LINE__, __func__, "unreachable")

#endif /* MYB_SRC_CORE_CHECK_H_ */
