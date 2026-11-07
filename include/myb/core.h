/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/core.h — version, status codes and attribute macros shared by every myb header.
 *
 * Public headers of myb are valid C17, C23 and C++17: the library itself is built with
 * Clang in C23 mode, but it can be used from any compiler. Headers contain declarations
 * only — no inline arithmetic — so every floating-point operation of myb is compiled with
 * the library's own flags (no contraction, no fast-math), never with the caller's flags.
 *
 * Numerical environment contract: myb assumes IEEE-754 binary64 arithmetic in the
 * default rounding mode (round to nearest, ties to even). It never changes the
 * floating-point environment; callers must not call it with another rounding mode. */
#ifndef MYB_CORE_H_
#define MYB_CORE_H_

#include <stddef.h>
#include <stdint.h>

/* ---- Version (single source of truth; CMake reads it from here) ---------------------- */
#define MYB_VERSION_MAJOR 0
#define MYB_VERSION_MINOR 1
#define MYB_VERSION_PATCH 0
#define MYB_VERSION_NUMBER \
  ((MYB_VERSION_MAJOR * 10000L) + (MYB_VERSION_MINOR * 100L) + MYB_VERSION_PATCH)

/* ---- Symbol visibility ------------------------------------------------------------------
 * MYB_STATIC  — the library is linked statically (set by the CMake target).
 * MYB_BUILDING — defined only while compiling the shared library itself. */
#if defined(MYB_STATIC)
#define MYB_API
#elif defined(_WIN32)
#if defined(MYB_BUILDING)
#define MYB_API __declspec(dllexport)
#else
#define MYB_API __declspec(dllimport)
#endif
#elif defined(__GNUC__) || defined(__clang__)
#define MYB_API __attribute__((visibility("default")))
#else
#define MYB_API
#endif

/* ---- Attributes ----------------------------------------------------------------------- */
#if defined(__cplusplus)
#if __cplusplus >= 201703L || (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L)
#define MYB_NODISCARD [[nodiscard]]
#endif
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ > 201710L
#define MYB_NODISCARD [[nodiscard]]
#endif
#if !defined(MYB_NODISCARD)
#if defined(__GNUC__) || defined(__clang__)
#define MYB_NODISCARD __attribute__((warn_unused_result))
#else
#define MYB_NODISCARD
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Status codes -----------------------------------------------------------------------
 * Expected errors (bad input data, coverage limits, non-convergence) are returned as
 * myb_status. A violated precondition of the API (NULL where an object is required,
 * inconsistent sizes) is a bug in the caller: myb terminates the process immediately with a
 * one-line message on stderr instead of continuing with corrupted state.
 * Values are stable and never reused. */
typedef enum myb_status {
  MYB_OK = 0,
  MYB_E_DOMAIN = 1,   /* argument outside the mathematical domain of the model */
  MYB_E_RANGE = 2,    /* epoch or argument outside the coverage of the data */
  MYB_E_NOCONV = 3,   /* iteration budget exhausted without convergence */
  MYB_E_STEP_MIN = 4, /* integrator step fell below the allowed minimum */
  MYB_E_NOMEM = 5,    /* allocation failed */
  MYB_E_OVERFLOW = 6, /* size or index computation overflowed */
  MYB_E_IO = 7,       /* file could not be read */
  MYB_E_FORMAT = 8,   /* malformed input data */
  MYB_E_NOT_FOUND = 9 /* named item (data variable, time scale, ...) does not exist */
} myb_status;

/* Stable ASCII name of a status code ("MYB_E_RANGE"); "MYB_E_UNKNOWN" for other values. */
MYB_API const char* myb_status_name(myb_status status);

/* ---- Diagnostics -------------------------------------------------------------------------
 * Functions that read external data accept an optional myb_diag*: on an error they record
 * where and why (one line, ASCII), so that a broken data file can be fixed without a
 * debugger. Passing NULL is always allowed. */
typedef struct myb_diag {
  long line;         /* 1-based line of the input, 0 if not applicable */
  char message[200]; /* NUL-terminated description, empty on success */
} myb_diag;

/* ---- Memory -------------------------------------------------------------------------------
 * Every object of myb that owns memory takes an allocator; NULL means the C runtime's
 * malloc/free. alloc returns memory aligned for any object type (like malloc) or NULL;
 * release accepts NULL. Allocation never happens inside numerical inner loops. */
typedef struct myb_allocator {
  void* (*alloc)(void* context, size_t size);
  void (*release)(void* context, void* pointer);
  void* context;
} myb_allocator;

/* Version of the linked library in the format of MYB_VERSION_NUMBER. Comparing it with the
 * macro detects a header/library mismatch. */
MYB_API long myb_version_number(void);

#ifdef __cplusplus
}
#endif

#endif /* MYB_CORE_H_ */
