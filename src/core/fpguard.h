/* SPDX-License-Identifier: BSD-3-Clause */
/* core/fpguard.h — compile-time guarantees of the numerical environment.
 *
 * Every translation unit of the library includes this header (through core/base.h). A build
 * that would silently break IEEE-754 semantics or cross-platform reproducibility stops here
 * with a message instead of producing different numbers. */
#ifndef MYB_SRC_CORE_FPGUARD_H_
#define MYB_SRC_CORE_FPGUARD_H_

#include <float.h>

#if !defined(__clang__)
#error "myb is built with Clang only (see docs/PLAN.md)"
#endif

#if !defined(__STDC_VERSION__) || __STDC_VERSION__ <= 201710L
#error "myb sources are C23: compile with -std=c23 (set by CMakeLists.txt)"
#endif

/* -ffast-math / -Ofast: reassociation, no NaN/Inf, flush to zero — wrong and
 * non-reproducible results. */
#if defined(__FAST_MATH__)
#error "myb: -ffast-math (or -Ofast) is forbidden: it breaks IEEE-754 semantics"
#endif
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__
#error "myb: -ffinite-math-only is forbidden: NaN and infinity must be honoured"
#endif

static_assert(FLT_RADIX == 2, "myb needs binary floating point");
static_assert(DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 && DBL_MIN_EXP == -1021,
              "myb needs IEEE-754 binary64 for double");
static_assert(sizeof(double) == 8, "myb needs a 64-bit double");
/* x87 extended evaluation (FLT_EVAL_METHOD 2) rounds twice and breaks error-free
 * transformations; SSE2 and AArch64 evaluate each operation in its own type. */
static_assert(FLT_EVAL_METHOD == 0, "myb needs FLT_EVAL_METHOD == 0 (no x87 evaluation)");
/* Bit-level code (myb_fma, hashing) reads a double through a uint64_t of the same byte
 * order. Clang has no target where they differ and does not define __FLOAT_WORD_ORDER__;
 * the check is kept for compilers that do. */
#if defined(__FLOAT_WORD_ORDER__)
static_assert(__FLOAT_WORD_ORDER__ == __BYTE_ORDER__,
              "myb needs the same byte order for double and uint64_t");
#endif

#endif /* MYB_SRC_CORE_FPGUARD_H_ */
