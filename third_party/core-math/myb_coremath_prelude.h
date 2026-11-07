/* SPDX-License-Identifier: BSD-3-Clause */
/* Included at the top of every vendored CORE-MATH source (the line added by import.cmake).
 *
 * 1. CORE-MATH assumes that __builtin_fma is correctly rounded; without hardware FMA the
 *    compiler turns it into a call to the C runtime's fma(), which on mingw-w64 is wrong for
 *    about 1 input in 6. The vendored files call MYB_CM_FMA instead: myb's exact fma.
 * 2. Debug output (MYB_CM_PRINTF) is dropped, and the "unexpected worst case" path of atan2,
 *    proved unreachable by the CORE-MATH authors, terminates through myb's fail-fast check
 *    instead of printing to stdout and calling exit(1) from inside a library.
 * 3. The public entry points are renamed (cr_sin -> myb_cm_sin, ...) so that the library can
 *    be linked together with an upstream CORE-MATH without symbol clashes. */
#ifndef MYB_THIRD_PARTY_CORE_MATH_PRELUDE_H_
#define MYB_THIRD_PARTY_CORE_MATH_PRELUDE_H_

#include "cmath/fma_soft.h"
#include "core/check.h"

#define MYB_CM_FMA(a, b, c) myb_internal_fma((a), (b), (c))
#define MYB_CM_PRINTF(...) ((void)0)
#define MYB_CM_EXIT(code)                                 \
  myb_internal_check_failed(__FILE__, __LINE__, __func__, \
                            "CORE-MATH: unexpected worst case (proved unreachable)")

#define cr_sin myb_cm_sin
#define cr_cos myb_cm_cos
#define cr_tan myb_cm_tan
#define cr_atan myb_cm_atan
#define cr_atan2 myb_cm_atan2
#define cr_asin myb_cm_asin
#define cr_acos myb_cm_acos
#define cr_exp myb_cm_exp
#define cr_expm1 myb_cm_expm1
#define cr_log myb_cm_log
#define cr_log1p myb_cm_log1p
#define cr_sinh myb_cm_sinh
#define cr_cosh myb_cm_cosh
#define cr_tanh myb_cm_tanh
#define cr_asinh myb_cm_asinh
#define cr_acosh myb_cm_acosh
#define cr_atanh myb_cm_atanh
#define cr_cbrt myb_cm_cbrt
#define cr_hypot myb_cm_hypot

#endif /* MYB_THIRD_PARTY_CORE_MATH_PRELUDE_H_ */
