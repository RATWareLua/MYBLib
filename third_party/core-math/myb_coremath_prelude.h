/* SPDX-License-Identifier: BSD-3-Clause */
/* Included at the top of every vendored CORE-MATH source (the one line added to each file).
 *
 * 1. CORE-MATH assumes that __builtin_fma is correctly rounded; without hardware FMA the
 *    compiler turns it into a call to the C runtime's fma(), which on mingw-w64 is wrong for
 *    about 1 input in 6. The vendored files therefore call MYB_CM_FMA, mapped here to myb's
 *    exact integer implementation.
 * 2. The public entry points are renamed (cr_sin -> myb_cm_sin, ...) so that the library can
 *    be linked together with an upstream CORE-MATH without symbol clashes. */
#ifndef MYB_THIRD_PARTY_CORE_MATH_PRELUDE_H_
#define MYB_THIRD_PARTY_CORE_MATH_PRELUDE_H_

#include "cmath/fma_soft.h"

#define MYB_CM_FMA(a, b, c) myb_internal_fma((a), (b), (c))

#define cr_sin myb_cm_sin
#define cr_cos myb_cm_cos

#endif /* MYB_THIRD_PARTY_CORE_MATH_PRELUDE_H_ */
