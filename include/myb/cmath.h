/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/cmath.h — correctly rounded elementary functions.
 *
 * Every function here returns RN(f(x)), the exactly rounded value, so its result is the
 * same bit pattern on every platform and every compiler. The system libm gives no such
 * guarantee (UCRT, glibc, musl and mingw-w64 differ in the last bits), and myb never calls
 * it for anything but exact operations.
 *
 * Implementations:
 *   myb_fma         own exact implementation (128-bit integer arithmetic). mingw-w64's fma()
 *                   is NOT correctly rounded (about 1 result in 6 is wrong on random input),
 *                   which is why myb does not use the CRT here.
 *   myb_sin/cos     CORE-MATH (Zimmermann, Sibidanov et al., MIT licence), vendored in
 *                   third_party/core-math, with its fma routed to myb_fma.
 *
 * Special values follow IEEE-754 / C Annex F (sin(+-0) = +-0, sin(inf) = NaN, ...).
 * Floating-point exception flags and errno are not part of the contract. */
#ifndef MYB_CMATH_H_
#define MYB_CMATH_H_

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

MYB_API double myb_fma(double a, double b, double c);
MYB_API double myb_sin(double x);
MYB_API double myb_cos(double x);

#ifdef __cplusplus
}
#endif

#endif /* MYB_CMATH_H_ */
