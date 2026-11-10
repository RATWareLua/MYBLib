/* SPDX-License-Identifier: BSD-3-Clause */
/* num/erf.h — error function and its complement (fdlibm algorithm, ~1 ulp, identical bits
 * on every platform; libm's erf is not used: it differs between C runtimes). */
#ifndef MYB_SRC_NUM_ERF_H_
#define MYB_SRC_NUM_ERF_H_

double myb_internal_erf(double x);
double myb_internal_erfc(double x);

#endif /* MYB_SRC_NUM_ERF_H_ */
