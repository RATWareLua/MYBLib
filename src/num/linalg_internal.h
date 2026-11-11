/* SPDX-License-Identifier: BSD-3-Clause */
/* num/linalg_internal.h — the dense solvers of myb/num.h for systems of any order: the
 * public ones stop at 64, which the invariant circles of cr3bp/torus.c exceed. */
#ifndef MYB_SRC_NUM_LINALG_INTERNAL_H_
#define MYB_SRC_NUM_LINALG_INTERNAL_H_

#include <stddef.h>

#include "myb/core.h"

/* As myb_lu_solve, without the limit on n (the work is n^3 / 3). */
[[nodiscard]] myb_status myb_internal_lu_solve(size_t n, double* a, double* b, size_t nrhs);

#endif /* MYB_SRC_NUM_LINALG_INTERNAL_H_ */
