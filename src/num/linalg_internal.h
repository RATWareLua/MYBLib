/* SPDX-License-Identifier: BSD-3-Clause */
/* num/linalg_internal.h — the dense solvers of myb/num.h for systems of any order: the
 * public ones stop at 64, which the invariant circles of cr3bp/torus.c (and the spectra of
 * their linearised maps) exceed. */
#ifndef MYB_SRC_NUM_LINALG_INTERNAL_H_
#define MYB_SRC_NUM_LINALG_INTERNAL_H_

#include <stddef.h>

#include "myb/core.h"

/* As myb_lu_solve, without the limit on n (the work is n^3 / 3). */
[[nodiscard]] myb_status myb_internal_lu_solve(size_t n, double* a, double* b, size_t nrhs);

/* As myb_eigenvalues, without the limit on n; work: n doubles (the reflectors of the
 * Hessenberg reduction). The same operations in the same order: for n <= 64 the results are
 * those of myb_eigenvalues bit for bit. */
[[nodiscard]] myb_status myb_internal_eigenvalues(size_t n, double* a, double* re, double* im,
                                                  double* work);

/* As myb_eigenvector_real, without the limit on n; work: n x n doubles. */
[[nodiscard]] myb_status myb_internal_eigenvector_real(size_t n, const double* a, double lambda,
                                                       double* v, double* work);

#endif /* MYB_SRC_NUM_LINALG_INTERNAL_H_ */
