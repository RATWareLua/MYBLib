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

/* An eigenvector v = vr + i vi of the real matrix a (n x n, row-major) for the eigenvalue
 * re + i im, by inverse iteration in complex arithmetic: (A - (lambda + delta) I) x = v three
 * times from a fixed start, delta = 1e-10 (|lambda| + max |a_ij|) keeping the system just
 * regular, one LU with partial pivoting (by |re| + |im|) for all three. v of unit Euclidean
 * norm, its largest component real and positive. work: 2 n^2 + n doubles. Errors: MYB_E_DOMAIN
 * (values not finite, a singular system), MYB_E_NOCONV (no vector came out). */
[[nodiscard]] myb_status myb_internal_eigenvector_complex(size_t n, const double* a, double re,
                                                          double im, double* vr, double* vi,
                                                          double* work);

#endif /* MYB_SRC_NUM_LINALG_INTERNAL_H_ */
