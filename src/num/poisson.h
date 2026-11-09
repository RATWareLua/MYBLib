/* SPDX-License-Identifier: BSD-3-Clause */
/* num/poisson.h — Poisson series in fundamental arguments, the form of the conventional series
 * of Earth orientation (IERS Conventions 2010, eq. 5.16; tidal corrections to the EOP, ch. 5
 * and 8):
 *
 *   Q(t) = scale * ( sum_i P_i t^i + sum_terms t^j (a_s sin ARG + a_c cos ARG) ),
 *   ARG = sum_k n_k F_k(t),   F_k(t) = sum_i c_ki t^i (degree <= 4),
 *
 * coefficients read from NAIF text kernels. A group holds several series over the same
 * arguments (X, Y and s + XY/2 share most of theirs). At load time the distinct argument
 * vectors are collected; at each date the phasors exp(i m F_k) of every multiple m that occurs
 * are computed once (correctly rounded sine and cosine), the phasor exp(i ARG) of each distinct
 * argument is the product of those of its nonzero n_k in a fixed order, and every term using
 * it adds its contribution. The relative error of a term is a few units of roundoff, far below
 * the 0.1 µas cut-off of the series, and the result is reproducible bit for bit. Terms are
 * summed per series and power j, the values with compensated summation. */
#ifndef MYB_SRC_NUM_POISSON_H_
#define MYB_SRC_NUM_POISSON_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"
#include "myb/io.h"

enum {
  MYB_POISSON_MAX_ARGUMENTS = 16,
  MYB_POISSON_MAX_MULTIPLE = 64, /* |n_k| */
  MYB_POISSON_MAX_POWER = 8,     /* j and the degree of the polynomial part */
  MYB_POISSON_ARGUMENT_DEGREE = 4,
  MYB_POISSON_MAX_SERIES = 8
};

typedef struct myb_internal_poisson_group {
  /* Fundamental arguments: F_k = (sum_i c_ki t^i mod turn_k) * scale_k radians. */
  int arguments;
  double argument_polynomial[MYB_POISSON_MAX_ARGUMENTS][MYB_POISSON_ARGUMENT_DEGREE + 1];
  double argument_scale[MYB_POISSON_MAX_ARGUMENTS];
  double argument_turn[MYB_POISSON_MAX_ARGUMENTS];
  int max_multiple[MYB_POISSON_MAX_ARGUMENTS];
  /* Series. */
  int series;
  double scale[MYB_POISSON_MAX_SERIES];
  int degree[MYB_POISSON_MAX_SERIES];
  double polynomial[MYB_POISSON_MAX_SERIES][MYB_POISSON_MAX_POWER + 1];
  /* Distinct arguments: nonzero (k, n_k) pairs, and the terms that use each. */
  size_t unique;
  uint32_t* pair_start;    /* unique + 1 */
  uint8_t* pair_argument;
  int8_t* pair_multiple;
  uint32_t* entry_start;   /* unique + 1 */
  uint8_t* entry_series;
  uint8_t* entry_power;
  double* entry_amplitude; /* a_s, a_c */
} myb_internal_poisson_group;

/* Reads the arguments <prefix>_ARGUMENT_COUNT, _ARGUMENT_POLYNOMIAL (5 per argument),
 * _ARGUMENT_SCALE, _ARGUMENT_TURN, and for each name N of `names` the series N_SCALE,
 * N_POLYNOMIAL and N_TERMS (3 + ARGUMENT_COUNT numbers per term: j, a_s, a_c, n_1 ..).
 * Errors: MYB_E_NOT_FOUND, MYB_E_FORMAT (diag names the variable), MYB_E_NOMEM. */
[[nodiscard]] myb_status myb_internal_poisson_load(const myb_kpool* pool, const char* prefix,
                                                   const char* const* names, int count,
                                                   const myb_allocator* allocator,
                                                   myb_internal_poisson_group* group,
                                                   myb_diag* diag);

void myb_internal_poisson_free(myb_internal_poisson_group* group,
                               const myb_allocator* allocator);

/* angles[k] = F_k(t) in radians, reduced to one turn. */
void myb_internal_poisson_arguments(const myb_internal_poisson_group* group, double t,
                                    double* angles);

/* values[i] = series i at t and rates[i] = its derivative per unit of t, both multiplied by
 * the series' scale. */
void myb_internal_poisson_eval(const myb_internal_poisson_group* group, double t,
                               double* values, double* rates);

#endif /* MYB_SRC_NUM_POISSON_H_ */
