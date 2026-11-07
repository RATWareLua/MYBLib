/* SPDX-License-Identifier: BSD-3-Clause */
/* num/series.h — Poisson series  S(t) = sum A_i * t^k_i * sin(w_i * t + phi_i).
 *
 * The general form of the conventional series of astronomy (TDB-TT, nutation, precession of
 * the CIP, Earth orientation corrections): coefficients come from data, never from code.
 * Terms are summed in their given order with compensated summation; sin is correctly
 * rounded, so the value is reproducible bit for bit. */
#ifndef MYB_SRC_NUM_SERIES_H_
#define MYB_SRC_NUM_SERIES_H_

#include <stddef.h>

enum { MYB_SERIES_MAX_POWER = 8 };

typedef struct myb_internal_series {
  size_t count;
  const double* terms; /* 4 numbers per term: k, A, w, phi */
} myb_internal_series;

/* 1 if every k is an integer in [0, MYB_SERIES_MAX_POWER] and all numbers are finite. */
int myb_internal_series_valid(const double* terms, size_t count);

double myb_internal_series_eval(const myb_internal_series* series, double t);

#endif /* MYB_SRC_NUM_SERIES_H_ */
