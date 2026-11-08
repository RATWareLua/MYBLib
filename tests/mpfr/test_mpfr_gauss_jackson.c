/* SPDX-License-Identifier: BSD-3-Clause */
/* The Gauss-Jackson tables (src/integ/gj_coefficients.h) recomputed in MPFR at 256 bits and
 * compared bit for bit: Bernoulli numbers from their recurrence, the Taylor coefficients of
 *   P1(D) = -sum_{m>=1} B_m D^(m-1) / m!,   P2(D) = sum_{m>=2} B_m (m-1) D^(m-2) / m!,
 * the Lagrange basis of the nodes -4..4 in monomial form, and
 *   V[j][k] = sum_q p1_q L_k^(q)(j),   X[j][k] = sum_q p2_q L_k^(q)(j). */
#include <mpfr.h>

#include "core/fpguard.h"
#include "integ/gj_coefficients.h"
#include "mtest.h"

enum { PREC = 256, NODES = 9, TERMS = 9 };

static double binomial(int n, int k) {
  double r = 1.0;
  for (int i = 0; i < k; ++i) {
    r = r * (double)(n - i) / (double)(i + 1); /* exact for these small values */
  }
  return r;
}

/* The double nearest to an MPFR value; values that are exactly zero as rationals come out of
 * the inexact divisions as residues below 2^-200 and are taken as +0 (the table's zero). */
static double nearest(mpfr_srcptr x) {
  const double d = mpfr_get_d(x, MPFR_RNDN);
  return (d < 0x1p-200 && d > -0x1p-200) ? 0.0 : d;
}

static double factorial(int n) {
  double r = 1.0;
  for (int i = 2; i <= n; ++i) {
    r *= (double)i;
  }
  return r;
}

int main(void) {
  mpfr_t bern[12];
  mpfr_t t;
  mpfr_init2(t, PREC);
  for (int m = 0; m < 12; ++m) {
    mpfr_init2(bern[m], PREC);
    /* sum_{k=0}^{m} C(m+1, k) B_k = 0 */
    if (m == 0) {
      mpfr_set_ui(bern[0], 1, MPFR_RNDN);
      continue;
    }
    mpfr_set_zero(bern[m], 1);
    for (int k = 0; k < m; ++k) {
      mpfr_mul_d(t, bern[k], binomial(m + 1, k), MPFR_RNDN);
      mpfr_add(bern[m], bern[m], t, MPFR_RNDN);
    }
    mpfr_div_si(bern[m], bern[m], -(m + 1), MPFR_RNDN);
  }
  /* Taylor coefficients of P1 and P2 */
  mpfr_t p1[TERMS];
  mpfr_t p2[TERMS];
  for (int q = 0; q < TERMS; ++q) {
    mpfr_inits2(PREC, p1[q], p2[q], (mpfr_ptr)0);
    mpfr_div_d(p1[q], bern[q + 1], -factorial(q + 1), MPFR_RNDN);
    mpfr_mul_ui(p2[q], bern[q + 2], (unsigned long)(q + 1), MPFR_RNDN);
    mpfr_div_d(p2[q], p2[q], factorial(q + 2), MPFR_RNDN);
  }
  /* Lagrange basis in monomial form */
  mpfr_t lag[NODES][NODES];
  for (int k = 0; k < NODES; ++k) {
    for (int m = 0; m < NODES; ++m) {
      mpfr_init2(lag[k][m], PREC);
      mpfr_set_zero(lag[k][m], 1);
    }
    mpfr_set_ui(lag[k][0], 1, MPFR_RNDN);
    int degree = 0;
    for (int i = 0; i < NODES; ++i) {
      if (i == k) {
        continue;
      }
      /* multiply by (u - node_i) / (node_k - node_i) */
      const double node_i = (double)(i - 4);
      const double denominator = (double)(k - i);
      for (int m = degree + 1; m >= 0; --m) {
        mpfr_set_zero(t, 1);
        if (m >= 1) {
          mpfr_set(t, lag[k][m - 1], MPFR_RNDN);
        }
        if (m <= degree) {
          mpfr_t prod;
          mpfr_init2(prod, PREC);
          mpfr_mul_d(prod, lag[k][m], -node_i, MPFR_RNDN);
          mpfr_add(t, t, prod, MPFR_RNDN);
          mpfr_clear(prod);
        }
        mpfr_div_d(lag[k][m], t, denominator, MPFR_RNDN);
      }
      ++degree;
    }
    for (int m = 0; m < NODES; ++m) {
      MT_CHECK_SAME(myb_gj_lagrange[k][m], nearest(lag[k][m]));
    }
  }
  /* V and X: sum_q p_q L_k^(q)(j) for j = -4..5 */
  mpfr_t dv;
  mpfr_t dx;
  mpfr_t deriv;
  mpfr_t power;
  mpfr_inits2(PREC, dv, dx, deriv, power, (mpfr_ptr)0);
  for (int row = 0; row < 10; ++row) {
    const double j = (double)(row - 4);
    for (int k = 0; k < NODES; ++k) {
      mpfr_set_zero(dv, 1);
      mpfr_set_zero(dx, 1);
      for (int q = 0; q < TERMS; ++q) {
        /* L_k^(q)(j) = sum_{m>=q} lag[k][m] m!/(m-q)! j^(m-q) */
        mpfr_set_zero(deriv, 1);
        for (int m = q; m < NODES; ++m) {
          mpfr_set_d(power, 1.0, MPFR_RNDN);
          for (int e = 0; e < m - q; ++e) {
            mpfr_mul_d(power, power, j, MPFR_RNDN);
          }
          mpfr_mul(power, power, lag[k][m], MPFR_RNDN);
          mpfr_mul_d(power, power, factorial(m) / factorial(m - q), MPFR_RNDN);
          mpfr_add(deriv, deriv, power, MPFR_RNDN);
        }
        mpfr_mul(t, deriv, p1[q], MPFR_RNDN);
        mpfr_add(dv, dv, t, MPFR_RNDN);
        mpfr_mul(t, deriv, p2[q], MPFR_RNDN);
        mpfr_add(dx, dx, t, MPFR_RNDN);
      }
      MT_CHECK_SAME(myb_gj_velocity[row][k], nearest(dv));
      MT_CHECK_SAME(myb_gj_position[row][k], nearest(dx));
    }
  }
  mpfr_clears(dv, dx, deriv, power, t, (mpfr_ptr)0);
  for (int q = 0; q < TERMS; ++q) {
    mpfr_clears(p1[q], p2[q], (mpfr_ptr)0);
  }
  for (int m = 0; m < 12; ++m) {
    mpfr_clear(bern[m]);
  }
  for (int k = 0; k < NODES; ++k) {
    for (int m = 0; m < NODES; ++m) {
      mpfr_clear(lag[k][m]);
    }
  }
  mpfr_free_cache();
  return mt_finish("test_mpfr_gauss_jackson");
}
