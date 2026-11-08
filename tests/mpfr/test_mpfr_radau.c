/* SPDX-License-Identifier: BSD-3-Clause */
/* The Gauss-Radau constants of IAS15 (src/integ/ias15_coefficients.h) recomputed in MPFR at
 * 512 bits and compared bit for bit with the tables: every entry must be the double nearest to
 * the exact value. Nodes: roots of P_7(2h - 1) + P_8(2h - 1), refined by bisection inside a
 * bracket around the tabulated value (which also proves there is a root next to it). */
#include <mpfr.h>

#include "core/fpguard.h"
#include "integ/ias15_coefficients.h"
#include "mtest.h"

enum { PREC = 512 };

/* P_7(x) + P_8(x) at x = 2h - 1 (three-term recurrence). */
static void radau_function(mpfr_srcptr h, mpfr_ptr out) {
  mpfr_t x;
  mpfr_t p0;
  mpfr_t p1;
  mpfr_t p2;
  mpfr_t t;
  mpfr_inits2(PREC, x, p0, p1, p2, t, (mpfr_ptr)0);
  mpfr_mul_ui(x, h, 2, MPFR_RNDN);
  mpfr_sub_ui(x, x, 1, MPFR_RNDN);
  mpfr_set_ui(p0, 1, MPFR_RNDN);
  mpfr_set(p1, x, MPFR_RNDN);
  mpfr_t p7;
  mpfr_init2(p7, PREC);
  for (unsigned k = 1; k < 8; ++k) { /* p2 = ((2k + 1) x p1 - k p0) / (k + 1) */
    mpfr_mul(t, x, p1, MPFR_RNDN);
    mpfr_mul_ui(t, t, 2 * k + 1, MPFR_RNDN);
    mpfr_mul_ui(p2, p0, k, MPFR_RNDN);
    mpfr_sub(p2, t, p2, MPFR_RNDN);
    mpfr_div_ui(p2, p2, k + 1, MPFR_RNDN);
    mpfr_set(p0, p1, MPFR_RNDN);
    mpfr_set(p1, p2, MPFR_RNDN);
    if (k == 6) {
      mpfr_set(p7, p1, MPFR_RNDN);
    }
  }
  mpfr_add(out, p7, p1, MPFR_RNDN);
  mpfr_clears(x, p0, p1, p2, t, p7, (mpfr_ptr)0);
}

/* Root inside [h - 1e-12, h + 1e-12]; 0 if the function does not change sign there. */
static int refine_node(double h, mpfr_ptr root) {
  mpfr_t a;
  mpfr_t b;
  mpfr_t m;
  mpfr_t fa;
  mpfr_t fm;
  mpfr_inits2(PREC, a, b, m, fa, fm, (mpfr_ptr)0);
  mpfr_set_d(a, h - 1e-12, MPFR_RNDN);
  mpfr_set_d(b, h + 1e-12, MPFR_RNDN);
  radau_function(a, fa);
  radau_function(b, fm);
  int ok = mpfr_sgn(fa) * mpfr_sgn(fm) < 0;
  /* bound: 600 halvings of a 2e-12 interval reach far below 2^-512 */
  for (int i = 0; ok && i < 600; ++i) {
    mpfr_add(m, a, b, MPFR_RNDN);
    mpfr_div_ui(m, m, 2, MPFR_RNDN);
    radau_function(m, fm);
    if (mpfr_sgn(fm) == 0) {
      mpfr_set(a, m, MPFR_RNDN);
      break;
    }
    if (mpfr_sgn(fm) == mpfr_sgn(fa)) {
      mpfr_set(a, m, MPFR_RNDN);
      mpfr_set(fa, fm, MPFR_RNDN);
    } else {
      mpfr_set(b, m, MPFR_RNDN);
    }
  }
  mpfr_set(root, a, MPFR_RNDN);
  mpfr_clears(a, b, m, fa, fm, (mpfr_ptr)0);
  return ok;
}

int main(void) {
  mpfr_t h[8];
  for (int i = 0; i < 8; ++i) {
    mpfr_init2(h[i], PREC);
  }
  mpfr_set_zero(h[0], 1);
  MT_CHECK(myb_ias15_h[0] == 0.0);
  for (int i = 1; i < 8; ++i) {
    MT_CHECK(refine_node(myb_ias15_h[i], h[i]));
    MT_CHECK_SAME(myb_ias15_h[i], mpfr_get_d(h[i], MPFR_RNDN));
  }
  mpfr_t t;
  mpfr_init2(t, PREC);
  /* C[k][j]: coefficient of h^(k+1) in h (h - h_1) ... (h - h_j). */
  mpfr_t poly[9];
  mpfr_t c[7][7];
  for (int i = 0; i < 9; ++i) {
    mpfr_init2(poly[i], PREC);
    mpfr_set_zero(poly[i], 1);
  }
  mpfr_set_ui(poly[0], 1, MPFR_RNDN);
  for (int j = 0; j < 7; ++j) {
    for (int d = 8; d >= 1; --d) { /* multiply by (h - h_j) */
      mpfr_mul(t, poly[d], h[j], MPFR_RNDN);
      mpfr_sub(poly[d], poly[d - 1], t, MPFR_RNDN);
    }
    mpfr_mul(t, poly[0], h[j], MPFR_RNDN);
    mpfr_neg(poly[0], t, MPFR_RNDN);
    for (int k = 0; k < 7; ++k) {
      mpfr_init2(c[k][j], PREC);
      mpfr_set(c[k][j], poly[k + 1], MPFR_RNDN);
      MT_CHECK_SAME(myb_ias15_c[k][j], mpfr_get_d(c[k][j], MPFR_RNDN));
    }
  }
  /* D = C^-1 by back substitution (C is unit upper triangular). */
  mpfr_t d[7][7];
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      mpfr_init2(d[i][j], PREC);
    }
  }
  for (int col = 0; col < 7; ++col) {
    for (int row = 6; row >= 0; --row) {
      mpfr_set_ui(d[row][col], row == col ? 1 : 0, MPFR_RNDN);
      for (int k = row + 1; k < 7; ++k) {
        mpfr_mul(t, c[row][k], d[k][col], MPFR_RNDN);
        mpfr_sub(d[row][col], d[row][col], t, MPFR_RNDN);
      }
      mpfr_div(d[row][col], d[row][col], c[row][row], MPFR_RNDN);
    }
  }
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      MT_CHECK_SAME(myb_ias15_d[i][j], mpfr_get_d(d[i][j], MPFR_RNDN));
    }
  }
  /* R[n][k] = 1 / (h_n - h_k) for k < n, 0 otherwise. */
  for (int n = 0; n < 8; ++n) {
    for (int k = 0; k < 8; ++k) {
      double expected = 0.0;
      if (k < n) {
        mpfr_sub(t, h[n], h[k], MPFR_RNDN);
        mpfr_ui_div(t, 1, t, MPFR_RNDN);
        expected = mpfr_get_d(t, MPFR_RNDN);
      }
      MT_CHECK_SAME(myb_ias15_r[n][k], expected);
    }
  }
  for (int i = 0; i < 9; ++i) {
    mpfr_clear(poly[i]);
  }
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      mpfr_clear(c[i][j]);
      mpfr_clear(d[i][j]);
    }
  }
  for (int i = 0; i < 8; ++i) {
    mpfr_clear(h[i]);
  }
  mpfr_clear(t);
  mpfr_free_cache();
  return mt_finish("test_mpfr_radau");
}
