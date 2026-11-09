/* SPDX-License-Identifier: BSD-3-Clause */
/* Tables of dyn/wh_coefficients.h against MPFR (256 bits):
 *   - every SABA/SBAB coefficient is the double nearest to its closed form;
 *   - the kicks integrate polynomials exactly to the degree of their quadrature (Gauss with
 *     n nodes: 2n - 1, Lobatto with n + 1 nodes: 2n - 1), drifts and kicks each sum to 1,
 *     and the schemes are symmetric — so the first-order error is O(e h^(2n));
 *   - the corrector a_i are the doubles nearest to i sqrt(7/40), the b of every order are the
 *     doubles nearest to the solution of sum_i b_i a_i^k = (2^(2n) - 2) B_(2n) / (4 2^(2n) 2n),
 *     k = 2n - 1, and the stored doubles satisfy the system to round-off; for order 3,
 *     b_1 = 1/(96 alpha). */
#include <math.h> /* fabs */
#include <mpfr.h>
#include <stdio.h>

#include "core/fpguard.h"
#include "dyn/wh_coefficients.h"
#include "mtest.h"
#include "myb/dyn.h"

enum { P = 256 };

static double nearest(const mpfr_t x) {
  return mpfr_get_d(x, MPFR_RNDN);
}

/* x = p/q + s sqrt(r)/t */
static void closed(mpfr_t x, long p, long q, long s, long r, long t) {
  mpfr_t y;
  mpfr_init2(y, P);
  mpfr_set_si(x, p, MPFR_RNDN);
  mpfr_div_si(x, x, q, MPFR_RNDN);
  mpfr_set_si(y, r, MPFR_RNDN);
  mpfr_sqrt(y, y, MPFR_RNDN);
  mpfr_mul_si(y, y, s, MPFR_RNDN);
  mpfr_div_si(y, y, t, MPFR_RNDN);
  mpfr_add(x, x, y, MPFR_RNDN);
  mpfr_clear(y);
}

static void check_closed(double stored, long p, long q, long s, long r, long t) {
  mpfr_t x;
  mpfr_init2(x, P);
  closed(x, p, q, s, r, t);
  MT_CHECK_SAME(stored, nearest(x));
  mpfr_clear(x);
}

static void test_closed_forms(void) {
  const myb_wh_table* t = myb_wh_tables;
  /* SABA2 */
  check_closed(t[MYB_WH_SABA2].drift[0], 1, 2, -1, 3, 6);
  check_closed(t[MYB_WH_SABA2].drift[1], 0, 1, 1, 3, 3);
  /* SABA3 */
  check_closed(t[MYB_WH_SABA3].drift[0], 1, 2, -1, 15, 10);
  check_closed(t[MYB_WH_SABA3].drift[1], 0, 1, 1, 15, 10);
  check_closed(t[MYB_WH_SABA3].kick[0], 5, 18, 0, 1, 1);
  check_closed(t[MYB_WH_SABA3].kick[1], 4, 9, 0, 1, 1);
  /* SABA4: u, v = sqrt(525 +- 70 sqrt30) */
  mpfr_t u;
  mpfr_t v;
  mpfr_t w;
  mpfr_t x;
  mpfr_inits2(P, u, v, w, x, (mpfr_ptr)0);
  mpfr_set_ui(w, 30, MPFR_RNDN);
  mpfr_sqrt(w, w, MPFR_RNDN);
  mpfr_mul_ui(u, w, 70, MPFR_RNDN);
  mpfr_ui_sub(v, 525, u, MPFR_RNDN);
  mpfr_add_ui(u, u, 525, MPFR_RNDN);
  mpfr_sqrt(u, u, MPFR_RNDN);
  mpfr_sqrt(v, v, MPFR_RNDN);
  mpfr_div_ui(x, u, 70, MPFR_RNDN);
  mpfr_ui_sub(x, 1, x, MPFR_RNDN);
  mpfr_sub_d(x, x, 0.5, MPFR_RNDN); /* 1/2 - u/70 */
  MT_CHECK_SAME(t[MYB_WH_SABA4].drift[0], nearest(x));
  mpfr_sub(x, u, v, MPFR_RNDN);
  mpfr_div_ui(x, x, 70, MPFR_RNDN);
  MT_CHECK_SAME(t[MYB_WH_SABA4].drift[1], nearest(x));
  mpfr_div_ui(x, v, 35, MPFR_RNDN);
  MT_CHECK_SAME(t[MYB_WH_SABA4].drift[2], nearest(x));
  mpfr_clears(u, v, w, x, (mpfr_ptr)0);
  check_closed(t[MYB_WH_SABA4].kick[0], 1, 4, -1, 30, 72);
  check_closed(t[MYB_WH_SABA4].kick[1], 1, 4, 1, 30, 72);
  /* SBAB */
  check_closed(t[MYB_WH_SBAB2].kick[0], 1, 6, 0, 1, 1);
  check_closed(t[MYB_WH_SBAB2].kick[1], 2, 3, 0, 1, 1);
  check_closed(t[MYB_WH_SBAB3].drift[0], 1, 2, -1, 5, 10);
  check_closed(t[MYB_WH_SBAB3].drift[1], 0, 1, 1, 5, 5);
  check_closed(t[MYB_WH_SBAB3].kick[0], 1, 12, 0, 1, 1);
  check_closed(t[MYB_WH_SBAB3].kick[1], 5, 12, 0, 1, 1);
  check_closed(t[MYB_WH_SBAB4].drift[0], 1, 2, -1, 21, 14);
  check_closed(t[MYB_WH_SBAB4].drift[1], 0, 1, 1, 21, 14);
  check_closed(t[MYB_WH_SBAB4].kick[0], 1, 20, 0, 1, 1);
  check_closed(t[MYB_WH_SBAB4].kick[1], 49, 180, 0, 1, 1);
  check_closed(t[MYB_WH_SBAB4].kick[2], 16, 45, 0, 1, 1);
}

/* Quadrature exactness and symmetry of every scheme, with the stored doubles (converted
 * exactly) in MPFR: |sum_k d_k s_k^p - 1/(p+1)| at round-off of the coefficients. */
static void test_quadrature(void) {
  for (int scheme = 0; scheme < MYB_WH_SCHEMES; ++scheme) {
    const myb_wh_table* t = &myb_wh_tables[scheme];
    const int drifts = t->kick_first ? t->kicks - 1 : t->kicks + 1;
    double drift_sum = 0.0;
    double kick_sum = 0.0;
    for (int k = 0; k < drifts; ++k) {
      drift_sum += t->drift[k];
      MT_CHECK_SAME(t->drift[k], t->drift[drifts - 1 - k]);
    }
    for (int k = 0; k < t->kicks; ++k) {
      kick_sum += t->kick[k];
      MT_CHECK_SAME(t->kick[k], t->kick[t->kicks - 1 - k]);
    }
    MT_CHECK(fabs(drift_sum - 1.0) <= 0x1p-52 && fabs(kick_sum - 1.0) <= 0x1p-52);
    /* node positions: cumulative drifts before each kick */
    mpfr_t s[MYB_WH_MAX_KICKS];
    mpfr_t term;
    mpfr_t sum;
    mpfr_inits2(P, term, sum, (mpfr_ptr)0);
    mpfr_t position;
    mpfr_init2(position, P);
    mpfr_set_ui(position, 0, MPFR_RNDN);
    for (int k = 0; k < t->kicks; ++k) {
      mpfr_init2(s[k], P);
      const int before = t->kick_first ? k : k + 1; /* drifts applied before kick k */
      mpfr_set_ui(position, 0, MPFR_RNDN);
      for (int d = 0; d < before; ++d) {
        mpfr_add_d(position, position, t->drift[d], MPFR_RNDN);
      }
      mpfr_set(s[k], position, MPFR_RNDN);
    }
    const int n = t->kick_first ? t->kicks - 1 : t->kicks; /* SABAn / SBABn */
    double worst = 0.0;
    for (int p = 0; p <= 2 * n - 1; ++p) {
      mpfr_set_ui(sum, 0, MPFR_RNDN);
      for (int k = 0; k < t->kicks; ++k) {
        mpfr_pow_ui(term, s[k], (unsigned long)p, MPFR_RNDN);
        mpfr_mul_d(term, term, t->kick[k], MPFR_RNDN);
        mpfr_add(sum, sum, term, MPFR_RNDN);
      }
      mpfr_set_ui(term, 1, MPFR_RNDN);
      mpfr_div_ui(term, term, (unsigned long)(p + 1), MPFR_RNDN);
      mpfr_sub(sum, sum, term, MPFR_RNDN);
      worst = fmax(worst, fabs(mpfr_get_d(sum, MPFR_RNDN)));
    }
    /* the next degree is not integrated exactly: the order is not higher than claimed */
    mpfr_set_ui(sum, 0, MPFR_RNDN);
    for (int k = 0; k < t->kicks; ++k) {
      mpfr_pow_ui(term, s[k], (unsigned long)(2 * n), MPFR_RNDN);
      mpfr_mul_d(term, term, t->kick[k], MPFR_RNDN);
      mpfr_add(sum, sum, term, MPFR_RNDN);
    }
    mpfr_set_ui(term, 1, MPFR_RNDN);
    mpfr_div_ui(term, term, (unsigned long)(2 * n + 1), MPFR_RNDN);
    mpfr_sub(sum, sum, term, MPFR_RNDN);
    const double next = fabs(mpfr_get_d(sum, MPFR_RNDN));
    printf("scheme %d (%s%d): quadrature exact to degree %d within %.1e, degree %d off by "
           "%.1e\n",
           scheme, t->kick_first ? "SBAB" : "SABA", n, 2 * n - 1, worst, 2 * n, next);
    MT_CHECK(worst < 4e-16);
    MT_CHECK(next > 1e-6);
    for (int k = 0; k < t->kicks; ++k) {
      mpfr_clear(s[k]);
    }
    mpfr_clears(term, sum, position, (mpfr_ptr)0);
  }
}

static void test_correctors(void) {
  static const long bn[8] = {1, -1, 1, -1, 5, -691, 7, -3617};
  static const long bd[8] = {6, 30, 42, 30, 66, 2730, 6, 510};
  mpfr_t alpha;
  mpfr_t t1;
  mpfr_t t2;
  mpfr_inits2(P, alpha, t1, t2, (mpfr_ptr)0);
  mpfr_set_ui(alpha, 7, MPFR_RNDN);
  mpfr_div_ui(alpha, alpha, 40, MPFR_RNDN);
  mpfr_sqrt(alpha, alpha, MPFR_RNDN);
  for (int i = 0; i < MYB_WH_CORRECTOR_STAGES; ++i) {
    mpfr_mul_ui(t1, alpha, (unsigned long)(i + 1), MPFR_RNDN);
    MT_CHECK_SAME(myb_wh_corrector_a[i], nearest(t1));
  }
  for (int m = 1; m <= MYB_WH_CORRECTOR_STAGES; ++m) {
    mpfr_t a[MYB_WH_CORRECTOR_STAGES][MYB_WH_CORRECTOR_STAGES + 1];
    for (int r = 0; r < m; ++r) {
      for (int c = 0; c <= m; ++c) {
        mpfr_init2(a[r][c], P);
      }
    }
    for (int r = 0; r < m; ++r) {
      const int n = r + 1;
      for (int c = 0; c < m; ++c) {
        mpfr_mul_ui(t1, alpha, (unsigned long)(c + 1), MPFR_RNDN);
        mpfr_pow_ui(a[r][c], t1, (unsigned long)(2 * n - 1), MPFR_RNDN);
      }
      mpfr_set_si(t1, bn[n - 1], MPFR_RNDN);
      mpfr_div_si(t1, t1, bd[n - 1], MPFR_RNDN);
      mpfr_set_ui(t2, 1, MPFR_RNDN);
      mpfr_mul_2ui(t2, t2, (unsigned long)(2 * n), MPFR_RNDN);
      mpfr_sub_ui(t2, t2, 2, MPFR_RNDN);
      mpfr_mul(t1, t1, t2, MPFR_RNDN);
      mpfr_div_2ui(t1, t1, (unsigned long)(2 * n + 2), MPFR_RNDN);
      mpfr_div_ui(a[r][m], t1, (unsigned long)(2 * n), MPFR_RNDN);
    }
    /* residual of the stored doubles, relative to the right-hand side */
    double worst = 0.0;
    for (int r = 0; r < m; ++r) {
      mpfr_set_ui(t1, 0, MPFR_RNDN);
      for (int c = 0; c < m; ++c) {
        mpfr_mul_d(t2, a[r][c], myb_wh_corrector_b[m - 1][c], MPFR_RNDN);
        mpfr_add(t1, t1, t2, MPFR_RNDN);
      }
      mpfr_sub(t1, t1, a[r][m], MPFR_RNDN);
      mpfr_div(t1, t1, a[r][m], MPFR_RNDN);
      worst = fmax(worst, fabs(mpfr_get_d(t1, MPFR_RNDN)));
    }
    /* Gauss-Jordan elimination for the exact solution */
    for (int p = 0; p < m; ++p) {
      int best = p;
      for (int r = p + 1; r < m; ++r) {
        if (mpfr_cmpabs(a[r][p], a[best][p]) > 0) {
          best = r;
        }
      }
      for (int c = 0; c <= m; ++c) {
        mpfr_swap(a[p][c], a[best][c]);
      }
      for (int r = 0; r < m; ++r) {
        if (r == p) {
          continue;
        }
        mpfr_div(t1, a[r][p], a[p][p], MPFR_RNDN);
        for (int c = p; c <= m; ++c) {
          mpfr_mul(t2, t1, a[p][c], MPFR_RNDN);
          mpfr_sub(a[r][c], a[r][c], t2, MPFR_RNDN);
        }
      }
    }
    for (int r = 0; r < m; ++r) {
      mpfr_div(t1, a[r][m], a[r][r], MPFR_RNDN);
      MT_CHECK_SAME(myb_wh_corrector_b[m - 1][r], nearest(t1));
    }
    printf("corrector of order %2d: system residual of the stored b %.1e (relative)\n", 2 * m + 1,
           worst);
    MT_CHECK(worst < 1e-14 * (double)(1 << m));
    for (int r = 0; r < m; ++r) {
      for (int c = 0; c <= m; ++c) {
        mpfr_clear(a[r][c]);
      }
    }
  }
  /* order 3: b_1 = 1/(96 alpha) = 0.024900596027799867..., Wisdom, Holman & Touma 1996 */
  mpfr_mul_ui(t1, alpha, 96, MPFR_RNDN);
  mpfr_ui_div(t1, 1, t1, MPFR_RNDN);
  MT_CHECK_SAME(myb_wh_corrector_b[0][0], nearest(t1));
  mpfr_clears(alpha, t1, t2, (mpfr_ptr)0);
}

int main(void) {
  test_closed_forms();
  test_quadrature();
  test_correctors();
  mpfr_free_cache();
  return mt_finish("test_mpfr_wh");
}
