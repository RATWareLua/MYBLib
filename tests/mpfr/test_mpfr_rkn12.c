/* SPDX-License-Identifier: BSD-3-Clause */
/* The RKN12(10) tables (src/integ/rkn12_coefficients.h) against the exact fractions of
 * tests/data/rkn12/dprkn12.txt (GMP rationals): every double must be the nearest to its
 * fraction, and the order conditions are re-checked exactly:
 *   sum_j a_ij = c_i^2 / 2 (rows 14-17: to 1e-45, their entries are 50-digit values),
 *   sum_i b_i c_i^k = 1/((k+1)(k+2)) for k <= 10, sum_i bp_i c_i^k = 1/(k+1) for k <= 11,
 *   b_i = bp_i (1 - c_i), and the error weights annihilate c^k for k <= 8. */
#include <math.h> /* fabs */
#include <mpfr.h>
#include <stdio.h>
#include <string.h>

#include "core/fpguard.h"
#include "integ/rkn12_coefficients.h"
#include "mtest.h"

enum { S = MYB_RKN12_STAGES, LINE = 512 };

static int read_section(FILE* file, const char* name, mpq_t* values, int count) {
  char line[LINE];
  if (fgets(line, LINE, file) == NULL || strncmp(line, name, strlen(name)) != 0) {
    return 0;
  }
  for (int i = 0; i < count; ++i) {
    char num[LINE];
    char den[LINE];
    if (fgets(line, LINE, file) == NULL || sscanf(line, "%511s %511s", num, den) != 2) {
      return 0;
    }
    char text[2 * LINE + 2];
    snprintf(text, sizeof text, "%s/%s", num, den);
    if (mpq_set_str(values[i], text, 10) != 0) {
      return 0;
    }
    mpq_canonicalize(values[i]);
  }
  return 1;
}

/* The double nearest to q: mpfr_set_q rounds the fraction correctly, once, to 53 bits. */
static double nearest(const mpq_t q) {
  mpfr_t y;
  mpfr_init2(y, 53);
  mpfr_set_q(y, q, MPFR_RNDN);
  const double d = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear(y);
  return d;
}

int main(void) {
  static mpq_t c[S];
  static mpq_t a[S * S];
  static mpq_t b[S];
  static mpq_t bp[S];
  static mpq_t bt[S];
  static mpq_t bpt[S];
  mpq_t* all[6] = {c, a, b, bp, bt, bpt};
  const int sizes[6] = {S, S * S, S, S, S, S};
  for (int k = 0; k < 6; ++k) {
    for (int i = 0; i < sizes[k]; ++i) {
      mpq_init(all[k][i]);
    }
  }
  FILE* file = fopen(MYB_TEST_DATA_DIR "/rkn12/dprkn12.txt", "r");
  MT_CHECK(file != NULL);
  if (file == NULL) {
    return mt_finish("test_mpfr_rkn12");
  }
  const char* names[6] = {"c", "a", "b", "bp", "bt", "bpt"};
  for (int k = 0; k < 6; ++k) {
    MT_CHECK(read_section(file, names[k], all[k], sizes[k]));
  }
  fclose(file);
  /* Bit for bit. */
  for (int i = 0; i < S; ++i) {
    MT_CHECK_SAME(myb_rkn12_c[i], nearest(c[i]));
    MT_CHECK_SAME(myb_rkn12_b[i], nearest(b[i]));
    MT_CHECK_SAME(myb_rkn12_bp[i], nearest(bp[i]));
    MT_CHECK_SAME(myb_rkn12_bt[i], nearest(bt[i]));
    MT_CHECK_SAME(myb_rkn12_bpt[i], nearest(bpt[i]));
    for (int j = 0; j < S; ++j) {
      MT_CHECK_SAME(myb_rkn12_a[i][j], nearest(a[i * S + j]));
    }
  }
  /* Order conditions, exactly. */
  mpq_t sum;
  mpq_t term;
  mpq_t power;
  mpq_t target;
  mpq_inits(sum, term, power, target, (mpq_ptr)0);
  int rows_exact = 0;
  for (int i = 0; i < S; ++i) {
    mpq_set_ui(sum, 0, 1);
    for (int j = 0; j < S; ++j) {
      mpq_add(sum, sum, a[i * S + j]);
    }
    mpq_mul(target, c[i], c[i]);
    mpq_div_2exp(target, target, 1);
    if (mpq_equal(sum, target)) {
      ++rows_exact;
    } else {
      mpq_sub(term, sum, target);
      MT_CHECK(i >= 13 && fabs(mpq_get_d(term)) < 1e-45);
    }
    mpq_set_ui(term, 1, 1);
    mpq_sub(term, term, c[i]);
    mpq_mul(term, term, bp[i]);
    MT_CHECK(mpq_equal(term, b[i]));
  }
  MT_CHECK(rows_exact == 13);
  for (int k = 0; k <= 11; ++k) {
    mpq_t qb;
    mpq_t qbp;
    mpq_t qbt;
    mpq_t qbpt;
    mpq_inits(qb, qbp, qbt, qbpt, (mpq_ptr)0);
    for (int i = 0; i < S; ++i) {
      mpq_set_ui(power, 1, 1);
      for (int e = 0; e < k; ++e) {
        mpq_mul(power, power, c[i]);
      }
      mpq_mul(term, b[i], power);
      mpq_add(qb, qb, term);
      mpq_mul(term, bp[i], power);
      mpq_add(qbp, qbp, term);
      mpq_mul(term, bt[i], power);
      mpq_add(qbt, qbt, term);
      mpq_mul(term, bpt[i], power);
      mpq_add(qbpt, qbpt, term);
    }
    if (k <= 10) {
      mpq_set_ui(target, 1, (unsigned long)((k + 1) * (k + 2)));
      MT_CHECK(mpq_equal(qb, target));
    }
    mpq_set_ui(target, 1, (unsigned long)(k + 1));
    MT_CHECK(mpq_equal(qbp, target));
    if (k <= 8) {
      MT_CHECK(mpq_sgn(qbt) == 0 && mpq_sgn(qbpt) == 0);
    }
    mpq_clears(qb, qbp, qbt, qbpt, (mpq_ptr)0);
  }
  mpq_clears(sum, term, power, target, (mpq_ptr)0);
  for (int k = 0; k < 6; ++k) {
    for (int i = 0; i < sizes[k]; ++i) {
      mpq_clear(all[k][i]);
    }
  }
  mpfr_free_cache();
  return mt_finish("test_mpfr_rkn12");
}
