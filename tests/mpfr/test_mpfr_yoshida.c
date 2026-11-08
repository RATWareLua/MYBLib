/* SPDX-License-Identifier: BSD-3-Clause */
/* The weights of the triple-jump compositions (src/integ/symplectic_coefficients.h) recomputed
 * in MPFR at 512 bits and compared bit for bit: every entry must be the double nearest to the
 * exact product of the factors z1 = 1 / (2 - 2^(1/(p+1))), z0 = 1 - 2 z1. */
#include <mpfr.h>

#include "core/fpguard.h"
#include "integ/symplectic_coefficients.h"
#include "mtest.h"

enum { PREC = 512 };

int main(void) {
  mpfr_t w[27];
  mpfr_t next[27];
  mpfr_t z1;
  mpfr_t z0;
  mpfr_t t;
  for (int i = 0; i < 27; ++i) {
    mpfr_init2(w[i], PREC);
    mpfr_init2(next[i], PREC);
  }
  mpfr_inits2(PREC, z1, z0, t, (mpfr_ptr)0);
  mpfr_set_ui(w[0], 1, MPFR_RNDN);
  MT_CHECK_SAME(myb_symplectic_order2[0], 1.0);
  const double* tables[3] = {myb_symplectic_order4, myb_symplectic_order6, myb_symplectic_order8};
  int count = 1;
  for (int level = 0; level < 3; ++level) {
    const unsigned long root = (unsigned long)(2 * level + 3); /* p + 1 for p = 2, 4, 6 */
    mpfr_set_ui(t, 2, MPFR_RNDN);
    mpfr_rootn_ui(t, t, root, MPFR_RNDN);
    mpfr_ui_sub(t, 2, t, MPFR_RNDN);
    mpfr_ui_div(z1, 1, t, MPFR_RNDN);
    mpfr_mul_ui(z0, z1, 2, MPFR_RNDN);
    mpfr_ui_sub(z0, 1, z0, MPFR_RNDN);
    for (int i = 0; i < count; ++i) {
      mpfr_mul(next[i], w[i], z1, MPFR_RNDN);
      mpfr_mul(next[count + i], w[i], z0, MPFR_RNDN);
      mpfr_mul(next[2 * count + i], w[i], z1, MPFR_RNDN);
    }
    count *= 3;
    for (int i = 0; i < count; ++i) {
      mpfr_set(w[i], next[i], MPFR_RNDN);
      MT_CHECK_SAME(tables[level][i], mpfr_get_d(w[i], MPFR_RNDN));
    }
  }
  for (int i = 0; i < 27; ++i) {
    mpfr_clear(w[i]);
    mpfr_clear(next[i]);
  }
  mpfr_clears(z1, z0, t, (mpfr_ptr)0);
  mpfr_free_cache();
  return mt_finish("test_mpfr_yoshida");
}
