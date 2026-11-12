/* SPDX-License-Identifier: BSD-3-Clause */
/* The compositions of 7 and 15 stages (src/integ/symplectic_coefficients.h). Their order
 * conditions are derived here: for a symmetric method of order 2, log S(h) = h X + h^3 Y3 +
 * h^5 Y5 + h^7 Y7 + ... with free X, Y3, Y5, Y7 (R. I. McLachlan, SIAM J. Sci. Comput. 16, 151,
 * 1995), and the composition S(g_s h) ... S(g_1 h) has order p when the logarithm of the
 * product (the formula of Baker, Campbell and Hausdorff, computed in the free associative
 * algebra truncated at weight 8: 55 words) is h X up to weight p. From the tabulated doubles
 * the free weights are refined by Gauss-Newton at 512 bits until every condition holds to
 * 2^-470; then every double must be the nearest to its refined weight, bit for bit. */
#include <math.h>
#include <mpfr.h>
#include <stdio.h>

#include "core/fpguard.h"
#include "integ/symplectic_coefficients.h"
#include "mtest.h"

enum { PREC = 512, MAXW = 8, MAXWORDS = 55, MAXLEN = 8, LETTERS = 4, MAXHALF = 7 };

static const int letter_weight[LETTERS] = {1, 3, 5, 7}; /* X, Y3, Y5, Y7 */

typedef struct word {
  int length;
  int weight;
  int letters[MAXLEN];
} word;

static word words[MAXWORDS];
static int word_count;
static int product[MAXWORDS][MAXWORDS]; /* the index of the concatenation, -1 beyond MAXW */

static int find_word(const int* letters, int length) {
  for (int k = 0; k < word_count; ++k) {
    if (words[k].length != length) {
      continue;
    }
    int same = 1;
    for (int i = 0; i < length; ++i) {
      same = same && words[k].letters[i] == letters[i];
    }
    if (same) {
      return k;
    }
  }
  return -1;
}

/* All words of weight <= MAXW, by increasing length (no recursion: layer by layer). */
static void make_words(void) {
  word_count = 1;
  words[0] = (word){0};
  int begin = 0;
  int end = 1;
  for (int layer = 0; layer < MAXLEN; ++layer) {
    for (int k = begin; k < end; ++k) {
      for (int l = 0; l < LETTERS; ++l) {
        if (words[k].weight + letter_weight[l] <= MAXW) {
          word w = words[k];
          w.letters[w.length++] = l;
          w.weight += letter_weight[l];
          MT_CHECK(word_count < MAXWORDS);
          if (word_count < MAXWORDS) {
            words[word_count++] = w;
          }
        }
      }
    }
    begin = end;
    end = word_count;
  }
  for (int a = 0; a < word_count; ++a) {
    for (int b = 0; b < word_count; ++b) {
      product[a][b] = -1;
      if (words[a].weight + words[b].weight <= MAXW) {
        int letters[MAXLEN];
        for (int i = 0; i < words[a].length; ++i) {
          letters[i] = words[a].letters[i];
        }
        for (int i = 0; i < words[b].length; ++i) {
          letters[words[a].length + i] = words[b].letters[i];
        }
        product[a][b] = find_word(letters, words[a].length + words[b].length);
      }
    }
  }
}

/* An element of the truncated algebra: a coefficient per word. */
typedef struct element {
  mpfr_t c[MAXWORDS];
} element;

static void element_init(element* e) {
  for (int k = 0; k < MAXWORDS; ++k) {
    mpfr_init2(e->c[k], PREC);
    mpfr_set_zero(e->c[k], 1);
  }
}

static void element_clear(element* e) {
  for (int k = 0; k < MAXWORDS; ++k) {
    mpfr_clear(e->c[k]);
  }
}

static void element_set(element* to, const element* from) {
  for (int k = 0; k < word_count; ++k) {
    mpfr_set(to->c[k], from->c[k], MPFR_RNDN);
  }
}

static mpfr_t scratch;

/* out = a b (out distinct from a and b) */
static void element_mul(element* out, const element* a, const element* b) {
  for (int k = 0; k < word_count; ++k) {
    mpfr_set_zero(out->c[k], 1);
  }
  for (int i = 0; i < word_count; ++i) {
    if (mpfr_zero_p(a->c[i])) {
      continue;
    }
    for (int j = 0; j < word_count; ++j) {
      const int k = product[i][j];
      if (k >= 0 && !mpfr_zero_p(b->c[j])) {
        mpfr_mul(scratch, a->c[i], b->c[j], MPFR_RNDN);
        mpfr_add(out->c[k], out->c[k], scratch, MPFR_RNDN);
      }
    }
  }
}

static element term;
static element power;
static element temporary;
static element exponent;
static element stage;
static element whole;
static element logarithm;

/* stage = exp(g X + g^3 Y3 + g^5 Y5 + g^7 Y7) */
static void stage_of(const mpfr_t g) {
  for (int k = 0; k < word_count; ++k) {
    mpfr_set_zero(exponent.c[k], 1);
  }
  for (int l = 0; l < LETTERS; ++l) {
    const int letters[1] = {l};
    mpfr_pow_ui(exponent.c[find_word(letters, 1)], g, (unsigned long)letter_weight[l],
                MPFR_RNDN);
  }
  for (int k = 0; k < word_count; ++k) {
    mpfr_set_zero(stage.c[k], 1);
    mpfr_set_zero(term.c[k], 1);
  }
  mpfr_set_ui(stage.c[0], 1, MPFR_RNDN);
  mpfr_set_ui(term.c[0], 1, MPFR_RNDN);
  for (int n = 1; n <= MAXW; ++n) {
    element_mul(&temporary, &term, &exponent);
    for (int k = 0; k < word_count; ++k) {
      mpfr_div_ui(term.c[k], temporary.c[k], (unsigned long)n, MPFR_RNDN);
      mpfr_add(stage.c[k], stage.c[k], term.c[k], MPFR_RNDN);
    }
  }
}

/* logarithm = log(whole), whole = 1 + q */
static void log_of_whole(void) {
  element_set(&term, &whole);
  mpfr_sub_ui(term.c[0], term.c[0], 1, MPFR_RNDN); /* q */
  for (int k = 0; k < word_count; ++k) {
    mpfr_set_zero(logarithm.c[k], 1);
    mpfr_set_zero(power.c[k], 1);
  }
  mpfr_set_ui(power.c[0], 1, MPFR_RNDN);
  for (int n = 1; n <= MAXW; ++n) {
    element_mul(&temporary, &power, &term);
    element_set(&power, &temporary);
    for (int k = 0; k < word_count; ++k) {
      mpfr_div_si(scratch, power.c[k], n % 2 == 1 ? n : -n, MPFR_RNDN);
      mpfr_add(logarithm.c[k], logarithm.c[k], scratch, MPFR_RNDN);
    }
  }
}

/* The residuals of the conditions of order p for the symmetric weights with the free half
 * `half` (outer to inner; the middle 1 - 2 sum): the coefficient of every word of weight
 * 1 .. p of the logarithm, less 1 for X. */
static int residuals(const mpfr_t* half, int m, int p, mpfr_t* r) {
  mpfr_t g[2 * MAXHALF + 1];
  const int s = 2 * m + 1;
  for (int k = 0; k < s; ++k) {
    mpfr_init2(g[k], PREC);
  }
  mpfr_set_ui(g[m], 1, MPFR_RNDN);
  for (int k = 0; k < m; ++k) {
    mpfr_set(g[k], half[k], MPFR_RNDN);
    mpfr_set(g[s - 1 - k], half[k], MPFR_RNDN);
    mpfr_mul_ui(scratch, half[k], 2, MPFR_RNDN);
    mpfr_sub(g[m], g[m], scratch, MPFR_RNDN);
  }
  for (int k = 0; k < word_count; ++k) {
    mpfr_set_zero(whole.c[k], 1);
  }
  mpfr_set_ui(whole.c[0], 1, MPFR_RNDN);
  for (int k = 0; k < s; ++k) {
    stage_of(g[k]);
    element_mul(&temporary, &stage, &whole);
    element_set(&whole, &temporary);
  }
  log_of_whole();
  int count = 0;
  for (int k = 1; k < word_count; ++k) {
    if (words[k].weight <= p) {
      mpfr_set(r[count], logarithm.c[k], MPFR_RNDN);
      if (words[k].length == 1 && words[k].letters[0] == 0) {
        mpfr_sub_ui(r[count], r[count], 1, MPFR_RNDN);
      }
      ++count;
    }
  }
  for (int k = 0; k < s; ++k) {
    mpfr_clear(g[k]);
  }
  return count;
}

static mpfr_t res[MAXWORDS];
static mpfr_t shifted[MAXWORDS];
static mpfr_t jac[MAXWORDS][MAXHALF];
static mpfr_t normal[MAXHALF][MAXHALF + 1];

/* Gauss-Newton from the doubles of a table; returns log2 of the largest residual at the end. */
static double refine(const double* table, int s, int p, mpfr_t* half) {
  const int m = (s - 1) / 2;
  for (int k = 0; k < m; ++k) {
    mpfr_set_d(half[k], table[k], MPFR_RNDN);
  }
  mpfr_t eps;
  mpfr_t largest;
  mpfr_init2(eps, PREC);
  mpfr_init2(largest, PREC);
  mpfr_set_ui_2exp(eps, 1, -256, MPFR_RNDN);
  int count = residuals(half, m, p, res);
  double first = 0.0;
  for (int iteration = 0; iteration < 12; ++iteration) {
    mpfr_set_zero(largest, 1);
    for (int i = 0; i < count; ++i) {
      mpfr_abs(scratch, res[i], MPFR_RNDN);
      mpfr_max(largest, largest, scratch, MPFR_RNDN);
    }
    if (iteration == 0) {
      first = mpfr_get_d(largest, MPFR_RNDN);
    }
    if (mpfr_cmp_ui_2exp(largest, 1, -470) < 0) {
      break;
    }
    /* the Jacobian by forward differences of 2^-256 */
    for (int j = 0; j < m; ++j) {
      mpfr_add(half[j], half[j], eps, MPFR_RNDN);
      residuals(half, m, p, shifted);
      mpfr_sub(half[j], half[j], eps, MPFR_RNDN);
      for (int i = 0; i < count; ++i) {
        mpfr_sub(jac[i][j], shifted[i], res[i], MPFR_RNDN);
        mpfr_div(jac[i][j], jac[i][j], eps, MPFR_RNDN);
      }
    }
    /* the normal equations J^T J d = J^T r, by elimination with partial pivoting */
    for (int a = 0; a < m; ++a) {
      for (int b = 0; b <= m; ++b) {
        mpfr_set_zero(normal[a][b], 1);
        for (int i = 0; i < count; ++i) {
          mpfr_mul(scratch, jac[i][a], b < m ? jac[i][b] : res[i], MPFR_RNDN);
          mpfr_add(normal[a][b], normal[a][b], scratch, MPFR_RNDN);
        }
      }
    }
    for (int col = 0; col < m; ++col) {
      int pivot = col;
      for (int row = col + 1; row < m; ++row) {
        if (mpfr_cmpabs(normal[row][col], normal[pivot][col]) > 0) {
          pivot = row;
        }
      }
      for (int b = 0; b <= m; ++b) {
        mpfr_swap(normal[col][b], normal[pivot][b]);
      }
      for (int row = col + 1; row < m; ++row) {
        mpfr_t factor;
        mpfr_init2(factor, PREC);
        mpfr_div(factor, normal[row][col], normal[col][col], MPFR_RNDN);
        for (int b = col; b <= m; ++b) {
          mpfr_mul(scratch, factor, normal[col][b], MPFR_RNDN);
          mpfr_sub(normal[row][b], normal[row][b], scratch, MPFR_RNDN);
        }
        mpfr_clear(factor);
      }
    }
    for (int row = m - 1; row >= 0; --row) {
      for (int b = row + 1; b < m; ++b) {
        mpfr_mul(scratch, normal[row][b], normal[b][m], MPFR_RNDN);
        mpfr_sub(normal[row][m], normal[row][m], scratch, MPFR_RNDN);
      }
      mpfr_div(normal[row][m], normal[row][m], normal[row][row], MPFR_RNDN);
      mpfr_sub(half[row], half[row], normal[row][m], MPFR_RNDN);
    }
    count = residuals(half, m, p, res);
  }
  const double last = mpfr_get_d(largest, MPFR_RNDN);
  printf("  %d stages, order %d: %d conditions; the doubles meet them to %.1e, the refined "
         "weights to 2^%.0f\n",
         s, p, count, first, last > 0.0 ? log2(last) : -9999.0);
  mpfr_clear(eps);
  mpfr_clear(largest);
  return last > 0.0 ? log2(last) : -9999.0;
}

int main(void) {
  make_words();
  MT_CHECK(word_count == 55);
  mpfr_init2(scratch, PREC);
  element* all[7] = {&term, &power, &temporary, &exponent, &stage, &whole, &logarithm};
  for (int k = 0; k < 7; ++k) {
    element_init(all[k]);
  }
  for (int i = 0; i < MAXWORDS; ++i) {
    mpfr_init2(res[i], PREC);
    mpfr_init2(shifted[i], PREC);
    for (int j = 0; j < MAXHALF; ++j) {
      mpfr_init2(jac[i][j], PREC);
    }
  }
  for (int i = 0; i < MAXHALF; ++i) {
    for (int j = 0; j <= MAXHALF; ++j) {
      mpfr_init2(normal[i][j], PREC);
    }
  }
  const struct {
    const double* table;
    int stages;
    int order;
  } sets[2] = {{myb_symplectic_order6_s7, 7, 6}, {myb_symplectic_order8_s15, 15, 8}};
  for (int k = 0; k < 2; ++k) {
    mpfr_t half[MAXHALF];
    mpfr_t middle;
    for (int i = 0; i < MAXHALF; ++i) {
      mpfr_init2(half[i], PREC);
    }
    mpfr_init2(middle, PREC);
    const double exponent_left = refine(sets[k].table, sets[k].stages, sets[k].order, half);
    MT_CHECK(exponent_left < -470.0);
    const int m = (sets[k].stages - 1) / 2;
    mpfr_set_ui(middle, 1, MPFR_RNDN);
    for (int i = 0; i < m; ++i) {
      MT_CHECK_SAME(sets[k].table[i], mpfr_get_d(half[i], MPFR_RNDN));
      MT_CHECK_SAME(sets[k].table[sets[k].stages - 1 - i], mpfr_get_d(half[i], MPFR_RNDN));
      mpfr_mul_ui(scratch, half[i], 2, MPFR_RNDN);
      mpfr_sub(middle, middle, scratch, MPFR_RNDN);
    }
    MT_CHECK_SAME(sets[k].table[m], mpfr_get_d(middle, MPFR_RNDN));
    for (int i = 0; i < MAXHALF; ++i) {
      mpfr_clear(half[i]);
    }
    mpfr_clear(middle);
  }
  for (int k = 0; k < 7; ++k) {
    element_clear(all[k]);
  }
  for (int i = 0; i < MAXWORDS; ++i) {
    mpfr_clear(res[i]);
    mpfr_clear(shifted[i]);
    for (int j = 0; j < MAXHALF; ++j) {
      mpfr_clear(jac[i][j]);
    }
  }
  for (int i = 0; i < MAXHALF; ++i) {
    for (int j = 0; j <= MAXHALF; ++j) {
      mpfr_clear(normal[i][j]);
    }
  }
  mpfr_clear(scratch);
  mpfr_free_cache();
  return mt_finish("test_mpfr_compositions");
}
