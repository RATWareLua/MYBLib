/* SPDX-License-Identifier: BSD-3-Clause */
/* Error-free transformations are exact, double-word operations stay within their proven
 * relative error bounds, Sum2/Dot2 within theirs — all measured against GNU MPFR. The
 * maximum observed errors are printed (in units of u^2 = 2^-106). */
#include <math.h>
#include <mpfr.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/num.h"

enum { EXACT_PREC = 2400, REF_PREC = 320 };

static const double u = 0x1p-53;

/* ---- Error-free transformations -------------------------------------------------------- */

static void test_eft_exact(void) {
  mpfr_t exact;
  mpfr_t got;
  mpfr_inits2(EXACT_PREC, exact, got, (mpfr_ptr)0);
  mt_rng rng = {.state = 31};
  long sum_wrong = 0;
  long prod_wrong = 0;
  for (int i = 0; i < 200000; ++i) {
    const double a = mt_double(&rng, -1000, 1000);
    const double b = mt_double(&rng, -1000, 1000);
    double err;
    const double s = myb_two_sum(a, b, &err);
    if (isfinite(s)) {
      mpfr_set_d(exact, a, MPFR_RNDN);
      mpfr_add_d(exact, exact, b, MPFR_RNDN); /* exact at EXACT_PREC */
      mpfr_set_d(got, s, MPFR_RNDN);
      mpfr_add_d(got, got, err, MPFR_RNDN);
      mtest.checks++;
      if (mpfr_cmp(exact, got) != 0) {
        mtest.failures++;
        sum_wrong++;
      }
    }
    const double a2 = mt_double(&rng, -480, 480);
    const double b2 = mt_double(&rng, -480, 480);
    const double p = myb_two_prod(a2, b2, &err);
    mpfr_set_d(exact, a2, MPFR_RNDN);
    mpfr_mul_d(exact, exact, b2, MPFR_RNDN);
    mpfr_set_d(got, p, MPFR_RNDN);
    mpfr_add_d(got, got, err, MPFR_RNDN);
    mtest.checks++;
    if (mpfr_cmp(exact, got) != 0) {
      mtest.failures++;
      prod_wrong++;
    }
  }
  mpfr_clears(exact, got, (mpfr_ptr)0);
  printf("two_sum  200000 inputs, %ld not exact\n", sum_wrong);
  printf("two_prod 200000 inputs, %ld not exact\n", prod_wrong);
}

/* ---- Double-word operations ------------------------------------------------------------- */

typedef enum op { OP_ADD, OP_SUB, OP_ADD_D, OP_MUL, OP_MUL_D, OP_DIV, OP_DIV_D, OP_SQRT } op;

typedef struct op_info {
  op kind;
  const char* name;
  double bound_u2; /* proven bound in units of u^2 (leading term; higher terms below 1e-12) */
} op_info;

/* Random normalised double-word number: hi with exponent in [emin, emax], lo up to ulp/2. */
static myb_dd random_dd(mt_rng* rng, int emin, int emax) {
  const double hi = mt_double(rng, emin, emax);
  const double lo = mt_double(rng, -1, -1) * hi * u * 0.5;
  double err;
  const double h = myb_two_sum(hi, lo, &err);
  return (myb_dd){.hi = h, .lo = err};
}

static void set_dd(mpfr_ptr target, myb_dd x) {
  mpfr_set_d(target, x.hi, MPFR_RNDN);
  mpfr_add_d(target, target, x.lo, MPFR_RNDN);
}

/* Relative error |z - ref| / |ref| in units of u^2. */
static double relative_error_u2(myb_dd z, mpfr_srcptr ref) {
  mpfr_t diff;
  mpfr_init2(diff, REF_PREC);
  set_dd(diff, z);
  mpfr_sub(diff, diff, ref, MPFR_RNDN);
  mpfr_div(diff, diff, ref, MPFR_RNDN);
  mpfr_abs(diff, diff, MPFR_RNDN);
  mpfr_mul_2si(diff, diff, 106, MPFR_RNDN);
  const double e = mpfr_get_d(diff, MPFR_RNDU);
  mpfr_clear(diff);
  return e;
}

static myb_dd apply(op kind, myb_dd x, myb_dd y) {
  switch (kind) {
    case OP_ADD:
      return myb_dd_add(x, y);
    case OP_SUB:
      return myb_dd_sub(x, y);
    case OP_ADD_D:
      return myb_dd_add_d(x, y.hi);
    case OP_MUL:
      return myb_dd_mul(x, y);
    case OP_MUL_D:
      return myb_dd_mul_d(x, y.hi);
    case OP_DIV:
      return myb_dd_div(x, y);
    case OP_DIV_D:
      return myb_dd_div_d(x, y.hi);
    case OP_SQRT:
      return myb_dd_sqrt(x);
  }
  return x;
}

static void reference(op kind, mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr y, double y_hi) {
  switch (kind) {
    case OP_ADD:
      mpfr_add(r, x, y, MPFR_RNDN);
      return;
    case OP_SUB:
      mpfr_sub(r, x, y, MPFR_RNDN);
      return;
    case OP_ADD_D:
      mpfr_add_d(r, x, y_hi, MPFR_RNDN);
      return;
    case OP_MUL:
      mpfr_mul(r, x, y, MPFR_RNDN);
      return;
    case OP_MUL_D:
      mpfr_mul_d(r, x, y_hi, MPFR_RNDN);
      return;
    case OP_DIV:
      mpfr_div(r, x, y, MPFR_RNDN);
      return;
    case OP_DIV_D:
      mpfr_div_d(r, x, y_hi, MPFR_RNDN);
      return;
    case OP_SQRT:
      mpfr_sqrt(r, x, MPFR_RNDN);
      return;
  }
}

static void test_dd_op(const op_info* info, uint64_t seed) {
  mpfr_t x;
  mpfr_t y;
  mpfr_t r;
  mpfr_inits2(REF_PREC, x, y, r, (mpfr_ptr)0);
  mt_rng rng = {.state = seed};
  double max_error = 0.0;
  long violations = 0;
  const long count = 300000;
  for (long i = 0; i < count; ++i) {
    myb_dd a = random_dd(&rng, -40, 40);
    myb_dd b = random_dd(&rng, -40, 40);
    if (info->kind == OP_SQRT && a.hi < 0.0) {
      a = myb_dd_neg(a);
    }
    if ((info->kind == OP_ADD || info->kind == OP_SUB) && (i % 4) == 0) {
      /* Deep cancellation: b is -a (or a) perturbed in its last bits. */
      const double perturb = mt_double(&rng, -1, -1) * a.hi * 0x1p-80;
      b = myb_dd_add_d(info->kind == OP_ADD ? myb_dd_neg(a) : a, perturb);
    }
    if (info->kind == OP_ADD_D || info->kind == OP_MUL_D || info->kind == OP_DIV_D) {
      b.lo = 0.0;
    }
    const myb_dd z = apply(info->kind, a, b);
    set_dd(x, a);
    set_dd(y, b);
    reference(info->kind, r, x, y, b.hi);
    if (mpfr_zero_p(r)) {
      continue;
    }
    const double e = relative_error_u2(z, r);
    if (e > max_error) {
      max_error = e;
    }
    mtest.checks++;
    if (!(e <= info->bound_u2 * (1.0 + 1e-12) + 1e-12)) {
      mtest.failures++;
      if (++violations <= 3) {
        fprintf(stderr, "%s: (%a + %a) op (%a + %a): error %.3f u^2 > bound %.3f\n", info->name,
                a.hi, a.lo, b.hi, b.lo, e, info->bound_u2);
      }
    }
  }
  mpfr_clears(x, y, r, (mpfr_ptr)0);
  printf("dd_%-6s %ld inputs, max error %.4f u^2 (bound %.2f u^2)\n", info->name, count, max_error,
         info->bound_u2);
}

/* ---- Sum2, Dot2 ---------------------------------------------------------------------------- */

static void test_sum_dot(void) {
  enum { N = 1000 };

  static double x[N];
  static double y[N];
  mpfr_t exact;
  mpfr_t term;
  mpfr_t abs_sum;
  mpfr_inits2(EXACT_PREC, exact, term, abs_sum, (mpfr_ptr)0);
  mt_rng rng = {.state = 41};
  const double gamma_n = (double)N * u / (1.0 - (double)N * u);
  long wrong = 0;
  for (int trial = 0; trial < 200; ++trial) {
    /* Ill-conditioned: pairs r, -(r + s) with s much smaller than r. */
    for (int i = 0; i < N; i += 2) {
      const double r = mt_double(&rng, 0, 60);
      const double s = mt_double(&rng, -40, -10);
      x[i] = r;
      x[i + 1] = -(r + s);
      y[i] = mt_double(&rng, -2, 2);
      y[i + 1] = y[i];
    }
    mpfr_set_zero(exact, 1);
    mpfr_set_zero(abs_sum, 1);
    for (int i = 0; i < N; ++i) {
      mpfr_add_d(exact, exact, x[i], MPFR_RNDN);
      mpfr_set_d(term, fabs(x[i]), MPFR_RNDN);
      mpfr_add(abs_sum, abs_sum, term, MPFR_RNDN);
    }
    const double sum = myb_sum2(x, N);
    const double s_exact = mpfr_get_d(exact, MPFR_RNDN);
    const double bound_sum = u * fabs(s_exact) + gamma_n * gamma_n * mpfr_get_d(abs_sum, MPFR_RNDU);
    mpfr_sub_d(term, exact, sum, MPFR_RNDN);
    mpfr_abs(term, term, MPFR_RNDN);
    mtest.checks++;
    if (!(mpfr_get_d(term, MPFR_RNDU) <= bound_sum)) {
      mtest.failures++;
      wrong++;
    }

    mpfr_set_zero(exact, 1);
    mpfr_set_zero(abs_sum, 1);
    for (int i = 0; i < N; ++i) {
      mpfr_set_d(term, x[i], MPFR_RNDN);
      mpfr_mul_d(term, term, y[i], MPFR_RNDN);
      mpfr_add(exact, exact, term, MPFR_RNDN);
      mpfr_abs(term, term, MPFR_RNDN);
      mpfr_add(abs_sum, abs_sum, term, MPFR_RNDN);
    }
    const double dot = myb_dot2(x, y, N);
    const double d_exact = mpfr_get_d(exact, MPFR_RNDN);
    const double bound_dot = u * fabs(d_exact) + gamma_n * gamma_n * mpfr_get_d(abs_sum, MPFR_RNDU);
    mpfr_sub_d(term, exact, dot, MPFR_RNDN);
    mpfr_abs(term, term, MPFR_RNDN);
    mtest.checks++;
    if (!(mpfr_get_d(term, MPFR_RNDU) <= bound_dot)) {
      mtest.failures++;
      wrong++;
    }
  }
  mpfr_clears(exact, term, abs_sum, (mpfr_ptr)0);
  printf("sum2/dot2 200 ill-conditioned trials each, %ld bound violations\n", wrong);
}

int main(void) {
  test_eft_exact();
  const op_info ops[] = {
      {OP_ADD, "add", 3.0},     {OP_SUB, "sub", 3.0},     {OP_ADD_D, "add_d", 2.0},
      {OP_MUL, "mul", 7.0},     {OP_MUL_D, "mul_d", 1.5}, {OP_DIV, "div", 15.0},
      {OP_DIV_D, "div_d", 3.5}, {OP_SQRT, "sqrt", 4.0},
  };
  for (size_t i = 0; i < sizeof ops / sizeof ops[0]; ++i) {
    test_dd_op(&ops[i], 1000 + i);
  }
  test_sum_dot();
  return mt_finish("test_mpfr_num");
}
