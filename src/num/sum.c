/* SPDX-License-Identifier: BSD-3-Clause */
/* Compensated summation and dot product: Sum2 and Dot2 of T. Ogita, S. M. Rump, S. Oishi,
 * "Accurate sum and dot product", SIAM J. Sci. Comput. 26(6), 2005. */
#include "core/base.h"
#include "myb/num.h"
#include "num/eft.h"

void myb_sum_acc_add(myb_sum_acc* acc, double x) {
  MYB_CHECK(acc != NULL);
  double err;
  acc->sum = myb_eft_two_sum(acc->sum, x, &err);
  acc->comp += err;
}

double myb_sum_acc_value(const myb_sum_acc* acc) {
  MYB_CHECK(acc != NULL);
  return acc->sum + acc->comp;
}

double myb_sum2(const double* x, size_t n) {
  if (n == 0) {
    return 0.0;
  }
  MYB_CHECK(x != NULL);
  double sum = x[0];
  double comp = 0.0;
  for (size_t i = 1; i < n; ++i) {
    double err;
    sum = myb_eft_two_sum(sum, x[i], &err);
    comp += err;
  }
  return sum + comp;
}

double myb_dot2(const double* x, const double* y, size_t n) {
  if (n == 0) {
    return 0.0;
  }
  MYB_CHECK(x != NULL && y != NULL);
  double comp;
  double sum = myb_eft_two_prod(x[0], y[0], &comp);
  for (size_t i = 1; i < n; ++i) {
    double prod_err;
    const double prod = myb_eft_two_prod(x[i], y[i], &prod_err);
    double sum_err;
    sum = myb_eft_two_sum(sum, prod, &sum_err);
    comp += sum_err + prod_err;
  }
  return sum + comp;
}
