/* SPDX-License-Identifier: BSD-3-Clause */
#include "num/series.h"

#include "core/base.h"
#include "myb/cmath.h"
#include "num/eft.h"

static int finite(double x) {
  return x - x == 0.0; /* false for infinities and NaN */
}

int myb_internal_series_valid(const double* terms, size_t count) {
  if (count != 0 && terms == NULL) {
    return 0;
  }
  for (size_t i = 0; i < count; ++i) {
    const double* term = terms + 4 * i;
    const double k = term[0];
    if (!(k >= 0.0 && k <= (double)MYB_SERIES_MAX_POWER && k == (double)(int)k)) {
      return 0;
    }
    if (!finite(term[1]) || !finite(term[2]) || !finite(term[3])) {
      return 0;
    }
  }
  return 1;
}

double myb_internal_series_eval(const myb_internal_series* series, double t) {
  MYB_CHECK(series != NULL);
  double powers[MYB_SERIES_MAX_POWER + 1];
  powers[0] = 1.0;
  for (int k = 1; k <= MYB_SERIES_MAX_POWER; ++k) {
    powers[k] = powers[k - 1] * t;
  }
  double sum = 0.0;
  double comp = 0.0;
  for (size_t i = 0; i < series->count; ++i) {
    const double* term = series->terms + 4 * i;
    const int k = (int)term[0];
    const double value = term[1] * powers[k] * myb_sin(term[2] * t + term[3]);
    double err;
    sum = myb_eft_two_sum(sum, value, &err);
    comp += err;
  }
  return sum + comp;
}
