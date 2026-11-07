/* SPDX-License-Identifier: BSD-3-Clause */
/* Compiled with -mfma in addition to the library flags: a*b + c must still be two roundings,
 * because -ffp-contract=off comes last on the command line. */
#include "core/fpguard.h"

double myb_test_probe_mul_add(double a, double b, double c);

double myb_test_probe_mul_add(double a, double b, double c) {
  return a * b + c;
}
