/* SPDX-License-Identifier: BSD-3-Clause */
/* A failed MYB_CHECK must print one line "file:line: function: CHECK failed: condition" and
 * terminate. CTest matches the line (PASS_REGULAR_EXPRESSION). */
#include "core/check.h"

int main(void) {
  volatile int two = 2;
  MYB_CHECK(1 + 1 == 3 + (two - 2));
  return 0; /* not reached */
}
