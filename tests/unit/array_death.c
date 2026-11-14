/* SPDX-License-Identifier: BSD-3-Clause */
/* An array form given NULL for an array of count > 0 elements breaks its contract: the
 * library prints one line naming the function and terminates (CTest matches the line). */
#include "myb/kepler.h"

int main(void) {
  const myb_pv pv0 = {{7.0e6, 0.0, 0.0}, {0.0, 7.5e3, 0.0}};
  myb_pv pv[2];
  volatile size_t count = 2;
  const myb_status status = myb_kepler_propagate_many(3.986004418e14, &pv0, NULL, count, pv, NULL);
  return status == MYB_OK ? 0 : 2; /* not reached */
}
