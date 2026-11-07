// SPDX-License-Identifier: BSD-3-Clause
// Compiled as C++17: the public headers must be usable from C++.
#include "myb/myb.h"

static_assert(sizeof(myb_dd) == 2 * sizeof(double), "myb_dd layout");

long myb_headers_cxx17_probe();

long myb_headers_cxx17_probe() {
  return MYB_VERSION_NUMBER;
}
