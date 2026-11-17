/* SPDX-License-Identifier: BSD-3-Clause */
/* Compiled as C23 by another compiler than the library's (GCC, tests/tools/check_gcc.cmake):
 * the public headers, those of ext/ among them, with the attributes of the shared library. */
#include "myb/myb.h"
#include "myb_ionosphere.h"
#include "myb_jb2008.h"
#include "myb_nrlmsise00.h"
#include "myb_sgp4.h"
#include "myb_troposphere.h"

#if !defined(__STDC_VERSION__) || __STDC_VERSION__ <= 201710L
#error "this file must be compiled as C23"
#endif

long myb_headers_c23_probe(void);

long myb_headers_c23_probe(void) {
  return MYB_VERSION_NUMBER;
}
