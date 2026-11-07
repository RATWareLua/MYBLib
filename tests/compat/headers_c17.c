/* SPDX-License-Identifier: BSD-3-Clause */
/* Compiled as C17 (not C23): the public headers must not require the library's own standard. */
#include "myb/myb.h"

#if !defined(__STDC_VERSION__) || __STDC_VERSION__ != 201710L
#error "this file must be compiled as C17"
#endif

long myb_headers_c17_probe(void);

long myb_headers_c17_probe(void) {
  return MYB_VERSION_NUMBER;
}
