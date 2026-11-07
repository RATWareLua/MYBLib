/* SPDX-License-Identifier: BSD-3-Clause */
#include "core/alloc.h"

#include <stdarg.h>
#include <stdckdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/base.h"

static void* default_alloc(void* context, size_t size) {
  (void)context;
  return malloc(size);
}

static void default_release(void* context, void* pointer) {
  (void)context;
  free(pointer);
}

myb_allocator myb_internal_allocator(const myb_allocator* requested) {
  if (requested != NULL) {
    MYB_CHECK(requested->alloc != NULL && requested->release != NULL);
    return *requested;
  }
  return (myb_allocator){.alloc = default_alloc, .release = default_release, .context = NULL};
}

myb_status myb_internal_alloc_array(const myb_allocator* allocator, size_t count, size_t size,
                                    void** out) {
  MYB_CHECK(allocator != NULL && out != NULL);
  size_t bytes;
  if (ckd_mul(&bytes, count, size)) {
    return MYB_E_OVERFLOW;
  }
  void* block = allocator->alloc(allocator->context, bytes == 0 ? 1 : bytes);
  if (block == NULL) {
    return MYB_E_NOMEM;
  }
  *out = block;
  return MYB_OK;
}

void myb_internal_release(const myb_allocator* allocator, void* pointer) {
  MYB_CHECK(allocator != NULL);
  if (pointer != NULL) {
    allocator->release(allocator->context, pointer);
  }
}

void myb_internal_diag(myb_diag* diag, long line, const char* format, ...) {
  if (diag == NULL) {
    return;
  }
  diag->line = line;
  va_list args;
  va_start(args, format);
  const int written = vsnprintf(diag->message, sizeof diag->message, format, args);
  va_end(args);
  if (written < 0) {
    diag->message[0] = '\0';
  }
  /* Some C runtimes (msvcrt) leave a truncated buffer unterminated. Messages contain only
   * integers and strings: floating-point formatting would depend on the locale. */
  diag->message[sizeof diag->message - 1] = '\0';
}

void myb_internal_diag_clear(myb_diag* diag) {
  if (diag != NULL) {
    diag->line = 0;
    diag->message[0] = '\0';
  }
}
