/* SPDX-License-Identifier: BSD-3-Clause */
#include "core/check.h"

#include <stddef.h>

#include "core/fpguard.h"
#include "platform/platform.h"

/* The message is assembled on the stack: the heap may already be corrupted when a check
 * fails, and stdio could block on a lock held by another thread. */
enum { MESSAGE_CAPACITY = 1024 };

typedef struct message {
  char data[MESSAGE_CAPACITY];
  size_t size;
} message;

static void message_put(message* m, const char* text) {
  if (text == NULL) {
    text = "?";
  }
  /* bound: length of text, at most MESSAGE_CAPACITY - 1 characters are stored */
  for (size_t i = 0; text[i] != '\0' && m->size < MESSAGE_CAPACITY - 1; ++i) {
    m->data[m->size++] = text[i];
  }
}

static void message_put_int(message* m, int value) {
  char digits[16];
  size_t count = 0;
  unsigned magnitude = value < 0 ? 0u - (unsigned)value : (unsigned)value;
  /* bound: an unsigned int has at most 10 decimal digits */
  for (size_t i = 0; i < sizeof digits; ++i) {
    digits[count++] = (char)('0' + (int)(magnitude % 10u));
    magnitude /= 10u;
    if (magnitude == 0u) {
      break;
    }
  }
  if (value < 0) {
    message_put(m, "-");
  }
  /* bound: count <= sizeof digits */
  for (size_t i = count; i > 0 && m->size < MESSAGE_CAPACITY - 1; --i) {
    m->data[m->size++] = digits[i - 1];
  }
}

void myb_internal_check_failed(const char* file, int line, const char* func, const char* what) {
  message m = {};
  message_put(&m, file);
  message_put(&m, ":");
  message_put_int(&m, line);
  message_put(&m, ": ");
  message_put(&m, func);
  message_put(&m, ": CHECK failed: ");
  message_put(&m, what);
  m.data[m.size++] = '\n';
  myb_platform_write_stderr(m.data, m.size);
  myb_platform_fatal_exit();
}
