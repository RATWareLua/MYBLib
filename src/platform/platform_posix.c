/* SPDX-License-Identifier: BSD-3-Clause */
#include <fcntl.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "core/alloc.h"
#include "core/base.h"
#include "platform/platform.h"

static myb_status read_all(int fd, char* buffer, size_t size) {
  size_t done = 0;
  /* bound: every iteration reads at least one byte or returns, so at most size + 1 */
  for (size_t calls = 0; done < size && calls <= size; ++calls) {
    const ssize_t got = read(fd, buffer + done, size - done);
    if (got <= 0) {
      return MYB_E_IO;
    }
    done += (size_t)got;
  }
  return done == size ? MYB_OK : MYB_E_IO;
}

myb_status myb_platform_read_file(const char* path, const myb_allocator* allocator, size_t max_size,
                                  char** data, size_t* size) {
  MYB_CHECK(path != NULL && allocator != NULL && data != NULL && size != NULL);
  const int fd = open(path, O_RDONLY | O_CLOEXEC);
  if (fd < 0) {
    return MYB_E_IO;
  }
  myb_status status = MYB_E_IO;
  char* buffer = NULL;
  struct stat info;
  if (fstat(fd, &info) != 0 || !S_ISREG(info.st_mode) || info.st_size < 0) {
    goto out;
  }
  if ((uint64_t)info.st_size > (uint64_t)max_size) {
    status = MYB_E_OVERFLOW;
    goto out;
  }
  const size_t bytes = (size_t)info.st_size;
  void* block = NULL;
  status = myb_internal_alloc_array(allocator, bytes + 1, 1, &block);
  if (status != MYB_OK) {
    goto out;
  }
  buffer = block;
  status = read_all(fd, buffer, bytes);
  if (status == MYB_OK) {
    buffer[bytes] = '\0';
    *data = buffer;
    *size = bytes;
    buffer = NULL;
  }
out:
  myb_internal_release(allocator, buffer);
  close(fd);
  return status;
}

void myb_platform_write_stderr(const char* data, size_t size) {
  /* bound: at most 8 partial writes; the fatal path must not hang on a broken pipe */
  for (int attempt = 0; attempt < 8 && size > 0; ++attempt) {
    const ssize_t written = write(2, data, size);
    if (written <= 0) {
      return;
    }
    data += written;
    size -= (size_t)written;
  }
}

void myb_platform_fatal_exit(void) {
  abort();
}
