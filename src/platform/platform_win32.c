/* SPDX-License-Identifier: BSD-3-Clause */
#define WIN32_LEAN_AND_MEAN
#include <intrin.h> /* __fastfail (declared noreturn) */
#include <io.h>     /* _write: unbuffered CRT write, no FILE lock */
#include <stdint.h>
#include <windows.h>

#include "core/alloc.h"
#include "core/base.h"
#include "platform/platform.h"

/* UTF-8 path -> NUL-terminated UTF-16 path in allocator memory. */
static myb_status wide_path(const char* path, const myb_allocator* allocator, wchar_t** out) {
  const int count = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, NULL, 0);
  if (count <= 0) {
    return MYB_E_IO;
  }
  void* block = NULL;
  const myb_status status =
      myb_internal_alloc_array(allocator, (size_t)count, sizeof(wchar_t), &block);
  if (status != MYB_OK) {
    return status;
  }
  wchar_t* wide = block;
  if (MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, wide, count) != count) {
    myb_internal_release(allocator, wide);
    return MYB_E_IO;
  }
  *out = wide;
  return MYB_OK;
}

/* Reads exactly size bytes; ReadFile transfers at most 2^32 - 1 bytes per call. */
static myb_status read_all(HANDLE file, char* buffer, size_t size) {
  enum { CHUNK = 1 << 30 };

  size_t done = 0;
  /* bound: each iteration reads a positive number of bytes or returns */
  for (size_t calls = 0; done < size && calls <= size / CHUNK + 1; ++calls) {
    const size_t remaining = size - done;
    const DWORD request = (DWORD)(remaining < CHUNK ? remaining : CHUNK);
    DWORD got = 0;
    if (!ReadFile(file, buffer + done, request, &got, NULL) || got == 0) {
      return MYB_E_IO;
    }
    done += got;
  }
  return done == size ? MYB_OK : MYB_E_IO;
}

myb_status myb_platform_read_file(const char* path, const myb_allocator* allocator, size_t max_size,
                                  char** data, size_t* size) {
  MYB_CHECK(path != NULL && allocator != NULL && data != NULL && size != NULL);
  wchar_t* wide = NULL;
  myb_status status = wide_path(path, allocator, &wide);
  if (status != MYB_OK) {
    return status;
  }
  HANDLE file = CreateFileW(wide, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                            FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN, NULL);
  myb_internal_release(allocator, wide);
  if (file == INVALID_HANDLE_VALUE) {
    return MYB_E_IO;
  }
  char* buffer = NULL;
  LARGE_INTEGER length;
  if (!GetFileSizeEx(file, &length) || length.QuadPart < 0) {
    status = MYB_E_IO;
    goto out;
  }
  if ((uint64_t)length.QuadPart > (uint64_t)max_size) {
    status = MYB_E_OVERFLOW;
    goto out;
  }
  const size_t bytes = (size_t)length.QuadPart;
  void* block = NULL;
  status = myb_internal_alloc_array(allocator, bytes + 1, 1, &block);
  if (status != MYB_OK) {
    goto out;
  }
  buffer = block;
  status = read_all(file, buffer, bytes);
  if (status == MYB_OK) {
    buffer[bytes] = '\0';
    *data = buffer;
    *size = bytes;
    buffer = NULL;
  }
out:
  myb_internal_release(allocator, buffer);
  CloseHandle(file);
  return status;
}

void myb_platform_write_stderr(const char* data, size_t size) {
  /* The fatal message is at most 1 KiB, so the conversion to unsigned is exact. */
  (void)_write(2, data, (unsigned)size);
}

void myb_platform_fatal_exit(void) {
  /* FAST_FAIL_FATAL_APP_EXIT: immediate termination without SEH handlers or dialogs; the
   * process exit code is 0xC0000409. */
  __fastfail(7);
}
