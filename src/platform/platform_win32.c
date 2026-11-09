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

myb_status myb_platform_map_file(const char* path, const myb_allocator* allocator,
                                 myb_platform_mapping* mapping) {
  MYB_CHECK(path != NULL && allocator != NULL && mapping != NULL);
  wchar_t* wide = NULL;
  myb_status status = wide_path(path, allocator, &wide);
  if (status != MYB_OK) {
    return status;
  }
  HANDLE file = CreateFileW(wide, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                            FILE_ATTRIBUTE_NORMAL | FILE_FLAG_RANDOM_ACCESS, NULL);
  myb_internal_release(allocator, wide);
  if (file == INVALID_HANDLE_VALUE) {
    return MYB_E_IO;
  }
  HANDLE section = NULL;
  LARGE_INTEGER length;
  status = MYB_E_IO;
  if (!GetFileSizeEx(file, &length) || length.QuadPart <= 0) {
    goto out;
  }
  if ((uint64_t)length.QuadPart > (uint64_t)SIZE_MAX) {
    status = MYB_E_OVERFLOW;
    goto out;
  }
  section = CreateFileMappingW(file, NULL, PAGE_READONLY, 0, 0, NULL);
  if (section == NULL) {
    goto out;
  }
  /* The view keeps the section and the file open after their handles are closed. */
  const void* view = MapViewOfFile(section, FILE_MAP_READ, 0, 0, 0);
  if (view != NULL) {
    mapping->data = view;
    mapping->size = (size_t)length.QuadPart;
    status = MYB_OK;
  }
out:
  if (section != NULL) {
    CloseHandle(section);
  }
  CloseHandle(file);
  return status;
}

void myb_platform_unmap_file(myb_platform_mapping* mapping) {
  MYB_CHECK(mapping != NULL);
  if (mapping->data != NULL) {
    MYB_CHECK(UnmapViewOfFile(mapping->data));
  }
  mapping->data = NULL;
  mapping->size = 0;
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

/* ---- Threads ----------------------------------------------------------------------------------
 * Slim reader/writer lock and condition variables (Windows Vista and later): no kernel object
 * per wait, and the workers stay parked between jobs. */
typedef struct team_worker {
  myb_platform_team* team;
  int member;
  HANDLE thread;
} team_worker;

struct myb_platform_team {
  myb_allocator allocator;
  int members;
  team_worker* workers; /* members - 1 */
  SRWLOCK lock;
  CONDITION_VARIABLE wake;
  CONDITION_VARIABLE done;
  uint64_t generation; /* incremented once per job */
  int pending;         /* workers that have not finished the current job */
  int busy;
  int quit;
  myb_platform_job* job;
  void* context;
};

static DWORD WINAPI team_worker_main(LPVOID parameter) {
  team_worker* w = parameter;
  myb_platform_team* t = w->team;
  uint64_t seen = 0;
  AcquireSRWLockExclusive(&t->lock);
  /* bound: one iteration per job; ends when destroy sets quit */
  for (;;) {
    /* bound: leaves on a new job or quit; a spurious wake-up re-tests the condition */
    while (t->generation == seen && !t->quit) {
      SleepConditionVariableSRW(&t->wake, &t->lock, INFINITE, 0);
    }
    if (t->quit) {
      break;
    }
    seen = t->generation;
    myb_platform_job* job = t->job;
    void* context = t->context;
    ReleaseSRWLockExclusive(&t->lock);
    job(context, w->member);
    AcquireSRWLockExclusive(&t->lock);
    t->pending--;
    if (t->pending == 0) {
      WakeConditionVariable(&t->done);
    }
  }
  ReleaseSRWLockExclusive(&t->lock);
  return 0;
}

/* Stops the first count workers and waits for them. */
static void team_stop(myb_platform_team* t, int count) {
  AcquireSRWLockExclusive(&t->lock);
  t->quit = 1;
  WakeAllConditionVariable(&t->wake);
  ReleaseSRWLockExclusive(&t->lock);
  for (int k = 0; k < count; ++k) {
    WaitForSingleObject(t->workers[k].thread, INFINITE);
    CloseHandle(t->workers[k].thread);
  }
}

myb_status myb_platform_team_create(int members, const myb_allocator* allocator,
                                    myb_platform_team** team) {
  MYB_CHECK(members >= 1 && allocator != NULL && team != NULL);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(allocator, 1, sizeof(myb_platform_team), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_platform_team* t = block;
  *t = (myb_platform_team){.allocator = *allocator, .members = members};
  InitializeSRWLock(&t->lock);
  InitializeConditionVariable(&t->wake);
  InitializeConditionVariable(&t->done);
  int started = 0;
  if (members > 1) {
    status = myb_internal_alloc_array(allocator, (size_t)(members - 1), sizeof(team_worker), &block);
    if (status != MYB_OK) {
      goto fail;
    }
    t->workers = block;
    for (int k = 1; k < members; ++k) {
      team_worker* w = &t->workers[k - 1];
      *w = (team_worker){.team = t, .member = k};
      w->thread = CreateThread(NULL, 0, team_worker_main, w, 0, NULL);
      if (w->thread == NULL) {
        status = MYB_E_IO;
        goto fail;
      }
      ++started;
    }
  }
  *team = t;
  return MYB_OK;

fail:
  team_stop(t, started);
  myb_internal_release(allocator, t->workers);
  myb_internal_release(allocator, t);
  return status;
}

void myb_platform_team_destroy(myb_platform_team* team) {
  if (team == NULL) {
    return;
  }
  MYB_CHECK(!team->busy);
  const myb_allocator allocator = team->allocator;
  team_stop(team, team->members - 1);
  myb_internal_release(&allocator, team->workers);
  myb_internal_release(&allocator, team);
}

void myb_platform_team_run(myb_platform_team* team, myb_platform_job* job, void* context) {
  MYB_CHECK(team != NULL && job != NULL);
  myb_platform_team* t = team;
  AcquireSRWLockExclusive(&t->lock);
  MYB_CHECK(!t->busy);
  t->busy = 1;
  t->job = job;
  t->context = context;
  t->pending = t->members - 1;
  t->generation++;
  WakeAllConditionVariable(&t->wake);
  ReleaseSRWLockExclusive(&t->lock);
  job(context, 0);
  AcquireSRWLockExclusive(&t->lock);
  /* bound: leaves when the members - 1 workers have finished this job */
  while (t->pending != 0) {
    SleepConditionVariableSRW(&t->done, &t->lock, INFINITE, 0);
  }
  t->busy = 0;
  ReleaseSRWLockExclusive(&t->lock);
}

int myb_platform_processor_count(void) {
  const DWORD count = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
  if (count == 0) {
    return 1;
  }
  return count > 4096 ? 4096 : (int)count;
}
