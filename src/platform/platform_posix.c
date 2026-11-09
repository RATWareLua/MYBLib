/* SPDX-License-Identifier: BSD-3-Clause */
/* POSIX.1-2008 (O_CLOEXEC, pthread, mmap): the library is built as strict C23, without the
 * C library's extensions. */
#define _POSIX_C_SOURCE 200809L

#include <fcntl.h>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/mman.h>
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

myb_status myb_platform_map_file(const char* path, const myb_allocator* allocator,
                                 myb_platform_mapping* mapping) {
  MYB_CHECK(path != NULL && allocator != NULL && mapping != NULL);
  const int fd = open(path, O_RDONLY | O_CLOEXEC);
  if (fd < 0) {
    return MYB_E_IO;
  }
  myb_status status = MYB_E_IO;
  struct stat info;
  if (fstat(fd, &info) != 0 || !S_ISREG(info.st_mode) || info.st_size <= 0) {
    goto out;
  }
  if ((uint64_t)info.st_size > (uint64_t)SIZE_MAX) {
    status = MYB_E_OVERFLOW;
    goto out;
  }
  const size_t bytes = (size_t)info.st_size;
  /* The mapping keeps the file referenced after the descriptor is closed. */
  void* view = mmap(NULL, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
  if (view != MAP_FAILED) {
    mapping->data = view;
    mapping->size = bytes;
    status = MYB_OK;
  }
out:
  close(fd);
  return status;
}

void myb_platform_unmap_file(myb_platform_mapping* mapping) {
  MYB_CHECK(mapping != NULL);
  if (mapping->data != NULL) {
    MYB_CHECK(munmap((void*)(uintptr_t)mapping->data, mapping->size) == 0);
  }
  mapping->data = NULL;
  mapping->size = 0;
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

/* ---- Threads --------------------------------------------------------------------------------- */
typedef struct team_worker {
  myb_platform_team* team;
  int member;
  pthread_t thread;
} team_worker;

struct myb_platform_team {
  myb_allocator allocator;
  int members;
  team_worker* workers; /* members - 1 */
  pthread_mutex_t lock;
  pthread_cond_t wake;
  pthread_cond_t done;
  uint64_t generation; /* incremented once per job */
  int pending;         /* workers that have not finished the current job */
  int busy;
  int quit;
  myb_platform_job* job;
  void* context;
};

static void* team_worker_main(void* parameter) {
  team_worker* w = parameter;
  myb_platform_team* t = w->team;
  uint64_t seen = 0;
  pthread_mutex_lock(&t->lock);
  /* bound: one iteration per job; ends when destroy sets quit */
  for (;;) {
    /* bound: leaves on a new job or quit; a spurious wake-up re-tests the condition */
    while (t->generation == seen && !t->quit) {
      pthread_cond_wait(&t->wake, &t->lock);
    }
    if (t->quit) {
      break;
    }
    seen = t->generation;
    myb_platform_job* job = t->job;
    void* context = t->context;
    pthread_mutex_unlock(&t->lock);
    job(context, w->member);
    pthread_mutex_lock(&t->lock);
    t->pending--;
    if (t->pending == 0) {
      pthread_cond_signal(&t->done);
    }
  }
  pthread_mutex_unlock(&t->lock);
  return NULL;
}

/* Stops the first count workers and waits for them. */
static void team_stop(myb_platform_team* t, int count) {
  pthread_mutex_lock(&t->lock);
  t->quit = 1;
  pthread_cond_broadcast(&t->wake);
  pthread_mutex_unlock(&t->lock);
  for (int k = 0; k < count; ++k) {
    pthread_join(t->workers[k].thread, NULL);
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
  int objects = 0; /* synchronization objects initialized: lock, wake, done */
  int started = 0;
  if (pthread_mutex_init(&t->lock, NULL) != 0) {
    status = MYB_E_IO;
    goto fail;
  }
  objects = 1;
  if (pthread_cond_init(&t->wake, NULL) != 0) {
    status = MYB_E_IO;
    goto fail;
  }
  objects = 2;
  if (pthread_cond_init(&t->done, NULL) != 0) {
    status = MYB_E_IO;
    goto fail;
  }
  objects = 3;
  if (members > 1) {
    status = myb_internal_alloc_array(allocator, (size_t)(members - 1), sizeof(team_worker), &block);
    if (status != MYB_OK) {
      goto fail;
    }
    t->workers = block;
    for (int k = 1; k < members; ++k) {
      team_worker* w = &t->workers[k - 1];
      *w = (team_worker){.team = t, .member = k};
      if (pthread_create(&w->thread, NULL, team_worker_main, w) != 0) {
        status = MYB_E_IO;
        goto fail;
      }
      ++started;
    }
  }
  *team = t;
  return MYB_OK;

fail:
  if (objects == 3) {
    team_stop(t, started);
    pthread_cond_destroy(&t->done);
  }
  if (objects >= 2) {
    pthread_cond_destroy(&t->wake);
  }
  if (objects >= 1) {
    pthread_mutex_destroy(&t->lock);
  }
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
  pthread_cond_destroy(&team->done);
  pthread_cond_destroy(&team->wake);
  pthread_mutex_destroy(&team->lock);
  myb_internal_release(&allocator, team->workers);
  myb_internal_release(&allocator, team);
}

void myb_platform_team_run(myb_platform_team* team, myb_platform_job* job, void* context) {
  MYB_CHECK(team != NULL && job != NULL);
  myb_platform_team* t = team;
  pthread_mutex_lock(&t->lock);
  MYB_CHECK(!t->busy);
  t->busy = 1;
  t->job = job;
  t->context = context;
  t->pending = t->members - 1;
  t->generation++;
  pthread_cond_broadcast(&t->wake);
  pthread_mutex_unlock(&t->lock);
  job(context, 0);
  pthread_mutex_lock(&t->lock);
  /* bound: leaves when the members - 1 workers have finished this job */
  while (t->pending != 0) {
    pthread_cond_wait(&t->done, &t->lock);
  }
  t->busy = 0;
  pthread_mutex_unlock(&t->lock);
}

int myb_platform_processor_count(void) {
  const long count = sysconf(_SC_NPROCESSORS_ONLN);
  if (count < 1) {
    return 1;
  }
  return count > 4096 ? 4096 : (int)count;
}
