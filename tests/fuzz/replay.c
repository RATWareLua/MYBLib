/* SPDX-License-Identifier: BSD-3-Clause */
/* replay.c — the main program of a fuzz harness where libFuzzer is not linked: every ordinary
 * build, and Windows, whose Clang has no libFuzzer runtime for the mingw-w64 target. The files
 * named on the command line run through LLVMFuzzerTestOneInput once each, as libFuzzer runs
 * the files it is given; with -mutate, also variants of each made by a fixed sequence of
 * small damages, so that a build without libFuzzer still feeds the readers broken input:
 *
 *   fuzz_<name> [-mutate=N] [-seed=S] [-variant=K] [-save=PATH] [-verbose] FILE... [@LIST]...
 *
 *   -mutate=N   for each file also its truncations (at 32 places through it and at each of
 *               its last 8 bytes) and N variants with one to four damages each: a byte
 *               replaced by a character that means something to a reader ('\n', '-', '=',
 *               '<', a digit, NUL) or by any byte, a span deleted or repeated, a number
 *               inserted (an overflow, a NaN, a subnormal), a word of a binary file replaced
 *               by an extreme integer or double. Variant 0 is the file itself.
 *   -seed=S     the seed of the variants (default 1).
 *   -variant=K  only variant K of each file (to reproduce a failure).
 *   -save=PATH  writes the last variant to PATH before running it (with -variant: the input
 *               that failed, to keep as a regression case in tests/fuzz/corpus/<name>/).
 *   -verbose    names every variant on stderr before it runs (the last name printed is the
 *               one that failed).
 *   @LIST       a file of input paths, one per line.
 *
 * Each input is copied into a block of exactly its size, so that a reader that reads past
 * the end is caught where the build has a sanitizer of addresses (linux-asan) and crashes
 * nowhere else by luck of a terminating NUL. A failed check, a crash or a hang (the timeout
 * of the test) ends the process: the exit code is 0 when every input ran. Files are read
 * through the platform layer of the library (UTF-8 paths everywhere). */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/alloc.h"
#include "fuzz.h"
#include "platform/platform.h"

enum {
  MAX_INPUT = 256 * 1024 * 1024,
  TRUNCATIONS = 32 + 8,
  MAX_GROWTH = 4 * 64 + 64, /* the most a variant adds: four insertions */
  MAX_ARGUMENT_FILES = 4096
};

typedef struct options {
  long mutate;
  uint64_t seed;
  long variant; /* -1: all */
  const char* save;
  int verbose;
} options;

/* ---- Deterministic variants ------------------------------------------------------------------ */

/* SplitMix64 (Steele, Lea, Flood, OOPSLA 2014). */
static uint64_t next_random(uint64_t* state) {
  uint64_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}

static size_t below(uint64_t* state, size_t n) {
  return n == 0 ? 0 : (size_t)(next_random(state) % n);
}

static const char SPECIAL[] = {'\0', '\n', '\r', ' ', '\t', '-', '+', '.', '0', '9', 'e', 'E',
                               'D', '=', '(', ')', '\'', '<', '>', '/', '&', ';', ',', '#', '$',
                               '@', 'x', 'Z', '[', ']', '"', '!', '\\', '\x80', '\xff'};

static const char* const NUMBERS[] = {"0", "-1", "1e308", "-1e309", "nan", "inf", "-0",
                                      "99999999999999999999", "2147483648", "-2147483649",
                                      "4.9e-324", "0x1p-1074", "1e-400", "18446744073709551616",
                                      "1.", ".5e", "3e+", "1D+99", "7777777777"};

/* Extreme words of a binary file: 32-bit integers and doubles, little-endian. */
static const uint64_t WORDS[] = {0,
                                 UINT64_C(0xFFFFFFFFFFFFFFFF),
                                 UINT64_C(0x000000007FFFFFFF),
                                 UINT64_C(0x0000000080000000),
                                 UINT64_C(0x0000000000000001),
                                 UINT64_C(0x0000000000000400),
                                 UINT64_C(0x7FF8000000000000), /* NaN */
                                 UINT64_C(0x7FF0000000000000), /* inf */
                                 UINT64_C(0xFFF0000000000000), /* -inf */
                                 UINT64_C(0x7FEFFFFFFFFFFFFF), /* the largest double */
                                 UINT64_C(0x8000000000000000), /* -0 */
                                 UINT64_C(0x0000000000000001), /* the smallest subnormal */
                                 UINT64_C(0x4340000000000000), /* 2^53 */
                                 UINT64_C(0x41DFFFFFFFC00000), /* 2^31 - 1 as a double */
                                 UINT64_C(0xC1E0000000000000), /* -2^31 as a double */
                                 UINT64_C(0x3FF0000000000000)}; /* 1 */

/* One damage applied to buffer[0 .. *size) in place; the buffer holds *size + MAX_GROWTH. */
static void damage(uint64_t* random, unsigned char* buffer, size_t* size) {
  const size_t n = *size;
  const size_t kind = below(random, 8);
  if (n == 0 && kind != 4) {
    buffer[0] = (unsigned char)SPECIAL[below(random, sizeof SPECIAL)];
    *size = 1;
    return;
  }
  const size_t at = below(random, n);
  switch (kind) {
    case 0: /* a meaningful character */
      buffer[at] = (unsigned char)SPECIAL[below(random, sizeof SPECIAL)];
      break;
    case 1: /* any byte */
      buffer[at] = (unsigned char)next_random(random);
      break;
    case 2: { /* a span deleted */
      const size_t length = 1 + below(random, n - at < 64 ? n - at : 64);
      memmove(buffer + at, buffer + at + length, n - at - length);
      *size = n - length;
      break;
    }
    case 3: { /* a span repeated */
      const size_t length = 1 + below(random, n - at < 64 ? n - at : 64);
      memmove(buffer + at + length, buffer + at, n - at);
      *size = n + length;
      break;
    }
    case 4: { /* a number inserted */
      const char* number = NUMBERS[below(random, sizeof NUMBERS / sizeof NUMBERS[0])];
      const size_t length = strlen(number);
      const size_t where = n == 0 ? 0 : below(random, n + 1);
      memmove(buffer + where + length, buffer + where, n - where);
      memcpy(buffer + where, number, length);
      *size = n + length;
      break;
    }
    case 5: /* a bit */
      buffer[at] ^= (unsigned char)(1u << below(random, 8));
      break;
    case 6: { /* an extreme word at an aligned place: 4 or 8 bytes */
      const uint64_t word = WORDS[below(random, sizeof WORDS / sizeof WORDS[0])];
      const size_t width = below(random, 2) == 0 ? 4 : 8;
      const size_t place = at & ~(width - 1);
      for (size_t i = 0; i < width && place + i < n; ++i) {
        buffer[place + i] = (unsigned char)(word >> (8 * i));
      }
      break;
    }
    default: /* eight digits inserted: a field wider than its reader expects */
      memmove(buffer + at + 8, buffer + at, n - at);
      memset(buffer + at, '9', 8);
      *size = n + 8;
      break;
  }
}

/* Builds variant k (>= 1) of input[0 .. size) into buffer (size + MAX_GROWTH bytes); returns
 * its size. Variants 1 .. TRUNCATIONS are truncations, the next ones damaged copies. */
static size_t make_variant(const unsigned char* input, size_t size, long k, uint64_t seed,
                           unsigned char* buffer) {
  if (k <= TRUNCATIONS) {
    size_t length = 0;
    if (k <= 32) {
      length = (size_t)(((uint64_t)size * (uint64_t)(k - 1)) / 32u);
    } else {
      const size_t back = (size_t)(k - 32);
      length = size > back ? size - back : 0;
    }
    memcpy(buffer, input, length);
    return length;
  }
  uint64_t random = seed ^ (UINT64_C(0xD1B54A32D192ED03) * (uint64_t)k);
  memcpy(buffer, input, size);
  size_t length = size;
  const size_t count = 1 + below(&random, 4);
  for (size_t i = 0; i < count; ++i) {
    damage(&random, buffer, &length);
  }
  return length;
}

/* ---- Running ------------------------------------------------------------------------------------ */

/* Copies [data, data + size) into a block of exactly that size and runs it. */
static void run_one(const unsigned char* data, size_t size) {
  unsigned char* copy = malloc(size == 0 ? 1 : size);
  MYB_CHECK(copy != NULL);
  memcpy(copy, data, size);
  (void)LLVMFuzzerTestOneInput(copy, size);
  free(copy);
}

static int save_variant(const char* path, const unsigned char* data, size_t size) {
  const myb_allocator a = myb_internal_allocator(NULL);
  const myb_status status = myb_platform_write_file(path, &a, data, size);
  if (status != MYB_OK) {
    fprintf(stderr, "cannot write %s (%s)\n", path, myb_status_name(status));
    return 1;
  }
  fprintf(stderr, "wrote %s (%zu bytes)\n", path, size);
  return 0;
}

/* Runs one file and its variants; returns 0, or 1 if it cannot be read. */
static int run_file(const char* path, const options* o, long* runs) {
  const myb_allocator a = myb_internal_allocator(NULL);
  char* text = NULL;
  size_t size = 0;
  const myb_status status = myb_platform_read_file(path, &a, MAX_INPUT, &text, &size);
  if (status != MYB_OK) {
    fprintf(stderr, "cannot read %s (%s)\n", path, myb_status_name(status));
    return 1;
  }
  const unsigned char* input = (const unsigned char*)text;
  unsigned char* buffer = malloc(size + MAX_GROWTH);
  MYB_CHECK(buffer != NULL);
  const long last = o->mutate > 0 ? TRUNCATIONS + o->mutate : 0;
  int result = 0;
  /* bound: the variants asked for */
  for (long k = 0; k <= last; ++k) {
    if (o->variant >= 0 && k != o->variant) {
      continue;
    }
    size_t length = size;
    const unsigned char* data = input;
    if (k > 0) {
      length = make_variant(input, size, k, o->seed, buffer);
      data = buffer;
    }
    if (o->verbose) {
      fprintf(stderr, "%s: variant %ld (%zu bytes)\n", path, k, length);
      fflush(stderr);
    }
    /* saved before it runs: a variant that ends the process is kept */
    if (o->save != NULL && (o->variant >= 0 || k == last)) {
      result |= save_variant(o->save, data, length);
    }
    run_one(data, length);
    ++*runs;
  }
  free(buffer);
  myb_internal_release(&a, text);
  return result;
}

/* Runs every path listed in a file, one per line (blank lines skipped). */
static int run_list(const char* path, const options* o, long* runs, long* files) {
  const myb_allocator a = myb_internal_allocator(NULL);
  char* text = NULL;
  size_t size = 0;
  const myb_status status = myb_platform_read_file(path, &a, MAX_INPUT, &text, &size);
  if (status != MYB_OK) {
    fprintf(stderr, "cannot read the list %s (%s)\n", path, myb_status_name(status));
    return 1;
  }
  int result = 0;
  size_t at = 0;
  /* bound: one line per iteration */
  while (at < size) {
    size_t end = at;
    /* bound: the rest of the text */
    while (end < size && text[end] != '\n') {
      ++end;
    }
    size_t stop = end;
    /* bound: the line */
    while (stop > at && (text[stop - 1] == '\r' || text[stop - 1] == ' ')) {
      --stop;
    }
    if (stop > at) {
      text[stop] = '\0';
      result |= run_file(text + at, o, runs);
      ++*files;
    }
    at = end + 1;
  }
  myb_internal_release(&a, text);
  return result;
}

static int numeric_option(const char* arg, const char* name, long* value) {
  const size_t n = strlen(name);
  if (strncmp(arg, name, n) != 0) {
    return 0;
  }
  char* end = NULL;
  *value = strtol(arg + n, &end, 10);
  return end != arg + n && *end == '\0';
}

int main(int argc, char** argv) {
  (void)LLVMFuzzerInitialize(&argc, &argv);
  options o = {.mutate = 0, .seed = 1, .variant = -1, .save = NULL, .verbose = 0};
  long runs = 0;
  long files = 0;
  int result = 0;
  if (argc > MAX_ARGUMENT_FILES) {
    fprintf(stderr, "too many arguments\n");
    return 1;
  }
  /* bound: the arguments */
  for (int i = 1; i < argc; ++i) {
    const char* arg = argv[i];
    long value = 0;
    if (numeric_option(arg, "-mutate=", &value) && value >= 0) {
      o.mutate = value;
    } else if (numeric_option(arg, "-seed=", &value)) {
      o.seed = (uint64_t)value;
    } else if (numeric_option(arg, "-variant=", &value) && value >= 0) {
      o.variant = value;
    } else if (strncmp(arg, "-save=", 6) == 0) {
      o.save = arg + 6;
    } else if (strcmp(arg, "-verbose") == 0) {
      o.verbose = 1;
    } else if (arg[0] == '-') {
      fprintf(stderr, "unknown option %s\n", arg);
      return 1;
    }
  }
  /* bound: the arguments */
  for (int i = 1; i < argc; ++i) {
    const char* arg = argv[i];
    if (arg[0] == '-') {
      continue;
    }
    if (arg[0] == '@') {
      result |= run_list(arg + 1, &o, &runs, &files);
    } else {
      result |= run_file(arg, &o, &runs);
      ++files;
    }
  }
  printf("%s: %ld inputs from %ld files\n", argc > 0 ? argv[0] : "replay", runs, files);
  return result;
}
