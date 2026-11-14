/* SPDX-License-Identifier: BSD-3-Clause */
/* site_names.h — the names of points a file of loads may hold, for the harnesses of the
 * readers that look a point up by its name: of the first lines that are not comments ($$),
 * the first word and the whole line without the blanks around it. */
#ifndef MYB_TESTS_FUZZ_SITE_NAMES_H_
#define MYB_TESTS_FUZZ_SITE_NAMES_H_

#include "fuzz.h"

enum { SITE_NAMES = 8, SITE_NAME_CHARS = 48 };

/* Fills names with up to SITE_NAMES candidates; returns how many. */
static inline int site_names(const uint8_t* data, size_t size,
                             char names[SITE_NAMES][SITE_NAME_CHARS]) {
  int count = 0;
  char line[SITE_NAME_CHARS];
  /* bound: at most 64 lines, SITE_NAMES names */
  for (size_t k = 0; k < 64 && count + 2 <= SITE_NAMES && fuzz_line(data, size, k, line, sizeof line);
       ++k) {
    if (line[0] == '$' && line[1] == '$') {
      continue;
    }
    size_t start = 0;
    /* bound: the line */
    while (line[start] == ' ' || line[start] == '\t') {
      ++start;
    }
    size_t end = strlen(line);
    /* bound: the line */
    while (end > start && (line[end - 1] == ' ' || line[end - 1] == '\t')) {
      --end;
    }
    if (end == start) {
      continue;
    }
    size_t word = start;
    /* bound: the line */
    while (word < end && line[word] != ' ' && line[word] != '\t') {
      ++word;
    }
    memcpy(names[count], line + start, word - start);
    names[count][word - start] = '\0';
    ++count;
    if (word < end) {
      memcpy(names[count], line + start, end - start);
      names[count][end - start] = '\0';
      ++count;
    }
  }
  return count;
}

#endif /* MYB_TESTS_FUZZ_SITE_NAMES_H_ */
