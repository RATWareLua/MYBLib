/* SPDX-License-Identifier: BSD-3-Clause */
/* io/daf.h — NAIF Double precision Array Files: the container of binary SPK and PCK kernels.
 *
 * Format ("DAF Required Reading", NAIF, JPL): 1024-byte records; record 1 is the file record
 * (identification word, ND and NI, the first and last summary records, the byte order); a
 * doubly linked list of summary records holds the descriptors ("summaries") of the arrays, ND
 * doubles and NI 32-bit integers each, packed into ND + (NI + 1) / 2 double words, and the
 * record after each summary record holds the arrays' names. Arrays are addressed in double
 * words from 1. Files of either IEEE byte order are read on any host: every value is decoded
 * from its bytes, never reinterpreted in place. */
#ifndef MYB_SRC_IO_DAF_H_
#define MYB_SRC_IO_DAF_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"

enum {
  MYB_DAF_RECORD = 1024,     /* bytes */
  MYB_DAF_MAX_ND = 124,      /* limits of the format */
  MYB_DAF_MAX_NI = 250,
  MYB_DAF_MAX_SUMMARY = 125, /* double words in a summary */
  MYB_DAF_MAX_NAME = 1000    /* characters in an array name */
};

typedef struct myb_internal_daf {
  const unsigned char* data;
  uint64_t words;  /* whole double words in the file */
  int big_endian;  /* byte order of the file */
  int nd;          /* doubles per summary */
  int ni;          /* integers per summary */
  char kind[5];    /* "SPK", "PCK", ... from the identification word "DAF/xxxx"; "" if old */
  uint64_t first;  /* first summary record (1-based record number) */
  uint64_t records; /* whole records in the file */
} myb_internal_daf;

/* Validates the file record of data[0 .. size). The memory must outlive the handle.
 * Errors: MYB_E_FORMAT (not a DAF, VAX formats, damaged by a text-mode transfer, ...). */
[[nodiscard]] myb_status myb_internal_daf_open(const unsigned char* data, size_t size,
                                               myb_internal_daf* daf, myb_diag* diag);

/* The double at a 1-based word address, 1 <= address <= daf->words (checked by the caller). */
double myb_internal_daf_double(const myb_internal_daf* daf, uint64_t address);

/* A position in the list of summaries. Zero-initialise before the first call. */
typedef struct myb_internal_daf_cursor {
  uint64_t record;  /* current summary record, 0 before the start */
  uint64_t visited; /* summary records seen (cycle guard) */
  int index;        /* next summary in the record */
  int count;        /* summaries in the record */
} myb_internal_daf_cursor;

/* One summary: nd doubles, ni integers and the name (NUL-terminated, trailing blanks removed). */
typedef struct myb_internal_daf_summary {
  double dc[MYB_DAF_MAX_ND];
  int32_t ic[MYB_DAF_MAX_NI];
  char name[MYB_DAF_MAX_NAME + 1];
} myb_internal_daf_summary;

/* Next summary in file order: MYB_OK with *found = 1, or *found = 0 at the end.
 * Errors: MYB_E_FORMAT (a broken or cyclic list). */
[[nodiscard]] myb_status myb_internal_daf_next(const myb_internal_daf* daf,
                                               myb_internal_daf_cursor* cursor,
                                               myb_internal_daf_summary* summary, int* found,
                                               myb_diag* diag);

#endif /* MYB_SRC_IO_DAF_H_ */
