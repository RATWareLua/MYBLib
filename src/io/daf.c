/* SPDX-License-Identifier: BSD-3-Clause */
#include "io/daf.h"

#include <string.h>

#include "core/alloc.h"
#include "core/base.h"

/* The FTP validation string: a text-mode transfer alters its CR, LF and 8-bit characters. */
static const unsigned char ftp_string[28] = {'F', 'T', 'P', 'S', 'T', 'R', ':', '\r', ':', '\n',
                                             ':', '\r', '\n', ':', '\r', 0x00, ':', 0x81, ':',
                                             0x10, 0xce, ':', 'E', 'N', 'D', 'F', 'T', 'P'};

enum {
  FILE_ND = 8,
  FILE_NI = 12,
  FILE_FWARD = 76,
  FILE_FORMAT = 88,
  FILE_FTP = 699,
  RECORD_WORDS = MYB_DAF_RECORD / 8
};

static uint64_t load64(const unsigned char* p, int big_endian) {
  uint64_t u = 0;
  for (int i = 0; i < 8; ++i) {
    const int k = big_endian ? i : 7 - i;
    u = (u << 8) | p[k];
  }
  return u;
}

static int32_t load32(const unsigned char* p, int big_endian) {
  uint32_t u = 0;
  for (int i = 0; i < 4; ++i) {
    const int k = big_endian ? i : 3 - i;
    u = (u << 8) | p[k];
  }
  int32_t s;
  memcpy(&s, &u, sizeof s);
  return s;
}

static double load_double(const unsigned char* p, int big_endian) {
  const uint64_t u = load64(p, big_endian);
  double d;
  memcpy(&d, &u, sizeof d);
  return d;
}

static int valid_counts(int32_t nd, int32_t ni) {
  return nd >= 0 && nd <= MYB_DAF_MAX_ND && ni >= 2 && ni <= MYB_DAF_MAX_NI &&
         nd + (ni + 1) / 2 <= MYB_DAF_MAX_SUMMARY;
}

static int all_zero(const unsigned char* p, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    if (p[i] != 0) {
      return 0;
    }
  }
  return 1;
}

static int all_blank(const unsigned char* p, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    if (p[i] != ' ' && p[i] != 0) {
      return 0;
    }
  }
  return 1;
}

myb_status myb_internal_daf_open(const unsigned char* data, size_t size, myb_internal_daf* daf,
                                 myb_diag* diag) {
  MYB_CHECK(daf != NULL && (data != NULL || size == 0));
  myb_internal_diag_clear(diag);
  if (size < 2 * (size_t)MYB_DAF_RECORD) {
    myb_internal_diag(diag, 0, "DAF: shorter than two records");
    return MYB_E_FORMAT;
  }
  char kind[5] = {0};
  if (memcmp(data, "DAF/", 4) == 0) {
    size_t n = 0;
    /* bound: 4 characters */
    for (size_t i = 4; i < 8 && data[i] != ' ' && data[i] != 0; ++i) {
      kind[n++] = (char)data[i];
    }
  } else if (memcmp(data, "NAIF/DAF", 8) != 0) {
    myb_internal_diag(diag, 0, "DAF: unknown identification word");
    return MYB_E_FORMAT;
  }
  int big_endian = -1;
  if (memcmp(data + FILE_FORMAT, "LTL-IEEE", 8) == 0) {
    big_endian = 0;
  } else if (memcmp(data + FILE_FORMAT, "BIG-IEEE", 8) == 0) {
    big_endian = 1;
  } else if (all_blank(data + FILE_FORMAT, 8)) {
    /* Files older than the format word: the byte order that gives valid ND and NI. */
    if (valid_counts(load32(data + FILE_ND, 0), load32(data + FILE_NI, 0))) {
      big_endian = 0;
    } else if (valid_counts(load32(data + FILE_ND, 1), load32(data + FILE_NI, 1))) {
      big_endian = 1;
    }
  }
  if (big_endian < 0) {
    myb_internal_diag(diag, 0, "DAF: unsupported binary format (only IEEE files are read)");
    return MYB_E_FORMAT;
  }
  if (!all_zero(data + FILE_FTP, sizeof ftp_string) &&
      memcmp(data + FILE_FTP, ftp_string, sizeof ftp_string) != 0) {
    myb_internal_diag(diag, 0, "DAF: damaged by a text-mode transfer (FTP string altered)");
    return MYB_E_FORMAT;
  }
  const int32_t nd = load32(data + FILE_ND, big_endian);
  const int32_t ni = load32(data + FILE_NI, big_endian);
  if (!valid_counts(nd, ni)) {
    myb_internal_diag(diag, 0, "DAF: invalid ND = %ld, NI = %ld", (long)nd, (long)ni);
    return MYB_E_FORMAT;
  }
  const int32_t fward = load32(data + FILE_FWARD, big_endian);
  const uint64_t records = (uint64_t)size / MYB_DAF_RECORD;
  if (fward < 2 || (uint64_t)fward + 1 > records) {
    myb_internal_diag(diag, 0, "DAF: first summary record %ld outside the file", (long)fward);
    return MYB_E_FORMAT;
  }
  daf->data = data;
  daf->words = (uint64_t)size / 8;
  daf->big_endian = big_endian;
  daf->nd = nd;
  daf->ni = ni;
  memcpy(daf->kind, kind, sizeof kind);
  daf->first = (uint64_t)fward;
  daf->records = records;
  return MYB_OK;
}

double myb_internal_daf_double(const myb_internal_daf* daf, uint64_t address) {
  MYB_DCHECK(address >= 1 && address <= daf->words);
  return load_double(daf->data + (address - 1) * 8, daf->big_endian);
}

/* Reads the control words of summary record `record`; MYB_E_FORMAT if they are invalid. */
static myb_status open_record(const myb_internal_daf* daf, uint64_t record, int* count,
                              myb_diag* diag) {
  if (record < 2 || record + 1 > daf->records) {
    myb_internal_diag(diag, 0, "DAF: summary record %llu outside the file",
                      (unsigned long long)record);
    return MYB_E_FORMAT;
  }
  const unsigned char* p = daf->data + (record - 1) * MYB_DAF_RECORD;
  const double nsum = load_double(p + 16, daf->big_endian);
  const int size = daf->nd + (daf->ni + 1) / 2;
  const int capacity = (RECORD_WORDS - 3) / size;
  if (!(nsum >= 0.0 && nsum <= capacity) || nsum != (double)(int)nsum) {
    myb_internal_diag(diag, 0, "DAF: summary record %llu: invalid count",
                      (unsigned long long)record);
    return MYB_E_FORMAT;
  }
  *count = (int)nsum;
  return MYB_OK;
}

myb_status myb_internal_daf_next(const myb_internal_daf* daf, myb_internal_daf_cursor* cursor,
                                 myb_internal_daf_summary* summary, int* found, myb_diag* diag) {
  MYB_CHECK(daf != NULL && cursor != NULL && summary != NULL && found != NULL);
  *found = 0;
  if (cursor->record == 0) {
    const myb_status status = open_record(daf, daf->first, &cursor->count, diag);
    if (status != MYB_OK) {
      return status;
    }
    cursor->record = daf->first;
    cursor->visited = 1;
    cursor->index = 0;
  }
  /* bound: each iteration moves to the next record; at most daf->records records exist */
  while (cursor->index >= cursor->count) {
    const unsigned char* p = daf->data + (cursor->record - 1) * MYB_DAF_RECORD;
    const double next = load_double(p, daf->big_endian);
    if (next == 0.0) {
      return MYB_OK;
    }
    if (!(next >= 2.0 && next <= (double)daf->records) || next != (double)(uint64_t)next) {
      myb_internal_diag(diag, 0, "DAF: invalid link in summary record %llu",
                        (unsigned long long)cursor->record);
      return MYB_E_FORMAT;
    }
    if (++cursor->visited > daf->records) {
      myb_internal_diag(diag, 0, "DAF: the list of summary records is cyclic");
      return MYB_E_FORMAT;
    }
    const uint64_t record = (uint64_t)next;
    const myb_status status = open_record(daf, record, &cursor->count, diag);
    if (status != MYB_OK) {
      return status;
    }
    cursor->record = record;
    cursor->index = 0;
  }
  const int size = daf->nd + (daf->ni + 1) / 2;
  const unsigned char* record = daf->data + (cursor->record - 1) * MYB_DAF_RECORD;
  const unsigned char* p = record + 24 + (size_t)cursor->index * (size_t)size * 8;
  for (int i = 0; i < daf->nd; ++i) {
    summary->dc[i] = load_double(p + 8 * i, daf->big_endian);
  }
  for (int i = 0; i < daf->ni; ++i) {
    summary->ic[i] = load32(p + 8 * daf->nd + 4 * i, daf->big_endian);
  }
  const size_t name_size = (size_t)size * 8;
  const unsigned char* name = record + MYB_DAF_RECORD + (size_t)cursor->index * name_size;
  size_t length = name_size;
  /* bound: at most name_size characters */
  while (length > 0 && (name[length - 1] == ' ' || name[length - 1] == 0)) {
    --length;
  }
  memcpy(summary->name, name, length);
  summary->name[length] = '\0';
  ++cursor->index;
  *found = 1;
  return MYB_OK;
}
