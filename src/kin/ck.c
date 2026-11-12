/* SPDX-License-Identifier: BSD-3-Clause */
/* Reading NAIF CK files (myb/ck.h): segments of types 1 to 6 in a DAF, after "CK Required
 * Reading" and "Generic Segments Required Reading" (NAIF, JPL). Which record answers a
 * request — inside an interval, at its ends, within the tolerance — and how the pointing is
 * carried between records were checked against ckgpav_c of CSPICE N0067 (tests/data/ck); the
 * rules of selection of types 4, 5 and 6 follow the descriptions of their readers in the
 * CSPICE sources (ckr04, ckr05, ckr06), which serve as the description of the format and are
 * not copied.
 *
 * Layouts (words of a segment, N pointings, R = 7 with angular velocities and 4 without):
 *   type 1: N records of R; N ticks; (N - 1) / 100 directory ticks; N
 *   type 2: N records of 8 (quaternion, angular velocity, seconds per tick); N starts;
 *           N stops; (N - 1) / 100 directory ticks
 *   type 3: N records of R; N ticks; (N - 1) / 100 directory ticks; M interval starts;
 *           (M - 1) / 100 directory ticks; M; N
 *   type 4: a generic segment of packets of variable size (mid, radius, the counts of
 *           coefficients of the seven components packed in one number, the coefficients),
 *           each with its start as the reference value; the metadata at the end
 *   type 5: N packets of P (by the subtype); N ticks; (N - 1) / 100 directory ticks;
 *           M interval starts; (M - 1) / 100 directory ticks; seconds per tick, subtype,
 *           window, M, N
 *   type 6: K mini-segments (each: packets, ticks, directory, seconds per tick, subtype,
 *           window, count); K + 1 boundaries; K / 100 directory boundaries; K + 1 pointers
 *           to the mini-segments (1-based offsets); the selection flag; K */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "io/daf.h"
#include "kin/spk_segment.h"
#include "myb/ck.h"
#include "myb/cmath.h"
#include "num/vec3.h"
#include "platform/platform.h"

typedef struct ck_segment {
  myb_internal_daf daf;
  int32_t instrument;
  int32_t frame;
  int32_t type;
  int has_rates;
  double start;
  double stop;
  uint64_t begin;
  uint64_t end;
  size_t file;
  char name[41];
  int supported;
  uint64_t count;     /* pointings or intervals */
  uint64_t record;    /* words in a record */
  uint64_t ticks;     /* address of the first tick (type 2: of the first start) */
  uint64_t stops;     /* type 2: address of the first stop */
  uint64_t intervals; /* types 3 and 5: how many, and the address of the first start */
  uint64_t starts;
  /* type 4: the partitions of the generic segment (addresses before their first items) */
  uint64_t packet_base;
  uint64_t directory_base;
  uint64_t reference_base;
  uint64_t packet_offset;
  /* type 5 */
  int subtype;
  int window;
  double clock_rate; /* seconds per tick */
  /* type 6: count mini-segments, the address of the first of count + 1 boundaries and of
   * count + 1 pointers */
  uint64_t bounds;
  uint64_t pointers;
  int select_last;
} ck_segment;

/* The packets of the subtypes of types 5 and 6 and their largest windows. */
static const uint64_t packet_size[4] = {8, 4, 14, 7};
static const uint64_t largest_window[4] = {12, 24, 12, 24};

struct myb_ck {
  myb_allocator allocator;
  myb_platform_mapping* files;
  size_t file_count;
  ck_segment* segments; /* load order */
  size_t count;
};

static double word(const ck_segment* s, uint64_t address) {
  return myb_internal_daf_double(&s->daf, address);
}

static myb_status format_error(myb_diag* diag, const ck_segment* s, const char* what) {
  myb_internal_diag(diag, 0, "segment '%s' (type %d): %s", s->name, (int)s->type, what);
  return MYB_E_FORMAT;
}

/* A double that holds a count below `limit`. */
static int as_count(double value, uint64_t limit, uint64_t* out) {
  if (!(value >= 0.0) || !(value <= (double)limit) || value != floor(value)) {
    return 0;
  }
  *out = (uint64_t)value;
  return 1;
}

static int increasing(const ck_segment* s, uint64_t first, uint64_t n) {
  for (uint64_t i = 0; i < n; ++i) {
    const double v = word(s, first + i);
    if (!isfinite(v) || (i > 0 && !(v > word(s, first + i - 1)))) {
      return 0;
    }
  }
  return 1;
}

/* Every value of b at or after a, each b found among the n increasing values from `first`
 * (a merge of two increasing sequences). */
static int among(const ck_segment* s, uint64_t first, uint64_t n, uint64_t from, uint64_t m) {
  uint64_t at = 0;
  for (uint64_t i = 0; i < m; ++i) {
    const double v = word(s, from + i);
    /* bound: the n values, each passed once over all the m */
    while (at < n && word(s, first + at) < v) {
      ++at;
    }
    if (at == n || word(s, first + at) != v) {
      return 0;
    }
  }
  return 1;
}

/* Type 4: the metadata of the generic segment and every packet. */
static myb_status parse_type4(ck_segment* s, myb_diag* diag) {
  const uint64_t length = s->end - s->begin + 1;
  uint64_t stored = 0;
  if (length < 17 || !as_count(word(s, s->end), length, &stored) || stored < 15) {
    return format_error(diag, s, "generic segment metadata");
  }
  /* 15 items: an old segment that does not count the count itself; more than 17: items
   * after those known */
  uint64_t known = stored;
  uint64_t span = stored;
  if (stored == 15) {
    known = 16;
    span = 16;
  } else if (stored > 17) {
    known = 17;
  }
  uint64_t meta[17] = {0};
  const uint64_t first = s->end - span + 1;
  for (uint64_t i = 0; i < known; ++i) {
    if (!as_count(word(s, first + i), 2 * length, &meta[i])) {
      return format_error(diag, s, "generic segment metadata");
    }
  }
  meta[16] = known;
  for (uint64_t i = known - 1; i < 16; ++i) {
    meta[i] = 0;
  }
  const uint64_t n = meta[11];
  const uint64_t base = s->begin - 1;
  /* variable packets (packet directory type 1), explicit references "at or before" (3) */
  if (meta[9] != 1 || meta[4] != 3 || n == 0 || meta[6] != n || meta[8] != n + 1 ||
      meta[5] + n > length || meta[7] + n + 1 > length || meta[10] >= length) {
    return format_error(diag, s, "not a segment of variable packets with references");
  }
  s->count = n;
  s->packet_base = base + meta[10];
  s->directory_base = base + meta[7];
  s->reference_base = base + meta[5];
  s->packet_offset = meta[15];
  if (!increasing(s, s->reference_base + 1, n)) {
    return format_error(diag, s, "references not increasing");
  }
  for (uint64_t k = 1; k <= n; ++k) {
    uint64_t b1 = 0;
    uint64_t b2 = 0;
    uint64_t code = 0;
    if (!as_count(word(s, s->directory_base + k), length, &b1) ||
        !as_count(word(s, s->directory_base + k + 1), length, &b2) ||
        b2 < b1 + s->packet_offset + 7) {
      return format_error(diag, s, "packet directory");
    }
    const uint64_t size = b2 - b1 - s->packet_offset;
    const uint64_t at = s->packet_base + b1;
    if (at + size - 1 > s->end ||
        !as_count(word(s, at + 2), (uint64_t)1 << 50, &code)) {
      return format_error(diag, s, "packet outside the segment");
    }
    uint64_t total = 0;
    for (int i = 0; i < 7; ++i) {
      const uint64_t c = code % 128;
      code /= 128;
      if (c > 19 || (i < 4 && c == 0) || (i >= 4 && s->has_rates && c == 0)) {
        return format_error(diag, s, "counts of coefficients");
      }
      total += c;
    }
    const double mid = word(s, at);
    const double radius = word(s, at + 1);
    if (code != 0 || total + 3 != size || !isfinite(mid) || !(radius > 0.0) ||
        !isfinite(radius)) {
      return format_error(diag, s, "packet");
    }
  }
  return MYB_OK;
}

/* The trailer of a set of packets of type 5 or of a mini-segment of type 6: seconds per
 * tick, subtype, window. */
static int window_parameters(const ck_segment* s, uint64_t at, double* rate, int* subtype,
                             int* window) {
  uint64_t sub_type = 0;
  uint64_t size = 0;
  *rate = word(s, at);
  if (!(*rate > 0.0) || !isfinite(*rate) || !as_count(word(s, at + 1), 3, &sub_type) ||
      !as_count(word(s, at + 2), MYB_CK_MAX_WINDOW, &size) || size < 2 || size % 2 != 0 ||
      size > largest_window[sub_type]) {
    return 0;
  }
  *subtype = (int)sub_type;
  *window = (int)size;
  return 1;
}

static myb_status parse_type5(ck_segment* s, myb_diag* diag) {
  const uint64_t length = s->end - s->begin + 1;
  uint64_t n = 0;
  uint64_t m = 0;
  if (length < 5 || !window_parameters(s, s->end - 4, &s->clock_rate, &s->subtype, &s->window) ||
      !as_count(word(s, s->end - 1), length, &m) || !as_count(word(s, s->end), length, &n) ||
      n == 0 || m == 0 || m > n) {
    return format_error(diag, s, "inconsistent trailer");
  }
  const uint64_t size = packet_size[s->subtype];
  if (n * (size + 1) + (n - 1) / 100 + m + (m - 1) / 100 + 5 != length) {
    return format_error(diag, s, "inconsistent length");
  }
  s->count = n;
  s->record = size;
  s->ticks = s->begin + n * size;
  s->intervals = m;
  s->starts = s->ticks + n + (n - 1) / 100;
  if (!increasing(s, s->ticks, n) || !increasing(s, s->starts, m) ||
      word(s, s->starts) != word(s, s->ticks) || !among(s, s->ticks, n, s->starts, m)) {
    return format_error(diag, s, "ticks or interval starts");
  }
  return MYB_OK;
}

/* A mini-segment of type 6: its first word and its last. */
static void mini_bounds(const ck_segment* s, uint64_t k, uint64_t* first, uint64_t* last) {
  *first = s->begin + (uint64_t)word(s, s->pointers + k) - 1;
  *last = s->begin + (uint64_t)word(s, s->pointers + k + 1) - 2;
}

static myb_status parse_type6(ck_segment* s, myb_diag* diag) {
  const uint64_t length = s->end - s->begin + 1;
  uint64_t n = 0;
  uint64_t select = 0;
  if (length < 2 || !as_count(word(s, s->end - 1), 1, &select) ||
      !as_count(word(s, s->end), length / 2, &n) || n == 0 ||
      2 * (n + 1) + n / 100 + 2 > length) {
    return format_error(diag, s, "inconsistent trailer");
  }
  s->count = n;
  s->select_last = (int)select;
  s->pointers = s->end - 2 - n;
  s->bounds = s->pointers - n / 100 - (n + 1);
  if (!increasing(s, s->bounds, n + 1)) {
    return format_error(diag, s, "boundaries not increasing");
  }
  uint64_t previous = 0;
  for (uint64_t k = 0; k <= n; ++k) {
    uint64_t p = 0;
    if (!as_count(word(s, s->pointers + k), length, &p) || (k == 0 && p != 1) ||
        (k > 0 && p < previous + 4) || (k == n && s->begin + p - 1 != s->bounds)) {
      return format_error(diag, s, "mini-segment pointers");
    }
    previous = p;
  }
  for (uint64_t k = 0; k < n; ++k) {
    uint64_t first = 0;
    uint64_t last = 0;
    mini_bounds(s, k, &first, &last);
    double rate = 0.0;
    int sub_type = 0;
    int window = 0;
    uint64_t count = 0;
    if (!window_parameters(s, last - 3, &rate, &sub_type, &window) ||
        !as_count(word(s, last), last - first + 1, &count) || count < 2) {
      return format_error(diag, s, "mini-segment trailer");
    }
    const uint64_t size = packet_size[sub_type];
    const uint64_t ticks = first + count * size;
    if (count * (size + 1) + (count - 1) / 100 + 4 != last - first + 1 ||
        !increasing(s, ticks, count) || !(word(s, ticks) <= word(s, s->bounds + k)) ||
        !(word(s, ticks + count - 1) >= word(s, s->bounds + k))) {
      return format_error(diag, s, "mini-segment");
    }
  }
  return MYB_OK;
}

/* The layout of a segment; `supported` is cleared for the types that are not read. */
static myb_status parse(ck_segment* s, myb_diag* diag) {
  const uint64_t length = s->end - s->begin + 1;
  const uint64_t r = s->has_rates ? 7 : 4;
  s->supported = 1;
  if (s->type == 1) {
    uint64_t n = 0;
    if (!as_count(word(s, s->end), length, &n) || n == 0 ||
        n * (r + 1) + (n - 1) / 100 + 1 != length) {
      return format_error(diag, s, "inconsistent directory");
    }
    s->count = n;
    s->record = r;
    s->ticks = s->begin + n * r;
  } else if (s->type == 2) {
    /* 10 N + (N - 1) / 100 words: N is the least whole number at or above
     * (100 length + 1) / 1001 */
    const uint64_t n = (100 * length + 1001) / 1001;
    if (n == 0 || 10 * n + (n - 1) / 100 != length) {
      return format_error(diag, s, "inconsistent length");
    }
    s->count = n;
    s->record = 8;
    s->ticks = s->begin + 8 * n;
    s->stops = s->ticks + n;
    for (uint64_t i = 0; i < n; ++i) {
      if (!(word(s, s->stops + i) >= word(s, s->ticks + i)) ||
          (i > 0 && !(word(s, s->ticks + i) >= word(s, s->stops + i - 1)))) {
        return format_error(diag, s, "intervals out of order");
      }
    }
  } else if (s->type == 3) {
    uint64_t n = 0;
    uint64_t m = 0;
    if (length < 2 || !as_count(word(s, s->end), length, &n) ||
        !as_count(word(s, s->end - 1), length, &m) || n == 0 || m == 0 || m > n ||
        n * (r + 1) + (n - 1) / 100 + m + (m - 1) / 100 + 2 != length) {
      return format_error(diag, s, "inconsistent directory");
    }
    s->count = n;
    s->record = r;
    s->ticks = s->begin + n * r;
    s->intervals = m;
    s->starts = s->ticks + n + (n - 1) / 100;
    if (!increasing(s, s->starts, m)) {
      return format_error(diag, s, "interval starts not increasing");
    }
  } else if (s->type == 4) {
    return parse_type4(s, diag);
  } else if (s->type == 5) {
    return parse_type5(s, diag);
  } else if (s->type == 6) {
    return parse_type6(s, diag);
  } else {
    s->supported = 0;
    return MYB_OK;
  }
  if (s->type != 2 && !increasing(s, s->ticks, s->count)) {
    return format_error(diag, s, "ticks not increasing");
  }
  return MYB_OK;
}

myb_status myb_ck_create(const myb_allocator* allocator, myb_ck** ck) {
  MYB_CHECK(ck != NULL);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_ck), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_ck* c = block;
  memset(c, 0, sizeof *c);
  c->allocator = a;
  *ck = c;
  return MYB_OK;
}

void myb_ck_destroy(myb_ck* ck) {
  if (ck == NULL) {
    return;
  }
  const myb_allocator a = ck->allocator;
  for (size_t i = 0; i < ck->file_count; ++i) {
    myb_platform_unmap_file(&ck->files[i]);
  }
  myb_internal_release(&a, ck->files);
  myb_internal_release(&a, ck->segments);
  myb_internal_release(&a, ck);
}

/* Appends the segments of a DAF image; on error the set is unchanged. */
static myb_status load(myb_ck* ck, const unsigned char* data, size_t size,
                       myb_platform_mapping mapping, myb_diag* diag) {
  myb_internal_diag_clear(diag);
  myb_internal_daf daf;
  myb_status status = myb_internal_daf_open(data, size, &daf, diag);
  if (status != MYB_OK) {
    return status;
  }
  if ((daf.kind[0] != '\0' && strcmp(daf.kind, "CK") != 0) || daf.nd != 2 || daf.ni != 6) {
    myb_internal_diag(diag, 0, "not a CK file (DAF/%s, ND = %d, NI = %d)", daf.kind, daf.nd,
                      daf.ni);
    return MYB_E_FORMAT;
  }
  const myb_allocator* a = &ck->allocator;
  ck_segment* segments = NULL;
  size_t count = ck->count;
  size_t capacity = 0;
  myb_internal_daf_cursor cursor = {0};
  myb_internal_daf_summary summary;
  /* bound: a summary occupies 5 words of the file, so at most words / 5 summaries */
  for (uint64_t guard = 0; guard <= daf.words / 5; ++guard) {
    int found = 0;
    status = myb_internal_daf_next(&daf, &cursor, &summary, &found, diag);
    if (status != MYB_OK || !found) {
      break;
    }
    if (segments == NULL || count == capacity) {
      const size_t grown = capacity == 0 ? ck->count + 16 : 2 * capacity;
      void* block = NULL;
      status = myb_internal_alloc_array(a, grown, sizeof *segments, &block);
      if (status != MYB_OK) {
        break;
      }
      if (segments != NULL) {
        memcpy(block, segments, count * sizeof *segments);
      } else if (ck->count > 0) {
        memcpy(block, ck->segments, ck->count * sizeof *segments);
      }
      myb_internal_release(a, segments);
      segments = block;
      capacity = grown;
    }
    ck_segment* s = &segments[count];
    memset(s, 0, sizeof *s);
    s->daf = daf;
    s->start = summary.dc[0];
    s->stop = summary.dc[1];
    s->instrument = summary.ic[0];
    s->frame = summary.ic[1];
    s->type = summary.ic[2];
    s->has_rates = summary.ic[3] != 0;
    s->file = ck->file_count;
    const size_t name_length = strlen(summary.name);
    memcpy(s->name, summary.name, name_length < 40 ? name_length : 40);
    const int32_t begin = summary.ic[4];
    const int32_t end = summary.ic[5];
    if (!(s->start <= s->stop) || !isfinite(s->start) || !isfinite(s->stop) || begin < 1 ||
        end < begin || (uint64_t)end > daf.words) {
      myb_internal_diag(diag, 0, "segment '%s': invalid descriptor", s->name);
      status = MYB_E_FORMAT;
      break;
    }
    s->begin = (uint64_t)begin;
    s->end = (uint64_t)end;
    status = parse(s, diag);
    if (status != MYB_OK) {
      break;
    }
    ++count;
  }
  myb_platform_mapping* files = NULL;
  if (status == MYB_OK) {
    void* block = NULL;
    status = myb_internal_alloc_array(a, ck->file_count + 1, sizeof *files, &block);
    files = block;
  }
  if (status != MYB_OK) {
    myb_internal_release(a, segments);
    return status;
  }
  if (ck->file_count > 0) {
    memcpy(files, ck->files, ck->file_count * sizeof *files);
  }
  files[ck->file_count] = mapping;
  myb_internal_release(a, ck->files);
  ck->files = files;
  ck->file_count += 1;
  if (segments != NULL) {
    myb_internal_release(a, ck->segments);
    ck->segments = segments;
    ck->count = count;
  }
  return MYB_OK;
}

myb_status myb_ck_load_memory(myb_ck* ck, const void* data, size_t size, myb_diag* diag) {
  MYB_CHECK(ck != NULL && (data != NULL || size == 0));
  const myb_platform_mapping none = {0};
  return load(ck, data, size, none, diag);
}

myb_status myb_ck_load_file(myb_ck* ck, const char* path, myb_diag* diag) {
  MYB_CHECK(ck != NULL && path != NULL);
  myb_internal_diag_clear(diag);
  myb_platform_mapping mapping = {0};
  myb_status status = myb_platform_map_file(path, &ck->allocator, &mapping);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  status = load(ck, mapping.data, mapping.size, mapping, diag);
  if (status != MYB_OK) {
    myb_platform_unmap_file(&mapping);
  }
  return status;
}

size_t myb_ck_segment_count(const myb_ck* ck) {
  MYB_CHECK(ck != NULL);
  return ck->count;
}

void myb_ck_segment(const myb_ck* ck, size_t index, myb_ck_segment_info* info) {
  MYB_CHECK(ck != NULL && info != NULL && index < ck->count);
  const ck_segment* s = &ck->segments[index];
  info->instrument = s->instrument;
  info->frame = s->frame;
  info->type = s->type;
  info->has_rates = s->has_rates;
  info->start = s->start;
  info->stop = s->stop;
  info->file = s->file;
  memcpy(info->name, s->name, sizeof info->name);
}

/* ---- Quaternions (SPICE's: the scalar first, the product in the order of the matrices) ----- */

static void q_multiply(const double a[4], const double b[4], double out[4]) {
  const double r[4] = {a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
                       a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
                       a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
                       a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]};
  memcpy(out, r, sizeof r);
}

static int q_unit(double q[4]) {
  const double n = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  if (!(n > 0.0) || !isfinite(n)) {
    return 0;
  }
  for (int i = 0; i < 4; ++i) {
    q[i] /= n;
  }
  return 1;
}

/* The quaternion of the turn by `angle` about the unit axis. */
static void q_turn(const double axis[3], double angle, double q[4]) {
  const double sn = myb_sin(0.5 * angle);
  q[0] = myb_cos(0.5 * angle);
  q[1] = sn * axis[0];
  q[2] = sn * axis[1];
  q[3] = sn * axis[2];
}

/* The record k of a segment: the quaternion made unit, the angular velocity (zeros without
 * one) and, for type 2, the seconds a tick lasts. 0 for a record that is not finite. */
static int read_record(const ck_segment* s, uint64_t k, double q[4], double rate[3],
                       double* seconds_per_tick) {
  const uint64_t at = s->begin + k * s->record;
  for (uint64_t i = 0; i < 4; ++i) {
    q[i] = word(s, at + i);
  }
  rate[0] = rate[1] = rate[2] = 0.0;
  if (s->record >= 7) {
    for (uint64_t i = 0; i < 3; ++i) {
      rate[i] = word(s, at + 4 + i);
    }
  }
  if (seconds_per_tick != NULL) {
    *seconds_per_tick = s->record == 8 ? word(s, at + 7) : 0.0;
  }
  return q_unit(q) && myb_vec3_finite(rate);
}

/* How many of the n increasing values from `first` are at or below x. */
static uint64_t count_le(const ck_segment* s, uint64_t first, uint64_t n, double x) {
  uint64_t lo = 0;
  uint64_t hi = n;
  /* bound: the interval halves */
  while (lo < hi) {
    const uint64_t mid = lo + (hi - lo) / 2;
    if (word(s, first + mid) <= x) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

/* What a segment answers: found = 0 if nothing is within the tolerance. */
typedef struct answer {
  int found;
  double q[4];
  double rate[3];
  double ticks;
} answer;

/* The record k alone, at its own instant. */
static myb_status at_record(const ck_segment* s, uint64_t k, answer* a) {
  if (!read_record(s, k, a->q, a->rate, NULL)) {
    return MYB_E_FORMAT;
  }
  a->ticks = word(s, s->ticks + k);
  a->found = 1;
  return MYB_OK;
}

/* Type 1: the pointing nearest to the instant, if within the tolerance. */
static myb_status answer_type1(const ck_segment* s, double ticks, double tolerance, answer* a) {
  const uint64_t le = count_le(s, s->ticks, s->count, ticks);
  uint64_t k = 0;
  if (le == 0) {
    k = 0;
  } else if (le == s->count) {
    k = s->count - 1;
  } else {
    /* the nearer of the two neighbours; the later one when they are as near */
    k = ticks - word(s, s->ticks + le - 1) < word(s, s->ticks + le) - ticks ? le - 1 : le;
  }
  if (fabs(word(s, s->ticks + k) - ticks) > tolerance) {
    return MYB_OK;
  }
  return at_record(s, k, a);
}

/* Type 2: inside an interval the frame has turned about the angular velocity since the
 * start of the interval; outside, the nearest end of an interval within the tolerance. */
static myb_status answer_type2(const ck_segment* s, double ticks, double tolerance, answer* a) {
  const uint64_t le = count_le(s, s->ticks, s->count, ticks);
  uint64_t k = 0;
  double at = ticks;
  if (le > 0 && ticks <= word(s, s->stops + le - 1)) {
    k = le - 1; /* inside */
  } else {
    /* between the stop of interval le - 1 and the start of interval le */
    const double before = le > 0 ? ticks - word(s, s->stops + le - 1) : HUGE_VAL;
    const double after = le < s->count ? word(s, s->ticks + le) - ticks : HUGE_VAL;
    if (before < after) { /* as near: the later one */
      k = le - 1;
      at = word(s, s->stops + k);
      if (before > tolerance) {
        return MYB_OK;
      }
    } else {
      k = le;
      at = word(s, s->ticks + k);
      if (after > tolerance) {
        return MYB_OK;
      }
    }
  }
  double seconds_per_tick = 0.0;
  if (!read_record(s, k, a->q, a->rate, &seconds_per_tick) || !isfinite(seconds_per_tick)) {
    return MYB_E_FORMAT;
  }
  const double speed = myb_vec3_norm(a->rate);
  const double angle = speed * seconds_per_tick * (at - word(s, s->ticks + k));
  if (speed > 0.0 && angle != 0.0) {
    /* the axes of the instrument turn about the angular velocity: C = C0 R^T */
    const double axis[3] = {a->rate[0] / speed, a->rate[1] / speed, a->rate[2] / speed};
    double turn[4];
    q_turn(axis, -angle, turn);
    q_multiply(a->q, turn, a->q);
  }
  a->ticks = at;
  a->found = 1;
  return MYB_OK;
}

/* Type 3: between two pointings of one interval the frame turns at a constant rate from
 * the one to the other, and the angular velocity goes linearly; elsewhere the nearest
 * pointing within the tolerance. */
static myb_status answer_type3(const ck_segment* s, double ticks, double tolerance, answer* a) {
  const uint64_t le = count_le(s, s->ticks, s->count, ticks);
  if (le == 0) {
    return word(s, s->ticks) - ticks <= tolerance ? at_record(s, 0, a) : MYB_OK;
  }
  const double t1 = word(s, s->ticks + le - 1);
  if (t1 == ticks) {
    return at_record(s, le - 1, a);
  }
  if (le == s->count) {
    return ticks - t1 <= tolerance ? at_record(s, le - 1, a) : MYB_OK;
  }
  const double t2 = word(s, s->ticks + le);
  /* the two are in one interval unless an interval starts at the later one */
  const uint64_t starts = count_le(s, s->starts, s->intervals, t2);
  const int apart = starts > 0 && word(s, s->starts + starts - 1) == t2;
  if (apart) {
    const int earlier = ticks - t1 < t2 - ticks; /* as near: the later one */
    if ((earlier ? ticks - t1 : t2 - ticks) > tolerance) {
      return MYB_OK;
    }
    return at_record(s, earlier ? le - 1 : le, a);
  }
  double q1[4];
  double q2[4];
  double r1[3];
  double r2[3];
  if (!read_record(s, le - 1, q1, r1, NULL) || !read_record(s, le, q2, r2, NULL)) {
    return MYB_E_FORMAT;
  }
  const double part = (ticks - t1) / (t2 - t1);
  /* q1^-1 q2: the turn from the first to the second, the short way */
  const double inverse[4] = {q1[0], -q1[1], -q1[2], -q1[3]};
  double between[4];
  q_multiply(inverse, q2, between);
  if (between[0] < 0.0) {
    for (int i = 0; i < 4; ++i) {
      between[i] = -between[i];
    }
  }
  const double sine = myb_vec3_norm(between + 1);
  memcpy(a->q, q1, sizeof q1);
  if (sine > 0.0) {
    const double angle = 2.0 * myb_atan2(sine, between[0]);
    const double axis[3] = {between[1] / sine, between[2] / sine, between[3] / sine};
    double turn[4];
    q_turn(axis, part * angle, turn);
    q_multiply(q1, turn, a->q);
  }
  for (int i = 0; i < 3; ++i) {
    a->rate[i] = r1[i] + part * (r2[i] - r1[i]);
  }
  a->ticks = ticks;
  a->found = 1;
  return MYB_OK;
}

/* How many of the n increasing values from `first` are below x. */
static uint64_t count_lt(const ck_segment* s, uint64_t first, uint64_t n, double x) {
  uint64_t lo = 0;
  uint64_t hi = n;
  /* bound: the interval halves */
  while (lo < hi) {
    const uint64_t mid = lo + (hi - lo) / 2;
    if (word(s, first + mid) < x) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

/* ---- Type 4 ---------------------------------------------------------------------------------- */

/* The k-th packet (0-based): the address of its first word (the mid). */
static uint64_t packet4(const ck_segment* s, uint64_t k) {
  return s->packet_base + (uint64_t)word(s, s->directory_base + k + 1);
}

/* A Chebyshev series of n coefficients from `at` at s in [-1, 1], with the full weight of
 * the term T0 (the sum of chbval in SPICE). */
static double chebyshev_at(const ck_segment* s, uint64_t at, uint64_t n, double x) {
  const double x2 = 2.0 * x;
  double w0 = 0.0;
  double w1 = 0.0;
  for (uint64_t j = n; j > 1; --j) {
    const double w2 = w1;
    w1 = w0;
    w0 = word(s, at + j - 1) + (x2 * w1 - w2);
  }
  return (x * w0 - w1) + word(s, at);
}

/* The pointing of packet k at the instant `at`. */
static myb_status evaluate4(const ck_segment* s, uint64_t k, double at, answer* a) {
  const uint64_t p = packet4(s, k);
  const double mid = word(s, p);
  const double radius = word(s, p + 1);
  uint64_t code = (uint64_t)word(s, p + 2);
  uint64_t counts[7];
  for (int i = 0; i < 7; ++i) {
    counts[i] = code % 128;
    code /= 128;
  }
  const double x = (at - mid) / radius;
  uint64_t c = p + 3;
  for (int i = 0; i < 4; ++i) {
    a->q[i] = chebyshev_at(s, c, counts[i], x);
    c += counts[i];
  }
  for (int i = 0; i < 3; ++i) {
    a->rate[i] = s->has_rates ? chebyshev_at(s, c, counts[4 + i], x) : 0.0;
    c += counts[4 + i];
  }
  if (!q_unit(a->q) || !myb_vec3_finite(a->rate)) {
    return MYB_E_FORMAT;
  }
  a->ticks = at;
  a->found = 1;
  return MYB_OK;
}

/* Inside an interval at the instant; in a gap or beyond the ends the nearer end of an
 * interval within the tolerance (as near: the earlier one). */
static myb_status answer_type4(const ck_segment* s, double ticks, double tolerance, answer* a) {
  const uint64_t n = s->count;
  const uint64_t le = count_le(s, s->reference_base + 1, n, ticks);
  if (le == 0) {
    const uint64_t p = packet4(s, 0);
    const double left = word(s, p) - word(s, p + 1);
    return ticks >= left - tolerance ? evaluate4(s, 0, left, a) : MYB_OK;
  }
  const uint64_t p = packet4(s, le - 1);
  const double right = word(s, p) + word(s, p + 1);
  if (le == n) {
    if (ticks > right + tolerance) {
      return MYB_OK;
    }
    return evaluate4(s, le - 1, ticks >= right ? right : ticks, a);
  }
  if (ticks <= right) {
    return evaluate4(s, le - 1, ticks, a);
  }
  const uint64_t q = packet4(s, le);
  const double left = word(s, q) - word(s, q + 1);
  if (ticks - right <= left - ticks) {
    return ticks <= right + tolerance ? evaluate4(s, le - 1, right, a) : MYB_OK;
  }
  return ticks >= left - tolerance ? evaluate4(s, le, left, a) : MYB_OK;
}

/* ---- Types 5 and 6: windows of packets ----------------------------------------------------- */

/* The packets an instant is interpolated from. */
typedef struct ck_window {
  double t;
  int subtype;
  int count;
  double rate; /* seconds per tick */
  double ticks[MYB_CK_MAX_WINDOW];
  double packets[MYB_CK_MAX_WINDOW * 7]; /* at most 24 of 7 or 12 of 14 */
} ck_window;

/* Copies `count` packets from index `first` of the packets at `packets` (ticks at `ticks`). */
static void fill_window(const ck_segment* s, uint64_t packets, uint64_t ticks, uint64_t first,
                        int count, ck_window* w) {
  const uint64_t size = packet_size[w->subtype];
  w->count = count;
  for (int i = 0; i < count; ++i) {
    w->ticks[i] = word(s, ticks + first + (uint64_t)i);
    for (uint64_t j = 0; j < size; ++j) {
      w->packets[(uint64_t)i * size + j] = word(s, packets + (first + (uint64_t)i) * size + j);
    }
  }
}

/* Lagrange interpolation at 0 of the m points (d_i, y_i), Neville's algorithm carried with
 * the derivative. */
static void lagrange(int m, const double* d, const double* y, double* value,
                     double* derivative) {
  double p[MYB_CK_MAX_WINDOW];
  double dp[MYB_CK_MAX_WINDOW];
  for (int i = 0; i < m; ++i) {
    p[i] = y[i];
    dp[i] = 0.0;
  }
  for (int k = 1; k < m; ++k) {
    for (int i = 0; i + k < m; ++i) {
      const double span = d[i + k] - d[i];
      dp[i] = (p[i + 1] - d[i] * dp[i + 1] - p[i] + d[i + k] * dp[i]) / span;
      p[i] = (d[i + k] * p[i] - d[i] * p[i + 1]) / span;
    }
  }
  *value = p[0];
  *derivative = dp[0];
}

/* The angular velocity from a unit quaternion and its derivative: -2 Im(q* dq). */
static void quaternion_rate(const double q[4], const double dq[4], double rate[3]) {
  const double conjugate[4] = {q[0], -q[1], -q[2], -q[3]};
  double product[4];
  q_multiply(conjugate, dq, product);
  for (int i = 0; i < 3; ++i) {
    rate[i] = -2.0 * product[1 + i];
  }
}

/* The quaternion and the angular velocity of a window at its instant. */
static myb_status evaluate_window(ck_window* w, answer* a) {
  const int m = w->count;
  const uint64_t size = packet_size[w->subtype];
  const int lagrangian = w->subtype == 1 || w->subtype == 3;
  /* each quaternion the nearer of q and -q to the one before it */
  for (int i = 1; i < m; ++i) {
    const double* before = w->packets + (uint64_t)(i - 1) * size;
    double* q = w->packets + (uint64_t)i * size;
    double same = 0.0;
    double opposite = 0.0;
    for (int j = 0; j < 4; ++j) {
      same += (before[j] - q[j]) * (before[j] - q[j]);
      opposite += (before[j] + q[j]) * (before[j] + q[j]);
    }
    if (sqrt(opposite) < sqrt(same)) {
      if (!lagrangian) {
        return MYB_E_FORMAT; /* Hermite data must come with consistent signs */
      }
      for (int j = 0; j < 4; ++j) {
        q[j] = -q[j];
      }
    }
  }
  double d[MYB_CK_MAX_WINDOW];
  for (int i = 0; i < m; ++i) {
    d[i] = w->ticks[i] - w->t;
  }
  double s[4];
  double ds[4];
  double y[MYB_CK_MAX_WINDOW];
  double yd[MYB_CK_MAX_WINDOW];
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < m; ++i) {
      y[i] = w->packets[(uint64_t)i * size + (uint64_t)j];
      yd[i] = lagrangian ? 0.0 : w->packets[(uint64_t)i * size + 4 + (uint64_t)j] * w->rate;
    }
    if (lagrangian) {
      lagrange(m, d, y, &s[j], &ds[j]);
    } else {
      myb_internal_hermite(m, d, y, yd, &s[j], &ds[j]);
    }
  }
  const double norm = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + s[3] * s[3]);
  if (!(norm > 0.0) || !isfinite(norm)) {
    return MYB_E_FORMAT;
  }
  for (int j = 0; j < 4; ++j) {
    a->q[j] = s[j] / norm;
  }
  if (w->subtype == 0 || w->subtype == 1) {
    /* d(s / |s|) = ds / |s| - (q . ds) / |s| q, per tick; then per second */
    const double along = (a->q[0] * ds[0] + a->q[1] * ds[1] + a->q[2] * ds[2] +
                          a->q[3] * ds[3]) / norm;
    double dq[4];
    for (int j = 0; j < 4; ++j) {
      dq[j] = ds[j] / norm - along * a->q[j];
    }
    quaternion_rate(a->q, dq, a->rate);
    for (int j = 0; j < 3; ++j) {
      a->rate[j] /= w->rate;
    }
  } else {
    const uint64_t from = w->subtype == 2 ? 8 : 4;
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < m; ++i) {
        y[i] = w->packets[(uint64_t)i * size + from + (uint64_t)j];
        yd[i] = w->subtype == 2 ? w->packets[(uint64_t)i * size + 11 + (uint64_t)j] * w->rate
                                : 0.0;
      }
      double unused = 0.0;
      if (w->subtype == 2) {
        myb_internal_hermite(m, d, y, yd, &a->rate[j], &unused);
      } else {
        lagrange(m, d, y, &a->rate[j], &unused);
      }
    }
  }
  if (!myb_vec3_finite(a->rate)) {
    return MYB_E_FORMAT;
  }
  a->ticks = w->t;
  a->found = 1;
  return MYB_OK;
}

/* Type 5: the window centred on the instant inside its interval, made smaller near the ends
 * of the interval; in a gap the nearer end within the tolerance (as near: the earlier). */
static myb_status answer_type5(const ck_segment* s, double ticks, double tolerance, answer* a) {
  const uint64_t n = s->count;
  double t = fmin(fmax(ticks, s->start), s->stop);
  const double first_tick = word(s, s->ticks);
  const double last_tick = word(s, s->ticks + n - 1);
  if (t < first_tick) {
    if (ticks + tolerance < first_tick) {
      return MYB_OK;
    }
    t = first_tick;
  } else if (t > last_tick) {
    if (ticks - tolerance > last_tick) {
      return MYB_OK;
    }
    t = last_tick;
  }
  /* 1-based indices of the two ticks around t */
  const uint64_t below = count_lt(s, s->ticks, n, t);
  uint64_t low = below == 0 ? 1 : below;
  uint64_t high = below == 0 ? (n == 1 ? 1 : 2) : below + 1;
  double low_tick = word(s, s->ticks + low - 1);
  double high_tick = word(s, s->ticks + high - 1);
  /* the interval: its start at or before t, the starts of the next two */
  const uint64_t j = count_le(s, s->starts, s->intervals, t);
  double start = word(s, s->starts + j - 1);
  double next = j < s->intervals ? word(s, s->starts + j) : HUGE_VAL;
  const double after_next = j + 1 < s->intervals ? word(s, s->starts + j + 1) : HUGE_VAL;
  if (high_tick == next) {
    /* t lies between two intervals */
    if (fabs(t - low_tick) <= fabs(high_tick - t)) {
      if (fabs(t - low_tick) > tolerance) {
        return MYB_OK;
      }
      t = low_tick;
      high = low;
      high_tick = low_tick;
    } else {
      if (fabs(high_tick - t) > tolerance) {
        return MYB_OK;
      }
      t = high_tick;
      low = high;
      low_tick = high_tick;
      start = next;
      next = after_next;
    }
  }
  const uint64_t half = (uint64_t)s->window / 2;
  const uint64_t left = half < low ? half : low;
  const uint64_t right = half < n - high + 1 ? half : n - high + 1;
  uint64_t first = low - left + 1;
  const uint64_t size = left + right;
  /* only the ticks of the interval */
  uint64_t i = 0;
  uint64_t k = 0;
  for (uint64_t m = 0; m < size; ++m) {
    const double tick = word(s, s->ticks + first - 1 + m);
    i += tick < start ? 1 : 0;
    k += tick < next ? 1 : 0;
  }
  if (i == k) {
    return MYB_OK;
  }
  first += i;
  ck_window w = {.t = t, .subtype = s->subtype, .rate = s->clock_rate};
  fill_window(s, s->begin, s->ticks, first - 1, (int)(k - i), &w);
  return evaluate_window(&w, a);
}

/* Type 6: the mini-segment whose interval holds the instant (on a common boundary the
 * earlier or the later by the flag of the segment); in a gap at the end of a mini-segment
 * its last tick or the start of the next one, whichever is nearer, within the tolerance.
 * The window is centred on the instant and made smaller only at the ends of the ticks. */
static myb_status answer_type6(const ck_segment* s, double ticks, double tolerance, answer* a) {
  const uint64_t n = s->count;
  double t = fmin(fmax(ticks, s->start), s->stop);
  uint64_t k = 0; /* 0-based mini-segment */
  if (s->select_last) {
    const uint64_t end = count_le(s, s->bounds, n + 1, t) + 1;
    const uint64_t clipped = end < n + 1 ? end : n + 1;
    if (clipped < 2) {
      return MYB_E_FORMAT;
    }
    k = clipped - 2;
  } else {
    uint64_t begin = count_lt(s, s->bounds, n + 1, t);
    begin = begin > 1 ? begin : 1;
    if (begin > n) {
      return MYB_E_FORMAT;
    }
    k = begin - 1;
  }
  uint64_t first = 0;
  uint64_t last = 0;
  mini_bounds(s, k, &first, &last);
  uint64_t count = (uint64_t)word(s, last);
  double last_tick = word(s, last - 4 - (count - 1) / 100);
  if (t > last_tick) {
    /* in the gap at the end of the mini-segment */
    const double stop = word(s, s->bounds + k + 1);
    if (tolerance == 0.0 || (t - last_tick > tolerance && stop - t > tolerance)) {
      return MYB_OK;
    }
    if (k + 1 == n) {
      if (t - last_tick > tolerance) {
        return MYB_OK;
      }
      t = last_tick;
    } else if (t - last_tick <= (stop - last_tick) / 2) {
      t = last_tick;
    } else {
      ++k;
      mini_bounds(s, k, &first, &last);
      count = (uint64_t)word(s, last);
      t = stop;
    }
  }
  ck_window w = {.t = t};
  if (!window_parameters(s, last - 3, &w.rate, &w.subtype, &w.count)) {
    return MYB_E_FORMAT;
  }
  const uint64_t window = (uint64_t)w.count;
  const uint64_t size = packet_size[w.subtype];
  const uint64_t tick0 = first + count * size;
  const uint64_t below = count_lt(s, tick0, count, t);
  const uint64_t low = below == 0 ? 1 : below;
  const uint64_t high = low + 1;
  const uint64_t half = window / 2;
  const uint64_t left = half < low ? half : low;
  const uint64_t right = high <= count ? (half < count - high + 1 ? half : count - high + 1) : 0;
  fill_window(s, first, tick0, low - left, (int)(left + right), &w);
  return evaluate_window(&w, a);
}

myb_status myb_ck_pointing_at(const myb_ck* ck, int32_t instrument, double ticks,
                              double tolerance, int need_rate, myb_ck_pointing* pointing) {
  MYB_CHECK(ck != NULL && pointing != NULL);
  if (!(tolerance >= 0.0) || !isfinite(tolerance) || !isfinite(ticks)) {
    return MYB_E_DOMAIN;
  }
  /* bound: the segments, the last loaded first */
  for (size_t k = ck->count; k-- > 0;) {
    const ck_segment* s = &ck->segments[k];
    if (s->instrument != instrument || !s->supported || (need_rate && !s->has_rates) ||
        ticks + tolerance < s->start || ticks - tolerance > s->stop) {
      continue;
    }
    answer a = {0};
    myb_status status = MYB_E_FORMAT;
    if (s->type == 1) {
      status = answer_type1(s, ticks, tolerance, &a);
    } else if (s->type == 2) {
      status = answer_type2(s, ticks, tolerance, &a);
    } else if (s->type == 3) {
      status = answer_type3(s, ticks, tolerance, &a);
    } else if (s->type == 4) {
      status = answer_type4(s, ticks, tolerance, &a);
    } else if (s->type == 5) {
      status = answer_type5(s, ticks, tolerance, &a);
    } else {
      status = answer_type6(s, ticks, tolerance, &a);
    }
    if (status != MYB_OK) {
      return status;
    }
    if (!a.found) {
      continue;
    }
    memcpy(pointing->q, a.q, sizeof a.q);
    memcpy(pointing->rate, a.rate, sizeof a.rate);
    pointing->has_rate = s->has_rates;
    pointing->ticks = a.ticks;
    pointing->frame = s->frame;
    const double c = a.q[0];
    const double x = a.q[1];
    const double y = a.q[2];
    const double z = a.q[3];
    const double m[9] = {1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - c * z), 2.0 * (x * z + c * y),
                         2.0 * (x * y + c * z), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - c * x),
                         2.0 * (x * z - c * y), 2.0 * (y * z + c * x), 1.0 - 2.0 * (x * x + y * y)};
    memcpy(pointing->rotation.m, m, sizeof m);
    /* the rows of the matrix are the axes of the instrument, which turn with the angular
     * velocity: the derivative of a row is rate x row */
    for (int i = 0; i < 3; ++i) {
      myb_vec3_cross(a.rate, m + 3 * i, pointing->rotation.dm + 3 * i);
    }
    return MYB_OK;
  }
  return MYB_E_RANGE;
}

myb_status myb_orientation_ck(void* context, myb_time t, myb_orientation* orientation) {
  const myb_ck_orientation* c = context;
  MYB_CHECK(c != NULL && c->ck != NULL && c->sclk != NULL && orientation != NULL);
  double ticks = 0.0;
  myb_status status = myb_sclk_from_time(c->sclk, t, &ticks);
  if (status != MYB_OK) {
    return status;
  }
  myb_ck_pointing pointing;
  status = myb_ck_pointing_at(c->ck, c->instrument, ticks, c->tolerance, 0, &pointing);
  if (status != MYB_OK) {
    return status;
  }
  if (c->base != 0 && pointing.frame != c->base) {
    return MYB_E_FORMAT;
  }
  orientation->rotation = pointing.rotation;
  orientation->ut = t;
  orientation->xp = 0.0;
  orientation->yp = 0.0;
  return MYB_OK;
}
