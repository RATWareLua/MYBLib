/* SPDX-License-Identifier: BSD-3-Clause */
/* Reading NAIF CK files (myb/ck.h): segments of types 1, 2 and 3 in a DAF, after "CK
 * Required Reading" (NAIF, JPL). Which record answers a request — inside an interval, at
 * its ends, within the tolerance — and how the pointing is carried between records were
 * checked against ckgpav_c of CSPICE N0067 (tests/data/ck).
 *
 * Layouts (words of a segment, N pointings, R = 7 with angular velocities and 4 without):
 *   type 1: N records of R; N ticks; (N - 1) / 100 directory ticks; N
 *   type 2: N records of 8 (quaternion, angular velocity, seconds per tick); N starts;
 *           N stops; (N - 1) / 100 directory ticks
 *   type 3: N records of R; N ticks; (N - 1) / 100 directory ticks; M interval starts;
 *           (M - 1) / 100 directory ticks; M; N */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "io/daf.h"
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
  uint64_t intervals; /* type 3: how many, and the address of the first start */
  uint64_t starts;
} ck_segment;

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
    } else {
      status = answer_type3(s, ticks, tolerance, &a);
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
