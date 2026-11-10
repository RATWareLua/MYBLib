/* SPDX-License-Identifier: BSD-3-Clause */
/* Writing SPK files: segments of discrete states (types 9 and 13) and of Chebyshev
 * polynomials (types 2 and 3) in a NAIF DAF, little-endian, as "SPK Required Reading" and "DAF
 * Required Reading" (NAIF, JPL) lay them out.
 *
 * File: record 1 is the file record; then come the summary records, each followed by the
 * record of the names of its (at most 25) segments, chained in a list; then the data. A
 * summary is 2 doubles (the coverage) and 6 integers (target, centre, frame, type, first and
 * last word of the data) packed into 5 double words; words are addressed from 1. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "myb/frames.h"
#include "myb/kin.h"
#include "platform/platform.h"

enum {
  RECORD_BYTES = 1024,
  RECORD_WORDS = 128,
  SUMMARY_WORDS = 5,
  SUMMARIES_PER_RECORD = 25,
  NAME_CHARS = 40,
  INTERNAL_NAME_CHARS = 60,
  MAX_DEGREE = 27,       /* the limit of the SPICE writers of these types */
  DIRECTORY_SPACING = 100 /* every 100th epoch of types 9 and 13 */
};

static const double metres_per_km = 1000.0;
static const double pi = 3.141592653589793;

typedef struct segment {
  int32_t target;
  int32_t center;
  int32_t frame;
  int32_t type;
  double start;
  double stop;
  char name[NAME_CHARS + 1];
  size_t offset; /* first word of its data in writer->data */
  size_t length;
} segment;

enum { KIND_SPK = 0, KIND_PCK = 1 };

struct myb_spk_writer {
  myb_allocator allocator;
  int kind;
  char internal_name[INTERNAL_NAME_CHARS + 1];
  segment* segments;
  size_t count;
  size_t capacity;
  double* data;
  size_t words;
  size_t word_capacity;
};

static int printable(const char* s, size_t max) {
  const size_t n = strlen(s);
  if (n == 0 || n > max) {
    return 0;
  }
  for (size_t i = 0; i < n; ++i) {
    if (s[i] < 32 || s[i] > 126) {
      return 0;
    }
  }
  return 1;
}

static myb_status create(const char* internal_name, const myb_allocator* allocator, int kind,
                         myb_spk_writer** writer) {
  const char* name = internal_name != NULL ? internal_name : "myb";
  if (!printable(name, INTERNAL_NAME_CHARS)) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_spk_writer), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_spk_writer* w = block;
  memset(w, 0, sizeof *w);
  w->allocator = a;
  w->kind = kind;
  memcpy(w->internal_name, name, strlen(name) + 1);
  *writer = w;
  return MYB_OK;
}

myb_status myb_spk_writer_create(const char* internal_name, const myb_allocator* allocator,
                                 myb_spk_writer** writer) {
  MYB_CHECK(writer != NULL);
  return create(internal_name, allocator, KIND_SPK, writer);
}

void myb_spk_writer_destroy(myb_spk_writer* writer) {
  if (writer == NULL) {
    return;
  }
  const myb_allocator a = writer->allocator;
  myb_internal_release(&a, writer->segments);
  myb_internal_release(&a, writer->data);
  myb_internal_release(&a, writer);
}

/* Room for `extra` more words of data and one more segment. */
static myb_status reserve(myb_spk_writer* w, size_t extra) {
  if (extra > SIZE_MAX / 16 - w->words) {
    return MYB_E_OVERFLOW;
  }
  if (w->words + extra > w->word_capacity) {
    size_t capacity = w->word_capacity < 1024 ? 1024 : w->word_capacity;
    /* bound: doubling reaches any size below SIZE_MAX / 16 in 64 steps */
    for (int step = 0; step < 64 && capacity < w->words + extra; ++step) {
      capacity *= 2;
    }
    void* block = NULL;
    const myb_status status =
        myb_internal_alloc_array(&w->allocator, capacity, sizeof(double), &block);
    if (status != MYB_OK) {
      return status;
    }
    if (w->words > 0) {
      memcpy(block, w->data, w->words * sizeof(double));
    }
    myb_internal_release(&w->allocator, w->data);
    w->data = block;
    w->word_capacity = capacity;
  }
  if (w->count == w->capacity) {
    const size_t capacity = w->capacity == 0 ? 8 : 2 * w->capacity;
    void* block = NULL;
    const myb_status status =
        myb_internal_alloc_array(&w->allocator, capacity, sizeof(segment), &block);
    if (status != MYB_OK) {
      return status;
    }
    if (w->count > 0) {
      memcpy(block, w->segments, w->count * sizeof(segment));
    }
    myb_internal_release(&w->allocator, w->segments);
    w->segments = block;
    w->capacity = capacity;
  }
  return MYB_OK;
}

static int valid_header(const myb_spk_segment_header* h) {
  return h->name != NULL && printable(h->name, NAME_CHARS);
}

/* The segment whose `length` words were just appended to the data. */
static void commit(myb_spk_writer* w, const myb_spk_segment_header* h, int type, double start,
                   double stop, size_t length) {
  segment* s = &w->segments[w->count++];
  memset(s, 0, sizeof *s);
  s->target = h->target;
  s->center = h->center;
  s->frame = h->frame;
  s->type = type;
  s->start = start;
  s->stop = stop;
  memcpy(s->name, h->name, strlen(h->name) + 1);
  s->offset = w->words;
  s->length = length;
  w->words += length;
}

myb_status myb_spk_writer_add_states(myb_spk_writer* writer,
                                     const myb_spk_segment_header* header, int type, int degree,
                                     const double* epochs, const double* states, size_t count) {
  MYB_CHECK(writer != NULL && header != NULL && epochs != NULL && states != NULL);
  if (!valid_header(header) || (type != 9 && type != 13) || degree < 1 || degree > MAX_DEGREE ||
      (type == 13 && degree % 2 == 0) || count < 2 || count > SIZE_MAX / 64) {
    return MYB_E_DOMAIN;
  }
  /* type 9: a polynomial through degree + 1 states; type 13 (Hermite): through (degree + 1) / 2
   * states with their velocities */
  const size_t window = type == 9 ? (size_t)degree + 1 : ((size_t)degree + 1) / 2;
  if (window > count) {
    return MYB_E_DOMAIN;
  }
  for (size_t k = 0; k < count; ++k) {
    if (!isfinite(epochs[k]) || (k > 0 && !(epochs[k] > epochs[k - 1]))) {
      return MYB_E_DOMAIN;
    }
    for (int i = 0; i < 6; ++i) {
      if (!isfinite(states[6 * k + (size_t)i])) {
        return MYB_E_DOMAIN;
      }
    }
  }
  const size_t directory = (count - 1) / DIRECTORY_SPACING;
  const size_t length = 7 * count + directory + 2;
  const myb_status status = reserve(writer, length);
  if (status != MYB_OK) {
    return status;
  }
  double* d = writer->data + writer->words;
  for (size_t k = 0; k < 6 * count; ++k) {
    *d++ = states[k] / metres_per_km;
  }
  for (size_t k = 0; k < count; ++k) {
    *d++ = epochs[k];
  }
  for (size_t k = 1; k <= directory; ++k) {
    *d++ = epochs[k * DIRECTORY_SPACING - 1];
  }
  *d++ = (double)(window - 1);
  *d++ = (double)count;
  commit(writer, header, type, epochs[0], epochs[count - 1], length);
  return MYB_OK;
}

/* sum c_k T_k(x) by Clenshaw's recurrence. */
myb_status myb_spk_writer_add_raw(myb_spk_writer* writer, const myb_spk_segment_header* header,
                                  int type, double start, double stop, const double* words,
                                  size_t count) {
  MYB_CHECK(writer != NULL && header != NULL && words != NULL);
  if (!valid_header(header) || type < 1 || !isfinite(start) || !isfinite(stop) ||
      !(stop >= start) || count < 1) {
    return MYB_E_DOMAIN;
  }
  const myb_status status = reserve(writer, count);
  if (status != MYB_OK) {
    return status;
  }
  memcpy(writer->data + writer->words, words, count * sizeof(double));
  commit(writer, header, type, start, stop, count);
  return MYB_OK;
}

static double chebyshev_value(const double* c, int n, double x) {
  double b1 = 0.0;
  double b2 = 0.0;
  for (int k = n - 1; k >= 1; --k) {
    const double b0 = 2.0 * x * b1 - b2 + c[k];
    b2 = b1;
    b1 = b0;
  }
  return x * b1 - b2 + c[0];
}

myb_status myb_spk_writer_add_chebyshev(myb_spk_writer* writer,
                                        const myb_spk_segment_header* header, int type,
                                        myb_state_function state, void* context, double start,
                                        double stop, size_t intervals, int degree,
                                        double* max_error) {
  MYB_CHECK(writer != NULL && header != NULL && state != NULL);
  if (!valid_header(header) || (type != 2 && type != 3) || degree < 1 || degree > MAX_DEGREE ||
      !isfinite(start) || !isfinite(stop) || !(stop > start) || intervals < 1 ||
      intervals > SIZE_MAX / 4096) {
    return MYB_E_DOMAIN;
  }
  const int n = degree + 1;             /* coefficients per component, and nodes */
  const int components = type == 2 ? 3 : 6;
  const size_t rsize = 2 + (size_t)components * (size_t)n;
  const size_t length = rsize * intervals + 4;
  myb_status status = reserve(writer, length);
  if (status != MYB_OK) {
    return status;
  }
  /* cos(pi j / (2 n)), j = 0 .. 4 n - 1: every cosine the transform needs */
  double cosine[4 * (MAX_DEGREE + 1)];
  for (int j = 0; j < 4 * n; ++j) {
    cosine[j] = myb_cos(pi * (double)j / (2.0 * (double)n));
  }
  const double step = (stop - start) / (double)intervals;
  double* d = writer->data + writer->words;
  double worst = 0.0;
  for (size_t interval = 0; interval < intervals; ++interval) {
    const double a = start + (stop - start) * ((double)interval / (double)intervals);
    const double radius = 0.5 * step;
    const double mid = a + radius;
    double f[MAX_DEGREE + 1][6];
    /* the states at the nodes x_i = cos(pi (2 i + 1) / (2 n)), in km and km/s */
    for (int i = 0; i < n; ++i) {
      double y[6];
      status = state(context, myb_time_from_seconds(mid + radius * cosine[2 * i + 1]), y);
      if (status != MYB_OK) {
        return status;
      }
      for (int c = 0; c < 6; ++c) {
        if (!isfinite(y[c])) {
          return MYB_E_DOMAIN;
        }
        f[i][c] = y[c] / metres_per_km;
      }
    }
    double* record = d;
    *d++ = mid;
    *d++ = radius;
    /* c_k = (2 / n) sum_i f_i cos(k pi (2 i + 1) / (2 n)), half of it for k = 0: the
     * polynomial of degree n - 1 through the nodes */
    for (int c = 0; c < components; ++c) {
      for (int k = 0; k < n; ++k) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
          sum += f[i][c] * cosine[(k * (2 * i + 1)) % (4 * n)];
        }
        *d++ = (k == 0 ? 1.0 : 2.0) * sum / (double)n;
      }
    }
    if (max_error != NULL) {
      /* between the nodes: at the extrema of T_n, x = cos(pi j / n) */
      for (int j = 0; j <= n; ++j) {
        const double x = cosine[2 * j];
        double y[6];
        status = state(context, myb_time_from_seconds(mid + radius * x), y);
        if (status != MYB_OK) {
          return status;
        }
        double square = 0.0;
        for (int c = 0; c < 3; ++c) {
          const double difference =
              metres_per_km * chebyshev_value(record + 2 + c * n, n, x) - y[c];
          square += difference * difference;
        }
        const double error = sqrt(square);
        worst = error > worst ? error : worst;
      }
    }
  }
  *d++ = start;
  *d++ = step;
  *d++ = (double)rsize;
  *d++ = (double)intervals;
  commit(writer, header, type, start, stop, length);
  if (max_error != NULL) {
    *max_error = worst;
  }
  return MYB_OK;
}

/* ---- The file ------------------------------------------------------------------------------- */

static size_t summary_records(const myb_spk_writer* w) {
  return (w->count + SUMMARIES_PER_RECORD - 1) / SUMMARIES_PER_RECORD;
}

/* First word (1-based) of the data area: after the file record and the summary and name
 * records. */
static size_t data_start(const myb_spk_writer* w) {
  return (1 + 2 * summary_records(w)) * RECORD_WORDS + 1;
}

size_t myb_spk_writer_size(const myb_spk_writer* writer) {
  MYB_CHECK(writer != NULL);
  if (writer->count == 0) {
    return 0;
  }
  const size_t words = data_start(writer) - 1 + writer->words;
  const size_t records = (words + RECORD_WORDS - 1) / RECORD_WORDS;
  return records * RECORD_BYTES;
}

static void store32(unsigned char* p, int32_t value) {
  const uint32_t u = (uint32_t)value;
  for (int i = 0; i < 4; ++i) {
    p[i] = (unsigned char)(u >> (8 * i));
  }
}

static void store_double(unsigned char* p, double value) {
  uint64_t u = 0;
  memcpy(&u, &value, sizeof u);
  for (int i = 0; i < 8; ++i) {
    p[i] = (unsigned char)(u >> (8 * i));
  }
}

myb_status myb_spk_writer_write(const myb_spk_writer* writer, void* buffer, size_t size) {
  MYB_CHECK(writer != NULL && buffer != NULL);
  const myb_spk_writer* w = writer;
  const size_t total = myb_spk_writer_size(w);
  if (w->count == 0) {
    return MYB_E_DOMAIN;
  }
  if (size < total) {
    return MYB_E_OVERFLOW;
  }
  const size_t first_data = data_start(w);
  if (first_data + w->words > (size_t)INT32_MAX) {
    return MYB_E_OVERFLOW; /* addresses are 32-bit integers */
  }
  unsigned char* out = buffer;
  memset(out, 0, total);
  /* the file record */
  static const unsigned char ftp[28] = {'F', 'T', 'P', 'S', 'T', 'R', ':', '\r', ':', '\n',
                                        ':', '\r', '\n', ':', '\r', 0x00, ':', 0x81, ':', 0x10,
                                        0xCE, ':', 'E', 'N', 'D', 'F', 'T', 'P'};
  const size_t groups = summary_records(w);
  memcpy(out, w->kind == KIND_SPK ? "DAF/SPK " : "DAF/PCK ", 8);
  store32(out + 8, 2);
  store32(out + 12, w->kind == KIND_SPK ? 6 : 5);
  memset(out + 16, ' ', INTERNAL_NAME_CHARS);
  memcpy(out + 16, w->internal_name, strlen(w->internal_name));
  store32(out + 76, 2);                                  /* first summary record */
  store32(out + 80, (int32_t)(2 * groups));              /* last summary record */
  store32(out + 84, (int32_t)(first_data + w->words));   /* first free word */
  memcpy(out + 88, "LTL-IEEE", 8);
  memcpy(out + 699, ftp, sizeof ftp);
  /* summary records with their name records */
  for (size_t g = 0; g < groups; ++g) {
    unsigned char* summary = out + (1 + 2 * g) * RECORD_BYTES;
    unsigned char* names = summary + RECORD_BYTES;
    const size_t first = g * SUMMARIES_PER_RECORD;
    const size_t in_group =
        w->count - first < SUMMARIES_PER_RECORD ? w->count - first : SUMMARIES_PER_RECORD;
    store_double(summary, g + 1 < groups ? (double)(2 * (g + 1) + 2) : 0.0); /* next */
    store_double(summary + 8, g > 0 ? (double)(2 * (g - 1) + 2) : 0.0);     /* previous */
    store_double(summary + 16, (double)in_group);
    memset(names, ' ', RECORD_BYTES);
    for (size_t k = 0; k < in_group; ++k) {
      const segment* s = &w->segments[first + k];
      unsigned char* p = summary + 24 + k * SUMMARY_WORDS * 8;
      const size_t begin = first_data + s->offset;
      store_double(p, s->start);
      store_double(p + 8, s->stop);
      if (w->kind == KIND_SPK) {
        store32(p + 16, s->target);
        store32(p + 20, s->center);
        store32(p + 24, s->frame);
        store32(p + 28, s->type);
        store32(p + 32, (int32_t)begin);
        store32(p + 36, (int32_t)(begin + s->length - 1));
      } else { /* the frame, its base, the type and the addresses: five integers */
        store32(p + 16, s->target);
        store32(p + 20, s->frame);
        store32(p + 24, s->type);
        store32(p + 28, (int32_t)begin);
        store32(p + 32, (int32_t)(begin + s->length - 1));
      }
      memcpy(names + k * NAME_CHARS, s->name, strlen(s->name));
    }
  }
  unsigned char* data = out + (first_data - 1) * 8;
  for (size_t k = 0; k < w->words; ++k) {
    store_double(data + 8 * k, w->data[k]);
  }
  return MYB_OK;
}

myb_status myb_spk_writer_write_file(const myb_spk_writer* writer, const char* path) {
  MYB_CHECK(writer != NULL && path != NULL);
  const size_t size = myb_spk_writer_size(writer);
  if (size == 0) {
    return MYB_E_DOMAIN;
  }
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&writer->allocator, size, 1, &block);
  if (status != MYB_OK) {
    return status;
  }
  status = myb_spk_writer_write(writer, block, size);
  if (status == MYB_OK) {
    status = myb_platform_write_file(path, &writer->allocator, block, size);
  }
  myb_internal_release(&writer->allocator, block);
  return status;
}

/* ---- Binary PCK (myb/frames.h) ------------------------------------------------------------------
 * The same file with another label and summaries of five integers: the writer is the same
 * object. */

struct myb_pck_writer {
  myb_spk_writer w;
};

myb_status myb_pck_writer_create(const char* internal_name, const myb_allocator* allocator,
                                 myb_pck_writer** writer) {
  MYB_CHECK(writer != NULL);
  myb_spk_writer* w = NULL;
  const myb_status status = create(internal_name, allocator, KIND_PCK, &w);
  if (status == MYB_OK) {
    *writer = (myb_pck_writer*)w; /* its only member */
  }
  return status;
}

void myb_pck_writer_destroy(myb_pck_writer* writer) {
  if (writer != NULL) {
    myb_spk_writer_destroy(&writer->w);
  }
}

/* The angles of M = R3(a[2]) R1(a[1]) R3(a[0]): the third row of M is
 * (sin a1 sin a0, -sin a1 cos a0, cos a1), its third column (sin a2 sin a1, cos a2 sin a1,
 * cos a1). 0 if they are not defined. */
static int euler313(const double m[9], double a[3]) {
  for (int i = 0; i < 9; ++i) {
    if (!isfinite(m[i])) {
      return 0;
    }
  }
  const double sine = myb_hypot(m[6], m[7]);
  if (!(sine > 1e-9)) {
    return 0;
  }
  a[0] = myb_atan2(m[6], -m[7]);
  a[1] = myb_atan2(sine, m[8]);
  a[2] = myb_atan2(m[2], m[5]);
  return 1;
}

myb_status myb_pck_writer_add_chebyshev(myb_pck_writer* writer, int32_t frame, int32_t base,
                                        const char* name, myb_orientation_function orientation,
                                        void* context, double start, double stop,
                                        size_t intervals, int degree, double* max_error) {
  MYB_CHECK(writer != NULL && orientation != NULL);
  myb_spk_writer* w = &writer->w;
  const myb_spk_segment_header header = {.target = frame, .center = 0, .frame = base, .name = name};
  if (!valid_header(&header) || degree < 1 || degree > MAX_DEGREE || !isfinite(start) ||
      !isfinite(stop) || !(stop > start) || intervals < 1 || intervals > SIZE_MAX / 4096) {
    return MYB_E_DOMAIN;
  }
  const int n = degree + 1;
  const size_t rsize = 2 + 3 * (size_t)n;
  const size_t length = rsize * intervals + 4;
  myb_status status = reserve(w, length);
  if (status != MYB_OK) {
    return status;
  }
  double cosine[4 * (MAX_DEGREE + 1)];
  for (int j = 0; j < 4 * n; ++j) {
    cosine[j] = myb_cos(pi * (double)j / (2.0 * (double)n));
  }
  const double step = (stop - start) / (double)intervals;
  double* d = w->data + w->words;
  double worst = 0.0;
  double last[3] = {0.0, 0.0, 0.0};
  int following = 0;
  for (size_t interval = 0; interval < intervals; ++interval) {
    const double a = start + (stop - start) * ((double)interval / (double)intervals);
    const double radius = 0.5 * step;
    const double mid = a + radius;
    double f[MAX_DEGREE + 1][3];
    /* the angles at the nodes x_i = cos(pi (2 i + 1) / (2 n)), in the order of time, each
     * the turn nearest to the one before */
    for (int i = n - 1; i >= 0; --i) {
      myb_orientation o;
      status = orientation(context, myb_time_from_seconds(mid + radius * cosine[2 * i + 1]), &o);
      if (status != MYB_OK) {
        return status;
      }
      if (!euler313(o.rotation.m, f[i])) {
        return MYB_E_DOMAIN;
      }
      for (int c = 0; c < 3 && following; ++c) {
        f[i][c] += 2.0 * pi * floor((last[c] - f[i][c]) / (2.0 * pi) + 0.5);
      }
      memcpy(last, f[i], sizeof last);
      following = 1;
    }
    double* record = d;
    *d++ = mid;
    *d++ = radius;
    for (int c = 0; c < 3; ++c) {
      for (int k = 0; k < n; ++k) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
          sum += f[i][c] * cosine[(k * (2 * i + 1)) % (4 * n)];
        }
        *d++ = (k == 0 ? 1.0 : 2.0) * sum / (double)n;
      }
    }
    if (max_error != NULL) {
      /* between the nodes: at the extrema of T_n */
      for (int j = 0; j <= n; ++j) {
        const double x = cosine[2 * j];
        myb_orientation o;
        status = orientation(context, myb_time_from_seconds(mid + radius * x), &o);
        if (status != MYB_OK) {
          return status;
        }
        const double angles[3] = {chebyshev_value(record + 2, n, x),
                                  chebyshev_value(record + 2 + n, n, x),
                                  chebyshev_value(record + 2 + 2 * n, n, x)};
        const double still[3] = {0.0, 0.0, 0.0};
        myb_rotation fitted;
        myb_rotation_euler313(angles, still, &fitted);
        /* a small turn by the angle e changes the matrix by sqrt(2) e (Frobenius) */
        double square = 0.0;
        for (int i = 0; i < 9; ++i) {
          const double difference = fitted.m[i] - o.rotation.m[i];
          square += difference * difference;
        }
        const double error = sqrt(0.5 * square);
        worst = error > worst ? error : worst;
      }
    }
  }
  *d++ = start;
  *d++ = step;
  *d++ = (double)rsize;
  *d++ = (double)intervals;
  commit(w, &header, 2, start, stop, length);
  if (max_error != NULL) {
    *max_error = worst;
  }
  return MYB_OK;
}

myb_status myb_pck_writer_add_raw(myb_pck_writer* writer, int32_t frame, int32_t base,
                                  const char* name, int type, double start, double stop,
                                  const double* words, size_t count) {
  MYB_CHECK(writer != NULL);
  const myb_spk_segment_header header = {.target = frame, .center = 0, .frame = base, .name = name};
  return myb_spk_writer_add_raw(&writer->w, &header, type, start, stop, words, count);
}

size_t myb_pck_writer_size(const myb_pck_writer* writer) {
  MYB_CHECK(writer != NULL);
  return myb_spk_writer_size(&writer->w);
}

myb_status myb_pck_writer_write(const myb_pck_writer* writer, void* buffer, size_t size) {
  MYB_CHECK(writer != NULL);
  return myb_spk_writer_write(&writer->w, buffer, size);
}

myb_status myb_pck_writer_write_file(const myb_pck_writer* writer, const char* path) {
  MYB_CHECK(writer != NULL);
  return myb_spk_writer_write_file(&writer->w, path);
}
