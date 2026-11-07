/* SPDX-License-Identifier: BSD-3-Clause */
/* Bit-exact reproducibility: a fixed workload over the whole public API is hashed and the hash
 * must be the same on every platform, compiler version, optimisation level and ISA level.
 * The golden value was produced on Windows 10 / Clang 23 / x86-64 (windows-dev preset); Linux
 * builds must reproduce it unchanged. A different hash means that some result changed in at
 * least one bit: find which one before touching the golden value.
 *
 * Inputs come from SplitMix64 and are built from integer bits, so the workload itself does
 * not depend on any library function. */
#include "core/fpguard.h"
#include "mtest.h"
#include "myb/cmath.h"
#include "myb/io.h"
#include "myb/kepler.h"
#include "myb/lambert.h"
#include "myb/num.h"
#include "myb/time.h"

static const uint64_t golden_hash = UINT64_C(0xd1b4b6c88d9abd1c);

/* FNV-1a, 64 bit, over the IEEE-754 encoding in little-endian byte order. */
typedef struct hasher {
  uint64_t state;
  long values;
} hasher;

static void hash_double(hasher* h, double x) {
  const uint64_t bits = mt_bits(x);
  for (int byte = 0; byte < 8; ++byte) {
    h->state ^= (bits >> (8 * byte)) & 0xFF;
    h->state *= UINT64_C(0x100000001B3);
  }
  h->values++;
}

static void hash_dd(hasher* h, myb_dd x) {
  hash_double(h, x.hi);
  hash_double(h, x.lo);
}

static void workload_cmath(hasher* h, mt_rng* rng) {
  for (int i = 0; i < 20000; ++i) {
    const double moderate = mt_double(rng, -30, 20);
    const double large = mt_double(rng, 20, 1023);
    hash_double(h, myb_sin(moderate));
    hash_double(h, myb_cos(moderate));
    hash_double(h, myb_sin(large));
    hash_double(h, myb_cos(large));
    /* One generator call per statement: the evaluation order of function arguments is
     * unspecified in C, and the input sequence must not depend on the compiler. */
    const double a = mt_double(rng, -40, 40);
    const double b = mt_double(rng, -40, 40);
    const double c = mt_double(rng, -80, 80);
    hash_double(h, myb_fma(a, b, c));
    const double any_a = mt_any_double(rng);
    const double any_b = mt_any_double(rng);
    const double any_c = mt_any_double(rng);
    hash_double(h, myb_fma(any_a, any_b, any_c));
  }
}

static void workload_dd(hasher* h, mt_rng* rng) {
  myb_dd chain = myb_dd_from_double(1.0);
  for (int i = 0; i < 20000; ++i) {
    const double a = mt_double(rng, -10, 10);
    double err;
    const double s = myb_two_sum(a, mt_double(rng, -40, 10), &err);
    hash_double(h, s);
    hash_double(h, err);
    const double p = myb_two_prod(a, mt_double(rng, -40, 40), &err);
    hash_double(h, p);
    hash_double(h, err);
    const myb_dd x = {.hi = s, .lo = 0.0};
    const myb_dd y = myb_dd_add_d(myb_dd_from_double(mt_double(rng, -5, 5)), err);
    hash_dd(h, myb_dd_add(x, y));
    hash_dd(h, myb_dd_mul(x, y));
    hash_dd(h, myb_dd_div(x, y));
    hash_dd(h, myb_dd_div_d(x, a));
    hash_dd(h, myb_dd_mul_d(y, a));
    hash_dd(h, myb_dd_sqrt(myb_dd_mul(y, y)));
    chain = myb_dd_div(myb_dd_add(myb_dd_mul(chain, y), x), myb_dd_add_d(y, 7.0));
    hash_dd(h, chain);
  }
}

static void workload_sums(hasher* h, mt_rng* rng) {
  enum { N = 1000 };

  static double x[N];
  static double y[N];
  for (int trial = 0; trial < 100; ++trial) {
    for (int i = 0; i < N; ++i) {
      x[i] = mt_double(rng, -60, 60);
      y[i] = mt_double(rng, -60, 60);
    }
    hash_double(h, myb_sum2(x, N));
    hash_double(h, myb_dot2(x, y, N));
  }
}

/* Decimal strings built digit by digit (no printf: its rounding and locale are not ours). */
static void workload_parse(hasher* h, mt_rng* rng) {
  char text[64];
  for (int i = 0; i < 20000; ++i) {
    size_t n = 0;
    const int digits = mt_int(rng, 1, 30);
    for (int d = 0; d < digits; ++d) {
      text[n++] = (char)('0' + (int)(mt_next(rng) % 10));
      if (d == 0) {
        text[n++] = '.';
      }
    }
    text[n++] = 'e';
    int exponent = mt_int(rng, -330, 310);
    if (exponent < 0) {
      text[n++] = '-';
      exponent = -exponent;
    }
    text[n++] = (char)('0' + exponent / 100);
    text[n++] = (char)('0' + exponent / 10 % 10);
    text[n++] = (char)('0' + exponent % 10);
    double value = 0.0;
    size_t used = 0;
    MT_CHECK(myb_parse_double(text, n, &value, &used) == MYB_OK && used == n);
    hash_double(h, value);
  }
}

static void hash_time(hasher* h, myb_time t) {
  hash_double(h, (double)t.s);
  hash_double(h, t.f);
}

static void workload_time(hasher* h, mt_rng* rng) {
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/time-scales.tk", NULL) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/tdb-tt-fb1990.tk", NULL) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/utc-tai.tk", NULL) == MYB_OK);
  myb_timesys* system = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &system, NULL) == MYB_OK);
  myb_kpool_destroy(pool);
  if (system == NULL) {
    return;
  }
  const char* names[] = {"TAI", "TT", "TCG", "TDB", "TCB", "GPS"};
  myb_scale scales[6];
  for (int i = 0; i < 6; ++i) {
    MT_CHECK(myb_timesys_scale(system, names[i], &scales[i]) == MYB_OK);
  }
  for (int i = 0; i < 5000; ++i) {
    const myb_time t = myb_time_from_seconds(mt_double(rng, 10, 32));
    const int from = mt_int(rng, 0, 5);
    const int to = mt_int(rng, 0, 5);
    myb_time out = {};
    MT_CHECK(myb_timesys_convert(system, t, scales[from], scales[to], &out) == MYB_OK);
    hash_time(h, out);
    myb_calendar utc = {};
    if (myb_timesys_to_utc(system, out, scales[to], &utc) == MYB_OK) {
      hash_double(h, (double)(utc.year * 100000000 + utc.month * 1000000 + utc.day * 10000 +
                              utc.hour * 100 + utc.minute));
      hash_double(h, utc.second);
    }
    double jd1 = 0.0;
    double jd2 = 0.0;
    myb_time_to_jd(out, &jd1, &jd2);
    hash_double(h, jd1);
    hash_double(h, jd2);
  }
  myb_timesys_destroy(system);
}

/* Every other elementary function, through a fixed set of arguments per function. */
static void workload_functions(hasher* h, mt_rng* rng) {
  for (int i = 0; i < 5000; ++i) {
    const double x = mt_double(rng, -20, 9);
    const double y = mt_double(rng, -20, 9);
    const double unit = mt_double(rng, -30, -1);
    hash_double(h, myb_tan(x));
    hash_double(h, myb_atan(x));
    hash_double(h, myb_atan2(y, x));
    hash_double(h, myb_asin(unit));
    hash_double(h, myb_acos(unit));
    hash_double(h, myb_exp(x));
    hash_double(h, myb_expm1(x));
    hash_double(h, myb_log(__builtin_fabs(x)));
    hash_double(h, myb_log1p(__builtin_fabs(x)));
    hash_double(h, myb_sinh(x));
    hash_double(h, myb_cosh(x));
    hash_double(h, myb_tanh(x));
    hash_double(h, myb_asinh(x));
    hash_double(h, myb_acosh(1.0 + __builtin_fabs(x)));
    hash_double(h, myb_atanh(unit));
    hash_double(h, myb_cbrt(x));
    hash_double(h, myb_hypot(x, y));
  }
}

static void hash_pv(hasher* h, const myb_pv* pv) {
  for (int k = 0; k < 3; ++k) {
    hash_double(h, pv->r[k]);
    hash_double(h, pv->v[k]);
  }
}

static void workload_kepler(hasher* h, mt_rng* rng) {
  const double mu = 3.986004418e14;
  for (int i = 0; i < 2000; ++i) {
    myb_elements el = {};
    el.p = 7e6 + 3e7 * (double)(mt_next(rng) >> 11) * 0x1p-53;
    const double e_draw = (double)(mt_next(rng) >> 11) * 0x1p-53;
    el.e = (i % 3 == 0) ? 1.0 + 3.0 * e_draw : 0.99 * e_draw;
    el.i = 3.0 * (double)(mt_next(rng) >> 11) * 0x1p-53;
    el.raan = 6.0 * (double)(mt_next(rng) >> 11) * 0x1p-53;
    el.argp = 6.0 * (double)(mt_next(rng) >> 11) * 0x1p-53;
    el.nu = (el.e > 1.0 ? 1.0 : 3.0) * (2.0 * (double)(mt_next(rng) >> 11) * 0x1p-53 - 1.0);
    myb_pv pv = {};
    MT_CHECK(myb_pv_from_elements(mu, &el, &pv) == MYB_OK);
    hash_pv(h, &pv);
    myb_elements back = {};
    MT_CHECK(myb_elements_from_pv(mu, &pv, &back) == MYB_OK);
    hash_double(h, back.p);
    hash_double(h, back.e);
    hash_double(h, back.nu);
    myb_equinoctial eq = {};
    MT_CHECK(myb_equinoctial_from_pv(mu, &pv, &eq) == MYB_OK);
    hash_double(h, eq.f);
    hash_double(h, eq.l);
    const double dt = mt_double(rng, 0, 14);
    myb_pv later = {};
    MT_CHECK(myb_kepler_propagate(mu, &pv, dt, &later) == MYB_OK);
    hash_pv(h, &later);
    myb_lambert_solution solutions[5];
    size_t count = 0;
    if (myb_lambert(mu, pv.r, later.r, dt, 0, 2, solutions, 5, &count) == MYB_OK) {
      for (size_t j = 0; j < count; ++j) {
        for (int k = 0; k < 3; ++k) {
          hash_double(h, solutions[j].v1[k]);
          hash_double(h, solutions[j].v2[k]);
        }
      }
    }
    double ecc = 0.0;
    MT_CHECK(myb_kepler_elliptic(mt_double(rng, -5, 3), 0.999 * e_draw, &ecc) == MYB_OK);
    hash_double(h, ecc);
  }
}

int main(void) {
  hasher h = {.state = UINT64_C(0xCBF29CE484222325), .values = 0};
  mt_rng rng = {.state = UINT64_C(20261107)};
  workload_cmath(&h, &rng);
  workload_dd(&h, &rng);
  workload_sums(&h, &rng);
  workload_parse(&h, &rng);
  workload_time(&h, &rng);
  workload_functions(&h, &rng);
  workload_kepler(&h, &rng);
  printf("repro hash: 0x%016" PRIx64 " over %ld values (golden 0x%016" PRIx64 ")\n", h.state,
         h.values, golden_hash);
  MT_CHECK(h.state == golden_hash);
  return mt_finish("test_repro");
}
