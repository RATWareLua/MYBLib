/* SPDX-License-Identifier: BSD-3-Clause */
/* Bit-exact reproducibility: a fixed workload over the whole public API is hashed and the hash
 * must be the same on every platform, compiler version, optimisation level and ISA level.
 * The golden value was produced on Windows 10 / Clang 23 / x86-64 (windows-dev preset); Linux
 * builds must reproduce it unchanged. A different hash means that some result changed in at
 * least one bit: find which one before touching the golden value.
 *
 * Inputs come from SplitMix64 and are built from integer bits, so the workload itself does
 * not depend on any library function. */
#include <stdlib.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/cmath.h"
#include "myb/dyn.h"
#include "myb/integ.h"
#include "myb/io.h"
#include "myb/kepler.h"
#include "myb/lambert.h"
#include "myb/num.h"
#include "myb/pool.h"
#include "myb/time.h"

static const uint64_t golden_hash = UINT64_C(0xde3d79ea60027aae);

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
    myb_pv same = {};
    double stm[36];
    MT_CHECK(myb_kepler_stm(mu, &pv, dt, &same, stm) == MYB_OK);
    for (int k = 0; k < 36; ++k) {
      hash_double(h, stm[k]);
    }
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

/* Earth with J2 (DOP853) and a planar three-body system (IAS15). */
static myb_status j2_rhs(void* context, double t, const double* y, double* dydt) {
  (void)context;
  (void)t;
  const double mu = 3.986004418e14;
  const double j2 = 1.08262668e-3;
  const double radius = 6378137.0;
  const double r2 = y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
  const double r = __builtin_sqrt(r2);
  const double k = -mu / (r2 * r);
  const double q = 1.5 * j2 * radius * radius / r2;
  const double z2 = y[2] * y[2] / r2;
  dydt[0] = y[3];
  dydt[1] = y[4];
  dydt[2] = y[5];
  dydt[3] = k * y[0] * (1.0 + q * (1.0 - 5.0 * z2));
  dydt[4] = k * y[1] * (1.0 + q * (1.0 - 5.0 * z2));
  dydt[5] = k * y[2] * (1.0 + q * (3.0 - 5.0 * z2));
  return MYB_OK;
}

static myb_status three_body(void* context, double t, const double* x, const double* v,
                             double* a) {
  (void)context;
  (void)t;
  (void)v;
  const double gm[3] = {1.0, 1e-3, 1e-6};
  for (int i = 0; i < 9; ++i) {
    a[i] = 0.0;
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i == j) {
        continue;
      }
      const double d[3] = {x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1],
                           x[3 * j + 2] - x[3 * i + 2]};
      const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
      const double k = gm[j] / (r2 * __builtin_sqrt(r2));
      for (int c = 0; c < 3; ++c) {
        a[3 * i + c] += k * d[c];
      }
    }
  }
  return MYB_OK;
}

static myb_status three_body_x(void* context, double t, const double* x, double* a) {
  return three_body(context, t, x, NULL, a);
}

static double radial_velocity(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[0] * y[3] + y[1] * y[4] + y[2] * y[5];
}

static void hash_event(hasher* h, int found, double t_event) {
  if (found) {
    hash_double(h, t_event);
  }
}

static void workload_integ(hasher* h) {
  const myb_dop853_options dop_options = {.rtol = 1e-11, .atol = 1e-6};
  myb_dop853* dop = NULL;
  MT_CHECK(myb_dop853_create(6, &dop_options, NULL, &dop) == MYB_OK);
  if (dop != NULL) {
    const double y0[6] = {7.0e6, 0.0, 1.0e5, 0.0, 6.0e3, 4.5e3};
    MT_CHECK(myb_dop853_start(dop, j2_rhs, NULL, 0.0, y0) == MYB_OK);
    double t_previous = 0.0;
    /* bound: 5 periods of about 6000 s; at most 10^5 steps */
    for (int i = 0; i < 100000 && myb_dop853_time(dop) < 30000.0; ++i) {
      MT_CHECK(myb_dop853_step(dop, 30000.0) == MYB_OK);
      double y[6] = {};
      const double t_now = myb_dop853_time(dop);
      MT_CHECK(myb_dop853_dense(dop, t_previous + 0.3 * (t_now - t_previous), y) == MYB_OK);
      t_previous = t_now;
      for (int k = 0; k < 6; ++k) {
        hash_double(h, myb_dop853_state(dop)[k]);
        hash_double(h, y[k]);
      }
      int found = 0;
      double t_event = 0.0;
      MT_CHECK(myb_dop853_locate(dop, radial_velocity, NULL, &found, &t_event, NULL) == MYB_OK);
      hash_event(h, found, t_event);
    }
    const myb_dop853_stats st = myb_dop853_statistics(dop);
    hash_double(h, (double)st.accepted);
    hash_double(h, (double)st.rejected);
    myb_dop853_destroy(dop);
  }
  const myb_ias15_options ias_options = {.block = 3};
  myb_ias15* ias = NULL;
  MT_CHECK(myb_ias15_create(9, &ias_options, NULL, &ias) == MYB_OK);
  if (ias != NULL) {
    /* Star, planet at 1 on a slightly eccentric orbit, test-like body at 1.6 (units G = 1). */
    const double y0[18] = {0.0, 0.0, 0.0, 1.0,  0.0, 0.01, 1.6, 0.0, -0.02,
                           0.0, -1.1e-3, 0.0, 0.0, 1.1, 0.0, 0.0, 0.79, 0.01};
    MT_CHECK(myb_ias15_start(ias, three_body, NULL, 0.0, y0) == MYB_OK);
    const double t_end = 40.0;
    double t_previous = 0.0;
    /* bound: about 6 orbits of the planet; at most 10^5 steps */
    for (int i = 0; i < 100000 && myb_ias15_time(ias) < t_end; ++i) {
      MT_CHECK(myb_ias15_step(ias, t_end) == MYB_OK);
      double y[18] = {};
      const double t_now = myb_ias15_time(ias);
      MT_CHECK(myb_ias15_dense(ias, t_previous + 0.3 * (t_now - t_previous), y) == MYB_OK);
      t_previous = t_now;
      for (int k = 0; k < 18; ++k) {
        hash_double(h, myb_ias15_state(ias)[k]);
        hash_double(h, y[k]);
      }
    }
    const myb_ias15_stats st = myb_ias15_statistics(ias);
    hash_double(h, (double)st.accepted);
    hash_double(h, (double)st.rejected);
    hash_double(h, (double)st.iterations);
    myb_ias15_destroy(ias);
  }
  for (int order = 2; order <= 8; order += 2) {
    myb_symplectic* sym = NULL;
    MT_CHECK(myb_symplectic_create(9, order, NULL, &sym) == MYB_OK);
    if (sym == NULL) {
      continue;
    }
    const double y0[18] = {0.0, 0.0, 0.0, 1.0,  0.0, 0.01, 1.6, 0.0, -0.02,
                           0.0, -1.1e-3, 0.0, 0.0, 1.1, 0.0, 0.0, 0.79, 0.01};
    MT_CHECK(myb_symplectic_start(sym, three_body_x, NULL, 0.0, y0) == MYB_OK);
    for (int k = 1; k <= 10; ++k) {
      MT_CHECK(myb_symplectic_integrate(sym, 1.0 * k, 0.05) == MYB_OK);
      for (int i = 0; i < 18; ++i) {
        hash_double(h, myb_symplectic_state(sym)[i]);
      }
    }
    myb_symplectic_destroy(sym);
  }
  const myb_rkn12_options rkn_options = {.rtol = 1e-12, .atol = 1e-12};
  myb_rkn12* rkn = NULL;
  MT_CHECK(myb_rkn12_create(9, &rkn_options, NULL, &rkn) == MYB_OK);
  if (rkn != NULL) {
    const double y0[18] = {0.0, 0.0, 0.0, 1.0,  0.0, 0.01, 1.6, 0.0, -0.02,
                           0.0, -1.1e-3, 0.0, 0.0, 1.1, 0.0, 0.0, 0.79, 0.01};
    MT_CHECK(myb_rkn12_start(rkn, three_body_x, NULL, 0.0, y0) == MYB_OK);
    double t_previous = 0.0;
    /* bound: about 6 orbits of the planet; at most 10^5 steps */
    for (int i = 0; i < 100000 && myb_rkn12_time(rkn) < 40.0; ++i) {
      MT_CHECK(myb_rkn12_step(rkn, 40.0) == MYB_OK);
      const double t_now = myb_rkn12_time(rkn);
      double y[18] = {};
      if (i % 10 == 0) {
        MT_CHECK(myb_rkn12_dense(rkn, t_previous + 0.3 * (t_now - t_previous), y) == MYB_OK);
      }
      for (int k = 0; k < 18; ++k) {
        hash_double(h, myb_rkn12_state(rkn)[k]);
        hash_double(h, y[k]);
      }
      t_previous = t_now;
    }
    const myb_rkn12_stats st = myb_rkn12_statistics(rkn);
    hash_double(h, (double)st.accepted);
    hash_double(h, (double)st.rejected);
    myb_rkn12_destroy(rkn);
  }
  const myb_gauss_jackson_options gj_options = {.h = 0.02};
  myb_gauss_jackson* gj = NULL;
  MT_CHECK(myb_gauss_jackson_create(9, &gj_options, NULL, &gj) == MYB_OK);
  if (gj != NULL) {
    const double y0[18] = {0.0, 0.0, 0.0, 1.0,  0.0, 0.01, 1.6, 0.0, -0.02,
                           0.0, -1.1e-3, 0.0, 0.0, 1.1, 0.0, 0.0, 0.79, 0.01};
    MT_CHECK(myb_gauss_jackson_start(gj, three_body, NULL, 0.0, y0) == MYB_OK);
    for (int k = 0; k < 500; ++k) {
      MT_CHECK(myb_gauss_jackson_step(gj) == MYB_OK);
      double y[18] = {};
      MT_CHECK(myb_gauss_jackson_dense(gj, myb_gauss_jackson_time(gj) - 0.013, y) == MYB_OK);
      for (int i = 0; i < 18; ++i) {
        hash_double(h, myb_gauss_jackson_state(gj)[i]);
        hash_double(h, y[i]);
      }
    }
    const myb_gauss_jackson_stats st = myb_gauss_jackson_statistics(gj);
    hash_double(h, (double)st.evaluations);
    hash_double(h, st.max_difference);
    myb_gauss_jackson_destroy(gj);
  }
}

/* N-body accelerations with 1 and 4 threads (must agree bit for bit, then hashed once), the
 * invariants, the figure eight by IAS15, and Wisdom-Holman schemes on a two-planet system. */
static void workload_dyn(hasher* h, mt_rng* rng) {
  enum { MASSIVE = 300, TEST = 21 };
  double* gm = malloc(MASSIVE * sizeof(double));
  double* x = malloc(3 * (MASSIVE + TEST) * sizeof(double));
  double* v = malloc(3 * (MASSIVE + TEST) * sizeof(double));
  double* a1 = malloc(3 * (MASSIVE + TEST) * sizeof(double));
  double* a4 = malloc(3 * (MASSIVE + TEST) * sizeof(double));
  MT_CHECK(gm != NULL && x != NULL && v != NULL && a1 != NULL && a4 != NULL);
  for (size_t i = 0; i < MASSIVE; ++i) {
    gm[i] = 1e8 * (1.0 + (double)(mt_next(rng) >> 11) * 0x1p-53);
  }
  for (size_t i = 0; i < 3 * (MASSIVE + TEST); ++i) {
    x[i] = 1e9 * ((double)(mt_next(rng) >> 11) * 0x1p-52 - 1.0);
    v[i] = 1e2 * ((double)(mt_next(rng) >> 11) * 0x1p-52 - 1.0);
  }
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(MASSIVE, TEST, gm, NULL, &s) == MYB_OK);
  myb_pool* pool = NULL;
  MT_CHECK(myb_pool_create(4, NULL, &pool) == MYB_OK);
  if (s != NULL && pool != NULL) {
    MT_CHECK(myb_nbody_set_interaction(s, 7, 3, MYB_INTERACTION_NONE) == MYB_OK);
    MT_CHECK(myb_nbody_acceleration(s, 0.0, x, a1) == MYB_OK);
    myb_nbody_set_pool(s, pool);
    MT_CHECK(myb_nbody_acceleration(s, 0.0, x, a4) == MYB_OK);
    myb_nbody_set_pool(s, NULL);
    int same = 1;
    for (size_t i = 0; i < 3 * (MASSIVE + TEST); ++i) {
      same = same && mt_same(a1[i], a4[i]);
      hash_double(h, a1[i]);
    }
    MT_CHECK(same);
    myb_nbody_invariants inv = {};
    MT_CHECK(myb_nbody_invariants_of(s, x, v, &inv) == MYB_OK);
    hash_dd(h, inv.energy);
    for (int k = 0; k < 3; ++k) {
      hash_double(h, inv.momentum[k]);
      hash_double(h, inv.angular_momentum[k]);
      hash_double(h, inv.center[k]);
    }
  }
  myb_nbody_destroy(s);
  myb_pool_destroy(pool);
  free(gm);
  free(x);
  free(v);
  free(a1);
  free(a4);
  /* Figure eight (start refined by tests/unit/test_nbody.c), one period by IAS15. */
  const double ones[3] = {1.0, 1.0, 1.0};
  MT_CHECK(myb_nbody_create(3, 0, ones, NULL, &s) == MYB_OK);
  const double p1 = 0.3471168881189257;
  const double p2 = 0.5327249453880308;
  const double y0[18] = {-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         p1,   p2,  0.0, p1,  p2,  0.0, -2.0 * p1, -2.0 * p2, 0.0};
  const myb_ias15_options o = {.block = 3};
  myb_ias15* ias = NULL;
  MT_CHECK(myb_ias15_create(9, &o, NULL, &ias) == MYB_OK);
  if (ias != NULL && s != NULL) {
    MT_CHECK(myb_ias15_start(ias, myb_nbody_acceleration_xv, s, 0.0, y0) == MYB_OK);
    MT_CHECK(myb_ias15_integrate(ias, 6.3259139829262061) == MYB_OK);
    for (int k = 0; k < 18; ++k) {
      hash_double(h, myb_ias15_state(ias)[k]);
    }
    hash_double(h, (double)myb_ias15_statistics(ias).accepted);
  }
  myb_ias15_destroy(ias);
  myb_nbody_destroy(s);
  /* Star, two planets, a test body: Wisdom-Holman schemes. */
  const double gm3[3] = {1.0, 9.5e-4, 2.85e-4};
  MT_CHECK(myb_nbody_create(3, 1, gm3, NULL, &s) == MYB_OK);
  const double xs[12] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.02, -1.4, 1.1, 0.0, 0.0, 2.6, 0.1};
  const double vs[12] = {0.0, -1.2e-3, 0.0, 0.0, 1.0, 0.0, -0.45, -0.57, 0.01, -0.62, 0.0, 0.0};
  const myb_wh_options options[3] = {{.h = 0.05, .scheme = MYB_WH_DKD, .corrector = 11},
                                     {.h = 0.08, .scheme = MYB_WH_SABA3},
                                     {.h = 0.08, .scheme = MYB_WH_SBAB2}};
  for (int k = 0; k < 3 && s != NULL; ++k) {
    myb_wh* w = NULL;
    MT_CHECK(myb_wh_create(s, &options[k], NULL, &w) == MYB_OK);
    if (w == NULL) {
      continue;
    }
    MT_CHECK(myb_wh_start(w, 0.0, xs, vs) == MYB_OK);
    MT_CHECK(myb_wh_integrate(w, 40.0) == MYB_OK);
    double xo[12];
    double vo[12];
    MT_CHECK(myb_wh_state(w, xo, vo) == MYB_OK);
    for (int i = 0; i < 12; ++i) {
      hash_double(h, xo[i]);
      hash_double(h, vo[i]);
    }
    MT_CHECK(myb_wh_dense(w, myb_wh_time(w) + 0.3 * options[k].h, xo, vo) == MYB_OK);
    for (int i = 0; i < 12; ++i) {
      hash_double(h, xo[i]);
      hash_double(h, vo[i]);
    }
    hash_double(h, (double)myb_wh_statistics(w).evaluations);
    myb_wh_destroy(w);
  }
  myb_nbody_destroy(s);
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
  workload_integ(&h);
  workload_dyn(&h, &rng);
  printf("repro hash: 0x%016" PRIx64 " over %ld values (golden 0x%016" PRIx64 ")\n", h.state,
         h.values, golden_hash);
  MT_CHECK(h.state == golden_hash);
  return mt_finish("test_repro");
}
