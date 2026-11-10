/* SPDX-License-Identifier: BSD-3-Clause */
/* Bit-exact reproducibility: a fixed workload over the whole public API is hashed and the hash
 * must be the same on every platform, compiler version, optimisation level and ISA level.
 * The golden value was produced on Windows 10 / Clang 23 / x86-64 (windows-dev preset); Linux
 * builds must reproduce it unchanged. A different hash means that some result changed in at
 * least one bit: find which one before touching the golden value.
 *
 * Inputs come from SplitMix64 and are built from integer bits, so the workload itself does
 * not depend on any library function. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/attitude.h"
#include "myb/budget.h"
#include "myb/cmath.h"
#include "myb/cr3bp.h"
#include "myb/dyn.h"
#include "myb/force.h"
#include "myb/frames.h"
#include "myb/grav.h"
#include "myb/integ.h"
#include "myb/io.h"
#include "myb/kepler.h"
#include "myb/kin.h"
#include "myb/lambert.h"
#include "myb/model.h"
#include "myb/num.h"
#include "myb/pool.h"
#include "myb/regular.h"
#include "myb/rel.h"
#include "myb/thrust.h"
#include "myb/time.h"
#include "myb/transfer.h"
#include "myb/uncert.h"
#include "myb_sgp4.h"

static const uint64_t golden_hash = UINT64_C(0x827497389931b8a9);

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

/* Phase 5: every SPK type (tests/data/spk), the IAU 2006/2000A CIP model, ERA, EOP with tidal
 * corrections, WGCCRE rotation, geodetic coordinates and orbital frames. */
static void workload_frames(hasher* h, mt_rng* rng) {
  const char* files[] = {"t01.bsp", "t02_03.bsp", "t05.bsp", "t08_09_12_13.bsp", "t14.bsp",
                         "t17.bsp", "t18_19.bsp", "t20.bsp", "eros.bsp"};
  myb_spk* spk = NULL;
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  for (size_t i = 0; i < sizeof files / sizeof files[0]; ++i) {
    char path[512];
    snprintf(path, sizeof path, "%s/spk/%s", MYB_TEST_DATA_DIR, files[i]);
    MT_CHECK(myb_spk_load_file(spk, path, NULL) == MYB_OK);
  }
  for (size_t i = 0; i < myb_spk_segment_count(spk); ++i) {
    myb_spk_segment_info info;
    myb_spk_segment(spk, i, &info);
    for (int k = 0; k < 8; ++k) {
      const double u = (double)mt_int(rng, 0, 1 << 30) / (double)(1 << 30);
      const double t = info.start + u * (info.stop - info.start);
      double state[6];
      MT_CHECK(myb_spk_state(spk, info.target, info.center, myb_time_from_seconds(t), state,
                             NULL) == MYB_OK);
      for (int c = 0; c < 6; ++c) {
        hash_double(h, state[c]);
      }
    }
  }
  myb_spk_destroy(spk);

  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const char* kernels[] = {"/time/time-scales.tk", "/time/tdb-tt-fb1990.tk", "/time/utc-tai.tk",
                           "/frames/iau2006-cip.tk", "/frames/eop-tides.tk"};
  for (size_t i = 0; i < sizeof kernels / sizeof kernels[0]; ++i) {
    char path[512];
    snprintf(path, sizeof path, "%s%s", MYB_DATA_DIR, kernels[i]);
    MT_CHECK(myb_kpool_load_file(pool, path, NULL) == MYB_OK);
  }
  MT_CHECK(myb_kpool_load_file(pool, MYB_TEST_DATA_DIR "/wgccre/pck00011.tpc", NULL) == MYB_OK);
  myb_timesys* sys = NULL;
  myb_cip* cip = NULL;
  myb_eop* eop = NULL;
  myb_wgccre* mars = NULL;
  myb_scale tt = {0};
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, NULL) == MYB_OK);
  MT_CHECK(myb_timesys_scale(sys, "TT", &tt) == MYB_OK);
  MT_CHECK(myb_cip_create(pool, "IAU2006A", NULL, &cip, NULL) == MYB_OK);
  MT_CHECK(myb_wgccre_create(pool, 499, NULL, &mars, NULL) == MYB_OK);
  MT_CHECK(myb_eop_create(NULL, &eop) == MYB_OK);
  MT_CHECK(myb_eop_load_c04_file(eop, MYB_TEST_DATA_DIR "/eop/eopc04-excerpt.txt", NULL) ==
           MYB_OK);
  MT_CHECK(myb_eop_set_tides(eop, pool, "IERS10_EOPT", NULL) == MYB_OK);
  const myb_ellipsoid wgs84 = {6378137.0, 1.0 / 298.257223563};
  for (int k = 0; k < 40; ++k) {
    /* 2020-01-03 .. 2020-12-28 TT */
    const myb_time t = {631411200 + (int64_t)mt_int(rng, 0, 31000000), 0.125};
    myb_eop_values v;
    myb_time ut1;
    MT_CHECK(myb_eop_at(eop, sys, t, tt, &v, &ut1) == MYB_OK);
    hash_double(h, v.xp);
    hash_double(h, v.dut1);
    hash_double(h, v.lod);
    myb_rotation m;
    myb_cip_terrestrial(cip, t, ut1, &v, &m);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, m.m[i]);
      hash_double(h, m.dm[i]);
    }
    myb_wgccre_rotation(mars, t, &m);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, m.m[i]);
      hash_double(h, m.dm[i]);
    }
    const double xyz[3] = {mt_double(rng, 20, 25), -mt_double(rng, 18, 24),
                           mt_double(rng, 10, 24)};
    double lon;
    double lat;
    double height;
    MT_CHECK(myb_geodetic_from_cartesian(&wgs84, xyz, &lon, &lat, &height) == MYB_OK);
    hash_double(h, lon);
    hash_double(h, lat);
    hash_double(h, height);
    const double state[6] = {xyz[0], xyz[1], xyz[2], 1.1e3, 7.4e3, -0.6e3};
    const double acc[3] = {-1.0, 0.5, 0.25};
    MT_CHECK(myb_orbital_rotation(MYB_ORBITAL_LVLH, state, acc, &m) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, m.dm[i]);
    }
  }
  myb_eop_destroy(eop);
  myb_wgccre_destroy(mars);
  myb_cip_destroy(cip);
  myb_timesys_destroy(sys);
  myb_kpool_destroy(pool);
}

/* Phase 6: EIH accelerations (1 and 4 threads compared bit for bit), 1PN invariants, the
 * satellite terms of IERS 2010 eq. 10.12, the BCRS -> GCRS transformation. */
static void workload_rel(hasher* h, mt_rng* rng) {
  enum { M = 40, T = 5, N = M + T };
  static double gm[M];
  static double x[3 * N];
  static double v[3 * N];
  static double a1[3 * N];
  static double a4[3 * N];
  for (int i = 0; i < M; ++i) {
    gm[i] = fabs(mt_double(rng, 10, 20));
  }
  for (int i = 0; i < 3 * N; ++i) {
    x[i] = mt_double(rng, 30, 34) * (mt_int(rng, 0, 1) ? 1.0 : -1.0);
    v[i] = mt_double(rng, 10, 14) * (mt_int(rng, 0, 1) ? 1.0 : -1.0);
  }
  myb_nbody* s = NULL;
  myb_pool* pool = NULL;
  MT_CHECK(myb_nbody_create(M, T, gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_pool_create(4, NULL, &pool) == MYB_OK);
  const myb_ppn ppn = {299792458.0, 1.0, 1.0};
  MT_CHECK(myb_nbody_set_relativity(s, &ppn) == MYB_OK);
  MT_CHECK(myb_nbody_set_interaction(s, 7, 3, MYB_INTERACTION_NEWTON) == MYB_OK);
  MT_CHECK(myb_nbody_set_interaction(s, 3, 7, MYB_INTERACTION_NONE) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration_xv(s, 0.0, x, v, a1) == MYB_OK);
  myb_nbody_set_pool(s, pool);
  MT_CHECK(myb_nbody_acceleration_xv(s, 0.0, x, v, a4) == MYB_OK);
  for (int i = 0; i < 3 * N; ++i) {
    MT_CHECK(mt_bits(a1[i]) == mt_bits(a4[i]));
    hash_double(h, a1[i]);
  }
  myb_nbody_pn_invariants inv;
  MT_CHECK(myb_nbody_pn_invariants_of(s, x, v, &inv) == MYB_OK);
  hash_dd(h, inv.energy);
  for (int k = 0; k < 3; ++k) {
    hash_double(h, inv.momentum[k]);
    hash_double(h, inv.center[k]);
  }
  myb_nbody_destroy(s);
  myb_pool_destroy(pool);
  const double gm_e = 3.986004418e14;
  const double j[3] = {0.0, 0.0, 9.8e8};
  myb_rel_local_frame f = {.c = 299792458.0, .xe = {-2.6e10, 1.3e11, 5.7e10},
                           .ve = {-2.98e4, -5.0e3, -2.2e3}, .ae = {1.5e-3, -5.6e-3, -2.4e-3},
                           .w_ext = 8.9e8};
  for (int k = 0; k < 50; ++k) {
    const double r[3] = {mt_double(rng, 22, 25), -mt_double(rng, 22, 25), mt_double(rng, 20, 24)};
    const double vel[3] = {mt_double(rng, 11, 13), mt_double(rng, 11, 13), -mt_double(rng, 9, 12)};
    double out[3];
    MT_CHECK(myb_rel_schwarzschild(&ppn, gm_e, r, vel, out) == MYB_OK);
    hash_double(h, out[0]);
    hash_double(h, out[2]);
    MT_CHECK(myb_rel_lense_thirring(&ppn, gm_e, j, r, vel, out) == MYB_OK);
    hash_double(h, out[1]);
    MT_CHECK(myb_rel_de_sitter(&ppn, 1.32712440041279419e20, f.xe, f.ve, vel, out) == MYB_OK);
    hash_double(h, out[0]);
    const double xb[3] = {f.xe[0] + r[0], f.xe[1] + r[1], f.xe[2] + r[2]};
    myb_rel_global_to_local(&f, xb, out);
    hash_double(h, out[0] - r[0]);
    myb_rel_local_to_global(&f, out, out);
    hash_double(h, out[1]);
    hash_double(h, myb_rel_local_time_shift(&f, xb));
  }
}

/* Phase 7: gravity fields (spherical and ellipsoidal harmonics, polyhedra, mascons, two
 * figures), tides, forces, SP3. */
static void workload_grav(hasher* h, mt_rng* rng) {
  myb_diag diag = {0};
  /* spherical harmonics, degree 40 with random coefficients */
  enum { N = 40 };
  static double c[(N + 1) * (N + 2) / 2];
  static double s[(N + 1) * (N + 2) / 2];
  for (int k = 0; k < (N + 1) * (N + 2) / 2; ++k) {
    c[k] = mt_double(rng, -30, -12) * (mt_int(rng, 0, 1) ? 1.0 : -1.0);
    s[k] = mt_double(rng, -30, -12) * (mt_int(rng, 0, 1) ? 1.0 : -1.0);
  }
  c[0] = 1.0;
  myb_gravity_field* f = NULL;
  myb_gravity_work* w = NULL;
  MT_CHECK(myb_gravity_field_create(3.986004415e14, 6378136.3, N, c, s, NULL, &f) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(N, NULL, &w) == MYB_OK);
  for (int k = 0; k < 20; ++k) {
    double r[3];
    for (int i = 0; i < 3; ++i) {
      r[i] = mt_double(rng, 22, 23) * (mt_int(rng, 0, 1) ? 1.0 : -1.0);
    }
    double v = 0.0;
    double a[3];
    MT_CHECK(myb_gravity_field_eval(f, w, r, k % 3, -1, &v, a) == MYB_OK);
    hash_double(h, v);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, a[i]);
    }
  }
  /* time-variable terms (ICGEM 2.0) and tides of IERS 2010 on a degree-4 field */
  myb_gravity_field* tv = NULL;
  MT_CHECK(myb_gravity_field_load_icgem_file(MYB_TEST_DATA_DIR "/grav/synthetic2.gfc", 0, NULL,
                                             &tv, &diag) == MYB_OK);
  myb_kpool* pool = NULL;
  myb_tides* tides = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/grav/iers2010-tides.tk", &diag) == MYB_OK);
  MT_CHECK(myb_tides_create(pool, "IERS10_TIDE", NULL, &tides, &diag) == MYB_OK);
  MT_CHECK(myb_tides_load_ocean_file(tides, MYB_TEST_DATA_DIR "/grav/fes2004_excerpt.dat", 3,
                                     &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  for (int k = 0; k < 4; ++k) {
    const myb_time t = myb_time_from_seconds(8.0e8 + 3.1e7 * k);
    myb_gravity_field_set_time(tv, t);
    const myb_tide_raiser raisers[2] = {{4.9028e12, {3.0e8, -2.2e8, 1.1e8 * k}},
                                        {1.32712440018e20, {1.4e11, 4.0e10, -2.0e10}}};
    const myb_tides_input in = {.raisers = raisers, .count = 2, .t = t, .xp = 7e-7,
                                .yp = 1.6e-6, .flags = MYB_TIDES_ALL};
    MT_CHECK(myb_tides_apply(tides, tv, &in) == MYB_OK);
    for (int n = 0; n <= 3; ++n) {
      for (int m = 0; m <= n; ++m) {
        double cc = 0.0;
        double ss = 0.0;
        MT_CHECK(myb_gravity_field_coefficient(tv, n, m, &cc, &ss) == MYB_OK);
        hash_double(h, cc);
        hash_double(h, ss);
      }
    }
  }
  myb_tides_destroy(tides);
  myb_gravity_field_destroy(tv);
  /* polyhedron, mascons, two figures */
  double vert[8][3];
  for (int k = 0; k < 8; ++k) {
    vert[k][0] = (k & 1) ? 2.0 : -1.0;
    vert[k][1] = (k & 2) ? 1.5 : -0.5;
    vert[k][2] = (k & 4) ? 0.7 : -2.0;
  }
  static const size_t faces[12][3] = {{0, 2, 3}, {0, 3, 1}, {4, 5, 7}, {4, 7, 6},
                                      {0, 1, 5}, {0, 5, 4}, {2, 6, 7}, {2, 7, 3},
                                      {0, 4, 6}, {0, 6, 2}, {1, 3, 7}, {1, 7, 5}};
  myb_polyhedron* poly = NULL;
  MT_CHECK(myb_polyhedron_create(&vert[0][0], 8, &faces[0][0], 12, 6.674e-8, NULL, &poly) ==
           MYB_OK);
  double gm[8];
  for (int k = 0; k < 8; ++k) {
    gm[k] = fabs(mt_double(rng, -2, 2));
  }
  for (int k = 0; k < 10; ++k) {
    double r[3];
    for (int i = 0; i < 3; ++i) {
      r[i] = mt_double(rng, 0, 3) * (mt_int(rng, 0, 1) ? 1.0 : -1.0);
    }
    double v = 0.0;
    double a[3];
    double lap = 0.0;
    MT_CHECK(myb_polyhedron_eval(poly, r, &v, a, &lap) == MYB_OK);
    hash_double(h, v);
    hash_double(h, a[k % 3]);
    hash_double(h, lap);
    MT_CHECK(myb_mascons_eval(gm, &vert[0][0], 8, r, &v, a) == MYB_OK);
    hash_double(h, v);
    hash_double(h, a[2]);
  }
  myb_polyhedron_destroy(poly);
  double qa[9];
  double qb[9];
  myb_quadrupole_from_coefficients(3.986e14, 6.378e6, (const double[3]){-4.8e-4, 1e-10, 2.4e-6},
                                   (const double[2]){1.4e-9, -1.4e-6}, qa);
  myb_quadrupole_from_coefficients(4.9e12, 1.738e6, (const double[3]){-9.1e-5, 0.0, 3.5e-5},
                                   (const double[2]){0.0, 0.0}, qb);
  const double rot[9] = {0.8, 0.6, 0.0, -0.6, 0.8, 0.0, 0.0, 0.0, 1.0};
  double qr[9];
  myb_quadrupole_rotate(rot, qb, qr);
  double aa[3] = {0.0, 0.0, 0.0};
  double ab[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_figure_figure(3.986e14, qa, 4.9e12, qr, (const double[3]){3.1e8, -2.2e8, 0.4e8},
                             aa, ab) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    hash_double(h, aa[i]);
    hash_double(h, ab[i]);
  }
  /* ellipsoidal harmonics, degree 6 */
  double alpha[49];
  for (int k = 0; k < 49; ++k) {
    alpha[k] = mt_double(rng, -12, -2) * (mt_int(rng, 0, 1) ? 1.0 : -1.0);
  }
  alpha[0] = 1.0;
  myb_ellipsoidal_field* e = NULL;
  MT_CHECK(myb_ellipsoidal_field_create(4.46e5, (const double[3]){20e3, 7e3, 5.5e3}, 6, alpha,
                                        NULL, &e) == MYB_OK);
  for (int k = 0; k < 8; ++k) {
    const double r[3] = {(21e3 + 1e3 * k) * (k % 2 ? 1.0 : -1.0), 3e3 * k, -2e3 + 500.0 * k};
    double v = 0.0;
    double a[3];
    MT_CHECK(myb_ellipsoidal_field_eval(e, r, &v, a) == MYB_OK);
    hash_double(h, v);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, a[i]);
    }
  }
  myb_ellipsoidal_field_destroy(e);
  /* forces */
  for (int k = 0; k < 10; ++k) {
    double r[3];
    double s3[3];
    for (int i = 0; i < 3; ++i) {
      r[i] = mt_double(rng, 22, 24) * (mt_int(rng, 0, 1) ? 1.0 : -1.0);
      s3[i] = mt_double(rng, 36, 37) * (mt_int(rng, 0, 1) ? 1.0 : -1.0);
    }
    double a[3] = {0.0, 0.0, 0.0};
    myb_force_third_body(1.32712440018e20, s3, r, a);
    myb_force_drag(1e-13, 0.02, r, a);
    MT_CHECK(myb_force_empirical_rtn(r, s3, (const double[3]){1e-9, -2e-9, 3e-9}, a) == MYB_OK);
    const double occ[3] = {-r[0], -r[1], -r[2]};
    const double src[3] = {s3[0] - r[0], s3[1] - r[1], s3[2] - r[2]};
    hash_double(h, myb_shadow_fraction(src, 6.957e8, occ, 6.378e6, (const double[2]){0.6, 0.1}));
    const myb_plate plates[2] = {{1.5, {0.0, 0.0, 1.0}, 0.3, 0.2, 1},
                                 {0.7, {0.6, 0.8, 0.0}, 0.1, 0.5, 0}};
    const double nr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    const double u[3] = {-r[0] / nr, -r[1] / nr, -r[2] / nr};
    MT_CHECK(myb_force_radiation_plates(4.56e-6, 500.0, plates, 2, u, a) == MYB_OK);
    const myb_planet_radiation knocke = {.radius = 6378137.0,
                                         .albedo = {0.34, 0.0, 0.10, 0.0, 0.29},
                                         .emissivity = {0.68, 0.0, -0.07, 0.0, -0.18},
                                         .omega = 1.99e-7};
    const double rr[3] = {r[0], r[1], r[2] + (r[2] < 0.0 ? -7e6 : 7e6)};
    MT_CHECK(myb_force_planet_radiation(&knocke, rr, s3, 1361.0, 0.5, 7.8e-4, 3 + k, a) ==
             MYB_OK);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, a[i]);
    }
  }
  /* SP3 */
  myb_sp3* sp3 = NULL;
  MT_CHECK(myb_sp3_load_file(MYB_TEST_DATA_DIR "/sp3/ilrsa.lageos1.250607.excerpt.sp3", NULL, &sp3,
                             &diag) == MYB_OK);
  if (sp3 != NULL) {
    for (size_t k = 0; k < myb_sp3_epoch_count(sp3); k += 97) {
      double p[3];
      double v[3];
      MT_CHECK(myb_sp3_state(sp3, 0, k, p, v) == 3);
      hash_double(h, p[k % 3]);
      hash_double(h, v[(k + 1) % 3]);
    }
    myb_sp3_destroy(sp3);
  }
  myb_gravity_work_destroy(w);
  myb_gravity_field_destroy(f);
}

/* Phase 8: linear algebra, force models and variational equations, random numbers, the
 * restricted three-body problem, B-plane, shooting, effect budgets. */
typedef struct repro_central {
  double gm;
} repro_central;

static myb_status repro_central_accel(void* c, double t, const double x[3], const double v[3],
                                      double a[3]) {
  (void)t;
  (void)v;
  const double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  const double k = -((const repro_central*)c)->gm / (r * r * r);
  for (int i = 0; i < 3; ++i) {
    a[i] = k * x[i];
  }
  return MYB_OK;
}

typedef struct repro_push {
  double k;
} repro_push;

static myb_status repro_push_accel(void* c, double t, const double x[3], const double v[3],
                                   double a[3]) {
  (void)t;
  (void)x;
  const double s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  const double k = ((const repro_push*)c)->k;
  a[0] = k - 1e-9 * s * v[0];
  a[1] = -0.5 * k - 1e-9 * s * v[1];
  a[2] = -1e-9 * s * v[2];
  return MYB_OK;
}

typedef struct repro_shot {
  double r1[3];
  double r2[3];
} repro_shot;

static myb_status repro_shot_f(void* c, const double* u, double* f, double* jac) {
  const repro_shot* s = c;
  const myb_pv p0 = {{s->r1[0], s->r1[1], s->r1[2]}, {u[0], u[1], u[2]}};
  myb_pv p;
  double phi[36];
  const myb_status st = myb_kepler_stm(3.986004418e14, &p0, 4000.0, &p, phi);
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 3; ++i) {
    f[i] = p.r[i] - s->r2[i];
    if (jac != NULL) {
      for (int j = 0; j < 3; ++j) {
        jac[3 * i + j] = phi[i * 6 + 3 + j];
      }
    }
  }
  return MYB_OK;
}

static void workload_phase8(hasher* h, mt_rng* rng) {
  /* linear algebra */
  double a[64];
  double b[8];
  for (int i = 0; i < 64; ++i) {
    a[i] = mt_double(rng, -2, 2) * (mt_int(rng, 0, 1) ? 1.0 : -1.0) + (i % 9 == 0 ? 4.0 : 0.0);
  }
  for (int i = 0; i < 8; ++i) {
    b[i] = mt_double(rng, -2, 2);
  }
  double e[64];
  memcpy(e, a, sizeof e);
  MT_CHECK(myb_lu_solve(8, a, b, 1) == MYB_OK);
  double re[8];
  double im[8];
  MT_CHECK(myb_eigenvalues(8, e, re, im) == MYB_OK);
  for (int i = 0; i < 8; ++i) {
    hash_double(h, b[i]);
    hash_double(h, re[i] * re[i] + im[i] * im[i]); /* order-independent per pair... */
  }
  /* random numbers */
  myb_rng r;
  myb_rng_seed(&r, 20261110);
  for (int i = 0; i < 16; ++i) {
    hash_double(h, myb_rng_normal(&r));
  }
  /* variational equations with a parameter */
  repro_central c = {3.986004418e14};
  repro_push p = {3e-7};
  myb_force_model* m = NULL;
  MT_CHECK(myb_force_model_create(NULL, &m) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "earth", repro_central_accel, NULL, &c, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "push", repro_push_accel, NULL, &p, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(m, "k", &p.k, 1e-9, NULL) == MYB_OK);
  myb_variational* var = NULL;
  MT_CHECK(myb_variational_create(m, NULL, NULL, &var) == MYB_OK);
  const double y0[6] = {7.1e6, -2e5, 1e5, 100.0, 7.4e3, 1.2e3};
  MT_CHECK(myb_variational_start(var, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_variational_integrate(var, 9000.0) == MYB_OK);
  double y[6];
  double phi[36];
  double sens[6];
  myb_variational_get(var, y, phi, sens);
  for (int i = 0; i < 6; ++i) {
    hash_double(h, y[i]);
    hash_double(h, phi[i * 7]);
    hash_double(h, sens[i]);
  }
  myb_variational_destroy(var);
  /* budget */
  myb_budget_entry entries[2];
  const myb_budget_options bo = {.tolerance = 1.0, .samples = 8};
  MT_CHECK(myb_budget(m, 0.0, y0, 7200.0, &bo, entries) == MYB_OK);
  hash_double(h, entries[1].effect);
  hash_double(h, entries[1].max_acceleration);
  myb_force_model_destroy(m);
  /* restricted three-body problem */
  const double mu = 1.215058560962404e-02;
  double l[5][3];
  MT_CHECK(myb_cr3bp_lagrange(mu, l) == MYB_OK);
  for (int k = 0; k < 3; ++k) {
    hash_double(h, l[k][0]);
  }
  double halo[6];
  double period = 0.0;
  MT_CHECK(myb_cr3bp_richardson(mu, 2, 0.02, 1, halo, &period) == MYB_OK);
  hash_double(h, halo[0]);
  hash_double(h, period);
  const myb_cr3bp_correction cc = {.symmetry = MYB_CR3BP_SYMMETRY_XZ, .fixed = 2};
  MT_CHECK(myb_cr3bp_correct(mu, &cc, halo, &period) == MYB_OK);
  hash_double(h, halo[0]);
  hash_double(h, halo[4]);
  hash_double(h, period);
  double mono[36];
  double index = 0.0;
  MT_CHECK(myb_cr3bp_monodromy(mu, halo, period, mono) == MYB_OK);
  MT_CHECK(myb_cr3bp_stability(mono, re, im, &index) == MYB_OK);
  hash_double(h, index);
  hash_double(h, myb_cr3bp_jacobi(mu, halo));
  /* B-plane and shooting */
  const double rr[3] = {-5e8, 3e7, 1e7};
  const double vv[3] = {5000.0, 100.0, -40.0};
  const double pole[3] = {0.0, 0.0, 1.0};
  myb_bplane bp;
  MT_CHECK(myb_bplane_compute(3.986004418e14, rr, vv, pole, &bp) == MYB_OK);
  hash_double(h, bp.b_t);
  hash_double(h, bp.b_r);
  hash_double(h, bp.time_to_periapsis);
  repro_shot shot = {{7e6, 1e6, -2e5}, {-3e6, 9e6, 2e6}};
  double u[3] = {0.0, 7000.0, 1000.0};
  const myb_shoot_options so = {.tolerance = 1e-6};
  double shoot_work[64];
  MT_CHECK(myb_shoot(repro_shot_f, &shot, 3, 3, u, &so, shoot_work, NULL, NULL) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    hash_double(h, u[i]);
  }
}

/* ---- Phase 9 ------------------------------------------------------------------------------ */

static double repro_uniform(mt_rng* rng) { /* [0, 1) from integer bits */
  return (double)(mt_next(rng) >> 11) * 0x1p-53;
}

static myb_status repro_gravity_gradient(void* context, double t, const double q[4],
                                         double torque[3]) {
  (void)context;
  (void)t;
  static const double moments[3] = {1.0, 1.6, 2.1};
  static const double r[3] = {2.0, 1.0, -0.5};
  double rb[3];
  myb_quat_rotate(q, r, rb);
  myb_torque_point_mass(0.8, rb, moments, torque);
  return MYB_OK;
}

static myb_status repro_j2(void* context, double t, const double x[3], const double v[3],
                           double a[3]) {
  (void)context;
  (void)t;
  (void)v;
  const double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  const double r = sqrt(r2);
  const double z2 = x[2] * x[2] / r2;
  const double f = -1.5e-3 / (r2 * r2 * r);
  a[0] = f * x[0] * (1.0 - 5.0 * z2);
  a[1] = f * x[1] * (1.0 - 5.0 * z2);
  a[2] = f * x[2] * (3.0 - 5.0 * z2);
  return MYB_OK;
}

static void workload_phase9(hasher* h, mt_rng* rng) {
  /* equinox-based frames */
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const char* kernels[] = {"/frames/iau1980.tk", "/time/time-scales.tk",
                           "/time/tdb-tt-fb1990.tk", "/time/utc-tai.tk"};
  for (size_t i = 0; i < sizeof kernels / sizeof kernels[0]; ++i) {
    char path[512];
    snprintf(path, sizeof path, "%s%s", MYB_DATA_DIR, kernels[i]);
    MT_CHECK(myb_kpool_load_file(pool, path, NULL) == MYB_OK);
  }
  myb_equinox* eq = NULL;
  myb_timesys* sys = NULL;
  MT_CHECK(myb_equinox_create(pool, "IAU1980", NULL, &eq, NULL) == MYB_OK);
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, NULL) == MYB_OK);
  myb_kpool_destroy(pool);
  for (int k = 0; k < 20; ++k) {
    const myb_time t = {(int64_t)mt_int(rng, -1000000, 1000000) * 3011, repro_uniform(rng)};
    for (int f = 0; f < 3; ++f) {
      myb_rotation r;
      myb_equinox_rotation(eq, (myb_equinox_frame)f, t, &r);
      for (int i = 0; i < 9; ++i) {
        hash_double(h, r.m[i]);
        hash_double(h, r.dm[i]);
      }
    }
    double rate;
    hash_double(h, myb_equinox_gmst(eq, t, &rate));
    hash_double(h, rate);
  }

  /* elliptic integrals and Jacobi functions */
  for (int k = 0; k < 40; ++k) {
    const double x = repro_uniform(rng) * 10.0;
    const double y = repro_uniform(rng) * 10.0 + 0.1;
    const double z = repro_uniform(rng) * 10.0 + 0.1;
    const double p = repro_uniform(rng) * 10.0 + 0.1;
    const double m = repro_uniform(rng);
    const double n = repro_uniform(rng) * 4.0 - 3.0;
    const double phi = repro_uniform(rng) * 20.0 - 10.0;
    hash_double(h, myb_carlson_rf(x, y, z));
    hash_double(h, myb_carlson_rd(x, y, z));
    hash_double(h, myb_carlson_rj(x, y, z, p));
    hash_double(h, myb_carlson_rc(x, y));
    hash_double(h, myb_elliptic_k(m));
    hash_double(h, myb_elliptic_e(m));
    hash_double(h, myb_elliptic_pi(n, m));
    hash_double(h, myb_elliptic_f_inc(phi, m));
    hash_double(h, myb_elliptic_e_inc(phi, m));
    hash_double(h, myb_elliptic_pi_inc(n, phi, m));
    double am;
    double sn;
    double cn;
    double dn;
    myb_jacobi(phi * 7.0, m, &am, &sn, &cn, &dn);
    hash_double(h, am);
    hash_double(h, sn);
    hash_double(h, cn);
    hash_double(h, dn);
  }

  /* attitude */
  for (int k = 0; k < 12; ++k) {
    double q[4];
    for (int i = 0; i < 4; ++i) {
      q[i] = repro_uniform(rng) - 0.5;
    }
    myb_quat_normalize(q);
    double m[9];
    myb_quat_to_matrix(q, m);
    double back[4];
    myb_quat_from_matrix(m, back);
    double phi[3];
    myb_quat_to_rotation_vector(q, phi);
    for (int i = 0; i < 4; ++i) {
      hash_double(h, back[i]);
    }
    for (int i = 0; i < 3; ++i) {
      hash_double(h, phi[i]);
    }
    const int axes[3] = {1 + k % 3, 1 + (k + 1) % 3, 1 + (k % 2 == 0 ? k % 3 : (k + 2) % 3)};
    double ang[3];
    MT_CHECK(myb_euler_from_matrix(axes, m, ang) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, ang[i]);
    }
    const double moments[3] = {1.0 + repro_uniform(rng), 2.0 + repro_uniform(rng),
                               3.0 + repro_uniform(rng)};
    const double l[3] = {repro_uniform(rng) - 0.5, repro_uniform(rng) - 0.5, 1.0};
    myb_free_rotation f;
    if (myb_free_rotation_init(moments, q, l, &f) == MYB_OK) {
      double qt[4];
      double lt[3];
      myb_free_rotation_state(&f, 17.0 * repro_uniform(rng), qt, lt);
      for (int i = 0; i < 4; ++i) {
        hash_double(h, qt[i]);
      }
      for (int i = 0; i < 3; ++i) {
        hash_double(h, lt[i]);
      }
    }
    const double in[9] = {moments[0], 0.1, 0.05, 0.1, moments[1], -0.02, 0.05, -0.02, moments[2]};
    double pm[3];
    double axes9[9];
    MT_CHECK(myb_principal_axes(in, pm, axes9) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, pm[i]);
    }
  }
  {
    const double moments[3] = {1.0, 1.6, 2.1};
    double q[4] = {0.9, 0.1, -0.3, 0.2};
    myb_quat_normalize(q);
    double l[3] = {0.3, 0.5, 1.1};
    MT_CHECK(myb_rotation_split(moments, repro_gravity_gradient, NULL, 6, 0.0, 0.1, 30, q, l) ==
             MYB_OK);
    for (int i = 0; i < 4; ++i) {
      hash_double(h, q[i]);
    }
    for (int i = 0; i < 3; ++i) {
      hash_double(h, l[i]);
    }
    const double im[9] = {0.39, 1e-5, 0, 1e-5, 0.3902, 0, 0, 0, 0.3904};
    const double core[2] = {1.6e-4, 1.6e-4 * 1.0002};
    const double wm[3] = {1e-7, -2e-8, 2.66e-6};
    const double wc[3] = {-7e-8, -1e-8, 2.66e-6};
    const double tq[3] = {1e-17, 2e-17, -3e-18};
    double dwm[3];
    double dwc[3];
    MT_CHECK(myb_core_mantle_rates(im, NULL, core, 4e-14, wm, wc, tq, dwm, dwc) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, dwm[i]);
      hash_double(h, dwc[i]);
    }
  }

  /* regularisation */
  {
    const myb_regular_problem p = {.mu = 1.0, .perturbation = repro_j2, .alpha = 1.5};
    const myb_dop853_options o = {.rtol = 1e-12, .atol = 1e-14};
    const double x0[3] = {1.8, 0.2, 0.1};
    const double v0[3] = {-0.05, 0.12, 0.08};
    for (int kind = 0; kind < 2; ++kind) {
      double x[3];
      double v[3];
      MT_CHECK(myb_regular_propagate((myb_regular_kind)kind, &p, &o, 0.0, x0, v0, 40.0, x, v,
                                     NULL) == MYB_OK);
      for (int i = 0; i < 3; ++i) {
        hash_double(h, x[i]);
        hash_double(h, v[i]);
      }
    }
  }

  /* thrust */
  {
    const myb_thrust_arc arcs[2] = {
        {.start = 0.0, .stop = 5000.0, .thrust = 0.3, .exhaust = 28000.0,
         .steering = MYB_STEER_RTN, .direction = {0.1, 1.0, 0.2}},
        {.start = 6000.0, .stop = 9000.0, .thrust = 0.2, .exhaust = 30000.0,
         .steering = MYB_STEER_EDELBAUM, .edelbaum_v0 = 7546.0, .edelbaum_beta0 = 1.08,
         .edelbaum_sense = -1.0}};
    myb_thrust th;
    MT_CHECK(myb_thrust_init(&th, 120.0, arcs, 2) == MYB_OK);
    const double x[3] = {7000e3, 100e3, 50e3};
    const double v[3] = {-10.0, 6600.0, 3600.0};
    for (int k = 0; k < 10; ++k) {
      const double t = 1000.0 * k;
      double a[3];
      MT_CHECK(myb_thrust_acceleration(&th, t, x, v, a) == MYB_OK);
      for (int i = 0; i < 3; ++i) {
        hash_double(h, a[i]);
      }
      hash_double(h, myb_thrust_mass(&th, t));
    }
    const myb_equinoctial e = {7e6, 0.01, -0.02, 0.2, 0.1, 1.3};
    const double artn[3] = {1e-4, 2e-3, -5e-4};
    double rates[6];
    myb_equinoctial_rates(3.986004418e14, &e, artn, rates);
    for (int i = 0; i < 6; ++i) {
      hash_double(h, rates[i]);
    }
    double dv;
    double beta0;
    MT_CHECK(myb_edelbaum(7546.0, 7440.0, 0.0174, &dv, &beta0) == MYB_OK);
    hash_double(h, dv);
    hash_double(h, beta0);
  }

  /* SGP4 (ext/sgp4) and SPK type 10 */
  {
    const char* lines[2][2] = {
        {"1 09880U 77021A   06176.56157475  .00000421  00000-0  10000-3 0  9814",
         "2 09880  64.5968 349.3786 7069051 270.0229  16.3320  2.00813614112380"},
        {"1 06251U 62025E   06176.82412014  .00008885  00000-0  12808-3 0  3985",
         "2 06251  58.0579  54.0425 0030035 139.1568 221.1854 15.56387291  6774"}};
    for (int s = 0; s < 2; ++s) {
      myb_tle tle;
      myb_sgp4 sat;
      MT_CHECK(myb_tle_parse(lines[s][0], lines[s][1], 0, &tle, NULL) == MYB_OK);
      MT_CHECK(myb_sgp4_init(&tle, MYB_SGP4_WGS72, MYB_SGP4_IMPROVED, &sat) == MYB_OK);
      for (int k = 0; k < 6; ++k) {
        double r[3];
        double v[3];
        MT_CHECK(myb_sgp4_propagate(&sat, 977.0 * k - 1000.0, r, v, NULL) == MYB_OK);
        for (int i = 0; i < 3; ++i) {
          hash_double(h, r[i]);
          hash_double(h, v[i]);
        }
      }
    }
    myb_spk* spk = NULL;
    MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
    MT_CHECK(myb_spk_load_file(spk, MYB_TEST_DATA_DIR "/spk/t10.bsp", NULL) == MYB_OK);
    for (size_t seg = 0; seg < myb_spk_segment_count(spk); ++seg) {
      myb_spk_segment_info info;
      myb_spk_segment(spk, seg, &info);
      for (int k = 0; k < 5; ++k) {
        const double t = info.start + (info.stop - info.start) * repro_uniform(rng);
        double st[6];
        MT_CHECK(myb_sgp4_spk10_state(spk, seg, sys, eq, myb_time_from_seconds(t), st) ==
                 MYB_OK);
        for (int i = 0; i < 6; ++i) {
          hash_double(h, st[i]);
        }
      }
    }
    myb_spk_destroy(spk);
  }
  myb_timesys_destroy(sys);
  myb_equinox_destroy(eq);
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
  workload_frames(&h, &rng); /* phase 5 */
  /* MYB_REPRO_PREFIX=<phase>: hash only the workloads up to that phase and compare with its
   * golden value, to show that adding a phase changed no earlier bit */
  const char* prefix = getenv("MYB_REPRO_PREFIX");
  const int upto = prefix != NULL ? atoi(prefix) : 99;
  if (upto >= 6) {
    workload_rel(&h, &rng);
  }
  if (upto >= 7) {
    workload_grav(&h, &rng);
  }
  if (upto >= 8) {
    workload_phase8(&h, &rng);
  }
  if (upto >= 9) {
    workload_phase9(&h, &rng);
  }
  uint64_t expected = golden_hash;
  if (upto == 5) {
    expected = UINT64_C(0x38f060f6206c4bf6);
  } else if (upto == 6) {
    expected = UINT64_C(0x9cb768b2ccddf18b);
  } else if (upto == 7) {
    expected = UINT64_C(0xbd6bc51c0c5d4fda);
  } else if (upto == 8) {
    expected = UINT64_C(0xb1b77fd24243268f);
  }
  printf("repro hash: 0x%016" PRIx64 " over %ld values (golden 0x%016" PRIx64 ")\n", h.state,
         h.values, expected);
  MT_CHECK(h.state == expected);
  return mt_finish("test_repro");
}
