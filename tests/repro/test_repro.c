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
#include "myb/ccsds.h"
#include "myb/cmath.h"
#include "myb/cr3bp.h"
#include "myb/dyn.h"
#include "myb/env.h"
#include "myb/estim.h"
#include "myb/events.h"
#include "myb/force.h"
#include "myb/frames.h"
#include "myb/geometry.h"
#include "myb/grav.h"
#include "myb/integ.h"
#include "myb/io.h"
#include "myb/kepler.h"
#include "myb/kin.h"
#include "myb/lambert.h"
#include "myb/model.h"
#include "myb/num.h"
#include "myb/observe.h"
#include "myb/orbit.h"
#include "myb/pool.h"
#include "myb/regular.h"
#include "myb/rel.h"
#include "myb/surface.h"
#include "myb/thrust.h"
#include "myb/time.h"
#include "myb/transfer.h"
#include "myb/uncert.h"
#include "myb_sgp4.h"

static const uint64_t golden_hash = UINT64_C(0xb76153a8dca45827);

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

/* ---- Phase 10 ----------------------------------------------------------------------------- */

typedef struct repro_mover {
  double x0[3];
  double v[3];
} repro_mover;

static myb_status repro_mover_state(void* context, myb_time t, double state[6]) {
  const repro_mover* m = context;
  const double s = myb_time_seconds(t);
  for (int i = 0; i < 3; ++i) {
    state[i] = m->x0[i] + m->v[i] * s;
    state[3 + i] = m->v[i];
  }
  return MYB_OK;
}

static void workload_phase10(hasher* h, mt_rng* rng) {
  /* compact binaries */
  for (int k = 0; k < 20; ++k) {
    const double x[3] = {1e3 + 100.0 * repro_uniform(rng), 50.0 * repro_uniform(rng), 3.0};
    const double v[3] = {0.001 * repro_uniform(rng), 0.03, 0.002 * repro_uniform(rng)};
    double a[3];
    MT_CHECK(myb_pn_binary_acceleration(0.6, 0.4, 1.0, MYB_PN_1 | MYB_PN_2 | MYB_PN_25, x, v, a) ==
             MYB_OK);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, a[i]);
    }
  }
  hash_double(h, myb_peters_period_derivative(1.9e20, 1.85e20, 299792458.0, 27907.0, 0.617));

  /* observables */
  repro_mover target = {{3e11, -1e11, 5e10}, {2.5e4, 1.8e4, -3e3}};
  repro_mover sun = {{1e8, -2e8, 3e7}, {12.0, -9.0, 0.5}};
  repro_mover station = {{1.4e11, 5e10, 2e10}, {-1e4, 2.9e4, 1e3}};
  const myb_light_body body = {.gm = 1.32712440041e20, .state = repro_mover_state, .context = &sun};
  const myb_light_model model = {.c = 299792458.0, .gamma = 1.0, .bodies = &body, .count = 1};
  for (int k = 0; k < 6; ++k) {
    const myb_time t = myb_time_from_seconds(1e5 * k + 1e3 * repro_uniform(rng));
    double obs[6];
    MT_CHECK(repro_mover_state(&station, t, obs) == MYB_OK);
    myb_light_path path;
    MT_CHECK(myb_light_time(&model, repro_mover_state, &target, t, obs,
                            (myb_light_direction)(k % 2), &path) == MYB_OK);
    hash_double(h, path.light_time);
    hash_double(h, path.range);
    hash_double(h, path.delay);
    hash_double(h, path.range_rate);
    double total;
    MT_CHECK(myb_round_trip(&model, repro_mover_state, &station, repro_mover_state, &target, t,
                            &total, NULL, NULL) == MYB_OK);
    hash_double(h, total);
    const double p[3] = {0.6, -0.48, 0.64};
    double ab[3];
    myb_aberration(p, obs + 3, 299792458.0, ab);
    double ld[3];
    const double e[3] = {0.0, 1.0, 0.0};
    myb_light_deflection(1.32712440041e20, 299792458.0, 1.0, p, p, e, 1.5e11, 1e-9, ld);
    double ra;
    double dec;
    myb_direction_radec(ab, &ra, &dec);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, ab[i]);
      hash_double(h, ld[i]);
    }
    hash_double(h, ra);
    hash_double(h, dec);
  }

  /* magnetic field, Lorentz force, co-rotation, plasma stream */
  double g[10] = {0.0, -2.9e-5, -1.5e-6, -2.5e-6, 3.0e-6, 1.7e-6, 1.3e-6, -2.3e-6, 1.2e-6, 6e-7};
  double hh[10] = {0.0, 0.0, 4.6e-6, 0.0, -2.9e-6, -6e-7, 0.0, -1.6e-7, 2.4e-7, -5.4e-7};
  myb_magnetic_field* field = NULL;
  MT_CHECK(myb_magnetic_field_create(6371.2e3, 3, g, hh, NULL, &field) == MYB_OK);
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_work_create(3, NULL, &work) == MYB_OK);
  for (int k = 0; k < 10; ++k) {
    const double r[3] = {7e6 + 1e6 * repro_uniform(rng), -3e6 * repro_uniform(rng),
                         2e6 * repro_uniform(rng)};
    double b[3];
    MT_CHECK(myb_magnetic_field_eval(field, work, r, b) == MYB_OK);
    const double omega[3] = {0.0, 0.0, 7.292115e-5};
    double e[3];
    myb_corotation_field(omega, r, b, e);
    const double v[3] = {1e3, 7e3, -2e3};
    double a[3];
    myb_lorentz_acceleration(1e-3, e, v, b, a);
    double w[3];
    myb_force_plasma_stream(8e-21, 1.496e11, 4e5, 0.02, r, v, w);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, b[i]);
      hash_double(h, a[i]);
      hash_double(h, w[i]);
    }
  }
  myb_gravity_work_destroy(work);
  myb_magnetic_field_destroy(field);
}

static void workload_phase11(hasher* h, mt_rng* rng) {
  /* compact binaries at 3.5PN with spins */
  myb_pn_binary b = {.gm1 = 0.6, .gm2 = 0.4, .c = 1.0, .orders = MYB_PN_ALL, .r0 = 3.0,
                     .kappa1 = 1.0, .kappa2 = 4.5};
  for (int k = 0; k < 20; ++k) {
    const double x[3] = {30.0 + 10.0 * repro_uniform(rng), 5.0 * repro_uniform(rng), 1.0};
    const double v[3] = {0.01 * repro_uniform(rng), 0.17, 0.02 * repro_uniform(rng)};
    const double c1[3] = {0.3, -0.5 * repro_uniform(rng), 0.6};
    const double c2[3] = {0.7 * repro_uniform(rng), 0.2, -0.4};
    double a[3];
    double d1[3];
    double d2[3];
    MT_CHECK(myb_pn_binary_evaluate(&b, x, v, c1, c2, a, d1, d2) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, a[i]);
      hash_double(h, d1[i]);
      hash_double(h, d2[i]);
    }
  }
  /* electron columns and dispersive delay */
  const myb_plasma_profile prof = {.radius = 6.957e8, .count = 3,
                                   .density = {2.99e14, 1.55e14, 3.44e11},
                                   .power = {16.0, 6.0, 2.3}};
  for (int k = 0; k < 6; ++k) {
    const double p1[3] = {1.496e11, 1e9 * repro_uniform(rng), 0.0};
    const double p2[3] = {-2e11, (3.0 + 20.0 * repro_uniform(rng)) * 6.957e8, 1e8};
    double n = 0.0;
    MT_CHECK(myb_plasma_column(&prof, p1, p2, &n) == MYB_OK);
    hash_double(h, n);
    hash_double(h, myb_plasma_group_delay(n, 8.4e9, 1.602176634e-19, 9.1093837139e-31,
                                          8.8541878188e-12));
  }
  /* charging */
  const myb_plasma_species wind[2] = {
      {.charge = -1.602176634e-19, .mass = 9.1093837139e-31, .density = 5e6,
       .temperature = 12.0 * 1.602176634e-19, .drift = 4e5},
      {.charge = 1.602176634e-19, .mass = 1.67262192595e-27, .density = 5e6,
       .temperature = 8.0 * 1.602176634e-19, .drift = 4e5}};
  const myb_photoemission photo = {.current_density = 4e-5,
                                   .temperature = 2.0 * 1.602176634e-19,
                                   .charge = 1.602176634e-19};
  for (int k = 0; k < 6; ++k) {
    double i = 0.0;
    MT_CHECK(myb_charging_current(wind, 2, &photo, 1e-6, -20.0 + 40.0 * repro_uniform(rng), &i) ==
             MYB_OK);
    hash_double(h, i);
  }
  double phi = 0.0;
  MT_CHECK(myb_charging_equilibrium(wind, 2, &photo, 1e-6, &phi) == MYB_OK);
  hash_double(h, phi);
  MT_CHECK(myb_charging_equilibrium(wind, 2, NULL, 1e-6, &phi) == MYB_OK);
  hash_double(h, phi);
}

/* ---- Phase 12: the environment, the ready-made terms, the propagator ---------------------- */

typedef struct p12_mover {
  double radius;
  double rate;
  double phase;
  double inclination;
} p12_mover;

static myb_status p12_mover_state(void* context, myb_time t, double y[6]) {
  const p12_mover* m = context;
  const double u = m->phase + m->rate * myb_time_seconds(t);
  const double cu = myb_cos(u);
  const double su = myb_sin(u);
  const double ci = myb_cos(m->inclination);
  const double si = myb_sin(m->inclination);
  y[0] = m->radius * cu;
  y[1] = m->radius * su * ci;
  y[2] = m->radius * su * si;
  y[3] = -m->radius * m->rate * su;
  y[4] = m->radius * m->rate * cu * ci;
  y[5] = m->radius * m->rate * cu * si;
  return MYB_OK;
}

typedef struct p12_spinner {
  double angles[3];
  double rate;
} p12_spinner;

static myb_status p12_spinner_orientation(void* context, myb_time t, myb_orientation* o) {
  const p12_spinner* s = context;
  const double angles[3] = {s->angles[0], s->angles[1],
                            s->angles[2] + s->rate * myb_time_seconds(t)};
  const double rates[3] = {0.0, 0.0, s->rate};
  myb_rotation_euler313(angles, rates, &o->rotation);
  o->ut = myb_time_add(t, -0.2);
  o->xp = 4e-7;
  o->yp = 1.5e-6;
  return MYB_OK;
}

static double p12_apsis(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  return y[0] * y[3] + y[1] * y[4] + y[2] * y[5];
}

static void workload_phase12(hasher* h, mt_rng* rng) {
  /* scales of the providers from the table */
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/time-scales.tk", NULL) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/time/tdb-tt-fb1990.tk", NULL) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/grav/iers2010-tides.tk", NULL) == MYB_OK);
  myb_timesys* sys = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, NULL) == MYB_OK);
  myb_tides* tides = NULL;
  MT_CHECK(myb_tides_create(pool, "IERS10_TIDE", NULL, &tides, NULL) == MYB_OK);
  myb_kpool_destroy(pool);
  if (sys == NULL || tides == NULL) {
    return;
  }
  const myb_time epoch = {801964800, 0.25};
  p12_mover sun = {.radius = 1.496e11, .rate = 1.99e-7, .phase = 1.0, .inclination = 0.4};
  p12_mover moon = {.radius = 3.84e8, .rate = 2.66e-6, .phase = 0.4, .inclination = 0.3};
  p12_spinner spin = {.angles = {0.2, 0.4, 1.0}, .rate = 7.292115e-5};
  const myb_env_options options = {.epoch = epoch, .timesys = sys, .scale = "TT",
                                   .table_step = 1800.0};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  size_t planet = 0;
  size_t i_sun = 0;
  size_t i_moon = 0;
  const myb_env_body b_planet = {.id = 399, .gm = 3.986004415e14,
                                 .orientation = p12_spinner_orientation,
                                 .orientation_context = &spin, .radius = 6378137.0,
                                 .spin = {0.0, 0.0, 9.8e8}};
  const myb_env_body b_sun = {.id = 10, .gm = 1.32712440018e20, .state = p12_mover_state,
                              .state_context = &sun, .state_scale = "TDB", .radius = 6.957e8,
                              .luminosity = 3.828e26, .limb = {0.6, 0.0}};
  const myb_env_body b_moon = {.id = 301, .gm = 4.9028e12, .state = p12_mover_state,
                               .state_context = &moon, .state_scale = "TCB", .radius = 1737400.0};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_sun, &i_sun) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_moon, &i_moon) == MYB_OK);
  for (int k = 0; k < 40; ++k) {
    const double t = 2.0e6 * (repro_uniform(rng) - 0.5);
    myb_time out = {};
    MT_CHECK(myb_env_time(env, t, k % 2 == 0 ? "TDB" : "TCB", &out) == MYB_OK);
    hash_time(h, out);
  }

  /* a model of every ready-made term */
  double c[15] = {0.0};
  double s[15] = {0.0};
  c[0] = 1.0;
  c[3] = -4.84165e-4;
  c[5] = 2.439e-6;
  s[5] = -1.400e-6;
  c[6] = 9.57e-7;
  c[7] = 2.03e-6;
  s[7] = 2.48e-7;
  c[10] = 5.4e-7;
  s[12] = 6.6e-7;
  myb_gravity_field* field = NULL;
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_field_create(b_planet.gm, 6378136.3, 4, c, s, NULL, &field) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(4, NULL, &work) == MYB_OK);
  myb_object object = {.mass = 400.0, .area = 0.28, .cr = 1.13, .cd = 2.2};
  myb_field_term gravity = {.env = env, .body = planet, .field = field, .work = work,
                            .n_min = 0, .n_max = -1, .update = 1, .tides = tides,
                            .tide_flags = MYB_TIDES_SOLID | MYB_TIDES_FREQUENCY |
                                          MYB_TIDES_POLE_SOLID,
                            .permanent = MYB_PERMANENT_TIDE_KEEP, .raisers = {i_sun, i_moon},
                            .raiser_count = 2};
  myb_third_body_term third_sun = {env, i_sun};
  myb_third_body_term third_moon = {env, i_moon};
  myb_relativity_term relativity = {.env = env, .ppn = {299792458.0, 1.0, 1.0},
                                    .central = planet, .external = i_sun};
  myb_radiation_term pressure = {.env = env, .source = i_sun, .object = &object,
                                 .c = 299792458.0, .occulters = {planet, i_moon},
                                 .occulter_count = 2};
  const myb_planet_radiation knocke = {.radius = 6378137.0,
                                       .albedo = {0.34, 0.0, 0.10, 0.0, 0.29},
                                       .emissivity = {0.68, 0.0, -0.07, 0.0, -0.18},
                                       .omega = 1.99106e-7};
  myb_planet_radiation_term earthlight = {.env = env, .planet = planet, .source = i_sun,
                                          .model = &knocke, .object = &object, .phase0 = 0.7,
                                          .rings = 5};
  myb_exponential_atmosphere air = {.ellipsoid = {6378137.0, 1.0 / 298.257223563},
                                    .count = 2,
                                    .base = {200000.0, 600000.0},
                                    .density = {2.5e-10, 1.1e-13},
                                    .scale = {37000.0, 71000.0}};
  myb_drag_term drag = {.env = env, .body = planet, .density = myb_density_exponential,
                        .density_context = &air, .object = &object};
  myb_empirical_term empirical = {.constant = {1e-10, -2e-10, 3e-11},
                                  .cosine = {4e-11, 5e-10, -6e-10},
                                  .sine = {-7e-11, 8e-11, 9e-10}};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "field", myb_term_field, NULL, &gravity, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "sun", myb_term_third_body, myb_term_third_body_jacobian,
                               &third_sun, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "moon", myb_term_third_body, myb_term_third_body_jacobian,
                               &third_moon, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "schwarzschild", myb_term_schwarzschild, NULL,
                               &relativity, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "lense-thirring", myb_term_lense_thirring, NULL,
                               &relativity, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "de sitter", myb_term_de_sitter, NULL, &relativity,
                               NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "sunlight", myb_term_radiation, NULL, &pressure, NULL) ==
           MYB_OK);
  MT_CHECK(myb_force_model_add(model, "earthlight", myb_term_planet_radiation, NULL,
                               &earthlight, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "drag", myb_term_drag, NULL, &drag, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "empirical", myb_term_empirical, NULL, &empirical,
                               NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(model, "cr", &object.cr, 0.01, NULL) == MYB_OK);
  for (int k = 0; k < 12; ++k) {
    const double t = 86400.0 * (repro_uniform(rng) - 0.5);
    const double x[3] = {6.9e6 + 2e5 * repro_uniform(rng), 1e6 * (repro_uniform(rng) - 0.5),
                         2e6 * (repro_uniform(rng) - 0.5)};
    const double v[3] = {300.0 * (repro_uniform(rng) - 0.5), 6.5e3 + 100.0 * repro_uniform(rng),
                         3.8e3};
    double a[3];
    MT_CHECK(myb_force_model_acceleration(model, t, x, v, a) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, a[i]);
    }
    for (size_t term = 0; term < myb_force_model_count(model); term += 3) {
      MT_CHECK(myb_force_model_term(model, term, t, x, v, a) == MYB_OK);
      for (int i = 0; i < 3; ++i) {
        hash_double(h, a[i]);
      }
    }
    if (k < 2) {
      double jx[9];
      double jv[9];
      double jp[3];
      MT_CHECK(myb_force_model_jacobian(model, t, x, v, jx, jv, jp) == MYB_OK);
      for (int i = 0; i < 9; ++i) {
        hash_double(h, jx[i]);
        hash_double(h, jv[i]);
      }
      for (int i = 0; i < 3; ++i) {
        hash_double(h, jp[i]);
      }
    }
  }

  /* the propagator: states, an event, partials */
  const double y0[6] = {7.0e6, 1.0e5, -2.0e5, -150.0, 6.6e3, 3.7e3};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, NULL, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 100.0, y0, 9100.0) == MYB_OK);
  const double times[4] = {100.0, 1234.5, 5000.0, 7300.25};
  double states[24];
  MT_CHECK(myb_orbit_states(orbit, times, 4, states) == MYB_OK);
  for (int i = 0; i < 24; ++i) {
    hash_double(h, states[i]);
  }
  int found = 0;
  double te = 0.0;
  double ye[6];
  MT_CHECK(myb_orbit_event(orbit, p12_apsis, NULL, 9100.0, 0, &found, &te, ye) == MYB_OK);
  hash_double(h, (double)found);
  hash_double(h, te);
  for (int i = 0; i < 6 && found; ++i) {
    hash_double(h, ye[i]);
  }
  hash_double(h, (double)myb_orbit_steps(orbit));
  myb_orbit_destroy(orbit);
  const myb_orbit_options variational = {.variational = 1};
  MT_CHECK(myb_orbit_create(model, &variational, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 100.0, y0, -1900.0) == MYB_OK);
  double y[6];
  double phi[36];
  double sensitivity[6];
  MT_CHECK(myb_orbit_partials(orbit, -1500.0, y, phi, sensitivity) == MYB_OK);
  for (int i = 0; i < 6; ++i) {
    hash_double(h, y[i]);
    hash_double(h, sensitivity[i]);
  }
  for (int i = 0; i < 36; ++i) {
    hash_double(h, phi[i]);
  }
  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
  myb_env_destroy(env);
  myb_tides_destroy(tides);
  myb_timesys_destroy(sys);
}

/* ---- Phase 14: gradient of the field, analytic Jacobians, estimation ------------------------ */

static void workload_phase14(hasher* h, mt_rng* rng) {
  /* the gradient of a field of degree 8 */
  double c[45] = {0.0};
  double s[45] = {0.0};
  c[0] = 1.0;
  for (int n = 2; n <= 8; ++n) {
    for (int m = 0; m <= n; ++m) {
      const double size = (n == 2 && m == 0 ? 4.8e-4 : 2e-6) / (n - 1);
      c[n * (n + 1) / 2 + m] = size * (2.0 * repro_uniform(rng) - 1.0);
      s[n * (n + 1) / 2 + m] = m == 0 ? 0.0 : size * (2.0 * repro_uniform(rng) - 1.0);
    }
  }
  const double gm = 3.986004415e14;
  myb_gravity_field* field = NULL;
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_field_create(gm, 6378136.3, 8, c, s, NULL, &field) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(8, NULL, &work) == MYB_OK);
  for (int k = 0; k < 12; ++k) {
    const double r[3] = {k < 2 ? 0.0 : 1.4e7 * (repro_uniform(rng) - 0.5),
                         k < 2 ? 0.0 : 1.4e7 * (repro_uniform(rng) - 0.5),
                         (k == 1 ? -1.0 : 1.0) * (6.6e6 + 1.0e6 * repro_uniform(rng))};
    double g[9];
    MT_CHECK(myb_gravity_field_gradient(field, work, r, k % 3 == 2 ? 2 : 0, -1, g) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, g[i]);
    }
  }
  /* the Jacobians of the relativistic accelerations */
  const myb_ppn ppn = {299792458.0, 1.0, 1.0};
  for (int k = 0; k < 6; ++k) {
    const double r[3] = {7.0e6 + 1e6 * repro_uniform(rng), 2e6 * (repro_uniform(rng) - 0.5), 1e6};
    const double v[3] = {300.0 * (repro_uniform(rng) - 0.5), 7.0e3, 2.0e3 * repro_uniform(rng)};
    const double j[3] = {1e8, -2e8 * repro_uniform(rng), 9.8e8};
    const double big_r[3] = {1.4e11, 4e10 * repro_uniform(rng), 1e10};
    const double big_v[3] = {-9e3, 2.8e4, 1e3 * repro_uniform(rng)};
    double dr[9];
    double dv[9];
    MT_CHECK(myb_rel_schwarzschild_jacobian(&ppn, gm, r, v, dr, dv) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, dr[i]);
      hash_double(h, dv[i]);
    }
    MT_CHECK(myb_rel_lense_thirring_jacobian(&ppn, gm, j, r, v, dr, dv) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, dr[i]);
      hash_double(h, dv[i]);
    }
    MT_CHECK(myb_rel_de_sitter_jacobian(&ppn, 1.32712440018e20, big_r, big_v, dv) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, dv[i]);
    }
  }

  /* a planet with a field, a rotation and an atmosphere; a station; the drag coefficient as
   * the parameter */
  p12_spinner spin = {.angles = {0.2, 0.4, 1.0}, .rate = 7.292115e-5};
  p12_mover sun = {.radius = 1.496e11, .rate = 1.99e-7, .phase = 1.0, .inclination = 0.4};
  const myb_env_options options = {.epoch = {801964800, 0.5}};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  size_t planet = 0;
  size_t i_sun = 0;
  const myb_env_body b_planet = {.gm = gm, .orientation = p12_spinner_orientation,
                                 .orientation_context = &spin, .radius = 6378137.0,
                                 .spin = {0.0, 0.0, 9.8e8}};
  const myb_env_body b_sun = {.gm = 1.32712440018e20, .state = p12_mover_state,
                              .state_context = &sun, .radius = 6.957e8, .luminosity = 3.828e26};
  MT_CHECK(myb_env_add_body(env, &b_planet, &planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &b_sun, &i_sun) == MYB_OK);
  myb_object object = {.mass = 200.0, .area = 4.0, .cr = 1.3, .cd = 2.2};
  myb_field_term gravity = {.env = env, .body = planet, .field = field, .work = work,
                            .n_min = 0, .n_max = -1};
  myb_exponential_atmosphere air = {.ellipsoid = {6378137.0, 1.0 / 298.257223563},
                                    .count = 1,
                                    .base = {400000.0},
                                    .density = {3.0e-12},
                                    .scale = {58000.0}};
  myb_drag_term drag = {.env = env, .body = planet, .density = myb_density_exponential,
                        .density_context = &air, .object = &object};
  myb_radiation_term light = {.env = env, .source = i_sun, .object = &object,
                              .c = 299792458.0, .occulters = {planet}, .occulter_count = 1};
  myb_relativity_term relativity = {.env = env, .ppn = ppn, .central = planet,
                                    .external = i_sun};
  myb_force_model* model = NULL;
  size_t i_drag = 0;
  size_t i_cd = 0;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "field", myb_term_field, myb_term_field_jacobian,
                               &gravity, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "drag", myb_term_drag, myb_term_drag_jacobian, &drag,
                               &i_drag) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "light", myb_term_radiation, myb_term_radiation_jacobian,
                               &light, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "schwarzschild", myb_term_schwarzschild,
                               myb_term_schwarzschild_jacobian, &relativity, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "lense-thirring", myb_term_lense_thirring,
                               myb_term_lense_thirring_jacobian, &relativity, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "de sitter", myb_term_de_sitter,
                               myb_term_de_sitter_jacobian, &relativity, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(model, "cd", &object.cd, 0.01, &i_cd) == MYB_OK);
  MT_CHECK(myb_force_model_bind_parameter(model, i_cd, i_drag) == MYB_OK);
  const double y0[6] = {6.8e6, 1.0e5, -2.0e5, -150.0, 6.5e3, 4.1e3};
  for (int k = 0; k < 3; ++k) {
    const double x[3] = {y0[0] + 1e5 * repro_uniform(rng), y0[1], y0[2] - 3e5 * repro_uniform(rng)};
    double jx[9];
    double jv[9];
    double jp[3];
    MT_CHECK(myb_force_model_jacobian(model, 100.0 * k, x, y0 + 3, jx, jv, jp) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, jx[i]);
      hash_double(h, jv[i]);
    }
    for (int i = 0; i < 3; ++i) {
      hash_double(h, jp[i]);
    }
  }

  /* observations of the truth with noise; least squares and the filter */
  myb_station station;
  MT_CHECK(myb_station_init(&station, env, planet, &air.ellipsoid, 0.3, 0.9, 150.0) == MYB_OK);
  enum { OBSERVATIONS = 48 };
  static myb_observation obs[OBSERVATIONS];
  const myb_measure_function measures[4] = {myb_measure_range, myb_measure_range_rate,
                                            myb_measure_radec, myb_measure_azel};
  const size_t counts[4] = {1, 1, 2, 2};
  const double sigmas[4] = {1.0, 1e-3, 2e-6, 2e-6};
  myb_rng noise;
  myb_rng_seed(&noise, 14);
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, NULL, NULL, &orbit) == MYB_OK);
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, 5000.0) == MYB_OK);
  for (int k = 0; k < OBSERVATIONS; ++k) {
    myb_observation* o = &obs[k];
    const int kind = k % 4;
    memset(o, 0, sizeof *o);
    o->t = 100.0 * (k + 1);
    o->measure = measures[kind];
    o->context = &station;
    o->count = counts[kind];
    o->periodic = kind >= 2 ? 1u : 0u;
    double y[6];
    MT_CHECK(myb_orbit_state(orbit, o->t, y) == MYB_OK);
    MT_CHECK(o->measure(o->context, o->t, y, o->value, NULL) == MYB_OK);
    for (size_t i = 0; i < o->count; ++i) {
      o->sigma[i] = sigmas[kind];
      o->value[i] += sigmas[kind] * myb_rng_normal(&noise);
      hash_double(h, o->value[i]);
    }
  }
  myb_orbit_destroy(orbit);
  double y[6] = {y0[0] + 120.0, y0[1] - 80.0, y0[2] + 60.0, y0[3] + 0.1, y0[4] - 0.05, y0[5]};
  object.cd = 1.8;
  myb_fit_result result;
  double covariance[49];
  static double residuals[2 * OBSERVATIONS];
  MT_CHECK(myb_fit_orbit(model, 0.0, y, obs, OBSERVATIONS, NULL, NULL, &result, covariance,
                         residuals) == MYB_OK);
  for (int i = 0; i < 6; ++i) {
    hash_double(h, y[i]);
  }
  hash_double(h, object.cd);
  hash_double(h, result.rms);
  hash_double(h, (double)result.iterations);
  for (int i = 0; i < 49; ++i) {
    hash_double(h, covariance[i]);
  }
  for (size_t i = 0; i < result.values; ++i) {
    hash_double(h, residuals[i]);
  }
  double p0[49] = {0.0};
  for (int i = 0; i < 3; ++i) {
    p0[8 * i] = 1.0e4;
    p0[8 * (3 + i)] = 1.0e-2;
  }
  p0[48] = 0.04;
  double yf[6] = {y0[0] + 50.0, y0[1] - 30.0, y0[2] + 20.0, y0[3] + 0.02, y0[4], y0[5] - 0.03};
  object.cd = 2.0;
  const myb_filter_options filter_options = {.process_noise = 1e-10};
  myb_filter* filter = NULL;
  MT_CHECK(myb_filter_create(model, &filter_options, NULL, &filter) == MYB_OK);
  MT_CHECK(myb_filter_start(filter, 0.0, yf, p0, 5000.0) == MYB_OK);
  for (int k = 0; k < OBSERVATIONS; ++k) {
    double innovation[MYB_MEASURE_MAX];
    MT_CHECK(myb_filter_update(filter, &obs[k], innovation) == MYB_OK);
    hash_double(h, innovation[0]);
  }
  MT_CHECK(myb_filter_predict(filter, 5000.0) == MYB_OK);
  myb_filter_get(filter, NULL, yf, covariance);
  for (int i = 0; i < 6; ++i) {
    hash_double(h, yf[i]);
  }
  hash_double(h, object.cd);
  for (int i = 0; i < 49; ++i) {
    hash_double(h, covariance[i]);
  }
  myb_filter_destroy(filter);
  myb_force_model_destroy(model);
  myb_env_destroy(env);
  myb_gravity_work_destroy(work);
  myb_gravity_field_destroy(field);
}

/* ---- Phase 15: output ----------------------------------------------------------------------- */

static void hash_bytes(hasher* h, const void* data, size_t size) {
  const unsigned char* p = data;
  for (size_t k = 0; k < size; ++k) {
    h->state ^= p[k];
    h->state *= UINT64_C(0x100000001B3);
  }
  h->values++;
}

static void workload_phase15(hasher* h, mt_rng* rng) {
  static char text[32768];
  size_t length = 0;
  /* numbers as text */
  for (int k = 0; k < 300; ++k) {
    const double x = mt_double(rng, -300, 300);
    hash_bytes(h, text, myb_format_double(x, text));
    MT_CHECK(myb_format_scientific(x, 1 + k % 20, text, sizeof text, &length) == MYB_OK);
    hash_bytes(h, text, length);
    const double y = (repro_uniform(rng) - 0.5) * 1.0e6;
    MT_CHECK(myb_format_fixed(y, k % 25, text, sizeof text, &length) == MYB_OK);
    hash_bytes(h, text, length);
  }
  /* an ephemeris file: tabulated states of both kinds and polynomials fitted to a motion */
  p12_mover mover = {.radius = 7.2e6, .rate = 1.03e-3, .phase = 0.4, .inclination = 0.9};
  enum { STATES = 40 };
  static double epochs[STATES];
  static double states[6 * STATES];
  static myb_calendar labels[STATES];
  for (int k = 0; k < STATES; ++k) {
    epochs[k] = 100.0 + 60.0 * k + 3.0 * repro_uniform(rng);
    const myb_time t = myb_time_from_seconds(epochs[k]);
    MT_CHECK(p12_mover_state(&mover, t, states + 6 * k) == MYB_OK);
    myb_time_to_calendar(t, &labels[k]);
  }
  myb_spk_writer* writer = NULL;
  MT_CHECK(myb_spk_writer_create("repro", NULL, &writer) == MYB_OK);
  const myb_spk_segment_header tabulated = {.target = -15, .center = 399, .frame = 1,
                                            .name = "states"};
  MT_CHECK(myb_spk_writer_add_states(writer, &tabulated, 9, 7, epochs, states, STATES) == MYB_OK);
  MT_CHECK(myb_spk_writer_add_states(writer, &tabulated, 13, 5, epochs, states, STATES) ==
           MYB_OK);
  const myb_spk_segment_header fitted = {.target = -16, .center = 399, .frame = 1, .name = "fit"};
  double fit_error = 0.0;
  MT_CHECK(myb_spk_writer_add_chebyshev(writer, &fitted, 3, p12_mover_state, &mover, 0.0, 6000.0,
                                        12, 11, &fit_error) == MYB_OK);
  hash_double(h, fit_error);
  MT_CHECK(myb_spk_writer_add_chebyshev(writer, &fitted, 2, p12_mover_state, &mover, 6000.0,
                                        9000.0, 5, 9, &fit_error) == MYB_OK);
  hash_double(h, fit_error);
  static unsigned char image[65536];
  const size_t size = myb_spk_writer_size(writer);
  MT_CHECK(size <= sizeof image);
  if (size <= sizeof image) {
    MT_CHECK(myb_spk_writer_write(writer, image, sizeof image) == MYB_OK);
    hash_bytes(h, image, size);
  }
  myb_spk_writer_destroy(writer);

  /* the messages of the CCSDS: written, read, the numbers read hashed */
  const myb_oem_segment_data segment = {.object_name = "REPRO",
                                        .object_id = "2026-001A",
                                        .center_name = "EARTH",
                                        .ref_frame = "GCRF",
                                        .time_system = "TT",
                                        .interpolation = "HERMITE",
                                        .interpolation_degree = 5,
                                        .count = STATES,
                                        .epochs = labels,
                                        .states = states};
  const myb_oem_data ephemeris = {.originator = "MYB",
                                  .creation_date = {2026, 11, 10, 12, 0, 0.5},
                                  .segment_count = 1,
                                  .segments = &segment};
  MT_CHECK(myb_oem_write(&ephemeris, text, sizeof text, &length) == MYB_OK);
  hash_bytes(h, text, length);
  myb_oem* oem = NULL;
  MT_CHECK(myb_oem_load(text, length, NULL, &oem, NULL) == MYB_OK);
  if (oem != NULL) {
    for (size_t k = 0; k < myb_oem_state_count(oem, 0); ++k) {
      myb_calendar epoch;
      double y[6];
      (void)myb_oem_state(oem, 0, k, &epoch, y, NULL);
      hash_double(h, epoch.second);
      for (int i = 0; i < 6; ++i) {
        hash_double(h, y[i]);
      }
    }
    myb_oem_destroy(oem);
  }
  static myb_opm opm;
  memset(&opm, 0, sizeof opm);
  memcpy(opm.originator, "MYB", 4);
  memcpy(opm.object_name, "REPRO", 6);
  memcpy(opm.object_id, "2026-001A", 10);
  memcpy(opm.center_name, "EARTH", 6);
  memcpy(opm.ref_frame, "GCRF", 5);
  memcpy(opm.time_system, "TT", 3);
  opm.creation_date = ephemeris.creation_date;
  opm.epoch = labels[3];
  memcpy(opm.state, states + 18, sizeof opm.state);
  opm.has_elements = 1;
  opm.semi_major_axis = 7.2e6 * (1.0 + 0.01 * repro_uniform(rng));
  opm.eccentricity = 0.01 * repro_uniform(rng);
  opm.inclination = 0.9 + 0.1 * repro_uniform(rng);
  opm.ra_of_asc_node = 6.0 * repro_uniform(rng);
  opm.arg_of_pericenter = 6.0 * repro_uniform(rng);
  opm.anomaly = 6.0 * repro_uniform(rng);
  opm.gm = 3.986004418e14;
  opm.has_spacecraft = 1;
  opm.mass = 500.0 + repro_uniform(rng);
  opm.solar_rad_area = 4.0;
  opm.solar_rad_coeff = 1.3;
  opm.drag_area = 3.5;
  opm.drag_coeff = 2.2;
  opm.maneuver_count = 2;
  for (int k = 0; k < 2; ++k) {
    opm.maneuvers[k].epoch_ignition = labels[10 + 5 * k];
    opm.maneuvers[k].duration = 30.0 * repro_uniform(rng);
    opm.maneuvers[k].delta_mass = -repro_uniform(rng);
    memcpy(opm.maneuvers[k].ref_frame, "RTN", 4);
    for (int i = 0; i < 3; ++i) {
      opm.maneuvers[k].dv[i] = 2.0 * (repro_uniform(rng) - 0.5);
    }
  }
  MT_CHECK(myb_opm_write(&opm, text, sizeof text, &length) == MYB_OK);
  hash_bytes(h, text, length);
  static myb_opm read_back;
  MT_CHECK(myb_opm_parse(text, length, &read_back, NULL) == MYB_OK);
  for (int i = 0; i < 6; ++i) {
    hash_double(h, read_back.state[i]);
  }
  hash_double(h, read_back.inclination);
  hash_double(h, read_back.anomaly);
  hash_double(h, read_back.gm);
  hash_double(h, read_back.maneuvers[1].dv[2]);

  /* an element set: written, as a message, and back */
  myb_tle tle;
  memset(&tle, 0, sizeof tle);
  tle.catalog_number = 123456;
  tle.classification = 'U';
  memcpy(tle.designator, "26001A", 7);
  tle.epoch_year = 2026;
  tle.epoch_day = 1.0 + 364.0 * repro_uniform(rng);
  tle.ndot = 1.0e-5 * (repro_uniform(rng) - 0.5);
  tle.nddot = 0.0;
  tle.bstar = 1.0e-4 * repro_uniform(rng);
  tle.element_number = 42;
  tle.inclination = 98.0 * repro_uniform(rng);
  tle.node = 359.0 * repro_uniform(rng);
  tle.eccentricity = 0.01 * repro_uniform(rng);
  tle.argument_of_perigee = 359.0 * repro_uniform(rng);
  tle.mean_anomaly = 359.0 * repro_uniform(rng);
  tle.mean_motion = 14.0 + repro_uniform(rng);
  tle.revolution = 4321;
  char line1[70];
  char line2[70];
  MT_CHECK(myb_tle_format(&tle, line1, line2) == MYB_OK);
  hash_bytes(h, line1, 69);
  hash_bytes(h, line2, 69);
  myb_tle parsed;
  MT_CHECK(myb_tle_parse(line1, line2, 0, &parsed, NULL) == MYB_OK);
  myb_omm omm;
  MT_CHECK(myb_omm_from_tle(&parsed, &omm) == MYB_OK);
  MT_CHECK(myb_omm_write(&omm, text, sizeof text, &length) == MYB_OK);
  hash_bytes(h, text, length);
  MT_CHECK(myb_omm_parse(text, length, &omm, NULL) == MYB_OK);
  MT_CHECK(myb_tle_from_omm(&omm, &parsed) == MYB_OK);
  hash_double(h, parsed.epoch_day);
  hash_double(h, parsed.mean_motion);
  hash_double(h, parsed.bstar);

  /* an SP3 file and a table of numbers */
  static const char* const ids[2] = {"L51", "L52"};
  static double positions[6 * 8];
  static double velocities[6 * 8];
  for (int e = 0; e < 8; ++e) {
    for (int s = 0; s < 2; ++s) {
      const double* y = states + 6 * (2 * e + s);
      memcpy(positions + 3 * (2 * e + s), y, 3 * sizeof(double));
      memcpy(velocities + 3 * (2 * e + s), y + 3, 3 * sizeof(double));
    }
  }
  const myb_sp3_data orbit = {.version = 'd',
                              .file_type = 'L',
                              .time_system = "TAI",
                              .frame = "GCRF",
                              .orbit_type = "FIT",
                              .agency = "MYB",
                              .data_used = "ORBIT",
                              .satellite_count = 2,
                              .satellites = ids,
                              .epoch_count = 8,
                              .epochs = labels,
                              .interval = 60.0,
                              .positions = positions,
                              .velocities = velocities};
  MT_CHECK(myb_sp3_write(&orbit, text, sizeof text, &length) == MYB_OK);
  hash_bytes(h, text, length);
  MT_CHECK(myb_table_write(NULL, 6, states, STATES, text, sizeof text, &length) == MYB_OK);
  hash_bytes(h, text, length);
  myb_table* table = NULL;
  MT_CHECK(myb_table_load(text, length, NULL, &table, NULL) == MYB_OK);
  if (table != NULL) {
    for (size_t r = 0; r < myb_table_rows(table); ++r) {
      for (int i = 0; i < 6; ++i) {
        hash_double(h, myb_table_row(table, r)[i]);
      }
    }
    myb_table_destroy(table);
  }
}

/* ---- Phase 16: events, breaks, the geometry of surfaces -------------------------------------- */

static void workload_phase16(hasher* h, mt_rng* rng) {
  /* rays and nearest points on an ellipsoid, the angles of illumination */
  const double radii[3] = {6378137.0, 6300000.0, 6356752.3};
  for (int k = 0; k < 40; ++k) {
    double origin[3];
    double look[3];
    const double reach = k % 5 == 4 ? 0.9 : 1.2 + 30.0 * repro_uniform(rng);
    for (int i = 0; i < 3; ++i) {
      origin[i] = reach * radii[i] * (2.0 * repro_uniform(rng) - 1.0);
    }
    for (int i = 0; i < 3; ++i) {
      look[i] = 0.7 * radii[i] * (2.0 * repro_uniform(rng) - 1.0) - origin[i];
    }
    double point[3] = {0.0, 0.0, 0.0};
    double distance = 0.0;
    const int hit = myb_ray_ellipsoid(radii, origin, look, &distance, point);
    hash_double(h, (double)hit);
    hash_double(h, distance);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, point[i]);
    }
    double near[3];
    double height = 0.0;
    MT_CHECK(myb_ellipsoid_nearest(radii, origin, near, &height) == MYB_OK);
    double normal[3];
    myb_ellipsoid_normal(radii, near, normal);
    const double sun[3] = {1.5e11 * repro_uniform(rng), 1.0e11, -3.0e10};
    double angles[3];
    MT_CHECK(myb_illumination_angles(near, normal, sun, origin, &angles[0], &angles[1],
                                     &angles[2]) == MYB_OK);
    hash_double(h, height);
    for (int i = 0; i < 3; ++i) {
      hash_double(h, near[i]);
      hash_double(h, normal[i]);
      hash_double(h, angles[i]);
    }
    const double to_sun[3] = {sun[0] - origin[0], sun[1] - origin[1], sun[2] - origin[2]};
    const double to_body[3] = {-origin[0], -origin[1], -origin[2]};
    hash_double(h, myb_shadow_margin(to_sun, 6.957e8, to_body, 6.3e6, k % 2));
  }
  /* a mesh: rays and solid angles */
  const double size = 1.2e7;
  const double vertices[18] = {size, 0, 0, -size, 0, 0, 0, size, 0, 0, -size, 0, 0, 0, size,
                               0,    0, -size};
  const size_t faces[24] = {0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4,
                            2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5};
  myb_polyhedron* shape = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, 6, faces, 8, 1.0, NULL, &shape) == MYB_OK);
  for (int k = 0; k < 30; ++k) {
    double origin[3];
    double look[3];
    for (int i = 0; i < 3; ++i) {
      origin[i] = (k % 3 == 2 ? 0.2 : 3.0) * size * (2.0 * repro_uniform(rng) - 1.0);
      look[i] = 0.5 * size * (2.0 * repro_uniform(rng) - 1.0) - origin[i];
    }
    double point[3] = {0.0, 0.0, 0.0};
    double distance = 0.0;
    size_t face = 0;
    const int hit = myb_polyhedron_ray(shape, origin, look, &distance, point, &face);
    hash_double(h, (double)hit);
    hash_double(h, distance);
    hash_double(h, (double)face);
    hash_double(h, point[0]);
    hash_double(h, myb_polyhedron_solid_angle(shape, origin));
  }

  /* an orbit through the shadow of the planet: breaks, partials across them, events */
  p12_mover sun_motion = {.radius = 1.495978707e11, .rate = 1.99e-7, .phase = 0.3,
                          .inclination = 0.4};
  p12_mover moon_motion = {.radius = 3.844e8, .rate = 2.66e-6, .phase = 2.0, .inclination = 0.3};
  p12_spinner spin = {.angles = {0.1, 0.3, 0.5}, .rate = 7.29e-5};
  const myb_env_options env_options = {0};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&env_options, NULL, &env) == MYB_OK);
  const myb_env_body planet = {.id = 399, .gm = 3.986004418e14,
                               .orientation = p12_spinner_orientation,
                               .orientation_context = &spin, .radius = 6378136.6};
  const myb_env_body sun = {.id = 10, .gm = 1.32712440018e20, .state = p12_mover_state,
                            .state_context = &sun_motion, .radius = 6.957e8,
                            .luminosity = 3.828e26};
  const myb_env_body moon = {.id = 301, .gm = 4.9028e12, .state = p12_mover_state,
                             .state_context = &moon_motion, .radius = 1.7374e6};
  size_t i_planet = 0;
  size_t i_sun = 0;
  size_t i_moon = 0;
  MT_CHECK(myb_env_add_body(env, &planet, &i_planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &moon, &i_moon) == MYB_OK);
  myb_object object = {.mass = 100.0, .area = 20.0, .cr = 1.3, .cd = 2.2};
  myb_point_mass_term central = {.gm = 3.986004418e14};
  myb_radiation_term light = {.env = env, .source = i_sun, .object = &object, .c = 299792458.0,
                              .occulters = {i_planet, i_moon}, .occulter_count = 2};
  myb_eclipse_event edges[4];
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "central", myb_term_point_mass,
                               myb_term_point_mass_jacobian, &central, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "light", myb_term_radiation, myb_term_radiation_jacobian,
                               &light, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(model, "cr", &object.cr, 1e-3, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_bind_parameter(model, 0, 1) == MYB_OK);
  MT_CHECK(myb_term_radiation_breaks(&light, model, edges) == MYB_OK);
  const double a = 1.2e7;
  const myb_elements el = {.p = a * (1.0 - 0.3 * 0.3), .e = 0.3, .i = 0.2 + 0.1 * repro_uniform(rng),
                           .raan = 0.3, .argp = 1.0, .nu = 0.5};
  myb_pv pv;
  MT_CHECK(myb_pv_from_elements(3.986004418e14, &el, &pv) == MYB_OK);
  double y0[6];
  memcpy(y0, pv.r, 3 * sizeof(double));
  memcpy(y0 + 3, pv.v, 3 * sizeof(double));
  const double span = 40000.0;
  for (int variational = 0; variational < 2; ++variational) {
    const myb_orbit_options options = {.epsilon = 1e-10, .variational = variational};
    myb_orbit* orbit = NULL;
    MT_CHECK(myb_orbit_create(model, &options, NULL, &orbit) == MYB_OK);
    MT_CHECK(myb_orbit_start(orbit, 0.0, y0, span) == MYB_OK);
    double y[6];
    double phi[36];
    double sens[6];
    for (int k = 1; k <= 8; ++k) {
      const double t = span * k / 8.0;
      if (variational) {
        MT_CHECK(myb_orbit_partials(orbit, t, y, phi, sens) == MYB_OK);
        for (int i = 0; i < 36; ++i) {
          hash_double(h, phi[i]);
        }
        for (int i = 0; i < 6; ++i) {
          hash_double(h, sens[i]);
        }
      } else {
        MT_CHECK(myb_orbit_state(orbit, t, y) == MYB_OK);
      }
      for (int i = 0; i < 6; ++i) {
        hash_double(h, y[i]);
      }
    }
    hash_double(h, (double)myb_orbit_steps(orbit));
    hash_double(h, (double)myb_orbit_breaks(orbit));
    myb_orbit_destroy(orbit);
  }
  /* the ready-made events along the same orbit */
  myb_range_event to_moon = {.env = env, .body = i_moon};
  myb_range_event sphere = {.body = MYB_EVENT_ORIGIN, .range = 1.1e7};
  myb_plane_event plane = {.env = env, .body = i_planet, .normal = {0.1, -0.2, 1.0}};
  myb_altitude_event level = {.env = env, .body = i_planet,
                              .ellipsoid = {.a = 6378137.0, .f = 1.0 / 298.257223563},
                              .altitude = 3.0e6};
  myb_station station;
  MT_CHECK(myb_station_init(&station, env, i_planet, &level.ellipsoid, 0.5, 0.2, 100.0) ==
           MYB_OK);
  myb_elevation_event horizon = {.station = &station, .elevation = 0.1};
  myb_shape_event surface = {.env = env, .body = i_planet, .shape = shape};
  const struct {
    myb_event_function g;
    void* context;
    int direction;
  } events[8] = {{myb_event_range_rate, &to_moon, +1}, {myb_event_range_rate, &sphere, -1},
                 {myb_event_range, &sphere, +1},       {myb_event_plane, &plane, +1},
                 {myb_event_altitude, &level, -1},     {myb_event_elevation, &horizon, 0},
                 {myb_event_eclipse, &edges[0], -1},   {myb_event_shape, &surface, -1}};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, NULL, NULL, &orbit) == MYB_OK);
  for (int k = 0; k < 8; ++k) {
    MT_CHECK(myb_orbit_start(orbit, 0.0, y0, span) == MYB_OK);
    /* bound: three events of a kind */
    for (int count = 0; count < 3; ++count) {
      int found = 0;
      double t_event = 0.0;
      double y_event[6] = {0.0};
      MT_CHECK(myb_orbit_event(orbit, events[k].g, events[k].context, span, events[k].direction,
                               &found, &t_event, y_event) == MYB_OK);
      hash_event(h, found, t_event);
      if (!found) {
        break;
      }
      hash_double(h, y_event[0]);
      hash_double(h, y_event[4]);
    }
  }
  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);
  myb_env_destroy(env);
  myb_polyhedron_destroy(shape);
}

/* ---- Phase 17: the small forces, the surface of the object, impulses ------------------------- */

static myb_status p17_density(void* context, myb_time t, const double r[3], double* density) {
  (void)context;
  (void)t;
  *density = 4.0e-12 * (1.0 + 1.0e-8 * r[2]);
  return MYB_OK;
}

static myb_status p17_wind(void* context, myb_time t, const double r[3], double w[3]) {
  (void)context;
  (void)t;
  w[0] = -1.0e-5 * r[1];
  w[1] = 1.0e-5 * r[0];
  w[2] = 3.0;
  return MYB_OK;
}

static void hash3(hasher* h, const double v[3]) {
  for (int i = 0; i < 3; ++i) {
    hash_double(h, v[i]);
  }
}

static void workload_phase17(hasher* h, mt_rng* rng) {
  const double c = 299792458.0;
  const double au = 1.495978707e11;
  /* the functions */
  myb_nongrav comet;
  myb_nongrav asteroid;
  myb_nongrav_marsden(1.0e-8, 2.0e-9, -3.0e-10, &comet);
  myb_nongrav_power(0.0, 1.0e-13, 0.0, au, 2.0, &asteroid);
  for (int k = 0; k < 12; ++k) {
    const double r[3] = {au * (0.3 + 4.0 * repro_uniform(rng)), au * (repro_uniform(rng) - 0.5),
                         0.1 * au * repro_uniform(rng)};
    const double v[3] = {3.0e4 * (repro_uniform(rng) - 0.5), 3.0e4 * repro_uniform(rng), 1.0e3};
    const double rn = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    const double u[3] = {r[0] / rn, r[1] / rn, r[2] / rn};
    double a[3] = {0.0, 0.0, 0.0};
    myb_force_radiation_moving(4.56e-6 * (au / rn) * (au / rn), 1.3, 0.02, u, v, c, a);
    hash3(h, a);
    MT_CHECK(myb_force_nongrav(k % 2 ? &comet : &asteroid, r, v, a) == MYB_OK);
    hash3(h, a);
    hash_double(h, myb_nongrav_g(&comet, rn));
    const double omega[3] = {1.0e-6, -2.0e-6, 7.0e-5 * repro_uniform(rng)};
    const double moon_r[3] = {r[0] * 2.5e-3, r[1] * 2.5e-3, r[2] * 2.5e-3};
    const double moon_v[3] = {v[0] * 0.03, v[1] * 0.03, v[2] * 0.03};
    MT_CHECK(myb_force_tidal_lag(3.986004418e14, 4.9028e12, 0.3, 6378136.6, 600.0, omega, moon_r,
                                 moon_v, k % 2, a) == MYB_OK);
    hash3(h, a);
    myb_force_emission(60.0 + 10.0 * repro_uniform(rng), 0.9, 240.0, c, u, a);
    hash3(h, a);
    const myb_yarkovsky body = {.radius = 100.0 + 900.0 * repro_uniform(rng), .density = 2500.0,
                                .absorptivity = 0.9, .emissivity = 0.9,
                                .thermal_inertia = 50.0 + 400.0 * repro_uniform(rng),
                                .rotation_rate = 3.0e-4, .obliquity = 3.0 * repro_uniform(rng)};
    double diurnal = 0.0;
    double seasonal = 0.0;
    MT_CHECK(myb_yarkovsky_drift(&body, 1361.0 * (au / rn) * (au / rn), 2.0e-7, c, &diurnal,
                                 &seasonal) == MYB_OK);
    hash_double(h, diurnal);
    hash_double(h, seasonal);
  }
  /* a gas on plates, light on a mesh */
  const myb_gas gas = {.temperature = 900.0, .molar_mass = 0.016, .wall_temperature = 300.0,
                       .accommodation = 0.85};
  myb_plate plates[6];
  memset(plates, 0, sizeof plates);
  for (int k = 0; k < 6; ++k) {
    plates[k].area = 1.0 + 0.5 * k;
    plates[k].normal[k / 2] = k % 2 ? 1.0 : -1.0;
    plates[k].specular = 0.1 * k;
    plates[k].diffuse = 0.05 * k;
    plates[k].thermal = k % 2;
  }
  for (int k = 0; k < 8; ++k) {
    const double speed = 200.0 + 8000.0 * repro_uniform(rng);
    const double v[3] = {speed * 0.6, -speed * 0.48, speed * 0.64};
    double a[3] = {0.0, 0.0, 0.0};
    double cd = 0.0;
    MT_CHECK(myb_force_gas_plates(3.0e-12, &gas, 100.0, plates, 6, v, a) == MYB_OK);
    MT_CHECK(myb_gas_sphere_cd(&gas, speed, &cd) == MYB_OK);
    hash3(h, a);
    hash_double(h, cd);
  }
  /* two boxes: a cube before a slab */
  double vertices[48];
  size_t faces[72];
  static const size_t quads[6][4] = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4},
                                     {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
  const double lo[2][3] = {{-0.5, -0.5, -0.5}, {-3.2, -1.5, -1.5}};
  const double hi[2][3] = {{0.5, 0.5, 0.5}, {-3.0, 1.5, 1.5}};
  for (size_t box = 0; box < 2; ++box) {
    for (size_t k = 0; k < 8; ++k) {
      for (size_t i = 0; i < 3; ++i) {
        vertices[3 * (8 * box + k) + i] = (k >> i) & 1u ? hi[box][i] : lo[box][i];
      }
    }
    for (size_t q = 0; q < 6; ++q) {
      size_t* f = faces + 6 * (6 * box + q);
      f[0] = 8 * box + quads[q][0];
      f[1] = 8 * box + quads[q][1];
      f[2] = 8 * box + quads[q][2];
      f[3] = 8 * box + quads[q][0];
      f[4] = 8 * box + quads[q][2];
      f[5] = 8 * box + quads[q][3];
    }
  }
  myb_polyhedron* craft = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, 16, faces, 24, 1.0, NULL, &craft) == MYB_OK);
  const myb_plate skin = {.specular = 0.2, .diffuse = 0.3, .thermal = 1};
  for (int k = 0; k < 6; ++k) {
    const double raw[3] = {-1.0, 0.4 * (repro_uniform(rng) - 0.5), 0.4 * (repro_uniform(rng) - 0.5)};
    const double n = sqrt(raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]);
    const double u[3] = {raw[0] / n, raw[1] / n, raw[2] / n};
    double force[3];
    double torque[3];
    MT_CHECK(myb_force_radiation_mesh(craft, &skin, 4.5e-6, u, k % 2 ? 4 : 0, force, torque) ==
             MYB_OK);
    hash3(h, force);
    hash3(h, torque);
  }

  /* the terms and the attitude laws in an environment */
  p12_mover sun_motion = {.radius = au, .rate = 1.99e-7, .phase = 0.3, .inclination = 0.4};
  p12_spinner spin = {.angles = {0.1, 0.3, 0.5}, .rate = 7.29e-5};
  const myb_env_options env_options = {0};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&env_options, NULL, &env) == MYB_OK);
  const myb_env_body planet = {.id = 399, .gm = 3.986004418e14,
                               .orientation = p12_spinner_orientation,
                               .orientation_context = &spin, .radius = 6378136.6};
  const myb_env_body sun = {.id = 10, .gm = 1.32712440018e20, .state = p12_mover_state,
                            .state_context = &sun_motion, .radius = 6.957e8,
                            .luminosity = 3.828e26};
  size_t i_planet = 0;
  size_t i_sun = 0;
  MT_CHECK(myb_env_add_body(env, &planet, &i_planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  const myb_thrust_arc burn = {.start = 500.0, .stop = 900.0, .thrust = 20.0, .exhaust = 3000.0,
                               .steering = MYB_STEER_RTN, .direction = {0.0, 1.0, 0.0}};
  myb_thrust thrust;
  MT_CHECK(myb_thrust_init(&thrust, 400.0, &burn, 1) == MYB_OK);
  myb_object object = {.mass = 400.0, .area = 6.0, .cr = 1.3, .cd = 2.2,
                       .mass_function = myb_thrust_mass_at, .mass_context = &thrust};
  myb_orbital_attitude orbital = {.env = env, .body = i_planet};
  myb_sun_attitude yaw = {.env = env, .body = i_planet, .source = i_sun};
  myb_spin_attitude turning = {.m0 = {0, 1, 0, 0, 0, 1, 1, 0, 0}, .rate = 0.01, .phase = 0.2};
  myb_radiation_term moving = {.env = env, .source = i_sun, .object = &object, .c = c,
                               .occulters = {i_planet}, .occulter_count = 1, .moving = 1};
  myb_surface_radiation_term panels = {.env = env, .source = i_sun, .object = &object, .c = c,
                                       .occulters = {i_planet}, .occulter_count = 1,
                                       .attitude = myb_attitude_sun_nadir,
                                       .attitude_context = &yaw, .plates = plates,
                                       .plate_count = 6};
  myb_surface_radiation_term hull = panels;
  hull.plates = NULL;
  hull.plate_count = 0;
  hull.shape = craft;
  hull.surface = &skin;
  hull.shadows = 3;
  hull.attitude = myb_attitude_orbital;
  hull.attitude_context = &orbital;
  myb_surface_drag_term flow = {.env = env, .body = i_planet, .density = p17_density,
                                .wind = p17_wind, .object = &object, .gas = gas,
                                .attitude = myb_attitude_spin, .attitude_context = &turning,
                                .plates = plates, .plate_count = 6};
  myb_drag_term drag = {.env = env, .body = i_planet, .density = p17_density, .object = &object,
                        .wind = p17_wind};
  myb_recoil_term recoil = {.object = &object, .c = c, .power = 300.0, .efficiency = 2.0 / 3.0,
                            .direction = {0.0, 0.0, 1.0}, .epoch = 100.0, .decay = 5000.0,
                            .attitude = myb_attitude_orbital, .attitude_context = &orbital};
  myb_nongrav_term outgassing = {.env = env, .source = i_sun, .model = comet};
  myb_tidal_lag_term tide = {.env = env, .body = i_planet, .gm_object = 4.9e12, .k2 = 0.3,
                             .lag = 600.0, .with_static = 1};
  for (int k = 0; k < 10; ++k) {
    const double t = 300.0 * k;
    const double x[3] = {7.0e6 * (repro_uniform(rng) - 0.5), 6.9e6 + 1.0e5 * repro_uniform(rng),
                         2.0e6 * (repro_uniform(rng) - 0.5)};
    const double v[3] = {-7.4e3, 500.0 * (repro_uniform(rng) - 0.5), 1.5e3 * repro_uniform(rng)};
    double a[3];
    double dadx[9];
    double dadv[9];
    double m[9];
    MT_CHECK(myb_term_radiation(&moving, t, x, v, a) == MYB_OK);
    hash3(h, a);
    MT_CHECK(myb_term_radiation_jacobian(&moving, t, x, v, dadx, dadv) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, dadx[i]);
      hash_double(h, dadv[i]);
    }
    MT_CHECK(myb_term_surface_radiation(&panels, t, x, v, a) == MYB_OK);
    hash3(h, a);
    MT_CHECK(myb_term_surface_radiation(&hull, t, x, v, a) == MYB_OK);
    hash3(h, a);
    MT_CHECK(myb_term_surface_drag(&flow, t, x, v, a) == MYB_OK);
    hash3(h, a);
    MT_CHECK(myb_term_drag(&drag, t, x, v, a) == MYB_OK);
    hash3(h, a);
    MT_CHECK(myb_term_drag_jacobian(&drag, t, x, v, dadx, dadv) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, dadv[i]);
    }
    MT_CHECK(myb_term_recoil(&recoil, t, x, v, a) == MYB_OK);
    hash3(h, a);
    MT_CHECK(myb_term_nongrav(&outgassing, t, x, v, a) == MYB_OK);
    hash3(h, a);
    MT_CHECK(myb_term_tidal_lag(&tide, t, x, v, a) == MYB_OK);
    hash3(h, a);
    MT_CHECK(myb_attitude_sun_nadir(&yaw, t, x, v, m) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, m[i]);
    }
    MT_CHECK(myb_attitude_spin(&turning, t, x, v, m) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, m[i]);
    }
  }
  /* an orbit under some of them, with two impulses */
  myb_point_mass_term central = {.gm = 3.986004418e14};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "central", myb_term_point_mass,
                               myb_term_point_mass_jacobian, &central, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "light", myb_term_radiation, myb_term_radiation_jacobian,
                               &moving, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "panels", myb_term_surface_radiation, NULL, &panels,
                               NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "recoil", myb_term_recoil, NULL, &recoil, NULL) == MYB_OK);
  const double y0[6] = {7.2e6, 0.0, 1.0e5, 0.0, 7.3e3, 1.2e3};
  for (int variational = 0; variational < 2; ++variational) {
    const myb_orbit_options options = {.epsilon = 1e-10, .variational = variational};
    myb_orbit* orbit = NULL;
    MT_CHECK(myb_orbit_create(model, &options, NULL, &orbit) == MYB_OK);
    MT_CHECK(myb_orbit_start(orbit, 0.0, y0, 9000.0) == MYB_OK);
    double y[6];
    double phi[36];
    const double first[3] = {0.0, 25.0, -3.0};
    const double second[3] = {-4.0, 10.0, 1.0};
    MT_CHECK(myb_orbit_kick(orbit, first) == MYB_OK);
    MT_CHECK(myb_orbit_state(orbit, 3100.0, y) == MYB_OK);
    hash_double(h, y[1]);
    MT_CHECK(myb_orbit_kick(orbit, second) == MYB_OK);
    if (variational) {
      MT_CHECK(myb_orbit_partials(orbit, 9000.0, y, phi, NULL) == MYB_OK);
      for (int i = 0; i < 36; ++i) {
        hash_double(h, phi[i]);
      }
    } else {
      MT_CHECK(myb_orbit_state(orbit, 9000.0, y) == MYB_OK);
    }
    for (int i = 0; i < 6; ++i) {
      hash_double(h, y[i]);
    }
    hash_double(h, (double)myb_orbit_steps(orbit));
    myb_orbit_destroy(orbit);
  }
  myb_force_model_destroy(model);
  myb_env_destroy(env);
  myb_polyhedron_destroy(craft);
}

/* ---- What phases 12-17 had left undone, closed before phase 18 ------------------------------- */

static myb_status debts_turning(void* context, myb_time t, myb_orientation* o) {
  const double* rate = context;
  const double s = myb_time_seconds(t);
  const double angles[3] = {0.3 + 1.0e-7 * s, 0.4 + 0.01 * myb_sin(4.0e-5 * s), 1.0 + *rate * s};
  const double rates[3] = {1.0e-7, 0.01 * 4.0e-5 * myb_cos(4.0e-5 * s), *rate};
  memset(o, 0, sizeof *o);
  myb_rotation_euler313(angles, rates, &o->rotation);
  o->ut = t;
  return MYB_OK;
}

static void workload_debts(hasher* h, mt_rng* rng) {
  const double c = 299792458.0;
  const double au = 1.495978707e11;
  const double gm_sun = 1.32712440018e20;
  const double pi = 3.141592653589793;
  static unsigned char image[16384];
  static char text[16384];

  /* forces: an outgassing that lags, the Yarkovsky effect of a conducting body — its answer,
   * the diurnal force, the seasonal one on an eccentric orbit —, the torque of a tide */
  myb_nongrav comet;
  myb_nongrav_marsden(1.0e-8, 2.0e-9, -3.0e-10, &comet);
  for (int k = 0; k < 8; ++k) {
    const double r[3] = {au * (0.8 + 2.0 * repro_uniform(rng)), au * (repro_uniform(rng) - 0.5),
                         0.1 * au * repro_uniform(rng)};
    const double rn = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    const double speed = sqrt(gm_sun / rn);
    const double v[3] = {-0.9 * speed * r[1] / rn, 0.9 * speed * r[0] / rn, 0.05 * speed};
    double a[3] = {0.0, 0.0, 0.0};
    MT_CHECK(myb_force_nongrav_delayed(&comet, gm_sun, 2.6e6 * (repro_uniform(rng) - 0.3), r, v,
                                       a) == MYB_OK);
    hash3(h, a);
    double re = 0.0;
    double im = 0.0;
    myb_yarkovsky_response(k < 2 ? 0.0 : 0.05 * k * k * (1.0 + repro_uniform(rng)),
                           0.2 + 3.0 * repro_uniform(rng), &re, &im);
    hash_double(h, re);
    hash_double(h, im);
    const myb_yarkovsky body = {.radius = 0.5 + 20.0 * repro_uniform(rng), .density = 2500.0,
                                .absorptivity = 0.9, .emissivity = 0.9,
                                .thermal_inertia = 100.0 + 300.0 * repro_uniform(rng),
                                .rotation_rate = 1.0e-3, .conductivity = k % 2 ? 1.5 : 0.0};
    const double spin[3] = {0.48, -0.6, 0.64};
    const double u[3] = {-r[0] / rn, -r[1] / rn, -r[2] / rn};
    MT_CHECK(myb_force_yarkovsky(&body, spin, 1361.0 * (au / rn) * (au / rn), c, u, a) == MYB_OK);
    hash3(h, a);
    const double state[6] = {r[0], r[1], r[2], v[0], v[1], v[2]};
    myb_yarkovsky_season* season = NULL;
    MT_CHECK(myb_yarkovsky_season_create(&body, spin, gm_sun, state, 1361.0, au, c, NULL,
                                         &season) == MYB_OK);
    if (season != NULL) {
      myb_force_yarkovsky_season(season, r, a);
      hash3(h, a);
      hash_double(h, myb_yarkovsky_season_drift(season));
      hash_double(h, (double)myb_yarkovsky_season_harmonics(season));
      myb_yarkovsky_season_destroy(season);
    }
    const double omega[3] = {1.0e-6, -2.0e-6, 7.0e-5};
    const double moon_r[3] = {r[0] * 2.5e-3, r[1] * 2.5e-3, r[2] * 2.5e-3};
    const double moon_v[3] = {v[0] * 0.03, v[1] * 0.03, v[2] * 0.03};
    double torque[3] = {0.0, 0.0, 0.0};
    MT_CHECK(myb_torque_tidal_lag(3.986004418e14, 4.9028e12, 0.3, 6378136.6, 600.0, omega,
                                  moon_r, moon_v, torque) == MYB_OK);
    hash3(h, torque);
  }

  /* shadows: three bodies on the disk, the air of the occulter */
  const double limb[2] = {0.4, 0.2};
  for (int k = 0; k < 8; ++k) {
    const double source[3] = {-au, 1.0e8 * (repro_uniform(rng) - 0.5), 0.0};
    const double bodies[9] = {-4.0e8, 6.0e6 * repro_uniform(rng), 1.0e6 * (repro_uniform(rng) - 0.5),
                              -1.0e8, -2.0e6 * repro_uniform(rng), 0.0,
                              -3.0e7, 3.0e5 * (repro_uniform(rng) - 0.5), 1.0e5};
    const double radii[3] = {6.4e6, 1.7e6, 1.0e5};
    hash_double(h, myb_shadow_fraction_many(source, 6.957e8, bodies, radii, 3, limb));
    const myb_refraction air = {.refractivity = 2.9e-4, .scale_height = 8.0e3,
                                .optical_depth = 0.1 * k};
    hash_double(h, myb_refraction_bending(&air, 6.4e6, 2.0e3 * k));
    const double behind[3] = {-1.4e9, 2.0e6 + 1.5e6 * k, 0.0};
    hash_double(h, myb_shadow_fraction_refracted(source, 6.957e8, behind, 6.4e6, &air, limb));
  }
  const myb_planet_radiation knocke = {.radius = 6378136.6,
                                       .albedo = {0.34, 0.0, 0.10, 0.0, 0.29},
                                       .emissivity = {0.68, 0.0, -0.07, 0.0, -0.18},
                                       .omega = 1.99e-7};
  for (int k = 0; k < 4; ++k) {
    const double r[3] = {7.0e6 + 1.0e6 * repro_uniform(rng), 2.0e6 * (repro_uniform(rng) - 0.5),
                         3.0e6 * repro_uniform(rng)};
    const double sun[3] = {0.8 * au, 0.6 * au, 0.0};
    double dadr[9];
    MT_CHECK(myb_force_planet_radiation_gradient(&knocke, r, sun, 1361.0, 0.7 * k, 0.02, 6,
                                                 dadr) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, dadr[i]);
    }
  }
  for (int k = 0; k < 6; ++k) {
    const double knudsen = myb_gas_knudsen(1.0e-12 * myb_exp(3.0 * k), 0.016, 3.7e-10, 2.0);
    hash_double(h, knudsen);
    hash_double(h, myb_gas_bridge(knudsen));
  }

  /* a lumpy body of 224 faces: rays through the tree of boxes, YORP */
  enum { RINGS = 8, SECTORS = 16, NV = 2 + (RINGS - 1) * SECTORS, NF = 2 * SECTORS * (RINGS - 1) };
  static double vertices[3 * NV];
  static size_t faces[3 * NF];
  size_t nv = 0;
  size_t nf = 0;
  vertices[0] = 0.0;
  vertices[1] = 0.0;
  vertices[2] = 100.0;
  nv = 1;
  for (int i = 1; i < RINGS; ++i) {
    const double theta = pi * (double)i / (double)RINGS;
    for (int j = 0; j < SECTORS; ++j) {
      const double phi = 2.0 * pi * (double)j / (double)SECTORS;
      const double radius = 100.0 * (1.0 + 0.25 * myb_sin(3.0 * theta) * myb_cos(2.0 * phi + 0.5));
      vertices[3 * nv] = radius * myb_sin(theta) * myb_cos(phi);
      vertices[3 * nv + 1] = radius * myb_sin(theta) * myb_sin(phi);
      vertices[3 * nv + 2] = radius * myb_cos(theta);
      ++nv;
    }
  }
  vertices[3 * nv] = 0.0;
  vertices[3 * nv + 1] = 0.0;
  vertices[3 * nv + 2] = -100.0;
  const size_t south = nv++;
  for (size_t j = 0; j < SECTORS; ++j) {
    const size_t next = (j + 1) % SECTORS;
    const size_t last = 1 + (size_t)(RINGS - 2) * SECTORS;
    faces[3 * nf] = 0;
    faces[3 * nf + 1] = 1 + j;
    faces[3 * nf + 2] = 1 + next;
    ++nf;
    faces[3 * nf] = south;
    faces[3 * nf + 1] = last + next;
    faces[3 * nf + 2] = last + j;
    ++nf;
  }
  for (size_t i = 0; i + 2 < RINGS; ++i) {
    for (size_t j = 0; j < SECTORS; ++j) {
      const size_t next = (j + 1) % SECTORS;
      const size_t a0 = 1 + i * SECTORS + j;
      const size_t b0 = 1 + i * SECTORS + next;
      const size_t c0 = 1 + (i + 1) * SECTORS + j;
      const size_t d0 = 1 + (i + 1) * SECTORS + next;
      faces[3 * nf] = a0;
      faces[3 * nf + 1] = c0;
      faces[3 * nf + 2] = d0;
      ++nf;
      faces[3 * nf] = a0;
      faces[3 * nf + 1] = d0;
      faces[3 * nf + 2] = b0;
      ++nf;
    }
  }
  MT_CHECK(nv == NV && nf == NF);
  myb_polyhedron* rock = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, NV, faces, NF, 1.0, NULL, &rock) == MYB_OK);
  for (int k = 0; k < 12; ++k) {
    double origin[3];
    double direction[3];
    for (int i = 0; i < 3; ++i) {
      origin[i] = 600.0 * (repro_uniform(rng) - 0.5);
    }
    for (int i = 0; i < 3; ++i) {
      direction[i] = 140.0 * (repro_uniform(rng) - 0.5) - origin[i];
    }
    double distance = 0.0;
    double point[3] = {0.0, 0.0, 0.0};
    size_t face = 0;
    const int hit = myb_polyhedron_ray(rock, origin, direction, &distance, point, &face);
    hash_double(h, hit ? distance : -1.0);
    hash_double(h, hit ? (double)face : -1.0);
    hash3(h, point);
  }
  double yorp[3];
  MT_CHECK(myb_force_yorp(rock, 4.5e-6, 1.1, 0.7, 1, 12, yorp) == MYB_OK);
  hash3(h, yorp);
  MT_CHECK(myb_force_yorp(rock, 4.5e-6, 0.4, 0.0, 0, 8, yorp) == MYB_OK);
  hash3(h, yorp);
  myb_polyhedron_destroy(rock);

  /* the limb and the terminator of an ellipsoid; what an observer sees with the light time */
  const double radii[3] = {6378136.6, 6350000.0, 6300000.0};
  for (int k = 0; k < 4; ++k) {
    const double observer[3] = {2.0e7 * (repro_uniform(rng) + 0.5), -1.5e7 * repro_uniform(rng),
                                9.0e6 * (repro_uniform(rng) - 0.5)};
    double centre[3];
    double major[3];
    double minor[3];
    MT_CHECK(myb_ellipsoid_limb(radii, observer, centre, major, minor) == MYB_OK);
    hash3(h, centre);
    hash3(h, major);
    hash3(h, minor);
    const double sun[3] = {0.6 * au, au * (0.8 - 0.3 * repro_uniform(rng)), 0.1 * au};
    double points[9];
    MT_CHECK(myb_ellipsoid_terminator(radii, sun, 6.957e8, k % 2, 3, points) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      hash_double(h, points[i]);
    }
  }
  p12_mover sun_motion = {.radius = au, .rate = 1.99e-7, .phase = 0.3, .inclination = 0.4};
  p12_spinner spinner = {.angles = {0.1, 0.3, 0.5}, .rate = 7.29e-5};
  const myb_env_options env_options = {0};
  myb_env* env = NULL;
  MT_CHECK(myb_env_create(&env_options, NULL, &env) == MYB_OK);
  const myb_env_body planet = {.id = 399, .gm = 3.986004418e14,
                               .orientation = p12_spinner_orientation,
                               .orientation_context = &spinner, .radius = 6378136.6};
  const myb_env_body sun_body = {.id = 10, .gm = gm_sun, .state = p12_mover_state,
                                 .state_context = &sun_motion, .radius = 6.957e8,
                                 .luminosity = 3.828e26};
  size_t i_planet = 0;
  size_t i_sun = 0;
  MT_CHECK(myb_env_add_body(env, &planet, &i_planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &sun_body, &i_sun) == MYB_OK);
  for (int k = 0; k < 4; ++k) {
    const double y[6] = {7.0e6 + 3.0e7 * repro_uniform(rng), 2.0e6 * (repro_uniform(rng) - 0.5),
                         1.0e6 * repro_uniform(rng), -1.0e3, 7.0e3 * repro_uniform(rng), 2.0e3};
    double position[3];
    myb_orientation o;
    double light_time = 0.0;
    MT_CHECK(myb_view_body(env, 100.0 * k, i_sun, y, c, 1, position, &o, &light_time) == MYB_OK);
    hash3(h, position);
    hash_double(h, light_time);
    const double look[3] = {-y[0], -y[1], 1.0e6 - y[2]};
    int hit = 0;
    double point[3] = {0.0, 0.0, 0.0};
    double tau = 0.0;
    MT_CHECK(myb_view_intercept(env, 100.0 * k, i_planet, radii, y, look, c, &hit, point, &tau) ==
             MYB_OK);
    hash_double(h, (double)hit);
    hash3(h, point);
    hash_double(h, tau);
  }
  /* the Jacobians of the terms that turn with the object, and of the Yarkovsky term */
  {
    myb_plate plates[3];
    memset(plates, 0, sizeof plates);
    plates[0] = (myb_plate){.area = 2.0, .normal = {1.0, 0.0, 0.0}, .specular = 0.2,
                            .diffuse = 0.1, .thermal = 1};
    plates[1] = (myb_plate){.area = 3.0, .normal = {0.0, -1.0, 0.0}, .specular = 0.4};
    plates[2] = (myb_plate){.area = 8.0, .normal = {0.6, 0.0, 0.8}, .diffuse = 0.3};
    const myb_object object = {.mass = 400.0, .area = 6.0, .cr = 1.3, .cd = 2.2};
    myb_orbital_attitude orbital = {.env = env, .body = i_planet};
    myb_sun_attitude yaw = {.env = env, .body = i_planet, .source = i_sun};
    myb_surface_radiation_term light = {.env = env, .source = i_sun, .object = &object, .c = c,
                                        .attitude = myb_attitude_orbital,
                                        .attitude_context = &orbital, .plates = plates,
                                        .plate_count = 3};
    myb_surface_drag_term drag = {.env = env, .body = i_planet, .density = p17_density,
                                  .wind = p17_wind, .object = &object,
                                  .gas = {.temperature = 900.0, .molar_mass = 0.016,
                                          .wall_temperature = 300.0, .accommodation = 0.85},
                                  .attitude = myb_attitude_sun_nadir, .attitude_context = &yaw,
                                  .plates = plates, .plate_count = 3, .length = 2.0,
                                  .molecule = 3.7e-10, .continuum_cd = 1.1,
                                  .continuum_area = 6.0};
    myb_recoil_term recoil = {.object = &object, .c = c, .power = 60.0, .efficiency = 0.9,
                              .direction = {0.6, 0.0, 0.8}, .attitude = myb_attitude_orbital,
                              .attitude_context = &orbital};
    const myb_yarkovsky stone = {.radius = 3.0, .density = 2500.0, .absorptivity = 0.9,
                                .emissivity = 0.9, .thermal_inertia = 300.0,
                                .rotation_rate = 1.0e-3, .conductivity = 1.5};
    myb_yarkovsky_term warm = {.env = env, .source = i_sun, .body = stone,
                               .spin = {0.48, -0.6, 0.64}, .c = c};
    for (int k = 0; k < 3; ++k) {
      const double x[3] = {7.0e6 * (0.6 + 0.3 * repro_uniform(rng)),
                           7.0e6 * (0.5 + 0.3 * repro_uniform(rng)),
                           7.0e6 * 0.4 * (repro_uniform(rng) - 0.5)};
      const double v[3] = {-5.0e3 * (0.5 + repro_uniform(rng)), 5.0e3, 2.0e3 * repro_uniform(rng)};
      double dadx[9];
      double dadv[9];
      MT_CHECK(myb_term_surface_radiation_jacobian(&light, 50.0 * k, x, v, dadx, dadv) == MYB_OK);
      for (int i = 0; i < 9; ++i) {
        hash_double(h, dadx[i]);
        hash_double(h, dadv[i]);
      }
      MT_CHECK(myb_term_surface_drag_jacobian(&drag, 50.0 * k, x, v, dadx, dadv) == MYB_OK);
      for (int i = 0; i < 9; ++i) {
        hash_double(h, dadx[i]);
        hash_double(h, dadv[i]);
      }
      MT_CHECK(myb_term_recoil_jacobian(&recoil, 50.0 * k, x, v, dadx, dadv) == MYB_OK);
      for (int i = 0; i < 9; ++i) {
        hash_double(h, dadx[i]);
        hash_double(h, dadv[i]);
      }
      /* the season of the orbit about the Sun that this state is on */
      double sun_state[6];
      MT_CHECK(myb_env_state(env, 50.0 * k, i_sun, sun_state) == MYB_OK);
      const double about_sun[6] = {x[0] - sun_state[0], x[1] - sun_state[1],
                                   x[2] - sun_state[2], v[0] - sun_state[3],
                                   v[1] - sun_state[4], v[2] - sun_state[5]};
      myb_yarkovsky_season* season = NULL;
      MT_CHECK(myb_yarkovsky_season_create(&stone, warm.spin, gm_sun, about_sun, 1361.0, au, c,
                                           NULL, &season) == MYB_OK);
      warm.season = season;
      MT_CHECK(myb_term_yarkovsky_jacobian(&warm, 50.0 * k, x, v, dadx, dadv) == MYB_OK);
      for (int i = 0; i < 9; ++i) {
        hash_double(h, dadx[i]);
      }
      warm.season = NULL;
      myb_yarkovsky_season_destroy(season);
    }
  }
  myb_env* twin = NULL;
  MT_CHECK(myb_env_clone(env, NULL, &twin) == MYB_OK);
  if (twin != NULL) {
    double y[6];
    MT_CHECK(myb_env_state(twin, 250.0, i_sun, y) == MYB_OK);
    hash3(h, y);
    hash3(h, y + 3);
    myb_env_destroy(twin);
  }
  myb_env_destroy(env);

  /* an orbit by DOP853 */
  myb_point_mass_term central = {.gm = 3.986004418e14};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "planet", myb_term_point_mass, myb_term_point_mass_jacobian,
                               &central, NULL) == MYB_OK);
  const myb_orbit_options orbit_options = {.integrator = MYB_ORBIT_DOP853, .rtol = 1e-11,
                                           .atol = 1e-6};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, &orbit_options, NULL, &orbit) == MYB_OK);
  const double y0[6] = {7.0e6, 0.0, 0.0, 0.0, 7.0e3, 2.0e3};
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, 6000.0) == MYB_OK);
  for (int k = 1; k <= 3; ++k) {
    double y[6];
    MT_CHECK(myb_orbit_state(orbit, 2000.0 * k, y) == MYB_OK);
    hash3(h, y);
    hash3(h, y + 3);
  }
  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);

  /* files: SP3 with clocks and accuracy codes */
  {
    static const char* const ids[2] = {"G01", "L51"};
    const myb_calendar labels[2] = {{2026, 1, 1, 0, 0, 0.0}, {2026, 1, 1, 0, 15, 0.0}};
    double positions[12];
    double velocities[12];
    double clocks[4];
    double rates[4];
    double sigmas[32];
    unsigned char flags[4];
    for (int k = 0; k < 12; ++k) {
      positions[k] = 4.0e7 * (repro_uniform(rng) - 0.5);
      velocities[k] = 6.0e3 * (repro_uniform(rng) - 0.5);
    }
    for (int k = 0; k < 4; ++k) {
      clocks[k] = 2.0e-4 * (repro_uniform(rng) - 0.5);
      rates[k] = 2.0e-10 * (repro_uniform(rng) - 0.5);
      flags[k] = (unsigned char)(k * 5 % 16);
      for (int i = 0; i < 3; ++i) {
        sigmas[8 * k + i] = 0.01 + 0.1 * repro_uniform(rng);
        sigmas[8 * k + 4 + i] = 1.0e-5 * (1.0 + repro_uniform(rng));
      }
      sigmas[8 * k + 3] = 1.0e-10 * (1.0 + repro_uniform(rng));
      sigmas[8 * k + 7] = 1.0e-14 * (1.0 + repro_uniform(rng));
    }
    const double accuracy[2] = {0.05, 0.0};
    const myb_sp3_data data = {.version = 'd', .file_type = 'M', .time_system = "GPS",
                               .frame = "IGS20", .orbit_type = "FIT", .agency = "MYB",
                               .data_used = "ORBIT", .satellite_count = 2, .satellites = ids,
                               .epoch_count = 2, .epochs = labels, .interval = 900.0,
                               .positions = positions, .velocities = velocities,
                               .clocks = clocks, .clock_rates = rates, .accuracy = accuracy,
                               .sigmas = sigmas, .flags = flags};
    size_t length = 0;
    MT_CHECK(myb_sp3_write(&data, text, sizeof text, &length) == MYB_OK);
    hash_bytes(h, text, length);
    myb_sp3* sp3 = NULL;
    MT_CHECK(myb_sp3_load(text, length, NULL, &sp3, NULL) == MYB_OK);
    if (sp3 != NULL) {
      double clock = 0.0;
      double rate = 0.0;
      double sigma[8];
      hash_double(h, (double)myb_sp3_clock(sp3, 1, 1, &clock, &rate));
      hash_double(h, clock);
      hash_double(h, rate);
      myb_sp3_sigma(sp3, 0, 1, sigma);
      for (int i = 0; i < 8; ++i) {
        hash_double(h, sigma[i]);
      }
      hash_double(h, myb_sp3_accuracy(sp3, 0));
      hash_double(h, (double)myb_sp3_flags(sp3, 1, 0));
      myb_sp3_destroy(sp3);
    }
  }
  /* a binary PCK from a provider; an SPK segment of two-line elements, word for word */
  {
    double rate = 7.2921e-5;
    myb_pck_writer* writer = NULL;
    MT_CHECK(myb_pck_writer_create("repro", NULL, &writer) == MYB_OK);
    double error = 0.0;
    MT_CHECK(myb_pck_writer_add_chebyshev(writer, 31099, 17, "turning", debts_turning, &rate, 0.0,
                                          86400.0, 6, 11, &error) == MYB_OK);
    hash_double(h, error);
    const size_t bytes = myb_pck_writer_size(writer);
    MT_CHECK(bytes > 0 && bytes <= sizeof image);
    MT_CHECK(myb_pck_writer_write(writer, image, sizeof image) == MYB_OK);
    hash_bytes(h, image, bytes);
    myb_pck_writer_destroy(writer);

    myb_sgp4_constants wgs72;
    myb_sgp4_constants_of(MYB_SGP4_WGS72, &wgs72);
    double constants[8];
    myb_sgp4_spk10_constants(&wgs72, constants);
    double packets[3 * 14];
    for (int k = 0; k < 3; ++k) {
      for (int i = 0; i < 14; ++i) {
        packets[14 * k + i] = repro_uniform(rng);
      }
      packets[14 * k + 9] = 1.0e8 + 5.0e4 * k + repro_uniform(rng);
    }
    double words[128];
    const size_t count = myb_sgp4_spk10_size(3);
    MT_CHECK(count == 73);
    MT_CHECK(myb_sgp4_spk10_words(constants, packets, 3, words, 128) == MYB_OK);
    myb_spk_writer* spk = NULL;
    MT_CHECK(myb_spk_writer_create("repro", NULL, &spk) == MYB_OK);
    const myb_spk_segment_header header = {.target = -77, .center = 399, .frame = 1,
                                           .name = "sets"};
    MT_CHECK(myb_spk_writer_add_raw(spk, &header, 10, 1.0e8, 1.0e8 + 1.5e5, words, count) ==
             MYB_OK);
    const size_t spk_bytes = myb_spk_writer_size(spk);
    MT_CHECK(spk_bytes > 0 && spk_bytes <= sizeof image);
    MT_CHECK(myb_spk_writer_write(spk, image, sizeof image) == MYB_OK);
    hash_bytes(h, image, spk_bytes);
    myb_spk_writer_destroy(spk);
  }
  /* the messages of the CCSDS: covariances, the two notations, the comprehensive message */
  {
    const myb_calendar epochs[2] = {{2026, 3, 1, 12, 0, 0.0}, {2026, 3, 1, 12, 1, 0.5}};
    double states[12];
    double covariance[36];
    for (int k = 0; k < 12; ++k) {
      states[k] = (k % 6 < 3 ? 7.0e6 : 7.0e3) * (repro_uniform(rng) - 0.5);
    }
    for (int row = 0; row < 6; ++row) {
      for (int column = 0; column <= row; ++column) {
        const double value = (row == column ? 100.0 : 1.0) * repro_uniform(rng) *
                             (row < 3 ? 1.0 : 1.0e-3) * (column < 3 ? 1.0 : 1.0e-3);
        covariance[6 * row + column] = value;
        covariance[6 * column + row] = value;
      }
    }
    const myb_oem_segment_data segment = {.object_name = "REPRO <1>", .object_id = "2026-001A",
                                          .center_name = "EARTH", .ref_frame = "EME2000",
                                          .time_system = "UTC", .count = 2, .epochs = epochs,
                                          .states = states, .covariance_count = 1,
                                          .covariance_epochs = epochs, .covariance_frame = "RTN",
                                          .covariances = covariance};
    const myb_oem_data data = {.originator = "myb", .creation_date = {2026, 3, 2, 0, 0, 0.0},
                               .segment_count = 1, .segments = &segment};
    for (int notation = 0; notation < 2; ++notation) {
      size_t length = 0;
      MT_CHECK(myb_oem_write_as(&data, (myb_odm_notation)notation, text, sizeof text, &length) ==
               MYB_OK);
      hash_bytes(h, text, length);
      myb_oem* oem = NULL;
      MT_CHECK(myb_oem_load(text, length, NULL, &oem, NULL) == MYB_OK);
      if (oem != NULL) {
        double p[36];
        double y[6];
        myb_oem_covariance(oem, 0, 0, NULL, NULL, p);
        (void)myb_oem_state(oem, 0, 1, NULL, y, NULL);
        for (int i = 0; i < 36; ++i) {
          hash_double(h, p[i]);
        }
        hash3(h, y);
        hash3(h, y + 3);
        myb_oem_destroy(oem);
      }
    }
    myb_opm opm;
    memset(&opm, 0, sizeof opm);
    memcpy(opm.originator, "myb", 4);
    memcpy(opm.object_name, "REPRO", 6);
    memcpy(opm.object_id, "2026-001A", 10);
    memcpy(opm.center_name, "EARTH", 6);
    memcpy(opm.ref_frame, "EME2000", 8);
    memcpy(opm.time_system, "UTC", 4);
    memcpy(opm.cov_ref_frame, "RTN", 4);
    opm.creation_date = data.creation_date;
    opm.epoch = epochs[1];
    memcpy(opm.state, states + 6, sizeof opm.state);
    opm.has_covariance = 1;
    memcpy(opm.covariance, covariance, sizeof covariance);
    size_t length = 0;
    MT_CHECK(myb_opm_write_as(&opm, MYB_ODM_XML, text, sizeof text, &length) == MYB_OK);
    hash_bytes(h, text, length);
    myb_opm back;
    MT_CHECK(myb_opm_parse(text, length, &back, NULL) == MYB_OK);
    for (int i = 0; i < 36; ++i) {
      hash_double(h, back.covariance[i]);
    }
    char lines[2][256];
    const double first[6] = {7000.0, -0.125, 1.0e-9, 7.5, 0.0, -3.25};
    MT_CHECK(myb_ocm_format_line(NULL, 0.0, first, 6, lines[0], sizeof lines[0], NULL) == MYB_OK);
    MT_CHECK(myb_ocm_format_line(&epochs[1], 0.0, states, 6, lines[1], sizeof lines[1], NULL) ==
             MYB_OK);
    const char* const trajectory[2] = {lines[0], lines[1]};
    const myb_ocm_pair meta[2] = {{"TIME_SYSTEM", "UTC"}, {"EPOCH_TZERO", "2026-03-01T12:00:00"}};
    const myb_ocm_pair traj_pairs[2] = {{"TRAJ_TYPE", "CARTPV"}, {"CENTER_NAME", "EARTH"}};
    const myb_ocm_pair user_pairs[1] = {{"USER_DEFINED_NOTE", "a < b"}};
    const myb_ocm_block_data blocks[2] = {
        {.kind = MYB_OCM_TRAJ, .pair_count = 2, .pairs = traj_pairs, .line_count = 2,
         .lines = trajectory},
        {.kind = MYB_OCM_USER, .pair_count = 1, .pairs = user_pairs}};
    const myb_ocm_data message = {.originator = "myb", .creation_date = data.creation_date,
                                  .meta_count = 2, .meta = meta, .block_count = 2,
                                  .blocks = blocks};
    for (int notation = 0; notation < 2; ++notation) {
      MT_CHECK(myb_ocm_write(&message, (myb_odm_notation)notation, text, sizeof text, &length) ==
               MYB_OK);
      hash_bytes(h, text, length);
      myb_ocm* ocm = NULL;
      MT_CHECK(myb_ocm_load(text, length, NULL, &ocm, NULL) == MYB_OK);
      if (ocm != NULL) {
        double y[9];
        int absolute = 0;
        int count = 0;
        double offset = 0.0;
        MT_CHECK(myb_ocm_state(ocm, 0, 1, &absolute, NULL, &offset, y, &count) == MYB_OK);
        hash_double(h, (double)absolute);
        hash_double(h, (double)count);
        hash3(h, y);
        hash3(h, y + 3);
        myb_ocm_destroy(ocm);
      }
    }
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
  if (upto >= 10) {
    workload_phase10(&h, &rng);
  }
  if (upto >= 11) {
    workload_phase11(&h, &rng);
  }
  if (upto >= 12) {
    workload_phase12(&h, &rng);
  }
  if (upto >= 14) {
    workload_phase14(&h, &rng);
  }
  if (upto >= 15) {
    workload_phase15(&h, &rng);
  }
  if (upto >= 16) {
    workload_phase16(&h, &rng);
  }
  if (upto >= 17) {
    workload_phase17(&h, &rng);
  }
  if (upto >= 18) { /* closed before phase 18 and counted with it */
    workload_debts(&h, &rng);
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
  } else if (upto == 9) {
    expected = UINT64_C(0x827497389931b8a9);
  } else if (upto == 10) {
    expected = UINT64_C(0xd1bce261642c87ac);
  } else if (upto == 11) {
    expected = UINT64_C(0xc9b659a710498026);
  } else if (upto == 12 || upto == 13) {
    expected = UINT64_C(0x6e21db689ceda9e6);
  } else if (upto == 14) {
    expected = UINT64_C(0xb110642ae757e6e0);
  } else if (upto == 15) {
    expected = UINT64_C(0x9128edc3c0c93fb3);
  } else if (upto == 16) {
    expected = UINT64_C(0x66d97d643b636e28);
  } else if (upto == 17) {
    expected = UINT64_C(0x99a2ea1ff9d038a9);
  }
  printf("repro hash: 0x%016" PRIx64 " over %ld values (golden 0x%016" PRIx64 ")\n", h.state,
         h.values, expected);
  MT_CHECK(h.state == expected);
  return mt_finish("test_repro");
}
