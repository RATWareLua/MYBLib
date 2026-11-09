/* SPDX-License-Identifier: BSD-3-Clause */
/* Wisdom-Holman integrators: exact two-body motion, observed orders of every scheme and
 * corrector in the regime e << 1 (error O(e h^p)), agreement with IAS15 for a Jupiter-Saturn
 * system, energy over long runs, time reversibility, dense output, test bodies, independence
 * of the thread count, error statuses. */
#include <math.h>
#include <stdlib.h>

#include "core/fpguard.h"
#include "mtest.h"
#include "myb/dyn.h"
#include "myb/integ.h"
#include "myb/kepler.h"
#include "myb/num.h"
#include "myb/pool.h"

static const double pi = 3.141592653589793;

typedef struct planet {
  double gm;
  double a;
  double e;
  double i;
  double node;
  double argp;
  double nu;
} planet;

/* Star (GM 1) and planets on heliocentric Kepler orbits, then moved to the barycentre. */
static void build(const planet* planets, size_t count, double* gm, double* x, double* v) {
  gm[0] = 1.0;
  for (int c = 0; c < 3; ++c) {
    x[c] = 0.0;
    v[c] = 0.0;
  }
  for (size_t k = 0; k < count; ++k) {
    const planet* p = &planets[k];
    gm[k + 1] = p->gm;
    const myb_elements el = {.p = p->a * (1.0 - p->e * p->e),
                             .e = p->e,
                             .i = p->i,
                             .raan = p->node,
                             .argp = p->argp,
                             .nu = p->nu};
    myb_pv pv = {};
    MT_CHECK(myb_pv_from_elements(1.0 + p->gm, &el, &pv) == MYB_OK);
    for (int c = 0; c < 3; ++c) {
      x[3 * (k + 1) + (size_t)c] = pv.r[c];
      v[3 * (k + 1) + (size_t)c] = pv.v[c];
    }
  }
  double total = 0.0;
  double cx[3] = {0.0, 0.0, 0.0};
  double cv[3] = {0.0, 0.0, 0.0};
  for (size_t k = 0; k <= count; ++k) {
    total += gm[k];
    for (int c = 0; c < 3; ++c) {
      cx[c] += gm[k] * x[3 * k + (size_t)c];
      cv[c] += gm[k] * v[3 * k + (size_t)c];
    }
  }
  for (size_t k = 0; k <= count; ++k) {
    for (int c = 0; c < 3; ++c) {
      x[3 * k + (size_t)c] -= cx[c] / total;
      v[3 * k + (size_t)c] -= cv[c] / total;
    }
  }
}

/* IAS15 reference state at t_end (x, v, 3n each). */
static void reference(myb_nbody* s, const double* x, const double* v, double t_end, double* xr,
                      double* vr) {
  const size_t n = 3 * myb_nbody_count(s);
  const myb_ias15_options o = {.block = 3};
  myb_ias15* in = NULL;
  MT_CHECK(myb_ias15_create(n, &o, NULL, &in) == MYB_OK);
  double* y = malloc(2 * n * sizeof(double));
  for (size_t i = 0; i < n; ++i) {
    y[i] = x[i];
    y[n + i] = v[i];
  }
  MT_CHECK(myb_ias15_start(in, myb_nbody_acceleration_xv, s, 0.0, y) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(in, t_end) == MYB_OK);
  const double* z = myb_ias15_state(in);
  for (size_t i = 0; i < n; ++i) {
    xr[i] = z[i];
    vr[i] = z[n + i];
  }
  free(y);
  myb_ias15_destroy(in);
}

static myb_wh* make_wh(myb_nbody* s, double h, myb_wh_scheme scheme, int corrector) {
  const myb_wh_options o = {.h = h, .scheme = scheme, .corrector = corrector};
  myb_wh* w = NULL;
  MT_CHECK(myb_wh_create(s, &o, NULL, &w) == MYB_OK);
  return w;
}

static double max_difference(const double* a, const double* b, size_t count) {
  double worst = 0.0;
  for (size_t i = 0; i < count; ++i) {
    worst = fmax(worst, fabs(a[i] - b[i]));
  }
  return worst;
}

/* Position error at t_end (largest component, units of the inner semi-major axis). */
static double wh_error(myb_nbody* s, const double* x, const double* v, double t_end, int steps,
                       myb_wh_scheme scheme, int corrector, const double* xr) {
  const size_t n = 3 * myb_nbody_count(s);
  myb_wh* w = make_wh(s, t_end / steps, scheme, corrector);
  MT_CHECK(myb_wh_start(w, 0.0, x, v) == MYB_OK);
  MT_CHECK(myb_wh_integrate(w, t_end) == MYB_OK);
  double* xw = malloc(n * sizeof(double));
  double* vw = malloc(n * sizeof(double));
  MT_CHECK(myb_wh_state(w, xw, vw) == MYB_OK);
  const double error = max_difference(xw, xr, n);
  free(xw);
  free(vw);
  myb_wh_destroy(w);
  return error;
}

static void test_two_body(void) {
  /* One planet and one test body: the kick vanishes to round-off, WH is the Kepler flow. */
  const planet p[1] = {{.gm = 1e-3, .a = 1.0, .e = 0.4, .i = 0.3, .node = 1.0, .argp = 0.5}};
  double gm[2];
  double x[9];
  double v[9];
  build(p, 1, gm, x, v);
  /* test body: circular at 2.5 around the star, i.e. around the barycentre of both */
  const myb_elements et = {.p = 2.5, .e = 0.0, .i = 0.1};
  myb_pv pt = {};
  MT_CHECK(myb_pv_from_elements(1.0 + 1e-3, &et, &pt) == MYB_OK);
  for (int c = 0; c < 3; ++c) {
    x[6 + c] = pt.r[c];
    v[6 + c] = pt.v[c];
  }
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(2, 1, gm, NULL, &s) == MYB_OK);
  /* The test body feels only the star here, so its Jacobi orbit is exactly Keplerian. */
  MT_CHECK(myb_nbody_set_interaction(s, 2, 1, MYB_INTERACTION_NONE) == MYB_OK);
  const double period = 2.0 * pi / sqrt(1.0 + 1e-3);
  const int steps = 2000;
  myb_wh* w = make_wh(s, period / 20.0, MYB_WH_DKD, 0);
  MT_CHECK(myb_wh_start(w, 0.0, x, v) == MYB_OK);
  for (int k = 0; k < steps; ++k) {
    MT_CHECK(myb_wh_step(w) == MYB_OK);
  }
  double xw[9];
  double vw[9];
  MT_CHECK(myb_wh_state(w, xw, vw) == MYB_OK);
  const double t = myb_wh_time(w);
  /* Relative orbit of the pair (GM 1 + 1e-3). */
  myb_pv rel0 = {};
  for (int c = 0; c < 3; ++c) {
    rel0.r[c] = x[3 + c] - x[c];
    rel0.v[c] = v[3 + c] - v[c];
  }
  myb_pv rel = {};
  MT_CHECK(myb_kepler_propagate(1.0 + 1e-3, &rel0, t, &rel) == MYB_OK);
  double error = 0.0;
  for (int c = 0; c < 3; ++c) {
    error = fmax(error, fabs(xw[3 + c] - xw[c] - rel.r[c]));
  }
  /* The star (and the barycentre of the pair) drives the test body: star position + Kepler. */
  myb_pv test0 = {};
  for (int c = 0; c < 3; ++c) {
    test0.r[c] = x[6 + c] - x[c];
    test0.v[c] = v[6 + c] - v[c];
  }
  /* Round-off only: the phase drifts with the energy random walk of 4000 Kepler drifts, each
   * accurate to a few ulp of the velocity (measured 1.0e-11, IAS15-like 2e-12 per 10 orbits). */
  printf("two bodies, %d WH steps (100 orbits): relative orbit off the Kepler propagator by "
         "%.2e of a\n",
         steps, error);
  MT_CHECK(error < 5e-11);
  const myb_wh_stats st = myb_wh_statistics(w);
  MT_CHECK(st.steps == steps && st.evaluations == steps && st.drifts == 2 * steps);
  myb_wh_destroy(w);
  myb_nbody_destroy(s);
}

/* Two planets far from resonance with small masses e; orders from the errors at 8, 16, 32
 * steps per inner period. */
static void test_orders(void) {
  const double eps = 1e-7;
  const planet p[2] = {
      {.gm = eps, .a = 1.0, .e = 0.05, .i = 0.02, .node = 0.3, .argp = 1.0, .nu = 0.0},
      {.gm = eps, .a = 1.71, .e = 0.04, .i = 0.03, .node = 2.0, .argp = 4.0, .nu = 2.5}};
  double gm[3];
  double x[9];
  double v[9];
  build(p, 2, gm, x, v);
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(3, 0, gm, NULL, &s) == MYB_OK);
  const double t_end = 4.0 * 2.0 * pi;
  double xr[9];
  double vr[9];
  reference(s, x, v, t_end, xr, vr);
  typedef struct case_ {
    const char* name;
    myb_wh_scheme scheme;
    int corrector;
    double order;
  } case_;
  const case_ cases[] = {
      {"DKD", MYB_WH_DKD, 0, 2.0},           {"DKD + corrector 3", MYB_WH_DKD, 3, 4.0},
      {"DKD + corrector 5", MYB_WH_DKD, 5, 6.0}, {"DKD + corrector 7", MYB_WH_DKD, 7, 8.0},
      {"SABA2", MYB_WH_SABA2, 0, 4.0},       {"SABA3", MYB_WH_SABA3, 0, 6.0},
      {"SABA4", MYB_WH_SABA4, 0, 8.0},       {"SBAB1 (KDK)", MYB_WH_SBAB1, 0, 2.0},
      {"SBAB2", MYB_WH_SBAB2, 0, 4.0},       {"SBAB3", MYB_WH_SBAB3, 0, 6.0},
      {"SBAB4", MYB_WH_SBAB4, 0, 8.0},
  };
  for (size_t c = 0; c < sizeof cases / sizeof cases[0]; ++c) {
    double e[3];
    for (int k = 0; k < 3; ++k) {
      e[k] = wh_error(s, x, v, t_end, 4 * (8 << k), cases[c].scheme, cases[c].corrector, xr);
    }
    /* 8 -> 16 steps: at least the order (higher terms still large there); 16 -> 32: the
     * asymptotic order, where the error at 32 steps is still above round-off (1e-13); the
     * sixth- and eighth-order schemes reach round-off at 16 steps already (their tables are
     * checked exactly in tests/mpfr/test_mpfr_wh.c). */
    const double order = log2(e[0] / e[1]);
    const double order_fine = log2(e[1] / e[2]);
    printf("%-18s errors %.2e %.2e %.2e (8/16/32 steps per orbit): observed order %.2f, %.2f\n",
           cases[c].name, e[0], e[1], e[2], order, order_fine);
    MT_CHECK(order > cases[c].order - 0.5);
    if (e[2] > 1e-13 && cases[c].order <= 6.0) {
      MT_CHECK(fabs(order_fine - cases[c].order) < 0.4);
    }
  }
  myb_nbody_destroy(s);
}

/* Jupiter- and Saturn-like planets: WH variants against IAS15 over 100 outer orbits. */
static void test_jupiter_saturn(void) {
  const planet p[2] = {
      {.gm = 9.5e-4, .a = 1.0, .e = 0.048, .i = 0.022, .node = 1.75, .argp = 4.78, .nu = 0.3},
      {.gm = 2.85e-4, .a = 1.835, .e = 0.054, .i = 0.043, .node = 1.98, .argp = 5.92, .nu = 2.1}};
  double gm[3];
  double x[9];
  double v[9];
  build(p, 2, gm, x, v);
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(3, 0, gm, NULL, &s) == MYB_OK);
  const double inner = 2.0 * pi;
  const double t_end = 100.0 * 2.0 * pi * pow(1.835, 1.5);
  double xr[9];
  double vr[9];
  reference(s, x, v, t_end, xr, vr);
  myb_nbody_invariants inv0 = {};
  MT_CHECK(myb_nbody_invariants_of(s, x, v, &inv0) == MYB_OK);
  const int per_orbit = 50;
  const int steps = (int)(t_end / inner * per_orbit);
  double errors[4];
  double energies[4];
  const myb_wh_scheme schemes[4] = {MYB_WH_DKD, MYB_WH_DKD, MYB_WH_SABA4, MYB_WH_SBAB4};
  const int correctors[4] = {0, 17, 0, 0};
  const char* names[4] = {"DKD", "DKD + corrector 17", "SABA4", "SBAB4"};
  for (int k = 0; k < 4; ++k) {
    myb_wh* w = make_wh(s, t_end / steps, schemes[k], correctors[k]);
    MT_CHECK(myb_wh_start(w, 0.0, x, v) == MYB_OK);
    double worst_energy = 0.0;
    for (int part = 1; part <= 20; ++part) {
      MT_CHECK(myb_wh_integrate(w, t_end * part / 20.0) == MYB_OK);
      double xw[9];
      double vw[9];
      MT_CHECK(myb_wh_state(w, xw, vw) == MYB_OK);
      myb_nbody_invariants inv = {};
      MT_CHECK(myb_nbody_invariants_of(s, xw, vw, &inv) == MYB_OK);
      worst_energy = fmax(worst_energy, fabs(myb_dd_to_double(myb_dd_sub(inv.energy, inv0.energy)) /
                                             myb_dd_to_double(inv0.energy)));
      if (part == 20) {
        errors[k] = max_difference(xw, xr, 9);
      }
    }
    const myb_wh_stats st = myb_wh_statistics(w);
    printf("Jupiter-Saturn, 100 Saturn orbits, %d steps per Jupiter orbit, %-18s: position off "
           "IAS15 by %.2e, energy within %.2e (%lld evaluations)\n",
           per_orbit, names[k], errors[k], worst_energy, (long long)st.evaluations);
    MT_CHECK(errors[k] < 1e-3);
    energies[k] = worst_energy;
    myb_wh_destroy(w);
  }
  /* The corrector removes the O(e h^2) energy error; the position error keeps the e^2 h^2
   * part, which no corrector touches, while SABA4/SBAB4 shrink both. */
  MT_CHECK(energies[1] < 0.01 * energies[0]);
  MT_CHECK(errors[1] < 0.5 * errors[0]);
  MT_CHECK(errors[2] < 0.05 * errors[0] && errors[3] < 0.05 * errors[0]);
  myb_nbody_destroy(s);
}

static void test_reversibility_and_dense(void) {
  const planet p[3] = {
      {.gm = 3e-6, .a = 1.0, .e = 0.0167, .i = 0.0, .node = 0.0, .argp = 1.8, .nu = 0.1},
      {.gm = 3.2e-7, .a = 1.524, .e = 0.093, .i = 0.032, .node = 0.86, .argp = 5.0, .nu = 3.0},
      {.gm = 9.5e-4, .a = 5.2, .e = 0.048, .i = 0.022, .node = 1.75, .argp = 4.78, .nu = 1.0}};
  double gm[4];
  double x[15];
  double v[15];
  build(p, 3, gm, x, v);
  /* a test body between Mars and Jupiter */
  const myb_elements et = {.p = 2.7 * (1 - 0.01), .e = 0.1, .i = 0.15, .raan = 1.0, .nu = 2.0};
  myb_pv pt = {};
  MT_CHECK(myb_pv_from_elements(1.0, &et, &pt) == MYB_OK);
  for (int c = 0; c < 3; ++c) {
    x[12 + c] = pt.r[c] + x[c];
    v[12 + c] = pt.v[c] + v[c];
  }
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(4, 1, gm, NULL, &s) == MYB_OK);
  const double h = 2.0 * pi / 60.0;
  myb_wh* w = make_wh(s, h, MYB_WH_DKD, 11);
  MT_CHECK(myb_wh_start(w, 0.0, x, v) == MYB_OK);
  MT_CHECK(myb_wh_integrate(w, 3000.0 * h) == MYB_OK);
  MT_CHECK(myb_wh_time(w) == 3000.0 * h);
  /* Dense output against IAS15 at a time inside a step. */
  const double t_mid = 3000.0 * h + 0.37 * h;
  double xd[15];
  double vd[15];
  MT_CHECK(myb_wh_dense(w, t_mid, xd, vd) == MYB_OK);
  double xr[15];
  double vr[15];
  reference(s, x, v, t_mid, xr, vr);
  const double dense_error = max_difference(xd, xr, 15);
  double xn[15];
  double vn[15];
  MT_CHECK(myb_wh_state(w, xn, vn) == MYB_OK);
  double xr_node[15];
  reference(s, x, v, 3000.0 * h, xr_node, vr);
  const double node_error = max_difference(xn, xr_node, 15);
  printf("Earth-Mars-Jupiter + test body, corrector 11, h = 1/60 year: after 50 years node off "
         "IAS15 by %.2e, dense output at 0.37 h by %.2e\n",
         node_error, dense_error);
  MT_CHECK(dense_error < 4.0 * node_error + 1e-12);
  MT_CHECK(myb_wh_dense(w, t_mid + 2.0 * h, xd, vd) == MYB_E_RANGE);
  /* Back to the start. */
  MT_CHECK(myb_wh_integrate(w, 0.0) == MYB_OK);
  MT_CHECK(myb_wh_time(w) == 0.0);
  MT_CHECK(myb_wh_state(w, xn, vn) == MYB_OK);
  const double back = fmax(max_difference(xn, x, 15), max_difference(vn, v, 15));
  printf("3000 steps forward and back: state returns within %.2e\n", back);
  MT_CHECK(back < 1e-12);
  const myb_wh_stats st = myb_wh_statistics(w);
  MT_CHECK(st.steps == 6000);
  myb_wh_destroy(w);
  myb_nbody_destroy(s);
}

/* Long run: energy error of the corrected kernel stays bounded (no drift). */
static void test_energy(void) {
  const planet p[2] = {
      {.gm = 9.5e-4, .a = 1.0, .e = 0.048, .i = 0.022, .node = 1.75, .argp = 4.78, .nu = 0.3},
      {.gm = 2.85e-4, .a = 1.835, .e = 0.054, .i = 0.043, .node = 1.98, .argp = 5.92, .nu = 2.1}};
  double gm[3];
  double x[9];
  double v[9];
  build(p, 2, gm, x, v);
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(3, 0, gm, NULL, &s) == MYB_OK);
  myb_nbody_invariants inv0 = {};
  MT_CHECK(myb_nbody_invariants_of(s, x, v, &inv0) == MYB_OK);
  const double h = 2.0 * pi / 40.0;
  myb_wh* w = make_wh(s, h, MYB_WH_DKD, 17);
  MT_CHECK(myb_wh_start(w, 0.0, x, v) == MYB_OK);
  double early = 0.0;
  double late = 0.0;
  for (int part = 1; part <= 40; ++part) {
    MT_CHECK(myb_wh_integrate(w, part * 5000.0 * h) == MYB_OK);
    double xw[9];
    double vw[9];
    MT_CHECK(myb_wh_state(w, xw, vw) == MYB_OK);
    myb_nbody_invariants inv = {};
    MT_CHECK(myb_nbody_invariants_of(s, xw, vw, &inv) == MYB_OK);
    const double de = fabs(myb_dd_to_double(myb_dd_sub(inv.energy, inv0.energy)) /
                           myb_dd_to_double(inv0.energy));
    if (part <= 4) {
      early = fmax(early, de);
    } else {
      late = fmax(late, de);
    }
  }
  printf("Jupiter-Saturn, corrector 17, 2e5 steps (5000 Jupiter orbits): relative energy error "
         "%.2e in the first 500 orbits, %.2e after\n",
         early, late);
  MT_CHECK(late < 1e-8);
  myb_wh_destroy(w);
  myb_nbody_destroy(s);
}

static void test_threads(void) {
  /* A star and 160 small bodies: enough pairs for the pool to be used. */
  const size_t m = 161;
  double* gm = malloc(m * sizeof(double));
  double* x = malloc(3 * m * sizeof(double));
  double* v = malloc(3 * m * sizeof(double));
  mt_rng rng = {5};
  planet* p = malloc((m - 1) * sizeof(planet));
  for (size_t k = 0; k < m - 1; ++k) {
    const double u = (double)(mt_next(&rng) >> 11) * 0x1p-53;
    p[k] = (planet){.gm = 1e-9 * (1.0 + u),
                    .a = 1.0 + 2.0 * u,
                    .e = 0.1 * u,
                    .i = 0.05 * u,
                    .node = 6.0 * u,
                    .argp = 3.0 * u,
                    .nu = 6.28 * (double)(mt_next(&rng) >> 11) * 0x1p-53};
  }
  build(p, m - 1, gm, x, v);
  double* finals[2] = {malloc(3 * m * sizeof(double)), malloc(3 * m * sizeof(double))};
  double* vel = malloc(3 * m * sizeof(double));
  for (int run = 0; run < 2; ++run) {
    myb_nbody* s = NULL;
    MT_CHECK(myb_nbody_create(m, 0, gm, NULL, &s) == MYB_OK);
    myb_pool* pool = NULL;
    if (run == 1) {
      MT_CHECK(myb_pool_create(4, NULL, &pool) == MYB_OK);
      myb_nbody_set_pool(s, pool);
    }
    myb_wh* w = make_wh(s, 0.05, MYB_WH_SABA2, 0);
    MT_CHECK(myb_wh_start(w, 0.0, x, v) == MYB_OK);
    MT_CHECK(myb_wh_integrate(w, 10.0) == MYB_OK);
    MT_CHECK(myb_wh_state(w, finals[run], vel) == MYB_OK);
    myb_wh_destroy(w);
    myb_nbody_destroy(s);
    myb_pool_destroy(pool);
  }
  int same = 1;
  for (size_t i = 0; i < 3 * m; ++i) {
    same = same && mt_same(finals[0][i], finals[1][i]);
  }
  MT_CHECK(same);
  free(finals[0]);
  free(finals[1]);
  free(vel);
  free(p);
  free(gm);
  free(x);
  free(v);
}

static void test_errors(void) {
  const double gm[2] = {0.0, 1.0};
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(2, 0, gm, NULL, &s) == MYB_OK);
  myb_wh* w = NULL;
  const myb_wh_options o = {.h = 0.1};
  MT_CHECK(myb_wh_create(s, &o, NULL, &w) == MYB_E_DOMAIN); /* GM_0 = 0 */
  myb_nbody_destroy(s);
  const double gm2[2] = {1.0, 1e-3};
  MT_CHECK(myb_nbody_create(2, 0, gm2, NULL, &s) == MYB_OK);
  const myb_wh_options bad[] = {{.h = 0.0},
                                {.h = 1.0 / 0.0},
                                {.h = 0.1, .scheme = (myb_wh_scheme)8},
                                {.h = 0.1, .corrector = 4},
                                {.h = 0.1, .corrector = 19},
                                {.h = 0.1, .scheme = MYB_WH_SABA2, .corrector = 3}};
  for (size_t k = 0; k < sizeof bad / sizeof bad[0]; ++k) {
    MT_CHECK(myb_wh_create(s, &bad[k], NULL, &w) == MYB_E_DOMAIN);
  }
  w = make_wh(s, 0.1, MYB_WH_DKD, 3);
  const double x[6] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  const double v[6] = {0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
  const double bad_x[6] = {0.0, 0.0, 0.0, 0.0 / 0.0, 0.0, 0.0};
  MT_CHECK(myb_wh_start(w, 0.0, bad_x, v) == MYB_E_DOMAIN);
  MT_CHECK(myb_wh_start(w, 0.0, x, v) == MYB_OK);
  MT_CHECK(myb_wh_integrate(w, 1.0 / 0.0) == MYB_E_DOMAIN);
  myb_wh_destroy(w);
  myb_wh_destroy(NULL);
  myb_nbody_destroy(s);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_errors();
  test_two_body();
  test_orders();
  test_jupiter_saturn();
  test_reversibility_and_dense();
  test_energy();
  test_threads();
  return mt_finish("test_whfast");
}
