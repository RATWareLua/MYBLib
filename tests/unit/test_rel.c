/* SPDX-License-Identifier: BSD-3-Clause */
/* Relativity at 1PN:
 *   - the EIH N-body accelerations against the Euler-Lagrange equations of the global PPN
 *     Lagrangian, derived symbolically and evaluated with 50 digits (tests/data/rel), for
 *     several (beta, gamma), 2 .. 5 bodies; the same bits on every path and thread count;
 *   - the conservation of the 1PN energy and momentum along an integration, at O(c^-4);
 *   - the perihelion advance of Mercury, 6 pi GM / (c^2 a (1 - e^2)) (2 + 2 gamma - beta) / 3
 *     per orbit, 42.98" per century in general relativity;
 *   - the satellite terms of IERS 2010 eq. 10.12: Schwarzschild perigee advance,
 *     Lense-Thirring node rate of LAGEOS (31 mas/yr), de Sitter precession (19.2 mas/yr);
 *   - the BCRS -> GCRS transformation of IAU 2000 B1.3 (Lorentz limit, round trip) and the
 *     proper time of a GPS clock (IERS 2010: 4.46e-10 faster than TT). */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dyn/nbody_internal.h"
#include "mtest.h"
#include "myb/dyn.h"
#include "myb/integ.h"
#include "myb/kepler.h"
#include "myb/pool.h"
#include "myb/rel.h"

static double norm3(const double* v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/* ---- EIH against the Lagrangian ------------------------------------------------------------ */

static void test_oracle(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/rel/eih_reference.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  myb_pool* pool = NULL;
  MT_CHECK(myb_pool_create(5, NULL, &pool) == MYB_OK);
  char line[1024];
  int cases = 0;
  double worst = 0.0;    /* |a - a_ref| / max |1PN part| */
  double smallest = 1.0; /* smallest 1PN part relative to the Newtonian one */
  /* bound: one case per iteration, at most the lines of the file */
  while (fgets(line, sizeof line, f) != NULL) {
    int n = 0;
    double beta = 0.0;
    double gamma = 0.0;
    double c = 0.0;
    if (sscanf(line, "case %d %lf %lf %lf", &n, &beta, &gamma, &c) != 4) {
      continue;
    }
    MT_CHECK(n >= 2 && n <= 8);
    double gm[8];
    double x[24];
    double v[24];
    double ref[24];
    for (int i = 0; i < n; ++i) {
      MT_CHECK(fgets(line, sizeof line, f) != NULL);
      MT_CHECK(sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &gm[i], &x[3 * i],
                      &x[3 * i + 1], &x[3 * i + 2], &v[3 * i], &v[3 * i + 1], &v[3 * i + 2],
                      &ref[3 * i], &ref[3 * i + 1], &ref[3 * i + 2]) == 10);
    }
    myb_nbody* s = NULL;
    myb_nbody* newton = NULL;
    MT_CHECK(myb_nbody_create((size_t)n, 0, gm, NULL, &s) == MYB_OK);
    MT_CHECK(myb_nbody_create((size_t)n, 0, gm, NULL, &newton) == MYB_OK);
    const myb_ppn ppn = {c, beta, gamma};
    MT_CHECK(myb_nbody_set_relativity(s, &ppn) == MYB_OK);
    double a[24];
    double an[24];
    MT_CHECK(myb_nbody_acceleration_xv(s, 0.0, x, v, a) == MYB_OK);
    MT_CHECK(myb_nbody_acceleration_xv(newton, 0.0, x, v, an) == MYB_OK);
    double pn = 0.0;
    double nn = 0.0;
    double err = 0.0;
    for (int i = 0; i < 3 * n; ++i) {
      pn = fmax(pn, fabs(ref[i] - an[i]));
      nn = fmax(nn, fabs(an[i]));
      err = fmax(err, fabs(a[i] - ref[i]));
    }
    worst = fmax(worst, err / pn);
    smallest = fmin(smallest, pn / nn);
    /* every path, with and without threads, gives the same bits */
    const myb_internal_nbody_path paths[3] = {MYB_NBODY_PATH_SYMMETRIC, MYB_NBODY_PATH_RECEIVER,
                                              MYB_NBODY_PATH_RECEIVER_SCALAR};
    for (int k = 0; k < 3; ++k) {
      myb_nbody_set_pool(s, k == 1 ? pool : NULL);
      double b[24];
      MT_CHECK(myb_internal_nbody_acceleration_xv_path(s, x, v, b, paths[k]) == MYB_OK);
      for (int i = 0; i < 3 * n; ++i) {
        MT_CHECK_SAME(b[i], a[i]);
      }
    }
    myb_nbody_destroy(s);
    myb_nbody_destroy(newton);
    ++cases;
  }
  fclose(f);
  myb_pool_destroy(pool);
  printf("EIH vs Euler-Lagrange of the PPN Lagrangian: %d cases, max error %.2e of the 1PN part "
         "(1PN part >= %.1e of Newton)\n",
         cases, worst, smallest);
  MT_CHECK(cases == 48);
  /* rounding of the double evaluation against the 1PN part: c = 1000 makes it 1e-6 of Newton */
  MT_CHECK(worst < 1e-9);
}

/* ---- Graph, errors -------------------------------------------------------------------------- */

static void test_graph(void) {
  const double gm[4] = {1.0, 0.3, 0.2, 0.0};
  const double x[12] = {0, 0, 0, 1, 0.1, 0, -0.4, 1.3, 0.2, 2.0, -1.0, 0.5};
  const double v[12] = {0, 0, 0, 0.1, 0.9, 0, -0.7, -0.2, 0.1, 0.3, 0.6, -0.1};
  myb_nbody* s = NULL;
  myb_nbody* plain = NULL;
  MT_CHECK(myb_nbody_create(3, 1, gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_nbody_create(3, 1, gm, NULL, &plain) == MYB_OK);
  double a[12];
  double b[12];
  MT_CHECK(myb_nbody_set_interaction(s, 0, 1, MYB_INTERACTION_EIH) == MYB_E_DOMAIN);
  const myb_ppn bad = {0.0, 1.0, 1.0};
  MT_CHECK(myb_nbody_set_relativity(s, &bad) == MYB_E_DOMAIN);
  const myb_ppn ppn = {30.0, 1.0, 1.0};
  MT_CHECK(myb_nbody_set_relativity(s, &ppn) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, x, a) == MYB_E_DOMAIN); /* needs velocities */
  MT_CHECK(myb_nbody_acceleration_xv(s, 0.0, x, v, a) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration_xv(plain, 0.0, x, v, b) == MYB_OK);
  int differ = 0;
  for (int i = 0; i < 12; ++i) {
    differ |= a[i] != b[i];
  }
  MT_CHECK(differ);
  /* every pair explicitly Newtonian: the Newtonian bits */
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      if (i != j) {
        MT_CHECK(myb_nbody_set_interaction(s, i, j, MYB_INTERACTION_NEWTON) == MYB_OK);
      }
    }
  }
  MT_CHECK(myb_nbody_acceleration_xv(s, 0.0, x, v, a) == MYB_OK);
  for (int i = 0; i < 12; ++i) {
    MT_CHECK_SAME(a[i], b[i]);
  }
  /* one pair back to EIH changes only its receiver */
  MT_CHECK(myb_nbody_set_interaction(s, 3, 0, MYB_INTERACTION_EIH) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration_xv(s, 0.0, x, v, a) == MYB_OK);
  for (int i = 0; i < 9; ++i) {
    MT_CHECK_SAME(a[i], b[i]);
  }
  MT_CHECK(a[9] != b[9] || a[10] != b[10] || a[11] != b[11]);
  /* an excluded pair stays excluded in the 1PN pass; NULL relativity: Newtonian again */
  MT_CHECK(myb_nbody_set_interaction(s, 3, 0, MYB_INTERACTION_NONE) == MYB_OK);
  MT_CHECK(myb_nbody_set_interaction(plain, 3, 0, MYB_INTERACTION_NONE) == MYB_OK);
  MT_CHECK(myb_nbody_set_relativity(s, NULL) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration(s, 0.0, x, a) == MYB_OK);
  MT_CHECK(myb_nbody_acceleration(plain, 0.0, x, b) == MYB_OK);
  for (int i = 0; i < 12; ++i) {
    MT_CHECK_SAME(a[i], b[i]);
  }
  /* Wisdom-Holman refuses velocity-dependent forces */
  MT_CHECK(myb_nbody_set_relativity(s, &ppn) == MYB_OK);
  myb_wh* wh = NULL;
  const myb_wh_options o = {.h = 0.01, .scheme = MYB_WH_DKD, .corrector = 0};
  MT_CHECK(myb_wh_create(s, &o, NULL, &wh) == MYB_OK);
  MT_CHECK(myb_wh_start(wh, 0.0, x, v) == MYB_E_DOMAIN);
  myb_wh_destroy(wh);
  myb_nbody_pn_invariants inv;
  MT_CHECK(myb_nbody_pn_invariants_of(plain, x, v, &inv) == MYB_E_DOMAIN);
  myb_nbody_destroy(s);
  myb_nbody_destroy(plain);
}

/* ---- Conservation ------------------------------------------------------------------------- */

typedef struct drift {
  double energy;  /* |h(T) - h(0)| / |h(0)| */
  double newton;  /* the same for the Newtonian energy */
  double momentum; /* |P(T) - P(0)| / |P scale| */
  double center;   /* barycentre off its uniform motion / size of the system */
} drift;

static drift conservation(double c) {
  const double gm[3] = {1.0, 0.5, 0.01};
  const double x[9] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.1, 0.0, -2.5, 0.0};
  const double v[9] = {0.0, -0.35, 0.0, 0.0, 0.85, 0.05, 0.8, 0.0, -0.05};
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(3, 0, gm, NULL, &s) == MYB_OK);
  const myb_ppn ppn = {c, 1.0, 1.0};
  MT_CHECK(myb_nbody_set_relativity(s, &ppn) == MYB_OK);
  const myb_ias15_options o = {.epsilon = 1e-10, .block = 3};
  myb_ias15* in = NULL;
  MT_CHECK(myb_ias15_create(9, &o, NULL, &in) == MYB_OK);
  double y[18];
  memcpy(y, x, sizeof x);
  memcpy(y + 9, v, sizeof v);
  MT_CHECK(myb_ias15_start(in, myb_nbody_acceleration_xv, s, 0.0, y) == MYB_OK);
  myb_nbody_pn_invariants i0;
  myb_nbody_invariants n0;
  MT_CHECK(myb_nbody_pn_invariants_of(s, x, v, &i0) == MYB_OK);
  MT_CHECK(myb_nbody_invariants_of(s, x, v, &n0) == MYB_OK);
  const double t_end = 60.0;
  MT_CHECK(myb_ias15_integrate(in, t_end) == MYB_OK);
  const double* z = myb_ias15_state(in);
  myb_nbody_pn_invariants i1;
  myb_nbody_invariants n1;
  MT_CHECK(myb_nbody_pn_invariants_of(s, z, z + 9, &i1) == MYB_OK);
  MT_CHECK(myb_nbody_invariants_of(s, z, z + 9, &n1) == MYB_OK);
  drift d;
  d.energy = fabs((i1.energy.hi - i0.energy.hi) + (i1.energy.lo - i0.energy.lo)) /
             fabs(i0.energy.hi);
  d.newton = fabs((n1.energy.hi - n0.energy.hi) + (n1.energy.lo - n0.energy.lo)) /
             fabs(n0.energy.hi);
  double dp[3];
  double dc[3];
  for (int k = 0; k < 3; ++k) {
    dp[k] = i1.momentum[k] - i0.momentum[k];
    dc[k] = i1.center[k] - (i0.center[k] + i0.center_velocity[k] * t_end);
  }
  d.momentum = norm3(dp) / 0.85;
  d.center = norm3(dc) / 2.5;
  myb_ias15_destroy(in);
  myb_nbody_destroy(s);
  return d;
}

static void test_conservation(void) {
  const drift a = conservation(100.0);
  const drift b = conservation(200.0);
  printf("1PN conservation over 60 time units, v/c ~ 1e-2: energy %.2e (Newtonian energy "
         "%.2e), momentum %.2e, barycentre %.2e; c doubled: %.2e, %.2e, %.2e, %.2e\n",
         a.energy, a.newton, a.momentum, a.center, b.energy, b.newton, b.momentum, b.center);
  /* the 1PN quantities change at O(c^-4), the Newtonian energy at O(c^-2) */
  MT_CHECK(a.energy < 1e-3 * a.newton && b.energy < 1e-3 * b.newton);
  MT_CHECK(a.energy / b.energy > 8.0 && a.newton / b.newton > 3.0);
  MT_CHECK(a.momentum < 1e-6 && a.center < 1e-6);
}

/* ---- Perihelion advance ------------------------------------------------------------------ */

static double radial_velocity(void* context, double t, const double* y) {
  (void)context;
  (void)t;
  double r[3];
  double v[3];
  for (int k = 0; k < 3; ++k) {
    r[k] = y[3 + k] - y[k];
    v[k] = y[9 + k] - y[6 + k];
  }
  return r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
}

/* Perihelion advance per orbit of a planet around its star (rad), measured over `orbits`. */
static double perihelion_rate(double gm_star, double gm_planet, double a, double e,
                              const myb_ppn* ppn, int orbits, double* predicted) {
  const double gm[2] = {gm_star, gm_planet};
  const double mu = gm_star + gm_planet;
  /* start at aphelion on the x axis, moving in the y direction */
  const double ra = a * (1.0 + e);
  const double va = sqrt(mu / a * (1.0 - e) / (1.0 + e));
  const double f = gm_planet / mu;
  const double x[6] = {-f * ra, 0.0, 0.0, (1.0 - f) * ra, 0.0, 0.0};
  const double v[6] = {0.0, -f * va, 0.0, 0.0, (1.0 - f) * va, 0.0};
  myb_nbody* s = NULL;
  MT_CHECK(myb_nbody_create(2, 0, gm, NULL, &s) == MYB_OK);
  MT_CHECK(myb_nbody_set_relativity(s, ppn) == MYB_OK);
  const myb_ias15_options o = {.epsilon = 1e-9, .block = 3};
  myb_ias15* in = NULL;
  MT_CHECK(myb_ias15_create(6, &o, NULL, &in) == MYB_OK);
  double y[12];
  memcpy(y, x, sizeof x);
  memcpy(y + 6, v, sizeof v);
  MT_CHECK(myb_ias15_start(in, myb_nbody_acceleration_xv, s, 0.0, y) == MYB_OK);
  const double period = 2.0 * 3.141592653589793 * sqrt(a * a * a / mu);
  double angles[400];
  int count = 0;
  double previous = 0.0;
  /* bound: the steps of `orbits` revolutions, a few hundred per revolution */
  for (long step = 0; step < 2000L * orbits && count < orbits && count < 400; ++step) {
    MT_CHECK(myb_ias15_step(in, (orbits + 0.75) * period) == MYB_OK);
    int found = 0;
    double t_event = 0.0;
    double ye[12];
    MT_CHECK(myb_ias15_locate(in, radial_velocity, NULL, &found, &t_event, ye) == MYB_OK);
    /* a sign change inside the step ending with r.v > 0: from - to +, a perihelion */
    if (found && radial_velocity(NULL, myb_ias15_time(in), myb_ias15_state(in)) > 0.0) {
      double angle = atan2(ye[4] - ye[1], ye[3] - ye[0]);
      /* unwrap */
      if (count > 0) {
        angle += 6.283185307179586 * round((previous - angle) / 6.283185307179586);
      }
      angles[count++] = angle;
      previous = angle;
    }
  }
  myb_ias15_destroy(in);
  myb_nbody_destroy(s);
  /* least-squares slope of the perihelion longitude per orbit */
  double sx = 0.0;
  double sy = 0.0;
  double sxx = 0.0;
  double sxy = 0.0;
  for (int k = 0; k < count; ++k) {
    sx += k;
    sy += angles[k];
    sxx += (double)k * k;
    sxy += k * angles[k];
  }
  const double slope = (count * sxy - sx * sy) / (count * sxx - sx * sx);
  *predicted = 6.0 * 3.141592653589793 * mu / (ppn->c * ppn->c * a * (1.0 - e * e)) *
               (2.0 + 2.0 * ppn->gamma - ppn->beta) / 3.0;
  MT_CHECK(count == orbits);
  return slope;
}

static void test_perihelion(void) {
  /* Sun and Mercury (GM from DE440, orbit of J2000) */
  const double gm_sun = 1.32712440041279419e20;
  const double gm_mercury = 2.2031868551e13;
  const double a = 5.7909226e10;
  const double e = 0.20563069;
  const double c = 299792458.0;
  const double params[3][2] = {{1.0, 1.0}, {1.5, 0.7}, {0.6, 1.3}};
  for (int k = 0; k < 3; ++k) {
    const myb_ppn ppn = {c, params[k][0], params[k][1]};
    double predicted = 0.0;
    const double rate = perihelion_rate(gm_sun, gm_mercury, a, e, &ppn, 60, &predicted);
    const double orbits_per_century = 36525.0 * 86400.0 /
                                      (2.0 * 3.141592653589793 *
                                       sqrt(a * a * a / (gm_sun + gm_mercury)));
    printf("Mercury, beta %.1f gamma %.1f: perihelion %.4f\"/century (predicted %.4f\"), "
           "relative difference %.1e\n",
           ppn.beta, ppn.gamma, rate * orbits_per_century * 206264.80624709636,
           predicted * orbits_per_century * 206264.80624709636, rate / predicted - 1.0);
    MT_CHECK(fabs(rate / predicted - 1.0) < 1e-6); /* mass ratio of Mercury: 1.7e-7 */
    if (k == 0) {
      MT_CHECK(fabs(rate * orbits_per_century * 206264.80624709636 - 42.98) < 0.01);
    }
  }
}

/* ---- Satellite terms (IERS 2010 eq. 10.12) ------------------------------------------------- */

typedef struct satellite {
  myb_ppn ppn;
  double gm;
  double spin[3];
  int terms; /* 1 Schwarzschild, 2 Lense-Thirring, 4 de Sitter */
  double gm_sun;
  double au;
  double n_sun; /* rad/s of the circular heliocentric orbit */
} satellite;

static myb_status satellite_acceleration(void* context, double t, const double* x,
                                         const double* v, double* a) {
  const satellite* s = context;
  const double r = norm3(x);
  for (int k = 0; k < 3; ++k) {
    a[k] = -s->gm * x[k] / (r * r * r);
  }
  double d[3];
  if (s->terms & 1) {
    MT_CHECK(myb_rel_schwarzschild(&s->ppn, s->gm, x, v, d) == MYB_OK);
    for (int k = 0; k < 3; ++k) {
      a[k] += d[k];
    }
  }
  if (s->terms & 2) {
    MT_CHECK(myb_rel_lense_thirring(&s->ppn, s->gm, s->spin, x, v, d) == MYB_OK);
    for (int k = 0; k < 3; ++k) {
      a[k] += d[k];
    }
  }
  if (s->terms & 4) {
    /* the central body on a circular orbit about the external one, in the x-y plane */
    const double R[3] = {s->au * cos(s->n_sun * t), s->au * sin(s->n_sun * t), 0.0};
    const double V[3] = {-s->au * s->n_sun * sin(s->n_sun * t),
                         s->au * s->n_sun * cos(s->n_sun * t), 0.0};
    MT_CHECK(myb_rel_de_sitter(&s->ppn, s->gm_sun, R, V, v, d) == MYB_OK);
    for (int k = 0; k < 3; ++k) {
      a[k] += d[k];
    }
  }
  return MYB_OK;
}

/* Integrates whole orbits and returns the angular momentum direction change about z: the
 * rotation of its projection onto the x-y plane (rad), and the elapsed time. */
static double node_shift(satellite* s, const myb_pv* start, int orbits, double* elapsed) {
  const double a = 12270e3;
  const double period = 2.0 * 3.141592653589793 * sqrt(a * a * a / s->gm);
  const myb_ias15_options o = {.epsilon = 1e-10};
  myb_ias15* in = NULL;
  MT_CHECK(myb_ias15_create(3, &o, NULL, &in) == MYB_OK);
  double y[6] = {start->r[0], start->r[1], start->r[2], start->v[0], start->v[1], start->v[2]};
  MT_CHECK(myb_ias15_start(in, satellite_acceleration, s, 0.0, y) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(in, orbits * period) == MYB_OK);
  const double* z = myb_ias15_state(in);
  double l0[3];
  double l1[3];
  const double r0[3] = {start->r[0], start->r[1], start->r[2]};
  const double v0[3] = {start->v[0], start->v[1], start->v[2]};
  l0[0] = r0[1] * v0[2] - r0[2] * v0[1];
  l0[1] = r0[2] * v0[0] - r0[0] * v0[2];
  l0[2] = r0[0] * v0[1] - r0[1] * v0[0];
  l1[0] = z[1] * z[5] - z[2] * z[4];
  l1[1] = z[2] * z[3] - z[0] * z[5];
  l1[2] = z[0] * z[4] - z[1] * z[3];
  *elapsed = orbits * period;
  myb_ias15_destroy(in);
  return atan2(l1[1], l1[0]) - atan2(l0[1], l0[0]);
}

static void test_satellite(void) {
  const double year = 365.25 * 86400.0;
  const double mas = 3.141592653589793 / 648000.0 / 1000.0;
  /* LAGEOS: a = 12270 km, e = 0.0045, i = 109.84 deg */
  myb_elements el = {.p = 12270e3 * (1.0 - 0.0045 * 0.0045), .e = 0.0045,
                     .i = 109.84 * 3.141592653589793 / 180.0, .raan = 0.6, .argp = 1.1,
                     .nu = 0.3};
  satellite s = {.ppn = {299792458.0, 1.0, 1.0}, .gm = 3.986004418e14,
                 .spin = {0.0, 0.0, 9.8e8}, .gm_sun = 1.32712440041279419e20,
                 .au = 1.495978707e11};
  s.n_sun = sqrt(s.gm_sun / (s.au * s.au * s.au));
  myb_pv start;
  MT_CHECK(myb_pv_from_elements(s.gm, &el, &start) == MYB_OK);
  double elapsed = 0.0;
  /* Lense-Thirring: node rate (1 + gamma) GM j / (c^2 a^3 (1 - e^2)^(3/2)) */
  s.terms = 2;
  const double lt = node_shift(&s, &start, 600, &elapsed) / elapsed;
  const double a = 12270e3;
  const double c2 = s.ppn.c * s.ppn.c;
  const double lt_expected =
      2.0 * s.gm * s.spin[2] / (c2 * a * a * a * pow(1.0 - 0.0045 * 0.0045, 1.5));
  /* de Sitter: the orbit turns about the pole of the heliocentric orbit at
   * (1 + 2 gamma)/2 GM_sun n / (c^2 au) */
  s.terms = 4;
  const double ds = node_shift(&s, &start, 600, &elapsed) / elapsed;
  const double ds_expected = 1.5 * s.gm_sun * s.n_sun / (c2 * s.au);
  printf("LAGEOS node: Lense-Thirring %.3f mas/yr (expected %.3f), de Sitter %.3f mas/yr "
         "(expected %.3f)\n",
         lt * year / mas, lt_expected * year / mas, ds * year / mas, ds_expected * year / mas);
  MT_CHECK(fabs(lt / lt_expected - 1.0) < 1e-3);
  MT_CHECK(fabs(ds / ds_expected - 1.0) < 1e-3);
  MT_CHECK(fabs(lt_expected * year / mas - 31.0) < 1.0 && fabs(ds_expected * year / mas - 19.2) < 0.1);
}

static double perigee_angle(const double* y, double gm) {
  /* direction of the eccentricity vector in the orbital plane x-y */
  const double r = norm3(y);
  const double* v = y + 3;
  const double rv = y[0] * v[0] + y[1] * v[1] + y[2] * v[2];
  const double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  double e[3];
  for (int k = 0; k < 3; ++k) {
    e[k] = ((v2 - gm / r) * y[k] - rv * v[k]) / gm;
  }
  return atan2(e[1], e[0]);
}

static void test_schwarzschild(void) {
  /* an eccentric orbit in the x-y plane: perigee advance per orbit as for Mercury */
  satellite s = {.ppn = {299792458.0, 1.0, 1.0}, .gm = 3.986004418e14, .terms = 1};
  const double params[2][2] = {{1.0, 1.0}, {1.4, 0.8}};
  for (int k = 0; k < 2; ++k) {
    s.ppn.beta = params[k][0];
    s.ppn.gamma = params[k][1];
    const double a = 12000e3;
    const double e = 0.3;
    const myb_elements el = {.p = a * (1.0 - e * e), .e = e, .i = 0.0, .raan = 0.0, .argp = 0.0,
                             .nu = 3.141592653589793};
    myb_pv start;
    MT_CHECK(myb_pv_from_elements(s.gm, &el, &start) == MYB_OK);
    const double period = 2.0 * 3.141592653589793 * sqrt(a * a * a / s.gm);
    const myb_ias15_options o = {.epsilon = 1e-10};
    myb_ias15* in = NULL;
    MT_CHECK(myb_ias15_create(3, &o, NULL, &in) == MYB_OK);
    double y[6] = {start.r[0], start.r[1], start.r[2], start.v[0], start.v[1], start.v[2]};
    MT_CHECK(myb_ias15_start(in, satellite_acceleration, &s, 0.0, y) == MYB_OK);
    const int orbits = 200;
    MT_CHECK(myb_ias15_integrate(in, orbits * period) == MYB_OK);
    /* the osculating eccentricity vector at apogee (same phase as the start) */
    const double shift = perigee_angle(myb_ias15_state(in), s.gm) - perigee_angle(y, s.gm);
    const double expected = orbits * 6.0 * 3.141592653589793 * s.gm /
                            (s.ppn.c * s.ppn.c * a * (1.0 - e * e)) *
                            (2.0 + 2.0 * s.ppn.gamma - s.ppn.beta) / 3.0;
    printf("Schwarzschild perigee advance (beta %.1f, gamma %.1f): %.6e rad in %d orbits "
           "(expected %.6e)\n",
           s.ppn.beta, s.ppn.gamma, shift, orbits, expected);
    MT_CHECK(fabs(shift / expected - 1.0) < 1e-3);
    myb_ias15_destroy(in);
  }
}

/* ---- Local frame, proper time ------------------------------------------------------------- */

static void test_local_frame(void) {
  const double c = 299792458.0;
  /* Lorentz limit: no external field, no acceleration */
  myb_rel_local_frame f = {.c = c, .xe = {1.4e11, -2.0e10, 7.0e9},
                           .ve = {2.1e4, 1.9e4, -8.0e3}, .ae = {0, 0, 0}, .w_ext = 0.0};
  const double x[3] = {1.4e11 + 6.0e6, -2.0e10 - 3.0e6, 7.0e9 + 1.0e6};
  double X[3];
  myb_rel_global_to_local(&f, x, X);
  double r[3];
  for (int k = 0; k < 3; ++k) {
    r[k] = x[k] - f.xe[k];
  }
  const double v2 = f.ve[0] * f.ve[0] + f.ve[1] * f.ve[1] + f.ve[2] * f.ve[2];
  const double lorentz = 1.0 / sqrt(1.0 - v2 / (c * c));
  const double vr = f.ve[0] * r[0] + f.ve[1] * r[1] + f.ve[2] * r[2];
  double worst = 0.0;
  for (int k = 0; k < 3; ++k) {
    const double exact = r[k] + (lorentz - 1.0) / v2 * vr * f.ve[k];
    worst = fmax(worst, fabs(X[k] - exact));
  }
  MT_CHECK(worst < 1e-9); /* m, for |r| = 6.8e6 m: the O(c^-4) terms are 1e-16 r */
  /* with the external potential and acceleration of the Earth: round trip to O(c^-4) */
  f.w_ext = 8.87e8 + 1.2e5;
  f.ae[0] = -5.9e-3;
  f.ae[1] = 1.0e-3;
  f.ae[2] = 4.0e-4;
  double back[3];
  myb_rel_global_to_local(&f, x, X);
  myb_rel_local_to_global(&f, X, back);
  double trip = 0.0;
  double scale = 0.0;
  for (int k = 0; k < 3; ++k) {
    trip = fmax(trip, fabs(back[k] - x[k]));
    scale = fmax(scale, fabs(X[k] - r[k]));
  }
  printf("BCRS -> GCRS at 6800 km: correction %.3f mm, round trip %.1e m, Lorentz limit %.1e m\n",
         scale * 1e3, trip, worst);
  MT_CHECK(trip < 1e-6 && scale > 5e-2 && scale < 1e-1);
  MT_CHECK(fabs(myb_rel_local_time_shift(&f, x) + vr / (c * c)) < 1e-20);
  /* GPS clock: dtau/dTT = 1 + L_G - (v^2/2 + U_E)/c^2 = 1 + 4.46e-10 (IERS 2010 sect. 10.2) */
  const double lg = 6.969290134e-10;
  const double rgps = 26561.75e3;
  const double gm = 3.986004418e14;
  const double vgps = sqrt(gm / rgps);
  const double rate = myb_rel_proper_time_rate(c, vgps * vgps, gm / rgps) + lg;
  printf("GPS clock: dtau/dTT - 1 = %.4e\n", rate - 1.0);
  MT_CHECK(fabs((rate - 1.0) - 4.465e-10) < 2e-12);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_oracle();
  test_graph();
  test_conservation();
  test_perihelion();
  test_satellite();
  test_schwarzschild();
  test_local_frame();
  return mt_finish("test_rel");
}
