/* SPDX-License-Identifier: BSD-3-Clause */
/* The circular restricted three-body problem against the NASA/JPL Three-Body Periodic Orbits
 * catalogue (tests/data/cr3bp/jpl_periodic_orbits.txt: halo L1/L2, Lyapunov L1/L2, vertical,
 * axial, distant retrograde orbits of the Earth-Moon system, halo L2 of the Sun-Earth system):
 * the Jacobi constant, the correction to periodicity (the catalogue state must be a fixed point
 * of the corrector), the period and the stability index; the Lagrange points of the catalogue;
 * the analytic Jacobian; Richardson's halo orbits as first guesses (error falling with the
 * amplitude); invariant manifolds (growth by the eigenvalue per period); the pulsating frame of
 * two real bodies. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/cr3bp.h"
#include "myb/integ.h"
#include "myb/io.h"
#include "myb/model.h"

static int read_line(FILE* f, char* buf, size_t size) {
  return fgets(buf, (int)size, f) != NULL;
}

static int parse(const char* s, double* v, int count, char* name, char* sym) {
  int used = 0;
  if (sscanf(s, "%31s %3s%n", name, sym, &used) != 2) {
    return 0;
  }
  const char* p = s + used;
  for (int k = 0; k < count; ++k) {
    while (*p == ' ') {
      ++p;
    }
    size_t n = 0;
    while (p[n] != ' ' && p[n] != '\n' && p[n] != '\r' && p[n] != '\0') {
      ++n;
    }
    size_t consumed = 0;
    if (myb_parse_double(p, n, &v[k], &consumed) != MYB_OK || consumed != n) {
      return 0;
    }
    p += n;
  }
  return 1;
}

static void test_catalogue(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/cr3bp/jpl_periodic_orbits.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  char line[1024];
  double worst_jacobi = 0.0;
  double worst_state = 0.0;
  double worst_period = 0.0;
  double worst_stability = 0.0;
  int rows = 0;
  while (read_line(f, line, sizeof line)) {
    if (line[0] == '#') {
      continue;
    }
    char name[32];
    char sym[4];
    double v[10];
    if (!parse(line, v, 10, name, sym)) {
      MT_CHECK(0);
      continue;
    }
    ++rows;
    const double mu = v[0];
    double y[6] = {v[1], v[2], v[3], v[4], v[5], v[6]};
    const double jacobi = v[7];
    const double period = v[8];
    const double stability = v[9];
    worst_jacobi = fmax(worst_jacobi, fabs(myb_cr3bp_jacobi(mu, y) - jacobi) / fabs(jacobi));
    /* symmetric form exactly (the catalogue gives 1e-13 residues) */
    const int xz = strcmp(sym, "xz") == 0;
    double y0[6] = {y[0], 0.0, xz ? y[2] : 0.0, 0.0, y[4], xz ? 0.0 : y[5]};
    const int planar = xz && fabs(y[2]) < 1e-12;
    myb_cr3bp_correction c = {.symmetry = xz ? MYB_CR3BP_SYMMETRY_XZ : MYB_CR3BP_SYMMETRY_X,
                              .fixed = xz ? (planar ? 0 : 2) : 0,
                              .tolerance = 1e-12};
    double t = 0.0;
    int ok = 0;
    for (int crossing = 1; crossing <= 4 && !ok; ++crossing) {
      double trial[6];
      memcpy(trial, y0, sizeof trial);
      c.crossing = crossing;
      if (myb_cr3bp_correct(mu, &c, trial, &t) == MYB_OK && fabs(t - period) < 1e-6 * period) {
        ok = 1;
        for (int i = 0; i < 6; ++i) {
          worst_state = fmax(worst_state, fabs(trial[i] - y0[i]));
        }
        worst_period = fmax(worst_period, fabs(t - period) / period);
        double mono[36];
        double re[6];
        double im[6];
        double index = 0.0;
        MT_CHECK(myb_cr3bp_monodromy(mu, trial, t, mono) == MYB_OK);
        MT_CHECK(myb_cr3bp_stability(mono, re, im, &index) == MYB_OK);
        worst_stability = fmax(worst_stability, fabs(index - stability) / stability);
        if (fabs(index - stability) > 1e-5 * stability) {
          MT_CHECK(strcmp(name, "em_lyap_l2") == 0);
          printf("  %s: stability %.10g vs %.10g\n", name, index, stability);
        }
      }
    }
    if (!ok) {
      printf("  %s (x0 %.6f): not corrected to the catalogue period %.6f\n", name, y0[0], period);
    }
    MT_CHECK(ok);
  }
  fclose(f);
  printf("JPL periodic orbits, %d of 8 families: Jacobi %.1e, corrected state %.1e, period %.1e, "
         "stability index %.1e\n",
         rows, worst_jacobi, worst_state, worst_period, worst_stability);
  MT_CHECK(rows == 56);
  MT_CHECK(worst_jacobi < 1e-12 && worst_state < 1e-9 && worst_period < 1e-10);
  /* measured 1.1e-4, from two L2 Lyapunov orbits passing 500 km above the Moon, where the
   * monodromy is most sensitive; ours moves by 1e-6 only when the integration is tightened
   * from 1e-11 to 1e-13 (and the step limit from 0.05 to 0.01); every other orbit agrees to
   * better than 1e-5 */
  MT_CHECK(worst_stability < 2e-4);
}

static void test_lagrange_and_jacobian(void) {
  const double mu = 1.215058560962404e-02;
  double l[5][3];
  MT_CHECK(myb_cr3bp_lagrange(mu, l) == MYB_OK);
  /* the catalogue's system block */
  MT_CHECK(fabs(l[0][0] - 0.836915125772357) < 1e-14);
  MT_CHECK(fabs(l[1][0] - 1.15568216544488) < 1e-14);
  MT_CHECK(fabs(l[2][0] - -1.00506264581028) < 1e-14);
  MT_CHECK(fabs(l[3][0] - 0.487849414390376) < 1e-15 && fabs(l[3][1] - 0.866025403784439) < 1e-15);
  for (int k = 0; k < 5; ++k) {
    double a[3];
    const double zero[3] = {0.0, 0.0, 0.0};
    myb_cr3bp_acceleration(mu, l[k], zero, a);
    MT_CHECK(fabs(a[0]) < 1e-14 && fabs(a[1]) < 1e-14);
  }
  double dummy[5][3];
  MT_CHECK(myb_cr3bp_lagrange(0.0, dummy) == MYB_E_DOMAIN);
  /* analytic Jacobian vs central differences of the force model */
  double m = mu;
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "cr3bp", myb_cr3bp_term, NULL, &m, NULL) == MYB_OK);
  const double x[3] = {0.83, 0.05, 0.07};
  const double v[3] = {0.01, 0.2, -0.03};
  double nx[9];
  double nv[9];
  double ax[9];
  double av[9];
  MT_CHECK(myb_force_model_jacobian(model, 0.0, x, v, nx, nv, NULL) == MYB_OK);
  myb_cr3bp_jacobian(mu, x, ax, av);
  double worst = 0.0;
  for (int i = 0; i < 9; ++i) {
    worst = fmax(worst, fmax(fabs(nx[i] - ax[i]), fabs(nv[i] - av[i])));
  }
  MT_CHECK(worst < 1e-7);
  myb_force_model_destroy(model);
}

static void test_richardson(void) {
  static const struct {
    double mu;
    int point;
    double az;
  } cases[] = {{1.215058560962404e-02, 1, 0.02}, {1.215058560962404e-02, 1, 0.01},
               {1.215058560962404e-02, 2, 0.02}, {1.215058560962404e-02, 2, 0.01},
               {3.0542e-06, 2, 8e-4},            {3.0542e-06, 2, 4e-4}};
  double err[6];
  for (int k = 0; k < 6; ++k) {
    double y0[6];
    double t = 0.0;
    MT_CHECK(myb_cr3bp_richardson(cases[k].mu, cases[k].point, cases[k].az, 1, y0, &t) ==
             MYB_OK);
    double y[6];
    memcpy(y, y0, sizeof y);
    const myb_cr3bp_correction c = {.symmetry = MYB_CR3BP_SYMMETRY_XZ, .fixed = 2};
    double period = 0.0;
    const myb_status st = myb_cr3bp_correct(cases[k].mu, &c, y, &period);
    MT_CHECK(st == MYB_OK);
    /* the error of the guess in units of the distance gamma of the point to the smaller
     * primary: halo orbits keep an in-plane amplitude of about 0.2 gamma however small Az is,
     * and that amplitude sets the third-order error */
    double l[5][3];
    MT_CHECK(myb_cr3bp_lagrange(cases[k].mu, l) == MYB_OK);
    const double gamma = fabs(l[cases[k].point - 1][0] - (1.0 - cases[k].mu));
    err[k] = fmax(fabs(y[0] - y0[0]), fabs(y[4] - y0[4])) / gamma;
    printf("Richardson halo L%d, Az %.4g: guess off by %.2e gamma; period %.6f (guess %.6f)\n",
           cases[k].point, cases[k].az, err[k], period, t);
    MT_CHECK(err[k] < 0.03);
    MT_CHECK(fabs(period - t) < 0.003 * t);
  }
  double y0[6];
  double t = 0.0;
  MT_CHECK(myb_cr3bp_richardson(0.0121, 3, 0.01, 1, y0, &t) == MYB_E_DOMAIN);
}

typedef struct cr3bp_ctx {
  double mu;
} cr3bp_ctx;

static myb_status rhs(void* context, double t, const double* x, const double* v, double* a) {
  (void)t;
  myb_cr3bp_acceleration(((const cr3bp_ctx*)context)->mu, x, v, a);
  return MYB_OK;
}

static double departure(double mu, const double y_start[6], const double orbit[6], double t) {
  cr3bp_ctx ctx = {mu};
  myb_ias15* in = NULL;
  const myb_ias15_options o = {.epsilon = 1e-12};
  MT_CHECK(myb_ias15_create(3, &o, NULL, &in) == MYB_OK);
  MT_CHECK(myb_ias15_start(in, rhs, &ctx, 0.0, y_start) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(in, t) == MYB_OK);
  const double* z = myb_ias15_state(in);
  const double d = sqrt((z[0] - orbit[0]) * (z[0] - orbit[0]) + (z[1] - orbit[1]) * (z[1] - orbit[1]) +
                        (z[2] - orbit[2]) * (z[2] - orbit[2]));
  myb_ias15_destroy(in);
  return d;
}

static void test_manifold(void) {
  /* an Earth-Moon L1 Lyapunov orbit of the JPL catalogue (stability index 1223.16) */
  const double mu = 1.215058560962404e-02;
  double y0[6] = {8.2522726142834313e-01, 0.0, 0.0, 0.0, 1.0748883813645686e-01, 0.0};
  double period = 0.0;
  const myb_cr3bp_correction c = {.symmetry = MYB_CR3BP_SYMMETRY_XZ, .fixed = 0};
  MT_CHECK(myb_cr3bp_correct(mu, &c, y0, &period) == MYB_OK);
  MT_CHECK(fabs(period - 2.7281356678175590) < 1e-10);
  double mono[36];
  double re[6];
  double im[6];
  double index = 0.0;
  MT_CHECK(myb_cr3bp_monodromy(mu, y0, period, mono) == MYB_OK);
  MT_CHECK(myb_cr3bp_stability(mono, re, im, &index) == MYB_OK);
  const double lambda = index + sqrt(index * index - 1.0);
  /* unstable manifold: forward one period the displacement grows by lambda; stable manifold:
   * backward one period likewise */
  const double eps = 1e-9;
  double yu[6];
  double ys[6];
  MT_CHECK(myb_cr3bp_manifold(mu, y0, period, 0.0, 1, 1, eps, yu) == MYB_OK);
  MT_CHECK(myb_cr3bp_manifold(mu, y0, period, 0.0, 0, 1, eps, ys) == MYB_OK);
  const double du = departure(mu, yu, y0, period);
  const double ds = departure(mu, ys, y0, -period);
  printf("L1 Lyapunov orbit (index %.2f, lambda %.1f): unstable manifold after one period "
         "%.3g, stable backwards %.3g (expected %.3g)\n",
         index, lambda, du, ds, eps * lambda);
  MT_CHECK(fabs(du / (eps * lambda) - 1.0) < 0.05 && fabs(ds / (eps * lambda) - 1.0) < 0.05);
  /* a point of the orbit itself departs by its own rounding only */
  MT_CHECK(departure(mu, y0, y0, period) < 1e-9);
  double yy[6];
  MT_CHECK(myb_cr3bp_manifold(mu, y0, period, 1.5, 1, 1, eps, yy) == MYB_E_DOMAIN);
}

static void test_frames(void) {
  /* Earth and Moon on a circular orbit: L4 is fixed in the rotating frame */
  const double gm1 = 3.986004418e14;
  const double gm2 = 4.9028e12;
  const double d = 3.844e8;
  const double n = sqrt((gm1 + gm2) / (d * d * d));
  const double mu = gm2 / (gm1 + gm2);
  const double th = 0.7;
  /* bodies about their barycentre, plus an arbitrary common velocity of the origin */
  const double vb[3] = {1000.0, -2000.0, 300.0};
  double b1[6];
  double b2[6];
  const double c = cos(th);
  const double s = sin(th);
  b1[0] = -mu * d * c;
  b1[1] = -mu * d * s;
  b1[2] = 0.0;
  b1[3] = mu * d * n * s + vb[0];
  b1[4] = -mu * d * n * c + vb[1];
  b1[5] = vb[2];
  b2[0] = (1.0 - mu) * d * c;
  b2[1] = (1.0 - mu) * d * s;
  b2[2] = 0.0;
  b2[3] = -(1.0 - mu) * d * n * s + vb[0];
  b2[4] = (1.0 - mu) * d * n * c + vb[1];
  b2[5] = vb[2];
  const double l4[6] = {0.5 - mu, 0.8660254037844386, 0.0, 0.0, 0.0, 0.0};
  double state[6];
  MT_CHECK(myb_cr3bp_to_ephemeris(b1, b2, gm1, gm2, l4, state) == MYB_OK);
  double back[6];
  MT_CHECK(myb_cr3bp_from_ephemeris(b1, b2, gm1, gm2, state, back) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 6; ++i) {
    worst = fmax(worst, fabs(back[i] - l4[i]));
  }
  MT_CHECK(worst < 1e-14);
  /* the inertial velocity of L4 is that of rigid rotation with the bodies */
  const double r4[3] = {state[0], state[1], state[2]};
  const double expect_v[3] = {-n * r4[1] + vb[0], n * r4[0] + vb[1], vb[2]};
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(state[3 + i] - expect_v[i]) < 1e-9);
  }
  /* a moving point: round trip */
  const double p[6] = {0.84, -0.01, 0.05, 0.02, -0.1, 0.03};
  MT_CHECK(myb_cr3bp_to_ephemeris(b1, b2, gm1, gm2, p, state) == MYB_OK);
  MT_CHECK(myb_cr3bp_from_ephemeris(b1, b2, gm1, gm2, state, back) == MYB_OK);
  worst = 0.0;
  for (int i = 0; i < 6; ++i) {
    worst = fmax(worst, fabs(back[i] - p[i]));
  }
  MT_CHECK(worst < 1e-14);
  printf("pulsating frame of two bodies: L4 fixed, round trips %.1e\n", worst);
}

/* How far the orbit through the first cut misses itself after the period, every component. */
static double closure(double mu, const myb_cr3bp_orbit* orbit) {
  cr3bp_ctx ctx = {mu};
  myb_ias15* in = NULL;
  const myb_ias15_options o = {.epsilon = 1e-12};
  MT_CHECK(myb_ias15_create(3, &o, NULL, &in) == MYB_OK);
  MT_CHECK(myb_ias15_start(in, rhs, &ctx, 0.0, orbit->states) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(in, orbit->period) == MYB_OK);
  const double* z = myb_ias15_state(in);
  double worst = 0.0;
  for (int i = 0; i < 6; ++i) {
    worst = fmax(worst, fabs(z[i] - orbit->states[i]));
  }
  myb_ias15_destroy(in);
  return worst;
}

static void test_shooting(void) {
  const double mu = 1.215058560962404e-02;
  /* the L1 Lyapunov orbit of the catalogue, cut into four */
  const double y0[6] = {8.2522726142834313e-01, 0.0, 0.0, 0.0, 1.0748883813645686e-01, 0.0};
  const double period = 2.7281356678175590;
  myb_cr3bp_orbit orbit;
  MT_CHECK(myb_cr3bp_orbit_from_state(mu, y0, period, 4, &orbit) == MYB_OK);
  MT_CHECK(orbit.patches == 4 && !orbit.has_tangent);
  MT_CHECK(myb_cr3bp_orbit_correct(mu, &orbit, 0.0, 0) == MYB_OK);
  double moved = 0.0;
  for (int i = 0; i < 6; ++i) {
    moved = fmax(moved, fabs(orbit.states[i] - y0[i]));
  }
  printf("multiple shooting, 4 cuts: the catalogue orbit is a fixed point to %.1e (period "
         "%.1e); it closes to %.1e\n",
         moved, fabs(orbit.period - period), closure(mu, &orbit));
  MT_CHECK(moved < 1e-9 && fabs(orbit.period - period) < 1e-9 && closure(mu, &orbit) < 1e-9);
  MT_CHECK(orbit.has_tangent);
  /* every cut pushed off, without any symmetry: back to a periodic orbit of the family */
  const double c0 = myb_cr3bp_jacobi(mu, orbit.states);
  myb_cr3bp_orbit pushed = orbit;
  pushed.has_tangent = 0;
  for (size_t k = 0; k < 24; ++k) {
    pushed.states[k] += 2.0e-5 * sin(1.0 + 2.7 * (double)k);
  }
  pushed.period *= 1.00002;
  const double open = closure(mu, &pushed);
  MT_CHECK(myb_cr3bp_orbit_correct(mu, &pushed, 0.0, 0) == MYB_OK);
  printf("  cuts pushed off by 2e-5: from %.1e to %.1e; the Jacobi constant %.9f (%.9f)\n",
         open, closure(mu, &pushed), myb_cr3bp_jacobi(mu, pushed.states), c0);
  MT_CHECK(open > 1e-4 && closure(mu, &pushed) < 1e-9);
  MT_CHECK(fabs(myb_cr3bp_jacobi(mu, pushed.states) - c0) < 1e-3);
  /* one cut: single shooting without symmetry */
  myb_cr3bp_orbit single;
  MT_CHECK(myb_cr3bp_orbit_from_state(mu, pushed.states, pushed.period * 1.00001, 1, &single) ==
           MYB_OK);
  single.states[0] += 1.0e-6;
  MT_CHECK(myb_cr3bp_orbit_correct(mu, &single, 0.0, 0) == MYB_OK);
  MT_CHECK(closure(mu, &single) < 1e-9);

  /* the family followed by its arclength: every member closes, and the corrector that
   * assumes the symmetry of the family, started from its first cut, finds the same orbit */
  myb_cr3bp_orbit member = orbit;
  double worst_closure = 0.0;
  double worst_symmetric = 0.0;
  double worst_period = 0.0;
  const myb_cr3bp_correction symmetric = {.symmetry = MYB_CR3BP_SYMMETRY_XZ, .fixed = 0};
  double last_period = member.period;
  int growing = 1;
  for (int step = 0; step < 25; ++step) {
    MT_CHECK(myb_cr3bp_orbit_continue(mu, &member, 0.01, 0.0, 0) == MYB_OK);
    worst_closure = fmax(worst_closure, closure(mu, &member));
    double ys[6] = {member.states[0], 0.0, 0.0, 0.0, member.states[4], 0.0};
    worst_symmetric = fmax(worst_symmetric, fabs(member.states[1]));
    worst_symmetric = fmax(worst_symmetric, fabs(member.states[3]));
    double ps = 0.0;
    MT_CHECK(myb_cr3bp_correct(mu, &symmetric, ys, &ps) == MYB_OK);
    worst_period = fmax(worst_period, fabs(ps - member.period));
    worst_symmetric = fmax(worst_symmetric, fabs(ys[4] - member.states[4]));
    growing = growing && member.period > last_period;
    last_period = member.period;
  }
  printf("  the family by its arclength, 25 steps of 0.01: x0 from %.6f to %.6f, period to "
         "%.6f; every member closes to %.1e; against the symmetric corrector: state %.1e, "
         "period %.1e\n",
         y0[0], member.states[0], member.period, worst_closure, worst_symmetric, worst_period);
  MT_CHECK(worst_closure < 1e-8);
  MT_CHECK(worst_symmetric < 1e-9 && worst_period < 1e-9);
  MT_CHECK(growing && member.period > period + 0.05);
  /* the other way from the same orbit: the period falls */
  member = orbit;
  MT_CHECK(myb_cr3bp_orbit_continue(mu, &member, -0.01, 0.0, 0) == MYB_OK);
  MT_CHECK(member.period < period && closure(mu, &member) < 1e-9);
  /* refused */
  myb_cr3bp_orbit raw;
  MT_CHECK(myb_cr3bp_orbit_from_state(mu, y0, period, 4, &raw) == MYB_OK);
  MT_CHECK(myb_cr3bp_orbit_continue(mu, &raw, 0.01, 0.0, 0) == MYB_E_DOMAIN); /* no tangent */
  MT_CHECK(myb_cr3bp_orbit_from_state(mu, y0, period, 11, &raw) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_orbit_from_state(mu, y0, -1.0, 4, &raw) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_orbit_continue(mu, &orbit, 0.0, 0.0, 0) == MYB_E_DOMAIN);
}

/* ---- Quasi-periodic orbits --------------------------------------------------------------------- */

static const double PI = 3.141592653589793;

static myb_status flow_rhs(void* context, double t, const double* y, double* dydt) {
  (void)t;
  dydt[0] = y[3];
  dydt[1] = y[4];
  dydt[2] = y[5];
  myb_cr3bp_acceleration(*(const double*)context, y, y + 3, dydt + 3);
  return MYB_OK;
}

/* The point of the circle at the angle 0 carried by the flow, by an integrator that the
 * solver does not use: after k T it must be on the circle at the angle k rho. The largest
 * distance over `turns` returns. */
static double returns_miss(double mu, const myb_cr3bp_torus* torus, int turns) {
  const myb_dop853_options options = {.rtol = 1e-13, .atol = 1e-13};
  myb_dop853* integrator = NULL;
  MT_CHECK(myb_dop853_create(6, &options, NULL, &integrator) == MYB_OK);
  double m = mu;
  MT_CHECK(myb_dop853_start(integrator, flow_rhs, &m, 0.0, torus->states) == MYB_OK);
  double worst = 0.0;
  for (int k = 1; k <= turns; ++k) {
    MT_CHECK(myb_dop853_integrate(integrator, k * torus->period) == MYB_OK);
    const double* y = myb_dop853_state(integrator);
    double on_circle[6];
    myb_cr3bp_torus_state(torus, k * torus->rotation, on_circle);
    for (int i = 0; i < 6; ++i) {
      worst = fmax(worst, fabs(y[i] - on_circle[i]));
    }
  }
  myb_dop853_destroy(integrator);
  return worst;
}

/* The spread of the Jacobi constant over the points, and the root mean square distance of
 * the circle from a state. */
static double jacobi_spread(double mu, const myb_cr3bp_torus* torus, double* mean) {
  double lo = 1e300;
  double hi = -1e300;
  double sum = 0.0;
  for (size_t i = 0; i < torus->points; ++i) {
    const double c = myb_cr3bp_jacobi(mu, torus->states + 6 * i);
    lo = fmin(lo, c);
    hi = fmax(hi, c);
    sum += c;
  }
  *mean = sum / (double)torus->points;
  return hi - lo;
}

static double circle_size(const myb_cr3bp_torus* torus, const double y0[6]) {
  double sum = 0.0;
  for (size_t i = 0; i < torus->points; ++i) {
    for (int k = 0; k < 6; ++k) {
      const double d = torus->states[6 * i + (size_t)k] - y0[k];
      sum += d * d;
    }
  }
  return sqrt(sum / (double)torus->points);
}

static void test_torus(void) {
  const double mu = 1.215058560962404e-02;
  const myb_cr3bp_correction c = {.symmetry = MYB_CR3BP_SYMMETRY_XZ, .fixed = 0};
  /* a distant retrograde orbit of the catalogue: stable, so the motion about it can be
   * followed for many turns */
  double dro[6] = {8.2915209119468902e-01, 0.0, 0.0, 0.0, 4.9478973105640744e-01, 0.0};
  double period = 0.0;
  MT_CHECK(myb_cr3bp_correct(mu, &c, dro, &period) == MYB_OK);
  static myb_cr3bp_torus torus;
  static myb_cr3bp_torus keep;
  MT_CHECK(myb_cr3bp_torus_from_orbit(mu, dro, period, 21, 0, 1e-3, &torus) == MYB_OK);
  const double rho0 = torus.rotation;
  const double linear_miss = returns_miss(mu, &torus, 1);
  MT_CHECK(myb_cr3bp_torus_continue(mu, &torus, MYB_CR3BP_TORUS_JACOBI, 0.0, 0.0, 0) == MYB_OK);
  double mean = 0.0;
  const double spread0 = jacobi_spread(mu, &torus, &mean);
  const double miss0 = returns_miss(mu, &torus, 20);
  printf("quasi-periodic orbits\n  about a distant retrograde orbit (T = %.6f), 21 points: "
         "rho = %.9f (%.9f of the monodromy matrix); the linear circle of size 1e-3 returns "
         "to %.1e, the corrected one over 20 turns to %.1e; Jacobi constant over the circle "
         "to %.1e; multipliers left %.1e\n",
         period, torus.rotation, rho0, linear_miss, miss0, spread0, torus.defect);
  MT_CHECK(linear_miss < 1e-4 && linear_miss > 1e-9); /* of the order of size^2 */
  MT_CHECK(miss0 < 1e-9);
  MT_CHECK(spread0 < 1e-10);
  MT_CHECK(torus.defect < 1e-9);
  MT_CHECK(fabs(torus.rotation - rho0) < 1e-4 && fabs(torus.period - period) < 1e-4);
  MT_CHECK(fabs(circle_size(&torus, dro) - 1e-3) < 1e-5);
  /* outwards along the family of one Jacobi constant */
  double mean_before = mean;
  for (int k = 0; k < 4; ++k) {
    MT_CHECK(myb_cr3bp_torus_continue(mu, &torus, MYB_CR3BP_TORUS_JACOBI, 5e-3, 0.0, 0) ==
             MYB_OK);
  }
  const double spread = jacobi_spread(mu, &torus, &mean);
  const double miss = returns_miss(mu, &torus, 20);
  printf("  four steps of 5e-3 along the family of one Jacobi constant: size %.5f, T = %.6f, "
         "rho = %.6f; over 20 turns to %.1e; the constant kept to %.1e, over the circle to "
         "%.1e; multipliers left %.1e\n",
         circle_size(&torus, dro), torus.period, torus.rotation, miss, fabs(mean - mean_before),
         spread, torus.defect);
  MT_CHECK(circle_size(&torus, dro) > 0.015);
  MT_CHECK(miss < 1e-8);
  MT_CHECK(fabs(mean - mean_before) < 1e-10 && spread < 1e-9);
  MT_CHECK(torus.defect < 1e-8);
  /* the two other families: T kept, rho kept; a step back returns */
  keep = torus;
  MT_CHECK(myb_cr3bp_torus_continue(mu, &torus, MYB_CR3BP_TORUS_PERIOD, 3e-3, 0.0, 0) == MYB_OK);
  MT_CHECK(fabs(torus.period - keep.period) < 1e-11 && fabs(torus.rotation - keep.rotation) > 1e-7);
  MT_CHECK(returns_miss(mu, &torus, 10) < 1e-8);
  /* the tangent kept is that of the family of the last step: a step of 0 turns it to the
   * new family */
  torus = keep;
  MT_CHECK(myb_cr3bp_torus_continue(mu, &torus, MYB_CR3BP_TORUS_ROTATION, 0.0, 0.0, 0) ==
           MYB_OK);
  static myb_cr3bp_torus base;
  base = torus;
  MT_CHECK(myb_cr3bp_torus_continue(mu, &torus, MYB_CR3BP_TORUS_ROTATION, 3e-3, 0.0, 0) ==
           MYB_OK);
  MT_CHECK(fabs(torus.rotation - keep.rotation) < 1e-11 && fabs(torus.period - keep.period) > 1e-7);
  MT_CHECK(returns_miss(mu, &torus, 10) < 1e-8);
  MT_CHECK(myb_cr3bp_torus_continue(mu, &torus, MYB_CR3BP_TORUS_ROTATION, -3e-3, 0.0, 0) ==
           MYB_OK);
  double back = 0.0;
  for (size_t k = 0; k < 6 * torus.points; ++k) {
    back = fmax(back, fabs(torus.states[k] - base.states[k]));
  }
  printf("  a step along the family of one rho and back: the circle returns to %.1e\n", back);
  MT_CHECK(back < 1e-4); /* of the order of the step squared */
  /* the other pair of the orbit: out of its plane */
  MT_CHECK(myb_cr3bp_torus_from_orbit(mu, dro, period, 21, 1, 2e-3, &torus) == MYB_OK);
  MT_CHECK(myb_cr3bp_torus_continue(mu, &torus, MYB_CR3BP_TORUS_JACOBI, 0.0, 0.0, 0) == MYB_OK);
  double flat = 0.0;
  double steep = 0.0;
  for (size_t i = 0; i < torus.points; ++i) {
    flat = fmax(flat, fabs(keep.states[6 * i + 2]));
    steep = fmax(steep, fabs(torus.states[6 * i + 2]));
  }
  printf("  the two pairs of the orbit: rho = %.6f, out of its plane (z up to %.1e), and "
         "%.6f, in it (z up to %.1e); the second over 20 turns to %.1e\n",
         keep.rotation, flat, torus.rotation, steep, returns_miss(mu, &torus, 20));
  MT_CHECK((flat < 1e-12) != (steep < 1e-12)); /* one of the two is planar */
  MT_CHECK(returns_miss(mu, &torus, 20) < 1e-9);
  MT_CHECK(myb_cr3bp_torus_from_orbit(mu, dro, period, 21, 2, 1e-3, &torus) == MYB_E_DOMAIN);

  /* a Lissajous orbit: about a small planar Lyapunov orbit of L1 the turn across it is that
   * of the linear motion out of the plane, 2 pi (1 - w_z / w_xy) */
  double points[5][3];
  MT_CHECK(myb_cr3bp_lagrange(mu, points) == MYB_OK);
  const double gamma = 1.0 - mu - points[0][0];
  const double c2 = (mu + (1.0 - mu) * gamma * gamma * gamma /
                              ((1.0 - gamma) * (1.0 - gamma) * (1.0 - gamma))) /
                    (gamma * gamma * gamma);
  const double in_plane = sqrt(0.5 * (2.0 - c2 + sqrt(9.0 * c2 * c2 - 8.0 * c2)));
  const double out_of_plane = sqrt(c2);
  const double kappa = (in_plane * in_plane + 1.0 + 2.0 * c2) / (2.0 * in_plane);
  const double amplitude = 1e-4;
  double small[6] = {points[0][0] - amplitude, 0.0, 0.0, 0.0, kappa * in_plane * amplitude, 0.0};
  double small_period = 0.0;
  MT_CHECK(myb_cr3bp_correct(mu, &c, small, &small_period) == MYB_OK);
  MT_CHECK(myb_cr3bp_torus_from_orbit(mu, small, small_period, 15, 0, 1e-5, &torus) == MYB_OK);
  const double expected = 2.0 * PI * (1.0 - out_of_plane / in_plane);
  printf("  about a Lyapunov orbit of L1 of 38 km: rho = %.7f against %.7f of the linear "
         "motion (frequencies %.6f in the plane, %.6f out of it)\n",
         torus.rotation, expected, in_plane, out_of_plane);
  MT_CHECK(fabs(small_period - 2.0 * PI / in_plane) < 1e-5);
  MT_CHECK(fabs(torus.rotation - expected) < 1e-5);
  /* about the Lyapunov orbit of the catalogue (4500 km, unstable 2400-fold in a turn): a
   * Lissajous orbit reaching 3000 km out of the plane */
  double lyapunov[6] = {8.2522726142834313e-01, 0.0, 0.0, 0.0, 1.0748883813645686e-01, 0.0};
  double lyapunov_period = 0.0;
  MT_CHECK(myb_cr3bp_correct(mu, &c, lyapunov, &lyapunov_period) == MYB_OK);
  MT_CHECK(myb_cr3bp_torus_from_orbit(mu, lyapunov, lyapunov_period, 25, 0, 1e-3, &torus) ==
           MYB_OK);
  MT_CHECK(myb_cr3bp_torus_continue(mu, &torus, MYB_CR3BP_TORUS_JACOBI, 0.0, 0.0, 0) == MYB_OK);
  for (int k = 0; k < 3; ++k) {
    MT_CHECK(myb_cr3bp_torus_continue(mu, &torus, MYB_CR3BP_TORUS_JACOBI, 4e-3, 0.0, 0) ==
             MYB_OK);
  }
  double height = 0.0;
  for (size_t i = 0; i < torus.points; ++i) {
    height = fmax(height, fabs(torus.states[6 * i + 2]));
  }
  const double lissajous_miss = returns_miss(mu, &torus, 1);
  const double lissajous_spread = jacobi_spread(mu, &torus, &mean);
  printf("  a Lissajous orbit about the Lyapunov orbit of 4500 km, 25 points: %.0f km out of "
         "the plane, T = %.6f, rho = %.6f; after one turn on the circle to %.1e; Jacobi "
         "constant over the circle to %.1e; multipliers left %.1e\n",
         height * 384400.0, torus.period, torus.rotation, lissajous_miss, lissajous_spread,
         torus.defect);
  MT_CHECK(height * 384400.0 > 1500.0);
  MT_CHECK(lissajous_miss < 1e-6);
  MT_CHECK(lissajous_spread < 1e-9);
  /* the circle between its points */
  double at[6];
  myb_cr3bp_torus_state(&torus, 2.0 * PI * 3.0 / 25.0, at);
  for (int k = 0; k < 6; ++k) {
    MT_CHECK(fabs(at[k] - torus.states[18 + k]) < 1e-14);
  }
  /* refused */
  MT_CHECK(myb_cr3bp_torus_from_orbit(mu, dro, period, 20, 0, 1e-3, &keep) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_torus_from_orbit(mu, dro, period, 21, 0, 0.0, &keep) == MYB_E_DOMAIN);
  keep = torus;
  keep.has_tangent = 0;
  MT_CHECK(myb_cr3bp_torus_continue(mu, &keep, MYB_CR3BP_TORUS_JACOBI, 0.0, 0.0, 0) ==
           MYB_E_DOMAIN);
  /* a step too long for the method is refused and leaves the torus as it was */
  keep = torus;
  MT_CHECK(myb_cr3bp_torus_continue(mu, &keep, MYB_CR3BP_TORUS_JACOBI, 0.5, 0.0, 0) ==
           MYB_E_NOCONV);
  MT_CHECK(memcmp(&keep, &torus, sizeof keep) == 0);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_lagrange_and_jacobian();
  test_frames();
  test_catalogue();
  test_richardson();
  test_manifold();
  test_shooting();
  test_torus();
  return mt_finish("test_cr3bp");
}
