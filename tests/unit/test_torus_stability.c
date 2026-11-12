/* SPDX-License-Identifier: BSD-3-Clause */
/* The stability of the invariant tori of the circular restricted three-body problem and their
 * manifolds (myb/cr3bp.h).
 *
 * The spectrum of the linearised map of the invariant circle against the monodromy matrix of
 * the periodic orbit that a small torus surrounds: about the L1 Lyapunov orbit of the JPL
 * catalogue (Lissajous tori) and two L2 halo orbits (quasi-halo tori) the unstable eigenvalue
 * of the torus is that of the orbit, the difference falling as the square of the size of the
 * torus; about a distant retrograde orbit every resolved eigenvalue is on the unit circle,
 * and the spurious ones of a larger torus are not taken for instability. The eigenvectors
 * against the map itself, by another integrator (DOP853) and finite differences. The
 * manifolds: their distance from the torus grows by the eigenvalue in a turn; carried to the
 * section of the Moon, x = 1 - mu, the manifolds of a Lissajous torus cross it where those of
 * its Lyapunov orbit do, out of the plane by the size of the torus and in it by its square,
 * the stable ones the mirror images of the unstable ones (the symmetry of the problem with
 * time reversed), the Jacobi constant kept. Propagation to a time and to a section against
 * DOP853 and the period of the catalogue orbit. */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "mtest.h"
#include "myb/cr3bp.h"
#include "myb/integ.h"

static const double MU = 1.215058560962404e-02;
static const double PI = 3.141592653589793;

/* the L1 Lyapunov orbit of the JPL catalogue (4500 km, stability index 1223.16) */
static const double LYAPUNOV[6] = {8.2522726142834313e-01, 0.0, 0.0, 0.0,
                                   1.0748883813645686e-01, 0.0};

static myb_status flow_rhs(void* context, double t, const double* y, double* dydt) {
  (void)t;
  dydt[0] = y[3];
  dydt[1] = y[4];
  dydt[2] = y[5];
  myb_cr3bp_acceleration(*(const double*)context, y, y + 3, dydt + 3);
  return MYB_OK;
}

/* The flow over t by DOP853, an integrator that the functions under test do not use. */
static void flow(double t, const double y0[6], double y[6]) {
  const myb_dop853_options options = {.rtol = 1e-13, .atol = 1e-13};
  myb_dop853* integrator = NULL;
  double mu = MU;
  MT_CHECK(myb_dop853_create(6, &options, NULL, &integrator) == MYB_OK);
  MT_CHECK(myb_dop853_start(integrator, flow_rhs, &mu, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_dop853_integrate(integrator, t) == MYB_OK);
  memcpy(y, myb_dop853_state(integrator), 6 * sizeof *y);
  myb_dop853_destroy(integrator);
}

/* The periodic orbit corrected (symmetric to the xz-plane, `fixed` held), its period, the
 * stability index and the largest eigenvalue of its monodromy matrix. */
static double orbit_lambda(double y0[6], int fixed, double* period, double* index) {
  const myb_cr3bp_correction c = {.symmetry = MYB_CR3BP_SYMMETRY_XZ, .fixed = fixed};
  MT_CHECK(myb_cr3bp_correct(MU, &c, y0, period) == MYB_OK);
  double mono[36];
  double re[6];
  double im[6];
  MT_CHECK(myb_cr3bp_monodromy(MU, y0, *period, mono) == MYB_OK);
  MT_CHECK(myb_cr3bp_stability(mono, re, im, index) == MYB_OK);
  return *index + sqrt(*index * *index - 1.0);
}

/* The corrected torus about the orbit (the first pair on the unit circle). */
static void torus_about(const double y0[6], double period, size_t points, double size,
                        myb_cr3bp_torus* torus) {
  MT_CHECK(myb_cr3bp_torus_from_orbit(MU, y0, period, points, 0, size, torus) == MYB_OK);
  MT_CHECK(myb_cr3bp_torus_continue(MU, torus, MYB_CR3BP_TORUS_JACOBI, 0.0, 0.0, 0) == MYB_OK);
}

static void print_spectrum(const char* name, const myb_cr3bp_torus* torus,
                           const myb_cr3bp_torus_spectrum* s) {
  printf(
      "  %s: %zu points, T = %.9f, rho = %.9f; %zu eigenvalues: %zu trivial, %zu elliptic, "
      "%zu hyperbolic, %zu unresolved; unstable %.9g, stable %.9g (product - 1 = %.1e), "
      "index %.6f; residual %.1e, tail %.1e\n",
      name, torus->points, torus->period, torus->rotation, s->count, s->trivial, s->elliptic,
      s->hyperbolic, s->unresolved, s->unstable, s->stable, s->unstable * s->stable - 1.0, s->index,
      s->residual, s->tail);
}

/* The Lissajous tori of 25 points about the Lyapunov orbit, of sizes 1e-3, 2e-3, 4e-3, and
 * the spectra of the first two: built once, used by the tests after the first. */
static double g_period;
static myb_cr3bp_torus g_tori[3];
static myb_cr3bp_torus_spectrum g_spectra[2];
static myb_cr3bp_torus g_scratch;
static myb_cr3bp_torus_spectrum g_spectrum;

/* ---- Propagation --------------------------------------------------------------------------- */

static void test_propagate(void) {
  double y0[6];
  memcpy(y0, LYAPUNOV, sizeof y0);
  double period = 0.0;
  double index = 0.0;
  (void)orbit_lambda(y0, 0, &period, &index);
  double y[6];
  double z[6];
  double t = 0.0;
  int crossed = 1;
  MT_CHECK(myb_cr3bp_propagate(MU, y0, 0.7 * period, NULL, y, &t, &crossed) == MYB_OK);
  MT_CHECK(t == 0.7 * period && crossed == 0);
  flow(0.7 * period, y0, z);
  double worst = 0.0;
  for (int i = 0; i < 6; ++i) {
    worst = fmax(worst, fabs(y[i] - z[i]));
  }
  /* the orbit starts on y = 0 with vy > 0 (not counted): it crosses y = 0 falling at T / 2,
   * rising at T; backward, the falling crossing is met at -T / 2 */
  const myb_cr3bp_section falling = {.axis = 1, .value = 0.0, .direction = -1};
  const myb_cr3bp_section rising = {.axis = 1, .value = 0.0, .direction = 1};
  double half = 0.0;
  double full = 0.0;
  double back = 0.0;
  double yh[6];
  double yf[6];
  double yb[6];
  int ch = 0;
  int cf = 0;
  int cb = 0;
  MT_CHECK(myb_cr3bp_propagate(MU, y0, 2.0 * period, &falling, yh, &half, &ch) == MYB_OK);
  MT_CHECK(myb_cr3bp_propagate(MU, y0, 2.0 * period, &rising, yf, &full, &cf) == MYB_OK);
  MT_CHECK(myb_cr3bp_propagate(MU, y0, -2.0 * period, &falling, yb, &back, &cb) == MYB_OK);
  double returned = 0.0;
  for (int i = 0; i < 6; ++i) {
    returned = fmax(returned, fabs(yf[i] - y0[i]));
  }
  printf(
      "propagation: against DOP853 over 0.7 T %.1e; y = 0 crossed falling at T / 2 %+.1e, "
      "rising at T %+.1e (the state returns to %.1e), falling backward at -T / 2 %+.1e\n",
      worst, half - 0.5 * period, full - period, returned, back + 0.5 * period);
  MT_CHECK(worst < 1e-11);
  MT_CHECK(ch && cf && cb);
  MT_CHECK(fabs(half - 0.5 * period) < 1e-10 && fabs(yh[1]) < 1e-14);
  MT_CHECK(fabs(full - period) < 1e-10 && returned < 1e-10);
  MT_CHECK(fabs(back + 0.5 * period) < 1e-10);
  /* in place */
  memcpy(y, y0, sizeof y);
  MT_CHECK(myb_cr3bp_propagate(MU, y, 0.7 * period, NULL, y, NULL, NULL) == MYB_OK);
  for (int i = 0; i < 6; ++i) {
    MT_CHECK(fabs(y[i] - z[i]) < 1e-11);
  }
  /* refused */
  const myb_cr3bp_section axis = {.axis = 6, .value = 0.0, .direction = 0};
  const myb_cr3bp_section direction = {.axis = 0, .value = 0.0, .direction = 2};
  MT_CHECK(myb_cr3bp_propagate(MU, y0, 1.0, &axis, y, NULL, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_propagate(MU, y0, 1.0, &direction, y, NULL, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_propagate(0.7, y0, 1.0, NULL, y, NULL, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_propagate(MU, y0, HUGE_VAL, NULL, y, NULL, NULL) == MYB_E_DOMAIN);
}

/* ---- The spectrum against the generating orbit ---------------------------------------------- */

static void test_lissajous(void) {
  double y0[6];
  memcpy(y0, LYAPUNOV, sizeof y0);
  double index = 0.0;
  const double lambda = orbit_lambda(y0, 0, &g_period, &index);
  printf("Lissajous tori about the L1 Lyapunov orbit (T = %.9f, lambda = %.6f, index %.6f)\n",
         g_period, lambda, index);
  const double sizes[3] = {1e-3, 2e-3, 4e-3};
  double differences[3];
  for (int k = 0; k < 3; ++k) {
    myb_cr3bp_torus* torus = &g_tori[k];
    myb_cr3bp_torus_spectrum* s = k < 2 ? &g_spectra[k] : &g_spectrum;
    torus_about(y0, g_period, 25, sizes[k], torus);
    MT_CHECK(myb_cr3bp_torus_stability(MU, torus, 0.0, s) == MYB_OK);
    print_spectrum("torus", torus, s);
    differences[k] = s->unstable / lambda - 1.0;
    printf(
        "    size %.0e: the unstable eigenvalue against that of the orbit %+.3e, the index "
        "%+.3e\n",
        sizes[k], differences[k], s->index / index - 1.0);
    /* every member of the two hyperbolic circles resolved, none elliptic */
    MT_CHECK(s->count == 150 && s->hyperbolic == 50 && s->elliptic == 0);
    MT_CHECK(s->trivial + s->unresolved == 100 && s->trivial >= 84);
    MT_CHECK(s->unstable > 1.0 && s->stable > 0.0 && s->stable < 1.0);
    MT_CHECK(fabs(s->unstable * s->stable - 1.0) < 1e-8);
    MT_CHECK(fabs(s->index / index - 1.0) < 2e-3);
    MT_CHECK(s->residual < 1e-7 && s->tail < 1e-12);
  }
  printf("    the differences grow as the size squared: %.3f, %.3f\n",
         differences[1] / differences[0], differences[2] / differences[1]);
  MT_CHECK(fabs(differences[0]) < 1e-4);
  MT_CHECK(fabs(differences[1] / differences[0] - 4.0) < 0.1);
  MT_CHECK(fabs(differences[2] / differences[1] - 4.0) < 0.1);
  /* the largest circle, 35 points (a matrix of 210): the same eigenvalue */
  torus_about(y0, g_period, MYB_CR3BP_MAX_CIRCLE, 1e-3, &g_scratch);
  const clock_t start = clock();
  MT_CHECK(myb_cr3bp_torus_stability(MU, &g_scratch, 0.0, &g_spectrum) == MYB_OK);
  const double seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
  print_spectrum("size 1e-3", &g_scratch, &g_spectrum);
  printf("    against 25 points %+.1e; %.2f s\n", g_spectrum.unstable / g_spectra[0].unstable - 1.0,
         seconds);
  MT_CHECK(g_spectrum.count == 210 && g_spectrum.hyperbolic == 70 && g_spectrum.elliptic == 0);
  MT_CHECK(fabs(g_spectrum.unstable / g_spectra[0].unstable - 1.0) < 1e-6);
}

static void test_halo(void) {
  /* two southern L2 halo orbits of the catalogue: indices 12.17 and 74.86 */
  const double halos[2][3] = {
      {1.1104372933451314e+00, -1.9396147171787578e-01, -2.2162607052127681e-01},
      {1.1430328046132601e+00, -1.5867031222543088e-01, -2.2220376826939378e-01}};
  for (int h = 0; h < 2; ++h) {
    double halo[6] = {halos[h][0], 0.0, halos[h][1], 0.0, halos[h][2], 0.0};
    double period = 0.0;
    double index = 0.0;
    const double lambda = orbit_lambda(halo, 2, &period, &index);
    torus_about(halo, period, 25, 1e-3, &g_scratch);
    MT_CHECK(myb_cr3bp_torus_stability(MU, &g_scratch, 0.0, &g_spectrum) == MYB_OK);
    printf("a quasi-halo torus of size 1e-3 about an L2 halo orbit (lambda %.6f, index %.6f)\n",
           lambda, index);
    print_spectrum("torus", &g_scratch, &g_spectrum);
    printf("    the unstable eigenvalue against that of the orbit %+.2e\n",
           g_spectrum.unstable / lambda - 1.0);
    MT_CHECK(g_spectrum.hyperbolic == 50 && g_spectrum.elliptic == 0);
    MT_CHECK(fabs(g_spectrum.unstable / lambda - 1.0) < 2e-5);
    MT_CHECK(fabs(g_spectrum.unstable * g_spectrum.stable - 1.0) < 1e-10);
  }
}

/* ---- A stable torus ------------------------------------------------------------------------ */

static void test_stable(void) {
  double dro[6] = {8.2915209119468902e-01, 0.0, 0.0, 0.0, 4.9478973105640744e-01, 0.0};
  double period = 0.0;
  double index = 0.0;
  (void)orbit_lambda(dro, 0, &period, &index);
  printf("a torus about a distant retrograde orbit (index %.12f)\n", index);
  torus_about(dro, period, 21, 1e-3, &g_scratch);
  MT_CHECK(myb_cr3bp_torus_stability(MU, &g_scratch, 0.0, &g_spectrum) == MYB_OK);
  print_spectrum("size 1e-3", &g_scratch, &g_spectrum);
  /* the normal pair exp(+-i sigma) on the unit circle: two elliptic circles */
  MT_CHECK(g_spectrum.hyperbolic == 0 && g_spectrum.elliptic == 42);
  MT_CHECK(g_spectrum.trivial + g_spectrum.unresolved == 84);
  MT_CHECK(g_spectrum.unstable == 0.0 && g_spectrum.stable == 0.0 && g_spectrum.index == 1.0);
  double y[6];
  MT_CHECK(myb_cr3bp_torus_manifold(&g_scratch, &g_spectrum, 0.0, 1, 1, 1e-6, y) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_torus_manifold(&g_scratch, &g_spectrum, 0.0, 0, 1, 1e-6, y) == MYB_E_DOMAIN);
  /* a larger one: some eigenvalues of high frequency leave the unit circle (by 5 %), but no
   * circle does */
  for (int k = 0; k < 4; ++k) {
    MT_CHECK(myb_cr3bp_torus_continue(MU, &g_scratch, MYB_CR3BP_TORUS_JACOBI, 5e-3, 0.0, 0) ==
             MYB_OK);
  }
  MT_CHECK(myb_cr3bp_torus_stability(MU, &g_scratch, 0.0, &g_spectrum) == MYB_OK);
  double off = 0.0;
  for (size_t k = 0; k < g_spectrum.count; ++k) {
    off = fmax(off, fabs(hypot(g_spectrum.re[k], g_spectrum.im[k]) - 1.0));
  }
  print_spectrum("size 0.02", &g_scratch, &g_spectrum);
  printf("    the farthest eigenvalue from the unit circle %.1e, unresolved\n", off);
  MT_CHECK(off > 1e-2);
  MT_CHECK(g_spectrum.hyperbolic == 0 && g_spectrum.elliptic > 0 && g_spectrum.index == 1.0);
}

/* ---- The eigenvectors against the nonlinear map ------------------------------------------- */

/* The map F(u) = R(-rho) phi_T(u) of the circle (or its inverse F^-1(u) = phi_-T(R(rho) u)),
 * the flow by DOP853, at the points of the circle u + eps v. */
static void map_of(const myb_cr3bp_torus* torus, const double* v, double eps, int inverse,
                   double* image) {
  const size_t n = torus->points;
  g_scratch = *torus;
  for (size_t k = 0; k < 6 * n; ++k) {
    g_scratch.states[k] = torus->states[k] + eps * v[k];
  }
  if (!inverse) {
    for (size_t j = 0; j < n; ++j) {
      flow(torus->period, g_scratch.states + 6 * j, image + 6 * j);
    }
    memcpy(g_scratch.states, image, 6 * n * sizeof *image);
    for (size_t i = 0; i < n; ++i) {
      myb_cr3bp_torus_state(&g_scratch, 2.0 * PI * (double)i / (double)n - torus->rotation,
                            image + 6 * i);
    }
  } else {
    double turned[6];
    for (size_t i = 0; i < n; ++i) {
      myb_cr3bp_torus_state(&g_scratch, 2.0 * PI * (double)i / (double)n + torus->rotation, turned);
      flow(-torus->period, turned, image + 6 * i);
    }
  }
}

/* (F(u + eps v) - F(u)) / eps = A v + O(eps) against lambda v (for the stable direction F^-1
 * and v / lambda): the largest difference relative to the largest component. */
static double map_miss(const myb_cr3bp_torus* torus, const double* v, double lambda, int inverse,
                       double eps) {
  static double base[6 * MYB_CR3BP_MAX_CIRCLE];
  static double moved[6 * MYB_CR3BP_MAX_CIRCLE];
  map_of(torus, v, 0.0, inverse, base);
  map_of(torus, v, eps, inverse, moved);
  const double expected = inverse ? 1.0 / lambda : lambda;
  double worst = 0.0;
  double size = 0.0;
  for (size_t k = 0; k < 6 * torus->points; ++k) {
    worst = fmax(worst, fabs((moved[k] - base[k]) / eps - expected * v[k]));
    size = fmax(size, fabs(expected * v[k]));
  }
  return worst / size;
}

static void test_directions(void) {
  const myb_cr3bp_torus* torus = &g_tori[0];
  const myb_cr3bp_torus_spectrum* s = &g_spectra[0];
  const double eps[2] = {1e-8, 1e-9};
  double unstable[2];
  double stable[2];
  for (int k = 0; k < 2; ++k) {
    unstable[k] = map_miss(torus, s->unstable_direction, s->unstable, 0, eps[k]);
    stable[k] = map_miss(torus, s->stable_direction, s->stable, 1, eps[k]);
  }
  printf(
      "the eigenvectors against the map by DOP853, (F(u + eps v) - F(u)) / eps - lambda v: "
      "eps 1e-8: unstable %.1e, stable %.1e; eps 1e-9: %.1e, %.1e\n",
      unstable[0], stable[0], unstable[1], stable[1]);
  MT_CHECK(unstable[0] < 1e-4 && stable[0] < 1e-4);
  MT_CHECK(unstable[1] < 2e-5 && stable[1] < 2e-5);
  /* the directions are unit in the mean square over the points, and turned to growing x */
  double norm = 0.0;
  double along = 0.0;
  for (size_t j = 0; j < torus->points; ++j) {
    along += s->unstable_direction[6 * j];
    for (int c = 0; c < 6; ++c) {
      norm += s->unstable_direction[6 * j + (size_t)c] * s->unstable_direction[6 * j + (size_t)c];
    }
  }
  MT_CHECK(fabs(norm / (double)torus->points - 1.0) < 1e-12 && along > 0.0);
}

/* ---- Manifolds ----------------------------------------------------------------------------- */

static double position_norm(const double y[6]) {
  return sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
}

/* A point of the manifold carried k turns (forward on the unstable one, backward on the
 * stable one) against the point of the circle it left carried as far, the circle's own error
 * (1e-11, grown by lambda^k) cancelling: their distance must be eps lambda^k |psi(theta +-
 * k rho)| / |psi(theta)| in position, psi the eigenfunction between the points by the
 * interpolant of the circle. The largest relative difference over seven angles; *rate: the
 * largest growth per turn (distance / eps)^(1/k). */
static double growth_miss(int unstable, int turns, double eps, double* rate) {
  const myb_cr3bp_torus* torus = &g_tori[0];
  const myb_cr3bp_torus_spectrum* s = &g_spectra[0];
  g_scratch = *torus;
  memcpy(g_scratch.states, unstable ? s->unstable_direction : s->stable_direction,
         sizeof g_scratch.states);
  const double sign = unstable ? 1.0 : -1.0;
  const double lambda = unstable ? s->unstable : 1.0 / s->stable;
  double worst = 0.0;
  *rate = 0.0;
  for (int k = 0; k < 7; ++k) {
    const double angle = 2.0 * PI * k / 7.0 + 0.1;
    double start[6];
    double on[6];
    double d0[6];
    double d1[6];
    double end[6];
    double elapsed = 0.0;
    int crossed = 1;
    MT_CHECK(myb_cr3bp_torus_manifold(torus, s, angle, unstable, 1, eps, start) == MYB_OK);
    MT_CHECK(myb_cr3bp_propagate(MU, start, sign * turns * torus->period, NULL, end, &elapsed,
                                 &crossed) == MYB_OK);
    MT_CHECK(crossed == 0 && elapsed == sign * turns * torus->period);
    myb_cr3bp_torus_state(torus, angle, on);
    MT_CHECK(myb_cr3bp_propagate(MU, on, sign * turns * torus->period, NULL, on, NULL, NULL) ==
             MYB_OK);
    myb_cr3bp_torus_state(&g_scratch, angle, d0);
    myb_cr3bp_torus_state(&g_scratch, angle + sign * turns * torus->rotation, d1);
    double gap[6];
    for (int i = 0; i < 6; ++i) {
      gap[i] = end[i] - on[i];
    }
    const double expected = eps * pow(lambda, turns) * position_norm(d1) / position_norm(d0);
    worst = fmax(worst, fabs(position_norm(gap) / expected - 1.0));
    *rate = fmax(*rate, pow(position_norm(gap) / eps, 1.0 / turns));
  }
  return worst;
}

static void test_growth(void) {
  const myb_cr3bp_torus* torus = &g_tori[0];
  const myb_cr3bp_torus_spectrum* s = &g_spectra[0];
  double rate_u1 = 0.0;
  double rate_s1 = 0.0;
  double rate_u2 = 0.0;
  double rate_s2 = 0.0;
  const double u1 = growth_miss(1, 1, 1e-9, &rate_u1);
  const double s1 = growth_miss(0, 1, 1e-9, &rate_s1);
  const double u2 = growth_miss(1, 2, 1e-12, &rate_u2);
  const double s2 = growth_miss(0, 2, 1e-12, &rate_s2);
  printf(
      "the manifolds of the torus of size 1e-3 (lambda %.4f): the distance from the torus "
      "grows per turn by up to %.4f (unstable, forward) and %.4f (stable, backward) over one "
      "turn from 1e-9, by %.4f and %.4f over two from 1e-12; against eps lambda^k |psi(theta "
      "+- k rho)| / |psi(theta)|: %.1e, %.1e and %.1e, %.1e\n",
      s->unstable, rate_u1, rate_s1, rate_u2, rate_s2, u1, s1, u2, s2);
  MT_CHECK(u1 < 2e-5 && s1 < 2e-5 && u2 < 2e-4 && s2 < 2e-4);
  MT_CHECK(fabs(rate_u1 / s->unstable - 1.0) < 1e-4 && fabs(rate_s1 * s->stable - 1.0) < 1e-4);
  /* the side +1 is that of growing x, on the mean over the circle */
  double along = 0.0;
  for (int k = 0; k < 12; ++k) {
    double start[6];
    double on[6];
    MT_CHECK(myb_cr3bp_torus_manifold(torus, s, 2.0 * PI * k / 12.0, 1, 1, 1e-6, start) == MYB_OK);
    myb_cr3bp_torus_state(torus, 2.0 * PI * k / 12.0, on);
    along += start[0] - on[0];
    MT_CHECK(fabs(position_norm(start) - position_norm(on)) <= 1e-6 * (1.0 + 1e-12));
  }
  MT_CHECK(along > 0.0);
}

/* The manifolds carried to the section of the Moon, x = 1 - mu: a point of the unstable
 * manifold of the torus at the angle theta and the displacement eps against the point where
 * the unstable manifold of the Lyapunov orbit, from its first state and with the same eps,
 * crosses it. The circle surrounds that state in the directions out of the plane (the
 * vertical pair of the planar orbit), so the two differ out of the plane as the size of the
 * torus and in it (y, vy) as its square, the motion in the plane being even in z. The stable
 * manifolds, backward: the mirror images (y, vx, vz -> -y, -vx, -vz, the symmetry of the
 * problem with time reversed) of the unstable ones. */
static void test_section(void) {
  double y0[6];
  memcpy(y0, LYAPUNOV, sizeof y0);
  const myb_cr3bp_section moon = {.axis = 0, .value = 1.0 - MU, .direction = 0};

  /* eps from 1e-4 down to 1e-6: the crossings spread over 0.6 of a turn (the growth of
   * a turn is lambda), and the error of the circle itself (1e-11, which grows as the
   * displacement does) is no more than 1e-5 of the displacement */
  enum { STEPS = 6, ANGLES = 12 };

  double orbit_u[STEPS][6];
  double orbit_s[STEPS][6];
  double y[6];
  MT_CHECK(myb_cr3bp_manifold(MU, y0, g_period, 0.0, 1, 1, 1e-6, y) == MYB_OK);
  const int side_u = y[0] > y0[0] ? 1 : -1;
  MT_CHECK(myb_cr3bp_manifold(MU, y0, g_period, 0.0, 0, 1, 1e-6, y) == MYB_OK);
  const int side_s = y[0] > y0[0] ? 1 : -1;
  double mirror = 0.0;
  double t_lo = 1e300;
  double t_hi = 0.0;
  for (int j = 0; j < STEPS; ++j) {
    const double eps = 1e-4 * pow(100.0, -j / (double)(STEPS - 1));
    double t = 0.0;
    double tb = 0.0;
    int crossed = 0;
    MT_CHECK(myb_cr3bp_manifold(MU, y0, g_period, 0.0, 1, side_u, eps, y) == MYB_OK);
    MT_CHECK(myb_cr3bp_propagate(MU, y, 4.0 * PI, &moon, orbit_u[j], &t, &crossed) == MYB_OK);
    MT_CHECK(crossed);
    MT_CHECK(myb_cr3bp_manifold(MU, y0, g_period, 0.0, 0, side_s, eps, y) == MYB_OK);
    MT_CHECK(myb_cr3bp_propagate(MU, y, -4.0 * PI, &moon, orbit_s[j], &tb, &crossed) == MYB_OK);
    MT_CHECK(crossed);
    mirror = fmax(mirror, fabs(orbit_s[j][1] + orbit_u[j][1]) +
                              fabs(orbit_s[j][4] - orbit_u[j][4]) + fabs(t + tb));
    t_lo = fmin(t_lo, t);
    t_hi = fmax(t_hi, t);
  }
  printf(
      "the unstable manifold of the Lyapunov orbit from its first state reaches x = 1 - mu "
      "after %.3f .. %.3f (vy %.4f .. %.4f); the stable one backward is its mirror image to "
      "%.0e\n",
      t_lo, t_hi, orbit_u[0][4], orbit_u[STEPS - 1][4], mirror);
  MT_CHECK(t_hi - t_lo > 0.5 * g_period && mirror < 1e-9);
  double planar[2];
  double heights[2];
  for (int s = 0; s < 2; ++s) {
    const myb_cr3bp_torus* torus = &g_tori[s];
    const myb_cr3bp_torus_spectrum* spectrum = &g_spectra[s];
    double in_plane = 0.0;
    double height = 0.0;
    double jacobi = 0.0;
    for (int k = 0; k < ANGLES; ++k) {
      for (int j = 0; j < STEPS; ++j) {
        const double angle = 2.0 * PI * k / (double)ANGLES;
        const double eps = 1e-4 * pow(100.0, -j / (double)(STEPS - 1));
        for (int unstable = 0; unstable <= 1; ++unstable) {
          const double* orbit = unstable ? orbit_u[j] : orbit_s[j];
          double at[6];
          double t = 0.0;
          int crossed = 0;
          MT_CHECK(myb_cr3bp_torus_manifold(torus, spectrum, angle, unstable, 1, eps, y) == MYB_OK);
          MT_CHECK(myb_cr3bp_propagate(MU, y, (unstable ? 4.0 : -4.0) * PI, &moon, at, &t,
                                       &crossed) == MYB_OK);
          MT_CHECK(crossed);
          in_plane = fmax(in_plane, hypot(at[1] - orbit[1], at[4] - orbit[4]));
          height = fmax(height, fabs(at[2]));
          jacobi = fmax(jacobi, fabs(myb_cr3bp_jacobi(MU, at) - myb_cr3bp_jacobi(MU, y)));
        }
      }
    }
    planar[s] = in_plane;
    heights[s] = height;
    printf(
        "  the torus of size %.0e: its manifolds (12 angles, 6 displacements, unstable and "
        "stable) cross the section within %.1e in (y, vy) of the orbit's, out of the plane "
        "up to z = %.2e; the Jacobi constant kept to %.0e\n",
        s == 0 ? 1e-3 : 2e-3, in_plane, height, jacobi);
    MT_CHECK(jacobi < 1e-13);
  }
  printf(
      "    from the size 1e-3 to 2e-3: in the plane %.3f times (the square), out of it %.4f "
      "(the size)\n",
      planar[1] / planar[0], heights[1] / heights[0]);
  MT_CHECK(fabs(heights[1] / heights[0] - 2.0) < 0.02);
  MT_CHECK(planar[1] / planar[0] > 3.5 && planar[1] / planar[0] < 4.5);
  MT_CHECK(planar[0] < 1e-3);
}

/* ---- Refused ------------------------------------------------------------------------------- */

static void test_refused(void) {
  const myb_cr3bp_torus* torus = &g_tori[0];
  const myb_cr3bp_torus_spectrum* s = &g_spectra[0];
  double y[6];
  MT_CHECK(myb_cr3bp_torus_stability(MU, torus, -1.0, &g_spectrum) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_torus_stability(0.0, torus, 0.0, &g_spectrum) == MYB_E_DOMAIN);
  /* a tolerance that the turn rho = 0.094 does not stand out of */
  MT_CHECK(myb_cr3bp_torus_stability(MU, torus, 0.01, &g_spectrum) == MYB_E_DOMAIN);
  g_scratch = *torus;
  g_scratch.points = 24;
  MT_CHECK(myb_cr3bp_torus_stability(MU, &g_scratch, 0.0, &g_spectrum) == MYB_E_DOMAIN);
  g_scratch = *torus;
  g_scratch.states[3] = (double)NAN;
  MT_CHECK(myb_cr3bp_torus_stability(MU, &g_scratch, 0.0, &g_spectrum) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_torus_manifold(torus, s, 0.0, 1, 0, 1e-6, y) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_torus_manifold(torus, s, 0.0, 1, 1, 0.0, y) == MYB_E_DOMAIN);
  MT_CHECK(myb_cr3bp_torus_manifold(torus, s, (double)NAN, 1, 1, 1e-6, y) == MYB_E_DOMAIN);
  /* a spectrum of another number of points */
  g_scratch = *torus;
  g_scratch.points = 23;
  MT_CHECK(myb_cr3bp_torus_manifold(&g_scratch, s, 0.0, 1, 1, 1e-6, y) == MYB_E_DOMAIN);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_propagate();
  test_lissajous();
  test_directions();
  test_growth();
  test_section();
  test_halo();
  test_stable();
  test_refused();
  return mt_finish("test_torus_stability");
}
