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

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_lagrange_and_jacobian();
  test_frames();
  test_catalogue();
  test_richardson();
  test_manifold();
  return mt_finish("test_cr3bp");
}
