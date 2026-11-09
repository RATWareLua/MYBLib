/* SPDX-License-Identifier: BSD-3-Clause */
/* Ellipsoidal harmonics:
 *  - every basis function is harmonic (4th-order finite-difference Laplacian), which holds only
 *    if the Lamé functions are solutions of Lamé's equation and R(lambda) is right;
 *  - on the reference ellipsoid a basis function equals E(mu) E(nu) / a^(2n) (Niven's Cartesian
 *    form against the ellipsoidal coordinates found from the confocal cubic);
 *  - the acceleration is the gradient of the potential (finite differences);
 *  - a homogeneous box (exact field: myb_polyhedron) fitted on its circumscribed ellipsoid:
 *    inside the Brillouin sphere, outside the ellipsoid, the error falls with the degree, while
 *    the exact spherical-harmonic series of the same box (volume moments by Gauss-Legendre
 *    quadrature, exact for the polynomials) does not converge there. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/grav.h"

static const double box_half[3] = {2.0, 1.2, 0.8};

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/* The roots s = lambda^2 > mu^2 > nu^2 of x^2/s + y^2/(s-h^2) + z^2/(s-k^2) = 1 by bisection on
 * (k^2, inf), (h^2, k^2), (0, h^2). */
static double confocal(double s, const double x[3], double h2, double k2) {
  return s * (s - h2) * (s - k2) - x[0] * x[0] * (s - h2) * (s - k2) - x[1] * x[1] * s * (s - k2) -
         x[2] * x[2] * s * (s - h2);
}

static double bisect(double lo, double hi, const double x[3], double h2, double k2) {
  double flo = confocal(lo, x, h2, k2);
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (mid <= lo || mid >= hi) {
      break;
    }
    const double fm = confocal(mid, x, h2, k2);
    if ((fm < 0.0) == (flo < 0.0)) {
      lo = mid;
      flo = fm;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

static void test_basis(void) {
  const double axes[3] = {3.0, 2.2, 1.5};
  const int N = 16;
  myb_ellipsoidal_field* f = NULL;
  MT_CHECK(myb_ellipsoidal_field_create(1.0, axes, N, NULL, NULL, &f) == MYB_OK);
  if (f == NULL) {
    return;
  }
  const int count = (N + 1) * (N + 1);
  const double h2 = axes[0] * axes[0] - axes[1] * axes[1];
  const double k2 = axes[0] * axes[0] - axes[2] * axes[2];
  mt_rng rng = {17};
  double worst_lap = 0.0;
  double worst_niven = 0.0;
  for (int trial = 0; trial < 12; ++trial) {
    double d[3];
    for (int i = 0; i < 3; ++i) {
      d[i] = (double)mt_int(&rng, -1000, 1000) / 1000.0;
    }
    /* a point on the reference ellipsoid along d, then 1.3 times farther for the Laplacian */
    const double q = sqrt(d[0] * d[0] / (axes[0] * axes[0]) + d[1] * d[1] / (axes[1] * axes[1]) +
                          d[2] * d[2] / (axes[2] * axes[2]));
    double s[3];
    double p[3];
    for (int i = 0; i < 3; ++i) {
      s[i] = d[i] / q;
      p[i] = 1.3 * s[i];
    }
    double base[289];
    double plus[289];
    double minus[289];
    double lap[289];
    double stencil_max[289];
    memset(lap, 0, sizeof lap);
    memset(stencil_max, 0, sizeof stencil_max);
    MT_CHECK(myb_ellipsoidal_field_basis(f, p, base) == MYB_OK);
    const double step = 2e-3;
    double plus2[289];
    double minus2[289];
    for (int c = 0; c < 3; ++c) {
      double x[3] = {p[0], p[1], p[2]};
      x[c] = p[c] + step;
      MT_CHECK(myb_ellipsoidal_field_basis(f, x, plus) == MYB_OK);
      x[c] = p[c] - step;
      MT_CHECK(myb_ellipsoidal_field_basis(f, x, minus) == MYB_OK);
      x[c] = p[c] + 2.0 * step;
      MT_CHECK(myb_ellipsoidal_field_basis(f, x, plus2) == MYB_OK);
      x[c] = p[c] - 2.0 * step;
      MT_CHECK(myb_ellipsoidal_field_basis(f, x, minus2) == MYB_OK);
      for (int k = 0; k < count; ++k) {
        lap[k] += (-plus2[k] + 16.0 * plus[k] - 30.0 * base[k] + 16.0 * minus[k] - minus2[k]) /
                  (12.0 * step * step);
        stencil_max[k] = fmax(stencil_max[k], fmax(fmax(fabs(plus[k]), fabs(minus[k])),
                                                   fmax(fabs(plus2[k]), fabs(minus2[k]))));
      }
    }
    const double r = norm3(p);
    for (int k = 0; k < count; ++k) {
      /* against the size of the second derivatives, |V| (n + 1)(n + 2) / r^2, |V| the largest
       * value of the stencil (a point may lie near a nodal surface of the term) */
      const int n = (int)floor(sqrt((double)k));
      const double vmax = fmax(fabs(base[k]), fmax(stencil_max[k], 1e-300));
      const double scale = vmax * (n + 1.0) * (n + 2.0) / (r * r);
      worst_lap = fmax(worst_lap, fabs(lap[k]) / scale);
    }
    /* Niven: on the surface, basis = sign * E(mu) E(nu) / a^(2n) */
    const double x0 = s[0];
    const double y0 = s[1];
    const double z0 = s[2];
    const double mu2 = bisect(h2, k2, s, h2, k2);
    const double nu2 = bisect(0.0, h2, s, h2, k2);
    double surf[289];
    MT_CHECK(myb_ellipsoidal_field_basis(f, s, surf) == MYB_OK);
    for (int n = 0; n <= N; ++n) {
      for (int pp = 0; pp <= 2 * n; ++pp) {
        double emu = 0.0;
        double enu = 0.0;
        MT_CHECK(myb_ellipsoidal_lame(f, n, pp, sqrt(mu2), &emu) == MYB_OK);
        MT_CHECK(myb_ellipsoidal_lame(f, n, pp, sqrt(nu2), &enu) == MYB_OK);
        /* parity of the classes: x^ex via the sign of nu, y for L and N, z for M and N */
        const int cls_index = pp; /* find the class from the counts K, L, M, N */
        const int nk = n / 2 + 1;
        const int nl = n >= 1 ? (n - 1) / 2 + 1 : 0;
        const int nm = nl;
        int cls = 0;
        if (cls_index >= nk + nl + nm) {
          cls = 3;
        } else if (cls_index >= nk + nl) {
          cls = 2;
        } else if (cls_index >= nk) {
          cls = 1;
        }
        const int o = cls == 0 ? 0 : (cls == 3 ? 2 : 1);
        const int ex = (n - o) % 2;
        double sign = 1.0;
        if (ex && x0 < 0.0) {
          sign = -sign;
        }
        if ((cls == 1 || cls == 3) && y0 < 0.0) {
          sign = -sign;
        }
        if ((cls == 2 || cls == 3) && z0 < 0.0) {
          sign = -sign;
        }
        double a2n = 1.0;
        for (int i = 0; i < 2 * n; ++i) {
          a2n *= axes[0];
        }
        const double expected = sign * emu * enu / a2n;
        const double v = surf[n * n + pp];
        worst_niven = fmax(worst_niven, fabs(v - expected) / (fabs(expected) + 1e-12));
      }
    }
  }
  printf("ellipsoidal basis to degree %d: Laplacian %.1e of the second derivatives, Niven's form "
         "vs E(mu) E(nu) %.1e\n",
         N, worst_lap, worst_niven);
  MT_CHECK(worst_lap < 1e-6);
  MT_CHECK(worst_niven < 1e-9);
  /* errors */
  double v = 0.0;
  MT_CHECK(myb_ellipsoidal_field_eval(f, (const double[3]){1.0, 0.5, 0.2}, &v, NULL) ==
           MYB_E_DOMAIN);
  myb_ellipsoidal_field* g = NULL;
  MT_CHECK(myb_ellipsoidal_field_create(1.0, (const double[3]){2.0, 2.0, 1.0}, 2, NULL, NULL, &g) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_ellipsoidal_field_create(1.0, axes, 21, NULL, NULL, &g) == MYB_E_DOMAIN);
  myb_ellipsoidal_field_destroy(f);
  myb_ellipsoidal_field_destroy(NULL);
}

/* ---- The box --------------------------------------------------------------------------------- */

static myb_polyhedron* make_box(void) {
  double v[8][3];
  for (int k = 0; k < 8; ++k) {
    v[k][0] = (k & 1) ? box_half[0] : -box_half[0];
    v[k][1] = (k & 2) ? box_half[1] : -box_half[1];
    v[k][2] = (k & 4) ? box_half[2] : -box_half[2];
  }
  static const size_t faces[12][3] = {{0, 2, 3}, {0, 3, 1}, {4, 5, 7}, {4, 7, 6},
                                      {0, 1, 5}, {0, 5, 4}, {2, 6, 7}, {2, 7, 3},
                                      {0, 4, 6}, {0, 6, 2}, {1, 3, 7}, {1, 7, 5}};
  myb_polyhedron* p = NULL;
  MT_CHECK(myb_polyhedron_create(&v[0][0], 8, &faces[0][0], 12, 1.0, NULL, &p) == MYB_OK);
  return p;
}

static int solve(int n, double* a, double* b) { /* Cholesky, a n x n row-major */
  for (int j = 0; j < n; ++j) {
    double d = a[j * n + j];
    for (int k = 0; k < j; ++k) {
      d -= a[j * n + k] * a[j * n + k];
    }
    if (!(d > 0.0)) {
      return 0;
    }
    a[j * n + j] = sqrt(d);
    for (int i = j + 1; i < n; ++i) {
      double s = a[i * n + j];
      for (int k = 0; k < j; ++k) {
        s -= a[i * n + k] * a[j * n + k];
      }
      a[i * n + j] = s / a[j * n + j];
    }
  }
  for (int i = 0; i < n; ++i) {
    double s = b[i];
    for (int k = 0; k < i; ++k) {
      s -= a[i * n + k] * b[k];
    }
    b[i] = s / a[i * n + i];
  }
  for (int i = n - 1; i >= 0; --i) {
    double s = b[i];
    for (int k = i + 1; k < n; ++k) {
      s -= a[k * n + i] * b[k];
    }
    b[i] = s / a[i * n + i];
  }
  return 1;
}

/* Gauss-Legendre nodes on [-1, 1] (n <= 32) by Newton on P_n. */
static void gauss(int n, double* x, double* w) {
  for (int k = 0; k < n; ++k) {
    double z = cos(3.141592653589793 * (k + 0.75) / (n + 0.5));
    double dp = 1.0;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0;
      double p1 = z;
      for (int j = 2; j <= n; ++j) {
        const double p2 = ((2.0 * j - 1.0) * z * p1 - (j - 1.0) * p0) / j;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (z * p1 - p0) / (z * z - 1.0);
      const double dz = p1 / dp;
      z -= dz;
      if (fabs(dz) < 1e-16) {
        break;
      }
    }
    x[k] = z;
    w[k] = 2.0 / ((1.0 - z * z) * dp * dp);
  }
}

/* Spherical-harmonic coefficients of the box (G rho = 1) by Gauss-Legendre in its volume. */
static myb_gravity_field* box_harmonics(int N, double R) {
  enum { Q = 16 };
  double gx[Q];
  double gw[Q];
  gauss(Q, gx, gw);
  static double c[(25 + 1) * (25 + 2) / 2];
  static double s[(25 + 1) * (25 + 2) / 2];
  memset(c, 0, sizeof c);
  memset(s, 0, sizeof s);
  const double mass = 8.0 * box_half[0] * box_half[1] * box_half[2];
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j) {
      for (int l = 0; l < Q; ++l) {
        const double p[3] = {box_half[0] * gx[i], box_half[1] * gx[j], box_half[2] * gx[l]};
        const double wm = gw[i] * gw[j] * gw[l] * box_half[0] * box_half[1] * box_half[2];
        const double r = norm3(p);
        const double t = p[2] / r;
        const double u = sqrt(p[0] * p[0] + p[1] * p[1]) / r;
        const double lam = atan2(p[1], p[0]);
        double pmm = 1.0;
        for (int m = 0; m <= N; ++m) {
          if (m == 1) {
            pmm = sqrt(3.0) * u;
          } else if (m > 1) {
            pmm *= sqrt((2.0 * m + 1.0) / (2.0 * m)) * u;
          }
          double p0 = pmm;
          double p1 = m < N ? sqrt(2.0 * m + 3.0) * t * pmm : 0.0;
          double rn = pow(r / R, m);
          for (int n = m; n <= N; ++n) {
            double pn;
            if (n == m) {
              pn = p0;
            } else if (n == m + 1) {
              pn = p1;
            } else {
              const double aa =
                  sqrt((2.0 * n - 1.0) * (2.0 * n + 1.0) / ((double)(n - m) * (n + m)));
              const double bb = sqrt((2.0 * n + 1.0) * (n + m - 1.0) * (n - m - 1.0) /
                                     ((double)(n - m) * (n + m) * (2.0 * n - 3.0)));
              pn = aa * t * p1 - bb * p0;
              p0 = p1;
              p1 = pn;
            }
            const double wt = wm / mass * rn / (2.0 * n + 1.0) * pn;
            c[n * (n + 1) / 2 + m] += wt * cos(m * lam);
            s[n * (n + 1) / 2 + m] += wt * sin(m * lam);
            rn *= r / R;
          }
        }
      }
    }
  }
  myb_gravity_field* f = NULL;
  MT_CHECK(myb_gravity_field_create(mass, R, N, c, s, NULL, &f) == MYB_OK);
  return f;
}

static void test_box(void) {
  myb_polyhedron* box = make_box();
  if (box == NULL) {
    return;
  }
  /* the ellipsoid through the corners, enlarged by 5 %: the body strictly inside */
  const double axes[3] = {1.05 * sqrt(3.0) * box_half[0], 1.05 * sqrt(3.0) * box_half[1],
                          1.05 * sqrt(3.0) * box_half[2]};
  const double brillouin = norm3(box_half);
  /* outside the ellipsoid, inside the Brillouin sphere; and one point outside the sphere */
  const double inner[3][3] = {{0.0, 0.0, 1.6}, {0.0, 2.25, 0.0}, {1.0, 1.5, 1.0}};
  const double outer[3] = {3.6, 0.3, 0.2};
  for (int k = 0; k < 3; ++k) {
    MT_CHECK(norm3(inner[k]) < brillouin);
  }
  const double mass = 8.0 * box_half[0] * box_half[1] * box_half[2];
  const int degrees[4] = {4, 8, 12, 16};
  double err_e[4];
  double err_s[4];
  double err_out[4];
  double err_acc[4];
  for (int di = 0; di < 4; ++di) {
    const int N = degrees[di];
    const int count = (N + 1) * (N + 1);
    myb_ellipsoidal_field* f = NULL;
    MT_CHECK(myb_ellipsoidal_field_create(mass, axes, N, NULL, NULL, &f) == MYB_OK);
    if (f == NULL) {
      break;
    }
    /* least squares on the reference ellipsoid */
    const int nt = 2 * N + 6;
    const int np = 2 * nt;
    double gx[64];
    double gw[64];
    gauss(nt, gx, gw);
    static double nmat[289 * 289];
    static double rhs[289];
    static double basis[289];
    memset(nmat, 0, sizeof(double) * (size_t)count * (size_t)count);
    memset(rhs, 0, sizeof(double) * (size_t)count);
    for (int i = 0; i < nt; ++i) {
      const double ct = gx[i];
      const double st = sqrt(1.0 - ct * ct);
      for (int j = 0; j < np; ++j) {
        const double ph = 2.0 * 3.141592653589793 * (j + 0.5) / np;
        const double x[3] = {axes[0] * st * cos(ph), axes[1] * st * sin(ph), axes[2] * ct};
        double u = 0.0;
        MT_CHECK(myb_polyhedron_eval(box, x, &u, NULL, NULL) == MYB_OK);
        MT_CHECK(myb_ellipsoidal_field_basis(f, x, basis) == MYB_OK);
        const double w = gw[i];
        const double scale = mass / axes[0];
        for (int a = 0; a < count; ++a) {
          rhs[a] += w * basis[a] * u / scale;
          for (int b = 0; b <= a; ++b) {
            nmat[a * count + b] += w * basis[a] * basis[b];
          }
        }
      }
    }
    for (int a = 0; a < count; ++a) {
      for (int b = 0; b < a; ++b) {
        nmat[b * count + a] = nmat[a * count + b];
      }
    }
    MT_CHECK(solve(count, nmat, rhs));
    MT_CHECK(myb_ellipsoidal_field_set_coefficients(f, rhs) == MYB_OK);
    myb_gravity_field* sh = box_harmonics(N, brillouin);
    myb_gravity_work* w = NULL;
    MT_CHECK(myb_gravity_work_create(N, NULL, &w) == MYB_OK);
    err_e[di] = 0.0;
    err_s[di] = 0.0;
    err_acc[di] = 0.0;
    for (int k = 0; k < 3; ++k) {
      double u = 0.0;
      double a[3];
      double ue = 0.0;
      double ae[3];
      double us = 0.0;
      MT_CHECK(myb_polyhedron_eval(box, inner[k], &u, a, NULL) == MYB_OK);
      MT_CHECK(myb_ellipsoidal_field_eval(f, inner[k], &ue, ae) == MYB_OK);
      MT_CHECK(myb_gravity_field_eval(sh, w, inner[k], 0, -1, &us, NULL) == MYB_OK);
      err_e[di] = fmax(err_e[di], fabs(ue - u) / u);
      err_s[di] = fmax(err_s[di], fabs(us - u) / u);
      const double da[3] = {ae[0] - a[0], ae[1] - a[1], ae[2] - a[2]};
      err_acc[di] = fmax(err_acc[di], norm3(da) / norm3(a));
    }
    double u = 0.0;
    double ue = 0.0;
    MT_CHECK(myb_polyhedron_eval(box, outer, &u, NULL, NULL) == MYB_OK);
    MT_CHECK(myb_ellipsoidal_field_eval(f, outer, &ue, NULL) == MYB_OK);
    err_out[di] = fabs(ue - u) / u;
    if (di == 1) {
      /* acceleration = gradient of the potential */
      const double* p = inner[2];
      double ae[3];
      double ue0 = 0.0;
      MT_CHECK(myb_ellipsoidal_field_eval(f, p, &ue0, ae) == MYB_OK);
      double worst = 0.0;
      for (int c = 0; c < 3; ++c) {
        double xp[3] = {p[0], p[1], p[2]};
        double xm[3] = {p[0], p[1], p[2]};
        xp[c] += 1e-5;
        xm[c] -= 1e-5;
        double vp = 0.0;
        double vm = 0.0;
        MT_CHECK(myb_ellipsoidal_field_eval(f, xp, &vp, NULL) == MYB_OK);
        MT_CHECK(myb_ellipsoidal_field_eval(f, xm, &vm, NULL) == MYB_OK);
        worst = fmax(worst, fabs((vp - vm) / 2e-5 - ae[c]) / norm3(ae));
      }
      printf("ellipsoidal acceleration vs differences of the potential: %.1e\n", worst);
      MT_CHECK(worst < 1e-8);
    }
    myb_gravity_work_destroy(w);
    myb_gravity_field_destroy(sh);
    myb_ellipsoidal_field_destroy(f);
  }
  printf("box 4 x 2.4 x 1.6 inside its Brillouin sphere, outside the ellipsoid:\n"
         "  degree   ellipsoidal (potential, acceleration)   spherical   | outside the sphere\n");
  for (int di = 0; di < 4; ++di) {
    printf("  %6d   %9.1e %9.1e                       %9.1e   | %9.1e\n", degrees[di], err_e[di],
           err_acc[di], err_s[di], err_out[di]);
  }
  /* measured: 8.9e-3 -> 6.5e-5 (ellipsoidal) while the spherical series grows 3.8e-2 -> 0.26;
   * the corners of the box, close to the reference ellipsoid, set the slow rate */
  MT_CHECK(err_e[3] < 0.02 * err_e[0] && err_e[3] < 2e-4);
  MT_CHECK(err_acc[3] < 0.1 * err_acc[0]);
  MT_CHECK(err_s[3] > err_s[0] && err_s[3] > 100.0 * err_e[3]);
  MT_CHECK(err_out[3] < 0.05 * err_out[0] && err_out[3] < 1e-3);
  myb_polyhedron_destroy(box);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_basis();
  test_box();
  return mt_finish("test_ellipsoidal");
}
