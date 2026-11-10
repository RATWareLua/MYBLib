/* SPDX-License-Identifier: BSD-3-Clause */
/* Small dense linear algebra: LU, Cholesky, Householder least squares, eigenvalues (companion
 * matrices with known roots, rotations, a symplectic matrix with the spread of a halo orbit's
 * monodromy, traces and determinants of random matrices), eigenvectors. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/num.h"

static double frand(mt_rng* rng) {
  return (double)mt_int(rng, -1000000, 1000000) / 1e6;
}

static void test_solvers(void) {
  mt_rng rng = {3};
  enum { N = 10 };
  double a[N * N];
  double a0[N * N];
  double x[N];
  double b[N];
  for (int i = 0; i < N * N; ++i) {
    a[i] = frand(&rng);
  }
  for (int i = 0; i < N; ++i) {
    a[i * N + i] += 4.0;
    x[i] = frand(&rng);
  }
  memcpy(a0, a, sizeof a);
  for (int i = 0; i < N; ++i) {
    b[i] = 0.0;
    for (int j = 0; j < N; ++j) {
      b[i] += a[i * N + j] * x[j];
    }
  }
  double y[N];
  memcpy(y, b, sizeof y);
  MT_CHECK(myb_lu_solve(N, a, y, 1) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < N; ++i) {
    worst = fmax(worst, fabs(y[i] - x[i]));
  }
  MT_CHECK(worst < 1e-14);
  /* Cholesky of A^T A */
  double s[N * N];
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      double t = 0.0;
      for (int k = 0; k < N; ++k) {
        t += a0[k * N + i] * a0[k * N + j];
      }
      s[i * N + j] = t;
    }
  }
  double rhs[N];
  for (int i = 0; i < N; ++i) {
    rhs[i] = 0.0;
    for (int j = 0; j < N; ++j) {
      rhs[i] += s[i * N + j] * x[j];
    }
  }
  MT_CHECK(myb_cholesky(N, s) == MYB_OK);
  myb_cholesky_solve(N, s, rhs);
  worst = 0.0;
  for (int i = 0; i < N; ++i) {
    worst = fmax(worst, fabs(rhs[i] - x[i]));
  }
  MT_CHECK(worst < 1e-13);
  /* least squares: 40 x 6 consistent, then the normal equations for a noisy one */
  enum { M = 40, P = 6 };
  double ls[M * P];
  double lb[M];
  double xp[P];
  for (int j = 0; j < P; ++j) {
    xp[j] = frand(&rng);
  }
  for (int i = 0; i < M; ++i) {
    lb[i] = 0.0;
    for (int j = 0; j < P; ++j) {
      ls[i * P + j] = frand(&rng);
      lb[i] += ls[i * P + j] * xp[j];
    }
  }
  double ls0[M * P];
  memcpy(ls0, ls, sizeof ls);
  double r = 1.0;
  MT_CHECK(myb_least_squares(M, P, ls, lb, &r) == MYB_OK);
  worst = 0.0;
  for (int j = 0; j < P; ++j) {
    worst = fmax(worst, fabs(lb[j] - xp[j]));
  }
  MT_CHECK(worst < 1e-14 && r < 1e-13);
  double noisy[M];
  for (int i = 0; i < M; ++i) {
    noisy[i] = frand(&rng);
  }
  double nb[M];
  memcpy(nb, noisy, sizeof nb);
  memcpy(ls, ls0, sizeof ls);
  MT_CHECK(myb_least_squares(M, P, ls, nb, NULL) == MYB_OK);
  double ne[P * P];
  double nr[P];
  for (int i = 0; i < P; ++i) {
    nr[i] = 0.0;
    for (int k = 0; k < M; ++k) {
      nr[i] += ls0[k * P + i] * noisy[k];
    }
    for (int j = 0; j < P; ++j) {
      double t = 0.0;
      for (int k = 0; k < M; ++k) {
        t += ls0[k * P + i] * ls0[k * P + j];
      }
      ne[i * P + j] = t;
    }
  }
  MT_CHECK(myb_cholesky(P, ne) == MYB_OK);
  myb_cholesky_solve(P, ne, nr);
  worst = 0.0;
  for (int j = 0; j < P; ++j) {
    worst = fmax(worst, fabs(nb[j] - nr[j]));
  }
  MT_CHECK(worst < 1e-13);
  /* errors */
  double sing[4] = {1.0, 2.0, 2.0, 4.0};
  double sb[2] = {1.0, 1.0};
  MT_CHECK(myb_lu_solve(2, sing, sb, 1) == MYB_E_DOMAIN);
  double neg[4] = {1.0, 2.0, 2.0, 1.0};
  MT_CHECK(myb_cholesky(2, neg) == MYB_E_DOMAIN);
  printf("LU, Cholesky, least squares: exact solutions recovered to 1e-14\n");
}

static int cmp_found(const double* re, const double* im, int n, double r, double i, double tol) {
  for (int k = 0; k < n; ++k) {
    if (fabs(re[k] - r) <= tol * fmax(1.0, fabs(r)) && fabs(im[k] - i) <= tol * fmax(1.0, fabs(r))) {
      return 1;
    }
  }
  return 0;
}

static void test_eigen(void) {
  /* companion matrix of (x-1)(x-2)...(x-6) */
  double c[7] = {1.0}; /* coefficients, highest first */
  int deg = 0;
  for (int root = 1; root <= 6; ++root) {
    for (int k = deg + 1; k >= 1; --k) {
      c[k] -= root * c[k - 1];
    }
    ++deg;
  }
  double comp[36];
  memset(comp, 0, sizeof comp);
  for (int j = 0; j < 6; ++j) {
    comp[j] = -c[j + 1];
  }
  for (int i = 1; i < 6; ++i) {
    comp[i * 6 + i - 1] = 1.0;
  }
  double re[12];
  double im[12];
  MT_CHECK(myb_eigenvalues(6, comp, re, im) == MYB_OK);
  for (int root = 1; root <= 6; ++root) {
    MT_CHECK(cmp_found(re, im, 6, root, 0.0, 1e-11));
  }
  /* rotation about an axis by 0.7 rad */
  const double t = 0.7;
  double rot[9] = {cos(t), -sin(t), 0.0, sin(t), cos(t), 0.0, 0.0, 0.0, 1.0};
  MT_CHECK(myb_eigenvalues(3, rot, re, im) == MYB_OK);
  MT_CHECK(cmp_found(re, im, 3, cos(t), sin(t), 1e-14) && cmp_found(re, im, 3, cos(t), -sin(t), 1e-14) &&
           cmp_found(re, im, 3, 1.0, 0.0, 1e-14));
  /* a symplectic-like spectrum: lambda, 1/lambda (lambda = 1500), mu, 1/mu, e^(+-i phi),
   * similar to a random matrix */
  mt_rng rng = {11};
  double tm[36];
  double d[36];
  memset(d, 0, sizeof d);
  const double lam = 1500.0;
  const double mu = 1.8;
  const double phi = 0.3;
  d[0] = lam;
  d[7] = 1.0 / lam;
  d[14] = mu;
  d[21] = 1.0 / mu;
  d[28] = cos(phi);
  d[29] = sin(phi);
  d[34] = -sin(phi);
  d[35] = cos(phi);
  for (int i = 0; i < 36; ++i) {
    tm[i] = frand(&rng) + (i % 7 == 0 ? 2.0 : 0.0);
  }
  /* A = T D T^-1: solve T^T X^T = (T D)^T, i.e. X = (T D) T^-1 */
  double td[36];
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double s = 0.0;
      for (int k = 0; k < 6; ++k) {
        s += tm[i * 6 + k] * d[k * 6 + j];
      }
      td[i * 6 + j] = s;
    }
  }
  double tt[36];
  double xt[36];
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      tt[i * 6 + j] = tm[j * 6 + i];
      xt[i * 6 + j] = td[j * 6 + i];
    }
  }
  MT_CHECK(myb_lu_solve(6, tt, xt, 6) == MYB_OK);
  double a[36];
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      a[i * 6 + j] = xt[j * 6 + i];
    }
  }
  double a0[36];
  memcpy(a0, a, sizeof a);
  MT_CHECK(myb_eigenvalues(6, a, re, im) == MYB_OK);
  MT_CHECK(cmp_found(re, im, 6, lam, 0.0, 1e-11));
  MT_CHECK(cmp_found(re, im, 6, 1.0 / lam, 0.0, 1e-8));
  MT_CHECK(cmp_found(re, im, 6, mu, 0.0, 1e-11));
  MT_CHECK(cmp_found(re, im, 6, 1.0 / mu, 0.0, 1e-11));
  MT_CHECK(cmp_found(re, im, 6, cos(phi), sin(phi), 1e-11));
  MT_CHECK(cmp_found(re, im, 6, cos(phi), -sin(phi), 1e-11));
  /* eigenvector of lambda: A v = lambda v */
  double v[6];
  double work[36];
  MT_CHECK(myb_eigenvector_real(6, a0, lam, v, work) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 6; ++i) {
    double s = 0.0;
    for (int j = 0; j < 6; ++j) {
      s += a0[i * 6 + j] * v[j];
    }
    worst = fmax(worst, fabs(s - lam * v[i]));
  }
  MT_CHECK(worst < 1e-10 * lam);
  /* random matrices: trace and determinant */
  for (int trial = 0; trial < 50; ++trial) {
    const int n = 2 + trial % 11;
    double m[144];
    double m2[144];
    for (int i = 0; i < n * n; ++i) {
      m[i] = frand(&rng);
    }
    memcpy(m2, m, sizeof(double) * (size_t)(n * n));
    double trace = 0.0;
    for (int i = 0; i < n; ++i) {
      trace += m[i * n + i];
    }
    MT_CHECK(myb_eigenvalues((size_t)n, m, re, im) == MYB_OK);
    double sr = 0.0;
    double si = 0.0;
    double pr = 1.0;
    double pi = 0.0;
    for (int k = 0; k < n; ++k) {
      sr += re[k];
      si += im[k];
      const double nr = pr * re[k] - pi * im[k];
      pi = pr * im[k] + pi * re[k];
      pr = nr;
    }
    /* determinant by LU (product of the pivots with the sign of the permutation is not
     * returned: compare |det|) */
    double lu[144];
    memcpy(lu, m2, sizeof(double) * (size_t)(n * n));
    double ones[12];
    for (int i = 0; i < n; ++i) {
      ones[i] = 1.0;
    }
    MT_CHECK(myb_lu_solve((size_t)n, lu, ones, 1) == MYB_OK);
    double det = 1.0;
    for (int i = 0; i < n; ++i) {
      det *= lu[i * n + i];
    }
    MT_CHECK(fabs(sr - trace) < 1e-12 * n && fabs(si) < 1e-12 * n);
    MT_CHECK(fabs(fabs(pr) - fabs(det)) < 1e-11 * fmax(1.0, fabs(det)) && fabs(pi) < 1e-11);
  }
  /* Jordan block: double eigenvalue 1 */
  double j2[4] = {1.0, 1.0, 0.0, 1.0};
  MT_CHECK(myb_eigenvalues(2, j2, re, im) == MYB_OK);
  MT_CHECK(fabs(re[0] - 1.0) < 1e-15 && fabs(re[1] - 1.0) < 1e-15);
  double bad[4] = {1.0, (double)NAN, 0.0, 1.0};
  MT_CHECK(myb_eigenvalues(2, bad, re, im) == MYB_E_DOMAIN);
  printf("eigenvalues: companion roots, rotation, spectrum 1500 .. 1/1500, 50 random matrices\n");
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_solvers();
  test_eigen();
  return mt_finish("test_linalg");
}
