/* SPDX-License-Identifier: BSD-3-Clause */
/* Small dense linear algebra for the solvers of myb (differential correction, estimation,
 * stability of periodic orbits): LU with partial pivoting, Cholesky, Householder QR least
 * squares, eigenvalues of real non-symmetric matrices by balancing, Householder reduction to
 * Hessenberg form and the implicit double-shift QR iteration of Francis (G. H. Golub,
 * C. F. Van Loan, Matrix Computations, 4th ed., 2013, Alg. 7.4.2, 7.5.1, 7.5.2; balancing of
 * B. N. Parlett, C. Reinsch, Numer. Math. 13, 293, 1969), eigenvectors by inverse iteration.
 * Every loop runs a fixed sequence of operations: results are reproducible bit for bit. */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "myb/num.h"
#include "num/linalg_internal.h"

enum { MAX_N = 64 };

static const double eps = 0x1p-52;

myb_status myb_lu_solve(size_t n, double* a, double* b, size_t nrhs) {
  MYB_CHECK(a != NULL && b != NULL && n <= MAX_N);
  return myb_internal_lu_solve(n, a, b, nrhs);
}

myb_status myb_internal_lu_solve(size_t n, double* a, double* b, size_t nrhs) {
  for (size_t k = 0; k < n; ++k) {
    /* pivot */
    size_t p = k;
    double best = fabs(a[k * n + k]);
    for (size_t i = k + 1; i < n; ++i) {
      if (fabs(a[i * n + k]) > best) {
        best = fabs(a[i * n + k]);
        p = i;
      }
    }
    if (!(best > 0.0) || !isfinite(best)) {
      return MYB_E_DOMAIN;
    }
    if (p != k) {
      for (size_t j = 0; j < n; ++j) {
        const double t = a[k * n + j];
        a[k * n + j] = a[p * n + j];
        a[p * n + j] = t;
      }
      for (size_t j = 0; j < nrhs; ++j) {
        const double t = b[k * nrhs + j];
        b[k * nrhs + j] = b[p * nrhs + j];
        b[p * nrhs + j] = t;
      }
    }
    for (size_t i = k + 1; i < n; ++i) {
      const double f = a[i * n + k] / a[k * n + k];
      a[i * n + k] = f;
      for (size_t j = k + 1; j < n; ++j) {
        a[i * n + j] -= f * a[k * n + j];
      }
      for (size_t j = 0; j < nrhs; ++j) {
        b[i * nrhs + j] -= f * b[k * nrhs + j];
      }
    }
  }
  for (size_t jj = 0; jj < nrhs; ++jj) {
    for (size_t ii = n; ii-- > 0;) {
      double s = b[ii * nrhs + jj];
      for (size_t j = ii + 1; j < n; ++j) {
        s -= a[ii * n + j] * b[j * nrhs + jj];
      }
      b[ii * nrhs + jj] = s / a[ii * n + ii];
    }
  }
  return MYB_OK;
}

myb_status myb_cholesky(size_t n, double* a) {
  MYB_CHECK(a != NULL && n <= MAX_N);
  for (size_t j = 0; j < n; ++j) {
    double d = a[j * n + j];
    for (size_t k = 0; k < j; ++k) {
      d -= a[j * n + k] * a[j * n + k];
    }
    if (!(d > 0.0) || !isfinite(d)) {
      return MYB_E_DOMAIN;
    }
    const double l = sqrt(d);
    a[j * n + j] = l;
    for (size_t i = j + 1; i < n; ++i) {
      double s = a[i * n + j];
      for (size_t k = 0; k < j; ++k) {
        s -= a[i * n + k] * a[j * n + k];
      }
      a[i * n + j] = s / l;
    }
    for (size_t i = 0; i < j; ++i) {
      a[i * n + j] = 0.0; /* the upper triangle of L */
    }
  }
  return MYB_OK;
}

void myb_cholesky_solve(size_t n, const double* l, double* b) {
  MYB_CHECK(l != NULL && b != NULL && n <= MAX_N);
  for (size_t i = 0; i < n; ++i) {
    double s = b[i];
    for (size_t k = 0; k < i; ++k) {
      s -= l[i * n + k] * b[k];
    }
    b[i] = s / l[i * n + i];
  }
  for (size_t i = n; i-- > 0;) {
    double s = b[i];
    for (size_t k = i + 1; k < n; ++k) {
      s -= l[k * n + i] * b[k];
    }
    b[i] = s / l[i * n + i];
  }
}

myb_status myb_least_squares(size_t m, size_t n, double* a, double* b, double* residual_norm) {
  MYB_CHECK(a != NULL && b != NULL && n <= MAX_N && m >= n);
  double diag[MAX_N];
  for (size_t k = 0; k < n; ++k) {
    /* Householder vector for column k, rows k .. m-1 */
    double scale = 0.0;
    for (size_t i = k; i < m; ++i) {
      scale = fmax(scale, fabs(a[i * n + k]));
    }
    if (!(scale > 0.0) || !isfinite(scale)) {
      return MYB_E_DOMAIN; /* rank deficient */
    }
    double norm2 = 0.0;
    for (size_t i = k; i < m; ++i) {
      const double t = a[i * n + k] / scale;
      norm2 += t * t;
    }
    double alpha = scale * sqrt(norm2);
    if (a[k * n + k] > 0.0) {
      alpha = -alpha;
    }
    /* v = x - alpha e1, stored in place; beta = 1 / (v.x - ...) */
    a[k * n + k] -= alpha;
    double vv = 0.0;
    for (size_t i = k; i < m; ++i) {
      vv += a[i * n + k] * a[i * n + k];
    }
    diag[k] = alpha;
    if (vv > 0.0) {
      for (size_t j = k + 1; j < n; ++j) {
        double s = 0.0;
        for (size_t i = k; i < m; ++i) {
          s += a[i * n + k] * a[i * n + j];
        }
        const double f = 2.0 * s / vv;
        for (size_t i = k; i < m; ++i) {
          a[i * n + j] -= f * a[i * n + k];
        }
      }
      double s = 0.0;
      for (size_t i = k; i < m; ++i) {
        s += a[i * n + k] * b[i];
      }
      const double f = 2.0 * s / vv;
      for (size_t i = k; i < m; ++i) {
        b[i] -= f * a[i * n + k];
      }
    }
  }
  if (residual_norm != NULL) {
    double r = 0.0;
    for (size_t i = n; i < m; ++i) {
      r += b[i] * b[i];
    }
    *residual_norm = sqrt(r);
  }
  /* R x = (Q^T b)[0..n-1], R upper with diagonal diag */
  for (size_t i = n; i-- > 0;) {
    double s = b[i];
    for (size_t j = i + 1; j < n; ++j) {
      s -= a[i * n + j] * b[j];
    }
    if (!(fabs(diag[i]) > 0.0)) {
      return MYB_E_DOMAIN;
    }
    b[i] = s / diag[i];
  }
  return MYB_OK;
}

/* ---- Eigenvalues ------------------------------------------------------------------------------ */

/* Balancing by powers of 2 (no rounding): rows and columns of comparable norms. */
static void balance(size_t n, double* a) {
  /* bound: at most 64 sweeps; each sweep changes scales by factors of 2 */
  for (int sweep = 0; sweep < 64; ++sweep) {
    int done = 1;
    for (size_t i = 0; i < n; ++i) {
      double c = 0.0;
      double r = 0.0;
      for (size_t j = 0; j < n; ++j) {
        if (j != i) {
          c += fabs(a[j * n + i]);
          r += fabs(a[i * n + j]);
        }
      }
      if (c == 0.0 || r == 0.0) {
        continue;
      }
      double f = 1.0;
      const double s = c + r;
      double cc = c;
      /* bound: the exponent range, at most 2100 halvings or doublings */
      for (int it = 0; it < 2100 && cc < r * 0.5; ++it) {
        cc *= 4.0;
        f *= 2.0;
      }
      for (int it = 0; it < 2100 && cc >= r * 2.0; ++it) {
        cc *= 0.25;
        f *= 0.5;
      }
      if ((c * f * f + r) / f < 0.95 * s && f != 1.0) {
        done = 0;
        for (size_t j = 0; j < n; ++j) {
          a[i * n + j] /= f;
          a[j * n + i] *= f;
        }
      }
    }
    if (done) {
      break;
    }
  }
}

/* Householder reduction to upper Hessenberg form (similarity); v: n doubles of work. */
static void hessenberg(size_t n, double* a, double* v) {
  for (size_t k = 0; k + 2 < n; ++k) {
    double scale = 0.0;
    for (size_t i = k + 1; i < n; ++i) {
      scale = fmax(scale, fabs(a[i * n + k]));
    }
    if (scale == 0.0) {
      continue;
    }
    double norm2 = 0.0;
    for (size_t i = k + 1; i < n; ++i) {
      v[i] = a[i * n + k] / scale;
      norm2 += v[i] * v[i];
    }
    double alpha = sqrt(norm2);
    if (v[k + 1] > 0.0) {
      alpha = -alpha;
    }
    v[k + 1] -= alpha;
    double vv = 0.0;
    for (size_t i = k + 1; i < n; ++i) {
      vv += v[i] * v[i];
    }
    if (vv == 0.0) {
      continue;
    }
    /* A = P A P, P = I - 2 v v^T / v.v on indices k+1 .. n-1 */
    for (size_t j = 0; j < n; ++j) {
      double s = 0.0;
      for (size_t i = k + 1; i < n; ++i) {
        s += v[i] * a[i * n + j];
      }
      const double f = 2.0 * s / vv;
      for (size_t i = k + 1; i < n; ++i) {
        a[i * n + j] -= f * v[i];
      }
    }
    for (size_t i = 0; i < n; ++i) {
      double s = 0.0;
      for (size_t j = k + 1; j < n; ++j) {
        s += a[i * n + j] * v[j];
      }
      const double f = 2.0 * s / vv;
      for (size_t j = k + 1; j < n; ++j) {
        a[i * n + j] -= f * v[j];
      }
    }
    for (size_t i = k + 2; i < n; ++i) {
      a[i * n + k] = 0.0;
    }
  }
}

/* Eigenvalues of the 2 x 2 block [[p, q], [r, s]]. */
static void eig2(double p, double q, double r, double s, double* re, double* im) {
  const double half = 0.5 * (p - s);
  const double disc = half * half + q * r;
  if (disc >= 0.0) {
    const double root = sqrt(disc);
    const double z = half + (half >= 0.0 ? root : -root);
    re[0] = s + z;
    re[1] = z != 0.0 ? s - q * r / z : s + z;
    im[0] = 0.0;
    im[1] = 0.0;
  } else {
    re[0] = s + half;
    re[1] = s + half;
    im[0] = sqrt(-disc);
    im[1] = -im[0];
  }
}

/* Applies the reflector I - 2 u u^T / u.u (u of length len at rows/columns k ..) on the
 * active block [lo, hi] of the Hessenberg matrix from both sides. */
static void reflect(size_t n, double* h, size_t k, int len, const double* u, size_t lo,
                    size_t hi) {
  double uu = 0.0;
  for (int i = 0; i < len; ++i) {
    uu += u[i] * u[i];
  }
  if (uu == 0.0) {
    return;
  }
  const size_t col0 = k > lo ? k - 1 : lo;
  for (size_t j = col0; j <= hi; ++j) {
    double s = 0.0;
    for (int i = 0; i < len; ++i) {
      s += u[i] * h[(k + (size_t)i) * n + j];
    }
    const double f = 2.0 * s / uu;
    for (int i = 0; i < len; ++i) {
      h[(k + (size_t)i) * n + j] -= f * u[i];
    }
  }
  const size_t row1 = k + (size_t)len < hi ? k + (size_t)len : hi;
  for (size_t i = lo; i <= row1; ++i) {
    double s = 0.0;
    for (int j = 0; j < len; ++j) {
      s += h[i * n + k + (size_t)j] * u[j];
    }
    const double f = 2.0 * s / uu;
    for (int j = 0; j < len; ++j) {
      h[i * n + k + (size_t)j] -= f * u[j];
    }
  }
}

static void householder3(double x, double y, double z, double* u) {
  const double norm = sqrt(x * x + y * y + z * z);
  u[0] = x + (x >= 0.0 ? norm : -norm);
  u[1] = y;
  u[2] = z;
}

/* The eigenvalues of A of any order; work: n doubles for the reflectors of the Hessenberg
 * reduction. */
static myb_status eigenvalues(size_t n, double* a, double* re, double* im, double* work) {
  for (size_t i = 0; i < n * n; ++i) {
    if (!isfinite(a[i])) {
      return MYB_E_DOMAIN;
    }
  }
  balance(n, a);
  hessenberg(n, a, work);
  double* h = a;
  size_t hi = n - 1;
  int iterations = 0;
  const int budget = 30 * (int)n;
  /* bound: every iteration either deflates or counts against 30 n sweeps */
  while (1) {
    if (hi == 0) {
      re[0] = h[0];
      im[0] = 0.0;
      break;
    }
    /* find the lowest negligible subdiagonal element in [0, hi] */
    size_t lo = hi;
    /* bound: lo decreases to 0 */
    while (lo > 0) {
      const double s = fabs(h[(lo - 1) * n + lo - 1]) + fabs(h[lo * n + lo]);
      if (fabs(h[lo * n + lo - 1]) <= eps * (s > 0.0 ? s : 1.0)) {
        h[lo * n + lo - 1] = 0.0;
        break;
      }
      --lo;
    }
    if (lo == hi) { /* one eigenvalue */
      re[hi] = h[hi * n + hi];
      im[hi] = 0.0;
      --hi;
      iterations = 0;
      continue;
    }
    if (lo + 1 == hi) { /* a 2 x 2 block */
      double r2[2];
      double i2[2];
      eig2(h[(hi - 1) * n + hi - 1], h[(hi - 1) * n + hi], h[hi * n + hi - 1], h[hi * n + hi],
           r2, i2);
      re[hi - 1] = r2[0];
      im[hi - 1] = i2[0];
      re[hi] = r2[1];
      im[hi] = i2[1];
      if (hi < 2) {
        break;
      }
      hi -= 2;
      iterations = 0;
      continue;
    }
    if (++iterations > budget) {
      return MYB_E_NOCONV;
    }
    /* double shift from the trailing 2 x 2 block (exceptional shifts at 10 and 20) */
    double s = h[(hi - 1) * n + hi - 1] + h[hi * n + hi];
    double t = h[(hi - 1) * n + hi - 1] * h[hi * n + hi] - h[(hi - 1) * n + hi] * h[hi * n + hi - 1];
    if (iterations == 10 || iterations == 20) {
      const double w = fabs(h[hi * n + hi - 1]) + fabs(h[(hi - 1) * n + hi - 2]);
      s = 1.5 * w;
      t = w * w;
    }
    double x = h[lo * n + lo] * h[lo * n + lo] + h[lo * n + lo + 1] * h[(lo + 1) * n + lo] -
               s * h[lo * n + lo] + t;
    double y = h[(lo + 1) * n + lo] * (h[lo * n + lo] + h[(lo + 1) * n + lo + 1] - s);
    double z = h[(lo + 1) * n + lo] * h[(lo + 2) * n + lo + 1];
    for (size_t k = lo; k + 2 <= hi; ++k) {
      double u[3];
      householder3(x, y, z, u);
      reflect(n, h, k, 3, u, lo, hi);
      x = h[(k + 1) * n + k];
      y = h[(k + 2) * n + k];
      if (k + 3 <= hi) {
        z = h[(k + 3) * n + k];
      }
    }
    double u[3];
    const double norm = sqrt(x * x + y * y);
    u[0] = x + (x >= 0.0 ? norm : -norm);
    u[1] = y;
    reflect(n, h, hi - 1, 2, u, lo, hi);
  }
  return MYB_OK;
}

myb_status myb_eigenvalues(size_t n, double* a, double* re, double* im) {
  MYB_CHECK(a != NULL && re != NULL && im != NULL && n >= 1 && n <= MAX_N);
  double v[MAX_N];
  return eigenvalues(n, a, re, im, v);
}

myb_status myb_internal_eigenvalues(size_t n, double* a, double* re, double* im, double* work) {
  MYB_CHECK(a != NULL && re != NULL && im != NULL && work != NULL && n >= 1);
  return eigenvalues(n, a, re, im, work);
}

myb_status myb_eigenvector_real(size_t n, const double* a, double lambda, double* v,
                                double* work) {
  MYB_CHECK(a != NULL && v != NULL && work != NULL && n >= 1 && n <= MAX_N);
  return myb_internal_eigenvector_real(n, a, lambda, v, work);
}

myb_status myb_internal_eigenvector_real(size_t n, const double* a, double lambda, double* v,
                                         double* work) {
  MYB_CHECK(a != NULL && v != NULL && work != NULL && n >= 1);
  double norm = 0.0;
  for (size_t i = 0; i < n * n; ++i) {
    norm = fmax(norm, fabs(a[i]));
  }
  if (!isfinite(lambda) || !isfinite(norm)) {
    return MYB_E_DOMAIN;
  }
  /* (A - (lambda + delta) I) x = v, delta keeping the system just regular */
  const double delta = 1e-10 * (fabs(lambda) + norm + 1e-300);
  double* m = work; /* n x n */
  for (size_t i = 0; i < n; ++i) {
    v[i] = 1.0 / sqrt((double)n) * (i % 2 == 0 ? 1.0 : 0.9);
  }
  for (int it = 0; it < 4; ++it) {
    memcpy(m, a, n * n * sizeof(double));
    for (size_t i = 0; i < n; ++i) {
      m[i * n + i] -= lambda + delta;
    }
    const myb_status st = myb_internal_lu_solve(n, m, v, 1);
    if (st != MYB_OK) {
      return st;
    }
    double s = 0.0;
    size_t big = 0;
    for (size_t i = 0; i < n; ++i) {
      s += v[i] * v[i];
      if (fabs(v[i]) > fabs(v[big])) {
        big = i;
      }
    }
    s = sqrt(s);
    if (!(s > 0.0) || !isfinite(s)) {
      return MYB_E_NOCONV;
    }
    const double sign = v[big] < 0.0 ? -1.0 : 1.0; /* largest component positive */
    for (size_t i = 0; i < n; ++i) {
      v[i] *= sign / s;
    }
  }
  return MYB_OK;
}
