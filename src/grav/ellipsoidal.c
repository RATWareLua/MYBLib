/* SPDX-License-Identifier: BSD-3-Clause */
/* Exterior ellipsoidal harmonics (G. Romain, B. Jean-Pierre, Celest. Mech. Dyn. Astron. 79, 235,
 * 2001; R. Garmier, J.-P. Barriot, Celest. Mech. Dyn. Astron. 79, 235, 2001; E. W. Hobson,
 * The Theory of Spherical and Ellipsoidal Harmonics, 1931, ch. XI).
 *
 * Units: lengths are scaled by the semi-major axis a of the reference ellipsoid, so h^2 =
 * 1 - b^2/a^2, k^2 = 1 - c^2/a^2 and lambda >= 1 outside it.
 *
 * Lamé functions of degree n: E(x) = psi(x) x^e Q(x^2), Q monic of degree r, in four classes
 *   K: psi = 1,                       e = n mod 2,       r = floor(n / 2)
 *   L: psi = sqrt|x^2 - h^2|,         e = (n - 1) mod 2, r = floor((n - 1) / 2)
 *   M: psi = sqrt|x^2 - k^2|,         e = (n - 1) mod 2, r = floor((n - 1) / 2)
 *   N: psi = sqrt|(x^2-h^2)(x^2-k^2)|, e = n mod 2,      r = floor((n - 2) / 2)
 * (2n + 1 functions per degree). Substituting P(x) = x^e Q(x^2) = sum_j b_j x^(m_j),
 * m_j = n - o - 2j (o = 0, 1, 1, 2 for K, L, M, N), into Lamé's equation
 *   (x^2 - h^2)(x^2 - k^2) E'' + x (2x^2 - h^2 - k^2) E' + (Lambda - n(n+1) x^2) E = 0
 * gives the three-term recurrence
 *   [(m_{j+1} + o)(m_{j+1} + o + 1) - n(n+1)] b_{j+1} + (Lambda - d(m_j)) b_j
 *     + h^2 k^2 m_{j-1}(m_{j-1} - 1) b_{j-1} = 0,
 *   d(m) = (h^2+k^2) m^2 (K), h^2 m^2 + k^2 (m+1)^2 (L), k^2 m^2 + h^2 (m+1)^2 (M),
 *          (h^2+k^2)(m+1)^2 (N),
 * whose eigenvalues Lambda are real and distinct (Sturm-Liouville); they are found as the roots
 * of the continuant, the vector b by the recurrence, then the r roots t_i of Q, which lie in
 * (0, k^2) (Stieltjes).
 *
 * The product E(lambda) E(mu) E(nu) of the three ellipsoidal coordinates is a solid harmonic
 * polynomial in x, y, z (F. Niven): with x^2/s + y^2/(s - h^2) + z^2/(s - k^2) = 1 for
 * s = lambda^2, mu^2, nu^2,
 *   lambda mu nu = h k x,  prod sqrt|c^2 - h^2| = h sqrt(k^2 - h^2) |y|,
 *   prod sqrt|c^2 - k^2| = k sqrt(k^2 - h^2) |z|,
 *   prod_c (c^2 - t) = (t - h^2)(t - k^2)(x^2 - t) + y^2 t (t - k^2) + z^2 t (t - h^2),
 * so H(x, y, z) = (h k x)^ex (h sqrt(k^2-h^2) y)^ey (k sqrt(k^2-h^2) z)^ez prod_i M_i(x, y, z)
 * with the signs of x, y, z — smooth everywhere, without the singular sets of the coordinates.
 * The exterior term is F(lambda) E(mu) E(nu) = R(lambda) H(x, y, z), with
 *   R(lambda) = integral_lambda^inf du / (E(u)^2 sqrt((u^2 - h^2)(u^2 - k^2)))
 *             = lambda^(-2n-1) integral_0^1 t^(2n) dt / (Ehat(t)^2 Shat(t)),
 * Ehat(t) = psihat(t) prod_i (1 - t_i t^2 / lambda^2), Shat(t) = sqrt((1 - h^2 t^2 / lambda^2)
 * (1 - k^2 t^2 / lambda^2)) bounded and positive for lambda >= 1 (Gauss-Legendre, 3 panels);
 * R'(lambda) = -1 / (E(lambda)^2 sqrt(...)). The gradient follows from d(lambda^2)/dx_i of the
 * confocal equation. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/grav.h"
#include "num/gauss24.h"

enum { MAX_DEGREE = 20, MAX_ROOTS = 10, CLASS_K = 0, CLASS_L = 1, CLASS_M = 2, CLASS_N = 3 };

typedef struct lame_term {
  int n;
  int cls;
  int ex;
  int ey;
  int ez;
  int roots;
  double t[MAX_ROOTS]; /* roots of Q (scaled by a^2) */
  double factor;       /* (h k)^ex (h sqrt(k^2-h^2))^ey (k sqrt(k^2-h^2))^ez / E(1) */
  double i_ref;        /* the integral of R at lambda = 1 */
} lame_term;

struct myb_ellipsoidal_field {
  myb_allocator allocator;
  double gm;
  double a;
  double h2;
  double k2;
  int degree;
  size_t count;
  lame_term* terms;
  double* alpha;
};

/* ---- Lamé functions -------------------------------------------------------------------------- */

static int class_offset(int cls) {
  return cls == CLASS_K ? 0 : (cls == CLASS_N ? 2 : 1);
}

static int class_size(int n, int cls) {
  const int o = class_offset(cls);
  return n < o ? 0 : (n - o) / 2 + 1;
}

static double diag_term(int cls, double m, double h2, double k2) {
  switch (cls) {
    case CLASS_K:
      return (h2 + k2) * m * m;
    case CLASS_L:
      return h2 * m * m + k2 * (m + 1.0) * (m + 1.0);
    case CLASS_M:
      return k2 * m * m + h2 * (m + 1.0) * (m + 1.0);
    default:
      return (h2 + k2) * (m + 1.0) * (m + 1.0);
  }
}

typedef struct tridiagonal {
  int size;
  double d[MAX_ROOTS + 2];
  double up[MAX_ROOTS + 2];  /* row j, column j + 1 */
  double low[MAX_ROOTS + 2]; /* row j, column j - 1 */
} tridiagonal;

static void lame_matrix(int n, int cls, double h2, double k2, tridiagonal* t) {
  const int o = class_offset(cls);
  const int size = class_size(n, cls);
  t->size = size;
  const double nn = (double)n * (n + 1);
  for (int j = 0; j < size; ++j) {
    const double m = (double)(n - o - 2 * j);
    t->d[j] = diag_term(cls, m, h2, k2);
    const double m_up = m - 2.0;
    t->up[j] = j + 1 < size ? -((m_up + o) * (m_up + o + 1.0) - nn) : 0.0;
    const double m_low = m + 2.0;
    t->low[j] = j > 0 ? -h2 * k2 * m_low * (m_low - 1.0) : 0.0;
  }
}

/* det(T - x I) by the continuant, scaled to stay in range. */
static double continuant(const tridiagonal* t, double x) {
  double f0 = 1.0;
  double f1 = t->d[0] - x;
  for (int j = 1; j < t->size; ++j) {
    const double f2 = (t->d[j] - x) * f1 - t->up[j - 1] * t->low[j] * f0;
    f0 = f1;
    f1 = f2;
    const double big = fabs(f1) > fabs(f0) ? fabs(f1) : fabs(f0);
    if (big > 1e100) {
      f0 *= 1e-100;
      f1 *= 1e-100;
    }
  }
  return f1;
}

/* Roots of f on [lo, hi], expected count; sampling with sign changes, then bisection. */
typedef double (*scalar_function)(const void* context, double x);

static int real_roots(scalar_function f, const void* context, double lo, double hi, int expected,
                      double* roots) {
  if (expected == 0) {
    return 1;
  }
  /* bound: 8 refinements of the sampling */
  for (int level = 0; level < 8; ++level) {
    const int samples = (64 * expected) << level;
    int found = 0;
    double x0 = lo;
    double f0 = f(context, x0);
    for (int k = 1; k <= samples && found <= expected; ++k) {
      const double x1 = lo + (hi - lo) * k / samples;
      const double f1 = f(context, x1);
      if (f0 == 0.0) {
        if (found < expected) {
          roots[found] = x0;
        }
        ++found;
      } else if ((f0 < 0.0) != (f1 < 0.0) && f1 != 0.0) {
        double a = x0;
        double b = x1;
        double fa = f0;
        /* bound: bisection to the resolution of doubles, at most 200 halvings */
        for (int it = 0; it < 200 && b - a > 0.0; ++it) {
          const double mid = 0.5 * (a + b);
          if (mid <= a || mid >= b) {
            break;
          }
          const double fm = f(context, mid);
          if ((fm < 0.0) == (fa < 0.0) && fm != 0.0) {
            a = mid;
            fa = fm;
          } else {
            b = mid;
          }
        }
        if (found < expected) {
          roots[found] = 0.5 * (a + b);
        }
        ++found;
      }
      x0 = x1;
      f0 = f1;
    }
    if (found == expected) {
      return 1;
    }
  }
  return 0;
}

static double continuant_function(const void* context, double x) {
  return continuant(context, x);
}

typedef struct polynomial {
  int degree;
  double c[MAX_ROOTS + 1]; /* c[0] t^degree + ... + c[degree] */
} polynomial;

static double polynomial_function(const void* context, double x) {
  const polynomial* p = context;
  double v = 0.0;
  for (int i = 0; i <= p->degree; ++i) {
    v = v * x + p->c[i];
  }
  return v;
}

/* E(x) / psi(x)... evaluated as psihat and products: |E(x)| for x >= 0 (scaled). */
static double lame_value(const lame_term* term, double x, double h2, double k2) {
  double v = term->ex ? x : 1.0;
  for (int i = 0; i < term->roots; ++i) {
    v *= x * x - term->t[i];
  }
  if (term->cls == CLASS_L || term->cls == CLASS_N) {
    v *= sqrt(fabs(x * x - h2));
  }
  if (term->cls == CLASS_M || term->cls == CLASS_N) {
    v *= sqrt(fabs(x * x - k2));
  }
  return v;
}

/* integral_0^1 t^(2n) / (Ehat(t)^2 Shat(t)) dt at lambda (scaled), three panels. */
static double radial_integral(const lame_term* term, double lambda, double h2, double k2) {
  static const double edges[4] = {0.0, 0.5, 0.85, 1.0};
  const double inv2 = 1.0 / (lambda * lambda);
  double sum = 0.0;
  for (int panel = 0; panel < 3; ++panel) {
    const double lo = edges[panel];
    const double half = 0.5 * (edges[panel + 1] - lo);
    for (int k = 0; k < 24; ++k) {
      const double x = k < 12 ? myb_gauss24_x[k] : -myb_gauss24_x[k - 12];
      const double w = myb_gauss24_w[k < 12 ? k : k - 12];
      const double t = lo + half * (x + 1.0);
      const double t2 = t * t;
      const double ph = 1.0 - h2 * inv2 * t2;
      const double pk = 1.0 - k2 * inv2 * t2;
      double e = 1.0;
      for (int i = 0; i < term->roots; ++i) {
        e *= 1.0 - term->t[i] * inv2 * t2;
      }
      double e2 = e * e;
      if (term->cls == CLASS_L || term->cls == CLASS_N) {
        e2 *= ph;
      }
      if (term->cls == CLASS_M || term->cls == CLASS_N) {
        e2 *= pk;
      }
      double tn = 1.0;
      for (int p = 0; p < term->n; ++p) {
        tn *= t2;
      }
      sum += w * half * tn / (e2 * sqrt(ph * pk));
    }
  }
  return sum;
}

static int build_term(lame_term* term, int n, int cls, int index, double h2, double k2) {
  memset(term, 0, sizeof *term);
  term->n = n;
  term->cls = cls;
  const int o = class_offset(cls);
  term->ex = (n - o) % 2;
  term->ey = cls == CLASS_L || cls == CLASS_N;
  term->ez = cls == CLASS_M || cls == CLASS_N;
  tridiagonal t;
  lame_matrix(n, cls, h2, k2, &t);
  const int size = t.size;
  term->roots = size - 1;
  double eig[MAX_ROOTS + 2];
  /* Gershgorin bounds */
  double lo = 0.0;
  double hi = 0.0;
  for (int j = 0; j < size; ++j) {
    const double r = fabs(t.up[j]) + fabs(t.low[j]);
    lo = j == 0 || t.d[j] - r < lo ? t.d[j] - r : lo;
    hi = j == 0 || t.d[j] + r > hi ? t.d[j] + r : hi;
  }
  lo -= 1e-9 * (fabs(lo) + fabs(hi) + 1.0);
  hi += 1e-9 * (fabs(lo) + fabs(hi) + 1.0);
  if (!real_roots(continuant_function, &t, lo, hi, size, eig)) {
    return 0;
  }
  const double lambda = eig[index];
  /* eigenvector by the recurrence, b_0 = 1 */
  double b[MAX_ROOTS + 2];
  b[0] = 1.0;
  if (size > 1) {
    b[1] = -(t.d[0] - lambda) / t.up[0];
    for (int j = 1; j + 1 < size; ++j) {
      b[j + 1] = -(t.low[j] * b[j - 1] + (t.d[j] - lambda) * b[j]) / t.up[j];
    }
  }
  /* roots of Q(t) = sum b_j t^(r - j) in (0, k^2) */
  polynomial q = {.degree = term->roots};
  for (int j = 0; j <= term->roots; ++j) {
    q.c[j] = b[j];
  }
  if (!real_roots(polynomial_function, &q, 0.0, k2, term->roots, term->t)) {
    return 0;
  }
  const double hk = sqrt(h2 * k2);
  const double kh = sqrt(k2 - h2);
  double factor = 1.0;
  if (term->ex) {
    factor *= hk;
  }
  if (term->ey) {
    factor *= sqrt(h2) * kh;
  }
  if (term->ez) {
    factor *= sqrt(k2) * kh;
  }
  term->factor = factor / lame_value(term, 1.0, h2, k2);
  term->i_ref = radial_integral(term, 1.0, h2, k2);
  return 1;
}

/* ---- Field ------------------------------------------------------------------------------------ */

void myb_ellipsoidal_field_destroy(myb_ellipsoidal_field* field) {
  if (field == NULL) {
    return;
  }
  const myb_allocator a = field->allocator;
  myb_internal_release(&a, field->terms);
  myb_internal_release(&a, field->alpha);
  myb_internal_release(&a, field);
}

myb_status myb_ellipsoidal_field_create(double gm, const double axes[3], int degree,
                                        const double* alpha, const myb_allocator* allocator,
                                        myb_ellipsoidal_field** field) {
  MYB_CHECK(axes != NULL && field != NULL);
  const double a = axes[0];
  const double b = axes[1];
  const double c = axes[2];
  if (!(gm > 0.0) || !isfinite(gm) || !(a > b) || !(b > c) || !(c > 0.0) || !isfinite(a) ||
      degree < 0 || degree > MAX_DEGREE) {
    return MYB_E_DOMAIN;
  }
  const size_t count = (size_t)(degree + 1) * (size_t)(degree + 1);
  if (alpha != NULL) {
    for (size_t k = 0; k < count; ++k) {
      if (!isfinite(alpha[k])) {
        return MYB_E_DOMAIN;
      }
    }
  }
  const myb_allocator al = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&al, 1, sizeof(myb_ellipsoidal_field), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_ellipsoidal_field* f = block;
  memset(f, 0, sizeof *f);
  f->allocator = al;
  f->gm = gm;
  f->a = a;
  f->h2 = 1.0 - (b / a) * (b / a);
  f->k2 = 1.0 - (c / a) * (c / a);
  f->degree = degree;
  f->count = count;
  status = myb_internal_alloc_array(&al, count, sizeof(lame_term), &block);
  if (status == MYB_OK) {
    f->terms = block;
    status = myb_internal_alloc_array(&al, count, sizeof(double), &block);
  }
  if (status != MYB_OK) {
    myb_ellipsoidal_field_destroy(f);
    return status;
  }
  f->alpha = block;
  for (size_t k = 0; k < count; ++k) {
    f->alpha[k] = alpha != NULL ? alpha[k] : 0.0;
  }
  size_t k = 0;
  for (int n = 0; n <= degree && status == MYB_OK; ++n) {
    for (int cls = CLASS_K; cls <= CLASS_N && status == MYB_OK; ++cls) {
      const int size = class_size(n, cls);
      for (int i = 0; i < size; ++i) {
        if (!build_term(&f->terms[k], n, cls, i, f->h2, f->k2)) {
          status = MYB_E_NOCONV;
          break;
        }
        ++k;
      }
    }
  }
  if (status == MYB_OK && k != count) {
    status = MYB_E_NOCONV;
  }
  if (status != MYB_OK) {
    myb_ellipsoidal_field_destroy(f);
    return status;
  }
  *field = f;
  return MYB_OK;
}

int myb_ellipsoidal_field_degree(const myb_ellipsoidal_field* field) {
  MYB_CHECK(field != NULL);
  return field->degree;
}

/* lambda^2 (scaled) of the point: the root s > k^2 of x^2/s + y^2/(s-h^2) + z^2/(s-k^2) = 1. */
static double lambda_squared(const myb_ellipsoidal_field* f, const double x[3]) {
  const double xx = x[0] * x[0];
  const double yy = x[1] * x[1];
  const double zz = x[2] * x[2];
  const double r2 = xx + yy + zz;
  if (!(r2 > 0.0)) {
    return 0.0;
  }
  double s = r2 + f->k2; /* G(s) <= 0: at or right of the root */
  /* bound: Newton on a convex decreasing function, monotone after the first step; 100 at most */
  for (int it = 0; it < 100; ++it) {
    const double a = s;
    const double b = s - f->h2;
    const double c = s - f->k2;
    const double g = xx / a + yy / b + zz / c - 1.0;
    const double dg = -(xx / (a * a) + yy / (b * b) + zz / (c * c));
    double next = s - g / dg;
    if (!(next > f->k2)) {
      next = 0.5 * (s + f->k2);
    }
    if (next == s || (it > 1 && fabs(next - s) <= 4e-16 * s)) {
      s = next;
      break;
    }
    s = next;
  }
  return s;
}

/* The basis function of term k (potential per unit alpha, in units of GM / a) and its gradient
 * with respect to the scaled coordinates. */
static void term_value(const myb_ellipsoidal_field* f, const lame_term* term, const double x[3],
                       double lambda, const double dlambda[3], double* value, double grad[3]) {
  /* H(x) / E(1) and its gradient */
  double prod = term->factor;
  double dprod[3] = {0.0, 0.0, 0.0};
  /* product rule over the factors: x, y, z powers and the quadrics M_i */
  double factors[MAX_ROOTS + 3];
  double dfactors[MAX_ROOTS + 3][3];
  int count = 0;
  if (term->ex) {
    factors[count] = x[0];
    dfactors[count][0] = 1.0;
    dfactors[count][1] = 0.0;
    dfactors[count][2] = 0.0;
    ++count;
  }
  if (term->ey) {
    factors[count] = x[1];
    dfactors[count][0] = 0.0;
    dfactors[count][1] = 1.0;
    dfactors[count][2] = 0.0;
    ++count;
  }
  if (term->ez) {
    factors[count] = x[2];
    dfactors[count][0] = 0.0;
    dfactors[count][1] = 0.0;
    dfactors[count][2] = 1.0;
    ++count;
  }
  for (int i = 0; i < term->roots; ++i) {
    const double t = term->t[i];
    const double th = t - f->h2;
    const double tk = t - f->k2;
    factors[count] = th * tk * (x[0] * x[0] - t) + x[1] * x[1] * t * tk + x[2] * x[2] * t * th;
    dfactors[count][0] = 2.0 * x[0] * th * tk;
    dfactors[count][1] = 2.0 * x[1] * t * tk;
    dfactors[count][2] = 2.0 * x[2] * t * th;
    ++count;
  }
  /* prefix and suffix products for the gradient without division */
  double prefix[MAX_ROOTS + 4];
  double suffix[MAX_ROOTS + 4];
  prefix[0] = 1.0;
  for (int i = 0; i < count; ++i) {
    prefix[i + 1] = prefix[i] * factors[i];
  }
  suffix[count] = 1.0;
  for (int i = count - 1; i >= 0; --i) {
    suffix[i] = suffix[i + 1] * factors[i];
  }
  const double h = prod * prefix[count];
  for (int i = 0; i < count; ++i) {
    const double others = prod * prefix[i] * suffix[i + 1];
    for (int c = 0; c < 3; ++c) {
      dprod[c] += others * dfactors[i][c];
    }
  }
  /* radial factor W = R(lambda) / R(1) and W' */
  const int n = term->n;
  double inv = 1.0;
  for (int p = 0; p < 2 * n + 1; ++p) {
    inv /= lambda;
  }
  const double w = inv * radial_integral(term, lambda, f->h2, f->k2) / term->i_ref;
  const double e = lame_value(term, lambda, f->h2, f->k2);
  const double s = sqrt((lambda * lambda - f->h2) * (lambda * lambda - f->k2));
  const double dw = -1.0 / (e * e * s) / term->i_ref;
  *value = w * h;
  if (grad != NULL) {
    for (int c = 0; c < 3; ++c) {
      grad[c] = dw * dlambda[c] * h + w * dprod[c];
    }
  }
}

static myb_status coordinates(const myb_ellipsoidal_field* f, const double r[3], double x[3],
                              double* lambda, double dlambda[3]) {
  if (!isfinite(r[0]) || !isfinite(r[1]) || !isfinite(r[2])) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    x[i] = r[i] / f->a;
  }
  const double s = lambda_squared(f, x);
  if (!(s >= 1.0 - 1e-12)) {
    return MYB_E_DOMAIN; /* inside the reference ellipsoid: the series need not converge */
  }
  *lambda = sqrt(s);
  const double a = s;
  const double b = s - f->h2;
  const double c = s - f->k2;
  const double gs = -(x[0] * x[0] / (a * a) + x[1] * x[1] / (b * b) + x[2] * x[2] / (c * c));
  const double gx[3] = {2.0 * x[0] / a, 2.0 * x[1] / b, 2.0 * x[2] / c};
  for (int i = 0; i < 3; ++i) {
    dlambda[i] = -gx[i] / gs / (2.0 * *lambda);
  }
  return MYB_OK;
}

myb_status myb_ellipsoidal_field_eval(const myb_ellipsoidal_field* field, const double r[3],
                                      double* potential, double a[3]) {
  MYB_CHECK(field != NULL && r != NULL);
  double x[3];
  double lambda = 0.0;
  double dlambda[3];
  const myb_status st = coordinates(field, r, x, &lambda, dlambda);
  if (st != MYB_OK) {
    return st;
  }
  double v = 0.0;
  double g[3] = {0.0, 0.0, 0.0};
  for (size_t k = 0; k < field->count; ++k) {
    if (field->alpha[k] == 0.0) {
      continue;
    }
    double tv = 0.0;
    double tg[3];
    term_value(field, &field->terms[k], x, lambda, dlambda, &tv, a != NULL ? tg : NULL);
    v += field->alpha[k] * tv;
    if (a != NULL) {
      for (int c = 0; c < 3; ++c) {
        g[c] += field->alpha[k] * tg[c];
      }
    }
  }
  const double scale = field->gm / field->a;
  if (potential != NULL) {
    *potential = scale * v;
  }
  if (a != NULL) {
    for (int c = 0; c < 3; ++c) {
      a[c] = scale / field->a * g[c];
    }
  }
  return MYB_OK;
}

myb_status myb_ellipsoidal_field_basis(const myb_ellipsoidal_field* field, const double r[3],
                                       double* values) {
  MYB_CHECK(field != NULL && r != NULL && values != NULL);
  double x[3];
  double lambda = 0.0;
  double dlambda[3];
  const myb_status st = coordinates(field, r, x, &lambda, dlambda);
  if (st != MYB_OK) {
    return st;
  }
  for (size_t k = 0; k < field->count; ++k) {
    term_value(field, &field->terms[k], x, lambda, dlambda, &values[k], NULL);
  }
  return MYB_OK;
}

myb_status myb_ellipsoidal_field_set_coefficients(myb_ellipsoidal_field* field,
                                                  const double* alpha) {
  MYB_CHECK(field != NULL && alpha != NULL);
  for (size_t k = 0; k < field->count; ++k) {
    if (!isfinite(alpha[k])) {
      return MYB_E_DOMAIN;
    }
  }
  memcpy(field->alpha, alpha, field->count * sizeof(double));
  return MYB_OK;
}

myb_status myb_ellipsoidal_lame(const myb_ellipsoidal_field* field, int n, int p, double x,
                                double* value) {
  MYB_CHECK(field != NULL && value != NULL);
  if (n < 0 || n > field->degree || p < 0 || p > 2 * n || !isfinite(x)) {
    return MYB_E_DOMAIN;
  }
  const lame_term* term = &field->terms[(size_t)n * (size_t)n + (size_t)p];
  const double s = x / field->a;
  double v = term->ex ? s : 1.0;
  for (int i = 0; i < term->roots; ++i) {
    v *= s * s - term->t[i];
  }
  if (term->ey) {
    v *= sqrt(fabs(s * s - field->h2));
  }
  if (term->ez) {
    v *= sqrt(fabs(s * s - field->k2));
  }
  /* E(x) in m^n: the scaled function times a^n */
  double an = 1.0;
  for (int i = 0; i < n; ++i) {
    an *= field->a;
  }
  *value = v * an;
  return MYB_OK;
}
