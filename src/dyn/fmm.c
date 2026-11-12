/* SPDX-License-Identifier: BSD-3-Clause */
/* The fast multipole method for the Newtonian accelerations of point masses: the approximate
 * level of the N-body graph (myb_nbody_set_approximation, myb/dyn.h).
 *
 * Expansions. With the spherical coordinates (r, theta, phi) of a vector and the associated
 * Legendre functions P_n^m with the Condon-Shortley phase, the scaled solid harmonics
 *   R_n^m = r^n P_n^m(cos theta) e^(i m phi) / (n + m)!          (regular)
 *   I_n^m = (n - m)! P_n^m(cos theta) e^(i m phi) / r^(n + 1)     (irregular)
 * with X_n^(-m) = (-1)^m conj(X_n^m) turn the translation theorems into plain convolutions
 * (C. A. White, M. Head-Gordon, J. Chem. Phys. 101, 6593, 1994; the normalisation of the
 * harmonics Upsilon and Theta of W. Dehnen, Comput. Astrophys. Cosmol. 1, 1, 2014):
 *   1 / |r - a|  = sum_n sum_m conj(R_n^m(a)) I_n^m(r)                       (|a| < |r|)
 *   R_n^m(a + b) = sum_(k <= n) sum_l R_k^l(a) R_(n-k)^(m-l)(b)
 *   I_n^m(r - a) = sum_k sum_l conj(R_k^l(a)) I_(n+k)^(m+l)(r)               (|a| < |r|)
 * (checked numerically before use: the addition of R to 3e-18, the whole chain below to its
 * truncation error). They come from the Cartesian components by recurrences — no angle and
 * no function of libm, one square root per vector:
 *   R_0^0 = 1, R_m^m = -(x + i y) R_(m-1)^(m-1) / (2m), R_(m+1)^m = z R_m^m,
 *   R_n^m = ((2n - 1) z R_(n-1)^m - r^2 R_(n-2)^m) / ((n + m)(n - m));
 *   I_0^0 = 1/r, I_m^m = -(2m - 1)(x + i y) I_(m-1)^(m-1) / r^2,
 *   I_(m+1)^m = (2m + 1) z I_m^m / r^2,
 *   I_n^m = ((2n - 1) z I_(n-1)^m - ((n - 1)^2 - m^2) I_(n-2)^m) / r^2.
 * For sources mu_b (their GM) at x_b, a cell with the centre z has the multipoles
 *   M_n^m = sum_b mu_b R_n^m(x_b - z)                                          (P2M)
 * and far from it the potential Phi(x) = sum_b mu_b / |x - x_b| = sum conj(M_n^m) I_n^m(x - z);
 * the local expansion about w, Phi(x) = sum conj(R_n^m(x - w)) F_n^m, follows from
 *   F_n^m = sum_k sum_l (-1)^k conj(M_k^l) I_(n+k)^(m+l)(z - w)                (M2L)
 * and the expansions move by
 *   M_n^m(z') = sum_k sum_l R_k^l(z - z') M_(n-k)^(m-l)(z)                     (M2M)
 *   F_n^m(w') = sum_k sum_l conj(R_k^l(w' - w)) F_(n+k)^(m+l)(w)               (L2L)
 * The acceleration grad Phi at x = w + y is read from the terms of degree 1 moved to x:
 *   G^m = sum_k sum_l conj(R_k^l(y)) F_(1+k)^(m+l),  a = (-Re G^1, -Im G^1, Re G^0)  (L2P)
 * Order p: every sum stops at the total degree n + k <= p, so the potential is kept to degree
 * p and the force to degree p - 1 in the ratio of the cell sizes to their separation. The
 * coefficients are complex and kept for m >= 0 only. Before a translation they are unfolded
 * into a "full form" (-n <= m <= n): a plane of real parts and one of imaginary parts, one row
 * per degree n with zeros beyond |m| = n, so that the sums over l run four lanes at a time
 * without a branch (M2L and L2L are sum_k sum_l a_k^l c_(n+k)^(m+l), the same loop).
 *
 * Scales. Lengths are measured in L, a power of two at least the half side of the root box,
 * and GM in G, a power of two at least the largest GM: the scaling is exact, the scaled
 * positions within the root box are <= sqrt(3) and the scaled weights <= 1. The direct sums
 * work in metres, on the positions as given (a close pair loses no digit to the scaling). The
 * irregular harmonics grow as 1/r^(n+1), so two cells interact through their expansions only
 * if their centres are at least 2^-30 L apart: with p <= 20 the largest I_n^m stays below
 * 2^900, the coefficients far from overflow whatever the clustering.
 *
 * The tree. An octree, built breadth first without recursion: a cell with more than `leaf`
 * points is split into its occupied octants (a counting sort, stable), a cell whose points
 * all fall into one octant shrinks to it instead (so every split cell has two children or
 * more, the cells are at most 2 N - 1, and the leaf size adapts to the clustering), down to
 * 32 halvings of the root box. The cells are stored level by level: the upward pass (P2M in
 * the leaves, M2M) runs over them backwards, the downward pass (L2L, L2P) forwards. The
 * centre of a cell is the centre of its GM (the centroid of its points if it holds test
 * bodies only); its radius r is the largest distance of its points from the centre.
 *
 * The interactions: the dual tree walk of W. Dehnen (J. Comput. Phys. 179, 27, 2002; the same
 * 2014) with an explicit stack. A pair of cells A, B is taken
 *   - directly (both ways, in metres) if it has no more than `direct_pairs` pairs of points,
 *     whatever its geometry — cheaper there than the expansions, and exact;
 *   - through the expansions, both ways from one table of I_n^m(z_B - z_A), if
 *       r_A + r_B < theta |z_A - z_B|        (the acceptance criterion, theta < 1);
 *   - otherwise directly if both are leaves, or by its children: the larger cell is split.
 * A cell with itself is split into its children and their pairs. `direct_pairs` is half the
 * number of complex products of a mutual M2L (one vectorised pair of points costs about two
 * of them: 825 pairs at p = 8), `leaf` 32 points (64 above p = 10): the best of the leaf sizes
 * 8 .. 128 and of the thresholds 1/16 .. 1 of the products, measured for 10^4 bodies (within
 * the noise of a shared machine, the leaf 32 and a half beat the others). Every interaction is
 * mutual: the force of B on A and that of A on B come from the same truncated expansion, and
 * the total force between two cells reduces to sums over their multipoles that are equal and
 * opposite. The momentum sum GM_i a_i is therefore zero up to rounding, whatever p and theta.
 * The direct sums take one receiver against four sources at a time (Clang's portable vectors,
 * as in dyn/nbody.c; no FMA), the sources' reactions stored back four at a time.
 *
 * The pairs pending on the stack are the unprocessed siblings of the pairs on the current path
 * of the walk: at most 35 per step, and a path has at most 2 x 32 steps (a split goes one
 * level deeper in one of the two cells), so the stack has a fixed size.
 *
 * Error and cost. See the measurements in myb/dyn.h (myb_nbody_set_approximation) and
 * tests/unit/test_fmm.c. The evaluation is serial and its order of operations is fixed (the
 * lanes of a vector are summed in a fixed order at the end): the same inputs give the same
 * bits on every machine. */
#include "dyn/fmm.h"

#include <math.h> /* sqrt */
#include <stdckdint.h>
#include <stdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/dyn.h"

enum {
  MAX_ORDER = MYB_NBODY_FMM_MAX_ORDER,
  MAX_DEPTH = 32,                         /* halvings of the root box */
  PUSH_MAX = 36,                          /* pairs pushed by one step: 8 + 28 for a cell itself */
  STACK = PUSH_MAX * (2 * MAX_DEPTH + 2), /* see the bound in the header comment */
  LANES = 4
};

/* Expansions only between centres at least this far apart, in units of L. */
static const double min_separation = 0x1p-30;

/* Vectors never cross a function boundary by value (-Wpsabi: the ABI of a 32-byte vector
 * differs with and without AVX); loads and stores go through memcpy. */
typedef double v4d __attribute__((vector_size(32)));
typedef long long v4i __attribute__((vector_size(32)));
#define MYB_SPLAT(x) ((v4d){(x), (x), (x), (x)})
#define MYB_LOAD(v, p) memcpy(&(v), (p), sizeof(v))
#define MYB_STORE(p, v) memcpy((p), &(v), sizeof(v))
#define MYB_HSUM(v) (((v)[0] + (v)[1]) + ((v)[2] + (v)[3]))

typedef struct point {
  double x; /* m */
  double y;
  double z;
  double w;     /* GM (m^3 s^-2); 0 for a test body */
  size_t index; /* in the caller's arrays */
  int test;     /* a test body: receives only */
} point;

typedef struct cell {
  double box[3]; /* centre of the box (m) */
  double half;   /* half side of the box */
  double z[3];   /* centre of the expansions */
  double radius; /* largest distance of a point from z */
  double weight; /* sum of GM */
  size_t begin;  /* points [begin, end) */
  size_t end;
  size_t parent;
  size_t child;      /* first child: the children are consecutive */
  unsigned children; /* 0 for a leaf */
  unsigned level;    /* in the tree */
  unsigned depth;    /* halvings of the root box */
  int local;         /* the local expansion is not zero */
} cell;

typedef struct cell_pair {
  size_t a;
  size_t b;
} cell_pair;

struct myb_internal_fmm {
  myb_allocator allocator;
  size_t count;
  int order;
  double theta;
  size_t leaf;
  size_t direct_pairs;
  size_t ncoef;  /* (p + 1)(p + 2) / 2: complex coefficients with m >= 0 */
  size_t stride; /* of a row of a full form: 2p + 8 */
  size_t origin; /* column of m = 0: p + 4 */
  size_t plane;  /* (p + 1) stride: the real plane, then the imaginary one */
  point* points; /* count: sorted into the cells */
  point* spare;  /* count: the partition */
  double* near;  /* 8 count: x, y, z, GM, test mask, ax, ay, az of the sorted points */
  double* nx;
  double* ny;
  double* nz;
  double* nw;
  double* nmask; /* all bits set for a massive point, zero for a test body */
  double* ax;    /* accelerations in the sorted order (m s^-2) */
  double* ay;
  double* az;
  cell* cells;
  size_t ncells;
  size_t capacity;   /* of cells, multipole and local */
  double* multipole; /* 2 ncoef per cell */
  double* local;     /* 2 ncoef per cell */
  cell_pair* stack;  /* STACK */
  double* scratch;   /* table, inv and the full forms in one block */
  double* table;     /* 2 ncoef: the harmonics of one vector */
  double* inv;       /* ncoef: 1 / ((n + m)(n - m)) of the regular recurrence */
  double* full_a;    /* 2 plane each */
  double* full_b;
  double* full_c;
  double inv_length; /* 1 / L */
  double inv_weight; /* 1 / G */
  double far_scale;  /* G / L^2: from the scaled field to m s^-2 */
  myb_internal_fmm_stats stats;
};

static int finite(double x) {
  return x - x == 0.0;
}

static inline size_t tri(int n, int m) {
  return (size_t)(n * (n + 1) / 2 + m);
}

/* The real part of row n of a full form at m = 0 (the imaginary part is f->plane further). */
static inline double* row(const myb_internal_fmm* f, double* full, int n) {
  return full + (size_t)n * f->stride + f->origin;
}

static inline const double* row_c(const myb_internal_fmm* f, const double* full, int n) {
  return full + (size_t)n * f->stride + f->origin;
}

/* The least power of two >= x for x > 0 (2^1023 for x above it, 2^-1022 below that). */
static double power_of_two_at_least(double x) {
  if (x >= 0x1p1023) {
    return 0x1p1023;
  }
  if (x <= 0x1p-1022) {
    return 0x1p-1022;
  }
  uint64_t bits = 0;
  memcpy(&bits, &x, sizeof bits);
  bits &= UINT64_C(0x7ff0000000000000);
  double p = 0.0;
  memcpy(&p, &bits, sizeof p);
  return p < x ? 2.0 * p : p;
}

/* ---- Harmonics and translations ----------------------------------------------------------- */

/* R_n^m(x, y, z) for 0 <= m <= n <= p into t (2 tri(p + 1, 0) values). */
static void regular(const myb_internal_fmm* f, int p, double x, double y, double z, double* t) {
  const double r2 = x * x + y * y + z * z;
  double re = 1.0;
  double im = 0.0;
  for (int m = 0; m <= p; ++m) {
    if (m > 0) {
      const double h = -0.5 / (double)m;
      const double re1 = h * (x * re - y * im);
      const double im1 = h * (x * im + y * re);
      re = re1;
      im = im1;
    }
    double* tm = t + 2 * tri(m, m);
    tm[0] = re;
    tm[1] = im;
    if (m < p) {
      double* t1 = t + 2 * tri(m + 1, m);
      t1[0] = z * re;
      t1[1] = z * im;
    }
    for (int n = m + 2; n <= p; ++n) {
      const size_t k = tri(n, m);
      const double* a = t + 2 * tri(n - 1, m);
      const double* b = t + 2 * tri(n - 2, m);
      const double c = (double)(2 * n - 1) * z;
      t[2 * k] = (c * a[0] - r2 * b[0]) * f->inv[k];
      t[2 * k + 1] = (c * a[1] - r2 * b[1]) * f->inv[k];
    }
  }
}

/* I_n^m(x, y, z) for 0 <= m <= n <= p into t. */
static void irregular(int p, double x, double y, double z, double* t) {
  const double r2 = x * x + y * y + z * z;
  const double s = 1.0 / r2;
  double re = 1.0 / sqrt(r2);
  double im = 0.0;
  for (int m = 0; m <= p; ++m) {
    if (m > 0) {
      const double h = -(double)(2 * m - 1) * s;
      const double re1 = h * (x * re - y * im);
      const double im1 = h * (x * im + y * re);
      re = re1;
      im = im1;
    }
    double* tm = t + 2 * tri(m, m);
    tm[0] = re;
    tm[1] = im;
    if (m < p) {
      const double c = (double)(2 * m + 1) * z * s;
      double* t1 = t + 2 * tri(m + 1, m);
      t1[0] = c * re;
      t1[1] = c * im;
    }
    for (int n = m + 2; n <= p; ++n) {
      const size_t k = tri(n, m);
      const double* a = t + 2 * tri(n - 1, m);
      const double* b = t + 2 * tri(n - 2, m);
      const double c1 = (double)(2 * n - 1) * z;
      const double c2 = (double)((n - 1) * (n - 1) - m * m);
      t[2 * k] = (c1 * a[0] - c2 * b[0]) * s;
      t[2 * k + 1] = (c1 * a[1] - c2 * b[1]) * s;
    }
  }
}

/* The full form out_n^m, -n <= m <= n <= p, of coefficients t kept for m >= 0 (X_n^(-m) =
 * (-1)^m conj(X_n^m)): s_n X_n^m, or s_n conj(X_n^m) with `conjugate`; s_n = (-1)^n with
 * `alternate`, 1 otherwise. Columns |m| > n are never written: they stay zero. */
static void expand(const myb_internal_fmm* f, int p, const double* t, int conjugate, int alternate,
                   double* out) {
  for (int n = 0; n <= p; ++n) {
    const double sn = (alternate && (n & 1) != 0) ? -1.0 : 1.0;
    const double* x = t + 2 * tri(n, 0);
    double* re = row(f, out, n);
    double* im = re + f->plane;
    for (int m = 0; m <= n; ++m) {
      const double xr = x[2 * m];
      const double xi = x[2 * m + 1];
      re[m] = sn * xr;
      im[m] = conjugate ? -(sn * xi) : sn * xi;
      if (m > 0) {
        const double sm = (m & 1) != 0 ? -sn : sn;
        re[-m] = sm * xr;
        im[-m] = conjugate ? sm * xi : -(sm * xi);
      }
    }
  }
}

/* out_n^m += s_n sum_(k <= p - n) sum_(l = -k .. k) a_k^l c_(n+k)^(m+l), 0 <= m <= n <= p
 * (out kept for m >= 0; a, c full forms; s_n = (-1)^n with `alternate`): M2L and L2L. The
 * lanes cover l = -k .. k + 3, where a is zero beyond k. */
static void convolve(const myb_internal_fmm* f, int p, const double* a, const double* c,
                     int alternate, double* out) {
  const size_t plane = f->plane;
  for (int n = 0; n <= p; ++n) {
    const double sign = (alternate && (n & 1) != 0) ? -1.0 : 1.0;
    for (int m = 0; m <= n; ++m) {
      v4d re = {};
      v4d im = {};
      for (int k = 0; k <= p - n; ++k) {
        const double* ar = row_c(f, a, k);
        const double* ai = ar + plane;
        const double* cr = row_c(f, c, n + k) + m;
        const double* ci = cr + plane;
        /* bound: (2k + 4) / 4 blocks */
        for (int l = -k; l <= k; l += LANES) {
          v4d xr;
          v4d xi;
          v4d yr;
          v4d yi;
          MYB_LOAD(xr, ar + l);
          MYB_LOAD(xi, ai + l);
          MYB_LOAD(yr, cr + l);
          MYB_LOAD(yi, ci + l);
          re += xr * yr - xi * yi;
          im += xr * yi + xi * yr;
        }
      }
      double* o = out + 2 * tri(n, m);
      o[0] += sign * MYB_HSUM(re);
      o[1] += sign * MYB_HSUM(im);
    }
  }
}

/* out_n^m += sum_(k <= n) sum_l r_k^l x_(n-k)^(m-l), |m - l| <= n - k (r, x full forms): M2M. */
static void shift_up(const myb_internal_fmm* f, int p, const double* r, const double* x,
                     double* out) {
  const size_t plane = f->plane;
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      double re = 0.0;
      double im = 0.0;
      for (int k = 0; k <= n; ++k) {
        const int j = n - k;
        const int lo = -k > m - j ? -k : m - j;
        const int hi = k < m + j ? k : m + j;
        const double* rr = row_c(f, r, k);
        const double* ri = rr + plane;
        const double* xr = row_c(f, x, j) + m; /* x_j^(m - l) = xr[-l] */
        const double* xi = xr + plane;
        for (int l = lo; l <= hi; ++l) {
          re += rr[l] * xr[-l] - ri[l] * xi[-l];
          im += rr[l] * xi[-l] + ri[l] * xr[-l];
        }
      }
      double* o = out + 2 * tri(n, m);
      o[0] += re;
      o[1] += im;
    }
  }
}

static double* multipole_of(const myb_internal_fmm* f, size_t c) {
  return f->multipole + 2 * f->ncoef * c;
}

static double* local_of(const myb_internal_fmm* f, size_t c) {
  return f->local + 2 * f->ncoef * c;
}

/* ---- Memory ------------------------------------------------------------------------------- */

/* Room for `needed` cells (never more than 2 count are needed: see the tree above). */
static myb_status reserve(myb_internal_fmm* f, size_t needed) {
  if (needed <= f->capacity) {
    return MYB_OK;
  }
  size_t most = 0;
  if (ckd_mul(&most, f->count, (size_t)2)) {
    return MYB_E_OVERFLOW;
  }
  MYB_CHECK(needed <= most);
  size_t capacity = f->capacity;
  if (ckd_mul(&capacity, capacity, (size_t)2) || capacity < needed) {
    capacity = needed;
  }
  if (capacity > most) {
    capacity = most;
  }
  size_t coefficients = 0;
  if (ckd_mul(&coefficients, capacity, 4 * f->ncoef)) {
    return MYB_E_OVERFLOW;
  }
  void* cells = NULL;
  myb_status status = myb_internal_alloc_array(&f->allocator, capacity, sizeof(cell), &cells);
  if (status != MYB_OK) {
    return status;
  }
  void* coef = NULL;
  status = myb_internal_alloc_array(&f->allocator, coefficients, sizeof(double), &coef);
  if (status != MYB_OK) {
    myb_internal_release(&f->allocator, cells);
    return status;
  }
  if (f->ncells != 0) {
    memcpy(cells, f->cells, f->ncells * sizeof(cell));
  }
  myb_internal_release(&f->allocator, f->cells);
  myb_internal_release(&f->allocator, f->multipole);
  f->cells = cells;
  f->multipole = coef;
  f->local = f->multipole + 2 * f->ncoef * capacity;
  f->capacity = capacity;
  return MYB_OK;
}

/* The mutual M2L of a pair of cells: 2 sum_(n <= p) (n + 1)(p - n + 1)^2 complex products. */
static size_t interaction_cost(int p) {
  size_t cost = 0;
  for (int n = 0; n <= p; ++n) {
    cost += (size_t)((n + 1) * (p - n + 1) * (p - n + 1));
  }
  return 2 * cost;
}

static size_t default_leaf(int p) {
  return p <= 10 ? 32 : 64;
}

void myb_internal_fmm_tune(myb_internal_fmm* fmm, size_t leaf_size, size_t direct_pairs) {
  MYB_CHECK(fmm != NULL);
  fmm->leaf = leaf_size != 0 ? leaf_size : default_leaf(fmm->order);
  fmm->direct_pairs = direct_pairs != 0 ? direct_pairs : interaction_cost(fmm->order) / 2;
  if (fmm->direct_pairs > SIZE_MAX / 4) {
    fmm->direct_pairs = SIZE_MAX / 4; /* the walk doubles it */
  }
}

void myb_internal_fmm_destroy(myb_internal_fmm* fmm) {
  if (fmm == NULL) {
    return;
  }
  const myb_allocator a = fmm->allocator;
  myb_internal_release(&a, fmm->scratch);
  myb_internal_release(&a, fmm->stack);
  myb_internal_release(&a, fmm->multipole);
  myb_internal_release(&a, fmm->cells);
  myb_internal_release(&a, fmm->near);
  myb_internal_release(&a, fmm->spare);
  myb_internal_release(&a, fmm->points);
  myb_internal_release(&a, fmm);
}

myb_status myb_internal_fmm_create(size_t count, int order, double theta, size_t leaf_size,
                                   const myb_allocator* allocator, myb_internal_fmm** fmm) {
  MYB_CHECK(fmm != NULL && count >= 1 && order >= 1 && order <= MAX_ORDER && theta > 0.0 &&
            theta < 1.0);
  const myb_allocator al = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&al, 1, sizeof(myb_internal_fmm), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_internal_fmm* f = block;
  const size_t p1 = (size_t)order + 1;
  const size_t stride = 2 * (size_t)order + 8;
  *f = (myb_internal_fmm){.allocator = al,
                          .count = count,
                          .order = order,
                          .theta = theta,
                          .ncoef = p1 * (p1 + 1) / 2,
                          .stride = stride,
                          .origin = (size_t)order + 4,
                          .plane = p1 * stride};
  myb_internal_fmm_tune(f, leaf_size, 0);
  size_t eight = 0;
  if (ckd_mul(&eight, count, (size_t)8)) {
    status = MYB_E_OVERFLOW;
    goto fail;
  }
  status = myb_internal_alloc_array(&al, count, sizeof(point), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  f->points = block;
  status = myb_internal_alloc_array(&al, count, sizeof(point), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  f->spare = block;
  status = myb_internal_alloc_array(&al, eight, sizeof(double), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  f->near = block;
  f->nx = f->near;
  f->ny = f->nx + count;
  f->nz = f->ny + count;
  f->nw = f->nz + count;
  f->nmask = f->nw + count;
  f->ax = f->nmask + count;
  f->ay = f->ax + count;
  f->az = f->ay + count;
  status = myb_internal_alloc_array(&al, STACK, sizeof(cell_pair), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  f->stack = block;
  const size_t scratch = 3 * f->ncoef + 6 * f->plane;
  status = myb_internal_alloc_array(&al, scratch, sizeof(double), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  f->scratch = block;
  memset(f->scratch, 0, scratch * sizeof(double)); /* the columns beyond |m| = n stay zero */
  f->table = f->scratch;
  f->inv = f->table + 2 * f->ncoef;
  f->full_a = f->inv + f->ncoef;
  f->full_b = f->full_a + 2 * f->plane;
  f->full_c = f->full_b + 2 * f->plane;
  for (int n = 0; n <= order; ++n) {
    for (int m = 0; m <= n; ++m) {
      f->inv[tri(n, m)] = n >= m + 2 ? 1.0 / ((double)(n + m) * (double)(n - m)) : 0.0;
    }
  }
  const size_t guess = 16 + 2 * (count / f->leaf);
  status = reserve(f, guess < 2 * count ? guess : 2 * count);
  if (status != MYB_OK) {
    goto fail;
  }
  *fmm = f;
  return MYB_OK;

fail:
  myb_internal_fmm_destroy(f);
  return status;
}

int myb_internal_fmm_order(const myb_internal_fmm* fmm) {
  MYB_CHECK(fmm != NULL);
  return fmm->order;
}

double myb_internal_fmm_theta(const myb_internal_fmm* fmm) {
  MYB_CHECK(fmm != NULL);
  return fmm->theta;
}

myb_internal_fmm_stats myb_internal_fmm_statistics(const myb_internal_fmm* fmm) {
  MYB_CHECK(fmm != NULL);
  return fmm->stats;
}

/* ---- The tree ----------------------------------------------------------------------------- */

static unsigned octant(const point* p, const double box[3]) {
  return (unsigned)(p->x >= box[0]) | (unsigned)(p->y >= box[1]) << 1 |
         (unsigned)(p->z >= box[2]) << 2;
}

static myb_status build(myb_internal_fmm* f) {
  const size_t n = f->count;
  point* pts = f->points;
  double lo[3] = {pts[0].x, pts[0].y, pts[0].z};
  double hi[3] = {lo[0], lo[1], lo[2]};
  for (size_t i = 1; i < n; ++i) {
    const double c[3] = {pts[i].x, pts[i].y, pts[i].z};
    for (int k = 0; k < 3; ++k) {
      lo[k] = c[k] < lo[k] ? c[k] : lo[k];
      hi[k] = c[k] > hi[k] ? c[k] : hi[k];
    }
  }
  double half = 0.0;
  double centre[3];
  for (int k = 0; k < 3; ++k) {
    const double h = 0.5 * (hi[k] - lo[k]);
    centre[k] = lo[k] + h;
    half = h > half ? h : half;
  }
  if (!finite(half)) {
    return MYB_E_DOMAIN;
  }
  /* a little more than the extent, so that the rounded centre leaves no point outside */
  half = half * (1.0 + 0x1p-20);
  if (!(half >= 0x1p-1000)) {
    half = 0x1p-1000;
  }
  f->inv_length = 1.0 / power_of_two_at_least(half);
  f->cells[0] = (cell){.box = {centre[0], centre[1], centre[2]}, .half = half, .end = n};
  f->ncells = 1;
  /* bound: the cells are at most 2 n - 1 (every split cell has two children or more) */
  for (size_t c = 0; c < f->ncells; ++c) {
    const size_t begin = f->cells[c].begin;
    const size_t end = f->cells[c].end;
    if (end - begin <= f->leaf) {
      continue;
    }
    double box[3] = {f->cells[c].box[0], f->cells[c].box[1], f->cells[c].box[2]};
    double h = f->cells[c].half;
    unsigned depth = f->cells[c].depth;
    size_t counts[8] = {};
    unsigned occupied = 0;
    /* bound: depth grows by one per pass, to MAX_DEPTH */
    while (depth < MAX_DEPTH) {
      memset(counts, 0, sizeof counts);
      for (size_t i = begin; i < end; ++i) {
        counts[octant(&pts[i], box)]++;
      }
      occupied = 0;
      unsigned only = 0;
      for (unsigned o = 0; o < 8; ++o) {
        if (counts[o] != 0) {
          occupied++;
          only = o;
        }
      }
      if (occupied > 1) {
        break;
      }
      /* all points in one octant: the cell shrinks to it */
      h *= 0.5;
      for (int k = 0; k < 3; ++k) {
        box[k] += ((only >> k) & 1u) != 0 ? h : -h;
      }
      depth++;
      occupied = 0;
    }
    cell* cl = &f->cells[c];
    cl->box[0] = box[0];
    cl->box[1] = box[1];
    cl->box[2] = box[2];
    cl->half = h;
    cl->depth = depth;
    if (occupied <= 1) {
      continue; /* a leaf at the depth limit */
    }
    const myb_status status = reserve(f, f->ncells + occupied);
    if (status != MYB_OK) {
      return status;
    }
    size_t offset[8];
    size_t at = begin;
    for (unsigned o = 0; o < 8; ++o) {
      offset[o] = at;
      at += counts[o];
    }
    for (size_t i = begin; i < end; ++i) {
      f->spare[offset[octant(&pts[i], box)]++] = pts[i];
    }
    memcpy(&pts[begin], &f->spare[begin], (end - begin) * sizeof(point));
    const size_t first = f->ncells;
    const unsigned level = f->cells[c].level + 1;
    size_t start = begin;
    const double q = 0.5 * h;
    for (unsigned o = 0; o < 8; ++o) {
      if (counts[o] == 0) {
        continue;
      }
      f->cells[f->ncells++] =
          (cell){.box = {box[0] + ((o & 1u) != 0 ? q : -q), box[1] + ((o & 2u) != 0 ? q : -q),
                         box[2] + ((o & 4u) != 0 ? q : -q)},
                 .half = q,
                 .begin = start,
                 .end = start + counts[o],
                 .parent = c,
                 .level = level,
                 .depth = depth + 1};
      start += counts[o];
    }
    f->cells[c].child = first;
    f->cells[c].children = occupied;
  }
  return MYB_OK;
}

/* Centres, radii, weights and multipoles, from the leaves up. */
static void upward(myb_internal_fmm* f) {
  const int p = f->order;
  const double s = f->inv_length;
  /* bound: ncells iterations */
  for (size_t c = f->ncells; c-- > 0;) {
    cell* cl = &f->cells[c];
    double* mc = multipole_of(f, c);
    memset(mc, 0, 2 * f->ncoef * sizeof(double));
    const double* b = cl->box;
    double weight = 0.0;
    double d[3] = {0.0, 0.0, 0.0};
    if (cl->children == 0) {
      for (size_t i = cl->begin; i < cl->end; ++i) {
        const point* q = &f->points[i];
        weight += q->w;
        d[0] += q->w * (q->x - b[0]);
        d[1] += q->w * (q->y - b[1]);
        d[2] += q->w * (q->z - b[2]);
      }
      if (!(weight > 0.0)) {
        for (size_t i = cl->begin; i < cl->end; ++i) {
          d[0] += f->points[i].x - b[0];
          d[1] += f->points[i].y - b[1];
          d[2] += f->points[i].z - b[2];
        }
        weight = (double)(cl->end - cl->begin);
      }
    } else {
      for (unsigned k = 0; k < cl->children; ++k) {
        const cell* ch = &f->cells[cl->child + k];
        weight += ch->weight;
        d[0] += ch->weight * (ch->z[0] - b[0]);
        d[1] += ch->weight * (ch->z[1] - b[1]);
        d[2] += ch->weight * (ch->z[2] - b[2]);
      }
      if (!(weight > 0.0)) {
        for (unsigned k = 0; k < cl->children; ++k) {
          const cell* ch = &f->cells[cl->child + k];
          const double nc = (double)(ch->end - ch->begin);
          d[0] += nc * (ch->z[0] - b[0]);
          d[1] += nc * (ch->z[1] - b[1]);
          d[2] += nc * (ch->z[2] - b[2]);
        }
        weight = (double)(cl->end - cl->begin);
      }
    }
    for (int k = 0; k < 3; ++k) {
      cl->z[k] = b[k] + d[k] / weight;
    }
    double r2 = 0.0;
    double mass = 0.0;
    for (size_t i = cl->begin; i < cl->end; ++i) {
      const point* q = &f->points[i];
      const double dx = q->x - cl->z[0];
      const double dy = q->y - cl->z[1];
      const double dz = q->z - cl->z[2];
      const double e = dx * dx + dy * dy + dz * dz;
      r2 = e > r2 ? e : r2;
      mass += q->w;
    }
    cl->radius = sqrt(r2);
    cl->weight = mass;
    cl->local = 0;
    if (!(mass > 0.0)) {
      continue;
    }
    if (cl->children == 0) {
      for (size_t i = cl->begin; i < cl->end; ++i) {
        const point* q = &f->points[i];
        if (q->w == 0.0) {
          continue;
        }
        regular(f, p, (q->x - cl->z[0]) * s, (q->y - cl->z[1]) * s, (q->z - cl->z[2]) * s,
                f->table);
        const double w = q->w * f->inv_weight;
        for (size_t k = 0; k < 2 * f->ncoef; ++k) {
          mc[k] += w * f->table[k];
        }
      }
    } else {
      for (unsigned k = 0; k < cl->children; ++k) {
        const size_t ci = cl->child + k;
        const cell* ch = &f->cells[ci];
        if (!(ch->weight > 0.0)) {
          continue;
        }
        regular(f, p, (ch->z[0] - cl->z[0]) * s, (ch->z[1] - cl->z[1]) * s,
                (ch->z[2] - cl->z[2]) * s, f->table);
        expand(f, p, f->table, 0, 0, f->full_a);
        expand(f, p, multipole_of(f, ci), 0, 0, f->full_b);
        shift_up(f, p, f->full_a, f->full_b, mc);
      }
    }
  }
}

/* ---- The interactions --------------------------------------------------------------------- */

/* Receiver i (scalar) against the sources [j0, j1): the vector part, j0 .. j0 + 4k. Returns
 * the first source left. The terms of test-test pairs are masked out (they could be 0 * inf if
 * the two bodies coincide). */
static size_t direct_lanes(myb_internal_fmm* f, size_t i, size_t j0, size_t j1, double out[3]) {
  const v4d xi = MYB_SPLAT(f->nx[i]);
  const v4d yi = MYB_SPLAT(f->ny[i]);
  const v4d zi = MYB_SPLAT(f->nz[i]);
  const v4d wi = MYB_SPLAT(f->nw[i]);
  const int test = f->points[i].test;
  v4d sx = {};
  v4d sy = {};
  v4d sz = {};
  size_t j = j0;
  /* bound: (j1 - j0) / 4 blocks */
  for (; j + LANES <= j1; j += LANES) {
    v4d xj;
    v4d yj;
    v4d zj;
    v4d wj;
    MYB_LOAD(xj, f->nx + j);
    MYB_LOAD(yj, f->ny + j);
    MYB_LOAD(zj, f->nz + j);
    MYB_LOAD(wj, f->nw + j);
    const v4d dx = xj - xi;
    const v4d dy = yj - yi;
    const v4d dz = zj - zi;
    const v4d r2 = dx * dx + dy * dy + dz * dz;
    const v4d q = MYB_SPLAT(1.0) / (r2 * __builtin_elementwise_sqrt(r2));
    const v4d sj = wj * q;
    const v4d si = wi * q;
    v4d t[3] = {sj * dx, sj * dy, sj * dz};
    v4d u[3] = {si * dx, si * dy, si * dz};
    if (test) {
      v4d mask;
      MYB_LOAD(mask, f->nmask + j);
      for (int k = 0; k < 3; ++k) {
        t[k] = (v4d)((v4i)t[k] & (v4i)mask);
        u[k] = (v4d)((v4i)u[k] & (v4i)mask);
      }
    }
    sx += t[0];
    sy += t[1];
    sz += t[2];
    v4d ax;
    v4d ay;
    v4d az;
    MYB_LOAD(ax, f->ax + j);
    MYB_LOAD(ay, f->ay + j);
    MYB_LOAD(az, f->az + j);
    ax -= u[0];
    ay -= u[1];
    az -= u[2];
    MYB_STORE(f->ax + j, ax);
    MYB_STORE(f->ay + j, ay);
    MYB_STORE(f->az + j, az);
  }
  out[0] = MYB_HSUM(sx);
  out[1] = MYB_HSUM(sy);
  out[2] = MYB_HSUM(sz);
  return j;
}

/* Points [i0, i1) with points [j0, j1), or the points [i0, i1) among themselves (`self`):
 * each pair once, both ways, in metres. */
static void direct(myb_internal_fmm* f, size_t i0, size_t i1, size_t j0, size_t j1, int self) {
  for (size_t i = i0; i < i1; ++i) {
    const size_t first = self ? i + 1 : j0;
    double s[3];
    size_t j = direct_lanes(f, i, first, j1, s);
    const double xi = f->nx[i];
    const double yi = f->ny[i];
    const double zi = f->nz[i];
    const double wi = f->nw[i];
    const int ti = f->points[i].test;
    for (; j < j1; ++j) {
      if (ti && f->points[j].test) {
        continue; /* two test bodies do not interact */
      }
      const double dx = f->nx[j] - xi;
      const double dy = f->ny[j] - yi;
      const double dz = f->nz[j] - zi;
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double q = 1.0 / (r2 * sqrt(r2));
      const double sj = f->nw[j] * q;
      const double si = wi * q;
      s[0] += sj * dx;
      s[1] += sj * dy;
      s[2] += sj * dz;
      f->ax[j] -= si * dx;
      f->ay[j] -= si * dy;
      f->az[j] -= si * dz;
    }
    f->ax[i] += s[0];
    f->ay[i] += s[1];
    f->az[i] += s[2];
  }
}

/* The local expansions of a and b from the multipoles of b and a: one table of I_n^m. */
static void interact(myb_internal_fmm* f, size_t a, size_t b) {
  const int p = f->order;
  const double s = f->inv_length;
  cell* ca = &f->cells[a];
  cell* cb = &f->cells[b];
  irregular(p, (cb->z[0] - ca->z[0]) * s, (cb->z[1] - ca->z[1]) * s, (cb->z[2] - ca->z[2]) * s,
            f->table);
  expand(f, p, f->table, 0, 0, f->full_c); /* I_n^m(z_b - z_a) */
  if (cb->weight > 0.0) {
    /* F_a = sum (-1)^k conj(M_b) I(z_b - z_a) */
    expand(f, p, multipole_of(f, b), 1, 1, f->full_a);
    convolve(f, p, f->full_a, f->full_c, 0, local_of(f, a));
    ca->local = 1;
  }
  if (ca->weight > 0.0) {
    /* F_b = sum (-1)^k conj(M_a) I(z_a - z_b) = (-1)^n sum conj(M_a) I(z_b - z_a) */
    expand(f, p, multipole_of(f, a), 1, 0, f->full_b);
    convolve(f, p, f->full_b, f->full_c, 1, local_of(f, b));
    cb->local = 1;
  }
}

/* na nb <= limit, without overflow. */
static int few_pairs(size_t na, size_t nb, size_t limit) {
  return na == 0 || nb <= limit / na;
}

static void walk(myb_internal_fmm* f) {
  cell_pair* stack = f->stack;
  size_t top = 0;
  stack[top++] = (cell_pair){0, 0};
  /* bound: a pair taken is either resolved or replaced by pairs one level deeper in one of
   * its cells; the tree is finite, so the loop ends, and the stack holds at most STACK pairs
   * (the header comment) */
  while (top > 0) {
    const cell_pair q = stack[--top];
    const cell* a = &f->cells[q.a];
    const size_t na = a->end - a->begin;
    if (q.a == q.b) {
      if (a->children == 0 || few_pairs(na, na - 1, 2 * f->direct_pairs)) {
        direct(f, a->begin, a->end, a->begin, a->end, 1);
        f->stats.p2p += na * (na - 1) / 2;
        continue;
      }
      MYB_CHECK(top + PUSH_MAX <= STACK);
      for (unsigned i = 0; i < a->children; ++i) {
        stack[top++] = (cell_pair){a->child + i, a->child + i};
        for (unsigned j = i + 1; j < a->children; ++j) {
          stack[top++] = (cell_pair){a->child + i, a->child + j};
        }
      }
      continue;
    }
    const cell* b = &f->cells[q.b];
    if (!(a->weight > 0.0) && !(b->weight > 0.0)) {
      continue; /* test bodies only */
    }
    const size_t nb = b->end - b->begin;
    if (few_pairs(na, nb, f->direct_pairs)) {
      direct(f, a->begin, a->end, b->begin, b->end, 0);
      f->stats.p2p += na * nb;
      continue;
    }
    const double dx = b->z[0] - a->z[0];
    const double dy = b->z[1] - a->z[1];
    const double dz = b->z[2] - a->z[2];
    const double r = sqrt(dx * dx + dy * dy + dz * dz);
    if (a->radius + b->radius < f->theta * r && r * f->inv_length >= min_separation) {
      interact(f, q.a, q.b);
      f->stats.m2l++;
      continue;
    }
    if (a->children == 0 && b->children == 0) {
      direct(f, a->begin, a->end, b->begin, b->end, 0);
      f->stats.p2p += na * nb;
      continue;
    }
    const int split_a = b->children == 0 || (a->children != 0 && a->radius >= b->radius);
    const cell* split = split_a ? a : b;
    const size_t other = split_a ? q.b : q.a;
    MYB_CHECK(top + PUSH_MAX <= STACK);
    for (unsigned i = 0; i < split->children; ++i) {
      stack[top++] = (cell_pair){split->child + i, other};
    }
  }
}

/* grad Phi (scaled) at x (scaled, relative to the centre) of a leaf whose local expansion is
 * the full form `full`. */
static void gradient(myb_internal_fmm* f, const double* full, double x, double y, double z,
                     double g[3]) {
  const int q = f->order - 1;
  const size_t plane = f->plane;
  regular(f, q, x, y, z, f->table);
  expand(f, q, f->table, 1, 0, f->full_a);
  v4d g0 = {};
  v4d g1r = {};
  v4d g1i = {};
  for (int k = 0; k <= q; ++k) {
    const double* ar = row_c(f, f->full_a, k);
    const double* ai = ar + plane;
    const double* c0r = row_c(f, full, k + 1);
    const double* c0i = c0r + plane;
    const double* c1r = c0r + 1;
    const double* c1i = c0i + 1;
    /* bound: (2k + 4) / 4 blocks */
    for (int l = -k; l <= k; l += LANES) {
      v4d xr;
      v4d xi;
      v4d y0r;
      v4d y0i;
      v4d y1r;
      v4d y1i;
      MYB_LOAD(xr, ar + l);
      MYB_LOAD(xi, ai + l);
      MYB_LOAD(y0r, c0r + l);
      MYB_LOAD(y0i, c0i + l);
      MYB_LOAD(y1r, c1r + l);
      MYB_LOAD(y1i, c1i + l);
      g0 += xr * y0r - xi * y0i;
      g1r += xr * y1r - xi * y1i;
      g1i += xr * y1i + xi * y1r;
    }
  }
  g[0] = -MYB_HSUM(g1r);
  g[1] = -MYB_HSUM(g1i);
  g[2] = MYB_HSUM(g0);
}

/* L2L from the root down, L2P in the leaves. */
static void downward(myb_internal_fmm* f) {
  const int p = f->order;
  const double s = f->inv_length;
  /* bound: ncells iterations; a parent comes before its children */
  for (size_t c = 0; c < f->ncells; ++c) {
    cell* cl = &f->cells[c];
    double* fc = local_of(f, c);
    if (c != 0 && f->cells[cl->parent].local) {
      const cell* pa = &f->cells[cl->parent];
      regular(f, p, (cl->z[0] - pa->z[0]) * s, (cl->z[1] - pa->z[1]) * s, (cl->z[2] - pa->z[2]) * s,
              f->table);
      expand(f, p, f->table, 1, 0, f->full_a);
      expand(f, p, local_of(f, cl->parent), 0, 0, f->full_b);
      convolve(f, p, f->full_a, f->full_b, 0, fc);
      cl->local = 1;
    }
    if (cl->children != 0 || !cl->local) {
      continue;
    }
    expand(f, p, fc, 0, 0, f->full_c);
    for (size_t i = cl->begin; i < cl->end; ++i) {
      double g[3];
      gradient(f, f->full_c, (f->nx[i] - cl->z[0]) * s, (f->ny[i] - cl->z[1]) * s,
               (f->nz[i] - cl->z[2]) * s, g);
      f->ax[i] += f->far_scale * g[0];
      f->ay[i] += f->far_scale * g[1];
      f->az[i] += f->far_scale * g[2];
    }
  }
}

myb_status myb_internal_fmm_evaluate(myb_internal_fmm* fmm, size_t massive, const double* gm,
                                     const double* px, const double* py, const double* pz,
                                     double* a) {
  MYB_CHECK(fmm != NULL && gm != NULL && px != NULL && py != NULL && pz != NULL && a != NULL &&
            massive <= fmm->count);
  myb_internal_fmm* f = fmm;
  const size_t n = f->count;
  double largest = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double w = i < massive ? gm[i] : 0.0;
    if (!finite(px[i]) || !finite(py[i]) || !finite(pz[i]) || !finite(w) || w < 0.0) {
      return MYB_E_DOMAIN;
    }
    f->points[i] =
        (point){.x = px[i], .y = py[i], .z = pz[i], .w = w, .index = i, .test = i >= massive};
    largest = w > largest ? w : largest;
  }
  f->stats = (myb_internal_fmm_stats){.leaf_size = f->leaf, .direct_pairs = f->direct_pairs};
  const myb_status status = build(f);
  if (status != MYB_OK) {
    return status;
  }
  const uint64_t ones = ~UINT64_C(0);
  double all = 0.0;
  memcpy(&all, &ones, sizeof all);
  for (size_t i = 0; i < n; ++i) {
    const point* q = &f->points[i];
    f->nx[i] = q->x;
    f->ny[i] = q->y;
    f->nz[i] = q->z;
    f->nw[i] = q->w;
    f->nmask[i] = q->test ? 0.0 : all;
  }
  const double weight = largest > 0.0 ? power_of_two_at_least(largest) : 1.0;
  f->inv_weight = 1.0 / weight;
  f->far_scale = weight * f->inv_length * f->inv_length;
  upward(f);
  memset(f->local, 0, f->ncells * 2 * f->ncoef * sizeof(double));
  memset(f->ax, 0, 3 * n * sizeof(double));
  walk(f);
  downward(f);
  for (size_t i = 0; i < n; ++i) {
    const size_t k = f->points[i].index;
    a[3 * k] = f->ax[i];
    a[3 * k + 1] = f->ay[i];
    a[3 * k + 2] = f->az[i];
  }
  f->stats.cells = f->ncells;
  for (size_t c = 0; c < f->ncells; ++c) {
    const cell* cl = &f->cells[c];
    if (cl->children == 0) {
      f->stats.leaves++;
    }
    f->stats.depth = cl->level > f->stats.depth ? cl->level : f->stats.depth;
  }
  return MYB_OK;
}
