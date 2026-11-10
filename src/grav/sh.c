/* SPDX-License-Identifier: BSD-3-Clause */
/* Spherical-harmonic gravity fields: storage, ICGEM files, time-variable terms, evaluation. */
#include <math.h>
#include <stdckdint.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "grav/field_internal.h"
#include "myb/grav.h"
#include "myb/io.h"
#include "num/vec3.h"
#include "platform/platform.h"

enum { MAX_DEGREE = 2700, GRADIENT_SHIFT = 2 };

static const double scale_up = 0x1p930;    /* seeds of the recursion */
static const double scale_down = 0x1p-930; /* exact inverse */
static const double seconds_per_year = 365.25 * 86400.0;

typedef enum term_kind { TERM_TREND = 0, TERM_COS = 1, TERM_SIN = 2 } term_kind;

typedef struct tv_term {
  size_t index;
  int kind;
  double c;
  double s;
  double period; /* s (periodic terms) */
} tv_term;

struct myb_gravity_field {
  myb_allocator allocator;
  double gm;
  double radius;
  int degree;
  size_t count; /* coefficients: (degree + 1)(degree + 2) / 2 */
  double* c_static;
  double* s_static;
  double* c; /* effective */
  double* s;
  double* a;    /* recursion: Pbar_nm = a_nm t Pbar_n-1,m - b_nm Pbar_n-2,m */
  double* b;
  double* e;    /* dPbar_nm/dphi = e_nm Pbar_n,m+1 - m t Pbar_nm / u */
  double* sect; /* sqrt((2m + 1) / (2m)) */
  tv_term* terms;
  size_t term_count;
  double* t0; /* per coefficient: reference epoch (s past J2000 TT), NaN if none */
  int dirty_degree; /* highest degree touched by corrections since set_time, -1 none */
  char tide[16];
};

struct myb_gravity_work {
  myb_allocator allocator;
  int degree;
  double* p;     /* column m */
  double* pnext; /* column m + 1 */
  double* q;     /* column m divided by cos(phi) */
  double* seed;
  double* seed_u;
  double* cm;
  double* sm;
  double* pw;
};

static size_t index_of(int n, int m) {
  return (size_t)n * (size_t)(n + 1) / 2 + (size_t)m;
}

static int finite(double x) {
  return x - x == 0.0;
}

/* ---- Creation ------------------------------------------------------------------------------ */

static myb_status allocate(double gm, double radius, int degree, const myb_allocator* allocator,
                           myb_gravity_field** out) {
  if (!(gm > 0.0) || !finite(gm) || !(radius > 0.0) || !finite(radius) || degree < 0 ||
      degree > MAX_DEGREE) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_gravity_field), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_gravity_field* f = block;
  memset(f, 0, sizeof *f);
  f->allocator = a;
  f->gm = gm;
  f->radius = radius;
  f->degree = degree;
  f->dirty_degree = -1;
  f->count = index_of(degree + 1, 0);
  /* the recursion of the Legendre functions reaches two degrees above the field: the second
   * derivatives of a harmonic of degree n are harmonics of degree n + 2 */
  const size_t table = index_of(degree + GRADIENT_SHIFT + 1, 0);
  double** arrays[8] = {&f->c_static, &f->s_static, &f->c, &f->s, &f->e, &f->t0, &f->a, &f->b};
  for (int k = 0; k < 8 && status == MYB_OK; ++k) {
    const size_t n = k < 6 ? f->count : table;
    status = myb_internal_alloc_array(&a, n, sizeof(double), &block);
    if (status == MYB_OK) {
      *arrays[k] = block;
      memset(block, 0, n * sizeof(double));
    }
  }
  if (status == MYB_OK) {
    status = myb_internal_alloc_array(&a, (size_t)degree + GRADIENT_SHIFT + 2, sizeof(double),
                                      &block);
    if (status == MYB_OK) {
      f->sect = block;
    }
  }
  if (status != MYB_OK) {
    myb_gravity_field_destroy(f);
    return status;
  }
  for (size_t k = 0; k < f->count; ++k) {
    f->t0[k] = (double)NAN;
  }
  for (int n = 0; n <= degree + GRADIENT_SHIFT; ++n) {
    for (int m = 0; m <= n; ++m) {
      const double dn = n;
      const double dm = m;
      const size_t k = index_of(n, m);
      if (n >= m + 2) {
        f->a[k] = sqrt((2.0 * dn - 1.0) * (2.0 * dn + 1.0) / ((dn - dm) * (dn + dm)));
        f->b[k] = sqrt((2.0 * dn + 1.0) * (dn + dm - 1.0) * (dn - dm - 1.0) /
                       ((dn - dm) * (dn + dm) * (2.0 * dn - 3.0)));
      }
      if (n <= degree) {
        f->e[k] = sqrt((dn - dm) * (dn + dm + 1.0) * (m == 0 ? 0.5 : 1.0));
      }
    }
  }
  f->sect[0] = 1.0;
  for (int m = 1; m <= degree + GRADIENT_SHIFT + 1; ++m) {
    f->sect[m] = sqrt((2.0 * m + 1.0) / (2.0 * m));
  }
  *out = f;
  return MYB_OK;
}

myb_status myb_gravity_field_create(double gm, double radius, int degree, const double* c,
                                    const double* s, const myb_allocator* allocator,
                                    myb_gravity_field** field) {
  MYB_CHECK(c != NULL && s != NULL && field != NULL);
  myb_gravity_field* f = NULL;
  const myb_status status = allocate(gm, radius, degree, allocator, &f);
  if (status != MYB_OK) {
    return status;
  }
  for (size_t k = 0; k < f->count; ++k) {
    if (!finite(c[k]) || !finite(s[k])) {
      myb_gravity_field_destroy(f);
      return MYB_E_DOMAIN;
    }
  }
  memcpy(f->c_static, c, f->count * sizeof(double));
  memcpy(f->s_static, s, f->count * sizeof(double));
  memcpy(f->c, c, f->count * sizeof(double));
  memcpy(f->s, s, f->count * sizeof(double));
  *field = f;
  return MYB_OK;
}

void myb_gravity_field_destroy(myb_gravity_field* field) {
  if (field == NULL) {
    return;
  }
  const myb_allocator a = field->allocator;
  myb_internal_release(&a, field->c_static);
  myb_internal_release(&a, field->s_static);
  myb_internal_release(&a, field->c);
  myb_internal_release(&a, field->s);
  myb_internal_release(&a, field->a);
  myb_internal_release(&a, field->b);
  myb_internal_release(&a, field->e);
  myb_internal_release(&a, field->t0);
  myb_internal_release(&a, field->sect);
  myb_internal_release(&a, field->terms);
  myb_internal_release(&a, field);
}

int myb_gravity_field_degree(const myb_gravity_field* field) {
  MYB_CHECK(field != NULL);
  return field->degree;
}

double myb_gravity_field_gm(const myb_gravity_field* field) {
  MYB_CHECK(field != NULL);
  return field->gm;
}

double myb_gravity_field_radius(const myb_gravity_field* field) {
  MYB_CHECK(field != NULL);
  return field->radius;
}

const char* myb_gravity_field_tide_system(const myb_gravity_field* field) {
  MYB_CHECK(field != NULL);
  return field->tide;
}

/* ---- Time-variable terms, corrections ------------------------------------------------------ */

void myb_gravity_field_set_time(myb_gravity_field* field, myb_time t) {
  MYB_CHECK(field != NULL);
  myb_gravity_field* f = field;
  /* undo the corrections */
  if (f->dirty_degree >= 0) {
    const size_t n = index_of(f->dirty_degree + 1, 0);
    memcpy(f->c, f->c_static, n * sizeof(double));
    memcpy(f->s, f->s_static, n * sizeof(double));
    f->dirty_degree = -1;
  }
  if (f->term_count == 0) {
    return;
  }
  const double now = myb_time_seconds(t);
  /* every coefficient with terms restarts from its static value */
  for (size_t k = 0; k < f->term_count; ++k) {
    const size_t i = f->terms[k].index;
    f->c[i] = f->c_static[i];
    f->s[i] = f->s_static[i];
  }
  double last_argument = (double)NAN;
  double cos_arg = 1.0;
  double sin_arg = 0.0;
  for (size_t k = 0; k < f->term_count; ++k) {
    const tv_term* term = &f->terms[k];
    const double dt = now - f->t0[term->index];
    double factor;
    if (term->kind == TERM_TREND) {
      factor = dt / seconds_per_year;
      f->c[term->index] += term->c * factor;
      f->s[term->index] += term->s * factor;
      continue;
    }
    const double argument = 6.283185307179586 * (dt / term->period);
    if (argument != last_argument) {
      cos_arg = myb_cos(argument);
      sin_arg = myb_sin(argument);
      last_argument = argument;
    }
    factor = term->kind == TERM_COS ? cos_arg : sin_arg;
    f->c[term->index] += term->c * factor;
    f->s[term->index] += term->s * factor;
  }
}

myb_status myb_gravity_field_add_correction(myb_gravity_field* field, int n, int m, double dc,
                                            double ds) {
  MYB_CHECK(field != NULL);
  if (n < 0 || m < 0 || m > n || n > field->degree || !finite(dc) || !finite(ds)) {
    return MYB_E_DOMAIN;
  }
  const size_t k = index_of(n, m);
  field->c[k] += dc;
  field->s[k] += ds;
  if (n > field->dirty_degree) {
    field->dirty_degree = n;
  }
  return MYB_OK;
}

myb_status myb_gravity_field_coefficient(const myb_gravity_field* field, int n, int m, double* c,
                                         double* s) {
  MYB_CHECK(field != NULL && c != NULL && s != NULL);
  if (n < 0 || m < 0 || m > n || n > field->degree) {
    return MYB_E_DOMAIN;
  }
  *c = field->c[index_of(n, m)];
  *s = field->s[index_of(n, m)];
  return MYB_OK;
}

/* ---- Work memory ----------------------------------------------------------------------------- */

myb_status myb_gravity_work_create(int degree, const myb_allocator* allocator,
                                   myb_gravity_work** work) {
  MYB_CHECK(work != NULL);
  if (degree < 0 || degree > MAX_DEGREE) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_gravity_work), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_gravity_work* w = block;
  memset(w, 0, sizeof *w);
  w->allocator = a;
  w->degree = degree;
  const size_t n = (size_t)degree + GRADIENT_SHIFT + 2;
  status = myb_internal_alloc_array(&a, 8 * n, sizeof(double), &block);
  if (status != MYB_OK) {
    myb_internal_release(&a, w);
    return status;
  }
  double* d = block;
  w->p = d;
  w->pnext = d + n;
  w->q = d + 2 * n;
  w->seed = d + 3 * n;
  w->seed_u = d + 4 * n;
  w->cm = d + 5 * n;
  w->sm = d + 6 * n;
  w->pw = d + 7 * n;
  *work = w;
  return MYB_OK;
}

void myb_gravity_work_destroy(myb_gravity_work* work) {
  if (work == NULL) {
    return;
  }
  const myb_allocator a = work->allocator;
  myb_internal_release(&a, work->p);
  myb_internal_release(&a, work);
}

/* ---- Evaluation ------------------------------------------------------------------------------ */

myb_status myb_gravity_field_eval(const myb_gravity_field* field, myb_gravity_work* work,
                                  const double r[3], int n_min, int n_max, double* potential,
                                  double a[3]) {
  MYB_CHECK(field != NULL && work != NULL && r != NULL);
  const myb_gravity_field* f = field;
  const int top = n_max < 0 ? f->degree : n_max;
  const double rr = myb_vec3_norm(r);
  if (!myb_vec3_finite(r) || !(rr > 0.0) || n_min < 0 || top > f->degree ||
      top > work->degree || n_min > top + 1) {
    return MYB_E_DOMAIN;
  }
  const double rho = myb_hypot(r[0], r[1]);
  const double t = r[2] / rr; /* sin(latitude) */
  const double u = rho / rr;  /* cos(latitude) */
  const double cl = rho > 0.0 ? r[0] / rho : 1.0;
  const double sl = rho > 0.0 ? r[1] / rho : 0.0;
  const double ratio = f->radius / rr;
  double* pw = work->pw;
  double* cm = work->cm;
  double* sm = work->sm;
  double* seed = work->seed;
  double* seed_u = work->seed_u;
  pw[0] = 1.0;
  cm[0] = 1.0;
  sm[0] = 0.0;
  for (int n = 1; n <= top; ++n) {
    pw[n] = pw[n - 1] * ratio;
    cm[n] = cm[n - 1] * cl - sm[n - 1] * sl;
    sm[n] = sm[n - 1] * cl + cm[n - 1] * sl;
  }
  /* sectorial seeds Pbar_mm and Pbar_mm / u, scaled by 2^930 */
  seed[0] = scale_up;
  seed_u[0] = 0.0;
  if (top >= 1) {
    seed[1] = 1.7320508075688772 * u * scale_up; /* Pbar_11 = sqrt(3) u */
    seed_u[1] = 1.7320508075688772 * scale_up;
  }
  for (int m = 2; m <= top; ++m) {
    seed[m] = f->sect[m] * u * seed[m - 1];
    seed_u[m] = f->sect[m] * u * seed_u[m - 1];
  }
  double* p = work->p;
  double* pnext = work->pnext;
  double* q = work->q;
  for (int n = 0; n <= top + 1; ++n) {
    pnext[n] = 0.0;
  }
  double v_sum = 0.0;
  double r_sum = 0.0;
  double phi_sum = 0.0;
  double lambda_sum = 0.0;
  for (int m = top; m >= 0; --m) {
    /* column m of Pbar (scaled) and of Pbar / u */
    p[m] = seed[m];
    q[m] = seed_u[m];
    if (m + 1 <= top) {
      const double k = sqrt(2.0 * m + 3.0) * t;
      p[m + 1] = k * seed[m];
      q[m + 1] = k * seed_u[m];
    }
    for (int n = m + 2; n <= top; ++n) {
      const size_t i = index_of(n, m);
      p[n] = f->a[i] * t * p[n - 1] - f->b[i] * p[n - 2];
      q[n] = f->a[i] * t * q[n - 1] - f->b[i] * q[n - 2];
    }
    double a0 = 0.0;
    double b0 = 0.0;
    double a1 = 0.0;
    double b1 = 0.0;
    double a2 = 0.0;
    double b2 = 0.0;
    double a3 = 0.0;
    double b3 = 0.0;
    const int first = m > n_min ? m : n_min;
    for (int n = first; n <= top; ++n) {
      const size_t i = index_of(n, m);
      const double wc = pw[n] * f->c[i];
      const double ws = pw[n] * f->s[i];
      const double dp = f->e[i] * pnext[n] - m * t * q[n];
      a0 += wc * p[n];
      b0 += ws * p[n];
      a1 += (n + 1) * wc * p[n];
      b1 += (n + 1) * ws * p[n];
      a2 += wc * dp;
      b2 += ws * dp;
      a3 += wc * q[n];
      b3 += ws * q[n];
    }
    v_sum += (a0 * cm[m] + b0 * sm[m]) * scale_down;
    r_sum += (a1 * cm[m] + b1 * sm[m]) * scale_down;
    phi_sum += (a2 * cm[m] + b2 * sm[m]) * scale_down;
    lambda_sum += m * (b3 * cm[m] - a3 * sm[m]) * scale_down;
    /* the current column becomes column m + 1 of the next step */
    double* swap = pnext;
    pnext = p;
    p = swap;
    for (int n = 0; n < m; ++n) {
      pnext[n] = 0.0;
    }
  }
  const double k = f->gm / rr;
  if (potential != NULL) {
    *potential = k * v_sum;
  }
  if (a != NULL) {
    const double g = k / rr;
    const double gr = -g * r_sum;
    const double gphi = g * phi_sum;
    const double glambda = g * lambda_sum;
    a[0] = gr * u * cl - gphi * t * cl - glambda * sl;
    a[1] = gr * u * sl - gphi * t * sl + glambda * cl;
    a[2] = gr * t + gphi * u;
  }
  return MYB_OK;
}

/* ---- Gradient -------------------------------------------------------------------------------- */

/* The second derivatives of a solid harmonic of degree n are solid harmonics of degree n + 2:
 * with Z_nm = (R / r)^(n+1) P_nm(sin phi) e^(i m lambda) (P_nm unnormalised) and the operators
 * D+ = d/dx + i d/dy, D- = d/dx - i d/dy, Dz = d/dz (L. E. Cunningham, Celest. Mech. 2, 207,
 * 1970; O. Montenbruck, E. Gill, Satellite Orbits, 2000, eq. 3.33),
 *   R D+ Z_nm = -Z_n+1,m+1,   R D- Z_nm = (n-m+2)(n-m+1) Z_n+1,m-1,   R Dz Z_nm = -(n-m+1) Z_n+1,m
 * (D- Z_n0 = conj(D+ Z_n0)). Applied twice to V = (GM / R) Re sum K_nm Z_nm, K = C - i S, each
 * of the five independent components
 *   V_zz,  V_xz = Re (D+ + D-) Dz V / 2,  V_yz = Re (D+ - D-) Dz V / (2 i),
 *   V_xx - V_yy = Re (D+ D+ + D- D-) V / 2,  V_xy = Re (D+ D+ - D- D-) V / (4 i)
 * (V_xx + V_yy = -V_zz) is a series of the same form with the coefficients of degree n + 2 and
 * orders m, m +- 1, m +- 2 below, where the factors of the fully normalised functions are
 *   sqrt((2n+1)/(2n+5)) sqrt(product of four of the numbers n -+ m + 1 .. n -+ m + 4)
 * and sqrt(2) appears where an order 0 meets an order above it. Nothing is divided by
 * cos(phi): the series are evaluated by the column recursion of the potential, so the poles
 * are regular points and the accuracy is that of the potential. */

/* sqrt((2n + 1) / (2n + 5) a b c d) for the small integers a .. d (their product is exact). */
static double shift_factor(int n, int a, int b, int c, int d) {
  const double product = (double)a * (double)b * (double)c * (double)d;
  return sqrt(product * (2.0 * n + 1.0) / (2.0 * n + 5.0));
}

/* The five derived coefficients (C', S') of degree n + 2 and order m (K' = C' - i S'), in the
 * order zz, xz, yz, xx - yy, xy. */
static void shifted_coefficients(const myb_gravity_field* f, int n, int m, double cs[5][2]) {
  static const double root_half = 0.70710678118654752;
  static const double root_two = 1.4142135623730951;
  for (int k = 0; k < 5; ++k) {
    cs[k][0] = 0.0;
    cs[k][1] = 0.0;
  }
  if (m <= n) { /* the same order: Dz Dz */
    const size_t i = index_of(n, m);
    const double w = shift_factor(n, n - m + 1, n - m + 2, n + m + 1, n + m + 2);
    cs[0][0] += w * f->c[i];
    cs[0][1] += w * f->s[i];
  }
  if (m >= 1 && m - 1 <= n) { /* from the order below: D+ Dz */
    const int j = m - 1;
    const size_t i = index_of(n, j);
    const double w = shift_factor(n, n - j + 1, n + j + 1, n + j + 2, n + j + 3);
    if (j == 0) {
      cs[1][0] += root_half * w * f->c[i];
      cs[2][1] += root_half * w * f->c[i];
    } else {
      cs[1][0] += 0.5 * w * f->c[i];
      cs[1][1] += 0.5 * w * f->s[i];
      cs[2][0] -= 0.5 * w * f->s[i];
      cs[2][1] += 0.5 * w * f->c[i];
    }
  }
  if (m + 1 <= n) { /* from the order above: D- Dz */
    const int j = m + 1;
    const size_t i = index_of(n, j);
    const double w = (m == 0 ? root_two : 1.0) *
                     shift_factor(n, n + j + 1, n - j + 1, n - j + 2, n - j + 3);
    cs[1][0] -= 0.5 * w * f->c[i];
    cs[1][1] -= 0.5 * w * f->s[i];
    cs[2][0] -= 0.5 * w * f->s[i];
    cs[2][1] += 0.5 * w * f->c[i];
  }
  if (m >= 2 && m - 2 <= n) { /* from two orders below: D+ D+ */
    const int j = m - 2;
    const size_t i = index_of(n, j);
    const double w = shift_factor(n, n + j + 1, n + j + 2, n + j + 3, n + j + 4);
    if (j == 0) {
      cs[3][0] += root_half * w * f->c[i];
      cs[4][1] += 0.5 * root_half * w * f->c[i];
    } else {
      cs[3][0] += 0.5 * w * f->c[i];
      cs[3][1] += 0.5 * w * f->s[i];
      cs[4][0] -= 0.25 * w * f->s[i];
      cs[4][1] += 0.25 * w * f->c[i];
    }
  }
  if (m + 2 <= n) { /* from two orders above: D- D- */
    const int j = m + 2;
    const size_t i = index_of(n, j);
    const double w = (m == 0 ? root_two : 1.0) *
                     shift_factor(n, n - j + 1, n - j + 2, n - j + 3, n - j + 4);
    cs[3][0] += 0.5 * w * f->c[i];
    cs[3][1] += 0.5 * w * f->s[i];
    cs[4][0] += 0.25 * w * f->s[i];
    cs[4][1] -= 0.25 * w * f->c[i];
  }
  if (m == 1 && n >= 1) { /* D- D- of the order 1 passes through 0 to the conjugate of order 1 */
    const size_t i = index_of(n, 1);
    const double w = shift_factor(n, n, n + 1, n + 2, n + 3);
    cs[3][0] -= 0.5 * w * f->c[i];
    cs[3][1] += 0.5 * w * f->s[i];
    cs[4][0] -= 0.25 * w * f->s[i];
    cs[4][1] -= 0.25 * w * f->c[i];
  }
}

myb_status myb_gravity_field_gradient(const myb_gravity_field* field, myb_gravity_work* work,
                                      const double r[3], int n_min, int n_max,
                                      double gradient[9]) {
  MYB_CHECK(field != NULL && work != NULL && r != NULL && gradient != NULL);
  const myb_gravity_field* f = field;
  const int source_top = n_max < 0 ? f->degree : n_max;
  const double rr = myb_vec3_norm(r);
  if (!myb_vec3_finite(r) || !(rr > 0.0) || n_min < 0 || source_top > f->degree ||
      source_top > work->degree || n_min > source_top + 1 ||
      source_top + GRADIENT_SHIFT > MAX_DEGREE) {
    return MYB_E_DOMAIN;
  }
  const int top = source_top + GRADIENT_SHIFT;
  const double rho = myb_hypot(r[0], r[1]);
  const double t = r[2] / rr;
  const double u = rho / rr;
  const double cl = rho > 0.0 ? r[0] / rho : 1.0;
  const double sl = rho > 0.0 ? r[1] / rho : 0.0;
  const double ratio = f->radius / rr;
  double* pw = work->pw;
  double* cm = work->cm;
  double* sm = work->sm;
  double* seed = work->seed;
  pw[0] = 1.0;
  cm[0] = 1.0;
  sm[0] = 0.0;
  for (int n = 1; n <= top; ++n) {
    pw[n] = pw[n - 1] * ratio;
    cm[n] = cm[n - 1] * cl - sm[n - 1] * sl;
    sm[n] = sm[n - 1] * cl + cm[n - 1] * sl;
  }
  seed[0] = scale_up;
  seed[1] = 1.7320508075688772 * u * scale_up;
  for (int m = 2; m <= top; ++m) {
    seed[m] = f->sect[m] * u * seed[m - 1];
  }
  double* p = work->p;
  double sum[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  for (int m = top; m >= 0; --m) {
    p[m] = seed[m];
    if (m + 1 <= top) {
      p[m + 1] = sqrt(2.0 * m + 3.0) * t * seed[m];
    }
    for (int n = m + 2; n <= top; ++n) {
      const size_t i = index_of(n, m);
      p[n] = f->a[i] * t * p[n - 1] - f->b[i] * p[n - 2];
    }
    double part[5][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    const int lowest = n_min + GRADIENT_SHIFT;
    const int first = m > lowest ? m : lowest;
    for (int n = first; n <= top; ++n) {
      double cs[5][2];
      shifted_coefficients(f, n - GRADIENT_SHIFT, m, cs);
      const double w = pw[n] * p[n];
      for (int k = 0; k < 5; ++k) {
        part[k][0] += w * cs[k][0];
        part[k][1] += w * cs[k][1];
      }
    }
    for (int k = 0; k < 5; ++k) {
      sum[k] += (part[k][0] * cm[m] + part[k][1] * sm[m]) * scale_down;
    }
  }
  /* the series are potentials of GM / R^2: (GM / R^2) / r times the sums */
  const double scale = f->gm / (f->radius * f->radius * rr);
  const double zz = scale * sum[0];
  const double xz = scale * sum[1];
  const double yz = scale * sum[2];
  const double q = scale * sum[3];
  const double xy = scale * sum[4];
  gradient[0] = 0.5 * (q - zz);
  gradient[1] = xy;
  gradient[2] = xz;
  gradient[3] = xy;
  gradient[4] = -0.5 * (q + zz);
  gradient[5] = yz;
  gradient[6] = xz;
  gradient[7] = yz;
  gradient[8] = zz;
  return MYB_OK;
}

/* ---- ICGEM files ----------------------------------------------------------------------------- */

typedef struct line_cursor {
  const char* text;
  size_t size;
  size_t at;
  long number;
} line_cursor;

/* Next line [*begin, *end); 0 at the end of the text. */
static int next_line(line_cursor* c, const char** begin, size_t* length) {
  if (c->at >= c->size) {
    return 0;
  }
  const size_t start = c->at;
  size_t end = start;
  /* bound: at most the size of the text */
  while (end < c->size && c->text[end] != '\n') {
    ++end;
  }
  c->at = end + 1;
  c->number++;
  size_t n = end - start;
  if (n > 0 && c->text[start + n - 1] == '\r') {
    --n;
  }
  *begin = c->text + start;
  *length = n;
  return 1;
}

/* Splits a line into at most `max` whitespace-separated tokens. */
static int tokens(const char* s, size_t n, const char** tok, size_t* len, int max) {
  int count = 0;
  size_t i = 0;
  /* bound: every iteration consumes at least one character */
  while (i < n && count < max) {
    /* bound: skips blanks */
    while (i < n && (s[i] == ' ' || s[i] == '\t')) {
      ++i;
    }
    if (i == n) {
      break;
    }
    const size_t start = i;
    /* bound: a token */
    while (i < n && s[i] != ' ' && s[i] != '\t') {
      ++i;
    }
    tok[count] = s + start;
    len[count] = i - start;
    ++count;
  }
  return count;
}

static int number(const char* s, size_t n, double* value) {
  size_t used = 0;
  return myb_parse_double(s, n, value, &used) == MYB_OK && used == n && finite(*value);
}

static int is(const char* s, size_t n, const char* word) {
  return n == strlen(word) && memcmp(s, word, n) == 0;
}

/* yyyymmdd or yyyymmdd.hhmm -> seconds past J2000 (TT). */
static int epoch(const char* s, size_t n, double* seconds) {
  double v = 0.0;
  if (!number(s, n, &v) || v < 10000101.0 || v > 99991231.2400) {
    return 0;
  }
  const double day_part = floor(v);
  const double fraction = v - day_part;
  const long ymd = (long)day_part;
  const int hhmm = (int)floor(fraction * 10000.0 + 0.5);
  const myb_calendar c = {ymd / 10000, (int)(ymd / 100 % 100), (int)(ymd % 100), hhmm / 100,
                          hhmm % 100, 0.0};
  myb_time t;
  if (myb_calendar_to_time(&c, &t) != MYB_OK) {
    return 0;
  }
  *seconds = myb_time_seconds(t);
  return 1;
}

static myb_status fail(myb_diag* diag, long line, const char* what) {
  myb_internal_diag(diag, line, "ICGEM: %s", what);
  return MYB_E_FORMAT;
}

myb_status myb_gravity_field_load_icgem(const char* text, size_t size, int max_degree,
                                        const myb_allocator* allocator, myb_gravity_field** field,
                                        myb_diag* diag) {
  MYB_CHECK(text != NULL && field != NULL);
  myb_internal_diag_clear(diag);
  line_cursor cur = {text, size, 0, 0};
  const char* line = NULL;
  size_t length = 0;
  double gm = 0.0;
  double radius = 0.0;
  double degree_value = -1.0;
  int normalized = 0;
  int in_head = 0;
  char tide[16] = "";
  const char* tok[12];
  size_t len[12];
  /* header: keys until end_of_head */
  int header_done = 0;
  /* bound: one line per iteration */
  while (!header_done && next_line(&cur, &line, &length)) {
    const int k = tokens(line, length, tok, len, 12);
    if (k == 0) {
      continue;
    }
    if (len[0] >= 11 && memcmp(tok[0], "end_of_head", 11) == 0) {
      header_done = 1;
      break;
    }
    if (len[0] >= 13 && memcmp(tok[0], "begin_of_head", 13) == 0) {
      in_head = 1;
      continue;
    }
    if (k < 2) {
      continue;
    }
    if (is(tok[0], len[0], "earth_gravity_constant") || is(tok[0], len[0], "gravity_constant")) {
      if (!number(tok[1], len[1], &gm)) {
        return fail(diag, cur.number, "invalid gravity constant");
      }
    } else if (is(tok[0], len[0], "radius")) {
      if (!number(tok[1], len[1], &radius)) {
        return fail(diag, cur.number, "invalid radius");
      }
    } else if (is(tok[0], len[0], "max_degree")) {
      if (!number(tok[1], len[1], &degree_value)) {
        return fail(diag, cur.number, "invalid max_degree");
      }
    } else if (is(tok[0], len[0], "norm")) {
      normalized = is(tok[1], len[1], "fully_normalized");
      if (!normalized) {
        return fail(diag, cur.number, "only fully_normalized coefficients are supported");
      }
    } else if (is(tok[0], len[0], "tide_system")) {
      const size_t n = len[1] < sizeof tide - 1 ? len[1] : sizeof tide - 1;
      memcpy(tide, tok[1], n);
      tide[n] = '\0';
    }
  }
  (void)in_head;
  if (!header_done || !(gm > 0.0) || !(radius > 0.0) || degree_value < 0.0 ||
      degree_value > (double)MAX_DEGREE || degree_value != floor(degree_value)) {
    return fail(diag, cur.number, "header incomplete (gravity constant, radius, max_degree)");
  }
  if (!normalized) {
    return fail(diag, cur.number, "norm missing (fully_normalized required)");
  }
  int degree = (int)degree_value;
  if (max_degree > 0 && max_degree < degree) {
    degree = max_degree;
  }
  myb_gravity_field* f = NULL;
  myb_status status = allocate(gm, radius, degree, allocator, &f);
  if (status != MYB_OK) {
    return status;
  }
  memcpy(f->tide, tide, sizeof tide);
  /* data: first pass counts the time-variable terms */
  const line_cursor data_start = cur;
  size_t term_count = 0;
  /* bound: one line per iteration */
  while (next_line(&cur, &line, &length)) {
    const int k = tokens(line, length, tok, len, 1);
    if (k == 1 && (is(tok[0], len[0], "trnd") || is(tok[0], len[0], "acos") ||
                   is(tok[0], len[0], "asin"))) {
      ++term_count;
    }
  }
  if (term_count > 0) {
    void* block = NULL;
    status = myb_internal_alloc_array(&f->allocator, term_count, sizeof(tv_term), &block);
    if (status != MYB_OK) {
      myb_gravity_field_destroy(f);
      return status;
    }
    f->terms = block;
  }
  cur = data_start;
  /* bound: one line per iteration */
  while (next_line(&cur, &line, &length)) {
    const int k = tokens(line, length, tok, len, 12);
    if (k == 0) {
      continue;
    }
    const int gfc = is(tok[0], len[0], "gfc");
    const int gfct = is(tok[0], len[0], "gfct");
    const int trnd = is(tok[0], len[0], "trnd");
    const int acos_line = is(tok[0], len[0], "acos");
    const int asin_line = is(tok[0], len[0], "asin");
    if (!gfc && !gfct && !trnd && !acos_line && !asin_line) {
      continue;
    }
    double nd = 0.0;
    double md = 0.0;
    double cv = 0.0;
    double sv = 0.0;
    if (k < 5 || !number(tok[1], len[1], &nd) || !number(tok[2], len[2], &md) ||
        !number(tok[3], len[3], &cv) || !number(tok[4], len[4], &sv) || nd != floor(nd) ||
        md != floor(md) || md < 0.0 || md > nd || nd > (double)MAX_DEGREE) {
      myb_gravity_field_destroy(f);
      return fail(diag, cur.number, "malformed coefficient line");
    }
    const int n = (int)nd;
    const int m = (int)md;
    if (n > degree) {
      continue;
    }
    const size_t i = index_of(n, m);
    if (gfc || gfct) {
      f->c_static[i] = cv;
      f->s_static[i] = sv;
      if (gfct) {
        double t0 = 0.0;
        if (k < 8 || !epoch(tok[7], len[7], &t0)) {
          myb_gravity_field_destroy(f);
          return fail(diag, cur.number, "gfct without a valid epoch (yyyymmdd)");
        }
        if (finite(f->t0[i])) {
          myb_gravity_field_destroy(f);
          return fail(diag, cur.number, "piecewise models (several gfct) are not supported");
        }
        f->t0[i] = t0;
      }
      continue;
    }
    tv_term* term = &f->terms[f->term_count];
    term->index = i;
    term->c = cv;
    term->s = sv;
    term->kind = trnd ? TERM_TREND : (acos_line ? TERM_COS : TERM_SIN);
    term->period = 0.0;
    if (!trnd) {
      double years = 0.0;
      if (k < 8 || !number(tok[7], len[7], &years) || !(years > 0.0)) {
        myb_gravity_field_destroy(f);
        return fail(diag, cur.number, "periodic term without a valid period (years)");
      }
      term->period = years * seconds_per_year;
    }
    f->term_count++;
  }
  for (size_t k = 0; k < f->term_count; ++k) {
    if (!finite(f->t0[f->terms[k].index])) {
      myb_gravity_field_destroy(f);
      return fail(diag, 0, "time-variable term without the epoch of its gfct line");
    }
  }
  memcpy(f->c, f->c_static, f->count * sizeof(double));
  memcpy(f->s, f->s_static, f->count * sizeof(double));
  *field = f;
  return MYB_OK;
}

myb_status myb_gravity_field_load_icgem_file(const char* path, int max_degree,
                                             const myb_allocator* allocator,
                                             myb_gravity_field** field, myb_diag* diag) {
  MYB_CHECK(path != NULL && field != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &a, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded = myb_gravity_field_load_icgem((const char*)map.data, map.size,
                                                         max_degree, allocator, field, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

/* ---- Internal access ------------------------------------------------------------------------- */

double* myb_internal_field_c(myb_gravity_field* field) {
  MYB_CHECK(field != NULL);
  return field->c;
}

double* myb_internal_field_s(myb_gravity_field* field) {
  MYB_CHECK(field != NULL);
  return field->s;
}

void myb_internal_field_touched(myb_gravity_field* field, int degree) {
  MYB_CHECK(field != NULL && degree >= 0 && degree <= field->degree);
  if (degree > field->dirty_degree) {
    field->dirty_degree = degree;
  }
}

myb_status myb_gravity_field_clone(const myb_gravity_field* field, const myb_allocator* allocator,
                                   myb_gravity_field** copy) {
  MYB_CHECK(field != NULL && copy != NULL);
  myb_gravity_field* f = NULL;
  myb_status status = allocate(field->gm, field->radius, field->degree, allocator, &f);
  if (status != MYB_OK) {
    return status;
  }
  const size_t bytes = field->count * sizeof(double);
  memcpy(f->c_static, field->c_static, bytes);
  memcpy(f->s_static, field->s_static, bytes);
  memcpy(f->c, field->c, bytes);
  memcpy(f->s, field->s, bytes);
  memcpy(f->t0, field->t0, bytes);
  memcpy(f->tide, field->tide, sizeof f->tide);
  f->dirty_degree = field->dirty_degree;
  if (field->term_count > 0) {
    void* block = NULL;
    status = myb_internal_alloc_array(&f->allocator, field->term_count, sizeof(tv_term), &block);
    if (status != MYB_OK) {
      myb_gravity_field_destroy(f);
      return status;
    }
    f->terms = block;
    memcpy(f->terms, field->terms, field->term_count * sizeof(tv_term));
    f->term_count = field->term_count;
  }
  *copy = f;
  return MYB_OK;
}
