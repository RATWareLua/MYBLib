/* SPDX-License-Identifier: BSD-3-Clause */
/* A small dense solver of nonlinear programs (myb/optimal.h): the method of multipliers
 * outside, a projected quasi-Newton iteration in the box inside. Every loop runs a fixed
 * sequence of operations: the result is the same on every platform. */
#include <math.h> /* sqrt, fabs, fmin, fmax */
#include <string.h>

#include "core/alloc.h"
#include "core/check.h"
#include "myb/optimal.h"

enum {
  HALVINGS = 40,     /* of a step that does not go down */
  OUTER = 50,
  INNER = 1000
};

static const double armijo = 1e-4;
static const double weight_factor = 10.0;
static const double weight_max = 1e12;

typedef struct solver {
  const myb_nlp_problem* problem;
  size_t n;
  size_t m;
  double weight;     /* rho */
  double* multiplier; /* m */
  double* c;         /* m */
  double* jacobian;  /* m n */
  double* gradient_f; /* n */
  double* b;         /* n n: the quasi-Newton matrix */
  double* factor;    /* n n: its free block, factored */
  double* g;         /* n: the gradient of the augmented Lagrangian */
  double* g_new;
  double* d;
  double* x_new;
  double* s;
  double* y;
  double* bs;
  double* rhs;
  size_t* free_index;
  int fresh;         /* the matrix is a multiple of the unit one */
  int evaluations;
  double objective;  /* at the last point evaluated */
} solver;

static double lower_of(const solver* s, size_t i) {
  return s->problem->lower != NULL ? s->problem->lower[i] : -HUGE_VAL;
}

static double upper_of(const solver* s, size_t i) {
  return s->problem->upper != NULL ? s->problem->upper[i] : HUGE_VAL;
}

static double into_box(const solver* s, size_t i, double value) {
  const double lo = lower_of(s, i);
  const double hi = upper_of(s, i);
  return value < lo ? lo : (value > hi ? hi : value);
}

/* The augmented Lagrangian at x, with its gradient when `gradient` is not NULL. */
static myb_status lagrangian(solver* s, const double* x, double* value, double* gradient) {
  const myb_nlp_problem* p = s->problem;
  double f = 0.0;
  const myb_status st =
      p->function(p->context, x, &f, s->c, gradient != NULL ? s->gradient_f : NULL,
                  gradient != NULL ? s->jacobian : NULL);
  ++s->evaluations;
  if (st != MYB_OK) {
    return st;
  }
  if (!isfinite(f)) {
    return MYB_E_DOMAIN;
  }
  s->objective = f;
  double sum = f;
  if (gradient != NULL) {
    memcpy(gradient, s->gradient_f, s->n * sizeof *gradient);
  }
  for (size_t i = 0; i < s->m; ++i) {
    const double c = s->c[i];
    if (!isfinite(c)) {
      return MYB_E_DOMAIN;
    }
    const double shifted = s->multiplier[i] + s->weight * c;
    double pull = shifted; /* the multiplier this point would have */
    if (i < p->equalities) {
      sum += c * (s->multiplier[i] + 0.5 * s->weight * c);
    } else {
      pull = shifted > 0.0 ? shifted : 0.0;
      sum += (pull * pull - s->multiplier[i] * s->multiplier[i]) / (2.0 * s->weight);
    }
    if (gradient != NULL && pull != 0.0) {
      const double* row = s->jacobian + i * s->n;
      for (size_t j = 0; j < s->n; ++j) {
        gradient[j] += pull * row[j];
      }
    }
  }
  *value = sum;
  return MYB_OK;
}

/* How far x is from a point where no feasible step goes down: the largest component of
 * x - P(x - g). */
static double critical(const solver* s, const double* x, const double* g) {
  double worst = 0.0;
  for (size_t i = 0; i < s->n; ++i) {
    worst = fmax(worst, fabs(x[i] - into_box(s, i, x[i] - g[i])));
  }
  return worst;
}

static void reset_matrix(solver* s, double scale) {
  memset(s->b, 0, s->n * s->n * sizeof *s->b);
  for (size_t i = 0; i < s->n; ++i) {
    s->b[i * s->n + i] = scale;
  }
  s->fresh = 1;
}

/* A = L L^T in place (lower part); 0 if A is not positive definite. */
static int factorise(size_t n, double* a) {
  for (size_t j = 0; j < n; ++j) {
    double d = a[j * n + j];
    for (size_t k = 0; k < j; ++k) {
      d -= a[j * n + k] * a[j * n + k];
    }
    if (!(d > 0.0) || !isfinite(d)) {
      return 0;
    }
    const double l = sqrt(d);
    a[j * n + j] = l;
    for (size_t i = j + 1; i < n; ++i) {
      double sum = a[i * n + j];
      for (size_t k = 0; k < j; ++k) {
        sum -= a[i * n + k] * a[j * n + k];
      }
      a[i * n + j] = sum / l;
    }
  }
  return 1;
}

static void factor_solve(size_t n, const double* l, double* b) {
  for (size_t i = 0; i < n; ++i) {
    double sum = b[i];
    for (size_t k = 0; k < i; ++k) {
      sum -= l[i * n + k] * b[k];
    }
    b[i] = sum / l[i * n + i];
  }
  for (size_t i = n; i-- > 0;) {
    double sum = b[i];
    for (size_t k = i + 1; k < n; ++k) {
      sum -= l[k * n + i] * b[k];
    }
    b[i] = sum / l[i * n + i];
  }
}

/* The step d from x: Newton's with the matrix for the variables that are free, the gradient
 * for those held at a bound that the gradient pushes them against (D. P. Bertsekas, 1982).
 * 0 if the matrix gives no way down on the free ones. */
static int direction(solver* s, const double* x, const double* g, double hold) {
  size_t free_count = 0;
  for (size_t i = 0; i < s->n; ++i) {
    const int held = (x[i] <= lower_of(s, i) + hold && g[i] > 0.0) ||
                     (x[i] >= upper_of(s, i) - hold && g[i] < 0.0);
    if (!held) {
      s->free_index[free_count++] = i;
    }
    s->d[i] = -g[i];
  }
  if (free_count == 0) {
    return 1;
  }
  for (size_t a = 0; a < free_count; ++a) {
    for (size_t b = 0; b < free_count; ++b) {
      s->factor[a * free_count + b] = s->b[s->free_index[a] * s->n + s->free_index[b]];
    }
    s->rhs[a] = -g[s->free_index[a]];
  }
  if (!factorise(free_count, s->factor)) {
    return 0;
  }
  factor_solve(free_count, s->factor, s->rhs);
  double slope = 0.0;
  for (size_t a = 0; a < free_count; ++a) {
    slope += g[s->free_index[a]] * s->rhs[a];
  }
  if (!(slope <= 0.0)) {
    return 0; /* also for a gradient that is not finite */
  }
  for (size_t a = 0; a < free_count; ++a) {
    s->d[s->free_index[a]] = s->rhs[a];
  }
  return 1;
}

/* The BFGS update of the matrix with the step s->s and the change of the gradient s->y,
 * damped so that it stays positive definite (M. J. D. Powell, 1978). */
static void update_matrix(solver* s) {
  const size_t n = s->n;
  double sy = 0.0;
  double ss = 0.0;
  double yy = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sy += s->s[i] * s->y[i];
    ss += s->s[i] * s->s[i];
    yy += s->y[i] * s->y[i];
  }
  if (!(ss > 0.0)) {
    return;
  }
  if (s->fresh && sy > 0.0 && yy > 0.0) {
    reset_matrix(s, yy / sy); /* the size of the curvature along the first step */
  }
  s->fresh = 0;
  double sbs = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double sum = 0.0;
    for (size_t j = 0; j < n; ++j) {
      sum += s->b[i * n + j] * s->s[j];
    }
    s->bs[i] = sum;
    sbs += s->s[i] * sum;
  }
  if (!(sbs > 0.0) || !isfinite(sbs)) {
    reset_matrix(s, 1.0);
    return;
  }
  const double theta = sy >= 0.2 * sbs ? 1.0 : 0.8 * sbs / (sbs - sy);
  double sr = 0.0;
  for (size_t i = 0; i < n; ++i) {
    s->y[i] = theta * s->y[i] + (1.0 - theta) * s->bs[i];
    sr += s->s[i] * s->y[i];
  }
  if (!(sr > 0.0) || !isfinite(sr)) {
    reset_matrix(s, 1.0);
    return;
  }
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      s->b[i * n + j] += s->y[i] * s->y[j] / sr - s->bs[i] * s->bs[j] / sbs;
    }
  }
}

/* Minimises the augmented Lagrangian in the box from x until the criticality is within the
 * tolerance, the iterations run out, or no step goes down any more. */
static myb_status minimise(solver* s, double* x, double tolerance, int max_iterations,
                           double* criticality, int* iterations) {
  const size_t n = s->n;
  double value = 0.0;
  myb_status st = lagrangian(s, x, &value, s->g);
  if (st != MYB_OK) {
    return st;
  }
  *criticality = critical(s, x, s->g);
  /* bound: max_iterations */
  for (int it = 0; it < max_iterations && *criticality > tolerance; ++it) {
    ++*iterations;
    if (!direction(s, x, s->g, fmin(1e-8, *criticality))) {
      reset_matrix(s, 1.0);
      if (!direction(s, x, s->g, fmin(1e-8, *criticality))) {
        return MYB_E_DOMAIN; /* a gradient that is not finite */
      }
    }
    double scale = 1.0;
    double value_new = 0.0;
    int accepted = 0;
    int with_gradient = 0;
    /* bound: HALVINGS */
    for (int half = 0; half < HALVINGS && !accepted; ++half) {
      double slope = 0.0;
      for (size_t i = 0; i < n; ++i) {
        s->x_new[i] = into_box(s, i, x[i] + scale * s->d[i]);
        slope += s->g[i] * (s->x_new[i] - x[i]);
      }
      if (slope < 0.0) {
        /* the full step usually stands: its gradient is asked for at once */
        st = lagrangian(s, s->x_new, &value_new, half == 0 ? s->g_new : NULL);
        if (st == MYB_OK && value_new <= value + armijo * slope) {
          accepted = 1;
          with_gradient = half == 0;
        }
      }
      scale *= 0.5;
    }
    if (!accepted) {
      if (s->fresh) {
        break; /* not even along the gradient: the function resolves no further */
      }
      reset_matrix(s, 1.0);
      continue;
    }
    if (!with_gradient) {
      st = lagrangian(s, s->x_new, &value_new, s->g_new);
      if (st != MYB_OK) {
        return st;
      }
    }
    for (size_t i = 0; i < n; ++i) {
      s->s[i] = s->x_new[i] - x[i];
      s->y[i] = s->g_new[i] - s->g[i];
    }
    update_matrix(s);
    memcpy(x, s->x_new, n * sizeof *x);
    memcpy(s->g, s->g_new, n * sizeof *s->g);
    value = value_new;
    *criticality = critical(s, x, s->g);
  }
  /* leave s->c and the objective at x */
  return lagrangian(s, x, &value, s->g);
}

static int problem_valid(const myb_nlp_problem* p) {
  if (p->function == NULL || p->n < 1 || p->n > MYB_NLP_MAX_UNKNOWNS ||
      p->equalities > MYB_NLP_MAX_CONSTRAINTS || p->inequalities > MYB_NLP_MAX_CONSTRAINTS ||
      p->equalities + p->inequalities > MYB_NLP_MAX_CONSTRAINTS) {
    return 0;
  }
  for (size_t i = 0; i < p->n; ++i) {
    const double lo = p->lower != NULL ? p->lower[i] : -HUGE_VAL;
    const double hi = p->upper != NULL ? p->upper[i] : HUGE_VAL;
    if (!(lo <= hi)) {
      return 0;
    }
  }
  return 1;
}

myb_status myb_nlp_solve(const myb_nlp_problem* problem, const myb_nlp_options* options,
                         double* x, double* multipliers, const myb_allocator* allocator,
                         myb_nlp_result* result) {
  MYB_CHECK(problem != NULL && x != NULL && result != NULL);
  memset(result, 0, sizeof *result);
  const myb_nlp_options defaults = {0};
  const myb_nlp_options* o = options != NULL ? options : &defaults;
  if (!problem_valid(problem) || !(o->feasibility >= 0.0) || !(o->optimality >= 0.0) ||
      o->max_outer < 0 || o->max_inner < 0) {
    return MYB_E_DOMAIN;
  }
  const size_t n = problem->n;
  const size_t m = problem->equalities + problem->inequalities;
  const double feasibility = o->feasibility > 0.0 ? o->feasibility : 1e-8;
  const double optimality = o->optimality > 0.0 ? o->optimality : 1e-6;
  const int max_outer = o->max_outer > 0 ? o->max_outer : OUTER;
  const int max_inner = o->max_inner > 0 ? o->max_inner : INNER;
  const myb_allocator a = myb_internal_allocator(allocator);
  const size_t doubles = 2 * n * n + 9 * n + 2 * m + m * n;
  void* block = NULL;
  void* indices = NULL;
  myb_status st = myb_internal_alloc_array(&a, doubles, sizeof(double), &block);
  if (st == MYB_OK) {
    st = myb_internal_alloc_array(&a, n, sizeof(size_t), &indices);
  }
  if (st != MYB_OK) {
    myb_internal_release(&a, block);
    return st;
  }
  double* w = block;
  solver s = {.problem = problem, .n = n, .m = m};
  s.b = w;
  s.factor = s.b + n * n;
  s.jacobian = s.factor + n * n;
  s.multiplier = s.jacobian + m * n;
  s.c = s.multiplier + m;
  s.gradient_f = s.c + m;
  s.g = s.gradient_f + n;
  s.g_new = s.g + n;
  s.d = s.g_new + n;
  s.x_new = s.d + n;
  s.s = s.x_new + n;
  s.y = s.s + n;
  s.bs = s.y + n;
  s.rhs = s.bs + n;
  s.free_index = indices;
  memset(s.multiplier, 0, m * sizeof *s.multiplier);
  for (size_t i = 0; i < n; ++i) {
    x[i] = into_box(&s, i, x[i]);
  }
  /* the first weight from the sizes of the objective and of the violation at the start
   * (E. G. Birgin, J. M. Martinez, Practical Augmented Lagrangian Methods, 2014, §12.4) */
  s.weight = 1.0;
  double value = 0.0;
  st = lagrangian(&s, x, &value, NULL);
  if (st == MYB_OK) {
    double violation = 0.0;
    for (size_t i = 0; i < m; ++i) {
      const double c = i < problem->equalities ? s.c[i] : fmax(s.c[i], 0.0);
      violation += 0.5 * c * c;
    }
    s.weight = fmin(1e8, fmax(1e-8, 10.0 * fmax(1.0, fabs(s.objective)) / fmax(1.0, violation)));
    reset_matrix(&s, 1.0);
  }
  double tolerance = 1e-3;
  double measure_before = HUGE_VAL;
  double measure = HUGE_VAL;
  double criticality = HUGE_VAL;
  int converged = 0;
  /* bound: max_outer */
  for (int outer = 0; outer < max_outer && st == MYB_OK && !converged; ++outer) {
    ++result->outer;
    tolerance = fmax(optimality, tolerance);
    st = minimise(&s, x, tolerance, max_inner, &criticality, &result->inner);
    if (st != MYB_OK) {
      break;
    }
    /* how far the constraints are from holding, the inequalities with their multipliers:
     * either the constraint is at zero or its multiplier is */
    measure = 0.0;
    for (size_t i = 0; i < m; ++i) {
      const double c = s.c[i];
      if (i < problem->equalities) {
        measure = fmax(measure, fabs(c));
        s.multiplier[i] += s.weight * c;
      } else {
        measure = fmax(measure, fabs(fmax(c, -s.multiplier[i] / s.weight)));
        s.multiplier[i] = fmax(0.0, s.multiplier[i] + s.weight * c);
      }
    }
    converged = measure <= feasibility && criticality <= optimality;
    if (!converged) {
      if (measure > 0.5 * measure_before && s.weight < weight_max) {
        s.weight *= weight_factor;
        reset_matrix(&s, 1.0);
      }
      measure_before = measure;
      tolerance *= 0.1;
    }
  }
  result->objective = s.objective;
  result->violation = measure;
  result->stationarity = criticality;
  result->evaluations = s.evaluations;
  if (multipliers != NULL) {
    memcpy(multipliers, s.multiplier, m * sizeof *multipliers);
  }
  myb_internal_release(&a, indices);
  myb_internal_release(&a, block);
  if (st != MYB_OK) {
    return st;
  }
  return converged ? MYB_OK : MYB_E_NOCONV;
}
