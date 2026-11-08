/* SPDX-License-Identifier: BSD-3-Clause */
/* Gauss-Jackson of order 8 in summed ordinate form (M. M. Berry, L. M. Healy, J. Astronaut.
 * Sci. 52, 331, 2004), with the coefficients of integ/gj_coefficients.h.
 *
 * Grid t_k = t0 + k h. With the first and second sums (compensated)
 *   s(k + 1) = s(k) + h a(k),   S(k + 1) = S(k) + h s(k + 1),
 * the state at window position j around the centre c is
 *   v(c + j) = s(c + j) + h   sum_k V[j + 4][k] a(c + k),
 *   x(c + j) = S(c + j) + h^2 sum_k X[j + 4][k] a(c + k),     k = -4..4.
 * A step predicts node c + 5 (j = 5), evaluates a there, moves the centre to c + 1 and corrects
 * the node as j = 4 until a stops changing. The start finds the nodes -4..4 around t0 by
 * iterating the same formulas from a Taylor guess. */
#include <math.h> /* fabs, ceil */
#include <stdckdint.h>
#include <stdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "integ/event.h"
#include "integ/gj_coefficients.h"
#include "myb/integ.h"
#include "num/eft.h"

enum {
  WINDOW = 9,
  HALF = 4,
  START_ITERATIONS = 60,
  DEFAULT_CORRECTIONS = 10,
};

static const double default_tolerance = 1e-15;
/* A stalled iteration counts as converged only at round-off level: a growing change means
 * the fixed point iteration diverges (the step is too long). */
static const double stall_limit = 1e-10;

struct myb_gauss_jackson {
  myb_allocator allocator;
  size_t n;
  double h;
  double tolerance;
  int max_corrections;
  myb_accel_function f;
  void* context;
  double t0;
  int64_t centre;  /* index of the window centre */
  int64_t current; /* index of the node reported as the state */
  int started;
  myb_gauss_jackson_stats stats;
  double* block; /* the single allocation of all vectors (the window pointers rotate) */
  double* window[WINDOW + 1]; /* a at nodes centre - 4 .. centre + 4, and the node in work */
  double* s;      /* first sum at the centre, n values ... */
  double* s_comp; /* ... and its compensation */
  double* big_s;  /* second sum at the centre */
  double* big_s_comp;
  double* y;      /* 2n: state at the current node */
  double* work;   /* 2n: state at the node in work */
  double* scratch; /* 2n: events */
  double* previous; /* n: a of the previous corrector iteration */
};

static int finite(double x) {
  return x - x == 0.0;
}

static double larger_of(double a, double b) {
  return a > b ? a : b;
}

myb_status myb_gauss_jackson_create(size_t n, const myb_gauss_jackson_options* options,
                                    const myb_allocator* allocator,
                                    myb_gauss_jackson** integrator) {
  MYB_CHECK(options != NULL && integrator != NULL);
  const myb_gauss_jackson_options o = *options;
  if (n == 0 || !finite(o.h) || o.h == 0.0 || !(o.tolerance >= 0.0) || o.max_corrections < 0) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_gauss_jackson), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_gauss_jackson* s = block;
  *s = (myb_gauss_jackson){.allocator = a,
                           .n = n,
                           .h = o.h,
                           .tolerance = o.tolerance > 0.0 ? o.tolerance : default_tolerance,
                           .max_corrections =
                               o.max_corrections > 0 ? o.max_corrections : DEFAULT_CORRECTIONS};
  /* 10 acceleration vectors, 4 sums, 3 states of 2n, previous: 21 n values. */
  size_t count = 0;
  if (ckd_mul(&count, (size_t)(WINDOW + 1 + 4 + 3 * 2 + 1), n)) {
    status = MYB_E_OVERFLOW;
    goto fail;
  }
  status = myb_internal_alloc_array(&a, count, sizeof(double), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  double* p = block;
  memset(p, 0, count * sizeof(double));
  s->block = p;
  for (int k = 0; k <= WINDOW; ++k) {
    s->window[k] = p;
    p += n;
  }
  double** single[5] = {&s->s, &s->s_comp, &s->big_s, &s->big_s_comp, &s->previous};
  for (int k = 0; k < 5; ++k) {
    *single[k] = p;
    p += n;
  }
  double** wide[3] = {&s->y, &s->work, &s->scratch};
  for (int k = 0; k < 3; ++k) {
    *wide[k] = p;
    p += 2 * n;
  }
  *integrator = s;
  return MYB_OK;

fail:
  myb_internal_release(&a, s);
  return status;
}

void myb_gauss_jackson_destroy(myb_gauss_jackson* integrator) {
  if (integrator == NULL) {
    return;
  }
  const myb_allocator a = integrator->allocator;
  myb_internal_release(&a, integrator->block);
  myb_internal_release(&a, integrator);
}

static double node_time(const myb_gauss_jackson* s, int64_t index) {
  return s->t0 + (double)index * s->h;
}

static myb_status evaluate(myb_gauss_jackson* s, int64_t index, const double* y, double* a) {
  s->stats.evaluations++;
  return s->f(s->context, node_time(s, index), y, y + s->n, a);
}

/* First and second sums of component i at window position j (-4..5) from those at the centre,
 * accumulated with the compensations; the result is rounded once. */
static void sums_at(const myb_gauss_jackson* s, size_t i, int j, double* first, double* second) {
  const double h = s->h;
  double v = s->s[i];
  double vc = s->s_comp[i];
  double x = s->big_s[i];
  double xc = s->big_s_comp[i];
  if (j >= 0) {
    for (int m = 0; m < j; ++m) { /* s(c+m+1) = s(c+m) + h a(c+m); S += h s(c+m+1) */
      v = myb_eft_two_sum(v, h * s->window[HALF + m][i] + vc, &vc);
      x = myb_eft_two_sum(x, h * (v + vc) + xc, &xc);
    }
  } else {
    for (int m = 0; m < -j; ++m) { /* S(c-m-1) = S(c-m) - h s(c-m); s(c-m-1) = s(c-m) - h a */
      x = myb_eft_two_sum(x, -(h * (v + vc)) + xc, &xc);
      v = myb_eft_two_sum(v, -(h * s->window[HALF - 1 - m][i]) + vc, &vc);
    }
  }
  *first = v + vc;
  *second = x + xc;
}

/* State at window position j (-4..5) into y (2n). */
static void state_at_node(const myb_gauss_jackson* s, int j, double* y) {
  const size_t n = s->n;
  const double h = s->h;
  const double* vrow = myb_gj_velocity[j + HALF];
  const double* xrow = myb_gj_position[j + HALF];
  for (size_t i = 0; i < n; ++i) {
    double first = 0.0;
    double second = 0.0;
    sums_at(s, i, j, &first, &second);
    double dv = 0.0;
    double dx = 0.0;
    for (int k = 0; k < WINDOW; ++k) {
      dv += vrow[k] * s->window[k][i];
      dx += xrow[k] * s->window[k][i];
    }
    y[i] = second + h * h * dx;
    y[n + i] = first + h * dv;
  }
}

/* Sums at the centre from the state y0 there: s = v - h sum V a, S = x - h^2 sum X a. */
static void anchor_sums(myb_gauss_jackson* s, const double* y0) {
  const size_t n = s->n;
  const double h = s->h;
  for (size_t i = 0; i < n; ++i) {
    double dv = 0.0;
    double dx = 0.0;
    for (int k = 0; k < WINDOW; ++k) {
      dv += myb_gj_velocity[HALF][k] * s->window[k][i];
      dx += myb_gj_position[HALF][k] * s->window[k][i];
    }
    s->s[i] = myb_eft_two_sum(y0[n + i], -(h * dv), &s->s_comp[i]);
    s->big_s[i] = myb_eft_two_sum(y0[i], -(h * h * dx), &s->big_s_comp[i]);
  }
}

/* max |a - b| / max |a| over n components (0 when both vanish). */
static double relative_change(size_t n, const double* a, const double* b) {
  double change = 0.0;
  double scale = 0.0;
  for (size_t i = 0; i < n; ++i) {
    change = larger_of(change, fabs(a[i] - b[i]));
    scale = larger_of(scale, fabs(a[i]));
  }
  if (change != change) {
    return change;
  }
  return scale > 0.0 ? change / scale : (change > 0.0 ? (double)INFINITY : 0.0);
}

myb_status myb_gauss_jackson_start(myb_gauss_jackson* integrator, myb_accel_function a,
                                   void* context, double t0, const double* y0) {
  MYB_CHECK(integrator != NULL && a != NULL && y0 != NULL);
  myb_gauss_jackson* s = integrator;
  const size_t n = s->n;
  if (!finite(t0)) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < 2 * n; ++i) {
    if (!finite(y0[i])) {
      return MYB_E_DOMAIN;
    }
  }
  s->f = a;
  s->context = context;
  s->t0 = t0;
  s->centre = 0;
  s->current = 0;
  s->started = 0;
  s->stats = (myb_gauss_jackson_stats){};
  myb_status status = evaluate(s, 0, y0, s->window[HALF]);
  if (status != MYB_OK) {
    return status;
  }
  /* Taylor guess x0 + v0 (jh) + a0 (jh)^2 / 2 at the other nodes. */
  const double h = s->h;
  for (int j = -HALF; j <= HALF; ++j) {
    if (j == 0) {
      continue;
    }
    const double dt = (double)j * h;
    for (size_t i = 0; i < n; ++i) {
      s->work[i] = y0[i] + dt * (y0[n + i] + 0.5 * dt * s->window[HALF][i]);
      s->work[n + i] = y0[n + i] + dt * s->window[HALF][i];
    }
    status = evaluate(s, j, s->work, s->window[HALF + j]);
    if (status != MYB_OK) {
      return status;
    }
  }
  /* Iterate the formulas at the nodes until the accelerations settle. */
  double change_last = (double)INFINITY;
  /* bound: START_ITERATIONS */
  for (int iteration = 0; iteration < START_ITERATIONS; ++iteration) {
    anchor_sums(s, y0);
    double change = 0.0;
    for (int j = -HALF; j <= HALF; ++j) {
      if (j == 0) {
        continue;
      }
      state_at_node(s, j, s->work);
      memcpy(s->previous, s->window[HALF + j], n * sizeof(double));
      status = evaluate(s, j, s->work, s->window[HALF + j]);
      if (status != MYB_OK) {
        return status;
      }
      change = larger_of(change, relative_change(n, s->window[HALF + j], s->previous));
    }
    if (change != change) {
      return MYB_E_NOCONV;
    }
    if (change <= s->tolerance ||
        (iteration >= 3 && !(change < change_last) && change <= stall_limit)) {
      anchor_sums(s, y0);
      memcpy(s->y, y0, 2 * n * sizeof(double));
      s->started = 1;
      return MYB_OK;
    }
    change_last = change;
  }
  return MYB_E_NOCONV; /* the step is too long for the start iteration to converge */
}

/* Moves the centre one node forward: the sums absorb a(centre), the window rotates. */
static void shift(myb_gauss_jackson* s) {
  const double h = s->h;
  for (size_t i = 0; i < s->n; ++i) {
    s->s[i] = myb_eft_two_sum(s->s[i], h * s->window[HALF][i] + s->s_comp[i], &s->s_comp[i]);
    s->big_s[i] = myb_eft_two_sum(s->big_s[i], h * (s->s[i] + s->s_comp[i]) + s->big_s_comp[i],
                                  &s->big_s_comp[i]);
  }
  double* oldest = s->window[0];
  for (int k = 0; k < WINDOW; ++k) {
    s->window[k] = s->window[k + 1];
  }
  s->window[WINDOW] = oldest;
  s->centre++;
}

/* Predicts and corrects the node centre + 5; the centre moves to centre + 1. */
static myb_status advance(myb_gauss_jackson* s) {
  const size_t n = s->n;
  state_at_node(s, HALF + 1, s->work); /* predictor */
  const int64_t index = s->centre + HALF + 1;
  myb_status status = evaluate(s, index, s->work, s->window[WINDOW]);
  if (status != MYB_OK) {
    return status;
  }
  memcpy(s->scratch, s->work, n * sizeof(double)); /* predicted positions */
  shift(s);
  double change_last = (double)INFINITY;
  int converged = 0;
  /* bound: max_corrections */
  for (int iteration = 0; iteration < s->max_corrections; ++iteration) {
    state_at_node(s, HALF, s->work); /* corrector */
    memcpy(s->previous, s->window[WINDOW - 1], n * sizeof(double));
    status = evaluate(s, index, s->work, s->window[WINDOW - 1]);
    if (status != MYB_OK) {
      return status;
    }
    const double change = relative_change(n, s->window[WINDOW - 1], s->previous);
    if (change != change) {
      return MYB_E_STEP_MIN; /* the solution left the doubles */
    }
    if (change <= s->tolerance || (iteration >= 2 && !(change < change_last))) {
      converged = change <= stall_limit;
      break;
    }
    change_last = change;
  }
  if (!converged) {
    s->stats.unconverged++;
  }
  /* Predictor-corrector difference of the positions: an estimate of the local error. */
  state_at_node(s, HALF, s->work);
  double difference = 0.0;
  double scale = 0.0;
  for (size_t i = 0; i < n; ++i) {
    difference = larger_of(difference, fabs(s->work[i] - s->scratch[i]));
    scale = larger_of(scale, fabs(s->work[i]));
  }
  if (scale > 0.0) {
    s->stats.max_difference = larger_of(s->stats.max_difference, difference / scale);
  }
  return MYB_OK;
}

myb_status myb_gauss_jackson_step(myb_gauss_jackson* integrator) {
  MYB_CHECK(integrator != NULL && integrator->started);
  myb_gauss_jackson* s = integrator;
  if (s->current == s->centre + HALF) {
    const myb_status status = advance(s);
    if (status != MYB_OK) {
      return status;
    }
  }
  s->current++;
  state_at_node(s, (int)(s->current - s->centre), s->y);
  s->stats.steps++;
  return MYB_OK;
}

myb_status myb_gauss_jackson_integrate(myb_gauss_jackson* integrator, double t_end) {
  MYB_CHECK(integrator != NULL && integrator->started);
  myb_gauss_jackson* s = integrator;
  if (!finite(t_end)) {
    return MYB_E_DOMAIN;
  }
  const double nodes = (t_end - node_time(s, s->current)) / s->h;
  if (!(nodes > -1e-9)) {
    return MYB_E_DOMAIN; /* behind the current node: the grid runs one way */
  }
  if (!(nodes < 0x1p62)) {
    return MYB_E_DOMAIN;
  }
  /* bound: the number of nodes to t_end, at most 2^62 */
  const int64_t target = s->current + (int64_t)ceil(nodes - 1e-9);
  for (int64_t i = s->current; i < target; ++i) {
    const myb_status status = myb_gauss_jackson_step(s);
    if (status != MYB_OK) {
      return status;
    }
  }
  return MYB_OK;
}

double myb_gauss_jackson_time(const myb_gauss_jackson* integrator) {
  MYB_CHECK(integrator != NULL);
  return node_time(integrator, integrator->current);
}

const double* myb_gauss_jackson_state(const myb_gauss_jackson* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->y;
}

myb_gauss_jackson_stats myb_gauss_jackson_statistics(const myb_gauss_jackson* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->stats;
}

myb_status myb_gauss_jackson_dense(myb_gauss_jackson* integrator, double t, double* y) {
  MYB_CHECK(integrator != NULL && y != NULL);
  const myb_gauss_jackson* s = integrator;
  if (!s->started) {
    return MYB_E_RANGE;
  }
  /* Position in steps from the centre; the rounding of t at the window ends is tolerated. */
  double u = (t - node_time(s, s->centre)) / s->h;
  if (!(u >= -(double)HALF - 1e-9 && u <= (double)HALF + 1e-9)) {
    return MYB_E_RANGE;
  }
  u = u > (double)HALF ? (double)HALF : (u < -(double)HALF ? -(double)HALF : u);
  /* x(u) = x_c + u h v_c + h^2 sum_k a_k I2_k(u), v(u) = v_c + h sum_k a_k I1_k(u), where
   * I1_k = int_0^u L_k, I2_k = int_0^u (u - w) L_k(w) dw: exact for a of degree 8. */
  const size_t n = s->n;
  const double h = s->h;
  double i1[WINDOW];
  double i2[WINDOW];
  for (int k = 0; k < WINDOW; ++k) {
    double p1 = 0.0;
    double p2 = 0.0;
    for (int m = WINDOW - 1; m >= 0; --m) {
      p1 = p1 * u + myb_gj_lagrange[k][m] / (double)(m + 1);
      p2 = p2 * u + myb_gj_lagrange[k][m] / (double)((m + 1) * (m + 2));
    }
    i1[k] = p1 * u;
    i2[k] = p2 * u * u;
  }
  for (size_t i = 0; i < n; ++i) {
    double first = 0.0;
    double second = 0.0;
    sums_at(s, i, 0, &first, &second);
    double dv = 0.0;
    double dx = 0.0;
    double iv = 0.0;
    double ix = 0.0;
    for (int k = 0; k < WINDOW; ++k) {
      const double a = s->window[k][i];
      dv += myb_gj_velocity[HALF][k] * a;
      dx += myb_gj_position[HALF][k] * a;
      iv += i1[k] * a;
      ix += i2[k] * a;
    }
    const double vc = first + h * dv;
    const double xc = second + h * h * dx;
    y[i] = xc + (u * h * vc + h * h * ix);
    y[n + i] = vc + h * iv;
  }
  return MYB_OK;
}

static myb_status dense_view(void* integrator, double t, double* y) {
  return myb_gauss_jackson_dense(integrator, t, y);
}

myb_status myb_gauss_jackson_locate(myb_gauss_jackson* integrator, myb_event_function g,
                                    void* context, int* found, double* t_event,
                                    double* y_event) {
  MYB_CHECK(integrator != NULL && g != NULL && found != NULL && t_event != NULL);
  myb_gauss_jackson* s = integrator;
  *found = 0;
  if (!s->started || s->current == 0) {
    return MYB_E_RANGE;
  }
  const double t1 = node_time(s, s->current);
  const double t0 = node_time(s, s->current - 1);
  double* y0 = s->work;
  const myb_status status = myb_gauss_jackson_dense(s, t0, y0);
  if (status != MYB_OK) {
    return status;
  }
  const myb_internal_step_view view = {.dense = dense_view,
                                       .integrator = s,
                                       .m = 2 * s->n,
                                       .t0 = t0,
                                       .y0 = y0,
                                       .t1 = t1,
                                       .y1 = s->y,
                                       .scratch = s->scratch};
  return myb_internal_locate_event(&view, g, context, found, t_event, y_event);
}
