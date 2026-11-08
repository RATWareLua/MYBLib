/* SPDX-License-Identifier: BSD-3-Clause */
/* Symmetric compositions of the leapfrog for x'' = a(t, x) (H. Yoshida, Phys. Lett. A 150, 262,
 * 1990). One stage with weight c is the kick-drift-kick leapfrog of length c h:
 *   v += (c h / 2) a(t_s, x);  x += c h v;  t_s += c h;  v += (c h / 2) a(t_s, x).
 * The acceleration after a drift serves the next kick as well (first same as last), so a step
 * costs one evaluation per stage. Time enters as a coordinate advanced by the drifts (extended
 * phase space), which keeps the map symplectic for a time-dependent potential too. Positions,
 * velocities and time are updated with compensated summation. */
#include <math.h> /* fabs, ceil */
#include <stdckdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "integ/symplectic_coefficients.h"
#include "myb/integ.h"
#include "num/eft.h"

struct myb_symplectic {
  myb_allocator allocator;
  size_t n;
  int stages;
  const double* weights;
  myb_position_accel_function f;
  void* context;
  double t;
  double t_comp;
  int started;
  int a_valid;
  myb_symplectic_stats stats;
  double* y;      /* 2n: x, v */
  double* y_comp; /* 2n */
  double* a;      /* n: a(t_s, x) after the last drift */
};

static int finite(double x) {
  return x - x == 0.0;
}

myb_status myb_symplectic_create(size_t n, int order, const myb_allocator* allocator,
                                 myb_symplectic** integrator) {
  MYB_CHECK(integrator != NULL);
  const double* weights = NULL;
  int stages = 0;
  switch (order) {
    case 2:
      weights = myb_symplectic_order2;
      stages = 1;
      break;
    case 4:
      weights = myb_symplectic_order4;
      stages = 3;
      break;
    case 6:
      weights = myb_symplectic_order6;
      stages = 9;
      break;
    case 8:
      weights = myb_symplectic_order8;
      stages = 27;
      break;
    default:
      return MYB_E_DOMAIN;
  }
  if (n == 0) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_symplectic), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_symplectic* s = block;
  *s = (myb_symplectic){.allocator = a, .n = n, .stages = stages, .weights = weights};
  size_t count = 0;
  if (ckd_mul(&count, (size_t)5, n)) {
    status = MYB_E_OVERFLOW;
    goto fail;
  }
  status = myb_internal_alloc_array(&a, count, sizeof(double), &block);
  if (status != MYB_OK) {
    goto fail;
  }
  double* p = block;
  memset(p, 0, count * sizeof(double));
  s->y = p;
  s->y_comp = p + 2 * n;
  s->a = p + 4 * n;
  *integrator = s;
  return MYB_OK;

fail:
  myb_internal_release(&a, s);
  return status;
}

void myb_symplectic_destroy(myb_symplectic* integrator) {
  if (integrator == NULL) {
    return;
  }
  const myb_allocator a = integrator->allocator;
  myb_internal_release(&a, integrator->y);
  myb_internal_release(&a, integrator);
}

static myb_status evaluate(myb_symplectic* s, double t) {
  s->stats.evaluations++;
  return s->f(s->context, t, s->y, s->a);
}

myb_status myb_symplectic_start(myb_symplectic* integrator, myb_position_accel_function a,
                                void* context, double t0, const double* y0) {
  MYB_CHECK(integrator != NULL && a != NULL && y0 != NULL);
  myb_symplectic* s = integrator;
  if (!finite(t0)) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < 2 * s->n; ++i) {
    if (!finite(y0[i])) {
      return MYB_E_DOMAIN;
    }
    s->y[i] = y0[i];
    s->y_comp[i] = 0.0;
  }
  s->f = a;
  s->context = context;
  s->t = t0;
  s->t_comp = 0.0;
  s->stats = (myb_symplectic_stats){};
  const myb_status status = evaluate(s, t0);
  s->a_valid = status == MYB_OK;
  s->started = status == MYB_OK;
  return status;
}

/* y[i] += delta with the compensation carried in comp[i]. */
static void add_compensated(double* y, double* comp, size_t i, double delta) {
  y[i] = myb_eft_two_sum(y[i], delta + comp[i], &comp[i]);
}

myb_status myb_symplectic_step(myb_symplectic* integrator, double h) {
  MYB_CHECK(integrator != NULL && integrator->started);
  myb_symplectic* s = integrator;
  if (!finite(h) || h == 0.0) {
    return MYB_E_DOMAIN;
  }
  const size_t n = s->n;
  double* x = s->y;
  double* v = s->y + n;
  double* xc = s->y_comp;
  double* vc = s->y_comp + n;
  if (!s->a_valid) {
    const myb_status status = evaluate(s, s->t + s->t_comp);
    if (status != MYB_OK) {
      return status;
    }
    s->a_valid = 1;
  }
  double elapsed = 0.0; /* sum of the weights so far, times h at use */
  for (int stage = 0; stage < s->stages; ++stage) {
    const double hc = s->weights[stage] * h;
    for (size_t i = 0; i < n; ++i) {
      add_compensated(v, vc, i, 0.5 * hc * s->a[i]);
    }
    for (size_t i = 0; i < n; ++i) {
      add_compensated(x, xc, i, hc * v[i]);
    }
    elapsed += s->weights[stage];
    s->a_valid = 0;
    const myb_status status = evaluate(s, s->t + (s->t_comp + elapsed * h));
    if (status != MYB_OK) {
      return status; /* the state is mid-step: restart with myb_symplectic_start */
    }
    s->a_valid = 1;
    for (size_t i = 0; i < n; ++i) {
      add_compensated(v, vc, i, 0.5 * hc * s->a[i]);
    }
  }
  s->t = myb_eft_two_sum(s->t, h + s->t_comp, &s->t_comp);
  s->stats.steps++;
  return MYB_OK;
}

myb_status myb_symplectic_integrate(myb_symplectic* integrator, double t_end, double h_max) {
  MYB_CHECK(integrator != NULL && integrator->started);
  myb_symplectic* s = integrator;
  if (!finite(t_end) || !(h_max > 0.0) || !finite(h_max)) {
    return MYB_E_DOMAIN;
  }
  const double span = (t_end - s->t) - s->t_comp;
  if (span == 0.0) {
    return MYB_OK;
  }
  const double steps = ceil(fabs(span) / h_max);
  if (!(steps <= 1e15)) {
    return MYB_E_DOMAIN; /* more steps than a counter of doubles can resolve */
  }
  const double h = span / steps;
  /* bound: steps <= 1e15, checked above */
  for (int64_t i = 0; i < (int64_t)steps; ++i) {
    const myb_status status = myb_symplectic_step(s, h);
    if (status != MYB_OK) {
      return status;
    }
  }
  s->t = t_end; /* the equal steps add up to the span up to rounding: land exactly */
  s->t_comp = 0.0;
  return MYB_OK;
}

double myb_symplectic_time(const myb_symplectic* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->t + integrator->t_comp;
}

const double* myb_symplectic_state(const myb_symplectic* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->y;
}

myb_symplectic_stats myb_symplectic_statistics(const myb_symplectic* integrator) {
  MYB_CHECK(integrator != NULL);
  return integrator->stats;
}
