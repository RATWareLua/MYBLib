/* SPDX-License-Identifier: BSD-3-Clause */
/* The budget of effects: each term switched against the reference propagation. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/budget.h"
#include "myb/integ.h"
#include "myb/num.h"

enum { MAX_SAMPLES = 10000 };

/* States (position, velocity) at the sample times t0 + (k + 1)(t1 - t0) / samples. */
static myb_status propagate(myb_force_model* model, double t0, const double y0[6], double t1,
                            double epsilon, size_t samples, double* states) {
  myb_ias15* in = NULL;
  const myb_ias15_options o = {.epsilon = epsilon};
  myb_status st = myb_ias15_create(3, &o, NULL, &in);
  if (st != MYB_OK) {
    return st;
  }
  st = myb_ias15_start(in, myb_force_model_accel, model, t0, y0);
  for (size_t k = 0; k < samples && st == MYB_OK; ++k) {
    const double t = t0 + (t1 - t0) * (double)(k + 1) / (double)samples;
    st = myb_ias15_integrate(in, t);
    if (st == MYB_OK) {
      const double* y = myb_ias15_state(in);
      memcpy(states + 6 * k, y, 3 * sizeof(double));
      memcpy(states + 6 * k + 3, y + 3, 3 * sizeof(double));
    }
  }
  myb_ias15_destroy(in);
  return st;
}

myb_status myb_budget(myb_force_model* model, double t0, const double y0[6], double t1,
                      const myb_budget_options* options, myb_budget_entry* entries) {
  MYB_CHECK(model != NULL && y0 != NULL && options != NULL && entries != NULL);
  const size_t samples = options->samples > 0 ? options->samples : 100;
  const double epsilon = options->epsilon > 0.0 ? options->epsilon : 1e-11;
  if (!(options->tolerance > 0.0) || options->epsilon < 0.0 || samples > MAX_SAMPLES ||
      !(t1 != t0) || !isfinite(t0) || !isfinite(t1)) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(NULL);
  void* block = NULL;
  myb_status st = myb_internal_alloc_array(&a, 12 * samples, sizeof(double), &block);
  if (st != MYB_OK) {
    return st;
  }
  double* reference = block;
  double* switched = reference + 6 * samples;
  st = propagate(model, t0, y0, t1, epsilon, samples, reference);
  const size_t count = myb_force_model_count(model);
  for (size_t k = 0; k < count && st == MYB_OK; ++k) {
    myb_budget_entry* e = &entries[k];
    memset(e, 0, sizeof *e);
    const char* name = myb_force_model_name(model, k);
    memcpy(e->name, name, strlen(name) + 1);
    e->enabled = myb_force_model_enabled(model, k);
    /* the largest acceleration of the term along the reference trajectory */
    for (size_t s = 0; s < samples && st == MYB_OK; ++s) {
      double acc[3];
      const double t = t0 + (t1 - t0) * (double)(s + 1) / (double)samples;
      st = myb_force_model_term(model, k, t, reference + 6 * s, reference + 6 * s + 3, acc);
      if (st == MYB_OK) {
        const double m = sqrt(acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]);
        e->max_acceleration = fmax(e->max_acceleration, m);
      }
    }
    if (st != MYB_OK) {
      break;
    }
    myb_force_model_set_enabled(model, k, !e->enabled);
    st = propagate(model, t0, y0, t1, epsilon, samples, switched);
    myb_force_model_set_enabled(model, k, e->enabled);
    if (st != MYB_OK) {
      break;
    }
    for (size_t s = 0; s < samples; ++s) {
      const double* p = reference + 6 * s;
      const double* q = switched + 6 * s;
      const double d = sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) +
                            (p[2] - q[2]) * (p[2] - q[2]));
      e->effect = fmax(e->effect, d);
    }
    if (!e->enabled && e->effect > options->tolerance) {
      e->flag = MYB_BUDGET_MISSING;
    } else if (e->enabled && e->effect < options->tolerance) {
      e->flag = MYB_BUDGET_NEGLIGIBLE;
    }
  }
  myb_internal_release(&a, reference);
  return st;
}

/* ---- The linear estimate ---------------------------------------------------------------------- */

typedef struct linear_context {
  myb_force_model* model;
  size_t count;
} linear_context;

/* x = (state, d_0 .. d_{count-1}), likewise v. The acceleration of the state is the sum of the
 * enabled terms taken as the model takes it, each term being computed once. */
static myb_status linear_accel(void* context, double t, const double* x, const double* v,
                               double* a) {
  const linear_context* c = context;
  double jx[9];
  double jv[9];
  myb_status st = myb_force_model_jacobian(c->model, t, x, v, jx, jv, NULL);
  if (st != MYB_OK) {
    return st;
  }
  myb_sum_acc sum[3] = {{0}};
  for (size_t k = 0; k < c->count; ++k) {
    const double* dx = x + 3 + 3 * k;
    const double* dv = v + 3 + 3 * k;
    double* da = a + 3 + 3 * k;
    double ak[3];
    st = myb_force_model_term(c->model, k, t, x, v, ak);
    if (st != MYB_OK) {
      return st;
    }
    const int enabled = myb_force_model_enabled(c->model, k);
    for (int i = 0; i < 3; ++i) {
      if (enabled) {
        myb_sum_acc_add(&sum[i], ak[i]);
      }
      da[i] = ak[i] + jx[3 * i] * dx[0] + jx[3 * i + 1] * dx[1] + jx[3 * i + 2] * dx[2] +
              jv[3 * i] * dv[0] + jv[3 * i + 1] * dv[1] + jv[3 * i + 2] * dv[2];
    }
  }
  for (int i = 0; i < 3; ++i) {
    a[i] = myb_sum_acc_value(&sum[i]);
  }
  return MYB_OK;
}

myb_status myb_budget_linear(myb_force_model* model, double t0, const double y0[6], double t1,
                             const myb_budget_options* options, myb_budget_entry* entries) {
  MYB_CHECK(model != NULL && y0 != NULL && options != NULL && entries != NULL);
  const size_t samples = options->samples > 0 ? options->samples : 100;
  const double epsilon = options->epsilon > 0.0 ? options->epsilon : 1e-11;
  if (!(options->tolerance > 0.0) || options->epsilon < 0.0 || samples > MAX_SAMPLES ||
      !(t1 != t0) || !isfinite(t0) || !isfinite(t1)) {
    return MYB_E_DOMAIN;
  }
  const size_t count = myb_force_model_count(model);
  const size_t n = 3 * (1 + count);
  const myb_allocator a = myb_internal_allocator(NULL);
  void* block = NULL;
  myb_status st = myb_internal_alloc_array(&a, 2 * n, sizeof(double), &block);
  if (st != MYB_OK) {
    return st;
  }
  double* start = block;
  memset(start, 0, 2 * n * sizeof(double));
  memcpy(start, y0, 3 * sizeof(double));
  memcpy(start + n, y0 + 3, 3 * sizeof(double));
  for (size_t k = 0; k < count; ++k) {
    myb_budget_entry* e = &entries[k];
    memset(e, 0, sizeof *e);
    const char* name = myb_force_model_name(model, k);
    memcpy(e->name, name, strlen(name) + 1);
    e->enabled = myb_force_model_enabled(model, k);
  }
  linear_context context = {.model = model, .count = count};
  myb_ias15* in = NULL;
  const myb_ias15_options o = {.epsilon = epsilon, .block = 3};
  st = myb_ias15_create(n, &o, NULL, &in);
  if (st == MYB_OK) {
    st = myb_ias15_start(in, linear_accel, &context, t0, start);
  }
  for (size_t s = 0; s < samples && st == MYB_OK; ++s) {
    const double t = t0 + (t1 - t0) * (double)(s + 1) / (double)samples;
    st = myb_ias15_integrate(in, t);
    if (st != MYB_OK) {
      break;
    }
    const double* x = myb_ias15_state(in);
    const double* v = x + n;
    for (size_t k = 0; k < count && st == MYB_OK; ++k) {
      myb_budget_entry* e = &entries[k];
      const double* d = x + 3 + 3 * k;
      e->effect = fmax(e->effect, sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
      double acc[3];
      st = myb_force_model_term(model, k, t, x, v, acc);
      if (st == MYB_OK) {
        const double m = sqrt(acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]);
        e->max_acceleration = fmax(e->max_acceleration, m);
      }
    }
  }
  myb_ias15_destroy(in);
  myb_internal_release(&a, start);
  for (size_t k = 0; k < count && st == MYB_OK; ++k) {
    myb_budget_entry* e = &entries[k];
    if (!e->enabled && e->effect > options->tolerance) {
      e->flag = MYB_BUDGET_MISSING;
    } else if (e->enabled && e->effect < options->tolerance) {
      e->flag = MYB_BUDGET_NEGLIGIBLE;
    }
  }
  return st;
}
