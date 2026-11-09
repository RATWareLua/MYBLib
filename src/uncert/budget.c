/* SPDX-License-Identifier: BSD-3-Clause */
/* The budget of effects: each term switched against the reference propagation. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/budget.h"
#include "myb/integ.h"

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
