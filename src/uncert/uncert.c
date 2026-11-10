/* SPDX-License-Identifier: BSD-3-Clause */
/* Variational equations, covariances, deterministic random numbers. */
#include <math.h>
#include <string.h>

#include "integ/ias15_internal.h"
#include "uncert/variational_internal.h"
#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "myb/num.h"
#include "myb/uncert.h"

struct myb_variational {
  myb_allocator allocator;
  myb_force_model* model;
  size_t np;
  size_t columns; /* 6 + np */
  size_t n;       /* coordinates: 3 (1 + columns) */
  myb_ias15* ias;
  double* scratch; /* 2 n state + dadp */
};

/* x = (state, dx_0 .. dx_{columns-1}), likewise v. */
static myb_status variational_accel(void* context, double t, const double* x, const double* v,
                                    double* a) {
  myb_variational* var = context;
  myb_status st = myb_force_model_acceleration(var->model, t, x, v, a);
  if (st != MYB_OK) {
    return st;
  }
  double jx[9];
  double jv[9];
  double jp[3 * MYB_MODEL_MAX_PARAMETERS];
  st = myb_force_model_jacobian(var->model, t, x, v, jx, jv, var->np > 0 ? jp : NULL);
  if (st != MYB_OK) {
    return st;
  }
  for (size_t c = 0; c < var->columns; ++c) {
    const double* dx = x + 3 + 3 * c;
    const double* dv = v + 3 + 3 * c;
    double* da = a + 3 + 3 * c;
    for (int i = 0; i < 3; ++i) {
      double s = jx[3 * i] * dx[0] + jx[3 * i + 1] * dx[1] + jx[3 * i + 2] * dx[2] +
                 jv[3 * i] * dv[0] + jv[3 * i + 1] * dv[1] + jv[3 * i + 2] * dv[2];
      if (c >= 6) {
        s += jp[(size_t)i * var->np + (c - 6)];
      }
      da[i] = s;
    }
  }
  return MYB_OK;
}

void myb_variational_destroy(myb_variational* variational) {
  if (variational == NULL) {
    return;
  }
  const myb_allocator a = variational->allocator;
  myb_ias15_destroy(variational->ias);
  myb_internal_release(&a, variational->scratch);
  myb_internal_release(&a, variational);
}

myb_status myb_variational_create(myb_force_model* model, const myb_variational_options* options,
                                  const myb_allocator* allocator,
                                  myb_variational** variational) {
  MYB_CHECK(model != NULL && variational != NULL);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_variational), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_variational* var = block;
  memset(var, 0, sizeof *var);
  var->allocator = a;
  var->model = model;
  var->np = myb_force_model_parameter_count(model);
  var->columns = 6 + var->np;
  var->n = 3 * (1 + var->columns);
  myb_ias15_options o = {.epsilon = 1e-10, .block = 3};
  if (options != NULL) {
    if (options->epsilon < 0.0 || options->h_max < 0.0 || options->max_steps < 0) {
      myb_variational_destroy(var);
      return MYB_E_DOMAIN;
    }
    if (options->epsilon > 0.0) {
      o.epsilon = options->epsilon;
    }
    o.h_max = options->h_max;
    o.max_steps = options->max_steps;
  }
  status = myb_ias15_create(var->n, &o, allocator, &var->ias);
  if (status == MYB_OK) {
    status = myb_internal_alloc_array(&a, 2 * var->n, sizeof(double), &block);
    if (status == MYB_OK) {
      var->scratch = block;
    }
  }
  if (status != MYB_OK) {
    myb_variational_destroy(var);
    return status;
  }
  *variational = var;
  return MYB_OK;
}

myb_status myb_variational_start(myb_variational* variational, double t0, const double y0[6]) {
  MYB_CHECK(variational != NULL && y0 != NULL);
  myb_variational* var = variational;
  double* y = var->scratch;
  memset(y, 0, 2 * var->n * sizeof(double));
  double* x = y;
  double* v = y + var->n;
  for (int i = 0; i < 3; ++i) {
    x[i] = y0[i];
    v[i] = y0[3 + i];
  }
  for (size_t c = 0; c < 6; ++c) {
    if (c < 3) {
      x[3 + 3 * c + c] = 1.0;
    } else {
      v[3 + 3 * c + (c - 3)] = 1.0;
    }
  }
  return myb_ias15_start(var->ias, variational_accel, var, t0, y);
}

myb_status myb_variational_step(myb_variational* variational, double t_limit) {
  MYB_CHECK(variational != NULL);
  return myb_ias15_step(variational->ias, t_limit);
}

myb_status myb_internal_variational_restart(myb_variational* variational, double t,
                                            const double y[6], const double phi[36],
                                            const double* sensitivity) {
  MYB_CHECK(variational != NULL && y != NULL && phi != NULL &&
            (sensitivity != NULL || variational->np == 0));
  myb_variational* var = variational;
  double* full = var->scratch;
  memset(full, 0, 2 * var->n * sizeof(double));
  double* x = full;
  double* v = full + var->n;
  for (int i = 0; i < 3; ++i) {
    x[i] = y[i];
    v[i] = y[3 + i];
  }
  for (size_t c = 0; c < var->columns; ++c) {
    double* dx = x + 3 + 3 * c;
    double* dv = v + 3 + 3 * c;
    for (int i = 0; i < 3; ++i) {
      if (c < 6) {
        dx[i] = phi[(size_t)i * 6 + c];
        dv[i] = phi[(size_t)(3 + i) * 6 + c];
      } else {
        dx[i] = sensitivity[(size_t)i * var->np + (c - 6)];
        dv[i] = sensitivity[(size_t)(3 + i) * var->np + (c - 6)];
      }
    }
  }
  /* the step the integration was making goes on: the automatic first step is made for
   * Phi = 1, S = 0 */
  const double planned = myb_internal_ias15_planned_step(var->ias);
  const myb_status status = myb_ias15_start(var->ias, variational_accel, var, t, full);
  if (status == MYB_OK && planned > 0.0 && isfinite(planned)) {
    myb_internal_ias15_plan_step(var->ias, planned);
  }
  return status;
}

myb_status myb_internal_variational_retake(myb_variational* variational, double t_limit) {
  MYB_CHECK(variational != NULL);
  return myb_internal_ias15_retake(variational->ias, t_limit);
}

myb_status myb_variational_integrate(myb_variational* variational, double t_end) {
  MYB_CHECK(variational != NULL);
  return myb_ias15_integrate(variational->ias, t_end);
}

double myb_variational_time(const myb_variational* variational) {
  MYB_CHECK(variational != NULL);
  return myb_ias15_time(variational->ias);
}

static void unpack(const myb_variational* var, const double* full, double y[6], double phi[36],
                   double* sensitivity) {
  const double* x = full;
  const double* v = full + var->n;
  if (y != NULL) {
    for (int i = 0; i < 3; ++i) {
      y[i] = x[i];
      y[3 + i] = v[i];
    }
  }
  for (size_t c = 0; c < var->columns; ++c) {
    const double* dx = x + 3 + 3 * c;
    const double* dv = v + 3 + 3 * c;
    for (int i = 0; i < 3; ++i) {
      if (c < 6 && phi != NULL) {
        phi[(size_t)i * 6 + c] = dx[i];
        phi[(size_t)(3 + i) * 6 + c] = dv[i];
      } else if (c >= 6 && sensitivity != NULL) {
        sensitivity[(size_t)i * var->np + (c - 6)] = dx[i];
        sensitivity[(size_t)(3 + i) * var->np + (c - 6)] = dv[i];
      }
    }
  }
}

void myb_variational_get(const myb_variational* variational, double y[6], double phi[36],
                         double* sensitivity) {
  MYB_CHECK(variational != NULL);
  unpack(variational, myb_ias15_state(variational->ias), y, phi, sensitivity);
}

myb_status myb_variational_dense(myb_variational* variational, double t, double y[6],
                                 double phi[36], double* sensitivity) {
  MYB_CHECK(variational != NULL);
  const myb_status st = myb_ias15_dense(variational->ias, t, variational->scratch);
  if (st != MYB_OK) {
    return st;
  }
  unpack(variational, variational->scratch, y, phi, sensitivity);
  return MYB_OK;
}

/* ---- Covariances ------------------------------------------------------------------------------ */

void myb_covariance_transform(size_t n, const double* m, const double* p, double* out) {
  MYB_CHECK(m != NULL && p != NULL && out != NULL && n <= 64);
  /* row by row: t = (M P)_i, then out_ij = t . M_j; one row of n doubles of scratch */
  double t[64];
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double s = 0.0;
      for (size_t k = 0; k < n; ++k) {
        s += m[i * n + k] * p[k * n + j];
      }
      t[j] = s;
    }
    for (size_t j = 0; j < n; ++j) {
      double s = 0.0;
      for (size_t k = 0; k < n; ++k) {
        s += t[k] * m[j * n + k];
      }
      out[i * n + j] = s;
    }
  }
}

myb_status myb_covariance_rtn(const double r[3], const double v[3], const double p[36],
                              double out[36]) {
  MYB_CHECK(r != NULL && v != NULL && p != NULL && out != NULL);
  const double rr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  const double h[3] = {r[1] * v[2] - r[2] * v[1], r[2] * v[0] - r[0] * v[2],
                       r[0] * v[1] - r[1] * v[0]};
  const double hh = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
  if (!(rr > 0.0) || !(hh > 0.0) || !isfinite(rr) || !isfinite(hh)) {
    return MYB_E_DOMAIN;
  }
  double rot[9];
  for (int i = 0; i < 3; ++i) {
    rot[i] = r[i] / rr;
    rot[6 + i] = h[i] / hh;
  }
  rot[3] = rot[7] * rot[2] - rot[8] * rot[1];
  rot[4] = rot[8] * rot[0] - rot[6] * rot[2];
  rot[5] = rot[6] * rot[1] - rot[7] * rot[0];
  double m[36];
  memset(m, 0, sizeof m);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      m[i * 6 + j] = rot[3 * i + j];
      m[(3 + i) * 6 + 3 + j] = rot[3 * i + j];
    }
  }
  myb_covariance_transform(6, m, p, out);
  return MYB_OK;
}

/* ---- Random numbers ---------------------------------------------------------------------------- */

static uint64_t splitmix64(uint64_t* x) {
  uint64_t z = (*x += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}

static uint64_t rotl(uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

void myb_rng_seed(myb_rng* rng, uint64_t seed) {
  MYB_CHECK(rng != NULL);
  uint64_t x = seed;
  for (int i = 0; i < 4; ++i) {
    rng->s[i] = splitmix64(&x);
  }
  rng->has_spare = 0;
  rng->spare = 0.0;
}

uint64_t myb_rng_next(myb_rng* rng) {
  MYB_CHECK(rng != NULL);
  uint64_t* s = rng->s;
  const uint64_t result = rotl(s[1] * 5, 7) * 9;
  const uint64_t t = s[1] << 17;
  s[2] ^= s[0];
  s[3] ^= s[1];
  s[1] ^= s[2];
  s[0] ^= s[3];
  s[2] ^= t;
  s[3] = rotl(s[3], 45);
  return result;
}

double myb_rng_uniform(myb_rng* rng) {
  return (double)(myb_rng_next(rng) >> 11) * 0x1p-53;
}

double myb_rng_normal(myb_rng* rng) {
  MYB_CHECK(rng != NULL);
  if (rng->has_spare) {
    rng->has_spare = 0;
    return rng->spare;
  }
  /* u1 in (0, 1]: log never sees 0 */
  const double u1 = 1.0 - myb_rng_uniform(rng);
  const double u2 = myb_rng_uniform(rng);
  const double radius = sqrt(-2.0 * myb_log(u1));
  const double angle = 6.283185307179586 * u2;
  rng->spare = radius * myb_sin(angle);
  rng->has_spare = 1;
  return radius * myb_cos(angle);
}

void myb_rng_gaussian(myb_rng* rng, size_t n, const double* mean, const double* l, double* out) {
  MYB_CHECK(rng != NULL && mean != NULL && l != NULL && out != NULL && n <= 64);
  double z[64];
  for (size_t i = 0; i < n; ++i) {
    z[i] = myb_rng_normal(rng);
  }
  for (size_t i = 0; i < n; ++i) {
    double s = mean[i];
    for (size_t k = 0; k <= i; ++k) {
      s += l[i * n + k] * z[k];
    }
    out[i] = s;
  }
}

myb_status myb_sample_statistics(size_t count, size_t n, const double* samples, double* mean,
                                 double* covariance) {
  MYB_CHECK(samples != NULL && mean != NULL && n <= 64);
  if (count < 2) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < n; ++i) {
    myb_sum_acc s = {0};
    for (size_t k = 0; k < count; ++k) {
      myb_sum_acc_add(&s, samples[k * n + i]);
    }
    mean[i] = myb_sum_acc_value(&s) / (double)count;
  }
  if (covariance != NULL) {
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j <= i; ++j) {
        myb_sum_acc s = {0};
        for (size_t k = 0; k < count; ++k) {
          myb_sum_acc_add(&s, (samples[k * n + i] - mean[i]) * (samples[k * n + j] - mean[j]));
        }
        covariance[i * n + j] = myb_sum_acc_value(&s) / (double)(count - 1);
        covariance[j * n + i] = covariance[i * n + j];
      }
    }
  }
  return MYB_OK;
}
