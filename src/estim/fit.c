/* SPDX-License-Identifier: BSD-3-Clause */
/* Weighted least squares over an arc: Gauss-Newton on the variational equations. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/estim.h"
#include "myb/num.h"
#include "myb/orbit.h"

enum {
  MAX_UNKNOWNS = 6 + MYB_MODEL_MAX_PARAMETERS + MYB_FIT_MAX_BIASES,
  DEFAULT_ITERATIONS = 10
};

static const double default_tolerance = 1e-6;
static const double pi = 3.141592653589793;

typedef struct fit {
  myb_force_model* model;
  myb_orbit* orbit;
  const myb_observation* observations;
  size_t count;
  double t0;
  double t_end;
  size_t np;
  size_t nb;
  double* bias;
  size_t nfree;
  size_t free_index[MAX_UNKNOWNS];
  const double* prior_sigma;
  double prior[MAX_UNKNOWNS]; /* the values the a priori information is about */
  double* normal;             /* nfree x nfree */
  double* rhs;                /* nfree */
  double* scaled;             /* nfree x nfree: the factor of the scaled normal matrix */
  double* scale;              /* nfree */
  double* column;             /* nfree */
} fit;

/* MYB_OK for an observation the estimators can use. */
static myb_status valid_observation(const myb_observation* o, size_t nb) {
  if (o->measure == NULL || o->count < 1 || o->count > MYB_MEASURE_MAX || !isfinite(o->t) ||
      o->bias < 0 || (o->bias > 0 && (size_t)o->bias - 1 + o->count > nb)) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < o->count; ++i) {
    if (!isfinite(o->value[i]) || !(o->sigma[i] > 0.0) || !isfinite(o->sigma[i])) {
      return MYB_E_DOMAIN;
    }
  }
  return MYB_OK;
}

/* observed - computed, an angle reduced to [-pi, pi). */
static double residual_of(const myb_observation* o, size_t i, double computed) {
  double r = o->value[i] - computed;
  if ((o->periodic >> i) & 1u) {
    r -= 2.0 * pi * floor((r + pi) / (2.0 * pi));
  }
  return r;
}

/* One pass over the arc from y0: the normal equations of the free unknowns, the residuals
 * (may be NULL) and the sum of the squares of the weighted residuals. */
static myb_status accumulate(fit* f, const double y0[6], double* residuals, double* chi2,
                             size_t* values) {
  const size_t nf = f->nfree;
  const size_t np = f->np;
  memset(f->normal, 0, nf * nf * sizeof(double));
  memset(f->rhs, 0, nf * sizeof(double));
  *chi2 = 0.0;
  *values = 0;
  myb_status status = myb_orbit_start(f->orbit, f->t0, y0, f->t_end);
  for (size_t k = 0; k < f->count && status == MYB_OK; ++k) {
    const myb_observation* o = &f->observations[k];
    double y[6];
    double phi[36];
    double sens[6 * MYB_MODEL_MAX_PARAMETERS];
    double computed[MYB_MEASURE_MAX];
    double partial[6 * MYB_MEASURE_MAX];
    status = myb_orbit_partials(f->orbit, o->t, y, phi, np > 0 ? sens : NULL);
    if (status == MYB_OK) {
      status = o->measure(o->context, o->t, y, computed, partial);
    }
    for (size_t i = 0; i < o->count && status == MYB_OK; ++i) {
      /* the shift of this value, if it has one, and its place among the unknowns */
      const size_t shift = o->bias > 0 ? (size_t)o->bias - 1 + i : f->nb;
      if (shift < f->nb) {
        computed[i] += f->bias[shift];
      }
      const double r = residual_of(o, i, computed[i]);
      if (!isfinite(r)) {
        status = MYB_E_DOMAIN;
        break;
      }
      if (residuals != NULL) {
        residuals[*values] = r;
      }
      ++*values;
      /* d value / d unknowns = (d value / d y) [Phi | S] */
      double row[MAX_UNKNOWNS];
      const double* h = partial + 6 * i;
      for (size_t j = 0; j < 6; ++j) {
        double sum = 0.0;
        for (size_t m = 0; m < 6; ++m) {
          sum += h[m] * phi[6 * m + j];
        }
        row[j] = sum;
      }
      for (size_t q = 0; q < np; ++q) {
        double sum = 0.0;
        for (size_t m = 0; m < 6; ++m) {
          sum += h[m] * sens[m * np + q];
        }
        row[6 + q] = sum;
      }
      for (size_t q = 0; q < f->nb; ++q) {
        row[6 + np + q] = q == shift ? 1.0 : 0.0;
      }
      const double w = 1.0 / o->sigma[i];
      const double wr = w * r;
      *chi2 += wr * wr;
      for (size_t p = 0; p < nf; ++p) {
        const double ap = w * row[f->free_index[p]];
        f->rhs[p] += ap * wr;
        for (size_t q = 0; q <= p; ++q) {
          f->normal[p * nf + q] += ap * (w * row[f->free_index[q]]);
        }
      }
    }
  }
  if (status != MYB_OK) {
    return status;
  }
  for (size_t p = 0; p < nf; ++p) {
    for (size_t q = 0; q < p; ++q) {
      f->normal[q * nf + p] = f->normal[p * nf + q];
    }
  }
  if (f->prior_sigma != NULL) {
    for (size_t p = 0; p < nf; ++p) {
      const size_t index = f->free_index[p];
      const double sigma = f->prior_sigma[index];
      if (sigma > 0.0) {
        const double now = index < 6 ? y0[index]
                           : index < 6 + np ? *myb_force_model_parameter(f->model, index - 6)
                                            : f->bias[index - 6 - np];
        f->normal[p * nf + p] += 1.0 / (sigma * sigma);
        f->rhs[p] += (f->prior[index] - now) / (sigma * sigma);
      }
    }
  }
  return MYB_OK;
}

/* Factorises the normal matrix scaled by its diagonal (f->scaled, f->scale). */
static myb_status factorise(fit* f) {
  const size_t nf = f->nfree;
  for (size_t p = 0; p < nf; ++p) {
    const double d = f->normal[p * nf + p];
    if (!(d > 0.0) || !isfinite(d)) {
      return MYB_E_DOMAIN;
    }
    f->scale[p] = sqrt(d);
  }
  for (size_t p = 0; p < nf; ++p) {
    for (size_t q = 0; q < nf; ++q) {
      f->scaled[p * nf + q] = f->normal[p * nf + q] / (f->scale[p] * f->scale[q]);
    }
  }
  return myb_cholesky(nf, f->scaled);
}

myb_status myb_fit_orbit(myb_force_model* model, double t0, double y0[6],
                         const myb_observation* observations, size_t count,
                         const myb_fit_options* options, const myb_allocator* allocator,
                         myb_fit_result* result, double* covariance, double* residuals) {
  MYB_CHECK(model != NULL && y0 != NULL && result != NULL &&
            (count == 0 || observations != NULL));
  const myb_fit_options none = {0};
  const myb_fit_options* opt = options != NULL ? options : &none;
  memset(result, 0, sizeof *result);
  if (count == 0 || !isfinite(t0) || opt->max_iterations < 0 || !(opt->tolerance >= 0.0) ||
      !(opt->epsilon >= 0.0) || opt->bias_count > MYB_FIT_MAX_BIASES ||
      (opt->bias_count > 0 && opt->bias == NULL)) {
    return MYB_E_DOMAIN;
  }
  for (size_t k = 0; k < opt->bias_count; ++k) {
    if (!isfinite(opt->bias[k])) {
      return MYB_E_DOMAIN;
    }
  }
  /* the observations: valid, ordered away from t0 */
  const double direction = observations[count - 1].t >= t0 ? 1.0 : -1.0;
  double previous = t0;
  for (size_t k = 0; k < count; ++k) {
    const myb_status valid = valid_observation(&observations[k], opt->bias_count);
    if (valid != MYB_OK) {
      return valid;
    }
    if ((observations[k].t - previous) * direction < 0.0) {
      return MYB_E_DOMAIN;
    }
    previous = observations[k].t;
  }
  fit f = {.model = model, .observations = observations, .count = count, .t0 = t0};
  /* an arc of no length still needs an interval to start on */
  f.t_end = observations[count - 1].t != t0 ? observations[count - 1].t : t0 + 1.0;
  f.np = myb_force_model_parameter_count(model);
  f.nb = opt->bias_count;
  f.bias = opt->bias;
  f.prior_sigma = opt->prior_sigma;
  const size_t nu = 6 + f.np + f.nb;
  for (size_t j = 0; j < nu; ++j) {
    f.prior[j] = j < 6 ? y0[j]
                 : j < 6 + f.np ? *myb_force_model_parameter(model, j - 6)
                                : f.bias[j - 6 - f.np];
    if (opt->prior_sigma != NULL && !(opt->prior_sigma[j] >= 0.0)) {
      return MYB_E_DOMAIN;
    }
    if (opt->fixed == NULL || !opt->fixed[j]) {
      f.free_index[f.nfree++] = j;
    }
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  const size_t nf = f.nfree;
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 2 * nf * nf + 3 * nf + 1, sizeof(double), &block);
  if (status != MYB_OK) {
    return status;
  }
  f.normal = block;
  f.scaled = f.normal + nf * nf;
  f.rhs = f.scaled + nf * nf;
  f.scale = f.rhs + nf;
  f.column = f.scale + nf;
  const myb_orbit_options orbit_options = {.epsilon = opt->epsilon, .variational = 1};
  status = myb_orbit_create(model, &orbit_options, allocator, &f.orbit);
  if (status != MYB_OK) {
    goto out;
  }
  const int max_iterations = opt->max_iterations > 0 ? opt->max_iterations : DEFAULT_ITERATIONS;
  const double tolerance = opt->tolerance > 0.0 ? opt->tolerance : default_tolerance;
  double rms_last = 0.0;
  for (int iteration = 0; iteration <= max_iterations; ++iteration) {
    double chi2 = 0.0;
    status = accumulate(&f, y0, residuals, &chi2, &result->values);
    if (status != MYB_OK) {
      goto out;
    }
    const double rms = sqrt(chi2 / (double)result->values);
    result->rms = rms;
    result->iterations = iteration;
    if (iteration == 0) {
      result->rms_before = rms;
    } else if (fabs(rms - rms_last) <= tolerance * rms_last) {
      result->converged = 1;
    }
    if (nf == 0) {
      result->converged = 1;
    }
    status = nf > 0 ? factorise(&f) : MYB_OK;
    if (status != MYB_OK) {
      goto out;
    }
    if (result->converged || iteration == max_iterations) {
      break;
    }
    /* the correction of the unknowns */
    for (size_t p = 0; p < nf; ++p) {
      f.column[p] = f.rhs[p] / f.scale[p];
    }
    myb_cholesky_solve(nf, f.scaled, f.column);
    for (size_t p = 0; p < nf; ++p) {
      const size_t index = f.free_index[p];
      const double step = f.column[p] / f.scale[p];
      if (index < 6) {
        y0[index] += step;
      } else if (index < 6 + f.np) {
        *myb_force_model_parameter(model, index - 6) += step;
      } else {
        f.bias[index - 6 - f.np] += step;
      }
    }
    rms_last = rms;
  }
  if (covariance != NULL) {
    memset(covariance, 0, nu * nu * sizeof(double));
    for (size_t p = 0; p < nf; ++p) {
      memset(f.column, 0, nf * sizeof(double));
      f.column[p] = 1.0;
      myb_cholesky_solve(nf, f.scaled, f.column);
      for (size_t q = 0; q < nf; ++q) {
        covariance[f.free_index[q] * nu + f.free_index[p]] =
            f.column[q] / (f.scale[p] * f.scale[q]);
      }
    }
  }
out:
  myb_orbit_destroy(f.orbit);
  myb_internal_release(&a, block);
  return status;
}
