/* SPDX-License-Identifier: BSD-3-Clause */
/* The force model of one object: named terms, compensated sum, Jacobians. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/model.h"
#include "myb/num.h"

static const double cbrt_eps = 6.0554544523933395e-6; /* (2^-52)^(1/3) */

typedef struct term {
  char name[MYB_MODEL_NAME];
  myb_term_function function;
  myb_term_jacobian jacobian;
  void* context;
  int enabled;
} term;

typedef struct parameter {
  char name[MYB_MODEL_NAME];
  double* value;
  double step;
} parameter;

struct myb_force_model {
  myb_allocator allocator;
  term terms[MYB_MODEL_MAX_TERMS];
  size_t count;
  parameter parameters[MYB_MODEL_MAX_PARAMETERS];
  size_t parameter_count;
  double length_scale;
  double velocity_scale;
};

myb_status myb_force_model_create(const myb_allocator* allocator, myb_force_model** model) {
  MYB_CHECK(model != NULL);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_force_model), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_force_model* m = block;
  memset(m, 0, sizeof *m);
  m->allocator = a;
  m->length_scale = 1.0;
  m->velocity_scale = 1.0;
  *model = m;
  return MYB_OK;
}

void myb_force_model_destroy(myb_force_model* model) {
  if (model == NULL) {
    return;
  }
  const myb_allocator a = model->allocator;
  myb_internal_release(&a, model);
}

static int valid_name(const char* name) {
  const size_t n = strlen(name);
  return n > 0 && n < MYB_MODEL_NAME;
}

myb_status myb_force_model_find(const myb_force_model* model, const char* name, size_t* index) {
  MYB_CHECK(model != NULL && name != NULL && index != NULL);
  for (size_t k = 0; k < model->count; ++k) {
    if (strcmp(model->terms[k].name, name) == 0) {
      *index = k;
      return MYB_OK;
    }
  }
  return MYB_E_NOT_FOUND;
}

myb_status myb_force_model_add(myb_force_model* model, const char* name, myb_term_function function,
                               myb_term_jacobian jacobian, void* context, size_t* index) {
  MYB_CHECK(model != NULL && name != NULL && function != NULL);
  size_t existing = 0;
  if (!valid_name(name) || myb_force_model_find(model, name, &existing) == MYB_OK) {
    return MYB_E_DOMAIN;
  }
  if (model->count == MYB_MODEL_MAX_TERMS) {
    return MYB_E_OVERFLOW;
  }
  term* t = &model->terms[model->count];
  memset(t, 0, sizeof *t);
  memcpy(t->name, name, strlen(name) + 1);
  t->function = function;
  t->jacobian = jacobian;
  t->context = context;
  t->enabled = 1;
  if (index != NULL) {
    *index = model->count;
  }
  model->count++;
  return MYB_OK;
}

size_t myb_force_model_count(const myb_force_model* model) {
  MYB_CHECK(model != NULL);
  return model->count;
}

const char* myb_force_model_name(const myb_force_model* model, size_t index) {
  MYB_CHECK(model != NULL && index < model->count);
  return model->terms[index].name;
}

int myb_force_model_enabled(const myb_force_model* model, size_t index) {
  MYB_CHECK(model != NULL && index < model->count);
  return model->terms[index].enabled;
}

void myb_force_model_set_enabled(myb_force_model* model, size_t index, int enabled) {
  MYB_CHECK(model != NULL && index < model->count);
  model->terms[index].enabled = enabled != 0;
}

myb_status myb_force_model_add_parameter(myb_force_model* model, const char* name, double* value,
                                         double step, size_t* index) {
  MYB_CHECK(model != NULL && name != NULL && value != NULL);
  if (!valid_name(name) || !(step > 0.0) || !isfinite(step)) {
    return MYB_E_DOMAIN;
  }
  for (size_t k = 0; k < model->parameter_count; ++k) {
    if (strcmp(model->parameters[k].name, name) == 0) {
      return MYB_E_DOMAIN;
    }
  }
  if (model->parameter_count == MYB_MODEL_MAX_PARAMETERS) {
    return MYB_E_OVERFLOW;
  }
  parameter* p = &model->parameters[model->parameter_count];
  memcpy(p->name, name, strlen(name) + 1);
  p->value = value;
  p->step = step;
  if (index != NULL) {
    *index = model->parameter_count;
  }
  model->parameter_count++;
  return MYB_OK;
}

size_t myb_force_model_parameter_count(const myb_force_model* model) {
  MYB_CHECK(model != NULL);
  return model->parameter_count;
}

const char* myb_force_model_parameter_name(const myb_force_model* model, size_t index) {
  MYB_CHECK(model != NULL && index < model->parameter_count);
  return model->parameters[index].name;
}

double* myb_force_model_parameter(const myb_force_model* model, size_t index) {
  MYB_CHECK(model != NULL && index < model->parameter_count);
  return model->parameters[index].value;
}

myb_status myb_force_model_set_scales(myb_force_model* model, double length_scale,
                                      double velocity_scale) {
  MYB_CHECK(model != NULL);
  if (!(length_scale > 0.0) || !(velocity_scale > 0.0) || !isfinite(length_scale) ||
      !isfinite(velocity_scale)) {
    return MYB_E_DOMAIN;
  }
  model->length_scale = length_scale;
  model->velocity_scale = velocity_scale;
  return MYB_OK;
}

myb_status myb_force_model_term(myb_force_model* model, size_t index, double t, const double x[3],
                                const double v[3], double a[3]) {
  MYB_CHECK(model != NULL && index < model->count && x != NULL && v != NULL && a != NULL);
  const term* tm = &model->terms[index];
  return tm->function(tm->context, t, x, v, a);
}

myb_status myb_force_model_acceleration(myb_force_model* model, double t, const double x[3],
                                        const double v[3], double a[3]) {
  MYB_CHECK(model != NULL && x != NULL && v != NULL && a != NULL);
  myb_sum_acc sum[3] = {{0}};
  for (size_t k = 0; k < model->count; ++k) {
    const term* tm = &model->terms[k];
    if (!tm->enabled) {
      continue;
    }
    double ak[3];
    const myb_status st = tm->function(tm->context, t, x, v, ak);
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 3; ++i) {
      myb_sum_acc_add(&sum[i], ak[i]);
    }
  }
  for (int i = 0; i < 3; ++i) {
    a[i] = myb_sum_acc_value(&sum[i]);
  }
  return MYB_OK;
}

myb_status myb_force_model_accel(void* model, double t, const double* x, const double* v,
                                 double* a) {
  return myb_force_model_acceleration(model, t, x, v, a);
}

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/* Central differences of one term with respect to x and v. */
static myb_status numeric_term(const term* tm, double t, const double x[3], const double v[3],
                               double hx, double hv, double dadx[9], double dadv[9]) {
  for (int j = 0; j < 3; ++j) {
    double xp[3] = {x[0], x[1], x[2]};
    double xm[3] = {x[0], x[1], x[2]};
    xp[j] += hx;
    xm[j] -= hx;
    const double h2 = xp[j] - xm[j]; /* the step actually taken */
    double ap[3];
    double am[3];
    myb_status st = tm->function(tm->context, t, xp, v, ap);
    if (st == MYB_OK) {
      st = tm->function(tm->context, t, xm, v, am);
    }
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 3; ++i) {
      dadx[3 * i + j] = (ap[i] - am[i]) / h2;
    }
    double vp[3] = {v[0], v[1], v[2]};
    double vm[3] = {v[0], v[1], v[2]};
    vp[j] += hv;
    vm[j] -= hv;
    const double k2 = vp[j] - vm[j];
    st = tm->function(tm->context, t, x, vp, ap);
    if (st == MYB_OK) {
      st = tm->function(tm->context, t, x, vm, am);
    }
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 3; ++i) {
      dadv[3 * i + j] = (ap[i] - am[i]) / k2;
    }
  }
  return MYB_OK;
}

myb_status myb_force_model_jacobian(myb_force_model* model, double t, const double x[3],
                                    const double v[3], double dadx[9], double dadv[9],
                                    double* dadp) {
  MYB_CHECK(model != NULL && x != NULL && v != NULL && dadx != NULL && dadv != NULL);
  const double hx = cbrt_eps * fmax(norm3(x), model->length_scale);
  const double hv = cbrt_eps * fmax(norm3(v), model->velocity_scale);
  myb_sum_acc sx[9] = {{0}};
  myb_sum_acc sv[9] = {{0}};
  for (size_t k = 0; k < model->count; ++k) {
    const term* tm = &model->terms[k];
    if (!tm->enabled) {
      continue;
    }
    double jx[9];
    double jv[9];
    const myb_status st = tm->jacobian != NULL
                              ? tm->jacobian(tm->context, t, x, v, jx, jv)
                              : numeric_term(tm, t, x, v, hx, hv, jx, jv);
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 9; ++i) {
      myb_sum_acc_add(&sx[i], jx[i]);
      myb_sum_acc_add(&sv[i], jv[i]);
    }
  }
  for (int i = 0; i < 9; ++i) {
    dadx[i] = myb_sum_acc_value(&sx[i]);
    dadv[i] = myb_sum_acc_value(&sv[i]);
  }
  if (dadp != NULL) {
    const size_t np = model->parameter_count;
    for (size_t k = 0; k < np; ++k) {
      parameter* p = &model->parameters[k];
      const double saved = *p->value;
      double ap[3];
      double am[3];
      *p->value = saved + p->step;
      const double up = *p->value;
      myb_status st = myb_force_model_acceleration(model, t, x, v, ap);
      *p->value = saved - p->step;
      const double down = *p->value;
      if (st == MYB_OK) {
        st = myb_force_model_acceleration(model, t, x, v, am);
      }
      *p->value = saved;
      if (st != MYB_OK) {
        return st;
      }
      for (int i = 0; i < 3; ++i) {
        dadp[(size_t)i * np + k] = (ap[i] - am[i]) / (up - down);
      }
    }
  }
  return MYB_OK;
}

void myb_point_mass_jacobian(double gm, const double r[3], double dadx[9]) {
  MYB_CHECK(r != NULL && dadx != NULL);
  const double rr = norm3(r);
  const double r3 = rr * rr * rr;
  const double r5 = r3 * rr * rr;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      dadx[3 * i + j] = gm * (3.0 * r[i] * r[j] / r5 - (i == j ? 1.0 / r3 : 0.0));
    }
  }
}
