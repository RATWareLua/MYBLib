/* SPDX-License-Identifier: BSD-3-Clause */
/* The environment of a force model: bodies with their providers, the time scales of the
 * providers, and the values at the current instant kept until the instant changes. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "env/env_internal.h"
#include "myb/env.h"
#include "myb/num.h"

static const double default_step = 3600.0;
static const double max_step = 1e15;              /* s: node times stay within myb_time */
static const double max_time = 1e17;              /* s: 3e9 years either way */
static const double max_nodes = 4503599627370496.0; /* 2^52: node indices are exact */

/* A scale other than the model's: (scale - model) at the model instants
 * epoch + (base + i) step, i = 0 .. 3. */
typedef struct link {
  myb_scale scale;
  int filled;
  int64_t base;
  double node[4];
} link;

typedef struct slot {
  myb_env_body body;
  int state_link; /* index of the provider's scale in links, -1: the model's */
  int orientation_link;
  int have_state;
  double state[6];
  int have_orientation;
  myb_orientation orientation;
  myb_rotation inverse;
} slot;

struct myb_env {
  myb_allocator allocator;
  myb_time epoch;
  const myb_timesys* timesys;
  myb_scale scale;
  double step;
  size_t link_count;
  link links[MYB_ENV_MAX_SCALES];
  size_t count;
  int origin;
  /* the instant the values below belong to */
  int cached;
  double cached_t;
  myb_time model_time;
  int have_time[MYB_ENV_MAX_SCALES];
  myb_time time[MYB_ENV_MAX_SCALES];
  slot slots[MYB_ENV_MAX_BODIES];
};

static int finite_all(const double* values, size_t count) {
  for (size_t i = 0; i < count; ++i) {
    if (!isfinite(values[i])) {
      return 0;
    }
  }
  return 1;
}

myb_status myb_env_create(const myb_env_options* options, const myb_allocator* allocator,
                          myb_env** env) {
  MYB_CHECK(options != NULL && env != NULL);
  double step = options->table_step;
  if (step == 0.0) {
    step = default_step;
  }
  if (!(step > 0.0) || !(step <= max_step)) {
    return MYB_E_DOMAIN;
  }
  myb_scale scale = {0};
  if (options->timesys != NULL) {
    if (options->scale == NULL) {
      return MYB_E_DOMAIN;
    }
    const myb_status status = myb_timesys_scale(options->timesys, options->scale, &scale);
    if (status != MYB_OK) {
      return status;
    }
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_env), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_env* e = block;
  memset(e, 0, sizeof *e);
  e->allocator = a;
  e->epoch = options->epoch;
  e->timesys = options->timesys;
  e->scale = scale;
  e->step = step;
  e->origin = -1;
  *env = e;
  return MYB_OK;
}

myb_status myb_env_clone(const myb_env* env, const myb_allocator* allocator, myb_env** clone) {
  MYB_CHECK(env != NULL && clone != NULL);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_env), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_env* e = block;
  memcpy(e, env, sizeof *e);
  e->allocator = a;
  /* nothing computed yet */
  e->cached = 0;
  for (size_t k = 0; k < MYB_ENV_MAX_SCALES; ++k) {
    e->links[k].filled = 0;
    e->have_time[k] = 0;
  }
  for (size_t k = 0; k < MYB_ENV_MAX_BODIES; ++k) {
    e->slots[k].have_state = 0;
    e->slots[k].have_orientation = 0;
  }
  *clone = e;
  return MYB_OK;
}

void myb_env_destroy(myb_env* env) {
  if (env == NULL) {
    return;
  }
  const myb_allocator a = env->allocator;
  myb_internal_release(&a, env);
}

/* The link of the scale `name` (NULL or the model's own: -1), added if new. */
static myb_status resolve_link(myb_env* env, const char* name, int* index) {
  *index = -1;
  if (name == NULL) {
    return MYB_OK;
  }
  if (env->timesys == NULL) {
    return MYB_E_DOMAIN;
  }
  myb_scale scale;
  const myb_status status = myb_timesys_scale(env->timesys, name, &scale);
  if (status != MYB_OK) {
    return status;
  }
  if (scale.index == env->scale.index) {
    return MYB_OK;
  }
  for (size_t k = 0; k < env->link_count; ++k) {
    if (env->links[k].scale.index == scale.index) {
      *index = (int)k;
      return MYB_OK;
    }
  }
  if (env->link_count == MYB_ENV_MAX_SCALES) {
    return MYB_E_OVERFLOW;
  }
  link* l = &env->links[env->link_count];
  memset(l, 0, sizeof *l);
  l->scale = scale;
  *index = (int)env->link_count;
  env->link_count++;
  return MYB_OK;
}

myb_status myb_env_add_body(myb_env* env, const myb_env_body* body, size_t* index) {
  MYB_CHECK(env != NULL && body != NULL);
  if (!isfinite(body->gm) || !(body->gm >= 0.0) || !isfinite(body->radius) ||
      !(body->radius >= 0.0) || !isfinite(body->luminosity) || !(body->luminosity >= 0.0) ||
      !finite_all(body->limb, 2) || !finite_all(body->spin, 3)) {
    return MYB_E_DOMAIN;
  }
  if (body->state == NULL && env->origin >= 0) {
    return MYB_E_DOMAIN;
  }
  if (env->count == MYB_ENV_MAX_BODIES) {
    return MYB_E_OVERFLOW;
  }
  const size_t links_before = env->link_count;
  int state_link = -1;
  int orientation_link = -1;
  myb_status status = body->state != NULL ? resolve_link(env, body->state_scale, &state_link)
                                          : MYB_OK;
  if (status == MYB_OK && body->orientation != NULL) {
    status = resolve_link(env, body->orientation_scale, &orientation_link);
  }
  if (status != MYB_OK) {
    env->link_count = links_before;
    return status;
  }
  slot* s = &env->slots[env->count];
  memset(s, 0, sizeof *s);
  s->body = *body;
  s->state_link = state_link;
  s->orientation_link = orientation_link;
  if (body->state == NULL) {
    env->origin = (int)env->count;
  }
  if (index != NULL) {
    *index = env->count;
  }
  env->count++;
  env->cached = 0;
  return MYB_OK;
}

size_t myb_env_body_count(const myb_env* env) {
  MYB_CHECK(env != NULL);
  return env->count;
}

double myb_env_gm(const myb_env* env, size_t body) {
  MYB_CHECK(env != NULL && body < env->count);
  return env->slots[body].body.gm;
}

myb_time myb_env_epoch(const myb_env* env) {
  MYB_CHECK(env != NULL);
  return env->epoch;
}

myb_status myb_env_find(const myb_env* env, int64_t id, size_t* index) {
  MYB_CHECK(env != NULL && index != NULL);
  for (size_t k = 0; k < env->count; ++k) {
    if (env->slots[k].body.id == id) {
      *index = k;
      return MYB_OK;
    }
  }
  return MYB_E_NOT_FOUND;
}

const myb_env_body* myb_internal_env_body(const myb_env* env, size_t body) {
  MYB_CHECK(env != NULL && body < env->count);
  return &env->slots[body].body;
}

int myb_internal_env_origin(const myb_env* env) {
  MYB_CHECK(env != NULL);
  return env->origin;
}

/* ---- Time ------------------------------------------------------------------------------------ */

/* Makes t the instant of the kept values, dropping those of another instant. */
static myb_status retarget(myb_env* env, double t) {
  if (!isfinite(t)) {
    return MYB_E_DOMAIN;
  }
  if (fabs(t) > max_time) {
    return MYB_E_RANGE;
  }
  if (env->cached && env->cached_t == t) {
    return MYB_OK;
  }
  env->cached = 1;
  env->cached_t = t;
  env->model_time = myb_time_add(env->epoch, t);
  for (size_t k = 0; k < env->link_count; ++k) {
    env->have_time[k] = 0;
  }
  for (size_t k = 0; k < env->count; ++k) {
    env->slots[k].have_state = 0;
    env->slots[k].have_orientation = 0;
  }
  return MYB_OK;
}

static myb_status node_value(const myb_env* env, const link* l, int64_t k, double* value) {
  const myb_time model = myb_time_add(env->epoch, (double)k * env->step);
  myb_time other;
  const myb_status status = myb_timesys_convert(env->timesys, model, env->scale, l->scale, &other);
  if (status != MYB_OK) {
    return status;
  }
  *value = myb_dd_to_double(myb_time_diff(other, model));
  return MYB_OK;
}

/* Makes the nodes k0 .. k0 + 3 those of the link; a node is a function of its index only, so
 * the ones already held are kept. */
static myb_status fill(const myb_env* env, link* l, int64_t k0) {
  if (l->filled && l->base == k0) {
    return MYB_OK;
  }
  double node[4];
  for (int i = 0; i < 4; ++i) {
    const int64_t k = k0 + i;
    if (l->filled && k >= l->base && k < l->base + 4) {
      node[i] = l->node[k - l->base];
      continue;
    }
    const myb_status status = node_value(env, l, k, &node[i]);
    if (status != MYB_OK) {
      return status;
    }
  }
  memcpy(l->node, node, sizeof node);
  l->base = k0;
  l->filled = 1;
  return MYB_OK;
}

/* (scale - model) at the model time t: Lagrange interpolation on the four nodes around t,
 * written about the node before t so that a large constant part of the difference does not
 * cancel. */
static myb_status link_offset(const myb_env* env, link* l, double t, double* offset) {
  const double u = t / env->step;
  if (!(fabs(u) < max_nodes)) {
    return MYB_E_RANGE;
  }
  const double whole = floor(u);
  const int64_t k0 = (int64_t)whole - 1;
  const double w = u - (whole - 1.0); /* in [1, 2): the abscissa counted from node k0 */
  const myb_status status = fill(env, l, k0);
  if (status != MYB_OK) {
    return status;
  }
  const double* f = l->node;
  const double l0 = -(w - 1.0) * (w - 2.0) * (w - 3.0) / 6.0;
  const double l2 = -w * (w - 1.0) * (w - 3.0) / 2.0;
  const double l3 = w * (w - 1.0) * (w - 2.0) / 6.0;
  *offset = f[1] + (l0 * (f[0] - f[1]) + l2 * (f[2] - f[1]) + l3 * (f[3] - f[1]));
  return MYB_OK;
}

/* The instant of the kept values on a link's scale. */
static myb_status link_time(myb_env* env, int index, myb_time* out) {
  if (index < 0) {
    *out = env->model_time;
    return MYB_OK;
  }
  if (!env->have_time[index]) {
    double offset = 0.0;
    const myb_status status = link_offset(env, &env->links[index], env->cached_t, &offset);
    if (status != MYB_OK) {
      return status;
    }
    env->time[index] = myb_time_add(env->model_time, offset);
    env->have_time[index] = 1;
  }
  *out = env->time[index];
  return MYB_OK;
}

myb_status myb_env_time(myb_env* env, double t, const char* scale, myb_time* out) {
  MYB_CHECK(env != NULL && out != NULL);
  int index = -1;
  if (scale != NULL) {
    if (env->timesys == NULL) {
      return MYB_E_NOT_FOUND;
    }
    myb_scale wanted;
    const myb_status found = myb_timesys_scale(env->timesys, scale, &wanted);
    if (found != MYB_OK) {
      return found;
    }
    if (wanted.index != env->scale.index) {
      for (size_t k = 0; k < env->link_count && index < 0; ++k) {
        if (env->links[k].scale.index == wanted.index) {
          index = (int)k;
        }
      }
      if (index < 0) {
        return MYB_E_NOT_FOUND;
      }
    }
  }
  const myb_status status = retarget(env, t);
  if (status != MYB_OK) {
    return status;
  }
  return link_time(env, index, out);
}

/* ---- Bodies ---------------------------------------------------------------------------------- */

myb_status myb_internal_env_state(myb_env* env, double t, size_t body, const double** state) {
  MYB_CHECK(env != NULL && body < env->count && state != NULL);
  myb_status status = retarget(env, t);
  if (status != MYB_OK) {
    return status;
  }
  slot* s = &env->slots[body];
  if (!s->have_state) {
    if (s->body.state == NULL) {
      memset(s->state, 0, sizeof s->state);
    } else {
      myb_time when;
      status = link_time(env, s->state_link, &when);
      if (status == MYB_OK) {
        status = s->body.state(s->body.state_context, when, s->state);
      }
      if (status != MYB_OK) {
        return status;
      }
      if (!finite_all(s->state, 6)) {
        return MYB_E_DOMAIN;
      }
    }
    s->have_state = 1;
  }
  *state = s->state;
  return MYB_OK;
}

myb_status myb_internal_env_orientation_time(myb_env* env, double t, size_t body, myb_time* out) {
  MYB_CHECK(env != NULL && body < env->count && out != NULL);
  const myb_status status = retarget(env, t);
  if (status != MYB_OK) {
    return status;
  }
  return link_time(env, env->slots[body].orientation_link, out);
}

myb_status myb_internal_env_orientation(myb_env* env, double t, size_t body,
                                        const myb_orientation** orientation,
                                        const myb_rotation** inverse) {
  MYB_CHECK(env != NULL && body < env->count && orientation != NULL);
  myb_status status = retarget(env, t);
  if (status != MYB_OK) {
    return status;
  }
  slot* s = &env->slots[body];
  if (!s->have_orientation) {
    myb_orientation* o = &s->orientation;
    if (s->body.orientation == NULL) {
      memset(o, 0, sizeof *o);
      o->rotation.m[0] = 1.0;
      o->rotation.m[4] = 1.0;
      o->rotation.m[8] = 1.0;
      o->ut = env->model_time;
    } else {
      myb_time when;
      status = link_time(env, s->orientation_link, &when);
      if (status == MYB_OK) {
        status = s->body.orientation(s->body.orientation_context, when, o);
      }
      if (status != MYB_OK) {
        return status;
      }
      if (!finite_all(o->rotation.m, 9) || !finite_all(o->rotation.dm, 9) || !isfinite(o->xp) ||
          !isfinite(o->yp)) {
        return MYB_E_DOMAIN;
      }
    }
    myb_rotation_invert(&o->rotation, &s->inverse);
    s->have_orientation = 1;
  }
  *orientation = &s->orientation;
  if (inverse != NULL) {
    *inverse = &s->inverse;
  }
  return MYB_OK;
}

myb_status myb_env_state(myb_env* env, double t, size_t body, double state[6]) {
  MYB_CHECK(state != NULL);
  const double* s = NULL;
  const myb_status status = myb_internal_env_state(env, t, body, &s);
  if (status == MYB_OK) {
    memcpy(state, s, 6 * sizeof(double));
  }
  return status;
}

myb_status myb_env_orientation(myb_env* env, double t, size_t body, myb_orientation* orientation) {
  MYB_CHECK(orientation != NULL);
  const myb_orientation* o = NULL;
  const myb_status status = myb_internal_env_orientation(env, t, body, &o, NULL);
  if (status == MYB_OK) {
    *orientation = *o;
  }
  return status;
}
