/* SPDX-License-Identifier: BSD-3-Clause */
/* glue.c — the C half of the Go binding (glue.h): the records of a description copied into
 * the structures of the library, and the calls a C program would make, in its order. No
 * arithmetic on numbers: only copies, comparisons and the indices of arrays. */
#include "glue.h"

#include <math.h> /* isfinite: a classification, no arithmetic */
#include <stdlib.h>
#include <string.h>

#include "glue_internal.h"

/* ---- The model ---------------------------------------------------------------------------------
 */

typedef struct body_store {
  mybgo_provider state;
  mybgo_orientation_store orientation;
} body_store;

typedef struct term_store {
  int kind;
  size_t index; /* of the term in the model */
  myb_point_mass_term point_mass;
  myb_third_body_term third_body;
  myb_field_term field;
  myb_gravity_field* own_field; /* the copy of a field that follows time */
  myb_gravity_work* work;
  myb_drag_term drag;
  myb_radiation_term radiation;
  myb_relativity_term relativity;
  myb_empirical_term empirical;
  myb_object object;
  myb_exponential_atmosphere atmosphere;
  myb_eclipse_event breaks[2 * MYB_ENV_MAX_OCCULTERS];
  mybgo_go_context go;
} term_store;

struct mybgo_model {
  myb_env* env;
  myb_force_model* model;
  size_t body_count;
  body_store* bodies;
  size_t term_count;
  term_store* terms;
};

static int body_index(const mybgo_model* m, size_t body) {
  return body < m->body_count;
}

static myb_status add_body(mybgo_model* m, const mybgo_body* d, body_store* s) {
  myb_env_body b = {0};
  b.id = d->id;
  b.gm = d->gm;
  b.radius = d->radius;
  b.luminosity = d->luminosity;
  b.limb[0] = d->limb[0];
  b.limb[1] = d->limb[1];
  b.spin[0] = d->spin[0];
  b.spin[1] = d->spin[1];
  b.spin[2] = d->spin[2];
  if (d->state.kind != MYBGO_PROVIDER_NONE) {
    s->state = d->state;
    b.state = mybgo_provider_state;
    b.state_context = &s->state;
    b.state_scale = d->state_scale;
  }
  myb_orientation_function orientation = NULL;
  void* context = NULL;
  const myb_status status =
      mybgo_orientation_of(&d->orientation, &s->orientation, &orientation, &context);
  if (status != MYB_OK) {
    return status;
  }
  if (orientation != NULL) {
    b.orientation = orientation;
    b.orientation_context = context;
    b.orientation_scale = d->orientation_scale;
  }
  size_t index = 0;
  return myb_env_add_body(m->env, &b, &index);
}

static myb_status add_field(mybgo_model* m, const mybgo_term* d, term_store* s) {
  if (d->field == NULL || !body_index(m, d->body)) {
    return MYB_E_DOMAIN;
  }
  myb_gravity_field* field = d->field;
  if (d->update != 0) {
    const myb_status status = myb_gravity_field_clone(d->field, NULL, &s->own_field);
    if (status != MYB_OK) {
      return status;
    }
    field = s->own_field;
  }
  const myb_status status =
      myb_gravity_work_create(myb_gravity_field_degree(field), NULL, &s->work);
  if (status != MYB_OK) {
    return status;
  }
  s->field = (myb_field_term){.env = m->env,
                              .body = d->body,
                              .field = field,
                              .work = s->work,
                              .n_min = d->n_min,
                              .n_max = d->n_max,
                              .update = d->update};
  return myb_force_model_add(m->model, d->name, myb_term_field, myb_term_field_jacobian, &s->field,
                             &s->index);
}

static myb_status add_relativity(mybgo_model* m, const mybgo_term* d, term_store* s) {
  if (!body_index(m, d->body) || (d->kind == MYBGO_TERM_DE_SITTER && !body_index(m, d->external))) {
    return MYB_E_DOMAIN;
  }
  s->relativity = (myb_relativity_term){
      .env = m->env, .ppn = d->ppn, .central = d->body, .external = d->external};
  myb_term_function function = myb_term_schwarzschild;
  myb_term_jacobian jacobian = myb_term_schwarzschild_jacobian;
  if (d->kind == MYBGO_TERM_LENSE_THIRRING) {
    function = myb_term_lense_thirring;
    jacobian = myb_term_lense_thirring_jacobian;
  } else if (d->kind == MYBGO_TERM_DE_SITTER) {
    function = myb_term_de_sitter;
    jacobian = myb_term_de_sitter_jacobian;
  }
  return myb_force_model_add(m->model, d->name, function, jacobian, &s->relativity, &s->index);
}

static myb_status add_radiation(mybgo_model* m, const mybgo_term* d, term_store* s) {
  if (!body_index(m, d->body) || d->occulter_count > MYB_ENV_MAX_OCCULTERS) {
    return MYB_E_DOMAIN;
  }
  s->object = d->object;
  s->object.mass_function = NULL;
  s->object.mass_context = NULL;
  s->radiation = (myb_radiation_term){.env = m->env,
                                      .source = d->body,
                                      .object = &s->object,
                                      .c = d->c,
                                      .occulter_count = d->occulter_count,
                                      .moving = d->moving,
                                      .atmospheres = NULL};
  for (size_t k = 0; k < d->occulter_count; ++k) {
    if (!body_index(m, d->occulters[k])) {
      return MYB_E_DOMAIN;
    }
    s->radiation.occulters[k] = d->occulters[k];
  }
  myb_status status = myb_force_model_add(m->model, d->name, myb_term_radiation,
                                          myb_term_radiation_jacobian, &s->radiation, &s->index);
  if (status == MYB_OK && d->breaks != 0) {
    status = myb_term_radiation_breaks(&s->radiation, m->model, s->breaks);
  }
  return status;
}

static myb_status add_term(mybgo_model* m, const mybgo_term* d, term_store* s) {
  if (d->name == NULL) {
    return MYB_E_DOMAIN;
  }
  s->kind = d->kind;
  myb_status status = MYB_OK;
  switch (d->kind) {
    case MYBGO_TERM_POINT_MASS:
      s->point_mass.gm = d->gm;
      status = myb_force_model_add(m->model, d->name, myb_term_point_mass,
                                   myb_term_point_mass_jacobian, &s->point_mass, &s->index);
      break;
    case MYBGO_TERM_THIRD_BODY:
      if (!body_index(m, d->body)) {
        return MYB_E_DOMAIN;
      }
      s->third_body = (myb_third_body_term){.env = m->env, .body = d->body};
      status = myb_force_model_add(m->model, d->name, myb_term_third_body,
                                   myb_term_third_body_jacobian, &s->third_body, &s->index);
      break;
    case MYBGO_TERM_FIELD:
      status = add_field(m, d, s);
      break;
    case MYBGO_TERM_DRAG:
      if (!body_index(m, d->body)) {
        return MYB_E_DOMAIN;
      }
      s->object = d->object;
      s->object.mass_function = NULL;
      s->object.mass_context = NULL;
      s->atmosphere = d->atmosphere;
      s->drag = (myb_drag_term){.env = m->env,
                                .body = d->body,
                                .density = myb_density_exponential,
                                .density_context = &s->atmosphere,
                                .object = &s->object};
      status = myb_force_model_add(m->model, d->name, myb_term_drag, myb_term_drag_jacobian,
                                   &s->drag, &s->index);
      break;
    case MYBGO_TERM_RADIATION:
      status = add_radiation(m, d, s);
      break;
    case MYBGO_TERM_SCHWARZSCHILD:
    case MYBGO_TERM_LENSE_THIRRING:
    case MYBGO_TERM_DE_SITTER:
      status = add_relativity(m, d, s);
      break;
    case MYBGO_TERM_EMPIRICAL:
      memcpy(s->empirical.constant, d->constant, sizeof d->constant);
      memcpy(s->empirical.cosine, d->cosine, sizeof d->cosine);
      memcpy(s->empirical.sine, d->sine, sizeof d->sine);
      status = myb_force_model_add(m->model, d->name, myb_term_empirical,
                                   myb_term_empirical_jacobian, &s->empirical, &s->index);
      break;
    case MYBGO_TERM_GO:
      s->go.handle = d->handle;
      status = myb_force_model_add(m->model, d->name, mybgo_go_term, NULL, &s->go, &s->index);
      break;
    default:
      return MYB_E_DOMAIN;
  }
  if (status == MYB_OK && d->position_only != 0) {
    myb_force_model_set_position_only(m->model, s->index, 1);
  }
  return status;
}

static myb_status add_parameter(mybgo_model* m, const mybgo_parameter* d) {
  if (d->name == NULL || d->term >= m->term_count) {
    return MYB_E_DOMAIN;
  }
  term_store* s = &m->terms[d->term];
  double* value = NULL;
  const int in_range = d->index >= 0 && d->index < 3;
  if (d->what == MYBGO_PARAMETER_CD && s->kind == MYBGO_TERM_DRAG) {
    value = &s->object.cd;
  } else if (d->what == MYBGO_PARAMETER_CR && s->kind == MYBGO_TERM_RADIATION) {
    value = &s->object.cr;
  } else if (s->kind == MYBGO_TERM_EMPIRICAL && in_range) {
    if (d->what == MYBGO_PARAMETER_CONSTANT) {
      value = &s->empirical.constant[d->index];
    } else if (d->what == MYBGO_PARAMETER_COSINE) {
      value = &s->empirical.cosine[d->index];
    } else if (d->what == MYBGO_PARAMETER_SINE) {
      value = &s->empirical.sine[d->index];
    }
  }
  if (value == NULL) {
    return MYB_E_DOMAIN;
  }
  size_t index = 0;
  const myb_status status =
      myb_force_model_add_parameter(m->model, d->name, value, d->step, &index);
  if (status != MYB_OK) {
    return status;
  }
  return myb_force_model_bind_parameter(m->model, index, s->index);
}

void mybgo_model_destroy(mybgo_model* model) {
  if (model == NULL) {
    return;
  }
  myb_force_model_destroy(model->model);
  myb_env_destroy(model->env);
  for (size_t k = 0; k < model->term_count && model->terms != NULL; ++k) {
    myb_gravity_work_destroy(model->terms[k].work);
    myb_gravity_field_destroy(model->terms[k].own_field);
  }
  free(model->terms);
  free(model->bodies);
  free(model);
}

myb_status mybgo_model_create(const mybgo_model_options* options, const mybgo_body* bodies,
                              size_t body_count, const mybgo_term* terms, size_t term_count,
                              const mybgo_parameter* parameters, size_t parameter_count,
                              mybgo_model** model, int* stage, size_t* item) {
  *model = NULL;
  *stage = MYBGO_STAGE_ENVIRONMENT;
  *item = 0;
  if (options == NULL || (body_count > 0 && bodies == NULL) || (term_count > 0 && terms == NULL) ||
      (parameter_count > 0 && parameters == NULL)) {
    return MYB_E_DOMAIN;
  }
  if (body_count > MYB_ENV_MAX_BODIES || term_count > MYB_MODEL_MAX_TERMS ||
      parameter_count > MYB_MODEL_MAX_PARAMETERS) {
    return MYB_E_OVERFLOW;
  }
  mybgo_model* m = calloc(1, sizeof *m);
  if (m == NULL) {
    return MYB_E_NOMEM;
  }
  myb_status status = MYB_E_NOMEM;
  m->bodies = calloc(body_count > 0 ? body_count : 1, sizeof *m->bodies);
  m->terms = calloc(term_count > 0 ? term_count : 1, sizeof *m->terms);
  if (m->bodies == NULL || m->terms == NULL) {
    goto fail;
  }
  const myb_env_options env_options = {.epoch = options->epoch,
                                       .timesys = options->timesys,
                                       .scale = options->scale,
                                       .table_step = options->table_step};
  status = myb_env_create(&env_options, NULL, &m->env);
  if (status == MYB_OK) {
    status = myb_force_model_create(NULL, &m->model);
  }
  if (status != MYB_OK) {
    goto fail;
  }
  *stage = MYBGO_STAGE_BODY;
  for (size_t i = 0; i < body_count; ++i) {
    *item = i;
    status = add_body(m, &bodies[i], &m->bodies[i]);
    if (status != MYB_OK) {
      goto fail;
    }
    m->body_count = i + 1;
  }
  *stage = MYBGO_STAGE_TERM;
  for (size_t j = 0; j < term_count; ++j) {
    *item = j;
    m->term_count = j + 1; /* the work and the field of the term are released on failure too */
    status = add_term(m, &terms[j], &m->terms[j]);
    if (status != MYB_OK) {
      goto fail;
    }
  }
  *stage = MYBGO_STAGE_PARAMETER;
  for (size_t p = 0; p < parameter_count; ++p) {
    *item = p;
    status = add_parameter(m, &parameters[p]);
    if (status != MYB_OK) {
      goto fail;
    }
  }
  *stage = 0;
  *item = 0;
  *model = m;
  return MYB_OK;
fail:
  mybgo_model_destroy(m);
  return status;
}

myb_force_model* mybgo_model_force(mybgo_model* model) {
  return model->model;
}

myb_env* mybgo_model_env(mybgo_model* model) {
  return model->env;
}

size_t mybgo_model_bodies(const mybgo_model* model) {
  return model->body_count;
}

size_t mybgo_model_parameter_count(const mybgo_model* model) {
  return myb_force_model_parameter_count(model->model);
}

void mybgo_model_get_parameters(const mybgo_model* model, double* values) {
  const size_t count = myb_force_model_parameter_count(model->model);
  for (size_t k = 0; k < count; ++k) {
    values[k] = *myb_force_model_parameter(model->model, k);
  }
}

myb_status mybgo_model_set_parameters(mybgo_model* model, const double* values) {
  const size_t count = myb_force_model_parameter_count(model->model);
  for (size_t k = 0; k < count; ++k) {
    if (!isfinite(values[k])) {
      return MYB_E_DOMAIN;
    }
  }
  for (size_t k = 0; k < count; ++k) {
    *myb_force_model_parameter(model->model, k) = values[k];
  }
  return MYB_OK;
}

myb_status mybgo_model_accelerations(mybgo_model* model, const double* t, const double* y,
                                     size_t count, double* a, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_force_model_acceleration(model->model, t[k], y + 6 * k, y + 6 * k + 3, a + 3 * k);
    if (status != MYB_OK) {
      break;
    }
  }
  *failed = k;
  return status;
}

/* ---- Propagation --------------------------------------------------------------------------------
 */

/* An impulse at kick_t comes before the state at t when it is not later in the direction of
 * the propagation (a state at the time of an impulse is the one after it). */
static int kick_first(int forward, double kick_t, double t) {
  return forward ? kick_t <= t : kick_t >= t;
}

myb_status mybgo_propagate(mybgo_model* model, const myb_orbit_options* options, double t0,
                           const double* y0, double t_end, const double* times, size_t count,
                           double* states, double* phi, double* sensitivity, const double* kicks,
                           size_t kick_count, int64_t* steps, size_t* failed) {
  *steps = 0;
  *failed = 0;
  myb_orbit* orbit = NULL;
  myb_status status = myb_orbit_create(model->model, options, NULL, &orbit);
  if (status == MYB_OK) {
    status = myb_orbit_start(orbit, t0, y0, t_end);
  }
  const int forward = t_end > t0;
  const size_t np = myb_force_model_parameter_count(model->model);
  const int partials = phi != NULL || sensitivity != NULL;
  size_t k = 0;
  size_t j = 0;
  /* bound: every pass takes one time or one impulse */
  while (status == MYB_OK && (k < count || j < kick_count)) {
    if (j < kick_count && (k == count || kick_first(forward, kicks[4 * j], times[k]))) {
      double y[6];
      status = myb_orbit_state(orbit, kicks[4 * j], y);
      if (status == MYB_OK) {
        status = myb_orbit_kick(orbit, kicks + 4 * j + 1);
      }
      ++j;
    } else {
      status = partials ? myb_orbit_partials(orbit, times[k], states + 6 * k,
                                             phi != NULL ? phi + 36 * k : NULL,
                                             sensitivity != NULL ? sensitivity + 6 * np * k : NULL)
                        : myb_orbit_state(orbit, times[k], states + 6 * k);
      if (status == MYB_OK) {
        ++k;
      }
    }
  }
  if (orbit != NULL) {
    *steps = myb_orbit_steps(orbit);
  }
  *failed = k;
  myb_orbit_destroy(orbit);
  return status;
}

/* ---- Events ---------------------------------------------------------------------------------- */

typedef struct event_context {
  myb_plane_event plane;
  myb_range_event range;
  myb_altitude_event altitude;
  myb_station station;
  myb_elevation_event elevation;
  myb_eclipse_event eclipse;
  mybgo_go_context go;
} event_context;

/* The function and the context of an event description. */
static myb_status event_of(mybgo_model* m, const mybgo_event* e, event_context* c,
                           myb_event_function* g, void** context) {
  const int body_or_origin = e->body == MYB_EVENT_ORIGIN || body_index(m, e->body);
  if (e->direction < -1 || e->direction > 1 || e->limit == 0) {
    return MYB_E_DOMAIN;
  }
  switch (e->kind) {
    case MYBGO_EVENT_PLANE:
      if (!body_or_origin) {
        return MYB_E_DOMAIN;
      }
      c->plane = (myb_plane_event){.env = m->env, .body = e->body};
      memcpy(c->plane.normal, e->normal, sizeof e->normal);
      *g = myb_event_plane;
      *context = &c->plane;
      return MYB_OK;
    case MYBGO_EVENT_RANGE:
    case MYBGO_EVENT_RANGE_RATE:
      if (!body_or_origin) {
        return MYB_E_DOMAIN;
      }
      c->range = (myb_range_event){.env = m->env, .body = e->body, .range = e->range};
      *g = e->kind == MYBGO_EVENT_RANGE ? myb_event_range : myb_event_range_rate;
      *context = &c->range;
      return MYB_OK;
    case MYBGO_EVENT_ALTITUDE:
      if (!body_index(m, e->body)) {
        return MYB_E_DOMAIN;
      }
      c->altitude = (myb_altitude_event){
          .env = m->env, .body = e->body, .ellipsoid = e->ellipsoid, .altitude = e->altitude};
      *g = myb_event_altitude;
      *context = &c->altitude;
      return MYB_OK;
    case MYBGO_EVENT_ELEVATION: {
      if (!body_index(m, e->body)) {
        return MYB_E_DOMAIN;
      }
      const myb_status status = myb_station_init(&c->station, m->env, e->body, &e->ellipsoid,
                                                 e->longitude, e->latitude, e->height);
      if (status != MYB_OK) {
        return status;
      }
      c->elevation = (myb_elevation_event){.station = &c->station, .elevation = e->elevation};
      *g = myb_event_elevation;
      *context = &c->elevation;
      return MYB_OK;
    }
    case MYBGO_EVENT_ECLIPSE:
      if (!body_index(m, e->source) || !body_index(m, e->occulter)) {
        return MYB_E_DOMAIN;
      }
      c->eclipse = (myb_eclipse_event){
          .env = m->env, .source = e->source, .occulter = e->occulter, .inner = e->inner};
      *g = myb_event_eclipse;
      *context = &c->eclipse;
      return MYB_OK;
    case MYBGO_EVENT_GO:
      c->go.handle = e->handle;
      *g = mybgo_go_event;
      *context = &c->go;
      return MYB_OK;
    default:
      return MYB_E_DOMAIN;
  }
}

typedef struct found_list {
  mybgo_found* items;
  size_t count;
  size_t capacity;
} found_list;

static myb_status append(found_list* list, size_t event, double t, const double y[6]) {
  if (list->count == list->capacity) {
    const size_t capacity = list->capacity > 0 ? 2 * list->capacity : 64;
    if (capacity > SIZE_MAX / sizeof(mybgo_found)) {
      return MYB_E_OVERFLOW;
    }
    mybgo_found* items = realloc(list->items, capacity * sizeof(mybgo_found));
    if (items == NULL) {
      return MYB_E_NOMEM;
    }
    list->items = items;
    list->capacity = capacity;
  }
  mybgo_found* f = &list->items[list->count++];
  f->event = event;
  f->t = t;
  memcpy(f->y, y, sizeof f->y);
  return MYB_OK;
}

/* The roots of one event (all of them when g is NULL: the propagation alone), then the state at
 * t_end and the steps when asked. */
static myb_status search(mybgo_model* m, const myb_orbit_options* options, double t0,
                         const double* y0, double t_end, myb_event_function g, void* context,
                         const mybgo_event* e, size_t index, found_list* list, double* y_end,
                         int64_t* steps) {
  myb_orbit* orbit = NULL;
  myb_status status = myb_orbit_create(m->model, options, NULL, &orbit);
  if (status == MYB_OK) {
    status = myb_orbit_start(orbit, t0, y0, t_end);
  }
  /* bound: at most e->limit roots */
  for (size_t r = 0; status == MYB_OK && g != NULL && r < e->limit; ++r) {
    int found = 0;
    double t = 0.0;
    double y[6];
    status = myb_orbit_event(orbit, g, context, t_end, e->direction, &found, &t, y);
    if (status != MYB_OK || !found) {
      break;
    }
    status = append(list, index, t, y);
  }
  if (status == MYB_OK && y_end != NULL) {
    status = myb_orbit_state(orbit, t_end, y_end);
    *steps = myb_orbit_steps(orbit);
  }
  myb_orbit_destroy(orbit);
  return status;
}

myb_status mybgo_events(mybgo_model* model, const myb_orbit_options* options, double t0,
                        const double* y0, double t_end, const mybgo_event* events,
                        size_t event_count, mybgo_found** found, size_t* found_count, double* y_end,
                        int64_t* steps, size_t* failed) {
  *found = NULL;
  *found_count = 0;
  *steps = 0;
  *failed = 0;
  found_list list = {0};
  myb_status status = MYB_OK;
  event_context* context = calloc(1, sizeof *context);
  if (context == NULL) {
    return MYB_E_NOMEM;
  }
  if (event_count == 0) {
    const mybgo_event none = {0};
    status = search(model, options, t0, y0, t_end, NULL, NULL, &none, 0, &list, y_end, steps);
  }
  for (size_t k = 0; k < event_count && status == MYB_OK; ++k) {
    *failed = k;
    myb_event_function g = NULL;
    void* c = NULL;
    status = event_of(model, &events[k], context, &g, &c);
    if (status == MYB_OK) {
      status = search(model, options, t0, y0, t_end, g, c, &events[k], k, &list,
                      k == 0 ? y_end : NULL, steps);
    }
  }
  free(context);
  if (status != MYB_OK) {
    free(list.items);
    return status;
  }
  *failed = event_count;
  *found = list.items;
  *found_count = list.count;
  return MYB_OK;
}

void mybgo_free(void* pointer) {
  free(pointer);
}

/* ---- Estimation ---------------------------------------------------------------------------- */

myb_status mybgo_stations_of(mybgo_model* model, const mybgo_station* stations, size_t count,
                             const mybgo_light* light, myb_station* st, myb_light_measure* lm) {
  for (size_t s = 0; s < count; ++s) {
    const mybgo_station* d = &stations[s];
    if (!body_index(model, d->body)) {
      return MYB_E_DOMAIN;
    }
    const myb_status status = myb_station_init(&st[s], model->env, d->body, &d->ellipsoid,
                                               d->longitude, d->latitude, d->height);
    if (status != MYB_OK) {
      return status;
    }
    if (light != NULL) {
      lm[s] = (myb_light_measure){.station = &st[s],
                                  .c = light->c,
                                  .gm = light->gm,
                                  .two_way = light->two_way,
                                  .aberration = light->aberration};
    }
  }
  return MYB_OK;
}

myb_status mybgo_measure_of(const mybgo_observation* o, myb_station* stations,
                            myb_light_measure* light, size_t station_count, myb_observation* out) {
  const int with_light = (o->kind & MYBGO_MEASURE_LIGHT) != 0;
  const int kind = o->kind & ~MYBGO_MEASURE_LIGHT;
  out->t = o->t;
  out->bias = o->bias;
  memcpy(out->value, o->value, sizeof out->value);
  memcpy(out->sigma, o->sigma, sizeof out->sigma);
  if (kind == MYBGO_MEASURE_POSITION || kind == MYBGO_MEASURE_STATE) {
    if (with_light) {
      return MYB_E_DOMAIN;
    }
    out->measure = kind == MYBGO_MEASURE_POSITION ? myb_measure_position : myb_measure_state;
    out->context = NULL;
    out->count = kind == MYBGO_MEASURE_POSITION ? 3 : 6;
    return MYB_OK;
  }
  if (o->station >= station_count || (with_light && light == NULL)) {
    return MYB_E_DOMAIN;
  }
  out->context = with_light ? (void*)&light[o->station] : (void*)&stations[o->station];
  switch (kind) {
    case MYBGO_MEASURE_RANGE:
      out->measure = with_light ? myb_measure_light_range : myb_measure_range;
      out->count = 1;
      return MYB_OK;
    case MYBGO_MEASURE_RANGE_RATE:
      out->measure = with_light ? myb_measure_light_range_rate : myb_measure_range_rate;
      out->count = 1;
      return MYB_OK;
    case MYBGO_MEASURE_RADEC:
      out->measure = with_light ? myb_measure_light_radec : myb_measure_radec;
      out->count = 2;
      out->periodic = 1;
      return MYB_OK;
    case MYBGO_MEASURE_AZEL:
      out->measure = with_light ? myb_measure_light_azel : myb_measure_azel;
      out->count = 2;
      out->periodic = 1;
      return MYB_OK;
    default:
      return MYB_E_DOMAIN;
  }
}

myb_status mybgo_fit(mybgo_model* model, double t0, double* y0, const mybgo_station* stations,
                     size_t station_count, const mybgo_light* light,
                     const mybgo_observation* observations, size_t count,
                     const mybgo_fit_options* options, const unsigned char* fixed,
                     const double* prior_sigma, double* bias, myb_fit_result* result,
                     double* covariance, double* residuals, size_t* failed) {
  *failed = count;
  myb_status status = MYB_E_NOMEM;
  myb_station* st = calloc(station_count > 0 ? station_count : 1, sizeof *st);
  myb_light_measure* lm = calloc(station_count > 0 ? station_count : 1, sizeof *lm);
  myb_observation* obs = calloc(count > 0 ? count : 1, sizeof *obs);
  if (st == NULL || lm == NULL || obs == NULL) {
    goto cleanup;
  }
  status = mybgo_stations_of(model, stations, station_count, light, st, lm);
  if (status != MYB_OK) {
    goto cleanup;
  }
  for (size_t k = 0; k < count; ++k) {
    status =
        mybgo_measure_of(&observations[k], st, light != NULL ? lm : NULL, station_count, &obs[k]);
    if (status != MYB_OK) {
      *failed = k;
      goto cleanup;
    }
  }
  const myb_fit_options fo = {.max_iterations = options->max_iterations,
                              .tolerance = options->tolerance,
                              .epsilon = options->epsilon,
                              .fixed = fixed,
                              .prior_sigma = prior_sigma,
                              .bias_count = options->bias_count,
                              .bias = bias};
  status =
      myb_fit_orbit(model->model, t0, y0, obs, count, &fo, NULL, result, covariance, residuals);
cleanup:
  free(obs);
  free(lm);
  free(st);
  return status;
}

/* ---- Files ------------------------------------------------------------------------------------
 */

myb_status mybgo_oem_write(const char* originator, const myb_calendar* created, int notation,
                           size_t segment_count, const char* const* names, const int* degrees,
                           const size_t* counts, const myb_calendar* epochs, const double* states,
                           char** text, size_t* length) {
  *text = NULL;
  *length = 0;
  myb_oem_segment_data* segments = calloc(segment_count > 0 ? segment_count : 1, sizeof *segments);
  if (segments == NULL) {
    return MYB_E_NOMEM;
  }
  size_t first = 0;
  for (size_t j = 0; j < segment_count; ++j) {
    myb_oem_segment_data* s = &segments[j];
    s->object_name = names[6 * j];
    s->object_id = names[6 * j + 1];
    s->center_name = names[6 * j + 2];
    s->ref_frame = names[6 * j + 3];
    s->time_system = names[6 * j + 4];
    s->interpolation = names[6 * j + 5];
    s->interpolation_degree = degrees[j];
    s->count = counts[j];
    s->epochs = counts[j] > 0 ? epochs + first : NULL;
    s->states = counts[j] > 0 ? states + 6 * first : NULL;
    first += counts[j];
  }
  const myb_oem_data data = {.originator = originator,
                             .creation_date = *created,
                             .segment_count = segment_count,
                             .segments = segment_count > 0 ? segments : NULL};
  const myb_odm_notation as = notation == MYB_ODM_XML ? MYB_ODM_XML : MYB_ODM_KVN;
  size_t needed = 0;
  myb_status status = myb_oem_write_as(&data, as, NULL, 0, &needed);
  if (status == MYB_E_OVERFLOW) {
    char* buffer = malloc(needed + 1);
    status = buffer != NULL ? myb_oem_write_as(&data, as, buffer, needed + 1, length) : MYB_E_NOMEM;
    if (status == MYB_OK) {
      *text = buffer;
    } else {
      free(buffer);
    }
  }
  free(segments);
  return status;
}

void mybgo_oem_meta(const myb_oem* oem, size_t segment, const char* const* keys, size_t n,
                    const char** values) {
  for (size_t k = 0; k < n; ++k) {
    values[k] = myb_oem_meta(oem, segment, keys[k]);
  }
}

size_t mybgo_oem_lines(const myb_oem* oem, size_t segment, myb_calendar* epochs, double* states,
                       double* accelerations) {
  const size_t count = myb_oem_state_count(oem, segment);
  size_t with = 0;
  for (size_t k = 0; k < count; ++k) {
    const int flags = myb_oem_state(oem, segment, k, &epochs[k], states + 6 * k,
                                    accelerations != NULL ? accelerations + 3 * k : NULL);
    if (flags == 3) {
      ++with;
    }
  }
  return with;
}

myb_status mybgo_spk_write(const char* path, const char* internal_name, size_t segment_count,
                           const myb_spk_segment_header* headers, const int* types,
                           const int* degrees, const size_t* counts, const double* epochs,
                           const double* states, void** data, size_t* size, size_t* failed) {
  *failed = segment_count;
  if (path == NULL && (data == NULL || size == NULL)) {
    return MYB_E_DOMAIN;
  }
  if (data != NULL) {
    *data = NULL;
  }
  if (size != NULL) {
    *size = 0;
  }
  myb_spk_writer* writer = NULL;
  myb_status status = myb_spk_writer_create(internal_name, NULL, &writer);
  size_t first = 0;
  for (size_t j = 0; j < segment_count && status == MYB_OK; ++j) {
    status = myb_spk_writer_add_states(writer, &headers[j], types[j], degrees[j], epochs + first,
                                       states + 6 * first, counts[j]);
    if (status != MYB_OK) {
      *failed = j;
    }
    first += counts[j];
  }
  if (status == MYB_OK && path != NULL) {
    status = myb_spk_writer_write_file(writer, path);
  } else if (status == MYB_OK) {
    const size_t bytes = myb_spk_writer_size(writer);
    void* buffer = malloc(bytes > 0 ? bytes : 1);
    status = buffer != NULL ? myb_spk_writer_write(writer, buffer, bytes) : MYB_E_NOMEM;
    if (status == MYB_OK) {
      *data = buffer;
      *size = bytes;
    } else {
      free(buffer);
    }
  }
  myb_spk_writer_destroy(writer);
  return status;
}
