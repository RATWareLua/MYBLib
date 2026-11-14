/* SPDX-License-Identifier: BSD-3-Clause */
/* glue_system.c — systems of bodies for the Go binding (glue_tasks.h): a myb_system made from
 * records, its runs through contacts and approaches, its invariants and derivatives, and runs
 * of the Wisdom-Holman integrators. */
#include <stdlib.h>
#include <string.h>

#include "glue_internal.h"

struct mybgo_system {
  myb_system* system;
  size_t bodies;
  mybgo_orientation_store* orientations; /* one a body */
  mybgo_provider* guided;
  size_t guided_count;
  myb_contact_point* points;
  mybgo_go_context* actions; /* one a link */
};

void mybgo_system_destroy(mybgo_system* s) {
  if (s == NULL) {
    return;
  }
  myb_system_destroy(s->system);
  free(s->actions);
  free(s->points);
  free(s->guided);
  free(s->orientations);
  free(s);
}

static myb_status set_body(mybgo_system* s, size_t i, const mybgo_system_body* d) {
  myb_orientation_function orientation = NULL;
  void* context = NULL;
  myb_status status =
      mybgo_orientation_of(&d->orientation, &s->orientations[i], &orientation, &context);
  if (status == MYB_OK && orientation != NULL) {
    status = myb_system_set_orientation(s->system, i, orientation, context);
  }
  if (status == MYB_OK && d->field != NULL) {
    status = myb_system_set_field(s->system, i, d->field, d->degree);
  }
  if (status == MYB_OK && d->has_rotation) {
    status = myb_system_set_rotation(s->system, i, &d->rotation);
  }
  if (status == MYB_OK && d->has_tide) {
    status = myb_system_set_tide(s->system, i, &d->tide);
  }
  if (status == MYB_OK && d->radius != 0.0) {
    status = myb_system_set_radius(s->system, i, d->radius);
  }
  if (status == MYB_OK && d->has_mass_rate) {
    status = myb_system_set_mass_rate(s->system, i, d->mass_rate, d->mass_epoch);
  }
  return status;
}

static myb_status add_link(mybgo_system* s, size_t j, const mybgo_system_link* d,
                           size_t point_count) {
  const int indices = d->a < s->bodies && d->b < s->bodies;
  switch (d->kind) {
    case MYBGO_LINK_INTERACTION:
      return indices ? myb_nbody_set_interaction(myb_system_nbody(s->system), d->a, d->b,
                                                 (myb_interaction)d->level)
                     : MYB_E_DOMAIN;
    case MYBGO_LINK_FIELD:
      return myb_system_add_field(s->system, d->a, d->b);
    case MYBGO_LINK_FIGURES:
      return myb_system_add_figures(s->system, d->a, d->b);
    case MYBGO_LINK_TIDE:
      return myb_system_add_tide(s->system, d->a, d->b, d->c);
    case MYBGO_LINK_DISTORTION:
      return myb_system_add_distortion(s->system, d->a, d->b);
    case MYBGO_LINK_PN_PAIR:
      return myb_system_add_pn_pair(s->system, d->a, d->b, &d->pn);
    case MYBGO_LINK_DRAGGING:
      return myb_system_add_frame_dragging(s->system, d->a, d->b, d->spin);
    case MYBGO_LINK_CONTACT: {
      if (d->point_count == 0 || d->first_point > point_count ||
          d->point_count > point_count - d->first_point) {
        return MYB_E_DOMAIN;
      }
      const myb_system_contact contact = {.surface = d->surface,
                                          .law = d->law,
                                          .points = s->points + d->first_point,
                                          .point_count = d->point_count};
      return myb_system_add_contact(s->system, d->a, d->b, &contact);
    }
    case MYBGO_LINK_ACTION:
      s->actions[j].handle = d->handle;
      return myb_system_add_action(s->system, d->a, mybgo_go_action, &s->actions[j]);
    default:
      return MYB_E_DOMAIN;
  }
}

myb_status mybgo_system_create(const mybgo_system_options* options, const mybgo_system_body* bodies,
                               size_t body_count, const mybgo_system_link* links, size_t link_count,
                               const myb_contact_point* points, size_t point_count,
                               const mybgo_provider* guided, const double* guided_gm,
                               size_t guided_count, mybgo_system** system, int* stage,
                               size_t* item) {
  *system = NULL;
  *stage = MYBGO_STAGE_ENVIRONMENT;
  *item = 0;
  if (body_count == 0) {
    return MYB_E_DOMAIN;
  }
  mybgo_system* s = calloc(1, sizeof *s);
  if (s == NULL) {
    return MYB_E_NOMEM;
  }
  myb_status status = MYB_E_NOMEM;
  double* gm = calloc(body_count, sizeof *gm);
  s->bodies = body_count;
  s->orientations = calloc(body_count, sizeof *s->orientations);
  s->guided = calloc(guided_count > 0 ? guided_count : 1, sizeof *s->guided);
  s->points = calloc(point_count > 0 ? point_count : 1, sizeof *s->points);
  s->actions = calloc(link_count > 0 ? link_count : 1, sizeof *s->actions);
  if (gm == NULL || s->orientations == NULL || s->guided == NULL || s->points == NULL ||
      s->actions == NULL) {
    goto fail;
  }
  for (size_t i = 0; i < body_count; ++i) {
    gm[i] = bodies[i].gm;
  }
  if (point_count > 0) {
    memcpy(s->points, points, point_count * sizeof *points);
  }
  status = myb_system_create(body_count, gm, NULL, &s->system);
  if (status == MYB_OK && options->has_relativity) {
    status = myb_nbody_set_relativity(myb_system_nbody(s->system), &options->relativity);
  }
  if (status == MYB_OK && options->fmm_order != 0) {
    status = myb_nbody_set_approximation(myb_system_nbody(s->system), options->fmm_order,
                                         options->fmm_theta);
  }
  if (status != MYB_OK) {
    goto fail;
  }
  *stage = MYBGO_STAGE_BODY;
  for (size_t i = 0; i < body_count; ++i) {
    *item = i;
    status = set_body(s, i, &bodies[i]);
    if (status != MYB_OK) {
      goto fail;
    }
  }
  *stage = MYBGO_STAGE_TERM;
  for (size_t j = 0; j < link_count; ++j) {
    *item = j;
    status = add_link(s, j, &links[j], point_count);
    if (status != MYB_OK) {
      goto fail;
    }
  }
  *stage = MYBGO_STAGE_PARAMETER;
  for (size_t g = 0; g < guided_count; ++g) {
    *item = g;
    s->guided[g] = guided[g];
    status =
        myb_system_add_guided(s->system, guided_gm[g], mybgo_provider_state, &s->guided[g], NULL);
    if (status != MYB_OK) {
      goto fail;
    }
  }
  s->guided_count = guided_count;
  free(gm);
  *stage = 0;
  *item = 0;
  *system = s;
  return MYB_OK;
fail:
  free(gm);
  mybgo_system_destroy(s);
  return status;
}

size_t mybgo_system_coordinates(const mybgo_system* system) {
  return myb_system_coordinates(system->system);
}

size_t mybgo_system_bodies(const mybgo_system* system) {
  return system->bodies;
}

int mybgo_system_present(const mybgo_system* system, size_t body) {
  return body < system->bodies && myb_system_present(system->system, body);
}

double mybgo_system_gm(const mybgo_system* system, size_t body) {
  return body < system->bodies ? myb_nbody_gm(myb_system_nbody(system->system), body) : 0.0;
}

myb_status mybgo_system_pack(const mybgo_system* s, const double* bodies,
                             const mybgo_rotation_state* rotations, size_t rotation_count,
                             const mybgo_spin_state* spins, size_t spin_count, double* y) {
  const size_t n = myb_system_coordinates(s->system);
  memset(y, 0, 2 * n * sizeof(double));
  for (size_t i = 0; i < s->bodies; ++i) {
    memcpy(y + 3 * i, bodies + 6 * i, 3 * sizeof(double));
    memcpy(y + n + 3 * i, bodies + 6 * i + 3, 3 * sizeof(double));
  }
  myb_status status = MYB_OK;
  for (size_t k = 0; k < rotation_count && status == MYB_OK; ++k) {
    const mybgo_rotation_state* r = &rotations[k];
    status = myb_system_pack_rotation(s->system, r->body, r->q, r->omega,
                                      r->has_core ? r->omega_core : NULL, y);
  }
  for (size_t k = 0; k < spin_count && status == MYB_OK; ++k) {
    status = myb_system_pack_spins(s->system, spins[k].pair, spins[k].chi1, spins[k].chi2, y);
  }
  return status;
}

myb_status mybgo_system_unpack(const mybgo_system* s, const double* y, double* bodies,
                               mybgo_rotation_state* rotations, size_t rotation_count,
                               mybgo_spin_state* spins, size_t spin_count) {
  const size_t n = myb_system_coordinates(s->system);
  for (size_t i = 0; i < s->bodies; ++i) {
    memcpy(bodies + 6 * i, y + 3 * i, 3 * sizeof(double));
    memcpy(bodies + 6 * i + 3, y + n + 3 * i, 3 * sizeof(double));
  }
  myb_status status = MYB_OK;
  for (size_t k = 0; k < rotation_count && status == MYB_OK; ++k) {
    mybgo_rotation_state* r = &rotations[k];
    status = myb_system_unpack_rotation(s->system, r->body, y, r->q, r->omega,
                                        r->has_core ? r->omega_core : NULL);
  }
  for (size_t k = 0; k < spin_count && status == MYB_OK; ++k) {
    status = myb_system_unpack_spins(s->system, spins[k].pair, y, spins[k].chi1, spins[k].chi2);
  }
  return status;
}

/* ---- Runs ------------------------------------------------------------------------------------ */

typedef struct event_list {
  mybgo_system_event* items;
  size_t count;
  size_t capacity;
} event_list;

static myb_status note(event_list* list, const myb_system_event* e, int action) {
  if (list->count == list->capacity) {
    const size_t capacity = list->capacity > 0 ? 2 * list->capacity : 16;
    mybgo_system_event* items = realloc(list->items, capacity * sizeof *items);
    if (items == NULL) {
      return MYB_E_NOMEM;
    }
    list->items = items;
    list->capacity = capacity;
  }
  mybgo_system_event* out = &list->items[list->count++];
  memset(out, 0, sizeof *out);
  out->kind = (int)e->kind;
  out->a = e->a;
  out->b = e->b;
  out->t = e->t;
  out->distance = e->distance;
  out->action = action;
  return MYB_OK;
}

/* What the policy does with a contact; the event noted last receives the outcome. */
static myb_status act(mybgo_system* s, myb_ias15* integrator, const mybgo_system_run_options* o,
                      const myb_system_event* e, mybgo_system_event* noted, int* stop) {
  switch (o->policy) {
    case MYBGO_CONTACT_STOP:
      *stop = 1;
      return MYB_OK;
    case MYBGO_CONTACT_MERGE: {
      const myb_nbody* nbody = myb_system_nbody(s->system);
      const int b_heavier = myb_nbody_gm(nbody, e->b) > myb_nbody_gm(nbody, e->a);
      const size_t into = b_heavier ? e->b : e->a;
      const size_t from = b_heavier ? e->a : e->b;
      noted->survivor = into;
      return myb_system_merge(s->system, integrator, into, from, &noted->merger);
    }
    case MYBGO_CONTACT_BOUNCE:
      return myb_system_bounce(s->system, integrator, e->a, e->b, o->restitution);
    case MYBGO_CONTACT_CONTINUE:
      return MYB_OK;
    default:
      return MYB_E_DOMAIN;
  }
}

myb_status mybgo_system_run(mybgo_system* s, const mybgo_system_run_options* options, double t0,
                            const double* y0, const double* times, size_t count, double* states,
                            mybgo_system_event** events, size_t* event_count, double* t_end,
                            double* y_end, myb_ias15_stats* stats, size_t* failed) {
  *events = NULL;
  *event_count = 0;
  *t_end = t0;
  memset(stats, 0, sizeof *stats);
  const size_t n = myb_system_coordinates(s->system);
  myb_ias15* integrator = NULL;
  myb_status status = myb_system_integrator(s->system, &options->ias15, NULL, &integrator);
  if (status == MYB_OK) {
    status = myb_ias15_start(integrator, myb_system_acceleration, s->system, t0, y0);
  }
  event_list list = {0};
  size_t k = 0;
  int stop = 0;
  /* bound: every pass reaches a time or notes one of at most `limit` events */
  while (status == MYB_OK && !stop && k < count) {
    myb_system_event e;
    status = myb_system_advance(s->system, integrator, times[k], options->watch, &e);
    if (status != MYB_OK) {
      break;
    }
    if (e.kind == MYB_SYSTEM_ARRIVED) {
      memcpy(states + 2 * n * k, myb_ias15_state(integrator), 2 * n * sizeof(double));
      ++k;
      continue;
    }
    if (list.count == options->limit) {
      status = MYB_E_OVERFLOW;
      break;
    }
    const int action = e.kind == MYB_SYSTEM_CONTACT ? options->policy : MYBGO_CONTACT_CONTINUE;
    status = note(&list, &e, action);
    if (status == MYB_OK && e.kind == MYB_SYSTEM_CONTACT) {
      status = act(s, integrator, options, &e, &list.items[list.count - 1], &stop);
    }
  }
  if (integrator != NULL) {
    *t_end = myb_ias15_time(integrator);
    memcpy(y_end, myb_ias15_state(integrator), 2 * n * sizeof(double));
    *stats = myb_ias15_statistics(integrator);
  }
  myb_ias15_destroy(integrator);
  *events = list.items;
  *event_count = list.count;
  return mybgo_end(status, k, failed);
}

myb_status mybgo_system_accelerations(mybgo_system* s, const double* t, const double* y,
                                      size_t count, double* a, size_t* failed) {
  const size_t n = myb_system_coordinates(s->system);
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_system_acceleration(s->system, t[k], y + 2 * n * k, y + 2 * n * k + n, a + n * k);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_system_invariants(mybgo_system* s, const double* y, size_t count,
                                   myb_nbody_invariants* invariants, size_t* failed) {
  const size_t n = myb_system_coordinates(s->system);
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_nbody_invariants_of(myb_system_nbody(s->system), y + 2 * n * k, y + 2 * n * k + n,
                                     &invariants[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_system_pn_invariants(mybgo_system* s, const double* y, size_t count,
                                      myb_nbody_pn_invariants* invariants, size_t* failed) {
  const size_t n = myb_system_coordinates(s->system);
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_nbody_pn_invariants_of(myb_system_nbody(s->system), y + 2 * n * k,
                                        y + 2 * n * k + n, &invariants[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_system_partials(mybgo_system* s, const myb_ias15_options* options, double t0,
                                 const double* y0, double t1, const size_t* bodies, size_t count,
                                 double* y1, double* phi, double* sensitivity) {
  return myb_system_partials(s->system, options, t0, y0, t1, bodies, count, y1, phi, sensitivity);
}

myb_status mybgo_wh_run(const double* gm, size_t massive, size_t test,
                        const myb_wh_options* options, double t0, const double* x, const double* v,
                        const double* times, size_t count, double* xs, double* vs,
                        myb_wh_stats* stats, size_t* failed) {
  memset(stats, 0, sizeof *stats);
  myb_nbody* nbody = NULL;
  myb_wh* wh = NULL;
  myb_status status = myb_nbody_create(massive, test, gm, NULL, &nbody);
  if (status == MYB_OK) {
    status = myb_wh_create(nbody, options, NULL, &wh);
  }
  if (status == MYB_OK) {
    status = myb_wh_start(wh, t0, x, v);
  }
  const size_t values = 3 * (massive + test);
  size_t k = 0;
  for (; k < count && status == MYB_OK; ++k) {
    status = myb_wh_integrate(wh, times[k]);
    if (status == MYB_OK) {
      status = myb_wh_dense(wh, times[k], xs + values * k, vs + values * k);
    }
    if (status != MYB_OK) {
      break;
    }
  }
  if (wh != NULL) {
    *stats = myb_wh_statistics(wh);
  }
  myb_wh_destroy(wh);
  myb_nbody_destroy(nbody);
  return mybgo_end(status, k, failed);
}
