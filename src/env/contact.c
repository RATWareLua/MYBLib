/* SPDX-License-Identifier: BSD-3-Clause */
/* The contact of the object of a force model with the surface of a body of its environment
 * (myb/contact.h): the term and its event. */
#include "myb/contact.h"

#include <math.h>
#include <string.h>

#include "core/base.h"
#include "env/env_internal.h"
#include "num/vec3.h"

/* The object as a sphere, in the frame of the surface: the force of the contact on it in the
 * axes of the model, and what myb_contact_force says of it. */
static myb_status contact_at(const myb_contact_term* c, double t, const double x[3],
                             const double v[3], double a[3], myb_contact_result* result) {
  const myb_contact_point sphere = {.offset = {0.0, 0.0, 0.0}, .radius = c->radius};
  if (c->env == NULL) {
    const myb_status st =
        myb_contact_force(&c->surface, &c->law, &sphere, 1, x, v, NULL, NULL, result);
    if (st == MYB_OK && a != NULL) {
      memcpy(a, result->force, 3 * sizeof *a);
    }
    return st;
  }
  const double* s = NULL;
  myb_status st = myb_internal_env_state(c->env, t, c->body, &s);
  if (st != MYB_OK) {
    return st;
  }
  const myb_orientation* o = NULL;
  st = myb_internal_env_orientation(c->env, t, c->body, &o, NULL);
  if (st != MYB_OK) {
    return st;
  }
  /* (M d, M w + dM d): the position and the velocity relative to the turning surface */
  double state[6];
  for (int i = 0; i < 3; ++i) {
    state[i] = x[i] - s[i];
    state[3 + i] = v[i] - s[3 + i];
  }
  myb_rotation_apply_state(&o->rotation, state, state);
  st = myb_contact_force(&c->surface, &c->law, &sphere, 1, state, state + 3, NULL, NULL, result);
  if (st == MYB_OK && a != NULL) {
    const double* m = o->rotation.m;
    const double* f = result->force;
    for (int i = 0; i < 3; ++i) {
      a[i] = m[i] * f[0] + m[3 + i] * f[1] + m[6 + i] * f[2];
    }
  }
  return st;
}

myb_status myb_term_contact(void* term, double t, const double x[3], const double v[3],
                            double a[3]) {
  const myb_contact_term* c = term;
  MYB_CHECK(c != NULL && x != NULL && v != NULL && a != NULL);
  myb_contact_result result;
  return contact_at(c, t, x, v, a, &result);
}

double myb_event_contact(void* term, double t, const double* y) {
  const myb_contact_term* c = term;
  MYB_CHECK(c != NULL && y != NULL);
  myb_contact_result result;
  if (contact_at(c, t, y, y + 3, NULL, &result) != MYB_OK) {
    return (double)NAN;
  }
  return result.gap;
}
