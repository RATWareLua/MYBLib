/* SPDX-License-Identifier: BSD-3-Clause */
/* A system of bodies with their models (myb/system.h): the levels of the interaction graph
 * beyond the point masses — fields, figures, tides, distortion, post-Newtonian pairs —, the
 * rotation of bodies integrated with the motion, and encounters.
 *
 * One evaluation goes in a fixed order, every sum in the order in which the levels were added:
 * the point masses (myb_nbody), the frames of the bodies, the inertia tensors of the integrated
 * rotations, then fields, figures, tides, pairs, and at the end the equations of rotation,
 * which need the torques gathered on the way. */
#include <math.h> /* sqrt, fabs */
#include <stdckdint.h>
#include <stdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "dyn/nbody_internal.h"
#include "frames/rotation.h"
#include "integ/ias15_internal.h"
#include "myb/attitude.h"
#include "myb/cmath.h"
#include "myb/kin.h"
#include "myb/num.h"
#include "myb/rel.h"
#include "myb/surface.h"
#include "myb/system.h"
#include "myb/time.h"
#include "num/vec3.h"
#include "rel/rel_internal.h"

enum {
  ROTATION_COORDINATES = 9, /* q (4), unused (2), the core (3) */
  SPIN_COORDINATES = 6,
  CORE_OFFSET = 6,
  ROOT_STEPS = 200, /* bisections of an event on a step: more than the bits of a double */
  NO_ROTATION = -1
};

typedef struct body_model {
  myb_orientation_function orientation;
  void* orientation_context;
  const myb_gravity_field* field;
  int degree;   /* of the field in use */
  int rotation; /* index of its integrated rotation, or NO_ROTATION */
  int has_tide;
  myb_system_tide tide;
  double radius; /* of the sphere of encounters; 0: none */
  int present;
  /* a mass that changes with the time: gm = gm_base (1 + gm_rate (t - gm_epoch)) */
  double gm_rate;
  double gm_epoch;
  double gm_base;
  /* of the evaluation in hand */
  int framed;
  double m[9]; /* inertial -> body */
} body_model;

typedef struct rotation_state {
  size_t body;
  myb_system_rotation model;
  /* of the evaluation in hand */
  double norm;     /* |q| of the state */
  double q[4];     /* unit */
  double dq[4];    /* rate of the unit quaternion */
  double omega[3]; /* mantle, in its frame */
  double core[3];  /* core, in the mantle frame */
  double mantle[9];
  double mantle_rate[9]; /* of the tides only */
  double total[9];
  double torque[3];
} rotation_state;

typedef struct pair {
  size_t a;
  size_t b;
} pair;

typedef struct triple {
  size_t body;
  size_t raiser;
  size_t point;
} triple;

typedef struct pn_state {
  size_t a;
  size_t b;
  myb_system_pn_pair model;
  int spins;
} pn_state;

typedef struct action_entry {
  size_t body;
  myb_system_action action;
  void* context;
} action_entry;

/* The spin of a body that drags a point. */
typedef struct drag_entry {
  size_t body;
  size_t point;
  double spin[3];
} drag_entry;

/* A body whose motion is given. */
typedef struct guided_body {
  double gm;
  myb_state_function state;
  void* context;
} guided_body;

enum { GUIDED_SLOTS = 8 }; /* instants whose places are kept: the nodes of a step and one */

struct myb_system {
  myb_allocator allocator;
  size_t bodies;
  myb_nbody* nbody;
  body_model* models;
  double* gm;
  double* point; /* 3 per body: the accelerations of the point masses of the evaluation */
  double* scratch; /* two states, for the encounters: the start of a step, and one to work on */
  size_t scratch_size;
  myb_gravity_work* work;
  int work_degree;
  size_t rotation_count;
  rotation_state rotations[MYB_SYSTEM_MAX_ROTATIONS];
  size_t field_count;
  pair fields[MYB_SYSTEM_MAX_FIELDS]; /* a: the body, b: the point */
  size_t figure_count;
  pair figures[MYB_SYSTEM_MAX_FIGURES];
  size_t tide_count;
  triple tides[MYB_SYSTEM_MAX_TIDES];
  size_t distortion_count;
  pair distortions[MYB_SYSTEM_MAX_DISTORTIONS]; /* a: the body, b: the raiser */
  size_t pn_count;
  pn_state pn[MYB_SYSTEM_MAX_PN_PAIRS];
  size_t action_count;
  action_entry actions[MYB_SYSTEM_MAX_ACTIONS];
  size_t drag_count;
  drag_entry drags[MYB_SYSTEM_MAX_DRAGS];
  /* the bodies with a given motion, and their places at the last instants asked for: an
   * integrator comes back to the same instants several times within a step */
  size_t guided_count;
  size_t guided_capacity;
  guided_body* guided;
  double* guided_places; /* GUIDED_SLOTS x 3 guided_capacity */
  double guided_time[GUIDED_SLOTS];
  int guided_filled[GUIDED_SLOTS];
  int guided_next;
  /* the closest approach reported last: not to be found again at the start of the next step */
  int reported;
  size_t reported_a;
  size_t reported_b;
  double reported_t;
};

/* ---- Creation -------------------------------------------------------------------------------- */

void myb_system_destroy(myb_system* system) {
  if (system == NULL) {
    return;
  }
  const myb_allocator a = system->allocator;
  myb_gravity_work_destroy(system->work);
  myb_nbody_destroy(system->nbody);
  myb_internal_release(&a, system->guided);
  myb_internal_release(&a, system->guided_places);
  myb_internal_release(&a, system->scratch);
  myb_internal_release(&a, system->point);
  myb_internal_release(&a, system->gm);
  myb_internal_release(&a, system->models);
  myb_internal_release(&a, system);
}

myb_status myb_system_create(size_t bodies, const double* gm, const myb_allocator* allocator,
                             myb_system** system) {
  MYB_CHECK(gm != NULL && system != NULL);
  if (bodies == 0) {
    return MYB_E_DOMAIN;
  }
  size_t triple_count;
  if (ckd_mul(&triple_count, bodies, (size_t)3)) {
    return MYB_E_OVERFLOW;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status st = myb_internal_alloc_array(&a, 1, sizeof(myb_system), &block);
  if (st != MYB_OK) {
    return st;
  }
  myb_system* s = block;
  memset(s, 0, sizeof *s);
  s->allocator = a;
  s->bodies = bodies;
  st = myb_nbody_create(bodies, 0, gm, allocator, &s->nbody);
  if (st == MYB_OK) {
    st = myb_internal_alloc_array(&a, bodies, sizeof(body_model), &block);
  }
  if (st == MYB_OK) {
    s->models = block;
    memset(s->models, 0, bodies * sizeof(body_model));
    st = myb_internal_alloc_array(&a, bodies, sizeof(double), &block);
  }
  if (st == MYB_OK) {
    s->gm = block;
    memcpy(s->gm, gm, bodies * sizeof(double));
    st = myb_internal_alloc_array(&a, triple_count, sizeof(double), &block);
  }
  if (st == MYB_OK) {
    s->point = block;
    memset(s->point, 0, triple_count * sizeof(double));
  }
  if (st != MYB_OK) {
    myb_system_destroy(s);
    return st;
  }
  for (size_t b = 0; b < bodies; ++b) {
    s->models[b].rotation = NO_ROTATION;
    s->models[b].present = 1;
  }
  *system = s;
  return MYB_OK;
}

size_t myb_system_bodies(const myb_system* system) {
  MYB_CHECK(system != NULL);
  return system->bodies;
}

myb_nbody* myb_system_nbody(myb_system* system) {
  MYB_CHECK(system != NULL);
  return system->nbody;
}

static size_t spin_pairs(const myb_system* s) {
  size_t n = 0;
  for (size_t k = 0; k < s->pn_count; ++k) {
    n += s->pn[k].spins ? 1u : 0u;
  }
  return n;
}

size_t myb_system_coordinates(const myb_system* system) {
  MYB_CHECK(system != NULL);
  return 3 * system->bodies + ROTATION_COORDINATES * system->rotation_count +
         SPIN_COORDINATES * spin_pairs(system);
}

static size_t rotation_offset(const myb_system* s, int rotation) {
  return 3 * s->bodies + ROTATION_COORDINATES * (size_t)rotation;
}

/* The offset of the spins of pair k, which has them. */
static size_t spin_offset(const myb_system* s, size_t k) {
  size_t before = 0;
  for (size_t i = 0; i < k; ++i) {
    before += s->pn[i].spins ? 1u : 0u;
  }
  return 3 * s->bodies + ROTATION_COORDINATES * s->rotation_count + SPIN_COORDINATES * before;
}

/* ---- The model of a body --------------------------------------------------------------------- */

static int usable(const myb_system* s, size_t body) {
  return body < s->bodies && s->models[body].present;
}

myb_status myb_system_set_orientation(myb_system* system, size_t body,
                                      myb_orientation_function orientation, void* context) {
  MYB_CHECK(system != NULL);
  if (!usable(system, body) || system->models[body].rotation != NO_ROTATION) {
    return MYB_E_DOMAIN;
  }
  system->models[body].orientation = orientation;
  system->models[body].orientation_context = context;
  return MYB_OK;
}

myb_status myb_system_set_field(myb_system* system, size_t body, const myb_gravity_field* field,
                                int degree) {
  MYB_CHECK(system != NULL);
  if (!usable(system, body)) {
    return MYB_E_DOMAIN;
  }
  body_model* m = &system->models[body];
  if (field == NULL) {
    m->field = NULL;
    m->degree = 0;
    return MYB_OK;
  }
  const int all = myb_gravity_field_degree(field);
  if (degree < 0 || degree > all || (degree != 0 && degree < 2)) {
    return MYB_E_DOMAIN;
  }
  const int used = degree == 0 ? all : degree;
  if (used > system->work_degree || system->work == NULL) {
    myb_gravity_work* work = NULL;
    const myb_status st = myb_gravity_work_create(used, &system->allocator, &work);
    if (st != MYB_OK) {
      return st;
    }
    myb_gravity_work_destroy(system->work);
    system->work = work;
    system->work_degree = used;
  }
  m->field = field;
  m->degree = used;
  return MYB_OK;
}

static int names_body_in_tides(const myb_system* s, size_t body) {
  for (size_t k = 0; k < s->tide_count; ++k) {
    if (s->tides[k].body == body) {
      return 1;
    }
  }
  return 0;
}

myb_status myb_system_set_rotation(myb_system* system, size_t body,
                                   const myb_system_rotation* rotation) {
  MYB_CHECK(system != NULL && rotation != NULL);
  const myb_system_rotation* r = rotation;
  if (!usable(system, body) || system->models[body].orientation != NULL ||
      !(system->gm[body] > 0.0) || names_body_in_tides(system, body)) {
    return MYB_E_DOMAIN;
  }
  double sum = r->radius + r->k2 + r->lag + r->mean_motion + r->core[0] + r->core[1] +
               r->core_friction;
  for (int i = 0; i < 3; ++i) {
    sum += r->moments[i];
    if (!(r->moments[i] > 0.0)) {
      return MYB_E_DOMAIN;
    }
  }
  const int core = r->core[0] != 0.0 || r->core[1] != 0.0;
  if (!isfinite(sum) || (r->k2 != 0.0 && !(r->radius > 0.0)) ||
      (core && (!(r->core[0] > 0.0) || !(r->core[1] > 0.0) || !(r->core[0] < r->moments[0]) ||
                !(r->core[1] < r->moments[2])))) {
    return MYB_E_DOMAIN;
  }
  int index = system->models[body].rotation;
  if (index == NO_ROTATION) {
    if (system->rotation_count >= MYB_SYSTEM_MAX_ROTATIONS) {
      return MYB_E_OVERFLOW;
    }
    index = (int)system->rotation_count++;
    system->models[body].rotation = index;
  }
  memset(&system->rotations[index], 0, sizeof system->rotations[index]);
  system->rotations[index].body = body;
  system->rotations[index].model = *r;
  return MYB_OK;
}

myb_status myb_system_set_tide(myb_system* system, size_t body, const myb_system_tide* tide) {
  MYB_CHECK(system != NULL && tide != NULL);
  if (!usable(system, body) || system->models[body].rotation != NO_ROTATION) {
    return MYB_E_DOMAIN;
  }
  double sum = tide->radius + tide->rotation_rate;
  for (int j = 0; j < 3; ++j) {
    sum += tide->k2[j] + tide->orbital_lag[j] + tide->rotation_lag[j];
  }
  if (!isfinite(sum) || !(tide->radius > 0.0)) {
    return MYB_E_DOMAIN;
  }
  system->models[body].tide = *tide;
  system->models[body].has_tide = 1;
  return MYB_OK;
}

myb_status myb_system_set_radius(myb_system* system, size_t body, double radius) {
  MYB_CHECK(system != NULL);
  if (!usable(system, body) || !(radius >= 0.0) || !isfinite(radius)) {
    return MYB_E_DOMAIN;
  }
  system->models[body].radius = radius;
  return MYB_OK;
}

double myb_system_radius(const myb_system* system, size_t body) {
  MYB_CHECK(system != NULL && body < system->bodies);
  return system->models[body].radius;
}

int myb_system_present(const myb_system* system, size_t body) {
  MYB_CHECK(system != NULL && body < system->bodies);
  return system->models[body].present;
}

/* ---- The levels ------------------------------------------------------------------------------ */

static int two_bodies(const myb_system* s, size_t a, size_t b) {
  return usable(s, a) && usable(s, b) && a != b;
}

myb_status myb_system_add_field(myb_system* system, size_t body, size_t point) {
  MYB_CHECK(system != NULL);
  if (!two_bodies(system, body, point) || !(system->gm[body] > 0.0)) {
    return MYB_E_DOMAIN;
  }
  if (system->field_count >= MYB_SYSTEM_MAX_FIELDS) {
    return MYB_E_OVERFLOW;
  }
  system->fields[system->field_count++] = (pair){body, point};
  return MYB_OK;
}

myb_status myb_system_add_figures(myb_system* system, size_t a, size_t b) {
  MYB_CHECK(system != NULL);
  if (!two_bodies(system, a, b) || !(system->gm[a] > 0.0) || !(system->gm[b] > 0.0)) {
    return MYB_E_DOMAIN;
  }
  if (system->figure_count >= MYB_SYSTEM_MAX_FIGURES) {
    return MYB_E_OVERFLOW;
  }
  system->figures[system->figure_count++] = (pair){a, b};
  return MYB_OK;
}

myb_status myb_system_add_tide(myb_system* system, size_t body, size_t raiser, size_t point) {
  MYB_CHECK(system != NULL);
  if (!two_bodies(system, body, raiser) || !two_bodies(system, body, point) ||
      !(system->gm[body] > 0.0) || system->models[body].rotation != NO_ROTATION) {
    return MYB_E_DOMAIN;
  }
  if (system->tide_count >= MYB_SYSTEM_MAX_TIDES) {
    return MYB_E_OVERFLOW;
  }
  system->tides[system->tide_count++] = (triple){body, raiser, point};
  return MYB_OK;
}

myb_status myb_system_add_distortion(myb_system* system, size_t body, size_t raiser) {
  MYB_CHECK(system != NULL);
  if (!two_bodies(system, body, raiser) || system->models[body].rotation == NO_ROTATION) {
    return MYB_E_DOMAIN;
  }
  if (system->distortion_count >= MYB_SYSTEM_MAX_DISTORTIONS) {
    return MYB_E_OVERFLOW;
  }
  system->distortions[system->distortion_count++] = (pair){body, raiser};
  return MYB_OK;
}

myb_status myb_system_add_pn_pair(myb_system* system, size_t a, size_t b,
                                  const myb_system_pn_pair* pair_model) {
  MYB_CHECK(system != NULL && pair_model != NULL);
  if (!two_bodies(system, a, b) || !(system->gm[a] > 0.0) || !(system->gm[b] > 0.0) ||
      (pair_model->orders & ~(unsigned)MYB_PN_ALL) != 0u || pair_model->orders == 0u) {
    return MYB_E_DOMAIN;
  }
  myb_system_pn_pair model = *pair_model;
  if (!(model.c > 0.0)) {
    myb_ppn ppn;
    if (model.c != 0.0 || !myb_internal_nbody_relativity(system->nbody, &ppn)) {
      return MYB_E_DOMAIN;
    }
    model.c = ppn.c;
  }
  if (!isfinite(model.c + model.r0 + model.kappa1 + model.kappa2) ||
      ((model.orders & MYB_PN_3) && !(model.r0 > 0.0)) ||
      ((model.orders & MYB_PN_SS) && (!(model.kappa1 >= 0.0) || !(model.kappa2 >= 0.0)))) {
    return MYB_E_DOMAIN;
  }
  if (system->pn_count >= MYB_SYSTEM_MAX_PN_PAIRS) {
    return MYB_E_OVERFLOW;
  }
  pn_state* p = &system->pn[system->pn_count++];
  p->a = a;
  p->b = b;
  p->model = model;
  p->spins = (model.orders & (MYB_PN_SO | MYB_PN_SS)) != 0u;
  return MYB_OK;
}

myb_status myb_system_add_action(myb_system* system, size_t body, myb_system_action action,
                                 void* context) {
  MYB_CHECK(system != NULL && action != NULL);
  if (!usable(system, body)) {
    return MYB_E_DOMAIN;
  }
  if (system->action_count >= MYB_SYSTEM_MAX_ACTIONS) {
    return MYB_E_OVERFLOW;
  }
  system->actions[system->action_count++] = (action_entry){body, action, context};
  return MYB_OK;
}

myb_status myb_system_light_action(void* light, double t, const myb_system_view* view,
                                   double acceleration[3], double torque[3]) {
  (void)t;
  MYB_CHECK(light != NULL && view != NULL && acceleration != NULL && torque != NULL);
  const myb_system_light* l = light;
  if (view->frame == NULL || l->shape == NULL || l->source >= view->bodies ||
      l->source == view->body || !(l->flux >= 0.0) || !(l->distance > 0.0) || !(l->c > 0.0) ||
      !(l->mass > 0.0) || !isfinite(l->flux + l->distance + l->c + l->mass)) {
    return MYB_E_DOMAIN;
  }
  double d[3];
  for (int i = 0; i < 3; ++i) {
    d[i] = view->x[3 * view->body + (size_t)i] - view->x[3 * l->source + (size_t)i];
  }
  const double r = myb_vec3_norm(d);
  if (!(r > 0.0) || !isfinite(r)) {
    return MYB_E_DOMAIN;
  }
  double u[3];
  myb_internal_mat3_vec(view->frame, d, u);
  const double length = myb_vec3_norm(u);
  for (int i = 0; i < 3; ++i) {
    u[i] /= length;
  }
  const double pressure = l->flux * (l->distance / r) * (l->distance / r) / l->c;
  double force[3];
  double moment[3];
  const myb_status st =
      myb_force_radiation_mesh(l->shape, &l->surface, pressure, u, l->shadows, force, moment);
  if (st != MYB_OK) {
    return st;
  }
  /* the force into the inertial axes: the transpose of the frame */
  for (int i = 0; i < 3; ++i) {
    acceleration[i] = (view->frame[i] * force[0] + view->frame[3 + i] * force[1] +
                       view->frame[6 + i] * force[2]) / l->mass;
    torque[i] = moment[i] / l->mass;
  }
  return MYB_OK;
}

/* ---- The state ------------------------------------------------------------------------------- */

myb_status myb_system_pack_rotation(const myb_system* system, size_t body, const double q[4],
                                    const double omega[3], const double omega_core[3],
                                    double* y) {
  MYB_CHECK(system != NULL && q != NULL && omega != NULL && y != NULL);
  if (body >= system->bodies || system->models[body].rotation == NO_ROTATION) {
    return MYB_E_DOMAIN;
  }
  double unit[4] = {q[0], q[1], q[2], q[3]};
  if (!(myb_quat_normalize(unit) > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const size_t n = myb_system_coordinates(system);
  const size_t at = rotation_offset(system, system->models[body].rotation);
  double* x = y + at;
  double* v = y + n + at;
  memset(x, 0, ROTATION_COORDINATES * sizeof *x);
  memset(v, 0, ROTATION_COORDINATES * sizeof *v);
  memcpy(x, unit, sizeof unit);
  myb_quat_derivative(unit, omega, v);
  if (omega_core != NULL) {
    memcpy(v + CORE_OFFSET, omega_core, 3 * sizeof *v);
  }
  return MYB_OK;
}

/* The unit quaternion, its rate and the angular velocity from the state of a rotation. */
static int read_rotation(const double* x, const double* v, double q[4], double dq[4],
                         double omega[3], double* norm) {
  memcpy(q, x, 4 * sizeof *q);
  const double n = myb_quat_normalize(q);
  if (!(n > 0.0) || !isfinite(n)) {
    return 0;
  }
  for (int i = 0; i < 4; ++i) {
    dq[i] = v[i] / n;
  }
  /* (0, omega) = 2 conj(q) dq */
  double conj[4];
  double w[4];
  myb_quat_conjugate(q, conj);
  myb_quat_multiply(conj, dq, w);
  for (int i = 0; i < 3; ++i) {
    omega[i] = 2.0 * w[1 + i];
  }
  *norm = n;
  return 1;
}

myb_status myb_system_unpack_rotation(const myb_system* system, size_t body, const double* y,
                                      double q[4], double omega[3], double omega_core[3]) {
  MYB_CHECK(system != NULL && y != NULL);
  if (body >= system->bodies || system->models[body].rotation == NO_ROTATION) {
    return MYB_E_DOMAIN;
  }
  const size_t n = myb_system_coordinates(system);
  const size_t at = rotation_offset(system, system->models[body].rotation);
  double unit[4];
  double dq[4];
  double w[3];
  double norm = 0.0;
  if (!read_rotation(y + at, y + n + at, unit, dq, w, &norm)) {
    return MYB_E_DOMAIN;
  }
  if (q != NULL) {
    memcpy(q, unit, sizeof unit);
  }
  if (omega != NULL) {
    memcpy(omega, w, sizeof w);
  }
  if (omega_core != NULL) {
    memcpy(omega_core, y + n + at + CORE_OFFSET, 3 * sizeof *omega_core);
  }
  return MYB_OK;
}

myb_status myb_system_pack_spins(const myb_system* system, size_t pair_index,
                                 const double chi1[3], const double chi2[3], double* y) {
  MYB_CHECK(system != NULL && chi1 != NULL && chi2 != NULL && y != NULL);
  if (pair_index >= system->pn_count || !system->pn[pair_index].spins) {
    return MYB_E_DOMAIN;
  }
  const size_t n = myb_system_coordinates(system);
  const size_t at = spin_offset(system, pair_index);
  memset(y + at, 0, SPIN_COORDINATES * sizeof *y);
  memcpy(y + n + at, chi1, 3 * sizeof *y);
  memcpy(y + n + at + 3, chi2, 3 * sizeof *y);
  return MYB_OK;
}

myb_status myb_system_unpack_spins(const myb_system* system, size_t pair_index, const double* y,
                                   double chi1[3], double chi2[3]) {
  MYB_CHECK(system != NULL && chi1 != NULL && chi2 != NULL && y != NULL);
  if (pair_index >= system->pn_count || !system->pn[pair_index].spins) {
    return MYB_E_DOMAIN;
  }
  const size_t n = myb_system_coordinates(system);
  const size_t at = spin_offset(system, pair_index);
  memcpy(chi1, y + n + at, 3 * sizeof *y);
  memcpy(chi2, y + n + at + 3, 3 * sizeof *y);
  return MYB_OK;
}

myb_status myb_system_integrator(const myb_system* system, const myb_ias15_options* options,
                                 const myb_allocator* allocator, myb_ias15** integrator) {
  MYB_CHECK(system != NULL && integrator != NULL);
  size_t groups[1 + MYB_SYSTEM_MAX_ROTATIONS + MYB_SYSTEM_MAX_PN_PAIRS];
  size_t count = 0;
  groups[count++] = 3 * system->bodies;
  for (size_t k = 0; k < system->rotation_count; ++k) {
    groups[count++] = ROTATION_COORDINATES;
  }
  for (size_t k = 0; k < system->pn_count; ++k) {
    if (system->pn[k].spins) {
      groups[count++] = SPIN_COORDINATES;
    }
  }
  myb_ias15_options o = {0};
  if (options != NULL) {
    o = *options;
  }
  o.block = 3;
  o.groups = groups;
  o.group_count = count;
  return myb_ias15_create(myb_system_coordinates(system), &o, allocator, integrator);
}

/* ---- One evaluation: frames ------------------------------------------------------------------ */

static void transposed_vec(const double m[9], const double v[3], double out[3]) {
  const double t[3] = {m[0] * v[0] + m[3] * v[1] + m[6] * v[2],
                       m[1] * v[0] + m[4] * v[1] + m[7] * v[2],
                       m[2] * v[0] + m[5] * v[1] + m[8] * v[2]};
  memcpy(out, t, sizeof t);
}

static void relative(const double* x, size_t from, size_t to, double d[3]) {
  for (int i = 0; i < 3; ++i) {
    d[i] = x[3 * to + (size_t)i] - x[3 * from + (size_t)i];
  }
}

static myb_status frames(myb_system* s, double t, const double* x, const double* v) {
  for (size_t b = 0; b < s->bodies; ++b) {
    body_model* m = &s->models[b];
    m->framed = 0;
    if (m->orientation != NULL && m->present) {
      myb_orientation o;
      const myb_status st = m->orientation(m->orientation_context, myb_time_from_seconds(t), &o);
      if (st != MYB_OK) {
        return st;
      }
      memcpy(m->m, o.rotation.m, sizeof m->m);
      m->framed = 1;
    }
  }
  for (size_t k = 0; k < s->rotation_count; ++k) {
    rotation_state* r = &s->rotations[k];
    const size_t at = rotation_offset(s, (int)k);
    if (!read_rotation(x + at, v + at, r->q, r->dq, r->omega, &r->norm)) {
      return MYB_E_DOMAIN;
    }
    memcpy(r->core, v + at + CORE_OFFSET, sizeof r->core);
    memset(r->torque, 0, sizeof r->torque);
    body_model* m = &s->models[r->body];
    myb_quat_to_matrix(r->q, m->m);
    m->framed = 1;
  }
  return MYB_OK;
}

/* The position and velocity of body `of` relative to body `from` a time `lag` before the
 * evaluation, from the motion at it: a Taylor series to the third order with the acceleration
 * of the point masses and the jerk of the two bodies alone. 0 if they coincide. */
static int earlier(const myb_system* s, const double* x, const double* v, size_t from, size_t of,
                   double lag, double r[3], double rv[3]) {
  double d[3];
  double dv[3];
  double da[3];
  relative(x, from, of, d);
  relative(v, from, of, dv);
  relative(s->point, from, of, da);
  const double r2 = myb_vec3_dot(d, d);
  if (!(r2 > 0.0) || !isfinite(r2)) {
    return 0;
  }
  const double rr = sqrt(r2);
  const double mu = s->gm[from] + s->gm[of];
  const double k3 = mu / (r2 * rr);
  const double k5 = 3.0 * mu * myb_vec3_dot(d, dv) / (r2 * r2 * rr);
  const double h2 = 0.5 * lag * lag;
  const double h3 = lag * lag * lag / 6.0;
  for (int i = 0; i < 3; ++i) {
    const double jerk = -(k3 * dv[i] - k5 * d[i]);
    r[i] = d[i] - lag * dv[i] + h2 * da[i] - h3 * jerk;
    if (rv != NULL) {
      rv[i] = dv[i] - lag * da[i] + h2 * jerk;
    }
  }
  return 1;
}

/* ---- One evaluation: the inertia of an integrated rotation ----------------------------------- */

static double distortion_scale(const rotation_state* r) {
  const double radius = r->model.radius;
  return r->model.k2 * radius * radius * radius * radius * radius;
}

static myb_status distort(myb_system* s, rotation_state* r, const double* x, const double* v) {
  const myb_system_rotation* model = &r->model;
  memset(r->mantle, 0, sizeof r->mantle);
  memset(r->mantle_rate, 0, sizeof r->mantle_rate);
  for (int i = 0; i < 3; ++i) {
    r->mantle[4 * i] = model->moments[i] - (i == 2 ? model->core[1] : model->core[0]);
  }
  if (model->k2 != 0.0) {
    const double k2r5 = distortion_scale(r);
    /* the mantle frame a lag ago: turned back by the angle of the lag */
    const double back[3] = {-r->omega[0] * model->lag, -r->omega[1] * model->lag,
                            -r->omega[2] * model->lag};
    double qe[4];
    double qd[4];
    double md[9];
    myb_quat_from_rotation_vector(back, qe);
    myb_quat_multiply(r->q, qe, qd);
    myb_quat_to_matrix(qd, md);
    for (size_t k = 0; k < s->distortion_count; ++k) {
      if (s->distortions[k].a != r->body) {
        continue;
      }
      const size_t raiser = s->distortions[k].b;
      double ri[3];
      double rvi[3];
      if (!earlier(s, x, v, r->body, raiser, model->lag, ri, rvi)) {
        return MYB_E_DOMAIN;
      }
      double rb[3];
      double rvb[3];
      double wxr[3];
      myb_internal_mat3_vec(md, ri, rb);
      myb_internal_mat3_vec(md, rvi, rvb);
      myb_vec3_cross(r->omega, rb, wxr);
      for (int i = 0; i < 3; ++i) {
        rvb[i] -= wxr[i]; /* the rate of the components in the turning frame */
      }
      const double r2 = myb_vec3_dot(rb, rb);
      const double r5 = r2 * r2 * sqrt(r2);
      const double rdot = myb_vec3_dot(rb, rvb);
      const double ft = -k2r5 * (s->gm[raiser] / s->gm[r->body]);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          const double d = i == j ? 1.0 : 0.0;
          r->mantle[3 * i + j] += ft * (rb[i] * rb[j] - d * r2 / 3.0) / r5;
          r->mantle_rate[3 * i + j] +=
              ft * ((rvb[i] * rb[j] + rb[i] * rvb[j]) / r5 -
                    5.0 * rdot * rb[i] * rb[j] / (r5 * r2) + d * rdot / r5);
        }
      }
    }
    const double fs = k2r5 / (3.0 * s->gm[r->body]);
    const double w2 = myb_vec3_dot(r->omega, r->omega);
    const double n2 = model->mean_motion * model->mean_motion;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const double d = i == j ? 1.0 : 0.0;
        r->mantle[3 * i + j] += fs * (r->omega[i] * r->omega[j] - d * (w2 - n2) / 3.0 -
                                      (i == 2 && j == 2 ? n2 : 0.0));
      }
    }
  }
  memcpy(r->total, r->mantle, sizeof r->total);
  r->total[0] += model->core[0];
  r->total[4] += model->core[0];
  r->total[8] += model->core[1];
  return MYB_OK;
}

/* G Q of a body in its own frame: from the inertia tensor of an integrated rotation, or from
 * the degree 2 of its field. */
static myb_status quadrupole_of(const myb_system* s, size_t body, double gq[9]) {
  const body_model* m = &s->models[body];
  if (m->rotation != NO_ROTATION) {
    const double* it = s->rotations[m->rotation].total;
    const double trace = it[0] + it[4] + it[8];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        gq[3 * i + j] = s->gm[body] * ((i == j ? trace : 0.0) - 3.0 * it[3 * i + j]);
      }
    }
    return MYB_OK;
  }
  if (m->field == NULL) {
    return MYB_E_DOMAIN;
  }
  double c[3];
  double sn[3];
  for (int k = 0; k < 3; ++k) {
    const myb_status st = myb_gravity_field_coefficient(m->field, 2, k, &c[k], &sn[k]);
    if (st != MYB_OK) {
      return st;
    }
  }
  myb_quadrupole_from_coefficients(s->gm[body], myb_gravity_field_radius(m->field), c, sn + 1,
                                   gq);
  return MYB_OK;
}

/* ---- One evaluation: the levels -------------------------------------------------------------- */

/* The acceleration (body frame) that the figure of `body` gives a point at rb. */
static myb_status figure_acceleration(myb_system* s, size_t body, const double rb[3],
                                      double ab[3]) {
  const body_model* m = &s->models[body];
  const double r2 = myb_vec3_dot(rb, rb);
  if (!(r2 > 0.0) || !isfinite(r2)) {
    return MYB_E_DOMAIN;
  }
  if (m->rotation == NO_ROTATION) {
    if (m->field == NULL) {
      return MYB_E_DOMAIN;
    }
    return myb_gravity_field_eval(m->field, s->work, rb, 2, m->degree, NULL, ab);
  }
  /* the degree 2 from the inertia tensor: U = G m (n . q n) / (2 r^3), q = tr(I) 1 - 3 I */
  double gq[9];
  myb_status st = quadrupole_of(s, body, gq);
  if (st != MYB_OK) {
    return st;
  }
  double qr[3];
  myb_internal_mat3_vec(gq, rb, qr);
  const double r5 = r2 * r2 * sqrt(r2);
  const double radial = 2.5 * myb_vec3_dot(rb, qr) / (r5 * r2);
  for (int i = 0; i < 3; ++i) {
    ab[i] = qr[i] / r5 - radial * rb[i];
  }
  if (m->field != NULL && m->degree >= 3) {
    double higher[3];
    st = myb_gravity_field_eval(m->field, s->work, rb, 3, m->degree, NULL, higher);
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 3; ++i) {
      ab[i] += higher[i];
    }
  }
  return MYB_OK;
}

static myb_status apply_fields(myb_system* s, const double* x, double* a) {
  for (size_t k = 0; k < s->field_count; ++k) {
    const size_t body = s->fields[k].a;
    const size_t point = s->fields[k].b;
    const body_model* m = &s->models[body];
    if (!m->framed) {
      return MYB_E_DOMAIN;
    }
    double d[3];
    double rb[3];
    double ab[3];
    double ai[3];
    relative(x, body, point, d);
    myb_internal_mat3_vec(m->m, d, rb);
    const myb_status st = figure_acceleration(s, body, rb, ab);
    if (st != MYB_OK) {
      return st;
    }
    transposed_vec(m->m, ab, ai);
    const double ratio = s->gm[point] / s->gm[body];
    for (int i = 0; i < 3; ++i) {
      a[3 * point + (size_t)i] += ai[i];
      a[3 * body + (size_t)i] -= ratio * ai[i];
    }
    if (m->rotation != NO_ROTATION) {
      double torque[3];
      myb_torque_field(ratio, rb, ab, torque);
      for (int i = 0; i < 3; ++i) {
        s->rotations[m->rotation].torque[i] += torque[i];
      }
    }
  }
  return MYB_OK;
}

/* out = m q m^T: a tensor of the inertial frame in the frame of m. */
static void into_frame(const double m[9], const double q[9], double out[9]) {
  double mt[9];
  myb_internal_mat3_transpose(m, mt);
  myb_quadrupole_rotate(mt, q, out);
}

static myb_status apply_figures(myb_system* s, const double* x, double* a) {
  for (size_t k = 0; k < s->figure_count; ++k) {
    const size_t ia = s->figures[k].a;
    const size_t ib = s->figures[k].b;
    const body_model* ma = &s->models[ia];
    const body_model* mb = &s->models[ib];
    if (!ma->framed || !mb->framed) {
      return MYB_E_DOMAIN;
    }
    double qa_body[9];
    double qb_body[9];
    double qa[9];
    double qb[9];
    myb_status st = quadrupole_of(s, ia, qa_body);
    if (st == MYB_OK) {
      st = quadrupole_of(s, ib, qb_body);
    }
    if (st != MYB_OK) {
      return st;
    }
    myb_quadrupole_rotate(ma->m, qa_body, qa);
    myb_quadrupole_rotate(mb->m, qb_body, qb);
    double r[3];
    relative(x, ia, ib, r);
    st = myb_figure_figure(s->gm[ia], qa, s->gm[ib], qb, r, a + 3 * ia, a + 3 * ib);
    if (st != MYB_OK) {
      return st;
    }
    /* the torque on each whose rotation is integrated, in its own frame */
    for (int side = 0; side < 2; ++side) {
      const size_t self = side == 0 ? ia : ib;
      const body_model* m = side == 0 ? ma : mb;
      if (m->rotation == NO_ROTATION) {
        continue;
      }
      double other[9];
      double rs[3];
      double towards[3];
      double torque[3];
      into_frame(m->m, side == 0 ? qb : qa, other);
      for (int i = 0; i < 3; ++i) {
        towards[i] = side == 0 ? r[i] : -r[i];
      }
      myb_internal_mat3_vec(m->m, towards, rs);
      st = myb_torque_figure_figure(s->gm[self], side == 0 ? qa_body : qb_body, other, rs,
                                    torque);
      if (st != MYB_OK) {
        return st;
      }
      for (int i = 0; i < 3; ++i) {
        s->rotations[m->rotation].torque[i] += torque[i];
      }
    }
  }
  return MYB_OK;
}

/* out = v turned about the unit axis p by the angle (right-handed). */
static void turn(const double p[3], double angle, const double v[3], double out[3]) {
  const double c = myb_cos(angle);
  const double sn = myb_sin(angle);
  double pxv[3];
  myb_vec3_cross(p, v, pxv);
  const double along = myb_vec3_dot(p, v) * (1.0 - c);
  for (int i = 0; i < 3; ++i) {
    out[i] = c * v[i] + sn * pxv[i] + along * p[i];
  }
}

static myb_status apply_tides(myb_system* s, const double* x, const double* v, double* a) {
  for (size_t k = 0; k < s->tide_count; ++k) {
    const triple* e = &s->tides[k];
    const body_model* m = &s->models[e->body];
    if (!m->framed || !m->has_tide) {
      return MYB_E_DOMAIN;
    }
    const myb_system_tide* tide = &m->tide;
    const double pole[3] = {m->m[6], m->m[7], m->m[8]};
    double d[3];
    relative(x, e->body, e->point, d);
    const double r2 = myb_vec3_dot(d, d);
    if (!(r2 > 0.0) || !isfinite(r2)) {
      return MYB_E_DOMAIN;
    }
    const double r5 = r2 * r2 * sqrt(r2);
    const double z = myb_vec3_dot(d, pole);
    double rho[3];
    for (int i = 0; i < 3; ++i) {
      rho[i] = d[i] - z * pole[i];
    }
    const double radius5 =
        tide->radius * tide->radius * tide->radius * tide->radius * tide->radius;
    double acc[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < 3; ++j) {
      if (tide->k2[j] == 0.0) {
        continue;
      }
      double rs[3];
      if (!earlier(s, x, v, e->body, e->raiser, tide->orbital_lag[j], rs, NULL)) {
        return MYB_E_DOMAIN;
      }
      if (j > 0) {
        turn(pole, tide->rotation_rate * tide->rotation_lag[j], rs, rs);
      }
      const double s2 = myb_vec3_dot(rs, rs);
      const double s5 = s2 * s2 * sqrt(s2);
      const double zs = myb_vec3_dot(rs, pole);
      double rhos[3];
      for (int i = 0; i < 3; ++i) {
        rhos[i] = rs[i] - zs * pole[i];
      }
      const double cross = myb_vec3_dot(rho, rhos);
      double w = 0.0;
      double grad[3];
      if (j == 0) {
        const double f = 3.0 * zs * zs - s2;
        w = 0.5 * (3.0 * z * z - r2) * f;
        for (int i = 0; i < 3; ++i) {
          grad[i] = f * (3.0 * z * pole[i] - d[i]);
        }
      } else if (j == 1) {
        w = 6.0 * z * zs * cross;
        for (int i = 0; i < 3; ++i) {
          grad[i] = 6.0 * zs * (cross * pole[i] + z * rhos[i]);
        }
      } else {
        const double rho2 = myb_vec3_dot(rho, rho);
        const double rhos2 = myb_vec3_dot(rhos, rhos);
        w = 1.5 * (2.0 * cross * cross - rho2 * rhos2);
        for (int i = 0; i < 3; ++i) {
          grad[i] = 3.0 * (2.0 * cross * rhos[i] - rhos2 * rho[i]);
        }
      }
      const double scale = 0.5 * s->gm[e->raiser] * radius5 * tide->k2[j] / (s5 * r5);
      const double radial = 5.0 * w / r2;
      for (int i = 0; i < 3; ++i) {
        acc[i] += scale * (grad[i] - radial * d[i]);
      }
    }
    const double ratio = s->gm[e->point] / s->gm[e->body];
    for (int i = 0; i < 3; ++i) {
      a[3 * e->point + (size_t)i] += acc[i];
      a[3 * e->body + (size_t)i] -= ratio * acc[i];
    }
  }
  return MYB_OK;
}

static myb_status apply_pairs(myb_system* s, const double* x, const double* v, double* a) {
  for (size_t k = 0; k < s->pn_count; ++k) {
    const pn_state* p = &s->pn[k];
    const myb_pn_binary binary = {.gm1 = s->gm[p->a], .gm2 = s->gm[p->b], .c = p->model.c,
                                  .orders = p->model.orders, .r0 = p->model.r0,
                                  .kappa1 = p->model.kappa1, .kappa2 = p->model.kappa2};
    double xr[3];
    double vr[3];
    relative(x, p->b, p->a, xr); /* of a with respect to b */
    relative(v, p->b, p->a, vr);
    double acc[3];
    myb_status st = MYB_OK;
    if (p->spins) {
      const size_t at = spin_offset(s, k);
      st = myb_internal_pn_binary_terms(&binary, xr, vr, v + at, v + at + 3, acc, a + at,
                                        a + at + 3);
    } else {
      st = myb_internal_pn_binary_terms(&binary, xr, vr, NULL, NULL, acc, NULL, NULL);
    }
    if (st != MYB_OK) {
      return st;
    }
    const double total = binary.gm1 + binary.gm2;
    for (int i = 0; i < 3; ++i) {
      a[3 * p->a + (size_t)i] += binary.gm2 / total * acc[i];
      a[3 * p->b + (size_t)i] -= binary.gm1 / total * acc[i];
    }
  }
  return MYB_OK;
}

static myb_status apply_actions(myb_system* s, double t, const double* x, const double* v,
                                double* a) {
  for (size_t k = 0; k < s->action_count; ++k) {
    const action_entry* e = &s->actions[k];
    const body_model* m = &s->models[e->body];
    if (!m->present) {
      continue;
    }
    const myb_system_view view = {.body = e->body, .bodies = s->bodies, .x = x, .v = v,
                                  .gm = s->gm, .frame = m->framed ? m->m : NULL};
    double acceleration[3] = {0.0, 0.0, 0.0};
    double torque[3] = {0.0, 0.0, 0.0};
    const myb_status st = e->action(e->context, t, &view, acceleration, torque);
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 3; ++i) {
      a[3 * e->body + (size_t)i] += acceleration[i];
      if (m->rotation != NO_ROTATION) {
        s->rotations[m->rotation].torque[i] += torque[i];
      }
    }
  }
  return MYB_OK;
}

myb_status myb_system_add_frame_dragging(myb_system* system, size_t body, size_t point,
                                         const double spin[3]) {
  MYB_CHECK(system != NULL && spin != NULL);
  myb_ppn ppn;
  if (!usable(system, body) || !usable(system, point) || body == point ||
      !myb_vec3_finite(spin) || !myb_internal_nbody_relativity(system->nbody, &ppn)) {
    return MYB_E_DOMAIN;
  }
  if (system->drag_count >= MYB_SYSTEM_MAX_DRAGS) {
    return MYB_E_OVERFLOW;
  }
  drag_entry* e = &system->drags[system->drag_count++];
  e->body = body;
  e->point = point;
  memcpy(e->spin, spin, sizeof e->spin);
  return MYB_OK;
}

/* The points dragged by the spins of bodies, and the bodies pushed back. */
static myb_status apply_drags(myb_system* s, const double* x, const double* v, double* a) {
  if (s->drag_count == 0) {
    return MYB_OK;
  }
  myb_ppn ppn;
  if (!myb_internal_nbody_relativity(s->nbody, &ppn)) {
    return MYB_E_DOMAIN;
  }
  for (size_t k = 0; k < s->drag_count; ++k) {
    const drag_entry* e = &s->drags[k];
    if (!s->models[e->body].present || !s->models[e->point].present) {
      continue;
    }
    double r[3];
    double w[3];
    double pull[3];
    relative(x, e->body, e->point, r);
    relative(v, e->body, e->point, w);
    const myb_status st = myb_rel_lense_thirring(&ppn, s->gm[e->body], e->spin, r, w, pull);
    if (st != MYB_OK) {
      return st;
    }
    const double ratio = s->gm[e->body] > 0.0 ? s->gm[e->point] / s->gm[e->body] : 0.0;
    for (int i = 0; i < 3; ++i) {
      a[3 * e->point + (size_t)i] += pull[i];
      a[3 * e->body + (size_t)i] -= ratio * pull[i];
    }
  }
  return MYB_OK;
}

myb_status myb_system_set_mass_rate(myb_system* system, size_t body, double rate,
                                    double epoch) {
  MYB_CHECK(system != NULL);
  if (body >= system->bodies || !system->models[body].present || !isfinite(rate) ||
      !isfinite(epoch)) {
    return MYB_E_DOMAIN;
  }
  body_model* m = &system->models[body];
  if (m->gm_rate == 0.0) {
    m->gm_base = system->gm[body]; /* what the body was given is its mass at the epoch */
  }
  m->gm_rate = rate;
  m->gm_epoch = epoch;
  return MYB_OK;
}

/* ---- Bodies with a given motion ----------------------------------------------------------------- */

myb_status myb_system_add_guided(myb_system* system, double gm, myb_state_function state,
                                 void* context, size_t* index) {
  MYB_CHECK(system != NULL && state != NULL);
  myb_system* s = system;
  if (!(gm >= 0.0) || !isfinite(gm)) {
    return MYB_E_DOMAIN;
  }
  if (s->guided_count == s->guided_capacity) {
    const size_t capacity = s->guided_capacity == 0 ? 64 : 2 * s->guided_capacity;
    size_t places;
    if (ckd_mul(&places, capacity, (size_t)(3 * GUIDED_SLOTS))) {
      return MYB_E_OVERFLOW;
    }
    void* bodies = NULL;
    void* room = NULL;
    myb_status st = myb_internal_alloc_array(&s->allocator, capacity, sizeof(guided_body), &bodies);
    if (st == MYB_OK) {
      st = myb_internal_alloc_array(&s->allocator, places, sizeof(double), &room);
    }
    if (st != MYB_OK) {
      myb_internal_release(&s->allocator, bodies);
      return st;
    }
    if (s->guided_count > 0) {
      memcpy(bodies, s->guided, s->guided_count * sizeof(guided_body));
    }
    myb_internal_release(&s->allocator, s->guided);
    myb_internal_release(&s->allocator, s->guided_places);
    s->guided = bodies;
    s->guided_places = room;
    s->guided_capacity = capacity;
  }
  s->guided[s->guided_count] = (guided_body){gm, state, context};
  if (index != NULL) {
    *index = s->guided_count;
  }
  ++s->guided_count;
  memset(s->guided_filled, 0, sizeof s->guided_filled); /* the places kept lack the new body */
  return MYB_OK;
}

size_t myb_system_guided(const myb_system* system) {
  MYB_CHECK(system != NULL);
  return system->guided_count;
}

/* The pull of the bodies with a given motion on every body of the system. */
static myb_status apply_guided(myb_system* s, double t, const double* x, double* a) {
  if (s->guided_count == 0) {
    return MYB_OK;
  }
  int slot = -1;
  for (int i = 0; i < GUIDED_SLOTS; ++i) {
    if (s->guided_filled[i] && s->guided_time[i] == t) {
      slot = i;
    }
  }
  if (slot < 0) {
    slot = s->guided_next;
    s->guided_next = (slot + 1) % GUIDED_SLOTS;
    s->guided_filled[slot] = 0;
    double* places = s->guided_places + (size_t)slot * 3 * s->guided_capacity;
    for (size_t k = 0; k < s->guided_count; ++k) {
      double state[6];
      const myb_status st =
          s->guided[k].state(s->guided[k].context, myb_time_from_seconds(t), state);
      if (st != MYB_OK) {
        return st;
      }
      if (!myb_vec3_finite(state)) {
        return MYB_E_DOMAIN;
      }
      memcpy(places + 3 * k, state, 3 * sizeof *places);
    }
    s->guided_time[slot] = t;
    s->guided_filled[slot] = 1;
  }
  const double* places = s->guided_places + (size_t)slot * 3 * s->guided_capacity;
  for (size_t b = 0; b < s->bodies; ++b) {
    if (!s->models[b].present) {
      continue;
    }
    double pull[3] = {0.0, 0.0, 0.0};
    for (size_t k = 0; k < s->guided_count; ++k) {
      const double d[3] = {places[3 * k] - x[3 * b], places[3 * k + 1] - x[3 * b + 1],
                           places[3 * k + 2] - x[3 * b + 2]};
      const double r2 = myb_vec3_dot(d, d);
      if (!(r2 > 0.0)) {
        return MYB_E_DOMAIN;
      }
      const double f = s->guided[k].gm / (r2 * sqrt(r2));
      for (int i = 0; i < 3; ++i) {
        pull[i] += f * d[i];
      }
    }
    for (int i = 0; i < 3; ++i) {
      a[3 * b + (size_t)i] += pull[i];
    }
  }
  return MYB_OK;
}

/* ---- One evaluation: the equations of rotation ----------------------------------------------- */

static myb_status rotate(const myb_system* s, const rotation_state* r, double* a) {
  const myb_system_rotation* model = &r->model;
  const double* wm = r->omega;
  const double* wc = r->core;
  const double ac = model->core[0];
  const double cc = model->core[1];
  const int core = cc > 0.0;
  double coupling[3] = {0.0, 0.0, 0.0}; /* of the core on the mantle */
  if (core) {
    const double zxw[3] = {-wc[1], wc[0], 0.0};
    for (int i = 0; i < 3; ++i) {
      coupling[i] = model->core_friction * (wc[i] - wm[i]) + (cc - ac) * wc[2] * zxw[i];
    }
  }
  double l[3];
  double gyro[3];
  double flow[3];
  myb_internal_mat3_vec(r->mantle, wm, l);
  myb_vec3_cross(wm, l, gyro);
  myb_internal_mat3_vec(r->mantle_rate, wm, flow);
  double rate[3];
  for (int i = 0; i < 3; ++i) {
    rate[i] = r->torque[i] + coupling[i] - flow[i] - gyro[i];
  }
  /* the rate of the spin distortion is linear in the rate looked for:
   * d/dt [fs (w w^T - |w|^2 1 / 3)] w = fs (|w|^2 w' + w (w . w') / 3) */
  double matrix[9];
  memcpy(matrix, r->mantle, sizeof matrix);
  if (model->k2 != 0.0) {
    const double fs = distortion_scale(r) / (3.0 * s->gm[r->body]);
    const double w2 = myb_vec3_dot(wm, wm);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        matrix[3 * i + j] += fs * ((i == j ? w2 : 0.0) + wm[i] * wm[j] / 3.0);
      }
    }
  }
  if (myb_lu_solve(3, matrix, rate, 1) != MYB_OK) {
    return MYB_E_DOMAIN;
  }
  /* q'' = (q' (0, w) + q (0, w')) / 2, for the quaternion of the state (of norm r->norm) */
  const double w4[4] = {0.0, wm[0], wm[1], wm[2]};
  const double dw4[4] = {0.0, rate[0], rate[1], rate[2]};
  double first[4];
  double second[4];
  myb_quat_multiply(r->dq, w4, first);
  myb_quat_multiply(r->q, dw4, second);
  for (int i = 0; i < 4; ++i) {
    a[i] = 0.5 * r->norm * (first[i] + second[i]);
  }
  if (core) {
    const double lc[3] = {ac * wc[0], ac * wc[1], cc * wc[2]};
    double gc[3];
    myb_vec3_cross(wm, lc, gc);
    a[CORE_OFFSET] = (-gc[0] - coupling[0]) / ac;
    a[CORE_OFFSET + 1] = (-gc[1] - coupling[1]) / ac;
    a[CORE_OFFSET + 2] = (-gc[2] - coupling[2]) / cc;
  }
  return MYB_OK;
}

/* The point masses, the frames and the inertia tensors: what every level needs. */
static myb_status prepare(myb_system* s, double t, const double* x, const double* v,
                          double* point) {
  for (size_t b = 0; b < s->bodies; ++b) {
    const body_model* m = &s->models[b];
    if (m->gm_rate != 0.0 && m->present) {
      const double gm = m->gm_base * (1.0 + m->gm_rate * (t - m->gm_epoch));
      if (!(gm >= 0.0) || !isfinite(gm)) {
        return MYB_E_DOMAIN; /* the mass has run out */
      }
      s->gm[b] = gm;
      myb_internal_nbody_set_gm(s->nbody, b, gm);
    }
  }
  myb_status st = myb_nbody_acceleration_xv(s->nbody, t, x, v, point);
  if (st != MYB_OK) {
    return st;
  }
  if (point != s->point) {
    memcpy(s->point, point, 3 * s->bodies * sizeof *point);
  }
  st = frames(s, t, x, v);
  for (size_t k = 0; k < s->rotation_count && st == MYB_OK; ++k) {
    st = distort(s, &s->rotations[k], x, v);
  }
  return st;
}

myb_status myb_system_acceleration(void* system, double t, const double* x, const double* v,
                                   double* a) {
  MYB_CHECK(system != NULL && x != NULL && v != NULL && a != NULL);
  myb_system* s = system;
  const size_t translational = 3 * s->bodies;
  myb_status st = prepare(s, t, x, v, a);
  if (st != MYB_OK) {
    return st;
  }
  memset(a + translational, 0, (myb_system_coordinates(s) - translational) * sizeof *a);
  st = apply_fields(s, x, a);
  if (st == MYB_OK) {
    st = apply_figures(s, x, a);
  }
  if (st == MYB_OK) {
    st = apply_tides(s, x, v, a);
  }
  if (st == MYB_OK) {
    st = apply_pairs(s, x, v, a);
  }
  if (st == MYB_OK) {
    st = apply_actions(s, t, x, v, a);
  }
  if (st == MYB_OK) {
    st = apply_drags(s, x, v, a);
  }
  if (st == MYB_OK) {
    st = apply_guided(s, t, x, a);
  }
  for (size_t k = 0; k < s->rotation_count && st == MYB_OK; ++k) {
    st = rotate(s, &s->rotations[k], a + rotation_offset(s, (int)k));
  }
  return st;
}

myb_status myb_system_inertia(myb_system* system, size_t body, double t, const double* y,
                              double inertia[9]) {
  MYB_CHECK(system != NULL && y != NULL && inertia != NULL);
  if (body >= system->bodies || system->models[body].rotation == NO_ROTATION) {
    return MYB_E_DOMAIN;
  }
  const size_t n = myb_system_coordinates(system);
  const myb_status st = prepare(system, t, y, y + n, system->point);
  if (st != MYB_OK) {
    return st;
  }
  memcpy(inertia, system->rotations[system->models[body].rotation].total, 9 * sizeof *inertia);
  return MYB_OK;
}

/* ---- Encounters ------------------------------------------------------------------------------ */

static myb_status need_scratch(myb_system* s) {
  size_t size;
  if (ckd_mul(&size, myb_system_coordinates(s), (size_t)4)) {
    return MYB_E_OVERFLOW;
  }
  if (s->scratch != NULL && s->scratch_size >= size) {
    return MYB_OK;
  }
  void* block = NULL;
  const myb_status st = myb_internal_alloc_array(&s->allocator, size, sizeof(double), &block);
  if (st != MYB_OK) {
    return st;
  }
  myb_internal_release(&s->allocator, s->scratch);
  s->scratch = block;
  s->scratch_size = size;
  return MYB_OK;
}

/* The distance of two bodies and its rate (the radial velocity) in the state y. */
static double separation(const myb_system* s, const double* y, size_t a, size_t b,
                         double* rate) {
  const size_t n = myb_system_coordinates(s);
  double d[3];
  double dv[3];
  relative(y, a, b, d);
  relative(y + n, a, b, dv);
  const double r = myb_vec3_norm(d);
  *rate = r > 0.0 ? myb_vec3_dot(d, dv) / r : 0.0;
  return r;
}

typedef enum search_kind { SEARCH_CONTACT = 0, SEARCH_CLOSEST = 1 } search_kind;

/* The function whose root is looked for, at the time t of the last step. */
static myb_status search_value(myb_system* s, myb_ias15* in, search_kind kind, size_t a,
                               size_t b, double t, double* value) {
  double* work = s->scratch + 2 * myb_system_coordinates(s);
  const myb_status st = myb_ias15_dense(in, t, work);
  if (st != MYB_OK) {
    return st;
  }
  double rate = 0.0;
  const double r = separation(s, work, a, b, &rate);
  *value = kind == SEARCH_CONTACT ? r - (s->models[a].radius + s->models[b].radius) : rate;
  return MYB_OK;
}

/* The root in [lo, hi], where the function is negative at one end and not at the other:
 * bisection to neighbouring doubles; the end returned is the one where the function has the
 * sign `after` takes — past the event. */
static myb_status bisect(myb_system* s, myb_ias15* in, search_kind kind, size_t a, size_t b,
                         double lo, double hi, double* root) {
  double f_lo = 0.0;
  myb_status st = search_value(s, in, kind, a, b, lo, &f_lo);
  /* bound: ROOT_STEPS halvings of an interval of doubles */
  for (int k = 0; k < ROOT_STEPS && st == MYB_OK; ++k) {
    const double mid = lo + 0.5 * (hi - lo);
    if (mid == lo || mid == hi) {
      break;
    }
    double f_mid = 0.0;
    st = search_value(s, in, kind, a, b, mid, &f_mid);
    if ((f_mid < 0.0) == (f_lo < 0.0)) {
      lo = mid;
      f_lo = f_mid;
    } else {
      hi = mid;
    }
  }
  *root = hi;
  return st;
}

myb_status myb_system_advance(myb_system* system, myb_ias15* integrator, double t_end,
                              double watch, myb_system_event* event) {
  MYB_CHECK(system != NULL && integrator != NULL && event != NULL);
  myb_system* s = system;
  if (!(watch >= 0.0) || !isfinite(watch) || !isfinite(t_end)) {
    return MYB_E_DOMAIN;
  }
  myb_status st = need_scratch(s);
  if (st != MYB_OK) {
    return st;
  }
  double* start = s->scratch;
  double* work = s->scratch + 2 * myb_system_coordinates(s);
  memset(event, 0, sizeof *event);
  /* bound: the steps of the integrator, which has its own budget per step and ends at t_end */
  for (int64_t step = 0; step < INT64_MAX; ++step) {
    const double t0 = myb_ias15_time(integrator);
    if (t0 == t_end) {
      event->kind = MYB_SYSTEM_ARRIVED;
      event->t = t0;
      return MYB_OK;
    }
    st = myb_ias15_step(integrator, t_end);
    if (st != MYB_OK) {
      return st;
    }
    const double t1 = myb_ias15_time(integrator);
    const double* y1 = myb_ias15_state(integrator);
    st = myb_ias15_dense(integrator, t0, start);
    if (st != MYB_OK) {
      return st;
    }
    /* the earliest event of the step, pair by pair */
    int found = 0;
    double best = t1;
    myb_system_event first = {0};
    for (size_t a = 0; a < s->bodies; ++a) {
      for (size_t b = a + 1; b < s->bodies; ++b) {
        const body_model* ma = &s->models[a];
        const body_model* mb = &s->models[b];
        const int spheres = ma->radius > 0.0 && mb->radius > 0.0;
        if (!ma->present || !mb->present || (!spheres && watch == 0.0)) {
          continue;
        }
        double rate0 = 0.0;
        double rate1 = 0.0;
        const double r0 = separation(s, start, a, b, &rate0);
        const double r1 = separation(s, y1, a, b, &rate1);
        const double sum = ma->radius + mb->radius;
        /* the closest approach inside the step, if the bodies turned from nearing to parting */
        double t_closest = t1;
        int turned = 0;
        const int again = s->reported && s->reported_a == a && s->reported_b == b &&
                          s->reported_t == t0;
        if (rate0 < 0.0 && rate1 >= 0.0 && !again) {
          st = bisect(s, integrator, SEARCH_CLOSEST, a, b, t0, t1, &t_closest);
          if (st != MYB_OK) {
            return st;
          }
          turned = 1;
        }
        double t_event = t1;
        myb_system_event_kind kind = MYB_SYSTEM_ARRIVED;
        if (spheres && r0 > sum) {
          /* contact: by the end of the step, or before the closest approach inside it */
          double r_closest = r1;
          if (turned) {
            double unused = 0.0;
            st = myb_ias15_dense(integrator, t_closest, work);
            if (st != MYB_OK) {
              return st;
            }
            r_closest = separation(s, work, a, b, &unused);
          }
          if (r_closest <= sum || r1 <= sum) {
            st = bisect(s, integrator, SEARCH_CONTACT, a, b, t0,
                        turned && r_closest <= sum ? t_closest : t1, &t_event);
            if (st != MYB_OK) {
              return st;
            }
            kind = MYB_SYSTEM_CONTACT;
          }
        }
        if (kind == MYB_SYSTEM_ARRIVED && turned && watch > 0.0) {
          double unused = 0.0;
          st = myb_ias15_dense(integrator, t_closest, work);
          if (st != MYB_OK) {
            return st;
          }
          if (separation(s, work, a, b, &unused) < watch) {
            kind = MYB_SYSTEM_APPROACH;
            t_event = t_closest;
          }
        }
        /* the earliest in the direction of the integration */
        if (kind != MYB_SYSTEM_ARRIVED &&
            (!found || fabs(t_event - t0) < fabs(best - t0))) {
          found = 1;
          best = t_event;
          first.kind = kind;
          first.a = a;
          first.b = b;
        }
      }
    }
    if (!found) {
      continue;
    }
    if (best != t1) {
      st = myb_internal_ias15_retake(integrator, best);
      if (st != MYB_OK) {
        return st;
      }
    }
    double rate = 0.0;
    first.t = myb_ias15_time(integrator);
    first.distance = separation(s, myb_ias15_state(integrator), first.a, first.b, &rate);
    *event = first;
    s->reported = first.kind == MYB_SYSTEM_APPROACH;
    s->reported_a = first.a;
    s->reported_b = first.b;
    s->reported_t = first.t;
    return MYB_OK;
  }
  return MYB_E_NOCONV;
}

/* Drops the levels that name `body`. */
static size_t drop_pairs(pair* list, size_t count, size_t body) {
  size_t kept = 0;
  for (size_t k = 0; k < count; ++k) {
    if (list[k].a != body && list[k].b != body) {
      list[kept++] = list[k];
    }
  }
  return kept;
}

static int in_pn_pair(const myb_system* s, size_t body) {
  for (size_t k = 0; k < s->pn_count; ++k) {
    if (s->pn[k].a == body || s->pn[k].b == body) {
      return 1;
    }
  }
  return 0;
}

/* Copies the state of the integrator and starts it again from the copy. */
static myb_status restart(myb_system* s, myb_ias15* in) {
  return myb_ias15_start(in, myb_system_acceleration, s, myb_ias15_time(in), s->scratch);
}

myb_status myb_system_merge(myb_system* system, myb_ias15* integrator, size_t a, size_t b,
                            myb_system_merger* merger) {
  MYB_CHECK(system != NULL && integrator != NULL);
  myb_system* s = system;
  if (!two_bodies(s, a, b) || s->models[a].rotation != NO_ROTATION ||
      s->models[b].rotation != NO_ROTATION || in_pn_pair(s, a) || in_pn_pair(s, b)) {
    return MYB_E_DOMAIN;
  }
  myb_status st = need_scratch(s);
  if (st != MYB_OK) {
    return st;
  }
  const size_t n = myb_system_coordinates(s);
  memcpy(s->scratch, myb_ias15_state(integrator), 2 * n * sizeof(double));
  double* x = s->scratch;
  double* v = s->scratch + n;
  const double ga = s->gm[a];
  const double gb = s->gm[b];
  const double total = ga + gb;
  double d[3];
  double dv[3];
  relative(x, a, b, d);
  relative(v, a, b, dv);
  const double r = myb_vec3_norm(d);
  if (!(r > 0.0) || !isfinite(r)) {
    return MYB_E_DOMAIN;
  }
  if (merger != NULL) {
    const double reduced = total > 0.0 ? ga * gb / total : 0.0;
    double l[3];
    myb_vec3_cross(d, dv, l);
    merger->energy = 0.5 * reduced * myb_vec3_dot(dv, dv) - ga * gb / r;
    for (int i = 0; i < 3; ++i) {
      merger->angular_momentum[i] = reduced * l[i];
    }
  }
  /* the centre of mass takes the place of a; two bodies without mass meet half-way */
  const double share = total > 0.0 ? gb / total : 0.5;
  for (int i = 0; i < 3; ++i) {
    x[3 * a + (size_t)i] += share * d[i];
    v[3 * a + (size_t)i] += share * dv[i];
  }
  const double ra = s->models[a].radius;
  const double rb = s->models[b].radius;
  s->models[a].radius = myb_cbrt(ra * ra * ra + rb * rb * rb);
  s->gm[a] = total;
  s->gm[b] = 0.0;
  /* a mass that changes goes on changing at its relative rate from what it is now */
  s->models[b].gm_rate = 0.0;
  if (s->models[a].gm_rate != 0.0) {
    s->models[a].gm_base = total;
    s->models[a].gm_epoch = myb_ias15_time(integrator);
  }
  myb_internal_nbody_set_gm(s->nbody, a, total);
  myb_internal_nbody_set_gm(s->nbody, b, 0.0);
  s->models[b].present = 0;
  s->models[b].radius = 0.0;
  for (size_t k = 0; k < s->bodies && st == MYB_OK; ++k) {
    if (k == b) {
      continue;
    }
    st = myb_nbody_set_interaction(s->nbody, k, b, MYB_INTERACTION_NONE);
    if (st == MYB_OK) {
      st = myb_nbody_set_interaction(s->nbody, b, k, MYB_INTERACTION_NONE);
    }
  }
  if (st != MYB_OK) {
    return st;
  }
  s->field_count = drop_pairs(s->fields, s->field_count, b);
  s->figure_count = drop_pairs(s->figures, s->figure_count, b);
  s->distortion_count = drop_pairs(s->distortions, s->distortion_count, b);
  size_t kept = 0;
  for (size_t k = 0; k < s->tide_count; ++k) {
    const triple* e = &s->tides[k];
    if (e->body != b && e->raiser != b && e->point != b) {
      s->tides[kept++] = *e;
    }
  }
  s->tide_count = kept;
  kept = 0;
  for (size_t k = 0; k < s->drag_count; ++k) {
    if (s->drags[k].body != b && s->drags[k].point != b) {
      s->drags[kept++] = s->drags[k];
    }
  }
  s->drag_count = kept;
  s->reported = 0;
  return restart(s, integrator);
}

myb_status myb_system_bounce(myb_system* system, myb_ias15* integrator, size_t a, size_t b,
                             double restitution) {
  MYB_CHECK(system != NULL && integrator != NULL);
  myb_system* s = system;
  if (!two_bodies(s, a, b) || !(restitution >= 0.0) || !(restitution <= 1.0)) {
    return MYB_E_DOMAIN;
  }
  myb_status st = need_scratch(s);
  if (st != MYB_OK) {
    return st;
  }
  const size_t n = myb_system_coordinates(s);
  memcpy(s->scratch, myb_ias15_state(integrator), 2 * n * sizeof(double));
  double* x = s->scratch;
  double* v = s->scratch + n;
  double d[3];
  double dv[3];
  relative(x, a, b, d);
  relative(v, a, b, dv);
  const double r = myb_vec3_norm(d);
  if (!(r > 0.0) || !isfinite(r)) {
    return MYB_E_DOMAIN;
  }
  const double along = myb_vec3_dot(d, dv) / r;
  if (along >= 0.0) {
    return MYB_OK; /* parting already */
  }
  const double total = s->gm[a] + s->gm[b];
  const double share_a = total > 0.0 ? s->gm[b] / total : 0.5;
  const double share_b = total > 0.0 ? s->gm[a] / total : 0.5;
  const double change = (1.0 + restitution) * along;
  for (int i = 0; i < 3; ++i) {
    const double unit = d[i] / r;
    v[3 * a + (size_t)i] += share_a * change * unit;
    v[3 * b + (size_t)i] -= share_b * change * unit;
  }
  s->reported = 0;
  return restart(s, integrator);
}
