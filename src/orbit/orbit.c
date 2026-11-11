/* SPDX-License-Identifier: BSD-3-Clause */
/* Propagation of one object under a force model: IAS15 with its own steps, every state from
 * the dense output of the step that contains its time.
 *
 * Breaks of the model (myb_force_model_add_break): after a step the break functions are
 * examined on its dense output as events are. A step that contains a sign change is
 * withdrawn and taken again from its beginning up to the change
 * (myb_internal_ias15_retake), so that every step kept lies on one side of every break. The
 * integrator is not started again: the state, the transition matrix and the sensitivities
 * simply continue, and the step after the break starts from the polynomial of the step
 * withdrawn.
 */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "integ/event.h"
#include "integ/ias15_internal.h"
#include "myb/orbit.h"
#include "myb/uncert.h"
#include "uncert/variational_internal.h"

enum { EVENT_PARTS = 4, EVENT_NUDGES = 3 };

static const double event_parts = 4.0; /* EVENT_PARTS */

static const double default_epsilon = 1e-10;
static const int64_t default_max_steps = 1000000000;

struct myb_orbit {
  myb_allocator allocator;
  myb_force_model* model;
  myb_ias15* ias;       /* without the variational equations */
  myb_variational* var; /* with them */
  myb_dop853* dop;      /* the other integrator */
  int approaching;      /* DOP853: a step that held a break is being redone up to it */
  double break_time;
  int64_t max_steps;
  int started;
  double t_end;
  double direction;  /* +1 or -1 */
  double reached;    /* always inside the last step (or t0 before any) */
  int have_step;
  double step_begin; /* the last step */
  double step_end;
  double y0[6];
  int64_t steps;
  double scratch[6];
  /* breaks */
  int64_t breaks;
  signed char break_sign[MYB_MODEL_MAX_BREAKS]; /* at the time integrated to; 0 = not known */
  double y_begin[6]; /* at the beginning of the step being taken */
  /* DOP853 with the variational equations: the state and its variations as one system of
   * the first order, `wide` values — the state (6), then for every column of the transition
   * matrix and of the sensitivities its part of the position (3) and of the velocity (3) */
  size_t wide;       /* 0 without them */
  size_t parameters; /* of the model, when the object was made */
  double* wide_start;   /* 3 * wide: the state the integration (re)started from, */
  double* wide_begin;   /* that at the beginning of the step being taken, */
  double* wide_scratch; /* and one to read into */
};

void myb_orbit_destroy(myb_orbit* orbit) {
  if (orbit == NULL) {
    return;
  }
  const myb_allocator a = orbit->allocator;
  myb_ias15_destroy(orbit->ias);
  myb_variational_destroy(orbit->var);
  myb_dop853_destroy(orbit->dop);
  myb_internal_release(&a, orbit->wide_start);
  myb_internal_release(&a, orbit);
}

myb_status myb_orbit_create(myb_force_model* model, const myb_orbit_options* options,
                            const myb_allocator* allocator, myb_orbit** orbit) {
  MYB_CHECK(model != NULL && orbit != NULL);
  const myb_orbit_options none = {0};
  const myb_orbit_options* opt = options != NULL ? options : &none;
  if (!(opt->epsilon >= 0.0) || !(opt->h_max >= 0.0) || opt->max_steps < 0 ||
      !isfinite(opt->epsilon) || !isfinite(opt->h_max) || !(opt->rtol >= 0.0) ||
      !(opt->atol >= 0.0) || !isfinite(opt->rtol) || !isfinite(opt->atol) ||
      (opt->integrator != MYB_ORBIT_IAS15 && opt->integrator != MYB_ORBIT_DOP853)) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_orbit), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_orbit* o = block;
  memset(o, 0, sizeof *o);
  o->allocator = a;
  o->model = model;
  o->max_steps = opt->max_steps > 0 ? opt->max_steps : default_max_steps;
  const double epsilon = opt->epsilon > 0.0 ? opt->epsilon : default_epsilon;
  if (opt->integrator == MYB_ORBIT_DOP853) {
    const myb_dop853_options options853 = {.rtol = opt->rtol > 0.0 ? opt->rtol : 1e-12,
                                           .atol = opt->atol > 0.0 ? opt->atol : 1e-9,
                                           .h_max = opt->h_max};
    size_t n = 6;
    if (opt->variational) {
      o->parameters = myb_force_model_parameter_count(model);
      o->wide = 6 + 6 * (6 + o->parameters);
      n = o->wide;
      status = myb_internal_alloc_array(&a, 3 * o->wide, sizeof(double), &block);
      if (status == MYB_OK) {
        o->wide_start = block;
        o->wide_begin = o->wide_start + o->wide;
        o->wide_scratch = o->wide_begin + o->wide;
        memset(o->wide_start, 0, 3 * o->wide * sizeof(double));
      }
    }
    if (status == MYB_OK) {
      status = myb_dop853_create(n, &options853, allocator, &o->dop);
    }
  } else if (opt->variational) {
    const myb_variational_options vo = {.epsilon = epsilon, .h_max = opt->h_max};
    status = myb_variational_create(model, &vo, allocator, &o->var);
  } else {
    const myb_ias15_options io = {.epsilon = epsilon, .h_max = opt->h_max};
    status = myb_ias15_create(3, &io, allocator, &o->ias);
  }
  if (status != MYB_OK) {
    myb_orbit_destroy(o);
    return status;
  }
  *orbit = o;
  return MYB_OK;
}

/* The equations of motion as a system of the first order, for DOP853: context = the model. */
static myb_status first_order(void* model, double t, const double* y, double* dydt) {
  dydt[0] = y[3];
  dydt[1] = y[4];
  dydt[2] = y[5];
  return myb_force_model_acceleration(model, t, y, y + 3, dydt + 3);
}

/* The same with the variations (context = the orbit): for every column d of the transition
 * matrix and of the sensitivities, d(position)' = d(velocity) and
 * d(velocity)' = (da/dx) d(position) + (da/dv) d(velocity) (+ da/dp for a parameter). */
static myb_status first_order_wide(void* orbit, double t, const double* y, double* dydt) {
  const myb_orbit* o = orbit;
  dydt[0] = y[3];
  dydt[1] = y[4];
  dydt[2] = y[5];
  myb_status status = myb_force_model_acceleration(o->model, t, y, y + 3, dydt + 3);
  if (status != MYB_OK) {
    return status;
  }
  double jx[9];
  double jv[9];
  double jp[3 * MYB_MODEL_MAX_PARAMETERS];
  status = myb_force_model_jacobian(o->model, t, y, y + 3, jx, jv, o->parameters > 0 ? jp : NULL);
  if (status != MYB_OK) {
    return status;
  }
  const size_t columns = 6 + o->parameters;
  for (size_t col = 0; col < columns; ++col) {
    const double* d = y + 6 + 6 * col;
    double* dd = dydt + 6 + 6 * col;
    for (int i = 0; i < 3; ++i) {
      dd[i] = d[3 + i];
      double s = jx[3 * i] * d[0] + jx[3 * i + 1] * d[1] + jx[3 * i + 2] * d[2] +
                 jv[3 * i] * d[3] + jv[3 * i + 1] * d[4] + jv[3 * i + 2] * d[5];
      if (col >= 6) {
        s += jp[(size_t)i * o->parameters + (col - 6)];
      }
      dd[3 + i] = s;
    }
  }
  return MYB_OK;
}

/* The state, the transition matrix and the sensitivities out of the wide vector. */
static void unpack_wide(const myb_orbit* o, const double* wide, double y[6], double phi[36],
                        double* sensitivity) {
  if (y != NULL) {
    memcpy(y, wide, 6 * sizeof *wide);
  }
  const size_t columns = 6 + o->parameters;
  for (size_t col = 0; col < columns; ++col) {
    const double* d = wide + 6 + 6 * col;
    for (size_t i = 0; i < 6; ++i) {
      if (col < 6 && phi != NULL) {
        phi[i * 6 + col] = d[i];
      } else if (col >= 6 && sensitivity != NULL) {
        sensitivity[i * o->parameters + (col - 6)] = d[i];
      }
    }
  }
}

/* (Re)starts DOP853 from o->wide_start, or from the six numbers y without the variations. */
static myb_status start_dop853(myb_orbit* o, double t, const double* y) {
  if (o->wide == 0) {
    return myb_dop853_start(o->dop, first_order, o->model, t, y);
  }
  return myb_dop853_start(o->dop, first_order_wide, o, t, o->wide_start);
}

static signed char sign_of(double g) {
  return g > 0.0 ? 1 : (g < 0.0 ? -1 : 0);
}

myb_status myb_orbit_start(myb_orbit* orbit, double t0, const double y0[6], double t_end) {
  MYB_CHECK(orbit != NULL && y0 != NULL);
  myb_orbit* o = orbit;
  o->started = 0;
  if (!isfinite(t0) || !isfinite(t_end) || t0 == t_end) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 6; ++i) {
    if (!isfinite(y0[i])) {
      return MYB_E_DOMAIN;
    }
  }
  if (o->wide != 0) { /* the state, the unit matrix, no sensitivity yet */
    memset(o->wide_start, 0, o->wide * sizeof(double));
    memcpy(o->wide_start, y0, 6 * sizeof(double));
    for (size_t col = 0; col < 6; ++col) {
      o->wide_start[6 + 6 * col + col] = 1.0;
    }
  }
  const myb_status status =
      o->dop != NULL   ? start_dop853(o, t0, y0)
      : o->var != NULL ? myb_variational_start(o->var, t0, y0)
                       : myb_ias15_start(o->ias, myb_force_model_accel, o->model, t0, y0);
  if (status != MYB_OK) {
    return status;
  }
  memcpy(o->y0, y0, sizeof o->y0);
  o->t_end = t_end;
  o->direction = t_end > t0 ? 1.0 : -1.0;
  o->reached = t0;
  o->have_step = 0;
  o->step_begin = t0;
  o->step_end = t0;
  o->steps = 0;
  o->breaks = 0;
  o->approaching = 0;
  const size_t count = myb_force_model_break_count(o->model);
  for (size_t k = 0; k < count; ++k) {
    myb_break_function g = NULL;
    void* context = NULL;
    myb_force_model_break(o->model, k, &g, &context);
    o->break_sign[k] = sign_of(g(context, t0, y0));
  }
  o->started = 1;
  return MYB_OK;
}

double myb_orbit_time(const myb_orbit* orbit) {
  MYB_CHECK(orbit != NULL && orbit->started);
  return orbit->reached;
}

int64_t myb_orbit_steps(const myb_orbit* orbit) {
  MYB_CHECK(orbit != NULL);
  return orbit->steps;
}

int64_t myb_orbit_breaks(const myb_orbit* orbit) {
  MYB_CHECK(orbit != NULL);
  return orbit->breaks;
}

static double integrator_time(const myb_orbit* o) {
  if (o->dop != NULL) {
    return myb_dop853_time(o->dop);
  }
  return o->var != NULL ? myb_variational_time(o->var) : myb_ias15_time(o->ias);
}

/* The state at t inside the last step (t0 before any step). */
static myb_status dense(myb_orbit* o, double t, double y[6]) {
  if (!o->have_step) {
    memcpy(y, o->y0, sizeof o->y0);
    return MYB_OK;
  }
  if (o->dop != NULL && o->wide != 0) {
    const myb_status status = myb_dop853_dense(o->dop, t, o->wide_scratch);
    if (status == MYB_OK) {
      memcpy(y, o->wide_scratch, 6 * sizeof *y);
    }
    return status;
  }
  if (o->dop != NULL) {
    return myb_dop853_dense(o->dop, t, y);
  }
  return o->var != NULL ? myb_variational_dense(o->var, t, y, NULL, NULL)
                        : myb_ias15_dense(o->ias, t, y);
}

static myb_status dense_view(void* orbit, double t, double* y) {
  return dense(orbit, t, y);
}

/* The first sign change of a break of the model in the step [ta, tb] just taken, ya the state
 * at ta: examined, as events are, at the ends of four parts of the step and located on its
 * dense output. The signs recorded are brought up to what was examined. */
static myb_status first_break(myb_orbit* o, size_t count, double ta, const double ya[6],
                              double tb, int* found, double* t_break) {
  *found = 0;
  double t0 = ta;
  double y0[6];
  memcpy(y0, ya, sizeof y0);
  for (int part = 1; part <= EVENT_PARTS; ++part) {
    const double t1 = part < EVENT_PARTS ? ta + (tb - ta) * ((double)part / event_parts) : tb;
    double y1[6];
    myb_status status = dense(o, t1, y1);
    if (status != MYB_OK) {
      return status;
    }
    size_t which = 0;
    for (size_t k = 0; k < count; ++k) {
      myb_break_function g = NULL;
      void* context = NULL;
      myb_force_model_break(o->model, k, &g, &context);
      const signed char s0 = o->break_sign[k];
      const signed char s1 = sign_of(g(context, t1, y1));
      if (s0 == 0 || s1 == s0) {
        o->break_sign[k] = s1;
        continue;
      }
      double candidate = t1; /* a zero exactly at the end of the part */
      if (s1 != 0) {
        const myb_internal_step_view view = {.dense = dense_view,
                                             .integrator = o,
                                             .m = 6,
                                             .t0 = t0,
                                             .y0 = y0,
                                             .t1 = t1,
                                             .y1 = y1,
                                             .scratch = o->scratch};
        int located = 0;
        status = myb_internal_locate_event(&view, g, context, &located, &candidate, NULL);
        if (status != MYB_OK) {
          return status;
        }
        if (!located) { /* the sign at t0 is no longer the one recorded */
          o->break_sign[k] = s1;
          continue;
        }
      }
      if (!*found || (candidate - *t_break) * o->direction < 0.0) {
        *found = 1;
        *t_break = candidate;
        which = k;
      }
    }
    if (*found) {
      o->break_sign[which] = (signed char)-o->break_sign[which];
      return MYB_OK;
    }
    t0 = t1;
    memcpy(y0, y1, sizeof y0);
  }
  return MYB_OK;
}

/* The same for DOP853. A method of one step has nothing to carry over: the step that held a
 * break is taken again from its beginning — in as many steps as the integrator likes — up to
 * the break, and the next one starts there. */
static myb_status advance_dop853(myb_orbit* o) {
  const size_t count = o->approaching ? 0 : myb_force_model_break_count(o->model);
  const double before = myb_dop853_time(o->dop);
  if (count > 0) {
    memcpy(o->y_begin, myb_dop853_state(o->dop), sizeof o->y_begin);
    if (o->wide != 0) {
      memcpy(o->wide_begin, myb_dop853_state(o->dop), o->wide * sizeof(double));
    }
  }
  myb_status status = myb_dop853_step(o->dop, o->approaching ? o->break_time : o->t_end);
  if (status != MYB_OK) {
    return status;
  }
  double after = myb_dop853_time(o->dop);
  if (after == before) {
    return MYB_OK;
  }
  o->step_begin = before;
  o->step_end = after;
  o->have_step = 1;
  o->steps++;
  if (o->approaching) {
    o->approaching = after != o->break_time;
    return MYB_OK;
  }
  if (count == 0) {
    return MYB_OK;
  }
  int found = 0;
  double t_break = 0.0;
  status = first_break(o, count, before, o->y_begin, after, &found, &t_break);
  if (status != MYB_OK || !found) {
    return status;
  }
  o->breaks++;
  if (t_break == before || t_break == after) {
    return MYB_OK;
  }
  o->have_step = 0;
  if (o->wide != 0) {
    memcpy(o->wide_start, o->wide_begin, o->wide * sizeof(double));
  }
  status = start_dop853(o, before, o->y_begin);
  if (status == MYB_OK) {
    status = myb_dop853_step(o->dop, t_break);
  }
  if (status != MYB_OK) {
    return status;
  }
  after = myb_dop853_time(o->dop);
  o->step_end = after;
  o->have_step = 1;
  o->approaching = after != t_break;
  o->break_time = t_break;
  return MYB_OK;
}

/* One step of the integrator towards t_end, ended at the first break of the model in it. */
static myb_status advance(myb_orbit* o) {
  if (o->dop != NULL) {
    return advance_dop853(o);
  }
  const size_t count = myb_force_model_break_count(o->model);
  const double before = integrator_time(o);
  if (count > 0) {
    if (o->var != NULL) {
      myb_variational_get(o->var, o->y_begin, NULL, NULL);
    } else {
      memcpy(o->y_begin, myb_ias15_state(o->ias), sizeof o->y_begin);
    }
  }
  myb_status status = o->var != NULL ? myb_variational_step(o->var, o->t_end)
                                     : myb_ias15_step(o->ias, o->t_end);
  if (status != MYB_OK) {
    return status;
  }
  const double after = integrator_time(o);
  if (after == before) {
    return MYB_OK;
  }
  o->step_begin = before;
  o->step_end = after;
  o->have_step = 1;
  o->steps++;
  if (count == 0) {
    return MYB_OK;
  }
  int found = 0;
  double t_break = 0.0;
  status = first_break(o, count, before, o->y_begin, after, &found, &t_break);
  if (status != MYB_OK || !found) {
    return status;
  }
  o->breaks++;
  /* a change at an end of the step leaves all of it on one side */
  if (t_break == before || t_break == after) {
    return MYB_OK;
  }
  status = o->var != NULL ? myb_internal_variational_retake(o->var, t_break)
                          : myb_internal_ias15_retake(o->ias, t_break);
  if (status == MYB_E_DOMAIN) {
    return MYB_OK; /* closer to the beginning than the time resolves: the step stays */
  }
  if (status != MYB_OK) {
    o->have_step = 0;
    return status;
  }
  o->step_end = integrator_time(o);
  return MYB_OK;
}

/* Steps until the integrator has reached t (already validated to lie before t_end). */
static myb_status reach(myb_orbit* o, double t) {
  int64_t budget = o->max_steps;
  /* bound: the step budget */
  while ((t - integrator_time(o)) * o->direction > 0.0) {
    if (budget == 0) {
      return MYB_E_NOCONV;
    }
    --budget;
    const myb_status status = advance(o);
    if (status != MYB_OK) {
      return status;
    }
  }
  return MYB_OK;
}

/* MYB_E_DOMAIN unless t lies between the time reached and t_end. */
static myb_status in_order(const myb_orbit* o, double t) {
  if (!isfinite(t) || (t - o->reached) * o->direction < 0.0 ||
      (o->t_end - t) * o->direction < 0.0) {
    return MYB_E_DOMAIN;
  }
  return MYB_OK;
}

myb_status myb_orbit_state(myb_orbit* orbit, double t, double y[6]) {
  MYB_CHECK(orbit != NULL && orbit->started && y != NULL);
  myb_status status = in_order(orbit, t);
  if (status != MYB_OK) {
    return status;
  }
  status = reach(orbit, t);
  if (status != MYB_OK) {
    return status;
  }
  status = dense(orbit, t, y);
  if (status == MYB_OK) {
    orbit->reached = t;
  }
  return status;
}

myb_status myb_orbit_states(myb_orbit* orbit, const double* times, size_t count,
                            double* states) {
  MYB_CHECK(orbit != NULL && (count == 0 || (times != NULL && states != NULL)));
  for (size_t k = 0; k < count; ++k) {
    const myb_status status = myb_orbit_state(orbit, times[k], states + 6 * k);
    if (status != MYB_OK) {
      return status;
    }
  }
  return MYB_OK;
}

myb_status myb_orbit_partials(myb_orbit* orbit, double t, double y[6], double phi[36],
                              double* sensitivity) {
  MYB_CHECK(orbit != NULL && orbit->started);
  myb_orbit* o = orbit;
  if (o->var == NULL && o->wide == 0) {
    return MYB_E_DOMAIN;
  }
  myb_status status = in_order(o, t);
  if (status != MYB_OK) {
    return status;
  }
  status = reach(o, t);
  if (status != MYB_OK) {
    return status;
  }
  if (o->wide != 0) {
    const double* wide = o->wide_start;
    if (o->have_step) {
      status = myb_dop853_dense(o->dop, t, o->wide_scratch);
      wide = o->wide_scratch;
    }
    if (status == MYB_OK) {
      unpack_wide(o, wide, y, phi, sensitivity);
    }
  } else if (o->have_step) {
    status = myb_variational_dense(o->var, t, y, phi, sensitivity);
  } else {
    myb_variational_get(o->var, y, phi, sensitivity);
  }
  if (status == MYB_OK) {
    o->reached = t;
  }
  return status;
}

myb_status myb_orbit_kick(myb_orbit* orbit, const double dv[3]) {
  MYB_CHECK(orbit != NULL && orbit->started && dv != NULL);
  myb_orbit* o = orbit;
  if (!isfinite(dv[0]) || !isfinite(dv[1]) || !isfinite(dv[2])) {
    return MYB_E_DOMAIN;
  }
  const double t = o->reached;
  double y[6];
  myb_status status = MYB_OK;
  if (o->var != NULL) {
    double phi[36];
    double sensitivity[6 * MYB_MODEL_MAX_PARAMETERS];
    if (o->have_step) {
      status = myb_variational_dense(o->var, t, y, phi, sensitivity);
    } else {
      myb_variational_get(o->var, y, phi, sensitivity);
    }
    if (status != MYB_OK) {
      return status;
    }
    for (int i = 0; i < 3; ++i) {
      y[3 + i] += dv[i];
    }
    status = myb_internal_variational_restart(o->var, t, y, phi, sensitivity);
  } else {
    status = dense(o, t, y);
    if (status != MYB_OK) {
      return status;
    }
    for (int i = 0; i < 3; ++i) {
      y[3 + i] += dv[i];
    }
    if (o->wide != 0) { /* the variations go on from where they are: an impulse changes none */
      if (o->have_step) {
        status = myb_dop853_dense(o->dop, t, o->wide_scratch);
        if (status != MYB_OK) {
          return status;
        }
        memcpy(o->wide_start, o->wide_scratch, o->wide * sizeof(double));
      }
      memcpy(o->wide_start, y, 6 * sizeof(double));
    }
    status = o->dop != NULL
                 ? start_dop853(o, t, y)
                 : myb_ias15_start(o->ias, myb_force_model_accel, o->model, t, y);
    o->approaching = 0;
  }
  if (status != MYB_OK) {
    o->started = 0;
    return status;
  }
  memcpy(o->y0, y, sizeof o->y0);
  o->have_step = 0;
  o->step_begin = t;
  o->step_end = t;
  const size_t count = myb_force_model_break_count(o->model);
  for (size_t k = 0; k < count; ++k) {
    myb_break_function g = NULL;
    void* context = NULL;
    myb_force_model_break(o->model, k, &g, &context);
    o->break_sign[k] = sign_of(g(context, t, y));
  }
  return MYB_OK;
}

/* A root of g between (ta, ga) and (tb, gb)? Counted when the signs are strictly opposite or
 * the far end is a zero, and g changes with t in the direction wanted. */
static int crossing(const myb_orbit* o, double ga, double gb, int direction) {
  const int changes = (ga < 0.0 && gb > 0.0) || (ga > 0.0 && gb < 0.0) || (ga != 0.0 && gb == 0.0);
  if (!changes || direction == 0) {
    return changes;
  }
  const int grows = o->direction > 0.0 ? ga < 0.0 : ga > 0.0;
  return direction > 0 ? grows : !grows;
}

/* Locates the root in [ta, tb] and moves it, if needed, to the first representable time where
 * g has left the sign it had at ta: a later search from the root then does not find it
 * again. */
static myb_status locate(myb_orbit* o, myb_event_function g, void* context, double ta,
                         const double ya[6], double ga, double tb, const double yb[6],
                         int* found, double* t_event, double y_event[6]) {
  const myb_internal_step_view view = {.dense = dense_view,
                                       .integrator = o,
                                       .m = 6,
                                       .t0 = ta,
                                       .y0 = ya,
                                       .t1 = tb,
                                       .y1 = yb,
                                       .scratch = o->scratch};
  double te = 0.0;
  double ye[6];
  myb_status status = myb_internal_locate_event(&view, g, context, found, &te, ye);
  if (status != MYB_OK || !*found) {
    return status;
  }
  for (int nudge = 0; nudge < EVENT_NUDGES && te != tb; ++nudge) {
    const double ge = g(context, te, ye);
    if (ge == 0.0 || (ge < 0.0) != (ga < 0.0)) {
      break;
    }
    te = nextafter(te, tb);
    status = dense(o, te, ye);
    if (status != MYB_OK) {
      return status;
    }
  }
  *t_event = te;
  memcpy(y_event, ye, sizeof ye);
  return MYB_OK;
}

myb_status myb_orbit_event(myb_orbit* orbit, myb_event_function g, void* context, double t_limit,
                           int direction, int* found, double* t_event, double y_event[6]) {
  MYB_CHECK(orbit != NULL && orbit->started && g != NULL && found != NULL && t_event != NULL);
  myb_orbit* o = orbit;
  *found = 0;
  myb_status status = in_order(o, t_limit);
  if (status != MYB_OK) {
    return status;
  }
  double ta = o->reached;
  double ya[6];
  status = dense(o, ta, ya);
  if (status != MYB_OK) {
    return status;
  }
  double ga = g(context, ta, ya);
  int64_t budget = o->max_steps;
  /* bound: the step budget; every pass examines what is left of a step, then takes the next */
  for (;;) {
    const int last = (t_limit - o->step_end) * o->direction <= 0.0;
    const double stop = last ? t_limit : o->step_end;
    for (int part = 1; part <= EVENT_PARTS && (stop - ta) * o->direction > 0.0; ++part) {
      double tb = part < EVENT_PARTS
                      ? o->step_begin + (o->step_end - o->step_begin) * ((double)part / event_parts)
                      : o->step_end;
      if ((tb - ta) * o->direction <= 0.0) {
        continue;
      }
      if ((tb - stop) * o->direction > 0.0) {
        tb = stop;
      }
      double yb[6];
      status = dense(o, tb, yb);
      if (status != MYB_OK) {
        return status;
      }
      const double gb = g(context, tb, yb);
      if (crossing(o, ga, gb, direction)) {
        double ye[6];
        status = locate(o, g, context, ta, ya, ga, tb, yb, found, t_event, ye);
        if (status != MYB_OK) {
          return status;
        }
        if (*found) {
          o->reached = *t_event;
          if (y_event != NULL) {
            memcpy(y_event, ye, sizeof ye);
          }
          return MYB_OK;
        }
      }
      ta = tb;
      memcpy(ya, yb, sizeof ya);
      ga = gb;
    }
    if (last) {
      o->reached = t_limit;
      return MYB_OK;
    }
    if (budget == 0) {
      return MYB_E_NOCONV;
    }
    --budget;
    status = advance(o);
    if (status != MYB_OK) {
      return status;
    }
  }
}
