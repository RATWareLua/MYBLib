/* SPDX-License-Identifier: BSD-3-Clause */
/* Ready-made terms of a force model: the functions of myb/grav.h, myb/force.h and myb/rel.h
 * fed from an environment. */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "env/env_internal.h"
#include "myb/env.h"
#include "myb/model.h"
#include "num/dual.h"
#include "num/vec3.h"

static void sub3(const double a[3], const double b[3], double out[3]) {
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}

static void zero3(double a[3]) {
  a[0] = 0.0;
  a[1] = 0.0;
  a[2] = 0.0;
}

static int valid_object(const myb_object* object) {
  return isfinite(object->mass) && object->mass > 0.0 && isfinite(object->area) &&
         isfinite(object->cr) && isfinite(object->cd);
}

/* The mass of the object at t; 0 if its function gives nothing usable. */
static double mass_of(const myb_object* object, double t) {
  if (object->mass_function == NULL) {
    return object->mass;
  }
  const double mass = object->mass_function(object->mass_context, t);
  return mass > 0.0 && isfinite(mass) ? mass : 0.0;
}

double myb_internal_object_mass(const myb_object* object, double t) {
  return isfinite(object->mass) && object->mass > 0.0 ? mass_of(object, t) : 0.0;
}

/* ---- Point masses ---------------------------------------------------------------------------- */

myb_status myb_term_point_mass(void* term, double t, const double x[3], const double v[3],
                               double a[3]) {
  const myb_point_mass_term* p = term;
  MYB_CHECK(p != NULL && x != NULL && v != NULL && a != NULL);
  (void)t;
  const double r = myb_vec3_norm(x);
  if (!(r > 0.0) || !isfinite(r) || !isfinite(p->gm)) {
    return MYB_E_DOMAIN;
  }
  const double k = -p->gm / (r * r * r);
  for (int i = 0; i < 3; ++i) {
    a[i] = k * x[i];
  }
  return MYB_OK;
}

myb_status myb_term_point_mass_jacobian(void* term, double t, const double x[3],
                                        const double v[3], double dadx[9], double dadv[9]) {
  const myb_point_mass_term* p = term;
  MYB_CHECK(p != NULL && x != NULL && v != NULL && dadx != NULL && dadv != NULL);
  (void)t;
  const double r = myb_vec3_norm(x);
  if (!(r > 0.0) || !isfinite(r) || !isfinite(p->gm)) {
    return MYB_E_DOMAIN;
  }
  myb_point_mass_jacobian(p->gm, x, dadx);
  memset(dadv, 0, 9 * sizeof(double));
  return MYB_OK;
}

/* The body of a third-body term and its position; MYB_E_DOMAIN for the body at the origin. */
static myb_status third_body(const myb_third_body_term* tb, double t, int* origin,
                             const double** s, double* gm) {
  *origin = myb_internal_env_origin(tb->env);
  *gm = myb_internal_env_body(tb->env, tb->body)->gm;
  if (*origin >= 0 && (size_t)*origin == tb->body) {
    return MYB_E_DOMAIN;
  }
  return myb_internal_env_state(tb->env, t, tb->body, s);
}

myb_status myb_term_third_body(void* term, double t, const double x[3], const double v[3],
                               double a[3]) {
  const myb_third_body_term* tb = term;
  MYB_CHECK(tb != NULL && tb->env != NULL && x != NULL && v != NULL && a != NULL);
  int origin = -1;
  const double* s = NULL;
  double gm = 0.0;
  const myb_status status = third_body(tb, t, &origin, &s, &gm);
  if (status != MYB_OK) {
    return status;
  }
  double d[3];
  sub3(s, x, d);
  const double dd = myb_vec3_norm(d);
  if (!(dd > 0.0) || !isfinite(dd)) {
    return MYB_E_DOMAIN;
  }
  zero3(a);
  if (origin >= 0) {
    if (!(myb_vec3_norm(s) > 0.0)) {
      return MYB_E_DOMAIN;
    }
    myb_force_third_body(gm, s, x, a);
    return MYB_OK;
  }
  const double k = gm / (dd * dd * dd);
  for (int i = 0; i < 3; ++i) {
    a[i] = k * d[i];
  }
  return MYB_OK;
}

myb_status myb_term_third_body_jacobian(void* term, double t, const double x[3],
                                        const double v[3], double dadx[9], double dadv[9]) {
  const myb_third_body_term* tb = term;
  MYB_CHECK(tb != NULL && tb->env != NULL && x != NULL && v != NULL && dadx != NULL &&
            dadv != NULL);
  int origin = -1;
  const double* s = NULL;
  double gm = 0.0;
  const myb_status status = third_body(tb, t, &origin, &s, &gm);
  if (status != MYB_OK) {
    return status;
  }
  double r[3];
  sub3(x, s, r);
  const double rr = myb_vec3_norm(r);
  if (!(rr > 0.0) || !isfinite(rr)) {
    return MYB_E_DOMAIN;
  }
  /* the indirect part does not depend on the object */
  myb_point_mass_jacobian(gm, r, dadx);
  memset(dadv, 0, 9 * sizeof(double));
  return MYB_OK;
}

/* ---- Gravity field ----------------------------------------------------------------------------- */

/* Brings the coefficients of the field to the instant t: time-variable terms and tides. */
static myb_status update_field(myb_field_term* f, double t, const double* s,
                               const myb_orientation* o) {
  myb_time when;
  myb_status status = myb_internal_env_orientation_time(f->env, t, f->body, &when);
  if (status != MYB_OK) {
    return status;
  }
  myb_gravity_field_set_time(f->field, when);
  if (f->tides == NULL) {
    return MYB_OK;
  }
  myb_tide_raiser raisers[MYB_ENV_MAX_RAISERS];
  for (size_t k = 0; k < f->raiser_count; ++k) {
    const double* sr = NULL;
    status = myb_internal_env_state(f->env, t, f->raisers[k], &sr);
    if (status != MYB_OK) {
      return status;
    }
    double d[3];
    sub3(sr, s, d);
    raisers[k].gm = myb_internal_env_body(f->env, f->raisers[k])->gm;
    myb_rotation_apply_position(&o->rotation, d, raisers[k].r);
  }
  const myb_tides_input input = {.raisers = raisers,
                                 .count = f->raiser_count,
                                 .t = o->ut,
                                 .xp = o->xp,
                                 .yp = o->yp,
                                 .flags = f->tide_flags,
                                 .permanent = f->permanent};
  return myb_tides_apply(f->tides, f->field, &input);
}

/* The field at the inertial position r relative to its body, in inertial axes. */
static myb_status field_at(const myb_field_term* f, const myb_orientation* o,
                           const myb_rotation* inverse, const double r[3], double g[3]) {
  double r_fixed[3];
  myb_rotation_apply_position(&o->rotation, r, r_fixed);
  double g_fixed[3];
  const myb_status status =
      myb_gravity_field_eval(f->field, f->work, r_fixed, f->n_min, f->n_max, NULL, g_fixed);
  if (status != MYB_OK) {
    return status;
  }
  myb_rotation_apply_position(inverse, g_fixed, g);
  return MYB_OK;
}

/* The state and orientation of the field's body at t, with the coefficients brought to t. */
static myb_status field_prepare(myb_field_term* f, double t, const double** s,
                                const myb_orientation** o, const myb_rotation** inverse) {
  myb_status status = myb_internal_env_state(f->env, t, f->body, s);
  if (status != MYB_OK) {
    return status;
  }
  status = myb_internal_env_orientation(f->env, t, f->body, o, inverse);
  if (status != MYB_OK) {
    return status;
  }
  if ((f->update || f->tides != NULL) && !(f->cached && f->cached_t == t)) {
    f->cached = 0;
    status = update_field(f, t, *s, *o);
    if (status != MYB_OK) {
      return status;
    }
    f->cached = 1;
    f->cached_t = t;
  }
  return MYB_OK;
}

myb_status myb_term_field(void* term, double t, const double x[3], const double v[3],
                          double a[3]) {
  myb_field_term* f = term;
  MYB_CHECK(f != NULL && f->env != NULL && f->field != NULL && f->work != NULL && x != NULL &&
            v != NULL && a != NULL && f->raiser_count <= MYB_ENV_MAX_RAISERS);
  const double* s = NULL;
  const myb_orientation* o = NULL;
  const myb_rotation* inverse = NULL;
  myb_status status = field_prepare(f, t, &s, &o, &inverse);
  if (status != MYB_OK) {
    return status;
  }
  double r[3];
  sub3(x, s, r);
  status = field_at(f, o, inverse, r, a);
  if (status != MYB_OK) {
    return status;
  }
  const int origin = myb_internal_env_origin(f->env);
  if (origin >= 0 && (size_t)origin != f->body) {
    const double r0[3] = {-s[0], -s[1], -s[2]};
    double g0[3];
    status = field_at(f, o, inverse, r0, g0);
    if (status != MYB_OK) {
      return status;
    }
    sub3(a, g0, a);
  }
  return MYB_OK;
}

myb_status myb_term_field_jacobian(void* term, double t, const double x[3], const double v[3],
                                   double dadx[9], double dadv[9]) {
  myb_field_term* f = term;
  MYB_CHECK(f != NULL && f->env != NULL && f->field != NULL && f->work != NULL && x != NULL &&
            v != NULL && dadx != NULL && dadv != NULL &&
            f->raiser_count <= MYB_ENV_MAX_RAISERS);
  const double* s = NULL;
  const myb_orientation* o = NULL;
  const myb_rotation* inverse = NULL;
  myb_status status = field_prepare(f, t, &s, &o, &inverse);
  if (status != MYB_OK) {
    return status;
  }
  double r[3];
  double r_fixed[3];
  sub3(x, s, r);
  myb_rotation_apply_position(&o->rotation, r, r_fixed);
  double g[9];
  status = myb_gravity_field_gradient(f->field, f->work, r_fixed, f->n_min, f->n_max, g);
  if (status != MYB_OK) {
    return status;
  }
  /* a = M^T g(M r): d a / d x = M^T G M; the indirect term does not depend on the object */
  const double* m = o->rotation.m;
  double gm[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      gm[3 * i + j] = g[3 * i] * m[j] + g[3 * i + 1] * m[3 + j] + g[3 * i + 2] * m[6 + j];
    }
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      dadx[3 * i + j] = m[i] * gm[j] + m[3 + i] * gm[3 + j] + m[6 + i] * gm[6 + j];
    }
  }
  memset(dadv, 0, 9 * sizeof(double));
  return MYB_OK;
}

/* ---- Radiation --------------------------------------------------------------------------------- */

/* The part of the disk of the source left visible at x by the occulters; d = x - (the
 * source). One occulter at a time gives the product of the separate fractions, as before;
 * where two or more airless ones cover parts of the disk at once, their union. */
static myb_status shadowing(myb_env* env, double t, size_t source_index, const size_t* occulters,
                            size_t count, const myb_refraction* atmospheres, const double x[3],
                            const double d[3], double* visible) {
  *visible = 1.0;
  if (count == 0) {
    return MYB_OK;
  }
  const myb_env_body* source = myb_internal_env_body(env, source_index);
  if (!(source->radius > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double to_source[3] = {-d[0], -d[1], -d[2]};
  const double* limb = (source->limb[0] != 0.0 || source->limb[1] != 0.0) ? source->limb : NULL;
  double product = 1.0;
  double with_air = 1.0; /* of the occulters with an atmosphere */
  double towards[3 * MYB_ENV_MAX_OCCULTERS];
  double radii[MYB_ENV_MAX_OCCULTERS];
  size_t partial = 0;
  for (size_t k = 0; k < count; ++k) {
    const myb_env_body* occulter = myb_internal_env_body(env, occulters[k]);
    if (!(occulter->radius > 0.0) || occulters[k] == source_index) {
      return MYB_E_DOMAIN;
    }
    const double* so = NULL;
    const myb_status status = myb_internal_env_state(env, t, occulters[k], &so);
    if (status != MYB_OK) {
      return status;
    }
    double to_occulter[3];
    sub3(so, x, to_occulter);
    if (atmospheres != NULL && atmospheres[k].scale_height > 0.0) {
      const double f = myb_shadow_fraction_refracted(to_source, source->radius, to_occulter,
                                                     occulter->radius, &atmospheres[k], limb);
      product *= f;
      with_air *= f;
      continue;
    }
    const double f =
        myb_shadow_fraction(to_source, source->radius, to_occulter, occulter->radius, limb);
    product *= f;
    if (f > 0.0 && f < 1.0) {
      memcpy(towards + 3 * partial, to_occulter, sizeof to_occulter);
      radii[partial] = occulter->radius;
      ++partial;
    }
  }
  if (partial >= 2 && product > 0.0) {
    product = with_air *
              myb_shadow_fraction_many(to_source, source->radius, towards, radii, partial, limb);
  }
  *visible = product;
  return MYB_OK;
}

myb_status myb_term_radiation(void* term, double t, const double x[3], const double v[3],
                              double a[3]) {
  const myb_radiation_term* r = term;
  MYB_CHECK(r != NULL && r->env != NULL && r->object != NULL && x != NULL && v != NULL &&
            a != NULL && r->occulter_count <= MYB_ENV_MAX_OCCULTERS);
  const myb_object* object = r->object;
  if (!valid_object(object) || !(r->c > 0.0) || !isfinite(r->c)) {
    return MYB_E_DOMAIN;
  }
  const myb_env_body* source = myb_internal_env_body(r->env, r->source);
  const double* s = NULL;
  myb_status status = myb_internal_env_state(r->env, t, r->source, &s);
  if (status != MYB_OK) {
    return status;
  }
  double d[3];
  sub3(x, s, d);
  const double dn = myb_vec3_norm(d);
  if (!(dn > 0.0) || !isfinite(dn)) {
    return MYB_E_DOMAIN;
  }
  const double u[3] = {d[0] / dn, d[1] / dn, d[2] / dn};
  double visible = 1.0;
  status = shadowing(r->env, t, r->source, r->occulters, r->occulter_count, r->atmospheres, x,
                     d, &visible);
  if (status != MYB_OK) {
    return status;
  }
  const double pressure = myb_radiation_flux(source->luminosity, dn) / r->c;
  const double mass = mass_of(object, t);
  if (!(mass > 0.0)) {
    return MYB_E_DOMAIN;
  }
  zero3(a);
  if (r->moving) {
    double w[3];
    sub3(v, s + 3, w);
    myb_force_radiation_moving(visible * pressure, object->cr, object->area / mass, u, w, r->c,
                               a);
  } else {
    myb_force_radiation_sphere(visible * pressure, object->cr, object->area / mass, u, a);
  }
  return MYB_OK;
}

/* Central differences of a term with respect to the position, steps h. */
static myb_status position_differences(myb_term_function function, void* term, double t,
                                       const double x[3], const double v[3], double h,
                                       double dadx[9]) {
  for (int j = 0; j < 3; ++j) {
    double plus[3] = {x[0], x[1], x[2]};
    double minus[3] = {x[0], x[1], x[2]};
    plus[j] += h;
    minus[j] -= h;
    double ap[3];
    double am[3];
    myb_status status = function(term, t, plus, v, ap);
    if (status == MYB_OK) {
      status = function(term, t, minus, v, am);
    }
    if (status != MYB_OK) {
      return status;
    }
    const double step = plus[j] - minus[j];
    for (int i = 0; i < 3; ++i) {
      dadx[3 * i + j] = (ap[i] - am[i]) / step;
    }
  }
  return MYB_OK;
}

myb_status myb_term_radiation_jacobian(void* term, double t, const double x[3],
                                       const double v[3], double dadx[9], double dadv[9]) {
  const myb_radiation_term* r = term;
  MYB_CHECK(r != NULL && r->env != NULL && r->object != NULL && x != NULL && v != NULL &&
            dadx != NULL && dadv != NULL && r->occulter_count <= MYB_ENV_MAX_OCCULTERS);
  const myb_object* object = r->object;
  if (!valid_object(object) || !(r->c > 0.0) || !isfinite(r->c)) {
    return MYB_E_DOMAIN;
  }
  const myb_env_body* source = myb_internal_env_body(r->env, r->source);
  const double* s = NULL;
  myb_status status = myb_internal_env_state(r->env, t, r->source, &s);
  if (status != MYB_OK) {
    return status;
  }
  double d[3];
  sub3(x, s, d);
  const double dn = myb_vec3_norm(d);
  if (!(dn > 0.0) || !isfinite(dn)) {
    return MYB_E_DOMAIN;
  }
  double visible = 1.0;
  status = shadowing(r->env, t, r->source, r->occulters, r->occulter_count, r->atmospheres, x,
                     d, &visible);
  if (status != MYB_OK) {
    return status;
  }
  memset(dadv, 0, 9 * sizeof(double));
  const double mass = mass_of(object, t);
  if (!(mass > 0.0)) {
    return MYB_E_DOMAIN;
  }
  /* a = k d / |d|^3 with k = visible L / (4 pi c) cr A / m */
  const double k = visible * (myb_radiation_flux(source->luminosity, dn) * dn * dn / r->c) *
                   object->cr * (object->area / mass);
  double w[3] = {0.0, 0.0, 0.0};
  const double u[3] = {d[0] / dn, d[1] / dn, d[2] / dn};
  if (r->moving) {
    /* the part in the velocity: -(k / (c |d|^2)) (u u^T + 1) */
    sub3(v, s + 3, w);
    const double scale = k / (r->c * dn * dn);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        dadv[3 * i + j] = -scale * (u[i] * u[j] + (i == j ? 1.0 : 0.0));
      }
    }
  }
  if (visible > 0.0 && visible < 1.0) {
    const double reach = myb_vec3_norm(x);
    const double h = 6.0554544523933395e-6 * (reach > 0.0 ? reach : dn);
    return position_differences(myb_term_radiation, term, t, x, v, h, dadx);
  }
  /* minus the Jacobian of a point mass k at the source */
  myb_point_mass_jacobian(-k, d, dadx);
  if (r->moving) {
    /* of -(k / c) [(u.w) u + w] / |d|^2:
     * -(k / (c |d|^3)) [u w^T + (u.w) (1 - 4 u u^T)] + (2 k / (c |d|^3)) w u^T */
    const double scale = k / (r->c * dn * dn * dn);
    const double uw = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        dadx[3 * i + j] += scale * (2.0 * w[i] * u[j] - u[i] * w[j] -
                                    uw * ((i == j ? 1.0 : 0.0) - 4.0 * u[i] * u[j]));
      }
    }
  }
  return MYB_OK;
}

myb_status myb_term_planet_radiation(void* term, double t, const double x[3], const double v[3],
                                     double a[3]) {
  const myb_planet_radiation_term* p = term;
  MYB_CHECK(p != NULL && p->env != NULL && p->model != NULL && p->object != NULL && x != NULL &&
            v != NULL && a != NULL);
  const myb_object* object = p->object;
  if (!valid_object(object) || !isfinite(p->phase0)) {
    return MYB_E_DOMAIN;
  }
  const double* sp = NULL;
  myb_status status = myb_internal_env_state(p->env, t, p->planet, &sp);
  if (status != MYB_OK) {
    return status;
  }
  const double* ss = NULL;
  status = myb_internal_env_state(p->env, t, p->source, &ss);
  if (status != MYB_OK) {
    return status;
  }
  const myb_orientation* o = NULL;
  const myb_rotation* inverse = NULL;
  status = myb_internal_env_orientation(p->env, t, p->planet, &o, &inverse);
  if (status != MYB_OK) {
    return status;
  }
  double r[3];
  double sun[3];
  sub3(x, sp, r);
  sub3(ss, sp, sun);
  const double rs = myb_vec3_norm(sun);
  if (!(rs > 0.0) || !isfinite(rs)) {
    return MYB_E_DOMAIN;
  }
  double r_fixed[3];
  double sun_fixed[3];
  myb_rotation_apply_position(&o->rotation, r, r_fixed);
  myb_rotation_apply_position(&o->rotation, sun, sun_fixed);
  const double flux =
      myb_radiation_flux(myb_internal_env_body(p->env, p->source)->luminosity, rs);
  const double phase = p->phase0 + p->model->omega * t;
  double a_fixed[3] = {0.0, 0.0, 0.0};
  const double mass = mass_of(object, t);
  if (!(mass > 0.0)) {
    return MYB_E_DOMAIN;
  }
  status = myb_force_planet_radiation(p->model, r_fixed, sun_fixed, flux, phase,
                                      object->cr * (object->area / mass), p->rings,
                                      a_fixed);
  if (status != MYB_OK) {
    return status;
  }
  myb_rotation_apply_position(inverse, a_fixed, a);
  return MYB_OK;
}

/* ---- Atmosphere --------------------------------------------------------------------------------- */

myb_status myb_term_drag(void* term, double t, const double x[3], const double v[3],
                         double a[3]) {
  const myb_drag_term* d = term;
  MYB_CHECK(d != NULL && d->env != NULL && d->density != NULL && d->object != NULL &&
            x != NULL && v != NULL && a != NULL);
  const myb_object* object = d->object;
  if (!valid_object(object)) {
    return MYB_E_DOMAIN;
  }
  const double* s = NULL;
  myb_status status = myb_internal_env_state(d->env, t, d->body, &s);
  if (status != MYB_OK) {
    return status;
  }
  const myb_orientation* o = NULL;
  const myb_rotation* inverse = NULL;
  status = myb_internal_env_orientation(d->env, t, d->body, &o, &inverse);
  if (status != MYB_OK) {
    return status;
  }
  myb_time when;
  status = myb_internal_env_orientation_time(d->env, t, d->body, &when);
  if (status != MYB_OK) {
    return status;
  }
  double y[6];
  sub3(x, s, y);
  sub3(v, s + 3, y + 3);
  /* the velocity in the body-fixed frame is the velocity relative to the co-rotating gas */
  double fixed[6];
  myb_rotation_apply_state(&o->rotation, y, fixed);
  double density = 0.0;
  status = d->density(d->density_context, when, fixed, &density);
  if (status != MYB_OK) {
    return status;
  }
  if (!(density >= 0.0) || !isfinite(density)) {
    return MYB_E_DOMAIN;
  }
  const double mass = mass_of(object, t);
  if (!(mass > 0.0)) {
    return MYB_E_DOMAIN;
  }
  if (d->wind != NULL) {
    double wind[3];
    status = d->wind(d->wind_context, when, fixed, wind);
    if (status != MYB_OK) {
      return status;
    }
    if (!myb_vec3_finite(wind)) {
      return MYB_E_DOMAIN;
    }
    sub3(fixed + 3, wind, fixed + 3);
  }
  double w[3];
  myb_rotation_apply_position(inverse, fixed + 3, w);
  zero3(a);
  myb_force_drag(density, object->cd * (object->area / mass), w, a);
  return MYB_OK;
}

myb_status myb_term_drag_jacobian(void* term, double t, const double x[3], const double v[3],
                                  double dadx[9], double dadv[9]) {
  const myb_drag_term* d = term;
  MYB_CHECK(d != NULL && d->env != NULL && d->density != NULL && d->object != NULL &&
            x != NULL && v != NULL && dadx != NULL && dadv != NULL);
  const myb_object* object = d->object;
  if (!valid_object(object)) {
    return MYB_E_DOMAIN;
  }
  const double* s = NULL;
  myb_status status = myb_internal_env_state(d->env, t, d->body, &s);
  if (status != MYB_OK) {
    return status;
  }
  const myb_orientation* o = NULL;
  const myb_rotation* inverse = NULL;
  status = myb_internal_env_orientation(d->env, t, d->body, &o, &inverse);
  if (status != MYB_OK) {
    return status;
  }
  myb_time when;
  status = myb_internal_env_orientation_time(d->env, t, d->body, &when);
  if (status != MYB_OK) {
    return status;
  }
  double y[6];
  sub3(x, s, y);
  sub3(v, s + 3, y + 3);
  double fixed[6];
  myb_rotation_apply_state(&o->rotation, y, fixed);
  double density = 0.0;
  status = d->density(d->density_context, when, fixed, &density);
  if (status != MYB_OK) {
    return status;
  }
  if (!(density >= 0.0) || !isfinite(density)) {
    return MYB_E_DOMAIN;
  }
  const double mass = mass_of(object, t);
  if (!(mass > 0.0)) {
    return MYB_E_DOMAIN;
  }
  double relative[3] = {fixed[3], fixed[4], fixed[5]};
  if (d->wind != NULL) {
    double wind[3];
    status = d->wind(d->wind_context, when, fixed, wind);
    if (status != MYB_OK) {
      return status;
    }
    if (!myb_vec3_finite(wind)) {
      return MYB_E_DOMAIN;
    }
    sub3(relative, wind, relative);
  }
  double w[3];
  myb_rotation_apply_position(inverse, relative, w);
  const double speed = myb_vec3_norm(w);
  const double beta = object->cd * (object->area / mass);
  /* a = -1/2 rho beta |w| w with w = (v - s') + M^T dM (x - s):
   * d a / d w = -1/2 rho beta (|w| 1 + w w^T / |w|), d w / d v = 1, d w / d x = M^T dM */
  double daw[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double outer = speed > 0.0 ? w[i] * w[j] / speed : 0.0;
      daw[3 * i + j] = -0.5 * density * beta * ((i == j ? speed : 0.0) + outer);
    }
  }
  memcpy(dadv, daw, 9 * sizeof(double));
  const double* m = o->rotation.m;
  const double* dm = o->rotation.dm;
  double spin[9]; /* M^T dM */
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      spin[3 * i + j] = m[i] * dm[j] + m[3 + i] * dm[3 + j] + m[6 + i] * dm[6 + j];
    }
  }
  /* the gradient of the density in the body-fixed frame, turned into the inertial axes */
  const double reach = myb_vec3_norm(fixed);
  const double h = 1e-7 * reach;
  double gradient_fixed[3] = {0.0, 0.0, 0.0};
  if (h > 0.0) {
    for (int j = 0; j < 3; ++j) {
      double plus[3] = {fixed[0], fixed[1], fixed[2]};
      double minus[3] = {fixed[0], fixed[1], fixed[2]};
      plus[j] += h;
      minus[j] -= h;
      double rho_plus = 0.0;
      double rho_minus = 0.0;
      status = d->density(d->density_context, when, plus, &rho_plus);
      if (status == MYB_OK) {
        status = d->density(d->density_context, when, minus, &rho_minus);
      }
      if (status != MYB_OK) {
        return status;
      }
      gradient_fixed[j] = (rho_plus - rho_minus) / (plus[j] - minus[j]);
    }
  }
  double gradient[3];
  myb_rotation_apply_position(inverse, gradient_fixed, gradient);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      dadx[3 * i + j] = daw[3 * i] * spin[j] + daw[3 * i + 1] * spin[3 + j] +
                        daw[3 * i + 2] * spin[6 + j] - 0.5 * beta * speed * w[i] * gradient[j];
    }
  }
  return MYB_OK;
}

/* ---- Relativity ---------------------------------------------------------------------------------- */

/* Position and velocity of the object relative to the central body. */
static myb_status about_central(const myb_relativity_term* r, double t, const double x[3],
                                const double v[3], double rr[3], double vv[3]) {
  const double* s = NULL;
  const myb_status status = myb_internal_env_state(r->env, t, r->central, &s);
  if (status != MYB_OK) {
    return status;
  }
  sub3(x, s, rr);
  sub3(v, s + 3, vv);
  return MYB_OK;
}

myb_status myb_term_schwarzschild(void* term, double t, const double x[3], const double v[3],
                                  double a[3]) {
  const myb_relativity_term* r = term;
  MYB_CHECK(r != NULL && r->env != NULL && x != NULL && v != NULL && a != NULL);
  double rr[3];
  double vv[3];
  const myb_status status = about_central(r, t, x, v, rr, vv);
  if (status != MYB_OK) {
    return status;
  }
  return myb_rel_schwarzschild(&r->ppn, myb_internal_env_body(r->env, r->central)->gm, rr, vv, a);
}

myb_status myb_term_lense_thirring(void* term, double t, const double x[3], const double v[3],
                                   double a[3]) {
  const myb_relativity_term* r = term;
  MYB_CHECK(r != NULL && r->env != NULL && x != NULL && v != NULL && a != NULL);
  double rr[3];
  double vv[3];
  myb_status status = about_central(r, t, x, v, rr, vv);
  if (status != MYB_OK) {
    return status;
  }
  const myb_orientation* o = NULL;
  const myb_rotation* inverse = NULL;
  status = myb_internal_env_orientation(r->env, t, r->central, &o, &inverse);
  if (status != MYB_OK) {
    return status;
  }
  const myb_env_body* central = myb_internal_env_body(r->env, r->central);
  double j[3];
  myb_rotation_apply_position(inverse, central->spin, j);
  return myb_rel_lense_thirring(&r->ppn, central->gm, j, rr, vv, a);
}

myb_status myb_term_de_sitter(void* term, double t, const double x[3], const double v[3],
                              double a[3]) {
  const myb_relativity_term* r = term;
  MYB_CHECK(r != NULL && r->env != NULL && x != NULL && v != NULL && a != NULL);
  if (r->central == r->external) {
    return MYB_E_DOMAIN;
  }
  double rr[3];
  double vv[3];
  myb_status status = about_central(r, t, x, v, rr, vv);
  if (status != MYB_OK) {
    return status;
  }
  const double* sc = NULL;
  status = myb_internal_env_state(r->env, t, r->central, &sc);
  if (status != MYB_OK) {
    return status;
  }
  const double* se = NULL;
  status = myb_internal_env_state(r->env, t, r->external, &se);
  if (status != MYB_OK) {
    return status;
  }
  double big_r[3];
  double big_v[3];
  sub3(sc, se, big_r);
  sub3(sc + 3, se + 3, big_v);
  return myb_rel_de_sitter(&r->ppn, myb_internal_env_body(r->env, r->external)->gm, big_r, big_v,
                           vv, a);
}

myb_status myb_term_schwarzschild_jacobian(void* term, double t, const double x[3],
                                           const double v[3], double dadx[9], double dadv[9]) {
  const myb_relativity_term* r = term;
  MYB_CHECK(r != NULL && r->env != NULL && x != NULL && v != NULL && dadx != NULL &&
            dadv != NULL);
  double rr[3];
  double vv[3];
  const myb_status status = about_central(r, t, x, v, rr, vv);
  if (status != MYB_OK) {
    return status;
  }
  return myb_rel_schwarzschild_jacobian(&r->ppn, myb_internal_env_body(r->env, r->central)->gm,
                                        rr, vv, dadx, dadv);
}

myb_status myb_term_lense_thirring_jacobian(void* term, double t, const double x[3],
                                            const double v[3], double dadx[9], double dadv[9]) {
  const myb_relativity_term* r = term;
  MYB_CHECK(r != NULL && r->env != NULL && x != NULL && v != NULL && dadx != NULL &&
            dadv != NULL);
  double rr[3];
  double vv[3];
  myb_status status = about_central(r, t, x, v, rr, vv);
  if (status != MYB_OK) {
    return status;
  }
  const myb_orientation* o = NULL;
  const myb_rotation* inverse = NULL;
  status = myb_internal_env_orientation(r->env, t, r->central, &o, &inverse);
  if (status != MYB_OK) {
    return status;
  }
  const myb_env_body* central = myb_internal_env_body(r->env, r->central);
  double j[3];
  myb_rotation_apply_position(inverse, central->spin, j);
  return myb_rel_lense_thirring_jacobian(&r->ppn, central->gm, j, rr, vv, dadx, dadv);
}

myb_status myb_term_de_sitter_jacobian(void* term, double t, const double x[3],
                                       const double v[3], double dadx[9], double dadv[9]) {
  const myb_relativity_term* r = term;
  MYB_CHECK(r != NULL && r->env != NULL && x != NULL && v != NULL && dadx != NULL &&
            dadv != NULL);
  if (r->central == r->external) {
    return MYB_E_DOMAIN;
  }
  const double* sc = NULL;
  myb_status status = myb_internal_env_state(r->env, t, r->central, &sc);
  if (status != MYB_OK) {
    return status;
  }
  const double* se = NULL;
  status = myb_internal_env_state(r->env, t, r->external, &se);
  if (status != MYB_OK) {
    return status;
  }
  double big_r[3];
  double big_v[3];
  sub3(sc, se, big_r);
  sub3(sc + 3, se + 3, big_v);
  memset(dadx, 0, 9 * sizeof(double));
  return myb_rel_de_sitter_jacobian(&r->ppn, myb_internal_env_body(r->env, r->external)->gm,
                                    big_r, big_v, dadv);
}

/* ---- Empirical ----------------------------------------------------------------------------------- */

myb_status myb_term_empirical(void* term, double t, const double x[3], const double v[3],
                              double a[3]) {
  const myb_empirical_term* e = term;
  MYB_CHECK(e != NULL && x != NULL && v != NULL && a != NULL);
  (void)t;
  double h[3];
  myb_vec3_cross(x, v, h);
  const double hn = myb_vec3_norm(h);
  const double r = myb_vec3_norm(x);
  if (!(hn > 0.0) || !(r > 0.0) || !isfinite(hn) || !isfinite(r)) {
    return MYB_E_DOMAIN;
  }
  /* the ascending node on the xy plane: z x h; an orbit in that plane counts from x */
  double node[3] = {-h[1], h[0], 0.0};
  const double nn = sqrt(node[0] * node[0] + node[1] * node[1]);
  if (nn > 0.0) {
    node[0] /= nn;
    node[1] /= nn;
  } else {
    node[0] = 1.0;
    node[1] = 0.0;
  }
  const double pole[3] = {h[0] / hn, h[1] / hn, h[2] / hn};
  double across[3];
  myb_vec3_cross(pole, node, across);
  const double cu = myb_vec3_dot(node, x) / r;
  const double su = myb_vec3_dot(across, x) / r;
  double rtn[3];
  for (int i = 0; i < 3; ++i) {
    rtn[i] = e->constant[i] + e->cosine[i] * cu + e->sine[i] * su;
  }
  zero3(a);
  return myb_force_empirical_rtn(x, v, rtn, a);
}

/* ---- Small bodies, the lagging tide ------------------------------------------------------------ */

myb_status myb_term_nongrav(void* term, double t, const double x[3], const double v[3],
                            double a[3]) {
  const myb_nongrav_term* n = term;
  MYB_CHECK(n != NULL && n->env != NULL && x != NULL && v != NULL && a != NULL);
  const double* s = NULL;
  const myb_status status = myb_internal_env_state(n->env, t, n->source, &s);
  if (status != MYB_OK) {
    return status;
  }
  double r[3];
  double w[3];
  sub3(x, s, r);
  sub3(v, s + 3, w);
  zero3(a);
  if (n->delay != 0.0) {
    return myb_force_nongrav_delayed(&n->model, myb_internal_env_body(n->env, n->source)->gm,
                                     n->delay, r, w, a);
  }
  return myb_force_nongrav(&n->model, r, w, a);
}

myb_status myb_term_tidal_lag(void* term, double t, const double x[3], const double v[3],
                              double a[3]) {
  const myb_tidal_lag_term* l = term;
  MYB_CHECK(l != NULL && l->env != NULL && x != NULL && v != NULL && a != NULL);
  const myb_env_body* body = myb_internal_env_body(l->env, l->body);
  const double* s = NULL;
  myb_status status = myb_internal_env_state(l->env, t, l->body, &s);
  if (status != MYB_OK) {
    return status;
  }
  const myb_orientation* o = NULL;
  status = myb_internal_env_orientation(l->env, t, l->body, &o, NULL);
  if (status != MYB_OK) {
    return status;
  }
  /* the rotation vector in the inertial axes: [omega x] = dM^T M */
  const double* m = o->rotation.m;
  const double* dm = o->rotation.dm;
  const double omega[3] = {dm[2] * m[1] + dm[5] * m[4] + dm[8] * m[7],
                           dm[0] * m[2] + dm[3] * m[5] + dm[6] * m[8],
                           dm[1] * m[0] + dm[4] * m[3] + dm[7] * m[6]};
  double r[3];
  double w[3];
  sub3(x, s, r);
  sub3(v, s + 3, w);
  zero3(a);
  return myb_force_tidal_lag(body->gm, l->gm_object, l->k2, body->radius, l->lag, omega, r, w,
                             l->with_static, a);
}

/* ---- What the terms of the surface share (surface/surface.c) ----------------------------------- */

myb_status myb_internal_env_light(myb_env* env, double t, size_t source_index,
                                  const size_t* occulters, size_t occulter_count,
                                  const myb_refraction* atmospheres, const double x[3],
                                  double d[3], double* distance, double* flux) {
  const myb_env_body* source = myb_internal_env_body(env, source_index);
  const double* s = NULL;
  myb_status status = myb_internal_env_state(env, t, source_index, &s);
  if (status != MYB_OK) {
    return status;
  }
  sub3(x, s, d);
  const double dn = myb_vec3_norm(d);
  if (!(dn > 0.0) || !isfinite(dn)) {
    return MYB_E_DOMAIN;
  }
  double visible = 1.0;
  status = shadowing(env, t, source_index, occulters, occulter_count, atmospheres, x, d,
                     &visible);
  if (status != MYB_OK) {
    return status;
  }
  *distance = dn;
  *flux = visible * myb_radiation_flux(source->luminosity, dn);
  return MYB_OK;
}

/* ---- More Jacobians, more terms ---------------------------------------------------------------- */

/* Central differences of a term in the position and the velocity. */
static myb_status state_differences(myb_term_function function, void* term, double t,
                                    const double x[3], const double v[3], double hx, double hv,
                                    double dadx[9], double dadv[9]) {
  for (int j = 0; j < 6; ++j) {
    double q[6] = {x[0], x[1], x[2], v[0], v[1], v[2]};
    const double h = j < 3 ? hx : hv;
    double plus[3];
    double minus[3];
    q[j] += h;
    myb_status status = function(term, t, q, q + 3, plus);
    if (status != MYB_OK) {
      return status;
    }
    q[j] -= 2.0 * h;
    status = function(term, t, q, q + 3, minus);
    if (status != MYB_OK) {
      return status;
    }
    double* out = j < 3 ? dadx : dadv;
    for (int i = 0; i < 3; ++i) {
      out[3 * i + j % 3] = (plus[i] - minus[i]) / (2.0 * h);
    }
  }
  return MYB_OK;
}

myb_status myb_term_planet_radiation_jacobian(void* term, double t, const double x[3],
                                              const double v[3], double dadx[9],
                                              double dadv[9]) {
  const myb_planet_radiation_term* p = term;
  MYB_CHECK(p != NULL && p->env != NULL && p->model != NULL && p->object != NULL && x != NULL &&
            v != NULL && dadx != NULL && dadv != NULL);
  const myb_object* object = p->object;
  if (!valid_object(object) || !isfinite(p->phase0)) {
    return MYB_E_DOMAIN;
  }
  const double* sp = NULL;
  myb_status status = myb_internal_env_state(p->env, t, p->planet, &sp);
  if (status != MYB_OK) {
    return status;
  }
  const double* ss = NULL;
  status = myb_internal_env_state(p->env, t, p->source, &ss);
  if (status != MYB_OK) {
    return status;
  }
  const myb_orientation* o = NULL;
  status = myb_internal_env_orientation(p->env, t, p->planet, &o, NULL);
  if (status != MYB_OK) {
    return status;
  }
  double r[3];
  double sun[3];
  sub3(x, sp, r);
  sub3(ss, sp, sun);
  const double rs = myb_vec3_norm(sun);
  const double mass = mass_of(object, t);
  if (!(rs > 0.0) || !isfinite(rs) || !(mass > 0.0)) {
    return MYB_E_DOMAIN;
  }
  double r_fixed[3];
  double sun_fixed[3];
  myb_rotation_apply_position(&o->rotation, r, r_fixed);
  myb_rotation_apply_position(&o->rotation, sun, sun_fixed);
  const double flux =
      myb_radiation_flux(myb_internal_env_body(p->env, p->source)->luminosity, rs);
  double g[9];
  status = myb_force_planet_radiation_gradient(p->model, r_fixed, sun_fixed, flux,
                                               p->phase0 + p->model->omega * t,
                                               object->cr * (object->area / mass), p->rings, g);
  if (status != MYB_OK) {
    return status;
  }
  /* M^T G M */
  const double* m = o->rotation.m;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double sum = 0.0;
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          sum += m[3 * k + i] * g[3 * k + l] * m[3 * l + j];
        }
      }
      dadx[3 * i + j] = sum;
    }
  }
  memset(dadv, 0, 9 * sizeof(double));
  return MYB_OK;
}

/* r / |r|, the unit normal of the plane of (r, w) and the direction across r in that plane,
 * as myb_force_empirical_rtn and myb_force_nongrav take them. */
static int orbital_frame(const myb_dual r[3], const myb_dual w[3], myb_dual radial[3],
                         myb_dual transverse[3], myb_dual normal[3], myb_dual* rn) {
  myb_dual h[3];
  myb_dual_cross(r, w, h);
  const myb_dual hn = myb_dual_norm(h);
  *rn = myb_dual_norm(r);
  if (!(hn.v > 0.0) || !(rn->v > 0.0) || !isfinite(hn.v) || !isfinite(rn->v)) {
    return 0;
  }
  for (int i = 0; i < 3; ++i) {
    radial[i] = myb_dual_div(r[i], *rn);
    normal[i] = myb_dual_div(h[i], hn);
  }
  myb_dual_cross(normal, radial, transverse);
  return 1;
}

myb_status myb_term_empirical_jacobian(void* term, double t, const double x[3],
                                       const double v[3], double dadx[9], double dadv[9]) {
  const myb_empirical_term* e = term;
  MYB_CHECK(e != NULL && x != NULL && v != NULL && dadx != NULL && dadv != NULL);
  (void)t;
  myb_dual r[3];
  myb_dual w[3];
  myb_dual_state(x, v, NULL, r, w);
  myb_dual radial[3];
  myb_dual transverse[3];
  myb_dual normal[3];
  myb_dual rn;
  if (!orbital_frame(r, w, radial, transverse, normal, &rn)) {
    return MYB_E_DOMAIN;
  }
  /* the ascending node z x h and the direction a quarter of a turn on in the orbit; the
   * normal has the direction of h */
  myb_dual node[3] = {myb_dual_scale(normal[1], -1.0), normal[0], myb_dual_constant(0.0)};
  const myb_dual n2 =
      myb_dual_add(myb_dual_mul(node[0], node[0]), myb_dual_mul(node[1], node[1]));
  if (n2.v > 0.0) {
    const myb_dual nn = myb_dual_sqrt(n2);
    node[0] = myb_dual_div(node[0], nn);
    node[1] = myb_dual_div(node[1], nn);
  } else {
    node[0] = myb_dual_constant(1.0);
    node[1] = myb_dual_constant(0.0);
  }
  myb_dual across[3];
  myb_dual_cross(normal, node, across);
  const myb_dual cu = myb_dual_dot(node, radial);
  const myb_dual su = myb_dual_dot(across, radial);
  myb_dual a[3] = {myb_dual_constant(0.0), myb_dual_constant(0.0), myb_dual_constant(0.0)};
  const myb_dual* axes[3] = {radial, transverse, normal};
  for (int k = 0; k < 3; ++k) {
    const myb_dual size =
        myb_dual_add(myb_dual_constant(e->constant[k]),
                     myb_dual_add(myb_dual_scale(cu, e->cosine[k]),
                                  myb_dual_scale(su, e->sine[k])));
    for (int i = 0; i < 3; ++i) {
      a[i] = myb_dual_add(a[i], myb_dual_mul(size, axes[k][i]));
    }
  }
  myb_dual_jacobian(a, dadx, dadv);
  return MYB_OK;
}

myb_status myb_term_nongrav_jacobian(void* term, double t, const double x[3], const double v[3],
                                     double dadx[9], double dadv[9]) {
  const myb_nongrav_term* n = term;
  MYB_CHECK(n != NULL && n->env != NULL && x != NULL && v != NULL && dadx != NULL &&
            dadv != NULL);
  const double* s = NULL;
  const myb_status status = myb_internal_env_state(n->env, t, n->source, &s);
  if (status != MYB_OK) {
    return status;
  }
  if (n->delay != 0.0) {
    const double d[3] = {x[0] - s[0], x[1] - s[1], x[2] - s[2]};
    const double w[3] = {v[0] - s[3], v[1] - s[4], v[2] - s[5]};
    return state_differences(myb_term_nongrav, term, t, x, v,
                             6.0554544523933395e-6 * myb_vec3_norm(d),
                             6.0554544523933395e-6 * myb_vec3_norm(w), dadx, dadv);
  }
  myb_dual r[3];
  myb_dual w[3];
  myb_dual_state(x, v, s, r, w);
  myb_dual radial[3];
  myb_dual transverse[3];
  myb_dual normal[3];
  myb_dual rn;
  const myb_nongrav* m = &n->model;
  if (!orbital_frame(r, w, radial, transverse, normal, &rn) || !(m->r0 > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const myb_dual q = myb_dual_scale(rn, 1.0 / m->r0);
  myb_dual g = myb_dual_scale(myb_dual_power(q, -m->m), m->alpha);
  if (m->k != 0.0) {
    g = myb_dual_mul(g, myb_dual_power(myb_dual_add(myb_dual_constant(1.0),
                                                    myb_dual_power(q, m->n)),
                                       -m->k));
  }
  myb_dual a[3];
  for (int i = 0; i < 3; ++i) {
    const myb_dual direction =
        myb_dual_add(myb_dual_scale(radial[i], m->a[0]),
                     myb_dual_add(myb_dual_scale(transverse[i], m->a[1]),
                                  myb_dual_scale(normal[i], m->a[2])));
    a[i] = myb_dual_mul(g, direction);
  }
  myb_dual_jacobian(a, dadx, dadv);
  return MYB_OK;
}

myb_status myb_term_tidal_lag_jacobian(void* term, double t, const double x[3],
                                       const double v[3], double dadx[9], double dadv[9]) {
  const myb_tidal_lag_term* l = term;
  MYB_CHECK(l != NULL && l->env != NULL && x != NULL && v != NULL && dadx != NULL &&
            dadv != NULL);
  const myb_env_body* body = myb_internal_env_body(l->env, l->body);
  const double* s = NULL;
  myb_status status = myb_internal_env_state(l->env, t, l->body, &s);
  if (status != MYB_OK) {
    return status;
  }
  const myb_orientation* o = NULL;
  status = myb_internal_env_orientation(l->env, t, l->body, &o, NULL);
  if (status != MYB_OK) {
    return status;
  }
  const double* m = o->rotation.m;
  const double* dm = o->rotation.dm;
  const myb_dual omega[3] = {
      myb_dual_constant(dm[2] * m[1] + dm[5] * m[4] + dm[8] * m[7]),
      myb_dual_constant(dm[0] * m[2] + dm[3] * m[5] + dm[6] * m[8]),
      myb_dual_constant(dm[1] * m[0] + dm[4] * m[3] + dm[7] * m[6])};
  myb_dual r[3];
  myb_dual w[3];
  myb_dual_state(x, v, s, r, w);
  const myb_dual r2 = myb_dual_dot(r, r);
  if (!(body->gm > 0.0) || !(body->radius > 0.0) || !(r2.v > 0.0) || !isfinite(r2.v)) {
    return MYB_E_DOMAIN;
  }
  const double radius2 = body->radius * body->radius;
  const double strength = 3.0 * l->k2 * l->gm_object * (1.0 + l->gm_object / body->gm) *
                          (radius2 * radius2 * body->radius);
  const myb_dual r4 = myb_dual_mul(r2, r2);
  const myb_dual scale = myb_dual_div(myb_dual_constant(strength), myb_dual_mul(r4, r4));
  myb_dual spin[3];
  myb_dual_cross(r, omega, spin);
  const myb_dual rv = myb_dual_dot(r, w);
  myb_dual a[3];
  for (int i = 0; i < 3; ++i) {
    const myb_dual along = myb_dual_div(myb_dual_scale(myb_dual_mul(rv, r[i]), 2.0), r2);
    myb_dual sum = myb_dual_scale(myb_dual_add(along, myb_dual_add(spin[i], w[i])), l->lag);
    if (l->with_static) {
      sum = myb_dual_add(sum, r[i]);
    }
    a[i] = myb_dual_scale(myb_dual_mul(scale, sum), -1.0);
  }
  myb_dual_jacobian(a, dadx, dadv);
  return MYB_OK;
}

myb_status myb_term_yarkovsky(void* term, double t, const double x[3], const double v[3],
                              double a[3]) {
  const myb_yarkovsky_term* y = term;
  MYB_CHECK(y != NULL && y->env != NULL && x != NULL && v != NULL && a != NULL);
  const myb_env_body* source = myb_internal_env_body(y->env, y->source);
  const double* s = NULL;
  const myb_status status = myb_internal_env_state(y->env, t, y->source, &s);
  if (status != MYB_OK) {
    return status;
  }
  double d[3];
  sub3(s, x, d);
  const double dn = myb_vec3_norm(d);
  if (!(dn > 0.0) || !isfinite(dn)) {
    return MYB_E_DOMAIN;
  }
  const double to_sun[3] = {d[0] / dn, d[1] / dn, d[2] / dn};
  zero3(a);
  if (y->season != NULL) {
    const double from_sun[3] = {-d[0], -d[1], -d[2]};
    myb_force_yarkovsky_season(y->season, from_sun, a);
  }
  return myb_force_yarkovsky(&y->body, y->spin, myb_radiation_flux(source->luminosity, dn), y->c,
                             to_sun, a);
}

myb_status myb_term_field_indirect(void* term, double t, const double x[3], const double v[3],
                                   double a[3]) {
  myb_field_term* f = term;
  MYB_CHECK(f != NULL && f->env != NULL && f->field != NULL && f->work != NULL && x != NULL &&
            v != NULL && a != NULL && f->raiser_count <= MYB_ENV_MAX_RAISERS);
  const double* s = NULL;
  const myb_orientation* o = NULL;
  const myb_rotation* inverse = NULL;
  myb_status status = field_prepare(f, t, &s, &o, &inverse);
  if (status != MYB_OK) {
    return status;
  }
  const myb_env_body* body = myb_internal_env_body(f->env, f->body);
  const int origin = myb_internal_env_origin(f->env);
  if (origin < 0 || (size_t)origin != f->body || !(body->gm > 0.0)) {
    return MYB_E_DOMAIN;
  }
  zero3(a);
  const int n_min = f->n_min > 2 ? f->n_min : 2;
  for (size_t k = 0; k < f->raiser_count; ++k) {
    const myb_env_body* other = myb_internal_env_body(f->env, f->raisers[k]);
    const double* so = NULL;
    status = myb_internal_env_state(f->env, t, f->raisers[k], &so);
    if (status != MYB_OK) {
      return status;
    }
    double r_fixed[3];
    double g_fixed[3];
    double g[3];
    myb_rotation_apply_position(&o->rotation, so, r_fixed);
    status = myb_gravity_field_eval(f->field, f->work, r_fixed, n_min, f->n_max, NULL, g_fixed);
    if (status != MYB_OK) {
      return status;
    }
    myb_rotation_apply_position(inverse, g_fixed, g);
    const double ratio = other->gm / body->gm;
    for (int i = 0; i < 3; ++i) {
      a[i] += ratio * g[i];
    }
  }
  return MYB_OK;
}
