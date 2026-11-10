/* SPDX-License-Identifier: BSD-3-Clause */
/* The surface of the object (myb/surface.h): free-molecular flow on plates, attitude laws and
 * the terms that turn plates and meshes with them. The light on a mesh is with the mesh, in
 * grav/polyhedron.c. */
#include <math.h>
#include <string.h>

#include "core/base.h"
#include "env/env_internal.h"
#include "force/surface_internal.h"
#include "myb/cmath.h"
#include "myb/surface.h"
#include "num/dual.h"
#include "num/erf.h"
#include "num/vec3.h"

static const double gas_constant = 8.314462618; /* J mol^-1 K^-1, exact since 2019 */
static const double sqrt_pi = 1.7724538509055160273;

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

/* ---- Free-molecular flow ---------------------------------------------------------------------- */

static int valid_gas(const myb_gas* gas) {
  return gas->temperature > 0.0 && gas->molar_mass > 0.0 && gas->wall_temperature >= 0.0 &&
         gas->accommodation >= 0.0 && gas->accommodation <= 1.0 && isfinite(gas->temperature) &&
         isfinite(gas->molar_mass) && isfinite(gas->wall_temperature);
}

/* The speed ratio S and the ratio of the speeds of the molecules leaving and arriving. */
static void flow_ratios(const myb_gas* gas, double speed, double* s, double* leaving) {
  const double thermal2 = 2.0 * gas_constant * gas->temperature / gas->molar_mass;
  *s = speed / sqrt(thermal2);
  const double wall = 4.0 * gas_constant * gas->wall_temperature / (gas->molar_mass * speed * speed);
  *leaving = sqrt(0.5 * (1.0 + gas->accommodation * (wall - 1.0)));
}

myb_status myb_gas_sphere_cd(const myb_gas* gas, double speed, double* cd) {
  MYB_CHECK(gas != NULL && cd != NULL);
  if (!valid_gas(gas) || !(speed > 0.0) || !isfinite(speed)) {
    return MYB_E_DOMAIN;
  }
  double s = 0.0;
  double leaving = 0.0;
  flow_ratios(gas, speed, &s, &leaving);
  const double s2 = s * s;
  *cd = (2.0 * s2 + 1.0) / (sqrt_pi * s2 * s) * myb_exp(-s2) +
        (4.0 * s2 * s2 + 4.0 * s2 - 1.0) / (2.0 * s2 * s2) * myb_internal_erf(s) +
        2.0 * sqrt_pi / 3.0 * leaving;
  return MYB_OK;
}

myb_status myb_force_gas_plates(double density, const myb_gas* gas, double mass,
                                const myb_plate* plates, size_t count, const double v[3],
                                double a[3]) {
  MYB_CHECK(gas != NULL && (plates != NULL || count == 0) && v != NULL && a != NULL);
  const double speed = myb_vec3_norm(v);
  if (!valid_gas(gas) || !(mass > 0.0) || !(density >= 0.0) || !(speed > 0.0) ||
      !isfinite(mass) || !isfinite(density) || !isfinite(speed)) {
    return MYB_E_DOMAIN;
  }
  double s = 0.0;
  double leaving = 0.0;
  flow_ratios(gas, speed, &s, &leaving);
  const double g = 1.0 / (2.0 * s * s);
  const double q = 1.0 + g;
  /* the gas arrives along -v: the drag is along it */
  const double drag[3] = {-v[0] / speed, -v[1] / speed, -v[2] / speed};
  double sum[3] = {0.0, 0.0, 0.0};
  for (size_t k = 0; k < count; ++k) {
    const myb_plate* p = &plates[k];
    const double n2 = myb_vec3_dot(p->normal, p->normal);
    if (!(p->area >= 0.0) || !isfinite(p->area) || !(fabs(n2 - 1.0) <= 1e-12)) {
      return MYB_E_DOMAIN;
    }
    const double gamma = -myb_vec3_dot(drag, p->normal);
    /* the part of the normal across the flow: the lift pushes against it */
    double across[3];
    for (int i = 0; i < 3; ++i) {
      across[i] = p->normal[i] + gamma * drag[i];
    }
    const double l = myb_vec3_norm(across);
    const double pe = myb_exp(-gamma * gamma * s * s) / s;
    const double z = 1.0 + myb_internal_erf(gamma * s);
    const double reemitted = 0.5 * leaving * (gamma * sqrt_pi * z + pe);
    const double cd = pe / sqrt_pi + gamma * q * z + gamma * reemitted;
    const double cl = l * g * z + l * reemitted;
    for (int i = 0; i < 3; ++i) {
      sum[i] += p->area * (cd * drag[i] - (l > 0.0 ? cl * across[i] / l : 0.0));
    }
  }
  const double scale = 0.5 * density * speed * speed / mass;
  for (int i = 0; i < 3; ++i) {
    a[i] += scale * sum[i];
  }
  return MYB_OK;
}

double myb_gas_knudsen(double density, double molar_mass, double diameter, double length) {
  static const double avogadro = 6.02214076e23; /* mol^-1, exact since 2019 */
  if (!(density > 0.0) || !(molar_mass > 0.0) || !(diameter > 0.0) || !(length > 0.0)) {
    return (double)INFINITY;
  }
  const double path =
      molar_mass / (1.4142135623730951 * 3.141592653589793 * diameter * diameter * density *
                    avogadro);
  return path / length;
}

double myb_gas_bridge(double knudsen) {
  if (!(knudsen > 1.0e-3)) {
    return 0.0;
  }
  if (knudsen >= 10.0) {
    return 1.0;
  }
  const double decades = myb_log(knudsen) / 2.302585092994046;
  const double s = myb_sin(3.141592653589793 * (3.0 + decades) / 8.0);
  return s * s;
}

/* ---- Attitude --------------------------------------------------------------------------------- */

myb_status myb_attitude_fixed(void* matrix, double t, const double x[3], const double v[3],
                              double m[9]) {
  MYB_CHECK(matrix != NULL && m != NULL);
  (void)t;
  (void)x;
  (void)v;
  memcpy(m, matrix, 9 * sizeof(double));
  return MYB_OK;
}

/* Rows x, y, z of the frame with x along a and z along a x b. */
static myb_status frame_of(const double a[3], const double b[3], double m[9]) {
  const double an = myb_vec3_norm(a);
  double h[3];
  myb_vec3_cross(a, b, h);
  const double hn = myb_vec3_norm(h);
  if (!(an > 0.0) || !(hn > 0.0) || !isfinite(an) || !isfinite(hn)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    m[i] = a[i] / an;
    m[6 + i] = h[i] / hn;
  }
  myb_vec3_cross(m + 6, m, m + 3);
  return MYB_OK;
}

myb_status myb_attitude_orbital(void* context, double t, const double x[3], const double v[3],
                                double m[9]) {
  const myb_orbital_attitude* o = context;
  MYB_CHECK(o != NULL && o->env != NULL && x != NULL && v != NULL && m != NULL);
  const double* s = NULL;
  const myb_status status = myb_internal_env_state(o->env, t, o->body, &s);
  if (status != MYB_OK) {
    return status;
  }
  double r[3];
  double w[3];
  sub3(x, s, r);
  sub3(v, s + 3, w);
  return frame_of(r, w, m);
}

myb_status myb_attitude_sun_nadir(void* context, double t, const double x[3], const double v[3],
                                  double m[9]) {
  const myb_sun_attitude* a = context;
  MYB_CHECK(a != NULL && a->env != NULL && x != NULL && m != NULL);
  (void)v;
  const double* s = NULL;
  myb_status status = myb_internal_env_state(a->env, t, a->body, &s);
  if (status != MYB_OK) {
    return status;
  }
  double down[3];
  sub3(s, x, down);
  status = myb_internal_env_state(a->env, t, a->source, &s);
  if (status != MYB_OK) {
    return status;
  }
  double sun[3];
  sub3(s, x, sun);
  /* z down, y = z x sun, x = y x z */
  double z[3];
  double y[3];
  const double dn = myb_vec3_norm(down);
  myb_vec3_cross(down, sun, y);
  const double yn = myb_vec3_norm(y);
  /* in line to the rounding of the vectors: no axis for the panels */
  if (!(dn > 0.0) || !(yn > 1.0e-12 * dn * myb_vec3_norm(sun)) || !isfinite(dn) ||
      !isfinite(yn)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    z[i] = down[i] / dn;
    y[i] = y[i] / yn;
  }
  myb_vec3_cross(y, z, m);
  memcpy(m + 3, y, sizeof y);
  memcpy(m + 6, z, sizeof z);
  return MYB_OK;
}

myb_status myb_attitude_spin(void* context, double t, const double x[3], const double v[3],
                             double m[9]) {
  const myb_spin_attitude* s = context;
  MYB_CHECK(s != NULL && m != NULL);
  (void)x;
  (void)v;
  const double angle = s->phase + s->rate * t;
  if (!isfinite(angle)) {
    return MYB_E_DOMAIN;
  }
  const double c = myb_cos(angle);
  const double sn = myb_sin(angle);
  for (int j = 0; j < 3; ++j) {
    m[j] = c * s->m0[j] + sn * s->m0[3 + j];
    m[3 + j] = -sn * s->m0[j] + c * s->m0[3 + j];
    m[6 + j] = s->m0[6 + j];
  }
  return MYB_OK;
}

/* out = M in, out = M^T in. */
static void to_body(const double m[9], const double in[3], double out[3]) {
  for (int i = 0; i < 3; ++i) {
    out[i] = m[3 * i] * in[0] + m[3 * i + 1] * in[1] + m[3 * i + 2] * in[2];
  }
}

static void from_body(const double m[9], const double in[3], double out[3]) {
  for (int i = 0; i < 3; ++i) {
    out[i] = m[i] * in[0] + m[3 + i] * in[1] + m[6 + i] * in[2];
  }
}

static int finite9(const double m[9]) {
  for (int i = 0; i < 9; ++i) {
    if (!isfinite(m[i])) {
      return 0;
    }
  }
  return 1;
}

/* ---- Terms ------------------------------------------------------------------------------------ */

myb_status myb_term_surface_radiation(void* term, double t, const double x[3], const double v[3],
                                      double a[3]) {
  const myb_surface_radiation_term* r = term;
  MYB_CHECK(r != NULL && r->env != NULL && r->object != NULL && r->attitude != NULL &&
            x != NULL && v != NULL && a != NULL && r->occulter_count <= MYB_ENV_MAX_OCCULTERS &&
            (r->shape == NULL ? (r->plates != NULL || r->plate_count == 0)
                              : r->surface != NULL));
  const double mass = myb_internal_object_mass(r->object, t);
  if (!(mass > 0.0) || !(r->c > 0.0) || !isfinite(r->c)) {
    return MYB_E_DOMAIN;
  }
  double d[3];
  double distance = 0.0;
  double flux = 0.0;
  myb_status status = myb_internal_env_light(r->env, t, r->source, r->occulters,
                                             r->occulter_count, r->atmospheres, x, d, &distance,
                                             &flux);
  if (status != MYB_OK) {
    return status;
  }
  double m[9];
  status = r->attitude(r->attitude_context, t, x, v, m);
  if (status != MYB_OK) {
    return status;
  }
  if (!finite9(m)) {
    return MYB_E_DOMAIN;
  }
  const double u[3] = {d[0] / distance, d[1] / distance, d[2] / distance};
  double u_body[3];
  to_body(m, u, u_body);
  double a_body[3] = {0.0, 0.0, 0.0};
  if (r->shape != NULL) {
    double force[3];
    status = myb_force_radiation_mesh(r->shape, r->surface, flux / r->c, u_body, r->shadows,
                                      force, NULL);
    for (int i = 0; i < 3; ++i) {
      a_body[i] = force[i] / mass;
    }
  } else {
    status = myb_force_radiation_plates(flux / r->c, mass, r->plates, r->plate_count, u_body,
                                        a_body);
  }
  if (status != MYB_OK) {
    return status;
  }
  from_body(m, a_body, a);
  return MYB_OK;
}

myb_status myb_term_surface_drag(void* term, double t, const double x[3], const double v[3],
                                 double a[3]) {
  const myb_surface_drag_term* d = term;
  MYB_CHECK(d != NULL && d->env != NULL && d->density != NULL && d->object != NULL &&
            d->attitude != NULL && (d->plates != NULL || d->plate_count == 0) && x != NULL &&
            v != NULL && a != NULL);
  const double mass = myb_internal_object_mass(d->object, t);
  if (!(mass > 0.0)) {
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
  /* the velocity relative to the gas in the inertial axes, then in the body of the object */
  double w[3];
  myb_rotation_apply_position(inverse, fixed + 3, w);
  double m[9];
  status = d->attitude(d->attitude_context, t, x, v, m);
  if (status != MYB_OK) {
    return status;
  }
  if (!finite9(m)) {
    return MYB_E_DOMAIN;
  }
  double w_body[3];
  to_body(m, w, w_body);
  double a_body[3] = {0.0, 0.0, 0.0};
  zero3(a);
  if (density == 0.0) {
    return MYB_OK;
  }
  status = myb_force_gas_plates(density, &d->gas, mass, d->plates, d->plate_count, w_body,
                                a_body);
  if (status != MYB_OK) {
    return status;
  }
  from_body(m, a_body, a);
  if (d->length > 0.0) {
    if (!(d->molecule > 0.0) || !isfinite(d->continuum_cd) || !isfinite(d->continuum_area)) {
      return MYB_E_DOMAIN;
    }
    const double weight =
        myb_gas_bridge(myb_gas_knudsen(density, d->gas.molar_mass, d->molecule, d->length));
    const double k = -0.5 * density * d->continuum_cd * d->continuum_area / mass *
                     myb_vec3_norm(w);
    for (int i = 0; i < 3; ++i) {
      a[i] = weight * a[i] + (1.0 - weight) * k * w[i];
    }
  }
  return MYB_OK;
}

myb_status myb_term_recoil(void* term, double t, const double x[3], const double v[3],
                           double a[3]) {
  const myb_recoil_term* r = term;
  MYB_CHECK(r != NULL && r->object != NULL && r->attitude != NULL && x != NULL && v != NULL &&
            a != NULL);
  const double mass = myb_internal_object_mass(r->object, t);
  const double n2 = myb_vec3_dot(r->direction, r->direction);
  if (!(mass > 0.0) || !(r->c > 0.0) || !isfinite(r->c) || !isfinite(r->power) ||
      !isfinite(r->efficiency) || !(fabs(n2 - 1.0) <= 1e-12) || !(r->decay >= 0.0) ||
      !isfinite(r->epoch)) {
    return MYB_E_DOMAIN;
  }
  double m[9];
  const myb_status status = r->attitude(r->attitude_context, t, x, v, m);
  if (status != MYB_OK) {
    return status;
  }
  if (!finite9(m)) {
    return MYB_E_DOMAIN;
  }
  const double power = r->decay > 0.0 ? r->power * myb_exp(-(t - r->epoch) / r->decay) : r->power;
  double a_body[3] = {0.0, 0.0, 0.0};
  myb_force_emission(power, r->efficiency, mass, r->c, r->direction, a_body);
  from_body(m, a_body, a);
  return MYB_OK;
}

/* ---- Jacobians --------------------------------------------------------------------------------
 * The terms are written once more on numbers that carry their derivatives in the state of the
 * object (num/dual.h): a = M^T f(M w, ...) with the attitude M, the vector w that the force
 * in the body depends on, and the scalars (flux, density). */

static const double difference_step = 6.0554544523933395e-6; /* the cube root of 2^-52 */

/* out = M in and out = M^T in for an attitude with its derivatives. */
static void to_body_dual(const myb_dual m[9], const myb_dual in[3], myb_dual out[3]) {
  for (int i = 0; i < 3; ++i) {
    out[i] = myb_dual_dot(m + 3 * i, in);
  }
}

static void from_body_dual(const myb_dual m[9], const myb_dual in[3], myb_dual out[3]) {
  for (int i = 0; i < 3; ++i) {
    const myb_dual column[3] = {m[i], m[3 + i], m[6 + i]};
    out[i] = myb_dual_dot(column, in);
  }
}

/* Rows x, y, z of the frame with x along a and z along a x b, as frame_of. */
static myb_status frame_of_dual(const myb_dual a[3], const myb_dual b[3], myb_dual m[9]) {
  const myb_dual an = myb_dual_norm(a);
  myb_dual h[3];
  myb_dual_cross(a, b, h);
  const myb_dual hn = myb_dual_norm(h);
  if (!(an.v > 0.0) || !(hn.v > 0.0) || !isfinite(an.v) || !isfinite(hn.v)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 3; ++i) {
    m[i] = myb_dual_div(a[i], an);
    m[6 + i] = myb_dual_div(h[i], hn);
  }
  myb_dual_cross(m + 6, m, m + 3);
  return MYB_OK;
}

/* The attitude with its derivatives in the state. The laws of this file are differentiated
 * as they are written; any other law by central differences of the law. */
static myb_status attitude_dual(myb_attitude_function attitude, void* context, double t,
                                const double x[3], const double v[3], myb_dual m[9]) {
  if (attitude == myb_attitude_orbital) {
    const myb_orbital_attitude* o = context;
    MYB_CHECK(o != NULL && o->env != NULL);
    const double* s = NULL;
    const myb_status status = myb_internal_env_state(o->env, t, o->body, &s);
    if (status != MYB_OK) {
      return status;
    }
    myb_dual r[3];
    myb_dual w[3];
    myb_dual_state(x, v, s, r, w);
    return frame_of_dual(r, w, m);
  }
  if (attitude == myb_attitude_sun_nadir) {
    const myb_sun_attitude* a = context;
    MYB_CHECK(a != NULL && a->env != NULL);
    const double* s = NULL;
    myb_status status = myb_internal_env_state(a->env, t, a->body, &s);
    if (status != MYB_OK) {
      return status;
    }
    myb_dual down[3];
    for (int i = 0; i < 3; ++i) {
      down[i] = myb_dual_constant(s[i] - x[i]);
      down[i].d[i] = -1.0;
    }
    status = myb_internal_env_state(a->env, t, a->source, &s);
    if (status != MYB_OK) {
      return status;
    }
    myb_dual sun[3];
    for (int i = 0; i < 3; ++i) {
      sun[i] = myb_dual_constant(s[i] - x[i]);
      sun[i].d[i] = -1.0;
    }
    /* z down, y = z x sun, x = y x z */
    myb_dual y[3];
    myb_dual_cross(down, sun, y);
    const myb_dual dn = myb_dual_norm(down);
    const myb_dual yn = myb_dual_norm(y);
    const double sn = sqrt(sun[0].v * sun[0].v + sun[1].v * sun[1].v + sun[2].v * sun[2].v);
    if (!(dn.v > 0.0) || !(yn.v > 1.0e-12 * dn.v * sn) || !isfinite(dn.v) || !isfinite(yn.v)) {
      return MYB_E_DOMAIN;
    }
    for (int i = 0; i < 3; ++i) {
      m[6 + i] = myb_dual_div(down[i], dn);
      m[3 + i] = myb_dual_div(y[i], yn);
    }
    myb_dual_cross(m + 3, m + 6, m);
    return MYB_OK;
  }
  double now[9];
  myb_status status = attitude(context, t, x, v, now);
  if (status != MYB_OK) {
    return status;
  }
  if (!finite9(now)) {
    return MYB_E_DOMAIN;
  }
  for (int k = 0; k < 9; ++k) {
    m[k] = myb_dual_constant(now[k]);
  }
  if (attitude == myb_attitude_fixed || attitude == myb_attitude_spin) {
    return MYB_OK; /* they do not look at the state */
  }
  const double reach = myb_vec3_norm(x);
  const double speed = myb_vec3_norm(v);
  const double hx = difference_step * (reach > 0.0 ? reach : 1.0);
  const double hv = difference_step * (speed > 0.0 ? speed : 1.0);
  for (int j = 0; j < 6; ++j) {
    double q[6] = {x[0], x[1], x[2], v[0], v[1], v[2]};
    const double h = j < 3 ? hx : hv;
    double plus[9];
    double minus[9];
    q[j] += h;
    const double up = q[j];
    status = attitude(context, t, q, q + 3, plus);
    if (status != MYB_OK) {
      return status;
    }
    q[j] -= 2.0 * h;
    status = attitude(context, t, q, q + 3, minus);
    if (status != MYB_OK) {
      return status;
    }
    for (int k = 0; k < 9; ++k) {
      m[k].d[j] = (plus[k] - minus[k]) / (up - q[j]);
    }
  }
  return MYB_OK;
}

/* myb_force_radiation_plates on such numbers: pressure is the light's pressure over the mass. */
static myb_status radiation_plates_dual(myb_dual pressure, const myb_plate* plates, size_t count,
                                        const myb_dual u[3], myb_dual a[3]) {
  for (int i = 0; i < 3; ++i) {
    a[i] = myb_dual_constant(0.0);
  }
  for (size_t k = 0; k < count; ++k) {
    const myb_plate* p = &plates[k];
    const double nn = myb_vec3_dot(p->normal, p->normal);
    if (!(p->area >= 0.0) || !isfinite(p->area) || !(fabs(nn - 1.0) <= 2e-12) ||
        !(p->specular >= 0.0 && p->specular <= 1.0) ||
        !(p->diffuse >= 0.0 && p->diffuse <= 1.0) || p->specular + p->diffuse > 1.0) {
      return MYB_E_DOMAIN;
    }
    const myb_dual n[3] = {myb_dual_constant(p->normal[0]), myb_dual_constant(p->normal[1]),
                           myb_dual_constant(p->normal[2])};
    const myb_dual c = myb_dual_scale(myb_dual_dot(n, u), -1.0);
    if (!(c.v > 0.0)) {
      continue; /* lit from behind or edge-on */
    }
    const double absorbed = 1.0 - p->specular - p->diffuse;
    myb_dual normal = myb_dual_scale(c, 2.0 * p->specular);
    normal.v += 2.0 * p->diffuse / 3.0 + (p->thermal ? 2.0 / 3.0 * absorbed : 0.0);
    const myb_dual lit = myb_dual_scale(myb_dual_mul(pressure, c), p->area);
    const myb_dual k_u = myb_dual_scale(lit, 1.0 - p->specular);
    const myb_dual k_n = myb_dual_scale(myb_dual_mul(lit, normal), -1.0);
    for (int i = 0; i < 3; ++i) {
      a[i] = myb_dual_add(a[i], myb_dual_add(myb_dual_mul(k_u, u[i]), myb_dual_mul(k_n, n[i])));
    }
  }
  return MYB_OK;
}

myb_status myb_term_surface_radiation_jacobian(void* term, double t, const double x[3],
                                               const double v[3], double dadx[9],
                                               double dadv[9]) {
  const myb_surface_radiation_term* r = term;
  MYB_CHECK(r != NULL && r->env != NULL && r->object != NULL && r->attitude != NULL &&
            x != NULL && v != NULL && dadx != NULL && dadv != NULL &&
            r->occulter_count <= MYB_ENV_MAX_OCCULTERS &&
            (r->shape == NULL ? (r->plates != NULL || r->plate_count == 0)
                              : r->surface != NULL));
  const double mass = myb_internal_object_mass(r->object, t);
  if (!(mass > 0.0) || !(r->c > 0.0) || !isfinite(r->c)) {
    return MYB_E_DOMAIN;
  }
  double d[3];
  double distance = 0.0;
  double flux = 0.0;
  myb_status status = myb_internal_env_light(r->env, t, r->source, r->occulters,
                                             r->occulter_count, r->atmospheres, x, d, &distance,
                                             &flux);
  if (status != MYB_OK) {
    return status;
  }
  memset(dadx, 0, 9 * sizeof(double));
  memset(dadv, 0, 9 * sizeof(double));
  if (!(flux > 0.0)) {
    return MYB_OK; /* in the full shadow */
  }
  const double full =
      myb_radiation_flux(myb_internal_env_body(r->env, r->source)->luminosity, distance);
  if (flux < full) {
    /* a penumbra: the light changes over tens of kilometres */
    const double reach = myb_vec3_norm(x);
    const double speed = myb_vec3_norm(v);
    return myb_internal_state_differences(myb_term_surface_radiation, term, t, x, v,
                                          difference_step * (reach > 0.0 ? reach : distance),
                                          difference_step * (speed > 0.0 ? speed : 1.0), dadx,
                                          dadv);
  }
  /* the vector from the source, the direction of the light and its pressure over the mass */
  myb_dual from[3];
  for (int i = 0; i < 3; ++i) {
    from[i] = myb_dual_variable(d[i], i);
  }
  const myb_dual dn = myb_dual_norm(from);
  myb_dual u[3];
  for (int i = 0; i < 3; ++i) {
    u[i] = myb_dual_div(from[i], dn);
  }
  const myb_dual pressure = myb_dual_div(
      myb_dual_constant(flux * distance * distance / (r->c * mass)), myb_dual_mul(dn, dn));
  myb_dual m[9];
  status = attitude_dual(r->attitude, r->attitude_context, t, x, v, m);
  if (status != MYB_OK) {
    return status;
  }
  myb_dual u_body[3];
  to_body_dual(m, u, u_body);
  myb_dual a_body[3];
  if (r->shape != NULL) {
    const double direction[3] = {u_body[0].v, u_body[1].v, u_body[2].v};
    double force[3];
    double gradient[9];
    status = myb_internal_radiation_mesh_gradient(r->shape, r->surface, flux / r->c, direction,
                                                  r->shadows, force, gradient);
    if (status != MYB_OK) {
      return status;
    }
    /* the force follows the pressure in proportion and the direction through the gradient */
    for (int i = 0; i < 3; ++i) {
      a_body[i] = myb_dual_scale(pressure, force[i] / (mass * pressure.v));
      for (int j = 0; j < 3; ++j) {
        a_body[i] = myb_dual_add(
            a_body[i], myb_dual_scale(myb_dual_sub(u_body[j], myb_dual_constant(u_body[j].v)),
                                      gradient[3 * i + j] / mass));
      }
    }
  } else {
    status = radiation_plates_dual(pressure, r->plates, r->plate_count, u_body, a_body);
    if (status != MYB_OK) {
      return status;
    }
  }
  myb_dual a[3];
  from_body_dual(m, a_body, a);
  myb_dual_jacobian(a, dadx, dadv);
  return MYB_OK;
}

static myb_dual erf_dual(myb_dual a) {
  myb_dual r;
  r.v = myb_internal_erf(a.v);
  const double slope = 2.0 / sqrt_pi * myb_exp(-a.v * a.v);
  for (int i = 0; i < MYB_DUAL_VARIABLES; ++i) {
    r.d[i] = slope * a.d[i];
  }
  return r;
}

/* myb_force_gas_plates on such numbers. The lift is written as (g z + reemitted) times the
 * part of the normal across the flow — the same as the coefficient times its unit vector,
 * without the division that a plate facing the flow would not allow. */
static myb_status gas_plates_dual(myb_dual density, const myb_gas* gas, double mass,
                                  const myb_plate* plates, size_t count, const myb_dual v[3],
                                  myb_dual a[3]) {
  const myb_dual speed = myb_dual_norm(v);
  if (!valid_gas(gas) || !(speed.v > 0.0) || !isfinite(speed.v)) {
    return MYB_E_DOMAIN;
  }
  const double thermal2 = 2.0 * gas_constant * gas->temperature / gas->molar_mass;
  const myb_dual s = myb_dual_scale(speed, 1.0 / sqrt(thermal2));
  const myb_dual speed2 = myb_dual_mul(speed, speed);
  myb_dual inside = myb_dual_scale(
      myb_dual_div(myb_dual_constant(4.0 * gas_constant * gas->wall_temperature / gas->molar_mass),
                   speed2),
      0.5 * gas->accommodation);
  inside.v += 0.5 * (1.0 - gas->accommodation);
  const myb_dual leaving = myb_dual_sqrt(inside);
  const myb_dual g = myb_dual_div(myb_dual_constant(0.5), myb_dual_mul(s, s));
  myb_dual q = g;
  q.v += 1.0;
  myb_dual drag[3];
  for (int i = 0; i < 3; ++i) {
    drag[i] = myb_dual_scale(myb_dual_div(v[i], speed), -1.0);
  }
  myb_dual sum[3] = {myb_dual_constant(0.0), myb_dual_constant(0.0), myb_dual_constant(0.0)};
  for (size_t k = 0; k < count; ++k) {
    const myb_plate* p = &plates[k];
    const double n2 = myb_vec3_dot(p->normal, p->normal);
    if (!(p->area >= 0.0) || !isfinite(p->area) || !(fabs(n2 - 1.0) <= 1e-12)) {
      return MYB_E_DOMAIN;
    }
    const myb_dual n[3] = {myb_dual_constant(p->normal[0]), myb_dual_constant(p->normal[1]),
                           myb_dual_constant(p->normal[2])};
    const myb_dual gamma = myb_dual_scale(myb_dual_dot(drag, n), -1.0);
    const myb_dual gs = myb_dual_mul(gamma, s);
    const myb_dual pe =
        myb_dual_div(myb_dual_exp(myb_dual_scale(myb_dual_mul(gs, gs), -1.0)), s);
    myb_dual z = erf_dual(gs);
    z.v += 1.0;
    const myb_dual reemitted = myb_dual_scale(
        myb_dual_mul(leaving, myb_dual_add(myb_dual_scale(myb_dual_mul(gamma, z), sqrt_pi), pe)),
        0.5);
    const myb_dual cd =
        myb_dual_add(myb_dual_scale(pe, 1.0 / sqrt_pi),
                     myb_dual_mul(gamma, myb_dual_add(myb_dual_mul(q, z), reemitted)));
    const myb_dual lift = myb_dual_add(myb_dual_mul(g, z), reemitted);
    for (int i = 0; i < 3; ++i) {
      const myb_dual across = myb_dual_add(n[i], myb_dual_mul(gamma, drag[i]));
      sum[i] = myb_dual_add(
          sum[i], myb_dual_scale(myb_dual_sub(myb_dual_mul(cd, drag[i]),
                                              myb_dual_mul(lift, across)),
                                 p->area));
    }
  }
  const myb_dual scale = myb_dual_scale(myb_dual_mul(density, speed2), 0.5 / mass);
  for (int i = 0; i < 3; ++i) {
    a[i] = myb_dual_mul(scale, sum[i]);
  }
  return MYB_OK;
}

myb_status myb_term_surface_drag_jacobian(void* term, double t, const double x[3],
                                          const double v[3], double dadx[9], double dadv[9]) {
  const myb_surface_drag_term* d = term;
  MYB_CHECK(d != NULL && d->env != NULL && d->density != NULL && d->object != NULL &&
            d->attitude != NULL && (d->plates != NULL || d->plate_count == 0) && x != NULL &&
            v != NULL && dadx != NULL && dadv != NULL);
  const double mass = myb_internal_object_mass(d->object, t);
  if (!(mass > 0.0)) {
    return MYB_E_DOMAIN;
  }
  double rho = 0.0;
  double gradient[3];
  double flow[3];
  double dwdx[9];
  myb_status status =
      myb_internal_env_flow(d->env, t, d->body, d->density, d->density_context, d->wind,
                            d->wind_context, x, v, &rho, gradient, flow, dwdx);
  if (status != MYB_OK) {
    return status;
  }
  memset(dadx, 0, 9 * sizeof(double));
  memset(dadv, 0, 9 * sizeof(double));
  myb_dual m[9];
  status = attitude_dual(d->attitude, d->attitude_context, t, x, v, m);
  if (status != MYB_OK) {
    return status;
  }
  if (rho == 0.0) {
    return MYB_OK; /* above the gas */
  }
  myb_dual density = myb_dual_constant(rho);
  myb_dual w[3];
  for (int i = 0; i < 3; ++i) {
    density.d[i] = gradient[i];
    w[i] = myb_dual_variable(flow[i], 3 + i);
    for (int j = 0; j < 3; ++j) {
      w[i].d[j] = dwdx[3 * i + j];
    }
  }
  myb_dual w_body[3];
  to_body_dual(m, w, w_body);
  myb_dual a_body[3];
  status = gas_plates_dual(density, &d->gas, mass, d->plates, d->plate_count, w_body, a_body);
  if (status != MYB_OK) {
    return status;
  }
  myb_dual a[3];
  from_body_dual(m, a_body, a);
  if (d->length > 0.0) {
    if (!(d->molecule > 0.0) || !isfinite(d->continuum_cd) || !isfinite(d->continuum_area)) {
      return MYB_E_DOMAIN;
    }
    /* the weight of the free molecules: sin^2(pi (3 + lg Kn) / 8), Kn as 1 / density */
    const double knudsen = myb_gas_knudsen(rho, d->gas.molar_mass, d->molecule, d->length);
    myb_dual weight = myb_dual_constant(myb_gas_bridge(knudsen));
    if (knudsen > 1.0e-3 && knudsen < 10.0) {
      const double angle =
          3.141592653589793 * (3.0 + myb_log(knudsen) / 2.302585092994046) / 8.0;
      const double slope = 2.0 * myb_sin(angle) * myb_cos(angle) * 3.141592653589793 /
                           (8.0 * 2.302585092994046) * (-1.0 / rho);
      for (int j = 0; j < MYB_DUAL_VARIABLES; ++j) {
        weight.d[j] = slope * density.d[j];
      }
    }
    const myb_dual k = myb_dual_scale(myb_dual_mul(density, myb_dual_norm(w)),
                                      -0.5 * d->continuum_cd * d->continuum_area / mass);
    myb_dual rest = myb_dual_scale(weight, -1.0);
    rest.v += 1.0;
    for (int i = 0; i < 3; ++i) {
      a[i] = myb_dual_add(myb_dual_mul(weight, a[i]),
                          myb_dual_mul(rest, myb_dual_mul(k, w[i])));
    }
  }
  myb_dual_jacobian(a, dadx, dadv);
  return MYB_OK;
}

myb_status myb_term_recoil_jacobian(void* term, double t, const double x[3], const double v[3],
                                    double dadx[9], double dadv[9]) {
  const myb_recoil_term* r = term;
  MYB_CHECK(r != NULL && r->object != NULL && r->attitude != NULL && x != NULL && v != NULL &&
            dadx != NULL && dadv != NULL);
  const double mass = myb_internal_object_mass(r->object, t);
  const double n2 = myb_vec3_dot(r->direction, r->direction);
  if (!(mass > 0.0) || !(r->c > 0.0) || !isfinite(r->c) || !isfinite(r->power) ||
      !isfinite(r->efficiency) || !(fabs(n2 - 1.0) <= 1e-12) || !(r->decay >= 0.0) ||
      !isfinite(r->epoch)) {
    return MYB_E_DOMAIN;
  }
  myb_dual m[9];
  const myb_status status = attitude_dual(r->attitude, r->attitude_context, t, x, v, m);
  if (status != MYB_OK) {
    return status;
  }
  const double power = r->decay > 0.0 ? r->power * myb_exp(-(t - r->epoch) / r->decay) : r->power;
  double fixed[3] = {0.0, 0.0, 0.0};
  myb_force_emission(power, r->efficiency, mass, r->c, r->direction, fixed);
  const myb_dual a_body[3] = {myb_dual_constant(fixed[0]), myb_dual_constant(fixed[1]),
                              myb_dual_constant(fixed[2])};
  myb_dual a[3];
  from_body_dual(m, a_body, a);
  myb_dual_jacobian(a, dadx, dadv);
  return MYB_OK;
}
