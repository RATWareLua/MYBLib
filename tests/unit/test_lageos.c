/* SPDX-License-Identifier: BSD-3-Clause */
/* LAGEOS-1 over one week against the combined orbit of the International Laser Ranging Service
 * (ILRS, ilrsa.orb.lageos1.250607.v80.sp3: 2025-06-01 .. 06-07, SLRF2020; the excerpt
 * tests/data/sp3/ilrsa.lageos1.250607.excerpt.sp3 keeps every 3rd epoch, 6 min, 1680 positions;
 * the file of EDC/DGFI, ftp://edc.dgfi.tum.de/pub/slr/products/orbits/TRF/lageos1/2025/250607/,
 * SHA-256 of the .gz ba1d9933c870c44b4e0533f423c343340e09947bdbb6096f6410de86a8bb7ef5).
 *
 * The model is assembled from myb's pieces with every number from data:
 *  - Earth: GOCO06s to degree 30 (zero-tide, with its trends and annual terms), solid tides of
 *    IERS 2010 (steps 1-3), FES2004 ocean tides to degree 20, solid and ocean pole tides; the
 *    orientation of IAU 2006/2000A with EOP 20 C04 and its diurnal corrections;
 *  - Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn from DE440 (GM of header.440);
 *  - relativity of IERS 2010 eq. 10.12 (Schwarzschild, Lense-Thirring, de Sitter);
 *  - solar radiation pressure on a sphere (cr fitted) with the penumbra of the Earth and the
 *    Moon and the limb darkening of the Sun; Earth's albedo and infrared after Knocke et al.;
 *  - empirical accelerations as estimated by the ILRS analysis centres: a constant along-track
 *    term and once-per-revolution along- and cross-track terms.
 * The model is built twice: from the ready-made terms of myb/env.h in a force model, propagated
 * by myb/orbit.h — the one that is fitted — and by hand from the same functions (the function
 * `acceleration` below), kept as the reference the first must agree with at every state.
 * The initial state and the five empirical parameters are fitted to the positions by the
 * least squares of myb/estim.h: variational equations with the analytic Jacobians of the
 * field (its gradient with the tides), the bodies, relativity and sunlight; the gradients of
 * the Earth's radiation and of the empirical terms, 1e-9 of the field's, are left out. The
 * residual is the misfit of our dynamics to the ILRS orbit (itself fitted to laser ranges at
 * the cm level by models of the same kind).
 * Needs data/download (scripts/fetch-data.cmd: de440.bsp, header.440, eopc04.1962-now,
 * GOCO06s.gfc, fes2004_Cnm-Snm.dat); skipped otherwise. Measured contributions (RMS over the
 * week after the fit, the full model 2.0 cm): without the solid tides 21 cm, without the ocean
 * tides 6.3 cm, without relativity 2.06 cm (radial 5.3 mm against 2.6), without shadows 2.1 cm;
 * MYB_LAGEOS_OFF selects parts to leave out, MYB_LAGEOS_PROFILE times an evaluation. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mtest.h"
#include "myb/myb.h"
#include "num/vec3.h"

#define DOWNLOAD MYB_DATA_DIR "/download/"

enum { DEGREE = 30, OCEAN_DEGREE = 20, PARAMETERS = 12, STRIDE = 1, BODIES = 7 };

static const double c_light = 299792458.0;
static const double pi = 3.141592653589793;
static const double lageos_mass = 406.965;                         /* kg */
static const double lageos_area = 3.141592653589793 * 0.30 * 0.30; /* m^2 */
static const double area_over_mass = 3.141592653589793 * 0.30 * 0.30 / 406.965; /* LAGEOS-1 */
static const double solar_flux_1au = 1361.0;  /* W m^-2 (Kopp & Lean 2011) */
static const double au_m = 149597870700.0;
static const double earth_j = 9.8e8;           /* angular momentum per unit mass, m^2/s */
static const double sun_radius = 695700000.0;  /* IAU 2015 B3 */
static const double moon_radius = 1737400.0;
static const double earth_radius = 6378137.0;

/* albedo and infrared of the Earth (Knocke et al. 1988), seasonal phase from 1981-12-22 */
static const myb_planet_radiation knocke = {
    .radius = 6378137.0,
    .albedo = {0.34, 0.0, 0.10, 0.0, 0.29},
    .emissivity = {0.68, 0.0, -0.07, 0.0, -0.18},
    .omega = 2.0 * 3.141592653589793 / (365.25 * 86400.0)};
static const double knocke_epoch = 568915200.0; /* s from 1981-12-22 to J2000.0 */

typedef struct model {
  myb_timesys* sys;
  myb_scale tt;
  myb_scale tdb;
  myb_cip* cip;
  myb_eop* eop;
  myb_spk* spk;
  myb_gravity_field* field;
  myb_gravity_work* work;
  myb_tides* tides;
  myb_time t0; /* TT */
  int body[BODIES];
  double gm_body[BODIES];
  double gm_earth;
  double cr;
  double empirical[5]; /* along: constant, cos u, sin u; cross: cos u, sin u */
  /* what depends on t only, kept for the next call at the same t (IAS15 evaluates every
   * substep of a step several times while it iterates) */
  int off; /* MYB_LAGEOS_OFF: 1 solid tides, 2 ocean, 4 pole, 8 relativity, 16 albedo,
              32 planets, 128 shadows (for the table of contributions) */
  int cached;
  double cache_t;
  myb_rotation rot;
  myb_rotation inv;
  myb_eop_values eop_values;
  /* TDB - TT (s) every hour from t = -1 h: the series of 787 terms costs more than all the
   * rest of an evaluation; 4-point interpolation of the table is exact to 1e-13 s */
  double tdb_tt[200];
  double s[BODIES][6];
} model;

static int read_constant(const char* header, const char* name, double* value) {
  const char* g1040 = strstr(header, "GROUP   1040");
  const char* g1041 = strstr(header, "GROUP   1041");
  if (g1040 == NULL || g1041 == NULL) {
    return 0;
  }
  const size_t len = strlen(name);
  const char* p = g1040 + 12;
  double count = 0.0;
  size_t used = 0;
  /* bound: white space */
  while (*p == ' ' || *p == '\n' || *p == '\r') {
    ++p;
  }
  if (myb_parse_double(p, strlen(p), &count, &used) != MYB_OK) {
    return 0;
  }
  p += used;
  int index = -1;
  for (int k = 0; k < (int)count; ++k) {
    /* bound: white space */
    while (*p == ' ' || *p == '\n' || *p == '\r') {
      ++p;
    }
    size_t n = 0;
    /* bound: a name */
    while (p[n] != ' ' && p[n] != '\n' && p[n] != '\r' && p[n] != '\0') {
      ++n;
    }
    if (n == len && memcmp(p, name, len) == 0) {
      index = k;
    }
    p += n;
  }
  if (index < 0) {
    return 0;
  }
  p = g1041 + 12;
  /* bound: white space */
  while (*p == ' ' || *p == '\n' || *p == '\r') {
    ++p;
  }
  if (myb_parse_double(p, strlen(p), &count, &used) != MYB_OK) {
    return 0;
  }
  p += used;
  for (int k = 0; k <= index; ++k) {
    /* bound: white space */
    while (*p == ' ' || *p == '\n' || *p == '\r') {
      ++p;
    }
    if (myb_parse_double(p, strlen(p), value, &used) != MYB_OK) {
      return 0;
    }
    p += used;
  }
  return 1;
}

static int read_file(const char* path, char* buffer, size_t size) {
  FILE* f = fopen(path, "rb");
  if (f == NULL) {
    return 0;
  }
  const size_t n = fread(buffer, 1, size - 1, f);
  fclose(f);
  buffer[n] = '\0';
  return n > 0;
}

static void scale_add3(double* a, double k, const double* b) {
  a[0] += k * b[0];
  a[1] += k * b[1];
  a[2] += k * b[2];
}

/* Earth orientation at TT instant t: celestial -> terrestrial, EOP values, UT1. */
static myb_status orientation(const model* m, myb_time t, myb_rotation* rot,
                              myb_eop_values* eop, myb_time* ut1) {
  const myb_status st = myb_eop_at(m->eop, m->sys, t, m->tt, eop, ut1);
  if (st != MYB_OK) {
    return st;
  }
  myb_cip_terrestrial(m->cip, t, *ut1, eop, rot);
  return MYB_OK;
}

static myb_status acceleration(void* context, double t, const double* x, const double* v,
                               double* a) {
  model* m = context;
  const myb_time tt = myb_time_add(m->t0, t);
  myb_status st = MYB_OK;
  if (!m->cached || m->cache_t != t) {
    m->cached = 0;
    const double u = t / 3600.0 + 1.0;
    const int k0 = (int)floor(u) - 1;
    if (k0 < 0 || k0 + 3 >= 200) {
      return MYB_E_RANGE;
    }
    const double w = u - (double)k0; /* in [1, 2) */
    const double* f = m->tdb_tt + k0;
    const double dtdb = -f[0] * (w - 1.0) * (w - 2.0) * (w - 3.0) / 6.0 +
                        f[1] * w * (w - 2.0) * (w - 3.0) / 2.0 -
                        f[2] * w * (w - 1.0) * (w - 3.0) / 2.0 +
                        f[3] * w * (w - 1.0) * (w - 2.0) / 6.0;
    const myb_time tdb = myb_time_add(tt, dtdb);
    myb_time ut1;
    st = orientation(m, tt, &m->rot, &m->eop_values, &ut1);
    if (st != MYB_OK) {
      return st;
    }
    myb_rotation_invert(&m->rot, &m->inv);
    /* bodies, geocentric (GCRS axes = ICRF axes) */
    for (int b = 0; b < BODIES; ++b) {
      st = myb_spk_state(m->spk, m->body[b], 399, tdb, m->s[b], NULL);
      if (st != MYB_OK) {
        return st;
      }
    }
    /* the field with its time-variable terms and tides */
    myb_gravity_field_set_time(m->field, tt);
    myb_tide_raiser raisers[2];
    raisers[0].gm = m->gm_body[0];
    myb_rotation_apply_position(&m->rot, m->s[0], raisers[0].r);
    raisers[1].gm = m->gm_body[1];
    myb_rotation_apply_position(&m->rot, m->s[1], raisers[1].r);
    const myb_tides_input tin = {.raisers = raisers, .count = 2, .t = ut1, .xp = m->eop_values.xp,
                                 .yp = m->eop_values.yp,
                                 .flags = MYB_TIDES_ALL &
                                          ~((m->off & 1 ? MYB_TIDES_SOLID | MYB_TIDES_FREQUENCY : 0) |
                                            (m->off & 256 ? MYB_TIDES_FREQUENCY : 0) |
                                            (m->off & 2 ? MYB_TIDES_OCEAN : 0) |
                                            (m->off & 4 ? MYB_TIDES_POLE_SOLID |
                                                              MYB_TIDES_POLE_OCEAN
                                                        : 0))};
    st = myb_tides_apply(m->tides, m->field, &tin);
    if (st != MYB_OK) {
      return st;
    }
    m->cached = 1;
    m->cache_t = t;
  }
  const myb_rotation rot = m->rot;
  const myb_rotation inv = m->inv;
  const double (*s)[6] = (const double (*)[6])m->s;
  const double* sun = s[0];
  const double* moon = s[1];
  double r_itrs[3];
  myb_rotation_apply_position(&rot, x, r_itrs);
  double g_itrs[3];
  st = myb_gravity_field_eval(m->field, m->work, r_itrs, 0, -1, NULL, g_itrs);
  if (st != MYB_OK) {
    return st;
  }
  myb_rotation_apply_position(&inv, g_itrs, a);
  /* third bodies */
  for (int b = 0; b < ((m->off & 32) ? 2 : BODIES); ++b) {
    myb_force_third_body(m->gm_body[b], s[b], x, a);
  }
  /* relativity */
  const myb_ppn ppn = {c_light, 1.0, 1.0};
  double tmp[3];
  const double rel_on = (m->off & 8) ? 0.0 : 1.0;
  st = myb_rel_schwarzschild(&ppn, m->gm_earth, x, v, tmp);
  if (st != MYB_OK) {
    return st;
  }
  scale_add3(a, rel_on, tmp);
  const double j_itrs[3] = {0.0, 0.0, earth_j};
  double j[3];
  myb_rotation_apply_position(&inv, j_itrs, j);
  st = myb_rel_lense_thirring(&ppn, m->gm_earth, j, x, v, tmp);
  if (st != MYB_OK) {
    return st;
  }
  scale_add3(a, rel_on, tmp);
  const double R[3] = {-sun[0], -sun[1], -sun[2]};
  const double V[3] = {-sun[3], -sun[4], -sun[5]};
  st = myb_rel_de_sitter(&ppn, m->gm_body[0], R, V, v, tmp);
  if (st != MYB_OK) {
    return st;
  }
  scale_add3(a, rel_on, tmp);
  /* solar radiation pressure with shadows */
  const double d[3] = {x[0] - sun[0], x[1] - sun[1], x[2] - sun[2]};
  const double dn = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  const double u[3] = {d[0] / dn, d[1] / dn, d[2] / dn};
  const double source[3] = {-d[0], -d[1], -d[2]};
  const double earth[3] = {-x[0], -x[1], -x[2]};
  const double luna[3] = {moon[0] - x[0], moon[1] - x[1], moon[2] - x[2]};
  const double limb[2] = {0.6, 0.0};
  const double nu = (m->off & 128) ? 1.0 : myb_shadow_fraction(source, sun_radius, earth, earth_radius, limb) *
                    myb_shadow_fraction(source, sun_radius, luna, moon_radius, limb);
  const double pressure = solar_flux_1au * (au_m / dn) * (au_m / dn) / c_light;
  myb_force_radiation_sphere(nu * pressure, m->cr, area_over_mass, u, a);
  /* albedo and infrared of the Earth */
  const double rs = sqrt(sun[0] * sun[0] + sun[1] * sun[1] + sun[2] * sun[2]);
  const double phase = knocke.omega * (myb_time_seconds(tt) + knocke_epoch);
  double sun_itrs[3];
  myb_rotation_apply_position(&rot, sun, sun_itrs);
  double albedo[3] = {0.0, 0.0, 0.0};
  st = myb_force_planet_radiation(&knocke, r_itrs, sun_itrs,
                                  solar_flux_1au * (au_m / rs) * (au_m / rs), phase,
                                  m->cr * area_over_mass, 10, albedo);
  if (st != MYB_OK) {
    return st;
  }
  myb_rotation_apply_position(&inv, albedo, tmp);
  scale_add3(a, (m->off & 16) ? 0.0 : 1.0, tmp);
  /* empirical: argument of latitude u from the node */
  double h[3];
  myb_vec3_cross(x, v, h);
  const double hn = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
  const double nodal[3] = {-h[1], h[0], 0.0}; /* z x h */
  const double nn = hypot(nodal[0], nodal[1]);
  const double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  double hxn[3];
  const double hu[3] = {h[0] / hn, h[1] / hn, h[2] / hn};
  const double nu_[3] = {nodal[0] / nn, nodal[1] / nn, 0.0};
  myb_vec3_cross(hu, nu_, hxn);
  const double cu = (nu_[0] * x[0] + nu_[1] * x[1]) / r;
  const double su = (hxn[0] * x[0] + hxn[1] * x[1] + hxn[2] * x[2]) / r;
  const double* e = m->empirical;
  const double rtn[3] = {0.0, e[0] + e[1] * cu + e[2] * su, e[3] * cu + e[4] * su};
  return myb_force_empirical_rtn(x, v, rtn, a);
}

/* ---- The same model from the ready-made terms of myb/env.h ---------------------------------- */

typedef struct assembly {
  myb_env* env;
  myb_force_model* forces;
  myb_orbit* orbit;
  myb_gravity_field* field; /* its own: the term changes the coefficients */
  myb_gravity_work* work;
  myb_cip_orientation earth_orientation;
  myb_spk_body ephemeris[BODIES];
  myb_object object;
  myb_field_term gravity;
  myb_third_body_term third[BODIES];
  myb_relativity_term relativity;
  myb_radiation_term pressure;
  myb_eclipse_event edges[4]; /* of the shadows of the Earth and the Moon */
  myb_planet_radiation_term albedo;
  myb_empirical_term empirical;
} assembly;

static void disassemble(assembly* s) {
  myb_orbit_destroy(s->orbit);
  myb_force_model_destroy(s->forces);
  myb_gravity_work_destroy(s->work);
  myb_gravity_field_destroy(s->field);
  myb_env_destroy(s->env);
}

/* Builds the environment, the terms and the propagator from the data of m (shared, read
 * only); m->off leaves out the same parts as in the reference. */
static int assemble(const model* m, assembly* s, int estimate_cr) {
  memset(s, 0, sizeof *s);
  const myb_env_options options = {.epoch = m->t0, .timesys = m->sys, .scale = "TT"};
  int ok = myb_env_create(&options, NULL, &s->env) == MYB_OK &&
           myb_gravity_field_clone(m->field, NULL, &s->field) == MYB_OK &&
           myb_gravity_work_create(DEGREE, NULL, &s->work) == MYB_OK &&
           myb_force_model_create(NULL, &s->forces) == MYB_OK;
  if (!ok) {
    return 0;
  }
  /* bodies: the Earth at the origin, the others from DE440 on TDB */
  s->earth_orientation = (myb_cip_orientation){m->cip, m->eop, m->sys, m->tt};
  size_t earth = 0;
  const myb_env_body b_earth = {.id = 399,
                                .gm = m->gm_earth,
                                .orientation = myb_orientation_cip,
                                .orientation_context = &s->earth_orientation,
                                .radius = earth_radius,
                                .spin = {0.0, 0.0, earth_j}};
  ok = myb_env_add_body(s->env, &b_earth, &earth) == MYB_OK;
  size_t index[BODIES] = {0};
  for (int b = 0; b < BODIES && ok; ++b) {
    s->ephemeris[b] =
        (myb_spk_body){.spk = m->spk, .target = m->body[b], .center = 399, .frame = 1};
    myb_env_body body = {.id = m->body[b],
                         .gm = m->gm_body[b],
                         .state = myb_spk_body_state,
                         .state_context = &s->ephemeris[b],
                         .state_scale = "TDB"};
    if (b == 0) {
      body.radius = sun_radius;
      body.luminosity = 4.0 * pi * au_m * au_m * solar_flux_1au;
      body.limb[0] = 0.6;
    } else if (b == 1) {
      body.radius = moon_radius;
    }
    ok = myb_env_add_body(s->env, &body, &index[b]) == MYB_OK;
  }
  if (!ok) {
    return 0;
  }
  const int off = m->off;
  s->object = (myb_object){.mass = lageos_mass, .area = lageos_area, .cr = 1.13};
  /* the field with its trends, annual terms and tides */
  s->gravity = (myb_field_term){
      .env = s->env,
      .body = earth,
      .field = s->field,
      .work = s->work,
      .n_min = 0,
      .n_max = -1,
      .update = 1,
      .tides = m->tides,
      .tide_flags = MYB_TIDES_ALL & ~((off & 1 ? MYB_TIDES_SOLID | MYB_TIDES_FREQUENCY : 0) |
                                      (off & 256 ? MYB_TIDES_FREQUENCY : 0) |
                                      (off & 2 ? MYB_TIDES_OCEAN : 0) |
                                      (off & 4 ? MYB_TIDES_POLE_SOLID | MYB_TIDES_POLE_OCEAN : 0)),
      .permanent = MYB_PERMANENT_TIDE_AUTO,
      .raisers = {index[0], index[1]},
      .raiser_count = 2};
  ok = myb_force_model_add(s->forces, "earth", myb_term_field, myb_term_field_jacobian,
                           &s->gravity, NULL) == MYB_OK;
  static const char* const names[BODIES] = {"sun",  "moon",    "mercury", "venus",
                                            "mars", "jupiter", "saturn"};
  for (int b = 0; b < BODIES && ok; ++b) {
    s->third[b] = (myb_third_body_term){s->env, index[b]};
    size_t term = 0;
    ok = myb_force_model_add(s->forces, names[b], myb_term_third_body,
                             myb_term_third_body_jacobian, &s->third[b], &term) == MYB_OK;
    if (ok && b >= 2 && (off & 32)) {
      myb_force_model_set_enabled(s->forces, term, 0);
    }
  }
  /* relativity of IERS 2010 eq. 10.12 */
  s->relativity = (myb_relativity_term){
      .env = s->env, .ppn = {c_light, 1.0, 1.0}, .central = earth, .external = index[0]};
  const myb_term_function relativistic[3] = {myb_term_schwarzschild, myb_term_lense_thirring,
                                             myb_term_de_sitter};
  const myb_term_jacobian relativistic_jacobian[3] = {myb_term_schwarzschild_jacobian,
                                                      myb_term_lense_thirring_jacobian,
                                                      myb_term_de_sitter_jacobian};
  static const char* const rel_names[3] = {"schwarzschild", "lense-thirring", "de sitter"};
  for (int k = 0; k < 3 && ok; ++k) {
    size_t term = 0;
    ok = myb_force_model_add(s->forces, rel_names[k], relativistic[k], relativistic_jacobian[k],
                             &s->relativity, &term) == MYB_OK;
    if (ok && (off & 8)) {
      myb_force_model_set_enabled(s->forces, term, 0);
    }
  }
  /* sunlight with the shadows of the Earth and the Moon; the Earth's own radiation */
  s->pressure = (myb_radiation_term){.env = s->env,
                                     .source = index[0],
                                     .object = &s->object,
                                     .c = c_light,
                                     .occulters = {earth, index[1]},
                                     .occulter_count = (off & 128) ? 0 : 2};
  s->albedo = (myb_planet_radiation_term){
      .env = s->env,
      .planet = earth,
      .source = index[0],
      .model = &knocke,
      .object = &s->object,
      .phase0 = knocke.omega * (myb_time_seconds(m->t0) + knocke_epoch),
      .rings = 10};
  size_t term = 0;
  ok = ok &&
       myb_force_model_add(s->forces, "sunlight", myb_term_radiation,
                           myb_term_radiation_jacobian, &s->pressure, NULL) == MYB_OK &&
       myb_force_model_add(s->forces, "earthlight", myb_term_planet_radiation,
                           myb_term_no_jacobian, &s->albedo, &term) == MYB_OK;
  if (ok && (off & 16)) {
    myb_force_model_set_enabled(s->forces, term, 0);
  }
  /* the edges of the shadows end the steps of the integrator (MYB_LAGEOS_NO_BREAKS: not) */
  if (ok && getenv("MYB_LAGEOS_NO_BREAKS") == NULL) {
    ok = myb_term_radiation_breaks(&s->pressure, s->forces, s->edges) == MYB_OK;
  }
  size_t empirical = 0;
  ok = ok && myb_force_model_add(s->forces, "empirical", myb_term_empirical,
                                 myb_term_no_jacobian, &s->empirical, &empirical) == MYB_OK;
  /* the parameters of the fit, in the order of p[6 ..]: cr only on request (it enters two
   * terms), the five empirical ones bound to their term */
  if (estimate_cr) {
    ok = ok && myb_force_model_add_parameter(s->forces, "cr", &s->object.cr, 0.01, NULL) == MYB_OK;
  }
  double* const fitted[5] = {&s->empirical.constant[1], &s->empirical.cosine[1],
                             &s->empirical.sine[1], &s->empirical.cosine[2],
                             &s->empirical.sine[2]};
  static const char* const fitted_names[5] = {"along", "along cos", "along sin", "cross cos",
                                              "cross sin"};
  for (int k = 0; k < 5 && ok; ++k) {
    size_t parameter = 0;
    ok = myb_force_model_add_parameter(s->forces, fitted_names[k], fitted[k], 1e-11,
                                       &parameter) == MYB_OK &&
         myb_force_model_bind_parameter(s->forces, parameter, empirical) == MYB_OK;
  }
  const myb_orbit_options orbit_options = {.epsilon = 1e-10};
  return ok && myb_orbit_create(s->forces, &orbit_options, NULL, &s->orbit) == MYB_OK;
}

/* The parameters of the fit after the state: cr, then along-track constant, cos u, sin u and
 * cross-track cos u, sin u. */
static void set_parameters(assembly* s, const double* p) {
  s->object.cr = p[6];
  s->empirical.constant[1] = p[7];
  s->empirical.cosine[1] = p[8];
  s->empirical.sine[1] = p[9];
  s->empirical.cosine[2] = p[10];
  s->empirical.sine[2] = p[11];
}

/* MYB_LAGEOS_PROFILE: time of the pieces of one evaluation of the acceleration. */
static void profile(model* m, const double* x, const double* v) {
  enum { REPS = 300 };
  double t_rot = 0.0;
  double t_tdb = 0.0;
  double t_spk = 0.0;
  double t_set = 0.0;
  double t_tides = 0.0;
  double t_eval = 0.0;
  double t_alb = 0.0;
  double t_all = 0.0;
  myb_rotation rot;
  myb_eop_values eop;
  myb_time ut1;
  myb_time tdb;
  double s[6];
  double a[3];
  for (int k = 0; k < REPS; ++k) {
    const myb_time tt = myb_time_add(m->t0, 100.0 * k);
    clock_t c0 = clock();
    (void)orientation(m, tt, &rot, &eop, &ut1);
    clock_t c1 = clock();
    (void)myb_timesys_convert(m->sys, tt, m->tt, m->tdb, &tdb);
    clock_t c2 = clock();
    for (int b = 0; b < BODIES; ++b) {
      (void)myb_spk_state(m->spk, m->body[b], 399, tdb, s, NULL);
    }
    clock_t c3 = clock();
    myb_gravity_field_set_time(m->field, tt);
    clock_t c4 = clock();
    myb_tide_raiser raisers[2] = {{m->gm_body[0], {1.4e11, 2e10, 1e10}},
                                  {m->gm_body[1], {3e8, 1e8, 1e8}}};
    const myb_tides_input tin = {.raisers = raisers, .count = 2, .t = ut1, .xp = eop.xp,
                                 .yp = eop.yp, .flags = MYB_TIDES_ALL};
    (void)myb_tides_apply(m->tides, m->field, &tin);
    clock_t c5 = clock();
    (void)myb_gravity_field_eval(m->field, m->work, x, 0, -1, NULL, a);
    clock_t c6 = clock();
    (void)myb_force_planet_radiation(&knocke, x, (const double[3]){1.4e11, 2e10, 1e10}, 1361.0,
                                     0.3, 1e-3, 10, a);
    clock_t c7 = clock();
    (void)acceleration(m, 100.0 * k, x, v, a);
    clock_t c8 = clock();
    t_rot += (double)(c1 - c0);
    t_tdb += (double)(c2 - c1);
    t_spk += (double)(c3 - c2);
    t_set += (double)(c4 - c3);
    t_tides += (double)(c5 - c4);
    t_eval += (double)(c6 - c5);
    t_alb += (double)(c7 - c6);
    t_all += (double)(c8 - c7);
  }
  const double us = 1e6 / (double)CLOCKS_PER_SEC / (double)REPS;
  printf("  profile (us per call): orientation %.1f, TT->TDB %.1f, 7 bodies %.1f, set_time %.1f, "
         "tides %.1f, field %.1f, albedo %.1f; whole acceleration %.1f\n",
         t_rot * us, t_tdb * us, t_spk * us, t_set * us, t_tides * us, t_eval * us, t_alb * us,
         t_all * us);
}

typedef struct observations {
  size_t count;
  double t[2000];       /* TT seconds from t0 */
  double r[2000][3];    /* GCRS positions, m */
  double rtn[2000][9];  /* rows R, T, N of the orbital frame at the observation */
} observations;

/* The force model of the ready-made terms against the model assembled by hand, at states
 * along the orbit from p: the largest difference relative to the acceleration. */
static myb_status compare_models(model* m, assembly* s, const observations* obs, const double* p,
                                 double* worst, int* states) {
  /* every parameter in play */
  double q[PARAMETERS];
  memcpy(q, p, sizeof q);
  const double trial[5] = {1e-11, -2e-11, 3e-11, 4e-11, -5e-11};
  memcpy(q + 7, trial, sizeof trial);
  set_parameters(s, q);
  m->cr = q[6];
  memcpy(m->empirical, q + 7, 5 * sizeof(double));
  *worst = 0.0;
  *states = 0;
  myb_status st = myb_orbit_start(s->orbit, 0.0, q, obs->t[obs->count - 1]);
  for (size_t k = 0; k < obs->count && st == MYB_OK; k += 40) {
    double y[6];
    double a[3];
    double ref[3];
    st = myb_orbit_state(s->orbit, obs->t[k], y);
    if (st == MYB_OK) {
      st = myb_force_model_acceleration(s->forces, obs->t[k], y, y + 3, a);
    }
    if (st == MYB_OK) {
      st = acceleration(m, obs->t[k], y, y + 3, ref);
    }
    if (st == MYB_OK) {
      const double d[3] = {a[0] - ref[0], a[1] - ref[1], a[2] - ref[2]};
      const double e = myb_vec3_norm(d) / myb_vec3_norm(ref);
      *worst = e > *worst ? e : *worst;
      ++*states;
    }
  }
  return st;
}

static void statistics(const observations* obs, const double* res, double* rms, double rtn[3]) {
  double s = 0.0;
  double c[3] = {0.0, 0.0, 0.0};
  for (size_t k = 0; k < obs->count; ++k) {
    for (int i = 0; i < 3; ++i) {
      const double* row = obs->rtn[k] + 3 * i;
      const double q = row[0] * res[3 * k] + row[1] * res[3 * k + 1] + row[2] * res[3 * k + 2];
      c[i] += q * q;
      s += res[3 * k + (size_t)i] * res[3 * k + (size_t)i];
    }
  }
  *rms = sqrt(s / (double)obs->count);
  for (int i = 0; i < 3; ++i) {
    rtn[i] = sqrt(c[i] / (double)obs->count);
  }
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  static char header[32768];
  static model m;
  static observations obs;
  myb_diag diag = {0};
  myb_sp3* sp3 = NULL;
  if (!read_file(DOWNLOAD "header.440", header, sizeof header) ||
      myb_sp3_load_file(MYB_TEST_DATA_DIR "/sp3/ilrsa.lageos1.250607.excerpt.sp3", NULL, &sp3, &diag) !=
          MYB_OK) {
    printf("header.440 not downloaded (or the SP3 excerpt missing): skipped\n");
    return MTEST_SKIP;
  }
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const char* kernels[] = {MYB_DATA_DIR "/time/time-scales.tk",
                           MYB_DATA_DIR "/time/tdb-tt-fb1990.tk", MYB_DATA_DIR "/time/utc-tai.tk",
                           MYB_DATA_DIR "/frames/iau2006-cip.tk",
                           MYB_DATA_DIR "/frames/eop-tides.tk",
                           MYB_DATA_DIR "/grav/iers2010-tides.tk"};
  for (size_t i = 0; i < sizeof kernels / sizeof kernels[0]; ++i) {
    MT_CHECK(myb_kpool_load_file(pool, kernels[i], &diag) == MYB_OK);
  }
  MT_CHECK(myb_timesys_create(pool, NULL, &m.sys, &diag) == MYB_OK);
  MT_CHECK(myb_timesys_scale(m.sys, "TT", &m.tt) == MYB_OK);
  MT_CHECK(myb_timesys_scale(m.sys, "TDB", &m.tdb) == MYB_OK);
  MT_CHECK(myb_cip_create(pool, "IAU2006A", NULL, &m.cip, &diag) == MYB_OK);
  MT_CHECK(myb_eop_create(NULL, &m.eop) == MYB_OK);
  MT_CHECK(myb_eop_set_tides(m.eop, pool, "IERS10_EOPT", &diag) == MYB_OK);
  MT_CHECK(myb_tides_create(pool, "IERS10_TIDE", NULL, &m.tides, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  MT_CHECK(myb_spk_create(NULL, &m.spk) == MYB_OK);
  const int ok =
      myb_eop_load_c04_file(m.eop, DOWNLOAD "eopc04.1962-now", &diag) == MYB_OK &&
      myb_spk_load_file(m.spk, DOWNLOAD "de440.bsp", &diag) == MYB_OK &&
      myb_gravity_field_load_icgem_file(DOWNLOAD "GOCO06s.gfc", DEGREE, NULL, &m.field, &diag) ==
          MYB_OK &&
      myb_tides_load_ocean_file(m.tides, DOWNLOAD "fes2004_Cnm-Snm.dat", OCEAN_DEGREE,
                                &diag) == MYB_OK;
  if (!ok) {
    printf("data/download incomplete (%s): skipped\n", diag.message);
    myb_sp3_destroy(sp3);
    return MTEST_SKIP;
  }
  MT_CHECK(myb_gravity_work_create(DEGREE, NULL, &m.work) == MYB_OK);
  m.gm_earth = myb_gravity_field_gm(m.field);
  /* GM of the perturbing bodies from DE440 (au^3/day^2) */
  double au = 0.0;
  double emrat = 0.0;
  MT_CHECK(read_constant(header, "AU", &au) && read_constant(header, "EMRAT", &emrat));
  const double unit = (au * 1e3) * (au * 1e3) * (au * 1e3) / (86400.0 * 86400.0);
  const char* names[BODIES] = {"GMS", "GMB", "GM1", "GM2", "GM4", "GM5", "GM6"};
  const int ids[BODIES] = {10, 301, 1, 2, 4, 5, 6};
  for (int b = 0; b < BODIES; ++b) {
    double g = 0.0;
    MT_CHECK(read_constant(header, names[b], &g));
    m.gm_body[b] = g * unit;
    m.body[b] = ids[b];
  }
  m.gm_body[1] /= 1.0 + emrat; /* Moon from the Earth-Moon barycentre */
  /* observations: every STRIDE-th epoch, ITRS (SLRF2020) -> GCRS */
  MT_CHECK(strcmp(myb_sp3_time_system(sp3), "UTC") == 0);
  const size_t epochs = myb_sp3_epoch_count(sp3);
  double p[PARAMETERS] = {0.0};
  for (size_t e = 0; e < epochs && obs.count < 2000; e += STRIDE) {
    myb_calendar label;
    myb_sp3_epoch(sp3, e, &label);
    myb_time t;
    MT_CHECK(myb_timesys_from_utc(m.sys, &label, m.tt, &t) == MYB_OK);
    if (e == 0) {
      m.t0 = t;
    }
    double itrs[6];
    const int flags = myb_sp3_state(sp3, 0, e, itrs, itrs + 3);
    if (!(flags & 1)) {
      continue;
    }
    myb_rotation rot;
    myb_eop_values eop;
    myb_time ut1;
    MT_CHECK(orientation(&m, t, &rot, &eop, &ut1) == MYB_OK);
    myb_rotation inv;
    myb_rotation_invert(&rot, &inv);
    double gcrs[6];
    myb_rotation_apply_state(&inv, itrs, gcrs);
    const size_t k = obs.count++;
    obs.t[k] = myb_dd_to_double(myb_time_diff(t, m.t0));
    memcpy(obs.r[k], gcrs, 3 * sizeof(double));
    double rot_rtn[9];
    double n[3];
    myb_vec3_cross(gcrs, gcrs + 3, n);
    const double rr = myb_vec3_norm(gcrs);
    const double nn = myb_vec3_norm(n);
    for (int i = 0; i < 3; ++i) {
      rot_rtn[i] = gcrs[i] / rr;
      rot_rtn[6 + i] = n[i] / nn;
    }
    myb_vec3_cross(rot_rtn + 6, rot_rtn, rot_rtn + 3);
    memcpy(obs.rtn[k], rot_rtn, sizeof rot_rtn);
    if (e == 0) {
      MT_CHECK(flags == 3);
      memcpy(p, gcrs, 6 * sizeof(double));
    }
  }
  myb_sp3_destroy(sp3);
  for (int k = 0; k < 200; ++k) {
    const myb_time tt = myb_time_add(m.t0, 3600.0 * (k - 1));
    myb_time tdb;
    MT_CHECK(myb_timesys_convert(m.sys, tt, m.tt, m.tdb, &tdb) == MYB_OK);
    m.tdb_tt[k] = myb_dd_to_double(myb_time_diff(tdb, tt));
  }
  printf("LAGEOS-1, ILRS combined orbit 2025-06-01 .. 06-07: %zu positions (every %d min)\n",
         obs.count, 6 * STRIDE);
  p[6] = 1.13;
  const int estimate_cr = getenv("MYB_LAGEOS_CR") != NULL;
  if (getenv("MYB_LAGEOS_OFF") != NULL) {
    m.off = atoi(getenv("MYB_LAGEOS_OFF"));
    printf("  model without the parts %d\n", m.off);
  }
  if (getenv("MYB_LAGEOS_PROFILE") != NULL) {
    profile(&m, p, p + 3);
  }
  static double res[6000];
  static myb_observation fit_obs[2000];
  for (size_t k = 0; k < obs.count; ++k) {
    fit_obs[k] = (myb_observation){.t = obs.t[k],
                                   .measure = myb_measure_position,
                                   .count = 3,
                                   .value = {obs.r[k][0], obs.r[k][1], obs.r[k][2]},
                                   .sigma = {1.0, 1.0, 1.0}};
  }
  static assembly set;
  MT_CHECK(assemble(&m, &set, estimate_cr));
  double model_difference = 1.0;
  int model_states = 0;
  MT_CHECK(compare_models(&m, &set, &obs, p, &model_difference, &model_states) == MYB_OK);
  printf("  ready-made terms against the model assembled by hand: %.1e of the acceleration at "
         "%d states\n",
         model_difference, model_states);
  MT_CHECK(model_states >= 40);
  MT_CHECK(model_difference < 1e-13);
  printf("  the week: %lld steps of the integrator, %lld edges of shadows met\n",
         (long long)myb_orbit_steps(set.orbit), (long long)myb_orbit_breaks(set.orbit));
  set_parameters(&set, p);
  /* the residuals of the first guess: every unknown held */
  const unsigned char held[PARAMETERS] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  const double epsilon = getenv("MYB_LAGEOS_EPS") != NULL ? atof(getenv("MYB_LAGEOS_EPS")) : 0.0;
  myb_fit_options options = {.fixed = held, .epsilon = epsilon};
  myb_fit_result result;
  double y0[6];
  memcpy(y0, p, sizeof y0);
  MT_CHECK(myb_fit_orbit(set.forces, 0.0, y0, fit_obs, obs.count, &options, NULL, &result, NULL,
                         res) == MYB_OK);
  double rms = 0.0;
  double rtn[3];
  statistics(&obs, res, &rms, rtn);
  printf("  from the ILRS state at the first epoch, cr 1.13, no empirical: RMS %.4f m "
         "(R %.4f, T %.4f, N %.4f)\n",
         rms, rtn[0], rtn[1], rtn[2]);
  /* with the edges of the shadows declared this number does not depend on the precision of
   * the integrator (25.8887 m at epsilon 1e-10, 1e-11 and 1e-12); with steps over the edges
   * it is 25.8821 m at 1e-10 and 25.8882 m at 1e-12. Without the secondary waves of the
   * ocean tides (MYB_TIDES_ADMITTANCE) it was 25.9587 m, and the fit below 2.05 cm in place
   * of 1.57 cm */
  if (m.off == 0 && !estimate_cr && getenv("MYB_LAGEOS_NO_BREAKS") == NULL) {
    MT_CHECK(fabs(rms - 25.8887) < 5.0e-4);
  }
  /* cr is held at 1.13 as by the ILRS analysis centres (it is not a parameter of the model
   * unless MYB_LAGEOS_CR asks): over a week a force of fixed inertial direction is nearly a
   * once-per-revolution signal in R, T, N, so cr and the once-per-revolution terms cannot be
   * separated */
  static double covariance[PARAMETERS * PARAMETERS];
  options = (myb_fit_options){.max_iterations = 5, .tolerance = 1e-3, .epsilon = epsilon};
  const clock_t started = clock();
  MT_CHECK(myb_fit_orbit(set.forces, 0.0, y0, fit_obs, obs.count, &options, NULL, &result,
                         covariance, res) == MYB_OK);
  const double seconds = (double)(clock() - started) / (double)CLOCKS_PER_SEC;
  statistics(&obs, res, &rms, rtn);
  const size_t unknowns = 6 + myb_force_model_parameter_count(set.forces);
  printf("  least squares of the library (%zu unknowns, %zu values): RMS %.4f m (R %.4f, "
         "T %.4f, N %.4f) after %d corrections, %.1f s\n",
         unknowns, result.values, rms, rtn[0], rtn[1], rtn[2], result.iterations, seconds);
  const double* e = &set.empirical.constant[0];
  printf("  fitted: cr %.4f; along-track %.2f, %.2f cos u, %.2f sin u; cross-track %.2f cos u, "
         "%.2f sin u (pm/s^2)\n",
         set.object.cr, set.empirical.constant[1] * 1e12, set.empirical.cosine[1] * 1e12,
         set.empirical.sine[1] * 1e12, set.empirical.cosine[2] * 1e12,
         set.empirical.sine[2] * 1e12);
  (void)e;
  /* measured: 1.6 cm (R 2.8 mm, T 1.3 cm, N 8 mm; 2.0 cm before the secondary waves of the
   * ocean tides), constant along-track -2.2 pm/s^2 (the thermal drag of LAGEOS); the
   * once-per-revolution cross-track terms (~600 pm/s^2) absorb what differs from the ILRS
   * models (ocean tides EOT11a, S1/S2 atmospheric tides) — the terms the analysis centres
   * estimate for the same reason */
  if (m.off == 0 && !estimate_cr) {
    MT_CHECK(result.converged);
    MT_CHECK(rms < 0.02);
    MT_CHECK(fabs(set.empirical.constant[1]) < 5e-12);
  }
  disassemble(&set);
  myb_gravity_work_destroy(m.work);
  myb_gravity_field_destroy(m.field);
  myb_tides_destroy(m.tides);
  myb_spk_destroy(m.spk);
  myb_eop_destroy(m.eop);
  myb_cip_destroy(m.cip);
  myb_timesys_destroy(m.sys);
  return mt_finish("test_lageos");
}
