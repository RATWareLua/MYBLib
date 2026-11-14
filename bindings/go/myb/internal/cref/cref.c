/* SPDX-License-Identifier: BSD-3-Clause */
/* cref.c — the C reference of the tests of the Go binding: the tasks of the binding written as
 * a C program writes them, with the functions of the library called directly. The Go tests
 * compare what the binding gives with what these give, bit for bit. Not part of the binding. */
#include "cref.h"

#include <stdlib.h>
#include <string.h>

/* ---- A system of bodies with mergers --------------------------------------------------------- */

int cref_system_run(const double* gm, const double* radius, size_t bodies, const double* y0,
                    const double* times, size_t count, double* states, double* events,
                    size_t* event_count, double* y_end) {
  myb_system* system = NULL;
  myb_ias15* integrator = NULL;
  myb_status status = myb_system_create(bodies, gm, NULL, &system);
  for (size_t i = 0; i < bodies && status == MYB_OK; ++i) {
    if (radius[i] != 0.0) {
      status = myb_system_set_radius(system, i, radius[i]);
    }
  }
  if (status == MYB_OK) {
    status = myb_system_integrator(system, NULL, NULL, &integrator);
  }
  if (status == MYB_OK) {
    status = myb_ias15_start(integrator, myb_system_acceleration, system, 0.0, y0);
  }
  const size_t n = 2 * myb_system_coordinates(system);
  *event_count = 0;
  size_t k = 0;
  while (status == MYB_OK && k < count && *event_count < 64) {
    myb_system_event e;
    status = myb_system_advance(system, integrator, times[k], 0.0, &e);
    if (status != MYB_OK) {
      break;
    }
    if (e.kind == MYB_SYSTEM_ARRIVED) {
      memcpy(states + n * k, myb_ias15_state(integrator), n * sizeof(double));
      ++k;
      continue;
    }
    double* out = events + 7 * (*event_count)++;
    const myb_nbody* nbody = myb_system_nbody(system);
    const size_t into = myb_nbody_gm(nbody, e.b) > myb_nbody_gm(nbody, e.a) ? e.b : e.a;
    const size_t from = into == e.a ? e.b : e.a;
    myb_system_merger merger;
    status = myb_system_merge(system, integrator, into, from, &merger);
    out[0] = (double)e.kind;
    out[1] = (double)e.a;
    out[2] = (double)e.b;
    out[3] = e.t;
    out[4] = e.distance;
    out[5] = (double)into;
    out[6] = merger.energy;
  }
  if (integrator != NULL) {
    memcpy(y_end, myb_ias15_state(integrator), n * sizeof(double));
  }
  myb_ias15_destroy(integrator);
  myb_system_destroy(system);
  return (int)status;
}

int cref_invariants(const double* gm, size_t bodies, const double* y, size_t count,
                    double* energy_hi, double* momentum) {
  myb_nbody* nbody = NULL;
  myb_status status = myb_nbody_create(bodies, 0, gm, NULL, &nbody);
  const size_t n = 3 * bodies;
  for (size_t k = 0; k < count && status == MYB_OK; ++k) {
    myb_nbody_invariants inv;
    status = myb_nbody_invariants_of(nbody, y + 2 * n * k, y + 2 * n * k + n, &inv);
    energy_hi[k] = inv.energy.hi;
    memcpy(momentum + 3 * k, inv.momentum, sizeof inv.momentum);
  }
  myb_nbody_destroy(nbody);
  return (int)status;
}

int cref_system_partials(const double* gm, size_t bodies, const double* y0, double t1, double* y1,
                         double* phi, double* sensitivity) {
  myb_system* system = NULL;
  myb_status status = myb_system_create(bodies, gm, NULL, &system);
  const size_t all[1] = {0};
  if (status == MYB_OK) {
    status = myb_system_partials(system, NULL, 0.0, y0, t1, all, 1, y1, phi, sensitivity);
  }
  myb_system_destroy(system);
  return (int)status;
}

int cref_wh(const double* gm, size_t massive, const double* x, const double* v, double h,
            const double* times, size_t count, double* xs, double* vs) {
  myb_nbody* nbody = NULL;
  myb_wh* wh = NULL;
  const myb_wh_options options = {.h = h, .scheme = MYB_WH_SABA3, .corrector = 0};
  myb_status status = myb_nbody_create(massive, 0, gm, NULL, &nbody);
  if (status == MYB_OK) {
    status = myb_wh_create(nbody, &options, NULL, &wh);
  }
  if (status == MYB_OK) {
    status = myb_wh_start(wh, 0.0, x, v);
  }
  for (size_t k = 0; k < count && status == MYB_OK; ++k) {
    status = myb_wh_integrate(wh, times[k]);
    if (status == MYB_OK) {
      status = myb_wh_dense(wh, times[k], xs + 3 * massive * k, vs + 3 * massive * k);
    }
  }
  myb_wh_destroy(wh);
  myb_nbody_destroy(nbody);
  return (int)status;
}

/* ---- Estimation ---------------------------------------------------------------------------- */

int cref_filter(double gm, const double* y0, const double* covariance, double t_end,
                const double* times, const double* positions, size_t count, double sigma,
                double* innovations, double* step_t, double* step_unknowns, size_t* steps) {
  myb_force_model* model = NULL;
  myb_filter* filter = NULL;
  myb_point_mass_term term = {.gm = gm};
  myb_status status = myb_force_model_create(NULL, &model);
  if (status == MYB_OK) {
    status = myb_force_model_add(model, "planet", myb_term_point_mass, myb_term_point_mass_jacobian,
                                 &term, NULL);
  }
  const myb_filter_options options = {.history = count + 1};
  if (status == MYB_OK) {
    status = myb_filter_create(model, &options, NULL, &filter);
  }
  if (status == MYB_OK) {
    status = myb_filter_start(filter, 0.0, y0, covariance, t_end);
  }
  for (size_t k = 0; k < count && status == MYB_OK; ++k) {
    myb_observation o = {.t = times[k], .measure = myb_measure_position, .count = 3};
    for (int i = 0; i < 3; ++i) {
      o.value[i] = positions[3 * k + (size_t)i];
      o.sigma[i] = sigma;
    }
    status = myb_filter_update(filter, &o, innovations + 3 * k);
  }
  if (status == MYB_OK) {
    status = myb_filter_smooth(filter);
  }
  *steps = 0;
  if (status == MYB_OK) {
    *steps = myb_filter_steps(filter);
    for (size_t j = 0; j < *steps; ++j) {
      myb_filter_step(filter, j, &step_t[j], step_unknowns + 6 * j, NULL);
    }
  }
  myb_filter_destroy(filter);
  myb_force_model_destroy(model);
  return (int)status;
}

int cref_range(double gm, double longitude, double latitude, double radius, const double* t,
               const double* y, size_t count, double* values, double* partials) {
  myb_env_options options = {0};
  myb_env* env = NULL;
  myb_status status = myb_env_create(&options, NULL, &env);
  const myb_env_body body = {.id = 399, .gm = gm};
  size_t index = 0;
  if (status == MYB_OK) {
    status = myb_env_add_body(env, &body, &index);
  }
  myb_station station;
  const myb_ellipsoid ellipsoid = {.a = radius, .f = 0.0};
  if (status == MYB_OK) {
    status = myb_station_init(&station, env, 0, &ellipsoid, longitude, latitude, 0.0);
  }
  for (size_t k = 0; k < count && status == MYB_OK; ++k) {
    status = myb_measure_range(&station, t[k], y + 6 * k, values + k, partials + 6 * k);
  }
  myb_env_destroy(env);
  return (int)status;
}

/* ---- The restricted problem ---------------------------------------------------------------------
 */

int cref_cr3bp_family(double mu, const double* y0, double period, size_t patches, double step,
                      size_t count, double* periods, double* states) {
  myb_cr3bp_orbit* orbit = malloc(sizeof *orbit);
  if (orbit == NULL) {
    return MYB_E_NOMEM;
  }
  myb_status status = myb_cr3bp_orbit_from_state(mu, y0, period, patches, orbit);
  if (status == MYB_OK) {
    status = myb_cr3bp_orbit_correct(mu, orbit, 0.0, 0);
  }
  for (size_t k = 0; k < count && status == MYB_OK; ++k) {
    if (k > 0) {
      status = myb_cr3bp_orbit_continue(mu, orbit, step, 0.0, 0);
    }
    periods[k] = orbit->period;
    memcpy(states + 6 * patches * k, orbit->states, 6 * patches * sizeof(double));
  }
  free(orbit);
  return (int)status;
}

int cref_torus(double mu, const double* y0, double period, size_t points, double size, double step,
               double* rotations, double* normal) {
  myb_cr3bp_torus* torus = malloc(sizeof *torus);
  myb_cr3bp_torus_spectrum* spectrum = malloc(sizeof *spectrum);
  myb_status status = torus != NULL && spectrum != NULL ? MYB_OK : MYB_E_NOMEM;
  if (status == MYB_OK) {
    status = myb_cr3bp_torus_from_orbit(mu, y0, period, points, 0, size, torus);
  }
  if (status == MYB_OK) {
    status = myb_cr3bp_torus_continue(mu, torus, MYB_CR3BP_TORUS_JACOBI, 0.0, 0.0, 0);
    rotations[0] = torus->rotation;
  }
  if (status == MYB_OK) {
    status = myb_cr3bp_torus_continue(mu, torus, MYB_CR3BP_TORUS_JACOBI, step, 0.0, 0);
    rotations[1] = torus->rotation;
  }
  if (status == MYB_OK) {
    status = myb_cr3bp_torus_stability(mu, torus, 0.0, spectrum);
  }
  if (status == MYB_OK) {
    normal[0] = spectrum->normal_rotation;
    normal[1] = spectrum->normal_eigenvalue[0];
    normal[2] = spectrum->normal_eigenvalue[1];
    normal[3] = spectrum->normal_re[0];
    normal[4] = spectrum->normal_im[5];
    normal[5] = (double)spectrum->elliptic;
  }
  free(spectrum);
  free(torus);
  return (int)status;
}

int cref_manifold(double mu, const double* y0, double period, double phase, double* y) {
  return (int)myb_cr3bp_manifold(mu, y0, period, phase, 1, 1, 1e-6, y);
}

/* ---- Targeting and transfers ---------------------------------------------------------------- */

static myb_status kepler_body(const double* state, myb_time epoch, double mu, myb_time t,
                              double out[6]) {
  myb_pv pv0;
  memcpy(pv0.r, state, sizeof pv0.r);
  memcpy(pv0.v, state + 3, sizeof pv0.v);
  myb_pv pv;
  const myb_status status =
      myb_kepler_propagate(mu, &pv0, myb_dd_to_double(myb_time_diff(t, epoch)), &pv);
  memcpy(out, pv.r, sizeof pv.r);
  memcpy(out + 3, pv.v, sizeof pv.v);
  return status;
}

int cref_lambert_cell(double mu, const double* from, const double* to, double departure,
                      double arrival, double* v1, double* v2) {
  const myb_time zero = myb_time_from_seconds(0.0);
  const myb_time td = myb_time_from_seconds(departure);
  const myb_time ta = myb_time_from_seconds(arrival);
  double a[6];
  double b[6];
  myb_status status = kepler_body(from, zero, mu, td, a);
  if (status == MYB_OK) {
    status = kepler_body(to, zero, mu, ta, b);
  }
  myb_lambert_solution s;
  size_t n = 0;
  if (status == MYB_OK) {
    status = myb_lambert(mu, a, b, myb_dd_to_double(myb_time_diff(ta, td)), 0, 0, &s, 1, &n);
  }
  memcpy(v1, s.v1, sizeof s.v1);
  memcpy(v2, s.v2, sizeof s.v2);
  return (int)status;
}

int cref_conjunctions(double gm, const double* y1, const double* y2, const double* covariance,
                      double radius, double t1, double step, double* first_t,
                      double* first_probability, size_t* count) {
  myb_force_model* model = NULL;
  myb_point_mass_term term = {.gm = gm};
  myb_status status = myb_force_model_create(NULL, &model);
  if (status == MYB_OK) {
    status = myb_force_model_add(model, "planet", myb_term_point_mass, myb_term_point_mass_jacobian,
                                 &term, NULL);
  }
  const myb_orbit_options options = {.variational = 1};
  myb_orbit* a = NULL;
  myb_orbit* b = NULL;
  if (status == MYB_OK) {
    status = myb_orbit_create(model, &options, NULL, &a);
  }
  if (status == MYB_OK) {
    status = myb_orbit_create(model, &options, NULL, &b);
  }
  if (status == MYB_OK) {
    status = myb_orbit_start(a, 0.0, y1, 2.0 * t1);
  }
  if (status == MYB_OK) {
    status = myb_orbit_start(b, 0.0, y2, 2.0 * t1);
  }
  myb_conjunction found[16];
  *count = 0;
  if (status == MYB_OK) {
    status = myb_conjunction_find(a, covariance, b, covariance, radius, t1, step, found, 16, count);
  }
  if (status == MYB_OK && *count > 0) {
    *first_t = found[0].approach.t;
    *first_probability = found[0].probability;
  }
  myb_orbit_destroy(b);
  myb_orbit_destroy(a);
  myb_force_model_destroy(model);
  return (int)status;
}

int cref_qlaw_fly(double mu, const double* target, const double* weight, double thrust,
                  double exhaust, const double* y0, double t_limit, double step, double* y,
                  double* t) {
  const myb_qlaw law = {.mu = mu,
                        .target = {target[0], target[1], target[2], target[3], target[4]},
                        .weight = {weight[0], weight[1], weight[2], weight[3], weight[4]}};
  const myb_qlaw_term engine = {.law = &law,
                                .thrust = thrust,
                                .exhaust = exhaust,
                                .tolerance = {1e3, 1e-3, 1e-3, 1e-3, 1e-3}};
  const myb_qlaw_flight flight = {.engine = &engine};
  myb_qlaw_arrival arrival;
  const myb_status status = myb_qlaw_fly(&flight, 0.0, y0, t_limit, step, NULL, &arrival);
  memcpy(y, arrival.y, sizeof arrival.y);
  *t = arrival.t;
  return (int)status;
}

/* The Earth, Venus and Mars of test_chain: circles as the conics of the Sun from states at 0. */
typedef struct cref_body {
  double state[6];
  double mu;
} cref_body;

static myb_status cref_body_state(void* context, myb_time t, double state[6]) {
  const cref_body* b = context;
  return kepler_body(b->state, myb_time_from_seconds(0.0), b->mu, t, state);
}

int cref_chain(double mu, const double* states, const double* gm, const double* periapsis,
               const double* impulse_exhaust, const double* bounds, double* x, size_t unknowns,
               double* propellant, double* constraints, size_t constraint_count) {
  cref_body bodies[3];
  myb_chain chain = {.mu = mu,
                     .epoch = myb_time_from_seconds(0.0),
                     .thrust = 1.0,
                     .exhaust = 30000.0,
                     .mass0 = 1000.0,
                     .legs = 2,
                     .segments = {6, 6},
                     .departure = {bounds[0], bounds[1]},
                     .flight = {{bounds[2], bounds[3]}, {bounds[4], bounds[5]}},
                     .launch = 2500.0,
                     .arrival = 0.0};
  for (int k = 0; k < 3; ++k) {
    memcpy(bodies[k].state, states + 6 * k, sizeof bodies[k].state);
    bodies[k].mu = mu;
    chain.bodies[k] = (myb_chain_body){cref_body_state, &bodies[k], gm[k], periapsis[k]};
    chain.impulse_exhaust[k] = impulse_exhaust[k];
  }
  size_t u = 0;
  size_t e = 0;
  size_t i = 0;
  myb_status status = myb_chain_size(&chain, &u, &e, &i);
  if (status == MYB_OK && (u != unknowns || e + i != constraint_count)) {
    status = MYB_E_DOMAIN;
  }
  if (status == MYB_OK) {
    status = myb_chain_guess(&chain, x);
  }
  if (status == MYB_OK) {
    status = myb_chain_evaluate(&chain, x, propellant, constraints, NULL, NULL, NULL);
  }
  return (int)status;
}

int cref_sims_flanagan(double mu, const double* start, const double* target, double t1,
                       size_t segments, const double* throttles, double* mismatch, double* mass1) {
  myb_sims_flanagan leg = {.mu = mu,
                           .thrust = 0.5,
                           .exhaust = 30000.0,
                           .t0 = 0.0,
                           .mass0 = 1000.0,
                           .t1 = t1,
                           .segments = segments};
  memcpy(leg.start, start, sizeof leg.start);
  memcpy(leg.target, target, sizeof leg.target);
  return (int)myb_sims_flanagan_mismatch(&leg, throttles, mismatch, mass1, NULL);
}

/* ---- Attitude and observables -------------------------------------------------------------- */

int cref_sun_nadir(double gm, const double* sun, double sun_gm, const double* t, const double* y,
                   size_t count, double* m) {
  myb_env_options options = {0};
  myb_env* env = NULL;
  myb_status status = myb_env_create(&options, NULL, &env);
  cref_body source = {.mu = sun_gm};
  memcpy(source.state, sun, sizeof source.state);
  const myb_env_body earth = {.id = 399, .gm = gm};
  const myb_env_body body = {
      .id = 10, .gm = sun_gm, .state = cref_body_state, .state_context = &source};
  size_t index = 0;
  if (status == MYB_OK) {
    status = myb_env_add_body(env, &earth, &index);
  }
  if (status == MYB_OK) {
    status = myb_env_add_body(env, &body, &index);
  }
  myb_sun_attitude law = {.env = env, .body = 0, .source = 1};
  for (size_t k = 0; k < count && status == MYB_OK; ++k) {
    status = myb_attitude_sun_nadir(&law, t[k], y + 6 * k, y + 6 * k + 3, m + 9 * k);
  }
  myb_env_destroy(env);
  return (int)status;
}

int cref_light_time(double c, double mu, const double* target, const double* observer, double t,
                    double* light_time, double* range) {
  cref_body body = {.mu = mu};
  memcpy(body.state, target, sizeof body.state);
  const myb_light_model model = {.c = c, .gamma = 1.0};
  myb_light_path path;
  const myb_status status = myb_light_time(&model, cref_body_state, &body, myb_time_from_seconds(t),
                                           observer, MYB_LIGHT_RECEPTION, &path);
  *light_time = path.light_time;
  *range = path.range;
  return (int)status;
}

/* ---- Files -------------------------------------------------------------------------------- */

int cref_ck_pointing(const char* ck_path, int instrument, double ticks, double* q) {
  myb_ck* ck = NULL;
  myb_status status = myb_ck_create(NULL, &ck);
  if (status == MYB_OK) {
    status = myb_ck_load_file(ck, ck_path, NULL);
  }
  myb_ck_pointing p;
  if (status == MYB_OK) {
    status = myb_ck_pointing_at(ck, instrument, ticks, 0.0, 0, &p);
  }
  if (status == MYB_OK) {
    memcpy(q, p.q, sizeof p.q);
  }
  myb_ck_destroy(ck);
  return (int)status;
}

int cref_sclk_time(const char* kernel, long long id, double ticks, long long* s, double* f) {
  myb_kpool* pool = NULL;
  myb_sclk* sclk = NULL;
  myb_status status = myb_kpool_create(NULL, &pool);
  if (status == MYB_OK) {
    status = myb_kpool_load_file(pool, kernel, NULL);
  }
  if (status == MYB_OK) {
    status = myb_sclk_create(pool, id, NULL, &sclk, NULL);
  }
  myb_time t = {0};
  if (status == MYB_OK) {
    status = myb_sclk_to_time(sclk, ticks, &t);
  }
  *s = t.s;
  *f = t.f;
  myb_sclk_destroy(sclk);
  myb_kpool_destroy(pool);
  return (int)status;
}
