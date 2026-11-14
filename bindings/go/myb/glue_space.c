/* SPDX-License-Identifier: BSD-3-Clause */
/* glue_space.c — the restricted problem of three bodies, targeting, relative motion, close
 * approaches and low-thrust transfers for the Go binding (glue_tasks.h). */
#include <stdlib.h>
#include <string.h>

#include "glue_internal.h"

static void to_pv(const double* y, myb_pv* pv) {
  memcpy(pv->r, y, sizeof pv->r);
  memcpy(pv->v, y + 3, sizeof pv->v);
}

static void from_pv(const myb_pv* pv, double* y) {
  memcpy(y, pv->r, sizeof pv->r);
  memcpy(y + 3, pv->v, sizeof pv->v);
}

/* ---- The restricted problem of three bodies ----------------------------------------------------
 */

myb_status mybgo_cr3bp_family(double mu, const myb_cr3bp_orbit* orbit, double step, size_t count,
                              double tolerance, int max_iterations, myb_cr3bp_orbit* members,
                              size_t* made) {
  *made = 0;
  if (count == 0) {
    return MYB_OK;
  }
  members[0] = *orbit;
  myb_status status = myb_cr3bp_orbit_correct(mu, &members[0], tolerance, max_iterations);
  for (size_t k = 1; k < count && status == MYB_OK; ++k) {
    members[k] = members[k - 1];
    status = myb_cr3bp_orbit_continue(mu, &members[k], step, tolerance, max_iterations);
    if (status == MYB_OK) {
      *made = k;
    }
  }
  if (status == MYB_OK) {
    *made = count;
  }
  return status;
}

myb_status mybgo_cr3bp_torus_family(double mu, const myb_cr3bp_torus* torus,
                                    myb_cr3bp_torus_family family, double step, size_t count,
                                    double tolerance, int max_iterations, myb_cr3bp_torus* members,
                                    size_t* made) {
  *made = 0;
  if (count == 0) {
    return MYB_OK;
  }
  members[0] = *torus;
  myb_status status =
      myb_cr3bp_torus_continue(mu, &members[0], family, 0.0, tolerance, max_iterations);
  for (size_t k = 1; k < count && status == MYB_OK; ++k) {
    members[k] = members[k - 1];
    status = myb_cr3bp_torus_continue(mu, &members[k], family, step, tolerance, max_iterations);
    if (status == MYB_OK) {
      *made = k;
    }
  }
  if (status == MYB_OK) {
    *made = count;
  }
  return status;
}

void mybgo_cr3bp_jacobi(double mu, const double* y, size_t count, double* jacobi) {
  for (size_t k = 0; k < count; ++k) {
    jacobi[k] = myb_cr3bp_jacobi(mu, y + 6 * k);
  }
}

myb_status mybgo_cr3bp_manifold(double mu, const double* y0, double period, const double* phase,
                                size_t count, int unstable, int sign, double epsilon, double* y,
                                size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_cr3bp_manifold(mu, y0, period, phase[k], unstable, sign, epsilon, y + 6 * k);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

void mybgo_cr3bp_torus_states(const myb_cr3bp_torus* torus, const double* angle, size_t count,
                              double* y) {
  for (size_t k = 0; k < count; ++k) {
    myb_cr3bp_torus_state(torus, angle[k], y + 6 * k);
  }
}

myb_status mybgo_cr3bp_torus_manifold(const myb_cr3bp_torus* torus,
                                      const myb_cr3bp_torus_spectrum* spectrum, const double* angle,
                                      size_t count, int unstable, int sign, double epsilon,
                                      double* y, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status =
        myb_cr3bp_torus_manifold(torus, spectrum, angle[k], unstable, sign, epsilon, y + 6 * k);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_cr3bp_propagate(double mu, const double* y0, size_t count, double duration,
                                 const myb_cr3bp_section* section, double* y, double* elapsed,
                                 int* crossed, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status =
        myb_cr3bp_propagate(mu, y0 + 6 * k, duration, section, y + 6 * k, &elapsed[k], &crossed[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_cr3bp_frames(int to_ephemeris, const double* body1, const double* body2,
                              double gm1, double gm2, const double* in, size_t count, double* out,
                              size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = to_ephemeris ? myb_cr3bp_to_ephemeris(body1 + 6 * k, body2 + 6 * k, gm1, gm2,
                                                   in + 6 * k, out + 6 * k)
                          : myb_cr3bp_from_ephemeris(body1 + 6 * k, body2 + 6 * k, gm1, gm2,
                                                     in + 6 * k, out + 6 * k);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

/* ---- Targeting ------------------------------------------------------------------------------- */

myb_status mybgo_lambert_many(double mu, const double* r1, const double* r2, const double* tof,
                              size_t count, int retrograde, int max_revolutions, size_t capacity,
                              myb_lambert_solution* solutions, size_t* counts, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_lambert(mu, r1 + 3 * k, r2 + 3 * k, tof[k], retrograde, max_revolutions,
                         solutions + capacity * k, capacity, &counts[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_lambert_grid(double mu, const mybgo_provider* departure,
                              const mybgo_provider* arrival, const myb_time* departures,
                              size_t departure_count, const myb_time* arrivals,
                              size_t arrival_count, int retrograde, int max_revolutions,
                              size_t capacity, double* departure_states, double* arrival_states,
                              myb_lambert_solution* solutions, size_t* counts, size_t* failed) {
  myb_status status =
      mybgo_provider_states(departure, departures, departure_count, departure_states, failed);
  if (status == MYB_OK) {
    status = mybgo_provider_states(arrival, arrivals, arrival_count, arrival_states, failed);
  }
  if (status != MYB_OK) {
    return status;
  }
  size_t cell = 0;
  for (size_t i = 0; i < departure_count && status == MYB_OK; ++i) {
    for (size_t j = 0; j < arrival_count; ++j) {
      cell = i * arrival_count + j;
      counts[cell] = 0;
      if (myb_time_compare(arrivals[j], departures[i]) <= 0) {
        continue;
      }
      const double tof = myb_dd_to_double(myb_time_diff(arrivals[j], departures[i]));
      status = myb_lambert(mu, departure_states + 6 * i, arrival_states + 6 * j, tof, retrograde,
                           max_revolutions, solutions + capacity * cell, capacity, &counts[cell]);
      if (status != MYB_OK) {
        break;
      }
    }
  }
  return mybgo_end(status, status == MYB_OK ? departure_count * arrival_count : cell, failed);
}

myb_status mybgo_bplanes(double gm, const double* states, size_t count, const double pole[3],
                         myb_bplane* bplanes, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_bplane_compute(gm, states + 6 * k, states + 6 * k + 3, pole, &bplanes[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_shoot(uintptr_t handle, size_t n, size_t m, double* u,
                       const myb_shoot_options* options, int* iterations, double* residual) {
  if (n == 0 || m == 0 || n > 64 || m > 64) {
    return MYB_E_DOMAIN;
  }
  double* work = malloc(myb_shoot_work_size(n, m) * sizeof(double));
  if (work == NULL) {
    return MYB_E_NOMEM;
  }
  mybgo_go_shoot_context context = {.handle = handle, .n = n};
  const myb_status status =
      myb_shoot(mybgo_go_shoot, &context, n, m, u, options, work, iterations, residual);
  free(work);
  return status;
}

myb_status mybgo_relative(int kind, double mu, const double* chief, const double* other,
                          size_t count, double* out, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    myb_pv c;
    myb_pv o;
    to_pv(chief + 6 * k, &c);
    switch (kind) {
      case MYBGO_RELATIVE_FROM_INERTIAL:
        to_pv(other + 6 * k, &o);
        status = myb_relative_from_inertial(&c, &o, out + 6 * k);
        break;
      case MYBGO_RELATIVE_TO_INERTIAL:
        status = myb_relative_to_inertial(&c, other + 6 * k, &o);
        if (status == MYB_OK) {
          from_pv(&o, out + 6 * k);
        }
        break;
      case MYBGO_RELATIVE_ROE:
        to_pv(other + 6 * k, &o);
        status = myb_roe_from_states(mu, &c, &o, out + 6 * k);
        break;
      default:
        status = MYB_E_DOMAIN;
    }
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_relative_stm(double mu, const double* chief, double n, const double* dt,
                              size_t count, double* phi, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  myb_pv c;
  if (chief != NULL) {
    to_pv(chief, &c);
  } else if (n == 0.0) {
    return mybgo_end(MYB_E_DOMAIN, 0, failed);
  }
  for (; k < count; ++k) {
    if (chief != NULL) {
      status = myb_relative_stm(mu, &c, dt[k], phi + 36 * k);
    } else {
      myb_cw_stm(n, dt[k], phi + 36 * k);
    }
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

void mybgo_roe_positions(double a, const double* u, size_t count, const double roe[6],
                         double* positions) {
  for (size_t k = 0; k < count; ++k) {
    myb_roe_to_position(a, u[k], roe, positions + 3 * k);
  }
}

/* ---- Close approaches -------------------------------------------------------------------------
 */

myb_status mybgo_approach_find(const mybgo_provider* first, const mybgo_provider* second, double t0,
                               double t1, double step, myb_approach** approaches, size_t* count) {
  *approaches = NULL;
  *count = 0;
  mybgo_provider a = *first;
  mybgo_provider b = *second;
  size_t capacity = 16;
  myb_approach* list = NULL;
  myb_status status = MYB_OK;
  /* bound: the second pass has room for every approach the first one counted */
  for (int pass = 0; pass < 2; ++pass) {
    free(list);
    list = calloc(capacity, sizeof *list);
    if (list == NULL) {
      return MYB_E_NOMEM;
    }
    size_t n = 0;
    status = myb_approach_find(mybgo_provider_state, &a, mybgo_provider_state, &b, t0, t1, step,
                               list, capacity, &n);
    if (status != MYB_OK || n <= capacity) {
      *count = status == MYB_OK ? n : 0;
      break;
    }
    capacity = n;
  }
  if (status == MYB_OK) {
    *approaches = list;
  } else {
    free(list);
  }
  return status;
}

/* ---- Low thrust -------------------------------------------------------------------------------
 */

myb_status mybgo_qlaw_directions(const myb_qlaw* law, const myb_elements* elements, size_t count,
                                 double* direction, double* q, double* q_rate, double* effectivity,
                                 size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_qlaw_direction(law, &elements[k], direction + 3 * k, &q[k], &q_rate[k],
                                &effectivity[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_qlaw_fly(const myb_qlaw* law, const myb_qlaw_term* engine, myb_force_model* field,
                          double rtol, double atol, double t0, const double* y0, double t_limit,
                          double step, myb_qlaw_arrival* arrival) {
  myb_qlaw_term term = *engine;
  term.law = law;
  const myb_qlaw_flight flight = {.engine = &term, .field = field, .rtol = rtol, .atol = atol};
  return myb_qlaw_fly(&flight, t0, y0, t_limit, step, NULL, arrival);
}

/* The chain with its providers in place (copied: the contexts live for the call). */
typedef struct chain_copy {
  myb_chain chain;
  mybgo_provider providers[MYB_CHAIN_MAX_LEGS + 1];
} chain_copy;

static myb_status chain_of(const myb_chain* chain, const mybgo_provider* providers,
                           chain_copy* copy) {
  copy->chain = *chain;
  if (chain->legs == 0 || chain->legs > MYB_CHAIN_MAX_LEGS) {
    return MYB_E_DOMAIN;
  }
  for (size_t k = 0; k <= chain->legs; ++k) {
    copy->providers[k] = providers[k];
    copy->chain.bodies[k].state = mybgo_provider_state;
    copy->chain.bodies[k].context = &copy->providers[k];
  }
  return MYB_OK;
}

myb_status mybgo_chain_size(const myb_chain* chain, const mybgo_provider* providers,
                            size_t* unknowns, size_t* equalities, size_t* inequalities) {
  chain_copy* copy = malloc(sizeof *copy);
  if (copy == NULL) {
    return MYB_E_NOMEM;
  }
  myb_status status = chain_of(chain, providers, copy);
  if (status == MYB_OK) {
    status = myb_chain_size(&copy->chain, unknowns, equalities, inequalities);
  }
  free(copy);
  return status;
}

myb_status mybgo_chain_guess(const myb_chain* chain, const mybgo_provider* providers, double* x) {
  chain_copy* copy = malloc(sizeof *copy);
  if (copy == NULL) {
    return MYB_E_NOMEM;
  }
  myb_status status = chain_of(chain, providers, copy);
  if (status == MYB_OK) {
    status = myb_chain_guess(&copy->chain, x);
  }
  free(copy);
  return status;
}

myb_status mybgo_chain_evaluate(const myb_chain* chain, const mybgo_provider* providers,
                                const double* x, double* propellant, double* constraints,
                                double* gradient, double* jacobian) {
  chain_copy* copy = malloc(sizeof *copy);
  if (copy == NULL) {
    return MYB_E_NOMEM;
  }
  myb_status status = chain_of(chain, providers, copy);
  if (status == MYB_OK) {
    status = myb_chain_evaluate(&copy->chain, x, propellant, constraints, gradient, jacobian, NULL);
  }
  free(copy);
  return status;
}

myb_status mybgo_chain_solve(const myb_chain* chain, const mybgo_provider* providers,
                             const myb_nlp_options* options, double* x, myb_chain_result* result) {
  chain_copy* copy = malloc(sizeof *copy);
  if (copy == NULL) {
    return MYB_E_NOMEM;
  }
  myb_status status = chain_of(chain, providers, copy);
  if (status == MYB_OK) {
    status = myb_chain_solve(&copy->chain, options, x, NULL, result);
  }
  free(copy);
  return status;
}
