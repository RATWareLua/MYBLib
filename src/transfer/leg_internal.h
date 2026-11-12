/* SPDX-License-Identifier: BSD-3-Clause */
/* transfer/leg_internal.h — a leg of Sims and Flanagan in any field: the coasts between the
 * impulses with their transition matrices, and the sweep of a leg that gives its mismatch with
 * the derivatives by everything it depends on — the impulses, the states at its two ends, the
 * length of its segments and the time it starts (for sims_flanagan.c and chain.c). Everything
 * is in the units of the solver: a length, the time sqrt(length^3 / mu), mu = 1. */
#ifndef MYB_SRC_TRANSFER_LEG_INTERNAL_H_
#define MYB_SRC_TRANSFER_LEG_INTERNAL_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/integ.h"
#include "myb/model.h"

/* What carries a state between two impulses. */
typedef struct myb_internal_coast {
  myb_force_model* field; /* NULL: the conic of mu = 1 */
  double length;          /* m */
  double time;            /* s */
  double origin;          /* s: the model time of the time 0 of the solver */
  myb_dop853* integrator; /* with a field: the state and its transition matrix, 42 values */
} myb_internal_coast;

/* allocator may be NULL. Errors: MYB_E_NOMEM. */
[[nodiscard]] myb_status myb_internal_coast_open(myb_internal_coast* coast, myb_force_model* field,
                                                 double length, double time, double origin,
                                                 const myb_allocator* allocator);
void myb_internal_coast_close(myb_internal_coast* coast); /* accepts a coast never opened */

/* The state y after dt (either sign) from y0 at the time t, and the transition matrix phi
 * (row-major) when it is not NULL. */
[[nodiscard]] myb_status myb_internal_coast_fly(myb_internal_coast* coast, double t,
                                                const double y0[6], double dt, double y[6],
                                                double* phi);

/* The rate of the state, (v, g(t, x, v)). */
[[nodiscard]] myb_status myb_internal_coast_rate(myb_internal_coast* coast, double t,
                                                 const double y[6], double rate[6]);

/* A leg: n segments of length dt from the time t, the first `forward` of them flown from
 * `start`, the others backward from `end`, an impulse (a change of velocity) in the middle of
 * each. */
typedef struct myb_internal_leg {
  size_t n;
  size_t forward;
  double t;
  double dt;
  double start[6];
  double end[6];
  const double* impulses; /* 3 n */
} myb_internal_leg;

/* The mismatch c = forward - backward at the meeting point and, at fixed impulses, its
 * derivatives by the state at the start and at the end (6 x 6), by dt (the end and the
 * impulses moving with it, t fixed) and by t (the whole leg moved in time; zero in a field that
 * does not depend on time). */
typedef struct myb_internal_sweep {
  double c[6];
  double by_start[36];
  double by_end[36];
  double by_step[6];
  double by_shift[6];
} myb_internal_sweep;

/* Doubles of work a sweep of n segments needs. */
size_t myb_internal_sweep_work(size_t n);

/* The sweep of a leg. g (may be NULL): dc / d impulse, 6 x 3 n row-major — with it, and with
 * `full`, the other derivatives too. work: myb_internal_sweep_work(n) doubles (used only for
 * derivatives). */
[[nodiscard]] myb_status myb_internal_leg_sweep(myb_internal_coast* coast,
                                                const myb_internal_leg* leg, double* g, int full,
                                                double* work, myb_internal_sweep* out);

#endif /* MYB_SRC_TRANSFER_LEG_INTERNAL_H_ */
