/* SPDX-License-Identifier: BSD-3-Clause */
/* kepler/drift.h — the two-body flow as an increment of the state, for splitting integrators
 * that accumulate the state with compensated summation (dyn/whfast.c). */
#ifndef MYB_SRC_KEPLER_DRIFT_H_
#define MYB_SRC_KEPLER_DRIFT_H_

#include "myb/kepler.h"

/* delta such that pv0 + delta is the state of myb_kepler_propagate(mu, pv0, dt) (the same
 * routes and errors). On the universal route delta is computed directly, accurate relative to
 * its own size; on the hyperbola restarted from the periapsis (|H0| > 1) it is the difference
 * of the propagated state and pv0. */
[[nodiscard]] myb_status myb_internal_kepler_delta(double mu, const myb_pv* pv0, double dt,
                                                   myb_pv* delta);

#endif /* MYB_SRC_KEPLER_DRIFT_H_ */
