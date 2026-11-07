/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/lambert.h — Lambert's problem: the conic arcs that join two positions in a given time.
 *
 * Method: D. Izzo, "Revisiting Lambert's problem", Celest. Mech. Dyn. Astron. 121 (2015) 1-15:
 * the Lancaster-Blanchard time equation T(x) in the variable x, with Battin's hypergeometric
 * series near x = 1 (parabola). Instead of Izzo's Householder iteration each root is found by
 * Newton's method inside a guaranteed bracket: T decreases monotonically for zero revolutions,
 * and for M >= 1 revolutions it has a single minimum T_min that separates the left and right
 * branches. Each solution is checked in the tests by propagating (r1, v1) over the transfer
 * time with myb_kepler_propagate. */
#ifndef MYB_LAMBERT_H_
#define MYB_LAMBERT_H_

#include <stddef.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct myb_lambert_solution {
  double v1[3];    /* velocity at r1, m/s */
  double v2[3];    /* velocity at r2, m/s */
  int revolutions; /* complete revolutions before arrival */
  int branch;      /* 0 for zero revolutions; -1 left, +1 right branch for M >= 1 */
} myb_lambert_solution;

/* All solutions with 0..max_revolutions complete revolutions (at most 1 + 2 max_revolutions,
 * which `capacity` must allow). Direction: retrograde == 0 selects motion counter-clockwise
 * about +z (angular momentum with positive z component), retrograde != 0 the other sense.
 * tof > 0 in seconds, mu > 0 in m^3 s^-2.
 * Errors: MYB_E_DOMAIN (non-finite input, r1 or r2 zero, r1 and r2 collinear: the plane of the
 * transfer is undefined), MYB_E_NOCONV (iteration budget; not expected). Multi-revolution
 * solutions that do not exist for this tof are simply absent from the result. */
MYB_API MYB_NODISCARD myb_status myb_lambert(double mu, const double r1[3], const double r2[3],
                                             double tof, int retrograde, int max_revolutions,
                                             myb_lambert_solution* solutions, size_t capacity,
                                             size_t* count);

#ifdef __cplusplus
}
#endif

#endif /* MYB_LAMBERT_H_ */
