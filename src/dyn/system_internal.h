/* SPDX-License-Identifier: BSD-3-Clause */
/* dyn/system_internal.h — the derivatives of a system of bodies by central differences of whole
 * propagations, what myb_system_partials (myb/system.h) computed before it integrated the
 * variational equations: kept as the reference of the tests (tests/unit/test_system.c compares
 * the two on a system with every level; tests/repro hashes it, as it did). */
#ifndef MYB_SRC_DYN_SYSTEM_INTERNAL_H_
#define MYB_SRC_DYN_SYSTEM_INTERNAL_H_

#include "myb/system.h"

/* The arguments, the outputs and the errors of myb_system_partials. Each column of phi costs two
 * propagations of the system from y0 changed by +-h in one coordinate (h = 1e-6 of the largest
 * coordinate of its kind — the positions of the bodies, then each group of three of the
 * rotations and the spins — or of their rates), each gm two propagations with it changed by
 * +-1e-6 of itself; the error is that of the differences, a part in 10^9 .. 10^11 of a
 * derivative while the motion does not branch. */
[[nodiscard]] myb_status myb_internal_system_partials_differences(
    myb_system* system, const myb_ias15_options* options, double t0, const double* y0, double t1,
    const size_t* bodies, size_t count, double* y1, double* phi, double* sensitivity);

/* The Jacobian that the variational equations of myb_system_partials integrate, at the state
 * (t, y) (2 n values, n = myb_system_coordinates): jacobian (n x 2 n, row-major) = d a / d x and
 * then d a / d v of a = myb_system_acceleration, and dgm (n x count, may be NULL with count 0) =
 * d a / d gm_k of the bodies given — every level as myb_system_partials takes it, analytic or by
 * differences. For the tests, which check it level by level against differences of the
 * acceleration. Errors: as myb_system_partials and myb_system_acceleration. */
[[nodiscard]] myb_status myb_internal_system_jacobian(myb_system* system, double t, const double* y,
                                                      const size_t* bodies, size_t count,
                                                      double* jacobian, double* dgm);

#endif /* MYB_SRC_DYN_SYSTEM_INTERNAL_H_ */
