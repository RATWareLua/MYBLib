/* SPDX-License-Identifier: BSD-3-Clause */
/* dyn/nbody_internal.h — evaluation paths of the N-body accelerations (for the tests, which
 * prove that they agree bit for bit) and facts about the graph used by the integrators. */
#ifndef MYB_SRC_DYN_NBODY_INTERNAL_H_
#define MYB_SRC_DYN_NBODY_INTERNAL_H_

#include "myb/dyn.h"

typedef enum myb_internal_nbody_path {
  MYB_NBODY_PATH_AUTO = 0,      /* what myb_nbody_acceleration does */
  MYB_NBODY_PATH_SYMMETRIC = 1, /* each pair once, single thread (needs no exclusion) */
  MYB_NBODY_PATH_RECEIVER = 2,  /* each receiver summed on its own, parallel when a pool is set,
                                   four receivers per vector block */
  MYB_NBODY_PATH_RECEIVER_SCALAR = 3, /* the same one receiver at a time (the reference) */
} myb_internal_nbody_path;

[[nodiscard]] myb_status myb_internal_nbody_acceleration_path(myb_nbody* nbody, const double* x,
                                                              double* a,
                                                              myb_internal_nbody_path path);

/* The same with velocities (v may be NULL when no 1PN term is active). */
[[nodiscard]] myb_status myb_internal_nbody_acceleration_xv_path(myb_nbody* nbody,
                                                                 const double* x,
                                                                 const double* v, double* a,
                                                                 myb_internal_nbody_path path);

/* 1 if the accelerations depend on the velocities (1PN terms). */
int myb_internal_nbody_velocity_dependent(const myb_nbody* nbody);

/* 1 if the massive bodies interact pairwise both ways (every exclusion between two massive
 * bodies is mirrored), so that sum GM_i a_i = 0 over the massive bodies (Newton's third law). */
int myb_internal_nbody_reciprocal(myb_nbody* nbody);

#endif /* MYB_SRC_DYN_NBODY_INTERNAL_H_ */
