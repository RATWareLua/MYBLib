/* SPDX-License-Identifier: BSD-3-Clause */
/* dyn/fmm.h — the fast multipole method for the Newtonian accelerations of point masses: the
 * approximate level of the N-body graph (myb_nbody_set_approximation in myb/dyn.h; the method
 * and its measured error are described in dyn/fmm.c). */
#ifndef MYB_SRC_DYN_FMM_H_
#define MYB_SRC_DYN_FMM_H_

#include <stddef.h>

#include "myb/core.h"

typedef struct myb_internal_fmm myb_internal_fmm;

/* What the last evaluation did. */
typedef struct myb_internal_fmm_stats {
  size_t cells;        /* cells of the octree */
  size_t leaves;       /* of them leaves */
  size_t depth;        /* deepest level (the root is level 0) */
  size_t m2l;          /* pairs of cells that interacted through their expansions, each once */
  size_t p2p;          /* pairs of points summed directly, each once */
  size_t leaf_size;    /* at most this many points in a leaf (unless at the depth limit) */
  size_t direct_pairs; /* a pair of cells with at most this many pairs of points is summed
                          directly, whatever the opening angle */
} myb_internal_fmm_stats;

/* An evaluator for `count` points (count >= 1) of expansion order 1 <= order <=
 * MYB_NBODY_FMM_MAX_ORDER and opening angle 0 < theta < 1; leaf_size 0 chooses the leaf size
 * from the order. Errors: MYB_E_NOMEM, MYB_E_OVERFLOW. */
[[nodiscard]] myb_status myb_internal_fmm_create(size_t count, int order, double theta,
                                                 size_t leaf_size, const myb_allocator* allocator,
                                                 myb_internal_fmm** fmm);

/* Accepts NULL. */
void myb_internal_fmm_destroy(myb_internal_fmm* fmm);

/* For the tests and the tuning: the leaf size (0: from the order) and the number of point
 * pairs below which a pair of cells is summed directly (0: from the order). */
void myb_internal_fmm_tune(myb_internal_fmm* fmm, size_t leaf_size, size_t direct_pairs);

int myb_internal_fmm_order(const myb_internal_fmm* fmm);
double myb_internal_fmm_theta(const myb_internal_fmm* fmm);

/* a (3 per point) = the Newtonian accelerations of the `count` points at (px, py, pz) (one
 * array per component) from the first `massive` of them, of GM gm[0 .. massive-1] (finite,
 * >= 0); the others feel the field and attract nothing. Errors: MYB_E_DOMAIN (a position not
 * finite, the cloud too wide), MYB_E_NOMEM and MYB_E_OVERFLOW (the tree grows). Coincident
 * interacting points give non-finite accelerations, as the direct sum does. */
[[nodiscard]] myb_status myb_internal_fmm_evaluate(myb_internal_fmm* fmm, size_t massive,
                                                   const double* gm, const double* px,
                                                   const double* py, const double* pz, double* a);

myb_internal_fmm_stats myb_internal_fmm_statistics(const myb_internal_fmm* fmm);

#endif /* MYB_SRC_DYN_FMM_H_ */
