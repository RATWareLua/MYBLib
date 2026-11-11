/* SPDX-License-Identifier: BSD-3-Clause */
/* rel/rel_internal.h — the post-Newtonian terms of a binary apart from the Newtonian one, for
 * a pair inside a system whose point masses already carry it (dyn/system.c). */
#ifndef MYB_SRC_REL_REL_INTERNAL_H_
#define MYB_SRC_REL_REL_INTERNAL_H_

#include "myb/rel.h"

/* As myb_pn_binary_evaluate, without the term -G m x / r^3: the terms named by b->orders
 * only, to their own accuracy (not a difference of two nearly equal accelerations). */
[[nodiscard]] myb_status myb_internal_pn_binary_terms(const myb_pn_binary* b, const double x[3],
                                                      const double v[3], const double chi1[3],
                                                      const double chi2[3], double a[3],
                                                      double dchi1[3], double dchi2[3]);

#endif /* MYB_SRC_REL_REL_INTERNAL_H_ */
