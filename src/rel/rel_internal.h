/* SPDX-License-Identifier: BSD-3-Clause */
/* rel/rel_internal.h — the post-Newtonian terms of a binary apart from the Newtonian one, for
 * a pair inside a system whose point masses already carry it (dyn/system.c), and the
 * interface between rel.c and the higher terms of rel/pn_higher.c. */
#ifndef MYB_SRC_REL_REL_INTERNAL_H_
#define MYB_SRC_REL_REL_INTERNAL_H_

#include "myb/rel.h"

/* Every flag of myb_pn_binary.orders, and those that need the spins. */
#define MYB_PN_KNOWN                                                                    \
  ((unsigned)MYB_PN_ALL | (unsigned)MYB_PN_4 | (unsigned)MYB_PN_4RR |                    \
   (unsigned)MYB_PN_SO25 | (unsigned)MYB_PN_SO35 | (unsigned)MYB_PN_SORR)
#define MYB_PN_SPIN_TERMS                                                              \
  ((unsigned)MYB_PN_SO | (unsigned)MYB_PN_SS | (unsigned)MYB_PN_SO25 |                   \
   (unsigned)MYB_PN_SO35 | (unsigned)MYB_PN_SORR)
/* The flags whose terms are computed by rel/pn_higher.c */
#define MYB_PN_HIGHER                                                                  \
  ((unsigned)MYB_PN_4 | (unsigned)MYB_PN_4RR | (unsigned)MYB_PN_SO25 |                    \
   (unsigned)MYB_PN_SO35 | (unsigned)MYB_PN_SORR)

/* As myb_pn_binary_evaluate, without the term -G m x / r^3: the terms named by b->orders
 * only, to their own accuracy (not a difference of two nearly equal accelerations). */
[[nodiscard]] myb_status myb_internal_pn_binary_terms(const myb_pn_binary* b, const double x[3],
                                                      const double v[3], const double chi1[3],
                                                      const double chi2[3], double a[3],
                                                      double dchi1[3], double dchi2[3]);

/* A_4PN and B_4PN of the instantaneous 4PN acceleration (rel/pn_higher.c), before the
 * division by c^8: u = G m / r, rd = n.v, v2 = v.v (SI), log_r = ln(r / r0). */
void myb_internal_pn4_coefficients(double nu, double u, double rd, double v2, double log_r,
                                   double* a4, double* b4);

/* The terms of MYB_PN_4 (its tail; A_4PN and B_4PN come from the function above), MYB_PN_4RR,
 * MYB_PN_SO25, MYB_PN_SO35 and MYB_PN_SORR named by b->orders: accelerations added to a (m/s^2),
 * precession angular velocities of the spins to w1, w2 (1/s). The state is validated by the
 * caller. MYB_E_DOMAIN when a tail is asked for and the osculating Newtonian orbit is unbound
 * or more eccentric than about 0.83 (more than 1024 harmonics); nothing is added then. */
[[nodiscard]] myb_status myb_internal_pn_higher(const myb_pn_binary* b, const double x[3],
                                                const double v[3], const double chi1[3],
                                                const double chi2[3], double a[3],
                                                double w1[3], double w2[3]);

#endif /* MYB_SRC_REL_REL_INTERNAL_H_ */
