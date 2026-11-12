/* SPDX-License-Identifier: BSD-3-Clause */
/* integ/ias15_internal.h — taking the last step of IAS15 again, shorter: for a force that is
 * not smooth at an instant inside the step (myb_force_model_add_break). */
#ifndef MYB_SRC_INTEG_IAS15_INTERNAL_H_
#define MYB_SRC_INTEG_IAS15_INTERNAL_H_

#include "myb/integ.h"

/* Withdraws the last accepted step and takes it again from its beginning up to t_limit, which
 * must lie strictly inside it: one step that ends exactly at t_limit, started from the
 * polynomial of the step withdrawn (the same one, on the part kept), so that it converges in
 * a sweep or two. The step after it starts from that polynomial too, continued past t_limit,
 * and is at most as long as the step withdrawn. MYB_E_DOMAIN, with nothing changed, if
 * t_limit is not inside the last step; otherwise the errors of myb_ias15_step. */
[[nodiscard]] myb_status myb_internal_ias15_retake(myb_ias15* integrator, double t_limit);

/* The length of the step the integrator would take next (0 before its first step), and a
 * length to take instead of the one a start has chosen: for a start that continues an
 * integration, where the automatic first step — made for a fresh problem, in which the
 * components that are zero are exactly zero — would be misled by components that merely
 * pass through zero. A step that is too long is rejected and redone as usual. */
double myb_internal_ias15_planned_step(const myb_ias15* integrator);
void myb_internal_ias15_plan_step(myb_ias15* integrator, double h);

/* Bounds of the position coordinate `coordinate` (0 .. n-1) over the whole last step: the
 * step's own polynomial (degree 9 in the fraction of the step) written in the basis of
 * Bernstein, whose coefficients enclose it, widened by their rounding. MYB_E_RANGE before a
 * step. */
[[nodiscard]] myb_status myb_internal_ias15_position_box(const myb_ias15* integrator,
                                                         size_t coordinate, double* lo,
                                                         double* hi);

#endif /* MYB_SRC_INTEG_IAS15_INTERNAL_H_ */
