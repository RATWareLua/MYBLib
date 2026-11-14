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

/* As myb_internal_ias15_retake, but t_limit may also lie beyond the end of the last step, by at
 * most its own length: the step is taken again from its beginning to end exactly at t_limit,
 * from its own polynomial (continued past its end when it is lengthened), which converges in a
 * sweep or two when t_limit is close to the end. Unlike a retake, the step after it is
 * predicted as after any step, from the polynomial of the step reshaped, and it is planned as
 * the step reshaped had planned it (or shorter, if the reshaped step says so). For a step
 * aimed at an instant predicted before it was taken — the edge of a shadow, the beginning of a
 * contact — and brought to that instant once it is known better (myb_orbit). MYB_E_DOMAIN,
 * with nothing changed, if t_limit is not within those bounds or equals the end of the step;
 * otherwise the errors of myb_ias15_step. */
[[nodiscard]] myb_status myb_internal_ias15_reshape(myb_ias15* integrator, double t_limit);

/* The first `count` (<= n) coordinates and their rates at the time t from the polynomial of
 * the last accepted step, which may be carried beyond its end (an extrapolation, for a
 * prediction); before any step (just after a start) from the Taylor polynomial of the state
 * and its acceleration, x + v dt + a dt^2 / 2. y: x[0 .. count-1], then v[0 .. count-1]. */
void myb_internal_ias15_extrapolate(const myb_ias15* integrator, double t, size_t count,
                                    double* y);

/* The length of the step the integrator means to take next (unsigned; infinity when it has no
 * timescale), after a start as after a step. */
double myb_internal_ias15_next_step(const myb_ias15* integrator);

/* The timescale of the step (myb_ias15_options.timescale) of an integrator already made: for
 * the integrators made inside other objects (the variational equations of myb/uncert.h). */
void myb_internal_ias15_set_timescale(myb_ias15* integrator, myb_ias15_timescale timescale);

/* The length of the step the integrator would take next (0 before its first step), and a
 * length to take instead of the one a start has chosen: for a start that continues an
 * integration, where the automatic first step — made for a fresh problem, in which the
 * components that are zero are exactly zero — would be misled by components that merely
 * pass through zero. A step that is too long is rejected and redone as usual. */
double myb_internal_ias15_planned_step(const myb_ias15* integrator);
void myb_internal_ias15_plan_step(myb_ias15* integrator, double h);

/* An integrator of n coordinates whose last `followers` (0 <= followers < n) are carried along
 * without a say in the integration: the step (the automatic first one included) and the end of
 * the predictor-corrector sweeps are those of the first n - followers coordinates alone — the
 * groups of the options must cover exactly these, the block divide both counts. The followers are
 * the variational equations of the leading coordinates (dyn/system.c, myb_system_partials): their
 * timescale is that of the motion they belong to, and a variation whose acceleration passes
 * through zero in every component (as one of them may) cannot stall the step. The leading
 * coordinates then take the same steps, and get the same bits, as when they are integrated
 * alone. Errors: as myb_ias15_create; MYB_E_DOMAIN also for followers >= n. */
[[nodiscard]] myb_status myb_internal_ias15_create_followed(size_t n, size_t followers,
                                                            const myb_ias15_options* options,
                                                            const myb_allocator* allocator,
                                                            myb_ias15** integrator);

/* Bounds of the position coordinate `coordinate` (0 .. n-1) over the whole last step: the
 * step's own polynomial (degree 9 in the fraction of the step) written in the basis of
 * Bernstein, whose coefficients enclose it, widened by their rounding. MYB_E_RANGE before a
 * step. */
[[nodiscard]] myb_status myb_internal_ias15_position_box(const myb_ias15* integrator,
                                                         size_t coordinate, double* lo,
                                                         double* hi);

#endif /* MYB_SRC_INTEG_IAS15_INTERNAL_H_ */
