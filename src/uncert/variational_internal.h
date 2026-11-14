/* SPDX-License-Identifier: BSD-3-Clause */
/* uncert/variational_internal.h — as integ/ias15_internal.h for the variational equations. */
#ifndef MYB_SRC_UNCERT_VARIATIONAL_INTERNAL_H_
#define MYB_SRC_UNCERT_VARIATIONAL_INTERNAL_H_

#include "myb/uncert.h"

[[nodiscard]] myb_status myb_internal_variational_retake(myb_variational* variational,
                                                         double t_limit);
[[nodiscard]] myb_status myb_internal_variational_reshape(myb_variational* variational,
                                                          double t_limit);
/* The state (6 values) at t from the polynomial of the last step carried as far as asked
 * (myb_internal_ias15_extrapolate), and the length of the next step. */
void myb_internal_variational_extrapolate(const myb_variational* variational, double t,
                                          double y[6]);
double myb_internal_variational_next_step(const myb_variational* variational);
void myb_internal_variational_plan_step(myb_variational* variational, double h);
void myb_internal_variational_set_timescale(myb_variational* variational,
                                            myb_ias15_timescale timescale);

/* Starts again at t from the state y with the transition matrix and the sensitivities given
 * (instead of the identity and zeros): the equations are linear in them, so they continue as
 * if the integration had not stopped. */
[[nodiscard]] myb_status myb_internal_variational_restart(myb_variational* variational, double t,
                                                          const double y[6],
                                                          const double phi[36],
                                                          const double* sensitivity);

#endif /* MYB_SRC_UNCERT_VARIATIONAL_INTERNAL_H_ */
