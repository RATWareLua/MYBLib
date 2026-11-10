/* SPDX-License-Identifier: BSD-3-Clause */
/* uncert/variational_internal.h — as integ/ias15_internal.h for the variational equations. */
#ifndef MYB_SRC_UNCERT_VARIATIONAL_INTERNAL_H_
#define MYB_SRC_UNCERT_VARIATIONAL_INTERNAL_H_

#include "myb/uncert.h"

[[nodiscard]] myb_status myb_internal_variational_retake(myb_variational* variational,
                                                         double t_limit);

/* Starts again at t from the state y with the transition matrix and the sensitivities given
 * (instead of the identity and zeros): the equations are linear in them, so they continue as
 * if the integration had not stopped. */
[[nodiscard]] myb_status myb_internal_variational_restart(myb_variational* variational, double t,
                                                          const double y[6],
                                                          const double phi[36],
                                                          const double* sensitivity);

#endif /* MYB_SRC_UNCERT_VARIATIONAL_INTERNAL_H_ */
