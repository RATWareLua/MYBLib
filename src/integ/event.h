/* SPDX-License-Identifier: BSD-3-Clause */
/* integ/event.h — root of an event function inside one integration step, located on the
 * integrator's dense output by the Illinois variant of regula falsi (M. Dowell, P. Jarratt,
 * BIT 11, 168, 1971), with a bisection every third iteration.
 *
 * The root is reported when g(t0) and g(t1) have strictly opposite signs, or when g(t1) == 0
 * (a zero at t0 belongs to the previous step). The result is the bracket end with the smaller
 * |g| once the bracket reaches the resolution of t. */
#ifndef MYB_SRC_INTEG_EVENT_H_
#define MYB_SRC_INTEG_EVENT_H_

#include <stddef.h>

#include "myb/integ.h"

/* State of the integrator at time t (within the step), m values. */
typedef myb_status (*myb_internal_dense_function)(void* integrator, double t, double* y);

typedef struct myb_internal_step_view {
  myb_internal_dense_function dense;
  void* integrator;
  size_t m;          /* state length */
  double t0;         /* start of the step */
  const double* y0;  /* state at t0 */
  double t1;         /* end of the step */
  const double* y1;  /* state at t1 */
  double* scratch;   /* m values */
} myb_internal_step_view;

/* g must not call back into the integrator. y_event may be NULL. */
myb_status myb_internal_locate_event(const myb_internal_step_view* step, myb_event_function g,
                                     void* context, int* found, double* t_event,
                                     double* y_event);

#endif /* MYB_SRC_INTEG_EVENT_H_ */
