/* SPDX-License-Identifier: BSD-3-Clause */
/* glue_internal.h — what the C files of the Go binding share among themselves (not seen by
 * Go): the trampolines that call Go functions back, and small helpers. */
#ifndef MYBGO_GLUE_INTERNAL_H_
#define MYBGO_GLUE_INTERNAL_H_

#include "glue.h"

/* The context of a Go function: the handle of its runtime/cgo.Handle. */
typedef struct mybgo_go_context {
  uintptr_t handle;
} mybgo_go_context;

/* The context of a residual function of myb_shoot: its handle and the count of unknowns. */
typedef struct mybgo_go_shoot_context {
  uintptr_t handle;
  size_t n;
} mybgo_go_shoot_context;

/* Trampolines with the signatures of the library; context: const mybgo_go_context* (that of
 * mybgo_go_shoot: const mybgo_go_shoot_context*). */
myb_status mybgo_go_term(void* context, double t, const double x[3], const double v[3],
                         double a[3]);
double mybgo_go_event(void* context, double t, const double* y);
myb_status mybgo_go_action(void* context, double t, const myb_system_view* view,
                           double acceleration[3], double torque[3]);
myb_status mybgo_go_shoot(void* context, const double* u, double* f, double* jacobian);
myb_status mybgo_go_gaussian(void* context, double t, double mean[6], double covariance[36]);

/* The stations of a model's environment (and their measurements with the light time when
 * light is not NULL), and the measurement of an observation record (glue.c). */
myb_status mybgo_stations_of(mybgo_model* model, const mybgo_station* stations, size_t count,
                             const mybgo_light* light, myb_station* st, myb_light_measure* lm);
myb_status mybgo_measure_of(const mybgo_observation* o, myb_station* stations,
                            myb_light_measure* light, size_t station_count, myb_observation* out);

/* The status of a call that ended with *failed: the index passed out. */
static inline myb_status mybgo_end(myb_status status, size_t index, size_t* failed) {
  if (failed != NULL) {
    *failed = index;
  }
  return status;
}

#endif /* MYBGO_GLUE_INTERNAL_H_ */
