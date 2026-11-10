/* SPDX-License-Identifier: BSD-3-Clause */
/* env/env_internal.h — what the terms read from an environment without copies. The pointers
 * stay valid until the environment is asked about another instant. */
#ifndef MYB_SRC_ENV_ENV_INTERNAL_H_
#define MYB_SRC_ENV_ENV_INTERNAL_H_

#include <stddef.h>

#include "myb/env.h"

/* The description of a body. */
const myb_env_body* myb_internal_env_body(const myb_env* env, size_t body);

/* Index of the body at the origin, or -1 for an inertial origin. */
int myb_internal_env_origin(const myb_env* env);

/* State of a body relative to the origin at the model time t (zeros for the body at the
 * origin). */
[[nodiscard]] myb_status myb_internal_env_state(myb_env* env, double t, size_t body,
                                                const double** state);

/* Orientation of a body at t and the inverse of its rotation (inverse may be NULL). */
[[nodiscard]] myb_status myb_internal_env_orientation(myb_env* env, double t, size_t body,
                                                      const myb_orientation** orientation,
                                                      const myb_rotation** inverse);

/* The instant t on the scale of the body's orientation provider. */
[[nodiscard]] myb_status myb_internal_env_orientation_time(myb_env* env, double t, size_t body,
                                                           myb_time* out);

/* The light of a source at x: d = x - (the source), its length and the flux (W m^-2) left by
 * the shadows of the occulters (as myb_term_radiation). */
[[nodiscard]] myb_status myb_internal_env_light(myb_env* env, double t, size_t source,
                                                const size_t* occulters, size_t occulter_count,
                                                const myb_refraction* atmospheres,
                                                const double x[3], double d[3],
                                                double* distance, double* flux);

/* The mass of the object at t (its function, if it has one); 0 if it is not valid. */
double myb_internal_object_mass(const myb_object* object, double t);

#endif /* MYB_SRC_ENV_ENV_INTERNAL_H_ */
