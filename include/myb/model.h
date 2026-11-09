/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/model.h — the force model of one object: a list of named acceleration terms.
 *
 * A term is a function a(t, x, v) supplied by the caller (gravity of a field, a third body, a
 * pressure, a manoeuvre, ...), usually built on myb/grav.h, myb/force.h and myb/rel.h. The
 * model sums the enabled terms in the order they were added, with compensated summation per
 * component, so the result does not depend on how terms are grouped; each term can be switched
 * off, evaluated alone (effect budgets) and differentiated (variational equations: an analytic
 * Jacobian if the term gives one, central differences otherwise). Parameters are doubles inside
 * the terms' contexts, registered by address, with respect to which the variational equations
 * also differentiate (cr, drag coefficients, empirical accelerations). */
#ifndef MYB_MODEL_H_
#define MYB_MODEL_H_

#include <stddef.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The acceleration of a term (written, not added). Any status other than MYB_OK stops the
 * evaluation and is returned. */
typedef myb_status (*myb_term_function)(void* context, double t, const double x[3],
                                        const double v[3], double a[3]);

/* Optional analytic Jacobian of a term: dadx[3 i + j] = d a_i / d x_j, dadv likewise. */
typedef myb_status (*myb_term_jacobian)(void* context, double t, const double x[3],
                                        const double v[3], double dadx[9], double dadv[9]);

typedef struct myb_force_model myb_force_model;

enum { MYB_MODEL_MAX_TERMS = 64, MYB_MODEL_MAX_PARAMETERS = 16, MYB_MODEL_NAME = 32 };

MYB_API MYB_NODISCARD myb_status myb_force_model_create(const myb_allocator* allocator,
                                                        myb_force_model** model);
MYB_API void myb_force_model_destroy(myb_force_model* model); /* accepts NULL */

/* Adds an enabled term; jacobian may be NULL; name: at most 31 characters, unique.
 * Errors: MYB_E_DOMAIN (name empty, too long or taken), MYB_E_OVERFLOW (64 terms). */
MYB_API MYB_NODISCARD myb_status myb_force_model_add(myb_force_model* model, const char* name,
                                                     myb_term_function function,
                                                     myb_term_jacobian jacobian, void* context,
                                                     size_t* index);
MYB_API MYB_NODISCARD myb_status myb_force_model_find(const myb_force_model* model,
                                                      const char* name, size_t* index);
MYB_API size_t myb_force_model_count(const myb_force_model* model);
MYB_API const char* myb_force_model_name(const myb_force_model* model, size_t index);
MYB_API int myb_force_model_enabled(const myb_force_model* model, size_t index);
MYB_API void myb_force_model_set_enabled(myb_force_model* model, size_t index, int enabled);

/* A parameter: *value (inside a term's context) and the step of its central differences.
 * Errors: MYB_E_DOMAIN (step <= 0, name invalid or taken), MYB_E_OVERFLOW (16 parameters). */
MYB_API MYB_NODISCARD myb_status myb_force_model_add_parameter(myb_force_model* model,
                                                               const char* name, double* value,
                                                               double step, size_t* index);
MYB_API size_t myb_force_model_parameter_count(const myb_force_model* model);
MYB_API const char* myb_force_model_parameter_name(const myb_force_model* model, size_t index);
MYB_API double* myb_force_model_parameter(const myb_force_model* model, size_t index);

/* Scales of the central differences: steps cbrt(eps) * scale per component (defaults: the norm
 * of x, resp. v, at the point, at least length_scale, resp. velocity_scale; both 1 by default).
 * Errors: MYB_E_DOMAIN for non-positive scales. */
MYB_API MYB_NODISCARD myb_status myb_force_model_set_scales(myb_force_model* model,
                                                            double length_scale,
                                                            double velocity_scale);

/* Sum of the enabled terms. */
MYB_API MYB_NODISCARD myb_status myb_force_model_acceleration(myb_force_model* model, double t,
                                                              const double x[3],
                                                              const double v[3], double a[3]);
/* One term, enabled or not. */
MYB_API MYB_NODISCARD myb_status myb_force_model_term(myb_force_model* model, size_t index,
                                                      double t, const double x[3],
                                                      const double v[3], double a[3]);

/* Jacobian of the sum of the enabled terms: d a / d x, d a / d v (3 x 3, row-major) and, if
 * dadp is not NULL, d a / d p (3 x parameter_count). Central differences for the terms without
 * an analytic Jacobian (relative error about 1e-10) and for the parameters. */
MYB_API MYB_NODISCARD myb_status myb_force_model_jacobian(myb_force_model* model, double t,
                                                          const double x[3], const double v[3],
                                                          double dadx[9], double dadv[9],
                                                          double* dadp);

/* Adapter to the integrators (myb_accel_function of myb/integ.h, n = 3): context = model. */
MYB_API MYB_NODISCARD myb_status myb_force_model_accel(void* model, double t, const double* x,
                                                       const double* v, double* a);

/* Analytic Jacobian of a point mass: a = -gm r / |r|^3, d a / d r = gm (3 r r^T / |r|^5 -
 * 1 / |r|^3); r is the position relative to the mass. */
MYB_API void myb_point_mass_jacobian(double gm, const double r[3], double dadx[9]);

#ifdef __cplusplus
}
#endif

#endif /* MYB_MODEL_H_ */
