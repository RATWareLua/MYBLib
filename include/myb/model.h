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

MYB_NODISCARD MYB_API myb_status myb_force_model_create(const myb_allocator* allocator,
                                                        myb_force_model** model);
MYB_API void myb_force_model_destroy(myb_force_model* model); /* accepts NULL */

/* Adds an enabled term; jacobian may be NULL; name: at most 31 characters, unique.
 * Errors: MYB_E_DOMAIN (name empty, too long or taken), MYB_E_OVERFLOW (64 terms). */
MYB_NODISCARD MYB_API myb_status myb_force_model_add(myb_force_model* model, const char* name,
                                                     myb_term_function function,
                                                     myb_term_jacobian jacobian, void* context,
                                                     size_t* index);
MYB_NODISCARD MYB_API myb_status myb_force_model_find(const myb_force_model* model,
                                                      const char* name, size_t* index);
MYB_API size_t myb_force_model_count(const myb_force_model* model);

/* Declares term `index` a function of the time and the place alone, not of the velocity (a
 * term is taken to depend on it until declared otherwise): a model whose enabled terms all are
 * may be integrated by RKN12 (myb/orbit.h). Of the terms of myb/env.h: myb_term_point_mass,
 * myb_term_third_body, myb_term_field, myb_term_field_indirect, and myb_term_radiation without
 * `moving`. The declaration is the caller's; it is not checked. */
MYB_API void myb_force_model_set_position_only(myb_force_model* model, size_t index,
                                               int position_only);
/* 1 if every enabled term is declared a function of the place alone. */
MYB_API int myb_force_model_position_only(const myb_force_model* model);
MYB_API const char* myb_force_model_name(const myb_force_model* model, size_t index);
MYB_API int myb_force_model_enabled(const myb_force_model* model, size_t index);
MYB_API void myb_force_model_set_enabled(myb_force_model* model, size_t index, int enabled);

/* A parameter: *value (inside a term's context) and the step of its central differences.
 * Errors: MYB_E_DOMAIN (step <= 0, name invalid or taken), MYB_E_OVERFLOW (16 parameters). */
MYB_NODISCARD MYB_API myb_status myb_force_model_add_parameter(myb_force_model* model,
                                                               const char* name, double* value,
                                                               double step, size_t* index);
MYB_API size_t myb_force_model_parameter_count(const myb_force_model* model);
MYB_API const char* myb_force_model_parameter_name(const myb_force_model* model, size_t index);
MYB_API double* myb_force_model_parameter(const myb_force_model* model, size_t index);

/* Binds a parameter to the one term it enters: its partial derivative is then the difference
 * of that term alone instead of the whole model (the same value, a fraction of the cost), and
 * zero while the term is switched off. Errors: MYB_E_DOMAIN (indices). */
MYB_NODISCARD MYB_API myb_status myb_force_model_bind_parameter(myb_force_model* model,
                                                                size_t parameter, size_t term);

/* Where the force is continuous but not smooth — the edge of a shadow, a thrust that switches
 * — a step of an integrator that straddles the instant loses its order: the error no longer
 * falls with the precision asked for. A break is a function g(t, y), y = (x, v), that changes
 * sign there (the form of an event function of myb/integ.h). myb_orbit (myb/orbit.h), and with
 * it the estimators, ends a step at every sign change of a break of the model, so that no
 * step contains one. The context must outlive the model.
 * Errors: MYB_E_OVERFLOW (more than MYB_MODEL_MAX_BREAKS). */
typedef double (*myb_break_function)(void* context, double t, const double* y);

enum { MYB_MODEL_MAX_BREAKS = 32 };

MYB_NODISCARD MYB_API myb_status myb_force_model_add_break(myb_force_model* model,
                                                           myb_break_function g, void* context);
MYB_API size_t myb_force_model_break_count(const myb_force_model* model);
MYB_API void myb_force_model_break(const myb_force_model* model, size_t index,
                                   myb_break_function* g, void** context);

/* How myb_orbit with IAS15 (and the estimators through it) ends its steps at the breaks.
 *   MYB_BREAKS_RETAKE   (the default) a step that turns out to contain a sign change is
 *                       withdrawn and taken again from its beginning up to it. Exact, but the
 *                       step that ran past the break is work thrown away, and the sweeps of a
 *                       step that holds a kink converge badly. A step that cannot be taken at
 *                       all across the break — a force so stiff beyond it (a contact) that the
 *                       sweeps of a long step diverge — is shortened until it can, which may
 *                       take many attempts.
 *   MYB_BREAKS_PREDICT  before a step the break functions are examined on the polynomial of
 *                       the last step carried over the step to come (on the Taylor polynomial
 *                       of the state after a start); a sign change predicted there becomes the
 *                       end of the step. Once taken, the step's own polynomial — now on one
 *                       side of the break — locates it again: inside the step or just beyond
 *                       its end. The step is then reshaped to end there (taken again from its
 *                       beginning, warm, in a sweep or two), until the root of its own
 *                       polynomial is at its end to 2^-30 of the step; no step is ever taken
 *                       across a break and thrown away. The step after the break is planned at
 *                       most a quarter of the one planned before it (the force starts anew
 *                       there, not smoothly, and nothing in the polynomial before the break
 *                       knows of it), and the step after any aimed one at most 4 times it (it
 *                       is predicted from it). A break the prediction misses is met as by
 *                       MYB_BREAKS_RETAKE.
 * Both give the same results to the accuracy of the integration. Measured on an orbit of
 * 12 000 km with 0.2 m^2/kg through 54 edges of the Earth's shadow in two days
 * (tests/unit/test_events.c): the final place over epsilon 1e-9 .. 1e-12 within 8.2e-7 m
 * withdrawing, 8.4e-7 m predicting (4e-7 m apart at 1e-13); 517, 701, 953, 1308 steps against
 * 557, 742, 993, 1362 — a step more an edge, the quarter step after it. Without the quarter the
 * prediction lost three orders (1.2e-4 m at 1e-9): the predicted edges themselves were right
 * to 3e-7 s, the step after them too long. On the fit of LAGEOS aiming saves 1 to 5 % of the
 * evaluations of the model and no time that shows (myb/orbit.h): RETAKE stays the default.
 * DOP853 and RKN12 restart at a break in either case. */
typedef enum myb_break_scheme {
  MYB_BREAKS_RETAKE = 0,
  MYB_BREAKS_PREDICT = 1
} myb_break_scheme;

/* Takes effect from the next myb_orbit_start. An unknown scheme is a contract violation. */
MYB_API void myb_force_model_set_break_scheme(myb_force_model* model, myb_break_scheme scheme);
MYB_API myb_break_scheme myb_force_model_break_scheme(const myb_force_model* model);

/* Scales of the central differences: steps cbrt(eps) * scale per component (defaults: the norm
 * of x, resp. v, at the point, at least length_scale, resp. velocity_scale; both 1 by default).
 * Errors: MYB_E_DOMAIN for non-positive scales. */
MYB_NODISCARD MYB_API myb_status myb_force_model_set_scales(myb_force_model* model,
                                                            double length_scale,
                                                            double velocity_scale);

/* Sum of the enabled terms. */
MYB_NODISCARD MYB_API myb_status myb_force_model_acceleration(myb_force_model* model, double t,
                                                              const double x[3],
                                                              const double v[3], double a[3]);
/* One term, enabled or not. */
MYB_NODISCARD MYB_API myb_status myb_force_model_term(myb_force_model* model, size_t index,
                                                      double t, const double x[3],
                                                      const double v[3], double a[3]);

/* Jacobian of the sum of the enabled terms: d a / d x, d a / d v (3 x 3, row-major) and, if
 * dadp is not NULL, d a / d p (3 x parameter_count). Central differences for the terms without
 * an analytic Jacobian (relative error about 1e-10) and for the parameters. */
MYB_NODISCARD MYB_API myb_status myb_force_model_jacobian(myb_force_model* model, double t,
                                                          const double x[3], const double v[3],
                                                          double dadx[9], double dadv[9],
                                                          double* dadp);

/* Adapter to the integrators (myb_accel_function of myb/integ.h, n = 3): context = model. */
MYB_NODISCARD MYB_API myb_status myb_force_model_accel(void* model, double t, const double* x,
                                                       const double* v, double* a);

/* A Jacobian of zeros (myb_term_jacobian): for a term whose derivatives are left out of the
 * variational equations on purpose — a force so far below the others that its gradient does
 * not show in the partials, and whose differences would cost more than they bring. */
MYB_NODISCARD MYB_API myb_status myb_term_no_jacobian(void* context, double t, const double x[3],
                                                      const double v[3], double dadx[9],
                                                      double dadv[9]);

/* Analytic Jacobian of a point mass: a = -gm r / |r|^3, d a / d r = gm (3 r r^T / |r|^5 -
 * 1 / |r|^3); r is the position relative to the mass. */
MYB_API void myb_point_mass_jacobian(double gm, const double r[3], double dadx[9]);

#ifdef __cplusplus
}
#endif

#endif /* MYB_MODEL_H_ */
