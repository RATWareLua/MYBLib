/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/orbit.h — propagation of one object under a force model: its state at the times asked
 * for, the events along the way and, on request, the state transition matrix and the
 * sensitivities to the parameters of the model.
 *
 * The object is integrated by IAS15 (myb/integ.h) from t0 towards t_end with the integrator's
 * own steps, and every state comes from the dense output of the step that contains its time.
 * The steps never depend on the times asked for (only the last one is cut, at t_end), so the
 * state at a time t is the same, bit for bit, whatever else was asked — sampled every second
 * or not at all. Times are those of the force model (seconds from the epoch of its
 * environment, myb/env.h) and must be asked in the order of the propagation, forward or
 * backward; starting again restarts it.
 *
 * Where the model declares breaks (myb_force_model_add_break: the edges of shadows) a step
 * ends at each of them: a step that would contain one is taken again up to it. The result
 * then does not depend on where the steps happen to fall around the edge, and its error falls
 * with the precision asked for as it does for a smooth force. The transition matrix and the
 * sensitivities continue across. */
#ifndef MYB_ORBIT_H_
#define MYB_ORBIT_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"
#include "myb/integ.h"
#include "myb/model.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct myb_orbit myb_orbit;

typedef enum myb_orbit_integrator {
  MYB_ORBIT_IAS15 = 0, /* the default */
  MYB_ORBIT_DOP853 = 1, /* the explicit Runge-Kutta pair of myb/integ.h: cheaper per step at
                           moderate accuracies, for forces of any dependence on the velocity.
                           With variational != 0 the variations are carried as part of one
                           system of the first order, under the same tolerances as the state */
  MYB_ORBIT_RKN12 = 2   /* the Runge-Kutta-Nystrom pair RKN12(10) of myb/integ.h with its
                           interpolant: for very tight tolerances, forces that do not depend on
                           the velocity only — every enabled term of the model declared so
                           (myb_force_model_set_position_only), else MYB_E_DOMAIN at the start;
                           no variational equations */
} myb_orbit_integrator;

typedef struct myb_orbit_options {
  double epsilon;    /* IAS15 precision parameter; 0 = 1e-10 */
  double h_max;      /* largest step (s); 0 = unlimited */
  int64_t max_steps; /* budget of steps of one call; 0 = 10^9 */
  int variational;   /* != 0: the variational equations are integrated with the state
                        (myb/uncert.h), for myb_orbit_partials */
  myb_orbit_integrator integrator;
  double rtol; /* DOP853 and RKN12: relative tolerance; 0 = 1e-12 */
  double atol; /* DOP853 and RKN12: absolute tolerance (m and m/s); 0 = 1e-9 */
} myb_orbit_options;

/* The model must outlive the object; with variational != 0 its parameters at this moment are
 * those differentiated. options may be NULL (defaults). Errors: MYB_E_DOMAIN (negative
 * options, an integrator that is not one of the three, RKN12 with variational), MYB_E_NOMEM,
 * MYB_E_OVERFLOW. */
MYB_API MYB_NODISCARD myb_status myb_orbit_create(myb_force_model* model,
                                                  const myb_orbit_options* options,
                                                  const myb_allocator* allocator,
                                                  myb_orbit** orbit);
MYB_API void myb_orbit_destroy(myb_orbit* orbit); /* accepts NULL */

/* Starts at t0 from y0 = (x, v) towards t_end != t0, which sets the direction and is never
 * passed. Errors: MYB_E_DOMAIN (times or state not finite, t_end == t0; with RKN12 a model with
 * a term not declared a function of the place alone), the statuses of the model. */
MYB_API MYB_NODISCARD myb_status myb_orbit_start(myb_orbit* orbit, double t0, const double y0[6],
                                                 double t_end);

/* The time reached: that of the last state or event returned (t0 after the start). */
MYB_API double myb_orbit_time(const myb_orbit* orbit);

/* The state at t, between the time reached and t_end. Errors: MYB_E_DOMAIN (t behind the time
 * reached, beyond t_end or not finite), MYB_E_NOCONV (step budget), the statuses of the
 * integrator and the model. */
MYB_API MYB_NODISCARD myb_status myb_orbit_state(myb_orbit* orbit, double t, double y[6]);

/* The states at count times in the order of the propagation: states[6 k .. 6 k + 5] at
 * times[k]. Errors: as myb_orbit_state (the states before the failing time are filled). */
MYB_API MYB_NODISCARD myb_status myb_orbit_states(myb_orbit* orbit, const double* times,
                                                  size_t count, double* states);

/* The state at t with phi[6 i + j] = d y_i / d y0_j and sensitivity[i np + k] = d y_i / d p_k
 * for the np parameters of the model (any output may be NULL). Errors: as myb_orbit_state;
 * MYB_E_DOMAIN for an object created without variational. */
MYB_API MYB_NODISCARD myb_status myb_orbit_partials(myb_orbit* orbit, double t, double y[6],
                                                    double phi[36], double* sensitivity);

/* The next root of g(t, y), y = (x, v), from the time reached up to t_limit (between it and
 * t_end): the sign of g is examined at the ends of four equal parts of every step, and a
 * change is located on the dense output to the resolution of t (as myb_ias15_locate), so
 * roots closer together than a quarter of a step may be missed. direction: +1 only roots
 * where g grows with t, -1 only where it falls, 0 both. With *found != 0 the time reached
 * becomes the root (t_event; y_event may be NULL) — a root at the time reached itself is not
 * reported again; otherwise it becomes t_limit. Errors: as myb_orbit_state. */
MYB_API MYB_NODISCARD myb_status myb_orbit_event(myb_orbit* orbit, myb_event_function g,
                                                 void* context, double t_limit, int direction,
                                                 int* found, double* t_event, double y_event[6]);

/* An impulsive manoeuvre: at the time reached the velocity changes by dv (m/s, in the
 * inertial axes of the model), and the propagation goes on from there towards the same
 * t_end. The transition matrix and the sensitivities continue: they are those of the
 * trajectory with the impulse, the impulse itself taken as given. A state asked at the time
 * of the impulse is the one after it. Errors: MYB_E_DOMAIN (dv not finite), the statuses of
 * the model. */
MYB_API MYB_NODISCARD myb_status myb_orbit_kick(myb_orbit* orbit, const double dv[3]);

/* Accepted steps since the start, and the breaks of the model met among them. */
MYB_API int64_t myb_orbit_steps(const myb_orbit* orbit);
MYB_API int64_t myb_orbit_breaks(const myb_orbit* orbit);

#ifdef __cplusplus
}
#endif

#endif /* MYB_ORBIT_H_ */
