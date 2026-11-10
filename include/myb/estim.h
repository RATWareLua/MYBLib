/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/estim.h — estimation of an orbit from observations: weighted least squares over an arc
 * and a sequential filter, on the variational equations of a force model (myb/orbit.h,
 * myb/uncert.h), with ready-made measurements.
 *
 * The unknowns are the state y0 = (x, v) at the epoch and the parameters of the force model
 * (myb_force_model_add_parameter), in that order: 6 + np numbers. A measurement is a function
 * of the state of the object at the time of the observation, with its derivatives; the
 * derivatives with respect to the unknowns follow from the state transition matrix and the
 * sensitivities. Measurements here are geometric (the state at the instant of the
 * observation): the light time between the object and the observer is not part of them. */
#ifndef MYB_ESTIM_H_
#define MYB_ESTIM_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"
#include "myb/env.h"
#include "myb/frames.h"
#include "myb/model.h"
#include "myb/surface.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Measurements ---------------------------------------------------------------------------- */

enum { MYB_MEASURE_MAX = 6, MYB_FIT_MAX_BIASES = 32 };

/* `count` values computed from the state y = (x, v) of the object at the model time t, and
 * their derivatives partials[6 i + j] = d value_i / d y_j (partials may be NULL). The count is
 * that of the observation the function serves. Any status other than MYB_OK stops the
 * estimation and is returned. */
typedef myb_status (*myb_measure_function)(void* context, double t, const double y[6],
                                           double* values, double* partials);

typedef struct myb_observation {
  double t; /* model time, s */
  myb_measure_function measure;
  void* context;
  size_t count;                  /* values, 1 .. MYB_MEASURE_MAX */
  double value[MYB_MEASURE_MAX]; /* observed */
  double sigma[MYB_MEASURE_MAX]; /* standard deviations, > 0 */
  unsigned periodic;             /* bit i: value i is an angle, its residual is reduced to
                                    (-pi, pi] */
  int bias;                      /* 0: none; k > 0: what is computed for value i is taken with
                                    the systematic shift number k - 1 + i of the estimator
                                    added (the delay of a station, the offset of a clock, of
                                    a pointing), an unknown estimated with the orbit */
} myb_observation;

/* The position (3 values) and the whole state (6 values) of the object itself, e.g. an
 * orbit to be fitted by the model; context is not used. */
MYB_API MYB_NODISCARD myb_status myb_measure_position(void* context, double t, const double y[6],
                                                      double* values, double* partials);
MYB_API MYB_NODISCARD myb_status myb_measure_state(void* context, double t, const double y[6],
                                                   double* values, double* partials);

/* An observer fixed on a body of an environment: its body-fixed position and the axes east,
 * north, up there (rows of enu, in the body-fixed frame). */
typedef struct myb_station {
  myb_env* env;
  size_t body;
  double position[3];
  double enu[9];
} myb_station;

/* A station at geodetic longitude, latitude (rad) and height (m) on an ellipsoid of the body.
 * Errors: those of myb_geodetic_to_cartesian. */
MYB_API MYB_NODISCARD myb_status myb_station_init(myb_station* station, myb_env* env, size_t body,
                                                  const myb_ellipsoid* ellipsoid,
                                                  double longitude, double latitude,
                                                  double height);

/* State of the station at the model time t relative to the origin of the environment, in its
 * inertial axes (the body's state and orientation). Errors: those of the environment. */
MYB_API MYB_NODISCARD myb_status myb_station_state(const myb_station* station, double t,
                                                   double state[6]);

/* Measurements from a station (context: const myb_station*), all of the vector d from the
 * station to the object and its rate:
 *   range       |d| (m), 1 value;
 *   range rate  d |d| / dt (m/s), 1 value;
 *   radec       right ascension in [0, 2 pi) and declination of d in the inertial axes, 2
 *               values (periodic = 1 for the first);
 *   azel        azimuth from north towards east in [0, 2 pi) and elevation above the
 *               station's horizon, 2 values (periodic = 1).
 * Errors: MYB_E_DOMAIN (the object at the station; for the angles also exactly at the pole of
 * their axes), the statuses of the environment. */
MYB_API MYB_NODISCARD myb_status myb_measure_range(void* station, double t, const double y[6],
                                                   double* values, double* partials);
MYB_API MYB_NODISCARD myb_status myb_measure_range_rate(void* station, double t,
                                                        const double y[6], double* values,
                                                        double* partials);
MYB_API MYB_NODISCARD myb_status myb_measure_radec(void* station, double t, const double y[6],
                                                   double* values, double* partials);
MYB_API MYB_NODISCARD myb_status myb_measure_azel(void* station, double t, const double y[6],
                                                  double* values, double* partials);

/* A measurement of a point of the object that is not its centre of mass — a retroreflector,
 * the phase centre of an antenna: the measurement `measure` (with its context) is taken at
 * the state of that point, offset (m, in the body frame) from the centre of mass and turned
 * by the attitude of the object; its velocity includes the turning of the body (the attitude
 * a second apart). The derivatives are those of the measurement at the point: how the
 * attitude follows the orbit — a part in |offset| / |x| — is left out. */
typedef struct myb_offset_measure {
  myb_measure_function measure;
  void* context;
  myb_attitude_function attitude;
  void* attitude_context;
  double offset[3];
} myb_offset_measure;

MYB_API MYB_NODISCARD myb_status myb_measure_offset(void* offset, double t, const double y[6],
                                                    double* values, double* partials);

/* ---- Least squares over an arc ---------------------------------------------------------------
 * Gauss-Newton on the weighted residuals (observed - computed) / sigma: every iteration
 * propagates the orbit with its variational equations through the observations, forms the
 * normal equations of the free unknowns (scaled by their diagonal, solved by Cholesky) and
 * corrects y0 and the parameters in place. It stops when the weighted RMS changes by less than
 * `tolerance` of itself, or after max_iterations corrections. */

typedef struct myb_fit_options {
  int max_iterations;          /* corrections; 0 = 10 */
  double tolerance;            /* 0 = 1e-6 */
  double epsilon;              /* IAS15 precision parameter of the propagation; 0 = 1e-10 */
  const unsigned char* fixed;  /* 6 + np flags, non-zero: the unknown keeps its value; NULL:
                                  all are estimated */
  const double* prior_sigma;   /* 6 + np standard deviations of the initial values as a priori
                                  information (0: none for that unknown); NULL: none */
  size_t bias_count;           /* systematic shifts of the measurements (myb_observation.bias),
                                  at most MYB_FIT_MAX_BIASES: unknowns after the parameters,
                                  so that fixed, prior_sigma and the covariance then have
                                  6 + np + bias_count entries a side */
  double* bias;                /* bias_count values: in the first guess, out the estimate */
} myb_fit_options;

typedef struct myb_fit_result {
  int iterations;    /* corrections applied */
  int converged;
  size_t values;     /* scalar measurements */
  double rms_before; /* weighted RMS at the initial values */
  double rms;        /* weighted RMS at the final values (1 for errors that match the sigmas) */
} myb_fit_result;

/* Fits y0 (in: the first guess at t0; out: the estimate) and the parameters of the model (in
 * place) to the observations, which must be ordered in time away from t0, all on one side of
 * it. covariance (may be NULL): (6 + np + bias_count)^2 values, the inverse of the normal
 * matrix at the estimate (rows and columns of fixed unknowns are zero). residuals (may be NULL): observed -
 * computed of every scalar measurement in the order of the observations, at the estimate.
 * Errors: MYB_E_DOMAIN (no observations, an invalid one, times out of order, a normal matrix
 * that is not positive definite: the observations do not determine the unknowns),
 * MYB_E_NOMEM, the statuses of the propagation and of the measurements. On an error y0 and
 * the parameters hold the values of the last iteration. */
MYB_API MYB_NODISCARD myb_status myb_fit_orbit(myb_force_model* model, double t0, double y0[6],
                                               const myb_observation* observations,
                                               size_t count, const myb_fit_options* options,
                                               const myb_allocator* allocator,
                                               myb_fit_result* result, double* covariance,
                                               double* residuals);

/* ---- Sequential filter ------------------------------------------------------------------------
 * The extended Kalman filter on the same unknowns z = (y, p): between observations the state
 * is propagated by the model and the covariance by its transition matrix,
 *   P <- F P F^T + Q,  F = | Phi  S |,  Q = q | dt^3/3  dt^2/2 |  (in the block of y, for an
 *                          |  0   1 |          | dt^2/2  dt     |   unmodelled white
 * acceleration of spectral density q, m^2 s^-3, per axis), and each observation corrects both
 * (K = P H^T (H P H^T + R)^-1; the covariance in Joseph's form, which stays symmetric and
 * positive). The parameters of the model are corrected in place. */

typedef struct myb_filter myb_filter;

typedef struct myb_filter_options {
  double epsilon;       /* IAS15 precision parameter; 0 = 1e-10 */
  double process_noise; /* q, >= 0 */
  size_t bias_count;    /* systematic shifts of the measurements, unknowns after the parameters
                           (they start at zero: myb_filter_set_bias) */
  size_t history;       /* instants to keep for the smoother (the start and every later time
                           of an observation or a prediction); 0: none are kept */
} myb_filter_options;

/* The model must outlive the filter. Errors: MYB_E_DOMAIN (negative options), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_filter_create(myb_force_model* model,
                                                   const myb_filter_options* options,
                                                   const myb_allocator* allocator,
                                                   myb_filter** filter);
MYB_API void myb_filter_destroy(myb_filter* filter); /* accepts NULL */

/* Starts at t0 with the state y0, the present values of the parameters and the covariance of
 * the 6 + np unknowns (row-major, symmetric); observations follow in time order towards
 * t_end, which is never passed. Errors: MYB_E_DOMAIN (times, values not finite). */
MYB_API MYB_NODISCARD myb_status myb_filter_start(myb_filter* filter, double t0,
                                                  const double y0[6], const double* covariance,
                                                  double t_end);

/* Propagates to the observation and applies it. innovation (may be NULL): observed - predicted
 * before the correction, `count` values. Errors: MYB_E_DOMAIN (an invalid observation, a time
 * behind the filter's or beyond t_end, an innovation covariance that is not positive
 * definite), the statuses of the propagation and the measurement. */
MYB_API MYB_NODISCARD myb_status myb_filter_update(myb_filter* filter,
                                                   const myb_observation* observation,
                                                   double* innovation);

/* Propagates to t without an observation (the covariance grows). Errors as above. */
MYB_API MYB_NODISCARD myb_status myb_filter_predict(myb_filter* filter, double t);

/* The time, the state and the covariance ((6 + np + bias_count)^2 values) of the filter; any
 * output may be NULL. */
MYB_API void myb_filter_get(const myb_filter* filter, double* t, double y[6],
                            double* covariance);

/* The systematic shifts: their first values (after myb_filter_start, before the first
 * observation; MYB_E_DOMAIN if not finite) and their present estimates. */
MYB_API MYB_NODISCARD myb_status myb_filter_set_bias(myb_filter* filter, const double* bias);
MYB_API void myb_filter_bias(const myb_filter* filter, double* bias);

/* ---- Smoother ----------------------------------------------------------------------------------
 * The filter's estimate at an instant uses the observations up to it; the smoother (Rauch,
 * Tung, Striebel) goes back over the instants the filter kept (myb_filter_options.history) and
 * gives at each the estimate from all the observations, with its covariance. At the first
 * instant that is the orbit a fit over the whole arc would give. */

/* The instants kept so far. */
MYB_API size_t myb_filter_steps(const myb_filter* filter);

/* The backward pass over the instants kept; more observations may follow, and the pass be
 * repeated. Errors: MYB_E_DOMAIN (a filter without history, a predicted covariance that is
 * not positive definite). MYB_E_OVERFLOW comes from the filter itself when an observation
 * would need more instants than it keeps. */
MYB_API MYB_NODISCARD myb_status myb_filter_smooth(myb_filter* filter);

/* The estimate at an instant kept: its time, the 6 + np + bias_count unknowns (state,
 * parameters, shifts) and their covariance — the filter's own until myb_filter_smooth has
 * run, the smoothed ones after. Any output may be NULL. */
MYB_API void myb_filter_step(const myb_filter* filter, size_t step, double* t, double* unknowns,
                             double* covariance);

#ifdef __cplusplus
}
#endif

#endif /* MYB_ESTIM_H_ */
