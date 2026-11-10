/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/observe.h — observables: light time, relativistic delay, range and range rate, round
 * trips, aberration, light deflection, right ascension and declination.
 *
 * Positions are in one inertial frame (the BCRS for the solar system) with the time coordinate
 * of that frame (TDB): the caller supplies the states of the participants as functions of time
 * (e.g. through myb_spk_state). The light-time equation of a one-way signal received at t by
 * an observer at x_o(t) from a target x_T,
 *   c (t - t_e) = |x_o(t) - x_T(t_e)| + c sum_B Dt_B,
 * is solved by fixed-point iteration (contraction |v_T| / c) to the resolution of the time;
 * Dt_B is the delay in the field of body B (I. I. Shapiro, Phys. Rev. Lett. 13, 789, 1964;
 * T. D. Moyer, Formulation for Observed and Computed Values of Deep Space Network Data Types,
 * DESCANSO Monograph 2, 2003, eq. 8-38, 8-55):
 *   Dt_B = (1 + gamma) gm_B / c^3 ln[(r1 + r2 + r12 + (1 + gamma) gm_B / c^2) /
 *                                     (r1 + r2 - r12 + (1 + gamma) gm_B / c^2)],
 * r1, r2 the distances of the end points from the body and r12 between them; the body is taken
 * at its position when the ray passes it, approximated as t - |x_B - x_o| / c (Moyer §8.3.1). */
#ifndef MYB_OBSERVE_H_
#define MYB_OBSERVE_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/kin.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The participants are providers (myb_state_function of myb/kin.h): the state of a point at
 * time t in the common frame. */

/* A body delaying (and deflecting) light. */
typedef struct myb_light_body {
  double gm;
  myb_state_function state;
  void* context;
} myb_light_body;

typedef struct myb_light_model {
  double c;                     /* m/s */
  double gamma;                 /* PPN gamma (1 in general relativity) */
  const myb_light_body* bodies; /* may be NULL: Newtonian light time */
  size_t count;
} myb_light_model;

typedef enum myb_light_direction {
  MYB_LIGHT_RECEPTION = 0,   /* the observer receives at t what the target emitted at t_e < t */
  MYB_LIGHT_TRANSMISSION = 1 /* the observer emits at t, the target receives at t_a > t */
} myb_light_direction;

typedef struct myb_light_path {
  myb_time t_target;  /* emission (reception) or arrival (transmission) time at the target */
  double target[6];   /* the target's state then */
  double light_time;  /* s, |t - t_target| */
  double range;       /* m, |x_T(t_target) - x_o(t)| */
  double delay;       /* s, the relativistic delay included in light_time */
  double range_rate;  /* m/s, d range / dt at the observer (Newtonian light time) */
  int iterations;
} myb_light_path;

/* One-way light time between the observer state `observer` at t and a target. Errors:
 * MYB_E_DOMAIN (invalid model), MYB_E_NOCONV (no convergence in 32 iterations), statuses of the
 * state functions. */
MYB_API MYB_NODISCARD myb_status myb_light_time(const myb_light_model* model,
                                                myb_state_function target, void* context,
                                                myb_time t, const double observer[6],
                                                myb_light_direction direction,
                                                myb_light_path* path);

/* The relativistic delay (s) of a ray from x1 to x2 by a body of gravitational parameter gm at
 * xb (Moyer eq. 8-38 above). */
MYB_API double myb_shapiro_delay(double gm, double c, double gamma, const double x1[3],
                                 const double x2[3], const double xb[3]);

/* Round trip received at t by a station: up-link transmitted at t_transmit, reflected or
 * transponded by the target at t_bounce. *total = t - t_transmit (s); the other outputs may be
 * NULL. Errors as myb_light_time. */
MYB_API MYB_NODISCARD myb_status myb_round_trip(const myb_light_model* model,
                                                myb_state_function station,
                                                void* station_context,
                                                myb_state_function target,
                                                void* target_context, myb_time t, double* total,
                                                myb_time* t_bounce, myb_time* t_transmit);

/* ---- Directions ---------------------------------------------------------------------------- */

/* Aberration: the natural direction p (unit, observer to source) as seen by an observer moving
 * with velocity v (m/s) in the frame, by the Lorentz transformation (exact in v / c):
 *   p' = (p / g + (1 + p.b / (1 + 1/g)) b) / |..|,  b = v / c, 1/g = sqrt(1 - b^2). */
MYB_API void myb_aberration(const double p[3], const double v[3], double c, double out[3]);

/* Light deflection by a body (first order in gm, PPN): the direction p (unit, observer to
 * source) becomes
 *   p + (1 + gamma) gm / (c^2 em) / max(q.(q + e), dlim) (p x (e x q)),
 * q the unit vector body -> source, e body -> observer, em the body-observer distance (m),
 * dlim > 0 the limiter near the limb (ERFA eraLd, of which this is the PPN form; not
 * normalised). */
MYB_API void myb_light_deflection(double gm, double c, double gamma, const double p[3],
                                  const double q[3], const double e[3], double em, double dlim,
                                  double out[3]);

/* Right ascension in [0, 2 pi) and declination of a direction (any length). */
MYB_API void myb_direction_radec(const double p[3], double* ra, double* dec);

/* ---- Dispersive delay in a plasma ----------------------------------------------------------
 * A radio signal of frequency f (Hz) crossing a column of N free electrons per m^2 (the total
 * electron content along the path) is delayed, to first order in (f_plasma / f)^2,
 *   group delay = +K N / f^2,  phase delay = -K N / f^2 (metres; the phase advances),
 *   K = e^2 / (8 pi^2 eps0 m_e) = 40.308 m^3/s^2
 * (e, m_e, eps0: MYB_CONST_E, MYB_CONST_ME, MYB_CONST_EPS0 of data/constants/codata-2022.tk).
 * Ranging measures the group delay, carrier phase and Doppler the phase delay; two frequencies
 * remove it. Higher orders ((f_p / f)^4, the magnetic field) are below 1 mm above 1 GHz on
 * Earth paths. Returns K N / f^2, NaN for invalid input. */
MYB_API double myb_plasma_group_delay(double column, double frequency, double e, double m_e,
                                      double eps0);

/* Power-law electron density about the centre of a body, n_e(r) = sum_k density[k] (radius /
 * r)^power[k] (m^-3), as the models of the solar corona and wind (e.g. Muhleman & Anderson,
 * Astrophys. J. 247, 1093, 1981: ~ r^-6 near the Sun and ~ r^-2 in the wind). */
enum { MYB_PLASMA_TERMS = 8 };
typedef struct myb_plasma_profile {
  double radius;                    /* reference radius (m, > 0) */
  size_t count;                     /* terms used, 1..MYB_PLASMA_TERMS */
  double density[MYB_PLASMA_TERMS]; /* n_k (m^-3, >= 0) */
  double power[MYB_PLASMA_TERMS];   /* p_k (> 0) */
} myb_plasma_profile;

/* Electron column (m^-2) of the profile along the straight segment from p1 to p2 (positions
 * relative to the centre of the body, m): adaptive Gauss-Kronrod quadrature in the angle seen
 * from the centre, relative accuracy 1e-13. Errors: MYB_E_DOMAIN (invalid profile, a segment
 * through the centre, non-finite input), MYB_E_NOCONV (quadrature budget exhausted). */
MYB_API MYB_NODISCARD myb_status myb_plasma_column(const myb_plasma_profile* profile,
                                                   const double p1[3], const double p2[3],
                                                   double* column);

#ifdef __cplusplus
}
#endif

#endif /* MYB_OBSERVE_H_ */
