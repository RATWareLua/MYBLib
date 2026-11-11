/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/thrust.h — continuous (low) thrust: burn arcs with a steering law as a term of a force
 * model, the Gauss equations in modified equinoctial elements, Edelbaum's transfer.
 *
 * A burn arc [start, stop] has a constant thrust T (N) and effective exhaust velocity c (m/s;
 * c = Isp g0 with the conventional g0 = 9.80665 m/s^2), so the mass falls linearly at T / c and
 * the acceleration T / m(t) grows: the rocket equation dv = c ln(m0 / m1) holds exactly for the
 * acceleration integrated along the arc. The mass is a function of time, so the state of an
 * integrator stays (x, v); arcs may not overlap. */
#ifndef MYB_THRUST_H_
#define MYB_THRUST_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/kepler.h"
#include "myb/model.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum myb_steering {
  MYB_STEER_VELOCITY = 0,      /* along v (prograde) */
  MYB_STEER_ANTI_VELOCITY = 1, /* against v */
  MYB_STEER_INERTIAL = 2,      /* fixed direction in the frame of x */
  MYB_STEER_RTN = 3,           /* fixed components in the radial, transverse, normal frame */
  MYB_STEER_EDELBAUM = 4       /* Edelbaum's yaw law (below) */
} myb_steering;

typedef struct myb_thrust_arc {
  double start;           /* s, on the time of the force model */
  double stop;            /* > start */
  double thrust;          /* N, > 0 */
  double exhaust;         /* effective exhaust velocity, m/s, > 0 */
  myb_steering steering;
  double direction[3];    /* INERTIAL and RTN: the direction (normalised at init) */
  /* EDELBAUM (J. Edelbaum, ARS J. 31, 1079, 1961; the yaw law as in J. A. Kechichian, J. Guid.
   * Control Dyn. 20, 803, 1997): thrust at the angle beta out of the orbital plane, along the
   * velocity in the plane, the out-of-plane sign following cos(u) (u the argument of latitude)
   * so that the inclination changes in `sense` (+1 increase, -1 decrease):
   *   tan beta = v0 sin beta0 / (v0 cos beta0 - dv),  dv = c ln(m(start) / m(t)). */
  double edelbaum_v0;     /* initial circular speed, m/s */
  double edelbaum_beta0;  /* initial yaw, rad (myb_edelbaum gives it) */
  double edelbaum_sense;  /* +1 or -1 */
} myb_thrust_arc;

enum { MYB_THRUST_MAX_ARCS = 16 };

typedef struct myb_thrust {
  double mass0; /* kg at the start of the first arc */
  size_t count;
  myb_thrust_arc arcs[MYB_THRUST_MAX_ARCS];
} myb_thrust;

/* Errors: MYB_E_DOMAIN (mass <= 0, count > MYB_THRUST_MAX_ARCS, an invalid arc, overlapping or
 * unordered arcs, a zero direction, the propellant exceeding the mass). */
MYB_API MYB_NODISCARD myb_status myb_thrust_init(myb_thrust* thrust, double mass0,
                                                 const myb_thrust_arc* arcs, size_t count);

/* Mass at t (kg). */
MYB_API double myb_thrust_mass(const myb_thrust* thrust, double t);

/* The same as a mass function of the object (myb_object.mass_function, myb/env.h; context:
 * const myb_thrust*): the surface forces then act on the mass that is left. */
MYB_API double myb_thrust_mass_at(void* thrust, double t);

/* The acceleration of the thrust at t (zero between arcs), x and v relative to the central body:
 * a term of a force model (myb_term_function, context: const myb_thrust*). */
MYB_API MYB_NODISCARD myb_status myb_thrust_acceleration(void* thrust, double t, const double x[3],
                                                         const double v[3], double a[3]);

/* The same acceleration in the radial, transverse, normal frame of (x, v). */
MYB_API MYB_NODISCARD myb_status myb_thrust_rtn(const myb_thrust* thrust, double t,
                                                const double x[3], const double v[3],
                                                double a_rtn[3]);

/* ---- Gauss equations in modified equinoctial elements -----------------------------------------
 * (M. J. H. Walker, B. Ireland, J. Owens, Celest. Mech. 36, 409, 1985): with w = 1 + f cos L +
 * g sin L, s^2 = 1 + h^2 + k^2, q = sqrt(p / mu) and the perturbing acceleration (R, T, N),
 *   p' = 2 p q T / w
 *   f' = q [ R sin L + ((w + 1) cos L + f) T / w - (h sin L - k cos L) g N / w ]
 *   g' = q [-R cos L + ((w + 1) sin L + g) T / w + (h sin L - k cos L) f N / w ]
 *   h' = q s^2 cos L N / (2 w),   k' = q s^2 sin L N / (2 w)
 *   L' = sqrt(mu p) (w / p)^2 + q (h sin L - k cos L) N / w. */
MYB_API void myb_equinoctial_rates(double mu, const myb_equinoctial* elements,
                                   const double a_rtn[3], double rates[6]);

/* ---- Edelbaum's circle-to-circle transfer -----------------------------------------------------
 * Constant acceleration between circular orbits of speeds v0, v1 with the plane change di
 * (rad) in the orbit-averaged approximation (relative accuracy of the order of the ratio of the
 * thrust to the gravitational acceleration):
 *   dV = sqrt(v0^2 - 2 v0 v1 cos(pi di / 2) + v1^2),
 *   tan beta0 = sin(pi di / 2) / (v0 / v1 - cos(pi di / 2)).
 * Errors: MYB_E_DOMAIN (speeds not positive). */
MYB_API MYB_NODISCARD myb_status myb_edelbaum(double v0, double v1, double di, double* dv,
                                              double* beta0);

/* Speed v and inclination change di after the velocity increment dv along Edelbaum's law:
 *   v = sqrt(v0^2 - 2 v0 dv cos beta0 + dv^2),
 *   di = (2 / pi) [atan((dv - v0 cos beta0) / (v0 sin beta0)) + pi / 2 - beta0]. */
MYB_API void myb_edelbaum_state(double v0, double beta0, double dv, double* v, double* di);

/* ---- The Q-law ---------------------------------------------------------------------------------
 * A feedback law for a transfer between any two orbits (A. E. Petropoulos, AIAA 2004-5089):
 * Q measures how far the orbit is from the target,
 *   Q = (1 + W_P P) sum_k W_k S_k [(oe_k - target_k) / rate_k]^2,
 * the distance of each element over the largest rate at which thrust of unit acceleration can
 * change it anywhere on the orbit, so that Q is about the square of the time still needed;
 * the thrust points where Q falls fastest — the gradient of Q carried through the Gauss
 * equations at the true anomaly of the moment. P = exp(100 (1 - r_p / rp_min)) keeps the
 * periapsis up. No optimisation and no guess is needed, and what comes out is near the
 * optimum in time for the simple transfers (against Edelbaum's: the tests); it is the usual
 * first guess for an optimiser. The elements are a, e, i, raan, argp: the law is singular
 * where an angle that is targeted is not defined (argp on a circle, raan in the equator).
 *
 * The direction is the gradient made a unit vector, so it jumps where the gradient passes
 * through zero — places of the orbit where thrust gains nothing — and near the target, where
 * the terms pull against each other, an engine that always burns dithers there and stops
 * advancing. Petropoulos' remedy is the effectivity: the rate of Q at the place of the moment
 * over the best rate anywhere on the orbit; the engine burns only above a threshold (0.2 to
 * 0.5 is usual) and coasts through the rest, which also saves propellant. */
typedef struct myb_qlaw {
  double mu;        /* of the centre, m^3 s^-2 */
  double target[5]; /* a (m), e, i, raan, argp (rad) */
  double weight[5]; /* W_k >= 0; 0 leaves the element free; at least one above 0 */
  double rp_min;    /* the lowest periapsis radius wanted (m); 0: no limit */
  double rp_weight; /* W_P; 0 = 1 */
} myb_qlaw;

/* The unit direction of the thrust (radial, transverse, normal components) for the orbit and
 * the place on it given by `now`; *q (may be NULL) the value of Q for unit acceleration
 * (s^2: of the order of the square of the time still needed), *q_rate (may be NULL) its
 * rate per unit acceleration (s; never positive), *effectivity (may be NULL) that rate over
 * the best at any of 72 places of the orbit, in [0, 1]. A zero direction where Q does not
 * change with any thrust (at the target). Errors: MYB_E_DOMAIN (an invalid law, an orbit
 * that is not an ellipse, non-finite input). */
MYB_API MYB_NODISCARD myb_status myb_qlaw_direction(const myb_qlaw* law, const myb_elements* now,
                                                    double direction[3], double* q,
                                                    double* q_rate, double* effectivity);

/* 1 if every targeted element is within its tolerance (in the units of the target). */
MYB_API int myb_qlaw_arrived(const myb_qlaw* law, const myb_elements* now,
                             const double tolerance[5]);

/* The law as a term of a force model (myb_term_function; the context is this structure): an
 * engine of constant thrust that points as the law says, is off wherever the effectivity is
 * below the threshold and whenever the orbit is within the tolerances. The mass falls at
 * thrust / exhaust from t0 on, as if the engine had burned without a pause: exact up to the
 * arrival without a threshold, an overestimate of what is spent with one (a term has no
 * memory of when it burned: myb_qlaw_fly below carries the mass in the state).
 * The thrust goes on and off sharply, so an integrator with steps of its own choosing is the
 * one to use (MYB_ORBIT_DOP853 of myb/orbit.h). Errors: those of myb_qlaw_direction;
 * MYB_E_DOMAIN (the mass exhausted, a threshold outside [0, 1)). */
typedef struct myb_qlaw_term {
  const myb_qlaw* law;
  double thrust;       /* N */
  double exhaust;      /* effective exhaust velocity, m/s */
  double mass0;        /* kg at t0 */
  double t0;           /* s */
  double tolerance[5]; /* of the elements, as the target */
  double effectivity;  /* the threshold, in [0, 1); 0: the engine always burns */
} myb_qlaw_term;

MYB_API MYB_NODISCARD myb_status myb_term_qlaw(void* term, double t, const double x[3],
                                               const double v[3], double a[3]);

/* The same engine with the mass as a seventh component of the state, y = (x, v, m): the mass
 * falls at thrust / exhaust only while the engine burns, so that what the pauses save is
 * counted. mass0 and t0 of the term are not used. */
typedef struct myb_qlaw_flight {
  const myb_qlaw_term* engine;
  myb_force_model* field; /* everything but the engine, the centre included (myb/model.h);
                             NULL: the point mass of the law */
  double rtol;            /* of DOP853 in myb_qlaw_fly; 0 = 1e-10 */
  double atol;            /* 0 = 1e-4 (m, m/s, kg) */
} myb_qlaw_flight;

/* The right-hand side for an integrator of first-order systems (myb_ode_function of
 * myb/integ.h; the context is a myb_qlaw_flight; 7 components). Errors: MYB_E_DOMAIN (an
 * invalid engine, the mass exhausted, an orbit that is not an ellipse), those of the field. */
MYB_API MYB_NODISCARD myb_status myb_qlaw_flight_rhs(void* flight, double t, const double* y,
                                                     double* dydt);

typedef struct myb_qlaw_arrival {
  double t;    /* of the arrival, or t_limit (s) */
  double y[7]; /* x, v (m, m/s) and the mass (kg) there */
  int arrived; /* 1: every targeted element is within the tolerance of the engine's term */
} myb_qlaw_arrival;

/* Flies the law by DOP853 from y0 at t0 until the orbit is within the tolerances, looked at
 * every `step` seconds, or until t_limit. Errors: MYB_E_DOMAIN (t_limit <= t0, step <= 0, a
 * mass that is not positive, non-finite input), MYB_E_NOMEM, those of the right-hand side
 * and of the integrator; the state reached by then is in *arrival. allocator may be NULL. */
MYB_API MYB_NODISCARD myb_status myb_qlaw_fly(const myb_qlaw_flight* flight, double t0,
                                              const double y0[7], double t_limit, double step,
                                              const myb_allocator* allocator,
                                              myb_qlaw_arrival* arrival);

#ifdef __cplusplus
}
#endif

#endif /* MYB_THRUST_H_ */
