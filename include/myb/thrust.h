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

#ifdef __cplusplus
}
#endif

#endif /* MYB_THRUST_H_ */
