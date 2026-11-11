/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/iod.h — a first orbit from a few measurements (initial orbit determination), for the
 * two-body motion about a centre of gravitational parameter mu. Positions in metres,
 * velocities in m/s, times in seconds, all vectors in one inertial frame centred on the body.
 * What comes out is a state to start an estimator with (myb/estim.h), not an estimate: the
 * methods take the measurements as exact.
 *
 * Two positions and the time between them are Lambert's problem (myb/lambert.h). */
#ifndef MYB_IOD_H_
#define MYB_IOD_H_

#include <stddef.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Three positions of one orbit, in the order of the motion: the velocity at the second
 * (J. W. Gibbs, 1889; in the form of D. A. Vallado, Fundamentals of Astrodynamics and
 * Applications, 4th ed., algorithm 54). Exact for any three points of a conic: the times are
 * not used. *out_of_plane (may be NULL) receives the angle (rad) between r1 and the plane of
 * r2 and r3, which measures how far the three are from one orbit. Close positions (a few
 * degrees apart) lose digits in the cross products: use myb_iod_herrick_gibbs there. Errors:
 * MYB_E_DOMAIN (a zero or non-finite vector, mu <= 0, positions on one line). */
MYB_API MYB_NODISCARD myb_status myb_iod_gibbs(double mu, const double r1[3], const double r2[3],
                                               const double r3[3], double v2[3],
                                               double* out_of_plane);

/* Three positions close together in time, t1 < t2 < t3: the velocity at the second from the
 * Taylor series of the motion to the fourth order (S. Herrick; Vallado, algorithm 55). The
 * error falls as the fourth power of the spacing: three parts in a million for a tenth of a
 * radian of the orbit between neighbours, seven in 10^8 for a twenty-fifth. Errors:
 * MYB_E_DOMAIN (times not increasing, a zero or non-finite vector, mu <= 0). */
MYB_API MYB_NODISCARD myb_status myb_iod_herrick_gibbs(double mu, const double r1[3], double t1,
                                                       const double r2[3], double t2,
                                                       const double r3[3], double t3,
                                                       double v2[3]);

enum { MYB_IOD_MAX_SOLUTIONS = 3 };

/* An orbit from the directions alone. */
typedef struct myb_iod_solution {
  double r[3];      /* position at the second time */
  double v[3];      /* velocity there */
  double range[3];  /* distances from the three observers (m), all > 0 */
  int iterations;   /* Newton steps of the refinement */
} myb_iod_solution;

/* Three directions to the object (unit vectors `sight`) from three places (`observer`: the
 * positions of the observers relative to the centre) at the times t[0] < t[1] < t[2] — Gauss's
 * method (Vallado, algorithm 52; H. D. Curtis, Orbital Mechanics for Engineering Students,
 * algorithms 5.5 and 5.6). The distance at the middle time is a root of
 * x^8 + a x^6 + b x^3 + c = 0, which may have up to three positive roots: each that gives
 * positive ranges is a candidate orbit. Every candidate is then refined by Newton's method
 * on the six conditions that the three sights set (the object lies on each of them) for the
 * state at the middle time, with the transition matrix of the two-body motion, until the
 * step is at the rounding of the state (at most 50 steps): the result passes through all
 * three sights in exact two-body motion, over arcs that the series of Gauss alone cannot
 * span. Candidates that do not settle, or settle behind an observer, are dropped; two roots
 * that settle on one orbit are given once. *count receives how many
 * came out, at most MYB_IOD_MAX_SOLUTIONS, in the order of the roots they started from. Errors:
 * MYB_E_DOMAIN (times not increasing, a sight that is not a unit vector to 1e-9, coplanar
 * sights, mu <= 0, non-finite input), MYB_E_NOCONV (no candidate settled). */
MYB_API MYB_NODISCARD myb_status myb_iod_gauss(double mu, const double t[3],
                                               const double observer[3][3],
                                               const double sight[3][3],
                                               myb_iod_solution solutions[MYB_IOD_MAX_SOLUTIONS],
                                               size_t* count);

#ifdef __cplusplus
}
#endif

#endif /* MYB_IOD_H_ */
