/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/conjunction.h — how close two objects come, and how likely they are to touch.
 *
 * The closest approaches of two objects are the minima of their distance: the instants at
 * which the range rate turns from negative to positive. They are looked for between the ends
 * of equal steps (two approaches closer together than a step may be missed).
 *
 * The probability of collision is that of a short encounter (J. L. Foster, H. S. Estes,
 * NASA/JSC-25898, 1992; S. Alfano, J. Astronaut. Sci. 53, 103, 2005): the two pass each other
 * in a straight line, quickly against the time in which their uncertainty changes, so that
 * only the plane across the relative velocity matters; the position of one from the other is
 * Gaussian there, with the sum of the two covariances; they collide if the miss is less than
 * the sum of their radii. The probability is the integral of that Gaussian over a disk. */
#ifndef MYB_CONJUNCTION_H_
#define MYB_CONJUNCTION_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/kin.h"
#include "myb/orbit.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct myb_approach {
  double t;           /* the time of the closest approach (s) */
  double distance;    /* of the two there (m) */
  double position[3]; /* of the second from the first (m) */
  double velocity[3]; /* of the second relative to the first (m/s) */
} myb_approach;

/* The closest approaches of two objects given by providers (myb_state_function of myb/kin.h;
 * both in one frame, at the instant myb_time_from_seconds(t)) for t in [t0, t1], looked for in
 * steps of `step` and located to the resolution of the time. Up to `capacity` are written, in
 * the order of the time; *count receives how many were found (it may exceed the capacity).
 * Errors: MYB_E_DOMAIN (t1 <= t0, step <= 0, non-finite input), the statuses of the providers. */
MYB_API MYB_NODISCARD myb_status myb_approach_find(myb_state_function first, void* first_context,
                                                   myb_state_function second,
                                                   void* second_context, double t0, double t1,
                                                   double step, myb_approach* approaches,
                                                   size_t capacity, size_t* count);

/* The same for two orbits being propagated (myb/orbit.h), both started and standing at one
 * time: they are taken together to t1 in steps of `step`, and an approach inside a step is
 * located on the cubic that joins the relative states at its ends — within step^4 / 384 of
 * the fourth derivative of the relative position (a quarter of a millimetre for 10 s in a
 * low orbit, a third of a metre for 60 s). Both orbits are left at t1. Errors: MYB_E_DOMAIN
 * (the orbits at different times, t1 not ahead of them, step <= 0), those of myb_orbit_state. */
MYB_API MYB_NODISCARD myb_status myb_approach_orbits(myb_orbit* first, myb_orbit* second,
                                                     double t1, double step,
                                                     myb_approach* approaches, size_t capacity,
                                                     size_t* count);

/* ---- The probability of collision ------------------------------------------------------------ */

/* An encounter seen in the plane across the relative velocity, in the principal axes of the
 * covariance projected onto it. */
typedef struct myb_encounter {
  double miss[2];  /* the miss vector (m) */
  double sigma[2]; /* standard deviations along the axes (m), sigma[0] >= sigma[1] > 0 */
  double radius;   /* the sum of the radii of the two (m), > 0 */
} myb_encounter;

/* From the relative position and velocity at the closest approach and the covariance of the
 * relative position (3 x 3, row-major, m^2: the sum of the covariances of the two positions,
 * in the frame of the vectors). Errors: MYB_E_DOMAIN (zero relative velocity, a covariance
 * that is not symmetric or not positive in the plane, radius <= 0, non-finite input). */
MYB_API MYB_NODISCARD myb_status myb_encounter_plane(const double position[3],
                                                     const double velocity[3],
                                                     const double covariance[9], double radius,
                                                     myb_encounter* encounter);

/* The integral of the Gaussian over the disk: one integral across the disk of the difference
 * of two error functions, by Gauss-Legendre panels in the angle that removes the roots at the
 * rim, doubled until the value settles (relative 1e-13, or absolute 1e-300). Errors:
 * MYB_E_DOMAIN (an invalid encounter), MYB_E_NOCONV. */
MYB_API MYB_NODISCARD myb_status myb_collision_probability(const myb_encounter* encounter,
                                                           double* probability);

/* The series of F. K. Chan (Spacecraft Collision Probability, 2008, ch. 4): the ellipse of
 * the covariance made a circle of the same area, which turns the integral into the
 * distribution of Rice, summed (all terms positive) until they fall below 1e-17 of the sum.
 * Exact for equal sigmas; for unequal ones an approximation, good while the disk is small
 * against both. NaN for an invalid encounter, and where R^2 / (2 sigma_x sigma_y) or half the
 * squared miss in sigmas reaches 700 (the exponentials leave the doubles). */
MYB_API double myb_collision_probability_chan(const myb_encounter* encounter);

/* The largest probability that the encounter can have if the covariance is known in shape
 * only — every sigma multiplied by one unknown factor — and the factor at which it is
 * reached (S. Alfano, AAS 03-548). Errors: as myb_collision_probability. */
MYB_API MYB_NODISCARD myb_status myb_collision_probability_max(const myb_encounter* encounter,
                                                               double* probability,
                                                               double* factor);

/* ---- Approaches with the uncertainty carried to them ----------------------------------------- */

typedef struct myb_conjunction {
  myb_approach approach;
  double covariance[9];    /* of the position of the second from the first there (m^2) */
  myb_encounter encounter; /* in the plane across the relative velocity */
  double probability;      /* of collision (myb_collision_probability) */
} myb_conjunction;

/* myb_approach_orbits with the covariances of the two states (6 x 6, row-major, m and m/s)
 * at the time the orbits were started, carried to every approach by the transition matrices
 * of the orbits themselves — both created with `variational` (myb/orbit.h), so that whatever
 * the force models hold acts on the uncertainty as it does on the motion. The covariance of
 * the relative position is Phi_r P Phi_r^T of the first plus that of the second, Phi_r the
 * rows of the positions, taken inside a step on the same cubic as the states (its rows of
 * the velocities are the rates of those of the positions); then the plane of the encounter
 * and the probability for the sum of the radii `radius`. Errors: as myb_approach_orbits,
 * myb_orbit_partials, myb_encounter_plane and myb_collision_probability. */
MYB_API MYB_NODISCARD myb_status myb_conjunction_find(myb_orbit* first,
                                                      const double first_covariance[36],
                                                      myb_orbit* second,
                                                      const double second_covariance[36],
                                                      double radius, double t1, double step,
                                                      myb_conjunction* conjunctions,
                                                      size_t capacity, size_t* count);

#ifdef __cplusplus
}
#endif

#endif /* MYB_CONJUNCTION_H_ */
