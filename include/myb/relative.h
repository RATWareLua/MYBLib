/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/relative.h — the motion of one object about another (a deputy about a chief) in the
 * frame that turns with the chief: R along its position, N along its angular momentum,
 * T = N x R (the Hill frame; RTN of myb/frames.h). A relative state is six numbers: the
 * position of the deputy from the chief (m) and its velocity as seen in the turning frame
 * (m/s), components R, T, N.
 *
 * To first order in the separation the relative state at one time is a linear map of that at
 * another. For a chief on a circle the map is that of Clohessy and Wiltshire (W. H. Clohessy,
 * R. S. Wiltshire, J. Aerospace Sci. 27, 653, 1960; G. W. Hill, 1878); on an ellipse it is the
 * transition matrix of K. Yamanaka and F. Ankersen (J. Guid. Control Dyn. 25, 60, 2002). Here
 * it is built for any conic from the transition matrix of the two-body motion of the chief
 * (myb_kepler_stm) taken into the turning frame at both ends — the same linear map, without
 * a formula per kind of orbit. */
#ifndef MYB_RELATIVE_H_
#define MYB_RELATIVE_H_

#include "myb/core.h"
#include "myb/kepler.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The relative state of the deputy from the inertial states of both, exactly, and back.
 * Errors: MYB_E_DOMAIN (a chief at the centre or on a line through it, non-finite input). */
MYB_API MYB_NODISCARD myb_status myb_relative_from_inertial(const myb_pv* chief,
                                                            const myb_pv* deputy,
                                                            double relative[6]);
MYB_API MYB_NODISCARD myb_status myb_relative_to_inertial(const myb_pv* chief,
                                                          const double relative[6],
                                                          myb_pv* deputy);

/* relative(t0 + dt) = phi relative(t0) to first order, for the chief `chief` at t0 moving as a
 * two-body orbit of the centre mu (any conic); phi is 6 x 6, row-major. The error of the map
 * is of the second order of the separation: (separation / distance to the centre) of the
 * separation, growing with dt. Errors: those of myb_kepler_stm; MYB_E_DOMAIN. */
MYB_API MYB_NODISCARD myb_status myb_relative_stm(double mu, const myb_pv* chief, double dt,
                                                  double phi[36]);

/* The same for a chief on a circle with the mean motion n (rad/s), in closed form
 * (Clohessy-Wiltshire); n must not be 0. */
MYB_API void myb_cw_stm(double n, double dt, double phi[36]);

/* Two impulses of a rendezvous in the linear motion: the first, at t0, takes the deputy from
 * the relative state `relative` to the relative position `target` at t0 + dt; the second,
 * there, leaves it with the relative velocity `target_velocity` (NULL: at rest with respect
 * to the chief). dv0 and dv1 in the turning frame of the chief at the time of each (m/s).
 * Errors: as myb_relative_stm; MYB_E_DOMAIN (a dt for which the position does not depend on
 * the velocity at t0 in every direction, as after a whole revolution). */
MYB_API MYB_NODISCARD myb_status myb_relative_rendezvous(double mu, const myb_pv* chief,
                                                         double dt, const double relative[6],
                                                         const double target[3],
                                                         const double target_velocity[3],
                                                         double dv0[3], double dv1[3]);

/* Relative orbital elements of the deputy (S. D'Amico, PhD thesis, TU Delft, 2010), which
 * do not change in two-body motion but for the drift of the second:
 *   roe = (da, dlambda, dex, dey, dix, diy),
 *   da = (a_d - a) / a,        dlambda = (u_d - u) + (raan_d - raan) cos i,
 *   dex = e_d cos w_d - e cos w,  dey = e_d sin w_d - e sin w,
 *   dix = i_d - i,             diy = (raan_d - raan) sin i,
 * u the mean argument of latitude. Both orbits must be ellipses; the elements are singular
 * for a chief in the equator of the frame. Errors: MYB_E_DOMAIN. */
MYB_API MYB_NODISCARD myb_status myb_roe_from_states(double mu, const myb_pv* chief,
                                                     const myb_pv* deputy, double roe[6]);

/* The relative position (R, T, N) that the elements give at the mean argument of latitude u
 * of a chief of semi-major axis a on a nearly circular orbit, to first order in the elements
 * and in the eccentricity:
 *   R = a (da - dex cos u - dey sin u),  T = a (dlambda + 2 dex sin u - 2 dey cos u),
 *   N = a (dix sin u - diy cos u). */
MYB_API void myb_roe_to_position(double a, double u, const double roe[6], double position[3]);

#ifdef __cplusplus
}
#endif

#endif /* MYB_RELATIVE_H_ */
