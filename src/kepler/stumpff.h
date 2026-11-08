/* SPDX-License-Identifier: BSD-3-Clause */
/* kepler/stumpff.h — Stumpff functions and the angle reduction used by the two-body code.
 *
 *   c2(z) = (1 - cos sqrt(z)) / z        c3(z) = (sqrt(z) - sin sqrt(z)) / z^(3/2)
 * continued analytically to z < 0 (cosh, sinh) and to z = 0 (1/2, 1/6). The forms used avoid
 * every cancellation: c2 = 2 (sin(s/2) / s)^2, and c3 by its power series for |z| < 1, so
 * that x - sin x = x^3 c3(x^2) and sinh x - x = x^3 c3(-x^2) keep full relative accuracy
 * near zero (near-parabolic orbits). Source: K. Stumpff, Himmelsmechanik (1959); R. H. Battin,
 * An Introduction to the Mathematics and Methods of Astrodynamics (1999), ch. 4. */
#ifndef MYB_SRC_KEPLER_STUMPFF_H_
#define MYB_SRC_KEPLER_STUMPFF_H_

double myb_internal_stumpff_c2(double z);
double myb_internal_stumpff_c3(double z);

/* c4(z) = (1/2 - c2(z)) / z and c5(z) = (1/6 - c3(z)) / z (1/24 and 1/120 at z = 0), for the
 * derivatives of the universal functions with respect to alpha. */
double myb_internal_stumpff_c4(double z);
double myb_internal_stumpff_c5(double z);

/* cos(s) for z = s^2 > 0, cosh(s) for z = -s^2 < 0, 1 for z = 0. */
double myb_internal_stumpff_c0(double z);

/* sin(s)/s for z = s^2 > 0, sinh(s)/s for z = -s^2 < 0, 1 for z = 0. */
double myb_internal_stumpff_c1(double z);

/* x reduced modulo 2*pi to [-pi, pi] with a double-word 2*pi (accurate for |x| < 2^40);
 * *turns receives the number of whole turns removed. */
double myb_internal_reduce_angle(double x, double* turns);

/* turns * 2*pi + x, with the double-word 2*pi. */
double myb_internal_add_turns(double x, double turns);

#endif /* MYB_SRC_KEPLER_STUMPFF_H_ */
