/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/meteoroid.h — random impacts of micrometeoroids: a Poisson process of impulses.
 *
 * The flux is the interplanetary flux model of E. Grun, H. A. Zook, H. Fechtig, R. H. Giese
 * (Icarus 62, 244, 1985, eq. A3): the cumulative flux of the particles of mass m or more at
 * 1 au, in m^-2 s^-1 with m in grams,
 *   F(m) = (c4 m^g4 + c5)^g5 + c6 (m + c7 m^g6 + c8 m^g7)^g8 + c9 (m + c10 m^g9)^g10,
 *   c4 = 2.2e3, c5 = 15, c6 = 1.3e-9, c7 = 1e11, c8 = 1e27, c9 = 1.3e-16, c10 = 1e6,
 *   g4 = 0.306, g5 = -4.38, g6 = 2, g7 = 4, g8 = -0.36, g9 = 2, g10 = -0.85,
 * for masses from about 1e-18 g to 100 g (the library takes masses in kg and converts): per
 * m^2 and year, 1088 particles above 1e-15 kg, 95 above 1e-12 kg, 1.49 above 1e-9 kg, 6e-4
 * above 1e-6 kg. F is the flux on one side of a flat plate turned every way, a quarter of the
 * flux from all directions: an object of mean cross-section A (its projected area averaged
 * over directions: pi R^2 for a sphere, a quarter of the surface of any convex body) takes
 * 4 A F impacts per second. Away from 1 au the flux is scaled by (r / 1 au)^-p, the exponent
 * the caller's: the density of the zodiacal cloud falls as r^-1.3 (C. Leinert, I. Richter,
 * E. Pitz, B. Planck, Astron. Astrophys. 103, 177, 1981), and what an object meets depends
 * also on its own speed through the cloud.
 *
 * The process. For particles between m_min and m_max the impacts come at the rate
 *   lambda = 4 A (F(m_min) - F(m_max)) (r / 1 au)^-p,
 * independent of each other: the waits between them are exponential, -ln(1 - u) / lambda; the
 * mass of each is drawn from the distribution F gives, the mass at which F(m) = F(m_min) - u
 * (F(m_min) - F(m_max)) (Newton's iteration in ln m); the speed is the caller's — one value for
 * all, 20 km/s by default (the mean impact speed Grun et al. take at 1 au), or drawn by a
 * function; the direction of the particle is isotropic, (sqrt(1 - z^2) cos phi, sqrt(1 - z^2)
 * sin phi, z) with z = 2 u - 1, phi = 2 pi u', or drawn by a function (an apex, a stream). An
 * impact gives the object the momentum beta m v along the direction the particle moves — beta
 * >= 1 counts the ejecta thrown back out of the crater —, and so the change of velocity beta m
 * v / M, for myb_orbit_kick (myb/orbit.h).
 *
 * Reproducibility: everything is drawn from the generator of myb/uncert.h (xoshiro256**) in a
 * fixed order — for each impact the wait, the mass, then the speed and the direction (a
 * function of the caller draws from the same generator) — and computed with the correctly
 * rounded functions of the library: one seed gives the same impacts, bit for bit, on every
 * platform.
 *
 * Measured (tests/unit/test_meteoroid.c): F against the formula evaluated with the C library's
 * pow to 6e-15 from 1e-21 to 0.1 kg. On 1 m^2 between 1e-15 and 1e-9 kg (1.38e-4 impacts a
 * second) 50085 impacts against 50000 expected, 0.4 sigma; 2000 runs of a mean of 20: the
 * mean 20.06 (0.6 sigma) and the variance 19.26 (-1.2 sigma) of Poisson's statistics. The
 * masses: the inverse of the distribution to 3e-15 of F, and 50000 drawn give the statistic of
 * Kolmogorov and Smirnov sqrt(n) D = 1.17 (1.95 is the bound at 0.001). The total impulse
 * (beta 2, 20 km/s, 1e-12 to 1e-11 kg, 1000 runs of a mean of 100 impacts): the sum of |J|
 * against lambda T beta v E[m] 0.3 sigma, a component 0.9 sigma from zero and its variance
 * against lambda T (beta v)^2 E[m^2] / 3 at -0.5 sigma (E[m], E[m^2] by quadrature of F). 310
 * impacts applied with myb_orbit_kick to a free motion: the position and the velocity of the
 * free motion with the kicks to 7e-16. */
#ifndef MYB_METEOROID_H_
#define MYB_METEOROID_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/orbit.h"
#include "myb/uncert.h"

#ifdef __cplusplus
extern "C" {
#endif

/* F(m) of Grun et al. at 1 au: the impacts per m^2 per s of the particles of mass >= mass (kg)
 * on one side of a plate turned every way. NaN for a mass that is not positive and finite. */
MYB_API double myb_meteoroid_flux(double mass);

/* A speed (m/s, > 0) or a direction (any length, scaled to unit; the way the particle moves,
 * in the inertial axes of the model) drawn for an impact at the model time t, from the
 * generator given. Any status other than MYB_OK stops the drawing and is returned. */
typedef myb_status (*myb_meteoroid_speed_function)(void* context, myb_rng* rng, double t,
                                                   double* speed);
typedef myb_status (*myb_meteoroid_direction_function)(void* context, myb_rng* rng, double t,
                                                       double direction[3]);

typedef struct myb_meteoroid_options {
  double cross_section; /* m^2: the mean cross-section of the object, > 0 */
  double mass_min;      /* kg, > 0 */
  double mass_max;      /* kg, > mass_min */
  double distance;      /* m: from the Sun, for the scaling; 0: 1 au */
  double exponent;      /* p of (distance / 1 au)^-p */
  double speed;         /* m/s, > 0; 0: 20 km/s. Not used with a speed function */
  myb_meteoroid_speed_function speed_function; /* NULL: speed */
  void* speed_context;
  myb_meteoroid_direction_function direction_function; /* NULL: isotropic */
  void* direction_context;
  double enhancement; /* beta, > 0; 0: 1 (the particle stays in the object, no ejecta) */
  double object_mass; /* M (kg), > 0 */
} myb_meteoroid_options;

typedef struct myb_meteoroid_impact {
  double t;            /* model time (s) */
  double mass;         /* kg */
  double speed;        /* m/s */
  double direction[3]; /* unit: the way the particle moves */
  double impulse[3];   /* N s: beta m v direction */
  double dv[3];        /* m/s: impulse / M */
} myb_meteoroid_impact;

/* lambda: the impacts per second. Errors: MYB_E_DOMAIN (invalid options). */
MYB_API MYB_NODISCARD myb_status myb_meteoroid_rate(const myb_meteoroid_options* options,
                                                    double* rate);

/* The mass (kg) below which a fraction u in [0, 1) of the impacts fall: the inverse of the
 * distribution of the masses. Errors: MYB_E_DOMAIN (invalid options, u outside [0, 1)),
 * MYB_E_NOCONV. */
MYB_API MYB_NODISCARD myb_status myb_meteoroid_mass(const myb_meteoroid_options* options, double u,
                                                    double* mass);

/* The next impact after the model time t. Errors: MYB_E_DOMAIN (invalid options, t not
 * finite, a speed or a direction drawn that is not valid), MYB_E_NOCONV, the statuses of the
 * functions of the caller. */
MYB_API MYB_NODISCARD myb_status myb_meteoroid_next(const myb_meteoroid_options* options,
                                                    myb_rng* rng, double t,
                                                    myb_meteoroid_impact* impact);

/* Advances the orbit from the time reached to t_end (forward, within its own t_end) with the
 * impacts on the way: at each the state is asked (myb_orbit_state) and the velocity changed
 * (myb_orbit_kick); the time reached is t_end at the end. The impacts are written to
 * impacts[0 .. capacity - 1] (may be NULL with capacity 0) and all of them counted in *count.
 * A later call draws from its own start: the same statistics (the process keeps no memory),
 * not the same impacts as one call over both intervals. Errors: MYB_E_DOMAIN (t_end behind the
 * time reached or not finite, invalid options), MYB_E_NOCONV (far more impacts than lambda T
 * allows: 1000 + lambda T + 20 sqrt(lambda T)), those of myb_meteoroid_next and of the orbit. */
MYB_API MYB_NODISCARD myb_status myb_meteoroid_propagate(myb_orbit* orbit,
                                                         const myb_meteoroid_options* options,
                                                         myb_rng* rng, double t_end,
                                                         myb_meteoroid_impact* impacts,
                                                         size_t capacity, size_t* count);

#ifdef __cplusplus
}
#endif

#endif /* MYB_METEOROID_H_ */
