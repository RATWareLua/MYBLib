/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/cr3bp.h — the circular restricted three-body problem.
 *
 * Units: the distance of the primaries, their total mass and their mean motion are 1; the
 * rotating frame has its origin at the barycentre, the larger primary at (-mu, 0, 0), the
 * smaller at (1 - mu, 0, 0), mu = m2 / (m1 + m2). Equations of motion
 *   x'' - 2 y' = dOmega/dx,  y'' + 2 x' = dOmega/dy,  z'' = dOmega/dz,
 *   Omega = (x^2 + y^2) / 2 + (1 - mu) / r1 + mu / r2,
 * Jacobi constant C = 2 Omega - |v|^2 (V. Szebehely, Theory of Orbits, 1967).
 *
 * Periodic orbits symmetric with respect to the xz-plane (halo, Lyapunov, distant retrograde)
 * or to the x-axis (vertical, axial) are corrected by single shooting on the half period
 * (K. C. Howell, Celest. Mech. 32, 53, 1984): the state transition matrix comes from the
 * variational equations (myb/uncert.h) with the analytic Jacobian. Third-order halo orbits of
 * D. L. Richardson (Celest. Mech. 22, 241, 1980) give first guesses. Stability from the
 * eigenvalues of the monodromy matrix; invariant manifolds from its real eigenvectors (W. S.
 * Koon, M. W. Lo, J. E. Marsden, S. D. Ross, Dynamical Systems, the Three-Body Problem and
 * Space Mission Design, 2011). */
#ifndef MYB_CR3BP_H_
#define MYB_CR3BP_H_

#include <stddef.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Acceleration in the rotating frame; Jacobian d a / d x and d a / d v (row-major). */
MYB_API void myb_cr3bp_acceleration(double mu, const double x[3], const double v[3], double a[3]);
MYB_API void myb_cr3bp_jacobian(double mu, const double x[3], double dadx[9], double dadv[9]);
MYB_API double myb_cr3bp_jacobi(double mu, const double y[6]);

/* Terms for a force model (myb/model.h), context = a pointer to mu. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_term(void* mu, double t, const double x[3],
                                                const double v[3], double a[3]);
MYB_API MYB_NODISCARD myb_status myb_cr3bp_term_jacobian(void* mu, double t, const double x[3],
                                                         const double v[3], double dadx[9],
                                                         double dadv[9]);

/* The five Lagrange points (L1 between the primaries, L2 beyond the smaller, L3 beyond the
 * larger, L4 and L5 at +-60 degrees). MYB_E_DOMAIN: mu outside (0, 0.5]. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_lagrange(double mu, double points[5][3]);

typedef enum myb_cr3bp_symmetry {
  MYB_CR3BP_SYMMETRY_XZ = 1, /* y0 = (x, 0, z, 0, vy, 0); at T/2: y = 0, vx = vz = 0 */
  MYB_CR3BP_SYMMETRY_X = 2   /* y0 = (x, 0, 0, 0, vy, vz); at T/2: y = 0, z = 0, vx = 0 */
} myb_cr3bp_symmetry;

typedef struct myb_cr3bp_correction {
  myb_cr3bp_symmetry symmetry;
  int fixed;          /* the component of y0 held fixed: XZ: 0 (x), 2 (z) or 4 (vy);
                         X: 0 (x), 4 (vy) or 5 (vz) */
  int crossing;       /* the half period ends at this crossing of y = 0 (1 = the first) */
  double tolerance;   /* on the constraints; 0 = 1e-12 */
  int max_iterations; /* 0 = 30 */
} myb_cr3bp_correction;

/* Corrects y0 (in place) to a periodic orbit of the given symmetry and returns its period.
 * Errors: MYB_E_DOMAIN (y0 not of the symmetric form, invalid settings), MYB_E_NOCONV (no
 * convergence, or no crossing within 50 time units). */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_correct(double mu,
                                                   const myb_cr3bp_correction* correction,
                                                   double y0[6], double* period);

/* Monodromy matrix (state transition matrix over one period) and its eigenvalues; index =
 * (|lambda_max| + 1 / |lambda_max|) / 2 (1 for a linearly stable orbit), the stability index
 * of the JPL catalogue. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_monodromy(double mu, const double y0[6], double period,
                                                     double monodromy[36]);
MYB_API MYB_NODISCARD myb_status myb_cr3bp_stability(const double monodromy[36], double re[6],
                                                     double im[6], double* index);

/* A state on the unstable (unstable != 0) or stable manifold of the periodic orbit at the
 * fraction phase (0 .. 1) of the period: the orbit's state displaced by epsilon (in length
 * units, along the position part of the eigenvector of the monodromy carried to that phase)
 * on the side sign (+1 or -1). MYB_E_DOMAIN: the orbit is not unstable. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_manifold(double mu, const double y0[6], double period,
                                                    double phase, int unstable, int sign,
                                                    double epsilon, double y[6]);

/* Third-order halo orbit of Richardson about L1 (point = 1) or L2 (2) with out-of-plane
 * amplitude az (in units of the distance of the primaries), branch +1 (northern, z > 0 at
 * the start) or -1: y0 at the xz-plane crossing with vx = vz = 0 (approximately), and its
 * period. MYB_E_DOMAIN: invalid input or an amplitude without a halo orbit. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_richardson(double mu, int point, double az, int branch,
                                                      double y0[6], double* period);

/* ---- Frames ------------------------------------------------------------------------------------
 * The rotating frame of two real bodies: origin at their barycentre, x along body 1 -> body 2,
 * z along their orbital angular momentum, lengths in units of their instantaneous distance
 * and time in units of 1 / n, n = sqrt((gm1 + gm2) / distance^3) (the "pulsating" frame of
 * ephemeris-based libration-point studies; the rotation of the orbital plane is neglected in
 * the velocities). States in m, m/s relative to any common origin. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_from_ephemeris(const double body1[6],
                                                          const double body2[6], double gm1,
                                                          double gm2, const double state[6],
                                                          double rotating[6]);
MYB_API MYB_NODISCARD myb_status myb_cr3bp_to_ephemeris(const double body1[6],
                                                        const double body2[6], double gm1,
                                                        double gm2, const double rotating[6],
                                                        double state[6]);

/* ---- Periodic orbits of any shape, and their families -----------------------------------------
 * Multiple shooting: the orbit is cut into `patches` equal parts in time; the unknowns are
 * the states at the cuts and the period; the conditions are that every part ends where the
 * next begins and the last where the first begins — all components but one, which the Jacobi
 * integral gives: the velocity component that is largest at the first cut — and that the
 * first state stays in the plane across the flow where it started (a phase). That is one
 * condition less than the unknowns: periodic orbits come in one-parameter families. The
 * correction takes minimum-norm Newton steps (myb_shoot of myb/transfer.h), the transition
 * matrices of the parts from the variational equations: the member of the family nearest to
 * the guess, with no symmetry assumed.
 *
 * Continuation by pseudo-arclength (H. B. Keller, 1977; E. J. Doedel et al.): from an orbit
 * and the tangent of its family — the null vector of the Jacobian of the conditions — the
 * next is predicted a step along the tangent and corrected in the plane across it. The
 * family is followed through its turning points in the period or in any coordinate. */

enum { MYB_CR3BP_MAX_PATCHES = 10 };

typedef struct myb_cr3bp_orbit {
  size_t patches;                            /* 1 .. MYB_CR3BP_MAX_PATCHES */
  double states[6 * MYB_CR3BP_MAX_PATCHES];  /* at the times k period / patches */
  double period;
  int has_tangent;                           /* the tangent below is that of this orbit */
  double tangent[6 * MYB_CR3BP_MAX_PATCHES + 1]; /* unit, in (states, period) */
} myb_cr3bp_orbit;

/* The cuts of the orbit through y0 with the given period, by propagation. Errors:
 * MYB_E_DOMAIN (mu outside (0, 0.5], the number of patches, period <= 0), those of the
 * propagation. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_orbit_from_state(double mu, const double y0[6],
                                                            double period, size_t patches,
                                                            myb_cr3bp_orbit* orbit);

/* Corrects the orbit in place until every condition holds to `tolerance` (0: 1e-11), in at
 * most max_iterations steps (0: 30); its tangent is found as well. Errors: MYB_E_DOMAIN,
 * MYB_E_NOCONV, MYB_E_NOMEM, those of the propagation. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_orbit_correct(double mu, myb_cr3bp_orbit* orbit,
                                                         double tolerance, int max_iterations);

/* The next orbit of the family, a distance `step` along the tangent in the space of the
 * unknowns (positive: the way the tangent of the orbit points; the tangent of a corrected
 * orbit keeps its sense from one step to the next, the first one points towards a longer
 * period, or towards a larger first coordinate where the period does not change). The orbit
 * must have been corrected. Errors: as myb_cr3bp_orbit_correct; MYB_E_DOMAIN also for an
 * orbit without a tangent; on MYB_E_NOCONV the orbit is left as it was: take a shorter step. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_orbit_continue(double mu, myb_cr3bp_orbit* orbit,
                                                          double step, double tolerance,
                                                          int max_iterations);

/* ---- Quasi-periodic orbits ----------------------------------------------------------------------
 * About a periodic orbit whose monodromy matrix has a pair of eigenvalues exp(+-i rho) on the
 * unit circle the motion winds on two-dimensional tori — Lissajous orbits about the planar
 * and vertical Lyapunov orbits, quasi-halo orbits, quasi-periodic distant retrograde orbits.
 * A torus is found as an invariant circle of the flow over a time T (G. Gomez, J. M. Mondelo,
 * Physica D 157, 283, 2001; Z. P. Olikara, D. J. Scheeres, AAS 12-361, 2012): the curve that
 * the flow carries in the time T onto itself turned by the angle rho,
 *   R(-rho) phi_T(u) - u = 0,
 * u the curve at `points` equal angles, R the turn of its trigonometric interpolant. T is the
 * time of one turn along the orbit, rho the turn across it in that time.
 *
 * The tori of the problem form families of two parameters, and the equation above has two
 * conditions fewer than it seems: the flow keeps the Jacobi integral and the action of a
 * closed curve, the integral of p . dq round it. Two multipliers take the place of what the
 * two invariants leave undetermined; they vanish at a solution, and what is left of them —
 * `defect` — measures how far the points resolve the circle. Two phase conditions hold the
 * points in place along the circle and along the flow; one condition picks a family of one
 * parameter (the Jacobi constant, the time T or the turn rho kept), and along that family
 * the tori are followed by pseudo-arclength as the periodic orbits above are. Newton's
 * method on all of it, the transition matrices of the points from the variational
 * equations. */

enum { MYB_CR3BP_MAX_CIRCLE = 35 };

typedef struct myb_cr3bp_torus {
  size_t points;                             /* odd, 5 .. MYB_CR3BP_MAX_CIRCLE */
  double states[6 * MYB_CR3BP_MAX_CIRCLE];   /* the circle at the angles 2 pi k / points */
  double period;                             /* T */
  double rotation;                           /* rho (rad) */
  double defect;                             /* the larger of the two multipliers left */
  int has_tangent;
  double tangent[6 * MYB_CR3BP_MAX_CIRCLE + 2]; /* of the family, in (states, T, rho) */
} myb_cr3bp_torus;

typedef enum myb_cr3bp_torus_family {
  MYB_CR3BP_TORUS_JACOBI = 0,  /* the mean Jacobi constant of the circle is kept */
  MYB_CR3BP_TORUS_PERIOD = 1,  /* T is kept */
  MYB_CR3BP_TORUS_ROTATION = 2 /* rho is kept */
} myb_cr3bp_torus_family;

/* The smallest torus about the periodic orbit through y0: the circle of the linear motion,
 * y0 + size (a cos theta - b sin theta), a + i b the eigenvector of the monodromy matrix for
 * exp(i rho), scaled so that `size` is the root mean square distance of the circle from y0
 * (positions and velocities together, in the units of the problem). `pair` picks among the
 * pairs on the unit circle, in the order of rho from 0 to pi (0: the first); the pair at 1
 * that every periodic orbit has does not count. The circle is exact to the order of size^2:
 * myb_cr3bp_torus_continue with a step of 0 corrects it, with a positive one moves outwards.
 * Errors: MYB_E_DOMAIN (mu, the period, the points, size <= 0, no such pair), MYB_E_NOCONV
 * (the eigenvector), those of the propagation. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_torus_from_orbit(double mu, const double y0[6],
                                                            double period, size_t points,
                                                            int pair, double size,
                                                            myb_cr3bp_torus* torus);

/* The torus of the family a distance `step` along its tangent — the root mean square
 * displacement of the circle, with the changes of T and rho — corrected across the tangent
 * until the conditions hold to `tolerance` (0: 1e-10) in at most max_iterations steps (0:
 * 20); a step of 0 corrects the torus where it is. The value the family keeps is that of the
 * torus on entry. The tangent kept in the torus is that of the family of the last step: to
 * change the family take a step of 0 in the new one first, which turns the tangent to it.
 * Errors: MYB_E_DOMAIN (an invalid torus or one without a tangent),
 * MYB_E_NOCONV (the torus is left as it was: take a shorter step, or more points),
 * MYB_E_NOMEM, those of the propagation. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_torus_continue(double mu, myb_cr3bp_torus* torus,
                                                          myb_cr3bp_torus_family family,
                                                          double step, double tolerance,
                                                          int max_iterations);

/* The point of the circle at any angle (rad), by its trigonometric interpolant. */
MYB_API void myb_cr3bp_torus_state(const myb_cr3bp_torus* torus, double angle, double y[6]);

/* ---- Stability of a torus and its manifolds -----------------------------------------------------
 * The motion near a torus is that of a circle near the invariant one under the map that finds
 * it, F(u) = R(-rho) phi_T(u): a displacement d of the points goes to A d,
 *   A = R(-rho) diag(Phi_1 .. Phi_n),
 * Phi_j the transition matrix of the point j over T. In the continuous circle this is the
 * Floquet problem of an invariant curve (A. Jorba, Nonlinearity 14, 943, 2001): with an
 * eigenfunction psi, Phi(theta) psi(theta) = lambda psi(theta + rho), comes psi(theta)
 * exp(i k theta) of the eigenvalue lambda exp(-i k rho), for every k, so that each eigenvalue
 * of the torus is a circle of n eigenvalues of A (the values at the points and the Fourier
 * coefficients of their interpolant are one space in two bases: A is the matrix of the turned
 * linearised map in either). A is symplectic, as R and every Phi_j are: its 6 n eigenvalues
 * come as lambda, 1 / lambda and their conjugates.
 *
 * Of the six circles four are trivial, at the turns exp(i k rho): those of the directions along
 * the circle and along the flow, which the torus carries into themselves, and two of the
 * family of tori (a change of the Jacobi constant or of the action changes T and rho), each a
 * Jordan block of 1 with one of the first two. The last two are the normal behaviour: a pair
 * exp(+-i sigma) on the unit circle (an elliptic torus, about which tori of three frequencies
 * lie), or a real pair lambda, 1 / lambda (a hyperbolic torus, with stable and unstable
 * manifolds of three dimensions; in a phase space of six the normal pair is never complex off
 * the unit circle). The representative of a real pair is its member at k = 0, real since n is
 * odd, and lambda is that of the torus.
 *
 * The points resolve the members of a circle up to some |k| only: beyond it the eigenvalues
 * of A are no longer turns of each other, and some leave the unit circle. A member is
 * resolved when its circle holds more than half the points, a circle being the eigenvalues
 * joined by turns by +-rho within the tolerance; the trivial circles are those joined to 1.
 * The Jordan blocks split by the square root of the error of A: 2.6e-6 radially for the tori
 * below, well within the tolerance of 1e-4.
 *
 * Measured, 25 points: about the L1 Lyapunov orbit of the JPL catalogue (lambda = 2446.324 of
 * its monodromy matrix) Lissajous tori of the sizes 1e-3, 2e-3, 4e-3 (the root mean square
 * distance of the circle from the orbit; 1e-3 is 384 km) have lambda = 2446.191, 2445.794,
 * 2444.210: -5.4e-5, -2.2e-4, -8.6e-4 relative, the size squared (ratios 3.997 and 3.988);
 * every member of the two hyperbolic circles is resolved (50 of 150 eigenvalues), lambda times
 * 1 / lambda is 1 to 4e-10, A v = lambda v to 4e-9; the map itself by another integrator and
 * finite differences gives back lambda v to 3.7e-6; 35 points (a matrix of 210, 0.2 s) give
 * the lambda of 25 to 2e-14. About two L2 halo orbits (stability indices 12.17 and 74.86)
 * quasi-halo tori of size 1e-3 have the lambda of the orbit to 3.5e-6 and 1.8e-6. About a
 * distant retrograde orbit all is on the unit circle: 42 elliptic eigenvalues of 126 (21
 * points), the two circles of the normal pair, the rest trivial or unresolved; at the size
 * 0.02 a few unresolved eigenvalues are 5.6 % off the unit circle and the torus is still
 * elliptic. Manifolds: a point displaced by 1e-9 departs from the torus by lambda in a turn
 * to 3.4e-6 (forward on the unstable manifold, backward on the stable one, against the
 * eigenfunction interpolated at the angle reached); carried to the section of the Moon
 * (x = 1 - mu) with displacements of 1e-6 .. 1e-4, the manifolds of a Lissajous torus cross
 * it where those of its Lyapunov orbit do, out of the plane as the size of the torus and in
 * it as its square (4e-4 in y, vy at the size 1e-3, ratio 3.994 to the size 2e-3), the
 * Jacobi constant kept to 4e-15, and the stable one is the mirror image of the unstable one
 * (the symmetry y, vx, vz, t -> -y, -vx, -vz, -t). */

enum { MYB_CR3BP_MAX_SPECTRUM = 6 * MYB_CR3BP_MAX_CIRCLE };

typedef struct myb_cr3bp_torus_spectrum {
  size_t count;                      /* 6 points */
  double re[MYB_CR3BP_MAX_SPECTRUM]; /* every eigenvalue of A, by modulus from the largest */
  double im[MYB_CR3BP_MAX_SPECTRUM];
  /* the classes of the eigenvalues; count in all */
  size_t trivial;    /* resolved, on the circles joined to 1 */
  size_t elliptic;   /* resolved, on the unit circle, on others */
  size_t hyperbolic; /* resolved, off the unit circle */
  size_t unresolved; /* on circles of no more than half the points */
  double index;      /* (|lambda| + 1 / |lambda|) / 2 of the unstable eigenvalue; 1 without */
  double unstable;   /* the real eigenvalue of the unstable direction; 0: none */
  double stable;     /* that of the stable direction; 0: none */
  /* their eigenfunctions at the points, of unit root mean square over the points, turned so
   * that the mean of the x component is positive (or zero) */
  double unstable_direction[6 * MYB_CR3BP_MAX_CIRCLE];
  double stable_direction[6 * MYB_CR3BP_MAX_CIRCLE];
  double residual; /* |A v - lambda v| / (|lambda| |v|) of the two directions, the larger */
  double tail;     /* the share of their energy above half the highest frequency, the larger */
} myb_cr3bp_torus_spectrum;

/* The spectrum of the torus: its points carried over T with their transition matrices, the
 * eigenvalues of A, their classes, and the eigenvectors of the hyperbolic pair. `tolerance`
 * (0: 1e-4) is the distance within which an eigenvalue counts as on the unit circle, at 1, or
 * a turn of another (relative to its modulus); it must stay well above the splitting of the
 * Jordan blocks and below the turn itself (|1 - exp(i rho)| > 10 tolerance). The work grows
 * as (6 points)^3. On an error the spectrum is not usable. Errors: MYB_E_DOMAIN (an invalid
 * torus or tolerance, a turn rho too near 0), MYB_E_NOCONV (the eigenvalues; or a resolved
 * circle off the unit circle with no real member: no spectrum of a torus of this problem —
 * more points, or another tolerance), MYB_E_NOMEM, those of the propagation. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_torus_stability(double mu, const myb_cr3bp_torus* torus,
                                                           double tolerance,
                                                           myb_cr3bp_torus_spectrum* spectrum);

/* A state on the unstable (unstable != 0) or stable manifold of the torus: the point of the
 * circle at `angle` displaced by epsilon (in length units, along the position part of the
 * eigenfunction there, interpolated as the circle is) on the side sign (+1 or -1; +1 is the
 * side of growing x on the mean over the circle: towards the Moon from a torus about L1). The
 * unstable manifold leaves the torus forward in time, its distance growing by |lambda| in T
 * while it is near; the stable one backward. The displacement is that of the linear manifold,
 * and the circle's own error (1e-11 for the tori above) grows with it: an epsilon well above
 * that (1e-6 .. 1e-4) when the manifold is followed far. MYB_E_DOMAIN: invalid input, a
 * spectrum of another number of points, or no hyperbolic pair. */
MYB_API MYB_NODISCARD myb_status myb_cr3bp_torus_manifold(const myb_cr3bp_torus* torus,
                                                          const myb_cr3bp_torus_spectrum* spectrum,
                                                          double angle, int unstable, int sign,
                                                          double epsilon, double y[6]);

/* ---- Propagation to a time or to a section ---------------------------------------------------
 * The flow of the problem by IAS15 (myb/integ.h, epsilon 1e-12), from y0 for the time
 * `duration` (negative: backward), stopping at the first crossing of a plane when a section
 * is given: y[axis] = value for a component of the state (0 .. 5), crossed with that component
 * growing in time (direction +1), falling (-1) or either (0); a start on the plane does not
 * count. The crossing is located on the dense output of the step by the event locator of the
 * integrators. y (may be y0): the state at the crossing, or at the end; *elapsed (may be
 * NULL) the time from the start; *crossed (may be NULL) 1 if the plane was reached. For the
 * manifolds of periodic orbits and of tori carried to a surface of section. Measured: the L1
 * Lyapunov orbit of the catalogue crosses y = 0 at T / 2 to 1e-15 and at T to 2e-12, the
 * state returning to 2e-12; against DOP853 over 0.7 T, 9e-13. Errors: MYB_E_DOMAIN (mu,
 * duration or y0 not finite, an invalid section), MYB_E_NOCONV (more than 10^6 steps), those
 * of the integrator (MYB_E_STEP_MIN at a collision with a primary). */
typedef struct myb_cr3bp_section {
  int axis; /* 0 .. 5 */
  double value;
  int direction; /* +1, -1 or 0 */
} myb_cr3bp_section;

MYB_API MYB_NODISCARD myb_status myb_cr3bp_propagate(double mu, const double y0[6], double duration,
                                                     const myb_cr3bp_section* section, double y[6],
                                                     double* elapsed, int* crossed);

#ifdef __cplusplus
}
#endif

#endif /* MYB_CR3BP_H_ */
