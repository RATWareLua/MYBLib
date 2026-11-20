/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/contact.h — contact of a body with the surface of a shape: the landing, the rebound,
 * sliding, and coming to rest.
 *
 * The body touches the surface with contact points: spheres fixed in its frame, each an offset
 * from its centre of mass and a radius — one sphere about the centre for a ball or a probe,
 * the pads of the legs of a lander (a radius of 0 is a point). A sphere presses into the
 * surface by the depth d = r - h, h the height of its centre over the surface (the distance to
 * the nearest point of the surface, negative inside), and the surface pushes back along the
 * outward normal there.
 *
 * The normal force is that of K. H. Hunt and F. R. E. Crossley (J. Appl. Mech. 42, 440, 1975):
 *   f_n = k d^n (1 + chi d'),   d' the rate of the depth,
 * an elastic force (n = 3/2: two bodies of Hertz; n = 1: a linear spring) with a damping that
 * grows with the depth. It is preferred to the linear spring and dashpot k d + c d' because
 * its force is continuous at both ends of the contact — zero where the surfaces touch, so that
 * an integrator meets no jump — and because it never pulls: a dashpot pushes with c d' the
 * instant the surfaces meet and pulls before they part, and holding its force at zero there
 * changes the restitution it was tuned for. And its restitution has a closed form. Per unit
 * mass, along the normal, d'' = -k d^n (1 + chi d'): v dv / (1 + chi v) = -k d^n dd, which
 * integrates over the contact, from d = 0 back to d = 0, to
 *   chi v_in - ln(1 + chi v_in) = -chi v_out - ln(1 - chi v_out),     e = v_out / v_in,
 * whatever k and n (the relation of D. W. Marhefka, D. E. Orin, IEEE Trans. Syst. Man
 * Cybern. A 29, 566, 1999). e falls from 1 at low speeds (1 - e = 2 chi v_in / 3 to the first
 * order) to 1 / (chi v_in) at high ones; the separating speed stays below 1 / chi, so the
 * force never changes sign and a contact of this law ends exactly where the surfaces part.
 * myb_contact_restitution and myb_contact_damping turn one into the other. A force that would
 * still be negative (a point thrown off faster than 1 / chi by something else) is held at
 * zero: there is no adhesion.
 *
 * Friction is Coulomb's, made smooth across the stick:
 *   f_t = -mu f_n u / sqrt(|u|^2 + v_s^2),   u the sliding velocity of the point,
 * the whole mu f_n for |u| >> v_s (short of it by (v_s / u)^2 / 2), linear in u below v_s. A
 * switch between sticking and sliding needs a state of its own carried through the
 * integration — the stretch of a bristle (P. R. Dahl, 1968; C. Canudas de Wit et al., IEEE
 * Trans. Autom. Control 40, 419, 1995) — where these forces are functions of the state alone,
 * as a term of a force model is; and of the smooth forms this one is smooth at every order (a
 * ramp min(|u| / v_s, 1) has a kink at v_s that costs an integrator of high order its order).
 * The price: on a slope of angle a below the angle of friction a body does not stick but
 * creeps, at v_s t / sqrt(1 - t^2), t = tan(a) / mu; and near rest the friction is stiff, mu
 * f_n / v_s per unit of speed, which sets the step there.
 *
 * Friction that sticks (MYB_CONTACT_FRICTION_STICK): a spring between the point and an anchor
 * on the surface, which slides when the spring would pull harder than the limit — the
 * elasto-plastic element of Prandtl, and the "spring with a slip limit" of the discrete
 * element method (P. A. Cundall, O. D. L. Strack, Geotechnique 29, 47, 1979). Its state is the
 * deflection z of each point, the point less its anchor, carried with the integrated state;
 * the force on the body is
 *   f_t = -P_L(k_t z + c_t u),   L = mu f_n,   P_L(F) = F min(1, L / |F|),
 * the spring and a dashpot, never more than the limit; and the deflection follows the point
 *   z' = u - phi zhat max(0, zhat . u) - zhat max(0, k_t |z| - 2 L) / c_t,
 *   phi = 0 below 0.95 of the limit, rising linearly to 1 at it (and beyond),
 * z and u the tangential parts, of the point at the surface. Below 0.95 L the anchor does not
 * move: a body at rest on a slope below the angle of friction stays where it is, held by
 * springs deflected by its weight along the slope (g sin a / k_t of a single point) — no
 * creep, whatever the time. At the limit the anchor is dragged with the point and the force
 * is mu f_n along the spring, in steady sliding along the velocity: Coulomb's law exactly, not
 * to (v_s / u)^2 / 2. When the normal force falls (a point that lifts off) the force is held
 * to the limit by P_L and goes to zero with it; a spring left stretched beyond twice the limit
 * is drawn back at the rate k_t / c_t (between once and twice the limit it keeps its stretch,
 * the force being the limit: a point that stops springs back by two slip deflections at
 * most). A point off the surface keeps its deflection; it is released against the limit,
 * which starts from zero, when the point touches again.
 * Why a band and not a switch at the limit itself: a leg that slides as it presses in sees
 * its limit grow under it, and a deflection switched between holding (z' = u, out) and
 * sliding (z' tangent to the limit) at the limit chatters on it — the steps of DOP853 fell to a
 * nanosecond at the first touch of the lander below. With the band the rate is a continuous
 * function of the state. Why not draw an overstretched spring back to the limit itself: a
 * steady slide would then sit at the top of the band, where the rate changes over 0.05 L / k_t
 * of deflection, stiff (|u| k_t / (0.05 L)) for a point lightly loaded. And why not a band
 * reaching beyond the limit: the spring would store energy past the cap faster than the
 * friction spends it — measured, the energy of the sliding lander rose by 2.6e-5.
 * Energy. With V = sum k_t |z|^2 / 2 counted with the other energies, the friction gives
 * dE/dt = -(the work against the anchor) - (the work of the dashpot) <= 0 in every case — below
 * the limit -c_t |u|^2 at most, at it, with the anchor dragged, -L |zhat . u| and less: any drag
 * of the anchor towards the point spends energy, and at the limit the spring never stretches
 * on — so that a contact never adds energy with this friction either. (On a curved surface
 * the tangent plane turns under a point that slides and the tangential part of z with it: the
 * statement is exact on a plane, which a point that sticks sees anyway.)
 * The rate is not smooth — the band has corners, max(0, zhat . u) one where the point moves
 * square to its spring, and P_L one at the limit — but it and the force are continuous: DOP853
 * goes through, measuring its error. Keep c_t > 0: it
 * damps the springs where the anchor holds and draws back a spring overstretched. The state
 * has no place in a term of a force model, which is a function of the state of the object
 * alone: myb_term_contact refuses such a law; it is a level of a system
 * (myb_system_add_contact), whose state carries 3 values per point for it, or
 * myb_contact_force_stick for an integrator of the caller's.
 *
 * Energy. With E the kinetic energy plus the potential of the outside forces plus the energy
 * of the springs, sum k d^(n+1) / (n + 1),
 *   dE/dt = -sum (k chi d^n d'^2 + mu f_n |u|^2 / sqrt(|u|^2 + v_s^2)) <= 0:
 * a contact never adds energy, whatever the motion (the depth is the distance to the surface,
 * whose gradient is the normal, so k d^n is the force of the potential of the springs even on
 * an edge of a mesh).
 *
 * The surface: a plane, a triaxial ellipsoid (a sphere) or a closed triangular mesh (myb/grav.h)
 * at rest in its own frame. On a mesh the nearest point may be inside a face, on an edge or on
 * a corner (myb_polyhedron_nearest, myb/geometry.h), and the normal is then the line from it to
 * the centre of the sphere: a ball that comes down on an edge is thrown back along the line
 * from the edge to its centre. A point meets the surface once, at its nearest point: a sphere
 * caught in a concave corner of a mesh feels the nearer wall alone.
 *
 * Three ways in:
 *   myb_contact_force       the force and the torque on a body in a given state, for any
 *                           integrator;
 *   myb_term_contact        a term of a force model (myb/model.h): the object of myb_orbit as a
 *                           sphere about its centre, on the surface of a body of an
 *                           environment, which may turn — the translation alone. Declared as a
 *                           break (myb_force_model_add_break), myb_event_contact ends the steps
 *                           where a contact begins and ends;
 *   myb_system_add_contact  (myb/system.h) a body of a system on the surface of another, with the
 *                           torques on those whose rotation is integrated: a lander on its legs.
 *
 * Integration: with DOP853 (myb_orbit with MYB_ORBIT_DOP853; for a system,
 * myb_system_acceleration as a system of the first order under myb_dop853 of myb/integ.h),
 * which measures its error against the state and goes through everything a contact does; or
 * with IAS15 in myb_orbit, given three things. (1) The timescale MYB_IAS15_TIMESCALE_JERK
 * (myb_orbit_options.timescale): a contact that holds a body against gravity passes the
 * acceleration through zero twice in every bounce, and the timescale of Pham, Rein and
 * Spiegel, relative to |a|, rejects every step that would cross such a zero — a ball dropped
 * on a plane stops it with MYB_E_STEP_MIN within its first contact; the same option rejects
 * the steps whose sweeps diverge or stall (a step of free flight, which IAS15 makes as long as
 * it likes, plunging into the stiff surface; a block resting on its spring, whose fast normal
 * motion the timescale does not see) and does not shorten a step for an acceleration that is
 * the rounding of a balance (myb/integ.h). (2) The breaks at the touch and the parting
 * (myb_force_model_add_break with myb_event_contact), so that no step holds the kink of a
 * touch; either scheme does (myb_force_model_set_break_scheme): MYB_BREAKS_PREDICT ends the
 * step of the fall at the touch predicted from the flight before it is taken, the default
 * withdraws it — with the same results and steps to a percent (below). (3) Some other force
 * on the body during the contact, gravity: where nothing else acts, the acceleration leaves
 * zero as a power of the time at both ends of a contact, where no timescale is defined.
 *
 * Measured (tests/unit/test_contact.c, DOP853 at 1e-13): a ball (Hertz, k = 1e9, chi = 0.5 s/m)
 * thrown at a plane without gravity comes back with the restitution of the closed form to
 * 3.5e-13; dropped on it under 1.62 m/s^2 it is off by -2.2e-4, and under a tenth of that
 * gravity at the same speed by -2.2e-5 (ratio 10.003): the impulse of the gravity during the
 * contact, of the order g tau / v, which the closed form leaves out. Dropped 0.2 m on the Moon,
 * a ball bounces 24 times and rests at the depth (g / k)^(1/n) to 2e-16 m, the energy never
 * up by more than 1e-16 m^2/s^2 from one millisecond to the next. A block sliding at 2 m/s
 * with mu = 0.4 follows the closed form of the regularised law to 2e-15 s, halves its speed
 * 6e-8 s after Coulomb's block and slides 0.509685 m against its 0.509684 m. On a cube (a
 * mesh) the rebound from a face is that from a plane to 5e-13, from an edge and from a corner
 * it goes back along the line to the ball to 1.2e-15 with the restitution of the closed form
 * to 1e-12. A ball on the equator of a turning asteroid (a term with an environment) is
 * carried along by the friction to the speed of the ground, to 2e-11 m/s. A lander on four
 * legs (a system, its rotation integrated) settles on a slope of 15 degrees with mu = 0.6 and
 * creeps at the speed of the closed form to 6e-13, and slides down one of 40 degrees at
 * g (sin a - mu cos a) to 9e-8; its energy never grows. A spinning ball rubbing past a turning
 * ellipsoid of its own system leaves the momentum of the two to 2e-16 and their angular
 * momentum, orbits and spins, to 1.2e-13. The nearest point of a mesh: a cube to its closed
 * form exactly, a star-shaped mesh of 528 faces to 4e-16 of a search of every face, and inside
 * and outside as the solid angle says at 3000 points.
 * Friction that sticks (k_t = 1e4 s^-2, c_t = 100 s^-1, mu = 0.6, g = 1.62 m/s^2; a system
 * under DOP853): a block let go on 25 degrees rests at the deflection of its weight, 6.8464e-5
 * m, to 1e-9 of it (its overshoot, 0.90 of the limit, stays below the band), moves 7.5e-18 m in
 * the last 30 s of a minute, and rests at the same place to 1.6e-18 m at rtol 1e-9 .. 1e-12;
 * the slope turned slowly to a part in a thousand below the angle of friction, it holds, 1.11
 * slip deflections down (5.7e-15 m in its last 40 s); a part above, it slides at
 * g (sin a - mu cos a) to 2.5e-13. The lander below (k_t = 1000 s^-2, c_t = 16 s^-1 a leg) on
 * 15 degrees stands: 9.3e-14 m/s after 20 s where the smooth law creeps at 5.0e-4 m/s, 1.6e-14
 * m in its last 10 s, its legs bearing its weight to 1e-11, at the same place to 1.1e-9 m at
 * rtol 1e-9 .. 1e-11 (3402 steps, the smooth law 7653); on 40 degrees it slides at
 * g (sin a - mu cos a) to 8.3e-8 (80 854 steps: the direction of a spring of 0.2 mm follows a
 * velocity of metres a second). The energy never grows by more than 4.4e-16 m^2/s^2.
 * Under IAS15 (epsilon 1e-9, 1e-10, 1e-12, the timescale of a', the touches predicted) against
 * DOP853 at rtol 1e-12: the ball dropped 0.2 m on the Moon touches 24 times, at 7.0e-8, 1.4e-8,
 * 7.6e-9 s of DOP853's touches, rests at its height to 1.5e-16 m, the energy never up by more
 * than 1.1e-16, in 13 789, 18 778, 34 922 steps (DOP853: 6729; the touches withdrawn:
 * 5.0e-8, 2.6e-8, 3.3e-9 s in 14 039, 19 024, 35 180 steps); the sliding block reaches its
 * three speeds at 7.9e-11, 1.0e-11, 9.1e-14 s of DOP853's (the closed form to as much) in
 * about 2300 steps (DOP853: 489). With the timescale of a alone the first contact stops IAS15
 * with MYB_E_STEP_MIN. */
#ifndef MYB_CONTACT_H_
#define MYB_CONTACT_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/env.h"
#include "myb/grav.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- The surface ---------------------------------------------------------------------------- */

typedef enum myb_contact_shape {
  MYB_CONTACT_PLANE = 0,     /* point, normal */
  MYB_CONTACT_ELLIPSOID = 1, /* radii: the semi-axes along x, y, z, centre at the origin */
  MYB_CONTACT_MESH = 2       /* mesh: closed, in the frame of the surface */
} myb_contact_shape;

typedef struct myb_contact_surface {
  myb_contact_shape shape;
  double point[3];            /* plane: a point of it (m) */
  double normal[3];           /* plane: towards the side the body is on; any length */
  double radii[3];            /* ellipsoid (m), > 0; equal for a sphere */
  const myb_polyhedron* mesh; /* mesh */
} myb_contact_surface;

/* The point of the surface nearest to a position, the outward unit normal there (of a mesh:
 * as myb_polyhedron_nearest) and the height of the position over it (m, negative inside); any
 * output may be NULL. Errors: MYB_E_DOMAIN (an invalid surface, a position that is not
 * finite), MYB_E_NOCONV (the ellipsoid, myb_ellipsoid_nearest). */
MYB_NODISCARD MYB_API myb_status myb_contact_surface_nearest(const myb_contact_surface* surface,
                                                             const double position[3],
                                                             double point[3], double normal[3],
                                                             double* height);

/* ---- The law ------------------------------------------------------------------------------- */

/* The friction (above). */
typedef enum myb_contact_friction {
  MYB_CONTACT_FRICTION_SMOOTH = 0, /* Coulomb's made smooth: a function of the state */
  MYB_CONTACT_FRICTION_STICK = 1   /* a spring with a slip limit: a state of its own */
} myb_contact_friction;

/* Every force is per unit mass of the body (an acceleration), so k and the other numbers are
 * those of the body as a whole: for the legs of a lander of mass m, the stiffness of one leg
 * in N m^-n over m. */
typedef struct myb_contact_law {
  double stiffness; /* k, > 0: m^(1-n) s^-2 */
  double exponent;  /* n >= 1: 1 a linear spring, 1.5 the sphere of Hertz */
  double damping;   /* chi >= 0 (s/m); 0: elastic */
  double friction;  /* mu >= 0 */
  double slip;      /* v_s (m/s) of the smooth friction, > 0 where mu > 0 */
  myb_contact_friction kind;   /* 0: smooth */
  double tangential_stiffness; /* k_t (s^-2) of the friction that sticks, > 0 where mu > 0 */
  double tangential_damping;   /* c_t (1/s) of the friction that sticks, >= 0 */
} myb_contact_law;

/* A contact point: a sphere fixed in the body, its centre at `offset` from the centre of mass
 * in the body frame (m), of radius >= 0 (m). */
typedef struct myb_contact_point {
  double offset[3];
  double radius;
} myb_contact_point;

typedef struct myb_contact_result {
  double force[3];  /* m s^-2: the force per unit mass of the body, in the axes of the surface */
  double torque[3]; /* m^2 s^-2: per unit mass, about the centre of mass, in the body axes */
  double energy;    /* m^2 s^-2: the energy of the springs, sum k d^(n+1) / (n + 1), and with
                       the friction that sticks sum k_t |z|^2 / 2 over every point */
  double gap;       /* m: the least height of a sphere over the surface less its radius,
                       negative while one presses into it — its sign changes are the contacts */
  size_t touching;  /* the points that press into the surface */
} myb_contact_result;

/* The force and the torque of the surface on a body: x and v of its centre of mass in the
 * frame of the surface (m, m/s), attitude the rotation from the axes of the surface to the
 * body axes (9 values, rows: the body axes; NULL: the same axes) and omega the angular velocity
 * of the body relative to the surface in the body axes (rad/s; NULL: none). Each point that
 * presses in is pushed at the nearest point of the surface, with the velocity of that point of
 * the body, v + omega x (s - x). Errors: MYB_E_DOMAIN (an invalid surface, law or point, no
 * points, a state that is not finite, a friction that sticks — it needs its state:
 * myb_contact_force_stick), MYB_E_NOCONV (an ellipsoid). */
MYB_NODISCARD MYB_API myb_status myb_contact_force(const myb_contact_surface* surface,
                                                   const myb_contact_law* law,
                                                   const myb_contact_point* points, size_t count,
                                                   const double x[3], const double v[3],
                                                   const double attitude[9], const double omega[3],
                                                   myb_contact_result* result);

/* As myb_contact_force, with the state of a friction that sticks: deflection[3 k .. 3 k + 2]
 * the deflection z of point k (m, in the axes of the surface), and deflection_rate (may be
 * NULL) its rate z', to be integrated with the motion — of a point off the surface 0. With the
 * smooth friction both are ignored (may be NULL) and the result is that of myb_contact_force.
 * Errors: as myb_contact_force (but for the friction that sticks), MYB_E_DOMAIN for a
 * deflection that is not finite or missing. */
MYB_NODISCARD MYB_API myb_status myb_contact_force_stick(
    const myb_contact_surface* surface, const myb_contact_law* law,
    const myb_contact_point* points, size_t count, const double x[3], const double v[3],
    const double attitude[9], const double omega[3], const double* deflection,
    double* deflection_rate, myb_contact_result* result);

/* The restitution of the law, e = v_out / v_in, of a contact met at the speed v_in (m/s)
 * along the normal, with no other force during it: the relation above, solved by Newton's
 * iteration (series where chi v is small). damping 0 or speed 0: 1. Errors: MYB_E_DOMAIN (a
 * value negative or not finite). */
MYB_NODISCARD MYB_API myb_status myb_contact_restitution(double damping, double speed,
                                                         double* restitution);

/* The damping chi (s/m) that gives the restitution e at the approach speed v (m/s): 0 for
 * e = 1. Errors: MYB_E_DOMAIN (e outside (0, 1], v not positive, values not finite). */
MYB_NODISCARD MYB_API myb_status myb_contact_damping(double restitution, double speed,
                                                     double* damping);

/* ---- A term of a force model ------------------------------------------------------------------
 * The object as a sphere of `radius` about its centre, on the surface of a body of an
 * environment given in the body-fixed frame of the body about its centre: the object is taken
 * into that frame with the state and the orientation of the body (myb_env_state,
 * myb_env_orientation, with the rate of the rotation), so that a turning body carries the
 * surface along and the friction acts on the velocity relative to it. With env == NULL the
 * surface is fixed at the origin of the model, in its axes. The body is not pushed back. No
 * torque: the translation alone (a body that turns: myb_system_add_contact of myb/system.h). */
typedef struct myb_contact_term {
  myb_env* env;                /* NULL: a surface fixed at the origin */
  size_t body;                 /* of env */
  myb_contact_surface surface; /* in the body-fixed frame of the body */
  myb_contact_law law;         /* per unit mass of the object */
  double radius;               /* of the object (m), >= 0 */
} myb_contact_term;

/* myb_term_function. Errors: those of myb_contact_force (a friction that sticks among them)
 * and of the environment. */
MYB_NODISCARD MYB_API myb_status myb_term_contact(void* term, double t, const double x[3],
                                                  const double v[3], double a[3]);

/* The height of the object over the surface less its radius (m; myb_contact_result.gap): its
 * fall through zero is the touch, its rise the parting. An event function (myb_orbit_event)
 * and a break of a force model (myb_force_model_add_break) with the same term as its context;
 * NaN where the environment fails. */
MYB_API double myb_event_contact(void* term, double t, const double* y);

#ifdef __cplusplus
}
#endif

#endif /* MYB_CONTACT_H_ */
