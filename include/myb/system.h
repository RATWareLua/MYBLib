/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/system.h — a system of bodies with their models: the point masses of myb/dyn.h joined by
 * what the bodies are — a figure, a rotation, an answer to tides — as further levels of the
 * interaction graph, with the rotation of chosen bodies integrated together with the motion.
 *
 * The point masses, with their graph, the 1PN terms and the threads, stay with the myb_nbody
 * inside (myb_system_nbody). On top of them every level is asked for by name:
 *
 *   field      the figure of a body acts on another body taken as a point, which acts back on
 *              the body (momentum is kept) and turns it (myb_system_add_field);
 *   figures    the quadrupoles of two bodies act on each other (myb_system_add_figures);
 *   tide       the tide that a body raises on another acts on a third, or on itself
 *              (myb_system_add_tide): Love numbers and time lags of the long-period, diurnal
 *              and semi-diurnal bands;
 *   distortion a body with an integrated rotation is distorted by the tide of another and by
 *              its own spin (myb_system_add_distortion): its inertia tensor, and with it its
 *              degree-2 field, follow;
 *   PN pair    two bodies close enough for the orders beyond 1PN (myb_system_add_pn_pair);
 *   action     anything else that pushes a body or turns it, given as a function
 *              (myb_system_add_action); ready-made: light on the mesh of a body, with the
 *              torque it gives at the moment (myb_system_light_action).
 *
 * The body-fixed frame of a body is either given from outside (myb_system_set_orientation: any
 * provider of myb/frames.h) or integrated (myb_system_set_rotation): a mantle that may hold a
 * fluid core, the model of the lunar ephemerides (W. M. Folkner, J. G. Williams, D. H. Boggs,
 * R. S. Park, P. Kuchynka, IPN Progress Report 42-196, 2014, sections III and IV; R. S. Park
 * et al., Astron. J. 161, 105, 2021).
 *
 * Integration. myb_system_acceleration has the signature of myb_accel_function (myb/integ.h)
 * for n = myb_system_coordinates(system) coordinates in blocks of three, so that IAS15
 * integrates the motion and the rotation as one system (myb_system_integrator makes one, with
 * the positions, every rotation and every pair of spins as groups of their own: a quaternion
 * is not measured in metres):
 *   3 per body          position (m); velocity (m/s);
 *   9 per rotation      the attitude quaternion q of myb/attitude.h (4 values, then 2 unused);
 *                       3 values whose rate is the angular velocity of the core in the mantle
 *                       frame. The velocities are dq/dt and that angular velocity. The
 *                       quaternion obeys q'' = (q' w + q w') / 2 with w = 2 conj(q) q', the
 *                       angular velocity of the mantle in its own frame;
 *   6 per PN pair       with spins: values whose rates are the two dimensionless spins.
 * myb_system_pack_rotation and its relatives fill and read these parts of a state.
 *
 * A system is used by one thread at a time. */
#ifndef MYB_SYSTEM_H_
#define MYB_SYSTEM_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/dyn.h"
#include "myb/frames.h"
#include "myb/grav.h"
#include "myb/integ.h"
#include "myb/surface.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct myb_system myb_system;

enum {
  MYB_SYSTEM_MAX_ROTATIONS = 16,
  MYB_SYSTEM_MAX_FIELDS = 256,   /* pairs (body, point) */
  MYB_SYSTEM_MAX_FIGURES = 32,   /* pairs of figures */
  MYB_SYSTEM_MAX_TIDES = 64,     /* triples (body, raiser, point) */
  MYB_SYSTEM_MAX_DISTORTIONS = 32,
  MYB_SYSTEM_MAX_PN_PAIRS = 16,
  MYB_SYSTEM_MAX_ACTIONS = 64
};

/* `bodies` >= 1 massive bodies with gm[0 .. bodies-1] (finite, >= 0). Errors: MYB_E_DOMAIN,
 * MYB_E_OVERFLOW, MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_system_create(size_t bodies, const double* gm,
                                                   const myb_allocator* allocator,
                                                   myb_system** system);
MYB_API void myb_system_destroy(myb_system* system); /* accepts NULL */

MYB_API size_t myb_system_bodies(const myb_system* system);
/* The point masses: set the graph, the relativity and the threads there. */
MYB_API myb_nbody* myb_system_nbody(myb_system* system);

/* ---- The model of a body --------------------------------------------------------------------- */

/* The body-fixed frame from outside: orientation(context, t, &o) gives the rotation from the
 * inertial frame of the system at the time t of the integration (seconds, as
 * myb_time_from_seconds). The provider must outlive the system. Errors: MYB_E_DOMAIN (the
 * body, a body whose rotation is integrated). */
MYB_API MYB_NODISCARD myb_status myb_system_set_orientation(myb_system* system, size_t body,
                                                            myb_orientation_function orientation,
                                                            void* context);

/* The figure of a body as harmonics in its frame, used up to `degree` (0: all of the field).
 * The field must outlive the system and its GM should be that of the body. For a body whose
 * rotation is integrated the degrees from 3 are taken from the field and the degree 2 from
 * the inertia tensor. Errors: MYB_E_DOMAIN (the body, a degree above the field), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_system_set_field(myb_system* system, size_t body,
                                                      const myb_gravity_field* field, int degree);

/* An integrated rotation. Every moment is per unit mass of the body (m^2). The inertia tensor
 * of the mantle in its frame is
 *   I_m = diag(A - A_c, B - A_c, C - C_c)
 *         - k2 R^5 sum_T (m_T / m) (r r^T - |r|^2 1 / 3) / |r|^5            (tides)
 *         + k2 R^5 / (3 G m) (w w^T - (|w|^2 - n^2) 1 / 3 - n^2 z z^T)      (spin),
 * r the position of the distorting body T (myb_system_add_distortion) at t - lag in the mantle
 * frame of t - lag, w the angular velocity of the mantle, n the mean motion that defines the
 * undistorted body. The rates of both parts enter the equations: that of the tides from the
 * motion, that of the spin through w' itself (the equation is solved for it). The core has the
 * moments (A_c, A_c, C_c) in the mantle frame and is coupled to the mantle by
 *   N = k_v (w_c - w) + (C_c - A_c) (z . w_c) (z x w_c). */
typedef struct myb_system_rotation {
  double moments[3];    /* A <= B <= C of the whole undistorted body (m^2), > 0 */
  double radius;        /* R of the distortion (m); needed with k2 != 0 */
  double k2;            /* potential Love number of degree 2 (0: rigid) */
  double lag;           /* time lag of the distortion (s) */
  double mean_motion;   /* n (rad/s) */
  double core[2];       /* A_c, C_c (m^2); both 0: no core */
  double core_friction; /* k_v (m^2/s) */
} myb_system_rotation;

/* Errors: MYB_E_DOMAIN (the body, a body with an orientation from outside or without mass,
 * invalid values), MYB_E_OVERFLOW (more than MYB_SYSTEM_MAX_ROTATIONS). */
MYB_API MYB_NODISCARD myb_status myb_system_set_rotation(myb_system* system, size_t body,
                                                         const myb_system_rotation* rotation);

/* The answer of a body to tides, by band (order j = 0 long-period, 1 diurnal, 2 semi-diurnal,
 * about the pole of the body — the third axis of its frame). A body T at r_T relative to the
 * body raises the potential, at r,
 *   U = G m_T R^5 / (2 |r|^5) sum_j k2[j] W_j(r, r*_j) / |r*_j|^5,
 *   W_0 = (3 z^2 - |r|^2) (3 z*^2 - |r*|^2) / 2,   W_1 = 6 z z* (rho . rho*),
 *   W_2 = 3 (2 (rho . rho*)^2 - rho^2 rho*^2) / 2       (z along the pole, rho across it),
 * which is k2 G m_T R^5 P_2(cos) / (|r|^3 |r*|^3) when the three bands answer alike. The raiser
 * is taken where it was: r*_j = R_pole(rotation_rate rotation_lag[j]) r_T(t - orbital_lag[j]),
 * a turn about the pole by the angle that the body rotates in the lag (Folkner et al. 2014,
 * eqs. 31-32). Earlier positions come from the motion at t: a Taylor series to the third
 * order, the jerk being that of the two bodies alone. */
typedef struct myb_system_tide {
  double radius;          /* R (m), > 0 */
  double k2[3];           /* Love numbers of the bands */
  double orbital_lag[3];  /* s; may be negative */
  double rotation_lag[3]; /* s; [0] is not used: the long-period band does not turn */
  double rotation_rate;   /* rad/s */
} myb_system_tide;

/* The body must have a frame (an orientation or an integrated rotation) by the time the
 * accelerations are asked for. Errors: MYB_E_DOMAIN. */
MYB_API MYB_NODISCARD myb_status myb_system_set_tide(myb_system* system, size_t body,
                                                     const myb_system_tide* tide);

/* ---- The levels of the graph ----------------------------------------------------------------- */

/* The figure of `body` acts on `point`: the acceleration of the point in the field, its
 * reaction on the body, and the torque on a body whose rotation is integrated. The body
 * needs a field or an integrated rotation. Errors: MYB_E_DOMAIN (indices, the same body
 * twice, a body without mass), MYB_E_OVERFLOW. */
MYB_API MYB_NODISCARD myb_status myb_system_add_field(myb_system* system, size_t body,
                                                      size_t point);

/* The quadrupoles of two bodies on each other (myb_figure_figure): accelerations of both and
 * the torques on those whose rotation is integrated. Errors: as above. */
MYB_API MYB_NODISCARD myb_status myb_system_add_figures(myb_system* system, size_t a, size_t b);

/* The tide that `raiser` raises on `body` (myb_system_set_tide) acts on `point`, which acts
 * back on the body; the raiser and the point may be one. Errors: as above. */
MYB_API MYB_NODISCARD myb_status myb_system_add_tide(myb_system* system, size_t body,
                                                     size_t raiser, size_t point);

/* `raiser` distorts `body`, whose rotation is integrated. Errors: as above. */
MYB_API MYB_NODISCARD myb_status myb_system_add_distortion(myb_system* system, size_t body,
                                                           size_t raiser);

/* The terms beyond 1PN of two bodies (myb_pn_binary_evaluate of myb/rel.h: `orders` of
 * MYB_PN_2 ... MYB_PN_SS, r0 for MYB_PN_3, kappa for MYB_PN_SS), added to what the point
 * masses give: the relative acceleration is shared by the masses, a_a += m_b / (m_a + m_b) a,
 * a_b -= m_a / (m_a + m_b) a, which keeps the Newtonian momentum. Ask for MYB_PN_1 only if
 * the point masses do not carry the 1PN terms of this pair already (its level in the
 * myb_nbody is Newtonian): it would be counted twice. With MYB_PN_SO or MYB_PN_SS the two
 * spins are integrated. The speed of light is `c`, or that of myb_nbody_set_relativity when
 * `c` is 0. Errors: MYB_E_DOMAIN (indices, a body without mass, unknown orders, no speed of
 * light, MYB_PN_3 without r0, a negative kappa), MYB_E_OVERFLOW. */
typedef struct myb_system_pn_pair {
  unsigned orders;
  double c;      /* m/s; 0: that of the relativity of the point masses */
  double r0;     /* m */
  double kappa1; /* of body a */
  double kappa2; /* of body b */
} myb_system_pn_pair;

MYB_API MYB_NODISCARD myb_status myb_system_add_pn_pair(myb_system* system, size_t a, size_t b,
                                                        const myb_system_pn_pair* pair);

/* What the levels above do not give — light, thrust, a gas. The function is called once an
 * evaluation for its body and writes (both are zero when it is called) the acceleration of the
 * body (m/s^2, inertial axes) and the torque on it per unit of its mass (m^2/s^2, in its own
 * frame; used only if its rotation is integrated). */
typedef struct myb_system_view {
  size_t body;         /* the body acted on */
  size_t bodies;
  const double* x;     /* positions of all the bodies, 3 each (m) */
  const double* v;     /* velocities (m/s) */
  const double* gm;    /* as they are now (mergers change them) */
  const double* frame; /* the rotation inertial -> body (9 values, rows), NULL without a frame */
} myb_system_view;

typedef myb_status (*myb_system_action)(void* context, double t, const myb_system_view* view,
                                        double acceleration[3], double torque[3]);

/* Errors: MYB_E_DOMAIN (the body), MYB_E_OVERFLOW (more than MYB_SYSTEM_MAX_ACTIONS). */
MYB_API MYB_NODISCARD myb_status myb_system_add_action(myb_system* system, size_t body,
                                                       myb_system_action action, void* context);

/* Light of the body `source` on the mesh of a body (myb_force_radiation_mesh of
 * myb/surface.h): the force on its centre and the torque about it as they are at the moment,
 * the momentum of the light absorbed included — what turns a body from one instant to the
 * next, where myb_force_yorp gives the mean over a rotation of the re-emitted part alone. The
 * mesh is given in the frame of the body about its centre of mass; the heat is re-emitted at
 * once (surface.thermal), without conduction. An action: pass the structure as the context.
 * Errors: MYB_E_DOMAIN (a body without a frame, invalid values), those of the mesh. */
typedef struct myb_system_light {
  size_t source;               /* the body that shines */
  const myb_polyhedron* shape; /* of the lit body */
  myb_plate surface;           /* its optical properties (area and normal are not used) */
  double flux;                 /* W m^-2 at `distance` from the source */
  double distance;             /* m */
  double c;                    /* m/s */
  double mass;                 /* of the lit body, kg */
  int shadows;                 /* as myb_force_radiation_mesh */
} myb_system_light;

MYB_API MYB_NODISCARD myb_status myb_system_light_action(void* light, double t,
                                                         const myb_system_view* view,
                                                         double acceleration[3],
                                                         double torque[3]);

/* ---- The state ------------------------------------------------------------------------------- */

/* Coordinates of the state: 3 per body, 9 per integrated rotation, 6 per PN pair with spins;
 * a state y holds the coordinates and then their rates, 2 n values. */
MYB_API size_t myb_system_coordinates(const myb_system* system);

/* The rotation of `body` into a state y and out of it: the attitude q (unit; scaled to unit
 * on the way out), the angular velocity of the mantle and that of the core (may be NULL: no
 * core), both in the mantle frame (rad/s). Errors: MYB_E_DOMAIN (the rotation of the body is
 * not integrated, q = 0). */
MYB_API MYB_NODISCARD myb_status myb_system_pack_rotation(const myb_system* system, size_t body,
                                                          const double q[4],
                                                          const double omega[3],
                                                          const double omega_core[3], double* y);
MYB_API MYB_NODISCARD myb_status myb_system_unpack_rotation(const myb_system* system,
                                                            size_t body, const double* y,
                                                            double q[4], double omega[3],
                                                            double omega_core[3]);

/* The spins of the PN pair number `pair` (in the order of myb_system_add_pn_pair), chi1 of its
 * body a and chi2 of its body b. Errors: MYB_E_DOMAIN (no such pair, a pair without spins). */
MYB_API MYB_NODISCARD myb_status myb_system_pack_spins(const myb_system* system, size_t pair,
                                                       const double chi1[3],
                                                       const double chi2[3], double* y);
MYB_API MYB_NODISCARD myb_status myb_system_unpack_spins(const myb_system* system, size_t pair,
                                                         const double* y, double chi1[3],
                                                         double chi2[3]);

/* An IAS15 integrator for the system as it stands (levels that add coordinates — rotations,
 * pairs with spins — must be in place): n coordinates, blocks of three, and the groups of
 * myb_ias15_options — the positions of all the bodies, each rotation, each pair of spins —
 * so that each sets the step by its own timescale. Of `options` (may be NULL) the precision,
 * the step limits and the budget are used. Start it with myb_system_acceleration and the
 * system. Errors: those of myb_ias15_create. */
MYB_API MYB_NODISCARD myb_status myb_system_integrator(const myb_system* system,
                                                       const myb_ias15_options* options,
                                                       const myb_allocator* allocator,
                                                       myb_ias15** integrator);

/* The second derivatives of the coordinates (myb_accel_function; the context is the system).
 * Errors: those of myb_nbody_acceleration_xv, of the orientation providers and of the fields;
 * MYB_E_DOMAIN (a level that needs a frame the body does not have, coincident bodies, an
 * inertia tensor that cannot be inverted). */
MYB_API MYB_NODISCARD myb_status myb_system_acceleration(void* system, double t, const double* x,
                                                         const double* v, double* a);

/* The inertia tensor of the whole body (mantle and core, m^2, in the mantle frame) of an
 * integrated rotation at the state (t, y): what the degree 2 of its field is taken from.
 * Errors: as myb_system_acceleration. */
MYB_API MYB_NODISCARD myb_status myb_system_inertia(myb_system* system, size_t body, double t,
                                                    const double* y, double inertia[9]);

/* ---- Encounters ------------------------------------------------------------------------------
 * Bodies with a radius (myb_system_set_radius) are spheres that may touch. myb_system_advance
 * steps an IAS15 integrator of the system towards t_end and stops at the first contact — the
 * instant at which the distance of two such bodies falls to the sum of their radii, found on
 * the polynomial of the step and made the end of the step — or, if `watch` > 0, at the
 * closest approach of two bodies that pass within `watch` of each other. What happens then
 * is the caller's: myb_system_merge, myb_system_bounce, or going on. */

MYB_API MYB_NODISCARD myb_status myb_system_set_radius(myb_system* system, size_t body,
                                                       double radius);
MYB_API double myb_system_radius(const myb_system* system, size_t body);
/* 0 for a body that has been merged into another. */
MYB_API int myb_system_present(const myb_system* system, size_t body);

typedef enum myb_system_event_kind {
  MYB_SYSTEM_ARRIVED = 0,  /* t_end reached */
  MYB_SYSTEM_CONTACT = 1,  /* the surfaces of a and b touch, approaching */
  MYB_SYSTEM_APPROACH = 2  /* a and b are at their closest, within `watch` */
} myb_system_event_kind;

typedef struct myb_system_event {
  myb_system_event_kind kind;
  size_t a;
  size_t b;        /* a < b */
  double t;        /* the time of the integrator */
  double distance; /* of the centres (m) */
} myb_system_event;

/* The integrator must have been started with myb_system_acceleration and this system. It is
 * left at event->t. A pair already in contact when the call starts (after an event that was
 * left alone) is not reported again until it has separated. Errors: those of myb_ias15_step;
 * MYB_E_DOMAIN (watch negative or not finite). */
MYB_API MYB_NODISCARD myb_status myb_system_advance(myb_system* system, myb_ias15* integrator,
                                                    double t_end, double watch,
                                                    myb_system_event* event);

/* What a merger takes out of the motion, with G as a factor like the invariants of
 * myb/dyn.h: the energy of the relative motion, G m_a m_b (|v|^2 / (2 (m_a + m_b)) - 1 / |r|),
 * and its angular momentum, G m_a m_b / (m_a + m_b) r x v, for r = x_b - x_a, v = v_b - v_a.
 * They go into heat and spin, which the system does not follow. */
typedef struct myb_system_merger {
  double energy;              /* G times the energy of the relative motion (m^5 s^-4) */
  double angular_momentum[3]; /* G times its angular momentum (m^5 s^-3) */
} myb_system_merger;

/* Merges body b into body a at the current state of the integrator: a takes the mass of both,
 * the position and velocity of their centre of mass (the momentum is kept) and the radius of
 * a sphere of both volumes; it keeps its own model. b is absent from then on: it attracts
 * nothing, feels nothing and its part of the state means nothing. Every level that names b is
 * dropped. The integrator is started again from the new state. Errors: MYB_E_DOMAIN (indices,
 * an absent body, a body whose rotation is integrated or that is part of a PN pair),
 * the errors of myb_ias15_start. */
MYB_API MYB_NODISCARD myb_status myb_system_merge(myb_system* system, myb_ias15* integrator,
                                                  size_t a, size_t b,
                                                  myb_system_merger* merger);

/* Two hard spheres rebound: the component of the relative velocity along the line of centres
 * is reversed and multiplied by `restitution` (1: elastic, 0: they stay together until
 * something parts them), the momentum is kept. The integrator is started again. Errors:
 * MYB_E_DOMAIN (indices, restitution outside 0 .. 1, coincident centres). */
MYB_API MYB_NODISCARD myb_status myb_system_bounce(myb_system* system, myb_ias15* integrator,
                                                   size_t a, size_t b, double restitution);

#ifdef __cplusplus
}
#endif

#endif /* MYB_SYSTEM_H_ */
