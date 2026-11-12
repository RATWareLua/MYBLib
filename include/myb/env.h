/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/env.h — the environment of a force model and its ready-made terms.
 *
 * A force model (myb/model.h) sums terms a(t, x, v) with t in seconds from an epoch. The
 * environment gives that t its meaning and holds what the terms share: the bodies around the
 * object, each a set of optional properties and providers (GM, the state from an ephemeris, the
 * orientation, a radius, a luminosity, ...), evaluated once per instant and kept until t
 * changes — integrators evaluate every substep several times, and most terms need the same
 * Sun, the same rotation of the central body. The library still knows no body: every number
 * and every provider is the caller's.
 *
 * Frames and origin. All bodies and the object share one set of inertial axes (those of the
 * state providers, e.g. the ICRF). The origin is the body without a state provider — the
 * central body of the motion, whose own acceleration by the others then enters as indirect
 * terms — or, when every body has a provider, the inertial origin of those providers (a
 * barycentre: no indirect terms).
 *
 * Time. The model time t = 0 is the instant `epoch` on the scale named by `scale`; each
 * provider is called with the same instant on the scale it was registered with (TDB for the
 * JPL ephemerides, TT for the rotation of the Earth). The difference of two scales is read
 * from a table of its values every `table_step` seconds of the model time, at the nodes
 * epoch + k table_step, interpolated on four nodes: a pure function of t (the same bits
 * whatever was asked before), at a small part of the cost of the series of TDB - TT (measured
 * against the exact conversion over 40 days either way at the default step of one hour:
 * 4e-16 s for TDB - TT, 3e-15 s for TCB - TT, a difference of 20 s with a rate).
 *
 * Terms. Each term is a context structure filled by the caller and a function of the form
 * myb_term_function: add them to a model with myb_force_model_add (the context must outlive
 * the model). Parameters to estimate are doubles inside these structures or inside the
 * myb_object they point to: register their addresses with myb_force_model_add_parameter.
 * What a term keeps between calls depends on t only. An environment and the terms using it
 * serve one thread at a time. */
#ifndef MYB_ENV_H_
#define MYB_ENV_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"
#include "myb/force.h"
#include "myb/frames.h"
#include "myb/grav.h"
#include "myb/kin.h"
#include "myb/rel.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- The environment ------------------------------------------------------------------------ */

enum {
  MYB_ENV_MAX_BODIES = 32,
  MYB_ENV_MAX_SCALES = 4, /* time scales other than the model's */
  MYB_ENV_MAX_RAISERS = 8,
  MYB_ENV_MAX_OCCULTERS = 8
};

typedef struct myb_env myb_env;

typedef struct myb_env_options {
  myb_time epoch;             /* the instant of t = 0 */
  const myb_timesys* timesys; /* NULL: one time scale, every provider gets epoch + t */
  const char* scale;          /* the scale of epoch and t, e.g. "TT" (needed with timesys) */
  double table_step;          /* s, > 0; 0 = 3600 */
} myb_env_options;

/* A body as the forces see it. Everything is optional: zero means absent. */
typedef struct myb_env_body {
  int64_t id; /* the caller's identifier (e.g. a NAIF ID), for myb_env_find */
  double gm;  /* m^3 s^-2 */
  /* State relative to the origin of the environment in its inertial axes; NULL: the body is
   * the origin. state_scale names the scale of the provider's argument (NULL: the model's). */
  myb_state_function state;
  void* state_context;
  const char* state_scale;
  /* Orientation: inertial axes -> body-fixed frame; NULL: the body-fixed frame is the inertial
   * axes. orientation_scale as above (NULL: the model's). */
  myb_orientation_function orientation;
  void* orientation_context;
  const char* orientation_scale;
  double radius;     /* m: the sphere that casts shadows and bounds the source's disk */
  double luminosity; /* W: the body as a source of radiation */
  double limb[2];    /* limb darkening u1, u2 of the source's disk (myb_shadow_fraction) */
  double spin[3];    /* angular momentum per unit mass in the body-fixed frame, m^2/s */
} myb_env_body;

/* The time system and every provider context must outlive the environment. Errors:
 * MYB_E_DOMAIN (a step that is negative, not finite or above 1e15 s, timesys without a scale
 * name),
 * MYB_E_NOT_FOUND (the scale), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_env_create(const myb_env_options* options,
                                                const myb_allocator* allocator, myb_env** env);
MYB_API void myb_env_destroy(myb_env* env); /* accepts NULL */

/* A copy with the same bodies and providers and values of its own: an environment keeps
 * what it computed for the last instant asked, so one object serves one thread, and a
 * second thread takes a copy (the providers and their contexts are shared and must
 * themselves be safe to call from several threads, as the ready-made ones on kernels opened
 * for reading are). Errors: MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_env_clone(const myb_env* env,
                                               const myb_allocator* allocator, myb_env** clone);

/* Adds a body (the description is copied) and returns its index. Errors: MYB_E_DOMAIN (a
 * value that is negative or not finite, a second body without a state provider, a scale name
 * without a time system), MYB_E_NOT_FOUND (a scale name), MYB_E_OVERFLOW (more than
 * MYB_ENV_MAX_BODIES bodies or MYB_ENV_MAX_SCALES scales). */
MYB_API MYB_NODISCARD myb_status myb_env_add_body(myb_env* env, const myb_env_body* body,
                                                  size_t* index);

MYB_API size_t myb_env_body_count(const myb_env* env);
MYB_API double myb_env_gm(const myb_env* env, size_t body);
MYB_API myb_time myb_env_epoch(const myb_env* env);

/* The first body with the identifier id. Errors: MYB_E_NOT_FOUND. */
MYB_API MYB_NODISCARD myb_status myb_env_find(const myb_env* env, int64_t id, size_t* index);

/* The instant epoch + t on the scale `scale` (NULL: the model's). Errors: MYB_E_DOMAIN (t not
 * finite), MYB_E_RANGE (|t| beyond 1e17 s or 2^52 table steps), MYB_E_NOT_FOUND (a scale that
 * is neither the model's nor one of a provider). */
MYB_API MYB_NODISCARD myb_status myb_env_time(myb_env* env, double t, const char* scale,
                                              myb_time* out);

/* State of a body relative to the origin (zeros for the body at the origin) and its
 * orientation at the model time t, as the terms see them. Errors: MYB_E_DOMAIN (t or the
 * provider's result not finite), MYB_E_RANGE (as myb_env_time), the statuses of the
 * providers. */
MYB_API MYB_NODISCARD myb_status myb_env_state(myb_env* env, double t, size_t body,
                                               double state[6]);
MYB_API MYB_NODISCARD myb_status myb_env_orientation(myb_env* env, double t, size_t body,
                                                     myb_orientation* orientation);

/* ---- The object ------------------------------------------------------------------------------
 * What the surface forces need of the body whose motion is modelled, treated as a sphere. */
typedef double (*myb_mass_function)(void* context, double t);

typedef struct myb_object {
  double mass; /* kg, > 0 */
  double area; /* m^2: the cross-section */
  double cr;   /* radiation pressure coefficient, 1 + reflectivity */
  double cd;   /* drag coefficient */
  /* a mass that changes: if given, the terms take mass_function(mass_context, t) (kg) at the
   * model time t instead of mass — myb_thrust_mass_at with a myb_thrust (myb/thrust.h) for an
   * object that burns propellant, myb_ablation_mass with a myb_ablation (myb/surface.h) for
   * one that the heat of an entry burns away */
  myb_mass_function mass_function;
  void* mass_context;
} myb_object;

/* ---- Terms -----------------------------------------------------------------------------------
 * Every function below is a myb_term_function (the acceleration is written, not added) with
 * its context structure as the first argument; x and v are relative to the origin of the
 * environment. Common errors: MYB_E_DOMAIN (the object at the centre of a body, invalid values
 * in the structure), the statuses of the providers. */

/* A point mass at the origin: a = -gm x / |x|^3, with its analytic Jacobian. */
typedef struct myb_point_mass_term {
  double gm;
} myb_point_mass_term;

MYB_API MYB_NODISCARD myb_status myb_term_point_mass(void* term, double t, const double x[3],
                                                     const double v[3], double a[3]);
MYB_API MYB_NODISCARD myb_status myb_term_point_mass_jacobian(void* term, double t,
                                                              const double x[3],
                                                              const double v[3], double dadx[9],
                                                              double dadv[9]);

/* The point mass of a body of the environment: with a body at the origin the perturbation
 * direct minus indirect (myb_force_third_body, without cancellation), otherwise the direct
 * attraction. Analytic Jacobian. MYB_E_DOMAIN also for the body at the origin itself. */
typedef struct myb_third_body_term {
  myb_env* env;
  size_t body;
} myb_third_body_term;

MYB_API MYB_NODISCARD myb_status myb_term_third_body(void* term, double t, const double x[3],
                                                     const double v[3], double a[3]);
MYB_API MYB_NODISCARD myb_status myb_term_third_body_jacobian(void* term, double t,
                                                              const double x[3],
                                                              const double v[3], double dadx[9],
                                                              double dadv[9]);

/* The gravity field of an extended body (myb_gravity_field_eval, degrees n_min .. n_max;
 * n_max = -1: all) in its body-fixed frame, by the body's orientation. For a body that is not
 * at the origin of an environment with one, the same field at the origin is subtracted (the
 * indirect term).
 *
 * With update != 0 the coefficients follow time (myb_gravity_field_set_time: trends and
 * periodic terms of the file) and with tides != NULL the tidal changes are added
 * (myb_tides_apply) raised by the bodies `raisers`, with the rotation time and the pole of the
 * body's orientation: both once per instant, at the instant on the scale of the body's
 * orientation (TT for the Earth). The term modifies the field: it needs its own
 * (myb_gravity_field_clone), and a work for the degrees used. */
typedef struct myb_field_term {
  myb_env* env;
  size_t body;
  myb_gravity_field* field;
  myb_gravity_work* work;
  int n_min;
  int n_max;
  int update;
  const myb_tides* tides;
  int tide_flags; /* MYB_TIDES_* */
  myb_permanent_tide permanent;
  size_t raisers[MYB_ENV_MAX_RAISERS];
  size_t raiser_count;
  int cached;      /* state: zero-initialise (and zero again after changing the fields above) */
  double cached_t;
} myb_field_term;

MYB_API MYB_NODISCARD myb_status myb_term_field(void* term, double t, const double x[3],
                                                const double v[3], double a[3]);
/* What the figure of the body at the origin adds to the motion about it when other bodies
 * pull on that figure: each of them (the raisers of the term) attracts the oblateness of the
 * body and so moves its centre, the origin of the coordinates, by
 *   a = sum (GM_j / GM_body) g(r_j),
 * g the field of the body without its central part (degrees from 2) at the position of the
 * body j. For the Earth and the Moon it is 1e-11 m/s^2. The term belongs to the same
 * myb_field_term as myb_term_field of the body at the origin, and does not depend on the
 * object: its Jacobian is myb_term_no_jacobian. MYB_E_DOMAIN: the body is not at the
 * origin. */
MYB_API MYB_NODISCARD myb_status myb_term_field_indirect(void* term, double t, const double x[3],
                                                         const double v[3], double a[3]);

/* Analytic: the gradient of the field (myb_gravity_field_gradient) turned into the inertial
 * axes. */
MYB_API MYB_NODISCARD myb_status myb_term_field_jacobian(void* term, double t, const double x[3],
                                                         const double v[3], double dadx[9],
                                                         double dadv[9]);

/* Radiation pressure of a source (a body with a luminosity) on the object as a sphere:
 *   a = nu L / (4 pi d^2 c) cr (area / mass) u,
 * d and u the distance and the direction from the source, nu the product of the fractions of
 * its disk left visible by the occulters (myb_shadow_fraction with the radii of the bodies and
 * the limb darkening of the source; where several airless occulters cover parts of the disk
 * at once, myb_shadow_fraction_many: their overlap counts once). c: the speed of light
 * (m/s). MYB_E_DOMAIN also for an occulter or, with occulters, a source without a radius. */
typedef struct myb_radiation_term {
  myb_env* env;
  size_t source;
  const myb_object* object;
  double c;
  size_t occulters[MYB_ENV_MAX_OCCULTERS];
  size_t occulter_count;
  int moving; /* != 0: with the terms of the first order in the velocity relative to the source
                 over c (myb_force_radiation_moving): the Poynting-Robertson drag */
  const myb_refraction* atmospheres; /* NULL, or one per occulter: its air bends and dims the
                                        light (myb_shadow_fraction_refracted; scale_height 0:
                                        no air) */
} myb_radiation_term;

MYB_API MYB_NODISCARD myb_status myb_term_radiation(void* term, double t, const double x[3],
                                                    const double v[3], double a[3]);
/* Analytic in full light and in the umbra; in a penumbra, where the visible fraction changes
 * over tens of kilometres, central differences of the term itself (steps of 6e-6 of the
 * distance from the origin). */
MYB_API MYB_NODISCARD myb_status myb_term_radiation_jacobian(void* term, double t,
                                                             const double x[3],
                                                             const double v[3], double dadx[9],
                                                             double dadv[9]);

/* Sunlight reflected and heat emitted by a planet (myb_force_planet_radiation, the model of
 * Knocke et al.) lit by the source, in the planet's body-fixed frame by its orientation; the
 * flux of the source at the planet from its luminosity. phase0: the phase of the seasonal term
 * at t = 0 (rad), growing at model->omega; rings: 1 .. 50. */
typedef struct myb_planet_radiation_term {
  myb_env* env;
  size_t planet;
  size_t source;
  const myb_planet_radiation* model;
  const myb_object* object;
  double phase0;
  int rings;
} myb_planet_radiation_term;

MYB_API MYB_NODISCARD myb_status myb_term_planet_radiation(void* term, double t,
                                                           const double x[3], const double v[3],
                                                           double a[3]);
/* The gradient of the field of that radiation (myb_force_planet_radiation_gradient) turned
 * into the inertial axes. */
MYB_API MYB_NODISCARD myb_status myb_term_planet_radiation_jacobian(void* term, double t,
                                                                    const double x[3],
                                                                    const double v[3],
                                                                    double dadx[9],
                                                                    double dadv[9]);

/* Drag in the atmosphere of a body (myb_force_drag): the density at the body-fixed position
 * and the instant on the scale of the body's orientation; the atmosphere rotates with the
 * body, so the relative velocity is that in the body-fixed frame (v - omega x r for a uniform
 * rotation). MYB_E_DOMAIN also for a density that is negative or not finite. */
typedef struct myb_drag_term {
  myb_env* env;
  size_t body;
  myb_density_function density;
  void* density_context;
  const myb_object* object;
  myb_wind_function wind; /* may be NULL: the gas at rest in the body-fixed frame */
  void* wind_context;
} myb_drag_term;

MYB_API MYB_NODISCARD myb_status myb_term_drag(void* term, double t, const double x[3],
                                               const double v[3], double a[3]);
/* Analytic in the velocity and in the rotation of the atmosphere; the gradients of the density
 * and of the wind by central differences of their functions over 1e-7 of the distance from
 * the body's centre (about 1e-9 of the gradient for a scale height of tens of kilometres). */
MYB_API MYB_NODISCARD myb_status myb_term_drag_jacobian(void* term, double t, const double x[3],
                                                        const double v[3], double dadx[9],
                                                        double dadv[9]);

/* Relativity of the motion about the body `central` (IERS Conventions 2010, eq. 10.12; the
 * functions of myb/rel.h with the position and velocity relative to that body): its
 * Schwarzschild field, the Lense-Thirring term of its spin (the body's spin through its
 * orientation), and the de Sitter term of its motion in the field of the body `external`. */
typedef struct myb_relativity_term {
  myb_env* env;
  myb_ppn ppn;
  size_t central;
  size_t external; /* de Sitter only */
} myb_relativity_term;

MYB_API MYB_NODISCARD myb_status myb_term_schwarzschild(void* term, double t, const double x[3],
                                                        const double v[3], double a[3]);
MYB_API MYB_NODISCARD myb_status myb_term_lense_thirring(void* term, double t, const double x[3],
                                                         const double v[3], double a[3]);
MYB_API MYB_NODISCARD myb_status myb_term_de_sitter(void* term, double t, const double x[3],
                                                    const double v[3], double a[3]);
/* Analytic (the Jacobians of myb/rel.h). */
MYB_API MYB_NODISCARD myb_status myb_term_schwarzschild_jacobian(void* term, double t,
                                                                 const double x[3],
                                                                 const double v[3],
                                                                 double dadx[9], double dadv[9]);
MYB_API MYB_NODISCARD myb_status myb_term_lense_thirring_jacobian(void* term, double t,
                                                                  const double x[3],
                                                                  const double v[3],
                                                                  double dadx[9],
                                                                  double dadv[9]);
MYB_API MYB_NODISCARD myb_status myb_term_de_sitter_jacobian(void* term, double t,
                                                             const double x[3],
                                                             const double v[3], double dadx[9],
                                                             double dadv[9]);

/* The non-gravitational acceleration of a small body (myb_force_nongrav, myb/force.h) with
 * the position and the velocity relative to the body `source` — the Sun. */
typedef struct myb_nongrav_term {
  myb_env* env;
  size_t source;
  myb_nongrav model;
  double delay; /* s: the activity lags the heating (myb_force_nongrav_delayed with the GM of
                   the source); 0: none */
} myb_nongrav_term;

MYB_API MYB_NODISCARD myb_status myb_term_nongrav(void* term, double t, const double x[3],
                                                  const double v[3], double a[3]);
/* Analytic; with a delay, central differences of the term. */
MYB_API MYB_NODISCARD myb_status myb_term_nongrav_jacobian(void* term, double t,
                                                           const double x[3], const double v[3],
                                                           double dadx[9], double dadv[9]);

/* The Yarkovsky acceleration of a body lit by `source`, its spin axis (unit) fixed in the
 * inertial axes: the diurnal part (myb_force_yarkovsky) and, with `season`, the seasonal one
 * (myb_force_yarkovsky_season; made for the orbit of the body about the source). */
typedef struct myb_yarkovsky_term {
  myb_env* env;
  size_t source;
  myb_yarkovsky body;
  double spin[3];
  double c;
  const myb_yarkovsky_season* season; /* NULL: the diurnal part alone */
} myb_yarkovsky_term;

MYB_API MYB_NODISCARD myb_status myb_term_yarkovsky(void* term, double t, const double x[3],
                                                    const double v[3], double a[3]);
/* Analytic (myb_force_yarkovsky_gradient, myb_force_yarkovsky_season_gradient). */
MYB_API MYB_NODISCARD myb_status myb_term_yarkovsky_jacobian(void* term, double t,
                                                             const double x[3],
                                                             const double v[3], double dadx[9],
                                                             double dadv[9]);

/* The lagging tide that the object raises on a body of the environment
 * (myb_force_tidal_lag): the motion of a moon about its planet. The rotation of the body is
 * that of its orientation; gm_object is the GM of the object. */
typedef struct myb_tidal_lag_term {
  myb_env* env;
  size_t body;
  double gm_object;
  double k2;
  double lag;      /* s */
  int with_static; /* != 0: also the tide without delay */
} myb_tidal_lag_term;

MYB_API MYB_NODISCARD myb_status myb_term_tidal_lag(void* term, double t, const double x[3],
                                                    const double v[3], double a[3]);
/* Analytic. */
MYB_API MYB_NODISCARD myb_status myb_term_tidal_lag_jacobian(void* term, double t,
                                                             const double x[3],
                                                             const double v[3], double dadx[9],
                                                             double dadv[9]);

/* Empirical accelerations in the radial, transverse and normal directions of the orbit (x, v)
 * as orbit analysis estimates them: a constant and a once-per-revolution part,
 *   a_rtn = constant + cosine cos u + sine sin u,
 * u the argument of latitude, counted in the orbital plane from its ascending node on the xy
 * plane of the inertial axes (from the x axis for an orbit in that plane). The nine numbers
 * are parameters. MYB_E_DOMAIN: x and v parallel or zero. */
typedef struct myb_empirical_term {
  double constant[3];
  double cosine[3];
  double sine[3];
} myb_empirical_term;

MYB_API MYB_NODISCARD myb_status myb_term_empirical(void* term, double t, const double x[3],
                                                    const double v[3], double a[3]);
/* Analytic. */
MYB_API MYB_NODISCARD myb_status myb_term_empirical_jacobian(void* term, double t,
                                                             const double x[3],
                                                             const double v[3], double dadx[9],
                                                             double dadv[9]);

#ifdef __cplusplus
}
#endif

#endif /* MYB_ENV_H_ */
