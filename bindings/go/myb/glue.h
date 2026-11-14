/* SPDX-License-Identifier: BSD-3-Clause */
/* glue.h — the C half of the Go binding: what one crossing from Go into C does as a whole.
 *
 * A force model is described by plain records (bodies, terms, parameters) that hold numbers
 * and pointers to objects of the library only — never to memory of Go — and turned into an
 * environment and a model inside one call; a propagation, a search for events or a fit is
 * another single call over arrays. Nothing here computes: the records are copied into the
 * structures of the library and its functions are called in the order a C program would call
 * them, so a number that comes back has the bits of the same calls made in C.
 *
 * Declarations only: the Go file that exports the callbacks includes this header in its
 * preamble, which cgo allows to hold no definitions. */
#ifndef MYBGO_GLUE_H_
#define MYBGO_GLUE_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/myb.h"

/* ---- Layout of the records Go shares with C ------------------------------------------------
 * The Go types Time, Calendar, State, Elements, Equinoctial, Rotation, EOPValues and
 * Orientation are laid out as these C structures; the Go side checks the offsets against
 * these constants when it compiles. */
enum {
  MYBGO_SIZE_TIME = sizeof(myb_time),
  MYBGO_SIZE_CALENDAR = sizeof(myb_calendar),
  MYBGO_OFFSET_CALENDAR_MONTH = offsetof(myb_calendar, month),
  MYBGO_OFFSET_CALENDAR_SECOND = offsetof(myb_calendar, second),
  MYBGO_SIZE_PV = sizeof(myb_pv),
  MYBGO_SIZE_ELEMENTS = sizeof(myb_elements),
  MYBGO_SIZE_EQUINOCTIAL = sizeof(myb_equinoctial),
  MYBGO_SIZE_ROTATION = sizeof(myb_rotation),
  MYBGO_SIZE_EOP = sizeof(myb_eop_values),
  MYBGO_SIZE_ORIENTATION = sizeof(myb_orientation),
  MYBGO_OFFSET_ORIENTATION_UT = offsetof(myb_orientation, ut),
  MYBGO_OFFSET_ORIENTATION_XP = offsetof(myb_orientation, xp)
};

/* ---- The description of a model -------------------------------------------------------------- */

enum {
  MYBGO_NAME = 32, /* names of scales, terms and parameters, with the NUL */

  MYBGO_PROVIDER_NONE = 0,   /* no state: the body is the origin of the environment */
  MYBGO_PROVIDER_SPK = 1,    /* a body of an SPK set (myb_spk_body_state) */
  MYBGO_PROVIDER_KEPLER = 2, /* the two-body orbit through a state at an epoch */
  MYBGO_PROVIDER_GO = 3,     /* a Go function, called back through cgo */

  MYBGO_ORIENTATION_NONE = 0,   /* the body-fixed frame is the inertial axes */
  MYBGO_ORIENTATION_FIXED = 1,  /* a constant rotation */
  MYBGO_ORIENTATION_CIP = 2,    /* the CIO-based model with observed EOP */
  MYBGO_ORIENTATION_PCK = 3,    /* a frame of a binary PCK */
  MYBGO_ORIENTATION_WGCCRE = 4, /* a rotation model of the WGCCRE */
  MYBGO_ORIENTATION_CK = 5      /* the pointing of a CK through its clock (myb_orientation_ck) */
};

/* A provider of states (myb_state_function of myb/kin.h) described by data. A KEPLER provider
 * gives myb_kepler_propagate(mu, state, the seconds from epoch to t): the conic of the problem
 * of two bodies that a C program would write for a planet on its orbit. */
typedef struct mybgo_provider {
  int kind; /* MYBGO_PROVIDER_* */
  const myb_spk* spk;
  int32_t target;
  int32_t center;
  int32_t frame;
  double mu;
  double state[6];
  myb_time epoch;
  uintptr_t handle;
} mybgo_provider;

/* myb_state_function; context: const mybgo_provider* (kept where it is while it is used).
 * Errors: MYB_E_DOMAIN (no provider, an SPK provider without a set), those of the source. */
myb_status mybgo_provider_state(void* context, myb_time t, double state[6]);

/* The states of a provider at count instants (the array form of a C loop over it). */
myb_status mybgo_provider_states(const mybgo_provider* provider, const myb_time* t, size_t count,
                                 double* states, size_t* failed);

/* An orientation provider described by data (myb/frames.h). */
typedef struct mybgo_orientation {
  int kind;           /* MYBGO_ORIENTATION_* */
  myb_rotation fixed; /* FIXED */
  myb_rotation base;  /* PCK, WGCCRE: inertial axes -> the base of the frame */
  int has_base;       /* 0: the base is the inertial axes */
  const myb_cip* cip; /* CIP: the model, the EOP, the time system and the name of TT */
  const myb_eop* eop;
  const myb_timesys* timesys;
  const char* tt_scale;
  const myb_pck* pck; /* PCK: the set and the frame class */
  int32_t pck_id;
  const myb_wgccre* wgccre; /* WGCCRE */
  const myb_ck* ck;         /* CK: the set, the clock, the instrument, the base frame wanted */
  const myb_sclk* sclk;
  int32_t instrument;
  int32_t ck_base;
  double tolerance; /* ticks */
} mybgo_orientation;

/* Where the contexts of an orientation provider live while it is used. */
typedef struct mybgo_orientation_store {
  myb_cip_orientation cip;
  myb_pck_orientation pck;
  myb_wgccre_orientation wgccre;
  myb_ck_orientation ck;
  myb_rotation fixed;
  myb_rotation base;
} mybgo_orientation_store;

/* The provider of a record and its context in `store`; *function NULL for NONE. Errors:
 * MYB_E_DOMAIN (an unknown kind, a missing object), those of finding the scale TT. */
myb_status mybgo_orientation_of(const mybgo_orientation* orientation,
                                mybgo_orientation_store* store, myb_orientation_function* function,
                                void** context);

/* The orientation at count instants (myb_orientation_many). */
myb_status mybgo_orientation_sample(const mybgo_orientation* orientation, const myb_time* t,
                                    size_t count, myb_orientation* out, size_t* failed);

typedef struct mybgo_body {
  int64_t id;
  double gm;
  double radius;
  double luminosity;
  double limb[2];
  double spin[3];
  mybgo_provider state;          /* NONE: the body is the origin */
  const char* state_scale;       /* NULL: the model's scale */
  mybgo_orientation orientation; /* NONE: the inertial axes */
  const char* orientation_scale; /* NULL: the model's scale */
} mybgo_body;

enum {
  MYBGO_TERM_POINT_MASS = 1,
  MYBGO_TERM_THIRD_BODY = 2,
  MYBGO_TERM_FIELD = 3,
  MYBGO_TERM_DRAG = 4,
  MYBGO_TERM_RADIATION = 5,
  MYBGO_TERM_SCHWARZSCHILD = 6,
  MYBGO_TERM_LENSE_THIRRING = 7,
  MYBGO_TERM_DE_SITTER = 8,
  MYBGO_TERM_EMPIRICAL = 9,
  MYBGO_TERM_GO = 10 /* a Go function, called back through cgo */
};

typedef struct mybgo_term {
  int kind; /* MYBGO_TERM_* */
  const char* name;
  size_t body;     /* THIRD_BODY, FIELD, DRAG: the body; RADIATION: the source; the relativity
                      terms: the central body */
  size_t external; /* DE_SITTER */
  double gm;       /* POINT_MASS */
  myb_gravity_field* field; /* FIELD: not modified unless update != 0, and then copied first */
  int n_min;
  int n_max;
  int update;
  myb_object object;                     /* DRAG, RADIATION (mass_function is ignored) */
  myb_exponential_atmosphere atmosphere; /* DRAG */
  double c;                              /* RADIATION */
  size_t occulters[MYB_ENV_MAX_OCCULTERS];
  size_t occulter_count;
  int moving;
  int breaks;         /* RADIATION: the edges of the shadows declared as breaks of the model */
  myb_ppn ppn;        /* the relativity terms */
  double constant[3]; /* EMPIRICAL */
  double cosine[3];
  double sine[3];
  uintptr_t handle;  /* GO: the Go function's handle */
  int position_only; /* declared a function of the place alone (myb_force_model_set_position_only)
                      */
} mybgo_term;

enum {
  MYBGO_PARAMETER_CD = 1,       /* the drag coefficient of a DRAG term */
  MYBGO_PARAMETER_CR = 2,       /* the radiation coefficient of a RADIATION term */
  MYBGO_PARAMETER_CONSTANT = 3, /* EMPIRICAL: constant[index] */
  MYBGO_PARAMETER_COSINE = 4,   /* EMPIRICAL: cosine[index] */
  MYBGO_PARAMETER_SINE = 5      /* EMPIRICAL: sine[index] */
};

typedef struct mybgo_parameter {
  const char* name;
  size_t term; /* the index of the term in the description */
  int what;    /* MYBGO_PARAMETER_* */
  int index;   /* 0 .. 2 for the empirical ones */
  double step; /* of the central differences */
} mybgo_parameter;

typedef struct mybgo_model_options {
  myb_time epoch;
  const myb_timesys* timesys; /* NULL: one time scale */
  const char* scale;
  double table_step;
} mybgo_model_options;

typedef struct mybgo_model mybgo_model;

/* Where a description failed: what the status of mybgo_model_create refers to. */
enum {
  MYBGO_STAGE_ENVIRONMENT = 1,
  MYBGO_STAGE_BODY = 2,
  MYBGO_STAGE_TERM = 3,
  MYBGO_STAGE_PARAMETER = 4
};

/* The environment and the force model of a description. The objects of the library it points
 * to must outlive the model; everything else is copied. Errors: those of the library, with
 * *stage (MYBGO_STAGE_*) and *item (the index of the body, term or parameter) naming the
 * record; MYB_E_DOMAIN also for an index of a body that the description does not have, an
 * unknown kind, a missing object. */
myb_status mybgo_model_create(const mybgo_model_options* options, const mybgo_body* bodies,
                              size_t body_count, const mybgo_term* terms, size_t term_count,
                              const mybgo_parameter* parameters, size_t parameter_count,
                              mybgo_model** model, int* stage, size_t* item);
void mybgo_model_destroy(mybgo_model* model); /* accepts NULL */

/* The values of the parameters now (they change in a fit), and new ones (MYB_E_DOMAIN if not
 * finite). */
size_t mybgo_model_parameter_count(const mybgo_model* model);
void mybgo_model_get_parameters(const mybgo_model* model, double* values);
myb_status mybgo_model_set_parameters(mybgo_model* model, const double* values);

/* What the tasks of the other parts of the glue take from a model: its force model, its
 * environment and the count of its bodies. */
myb_force_model* mybgo_model_force(mybgo_model* model);
myb_env* mybgo_model_env(mybgo_model* model);
size_t mybgo_model_bodies(const mybgo_model* model);

/* The sum of the enabled terms at count points: a[3 k ..] at t[k] and the state y[6 k ..].
 * *failed: as for the array forms of the library. */
myb_status mybgo_model_accelerations(mybgo_model* model, const double* t, const double* y,
                                     size_t count, double* a, size_t* failed);

/* ---- Propagation ------------------------------------------------------------------------------
 */

/* From y0 at t0 towards t_end: the states at the count times (in the order of the propagation;
 * states: 6 each), with options->variational the transition matrices (phi: 36 each) and the
 * sensitivities (6 np each) — either may be NULL —, and the impulses kicks[4 j ..] = (t, dv)
 * applied on the way (at their times, in the same order; a state asked at the time of an
 * impulse is the one after it). *steps: the accepted steps when it ended; *failed: the index
 * of the time that failed (count if none did). */
myb_status mybgo_propagate(mybgo_model* model, const myb_orbit_options* options, double t0,
                           const double* y0, double t_end, const double* times, size_t count,
                           double* states, double* phi, double* sensitivity, const double* kicks,
                           size_t kick_count, int64_t* steps, size_t* failed);

/* ---- Events ------------------------------------------------------------------------------------
 */

enum {
  MYBGO_EVENT_PLANE = 1,      /* myb_event_plane */
  MYBGO_EVENT_RANGE = 2,      /* myb_event_range */
  MYBGO_EVENT_RANGE_RATE = 3, /* myb_event_range_rate: the apsides */
  MYBGO_EVENT_ALTITUDE = 4,   /* myb_event_altitude */
  MYBGO_EVENT_ELEVATION = 5,  /* myb_event_elevation over a station */
  MYBGO_EVENT_ECLIPSE = 6,    /* myb_event_eclipse */
  MYBGO_EVENT_GO = 7          /* a Go function, called back through cgo */
};

typedef struct mybgo_event {
  int kind;      /* MYBGO_EVENT_* */
  int direction; /* +1, -1, 0 */
  size_t limit;  /* the most roots to report */
  size_t body;   /* PLANE, RANGE, RANGE_RATE, ALTITUDE (MYB_EVENT_ORIGIN for the first three) */
  double normal[3];
  double range;
  myb_ellipsoid ellipsoid; /* ALTITUDE, and the ellipsoid of the station */
  double altitude;
  double longitude; /* ELEVATION: the station on `body` */
  double latitude;
  double height;
  double elevation;
  size_t source; /* ECLIPSE */
  size_t occulter;
  int inner;
  uintptr_t handle; /* GO */
} mybgo_event;

typedef struct mybgo_found {
  size_t event; /* the index of the event in the request */
  double t;
  double y[6];
} mybgo_found;

/* The roots of every event from t0 to t_end (each event searched on an orbit of its own, all
 * the same orbit: the steps do not depend on what is asked), event by event and in the order
 * of the propagation within one; *found: an array of *found_count records allocated here (free
 * it with mybgo_free; NULL when none were found), *y_end the state at t_end, *steps the
 * accepted steps of the propagation. *failed: the index of the event whose search failed. */
myb_status mybgo_events(mybgo_model* model, const myb_orbit_options* options, double t0,
                        const double* y0, double t_end, const mybgo_event* events,
                        size_t event_count, mybgo_found** found, size_t* found_count, double* y_end,
                        int64_t* steps, size_t* failed);

void mybgo_free(void* pointer);

/* ---- Estimation ---------------------------------------------------------------------------- */

enum {
  MYBGO_MEASURE_POSITION = 1, /* 3 values */
  MYBGO_MEASURE_STATE = 2,    /* 6 values */
  MYBGO_MEASURE_RANGE = 3,    /* from a station, 1 value */
  MYBGO_MEASURE_RANGE_RATE = 4,
  MYBGO_MEASURE_RADEC = 5, /* 2 values, the first an angle */
  MYBGO_MEASURE_AZEL = 6,
  MYBGO_MEASURE_LIGHT = 16 /* added to the four of a station: with the light time */
};

typedef struct mybgo_station {
  size_t body;
  myb_ellipsoid ellipsoid;
  double longitude;
  double latitude;
  double height;
} mybgo_station;

typedef struct mybgo_light {
  double c;
  double gm;
  int two_way;
  int aberration;
} mybgo_light;

typedef struct mybgo_observation {
  double t;
  int kind;       /* MYBGO_MEASURE_* */
  size_t station; /* the index of the station for the measurements from one */
  double value[MYB_MEASURE_MAX];
  double sigma[MYB_MEASURE_MAX];
  int bias;
} mybgo_observation;

typedef struct mybgo_fit_options {
  int max_iterations;
  double tolerance;
  double epsilon;
  size_t bias_count;
} mybgo_fit_options;

/* myb_fit_orbit over the observations: y0 (in the first guess, out the estimate) and the
 * parameters of the model, in place; fixed and prior_sigma (6 + np + bias_count each, may be
 * NULL), bias (bias_count values, in and out), covariance and residuals as myb_fit_orbit.
 * light may be NULL when no observation takes the light time. *failed: the index of the
 * observation that a description error refers to (count otherwise). */
myb_status mybgo_fit(mybgo_model* model, double t0, double* y0, const mybgo_station* stations,
                     size_t station_count, const mybgo_light* light,
                     const mybgo_observation* observations, size_t count,
                     const mybgo_fit_options* options, const unsigned char* fixed,
                     const double* prior_sigma, double* bias, myb_fit_result* result,
                     double* covariance, double* residuals, size_t* failed);

/* ---- Files ------------------------------------------------------------------------------------
 */

/* One OEM message from segments laid end to end: segment j has counts[j] epochs and states
 * (taken in turn from `epochs` and `states`, 6 each) and the texts names[6 j .. 6 j + 5] =
 * object name, object id, centre name, reference frame, time system, interpolation (the last
 * may be NULL) with degrees[j]. *text: allocated here (mybgo_free), *length its characters. */
myb_status mybgo_oem_write(const char* originator, const myb_calendar* created, int notation,
                           size_t segment_count, const char* const* names, const int* degrees,
                           const size_t* counts, const myb_calendar* epochs, const double* states,
                           char** text, size_t* length);

/* The metadata keywords keys[0 .. n) of a segment (NULL where it has none; pointers into the
 * message) and its lines: epochs, states (6 each) and accelerations (3 each, may be NULL);
 * returns how many lines carry an acceleration. */
void mybgo_oem_meta(const myb_oem* oem, size_t segment, const char* const* keys, size_t n,
                    const char** values);
size_t mybgo_oem_lines(const myb_oem* oem, size_t segment, myb_calendar* epochs, double* states,
                       double* accelerations);

/* One SPK file from segments laid end to end (types 9 or 13): segment j of counts[j] epochs
 * and states. Written to path, or with path NULL into *data (allocated here, mybgo_free) of
 * *size bytes. *failed: the index of the segment the writer refused. */
myb_status mybgo_spk_write(const char* path, const char* internal_name, size_t segment_count,
                           const myb_spk_segment_header* headers, const int* types,
                           const int* degrees, const size_t* counts, const double* epochs,
                           const double* states, void** data, size_t* size, size_t* failed);

#include "glue_tasks.h"

#endif /* MYBGO_GLUE_H_ */
