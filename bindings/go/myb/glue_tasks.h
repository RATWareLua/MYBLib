/* SPDX-License-Identifier: BSD-3-Clause */
/* glue_tasks.h — the tasks of the C half of the Go binding beyond the force model: systems of
 * bodies, estimation, the restricted problem of three bodies, targeting and optimal transfers,
 * attitude, observables and files. As in glue.h: records of numbers and of pointers to objects
 * of the library, one call a task, the functions of the library called in the order a C
 * program calls them; the arrays are those of the caller, nothing kept after the call unless
 * it is said. Included by glue.h; declarations only. */
#ifndef MYBGO_GLUE_TASKS_H_
#define MYBGO_GLUE_TASKS_H_

/* ---- Systems of bodies (myb/system.h, myb/dyn.h) ----------------------------------------------
 */

typedef struct mybgo_system_body {
  double gm;
  double radius;                 /* 0: no sphere (myb_system_set_radius) */
  mybgo_orientation orientation; /* given from outside; NONE: none */
  myb_gravity_field* field;      /* NULL: no figure */
  int degree;
  int has_rotation; /* an integrated rotation */
  myb_system_rotation rotation;
  int has_tide;
  myb_system_tide tide;
  int has_mass_rate;
  double mass_rate;
  double mass_epoch;
} mybgo_system_body;

enum {
  MYBGO_LINK_INTERACTION = 1, /* a receives from b at `level` (myb_nbody_set_interaction) */
  MYBGO_LINK_FIELD = 2,       /* the figure of a on the point b */
  MYBGO_LINK_FIGURES = 3,     /* the figures of a and b */
  MYBGO_LINK_TIDE = 4,        /* the tide raised on a by b acts on c */
  MYBGO_LINK_DISTORTION = 5,  /* b distorts a */
  MYBGO_LINK_PN_PAIR = 6,     /* the terms beyond 1PN of a and b */
  MYBGO_LINK_DRAGGING = 7,    /* the spin of a drags b */
  MYBGO_LINK_CONTACT = 8,     /* a touches the ground b with points */
  MYBGO_LINK_ACTION = 9       /* a Go function acts on a */
};

typedef struct mybgo_system_link {
  int kind; /* MYBGO_LINK_* */
  size_t a;
  size_t b;
  size_t c;
  int level;
  myb_system_pn_pair pn;
  double spin[3];
  myb_contact_surface surface;
  myb_contact_law law;
  size_t first_point; /* CONTACT: points[first_point .. first_point + point_count) */
  size_t point_count;
  uintptr_t handle; /* ACTION */
} mybgo_system_link;

typedef struct mybgo_system_options {
  int has_relativity;
  myb_ppn relativity;
  int fmm_order; /* 0: exact sums */
  double fmm_theta;
} mybgo_system_options;

typedef struct mybgo_system mybgo_system;

/* A system of the bodies with their levels, the guided bodies (gm and provider each) and the
 * contact points (copied, kept by the system). The objects pointed to must outlive it.
 * Errors: those of the library; *stage (MYBGO_STAGE_BODY: a body, MYBGO_STAGE_TERM: a link,
 * MYBGO_STAGE_PARAMETER: a guided body) and *item name the record. */
myb_status mybgo_system_create(const mybgo_system_options* options, const mybgo_system_body* bodies,
                               size_t body_count, const mybgo_system_link* links, size_t link_count,
                               const myb_contact_point* points, size_t point_count,
                               const mybgo_provider* guided, const double* guided_gm,
                               size_t guided_count, mybgo_system** system, int* stage,
                               size_t* item);
void mybgo_system_destroy(mybgo_system* system); /* accepts NULL */
size_t mybgo_system_coordinates(const mybgo_system* system);
size_t mybgo_system_bodies(const mybgo_system* system);
int mybgo_system_present(const mybgo_system* system, size_t body);
double mybgo_system_gm(const mybgo_system* system, size_t body);

typedef struct mybgo_rotation_state {
  size_t body;
  double q[4];
  double omega[3];
  double omega_core[3];
  int has_core;
} mybgo_rotation_state;

typedef struct mybgo_spin_state {
  size_t pair;
  double chi1[3];
  double chi2[3];
} mybgo_spin_state;

/* A state y (2 n values) from the states of the bodies (6 each), the rotations and the spins,
 * and back (the records name the bodies and pairs to read). */
myb_status mybgo_system_pack(const mybgo_system* system, const double* bodies,
                             const mybgo_rotation_state* rotations, size_t rotation_count,
                             const mybgo_spin_state* spins, size_t spin_count, double* y);
myb_status mybgo_system_unpack(const mybgo_system* system, const double* y, double* bodies,
                               mybgo_rotation_state* rotations, size_t rotation_count,
                               mybgo_spin_state* spins, size_t spin_count);

/* The deflections of the springs of a friction that sticks (myb_system_pack_contact, with pack
 * != 0, or myb_system_unpack_contact): the contact number contacts[k] among those in place
 * takes 3 values a point, in turn from `deflections` (capacity values). Packing, counts[k] is
 * the points given for it and must be its count; unpacking, counts[k] receives that count.
 * MYB_E_DOMAIN for no such contact or another count. *failed: the record that failed. */
myb_status mybgo_system_contacts(const mybgo_system* system, int pack, double* y,
                                 const size_t* contacts, size_t* counts, size_t count,
                                 double* deflections, size_t capacity, size_t* failed);

/* ---- The force of a contact ------------------------------------------------------------------ */

typedef struct mybgo_contact_sample {
  double x[3];
  double v[3];
  double attitude[9];
  double omega[3];
  int has_attitude;
  int has_omega;
} mybgo_contact_sample;

/* The force of the surface on the points at each of `count` states (myb_contact_force_stick):
 * with a friction that sticks, deflections holds 3 values a point a state and
 * deflection_rates (may be NULL) as many; with the smooth one both are ignored (may be NULL).
 * *failed: the state that failed. */
myb_status mybgo_contact_forces(const myb_contact_surface* surface, const myb_contact_law* law,
                                const myb_contact_point* points, size_t point_count,
                                const mybgo_contact_sample* samples, size_t count,
                                const double* deflections, double* deflection_rates,
                                myb_contact_result* results, size_t* failed);

enum {
  MYBGO_CONTACT_STOP = 0,    /* the run ends at the contact */
  MYBGO_CONTACT_MERGE = 1,   /* the lighter body is merged into the heavier (b into a if equal) */
  MYBGO_CONTACT_BOUNCE = 2,  /* the two rebound with the restitution */
  MYBGO_CONTACT_CONTINUE = 3 /* nothing is done */
};

typedef struct mybgo_system_event {
  int kind; /* myb_system_event_kind */
  size_t a;
  size_t b;
  double t;
  double distance;
  int action;      /* MYBGO_CONTACT_*: what was done (CONTINUE for an approach) */
  size_t survivor; /* MERGE: the body that remains */
  myb_system_merger merger;
} mybgo_system_event;

typedef struct mybgo_system_run_options {
  myb_ias15_options ias15; /* precision, step limits, budget (myb_system_integrator) */
  double watch;            /* myb_system_advance */
  int policy;              /* MYBGO_CONTACT_* */
  double restitution;      /* BOUNCE */
  size_t limit;            /* the most events */
} mybgo_system_run_options;

/* From y0 at t0 through the times (in the order of the integration): the integrator ends a step
 * at every time asked (myb_system_advance to each), and at the events: contacts acted on by the
 * policy, approaches noted. states: 2 n each; *events: allocated here (mybgo_free), NULL when
 * none; *t_end and y_end where the run stopped; *failed the index of the time not reached. */
myb_status mybgo_system_run(mybgo_system* system, const mybgo_system_run_options* options,
                            double t0, const double* y0, const double* times, size_t count,
                            double* states, mybgo_system_event** events, size_t* event_count,
                            double* t_end, double* y_end, myb_ias15_stats* stats, size_t* failed);

/* myb_system_acceleration at count states (a: n each); the invariants of the point masses
 * (myb_nbody_invariants_of, myb_nbody_pn_invariants_of) of count states. */
myb_status mybgo_system_accelerations(mybgo_system* system, const double* t, const double* y,
                                      size_t count, double* a, size_t* failed);
myb_status mybgo_system_invariants(mybgo_system* system, const double* y, size_t count,
                                   myb_nbody_invariants* invariants, size_t* failed);
myb_status mybgo_system_pn_invariants(mybgo_system* system, const double* y, size_t count,
                                      myb_nbody_pn_invariants* invariants, size_t* failed);
myb_status mybgo_system_partials(mybgo_system* system, const myb_ias15_options* options, double t0,
                                 const double* y0, double t1, const size_t* bodies, size_t count,
                                 double* y1, double* phi, double* sensitivity);

/* A Wisdom-Holman integration of point masses (myb_wh): massive + test bodies, x and v 3 each,
 * the states at the times (myb_wh_integrate, then myb_wh_dense). */
myb_status mybgo_wh_run(const double* gm, size_t massive, size_t test,
                        const myb_wh_options* options, double t0, const double* x, const double* v,
                        const double* times, size_t count, double* xs, double* vs,
                        myb_wh_stats* stats, size_t* failed);

/* ---- Estimation (myb/estim.h) ---------------------------------------------------------------- */

typedef struct mybgo_filter_options {
  double epsilon;
  double process_noise;
  size_t bias_count;
  int smooth; /* keep every instant and smooth at the end */
} mybgo_filter_options;

/* The filter from (t0, y0, covariance of the 6 + np unknowns) towards t_end over the
 * observations, then the smoother if asked: innovations (the values of each observation in
 * turn), the instants kept (*steps of them: their times, unknowns u = 6 + np + biases each and
 * covariances u^2 each — room for count + 1), the final state, covariance and biases.
 * *failed: the observation that failed. */
myb_status mybgo_filter_run(mybgo_model* model, const mybgo_filter_options* options, double t0,
                            const double* y0, const double* covariance, double t_end,
                            const double* bias, const mybgo_station* stations, size_t station_count,
                            const mybgo_light* light, const mybgo_observation* observations,
                            size_t count, double* innovations, size_t* steps, double* step_t,
                            double* step_unknowns, double* step_covariance, double* final_t,
                            double* final_y, double* final_covariance, double* final_bias,
                            size_t* failed);

/* A measurement (kind MYBGO_MEASURE_*, from `station` when it needs one) at count instants and
 * states: values (as many as the kind has, each) and partials (6 per value, may be NULL). */
myb_status mybgo_measure(mybgo_model* model, int kind, const mybgo_station* station,
                         const mybgo_light* light, const double* t, const double* y, size_t count,
                         double* values, double* partials, size_t* failed);

/* ---- The restricted problem of three bodies (myb/cr3bp.h) ------------------------------------ */

/* The family through an orbit: the orbit corrected, then count - 1 steps along the family
 * (myb_cr3bp_orbit_continue); members: count orbits; *made how many were found. */
myb_status mybgo_cr3bp_family(double mu, const myb_cr3bp_orbit* orbit, double step, size_t count,
                              double tolerance, int max_iterations, myb_cr3bp_orbit* members,
                              size_t* made);
/* The same for tori: a step of 0 in the family, then count - 1 steps. */
myb_status mybgo_cr3bp_torus_family(double mu, const myb_cr3bp_torus* torus,
                                    myb_cr3bp_torus_family family, double step, size_t count,
                                    double tolerance, int max_iterations, myb_cr3bp_torus* members,
                                    size_t* made);
/* Over slices: the Jacobi constants; states on a manifold at phases; points of a torus at
 * angles; states on the manifold of a torus at angles; propagations of count states (to a
 * section when it is given). */
void mybgo_cr3bp_jacobi(double mu, const double* y, size_t count, double* jacobi);
myb_status mybgo_cr3bp_manifold(double mu, const double* y0, double period, const double* phase,
                                size_t count, int unstable, int sign, double epsilon, double* y,
                                size_t* failed);
void mybgo_cr3bp_torus_states(const myb_cr3bp_torus* torus, const double* angle, size_t count,
                              double* y);
myb_status mybgo_cr3bp_torus_manifold(const myb_cr3bp_torus* torus,
                                      const myb_cr3bp_torus_spectrum* spectrum, const double* angle,
                                      size_t count, int unstable, int sign, double epsilon,
                                      double* y, size_t* failed);
myb_status mybgo_cr3bp_propagate(double mu, const double* y0, size_t count, double duration,
                                 const myb_cr3bp_section* section, double* y, double* elapsed,
                                 int* crossed, size_t* failed);
myb_status mybgo_cr3bp_frames(int to_ephemeris, const double* body1, const double* body2,
                              double gm1, double gm2, const double* in, size_t count, double* out,
                              size_t* failed);

/* ---- Targeting and transfers (myb/lambert.h, myb/transfer.h, myb/relative.h) --------------- */

/* Lambert's problem for count pairs (r1, r2, tof; 3, 3, 1 each): up to 1 + 2 max_revolutions
 * solutions each (`capacity` a problem), counts[k] of them. */
myb_status mybgo_lambert_many(double mu, const double* r1, const double* r2, const double* tof,
                              size_t count, int retrograde, int max_revolutions, size_t capacity,
                              myb_lambert_solution* solutions, size_t* counts, size_t* failed);
/* The grid of departures and arrivals between two providers: the states of the departure body
 * at the departures and of the arrival body at the arrivals, and for each pair (i, j) with
 * arrival after departure the solutions of Lambert between the two positions over the time
 * between (cell i * arrivals + j, `capacity` solutions a cell, counts 0 where the arrival is
 * not after the departure). */
myb_status mybgo_lambert_grid(double mu, const mybgo_provider* departure,
                              const mybgo_provider* arrival, const myb_time* departures,
                              size_t departure_count, const myb_time* arrivals,
                              size_t arrival_count, int retrograde, int max_revolutions,
                              size_t capacity, double* departure_states, double* arrival_states,
                              myb_lambert_solution* solutions, size_t* counts, size_t* failed);
myb_status mybgo_bplanes(double gm, const double* states, size_t count, const double pole[3],
                         myb_bplane* bplanes, size_t* failed);
/* myb_shoot with a Go residual function (handle): u updated in place. */
myb_status mybgo_shoot(uintptr_t handle, size_t n, size_t m, double* u,
                       const myb_shoot_options* options, int* iterations, double* residual);

enum {
  MYBGO_RELATIVE_FROM_INERTIAL = 1, /* in: chief, deputy (6 each); out: relative (6) */
  MYBGO_RELATIVE_TO_INERTIAL = 2,   /* in: chief, relative; out: deputy */
  MYBGO_RELATIVE_ROE = 3            /* in: chief, deputy; out: roe (6) */
};

myb_status mybgo_relative(int kind, double mu, const double* chief, const double* other,
                          size_t count, double* out, size_t* failed);
/* The transition matrices for count offsets: of a chief on any conic (myb_relative_stm), or
 * of the circle of mean motion n (myb_cw_stm) when chief is NULL. */
myb_status mybgo_relative_stm(double mu, const double* chief, double n, const double* dt,
                              size_t count, double* phi, size_t* failed);
void mybgo_roe_positions(double a, const double* u, size_t count, const double roe[6],
                         double* positions);

/* ---- Conjunctions and the uncertainty (myb/conjunction.h, myb/uncert.h) --------------------- */

/* The closest approaches of two providers: *approaches allocated here (mybgo_free). */
myb_status mybgo_approach_find(const mybgo_provider* first, const mybgo_provider* second, double t0,
                               double t1, double step, myb_approach** approaches, size_t* count);

/* An orbit of a model from (t0, y0) towards t_end, made with the variational equations. */
typedef struct mybgo_orbit_start {
  myb_orbit_options options;
  double t0;
  double y0[6];
  double t_end;
  double covariance[36];
} mybgo_orbit_start;

/* myb_conjunction_find of the two orbits (both started at one time) up to t1; *found allocated
 * here (mybgo_free). The long encounter (myb_conjunction_long) when probability is not NULL:
 * then *found is not used. */
myb_status mybgo_conjunctions(mybgo_model* first, const mybgo_orbit_start* first_start,
                              mybgo_model* second, const mybgo_orbit_start* second_start,
                              double radius, double t1, double step, myb_conjunction** found,
                              size_t* count, double* probability, double* error);

/* The covariance carried along an orbit: P(t) = Phi P0 Phi^T at the times (states 6 each,
 * covariances 36 each). */
myb_status mybgo_covariance_propagate(mybgo_model* model, const mybgo_orbit_start* start,
                                      const double* times, size_t count, double* states,
                                      double* covariances, size_t* failed);
/* myb_covariance_rtn of count states and covariances. */
myb_status mybgo_covariance_rtn(const double* states, const double* covariances, size_t count,
                                double* out, size_t* failed);
/* myb_collision_probability_long with a Go function giving the relative Gaussian. */
myb_status mybgo_collision_long(uintptr_t handle, double t0, double t1, double step, double radius,
                                double* probability, double* error);

/* ---- Low thrust and optimal transfers (myb/thrust.h, myb/optimal.h) ------------------------- */

/* The direction of the Q-law at count sets of elements (direction 3, q, rate, effectivity 1
 * each). */
myb_status mybgo_qlaw_directions(const myb_qlaw* law, const myb_elements* elements, size_t count,
                                 double* direction, double* q, double* q_rate, double* effectivity,
                                 size_t* failed);
/* myb_qlaw_fly: the law, the engine (its law pointer set here) and the field (NULL: the point
 * mass of the law). */
myb_status mybgo_qlaw_fly(const myb_qlaw* law, const myb_qlaw_term* engine, myb_force_model* field,
                          double rtol, double atol, double t0, const double* y0, double t_limit,
                          double step, myb_qlaw_arrival* arrival);

/* The chain with its bodies given by providers (bodies[k].state and context set here); the
 * functions of myb/optimal.h. */
myb_status mybgo_chain_size(const myb_chain* chain, const mybgo_provider* providers,
                            size_t* unknowns, size_t* equalities, size_t* inequalities);
myb_status mybgo_chain_guess(const myb_chain* chain, const mybgo_provider* providers, double* x);
myb_status mybgo_chain_evaluate(const myb_chain* chain, const mybgo_provider* providers,
                                const double* x, double* propellant, double* constraints,
                                double* gradient, double* jacobian);
myb_status mybgo_chain_solve(const myb_chain* chain, const mybgo_provider* providers,
                             const myb_nlp_options* options, double* x, myb_chain_result* result);

/* ---- Attitude (myb/surface.h, myb/attitude.h) -------------------------------------------------
 */

enum {
  MYBGO_ATTITUDE_FIXED = 1,     /* myb_attitude_fixed: m */
  MYBGO_ATTITUDE_ORBITAL = 2,   /* myb_attitude_orbital: body */
  MYBGO_ATTITUDE_SUN_NADIR = 3, /* myb_attitude_sun_nadir: body, source */
  MYBGO_ATTITUDE_SPIN = 4       /* myb_attitude_spin: m (M0), rate, phase */
};

typedef struct mybgo_attitude {
  int kind;
  double m[9];
  size_t body;
  size_t source;
  double rate;
  double phase;
} mybgo_attitude;

/* The law at count instants and states: m (9 each). model: the environment of the laws that
 * need one (may be NULL for the others). */
myb_status mybgo_attitudes(mybgo_model* model, const mybgo_attitude* law, const double* t,
                           const double* y, size_t count, double* m, size_t* failed);
/* The free rotation (myb_free_rotation_init, then _state) at count times: q 4, l 3 each. */
myb_status mybgo_free_rotation(const double moments[3], const double q0[4], const double l0[3],
                               const double* t, size_t count, double* q, double* l);

/* ---- Observables (myb/observe.h) ------------------------------------------------------------ */

typedef struct mybgo_light_model {
  double c;
  double gamma;
} mybgo_light_model;

/* myb_light_time of the target for count instants and observers (6 each); the bodies that
 * delay the light by their gm and providers. */
myb_status mybgo_light_times(const mybgo_light_model* model, const double* gm,
                             const mybgo_provider* bodies, size_t body_count,
                             const mybgo_provider* target, const myb_time* t,
                             const double* observers, size_t count, int direction,
                             myb_light_path* paths, size_t* failed);
myb_status mybgo_round_trips(const mybgo_light_model* model, const double* gm,
                             const mybgo_provider* bodies, size_t body_count,
                             const mybgo_provider* station, const mybgo_provider* target,
                             const myb_time* t, size_t count, double* total, myb_time* bounce,
                             myb_time* transmit, size_t* failed);
void mybgo_aberrations(const double* p, const double* v, size_t count, double c, double* out);
void mybgo_deflections(double gm, double c, double gamma, const double* p, const double* q,
                       const double* e, const double* em, size_t count, double dlim, double* out);
void mybgo_shapiro_delays(double gm, double c, double gamma, const double* x1, const double* x2,
                          const double* xb, size_t count, double* delay);
void mybgo_radecs(const double* p, size_t count, double* ra, double* dec);
myb_status mybgo_plasma_columns(const myb_plasma_profile* profile, const double* p1,
                                const double* p2, size_t count, double* column, size_t* failed);

/* ---- Files ---------------------------------------------------------------------------------- */

/* An OPM or OMM as text (myb_opm_write_as, myb_omm_write_as): *text allocated (mybgo_free). */
myb_status mybgo_opm_write(const myb_opm* opm, int notation, char** text, size_t* length);
myb_status mybgo_omm_write(const myb_omm* omm, int notation, char** text, size_t* length);

/* An SP3 file from its parts (myb_sp3_data; the optional arrays may be NULL). */
myb_status mybgo_sp3_write(const myb_sp3_data* head, const char* const* satellites,
                           const myb_calendar* epochs, const double* positions,
                           const double* velocities, const double* clocks,
                           const double* clock_rates, const double* accuracy, const double* sigmas,
                           const unsigned char* flags, char** text, size_t* length);
/* Every record of a loaded SP3: positions and velocities (3 each, NaN where missing), what
 * each gives (1 position, 2 velocity), clocks and rates, sigmas (8 each), flags. Every array
 * epochs x satellites. */
void mybgo_sp3_records(const myb_sp3* sp3, double* positions, double* velocities, int* given,
                       double* clocks, double* clock_rates, double* sigmas, unsigned char* flags);

enum {
  MYBGO_CK_POINTINGS = 1, /* types 1 and 3: ticks, quaternions, rates (or none), intervals */
  MYBGO_CK_INTERVALS = 2  /* type 2: starts, stops, quaternions, rates, seconds per tick */
};

typedef struct mybgo_ck_segment {
  int kind;
  int type; /* POINTINGS: 1 or 3 */
  myb_ck_segment_header header;
  size_t count;          /* pointings or intervals */
  int has_rates;         /* POINTINGS */
  size_t interval_count; /* POINTINGS, type 3: the starts of intervals; 0: one interval */
} mybgo_ck_segment;

/* One CK file from segments laid end to end: each takes count values from ticks (POINTINGS)
 * or from starts and stops (INTERVALS), quaternions 4 and rates 3 a value, seconds_per_tick
 * one an interval, interval starts interval_count. To path, or into *data (mybgo_free).
 * *failed: the segment refused. */
myb_status mybgo_ck_write(const char* path, const char* internal_name,
                          const mybgo_ck_segment* segments, size_t segment_count,
                          const double* ticks, const double* stops, const double* quaternions,
                          const double* rates, const double* seconds_per_tick,
                          const double* interval_starts, void** data, size_t* size, size_t* failed);
/* myb_ck_pointing_at for count ticks. */
myb_status mybgo_ck_pointings(const myb_ck* ck, int32_t instrument, const double* ticks,
                              size_t count, double tolerance, int need_rate,
                              myb_ck_pointing* pointings, size_t* failed);
/* The clock: ticks -> instants and back, for count values. */
myb_status mybgo_sclk_times(const myb_sclk* sclk, const double* ticks, size_t count, myb_time* t,
                            size_t* failed);
myb_status mybgo_sclk_ticks(const myb_sclk* sclk, const myb_time* t, size_t count, double* ticks,
                            size_t* failed);

#endif /* MYBGO_GLUE_TASKS_H_ */
