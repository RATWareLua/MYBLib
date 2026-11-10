/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/events.h — ready-made event functions for myb_orbit_event (myb/orbit.h): functions
 * g(t, y) of the model time and the state y = (x, v) of the object whose sign changes mark
 * what is looked for — an apsis, a close approach, the crossing of a plane, of a distance or
 * of an altitude, the rise and set over a station, the edges of a shadow, the surface of a
 * shape. Each takes a context structure that the caller fills.
 *
 * Bodies are those of an environment (myb/env.h), given by index; MYB_EVENT_ORIGIN stands for
 * the origin of the environment itself (the context's env may then be NULL). The geometry is
 * instantaneous: no light time. A function has no status to return: when a provider of the
 * environment fails it gives NaN, which is no sign change, and the event is not found. */
#ifndef MYB_EVENTS_H_
#define MYB_EVENTS_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/env.h"
#include "myb/estim.h"
#include "myb/frames.h"
#include "myb/grav.h"
#include "myb/model.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MYB_EVENT_ORIGIN ((size_t)-1)

/* ---- Distance ------------------------------------------------------------------------------- */

typedef struct myb_range_event {
  myb_env* env;
  size_t body;
  double range; /* m; not used by myb_event_range_rate */
} myb_range_event;

/* |x - x_body| - range: the entry into a sphere around the body (falling) and the exit. */
MYB_API double myb_event_range(void* event, double t, const double* y);

/* The rate of the distance to the body (m/s): zero and growing at the closest approach — the
 * pericentre of an orbit about the body —, zero and falling at the farthest point. */
MYB_API double myb_event_range_rate(void* event, double t, const double* y);

/* ---- A plane -------------------------------------------------------------------------------- */

typedef struct myb_plane_event {
  myb_env* env;
  size_t body;
  double normal[3]; /* in the inertial axes; any length */
} myb_plane_event;

/* normal . (x - x_body): the ascending node on the plane through the body (growing) and the
 * descending one. With the pole of the body as the normal, its equator. */
MYB_API double myb_event_plane(void* event, double t, const double* y);

/* ---- Altitude ------------------------------------------------------------------------------- */

typedef struct myb_altitude_event {
  myb_env* env;
  size_t body; /* must have an orientation */
  myb_ellipsoid ellipsoid;
  double altitude; /* m */
} myb_altitude_event;

/* The geodetic height over the ellipsoid of the body minus altitude: the entry into an
 * atmosphere, the impact (altitude 0). */
MYB_API double myb_event_altitude(void* event, double t, const double* y);

/* ---- A station ------------------------------------------------------------------------------ */

typedef struct myb_elevation_event {
  const myb_station* station;
  double elevation; /* rad */
} myb_elevation_event;

/* The elevation of the object over the horizon of the station minus elevation (rad): rise
 * (growing) and set. Geometric: no refraction. */
MYB_API double myb_event_elevation(void* event, double t, const double* y);

/* ---- Shadows -------------------------------------------------------------------------------- */

typedef struct myb_eclipse_event {
  myb_env* env;
  size_t source;   /* the luminous body; with its radius */
  size_t occulter; /* with its radius */
  int inner;       /* 0: the disks touch from outside; != 0: from inside */
} myb_eclipse_event;

/* myb_shadow_margin (myb/force.h) of the two bodies seen from the object (rad). With
 * inner == 0 it falls through zero where the object enters the penumbra and grows where it
 * leaves; with inner != 0 the same for the umbra (or for the annular phase, when the occulter
 * looks smaller than the source). These are the instants where the light on the object, as
 * myb_term_radiation takes it, stops being smooth. */
MYB_API double myb_event_eclipse(void* event, double t, const double* y);

/* Declares the edges of the shadows of a radiation term to a model as breaks
 * (myb_force_model_add_break): two per occulter. events: 2 * term->occulter_count structures
 * that this call fills and that must outlive the model. Errors: MYB_E_OVERFLOW. */
MYB_API MYB_NODISCARD myb_status myb_term_radiation_breaks(const myb_radiation_term* term,
                                                           myb_force_model* model,
                                                           myb_eclipse_event* events);

/* ---- The surface of a shape ----------------------------------------------------------------- */

typedef struct myb_shape_event {
  myb_env* env;
  size_t body; /* must have an orientation */
  const myb_polyhedron* shape; /* in the body-fixed frame */
} myb_shape_event;

/* Positive outside the shape and negative inside (2 pi minus the solid angle the surface
 * subtends, myb_polyhedron_solid_angle): its fall is the impact. The function jumps at the
 * surface, so the crossing is bracketed and halved rather than interpolated — some fifty
 * evaluations, each a pass over the faces. */
MYB_API double myb_event_shape(void* event, double t, const double* y);

#ifdef __cplusplus
}
#endif

#endif /* MYB_EVENTS_H_ */
