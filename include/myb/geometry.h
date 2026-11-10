/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/geometry.h — the geometry of the surface of a body: where a ray meets it, the point
 * nearest to an observer, the normal there and the angles of illumination.
 *
 * Everything is in the body-fixed frame of the body, centre at the origin (the caller turns
 * inertial vectors with the orientation of the body: myb_env_orientation, myb/env.h). Two
 * shapes: the triaxial ellipsoid of semi-axes radii[3] along x, y, z, and the closed
 * triangular mesh of myb/grav.h (myb_polyhedron). The geometry is instantaneous: light time
 * and aberration are the caller's (myb/observe.h). */
#ifndef MYB_GEOMETRY_H_
#define MYB_GEOMETRY_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/env.h"
#include "myb/frames.h"
#include "myb/grav.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Triaxial ellipsoid --------------------------------------------------------------------- */

/* The first point where the ray origin + s direction, s >= 0, meets the ellipsoid: returns 1
 * and gives the distance s (m; the direction may have any length) and the point (either may
 * be NULL), or returns 0 when the ray misses it (also for a direction of zero length and for
 * input that is not finite). From inside the ellipsoid it is the point where the ray leaves.
 * The radii must be positive. */
MYB_API int myb_ray_ellipsoid(const double radii[3], const double origin[3],
                              const double direction[3], double* distance, double point[3]);

/* The point of the ellipsoid nearest to a position and the height of the position over it
 * along the normal there (m, negative inside; either output may be NULL): for an observer,
 * the point under it. Deep inside, in a plane of symmetry, two points (or a circle of them)
 * may be equally near: the one on the positive side of the shortest axis is given. Errors:
 * MYB_E_DOMAIN (input not finite), MYB_E_NOCONV. */
MYB_API MYB_NODISCARD myb_status myb_ellipsoid_nearest(const double radii[3],
                                                       const double position[3],
                                                       double point[3], double* height);

/* The outward unit normal at a point of the ellipsoid. */
MYB_API void myb_ellipsoid_normal(const double radii[3], const double point[3],
                                  double normal[3]);

/* The limb of the ellipsoid seen from a point outside it — where the lines of sight touch
 * the surface. It is an ellipse in a plane: centre + cos(u) major + sin(u) minor, major and
 * minor its semi-axes, perpendicular, |major| >= |minor|. Errors: MYB_E_DOMAIN (the observer
 * on or inside the ellipsoid, input not finite). */
MYB_API MYB_NODISCARD myb_status myb_ellipsoid_limb(const double radii[3],
                                                    const double observer[3], double centre[3],
                                                    double major[3], double minor[3]);

/* The terminator of the ellipsoid lit by a spherical source of the given radius at the
 * position given: `count` points (3 values each) in half-planes that hinge on the line from
 * the centre of the ellipsoid to that of the source, equally spaced round it, the first one
 * on the side of the +z axis of the frame (of +x if the source is on the z axis).
 * penumbral == 0: the edge of the region that sees no part of the source (the plane tangent
 * to the surface there touches the source from the inside of the shadow); != 0: the edge of
 * the region that sees all of it. For a source of radius 0 both are where the source is on
 * the horizon. Errors: MYB_E_DOMAIN (the source not outside the ellipsoid, a negative
 * radius, input not finite). */
MYB_API MYB_NODISCARD myb_status myb_ellipsoid_terminator(const double radii[3],
                                                          const double source[3],
                                                          double source_radius, int penumbral,
                                                          size_t count, double* points);

/* ---- Triangular mesh ------------------------------------------------------------------------ */

MYB_API size_t myb_polyhedron_face_count(const myb_polyhedron* polyhedron);

/* The outward unit normal of a face. */
MYB_API void myb_polyhedron_face_normal(const myb_polyhedron* polyhedron, size_t face,
                                        double normal[3]);

/* As myb_ray_ellipsoid for the mesh: the nearest face met by the ray, from outside or from
 * inside, and its index (may be NULL). The faces are kept in a tree of boxes, halved until
 * four are left in each: a ray examines the boxes it passes through — tens of faces of a
 * mesh of a hundred thousand. Faces are taken a part in 10^12 larger than they are, so that
 * a ray through an edge shared by two faces meets one of them (the first in the mesh). */
MYB_API int myb_polyhedron_ray(const myb_polyhedron* polyhedron, const double origin[3],
                               const double direction[3], double* distance, double point[3],
                               size_t* face);

/* The solid angle (sr) that the surface of the mesh subtends at r: 4 pi inside, 0 outside.
 * On the surface itself it is one of the two or anything between: the face under the point
 * counts for 2 pi with the sign of a rounding error. */
MYB_API double myb_polyhedron_solid_angle(const myb_polyhedron* polyhedron, const double r[3]);

/* ---- What an observer sees ------------------------------------------------------------------
 * The functions above take the body where it is. An observer at the model time t sees it
 * where it was, and turned as it was, when the light left it: a light time tau earlier, tau
 * the distance then over c, found by iteration on the states of the environment (myb/env.h).
 * y = (x, v) is the observer relative to the origin of the environment; its velocity enters
 * the aberration and should then be that relative to an inertial origin. */

/* The vector (m, inertial axes) from the observer to the body at t - tau — with
 * aberration != 0 turned by the observer's own motion (myb_aberration, myb/observe.h: the
 * exact transformation; the classical formula of other software differs from it by
 * (v / c)^2) —, the orientation of the body at t - tau and tau (s). Any output may be NULL.
 * Errors: MYB_E_DOMAIN, the statuses of the environment. */
MYB_API MYB_NODISCARD myb_status myb_view_body(myb_env* env, double t, size_t body,
                                               const double y[6], double c, int aberration,
                                               double position[3], myb_orientation* orientation,
                                               double* light_time);

/* Where the line of sight `direction` (inertial axes, any length) of the observer meets the
 * ellipsoid of the body: the point in the body-fixed frame and the light time from it, the
 * body taken where it was and turned as it was that long before t. *hit = 0 if the line of
 * sight misses the body (point and light_time, which may be NULL, are then untouched).
 * Errors as above. */
MYB_API MYB_NODISCARD myb_status myb_view_intercept(myb_env* env, double t, size_t body,
                                                    const double radii[3], const double y[6],
                                                    const double direction[3], double c,
                                                    int* hit, double point[3],
                                                    double* light_time);

/* ---- Illumination --------------------------------------------------------------------------- */

/* The angles at a point of a surface with the outward normal there (any length), all vectors
 * in one frame: phase — between the directions to the source and to the observer; incidence —
 * between the normal and the direction to the source (more than pi/2: the source is below the
 * local horizon); emission — between the normal and the direction to the observer. Radians in
 * [0, pi]; any output may be NULL. Errors: MYB_E_DOMAIN (a vector of zero length, the source
 * or the observer at the point, input not finite). */
MYB_API MYB_NODISCARD myb_status myb_illumination_angles(const double point[3],
                                                         const double normal[3],
                                                         const double source[3],
                                                         const double observer[3],
                                                         double* phase, double* incidence,
                                                         double* emission);

#ifdef __cplusplus
}
#endif

#endif /* MYB_GEOMETRY_H_ */
