/* SPDX-License-Identifier: BSD-3-Clause */
/* What the Jacobians of the surface terms (force/surface.c) take from the mesh
 * (grav/polyhedron.c). */
#ifndef MYB_SRC_FORCE_SURFACE_INTERNAL_H_
#define MYB_SRC_FORCE_SURFACE_INTERNAL_H_

#include "myb/surface.h"
#include "num/dual.h"

/* As myb_force_radiation_mesh, with gradient[3 i + j] = d force_i / d u_j (u taken as three
 * free numbers) for the parts of the mesh that are lit at u. */
[[nodiscard]] myb_status myb_internal_radiation_mesh_gradient(const myb_polyhedron* shape,
                                                              const myb_plate* surface,
                                                              double pressure,
                                                              const double u[3], int shadows,
                                                              double force[3],
                                                              double gradient[9]);

/* The attitude m (rows: inertial axes -> body frame) with its derivatives in the state. The
 * laws of force/surface.c are differentiated as they are written; any other law by central
 * differences of the law. Errors: those of the law, MYB_E_DOMAIN. */
[[nodiscard]] myb_status myb_internal_attitude_dual(myb_attitude_function attitude,
                                                    void* context, double t, const double x[3],
                                                    const double v[3], myb_dual m[9]);

#endif /* MYB_SRC_FORCE_SURFACE_INTERNAL_H_ */
