/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/surface.h — the object as more than a sphere: a surface of plates or a mesh, turned by
 * an attitude law, under light and under a rarefied gas; and the recoil of what it emits.
 *
 * Light on plates is myb_force_radiation_plates (myb/force.h); here are the same for a mesh
 * with the shadows its parts cast on each other, the force and the lift of a free-molecular
 * flow (Sentman), the attitude laws and the terms of a force model (myb/model.h) that put
 * them together with an environment (myb/env.h). */
#ifndef MYB_SURFACE_H_
#define MYB_SURFACE_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/env.h"
#include "myb/force.h"
#include "myb/grav.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Light on a mesh -------------------------------------------------------------------------
 * Every face of a closed triangular mesh (myb/grav.h; several closed parts may make one mesh:
 * a body and its panels) takes the light as a plate with the properties of `surface`
 * (specular, diffuse, thermal; its area and normal are not used). pressure: flux / c
 * (N m^-2); u: the unit direction of propagation of the light in the frame of the mesh.
 * shadows == 0: every face turned to the source is lit. shadows = k >= 1: each such face is
 * divided into k^2 equal triangles, and a part counts only if nothing of the mesh stands
 * between its centre and the source (a ray per part, each a pass over the faces: the cost
 * grows as faces^2 k^2). force (N) and torque (N m, about the origin of the mesh) are written;
 * either may be NULL. With specular = diffuse = 0 and thermal != 0 this is the zero-conductivity
 * model of D. P. Rubincam (Icarus 148, 2, 2000) whose torque, averaged over the rotation and
 * the orbit, is the YORP effect. Errors: MYB_E_DOMAIN (u not unit to 1e-12, fractions outside
 * [0, 1], values not finite). */
MYB_API MYB_NODISCARD myb_status myb_force_radiation_mesh(const myb_polyhedron* shape,
                                                          const myb_plate* surface,
                                                          double pressure, const double u[3],
                                                          int shadows, double force[3],
                                                          double torque[3]);

/* The YORP torque of a body that turns about the z axis of its mesh, averaged over a
 * rotation, with the heat conducted into the surface: every face re-emits what it absorbed
 * as a Lambertian radiator, but later and smoother — the answer of a half-space to a periodic
 * heating, 1 / (1 + theta sqrt(k) (1 + i) / 2) for the harmonic k of the rotation, theta the
 * thermal parameter at the rotation frequency (0: at once, Rubincam's model). The Sun stands
 * at the colatitude given from the spin axis; pressure, shadows as above; `phases` equal
 * steps of the rotation (at least 8). torque (N m): [0] in the plane of the axis and the
 * Sun, across the axis, on the Sun's side; [1] completing; [2] along the axis. The mean torque
 * along the axis needs only the mean of what each face emits, which conduction does not
 * change: [2] does not depend on theta. The other two — they turn the axis — do (D. Capek,
 * D. Vokrouhlicky, Icarus 172, 526, 2004). The momentum of the light absorbed is left out:
 * over a rotation and an orbit it averages away. Errors: MYB_E_DOMAIN. */
MYB_API MYB_NODISCARD myb_status myb_force_yorp(const myb_polyhedron* shape, double pressure,
                                                double sun_colatitude, double theta,
                                                int shadows, int phases, double torque[3]);

/* ---- Free-molecular flow ---------------------------------------------------------------------
 * The force of a gas so thin that its molecules meet the surface and not each other, re-emitted
 * diffusely (L. H. Sentman, Lockheed LMSC-448514, 1961; in the form of E. Doornbos,
 * Thermospheric Density and Wind Determination from Satellite Dynamics, 2012, eqs. 3.50-3.53):
 * drag along the flow and lift across it, for every plate, the faces turned away included (the
 * thermal motion of the molecules reaches them). */

typedef struct myb_gas {
  double temperature;      /* of the gas, K */
  double molar_mass;       /* of the gas, kg/mol */
  double wall_temperature; /* of the surface, K */
  double accommodation;    /* of the energy: 0 specular-like .. 1 full */
} myb_gas;

/* a += (1/2) rho V^2 / mass * sum A_i (C_D,i d + C_L,i l_i) over the plates (area and normal
 * used), v the velocity of the object relative to the gas in the frame of the plates. Errors:
 * MYB_E_DOMAIN (a gas, a mass, a density or a plate that is not valid; v = 0). */
MYB_API MYB_NODISCARD myb_status myb_force_gas_plates(double density, const myb_gas* gas,
                                                      double mass, const myb_plate* plates,
                                                      size_t count, const double v[3],
                                                      double a[3]);

/* The drag coefficient of a sphere in the same flow (the closed form), referred to its
 * cross-section. Errors: MYB_E_DOMAIN. */
MYB_API MYB_NODISCARD myb_status myb_gas_sphere_cd(const myb_gas* gas, double speed, double* cd);

/* Between the free molecules and the continuum. The Knudsen number of a body of the given
 * length (m) in a gas of the given density: the mean free path M / (sqrt(2) pi d^2 rho N_A)
 * of molecules of diameter d (m) over that length; and the weight of the free-molecular value
 * at that number in the bridge used for entry vehicles (R. G. Wilmoth, R. A. Mitcheltree,
 * J. N. Moss, J. Spacecraft Rockets 36, 436, 1999),
 *   C = C_continuum + (C_free - C_continuum) sin^2(pi (3 + log10 Kn) / 8),
 * 0 below Kn = 1e-3 and 1 above 10. The continuum value is the caller's: the library has no
 * continuum aerodynamics. */
MYB_API double myb_gas_knudsen(double density, double molar_mass, double diameter,
                               double length);
MYB_API double myb_gas_bridge(double knudsen);

/* ---- Attitude --------------------------------------------------------------------------------
 * The attitude of the object: m (rows) turns the inertial axes of the model into the body
 * frame of the object, body = M inertial, at the model time t and the state (x, v) relative
 * to the origin of the environment. */
typedef myb_status (*myb_attitude_function)(void* context, double t, const double x[3],
                                            const double v[3], double m[9]);

/* context: const double[9], a constant matrix. */
MYB_API MYB_NODISCARD myb_status myb_attitude_fixed(void* matrix, double t, const double x[3],
                                                    const double v[3], double m[9]);

/* Body axes along the orbit about a body: x from the body to the object, z along the orbital
 * angular momentum, y completing (along the motion for a circular orbit). */
typedef struct myb_orbital_attitude {
  myb_env* env;
  size_t body;
} myb_orbital_attitude;

MYB_API MYB_NODISCARD myb_status myb_attitude_orbital(void* context, double t, const double x[3],
                                                      const double v[3], double m[9]);

/* The yaw steering of a navigation satellite: z towards the body, y across the plane of the
 * body, the object and the source (the axis of its panels), x in that plane on the side of
 * the source. Undefined (MYB_E_DOMAIN) when the three are in line. */
typedef struct myb_sun_attitude {
  myb_env* env;
  size_t body;
  size_t source;
} myb_sun_attitude;

MYB_API MYB_NODISCARD myb_status myb_attitude_sun_nadir(void* context, double t,
                                                        const double x[3], const double v[3],
                                                        double m[9]);

/* A uniform rotation: body = R3(phase + rate t) M0 inertial, with M0 (rows) the attitude at
 * phase + rate t = 0; the body z axis is the axis of rotation. */
typedef struct myb_spin_attitude {
  double m0[9];
  double rate;  /* rad/s */
  double phase; /* rad */
} myb_spin_attitude;

MYB_API MYB_NODISCARD myb_status myb_attitude_spin(void* context, double t, const double x[3],
                                                   const double v[3], double m[9]);

/* ---- Terms of a force model ------------------------------------------------------------------ */

/* The light of a source on the surface of the object: plates (myb_force_radiation_plates) or,
 * if shape is given, a mesh (myb_force_radiation_mesh) with the properties of *surface on
 * every face; the shadows of the occulters as in myb_radiation_term. The object gives the
 * mass. */
typedef struct myb_surface_radiation_term {
  myb_env* env;
  size_t source;
  const myb_object* object;
  double c;
  size_t occulters[MYB_ENV_MAX_OCCULTERS];
  size_t occulter_count;
  myb_attitude_function attitude;
  void* attitude_context;
  const myb_plate* plates;
  size_t plate_count;
  const myb_polyhedron* shape;
  const myb_plate* surface;
  int shadows;
  const myb_refraction* atmospheres; /* as in myb_radiation_term */
} myb_surface_radiation_term;

MYB_API MYB_NODISCARD myb_status myb_term_surface_radiation(void* term, double t,
                                                            const double x[3], const double v[3],
                                                            double a[3]);
/* The Jacobians of the three terms of this header are analytic: the force in the body of the
 * object is differentiated as it is written, and so is the attitude when its law is one of
 * the four of this header (it is known by its address). The derivatives of an attitude law
 * of the caller's own are central differences of that law over 6e-6 of |x| and of |v| —
 * 12 calls of the law, not of the term. Where the force is not smooth the derivative is
 * that of the side the object is on: a plate that the light or the flow just grazes, the
 * parts of a mesh that are lit now. In a penumbra the light changes over tens of kilometres
 * and the Jacobian of this term is central differences of the term, as that of
 * myb_term_radiation. */
MYB_API MYB_NODISCARD myb_status myb_term_surface_radiation_jacobian(void* term, double t,
                                                                     const double x[3],
                                                                     const double v[3],
                                                                     double dadx[9],
                                                                     double dadv[9]);

/* The free-molecular flow of the atmosphere of a body on the plates of the object
 * (myb_force_gas_plates): density and wind as in myb_drag_term. */
typedef struct myb_surface_drag_term {
  myb_env* env;
  size_t body;
  myb_density_function density;
  void* density_context;
  myb_wind_function wind; /* may be NULL */
  void* wind_context;
  const myb_object* object;
  myb_gas gas;
  myb_attitude_function attitude;
  void* attitude_context;
  const myb_plate* plates;
  size_t plate_count;
  /* the transition to the continuum (myb_gas_bridge), if length > 0: the length (m) of the
   * object and the diameter (m) of the molecules for the Knudsen number, and the drag of the
   * continuum, -1/2 rho |w| w continuum_cd continuum_area / mass */
  double length;
  double molecule;
  double continuum_cd;
  double continuum_area;
} myb_surface_drag_term;

MYB_API MYB_NODISCARD myb_status myb_term_surface_drag(void* term, double t, const double x[3],
                                                       const double v[3], double a[3]);
/* As above; the gradients of the density and of the wind are central differences of their
 * functions, as in myb_term_drag_jacobian. */
MYB_API MYB_NODISCARD myb_status myb_term_surface_drag_jacobian(void* term, double t,
                                                                const double x[3],
                                                                const double v[3],
                                                                double dadx[9], double dadv[9]);

/* The recoil of radiation that the object sends out in a direction fixed in its body: the
 * heat of a radioisotope generator or of the electronics, the beam of an antenna —
 *   a = - efficiency power / (mass c) direction,
 * efficiency the part of the power that leaves along the direction (1 for a beam, 2/3 for a
 * flat Lambertian radiator). With decay > 0 the power falls as exp(-(t - epoch) / decay): a
 * radioactive source, or a surface that stops outgassing. */
typedef struct myb_recoil_term {
  const myb_object* object;
  double c;
  double power;        /* W */
  double efficiency;
  double direction[3]; /* unit, in the body frame: where the radiation goes */
  double epoch;        /* model time (s) */
  double decay;        /* s; 0: constant */
  myb_attitude_function attitude;
  void* attitude_context;
} myb_recoil_term;

MYB_API MYB_NODISCARD myb_status myb_term_recoil(void* term, double t, const double x[3],
                                                 const double v[3], double a[3]);
/* As above: the recoil is fixed in the body, so its derivatives are those of the attitude. */
MYB_API MYB_NODISCARD myb_status myb_term_recoil_jacobian(void* term, double t,
                                                          const double x[3], const double v[3],
                                                          double dadx[9], double dadv[9]);

#ifdef __cplusplus
}
#endif

#endif /* MYB_SURFACE_H_ */
