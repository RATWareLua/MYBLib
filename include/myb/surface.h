/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/surface.h — the object as more than a sphere: a surface of plates or a mesh, turned by
 * an attitude law, under light and under a rarefied gas or a dense one; the recoil of what it
 * emits; the heat of an entry and the mass it burns away.
 *
 * Light on plates is myb_force_radiation_plates (myb/force.h); here are the same for a mesh
 * with the shadows its parts cast on each other, the force and the lift of a free-molecular
 * flow (Sentman), the hypersonic continuum of the modified Newtonian theory and the bridge
 * between the two, the attitude laws, the terms of a force model (myb/model.h) that put them
 * together with an environment (myb/env.h), and the heating and the ablation of an object
 * that enters an atmosphere. */
#ifndef MYB_SURFACE_H_
#define MYB_SURFACE_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/env.h"
#include "myb/force.h"
#include "myb/grav.h"
#include "myb/model.h"

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

/* The force and the torque of light on the mesh at the moment, as myb_force_radiation_mesh,
 * but with the heat conducted into the surface: the body turns at a steady rate about `axis`
 * (unit, in the frame of the mesh, in the sense of the rotation) and the source stands still
 * over a rotation; every face — every part of it with shadows — re-emits what it absorbs as a
 * Lambertian radiator later and smoother: its absorbed power over the rotation, taken at
 * `phases` equal steps (at least 8) from the moment on, is split into its harmonics k of the
 * rotation, each answered by 1 / (1 + theta sqrt(k) (1 + i) / 2) as in myb_force_yorp (theta
 * the thermal parameter at the rotation frequency), and what it emits now is the sum at the
 * phase of the moment (all the harmonics below that of Nyquist). The momentum of the light absorbed and of the light reflected is that
 * of the moment. With theta = 0 this is myb_force_radiation_mesh; over a rotation its torque
 * averages to that of myb_force_yorp. Errors: MYB_E_DOMAIN (as myb_force_radiation_mesh; an
 * axis that is not unit to 1e-12, theta negative or not finite, phases outside 8 .. 256). */
MYB_API MYB_NODISCARD myb_status myb_force_radiation_mesh_lagged(
    const myb_polyhedron* shape, const myb_plate* surface, double pressure, const double u[3],
    const double axis[3], double theta, int phases, int shadows, double force[3],
    double torque[3]);

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
 * 0 below Kn = 1e-3 and 1 above 10. The continuum value is the modified Newtonian flow below
 * or a number of the caller's (myb_surface_drag_term). */
MYB_API double myb_gas_knudsen(double density, double molar_mass, double diameter,
                               double length);
MYB_API double myb_gas_bridge(double knudsen);

/* ---- Continuum flow --------------------------------------------------------------------------
 * In air dense enough to flow round the body, the hypersonic pressure of the modified
 * Newtonian theory (L. Lees, Hypersonic flow, Fifth International Aeronautical Conference,
 * Los Angeles, 1955; J. D. Anderson, Hypersonic and High-Temperature Gas Dynamics, 2nd ed.,
 * AIAA, 2006, ch. 3): an element of the surface whose outward normal n makes the angle
 * 90 deg - theta with the motion v relative to the air takes
 *   p - p_inf = (1/2) rho V^2 C_p,max sin^2 theta,   sin theta = n . v / V,
 * along -n where it faces the flow (sin theta > 0) and nothing where it does not: the flow is
 * taken to leave the surface there, and the base pressure of a real body is near p_inf.
 * C_p,max is the stagnation pressure behind a normal shock in a perfect gas over the dynamic
 * pressure,
 *   C_p,max = 2 / (gamma M^2) (p_02 / p_inf - 1),
 *   p_02 / p_inf = [(gamma + 1)^2 M^2 / (4 gamma M^2 - 2 (gamma - 1))]^(gamma / (gamma - 1))
 *                  (1 - gamma + 2 gamma M^2) / (gamma + 1)
 * (Rayleigh's pitot formula; the tables of NACA Report 1135, 1953): for gamma = 1.4 it is
 * 1.8394 as M -> infinity, 1.8088 at M = 5, 1.6573 at M = 2 (Newton's own theory has 2). A
 * sphere then has C_D = C_p,max / 2 on its cross-section, a cone of half-angle delta along its
 * axis C_p,max sin^2 delta on its base, a flat plate at the incidence alpha the normal force
 * C_p,max sin^2 alpha. The theory is one of hypersonic speeds. Below M = 1 there is no shock
 * and C_p,max is the isentropic stagnation pressure, 2 / (gamma M^2) ((1 + (gamma - 1) M^2 /
 * 2)^(gamma / (gamma - 1)) - 1), which meets the other at M = 1 and tends to 1 as M -> 0: a
 * continuation, so that a descent through the speed of sound has a force at all — it does not
 * make the theory right there. */

/* The air at a place: what the continuum needs of it besides the density. The speed of sound
 * sqrt(gamma R T / molar mass) gives the Mach number, gamma the pressure behind the shock.
 * An atmosphere model gives the temperature and the composition (NRLMSISE-00 of ext/: the
 * temperature at the height and the number densities of its species). */
typedef struct myb_air {
  double temperature; /* K, > 0 */
  double molar_mass;  /* kg/mol, > 0: the mean of the mixture */
  double heat_ratio;  /* gamma = c_p / c_v, > 1: 7/5 for N2 and O2, 5/3 for O, N, He, Ar */
} myb_air;

/* The air at the body-fixed position r (m) and the instant t, in the form of a
 * myb_density_function. Any status other than MYB_OK stops the evaluation. */
typedef myb_status (*myb_air_function)(void* context, myb_time t, const double r[3],
                                       myb_air* air);

/* gamma of a mixture of ideal gases with the given amounts (in any one unit: mole fractions,
 * number densities) of monatomic and of diatomic molecules whose vibrations are not excited:
 * (5 n1 + 7 n2) / (3 n1 + 5 n2). 0 if the amounts are not valid (negative, not finite, no
 * gas at all). */
MYB_API double myb_air_heat_ratio(double monatomic, double diatomic);

/* C_p,max at the Mach number (>= 0) in a gas of the heat ratio gamma (> 1). Errors:
 * MYB_E_DOMAIN. */
MYB_API MYB_NODISCARD myb_status myb_newtonian_cp_max(double mach, double heat_ratio,
                                                      double* cp_max);

/* a += -(1/2) rho V^2 C_p,max / mass sum A_i sin^2 theta_i n_i over the plates that face the
 * flow (area and normal used), v the velocity of the object relative to the air in the frame
 * of the plates, M = V / sqrt(gamma R T / molar mass). Plates have no places: none of them
 * hides another. Errors: MYB_E_DOMAIN (an air, a mass, a density or a plate that is not
 * valid; v = 0). */
MYB_API MYB_NODISCARD myb_status myb_force_newtonian_plates(double density, const myb_air* air,
                                                            double mass, const myb_plate* plates,
                                                            size_t count, const double v[3],
                                                            double a[3]);

/* The same pressure on a closed mesh (myb/grav.h): force (N) and torque (N m, about the origin
 * of the mesh) are written, either may be NULL; v as above, in the frame of the mesh.
 * shadows == 0: every face turned to the flow takes it. shadows = k >= 1: each such face is
 * divided into k^2 equal triangles, and a part counts only if nothing of the mesh stands
 * upstream of its centre — the theory leaves the parts that the body hides from the flow
 * without pressure, as it leaves the faces turned away. The pressure is that of light of the
 * pressure (1/4) rho V^2 C_p,max on a mirror, 2 p sin^2 theta along -n, and the faces hidden
 * from the flow are those in the shadow: the computation is myb_force_radiation_mesh with a
 * fully specular surface, its rays and its cost. Errors: MYB_E_DOMAIN. */
MYB_API MYB_NODISCARD myb_status myb_force_newtonian_mesh(const myb_polyhedron* shape,
                                                          double density, const myb_air* air,
                                                          const double v[3], int shadows,
                                                          double force[3], double torque[3]);

/* The continuum of a myb_surface_drag_term in place of its one number continuum_cd: the
 * modified Newtonian flow on the plates of the term or on a mesh. The air is that of the
 * function if one is given — and it then serves the whole term: its temperature and molar
 * mass replace those of the term's gas for the free molecules and for the Knudsen number too
 * (the wall temperature and the accommodation stay) —, otherwise the temperature and the
 * molar mass of the term's gas with heat_ratio. */
typedef struct myb_continuum {
  double heat_ratio;           /* gamma of the term's gas, when air is NULL */
  myb_air_function air;        /* may be NULL */
  void* air_context;
  const myb_polyhedron* shape; /* NULL: the plates of the term take the pressure */
  int shadows;                 /* of the mesh, as in myb_force_newtonian_mesh */
} myb_continuum;

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
/* The Jacobians of the three terms of this header are analytic (but that of the drag with a
 * continuum, below): the force in the body of the object is differentiated as it is written,
 * and so is the attitude when its law is one of the four of this header (it is known by its
 * address). The derivatives of an attitude law
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
 * (myb_force_gas_plates): density and wind as in myb_drag_term. With length > 0 the
 * transition to the continuum, a = W a_free + (1 - W) a_continuum with the weight
 * W = myb_gas_bridge(Kn) of the free molecules, the continuum
 *  - without `continuum`: the drag -1/2 rho |w| w continuum_cd continuum_area / mass of the
 *    one number given (the term as it was before the continuum flow came, to the bit);
 *  - with it: the modified Newtonian flow (myb_force_newtonian_plates on the plates of the
 *    term, or myb_force_newtonian_mesh on continuum->shape) in the air of myb_continuum; the
 *    part whose weight is 0 is not computed.
 * The free molecules act on the plates only: a term whose continuum is a mesh and that has
 * no plates feels nothing in a free-molecular flow. */
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
  /* NULL: continuum_cd and continuum_area as above; otherwise the modified Newtonian flow
   * (myb_continuum), and continuum_cd and continuum_area are not used */
  const myb_continuum* continuum;
} myb_surface_drag_term;

MYB_API MYB_NODISCARD myb_status myb_term_surface_drag(void* term, double t, const double x[3],
                                                       const double v[3], double a[3]);
/* As above; the gradients of the density and of the wind are central differences of their
 * functions, as in myb_term_drag_jacobian. With a continuum the Jacobian is central
 * differences of the term itself, over 6e-6 of |x| and of |v| (12 calls of the term): its
 * force goes through the Mach number, the air of the place and, on a mesh, the parts the flow
 * reaches, which change in steps — a step of the differences across such a change gives the
 * slope of the step. */
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

/* ---- Entry: the heat and the ablation ----------------------------------------------------------
 * The convective heating of the stagnation point (K. Sutton, R. A. Graves, A general
 * stagnation-point convective-heating equation for arbitrary gas mixtures, NASA TR R-376,
 * 1971):
 *   q = k sqrt(rho / R_n) V^3   (W m^-2),
 * R_n the radius of the nose, V the speed relative to the air; k depends on the composition
 * of the gas alone, 1.7415e-4 kg^1/2 m^-1 for the air of the Earth. The heat that the shock
 * layer radiates, which takes over above some 10 km/s for a large nose, is not in it.
 * Errors: MYB_E_DOMAIN (a negative or non-finite value, a radius that is not positive). */
MYB_API MYB_NODISCARD myb_status myb_heating_sutton_graves(double k, double density,
                                                           double nose_radius, double speed,
                                                           double* q);

/* An object that the heat consumes loses mass at
 *   dm/dt = -q A / Q*,
 * q of Sutton and Graves at the state, A the area that takes it, Q* the effective heat of
 * ablation (J/kg: what a kilogram of the surface carries away — its heating, melting and
 * vaporisation, and the heat that the gases it gives off keep from the wall).
 *
 * The mass lost depends on the path — the density along it and the speed — so no function of
 * time alone can give it before the path is known, and a term cannot keep it (what a term
 * keeps depends on t only, and an integrator withdraws steps). The mass is therefore a seventh
 * component of the state, y = (x, v, m), as in myb_qlaw_flight_rhs (myb/thrust.h), and the
 * forces see it through the object: with myb_ablation_mass as the mass function of the object
 * (myb_object.mass_function, mass_context = this structure), every term of the field — the
 * drag, the free molecules and the continuum, the light — acts on the mass of the state that
 * myb_ablation_rhs evaluates, which it writes here before it asks the field. The radius of the
 * nose and the area stay as given while the mass falls. */
typedef struct myb_ablation {
  myb_env* env;
  size_t body; /* whose air */
  myb_density_function density;
  void* density_context;
  myb_wind_function wind; /* may be NULL, as in myb_drag_term */
  void* wind_context;
  double k;               /* of Sutton and Graves, kg^1/2 m^-1 */
  double nose_radius;     /* R_n, m */
  double area;            /* A, m^2 */
  double heat;            /* Q*, J/kg */
  myb_force_model* field; /* the forces on the object; NULL: none, the motion is a straight line */
  double mass;            /* kg: the mass of the state being evaluated (myb_ablation_rhs) */
} myb_ablation;

/* The heating (W m^-2) and the rate of the mass (kg/s, <= 0) at the state (x, v) relative to
 * the origin of the environment, at the model time t; either output may be NULL. Errors:
 * MYB_E_DOMAIN (values of the structure that are not valid), those of the environment and of
 * the density and wind functions. */
MYB_API MYB_NODISCARD myb_status myb_ablation_rate(const myb_ablation* ablation, double t,
                                                   const double x[3], const double v[3],
                                                   double* heating, double* rate);

/* The right-hand side for an integrator of first-order systems (myb_ode_function of
 * myb/integ.h; the context is a myb_ablation; 7 components): dx/dt = v, dv/dt = the field,
 * dm/dt = -q A / Q*. Errors: MYB_E_DOMAIN (a mass that is not positive or not finite: the
 * object is consumed), those of the field and of myb_ablation_rate. */
MYB_API MYB_NODISCARD myb_status myb_ablation_rhs(void* ablation, double t, const double* y,
                                                  double* dydt);

/* The mass function of an object that ablates (myb_mass_function; the context is the
 * myb_ablation): the mass of the state being evaluated, whatever t. */
MYB_API double myb_ablation_mass(void* ablation, double t);

#ifdef __cplusplus
}
#endif

#endif /* MYB_SURFACE_H_ */
