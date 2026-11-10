/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/force.h — forces on a body other than its own field's gravity: perturbations by point
 * masses, radiation of a source (direct, reflected and emitted by a planet) and its shadows,
 * atmospheric drag, empirical accelerations.
 *
 * Every function is a pure function of its arguments (positions and velocities in one frame,
 * SI units); none knows a particular body: luminosities, radii, coefficients of a planet's
 * albedo, densities of an atmosphere come from the caller. Accelerations are ADDED to a[3], so
 * that a model is the sum of its calls. */
#ifndef MYB_FORCE_H_
#define MYB_FORCE_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/grav.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Point masses -------------------------------------------------------------------------- */

/* Perturbation of a third body of gravitational parameter gm at s (m, relative to the central
 * body) on an object at r (relative to the central body): the direct term minus the indirect
 * one, a += gm ((s - r)/|s - r|^3 - s/|s|^3), evaluated in R. H. Battin's form (An
 * Introduction to the Mathematics and Methods of Astrodynamics, 1999, §8.4) without the
 * cancellation of the two terms when |r| << |s|. */
MYB_API void myb_force_third_body(double gm, const double s[3], const double r[3], double a[3]);

/* ---- Radiation ------------------------------------------------------------------------------- */

/* Flux (W m^-2) of a source of luminosity L (W) at distance d (m): L / (4 pi d^2). The pressure
 * on a perfect absorber is flux / c. */
MYB_API double myb_radiation_flux(double luminosity, double distance);

/* A sphere (or an object treated as one): a += pressure * cr * area_over_mass * u, with u the
 * unit direction of propagation (from the source to the object), cr = 1 + reflectivity. */
MYB_API void myb_force_radiation_sphere(double pressure, double cr, double area_over_mass,
                                        const double u[3], double a[3]);

/* A flat plate: area (m^2), outward unit normal, fractions of the incoming light reflected
 * specularly and diffusely (the rest is absorbed); thermal != 0 re-emits the absorbed part at
 * once as a Lambertian emitter on the lit side. */
typedef struct myb_plate {
  double area;
  double normal[3];
  double specular;
  double diffuse;
  int thermal;
} myb_plate;

/* Sum over plates lit from the front (n.u < 0) of
 *   F = pressure A cos(theta) [ (1 - specular) u - 2 (specular cos(theta) + diffuse / 3) n ]
 *       (- (2/3) absorbed A cos(theta) pressure n with thermal re-emission),
 * cos(theta) = -n.u (D. Vallado, Fundamentals of Astrodynamics, 2013, §9.6.3; Montenbruck &
 * Gill, Satellite Orbits, 2000, eq. 3.75); a += F / mass. Plates, u and a share one frame
 * (usually the body frame of the object). No self-shadowing. MYB_E_DOMAIN: mass <= 0, a plate
 * with a negative area, a normal that is not unit to 1e-12, fractions outside [0, 1] or summing
 * to more than 1. */
MYB_API MYB_NODISCARD myb_status myb_force_radiation_plates(double pressure, double mass,
                                                            const myb_plate* plates, size_t count,
                                                            const double u[3], double a[3]);

/* Fraction (0 .. 1) of the light of a spherical source still reaching an object when an
 * occulting sphere stands in front: source and occulter are positions relative to the object
 * (m), with their radii. Exact geometry on the sky sphere (the occulter's disk is a spherical
 * cap, not a plane circle, which matters close to a planet), and the brightness of the source's
 * disk I(mu) = 1 - u1 (1 - mu) - u2 (1 - mu)^2 (limb darkening; limb = NULL: uniform). Inside
 * the occulter: 0; the occulter not in front (beyond the source or behind the object): 1.
 * Integration over the source's disk by annuli, whose occulted arcs are exact, with Gauss-
 * Legendre rules on the intervals between the kinks of the integrand: relative error < 1e-12. */
MYB_API double myb_shadow_fraction(const double source[3], double source_radius,
                                   const double occulter[3], double occulter_radius,
                                   const double limb[2]);

/* Radiation of a planet: sunlight reflected by its surface (albedo, Lambertian) and its own
 * thermal emission (emissivity, the absorbed flux re-emitted uniformly over the sphere), both
 * coefficients zonal functions a = a0 + a1 P1(sin phi) + a2 P2(sin phi) of the latitude, with
 * a1 = c0 + c1 cos(w (t - t0)) + c2 sin(w (t - t0)) (seasonal) — the model of G. Knocke,
 * J. Ries, B. Tapley (AIAA 88-4292, 1988). The visible cap of the planet is divided into a
 * central disk and `rings` rings of 6 k elements of equal projected weight. */
typedef struct myb_planet_radiation {
  double radius;         /* m */
  double albedo[5];      /* a0, c0, c1, c2, a2 */
  double emissivity[5];  /* e0, k0, k1, k2, e2 */
  double omega;          /* rad s^-1 of the seasonal term */
} myb_planet_radiation;

/* r: object, sun: source, both in the planet's body-fixed frame (z = rotation axis), m;
 * solar_flux: the flux of the source at the planet (W m^-2); phase = omega (t - t0) (rad);
 * cr_area_over_mass: cr A / m of the object treated as a sphere; rings: 1 .. 50.
 * a += acceleration (same frame). MYB_E_DOMAIN: object inside the planet, invalid rings or
 * values. */
MYB_API MYB_NODISCARD myb_status myb_force_planet_radiation(const myb_planet_radiation* planet,
                                                            const double r[3], const double sun[3],
                                                            double solar_flux, double phase,
                                                            double cr_area_over_mass, int rings,
                                                            double a[3]);

/* ---- Atmosphere -------------------------------------------------------------------------------- */

/* Drag: a += -1/2 rho cd_area_over_mass |v| v, with v the velocity relative to the atmosphere
 * (e.g. v - w x r for an atmosphere rotating with the planet) and rho in kg m^-3. */
MYB_API void myb_force_drag(double density, double cd_area_over_mass, const double v[3],
                            double a[3]);

/* ---- Empirical ------------------------------------------------------------------------------- */

/* An acceleration given in the radial, transverse (along-track for a circular orbit) and normal
 * directions of the orbit (r, v): a += a_rtn[0] R + a_rtn[1] T + a_rtn[2] N. MYB_E_DOMAIN: r
 * and v parallel or zero. */
MYB_API MYB_NODISCARD myb_status myb_force_empirical_rtn(const double r[3], const double v[3],
                                                         const double a_rtn[3], double a[3]);

/* ---- Charge and plasma ------------------------------------------------------------------------
 * Magnetic fields of bodies as Gauss coefficients: the potential of the internal field
 *   V = a sum_{n>=1} (a/r)^(n+1) sum_m (g_nm cos m lambda + h_nm sin m lambda) P_n^m(sin phi),
 * P_n^m Schmidt semi-normalised (no Condon-Shortley phase), B = -grad V in the body-fixed frame
 * (a the reference radius, phi the geocentric latitude). Evaluated by the spherical-harmonic code
 * of myb/grav.h (P_n^m = Pbar_nm / sqrt(2n + 1)), so a work of myb_gravity_work_create serves
 * both. Coefficients are indexed n (n + 1) / 2 + m (n = 0 unused), in tesla. */

typedef struct myb_magnetic_field myb_magnetic_field;

/* Errors: MYB_E_DOMAIN (degree < 1, radius <= 0, non-finite coefficients), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_magnetic_field_create(double radius, int degree,
                                                           const double* g, const double* h,
                                                           const myb_allocator* allocator,
                                                           myb_magnetic_field** field);

/* From a coefficient table in the IGRF format (columns g/h n m, then the field at each epoch
 * of the header line "g/h n m 1900.0 ... 2025.0 2025-30" in nT, the last column the secular
 * variation in nT/year): linear in time between the epochs, the secular variation after the last
 * (the IGRF prescription), at the decimal year `year` (>= the first epoch). The reference radius
 * (6371.2 km for the IGRF) is the caller's. Errors: MYB_E_FORMAT (diag gives the line),
 * MYB_E_RANGE (year before the table), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_magnetic_field_load_table(const char* text, size_t size,
                                                               double radius, double year,
                                                               const myb_allocator* allocator,
                                                               myb_magnetic_field** field,
                                                               myb_diag* diag);

MYB_API void myb_magnetic_field_destroy(myb_magnetic_field* field); /* accepts NULL */
MYB_API int myb_magnetic_field_degree(const myb_magnetic_field* field);

/* B (T) at the body-fixed position r (m); work created for at least the degree. Errors:
 * MYB_E_DOMAIN (r = 0 or non-finite). */
MYB_API MYB_NODISCARD myb_status myb_magnetic_field_eval(const myb_magnetic_field* field,
                                                         myb_gravity_work* work,
                                                         const double r[3], double b[3]);

/* Lorentz force per unit mass on a charge-to-mass ratio q/m (C/kg) moving with velocity v in
 * the fields E (V/m) and B (T) of the same inertial frame: a = (q/m) (E + v x B). */
MYB_API void myb_lorentz_acceleration(double charge_to_mass, const double e[3], const double v[3],
                                      const double b[3], double a[3]);

/* The electric field of a plasma co-rotating with angular velocity omega (ideal MHD, as in a
 * planet's inner magnetosphere), in the inertial frame: E = -(omega x r) x B. */
MYB_API void myb_corotation_field(const double omega[3], const double r[3], const double b[3],
                                  double e[3]);

/* Momentum transfer by a radial plasma stream (the solar wind): density rho = rho0 (r0 / |r|)^2
 * (rho0 at the distance r0) flowing at the speed `speed` away from the source; r and v are the
 * object's position and velocity relative to the source. With the relative velocity
 * w = speed r/|r| - v: a = (cd / 2) (area / mass) rho |w| w (cd ~ 2 for full momentum
 * transfer). */
MYB_API void myb_force_plasma_stream(double rho0, double r0, double speed,
                                     double cd_area_over_mass, const double r[3],
                                     const double v[3], double a[3]);

/* ---- Charging in a plasma ------------------------------------------------------------------
 * A conducting sphere of radius a at the potential phi (V, relative to the undisturbed plasma)
 * collects the particles of each species s (charge q, mass m, density n, temperature T as
 * k T in joules, bulk speed u relative to the sphere) whose orbits meet it: the
 * orbit-motion-limited (OML) current (H. M. Mott-Smith, I. Langmuir, Phys. Rev. 28, 727, 1926;
 * for drifting plasmas E. C. Whipple, Rep. Prog. Phys. 44, 1197, 1981), with v_th = sqrt(2 k T
 * / m), M = u / v_th, chi = q phi / (k T):
 *   attracted (chi <= 0): I = q n pi a^2 v_th [(M + 1/(2M)) erf M + e^(-M^2) / sqrt(pi)
 *                               - chi erf(M) / M],
 *   repelled  (chi > 0):  I = q n pi a^2 v_th / (sqrt(pi) M) int_sqrt(chi)^inf (x^2 - chi)
 *                               [e^(-(x - M)^2) - e^(-(x + M)^2)] dx   (closed form in erfc),
 * which for M -> 0 are 4 pi a^2 q n sqrt(k T / (2 pi m)) times (1 - chi) and e^(-chi). Sunlight
 * releases photoelectrons from the cross-section: I = pi a^2 j (j the emitted current density,
 * A/m^2, at the body's distance), times (1 + chi_p) e^(-chi_p) when phi > 0 holds them back
 * (chi_p = e phi / (k T_p)). The total current makes the body float at the potential where it
 * vanishes; the charge is then Q = 4 pi eps0 a phi (a sphere far from others), and the
 * charge-to-mass ratio for myb_lorentz_acceleration Q / m. Secondary emission by electron
 * impact (important for some materials above ~100 eV) and the thick-sheath limit (a much
 * larger than the Debye length) are not modelled. */
typedef struct myb_plasma_species {
  double charge;      /* C (-e for electrons) */
  double mass;        /* kg */
  double density;     /* m^-3 */
  double temperature; /* k T, J */
  double drift;       /* bulk speed relative to the body, m/s (>= 0) */
} myb_plasma_species;

typedef struct myb_photoemission {
  double current_density; /* emitted current density of the sunlit cross-section, A/m^2 */
  double temperature;     /* k T of the photoelectrons, J */
  double charge;          /* elementary charge, C (MYB_CONST_E) */
} myb_photoemission;

/* Total current (A; positive charges the body positively) to a sphere of radius a at the
 * potential phi from `count` species and, if photo is not NULL, photoemission. Errors:
 * MYB_E_DOMAIN (invalid input). */
MYB_API MYB_NODISCARD myb_status myb_charging_current(const myb_plasma_species* species,
                                                      size_t count,
                                                      const myb_photoemission* photo,
                                                      double radius, double potential,
                                                      double* current);

/* Floating potential (V): the root of myb_charging_current, which decreases monotonically with
 * phi (bracketing outwards from 0 in steps of the largest k T / |q|, then bisection to the
 * resolution of the doubles). Errors: MYB_E_DOMAIN (invalid input, or currents of one sign
 * only, so that no potential balances them), MYB_E_NOCONV. */
MYB_API MYB_NODISCARD myb_status myb_charging_equilibrium(const myb_plasma_species* species,
                                                          size_t count,
                                                          const myb_photoemission* photo,
                                                          double radius, double* potential);

#ifdef __cplusplus
}
#endif

#endif /* MYB_FORCE_H_ */
