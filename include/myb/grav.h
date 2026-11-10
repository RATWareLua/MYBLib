/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/grav.h — gravity fields of extended bodies.
 *
 * Spherical harmonics. In the body-fixed frame, with geocentric latitude phi, longitude lambda
 * and fully normalised associated Legendre functions Pbar_nm (geodesy convention, without the
 * Condon-Shortley phase),
 *   V(r, phi, lambda) = GM/r sum_{n=0}^{N} (R/r)^n sum_{m=0}^{n} Pbar_nm(sin phi)
 *                       (Cbar_nm cos(m lambda) + Sbar_nm sin(m lambda)).
 * Evaluation: forward column recursion of Pbar_nm in the form of S. A. Holmes and
 * W. E. Featherstone (J. Geodesy 76, 279, 2002) with the sectorial seeds scaled by 2^930, so
 * that degree 2700 stays within range at every latitude; the gradient is written without any
 * division by cos(phi) (the columns Pbar_nm / cos(phi) are recursed from their own seeds, and
 * dPbar_nm/dphi = e_nm Pbar_n,m+1 - m sin(phi) Pbar_nm / cos(phi)), so the poles are regular
 * points. Coefficients come from data: ICGEM files (format 1.0 and 2.0, with the time-variable
 * terms gfct, trnd, acos, asin), or arrays.
 *
 * Polyhedra: the exact field of a homogeneous closed polyhedron (R. A. Werner and
 * D. J. Scheeres, Celest. Mech. Dyn. Astron. 65, 313, 1997), valid inside and outside. */
#ifndef MYB_GRAV_H_
#define MYB_GRAV_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/io.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Spherical harmonics ----------------------------------------------------------------- */

typedef struct myb_gravity_field myb_gravity_field;

/* Coefficient arrays are indexed by n (n + 1) / 2 + m, 0 <= m <= n <= degree. degree >= 0 and
 * <= 2700, gm > 0, radius > 0. Errors: MYB_E_DOMAIN, MYB_E_NOMEM, MYB_E_OVERFLOW. */
MYB_API MYB_NODISCARD myb_status myb_gravity_field_create(double gm, double radius, int degree,
                                                          const double* c, const double* s,
                                                          const myb_allocator* allocator,
                                                          myb_gravity_field** field);

/* From an ICGEM file (text): keys earth_gravity_constant (or gravity_constant), radius,
 * max_degree, norm (fully_normalized only), tide_system, format; lines gfc (static), and in
 * format 2.0 gfct (static part with its reference epoch), trnd (rate per year), acos and asin
 * (periodic terms with their period in years). Truncated to degree max_degree (0 = all).
 * Errors: MYB_E_FORMAT (diag gives the line), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_gravity_field_load_icgem(const char* text, size_t size,
                                                              int max_degree,
                                                              const myb_allocator* allocator,
                                                              myb_gravity_field** field,
                                                              myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_gravity_field_load_icgem_file(const char* path,
                                                                   int max_degree,
                                                                   const myb_allocator* allocator,
                                                                   myb_gravity_field** field,
                                                                   myb_diag* diag);

MYB_API void myb_gravity_field_destroy(myb_gravity_field* field); /* accepts NULL */

/* An independent copy (coefficients, time-variable terms, current corrections): fields are
 * modified by set_time and the tides, so each thread evaluating at its own time needs its own.
 * Errors: MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_gravity_field_clone(const myb_gravity_field* field,
                                                         const myb_allocator* allocator,
                                                         myb_gravity_field** copy);

MYB_API int myb_gravity_field_degree(const myb_gravity_field* field);
MYB_API double myb_gravity_field_gm(const myb_gravity_field* field);
MYB_API double myb_gravity_field_radius(const myb_gravity_field* field);

/* "tide_free", "zero_tide", "mean_tide" or "" as declared by the file. */
MYB_API const char* myb_gravity_field_tide_system(const myb_gravity_field* field);

/* Effective coefficients = static ones + the time-variable terms at t (TT, as used by the
 * models) + the corrections added since. myb_gravity_field_set_time recomputes them for t;
 * myb_gravity_field_add_correction adds dC, dS to one coefficient (tides, pole tide) until the
 * next set_time. Both modify the field: not concurrently with evaluations. Errors:
 * MYB_E_DOMAIN (n, m outside the field, non-finite values). */
MYB_API void myb_gravity_field_set_time(myb_gravity_field* field, myb_time t);
MYB_API MYB_NODISCARD myb_status myb_gravity_field_add_correction(myb_gravity_field* field,
                                                                  int n, int m, double dc,
                                                                  double ds);
MYB_API MYB_NODISCARD myb_status myb_gravity_field_coefficient(const myb_gravity_field* field,
                                                               int n, int m, double* c,
                                                               double* s);

/* Scratch memory for evaluations up to a degree; one per thread. */
typedef struct myb_gravity_work myb_gravity_work;
MYB_API MYB_NODISCARD myb_status myb_gravity_work_create(int degree, const myb_allocator* allocator,
                                                         myb_gravity_work** work);
MYB_API void myb_gravity_work_destroy(myb_gravity_work* work); /* accepts NULL */

/* Potential (m^2 s^-2, positive: GM/r for a point mass) and acceleration (m s^-2, its
 * gradient) at the body-fixed position r (m), using the degrees n_min .. n_max (n_max = -1:
 * all; n_min = 2 leaves out the central and the degree-1 terms). potential and a may be NULL.
 * The work must have been created for at least n_max. Errors: MYB_E_DOMAIN (r = 0, r non-finite,
 * degrees outside the field or the work). Thread-safe with one work per thread. */
MYB_API MYB_NODISCARD myb_status myb_gravity_field_eval(const myb_gravity_field* field,
                                                        myb_gravity_work* work, const double r[3],
                                                        int n_min, int n_max, double* potential,
                                                        double a[3]);

/* The gradient of the acceleration (the second derivatives of the potential, s^-2) at the
 * body-fixed position r, for the degrees n_min .. n_max as above: gradient[3 i + j] =
 * d a_i / d r_j, symmetric and, outside the body, free of trace. Exact to rounding like the
 * potential: every component is a series of harmonics two degrees higher (in the source),
 * with no division by cos(phi), so the poles are regular points. Errors: as
 * myb_gravity_field_eval; MYB_E_DOMAIN also for n_max above 2698. */
MYB_API MYB_NODISCARD myb_status myb_gravity_field_gradient(const myb_gravity_field* field,
                                                            myb_gravity_work* work,
                                                            const double r[3], int n_min,
                                                            int n_max, double gradient[9]);

/* ---- Tides ----------------------------------------------------------------------------------
 * Tidal changes of a field's coefficients after IERS Conventions 2010, chapter 6 (the formulas
 * are general; every number comes from data, e.g. data/grav/iers2010-tides.tk for the Earth):
 *   solid tides, step 1: dC_nm - i dS_nm = k_nm / (2n + 1) sum_j (GM_j / GM)(R / r_j)^(n+1)
 *     Pbar_nm(sin phi_j) exp(-i m lambda_j) for n = 2, 3, and degree 4 from the degree-2 tide
 *     with k_2m(+) (eq. 6.6, 6.7), GM and R those of the field;
 *   step 2: frequency-dependent corrections of C20, C21, S21, C22, S22 as Poisson series over
 *     chi = GMST + pi and the Delaunay arguments (eq. 6.8, Tables 6.5a-c);
 *   step 3: for a zero-tide field the permanent part A0 H0 Re(k20) is removed from C20
 *     (eq. 6.13, 6.14);
 *   pole tides: dC21 = c1 (m1 + a1 m2), dS21 = c2 (m2 + a2 m1), with the wobble m1 = xp - xs,
 *     m2 = -(yp - ys) about the secular pole (sections 6.4, 6.5; chapter 7, eq. 21);
 *   ocean tides: dC_nm - i dS_nm = sum_f (C+ - i S+) e^(i theta_f) + (C- + i S-) e^(-i theta_f)
 *     (eq. 6.15) from files of the form of the IERS fes2004_Cnm-Snm.dat (Doodson number, name,
 *     n, m, C+, S+, C-, S- in units of 1e-11), degree 1 skipped; S_n0 is not changed.
 * The kernel variables: <prefix>_ARGUMENT_* and the series <prefix>_C20 .. _S22 (form of
 * data/frames/iau2006-cip.tk), <prefix>_LOVE (rows n, m, Re k, Im k, k+), _PERMANENT (A0, H0),
 * _POLE_SOLID and _POLE_OCEAN (c1, a1, c2, a2), _SECULAR_POLE (xs0, xs1, ys0, ys1 in mas). */

typedef struct myb_tides myb_tides;

/* Errors: MYB_E_NOT_FOUND, MYB_E_FORMAT (diag names the variable), MYB_E_DOMAIN (prefix longer
 * than 20 characters), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_tides_create(const myb_kpool* pool, const char* prefix,
                                                  const myb_allocator* allocator,
                                                  myb_tides** tides, myb_diag* diag);
MYB_API void myb_tides_destroy(myb_tides* tides); /* accepts NULL */

/* Ocean tide coefficients up to max_degree (2 .. 360), replacing those loaded before.
 * Errors: MYB_E_FORMAT (diag gives the line), MYB_E_DOMAIN, MYB_E_NOMEM, MYB_E_IO. */
MYB_API MYB_NODISCARD myb_status myb_tides_load_ocean(myb_tides* tides, const char* text,
                                                      size_t size, int max_degree,
                                                      myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_tides_load_ocean_file(myb_tides* tides, const char* path,
                                                           int max_degree, myb_diag* diag);
MYB_API size_t myb_tides_ocean_waves(const myb_tides* tides);

/* A tide-raising body: its GM and its position in the body-fixed frame of the field (m). */
typedef struct myb_tide_raiser {
  double gm;
  double r[3];
} myb_tide_raiser;

enum {
  MYB_TIDES_SOLID = 1,      /* step 1 (and step 3) */
  MYB_TIDES_FREQUENCY = 2,  /* step 2 */
  MYB_TIDES_OCEAN = 4,
  MYB_TIDES_POLE_SOLID = 8,
  MYB_TIDES_POLE_OCEAN = 16,
  MYB_TIDES_ALL = 31
};

typedef enum myb_permanent_tide {
  MYB_PERMANENT_TIDE_AUTO = 0,   /* from the field: zero_tide -> remove, tide_free -> keep */
  MYB_PERMANENT_TIDE_KEEP = 1,
  MYB_PERMANENT_TIDE_REMOVE = 2
} myb_permanent_tide;

typedef struct myb_tides_input {
  const myb_tide_raiser* raisers;
  size_t count;
  myb_time t;     /* argument time: UT1 (chi = GMST + pi); also the date of the secular pole */
  double xp;      /* polar motion, rad (pole tides) */
  double yp;
  int flags;      /* MYB_TIDES_* */
  myb_permanent_tide permanent;
} myb_tides_input;

/* Adds the selected tidal changes to the field as corrections (myb_gravity_field_add_correction;
 * they last until the next myb_gravity_field_set_time); coefficients above the field's degree
 * are skipped. Errors: MYB_E_DOMAIN (a raiser at the origin, non-finite input, AUTO with a field
 * that is neither zero-tide nor tide-free). */
MYB_API MYB_NODISCARD myb_status myb_tides_apply(const myb_tides* tides, myb_gravity_field* field,
                                                 const myb_tides_input* input);

/* ---- Ellipsoidal harmonics ----------------------------------------------------------------------
 * Exterior expansion about a reference triaxial ellipsoid of semi-axes a > b > c, which unlike
 * the spherical-harmonic series converges down to the ellipsoid itself — inside the Brillouin
 * sphere of an elongated body (R. Garmier, J.-P. Barriot, Celest. Mech. Dyn. Astron. 79, 235,
 * 2001):
 *   U = (GM / a) sum_{n=0}^{N} sum_{p=0}^{2n} alpha_np [F_np(lambda) / F_np(a)] E_np(mu) E_np(nu)
 *       / a^(2n),
 * lambda >= a >= mu >= ... the ellipsoidal coordinates (h^2 = a^2 - b^2, k^2 = a^2 - c^2), E_np
 * the Lamé functions of the first kind with a monic polynomial factor, F_np those of the second
 * kind. Within a degree the 2n + 1 functions are ordered by class K, L, M, N (classes in the
 * order of their eigenvalue). The terms are evaluated as R(lambda) H(x, y, z) with H the
 * Cartesian solid-harmonic polynomial of Niven, so the gradient has no singular set. Valid on
 * and outside the reference ellipsoid (MYB_E_DOMAIN inside); degree <= 20. */

typedef struct myb_ellipsoidal_field myb_ellipsoidal_field;

/* axes: a > b > c > 0 (m); alpha: (degree + 1)^2 coefficients, index n^2 + p (NULL: zeros).
 * Errors: MYB_E_DOMAIN, MYB_E_NOMEM, MYB_E_NOCONV (Lamé functions not resolved). */
MYB_API MYB_NODISCARD myb_status myb_ellipsoidal_field_create(double gm, const double axes[3],
                                                              int degree, const double* alpha,
                                                              const myb_allocator* allocator,
                                                              myb_ellipsoidal_field** field);
MYB_API void myb_ellipsoidal_field_destroy(myb_ellipsoidal_field* field); /* accepts NULL */
MYB_API int myb_ellipsoidal_field_degree(const myb_ellipsoidal_field* field);
MYB_API MYB_NODISCARD myb_status myb_ellipsoidal_field_set_coefficients(
    myb_ellipsoidal_field* field, const double* alpha);

/* Potential (m^2 s^-2) and acceleration (m s^-2) at r (m, frame of the ellipsoid's axes);
 * either output may be NULL. Thread-safe. */
MYB_API MYB_NODISCARD myb_status myb_ellipsoidal_field_eval(const myb_ellipsoidal_field* field,
                                                            const double r[3], double* potential,
                                                            double a[3]);

/* The (degree + 1)^2 basis functions at r: the potential of each unit coefficient divided by
 * GM / a (for estimating alpha from data). */
MYB_API MYB_NODISCARD myb_status myb_ellipsoidal_field_basis(const myb_ellipsoidal_field* field,
                                                             const double r[3], double* values);

/* E_np(x) (m^n) for x >= 0, with |x^2 - h^2|, |x^2 - k^2| under the square roots. */
MYB_API MYB_NODISCARD myb_status myb_ellipsoidal_lame(const myb_ellipsoidal_field* field, int n,
                                                      int p, double x, double* value);

/* ---- Mass concentrations ----------------------------------------------------------------------
 * The field of count point masses gm[k] at positions[3 k .. 3 k + 2] (m, body frame): potential
 * sum gm / |p - r| and acceleration, summed with compensation in array order. potential and a
 * may be NULL. MYB_E_DOMAIN: r at a mascon or non-finite. */
MYB_API MYB_NODISCARD myb_status myb_mascons_eval(const double* gm, const double* positions,
                                                  size_t count, const double r[3],
                                                  double* potential, double a[3]);

/* ---- Two extended bodies --------------------------------------------------------------------
 * The mutual potential of two bodies A and B expands as U = G sum (-1)^l' / (l! l'!) M_A^L
 * M_B^L' d_LL'(1/r) over their STF multipole moments (T. Damour, M. Soffel, C. Xu, Phys. Rev. D
 * 43, 3273, 1991). The terms l' = 0 (resp. l = 0) are the field of A at the centre of mass of B
 * (resp. of B at A), given to any degree by myb_gravity_field_eval in the body frame; the first
 * term coupling the two figures, l = l' = 2, is
 *   U22 = G / (36 r^5) [105 (n.QA.n)(n.QB.n) - 60 n.QA QB.n + 6 tr(QA QB)],
 * with the trace-free quadrupoles Q = integral rho (3 x x^T - |x|^2 1) dV = tr(I) 1 - 3 I
 * (I: inertia tensor), r = x_B - x_A, n = r / r. Its relative size is (size_A size_B / r^2)^2
 * of the point-mass term; the next terms are of relative order (size / r)^6 for bodies
 * without octupoles, or 5 otherwise. Quadrupoles are passed as G Q (m^5 s^-2), row-major. */

/* G Q in the body frame from the fully normalised degree-2 coefficients of a field:
 * c = (C20, C21, C22), s = (S21, S22), with gm and the reference radius. */
MYB_API void myb_quadrupole_from_coefficients(double gm, double radius, const double c[3],
                                              const double s[2], double gq[9]);

/* G Q in an inertial frame from G Q in the body frame, m = the rotation inertial -> body
 * (rows, as myb_rotation.m): gq = m^T gq_body m. */
MYB_API void myb_quadrupole_rotate(const double m[9], const double gq_body[9], double gq[9]);

/* Accelerations of the centres of mass of A and B from U22 (both in the frame of r and the
 * quadrupoles), ADDED to a_a and a_b (either may be NULL). MYB_E_DOMAIN: r = 0, non-finite
 * input, gm <= 0. */
MYB_API MYB_NODISCARD myb_status myb_figure_figure(double gm_a, const double gq_a[9], double gm_b,
                                                   const double gq_b[9], const double r[3],
                                                   double a_a[3], double a_b[3]);

/* ---- Polyhedra ----------------------------------------------------------------------------- */

typedef struct myb_polyhedron myb_polyhedron;

/* A closed, consistently oriented triangular mesh (counter-clockwise seen from outside):
 * vertices (3 per vertex, m) and faces (3 vertex indices per face), of density so that
 * G rho = g_rho (m^3 kg^-1 s^-2 times kg m^-3 = s^-2). Errors: MYB_E_DOMAIN (an index out of
 * range, a degenerate face, a mesh that is not closed: every edge must be shared by exactly two
 * faces in opposite directions), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_polyhedron_create(const double* vertices, size_t vertex_count,
                                                       const size_t* faces, size_t face_count,
                                                       double g_rho,
                                                       const myb_allocator* allocator,
                                                       myb_polyhedron** polyhedron);
MYB_API void myb_polyhedron_destroy(myb_polyhedron* polyhedron); /* accepts NULL */

/* Volume (m^3) and centre of mass (m) of the mesh. */
MYB_API double myb_polyhedron_volume(const myb_polyhedron* polyhedron, double center[3]);

/* Potential (positive), acceleration and the Laplacian (-4 pi G rho inside, 0 outside, from the
 * solid angles) at r (m) in the frame of the mesh; NULL outputs are skipped. The potential and
 * the acceleration are continuous everywhere, the surface included (edge and face terms take
 * their limits); on the surface the Laplacian is that of the local solid angle. The terms grow
 * like r^2 while the potential falls like 1/r: the relative rounding error is about
 * 1e-16 (r / size)^2, so far from the body use its spherical-harmonic expansion.
 * Errors: MYB_E_DOMAIN (non-finite r). */
MYB_API MYB_NODISCARD myb_status myb_polyhedron_eval(const myb_polyhedron* polyhedron,
                                                     const double r[3], double* potential,
                                                     double a[3], double* laplacian);

#ifdef __cplusplus
}
#endif

#endif /* MYB_GRAV_H_ */
