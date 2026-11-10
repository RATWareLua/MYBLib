/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/frames.h — reference frames: rotations with their time derivatives and the providers of
 * the orientation of bodies.
 *
 * A rotation maps coordinates in frame A to frame B: x_B = M x_A, and states with
 * v_B = M v_A + dM/dt x_A. M is a proper orthogonal 3x3 matrix stored by rows; dM/dt is in
 * s^-1. Composition and inversion keep the derivative exact to rounding.
 *
 * Binary PCK files (NAIF "PCK Required Reading"): the orientation of a body-fixed frame
 * relative to a base frame as three Euler angles (phi, delta, w) in radians, such that
 *   M = R3(w) R1(delta) R3(phi),      R1, R3: rotations of the frame about x and z,
 * read from DAF segments of types 2, 3 and 20 (high-precision Earth orientation, the lunar
 * principal-axis frame of the JPL ephemerides). The identifiers are the file's (frame class
 * IDs, NAIF frame codes); the library attaches no meaning to them. */
#ifndef MYB_FRAMES_H_
#define MYB_FRAMES_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"
#include "myb/io.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Rotations --------------------------------------------------------------------------- */

typedef struct myb_rotation {
  double m[9];  /* rows of M */
  double dm[9]; /* rows of dM/dt (s^-1) */
} myb_rotation;

/* out = M in (positions, 3 values) — in and out may alias. */
MYB_API void myb_rotation_apply_position(const myb_rotation* rotation, const double in[3],
                                         double out[3]);

/* out = (M r, M v + dM/dt r) for a state in = (r, v) — in and out may alias. */
MYB_API void myb_rotation_apply_state(const myb_rotation* rotation, const double in[6],
                                      double out[6]);

/* out = second after first (A -> B by first, B -> C by second); out may alias either. */
MYB_API void myb_rotation_compose(const myb_rotation* second, const myb_rotation* first,
                                  myb_rotation* out);

/* The inverse rotation B -> A: the transposes of M and dM/dt; out may alias rotation. */
MYB_API void myb_rotation_invert(const myb_rotation* rotation, myb_rotation* out);

/* M = R3(angles[2]) R1(angles[1]) R3(angles[0]) and its derivative for the angle rates
 * (rad/s). */
MYB_API void myb_rotation_euler313(const double angles[3], const double rates[3],
                                   myb_rotation* rotation);

/* ---- Celestial intermediate pole and origin ---------------------------------------------------
 * The CIO-based transformation between a celestial frame (GCRS) and the terrestrial frame of
 * a body (ITRS for the Earth), IERS Conventions 2010, ch. 5:
 *   [terrestrial] = W(t) R3(ERA) C(t) [celestial],
 *   C = R3(-s) B^T,  B = | 1 - aX^2   -aXY      X        |   a = 1 / (1 + sqrt(1 - X^2 - Y^2))
 *                        | -aXY      1 - aY^2   Y        |   (eq. 5.10)
 *                        | -X        -Y         1 - a(X^2 + Y^2) |
 *   W = R1(-y_p) R2(-x_p) R3(s')   (polar motion, eq. 5.3),
 * with X, Y the coordinates of the CIP and s the CIO locator given by Poisson series (eq. 5.16),
 * ERA the rotation angle (eq. 5.15) and s' the TIO locator (eq. 5.13). Every coefficient comes
 * from a kernel (data/frames/iau2006-cip.tk: IAU 2006/2000A, prefix "IAU2006A"). */

typedef struct myb_cip myb_cip;

/* Reads the model <prefix>_* (format in data/frames/iau2006-cip.tk); the pool may be destroyed
 * afterwards. Errors: MYB_E_NOT_FOUND, MYB_E_FORMAT (diag names the variable), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_cip_create(const myb_kpool* pool, const char* prefix,
                                                const myb_allocator* allocator, myb_cip** cip,
                                                myb_diag* diag);
MYB_API void myb_cip_destroy(myb_cip* cip); /* accepts NULL */

/* X, Y and s (rad) at the TT instant tt, and their rates (rad per TT second; rates may be
 * NULL). Thread-safe. */
MYB_API void myb_cip_xys(const myb_cip* cip, myb_time tt, double xys[3], double rates[3]);

/* Rotation angle (rad, in [0, 2 pi)) at the UT1 instant ut1, and its rate (rad per UT1 second;
 * may be NULL). The fraction of the day is kept apart, so the angle is exact to ~1e-16 rad. */
MYB_API double myb_cip_era(const myb_cip* cip, myb_time ut1, double* rate);

/* TIO locator s' (rad) at the TT instant tt. */
MYB_API double myb_cip_sprime(const myb_cip* cip, myb_time tt);

/* C (celestial -> intermediate) from X, Y, s and their rates (NULL: zero rates). */
MYB_API void myb_cip_c2i(const double xys[3], const double rates[3], myb_rotation* rotation);

/* Earth orientation parameters at one instant (IERS conventions): polar motion x_p, y_p (rad),
 * UT1 - UTC (s), excess length of day LOD (s), celestial pole offsets dX, dY (rad). */
typedef struct myb_eop_values {
  double xp;
  double yp;
  double dut1;
  double lod;
  double dx;
  double dy;
} myb_eop_values;

/* Daily EOP series of the IERS and their interpolation.
 *
 * Files: EOP 20 C04 (IERS Earth Orientation Centre, "eopc04.1962-now": x, y, UT1-UTC, dX, dY,
 * LOD at 0h UTC) and finals2000A (IERS Rapid Service / Prediction Centre, Bulletin A values
 * including predictions; blank LOD, dX, dY columns are kept as missing). dX, dY must refer to
 * the precession-nutation model in use (IAU 2006/2000A for data/frames/iau2006-cip.tk).
 *
 * Interpolation (myb_eop_at): 4-point Lagrange, as the IERS routine interp.f, on the rows
 * j-1 .. j+2 around the date, in the time of the caller's scale (continuous, unlike UTC
 * labels); UT1 is interpolated as UT1 minus that scale, not UT1 - UTC, so that a leap second
 * inside the 4 days causes no jump. A quantity missing at any node of the interpolation counts as zero (dX, dY, LOD of
 * predictions); missing x, y or UT1 - UTC is an error. Optionally the diurnal and semidiurnal
 * variations absent from daily series are added (ocean tides, libration;
 * data/frames/eop-tides.tk), the LOD variation following from the UT1 one. */
typedef struct myb_eop myb_eop;

MYB_API MYB_NODISCARD myb_status myb_eop_create(const myb_allocator* allocator, myb_eop** eop);
MYB_API void myb_eop_destroy(myb_eop* eop); /* accepts NULL */

/* Replace the table (atomically; on error it is unchanged). Errors: MYB_E_FORMAT (diag gives
 * the line), MYB_E_NOMEM, MYB_E_IO for files. Rows must have increasing MJD. */
MYB_API MYB_NODISCARD myb_status myb_eop_load_c04(myb_eop* eop, const char* text, size_t size,
                                                  myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_eop_load_c04_file(myb_eop* eop, const char* path,
                                                       myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_eop_load_finals(myb_eop* eop, const char* text, size_t size,
                                                     myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_eop_load_finals_file(myb_eop* eop, const char* path,
                                                          myb_diag* diag);

/* Diurnal/semidiurnal corrections from kernel variables <prefix>_* (series form of
 * data/frames/eop-tides.tk with prefix "IERS10_EOPT", series XP, YP, UT1; kernel names are
 * limited to 32 characters, so prefixes to 12). Errors as myb_cip_create. */
MYB_API MYB_NODISCARD myb_status myb_eop_set_tides(myb_eop* eop, const myb_kpool* pool,
                                                   const char* prefix, myb_diag* diag);

/* First and last MJD (UTC) of the table, and its number of rows (0 if empty). */
MYB_API size_t myb_eop_rows(const myb_eop* eop, double* first, double* last);

/* EOP at the instant t of `scale` (a scale of `system` that its UTC table converts to, e.g.
 * TT), and the same instant on UT1 (*ut1, seconds since J2000.0 UT1; may be NULL).
 * Errors: MYB_E_RANGE (outside the table or the UTC table, fewer than two rows on either side,
 * missing x, y or UT1 - UTC), MYB_E_NOT_FOUND (no UTC table). Thread-safe. */
MYB_API MYB_NODISCARD myb_status myb_eop_at(const myb_eop* eop, const myb_timesys* system,
                                            myb_time t, myb_scale scale, myb_eop_values* values,
                                            myb_time* ut1);

/* Celestial -> terrestrial rotation at one instant given as tt and ut1, with the EOP (NULL:
 * all zero). dX, dY are added to the X, Y of the series. The derivative is with respect to TT
 * and includes the rates of X, Y, s and of the angle, ERA rate * (1 - LOD / 86400 s). */
MYB_API void myb_cip_terrestrial(const myb_cip* cip, myb_time tt, myb_time ut1,
                                 const myb_eop_values* eop, myb_rotation* rotation);

/* ---- Equinox-based systems ---------------------------------------------------------------
 * The classical (FK5-era) systems of the Earth, from the mean equator and equinox of J2000.0
 * (the NAIF frame J2000, which agrees with the GCRS to the 23 mas frame bias):
 *   MOD  = P [J2000],  P = R3(-z_A) R2(theta_A) R3(-zeta_A)     (IAU 1976 precession)
 *   TOD  = N [MOD],    N = R1(-eps_A - Delta eps) R3(-Delta psi) R1(eps_A)   (IAU 1980 nutation)
 *   TEME: z = the pole of TOD, x = the mean equinox of date projected on the true equator
 *         (the frame of the SGP4 theory, as SPICE defines it for SPK type 10; it differs from
 *         R3(Delta psi cos eps_A) [TOD] by the second-order Delta psi Delta eps sin eps, ~1e-9 rad);
 *   PEF  = R3(GMST82) [TEME]   (pseudo-Earth-fixed: then polar motion gives the ITRS).
 * Every coefficient comes from a kernel (data/frames/iau1980.tk, prefix "IAU1980"). Times are
 * TT (TDB alike) except UT1 for the sidereal time. Derivatives are per second. */

typedef struct myb_equinox myb_equinox;

/* Reads the model <prefix>_* (format in data/frames/iau1980.tk); the pool may be destroyed
 * afterwards. Errors: MYB_E_NOT_FOUND, MYB_E_FORMAT (diag names the variable), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_equinox_create(const myb_kpool* pool, const char* prefix,
                                                    const myb_allocator* allocator,
                                                    myb_equinox** equinox, myb_diag* diag);
MYB_API void myb_equinox_destroy(myb_equinox* equinox); /* accepts NULL */

/* Precession angles zeta_A, z_A, theta_A (rad) and their rates (may be NULL). */
MYB_API void myb_equinox_precession(const myb_equinox* equinox, myb_time tt, double angles[3],
                                    double rates[3]);

/* Mean obliquity eps_A (rad) and its rate (may be NULL). */
MYB_API double myb_equinox_obliquity(const myb_equinox* equinox, myb_time tt, double* rate);

/* Nutation (Delta psi, Delta eps) in rad and their rates (may be NULL). */
MYB_API void myb_equinox_nutation(const myb_equinox* equinox, myb_time tt, double nutation[2],
                                  double rates[2]);

typedef enum myb_equinox_frame {
  MYB_EQUINOX_MOD = 0,
  MYB_EQUINOX_TOD = 1,
  MYB_EQUINOX_TEME = 2
} myb_equinox_frame;

/* The rotation J2000 -> frame at tt, with its derivative. */
MYB_API void myb_equinox_rotation(const myb_equinox* equinox, myb_equinox_frame frame,
                                  myb_time tt, myb_rotation* rotation);

/* Greenwich mean sidereal time of 1982 (rad, [0, 2 pi)) and its rate per UT1 second (may be
 * NULL), at ut1 (an instant on the UT1 scale). */
MYB_API double myb_equinox_gmst(const myb_equinox* equinox, myb_time ut1, double* rate);

/* The rotation TEME -> PEF = R3(GMST82) at ut1, with its derivative (the Earth's rotation). */
MYB_API void myb_equinox_teme_pef(const myb_equinox* equinox, myb_time ut1,
                                  myb_rotation* rotation);

/* ---- Rotation models of the IAU Working Group on Cartographic Coordinates ------------------
 * The orientation of a body as published by the WGCCRE (Archinal et al., Celest. Mech. Dyn.
 * Astron. 130, 22, 2018) in NAIF text PCK form (pck00011.tpc): with T the TDB Julian
 * centuries and d the days from the epoch (J2000.0 unless BODYn_CONSTANTS_JED_EPOCH),
 *   alpha = a0 + a1 T + a2 T^2 + sum_i A_i sin(theta_i)       BODYn_POLE_RA, _NUT_PREC_RA
 *   delta = d0 + d1 T + d2 T^2 + sum_i D_i cos(theta_i)       BODYn_POLE_DEC, _NUT_PREC_DEC
 *   W     = w0 + w1 d + w2 d^2 + sum_i W_i sin(theta_i)       BODYn_PM, _NUT_PREC_PM
 * (degrees), the phase angles theta_i polynomials in T of degree BODYb_MAX_PHASE_DEGREE
 * (default 1) from BODYb_NUT_PREC_ANGLES, b = n / 100 for 0 < n < 1000 (the system's
 * barycentre, NAIF's convention) and b = n otherwise. The rotation from the reference frame
 * (BODYn_CONSTANTS_REF_FRAME, default 1 = J2000) to the body-fixed frame is
 *   M = R3(W) R1(90 deg - delta) R3(90 deg + alpha),
 * Euler angles (90 deg + alpha, 90 deg - delta, W) in the form of myb_rotation_euler313. */

typedef struct myb_wgccre myb_wgccre;

/* The model of body `body` from the kernel variables. Errors: MYB_E_NOT_FOUND (no
 * BODYn_POLE_RA, or phase angles missing while terms need them), MYB_E_FORMAT (wrong counts,
 * more terms than phase angles), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_wgccre_create(const myb_kpool* pool, int64_t body,
                                                   const myb_allocator* allocator,
                                                   myb_wgccre** model, myb_diag* diag);
MYB_API void myb_wgccre_destroy(myb_wgccre* model); /* accepts NULL */

/* NAIF code of the reference frame of the model (1 = J2000/ICRF). */
MYB_API int32_t myb_wgccre_frame(const myb_wgccre* model);

/* Euler angles (rad) and rates (rad/s) at the TDB instant t. Thread-safe. */
MYB_API void myb_wgccre_euler(const myb_wgccre* model, myb_time t, double angles[3],
                              double rates[3]);

/* The rotation reference frame -> body-fixed frame at the TDB instant t. */
MYB_API void myb_wgccre_rotation(const myb_wgccre* model, myb_time t, myb_rotation* rotation);

/* ---- Geodetic coordinates and topocentric frames ------------------------------------------
 * On an oblate ellipsoid of revolution with equatorial radius a and flattening f (from data,
 * e.g. WGS-84 a = 6378137 m, f = 1/298.257223563). Cartesian -> geodetic is the closed form of
 * H. Vermeille, J. Geodesy 85, 105 (2011), valid everywhere including the interior; no
 * iteration. */

typedef struct myb_ellipsoid {
  double a; /* equatorial radius (m), > 0 */
  double f; /* flattening, 0 <= f < 1 */
} myb_ellipsoid;

/* Longitude in (-pi, pi], geodetic latitude in [-pi/2, pi/2] (rad) and height (m) of the
 * body-fixed point xyz (m). The centre gives the north pole at height -b. Errors: MYB_E_DOMAIN
 * (invalid ellipsoid, non-finite input). */
MYB_API MYB_NODISCARD myb_status myb_geodetic_from_cartesian(const myb_ellipsoid* ellipsoid,
                                                             const double xyz[3],
                                                             double* longitude,
                                                             double* latitude, double* height);

/* Body-fixed point of geodetic longitude, latitude (rad) and height (m). Errors: MYB_E_DOMAIN
 * (invalid ellipsoid or latitude, non-finite input). */
MYB_API MYB_NODISCARD myb_status myb_geodetic_to_cartesian(const myb_ellipsoid* ellipsoid,
                                                           double longitude, double latitude,
                                                           double height, double xyz[3]);

/* Rotation body-fixed -> local East-North-Up at geodetic longitude and latitude (constant:
 * dM/dt = 0). */
MYB_API void myb_topocentric_enu(double longitude, double latitude, myb_rotation* rotation);

/* ---- Orbital frames ---------------------------------------------------------------------- */

typedef enum myb_orbital_frame {
  MYB_ORBITAL_RTN = 1, /* radial, transverse (in-track), normal (= RSW, RIC) */
  MYB_ORBITAL_VNC = 2, /* velocity, normal (orbit pole), co-normal = V x N */
  MYB_ORBITAL_LVLH = 3 /* CCSDS local vertical local horizontal: z nadir, y = -pole, x = y x z */
} myb_orbital_frame;

/* Rotation from the frame of `state` (r, v) to the orbital frame, with the derivative that
 * follows from the acceleration (the rates of the axes need it; a central force gives a
 * constant pole). Errors: MYB_E_DOMAIN (r or v zero, r parallel to v, non-finite input). */
MYB_API MYB_NODISCARD myb_status myb_orbital_rotation(myb_orbital_frame kind,
                                                      const double state[6],
                                                      const double acceleration[3],
                                                      myb_rotation* rotation);

/* ---- Fixed frames of NAIF frame kernels ---------------------------------------------------
 * "TK frames" (NAIF Frames Required Reading): a constant rotation of frame `id` relative to
 * the frame named by TKFRAME_<id>_RELATIVE, given as
 *   TKFRAME_<id>_SPEC = 'ANGLES'      _ANGLES (a1 a2 a3), _AXES (i1 i2 i3), _UNITS
 *                                     (TK -> relative = [a1]_i1 [a2]_i2 [a3]_i3)
 *                       'MATRIX'      _MATRIX (TK -> relative, by columns)
 *                       'QUATERNION'  _Q (q0 q1 q2 q3, SPICE style: q0 scalar; TK -> relative).
 * data/frames/naif-inertial.tk defines this way the inertial frames NAIF builds into SPICE that
 * binary PCK and SPK files refer to (B1950, FK4, GALACTIC, ECLIPJ2000, ECLIPB1950). */

/* The rotation relative frame -> frame `id` (constant, dM/dt = 0) and the relative frame's
 * name (NUL-terminated, at most 32 characters). Errors: MYB_E_NOT_FOUND (no definition),
 * MYB_E_FORMAT (unknown SPEC or UNITS, wrong counts, a matrix that is not a rotation). */
MYB_API MYB_NODISCARD myb_status myb_tkframe(const myb_kpool* pool, int64_t id,
                                             myb_rotation* rotation, char relative[33],
                                             myb_diag* diag);

/* ---- Binary PCK -------------------------------------------------------------------------- */

typedef struct myb_pck myb_pck;

MYB_API MYB_NODISCARD myb_status myb_pck_create(const myb_allocator* allocator, myb_pck** pck);
MYB_API void myb_pck_destroy(myb_pck* pck); /* accepts NULL */

/* As myb_spk_load_file / myb_spk_load_memory (myb/kin.h), for binary PCK files. */
MYB_API MYB_NODISCARD myb_status myb_pck_load_file(myb_pck* pck, const char* path, myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_pck_load_memory(myb_pck* pck, const void* data, size_t size,
                                                     myb_diag* diag);

/* Euler angles (rad) and their rates (rad/s) of frame class `id` at t (seconds past J2000 on
 * the file's time scale, TDB for NAIF files) from the highest-priority covering segment, and
 * its base frame code (*base may be NULL). Errors: MYB_E_RANGE (no segment covers t),
 * MYB_E_FORMAT (unsupported type, invalid data). Thread-safe. */
MYB_API MYB_NODISCARD myb_status myb_pck_euler(const myb_pck* pck, int32_t id, myb_time t,
                                               double angles[3], double rates[3], int32_t* base);

/* The rotation from the base frame to frame class `id` at t; errors as myb_pck_euler. */
MYB_API MYB_NODISCARD myb_status myb_pck_rotation(const myb_pck* pck, int32_t id, myb_time t,
                                                  myb_rotation* rotation, int32_t* base);

/* ---- Orientation providers, part 1: see below ------------------------------------------------ */
typedef struct myb_orientation myb_orientation;
typedef myb_status (*myb_orientation_function)(void* context, myb_time t,
                                               myb_orientation* orientation);

/* ---- Writing binary PCK files -----------------------------------------------------------------
 * As the SPK writer of myb/kin.h: segments are collected in memory and laid out as one
 * little-endian file, which SPICE and myb_pck_load_file read. */

typedef struct myb_pck_writer myb_pck_writer;

/* internal_name: the label of the file, at most 60 printable ASCII characters (NULL: "myb").
 * Errors: MYB_E_DOMAIN, MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_pck_writer_create(const char* internal_name,
                                                       const myb_allocator* allocator,
                                                       myb_pck_writer** writer);
MYB_API void myb_pck_writer_destroy(myb_pck_writer* writer); /* accepts NULL */

/* A segment (type 2) of Chebyshev polynomials of the given degree (1 .. 27) on `intervals`
 * equal parts of [start, stop] (seconds past J2000, TDB) for the Euler angles of the frame of
 * class `frame` (e.g. 3000 for ITRF93, 31008 for the lunar principal axes) relative to the
 * inertial frame `base` (NAIF code: 1 = J2000, 17 = ECLIPJ2000): the polynomials through the
 * angles of the provider's rotation at the degree + 1 Chebyshev nodes of each part; its ut,
 * xp, yp are not used. The angles (myb_rotation_euler313) are followed from node to node:
 * the frame must turn by less than half a revolution between two nodes, and its pole must
 * stay away from the z axis of the base frame, where the angles are not defined (for the
 * Earth the ecliptic frame 17 is the usual base). name: 1 .. 40 printable ASCII characters.
 * *max_error (may be NULL; asking costs degree + 2 more calls a part): the largest angle
 * (rad) between the rotation of the provider and that of the polynomials half-way between
 * the nodes. Errors: MYB_E_DOMAIN (the arguments, a rotation that is not finite, the pole on
 * the z axis), MYB_E_NOMEM, MYB_E_OVERFLOW, the statuses of the provider. */
MYB_API MYB_NODISCARD myb_status myb_pck_writer_add_chebyshev(
    myb_pck_writer* writer, int32_t frame, int32_t base, const char* name,
    myb_orientation_function orientation, void* context, double start, double stop,
    size_t intervals, int degree, double* max_error);

/* As myb_spk_writer_add_raw: a segment of any type from its words. */
MYB_API MYB_NODISCARD myb_status myb_pck_writer_add_raw(myb_pck_writer* writer, int32_t frame,
                                                        int32_t base, const char* name,
                                                        int type, double start, double stop,
                                                        const double* words, size_t count);

/* As myb_spk_writer_size, _write and _write_file. */
MYB_API size_t myb_pck_writer_size(const myb_pck_writer* writer);
MYB_API MYB_NODISCARD myb_status myb_pck_writer_write(const myb_pck_writer* writer, void* buffer,
                                                      size_t size);
MYB_API MYB_NODISCARD myb_status myb_pck_writer_write_file(const myb_pck_writer* writer,
                                                           const char* path);

/* ---- Orientation providers -----------------------------------------------------------------
 * The orientation of a body as a function of time, in the form the environment of a force
 * model takes it (myb/env.h): the rotation from inertial axes to the body-fixed frame and, for
 * a body whose rotation is observed (the Earth), what its tides need — the body's own rotation
 * time (UT1) and the coordinates of its pole (polar motion). Providers without those give the
 * argument t and zeros. The time scale of t is the provider's: named where it is registered. */

struct myb_orientation {
  myb_rotation rotation; /* inertial axes -> body-fixed frame */
  myb_time ut;           /* rotation time of the body (UT1 for the Earth); t if it has none */
  double xp;             /* pole coordinates (rad); 0 if not modelled */
  double yp;
};

/* myb_orientation_function (declared above): any status other than MYB_OK is returned to the
 * caller unchanged. */

/* The CIO-based model with observed EOP (myb_eop_at, then myb_cip_terrestrial): t is TT, given
 * as the instant on the scale `tt` of `timesys`; the rotation is celestial -> terrestrial, ut is
 * UT1, xp and yp the polar motion. Errors: those of myb_eop_at. */
typedef struct myb_cip_orientation {
  const myb_cip* cip;
  const myb_eop* eop;
  const myb_timesys* timesys;
  myb_scale tt;
} myb_cip_orientation;

/* myb_orientation_function; context: const myb_cip_orientation*. */
MYB_API MYB_NODISCARD myb_status myb_orientation_cip(void* context, myb_time t,
                                                     myb_orientation* orientation);

/* A frame of a binary PCK (myb_pck_rotation; t on the file's scale, TDB for NAIF files). The
 * segments give the frame relative to their base frame: `base` is the rotation from the
 * inertial axes in use to that base (e.g. myb_tkframe of ECLIPJ2000 for ITRF93), NULL when the
 * two coincide. Errors: those of myb_pck_rotation. */
typedef struct myb_pck_orientation {
  const myb_pck* pck;
  int32_t id;
  const myb_rotation* base;
} myb_pck_orientation;

/* myb_orientation_function; context: const myb_pck_orientation*. */
MYB_API MYB_NODISCARD myb_status myb_orientation_pck(void* context, myb_time t,
                                                     myb_orientation* orientation);

/* A rotation model of the WGCCRE (myb_wgccre_rotation; t is TDB); base as above, relative to
 * the model's reference frame (myb_wgccre_frame). */
typedef struct myb_wgccre_orientation {
  const myb_wgccre* model;
  const myb_rotation* base;
} myb_wgccre_orientation;

/* myb_orientation_function; context: const myb_wgccre_orientation*. */
MYB_API MYB_NODISCARD myb_status myb_orientation_wgccre(void* context, myb_time t,
                                                        myb_orientation* orientation);

/* A constant orientation. myb_orientation_function; context: const myb_rotation*. */
MYB_API MYB_NODISCARD myb_status myb_orientation_fixed(void* rotation, myb_time t,
                                                       myb_orientation* orientation);

#ifdef __cplusplus
}
#endif

#endif /* MYB_FRAMES_H_ */
