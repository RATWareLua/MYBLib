/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/site.h — the motion of a point of the crust of a body: what its tides and the loads on
 * it add to the place the point has in the body-fixed frame.
 *
 * The formulas are those of the IERS Conventions 2010, chapter 7 (G. Petit, B. Luzum (eds.),
 * IERS Technical Note 36, with the update of 2018), written from the published equations;
 * every number comes from data (data/frames/iers2010-site.tk for the Earth):
 *   the tide of the solid body, step 1: degrees 2 and 3 with nominal Love and Shida numbers
 *     that depend on the latitude (eq. 5, 6), the l(1) terms of the diurnal and semidiurnal
 *     bands (eq. 8, 9) and the out-of-phase parts of the two bands (eq. 10, 11), for every
 *     raising body;
 *   step 2: the frequency dependence of the Love and Shida numbers in the diurnal and the
 *     long-period bands, as Poisson series over chi = GMST + pi and the Delaunay arguments
 *     (eq. 12, 13, Tables 7.3a, 7.3b);
 *   the permanent part of the tide (eq. 14), left in the displacement for a place given as
 *     "conventional tide free" (the places of the ITRF) and taken out for a mean-tide place;
 *   the pole tide: the rotational deformation by the wobble m1 = xp - xs, m2 = -(yp - ys)
 *     about the secular pole (eq. 21, 26);
 *   the load of the ocean tides (section 7.1.2): the amplitudes and phases of a few main
 *     waves at the point (a table in the BLQ form) are carried to every wave of the
 *     tide-generating potential by its admittance — the response A exp(-i phi) / |H| of the
 *     given waves, interpolated in frequency within each band (long-period, diurnal,
 *     semidiurnal) by a cubic spline whose end slopes are those of the parabola through the
 *     three end points (straight lines for three points or fewer), constant outside the
 *     given waves;
 *   the load of the atmospheric tides S1 and S2 (eq. 19) from the twelve coefficients of
 *     the point;
 *   the load of the ocean pole tide (eq. 29) from the six coefficients of the point (a grid
 *     of them: myb_site_grid).
 * A displacement is a vector in the body-fixed frame, metres. The tide of the solid body
 * and the pole tide are along the radius and across it (geocentric latitude); the loads are
 * given along the axes east, north, up of the point (myb_topocentric_enu, myb/frames.h). */
#ifndef MYB_SITE_H_
#define MYB_SITE_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/grav.h"
#include "myb/io.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The numbers of the model. The kernel variables: <prefix>_ARGUMENT_* and the series
 * <prefix>_BAND1_RS, _RC, _TS, _TC, <prefix>_BAND0_R, _T (form of
 * data/frames/iau2006-cip.tk), _RADIUS, _LOVE, _L1, _IMAGINARY, _PERMANENT, _POLE,
 * _SECULAR_POLE, _OCEAN_POLE, _WAVES, _BLQ, _WAVE_SHIFT (described in
 * data/frames/iers2010-site.tk). */
typedef struct myb_site_tides myb_site_tides;

/* Errors: MYB_E_NOT_FOUND, MYB_E_FORMAT (diag names the variable), MYB_E_DOMAIN (prefix longer
 * than 20 characters), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_site_tides_create(const myb_kpool* pool, const char* prefix,
                                                       const myb_allocator* allocator,
                                                       myb_site_tides** tides, myb_diag* diag);
MYB_API void myb_site_tides_destroy(myb_site_tides* tides); /* accepts NULL */

/* The waves of the potential and those a table of ocean loading gives. */
MYB_API size_t myb_site_tides_waves(const myb_site_tides* tides);
MYB_API size_t myb_site_tides_given_waves(const myb_site_tides* tides);

/* ---- The load of the ocean tides at a point ------------------------------------------------- */

enum { MYB_SITE_MAX_GIVEN_WAVES = 16 };

typedef struct myb_ocean_loading myb_ocean_loading;

/* From the amplitudes (m) and the phase lags (degrees) of the given waves, in the order of
 * <prefix>_BLQ: three rows of myb_site_tides_given_waves numbers each — along the radius
 * (up), to the west, to the south. Errors: MYB_E_DOMAIN (a value that is negative or not
 * finite), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_ocean_loading_create(const myb_site_tides* tides,
                                                          const double* amplitudes,
                                                          const double* phases,
                                                          const myb_allocator* allocator,
                                                          myb_ocean_loading** loading);

/* From the record of the point `site` in a file of the BLQ form (the ocean loading service
 * of Onsala; lines beginning with $$ are comments; a record is the name of the point on a
 * line of its own, then three lines of amplitudes and three of phases). The name is compared
 * without regard to case and to the blanks around it. Errors: MYB_E_NOT_FOUND (no such
 * point), MYB_E_FORMAT (diag gives the line), MYB_E_NOMEM, MYB_E_IO. */
MYB_API MYB_NODISCARD myb_status myb_ocean_loading_load(const myb_site_tides* tides,
                                                        const char* text, size_t size,
                                                        const char* site,
                                                        const myb_allocator* allocator,
                                                        myb_ocean_loading** loading,
                                                        myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_ocean_loading_load_file(const myb_site_tides* tides,
                                                             const char* path, const char* site,
                                                             const myb_allocator* allocator,
                                                             myb_ocean_loading** loading,
                                                             myb_diag* diag);
MYB_API void myb_ocean_loading_destroy(myb_ocean_loading* loading); /* accepts NULL */

/* The displacement up, to the west and to the south (m) at the argument time t (UT1) with
 * the difference `dynamical` (both as in myb_site_input). */
MYB_API void myb_ocean_loading_at(const myb_site_tides* tides, const myb_ocean_loading* loading,
                                  myb_time t, double dynamical, double displacement[3]);

/* ---- A grid of coefficients ------------------------------------------------------------------
 * Numbers given at the nodes of a grid of equal steps in longitude and latitude, in the form
 * of the file of the ocean pole tide load of the Conventions (opoleloadcoefcmcor.txt): a
 * header of lines "name = value" with Number_longitude_Grid_Points, First_longitude_degrees,
 * Longitude_step_degrees, Number_latitude_grid_points, First_latitude_degrees,
 * Latitude_step_degrees (names compared without regard to case), a line of dashes, then one
 * row per node: longitude, latitude (degrees) and the values, the longitude running first. */
enum { MYB_SITE_GRID_MAX_VALUES = 12 };

typedef struct myb_site_grid myb_site_grid;

/* Errors: MYB_E_FORMAT (diag gives the line), MYB_E_NOMEM, MYB_E_IO. */
MYB_API MYB_NODISCARD myb_status myb_site_grid_load(const char* text, size_t size,
                                                    const myb_allocator* allocator,
                                                    myb_site_grid** grid, myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_site_grid_load_file(const char* path,
                                                         const myb_allocator* allocator,
                                                         myb_site_grid** grid, myb_diag* diag);
MYB_API void myb_site_grid_destroy(myb_site_grid* grid); /* accepts NULL */
MYB_API size_t myb_site_grid_values(const myb_site_grid* grid); /* per node */

/* The values at a longitude and a latitude (rad), linear between the four nodes around; a
 * grid that goes round the body is continued across its first meridian. Errors: MYB_E_RANGE
 * (outside the grid), MYB_E_DOMAIN (not finite). */
MYB_API MYB_NODISCARD myb_status myb_site_grid_at(const myb_site_grid* grid, double longitude,
                                                  double latitude, double* values);

/* ---- The displacement ------------------------------------------------------------------------ */

/* The loads of a point: what is absent is NULL or zeros. */
typedef struct myb_site_loads {
  const myb_ocean_loading* ocean;
  double atmosphere[12]; /* S1 and S2 (eq. 19): for up, east, north the four numbers A_d1,
                            B_d1, A_d2, B_d2 (m) of cos and sin of one and two cycles per day
                            of UT1 */
  double ocean_pole[6];  /* u_r^R, u_r^I, u_n^R, u_n^I, u_e^R, u_e^I of eq. 29, in the order
                            of the columns of the grid */
} myb_site_loads;

enum {
  MYB_SITE_SOLID = 1,       /* step 1 */
  MYB_SITE_FREQUENCY = 2,   /* step 2 */
  MYB_SITE_POLE = 4,
  MYB_SITE_OCEAN = 8,
  MYB_SITE_ATMOSPHERE = 16,
  MYB_SITE_OCEAN_POLE = 32,
  MYB_SITE_ALL = 63,
  MYB_SITE_MEAN_TIDE = 64   /* the place is a mean-tide one: the permanent part of the tide
                               is taken out of the displacement */
};

typedef struct myb_site_input {
  double gm;                      /* of the body, m^3 s^-2 */
  const myb_tide_raiser* raisers; /* GM and place in the body-fixed frame (myb/grav.h) */
  size_t count;
  myb_time t;                     /* argument time: UT1 (chi = GMST + pi; the day of S1, S2);
                                     also the date of the secular pole */
  double dynamical;               /* TT - UT1 at that instant, s: the arguments of the Moon
                                     and the Sun are taken on TT. With 0 they are taken on
                                     UT1 — 70 s late, 4e-4 rad in the phase of a lunar wave,
                                     a few micrometres */
  double xp;                      /* polar motion, rad */
  double yp;
  int flags;                      /* MYB_SITE_* */
} myb_site_input;

/* The displacement of the point at `position` (body-fixed, m) by what the flags select.
 * `enu` — the axes east, north, up of the point as rows, in the body-fixed frame — and
 * `loads` are needed by the loads only and may be NULL without them. Errors: MYB_E_DOMAIN (a
 * point or a raiser at the centre, a point on the axis of rotation — its longitude is not
 * defined —, gm that is not positive, input that is not finite, loads without the axes). */
MYB_API MYB_NODISCARD myb_status myb_site_displacement(const myb_site_tides* tides,
                                                       const double position[3],
                                                       const double* enu,
                                                       const myb_site_loads* loads,
                                                       const myb_site_input* input,
                                                       double displacement[3]);

#ifdef __cplusplus
}
#endif

#endif /* MYB_SITE_H_ */
