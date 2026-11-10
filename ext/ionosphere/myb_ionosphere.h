/* myb_ionosphere.h — the Earth's ionosphere for radio signals as an optional module of myb (target
 * myb_ionosphere; see NOTICE): global maps of the vertical total electron content in the IONEX
 * format (S. Schaer, W. Gurtner, J. Feltens, "IONEX: The IONosphere Map EXchange Format
 * Version 1", 1998; version 1.1, 2015) — the products of the IGS analysis centres (CODE, JPL,
 * ESA, UPC, ...) — and the single-layer (thin shell) model that turns them into the slant
 * content of a ray.
 *
 * The slant content N = M(z) VTEC(pierce point) gives the delay with myb_plasma_group_delay
 * (myb/observe.h): 40.3 N / f^2 m, about 16 cm per TECU (10^16 electrons/m^2) at L1. Angles in
 * radians, heights and radii in metres, time as a modified Julian date in UT (the maps are
 * hourly to two-hourly; a microsecond is far below their resolution). */
#ifndef MYB_IONOSPHERE_H_
#define MYB_IONOSPHERE_H_

#include <stddef.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Single-layer model ----------------------------------------------------------------------
 * All electrons in a shell of height H above a sphere of radius R. For a station at geodetic
 * latitude and longitude and a ray at azimuth (from north towards east) and elevation e, with
 * z = pi/2 - e:
 *   sin z' = R / (R + H) sin(alpha z),  M = 1 / cos z'           (the mapping function),
 * and the ray pierces the shell at the central angle psi = z - asin(R / (R + H) sin z) from the
 * station. alpha = 1 is the standard model (with R = 6371 km, H = 450 km of the IGS maps; the
 * values of a file: myb_ionex_shell); the modified model of CODE uses H = 506.7 km,
 * alpha = 0.9782. */
MYB_EXT_API MYB_NODISCARD myb_status myb_iono_pierce_point(double latitude, double longitude,
                                                       double azimuth, double elevation,
                                                       double radius, double height,
                                                       double alpha, double* pierce_latitude,
                                                       double* pierce_longitude,
                                                       double* mapping);

/* ---- IONEX maps ---------------------------------------------------------------------------- */

typedef struct myb_ionex myb_ionex;

/* Reads the two-dimensional TEC maps of an IONEX file (text of `size` bytes; RMS and height maps
 * are skipped). Errors: MYB_E_FORMAT (with the line in diag), MYB_E_NOMEM, MYB_E_DOMAIN
 * (three-dimensional maps). */
MYB_EXT_API MYB_NODISCARD myb_status myb_ionex_load(const char* text, size_t size,
                                                const myb_allocator* allocator,
                                                myb_ionex** maps, myb_diag* diag);
MYB_EXT_API void myb_ionex_destroy(myb_ionex* maps); /* NULL is allowed */

/* Base radius and shell height of the maps (m). */
MYB_EXT_API void myb_ionex_shell(const myb_ionex* maps, double* radius, double* height);

/* Vertical total electron content (electrons/m^2) at the modified Julian date (UT), geodetic
 * latitude and longitude of the shell: bilinear in latitude and longitude, linear in time
 * between the maps rotated with the Sun (the recommended interpolation of the IONEX
 * description, eq. 3: E_i at longitude + (t - T_i) 360 deg / day). Errors: MYB_E_RANGE (time
 * outside the maps, latitude outside the grid, a missing value of the map). */
MYB_EXT_API MYB_NODISCARD myb_status myb_ionex_vtec(const myb_ionex* maps, double mjd,
                                                double latitude, double longitude, double* vtec);

#ifdef __cplusplus
}
#endif

#endif /* MYB_IONOSPHERE_H_ */
