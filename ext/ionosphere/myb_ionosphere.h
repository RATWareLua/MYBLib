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
 * hourly to two-hourly; a microsecond is far below their resolution).
 *
 * Beyond the first order (IERS Conventions 2010, section 9.4, eq. 9.19 to 9.23 and 9.40; M.
 * Fritsche et al., Geophys. Res. Lett. 32, L23311, 2005) the delay of a signal of frequency f
 * is, for the carrier phase and for the code,
 *   -s1 / f^2 - s2 / f^3 - s3 / f^4   and   s1 / f^2 + 2 s2 / f^3 + 3 s3 / f^4,
 *   s1 = K1 int N dl,  s2 = K2 int N B cos(theta) dl,
 *   s3 = K3 int N^2 dl + K4 int N B^2 (1 + cos^2(theta)) dl,
 * N the density of free electrons, B the magnetic field and theta its angle to the direction
 * of propagation, K1 = e^2 / (8 pi^2 m eps0) = 40.308, K2 = e^3 / (16 pi^3 m^2 eps0) =
 * 1.1283e12, K3 = e^4 / (128 pi^4 m^2 eps0^2) = 812.38, K4 = e^4 / (64 pi^4 m^3 eps0) =
 * 1.5792e22 (SI). The integrals along the ray need the density in space: a function of the
 * place, here a layer over the maps or a grid of densities; the field is one of myb/force.h
 * (the IGRF). */
#ifndef MYB_IONOSPHERE_H_
#define MYB_IONOSPHERE_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/force.h"
#include "myb/grav.h"

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

/* ---- The density of electrons in space --------------------------------------------------------
 * The density of free electrons (m^-3) at the modified Julian date (UT) and the place r (m)
 * in the Earth-fixed frame. Any status other than MYB_OK stops the integration and is
 * returned. */
typedef myb_status (*myb_electron_density)(void* context, double mjd, const double r[3],
                                           double* density);

/* A layer over the maps: the vertical content of the maps under the place (at its geocentric
 * latitude and longitude) spread over the height h above the sphere of the maps' base radius
 * as a Chapman layer,
 *   N(h) = N_max exp((1 - z - exp(-z)) / 2),  z = (h - peak) / scale,
 * N_max = VTEC / (scale sqrt(2 pi e)): its integral over the height is the content of the map.
 * Its integral of N^2 is N_max VTEC sqrt(e / (2 pi)): the shape factor 0.66 of the
 * Conventions (eq. 9.47). myb_iono_layer_density is a myb_electron_density (context: const
 * myb_iono_layer*). Errors: those of myb_ionex_vtec, MYB_E_DOMAIN (scale not positive). */
typedef struct myb_iono_layer {
  const myb_ionex* maps;
  double peak;  /* the height of the largest density, m (300 .. 450 km) */
  double scale; /* the scale height, m (60 .. 100 km; the peak several of them up) */
} myb_iono_layer;

MYB_EXT_API MYB_NODISCARD myb_status myb_iono_layer_density(void* layer, double mjd,
                                                        const double r[3], double* density);

/* Densities at the nodes of a grid of equal steps in the geocentric latitude, the longitude
 * (rad), the height above a sphere (m) and at given epochs — three-dimensional maps, in
 * arrays of the caller: linear between the eight nodes around the place and between the two
 * epochs around the date, each rotated with the Sun as the maps of myb_ionex_vtec; nothing
 * below the lowest and above the highest node. myb_iono_grid_density is a
 * myb_electron_density (context: const myb_iono_grid*). Errors: MYB_E_RANGE (the date or the
 * latitude outside the grid), MYB_E_DOMAIN (a grid that is none). */
typedef struct myb_iono_grid {
  double radius;         /* of the sphere the heights are counted from, m */
  double latitude;       /* of the first node, and the step (either sign) */
  double latitude_step;
  int latitudes;         /* >= 2 */
  double longitude;
  double longitude_step; /* > 0; the grid goes round when longitudes * step is a turn */
  int longitudes;        /* >= 2 */
  double height;
  double height_step;    /* > 0 */
  int heights;           /* >= 2 */
  const double* epochs;  /* MJD, increasing */
  int maps;              /* >= 1 */
  const double* density; /* [map][height][latitude][longitude], m^-3 */
} myb_iono_grid;

MYB_EXT_API MYB_NODISCARD myb_status myb_iono_grid_density(void* grid, double mjd,
                                                       const double r[3], double* density);

/* ---- The integrals along a ray and the orders of the delay ---------------------------------- */

typedef struct myb_iono_ray {
  double content;      /* int N dl: the slant total electron content, m^-2 */
  double field;        /* int N B cos(theta) dl, T m^-2 */
  double square;       /* int N^2 dl, m^-5 */
  double field_square; /* int N B^2 (1 + cos^2(theta)) dl, T^2 m^-2 */
} myb_iono_ray;

/* The integrals along the straight line between the receiver and the transmitter (places in
 * the Earth-fixed frame, m), over its parts between the radii `bottom` and `top` of the
 * ionosphere (e.g. the base radius + 50 km and + 2000 km): on each part `panels` equal panels
 * of the four-point Gauss rule (1 .. 4096; panels as wide as the scale height of the layer
 * give 4e-8 of the content, twice as wide 3e-6). theta is the angle of the field to the direction transmitter -> receiver. `field`
 * may be NULL (the two integrals with the field are then zero); `work` as for
 * myb_magnetic_field_eval. Errors: MYB_E_DOMAIN (the ends at one place, radii or panels that
 * are none, non-finite input), the statuses of the density and of the field. */
MYB_EXT_API MYB_NODISCARD myb_status myb_iono_ray_integrals(myb_electron_density density,
                                                        void* context,
                                                        const myb_magnetic_field* field,
                                                        myb_gravity_work* work, double mjd,
                                                        const double receiver[3],
                                                        const double transmitter[3],
                                                        double bottom, double top, int panels,
                                                        myb_iono_ray* ray);

/* The coefficients s1 (m Hz^2), s2 (m Hz^3), s3 (m Hz^4) of the delay from the integrals and
 * the constants e, m_e, eps0 (MYB_CONST_E, MYB_CONST_ME, MYB_CONST_EPS0 of
 * data/constants/codata-2022.tk). */
typedef struct myb_iono_orders {
  double s1;
  double s2;
  double s3;
} myb_iono_orders;

MYB_EXT_API void myb_iono_orders_of(const myb_iono_ray* ray, double e, double m_e, double eps0,
                                    myb_iono_orders* orders);

/* The delay (m) of the carrier phase, -s1 / f^2 - s2 / f^3 - s3 / f^4, and of the code,
 * s1 / f^2 + 2 s2 / f^3 + 3 s3 / f^4, at the frequency f (Hz). The sign of s2 is that of a
 * right-hand circularly polarised signal (the GNSS). NaN for a frequency that is not
 * positive. */
MYB_EXT_API double myb_iono_phase_delay(const myb_iono_orders* orders, double frequency);
MYB_EXT_API double myb_iono_code_delay(const myb_iono_orders* orders, double frequency);

/* What the bending of the ray adds to s3 (m Hz^4) for the slant content (m^-2) and the
 * elevation (rad), by the fit of N. Jakowski et al. (1994) to ray tracing at the frequencies
 * of GPS (eq. 9.37 of the Conventions, which name no bending model conventional):
 * 2.495e8 ((1 - 0.8592 cos^2 E)^(-1/2) - 1) S^2 in mm MHz^4 with S in TECU. NaN for invalid
 * input. */
MYB_EXT_API double myb_iono_bending(double content, double elevation);

#ifdef __cplusplus
}
#endif

#endif /* MYB_IONOSPHERE_H_ */
