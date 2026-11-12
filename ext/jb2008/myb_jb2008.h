/* myb_jb2008.h — the Jacchia-Bowman 2008 model of the density of the Earth's thermosphere as an
 * optional module of myb (target myb_jb2008; see NOTICE): a clean implementation from the
 * published papers, not a translation of the program of Space Environment Technologies (whose
 * licence forbids that).
 *
 * Sources: B. R. Bowman, W. K. Tobiska, F. A. Marcos, C. Y. Huang, C. S. Lin, W. J. Burke, "A
 * New Empirical Thermospheric Density Model JB2008 Using New Solar and Geomagnetic Indices",
 * AIAA 2008-6438 (the exospheric temperature from F10, S10, M10, Y10; the semiannual
 * variation); B. R. Bowman et al., "The JB2006 empirical thermospheric density model", J. Atm.
 * Solar-Terr. Phys. 70, 774 (2008) (the corrections of the temperature by local time, latitude
 * and height; the factor above 1000 km); L. G. Jacchia, SAO Special Reports 313 (1970) and 332
 * (1971) (the static model: the temperature profile, the barometric and diffusion equations,
 * the composition at 90 km, hydrogen; the diurnal and the seasonal-latitudinal variations).
 *
 * Where the papers leave a choice open, it was settled by comparing the outputs of the program
 * of SET (run unchanged, as a black box, never read) for the same inputs: the J70 constants of
 * the static model and the diurnal variation (m = 2.5, R = 0.31); the corrections of JB2006
 * equations 7 and 8 in their bands of height (200-240 and 300-600 km), joined by cubics of
 * Hermite to zero at 120 and 800 km and to each other between 240 and 300 km, with the daily
 * F10; the seasonal-latitudinal term of J70 with the time argument of J71; no
 * seasonal-latitudinal variation of helium; the semiannual F(z) not below zero. The exospheric
 * and the local temperatures then agree with that program to 1e-6 and 2e-5 K; the density to
 * 9e-4 in log10 (0.2 %) rms and 4.8e-3 (1.1 %) at worst over 655 inputs from 90 to 6000 km,
 * set by the three significant digits of the semiannual coefficients and the four of the
 * factor above 1000 km as the papers print them, and by the integration (here Gauss-Legendre;
 * the program's own Newton-Cotes steps are not published). The static model reproduces the
 * tables of Jacchia (1970) to their rounding.
 *
 * The indices come from the files of SET (SOLFSMY.TXT, DTCFILE.TXT: obtained under their
 * licence, not distributed with myb), read here by their documented formats. */
#ifndef MYB_JB2008_H_
#define MYB_JB2008_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/frames.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The indices of one instant, already lagged: F10 and S10 of the day before, M10 of two days
 * before, Y10 of five days before, each with its 81-day centred mean of the same day (sfu); and
 * the change of the exospheric temperature by the geomagnetic storm (K). */
typedef struct myb_jb2008_indices {
  double f10;
  double f10_average;
  double s10;
  double s10_average;
  double m10;
  double m10_average;
  double y10;
  double y10_average;
  double dtc;
} myb_jb2008_indices;

typedef struct myb_jb2008_output {
  double density;                /* kg m^-3 */
  double exospheric_temperature; /* K, above the place */
  double temperature;            /* K, at the place */
} myb_jb2008_output;

/* The model at the modified Julian date `mjd` (UT), for the Sun at the right ascension and the
 * declination given (rad) and the place at the right ascension `right_ascension` (rad: the
 * sidereal angle plus the east longitude), the geocentric latitude `latitude` (rad) and the
 * height `height` (km, from 90 km up; checked against the program to 4000 km). Errors:
 * MYB_E_DOMAIN (non-finite input, height below 90 km, negative indices). */
MYB_EXT_API MYB_NODISCARD myb_status myb_jb2008(double mjd, double sun_right_ascension,
                                            double sun_declination, double right_ascension,
                                            double latitude, double height,
                                            const myb_jb2008_indices* indices,
                                            myb_jb2008_output* output);

/* ---- The files of indices ---------------------------------------------------------------------
 * SOLFSMY.TXT: comment lines begin with '#'; a record is "YYYY DDD JulianDay F10 F81c S10 S81c M10
 * M81c Y10 Y81c Ssrc", the values of one UTC day. DTCFILE.TXT: "DTC YYYY DDD" and 24 values of
 * the hours of that UTC day (K). The days must follow one another without gaps. */
typedef struct myb_jb2008_solar myb_jb2008_solar;
typedef struct myb_jb2008_storm myb_jb2008_storm;

/* Errors: MYB_E_FORMAT (diag gives the line), MYB_E_NOMEM, MYB_E_IO. */
MYB_EXT_API MYB_NODISCARD myb_status myb_jb2008_solar_load(const char* text, size_t size,
                                                       const myb_allocator* allocator,
                                                       myb_jb2008_solar** table, myb_diag* diag);
MYB_EXT_API MYB_NODISCARD myb_status myb_jb2008_solar_load_file(const char* path,
                                                            const myb_allocator* allocator,
                                                            myb_jb2008_solar** table,
                                                            myb_diag* diag);
MYB_EXT_API void myb_jb2008_solar_destroy(myb_jb2008_solar* table); /* NULL is allowed */

/* The eight solar indices of an instant on the UTC day `mjd` (the integer part is the day),
 * lagged as the model wants them (the dtc of `indices` is left as it is). Errors: MYB_E_RANGE
 * (a day the table does not hold). */
MYB_EXT_API MYB_NODISCARD myb_status myb_jb2008_solar_at(const myb_jb2008_solar* table,
                                                     double mjd, myb_jb2008_indices* indices);

MYB_EXT_API MYB_NODISCARD myb_status myb_jb2008_storm_load(const char* text, size_t size,
                                                       const myb_allocator* allocator,
                                                       myb_jb2008_storm** table, myb_diag* diag);
MYB_EXT_API MYB_NODISCARD myb_status myb_jb2008_storm_load_file(const char* path,
                                                            const myb_allocator* allocator,
                                                            myb_jb2008_storm** table,
                                                            myb_diag* diag);
MYB_EXT_API void myb_jb2008_storm_destroy(myb_jb2008_storm* table); /* NULL is allowed */

/* The value of DTCFILE for the UTC hour that holds the instant `mjd`. Errors: MYB_E_RANGE. */
MYB_EXT_API MYB_NODISCARD myb_status myb_jb2008_storm_at(const myb_jb2008_storm* table,
                                                     double mjd, double* dtc);

/* ---- The atmosphere of a force model ----------------------------------------------------------
 * myb_density_function (myb/force.h; the context is this structure): the density at the
 * body-fixed place r (m) at the instant t. The UTC day and hour of t (on `scale`) choose the
 * indices of the tables (or, without them, those of `fixed` at all times); the place has its
 * geocentric latitude and its height over the ellipsoid; the Sun is the low-precision Sun of
 * the Astronomical Almanac (0.01 degree) and the right ascension of the place the Greenwich
 * mean sidereal time of 1982 on UTC plus the longitude: the hour angle of the Sun, which is what
 * the model uses, to about a hundredth of a degree. Errors: those of myb_timesys_to_utc,
 * myb_geodetic_from_cartesian, the tables (MYB_E_RANGE) and myb_jb2008. */
typedef struct myb_jb2008_atmosphere {
  const myb_jb2008_solar* solar; /* NULL: fixed */
  const myb_jb2008_storm* storm; /* NULL: fixed.dtc */
  myb_jb2008_indices fixed;
  const myb_timesys* timesys;
  myb_scale scale;          /* of t */
  myb_ellipsoid ellipsoid;  /* of the Earth */
} myb_jb2008_atmosphere;

MYB_EXT_API MYB_NODISCARD myb_status myb_jb2008_density(void* atmosphere, myb_time t,
                                                    const double r[3], double* density);

#ifdef __cplusplus
}
#endif

#endif /* MYB_JB2008_H_ */
