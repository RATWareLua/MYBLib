/* myb_sgp4.h — the SGP4/SDP4 theory of NORAD two-line element sets (TLE) as an optional module of
 * myb (target myb_sgp4; see NOTICE).
 *
 * TLE are mean elements of the US Space Force catalogue, meaningful only with the theory that
 * produced them: F. R. Hoots, R. L. Roehrich, Spacetrack Report #3 (1980), in the revision of
 * D. A. Vallado, P. Crawford, R. Hujsak, T. S. Kelso, "Revisiting Spacetrack Report #3",
 * AIAA 2006-6753, and its code of 2020. The theory is specific to the Earth (WGS-72 constants,
 * lunar and solar terms with fixed coefficients), hence outside the library proper
 * (docs/PLAN.en.md §3). Its accuracy is that of the elements, about a kilometre at epoch growing
 * by kilometres per day — for precise work fit a numerical orbit to its states.
 *
 * This port follows the Vallado code operation by operation, with correctly rounded elementary
 * functions and an exact fmod: it reproduces the verification output of the code (tcppver.out)
 * to its printed digits on every platform, and propagation is a pure function of the
 * initialised record (the code's deep-space integrator cache is not kept; it never changed the
 * results), so one record serves any number of threads.
 *
 * States are in TEME (true equator, mean equinox of the epoch of the state; see myb/frames.h,
 * myb_equinox_teme and myb_teme_pef), in m and m/s. */
#ifndef MYB_SGP4_H_
#define MYB_SGP4_H_

#include <stdint.h>

#include "myb/ccsds.h"
#include "myb/core.h"
#include "myb/frames.h"
#include "myb/kin.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Two-line element sets ----------------------------------------------------------------- */

typedef struct myb_tle {
  char satnum[6];         /* catalogue number as written (5 characters, Alpha-5 allowed) */
  long catalog_number;    /* the same decoded (Alpha-5: A = 10 ... Z = 33, without I and O) */
  char classification;    /* 'U', 'C', 'S' */
  char designator[9];     /* international designator, trailing blanks removed */
  int epoch_year;         /* four digits (two-digit years 57..99 -> 19xx, 00..56 -> 20xx) */
  double epoch_day;       /* day of the year and fraction, 1.0 = 1 January 0h UTC */
  double ndot;            /* first derivative of the mean motion / 2, rev/day^2 */
  double nddot;           /* second derivative / 6, rev/day^3 */
  double bstar;           /* drag term, 1 / Earth radii */
  int ephemeris_type;
  int element_number;
  double inclination;     /* deg */
  double node;            /* right ascension of the ascending node, deg */
  double eccentricity;
  double argument_of_perigee; /* deg */
  double mean_anomaly;    /* deg */
  double mean_motion;     /* rev/day (Kozai mean motion) */
  long revolution;        /* revolution number at epoch */
} myb_tle;

enum {
  MYB_TLE_NO_CHECKSUM = 1 /* accept lines whose checksum does not match (edited test sets) */
};

/* Parses the two lines (at least 68 characters each; a trailing newline and anything after
 * column 69 are ignored). When column 69 holds a digit it must equal the modulo-10 checksum
 * (digits, '-' = 1) unless flags has MYB_TLE_NO_CHECKSUM. Errors: MYB_E_FORMAT (diag tells the
 * line and column). */
MYB_EXT_API MYB_NODISCARD myb_status myb_tle_parse(const char* line1, const char* line2,
                                               unsigned flags, myb_tle* tle, myb_diag* diag);

/* The epoch as a UTC calendar date (the fraction of the day converted to seconds with one
 * rounding; the theory itself counts days since 1949 December 31 0h UTC). */
MYB_EXT_API void myb_tle_epoch(const myb_tle* tle, myb_calendar* utc);

/* The two lines of an element set, 69 characters each and a NUL, with their checksums: what
 * myb_tle_parse reads back to the same values. The catalogue number is written from satnum
 * when that holds 5 characters, otherwise from catalog_number (Alpha-5 above 99 999). Errors:
 * MYB_E_DOMAIN (a value that does not fit its field: an eccentricity outside [0, 1), an angle
 * outside [0, 1000), a derivative of the mean motion of 1 or more, a drag term whose decimal
 * exponent has two digits, a number above the width of its field, values not finite). */
MYB_EXT_API MYB_NODISCARD myb_status myb_tle_format(const myb_tle* tle, char line1[70],
                                                    char line2[70]);

/* The element set of a mean-elements message of the CCSDS (myb/ccsds.h; the message must give
 * MEAN_MOTION and the TLE-related parameters) and the message of an element set: the numbers
 * pass unchanged, in the units of the catalogue. The epoch of a message is a calendar date,
 * that of an element set a day of the year: the conversion is exact to the rounding of one
 * division. From an element set the message gets the object's identifiers, centre EARTH,
 * frame TEME, time UTC and theory SGP4; the originator, the name and the creation date are
 * the caller's to fill. Errors: MYB_E_DOMAIN (a message without those values, a catalogue
 * number above 339 999). */
MYB_EXT_API MYB_NODISCARD myb_status myb_tle_from_omm(const myb_omm* omm, myb_tle* tle);
MYB_EXT_API MYB_NODISCARD myb_status myb_omm_from_tle(const myb_tle* tle, myb_omm* omm);

/* ---- The theory ---------------------------------------------------------------------------- */

typedef enum myb_sgp4_gravity {
  MYB_SGP4_WGS72_OLD = 0, /* xke = 0.0743669161 as in the original code */
  MYB_SGP4_WGS72 = 1,     /* the constants of the catalogue (default) */
  MYB_SGP4_WGS84 = 2
} myb_sgp4_gravity;

typedef enum myb_sgp4_mode {
  MYB_SGP4_IMPROVED = 0, /* Vallado's corrections ('i') */
  MYB_SGP4_AFSPC = 1     /* the operational code of the Air Force ('a') */
} myb_sgp4_mode;

/* Geophysical constants of the theory. */
typedef struct myb_sgp4_constants {
  double radius; /* equatorial radius of the Earth, km */
  double xke;    /* sqrt(GM) in Earth radii^1.5 per minute */
  double j2, j3, j4;
} myb_sgp4_constants;

/* The constants of a standard model (WGS-72 old: xke = 0.0743669161; WGS-72 and WGS-84:
 * xke = 60 / sqrt(radius^3 / GM) with GM = 398600.8 and 398600.5 km^3/s^2). */
MYB_EXT_API void myb_sgp4_constants_of(myb_sgp4_gravity gravity, myb_sgp4_constants* constants);

/* Initialised record; fields are internal (read-only for the caller). */
typedef struct myb_sgp4 {
  /* constants */
  double tumin, radius, xke, j2, j3, j4, j3oj2;
  myb_sgp4_mode mode;
  int deep;     /* 1: deep space (period >= 225 min) */
  int isimp;
  int irez;     /* resonance: 0 none, 1 one day, 2 half day */
  double epoch; /* days since 1949 December 31 0h UTC */
  /* mean elements */
  double bstar, ecco, argpo, inclo, mo, no_kozai, nodeo, no_unkozai;
  /* near-Earth */
  double aycof, con41, cc1, cc4, cc5, d2, d3, d4, delmo, eta, argpdot, omgcof, sinmao, t2cof,
      t3cof, t4cof, t5cof, x1mth2, x7thm1, mdot, nodedot, xlcof, xmcof, nodecf;
  /* deep space */
  double d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433, dedt, del1, del2,
      del3, didt, dmdt, dnodt, domdt, e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2, sgh3,
      sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, gsto, xfact, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3,
      xl2, xl3, xl4, xlamo, zmol, zmos;
} myb_sgp4;

/* Initialises from a parsed TLE. Errors: MYB_E_DOMAIN (eccentricity outside [0, 1), mean motion
 * <= 0, perigee below the surface). */
MYB_EXT_API MYB_NODISCARD myb_status myb_sgp4_init(const myb_tle* tle, myb_sgp4_gravity gravity,
                                               myb_sgp4_mode mode, myb_sgp4* sgp4);

/* The same from elements in the units of the theory (epoch: days since 1949 December 31 0h
 * UTC; angles rad, Kozai mean motion rad/min; the derivatives of the mean motion in a TLE are
 * not used by SGP4). */
MYB_EXT_API MYB_NODISCARD myb_status myb_sgp4_init_elements(
    const myb_sgp4_constants* constants, myb_sgp4_mode mode, double epoch, double bstar,
    double eccentricity, double argument_of_perigee, double inclination, double mean_anomaly,
    double mean_motion, double node, myb_sgp4* sgp4);

/* Errors of propagation, as the codes of the theory. */
typedef enum myb_sgp4_error {
  MYB_SGP4_OK = 0,
  MYB_SGP4_MEAN_ECCENTRICITY = 1,      /* mean eccentricity outside [-0.001, 1) */
  MYB_SGP4_MEAN_MOTION = 2,            /* mean motion <= 0 */
  MYB_SGP4_PERTURBED_ECCENTRICITY = 3, /* outside [0, 1] */
  MYB_SGP4_SEMILATUS_RECTUM = 4,       /* < 0 */
  MYB_SGP4_DECAYED = 6                 /* radius below the Earth's (state still returned) */
} myb_sgp4_error;

/* TEME position (m) and velocity (m/s) at the given minutes since the epoch of the record.
 * *error (may be NULL) receives the code; codes 1..4 return MYB_E_DOMAIN and leave r, v
 * unchanged, code 6 returns MYB_OK with the state. Thread-safe. */
MYB_EXT_API MYB_NODISCARD myb_status myb_sgp4_propagate(const myb_sgp4* sgp4, double minutes,
                                                    double r[3], double v[3],
                                                    myb_sgp4_error* error);

/* ---- Fitting an element set ------------------------------------------------------------------
 * The mean elements of the theory are not the osculating elements of an orbit: an element set
 * for an object whose motion is known (an ephemeris, a numerical propagation) is the one whose
 * SGP4 orbit passes closest to it. */

typedef struct myb_tle_fit_options {
  int bstar;                /* != 0: the drag term B* is fitted too (the arc must feel it) */
  int max_iterations;       /* 0 = 40 */
  double tolerance;         /* the fit stops when a step lowers the root mean square of the
                               residuals by less than this part of it; 0 = 1e-9 */
  myb_sgp4_gravity gravity; /* the constants to fit for; the catalogue uses MYB_SGP4_WGS72 */
  myb_sgp4_mode mode;
} myb_tle_fit_options;

/* Fits the elements of *tle — mean motion, eccentricity, inclination, node, argument of
 * perigee, mean anomaly and, if asked for, B* — to `count` >= 3 TEME positions (m, 3 each) at
 * the times `minutes` from the epoch of the set, in the sense of least squares: Gauss-Newton
 * steps in unknowns without the singularities of the classical angles, the derivatives by
 * differences of the theory. *tle gives the first guess (myb_tle_from_state makes one), the
 * epoch and what identifies the object, and receives the result; *rms (may be NULL) the root
 * mean square of the residuals of the positions (m), *iterations (may be NULL) the steps
 * taken. options may be NULL (WGS-72, improved mode, no B*). Errors: MYB_E_DOMAIN (fewer than
 * three positions, invalid options, elements that are no orbit of the theory, unknowns that
 * the positions do not determine), MYB_E_NOCONV (the budget of steps; *tle holds the best so
 * far), MYB_E_NOMEM, those of myb_sgp4_propagate. */
MYB_EXT_API MYB_NODISCARD myb_status myb_tle_fit(const double* minutes, const double* positions,
                                             size_t count, const myb_tle_fit_options* options,
                                             myb_tle* tle, double* rms, int* iterations);

/* A first guess for the fit: the osculating elements of the TEME state (m, m/s) at the epoch
 * of the set written as its mean ones — off by the short-period terms of the oblateness, some
 * kilometres. The epoch, B* and the identification of *tle are left as they are. Errors:
 * MYB_E_DOMAIN (not an ellipse), those of myb_elements_from_pv. */
MYB_EXT_API MYB_NODISCARD myb_status myb_tle_from_state(const double r[3], const double v[3],
                                                    myb_sgp4_gravity gravity, myb_tle* tle);

/* The same in the units of the theory (km, km/s), as printed by its verification program. */
MYB_EXT_API MYB_NODISCARD myb_status myb_sgp4_propagate_km(const myb_sgp4* sgp4, double minutes,
                                                       double r[3], double v[3],
                                                       myb_sgp4_error* error);

/* ---- SPK type 10 --------------------------------------------------------------------------
 * Segments of two-line element sets (NAIF SPK Required Reading, type 10; written by SPKW10):
 * eight constants J2, J3, J4, KE, QO, SO, ER, AE (only J2, J3, J4, KE, ER are used, as by
 * SPICE), then for each epoch a packet NDT20, NDD60, BSTAR, INCL, NODE0, ECC, OMEGA, MO, NO
 * (rad, rad/min), EPOCH (TDB seconds past J2000) and four nutation values (unused). The state at
 * et blends the SGP4 states s1, s2 of the two element sets bracketing et (epochs t1 <= et <= t2)
 * with W = (1 + cos(pi (et - t1) / (t2 - t1))) / 2:
 *   p = W p1 + (1 - W) p2,   v = W v1 + (1 - W) v2 + W' (p1 - p2),
 * evaluated in the AFSPC mode with the epoch of each set converted to a UTC Julian date (SPICE
 * rounds it to one double, 40 microseconds; so do we), and rotated from TEME of et to J2000
 * (myb_equinox_rotation; et is used as TT, as SPICE does). Before the first and after the last
 * epoch the nearest set alone is used. */

/* State (m, m/s) of segment `index` of spk, which must be of type 10, at et (TDB), in J2000
 * relative to the segment's centre. system converts TDB to UTC (scale names "TDB"); equinox
 * gives TEME. Errors: MYB_E_FORMAT (not type 10, inconsistent layout), MYB_E_RANGE (et outside
 * the UTC table), statuses of the propagation (MYB_E_DOMAIN). */
MYB_EXT_API MYB_NODISCARD myb_status myb_sgp4_spk10_state(const myb_spk* spk, size_t index,
                                                      const myb_timesys* system,
                                                      const myb_equinox* equinox, myb_time et,
                                                      double state[6]);

/* Writing such a segment: its words for myb_spk_writer_add_raw (myb/kin.h; type 10, frame 1,
 * the centre 399, start and stop the span the segment is to answer for), laid out as SPKW10
 * lays them out — the same words for the same element sets. */

/* The eight constants: J2, J3, J4, KE of the model, QO = 120 and SO = 78 km (the density
 * function), ER the radius in km, AE = 1. */
MYB_EXT_API void myb_sgp4_spk10_constants(const myb_sgp4_constants* model, double constants[8]);

/* The packet of an element set, in the units of the format: NDT20 (rad/min^2), NDD60
 * (rad/min^3), BSTAR, the angles in rad, the mean motion in rad/min, the epoch in TDB seconds
 * past J2000 (system gives the UTC table and the scale "TDB") and the nutation in obliquity
 * and in longitude with their rates at the epoch (IAU 1980; rad, rad/s) from equinox — or
 * zeros if it is NULL: SPICE since N0066 and myb_sgp4_spk10_state do not read them.
 * Errors: those of myb_timesys_from_utc. */
MYB_EXT_API MYB_NODISCARD myb_status myb_sgp4_spk10_packet(const myb_tle* tle,
                                                       const myb_timesys* system,
                                                       const myb_equinox* equinox,
                                                       double packet[14]);

/* The number of words of a segment of `count` element sets (0 if there are none or too
 * many), and the words from the packets (14 each, in the order of their epochs, which must
 * increase). Errors: MYB_E_DOMAIN (no sets, values not finite, epochs not increasing),
 * MYB_E_OVERFLOW (size < myb_sgp4_spk10_size(count)). */
MYB_EXT_API size_t myb_sgp4_spk10_size(size_t count);
MYB_EXT_API MYB_NODISCARD myb_status myb_sgp4_spk10_words(const double constants[8],
                                                      const double* packets, size_t count,
                                                      double* words, size_t size);

/* Greenwich mean sidereal time of IAU 1982 (rad, [0, 2 pi)) at the UT1 Julian date jd1 + jd2,
 * as the theory computes it (gstime). */
MYB_EXT_API double myb_sgp4_gmst(double jd1, double jd2);

#ifdef __cplusplus
}
#endif

#endif /* MYB_SGP4_H_ */
