/* myb_sgp4.h — the SGP4/SDP4 theory of NORAD two-line element sets (TLE) as an optional module of
 * myb (target myb_sgp4; see NOTICE).
 *
 * TLE are mean elements of the US Space Force catalogue, meaningful only with the theory that
 * produced them: F. R. Hoots, R. L. Roehrich, Spacetrack Report #3 (1980), in the revision of
 * D. A. Vallado, P. Crawford, R. Hujsak, T. S. Kelso, "Revisiting Spacetrack Report #3",
 * AIAA 2006-6753, and its code of 2020. The theory is specific to the Earth (WGS-72 constants,
 * lunar and solar terms with fixed coefficients), hence outside the library proper
 * (docs/PLAN.md §3). Its accuracy is that of the elements, about a kilometre at epoch growing
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
MYB_API MYB_NODISCARD myb_status myb_tle_parse(const char* line1, const char* line2,
                                               unsigned flags, myb_tle* tle, myb_diag* diag);

/* The epoch as a UTC calendar date (the fraction of the day converted to seconds with one
 * rounding; the theory itself counts days since 1949 December 31 0h UTC). */
MYB_API void myb_tle_epoch(const myb_tle* tle, myb_calendar* utc);

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
MYB_API void myb_sgp4_constants_of(myb_sgp4_gravity gravity, myb_sgp4_constants* constants);

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
MYB_API MYB_NODISCARD myb_status myb_sgp4_init(const myb_tle* tle, myb_sgp4_gravity gravity,
                                               myb_sgp4_mode mode, myb_sgp4* sgp4);

/* The same from elements in the units of the theory (epoch: days since 1949 December 31 0h
 * UTC; angles rad, Kozai mean motion rad/min; the derivatives of the mean motion in a TLE are
 * not used by SGP4). */
MYB_API MYB_NODISCARD myb_status myb_sgp4_init_elements(
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
MYB_API MYB_NODISCARD myb_status myb_sgp4_propagate(const myb_sgp4* sgp4, double minutes,
                                                    double r[3], double v[3],
                                                    myb_sgp4_error* error);

/* The same in the units of the theory (km, km/s), as printed by its verification program. */
MYB_API MYB_NODISCARD myb_status myb_sgp4_propagate_km(const myb_sgp4* sgp4, double minutes,
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
MYB_API MYB_NODISCARD myb_status myb_sgp4_spk10_state(const myb_spk* spk, size_t index,
                                                      const myb_timesys* system,
                                                      const myb_equinox* equinox, myb_time et,
                                                      double state[6]);

/* Greenwich mean sidereal time of IAU 1982 (rad, [0, 2 pi)) at the UT1 Julian date jd1 + jd2,
 * as the theory computes it (gstime). */
MYB_API double myb_sgp4_gmst(double jd1, double jd2);

#ifdef __cplusplus
}
#endif

#endif /* MYB_SGP4_H_ */
