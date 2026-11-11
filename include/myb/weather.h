/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/weather.h — space weather: the indices of the activity of the Sun and of the magnetic
 * field of the Earth that the models of the upper atmosphere take.
 *
 * The indices are observations, kept in files that grow every day: the solar radio flux at
 * 10.7 cm, F10.7 (in solar flux units, 1e-22 W m^-2 Hz^-1: as observed at the distance of
 * the Earth and adjusted to 1 au), its averages over 81 days, the planetary indices Kp and
 * ap of every three hours and their daily mean Ap. The reader takes a table of days as text
 * with a first line of column names (the CSV files of CelesTrak, SW-All.csv and
 * SW-Last5Years.csv: T. S. Kelso, "Space Weather Data Documentation"), finds its columns by
 * their names and keeps what it knows:
 *   DATE (YYYY-MM-DD, the UTC day),
 *   KP1 .. KP8 (Kp of the eight three-hour spans, in tenths), AP1 .. AP8, AP_AVG,
 *   F10.7_OBS, F10.7_ADJ, F10.7_DATA_TYPE (OBS observed, INT interpolated, PRD and PRM
 *   predicted), F10.7_OBS_CENTER81, F10.7_OBS_LAST81, F10.7_ADJ_CENTER81,
 *   F10.7_ADJ_LAST81.
 * Only DATE is required; a column that is absent, or a field left empty (as in the rows of
 * predictions), is a value that is missing. */
#ifndef MYB_WEATHER_H_
#define MYB_WEATHER_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"
#include "myb/frames.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct myb_space_weather myb_space_weather;

/* Errors: MYB_E_FORMAT (diag gives the line: no DATE column, a date or a number that is not
 * one, days out of order), MYB_E_NOMEM, MYB_E_OVERFLOW; MYB_E_IO (files). */
MYB_API MYB_NODISCARD myb_status myb_space_weather_load(const char* text, size_t size,
                                                        const myb_allocator* allocator,
                                                        myb_space_weather** weather,
                                                        myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_space_weather_load_file(const char* path,
                                                             const myb_allocator* allocator,
                                                             myb_space_weather** weather,
                                                             myb_diag* diag);
MYB_API void myb_space_weather_destroy(myb_space_weather* weather); /* accepts NULL */

/* The number of days, and the Modified Julian Dates of the first and the last (either may
 * be NULL). */
MYB_API size_t myb_space_weather_days(const myb_space_weather* weather, int64_t* first,
                                      int64_t* last);

enum { MYB_SPACE_WEATHER_OBSERVED = 0, MYB_SPACE_WEATHER_INTERPOLATED = 1,
       MYB_SPACE_WEATHER_PREDICTED = 2 };

/* A day of the table; a value that the file does not give is a NaN. */
typedef struct myb_space_weather_day {
  int64_t mjd;                /* of the UTC day */
  double kp[8];               /* 0 .. 9, in thirds */
  double ap[8];
  double ap_daily;
  double f107_observed;       /* at the Earth */
  double f107_adjusted;       /* to 1 au */
  double f107_observed_81[2]; /* the mean of 81 days: centred on the day, ending with it */
  double f107_adjusted_81[2];
  int kind;                   /* of the flux: MYB_SPACE_WEATHER_* */
} myb_space_weather_day;

/* Errors: MYB_E_RANGE (a day outside the table or in a gap of it). */
MYB_API MYB_NODISCARD myb_status myb_space_weather_day_at(const myb_space_weather* weather,
                                                          int64_t mjd,
                                                          myb_space_weather_day* day);

/* What the models of the thermosphere of the MSIS family take at an instant of UTC — the
 * day `mjd` and the seconds gone of it:
 *   f107          the flux of the day before;
 *   f107_average  its mean over 81 days centred on the day;
 *   ap            the daily Ap;
 *   ap_history    the daily Ap; ap of the three hours of the instant, of those 3, 6 and 9
 *                 hours before; the means of the eight values from 12 to 33 and from 36 to
 *                 57 hours before (the seven numbers of A. E. Hedin's models);
 *   kp            Kp of the three hours of the instant.
 * adjusted != 0 takes the flux adjusted to 1 au in place of the observed one (the MSIS
 * models were built on the observed). has_history = 0 where a three-hour value of the 57
 * hours is missing: ap_history then holds the daily Ap seven times. Errors: MYB_E_RANGE (a
 * day outside the table, a flux or a daily Ap that the table does not give), MYB_E_DOMAIN
 * (seconds outside 0 .. 86401). */
typedef struct myb_space_weather_inputs {
  double f107;
  double f107_average;
  double ap;
  double ap_history[7];
  int has_history;
  double kp;
} myb_space_weather_inputs;

MYB_API MYB_NODISCARD myb_status myb_space_weather_inputs_at(const myb_space_weather* weather,
                                                             int64_t mjd, double seconds,
                                                             int adjusted,
                                                             myb_space_weather_inputs* inputs);

/* ---- What a model of the upper atmosphere is asked with --------------------------------------
 * The empirical models of the thermosphere (the MSIS family, the Jacchia-Bowman and DTM
 * models, the wind model HWM) take the day of the year and the time of the day on UTC, the
 * geodetic place, the local solar time and the indices. A model itself is not part of the
 * library unless its licence allows it (ext/nrlmsise00); with this a copy of any other one
 * becomes a myb_density_function or a myb_wind_function (myb/force.h) in a few lines:
 *
 *   myb_status my_density(void* context, myb_time t, const double r[3], double* density) {
 *     myb_thermosphere_point p;
 *     const myb_status status = myb_thermosphere_point_at(context, t, r, &p);
 *     if (status == MYB_OK) {
 *       *density = the_model(p.day_of_year, p.seconds, p.height, p.latitude, ...);
 *     }
 *     return status;
 *   }
 */
typedef struct myb_thermosphere_setting {
  const myb_space_weather* weather; /* NULL: the three indices below, at all times */
  double f107;                      /* daily F10.7 of the day before */
  double f107_average;              /* its mean over 81 days */
  double ap;                        /* daily Ap */
  int adjusted;                     /* the flux of the table: 0 observed, 1 adjusted to 1 au */
  const myb_timesys* timesys;
  myb_scale scale;                  /* of the instants asked */
  myb_ellipsoid ellipsoid;          /* of the body */
} myb_thermosphere_setting;

typedef struct myb_thermosphere_point {
  int64_t mjd;             /* the UTC day */
  int day_of_year;         /* 1 .. 366 */
  double seconds;          /* of the day, UTC; below 86400 (a leap second: the day's end) */
  double height;           /* geodetic, m */
  double latitude;         /* geodetic, rad */
  double longitude;        /* rad, (-pi, pi] */
  double local_solar_time; /* hours: the time of the day plus the longitude */
  myb_space_weather_inputs indices;
} myb_thermosphere_point;

/* The point of the instant t and the body-fixed place r (m). Errors: those of
 * myb_timesys_to_utc, myb_geodetic_from_cartesian and myb_space_weather_inputs_at
 * (MYB_E_RANGE outside the table). */
MYB_API MYB_NODISCARD myb_status myb_thermosphere_point_at(
    const myb_thermosphere_setting* setting, myb_time t, const double r[3],
    myb_thermosphere_point* point);

#ifdef __cplusplus
}
#endif

#endif /* MYB_WEATHER_H_ */
