/* myb_nrlmsise00.h — the NRLMSISE-00 model of the Earth's neutral atmosphere as an optional
 * module of myb (target myb_nrlmsise00; see NOTICE for the status of the material).
 *
 * The model of J. M. Picone, A. E. Hedin, D. P. Drob (J. Geophys. Res. 107, 1468, 2002) is
 * specific to the Earth and not reducible to data, hence outside the library proper (docs/PLAN.en.md
 * §3). This module evaluates it reproducibly (correctly rounded elementary functions) and
 * thread-safely (no shared state). */
#ifndef MYB_NRLMSISE00_H_
#define MYB_NRLMSISE00_H_

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct myb_nrlmsise00_input {
  int day_of_year;          /* 1 .. 366 */
  double seconds;           /* seconds of the day, UT */
  double altitude;          /* geodetic altitude, km */
  double latitude;          /* geodetic latitude, degrees */
  double longitude;         /* degrees */
  double local_solar_time;  /* hours; the model's convention: seconds / 3600 + longitude / 15 */
  double f107_average;      /* 81-day average of F10.7, centred on the day */
  double f107;              /* daily F10.7 of the previous day */
  double ap;                /* daily magnetic index */
  const double* ap_history; /* NULL, or 7 values (see the model's documentation) used when
                               switch 9 is -1: daily Ap, 3-hour ap now, 3, 6, 9 hours before,
                               averages of 8 3-hour values 12-33 and 36-57 hours before */
} myb_nrlmsise00_input;

typedef struct myb_nrlmsise00_output {
  /* number densities of He, O, N2, O2, Ar (index 0 .. 4), H, N, anomalous O (6 .. 8) in m^-3,
   * and the total mass density (index 5) in kg m^-3 */
  double density[9];
  double exospheric_temperature; /* K */
  double temperature;            /* K at the altitude */
} myb_nrlmsise00_output;

/* switches: 24 values as in the model (NULL: all 1 = standard model, ap_history used if given);
 * switch 0 is ignored (SI units always). anomalous_oxygen != 0 (GTD7D): the total mass density
 * includes anomalous oxygen, as appropriate for drag above 500 km. Errors: MYB_E_DOMAIN
 * (non-finite input, altitude < 0 or > 20000 km, day outside 1 .. 366, latitude outside
 * [-90, 90], negative indices, switch 9 = -1 without ap_history). */
MYB_EXT_API MYB_NODISCARD myb_status myb_nrlmsise00(const myb_nrlmsise00_input* input,
                                                const int* switches, int anomalous_oxygen,
                                                myb_nrlmsise00_output* output);

#ifdef __cplusplus
}
#endif

#endif /* MYB_NRLMSISE00_H_ */
