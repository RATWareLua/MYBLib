/* SPDX-License-Identifier: BSD-3-Clause */
/* Orbit data messages of the CCSDS (myb/ccsds.h) on messages in the style of the examples of
 * the Blue Book (CCSDS 502.0-B-2, figures 3-1, 3-3, 4-1, 5-1 and annex examples):
 *  - OEM: header and metadata, two segments, both forms of the epoch, accelerations, a
 *    covariance section passed over; the states in SI; written and read back; errors with
 *    their lines;
 *  - OPM: state, elements, spacecraft, manoeuvres, units in brackets; written and read back;
 *  - OMM: mean elements with the TLE-related parameters; written and read back. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/ccsds.h"

static const char oem_text[] =
    "CCSDS_OEM_VERS = 2.0\n"
    "COMMENT a sample in the style of CCSDS 502.0-B-2, figure 5-1\n"
    "CREATION_DATE = 1996-11-04T17:22:31\n"
    "ORIGINATOR = NASA/JPL\n"
    "\n"
    "META_START\n"
    "OBJECT_NAME = MARS GLOBAL SURVEYOR\n"
    "OBJECT_ID = 1996-062A\n"
    "CENTER_NAME = MARS BARYCENTER\n"
    "REF_FRAME = EME2000\n"
    "TIME_SYSTEM = UTC\n"
    "START_TIME = 1996-12-18T12:00:00.331\n"
    "USEABLE_START_TIME = 1996-12-18T12:10:00.331\n"
    "USEABLE_STOP_TIME = 1996-12-28T21:23:00.331\n"
    "STOP_TIME = 1996-12-28T21:28:00.331\n"
    "INTERPOLATION = HERMITE\n"
    "INTERPOLATION_DEGREE = 7\n"
    "META_STOP\n"
    "\n"
    "COMMENT This file was produced by M.R. Somebody\n"
    "1996-12-18T12:00:00.331 2789.619 -280.045 -1746.755 4.73372 -2.49586 -1.04195\n"
    "1996-12-18T12:01:00.331 2783.419 -308.143 -1877.071 5.18604 -2.42124 -1.99608\r\n"
    "  1996-12-18T12:02:00.331   2776.033 -336.859 -2008.682 5.63678 -2.33951 -1.94687\n"
    "\n"
    "META_START\n"
    "OBJECT_NAME = MARS GLOBAL SURVEYOR\n"
    "OBJECT_ID = 1996-062A\n"
    "CENTER_NAME = MARS BARYCENTER\n"
    "REF_FRAME = EME2000\n"
    "TIME_SYSTEM = UTC\n"
    "START_TIME = 1996-363T21:29:00.331Z\n"
    "STOP_TIME = 1996-12-28T21:30:00.331\n"
    "META_STOP\n"
    "1996-363T21:29:00.331 -2432.166 -063.042 1742.754 7.33702 -3.495867 -1.041945 0.001 "
    "-0.002 3e-4\n"
    "1996-12-28T21:30:00.331 -2445.234 -878.141 1873.073 1.86043 -3.421256 -0.996366 0.0011 "
    "-0.0021 3.1E-04\n"
    "\n"
    "COVARIANCE_START\n"
    "EPOCH = 1996-12-28T21:29:07.267\n"
    "COV_REF_FRAME = EME2000\n"
    " 3.3313494e-04\n"
    " 4.6189273e-04  6.7824216e-04\n"
    "-3.0700078e-04 -4.2212341e-04  3.2319319e-04\n"
    "-3.3493650e-07 -4.6860842e-07  2.4849495e-07  4.2960228e-10\n"
    "-2.2118325e-07 -2.8641868e-07  1.7980986e-07  2.6088992e-10  1.7675147e-10\n"
    "-3.0413460e-07 -4.9894969e-07  3.5403109e-07  1.8692631e-10  1.0088625e-10  6.2244443e-10\n"
    "EPOCH = 1996-12-29T21:00:00\n"
    "1\n"
    "0 2\n"
    "0 0 3\n"
    "0 0 0 4e-6\n"
    "0 0 0 0 5e-6\n"
    "0.5 0 0 0 0 6e-6\n"
    "COVARIANCE_STOP\n";

static myb_status load(const char* text, myb_oem** oem, myb_diag* diag) {
  return myb_oem_load(text, strlen(text), NULL, oem, diag);
}

static void test_oem(void) {
  myb_oem* oem = NULL;
  myb_diag diag = {0};
  MT_CHECK(load(oem_text, &oem, &diag) == MYB_OK);
  if (oem == NULL) {
    fprintf(stderr, "line %ld: %s\n", diag.line, diag.message);
    return;
  }
  MT_CHECK(strcmp(myb_oem_header(oem, "CCSDS_OEM_VERS"), "2.0") == 0);
  MT_CHECK(strcmp(myb_oem_header(oem, "ORIGINATOR"), "NASA/JPL") == 0);
  MT_CHECK(myb_oem_header(oem, "OBJECT_NAME") == NULL);
  MT_CHECK(myb_oem_segment_count(oem) == 2);
  MT_CHECK(strcmp(myb_oem_meta(oem, 0, "OBJECT_NAME"), "MARS GLOBAL SURVEYOR") == 0);
  MT_CHECK(strcmp(myb_oem_meta(oem, 0, "INTERPOLATION_DEGREE"), "7") == 0);
  MT_CHECK(myb_oem_meta(oem, 1, "INTERPOLATION") == NULL);
  MT_CHECK(strcmp(myb_oem_meta(oem, 1, "START_TIME"), "1996-363T21:29:00.331Z") == 0);
  MT_CHECK(myb_oem_state_count(oem, 0) == 3 && myb_oem_state_count(oem, 1) == 2);
  myb_calendar epoch;
  double state[6];
  double acceleration[3];
  MT_CHECK(myb_oem_state(oem, 0, 2, &epoch, state, acceleration) == 1);
  MT_CHECK(epoch.year == 1996 && epoch.month == 12 && epoch.day == 18 && epoch.hour == 12 &&
           epoch.minute == 2);
  MT_CHECK_SAME(epoch.second, 0.331);
  MT_CHECK_SAME(state[0], 2776.033 * 1000.0);
  MT_CHECK_SAME(state[2], -2008.682 * 1000.0);
  MT_CHECK_SAME(state[5], -1.94687 * 1000.0);
  MT_CHECK(acceleration[0] == 0.0 && acceleration[2] == 0.0);
  /* the day-of-year form of the epoch; accelerations */
  MT_CHECK(myb_oem_state(oem, 1, 0, &epoch, state, acceleration) == 3);
  MT_CHECK(epoch.year == 1996 && epoch.month == 12 && epoch.day == 28 && epoch.hour == 21 &&
           epoch.minute == 29);
  MT_CHECK_SAME(state[1], -63.042 * 1000.0);
  MT_CHECK_SAME(acceleration[2], 3e-4 * 1000.0);
  MT_CHECK(myb_oem_state(oem, 1, 1, NULL, NULL, acceleration) == 3);
  MT_CHECK_SAME(acceleration[2], 3.1e-4 * 1000.0);

  /* written from what was read, and read again */
  static myb_calendar epochs[2][3];
  static double states[2][18];
  myb_oem_segment_data segments[2];
  for (size_t s = 0; s < 2; ++s) {
    const size_t n = myb_oem_state_count(oem, s);
    for (size_t k = 0; k < n; ++k) {
      (void)myb_oem_state(oem, s, k, &epochs[s][k], states[s] + 6 * k, NULL);
    }
    segments[s] = (myb_oem_segment_data){
        .object_name = myb_oem_meta(oem, s, "OBJECT_NAME"),
        .object_id = myb_oem_meta(oem, s, "OBJECT_ID"),
        .center_name = myb_oem_meta(oem, s, "CENTER_NAME"),
        .ref_frame = myb_oem_meta(oem, s, "REF_FRAME"),
        .time_system = myb_oem_meta(oem, s, "TIME_SYSTEM"),
        .interpolation = s == 0 ? "HERMITE" : NULL,
        .interpolation_degree = 7,
        .count = n,
        .epochs = epochs[s],
        .states = states[s]};
  }
  const myb_oem_data data = {.originator = "MYB",
                             .creation_date = {2026, 11, 10, 9, 8, 7.25},
                             .segment_count = 2,
                             .segments = segments};
  size_t length = 0;
  MT_CHECK(myb_oem_write(&data, NULL, 0, &length) == MYB_E_OVERFLOW);
  MT_CHECK(length > 500);
  char* text = malloc(length + 1);
  char* again = malloc(length + 1);
  MT_CHECK(text != NULL && again != NULL);
  if (text == NULL || again == NULL) {
    return;
  }
  size_t written = 0;
  MT_CHECK(myb_oem_write(&data, text, length, &written) == MYB_E_OVERFLOW); /* no room for NUL */
  MT_CHECK(myb_oem_write(&data, text, length + 1, &written) == MYB_OK);
  MT_CHECK(written == length && strlen(text) == length);
  MT_CHECK(strstr(text, "CREATION_DATE = 2026-11-10T09:08:07.250\n") != NULL);
  MT_CHECK(strstr(text, "START_TIME = 1996-12-28T21:29:00.331\n") != NULL);
  MT_CHECK(strstr(text, "INTERPOLATION_DEGREE = 7\n") != NULL);
  MT_CHECK(strstr(text, "1996-12-18T12:00:00.331 2789.619 -280.045 -1746.755 4.73372 -2.49586 "
                        "-1.04195\n") != NULL);
  myb_oem* back = NULL;
  MT_CHECK(load(text, &back, &diag) == MYB_OK);
  if (back != NULL) {
    MT_CHECK(myb_oem_segment_count(back) == 2);
    MT_CHECK(strcmp(myb_oem_header(back, "ORIGINATOR"), "MYB") == 0);
    static myb_calendar epochs2[2][3];
    static double states2[2][18];
    myb_oem_segment_data segments2[2];
    double worst = 0.0;
    for (size_t s = 0; s < 2; ++s) {
      const size_t n = myb_oem_state_count(back, s);
      MT_CHECK(n == segments[s].count);
      for (size_t k = 0; k < n; ++k) {
        MT_CHECK(myb_oem_state(back, s, k, &epochs2[s][k], states2[s] + 6 * k, NULL) == 1);
        MT_CHECK(memcmp(&epochs2[s][k], &epochs[s][k], sizeof(myb_calendar)) == 0);
        for (int i = 0; i < 6; ++i) {
          const double a = states[s][6 * k + (size_t)i];
          const double b = states2[s][6 * k + (size_t)i];
          worst = fmax(worst, fabs(a - b) / (fabs(a) * 2.220446049250313e-16));
        }
      }
      segments2[s] = segments[s];
      segments2[s].epochs = epochs2[s];
      segments2[s].states = states2[s];
    }
    printf("OEM written and read back: states within %.2f ulp\n", worst);
    MT_CHECK(worst <= 1.0);
    /* and written again: the same text */
    myb_oem_data data2 = data;
    data2.segments = segments2;
    MT_CHECK(myb_oem_write(&data2, again, length + 1, &written) == MYB_OK);
    MT_CHECK(written == length && memcmp(again, text, length) == 0);
    myb_oem_destroy(back);
  }
  free(text);
  free(again);
  /* invalid descriptions */
  myb_oem_data bad = data;
  bad.segment_count = 0;
  MT_CHECK(myb_oem_write(&bad, NULL, 0, &length) == MYB_E_DOMAIN);
  bad = data;
  bad.originator = "";
  MT_CHECK(myb_oem_write(&bad, NULL, 0, &length) == MYB_E_DOMAIN);
  bad = data;
  bad.creation_date.month = 13;
  MT_CHECK(myb_oem_write(&bad, NULL, 0, &length) == MYB_E_DOMAIN);
  myb_oem_destroy(oem);
  myb_oem_destroy(NULL);

  /* damaged messages: the line is reported */
  static const struct {
    const char* text;
    long line;
  } broken[] = {
      {"CREATION_DATE = 1996-11-04T17:22:31\n", 1}, /* does not start with the version */
      {"CCSDS_OEM_VERS = 2.0\nMETA_STOP\n", 2},
      {"CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = A\n1996-12-18T12:00:00 1 2 3 4 5 6\n", 4},
      {"CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = A\nMETA_STOP\n1996-12-18T12:00:00 1 2 3 "
       "4 5\n",
       5},
      {"CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = A\nMETA_STOP\n1996-13-18T12:00:00 1 2 3 "
       "4 5 6\n",
       5},
      {"CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = A\nMETA_STOP\n1996-12-18T12:00:00 1 2 3 "
       "4 5 six\n",
       5},
      {"CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = A\nMETA_STOP\nKEY = 1\n", 5},
      {"CCSDS_OEM_VERS = 2.0\nORIGINATOR = X\n", 2},                           /* no segment */
      {"CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = A\n", 3},              /* not closed */
      {"CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = A\nMETA_STOP\n", 0},   /* no OBJECT_ID */
  };
  for (size_t k = 0; k < sizeof broken / sizeof broken[0]; ++k) {
    myb_oem* none = NULL;
    MT_CHECK(load(broken[k].text, &none, &diag) == MYB_E_FORMAT);
    MT_CHECK(none == NULL && diag.message[0] != '\0');
    MT_CHECK(diag.line == broken[k].line);
  }
  myb_oem* none = NULL;
  MT_CHECK(myb_oem_load_file("no-such-file.oem", NULL, &none, &diag) == MYB_E_IO);
}

static const char opm_text[] =
    "CCSDS_OPM_VERS = 2.0\n"
    "CREATION_DATE = 2000-06-03T05:33:00.123\n"
    "ORIGINATOR = GSOC\n"
    "COMMENT GEOCENTRIC, CARTESIAN, EARTH FIXED\n"
    "OBJECT_NAME = EUTELSAT W4\n"
    "OBJECT_ID = 2000-028A\n"
    "CENTER_NAME = EARTH\n"
    "REF_FRAME = TOD\n"
    "TIME_SYSTEM = UTC\n"
    "EPOCH = 2006-06-03T00:00:00.000\n"
    "X = 6655.9942 [km]\n"
    "Y = -40218.5751 [km]\n"
    "Z = -82.9177 [km]\n"
    "X_DOT = 3.11548208 [km/s]\n"
    "Y_DOT = 0.47042605 [km/s]\n"
    "Z_DOT = -0.00101495 [km/s]\n"
    "SEMI_MAJOR_AXIS = 41399.5123 [km]\n"
    "ECCENTRICITY = 0.020842611\n"
    "INCLINATION = 0.117746 [deg]\n"
    "RA_OF_ASC_NODE = 17.604721 [deg]\n"
    "ARG_OF_PERICENTER = 218.242943 [deg]\n"
    "TRUE_ANOMALY = 41.922339 [deg]\n"
    "GM = 398600.4415 [km**3/s**2]\n"
    "MASS = 1913.000 [kg]\n"
    "SOLAR_RAD_AREA = 10.000 [m**2]\n"
    "SOLAR_RAD_COEFF = 1.300\n"
    "DRAG_AREA = 10.000 [m**2]\n"
    "DRAG_COEFF = 2.300\n"
    "MAN_EPOCH_IGNITION = 2000-06-03T09:00:34.1\n"
    "MAN_DURATION = 132.60 [s]\n"
    "MAN_DELTA_MASS = -18.418 [kg]\n"
    "MAN_REF_FRAME = EME2000\n"
    "MAN_DV_1 = -0.02325700 [km/s]\n"
    "MAN_DV_2 = 0.01683160 [km/s]\n"
    "MAN_DV_3 = -0.00893444 [km/s]\n"
    "MAN_EPOCH_IGNITION = 2000-06-05T18:59:21.0\n"
    "MAN_DURATION = 0.0 [s]\n"
    "MAN_DELTA_MASS = -1.469 [kg]\n"
    "MAN_REF_FRAME = RTN\n"
    "MAN_DV_1 = 0.00101500 [km/s]\n"
    "MAN_DV_2 = -0.00187300 [km/s]\n"
    "MAN_DV_3 = 0.00000000 [km/s]\n"
    "USER_DEFINED_EARTH_MODEL = WGS-84\n";

static void test_opm(void) {
  static myb_opm opm;
  static myb_opm back;
  myb_diag diag = {0};
  MT_CHECK(myb_opm_parse(opm_text, strlen(opm_text), &opm, &diag) == MYB_OK);
  MT_CHECK(strcmp(opm.object_name, "EUTELSAT W4") == 0 && strcmp(opm.ref_frame, "TOD") == 0 &&
           strcmp(opm.originator, "GSOC") == 0);
  MT_CHECK(opm.epoch.year == 2006 && opm.epoch.month == 6 && opm.epoch.day == 3);
  MT_CHECK_SAME(opm.creation_date.second, 0.123);
  MT_CHECK_SAME(opm.state[1], -40218.5751 * 1000.0);
  MT_CHECK_SAME(opm.state[3], 3.11548208 * 1000.0);
  MT_CHECK(opm.has_elements && !opm.anomaly_is_mean);
  MT_CHECK_SAME(opm.semi_major_axis, 41399.5123 * 1000.0);
  MT_CHECK_SAME(opm.inclination, 0.117746 * 0.017453292519943295);
  MT_CHECK_SAME(opm.gm, 398600.4415 * 1.0e9);
  MT_CHECK(opm.has_spacecraft);
  MT_CHECK_SAME(opm.mass, 1913.0);
  MT_CHECK_SAME(opm.drag_coeff, 2.3);
  MT_CHECK(opm.maneuver_count == 2);
  MT_CHECK(strcmp(opm.maneuvers[1].ref_frame, "RTN") == 0);
  MT_CHECK_SAME(opm.maneuvers[0].duration, 132.6);
  MT_CHECK_SAME(opm.maneuvers[0].dv[0], -0.023257 * 1000.0);
  MT_CHECK_SAME(opm.maneuvers[1].delta_mass, -1.469);
  MT_CHECK(opm.maneuvers[0].epoch_ignition.hour == 9 && opm.maneuvers[0].epoch_ignition.minute == 0);
  MT_CHECK_SAME(opm.maneuvers[0].epoch_ignition.second, 34.1);
  /* the state vector alone is a whole message */
  char text[4096];
  size_t length = 0;
  MT_CHECK(myb_opm_write(&opm, text, sizeof text, &length) == MYB_OK && strlen(text) == length);
  MT_CHECK(strstr(text, "X = 6655.9942\n") != NULL);
  MT_CHECK(strstr(text, "MAN_EPOCH_IGNITION = 2000-06-03T09:00:34.100\n") != NULL);
  MT_CHECK(myb_opm_parse(text, length, &back, &diag) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 6; ++i) {
    worst = fmax(worst, fabs(back.state[i] - opm.state[i]) /
                            (fabs(opm.state[i]) * 2.220446049250313e-16));
  }
  worst = fmax(worst, fabs(back.inclination - opm.inclination) /
                          (opm.inclination * 2.220446049250313e-16));
  worst = fmax(worst, fabs(back.gm - opm.gm) / (opm.gm * 2.220446049250313e-16));
  printf("OPM written and read back: values within %.2f ulp\n", worst);
  MT_CHECK(worst <= 1.5);
  MT_CHECK(back.maneuver_count == 2 && back.has_elements && back.has_spacecraft);
  MT_CHECK(memcmp(&back.epoch, &opm.epoch, sizeof opm.epoch) == 0);
  MT_CHECK(memcmp(&back.maneuvers[1].epoch_ignition, &opm.maneuvers[1].epoch_ignition,
                  sizeof(myb_calendar)) == 0);
  MT_CHECK_SAME(back.maneuvers[1].delta_mass, opm.maneuvers[1].delta_mass);
  /* errors */
  MT_CHECK(myb_opm_parse("CCSDS_OPM_VERS = 2.0\nX = 1\n", 27, &back, &diag) == MYB_E_FORMAT);
  MT_CHECK(strstr(diag.message, "CREATION_DATE") != NULL);
  MT_CHECK(myb_opm_parse("X = 1\n", 6, &back, &diag) == MYB_E_FORMAT && diag.line == 1);
  static const char bad_number[] = "CCSDS_OPM_VERS = 2.0\nCREATION_DATE = 2000-06-03T05:33:00\n"
                                   "X = seven [km]\n";
  MT_CHECK(myb_opm_parse(bad_number, strlen(bad_number), &back, &diag) == MYB_E_FORMAT &&
           diag.line == 3);
  static const char half[] =
      "CCSDS_OPM_VERS = 2.0\nCREATION_DATE = 2000-06-03T05:33:00\nORIGINATOR = A\n"
      "OBJECT_NAME = B\nOBJECT_ID = C\nCENTER_NAME = D\nREF_FRAME = E\nTIME_SYSTEM = UTC\n"
      "EPOCH = 2000-06-03T05:33:00\nX = 1\nY = 2\nZ = 3\nX_DOT = 4\nY_DOT = 5\nZ_DOT = 6\n";
  MT_CHECK(myb_opm_parse(half, strlen(half), &back, &diag) == MYB_OK);
  MT_CHECK(!back.has_elements && !back.has_spacecraft && back.maneuver_count == 0);
  char with_elements[1024];
  snprintf(with_elements, sizeof with_elements, "%sSEMI_MAJOR_AXIS = 7000\n", half);
  MT_CHECK(myb_opm_parse(with_elements, strlen(with_elements), &back, &diag) == MYB_E_FORMAT);
  snprintf(with_elements, sizeof with_elements, "%sMAN_DURATION = 1\n", half);
  MT_CHECK(myb_opm_parse(with_elements, strlen(with_elements), &back, &diag) == MYB_E_FORMAT);
  opm.epoch.day = 31; /* June has 30 */
  MT_CHECK(myb_opm_write(&opm, text, sizeof text, &length) == MYB_E_DOMAIN);
}

static const char omm_text[] =
    "CCSDS_OMM_VERS = 2.0\n"
    "CREATION_DATE = 2020-06-01T10:00:00\n"
    "ORIGINATOR = 18 SPCS\n"
    "OBJECT_NAME = ISS (ZARYA)\n"
    "OBJECT_ID = 1998-067A\n"
    "CENTER_NAME = EARTH\n"
    "REF_FRAME = TEME\n"
    "TIME_SYSTEM = UTC\n"
    "MEAN_ELEMENT_THEORY = SGP4\n"
    "EPOCH = 2020-06-01T03:22:38.887392\n"
    "MEAN_MOTION = 15.49438371\n"
    "ECCENTRICITY = 0.0002335\n"
    "INCLINATION = 51.6443\n"
    "RA_OF_ASC_NODE = 93.0000\n"
    "ARG_OF_PERICENTER = 45.4201\n"
    "MEAN_ANOMALY = 314.6988\n"
    "EPHEMERIS_TYPE = 0\n"
    "CLASSIFICATION_TYPE = U\n"
    "NORAD_CAT_ID = 25544\n"
    "ELEMENT_SET_NO = 999\n"
    "REV_AT_EPOCH = 22954\n"
    "BSTAR = 0.000024949\n"
    "MEAN_MOTION_DOT = 0.00000907\n"
    "MEAN_MOTION_DDOT = 0\n";

static void test_omm(void) {
  myb_omm omm;
  myb_omm back;
  myb_diag diag = {0};
  MT_CHECK(myb_omm_parse(omm_text, strlen(omm_text), &omm, &diag) == MYB_OK);
  MT_CHECK(strcmp(omm.object_name, "ISS (ZARYA)") == 0 &&
           strcmp(omm.mean_element_theory, "SGP4") == 0 && strcmp(omm.ref_frame, "TEME") == 0);
  MT_CHECK(omm.has_mean_motion && !omm.has_gm && omm.has_tle);
  MT_CHECK_SAME(omm.mean_motion, 15.49438371);
  MT_CHECK_SAME(omm.eccentricity, 0.0002335);
  MT_CHECK_SAME(omm.inclination, 51.6443);
  MT_CHECK_SAME(omm.epoch.second, 38.887392);
  MT_CHECK(omm.norad_cat_id == 25544 && omm.element_set_no == 999 && omm.rev_at_epoch == 22954 &&
           omm.classification_type == 'U' && omm.ephemeris_type == 0);
  MT_CHECK_SAME(omm.bstar, 0.000024949);
  MT_CHECK_SAME(omm.mean_motion_dot, 0.00000907);
  char text[2048];
  size_t length = 0;
  MT_CHECK(myb_omm_write(&omm, text, sizeof text, &length) == MYB_OK && strlen(text) == length);
  MT_CHECK(strstr(text, "EPOCH = 2020-06-01T03:22:38.887392\n") != NULL);
  MT_CHECK(strstr(text, "BSTAR = 2.4949e-05\n") != NULL);
  MT_CHECK(myb_omm_parse(text, length, &back, &diag) == MYB_OK);
  /* the numbers of the catalogue pass unchanged: every bit */
  MT_CHECK(memcmp(&back, &omm, sizeof omm) == 0);
  /* a semi-major axis in place of the mean motion; errors */
  static const char axis[] =
      "CCSDS_OMM_VERS = 2.0\nCREATION_DATE = 2020-06-01T10:00:00\nORIGINATOR = A\n"
      "OBJECT_NAME = B\nOBJECT_ID = C\nCENTER_NAME = EARTH\nREF_FRAME = TEME\nTIME_SYSTEM = UTC\n"
      "MEAN_ELEMENT_THEORY = DSST\nEPOCH = 2020-06-01T03:22:38\nSEMI_MAJOR_AXIS = 6800.5\n"
      "ECCENTRICITY = 0.001\nINCLINATION = 51\nRA_OF_ASC_NODE = 1\nARG_OF_PERICENTER = 2\n"
      "MEAN_ANOMALY = 3\nGM = 398600.4418\n";
  MT_CHECK(myb_omm_parse(axis, strlen(axis), &back, &diag) == MYB_OK);
  MT_CHECK(!back.has_mean_motion && back.has_gm && !back.has_tle);
  MT_CHECK_SAME(back.semi_major_axis, 6800.5);
  MT_CHECK(myb_omm_write(&back, text, sizeof text, &length) == MYB_OK);
  MT_CHECK(strstr(text, "SEMI_MAJOR_AXIS = 6800.5\n") != NULL && strstr(text, "BSTAR") == NULL);
  static const char no_motion[] =
      "CCSDS_OMM_VERS = 2.0\nCREATION_DATE = 2020-06-01T10:00:00\nORIGINATOR = A\n"
      "OBJECT_NAME = B\nOBJECT_ID = C\nCENTER_NAME = EARTH\nREF_FRAME = TEME\nTIME_SYSTEM = UTC\n"
      "MEAN_ELEMENT_THEORY = SGP4\nEPOCH = 2020-06-01T03:22:38\n"
      "ECCENTRICITY = 0.001\nINCLINATION = 51\nRA_OF_ASC_NODE = 1\nARG_OF_PERICENTER = 2\n"
      "MEAN_ANOMALY = 3\n";
  MT_CHECK(myb_omm_parse(no_motion, strlen(no_motion), &back, &diag) == MYB_E_FORMAT);
  MT_CHECK(myb_omm_parse("CCSDS_OPM_VERS = 2.0\n", 21, &back, &diag) == MYB_E_FORMAT);
}

static uint64_t next_random(uint64_t* state) {
  *state = *state * 6364136223846793005ull + 1442695040888963407ull;
  uint64_t z = *state;
  z ^= z >> 33;
  z *= 0xff51afd7ed558ccdull;
  z ^= z >> 33;
  return z;
}

/* States that are not decimals of a kilometre: what a propagation gives. Written and read
 * they come back exactly, but for the few that no number of kilometres times 1000 reaches. */
static void test_random(void) {
  enum { COUNT = 2000 };
  static myb_calendar epochs[COUNT];
  static double states[6 * COUNT];
  uint64_t state = 15;
  for (size_t k = 0; k < COUNT; ++k) {
    const double seconds = 0.001 * (double)(next_random(&state) % 86400000ull) +
                           (double)(next_random(&state) >> 11) * 0x1p-53 * 1.0e-3;
    epochs[k] = (myb_calendar){2026, 11, 10, (int)(seconds / 3600.0),
                               (int)(seconds / 60.0) % 60, 0.0};
    epochs[k].second = seconds - 3600.0 * epochs[k].hour - 60.0 * epochs[k].minute;
    for (int i = 0; i < 6; ++i) {
      const double u = (double)(next_random(&state) >> 11) * 0x1p-53;
      states[6 * k + (size_t)i] = (i < 3 ? 1.0e6 + 4.9e7 * u : 10.0 + 7990.0 * u) *
                                  (next_random(&state) & 1u ? 1.0 : -1.0);
    }
  }
  const myb_oem_segment_data segment = {.object_name = "TEST",
                                        .object_id = "2026-000A",
                                        .center_name = "EARTH",
                                        .ref_frame = "GCRF",
                                        .time_system = "TT",
                                        .count = COUNT,
                                        .epochs = epochs,
                                        .states = states};
  const myb_oem_data data = {.originator = "MYB",
                             .creation_date = {2026, 11, 10, 0, 0, 0.0},
                             .segment_count = 1,
                             .segments = &segment};
  size_t length = 0;
  MT_CHECK(myb_oem_write(&data, NULL, 0, &length) == MYB_E_OVERFLOW);
  char* text = malloc(length + 1);
  MT_CHECK(text != NULL);
  if (text == NULL) {
    return;
  }
  MT_CHECK(myb_oem_write(&data, text, length + 1, NULL) == MYB_OK);
  myb_oem* oem = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_oem_load(text, length, NULL, &oem, &diag) == MYB_OK);
  if (oem != NULL) {
    MT_CHECK(myb_oem_state_count(oem, 0) == COUNT);
    size_t inexact = 0;
    size_t epochs_changed = 0;
    double worst = 0.0;
    for (size_t k = 0; k < COUNT; ++k) {
      myb_calendar epoch;
      double y[6];
      (void)myb_oem_state(oem, 0, k, &epoch, y, NULL);
      epochs_changed += memcmp(&epoch, &epochs[k], sizeof epoch) != 0;
      for (int i = 0; i < 6; ++i) {
        const double v = states[6 * k + (size_t)i];
        if (y[i] != v) {
          ++inexact;
          const double ulp = nextafter(fabs(v), (double)INFINITY) - fabs(v);
          worst = fmax(worst, fabs(y[i] - v) / ulp);
        }
      }
    }
    printf("OEM of %d random states: %zu of %d numbers not exact (%.2f%%, at most %.0f ulp), "
           "%zu epochs changed\n",
           COUNT, inexact, 6 * COUNT, 100.0 * (double)inexact / (double)(6 * COUNT), worst,
           epochs_changed);
    MT_CHECK(epochs_changed == 0);
    MT_CHECK(worst <= 1.0 && inexact <= 6 * COUNT * 3 / 100);
    myb_oem_destroy(oem);
  } else {
    printf("line %ld: %s\n", diag.line, diag.message);
  }
  free(text);
}

/* ---- Covariances, the XML notation, the comprehensive message --------------------------------- */

/* The message of CCSDS 505.0-B in the style of its annex: units as attributes, comments,
 * an entity, a state with accelerations, a covariance. */
static const char oem_xml[] =
    "\xEF\xBB\xBF<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<!-- an ephemeris in XML -->\n"
    "<oem xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
    "     xsi:noNamespaceSchemaLocation=\"http://sanaregistry.org/r/ndmxml/ndmxml-1.0-master.xsd\"\n"
    "     id=\"CCSDS_OEM_VERS\" version=\"2.0\">\n"
    "  <header>\n"
    "    <COMMENT>a > b</COMMENT>\n"
    "    <CREATION_DATE>1996-11-04T17:22:31</CREATION_DATE>\n"
    "    <ORIGINATOR>NASA &amp; JPL</ORIGINATOR>\n"
    "  </header>\n"
    "  <body>\n"
    "    <segment>\n"
    "      <metadata>\n"
    "        <OBJECT_NAME>MARS GLOBAL SURVEYOR</OBJECT_NAME>\n"
    "        <OBJECT_ID>1996-062A</OBJECT_ID>\n"
    "        <CENTER_NAME>MARS BARYCENTER</CENTER_NAME>\n"
    "        <REF_FRAME>EME2000</REF_FRAME>\n"
    "        <TIME_SYSTEM>UTC</TIME_SYSTEM>\n"
    "        <START_TIME>1996-12-18T12:00:00.331</START_TIME>\n"
    "        <STOP_TIME>1996-12-18T12:01:00.331</STOP_TIME>\n"
    "      </metadata>\n"
    "      <data>\n"
    "        <stateVector>\n"
    "          <EPOCH>1996-12-18T12:00:00.331</EPOCH>\n"
    "          <X units=\"km\">2789.619</X><Y units=\"km\">-280.045</Y><Z>-1746.755</Z>\n"
    "          <X_DOT units=\"km/s\">4.73372</X_DOT><Y_DOT>-2.49586</Y_DOT><Z_DOT>-1.04195</Z_DOT>\n"
    "        </stateVector>\n"
    "        <stateVector>\n"
    "          <EPOCH>1996-12-18T12:01:00.331</EPOCH>\n"
    "          <X>2783.419</X><Y>-308.143</Y><Z>-1877.071</Z>\n"
    "          <X_DOT>5.18604</X_DOT><Y_DOT>-2.42124</Y_DOT><Z_DOT>-1.99608</Z_DOT>\n"
    "          <X_DDOT>0.001</X_DDOT><Y_DDOT>-0.002</Y_DDOT><Z_DDOT>3e-4</Z_DDOT>\n"
    "        </stateVector>\n"
    "        <covarianceMatrix>\n"
    "          <EPOCH>1996-12-18T12:00:30</EPOCH>\n"
    "          <COV_REF_FRAME>RTN</COV_REF_FRAME>\n"
    "          <CX_X>1</CX_X><CY_X>2</CY_X><CY_Y>3</CY_Y>\n"
    "          <CZ_X>4</CZ_X><CZ_Y>5</CZ_Y><CZ_Z>6</CZ_Z>\n"
    "          <CX_DOT_X>7</CX_DOT_X><CX_DOT_Y>8</CX_DOT_Y><CX_DOT_Z>9</CX_DOT_Z>\n"
    "          <CX_DOT_X_DOT>10</CX_DOT_X_DOT>\n"
    "          <CY_DOT_X>11</CY_DOT_X><CY_DOT_Y>12</CY_DOT_Y><CY_DOT_Z>13</CY_DOT_Z>\n"
    "          <CY_DOT_X_DOT>14</CY_DOT_X_DOT><CY_DOT_Y_DOT>15</CY_DOT_Y_DOT>\n"
    "          <CZ_DOT_X>16</CZ_DOT_X><CZ_DOT_Y>17</CZ_DOT_Y><CZ_DOT_Z>18</CZ_DOT_Z>\n"
    "          <CZ_DOT_X_DOT>19</CZ_DOT_X_DOT><CZ_DOT_Y_DOT>20</CZ_DOT_Y_DOT>\n"
    "          <CZ_DOT_Z_DOT>21</CZ_DOT_Z_DOT>\n"
    "        </covarianceMatrix>\n"
    "      </data>\n"
    "    </segment>\n"
    "  </body>\n"
    "</oem>\n";

/* The same message as a document with a type declaration (an internal subset whose comment,
 * entity values and processing instruction hold every character that could end it early) and
 * values given as CDATA sections. */
#define OEM_XML_HEAD                                                                             \
  "<?xml version=\"1.0\"?>\n"                                                                    \
  "<!DOCTYPE oem SYSTEM \"oem.dtd\" [\n"                                                         \
  "  <!-- it's a comment with ] and > and \" -->\n"                                              \
  "  <!ENTITY agency \"NASA > JPL ]\">\n"                                                        \
  "  <!ENTITY trap '<ORIGINATOR>WRONG</ORIGINATOR>'>\n"                                          \
  "  <?note ]> is not the end ?>\n"                                                             \
  "]>\n"                                                                                        \
  "<oem id=\"CCSDS_OEM_VERS\" version=\"2.0\">\n"                                               \
  "  <header>\n"                                                                                \
  "    <CREATION_DATE><![CDATA[1996-11-04T17:22:31]]></CREATION_DATE>\n"
#define OEM_XML_TAIL                                                                             \
  "  </header>\n"                                                                               \
  "  <body><segment>\n"                                                                         \
  "    <metadata>\n"                                                                            \
  "      <OBJECT_NAME>\n        <![CDATA[ MGS <1> &amp; ]]>\n      </OBJECT_NAME>\n"             \
  "      <OBJECT_ID>1996-062A</OBJECT_ID>\n"                                                    \
  "      <CENTER_NAME>MARS &agency;</CENTER_NAME>\n"                                            \
  "      <REF_FRAME>EME2000</REF_FRAME>\n"                                                      \
  "      <TIME_SYSTEM><![CDATA[UTC]]></TIME_SYSTEM>\n"                                          \
  "      <START_TIME>1996-12-18T12:00:00.331</START_TIME>\n"                                    \
  "      <STOP_TIME>1996-12-18T12:01:00.331</STOP_TIME>\n"                                      \
  "    </metadata>\n"                                                                           \
  "    <data>\n"                                                                                \
  "      <![CDATA[ between the elements: <stateVector> nothing reads it ]]>\n"                  \
  "      <stateVector>\n"                                                                       \
  "        <EPOCH><![CDATA[1996-12-18T12:00:00.331]]></EPOCH>\n"                                \
  "        <X><![CDATA[2789.619]]></X><Y> <![CDATA[ -280.045 ]]> </Y><Z>-1746.755</Z>\n"        \
  "        <X_DOT>4.73372</X_DOT><Y_DOT>-2.49586</Y_DOT><Z_DOT><![CDATA[-1.04195]]></Z_DOT>\n"  \
  "      </stateVector>\n"                                                                      \
  "    </data>\n"                                                                               \
  "  </segment></body>\n"                                                                       \
  "</oem>\n"

static void test_xml_sections(void) {
  myb_diag diag = {0};
  myb_oem* oem = NULL;
  static const char good[] =
      OEM_XML_HEAD "    <ORIGINATOR><![CDATA[A <b> &amp; C]]></ORIGINATOR>\n" OEM_XML_TAIL;
  MT_CHECK(myb_oem_load(good, sizeof good - 1, NULL, &oem, &diag) == MYB_OK);
  if (oem != NULL) {
    /* a section is taken as written: no entity in it, and its blanks at the ends dropped */
    MT_CHECK(strcmp(myb_oem_header(oem, "ORIGINATOR"), "A <b> &amp; C") == 0);
    MT_CHECK(strcmp(myb_oem_header(oem, "CREATION_DATE"), "1996-11-04T17:22:31") == 0);
    MT_CHECK(strcmp(myb_oem_meta(oem, 0, "OBJECT_NAME"), "MGS <1> &amp;") == 0);
    MT_CHECK(strcmp(myb_oem_meta(oem, 0, "TIME_SYSTEM"), "UTC") == 0);
    /* an entity of the internal subset is not known: as written */
    MT_CHECK(strcmp(myb_oem_meta(oem, 0, "CENTER_NAME"), "MARS &agency;") == 0);
    MT_CHECK(myb_oem_segment_count(oem) == 1 && myb_oem_state_count(oem, 0) == 1);
    myb_calendar epoch;
    double y[6];
    double unused[3];
    MT_CHECK(myb_oem_state(oem, 0, 0, &epoch, y, unused) == 1);
    MT_CHECK_SAME(epoch.second, 0.331);
    MT_CHECK_SAME(y[0], 2789.619 * 1000.0);
    MT_CHECK_SAME(y[1], -280.045 * 1000.0);
    MT_CHECK_SAME(y[5], -1.04195 * 1000.0);
    myb_oem_destroy(oem);
  } else {
    fprintf(stderr, "line %ld: %s\n", diag.line, diag.message);
  }
  /* a section among other characters, two sections, a section that does not end: no value */
  static const char* const bad[4] = {
      OEM_XML_HEAD "    <ORIGINATOR>A<![CDATA[B]]></ORIGINATOR>\n" OEM_XML_TAIL,
      OEM_XML_HEAD "    <ORIGINATOR><![CDATA[B]]>C</ORIGINATOR>\n" OEM_XML_TAIL,
      OEM_XML_HEAD "    <ORIGINATOR><![CDATA[A]]><![CDATA[B]]></ORIGINATOR>\n" OEM_XML_TAIL,
      OEM_XML_HEAD "    <ORIGINATOR><![CDATA[A</ORIGINATOR>\n" OEM_XML_TAIL};
  for (int k = 0; k < 4; ++k) {
    oem = NULL;
    MT_CHECK(myb_oem_load(bad[k], strlen(bad[k]), NULL, &oem, &diag) == MYB_E_FORMAT);
    MT_CHECK(oem == NULL);
    if (k < 3) {
      MT_CHECK(diag.line == 11); /* the line of the element */
    }
  }
  /* a line of data of the comprehensive message as one section */
  static const char ocm_text[] =
      "<?xml version=\"1.0\"?>\n"
      "<!DOCTYPE ocm [ <!ENTITY x \"]>\"> ]>\n"
      "<ocm id=\"CCSDS_OCM_VERS\" version=\"3.0\">\n"
      "<header><CREATION_DATE>2022-11-06T09:23:57</CREATION_DATE>"
      "<ORIGINATOR><![CDATA[A&B]]></ORIGINATOR></header>\n"
      "<body><segment><metadata><TIME_SYSTEM>UTC</TIME_SYSTEM>"
      "<EPOCH_TZERO>2022-12-18T14:28:15.1172</EPOCH_TZERO></metadata>\n"
      "<data><traj><CENTER_NAME>EARTH</CENTER_NAME><TRAJ_REF_FRAME>ICRF</TRAJ_REF_FRAME>"
      "<TRAJ_TYPE>CARTPV</TRAJ_TYPE>\n"
      "<trajLine> <![CDATA[ 0.0 2789.6 -280.0 -1746.8 4.73 -2.50 -1.04 ]]> </trajLine>\n"
      "<trajLine>10.0 2836.0 -305.0 -1757.2 4.72 -2.51 -1.05</trajLine>\n"
      "</traj></data></segment></body></ocm>\n";
  myb_ocm* ocm = NULL;
  MT_CHECK(myb_ocm_load(ocm_text, sizeof ocm_text - 1, NULL, &ocm, &diag) == MYB_OK);
  if (ocm != NULL) {
    MT_CHECK(strcmp(myb_ocm_header(ocm, "ORIGINATOR"), "A&B") == 0);
    MT_CHECK(myb_ocm_block_count(ocm, MYB_OCM_TRAJ) == 1 &&
             myb_ocm_line_count(ocm, MYB_OCM_TRAJ, 0) == 2);
    MT_CHECK(strcmp(myb_ocm_line_text(ocm, MYB_OCM_TRAJ, 0, 0),
                    "0.0 2789.6 -280.0 -1746.8 4.73 -2.50 -1.04") == 0);
    myb_ocm_destroy(ocm);
  } else {
    fprintf(stderr, "line %ld: %s\n", diag.line, diag.message);
  }
  printf("XML: values as CDATA sections, a document type with an internal subset\n");
}

static void test_covariance_and_xml(void) {
  myb_diag diag = {0};
  /* the covariances of the sample: the figure of the standard and one of whole numbers */
  myb_oem* oem = NULL;
  MT_CHECK(load(oem_text, &oem, &diag) == MYB_OK);
  if (oem == NULL) {
    return;
  }
  MT_CHECK(myb_oem_covariance_count(oem, 0) == 0 && myb_oem_covariance_count(oem, 1) == 2);
  myb_calendar epoch;
  const char* frame = "";
  double p[36];
  myb_oem_covariance(oem, 1, 0, &epoch, &frame, p);
  MT_CHECK(epoch.day == 28 && epoch.hour == 21 && epoch.minute == 29);
  MT_CHECK_SAME(epoch.second, 7.267);
  MT_CHECK(frame != NULL && strcmp(frame, "EME2000") == 0);
  MT_CHECK_SAME(p[0], 3.3313494e-04 * 1e6);
  MT_CHECK_SAME(p[6 * 1 + 0], 4.6189273e-04 * 1e6);
  MT_CHECK_SAME(p[6 * 0 + 1], 4.6189273e-04 * 1e6);
  MT_CHECK_SAME(p[6 * 3 + 2], 2.4849495e-07 * 1e6);
  MT_CHECK_SAME(p[6 * 2 + 3], 2.4849495e-07 * 1e6);
  MT_CHECK_SAME(p[6 * 5 + 5], 6.2244443e-10 * 1e6);
  MT_CHECK_SAME(p[6 * 5 + 0], -3.0413460e-07 * 1e6);
  myb_oem_covariance(oem, 1, 1, &epoch, &frame, p);
  MT_CHECK(frame == NULL && epoch.day == 29);
  MT_CHECK(p[0] == 1e6 && p[7] == 2e6 && p[14] == 3e6 && p[6 * 5] == 0.5e6 && p[5] == 0.5e6);
  MT_CHECK_SAME(p[35], 6e-6 * 1e6);

  /* written with the covariances in both notations and read back: every number */
  enum { STATES = 2, MATRICES = 2 };
  myb_calendar epochs[STATES];
  double states[6 * STATES];
  myb_calendar covariance_epochs[MATRICES];
  double covariances[36 * MATRICES];
  for (size_t k = 0; k < STATES; ++k) {
    (void)myb_oem_state(oem, 1, k, &epochs[k], states + 6 * k, NULL);
  }
  for (size_t k = 0; k < MATRICES; ++k) {
    myb_oem_covariance(oem, 1, k, &covariance_epochs[k], NULL, covariances + 36 * k);
  }
  const myb_oem_segment_data segment = {.object_name = "MGS <2>",
                                        .object_id = "1996-062A",
                                        .center_name = "MARS & CO",
                                        .ref_frame = "EME2000",
                                        .time_system = "UTC",
                                        .count = STATES,
                                        .epochs = epochs,
                                        .states = states,
                                        .covariance_count = MATRICES,
                                        .covariance_epochs = covariance_epochs,
                                        .covariance_frame = "RTN",
                                        .covariances = covariances};
  const myb_oem_data data = {.originator = "myb",
                             .creation_date = {2026, 11, 10, 9, 0, 0.0},
                             .segment_count = 1,
                             .segments = &segment};
  static char text[2][8192];
  for (int notation = 0; notation < 2; ++notation) {
    size_t length = 0;
    MT_CHECK(myb_oem_write_as(&data, (myb_odm_notation)notation, text[notation],
                              sizeof text[notation], &length) == MYB_OK);
    MT_CHECK((text[notation][0] == '<') == (notation == MYB_ODM_XML));
    myb_oem* back = NULL;
    MT_CHECK(myb_oem_load(text[notation], length, NULL, &back, &diag) == MYB_OK);
    if (back == NULL) {
      fprintf(stderr, "line %ld: %s\n", diag.line, diag.message);
      continue;
    }
    MT_CHECK(strcmp(myb_oem_header(back, "CCSDS_OEM_VERS"), "2.0") == 0);
    MT_CHECK(strcmp(myb_oem_meta(back, 0, "OBJECT_NAME"), "MGS <2>") == 0);
    MT_CHECK(strcmp(myb_oem_meta(back, 0, "CENTER_NAME"), "MARS & CO") == 0);
    MT_CHECK(myb_oem_state_count(back, 0) == STATES);
    MT_CHECK(myb_oem_covariance_count(back, 0) == MATRICES);
    for (size_t k = 0; k < STATES; ++k) {
      double y[6];
      (void)myb_oem_state(back, 0, k, &epoch, y, NULL);
      MT_CHECK(memcmp(y, states + 6 * k, sizeof y) == 0 && epoch.second == epochs[k].second);
    }
    for (size_t k = 0; k < MATRICES; ++k) {
      myb_oem_covariance(back, 0, k, &epoch, &frame, p);
      MT_CHECK(frame != NULL && strcmp(frame, "RTN") == 0);
      MT_CHECK(epoch.second == covariance_epochs[k].second);
      for (int i = 0; i < 36; ++i) {
        MT_CHECK(fabs(p[i] - covariances[36 * k + (size_t)i]) <= 2e-16 * fabs(p[i]));
      }
    }
    myb_oem_destroy(back);
  }
  MT_CHECK(strstr(text[1], "<OBJECT_NAME>MGS &lt;2&gt;</OBJECT_NAME>") != NULL);
  MT_CHECK(strstr(text[1], "<covarianceMatrix>") != NULL && strstr(text[1], "</oem>") != NULL);
  MT_CHECK(strstr(text[0], "\nCOVARIANCE_START\nEPOCH = 1996-12-28T21:29:07.267\n"
                           "COV_REF_FRAME = RTN\n0.00033313494\n0.00046189273 0.00067824216\n") !=
           NULL);
  myb_oem_destroy(oem);

  /* a message in XML as others write it */
  oem = NULL;
  MT_CHECK(load(oem_xml, &oem, &diag) == MYB_OK);
  if (oem != NULL) {
    MT_CHECK(strcmp(myb_oem_header(oem, "CCSDS_OEM_VERS"), "2.0") == 0);
    MT_CHECK(strcmp(myb_oem_header(oem, "ORIGINATOR"), "NASA & JPL") == 0);
    MT_CHECK(myb_oem_header(oem, "COMMENT") == NULL);
    MT_CHECK(myb_oem_segment_count(oem) == 1 && myb_oem_state_count(oem, 0) == 2);
    MT_CHECK(strcmp(myb_oem_meta(oem, 0, "CENTER_NAME"), "MARS BARYCENTER") == 0);
    double y[6];
    double a[3];
    MT_CHECK(myb_oem_state(oem, 0, 0, &epoch, y, a) == 1);
    MT_CHECK_SAME(y[0], 2789.619 * 1000.0);
    MT_CHECK_SAME(y[5], -1.04195 * 1000.0);
    MT_CHECK_SAME(epoch.second, 0.331);
    MT_CHECK(myb_oem_state(oem, 0, 1, &epoch, y, a) == 3);
    MT_CHECK_SAME(y[1], -308.143 * 1000.0);
    MT_CHECK_SAME(a[2], 3e-4 * 1000.0);
    MT_CHECK(myb_oem_covariance_count(oem, 0) == 1);
    myb_oem_covariance(oem, 0, 0, &epoch, &frame, p);
    MT_CHECK(strcmp(frame, "RTN") == 0 && epoch.second == 30.0);
    int k = 0;
    for (int row = 0; row < 6; ++row) {
      for (int column = 0; column <= row; ++column) {
        ++k;
        MT_CHECK(p[6 * row + column] == 1e6 * k && p[6 * column + row] == 1e6 * k);
      }
    }
    myb_oem_destroy(oem);
  } else {
    fprintf(stderr, "line %ld: %s\n", diag.line, diag.message);
  }

  /* what is refused */
  static const char* const bad[] = {
      /* a matrix of 20 numbers */
      "CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = A\nOBJECT_ID = B\nCENTER_NAME = C\n"
      "REF_FRAME = D\nTIME_SYSTEM = UTC\nSTART_TIME = 2020-01-01T00:00:00\n"
      "STOP_TIME = 2020-01-01T00:00:00\nMETA_STOP\n2020-01-01T00:00:00 1 2 3 4 5 6\n"
      "COVARIANCE_START\nEPOCH = 2020-01-01T00:00:00\n1\n2 3\n4 5 6\n7 8 9 10\n11 12 13 14 15\n"
      "16 17 18 19 20\nCOVARIANCE_STOP\n",
      /* numbers before the epoch */
      "CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = A\nOBJECT_ID = B\nCENTER_NAME = C\n"
      "REF_FRAME = D\nTIME_SYSTEM = UTC\nSTART_TIME = 2020-01-01T00:00:00\n"
      "STOP_TIME = 2020-01-01T00:00:00\nMETA_STOP\n2020-01-01T00:00:00 1 2 3 4 5 6\n"
      "COVARIANCE_START\n1 2\nCOVARIANCE_STOP\n",
      /* a root without the version, a tag that does not end */
      "<oem><header><CREATION_DATE>2020-01-01T00:00:00</CREATION_DATE></header></oem>",
      "<oem id=\"CCSDS_OEM_VERS\" version=\"2.0\"><header><ORIGINATOR",
      /* an OPM in the place of an OEM */
      "<opm id=\"CCSDS_OPM_VERS\" version=\"2.0\"></opm>"};
  for (size_t i = 0; i < sizeof bad / sizeof bad[0]; ++i) {
    myb_oem* none = NULL;
    MT_CHECK(load(bad[i], &none, &diag) == MYB_E_FORMAT && none == NULL);
  }
  myb_oem_segment_data wrong = segment;
  wrong.covariance_epochs = NULL;
  myb_oem_data wrong_data = data;
  wrong_data.segments = &wrong;
  size_t length = 0;
  MT_CHECK(myb_oem_write_as(&wrong_data, MYB_ODM_KVN, text[0], sizeof text[0], &length) ==
           MYB_E_DOMAIN);
  MT_CHECK(myb_oem_write_as(&data, (myb_odm_notation)7, text[0], sizeof text[0], &length) ==
           MYB_E_DOMAIN);

  /* the parameter messages: the covariance, the two notations */
  myb_opm opm;
  MT_CHECK(myb_opm_parse(opm_text, strlen(opm_text), &opm, &diag) == MYB_OK);
  opm.has_covariance = 1;
  snprintf(opm.cov_ref_frame, sizeof opm.cov_ref_frame, "%s", "RTN");
  for (int row = 0; row < 6; ++row) {
    for (int column = 0; column < 6; ++column) {
      opm.covariance[6 * row + column] =
          (row == column ? 100.0 : 0.25) * (1.0 + row + column) * (row < 3 ? 1.0 : 1e-3) *
          (column < 3 ? 1.0 : 1e-3);
    }
  }
  for (int notation = 0; notation < 2; ++notation) {
    size_t n = 0;
    MT_CHECK(myb_opm_write_as(&opm, (myb_odm_notation)notation, text[0], sizeof text[0], &n) ==
             MYB_OK);
    myb_opm back;
    MT_CHECK(myb_opm_parse(text[0], n, &back, &diag) == MYB_OK);
    MT_CHECK(back.has_covariance && strcmp(back.cov_ref_frame, "RTN") == 0);
    MT_CHECK(strcmp(back.object_name, opm.object_name) == 0 &&
             back.maneuver_count == opm.maneuver_count && back.has_elements == opm.has_elements);
    double worst = 0.0;
    for (int i = 0; i < 36; ++i) {
      worst = fmax(worst, fabs(back.covariance[i] - opm.covariance[i]) /
                              fabs(opm.covariance[i]));
    }
    for (int i = 0; i < 6; ++i) {
      worst = fmax(worst, fabs(back.state[i] - opm.state[i]) / fabs(opm.state[i]));
    }
    MT_CHECK(worst < 3e-16);
    MT_CHECK(back.semi_major_axis == opm.semi_major_axis && back.mass == opm.mass);
    if (opm.maneuver_count > 0) {
      MT_CHECK(back.maneuvers[0].dv[1] == opm.maneuvers[0].dv[1] &&
               strcmp(back.maneuvers[0].ref_frame, opm.maneuvers[0].ref_frame) == 0);
    }
    if (notation == MYB_ODM_XML) {
      MT_CHECK(strstr(text[0], "<opm ") != NULL && strstr(text[0], "<stateVector>") != NULL &&
               strstr(text[0], "<covarianceMatrix>\n<COV_REF_FRAME>RTN</COV_REF_FRAME>") != NULL);
    } else {
      MT_CHECK(strstr(text[0], "\nCOV_REF_FRAME = RTN\nCX_X = ") != NULL);
    }
  }
  /* an incomplete covariance */
  static const char partial[] = "\nCX_X = 1\nCY_X = 2\n";
  size_t n = 0;
  opm.has_covariance = 0;
  MT_CHECK(myb_opm_write(&opm, text[0], sizeof text[0] - sizeof partial, &n) == MYB_OK);
  memcpy(text[0] + n, partial, sizeof partial);
  myb_opm incomplete;
  MT_CHECK(myb_opm_parse(text[0], strlen(text[0]), &incomplete, &diag) == MYB_E_FORMAT);

  myb_omm omm;
  MT_CHECK(myb_omm_parse(omm_text, strlen(omm_text), &omm, &diag) == MYB_OK);
  omm.has_covariance = 1;
  for (int i = 0; i < 36; ++i) {
    omm.covariance[i] = 1.0 + (i / 6) * (i % 6); /* symmetric */
  }
  for (int notation = 0; notation < 2; ++notation) {
    MT_CHECK(myb_omm_write_as(&omm, (myb_odm_notation)notation, text[0], sizeof text[0], &n) ==
             MYB_OK);
    myb_omm back;
    MT_CHECK(myb_omm_parse(text[0], n, &back, &diag) == MYB_OK);
    MT_CHECK(back.mean_motion == omm.mean_motion && back.eccentricity == omm.eccentricity &&
             back.bstar == omm.bstar && back.norad_cat_id == omm.norad_cat_id &&
             back.has_tle == omm.has_tle && back.epoch.second == omm.epoch.second &&
             strcmp(back.mean_element_theory, omm.mean_element_theory) == 0);
    MT_CHECK(back.has_covariance && memcmp(back.covariance, omm.covariance, sizeof omm.covariance) == 0);
    MT_CHECK((strstr(text[0], "<tleParameters>") != NULL) == (notation == MYB_ODM_XML && omm.has_tle));
  }
  printf("covariances of OEM, OPM, OMM and the XML notation: written and read back in both\n");
}

static const char ocm_text[] =
    "CCSDS_OCM_VERS = 3.0\n"
    "COMMENT in the style of CCSDS 502.0-B-3, annex G\n"
    "CREATION_DATE = 2022-11-06T09:23:57\n"
    "ORIGINATOR = JAXA\n"
    "META_START\n"
    "OBJECT_NAME = OSPREY 5\n"
    "INTERNATIONAL_DESIGNATOR = 1998-999A\n"
    "TIME_SYSTEM = UTC\n"
    "EPOCH_TZERO = 2022-12-18T00:00:00.0000\n"
    "META_STOP\n"
    "TRAJ_START\n"
    "COMMENT geocentric\n"
    "TRAJ_ID = PREDICT\n"
    "CENTER_NAME = EARTH\n"
    "TRAJ_REF_FRAME = ICRF3\n"
    "TRAJ_TYPE = CARTPV\n"
    "0.0 2789.6 -280.0 -1746.8 4.73 -2.50 -1.04\n"
    "10.0 2836.9 -305.0 -1757.2 4.72 -2.49 -1.05\n"
    "2022-12-18T00:00:20 2884.1 -329.9 -1767.7 4.71 -2.48 -1.06\n"
    "TRAJ_STOP\n"
    "TRAJ_START\n"
    "TRAJ_TYPE = KEPLERIAN\n"
    "TRAJ_UNITS = [km,n/a,deg,deg,deg,deg]\n"
    "0 6778 0.001 51.6 10 20 30\n"
    "TRAJ_STOP\n"
    "PHYS_START\n"
    "WET_MASS = 100.0 [kg]\n"
    "DRAG_CONST_AREA = 2.5\n"
    "PHYS_STOP\n"
    "COV_START\n"
    "COV_REF_FRAME = TNW_INERTIAL\n"
    "COV_TYPE = CARTP\n"
    "0.0 1 0.1 2 0.2 0.3 3\n"
    "COV_STOP\n"
    "MAN_START\n"
    "MAN_ID = E_W_20160305B\n"
    "MAN_DEVICE_ID = THR_01\n"
    "MAN_COMPOSITION = TIME_RELATIVE,MAN_DURA,THR_INTERP,DV_X,DV_Y,DV_Z\n"
    "500.0 25.0 ON 0.001 0.002 -0.003\n"
    "MAN_STOP\n"
    "PERT_START\n"
    "GRAVITY_MODEL = EGM-96: 36D 36O\n"
    "PERT_STOP\n"
    "OD_START\n"
    "OD_ID = OD_#12\n"
    "OD_METHOD = SEQUENTIAL\n"
    "OD_STOP\n"
    "USER_START\n"
    "USER_DEFINED_CONSOLE_POC = MAXWELL SMART\n"
    "USER_STOP\n";

static void test_ocm(void) {
  myb_diag diag = {0};
  myb_ocm* ocm = NULL;
  MT_CHECK(myb_ocm_load(ocm_text, strlen(ocm_text), NULL, &ocm, &diag) == MYB_OK);
  if (ocm == NULL) {
    fprintf(stderr, "line %ld: %s\n", diag.line, diag.message);
    return;
  }
  /* the message is rebuilt from what was read, written in both notations and read again:
   * the three give the same answers */
  myb_ocm* versions[3] = {ocm, NULL, NULL};
  static char text[2][8192];
  {
    enum { MOST = 16 };
    static myb_ocm_pair pairs[64];
    static const char* lines[32];
    static myb_ocm_block_data blocks[MOST];
    static const char* const meta_keys[] = {"OBJECT_NAME", "INTERNATIONAL_DESIGNATOR",
                                            "TIME_SYSTEM", "EPOCH_TZERO"};
    static const char* const block_keys[] = {
        "TRAJ_ID", "CENTER_NAME", "TRAJ_REF_FRAME", "TRAJ_TYPE", "TRAJ_UNITS", "WET_MASS",
        "DRAG_CONST_AREA", "COV_REF_FRAME", "COV_TYPE", "MAN_ID", "MAN_DEVICE_ID",
        "MAN_COMPOSITION", "GRAVITY_MODEL", "OD_ID", "OD_METHOD", "USER_DEFINED_CONSOLE_POC"};
    size_t used = 0;
    size_t used_lines = 0;
    size_t count = 0;
    myb_ocm_pair meta[4];
    for (int i = 0; i < 4; ++i) {
      meta[i] = (myb_ocm_pair){meta_keys[i], myb_ocm_meta(ocm, meta_keys[i])};
    }
    for (int kind = 0; kind < 7; ++kind) {
      for (size_t b = 0; b < myb_ocm_block_count(ocm, (myb_ocm_kind)kind); ++b) {
        myb_ocm_block_data* block = &blocks[count++];
        *block = (myb_ocm_block_data){.kind = (myb_ocm_kind)kind, .pairs = pairs + used,
                                      .lines = lines + used_lines};
        for (size_t i = 0; i < sizeof block_keys / sizeof block_keys[0]; ++i) {
          const char* value = myb_ocm_value(ocm, (myb_ocm_kind)kind, b, block_keys[i]);
          if (value != NULL) {
            pairs[used++] = (myb_ocm_pair){block_keys[i], value};
            block->pair_count++;
          }
        }
        for (size_t i = 0; i < myb_ocm_line_count(ocm, (myb_ocm_kind)kind, b); ++i) {
          lines[used_lines++] = myb_ocm_line_text(ocm, (myb_ocm_kind)kind, b, i);
          block->line_count++;
        }
      }
    }
    MT_CHECK(count == 8);
    const myb_ocm_data data = {.originator = myb_ocm_header(ocm, "ORIGINATOR"),
                               .creation_date = {2022, 11, 6, 9, 23, 57.0},
                               .meta_count = 4,
                               .meta = meta,
                               .block_count = count,
                               .blocks = blocks};
    for (int notation = 0; notation < 2; ++notation) {
      size_t length = 0;
      MT_CHECK(myb_ocm_write(&data, (myb_odm_notation)notation, text[notation],
                             sizeof text[notation], &length) == MYB_OK);
      MT_CHECK(myb_ocm_load(text[notation], length, NULL, &versions[1 + notation], &diag) ==
               MYB_OK);
      if (versions[1 + notation] == NULL) {
        fprintf(stderr, "line %ld: %s\n", diag.line, diag.message);
      }
    }
    MT_CHECK(strstr(text[1], "<USER_DEFINED parameter=\"CONSOLE_POC\">MAXWELL SMART") != NULL);
    MT_CHECK(strstr(text[1], "<trajLine>0.0 2789.6 -280.0") != NULL);
    MT_CHECK(strstr(text[0], "\nTRAJ_START\nTRAJ_ID = PREDICT\n") != NULL);
    /* what cannot be written */
    size_t length = 0;
    myb_ocm_data wrong = data;
    wrong.meta_count = 2; /* no TIME_SYSTEM */
    MT_CHECK(myb_ocm_write(&wrong, MYB_ODM_KVN, text[0], sizeof text[0], &length) == MYB_E_DOMAIN);
    const myb_ocm_pair lower = {"lower", "x"};
    const myb_ocm_block_data bad_key = {.kind = MYB_OCM_PHYS, .pair_count = 1, .pairs = &lower};
    wrong = data;
    wrong.block_count = 1;
    wrong.blocks = &bad_key;
    MT_CHECK(myb_ocm_write(&wrong, MYB_ODM_KVN, text[0], sizeof text[0], &length) == MYB_E_DOMAIN);
    static const char* const one_line[1] = {"0 1"};
    const myb_ocm_block_data bad_lines = {.kind = MYB_OCM_PHYS, .line_count = 1, .lines = one_line};
    wrong.blocks = &bad_lines;
    MT_CHECK(myb_ocm_write(&wrong, MYB_ODM_XML, text[0], sizeof text[0], &length) == MYB_E_DOMAIN);
  }
  for (int v = 0; v < 3; ++v) {
    const myb_ocm* o = versions[v];
    if (o == NULL) {
      continue;
    }
    MT_CHECK(strcmp(myb_ocm_header(o, "CCSDS_OCM_VERS"), "3.0") == 0);
    MT_CHECK(strcmp(myb_ocm_header(o, "ORIGINATOR"), "JAXA") == 0);
    MT_CHECK(strcmp(myb_ocm_meta(o, "OBJECT_NAME"), "OSPREY 5") == 0);
    MT_CHECK(myb_ocm_meta(o, "TRAJ_ID") == NULL);
    static const size_t counts[7] = {2, 1, 1, 1, 1, 1, 1};
    for (int kind = 0; kind < 7; ++kind) {
      MT_CHECK(myb_ocm_block_count(o, (myb_ocm_kind)kind) == counts[kind]);
    }
    MT_CHECK(strcmp(myb_ocm_value(o, MYB_OCM_TRAJ, 0, "TRAJ_TYPE"), "CARTPV") == 0);
    MT_CHECK(strcmp(myb_ocm_value(o, MYB_OCM_TRAJ, 1, "TRAJ_TYPE"), "KEPLERIAN") == 0);
    MT_CHECK(myb_ocm_value(o, MYB_OCM_TRAJ, 1, "TRAJ_ID") == NULL);
    MT_CHECK(strcmp(myb_ocm_value(o, MYB_OCM_PHYS, 0, "WET_MASS"), "100.0 [kg]") == 0);
    MT_CHECK(strcmp(myb_ocm_value(o, MYB_OCM_PERT, 0, "GRAVITY_MODEL"), "EGM-96: 36D 36O") == 0);
    MT_CHECK(strcmp(myb_ocm_value(o, MYB_OCM_OD, 0, "OD_ID"), "OD_#12") == 0);
    MT_CHECK(strcmp(myb_ocm_value(o, MYB_OCM_USER, 0, "USER_DEFINED_CONSOLE_POC"),
                    "MAXWELL SMART") == 0);
    MT_CHECK(myb_ocm_line_count(o, MYB_OCM_TRAJ, 0) == 3 &&
             myb_ocm_line_count(o, MYB_OCM_TRAJ, 1) == 1 &&
             myb_ocm_line_count(o, MYB_OCM_PHYS, 0) == 0 &&
             myb_ocm_line_count(o, MYB_OCM_COV, 0) == 1 &&
             myb_ocm_line_count(o, MYB_OCM_MAN, 0) == 1);
    /* a trajectory in SI: relative and absolute times */
    int absolute = 1;
    myb_calendar epoch;
    double offset = -1.0;
    double state[9];
    int count = 0;
    MT_CHECK(myb_ocm_state(o, 0, 1, &absolute, &epoch, &offset, state, &count) == MYB_OK);
    MT_CHECK(!absolute && offset == 10.0 && count == 6 && state[6] == 0.0);
    MT_CHECK_SAME(state[0], 2836.9 * 1000.0);
    MT_CHECK_SAME(state[5], -1.05 * 1000.0);
    MT_CHECK(myb_ocm_state(o, 0, 2, &absolute, &epoch, &offset, state, &count) == MYB_OK);
    MT_CHECK(absolute && epoch.year == 2022 && epoch.second == 20.0);
    MT_CHECK(myb_ocm_state(o, 1, 0, &absolute, &epoch, &offset, state, &count) == MYB_E_DOMAIN);
    /* the lines as they are */
    double values[8];
    size_t found = 0;
    MT_CHECK(myb_ocm_line(o, MYB_OCM_TRAJ, 1, 0, &absolute, NULL, &offset, values, 8, &found) ==
             MYB_OK);
    MT_CHECK(!absolute && offset == 0.0 && found == 6 && values[0] == 6778.0 && values[5] == 30.0);
    MT_CHECK(myb_ocm_line(o, MYB_OCM_COV, 0, 0, NULL, NULL, NULL, values, 2, &found) == MYB_OK);
    MT_CHECK(found == 6 && values[0] == 1.0 && values[1] == 0.1);
    MT_CHECK(myb_ocm_line(o, MYB_OCM_MAN, 0, 0, NULL, NULL, NULL, values, 8, &found) ==
             MYB_E_FORMAT); /* a word among the numbers */
    MT_CHECK(strcmp(myb_ocm_line_text(o, MYB_OCM_MAN, 0, 0),
                    "500.0 25.0 ON 0.001 0.002 -0.003") == 0);
  }
  myb_ocm_destroy(versions[1]);
  myb_ocm_destroy(versions[2]);
  myb_ocm_destroy(ocm);

  /* a line from numbers */
  char line[256];
  size_t length = 0;
  const double numbers[3] = {7000.0, -0.1, 1e-9};
  const myb_calendar label = {2022, 12, 18, 0, 0, 20.5};
  MT_CHECK(myb_ocm_format_line(NULL, 12.5, numbers, 3, line, sizeof line, &length) == MYB_OK);
  printf("a line of an OCM from numbers: %s\n", line);
  MT_CHECK(strcmp(line, "12.5 7000.0 -0.1 1e-09") == 0 && length == strlen(line));
  MT_CHECK(myb_ocm_format_line(&label, 0.0, numbers, 1, line, sizeof line, &length) == MYB_OK);
  MT_CHECK(strcmp(line, "2022-12-18T00:00:20.500 7000.0") == 0);
  MT_CHECK(myb_ocm_format_line(NULL, 0.0, numbers, 3, line, 4, &length) == MYB_E_OVERFLOW);
  MT_CHECK(myb_ocm_format_line(NULL, (double)NAN, numbers, 3, line, sizeof line, &length) ==
           MYB_E_DOMAIN);

  /* what is refused */
  static const char* const bad[] = {
      "CCSDS_OEM_VERS = 2.0\nMETA_START\nTIME_SYSTEM = UTC\nEPOCH_TZERO = 2020-001T00:00:00\n"
      "META_STOP\n",
      "CCSDS_OCM_VERS = 3.0\nMETA_START\nTIME_SYSTEM = UTC\nMETA_STOP\n",
      "CCSDS_OCM_VERS = 3.0\nMETA_START\nTIME_SYSTEM = UTC\nEPOCH_TZERO = 2020-001T00:00:00\n",
      "CCSDS_OCM_VERS = 3.0\nMETA_START\nTIME_SYSTEM = UTC\nEPOCH_TZERO = 2020-001T00:00:00\n"
      "META_STOP\nTRAJ_START\nTRAJ_TYPE = CARTP\n0 1 2 3\n",
      "CCSDS_OCM_VERS = 3.0\nMETA_START\nTIME_SYSTEM = UTC\nEPOCH_TZERO = 2020-001T00:00:00\n"
      "META_STOP\nPHYS_START\n0 1 2 3\nPHYS_STOP\n",
      "CCSDS_OCM_VERS = 3.0\nMETA_START\nTIME_SYSTEM = UTC\nEPOCH_TZERO = 2020-001T00:00:00\n"
      "META_STOP\nTRAJ_START\nTRAJ_STOP\nWET_MASS = 1\n",
      "CCSDS_OCM_VERS = 3.0\nMETA_START\nTIME_SYSTEM = UTC\nEPOCH_TZERO = 2020-001T00:00:00\n"
      "META_STOP\nTRAJ_START\nCOV_STOP\n"};
  for (size_t i = 0; i < sizeof bad / sizeof bad[0]; ++i) {
    myb_ocm* none = NULL;
    MT_CHECK(myb_ocm_load(bad[i], strlen(bad[i]), NULL, &none, &diag) == MYB_E_FORMAT &&
             none == NULL);
  }
  printf("OCM: 8 blocks of the 7 kinds read, written in both notations and read back\n");
}

int main(void) {
  test_random();
  test_oem();
  test_opm();
  test_omm();
  test_covariance_and_xml();
  test_xml_sections();
  test_ocm();
  return mt_finish("test_ccsds");
}
