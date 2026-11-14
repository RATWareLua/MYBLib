/* SPDX-License-Identifier: BSD-3-Clause */
/* CCSDS parameter messages, KVN and XML (myb_opm_parse): a message read is written in both
 * notations (myb_opm_write_as) and read again: a message the writer accepts must read back
 * with the same texts and dates, the numbers within the unit of the last place that the
 * conversion between the units of the file and SI may cost, the same manoeuvres and flags. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fuzz.h"
#include "myb/ccsds.h"

MYB_FUZZ_NO_INITIALIZE

/* Within a few units of the last place (the file keeps km, degrees: value / scale, then
 * scale times what is read), or a part in 1e-300 of nothing where the quotient was
 * subnormal. */
static int close_to(double a, double b) {
  return a == b || fabs(a - b) <= 4.0 * 2.220446049250313e-16 * fabs(a) + 1e-290;
}

static int same_date(const myb_calendar* a, const myb_calendar* b) {
  return a->year == b->year && a->month == b->month && a->day == b->day && a->hour == b->hour &&
         a->minute == b->minute && a->second == b->second;
}

static void round_trip(const myb_opm* opm) {
  for (int notation = 0; notation < 2; ++notation) {
    size_t length = 0;
    const myb_status sized = myb_opm_write_as(opm, (myb_odm_notation)notation, NULL, 0, &length);
    if (sized != MYB_E_OVERFLOW) {
      MYB_CHECK(sized == MYB_E_DOMAIN);
      continue;
    }
    char* text = malloc(length + 1);
    MYB_CHECK(text != NULL);
    MYB_CHECK(myb_opm_write_as(opm, (myb_odm_notation)notation, text, length + 1, &length) ==
              MYB_OK);
    myb_opm* back = malloc(sizeof *back);
    MYB_CHECK(back != NULL);
    myb_diag diag = {0};
    MYB_CHECK(myb_opm_parse(text, length, back, &diag) == MYB_OK);
    MYB_CHECK(back->maneuver_count == opm->maneuver_count);
    MYB_CHECK(back->has_elements == opm->has_elements &&
              back->has_spacecraft == opm->has_spacecraft &&
              back->has_covariance == opm->has_covariance);
    MYB_CHECK(strcmp(back->originator, opm->originator) == 0 &&
              strcmp(back->object_name, opm->object_name) == 0 &&
              strcmp(back->object_id, opm->object_id) == 0 &&
              strcmp(back->center_name, opm->center_name) == 0 &&
              strcmp(back->ref_frame, opm->ref_frame) == 0 &&
              strcmp(back->time_system, opm->time_system) == 0);
    MYB_CHECK(same_date(&back->creation_date, &opm->creation_date) &&
              same_date(&back->epoch, &opm->epoch));
    for (int i = 0; i < 6; ++i) {
      MYB_CHECK(close_to(back->state[i], opm->state[i]));
    }
    if (opm->has_elements) {
      MYB_CHECK(close_to(back->semi_major_axis, opm->semi_major_axis) &&
                close_to(back->eccentricity, opm->eccentricity) &&
                close_to(back->inclination, opm->inclination) &&
                close_to(back->ra_of_asc_node, opm->ra_of_asc_node) &&
                close_to(back->arg_of_pericenter, opm->arg_of_pericenter) &&
                close_to(back->anomaly, opm->anomaly) && close_to(back->gm, opm->gm) &&
                back->anomaly_is_mean == opm->anomaly_is_mean);
    }
    if (opm->has_spacecraft) {
      MYB_CHECK(back->mass == opm->mass && back->solar_rad_area == opm->solar_rad_area &&
                back->solar_rad_coeff == opm->solar_rad_coeff &&
                back->drag_area == opm->drag_area && back->drag_coeff == opm->drag_coeff);
    }
    for (size_t m = 0; m < opm->maneuver_count; ++m) {
      const myb_opm_maneuver* a = &opm->maneuvers[m];
      const myb_opm_maneuver* b = &back->maneuvers[m];
      MYB_CHECK(same_date(&a->epoch_ignition, &b->epoch_ignition) &&
                strcmp(a->ref_frame, b->ref_frame) == 0 && b->delta_mass == a->delta_mass &&
                close_to(b->duration, a->duration) &&
                close_to(b->dv[0], a->dv[0]) && close_to(b->dv[1], a->dv[1]) &&
                close_to(b->dv[2], a->dv[2]));
    }
    if (opm->has_covariance) {
      MYB_CHECK(strcmp(back->cov_ref_frame, opm->cov_ref_frame) == 0);
      for (int k = 0; k < 36; ++k) {
        MYB_CHECK(close_to(back->covariance[k], opm->covariance[k]));
      }
    }
    free(back);
    free(text);
  }
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_opm* opm = malloc(sizeof *opm);
  MYB_CHECK(opm != NULL);
  myb_diag diag = {0};
  if (myb_opm_parse((const char*)data, size, opm, &diag) == MYB_OK) {
    MYB_CHECK(opm->maneuver_count <= MYB_OPM_MAX_MANEUVERS);
    for (int i = 0; i < 6; ++i) {
      MYB_CHECK(isfinite(opm->state[i]));
    }
    round_trip(opm);
  }
  free(opm);
  return 0;
}
