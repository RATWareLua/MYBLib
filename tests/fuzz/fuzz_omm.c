/* SPDX-License-Identifier: BSD-3-Clause */
/* CCSDS mean-elements messages, KVN and XML (myb_omm_parse): a message read is written in
 * both notations (myb_omm_write_as) and read again — the numbers of a catalogue pass
 * unchanged, so the message read back must equal the one written, bit for bit — and turned
 * into a two-line element set (ext/sgp4: myb_tle_from_omm) that is formatted, read back and
 * propagated. */
#include <math.h>
#include <stdlib.h>

#include "fuzz.h"
#include "myb/ccsds.h"
#include "myb_sgp4.h"

MYB_FUZZ_NO_INITIALIZE

static int same_date(const myb_calendar* a, const myb_calendar* b) {
  return a->year == b->year && a->month == b->month && a->day == b->day && a->hour == b->hour &&
         a->minute == b->minute && a->second == b->second;
}

/* What the writer wrote, as it was: the texts, the dates, every number (none is scaled). */
static void same_message(const myb_omm* a, const myb_omm* b) {
  MYB_CHECK(strcmp(a->originator, b->originator) == 0 &&
            strcmp(a->object_name, b->object_name) == 0 &&
            strcmp(a->object_id, b->object_id) == 0 &&
            strcmp(a->center_name, b->center_name) == 0 &&
            strcmp(a->ref_frame, b->ref_frame) == 0 &&
            strcmp(a->time_system, b->time_system) == 0 &&
            strcmp(a->mean_element_theory, b->mean_element_theory) == 0);
  MYB_CHECK(same_date(&a->creation_date, &b->creation_date) && same_date(&a->epoch, &b->epoch));
  MYB_CHECK(a->has_mean_motion ? a->mean_motion == b->mean_motion
                               : a->semi_major_axis == b->semi_major_axis);
  MYB_CHECK(a->eccentricity == b->eccentricity && a->inclination == b->inclination &&
            a->ra_of_asc_node == b->ra_of_asc_node &&
            a->arg_of_pericenter == b->arg_of_pericenter && a->mean_anomaly == b->mean_anomaly);
  MYB_CHECK(!a->has_gm || a->gm == b->gm);
  if (a->has_tle) {
    MYB_CHECK(a->ephemeris_type == b->ephemeris_type &&
              a->classification_type == b->classification_type &&
              a->norad_cat_id == b->norad_cat_id && a->element_set_no == b->element_set_no &&
              a->rev_at_epoch == b->rev_at_epoch && a->bstar == b->bstar &&
              a->mean_motion_dot == b->mean_motion_dot &&
              a->mean_motion_ddot == b->mean_motion_ddot);
  }
  if (a->has_covariance) {
    MYB_CHECK(strcmp(a->cov_ref_frame, b->cov_ref_frame) == 0);
    for (int k = 0; k < 36; ++k) {
      MYB_CHECK(a->covariance[k] == b->covariance[k]);
    }
  }
}

static void round_trip(const myb_omm* omm) {
  for (int notation = 0; notation < 2; ++notation) {
    size_t length = 0;
    const myb_status sized = myb_omm_write_as(omm, (myb_odm_notation)notation, NULL, 0, &length);
    if (sized != MYB_E_OVERFLOW) {
      MYB_CHECK(sized == MYB_E_DOMAIN);
      continue;
    }
    char* text = malloc(length + 1);
    MYB_CHECK(text != NULL);
    MYB_CHECK(myb_omm_write_as(omm, (myb_odm_notation)notation, text, length + 1, &length) ==
              MYB_OK);
    myb_omm back;
    myb_diag diag = {0};
    MYB_CHECK(myb_omm_parse(text, length, &back, &diag) == MYB_OK);
    MYB_CHECK(back.has_mean_motion == omm->has_mean_motion && back.has_tle == omm->has_tle &&
              back.has_gm == omm->has_gm && back.has_covariance == omm->has_covariance);
    same_message(omm, &back);
    free(text);
  }
}

static void as_tle(const myb_omm* omm) {
  myb_tle tle;
  if (myb_tle_from_omm(omm, &tle) != MYB_OK) {
    return;
  }
  char line1[70];
  char line2[70];
  if (myb_tle_format(&tle, line1, line2) == MYB_OK) {
    myb_tle back;
    myb_diag diag = {0};
    MYB_CHECK(myb_tle_parse(line1, line2, 0, &back, &diag) == MYB_OK);
  }
  myb_sgp4 sgp4;
  if (myb_sgp4_init(&tle, MYB_SGP4_WGS72, MYB_SGP4_IMPROVED, &sgp4) == MYB_OK) {
    double r[3];
    double v[3];
    if (myb_sgp4_propagate(&sgp4, 720.0, r, v, NULL) == MYB_OK) {
      for (int i = 0; i < 3; ++i) {
        MYB_CHECK(isfinite(r[i]) && isfinite(v[i]));
      }
    }
  }
  myb_omm again;
  (void)myb_omm_from_tle(&tle, &again);
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_omm omm;
  myb_diag diag = {0};
  if (myb_omm_parse((const char*)data, size, &omm, &diag) == MYB_OK) {
    round_trip(&omm);
    as_tle(&omm);
  }
  return 0;
}
