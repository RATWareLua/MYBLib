/* SPDX-License-Identifier: BSD-3-Clause */
/* CCSDS ephemeris messages, KVN and XML (myb_oem_load): the keywords of the header and of
 * each segment, every state with its epoch, every covariance; then the message rebuilt from
 * what was read, written in both notations (myb_oem_write_as) and read again: a message the
 * writer accepts must read back with the same texts, dates, states and covariances (these to
 * the unit of the last place the conversion between km and m may cost). */
#include <math.h>
#include <stdlib.h>

#include "fuzz.h"
#include "myb/ccsds.h"

MYB_FUZZ_NO_INITIALIZE

enum { MAX_SEGMENTS = 16, MAX_STATES = 512, MAX_COVARIANCES = 64 };

static const char* const header_keys[] = {"CCSDS_OEM_VERS", "CREATION_DATE", "ORIGINATOR",
                                          "MESSAGE_ID", "COMMENT"};
static const char* const meta_keys[] = {"OBJECT_NAME", "OBJECT_ID", "CENTER_NAME", "REF_FRAME",
                                        "REF_FRAME_EPOCH", "TIME_SYSTEM", "START_TIME",
                                        "USEABLE_START_TIME", "USEABLE_STOP_TIME", "STOP_TIME",
                                        "INTERPOLATION", "INTERPOLATION_DEGREE"};

static const char* or_empty(const char* text) {
  return text != NULL ? text : "";
}

/* Within a few units of the last place (the file keeps km: value / 1000, then 1000 times what
 * is read), or a part in 1e-300 of nothing where the quotient was subnormal. */
static int close_to(double a, double b) {
  return a == b || fabs(a - b) <= 4.0 * 2.220446049250313e-16 * fabs(a) + 1e-290;
}

static int same_date(const myb_calendar* a, const myb_calendar* b) {
  return a->year == b->year && a->month == b->month && a->day == b->day && a->hour == b->hour &&
         a->minute == b->minute && a->second == b->second;
}

/* The covariances of a segment, when they share one frame (the writer writes one for all). */
static size_t covariances_of(const myb_oem* oem, size_t k, myb_calendar* epochs, double* matrices,
                             const char** frame) {
  size_t count = myb_oem_covariance_count(oem, k);
  count = count < MAX_COVARIANCES ? count : MAX_COVARIANCES;
  *frame = NULL;
  for (size_t i = 0; i < count; ++i) {
    const char* this_frame = NULL;
    myb_oem_covariance(oem, k, i, &epochs[i], &this_frame, matrices + 36 * i);
    if (i == 0) {
      *frame = this_frame;
    } else if ((this_frame == NULL) != (*frame == NULL) ||
               (this_frame != NULL && strcmp(this_frame, *frame) != 0)) {
      return 0;
    }
  }
  return count;
}

/* Writes what was read back as a message and reads it again: the same texts, dates, states
 * and covariances. */
static void round_trip(const myb_oem* oem) {
  const size_t segments = myb_oem_segment_count(oem);
  if (segments == 0 || segments > MAX_SEGMENTS) {
    return;
  }
  myb_oem_segment_data parts[MAX_SEGMENTS];
  myb_calendar* epochs = calloc(segments * MAX_STATES, sizeof(myb_calendar));
  double* states = calloc(segments * MAX_STATES * 6, sizeof(double));
  myb_calendar* covariance_epochs = calloc(segments * MAX_COVARIANCES, sizeof(myb_calendar));
  double* covariances = calloc(segments * MAX_COVARIANCES * 36, sizeof(double));
  MYB_CHECK(epochs != NULL && states != NULL && covariance_epochs != NULL &&
            covariances != NULL);
  for (size_t k = 0; k < segments; ++k) {
    size_t count = myb_oem_state_count(oem, k);
    count = count < MAX_STATES ? count : MAX_STATES;
    for (size_t i = 0; i < count; ++i) {
      (void)myb_oem_state(oem, k, i, &epochs[k * MAX_STATES + i],
                          &states[(k * MAX_STATES + i) * 6], NULL);
    }
    const char* frame = NULL;
    const size_t matrices = covariances_of(oem, k, covariance_epochs + k * MAX_COVARIANCES,
                                           covariances + k * MAX_COVARIANCES * 36, &frame);
    parts[k] = (myb_oem_segment_data){
        .object_name = or_empty(myb_oem_meta(oem, k, "OBJECT_NAME")),
        .object_id = or_empty(myb_oem_meta(oem, k, "OBJECT_ID")),
        .center_name = or_empty(myb_oem_meta(oem, k, "CENTER_NAME")),
        .ref_frame = or_empty(myb_oem_meta(oem, k, "REF_FRAME")),
        .time_system = or_empty(myb_oem_meta(oem, k, "TIME_SYSTEM")),
        .interpolation = NULL,
        .count = count,
        .epochs = epochs + k * MAX_STATES,
        .states = states + k * MAX_STATES * 6,
        .covariance_count = matrices,
        .covariance_epochs = covariance_epochs + k * MAX_COVARIANCES,
        .covariances = covariances + k * MAX_COVARIANCES * 36,
        .covariance_frame = frame};
  }
  myb_calendar created = {2000, 1, 1, 0, 0, 0.0};
  const myb_oem_data message = {.originator = or_empty(myb_oem_header(oem, "ORIGINATOR")),
                                .creation_date = created,
                                .segment_count = segments,
                                .segments = parts};
  for (int notation = 0; notation < 2; ++notation) {
    size_t length = 0;
    const myb_status sized =
        myb_oem_write_as(&message, (myb_odm_notation)notation, NULL, 0, &length);
    if (sized != MYB_E_OVERFLOW) {
      MYB_CHECK(sized == MYB_E_DOMAIN);
      continue;
    }
    char* text = malloc(length + 1);
    MYB_CHECK(text != NULL);
    MYB_CHECK(myb_oem_write_as(&message, (myb_odm_notation)notation, text, length + 1, &length) ==
              MYB_OK);
    myb_oem* back = NULL;
    myb_diag diag = {0};
    const myb_status read = myb_oem_load(text, length, NULL, &back, &diag);
    MYB_CHECK(read == MYB_OK);
    MYB_CHECK(myb_oem_segment_count(back) == segments);
    MYB_CHECK(strcmp(myb_oem_header(back, "ORIGINATOR"), message.originator) == 0);
    for (size_t k = 0; k < segments; ++k) {
      const myb_oem_segment_data* part = &parts[k];
      MYB_CHECK(strcmp(myb_oem_meta(back, k, "OBJECT_NAME"), part->object_name) == 0 &&
                strcmp(myb_oem_meta(back, k, "OBJECT_ID"), part->object_id) == 0 &&
                strcmp(myb_oem_meta(back, k, "CENTER_NAME"), part->center_name) == 0 &&
                strcmp(myb_oem_meta(back, k, "REF_FRAME"), part->ref_frame) == 0 &&
                strcmp(myb_oem_meta(back, k, "TIME_SYSTEM"), part->time_system) == 0);
      MYB_CHECK(myb_oem_state_count(back, k) == part->count &&
                myb_oem_covariance_count(back, k) == part->covariance_count);
      for (size_t i = 0; i < part->count; ++i) {
        myb_calendar epoch;
        double state[6];
        (void)myb_oem_state(back, k, i, &epoch, state, NULL);
        MYB_CHECK(same_date(&epoch, &part->epochs[i]));
        for (int j = 0; j < 6; ++j) {
          MYB_CHECK(close_to(state[j], part->states[6 * i + (size_t)j]));
        }
      }
      for (size_t i = 0; i < part->covariance_count; ++i) {
        myb_calendar epoch;
        const char* frame = NULL;
        double matrix[36];
        myb_oem_covariance(back, k, i, &epoch, &frame, matrix);
        MYB_CHECK(same_date(&epoch, &part->covariance_epochs[i]));
        MYB_CHECK((frame == NULL) == (part->covariance_frame == NULL) &&
                  (frame == NULL || strcmp(frame, part->covariance_frame) == 0));
        for (int j = 0; j < 36; ++j) {
          MYB_CHECK(close_to(matrix[j], part->covariances[36 * i + (size_t)j]));
        }
      }
    }
    myb_oem_destroy(back);
    free(text);
  }
  free(covariances);
  free(covariance_epochs);
  free(states);
  free(epochs);
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_oem* oem = NULL;
  myb_diag diag = {0};
  if (myb_oem_load((const char*)data, size, NULL, &oem, &diag) != MYB_OK) {
    MYB_CHECK(oem == NULL);
    return 0;
  }
  for (size_t i = 0; i < sizeof header_keys / sizeof header_keys[0]; ++i) {
    (void)myb_oem_header(oem, header_keys[i]);
  }
  const size_t segments = myb_oem_segment_count(oem);
  for (size_t k = 0; k < segments && k < MAX_SEGMENTS; ++k) {
    for (size_t i = 0; i < sizeof meta_keys / sizeof meta_keys[0]; ++i) {
      (void)myb_oem_meta(oem, k, meta_keys[i]);
    }
    const size_t states = myb_oem_state_count(oem, k);
    for (size_t i = 0; i < states && i < MAX_STATES; ++i) {
      myb_calendar epoch;
      double state[6];
      double acceleration[3];
      const int given = myb_oem_state(oem, k, i, &epoch, state, acceleration);
      MYB_CHECK(given == 1 || given == 3);
      for (int j = 0; j < 6; ++j) {
        MYB_CHECK(isfinite(state[j]) && (j >= 3 || isfinite(acceleration[j])));
      }
    }
    const size_t covariances = myb_oem_covariance_count(oem, k);
    for (size_t i = 0; i < covariances && i < MAX_COVARIANCES; ++i) {
      myb_calendar epoch;
      const char* frame = NULL;
      double covariance[36];
      myb_oem_covariance(oem, k, i, &epoch, &frame, covariance);
    }
  }
  round_trip(oem);
  myb_oem_destroy(oem);
  return 0;
}
