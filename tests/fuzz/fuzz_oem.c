/* SPDX-License-Identifier: BSD-3-Clause */
/* CCSDS ephemeris messages, KVN and XML (myb_oem_load): the keywords of the header and of
 * each segment, every state with its epoch, every covariance; then the message rebuilt from
 * what was read, written in both notations (myb_oem_write_as) and read again: a message the
 * writer accepts must read back. */
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

/* Writes what was read back as a message and reads it again. */
static void round_trip(const myb_oem* oem) {
  const size_t segments = myb_oem_segment_count(oem);
  if (segments == 0 || segments > MAX_SEGMENTS) {
    return;
  }
  myb_oem_segment_data parts[MAX_SEGMENTS];
  myb_calendar* epochs = calloc(segments * MAX_STATES, sizeof(myb_calendar));
  double* states = calloc(segments * MAX_STATES * 6, sizeof(double));
  MYB_CHECK(epochs != NULL && states != NULL);
  for (size_t k = 0; k < segments; ++k) {
    size_t count = myb_oem_state_count(oem, k);
    count = count < MAX_STATES ? count : MAX_STATES;
    for (size_t i = 0; i < count; ++i) {
      (void)myb_oem_state(oem, k, i, &epochs[k * MAX_STATES + i],
                          &states[(k * MAX_STATES + i) * 6], NULL);
    }
    parts[k] = (myb_oem_segment_data){
        .object_name = or_empty(myb_oem_meta(oem, k, "OBJECT_NAME")),
        .object_id = or_empty(myb_oem_meta(oem, k, "OBJECT_ID")),
        .center_name = or_empty(myb_oem_meta(oem, k, "CENTER_NAME")),
        .ref_frame = or_empty(myb_oem_meta(oem, k, "REF_FRAME")),
        .time_system = or_empty(myb_oem_meta(oem, k, "TIME_SYSTEM")),
        .interpolation = NULL,
        .count = count,
        .epochs = epochs + k * MAX_STATES,
        .states = states + k * MAX_STATES * 6};
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
    myb_oem_destroy(back);
    free(text);
  }
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
