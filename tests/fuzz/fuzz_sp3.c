/* SPDX-License-Identifier: BSD-3-Clause */
/* SP3-c and SP3-d precise orbits (myb_sp3_load): the satellites, epochs, frame and time
 * system, and of every record the state, the clock, the standard deviations and the flags;
 * then the file rebuilt from what was read (myb_sp3_write) and read again: a description the
 * writer accepts must read back with the same satellites, labels to the 8 decimals of their
 * seconds, records present where they were and their numbers to the 6 decimals of the
 * fields, the same flags. */
#include <math.h>
#include <stdlib.h>

#include "fuzz.h"
#include "myb/kin.h"

MYB_FUZZ_NO_INITIALIZE

enum { MAX_RECORDS = 1 << 14 };

/* Within half a unit of the last decimal written (in the unit of the field: `unit` of the
 * library's), and the rounding of the conversion. */
static int within(double a, double b, double unit) {
  return fabs(a - b) <= 0.5 * unit * (1.0 + 1e-9) + 1e-15 * fabs(a);
}

static void round_trip(const myb_sp3* sp3) {
  const size_t satellites = myb_sp3_satellite_count(sp3);
  const size_t epochs = myb_sp3_epoch_count(sp3);
  const size_t cells = satellites * epochs;
  if (epochs == 0 || cells > MAX_RECORDS) {
    return;
  }
  const char** ids = malloc(satellites * sizeof *ids);
  myb_calendar* labels = malloc(epochs * sizeof *labels);
  double* positions = malloc(3 * cells * sizeof(double));
  double* velocities = malloc(3 * cells * sizeof(double));
  double* clocks = malloc(cells * sizeof(double));
  double* rates = malloc(cells * sizeof(double));
  unsigned char* flags = malloc(cells);
  MYB_CHECK(ids != NULL && labels != NULL && positions != NULL && velocities != NULL &&
            clocks != NULL && rates != NULL && flags != NULL);
  for (size_t s = 0; s < satellites; ++s) {
    ids[s] = myb_sp3_satellite(sp3, s);
  }
  for (size_t e = 0; e < epochs; ++e) {
    myb_sp3_epoch(sp3, e, &labels[e]);
    for (size_t s = 0; s < satellites; ++s) {
      const size_t cell = e * satellites + s;
      const int given = myb_sp3_state(sp3, s, e, positions + 3 * cell, velocities + 3 * cell);
      if (!(given & 1)) {
        positions[3 * cell] = (double)NAN;
      }
      if (!(given & 2)) {
        velocities[3 * cell] = (double)NAN;
      }
      const int timed = myb_sp3_clock(sp3, s, e, &clocks[cell], &rates[cell]);
      if (!(timed & 1)) {
        clocks[cell] = (double)NAN;
      }
      if (!(timed & 2)) {
        rates[cell] = (double)NAN;
      }
      flags[cell] = (unsigned char)myb_sp3_flags(sp3, s, e);
    }
  }
  const myb_sp3_data data = {.version = 'd',
                             .file_type = 'M',
                             .time_system = myb_sp3_time_system(sp3),
                             .frame = myb_sp3_frame(sp3),
                             .orbit_type = "FIT",
                             .agency = "MYB",
                             .data_used = "ORBIT",
                             .satellite_count = satellites,
                             .satellites = ids,
                             .epoch_count = epochs,
                             .epochs = labels,
                             .interval = 900.0,
                             .positions = positions,
                             .velocities = velocities,
                             .clocks = clocks,
                             .clock_rates = rates,
                             .flags = flags};
  size_t length = 0;
  const myb_status sized = myb_sp3_write(&data, NULL, 0, &length);
  if (sized == MYB_E_OVERFLOW) {
    char* text = malloc(length + 1);
    MYB_CHECK(text != NULL);
    MYB_CHECK(myb_sp3_write(&data, text, length + 1, &length) == MYB_OK);
    myb_sp3* back = NULL;
    myb_diag diag = {0};
    MYB_CHECK(myb_sp3_load(text, length, NULL, &back, &diag) == MYB_OK);
    MYB_CHECK(myb_sp3_satellite_count(back) == satellites && myb_sp3_epoch_count(back) == epochs);
    MYB_CHECK(strcmp(myb_sp3_time_system(back), data.time_system) == 0 &&
              strcmp(myb_sp3_frame(back), data.frame) == 0);
    for (size_t s = 0; s < satellites; ++s) {
      MYB_CHECK(strcmp(myb_sp3_satellite(back, s), ids[s]) == 0);
    }
    for (size_t e = 0; e < epochs; ++e) {
      myb_calendar label;
      myb_sp3_epoch(back, e, &label);
      MYB_CHECK(label.year == labels[e].year && label.month == labels[e].month &&
                label.day == labels[e].day && label.hour == labels[e].hour &&
                label.minute == labels[e].minute && within(label.second, labels[e].second, 1e-8));
      for (size_t s = 0; s < satellites; ++s) {
        const size_t cell = e * satellites + s;
        double r[3];
        double v[3];
        const int given = myb_sp3_state(back, s, e, r, v);
        MYB_CHECK(((given & 1) != 0) == isfinite(positions[3 * cell]));
        MYB_CHECK(((given & 2) != 0) == isfinite(velocities[3 * cell]));
        for (int i = 0; i < 3; ++i) {
          MYB_CHECK(!(given & 1) || within(r[i], positions[3 * cell + (size_t)i], 1e-3));
          MYB_CHECK(!(given & 2) || within(v[i], velocities[3 * cell + (size_t)i], 1e-7));
        }
        double clock = 0.0;
        double rate = 0.0;
        const int timed = myb_sp3_clock(back, s, e, &clock, &rate);
        MYB_CHECK(((timed & 1) != 0) == isfinite(clocks[cell]) &&
                  ((timed & 2) != 0) == isfinite(rates[cell]));
        MYB_CHECK(!(timed & 1) || within(clock, clocks[cell], 1e-12));
        MYB_CHECK(!(timed & 2) || within(rate, rates[cell], 1e-16));
        MYB_CHECK(myb_sp3_flags(back, s, e) == flags[cell]);
      }
    }
    myb_sp3_destroy(back);
    free(text);
  } else {
    MYB_CHECK(sized == MYB_E_DOMAIN);
  }
  free(flags);
  free(rates);
  free(clocks);
  free(velocities);
  free(positions);
  free(labels);
  free((void*)ids);
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  myb_sp3* sp3 = NULL;
  myb_diag diag = {0};
  if (myb_sp3_load((const char*)data, size, NULL, &sp3, &diag) != MYB_OK) {
    MYB_CHECK(sp3 == NULL);
    return 0;
  }
  const size_t satellites = myb_sp3_satellite_count(sp3);
  const size_t epochs = myb_sp3_epoch_count(sp3);
  MYB_CHECK(myb_sp3_time_system(sp3) != NULL && myb_sp3_frame(sp3) != NULL);
  size_t records = 0;
  for (size_t s = 0; s < satellites; ++s) {
    MYB_CHECK(myb_sp3_satellite(sp3, s) != NULL);
    (void)myb_sp3_accuracy(sp3, s);
  }
  for (size_t e = 0; e < epochs && records < MAX_RECORDS; ++e) {
    myb_calendar label;
    myb_sp3_epoch(sp3, e, &label);
    /* a label that was read: in the ranges of the format */
    MYB_CHECK(label.year >= 1900 && label.year <= 2999 && label.month >= 1 && label.month <= 12 &&
              label.day >= 1 && label.day <= 31 && label.second >= 0.0 && label.second < 61.0);
    for (size_t s = 0; s < satellites && records < MAX_RECORDS; ++s, ++records) {
      double position[3];
      double velocity[3];
      const int given = myb_sp3_state(sp3, s, e, position, velocity);
      MYB_CHECK(given >= 0 && given <= 3);
      double clock = 0.0;
      double rate = 0.0;
      const int clocks = myb_sp3_clock(sp3, s, e, &clock, &rate);
      MYB_CHECK(clocks >= 0 && clocks <= 3);
      double sigma[8];
      myb_sp3_sigma(sp3, s, e, sigma);
      for (int i = 0; i < 8; ++i) {
        MYB_CHECK(sigma[i] >= 0.0 && isfinite(sigma[i]));
      }
      (void)myb_sp3_flags(sp3, s, e);
    }
  }
  round_trip(sp3);
  myb_sp3_destroy(sp3);
  return 0;
}
