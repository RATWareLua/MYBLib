/* SPDX-License-Identifier: BSD-3-Clause */
/* glue_files.c — orbit data messages, SP3, CK and the clocks of SCLK kernels for the Go
 * binding (glue_tasks.h). */
#include <stdlib.h>
#include <string.h>

#include "glue_internal.h"

/* The mark of a missing record. */
static const double missing = (double)__builtin_nan("");

/* A writer that takes a buffer and tells the length it needs when the buffer is short. */
typedef myb_status (*text_writer)(const void* message, myb_odm_notation notation, char* buffer,
                                  size_t size, size_t* length);

static myb_status write_opm(const void* message, myb_odm_notation notation, char* buffer,
                            size_t size, size_t* length) {
  return myb_opm_write_as(message, notation, buffer, size, length);
}

static myb_status write_omm(const void* message, myb_odm_notation notation, char* buffer,
                            size_t size, size_t* length) {
  return myb_omm_write_as(message, notation, buffer, size, length);
}

static myb_status write_text(text_writer writer, const void* message, int notation, char** text,
                             size_t* length) {
  *text = NULL;
  *length = 0;
  const myb_odm_notation as = notation == MYB_ODM_XML ? MYB_ODM_XML : MYB_ODM_KVN;
  size_t needed = 0;
  myb_status status = writer(message, as, NULL, 0, &needed);
  if (status != MYB_E_OVERFLOW) {
    return status == MYB_OK ? MYB_E_DOMAIN : status;
  }
  char* buffer = malloc(needed + 1);
  if (buffer == NULL) {
    return MYB_E_NOMEM;
  }
  status = writer(message, as, buffer, needed + 1, length);
  if (status == MYB_OK) {
    *text = buffer;
  } else {
    free(buffer);
  }
  return status;
}

myb_status mybgo_opm_write(const myb_opm* opm, int notation, char** text, size_t* length) {
  return write_text(write_opm, opm, notation, text, length);
}

myb_status mybgo_omm_write(const myb_omm* omm, int notation, char** text, size_t* length) {
  return write_text(write_omm, omm, notation, text, length);
}

myb_status mybgo_sp3_write(const myb_sp3_data* head, const char* const* satellites,
                           const myb_calendar* epochs, const double* positions,
                           const double* velocities, const double* clocks,
                           const double* clock_rates, const double* accuracy, const double* sigmas,
                           const unsigned char* flags, char** text, size_t* length) {
  *text = NULL;
  *length = 0;
  myb_sp3_data data = *head;
  data.satellites = satellites;
  data.epochs = epochs;
  data.positions = positions;
  data.velocities = velocities;
  data.clocks = clocks;
  data.clock_rates = clock_rates;
  data.accuracy = accuracy;
  data.sigmas = sigmas;
  data.flags = flags;
  size_t needed = 0;
  myb_status status = myb_sp3_write(&data, NULL, 0, &needed);
  if (status != MYB_E_OVERFLOW) {
    return status == MYB_OK ? MYB_E_DOMAIN : status;
  }
  char* buffer = malloc(needed + 1);
  if (buffer == NULL) {
    return MYB_E_NOMEM;
  }
  status = myb_sp3_write(&data, buffer, needed + 1, length);
  if (status == MYB_OK) {
    *text = buffer;
  } else {
    free(buffer);
  }
  return status;
}

void mybgo_sp3_records(const myb_sp3* sp3, double* positions, double* velocities, int* given,
                       double* clocks, double* clock_rates, double* sigmas, unsigned char* flags) {
  const size_t satellites = myb_sp3_satellite_count(sp3);
  const size_t epochs = myb_sp3_epoch_count(sp3);
  for (size_t e = 0; e < epochs; ++e) {
    for (size_t s = 0; s < satellites; ++s) {
      const size_t k = e * satellites + s;
      double p[3] = {missing, missing, missing};
      double v[3] = {missing, missing, missing};
      given[k] = myb_sp3_state(sp3, s, e, p, v);
      if ((given[k] & 1) == 0) {
        p[0] = p[1] = p[2] = missing;
      }
      if ((given[k] & 2) == 0) {
        v[0] = v[1] = v[2] = missing;
      }
      memcpy(positions + 3 * k, p, sizeof p);
      memcpy(velocities + 3 * k, v, sizeof v);
      double clock = missing;
      double rate = missing;
      const int which = myb_sp3_clock(sp3, s, e, &clock, &rate);
      clocks[k] = (which & 1) != 0 ? clock : missing;
      clock_rates[k] = (which & 2) != 0 ? rate : missing;
      myb_sp3_sigma(sp3, s, e, sigmas + 8 * k);
      flags[k] = (unsigned char)myb_sp3_flags(sp3, s, e);
    }
  }
}

myb_status mybgo_ck_write(const char* path, const char* internal_name,
                          const mybgo_ck_segment* segments, size_t segment_count,
                          const double* ticks, const double* stops, const double* quaternions,
                          const double* rates, const double* seconds_per_tick,
                          const double* interval_starts, void** data, size_t* size,
                          size_t* failed) {
  *failed = segment_count;
  if (path == NULL && (data == NULL || size == NULL)) {
    return MYB_E_DOMAIN;
  }
  if (data != NULL) {
    *data = NULL;
    *size = 0;
  }
  myb_ck_writer* writer = NULL;
  myb_status status = myb_ck_writer_create(internal_name, NULL, &writer);
  size_t value = 0;    /* the first value of the segment */
  size_t interval = 0; /* the first interval start, the first seconds per tick */
  for (size_t j = 0; j < segment_count && status == MYB_OK; ++j) {
    const mybgo_ck_segment* s = &segments[j];
    if (s->kind == MYBGO_CK_POINTINGS) {
      status = myb_ck_writer_add_pointings(
          writer, &s->header, s->type, ticks + value, quaternions + 4 * value,
          s->has_rates ? rates + 3 * value : NULL, s->count,
          s->interval_count > 0 ? interval_starts + interval : NULL, s->interval_count);
      interval += s->interval_count;
    } else if (s->kind == MYBGO_CK_INTERVALS) {
      status = myb_ck_writer_add_intervals(writer, &s->header, ticks + value, stops + value,
                                           quaternions + 4 * value, rates + 3 * value,
                                           seconds_per_tick + value, s->count);
    } else {
      status = MYB_E_DOMAIN;
    }
    if (status != MYB_OK) {
      *failed = j;
    }
    value += s->count;
  }
  if (status == MYB_OK && path != NULL) {
    status = myb_ck_writer_write_file(writer, path);
  } else if (status == MYB_OK) {
    const size_t bytes = myb_ck_writer_size(writer);
    void* buffer = malloc(bytes > 0 ? bytes : 1);
    status = buffer != NULL ? myb_ck_writer_write(writer, buffer, bytes) : MYB_E_NOMEM;
    if (status == MYB_OK) {
      *data = buffer;
      *size = bytes;
    } else {
      free(buffer);
    }
  }
  myb_ck_writer_destroy(writer);
  return status;
}

myb_status mybgo_ck_pointings(const myb_ck* ck, int32_t instrument, const double* ticks,
                              size_t count, double tolerance, int need_rate,
                              myb_ck_pointing* pointings, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_ck_pointing_at(ck, instrument, ticks[k], tolerance, need_rate, &pointings[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_sclk_times(const myb_sclk* sclk, const double* ticks, size_t count, myb_time* t,
                            size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_sclk_to_time(sclk, ticks[k], &t[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_sclk_ticks(const myb_sclk* sclk, const myb_time* t, size_t count, double* ticks,
                            size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_sclk_from_time(sclk, t[k], &ticks[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}
