/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/ck.h — the orientation of spacecraft and their instruments from NAIF CK files, and the
 * clocks they are counted by.
 *
 * A CK file ("CK Required Reading", NAIF, JPL) holds, for an instrument code, the rotation
 * from a base frame to the instrument's frame — the C-matrix — as quaternions at instants of
 * the spacecraft clock, optionally with the angular velocity of the frame. The clock counts
 * ticks; a clock kernel ("SCLK Required Reading") ties them to a time scale by pieces of
 * constant rate. The codes are the files': the library attaches no meaning to them.
 *
 * Quaternions are SPICE's: q = (cos(a/2), sin(a/2) u) for the matrix that turns vectors by
 * the angle a about u, which here is the C-matrix. Angular velocities are those of the
 * instrument frame relative to the base frame, in the base frame, rad/s.
 *
 * Segment types read: 1 (pointing at discrete instants), 2 (intervals of constant angular
 * velocity), 3 (interpolation between instants inside intervals). Other types load, and a
 * query that would need one of their segments passes it over. */
#ifndef MYB_CK_H_
#define MYB_CK_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"
#include "myb/frames.h"
#include "myb/io.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Spacecraft clocks -------------------------------------------------------------------------
 * A clock of type 1: a count of several fields (each with a modulus and an offset, the last
 * one the tick), partitions — spans between resets of the count — and a table of (ticks,
 * time, rate): from each row on the time runs at `rate` seconds per count of the most
 * significant field. Ticks are counted through the partitions without a break ("encoded"),
 * as the files of the orientation give them. */

enum { MYB_SCLK_MAX_FIELDS = 10 };

typedef struct myb_sclk myb_sclk;

/* The clock whose kernel variables end in _<id> (the spacecraft code without its sign):
 * SCLK_DATA_TYPE, SCLK01_N_FIELDS, SCLK01_MODULI, SCLK01_OFFSETS, SCLK01_OUTPUT_DELIM,
 * SCLK_PARTITION_START, SCLK_PARTITION_END, SCLK01_COEFFICIENTS and, optionally,
 * SCLK01_TIME_SYSTEM. What is needed is copied. Errors: MYB_E_NOT_FOUND (a variable),
 * MYB_E_FORMAT (a type other than 1, inconsistent values; diag tells which), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_sclk_create(const myb_kpool* pool, int64_t id,
                                                 const myb_allocator* allocator, myb_sclk** sclk,
                                                 myb_diag* diag);
MYB_API void myb_sclk_destroy(myb_sclk* sclk); /* accepts NULL */

/* The scale of the times of the table: 1 = TDB, 2 = TT. */
MYB_API int myb_sclk_time_system(const myb_sclk* sclk);
MYB_API size_t myb_sclk_partitions(const myb_sclk* sclk);

/* Encoded ticks -> the instant on the clock's time scale, and back (a table whose times
 * are not continuous has no single inverse: the row is chosen by the time, as SPICE does).
 * Errors: MYB_E_RANGE (before the first row of the table or after the end of the last
 * partition). */
MYB_API MYB_NODISCARD myb_status myb_sclk_to_time(const myb_sclk* sclk, double ticks,
                                                  myb_time* t);
MYB_API MYB_NODISCARD myb_status myb_sclk_from_time(const myb_sclk* sclk, myb_time t,
                                                    double* ticks);

/* A clock string -> encoded ticks: "partition/field field ..." with the fields separated by
 * any of . : - , or blanks, whole numbers, missing fields on the right taken at their
 * offsets; without a partition the first one that holds the count is taken. Errors:
 * MYB_E_FORMAT (not such a string), MYB_E_DOMAIN (no such partition, a count outside it). */
MYB_API MYB_NODISCARD myb_status myb_sclk_encode(const myb_sclk* sclk, const char* text,
                                                 size_t length, double* ticks);

/* Encoded ticks (rounded to a whole tick) -> the string with the partition, every field at
 * the width of its largest value. *length (may be NULL): its characters without the NUL.
 * Errors: MYB_E_RANGE (ticks outside the partitions), MYB_E_OVERFLOW (the buffer). */
MYB_API MYB_NODISCARD myb_status myb_sclk_decode(const myb_sclk* sclk, double ticks, char* text,
                                                 size_t size, size_t* length);

/* ---- Orientation files ----------------------------------------------------------------------- */

typedef struct myb_ck myb_ck;

MYB_API MYB_NODISCARD myb_status myb_ck_create(const myb_allocator* allocator, myb_ck** ck);
MYB_API void myb_ck_destroy(myb_ck* ck); /* accepts NULL */

/* As myb_spk_load_file / myb_spk_load_memory (myb/kin.h): later files take precedence over
 * earlier ones, later segments of a file over earlier ones. */
MYB_API MYB_NODISCARD myb_status myb_ck_load_file(myb_ck* ck, const char* path, myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_ck_load_memory(myb_ck* ck, const void* data, size_t size,
                                                    myb_diag* diag);

typedef struct myb_ck_segment_info {
  int32_t instrument;
  int32_t frame;  /* the base frame (NAIF code) */
  int32_t type;
  int has_rates;  /* 1: the segment gives angular velocities */
  double start;   /* coverage, encoded ticks */
  double stop;
  size_t file;    /* 0-based load order of the file */
  char name[41];
} myb_ck_segment_info;

MYB_API size_t myb_ck_segment_count(const myb_ck* ck);
MYB_API void myb_ck_segment(const myb_ck* ck, size_t index, myb_ck_segment_info* info);

typedef struct myb_ck_pointing {
  double q[4];           /* the quaternion of the C-matrix */
  double rate[3];        /* the angular velocity, in the base frame (rad/s); zeros without */
  int has_rate;
  double ticks;          /* the instant of the clock the pointing is for */
  int32_t frame;         /* the base frame */
  myb_rotation rotation; /* base frame -> instrument frame; its derivative from the angular
                            velocity (zeros without one) */
} myb_ck_pointing;

/* The orientation of an instrument at `ticks` of its clock, or the nearest that the files
 * hold within `tolerance` ticks (the rule of SPICE's ckgpav_c): the segments of the
 * instrument are tried in the order of their precedence, and the first that gives a pointing
 * within the tolerance is taken — inside an interval of types 2 and 3 at the instant asked,
 * elsewhere at the nearest instant of the segment. need_rate != 0 passes over segments
 * without angular velocities. Errors: MYB_E_RANGE (no segment gives one), MYB_E_DOMAIN (a
 * negative or non-finite tolerance), MYB_E_FORMAT (data that turn out invalid).
 * Thread-safe. */
MYB_API MYB_NODISCARD myb_status myb_ck_pointing_at(const myb_ck* ck, int32_t instrument,
                                                    double ticks, double tolerance,
                                                    int need_rate, myb_ck_pointing* pointing);

/* An instrument frame as an orientation provider (myb_orientation_function; context: const
 * myb_ck_orientation*): the argument is the instant on the clock's time scale
 * (myb_sclk_time_system — the scale to register the provider with in a graph of frames).
 * Errors: those of myb_sclk_from_time and myb_ck_pointing_at; MYB_E_FORMAT for a pointing
 * given in another frame than `base`. */
typedef struct myb_ck_orientation {
  const myb_ck* ck;
  const myb_sclk* sclk;
  int32_t instrument;
  int32_t base;     /* the frame the pointing must be given in; 0: any */
  double tolerance; /* ticks; 0 asks for the instant itself */
} myb_ck_orientation;

MYB_API MYB_NODISCARD myb_status myb_orientation_ck(void* context, myb_time t,
                                                    myb_orientation* orientation);

/* ---- Writing orientation files ------------------------------------------------------------------
 * As the SPK writer of myb/kin.h: segments are collected in memory and laid out as one
 * little-endian file that SPICE and myb_ck_load_file read. */

typedef struct myb_ck_writer myb_ck_writer;

MYB_API MYB_NODISCARD myb_status myb_ck_writer_create(const char* internal_name,
                                                      const myb_allocator* allocator,
                                                      myb_ck_writer** writer);
MYB_API void myb_ck_writer_destroy(myb_ck_writer* writer); /* accepts NULL */

typedef struct myb_ck_segment_header {
  int32_t instrument;
  int32_t frame;    /* the base frame (NAIF code, 1 = J2000) */
  const char* name; /* 1 .. 40 printable ASCII characters */
} myb_ck_segment_header;

/* A segment of `count` pointings at increasing ticks: quaternions (4 a pointing; made unit
 * on the way), rates (3 a pointing, rad/s; NULL: none). type 1 gives them back at their
 * instants only; type 3 interpolates inside the intervals that begin at the ticks
 * `interval_starts` (increasing, each one of `ticks`, the first of them the first tick;
 * NULL: one interval of them all). Errors: MYB_E_DOMAIN (the header, the type, ticks not
 * increasing, a zero quaternion, values not finite, interval starts that are not among the
 * ticks), MYB_E_NOMEM, MYB_E_OVERFLOW. */
MYB_API MYB_NODISCARD myb_status myb_ck_writer_add_pointings(
    myb_ck_writer* writer, const myb_ck_segment_header* header, int type, const double* ticks,
    const double* quaternions, const double* rates, size_t count,
    const double* interval_starts, size_t interval_count);

/* A segment of type 2: `count` intervals [starts[k], stops[k]] that follow one another
 * without overlap, each with the quaternion at its start, a constant angular velocity (3,
 * rad/s) and the seconds a tick lasts in it. Errors as above. */
MYB_API MYB_NODISCARD myb_status myb_ck_writer_add_intervals(
    myb_ck_writer* writer, const myb_ck_segment_header* header, const double* starts,
    const double* stops, const double* quaternions, const double* rates,
    const double* seconds_per_tick, size_t count);

MYB_API size_t myb_ck_writer_size(const myb_ck_writer* writer);
MYB_API MYB_NODISCARD myb_status myb_ck_writer_write(const myb_ck_writer* writer, void* buffer,
                                                     size_t size);
MYB_API MYB_NODISCARD myb_status myb_ck_writer_write_file(const myb_ck_writer* writer,
                                                          const char* path);

#ifdef __cplusplus
}
#endif

#endif /* MYB_CK_H_ */
