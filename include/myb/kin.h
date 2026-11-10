/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/kin.h — kinematic providers: states of bodies given as functions of time.
 *
 * SPK ephemerides. myb reads NAIF SPK files (planetary ephemerides DE440/DE441 of JPL and
 * EPM of IAA RAS, small-body files of JPL Horizons, spacecraft trajectories) with its own
 * reader, after "SPK Required Reading" and "DAF Required Reading" (NAIF, JPL). A segment gives
 * the state of a target relative to a center over a time interval, in a reference frame named
 * by a NAIF frame code (1 = J2000, i.e. ICRF for the JPL ephemerides). The identifiers are the
 * file's: the library attaches no meaning to them.
 *
 * Segment types read: 1 and 21 (modified difference arrays), 2 and 3 (Chebyshev, position or
 * position and velocity), 5 (discrete states with two-body propagation), 8 and 9 (Lagrange),
 * 12 and 13 (Hermite), 14 (Chebyshev, unequal steps), 17 (precessing equinoctial elements),
 * 18 and 19 (ESOC/DDID Hermite and Lagrange), 20 (Chebyshev velocity, integrated position).
 * Other types load, and a query that needs one of their segments fails with MYB_E_FORMAT.
 *
 * Time is the ephemeris' argument, TDB for the JPL and IAA ephemerides (T_eph), as a myb_time
 * on that scale: its two parts keep the full resolution in the Chebyshev and interpolation
 * arguments, which a single double of seconds past J2000 (resolution 0.1 µs today) does not.
 * Results are in m and m/s. */
#ifndef MYB_KIN_H_
#define MYB_KIN_H_

#include <stddef.h>
#include <stdint.h>

#include "myb/core.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct myb_spk myb_spk;

/* An empty set of SPK files. allocator may be NULL (malloc/free). */
MYB_API MYB_NODISCARD myb_status myb_spk_create(const myb_allocator* allocator, myb_spk** spk);

/* Unmaps the files; accepts NULL. */
MYB_API void myb_spk_destroy(myb_spk* spk);

/* Adds a file (path in UTF-8), mapped into memory and read on demand. Later files take
 * precedence over earlier ones, and within a file later segments over earlier ones (the SPICE
 * rule). Loading is atomic: on error nothing is added and diag tells why. Errors: MYB_E_IO,
 * MYB_E_FORMAT (not a DAF/SPK file, an inconsistent segment), MYB_E_NOMEM, MYB_E_OVERFLOW.
 * Not to be called while another thread queries the set. */
MYB_API MYB_NODISCARD myb_status myb_spk_load_file(myb_spk* spk, const char* path, myb_diag* diag);

/* The same for an SPK image in memory, which the caller keeps alive and unchanged until the
 * set is destroyed. */
MYB_API MYB_NODISCARD myb_status myb_spk_load_memory(myb_spk* spk, const void* data, size_t size,
                                                     myb_diag* diag);

typedef struct myb_spk_segment_info {
  int32_t target;
  int32_t center;
  int32_t frame;     /* NAIF frame code */
  int32_t type;      /* SPK data type */
  double start;      /* coverage, seconds past J2000 on the ephemeris time scale */
  double stop;
  size_t file;       /* 0-based load order of the file */
  char name[41];     /* segment identifier, trailing blanks removed */
} myb_spk_segment_info;

/* Segments in load order (file by file, in file order). */
MYB_API size_t myb_spk_segment_count(const myb_spk* spk);
MYB_API void myb_spk_segment(const myb_spk* spk, size_t index, myb_spk_segment_info* info);

/* Number of double words of a segment, and count of them from 0-based offset (in the native
 * byte order) — for evaluators of segment types outside the library (SPK type 10, two-line
 * elements, in ext/sgp4). Errors: MYB_E_RANGE (offset + count beyond the segment). */
MYB_API size_t myb_spk_segment_size(const myb_spk* spk, size_t index);
MYB_API MYB_NODISCARD myb_status myb_spk_segment_data(const myb_spk* spk, size_t index,
                                                      size_t offset, size_t count,
                                                      double* values);

/* State (position m, velocity m/s) of target relative to center at t. The two bodies are
 * connected through the segments that cover t (each body's highest-priority covering segment
 * gives it relative to its center, up to a common body), and the states along the two branches
 * are summed with compensated summation. *frame (may be NULL) receives the frame code, which
 * all those segments must share. target == center gives zeros (frame 0 if no data).
 * Errors: MYB_E_RANGE (no chain of segments covering t joins the two bodies),
 * MYB_E_FORMAT (segments in different frames, an unsupported segment type, or data that turn
 * out invalid), MYB_E_DOMAIN (two-body propagation of type 5 failed). Thread-safe: queries do
 * not modify the set. */
MYB_API MYB_NODISCARD myb_status myb_spk_state(const myb_spk* spk, int32_t target, int32_t center,
                                               myb_time t, double state[6], int32_t* frame);

/* ---- Providers -----------------------------------------------------------------------------
 * A provider gives the state (position m, velocity m/s) of a point as a function of time, in
 * the frame and on the time scale agreed between it and its user: the form in which
 * ephemerides enter the light-time equation (myb/observe.h) and the environment of a force
 * model (myb/env.h). Any status other than MYB_OK is returned to the caller unchanged. */
typedef myb_status (*myb_state_function)(void* context, myb_time t, double state[6]);

/* A body of an SPK set as a provider: target relative to center at t on the ephemeris' time
 * scale (myb_spk_state). frame != 0 requires the segments to be in that frame (NAIF code,
 * 1 = J2000/ICRF) and fails with MYB_E_FORMAT otherwise; 0 accepts any. */
typedef struct myb_spk_body {
  const myb_spk* spk;
  int32_t target;
  int32_t center;
  int32_t frame;
} myb_spk_body;

/* myb_state_function; context: const myb_spk_body*. */
MYB_API MYB_NODISCARD myb_status myb_spk_body_state(void* body, myb_time t, double state[6]);

/* ---- Writing SPK files -----------------------------------------------------------------------
 * A writer collects segments in memory and lays them out as one little-endian SPK file, which
 * SPICE and myb_spk_load_file read. States are taken in m and m/s (the file holds km and
 * km/s: values pass through one division by 1000 on the way in and one multiplication on the
 * way out, so a state read back agrees with the one written to an ulp), epochs in seconds
 * past J2000 on the ephemeris' time scale (TDB). */

typedef struct myb_spk_writer myb_spk_writer;

/* internal_name: the label of the file, at most 60 printable ASCII characters (NULL: "myb").
 * Errors: MYB_E_DOMAIN (the name), MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_spk_writer_create(const char* internal_name,
                                                       const myb_allocator* allocator,
                                                       myb_spk_writer** writer);
MYB_API void myb_spk_writer_destroy(myb_spk_writer* writer); /* accepts NULL */

/* What a segment describes: target relative to center in the frame (NAIF codes, 1 = J2000);
 * name: its identifier, 1 .. 40 printable ASCII characters. */
typedef struct myb_spk_segment_header {
  int32_t target;
  int32_t center;
  int32_t frame;
  const char* name;
} myb_spk_segment_header;

/* A segment of count >= 2 states at increasing epochs, interpolated by Lagrange polynomials of
 * the given degree through degree + 1 states (type 9), or by Hermite polynomials of the given
 * odd degree through (degree + 1) / 2 states with their velocities (type 13). degree: 1 .. 27.
 * At the epochs themselves both give the states back. Errors: MYB_E_DOMAIN (the header, the
 * type, the degree, fewer states than the window, epochs not increasing, values not finite),
 * MYB_E_NOMEM, MYB_E_OVERFLOW. */
MYB_API MYB_NODISCARD myb_status myb_spk_writer_add_states(myb_spk_writer* writer,
                                                           const myb_spk_segment_header* header,
                                                           int type, int degree,
                                                           const double* epochs,
                                                           const double* states, size_t count);

/* A segment of Chebyshev polynomials of the given degree (1 .. 27) on `intervals` equal parts
 * of [start, stop], fitted to a provider (called with the epochs as myb_time of seconds past
 * J2000): the polynomials through the states at the degree + 1 Chebyshev nodes of each part.
 * Type 2 keeps the position (readers differentiate it for the velocity), type 3 the position
 * and the velocity. *max_error (may be NULL; asking costs degree + 2 more calls of the
 * provider per part): the largest distance (m) between the provider and the polynomials
 * half-way between the nodes — to choose the degree and the number of parts. Errors:
 * MYB_E_DOMAIN (as above, stop <= start, a state not finite), the statuses of the provider. */
MYB_API MYB_NODISCARD myb_status myb_spk_writer_add_chebyshev(
    myb_spk_writer* writer, const myb_spk_segment_header* header, int type,
    myb_state_function state, void* context, double start, double stop, size_t intervals,
    int degree, double* max_error);

/* A segment of any type from its words as a file holds them (myb_spk_segment_data gives
 * them): copies, parts and merges of files, and the types that are not built here — 1 and
 * 21, the difference arrays of an integrator, which only that integrator produces; 10,
 * two-line elements (myb_sgp4_spk10_words of ext/sgp4 lays them out); 5, 8, 12, 14, 17 - 20.
 * The words are written as given: nothing is converted or checked against the type.
 * Errors: MYB_E_DOMAIN (the header, type < 1, stop < start, no words), MYB_E_NOMEM,
 * MYB_E_OVERFLOW. */
MYB_API MYB_NODISCARD myb_status myb_spk_writer_add_raw(myb_spk_writer* writer,
                                                        const myb_spk_segment_header* header,
                                                        int type, double start, double stop,
                                                        const double* words, size_t count);

/* Bytes of the file with the segments added so far (0 without segments). */
MYB_API size_t myb_spk_writer_size(const myb_spk_writer* writer);

/* The file into a buffer of at least myb_spk_writer_size bytes. Errors: MYB_E_DOMAIN (no
 * segments), MYB_E_OVERFLOW (the buffer, or more than 2^31 words of data). */
MYB_API MYB_NODISCARD myb_status myb_spk_writer_write(const myb_spk_writer* writer, void* buffer,
                                                      size_t size);

/* The file at path (UTF-8), created or replaced. Additional errors: MYB_E_IO, MYB_E_NOMEM. */
MYB_API MYB_NODISCARD myb_status myb_spk_writer_write_file(const myb_spk_writer* writer,
                                                           const char* path);

/* ---- SP3 precise orbits ---------------------------------------------------------------------
 * SP3-c and SP3-d files (IGS, ILRS analysis centres): satellite positions (and velocities) at
 * tabulated epochs, usually in an Earth-fixed frame (ITRF, SLRF) named by the header. Epochs are
 * calendar labels in the file's time system ("GPS", "UTC", "TAI", ...), which the caller
 * converts (myb_timesys_from_utc for UTC). Positions in m, velocities in m/s; positions given
 * as 0, 0, 0 (bad or missing) are reported missing. */

typedef struct myb_sp3 myb_sp3;

/* Errors: MYB_E_FORMAT (diag gives the line), MYB_E_NOMEM, MYB_E_OVERFLOW, MYB_E_IO (files). */
MYB_API MYB_NODISCARD myb_status myb_sp3_load(const char* text, size_t size,
                                              const myb_allocator* allocator, myb_sp3** sp3,
                                              myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_sp3_load_file(const char* path,
                                                   const myb_allocator* allocator,
                                                   myb_sp3** sp3, myb_diag* diag);
MYB_API void myb_sp3_destroy(myb_sp3* sp3); /* accepts NULL */

MYB_API size_t myb_sp3_satellite_count(const myb_sp3* sp3);
MYB_API const char* myb_sp3_satellite(const myb_sp3* sp3, size_t index); /* e.g. "L51" */
MYB_API size_t myb_sp3_epoch_count(const myb_sp3* sp3);
MYB_API const char* myb_sp3_time_system(const myb_sp3* sp3);
MYB_API const char* myb_sp3_frame(const myb_sp3* sp3); /* coordinate system, e.g. "SLR20" */
MYB_API void myb_sp3_epoch(const myb_sp3* sp3, size_t epoch, myb_calendar* label);

/* Position and velocity of a satellite at an epoch (either output may be NULL); returns the
 * flags of what the file gives: 1 = position, 2 = velocity. */
MYB_API int myb_sp3_state(const myb_sp3* sp3, size_t satellite, size_t epoch, double position[3],
                          double velocity[3]);

/* The clock correction of a satellite at an epoch (s; the file keeps microseconds) and its
 * rate (s/s), either output may be NULL; returns what the file gives: 1 = the clock, 2 = the
 * rate. */
MYB_API int myb_sp3_clock(const myb_sp3* sp3, size_t satellite, size_t epoch, double* clock,
                          double* rate);

/* The standard deviations the records of a satellite at an epoch carry: [0..2] of the
 * position (m), [3] of the clock (s), [4..6] of the velocity (m/s), [7] of the clock rate
 * (s/s); 0 where the file gives none. The file keeps them as whole powers of the two bases
 * of its header (1.25 and 1.025 by habit): a value read is within 12 % (1.2 % for clocks)
 * of the one the producer had. */
MYB_API void myb_sp3_sigma(const myb_sp3* sp3, size_t satellite, size_t epoch, double sigma[8]);

/* The accuracy of the orbit of a satellite that the header states (m; a power of two of
 * millimetres), 0 if unknown. */
MYB_API double myb_sp3_accuracy(const myb_sp3* sp3, size_t satellite);

/* The flags of the position record of a satellite at an epoch. */
enum {
  MYB_SP3_CLOCK_EVENT = 1,     /* a discontinuity of the clock since the last epoch */
  MYB_SP3_CLOCK_PREDICTED = 2, /* the clock is a prediction */
  MYB_SP3_MANEUVER = 4,        /* a manoeuvre since the last epoch */
  MYB_SP3_ORBIT_PREDICTED = 8  /* the position is a prediction */
};
MYB_API unsigned myb_sp3_flags(const myb_sp3* sp3, size_t satellite, size_t epoch);

/* What an SP3 file is written from: positions[3 (epoch * satellite_count + satellite)] in m
 * (a value that is not finite marks a missing position), velocities likewise in m/s or NULL.
 * The format keeps six decimals of km and of dm/s: a millimetre and 0.1 micrometre per
 * second; of microseconds for the clocks (a picosecond; |clock| < 1 s). The members after
 * `velocities` may all be NULL: no clocks, no accuracies, no flags. */
typedef struct myb_sp3_data {
  char version;            /* 'c' (at most 85 satellites) or 'd' */
  char file_type;          /* 'G' GPS, 'M' mixed, 'L' low orbiters, ... */
  const char* time_system; /* 3 characters: "GPS", "UTC", "TAI", ... */
  const char* frame;       /* at most 5: "IGS20", "SLR20", ... */
  const char* orbit_type;  /* 3: "FIT", "EXT", "BCT", "HLM" */
  const char* agency;      /* at most 4 */
  const char* data_used;   /* at most 5: "ORBIT", "SLR", ... */
  size_t satellite_count;
  const char* const* satellites; /* identifiers of 3 characters: "G01", "L51" */
  size_t epoch_count;
  const myb_calendar* epochs; /* labels in the time system */
  double interval;            /* s between epochs, > 0 */
  const double* positions;
  const double* velocities;
  const double* clocks;      /* [epoch * satellite_count + satellite], s; not finite: none */
  const double* clock_rates; /* likewise, s/s; written with the velocities */
  const double* accuracy;    /* [satellite], m; 0: unknown */
  const double* sigmas;      /* 8 a record, as myb_sp3_sigma gives them; 0: none */
  const unsigned char* flags; /* a record: MYB_SP3_* */
} myb_sp3_data;

/* The file as text. *length (may be NULL): its characters without the NUL. Errors:
 * MYB_E_OVERFLOW (the buffer; call again with length + 1 bytes), MYB_E_DOMAIN (an invalid
 * description). */
MYB_API MYB_NODISCARD myb_status myb_sp3_write(const myb_sp3_data* data, char* buffer,
                                               size_t size, size_t* length);

#ifdef __cplusplus
}
#endif

#endif /* MYB_KIN_H_ */
