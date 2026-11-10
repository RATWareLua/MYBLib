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

#ifdef __cplusplus
}
#endif

#endif /* MYB_KIN_H_ */
