/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/ccsds.h — orbit data messages of the CCSDS: the ephemeris message OEM, the parameter
 * message OPM, the mean-elements message OMM and the comprehensive message OCM, read and
 * written ("Orbit Data Messages", CCSDS 502.0-B-2/B-3) in both notations — keyword = value
 * (KVN) and XML ("XML Specification for Navigation Data Messages", CCSDS 505.0-B). The
 * readers tell the notation by the first character of the text ('<' for XML); the writers
 * take it as an argument.
 *
 * Epochs are calendar labels (YYYY-MM-DDThh:mm:ss.sss or YYYY-DDDThh:mm:ss.sss) on the time
 * system the message names (TIME_SYSTEM: UTC, TAI, TT, TDB, GPS, ...): the caller converts
 * them (myb_timesys_from_utc for UTC, myb_calendar_to_time for the continuous scales).
 * Covariances are 6 x 6 matrices of a state, by rows, of which the messages carry the lower
 * triangle. Comments and, in OEM, OPM and OMM, user-defined keywords are passed over. Of XML
 * the elements and attributes of the messages are read; CDATA sections and document types
 * with an internal subset are not. */
#ifndef MYB_CCSDS_H_
#define MYB_CCSDS_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/time.h"

#ifdef __cplusplus
extern "C" {
#endif

enum { MYB_ODM_TEXT = 64, MYB_OPM_MAX_MANEUVERS = 16 };

typedef enum myb_odm_notation {
  MYB_ODM_KVN = 0, /* KEYWORD = value lines */
  MYB_ODM_XML = 1
} myb_odm_notation;

/* ---- OEM: ephemeris -------------------------------------------------------------------------
 * A header and one or more segments, each a block of metadata (object, centre, reference
 * frame, time system, ...) and the states at tabulated epochs, in the file in km and km/s
 * (accelerations, if present, in km/s^2), here in m, m/s and m/s^2. */

typedef struct myb_oem myb_oem;

/* Errors: MYB_E_FORMAT (diag gives the line), MYB_E_NOMEM, MYB_E_OVERFLOW; MYB_E_IO (files). */
MYB_API MYB_NODISCARD myb_status myb_oem_load(const char* text, size_t size,
                                              const myb_allocator* allocator, myb_oem** oem,
                                              myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_oem_load_file(const char* path,
                                                   const myb_allocator* allocator, myb_oem** oem,
                                                   myb_diag* diag);
MYB_API void myb_oem_destroy(myb_oem* oem); /* accepts NULL */

/* The value of a header keyword ("CCSDS_OEM_VERS", "CREATION_DATE", "ORIGINATOR", ...) or of
 * a metadata keyword of a segment ("OBJECT_NAME", "REF_FRAME", "TIME_SYSTEM", "START_TIME",
 * "INTERPOLATION", ...), as written; NULL if the message does not have it. */
MYB_API const char* myb_oem_header(const myb_oem* oem, const char* key);
MYB_API size_t myb_oem_segment_count(const myb_oem* oem);
MYB_API const char* myb_oem_meta(const myb_oem* oem, size_t segment, const char* key);

MYB_API size_t myb_oem_state_count(const myb_oem* oem, size_t segment);

/* The epoch, the state and the acceleration (any output may be NULL) of a line of a segment;
 * returns 1, or 3 when the line also gives the acceleration (zeros otherwise). */
MYB_API int myb_oem_state(const myb_oem* oem, size_t segment, size_t index,
                          myb_calendar* epoch, double state[6], double acceleration[3]);

/* The covariance matrices of a segment: how many, and the epoch, the reference frame
 * (COV_REF_FRAME; NULL if the message names none: then that of the segment) and the matrix
 * (36 numbers by rows: m^2, m^2/s, m^2/s^2) of one; any output may be NULL. */
MYB_API size_t myb_oem_covariance_count(const myb_oem* oem, size_t segment);
MYB_API void myb_oem_covariance(const myb_oem* oem, size_t segment, size_t index,
                                myb_calendar* epoch, const char** ref_frame,
                                double covariance[36]);

/* What a writer puts into a segment; interpolation may be NULL (then no interpolation
 * keywords are written), and so the covariances (covariance_count 0). */
typedef struct myb_oem_segment_data {
  const char* object_name;
  const char* object_id;
  const char* center_name;
  const char* ref_frame;
  const char* time_system;
  const char* interpolation; /* e.g. "HERMITE", "LAGRANGE" */
  int interpolation_degree;
  size_t count;
  const myb_calendar* epochs;
  const double* states; /* 6 per epoch: m, m/s */
  size_t covariance_count;
  const myb_calendar* covariance_epochs;
  const char* covariance_frame; /* COV_REF_FRAME of each, or NULL */
  const double* covariances;    /* 36 per epoch, by rows (the lower triangle is written) */
} myb_oem_segment_data;

typedef struct myb_oem_data {
  const char* originator;
  myb_calendar creation_date;
  size_t segment_count;
  const myb_oem_segment_data* segments;
} myb_oem_data;

/* The message (version 2.0) as text. Numbers are the shortest decimals in the file's units
 * that the reader turns back into the same values: a message written and read gives the
 * states back exactly, but for a few values in a hundred that the conversion between km and m
 * cannot reach and that come back one unit of the last place away; and a message read and
 * written keeps the numbers as they were written (7000 km stays 7000). Seconds of epochs get
 * as many decimals as they need, three at least. *length (may be NULL): the characters of the
 * whole text without the NUL. Errors: MYB_E_OVERFLOW (the buffer; call again with length + 1
 * bytes), MYB_E_DOMAIN (no segments, a segment without states, a missing text, values not
 * finite, an invalid date). */
MYB_API MYB_NODISCARD myb_status myb_oem_write(const myb_oem_data* data, char* buffer,
                                               size_t size, size_t* length);
/* The same in the notation given (myb_oem_write: KVN). */
MYB_API MYB_NODISCARD myb_status myb_oem_write_as(const myb_oem_data* data,
                                                  myb_odm_notation notation, char* buffer,
                                                  size_t size, size_t* length);

/* ---- OPM: state at one epoch ------------------------------------------------------------------
 * The state (file: km, km/s; here m, m/s), optionally the osculating elements (km, degrees,
 * km^3/s^2; here m, rad, m^3/s^2), the parameters of the spacecraft (kg, m^2) and manoeuvres. */

typedef struct myb_opm_maneuver {
  myb_calendar epoch_ignition;
  double duration;   /* s */
  double delta_mass; /* kg, negative for a burn */
  char ref_frame[MYB_ODM_TEXT];
  double dv[3]; /* m/s in ref_frame */
} myb_opm_maneuver;

typedef struct myb_opm {
  char originator[MYB_ODM_TEXT];
  char object_name[MYB_ODM_TEXT];
  char object_id[MYB_ODM_TEXT];
  char center_name[MYB_ODM_TEXT];
  char ref_frame[MYB_ODM_TEXT];
  char time_system[MYB_ODM_TEXT];
  myb_calendar creation_date;
  myb_calendar epoch;
  double state[6];
  int has_elements;
  double semi_major_axis;
  double eccentricity;
  double inclination;
  double ra_of_asc_node;
  double arg_of_pericenter;
  double anomaly;
  int anomaly_is_mean; /* MEAN_ANOMALY instead of TRUE_ANOMALY */
  double gm;
  int has_spacecraft;
  double mass;
  double solar_rad_area;
  double solar_rad_coeff;
  double drag_area;
  double drag_coeff;
  size_t maneuver_count;
  myb_opm_maneuver maneuvers[MYB_OPM_MAX_MANEUVERS];
  int has_covariance;
  char cov_ref_frame[MYB_ODM_TEXT]; /* empty: that of the state */
  double covariance[36];            /* by rows: m^2, m^2/s, m^2/s^2 */
} myb_opm;

/* Errors: MYB_E_FORMAT (diag gives the line: not an OPM, a missing mandatory keyword, a value
 * that is not a number or a date, more than MYB_OPM_MAX_MANEUVERS manoeuvres); the structure
 * is then unspecified. */
MYB_API MYB_NODISCARD myb_status myb_opm_parse(const char* text, size_t size, myb_opm* opm,
                                               myb_diag* diag);
/* As myb_oem_write and myb_oem_write_as. */
MYB_API MYB_NODISCARD myb_status myb_opm_write(const myb_opm* opm, char* buffer, size_t size,
                                               size_t* length);
MYB_API MYB_NODISCARD myb_status myb_opm_write_as(const myb_opm* opm, myb_odm_notation notation,
                                                  char* buffer, size_t size, size_t* length);

/* ---- OMM: mean elements -----------------------------------------------------------------------
 * Mean elements of a theory named by MEAN_ELEMENT_THEORY (SGP4 for the elements of a TLE),
 * meaningful only with that theory: they are kept in the units of the message — rev/day,
 * degrees, km, km^3/s^2 — exactly as written, so that the theory receives the numbers of the
 * catalogue (ext/sgp4 turns a message into a myb_tle and back). */

typedef struct myb_omm {
  char originator[MYB_ODM_TEXT];
  char object_name[MYB_ODM_TEXT];
  char object_id[MYB_ODM_TEXT];
  char center_name[MYB_ODM_TEXT];
  char ref_frame[MYB_ODM_TEXT];
  char time_system[MYB_ODM_TEXT];
  char mean_element_theory[MYB_ODM_TEXT];
  myb_calendar creation_date;
  myb_calendar epoch;
  int has_mean_motion;    /* MEAN_MOTION given (otherwise SEMI_MAJOR_AXIS) */
  double mean_motion;     /* rev/day */
  double semi_major_axis; /* km */
  double eccentricity;
  double inclination;       /* deg */
  double ra_of_asc_node;    /* deg */
  double arg_of_pericenter; /* deg */
  double mean_anomaly;      /* deg */
  int has_gm;
  double gm; /* km^3/s^2 */
  int has_tle; /* the TLE-related parameters below */
  int ephemeris_type;
  char classification_type;
  long norad_cat_id;
  long element_set_no;
  long rev_at_epoch;
  double bstar;            /* 1 / Earth radii */
  double mean_motion_dot;  /* rev/day^2 */
  double mean_motion_ddot; /* rev/day^3 */
  int has_covariance;
  char cov_ref_frame[MYB_ODM_TEXT];
  double covariance[36]; /* by rows, in the units of the message: km^2, km^2/s, km^2/s^2 */
} myb_omm;

MYB_API MYB_NODISCARD myb_status myb_omm_parse(const char* text, size_t size, myb_omm* omm,
                                               myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_omm_write(const myb_omm* omm, char* buffer, size_t size,
                                               size_t* length);
MYB_API MYB_NODISCARD myb_status myb_omm_write_as(const myb_omm* omm, myb_odm_notation notation,
                                                  char* buffer, size_t size, size_t* length);

/* ---- OCM: the comprehensive message ------------------------------------------------------------
 * Version 3 of the standard: a header, metadata (TIME_SYSTEM and EPOCH_TZERO are mandatory)
 * and blocks of seven kinds, any number of each: trajectories, the physical properties of the
 * object, covariances, manoeuvres, the perturbations of the orbit determination, its
 * description and user-defined values. A block is a list of keywords and — for trajectories,
 * covariances and manoeuvres — lines of data: a time, absolute or in seconds from
 * EPOCH_TZERO, and the values that the keywords of the block describe (TRAJ_TYPE and
 * TRAJ_UNITS, COV_TYPE, MAN_COMPOSITION). The standard defines some two hundred keywords and
 * dozens of types of trajectories and covariances: the message is given here as it is —
 * keywords by name, lines as text and as numbers in the units of the message — and the
 * Cartesian trajectories also in SI. */

typedef struct myb_ocm myb_ocm;

typedef enum myb_ocm_kind {
  MYB_OCM_TRAJ = 0,
  MYB_OCM_PHYS = 1,
  MYB_OCM_COV = 2,
  MYB_OCM_MAN = 3,
  MYB_OCM_PERT = 4,
  MYB_OCM_OD = 5,
  MYB_OCM_USER = 6
} myb_ocm_kind;

/* Errors: as myb_oem_load. */
MYB_API MYB_NODISCARD myb_status myb_ocm_load(const char* text, size_t size,
                                              const myb_allocator* allocator, myb_ocm** ocm,
                                              myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_ocm_load_file(const char* path,
                                                   const myb_allocator* allocator, myb_ocm** ocm,
                                                   myb_diag* diag);
MYB_API void myb_ocm_destroy(myb_ocm* ocm); /* accepts NULL */

/* The value of a keyword of the header, of the metadata, or of block `index` of a kind, as
 * written; NULL if it is not there. A user-defined value is asked for as USER_DEFINED_<name>
 * in both notations. */
MYB_API const char* myb_ocm_header(const myb_ocm* ocm, const char* key);
MYB_API const char* myb_ocm_meta(const myb_ocm* ocm, const char* key);
MYB_API size_t myb_ocm_block_count(const myb_ocm* ocm, myb_ocm_kind kind);
MYB_API const char* myb_ocm_value(const myb_ocm* ocm, myb_ocm_kind kind, size_t index,
                                  const char* key);

/* The lines of data of a block: how many, and one as text (its words parted by blanks). */
MYB_API size_t myb_ocm_line_count(const myb_ocm* ocm, myb_ocm_kind kind, size_t index);
MYB_API const char* myb_ocm_line_text(const myb_ocm* ocm, myb_ocm_kind kind, size_t index,
                                      size_t line);

/* A line as numbers: its time — *absolute != 0: the calendar label *epoch; 0: *offset seconds
 * from EPOCH_TZERO — and the values after it in the units of the message, the first `max` of
 * them into values (may be NULL with max 0); *count: how many the line has. Errors:
 * MYB_E_FORMAT (a word that is neither a time nor a number: lines of manoeuvres may carry
 * words — take those from the text). */
MYB_API MYB_NODISCARD myb_status myb_ocm_line(const myb_ocm* ocm, myb_ocm_kind kind,
                                              size_t index, size_t line, int* absolute,
                                              myb_calendar* epoch, double* offset,
                                              double* values, size_t max, size_t* count);

/* A line of a trajectory of TRAJ_TYPE CARTP, CARTPV or CARTPVA in m, m/s and m/s^2: *count is
 * 3, 6 or 9, the rest of state is zero. Errors: MYB_E_DOMAIN (another type), MYB_E_FORMAT. */
MYB_API MYB_NODISCARD myb_status myb_ocm_state(const myb_ocm* ocm, size_t index, size_t line,
                                               int* absolute, myb_calendar* epoch,
                                               double* offset, double state[9], int* count);

/* What a message is written from. Keywords are upper-case letters, digits and '_'; a
 * user-defined one is USER_DEFINED_<name>. Lines are given as text (myb_ocm_format_line
 * makes one from numbers) and belong to trajectories, covariances and manoeuvres only. */
typedef struct myb_ocm_pair {
  const char* key;
  const char* value;
} myb_ocm_pair;

typedef struct myb_ocm_block_data {
  myb_ocm_kind kind;
  size_t pair_count;
  const myb_ocm_pair* pairs;
  size_t line_count;
  const char* const* lines;
} myb_ocm_block_data;

typedef struct myb_ocm_data {
  const char* originator;
  myb_calendar creation_date;
  size_t meta_count;
  const myb_ocm_pair* meta; /* with TIME_SYSTEM and EPOCH_TZERO */
  size_t block_count;
  const myb_ocm_block_data* blocks;
} myb_ocm_data;

/* The message (version 3.0) as text, as myb_oem_write_as. Errors: MYB_E_OVERFLOW,
 * MYB_E_DOMAIN (a keyword or a value that cannot be written, metadata without TIME_SYSTEM or
 * EPOCH_TZERO, lines in a block that has none). */
MYB_API MYB_NODISCARD myb_status myb_ocm_write(const myb_ocm_data* data,
                                               myb_odm_notation notation, char* buffer,
                                               size_t size, size_t* length);

/* A line of data: the time — the label *epoch, or `offset` seconds if epoch is NULL — and
 * `count` values, each the shortest decimal that reads back. Errors: MYB_E_OVERFLOW (call
 * again with *length + 1 bytes), MYB_E_DOMAIN. */
MYB_API MYB_NODISCARD myb_status myb_ocm_format_line(const myb_calendar* epoch, double offset,
                                                     const double* values, size_t count,
                                                     char* buffer, size_t size, size_t* length);

#ifdef __cplusplus
}
#endif

#endif /* MYB_CCSDS_H_ */
