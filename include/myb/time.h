/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/time.h — instants, calendars and time scales.
 *
 * An instant is stored as whole seconds (int64) plus a fraction in [0, 1) (double), counted
 * from J2000.0 of its scale (2000-01-01T12:00:00, JD 2451545.0). Resolution is better than
 * 10^-16 s over ±2.9e11 years; a single double Julian date would resolve only ~40 µs today.
 *
 * Time scales are not built into the library: a myb_timesys is created from kernel
 * variables (data/time/time-scales.tk defines TAI, TT, TCG, TDB, TCB, GPS, GST, BDT with the
 * IAU/IERS constants) and is a tree of relations between scales. UTC is not a scale of the
 * tree: it labels an atomic scale with leap seconds and exists only as calendar labels
 * (data/time/utc-tai.tk). Calendars are proleptic Gregorian with astronomical year numbering
 * (year 0 = 1 BC). */
#ifndef MYB_TIME_H_
#define MYB_TIME_H_

#include <stdint.h>

#include "myb/core.h"
#include "myb/io.h"
#include "myb/num.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Instants ------------------------------------------------------------------------------ */

typedef struct myb_time {
  int64_t s; /* whole seconds since J2000.0 of the scale */
  double f;  /* fraction of a second, 0 <= f < 1 */
} myb_time;

/* Conversions and arithmetic. Results are normalised (0 <= f < 1); the only rounding is that
 * of the final fraction (at most 2^-53 s). Arguments must be finite with |seconds| < 2^62. */
MYB_API myb_time myb_time_from_seconds(double seconds);
MYB_API myb_time myb_time_from_dd(myb_dd seconds);
MYB_API myb_time myb_time_add(myb_time t, double seconds);
MYB_API myb_time myb_time_add_dd(myb_time t, myb_dd seconds);
MYB_API myb_dd myb_time_diff(myb_time a, myb_time b); /* a - b, relative error < 2^-103 */
MYB_API double myb_time_seconds(myb_time t);          /* RN(seconds since J2000.0) */
MYB_API int myb_time_compare(myb_time a, myb_time b); /* -1, 0, 1 */

/* Two-part Julian date jd1 + jd2 (days). From: any split. To: jd1 is 2451545.0 plus whole
 * days and 0 <= jd2 < 1, which keeps the full resolution in jd2. */
MYB_API myb_time myb_time_from_jd(double jd1, double jd2);
MYB_API void myb_time_to_jd(myb_time t, double* jd1, double* jd2);

/* ---- Calendar ------------------------------------------------------------------------------ */

typedef struct myb_calendar {
  int64_t year;  /* astronomical numbering */
  int month;     /* 1..12 */
  int day;       /* 1..31 */
  int hour;      /* 0..23 */
  int minute;    /* 0..59 */
  double second; /* [0, 60) on a continuous scale; UTC allows [0, 60 + leap) */
} myb_calendar;

/* Calendar date and time of a continuous scale -> instant of that scale. MYB_E_DOMAIN for an
 * invalid date or time (including second >= 60) or |year| > 10^9. */
MYB_API MYB_NODISCARD myb_status myb_calendar_to_time(const myb_calendar* calendar, myb_time* t);
MYB_API void myb_time_to_calendar(myb_time t, myb_calendar* calendar);

/* ---- Time scale systems --------------------------------------------------------------------- */

typedef struct myb_timesys myb_timesys;

typedef struct myb_scale {
  int index;
} myb_scale;

/* Builds the scale tree from kernel variables (format in data/time/time-scales.tk; the UTC
 * table, data/time/utc-tai.tk, is optional). The system copies what it needs: the pool may
 * be destroyed afterwards. Errors: MYB_E_NOT_FOUND (missing variable), MYB_E_FORMAT
 * (inconsistent definitions, e.g. a cycle or an unordered UTC table), MYB_E_NOMEM; diag
 * names the offending variable. */
MYB_API MYB_NODISCARD myb_status myb_timesys_create(const myb_kpool* pool,
                                                    const myb_allocator* allocator,
                                                    myb_timesys** system, myb_diag* diag);
MYB_API void myb_timesys_destroy(myb_timesys* system); /* NULL is allowed */

/* Scale by name ("TT"). MYB_E_NOT_FOUND if the system has no such scale. */
MYB_API MYB_NODISCARD myb_status myb_timesys_scale(const myb_timesys* system, const char* name,
                                                   myb_scale* scale);
MYB_API const char* myb_timesys_scale_name(const myb_timesys* system, myb_scale scale);

/* The same instant on another scale. Exact relations (offsets) add no error beyond the final
 * rounding of the fraction; rate relations are evaluated in double-word arithmetic; series
 * relations are iterated to below 10^-15 s. Errors: the statuses of a function set by
 * myb_timesys_set_difference. */
MYB_API MYB_NODISCARD myb_status myb_timesys_convert(const myb_timesys* system, myb_time t,
                                                     myb_scale from, myb_scale to, myb_time* out);

/* A relation of the kind SERIES given by a function in place of its series: the difference
 * scale - parent (s) at an instant of the scale — TDB - TT from the time ephemeris of a
 * planetary ephemeris (myb_spk_time_difference of myb/kin.h) in place of the analytic
 * series. function = NULL returns to the series. A table made before is dropped. The
 * context must outlive the system. Not to be called while another thread converts. Errors:
 * MYB_E_DOMAIN (a scale whose relation is not a series). */
typedef myb_status (*myb_scale_difference)(void* context, myb_time t, double* difference);

MYB_API MYB_NODISCARD myb_status myb_timesys_set_difference(myb_timesys* system, myb_scale scale,
                                                            myb_scale_difference function,
                                                            void* context);

/* The relation of the kind SERIES of a scale tabulated every `step` seconds (0: an hour)
 * between the instants t0 and t1 of the scale: conversions through it then take the cubic
 * through the four nodes around the instant in place of the series (787 sines for TDB - TT)
 * or of the function, and outside the table what defines the relation as before. With a
 * step of an hour the table of TDB - TT is within 4e-16 s of the series; its values depend
 * on the nodes only, not on what was asked before. A second call replaces the table.
 * Errors: MYB_E_DOMAIN (another kind of relation, t1 not after t0, a step that is negative
 * or not finite, more than 2^26 nodes), MYB_E_NOMEM, the statuses of the function. */
MYB_API MYB_NODISCARD myb_status myb_timesys_tabulate(myb_timesys* system, myb_scale scale,
                                                      myb_time t0, myb_time t1, double step);

/* UTC calendar label -> instant on scale `to`. MYB_E_DOMAIN for an invalid label (e.g.
 * 23:59:60 on a day without a leap second), MYB_E_RANGE before the first table row or from
 * the table's expiry date on, MYB_E_NOT_FOUND if the system has no UTC table. */
MYB_API MYB_NODISCARD myb_status myb_timesys_from_utc(const myb_timesys* system,
                                                      const myb_calendar* utc, myb_scale to,
                                                      myb_time* out);

/* Instant on scale `from` -> UTC calendar label (23:59:60.x during a positive leap second). */
MYB_API MYB_NODISCARD myb_status myb_timesys_to_utc(const myb_timesys* system, myb_time t,
                                                    myb_scale from, myb_calendar* utc);

#ifdef __cplusplus
}
#endif

#endif /* MYB_TIME_H_ */
