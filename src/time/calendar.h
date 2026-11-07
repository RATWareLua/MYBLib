/* SPDX-License-Identifier: BSD-3-Clause */
/* time/calendar.h — integer calendar helpers shared by the time module. */
#ifndef MYB_SRC_TIME_CALENDAR_H_
#define MYB_SRC_TIME_CALENDAR_H_

#include <stdint.h>

/* floor(a / b) for b > 0. */
int64_t myb_internal_floor_div(int64_t a, int64_t b);

/* Proleptic Gregorian date <-> days since 1970-01-01 (Hinnant's algorithms). */
int64_t myb_internal_days_from_civil(int64_t year, int month, int day);
void myb_internal_civil_from_days(int64_t days, int64_t* year, int* month, int* day);

int myb_internal_days_in_month(int64_t year, int month);

/* 1 if year/month/day is a valid date with |year| <= 10^9. */
int myb_internal_valid_date(int64_t year, int month, int day);

/* Modified Julian Date 0 is 1858-11-17, i.e. day -40587 from 1970-01-01. */
enum { MYB_MJD_OF_UNIX_EPOCH = 40587 };

#endif /* MYB_SRC_TIME_CALENDAR_H_ */
