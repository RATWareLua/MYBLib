/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/io.h — reading external data: numbers and NAIF text kernels.
 *
 * Numbers: myb_parse_double converts text to the correctly rounded binary64 value. It does not
 * depend on the C locale (the decimal separator is always '.') nor on the C runtime's strtod.
 *
 * Data files: myb reads every convention and catalogue (time scales, leap seconds, series
 * coefficients, body constants) from NAIF text kernels ("Kernel Required Reading", NAIF,
 * JPL): text outside \begindata ... \begintext blocks is commentary; inside, assignments
 *     NAME  = value                 NAME  = ( value value, value ... )
 *     NAME += ( values to append )
 * with numbers (Fortran 'D' exponents allowed), strings in single quotes ('' is a quote) and
 * @-dates (kept as strings). The same reader loads JPL's own kernels (gm_de440.tpc,
 * pck00011.tpc, ...), so the library never needs a copy of their data. */
#ifndef MYB_IO_H_
#define MYB_IO_H_

#include <stddef.h>

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Numbers ---------------------------------------------------------------------------- */

/* Parses one number at the start of text[0 .. length):
 *   [+-] digits [. digits] | [+-] . digits   followed by optional [eEdD] [+-] digits
 *   [+-] 0x hexdigits [. hexdigits] [pP [+-] digits]
 *   [+-] inf | infinity | nan                (case-insensitive)
 * No leading white space is skipped. On success *value is the correctly rounded binary64
 * value (round to nearest, ties to even; overflow gives infinity, underflow is gradual) and
 * *consumed (may be NULL) the number of characters used. MYB_E_FORMAT if no number starts
 * there; outputs are then untouched. */
MYB_API MYB_NODISCARD myb_status myb_parse_double(const char* text, size_t length, double* value,
                                                  size_t* consumed);

/* ---- Kernel pool -------------------------------------------------------------------------- */

typedef struct myb_kpool myb_kpool;

typedef enum myb_kvar_type {
  MYB_KVAR_NUMBERS = 1,
  MYB_KVAR_STRINGS = 2 /* strings and @-dates */
} myb_kvar_type;

/* Creates an empty pool. allocator may be NULL (malloc/free). */
MYB_API MYB_NODISCARD myb_status myb_kpool_create(const myb_allocator* allocator, myb_kpool** pool);
MYB_API void myb_kpool_destroy(myb_kpool* pool); /* NULL is allowed */

/* Loads kernel text. Loading is atomic: on any error the pool is unchanged and diag (may be
 * NULL) tells the line and the reason. Later assignments replace (=) or extend (+=) earlier
 * ones, also across files. Errors: MYB_E_FORMAT, MYB_E_NOMEM, MYB_E_OVERFLOW. */
MYB_API MYB_NODISCARD myb_status myb_kpool_load_text(myb_kpool* pool, const char* text, size_t size,
                                                     myb_diag* diag);

/* Loads a kernel file; path is UTF-8 on every platform. Additional error: MYB_E_IO. */
MYB_API MYB_NODISCARD myb_status myb_kpool_load_file(myb_kpool* pool, const char* path,
                                                     myb_diag* diag);

/* Type and number of values of a variable. MYB_E_NOT_FOUND if it does not exist. */
MYB_API MYB_NODISCARD myb_status myb_kpool_describe(const myb_kpool* pool, const char* name,
                                                    myb_kvar_type* type, size_t* count);

/* Numeric values of a variable. MYB_E_NOT_FOUND, or MYB_E_FORMAT if it holds strings.
 * The pointer stays valid until the next load into the pool or its destruction. */
MYB_API MYB_NODISCARD myb_status myb_kpool_numbers(const myb_kpool* pool, const char* name,
                                                   const double** values, size_t* count);

/* One string value (not NUL-terminated: use *length). MYB_E_NOT_FOUND, MYB_E_FORMAT if the
 * variable is numeric, MYB_E_RANGE if index >= count. Pointer lifetime as above. */
MYB_API MYB_NODISCARD myb_status myb_kpool_string(const myb_kpool* pool, const char* name,
                                                  size_t index, const char** value, size_t* length);

#ifdef __cplusplus
}
#endif

#endif /* MYB_IO_H_ */
