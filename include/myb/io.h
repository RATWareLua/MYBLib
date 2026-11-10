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

/* ---- Printing numbers ----------------------------------------------------------------------
 * The inverse of myb_parse_double, as exact: every digit follows from integer arithmetic on
 * the value, the same on every platform and independent of the C runtime and the locale.
 * Non-finite values print as "nan", "inf", "-inf". */

enum { MYB_DOUBLE_CHARS = 32 };

/* The shortest decimal that myb_parse_double reads back to the same value, the closest to it
 * among those of that length, in the form of Python's repr: positional for 1e-4 <= |x| < 1e16
 * ("0.001", "1234.5", "10.0"), otherwise with an exponent ("1.5e-07", "1e+300"); "0.0" and
 * "-0.0". Returns the number of characters written, without the terminating NUL (at most 24). */
MYB_API size_t myb_format_double(double value, char out[MYB_DOUBLE_CHARS]);

/* `decimals` (0 .. 60) digits after the point, as "%.Nf": the exact value rounded to nearest,
 * ties to even. *length (may be NULL) receives the number of characters of the full text
 * without the NUL; MYB_E_OVERFLOW if size cannot hold it (the text is then cut, still
 * terminated). Errors: MYB_E_DOMAIN (decimals outside the range). */
MYB_API MYB_NODISCARD myb_status myb_format_fixed(double value, int decimals, char* out,
                                                  size_t size, size_t* length);

/* `digits` (1 .. 60) significant digits with an exponent, as "%.(digits-1)e": "1.234e+05".
 * Otherwise as myb_format_fixed. */
MYB_API MYB_NODISCARD myb_status myb_format_scientific(double value, int digits, char* out,
                                                       size_t size, size_t* length);

/* ---- Tables of numbers ---------------------------------------------------------------------
 * Text of rows of numbers separated by commas, semicolons, tabs or blanks (CSV and the like):
 * every row the same number of fields; '#' starts a comment; a first row that is not numbers
 * names the columns; an empty field between separators is a NaN. */

typedef struct myb_table myb_table;

/* Errors: MYB_E_FORMAT (diag gives the line: a field that is not a number, a row of another
 * length, no values), MYB_E_NOMEM, MYB_E_OVERFLOW; MYB_E_IO (files). */
MYB_API MYB_NODISCARD myb_status myb_table_load(const char* text, size_t size,
                                                const myb_allocator* allocator,
                                                myb_table** table, myb_diag* diag);
MYB_API MYB_NODISCARD myb_status myb_table_load_file(const char* path,
                                                     const myb_allocator* allocator,
                                                     myb_table** table, myb_diag* diag);
MYB_API void myb_table_destroy(myb_table* table); /* accepts NULL */

MYB_API size_t myb_table_rows(const myb_table* table);
MYB_API size_t myb_table_columns(const myb_table* table);
/* The name of a column (at most 31 characters are kept); "" without a row of names. */
MYB_API const char* myb_table_name(const myb_table* table, size_t column);
/* The values of a row: myb_table_columns numbers. */
MYB_API const double* myb_table_row(const myb_table* table, size_t row);

/* rows x columns values (row by row) as comma-separated text, each the shortest decimal that
 * reads back to the same number, a NaN as an empty field; names (may be NULL): a first row of
 * column names. *length (may be NULL): the characters of the whole text. Errors:
 * MYB_E_OVERFLOW (the buffer; call again with length + 1 bytes), MYB_E_DOMAIN (no columns,
 * more than 1024, a missing name). */
MYB_API MYB_NODISCARD myb_status myb_table_write(const char* const* names, size_t columns,
                                                 const double* values, size_t rows,
                                                 char* buffer, size_t size, size_t* length);

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
