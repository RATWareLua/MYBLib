/* SPDX-License-Identifier: BSD-3-Clause */
/* kin/spk_segment.h — one SPK segment: validation of its layout when a file is loaded, and the
 * evaluation of its state, in the file's units (km, km/s).
 *
 * Layouts and evaluation rules follow "SPK Required Reading" (NAIF, JPL); the choice of records
 * and interpolation windows, and the conventions of types 14, 17, 19 and 20, were checked
 * against the evaluations of CSPICE N0067 (tests/data/spk). */
#ifndef MYB_SRC_KIN_SPK_SEGMENT_H_
#define MYB_SRC_KIN_SPK_SEGMENT_H_

#include <stddef.h>
#include <stdint.h>

#include "io/daf.h"
#include "myb/core.h"
#include "myb/time.h"

enum {
  MYB_SPK_MAX_WINDOW = 32,      /* states in an interpolation window (SPICE: at most 28) */
  MYB_SPK_MAX_COEFFICIENTS = 64, /* Chebyshev coefficients per component (SPICE: <= 51) */
  MYB_SPK_MAX_MDA = 25           /* MAXTRM of the difference lines of types 1 and 21 */
};

typedef struct myb_internal_spk_segment {
  myb_internal_daf daf; /* the file the segment lives in */
  int32_t target;
  int32_t center;
  int32_t frame;
  int32_t type;
  double start; /* coverage, seconds past J2000 */
  double stop;
  uint64_t begin; /* first and last word address of the segment */
  uint64_t end;
  size_t file;
  char name[41];
  int supported;

  /* Layout, validated at load time; meaning by type. */
  uint64_t count;    /* records, states, packets or intervals */
  uint64_t size;     /* words per record or packet */
  uint64_t epochs;   /* address of the first epoch (0: equally spaced) */
  uint64_t data;     /* type 14: address of the first packet */
  int window;        /* interpolation window, MDA dimension */
  int subtype;       /* type 19: 1 if the later interval is used on a boundary */
  double init;       /* first epoch or record start (s past J2000) */
  double step;       /* epoch spacing or record length (s) */
  double gm;         /* type 5 (km^3 s^-2) */
  double dscale;     /* type 20 (km) */
  double tscale;     /* type 20 (s) */
  myb_time init_time; /* type 20: start of the first record */
  int64_t offset;    /* type 14: words before the data of each packet */
} myb_internal_spk_segment;

/* Fills the layout of a segment whose summary fields (daf, target ... end) are set; sets
 * `supported`. MYB_E_FORMAT for an inconsistent layout of a supported type. */
[[nodiscard]] myb_status myb_internal_spk_segment_parse(myb_internal_spk_segment* segment,
                                                        myb_diag* diag);

/* State (km, km/s) of the segment at t (which the caller has checked to be covered).
 * Errors: MYB_E_FORMAT (unsupported type, invalid record data), MYB_E_DOMAIN (type 5 two-body
 * propagation failed). */
[[nodiscard]] myb_status myb_internal_spk_segment_state(const myb_internal_spk_segment* segment,
                                                        myb_time t, double state[6]);

/* Hermite interpolation at 0 of values y_i and derivatives yd_i at the m nodes d_i (m at
 * most MYB_SPK_MAX_WINDOW; Newton divided differences on the doubled nodes): the value and the
 * derivative. Shared with the CK reader. */
void myb_internal_hermite(int m, const double* d, const double* y, const double* yd,
                          double* value, double* derivative);

#endif /* MYB_SRC_KIN_SPK_SEGMENT_H_ */
