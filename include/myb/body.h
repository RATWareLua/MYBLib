/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/body.h — parameters with uncertainty and the model of a body.
 *
 * The library knows no particular body: a body is a set of parameters supplied by the user
 * (usually from kernels: GM from an ephemeris such as gm_de440.tpc, fields from gravity
 * models, ...). Every parameter carries its standard uncertainty, so that the uncertainty of
 * a result can be computed rather than assumed.
 *
 * Present components: GM and the kind of motion. Shape, gravity field, rotation, tides,
 * atmosphere and radiation are added with their modules (docs/PLAN.en.md, §3). */
#ifndef MYB_BODY_H_
#define MYB_BODY_H_

#include <stdint.h>

#include "myb/core.h"
#include "myb/io.h"

#ifdef __cplusplus
extern "C" {
#endif

/* A value and its standard uncertainty (1 sigma). sigma = 0: exact (defining) value. */
typedef struct myb_param {
  double value;
  double sigma;
} myb_param;

/* Product a * b with first-order uncertainty for independent errors:
 * sigma^2 = (b sigma_a)^2 + (a sigma_b)^2, evaluated without overflow.
 * MYB_E_DOMAIN for non-finite values or negative sigmas. */
MYB_API MYB_NODISCARD myb_status myb_param_product(myb_param a, myb_param b, myb_param* product);

/* GM = G * M for a body known by its mass (kg) and a value of G (m^3 kg^-1 s^-2, e.g.
 * CODATA 2022 in data/constants/codata-2022.tk). The relative uncertainty of G (2.2e-5) then
 * dominates; prefer a measured GM. MYB_E_DOMAIN for negative or non-finite inputs. */
MYB_API MYB_NODISCARD myb_status myb_gm_from_mass(myb_param mass, myb_param g, myb_param* gm);

/* Reads a parameter from kernel variables: the value from `name` (exactly one number), the
 * uncertainty from `sigma_name` (exactly one number; NULL means exact), both multiplied by
 * `scale` (unit conversion, e.g. 1e9 for km^3 s^-2 -> m^3 s^-2; the multiplication is the
 * only rounding). Errors: MYB_E_NOT_FOUND, MYB_E_FORMAT (not one number), MYB_E_DOMAIN
 * (non-finite result or negative sigma). */
MYB_API MYB_NODISCARD myb_status myb_param_from_kpool(const myb_kpool* pool, const char* name,
                                                      const char* sigma_name, double scale,
                                                      myb_param* param);

typedef enum myb_motion {
  MYB_MOTION_DYNAMIC = 1,  /* the state is integrated */
  MYB_MOTION_KINEMATIC = 2 /* the state is given by a provider (ephemeris, table, function) */
} myb_motion;

typedef struct myb_body {
  int64_t id;   /* opaque identifier chosen by the user (e.g. a NAIF ID) */
  myb_param gm; /* m^3 s^-2; 0 for a test particle */
  myb_motion motion;
} myb_body;

/* MYB_OK if the body is physically meaningful: finite GM >= 0, finite sigma >= 0, a known
 * kind of motion; MYB_E_DOMAIN otherwise. */
MYB_API MYB_NODISCARD myb_status myb_body_validate(const myb_body* body);

#ifdef __cplusplus
}
#endif

#endif /* MYB_BODY_H_ */
