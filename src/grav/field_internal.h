/* SPDX-License-Identifier: BSD-3-Clause */
/* grav/field_internal.h — direct access to the effective coefficients of a field for the
 * models that change them in bulk (tides): the same effect as myb_gravity_field_add_correction
 * without a call and a check per coefficient. */
#ifndef MYB_SRC_GRAV_FIELD_INTERNAL_H_
#define MYB_SRC_GRAV_FIELD_INTERNAL_H_

#include "myb/grav.h"

/* Effective C and S, indexed n (n + 1) / 2 + m. Writers must call
 * myb_internal_field_touched with the highest degree they changed, so that the next
 * myb_gravity_field_set_time restores the coefficients. */
double* myb_internal_field_c(myb_gravity_field* field);
double* myb_internal_field_s(myb_gravity_field* field);
void myb_internal_field_touched(myb_gravity_field* field, int degree);

#endif /* MYB_SRC_GRAV_FIELD_INTERNAL_H_ */
