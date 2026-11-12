/* SPDX-License-Identifier: BSD-3-Clause */
/* What the levels that carry contacts (dyn/system.c) check of them when they are added. */
#ifndef MYB_SRC_FORCE_CONTACT_INTERNAL_H_
#define MYB_SRC_FORCE_CONTACT_INTERNAL_H_

#include <stddef.h>

#include "myb/contact.h"

/* 1 if the surface, the law and the points (count >= 1) are valid for myb_contact_force. */
int myb_internal_contact_valid(const myb_contact_surface* surface, const myb_contact_law* law,
                               const myb_contact_point* points, size_t count);

#endif /* MYB_SRC_FORCE_CONTACT_INTERNAL_H_ */
