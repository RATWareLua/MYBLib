/* SPDX-License-Identifier: BSD-3-Clause */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "frames/rotation.h"
#include "myb/frames.h"

static const double pi = 0x1.921fb54442d18p+1;

static int string_is(const char* value, size_t length, const char* expected) {
  return length == strlen(expected) && memcmp(value, expected, length) == 0;
}

static myb_status get_string(const myb_kpool* pool, const char* name, const char** value,
                             size_t* length, myb_diag* diag) {
  const myb_status status = myb_kpool_string(pool, name, 0, value, length);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "%s: %s", name,
                      status == MYB_E_NOT_FOUND ? "not defined" : "not a string");
  }
  return status;
}

static myb_status get_numbers(const myb_kpool* pool, const char* name, size_t expected,
                              double* out, myb_diag* diag) {
  const double* v = NULL;
  size_t n = 0;
  const myb_status status = myb_kpool_numbers(pool, name, &v, &n);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "%s: %s", name,
                      status == MYB_E_NOT_FOUND ? "not defined" : "not numeric");
    return status;
  }
  if (n != expected) {
    myb_internal_diag(diag, 0, "%s: expected %zu numbers", name, expected);
    return MYB_E_FORMAT;
  }
  for (size_t i = 0; i < n; ++i) {
    if (!isfinite(v[i])) {
      myb_internal_diag(diag, 0, "%s: non-finite value", name);
      return MYB_E_FORMAT;
    }
    out[i] = v[i];
  }
  return MYB_OK;
}

/* Radians per unit, 0 for an unknown unit. */
static double unit_factor(const char* u, size_t n) {
  if (string_is(u, n, "RADIANS")) {
    return 1.0;
  }
  if (string_is(u, n, "DEGREES")) {
    return pi / 180.0;
  }
  if (string_is(u, n, "ARCMINUTES")) {
    return pi / 10800.0;
  }
  if (string_is(u, n, "ARCSECONDS")) {
    return pi / 648000.0;
  }
  if (string_is(u, n, "HOURANGLE")) {
    return pi / 12.0;
  }
  if (string_is(u, n, "MINUTEANGLE")) {
    return pi / 720.0;
  }
  if (string_is(u, n, "SECONDANGLE")) {
    return pi / 43200.0;
  }
  return 0.0;
}

myb_status myb_tkframe(const myb_kpool* pool, int64_t id, myb_rotation* rotation,
                       char relative[33], myb_diag* diag) {
  MYB_CHECK(pool != NULL && rotation != NULL && relative != NULL);
  myb_internal_diag_clear(diag);
  char name[64];
  const char* value = NULL;
  size_t length = 0;
  snprintf(name, sizeof name, "TKFRAME_%" PRId64 "_RELATIVE", id);
  myb_status status = get_string(pool, name, &value, &length, diag);
  if (status != MYB_OK) {
    return status;
  }
  if (length == 0 || length > 32) {
    myb_internal_diag(diag, 0, "%s: frame name of 1 .. 32 characters expected", name);
    return MYB_E_FORMAT;
  }
  char rel[33];
  memcpy(rel, value, length);
  rel[length] = '\0';
  snprintf(name, sizeof name, "TKFRAME_%" PRId64 "_SPEC", id);
  status = get_string(pool, name, &value, &length, diag);
  if (status != MYB_OK) {
    return status;
  }
  double to_relative[9]; /* TK -> relative */
  if (string_is(value, length, "ANGLES")) {
    double angles[3];
    double axes[3];
    snprintf(name, sizeof name, "TKFRAME_%" PRId64 "_ANGLES", id);
    status = get_numbers(pool, name, 3, angles, diag);
    if (status == MYB_OK) {
      snprintf(name, sizeof name, "TKFRAME_%" PRId64 "_AXES", id);
      status = get_numbers(pool, name, 3, axes, diag);
    }
    if (status != MYB_OK) {
      return status;
    }
    snprintf(name, sizeof name, "TKFRAME_%" PRId64 "_UNITS", id);
    status = get_string(pool, name, &value, &length, diag);
    if (status != MYB_OK) {
      return status;
    }
    const double factor = unit_factor(value, length);
    if (factor == 0.0) {
      myb_internal_diag(diag, 0, "%s: unknown unit", name);
      return MYB_E_FORMAT;
    }
    /* M = [a1]_i1 [a2]_i2 [a3]_i3 */
    double m[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    for (int k = 0; k < 3; ++k) {
      if (!(axes[k] == 1.0 || axes[k] == 2.0 || axes[k] == 3.0)) {
        myb_internal_diag(diag, 0, "TKFRAME_%" PRId64 "_AXES: axes must be 1, 2 or 3", id);
        return MYB_E_FORMAT;
      }
      double r[9];
      double unused[9];
      myb_internal_rot_axis((int)axes[k] - 1, angles[k] * factor, r, unused);
      myb_internal_mat3_mul(m, r, m);
    }
    memcpy(to_relative, m, sizeof m);
  } else if (string_is(value, length, "MATRIX")) {
    double columns[9];
    snprintf(name, sizeof name, "TKFRAME_%" PRId64 "_MATRIX", id);
    status = get_numbers(pool, name, 9, columns, diag);
    if (status != MYB_OK) {
      return status;
    }
    myb_internal_mat3_transpose(columns, to_relative); /* columns -> rows */
  } else if (string_is(value, length, "QUATERNION")) {
    double q[4];
    snprintf(name, sizeof name, "TKFRAME_%" PRId64 "_Q", id);
    status = get_numbers(pool, name, 4, q, diag);
    if (status != MYB_OK) {
      return status;
    }
    const double n = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    if (!(n > 0.0)) {
      myb_internal_diag(diag, 0, "%s: zero quaternion", name);
      return MYB_E_FORMAT;
    }
    const double c = q[0] / n;
    const double x = q[1] / n;
    const double y = q[2] / n;
    const double z = q[3] / n;
    /* SPICE q2m: the matrix rotating vectors by the quaternion's angle about its axis */
    const double m[9] = {1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - c * z), 2.0 * (x * z + c * y),
                         2.0 * (x * y + c * z), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - c * x),
                         2.0 * (x * z - c * y), 2.0 * (y * z + c * x), 1.0 - 2.0 * (x * x + y * y)};
    memcpy(to_relative, m, sizeof m);
  } else {
    myb_internal_diag(diag, 0, "%s: SPEC must be 'ANGLES', 'MATRIX' or 'QUATERNION'", name);
    return MYB_E_FORMAT;
  }
  /* a rotation: M M^T = I and det M = +1 */
  double mt[9];
  double product[9];
  myb_internal_mat3_transpose(to_relative, mt);
  myb_internal_mat3_mul(to_relative, mt, product);
  double worst = 0.0;
  for (int i = 0; i < 9; ++i) {
    worst = fmax(worst, fabs(product[i] - ((i % 4) == 0 ? 1.0 : 0.0)));
  }
  const double* a = to_relative;
  const double det = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6]) +
                     a[2] * (a[3] * a[7] - a[4] * a[6]);
  if (!(worst < 1e-10) || !(det > 0.0)) {
    myb_internal_diag(diag, 0, "TKFRAME_%" PRId64 ": not a rotation matrix", id);
    return MYB_E_FORMAT;
  }
  memcpy(rotation->m, mt, sizeof mt); /* relative -> TK */
  memset(rotation->dm, 0, sizeof rotation->dm);
  memcpy(relative, rel, sizeof rel);
  return MYB_OK;
}
