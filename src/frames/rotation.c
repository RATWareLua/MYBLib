/* SPDX-License-Identifier: BSD-3-Clause */
#include <string.h>

#include "core/base.h"
#include "frames/rotation.h"
#include "myb/cmath.h"
#include "myb/frames.h"

void myb_internal_mat3_mul(const double a[9], const double b[9], double out[9]) {
  double r[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      r[3 * i + j] = a[3 * i] * b[j] + a[3 * i + 1] * b[3 + j] + a[3 * i + 2] * b[6 + j];
    }
  }
  memcpy(out, r, sizeof r);
}

void myb_internal_mat3_vec(const double m[9], const double v[3], double out[3]) {
  const double x = m[0] * v[0] + m[1] * v[1] + m[2] * v[2];
  const double y = m[3] * v[0] + m[4] * v[1] + m[5] * v[2];
  const double z = m[6] * v[0] + m[7] * v[1] + m[8] * v[2];
  out[0] = x;
  out[1] = y;
  out[2] = z;
}

void myb_internal_mat3_transpose(const double m[9], double out[9]) {
  double r[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      r[3 * j + i] = m[3 * i + j];
    }
  }
  memcpy(out, r, sizeof r);
}

void myb_internal_rot_axis(int axis, double angle, double m[9], double dm[9]) {
  const double c = myb_cos(angle);
  const double s = myb_sin(angle);
  memset(m, 0, 9 * sizeof *m);
  memset(dm, 0, 9 * sizeof *dm);
  /* Frame rotation: the coordinates of a fixed vector in axes turned by +angle about `axis`. */
  const int i = (axis + 1) % 3;
  const int j = (axis + 2) % 3;
  m[3 * axis + axis] = 1.0;
  m[3 * i + i] = c;
  m[3 * i + j] = s;
  m[3 * j + i] = -s;
  m[3 * j + j] = c;
  dm[3 * i + i] = -s;
  dm[3 * i + j] = c;
  dm[3 * j + i] = -c;
  dm[3 * j + j] = -s;
}

void myb_rotation_apply_position(const myb_rotation* rotation, const double in[3],
                                 double out[3]) {
  MYB_CHECK(rotation != NULL && in != NULL && out != NULL);
  myb_internal_mat3_vec(rotation->m, in, out);
}

void myb_rotation_apply_state(const myb_rotation* rotation, const double in[6], double out[6]) {
  MYB_CHECK(rotation != NULL && in != NULL && out != NULL);
  double r[3];
  double v[3];
  double w[3];
  myb_internal_mat3_vec(rotation->m, in, r);
  myb_internal_mat3_vec(rotation->m, in + 3, v);
  myb_internal_mat3_vec(rotation->dm, in, w);
  for (int i = 0; i < 3; ++i) {
    out[i] = r[i];
    out[i + 3] = v[i] + w[i];
  }
}

void myb_rotation_compose(const myb_rotation* second, const myb_rotation* first,
                          myb_rotation* out) {
  MYB_CHECK(second != NULL && first != NULL && out != NULL);
  double m[9];
  double a[9];
  double b[9];
  myb_internal_mat3_mul(second->m, first->m, m);
  myb_internal_mat3_mul(second->dm, first->m, a);
  myb_internal_mat3_mul(second->m, first->dm, b);
  for (int i = 0; i < 9; ++i) {
    out->m[i] = m[i];
    out->dm[i] = a[i] + b[i];
  }
}

void myb_rotation_invert(const myb_rotation* rotation, myb_rotation* out) {
  MYB_CHECK(rotation != NULL && out != NULL);
  myb_internal_mat3_transpose(rotation->m, out->m);
  myb_internal_mat3_transpose(rotation->dm, out->dm);
}

void myb_rotation_euler313(const double angles[3], const double rates[3],
                           myb_rotation* rotation) {
  MYB_CHECK(angles != NULL && rates != NULL && rotation != NULL);
  double a[9];
  double da[9];
  double b[9];
  double db[9];
  double c[9];
  double dc[9];
  myb_internal_rot_axis(2, angles[2], a, da);
  myb_internal_rot_axis(0, angles[1], b, db);
  myb_internal_rot_axis(2, angles[0], c, dc);
  double bc[9];
  double ab[9];
  double t1[9];
  double t2[9];
  double t3[9];
  myb_internal_mat3_mul(b, c, bc);
  myb_internal_mat3_mul(a, b, ab);
  myb_internal_mat3_mul(a, bc, rotation->m);
  myb_internal_mat3_mul(da, bc, t1);
  double dbc[9];
  myb_internal_mat3_mul(db, c, dbc);
  myb_internal_mat3_mul(a, dbc, t2);
  myb_internal_mat3_mul(ab, dc, t3);
  for (int i = 0; i < 9; ++i) {
    rotation->dm[i] = t1[i] * rates[2] + t2[i] * rates[1] + t3[i] * rates[0];
  }
}
