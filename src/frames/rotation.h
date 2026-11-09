/* SPDX-License-Identifier: BSD-3-Clause */
/* frames/rotation.h — 3x3 matrix helpers (rows, fixed order of operations) and elementary
 * frame rotations. */
#ifndef MYB_SRC_FRAMES_ROTATION_H_
#define MYB_SRC_FRAMES_ROTATION_H_

/* out = a b; out may alias a or b. */
void myb_internal_mat3_mul(const double a[9], const double b[9], double out[9]);

/* out = m v; out may alias v. */
void myb_internal_mat3_vec(const double m[9], const double v[3], double out[3]);

/* out = m^T; out may alias m. */
void myb_internal_mat3_transpose(const double m[9], double out[9]);

/* Rotation of the frame by `angle` about axis 0 (x), 1 (y) or 2 (z), and its derivative with
 * respect to the angle. */
void myb_internal_rot_axis(int axis, double angle, double m[9], double dm[9]);

#endif /* MYB_SRC_FRAMES_ROTATION_H_ */
