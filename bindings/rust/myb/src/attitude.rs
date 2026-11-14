// SPDX-License-Identifier: BSD-3-Clause
//! Orientation and rotational dynamics of rigid bodies (`myb/attitude.h`): quaternions
//! (scalar first, Hamilton product, the frame rotation M from inertial to body axes), Euler
//! angles, principal axes, the equations of motion with and without a fluid core,
//! gravitational torques, the exact torque-free rotation and the splitting integrator.
use std::ffi::{c_int, c_void};

use crate::{check, domain, sys, Result, Vec3};

/// A quaternion (w, x, y, z).
pub type Quat = [f64; 4];

/// a * b.
pub fn quat_multiply(a: &Quat, b: &Quat) -> Quat {
    let mut out = [0.0; 4];
    // SAFETY: arrays of four numbers
    unsafe { sys::myb_quat_multiply(a.as_ptr(), b.as_ptr(), out.as_mut_ptr()) };
    out
}

pub fn quat_conjugate(q: &Quat) -> Quat {
    let mut out = [0.0; 4];
    // SAFETY: arrays of four numbers
    unsafe { sys::myb_quat_conjugate(q.as_ptr(), out.as_mut_ptr()) };
    out
}

/// q scaled to unit norm, and the former norm (q unchanged if it is 0).
pub fn quat_normalize(q: &Quat) -> (Quat, f64) {
    let mut out = *q;
    // SAFETY: an array of four numbers
    let norm = unsafe { sys::myb_quat_normalize(out.as_mut_ptr()) };
    (out, norm)
}

/// M(q), rows.
pub fn quat_to_matrix(q: &Quat) -> [f64; 9] {
    let mut m = [0.0; 9];
    // SAFETY: arrays of four and nine numbers
    unsafe { sys::myb_quat_to_matrix(q.as_ptr(), m.as_mut_ptr()) };
    m
}

/// The unit quaternion of a rotation matrix, w >= 0.
pub fn quat_from_matrix(m: &[f64; 9]) -> Quat {
    let mut q = [0.0; 4];
    // SAFETY: arrays of nine and four numbers
    unsafe { sys::myb_quat_from_matrix(m.as_ptr(), q.as_mut_ptr()) };
    q
}

/// M(q) v.
pub fn quat_rotate(q: &Quat, v: &Vec3) -> Vec3 {
    let mut out = [0.0; 3];
    // SAFETY: arrays of four and three numbers
    unsafe { sys::myb_quat_rotate(q.as_ptr(), v.as_ptr(), out.as_mut_ptr()) };
    out
}

/// The quaternion of a rotation vector angle * axis.
pub fn quat_from_rotation_vector(phi: &Vec3) -> Quat {
    let mut q = [0.0; 4];
    // SAFETY: arrays of three and four numbers
    unsafe { sys::myb_quat_from_rotation_vector(phi.as_ptr(), q.as_mut_ptr()) };
    q
}

/// The rotation vector of a quaternion (the angle in [0, pi]).
pub fn quat_to_rotation_vector(q: &Quat) -> Vec3 {
    let mut phi = [0.0; 3];
    // SAFETY: arrays of four and three numbers
    unsafe { sys::myb_quat_to_rotation_vector(q.as_ptr(), phi.as_mut_ptr()) };
    phi
}

/// Spherical linear interpolation along the shorter arc, t in [0, 1].
pub fn quat_slerp(a: &Quat, b: &Quat, t: f64) -> Quat {
    let mut out = [0.0; 4];
    // SAFETY: arrays of four numbers
    unsafe { sys::myb_quat_slerp(a.as_ptr(), b.as_ptr(), t, out.as_mut_ptr()) };
    out
}

/// dq/dt = q * (0, omega) / 2.
pub fn quat_derivative(q: &Quat, omega: &Vec3) -> Quat {
    let mut dq = [0.0; 4];
    // SAFETY: arrays of four and three numbers
    unsafe { sys::myb_quat_derivative(q.as_ptr(), omega.as_ptr(), dq.as_mut_ptr()) };
    dq
}

fn axes_of(axes: &[u32; 3]) -> Result<[c_int; 3]> {
    let mut out = [0; 3];
    for (slot, &axis) in out.iter_mut().zip(axes) {
        *slot = c_int::try_from(axis).or_else(|_| domain())?;
    }
    Ok(out)
}

/// M = R_k(c) R_j(b) R_i(a) for the sequence of axes (i, j, k), each 1, 2 or 3.
pub fn euler_to_matrix(axes: &[u32; 3], angles: &Vec3) -> Result<[f64; 9]> {
    let axes = axes_of(axes)?;
    let mut m = [0.0; 9];
    // SAFETY: arrays of three, three and nine numbers
    check(unsafe { sys::myb_euler_to_matrix(axes.as_ptr(), angles.as_ptr(), m.as_mut_ptr()) })?;
    Ok(m)
}

/// The angles of a matrix for the sequence of axes.
pub fn euler_from_matrix(axes: &[u32; 3], m: &[f64; 9]) -> Result<Vec3> {
    let axes = axes_of(axes)?;
    let mut angles = [0.0; 3];
    // SAFETY: arrays of three, nine and three numbers
    check(unsafe { sys::myb_euler_from_matrix(axes.as_ptr(), m.as_ptr(), angles.as_mut_ptr()) })?;
    Ok(angles)
}

/// The principal moments (ascending) and axes (rows) of a symmetric inertia tensor.
pub fn principal_axes(inertia: &[f64; 9]) -> Result<(Vec3, [f64; 9])> {
    let (mut moments, mut axes) = ([0.0; 3], [0.0; 9]);
    // SAFETY: arrays of nine, three and nine numbers
    check(unsafe { sys::myb_principal_axes(inertia.as_ptr(), moments.as_mut_ptr(), axes.as_mut_ptr()) })?;
    Ok((moments, axes))
}

/// dq/dt and dL/dt of a rigid body in principal axes (torque `None`: free).
pub fn rigid_body_rates(moments: &Vec3, q: &Quat, l: &Vec3, torque: Option<&Vec3>) -> (Quat, Vec3) {
    let (mut dq, mut dl) = ([0.0; 4], [0.0; 3]);
    // SAFETY: arrays of the sizes of the C interface
    unsafe {
        sys::myb_rigid_body_rates(
            moments.as_ptr(),
            q.as_ptr(),
            l.as_ptr(),
            torque.map_or(std::ptr::null(), |n| n.as_ptr()),
            dq.as_mut_ptr(),
            dl.as_mut_ptr(),
        )
    };
    (dq, dl)
}

/// The rates of the angular velocities of a mantle and its fluid core (`myb_core_mantle_rates`);
/// inertia_rate `None`: constant.
#[allow(clippy::too_many_arguments)]
pub fn core_mantle_rates(
    inertia_mantle: &[f64; 9],
    inertia_rate: Option<&[f64; 9]>,
    core: &[f64; 2],
    kv: f64,
    omega_mantle: &Vec3,
    omega_core: &Vec3,
    torque: &Vec3,
) -> Result<(Vec3, Vec3)> {
    let (mut rate_mantle, mut rate_core) = ([0.0; 3], [0.0; 3]);
    // SAFETY: arrays of the sizes of the C interface
    check(unsafe {
        sys::myb_core_mantle_rates(
            inertia_mantle.as_ptr(),
            inertia_rate.map_or(std::ptr::null(), |r| r.as_ptr()),
            core.as_ptr(),
            kv,
            omega_mantle.as_ptr(),
            omega_core.as_ptr(),
            torque.as_ptr(),
            rate_mantle.as_mut_ptr(),
            rate_core.as_mut_ptr(),
        )
    })?;
    Ok((rate_mantle, rate_core))
}

/// The torque of a point mass on a body of inertia diag(moments), r in the body frame.
pub fn torque_point_mass(gm: f64, r: &Vec3, moments: &Vec3) -> Vec3 {
    let mut torque = [0.0; 3];
    // SAFETY: arrays of three numbers
    unsafe { sys::myb_torque_point_mass(gm, r.as_ptr(), moments.as_ptr(), torque.as_mut_ptr()) };
    torque
}

/// -(m_other / m_body) r x a: the torque per unit mass from the field of the body.
pub fn torque_field(mass_ratio: f64, r: &Vec3, a: &Vec3) -> Vec3 {
    let mut torque = [0.0; 3];
    // SAFETY: arrays of three numbers
    unsafe { sys::myb_torque_field(mass_ratio, r.as_ptr(), a.as_ptr(), torque.as_mut_ptr()) };
    torque
}

/// The torque per unit mass of A from the quadrupole-quadrupole potential.
pub fn torque_figure_figure(gm_a: f64, gq_a: &[f64; 9], gq_b: &[f64; 9], r: &Vec3) -> Result<Vec3> {
    let mut torque = [0.0; 3];
    // SAFETY: arrays of the sizes of the C interface
    check(unsafe { sys::myb_torque_figure_figure(gm_a, gq_a.as_ptr(), gq_b.as_ptr(), r.as_ptr(), torque.as_mut_ptr()) })?;
    Ok(torque)
}

/// The Euler-Poinsot motion in closed form (`myb_free_rotation`).
#[derive(Debug, Clone, Copy)]
pub struct FreeRotation(sys::myb_free_rotation);

impl FreeRotation {
    /// From the principal moments, the attitude and the body angular momentum at t = 0.
    pub fn new(moments: &Vec3, q: &Quat, l: &Vec3) -> Result<Self> {
        let mut rotation = sys::myb_free_rotation::default();
        // SAFETY: arrays of the sizes of the C interface
        check(unsafe { sys::myb_free_rotation_init(moments.as_ptr(), q.as_ptr(), l.as_ptr(), &mut rotation) })?;
        Ok(FreeRotation(rotation))
    }

    /// The attitude and the body angular momentum at t.
    pub fn state(&self, t: f64) -> (Quat, Vec3) {
        let (mut q, mut l) = ([0.0; 4], [0.0; 3]);
        // SAFETY: the rotation is initialised; arrays of four and three numbers
        unsafe { sys::myb_free_rotation_state(&self.0, t, q.as_mut_ptr(), l.as_mut_ptr()) };
        (q, l)
    }
}

unsafe extern "C" fn call_torque<F: Fn(f64, &Quat) -> Result<Vec3>>(
    context: *mut c_void,
    t: f64,
    q: *const f64,
    torque: *mut f64,
) -> sys::myb_status {
    // SAFETY: the context is the closure of the call; q holds four numbers, torque three
    let function = unsafe { &*(context as *const F) };
    let q = unsafe { *(q as *const Quat) };
    match function(t, &q) {
        Ok(n) => {
            unsafe { *(torque as *mut Vec3) = n };
            sys::MYB_OK
        }
        Err(error) => error.status.raw(),
    }
}

/// Advances (q, l) from t0 by `steps` steps of h of the symplectic splitting of the order (2,
/// 4, 6, 8), with the torque of a potential the closure gives at (t, q) (`myb_rotation_split`).
#[allow(clippy::too_many_arguments)]
pub fn rotation_split<F: Fn(f64, &Quat) -> Result<Vec3>>(
    moments: &Vec3,
    torque: F,
    order: u32,
    t0: f64,
    h: f64,
    steps: u64,
    q: &Quat,
    l: &Vec3,
) -> Result<(Quat, Vec3)> {
    let order = c_int::try_from(order).or_else(|_| domain())?;
    let steps = i64::try_from(steps).or_else(|_| domain())?;
    let (mut q, mut l) = (*q, *l);
    // SAFETY: the closure lives for the call; arrays of four and three numbers
    check(unsafe {
        sys::myb_rotation_split(
            moments.as_ptr(),
            Some(call_torque::<F>),
            &torque as *const F as *mut c_void,
            order,
            t0,
            h,
            steps,
            q.as_mut_ptr(),
            l.as_mut_ptr(),
        )
    })?;
    Ok((q, l))
}
