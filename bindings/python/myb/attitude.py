# SPDX-License-Identifier: BSD-3-Clause
"""Orientation and rotational dynamics of rigid bodies (myb/attitude.h): quaternions
(w, x, y, z, scalar first), Euler angles, principal axes, the equations of motion, the
gravitational torques, the exact torque-free rotation and the splitting integrator.

A 3 x 3 matrix is nine numbers, the rows one after another, as in the C interface (inputs
may also be three rows); vectors are triples, quaternions quadruples."""
from ._core import _Calls, _check, _doubles, _optional, _settle
from ._native import ffi, lib

__all__ = [
    "FreeRotation", "core_mantle_rates", "euler_from_matrix", "euler_to_matrix",
    "principal_axes", "quat_conjugate", "quat_derivative", "quat_from_matrix",
    "quat_from_rotation_vector", "quat_multiply", "quat_normalize", "quat_rotate",
    "quat_slerp", "quat_to_matrix", "quat_to_rotation_vector", "rigid_body_rates",
    "rotation_split", "torque_field", "torque_figure_figure", "torque_point_mass",
]


def _q(q):
    return _doubles(q, 4, "a quaternion")


def _v(v, what="a vector"):
    return _doubles(v, 3, what)


def quat_multiply(a, b):
    out = ffi.new("double[4]")
    lib.myb_quat_multiply(_q(a), _q(b), out)
    return tuple(out)


def quat_conjugate(q):
    out = ffi.new("double[4]")
    lib.myb_quat_conjugate(_q(q), out)
    return tuple(out)


def quat_normalize(q):
    """(the quaternion scaled to unit norm, its former norm)."""
    made = _q(q)
    norm = lib.myb_quat_normalize(made)
    return tuple(made), norm


def quat_to_matrix(q):
    out = ffi.new("double[9]")
    lib.myb_quat_to_matrix(_q(q), out)
    return tuple(out)


def quat_from_matrix(m):
    out = ffi.new("double[4]")
    lib.myb_quat_from_matrix(_doubles(m, 9, "a matrix"), out)
    return tuple(out)


def quat_rotate(q, v):
    out = ffi.new("double[3]")
    lib.myb_quat_rotate(_q(q), _v(v), out)
    return tuple(out)


def quat_from_rotation_vector(phi):
    out = ffi.new("double[4]")
    lib.myb_quat_from_rotation_vector(_v(phi), out)
    return tuple(out)


def quat_to_rotation_vector(q):
    out = ffi.new("double[3]")
    lib.myb_quat_to_rotation_vector(_q(q), out)
    return tuple(out)


def quat_slerp(a, b, t):
    out = ffi.new("double[4]")
    lib.myb_quat_slerp(_q(a), _q(b), float(t), out)
    return tuple(out)


def quat_derivative(q, omega):
    out = ffi.new("double[4]")
    lib.myb_quat_derivative(_q(q), _v(omega), out)
    return tuple(out)


def _axes(axes):
    values = [int(a) for a in axes]
    if len(values) != 3:
        raise ValueError("a sequence of three axes")
    return ffi.new("int[3]", values)


def euler_to_matrix(axes, angles):
    """The matrix of the angles (a, b, c) about the axes (i, j, k), e.g. (3, 1, 3)."""
    out = ffi.new("double[9]")
    _check(lib.myb_euler_to_matrix(_axes(axes), _v(angles, "angles"), out))
    return tuple(out)


def euler_from_matrix(axes, m):
    out = ffi.new("double[3]")
    _check(lib.myb_euler_from_matrix(_axes(axes), _doubles(m, 9, "a matrix"), out))
    return tuple(out)


def principal_axes(inertia):
    """(moments ascending, axes: the rotation to the principal frame, nine numbers)."""
    moments = ffi.new("double[3]")
    axes = ffi.new("double[9]")
    _check(lib.myb_principal_axes(_doubles(inertia, 9, "inertia"), moments, axes))
    return tuple(moments), tuple(axes)


def rigid_body_rates(moments, q, l, torque=None):
    """(dq/dt, dL/dt) in principal axes."""
    dq = ffi.new("double[4]")
    dl = ffi.new("double[3]")
    lib.myb_rigid_body_rates(_v(moments, "moments"), _q(q), _v(l), _optional(torque, 3, "torque"),
                             dq, dl)
    return tuple(dq), tuple(dl)


def core_mantle_rates(inertia_mantle, inertia_rate, core, kv, omega_mantle, omega_core,
                      torque):
    """(rate of the mantle's angular velocity, rate of the core's) of a body with a fluid
    core; inertia_rate None for a constant inertia."""
    mantle = ffi.new("double[3]")
    rate_core = ffi.new("double[3]")
    _check(lib.myb_core_mantle_rates(_doubles(inertia_mantle, 9, "inertia_mantle"),
                                     _optional(inertia_rate, 9, "inertia_rate"),
                                     _doubles(core, 2, "core"), float(kv),
                                     _v(omega_mantle, "omega_mantle"),
                                     _v(omega_core, "omega_core"),
                                     _v(torque, "torque"), mantle, rate_core))
    return tuple(mantle), tuple(rate_core)


def torque_point_mass(gm, r, moments):
    out = ffi.new("double[3]")
    lib.myb_torque_point_mass(float(gm), _v(r), _v(moments, "moments"), out)
    return tuple(out)


def torque_field(mass_ratio, r, a):
    out = ffi.new("double[3]")
    lib.myb_torque_field(float(mass_ratio), _v(r), _v(a), out)
    return tuple(out)


def torque_figure_figure(gm_a, gq_a, gq_b, r):
    out = ffi.new("double[3]")
    _check(lib.myb_torque_figure_figure(float(gm_a), _doubles(gq_a, 9, "gq_a"),
                                        _doubles(gq_b, 9, "gq_b"), _v(r), out))
    return tuple(out)


class FreeRotation:
    """The exact torque-free rotation (myb_free_rotation) from principal moments, the
    attitude q and the body angular momentum l at t = 0."""

    def __init__(self, moments, q, l):
        self.raw = ffi.new("myb_free_rotation*")
        _check(lib.myb_free_rotation_init(_v(moments, "moments"), _q(q), _v(l), self.raw))

    def state(self, t):
        """(q, l) at t."""
        q = ffi.new("double[4]")
        l = ffi.new("double[3]")
        lib.myb_free_rotation_state(self.raw, float(t), q, l)
        return tuple(q), tuple(l)


def rotation_split(moments, torque, order, t0, h, steps, q, l):
    """(q, l) after `steps` steps of h of the splitting integrator of the given order (2, 4,
    6, 8); torque: a callable f(t, q) -> three numbers in the body frame (myb_rotation_split;
    the free motion between the kicks is exact)."""
    calls = _Calls()

    @ffi.callback("myb_status(void*, double, const double[4], double[3])")
    def function(context, t, attitude, out):
        try:
            out[0:3] = list(_v(torque(t, tuple(attitude[0:4])), "a torque"))
            return lib.MYB_OK
        except BaseException as error:  # noqa: B902 - nothing may cross the C code
            return calls.failed(error)

    attitude = _q(q)
    momentum = _v(l)
    _settle(lib.myb_rotation_split(_v(moments, "moments"), function, ffi.NULL, int(order),
                                   float(t0), float(h), int(steps), attitude, momentum), calls)
    return tuple(attitude), tuple(momentum)

