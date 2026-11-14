// SPDX-License-Identifier: BSD-3-Clause
//! Providers: what the library asks for as a function and a context.
//!
//! * [`StateSource`] — the state of a point at an instant (`myb_state_function` of
//!   `myb/kin.h`): a closure `Fn(Time) -> Result<State>` or a body of an SPK set
//!   ([`SpkBody`]);
//! * [`OrientationSource`] — the rotation of a body-fixed frame at an instant
//!   (`myb_orientation_function` of `myb/frames.h`): a closure, a constant rotation
//!   ([`FixedOrientation`]) or an instrument of a CK file ([`CkOrientation`]);
//! * [`AttitudeSource`] — the attitude of an object at a model time and state
//!   (`myb_attitude_function` of `myb/surface.h`): a closure or one of the laws of the library
//!   ([`Attitude`]).
//!
//! A closure reports a failure by `Err` with the status the library then returns; a panic in
//! it ends the process (it runs inside a function called from C).
//!
//! The traits are `unsafe` to implement: what `as_raw` gives is called by C. Implementing one
//! for a structure of `myb::sys` with its C function (a `myb_cip_orientation` with
//! `myb_orientation_cip`, say) is how a provider of the lower layer joins the upper one.
use std::ffi::c_void;
use std::rc::Rc;

use crate::{sys, Ck, Environment, Ephemeris, Result, Sclk, State, Time, Vec3};

/// The rotation of a frame with its derivative: rows of M and of dM/dt (`myb/frames.h`).
pub type Rotation = sys::myb_rotation;
/// The orientation of a body: the rotation from the inertial axes to its frame, its rotation
/// time and its pole coordinates.
pub type Orientation = sys::myb_orientation;

// ---- States ---------------------------------------------------------------------------------------

/// A provider of states: position (m) and velocity (m/s) at an instant, in the frame and on
/// the time scale agreed with its user.
///
/// # Safety
/// `as_raw` must give a function that may be called with the context, from any call of the
/// library that receives them, for as long as `self` lives and is not moved (any provider of
/// the C interface with its context qualifies: that is how one of `myb::sys` joins in).
pub unsafe trait StateSource {
    /// The function and its context for the C interface, valid as long as `self` is.
    fn as_raw(&self) -> (sys::myb_state_function, *mut c_void);
}

// SAFETY: call_state::<F> reads the closure behind the context, which is `self`
unsafe impl<F: Fn(Time) -> Result<State>> StateSource for F {
    fn as_raw(&self) -> (sys::myb_state_function, *mut c_void) {
        (Some(call_state::<F>), self as *const F as *mut c_void)
    }
}

pub(crate) unsafe extern "C" fn call_state<F: Fn(Time) -> Result<State>>(
    context: *mut c_void,
    t: Time,
    state: *mut f64,
) -> sys::myb_status {
    // SAFETY: the context is the closure as_raw lent; state holds six numbers
    let source = unsafe { &*(context as *const F) };
    match source(t) {
        Ok(y) => {
            unsafe { *(state as *mut State) = y };
            sys::MYB_OK
        }
        Err(error) => error.status.raw(),
    }
}

/// A body of an SPK set as a provider (`myb_spk_body_state`): target relative to center at
/// an instant of the files' scale; frame 0 accepts any frame of the segments.
pub struct SpkBody {
    raw: sys::myb_spk_body,
    _ephemeris: Rc<Ephemeris>,
}

impl SpkBody {
    pub fn new(ephemeris: &Rc<Ephemeris>, target: i32, center: i32, frame: i32) -> Self {
        SpkBody {
            raw: sys::myb_spk_body { spk: ephemeris.as_ptr(), target, center, frame, graph: std::ptr::null() },
            _ephemeris: Rc::clone(ephemeris),
        }
    }

    /// The state at an instant, as the library gives it to its users.
    pub fn state(&self, t: Time) -> Result<State> {
        state_from(self, t)
    }
}

// SAFETY: myb_spk_body_state reads the myb_spk_body inside self, whose set self keeps
unsafe impl StateSource for SpkBody {
    fn as_raw(&self) -> (sys::myb_state_function, *mut c_void) {
        (Some(sys::myb_spk_body_state), &self.raw as *const sys::myb_spk_body as *mut c_void)
    }
}

/// The state a provider gives at an instant, through its C function.
pub fn state_from(source: &(impl StateSource + ?Sized), t: Time) -> Result<State> {
    let (function, context) = source.as_raw();
    let function = function.ok_or(crate::Error::new(crate::Status::Domain))?;
    let mut y = [0.0; 6];
    // SAFETY: the function and context are valid while `source` is borrowed; y holds six numbers
    crate::check(unsafe { function(context, t, y.as_mut_ptr()) })?;
    Ok(y)
}

// ---- Orientations ---------------------------------------------------------------------------------

/// A provider of the orientation of a body-fixed frame at an instant.
/// # Safety
/// As [`StateSource`]: the function and context from `as_raw` must stay callable while `self`
/// lives and is not moved.
pub unsafe trait OrientationSource {
    /// The function and its context for the C interface, valid as long as `self` is.
    fn as_raw(&self) -> (sys::myb_orientation_function, *mut c_void);
}

// SAFETY: call_orientation::<F> reads the closure behind the context, which is `self`
unsafe impl<F: Fn(Time) -> Result<Orientation>> OrientationSource for F {
    fn as_raw(&self) -> (sys::myb_orientation_function, *mut c_void) {
        (Some(call_orientation::<F>), self as *const F as *mut c_void)
    }
}

unsafe extern "C" fn call_orientation<F: Fn(Time) -> Result<Orientation>>(
    context: *mut c_void,
    t: Time,
    orientation: *mut sys::myb_orientation,
) -> sys::myb_status {
    // SAFETY: the context is the closure as_raw lent; the output is valid
    let source = unsafe { &*(context as *const F) };
    match source(t) {
        Ok(o) => {
            unsafe { *orientation = o };
            sys::MYB_OK
        }
        Err(error) => error.status.raw(),
    }
}

/// A constant orientation (`myb_orientation_fixed`).
pub struct FixedOrientation(pub Rotation);

// SAFETY: myb_orientation_fixed reads the rotation inside self
unsafe impl OrientationSource for FixedOrientation {
    fn as_raw(&self) -> (sys::myb_orientation_function, *mut c_void) {
        (Some(sys::myb_orientation_fixed), &self.0 as *const Rotation as *mut c_void)
    }
}

/// An instrument frame of CK files on its spacecraft clock (`myb_orientation_ck`): the
/// argument is the instant on the clock's time scale.
pub struct CkOrientation {
    raw: sys::myb_ck_orientation,
    _ck: Rc<Ck>,
    _sclk: Rc<Sclk>,
}

impl CkOrientation {
    /// base: the frame the pointing must be given in (0: any); tolerance in ticks.
    pub fn new(ck: &Rc<Ck>, sclk: &Rc<Sclk>, instrument: i32, base: i32, tolerance: f64) -> Self {
        CkOrientation {
            raw: sys::myb_ck_orientation { ck: ck.as_ptr(), sclk: sclk.as_ptr(), instrument, base, tolerance },
            _ck: Rc::clone(ck),
            _sclk: Rc::clone(sclk),
        }
    }
}

// SAFETY: myb_orientation_ck reads the structure inside self, whose files and clock self keeps
unsafe impl OrientationSource for CkOrientation {
    fn as_raw(&self) -> (sys::myb_orientation_function, *mut c_void) {
        (Some(sys::myb_orientation_ck), &self.raw as *const sys::myb_ck_orientation as *mut c_void)
    }
}

/// The orientation a provider gives at an instant, through its C function.
pub fn orientation_from(source: &(impl OrientationSource + ?Sized), t: Time) -> Result<Orientation> {
    let (function, context) = source.as_raw();
    let function = function.ok_or(crate::Error::new(crate::Status::Domain))?;
    let mut o = Orientation::default();
    // SAFETY: the function and context are valid while `source` is borrowed
    crate::check(unsafe { function(context, t, &mut o) })?;
    Ok(o)
}

// ---- Attitudes ------------------------------------------------------------------------------------

/// A provider of the attitude of an object: the rows of the matrix that turns the inertial
/// axes of the model into the body frame, at the model time t and the state (x, v).
/// # Safety
/// As [`StateSource`]: the function and context from `as_raw` must stay callable while `self`
/// lives and is not moved.
pub unsafe trait AttitudeSource {
    /// The function and its context for the C interface, valid as long as `self` is.
    fn as_raw(&self) -> (sys::myb_attitude_function, *mut c_void);
}

// SAFETY: call_attitude::<F> reads the closure behind the context, which is `self`
unsafe impl<F: Fn(f64, &Vec3, &Vec3) -> Result<[f64; 9]>> AttitudeSource for F {
    fn as_raw(&self) -> (sys::myb_attitude_function, *mut c_void) {
        (Some(call_attitude::<F>), self as *const F as *mut c_void)
    }
}

unsafe extern "C" fn call_attitude<F: Fn(f64, &Vec3, &Vec3) -> Result<[f64; 9]>>(
    context: *mut c_void,
    t: f64,
    x: *const f64,
    v: *const f64,
    m: *mut f64,
) -> sys::myb_status {
    // SAFETY: the context is the closure as_raw lent; x, v hold 3 numbers, m 9
    let source = unsafe { &*(context as *const F) };
    let (x, v) = unsafe { (*(x as *const Vec3), *(v as *const Vec3)) };
    match source(t, &x, &v) {
        Ok(matrix) => {
            unsafe { *(m as *mut [f64; 9]) = matrix };
            sys::MYB_OK
        }
        Err(error) => error.status.raw(),
    }
}

/// The attitude laws of the library (`myb/surface.h`).
pub struct Attitude(Law);

enum Law {
    Fixed([f64; 9]),
    Orbital { raw: sys::myb_orbital_attitude, _env: Rc<Environment> },
    SunNadir { raw: sys::myb_sun_attitude, _env: Rc<Environment> },
    Spin(sys::myb_spin_attitude),
}

impl Attitude {
    /// A constant matrix, rows (`myb_attitude_fixed`).
    pub fn fixed(m: [f64; 9]) -> Self {
        Attitude(Law::Fixed(m))
    }

    /// x from a body of the environment to the object, z along the orbital angular momentum
    /// (`myb_attitude_orbital`).
    pub fn orbital(env: &Rc<Environment>, body: usize) -> Result<Self> {
        if body >= env.body_count() {
            return crate::domain();
        }
        Ok(Attitude(Law::Orbital { raw: sys::myb_orbital_attitude { env: env.as_ptr(), body }, _env: Rc::clone(env) }))
    }

    /// The yaw steering of a navigation satellite (`myb_attitude_sun_nadir`).
    pub fn sun_nadir(env: &Rc<Environment>, body: usize, source: usize) -> Result<Self> {
        if body >= env.body_count() || source >= env.body_count() {
            return crate::domain();
        }
        Ok(Attitude(Law::SunNadir { raw: sys::myb_sun_attitude { env: env.as_ptr(), body, source }, _env: Rc::clone(env) }))
    }

    /// A uniform rotation about the body z axis (`myb_attitude_spin`): m0 the attitude at
    /// phase + rate t = 0 (rows); rate in rad/s, phase in rad.
    pub fn spin(m0: [f64; 9], rate: f64, phase: f64) -> Self {
        Attitude(Law::Spin(sys::myb_spin_attitude { m0, rate, phase }))
    }
}

// SAFETY: the laws read the structures inside self, whose environments self keeps
unsafe impl AttitudeSource for Attitude {
    fn as_raw(&self) -> (sys::myb_attitude_function, *mut c_void) {
        match &self.0 {
            Law::Fixed(m) => (Some(sys::myb_attitude_fixed), m.as_ptr() as *mut c_void),
            Law::Orbital { raw, .. } => (Some(sys::myb_attitude_orbital), raw as *const _ as *mut c_void),
            Law::SunNadir { raw, .. } => (Some(sys::myb_attitude_sun_nadir), raw as *const _ as *mut c_void),
            Law::Spin(raw) => (Some(sys::myb_attitude_spin), raw as *const _ as *mut c_void),
        }
    }
}

/// The attitude a provider gives at the model time t and the state (x, v).
pub fn attitude_from(source: &(impl AttitudeSource + ?Sized), t: f64, x: &Vec3, v: &Vec3) -> Result<[f64; 9]> {
    let (function, context) = source.as_raw();
    let function = function.ok_or(crate::Error::new(crate::Status::Domain))?;
    let mut m = [0.0; 9];
    // SAFETY: the function and context are valid while `source` is borrowed; x, v hold three
    // numbers, m nine
    crate::check(unsafe { function(context, t, x.as_ptr(), v.as_ptr(), m.as_mut_ptr()) })?;
    Ok(m)
}
