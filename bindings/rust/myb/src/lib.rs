// SPDX-License-Identifier: BSD-3-Clause
//! myb for Rust: the motion of objects in space.
//!
//! Two layers over the C library:
//!
//! * [`sys`] (the crate `myb-sys`) declares every function and structure of the C interface,
//!   by the names of the C headers; it is `unsafe`, as C is.
//! * The types of this crate wrap the path most programs take — kernels and time,
//!   ephemerides, the environment and the force model, propagation, events, the text
//!   messages — with ownership (`Drop`), outcomes as [`Result`] with the status of the
//!   library, and closures for terms, providers and events. A force model keeps alive the
//!   environments and fields its terms point to (they are shared through `Rc`).
//!
//! A closure handed to the library reports a failure by returning `Err`. A panic is not
//! caught: the closure runs inside a function called from C, where a panic cannot unwind,
//! so it ends the program — a broken contract stops at once, as `MYB_CHECK` does in the
//! library.
//!
//! No arithmetic is done here: a number obtained through this crate has the bits of the same
//! calls made in C. States are `[f64; 6]` (position in m, velocity in m/s); the times of a
//! model are seconds from the epoch of its environment.
//!
//! The objects are not `Send`: the library's objects keep working memory, one thread at a
//! time (an environment for another thread is a clone made in C, `myb_env_clone`).
//!
//! ```no_run
//! use std::rc::Rc;
//! let gm = 3.986004418e14;
//! let model = Rc::new(myb::ForceModel::new()?);
//! model.add_point_mass("planet", gm)?;
//! let mut orbit = myb::Orbit::new(&model, &myb::OrbitOptions::default())?;
//! orbit.start(0.0, &[7.0e6, 0.0, 0.0, 0.0, 7.5e3, 0.0], 6000.0)?;
//! let later = orbit.state(6000.0)?;
//! # Ok::<(), myb::Error>(())
//! ```
use std::any::Any;
use std::cell::RefCell;
use std::ffi::{c_char, c_int, c_void, CStr, CString};
use std::fmt;
use std::path::Path;
use std::ptr::{self, NonNull};
use std::rc::Rc;

pub use myb_sys as sys;

/// Position (m) and velocity (m/s).
pub type State = [f64; 6];
pub type Vec3 = [f64; 3];
/// An instant: whole seconds and a fraction past J2000 of its time scale.
pub type Time = sys::myb_time;
/// A calendar label: year, month, day, hour, minute, second.
pub type Calendar = sys::myb_calendar;
/// Classical elements: p (m), e, i, raan, argp, nu (rad).
pub type Elements = sys::myb_elements;

// ---- Outcomes ---------------------------------------------------------------------------------

/// What the library answers when a call fails.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// An argument outside the mathematical domain of the model.
    Domain,
    /// An epoch or argument outside the coverage of the data.
    Range,
    /// The iteration budget was exhausted without convergence.
    NoConvergence,
    /// The step of an integrator fell below the allowed minimum.
    StepMin,
    /// An allocation failed.
    NoMemory,
    /// A size or index computation overflowed.
    Overflow,
    /// A file could not be read.
    Io,
    /// Malformed input data.
    Format,
    /// A named item does not exist.
    NotFound,
    /// A status this crate does not know (a newer library).
    Other(u32),
}

impl Status {
    fn from_raw(code: sys::myb_status) -> Option<Status> {
        Some(match code {
            sys::MYB_OK => return None,
            sys::MYB_E_DOMAIN => Status::Domain,
            sys::MYB_E_RANGE => Status::Range,
            sys::MYB_E_NOCONV => Status::NoConvergence,
            sys::MYB_E_STEP_MIN => Status::StepMin,
            sys::MYB_E_NOMEM => Status::NoMemory,
            sys::MYB_E_OVERFLOW => Status::Overflow,
            sys::MYB_E_IO => Status::Io,
            sys::MYB_E_FORMAT => Status::Format,
            sys::MYB_E_NOT_FOUND => Status::NotFound,
            other => Status::Other(other),
        })
    }

    /// The code of the C interface.
    pub fn raw(self) -> sys::myb_status {
        match self {
            Status::Domain => sys::MYB_E_DOMAIN,
            Status::Range => sys::MYB_E_RANGE,
            Status::NoConvergence => sys::MYB_E_NOCONV,
            Status::StepMin => sys::MYB_E_STEP_MIN,
            Status::NoMemory => sys::MYB_E_NOMEM,
            Status::Overflow => sys::MYB_E_OVERFLOW,
            Status::Io => sys::MYB_E_IO,
            Status::Format => sys::MYB_E_FORMAT,
            Status::NotFound => sys::MYB_E_NOT_FOUND,
            Status::Other(code) => code,
        }
    }

    /// The name the library gives the status: "MYB_E_DOMAIN", ...
    pub fn name(self) -> &'static str {
        // SAFETY: myb_status_name returns a static NUL-terminated string for any code
        unsafe { CStr::from_ptr(sys::myb_status_name(self.raw())) }
            .to_str()
            .unwrap_or("MYB_E_?")
    }
}

/// A failed call: the status and, from parsers, what was wrong and on which line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error {
    pub status: Status,
    pub message: String,
}

impl Error {
    pub fn new(status: Status) -> Self {
        Error { status, message: String::new() }
    }
}

impl From<Status> for Error {
    fn from(status: Status) -> Self {
        Error::new(status)
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.message.is_empty() {
            f.write_str(self.status.name())
        } else {
            write!(f, "{}: {}", self.status.name(), self.message)
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

fn check(code: sys::myb_status) -> Result<()> {
    match Status::from_raw(code) {
        None => Ok(()),
        Some(status) => Err(Error::new(status)),
    }
}

fn check_diag(code: sys::myb_status, diag: &sys::myb_diag) -> Result<()> {
    let Some(status) = Status::from_raw(code) else {
        return Ok(());
    };
    // SAFETY: the library keeps the message NUL-terminated within its 200 characters
    let text = unsafe { CStr::from_ptr(diag.message.as_ptr()) }.to_string_lossy();
    let message = if diag.line > 0 { format!("line {}: {}", diag.line, text) } else { text.into_owned() };
    Err(Error { status, message })
}

fn c_string(text: &str) -> Result<CString> {
    CString::new(text).map_err(|_| Error::new(Status::Domain))
}

fn c_path(path: &Path) -> Result<CString> {
    c_string(path.to_str().ok_or(Error::new(Status::Domain))?)
}

/// A text of the library's writers: they tell the length when the buffer is too small.
fn text(write: impl Fn(*mut c_char, usize, *mut usize) -> sys::myb_status) -> Result<String> {
    let mut length = 0usize;
    let code = write(ptr::null_mut(), 0, &mut length);
    if code != sys::MYB_OK && code != sys::MYB_E_OVERFLOW {
        check(code)?;
    }
    let mut buffer = vec![0u8; length + 1];
    check(write(buffer.as_mut_ptr().cast(), buffer.len(), &mut length))?;
    buffer.truncate(length);
    String::from_utf8(buffer).map_err(|_| Error::new(Status::Format))
}

/// What an object points to and must outlive.
#[derive(Default)]
struct Kept(RefCell<Vec<Box<dyn Any>>>);

impl Kept {
    /// Stores a value for the life of the owner; the pointer stays valid until then.
    fn keep<T: 'static>(&self, value: T) -> *mut T {
        let mut boxed = Box::new(value);
        let pointer: *mut T = &mut *boxed;
        self.0.borrow_mut().push(boxed);
        pointer
    }
}

// ---- Time ---------------------------------------------------------------------------------------

/// An instant from seconds past J2000 of its time scale.
pub fn time_from_seconds(seconds: f64) -> Time {
    // SAFETY: a function of numbers
    unsafe { sys::myb_time_from_seconds(seconds) }
}

/// The seconds past J2000 of an instant, as one number.
pub fn seconds(instant: Time) -> f64 {
    // SAFETY: a function of numbers
    unsafe { sys::myb_time_seconds(instant) }
}

/// A pool of text kernels (`myb/io.h`).
pub struct KernelPool {
    raw: NonNull<sys::myb_kpool>,
}

impl KernelPool {
    pub fn new() -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: a valid place for the pointer; the default allocator
        check(unsafe { sys::myb_kpool_create(ptr::null(), &mut raw) })?;
        Ok(KernelPool { raw: NonNull::new(raw).ok_or(Error::new(Status::NoMemory))? })
    }

    pub fn load_file(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let path = c_path(path.as_ref())?;
        let mut diag = sys::myb_diag::default();
        // SAFETY: the pool is alive; the path is NUL-terminated
        let code = unsafe { sys::myb_kpool_load_file(self.raw.as_ptr(), path.as_ptr(), &mut diag) };
        check_diag(code, &diag)
    }

    pub fn as_ptr(&self) -> *mut sys::myb_kpool {
        self.raw.as_ptr()
    }
}

impl Drop for KernelPool {
    fn drop(&mut self) {
        // SAFETY: the pool was created by myb_kpool_create and is dropped once
        unsafe { sys::myb_kpool_destroy(self.raw.as_ptr()) }
    }
}

/// The time scales of a pool (`myb/time.h`); it keeps nothing of the pool. Scales are named:
/// "TT", "TAI", "TDB", ...; UTC is a calendar label.
pub struct TimeSystem {
    raw: NonNull<sys::myb_timesys>,
}

impl TimeSystem {
    pub fn new(pool: &KernelPool) -> Result<Self> {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the pool is alive for the call
        let code = unsafe { sys::myb_timesys_create(pool.as_ptr(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(TimeSystem { raw: NonNull::new(raw).ok_or(Error::new(Status::NoMemory))? })
    }

    fn scale(&self, name: &str) -> Result<sys::myb_scale> {
        let name = c_string(name)?;
        let mut scale = sys::myb_scale::default();
        // SAFETY: the system is alive; the name is NUL-terminated
        check(unsafe { sys::myb_timesys_scale(self.raw.as_ptr(), name.as_ptr(), &mut scale) })?;
        Ok(scale)
    }

    /// The instant on another scale.
    pub fn convert(&self, instant: Time, from: &str, to: &str) -> Result<Time> {
        let (from, to) = (self.scale(from)?, self.scale(to)?);
        let mut out = Time::default();
        // SAFETY: the system is alive; the scales are its own
        check(unsafe { sys::myb_timesys_convert(self.raw.as_ptr(), instant, from, to, &mut out) })?;
        Ok(out)
    }

    pub fn from_utc(&self, label: &Calendar, to: &str) -> Result<Time> {
        let to = self.scale(to)?;
        let mut out = Time::default();
        // SAFETY: as above
        check(unsafe { sys::myb_timesys_from_utc(self.raw.as_ptr(), label, to, &mut out) })?;
        Ok(out)
    }

    pub fn to_utc(&self, instant: Time, from: &str) -> Result<Calendar> {
        let from = self.scale(from)?;
        let mut out = Calendar::default();
        // SAFETY: as above
        check(unsafe { sys::myb_timesys_to_utc(self.raw.as_ptr(), instant, from, &mut out) })?;
        Ok(out)
    }

    pub fn as_ptr(&self) -> *mut sys::myb_timesys {
        self.raw.as_ptr()
    }
}

impl Drop for TimeSystem {
    fn drop(&mut self) {
        // SAFETY: created by myb_timesys_create, dropped once
        unsafe { sys::myb_timesys_destroy(self.raw.as_ptr()) }
    }
}

// ---- Two bodies ---------------------------------------------------------------------------------

fn pv_of(state: &State) -> sys::myb_pv {
    sys::myb_pv { r: [state[0], state[1], state[2]], v: [state[3], state[4], state[5]] }
}

fn state_of(pv: &sys::myb_pv) -> State {
    [pv.r[0], pv.r[1], pv.r[2], pv.v[0], pv.v[1], pv.v[2]]
}

pub fn state_from_elements(mu: f64, elements: &Elements) -> Result<State> {
    let mut pv = sys::myb_pv::default();
    // SAFETY: valid places for the input and the output
    check(unsafe { sys::myb_pv_from_elements(mu, elements, &mut pv) })?;
    Ok(state_of(&pv))
}

pub fn elements_from_state(mu: f64, state: &State) -> Result<Elements> {
    let mut elements = Elements::default();
    // SAFETY: as above
    check(unsafe { sys::myb_elements_from_pv(mu, &pv_of(state), &mut elements) })?;
    Ok(elements)
}

/// The state a time dt later on the orbit of two bodies.
pub fn kepler_propagate(mu: f64, state: &State, dt: f64) -> Result<State> {
    let mut pv = sys::myb_pv::default();
    // SAFETY: as above
    check(unsafe { sys::myb_kepler_propagate(mu, &pv_of(state), dt, &mut pv) })?;
    Ok(state_of(&pv))
}

// ---- Ephemerides and fields ---------------------------------------------------------------------

/// A set of SPK files (`myb/kin.h`).
pub struct Ephemeris {
    raw: NonNull<sys::myb_spk>,
}

impl Ephemeris {
    pub fn new() -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: a valid place for the pointer
        check(unsafe { sys::myb_spk_create(ptr::null(), &mut raw) })?;
        Ok(Ephemeris { raw: NonNull::new(raw).ok_or(Error::new(Status::NoMemory))? })
    }

    pub fn load_file(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let path = c_path(path.as_ref())?;
        let mut diag = sys::myb_diag::default();
        // SAFETY: the set is alive; the path is NUL-terminated
        let code = unsafe { sys::myb_spk_load_file(self.raw.as_ptr(), path.as_ptr(), &mut diag) };
        check_diag(code, &diag)
    }

    /// Target relative to center at an instant of the files' time scale.
    pub fn state(&self, target: i32, center: i32, instant: Time) -> Result<State> {
        let mut y = [0.0; 6];
        // SAFETY: the set is alive; y holds six numbers
        check(unsafe {
            sys::myb_spk_state(self.raw.as_ptr(), target, center, instant, y.as_mut_ptr(), ptr::null_mut())
        })?;
        Ok(y)
    }

    /// The segments of the set, in the order of loading.
    pub fn segments(&self) -> Vec<sys::myb_spk_segment_info> {
        // SAFETY: the set is alive; the index is within the count
        unsafe {
            (0..sys::myb_spk_segment_count(self.raw.as_ptr()))
                .map(|index| {
                    let mut info = sys::myb_spk_segment_info::default();
                    sys::myb_spk_segment(self.raw.as_ptr(), index, &mut info);
                    info
                })
                .collect()
        }
    }

    pub fn as_ptr(&self) -> *mut sys::myb_spk {
        self.raw.as_ptr()
    }
}

impl Drop for Ephemeris {
    fn drop(&mut self) {
        // SAFETY: created by myb_spk_create, dropped once
        unsafe { sys::myb_spk_destroy(self.raw.as_ptr()) }
    }
}

/// The field of a body in spherical harmonics (`myb/grav.h`).
pub struct GravityField {
    raw: NonNull<sys::myb_gravity_field>,
}

impl GravityField {
    /// c, s: fully normalised coefficients at the index n (n + 1) / 2 + m, to the degree.
    pub fn new(gm: f64, radius: f64, degree: u32, c: &[f64], s: &[f64]) -> Result<Self> {
        let needed = (degree as usize + 1) * (degree as usize + 2) / 2;
        if c.len() < needed || s.len() < needed || degree > i32::MAX as u32 {
            return Err(Error::new(Status::Domain));
        }
        let mut raw = ptr::null_mut();
        // SAFETY: c and s hold the coefficients the degree asks for
        check(unsafe {
            sys::myb_gravity_field_create(gm, radius, degree as c_int, c.as_ptr(), s.as_ptr(), ptr::null(), &mut raw)
        })?;
        Ok(GravityField { raw: NonNull::new(raw).ok_or(Error::new(Status::NoMemory))? })
    }

    /// An ICGEM file (gfc) to a degree; `None`: all of it.
    pub fn load_icgem_file(path: impl AsRef<Path>, max_degree: Option<u32>) -> Result<Self> {
        let path = c_path(path.as_ref())?;
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        let degree = max_degree.map_or(-1, |d| d.min(i32::MAX as u32) as c_int);
        // SAFETY: the path is NUL-terminated; valid places for the outputs
        let code = unsafe {
            sys::myb_gravity_field_load_icgem_file(path.as_ptr(), degree, ptr::null(), &mut raw, &mut diag)
        };
        check_diag(code, &diag)?;
        Ok(GravityField { raw: NonNull::new(raw).ok_or(Error::new(Status::NoMemory))? })
    }

    pub fn degree(&self) -> i32 {
        // SAFETY: the field is alive
        unsafe { sys::myb_gravity_field_degree(self.raw.as_ptr()) }
    }

    pub fn gm(&self) -> f64 {
        // SAFETY: the field is alive
        unsafe { sys::myb_gravity_field_gm(self.raw.as_ptr()) }
    }

    pub fn radius(&self) -> f64 {
        // SAFETY: the field is alive
        unsafe { sys::myb_gravity_field_radius(self.raw.as_ptr()) }
    }

    pub fn as_ptr(&self) -> *mut sys::myb_gravity_field {
        self.raw.as_ptr()
    }
}

impl Drop for GravityField {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_gravity_field_destroy(self.raw.as_ptr()) }
    }
}

// ---- The environment ----------------------------------------------------------------------------

/// A body of an environment: where it is and what it is to the terms of a model.
#[derive(Default)]
pub struct Body {
    /// The caller's identifier (a NAIF code, say).
    pub id: i64,
    /// m^3 s^-2.
    pub gm: f64,
    /// m: the sphere that casts shadows and bounds the disk of a source.
    pub radius: f64,
    /// W: the body as a source of radiation.
    pub luminosity: f64,
    /// Limb darkening u1, u2 of the disk of a source.
    pub limb: [f64; 2],
}

type Provider = dyn Fn(f64) -> Result<State>;

/// The bodies of a scenario (`myb/env.h`): what the terms of a force model ask for.
pub struct Environment {
    raw: NonNull<sys::myb_env>,
    kept: Kept,
}

impl Environment {
    /// epoch: the instant of the model time 0.
    pub fn new(epoch: Time) -> Result<Self> {
        let options = sys::myb_env_options { epoch, ..Default::default() };
        let mut raw = ptr::null_mut();
        // SAFETY: valid options (no time system: one scale) and a place for the pointer
        check(unsafe { sys::myb_env_create(&options, ptr::null(), &mut raw) })?;
        Ok(Environment { raw: NonNull::new(raw).ok_or(Error::new(Status::NoMemory))?, kept: Kept::default() })
    }

    fn add(&self, body: &Body, state: sys::myb_state_function, context: *mut c_void) -> Result<usize> {
        let made = sys::myb_env_body {
            id: body.id,
            gm: body.gm,
            state,
            state_context: context,
            radius: body.radius,
            luminosity: body.luminosity,
            limb: body.limb,
            ..Default::default()
        };
        let mut index = 0usize;
        // SAFETY: the environment is alive; the provider and its context live as long as it
        // does (they are in self.kept)
        check(unsafe { sys::myb_env_add_body(self.raw.as_ptr(), &made, &mut index) })?;
        Ok(index)
    }

    /// The body at the origin of the environment. Returns its index.
    pub fn add_origin(&self, body: &Body) -> Result<usize> {
        self.add(body, None, ptr::null_mut())
    }

    /// A body whose state relative to the origin is that of `target` relative to `center`
    /// in an SPK set (frame 0: whatever frame the segments are in).
    pub fn add_from_ephemeris(
        &self,
        body: &Body,
        ephemeris: &Rc<Ephemeris>,
        target: i32,
        center: i32,
        frame: i32,
    ) -> Result<usize> {
        let context = self.kept.keep(sys::myb_spk_body {
            spk: ephemeris.as_ptr(),
            target,
            center,
            frame,
            graph: std::ptr::null(),
        });
        self.kept.keep(Rc::clone(ephemeris));
        self.add(body, Some(sys::myb_spk_body_state), context.cast())
    }

    /// A body whose state relative to the origin a closure gives, at an instant in seconds
    /// past J2000 of the model's scale.
    pub fn add_with(&self, body: &Body, state: impl Fn(f64) -> Result<State> + 'static) -> Result<usize> {
        let provider: Box<Provider> = Box::new(state);
        let context = self.kept.keep(provider);
        self.add(body, Some(call_provider), context.cast())
    }

    pub fn find(&self, id: i64) -> Result<usize> {
        let mut index = 0usize;
        // SAFETY: the environment is alive
        check(unsafe { sys::myb_env_find(self.raw.as_ptr(), id, &mut index) })?;
        Ok(index)
    }

    /// The state of a body relative to the origin at the model time t.
    pub fn state(&self, t: f64, body: usize) -> Result<State> {
        if body >= self.body_count() {
            return Err(Error::new(Status::Domain));
        }
        let mut y = [0.0; 6];
        // SAFETY: the environment is alive; the body is one of its own
        check(unsafe { sys::myb_env_state(self.raw.as_ptr(), t, body, y.as_mut_ptr()) })?;
        Ok(y)
    }

    pub fn body_count(&self) -> usize {
        // SAFETY: the environment is alive
        unsafe { sys::myb_env_body_count(self.raw.as_ptr()) }
    }

    pub fn as_ptr(&self) -> *mut sys::myb_env {
        self.raw.as_ptr()
    }
}

impl Drop for Environment {
    fn drop(&mut self) {
        // SAFETY: created by myb_env_create, dropped once; what it points to is dropped after
        unsafe { sys::myb_env_destroy(self.raw.as_ptr()) }
    }
}

unsafe extern "C" fn call_provider(context: *mut c_void, t: Time, state: *mut f64) -> sys::myb_status {
    // SAFETY: the context is the Box<Provider> the environment keeps; state holds six numbers
    let provider = unsafe { &*(context as *const Box<Provider>) };
    let seconds = unsafe { sys::myb_time_seconds(t) };
    match provider(seconds) {
        Ok(y) => {
            unsafe { *(state as *mut State) = y };
            sys::MYB_OK
        }
        Err(error) => error.status.raw(),
    }
}

// ---- The force model ----------------------------------------------------------------------------

/// The object the forces act on.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Object {
    /// kg.
    pub mass: f64,
    /// m^2: the cross-section.
    pub area: f64,
    /// The radiation pressure coefficient, 1 + reflectivity.
    pub cr: f64,
    /// The drag coefficient.
    pub cd: f64,
}

impl Object {
    fn raw(&self) -> sys::myb_object {
        sys::myb_object { mass: self.mass, area: self.area, cr: self.cr, cd: self.cd, ..Default::default() }
    }
}

/// A layered exponential atmosphere over an ellipsoid: (base height m, density kg/m^3 at the
/// base, scale height m) of each layer, by height.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ExponentialAtmosphere {
    pub equatorial_radius: f64,
    pub flattening: f64,
    pub layers: Vec<(f64, f64, f64)>,
}

type Term = dyn Fn(f64, &Vec3, &Vec3) -> Result<Vec3>;

/// The forces on one object (`myb/model.h`) with the ready-made terms of `myb/env.h`. The
/// model owns the structures of its terms and keeps what they point to alive.
pub struct ForceModel {
    raw: NonNull<sys::myb_force_model>,
    kept: Kept,
}

impl ForceModel {
    pub fn new() -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: a valid place for the pointer
        check(unsafe { sys::myb_force_model_create(ptr::null(), &mut raw) })?;
        Ok(ForceModel { raw: NonNull::new(raw).ok_or(Error::new(Status::NoMemory))?, kept: Kept::default() })
    }

    /// Any term of the C interface. Returns its index.
    ///
    /// # Safety
    /// `function` and `jacobian` must be callable with `context` for as long as the model
    /// lives, as the C interface asks.
    pub unsafe fn add_raw(
        &self,
        name: &str,
        function: sys::myb_term_function,
        jacobian: sys::myb_term_jacobian,
        context: *mut c_void,
    ) -> Result<usize> {
        let name = c_string(name)?;
        let mut index = 0usize;
        // SAFETY: the model is alive; the rest is the caller's promise
        check(unsafe { sys::myb_force_model_add(self.raw.as_ptr(), name.as_ptr(), function, jacobian, context, &mut index) })?;
        Ok(index)
    }

    /// The attraction of a point mass at the origin.
    pub fn add_point_mass(&self, name: &str, gm: f64) -> Result<usize> {
        let term = self.kept.keep(sys::myb_point_mass_term { gm });
        // SAFETY: the term lives in self.kept
        unsafe { self.add_raw(name, Some(sys::myb_term_point_mass), Some(sys::myb_term_point_mass_jacobian), term.cast()) }
    }

    fn body_of(env: &Environment, body: usize) -> Result<usize> {
        if body < env.body_count() {
            Ok(body)
        } else {
            Err(Error::new(Status::Domain))
        }
    }

    /// The attraction of a body of the environment away from the origin.
    pub fn add_third_body(&self, name: &str, env: &Rc<Environment>, body: usize) -> Result<usize> {
        let body = Self::body_of(env, body)?;
        let term = self.kept.keep(sys::myb_third_body_term { env: env.as_ptr(), body });
        self.kept.keep(Rc::clone(env));
        // SAFETY: the term and the environment live as long as the model
        unsafe { self.add_raw(name, Some(sys::myb_term_third_body), Some(sys::myb_term_third_body_jacobian), term.cast()) }
    }

    /// The field of a body in harmonics of the degrees n_min ..= n_max (`None`: all it has).
    pub fn add_field(
        &self,
        name: &str,
        env: &Rc<Environment>,
        body: usize,
        field: &Rc<GravityField>,
        n_min: u32,
        n_max: Option<u32>,
    ) -> Result<usize> {
        let body = Self::body_of(env, body)?;
        let mut work = ptr::null_mut();
        // SAFETY: a valid place for the pointer
        check(unsafe { sys::myb_gravity_work_create(field.degree(), ptr::null(), &mut work) })?;
        self.kept.keep(Work(work));
        let term = self.kept.keep(sys::myb_field_term {
            env: env.as_ptr(),
            body,
            field: field.as_ptr(),
            work,
            n_min: n_min.min(i32::MAX as u32) as c_int,
            n_max: n_max.map_or(-1, |n| n.min(i32::MAX as u32) as c_int),
            ..Default::default()
        });
        self.kept.keep(Rc::clone(env));
        self.kept.keep(Rc::clone(field));
        // SAFETY: the term, its work space, the environment and the field live with the model
        unsafe { self.add_raw(name, Some(sys::myb_term_field), Some(sys::myb_term_field_jacobian), term.cast()) }
    }

    /// The pressure of the light of a body of the environment on the object, in the shadows
    /// of the occulting bodies; c: the speed of light.
    pub fn add_radiation(
        &self,
        name: &str,
        env: &Rc<Environment>,
        source: usize,
        object: &Object,
        c: f64,
        occulters: &[usize],
    ) -> Result<usize> {
        let source = Self::body_of(env, source)?;
        let mut term = sys::myb_radiation_term { env: env.as_ptr(), source, c, ..Default::default() };
        if occulters.len() > term.occulters.len() {
            return Err(Error::new(Status::Domain));
        }
        for (slot, &body) in term.occulters.iter_mut().zip(occulters) {
            *slot = Self::body_of(env, body)?;
        }
        term.occulter_count = occulters.len();
        term.object = self.kept.keep(object.raw());
        let term = self.kept.keep(term);
        self.kept.keep(Rc::clone(env));
        // SAFETY: the term, the object and the environment live with the model
        unsafe { self.add_raw(name, Some(sys::myb_term_radiation), Some(sys::myb_term_radiation_jacobian), term.cast()) }
    }

    /// The drag of an exponential atmosphere that turns with a body of the environment.
    pub fn add_drag(
        &self,
        name: &str,
        env: &Rc<Environment>,
        body: usize,
        air: &ExponentialAtmosphere,
        object: &Object,
    ) -> Result<usize> {
        let body = Self::body_of(env, body)?;
        let mut made = sys::myb_exponential_atmosphere::default();
        if air.layers.is_empty() || air.layers.len() > made.base.len() {
            return Err(Error::new(Status::Domain));
        }
        made.ellipsoid = sys::myb_ellipsoid { a: air.equatorial_radius, f: air.flattening };
        made.count = air.layers.len();
        for (k, &(base, density, scale)) in air.layers.iter().enumerate() {
            made.base[k] = base;
            made.density[k] = density;
            made.scale[k] = scale;
        }
        let term = self.kept.keep(sys::myb_drag_term {
            env: env.as_ptr(),
            body,
            density: Some(sys::myb_density_exponential),
            density_context: self.kept.keep(made).cast(),
            object: self.kept.keep(object.raw()),
            ..Default::default()
        });
        self.kept.keep(Rc::clone(env));
        // SAFETY: the term, the atmosphere, the object and the environment live with the model
        unsafe { self.add_raw(name, Some(sys::myb_term_drag), Some(sys::myb_term_drag_jacobian), term.cast()) }
    }

    /// The post-Newtonian correction of a central body of the environment (Schwarzschild),
    /// in general relativity (beta = gamma = 1).
    pub fn add_schwarzschild(&self, name: &str, env: &Rc<Environment>, central: usize, c: f64) -> Result<usize> {
        let central = Self::body_of(env, central)?;
        let term = self.kept.keep(sys::myb_relativity_term {
            env: env.as_ptr(),
            ppn: sys::myb_ppn { c, beta: 1.0, gamma: 1.0 },
            central,
            external: central,
        });
        self.kept.keep(Rc::clone(env));
        // SAFETY: the term and the environment live with the model
        unsafe { self.add_raw(name, Some(sys::myb_term_schwarzschild), Some(sys::myb_term_schwarzschild_jacobian), term.cast()) }
    }

    /// A term of the caller's own: the acceleration at the model time t and the state
    /// (x, v). Its derivatives are taken by differences of the model.
    pub fn add_term(&self, name: &str, term: impl Fn(f64, &Vec3, &Vec3) -> Result<Vec3> + 'static) -> Result<usize> {
        let boxed: Box<Term> = Box::new(term);
        let context = self.kept.keep(boxed);
        // SAFETY: call_term reads the Box<Term> that self.kept holds
        unsafe { self.add_raw(name, Some(call_term), None, context.cast()) }
    }

    /// The acceleration of the object at the model time t and the state (x, v).
    pub fn acceleration(&self, t: f64, x: &Vec3, v: &Vec3) -> Result<Vec3> {
        let mut a = [0.0; 3];
        // SAFETY: the model is alive; the arrays hold three numbers each
        check(unsafe { sys::myb_force_model_acceleration(self.raw.as_ptr(), t, x.as_ptr(), v.as_ptr(), a.as_mut_ptr()) })?;
        Ok(a)
    }

    pub fn term_count(&self) -> usize {
        // SAFETY: the model is alive
        unsafe { sys::myb_force_model_count(self.raw.as_ptr()) }
    }

    pub fn term_name(&self, index: usize) -> Option<String> {
        if index >= self.term_count() {
            return None;
        }
        // SAFETY: the index is within the count; the name is a NUL-terminated string of the model
        Some(unsafe { CStr::from_ptr(sys::myb_force_model_name(self.raw.as_ptr(), index)) }.to_string_lossy().into_owned())
    }

    /// Switches a term on or off; false if there is no such term.
    pub fn set_enabled(&self, index: usize, enabled: bool) -> bool {
        if index >= self.term_count() {
            return false;
        }
        // SAFETY: the index is within the count
        unsafe { sys::myb_force_model_set_enabled(self.raw.as_ptr(), index, c_int::from(enabled)) };
        true
    }

    pub fn as_ptr(&self) -> *mut sys::myb_force_model {
        self.raw.as_ptr()
    }
}

impl Drop for ForceModel {
    fn drop(&mut self) {
        // SAFETY: created by myb_force_model_create, dropped once; the fields — what the
        // terms point to — are dropped after this
        unsafe { sys::myb_force_model_destroy(self.raw.as_ptr()) }
    }
}

/// The work space of a field term.
struct Work(*mut sys::myb_gravity_work);

impl Drop for Work {
    fn drop(&mut self) {
        // SAFETY: created by myb_gravity_work_create, dropped once
        unsafe { sys::myb_gravity_work_destroy(self.0) }
    }
}

unsafe extern "C" fn call_term(context: *mut c_void, t: f64, x: *const f64, v: *const f64, a: *mut f64) -> sys::myb_status {
    // SAFETY: the context is the Box<Term> the model keeps; x, v, a hold three numbers each
    let term = unsafe { &*(context as *const Box<Term>) };
    let (x, v) = unsafe { (*(x as *const Vec3), *(v as *const Vec3)) };
    match term(t, &x, &v) {
        Ok(acceleration) => {
            unsafe { *(a as *mut Vec3) = acceleration };
            sys::MYB_OK
        }
        Err(error) => error.status.raw(),
    }
}

// ---- Propagation --------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Integrator {
    #[default]
    Ias15,
    Dop853,
}

/// How an orbit is integrated; zeros stand for the defaults of the library.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct OrbitOptions {
    /// The precision of IAS15 (0: 1e-10).
    pub epsilon: f64,
    /// The largest step in s (0: unlimited).
    pub h_max: f64,
    /// The budget of steps of one call (0: 10^9).
    pub max_steps: u64,
    /// Integrate the transition matrix too, for [`Orbit::partials`].
    pub variational: bool,
    pub integrator: Integrator,
    /// DOP853: the relative tolerance (0: 1e-12) and the absolute one (0: 1e-9).
    pub rtol: f64,
    pub atol: f64,
}

/// Which zeros of an event function count.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Direction {
    #[default]
    Any,
    Rising,
    Falling,
}

/// One object under a force model (`myb/orbit.h`): states at the times asked for, events.
pub struct Orbit {
    raw: NonNull<sys::myb_orbit>,
    _model: Rc<ForceModel>,
}

impl Orbit {
    pub fn new(model: &Rc<ForceModel>, options: &OrbitOptions) -> Result<Self> {
        let made = sys::myb_orbit_options {
            epsilon: options.epsilon,
            h_max: options.h_max,
            max_steps: options.max_steps.min(i64::MAX as u64) as i64,
            variational: c_int::from(options.variational),
            integrator: match options.integrator {
                Integrator::Ias15 => sys::MYB_ORBIT_IAS15,
                Integrator::Dop853 => sys::MYB_ORBIT_DOP853,
            },
            rtol: options.rtol,
            atol: options.atol,
        };
        let mut raw = ptr::null_mut();
        // SAFETY: the model is alive and kept alive by the orbit
        check(unsafe { sys::myb_orbit_create(model.as_ptr(), &made, ptr::null(), &mut raw) })?;
        Ok(Orbit { raw: NonNull::new(raw).ok_or(Error::new(Status::NoMemory))?, _model: Rc::clone(model) })
    }

    /// The object is at `state` at t0 and will be asked about up to t_end (before t0:
    /// backwards).
    pub fn start(&mut self, t0: f64, state: &State, t_end: f64) -> Result<()> {
        // SAFETY: the orbit is alive; the state holds six numbers
        check(unsafe { sys::myb_orbit_start(self.raw.as_ptr(), t0, state.as_ptr(), t_end) })
    }

    pub fn state(&mut self, t: f64) -> Result<State> {
        let mut y = [0.0; 6];
        // SAFETY: as above
        check(unsafe { sys::myb_orbit_state(self.raw.as_ptr(), t, y.as_mut_ptr()) })?;
        Ok(y)
    }

    /// The states at times in the order of the propagation.
    pub fn states(&mut self, times: &[f64]) -> Result<Vec<State>> {
        let mut states = vec![[0.0; 6]; times.len()];
        // SAFETY: `states` is times.len() arrays of six numbers in a row
        check(unsafe { sys::myb_orbit_states(self.raw.as_ptr(), times.as_ptr(), times.len(), states.as_mut_ptr().cast()) })?;
        Ok(states)
    }

    /// The next zero of g(t, state) up to t_limit: its time and the state there, or `None`
    /// if there is none before t_limit.
    pub fn next_event<G: Fn(f64, &State) -> f64>(&mut self, g: G, t_limit: f64, direction: Direction) -> Result<Option<(f64, State)>> {
        let direction = match direction {
            Direction::Any => 0,
            Direction::Rising => 1,
            Direction::Falling => -1,
        };
        let mut found: c_int = 0;
        let mut t_event = 0.0;
        let mut y_event = [0.0; 6];
        // SAFETY: g outlives the call, which is the only user of the pointer to it
        check(unsafe {
            sys::myb_orbit_event(
                self.raw.as_ptr(),
                Some(call_event::<G>),
                &g as *const G as *mut c_void,
                t_limit,
                direction,
                &mut found,
                &mut t_event,
                y_event.as_mut_ptr(),
            )
        })?;
        Ok((found != 0).then_some((t_event, y_event)))
    }

    /// An impulse: the velocity changes by dv where the object is now.
    pub fn kick(&mut self, dv: &Vec3) -> Result<()> {
        // SAFETY: the orbit is alive; dv holds three numbers
        check(unsafe { sys::myb_orbit_kick(self.raw.as_ptr(), dv.as_ptr()) })
    }

    /// The state at t and the transition matrix from the start, by rows
    /// (`OrbitOptions::variational` must be set).
    pub fn partials(&mut self, t: f64) -> Result<(State, [f64; 36])> {
        let mut y = [0.0; 6];
        let mut phi = [0.0; 36];
        // SAFETY: the outputs hold 6 and 36 numbers; no parameters are asked for
        check(unsafe { sys::myb_orbit_partials(self.raw.as_ptr(), t, y.as_mut_ptr(), phi.as_mut_ptr(), ptr::null_mut()) })?;
        Ok((y, phi))
    }

    /// The time the propagation has reached.
    pub fn time(&self) -> f64 {
        // SAFETY: the orbit is alive
        unsafe { sys::myb_orbit_time(self.raw.as_ptr()) }
    }

    pub fn steps(&self) -> i64 {
        // SAFETY: the orbit is alive
        unsafe { sys::myb_orbit_steps(self.raw.as_ptr()) }
    }

    pub fn as_ptr(&self) -> *mut sys::myb_orbit {
        self.raw.as_ptr()
    }
}

impl Drop for Orbit {
    fn drop(&mut self) {
        // SAFETY: created by myb_orbit_create, dropped once, before the model it uses
        unsafe { sys::myb_orbit_destroy(self.raw.as_ptr()) }
    }
}

unsafe extern "C" fn call_event<G: Fn(f64, &State) -> f64>(context: *mut c_void, t: f64, y: *const f64) -> f64 {
    // SAFETY: the context is the closure next_event borrowed for the call; y holds a state
    let g = unsafe { &*(context as *const G) };
    let state = unsafe { *(y as *const State) };
    g(t, &state)
}

// ---- Messages -----------------------------------------------------------------------------------

/// The notation of a message of the CCSDS.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Notation {
    /// KEYWORD = value lines.
    #[default]
    Kvn,
    Xml,
}

/// A segment of an orbit ephemeris message: states in m and m/s at calendar labels.
#[derive(Debug, Clone, Default)]
pub struct OemSegment {
    pub object_name: String,
    pub object_id: String,
    pub center_name: String,
    pub ref_frame: String,
    pub time_system: String,
    /// The interpolation and its degree, if the message is to name them.
    pub interpolation: Option<(String, u32)>,
    pub epochs: Vec<Calendar>,
    pub states: Vec<State>,
}

/// An orbit ephemeris message as text.
pub fn write_oem(originator: &str, creation_date: &Calendar, segments: &[OemSegment], notation: Notation) -> Result<String> {
    struct Texts {
        fixed: [CString; 5],
        interpolation: Option<CString>,
    }
    let mut texts = Vec::with_capacity(segments.len());
    for segment in segments {
        if segment.epochs.len() != segment.states.len() {
            return Err(Error::new(Status::Domain));
        }
        texts.push(Texts {
            fixed: [
                c_string(&segment.object_name)?,
                c_string(&segment.object_id)?,
                c_string(&segment.center_name)?,
                c_string(&segment.ref_frame)?,
                c_string(&segment.time_system)?,
            ],
            interpolation: segment.interpolation.as_ref().map(|(name, _)| c_string(name)).transpose()?,
        });
    }
    let made: Vec<sys::myb_oem_segment_data> = segments
        .iter()
        .zip(&texts)
        .map(|(segment, text)| sys::myb_oem_segment_data {
            object_name: text.fixed[0].as_ptr(),
            object_id: text.fixed[1].as_ptr(),
            center_name: text.fixed[2].as_ptr(),
            ref_frame: text.fixed[3].as_ptr(),
            time_system: text.fixed[4].as_ptr(),
            interpolation: text.interpolation.as_ref().map_or(ptr::null(), |name| name.as_ptr()),
            interpolation_degree: segment.interpolation.as_ref().map_or(0, |&(_, degree)| degree.min(i32::MAX as u32) as c_int),
            count: segment.epochs.len(),
            epochs: segment.epochs.as_ptr(),
            states: segment.states.as_ptr().cast(),
            ..Default::default()
        })
        .collect();
    let originator = c_string(originator)?;
    let data = sys::myb_oem_data {
        originator: originator.as_ptr(),
        creation_date: *creation_date,
        segment_count: made.len(),
        segments: made.as_ptr(),
    };
    let notation = match notation {
        Notation::Kvn => sys::MYB_ODM_KVN,
        Notation::Xml => sys::MYB_ODM_XML,
    };
    // SAFETY: every pointer of `data` is into `texts`, `made` and `segments`, alive here
    text(|buffer, size, length| unsafe { sys::myb_oem_write_as(&data, notation, buffer, size, length) })
}

/// An orbit ephemeris message read from text in either notation (`myb/ccsds.h`).
pub struct Oem {
    raw: NonNull<sys::myb_oem>,
}

impl Oem {
    pub fn parse(text: &str) -> Result<Self> {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the text is text.len() bytes; valid places for the outputs
        let code = unsafe { sys::myb_oem_load(text.as_ptr().cast(), text.len(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Oem { raw: NonNull::new(raw).ok_or(Error::new(Status::NoMemory))? })
    }

    pub fn segment_count(&self) -> usize {
        // SAFETY: the message is alive
        unsafe { sys::myb_oem_segment_count(self.raw.as_ptr()) }
    }

    /// The number of states of a segment (0 if there is no such segment).
    pub fn state_count(&self, segment: usize) -> usize {
        if segment >= self.segment_count() {
            return 0;
        }
        // SAFETY: the segment is one of the message
        unsafe { sys::myb_oem_state_count(self.raw.as_ptr(), segment) }
    }

    /// A keyword of the metadata of a segment, as written.
    pub fn meta(&self, segment: usize, key: &str) -> Option<String> {
        if segment >= self.segment_count() {
            return None;
        }
        let key = CString::new(key).ok()?;
        // SAFETY: the segment is one of the message; the key is NUL-terminated
        let value = unsafe { sys::myb_oem_meta(self.raw.as_ptr(), segment, key.as_ptr()) };
        // SAFETY: a value, if any, is a NUL-terminated string that lives with the message
        (!value.is_null()).then(|| unsafe { CStr::from_ptr(value) }.to_string_lossy().into_owned())
    }

    /// The epoch and the state of a line of a segment.
    pub fn state(&self, segment: usize, index: usize) -> Option<(Calendar, State)> {
        if index >= self.state_count(segment) {
            return None;
        }
        let mut epoch = Calendar::default();
        let mut y = [0.0; 6];
        // SAFETY: the line is one of the segment; the outputs are valid
        unsafe { sys::myb_oem_state(self.raw.as_ptr(), segment, index, &mut epoch, y.as_mut_ptr(), ptr::null_mut()) };
        Some((epoch, y))
    }
}

impl Drop for Oem {
    fn drop(&mut self) {
        // SAFETY: created by myb_oem_load, dropped once
        unsafe { sys::myb_oem_destroy(self.raw.as_ptr()) }
    }
}
