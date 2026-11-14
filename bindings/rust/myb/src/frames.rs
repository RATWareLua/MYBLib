// SPDX-License-Identifier: BSD-3-Clause
//! Reference frames (`myb/frames.h`): rotations with their derivatives, the CIO-based model of
//! the Earth ([`Cip`]) with the observed Earth orientation ([`Eop`]), the equinox-based
//! systems ([`Equinox`]), the rotation models of the WGCCRE ([`Wgccre`]), binary PCK files
//! ([`Pck`]), fixed frames of frame kernels ([`tkframe`]) and the graph of frames
//! ([`FrameGraph`], with SPK states across frames: [`crate::Ephemeris::state_in`]).
//!
//! Every model has its array forms (`*_many`): one crossing into C for a whole series, the
//! same bits as the single calls. The models that give the orientation of a body are also
//! providers ([`CipOrientation`], [`PckOrientation`], [`WgccreOrientation`],
//! [`GraphOrientation`]) for a [`crate::System`], a graph or [`crate::orientation_many`].
use std::ffi::{c_char, c_void, CString};
use std::path::Path;
use std::ptr::{self, NonNull};
use std::rc::Rc;

use crate::many::check_at;
use crate::{
    array_text, c_path, c_string, check, check_diag, created, domain, sys, Ephemeris, KernelPool, Kept, OrientationSource, Result,
    Rotation, State, Time, TimeSystem, Vec3,
};

/// Earth orientation parameters at an instant: polar motion x_p, y_p (rad), UT1 - UTC (s),
/// the excess length of day (s), the celestial pole offsets dX, dY (rad).
pub type EopValues = sys::myb_eop_values;

fn times_check<T>(a: &[T], b: &[T]) -> Result<()> {
    if a.len() == b.len() {
        Ok(())
    } else {
        domain()
    }
}

// ---- Rotations ------------------------------------------------------------------------------------

/// M x for a position.
pub fn rotation_apply_position(rotation: &Rotation, position: &Vec3) -> Vec3 {
    let mut out = [0.0; 3];
    // SAFETY: arrays of three numbers
    unsafe { sys::myb_rotation_apply_position(rotation, position.as_ptr(), out.as_mut_ptr()) };
    out
}

/// (M r, M v + dM/dt r) for a state.
pub fn rotation_apply_state(rotation: &Rotation, state: &State) -> State {
    let mut out = [0.0; 6];
    // SAFETY: arrays of six numbers
    unsafe { sys::myb_rotation_apply_state(rotation, state.as_ptr(), out.as_mut_ptr()) };
    out
}

/// `second` after `first` (A -> B by first, B -> C by second).
pub fn rotation_compose(second: &Rotation, first: &Rotation) -> Rotation {
    let mut out = Rotation::default();
    // SAFETY: valid places for the inputs and the output
    unsafe { sys::myb_rotation_compose(second, first, &mut out) };
    out
}

/// The inverse rotation, with its derivative.
pub fn rotation_invert(rotation: &Rotation) -> Rotation {
    let mut out = Rotation::default();
    // SAFETY: valid places for the input and the output
    unsafe { sys::myb_rotation_invert(rotation, &mut out) };
    out
}

/// `M = R3(angles[2]) R1(angles[1]) R3(angles[0])` and its derivative for the angle rates.
pub fn rotation_euler313(angles: &Vec3, rates: &Vec3) -> Rotation {
    let mut out = Rotation::default();
    // SAFETY: arrays of three numbers
    unsafe { sys::myb_rotation_euler313(angles.as_ptr(), rates.as_ptr(), &mut out) };
    out
}

/// The rotation celestial -> intermediate from X, Y, s and their rates (`None`: zero rates).
pub fn cip_c2i(xys: &Vec3, rates: Option<&Vec3>) -> Rotation {
    let mut out = Rotation::default();
    // SAFETY: arrays of three numbers
    unsafe { sys::myb_cip_c2i(xys.as_ptr(), rates.map_or(ptr::null(), |r| r.as_ptr()), &mut out) };
    out
}

/// A fixed frame of a frames kernel: the rotation from its relative frame (constant) and the
/// name of that frame (`myb_tkframe`).
pub fn tkframe(pool: &KernelPool, id: i64) -> Result<(Rotation, String)> {
    let mut rotation = Rotation::default();
    let mut relative = [0 as c_char; 33];
    let mut diag = sys::myb_diag::default();
    // SAFETY: the pool is alive; relative holds 33 characters
    let code = unsafe { sys::myb_tkframe(pool.as_ptr(), id, &mut rotation, relative.as_mut_ptr(), &mut diag) };
    check_diag(code, &diag)?;
    Ok((rotation, array_text(&relative)))
}

// ---- The CIO-based model --------------------------------------------------------------------------

/// The celestial intermediate pole and origin of a body (`myb_cip`), from the kernel
/// variables under a prefix ("IAU2006A" of data/frames/iau2006-cip.tk).
pub struct Cip {
    raw: NonNull<sys::myb_cip>,
}

impl Cip {
    pub fn new(pool: &KernelPool, prefix: &str) -> Result<Self> {
        let prefix = c_string(prefix)?;
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the pool is alive for the call
        let code = unsafe { sys::myb_cip_create(pool.as_ptr(), prefix.as_ptr(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Cip { raw: created(raw)? })
    }

    /// X, Y, s (rad) at the TT instant and their rates (rad/s).
    pub fn xys(&self, tt: Time) -> (Vec3, Vec3) {
        let (mut xys, mut rates) = ([0.0; 3], [0.0; 3]);
        // SAFETY: arrays of three numbers
        unsafe { sys::myb_cip_xys(self.raw.as_ptr(), tt, xys.as_mut_ptr(), rates.as_mut_ptr()) };
        (xys, rates)
    }

    /// The same at many instants.
    pub fn xys_many(&self, tt: &[Time]) -> (Vec<Vec3>, Vec<Vec3>) {
        let (mut xys, mut rates) = (vec![[0.0; 3]; tt.len()], vec![[0.0; 3]; tt.len()]);
        if !tt.is_empty() {
            // SAFETY: the outputs hold three numbers per instant
            unsafe { sys::myb_cip_xys_many(self.raw.as_ptr(), tt.as_ptr(), tt.len(), xys.as_mut_ptr().cast(), rates.as_mut_ptr().cast()) };
        }
        (xys, rates)
    }

    /// Tabulates the series every `step` seconds of TT (0: an hour) from tt0 to tt1.
    pub fn tabulate(&mut self, tt0: Time, tt1: Time, step: f64) -> Result<()> {
        // SAFETY: the model is alive and not shared while it changes
        check(unsafe { sys::myb_cip_tabulate(self.raw.as_ptr(), tt0, tt1, step) })
    }

    /// The Earth rotation angle (rad) at the UT1 instant and its rate (rad per UT1 second).
    pub fn era(&self, ut1: Time) -> (f64, f64) {
        let mut rate = 0.0;
        // SAFETY: the model is alive
        let angle = unsafe { sys::myb_cip_era(self.raw.as_ptr(), ut1, &mut rate) };
        (angle, rate)
    }

    /// The TIO locator s' (rad) at the TT instant.
    pub fn sprime(&self, tt: Time) -> f64 {
        // SAFETY: the model is alive
        unsafe { sys::myb_cip_sprime(self.raw.as_ptr(), tt) }
    }

    /// The rotation celestial -> terrestrial at tt and ut1 with the EOP (`None`: all zero).
    pub fn terrestrial(&self, tt: Time, ut1: Time, eop: Option<&EopValues>) -> Rotation {
        let mut out = Rotation::default();
        // SAFETY: valid places for the inputs and the output
        unsafe { sys::myb_cip_terrestrial(self.raw.as_ptr(), tt, ut1, eop.map_or(ptr::null(), |e| e as *const EopValues), &mut out) };
        out
    }

    /// The same at many instants, with the EOP of each (`None`: all zero).
    pub fn terrestrial_many(&self, tt: &[Time], ut1: &[Time], eop: Option<&[EopValues]>) -> Result<Vec<Rotation>> {
        times_check(tt, ut1)?;
        if eop.is_some_and(|e| e.len() != tt.len()) {
            return domain();
        }
        let mut out = vec![Rotation::default(); tt.len()];
        if !tt.is_empty() {
            // SAFETY: the arrays hold one entry per instant
            unsafe {
                sys::myb_cip_terrestrial_many(
                    self.raw.as_ptr(),
                    tt.as_ptr(),
                    ut1.as_ptr(),
                    eop.map_or(ptr::null(), |e| e.as_ptr()),
                    tt.len(),
                    out.as_mut_ptr(),
                )
            };
        }
        Ok(out)
    }

    pub fn as_ptr(&self) -> *const sys::myb_cip {
        self.raw.as_ptr()
    }
}

impl Drop for Cip {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_cip_destroy(self.raw.as_ptr()) }
    }
}

// ---- Earth orientation parameters -----------------------------------------------------------------

/// The daily EOP series of the IERS and their interpolation (`myb_eop`).
pub struct Eop {
    raw: NonNull<sys::myb_eop>,
}

type EopLoad = unsafe extern "C" fn(*mut sys::myb_eop, *const c_char, usize, *mut sys::myb_diag) -> sys::myb_status;
type EopLoadFile = unsafe extern "C" fn(*mut sys::myb_eop, *const c_char, *mut sys::myb_diag) -> sys::myb_status;

impl Eop {
    /// An empty table.
    pub fn new() -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: a valid place for the pointer
        check(unsafe { sys::myb_eop_create(ptr::null(), &mut raw) })?;
        Ok(Eop { raw: created(raw)? })
    }

    fn load(&mut self, load: EopLoad, text: &str) -> Result<()> {
        let mut diag = sys::myb_diag::default();
        // SAFETY: the text is text.len() bytes
        let code = unsafe { load(self.raw.as_ptr(), text.as_ptr().cast(), text.len(), &mut diag) };
        check_diag(code, &diag)
    }

    fn load_file(&mut self, load: EopLoadFile, path: &Path) -> Result<()> {
        let path = c_path(path)?;
        let mut diag = sys::myb_diag::default();
        // SAFETY: the path is NUL-terminated
        let code = unsafe { load(self.raw.as_ptr(), path.as_ptr(), &mut diag) };
        check_diag(code, &diag)
    }

    /// Replaces the table by an EOP 20 C04 text.
    pub fn load_c04(&mut self, text: &str) -> Result<()> {
        self.load(sys::myb_eop_load_c04, text)
    }

    pub fn load_c04_file(&mut self, path: impl AsRef<Path>) -> Result<()> {
        self.load_file(sys::myb_eop_load_c04_file, path.as_ref())
    }

    /// Replaces the table by a finals2000A text.
    pub fn load_finals(&mut self, text: &str) -> Result<()> {
        self.load(sys::myb_eop_load_finals, text)
    }

    pub fn load_finals_file(&mut self, path: impl AsRef<Path>) -> Result<()> {
        self.load_file(sys::myb_eop_load_finals_file, path.as_ref())
    }

    /// The diurnal and semidiurnal corrections of the kernel variables under a prefix.
    pub fn set_tides(&mut self, pool: &KernelPool, prefix: &str) -> Result<()> {
        let prefix = c_string(prefix)?;
        let mut diag = sys::myb_diag::default();
        // SAFETY: the pool is alive for the call
        let code = unsafe { sys::myb_eop_set_tides(self.raw.as_ptr(), pool.as_ptr(), prefix.as_ptr(), &mut diag) };
        check_diag(code, &diag)
    }

    /// The rows of the table and its first and last MJD (UTC).
    pub fn rows(&self) -> (usize, f64, f64) {
        let (mut first, mut last) = (0.0, 0.0);
        // SAFETY: the table is alive
        let rows = unsafe { sys::myb_eop_rows(self.raw.as_ptr(), &mut first, &mut last) };
        (rows, first, last)
    }

    /// The EOP at an instant of the named scale of the time system, and the same instant on UT1.
    pub fn at(&self, system: &TimeSystem, t: Time, scale: &str) -> Result<(EopValues, Time)> {
        let scale = system.scale(scale)?;
        let (mut values, mut ut1) = (EopValues::default(), Time::default());
        // SAFETY: the table and the system are alive
        check(unsafe { sys::myb_eop_at(self.raw.as_ptr(), system.as_ptr(), t, scale, &mut values, &mut ut1) })?;
        Ok((values, ut1))
    }

    /// The same at many instants.
    pub fn at_many(&self, system: &TimeSystem, t: &[Time], scale: &str) -> Result<(Vec<EopValues>, Vec<Time>)> {
        let scale = system.scale(scale)?;
        let (mut values, mut ut1) = (vec![EopValues::default(); t.len()], vec![Time::default(); t.len()]);
        if t.is_empty() {
            return Ok((values, ut1));
        }
        let mut failed = 0usize;
        // SAFETY: the outputs hold one entry per instant
        check_at(
            unsafe {
                sys::myb_eop_at_many(self.raw.as_ptr(), system.as_ptr(), t.as_ptr(), t.len(), scale, values.as_mut_ptr(), ut1.as_mut_ptr(), &mut failed)
            },
            failed,
        )?;
        Ok((values, ut1))
    }

    pub fn as_ptr(&self) -> *const sys::myb_eop {
        self.raw.as_ptr()
    }
}

impl Drop for Eop {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_eop_destroy(self.raw.as_ptr()) }
    }
}

/// The CIO-based model with observed EOP as a provider (`myb_orientation_cip`): the argument
/// is TT on the named scale of the time system; the rotation celestial -> terrestrial, the
/// rotation time UT1 and the polar motion.
pub struct CipOrientation {
    raw: sys::myb_cip_orientation,
    _kept: (Rc<Cip>, Rc<Eop>, Rc<TimeSystem>),
}

impl CipOrientation {
    pub fn new(cip: &Rc<Cip>, eop: &Rc<Eop>, timesys: &Rc<TimeSystem>, tt: &str) -> Result<Self> {
        let tt = timesys.scale(tt)?;
        Ok(CipOrientation {
            raw: sys::myb_cip_orientation { cip: cip.as_ptr(), eop: eop.as_ptr(), timesys: timesys.as_ptr(), tt },
            _kept: (Rc::clone(cip), Rc::clone(eop), Rc::clone(timesys)),
        })
    }
}

// SAFETY: myb_orientation_cip reads the structure inside self, whose objects self keeps
unsafe impl OrientationSource for CipOrientation {
    fn as_raw(&self) -> (sys::myb_orientation_function, *mut c_void) {
        (Some(sys::myb_orientation_cip), &self.raw as *const sys::myb_cip_orientation as *mut c_void)
    }
}

// ---- Equinox-based systems ------------------------------------------------------------------------

/// The frames of the equinox-based systems.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquinoxFrame {
    /// Mean of date.
    Mod,
    /// True of date.
    Tod,
    /// True equator, mean equinox (the frame of SGP4).
    Teme,
}

/// The classical systems of the Earth: IAU 1976 precession, IAU 1980 nutation, GMST 1982
/// (`myb_equinox`), from the kernel variables under a prefix ("IAU1980" of data/frames/iau1980.tk).
pub struct Equinox {
    raw: NonNull<sys::myb_equinox>,
}

impl Equinox {
    pub fn new(pool: &KernelPool, prefix: &str) -> Result<Self> {
        let prefix = c_string(prefix)?;
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the pool is alive for the call
        let code = unsafe { sys::myb_equinox_create(pool.as_ptr(), prefix.as_ptr(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Equinox { raw: created(raw)? })
    }

    /// The precession angles zeta_A, z_A, theta_A (rad) and their rates.
    pub fn precession(&self, tt: Time) -> (Vec3, Vec3) {
        let (mut angles, mut rates) = ([0.0; 3], [0.0; 3]);
        // SAFETY: arrays of three numbers
        unsafe { sys::myb_equinox_precession(self.raw.as_ptr(), tt, angles.as_mut_ptr(), rates.as_mut_ptr()) };
        (angles, rates)
    }

    /// The mean obliquity (rad) and its rate.
    pub fn obliquity(&self, tt: Time) -> (f64, f64) {
        let mut rate = 0.0;
        // SAFETY: the model is alive
        let eps = unsafe { sys::myb_equinox_obliquity(self.raw.as_ptr(), tt, &mut rate) };
        (eps, rate)
    }

    /// The nutation (Delta psi, Delta eps) and their rates.
    pub fn nutation(&self, tt: Time) -> ([f64; 2], [f64; 2]) {
        let (mut nutation, mut rates) = ([0.0; 2], [0.0; 2]);
        // SAFETY: arrays of two numbers
        unsafe { sys::myb_equinox_nutation(self.raw.as_ptr(), tt, nutation.as_mut_ptr(), rates.as_mut_ptr()) };
        (nutation, rates)
    }

    fn frame_of(frame: EquinoxFrame) -> sys::myb_equinox_frame {
        match frame {
            EquinoxFrame::Mod => sys::MYB_EQUINOX_MOD,
            EquinoxFrame::Tod => sys::MYB_EQUINOX_TOD,
            EquinoxFrame::Teme => sys::MYB_EQUINOX_TEME,
        }
    }

    /// The rotation J2000 -> frame at tt.
    pub fn rotation(&self, frame: EquinoxFrame, tt: Time) -> Rotation {
        let mut out = Rotation::default();
        // SAFETY: a valid place for the output
        unsafe { sys::myb_equinox_rotation(self.raw.as_ptr(), Self::frame_of(frame), tt, &mut out) };
        out
    }

    /// The same at many instants.
    pub fn rotation_many(&self, frame: EquinoxFrame, tt: &[Time]) -> Vec<Rotation> {
        let mut out = vec![Rotation::default(); tt.len()];
        if !tt.is_empty() {
            // SAFETY: out holds one rotation per instant
            unsafe { sys::myb_equinox_rotation_many(self.raw.as_ptr(), Self::frame_of(frame), tt.as_ptr(), tt.len(), out.as_mut_ptr()) };
        }
        out
    }

    /// GMST 1982 (rad) at the UT1 instant and its rate.
    pub fn gmst(&self, ut1: Time) -> (f64, f64) {
        let mut rate = 0.0;
        // SAFETY: the model is alive
        let angle = unsafe { sys::myb_equinox_gmst(self.raw.as_ptr(), ut1, &mut rate) };
        (angle, rate)
    }

    /// The rotation TEME -> PEF at the UT1 instant.
    pub fn teme_pef(&self, ut1: Time) -> Rotation {
        let mut out = Rotation::default();
        // SAFETY: a valid place for the output
        unsafe { sys::myb_equinox_teme_pef(self.raw.as_ptr(), ut1, &mut out) };
        out
    }

    /// The same at many instants.
    pub fn teme_pef_many(&self, ut1: &[Time]) -> Vec<Rotation> {
        let mut out = vec![Rotation::default(); ut1.len()];
        if !ut1.is_empty() {
            // SAFETY: out holds one rotation per instant
            unsafe { sys::myb_equinox_teme_pef_many(self.raw.as_ptr(), ut1.as_ptr(), ut1.len(), out.as_mut_ptr()) };
        }
        out
    }

    pub fn as_ptr(&self) -> *const sys::myb_equinox {
        self.raw.as_ptr()
    }
}

impl Drop for Equinox {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_equinox_destroy(self.raw.as_ptr()) }
    }
}

// ---- WGCCRE ---------------------------------------------------------------------------------------

/// The rotation model of a body of the WGCCRE from a text PCK (`myb_wgccre`).
pub struct Wgccre {
    raw: NonNull<sys::myb_wgccre>,
}

impl Wgccre {
    pub fn new(pool: &KernelPool, body: i64) -> Result<Self> {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the pool is alive for the call
        let code = unsafe { sys::myb_wgccre_create(pool.as_ptr(), body, ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Wgccre { raw: created(raw)? })
    }

    /// The NAIF code of the reference frame of the model.
    pub fn frame(&self) -> i32 {
        // SAFETY: the model is alive
        unsafe { sys::myb_wgccre_frame(self.raw.as_ptr()) }
    }

    /// The Euler angles (rad) and rates (rad/s) at the TDB instant.
    pub fn euler(&self, t: Time) -> (Vec3, Vec3) {
        let (mut angles, mut rates) = ([0.0; 3], [0.0; 3]);
        // SAFETY: arrays of three numbers
        unsafe { sys::myb_wgccre_euler(self.raw.as_ptr(), t, angles.as_mut_ptr(), rates.as_mut_ptr()) };
        (angles, rates)
    }

    /// The rotation reference frame -> body-fixed frame at the TDB instant.
    pub fn rotation(&self, t: Time) -> Rotation {
        let mut out = Rotation::default();
        // SAFETY: a valid place for the output
        unsafe { sys::myb_wgccre_rotation(self.raw.as_ptr(), t, &mut out) };
        out
    }

    /// The same at many instants.
    pub fn rotation_many(&self, t: &[Time]) -> Vec<Rotation> {
        let mut out = vec![Rotation::default(); t.len()];
        if !t.is_empty() {
            // SAFETY: out holds one rotation per instant
            unsafe { sys::myb_wgccre_rotation_many(self.raw.as_ptr(), t.as_ptr(), t.len(), out.as_mut_ptr()) };
        }
        out
    }

    pub fn as_ptr(&self) -> *const sys::myb_wgccre {
        self.raw.as_ptr()
    }
}

impl Drop for Wgccre {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_wgccre_destroy(self.raw.as_ptr()) }
    }
}

/// A WGCCRE model as a provider (`myb_orientation_wgccre`); base: the rotation from the
/// inertial axes in use to the model's reference frame (`None`: they coincide).
pub struct WgccreOrientation {
    raw: sys::myb_wgccre_orientation,
    _base: Option<Box<Rotation>>,
    _model: Rc<Wgccre>,
}

impl WgccreOrientation {
    pub fn new(model: &Rc<Wgccre>, base: Option<&Rotation>) -> Self {
        let base = base.map(|b| Box::new(*b));
        WgccreOrientation {
            raw: sys::myb_wgccre_orientation { model: model.as_ptr(), base: base.as_deref().map_or(ptr::null(), |b| b as *const Rotation) },
            _base: base,
            _model: Rc::clone(model),
        }
    }
}

// SAFETY: myb_orientation_wgccre reads the structure inside self; the model and the boxed base
// it points to live with self
unsafe impl OrientationSource for WgccreOrientation {
    fn as_raw(&self) -> (sys::myb_orientation_function, *mut c_void) {
        (Some(sys::myb_orientation_wgccre), &self.raw as *const sys::myb_wgccre_orientation as *mut c_void)
    }
}

// ---- Binary PCK -----------------------------------------------------------------------------------

/// A set of binary PCK files (`myb_pck`).
pub struct Pck {
    raw: NonNull<sys::myb_pck>,
    kept: Kept,
}

impl Pck {
    pub fn new() -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: a valid place for the pointer
        check(unsafe { sys::myb_pck_create(ptr::null(), &mut raw) })?;
        Ok(Pck { raw: created(raw)?, kept: Kept::default() })
    }

    pub fn load_file(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let path = c_path(path.as_ref())?;
        let mut diag = sys::myb_diag::default();
        // SAFETY: the path is NUL-terminated
        let code = unsafe { sys::myb_pck_load_file(self.raw.as_ptr(), path.as_ptr(), &mut diag) };
        check_diag(code, &diag)
    }

    /// A PCK image in memory; the set keeps the bytes.
    pub fn load_memory(&mut self, image: Vec<u8>) -> Result<()> {
        let image = image.into_boxed_slice();
        let (data, size) = (image.as_ptr(), image.len());
        self.kept.keep(image); // the bytes stay where they are
        let mut diag = sys::myb_diag::default();
        // SAFETY: the image lives in self.kept, unchanged, as long as the set
        let code = unsafe { sys::myb_pck_load_memory(self.raw.as_ptr(), data.cast(), size, &mut diag) };
        check_diag(code, &diag)
    }

    /// The Euler angles (rad), their rates (rad/s) and the base frame of frame class `id` at t.
    pub fn euler(&self, id: i32, t: Time) -> Result<(Vec3, Vec3, i32)> {
        let (mut angles, mut rates, mut base) = ([0.0; 3], [0.0; 3], 0);
        // SAFETY: arrays of three numbers
        check(unsafe { sys::myb_pck_euler(self.raw.as_ptr(), id, t, angles.as_mut_ptr(), rates.as_mut_ptr(), &mut base) })?;
        Ok((angles, rates, base))
    }

    /// The rotation base frame -> frame class `id` at t, and the base frame.
    pub fn rotation(&self, id: i32, t: Time) -> Result<(Rotation, i32)> {
        let (mut rotation, mut base) = (Rotation::default(), 0);
        // SAFETY: valid places for the outputs
        check(unsafe { sys::myb_pck_rotation(self.raw.as_ptr(), id, t, &mut rotation, &mut base) })?;
        Ok((rotation, base))
    }

    /// The same at many instants.
    pub fn rotation_many(&self, id: i32, t: &[Time]) -> Result<(Vec<Rotation>, Vec<i32>)> {
        let (mut rotations, mut bases) = (vec![Rotation::default(); t.len()], vec![0i32; t.len()]);
        if t.is_empty() {
            return Ok((rotations, bases));
        }
        let mut failed = 0usize;
        // SAFETY: the outputs hold one entry per instant
        check_at(
            unsafe { sys::myb_pck_rotation_many(self.raw.as_ptr(), id, t.as_ptr(), t.len(), rotations.as_mut_ptr(), bases.as_mut_ptr(), &mut failed) },
            failed,
        )?;
        Ok((rotations, bases))
    }

    pub fn as_ptr(&self) -> *const sys::myb_pck {
        self.raw.as_ptr()
    }
}

impl Drop for Pck {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once; the images after it
        unsafe { sys::myb_pck_destroy(self.raw.as_ptr()) }
    }
}

/// A frame of a PCK set as a provider (`myb_orientation_pck`); base: the rotation from the
/// inertial axes in use to the base frame of the segments (`None`: they coincide).
pub struct PckOrientation {
    raw: sys::myb_pck_orientation,
    _base: Option<Box<Rotation>>,
    _pck: Rc<Pck>,
}

impl PckOrientation {
    pub fn new(pck: &Rc<Pck>, id: i32, base: Option<&Rotation>) -> Self {
        let base = base.map(|b| Box::new(*b));
        PckOrientation {
            raw: sys::myb_pck_orientation { pck: pck.as_ptr(), id, base: base.as_deref().map_or(ptr::null(), |b| b as *const Rotation) },
            _base: base,
            _pck: Rc::clone(pck),
        }
    }
}

// SAFETY: myb_orientation_pck reads the structure inside self; the set and the boxed base it
// points to live with self
unsafe impl OrientationSource for PckOrientation {
    fn as_raw(&self) -> (sys::myb_orientation_function, *mut c_void) {
        (Some(sys::myb_orientation_pck), &self.raw as *const sys::myb_pck_orientation as *mut c_void)
    }
}

// ---- The graph of frames --------------------------------------------------------------------------

/// Frames named by codes, each given in its base by a provider or a constant rotation
/// (`myb_frame_graph`). Build it, then share it through an `Rc` with what queries it.
pub struct FrameGraph {
    raw: NonNull<sys::myb_frame_graph>,
    kept: Kept,
    _timesys: Option<Rc<TimeSystem>>,
}

impl FrameGraph {
    /// One time scale: every provider gets the instant asked.
    pub fn new() -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: no options; a valid place for the pointer
        check(unsafe { sys::myb_frame_graph_create(ptr::null(), ptr::null(), &mut raw) })?;
        Ok(FrameGraph { raw: created(raw)?, kept: Kept::default(), _timesys: None })
    }

    /// The instants asked for on the named scale of the time system, which converts them for
    /// providers registered with other scales; the graph keeps the system.
    pub fn with_time(timesys: &Rc<TimeSystem>, scale: &str) -> Result<Self> {
        let scale = c_string(scale)?;
        let options = sys::myb_frame_graph_options { timesys: timesys.as_ptr(), scale: scale.as_ptr() };
        let mut raw = ptr::null_mut();
        // SAFETY: the system outlives the graph (it keeps it); the name is copied or looked up
        check(unsafe { sys::myb_frame_graph_create(&options, ptr::null(), &mut raw) })?;
        Ok(FrameGraph { raw: created(raw)?, kept: Kept::default(), _timesys: Some(Rc::clone(timesys)) })
    }

    /// Frame `id` given in `base` by a provider of the rotation base -> id; scale: that of the
    /// provider's argument (`None`: the graph's).
    pub fn add(&mut self, id: i32, base: i32, orientation: impl OrientationSource + 'static, scale: Option<&str>) -> Result<()> {
        let scale: Option<CString> = scale.map(c_string).transpose()?;
        let kept = self.kept.keep(orientation);
        // SAFETY: the provider lives in self.kept, as long as the graph
        let (function, context) = unsafe { (*kept).as_raw() };
        check(unsafe {
            sys::myb_frame_graph_add(self.raw.as_ptr(), id, base, function, context, scale.as_ref().map_or(ptr::null(), |s| s.as_ptr()))
        })
    }

    /// Frame `id` at a constant rotation base -> id.
    pub fn add_fixed(&mut self, id: i32, base: i32, rotation: &Rotation) -> Result<()> {
        // SAFETY: the rotation is copied
        check(unsafe { sys::myb_frame_graph_add_fixed(self.raw.as_ptr(), id, base, rotation) })
    }

    /// The fixed frame `id` of a frames kernel, in the frame its kernel names.
    pub fn add_tkframe(&mut self, pool: &KernelPool, id: i32) -> Result<()> {
        let mut diag = sys::myb_diag::default();
        // SAFETY: the pool is alive for the call
        let code = unsafe { sys::myb_frame_graph_add_tkframe(self.raw.as_ptr(), pool.as_ptr(), id, &mut diag) };
        check_diag(code, &diag)
    }

    /// The rotation from frame `from` to frame `to` at t.
    pub fn rotation(&self, from: i32, to: i32, t: Time) -> Result<Rotation> {
        let mut out = Rotation::default();
        // SAFETY: the graph and its providers are alive
        check(unsafe { sys::myb_frame_graph_rotation(self.raw.as_ptr(), from, to, t, &mut out) })?;
        Ok(out)
    }

    /// The same at many instants.
    pub fn rotation_many(&self, from: i32, to: i32, t: &[Time]) -> Result<Vec<Rotation>> {
        let mut out = vec![Rotation::default(); t.len()];
        if t.is_empty() {
            return Ok(out);
        }
        let mut failed = 0usize;
        // SAFETY: out holds one rotation per instant
        check_at(unsafe { sys::myb_frame_graph_rotation_many(self.raw.as_ptr(), from, to, t.as_ptr(), t.len(), out.as_mut_ptr(), &mut failed) }, failed)?;
        Ok(out)
    }

    pub fn as_ptr(&self) -> *const sys::myb_frame_graph {
        self.raw.as_ptr()
    }
}

impl Drop for FrameGraph {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once; the providers and the system after it
        unsafe { sys::myb_frame_graph_destroy(self.raw.as_ptr()) }
    }
}

/// Two frames of a graph as a provider (`myb_orientation_graph`): the rotation from -> to.
pub struct GraphOrientation {
    raw: sys::myb_graph_orientation,
    _graph: Rc<FrameGraph>,
}

impl GraphOrientation {
    pub fn new(graph: &Rc<FrameGraph>, from: i32, to: i32) -> Self {
        GraphOrientation { raw: sys::myb_graph_orientation { graph: graph.as_ptr(), from, to }, _graph: Rc::clone(graph) }
    }
}

// SAFETY: myb_orientation_graph reads the structure inside self, whose graph self keeps
unsafe impl OrientationSource for GraphOrientation {
    fn as_raw(&self) -> (sys::myb_orientation_function, *mut c_void) {
        (Some(sys::myb_orientation_graph), &self.raw as *const sys::myb_graph_orientation as *mut c_void)
    }
}

// ---- SPK states across frames ---------------------------------------------------------------------

impl Ephemeris {
    /// Target relative to center at t in `frame`, the segments turned into it by the graph
    /// (`myb_spk_state_in`).
    pub fn state_in(&self, graph: &FrameGraph, target: i32, center: i32, t: Time, frame: i32) -> Result<State> {
        let mut y = [0.0; 6];
        // SAFETY: the set and the graph are alive; y holds six numbers
        check(unsafe { sys::myb_spk_state_in(self.as_ptr(), graph.as_ptr(), target, center, t, frame, y.as_mut_ptr()) })?;
        Ok(y)
    }

    /// The same at many instants.
    pub fn state_in_many(&self, graph: &FrameGraph, target: i32, center: i32, t: &[Time], frame: i32) -> Result<Vec<State>> {
        let mut states = vec![[0.0; 6]; t.len()];
        if t.is_empty() {
            return Ok(states);
        }
        let mut failed = 0usize;
        // SAFETY: states holds six numbers per instant
        check_at(
            unsafe {
                sys::myb_spk_state_in_many(self.as_ptr(), graph.as_ptr(), target, center, t.as_ptr(), t.len(), frame, states.as_mut_ptr().cast(), &mut failed)
            },
            failed,
        )?;
        Ok(states)
    }
}
