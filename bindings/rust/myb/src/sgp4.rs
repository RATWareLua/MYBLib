// SPDX-License-Identifier: BSD-3-Clause
//! Two-line element sets and the SGP4/SDP4 theory (`ext/sgp4`, the feature `sgp4`): parsing
//! and formatting element sets, their mean-elements messages, propagation, the fit of an
//! element set to positions, and the words of SPK segments of type 10. States are in TEME, in
//! m and m/s.
use std::ffi::c_int;

use crate::{check, check_diag, domain, sys, Calendar, Omm, Result, TimeSystem, Vec3};

use sys::sgp4 as ext;

/// An element set as the catalogue gives it (`myb_tle`).
pub type Tle = ext::myb_tle;
/// The geophysical constants of the theory.
pub type Sgp4Constants = ext::myb_sgp4_constants;

/// The constants of a standard model.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Sgp4Gravity {
    Wgs72Old,
    /// The constants of the catalogue.
    #[default]
    Wgs72,
    Wgs84,
}

/// Vallado's corrections or the operational code of the Air Force.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Sgp4Mode {
    #[default]
    Improved,
    Afspc,
}

impl Sgp4Gravity {
    fn raw(self) -> ext::myb_sgp4_gravity {
        match self {
            Sgp4Gravity::Wgs72Old => ext::MYB_SGP4_WGS72_OLD,
            Sgp4Gravity::Wgs72 => ext::MYB_SGP4_WGS72,
            Sgp4Gravity::Wgs84 => ext::MYB_SGP4_WGS84,
        }
    }
}

impl Sgp4Mode {
    fn raw(self) -> ext::myb_sgp4_mode {
        match self {
            Sgp4Mode::Improved => ext::MYB_SGP4_IMPROVED,
            Sgp4Mode::Afspc => ext::MYB_SGP4_AFSPC,
        }
    }
}

/// Parses the two lines; `checksum` false accepts lines whose checksum does not match.
pub fn tle_parse(line1: &str, line2: &str, checksum: bool) -> Result<Tle> {
    let (line1, line2) = (crate::c_string(line1)?, crate::c_string(line2)?);
    let flags = if checksum { 0 } else { ext::MYB_TLE_NO_CHECKSUM };
    let mut tle = Tle::default();
    let mut diag = sys::myb_diag::default();
    // SAFETY: NUL-terminated lines; valid places for the outputs
    let code = unsafe { ext::myb_tle_parse(line1.as_ptr(), line2.as_ptr(), flags, &mut tle, &mut diag) };
    check_diag(code, &diag)?;
    Ok(tle)
}

/// The epoch of an element set as a UTC label.
pub fn tle_epoch(tle: &Tle) -> Calendar {
    let mut utc = Calendar::default();
    // SAFETY: valid places for the input and the output
    unsafe { ext::myb_tle_epoch(tle, &mut utc) };
    utc
}

/// The two lines of an element set, with their checksums.
pub fn tle_format(tle: &Tle) -> Result<(String, String)> {
    let (mut line1, mut line2) = ([0u8; 70], [0u8; 70]);
    // SAFETY: two buffers of 70 characters
    check(unsafe { ext::myb_tle_format(tle, line1.as_mut_ptr().cast(), line2.as_mut_ptr().cast()) })?;
    let text = |line: &[u8; 70]| {
        let end = line.iter().position(|&c| c == 0).unwrap_or(line.len());
        String::from_utf8_lossy(&line[..end]).into_owned()
    };
    Ok((text(&line1), text(&line2)))
}

/// The element set of a mean-elements message.
pub fn tle_from_omm(omm: &Omm) -> Result<Tle> {
    let mut tle = Tle::default();
    // SAFETY: valid places for the input and the output
    check(unsafe { ext::myb_tle_from_omm(omm, &mut tle) })?;
    Ok(tle)
}

/// The mean-elements message of an element set.
pub fn omm_from_tle(tle: &Tle) -> Result<Omm> {
    let mut omm = Omm::default();
    // SAFETY: valid places for the input and the output
    check(unsafe { ext::myb_omm_from_tle(tle, &mut omm) })?;
    Ok(omm)
}

/// A first guess for a fit: the osculating elements of a TEME state written as mean ones
/// into the element set (its epoch, B* and identification are kept).
pub fn tle_from_state(r: &Vec3, v: &Vec3, gravity: Sgp4Gravity, tle: &mut Tle) -> Result<()> {
    // SAFETY: arrays of three numbers
    check(unsafe { ext::myb_tle_from_state(r.as_ptr(), v.as_ptr(), gravity.raw(), tle) })
}

/// How an element set is fitted; zeros stand for the defaults.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct TleFitOptions {
    /// B* is fitted too.
    pub bstar: bool,
    pub max_iterations: u32,
    pub tolerance: f64,
    pub gravity: Sgp4Gravity,
    pub mode: Sgp4Mode,
}

/// The root mean square of the residuals (m) and the steps of a fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TleFit {
    pub rms: f64,
    pub iterations: i32,
}

/// Fits the elements of `tle` (the first guess) to TEME positions at minutes from its epoch.
pub fn tle_fit(minutes: &[f64], positions: &[Vec3], options: &TleFitOptions, tle: &mut Tle) -> Result<TleFit> {
    if minutes.len() != positions.len() {
        return domain();
    }
    let made = ext::myb_tle_fit_options {
        bstar: c_int::from(options.bstar),
        max_iterations: c_int::try_from(options.max_iterations).or_else(|_| domain())?,
        tolerance: options.tolerance,
        gravity: options.gravity.raw(),
        mode: options.mode.raw(),
    };
    let mut fit = TleFit { rms: 0.0, iterations: 0 };
    // SAFETY: one time and three numbers per position
    check(unsafe {
        ext::myb_tle_fit(minutes.as_ptr(), positions.as_ptr().cast(), minutes.len(), &made, tle, &mut fit.rms, &mut fit.iterations)
    })?;
    Ok(fit)
}

/// The constants of a standard model.
pub fn sgp4_constants(gravity: Sgp4Gravity) -> Sgp4Constants {
    let mut constants = Sgp4Constants::default();
    // SAFETY: a valid place for the output
    unsafe { ext::myb_sgp4_constants_of(gravity.raw(), &mut constants) };
    constants
}

/// The code of the theory for a propagation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sgp4Code {
    Ok,
    /// The radius below the Earth's: the state is still given.
    Decayed,
    /// Another code of the theory (the call then fails with Domain).
    Other(u32),
}

/// The initialised theory of an element set (`myb_sgp4`).
#[derive(Debug, Clone, Copy)]
pub struct Sgp4(ext::myb_sgp4);

impl Sgp4 {
    pub fn new(tle: &Tle, gravity: Sgp4Gravity, mode: Sgp4Mode) -> Result<Self> {
        // SAFETY: the record is plain numbers, written whole by the library
        let mut sgp4: ext::myb_sgp4 = unsafe { std::mem::zeroed() };
        // SAFETY: valid places for the input and the output
        check(unsafe { ext::myb_sgp4_init(tle, gravity.raw(), mode.raw(), &mut sgp4) })?;
        Ok(Sgp4(sgp4))
    }

    /// From elements in the units of the theory (epoch in days since 1949 December 31 0h UTC,
    /// angles in rad, the Kozai mean motion in rad/min).
    #[allow(clippy::too_many_arguments)]
    pub fn from_elements(
        constants: &Sgp4Constants,
        mode: Sgp4Mode,
        epoch: f64,
        bstar: f64,
        eccentricity: f64,
        argument_of_perigee: f64,
        inclination: f64,
        mean_anomaly: f64,
        mean_motion: f64,
        node: f64,
    ) -> Result<Self> {
        // SAFETY: as in new
        let mut sgp4: ext::myb_sgp4 = unsafe { std::mem::zeroed() };
        // SAFETY: valid places for the input and the output
        check(unsafe {
            ext::myb_sgp4_init_elements(
                constants,
                mode.raw(),
                epoch,
                bstar,
                eccentricity,
                argument_of_perigee,
                inclination,
                mean_anomaly,
                mean_motion,
                node,
                &mut sgp4,
            )
        })?;
        Ok(Sgp4(sgp4))
    }

    fn run(
        &self,
        minutes: f64,
        call: unsafe extern "C" fn(*const ext::myb_sgp4, f64, *mut f64, *mut f64, *mut ext::myb_sgp4_error) -> sys::myb_status,
    ) -> Result<(Vec3, Vec3, Sgp4Code)> {
        let (mut r, mut v, mut code) = ([0.0; 3], [0.0; 3], 0);
        // SAFETY: the record is initialised; arrays of three numbers
        check(unsafe { call(&self.0, minutes, r.as_mut_ptr(), v.as_mut_ptr(), &mut code) })?;
        let code = match code {
            ext::MYB_SGP4_OK => Sgp4Code::Ok,
            ext::MYB_SGP4_DECAYED => Sgp4Code::Decayed,
            other => Sgp4Code::Other(other),
        };
        Ok((r, v, code))
    }

    /// The TEME position (m) and velocity (m/s) at minutes from the epoch.
    pub fn propagate(&self, minutes: f64) -> Result<(Vec3, Vec3, Sgp4Code)> {
        self.run(minutes, ext::myb_sgp4_propagate)
    }

    /// The same in km and km/s, as the verification program prints.
    pub fn propagate_km(&self, minutes: f64) -> Result<(Vec3, Vec3, Sgp4Code)> {
        self.run(minutes, ext::myb_sgp4_propagate_km)
    }

    pub fn as_raw(&self) -> &ext::myb_sgp4 {
        &self.0
    }
}

/// The eight constants of an SPK segment of type 10 for a model.
pub fn sgp4_spk10_constants(model: &Sgp4Constants) -> [f64; 8] {
    let mut constants = [0.0; 8];
    // SAFETY: a valid place for the input; constants holds eight numbers
    unsafe { ext::myb_sgp4_spk10_constants(model, constants.as_mut_ptr()) };
    constants
}

/// The packet of an element set (14 numbers) for an SPK segment of type 10, the nutation
/// left zero (no equinox model).
pub fn sgp4_spk10_packet(tle: &Tle, system: &TimeSystem) -> Result<[f64; 14]> {
    let mut packet = [0.0; 14];
    // SAFETY: the system is alive; the packet holds 14 numbers; no equinox
    check(unsafe { ext::myb_sgp4_spk10_packet(tle, system.as_ptr(), std::ptr::null(), packet.as_mut_ptr()) })?;
    Ok(packet)
}

/// The words of a segment of type 10 from its constants and packets (in the order of their
/// epochs), for [`crate::SpkWriter::add_raw`].
pub fn sgp4_spk10_words(constants: &[f64; 8], packets: &[[f64; 14]]) -> Result<Vec<f64>> {
    // SAFETY: a function of numbers
    let size = unsafe { ext::myb_sgp4_spk10_size(packets.len()) };
    if size == 0 {
        return domain();
    }
    let mut words = vec![0.0; size];
    // SAFETY: packets holds 14 numbers per set; words the size the library asks for
    check(unsafe { ext::myb_sgp4_spk10_words(constants.as_ptr(), packets.as_ptr().cast(), packets.len(), words.as_mut_ptr(), size) })?;
    Ok(words)
}

/// The number of words of a segment of `count` element sets (0 if none or too many).
pub fn sgp4_spk10_size(count: usize) -> usize {
    // SAFETY: a function of numbers
    unsafe { ext::myb_sgp4_spk10_size(count) }
}

/// The Greenwich mean sidereal time of IAU 1982 at the UT1 Julian date jd1 + jd2, as the
/// theory computes it.
pub fn sgp4_gmst(jd1: f64, jd2: f64) -> f64 {
    // SAFETY: a function of numbers
    unsafe { ext::myb_sgp4_gmst(jd1, jd2) }
}
