// SPDX-License-Identifier: BSD-3-Clause
//! The array forms of the library (`*_many` of `myb/time.h`, `myb/kin.h`, `myb/frames.h`,
//! `myb/kepler.h`) as functions over slices: one crossing into C for a whole series, the
//! same bits as the single calls in a loop. A form that can fail stops at the first element
//! that fails: the error carries its status, and its index in the message ("element k").
use crate::{
    check, domain, pv_of, state_of, sys, Calendar, Elements, Ephemeris, Equinoctial, Error, Orientation, OrientationSource, Result,
    Rotation, State, Time, TimeSystem, Vec3,
};

pub(crate) fn check_at(code: sys::myb_status, failed: usize) -> Result<()> {
    match check(code) {
        Ok(()) => Ok(()),
        Err(error) => Err(Error { status: error.status, message: format!("element {failed}") }),
    }
}

fn pvs(states: &[State]) -> Vec<sys::myb_pv> {
    states.iter().map(pv_of).collect()
}

// ---- Time -----------------------------------------------------------------------------------------

/// Instants from seconds past J2000 of their scale.
pub fn time_from_seconds_many(seconds: &[f64]) -> Vec<Time> {
    let mut t = vec![Time::default(); seconds.len()];
    if !t.is_empty() {
        // SAFETY: both arrays hold seconds.len() entries
        unsafe { sys::myb_time_from_seconds_many(seconds.as_ptr(), seconds.len(), t.as_mut_ptr()) };
    }
    t
}

/// The seconds past J2000 of instants.
pub fn time_seconds_many(t: &[Time]) -> Vec<f64> {
    let mut seconds = vec![0.0; t.len()];
    if !t.is_empty() {
        // SAFETY: both arrays hold t.len() entries
        unsafe { sys::myb_time_seconds_many(t.as_ptr(), t.len(), seconds.as_mut_ptr()) };
    }
    seconds
}

/// The instant of a calendar label on a continuous scale (`myb_calendar_to_time`).
pub fn calendar_to_time(calendar: &Calendar) -> Result<Time> {
    let mut t = Time::default();
    // SAFETY: valid places for the input and the output
    check(unsafe { sys::myb_calendar_to_time(calendar, &mut t) })?;
    Ok(t)
}

/// The calendar label of an instant (`myb_time_to_calendar`).
pub fn time_to_calendar(t: Time) -> Calendar {
    let mut calendar = Calendar::default();
    // SAFETY: a valid place for the output
    unsafe { sys::myb_time_to_calendar(t, &mut calendar) };
    calendar
}

pub fn calendar_to_time_many(calendars: &[Calendar]) -> Result<Vec<Time>> {
    let mut t = vec![Time::default(); calendars.len()];
    if t.is_empty() {
        return Ok(t);
    }
    let mut failed = 0usize;
    // SAFETY: both arrays hold calendars.len() entries
    check_at(unsafe { sys::myb_calendar_to_time_many(calendars.as_ptr(), calendars.len(), t.as_mut_ptr(), &mut failed) }, failed)?;
    Ok(t)
}

pub fn time_to_calendar_many(t: &[Time]) -> Vec<Calendar> {
    let mut calendars = vec![Calendar::default(); t.len()];
    if !t.is_empty() {
        // SAFETY: both arrays hold t.len() entries
        unsafe { sys::myb_time_to_calendar_many(t.as_ptr(), t.len(), calendars.as_mut_ptr()) };
    }
    calendars
}

impl TimeSystem {
    /// The instants on another scale.
    pub fn convert_many(&self, t: &[Time], from: &str, to: &str) -> Result<Vec<Time>> {
        let (from, to) = (self.scale(from)?, self.scale(to)?);
        let mut out = vec![Time::default(); t.len()];
        if out.is_empty() {
            return Ok(out);
        }
        let mut failed = 0usize;
        // SAFETY: both arrays hold t.len() entries; the scales are the system's
        check_at(unsafe { sys::myb_timesys_convert_many(self.as_ptr(), t.as_ptr(), t.len(), from, to, out.as_mut_ptr(), &mut failed) }, failed)?;
        Ok(out)
    }

    pub fn from_utc_many(&self, labels: &[Calendar], to: &str) -> Result<Vec<Time>> {
        let to = self.scale(to)?;
        let mut out = vec![Time::default(); labels.len()];
        if out.is_empty() {
            return Ok(out);
        }
        let mut failed = 0usize;
        // SAFETY: both arrays hold labels.len() entries
        check_at(
            unsafe { sys::myb_timesys_from_utc_many(self.as_ptr(), labels.as_ptr(), labels.len(), to, out.as_mut_ptr(), &mut failed) },
            failed,
        )?;
        Ok(out)
    }

    pub fn to_utc_many(&self, t: &[Time], from: &str) -> Result<Vec<Calendar>> {
        let from = self.scale(from)?;
        let mut out = vec![Calendar::default(); t.len()];
        if out.is_empty() {
            return Ok(out);
        }
        let mut failed = 0usize;
        // SAFETY: both arrays hold t.len() entries
        check_at(unsafe { sys::myb_timesys_to_utc_many(self.as_ptr(), t.as_ptr(), t.len(), from, out.as_mut_ptr(), &mut failed) }, failed)?;
        Ok(out)
    }
}

// ---- Ephemerides ----------------------------------------------------------------------------------

impl Ephemeris {
    /// Target relative to center at the instants, with the frame of each state.
    pub fn state_many(&self, target: i32, center: i32, t: &[Time]) -> Result<(Vec<State>, Vec<i32>)> {
        let mut states = vec![[0.0; 6]; t.len()];
        let mut frames = vec![0i32; t.len()];
        if t.is_empty() {
            return Ok((states, frames));
        }
        let mut failed = 0usize;
        // SAFETY: states holds six numbers per instant, frames one
        check_at(
            unsafe {
                sys::myb_spk_state_many(
                    self.as_ptr(),
                    target,
                    center,
                    t.as_ptr(),
                    t.len(),
                    states.as_mut_ptr().cast(),
                    frames.as_mut_ptr(),
                    &mut failed,
                )
            },
            failed,
        )?;
        Ok((states, frames))
    }
}

// ---- Frames ---------------------------------------------------------------------------------------

/// Each rotation applied to the position of the same index.
pub fn rotation_apply_position_many(rotations: &[Rotation], positions: &[Vec3]) -> Result<Vec<Vec3>> {
    if rotations.len() != positions.len() {
        return domain();
    }
    let mut out = vec![[0.0; 3]; positions.len()];
    if !out.is_empty() {
        // SAFETY: the arrays hold one rotation and three numbers per element
        unsafe { sys::myb_rotation_apply_position_many(rotations.as_ptr(), positions.as_ptr().cast(), positions.len(), out.as_mut_ptr().cast()) };
    }
    Ok(out)
}

/// Each rotation applied to the state of the same index (the turning of the frame included).
pub fn rotation_apply_state_many(rotations: &[Rotation], states: &[State]) -> Result<Vec<State>> {
    if rotations.len() != states.len() {
        return domain();
    }
    let mut out = vec![[0.0; 6]; states.len()];
    if !out.is_empty() {
        // SAFETY: the arrays hold one rotation and six numbers per element
        unsafe { sys::myb_rotation_apply_state_many(rotations.as_ptr(), states.as_ptr().cast(), states.len(), out.as_mut_ptr().cast()) };
    }
    Ok(out)
}

/// What an orientation provider gives at the instants.
pub fn orientation_many(source: &(impl OrientationSource + ?Sized), t: &[Time]) -> Result<Vec<Orientation>> {
    let mut out = vec![Orientation::default(); t.len()];
    if out.is_empty() {
        return Ok(out);
    }
    let (function, context) = source.as_raw();
    let mut failed = 0usize;
    // SAFETY: the provider lives for the call; out holds one entry per instant
    check_at(unsafe { sys::myb_orientation_many(function, context, t.as_ptr(), t.len(), out.as_mut_ptr(), &mut failed) }, failed)?;
    Ok(out)
}

// ---- Two bodies -----------------------------------------------------------------------------------

pub fn elements_from_state_many(mu: f64, states: &[State]) -> Result<Vec<Elements>> {
    let pv = pvs(states);
    let mut out = vec![Elements::default(); pv.len()];
    if out.is_empty() {
        return Ok(out);
    }
    let mut failed = 0usize;
    // SAFETY: both arrays hold one entry per state
    check_at(unsafe { sys::myb_elements_from_pv_many(mu, pv.as_ptr(), pv.len(), out.as_mut_ptr(), &mut failed) }, failed)?;
    Ok(out)
}

pub fn state_from_elements_many(mu: f64, elements: &[Elements]) -> Result<Vec<State>> {
    let mut pv = vec![sys::myb_pv::default(); elements.len()];
    if pv.is_empty() {
        return Ok(Vec::new());
    }
    let mut failed = 0usize;
    // SAFETY: both arrays hold one entry per set of elements
    check_at(unsafe { sys::myb_pv_from_elements_many(mu, elements.as_ptr(), elements.len(), pv.as_mut_ptr(), &mut failed) }, failed)?;
    Ok(pv.iter().map(state_of).collect())
}

pub fn equinoctial_from_state_many(mu: f64, states: &[State]) -> Result<Vec<Equinoctial>> {
    let pv = pvs(states);
    let mut out = vec![Equinoctial::default(); pv.len()];
    if out.is_empty() {
        return Ok(out);
    }
    let mut failed = 0usize;
    // SAFETY: both arrays hold one entry per state
    check_at(unsafe { sys::myb_equinoctial_from_pv_many(mu, pv.as_ptr(), pv.len(), out.as_mut_ptr(), &mut failed) }, failed)?;
    Ok(out)
}

pub fn state_from_equinoctial_many(mu: f64, elements: &[Equinoctial]) -> Result<Vec<State>> {
    let mut pv = vec![sys::myb_pv::default(); elements.len()];
    if pv.is_empty() {
        return Ok(Vec::new());
    }
    let mut failed = 0usize;
    // SAFETY: both arrays hold one entry per set of elements
    check_at(unsafe { sys::myb_pv_from_equinoctial_many(mu, elements.as_ptr(), elements.len(), pv.as_mut_ptr(), &mut failed) }, failed)?;
    Ok(pv.iter().map(state_of).collect())
}

/// The states of the two-body orbit through `state` after each dt.
pub fn kepler_propagate_many(mu: f64, state: &State, dt: &[f64]) -> Result<Vec<State>> {
    let mut pv = vec![sys::myb_pv::default(); dt.len()];
    if pv.is_empty() {
        return Ok(Vec::new());
    }
    let mut failed = 0usize;
    // SAFETY: pv holds one entry per offset; pv0 is not an output
    check_at(unsafe { sys::myb_kepler_propagate_many(mu, &pv_of(state), dt.as_ptr(), dt.len(), pv.as_mut_ptr(), &mut failed) }, failed)?;
    Ok(pv.iter().map(state_of).collect())
}

/// The same with the transition matrix of each (rows).
pub fn kepler_stm_many(mu: f64, state: &State, dt: &[f64]) -> Result<Vec<(State, [f64; 36])>> {
    let mut pv = vec![sys::myb_pv::default(); dt.len()];
    let mut stm = vec![[0.0; 36]; dt.len()];
    if pv.is_empty() {
        return Ok(Vec::new());
    }
    let mut failed = 0usize;
    // SAFETY: pv holds one entry per offset, stm 36 numbers per offset
    check_at(
        unsafe { sys::myb_kepler_stm_many(mu, &pv_of(state), dt.as_ptr(), dt.len(), pv.as_mut_ptr(), stm.as_mut_ptr().cast(), &mut failed) },
        failed,
    )?;
    Ok(pv.iter().map(state_of).zip(stm).collect())
}

/// The state and transition matrix (rows) of a two-body orbit after dt (`myb_kepler_stm`).
pub fn kepler_stm(mu: f64, state: &State, dt: f64) -> Result<(State, [f64; 36])> {
    let mut pv = sys::myb_pv::default();
    let mut stm = [0.0; 36];
    // SAFETY: valid places for the outputs
    check(unsafe { sys::myb_kepler_stm(mu, &pv_of(state), dt, &mut pv, stm.as_mut_ptr()) })?;
    Ok((state_of(&pv), stm))
}
