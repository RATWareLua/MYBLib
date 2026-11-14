// SPDX-License-Identifier: BSD-3-Clause
//! Observables (`myb/observe.h`): the light time between an observer and a target with the
//! relativistic delay of the bodies on the way, round trips, the Shapiro delay, aberration,
//! light deflection, right ascension and declination, the dispersive delay of a plasma.
use std::rc::Rc;

use crate::{check, sys, Result, State, StateSource, Time, Vec3};

/// The times, the target's state, the light time, the range, the delay and the range rate of
/// a one-way signal.
pub type LightPath = sys::myb_light_path;
/// A power-law electron density about a body.
pub type PlasmaProfile = sys::myb_plasma_profile;

/// A body that delays and deflects light: its GM and a provider of its state.
#[derive(Clone)]
pub struct LightBody {
    pub gm: f64,
    pub state: Rc<dyn StateSource>,
}

/// Which way the signal goes.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum LightDirection {
    /// The observer receives at t what the target emitted earlier.
    #[default]
    Reception,
    /// The observer emits at t, the target receives later.
    Transmission,
}

/// The time it takes for one round trip, and when the target returned and the station sent.
#[derive(Debug, Clone, Copy)]
pub struct RoundTrip {
    pub total: f64,
    pub t_bounce: Time,
    pub t_transmit: Time,
}

/// The speed of light, PPN gamma and the bodies on the way (none: Newtonian light time).
#[derive(Clone, Default)]
pub struct LightModel {
    pub c: f64,
    pub gamma: f64,
    pub bodies: Vec<LightBody>,
}

impl LightModel {
    fn with_raw<T>(&self, call: impl FnOnce(&sys::myb_light_model) -> Result<T>) -> Result<T> {
        let bodies: Vec<sys::myb_light_body> = self
            .bodies
            .iter()
            .map(|body| {
                let (state, context) = body.state.as_raw();
                sys::myb_light_body { gm: body.gm, state, context }
            })
            .collect();
        let model = sys::myb_light_model {
            c: self.c,
            gamma: self.gamma,
            bodies: if bodies.is_empty() { std::ptr::null() } else { bodies.as_ptr() },
            count: bodies.len(),
        };
        call(&model)
    }

    /// The one-way light time between the observer state at t and a target
    /// (`myb_light_time`).
    pub fn light_time(&self, target: &(impl StateSource + ?Sized), t: Time, observer: &State, direction: LightDirection) -> Result<LightPath> {
        let (function, context) = target.as_raw();
        let direction = match direction {
            LightDirection::Reception => sys::MYB_LIGHT_RECEPTION,
            LightDirection::Transmission => sys::MYB_LIGHT_TRANSMISSION,
        };
        self.with_raw(|model| {
            let mut path = LightPath::default();
            // SAFETY: the providers live in self and in `target` for the call
            check(unsafe { sys::myb_light_time(model, function, context, t, observer.as_ptr(), direction, &mut path) })?;
            Ok(path)
        })
    }

    /// A round trip received by the station at t (`myb_round_trip`).
    pub fn round_trip(&self, station: &(impl StateSource + ?Sized), target: &(impl StateSource + ?Sized), t: Time) -> Result<RoundTrip> {
        let (sf, sc) = station.as_raw();
        let (tf, tc) = target.as_raw();
        self.with_raw(|model| {
            let mut out = RoundTrip { total: 0.0, t_bounce: Time::default(), t_transmit: Time::default() };
            // SAFETY: the providers live for the call; valid places for the outputs
            check(unsafe { sys::myb_round_trip(model, sf, sc, tf, tc, t, &mut out.total, &mut out.t_bounce, &mut out.t_transmit) })?;
            Ok(out)
        })
    }
}

/// The relativistic delay (s) of a ray from x1 to x2 by a body gm at xb.
pub fn shapiro_delay(gm: f64, c: f64, gamma: f64, x1: &Vec3, x2: &Vec3, xb: &Vec3) -> f64 {
    // SAFETY: arrays of three numbers
    unsafe { sys::myb_shapiro_delay(gm, c, gamma, x1.as_ptr(), x2.as_ptr(), xb.as_ptr()) }
}

/// The direction p (unit) as an observer moving with velocity v sees it.
pub fn aberration(p: &Vec3, v: &Vec3, c: f64) -> Vec3 {
    let mut out = [0.0; 3];
    // SAFETY: arrays of three numbers
    unsafe { sys::myb_aberration(p.as_ptr(), v.as_ptr(), c, out.as_mut_ptr()) };
    out
}

/// The direction p deflected by a body (first order, PPN; not normalised).
#[allow(clippy::too_many_arguments)]
pub fn light_deflection(gm: f64, c: f64, gamma: f64, p: &Vec3, q: &Vec3, e: &Vec3, em: f64, dlim: f64) -> Vec3 {
    let mut out = [0.0; 3];
    // SAFETY: arrays of three numbers
    unsafe { sys::myb_light_deflection(gm, c, gamma, p.as_ptr(), q.as_ptr(), e.as_ptr(), em, dlim, out.as_mut_ptr()) };
    out
}

/// Right ascension in [0, 2 pi) and declination of a direction.
pub fn direction_radec(p: &Vec3) -> (f64, f64) {
    let (mut ra, mut dec) = (0.0, 0.0);
    // SAFETY: an array of three numbers; valid places for the outputs
    unsafe { sys::myb_direction_radec(p.as_ptr(), &mut ra, &mut dec) };
    (ra, dec)
}

/// The group delay (m) K N / f^2 of an electron column (NaN for invalid input).
pub fn plasma_group_delay(column: f64, frequency: f64, e: f64, m_e: f64, eps0: f64) -> f64 {
    // SAFETY: a function of numbers
    unsafe { sys::myb_plasma_group_delay(column, frequency, e, m_e, eps0) }
}

/// The electron column (m^-2) of a profile along the segment from p1 to p2.
pub fn plasma_column(profile: &PlasmaProfile, p1: &Vec3, p2: &Vec3) -> Result<f64> {
    let mut column = 0.0;
    // SAFETY: valid places for the inputs and the output
    check(unsafe { sys::myb_plasma_column(profile, p1.as_ptr(), p2.as_ptr(), &mut column) })?;
    Ok(column)
}
