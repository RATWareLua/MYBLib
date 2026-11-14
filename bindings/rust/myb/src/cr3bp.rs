// SPDX-License-Identifier: BSD-3-Clause
//! The circular restricted three-body problem (`myb/cr3bp.h`) in its units: the distance of
//! the primaries, their total mass and their mean motion are 1. [`Cr3bp`] holds the mass
//! ratio mu and carries the functions that take it; periodic orbits and tori are the
//! structures of the library ([`Cr3bpOrbit`], [`Cr3bpTorus`]), boxed where they are large.
use std::ffi::c_int;
use std::ptr;

use crate::{check, domain, sys, ForceModel, Result, State, Vec3};

/// A periodic orbit cut into patches, with the tangent of its family.
pub type Cr3bpOrbit = sys::myb_cr3bp_orbit;
/// An invariant circle of a torus, with the tangent of its family.
pub type Cr3bpTorus = sys::myb_cr3bp_torus;
/// The spectrum of a torus: the eigenvalues of its map, their classes, the directions of the
/// hyperbolic pair and the normal pair of an elliptic torus.
pub type Cr3bpTorusSpectrum = sys::myb_cr3bp_torus_spectrum;

/// The symmetry of a periodic orbit for single shooting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cr3bpSymmetry {
    /// y0 = (x, 0, z, 0, vy, 0); at T/2 y = 0, vx = vz = 0.
    Xz,
    /// y0 = (x, 0, 0, 0, vy, vz); at T/2 y = 0, z = 0, vx = 0.
    X,
}

/// How [`Cr3bp::correct`] corrects; zeros stand for the defaults of the library.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cr3bpCorrection {
    pub symmetry: Cr3bpSymmetry,
    /// The component of y0 held fixed.
    pub fixed: u32,
    /// The half period ends at this crossing of y = 0 (1: the first).
    pub crossing: u32,
    pub tolerance: f64,
    pub max_iterations: u32,
}

/// The family of tori a continuation follows: what it keeps.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Cr3bpTorusFamily {
    #[default]
    Jacobi,
    Period,
    Rotation,
}

/// A plane `y[axis] = value` crossed with that component growing (+1), falling (-1) or either (0).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cr3bpSection {
    pub axis: u32,
    pub value: f64,
    pub direction: i32,
}

/// The eigenvalues of a monodromy matrix and the stability index.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cr3bpStability {
    pub re: [f64; 6],
    pub im: [f64; 6],
    pub index: f64,
}

/// Where [`Cr3bp::propagate`] stopped: the state, the time from the start, and whether the
/// plane was reached.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cr3bpPropagation {
    pub y: State,
    pub elapsed: f64,
    pub crossed: bool,
}

fn int(value: u32) -> Result<c_int> {
    c_int::try_from(value).or_else(|_| domain())
}

/// The problem of the mass ratio mu = m2 / (m1 + m2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cr3bp {
    pub mu: f64,
}

impl Cr3bp {
    pub fn new(mu: f64) -> Self {
        Cr3bp { mu }
    }

    /// The acceleration in the rotating frame.
    pub fn acceleration(&self, x: &Vec3, v: &Vec3) -> Vec3 {
        let mut a = [0.0; 3];
        // SAFETY: arrays of three numbers
        unsafe { sys::myb_cr3bp_acceleration(self.mu, x.as_ptr(), v.as_ptr(), a.as_mut_ptr()) };
        a
    }

    /// d a / d x and d a / d v (rows).
    pub fn jacobian(&self, x: &Vec3) -> ([f64; 9], [f64; 9]) {
        let (mut dadx, mut dadv) = ([0.0; 9], [0.0; 9]);
        // SAFETY: the outputs hold nine numbers each
        unsafe { sys::myb_cr3bp_jacobian(self.mu, x.as_ptr(), dadx.as_mut_ptr(), dadv.as_mut_ptr()) };
        (dadx, dadv)
    }

    /// The Jacobi constant C = 2 Omega - |v|^2.
    pub fn jacobi(&self, y: &State) -> f64 {
        // SAFETY: y holds six numbers
        unsafe { sys::myb_cr3bp_jacobi(self.mu, y.as_ptr()) }
    }

    /// L1 .. L5.
    pub fn lagrange(&self) -> Result<[Vec3; 5]> {
        let mut points = [[0.0; 3]; 5];
        // SAFETY: five points of three numbers
        check(unsafe { sys::myb_cr3bp_lagrange(self.mu, points.as_mut_ptr()) })?;
        Ok(points)
    }

    /// The symmetric periodic orbit near y0 and its period.
    pub fn correct(&self, correction: &Cr3bpCorrection, y0: &State) -> Result<(State, f64)> {
        let made = sys::myb_cr3bp_correction {
            symmetry: match correction.symmetry {
                Cr3bpSymmetry::Xz => sys::MYB_CR3BP_SYMMETRY_XZ,
                Cr3bpSymmetry::X => sys::MYB_CR3BP_SYMMETRY_X,
            },
            fixed: int(correction.fixed)?,
            crossing: int(correction.crossing)?,
            tolerance: correction.tolerance,
            max_iterations: int(correction.max_iterations)?,
        };
        let (mut y, mut period) = (*y0, 0.0);
        // SAFETY: y holds six numbers
        check(unsafe { sys::myb_cr3bp_correct(self.mu, &made, y.as_mut_ptr(), &mut period) })?;
        Ok((y, period))
    }

    /// The state transition matrix over one period (rows).
    pub fn monodromy(&self, y0: &State, period: f64) -> Result<[f64; 36]> {
        let mut m = [0.0; 36];
        // SAFETY: m holds 36 numbers
        check(unsafe { sys::myb_cr3bp_monodromy(self.mu, y0.as_ptr(), period, m.as_mut_ptr()) })?;
        Ok(m)
    }

    /// A state on the unstable (`unstable`) or stable manifold at the phase (0 .. 1).
    pub fn manifold(&self, y0: &State, period: f64, phase: f64, unstable: bool, sign: i32, epsilon: f64) -> Result<State> {
        let mut y = [0.0; 6];
        // SAFETY: y holds six numbers
        check(unsafe {
            sys::myb_cr3bp_manifold(self.mu, y0.as_ptr(), period, phase, c_int::from(unstable), sign, epsilon, y.as_mut_ptr())
        })?;
        Ok(y)
    }

    /// The third-order halo orbit of Richardson about L1 (1) or L2 (2): y0 and its period.
    pub fn richardson(&self, point: i32, az: f64, branch: i32) -> Result<(State, f64)> {
        let (mut y, mut period) = ([0.0; 6], 0.0);
        // SAFETY: y holds six numbers
        check(unsafe { sys::myb_cr3bp_richardson(self.mu, point, az, branch, y.as_mut_ptr(), &mut period) })?;
        Ok((y, period))
    }

    /// The cuts of the orbit through y0 with the period, by propagation.
    pub fn orbit_from_state(&self, y0: &State, period: f64, patches: usize) -> Result<Cr3bpOrbit> {
        let mut orbit = Cr3bpOrbit::default();
        // SAFETY: valid places for the input and the output
        check(unsafe { sys::myb_cr3bp_orbit_from_state(self.mu, y0.as_ptr(), period, patches, &mut orbit) })?;
        Ok(orbit)
    }

    /// Corrects the orbit in place by multiple shooting (0: the default tolerance and budget).
    pub fn orbit_correct(&self, orbit: &mut Cr3bpOrbit, tolerance: f64, max_iterations: u32) -> Result<()> {
        // SAFETY: the orbit is valid
        check(unsafe { sys::myb_cr3bp_orbit_correct(self.mu, orbit, tolerance, int(max_iterations)?) })
    }

    /// The next orbit of the family a step along the tangent (left as it was on failure).
    pub fn orbit_continue(&self, orbit: &mut Cr3bpOrbit, step: f64, tolerance: f64, max_iterations: u32) -> Result<()> {
        // SAFETY: the orbit is valid
        check(unsafe { sys::myb_cr3bp_orbit_continue(self.mu, orbit, step, tolerance, int(max_iterations)?) })
    }

    /// The smallest torus about the periodic orbit through y0 (`myb_cr3bp_torus_from_orbit`).
    pub fn torus_from_orbit(&self, y0: &State, period: f64, points: usize, pair: u32, size: f64) -> Result<Box<Cr3bpTorus>> {
        let mut torus = Box::<Cr3bpTorus>::default();
        // SAFETY: the torus is a valid place
        check(unsafe { sys::myb_cr3bp_torus_from_orbit(self.mu, y0.as_ptr(), period, points, int(pair)?, size, &mut *torus) })?;
        Ok(torus)
    }

    /// The torus of the family a step along its tangent (0: corrects it where it is).
    pub fn torus_continue(
        &self,
        torus: &mut Cr3bpTorus,
        family: Cr3bpTorusFamily,
        step: f64,
        tolerance: f64,
        max_iterations: u32,
    ) -> Result<()> {
        let family = match family {
            Cr3bpTorusFamily::Jacobi => sys::MYB_CR3BP_TORUS_JACOBI,
            Cr3bpTorusFamily::Period => sys::MYB_CR3BP_TORUS_PERIOD,
            Cr3bpTorusFamily::Rotation => sys::MYB_CR3BP_TORUS_ROTATION,
        };
        // SAFETY: the torus is valid
        check(unsafe { sys::myb_cr3bp_torus_continue(self.mu, torus, family, step, tolerance, int(max_iterations)?) })
    }

    /// The spectrum of the torus (tolerance 0: 1e-4).
    pub fn torus_stability(&self, torus: &Cr3bpTorus, tolerance: f64) -> Result<Box<Cr3bpTorusSpectrum>> {
        let mut spectrum = Box::<Cr3bpTorusSpectrum>::default();
        // SAFETY: valid places for the input and the output
        check(unsafe { sys::myb_cr3bp_torus_stability(self.mu, torus, tolerance, &mut *spectrum) })?;
        Ok(spectrum)
    }

    /// The flow from y0 for `duration`, stopping at the first crossing of a section.
    pub fn propagate(&self, y0: &State, duration: f64, section: Option<&Cr3bpSection>) -> Result<Cr3bpPropagation> {
        let made = section
            .map(|s| -> Result<sys::myb_cr3bp_section> {
                Ok(sys::myb_cr3bp_section { axis: int(s.axis)?, value: s.value, direction: s.direction })
            })
            .transpose()?;
        let (mut y, mut elapsed, mut crossed) = ([0.0; 6], 0.0, 0);
        // SAFETY: valid places for the outputs
        check(unsafe {
            sys::myb_cr3bp_propagate(
                self.mu,
                y0.as_ptr(),
                duration,
                made.as_ref().map_or(ptr::null(), |s| s as *const sys::myb_cr3bp_section),
                y.as_mut_ptr(),
                &mut elapsed,
                &mut crossed,
            )
        })?;
        Ok(Cr3bpPropagation { y, elapsed, crossed: crossed != 0 })
    }
}

/// The eigenvalues of a monodromy matrix and the stability index of the JPL catalogue.
pub fn cr3bp_stability(monodromy: &[f64; 36]) -> Result<Cr3bpStability> {
    let mut out = Cr3bpStability { re: [0.0; 6], im: [0.0; 6], index: 0.0 };
    // SAFETY: valid places for the outputs
    check(unsafe { sys::myb_cr3bp_stability(monodromy.as_ptr(), out.re.as_mut_ptr(), out.im.as_mut_ptr(), &mut out.index) })?;
    Ok(out)
}

/// The point of the circle of a torus at an angle (rad).
pub fn cr3bp_torus_state(torus: &Cr3bpTorus, angle: f64) -> State {
    let mut y = [0.0; 6];
    // SAFETY: y holds six numbers
    unsafe { sys::myb_cr3bp_torus_state(torus, angle, y.as_mut_ptr()) };
    y
}

/// A state on the unstable or stable manifold of a torus at an angle of its circle.
pub fn cr3bp_torus_manifold(
    torus: &Cr3bpTorus,
    spectrum: &Cr3bpTorusSpectrum,
    angle: f64,
    unstable: bool,
    sign: i32,
    epsilon: f64,
) -> Result<State> {
    let mut y = [0.0; 6];
    // SAFETY: y holds six numbers
    check(unsafe { sys::myb_cr3bp_torus_manifold(torus, spectrum, angle, c_int::from(unstable), sign, epsilon, y.as_mut_ptr()) })?;
    Ok(y)
}

/// A state (m, m/s) in the pulsating rotating frame of two real bodies.
pub fn cr3bp_from_ephemeris(body1: &State, body2: &State, gm1: f64, gm2: f64, state: &State) -> Result<State> {
    let mut rotating = [0.0; 6];
    // SAFETY: arrays of six numbers
    check(unsafe { sys::myb_cr3bp_from_ephemeris(body1.as_ptr(), body2.as_ptr(), gm1, gm2, state.as_ptr(), rotating.as_mut_ptr()) })?;
    Ok(rotating)
}

/// Back from the rotating frame to a state in m, m/s.
pub fn cr3bp_to_ephemeris(body1: &State, body2: &State, gm1: f64, gm2: f64, rotating: &State) -> Result<State> {
    let mut state = [0.0; 6];
    // SAFETY: arrays of six numbers
    check(unsafe { sys::myb_cr3bp_to_ephemeris(body1.as_ptr(), body2.as_ptr(), gm1, gm2, rotating.as_ptr(), state.as_mut_ptr()) })?;
    Ok(state)
}

impl ForceModel {
    /// The CR3BP of the mass ratio mu as a term (`myb_cr3bp_term` with its Jacobian).
    pub fn add_cr3bp(&self, name: &str, mu: f64) -> Result<usize> {
        let context = self.kept.keep(mu);
        // SAFETY: mu lives in self.kept, as long as the model
        unsafe { self.add_raw(name, Some(sys::myb_cr3bp_term), Some(sys::myb_cr3bp_term_jacobian), context.cast()) }
    }
}
