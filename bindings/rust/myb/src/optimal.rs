// SPDX-License-Identifier: BSD-3-Clause
//! Low thrust (`myb/optimal.h`, `myb/thrust.h`): the indirect method, the solver of nonlinear
//! programs, the direct method of Sims and Flanagan, chains of its legs joined by flybys
//! (unpowered or powered), burn arcs with steering laws as terms of a force model, the Gauss
//! equations in equinoctial elements, Edelbaum's transfer and the Q-law.
use std::ffi::{c_int, c_void};
use std::ptr;
use std::rc::Rc;

use crate::{check, domain, out_ptr, sys, Elements, Equinoctial, ForceModel, Result, State, StateSource, Time, Vec3};

/// A solution of the indirect method: the costates at t0, the arrival, the mass left.
pub type IndirectSolution = sys::myb_indirect_solution;
/// The solution at a time: state and mass, thrust, switching function, Hamiltonian, costates.
pub type IndirectPoint = sys::myb_indirect_point;
/// Tolerances, budgets and the first weight of the solver; zeros stand for the defaults.
pub type NlpOptions = sys::myb_nlp_options;
/// What the solver reached.
pub type NlpResult = sys::myb_nlp_result;
/// The throttles of a leg that spend the least: the mass left, the mismatch, the solver's result.
pub type SimsFlanaganResult = sys::myb_sims_flanagan_result;
/// A chain solved: the masses, excess speeds, turns, periapses, mismatches, impulses of the
/// powered flybys and the results of both passes.
pub type ChainResult = sys::myb_chain_result;
/// The orbital elements a Q-law targets, their weights and the periapsis limit.
pub type Qlaw = sys::myb_qlaw;
/// Where a flight of the Q-law ended: the time, the state with the mass, whether it arrived.
pub type QlawArrival = sys::myb_qlaw_arrival;

fn model_ptr(field: &Option<Rc<ForceModel>>) -> *mut sys::myb_force_model {
    field.as_ref().map_or(ptr::null_mut(), |model| model.as_ptr())
}

fn int(value: u32) -> Result<c_int> {
    c_int::try_from(value).or_else(|_| domain())
}

// ---- The indirect method --------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum IndirectGoal {
    /// The least propellant in a given time.
    #[default]
    Fuel,
    /// The least time.
    Time,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum IndirectArrival {
    /// The state `target` (Fuel) or the target that moves (Time).
    #[default]
    Rendezvous,
    /// Anywhere on the orbit through `target`: the phase free.
    Orbit,
}

/// A transfer for the indirect method (`myb_indirect_problem`); zeros stand for the defaults.
#[derive(Clone, Default)]
pub struct IndirectProblem {
    /// Of the centre (m^3 s^-2): the field without `field`, and the units of the solver.
    pub mu: f64,
    /// Everything but the engine, the centre included; `None`: the point mass mu.
    pub field: Option<Rc<ForceModel>>,
    /// N.
    pub thrust: f64,
    /// m/s.
    pub exhaust: f64,
    pub t0: f64,
    pub start: State,
    pub mass0: f64,
    /// Fuel: the time of arrival; Time: a first guess of it.
    pub t1: f64,
    pub target: State,
    pub goal: IndirectGoal,
    pub smoothing: f64,
    pub tolerance: f64,
    pub tries: u32,
    pub arrival: IndirectArrival,
    /// Orbit, with a guess of the costates: a first guess of the phase (s).
    pub phase: f64,
}

impl IndirectProblem {
    fn raw(&self) -> Result<sys::myb_indirect_problem> {
        Ok(sys::myb_indirect_problem {
            mu: self.mu,
            field: model_ptr(&self.field),
            thrust: self.thrust,
            exhaust: self.exhaust,
            t0: self.t0,
            start: self.start,
            mass0: self.mass0,
            t1: self.t1,
            target: self.target,
            goal: match self.goal {
                IndirectGoal::Fuel => sys::MYB_INDIRECT_FUEL,
                IndirectGoal::Time => sys::MYB_INDIRECT_TIME,
            },
            smoothing: self.smoothing,
            tolerance: self.tolerance,
            tries: int(self.tries)?,
            arrival: match self.arrival {
                IndirectArrival::Rendezvous => sys::MYB_INDIRECT_RENDEZVOUS,
                IndirectArrival::Orbit => sys::MYB_INDIRECT_ORBIT,
            },
            phase: self.phase,
        })
    }

    /// Solves the conditions of Pontryagin's principle (`myb_indirect_solve`); guess: the eight
    /// costates of a solution nearby.
    pub fn solve(&self, guess: Option<&[f64; 8]>) -> Result<IndirectSolution> {
        let made = self.raw()?;
        let mut solution = IndirectSolution::default();
        // SAFETY: the field, if any, is alive in self; the guess holds eight numbers
        check(unsafe { sys::myb_indirect_solve(&made, guess.map_or(ptr::null(), |g| g.as_ptr()), ptr::null(), &mut solution) })?;
        Ok(solution)
    }

    /// The solution at increasing times between t0 and the arrival (`myb_indirect_sample`).
    pub fn sample(&self, solution: &IndirectSolution, times: &[f64]) -> Result<Vec<IndirectPoint>> {
        let made = self.raw()?;
        let mut points = vec![IndirectPoint::default(); times.len()];
        // SAFETY: points holds one entry per time
        check(unsafe { sys::myb_indirect_sample(&made, solution, times.as_ptr(), times.len(), ptr::null(), points.as_mut_ptr()) })?;
        Ok(points)
    }
}

// ---- Nonlinear programs ---------------------------------------------------------------------------

/// The shape of a nonlinear program: the constraints (equalities first) and the bounds of
/// the unknowns (`None`: none), as many as the unknowns.
#[derive(Debug, Clone, Copy, Default)]
pub struct NlpProblem<'a> {
    pub equalities: usize,
    pub inequalities: usize,
    pub lower: Option<&'a [f64]>,
    pub upper: Option<&'a [f64]>,
}

struct Program<'a, F> {
    n: usize,
    m: usize,
    function: &'a F,
}

type Derivatives<'a> = Option<(&'a mut [f64], &'a mut [f64])>;

unsafe extern "C" fn call_nlp<F: Fn(&[f64], &mut [f64], Derivatives<'_>) -> Result<f64>>(
    context: *mut c_void,
    x: *const f64,
    f: *mut f64,
    c: *mut f64,
    gradient: *mut f64,
    jacobian: *mut f64,
) -> sys::myb_status {
    // SAFETY: the context is the Program of the call; x holds n numbers, c m, the gradient n
    // and the Jacobian m n
    let program = unsafe { &*(context as *const Program<F>) };
    let (n, m) = (program.n, program.m);
    let x = unsafe { std::slice::from_raw_parts(x, n) };
    let c: &mut [f64] = if m == 0 || c.is_null() { &mut [] } else { unsafe { std::slice::from_raw_parts_mut(c, m) } };
    let derivatives = (!gradient.is_null()).then(|| {
        let gradient = unsafe { std::slice::from_raw_parts_mut(gradient, n) };
        let jacobian: &mut [f64] =
            if m == 0 || jacobian.is_null() { &mut [] } else { unsafe { std::slice::from_raw_parts_mut(jacobian, m * n) } };
        (gradient, jacobian)
    });
    match (program.function)(x, c, derivatives) {
        Ok(value) => {
            unsafe { *f = value };
            sys::MYB_OK
        }
        Err(error) => error.status.raw(),
    }
}

/// Minimises f(x) under the constraints by the method of multipliers (`myb_nlp_solve`). The
/// closure returns f at x, fills the constraints c and, when asked, the gradient of f and the
/// Jacobian of c (constraints x unknowns, rows). x: the first guess, the solution on return
/// (also on `NoConvergence`: what was reached). Returns the result and the multipliers.
pub fn nlp_solve<F: Fn(&[f64], &mut [f64], Derivatives<'_>) -> Result<f64>>(
    problem: &NlpProblem,
    function: F,
    options: Option<&NlpOptions>,
    x: &mut [f64],
) -> Result<(NlpResult, Vec<f64>)> {
    let n = x.len();
    let m = problem.equalities + problem.inequalities;
    if problem.lower.is_some_and(|b| b.len() != n) || problem.upper.is_some_and(|b| b.len() != n) {
        return domain();
    }
    let program = Program { n, m, function: &function };
    let made = sys::myb_nlp_problem {
        n,
        equalities: problem.equalities,
        inequalities: problem.inequalities,
        lower: problem.lower.map_or(ptr::null(), |b| b.as_ptr()),
        upper: problem.upper.map_or(ptr::null(), |b| b.as_ptr()),
        function: Some(call_nlp::<F>),
        context: &program as *const Program<F> as *mut c_void,
    };
    let mut multipliers = vec![0.0; m];
    let mut result = NlpResult::default();
    // SAFETY: the program lives for the call; x holds n numbers, the multipliers m
    check(unsafe {
        sys::myb_nlp_solve(
            &made,
            options.map_or(ptr::null(), |o| o as *const NlpOptions),
            x.as_mut_ptr(),
            out_ptr(&mut multipliers),
            ptr::null(),
            &mut result,
        )
    })?;
    Ok((result, multipliers))
}

// ---- Sims and Flanagan ----------------------------------------------------------------------------

/// A leg of the direct method (`myb_sims_flanagan`).
#[derive(Clone, Default)]
pub struct SimsFlanagan {
    pub mu: f64,
    pub thrust: f64,
    pub exhaust: f64,
    pub t0: f64,
    pub start: State,
    pub mass0: f64,
    pub t1: f64,
    pub target: State,
    /// 2 ..= 64.
    pub segments: usize,
    /// Everything but the engine; `None`: the conics of mu.
    pub field: Option<Rc<ForceModel>>,
}

/// Where the halves of a leg meet: the mismatch (forward minus backward, m and m/s), the mass at
/// the arrival and the Jacobian of the mismatch by the throttles (6 x 3 segments, rows).
#[derive(Debug, Clone, PartialEq)]
pub struct SimsFlanaganMismatch {
    pub mismatch: [f64; 6],
    pub mass1: f64,
    pub jacobian: Vec<f64>,
}

impl SimsFlanagan {
    fn raw(&self) -> sys::myb_sims_flanagan {
        sys::myb_sims_flanagan {
            mu: self.mu,
            thrust: self.thrust,
            exhaust: self.exhaust,
            t0: self.t0,
            start: self.start,
            mass0: self.mass0,
            t1: self.t1,
            target: self.target,
            segments: self.segments,
            field: model_ptr(&self.field),
        }
    }

    /// The mismatch for the throttles (3 per segment).
    pub fn mismatch(&self, throttles: &[f64]) -> Result<SimsFlanaganMismatch> {
        if throttles.len() != 3 * self.segments {
            return domain();
        }
        let made = self.raw();
        let mut out = SimsFlanaganMismatch { mismatch: [0.0; 6], mass1: 0.0, jacobian: vec![0.0; 18 * self.segments] };
        // SAFETY: the outputs hold 6 and 18 segments numbers
        check(unsafe {
            sys::myb_sims_flanagan_mismatch(&made, throttles.as_ptr(), out.mismatch.as_mut_ptr(), &mut out.mass1, out_ptr(&mut out.jacobian))
        })?;
        Ok(out)
    }

    /// The throttles (3 per segment; all zeros: the built-in guess) that join the states with
    /// the least propellant, in place.
    pub fn solve(&self, options: Option<&NlpOptions>, throttles: &mut [f64]) -> Result<SimsFlanaganResult> {
        if throttles.len() != 3 * self.segments {
            return domain();
        }
        let made = self.raw();
        let mut result = SimsFlanaganResult::default();
        // SAFETY: throttles holds 3 numbers per segment
        check(unsafe {
            sys::myb_sims_flanagan_solve(
                &made,
                options.map_or(ptr::null(), |o| o as *const NlpOptions),
                throttles.as_mut_ptr(),
                ptr::null(),
                &mut result,
            )
        })?;
        Ok(result)
    }
}

// ---- Chains ---------------------------------------------------------------------------------------

/// A body of a chain: a provider of its state relative to the centre (at epoch + t), its GM
/// and least flyby distance (for a body flown by), and the effective exhaust velocity of the
/// impulse of a powered flyby (0: unpowered).
#[derive(Clone)]
pub struct ChainBody {
    pub state: Rc<dyn StateSource>,
    pub gm: f64,
    pub periapsis: f64,
    pub impulse_exhaust: f64,
}

/// A leg of a chain: its segments and the shortest and longest flight (s).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChainLeg {
    pub segments: usize,
    pub flight: [f64; 2],
}

/// Legs of Sims and Flanagan from body to body joined by flybys (`myb_chain`).
#[derive(Clone)]
pub struct Chain {
    pub mu: f64,
    pub field: Option<Rc<ForceModel>>,
    /// The instant of the time 0 of the chain.
    pub epoch: Time,
    pub thrust: f64,
    pub exhaust: f64,
    pub mass0: f64,
    pub legs: Vec<ChainLeg>,
    /// The departure, those flown by, the arrival: one more than the legs.
    pub bodies: Vec<ChainBody>,
    /// The earliest and the latest departure (s).
    pub departure: [f64; 2],
    /// The largest excess speed at the departure (m/s).
    pub launch: f64,
    /// The largest at the arrival (0: a rendezvous).
    pub arrival: f64,
}

/// The unknowns and the constraints of a chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChainSize {
    pub unknowns: usize,
    pub equalities: usize,
    pub inequalities: usize,
}

/// The propellant, the constraints, the gradient and the Jacobian (rows) of a chain at x.
#[derive(Debug, Clone, PartialEq)]
pub struct ChainEvaluation {
    pub propellant: f64,
    pub constraints: Vec<f64>,
    pub gradient: Vec<f64>,
    pub jacobian: Vec<f64>,
}

/// A powered flyby: the impulse at the periapsis and the turn, with their derivatives by
/// (v_in, v_out, periapsis) — those of dv, then of the turn.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChainFlyby {
    pub dv: f64,
    pub turn: f64,
    pub partials: [f64; 6],
}

impl Chain {
    fn raw(&self) -> Result<sys::myb_chain> {
        let mut made = sys::myb_chain {
            mu: self.mu,
            field: model_ptr(&self.field),
            epoch: self.epoch,
            thrust: self.thrust,
            exhaust: self.exhaust,
            mass0: self.mass0,
            legs: self.legs.len(),
            departure: self.departure,
            launch: self.launch,
            arrival: self.arrival,
            ..Default::default()
        };
        if self.legs.len() > made.segments.len() || self.bodies.len() != self.legs.len() + 1 {
            return domain();
        }
        for (k, leg) in self.legs.iter().enumerate() {
            made.segments[k] = leg.segments;
            made.flight[k] = leg.flight;
        }
        for (k, body) in self.bodies.iter().enumerate() {
            let (state, context) = body.state.as_raw();
            made.bodies[k] = sys::myb_chain_body { state, context, gm: body.gm, periapsis: body.periapsis };
            made.impulse_exhaust[k] = body.impulse_exhaust;
        }
        Ok(made)
    }

    /// The number of unknowns and constraints (`myb_chain_size`).
    pub fn size(&self) -> Result<ChainSize> {
        let made = self.raw()?;
        let mut size = ChainSize { unknowns: 0, equalities: 0, inequalities: 0 };
        // SAFETY: the providers live in self for the call
        check(unsafe { sys::myb_chain_size(&made, &mut size.unknowns, &mut size.equalities, &mut size.inequalities) })?;
        Ok(size)
    }

    /// The first guess built in (`myb_chain_guess`).
    pub fn guess(&self) -> Result<Vec<f64>> {
        let size = self.size()?;
        let made = self.raw()?;
        let mut x = vec![0.0; size.unknowns];
        // SAFETY: x holds the unknowns
        check(unsafe { sys::myb_chain_guess(&made, x.as_mut_ptr()) })?;
        Ok(x)
    }

    /// The propellant, the constraints and their derivatives at x (`myb_chain_evaluate`).
    pub fn evaluate(&self, x: &[f64]) -> Result<ChainEvaluation> {
        let size = self.size()?;
        if x.len() != size.unknowns {
            return domain();
        }
        let m = size.equalities + size.inequalities;
        let made = self.raw()?;
        let mut out = ChainEvaluation {
            propellant: 0.0,
            constraints: vec![0.0; m],
            gradient: vec![0.0; x.len()],
            jacobian: vec![0.0; m * x.len()],
        };
        // SAFETY: the outputs hold the sizes of the chain
        check(unsafe {
            sys::myb_chain_evaluate(
                &made,
                x.as_ptr(),
                &mut out.propellant,
                out_ptr(&mut out.constraints),
                out_ptr(&mut out.gradient),
                out_ptr(&mut out.jacobian),
                ptr::null(),
            )
        })?;
        Ok(out)
    }

    /// The chain with the least propellant from x (all zeros: the built-in guess), in place
    /// (`myb_chain_solve`).
    pub fn solve(&self, options: Option<&NlpOptions>, x: &mut [f64]) -> Result<ChainResult> {
        let size = self.size()?;
        if x.len() != size.unknowns {
            return domain();
        }
        let made = self.raw()?;
        let mut result = ChainResult::default();
        // SAFETY: x holds the unknowns
        check(unsafe {
            sys::myb_chain_solve(&made, options.map_or(ptr::null(), |o| o as *const NlpOptions), x.as_mut_ptr(), ptr::null(), &mut result)
        })?;
        Ok(result)
    }
}

/// The impulse and the turn of a powered flyby for the excess speeds and the periapsis
/// (`myb_chain_flyby`).
pub fn chain_flyby(gm: f64, v_in: f64, v_out: f64, periapsis: f64) -> Result<ChainFlyby> {
    let mut out = ChainFlyby { dv: 0.0, turn: 0.0, partials: [0.0; 6] };
    // SAFETY: valid places for the outputs
    check(unsafe { sys::myb_chain_flyby(gm, v_in, v_out, periapsis, &mut out.dv, &mut out.turn, out.partials.as_mut_ptr()) })?;
    Ok(out)
}

// ---- Burn arcs ------------------------------------------------------------------------------------

/// Where the thrust of an arc points.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Steering {
    #[default]
    Velocity,
    AntiVelocity,
    /// A fixed direction in the frame of x.
    Inertial,
    /// Fixed components in the radial, transverse, normal frame.
    Rtn,
    /// Edelbaum's yaw law.
    Edelbaum,
}

/// A burn arc of constant thrust (`myb_thrust_arc`).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ThrustArc {
    pub start: f64,
    pub stop: f64,
    /// N.
    pub thrust: f64,
    /// m/s.
    pub exhaust: f64,
    pub steering: Steering,
    pub direction: Vec3,
    pub edelbaum_v0: f64,
    pub edelbaum_beta0: f64,
    pub edelbaum_sense: f64,
}

/// Burn arcs with the mass that falls along them (`myb_thrust`).
#[derive(Debug, Clone, Copy)]
pub struct Thrust(sys::myb_thrust);

impl Thrust {
    pub fn new(mass0: f64, arcs: &[ThrustArc]) -> Result<Self> {
        let made: Vec<sys::myb_thrust_arc> = arcs
            .iter()
            .map(|arc| sys::myb_thrust_arc {
                start: arc.start,
                stop: arc.stop,
                thrust: arc.thrust,
                exhaust: arc.exhaust,
                steering: match arc.steering {
                    Steering::Velocity => sys::MYB_STEER_VELOCITY,
                    Steering::AntiVelocity => sys::MYB_STEER_ANTI_VELOCITY,
                    Steering::Inertial => sys::MYB_STEER_INERTIAL,
                    Steering::Rtn => sys::MYB_STEER_RTN,
                    Steering::Edelbaum => sys::MYB_STEER_EDELBAUM,
                },
                direction: arc.direction,
                edelbaum_v0: arc.edelbaum_v0,
                edelbaum_beta0: arc.edelbaum_beta0,
                edelbaum_sense: arc.edelbaum_sense,
            })
            .collect();
        let mut thrust = sys::myb_thrust::default();
        // SAFETY: made holds arcs.len() arcs
        check(unsafe { sys::myb_thrust_init(&mut thrust, mass0, made.as_ptr(), made.len()) })?;
        Ok(Thrust(thrust))
    }

    /// The mass at t (kg).
    pub fn mass(&self, t: f64) -> f64 {
        // SAFETY: the structure is initialised
        unsafe { sys::myb_thrust_mass(&self.0, t) }
    }

    /// The same as a mass function of an object (`myb_thrust_mass_at`).
    pub fn mass_at(&self, t: f64) -> f64 {
        // SAFETY: the function only reads the structure
        unsafe { sys::myb_thrust_mass_at(&self.0 as *const sys::myb_thrust as *mut c_void, t) }
    }

    /// The acceleration of the thrust at t (`myb_thrust_acceleration`).
    pub fn acceleration(&self, t: f64, x: &Vec3, v: &Vec3) -> Result<Vec3> {
        let mut a = [0.0; 3];
        // SAFETY: the function only reads the structure; arrays of three numbers
        check(unsafe {
            sys::myb_thrust_acceleration(&self.0 as *const sys::myb_thrust as *mut c_void, t, x.as_ptr(), v.as_ptr(), a.as_mut_ptr())
        })?;
        Ok(a)
    }

    /// The acceleration in the radial, transverse, normal frame of (x, v).
    pub fn rtn(&self, t: f64, x: &Vec3, v: &Vec3) -> Result<Vec3> {
        let mut a = [0.0; 3];
        // SAFETY: arrays of three numbers
        check(unsafe { sys::myb_thrust_rtn(&self.0, t, x.as_ptr(), v.as_ptr(), a.as_mut_ptr()) })?;
        Ok(a)
    }

    pub fn as_raw(&self) -> &sys::myb_thrust {
        &self.0
    }
}

/// The rates of the modified equinoctial elements under a perturbing acceleration (R, T, N).
pub fn equinoctial_rates(mu: f64, elements: &Equinoctial, a_rtn: &Vec3) -> [f64; 6] {
    let mut rates = [0.0; 6];
    // SAFETY: valid places for the inputs and the output
    unsafe { sys::myb_equinoctial_rates(mu, elements, a_rtn.as_ptr(), rates.as_mut_ptr()) };
    rates
}

/// Edelbaum's transfer between circles of speeds v0, v1 with the plane change di: dV and the
/// initial yaw beta0.
pub fn edelbaum(v0: f64, v1: f64, di: f64) -> Result<(f64, f64)> {
    let (mut dv, mut beta0) = (0.0, 0.0);
    // SAFETY: valid places for the outputs
    check(unsafe { sys::myb_edelbaum(v0, v1, di, &mut dv, &mut beta0) })?;
    Ok((dv, beta0))
}

/// The speed and the inclination change after dv along Edelbaum's law.
pub fn edelbaum_state(v0: f64, beta0: f64, dv: f64) -> (f64, f64) {
    let (mut v, mut di) = (0.0, 0.0);
    // SAFETY: valid places for the outputs
    unsafe { sys::myb_edelbaum_state(v0, beta0, dv, &mut v, &mut di) };
    (v, di)
}

// ---- The Q-law ------------------------------------------------------------------------------------

/// The direction of a Q-law (radial, transverse, normal), Q, its rate and the effectivity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QlawDirection {
    pub direction: Vec3,
    pub q: f64,
    pub q_rate: f64,
    pub effectivity: f64,
}

/// Where the law points for the orbit `now` (`myb_qlaw_direction`).
pub fn qlaw_direction(law: &Qlaw, now: &Elements) -> Result<QlawDirection> {
    let mut out = QlawDirection { direction: [0.0; 3], q: 0.0, q_rate: 0.0, effectivity: 0.0 };
    // SAFETY: valid places for the inputs and the outputs
    check(unsafe { sys::myb_qlaw_direction(law, now, out.direction.as_mut_ptr(), &mut out.q, &mut out.q_rate, &mut out.effectivity) })?;
    Ok(out)
}

/// Whether every targeted element is within its tolerance.
pub fn qlaw_arrived(law: &Qlaw, now: &Elements, tolerance: &[f64; 5]) -> bool {
    // SAFETY: valid places for the inputs
    unsafe { sys::myb_qlaw_arrived(law, now, tolerance.as_ptr()) != 0 }
}

/// An engine of constant thrust that points as a Q-law says (`myb_qlaw_term`).
#[derive(Debug, Clone, Copy, Default)]
pub struct QlawTerm {
    pub law: Qlaw,
    pub thrust: f64,
    pub exhaust: f64,
    pub mass0: f64,
    pub t0: f64,
    pub tolerance: [f64; 5],
    /// The threshold in [0, 1); 0: the engine always burns.
    pub effectivity: f64,
}

impl QlawTerm {
    fn raw(&self, law: *const sys::myb_qlaw) -> sys::myb_qlaw_term {
        sys::myb_qlaw_term {
            law,
            thrust: self.thrust,
            exhaust: self.exhaust,
            mass0: self.mass0,
            t0: self.t0,
            tolerance: self.tolerance,
            effectivity: self.effectivity,
        }
    }

    /// The acceleration of the engine at t (`myb_term_qlaw`).
    pub fn acceleration(&self, t: f64, x: &Vec3, v: &Vec3) -> Result<Vec3> {
        let mut term = self.raw(&self.law);
        let mut a = [0.0; 3];
        // SAFETY: the term and its law live for the call; arrays of three numbers
        check(unsafe { sys::myb_term_qlaw((&mut term as *mut sys::myb_qlaw_term).cast(), t, x.as_ptr(), v.as_ptr(), a.as_mut_ptr()) })?;
        Ok(a)
    }
}

/// The same engine with the mass in the state (`myb_qlaw_flight`).
#[derive(Clone, Default)]
pub struct QlawFlight {
    pub engine: QlawTerm,
    /// Everything but the engine; `None`: the point mass of the law.
    pub field: Option<Rc<ForceModel>>,
    pub rtol: f64,
    pub atol: f64,
}

impl QlawFlight {
    /// The right-hand side for y = (x, v, m) (`myb_qlaw_flight_rhs`).
    pub fn rhs(&self, t: f64, y: &[f64; 7]) -> Result<[f64; 7]> {
        let term = self.engine.raw(&self.engine.law);
        let mut flight = sys::myb_qlaw_flight { engine: &term, field: model_ptr(&self.field), rtol: self.rtol, atol: self.atol };
        let mut dydt = [0.0; 7];
        // SAFETY: the flight, its engine and its law live for the call; arrays of seven numbers
        check(unsafe { sys::myb_qlaw_flight_rhs((&mut flight as *mut sys::myb_qlaw_flight).cast(), t, y.as_ptr(), dydt.as_mut_ptr()) })?;
        Ok(dydt)
    }

    /// Flies the law from y0 at t0 until the orbit is within the tolerances, looked at every
    /// `step` seconds, or until t_limit (`myb_qlaw_fly`).
    pub fn fly(&self, t0: f64, y0: &[f64; 7], t_limit: f64, step: f64) -> Result<QlawArrival> {
        let term = self.engine.raw(&self.engine.law);
        let flight = sys::myb_qlaw_flight { engine: &term, field: model_ptr(&self.field), rtol: self.rtol, atol: self.atol };
        let mut arrival = QlawArrival::default();
        // SAFETY: the flight, its engine and its law live for the call
        check(unsafe { sys::myb_qlaw_fly(&flight, t0, y0.as_ptr(), t_limit, step, ptr::null(), &mut arrival) })?;
        Ok(arrival)
    }
}

impl ForceModel {
    /// Burn arcs as a term (`myb_thrust_acceleration`); the model keeps a copy.
    pub fn add_thrust(&self, name: &str, thrust: &Thrust) -> Result<usize> {
        let context = self.kept.keep(thrust.0);
        // SAFETY: the copy lives in self.kept, as long as the model
        unsafe { self.add_raw(name, Some(sys::myb_thrust_acceleration), None, context.cast()) }
    }

    /// A Q-law engine as a term (`myb_term_qlaw`); the model keeps a copy of the term and its law.
    pub fn add_qlaw(&self, name: &str, term: &QlawTerm) -> Result<usize> {
        let law = self.kept.keep(term.law);
        let context = self.kept.keep(term.raw(law));
        // SAFETY: the term and its law live in self.kept, as long as the model
        unsafe { self.add_raw(name, Some(sys::myb_term_qlaw), None, context.cast()) }
    }
}
