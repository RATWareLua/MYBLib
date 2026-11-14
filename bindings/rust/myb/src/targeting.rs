// SPDX-License-Identifier: BSD-3-Clause
//! Targeting (`myb/transfer.h`, `myb/lambert.h`, `myb/iod.h`, `myb/relative.h`,
//! `myb/conjunction.h`, `myb/uncert.h`): the B-plane and patched conics, a shooting solver,
//! Lambert's problem, first orbits from a few measurements, the motion of a deputy about a
//! chief, closest approaches and the probability of collision, the variational equations of a
//! force model, covariances and deterministic random numbers.
use std::ffi::{c_int, c_void};
use std::ptr::{self, NonNull};
use std::rc::Rc;

use crate::{check, created, domain, out_ptr, pv_of, state_of, sys, ForceModel, Orbit, Result, State, StateSource, Vec3};

/// The B-plane of an approach (`myb_bplane`).
pub type BPlane = sys::myb_bplane;
/// One arc of Lambert's problem: the velocities at both ends, the revolutions, the branch.
pub type LambertSolution = sys::myb_lambert_solution;
/// An orbit from the directions alone: the state at the middle time, the ranges, the steps.
pub type IodSolution = sys::myb_iod_solution;
/// A closest approach: its time, the distance, the relative position and velocity.
pub type Approach = sys::myb_approach;
/// An encounter in the plane across the relative velocity.
pub type Encounter = sys::myb_encounter;
/// An approach with the covariance carried to it and the probability of collision.
pub type Conjunction = sys::myb_conjunction;

// ---- Transfers ------------------------------------------------------------------------------------

/// The B-plane of the two-body hyperbola through (r, v) about a target of parameter gm, for a
/// reference pole.
pub fn bplane(gm: f64, r: &Vec3, v: &Vec3, pole: &Vec3) -> Result<BPlane> {
    let mut out = BPlane::default();
    // SAFETY: arrays of three numbers; a valid place for the output
    check(unsafe { sys::myb_bplane_compute(gm, r.as_ptr(), v.as_ptr(), pole.as_ptr(), &mut out) })?;
    Ok(out)
}

/// The sphere of influence a (m/M)^(2/5).
pub fn sphere_of_influence(a: f64, mass_ratio: f64) -> f64 {
    // SAFETY: a function of numbers
    unsafe { sys::myb_sphere_of_influence(a, mass_ratio) }
}

/// The Hill radius a (1 - e) (m / (3 M))^(1/3).
pub fn hill_radius(a: f64, e: f64, mass_ratio: f64) -> f64 {
    // SAFETY: a function of numbers
    unsafe { sys::myb_hill_radius(a, e, mass_ratio) }
}

/// The turn angle of a flyby.
pub fn flyby_turn(gm: f64, v_infinity: f64, periapsis: f64) -> f64 {
    // SAFETY: a function of numbers
    unsafe { sys::myb_flyby_turn(gm, v_infinity, periapsis) }
}

/// The impulse from a circular orbit of the radius to a hyperbola of excess v_infinity.
pub fn departure_impulse(gm: f64, radius: f64, v_infinity: f64) -> f64 {
    // SAFETY: a function of numbers
    unsafe { sys::myb_departure_impulse(gm, radius, v_infinity) }
}

/// The doubles of work memory of [`shoot`] (allocated by it).
pub fn shoot_work_size(n: usize, m: usize) -> usize {
    // SAFETY: a function of numbers
    unsafe { sys::myb_shoot_work_size(n, m) }
}

/// What [`shoot`] reached.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Shot {
    pub iterations: i32,
    /// max |f_i|.
    pub residual: f64,
}

struct Shooting<'a, F> {
    n: usize,
    m: usize,
    function: &'a F,
}

unsafe extern "C" fn call_shoot<F: Fn(&[f64], &mut [f64], Option<&mut [f64]>) -> Result<()>>(
    context: *mut c_void,
    u: *const f64,
    f: *mut f64,
    jacobian: *mut f64,
) -> sys::myb_status {
    // SAFETY: the context is the Shooting of the call; u holds n numbers, f m, the Jacobian m n
    let shooting = unsafe { &*(context as *const Shooting<F>) };
    let u = unsafe { std::slice::from_raw_parts(u, shooting.n) };
    let f = unsafe { std::slice::from_raw_parts_mut(f, shooting.m) };
    let jacobian = (!jacobian.is_null()).then(|| unsafe { std::slice::from_raw_parts_mut(jacobian, shooting.m * shooting.n) });
    match (shooting.function)(u, f, jacobian) {
        Ok(()) => sys::MYB_OK,
        Err(error) => error.status.raw(),
    }
}

/// Newton's method for f(u) = 0 with u.len() unknowns and m equations (`myb_shoot`): the
/// closure fills f and, when asked, the m x n Jacobian (rows). u is updated in place.
/// max_iterations 0: 20.
pub fn shoot<F: Fn(&[f64], &mut [f64], Option<&mut [f64]>) -> Result<()>>(
    function: F,
    m: usize,
    u: &mut [f64],
    tolerance: f64,
    max_iterations: u32,
) -> Result<Shot> {
    let n = u.len();
    let options = sys::myb_shoot_options { tolerance, max_iterations: c_int::try_from(max_iterations).or_else(|_| domain())? };
    let mut work = vec![0.0; shoot_work_size(n, m)];
    let shooting = Shooting { n, m, function: &function };
    let mut shot = Shot { iterations: 0, residual: 0.0 };
    // SAFETY: the context lives for the call; u holds n numbers, work what the library asks for
    check(unsafe {
        sys::myb_shoot(
            Some(call_shoot::<F>),
            &shooting as *const Shooting<F> as *mut c_void,
            n,
            m,
            u.as_mut_ptr(),
            &options,
            out_ptr(&mut work),
            &mut shot.iterations,
            &mut shot.residual,
        )
    })?;
    Ok(shot)
}

/// All the arcs with 0 ..= max_revolutions revolutions from r1 to r2 in tof seconds.
pub fn lambert(mu: f64, r1: &Vec3, r2: &Vec3, tof: f64, retrograde: bool, max_revolutions: u32) -> Result<Vec<LambertSolution>> {
    let revolutions = c_int::try_from(max_revolutions).or_else(|_| domain())?;
    let mut solutions = vec![LambertSolution::default(); 1 + 2 * max_revolutions as usize];
    let mut count = 0usize;
    // SAFETY: solutions holds the capacity given
    check(unsafe {
        sys::myb_lambert(
            mu,
            r1.as_ptr(),
            r2.as_ptr(),
            tof,
            c_int::from(retrograde),
            revolutions,
            solutions.as_mut_ptr(),
            solutions.len(),
            &mut count,
        )
    })?;
    solutions.truncate(count);
    Ok(solutions)
}

// ---- First orbits ---------------------------------------------------------------------------------

/// The velocity at the second of three positions of one orbit (Gibbs) and the angle of r1 out
/// of the plane of r2 and r3.
pub fn iod_gibbs(mu: f64, r1: &Vec3, r2: &Vec3, r3: &Vec3) -> Result<(Vec3, f64)> {
    let (mut v2, mut out_of_plane) = ([0.0; 3], 0.0);
    // SAFETY: arrays of three numbers
    check(unsafe { sys::myb_iod_gibbs(mu, r1.as_ptr(), r2.as_ptr(), r3.as_ptr(), v2.as_mut_ptr(), &mut out_of_plane) })?;
    Ok((v2, out_of_plane))
}

/// The velocity at the second of three close positions (Herrick-Gibbs).
pub fn iod_herrick_gibbs(mu: f64, r1: &Vec3, t1: f64, r2: &Vec3, t2: f64, r3: &Vec3, t3: f64) -> Result<Vec3> {
    let mut v2 = [0.0; 3];
    // SAFETY: arrays of three numbers
    check(unsafe { sys::myb_iod_herrick_gibbs(mu, r1.as_ptr(), t1, r2.as_ptr(), t2, r3.as_ptr(), t3, v2.as_mut_ptr()) })?;
    Ok(v2)
}

type Sights = unsafe extern "C" fn(f64, *const f64, *const [f64; 3], *const [f64; 3], *mut sys::myb_iod_solution, *mut usize) -> sys::myb_status;

fn iod_candidates(method: Sights, mu: f64, t: &[f64; 3], observer: &[Vec3; 3], sight: &[Vec3; 3]) -> Result<Vec<IodSolution>> {
    let mut solutions = [IodSolution::default(); sys::MYB_IOD_MAX_SOLUTIONS as usize];
    let mut count = 0usize;
    // SAFETY: the arrays have the sizes of the C interface
    check(unsafe { method(mu, t.as_ptr(), observer.as_ptr(), sight.as_ptr(), solutions.as_mut_ptr(), &mut count) })?;
    Ok(solutions[..count.min(solutions.len())].to_vec())
}

/// Gauss's method: three unit sights from three observers at `t[0] < t[1] < t[2]`; the
/// candidate orbits, refined.
pub fn iod_gauss(mu: f64, t: &[f64; 3], observer: &[Vec3; 3], sight: &[Vec3; 3]) -> Result<Vec<IodSolution>> {
    iod_candidates(sys::myb_iod_gauss, mu, t, observer, sight)
}

/// Laplace's method, as [`iod_gauss`].
pub fn iod_laplace(mu: f64, t: &[f64; 3], observer: &[Vec3; 3], sight: &[Vec3; 3]) -> Result<Vec<IodSolution>> {
    iod_candidates(sys::myb_iod_laplace, mu, t, observer, sight)
}

/// The double-r iteration from guesses of the distances at the first two times.
pub fn iod_double_r(mu: f64, t: &[f64; 3], observer: &[Vec3; 3], sight: &[Vec3; 3], distance: &[f64; 2]) -> Result<IodSolution> {
    let mut solution = IodSolution::default();
    // SAFETY: the arrays have the sizes of the C interface
    check(unsafe { sys::myb_iod_double_r(mu, t.as_ptr(), observer.as_ptr(), sight.as_ptr(), distance.as_ptr(), &mut solution) })?;
    Ok(solution)
}

/// Gooding's method from guesses of the first and third ranges.
pub fn iod_gooding(
    mu: f64,
    t: &[f64; 3],
    observer: &[Vec3; 3],
    sight: &[Vec3; 3],
    range: &[f64; 2],
    revolutions: u32,
    branch: i32,
) -> Result<IodSolution> {
    let revolutions = c_int::try_from(revolutions).or_else(|_| domain())?;
    let mut solution = IodSolution::default();
    // SAFETY: the arrays have the sizes of the C interface
    check(unsafe {
        sys::myb_iod_gooding(mu, t.as_ptr(), observer.as_ptr(), sight.as_ptr(), range.as_ptr(), revolutions, branch, &mut solution)
    })?;
    Ok(solution)
}

// ---- Relative motion ------------------------------------------------------------------------------

/// The state of the deputy from the chief in the turning frame (R, T, N).
pub fn relative_from_inertial(chief: &State, deputy: &State) -> Result<[f64; 6]> {
    let mut relative = [0.0; 6];
    // SAFETY: valid places for the inputs and the output
    check(unsafe { sys::myb_relative_from_inertial(&pv_of(chief), &pv_of(deputy), relative.as_mut_ptr()) })?;
    Ok(relative)
}

/// The inertial state of the deputy from its relative state.
pub fn relative_to_inertial(chief: &State, relative: &[f64; 6]) -> Result<State> {
    let mut deputy = sys::myb_pv::default();
    // SAFETY: valid places for the inputs and the output
    check(unsafe { sys::myb_relative_to_inertial(&pv_of(chief), relative.as_ptr(), &mut deputy) })?;
    Ok(state_of(&deputy))
}

/// The linear map of the relative state over dt for a chief on any conic of mu (rows).
pub fn relative_stm(mu: f64, chief: &State, dt: f64) -> Result<[f64; 36]> {
    let mut phi = [0.0; 36];
    // SAFETY: phi holds 36 numbers
    check(unsafe { sys::myb_relative_stm(mu, &pv_of(chief), dt, phi.as_mut_ptr()) })?;
    Ok(phi)
}

/// The Clohessy-Wiltshire map for the mean motion n (rows).
pub fn cw_stm(n: f64, dt: f64) -> [f64; 36] {
    let mut phi = [0.0; 36];
    // SAFETY: phi holds 36 numbers
    unsafe { sys::myb_cw_stm(n, dt, phi.as_mut_ptr()) };
    phi
}

/// The two impulses of a linear rendezvous (turning frame at each, m/s); target_velocity
/// `None`: at rest relative to the chief.
pub fn relative_rendezvous(
    mu: f64,
    chief: &State,
    dt: f64,
    relative: &[f64; 6],
    target: &Vec3,
    target_velocity: Option<&Vec3>,
) -> Result<(Vec3, Vec3)> {
    let (mut dv0, mut dv1) = ([0.0; 3], [0.0; 3]);
    // SAFETY: valid places for the inputs and the outputs
    check(unsafe {
        sys::myb_relative_rendezvous(
            mu,
            &pv_of(chief),
            dt,
            relative.as_ptr(),
            target.as_ptr(),
            target_velocity.map_or(ptr::null(), |v| v.as_ptr()),
            dv0.as_mut_ptr(),
            dv1.as_mut_ptr(),
        )
    })?;
    Ok((dv0, dv1))
}

/// The relative orbital elements (da, dlambda, dex, dey, dix, diy) of the deputy.
pub fn roe_from_states(mu: f64, chief: &State, deputy: &State) -> Result<[f64; 6]> {
    let mut roe = [0.0; 6];
    // SAFETY: valid places for the inputs and the output
    check(unsafe { sys::myb_roe_from_states(mu, &pv_of(chief), &pv_of(deputy), roe.as_mut_ptr()) })?;
    Ok(roe)
}

/// The relative position (R, T, N) the elements give at the mean argument of latitude u.
pub fn roe_to_position(a: f64, u: f64, roe: &[f64; 6]) -> Vec3 {
    let mut position = [0.0; 3];
    // SAFETY: arrays of the sizes of the C interface
    unsafe { sys::myb_roe_to_position(a, u, roe.as_ptr(), position.as_mut_ptr()) };
    position
}

// ---- Conjunctions ---------------------------------------------------------------------------------

/// The closest approaches of two providers (instants `myb_time_from_seconds(t)`) for t in
/// [t0, t1], looked for in steps: at most `capacity` of them, and how many there were.
pub fn approach_find(
    first: &(impl StateSource + ?Sized),
    second: &(impl StateSource + ?Sized),
    t0: f64,
    t1: f64,
    step: f64,
    capacity: usize,
) -> Result<(Vec<Approach>, usize)> {
    let (f1, c1) = first.as_raw();
    let (f2, c2) = second.as_raw();
    let mut approaches = vec![Approach::default(); capacity];
    let mut count = 0usize;
    // SAFETY: the providers live for the call; approaches holds `capacity` entries
    check(unsafe { sys::myb_approach_find(f1, c1, f2, c2, t0, t1, step, approaches.as_mut_ptr(), capacity, &mut count) })?;
    approaches.truncate(count.min(capacity));
    Ok((approaches, count))
}

/// The closest approaches of two orbits standing at one time, taken together to t1.
pub fn approach_orbits(first: &mut Orbit, second: &mut Orbit, t1: f64, step: f64, capacity: usize) -> Result<(Vec<Approach>, usize)> {
    let mut approaches = vec![Approach::default(); capacity];
    let mut count = 0usize;
    // SAFETY: both orbits are alive and distinct; approaches holds `capacity` entries
    check(unsafe {
        sys::myb_approach_orbits(first.as_ptr(), second.as_ptr(), t1, step, approaches.as_mut_ptr(), capacity, &mut count)
    })?;
    approaches.truncate(count.min(capacity));
    Ok((approaches, count))
}

/// The encounter from the relative position and velocity and the covariance of the relative
/// position (3 x 3, m^2), for the sum of the radii.
pub fn encounter_plane(position: &Vec3, velocity: &Vec3, covariance: &[f64; 9], radius: f64) -> Result<Encounter> {
    let mut out = Encounter::default();
    // SAFETY: arrays of the sizes of the C interface
    check(unsafe { sys::myb_encounter_plane(position.as_ptr(), velocity.as_ptr(), covariance.as_ptr(), radius, &mut out) })?;
    Ok(out)
}

/// The probability of collision of a short encounter.
pub fn collision_probability(encounter: &Encounter) -> Result<f64> {
    let mut p = 0.0;
    // SAFETY: valid places for the input and the output
    check(unsafe { sys::myb_collision_probability(encounter, &mut p) })?;
    Ok(p)
}

/// The series of Chan (NaN for an invalid encounter).
pub fn collision_probability_chan(encounter: &Encounter) -> f64 {
    // SAFETY: a valid place for the input
    unsafe { sys::myb_collision_probability_chan(encounter) }
}

/// The largest probability over a scale of the covariance, and the factor where it is.
pub fn collision_probability_max(encounter: &Encounter) -> Result<(f64, f64)> {
    let (mut p, mut factor) = (0.0, 0.0);
    // SAFETY: valid places for the input and the outputs
    check(unsafe { sys::myb_collision_probability_max(encounter, &mut p, &mut factor) })?;
    Ok((p, factor))
}

/// The approaches of two variational orbits with the covariances (6 x 6) of their starting
/// states carried to each, and the probability for the sum of the radii.
#[allow(clippy::too_many_arguments)]
pub fn conjunction_find(
    first: &mut Orbit,
    first_covariance: &[f64; 36],
    second: &mut Orbit,
    second_covariance: &[f64; 36],
    radius: f64,
    t1: f64,
    step: f64,
    capacity: usize,
) -> Result<(Vec<Conjunction>, usize)> {
    let mut conjunctions = vec![Conjunction::default(); capacity];
    let mut count = 0usize;
    // SAFETY: both orbits are alive and distinct; the arrays have the sizes asked for
    check(unsafe {
        sys::myb_conjunction_find(
            first.as_ptr(),
            first_covariance.as_ptr(),
            second.as_ptr(),
            second_covariance.as_ptr(),
            radius,
            t1,
            step,
            conjunctions.as_mut_ptr(),
            capacity,
            &mut count,
        )
    })?;
    conjunctions.truncate(count.min(capacity));
    Ok((conjunctions, count))
}

/// The rate of entry into the sphere of the radius for a Gaussian relative state.
pub fn collision_rate(mean: &[f64; 6], covariance: &[f64; 36], radius: f64) -> Result<f64> {
    let mut rate = 0.0;
    // SAFETY: arrays of the sizes of the C interface
    check(unsafe { sys::myb_collision_rate(mean.as_ptr(), covariance.as_ptr(), radius, &mut rate) })?;
    Ok(rate)
}

unsafe extern "C" fn call_gaussian<F: Fn(f64) -> Result<([f64; 6], [f64; 36])>>(
    context: *mut c_void,
    t: f64,
    mean: *mut f64,
    covariance: *mut f64,
) -> sys::myb_status {
    // SAFETY: the context is the closure of the call; the outputs hold 6 and 36 numbers
    let function = unsafe { &*(context as *const F) };
    match function(t) {
        Ok((m, p)) => {
            unsafe {
                *(mean as *mut [f64; 6]) = m;
                *(covariance as *mut [f64; 36]) = p;
            }
            sys::MYB_OK
        }
        Err(error) => error.status.raw(),
    }
}

/// The expected number of entries over [t0, t1] of a relative state the closure gives as a
/// Gaussian (mean, covariance) at each time; the probability and the error estimate.
pub fn collision_probability_long<F: Fn(f64) -> Result<([f64; 6], [f64; 36])>>(
    gaussian: F,
    t0: f64,
    t1: f64,
    step: f64,
    radius: f64,
) -> Result<(f64, f64)> {
    let (mut p, mut error) = (0.0, 0.0);
    // SAFETY: the closure lives for the call
    check(unsafe {
        sys::myb_collision_probability_long(
            Some(call_gaussian::<F>),
            &gaussian as *const F as *mut c_void,
            t0,
            t1,
            step,
            radius,
            &mut p,
            &mut error,
        )
    })?;
    Ok((p, error))
}

/// The same for two variational orbits with the covariances of their starting states.
#[allow(clippy::too_many_arguments)]
pub fn conjunction_long(
    first: &mut Orbit,
    first_covariance: &[f64; 36],
    second: &mut Orbit,
    second_covariance: &[f64; 36],
    radius: f64,
    t1: f64,
    step: f64,
) -> Result<(f64, f64)> {
    let (mut p, mut error) = (0.0, 0.0);
    // SAFETY: both orbits are alive and distinct; the arrays have the sizes asked for
    check(unsafe {
        sys::myb_conjunction_long(
            first.as_ptr(),
            first_covariance.as_ptr(),
            second.as_ptr(),
            second_covariance.as_ptr(),
            radius,
            t1,
            step,
            &mut p,
            &mut error,
        )
    })?;
    Ok((p, error))
}

// ---- Variational equations ------------------------------------------------------------------------

/// How the variational equations are integrated; zeros stand for the defaults.
pub type VariationalOptions = sys::myb_variational_options;

/// The state with its transition matrix (rows) and sensitivities (6 x np, rows).
#[derive(Debug, Clone, PartialEq)]
pub struct VariationalState {
    pub y: State,
    pub phi: [f64; 36],
    pub sensitivity: Vec<f64>,
}

/// The state, the transition matrix and the sensitivities to the parameters of a model
/// (`myb_variational`).
pub struct Variational {
    raw: NonNull<sys::myb_variational>,
    _model: Rc<ForceModel>,
    parameters: usize,
}

impl Variational {
    pub fn new(model: &Rc<ForceModel>, options: &VariationalOptions) -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: the model outlives the object (it keeps it)
        check(unsafe { sys::myb_variational_create(model.as_ptr(), options, ptr::null(), &mut raw) })?;
        Ok(Variational { raw: created(raw)?, _model: Rc::clone(model), parameters: model.parameter_count() })
    }

    /// Starts at t0 from y0 with Phi = I and S = 0.
    pub fn start(&mut self, t0: f64, y0: &State) -> Result<()> {
        // SAFETY: y0 holds six numbers
        check(unsafe { sys::myb_variational_start(self.raw.as_ptr(), t0, y0.as_ptr()) })
    }

    /// One step towards t_limit.
    pub fn step(&mut self, t_limit: f64) -> Result<()> {
        // SAFETY: the object is alive
        check(unsafe { sys::myb_variational_step(self.raw.as_ptr(), t_limit) })
    }

    /// Steps up to t_end.
    pub fn integrate(&mut self, t_end: f64) -> Result<()> {
        // SAFETY: the object is alive
        check(unsafe { sys::myb_variational_integrate(self.raw.as_ptr(), t_end) })
    }

    pub fn time(&self) -> f64 {
        // SAFETY: the object is alive
        unsafe { sys::myb_variational_time(self.raw.as_ptr()) }
    }

    fn empty(&self) -> VariationalState {
        VariationalState { y: [0.0; 6], phi: [0.0; 36], sensitivity: vec![0.0; 6 * self.parameters] }
    }

    /// Where the integration stands.
    pub fn get(&self) -> VariationalState {
        let mut out = self.empty();
        // SAFETY: the outputs hold 6, 36 and 6 np numbers
        unsafe { sys::myb_variational_get(self.raw.as_ptr(), out.y.as_mut_ptr(), out.phi.as_mut_ptr(), out_ptr(&mut out.sensitivity)) };
        out
    }

    /// At any t within the last step.
    pub fn dense(&mut self, t: f64) -> Result<VariationalState> {
        let mut out = self.empty();
        // SAFETY: the outputs hold 6, 36 and 6 np numbers
        check(unsafe {
            sys::myb_variational_dense(self.raw.as_ptr(), t, out.y.as_mut_ptr(), out.phi.as_mut_ptr(), out_ptr(&mut out.sensitivity))
        })?;
        Ok(out)
    }
}

impl Drop for Variational {
    fn drop(&mut self) {
        // SAFETY: created by myb_variational_create, dropped once, before its model
        unsafe { sys::myb_variational_destroy(self.raw.as_ptr()) }
    }
}

// ---- Covariances and random numbers ---------------------------------------------------------------

/// M P M^T for n x n matrices (rows).
pub fn covariance_transform(n: usize, m: &[f64], p: &[f64]) -> Result<Vec<f64>> {
    if m.len() != n * n || p.len() != n * n {
        return domain();
    }
    let mut out = vec![0.0; n * n];
    // SAFETY: the matrices hold n^2 numbers and do not alias
    unsafe { sys::myb_covariance_transform(n, m.as_ptr(), p.as_ptr(), out_ptr(&mut out)) };
    Ok(out)
}

/// A 6 x 6 covariance in the radial, transverse, normal axes of (r, v).
pub fn covariance_rtn(r: &Vec3, v: &Vec3, p: &[f64; 36]) -> Result<[f64; 36]> {
    let mut out = [0.0; 36];
    // SAFETY: arrays of the sizes of the C interface
    check(unsafe { sys::myb_covariance_rtn(r.as_ptr(), v.as_ptr(), p.as_ptr(), out.as_mut_ptr()) })?;
    Ok(out)
}

/// xoshiro256** with Box-Muller normals: the same sequence on every platform (`myb_rng`).
#[derive(Debug, Clone, Copy)]
pub struct Rng(sys::myb_rng);

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut rng = sys::myb_rng::default();
        // SAFETY: a valid place for the generator
        unsafe { sys::myb_rng_seed(&mut rng, seed) };
        Rng(rng)
    }

    pub fn next_u64(&mut self) -> u64 {
        // SAFETY: the generator is valid
        unsafe { sys::myb_rng_next(&mut self.0) }
    }

    /// In [0, 1), 53 random bits.
    pub fn uniform(&mut self) -> f64 {
        // SAFETY: the generator is valid
        unsafe { sys::myb_rng_uniform(&mut self.0) }
    }

    /// Standard normal.
    pub fn normal(&mut self) -> f64 {
        // SAFETY: the generator is valid
        unsafe { sys::myb_rng_normal(&mut self.0) }
    }

    /// A sample of N(mean, L L^T), L lower triangular (n x n, rows).
    pub fn gaussian(&mut self, mean: &[f64], l: &[f64]) -> Result<Vec<f64>> {
        let n = mean.len();
        if l.len() != n * n {
            return domain();
        }
        let mut out = vec![0.0; n];
        // SAFETY: the arrays hold n and n^2 numbers
        unsafe { sys::myb_rng_gaussian(&mut self.0, n, mean.as_ptr(), l.as_ptr(), out_ptr(&mut out)) };
        Ok(out)
    }

    pub fn as_raw(&mut self) -> &mut sys::myb_rng {
        &mut self.0
    }
}

/// The mean and covariance (divisor count - 1) of samples of dimension n stored in a row.
pub fn sample_statistics(n: usize, samples: &[f64]) -> Result<(Vec<f64>, Vec<f64>)> {
    if n == 0 || samples.len() % n != 0 {
        return domain();
    }
    let (mut mean, mut covariance) = (vec![0.0; n], vec![0.0; n * n]);
    // SAFETY: samples holds count n numbers; the outputs n and n^2
    check(unsafe { sys::myb_sample_statistics(samples.len() / n, n, samples.as_ptr(), mean.as_mut_ptr(), covariance.as_mut_ptr()) })?;
    Ok((mean, covariance))
}
