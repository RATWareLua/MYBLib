// SPDX-License-Identifier: BSD-3-Clause
//! Systems of bodies (`myb/dyn.h`, `myb/system.h`): the Newtonian N-body problem with its graph
//! of interactions, the 1PN terms, threads and the fast multipole method; the Wisdom-Holman
//! integrators; a system of bodies with their models — figures, tides, rotations, PN pairs,
//! contacts, guided bodies, actions — integrated by IAS15 with mergers and bounces.
//!
//! Positions, velocities and accelerations are slices of 3 numbers per body (x, y, z of body 0,
//! then body 1, ...). The state of a [`System`] is 2 n numbers: its n coordinates
//! ([`System::coordinates`]) and their rates.
use std::cell::RefCell;
use std::ffi::{c_int, c_void};
use std::ptr::{self, NonNull};
use std::rc::Rc;

use crate::{check, created, domain, sys, GravityField, Ias15Timescale, Kept, OrientationSource, Result, StateSource, Vec3};

// ---- Threads --------------------------------------------------------------------------------------

/// A pool of threads (`myb/pool.h`) for the accelerations of an [`NBody`]; the results do not
/// depend on it.
pub struct Pool {
    raw: NonNull<sys::myb_pool>,
}

impl Pool {
    /// threads: the workers (the calling thread works too).
    pub fn new(threads: u32) -> Result<Self> {
        let threads = c_int::try_from(threads).or_else(|_| domain())?;
        let mut raw = ptr::null_mut();
        // SAFETY: a valid place for the pointer; the default allocator
        check(unsafe { sys::myb_pool_create(threads, ptr::null(), &mut raw) })?;
        Ok(Pool { raw: created(raw)? })
    }

    pub fn threads(&self) -> i32 {
        // SAFETY: the pool is alive
        unsafe { sys::myb_pool_threads(self.raw.as_ptr()) }
    }

    pub fn as_ptr(&self) -> *mut sys::myb_pool {
        self.raw.as_ptr()
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        // SAFETY: created by myb_pool_create, dropped once
        unsafe { sys::myb_pool_destroy(self.raw.as_ptr()) }
    }
}

// ---- The N-body problem ---------------------------------------------------------------------------

/// How a source body acts on a receiver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interaction {
    None,
    Newton,
    /// With the 1PN terms of Einstein, Infeld and Hoffmann ([`NBody::set_relativity`] first).
    Eih,
}

/// The parameters of the post-Newtonian terms: c (m/s), beta, gamma (`myb/rel.h`).
pub type Ppn = sys::myb_ppn;
/// The conserved quantities of the massive bodies, each with one factor G.
pub type NBodyInvariants = sys::myb_nbody_invariants;
/// The same at the first post-Newtonian order.
pub type NBodyPnInvariants = sys::myb_nbody_pn_invariants;

/// Point masses (`myb/dyn.h`): massive bodies first, then test bodies that attract nothing.
pub struct NBody {
    raw: NonNull<sys::myb_nbody>,
    owned: bool,
    pool: RefCell<Option<Rc<Pool>>>,
}

impl NBody {
    /// gm: the massive bodies (at least one); `test` test bodies follow them.
    pub fn new(gm: &[f64], test: usize) -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: gm holds gm.len() numbers; a valid place for the pointer
        check(unsafe { sys::myb_nbody_create(gm.len(), test, gm.as_ptr(), ptr::null(), &mut raw) })?;
        Ok(NBody { raw: created(raw)?, owned: true, pool: RefCell::new(None) })
    }

    /// Massive and test bodies.
    pub fn count(&self) -> usize {
        // SAFETY: the system is alive
        unsafe { sys::myb_nbody_count(self.raw.as_ptr()) }
    }

    pub fn massive(&self) -> usize {
        // SAFETY: the system is alive
        unsafe { sys::myb_nbody_massive(self.raw.as_ptr()) }
    }

    /// The GM of a body (0 for a test body); `None` if there is no such body.
    pub fn gm(&self, body: usize) -> Option<f64> {
        // SAFETY: the body is one of the system
        (body < self.count()).then(|| unsafe { sys::myb_nbody_gm(self.raw.as_ptr(), body) })
    }

    /// Threads for the accelerations (`None`: the calling thread only); the system keeps the pool.
    pub fn set_pool(&self, pool: Option<&Rc<Pool>>) {
        // SAFETY: the pool, if any, is kept alive by the system below
        unsafe { sys::myb_nbody_set_pool(self.raw.as_ptr(), pool.map_or(ptr::null_mut(), |p| p.as_ptr())) };
        *self.pool.borrow_mut() = pool.cloned();
    }

    /// The level of the action of `source` on `receiver` (one direction of the pair).
    pub fn set_interaction(&self, receiver: usize, source: usize, level: Interaction) -> Result<()> {
        let level = match level {
            Interaction::None => sys::MYB_INTERACTION_NONE,
            Interaction::Newton => sys::MYB_INTERACTION_NEWTON,
            Interaction::Eih => sys::MYB_INTERACTION_EIH,
        };
        // SAFETY: the system is alive; the library checks the indices
        check(unsafe { sys::myb_nbody_set_interaction(self.raw.as_ptr(), receiver, source, level) })
    }

    /// The 1PN terms for every pair not set otherwise (`None`: back to Newtonian).
    pub fn set_relativity(&self, ppn: Option<&Ppn>) -> Result<()> {
        // SAFETY: the parameters are copied
        check(unsafe { sys::myb_nbody_set_relativity(self.raw.as_ptr(), ppn.map_or(ptr::null(), |p| p as *const Ppn)) })
    }

    fn sized(&self, values: &[f64]) -> Result<()> {
        if values.len() == 3 * self.count() {
            Ok(())
        } else {
            domain()
        }
    }

    /// The accelerations at the positions x (no 1PN terms).
    pub fn acceleration(&self, t: f64, x: &[f64]) -> Result<Vec<f64>> {
        self.sized(x)?;
        let mut a = vec![0.0; x.len()];
        // SAFETY: x and a hold 3 numbers per body
        check(unsafe { sys::myb_nbody_acceleration(self.raw.as_ptr().cast(), t, x.as_ptr(), a.as_mut_ptr()) })?;
        Ok(a)
    }

    /// The accelerations at the positions x and velocities v (with the 1PN terms, if set).
    pub fn acceleration_xv(&self, t: f64, x: &[f64], v: &[f64]) -> Result<Vec<f64>> {
        self.sized(x)?;
        self.sized(v)?;
        let mut a = vec![0.0; x.len()];
        // SAFETY: x, v and a hold 3 numbers per body
        check(unsafe {
            sys::myb_nbody_acceleration_xv(self.raw.as_ptr().cast(), t, x.as_ptr(), v.as_ptr(), a.as_mut_ptr())
        })?;
        Ok(a)
    }

    /// The Newtonian invariants of the massive bodies (`myb_nbody_invariants_of`).
    pub fn invariants(&self, x: &[f64], v: &[f64]) -> Result<NBodyInvariants> {
        self.sized(x)?;
        self.sized(v)?;
        let mut out = NBodyInvariants::default();
        // SAFETY: x and v hold 3 numbers per body
        check(unsafe { sys::myb_nbody_invariants_of(self.raw.as_ptr(), x.as_ptr(), v.as_ptr(), &mut out) })?;
        Ok(out)
    }

    /// The 1PN invariants (`myb_nbody_pn_invariants_of`; relativity must be set).
    pub fn pn_invariants(&self, x: &[f64], v: &[f64]) -> Result<NBodyPnInvariants> {
        self.sized(x)?;
        self.sized(v)?;
        let mut out = NBodyPnInvariants::default();
        // SAFETY: x and v hold 3 numbers per body
        check(unsafe { sys::myb_nbody_pn_invariants_of(self.raw.as_ptr(), x.as_ptr(), v.as_ptr(), &mut out) })?;
        Ok(out)
    }

    /// The fast multipole method of the order (1 ..= 20; 0: exact sums again) and the opening
    /// angle theta (0 < theta < 1).
    pub fn set_approximation(&self, order: u32, theta: f64) -> Result<()> {
        let order = c_int::try_from(order).or_else(|_| domain())?;
        // SAFETY: the system is alive
        check(unsafe { sys::myb_nbody_set_approximation(self.raw.as_ptr(), order, theta) })
    }

    /// The order and the opening angle of the approximation; (0, 0) for exact sums.
    pub fn approximation(&self) -> (i32, f64) {
        let mut theta = 0.0;
        // SAFETY: the system is alive
        let order = unsafe { sys::myb_nbody_approximation(self.raw.as_ptr(), &mut theta) };
        (order, theta)
    }

    pub fn as_ptr(&self) -> *mut sys::myb_nbody {
        self.raw.as_ptr()
    }
}

impl Drop for NBody {
    fn drop(&mut self) {
        if self.owned {
            // SAFETY: created by myb_nbody_create, dropped once; the pool after it
            unsafe { sys::myb_nbody_destroy(self.raw.as_ptr()) }
        }
    }
}

// ---- Wisdom-Holman --------------------------------------------------------------------------------

/// The splitting schemes of `myb/dyn.h`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum WhScheme {
    #[default]
    Dkd,
    Saba2,
    Saba3,
    Saba4,
    Sbab1,
    Sbab2,
    Sbab3,
    Sbab4,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WhOptions {
    /// The step (s), its sign the direction.
    pub h: f64,
    pub scheme: WhScheme,
    /// 0, or an odd order 3 ..= 17 of the symplectic corrector (Dkd only).
    pub corrector: u32,
}

/// Kernel steps, evaluations of the accelerations and Kepler drifts.
pub type WhStats = sys::myb_wh_stats;

/// A Wisdom-Holman integrator of an [`NBody`] in Jacobi coordinates.
pub struct Wh {
    raw: NonNull<sys::myb_wh>,
    nbody: Rc<NBody>,
}

impl Wh {
    pub fn new(system: &Rc<NBody>, options: &WhOptions) -> Result<Self> {
        let made = sys::myb_wh_options {
            h: options.h,
            scheme: match options.scheme {
                WhScheme::Dkd => sys::MYB_WH_DKD,
                WhScheme::Saba2 => sys::MYB_WH_SABA2,
                WhScheme::Saba3 => sys::MYB_WH_SABA3,
                WhScheme::Saba4 => sys::MYB_WH_SABA4,
                WhScheme::Sbab1 => sys::MYB_WH_SBAB1,
                WhScheme::Sbab2 => sys::MYB_WH_SBAB2,
                WhScheme::Sbab3 => sys::MYB_WH_SBAB3,
                WhScheme::Sbab4 => sys::MYB_WH_SBAB4,
            },
            corrector: c_int::try_from(options.corrector).or_else(|_| domain())?,
        };
        let mut raw = ptr::null_mut();
        // SAFETY: the system outlives the integrator (it keeps it)
        check(unsafe { sys::myb_wh_create(system.as_ptr(), &made, ptr::null(), &mut raw) })?;
        Ok(Wh { raw: created(raw)?, nbody: Rc::clone(system) })
    }

    /// Starts at t0 from the inertial positions and velocities.
    pub fn start(&mut self, t0: f64, x: &[f64], v: &[f64]) -> Result<()> {
        self.nbody.sized(x)?;
        self.nbody.sized(v)?;
        // SAFETY: x and v hold 3 numbers per body
        check(unsafe { sys::myb_wh_start(self.raw.as_ptr(), t0, x.as_ptr(), v.as_ptr()) })
    }

    pub fn step(&mut self) -> Result<()> {
        // SAFETY: the integrator is alive
        check(unsafe { sys::myb_wh_step(self.raw.as_ptr()) })
    }

    /// Steps to the node nearest to t_end.
    pub fn integrate(&mut self, t_end: f64) -> Result<()> {
        // SAFETY: the integrator is alive
        check(unsafe { sys::myb_wh_integrate(self.raw.as_ptr(), t_end) })
    }

    /// The time of the current node.
    pub fn time(&self) -> f64 {
        // SAFETY: the integrator is alive
        unsafe { sys::myb_wh_time(self.raw.as_ptr()) }
    }

    fn pair(&mut self, call: impl FnOnce(*mut sys::myb_wh, *mut f64, *mut f64) -> sys::myb_status) -> Result<(Vec<f64>, Vec<f64>)> {
        let n = 3 * self.nbody.count();
        let (mut x, mut v) = (vec![0.0; n], vec![0.0; n]);
        check(call(self.raw.as_ptr(), x.as_mut_ptr(), v.as_mut_ptr()))?;
        Ok((x, v))
    }

    /// The corrected inertial positions and velocities at the current node.
    pub fn state(&mut self) -> Result<(Vec<f64>, Vec<f64>)> {
        // SAFETY: the outputs hold 3 numbers per body
        self.pair(|wh, x, v| unsafe { sys::myb_wh_state(wh, x, v) })
    }

    /// The state at t within a step of the node.
    pub fn dense(&mut self, t: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        // SAFETY: the outputs hold 3 numbers per body
        self.pair(|wh, x, v| unsafe { sys::myb_wh_dense(wh, t, x, v) })
    }

    pub fn statistics(&self) -> WhStats {
        // SAFETY: the integrator is alive
        unsafe { sys::myb_wh_statistics(self.raw.as_ptr()) }
    }
}

impl Drop for Wh {
    fn drop(&mut self) {
        // SAFETY: created by myb_wh_create, dropped once, before its system
        unsafe { sys::myb_wh_destroy(self.raw.as_ptr()) }
    }
}

// ---- Polyhedra ------------------------------------------------------------------------------------

/// A closed triangular mesh (`myb/grav.h`): the shape of a body for light and contacts, and
/// its field for a constant density.
pub struct Polyhedron {
    raw: NonNull<sys::myb_polyhedron>,
}

impl Polyhedron {
    /// Vertices (m), faces counter-clockwise from outside, and G rho (s^-2).
    pub fn new(vertices: &[Vec3], faces: &[[usize; 3]], g_rho: f64) -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: the arrays hold 3 numbers per vertex and 3 indices per face, in a row
        check(unsafe {
            sys::myb_polyhedron_create(
                vertices.as_ptr().cast(),
                vertices.len(),
                faces.as_ptr().cast(),
                faces.len(),
                g_rho,
                ptr::null(),
                &mut raw,
            )
        })?;
        Ok(Polyhedron { raw: created(raw)? })
    }

    /// The volume (m^3) and the centre of mass (m).
    pub fn volume(&self) -> (f64, Vec3) {
        let mut center = [0.0; 3];
        // SAFETY: the mesh is alive; center holds three numbers
        let volume = unsafe { sys::myb_polyhedron_volume(self.raw.as_ptr(), center.as_mut_ptr()) };
        (volume, center)
    }

    /// The potential, the acceleration and the Laplacian at r in the frame of the mesh.
    pub fn eval(&self, r: &Vec3) -> Result<(f64, Vec3, f64)> {
        let (mut potential, mut a, mut laplacian) = (0.0, [0.0; 3], 0.0);
        // SAFETY: valid places for the outputs
        check(unsafe { sys::myb_polyhedron_eval(self.raw.as_ptr(), r.as_ptr(), &mut potential, a.as_mut_ptr(), &mut laplacian) })?;
        Ok((potential, a, laplacian))
    }

    pub fn as_ptr(&self) -> *const sys::myb_polyhedron {
        self.raw.as_ptr()
    }
}

impl Drop for Polyhedron {
    fn drop(&mut self) {
        // SAFETY: created by myb_polyhedron_create, dropped once
        unsafe { sys::myb_polyhedron_destroy(self.raw.as_ptr()) }
    }
}

// ---- A system of bodies with models ----------------------------------------------------------------

/// An integrated rotation: moments per unit mass, distortion, core (`myb_system_rotation`).
pub type SystemRotation = sys::myb_system_rotation;
/// The answer of a body to tides by band (`myb_system_tide`).
pub type SystemTide = sys::myb_system_tide;
/// The terms beyond 1PN of a pair: `orders` or-ed from the `PN_*` flags below (those of
/// `myb/rel.h`).
pub type SystemPnPair = sys::myb_system_pn_pair;

/// The orders of the post-Newtonian terms of a pair (`myb/rel.h`), for [`SystemPnPair::orders`].
pub const PN_1: u32 = sys::MYB_PN_1;
pub const PN_2: u32 = sys::MYB_PN_2;
pub const PN_25: u32 = sys::MYB_PN_25;
pub const PN_3: u32 = sys::MYB_PN_3;
pub const PN_35: u32 = sys::MYB_PN_35;
/// Spin-orbit, leading order.
pub const PN_SO: u32 = sys::MYB_PN_SO;
/// Spin-spin, leading order.
pub const PN_SS: u32 = sys::MYB_PN_SS;
/// PN_1 .. PN_SS.
pub const PN_ALL: u32 = sys::MYB_PN_ALL;
pub const PN_4: u32 = sys::MYB_PN_4;
/// Radiation reaction at 4.5PN.
pub const PN_4RR: u32 = sys::MYB_PN_4RR;
pub const PN_SO25: u32 = sys::MYB_PN_SO25;
pub const PN_SO35: u32 = sys::MYB_PN_SO35;
/// Spin-orbit radiation reaction.
pub const PN_SORR: u32 = sys::MYB_PN_SORR;
/// Spin-spin at next-to-leading order (1/c^6), with its precession.
pub const PN_SS3: u32 = sys::MYB_PN_SS3;
/// Spin-orbit radiation reaction at next-to-leading order (1/c^10).
pub const PN_SORR5: u32 = sys::MYB_PN_SORR5;
/// The spin-orbit part of the dissipative tail, quasi-circular orbits.
pub const PN_SOTAIL: u32 = sys::MYB_PN_SOTAIL;
/// The tails of PN_4 and PN_4RR on any osculating conic.
pub const PN_TAILX: u32 = sys::MYB_PN_TAILX;
/// What a merger takes out of the motion, with G as a factor.
pub type SystemMerger = sys::myb_system_merger;
/// The optical properties of a surface (`myb/force.h`).
pub type Plate = sys::myb_plate;
/// The law of a contact (`myb/contact.h`): the normal spring and damping, the friction —
/// smooth (`slip`) or sticking (`kind` = [`ContactFriction::Stick`]`.raw()`, with
/// `tangential_stiffness` and `tangential_damping`).
pub type ContactLaw = sys::myb_contact_law;
/// A contact point of a body: its offset in the body frame and its radius.
pub type ContactPoint = sys::myb_contact_point;
/// The force and torque per unit mass of a contact, the energy of its springs, the gap and
/// the points that press.
pub type ContactResult = sys::myb_contact_result;

/// The friction of a contact law.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum ContactFriction {
    /// Coulomb's made smooth: a function of the state.
    #[default]
    Smooth,
    /// A spring with a slip limit: a state of its own, the deflection of each point.
    Stick,
}

impl ContactFriction {
    /// The value of `ContactLaw::kind`.
    pub fn raw(self) -> sys::myb_contact_friction {
        match self {
            ContactFriction::Smooth => sys::MYB_CONTACT_FRICTION_SMOOTH,
            ContactFriction::Stick => sys::MYB_CONTACT_FRICTION_STICK,
        }
    }

    /// The friction of a law (`None`: a value this crate does not know).
    pub fn of(law: &ContactLaw) -> Option<Self> {
        match law.kind {
            sys::MYB_CONTACT_FRICTION_SMOOTH => Some(ContactFriction::Smooth),
            sys::MYB_CONTACT_FRICTION_STICK => Some(ContactFriction::Stick),
            _ => None,
        }
    }
}
/// Accepted and rejected steps, evaluations of IAS15.
pub type Ias15Stats = sys::myb_ias15_stats;

/// What an action sees: the body acted on and all the bodies as they are.
pub struct SystemView<'a> {
    pub body: usize,
    pub bodies: usize,
    /// Positions of all the bodies, 3 each (m).
    pub x: &'a [f64],
    /// Velocities (m/s).
    pub v: &'a [f64],
    /// The GM of the moment.
    pub gm: &'a [f64],
    /// The rotation inertial -> body (rows), if the body has a frame.
    pub frame: Option<&'a [f64; 9]>,
}

type Action = dyn Fn(f64, &SystemView) -> Result<(Vec3, Vec3)>;

unsafe extern "C" fn call_action(
    context: *mut c_void,
    t: f64,
    view: *const sys::myb_system_view,
    acceleration: *mut f64,
    torque: *mut f64,
) -> sys::myb_status {
    // SAFETY: the context is the Box<Action> the system keeps; the view and its arrays are the
    // library's, valid for the call
    let action = unsafe { &*(context as *const Box<Action>) };
    let raw = unsafe { &*view };
    let bodies = raw.bodies;
    let seen = unsafe {
        SystemView {
            body: raw.body,
            bodies,
            x: std::slice::from_raw_parts(raw.x, 3 * bodies),
            v: std::slice::from_raw_parts(raw.v, 3 * bodies),
            gm: std::slice::from_raw_parts(raw.gm, bodies),
            frame: (!raw.frame.is_null()).then(|| &*(raw.frame as *const [f64; 9])),
        }
    };
    match action(t, &seen) {
        Ok((a, n)) => {
            unsafe {
                *(acceleration as *mut Vec3) = a;
                *(torque as *mut Vec3) = n;
            }
            sys::MYB_OK
        }
        Err(error) => error.status.raw(),
    }
}

/// Light of a body on the mesh of another, with the torque of the moment
/// (`myb_system_light`).
pub struct SystemLight {
    /// The body that shines.
    pub source: usize,
    /// The mesh of the lit body, in its frame about its centre of mass.
    pub shape: Rc<Polyhedron>,
    /// Its optical properties (area and normal are not used).
    pub surface: Plate,
    /// W m^-2 at `distance` (m) from the source.
    pub flux: f64,
    pub distance: f64,
    /// m/s.
    pub c: f64,
    /// Of the lit body, kg.
    pub mass: f64,
    pub shadows: bool,
    /// The thermal parameter at the rotation frequency; 0: re-emitted at once.
    pub theta: f64,
    /// With theta > 0: unit, in the body frame, the sense of the rotation.
    pub spin_axis: Vec3,
    /// With theta > 0: steps of a turn (0: 64).
    pub phases: u32,
}

/// The surface of the ground of a contact, in its frame.
pub enum ContactSurface {
    /// A point of the plane and its normal towards the body.
    Plane { point: Vec3, normal: Vec3 },
    /// The semi-axes along x, y, z, centred at the origin.
    Ellipsoid { radii: Vec3 },
    /// A closed mesh.
    Mesh(Rc<Polyhedron>),
}

impl ContactSurface {
    /// The structure of the C interface; a mesh is pointed to, alive as long as `self`.
    fn raw(&self) -> sys::myb_contact_surface {
        let mut surface = sys::myb_contact_surface::default();
        match self {
            ContactSurface::Plane { point, normal } => {
                surface.shape = sys::MYB_CONTACT_PLANE;
                surface.point = *point;
                surface.normal = *normal;
            }
            ContactSurface::Ellipsoid { radii } => {
                surface.shape = sys::MYB_CONTACT_ELLIPSOID;
                surface.radii = *radii;
            }
            ContactSurface::Mesh(mesh) => {
                surface.shape = sys::MYB_CONTACT_MESH;
                surface.mesh = mesh.as_ptr();
            }
        }
        surface
    }
}

/// The force of a contact with the smooth friction on a body at x, v (in the axes of the
/// surface), turned by `attitude` (rows, surface -> body) at the angular velocity omega (body
/// axes) (`myb_contact_force`); a law whose friction sticks needs its state:
/// [`contact_force_stick`].
pub fn contact_force(
    surface: &ContactSurface,
    law: &ContactLaw,
    points: &[ContactPoint],
    x: &Vec3,
    v: &Vec3,
    attitude: &[f64; 9],
    omega: &Vec3,
) -> Result<ContactResult> {
    let made = surface.raw();
    let mut result = ContactResult::default();
    // SAFETY: the surface (and its mesh), the law and the points are alive for the call
    check(unsafe {
        sys::myb_contact_force(&made, law, points.as_ptr(), points.len(), x.as_ptr(), v.as_ptr(), attitude.as_ptr(), omega.as_ptr(), &mut result)
    })?;
    Ok(result)
}

/// The same with the state of a friction that sticks (`myb_contact_force_stick`): the
/// deflection of each point (3 numbers a point, axes of the surface; empty with the smooth
/// friction). Returns the result and the rates of the deflections.
#[allow(clippy::too_many_arguments)]
pub fn contact_force_stick(
    surface: &ContactSurface,
    law: &ContactLaw,
    points: &[ContactPoint],
    x: &Vec3,
    v: &Vec3,
    attitude: &[f64; 9],
    omega: &Vec3,
    deflection: &[f64],
) -> Result<(ContactResult, Vec<f64>)> {
    if !deflection.is_empty() && deflection.len() != 3 * points.len() {
        return domain();
    }
    let made = surface.raw();
    let mut result = ContactResult::default();
    let mut rate = vec![0.0; 3 * points.len()];
    // SAFETY: as above; the deflection, if any, and its rate hold three numbers a point
    check(unsafe {
        sys::myb_contact_force_stick(
            &made,
            law,
            points.as_ptr(),
            points.len(),
            x.as_ptr(),
            v.as_ptr(),
            attitude.as_ptr(),
            omega.as_ptr(),
            if deflection.is_empty() { ptr::null() } else { deflection.as_ptr() },
            crate::out_ptr(&mut rate),
            &mut result,
        )
    })?;
    Ok((result, rate))
}

/// A body on the surface of another (`myb_system_contact`).
pub struct SystemContact {
    pub surface: ContactSurface,
    /// Per unit mass of the body.
    pub law: ContactLaw,
    /// Of the body, in its frame.
    pub points: Vec<ContactPoint>,
}

/// The options of IAS15 that a system uses: precision, step limits, budget; zeros stand for
/// the defaults of the library.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Ias15Options {
    pub epsilon: f64,
    pub h_initial: f64,
    pub h_max: f64,
    pub h_min: f64,
    pub max_steps: u64,
    /// The timescale of the step (Jerk for an acceleration that passes through zero).
    pub timescale: Ias15Timescale,
}

impl Ias15Options {
    fn raw(&self) -> sys::myb_ias15_options {
        sys::myb_ias15_options {
            epsilon: self.epsilon,
            h_initial: self.h_initial,
            h_max: self.h_max,
            h_min: self.h_min,
            max_steps: self.max_steps.min(i64::MAX as u64) as i64,
            timescale: self.timescale.raw(),
            ..Default::default()
        }
    }
}

/// The derivatives of the state at t1: the state, d y1 / d y0 (2 n x 2 n, rows) and d y1 / d gm
/// of the bodies asked for (2 n x count, rows).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SystemPartials {
    pub y1: Vec<f64>,
    pub phi: Vec<f64>,
    pub sensitivity: Vec<f64>,
}

/// Bodies with their models (`myb/system.h`). Configure it, then integrate it with a
/// [`SystemIntegrator`] (which keeps it through an `Rc`).
pub struct System {
    raw: NonNull<sys::myb_system>,
    nbody: NBody,
    kept: Kept,
    /// The contacts in place, in their order: body, ground, points.
    contacts: RefCell<Vec<(usize, usize, usize)>>,
    /// The most points a contact has had: the size of the buffers of its deflections.
    most_points: std::cell::Cell<usize>,
}

impl System {
    /// gm: the bodies (at least one).
    pub fn new(gm: &[f64]) -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: gm holds gm.len() numbers; a valid place for the pointer
        check(unsafe { sys::myb_system_create(gm.len(), gm.as_ptr(), ptr::null(), &mut raw) })?;
        let raw = created(raw)?;
        // SAFETY: the point masses live inside the system, as long as it does
        let points = created(unsafe { sys::myb_system_nbody(raw.as_ptr()) })?;
        Ok(System {
            raw,
            nbody: NBody { raw: points, owned: false, pool: RefCell::new(None) },
            kept: Kept::default(),
            contacts: RefCell::new(Vec::new()),
            most_points: std::cell::Cell::new(0),
        })
    }

    pub fn bodies(&self) -> usize {
        // SAFETY: the system is alive
        unsafe { sys::myb_system_bodies(self.raw.as_ptr()) }
    }

    /// The point masses: their graph, relativity, threads and approximation.
    pub fn nbody(&self) -> &NBody {
        &self.nbody
    }

    /// The body-fixed frame of a body from a provider (on the time of the integration).
    pub fn set_orientation(&self, body: usize, orientation: impl OrientationSource + 'static) -> Result<()> {
        let kept = self.kept.keep(orientation);
        // SAFETY: the provider lives in self.kept, as long as the system
        let (function, context) = unsafe { (*kept).as_raw() };
        check(unsafe { sys::myb_system_set_orientation(self.raw.as_ptr(), body, function, context) })
    }

    /// The figure of a body up to a degree (`None`: all of the field).
    pub fn set_field(&self, body: usize, field: &Rc<GravityField>, degree: Option<u32>) -> Result<()> {
        let degree = degree.map_or(Ok(0), c_int::try_from).or_else(|_| domain())?;
        self.kept.keep(Rc::clone(field));
        // SAFETY: the field is kept alive by the system
        check(unsafe { sys::myb_system_set_field(self.raw.as_ptr(), body, field.as_ptr(), degree) })
    }

    /// The rotation of a body integrated with the motion.
    pub fn set_rotation(&self, body: usize, rotation: &SystemRotation) -> Result<()> {
        // SAFETY: the structure is copied
        check(unsafe { sys::myb_system_set_rotation(self.raw.as_ptr(), body, rotation) })
    }

    /// The answer of a body to tides.
    pub fn set_tide(&self, body: usize, tide: &SystemTide) -> Result<()> {
        // SAFETY: the structure is copied
        check(unsafe { sys::myb_system_set_tide(self.raw.as_ptr(), body, tide) })
    }

    /// The figure of `body` acts on `point`, which acts back.
    pub fn add_field(&self, body: usize, point: usize) -> Result<()> {
        // SAFETY: the library checks the indices
        check(unsafe { sys::myb_system_add_field(self.raw.as_ptr(), body, point) })
    }

    /// The quadrupoles of two bodies on each other.
    pub fn add_figures(&self, a: usize, b: usize) -> Result<()> {
        // SAFETY: the library checks the indices
        check(unsafe { sys::myb_system_add_figures(self.raw.as_ptr(), a, b) })
    }

    /// The tide `raiser` raises on `body` acts on `point`.
    pub fn add_tide(&self, body: usize, raiser: usize, point: usize) -> Result<()> {
        // SAFETY: the library checks the indices
        check(unsafe { sys::myb_system_add_tide(self.raw.as_ptr(), body, raiser, point) })
    }

    /// `raiser` distorts `body`, whose rotation is integrated.
    pub fn add_distortion(&self, body: usize, raiser: usize) -> Result<()> {
        // SAFETY: the library checks the indices
        check(unsafe { sys::myb_system_add_distortion(self.raw.as_ptr(), body, raiser) })
    }

    /// The terms beyond 1PN of a pair of bodies.
    pub fn add_pn_pair(&self, a: usize, b: usize, pair: &SystemPnPair) -> Result<()> {
        // SAFETY: the structure is copied
        check(unsafe { sys::myb_system_add_pn_pair(self.raw.as_ptr(), a, b, pair) })
    }

    /// Anything else that pushes or turns a body: the closure gives the acceleration (m/s^2,
    /// inertial axes) and the torque per unit mass (m^2/s^2, body frame).
    pub fn add_action(&self, body: usize, action: impl Fn(f64, &SystemView) -> Result<(Vec3, Vec3)> + 'static) -> Result<()> {
        let boxed: Box<Action> = Box::new(action);
        let context = self.kept.keep(boxed);
        // SAFETY: call_action reads the Box<Action> that self.kept holds
        check(unsafe { sys::myb_system_add_action(self.raw.as_ptr(), body, Some(call_action), context.cast()) })
    }

    /// The light of a source on the mesh of `body` (`myb_system_light_action`).
    pub fn add_light(&self, body: usize, light: &SystemLight) -> Result<()> {
        let made = sys::myb_system_light {
            source: light.source,
            shape: light.shape.as_ptr(),
            surface: light.surface,
            flux: light.flux,
            distance: light.distance,
            c: light.c,
            mass: light.mass,
            shadows: c_int::from(light.shadows),
            theta: light.theta,
            spin_axis: light.spin_axis,
            phases: c_int::try_from(light.phases).or_else(|_| domain())?,
        };
        self.kept.keep(Rc::clone(&light.shape));
        let context = self.kept.keep(made);
        // SAFETY: the structure and its mesh live as long as the system
        check(unsafe { sys::myb_system_add_action(self.raw.as_ptr(), body, Some(sys::myb_system_light_action), context.cast()) })
    }

    /// The spin of `body` (angular momentum per unit mass, m^2/s, axes of the system) drags
    /// `point` (Lense-Thirring; relativity must be set).
    pub fn add_frame_dragging(&self, body: usize, point: usize, spin: &Vec3) -> Result<()> {
        // SAFETY: the spin is copied
        check(unsafe { sys::myb_system_add_frame_dragging(self.raw.as_ptr(), body, point, spin.as_ptr()) })
    }

    /// gm (1 + rate (t - epoch)) for the body.
    pub fn set_mass_rate(&self, body: usize, rate: f64, epoch: f64) -> Result<()> {
        // SAFETY: the system is alive
        check(unsafe { sys::myb_system_set_mass_rate(self.raw.as_ptr(), body, rate, epoch) })
    }

    /// A body whose motion a provider gives (at seconds of the integration as instants); it
    /// pulls the bodies and feels nothing. Returns its number among the guided bodies.
    pub fn add_guided(&self, gm: f64, state: impl StateSource + 'static) -> Result<usize> {
        let kept = self.kept.keep(state);
        // SAFETY: the provider lives in self.kept, as long as the system
        let (function, context) = unsafe { (*kept).as_raw() };
        let mut index = 0usize;
        check(unsafe { sys::myb_system_add_guided(self.raw.as_ptr(), gm, function, context, &mut index) })?;
        Ok(index)
    }

    pub fn guided(&self) -> usize {
        // SAFETY: the system is alive
        unsafe { sys::myb_system_guided(self.raw.as_ptr()) }
    }

    /// `body` touches the surface of `ground` with its contact points. A friction that sticks
    /// adds 3 coordinates a point to the state, the deflections of its springs
    /// ([`System::pack_contact`]).
    pub fn add_contact(&self, body: usize, ground: usize, contact: &SystemContact) -> Result<()> {
        let surface = contact.surface.raw();
        if let ContactSurface::Mesh(mesh) = &contact.surface {
            self.kept.keep(Rc::clone(mesh));
        }
        let points = contact.points.clone().into_boxed_slice();
        let made = sys::myb_system_contact { surface, law: contact.law, points: points.as_ptr(), point_count: points.len() };
        self.kept.keep(points); // the points stay where they are
        // SAFETY: the points live in self.kept; the surface, law and count are copied
        check(unsafe { sys::myb_system_add_contact(self.raw.as_ptr(), body, ground, &made) })?;
        self.contacts.borrow_mut().push((body, ground, contact.points.len()));
        self.most_points.set(self.most_points.get().max(contact.points.len()));
        Ok(())
    }

    /// The points of the contact number `contact` in place, if there is one.
    fn contact_points(&self, contact: usize) -> Result<usize> {
        self.contacts.borrow().get(contact).map(|&(_, _, points)| points).ok_or(crate::Error::new(crate::Status::Domain))
    }

    /// Writes the deflections of the springs of a contact whose friction sticks (3 numbers a
    /// point, axes of its surface) into a state y.
    pub fn pack_contact(&self, contact: usize, deflection: &[f64], y: &mut [f64]) -> Result<()> {
        self.state_sized(y)?;
        if deflection.len() != 3 * self.contact_points(contact)? {
            return domain();
        }
        // a buffer for the most points any contact has: the library reads those of its own
        let mut buffer = vec![0.0; 3 * self.most_points.get()];
        buffer[..deflection.len()].copy_from_slice(deflection);
        // SAFETY: y holds the whole state; the buffer is as long as any contact needs
        check(unsafe { sys::myb_system_pack_contact(self.raw.as_ptr(), contact, buffer.as_ptr(), y.as_mut_ptr()) })
    }

    /// The deflections of the springs of a contact in y.
    pub fn unpack_contact(&self, contact: usize, y: &[f64]) -> Result<Vec<f64>> {
        self.state_sized(y)?;
        let points = self.contact_points(contact)?;
        let mut buffer = vec![0.0; 3 * self.most_points.get()];
        // SAFETY: y holds the whole state; the buffer is as long as any contact needs
        check(unsafe { sys::myb_system_unpack_contact(self.raw.as_ptr(), contact, y.as_ptr(), buffer.as_mut_ptr()) })?;
        buffer.truncate(3 * points);
        Ok(buffer)
    }

    /// The coordinates of the state: 3 per body, 9 per rotation, 6 per PN pair with spins, 3
    /// per point of a contact whose friction sticks.
    pub fn coordinates(&self) -> usize {
        // SAFETY: the system is alive
        unsafe { sys::myb_system_coordinates(self.raw.as_ptr()) }
    }

    fn state_sized(&self, y: &[f64]) -> Result<()> {
        if y.len() == 2 * self.coordinates() {
            Ok(())
        } else {
            domain()
        }
    }

    /// Writes the rotation of a body into a state y (2 n numbers).
    pub fn pack_rotation(&self, body: usize, q: &[f64; 4], omega: &Vec3, omega_core: Option<&Vec3>, y: &mut [f64]) -> Result<()> {
        self.state_sized(y)?;
        // SAFETY: y holds the whole state
        check(unsafe {
            sys::myb_system_pack_rotation(
                self.raw.as_ptr(),
                body,
                q.as_ptr(),
                omega.as_ptr(),
                omega_core.map_or(ptr::null(), |w| w.as_ptr()),
                y.as_mut_ptr(),
            )
        })
    }

    /// The attitude, the angular velocity of the mantle and that of the core of a body in y.
    pub fn unpack_rotation(&self, body: usize, y: &[f64]) -> Result<([f64; 4], Vec3, Vec3)> {
        self.state_sized(y)?;
        let (mut q, mut omega, mut core) = ([0.0; 4], [0.0; 3], [0.0; 3]);
        // SAFETY: y holds the whole state; the outputs are valid
        check(unsafe {
            sys::myb_system_unpack_rotation(self.raw.as_ptr(), body, y.as_ptr(), q.as_mut_ptr(), omega.as_mut_ptr(), core.as_mut_ptr())
        })?;
        Ok((q, omega, core))
    }

    /// The spins of a PN pair into y.
    pub fn pack_spins(&self, pair: usize, chi1: &Vec3, chi2: &Vec3, y: &mut [f64]) -> Result<()> {
        self.state_sized(y)?;
        // SAFETY: y holds the whole state
        check(unsafe { sys::myb_system_pack_spins(self.raw.as_ptr(), pair, chi1.as_ptr(), chi2.as_ptr(), y.as_mut_ptr()) })
    }

    /// The spins of a PN pair in y.
    pub fn unpack_spins(&self, pair: usize, y: &[f64]) -> Result<(Vec3, Vec3)> {
        self.state_sized(y)?;
        let (mut chi1, mut chi2) = ([0.0; 3], [0.0; 3]);
        // SAFETY: y holds the whole state
        check(unsafe { sys::myb_system_unpack_spins(self.raw.as_ptr(), pair, y.as_ptr(), chi1.as_mut_ptr(), chi2.as_mut_ptr()) })?;
        Ok((chi1, chi2))
    }

    /// The second derivatives of the n coordinates x with their rates v.
    pub fn acceleration(&self, t: f64, x: &[f64], v: &[f64]) -> Result<Vec<f64>> {
        let n = self.coordinates();
        if x.len() != n || v.len() != n {
            return domain();
        }
        let mut a = vec![0.0; n];
        // SAFETY: x, v and a hold n numbers
        check(unsafe { sys::myb_system_acceleration(self.raw.as_ptr().cast(), t, x.as_ptr(), v.as_ptr(), a.as_mut_ptr()) })?;
        Ok(a)
    }

    /// The inertia tensor (m^2, mantle frame) of an integrated rotation at the state (t, y).
    pub fn inertia(&self, body: usize, t: f64, y: &[f64]) -> Result<[f64; 9]> {
        self.state_sized(y)?;
        let mut inertia = [0.0; 9];
        // SAFETY: y holds the whole state
        check(unsafe { sys::myb_system_inertia(self.raw.as_ptr(), body, t, y.as_ptr(), inertia.as_mut_ptr()) })?;
        Ok(inertia)
    }

    /// The derivatives of the state at t1 by the state at t0 and by the GM of `bodies`.
    pub fn partials(&self, options: &Ias15Options, t0: f64, y0: &[f64], t1: f64, bodies: &[usize]) -> Result<SystemPartials> {
        self.state_sized(y0)?;
        let m = y0.len();
        let mut out = SystemPartials { y1: vec![0.0; m], phi: vec![0.0; m * m], sensitivity: vec![0.0; m * bodies.len()] };
        let made = options.raw();
        // SAFETY: the outputs hold 2n, (2n)^2 and 2n count numbers
        check(unsafe {
            sys::myb_system_partials(
                self.raw.as_ptr(),
                &made,
                t0,
                y0.as_ptr(),
                t1,
                bodies.as_ptr(),
                bodies.len(),
                out.y1.as_mut_ptr(),
                out.phi.as_mut_ptr(),
                crate::out_ptr(&mut out.sensitivity),
            )
        })?;
        Ok(out)
    }

    /// The radius of a body, which makes it a sphere that may touch others.
    pub fn set_radius(&self, body: usize, radius: f64) -> Result<()> {
        // SAFETY: the library checks the index
        check(unsafe { sys::myb_system_set_radius(self.raw.as_ptr(), body, radius) })
    }

    /// `None` if there is no such body.
    pub fn radius(&self, body: usize) -> Option<f64> {
        // SAFETY: the body is one of the system
        (body < self.bodies()).then(|| unsafe { sys::myb_system_radius(self.raw.as_ptr(), body) })
    }

    /// False for a body merged into another (or one the system does not have).
    pub fn present(&self, body: usize) -> bool {
        // SAFETY: the body is one of the system
        body < self.bodies() && unsafe { sys::myb_system_present(self.raw.as_ptr(), body) } != 0
    }

    pub fn as_ptr(&self) -> *mut sys::myb_system {
        self.raw.as_ptr()
    }
}

impl Drop for System {
    fn drop(&mut self) {
        // SAFETY: created by myb_system_create, dropped once; what it points to after it
        unsafe { sys::myb_system_destroy(self.raw.as_ptr()) }
    }
}

/// Why [`SystemIntegrator::advance`] stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SystemEventKind {
    /// t_end reached.
    Arrived,
    /// The surfaces of a and b touch, approaching.
    Contact,
    /// a and b are at their closest, within `watch`.
    Approach,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SystemEvent {
    pub kind: SystemEventKind,
    /// a < b.
    pub a: usize,
    pub b: usize,
    /// The time of the integrator.
    pub t: f64,
    /// Of the centres (m).
    pub distance: f64,
}

/// IAS15 over a [`System`] (`myb_system_integrator`), with the encounters of its bodies.
pub struct SystemIntegrator {
    raw: NonNull<sys::myb_ias15>,
    system: Rc<System>,
    n: usize,
}

impl SystemIntegrator {
    /// For the system as it stands (the levels that add coordinates in place).
    pub fn new(system: &Rc<System>, options: &Ias15Options) -> Result<Self> {
        let made = options.raw();
        let mut raw = ptr::null_mut();
        // SAFETY: a valid place for the pointer
        check(unsafe { sys::myb_system_integrator(system.as_ptr(), &made, ptr::null(), &mut raw) })?;
        Ok(SystemIntegrator { raw: created(raw)?, system: Rc::clone(system), n: system.coordinates() })
    }

    pub fn system(&self) -> &Rc<System> {
        &self.system
    }

    /// Starts at t0 from the state y0 (2 n numbers) with the accelerations of the system.
    pub fn start(&mut self, t0: f64, y0: &[f64]) -> Result<()> {
        if y0.len() != 2 * self.n {
            return domain();
        }
        // SAFETY: the system outlives the integrator (it keeps it); y0 holds 2n numbers
        check(unsafe {
            sys::myb_ias15_start(self.raw.as_ptr(), Some(sys::myb_system_acceleration), self.system.as_ptr().cast(), t0, y0.as_ptr())
        })
    }

    /// One accepted step towards t_limit.
    pub fn step(&mut self, t_limit: f64) -> Result<()> {
        // SAFETY: the integrator is alive
        check(unsafe { sys::myb_ias15_step(self.raw.as_ptr(), t_limit) })
    }

    /// Steps until t_end exactly (no encounters looked for).
    pub fn integrate(&mut self, t_end: f64) -> Result<()> {
        // SAFETY: the integrator is alive
        check(unsafe { sys::myb_ias15_integrate(self.raw.as_ptr(), t_end) })
    }

    pub fn time(&self) -> f64 {
        // SAFETY: the integrator is alive
        unsafe { sys::myb_ias15_time(self.raw.as_ptr()) }
    }

    /// The state (2 n numbers) where the integrator stands.
    pub fn state(&self) -> Vec<f64> {
        // SAFETY: the integrator holds 2n numbers of state
        unsafe { std::slice::from_raw_parts(sys::myb_ias15_state(self.raw.as_ptr()), 2 * self.n) }.to_vec()
    }

    /// The state at t within the last step.
    pub fn dense(&mut self, t: f64) -> Result<Vec<f64>> {
        let mut y = vec![0.0; 2 * self.n];
        // SAFETY: y holds 2n numbers
        check(unsafe { sys::myb_ias15_dense(self.raw.as_ptr(), t, y.as_mut_ptr()) })?;
        Ok(y)
    }

    pub fn statistics(&self) -> Ias15Stats {
        // SAFETY: the integrator is alive
        unsafe { sys::myb_ias15_statistics(self.raw.as_ptr()) }
    }

    /// Steps towards t_end, stopping at the first contact of two bodies with radii or, with
    /// watch > 0, at the closest approach of two bodies within `watch`.
    pub fn advance(&mut self, t_end: f64, watch: f64) -> Result<SystemEvent> {
        let mut event = sys::myb_system_event::default();
        // SAFETY: the integrator was made for this system
        check(unsafe { sys::myb_system_advance(self.system.as_ptr(), self.raw.as_ptr(), t_end, watch, &mut event) })?;
        Ok(SystemEvent {
            kind: match event.kind {
                sys::MYB_SYSTEM_CONTACT => SystemEventKind::Contact,
                sys::MYB_SYSTEM_APPROACH => SystemEventKind::Approach,
                _ => SystemEventKind::Arrived,
            },
            a: event.a,
            b: event.b,
            t: event.t,
            distance: event.distance,
        })
    }

    /// Merges body b into body a; the integrator starts again.
    pub fn merge(&mut self, a: usize, b: usize) -> Result<SystemMerger> {
        let mut merger = SystemMerger::default();
        // SAFETY: the integrator was made for this system
        check(unsafe { sys::myb_system_merge(self.system.as_ptr(), self.raw.as_ptr(), a, b, &mut merger) })?;
        // every level that names b is dropped: so are its contacts, and the rest renumbered
        self.system.contacts.borrow_mut().retain(|&(body, ground, _)| body != b && ground != b);
        Ok(merger)
    }

    /// Two hard spheres rebound with the restitution (0 ..= 1); the integrator starts again.
    pub fn bounce(&mut self, a: usize, b: usize, restitution: f64) -> Result<()> {
        // SAFETY: the integrator was made for this system
        check(unsafe { sys::myb_system_bounce(self.system.as_ptr(), self.raw.as_ptr(), a, b, restitution) })
    }

    pub fn as_ptr(&self) -> *mut sys::myb_ias15 {
        self.raw.as_ptr()
    }
}

impl Drop for SystemIntegrator {
    fn drop(&mut self) {
        // SAFETY: created by myb_system_integrator, dropped once, before its system
        unsafe { sys::myb_ias15_destroy(self.raw.as_ptr()) }
    }
}
