// SPDX-License-Identifier: BSD-3-Clause
//! The last of phase 21 through the safe layer against the same C calls, bit for bit: the
//! break schemes of a model, the timescale of IAS15, the friction that sticks (contact.h and
//! its state in a system), the new orders of a PN pair; failures as `Err` with the status.
use std::fmt::Debug;
use std::ptr;
use std::rc::Rc;

use myb::sys;
use myb::{
    BreakScheme, ContactFriction, ContactLaw, ContactPoint, ContactSurface, ForceModel, Ias15Options, Ias15Timescale, Orbit,
    OrbitOptions, Ppn, State, Status, System, SystemContact, SystemIntegrator, SystemPnPair,
};

const GM: f64 = 3.986004418e14;

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn same<T: Debug>(a: &T, b: &T) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

fn orbit_state() -> State {
    myb::state_from_elements(GM, &myb::Elements { p: 7.0e6, e: 0.05, i: 0.9, raan: 0.3, argp: 0.2, nu: 0.1 }).unwrap()
}

#[test]
fn break_schemes_and_timescales_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let model = Rc::new(ForceModel::new()?);
    model.add_point_mass("earth", GM)?;
    assert_eq!(model.break_scheme(), BreakScheme::Retake);
    model.set_break_scheme(BreakScheme::Predict);
    assert_eq!(model.break_scheme(), BreakScheme::Predict);
    // SAFETY: the model is alive
    assert_eq!(unsafe { sys::myb_force_model_break_scheme(model.as_ptr()) }, sys::MYB_BREAKS_PREDICT);

    let y0 = orbit_state();
    let options = OrbitOptions { timescale: Ias15Timescale::Jerk, ..Default::default() };
    let mut orbit = Orbit::new(&model, &options)?;
    orbit.start(0.0, &y0, 6000.0)?;
    let end = orbit.state(6000.0)?;
    let raw_options = sys::myb_orbit_options { timescale: sys::MYB_IAS15_TIMESCALE_JERK, ..Default::default() };
    let mut raw_end = [0.0; 6];
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    let raw_steps = unsafe {
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_orbit_create(model.as_ptr(), &raw_options, ptr::null(), &mut raw), sys::MYB_OK);
        assert_eq!(sys::myb_orbit_start(raw, 0.0, y0.as_ptr(), 6000.0), sys::MYB_OK);
        assert_eq!(sys::myb_orbit_state(raw, 6000.0, raw_end.as_mut_ptr()), sys::MYB_OK);
        let steps = sys::myb_orbit_steps(raw);
        sys::myb_orbit_destroy(raw);
        steps
    };
    assert_eq!((bits(&end), orbit.steps()), (bits(&raw_end), raw_steps));
    model.set_break_scheme(BreakScheme::Retake);
    assert_eq!(model.break_scheme(), BreakScheme::Retake);

    // the timescale of the integrator of a system
    let gm = [GM, 1.0];
    let mut y = vec![0.0; 12];
    y[3..6].copy_from_slice(&y0[..3]);
    y[9..12].copy_from_slice(&y0[3..]);
    let system = Rc::new(System::new(&gm)?);
    let jerk = Ias15Options { timescale: Ias15Timescale::Jerk, ..Default::default() };
    let mut integrator = SystemIntegrator::new(&system, &jerk)?;
    integrator.start(0.0, &y)?;
    integrator.integrate(3000.0)?;
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    let raw_state = unsafe {
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_system_create(2, gm.as_ptr(), ptr::null(), &mut raw), sys::MYB_OK);
        let options = sys::myb_ias15_options { timescale: sys::MYB_IAS15_TIMESCALE_JERK, ..Default::default() };
        let mut ias15 = ptr::null_mut();
        assert_eq!(sys::myb_system_integrator(raw, &options, ptr::null(), &mut ias15), sys::MYB_OK);
        assert_eq!(sys::myb_ias15_start(ias15, Some(sys::myb_system_acceleration), raw.cast(), 0.0, y.as_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_ias15_integrate(ias15, 3000.0), sys::MYB_OK);
        let state = std::slice::from_raw_parts(sys::myb_ias15_state(ias15), 12).to_vec();
        sys::myb_ias15_destroy(ias15);
        sys::myb_system_destroy(raw);
        state
    };
    assert_eq!(bits(&integrator.state()), bits(&raw_state));
    // failures: options the library refuses
    let wrong = OrbitOptions { timescale: Ias15Timescale::Jerk, epsilon: -1.0, ..Default::default() };
    assert_eq!(Orbit::new(&model, &wrong).err().unwrap().status, Status::Domain);
    let wrong = Ias15Options { timescale: Ias15Timescale::Jerk, h_max: -1.0, ..Default::default() };
    assert_eq!(SystemIntegrator::new(&system, &wrong).err().unwrap().status, Status::Domain);
    Ok(())
}

// ---- Friction that sticks -------------------------------------------------------------------------

fn law(kind: ContactFriction) -> ContactLaw {
    ContactLaw {
        stiffness: 1.0e4,
        exponent: 1.0,
        damping: 10.0,
        friction: 0.6,
        slip: 0.01,
        kind: kind.raw(),
        tangential_stiffness: 1.0e3,
        tangential_damping: 16.0,
    }
}

fn legs() -> Vec<ContactPoint> {
    [[1.0, 1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, -1.0], [1.0, -1.0, -1.0]]
        .map(|offset| ContactPoint { offset, radius: 0.1 })
        .to_vec()
}

const PLANE: ([f64; 3], [f64; 3]) = ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn raw_plane() -> sys::myb_contact_surface {
    sys::myb_contact_surface { shape: sys::MYB_CONTACT_PLANE, point: PLANE.0, normal: PLANE.1, ..Default::default() }
}

#[test]
fn the_friction_that_sticks_gives_the_bits_of_the_c_calls() -> myb::Result<()> {
    let surface = ContactSurface::Plane { point: PLANE.0, normal: PLANE.1 };
    let points = legs();
    let (x, v, omega) = ([0.0, 0.0, 1.05], [0.1, -0.02, -0.01], [0.0, 0.0, 0.01]);
    let deflection: Vec<f64> = (0..12).map(|k| 1.0e-5 * ((k as f64) * 0.9).sin()).collect();
    let smooth = law(ContactFriction::Smooth);
    let stick = law(ContactFriction::Stick);
    assert_eq!((ContactFriction::of(&stick), ContactFriction::of(&smooth)), (Some(ContactFriction::Stick), Some(ContactFriction::Smooth)));
    let free = myb::contact_force(&surface, &smooth, &points, &x, &v, &IDENTITY, &omega)?;
    let (held, rate) = myb::contact_force_stick(&surface, &stick, &points, &x, &v, &IDENTITY, &omega, &deflection)?;
    let made = raw_plane();
    let (mut raw_free, mut raw_held, mut raw_rate) = (sys::myb_contact_result::default(), sys::myb_contact_result::default(), [0.0; 12]);
    // SAFETY: arrays of the sizes the C interface asks for
    unsafe {
        assert_eq!(
            sys::myb_contact_force(&made, &smooth, points.as_ptr(), 4, x.as_ptr(), v.as_ptr(), IDENTITY.as_ptr(), omega.as_ptr(), &mut raw_free),
            sys::MYB_OK
        );
        assert_eq!(
            sys::myb_contact_force_stick(
                &made,
                &stick,
                points.as_ptr(),
                4,
                x.as_ptr(),
                v.as_ptr(),
                IDENTITY.as_ptr(),
                omega.as_ptr(),
                deflection.as_ptr(),
                raw_rate.as_mut_ptr(),
                &mut raw_held,
            ),
            sys::MYB_OK
        );
    }
    same(&free, &raw_free);
    same(&held, &raw_held);
    assert_eq!(bits(&rate), bits(&raw_rate));
    assert_eq!(held.touching, 4);
    // the smooth friction through the stick form: the same as contact_force, no state
    let (again, _) = myb::contact_force_stick(&surface, &smooth, &points, &x, &v, &IDENTITY, &omega, &[])?;
    same(&again, &free);
    // failures: the friction that sticks without its state; a state of another length
    assert_eq!(myb::contact_force(&surface, &stick, &points, &x, &v, &IDENTITY, &omega).unwrap_err().status, Status::Domain);
    assert_eq!(myb::contact_force_stick(&surface, &stick, &points, &x, &v, &IDENTITY, &omega, &[]).unwrap_err().status, Status::Domain);
    assert_eq!(myb::contact_force_stick(&surface, &stick, &points, &x, &v, &IDENTITY, &omega, &deflection[..6]).unwrap_err().status, Status::Domain);

    // in a system: the deflections are coordinates of the state
    let gm = [0.0, 1.0e-6];
    let system = System::new(&gm)?;
    system.add_contact(1, 0, &SystemContact { surface, law: stick, points: points.clone() })?;
    assert_eq!(system.coordinates(), 6 + 12);
    let mut y = vec![0.0; 36];
    y[3..6].copy_from_slice(&x);
    y[21..24].copy_from_slice(&v);
    system.pack_contact(0, &deflection, &mut y)?;
    assert_eq!(bits(&system.unpack_contact(0, &y)?), bits(&deflection));
    let a = system.acceleration(0.0, &y[..18], &y[18..])?;
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    unsafe {
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_system_create(2, gm.as_ptr(), ptr::null(), &mut raw), sys::MYB_OK);
        let contact = sys::myb_system_contact { surface: made, law: stick, points: points.as_ptr(), point_count: 4 };
        assert_eq!(sys::myb_system_add_contact(raw, 1, 0, &contact), sys::MYB_OK);
        let mut raw_y = vec![0.0; 36];
        raw_y[3..6].copy_from_slice(&x);
        raw_y[21..24].copy_from_slice(&v);
        assert_eq!(sys::myb_system_pack_contact(raw, 0, deflection.as_ptr(), raw_y.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(bits(&y), bits(&raw_y));
        let mut out = [0.0; 12];
        assert_eq!(sys::myb_system_unpack_contact(raw, 0, raw_y.as_ptr(), out.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(bits(&out), bits(&deflection));
        let mut raw_a = vec![0.0; 18];
        assert_eq!(sys::myb_system_acceleration(raw.cast(), 0.0, raw_y.as_ptr(), raw_y[18..].as_ptr(), raw_a.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(bits(&a), bits(&raw_a));
        sys::myb_system_destroy(raw);
    }
    assert!(a[5] > 0.0, "{a:?}"); // pushed up by the springs
    // failures: no such contact; a contact whose friction does not stick
    assert_eq!(system.pack_contact(1, &deflection, &mut y).unwrap_err().status, Status::Domain);
    assert_eq!(system.pack_contact(0, &deflection[..9], &mut y).unwrap_err().status, Status::Domain);
    let slippery = System::new(&gm)?;
    slippery.add_contact(1, 0, &SystemContact { surface: ContactSurface::Plane { point: PLANE.0, normal: PLANE.1 }, law: smooth, points })?;
    let mut y = vec![0.0; 2 * slippery.coordinates()];
    assert_eq!(slippery.pack_contact(0, &deflection, &mut y).unwrap_err().status, Status::Domain);
    assert_eq!(slippery.unpack_contact(0, &y).unwrap_err().status, Status::Domain);
    Ok(())
}

// ---- The new orders of a PN pair ------------------------------------------------------------------

const C: f64 = 299792458.0;
const BH: [f64; 2] = [1.327e21, 0.9e21];

/// A system of two black holes with a PN pair, through the safe layer and through C: the
/// accelerations at the state, or the status of the C calls.
fn pn_pair(orders: u32, speed: f64, spins: bool) -> (myb::Result<Vec<f64>>, Result<Vec<f64>, sys::myb_status>) {
    let pair = SystemPnPair { orders, c: 0.0, r0: 1.0e7, kappa1: 1.0, kappa2: 1.2 };
    let ppn = Ppn { c: C, beta: 1.0, gamma: 1.0 };
    let r = 2.0e7;
    let state = |n: usize| {
        let mut y = vec![0.0; 2 * n];
        y[3] = r;
        y[n + 4] = speed;
        y
    };
    let (chi1, chi2) = ([0.1, -0.2, 0.6], [0.3, 0.1, -0.4]);
    let safe = (|| -> myb::Result<Vec<f64>> {
        let system = System::new(&BH)?;
        system.nbody().set_relativity(Some(&ppn))?;
        system.add_pn_pair(0, 1, &pair)?;
        let n = system.coordinates();
        let mut y = state(n);
        if spins {
            system.pack_spins(0, &chi1, &chi2, &mut y)?;
        }
        system.acceleration(0.0, &y[..n], &y[n..])
    })();
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    let raw = unsafe {
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_system_create(2, BH.as_ptr(), ptr::null(), &mut raw), sys::MYB_OK);
        assert_eq!(sys::myb_nbody_set_relativity(sys::myb_system_nbody(raw), &ppn), sys::MYB_OK);
        let mut outcome = Err(sys::myb_system_add_pn_pair(raw, 0, 1, &pair));
        if outcome == Err(sys::MYB_OK) {
            let n = sys::myb_system_coordinates(raw);
            let mut y = state(n);
            if spins {
                assert_eq!(sys::myb_system_pack_spins(raw, 0, chi1.as_ptr(), chi2.as_ptr(), y.as_mut_ptr()), sys::MYB_OK);
            }
            let mut a = vec![0.0; n];
            let code = sys::myb_system_acceleration(raw.cast(), 0.0, y.as_ptr(), y[n..].as_ptr(), a.as_mut_ptr());
            outcome = if code == sys::MYB_OK { Ok(a) } else { Err(code) };
        }
        sys::myb_system_destroy(raw);
        outcome
    };
    (safe, raw)
}

#[test]
fn the_new_pn_orders_give_the_bits_of_the_c_calls() {
    use myb::{PN_2, PN_4, PN_4RR, PN_SO, PN_SORR5, PN_SOTAIL, PN_SS, PN_SS3, PN_TAILX};
    let circular = ((BH[0] + BH[1]) / 2.0e7_f64).sqrt();
    // the spin terms at the next orders, with the spins integrated
    let (safe, raw) = pn_pair(PN_2 | PN_SO | PN_SS | PN_SS3 | PN_SORR5 | PN_SOTAIL, circular, true);
    let (safe, raw) = (safe.unwrap(), raw.unwrap());
    assert_eq!(safe.len(), 12);
    assert_eq!(bits(&safe), bits(&raw));
    // the tails on an unbound conic
    let (safe, raw) = pn_pair(PN_4 | PN_4RR | PN_TAILX, 2.0 * circular, false);
    assert_eq!(bits(&safe.unwrap()), bits(&raw.unwrap()));
    // failures: the same tails without PN_TAILX; a negative kappa with PN_SS3
    let (safe, raw) = pn_pair(PN_4 | PN_4RR, 2.0 * circular, false);
    assert_eq!(safe.unwrap_err().status.raw(), raw.unwrap_err());
    let system = System::new(&BH).unwrap();
    system.nbody().set_relativity(Some(&Ppn { c: C, beta: 1.0, gamma: 1.0 })).unwrap();
    let negative = SystemPnPair { orders: PN_SO | PN_SS3, c: 0.0, r0: 0.0, kappa1: -1.0, kappa2: 1.0 };
    assert_eq!(system.add_pn_pair(0, 1, &negative).unwrap_err().status, Status::Domain);
    assert_eq!(myb::PN_TAILX, 32768);
}
