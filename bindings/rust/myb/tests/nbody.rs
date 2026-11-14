// SPDX-License-Identifier: BSD-3-Clause
//! N bodies (`myb/dyn.h`, `myb/system.h`) through the safe types against the same C calls
//! through `myb::sys`, bit for bit; failures as `Err` with the status of the library.
use std::ffi::c_void;
use std::fmt::Debug;
use std::ptr;
use std::rc::Rc;

use myb::sys;
use myb::{
    FixedOrientation, GravityField, Ias15Options, Interaction, NBody, Pool, Ppn, Rotation, Status, System, SystemEventKind,
    SystemIntegrator, SystemRotation, Wh, WhOptions, WhScheme,
};

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

/// Two structures of the library are the same when their exact printed forms are (the Debug
/// form of a float reads back to its bits).
fn same<T: Debug>(a: &T, b: &T) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

const GM: [f64; 3] = [1.32712440018e20, 1.26686534e17, 3.7931187e16];

fn bodies() -> (Vec<f64>, Vec<f64>) {
    let x = vec![0.0, 0.0, 0.0, 7.78e11, 0.0, 1.0e9, -1.0e11, 1.42e12, -2.0e10, 5.0e11, 2.0e10, 0.0];
    let v = vec![0.0, 0.0, 0.0, 0.0, 13.07e3, 10.0, -9.68e3, -0.68e3, 0.0, -1.0e3, 16.29e3, 100.0];
    (x, v)
}

fn raw_nbody(test: usize) -> *mut sys::myb_nbody {
    let mut raw = ptr::null_mut();
    // SAFETY: GM holds three numbers
    assert_eq!(unsafe { sys::myb_nbody_create(3, test, GM.as_ptr(), ptr::null(), &mut raw) }, sys::MYB_OK);
    raw
}

#[test]
fn point_masses_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let (x, v) = bodies();
    let nbody = NBody::new(&GM, 1)?;
    assert_eq!((nbody.count(), nbody.massive(), nbody.gm(3), nbody.gm(4)), (4, 3, Some(0.0), None));
    let raw = raw_nbody(1);
    let mut a = vec![0.0; 12];
    let mut invariants = sys::myb_nbody_invariants::default();
    // SAFETY: the arrays hold 3 numbers per body
    unsafe {
        assert_eq!(sys::myb_nbody_acceleration(raw.cast(), 0.0, x.as_ptr(), a.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_nbody_invariants_of(raw, x.as_ptr(), v.as_ptr(), &mut invariants), sys::MYB_OK);
    }
    assert_eq!(bits(&nbody.acceleration(0.0, &x)?), bits(&a));
    same(&nbody.invariants(&x, &v)?, &invariants);

    // threads give the same bits
    let pool = Rc::new(Pool::new(2)?);
    assert_eq!(pool.threads(), 2);
    nbody.set_pool(Some(&pool));
    assert_eq!(bits(&nbody.acceleration(0.0, &x)?), bits(&a));
    nbody.set_pool(None);

    // the graph: the test body does not feel Saturn
    nbody.set_interaction(3, 2, Interaction::None)?;
    // SAFETY: the indices are valid
    assert_eq!(unsafe { sys::myb_nbody_set_interaction(raw, 3, 2, sys::MYB_INTERACTION_NONE) }, sys::MYB_OK);
    unsafe { sys::myb_nbody_acceleration(raw.cast(), 0.0, x.as_ptr(), a.as_mut_ptr()) };
    assert_eq!(bits(&nbody.acceleration(0.0, &x)?), bits(&a));

    // the approximation
    nbody.set_approximation(8, 0.5)?;
    assert_eq!(nbody.approximation(), (8, 0.5));
    // SAFETY: valid settings
    unsafe {
        assert_eq!(sys::myb_nbody_set_approximation(raw, 8, 0.5), sys::MYB_OK);
        assert_eq!(sys::myb_nbody_acceleration(raw.cast(), 0.0, x.as_ptr(), a.as_mut_ptr()), sys::MYB_OK);
    }
    assert_eq!(bits(&nbody.acceleration(0.0, &x)?), bits(&a));
    nbody.set_approximation(0, 0.0)?;
    unsafe { sys::myb_nbody_set_approximation(raw, 0, 0.0) };

    // the 1PN terms
    let ppn = Ppn { c: 299792458.0, beta: 1.0, gamma: 1.0 };
    nbody.set_relativity(Some(&ppn))?;
    let mut pn = sys::myb_nbody_pn_invariants::default();
    // SAFETY: valid parameters and arrays
    unsafe {
        assert_eq!(sys::myb_nbody_set_relativity(raw, &ppn), sys::MYB_OK);
        assert_eq!(sys::myb_nbody_acceleration_xv(raw.cast(), 0.0, x.as_ptr(), v.as_ptr(), a.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_nbody_pn_invariants_of(raw, x.as_ptr(), v.as_ptr(), &mut pn), sys::MYB_OK);
    }
    assert_eq!(bits(&nbody.acceleration_xv(0.0, &x, &v)?), bits(&a));
    same(&nbody.pn_invariants(&x, &v)?, &pn);
    assert_eq!(nbody.acceleration(0.0, &x).unwrap_err().status, Status::Domain);
    unsafe { sys::myb_nbody_destroy(raw) };

    // failures
    assert_eq!(NBody::new(&[], 0).err().unwrap().status, Status::Domain);
    assert_eq!(NBody::new(&[-1.0], 0).err().unwrap().status, Status::Domain);
    assert_eq!(nbody.set_interaction(1, 1, Interaction::Newton).unwrap_err().status, Status::Domain);
    assert_eq!(nbody.set_approximation(21, 0.5).unwrap_err().status, Status::Domain);
    assert_eq!(nbody.acceleration_xv(0.0, &x[..9], &v[..9]).unwrap_err().status, Status::Domain);
    let coincident = vec![0.0; 12];
    assert_eq!(nbody.invariants(&coincident, &v).unwrap_err().status, Status::Domain);
    Ok(())
}

#[test]
fn wisdom_holman_gives_the_bits_of_the_c_calls() -> myb::Result<()> {
    let (x, v) = bodies();
    let options = WhOptions { h: 86400.0 * 10.0, scheme: WhScheme::Dkd, corrector: 5 };
    let nbody = Rc::new(NBody::new(&GM, 1)?);
    let mut wh = Wh::new(&nbody, &options)?;
    drop(nbody); // the integrator keeps its system
    wh.start(0.0, &x, &v)?;
    wh.integrate(86400.0 * 3650.0)?;
    let (xs, vs) = wh.state()?;
    let t_dense = wh.time() + 86400.0 * 3.0;
    let (xd, vd) = wh.dense(t_dense)?;

    let raw = raw_nbody(1);
    let made = sys::myb_wh_options { h: options.h, scheme: sys::MYB_WH_DKD, corrector: 5 };
    let (mut xr, mut vr, mut xdr, mut vdr) = (vec![0.0; 12], vec![0.0; 12], vec![0.0; 12], vec![0.0; 12]);
    let mut integrator = ptr::null_mut();
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    let stats = unsafe {
        assert_eq!(sys::myb_wh_create(raw, &made, ptr::null(), &mut integrator), sys::MYB_OK);
        assert_eq!(sys::myb_wh_start(integrator, 0.0, x.as_ptr(), v.as_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_wh_integrate(integrator, 86400.0 * 3650.0), sys::MYB_OK);
        assert_eq!(sys::myb_wh_state(integrator, xr.as_mut_ptr(), vr.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_wh_dense(integrator, t_dense, xdr.as_mut_ptr(), vdr.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_wh_time(integrator).to_bits(), wh.time().to_bits());
        let stats = sys::myb_wh_statistics(integrator);
        sys::myb_wh_destroy(integrator);
        sys::myb_nbody_destroy(raw);
        stats
    };
    assert_eq!((bits(&xs), bits(&vs)), (bits(&xr), bits(&vr)));
    assert_eq!((bits(&xd), bits(&vd)), (bits(&xdr), bits(&vdr)));
    same(&wh.statistics(), &stats);
    wh.step()?;

    // failures: no step; an approximation in the system
    let nbody = Rc::new(NBody::new(&GM, 1)?);
    assert_eq!(Wh::new(&nbody, &WhOptions::default()).err().unwrap().status, Status::Domain);
    nbody.set_approximation(4, 0.5)?;
    let mut refused = Wh::new(&nbody, &options)?;
    assert_eq!(refused.start(0.0, &x, &v).unwrap_err().status, Status::Domain);
    assert_eq!(refused.start(0.0, &x[..3], &v).unwrap_err().status, Status::Domain);
    Ok(())
}

// ---- A system of bodies ---------------------------------------------------------------------------

const EARTH: f64 = 3.986004418e14;
const RADIUS: f64 = 6378137.0;
const PUSH: [f64; 3] = [1.0e-9, -2.0e-9, 0.5e-9];
const FAR: [f64; 6] = [3.8e8, 1.0e8, 0.0, 0.0, 0.0, 0.0];

fn identity() -> Rotation {
    Rotation { m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], dm: [0.0; 9] }
}

fn j2() -> ([f64; 6], [f64; 6]) {
    ([1.0, 0.0, 0.0, -1.08263e-3 / 5.0_f64.sqrt(), 0.0, 0.0], [0.0; 6])
}

fn system_state() -> Vec<f64> {
    let r = RADIUS + 500.0e3;
    let speed = (EARTH / r).sqrt();
    let mut y = vec![0.0; 18];
    // positions: the Earth, the Moon, the satellite; then the velocities
    y[3..6].copy_from_slice(&[3.84e8, 0.0, 0.0]);
    y[6..9].copy_from_slice(&[r, 0.0, 0.0]);
    y[13] = 1.022e3;
    y[16] = speed * 0.8;
    y[17] = speed * 0.6;
    y
}

unsafe extern "C" fn push(_context: *mut c_void, _t: f64, _view: *const sys::myb_system_view, a: *mut f64, _n: *mut f64) -> sys::myb_status {
    // SAFETY: a holds three numbers
    unsafe { *(a as *mut [f64; 3]) = PUSH };
    sys::MYB_OK
}

unsafe extern "C" fn far(_context: *mut c_void, _t: sys::myb_time, y: *mut f64) -> sys::myb_status {
    // SAFETY: y holds six numbers
    unsafe { *(y as *mut [f64; 6]) = FAR };
    sys::MYB_OK
}

/// The integration and the derivatives of the system through the C calls.
fn system_raw(y0: &[f64], t1: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (c, s) = j2();
    let rotation = identity();
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    unsafe {
        let mut field = ptr::null_mut();
        assert_eq!(sys::myb_gravity_field_create(EARTH, RADIUS, 2, c.as_ptr(), s.as_ptr(), ptr::null(), &mut field), sys::MYB_OK);
        let gm = [EARTH, 4.9028e12, 1.0e3];
        let mut system = ptr::null_mut();
        assert_eq!(sys::myb_system_create(3, gm.as_ptr(), ptr::null(), &mut system), sys::MYB_OK);
        let context = &rotation as *const Rotation as *mut c_void;
        assert_eq!(sys::myb_system_set_orientation(system, 0, Some(sys::myb_orientation_fixed), context), sys::MYB_OK);
        assert_eq!(sys::myb_system_set_field(system, 0, field, 0), sys::MYB_OK);
        assert_eq!(sys::myb_system_add_field(system, 0, 2), sys::MYB_OK);
        assert_eq!(sys::myb_system_add_action(system, 2, Some(push), ptr::null_mut()), sys::MYB_OK);
        assert_eq!(sys::myb_system_add_guided(system, 1.0e10, Some(far), ptr::null_mut(), ptr::null_mut()), sys::MYB_OK);
        let mut integrator = ptr::null_mut();
        assert_eq!(sys::myb_system_integrator(system, ptr::null(), ptr::null(), &mut integrator), sys::MYB_OK);
        assert_eq!(sys::myb_ias15_start(integrator, Some(sys::myb_system_acceleration), system.cast(), 0.0, y0.as_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_ias15_integrate(integrator, t1), sys::MYB_OK);
        let state = std::slice::from_raw_parts(sys::myb_ias15_state(integrator), 18).to_vec();
        let options = sys::myb_ias15_options::default();
        let (mut y1, mut phi, mut sensitivity) = (vec![0.0; 18], vec![0.0; 324], vec![0.0; 18]);
        let chosen = [1usize];
        assert_eq!(
            sys::myb_system_partials(system, &options, 0.0, y0.as_ptr(), 600.0, chosen.as_ptr(), 1, y1.as_mut_ptr(), phi.as_mut_ptr(), sensitivity.as_mut_ptr()),
            sys::MYB_OK
        );
        sys::myb_ias15_destroy(integrator);
        sys::myb_system_destroy(system);
        sys::myb_gravity_field_destroy(field);
        phi.extend_from_slice(&y1);
        (state, phi, sensitivity)
    }
}

#[test]
fn a_system_of_bodies_gives_the_bits_of_the_c_calls() -> myb::Result<()> {
    let (c, s) = j2();
    let field = Rc::new(GravityField::new(EARTH, RADIUS, 2, &c, &s)?);
    let system = System::new(&[EARTH, 4.9028e12, 1.0e3])?;
    assert_eq!((system.bodies(), system.coordinates()), (3, 9));
    system.set_orientation(0, FixedOrientation(identity()))?;
    system.set_field(0, &field, None)?;
    drop(field); // the system keeps it
    system.add_field(0, 2)?;
    let seen = Rc::new(std::cell::Cell::new((0, 0, 0, 0)));
    let record = Rc::clone(&seen);
    system.add_action(2, move |_t, view| {
        record.set((view.body, view.bodies, view.x.len(), view.gm.len()));
        Ok((PUSH, [0.0; 3]))
    })?;
    let guided = system.add_guided(1.0e10, |_t: myb::Time| -> myb::Result<myb::State> { Ok(FAR) })?;
    assert_eq!((guided, system.guided()), (0, 1));
    let system = Rc::new(system);
    let y0 = system_state();
    let mut integrator = SystemIntegrator::new(&system, &Ias15Options::default())?;
    integrator.start(0.0, &y0)?;
    integrator.integrate(3600.0)?;
    let partials = system.partials(&Ias15Options::default(), 0.0, &y0, 600.0, &[1])?;
    let (state, mut phi_raw, sensitivity) = system_raw(&y0, 3600.0);
    let y1_raw = phi_raw.split_off(324);
    assert_eq!(bits(&integrator.state()), bits(&state));
    assert_eq!(integrator.time(), 3600.0);
    assert_eq!(seen.get(), (2, 3, 9, 3));
    assert_eq!(bits(&partials.y1), bits(&y1_raw));
    assert_eq!(bits(&partials.phi), bits(&phi_raw));
    assert_eq!(bits(&partials.sensitivity), bits(&sensitivity));
    let a = system.acceleration(0.0, &y0[..9], &y0[9..])?;
    assert!(a[6] < 0.0 && a.iter().all(|x| x.is_finite()));
    assert!(integrator.statistics().accepted > 0);
    let last = integrator.state();
    let dense = integrator.dense(integrator.time())?;
    assert_eq!(bits(&dense), bits(&last));

    // failures
    assert_eq!(system.add_field(0, 0).unwrap_err().status, Status::Domain);
    assert_eq!(system.add_figures(0, 7).unwrap_err().status, Status::Domain);
    let mut y = vec![0.0; 18];
    assert_eq!(system.pack_rotation(0, &[1.0, 0.0, 0.0, 0.0], &[0.0; 3], None, &mut y).unwrap_err().status, Status::Domain);
    assert_eq!(system.pack_rotation(0, &[1.0, 0.0, 0.0, 0.0], &[0.0; 3], None, &mut y[..3]).unwrap_err().status, Status::Domain);
    assert_eq!(integrator.start(0.0, &y0[..9]).unwrap_err().status, Status::Domain);
    assert_eq!(system.partials(&Ias15Options::default(), 0.0, &y0, 0.0, &[]).unwrap_err().status, Status::Domain);
    assert_eq!((system.radius(9), system.present(9)), (None, false));
    // an action that fails stops the evaluation with its status
    let broken = System::new(&[EARTH, 1.0])?;
    broken.add_action(1, |_, _| Err(myb::Error::new(Status::Range)))?;
    let y = [0.0, 0.0, 0.0, 7.0e6, 0.0, 0.0];
    assert_eq!(broken.acceleration(0.0, &y, &[0.0; 6]).unwrap_err().status, Status::Range);
    Ok(())
}

#[test]
fn rotations_and_encounters() -> myb::Result<()> {
    // an integrated rotation packed and unpacked as the C calls do it
    let rotation = SystemRotation { moments: [0.32, 0.33, 0.34].map(|k| k * RADIUS * RADIUS), mean_motion: 7.29e-5, ..Default::default() };
    let system = System::new(&[EARTH, 4.9028e12])?;
    system.set_rotation(0, &rotation)?;
    assert_eq!(system.coordinates(), 15);
    let q = [0.9, 0.1, -0.3, 0.2];
    let omega = [1.0e-6, -2.0e-6, 7.29e-5];
    let mut y = vec![0.0; 30];
    system.pack_rotation(0, &q, &omega, None, &mut y)?;
    let (q_out, omega_out, _) = system.unpack_rotation(0, &y)?;
    let mut raw_y = vec![0.0; 30];
    let (mut raw_q, mut raw_omega, mut raw_core) = ([0.0; 4], [0.0; 3], [0.0; 3]);
    // SAFETY: the calls of the C interface with valid arguments
    unsafe {
        let mut raw = ptr::null_mut();
        let gm = [EARTH, 4.9028e12];
        assert_eq!(sys::myb_system_create(2, gm.as_ptr(), ptr::null(), &mut raw), sys::MYB_OK);
        assert_eq!(sys::myb_system_set_rotation(raw, 0, &rotation), sys::MYB_OK);
        assert_eq!(sys::myb_system_pack_rotation(raw, 0, q.as_ptr(), omega.as_ptr(), ptr::null(), raw_y.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(
            sys::myb_system_unpack_rotation(raw, 0, raw_y.as_ptr(), raw_q.as_mut_ptr(), raw_omega.as_mut_ptr(), raw_core.as_mut_ptr()),
            sys::MYB_OK
        );
        sys::myb_system_destroy(raw);
    }
    assert_eq!(bits(&y), bits(&raw_y));
    assert_eq!((bits(&q_out), bits(&omega_out)), (bits(&raw_q), bits(&raw_omega)));
    assert_eq!(system.pack_spins(0, &[0.0; 3], &[0.0; 3], &mut y).unwrap_err().status, Status::Domain);

    // two spheres on a collision course: the contact, then a merger
    let gm = [1.0e10, 1.0e9];
    let y0 = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0];
    let system = Rc::new(System::new(&gm)?);
    system.set_radius(0, 1000.0)?;
    system.set_radius(1, 1000.0)?;
    assert_eq!(system.radius(1), Some(1000.0));
    let mut integrator = SystemIntegrator::new(&system, &Ias15Options::default())?;
    integrator.start(0.0, &y0)?;
    let event = integrator.advance(1.0e5, 0.0)?;
    assert_eq!((event.kind, event.a, event.b), (SystemEventKind::Contact, 0, 1));
    let merger = integrator.merge(0, 1)?;
    assert!(!system.present(1) && system.present(0));
    let (raw_event, raw_merger) = unsafe {
        // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_system_create(2, gm.as_ptr(), ptr::null(), &mut raw), sys::MYB_OK);
        assert_eq!(sys::myb_system_set_radius(raw, 0, 1000.0), sys::MYB_OK);
        assert_eq!(sys::myb_system_set_radius(raw, 1, 1000.0), sys::MYB_OK);
        let mut ias15 = ptr::null_mut();
        assert_eq!(sys::myb_system_integrator(raw, ptr::null(), ptr::null(), &mut ias15), sys::MYB_OK);
        assert_eq!(sys::myb_ias15_start(ias15, Some(sys::myb_system_acceleration), raw.cast(), 0.0, y0.as_ptr()), sys::MYB_OK);
        let mut event = sys::myb_system_event::default();
        assert_eq!(sys::myb_system_advance(raw, ias15, 1.0e5, 0.0, &mut event), sys::MYB_OK);
        let mut merger = sys::myb_system_merger::default();
        assert_eq!(sys::myb_system_merge(raw, ias15, 0, 1, &mut merger), sys::MYB_OK);
        sys::myb_ias15_destroy(ias15);
        sys::myb_system_destroy(raw);
        (event, merger)
    };
    assert_eq!((event.t.to_bits(), event.distance.to_bits()), (raw_event.t.to_bits(), raw_event.distance.to_bits()));
    same(&merger, &raw_merger);
    assert_eq!(integrator.merge(0, 1).unwrap_err().status, Status::Domain);

    // a bounce instead
    let system = Rc::new(System::new(&gm)?);
    system.set_radius(0, 1000.0)?;
    system.set_radius(1, 1000.0)?;
    let mut integrator = SystemIntegrator::new(&system, &Ias15Options::default())?;
    integrator.start(0.0, &y0)?;
    integrator.advance(1.0e5, 0.0)?;
    integrator.bounce(0, 1, 1.0)?;
    let after = integrator.state();
    assert!(after[9] > 0.0, "{after:?}");
    assert_eq!(integrator.bounce(0, 1, 2.0).unwrap_err().status, Status::Domain);
    assert_eq!(integrator.advance(2.0e5, -1.0).unwrap_err().status, Status::Domain);
    Ok(())
}
