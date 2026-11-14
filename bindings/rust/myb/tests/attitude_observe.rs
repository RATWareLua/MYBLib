// SPDX-License-Identifier: BSD-3-Clause
//! Attitude (`myb/attitude.h`) and observables (`myb/observe.h`) through the safe layer
//! against the same C calls, bit for bit; failures as `Err` with the status of the library.
use std::ffi::c_void;
use std::fmt::Debug;
use std::ptr;
use std::rc::Rc;

use myb::sys;
use myb::{FreeRotation, LightBody, LightDirection, LightModel, PlasmaProfile, State, Status, Vec3};

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn same<T: Debug>(a: &T, b: &T) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

#[test]
fn quaternions_and_rotation_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let (a, b) = ([0.9, 0.1, -0.3, 0.2], [0.5, -0.5, 0.5, 0.4]);
    let (v, omega, phi) = ([1.0, 2.0, 3.0], [0.01, -0.02, 0.3], [0.2, -0.1, 0.05]);
    let (normalized, norm) = myb::quat_normalize(&a);
    let m = myb::quat_to_matrix(&a);
    // SAFETY: arrays of the sizes of the C interface
    unsafe {
        let mut out = [0.0; 4];
        sys::myb_quat_multiply(a.as_ptr(), b.as_ptr(), out.as_mut_ptr());
        assert_eq!(bits(&myb::quat_multiply(&a, &b)), bits(&out));
        sys::myb_quat_conjugate(a.as_ptr(), out.as_mut_ptr());
        assert_eq!(bits(&myb::quat_conjugate(&a)), bits(&out));
        let mut q = a;
        let n = sys::myb_quat_normalize(q.as_mut_ptr());
        assert_eq!((bits(&normalized), norm.to_bits()), (bits(&q), n.to_bits()));
        let mut raw_m = [0.0; 9];
        sys::myb_quat_to_matrix(a.as_ptr(), raw_m.as_mut_ptr());
        assert_eq!(bits(&m), bits(&raw_m));
        sys::myb_quat_from_matrix(m.as_ptr(), out.as_mut_ptr());
        assert_eq!(bits(&myb::quat_from_matrix(&m)), bits(&out));
        let mut r = [0.0; 3];
        sys::myb_quat_rotate(a.as_ptr(), v.as_ptr(), r.as_mut_ptr());
        assert_eq!(bits(&myb::quat_rotate(&a, &v)), bits(&r));
        sys::myb_quat_from_rotation_vector(phi.as_ptr(), out.as_mut_ptr());
        assert_eq!(bits(&myb::quat_from_rotation_vector(&phi)), bits(&out));
        sys::myb_quat_to_rotation_vector(normalized.as_ptr(), r.as_mut_ptr());
        assert_eq!(bits(&myb::quat_to_rotation_vector(&normalized)), bits(&r));
        sys::myb_quat_slerp(normalized.as_ptr(), b.as_ptr(), 0.3, out.as_mut_ptr());
        assert_eq!(bits(&myb::quat_slerp(&normalized, &b, 0.3)), bits(&out));
        sys::myb_quat_derivative(a.as_ptr(), omega.as_ptr(), out.as_mut_ptr());
        assert_eq!(bits(&myb::quat_derivative(&a, &omega)), bits(&out));

        let axes = [3, 1, 3];
        let angles = [0.3, 0.4, 0.5];
        sys::myb_euler_to_matrix(axes.as_ptr(), angles.as_ptr(), raw_m.as_mut_ptr());
        let euler = myb::euler_to_matrix(&[3, 1, 3], &angles)?;
        assert_eq!(bits(&euler), bits(&raw_m));
        sys::myb_euler_from_matrix(axes.as_ptr(), raw_m.as_ptr(), r.as_mut_ptr());
        assert_eq!(bits(&myb::euler_from_matrix(&[3, 1, 3], &euler)?), bits(&r));

        let inertia = [10.0, 0.5, 0.1, 0.5, 12.0, -0.2, 0.1, -0.2, 15.0];
        let (moments, axes_out) = myb::principal_axes(&inertia)?;
        let (mut raw_moments, mut raw_axes) = ([0.0; 3], [0.0; 9]);
        assert_eq!(sys::myb_principal_axes(inertia.as_ptr(), raw_moments.as_mut_ptr(), raw_axes.as_mut_ptr()), sys::MYB_OK);
        assert_eq!((bits(&moments), bits(&axes_out)), (bits(&raw_moments), bits(&raw_axes)));

        let principal = [10.0, 12.0, 15.0];
        let l = [1.0, -2.0, 30.0];
        let (dq, dl) = myb::rigid_body_rates(&principal, &normalized, &l, Some(&[0.1, 0.0, -0.1]));
        let (mut raw_dq, mut raw_dl) = ([0.0; 4], [0.0; 3]);
        sys::myb_rigid_body_rates(principal.as_ptr(), normalized.as_ptr(), l.as_ptr(), [0.1, 0.0, -0.1].as_ptr(), raw_dq.as_mut_ptr(), raw_dl.as_mut_ptr());
        assert_eq!((bits(&dq), bits(&dl)), (bits(&raw_dq), bits(&raw_dl)));

        let core = [0.5, 0.6];
        let (wm, wc) = ([1e-6, 2e-6, 2.66e-6], [0.0, 0.0, 2.6e-6]);
        let (rm, rc) = myb::core_mantle_rates(&inertia, None, &core, 1e-8, &wm, &wc, &[0.0; 3])?;
        let (mut raw_rm, mut raw_rc) = ([0.0; 3], [0.0; 3]);
        assert_eq!(
            sys::myb_core_mantle_rates(inertia.as_ptr(), ptr::null(), core.as_ptr(), 1e-8, wm.as_ptr(), wc.as_ptr(), [0.0; 3].as_ptr(), raw_rm.as_mut_ptr(), raw_rc.as_mut_ptr()),
            sys::MYB_OK
        );
        assert_eq!((bits(&rm), bits(&rc)), (bits(&raw_rm), bits(&raw_rc)));

        let rel = [3.8e8, 1.0e7, -2.0e6];
        sys::myb_torque_point_mass(4.9e12, rel.as_ptr(), principal.as_ptr(), r.as_mut_ptr());
        assert_eq!(bits(&myb::torque_point_mass(4.9e12, &rel, &principal)), bits(&r));
        sys::myb_torque_field(0.0123, rel.as_ptr(), [1e-3, 2e-3, -1e-3].as_ptr(), r.as_mut_ptr());
        assert_eq!(bits(&myb::torque_field(0.0123, &rel, &[1e-3, 2e-3, -1e-3])), bits(&r));
        let gq = [1e9, 0.0, 0.0, 0.0, 2e9, 0.0, 0.0, 0.0, -3e9];
        assert_eq!(sys::myb_torque_figure_figure(3.986e14, gq.as_ptr(), gq.as_ptr(), rel.as_ptr(), r.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(bits(&myb::torque_figure_figure(3.986e14, &gq, &gq, &rel)?), bits(&r));

        // the exact free rotation
        let free = FreeRotation::new(&principal, &normalized, &l)?;
        let (q_t, l_t) = free.state(123.0);
        let mut raw_free = sys::myb_free_rotation::default();
        assert_eq!(sys::myb_free_rotation_init(principal.as_ptr(), normalized.as_ptr(), l.as_ptr(), &mut raw_free), sys::MYB_OK);
        let (mut raw_q, mut raw_l) = ([0.0; 4], [0.0; 3]);
        sys::myb_free_rotation_state(&raw_free, 123.0, raw_q.as_mut_ptr(), raw_l.as_mut_ptr());
        assert_eq!((bits(&q_t), bits(&l_t)), (bits(&raw_q), bits(&raw_l)));
    }

    // the splitting integrator with the torque of a closure, against a C function
    unsafe extern "C" fn raw_torque(_c: *mut c_void, t: f64, q: *const f64, n: *mut f64) -> sys::myb_status {
        // SAFETY: q holds four numbers, n three
        unsafe { *(n as *mut Vec3) = [1e-3 * *q.add(1), -1e-3 * *q, 1e-4 * t.sin()] };
        sys::MYB_OK
    }
    let principal = [10.0, 12.0, 15.0];
    let l = [1.0, -2.0, 30.0];
    let torque = |t: f64, q: &myb::Quat| -> myb::Result<Vec3> { Ok([1e-3 * q[1], -1e-3 * q[0], 1e-4 * t.sin()]) };
    let (q, l_end) = myb::rotation_split(&principal, torque, 4, 0.0, 0.5, 200, &normalized, &l)?;
    let (mut raw_q, mut raw_l) = (normalized, l);
    // SAFETY: arrays of four and three numbers
    assert_eq!(
        unsafe { sys::myb_rotation_split(principal.as_ptr(), Some(raw_torque), ptr::null_mut(), 4, 0.0, 0.5, 200, raw_q.as_mut_ptr(), raw_l.as_mut_ptr()) },
        sys::MYB_OK
    );
    assert_eq!((bits(&q), bits(&l_end)), (bits(&raw_q), bits(&raw_l)));

    // failures
    assert_eq!(myb::euler_to_matrix(&[3, 3, 1], &[0.0; 3]).unwrap_err().status, Status::Domain);
    assert_eq!(myb::principal_axes(&[1.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap_err().status, Status::Domain);
    assert_eq!(FreeRotation::new(&[0.0, 1.0, 1.0], &normalized, &l).err().unwrap().status, Status::Domain);
    assert_eq!(myb::rotation_split(&principal, torque, 3, 0.0, 0.5, 10, &normalized, &l).unwrap_err().status, Status::Domain);
    let failing = |_: f64, _: &myb::Quat| -> myb::Result<Vec3> { Err(myb::Error::new(Status::Range)) };
    assert_eq!(myb::rotation_split(&principal, failing, 2, 0.0, 0.5, 10, &normalized, &l).unwrap_err().status, Status::Range);
    assert_eq!(myb::torque_figure_figure(3.986e14, &[0.0; 9], &[0.0; 9], &[0.0; 3]).unwrap_err().status, Status::Domain);
    let _ = m;
    Ok(())
}

// ---- Observables ----------------------------------------------------------------------------------

const C: f64 = 299792458.0;

fn target_at(seconds: f64) -> State {
    [1.5e11 + 3.0e4 * seconds, 2.0e10, -1.0e9, 3.0e4, 0.0, 0.0]
}

fn sun_at(_seconds: f64) -> State {
    [0.0; 6]
}

fn station_at(seconds: f64) -> State {
    let w = 7.29e-5;
    [6.4e6 * (w * seconds).cos(), 6.4e6 * (w * seconds).sin(), 0.0, -6.4e6 * w * (w * seconds).sin(), 6.4e6 * w * (w * seconds).cos(), 0.0]
}

unsafe extern "C" fn raw_target(_c: *mut c_void, t: sys::myb_time, y: *mut f64) -> sys::myb_status {
    // SAFETY: y holds six numbers
    unsafe { *(y as *mut State) = target_at(sys::myb_time_seconds(t)) };
    sys::MYB_OK
}

unsafe extern "C" fn raw_sun(_c: *mut c_void, t: sys::myb_time, y: *mut f64) -> sys::myb_status {
    // SAFETY: y holds six numbers
    unsafe { *(y as *mut State) = sun_at(sys::myb_time_seconds(t)) };
    sys::MYB_OK
}

unsafe extern "C" fn raw_station(_c: *mut c_void, t: sys::myb_time, y: *mut f64) -> sys::myb_status {
    // SAFETY: y holds six numbers
    unsafe { *(y as *mut State) = station_at(sys::myb_time_seconds(t)) };
    sys::MYB_OK
}

#[test]
fn observables_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let sun = LightBody { gm: 1.32712440018e20, state: Rc::new(|t: myb::Time| -> myb::Result<State> { Ok(sun_at(myb::seconds(t))) }) };
    let model = LightModel { c: C, gamma: 1.0, bodies: vec![sun] };
    let target = |t: myb::Time| -> myb::Result<State> { Ok(target_at(myb::seconds(t))) };
    let station = |t: myb::Time| -> myb::Result<State> { Ok(station_at(myb::seconds(t))) };
    let t = myb::time_from_seconds(1000.0);
    let observer = station_at(1000.0);
    let path = model.light_time(&target, t, &observer, LightDirection::Reception)?;
    let forward = model.light_time(&target, t, &observer, LightDirection::Transmission)?;
    let trip = model.round_trip(&station, &target, t)?;
    let bodies = [sys::myb_light_body { gm: 1.32712440018e20, state: Some(raw_sun), context: ptr::null_mut() }];
    let raw_model = sys::myb_light_model { c: C, gamma: 1.0, bodies: bodies.as_ptr(), count: 1 };
    let (mut raw_path, mut raw_forward) = (sys::myb_light_path::default(), sys::myb_light_path::default());
    let (mut total, mut bounce, mut transmit) = (0.0, sys::myb_time::default(), sys::myb_time::default());
    // SAFETY: the providers and the model live for the calls
    unsafe {
        assert_eq!(sys::myb_light_time(&raw_model, Some(raw_target), ptr::null_mut(), t, observer.as_ptr(), sys::MYB_LIGHT_RECEPTION, &mut raw_path), sys::MYB_OK);
        assert_eq!(sys::myb_light_time(&raw_model, Some(raw_target), ptr::null_mut(), t, observer.as_ptr(), sys::MYB_LIGHT_TRANSMISSION, &mut raw_forward), sys::MYB_OK);
        assert_eq!(
            sys::myb_round_trip(&raw_model, Some(raw_station), ptr::null_mut(), Some(raw_target), ptr::null_mut(), t, &mut total, &mut bounce, &mut transmit),
            sys::MYB_OK
        );
    }
    same(&path, &raw_path);
    same(&forward, &raw_forward);
    assert_eq!(trip.total.to_bits(), total.to_bits());
    same(&(trip.t_bounce, trip.t_transmit), &(bounce, transmit));
    assert!(path.delay > 0.0 && path.light_time > 500.0);

    // directions and delays
    let (x1, x2, xb) = ([1.5e11, 0.0, 0.0], [-2.0e11, 1.0e9, 0.0], [0.0; 3]);
    let p = [0.6, 0.8, 0.0];
    let profile = PlasmaProfile { radius: 6.96e8, count: 2, density: [1e14, 1e12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], power: [6.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0] };
    let column = myb::plasma_column(&profile, &x1, &x2)?;
    // SAFETY: arrays of the sizes of the C interface
    unsafe {
        assert_eq!(myb::shapiro_delay(1.327e20, C, 1.0, &x1, &x2, &xb).to_bits(), sys::myb_shapiro_delay(1.327e20, C, 1.0, x1.as_ptr(), x2.as_ptr(), xb.as_ptr()).to_bits());
        let mut out = [0.0; 3];
        sys::myb_aberration(p.as_ptr(), [3.0e4, 0.0, 1.0e3].as_ptr(), C, out.as_mut_ptr());
        assert_eq!(bits(&myb::aberration(&p, &[3.0e4, 0.0, 1.0e3], C)), bits(&out));
        let (q, e) = ([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
        sys::myb_light_deflection(1.327e20, C, 1.0, p.as_ptr(), q.as_ptr(), e.as_ptr(), 1.5e11, 1e-6, out.as_mut_ptr());
        assert_eq!(bits(&myb::light_deflection(1.327e20, C, 1.0, &p, &q, &e, 1.5e11, 1e-6)), bits(&out));
        let (mut ra, mut dec) = (0.0, 0.0);
        sys::myb_direction_radec([1.0, -2.0, 0.5].as_ptr(), &mut ra, &mut dec);
        let (r, d) = myb::direction_radec(&[1.0, -2.0, 0.5]);
        assert_eq!((r.to_bits(), d.to_bits()), (ra.to_bits(), dec.to_bits()));
        let e_charge = 1.602176634e-19;
        assert_eq!(
            myb::plasma_group_delay(1e17, 2.2e9, e_charge, 9.1093837139e-31, 8.8541878188e-12).to_bits(),
            sys::myb_plasma_group_delay(1e17, 2.2e9, e_charge, 9.1093837139e-31, 8.8541878188e-12).to_bits()
        );
        let mut raw_column = 0.0;
        assert_eq!(sys::myb_plasma_column(&profile, x1.as_ptr(), x2.as_ptr(), &mut raw_column), sys::MYB_OK);
        assert_eq!(column.to_bits(), raw_column.to_bits());
    }

    // failures
    let invalid = LightModel { c: 0.0, gamma: 1.0, bodies: Vec::new() };
    assert_eq!(invalid.light_time(&target, t, &observer, LightDirection::Reception).unwrap_err().status, Status::Domain);
    let lost = |_t: myb::Time| -> myb::Result<State> { Err(myb::Error::new(Status::Range)) };
    assert_eq!(model.light_time(&lost, t, &observer, LightDirection::Reception).unwrap_err().status, Status::Range);
    assert_eq!(model.round_trip(&lost, &target, t).unwrap_err().status, Status::Range);
    assert_eq!(myb::plasma_column(&profile, &[0.0; 3], &x2).unwrap_err().status, Status::Domain);
    let empty = PlasmaProfile { count: 0, ..profile };
    assert_eq!(myb::plasma_column(&empty, &x1, &x2).unwrap_err().status, Status::Domain);
    Ok(())
}
