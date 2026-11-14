// SPDX-License-Identifier: BSD-3-Clause
//! Targeting (`myb/transfer.h`, `myb/lambert.h`, `myb/iod.h`, `myb/relative.h`,
//! `myb/conjunction.h`, `myb/uncert.h`) through the safe layer against the same C calls, bit
//! for bit; failures as `Err` with the status of the library.
use std::ffi::c_void;
use std::fmt::Debug;
use std::ptr;
use std::rc::Rc;

use myb::sys;
use myb::{Error, ForceModel, Orbit, OrbitOptions, Rng, State, Status, Variational, VariationalOptions};

const GM: f64 = 3.986004418e14;

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn same<T: Debug>(a: &T, b: &T) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

fn pv(y: &State) -> sys::myb_pv {
    sys::myb_pv { r: [y[0], y[1], y[2]], v: [y[3], y[4], y[5]] }
}

fn orbit(i: f64, nu: f64) -> State {
    let elements = myb::Elements { p: 7.0e6, e: 0.0, i, raan: 0.0, argp: 0.0, nu };
    myb::state_from_elements(GM, &elements).unwrap()
}

#[test]
fn transfers_and_first_orbits_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let (r, v, pole) = ([-1.0e9, 2.0e8, 1.0e7], [3.0e3, -1.0e2, 50.0], [0.0, 0.0, 1.0]);
    let bplane = myb::bplane(GM, &r, &v, &pole)?;
    let mut raw = sys::myb_bplane::default();
    // SAFETY: arrays of three numbers; a valid place for the output
    unsafe {
        assert_eq!(sys::myb_bplane_compute(GM, r.as_ptr(), v.as_ptr(), pole.as_ptr(), &mut raw), sys::MYB_OK);
        assert_eq!(myb::sphere_of_influence(1.5e11, 3e-6).to_bits(), sys::myb_sphere_of_influence(1.5e11, 3e-6).to_bits());
        assert_eq!(myb::hill_radius(1.5e11, 0.0167, 3e-6).to_bits(), sys::myb_hill_radius(1.5e11, 0.0167, 3e-6).to_bits());
        assert_eq!(myb::flyby_turn(GM, 3e3, 7e6).to_bits(), sys::myb_flyby_turn(GM, 3e3, 7e6).to_bits());
        assert_eq!(myb::departure_impulse(GM, 7e6, 3e3).to_bits(), sys::myb_departure_impulse(GM, 7e6, 3e3).to_bits());
    }
    same(&bplane, &raw);
    assert_eq!(myb::bplane(GM, &[7.0e6, 0.0, 0.0], &[0.0, 7.5e3, 0.0], &pole).unwrap_err().status, Status::Domain);

    // Lambert
    let (r1, r2) = ([7.0e6, 0.0, 0.0], [0.0, 8.0e6, 1.0e5]);
    let solutions = myb::lambert(GM, &r1, &r2, 9000.0, false, 1)?;
    let mut raw_solutions = [sys::myb_lambert_solution::default(); 3];
    let mut count = 0usize;
    // SAFETY: the capacity is that of the array
    assert_eq!(
        unsafe { sys::myb_lambert(GM, r1.as_ptr(), r2.as_ptr(), 9000.0, 0, 1, raw_solutions.as_mut_ptr(), 3, &mut count) },
        sys::MYB_OK
    );
    assert_eq!(solutions.len(), count);
    same(&solutions, &raw_solutions[..count].to_vec());
    assert_eq!(myb::lambert(GM, &r1, &[1.4e7, 0.0, 0.0], 3600.0, false, 0).unwrap_err().status, Status::Domain);

    // shooting with a closure, against a C function doing the same
    let residual = |u: &[f64], f: &mut [f64], jacobian: Option<&mut [f64]>| -> myb::Result<()> {
        f[0] = u[0] * u[0] - 2.0;
        f[1] = u[0] * u[1] - 3.0;
        if let Some(j) = jacobian {
            j.copy_from_slice(&[2.0 * u[0], 0.0, u[1], u[0]]);
        }
        Ok(())
    };
    unsafe extern "C" fn raw_residual(_c: *mut c_void, u: *const f64, f: *mut f64, j: *mut f64) -> sys::myb_status {
        // SAFETY: u holds two numbers, f two, j four
        unsafe {
            let (u0, u1) = (*u, *u.add(1));
            *f = u0 * u0 - 2.0;
            *f.add(1) = u0 * u1 - 3.0;
            if !j.is_null() {
                for (k, value) in [2.0 * u0, 0.0, u1, u0].into_iter().enumerate() {
                    *j.add(k) = value;
                }
            }
        }
        sys::MYB_OK
    }
    let mut u = [1.0, 1.0];
    let shot = myb::shoot(residual, 2, &mut u, 1e-14, 0)?;
    let mut raw_u = [1.0, 1.0];
    let mut work = vec![0.0; myb::shoot_work_size(2, 2)];
    let (mut iterations, mut left) = (0, 0.0);
    let options = sys::myb_shoot_options { tolerance: 1e-14, max_iterations: 0 };
    // SAFETY: u holds two numbers, work what the library asks for
    assert_eq!(
        unsafe { sys::myb_shoot(Some(raw_residual), ptr::null_mut(), 2, 2, raw_u.as_mut_ptr(), &options, work.as_mut_ptr(), &mut iterations, &mut left) },
        sys::MYB_OK
    );
    assert_eq!((bits(&u), shot.iterations, shot.residual.to_bits()), (bits(&raw_u), iterations, left.to_bits()));
    let failing = |_: &[f64], _: &mut [f64], _: Option<&mut [f64]>| -> myb::Result<()> { Err(Error::new(Status::Range)) };
    assert_eq!(myb::shoot(failing, 2, &mut [1.0, 1.0], 1e-12, 0).unwrap_err().status, Status::Range);

    // first orbits from three positions and from three sights
    let y0 = myb::state_from_elements(GM, &myb::Elements { p: 7.2e6, e: 0.05, i: 0.6, raan: 0.4, argp: 0.3, nu: 0.2 })?;
    let at = |t: f64| myb::kepler_propagate(GM, &y0, t).unwrap();
    let (p1, p2, p3) = (at(0.0), at(600.0), at(1200.0));
    let (r1, r2, r3) = ([p1[0], p1[1], p1[2]], [p2[0], p2[1], p2[2]], [p3[0], p3[1], p3[2]]);
    let (v2, out) = myb::iod_gibbs(GM, &r1, &r2, &r3)?;
    let (mut raw_v2, mut raw_out) = ([0.0; 3], 0.0);
    // SAFETY: arrays of three numbers
    assert_eq!(unsafe { sys::myb_iod_gibbs(GM, r1.as_ptr(), r2.as_ptr(), r3.as_ptr(), raw_v2.as_mut_ptr(), &mut raw_out) }, sys::MYB_OK);
    assert_eq!((bits(&v2), out.to_bits()), (bits(&raw_v2), raw_out.to_bits()));
    let v2 = myb::iod_herrick_gibbs(GM, &r1, 0.0, &r2, 600.0, &r3, 1200.0)?;
    // SAFETY: arrays of three numbers
    assert_eq!(
        unsafe { sys::myb_iod_herrick_gibbs(GM, r1.as_ptr(), 0.0, r2.as_ptr(), 600.0, r3.as_ptr(), 1200.0, raw_v2.as_mut_ptr()) },
        sys::MYB_OK
    );
    assert_eq!(bits(&v2), bits(&raw_v2));
    assert_eq!(myb::iod_gibbs(GM, &[0.0; 3], &r2, &r3).unwrap_err().status, Status::Domain);

    let t = [0.0, 300.0, 600.0];
    let observer = [0.0, 0.003, 0.006].map(|a: f64| [6378137.0 * a.cos(), 6378137.0 * a.sin(), 1000.0]);
    let mut sight = [[0.0; 3]; 3];
    for k in 0..3 {
        let p = at(t[k]);
        let d = [p[0] - observer[k][0], p[1] - observer[k][1], p[2] - observer[k][2]];
        let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        sight[k] = [d[0] / n, d[1] / n, d[2] / n];
    }
    let gauss = myb::iod_gauss(GM, &t, &observer, &sight)?;
    let laplace = myb::iod_laplace(GM, &t, &observer, &sight)?;
    // guesses near the truth: the distances from the centre, the ranges of Gauss's orbit
    let norm = |p: &State| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    let distance = [norm(&at(t[0])) * 1.01, norm(&at(t[1])) * 0.99];
    let range = [gauss[0].range[0] * 1.02, gauss[0].range[2] * 0.98];
    let double_r = myb::iod_double_r(GM, &t, &observer, &sight, &distance)?;
    let gooding = myb::iod_gooding(GM, &t, &observer, &sight, &range, 0, 0)?;
    let mut raw_solutions = [sys::myb_iod_solution::default(); 3];
    let mut raw_one = sys::myb_iod_solution::default();
    // SAFETY: the arrays have the sizes of the C interface
    unsafe {
        assert_eq!(sys::myb_iod_gauss(GM, t.as_ptr(), observer.as_ptr(), sight.as_ptr(), raw_solutions.as_mut_ptr(), &mut count), sys::MYB_OK);
        same(&gauss, &raw_solutions[..count].to_vec());
        assert_eq!(sys::myb_iod_laplace(GM, t.as_ptr(), observer.as_ptr(), sight.as_ptr(), raw_solutions.as_mut_ptr(), &mut count), sys::MYB_OK);
        same(&laplace, &raw_solutions[..count].to_vec());
        assert_eq!(sys::myb_iod_double_r(GM, t.as_ptr(), observer.as_ptr(), sight.as_ptr(), distance.as_ptr(), &mut raw_one), sys::MYB_OK);
        same(&double_r, &raw_one);
        assert_eq!(sys::myb_iod_gooding(GM, t.as_ptr(), observer.as_ptr(), sight.as_ptr(), range.as_ptr(), 0, 0, &mut raw_one), sys::MYB_OK);
        same(&gooding, &raw_one);
    }
    assert!(!gauss.is_empty());
    assert_eq!(myb::iod_gauss(GM, &[0.0, 0.0, 600.0], &observer, &sight).unwrap_err().status, Status::Domain);
    assert_eq!(myb::iod_gooding(GM, &t, &observer, &sight, &range, 9, 1).unwrap_err().status, Status::Domain);
    Ok(())
}

#[test]
fn relative_motion_gives_the_bits_of_the_c_calls() -> myb::Result<()> {
    let chief = orbit(0.9, 0.3);
    let mut deputy = chief;
    for (k, d) in [100.0, -50.0, 20.0, 0.1, 0.0, -0.05].iter().enumerate() {
        deputy[k] += d;
    }
    let relative = myb::relative_from_inertial(&chief, &deputy)?;
    let back = myb::relative_to_inertial(&chief, &relative)?;
    let stm = myb::relative_stm(GM, &chief, 600.0)?;
    let cw = myb::cw_stm(1.1e-3, 600.0);
    let (dv0, dv1) = myb::relative_rendezvous(GM, &chief, 1200.0, &relative, &[0.0, 10.0, 0.0], None)?;
    let roe = myb::roe_from_states(GM, &chief, &deputy)?;
    let position = myb::roe_to_position(7.0e6, 0.5, &roe);
    let (mut r, mut b, mut s, mut c, mut d0, mut d1, mut e, mut p) =
        ([0.0; 6], sys::myb_pv::default(), [0.0; 36], [0.0; 36], [0.0; 3], [0.0; 3], [0.0; 6], [0.0; 3]);
    let target = [0.0, 10.0, 0.0];
    // SAFETY: arrays of the sizes of the C interface
    unsafe {
        assert_eq!(sys::myb_relative_from_inertial(&pv(&chief), &pv(&deputy), r.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_relative_to_inertial(&pv(&chief), r.as_ptr(), &mut b), sys::MYB_OK);
        assert_eq!(sys::myb_relative_stm(GM, &pv(&chief), 600.0, s.as_mut_ptr()), sys::MYB_OK);
        sys::myb_cw_stm(1.1e-3, 600.0, c.as_mut_ptr());
        assert_eq!(
            sys::myb_relative_rendezvous(GM, &pv(&chief), 1200.0, r.as_ptr(), target.as_ptr(), ptr::null(), d0.as_mut_ptr(), d1.as_mut_ptr()),
            sys::MYB_OK
        );
        assert_eq!(sys::myb_roe_from_states(GM, &pv(&chief), &pv(&deputy), e.as_mut_ptr()), sys::MYB_OK);
        sys::myb_roe_to_position(7.0e6, 0.5, e.as_ptr(), p.as_mut_ptr());
    }
    assert_eq!(bits(&relative), bits(&r));
    assert_eq!(bits(&back), bits(&[b.r, b.v].concat()));
    assert_eq!((bits(&stm), bits(&cw)), (bits(&s), bits(&c)));
    assert_eq!((bits(&dv0), bits(&dv1)), (bits(&d0), bits(&d1)));
    assert_eq!((bits(&roe), bits(&position)), (bits(&e), bits(&p)));
    assert_eq!(myb::relative_from_inertial(&[0.0; 6], &deputy).unwrap_err().status, Status::Domain);
    let period = 2.0 * std::f64::consts::PI * (7.0e6_f64.powi(3) / GM).sqrt();
    assert_eq!(myb::relative_rendezvous(GM, &chief, period, &relative, &target, None).unwrap_err().status, Status::Domain);
    Ok(())
}

// ---- Conjunctions ---------------------------------------------------------------------------------

fn line(seconds: f64) -> State {
    [-1000.0 + 10.0 * seconds, 50.0, 0.0, 10.0, 0.0, 0.0]
}

unsafe extern "C" fn raw_still(_c: *mut c_void, _t: sys::myb_time, y: *mut f64) -> sys::myb_status {
    // SAFETY: y holds six numbers
    unsafe { *(y as *mut State) = [0.0; 6] };
    sys::MYB_OK
}

unsafe extern "C" fn raw_line(_c: *mut c_void, t: sys::myb_time, y: *mut f64) -> sys::myb_status {
    // SAFETY: y holds six numbers
    unsafe { *(y as *mut State) = line(sys::myb_time_seconds(t)) };
    sys::MYB_OK
}

fn raw_orbit(model: &ForceModel, y0: &State) -> *mut sys::myb_orbit {
    let options = sys::myb_orbit_options { variational: 1, ..Default::default() };
    let mut raw = ptr::null_mut();
    // SAFETY: the model outlives the orbit in the test
    unsafe {
        assert_eq!(sys::myb_orbit_create(model.as_ptr(), &options, ptr::null(), &mut raw), sys::MYB_OK);
        assert_eq!(sys::myb_orbit_start(raw, 0.0, y0.as_ptr(), 4000.0), sys::MYB_OK);
    }
    raw
}

fn diagonal(position: f64, velocity: f64) -> [f64; 36] {
    let mut p = [0.0; 36];
    for k in 0..6 {
        p[7 * k] = if k < 3 { position } else { velocity };
    }
    p
}

#[test]
fn conjunctions_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let still = |_t: myb::Time| -> myb::Result<State> { Ok([0.0; 6]) };
    let moving = |t: myb::Time| -> myb::Result<State> { Ok(line(myb::seconds(t))) };
    let (approaches, count) = myb::approach_find(&still, &moving, 0.0, 200.0, 7.0, 4)?;
    let mut raw = [sys::myb_approach::default(); 4];
    let mut raw_count = 0usize;
    // SAFETY: the capacity is that of the array
    assert_eq!(
        unsafe { sys::myb_approach_find(Some(raw_still), ptr::null_mut(), Some(raw_line), ptr::null_mut(), 0.0, 200.0, 7.0, raw.as_mut_ptr(), 4, &mut raw_count) },
        sys::MYB_OK
    );
    assert_eq!((count, approaches.len()), (raw_count, 1));
    same(&approaches, &raw[..1].to_vec());
    assert!((approaches[0].t - 100.0).abs() < 1e-6 && (approaches[0].distance - 50.0).abs() < 1e-6);
    assert_eq!(myb::approach_find(&still, &moving, 10.0, 0.0, 1.0, 4).unwrap_err().status, Status::Domain);
    let lost = |_t: myb::Time| -> myb::Result<State> { Err(Error::new(Status::Range)) };
    assert_eq!(myb::approach_find(&still, &lost, 0.0, 10.0, 1.0, 4).unwrap_err().status, Status::Range);

    // the plane of an encounter and its probability
    let (position, velocity, covariance) = ([0.0, 50.0, 10.0], [10.0, 0.0, 0.0], [100.0, 0.0, 0.0, 0.0, 400.0, 0.0, 0.0, 0.0, 900.0]);
    let encounter = myb::encounter_plane(&position, &velocity, &covariance, 20.0)?;
    let mut raw_encounter = sys::myb_encounter::default();
    let (mut p, mut p_max, mut factor) = (0.0, 0.0, 0.0);
    // SAFETY: arrays of the sizes of the C interface
    unsafe {
        assert_eq!(sys::myb_encounter_plane(position.as_ptr(), velocity.as_ptr(), covariance.as_ptr(), 20.0, &mut raw_encounter), sys::MYB_OK);
        assert_eq!(sys::myb_collision_probability(&raw_encounter, &mut p), sys::MYB_OK);
        assert_eq!(sys::myb_collision_probability_max(&raw_encounter, &mut p_max, &mut factor), sys::MYB_OK);
        assert_eq!(myb::collision_probability_chan(&encounter).to_bits(), sys::myb_collision_probability_chan(&raw_encounter).to_bits());
    }
    same(&encounter, &raw_encounter);
    assert_eq!(myb::collision_probability(&encounter)?.to_bits(), p.to_bits());
    let (q, f) = myb::collision_probability_max(&encounter)?;
    assert_eq!((q.to_bits(), f.to_bits()), (p_max.to_bits(), factor.to_bits()));
    assert_eq!(myb::encounter_plane(&position, &velocity, &covariance, 0.0).unwrap_err().status, Status::Domain);

    // the rate of entry and its integral over a constant Gaussian
    let mean = [0.0, 30.0, 0.0, 5.0, 0.0, 0.0];
    let cov = diagonal(400.0, 1.0);
    let rate = myb::collision_rate(&mean, &cov, 20.0)?;
    let mut raw_rate = 0.0;
    // SAFETY: arrays of the sizes of the C interface
    assert_eq!(unsafe { sys::myb_collision_rate(mean.as_ptr(), cov.as_ptr(), 20.0, &mut raw_rate) }, sys::MYB_OK);
    assert_eq!(rate.to_bits(), raw_rate.to_bits());
    unsafe extern "C" fn raw_gaussian(_c: *mut c_void, t: f64, mean: *mut f64, covariance: *mut f64) -> sys::myb_status {
        // SAFETY: the outputs hold 6 and 36 numbers
        unsafe {
            *(mean as *mut [f64; 6]) = [-5.0 * (t - 5.0), 30.0, 0.0, 5.0, 0.0, 0.0];
            *(covariance as *mut [f64; 36]) = diagonal(400.0, 1.0);
        }
        sys::MYB_OK
    }
    let gaussian = |t: f64| -> myb::Result<([f64; 6], [f64; 36])> { Ok(([-5.0 * (t - 5.0), 30.0, 0.0, 5.0, 0.0, 0.0], diagonal(400.0, 1.0))) };
    let (long, error) = myb::collision_probability_long(gaussian, 0.0, 10.0, 1.0, 20.0)?;
    let (mut raw_long, mut raw_error) = (0.0, 0.0);
    // SAFETY: valid places for the outputs
    assert_eq!(
        unsafe { sys::myb_collision_probability_long(Some(raw_gaussian), ptr::null_mut(), 0.0, 10.0, 1.0, 20.0, &mut raw_long, &mut raw_error) },
        sys::MYB_OK
    );
    assert_eq!((long.to_bits(), error.to_bits()), (raw_long.to_bits(), raw_error.to_bits()));
    assert_eq!(myb::collision_rate(&mean, &cov, -1.0).unwrap_err().status, Status::Domain);

    // two orbits that cross near their nodes, with their uncertainty
    let model = Rc::new(ForceModel::new()?);
    model.add_point_mass("earth", GM)?;
    let options = OrbitOptions { variational: true, ..Default::default() };
    let (a0, b0) = (orbit(0.0, 0.0), orbit(0.1, -1.0e-4));
    let (pa, pb) = (diagonal(100.0, 1.0e-4), diagonal(400.0, 4.0e-4));
    let mut a = Orbit::new(&model, &options)?;
    let mut b = Orbit::new(&model, &options)?;
    a.start(0.0, &a0, 4000.0)?;
    b.start(0.0, &b0, 4000.0)?;
    let (conjunctions, count) = myb::conjunction_find(&mut a, &pa, &mut b, &pb, 10.0, 3200.0, 10.0, 8)?;
    let (ra, rb) = (raw_orbit(&model, &a0), raw_orbit(&model, &b0));
    let mut raw_conjunctions = [sys::myb_conjunction::default(); 8];
    // SAFETY: two distinct orbits; the capacity is that of the array
    assert_eq!(
        unsafe { sys::myb_conjunction_find(ra, pa.as_ptr(), rb, pb.as_ptr(), 10.0, 3200.0, 10.0, raw_conjunctions.as_mut_ptr(), 8, &mut raw_count) },
        sys::MYB_OK
    );
    assert!(count >= 2, "{count}");
    assert_eq!(count, raw_count);
    same(&conjunctions, &raw_conjunctions[..count.min(8)].to_vec());
    let (probability, error) = myb::conjunction_long(&mut a, &pa, &mut b, &pb, 10.0, 3300.0, 20.0)?;
    let (mut raw_probability, mut raw_error) = (0.0, 0.0);
    // SAFETY: two distinct orbits standing at one time
    assert_eq!(
        unsafe { sys::myb_conjunction_long(ra, pa.as_ptr(), rb, pb.as_ptr(), 10.0, 3300.0, 20.0, &mut raw_probability, &mut raw_error) },
        sys::MYB_OK
    );
    assert_eq!((probability.to_bits(), error.to_bits()), (raw_probability.to_bits(), raw_error.to_bits()));
    // and the approaches of the two, without the uncertainty
    let (approaches, count) = myb::approach_orbits(&mut a, &mut b, 3900.0, 10.0, 4)?;
    // SAFETY: as above
    assert_eq!(unsafe { sys::myb_approach_orbits(ra, rb, 3900.0, 10.0, raw.as_mut_ptr(), 4, &mut raw_count) }, sys::MYB_OK);
    assert_eq!(count, raw_count);
    same(&approaches, &raw[..count.min(4)].to_vec());
    // SAFETY: created above, destroyed once
    unsafe {
        sys::myb_orbit_destroy(ra);
        sys::myb_orbit_destroy(rb);
    }
    assert_eq!(myb::approach_orbits(&mut a, &mut b, 3800.0, 10.0, 8).unwrap_err().status, Status::Domain);
    Ok(())
}

#[test]
fn uncertainty_gives_the_bits_of_the_c_calls() -> myb::Result<()> {
    let model = Rc::new(ForceModel::new()?);
    model.add_point_mass("earth", GM)?;
    let y0 = orbit(0.5, 0.2);
    let mut variational = Variational::new(&model, &VariationalOptions::default())?;
    variational.start(0.0, &y0)?;
    variational.integrate(600.0)?;
    let end = variational.get();
    let dense = variational.dense(variational.time() - 1.0)?;
    let (mut y, mut phi, mut yd, mut phid) = ([0.0; 6], [0.0; 36], [0.0; 6], [0.0; 36]);
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    unsafe {
        let options = sys::myb_variational_options::default();
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_variational_create(model.as_ptr(), &options, ptr::null(), &mut raw), sys::MYB_OK);
        assert_eq!(sys::myb_variational_start(raw, 0.0, y0.as_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_variational_integrate(raw, 600.0), sys::MYB_OK);
        sys::myb_variational_get(raw, y.as_mut_ptr(), phi.as_mut_ptr(), ptr::null_mut());
        assert_eq!(sys::myb_variational_dense(raw, sys::myb_variational_time(raw) - 1.0, yd.as_mut_ptr(), phid.as_mut_ptr(), ptr::null_mut()), sys::MYB_OK);
        sys::myb_variational_destroy(raw);
    }
    assert_eq!((bits(&end.y), bits(&end.phi)), (bits(&y), bits(&phi)));
    assert_eq!((bits(&dense.y), bits(&dense.phi)), (bits(&yd), bits(&phid)));
    assert!(end.sensitivity.is_empty());
    assert_eq!(variational.dense(-100.0).unwrap_err().status, Status::Range);
    variational.step(700.0)?;

    // covariances
    let p = {
        let mut p = [0.0; 36];
        for i in 0..6 {
            for j in 0..6 {
                p[6 * i + j] = if i == j { 10.0 + i as f64 } else { 0.5 / (1.0 + (i + j) as f64) };
            }
        }
        p
    };
    let m: Vec<f64> = (0..36).map(|k| (k as f64 * 0.37).sin()).collect();
    let transformed = myb::covariance_transform(6, &m, &p)?;
    let rtn = myb::covariance_rtn(&[y0[0], y0[1], y0[2]], &[y0[3], y0[4], y0[5]], &p)?;
    let (mut out, mut raw_rtn) = (vec![0.0; 36], [0.0; 36]);
    // SAFETY: matrices of 36 numbers
    unsafe {
        sys::myb_covariance_transform(6, m.as_ptr(), p.as_ptr(), out.as_mut_ptr());
        assert_eq!(sys::myb_covariance_rtn(y0.as_ptr(), y0[3..].as_ptr(), p.as_ptr(), raw_rtn.as_mut_ptr()), sys::MYB_OK);
    }
    assert_eq!((bits(&transformed), bits(&rtn)), (bits(&out), bits(&raw_rtn)));
    assert_eq!(myb::covariance_transform(6, &m[..35], &p).unwrap_err().status, Status::Domain);
    assert_eq!(myb::covariance_rtn(&[1.0, 0.0, 0.0], &[2.0, 0.0, 0.0], &p).unwrap_err().status, Status::Domain);

    // the generator and statistics
    let mut rng = Rng::new(42);
    let mut raw = sys::myb_rng::default();
    // SAFETY: a valid generator
    unsafe { sys::myb_rng_seed(&mut raw, 42) };
    let l = [2.0, 0.0, 0.5, 1.0];
    let mut samples = Vec::new();
    for _ in 0..50 {
        let mut g = [0.0; 2];
        // SAFETY: the generator is valid; arrays of the sizes asked for
        unsafe {
            assert_eq!(rng.next_u64(), sys::myb_rng_next(&mut raw));
            assert_eq!(rng.uniform().to_bits(), sys::myb_rng_uniform(&mut raw).to_bits());
            assert_eq!(rng.normal().to_bits(), sys::myb_rng_normal(&mut raw).to_bits());
            sys::myb_rng_gaussian(&mut raw, 2, [1.0, -1.0].as_ptr(), l.as_ptr(), g.as_mut_ptr());
        }
        let sample = rng.gaussian(&[1.0, -1.0], &l)?;
        assert_eq!(bits(&sample), bits(&g));
        samples.extend(sample);
    }
    let (mean, covariance) = myb::sample_statistics(2, &samples)?;
    let (mut raw_mean, mut raw_covariance) = ([0.0; 2], [0.0; 4]);
    // SAFETY: 50 samples of two numbers
    assert_eq!(unsafe { sys::myb_sample_statistics(50, 2, samples.as_ptr(), raw_mean.as_mut_ptr(), raw_covariance.as_mut_ptr()) }, sys::MYB_OK);
    assert_eq!((bits(&mean), bits(&covariance)), (bits(&raw_mean), bits(&raw_covariance)));
    assert_eq!(myb::sample_statistics(2, &samples[..2]).unwrap_err().status, Status::Domain);
    assert_eq!(rng.gaussian(&[0.0; 2], &l[..3]).unwrap_err().status, Status::Domain);
    Ok(())
}
