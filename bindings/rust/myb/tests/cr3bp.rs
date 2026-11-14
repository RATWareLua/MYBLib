// SPDX-License-Identifier: BSD-3-Clause
//! The circular restricted three-body problem (`myb/cr3bp.h`) through the safe types against
//! the same C calls, bit for bit; failures as `Err` with the status of the library.
#![allow(clippy::excessive_precision)] // the states of the JPL catalogue, as printed there
use std::fmt::Debug;
use std::rc::Rc;

use myb::sys;
use myb::{Cr3bp, Cr3bpCorrection, Cr3bpSection, Cr3bpSymmetry, Cr3bpTorusFamily, ForceModel, State, Status};

const MU: f64 = 1.215058560962404e-02;

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn same<T: Debug>(a: &T, b: &T) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

fn correction() -> Cr3bpCorrection {
    Cr3bpCorrection { symmetry: Cr3bpSymmetry::Xz, fixed: 0, crossing: 1, tolerance: 0.0, max_iterations: 0 }
}

fn raw_correction() -> sys::myb_cr3bp_correction {
    sys::myb_cr3bp_correction { symmetry: sys::MYB_CR3BP_SYMMETRY_XZ, fixed: 0, crossing: 1, tolerance: 0.0, max_iterations: 0 }
}

/// A periodic orbit corrected both ways: the state, the period.
fn corrected(problem: &Cr3bp, guess: State) -> (State, f64) {
    let (y0, period) = problem.correct(&correction(), &guess).unwrap();
    let (mut y, mut t) = (guess, 0.0);
    // SAFETY: y holds six numbers
    assert_eq!(unsafe { sys::myb_cr3bp_correct(MU, &raw_correction(), y.as_mut_ptr(), &mut t) }, sys::MYB_OK);
    assert_eq!((bits(&y0), period.to_bits()), (bits(&y), t.to_bits()));
    (y0, period)
}

#[test]
fn periodic_orbits_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let problem = Cr3bp::new(MU);
    let mut points = [[0.0; 3]; 5];
    // SAFETY: five points of three numbers
    assert_eq!(unsafe { sys::myb_cr3bp_lagrange(MU, points.as_mut_ptr()) }, sys::MYB_OK);
    assert_eq!(problem.lagrange()?, points);

    let (y0, period) = corrected(&problem, [8.2522726142834313e-01, 0.0, 0.0, 0.0, 1.0748883813645686e-01, 0.0]);
    // SAFETY: arrays of the sizes of the C interface
    unsafe {
        assert_eq!(problem.jacobi(&y0).to_bits(), sys::myb_cr3bp_jacobi(MU, y0.as_ptr()).to_bits());
        let (x, v) = ([y0[0], y0[1], y0[2]], [y0[3], y0[4], y0[5]]);
        let (mut a, mut dadx, mut dadv) = ([0.0; 3], [0.0; 9], [0.0; 9]);
        sys::myb_cr3bp_acceleration(MU, x.as_ptr(), v.as_ptr(), a.as_mut_ptr());
        sys::myb_cr3bp_jacobian(MU, x.as_ptr(), dadx.as_mut_ptr(), dadv.as_mut_ptr());
        assert_eq!(bits(&problem.acceleration(&x, &v)), bits(&a));
        let (jx, jv) = problem.jacobian(&x);
        assert_eq!((bits(&jx), bits(&jv)), (bits(&dadx), bits(&dadv)));
        // the term in a force model, against the same model made in C
        let model = ForceModel::new()?;
        model.add_cr3bp("cr3bp", MU)?;
        let mut mu = MU;
        let mut raw = std::ptr::null_mut();
        let mut in_model = [0.0; 3];
        assert_eq!(sys::myb_force_model_create(std::ptr::null(), &mut raw), sys::MYB_OK);
        let context = (&mut mu as *mut f64).cast();
        assert_eq!(
            sys::myb_force_model_add(raw, c"cr3bp".as_ptr(), Some(sys::myb_cr3bp_term), Some(sys::myb_cr3bp_term_jacobian), context, std::ptr::null_mut()),
            sys::MYB_OK
        );
        assert_eq!(sys::myb_force_model_acceleration(raw, 0.0, x.as_ptr(), v.as_ptr(), in_model.as_mut_ptr()), sys::MYB_OK);
        sys::myb_force_model_destroy(raw);
        assert_eq!(bits(&model.acceleration(0.0, &x, &v)?), bits(&in_model));
    }

    let monodromy = problem.monodromy(&y0, period)?;
    let stability = myb::cr3bp_stability(&monodromy)?;
    let manifold = problem.manifold(&y0, period, 0.25, true, 1, 1e-6)?;
    let (mut m, mut re, mut im, mut index, mut w) = ([0.0; 36], [0.0; 6], [0.0; 6], 0.0, [0.0; 6]);
    // SAFETY: arrays of the sizes of the C interface
    unsafe {
        assert_eq!(sys::myb_cr3bp_monodromy(MU, y0.as_ptr(), period, m.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_cr3bp_stability(m.as_ptr(), re.as_mut_ptr(), im.as_mut_ptr(), &mut index), sys::MYB_OK);
        assert_eq!(sys::myb_cr3bp_manifold(MU, y0.as_ptr(), period, 0.25, 1, 1, 1e-6, w.as_mut_ptr()), sys::MYB_OK);
    }
    assert_eq!(bits(&monodromy), bits(&m));
    assert_eq!((bits(&stability.re), bits(&stability.im), stability.index.to_bits()), (bits(&re), bits(&im), index.to_bits()));
    assert_eq!(bits(&manifold), bits(&w));

    // Richardson's halo
    let (halo, halo_period) = problem.richardson(1, 0.02, 1)?;
    let (mut h, mut ht) = ([0.0; 6], 0.0);
    // SAFETY: h holds six numbers
    assert_eq!(unsafe { sys::myb_cr3bp_richardson(MU, 1, 0.02, 1, h.as_mut_ptr(), &mut ht) }, sys::MYB_OK);
    assert_eq!((bits(&halo), halo_period.to_bits()), (bits(&h), ht.to_bits()));

    // half a period to the section y = 0
    let section = Cr3bpSection { axis: 1, value: 0.0, direction: 0 };
    let reached = problem.propagate(&y0, period, Some(&section))?;
    let raw_section = sys::myb_cr3bp_section { axis: 1, value: 0.0, direction: 0 };
    let (mut y, mut elapsed, mut crossed) = ([0.0; 6], 0.0, 0);
    // SAFETY: valid places for the outputs
    assert_eq!(
        unsafe { sys::myb_cr3bp_propagate(MU, y0.as_ptr(), period, &raw_section, y.as_mut_ptr(), &mut elapsed, &mut crossed) },
        sys::MYB_OK
    );
    assert!(reached.crossed && crossed == 1);
    assert_eq!((bits(&reached.y), reached.elapsed.to_bits()), (bits(&y), elapsed.to_bits()));
    let end = problem.propagate(&y0, 0.1, None)?;
    assert!(!end.crossed && end.elapsed == 0.1);

    // multiple shooting and continuation
    let mut orbit = problem.orbit_from_state(&y0, period, 4)?;
    problem.orbit_correct(&mut orbit, 0.0, 0)?;
    problem.orbit_continue(&mut orbit, 1e-3, 0.0, 0)?;
    let mut raw_orbit = sys::myb_cr3bp_orbit::default();
    // SAFETY: valid places for the orbit
    unsafe {
        assert_eq!(sys::myb_cr3bp_orbit_from_state(MU, y0.as_ptr(), period, 4, &mut raw_orbit), sys::MYB_OK);
        assert_eq!(sys::myb_cr3bp_orbit_correct(MU, &mut raw_orbit, 0.0, 0), sys::MYB_OK);
        assert_eq!(sys::myb_cr3bp_orbit_continue(MU, &mut raw_orbit, 1e-3, 0.0, 0), sys::MYB_OK);
    }
    same(&orbit, &raw_orbit);

    // the pulsating frame of two real bodies
    let earth = [-4.67e6, 0.0, 0.0, 0.0, -12.5, 0.0];
    let moon = [3.797e8, 0.0, 0.0, 0.0, 1.0e3, 10.0];
    let state = [3.2e8, 1.0e7, -2.0e6, 10.0, 800.0, -5.0];
    let rotating = myb::cr3bp_from_ephemeris(&earth, &moon, 3.986e14, 4.9e12, &state)?;
    let back = myb::cr3bp_to_ephemeris(&earth, &moon, 3.986e14, 4.9e12, &rotating)?;
    let (mut r, mut s) = ([0.0; 6], [0.0; 6]);
    // SAFETY: arrays of six numbers
    unsafe {
        assert_eq!(sys::myb_cr3bp_from_ephemeris(earth.as_ptr(), moon.as_ptr(), 3.986e14, 4.9e12, state.as_ptr(), r.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_cr3bp_to_ephemeris(earth.as_ptr(), moon.as_ptr(), 3.986e14, 4.9e12, r.as_ptr(), s.as_mut_ptr()), sys::MYB_OK);
    }
    assert_eq!((bits(&rotating), bits(&back)), (bits(&r), bits(&s)));

    // failures
    assert_eq!(Cr3bp::new(0.7).lagrange().unwrap_err().status, Status::Domain);
    assert_eq!(problem.richardson(3, 0.01, 1).unwrap_err().status, Status::Domain);
    let bad = Cr3bpSection { axis: 7, value: 0.0, direction: 0 };
    assert_eq!(problem.propagate(&y0, 1.0, Some(&bad)).unwrap_err().status, Status::Domain);
    assert_eq!(problem.orbit_from_state(&y0, period, 11).unwrap_err().status, Status::Domain);
    let asymmetric = [y0[0], 0.1, 0.0, 0.0, y0[4], 0.0];
    assert_eq!(problem.correct(&correction(), &asymmetric).unwrap_err().status, Status::Domain);
    assert_eq!(myb::cr3bp_from_ephemeris(&earth, &earth, 3.986e14, 4.9e12, &state).unwrap_err().status, Status::Domain);
    Ok(())
}

#[test]
fn tori_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let problem = Rc::new(Cr3bp::new(MU));
    // a Lissajous torus about the Lyapunov orbit: hyperbolic
    let (lyapunov, period) = corrected(&problem, [8.2522726142834313e-01, 0.0, 0.0, 0.0, 1.0748883813645686e-01, 0.0]);
    let mut torus = problem.torus_from_orbit(&lyapunov, period, 15, 0, 1e-3)?;
    problem.torus_continue(&mut torus, Cr3bpTorusFamily::Jacobi, 0.0, 0.0, 0)?;
    let spectrum = problem.torus_stability(&torus, 0.0)?;
    let point = myb::cr3bp_torus_state(&torus, 0.3);
    let manifold = myb::cr3bp_torus_manifold(&torus, &spectrum, 0.3, true, 1, 1e-6)?;
    let mut raw = Box::<sys::myb_cr3bp_torus>::default();
    let mut raw_spectrum = Box::<sys::myb_cr3bp_torus_spectrum>::default();
    let (mut y, mut w) = ([0.0; 6], [0.0; 6]);
    // SAFETY: valid places for the torus and its spectrum
    unsafe {
        assert_eq!(sys::myb_cr3bp_torus_from_orbit(MU, lyapunov.as_ptr(), period, 15, 0, 1e-3, &mut *raw), sys::MYB_OK);
        assert_eq!(sys::myb_cr3bp_torus_continue(MU, &mut *raw, sys::MYB_CR3BP_TORUS_JACOBI, 0.0, 0.0, 0), sys::MYB_OK);
        assert_eq!(sys::myb_cr3bp_torus_stability(MU, &*raw, 0.0, &mut *raw_spectrum), sys::MYB_OK);
        sys::myb_cr3bp_torus_state(&*raw, 0.3, y.as_mut_ptr());
        assert_eq!(sys::myb_cr3bp_torus_manifold(&*raw, &*raw_spectrum, 0.3, 1, 1, 1e-6, w.as_mut_ptr()), sys::MYB_OK);
    }
    same(&*torus, &*raw);
    same(&*spectrum, &*raw_spectrum);
    assert!(spectrum.unstable > 1.0 && spectrum.hyperbolic > 0);
    assert_eq!((bits(&point), bits(&manifold)), (bits(&y), bits(&w)));

    // a torus about a distant retrograde orbit: elliptic, with its normal rotation
    let (dro, dro_period) = corrected(&problem, [8.2915209119468902e-01, 0.0, 0.0, 0.0, 4.9478973105640744e-01, 0.0]);
    let mut elliptic = problem.torus_from_orbit(&dro, dro_period, 21, 0, 1e-3)?;
    problem.torus_continue(&mut elliptic, Cr3bpTorusFamily::Jacobi, 0.0, 0.0, 0)?;
    let spectrum = problem.torus_stability(&elliptic, 0.0)?;
    // SAFETY: valid places for the torus and its spectrum
    unsafe {
        assert_eq!(sys::myb_cr3bp_torus_from_orbit(MU, dro.as_ptr(), dro_period, 21, 0, 1e-3, &mut *raw), sys::MYB_OK);
        assert_eq!(sys::myb_cr3bp_torus_continue(MU, &mut *raw, sys::MYB_CR3BP_TORUS_JACOBI, 0.0, 0.0, 0), sys::MYB_OK);
        assert_eq!(sys::myb_cr3bp_torus_stability(MU, &*raw, 0.0, &mut *raw_spectrum), sys::MYB_OK);
    }
    same(&*elliptic, &*raw);
    same(&*spectrum, &*raw_spectrum);
    assert_eq!(
        (spectrum.normal_rotation.to_bits(), bits(&spectrum.normal_eigenvalue), spectrum.normal_frequency.to_bits()),
        (raw_spectrum.normal_rotation.to_bits(), bits(&raw_spectrum.normal_eigenvalue), raw_spectrum.normal_frequency.to_bits())
    );
    assert_eq!((bits(&spectrum.normal_re), bits(&spectrum.normal_im)), (bits(&raw_spectrum.normal_re), bits(&raw_spectrum.normal_im)));
    assert!(spectrum.elliptic > 0 && spectrum.normal_rotation != 0.0);
    // no hyperbolic pair: no manifold
    assert_eq!(myb::cr3bp_torus_manifold(&elliptic, &spectrum, 0.0, true, 1, 1e-6).unwrap_err().status, Status::Domain);
    assert_eq!(problem.torus_from_orbit(&dro, dro_period, 20, 0, 1e-3).unwrap_err().status, Status::Domain);
    assert_eq!(problem.torus_from_orbit(&dro, dro_period, 21, 0, 0.0).unwrap_err().status, Status::Domain);
    assert_eq!(problem.torus_stability(&elliptic, -1.0).unwrap_err().status, Status::Domain);
    Ok(())
}
