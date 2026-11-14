// SPDX-License-Identifier: BSD-3-Clause
//! Low thrust (`myb/optimal.h`, `myb/thrust.h`) through the safe layer against the same C
//! calls, bit for bit; failures as `Err` with the status of the library.
use std::ffi::c_void;
use std::fmt::Debug;
use std::ptr;
use std::rc::Rc;

use myb::sys;
use myb::{
    Chain, ChainBody, ChainLeg, ForceModel, IndirectProblem, NlpProblem, Qlaw, QlawFlight, QlawTerm, SimsFlanagan, State, Status,
    Steering, Thrust, ThrustArc,
};

const GM: f64 = 3.986004418e14;
const SUN: f64 = 1.32712440018e20;
const AU: f64 = 1.495978707e11;

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn same<T: Debug>(a: &T, b: &T) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

fn on_circle(mu: f64, r: f64, angle: f64, tilt: f64) -> State {
    let v = (mu / r).sqrt();
    [
        r * angle.cos(),
        r * angle.sin() * tilt.cos(),
        r * angle.sin() * tilt.sin(),
        -v * angle.sin(),
        v * angle.cos() * tilt.cos(),
        v * angle.cos() * tilt.sin(),
    ]
}

fn raw_model() -> *mut sys::myb_force_model {
    let mut raw = ptr::null_mut();
    // SAFETY: the calls of the C interface with valid arguments
    unsafe {
        assert_eq!(sys::myb_force_model_create(ptr::null(), &mut raw), sys::MYB_OK);
    }
    raw
}

#[test]
fn the_indirect_method_and_the_solver_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let problem = IndirectProblem {
        mu: GM,
        thrust: 300.0,
        exhaust: 3000.0,
        t0: 100.0,
        start: on_circle(GM, 7.0e6, 0.0, 0.0),
        mass0: 1000.0,
        t1: 4600.0,
        target: on_circle(GM, 8.0e6, 4.4, 0.06),
        smoothing: 1.0,
        ..Default::default()
    };
    let solution = problem.solve(None)?;
    let times = [100.0, 1000.0, 2500.0, 4000.0, solution.t1];
    let points = problem.sample(&solution, &times)?;
    let raw = sys::myb_indirect_problem {
        mu: GM,
        field: ptr::null_mut(),
        thrust: 300.0,
        exhaust: 3000.0,
        t0: 100.0,
        start: problem.start,
        mass0: 1000.0,
        t1: 4600.0,
        target: problem.target,
        goal: sys::MYB_INDIRECT_FUEL,
        smoothing: 1.0,
        tolerance: 0.0,
        tries: 0,
        arrival: sys::MYB_INDIRECT_RENDEZVOUS,
        phase: 0.0,
    };
    let mut raw_solution = sys::myb_indirect_solution::default();
    let mut raw_points = [sys::myb_indirect_point::default(); 5];
    // SAFETY: the calls of the C interface with valid arguments
    unsafe {
        assert_eq!(sys::myb_indirect_solve(&raw, ptr::null(), ptr::null(), &mut raw_solution), sys::MYB_OK);
        assert_eq!(sys::myb_indirect_sample(&raw, &raw_solution, times.as_ptr(), 5, ptr::null(), raw_points.as_mut_ptr()), sys::MYB_OK);
    }
    same(&solution, &raw_solution);
    same(&points, &raw_points.to_vec());
    // from the costates of the solution, at once
    let again = problem.solve(Some(&solution.costate))?;
    assert!((again.mass1 - solution.mass1).abs() < 1e-6);
    let broken = IndirectProblem { thrust: 0.0, ..problem.clone() };
    assert_eq!(broken.solve(None).unwrap_err().status, Status::Domain);
    assert_eq!(problem.sample(&solution, &[4000.0, 1000.0]).unwrap_err().status, Status::Domain);

    // a small nonlinear program: a closure against a C function doing the same
    let program = NlpProblem { equalities: 1, inequalities: 1, lower: Some(&[-10.0, -10.0]), upper: None };
    let function = |x: &[f64], c: &mut [f64], derivatives: Option<(&mut [f64], &mut [f64])>| -> myb::Result<f64> {
        c[0] = x[0] + x[1] - 2.0;
        c[1] = x[0] - 0.8;
        if let Some((gradient, jacobian)) = derivatives {
            gradient.copy_from_slice(&[2.0 * (x[0] - 1.0), 2.0 * (x[1] - 2.0)]);
            jacobian.copy_from_slice(&[1.0, 1.0, 1.0, 0.0]);
        }
        Ok((x[0] - 1.0) * (x[0] - 1.0) + (x[1] - 2.0) * (x[1] - 2.0))
    };
    unsafe extern "C" fn raw_function(_c: *mut c_void, x: *const f64, f: *mut f64, c: *mut f64, g: *mut f64, j: *mut f64) -> sys::myb_status {
        // SAFETY: x holds two numbers, c two, g two and j four
        unsafe {
            let (x0, x1) = (*x, *x.add(1));
            *c = x0 + x1 - 2.0;
            *c.add(1) = x0 - 0.8;
            if !g.is_null() {
                *g = 2.0 * (x0 - 1.0);
                *g.add(1) = 2.0 * (x1 - 2.0);
                for (k, value) in [1.0, 1.0, 1.0, 0.0].into_iter().enumerate() {
                    *j.add(k) = value;
                }
            }
            *f = (x0 - 1.0) * (x0 - 1.0) + (x1 - 2.0) * (x1 - 2.0);
        }
        sys::MYB_OK
    }
    let mut x = [0.0, 0.0];
    let (result, multipliers) = myb::nlp_solve(&program, function, None, &mut x)?;
    let lower = [-10.0, -10.0];
    let raw_program = sys::myb_nlp_problem {
        n: 2,
        equalities: 1,
        inequalities: 1,
        lower: lower.as_ptr(),
        upper: ptr::null(),
        function: Some(raw_function),
        context: ptr::null_mut(),
    };
    let (mut raw_x, mut raw_multipliers, mut raw_result) = ([0.0, 0.0], [0.0; 2], sys::myb_nlp_result::default());
    // SAFETY: x holds two numbers, the multipliers two
    assert_eq!(
        unsafe { sys::myb_nlp_solve(&raw_program, ptr::null(), raw_x.as_mut_ptr(), raw_multipliers.as_mut_ptr(), ptr::null(), &mut raw_result) },
        sys::MYB_OK
    );
    assert_eq!((bits(&x), bits(&multipliers)), (bits(&raw_x), bits(&raw_multipliers)));
    same(&result, &raw_result);
    assert!((x[0] - 0.5).abs() < 1e-5 && (x[1] - 1.5).abs() < 1e-5, "{x:?}");
    let short = NlpProblem { lower: Some(&[0.0]), ..program };
    assert_eq!(myb::nlp_solve(&short, function, None, &mut [0.0, 0.0]).unwrap_err().status, Status::Domain);
    let failing = |_: &[f64], _: &mut [f64], _: Option<(&mut [f64], &mut [f64])>| -> myb::Result<f64> { Err(myb::Error::new(Status::Range)) };
    assert_eq!(myb::nlp_solve(&program, failing, None, &mut [0.0, 0.0]).unwrap_err().status, Status::Range);
    Ok(())
}

fn raw_leg(leg: &SimsFlanagan, field: *mut sys::myb_force_model) -> sys::myb_sims_flanagan {
    sys::myb_sims_flanagan {
        mu: leg.mu,
        thrust: leg.thrust,
        exhaust: leg.exhaust,
        t0: leg.t0,
        start: leg.start,
        mass0: leg.mass0,
        t1: leg.t1,
        target: leg.target,
        segments: leg.segments,
        field,
    }
}

#[test]
fn sims_flanagan_gives_the_bits_of_the_c_calls() -> myb::Result<()> {
    let model = Rc::new(ForceModel::new()?);
    model.add_point_mass("earth", GM)?;
    let mut leg = SimsFlanagan {
        mu: GM,
        thrust: 0.5,
        exhaust: 30000.0,
        t0: 0.0,
        start: on_circle(GM, 7.0e6, 0.0, 0.0),
        mass0: 500.0,
        t1: 3000.0,
        target: on_circle(GM, 7.05e6, 2.6, 0.0),
        segments: 6,
        field: None,
    };
    let throttles: Vec<f64> = (0..18).map(|k| 0.1 * ((k as f64) * 0.7).sin()).collect();
    for with_field in [false, true] {
        leg.field = with_field.then(|| Rc::clone(&model));
        let mismatch = leg.mismatch(&throttles)?;
        let raw_field = raw_model();
        // SAFETY: the calls of the C interface with valid arguments
        let (m, mass1, jacobian) = unsafe {
            assert_eq!(sys::myb_force_model_add(raw_field, c"earth".as_ptr(), Some(sys::myb_term_point_mass), Some(sys::myb_term_point_mass_jacobian), (&GM as *const f64 as *mut f64).cast(), ptr::null_mut()), sys::MYB_OK);
            let made = raw_leg(&leg, if with_field { raw_field } else { ptr::null_mut() });
            let (mut m, mut mass1, mut jacobian) = ([0.0; 6], 0.0, vec![0.0; 108]);
            assert_eq!(sys::myb_sims_flanagan_mismatch(&made, throttles.as_ptr(), m.as_mut_ptr(), &mut mass1, jacobian.as_mut_ptr()), sys::MYB_OK);
            sys::myb_force_model_destroy(raw_field);
            (m, mass1, jacobian)
        };
        assert_eq!(bits(&mismatch.mismatch), bits(&m));
        assert_eq!((mismatch.mass1.to_bits(), bits(&mismatch.jacobian)), (mass1.to_bits(), bits(&jacobian)));
    }
    leg.field = None;
    let mut solved = vec![0.0; 18];
    let outcome = leg.solve(None, &mut solved);
    let mut raw_solved = vec![0.0; 18];
    let mut raw_result = sys::myb_sims_flanagan_result::default();
    let made = raw_leg(&leg, ptr::null_mut());
    // SAFETY: the throttles hold three numbers per segment
    let code = unsafe { sys::myb_sims_flanagan_solve(&made, ptr::null(), raw_solved.as_mut_ptr(), ptr::null(), &mut raw_result) };
    match outcome {
        Ok(result) => {
            assert_eq!(code, sys::MYB_OK);
            same(&result, &raw_result);
        }
        Err(error) => assert_eq!(error.status.raw(), code),
    }
    assert_eq!(bits(&solved), bits(&raw_solved));
    assert_eq!(leg.mismatch(&throttles[..17]).unwrap_err().status, Status::Domain);
    let heavy = vec![2.0; 18];
    assert_eq!(leg.mismatch(&heavy).unwrap_err().status, Status::Domain);
    Ok(())
}

// ---- Chains ---------------------------------------------------------------------------------------

fn planet(radius: f64, phase: f64, t: myb::Time) -> State {
    on_circle(SUN, radius, phase + (SUN / radius.powi(3)).sqrt() * myb::seconds(t), 0.0)
}

unsafe extern "C" fn raw_earth(_c: *mut c_void, t: sys::myb_time, y: *mut f64) -> sys::myb_status {
    // SAFETY: y holds six numbers
    unsafe { *(y as *mut State) = planet(AU, 0.0, t) };
    sys::MYB_OK
}

unsafe extern "C" fn raw_venus(_c: *mut c_void, t: sys::myb_time, y: *mut f64) -> sys::myb_status {
    // SAFETY: y holds six numbers
    unsafe { *(y as *mut State) = planet(0.723 * AU, 1.0, t) };
    sys::MYB_OK
}

unsafe extern "C" fn raw_mars(_c: *mut c_void, t: sys::myb_time, y: *mut f64) -> sys::myb_status {
    // SAFETY: y holds six numbers
    unsafe { *(y as *mut State) = planet(1.524 * AU, 0.9, t) };
    sys::MYB_OK
}

#[test]
fn chains_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let earth = ChainBody { state: Rc::new(|t: myb::Time| -> myb::Result<State> { Ok(planet(AU, 0.0, t)) }), gm: 0.0, periapsis: 0.0, impulse_exhaust: 0.0 };
    let venus = ChainBody {
        state: Rc::new(|t: myb::Time| -> myb::Result<State> { Ok(planet(0.723 * AU, 1.0, t)) }),
        gm: 3.24859e14,
        periapsis: 6.5e6,
        impulse_exhaust: 3000.0,
    };
    let mars = ChainBody { state: Rc::new(|t: myb::Time| -> myb::Result<State> { Ok(planet(1.524 * AU, 0.9, t)) }), gm: 0.0, periapsis: 0.0, impulse_exhaust: 0.0 };
    let day = 86400.0;
    let mut chain = Chain {
        mu: SUN,
        field: None,
        epoch: myb::time_from_seconds(0.0),
        thrust: 1.0,
        exhaust: 30000.0,
        mass0: 1000.0,
        legs: vec![ChainLeg { segments: 4, flight: [150.0 * day, 170.0 * day] }, ChainLeg { segments: 4, flight: [250.0 * day, 270.0 * day] }],
        bodies: vec![earth.clone(), venus, mars.clone()],
        departure: [0.0, 20.0 * day],
        launch: 2500.0,
        arrival: 0.0,
    };
    let mut raw = sys::myb_chain {
        mu: SUN,
        field: ptr::null_mut(),
        epoch: myb::time_from_seconds(0.0),
        thrust: 1.0,
        exhaust: 30000.0,
        mass0: 1000.0,
        legs: 2,
        departure: [0.0, 20.0 * day],
        launch: 2500.0,
        arrival: 0.0,
        ..Default::default()
    };
    raw.segments[..2].copy_from_slice(&[4, 4]);
    raw.flight[0] = [150.0 * day, 170.0 * day];
    raw.flight[1] = [250.0 * day, 270.0 * day];
    raw.bodies[0] = sys::myb_chain_body { state: Some(raw_earth), context: ptr::null_mut(), gm: 0.0, periapsis: 0.0 };
    raw.bodies[1] = sys::myb_chain_body { state: Some(raw_venus), context: ptr::null_mut(), gm: 3.24859e14, periapsis: 6.5e6 };
    raw.bodies[2] = sys::myb_chain_body { state: Some(raw_mars), context: ptr::null_mut(), gm: 0.0, periapsis: 0.0 };
    raw.impulse_exhaust[1] = 3000.0;
    let size = chain.size()?;
    let guess = chain.guess()?;
    let evaluation = chain.evaluate(&guess)?;
    let (mut unknowns, mut equalities, mut inequalities) = (0usize, 0usize, 0usize);
    // SAFETY: the calls of the C interface with valid arguments
    unsafe {
        assert_eq!(sys::myb_chain_size(&raw, &mut unknowns, &mut equalities, &mut inequalities), sys::MYB_OK);
        assert_eq!((size.unknowns, size.equalities, size.inequalities), (unknowns, equalities, inequalities));
        let mut x = vec![0.0; unknowns];
        assert_eq!(sys::myb_chain_guess(&raw, x.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(bits(&guess), bits(&x));
        let m = equalities + inequalities;
        let (mut propellant, mut constraints, mut gradient, mut jacobian) = (0.0, vec![0.0; m], vec![0.0; unknowns], vec![0.0; m * unknowns]);
        assert_eq!(
            sys::myb_chain_evaluate(&raw, x.as_ptr(), &mut propellant, constraints.as_mut_ptr(), gradient.as_mut_ptr(), jacobian.as_mut_ptr(), ptr::null()),
            sys::MYB_OK
        );
        assert_eq!(evaluation.propellant.to_bits(), propellant.to_bits());
        assert_eq!((bits(&evaluation.constraints), bits(&evaluation.gradient)), (bits(&constraints), bits(&gradient)));
        assert_eq!(bits(&evaluation.jacobian), bits(&jacobian));
    }
    // a powered flyby
    let flyby = myb::chain_flyby(3.24859e14, 3000.0, 3500.0, 7.0e6)?;
    let (mut dv, mut turn, mut partials) = (0.0, 0.0, [0.0; 6]);
    // SAFETY: valid places for the outputs
    assert_eq!(unsafe { sys::myb_chain_flyby(3.24859e14, 3000.0, 3500.0, 7.0e6, &mut dv, &mut turn, partials.as_mut_ptr()) }, sys::MYB_OK);
    assert_eq!((flyby.dv.to_bits(), flyby.turn.to_bits(), bits(&flyby.partials)), (dv.to_bits(), turn.to_bits(), bits(&partials)));
    assert_eq!(myb::chain_flyby(3.24859e14, -1.0, 3500.0, 7.0e6).unwrap_err().status, Status::Domain);

    // failures: a body short of the legs; unknowns of another number; a failing provider
    assert_eq!(chain.evaluate(&guess[1..]).unwrap_err().status, Status::Domain);
    assert_eq!(chain.solve(None, &mut [0.0; 3]).unwrap_err().status, Status::Domain);
    chain.bodies.pop();
    assert_eq!(chain.size().unwrap_err().status, Status::Domain);
    chain.bodies.push(ChainBody { state: Rc::new(|_t: myb::Time| -> myb::Result<State> { Err(myb::Error::new(Status::Range)) }), ..mars });
    assert_eq!(chain.guess().unwrap_err().status, Status::Range);
    let _ = earth;
    Ok(())
}

// ---- Burn arcs and the Q-law ----------------------------------------------------------------------

#[test]
fn thrust_and_the_qlaw_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let arc = ThrustArc {
        start: 100.0,
        stop: 2000.0,
        thrust: 0.2,
        exhaust: 20000.0,
        steering: Steering::Rtn,
        direction: [0.1, 1.0, 0.2],
        ..Default::default()
    };
    let thrust = Thrust::new(500.0, &[arc])?;
    let raw_arc = sys::myb_thrust_arc {
        start: 100.0,
        stop: 2000.0,
        thrust: 0.2,
        exhaust: 20000.0,
        steering: sys::MYB_STEER_RTN,
        direction: [0.1, 1.0, 0.2],
        edelbaum_v0: 0.0,
        edelbaum_beta0: 0.0,
        edelbaum_sense: 0.0,
    };
    let mut raw = sys::myb_thrust::default();
    let y = on_circle(GM, 7.0e6, 0.4, 0.3);
    let (x, v) = ([y[0], y[1], y[2]], [y[3], y[4], y[5]]);
    let (mut a, mut rtn) = ([0.0; 3], [0.0; 3]);
    // SAFETY: the calls of the C interface with valid arguments
    unsafe {
        assert_eq!(sys::myb_thrust_init(&mut raw, 500.0, &raw_arc, 1), sys::MYB_OK);
        same(thrust.as_raw(), &raw);
        assert_eq!(thrust.mass(1000.0).to_bits(), sys::myb_thrust_mass(&raw, 1000.0).to_bits());
        assert_eq!(thrust.mass_at(1000.0).to_bits(), sys::myb_thrust_mass(&raw, 1000.0).to_bits());
        let context = (&mut raw as *mut sys::myb_thrust).cast();
        assert_eq!(sys::myb_thrust_acceleration(context, 1000.0, x.as_ptr(), v.as_ptr(), a.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_thrust_rtn(&raw, 1000.0, x.as_ptr(), v.as_ptr(), rtn.as_mut_ptr()), sys::MYB_OK);
    }
    assert_eq!(bits(&thrust.acceleration(1000.0, &x, &v)?), bits(&a));
    assert_eq!(bits(&thrust.rtn(1000.0, &x, &v)?), bits(&rtn));
    let model = ForceModel::new()?;
    model.add_thrust("engine", &thrust)?;
    let raw_engine = raw_model();
    let mut in_model = [0.0; 3];
    // SAFETY: the calls of the C interface; the thrust outlives the model
    unsafe {
        let context = (&mut raw as *mut sys::myb_thrust).cast();
        assert_eq!(sys::myb_force_model_add(raw_engine, c"engine".as_ptr(), Some(sys::myb_thrust_acceleration), None, context, ptr::null_mut()), sys::MYB_OK);
        assert_eq!(sys::myb_force_model_acceleration(raw_engine, 1000.0, x.as_ptr(), v.as_ptr(), in_model.as_mut_ptr()), sys::MYB_OK);
        sys::myb_force_model_destroy(raw_engine);
    }
    assert_eq!(bits(&model.acceleration(1000.0, &x, &v)?), bits(&in_model));

    // Edelbaum and the Gauss equations
    let (dv, beta0) = myb::edelbaum(7.5e3, 7.0e3, 0.1)?;
    let (speed, di) = myb::edelbaum_state(7.5e3, beta0, 100.0);
    let equinoctial = myb::Equinoctial { p: 7.0e6, f: 0.01, g: -0.02, h: 0.1, k: 0.05, l: 1.0 };
    let rates = myb::equinoctial_rates(GM, &equinoctial, &[1e-4, 2e-4, -1e-4]);
    let (mut raw_dv, mut raw_beta, mut raw_speed, mut raw_di, mut raw_rates) = (0.0, 0.0, 0.0, 0.0, [0.0; 6]);
    // SAFETY: valid places for the outputs
    unsafe {
        assert_eq!(sys::myb_edelbaum(7.5e3, 7.0e3, 0.1, &mut raw_dv, &mut raw_beta), sys::MYB_OK);
        sys::myb_edelbaum_state(7.5e3, raw_beta, 100.0, &mut raw_speed, &mut raw_di);
        sys::myb_equinoctial_rates(GM, &equinoctial, [1e-4, 2e-4, -1e-4].as_ptr(), raw_rates.as_mut_ptr());
    }
    assert_eq!((dv.to_bits(), beta0.to_bits(), speed.to_bits(), di.to_bits()), (raw_dv.to_bits(), raw_beta.to_bits(), raw_speed.to_bits(), raw_di.to_bits()));
    assert_eq!(bits(&rates), bits(&raw_rates));

    // the Q-law
    let law = Qlaw { mu: GM, target: [7.5e6, 0.001, 0.35, 0.0, 0.0], weight: [1.0, 1.0, 1.0, 0.0, 0.0], rp_min: 6.6e6, rp_weight: 0.0 };
    let now = myb::elements_from_state(GM, &y)?;
    let direction = myb::qlaw_direction(&law, &now)?;
    let (mut d, mut q, mut q_rate, mut effectivity) = ([0.0; 3], 0.0, 0.0, 0.0);
    let tolerance = [1e3, 1e-3, 1e-3, 1.0, 1.0];
    // SAFETY: valid places for the inputs and the outputs
    unsafe {
        assert_eq!(sys::myb_qlaw_direction(&law, &now, d.as_mut_ptr(), &mut q, &mut q_rate, &mut effectivity), sys::MYB_OK);
        assert_eq!(myb::qlaw_arrived(&law, &now, &tolerance), sys::myb_qlaw_arrived(&law, &now, tolerance.as_ptr()) != 0);
    }
    assert_eq!(bits(&direction.direction), bits(&d));
    assert_eq!((direction.q.to_bits(), direction.q_rate.to_bits(), direction.effectivity.to_bits()), (q.to_bits(), q_rate.to_bits(), effectivity.to_bits()));
    let term = QlawTerm { law, thrust: 1.0, exhaust: 20000.0, mass0: 500.0, t0: 0.0, tolerance, effectivity: 0.0 };
    let mut raw_term = sys::myb_qlaw_term { law: &law, thrust: 1.0, exhaust: 20000.0, mass0: 500.0, t0: 0.0, tolerance, effectivity: 0.0 };
    let mut raw_a = [0.0; 3];
    // SAFETY: the term and its law live for the call
    assert_eq!(
        unsafe { sys::myb_term_qlaw((&mut raw_term as *mut sys::myb_qlaw_term).cast(), 10.0, x.as_ptr(), v.as_ptr(), raw_a.as_mut_ptr()) },
        sys::MYB_OK
    );
    assert_eq!(bits(&term.acceleration(10.0, &x, &v)?), bits(&raw_a));
    let model = ForceModel::new()?;
    model.add_qlaw("qlaw", &term)?;
    let qlaw_model = model.acceleration(10.0, &x, &v)?;
    let raw_qlaw = raw_model();
    // SAFETY: the term outlives the model
    unsafe {
        assert_eq!(sys::myb_force_model_add(raw_qlaw, c"qlaw".as_ptr(), Some(sys::myb_term_qlaw), None, (&mut raw_term as *mut sys::myb_qlaw_term).cast(), ptr::null_mut()), sys::MYB_OK);
        assert_eq!(sys::myb_force_model_acceleration(raw_qlaw, 10.0, x.as_ptr(), v.as_ptr(), raw_a.as_mut_ptr()), sys::MYB_OK);
        sys::myb_force_model_destroy(raw_qlaw);
    }
    assert_eq!(bits(&qlaw_model), bits(&raw_a));

    // the flight with the mass in the state
    let flight = QlawFlight { engine: term, field: None, rtol: 0.0, atol: 0.0 };
    let y7 = [y[0], y[1], y[2], y[3], y[4], y[5], 500.0];
    let rhs = flight.rhs(0.0, &y7)?;
    let arrival = flight.fly(0.0, &y7, 1200.0, 300.0)?;
    let raw_flight = sys::myb_qlaw_flight { engine: &raw_term, field: ptr::null_mut(), rtol: 0.0, atol: 0.0 };
    let (mut raw_rhs, mut raw_arrival) = ([0.0; 7], sys::myb_qlaw_arrival::default());
    // SAFETY: the flight, its engine and its law live for the calls
    unsafe {
        assert_eq!(sys::myb_qlaw_flight_rhs((&raw_flight as *const sys::myb_qlaw_flight as *mut sys::myb_qlaw_flight).cast(), 0.0, y7.as_ptr(), raw_rhs.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(sys::myb_qlaw_fly(&raw_flight, 0.0, y7.as_ptr(), 1200.0, 300.0, ptr::null(), &mut raw_arrival), sys::MYB_OK);
    }
    assert_eq!(bits(&rhs), bits(&raw_rhs));
    same(&arrival, &raw_arrival);

    // failures
    assert_eq!(Thrust::new(0.0, &[arc]).err().unwrap().status, Status::Domain);
    assert_eq!(myb::edelbaum(-1.0, 7.0e3, 0.1).unwrap_err().status, Status::Domain);
    let hyperbola = myb::Elements { p: 7.0e6, e: 1.5, ..now };
    assert_eq!(myb::qlaw_direction(&law, &hyperbola).unwrap_err().status, Status::Domain);
    assert_eq!(flight.fly(0.0, &y7, 0.0, 300.0).unwrap_err().status, Status::Domain);
    Ok(())
}
