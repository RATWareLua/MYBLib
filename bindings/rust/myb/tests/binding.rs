// SPDX-License-Identifier: BSD-3-Clause
//! The safe types against the C interface they wrap.
//!
//!  - the low orbit of examples/low_orbit.c built from the types gives the nodes and the
//!    final state of the same calls made through `myb::sys` — the C functions — bit for bit;
//!  - a failed call is an `Err` with the status of the library and the message of a parser;
//!  - terms, providers and events as closures; one that fails stops the call with its
//!    status;
//!  - a model keeps alive the environment and the field it was given, an orbit its model;
//!  - kernels and time scales, an SPK set and a body of it in an environment, a message
//!    written and read in both notations;
//!  - a panic in a closure that C calls ends the process: it is not caught and does not
//!    unwind through the frames of C (the test runs itself again as the process that dies).
use std::cell::Cell;
use std::path::PathBuf;
use std::ptr;
use std::rc::Rc;

use myb::sys;
use myb::{Body, Direction, Environment, Error, ForceModel, Orbit, OrbitOptions, State, Status};

const GM: f64 = 3.986004418e14;
const RADIUS: f64 = 6378137.0;
const J2: f64 = 1.08263e-3;
const DAY: f64 = 86400.0;

fn repository() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn data(relative: &str) -> PathBuf {
    std::env::var_os("MYB_DATA_DIR").map_or_else(|| repository().join("data"), PathBuf::from).join(relative)
}

fn test_data(relative: &str) -> PathBuf {
    std::env::var_os("MYB_TEST_DATA_DIR")
        .map_or_else(|| repository().join("tests/data"), PathBuf::from)
        .join(relative)
}

fn elements() -> myb::Elements {
    let e = 0.001;
    myb::Elements {
        p: (RADIUS + 400000.0) * (1.0 - e * e),
        e,
        i: 51.6 * std::f64::consts::PI / 180.0,
        raan: 0.5,
        argp: 1.0,
        nu: 0.0,
    }
}

fn coefficients() -> ([f64; 6], [f64; 6]) {
    ([1.0, 0.0, 0.0, -J2 / 5.0_f64.sqrt(), 0.0, 0.0], [0.0; 6])
}

unsafe extern "C" fn equator(_context: *mut std::ffi::c_void, _t: f64, y: *const f64) -> f64 {
    // SAFETY: y is the state the library passes
    unsafe { *y.add(2) }
}

/// The example call for call through the C functions.
fn low_orbit_raw() -> (Vec<f64>, State, i64) {
    let (c, s) = coefficients();
    // SAFETY: the calls of the example, each with valid arguments; everything is destroyed
    unsafe {
        let mut field = ptr::null_mut();
        let mut work = ptr::null_mut();
        assert_eq!(sys::myb_gravity_field_create(GM, RADIUS, 2, c.as_ptr(), s.as_ptr(), ptr::null(), &mut field), sys::MYB_OK);
        assert_eq!(sys::myb_gravity_work_create(2, ptr::null(), &mut work), sys::MYB_OK);
        let options = sys::myb_env_options::default();
        let mut env = ptr::null_mut();
        assert_eq!(sys::myb_env_create(&options, ptr::null(), &mut env), sys::MYB_OK);
        let body = sys::myb_env_body { id: 1, gm: GM, ..Default::default() };
        let mut planet = 0usize;
        assert_eq!(sys::myb_env_add_body(env, &body, &mut planet), sys::MYB_OK);
        let mut gravity = sys::myb_field_term { env, body: planet, field, work, n_max: -1, ..Default::default() };
        let mut air = sys::myb_exponential_atmosphere::default();
        air.ellipsoid.a = RADIUS;
        air.count = 1;
        air.base[0] = 400000.0;
        air.density[0] = 3.0e-12;
        air.scale[0] = 58000.0;
        let object = sys::myb_object { mass: 500.0, area: 4.0, cr: 0.0, cd: 2.2, ..Default::default() };
        let mut drag = sys::myb_drag_term {
            env,
            body: planet,
            density: Some(sys::myb_density_exponential),
            density_context: (&mut air as *mut sys::myb_exponential_atmosphere).cast(),
            object: &object,
            ..Default::default()
        };
        let mut model = ptr::null_mut();
        assert_eq!(sys::myb_force_model_create(ptr::null(), &mut model), sys::MYB_OK);
        assert_eq!(
            sys::myb_force_model_add(model, c"gravity".as_ptr(), Some(sys::myb_term_field), Some(sys::myb_term_field_jacobian), (&mut gravity as *mut sys::myb_field_term).cast(), ptr::null_mut()),
            sys::MYB_OK
        );
        assert_eq!(
            sys::myb_force_model_add(model, c"drag".as_ptr(), Some(sys::myb_term_drag), Some(sys::myb_term_drag_jacobian), (&mut drag as *mut sys::myb_drag_term).cast(), ptr::null_mut()),
            sys::MYB_OK
        );
        let mut pv = sys::myb_pv::default();
        assert_eq!(sys::myb_pv_from_elements(GM, &elements(), &mut pv), sys::MYB_OK);
        let y0 = [pv.r[0], pv.r[1], pv.r[2], pv.v[0], pv.v[1], pv.v[2]];
        let mut orbit = ptr::null_mut();
        assert_eq!(sys::myb_orbit_create(model, ptr::null(), ptr::null(), &mut orbit), sys::MYB_OK);
        assert_eq!(sys::myb_orbit_start(orbit, 0.0, y0.as_ptr(), DAY), sys::MYB_OK);
        let mut nodes = Vec::new();
        for _ in 0..32 {
            let mut found = 0;
            let mut t_node = 0.0;
            assert_eq!(sys::myb_orbit_event(orbit, Some(equator), ptr::null_mut(), DAY, 1, &mut found, &mut t_node, ptr::null_mut()), sys::MYB_OK);
            if found == 0 {
                break;
            }
            nodes.push(t_node);
        }
        let mut last = [0.0; 6];
        assert_eq!(sys::myb_orbit_state(orbit, DAY, last.as_mut_ptr()), sys::MYB_OK);
        let steps = sys::myb_orbit_steps(orbit);
        sys::myb_orbit_destroy(orbit);
        sys::myb_force_model_destroy(model);
        sys::myb_env_destroy(env);
        sys::myb_gravity_work_destroy(work);
        sys::myb_gravity_field_destroy(field);
        (nodes, last, steps)
    }
}

/// The model of the example from the types; what it is built from is dropped.
fn low_orbit_model() -> myb::Result<Rc<ForceModel>> {
    let (c, s) = coefficients();
    let field = Rc::new(myb::GravityField::new(GM, RADIUS, 2, &c, &s)?);
    let env = Rc::new(Environment::new(myb::time_from_seconds(0.0))?);
    let planet = env.add_origin(&Body { id: 1, gm: GM, ..Default::default() })?;
    let model = Rc::new(ForceModel::new()?);
    model.add_field("gravity", &env, planet, &field, 0, None)?;
    let air = myb::ExponentialAtmosphere { equatorial_radius: RADIUS, flattening: 0.0, layers: vec![(400000.0, 3.0e-12, 58000.0)] };
    let object = myb::Object { mass: 500.0, area: 4.0, cr: 0.0, cd: 2.2 };
    model.add_drag("drag", &env, planet, &air, &object)?;
    Ok(model)
}

fn bits(state: &State) -> [u64; 6] {
    state.map(f64::to_bits)
}

#[test]
fn the_types_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let (nodes_raw, last_raw, steps_raw) = low_orbit_raw();
    let model = low_orbit_model()?; // the field and the environment live in the model alone
    assert_eq!(model.term_count(), 2);
    assert_eq!(model.term_name(1).as_deref(), Some("drag"));
    assert_eq!(model.term_name(2), None);
    let mut orbit = Orbit::new(&model, &OrbitOptions::default())?;
    drop(model); // and the model in the orbit
    orbit.start(0.0, &myb::state_from_elements(GM, &elements())?, DAY)?;
    let mut nodes = Vec::new();
    while let Some((t, y)) = orbit.next_event(|_, y| y[2], DAY, Direction::Rising)? {
        assert!(y[2].abs() < 1e-6 && y[5] > 0.0);
        nodes.push(t);
        assert!(nodes.len() <= 32);
    }
    let last = orbit.state(DAY)?;
    assert!(nodes.len() >= 14);
    assert_eq!(nodes.iter().map(|t| t.to_bits()).collect::<Vec<_>>(), nodes_raw.iter().map(|t| t.to_bits()).collect::<Vec<_>>());
    assert_eq!(bits(&last), bits(&last_raw));
    assert_eq!(orbit.steps(), steps_raw);
    assert_eq!(orbit.time(), DAY);
    // the node drifts as the first-order theory says
    let a0 = RADIUS + 400000.0;
    let theory = -1.5 * (GM / a0.powi(3)).sqrt() * J2 * (RADIUS / elements().p).powi(2) * elements().i.cos();
    let after = myb::elements_from_state(GM, &last)?;
    assert!(((after.raan - elements().raan) / DAY / theory - 1.0).abs() < 0.02);
    Ok(())
}

#[test]
fn a_failure_is_an_error_with_the_status() {
    let wrong = myb::Elements { p: -1.0, ..Default::default() };
    let error = myb::state_from_elements(GM, &wrong).unwrap_err();
    assert_eq!(error.status, Status::Domain);
    assert_eq!(error.status.name(), "MYB_E_DOMAIN");
    assert_eq!(error.to_string(), "MYB_E_DOMAIN");
    // a parser names the line
    let error = myb::Oem::parse("CCSDS_OEM_VERS = 2.0\nMETA_START\n").err().unwrap();
    assert_eq!(error.status, Status::Format);
    assert!(error.message.starts_with("line 2"), "{error}");
    let mut pool = myb::KernelPool::new().unwrap();
    assert_eq!(pool.load_file(data("no-such-kernel.tk")).unwrap_err().status, Status::Io);
    assert_eq!(myb::GravityField::new(GM, RADIUS, 4, &[1.0], &[0.0]).err().unwrap().status, Status::Domain);
    // a body that the environment does not have
    let env = Rc::new(Environment::new(myb::time_from_seconds(0.0)).unwrap());
    let model = ForceModel::new().unwrap();
    assert_eq!(model.add_third_body("nobody", &env, 3).unwrap_err().status, Status::Domain);
    assert_eq!(env.state(0.0, 0).unwrap_err().status, Status::Domain);
    assert!(!model.set_enabled(0, false));
}

#[test]
fn terms_and_events_are_closures() -> myb::Result<()> {
    let x = [7.0e6, 1.0e5, -2.0e5];
    let v = [0.0, 7.5e3, 1.0e3];
    let model = Rc::new(ForceModel::new()?);
    model.add_point_mass("planet", GM)?;
    let push = [1.0e-6, -2.0e-6, 3.0e-6];
    let calls = Rc::new(Cell::new(0u32));
    let counted = Rc::clone(&calls);
    let index = model.add_term("push", move |_, _, _| {
        counted.set(counted.get() + 1);
        Ok(push)
    })?;
    let both = model.acceleration(0.0, &x, &v)?;
    assert!(model.set_enabled(index, false));
    let planet = model.acceleration(0.0, &x, &v)?;
    assert!(model.set_enabled(index, true));
    for i in 0..3 {
        assert!((both[i] - planet[i] - push[i]).abs() < 1e-14);
    }
    assert_eq!(calls.get(), 1);
    // an orbit with the push: an event, many states, the transition matrix, an impulse
    let options = OrbitOptions { variational: true, ..Default::default() };
    let mut orbit = Orbit::new(&model, &options)?;
    let y0 = myb::state_from_elements(GM, &elements())?;
    orbit.start(0.0, &y0, 6000.0)?;
    let crossing = orbit.next_event(|_, y| y[2], 6000.0, Direction::Falling)?;
    let (_, at) = crossing.expect("the descending node of the first revolution");
    assert!(at[2].abs() < 1e-6 && at[5] < 0.0);
    let (y, phi) = orbit.partials(5000.0)?;
    assert_eq!(bits(&y), bits(&orbit.state(5000.0)?));
    assert!(phi[0] != 0.0 && phi[35].is_finite());
    let many = orbit.states(&[5100.0, 5200.0, 5300.0])?;
    assert_eq!(many.len(), 3);
    orbit.kick(&[0.0, 1.0, 0.0])?;
    assert!(calls.get() > 100);
    // DOP853 agrees with IAS15 to its tolerance
    let plain = Rc::new(ForceModel::new()?);
    plain.add_point_mass("planet", GM)?;
    let mut ias15 = Orbit::new(&plain, &OrbitOptions::default())?;
    let dop = OrbitOptions { integrator: myb::Integrator::Dop853, rtol: 1e-11, atol: 1e-6, ..Default::default() };
    let mut dop853 = Orbit::new(&plain, &dop)?;
    ias15.start(0.0, &y0, 3000.0)?;
    dop853.start(0.0, &y0, 3000.0)?;
    let (a, b) = (ias15.state(3000.0)?, dop853.state(3000.0)?);
    assert!((0..3).all(|i| (a[i] - b[i]).abs() < 1e-2));
    let kepler = myb::kepler_propagate(GM, &y0, 3000.0)?;
    assert!((0..3).all(|i| (a[i] - kepler[i]).abs() < 1e-6));

    // a term that fails: its status comes back from the evaluation
    let failing = ForceModel::new()?;
    failing.add_term("refuses", |_, _, _| Err(Error::new(Status::Range)))?;
    assert_eq!(failing.acceleration(0.0, &x, &v).unwrap_err().status, Status::Range);
    Ok(())
}

#[test]
fn a_provider_is_a_closure() -> myb::Result<()> {
    let rate = 2.0e-7;
    let distance = 1.5e11;
    let sun = move |seconds: f64| -> myb::Result<State> {
        let angle = rate * seconds;
        Ok([distance * angle.cos(), distance * angle.sin(), 0.0, -distance * rate * angle.sin(), distance * rate * angle.cos(), 0.0])
    };
    let env = Rc::new(Environment::new(myb::time_from_seconds(0.0))?);
    env.add_origin(&Body { id: 399, gm: GM, ..Default::default() })?;
    let index = env.add_with(&Body { id: 10, gm: 1.32712440018e20, ..Default::default() }, sun)?;
    assert_eq!(env.body_count(), 2);
    assert_eq!(env.find(10)?, index);
    assert_eq!(env.find(11).unwrap_err().status, Status::NotFound);
    assert_eq!(bits(&env.state(1000.0, index)?), bits(&sun(1000.0)?));
    let model = ForceModel::new()?;
    model.add_point_mass("planet", GM)?;
    model.add_third_body("sun", &env, index)?;
    drop(env); // the model keeps the environment, and the environment its closure
    let x = [7.0e6, 1.0e5, -2.0e5];
    let v = [0.0, 7.5e3, 1.0e3];
    let with_sun = model.acceleration(0.0, &x, &v)?;
    model.set_enabled(1, false);
    let alone = model.acceleration(0.0, &x, &v)?;
    let tide = (0..3).map(|i| (with_sun[i] - alone[i]).powi(2)).sum::<f64>().sqrt();
    assert!(tide > 1e-7 && tide < 1e-5); // 2 GM r / d^3 at 7000 km is 5e-7 m/s^2
    // a provider that fails stops the evaluation with its status
    let broken = Environment::new(myb::time_from_seconds(0.0))?;
    broken.add_origin(&Body::default())?;
    let lost = broken.add_with(&Body { id: 5, ..Default::default() }, |_| Err(Error::new(Status::Range)))?;
    assert_eq!(broken.state(0.0, lost).unwrap_err().status, Status::Range);
    Ok(())
}

#[test]
fn time_scales_and_an_spk_set() -> myb::Result<()> {
    let mut pool = myb::KernelPool::new()?;
    for kernel in ["time/time-scales.tk", "time/tdb-tt-fb1990.tk", "time/utc-tai.tk"] {
        pool.load_file(data(kernel))?;
    }
    let system = myb::TimeSystem::new(&pool)?;
    drop(pool); // the system keeps nothing of the pool
    let label = myb::Calendar { year: 2026, month: 11, day: 11, hour: 12, minute: 30, second: 15.25 };
    let instant = system.from_utc(&label, "TT")?;
    let back = system.to_utc(instant, "TT")?;
    assert_eq!((back.year, back.month, back.day, back.hour, back.minute), (2026, 11, 11, 12, 30));
    assert!((back.second - 15.25).abs() < 1e-9);
    let tai = system.convert(instant, "TT", "TAI")?;
    assert!((myb::seconds(instant) - myb::seconds(tai) - 32.184).abs() < 1e-6);
    assert_eq!(system.convert(instant, "TT", "NO SUCH SCALE").unwrap_err().status, Status::NotFound);

    let mut loaded = myb::Ephemeris::new()?;
    loaded.load_file(test_data("spk/t02_03.bsp"))?;
    assert!(loaded.load_file(test_data("spk/missing.bsp")).is_err());
    let eph = Rc::new(loaded);
    let first = eph.segments()[0];
    let middle = myb::time_from_seconds(0.5 * (first.start + first.stop));
    let there = eph.state(first.target, first.center, middle)?;
    let mut raw = [0.0; 6];
    // SAFETY: the set is alive; raw holds six numbers
    assert_eq!(unsafe { sys::myb_spk_state(eph.as_ptr(), first.target, first.center, middle, raw.as_mut_ptr(), ptr::null_mut()) }, sys::MYB_OK);
    assert_eq!(bits(&there), bits(&raw));
    assert_eq!(eph.state(first.target, first.center, myb::time_from_seconds(first.stop + 1.0e9)).unwrap_err().status, Status::Range);
    let env = Environment::new(middle)?;
    env.add_origin(&Body { id: i64::from(first.center), gm: GM, ..Default::default() })?;
    let index = env.add_from_ephemeris(&Body { id: i64::from(first.target), gm: 1.0, ..Default::default() }, &eph, first.target, first.center, 0)?;
    drop(eph); // the environment keeps the ephemeris
    assert_eq!(bits(&env.state(0.0, index)?), bits(&there));
    Ok(())
}

#[test]
fn a_message_in_both_notations() -> myb::Result<()> {
    let segment = myb::OemSegment {
        object_name: "A & B".into(),
        object_id: "2026-001A".into(),
        center_name: "EARTH".into(),
        ref_frame: "EME2000".into(),
        time_system: "UTC".into(),
        interpolation: Some(("HERMITE".into(), 7)),
        epochs: vec![
            myb::Calendar { year: 2026, month: 3, day: 1, hour: 12, minute: 0, second: 0.0 },
            myb::Calendar { year: 2026, month: 3, day: 1, hour: 12, minute: 1, second: 0.5 },
        ],
        states: vec![[7.0e6, 1.0, -2.5, 0.25, 7.5e3, -1.0e3], [7.0e6 + 15.0, 4.5e5, -6.0e4, -480.0, 7.4e3, -990.0]],
    };
    let created = myb::Calendar { year: 2026, month: 3, day: 2, hour: 0, minute: 0, second: 0.0 };
    for notation in [myb::Notation::Kvn, myb::Notation::Xml] {
        let text = myb::write_oem("myb", &created, std::slice::from_ref(&segment), notation)?;
        assert_eq!(text.starts_with('<'), notation == myb::Notation::Xml);
        let message = myb::Oem::parse(&text)?;
        assert_eq!((message.segment_count(), message.state_count(0), message.state_count(1)), (1, 2, 0));
        assert_eq!(message.meta(0, "OBJECT_NAME").as_deref(), Some("A & B"));
        assert_eq!(message.meta(0, "INTERPOLATION_DEGREE").as_deref(), Some("7"));
        assert_eq!(message.meta(0, "NO_SUCH_KEY"), None);
        let (epoch, y) = message.state(0, 1).expect("the second line");
        assert_eq!((epoch.minute, epoch.second), (1, 0.5));
        assert_eq!(bits(&y), bits(&segment.states[1]));
        assert!(message.state(0, 2).is_none());
    }
    assert_eq!(myb::write_oem("myb", &created, &[], myb::Notation::Kvn).unwrap_err().status, Status::Domain);
    let uneven = myb::OemSegment { states: vec![[0.0; 6]], ..segment };
    assert_eq!(myb::write_oem("myb", &created, &[uneven], myb::Notation::Kvn).unwrap_err().status, Status::Domain);
    Ok(())
}

/// The name of the variable that makes this test binary the process that panics.
const DIES: &str = "MYB_TEST_PANIC_IN";

fn panic_in(what: &str) {
    let model = ForceModel::new().unwrap();
    model.add_point_mass("centre", GM).unwrap();
    match what {
        "term" => {
            model.add_term("broken", |_t, _x, _v| -> myb::Result<[f64; 3]> { panic!("in a term") }).unwrap();
            let _ = model.acceleration(0.0, &[7.0e6, 0.0, 0.0], &[0.0, 7.5e3, 0.0]);
        }
        _ => {
            let model = Rc::new(model);
            let mut orbit = Orbit::new(&model, &OrbitOptions::default()).unwrap();
            orbit.start(0.0, &[7.0e6, 0.0, 0.0, 0.0, 7.5e3, 1.0e3], 6000.0).unwrap();
            let _ = orbit.next_event(|_t, _y| -> f64 { panic!("in an event") }, 6000.0, Direction::Any);
        }
    }
}

#[test]
fn a_panic_in_a_closure_called_from_c_ends_the_process() {
    if let Ok(what) = std::env::var(DIES) {
        panic_in(&what);
        // not reached: the process is gone. Leaving normally is what the parent takes as failure.
        std::process::exit(0);
    }
    let this = std::env::current_exe().unwrap();
    for what in ["term", "event"] {
        let output = std::process::Command::new(&this)
            .args(["--exact", "a_panic_in_a_closure_called_from_c_ends_the_process", "--nocapture", "--test-threads", "1"])
            .env(DIES, what)
            .output()
            .unwrap();
        let said = String::from_utf8_lossy(&output.stderr);
        // the message of the panic, then the end: not the exit of a test that failed by
        // unwinding (101), not success
        assert!(said.contains(&format!("in {} {}", if what == "term" { "a" } else { "an" }, what)), "{what}: {said}");
        assert!(!output.status.success(), "{what}: the process went on");
        assert_ne!(output.status.code(), Some(101), "{what}: the panic unwound through C: {said}");
        assert!(said.contains("panic in a function that cannot unwind"), "{what}: {said}");
    }
}
