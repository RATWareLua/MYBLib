// SPDX-License-Identifier: BSD-3-Clause
//! Estimation (`myb/estim.h`) through the safe types against the same C calls, bit for bit;
//! failures as `Err` with the status of the library.
use std::ffi::c_void;
use std::fmt::Debug;
use std::ptr;
use std::rc::Rc;

use myb::sys;
use myb::{
    Attitude, Body, Ellipsoid, Environment, Error, Filter, FilterOptions, FitOptions, ForceModel, LightMeasure, Measurement, Observation,
    State, Station, Status,
};

const GM: f64 = 3.986004418e14;
const RADIUS: f64 = 6378137.0;

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn same<T: Debug>(a: &T, b: &T) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

fn orbit() -> State {
    let elements = myb::Elements { p: 7.0e6, e: 0.01, i: 0.9, raan: 0.3, argp: 0.2, nu: 0.1 };
    myb::state_from_elements(GM, &elements).unwrap()
}

fn ellipsoid() -> Ellipsoid {
    Ellipsoid { a: RADIUS, f: 1.0 / 298.257223563 }
}

type RawMeasure = unsafe extern "C" fn(*mut c_void, f64, *const f64, *mut f64, *mut f64) -> sys::myb_status;

fn raw_measure(function: RawMeasure, context: *mut c_void, count: usize, t: f64, y: &State) -> (Vec<f64>, Vec<f64>) {
    let (mut values, mut partials) = (vec![0.0; count], vec![0.0; 6 * count]);
    // SAFETY: the context is what the function asks for; the outputs hold count and 6 count numbers
    assert_eq!(unsafe { function(context, t, y.as_ptr(), values.as_mut_ptr(), partials.as_mut_ptr()) }, sys::MYB_OK);
    (values, partials)
}

#[test]
fn measurements_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let env = Rc::new(Environment::new(myb::time_from_seconds(0.0))?);
    env.add_origin(&Body { id: 399, gm: GM, ..Default::default() })?;
    let station = Rc::new(Station::new(&env, 0, &ellipsoid(), 0.3, 0.8, 100.0)?);
    let mut raw = sys::myb_station::default();
    // SAFETY: the environment is alive for the test
    assert_eq!(unsafe { sys::myb_station_init(&mut raw, env.as_ptr(), 0, &ellipsoid(), 0.3, 0.8, 100.0) }, sys::MYB_OK);
    assert_eq!(bits(&station.position()), bits(&raw.position));
    assert_eq!(bits(&station.enu()), bits(&raw.enu));
    let mut there = [0.0; 6];
    // SAFETY: the station is valid
    assert_eq!(unsafe { sys::myb_station_state(&raw, 120.0, there.as_mut_ptr()) }, sys::MYB_OK);
    assert_eq!(bits(&station.state(120.0)?), bits(&there));

    let (t, y) = (120.0, orbit());
    let context = &mut raw as *mut sys::myb_station as *mut c_void;
    let cases: [(Measurement, RawMeasure, *mut c_void, usize); 6] = [
        (Measurement::range(&station), sys::myb_measure_range, context, 1),
        (Measurement::range_rate(&station), sys::myb_measure_range_rate, context, 1),
        (Measurement::radec(&station), sys::myb_measure_radec, context, 2),
        (Measurement::azel(&station), sys::myb_measure_azel, context, 2),
        (Measurement::position(), sys::myb_measure_position, ptr::null_mut(), 3),
        (Measurement::state(), sys::myb_measure_state, ptr::null_mut(), 6),
    ];
    for (measurement, function, context, count) in cases {
        assert_eq!(measurement.count(), count);
        let (values, partials) = measurement.evaluate(t, &y)?;
        let (raw_values, raw_partials) = raw_measure(function, context, count, t, &y);
        assert_eq!((bits(&values), bits(&partials)), (bits(&raw_values), bits(&raw_partials)));
        assert_eq!(bits(&measurement.value(t, &y)?), bits(&raw_values));
    }

    // with the light time
    let light = LightMeasure { station: Rc::clone(&station), c: 299792458.0, gm: GM, two_way: true, aberration: false };
    let mut raw_light = sys::myb_light_measure { station: &raw, c: 299792458.0, gm: GM, two_way: 1, aberration: 0 };
    let light_context = &mut raw_light as *mut sys::myb_light_measure as *mut c_void;
    let cases: [(Measurement, RawMeasure, usize); 4] = [
        (Measurement::light_range(&light), sys::myb_measure_light_range, 1),
        (Measurement::light_range_rate(&light), sys::myb_measure_light_range_rate, 1),
        (Measurement::light_radec(&light), sys::myb_measure_light_radec, 2),
        (Measurement::light_azel(&light), sys::myb_measure_light_azel, 2),
    ];
    for (measurement, function, count) in cases {
        let (values, partials) = measurement.evaluate(t, &y)?;
        let (raw_values, raw_partials) = raw_measure(function, light_context, count, t, &y);
        assert_eq!((bits(&values), bits(&partials)), (bits(&raw_values), bits(&raw_partials)));
    }

    // a point off the centre of mass, turned by a fixed attitude
    let matrix = [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
    let offset = Measurement::offset(&Measurement::range(&station), Attitude::fixed(matrix), [1.0, 2.0, 3.0], GM);
    let mut raw_offset = sys::myb_offset_measure {
        measure: Some(sys::myb_measure_range),
        context,
        attitude: Some(sys::myb_attitude_fixed),
        attitude_context: matrix.as_ptr() as *mut c_void,
        offset: [1.0, 2.0, 3.0],
        gm: GM,
    };
    let (values, partials) = offset.evaluate(t, &y)?;
    let (raw_values, raw_partials) = raw_measure(sys::myb_measure_offset, (&mut raw_offset as *mut sys::myb_offset_measure).cast(), 1, t, &y);
    assert_eq!((bits(&values), bits(&partials)), (bits(&raw_values), bits(&raw_partials)));

    // a closure
    let custom = Measurement::custom(1, |_t, y, values, partials| {
        values[0] = y[0];
        if let Some(p) = partials {
            p[0] = 1.0;
        }
        Ok(())
    });
    let (values, partials) = custom.evaluate(t, &y)?;
    assert_eq!((values[0], partials[0], partials[1]), (y[0], 1.0, 0.0));
    let failing = Measurement::custom(1, |_, _, _, _| Err(Error::new(Status::Range)));
    assert_eq!(failing.value(t, &y).unwrap_err().status, Status::Range);
    // a station on a body the environment does not have; a station at the object
    assert_eq!(Station::new(&env, 3, &ellipsoid(), 0.0, 0.0, 0.0).err().unwrap().status, Status::Domain);
    let at_station = station.state(t)?;
    assert_eq!(Measurement::range(&station).value(t, &at_station).unwrap_err().status, Status::Domain);
    Ok(())
}

fn observations(measurement: &Measurement) -> Vec<Observation> {
    let truth = orbit();
    (1..=12)
        .map(|k| {
            let t = 60.0 * k as f64;
            let y = myb::kepler_propagate(GM, &truth, t).unwrap();
            Observation::new(t, measurement, &y[..3], &[1.0, 1.0, 1.0])
        })
        .collect()
}

fn raw_observations(observations: &[Observation]) -> Vec<sys::myb_observation> {
    observations
        .iter()
        .map(|o| {
            let mut made = sys::myb_observation { t: o.t, measure: Some(sys::myb_measure_position), count: 3, ..Default::default() };
            made.value[..3].copy_from_slice(&o.value);
            made.sigma[..3].copy_from_slice(&o.sigma);
            made
        })
        .collect()
}

fn guess() -> State {
    let mut y = orbit();
    for (k, d) in [100.0, -50.0, 20.0, 0.1, -0.05, 0.02].iter().enumerate() {
        y[k] += d;
    }
    y
}

#[test]
fn a_fit_and_a_filter_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let model = Rc::new(ForceModel::new()?);
    model.add_point_mass("earth", GM)?;
    let observed = observations(&Measurement::position());
    let fit = myb::fit_orbit(&model, 0.0, &guess(), &observed, &FitOptions::default())?;
    // exact positions: the weighted RMS falls to rounding (and keeps changing relative to itself)
    assert!(fit.result.rms < 1.0e-6 && fit.result.rms_before > 10.0, "{:?}", fit.result);

    let raws = raw_observations(&observed);
    let mut y = guess();
    let mut result = sys::myb_fit_result::default();
    let (mut covariance, mut residuals) = (vec![0.0; 36], vec![0.0; 36]);
    let options = sys::myb_fit_options::default();
    // SAFETY: the observations and outputs have the sizes the library asks for
    assert_eq!(
        unsafe {
            sys::myb_fit_orbit(
                model.as_ptr(),
                0.0,
                y.as_mut_ptr(),
                raws.as_ptr(),
                raws.len(),
                &options,
                ptr::null(),
                &mut result,
                covariance.as_mut_ptr(),
                residuals.as_mut_ptr(),
            )
        },
        sys::MYB_OK
    );
    assert_eq!(bits(&fit.y0), bits(&y));
    assert_eq!(bits(&fit.covariance), bits(&covariance));
    assert_eq!(bits(&fit.residuals), bits(&residuals));
    same(&fit.result, &result);

    // the filter and its smoother
    let mut p0 = vec![0.0; 36];
    for k in 0..6 {
        p0[7 * k] = if k < 3 { 1.0e4 } else { 1.0 };
    }
    let options = FilterOptions { history: 20, ..Default::default() };
    let mut filter = Filter::new(&model, &options)?;
    filter.start(0.0, &guess(), &p0, 800.0)?;
    let innovations: Vec<Vec<f64>> = observed.iter().map(|o| filter.update(o)).collect::<myb::Result<_>>()?;
    filter.predict(780.0)?;
    let estimate = filter.get();
    let steps = filter.steps();
    filter.smooth()?;
    let first = filter.step(0).expect("the first instant");
    assert!(filter.step(steps).is_none());
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    unsafe {
        let made = sys::myb_filter_options { history: 20, ..Default::default() };
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_filter_create(model.as_ptr(), &made, ptr::null(), &mut raw), sys::MYB_OK);
        assert_eq!(sys::myb_filter_start(raw, 0.0, guess().as_ptr(), p0.as_ptr(), 800.0), sys::MYB_OK);
        for (observation, innovation) in raws.iter().zip(&innovations) {
            let mut raw_innovation = [0.0; 3];
            assert_eq!(sys::myb_filter_update(raw, observation, raw_innovation.as_mut_ptr()), sys::MYB_OK);
            assert_eq!(bits(innovation), bits(&raw_innovation));
        }
        assert_eq!(sys::myb_filter_predict(raw, 780.0), sys::MYB_OK);
        let (mut t, mut y, mut covariance) = (0.0, [0.0; 6], vec![0.0; 36]);
        sys::myb_filter_get(raw, &mut t, y.as_mut_ptr(), covariance.as_mut_ptr());
        assert_eq!((t.to_bits(), bits(&estimate.y), bits(&estimate.covariance)), (estimate.t.to_bits(), bits(&y), bits(&covariance)));
        assert_eq!(sys::myb_filter_steps(raw), steps);
        assert_eq!(sys::myb_filter_smooth(raw), sys::MYB_OK);
        let (mut unknowns, mut smoothed) = (vec![0.0; 6], vec![0.0; 36]);
        sys::myb_filter_step(raw, 0, &mut t, unknowns.as_mut_ptr(), smoothed.as_mut_ptr());
        assert_eq!((bits(&first.unknowns), bits(&first.covariance)), (bits(&unknowns), bits(&smoothed)));
        sys::myb_filter_destroy(raw);
    }
    assert_eq!(filter.bias(), Vec::<f64>::new());

    // failures
    assert_eq!(myb::fit_orbit(&model, 0.0, &guess(), &[], &FitOptions::default()).unwrap_err().status, Status::Domain);
    let uneven = Observation::new(60.0, &Measurement::position(), &[1.0, 2.0], &[1.0, 1.0]);
    assert_eq!(myb::fit_orbit(&model, 0.0, &guess(), std::slice::from_ref(&uneven), &FitOptions::default()).unwrap_err().status, Status::Domain);
    let fixed = FitOptions { fixed: Some(vec![false; 5]), ..Default::default() };
    assert_eq!(myb::fit_orbit(&model, 0.0, &guess(), &observed, &fixed).unwrap_err().status, Status::Domain);
    let failing = Measurement::custom(3, |_, _, _, _| Err(Error::new(Status::Range)));
    let refused = observations(&failing);
    assert_eq!(myb::fit_orbit(&model, 0.0, &guess(), &refused, &FitOptions::default()).unwrap_err().status, Status::Range);
    let mut filter = Filter::new(&model, &FilterOptions::default())?;
    assert_eq!(filter.start(0.0, &guess(), &p0[..35], 800.0).unwrap_err().status, Status::Domain);
    filter.start(0.0, &guess(), &p0, 800.0)?;
    assert_eq!(filter.update(&uneven).unwrap_err().status, Status::Domain);
    assert_eq!(filter.predict(900.0).unwrap_err().status, Status::Domain);
    assert_eq!(filter.set_bias(&[1.0]).unwrap_err().status, Status::Domain);
    assert_eq!(filter.smooth().unwrap_err().status, Status::Domain);
    assert_eq!(Filter::new(&model, &FilterOptions { process_noise: -1.0, ..Default::default() }).err().unwrap().status, Status::Domain);
    Ok(())
}
