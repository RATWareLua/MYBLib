// SPDX-License-Identifier: BSD-3-Clause
//! The array forms (`*_many` of `myb/time.h`, `myb/kin.h`, `myb/frames.h`, `myb/kepler.h`)
//! over slices against the single C calls in a loop, bit for bit; a failing element is an
//! `Err` with its status and index.
use std::fmt::Debug;
use std::path::PathBuf;

use myb::sys;
use myb::{Elements, Ephemeris, FixedOrientation, KernelPool, Rotation, State, Status, Time, TimeSystem};

const GM: f64 = 3.986004418e14;

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn same<T: Debug>(a: &T, b: &T) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

fn data(relative: &str) -> PathBuf {
    std::env::var_os("MYB_DATA_DIR").map_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../data"), PathBuf::from).join(relative)
}

fn test_data(relative: &str) -> PathBuf {
    std::env::var_os("MYB_TEST_DATA_DIR")
        .map_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../tests/data"), PathBuf::from)
        .join(relative)
}

fn pv(y: &State) -> sys::myb_pv {
    sys::myb_pv { r: [y[0], y[1], y[2]], v: [y[3], y[4], y[5]] }
}

#[test]
fn time_forms_give_the_bits_of_the_single_calls() -> myb::Result<()> {
    let seconds = [0.0, 1.5, -86400.25, 8.0e8 + 0.125];
    let t = myb::time_from_seconds_many(&seconds);
    // SAFETY: functions of numbers
    let single: Vec<Time> = seconds.iter().map(|&s| unsafe { sys::myb_time_from_seconds(s) }).collect();
    same(&t, &single);
    assert_eq!(bits(&myb::time_seconds_many(&t)), bits(&seconds));
    assert!(myb::time_from_seconds_many(&[]).is_empty());
    let labels = myb::time_to_calendar_many(&t);
    same(&labels, &t.iter().map(|&x| myb::time_to_calendar(x)).collect::<Vec<_>>());
    same(&myb::calendar_to_time_many(&labels)?, &labels.iter().map(|c| myb::calendar_to_time(c).unwrap()).collect::<Vec<_>>());
    let mut wrong = labels.clone();
    wrong[2].month = 13;
    let error = myb::calendar_to_time_many(&wrong).unwrap_err();
    assert_eq!((error.status, error.message.as_str()), (Status::Domain, "element 2"));

    let mut pool = KernelPool::new()?;
    for kernel in ["time/time-scales.tk", "time/tdb-tt-fb1990.tk", "time/utc-tai.tk"] {
        pool.load_file(data(kernel))?;
    }
    let system = TimeSystem::new(&pool)?;
    let tdb = system.convert_many(&t, "TT", "TDB")?;
    let utc = system.to_utc_many(&t, "TT")?;
    let again = system.from_utc_many(&utc, "TT")?;
    for k in 0..t.len() {
        same(&tdb[k], &system.convert(t[k], "TT", "TDB")?);
        same(&utc[k], &system.to_utc(t[k], "TT")?);
        same(&again[k], &system.from_utc(&utc[k], "TT")?);
    }
    assert_eq!(system.convert_many(&t, "TT", "NO SUCH").unwrap_err().status, Status::NotFound);
    let mut bad = utc.clone();
    bad[1].second = 75.0;
    assert_eq!(system.from_utc_many(&bad, "TT").unwrap_err().message, "element 1");
    Ok(())
}

#[test]
fn ephemeris_and_frame_forms_give_the_bits_of_the_single_calls() -> myb::Result<()> {
    let mut ephemeris = Ephemeris::new()?;
    ephemeris.load_file(test_data("spk/t02_03.bsp"))?;
    let first = ephemeris.segments()[0];
    let t: Vec<Time> = (0..5).map(|k| myb::time_from_seconds(first.start + (first.stop - first.start) * (k as f64 + 0.5) / 5.0)).collect();
    let (states, frames) = ephemeris.state_many(first.target, first.center, &t)?;
    for k in 0..t.len() {
        assert_eq!(bits(&states[k]), bits(&ephemeris.state(first.target, first.center, t[k])?));
        assert_eq!(frames[k], first.frame);
    }
    let mut late = t.clone();
    late[3] = myb::time_from_seconds(first.stop + 1.0e9);
    let error = ephemeris.state_many(first.target, first.center, &late).unwrap_err();
    assert_eq!((error.status, error.message.as_str()), (Status::Range, "element 3"));

    // rotations applied element by element
    let rotations: Vec<Rotation> = (0..4)
        .map(|k| {
            let mut r = Rotation::default();
            let angles = [0.1 * k as f64, 0.2, 0.3];
            // SAFETY: arrays of three numbers
            unsafe { sys::myb_rotation_euler313(angles.as_ptr(), [1e-5, 0.0, 7e-5].as_ptr(), &mut r) };
            r
        })
        .collect();
    let positions: Vec<[f64; 3]> = (0..4).map(|k| [1.0e7 + k as f64, 2.0e6, -3.0e5]).collect();
    let moving: Vec<State> = (0..4).map(|k| [1.0e7, 2.0e6 * k as f64, -3.0e5, 10.0, 7.0e3, 1.0]).collect();
    let turned = myb::rotation_apply_position_many(&rotations, &positions)?;
    let carried = myb::rotation_apply_state_many(&rotations, &moving)?;
    for k in 0..4 {
        let (mut p, mut s) = ([0.0; 3], [0.0; 6]);
        // SAFETY: arrays of three and six numbers
        unsafe {
            sys::myb_rotation_apply_position(&rotations[k], positions[k].as_ptr(), p.as_mut_ptr());
            sys::myb_rotation_apply_state(&rotations[k], moving[k].as_ptr(), s.as_mut_ptr());
        }
        assert_eq!((bits(&turned[k]), bits(&carried[k])), (bits(&p), bits(&s)));
    }
    assert_eq!(myb::rotation_apply_state_many(&rotations[..3], &moving).unwrap_err().status, Status::Domain);

    // an orientation provider at many instants
    let fixed = FixedOrientation(rotations[2]);
    let orientations = myb::orientation_many(&fixed, &t)?;
    for k in 0..t.len() {
        same(&orientations[k], &myb::orientation_from(&fixed, t[k])?);
    }
    let failing = |x: Time| -> myb::Result<myb::Orientation> {
        if myb::seconds(x) > myb::seconds(t[1]) {
            Err(myb::Error::new(Status::Range))
        } else {
            Ok(myb::Orientation::default())
        }
    };
    assert_eq!(myb::orientation_many(&failing, &t).unwrap_err().message, "element 2");
    Ok(())
}

#[test]
fn kepler_forms_give_the_bits_of_the_single_calls() -> myb::Result<()> {
    let elements: Vec<Elements> = (0..5).map(|k| Elements { p: 7.0e6 + 1.0e5 * k as f64, e: 0.01 * k as f64, i: 0.9, raan: 0.3, argp: 0.2, nu: 0.5 * k as f64 }).collect();
    let states = myb::state_from_elements_many(GM, &elements)?;
    let back = myb::elements_from_state_many(GM, &states)?;
    let equinoctial = myb::equinoctial_from_state_many(GM, &states)?;
    let again = myb::state_from_equinoctial_many(GM, &equinoctial)?;
    for k in 0..elements.len() {
        assert_eq!(bits(&states[k]), bits(&myb::state_from_elements(GM, &elements[k])?));
        same(&back[k], &myb::elements_from_state(GM, &states[k])?);
        let (mut eq, mut pv_out) = (myb::Equinoctial::default(), sys::myb_pv::default());
        // SAFETY: valid places for the inputs and the outputs
        unsafe {
            assert_eq!(sys::myb_equinoctial_from_pv(GM, &pv(&states[k]), &mut eq), sys::MYB_OK);
            assert_eq!(sys::myb_pv_from_equinoctial(GM, &eq, &mut pv_out), sys::MYB_OK);
        }
        same(&equinoctial[k], &eq);
        assert_eq!(bits(&again[k]), bits(&[pv_out.r, pv_out.v].concat()));
    }
    let dt = [0.0, 60.0, -600.0, 5400.0];
    let propagated = myb::kepler_propagate_many(GM, &states[1], &dt)?;
    let with_stm = myb::kepler_stm_many(GM, &states[1], &dt)?;
    for k in 0..dt.len() {
        assert_eq!(bits(&propagated[k]), bits(&myb::kepler_propagate(GM, &states[1], dt[k])?));
        let (y, stm) = myb::kepler_stm(GM, &states[1], dt[k])?;
        assert_eq!((bits(&with_stm[k].0), bits(&with_stm[k].1)), (bits(&y), bits(&stm)));
    }
    let mut wrong = elements.clone();
    wrong[4].p = -1.0;
    let error = myb::state_from_elements_many(GM, &wrong).unwrap_err();
    assert_eq!((error.status, error.message.as_str()), (Status::Domain, "element 4"));
    assert_eq!(myb::kepler_propagate_many(-1.0, &states[1], &dt).unwrap_err().message, "element 0");
    Ok(())
}
