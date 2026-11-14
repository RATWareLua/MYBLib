// SPDX-License-Identifier: BSD-3-Clause
//! Reference frames (`myb/frames.h`): the models of the Earth, EOP, WGCCRE, binary PCK and the
//! graph of frames through the safe types against the same C calls through `myb::sys`, bit
//! for bit, the array forms against the single calls; failures as `Err` with the status of
//! the library (and the index of the element of an array form).
use std::ffi::{c_void, CString};
use std::fmt::Debug;
use std::path::{Path, PathBuf};
use std::ptr;
use std::rc::Rc;

use myb::sys;
use myb::{
    Calendar, Cip, CipOrientation, Eop, Ephemeris, Equinox, EquinoxFrame, FrameGraph, GraphOrientation, KernelPool, Pck,
    PckOrientation, Rotation, Status, Time, TimeSystem, Wgccre, WgccreOrientation,
};

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

fn c_path(path: &Path) -> CString {
    CString::new(path.to_str().unwrap()).unwrap()
}

fn pool(kernels: &[PathBuf]) -> myb::Result<KernelPool> {
    let mut pool = KernelPool::new()?;
    for kernel in kernels {
        pool.load_file(kernel)?;
    }
    Ok(pool)
}

fn time_kernels() -> Vec<PathBuf> {
    ["time/time-scales.tk", "time/tdb-tt-fb1990.tk", "time/utc-tai.tk"].iter().map(|k| data(k)).collect()
}

fn utc(year: i64, month: i32, day: i32, hour: i32) -> Calendar {
    Calendar { year, month, day, hour, minute: 17, second: 12.625 }
}

#[test]
fn rotations_and_the_models_of_the_earth_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    // rotations
    let a = myb::rotation_euler313(&[0.3, 0.4, 0.5], &[1e-5, 2e-5, 7e-5]);
    let b = myb::rotation_euler313(&[-0.2, 1.1, 2.0], &[0.0, 0.0, 1e-4]);
    let (state, position) = ([7.0e6, 1.0e5, -2.0e5, 10.0, 7.5e3, -1.0e3], [1.0, 2.0, 3.0]);
    // SAFETY: arrays of the sizes of the C interface
    unsafe {
        let mut r = Rotation::default();
        sys::myb_rotation_euler313([0.3, 0.4, 0.5].as_ptr(), [1e-5, 2e-5, 7e-5].as_ptr(), &mut r);
        same(&a, &r);
        sys::myb_rotation_compose(&b, &a, &mut r);
        same(&myb::rotation_compose(&b, &a), &r);
        sys::myb_rotation_invert(&a, &mut r);
        same(&myb::rotation_invert(&a), &r);
        let (mut p, mut s) = ([0.0; 3], [0.0; 6]);
        sys::myb_rotation_apply_position(&a, position.as_ptr(), p.as_mut_ptr());
        sys::myb_rotation_apply_state(&a, state.as_ptr(), s.as_mut_ptr());
        assert_eq!((bits(&myb::rotation_apply_position(&a, &position)), bits(&myb::rotation_apply_state(&a, &state))), (bits(&p), bits(&s)));
        sys::myb_cip_c2i([1e-4, -2e-5, 3e-8].as_ptr(), ptr::null(), &mut r);
        same(&myb::cip_c2i(&[1e-4, -2e-5, 3e-8], None), &r);
    }

    // the CIO-based model with the EOP
    let mut kernels = time_kernels();
    kernels.extend([data("frames/iau2006-cip.tk"), data("frames/eop-tides.tk")]);
    let kernel_pool = pool(&kernels)?;
    let system = Rc::new(TimeSystem::new(&kernel_pool)?);
    let cip = Cip::new(&kernel_pool, "IAU2006A")?;
    let mut eop = Eop::new()?;
    eop.load_c04_file(test_data("eop/eopc04-excerpt.txt"))?;
    eop.set_tides(&kernel_pool, "IERS10_EOPT")?;
    let labels = [utc(2017, 1, 3, 4), utc(2018, 6, 30, 23), utc(2019, 12, 1, 0), utc(2020, 7, 4, 12)];
    let tt = system.from_utc_many(&labels, "TT")?;
    let (values, ut1) = eop.at_many(&system, &tt, "TT")?;
    let (xys, rates) = cip.xys_many(&tt);
    let rotations = cip.terrestrial_many(&tt, &ut1, Some(&values))?;
    let (rows, first, last) = eop.rows();

    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    unsafe {
        let mut raw_pool = ptr::null_mut();
        assert_eq!(sys::myb_kpool_create(ptr::null(), &mut raw_pool), sys::MYB_OK);
        let mut diag = sys::myb_diag::default();
        for kernel in &kernels {
            assert_eq!(sys::myb_kpool_load_file(raw_pool, c_path(kernel).as_ptr(), &mut diag), sys::MYB_OK);
        }
        let (mut raw_cip, mut raw_eop) = (ptr::null_mut(), ptr::null_mut());
        assert_eq!(sys::myb_cip_create(raw_pool, c"IAU2006A".as_ptr(), ptr::null(), &mut raw_cip, &mut diag), sys::MYB_OK);
        assert_eq!(sys::myb_eop_create(ptr::null(), &mut raw_eop), sys::MYB_OK);
        assert_eq!(sys::myb_eop_load_c04_file(raw_eop, c_path(&test_data("eop/eopc04-excerpt.txt")).as_ptr(), &mut diag), sys::MYB_OK);
        assert_eq!(sys::myb_eop_set_tides(raw_eop, raw_pool, c"IERS10_EOPT".as_ptr(), &mut diag), sys::MYB_OK);
        let (mut raw_first, mut raw_last) = (0.0, 0.0);
        assert_eq!((sys::myb_eop_rows(raw_eop, &mut raw_first, &mut raw_last), raw_first, raw_last), (rows, first, last));
        let mut scale = sys::myb_scale::default();
        assert_eq!(sys::myb_timesys_scale(system.as_ptr(), c"TT".as_ptr(), &mut scale), sys::MYB_OK);
        for k in 0..tt.len() {
            let (mut v, mut u) = (sys::myb_eop_values::default(), sys::myb_time::default());
            assert_eq!(sys::myb_eop_at(raw_eop, system.as_ptr(), tt[k], scale, &mut v, &mut u), sys::MYB_OK);
            same(&(values[k], ut1[k]), &(v, u));
            let (single, single_ut1) = eop.at(&system, tt[k], "TT")?;
            same(&(single, single_ut1), &(v, u));
            let (mut x, mut r) = ([0.0; 3], [0.0; 3]);
            sys::myb_cip_xys(raw_cip, tt[k], x.as_mut_ptr(), r.as_mut_ptr());
            assert_eq!((bits(&xys[k]), bits(&rates[k])), (bits(&x), bits(&r)));
            assert_eq!((bits(&cip.xys(tt[k]).0), bits(&cip.xys(tt[k]).1)), (bits(&x), bits(&r)));
            let mut m = Rotation::default();
            sys::myb_cip_terrestrial(raw_cip, tt[k], u, &v, &mut m);
            same(&rotations[k], &m);
            same(&cip.terrestrial(tt[k], u, Some(&v)), &m);
            let mut rate = 0.0;
            let era = sys::myb_cip_era(raw_cip, u, &mut rate);
            assert_eq!((cip.era(u).0.to_bits(), cip.era(u).1.to_bits()), (era.to_bits(), rate.to_bits()));
            assert_eq!(cip.sprime(tt[k]).to_bits(), sys::myb_cip_sprime(raw_cip, tt[k]).to_bits());
        }
        // the provider of the model
        let context = sys::myb_cip_orientation { cip: raw_cip, eop: raw_eop, timesys: system.as_ptr(), tt: scale };
        let mut o = sys::myb_orientation::default();
        assert_eq!(sys::myb_orientation_cip((&context as *const sys::myb_cip_orientation as *mut sys::myb_cip_orientation).cast(), tt[1], &mut o), sys::MYB_OK);
        let provider = CipOrientation::new(&Rc::new(cip), &Rc::new(eop), &system, "TT")?;
        same(&myb::orientation_from(&provider, tt[1])?, &o);
        same(&myb::orientation_many(&provider, &tt)?[1], &o);
        // a tabulated model
        let mut tabulated = Cip::new(&kernel_pool, "IAU2006A")?;
        tabulated.tabulate(tt[0], tt[3], 3600.0)?;
        assert_eq!(sys::myb_cip_tabulate(raw_cip, tt[0], tt[3], 3600.0), sys::MYB_OK);
        let (x, _) = tabulated.xys(tt[2]);
        let mut raw_x = [0.0; 3];
        sys::myb_cip_xys(raw_cip, tt[2], raw_x.as_mut_ptr(), ptr::null_mut());
        assert_eq!(bits(&x), bits(&raw_x));
        assert_eq!(tabulated.tabulate(tt[3], tt[0], 3600.0).unwrap_err().status, Status::Domain);
        sys::myb_eop_destroy(raw_eop);
        sys::myb_cip_destroy(raw_cip);
        sys::myb_kpool_destroy(raw_pool);
    }
    let cip = Cip::new(&kernel_pool, "IAU2006A")?;
    assert_eq!(cip.terrestrial_many(&tt, &ut1[..3], None).unwrap_err().status, Status::Domain);
    let zero = cip.terrestrial_many(&tt[..1], &ut1[..1], None)?;
    same(&zero[0], &cip.terrestrial(tt[0], ut1[0], None));

    // the finals series; failures
    let mut finals = Eop::new()?;
    finals.load_finals_file(test_data("eop/finals2000A-excerpt.txt"))?;
    assert!(finals.rows().0 > 0);
    let mut outside = tt.clone();
    outside[2] = system.from_utc(&utc(1990, 1, 1, 0), "TT")?;
    let mut eop = Eop::new()?;
    eop.load_c04_file(test_data("eop/eopc04-excerpt.txt"))?;
    let error = eop.at_many(&system, &outside, "TT").unwrap_err();
    assert_eq!((error.status, error.message.as_str()), (Status::Range, "element 2"));
    assert_eq!(eop.at(&system, outside[2], "TT").unwrap_err().status, Status::Range);
    assert_eq!(eop.at(&system, tt[0], "NO SUCH").unwrap_err().status, Status::NotFound);
    let before = eop.rows();
    let error = eop.load_c04("2017 1 1 0 nonsense\n").unwrap_err();
    assert_eq!(error.status, Status::Format);
    assert!(error.message.starts_with("line 1"), "{error}");
    assert_eq!(eop.rows(), before); // the failed load left the table as it was
    assert_eq!(Cip::new(&kernel_pool, "NOPE").err().unwrap().status, Status::NotFound);

    // the equinox-based systems
    let equinox_pool = pool(&[data("frames/iau1980.tk")])?;
    let equinox = Equinox::new(&equinox_pool, "IAU1980")?;
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    unsafe {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        assert_eq!(sys::myb_equinox_create(equinox_pool.as_ptr(), c"IAU1980".as_ptr(), ptr::null(), &mut raw, &mut diag), sys::MYB_OK);
        for (frame, code) in [(EquinoxFrame::Mod, sys::MYB_EQUINOX_MOD), (EquinoxFrame::Tod, sys::MYB_EQUINOX_TOD), (EquinoxFrame::Teme, sys::MYB_EQUINOX_TEME)] {
            let many = equinox.rotation_many(frame, &tt);
            for k in 0..tt.len() {
                let mut r = Rotation::default();
                sys::myb_equinox_rotation(raw, code, tt[k], &mut r);
                same(&many[k], &r);
                same(&equinox.rotation(frame, tt[k]), &r);
            }
        }
        let teme_pef = equinox.teme_pef_many(&tt);
        for k in 0..tt.len() {
            let mut r = Rotation::default();
            sys::myb_equinox_teme_pef(raw, tt[k], &mut r);
            same(&teme_pef[k], &r);
            same(&equinox.teme_pef(tt[k]), &r);
        }
        let (mut angles, mut rates, mut nutation, mut nutation_rates, mut rate) = ([0.0; 3], [0.0; 3], [0.0; 2], [0.0; 2], 0.0);
        sys::myb_equinox_precession(raw, tt[1], angles.as_mut_ptr(), rates.as_mut_ptr());
        assert_eq!((bits(&equinox.precession(tt[1]).0), bits(&equinox.precession(tt[1]).1)), (bits(&angles), bits(&rates)));
        let eps = sys::myb_equinox_obliquity(raw, tt[1], &mut rate);
        assert_eq!((equinox.obliquity(tt[1]).0.to_bits(), equinox.obliquity(tt[1]).1.to_bits()), (eps.to_bits(), rate.to_bits()));
        sys::myb_equinox_nutation(raw, tt[1], nutation.as_mut_ptr(), nutation_rates.as_mut_ptr());
        assert_eq!((bits(&equinox.nutation(tt[1]).0), bits(&equinox.nutation(tt[1]).1)), (bits(&nutation), bits(&nutation_rates)));
        let gmst = sys::myb_equinox_gmst(raw, tt[1], &mut rate);
        assert_eq!((equinox.gmst(tt[1]).0.to_bits(), equinox.gmst(tt[1]).1.to_bits()), (gmst.to_bits(), rate.to_bits()));
        sys::myb_equinox_destroy(raw);
    }
    assert!(equinox.rotation_many(EquinoxFrame::Mod, &[]).is_empty());
    assert_eq!(Equinox::new(&equinox_pool, "NOPE").err().unwrap().status, Status::NotFound);
    Ok(())
}

// ---- WGCCRE, PCK and the graph --------------------------------------------------------------------

const TURNING: i32 = 1400099;
const TURNING_CLASS: i32 = 31099;

unsafe extern "C" fn witness(_c: *mut c_void, t: sys::myb_time, o: *mut sys::myb_orientation) -> sys::myb_status {
    // SAFETY: the output is valid
    unsafe {
        let s = sys::myb_time_seconds(t);
        let mut r = Rotation::default();
        sys::myb_rotation_euler313([0.1, 0.2, 7.29e-5 * s].as_ptr(), [0.0, 0.0, 7.29e-5].as_ptr(), &mut r);
        *o = sys::myb_orientation { rotation: r, ut: t, xp: 0.0, yp: 0.0 };
    }
    sys::MYB_OK
}

fn witness_closure(t: Time) -> myb::Result<myb::Orientation> {
    let r = myb::rotation_euler313(&[0.1, 0.2, 7.29e-5 * myb::seconds(t)], &[0.0, 0.0, 7.29e-5]);
    Ok(myb::Orientation { rotation: r, ut: t, xp: 0.0, yp: 0.0 })
}

#[test]
fn wgccre_pck_and_the_graph_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let t: Vec<Time> = (0..4).map(|k| myb::time_from_seconds(6.0e8 + 1000.0 + 1500.0 * k as f64)).collect();

    // WGCCRE
    let text_pck = pool(&[test_data("wgccre/pck00011.tpc")])?;
    let mars = Rc::new(Wgccre::new(&text_pck, 499)?);
    let many = mars.rotation_many(&t);
    let base = myb::rotation_euler313(&[0.0, 0.40909, 0.0], &[0.0; 3]);
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    unsafe {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        assert_eq!(sys::myb_wgccre_create(text_pck.as_ptr(), 499, ptr::null(), &mut raw, &mut diag), sys::MYB_OK);
        assert_eq!(mars.frame(), sys::myb_wgccre_frame(raw));
        for k in 0..t.len() {
            let mut r = Rotation::default();
            sys::myb_wgccre_rotation(raw, t[k], &mut r);
            same(&many[k], &r);
            same(&mars.rotation(t[k]), &r);
        }
        let (mut angles, mut rates) = ([0.0; 3], [0.0; 3]);
        sys::myb_wgccre_euler(raw, t[0], angles.as_mut_ptr(), rates.as_mut_ptr());
        assert_eq!((bits(&mars.euler(t[0]).0), bits(&mars.euler(t[0]).1)), (bits(&angles), bits(&rates)));
        let context = sys::myb_wgccre_orientation { model: raw, base: &base };
        let mut o = sys::myb_orientation::default();
        assert_eq!(sys::myb_orientation_wgccre((&context as *const sys::myb_wgccre_orientation as *mut sys::myb_wgccre_orientation).cast(), t[1], &mut o), sys::MYB_OK);
        let provider = WgccreOrientation::new(&mars, Some(&base));
        same(&myb::orientation_from(&provider, t[1])?, &o);
        sys::myb_wgccre_destroy(raw);
    }
    assert_eq!(Wgccre::new(&text_pck, 123456).err().unwrap().status, Status::NotFound);

    // binary PCK, from the file and from memory
    let mut loaded = Pck::new()?;
    loaded.load_file(test_data("frames/graph.bpc"))?;
    let mut in_memory = Pck::new()?;
    in_memory.load_memory(std::fs::read(test_data("frames/graph.bpc")).unwrap())?;
    let pck = Rc::new(loaded);
    let (rotations, bases) = pck.rotation_many(TURNING_CLASS, &t)?;
    let (raw_pck, raw_spk) = unsafe {
        // SAFETY: the calls of the C interface with valid arguments
        let (mut raw_pck, mut raw_spk) = (ptr::null_mut(), ptr::null_mut());
        let mut diag = sys::myb_diag::default();
        assert_eq!(sys::myb_pck_create(ptr::null(), &mut raw_pck), sys::MYB_OK);
        assert_eq!(sys::myb_pck_load_file(raw_pck, c_path(&test_data("frames/graph.bpc")).as_ptr(), &mut diag), sys::MYB_OK);
        assert_eq!(sys::myb_spk_create(ptr::null(), &mut raw_spk), sys::MYB_OK);
        assert_eq!(sys::myb_spk_load_file(raw_spk, c_path(&test_data("frames/graph.bsp")).as_ptr(), &mut diag), sys::MYB_OK);
        (raw_pck, raw_spk)
    };
    // SAFETY: the set is alive
    unsafe {
        for k in 0..t.len() {
            let (mut r, mut b) = (Rotation::default(), 0);
            assert_eq!(sys::myb_pck_rotation(raw_pck, TURNING_CLASS, t[k], &mut r, &mut b), sys::MYB_OK);
            same(&(rotations[k], bases[k]), &(r, b));
            same(&pck.rotation(TURNING_CLASS, t[k])?, &(r, b));
            same(&in_memory.rotation(TURNING_CLASS, t[k])?, &(r, b));
        }
        let (mut angles, mut rates, mut b) = ([0.0; 3], [0.0; 3], 0);
        assert_eq!(sys::myb_pck_euler(raw_pck, TURNING_CLASS, t[0], angles.as_mut_ptr(), rates.as_mut_ptr(), &mut b), sys::MYB_OK);
        let (e, r, eb) = pck.euler(TURNING_CLASS, t[0])?;
        assert_eq!((bits(&e), bits(&r), eb), (bits(&angles), bits(&rates), b));
        let context = sys::myb_pck_orientation { pck: raw_pck, id: TURNING_CLASS, base: ptr::null() };
        let mut o = sys::myb_orientation::default();
        assert_eq!(sys::myb_orientation_pck((&context as *const sys::myb_pck_orientation as *mut sys::myb_pck_orientation).cast(), t[2], &mut o), sys::MYB_OK);
        same(&myb::orientation_from(&PckOrientation::new(&pck, TURNING_CLASS, None), t[2])?, &o);
    }
    let mut late = t.clone();
    late[1] = myb::time_from_seconds(1.0e10);
    assert_eq!(pck.rotation_many(TURNING_CLASS, &late).unwrap_err().message, "element 1");
    assert_eq!(pck.rotation(TURNING_CLASS, late[1]).unwrap_err().status, Status::Range);
    assert_eq!(Pck::new()?.load_memory(vec![0u8; 64]).unwrap_err().status, Status::Format);

    // the graph of the inertial frames of the kernel and the turning frame of the PCK
    let inertial = pool(&[data("frames/naif-inertial.tk")])?;
    let mut graph = FrameGraph::new()?;
    for id in [2, 3, 13, 17, 18] {
        graph.add_tkframe(&inertial, id)?;
    }
    assert_eq!(graph.add_tkframe(&inertial, 17).unwrap_err().status, Status::Domain);
    assert_eq!(graph.add_tkframe(&inertial, 99).unwrap_err().status, Status::NotFound);
    graph.add(TURNING, 1, PckOrientation::new(&pck, TURNING_CLASS, None), None)?;
    graph.add_fixed(5000, 17, &base)?;
    assert_eq!(graph.add_fixed(5000, 17, &base).unwrap_err().status, Status::Domain);
    assert_eq!(graph.add(0, 1, PckOrientation::new(&pck, TURNING_CLASS, None), None).unwrap_err().status, Status::Domain);
    assert_eq!(graph.add(6000, 1, PckOrientation::new(&pck, TURNING_CLASS, None), Some("TT")).unwrap_err().status, Status::Domain);
    let graph = Rc::new(graph);
    let turned = graph.rotation_many(17, TURNING, &t)?;
    let mut ephemeris = Ephemeris::new()?;
    ephemeris.load_file(test_data("frames/graph.bsp"))?;
    let ephemeris = Rc::new(ephemeris);
    let states = ephemeris.state_in_many(&graph, 1602, 399, &t, 13)?;
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    unsafe {
        let mut raw_graph = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        assert_eq!(sys::myb_frame_graph_create(ptr::null(), ptr::null(), &mut raw_graph), sys::MYB_OK);
        for id in [2, 3, 13, 17, 18] {
            assert_eq!(sys::myb_frame_graph_add_tkframe(raw_graph, inertial.as_ptr(), id, &mut diag), sys::MYB_OK);
        }
        let turning = sys::myb_pck_orientation { pck: raw_pck, id: TURNING_CLASS, base: ptr::null() };
        let context = (&turning as *const sys::myb_pck_orientation as *mut sys::myb_pck_orientation).cast();
        assert_eq!(sys::myb_frame_graph_add(raw_graph, TURNING, 1, Some(sys::myb_orientation_pck), context, ptr::null()), sys::MYB_OK);
        assert_eq!(sys::myb_frame_graph_add_fixed(raw_graph, 5000, 17, &base), sys::MYB_OK);
        for k in 0..t.len() {
            let mut r = Rotation::default();
            assert_eq!(sys::myb_frame_graph_rotation(raw_graph, 17, TURNING, t[k], &mut r), sys::MYB_OK);
            same(&turned[k], &r);
            same(&graph.rotation(17, TURNING, t[k])?, &r);
            let mut y = [0.0; 6];
            assert_eq!(sys::myb_spk_state_in(raw_spk, raw_graph, 1602, 399, t[k], 13, y.as_mut_ptr()), sys::MYB_OK);
            assert_eq!(bits(&states[k]), bits(&y));
            assert_eq!(bits(&ephemeris.state_in(&graph, 1602, 399, t[k], 13)?), bits(&y));
        }
        let mut r = Rotation::default();
        assert_eq!(sys::myb_frame_graph_rotation(raw_graph, 5000, 2, t[0], &mut r), sys::MYB_OK);
        same(&graph.rotation(5000, 2, t[0])?, &r);
        let pair = sys::myb_graph_orientation { graph: raw_graph, from: 17, to: TURNING };
        let mut o = sys::myb_orientation::default();
        assert_eq!(sys::myb_orientation_graph((&pair as *const sys::myb_graph_orientation as *mut sys::myb_graph_orientation).cast(), t[3], &mut o), sys::MYB_OK);
        same(&myb::orientation_from(&GraphOrientation::new(&graph, 17, TURNING), t[3])?, &o);
        // the fixed frame of the kernel itself
        let (mut fixed, mut relative) = (Rotation::default(), [0 as std::ffi::c_char; 33]);
        assert_eq!(sys::myb_tkframe(inertial.as_ptr(), 17, &mut fixed, relative.as_mut_ptr(), &mut diag), sys::MYB_OK);
        let (rotation, name) = myb::tkframe(&inertial, 17)?;
        same(&rotation, &fixed);
        assert_eq!(name, "J2000");
        sys::myb_frame_graph_destroy(raw_graph);
        sys::myb_spk_destroy(raw_spk);
        sys::myb_pck_destroy(raw_pck);
    }
    // a body of the set in a frame of the graph
    let body = myb::SpkBody::in_frame(&ephemeris, &graph, 1602, 399, 13);
    assert_eq!(bits(&body.state(t[2])?), bits(&states[2]));
    // failures
    assert_eq!(graph.rotation(17, 777, t[0]).unwrap_err().status, Status::NotFound);
    assert_eq!(ephemeris.state_in(&graph, 1602, 399, t[0], 777).unwrap_err().status, Status::NotFound);
    let error = ephemeris.state_in_many(&graph, 1602, 399, &late, 13).unwrap_err();
    assert_eq!((error.status, error.message.as_str()), (Status::Range, "element 1"));
    assert_eq!(graph.rotation_many(17, TURNING, &late).unwrap_err().message, "element 1");
    assert_eq!(myb::tkframe(&inertial, 99).unwrap_err().status, Status::NotFound);

    // a graph on TDB with a provider on TT, a closure against a C function
    let system = Rc::new(TimeSystem::new(&pool(&time_kernels())?)?);
    let mut timed = FrameGraph::with_time(&system, "TDB")?;
    timed.add(10, 1, witness_closure, Some("TT"))?;
    timed.add(11, 10, |_t: Time| -> myb::Result<myb::Orientation> { Err(myb::Error::new(Status::Range)) }, None)?;
    let timed = Rc::new(timed);
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    unsafe {
        let options = sys::myb_frame_graph_options { timesys: system.as_ptr(), scale: c"TDB".as_ptr() };
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_frame_graph_create(&options, ptr::null(), &mut raw), sys::MYB_OK);
        assert_eq!(sys::myb_frame_graph_add(raw, 10, 1, Some(witness), ptr::null_mut(), c"TT".as_ptr()), sys::MYB_OK);
        let mut r = Rotation::default();
        assert_eq!(sys::myb_frame_graph_rotation(raw, 1, 10, t[1], &mut r), sys::MYB_OK);
        same(&timed.rotation(1, 10, t[1])?, &r);
        sys::myb_frame_graph_destroy(raw);
    }
    assert_eq!(timed.rotation(1, 11, t[1]).unwrap_err().status, Status::Range);
    assert_eq!(FrameGraph::with_time(&system, "XYZ").err().unwrap().status, Status::NotFound);
    Ok(())
}

// ---- SPK type 10 with an equinox model ------------------------------------------------------------

#[cfg(feature = "sgp4")]
#[test]
fn element_sets_in_spk_files_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    use myb::sgp4;
    use sys::sgp4 as ext;
    let mut kernels = time_kernels();
    kernels.push(data("frames/iau1980.tk"));
    let kernel_pool = pool(&kernels)?;
    let system = TimeSystem::new(&kernel_pool)?;
    let equinox = Equinox::new(&kernel_pool, "IAU1980")?;
    let mut ephemeris = Ephemeris::new()?;
    ephemeris.load_file(test_data("spk/t10.bsp"))?;
    let segment = ephemeris.segments()[0];
    let et = myb::time_from_seconds(0.5 * (segment.start + segment.stop));
    let state = sgp4::sgp4_spk10_state(&ephemeris, 0, &system, &equinox, et)?;
    let tle = sgp4::tle_parse(
        "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753",
        "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667",
        true,
    )?;
    let packet = sgp4::sgp4_spk10_packet(&tle, &system, Some(&equinox))?;
    let (mut raw_state, mut raw_packet) = ([0.0; 6], [0.0; 14]);
    // SAFETY: the set, the system and the model are alive; arrays of the sizes asked for
    unsafe {
        assert_eq!(ext::myb_sgp4_spk10_state(ephemeris.as_ptr(), 0, system.as_ptr(), equinox.as_ptr(), et, raw_state.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(ext::myb_sgp4_spk10_packet(&tle, system.as_ptr(), equinox.as_ptr(), raw_packet.as_mut_ptr()), sys::MYB_OK);
    }
    assert_eq!((bits(&state), bits(&packet)), (bits(&raw_state), bits(&raw_packet)));
    assert!(packet[10] != 0.0); // the nutation, which no equinox leaves zero
    assert_eq!(sgp4::sgp4_spk10_packet(&tle, &system, None)?[10], 0.0);
    // failures: no such segment; a segment of another type
    let segments = ephemeris.segments().len();
    assert_eq!(sgp4::sgp4_spk10_state(&ephemeris, segments, &system, &equinox, et).unwrap_err().status, Status::Domain);
    let mut other = Ephemeris::new()?;
    other.load_file(test_data("spk/t02_03.bsp"))?;
    assert_eq!(sgp4::sgp4_spk10_state(&other, 0, &system, &equinox, et).unwrap_err().status, Status::Format);
    Ok(())
}
