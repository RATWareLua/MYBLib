// SPDX-License-Identifier: BSD-3-Clause
//! Files (`myb/io.h`, `myb/kin.h`, `myb/ccsds.h`, `myb/ck.h`, `ext/sgp4`) through the safe
//! layer against the same C calls, bit for bit; failures as `Err` with the status of the
//! library.
use std::ffi::{c_char, c_void, CString};
use std::fmt::Debug;
use std::path::PathBuf;
use std::ptr;
use std::rc::Rc;

use myb::sys;
use myb::{
    Ck, CkPackets, CkSegmentHeader, CkWriter, Ephemeris, KernelPool, KernelVariable, Notation, OcmBlock, OcmData, OcmKind, OcmTime,
    Sclk, Sp3, Sp3Data, SpkSegmentHeader, SpkWriter, State, Status, Table,
};

const GM: f64 = 3.986004418e14;

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn same<T: Debug>(a: &T, b: &T) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

fn test_data(relative: &str) -> PathBuf {
    std::env::var_os("MYB_TEST_DATA_DIR")
        .map_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../tests/data"), PathBuf::from)
        .join(relative)
}

/// A text a writer of the C interface gives into a buffer large enough.
fn raw_text(write: impl Fn(*mut c_char, usize, *mut usize) -> sys::myb_status) -> String {
    let mut buffer = vec![0u8; 1 << 16];
    let mut length = 0usize;
    assert_eq!(write(buffer.as_mut_ptr().cast(), buffer.len(), &mut length), sys::MYB_OK);
    buffer.truncate(length);
    String::from_utf8(buffer).unwrap()
}

#[test]
fn numbers_tables_and_kernel_variables() -> myb::Result<()> {
    let (value, consumed) = myb::parse_double("1.25e3xyz")?;
    let (mut raw_value, mut raw_consumed) = (0.0, 0usize);
    // SAFETY: the text is nine bytes
    assert_eq!(unsafe { sys::myb_parse_double(c"1.25e3xyz".as_ptr(), 9, &mut raw_value, &mut raw_consumed) }, sys::MYB_OK);
    assert_eq!((value.to_bits(), consumed), (raw_value.to_bits(), raw_consumed));
    assert_eq!(myb::parse_double("abc").unwrap_err().status, Status::Format);
    for x in [0.1, 1.0 / 3.0, -2.5e-300, 1e16, f64::NAN, -0.0] {
        let mut out = [0u8; 32];
        // SAFETY: out holds 32 characters
        let length = unsafe { sys::myb_format_double(x, out.as_mut_ptr().cast()) };
        assert_eq!(myb::format_double(x).as_bytes(), &out[..length]);
        // SAFETY: the buffers are as large as said
        let fixed = raw_text(|b, s, l| unsafe { sys::myb_format_fixed(x, 20, b, s, l) });
        assert_eq!(myb::format_fixed(x, 20)?, fixed);
        let scientific = raw_text(|b, s, l| unsafe { sys::myb_format_scientific(x, 17, b, s, l) });
        assert_eq!(myb::format_scientific(x, 17)?, scientific);
    }
    assert_eq!(myb::format_fixed(1.0, 61).unwrap_err().status, Status::Domain);
    assert_eq!(myb::format_fixed(1e300, 60)?.len(), 362); // longer than the first buffer

    // tables
    let table = Table::parse("a,b\n1,2.5\n3,\n")?;
    assert_eq!((table.rows(), table.columns(), table.name(0).as_deref(), table.name(2)), (2, 2, Some("a"), None));
    assert_eq!((table.row(0).unwrap(), table.row(2)), (&[1.0, 2.5][..], None));
    assert!(table.row(1).unwrap()[1].is_nan());
    let values = [1.0, 2.5, f64::NAN, 4.0];
    let text = myb::write_table(Some(&["a", "b"]), 2, &values)?;
    let names = [c"a".as_ptr(), c"b".as_ptr()];
    // SAFETY: the names and values live for the call
    let raw = raw_text(|b, s, l| unsafe { sys::myb_table_write(names.as_ptr(), 2, values.as_ptr(), 2, b, s, l) });
    assert_eq!(text, raw);
    assert_eq!(Table::parse(&text)?.row(1).unwrap()[1], 4.0);
    let error = Table::parse("1,2\n3\n").err().unwrap();
    assert_eq!(error.status, Status::Format);
    assert!(error.message.starts_with("line 2"), "{error}");
    assert_eq!(myb::write_table(None, 3, &values).unwrap_err().status, Status::Domain);
    assert_eq!(Table::load_file(test_data("io/missing.csv")).err().unwrap().status, Status::Io);

    // kernel variables
    let mut pool = KernelPool::new()?;
    pool.load_text("\\begindata\nX = ( 1.5 2.5D3 )\nNAME = ( 'abc' 'it''s' )\n\\begintext\n")?;
    assert_eq!(pool.describe("X")?, KernelVariable::Numbers(2));
    assert_eq!(pool.describe("NAME")?, KernelVariable::Strings(2));
    assert_eq!(pool.numbers("X")?, vec![1.5, 2500.0]);
    assert_eq!((pool.string("NAME", 0)?, pool.string("NAME", 1)?), ("abc".to_string(), "it's".to_string()));
    assert_eq!(pool.describe("NONE").unwrap_err().status, Status::NotFound);
    assert_eq!(pool.numbers("NAME").unwrap_err().status, Status::Format);
    assert_eq!(pool.string("NAME", 2).unwrap_err().status, Status::Range);
    assert_eq!(pool.load_text("\\begindata\nX = ( 1\n").unwrap_err().status, Status::Format);
    Ok(())
}

// ---- SPK and SP3 ----------------------------------------------------------------------------------

fn kepler(seconds: f64) -> State {
    let y0 = myb::state_from_elements(GM, &myb::Elements { p: 7.0e6, e: 0.01, i: 0.9, raan: 0.3, argp: 0.2, nu: 0.1 }).unwrap();
    myb::kepler_propagate(GM, &y0, seconds).unwrap()
}

unsafe extern "C" fn raw_kepler(_c: *mut c_void, t: sys::myb_time, y: *mut f64) -> sys::myb_status {
    // SAFETY: y holds six numbers
    unsafe { *(y as *mut State) = kepler(sys::myb_time_seconds(t)) };
    sys::MYB_OK
}

#[test]
fn the_spk_writer_and_sp3_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let epochs: Vec<f64> = (0..=10).map(|k| 60.0 * k as f64).collect();
    let states: Vec<State> = epochs.iter().map(|&t| kepler(t)).collect();
    let mut writer = SpkWriter::new(Some("myb test"))?;
    let hermite = SpkSegmentHeader { target: -1001, center: 399, frame: 1, name: "HERMITE".into() };
    writer.add_states(&hermite, 13, 7, &epochs, &states)?;
    let chebyshev = SpkSegmentHeader { target: -1002, center: 399, frame: 1, name: "CHEBYSHEV".into() };
    let provider = |t: myb::Time| -> myb::Result<State> { Ok(kepler(myb::seconds(t))) };
    let error = writer.add_chebyshev(&chebyshev, 3, &provider, 0.0, 600.0, 4, 12, true)?;
    assert!(error.is_some_and(|e| e < 1e-3), "{error:?}");
    writer.add_raw(&SpkSegmentHeader { target: -1003, center: 399, frame: 1, name: "RAW".into() }, 2, 0.0, 600.0, &[300.0, 300.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 600.0, 11.0, 1.0])?;
    let image = writer.write()?;
    assert_eq!(image.len(), writer.size());

    let raw_image = unsafe {
        // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_spk_writer_create(c"myb test".as_ptr(), ptr::null(), &mut raw), sys::MYB_OK);
        let header = sys::myb_spk_segment_header { target: -1001, center: 399, frame: 1, name: c"HERMITE".as_ptr() };
        assert_eq!(sys::myb_spk_writer_add_states(raw, &header, 13, 7, epochs.as_ptr(), states.as_ptr().cast(), 11), sys::MYB_OK);
        let header = sys::myb_spk_segment_header { target: -1002, center: 399, frame: 1, name: c"CHEBYSHEV".as_ptr() };
        let mut max_error = 0.0;
        assert_eq!(sys::myb_spk_writer_add_chebyshev(raw, &header, 3, Some(raw_kepler), ptr::null_mut(), 0.0, 600.0, 4, 12, &mut max_error), sys::MYB_OK);
        assert_eq!(max_error.to_bits(), error.unwrap().to_bits());
        let header = sys::myb_spk_segment_header { target: -1003, center: 399, frame: 1, name: c"RAW".as_ptr() };
        let words = [300.0, 300.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 600.0, 11.0, 1.0];
        assert_eq!(sys::myb_spk_writer_add_raw(raw, &header, 2, 0.0, 600.0, words.as_ptr(), words.len()), sys::MYB_OK);
        let mut image = vec![0u8; sys::myb_spk_writer_size(raw)];
        assert_eq!(sys::myb_spk_writer_write(raw, image.as_mut_ptr().cast(), image.len()), sys::MYB_OK);
        sys::myb_spk_writer_destroy(raw);
        image
    };
    assert!(image == raw_image, "the files differ");
    let mut ephemeris = Ephemeris::new()?;
    ephemeris.load_memory(image)?;
    let t = myb::time_from_seconds(123.4);
    let mut raw_state = [0.0; 6];
    // SAFETY: the set and its image are alive
    unsafe {
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_spk_create(ptr::null(), &mut raw), sys::MYB_OK);
        let mut diag = sys::myb_diag::default();
        assert_eq!(sys::myb_spk_load_memory(raw, raw_image.as_ptr().cast(), raw_image.len(), &mut diag), sys::MYB_OK);
        assert_eq!(sys::myb_spk_state(raw, -1002, 399, t, raw_state.as_mut_ptr(), ptr::null_mut()), sys::MYB_OK);
        sys::myb_spk_destroy(raw);
    }
    assert_eq!(bits(&ephemeris.state(-1002, 399, t)?), bits(&raw_state));
    assert_eq!(bits(&ephemeris.state(-1001, 399, myb::time_from_seconds(120.0))?), bits(&states[2]));
    let body = myb::SpkBody::new(&Rc::new(ephemeris), -1001, 399, 0);
    assert_eq!(bits(&body.state(myb::time_from_seconds(120.0))?), bits(&states[2]));

    // failures
    let mut empty = SpkWriter::new(None)?;
    assert_eq!(empty.write().unwrap_err().status, Status::Domain);
    assert_eq!(empty.add_states(&hermite, 13, 7, &epochs[..3], &states).unwrap_err().status, Status::Domain);
    let long = SpkSegmentHeader { name: "N".repeat(41), ..hermite.clone() };
    assert_eq!(empty.add_states(&long, 13, 7, &epochs, &states).unwrap_err().status, Status::Domain);
    let lost = |_t: myb::Time| -> myb::Result<State> { Err(myb::Error::new(Status::Range)) };
    assert_eq!(empty.add_chebyshev(&chebyshev, 3, &lost, 0.0, 600.0, 4, 12, false).unwrap_err().status, Status::Range);
    assert_eq!(SpkWriter::new(Some("\u{1}bad")).err().unwrap().status, Status::Domain);
    assert_eq!(Ephemeris::new()?.load_memory(vec![0u8; 100]).unwrap_err().status, Status::Format);

    // SP3: read, and written again as the C writer writes it
    let sp3 = Sp3::load_file(test_data("sp3/ilrsa.lageos1.250607.excerpt.sp3"))?;
    let (satellites, epochs) = (sp3.satellite_count(), sp3.epoch_count());
    assert_eq!((satellites, sp3.satellite(0).as_deref(), sp3.time_system(), sp3.frame()), (1, Some("L51"), "UTC".to_string(), "SLR20".to_string()));
    let state = sp3.state(0, 0).unwrap();
    assert_eq!(state.given, 3);
    assert!(sp3.state(0, epochs).is_none() && sp3.satellite(1).is_none());
    let raw = sp3_raw(&test_data("sp3/ilrsa.lageos1.250607.excerpt.sp3"));
    assert_eq!((bits(&state.position), bits(&state.velocity)), (bits(&raw.0), bits(&raw.1)));
    same(&sp3.epoch(0).unwrap(), &raw.2);
    assert_eq!(sp3.flags(0, 0), Some(0));
    assert!(sp3.clock(0, 0).is_some() && sp3.sigma(0, 0).is_some() && sp3.accuracy(0).is_some());
    let n = epochs.min(4);
    let mut data = Sp3Data {
        version: b'd',
        file_type: b'L',
        time_system: "UTC".into(),
        frame: "SLR20".into(),
        orbit_type: "FIT".into(),
        agency: "MYB".into(),
        data_used: "SLR".into(),
        satellites: vec!["L51".into()],
        epochs: (0..n).map(|k| sp3.epoch(k).unwrap()).collect(),
        interval: 360.0,
        ..Default::default()
    };
    for k in 0..n {
        let s = sp3.state(0, k).unwrap();
        data.positions.extend(s.position);
        data.velocities.get_or_insert_with(Vec::new).extend(s.velocity);
    }
    let text = myb::write_sp3(&data)?;
    let texts: Vec<CString> = ["UTC", "SLR20", "FIT", "MYB", "SLR", "L51"].iter().map(|s| CString::new(*s).unwrap()).collect();
    let satellite = [texts[5].as_ptr()];
    let made = sys::myb_sp3_data {
        version: b'd' as c_char,
        file_type: b'L' as c_char,
        time_system: texts[0].as_ptr(),
        frame: texts[1].as_ptr(),
        orbit_type: texts[2].as_ptr(),
        agency: texts[3].as_ptr(),
        data_used: texts[4].as_ptr(),
        satellite_count: 1,
        satellites: satellite.as_ptr(),
        epoch_count: n,
        epochs: data.epochs.as_ptr(),
        interval: 360.0,
        positions: data.positions.as_ptr(),
        velocities: data.velocities.as_ref().unwrap().as_ptr(),
        clocks: ptr::null(),
        clock_rates: ptr::null(),
        accuracy: ptr::null(),
        sigmas: ptr::null(),
        flags: ptr::null(),
    };
    // SAFETY: every pointer of `made` is alive here
    assert_eq!(text, raw_text(|b, s, l| unsafe { sys::myb_sp3_write(&made, b, s, l) }));
    let back = Sp3::parse(&text)?;
    assert_eq!((back.epoch_count(), bits(&back.state(0, 1).unwrap().position)), (n, bits(&sp3.state(0, 1).unwrap().position)));
    data.positions.pop();
    assert_eq!(myb::write_sp3(&data).unwrap_err().status, Status::Domain);
    assert_eq!(Sp3::parse("not an sp3 file\n").err().unwrap().status, Status::Format);
    Ok(())
}

fn sp3_raw(path: &std::path::Path) -> ([f64; 3], [f64; 3], myb::Calendar) {
    let path = CString::new(path.to_str().unwrap()).unwrap();
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    unsafe {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        assert_eq!(sys::myb_sp3_load_file(path.as_ptr(), ptr::null(), &mut raw, &mut diag), sys::MYB_OK);
        let (mut p, mut v, mut epoch) = ([0.0; 3], [0.0; 3], myb::Calendar::default());
        assert_eq!(sys::myb_sp3_state(raw, 0, 0, p.as_mut_ptr(), v.as_mut_ptr()), 3);
        sys::myb_sp3_epoch(raw, 0, &mut epoch);
        sys::myb_sp3_destroy(raw);
        (p, v, epoch)
    }
}

// ---- Orbit data messages --------------------------------------------------------------------------

fn calendar(minute: i32, second: f64) -> myb::Calendar {
    myb::Calendar { year: 2026, month: 3, day: 1, hour: 12, minute, second }
}

#[test]
fn messages_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    // OPM
    let mut opm = myb::Opm::default();
    myb::set_odm_text(&mut opm.originator, "MYB")?;
    myb::set_odm_text(&mut opm.object_name, "SAT")?;
    myb::set_odm_text(&mut opm.object_id, "2026-001A")?;
    myb::set_odm_text(&mut opm.center_name, "EARTH")?;
    myb::set_odm_text(&mut opm.ref_frame, "EME2000")?;
    myb::set_odm_text(&mut opm.time_system, "UTC")?;
    opm.creation_date = calendar(0, 0.0);
    opm.epoch = calendar(5, 12.5);
    opm.state = kepler(0.0);
    for notation in [Notation::Kvn, Notation::Xml] {
        let text = myb::write_opm(&opm, notation)?;
        let raw = raw_text(|b, s, l| unsafe {
            // SAFETY: the buffers are as large as said
            sys::myb_opm_write_as(&opm, if notation == Notation::Xml { sys::MYB_ODM_XML } else { sys::MYB_ODM_KVN }, b, s, l)
        });
        assert_eq!(text, raw);
        let back = myb::parse_opm(&text)?;
        let mut raw_back = myb::Opm::default();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the text is text.len() bytes
        assert_eq!(unsafe { sys::myb_opm_parse(text.as_ptr().cast(), text.len(), &mut raw_back, &mut diag) }, sys::MYB_OK);
        same(&back, &raw_back);
        assert_eq!(myb::odm_text(&back.object_name), "SAT");
    }
    assert_eq!(myb::set_odm_text(&mut opm.object_name, &"X".repeat(64)).unwrap_err().status, Status::Domain);
    let error = myb::parse_opm("CCSDS_OPM_VERS = 2.0\nNOT A LINE\n").err().unwrap();
    assert_eq!(error.status, Status::Format);
    assert!(error.message.starts_with("line"), "{error}");

    // OMM
    let mut omm = myb::Omm::default();
    for (field, text) in [
        (&mut omm.originator, "MYB"),
        (&mut omm.object_name, "SAT"),
        (&mut omm.object_id, "2026-001A"),
        (&mut omm.center_name, "EARTH"),
        (&mut omm.ref_frame, "TEME"),
        (&mut omm.time_system, "UTC"),
        (&mut omm.mean_element_theory, "SGP4"),
    ] {
        myb::set_odm_text(field, text)?;
    }
    omm.creation_date = calendar(0, 0.0);
    omm.epoch = calendar(1, 0.25);
    omm.has_mean_motion = 1;
    omm.mean_motion = 15.5;
    omm.eccentricity = 0.001;
    omm.inclination = 51.6;
    omm.ra_of_asc_node = 10.0;
    omm.arg_of_pericenter = 20.0;
    omm.mean_anomaly = 30.0;
    for notation in [Notation::Kvn, Notation::Xml] {
        let text = myb::write_omm(&omm, notation)?;
        let raw = raw_text(|b, s, l| unsafe {
            // SAFETY: the buffers are as large as said
            if notation == Notation::Xml { sys::myb_omm_write_as(&omm, sys::MYB_ODM_XML, b, s, l) } else { sys::myb_omm_write(&omm, b, s, l) }
        });
        assert_eq!(text, raw);
        let back = myb::parse_omm(&text)?;
        assert_eq!((back.mean_motion, myb::odm_text(&back.mean_element_theory)), (15.5, "SGP4".to_string()));
    }

    // OEM with a covariance
    let mut covariance = [0.0; 36];
    for k in 0..6 {
        covariance[7 * k] = 1.0 + k as f64;
    }
    covariance[1] = 0.25;
    covariance[6] = 0.25;
    let segment = myb::OemSegment {
        object_name: "SAT".into(),
        object_id: "2026-001A".into(),
        center_name: "EARTH".into(),
        ref_frame: "EME2000".into(),
        time_system: "UTC".into(),
        interpolation: None,
        epochs: vec![calendar(0, 0.0), calendar(1, 0.0)],
        states: vec![kepler(0.0), kepler(60.0)],
        covariance_frame: Some("RTN".into()),
        covariances: vec![(calendar(0, 0.0), covariance)],
    };
    let text = myb::write_oem("MYB", &calendar(0, 0.0), std::slice::from_ref(&segment), Notation::Kvn)?;
    let message = myb::Oem::parse(&text)?;
    assert_eq!((message.covariance_count(0), message.covariance_count(1)), (1, 0));
    let read = message.covariance(0, 0).unwrap();
    assert_eq!(read.ref_frame.as_deref(), Some("RTN"));
    assert_eq!(message.header("ORIGINATOR").as_deref(), Some("MYB"));
    let (_, y, acceleration) = message.state_with_acceleration(0, 1).unwrap();
    assert!(acceleration.is_none() && (y[0] - segment.states[1][0]).abs() < 1e-6);
    let mut raw_matrix = [0.0; 36];
    let mut raw_frame = ptr::null();
    let mut raw_epoch = myb::Calendar::default();
    // SAFETY: the covariance is one of the segment of a message read in C
    unsafe {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        assert_eq!(sys::myb_oem_load(text.as_ptr().cast(), text.len(), ptr::null(), &mut raw, &mut diag), sys::MYB_OK);
        sys::myb_oem_covariance(raw, 0, 0, &mut raw_epoch, &mut raw_frame, raw_matrix.as_mut_ptr());
        sys::myb_oem_destroy(raw);
    }
    assert_eq!(bits(&read.matrix), bits(&raw_matrix));
    same(&read.epoch, &raw_epoch);
    assert!(message.covariance(0, 1).is_none());

    // OCM: written from keywords and lines, and read again
    let line = |offset: f64, y: &State| myb::ocm_format_line(None, offset, &y.map(|v| v / 1000.0));
    let data = OcmData {
        originator: "MYB".into(),
        creation_date: calendar(0, 0.0),
        meta: vec![("TIME_SYSTEM".into(), "UTC".into()), ("EPOCH_TZERO".into(), "2026-03-01T12:00:00".into())],
        blocks: vec![
            OcmBlock {
                kind: OcmKind::Trajectory,
                pairs: vec![("CENTER_NAME".into(), "EARTH".into()), ("TRAJ_REF_FRAME".into(), "ICRF3".into()), ("TRAJ_TYPE".into(), "CARTPV".into())],
                lines: vec![line(0.0, &kepler(0.0))?, line(60.0, &kepler(60.0))?],
            },
            OcmBlock { kind: OcmKind::Physical, pairs: vec![("WET_MASS".into(), "100.0".into())], lines: Vec::new() },
        ],
    };
    for notation in [Notation::Kvn, Notation::Xml] {
        let text = myb::write_ocm(&data, notation)?;
        let ocm = myb::Ocm::parse(&text)?;
        assert_eq!((ocm.block_count(OcmKind::Trajectory), ocm.block_count(OcmKind::Physical), ocm.line_count(OcmKind::Trajectory, 0)), (1, 1, 2));
        assert_eq!(ocm.meta("TIME_SYSTEM").as_deref(), Some("UTC"));
        assert_eq!(ocm.header("ORIGINATOR").as_deref(), Some("MYB"));
        assert_eq!(ocm.value(OcmKind::Physical, 0, "WET_MASS").as_deref(), Some("100.0"));
        let (time, state) = ocm.state(0, 1)?;
        assert!(matches!(time, OcmTime::Offset(t) if t == 60.0));
        let (mut absolute, mut epoch, mut offset, mut raw_state, mut count) = (0, myb::Calendar::default(), 0.0, [0.0; 9], 0);
        let (_, values) = ocm.line(OcmKind::Trajectory, 0, 1)?;
        // SAFETY: a message read in C; the line is one of its block
        unsafe {
            let mut raw = ptr::null_mut();
            let mut diag = sys::myb_diag::default();
            assert_eq!(sys::myb_ocm_load(text.as_ptr().cast(), text.len(), ptr::null(), &mut raw, &mut diag), sys::MYB_OK);
            assert_eq!(sys::myb_ocm_state(raw, 0, 1, &mut absolute, &mut epoch, &mut offset, raw_state.as_mut_ptr(), &mut count), sys::MYB_OK);
            let mut raw_values = [0.0; 6];
            let mut n = 0usize;
            assert_eq!(sys::myb_ocm_line(raw, sys::MYB_OCM_TRAJ, 0, 1, &mut absolute, &mut epoch, &mut offset, raw_values.as_mut_ptr(), 6, &mut n), sys::MYB_OK);
            assert_eq!(bits(&values), bits(&raw_values[..n]));
            sys::myb_ocm_destroy(raw);
        }
        assert_eq!(bits(&state), bits(&raw_state[..count as usize]));
        assert!(ocm.line_text(OcmKind::Trajectory, 0, 0).is_some() && ocm.line_text(OcmKind::Trajectory, 0, 2).is_none());
        assert_eq!(ocm.state(0, 2).unwrap_err().status, Status::Domain);
    }
    let without = OcmData { meta: vec![("TIME_SYSTEM".into(), "UTC".into())], ..data.clone() };
    assert_eq!(myb::write_ocm(&without, Notation::Kvn).unwrap_err().status, Status::Domain);
    let labelled = myb::ocm_format_line(Some(&calendar(2, 0.5)), 0.0, &[1.0, 2.0])?;
    assert!(labelled.starts_with("2026-03-01T12:02:00.5"), "{labelled}");
    Ok(())
}

// ---- CK and SCLK ----------------------------------------------------------------------------------

#[test]
fn ck_and_sclk_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    let mut pool = KernelPool::new()?;
    pool.load_file(test_data("ck/test.tsc"))?;
    let sclk = Rc::new(Sclk::new(&pool, 9099)?);
    let (raw_pool, raw_sclk) = unsafe {
        // SAFETY: the calls of the C interface with valid arguments
        let mut raw_pool = ptr::null_mut();
        assert_eq!(sys::myb_kpool_create(ptr::null(), &mut raw_pool), sys::MYB_OK);
        let path = CString::new(test_data("ck/test.tsc").to_str().unwrap()).unwrap();
        let mut diag = sys::myb_diag::default();
        assert_eq!(sys::myb_kpool_load_file(raw_pool, path.as_ptr(), &mut diag), sys::MYB_OK);
        let mut raw_sclk = ptr::null_mut();
        assert_eq!(sys::myb_sclk_create(raw_pool, 9099, ptr::null(), &mut raw_sclk, &mut diag), sys::MYB_OK);
        (raw_pool, raw_sclk)
    };
    let ticks = sclk.encode("1/1000.128.3")?;
    let t = sclk.to_time(ticks)?;
    let back = sclk.from_time(t)?;
    let text = sclk.decode(ticks)?;
    // SAFETY: the clock is valid
    unsafe {
        let mut raw_ticks = 0.0;
        assert_eq!(sys::myb_sclk_encode(raw_sclk, c"1/1000.128.3".as_ptr(), 12, &mut raw_ticks), sys::MYB_OK);
        assert_eq!(ticks.to_bits(), raw_ticks.to_bits());
        let mut raw_t = sys::myb_time::default();
        assert_eq!(sys::myb_sclk_to_time(raw_sclk, raw_ticks, &mut raw_t), sys::MYB_OK);
        same(&t, &raw_t);
        let mut raw_back = 0.0;
        assert_eq!(sys::myb_sclk_from_time(raw_sclk, raw_t, &mut raw_back), sys::MYB_OK);
        assert_eq!(back.to_bits(), raw_back.to_bits());
        assert_eq!(text, raw_text(|b, s, l| sys::myb_sclk_decode(raw_sclk, raw_ticks, b, s, l)));
        assert_eq!((sclk.time_system(), sclk.partitions()), (sys::myb_sclk_time_system(raw_sclk), sys::myb_sclk_partitions(raw_sclk)));
    }
    assert_eq!(sclk.encode("not a clock").unwrap_err().status, Status::Format);
    assert_eq!(sclk.to_time(-1.0e12).unwrap_err().status, Status::Range);
    assert_eq!(Sclk::new(&pool, 1234).err().unwrap().status, Status::NotFound);

    // a CK file read both ways
    let mut ck = Ck::new()?;
    ck.load_file(test_data("ck/test.bc"))?;
    let segments = ck.segments();
    let first = segments[0];
    let middle = 0.5 * (first.start + first.stop);
    let pointing = ck.pointing_at(first.instrument, middle, 1.0e6, false)?;
    // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
    let raw_pointing = unsafe {
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_ck_create(ptr::null(), &mut raw), sys::MYB_OK);
        let path = CString::new(test_data("ck/test.bc").to_str().unwrap()).unwrap();
        let mut diag = sys::myb_diag::default();
        assert_eq!(sys::myb_ck_load_file(raw, path.as_ptr(), &mut diag), sys::MYB_OK);
        assert_eq!(sys::myb_ck_segment_count(raw), segments.len());
        let mut pointing = sys::myb_ck_pointing::default();
        assert_eq!(sys::myb_ck_pointing_at(raw, first.instrument, middle, 1.0e6, 0, &mut pointing), sys::MYB_OK);
        sys::myb_ck_destroy(raw);
        pointing
    };
    same(&pointing, &raw_pointing);
    assert_eq!(ck.pointing_at(123456, middle, 0.0, false).unwrap_err().status, Status::Range);
    assert_eq!(ck.pointing_at(first.instrument, middle, -1.0, false).unwrap_err().status, Status::Domain);

    // the orientation of the instrument as a provider
    let ck = Rc::new(ck);
    let orientation = myb::CkOrientation::new(&ck, &sclk, first.instrument, 0, 1.0e6);
    let at = sclk.to_time(middle)?;
    let o = myb::orientation_from(&orientation, at)?;
    let raw_context = sys::myb_ck_orientation { ck: ck.as_ptr(), sclk: raw_sclk, instrument: first.instrument, base: 0, tolerance: 1.0e6 };
    let mut raw_o = sys::myb_orientation::default();
    // SAFETY: the context points to live objects
    assert_eq!(unsafe { sys::myb_orientation_ck((&raw_context as *const sys::myb_ck_orientation as *mut sys::myb_ck_orientation).cast(), at, &mut raw_o) }, sys::MYB_OK);
    same(&o, &raw_o);
    // SAFETY: created above, destroyed once
    unsafe {
        sys::myb_sclk_destroy(raw_sclk);
        sys::myb_kpool_destroy(raw_pool);
    }

    // a CK written and read again
    let ticks: Vec<f64> = (0..8).map(|k| 1000.0 + 100.0 * k as f64).collect();
    let quaternions: Vec<[f64; 4]> = (0..8).map(|k| {
        let a = 0.05 * k as f64;
        [a.cos(), 0.0, 0.0, a.sin()]
    }).collect();
    let header = CkSegmentHeader { instrument: -9099001, frame: 1, name: "WRITTEN".into() };
    let mut writer = CkWriter::new(Some("myb test"))?;
    writer.add_pointings(&header, 3, &ticks, &quaternions, None, None)?;
    let rates = vec![[0.0, 0.0, 1e-3]; 2];
    writer.add_intervals(&CkSegmentHeader { instrument: -9099002, ..header.clone() }, &[0.0, 500.0], &[400.0, 900.0], &quaternions[..2], &rates, &[1.0, 1.0])?;
    let packets = CkPackets { subtype: 1, window: 4, seconds_per_tick: 1.0, ticks: ticks.clone(), packets: quaternions.concat() };
    writer.add_packets(&CkSegmentHeader { instrument: -9099003, ..header.clone() }, &packets, false, 1000.0, 1700.0, None)?;
    let image = writer.write()?;
    let raw_image = unsafe {
        // SAFETY: the calls of the C interface with valid arguments; everything is destroyed
        let mut raw = ptr::null_mut();
        assert_eq!(sys::myb_ck_writer_create(c"myb test".as_ptr(), ptr::null(), &mut raw), sys::MYB_OK);
        let mut made = sys::myb_ck_segment_header { instrument: -9099001, frame: 1, name: c"WRITTEN".as_ptr() };
        assert_eq!(sys::myb_ck_writer_add_pointings(raw, &made, 3, ticks.as_ptr(), quaternions.as_ptr().cast(), ptr::null(), 8, ptr::null(), 0), sys::MYB_OK);
        made.instrument = -9099002;
        let (starts, stops, spt) = ([0.0, 500.0], [400.0, 900.0], [1.0, 1.0]);
        assert_eq!(sys::myb_ck_writer_add_intervals(raw, &made, starts.as_ptr(), stops.as_ptr(), quaternions.as_ptr().cast(), rates.as_ptr().cast(), spt.as_ptr(), 2), sys::MYB_OK);
        made.instrument = -9099003;
        let flat = quaternions.concat();
        let raw_packets = sys::myb_ck_packets { subtype: 1, window: 4, seconds_per_tick: 1.0, ticks: ticks.as_ptr(), packets: flat.as_ptr(), count: 8 };
        assert_eq!(sys::myb_ck_writer_add_packets(raw, &made, &raw_packets, 0, 1000.0, 1700.0, ptr::null(), 0), sys::MYB_OK);
        let mut image = vec![0u8; sys::myb_ck_writer_size(raw)];
        assert_eq!(sys::myb_ck_writer_write(raw, image.as_mut_ptr().cast(), image.len()), sys::MYB_OK);
        sys::myb_ck_writer_destroy(raw);
        image
    };
    assert!(image == raw_image, "the files differ");
    let mut written = Ck::new()?;
    written.load_memory(image)?;
    assert_eq!(written.segments().len(), 3);
    let read = written.pointing_at(-9099001, 1250.0, 0.0, false)?;
    assert!((read.q[0].abs() - (0.05 * 2.5_f64).cos()).abs() < 1e-6, "{:?}", read.q);
    // failures
    let wrong = CkPackets { subtype: 7, ..packets.clone() };
    assert_eq!(writer.add_packets(&header, &wrong, false, 1000.0, 1700.0, None).unwrap_err().status, Status::Domain);
    let short = CkPackets { packets: packets.packets[..30].to_vec(), ..packets };
    assert_eq!(writer.add_packets(&header, &short, false, 1000.0, 1700.0, None).unwrap_err().status, Status::Domain);
    assert_eq!(writer.add_pointings(&header, 3, &ticks[..3], &quaternions, None, None).unwrap_err().status, Status::Domain);
    assert_eq!(writer.add_chebyshev(&header, &[0.0], &[1.0], &[1; 6], &[1.0; 6], false).unwrap_err().status, Status::Domain);
    assert_eq!(CkWriter::new(None)?.write().unwrap_err().status, Status::Domain);
    assert_eq!(Ck::new()?.load_file(test_data("ck/test.tsc")).unwrap_err().status, Status::Format);
    Ok(())
}

// ---- Two-line element sets ------------------------------------------------------------------------

#[cfg(feature = "sgp4")]
#[test]
fn element_sets_give_the_bits_of_the_c_calls() -> myb::Result<()> {
    use myb::sgp4::{self, Sgp4, Sgp4Code, Sgp4Gravity, Sgp4Mode, TleFitOptions};
    use sys::sgp4 as ext;
    let line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    let line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
    let tle = sgp4::tle_parse(line1, line2, true)?;
    let mut raw = ext::myb_tle::default();
    let mut diag = sys::myb_diag::default();
    let (c1, c2) = (CString::new(line1).unwrap(), CString::new(line2).unwrap());
    // SAFETY: NUL-terminated lines
    assert_eq!(unsafe { ext::myb_tle_parse(c1.as_ptr(), c2.as_ptr(), 0, &mut raw, &mut diag) }, sys::MYB_OK);
    same(&tle, &raw);
    assert_eq!(sgp4::tle_format(&tle)?, (line1.to_string(), line2.to_string()));
    let epoch = sgp4::tle_epoch(&tle);
    assert_eq!((epoch.year, epoch.month, epoch.day), (2000, 6, 27));

    let theory = Sgp4::new(&tle, Sgp4Gravity::Wgs72, Sgp4Mode::Improved)?;
    let (r, v, code) = theory.propagate(360.0)?;
    let (rk, _, _) = theory.propagate_km(360.0)?;
    let (mut raw_r, mut raw_v, mut raw_code) = ([0.0; 3], [0.0; 3], 0);
    // SAFETY: the calls of the C interface with valid arguments
    unsafe {
        let mut record: ext::myb_sgp4 = std::mem::zeroed();
        assert_eq!(ext::myb_sgp4_init(&raw, ext::MYB_SGP4_WGS72, ext::MYB_SGP4_IMPROVED, &mut record), sys::MYB_OK);
        assert_eq!(ext::myb_sgp4_propagate(&record, 360.0, raw_r.as_mut_ptr(), raw_v.as_mut_ptr(), &mut raw_code), sys::MYB_OK);
    }
    assert_eq!((bits(&r), bits(&v), code), (bits(&raw_r), bits(&raw_v), Sgp4Code::Ok));
    assert_eq!(raw_code, ext::MYB_SGP4_OK);
    assert!((rk[0] * 1000.0 - r[0]).abs() < 1e-6);

    // the element set fitted to its own positions
    let minutes: Vec<f64> = (0..12).map(|k| 30.0 * k as f64).collect();
    let positions: Vec<[f64; 3]> = minutes.iter().map(|&m| theory.propagate(m).unwrap().0).collect();
    let (r0, v0, _) = theory.propagate(0.0)?;
    let mut guess = tle;
    sgp4::tle_from_state(&r0, &v0, Sgp4Gravity::Wgs72, &mut guess)?;
    let mut fitted = guess;
    let fit = sgp4::tle_fit(&minutes, &positions, &TleFitOptions::default(), &mut fitted)?;
    let mut raw_fitted = guess;
    let (mut rms, mut iterations) = (0.0, 0);
    // SAFETY: one time and three numbers per position
    assert_eq!(
        unsafe { ext::myb_tle_fit(minutes.as_ptr(), positions.as_ptr().cast(), 12, ptr::null(), &mut raw_fitted, &mut rms, &mut iterations) },
        sys::MYB_OK
    );
    same(&fitted, &raw_fitted);
    assert_eq!((fit.rms.to_bits(), fit.iterations), (rms.to_bits(), iterations));
    assert!(fit.rms < 1.0, "{fit:?}");

    // a mean-elements message and back
    let mut omm = sgp4::omm_from_tle(&tle)?;
    myb::set_odm_text(&mut omm.originator, "MYB")?;
    myb::set_odm_text(&mut omm.object_name, "VANGUARD 1")?;
    let text = myb::write_omm(&omm, Notation::Xml)?;
    same(&sgp4::tle_from_omm(&myb::parse_omm(&text)?)?, &tle);

    // the words of an SPK segment of type 10
    let constants = sgp4::sgp4_constants(Sgp4Gravity::Wgs72);
    let spk10 = sgp4::sgp4_spk10_constants(&constants);
    let mut time_pool = KernelPool::new()?;
    let data = std::env::var_os("MYB_DATA_DIR").map_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../data"), PathBuf::from);
    for kernel in ["time/time-scales.tk", "time/tdb-tt-fb1990.tk", "time/utc-tai.tk"] {
        time_pool.load_file(data.join(kernel))?;
    }
    let system = myb::TimeSystem::new(&time_pool)?;
    let packet = sgp4::sgp4_spk10_packet(&tle, &system, None)?;
    let words = sgp4::sgp4_spk10_words(&spk10, &[packet])?;
    assert_eq!(words.len(), sgp4::sgp4_spk10_size(1));
    let (mut raw_constants, mut raw_packet) = ([0.0; 8], [0.0; 14]);
    let mut raw_words = vec![0.0; words.len()];
    // SAFETY: arrays of the sizes the C interface asks for
    unsafe {
        ext::myb_sgp4_spk10_constants(&constants, raw_constants.as_mut_ptr());
        assert_eq!(ext::myb_sgp4_spk10_packet(&raw, system.as_ptr(), ptr::null(), raw_packet.as_mut_ptr()), sys::MYB_OK);
        assert_eq!(ext::myb_sgp4_spk10_words(raw_constants.as_ptr(), raw_packet.as_ptr(), 1, raw_words.as_mut_ptr(), raw_words.len()), sys::MYB_OK);
        assert_eq!(sgp4::sgp4_gmst(2451545.0, 0.25).to_bits(), ext::myb_sgp4_gmst(2451545.0, 0.25).to_bits());
    }
    assert_eq!((bits(&spk10), bits(&packet), bits(&words)), (bits(&raw_constants), bits(&raw_packet), bits(&raw_words)));

    // failures
    let error = sgp4::tle_parse(line1, "2 00005  34.2682", true).err().unwrap();
    assert_eq!(error.status, Status::Format);
    let edited = format!("{}0", &line1[..68]);
    assert_eq!(sgp4::tle_parse(&edited, line2, true).err().unwrap().status, Status::Format);
    sgp4::tle_parse(&edited, line2, false)?;
    let mut hyperbolic = tle;
    hyperbolic.eccentricity = 1.5;
    assert_eq!(Sgp4::new(&hyperbolic, Sgp4Gravity::Wgs72, Sgp4Mode::Improved).err().unwrap().status, Status::Domain);
    assert_eq!(sgp4::tle_fit(&minutes[..2], &positions[..2], &TleFitOptions::default(), &mut fitted).unwrap_err().status, Status::Domain);
    assert_eq!(sgp4::sgp4_spk10_words(&spk10, &[]).unwrap_err().status, Status::Domain);
    Ok(())
}
