// SPDX-License-Identifier: BSD-3-Clause
//! Files (`myb/io.h`, `myb/kin.h`, `myb/ccsds.h`, `myb/ck.h`): numbers read and printed
//! exactly, tables of numbers, the variables of a kernel pool, the SPK writer, SP3 precise
//! orbits, the orbit data messages OPM, OMM and OCM (the OEM is at the root of the crate),
//! CK orientation files with their writer, and spacecraft clocks.
use std::ffi::{c_char, c_int, CString};
use std::path::Path;
use std::ptr::{self, NonNull};

use crate::{
    array_text, c_path, c_string, c_text, check, check_diag, created, domain, flat, set_array_text, sys, Calendar, Error,
    KernelPool, Kept, Notation, Oem, Result, State, StateSource, Status, Time, Vec3,
};

/// A text of a writer of the library, with a buffer first and a larger one if it said so.
pub(crate) fn written(write: impl Fn(*mut c_char, usize, *mut usize) -> sys::myb_status) -> Result<String> {
    let mut buffer = vec![0u8; 256];
    let mut length = 0usize;
    let code = write(buffer.as_mut_ptr().cast(), buffer.len(), &mut length);
    if code == sys::MYB_E_OVERFLOW {
        buffer = vec![0u8; length + 1];
        check(write(buffer.as_mut_ptr().cast(), buffer.len(), &mut length))?;
    } else {
        check(code)?;
    }
    buffer.truncate(length);
    String::from_utf8(buffer).map_err(|_| Error::new(Status::Format))
}

fn notation_of(notation: Notation) -> sys::myb_odm_notation {
    match notation {
        Notation::Kvn => sys::MYB_ODM_KVN,
        Notation::Xml => sys::MYB_ODM_XML,
    }
}

// ---- Numbers --------------------------------------------------------------------------------------

/// The number at the start of the text, correctly rounded, and the characters it took.
pub fn parse_double(text: &str) -> Result<(f64, usize)> {
    let (mut value, mut consumed) = (0.0, 0usize);
    // SAFETY: the text is text.len() bytes
    check(unsafe { sys::myb_parse_double(text.as_ptr().cast(), text.len(), &mut value, &mut consumed) })?;
    Ok((value, consumed))
}

/// The shortest decimal that reads back to the same value ("0.001", "1e+300").
pub fn format_double(value: f64) -> String {
    let mut out = [0u8; sys::MYB_DOUBLE_CHARS as usize];
    // SAFETY: out holds MYB_DOUBLE_CHARS characters
    let length = unsafe { sys::myb_format_double(value, out.as_mut_ptr().cast()) };
    String::from_utf8_lossy(&out[..length.min(out.len())]).into_owned()
}

/// `decimals` digits after the point, as "%.Nf", exactly rounded.
pub fn format_fixed(value: f64, decimals: u32) -> Result<String> {
    let decimals = c_int::try_from(decimals).or_else(|_| domain())?;
    // SAFETY: the buffers are as large as said
    written(|buffer, size, length| unsafe { sys::myb_format_fixed(value, decimals, buffer, size, length) })
}

/// `digits` significant digits with an exponent, as "%.(digits-1)e".
pub fn format_scientific(value: f64, digits: u32) -> Result<String> {
    let digits = c_int::try_from(digits).or_else(|_| domain())?;
    // SAFETY: the buffers are as large as said
    written(|buffer, size, length| unsafe { sys::myb_format_scientific(value, digits, buffer, size, length) })
}

// ---- Tables ---------------------------------------------------------------------------------------

/// Rows of numbers read from CSV-like text (`myb_table`).
pub struct Table {
    raw: NonNull<sys::myb_table>,
}

impl Table {
    pub fn parse(text: &str) -> Result<Self> {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the text is text.len() bytes
        let code = unsafe { sys::myb_table_load(text.as_ptr().cast(), text.len(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Table { raw: created(raw)? })
    }

    pub fn load_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = c_path(path.as_ref())?;
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the path is NUL-terminated
        let code = unsafe { sys::myb_table_load_file(path.as_ptr(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Table { raw: created(raw)? })
    }

    pub fn rows(&self) -> usize {
        // SAFETY: the table is alive
        unsafe { sys::myb_table_rows(self.raw.as_ptr()) }
    }

    pub fn columns(&self) -> usize {
        // SAFETY: the table is alive
        unsafe { sys::myb_table_columns(self.raw.as_ptr()) }
    }

    /// The name of a column ("" without a row of names).
    pub fn name(&self, column: usize) -> Option<String> {
        // SAFETY: the column is one of the table
        (column < self.columns()).then(|| c_text(unsafe { sys::myb_table_name(self.raw.as_ptr(), column) }).unwrap_or_default())
    }

    /// The values of a row.
    pub fn row(&self, row: usize) -> Option<&[f64]> {
        if row >= self.rows() {
            return None;
        }
        // SAFETY: the row is one of the table and holds `columns` numbers that live with it
        Some(unsafe { std::slice::from_raw_parts(sys::myb_table_row(self.raw.as_ptr(), row), self.columns()) })
    }
}

impl Drop for Table {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_table_destroy(self.raw.as_ptr()) }
    }
}

/// Rows of `columns` values (row by row) as comma-separated text, with a first row of names.
pub fn write_table(names: Option<&[&str]>, columns: usize, values: &[f64]) -> Result<String> {
    if columns == 0 || values.len() % columns != 0 || names.is_some_and(|n| n.len() != columns) {
        return domain();
    }
    let names: Option<Vec<CString>> = names.map(|n| n.iter().map(|&name| c_string(name)).collect()).transpose()?;
    let pointers: Option<Vec<*const c_char>> = names.as_ref().map(|n| n.iter().map(|name| name.as_ptr()).collect());
    let rows = values.len() / columns;
    // SAFETY: the names and values live for the call; the buffers are as large as said
    written(|buffer, size, length| unsafe {
        sys::myb_table_write(
            pointers.as_ref().map_or(ptr::null(), |p| p.as_ptr()),
            columns,
            values.as_ptr(),
            rows,
            buffer,
            size,
            length,
        )
    })
}

// ---- Kernel variables -----------------------------------------------------------------------------

/// What a variable of a kernel pool holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelVariable {
    Numbers(usize),
    Strings(usize),
}

impl KernelPool {
    /// Loads kernel text (atomic: on an error the pool is unchanged).
    pub fn load_text(&mut self, text: &str) -> Result<()> {
        let mut diag = sys::myb_diag::default();
        // SAFETY: the text is text.len() bytes
        let code = unsafe { sys::myb_kpool_load_text(self.as_ptr(), text.as_ptr().cast(), text.len(), &mut diag) };
        check_diag(code, &diag)
    }

    /// The type and the number of values of a variable.
    pub fn describe(&self, name: &str) -> Result<KernelVariable> {
        let name = c_string(name)?;
        let (mut kind, mut count) = (0, 0usize);
        // SAFETY: the name is NUL-terminated
        check(unsafe { sys::myb_kpool_describe(self.as_ptr(), name.as_ptr(), &mut kind, &mut count) })?;
        Ok(if kind == sys::MYB_KVAR_STRINGS { KernelVariable::Strings(count) } else { KernelVariable::Numbers(count) })
    }

    /// The numbers of a variable.
    pub fn numbers(&self, name: &str) -> Result<Vec<f64>> {
        let name = c_string(name)?;
        let (mut values, mut count) = (ptr::null(), 0usize);
        // SAFETY: the name is NUL-terminated; the values live until the next load, and are copied
        check(unsafe { sys::myb_kpool_numbers(self.as_ptr(), name.as_ptr(), &mut values, &mut count) })?;
        Ok(if count == 0 { Vec::new() } else { unsafe { std::slice::from_raw_parts(values, count) }.to_vec() })
    }

    /// One string of a variable.
    pub fn string(&self, name: &str, index: usize) -> Result<String> {
        let name = c_string(name)?;
        let (mut value, mut length) = (ptr::null(), 0usize);
        // SAFETY: the name is NUL-terminated; the value lives until the next load, and is copied
        check(unsafe { sys::myb_kpool_string(self.as_ptr(), name.as_ptr(), index, &mut value, &mut length) })?;
        let bytes = if length == 0 { &[][..] } else { unsafe { std::slice::from_raw_parts(value.cast::<u8>(), length) } };
        Ok(String::from_utf8_lossy(bytes).into_owned())
    }
}

// ---- Writing SPK files ----------------------------------------------------------------------------

/// What a segment describes: target relative to center in a frame (NAIF codes), and its name
/// (1 ..= 40 printable ASCII characters).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SpkSegmentHeader {
    pub target: i32,
    pub center: i32,
    pub frame: i32,
    pub name: String,
}

/// Collects SPK segments in memory and lays them out as one file (`myb_spk_writer`).
pub struct SpkWriter {
    raw: NonNull<sys::myb_spk_writer>,
}

fn spk_header(header: &SpkSegmentHeader) -> Result<(sys::myb_spk_segment_header, CString)> {
    let name = c_string(&header.name)?;
    let made = sys::myb_spk_segment_header { target: header.target, center: header.center, frame: header.frame, name: name.as_ptr() };
    Ok((made, name))
}

impl SpkWriter {
    /// internal_name: the label of the file (`None`: "myb").
    pub fn new(internal_name: Option<&str>) -> Result<Self> {
        let name = internal_name.map(c_string).transpose()?;
        let mut raw = ptr::null_mut();
        // SAFETY: the name, if any, is NUL-terminated
        check(unsafe { sys::myb_spk_writer_create(name.as_ref().map_or(ptr::null(), |n| n.as_ptr()), ptr::null(), &mut raw) })?;
        Ok(SpkWriter { raw: created(raw)? })
    }

    /// A segment of states at increasing epochs (s past J2000, TDB): type 9 (Lagrange) or 13
    /// (Hermite) of the degree.
    pub fn add_states(&mut self, header: &SpkSegmentHeader, kind: i32, degree: u32, epochs: &[f64], states: &[State]) -> Result<()> {
        if epochs.len() != states.len() {
            return domain();
        }
        let degree = c_int::try_from(degree).or_else(|_| domain())?;
        let (made, _name) = spk_header(header)?;
        // SAFETY: the arrays hold one epoch and six numbers per state
        check(unsafe {
            sys::myb_spk_writer_add_states(self.raw.as_ptr(), &made, kind, degree, epochs.as_ptr(), flat(states).as_ptr(), states.len())
        })
    }

    /// A segment of Chebyshev polynomials (type 2 or 3) fitted to a provider on `intervals`
    /// parts of [start, stop]; with `measure_error`, the largest distance (m) between the
    /// provider and the polynomials half-way between the nodes.
    #[allow(clippy::too_many_arguments)]
    pub fn add_chebyshev(
        &mut self,
        header: &SpkSegmentHeader,
        kind: i32,
        state: &(impl StateSource + ?Sized),
        start: f64,
        stop: f64,
        intervals: usize,
        degree: u32,
        measure_error: bool,
    ) -> Result<Option<f64>> {
        let degree = c_int::try_from(degree).or_else(|_| domain())?;
        let (made, _name) = spk_header(header)?;
        let (function, context) = state.as_raw();
        let mut max_error = 0.0;
        // SAFETY: the provider lives for the call
        check(unsafe {
            sys::myb_spk_writer_add_chebyshev(
                self.raw.as_ptr(),
                &made,
                kind,
                function,
                context,
                start,
                stop,
                intervals,
                degree,
                if measure_error { &mut max_error } else { ptr::null_mut() },
            )
        })?;
        Ok(measure_error.then_some(max_error))
    }

    /// A segment of any type from its words as a file holds them.
    pub fn add_raw(&mut self, header: &SpkSegmentHeader, kind: i32, start: f64, stop: f64, words: &[f64]) -> Result<()> {
        let (made, _name) = spk_header(header)?;
        // SAFETY: words holds words.len() numbers
        check(unsafe { sys::myb_spk_writer_add_raw(self.raw.as_ptr(), &made, kind, start, stop, words.as_ptr(), words.len()) })
    }

    /// Bytes of the file so far.
    pub fn size(&self) -> usize {
        // SAFETY: the writer is alive
        unsafe { sys::myb_spk_writer_size(self.raw.as_ptr()) }
    }

    /// The file in memory.
    pub fn write(&self) -> Result<Vec<u8>> {
        let mut image = vec![0u8; self.size()];
        // SAFETY: image holds `size` bytes
        check(unsafe { sys::myb_spk_writer_write(self.raw.as_ptr(), image.as_mut_ptr().cast(), image.len()) })?;
        Ok(image)
    }

    pub fn write_file(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = c_path(path.as_ref())?;
        // SAFETY: the path is NUL-terminated
        check(unsafe { sys::myb_spk_writer_write_file(self.raw.as_ptr(), path.as_ptr()) })
    }
}

impl Drop for SpkWriter {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_spk_writer_destroy(self.raw.as_ptr()) }
    }
}

// ---- SP3 ------------------------------------------------------------------------------------------

/// SP3-c and SP3-d precise orbits (`myb_sp3`).
pub struct Sp3 {
    raw: NonNull<sys::myb_sp3>,
}

/// What a record gives of a satellite at an epoch: 1 position, 2 velocity (flags), and both.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sp3State {
    pub given: i32,
    pub position: Vec3,
    pub velocity: Vec3,
}

/// The clock of a record: what is given (1 the clock, 2 the rate), the clock (s) and its rate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sp3Clock {
    pub given: i32,
    pub clock: f64,
    pub rate: f64,
}

impl Sp3 {
    pub fn parse(text: &str) -> Result<Self> {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the text is text.len() bytes
        let code = unsafe { sys::myb_sp3_load(text.as_ptr().cast(), text.len(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Sp3 { raw: created(raw)? })
    }

    pub fn load_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = c_path(path.as_ref())?;
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the path is NUL-terminated
        let code = unsafe { sys::myb_sp3_load_file(path.as_ptr(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Sp3 { raw: created(raw)? })
    }

    pub fn satellite_count(&self) -> usize {
        // SAFETY: the file is alive
        unsafe { sys::myb_sp3_satellite_count(self.raw.as_ptr()) }
    }

    /// The identifier of a satellite ("L51").
    pub fn satellite(&self, index: usize) -> Option<String> {
        // SAFETY: the index is one of the file
        if index < self.satellite_count() {
            c_text(unsafe { sys::myb_sp3_satellite(self.raw.as_ptr(), index) })
        } else {
            None
        }
    }

    pub fn epoch_count(&self) -> usize {
        // SAFETY: the file is alive
        unsafe { sys::myb_sp3_epoch_count(self.raw.as_ptr()) }
    }

    pub fn time_system(&self) -> String {
        // SAFETY: the file is alive
        c_text(unsafe { sys::myb_sp3_time_system(self.raw.as_ptr()) }).unwrap_or_default()
    }

    pub fn frame(&self) -> String {
        // SAFETY: the file is alive
        c_text(unsafe { sys::myb_sp3_frame(self.raw.as_ptr()) }).unwrap_or_default()
    }

    fn record(&self, satellite: usize, epoch: usize) -> bool {
        satellite < self.satellite_count() && epoch < self.epoch_count()
    }

    pub fn epoch(&self, epoch: usize) -> Option<Calendar> {
        if epoch >= self.epoch_count() {
            return None;
        }
        let mut label = Calendar::default();
        // SAFETY: the epoch is one of the file
        unsafe { sys::myb_sp3_epoch(self.raw.as_ptr(), epoch, &mut label) };
        Some(label)
    }

    pub fn state(&self, satellite: usize, epoch: usize) -> Option<Sp3State> {
        if !self.record(satellite, epoch) {
            return None;
        }
        let mut out = Sp3State { given: 0, position: [0.0; 3], velocity: [0.0; 3] };
        // SAFETY: the record is one of the file
        out.given = unsafe { sys::myb_sp3_state(self.raw.as_ptr(), satellite, epoch, out.position.as_mut_ptr(), out.velocity.as_mut_ptr()) };
        Some(out)
    }

    pub fn clock(&self, satellite: usize, epoch: usize) -> Option<Sp3Clock> {
        if !self.record(satellite, epoch) {
            return None;
        }
        let mut out = Sp3Clock { given: 0, clock: 0.0, rate: 0.0 };
        // SAFETY: the record is one of the file
        out.given = unsafe { sys::myb_sp3_clock(self.raw.as_ptr(), satellite, epoch, &mut out.clock, &mut out.rate) };
        Some(out)
    }

    /// The standard deviations of a record: position (3), clock, velocity (3), clock rate.
    pub fn sigma(&self, satellite: usize, epoch: usize) -> Option<[f64; 8]> {
        if !self.record(satellite, epoch) {
            return None;
        }
        let mut sigma = [0.0; 8];
        // SAFETY: the record is one of the file
        unsafe { sys::myb_sp3_sigma(self.raw.as_ptr(), satellite, epoch, sigma.as_mut_ptr()) };
        Some(sigma)
    }

    /// The accuracy of the orbit of a satellite the header states (m), 0 if unknown.
    pub fn accuracy(&self, satellite: usize) -> Option<f64> {
        // SAFETY: the satellite is one of the file
        (satellite < self.satellite_count()).then(|| unsafe { sys::myb_sp3_accuracy(self.raw.as_ptr(), satellite) })
    }

    /// MYB_SP3_* flags of a record.
    pub fn flags(&self, satellite: usize, epoch: usize) -> Option<u32> {
        // SAFETY: the record is one of the file
        self.record(satellite, epoch).then(|| unsafe { sys::myb_sp3_flags(self.raw.as_ptr(), satellite, epoch) })
    }
}

impl Drop for Sp3 {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_sp3_destroy(self.raw.as_ptr()) }
    }
}

/// What an SP3 file is written from (`myb_sp3_data`). Records are [epoch * satellites +
/// satellite]: positions and velocities 3 each, sigmas 8 each, clocks one each.
#[derive(Debug, Clone, Default)]
pub struct Sp3Data {
    /// 'c' or 'd'.
    pub version: u8,
    /// 'G', 'M', 'L', ...
    pub file_type: u8,
    pub time_system: String,
    pub frame: String,
    pub orbit_type: String,
    pub agency: String,
    pub data_used: String,
    pub satellites: Vec<String>,
    pub epochs: Vec<Calendar>,
    /// s between epochs.
    pub interval: f64,
    pub positions: Vec<f64>,
    pub velocities: Option<Vec<f64>>,
    pub clocks: Option<Vec<f64>>,
    pub clock_rates: Option<Vec<f64>>,
    pub accuracy: Option<Vec<f64>>,
    pub sigmas: Option<Vec<f64>>,
    pub flags: Option<Vec<u8>>,
}

/// An SP3 file as text.
pub fn write_sp3(data: &Sp3Data) -> Result<String> {
    let records = data.satellites.len() * data.epochs.len();
    let sized = |values: &Option<Vec<f64>>, each: usize| values.as_ref().is_none_or(|v| v.len() == each * records);
    if data.positions.len() != 3 * records
        || !sized(&data.velocities, 3)
        || !sized(&data.clocks, 1)
        || !sized(&data.clock_rates, 1)
        || !sized(&data.sigmas, 8)
        || data.accuracy.as_ref().is_some_and(|a| a.len() != data.satellites.len())
        || data.flags.as_ref().is_some_and(|f| f.len() != records)
    {
        return domain();
    }
    let texts = [&data.time_system, &data.frame, &data.orbit_type, &data.agency, &data.data_used]
        .map(|text| c_string(text))
        .into_iter()
        .collect::<Result<Vec<_>>>()?;
    let satellites = data.satellites.iter().map(|s| c_string(s)).collect::<Result<Vec<_>>>()?;
    let pointers: Vec<*const c_char> = satellites.iter().map(|s| s.as_ptr()).collect();
    let optional = |values: &Option<Vec<f64>>| values.as_ref().map_or(ptr::null(), |v| v.as_ptr());
    let made = sys::myb_sp3_data {
        version: data.version as c_char,
        file_type: data.file_type as c_char,
        time_system: texts[0].as_ptr(),
        frame: texts[1].as_ptr(),
        orbit_type: texts[2].as_ptr(),
        agency: texts[3].as_ptr(),
        data_used: texts[4].as_ptr(),
        satellite_count: pointers.len(),
        satellites: pointers.as_ptr(),
        epoch_count: data.epochs.len(),
        epochs: data.epochs.as_ptr(),
        interval: data.interval,
        positions: data.positions.as_ptr(),
        velocities: optional(&data.velocities),
        clocks: optional(&data.clocks),
        clock_rates: optional(&data.clock_rates),
        accuracy: optional(&data.accuracy),
        sigmas: optional(&data.sigmas),
        flags: data.flags.as_ref().map_or(ptr::null(), |f| f.as_ptr()),
    };
    // SAFETY: every pointer of `made` is into `data`, `texts` and `satellites`, alive here
    written(|buffer, size, length| unsafe { sys::myb_sp3_write(&made, buffer, size, length) })
}

// ---- OEM: what the reader of the root does not give ------------------------------------------------

/// A covariance of an OEM segment: its epoch, its frame (`None`: that of the segment) and the
/// matrix (rows; m^2, m^2/s, m^2/s^2).
#[derive(Debug, Clone)]
pub struct OemCovariance {
    pub epoch: Calendar,
    pub ref_frame: Option<String>,
    pub matrix: [f64; 36],
}

impl Oem {
    pub fn load_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = c_path(path.as_ref())?;
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the path is NUL-terminated
        let code = unsafe { sys::myb_oem_load_file(path.as_ptr(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Oem { raw: created(raw)? })
    }

    /// A keyword of the header, as written.
    pub fn header(&self, key: &str) -> Option<String> {
        let key = CString::new(key).ok()?;
        // SAFETY: the key is NUL-terminated; a value lives with the message
        c_text(unsafe { sys::myb_oem_header(self.raw.as_ptr(), key.as_ptr()) })
    }

    /// The epoch, the state and, if the line gives it, the acceleration of a line.
    pub fn state_with_acceleration(&self, segment: usize, index: usize) -> Option<(Calendar, State, Option<Vec3>)> {
        if index >= self.state_count(segment) {
            return None;
        }
        let (mut epoch, mut y, mut a) = (Calendar::default(), [0.0; 6], [0.0; 3]);
        // SAFETY: the line is one of the segment
        let given = unsafe { sys::myb_oem_state(self.raw.as_ptr(), segment, index, &mut epoch, y.as_mut_ptr(), a.as_mut_ptr()) };
        Some((epoch, y, (given == 3).then_some(a)))
    }

    pub fn covariance_count(&self, segment: usize) -> usize {
        if segment >= self.segment_count() {
            return 0;
        }
        // SAFETY: the segment is one of the message
        unsafe { sys::myb_oem_covariance_count(self.raw.as_ptr(), segment) }
    }

    pub fn covariance(&self, segment: usize, index: usize) -> Option<OemCovariance> {
        if index >= self.covariance_count(segment) {
            return None;
        }
        let (mut epoch, mut frame, mut matrix) = (Calendar::default(), ptr::null(), [0.0; 36]);
        // SAFETY: the covariance is one of the segment; the frame lives with the message
        unsafe { sys::myb_oem_covariance(self.raw.as_ptr(), segment, index, &mut epoch, &mut frame, matrix.as_mut_ptr()) };
        Some(OemCovariance { epoch, ref_frame: c_text(frame), matrix })
    }
}

// ---- OPM and OMM ----------------------------------------------------------------------------------

/// A parameter message: the state at an epoch, elements, spacecraft, manoeuvres, covariance;
/// the texts are fixed arrays ([`odm_text`], [`set_odm_text`]).
pub type Opm = sys::myb_opm;
/// A manoeuvre of an OPM.
pub type OpmManeuver = sys::myb_opm_maneuver;
/// A mean-elements message, in the units of the message.
pub type Omm = sys::myb_omm;

/// A text field of a message.
pub fn odm_text(field: &[c_char]) -> String {
    array_text(field)
}

/// Sets a text field of a message (Domain if it does not fit with its NUL).
pub fn set_odm_text(field: &mut [c_char], text: &str) -> Result<()> {
    set_array_text(field, text)
}

/// An OPM read from text in either notation.
pub fn parse_opm(text: &str) -> Result<Opm> {
    let mut opm = Opm::default();
    let mut diag = sys::myb_diag::default();
    // SAFETY: the text is text.len() bytes
    let code = unsafe { sys::myb_opm_parse(text.as_ptr().cast(), text.len(), &mut opm, &mut diag) };
    check_diag(code, &diag)?;
    Ok(opm)
}

/// An OPM as text (`myb_opm_write` for KVN, `myb_opm_write_as` for XML).
pub fn write_opm(opm: &Opm, notation: Notation) -> Result<String> {
    // SAFETY: the buffers are as large as said
    written(|buffer, size, length| unsafe {
        match notation {
            Notation::Kvn => sys::myb_opm_write(opm, buffer, size, length),
            Notation::Xml => sys::myb_opm_write_as(opm, notation_of(notation), buffer, size, length),
        }
    })
}

/// An OMM read from text in either notation.
pub fn parse_omm(text: &str) -> Result<Omm> {
    let mut omm = Omm::default();
    let mut diag = sys::myb_diag::default();
    // SAFETY: the text is text.len() bytes
    let code = unsafe { sys::myb_omm_parse(text.as_ptr().cast(), text.len(), &mut omm, &mut diag) };
    check_diag(code, &diag)?;
    Ok(omm)
}

/// An OMM as text (`myb_omm_write` for KVN, `myb_omm_write_as` for XML).
pub fn write_omm(omm: &Omm, notation: Notation) -> Result<String> {
    // SAFETY: the buffers are as large as said
    written(|buffer, size, length| unsafe {
        match notation {
            Notation::Kvn => sys::myb_omm_write(omm, buffer, size, length),
            Notation::Xml => sys::myb_omm_write_as(omm, notation_of(notation), buffer, size, length),
        }
    })
}

// ---- OCM ------------------------------------------------------------------------------------------

/// The kinds of blocks of an OCM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OcmKind {
    Trajectory,
    Physical,
    Covariance,
    Maneuver,
    Perturbations,
    OrbitDetermination,
    User,
}

impl OcmKind {
    fn raw(self) -> sys::myb_ocm_kind {
        match self {
            OcmKind::Trajectory => sys::MYB_OCM_TRAJ,
            OcmKind::Physical => sys::MYB_OCM_PHYS,
            OcmKind::Covariance => sys::MYB_OCM_COV,
            OcmKind::Maneuver => sys::MYB_OCM_MAN,
            OcmKind::Perturbations => sys::MYB_OCM_PERT,
            OcmKind::OrbitDetermination => sys::MYB_OCM_OD,
            OcmKind::User => sys::MYB_OCM_USER,
        }
    }
}

/// The time of a line of data: a calendar label, or seconds from EPOCH_TZERO.
#[derive(Debug, Clone, Copy)]
pub enum OcmTime {
    Absolute(Calendar),
    Offset(f64),
}

fn ocm_time(absolute: c_int, epoch: Calendar, offset: f64) -> OcmTime {
    if absolute != 0 {
        OcmTime::Absolute(epoch)
    } else {
        OcmTime::Offset(offset)
    }
}

/// A comprehensive message read from text (`myb_ocm`).
pub struct Ocm {
    raw: NonNull<sys::myb_ocm>,
}

impl Ocm {
    pub fn parse(text: &str) -> Result<Self> {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the text is text.len() bytes
        let code = unsafe { sys::myb_ocm_load(text.as_ptr().cast(), text.len(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Ocm { raw: created(raw)? })
    }

    pub fn load_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = c_path(path.as_ref())?;
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the path is NUL-terminated
        let code = unsafe { sys::myb_ocm_load_file(path.as_ptr(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Ocm { raw: created(raw)? })
    }

    pub fn header(&self, key: &str) -> Option<String> {
        let key = CString::new(key).ok()?;
        // SAFETY: the key is NUL-terminated; a value lives with the message
        c_text(unsafe { sys::myb_ocm_header(self.raw.as_ptr(), key.as_ptr()) })
    }

    pub fn meta(&self, key: &str) -> Option<String> {
        let key = CString::new(key).ok()?;
        // SAFETY: as above
        c_text(unsafe { sys::myb_ocm_meta(self.raw.as_ptr(), key.as_ptr()) })
    }

    pub fn block_count(&self, kind: OcmKind) -> usize {
        // SAFETY: the message is alive
        unsafe { sys::myb_ocm_block_count(self.raw.as_ptr(), kind.raw()) }
    }

    /// A keyword of a block, as written.
    pub fn value(&self, kind: OcmKind, index: usize, key: &str) -> Option<String> {
        if index >= self.block_count(kind) {
            return None;
        }
        let key = CString::new(key).ok()?;
        // SAFETY: the block is one of the message
        c_text(unsafe { sys::myb_ocm_value(self.raw.as_ptr(), kind.raw(), index, key.as_ptr()) })
    }

    pub fn line_count(&self, kind: OcmKind, index: usize) -> usize {
        if index >= self.block_count(kind) {
            return 0;
        }
        // SAFETY: the block is one of the message
        unsafe { sys::myb_ocm_line_count(self.raw.as_ptr(), kind.raw(), index) }
    }

    /// A line of data as text.
    pub fn line_text(&self, kind: OcmKind, index: usize, line: usize) -> Option<String> {
        if line >= self.line_count(kind, index) {
            return None;
        }
        // SAFETY: the line is one of the block
        c_text(unsafe { sys::myb_ocm_line_text(self.raw.as_ptr(), kind.raw(), index, line) })
    }

    /// A line of data as its time and numbers in the units of the message.
    pub fn line(&self, kind: OcmKind, index: usize, line: usize) -> Result<(OcmTime, Vec<f64>)> {
        if line >= self.line_count(kind, index) {
            return domain();
        }
        let (mut absolute, mut epoch, mut offset, mut count) = (0, Calendar::default(), 0.0, 0usize);
        // SAFETY: the line is one of the block; no values are asked for
        check(unsafe {
            sys::myb_ocm_line(self.raw.as_ptr(), kind.raw(), index, line, &mut absolute, &mut epoch, &mut offset, ptr::null_mut(), 0, &mut count)
        })?;
        let mut values = vec![0.0; count];
        // SAFETY: values holds `count` numbers
        check(unsafe {
            sys::myb_ocm_line(
                self.raw.as_ptr(),
                kind.raw(),
                index,
                line,
                &mut absolute,
                &mut epoch,
                &mut offset,
                crate::out_ptr(&mut values),
                count,
                &mut count,
            )
        })?;
        Ok((ocm_time(absolute, epoch, offset), values))
    }

    /// A line of a Cartesian trajectory in m, m/s and m/s^2: 3, 6 or 9 numbers.
    pub fn state(&self, index: usize, line: usize) -> Result<(OcmTime, Vec<f64>)> {
        if line >= self.line_count(OcmKind::Trajectory, index) {
            return domain();
        }
        let (mut absolute, mut epoch, mut offset, mut state, mut count) = (0, Calendar::default(), 0.0, [0.0; 9], 0);
        // SAFETY: the line is one of the block; state holds nine numbers
        check(unsafe {
            sys::myb_ocm_state(self.raw.as_ptr(), index, line, &mut absolute, &mut epoch, &mut offset, state.as_mut_ptr(), &mut count)
        })?;
        Ok((ocm_time(absolute, epoch, offset), state[..(count.max(0) as usize).min(9)].to_vec()))
    }
}

impl Drop for Ocm {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_ocm_destroy(self.raw.as_ptr()) }
    }
}

/// A block of an OCM to write: its keywords and, for trajectories, covariances and
/// manoeuvres, its lines of data ([`ocm_format_line`] makes one).
#[derive(Debug, Clone, PartialEq)]
pub struct OcmBlock {
    pub kind: OcmKind,
    pub pairs: Vec<(String, String)>,
    pub lines: Vec<String>,
}

/// What an OCM is written from; the metadata must give TIME_SYSTEM and EPOCH_TZERO.
#[derive(Debug, Clone)]
pub struct OcmData {
    pub originator: String,
    pub creation_date: Calendar,
    pub meta: Vec<(String, String)>,
    pub blocks: Vec<OcmBlock>,
}

fn pairs_of(pairs: &[(String, String)], texts: &mut Vec<CString>) -> Result<Vec<sys::myb_ocm_pair>> {
    let mut out = Vec::with_capacity(pairs.len());
    for (key, value) in pairs {
        let (key, value) = (c_string(key)?, c_string(value)?);
        out.push(sys::myb_ocm_pair { key: key.as_ptr(), value: value.as_ptr() });
        texts.push(key);
        texts.push(value);
    }
    Ok(out)
}

/// An OCM (version 3.0) as text.
pub fn write_ocm(data: &OcmData, notation: Notation) -> Result<String> {
    let mut texts = Vec::new();
    let originator = c_string(&data.originator)?;
    let meta = pairs_of(&data.meta, &mut texts)?;
    let mut block_pairs = Vec::with_capacity(data.blocks.len());
    let mut block_lines = Vec::with_capacity(data.blocks.len());
    for block in &data.blocks {
        block_pairs.push(pairs_of(&block.pairs, &mut texts)?);
        let lines = block.lines.iter().map(|line| c_string(line)).collect::<Result<Vec<_>>>()?;
        let pointers: Vec<*const c_char> = lines.iter().map(|line| line.as_ptr()).collect();
        texts.extend(lines);
        block_lines.push(pointers);
    }
    let blocks: Vec<sys::myb_ocm_block_data> = data
        .blocks
        .iter()
        .zip(block_pairs.iter().zip(&block_lines))
        .map(|(block, (pairs, lines))| sys::myb_ocm_block_data {
            kind: block.kind.raw(),
            pair_count: pairs.len(),
            pairs: pairs.as_ptr(),
            line_count: lines.len(),
            lines: lines.as_ptr(),
        })
        .collect();
    let made = sys::myb_ocm_data {
        originator: originator.as_ptr(),
        creation_date: data.creation_date,
        meta_count: meta.len(),
        meta: meta.as_ptr(),
        block_count: blocks.len(),
        blocks: blocks.as_ptr(),
    };
    // SAFETY: every pointer of `made` is into the texts, pairs and blocks alive here
    written(|buffer, size, length| unsafe { sys::myb_ocm_write(&made, notation_of(notation), buffer, size, length) })
}

/// A line of data: the time (a label, or `offset` seconds when epoch is `None`) and the values.
pub fn ocm_format_line(epoch: Option<&Calendar>, offset: f64, values: &[f64]) -> Result<String> {
    // SAFETY: values holds values.len() numbers; the buffers are as large as said
    written(|buffer, size, length| unsafe {
        sys::myb_ocm_format_line(
            epoch.map_or(ptr::null(), |e| e as *const Calendar),
            offset,
            values.as_ptr(),
            values.len(),
            buffer,
            size,
            length,
        )
    })
}

// ---- Spacecraft clocks ----------------------------------------------------------------------------

/// A spacecraft clock of type 1 from the kernel variables of its code (`myb_sclk`).
pub struct Sclk {
    raw: NonNull<sys::myb_sclk>,
}

impl Sclk {
    pub fn new(pool: &KernelPool, id: i64) -> Result<Self> {
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the pool is alive for the call
        let code = unsafe { sys::myb_sclk_create(pool.as_ptr(), id, ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(Sclk { raw: created(raw)? })
    }

    /// The scale of the table: 1 TDB, 2 TT.
    pub fn time_system(&self) -> i32 {
        // SAFETY: the clock is alive
        unsafe { sys::myb_sclk_time_system(self.raw.as_ptr()) }
    }

    pub fn partitions(&self) -> usize {
        // SAFETY: the clock is alive
        unsafe { sys::myb_sclk_partitions(self.raw.as_ptr()) }
    }

    /// Encoded ticks to the instant on the clock's scale.
    pub fn to_time(&self, ticks: f64) -> Result<Time> {
        let mut t = Time::default();
        // SAFETY: the clock is alive
        check(unsafe { sys::myb_sclk_to_time(self.raw.as_ptr(), ticks, &mut t) })?;
        Ok(t)
    }

    pub fn from_time(&self, t: Time) -> Result<f64> {
        let mut ticks = 0.0;
        // SAFETY: the clock is alive
        check(unsafe { sys::myb_sclk_from_time(self.raw.as_ptr(), t, &mut ticks) })?;
        Ok(ticks)
    }

    /// A clock string ("partition/field field ...") to encoded ticks.
    pub fn encode(&self, text: &str) -> Result<f64> {
        let mut ticks = 0.0;
        // SAFETY: the text is text.len() bytes
        check(unsafe { sys::myb_sclk_encode(self.raw.as_ptr(), text.as_ptr().cast(), text.len(), &mut ticks) })?;
        Ok(ticks)
    }

    /// Encoded ticks to the clock string with its partition.
    pub fn decode(&self, ticks: f64) -> Result<String> {
        // SAFETY: the buffers are as large as said
        written(|buffer, size, length| unsafe { sys::myb_sclk_decode(self.raw.as_ptr(), ticks, buffer, size, length) })
    }

    pub fn as_ptr(&self) -> *const sys::myb_sclk {
        self.raw.as_ptr()
    }
}

impl Drop for Sclk {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_sclk_destroy(self.raw.as_ptr()) }
    }
}

// ---- CK files -------------------------------------------------------------------------------------

/// A segment of a CK set.
pub type CkSegmentInfo = sys::myb_ck_segment_info;
/// The orientation of an instrument at an instant of its clock.
pub type CkPointing = sys::myb_ck_pointing;

/// A set of CK files (`myb_ck`).
pub struct Ck {
    raw: NonNull<sys::myb_ck>,
    kept: Kept,
}

impl Ck {
    pub fn new() -> Result<Self> {
        let mut raw = ptr::null_mut();
        // SAFETY: a valid place for the pointer
        check(unsafe { sys::myb_ck_create(ptr::null(), &mut raw) })?;
        Ok(Ck { raw: created(raw)?, kept: Kept::default() })
    }

    pub fn load_file(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let path = c_path(path.as_ref())?;
        let mut diag = sys::myb_diag::default();
        // SAFETY: the path is NUL-terminated
        let code = unsafe { sys::myb_ck_load_file(self.raw.as_ptr(), path.as_ptr(), &mut diag) };
        check_diag(code, &diag)
    }

    /// A CK image in memory (what [`CkWriter::write`] gives); the set keeps the bytes.
    pub fn load_memory(&mut self, image: Vec<u8>) -> Result<()> {
        let image = image.into_boxed_slice();
        let (data, size) = (image.as_ptr(), image.len());
        self.kept.keep(image); // the bytes stay where they are
        let mut diag = sys::myb_diag::default();
        // SAFETY: the image lives in self.kept, unchanged, as long as the set
        let code = unsafe { sys::myb_ck_load_memory(self.raw.as_ptr(), data.cast(), size, &mut diag) };
        check_diag(code, &diag)
    }

    pub fn segments(&self) -> Vec<CkSegmentInfo> {
        // SAFETY: the set is alive; the index is within the count
        unsafe {
            (0..sys::myb_ck_segment_count(self.raw.as_ptr()))
                .map(|index| {
                    let mut info: CkSegmentInfo = std::mem::zeroed();
                    sys::myb_ck_segment(self.raw.as_ptr(), index, &mut info);
                    info
                })
                .collect()
        }
    }

    /// The orientation of an instrument at `ticks`, or the nearest within the tolerance.
    pub fn pointing_at(&self, instrument: i32, ticks: f64, tolerance: f64, need_rate: bool) -> Result<CkPointing> {
        let mut pointing = CkPointing::default();
        // SAFETY: the set is alive
        check(unsafe { sys::myb_ck_pointing_at(self.raw.as_ptr(), instrument, ticks, tolerance, c_int::from(need_rate), &mut pointing) })?;
        Ok(pointing)
    }

    pub fn as_ptr(&self) -> *const sys::myb_ck {
        self.raw.as_ptr()
    }
}

impl Drop for Ck {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once; the images after it
        unsafe { sys::myb_ck_destroy(self.raw.as_ptr()) }
    }
}

/// What a CK segment describes: an instrument in a base frame, and its name.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CkSegmentHeader {
    pub instrument: i32,
    pub frame: i32,
    pub name: String,
}

/// Packets for the sliding windows of types 5 and 6 (`myb_ck_packets`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CkPackets {
    /// 0 .. 3: packets of 8, 4, 14 or 7 numbers.
    pub subtype: u32,
    pub window: u32,
    pub seconds_per_tick: f64,
    pub ticks: Vec<f64>,
    pub packets: Vec<f64>,
}

impl CkPackets {
    fn raw(&self) -> Result<sys::myb_ck_packets> {
        let size = match self.subtype {
            0 => 8,
            1 => 4,
            2 => 14,
            3 => 7,
            _ => return domain(),
        };
        if self.packets.len() != size * self.ticks.len() {
            return domain();
        }
        Ok(sys::myb_ck_packets {
            subtype: self.subtype as c_int,
            window: c_int::try_from(self.window).or_else(|_| domain())?,
            seconds_per_tick: self.seconds_per_tick,
            ticks: self.ticks.as_ptr(),
            packets: self.packets.as_ptr(),
            count: self.ticks.len(),
        })
    }
}

fn ck_header(header: &CkSegmentHeader) -> Result<(sys::myb_ck_segment_header, CString)> {
    let name = c_string(&header.name)?;
    Ok((sys::myb_ck_segment_header { instrument: header.instrument, frame: header.frame, name: name.as_ptr() }, name))
}

fn optional(values: Option<&[f64]>) -> (*const f64, usize) {
    values.map_or((ptr::null(), 0), |v| (v.as_ptr(), v.len()))
}

/// Collects CK segments in memory and lays them out as one file (`myb_ck_writer`).
pub struct CkWriter {
    raw: NonNull<sys::myb_ck_writer>,
}

impl CkWriter {
    pub fn new(internal_name: Option<&str>) -> Result<Self> {
        let name = internal_name.map(c_string).transpose()?;
        let mut raw = ptr::null_mut();
        // SAFETY: the name, if any, is NUL-terminated
        check(unsafe { sys::myb_ck_writer_create(name.as_ref().map_or(ptr::null(), |n| n.as_ptr()), ptr::null(), &mut raw) })?;
        Ok(CkWriter { raw: created(raw)? })
    }

    /// Pointings at increasing ticks: type 1, or type 3 interpolated inside the intervals that
    /// begin at `interval_starts` (`None`: one interval).
    pub fn add_pointings(
        &mut self,
        header: &CkSegmentHeader,
        kind: i32,
        ticks: &[f64],
        quaternions: &[[f64; 4]],
        rates: Option<&[Vec3]>,
        interval_starts: Option<&[f64]>,
    ) -> Result<()> {
        if quaternions.len() != ticks.len() || rates.is_some_and(|r| r.len() != ticks.len()) {
            return domain();
        }
        let (made, _name) = ck_header(header)?;
        let (starts, start_count) = optional(interval_starts);
        // SAFETY: the arrays hold one entry per tick
        check(unsafe {
            sys::myb_ck_writer_add_pointings(
                self.raw.as_ptr(),
                &made,
                kind,
                ticks.as_ptr(),
                quaternions.as_ptr().cast(),
                rates.map_or(ptr::null(), |r| r.as_ptr().cast()),
                ticks.len(),
                starts,
                start_count,
            )
        })
    }

    /// Type 2: intervals [starts, stops] each with the quaternion at its start, a constant
    /// angular velocity and the seconds of a tick.
    pub fn add_intervals(
        &mut self,
        header: &CkSegmentHeader,
        starts: &[f64],
        stops: &[f64],
        quaternions: &[[f64; 4]],
        rates: &[Vec3],
        seconds_per_tick: &[f64],
    ) -> Result<()> {
        let n = starts.len();
        if stops.len() != n || quaternions.len() != n || rates.len() != n || seconds_per_tick.len() != n {
            return domain();
        }
        let (made, _name) = ck_header(header)?;
        // SAFETY: the arrays hold one entry per interval
        check(unsafe {
            sys::myb_ck_writer_add_intervals(
                self.raw.as_ptr(),
                &made,
                starts.as_ptr(),
                stops.as_ptr(),
                quaternions.as_ptr().cast(),
                rates.as_ptr().cast(),
                seconds_per_tick.as_ptr(),
                n,
            )
        })
    }

    /// Type 4: Chebyshev series on intervals; counts 7 per interval, the coefficients in a row.
    pub fn add_chebyshev(
        &mut self,
        header: &CkSegmentHeader,
        starts: &[f64],
        stops: &[f64],
        counts: &[i32],
        coefficients: &[f64],
        has_rates: bool,
    ) -> Result<()> {
        let n = starts.len();
        let total: i64 = counts.iter().map(|&c| i64::from(c.max(0))).sum();
        if stops.len() != n || counts.len() != 7 * n || coefficients.len() as i64 != total {
            return domain();
        }
        let (made, _name) = ck_header(header)?;
        // SAFETY: the arrays have the sizes checked above
        check(unsafe {
            sys::myb_ck_writer_add_chebyshev(
                self.raw.as_ptr(),
                &made,
                starts.as_ptr(),
                stops.as_ptr(),
                counts.as_ptr(),
                coefficients.as_ptr(),
                n,
                c_int::from(has_rates),
            )
        })
    }

    /// Type 5: packets covering [start, stop], interpolated inside the intervals.
    pub fn add_packets(
        &mut self,
        header: &CkSegmentHeader,
        packets: &CkPackets,
        has_rates: bool,
        start: f64,
        stop: f64,
        interval_starts: Option<&[f64]>,
    ) -> Result<()> {
        let raw = packets.raw()?;
        let (made, _name) = ck_header(header)?;
        let (starts, start_count) = optional(interval_starts);
        // SAFETY: the packets live in `packets` for the call
        check(unsafe {
            sys::myb_ck_writer_add_packets(self.raw.as_ptr(), &made, &raw, c_int::from(has_rates), start, stop, starts, start_count)
        })
    }

    /// Type 6: mini-segments between boundaries (one more than the mini-segments).
    #[allow(clippy::too_many_arguments)]
    pub fn add_minisegments(
        &mut self,
        header: &CkSegmentHeader,
        minisegments: &[CkPackets],
        boundaries: &[f64],
        has_rates: bool,
        select_last: bool,
        start: f64,
        stop: f64,
    ) -> Result<()> {
        if boundaries.len() != minisegments.len() + 1 {
            return domain();
        }
        let raws = minisegments.iter().map(CkPackets::raw).collect::<Result<Vec<_>>>()?;
        let (made, _name) = ck_header(header)?;
        // SAFETY: the packets live in `minisegments` for the call
        check(unsafe {
            sys::myb_ck_writer_add_minisegments(
                self.raw.as_ptr(),
                &made,
                raws.as_ptr(),
                boundaries.as_ptr(),
                raws.len(),
                c_int::from(has_rates),
                c_int::from(select_last),
                start,
                stop,
            )
        })
    }

    pub fn size(&self) -> usize {
        // SAFETY: the writer is alive
        unsafe { sys::myb_ck_writer_size(self.raw.as_ptr()) }
    }

    pub fn write(&self) -> Result<Vec<u8>> {
        let mut image = vec![0u8; self.size()];
        // SAFETY: image holds `size` bytes
        check(unsafe { sys::myb_ck_writer_write(self.raw.as_ptr(), image.as_mut_ptr().cast(), image.len()) })?;
        Ok(image)
    }

    pub fn write_file(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = c_path(path.as_ref())?;
        // SAFETY: the path is NUL-terminated
        check(unsafe { sys::myb_ck_writer_write_file(self.raw.as_ptr(), path.as_ptr()) })
    }
}

impl Drop for CkWriter {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_ck_writer_destroy(self.raw.as_ptr()) }
    }
}
