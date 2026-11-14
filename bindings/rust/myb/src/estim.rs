// SPDX-License-Identifier: BSD-3-Clause
//! Estimation of an orbit (`myb/estim.h`): stations on the bodies of an environment, the
//! measurements of the library and of the caller (values with their derivatives), the fit of
//! an orbit by least squares over an arc, and the sequential filter with its smoother.
//!
//! The unknowns are the state at the epoch, the parameters of the force model and the
//! systematic shifts of the measurements, in that order.
use std::any::Any;
use std::ffi::{c_int, c_void};
use std::ptr::{self, NonNull};
use std::rc::Rc;

use crate::{
    c_path, c_string, check, check_diag, created, domain, out_ptr, sys, AttitudeSource, Environment, ForceModel, KernelPool,
    Result, State, Vec3,
};

/// An ellipsoid of a body: equatorial radius a (m) and flattening f.
pub type Ellipsoid = sys::myb_ellipsoid;
/// What a fit reached: iterations, convergence, the weighted RMS before and after.
pub type FitResult = sys::myb_fit_result;

// ---- What moves a station -------------------------------------------------------------------------

/// The tides of a body for the displacements of its stations (`myb/site.h`), from the kernel
/// variables under a prefix.
pub struct SiteTides {
    raw: NonNull<sys::myb_site_tides>,
}

impl SiteTides {
    pub fn new(pool: &KernelPool, prefix: &str) -> Result<Self> {
        let prefix = c_string(prefix)?;
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the pool is alive for the call; valid places for the outputs
        let code = unsafe { sys::myb_site_tides_create(pool.as_ptr(), prefix.as_ptr(), ptr::null(), &mut raw, &mut diag) };
        check_diag(code, &diag)?;
        Ok(SiteTides { raw: created(raw)? })
    }

    /// The waves of the potential.
    pub fn waves(&self) -> usize {
        // SAFETY: the tides are alive
        unsafe { sys::myb_site_tides_waves(self.raw.as_ptr()) }
    }

    /// The waves a table of ocean loading gives.
    pub fn given_waves(&self) -> usize {
        // SAFETY: the tides are alive
        unsafe { sys::myb_site_tides_given_waves(self.raw.as_ptr()) }
    }

    pub fn as_ptr(&self) -> *const sys::myb_site_tides {
        self.raw.as_ptr()
    }
}

impl Drop for SiteTides {
    fn drop(&mut self) {
        // SAFETY: created by myb_site_tides_create, dropped once
        unsafe { sys::myb_site_tides_destroy(self.raw.as_ptr()) }
    }
}

/// The load of the ocean tides at a point (`myb/site.h`).
pub struct OceanLoading {
    raw: NonNull<sys::myb_ocean_loading>,
    _tides: Rc<SiteTides>,
}

impl OceanLoading {
    /// Amplitudes (m) and phase lags (degrees): three rows (up, west, south) of the given waves.
    pub fn new(tides: &Rc<SiteTides>, amplitudes: &[f64], phases: &[f64]) -> Result<Self> {
        let needed = 3 * tides.given_waves();
        if amplitudes.len() != needed || phases.len() != needed {
            return domain();
        }
        let mut raw = ptr::null_mut();
        // SAFETY: the arrays hold three rows of the given waves
        check(unsafe { sys::myb_ocean_loading_create(tides.as_ptr(), amplitudes.as_ptr(), phases.as_ptr(), ptr::null(), &mut raw) })?;
        Ok(OceanLoading { raw: created(raw)?, _tides: Rc::clone(tides) })
    }

    /// The block of a site in the text of a BLQ file.
    pub fn load(tides: &Rc<SiteTides>, text: &str, site: &str) -> Result<Self> {
        let site = c_string(site)?;
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: the text is text.len() bytes; valid places for the outputs
        let code = unsafe {
            sys::myb_ocean_loading_load(tides.as_ptr(), text.as_ptr().cast(), text.len(), site.as_ptr(), ptr::null(), &mut raw, &mut diag)
        };
        check_diag(code, &diag)?;
        Ok(OceanLoading { raw: created(raw)?, _tides: Rc::clone(tides) })
    }

    /// The block of a site in a BLQ file.
    pub fn load_file(tides: &Rc<SiteTides>, path: impl AsRef<std::path::Path>, site: &str) -> Result<Self> {
        let (path, site) = (c_path(path.as_ref())?, c_string(site)?);
        let mut raw = ptr::null_mut();
        let mut diag = sys::myb_diag::default();
        // SAFETY: NUL-terminated texts; valid places for the outputs
        let code = unsafe {
            sys::myb_ocean_loading_load_file(tides.as_ptr(), path.as_ptr(), site.as_ptr(), ptr::null(), &mut raw, &mut diag)
        };
        check_diag(code, &diag)?;
        Ok(OceanLoading { raw: created(raw)?, _tides: Rc::clone(tides) })
    }

    pub fn as_ptr(&self) -> *const sys::myb_ocean_loading {
        self.raw.as_ptr()
    }
}

impl Drop for OceanLoading {
    fn drop(&mut self) {
        // SAFETY: created by the library, dropped once
        unsafe { sys::myb_ocean_loading_destroy(self.raw.as_ptr()) }
    }
}

/// The loads of a point: the ocean tides, the atmospheric tides S1 and S2, the ocean pole tide.
#[derive(Clone, Default)]
pub struct SiteLoads {
    pub ocean: Option<Rc<OceanLoading>>,
    pub atmosphere: [f64; 12],
    pub ocean_pole: [f64; 6],
}

/// What moves a station on its body (`myb_station_motion`).
#[derive(Clone)]
pub struct StationMotion {
    pub tides: Rc<SiteTides>,
    pub loads: Option<SiteLoads>,
    /// Bodies of the environment that raise the tides (at most 8).
    pub raisers: Vec<usize>,
    /// MYB_SITE_* of `myb/site.h`.
    pub flags: i32,
    /// The velocity of the displacement too.
    pub rate: bool,
}

struct Motion {
    raw: sys::myb_station_motion,
    _loads: Option<Box<sys::myb_site_loads>>,
    _tides: Rc<SiteTides>,
    _ocean: Option<Rc<OceanLoading>>,
}

// ---- Stations -------------------------------------------------------------------------------------

/// An observer on a body of an environment (`myb_station`).
pub struct Station {
    raw: sys::myb_station,
    env: Rc<Environment>,
    _motion: Option<Box<Motion>>,
}

impl Station {
    /// At geodetic longitude, latitude (rad) and height (m) on an ellipsoid of the body.
    pub fn new(env: &Rc<Environment>, body: usize, ellipsoid: &Ellipsoid, longitude: f64, latitude: f64, height: f64) -> Result<Self> {
        if body >= env.body_count() {
            return domain();
        }
        let mut raw = sys::myb_station::default();
        // SAFETY: the environment is kept alive by the station
        check(unsafe { sys::myb_station_init(&mut raw, env.as_ptr(), body, ellipsoid, longitude, latitude, height) })?;
        Ok(Station { raw, env: Rc::clone(env), _motion: None })
    }

    /// The same station whose position is that of a terrestrial frame, moved at every instant
    /// by the tides and loads of the motion.
    pub fn with_motion(mut self, motion: &StationMotion) -> Result<Self> {
        let mut raw = sys::myb_station_motion::default();
        if motion.raisers.len() > raw.raisers.len() || motion.raisers.iter().any(|&b| b >= self.env.body_count()) {
            return domain();
        }
        raw.raisers[..motion.raisers.len()].copy_from_slice(&motion.raisers);
        raw.raiser_count = motion.raisers.len();
        raw.flags = motion.flags;
        raw.rate = c_int::from(motion.rate);
        raw.tides = motion.tides.as_ptr();
        let loads = motion.loads.as_ref().map(|loads| {
            Box::new(sys::myb_site_loads {
                ocean: loads.ocean.as_ref().map_or(ptr::null(), |o| o.as_ptr()),
                atmosphere: loads.atmosphere,
                ocean_pole: loads.ocean_pole,
            })
        });
        raw.loads = loads.as_ref().map_or(ptr::null(), |l| &**l as *const sys::myb_site_loads);
        let kept = Box::new(Motion {
            raw,
            _loads: loads,
            _tides: Rc::clone(&motion.tides),
            _ocean: motion.loads.as_ref().and_then(|l| l.ocean.clone()),
        });
        self.raw.motion = &kept.raw;
        self._motion = Some(kept);
        Ok(self)
    }

    /// The body-fixed position (m).
    pub fn position(&self) -> Vec3 {
        self.raw.position
    }

    /// The axes east, north, up (rows, body-fixed frame).
    pub fn enu(&self) -> [f64; 9] {
        self.raw.enu
    }

    pub fn body(&self) -> usize {
        self.raw.body
    }

    /// The state relative to the origin of the environment, inertial axes, at the model time t.
    pub fn state(&self, t: f64) -> Result<State> {
        let mut y = [0.0; 6];
        // SAFETY: the station, its environment and its motion are alive
        check(unsafe { sys::myb_station_state(&self.raw, t, y.as_mut_ptr()) })?;
        Ok(y)
    }

    pub fn as_ptr(&self) -> *const sys::myb_station {
        &self.raw
    }
}

/// The light-time measurements from a station (`myb_light_measure`).
#[derive(Clone)]
pub struct LightMeasure {
    pub station: Rc<Station>,
    /// The speed of light, m/s.
    pub c: f64,
    /// Of the body at the origin (m^3 s^-2); 0: straight motion over the light time.
    pub gm: f64,
    /// The station emits and the object returns.
    pub two_way: bool,
    /// The angles as the moving station sees them.
    pub aberration: bool,
}

// ---- Measurements ---------------------------------------------------------------------------------

type Measure = dyn Fn(f64, &State, &mut [f64], Option<&mut [f64]>) -> Result<()>;

struct Custom {
    count: usize,
    measure: Box<Measure>,
}

unsafe extern "C" fn call_measure(context: *mut c_void, t: f64, y: *const f64, values: *mut f64, partials: *mut f64) -> sys::myb_status {
    // SAFETY: the context is the Custom a Measurement keeps; y holds a state, values count
    // numbers and partials, if any, 6 count
    let custom = unsafe { &*(context as *const Custom) };
    let state = unsafe { *(y as *const State) };
    let values = unsafe { std::slice::from_raw_parts_mut(values, custom.count) };
    let partials = (!partials.is_null()).then(|| unsafe { std::slice::from_raw_parts_mut(partials, 6 * custom.count) });
    match (custom.measure)(t, &state, values, partials) {
        Ok(()) => sys::MYB_OK,
        Err(error) => error.status.raw(),
    }
}

struct Offset {
    raw: sys::myb_offset_measure,
    _measure: Measurement,
    _attitude: Box<dyn AttitudeSource>,
}

/// A measurement: `count` values computed from the state of the object at a model time, with
/// their derivatives by the state — one of the library's or a closure.
#[derive(Clone)]
pub struct Measurement {
    function: sys::myb_measure_function,
    context: *mut c_void,
    count: usize,
    _kept: Option<Rc<dyn Any>>,
}

impl Measurement {
    fn of_station(function: sys::myb_measure_function, station: &Rc<Station>, count: usize) -> Self {
        let context = station.as_ptr() as *mut c_void;
        Measurement { function, context, count, _kept: Some(Rc::clone(station) as Rc<dyn Any>) }
    }

    fn of_light(function: sys::myb_measure_function, light: &LightMeasure, count: usize) -> Self {
        let kept = Rc::new((
            sys::myb_light_measure {
                station: light.station.as_ptr(),
                c: light.c,
                gm: light.gm,
                two_way: c_int::from(light.two_way),
                aberration: c_int::from(light.aberration),
            },
            Rc::clone(&light.station),
        ));
        let context = &kept.0 as *const sys::myb_light_measure as *mut c_void;
        Measurement { function, context, count, _kept: Some(kept as Rc<dyn Any>) }
    }

    /// The position of the object itself (3 values).
    pub fn position() -> Self {
        Measurement { function: Some(sys::myb_measure_position), context: ptr::null_mut(), count: 3, _kept: None }
    }

    /// The whole state of the object (6 values).
    pub fn state() -> Self {
        Measurement { function: Some(sys::myb_measure_state), context: ptr::null_mut(), count: 6, _kept: None }
    }

    /// The range from a station (1 value, m).
    pub fn range(station: &Rc<Station>) -> Self {
        Self::of_station(Some(sys::myb_measure_range), station, 1)
    }

    /// The range rate (1 value, m/s).
    pub fn range_rate(station: &Rc<Station>) -> Self {
        Self::of_station(Some(sys::myb_measure_range_rate), station, 1)
    }

    /// Right ascension and declination in the inertial axes (2 values; the first periodic).
    pub fn radec(station: &Rc<Station>) -> Self {
        Self::of_station(Some(sys::myb_measure_radec), station, 2)
    }

    /// Azimuth and elevation at the station (2 values; the first periodic).
    pub fn azel(station: &Rc<Station>) -> Self {
        Self::of_station(Some(sys::myb_measure_azel), station, 2)
    }

    /// The range with the light time.
    pub fn light_range(light: &LightMeasure) -> Self {
        Self::of_light(Some(sys::myb_measure_light_range), light, 1)
    }

    /// The range rate with the light time.
    pub fn light_range_rate(light: &LightMeasure) -> Self {
        Self::of_light(Some(sys::myb_measure_light_range_rate), light, 1)
    }

    /// The direction the light comes from, as right ascension and declination.
    pub fn light_radec(light: &LightMeasure) -> Self {
        Self::of_light(Some(sys::myb_measure_light_radec), light, 2)
    }

    /// The direction the light comes from, as azimuth and elevation.
    pub fn light_azel(light: &LightMeasure) -> Self {
        Self::of_light(Some(sys::myb_measure_light_azel), light, 2)
    }

    /// `measure` taken at a point offset (m, body frame) from the centre of mass and turned by
    /// the attitude of the object (`myb_measure_offset`); gm of the body at the origin, or 0.
    pub fn offset(measure: &Measurement, attitude: impl AttitudeSource + 'static, offset: Vec3, gm: f64) -> Self {
        let attitude: Box<dyn AttitudeSource> = Box::new(attitude);
        let (function, context) = attitude.as_raw();
        let kept = Rc::new(Offset {
            raw: sys::myb_offset_measure {
                measure: measure.function,
                context: measure.context,
                attitude: function,
                attitude_context: context,
                offset,
                gm,
            },
            _measure: measure.clone(),
            _attitude: attitude,
        });
        let context = &kept.raw as *const sys::myb_offset_measure as *mut c_void;
        Measurement { function: Some(sys::myb_measure_offset), context, count: measure.count, _kept: Some(kept as Rc<dyn Any>) }
    }

    /// A measurement of the caller: the closure fills `count` values at (t, y) and, when asked,
    /// their derivatives partials[6 i + j] = d value_i / d y_j.
    pub fn custom(count: usize, measure: impl Fn(f64, &State, &mut [f64], Option<&mut [f64]>) -> Result<()> + 'static) -> Self {
        let kept = Rc::new(Custom { count, measure: Box::new(measure) });
        let context = Rc::as_ptr(&kept) as *mut c_void;
        Measurement { function: Some(call_measure), context, count, _kept: Some(kept as Rc<dyn Any>) }
    }

    /// The number of values.
    pub fn count(&self) -> usize {
        self.count
    }

    /// The values at the model time t and the state y.
    pub fn value(&self, t: f64, y: &State) -> Result<Vec<f64>> {
        let mut values = vec![0.0; self.count];
        let function = self.function.ok_or(crate::Error::new(crate::Status::Domain))?;
        // SAFETY: the context lives with self; values holds count numbers
        check(unsafe { function(self.context, t, y.as_ptr(), values.as_mut_ptr(), ptr::null_mut()) })?;
        Ok(values)
    }

    /// The values and their derivatives by the state (count x 6, rows).
    pub fn evaluate(&self, t: f64, y: &State) -> Result<(Vec<f64>, Vec<f64>)> {
        let mut values = vec![0.0; self.count];
        let mut partials = vec![0.0; 6 * self.count];
        let function = self.function.ok_or(crate::Error::new(crate::Status::Domain))?;
        // SAFETY: the context lives with self; the outputs hold count and 6 count numbers
        check(unsafe { function(self.context, t, y.as_ptr(), values.as_mut_ptr(), partials.as_mut_ptr()) })?;
        Ok((values, partials))
    }

    pub fn as_raw(&self) -> (sys::myb_measure_function, *mut c_void) {
        (self.function, self.context)
    }
}

/// An observation: what a measurement gave at a model time, with its standard deviations.
#[derive(Clone)]
pub struct Observation {
    /// Model time, s.
    pub t: f64,
    pub measurement: Measurement,
    /// Observed values, as many as the measurement gives.
    pub value: Vec<f64>,
    /// Standard deviations, > 0.
    pub sigma: Vec<f64>,
    /// Bit i: value i is an angle.
    pub periodic: u32,
    /// 0: none; k > 0: value i takes the systematic shift number k - 1 + i.
    pub bias: u32,
}

impl Observation {
    pub fn new(t: f64, measurement: &Measurement, value: &[f64], sigma: &[f64]) -> Self {
        Observation { t, measurement: measurement.clone(), value: value.to_vec(), sigma: sigma.to_vec(), periodic: 0, bias: 0 }
    }

    fn raw(&self) -> Result<sys::myb_observation> {
        let count = self.measurement.count;
        let mut made = sys::myb_observation {
            t: self.t,
            measure: self.measurement.function,
            context: self.measurement.context,
            count,
            periodic: self.periodic,
            bias: c_int::try_from(self.bias).or_else(|_| domain())?,
            ..Default::default()
        };
        if count > made.value.len() || self.value.len() != count || self.sigma.len() != count {
            return domain();
        }
        made.value[..count].copy_from_slice(&self.value);
        made.sigma[..count].copy_from_slice(&self.sigma);
        Ok(made)
    }
}

// ---- Least squares --------------------------------------------------------------------------------

/// How a fit goes; zeros stand for the defaults of the library.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FitOptions {
    pub max_iterations: u32,
    pub tolerance: f64,
    /// IAS15 precision of the propagation.
    pub epsilon: f64,
    /// For each unknown (6 + parameters + shifts): true keeps its value.
    pub fixed: Option<Vec<bool>>,
    /// For each unknown: the a priori standard deviation of its first value (0: none).
    pub prior_sigma: Option<Vec<f64>>,
    /// The first guesses of the systematic shifts; their number is that of the shifts.
    pub bias: Vec<f64>,
}

/// What a fit gives.
#[derive(Debug, Clone)]
pub struct Fit {
    /// The state at the epoch.
    pub y0: State,
    /// The estimated shifts.
    pub bias: Vec<f64>,
    pub result: FitResult,
    /// Of the unknowns, (6 + np + shifts)^2, rows.
    pub covariance: Vec<f64>,
    /// Observed - computed of every value, in the order of the observations.
    pub residuals: Vec<f64>,
}

/// Fits the state at t0 (from the first guess y0) and the parameters of the model (in place)
/// to the observations, ordered in time away from t0 (`myb_fit_orbit`).
pub fn fit_orbit(model: &ForceModel, t0: f64, y0: &State, observations: &[Observation], options: &FitOptions) -> Result<Fit> {
    let raws = observations.iter().map(Observation::raw).collect::<Result<Vec<_>>>()?;
    let unknowns = 6 + model.parameter_count() + options.bias.len();
    let fixed: Option<Vec<u8>> = options.fixed.as_ref().map(|f| f.iter().map(|&b| u8::from(b)).collect());
    if fixed.as_ref().is_some_and(|f| f.len() != unknowns) || options.prior_sigma.as_ref().is_some_and(|p| p.len() != unknowns) {
        return domain();
    }
    let mut bias = options.bias.clone();
    let made = sys::myb_fit_options {
        max_iterations: c_int::try_from(options.max_iterations).or_else(|_| domain())?,
        tolerance: options.tolerance,
        epsilon: options.epsilon,
        fixed: fixed.as_ref().map_or(ptr::null(), |f| f.as_ptr()),
        prior_sigma: options.prior_sigma.as_ref().map_or(ptr::null(), |p| p.as_ptr()),
        bias_count: bias.len(),
        bias: out_ptr(&mut bias),
    };
    let values: usize = raws.iter().map(|o| o.count).sum();
    let mut y = *y0;
    let mut result = FitResult::default();
    let mut covariance = vec![0.0; unknowns * unknowns];
    let mut residuals = vec![0.0; values];
    // SAFETY: the observations point to measurements alive in `observations`; the outputs hold
    // the sizes the library asks for
    check(unsafe {
        sys::myb_fit_orbit(
            model.as_ptr(),
            t0,
            y.as_mut_ptr(),
            raws.as_ptr(),
            raws.len(),
            &made,
            ptr::null(),
            &mut result,
            covariance.as_mut_ptr(),
            out_ptr(&mut residuals),
        )
    })?;
    Ok(Fit { y0: y, bias, result, covariance, residuals })
}

// ---- The filter -----------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct FilterOptions {
    /// IAS15 precision (0: 1e-10).
    pub epsilon: f64,
    /// The spectral density of an unmodelled white acceleration, m^2 s^-3.
    pub process_noise: f64,
    /// Systematic shifts of the measurements.
    pub bias_count: usize,
    /// Instants kept for the smoother (0: none).
    pub history: usize,
}

/// The time, the state and the covariance of the unknowns.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterEstimate {
    pub t: f64,
    pub y: State,
    pub covariance: Vec<f64>,
}

/// The estimate at an instant the filter kept: its unknowns (state, parameters, shifts) and
/// their covariance.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterStep {
    pub t: f64,
    pub unknowns: Vec<f64>,
    pub covariance: Vec<f64>,
}

/// The extended Kalman filter on the state and the parameters of a model (`myb_filter`).
pub struct Filter {
    raw: NonNull<sys::myb_filter>,
    _model: Rc<ForceModel>,
    parameters: usize,
    bias_count: usize,
}

impl Filter {
    pub fn new(model: &Rc<ForceModel>, options: &FilterOptions) -> Result<Self> {
        let made = sys::myb_filter_options {
            epsilon: options.epsilon,
            process_noise: options.process_noise,
            bias_count: options.bias_count,
            history: options.history,
        };
        let mut raw = ptr::null_mut();
        // SAFETY: the model outlives the filter (it keeps it)
        check(unsafe { sys::myb_filter_create(model.as_ptr(), &made, ptr::null(), &mut raw) })?;
        Ok(Filter { raw: created(raw)?, _model: Rc::clone(model), parameters: model.parameter_count(), bias_count: options.bias_count })
    }

    fn unknowns(&self) -> usize {
        6 + self.parameters + self.bias_count
    }

    /// Starts at t0 with the state, the present parameters and the covariance of the
    /// 6 + np unknowns (rows), towards t_end.
    pub fn start(&mut self, t0: f64, y0: &State, covariance: &[f64], t_end: f64) -> Result<()> {
        let n = 6 + self.parameters;
        if covariance.len() != n * n {
            return domain();
        }
        // SAFETY: the covariance holds (6 + np)^2 numbers
        check(unsafe { sys::myb_filter_start(self.raw.as_ptr(), t0, y0.as_ptr(), covariance.as_ptr(), t_end) })
    }

    /// Propagates to the observation and applies it; returns the innovation.
    pub fn update(&mut self, observation: &Observation) -> Result<Vec<f64>> {
        let made = observation.raw()?;
        let mut innovation = vec![0.0; made.count];
        // SAFETY: the measurement is alive in `observation`; innovation holds count numbers
        check(unsafe { sys::myb_filter_update(self.raw.as_ptr(), &made, innovation.as_mut_ptr()) })?;
        Ok(innovation)
    }

    /// Propagates to t without an observation.
    pub fn predict(&mut self, t: f64) -> Result<()> {
        // SAFETY: the filter is alive
        check(unsafe { sys::myb_filter_predict(self.raw.as_ptr(), t) })
    }

    /// The time, the state and the covariance ((6 + np + shifts)^2) of the filter.
    pub fn get(&self) -> FilterEstimate {
        let n = self.unknowns();
        let mut out = FilterEstimate { t: 0.0, y: [0.0; 6], covariance: vec![0.0; n * n] };
        // SAFETY: the outputs hold the sizes of the filter
        unsafe { sys::myb_filter_get(self.raw.as_ptr(), &mut out.t, out.y.as_mut_ptr(), out.covariance.as_mut_ptr()) };
        out
    }

    /// The first values of the systematic shifts (after start, before the first observation).
    pub fn set_bias(&mut self, bias: &[f64]) -> Result<()> {
        if bias.len() != self.bias_count {
            return domain();
        }
        // SAFETY: bias holds bias_count numbers
        check(unsafe { sys::myb_filter_set_bias(self.raw.as_ptr(), bias.as_ptr()) })
    }

    /// The present estimates of the shifts.
    pub fn bias(&self) -> Vec<f64> {
        let mut bias = vec![0.0; self.bias_count];
        if !bias.is_empty() {
            // SAFETY: bias holds bias_count numbers
            unsafe { sys::myb_filter_bias(self.raw.as_ptr(), bias.as_mut_ptr()) };
        }
        bias
    }

    /// The instants kept so far.
    pub fn steps(&self) -> usize {
        // SAFETY: the filter is alive
        unsafe { sys::myb_filter_steps(self.raw.as_ptr()) }
    }

    /// The backward pass of the smoother over the instants kept.
    pub fn smooth(&mut self) -> Result<()> {
        // SAFETY: the filter is alive
        check(unsafe { sys::myb_filter_smooth(self.raw.as_ptr()) })
    }

    /// The estimate at an instant kept (`None` if there is no such instant).
    pub fn step(&self, index: usize) -> Option<FilterStep> {
        if index >= self.steps() {
            return None;
        }
        let n = self.unknowns();
        let mut out = FilterStep { t: 0.0, unknowns: vec![0.0; n], covariance: vec![0.0; n * n] };
        // SAFETY: the instant is one kept; the outputs hold the sizes of the filter
        unsafe { sys::myb_filter_step(self.raw.as_ptr(), index, &mut out.t, out.unknowns.as_mut_ptr(), out.covariance.as_mut_ptr()) };
        Some(out)
    }
}

impl Drop for Filter {
    fn drop(&mut self) {
        // SAFETY: created by myb_filter_create, dropped once, before its model
        unsafe { sys::myb_filter_destroy(self.raw.as_ptr()) }
    }
}
