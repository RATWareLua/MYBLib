# Orbit estimation

**English** | [Русский](estimation.ru.md) · [All passports](index.en.md)

## What it is

Estimation of an orbit from observations: the state at an epoch, the parameters of the force
model and systematic shifts of the measurements, by weighted least squares over an arc
(Gauss–Newton on the variational equations) or by an extended Kalman filter with process noise,
followed if wanted by a Rauch–Tung–Striebel smoother. A measurement is any function of the state
of the object at the time of the observation with its derivatives; ready-made ones are the
position and the state of the object, range, range rate, right ascension and declination,
azimuth and elevation from a station on a body of the environment, and the measurement of a
point offset from the centre of mass through the attitude. The derivatives with respect to the
unknowns come from the transition matrix and the sensitivities of the
[propagation](propagation-uncertainty.en.md).

## Source

- Weighted least squares with normal equations scaled by their diagonal and solved by Cholesky;
  the extended Kalman filter with the covariance update in Joseph's form; process noise as an
  unmodelled white acceleration; the smoother of H. E. Rauch, F. Tung and C. T. Striebel —
  standard methods, implemented by the library (the header names no other reference).
- The attitude laws of `myb/surface.h` are differentiated as they are written, a law of the
  caller's by central differences of the law.
- No third-party code, no data files.

## Interface

Header [`myb/estim.h`](../../include/myb/estim.h):

- measurements (type `myb_measure_function`): `myb_measure_position`, `myb_measure_state`,
  `myb_measure_range`, `myb_measure_range_rate`, `myb_measure_radec`, `myb_measure_azel`,
  `myb_measure_offset`; with the light time `myb_measure_light_range`,
  `myb_measure_light_range_rate`, `myb_measure_light_radec`, `myb_measure_light_azel` (see
  [observables and light time](observables-light-time.en.md));
- structures `myb_observation`, `myb_station`, `myb_station_motion`, `myb_offset_measure`,
  `myb_light_measure`; `myb_station_init`, `myb_station_state` (the motion of a station:
  [site displacement](site-displacement.en.md));
- least squares: `myb_fit_orbit`, `myb_fit_options`, `myb_fit_result`;
- the filter: `myb_filter` (opaque), `myb_filter_options`, `myb_filter_create`,
  `myb_filter_destroy`, `myb_filter_start`, `myb_filter_update`, `myb_filter_predict`,
  `myb_filter_get`, `myb_filter_set_bias`, `myb_filter_bias`;
- the smoother: `myb_filter_steps`, `myb_filter_smooth`, `myb_filter_step`;
- constants `MYB_MEASURE_MAX` (6), `MYB_FIT_MAX_BIASES` (32).

## Inputs and outputs

- The unknowns, in order: y0 = (x, v) at the epoch t0, the np parameters of the force model
  (`myb_force_model_add_parameter`), the biases — 6 + np + nb numbers. Times are the model time
  in s; positions in m and velocities in m/s, relative to the origin of the environment in its
  inertial axes.
- An observation: its time, the measurement and its context, 1 … 6 observed values with their
  standard deviations (> 0); a bit per value marks an angle whose residual is reduced to
  (−π, π]; `bias = k > 0` adds the shift k − 1 + i of the estimator to value i (the delay of a
  station, the offset of a clock or of a pointing). Residuals are observed − computed.
- A station: geodetic longitude and latitude (rad) and height (m) on an ellipsoid of the body,
  the axes east, north, up; fixed on the body, or with a `myb_station_motion` (tides, pole tide,
  loads) added at every instant.
- Range (m) and range rate (m/s) of the vector from the station to the object; right ascension
  in [0, 2π) and declination in the inertial axes; azimuth from north towards east in [0, 2π)
  and elevation above the horizon of the station. These are geometric: the state at the instant
  of the observation.
- An offset point: the offset in m in the body frame, turned by the attitude; its velocity takes
  the turning of the body from the attitudes half a second before and after — at the states of
  a Kepler step about a point mass `gm` at the origin when gm > 0 (a law that follows the
  velocity turns with it), in a straight line when gm = 0.
- `myb_fit_orbit`: observations in time order away from t0, all on one side of it (an arc
  forward or backward). It stops when the weighted RMS changes by less than `tolerance`
  (10⁻⁶) of itself or after `max_iterations` (10) corrections; IAS15 precision 10⁻¹⁰; unknowns
  may be held fixed or given a priori sigmas. Out: the estimate (y0 and the parameters in
  place), the covariance (the inverse of the normal matrix; zero rows for fixed unknowns), the
  residuals, the weighted RMS before and after (1 when the errors match the sigmas).
- The filter: between observations the state goes by the model and the covariance by
  P ← F P Fᵀ + Q, F = [Φ S; 0 1], Q = q [Δt³/3 Δt²/2; Δt²/2 Δt] per axis in the block of y
  (q in m² s⁻³); each observation corrects both with K = P Hᵀ (H P Hᵀ + R)⁻¹, the covariance in
  Joseph's form. The parameters of the model are corrected in place; biases start at zero or at
  `myb_filter_set_bias`.
- The smoother goes back over the instants the filter kept (`history`) and gives each the
  estimate from all the observations with its covariance; at the first instant that is the orbit
  of a fit over the whole arc.

## Domain of validity

- Any force model of the library; at most 16 parameters (`MYB_MODEL_MAX_PARAMETERS`), 32 biases,
  6 values per observation.
- The derivatives of an offset point through the attitude are a part in |offset| / |x| of the
  whole.
- Refused with `MYB_E_DOMAIN`: no observations, an invalid one, times out of order, a normal
  matrix that is not positive definite (the observations do not determine the unknowns); in the
  filter negative options, a time behind the filter's or beyond t_end, an innovation covariance
  that is not positive definite; a smoother without history; the object at the station; angles
  exactly at the pole of their axes; a negative `gm`.
- `MYB_E_OVERFLOW`: an observation would need more instants than the filter keeps;
  `MYB_E_NOMEM`; the statuses of the propagation and of the measurements. On an error of the fit
  y0 and the parameters hold the values of the last iteration.

## Accuracy

`test_estim` ([tests/unit/test_estim.c](../../tests/unit/test_estim.c)) works on a made-up
planet (a field of degree 4, uniform rotation, a layer of exponential atmosphere, three
stations, an orbit of 420 km), so nothing is downloaded:

| What | Reference | Result |
|---|---|---|
| Derivatives of every measurement | central differences | up to 3·10⁻¹⁰ |
| Least squares without noise, the guess displaced by hundreds of metres, c_d off by 0.5 | the truth | 2·10⁻⁹ m, 6·10⁻¹¹ in c_d |
| Least squares with noise (σ 1 m, 1 mm/s, 2 µrad) | the truth and its own covariance | weighted RMS 1.035, largest error 1.7σ |
| The filter, the same noisy observations from a covariance of a kilometre | the truth | 0.68 m at the end with σ 0.43 m, largest error 1.9σ |
| Process noise over 600 s | q Δt³/3 | 2.166 m² against 2.16 |
| Three range delays and a range-rate offset (35, −12, 4.5 m; 2.0 mm/s put in) | the shifts | σ 0.19 m, largest error 1.1σ, weighted RMS 0.963 |
| The smoother, the first of 140 instants | a fit over the same arc | estimates 0.000σ apart, sigmas to 7·10⁻⁶; position σ 520 m before the pass, 0.94 m after |
| An offset point, orbital axes | the point in closed form | position exact, velocity 3·10⁻¹⁰ m/s |

`test_light_measure` ([tests/unit/test_light_measure.c](../../tests/unit/test_light_measure.c)):
an arm of 37 m with an attitude along the velocity — the velocity of the point against a
difference of the true orbit 5.7·10⁻⁹ m/s with `gm` (3.4·10⁻² m/s without); the attitude adds up
to 2·10⁻³ to the derivatives, which agree with differences within what a difference allows.

Real orbits through the same estimator:

| What | Reference | Result |
|---|---|---|
| LAGEOS-1, 2025-06-01 … 07, 1680 positions; the state and five empirical accelerations fitted, c_r = 1.13 held (`test_lageos`) | ILRS combined orbit, SLRF2020 | RMS 1.57 cm (R 2.8 mm, T 1.3 cm, N 8 mm); 2.05 cm before the secondary waves of the ocean tide were added |
| LAGEOS-1, residual of the first guess | itself at ε = 10⁻¹⁰, 10⁻¹¹, 10⁻¹² | 25.8887 m at all three (the edges of shadows declared as breaks) |
| Tesla Roadster, 2018-03 … 2025-03 (`test_spacecraft`) | JPL Horizons | RMS 306 m, largest 1.3 km |
| OSIRIS-REx through the Earth flyby of 2017-09-22, 6 days | JPL Horizons | RMS 12.6 m |
| JWST, the week 2023-02-01 … 08 | JPL Horizons | RMS 260 m |

([tests/unit/test_lageos.c](../../tests/unit/test_lageos.c),
[tests/unit/test_spacecraft.c](../../tests/unit/test_spacecraft.c).)

## Cost

- LAGEOS: the fit converges in two corrections; measured in phase 14 at 21 s (against ~35 s for
  the former fit by differences of whole propagations). Ending the steps at the edges of shadows
  (phase 16) costs +26 % of the time of that fit.

## Limitations

- The station measurements of the list above are geometric: no light time and no aberration —
  the `myb_measure_light_*` functions have them. Delays of the media (troposphere, ionosphere)
  and the relativistic delay are not part of any measurement: they are corrections the caller
  adds to the observed value.
- The fit needs the observations ordered in time on one side of the epoch.
- The filter linearises about its current estimate; process noise is a white acceleration of one
  density on every axis of the state block, none on the parameters or the biases.
- The smoother sees only the instants kept: `history` sets how many.
- The velocity of an offset point with gm = 0 misses the turning of an attitude that follows the
  velocity alone (3.4·10⁻² m/s on an arm of 37 m).
- Partials of a law of the caller's are central differences of the law.

## Data

None for the estimator. The tests on real orbits read, from the repository,
`tests/data/sp3/ilrsa.lageos1.250607.excerpt.sp3` (ILRS, free with acknowledgement) and
`tests/data/horizons/*.txt` (JPL Horizons), and the downloaded files (`cmake/FetchData.cmake`)
`de440.bsp`, `header.440`, `eopc04.1962-now`, `GOCO06s.gfc` (CC BY 4.0) and
`fes2004_Cnm-Snm.dat`; they are skipped without them.
