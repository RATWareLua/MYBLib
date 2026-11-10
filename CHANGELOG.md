# Changelog

**English** | [Русский](CHANGELOG.ru.md)

Versions are major.minor.patch.tweak (`MYB_VERSION_*` in `include/myb/core.h`). Within a major
version the API only grows. Details, checks and measured accuracies of every phase:
[`docs/PLAN.en.md`](docs/PLAN.en.md), §12.

## 1.2.0.0 — in progress (no tag)

### Added
- `myb/env.h`: the environment of a force model (bodies with providers of state and
  orientation, values once per instant, time scales of the providers from a table) and
  ready-made terms: point mass, third body, gravity field with tides, radiation pressure with
  shadows, radiation of a planet, drag, relativity, empirical accelerations.
- `myb/orbit.h`: propagation of one object — states at the times asked for (independent of the
  sampling, bit for bit), events with a direction filter, state transition matrix.
- Providers: `myb_spk_body` (`myb/kin.h`); `myb_orientation_function` with the CIP + EOP, PCK,
  WGCCRE and constant orientations (`myb/frames.h`).
- `myb_density_function` and the layered exponential atmosphere (`myb/force.h`).
- `myb_gravity_field_gradient`: the gradient of the acceleration of a field of harmonics; analytic
  Jacobians of the ready-made terms and of the relativistic accelerations
  (`myb_term_*_jacobian`, `myb_rel_*_jacobian`); `myb_term_no_jacobian`,
  `myb_force_model_bind_parameter`.
- `myb/estim.h`: least squares over an arc (`myb_fit_orbit`), an extended Kalman filter
  (`myb_filter`), measurements with their derivatives (position, range, range rate, RA/Dec,
  azimuth and elevation from a `myb_station`).
- Output: exact printing of numbers (`myb_format_double`, `myb_format_fixed`,
  `myb_format_scientific`); the SPK writer `myb_spk_writer` (types 9, 13, 2, 3); `myb/ccsds.h`
  — CCSDS OEM, OPM and OMM, read and written; `myb_sp3_write`; tables of numbers
  (`myb_table_*`); in `ext/sgp4`: `myb_tle_format`, `myb_tle_from_omm`, `myb_omm_from_tle`.
- Breaks of the step (`myb_force_model_add_break`, `myb_term_radiation_breaks`): `myb_orbit`
  and the estimators end a step at every edge of a shadow; `myb_orbit_breaks`.
- `myb/events.h`: ready-made events (`myb_event_range`, `myb_event_range_rate`,
  `myb_event_plane`, `myb_event_altitude`, `myb_event_elevation`, `myb_event_eclipse`,
  `myb_event_shape`); `myb_shadow_margin` (`myb/force.h`).
- `myb/geometry.h`: `myb_ray_ellipsoid`, `myb_ellipsoid_nearest`, `myb_ellipsoid_normal`,
  `myb_illumination_angles`, `myb_polyhedron_ray`, `myb_polyhedron_solid_angle`.
- Small forces (`myb/force.h`): `myb_force_radiation_moving` (Poynting–Robertson),
  `myb_nongrav` with `myb_nongrav_marsden` and `myb_nongrav_power`, `myb_yarkovsky_drift`,
  `myb_force_emission`, `myb_force_tidal_lag`, `myb_wind_function`; terms `myb_term_nongrav`,
  `myb_term_tidal_lag` (`myb/env.h`).
- `myb/surface.h`: attitude laws (`myb_attitude_fixed`, `_orbital`, `_sun_nadir`, `_spin`),
  free-molecular flow (`myb_force_gas_plates`, `myb_gas_sphere_cd`), light on a mesh with
  self-shadowing and its torque (`myb_force_radiation_mesh`), terms
  `myb_term_surface_radiation`, `myb_term_surface_drag`, `myb_term_recoil`.
- `myb_orbit_kick`: impulsive manoeuvres; `myb_thrust_mass_at`.
- Installation (`cmake --install`), the CMake package `find_package(myb)`, `myb.pc`, the
  optional shared library (`MYB_BUILD_SHARED`), `examples/`.
- Estimation: biases of the measurements as unknowns (`myb_fit_options.bias_count`,
  `myb_observation.bias`, `myb_filter_set_bias`, `myb_filter_bias`); a smoother
  (`myb_filter_smooth`, `myb_filter_steps`, `myb_filter_step`); a measured point away from
  the centre of mass (`myb_measure_offset`).
- Forces: `myb_force_nongrav_delayed`; `myb_yarkovsky_response`, `myb_force_yarkovsky`,
  `myb_yarkovsky_season` (the seasonal effect on an eccentric orbit) and the term
  `myb_term_yarkovsky`; `myb_torque_tidal_lag`; `myb_term_field_indirect`;
  `myb_shadow_fraction_many`; `myb_refraction`, `myb_refraction_bending`,
  `myb_shadow_fraction_refracted`; `myb_force_planet_radiation_gradient`; `myb_gas_knudsen`,
  `myb_gas_bridge`; `myb_force_yorp`; the Jacobians `myb_term_nongrav_jacobian`,
  `myb_term_tidal_lag_jacobian`, `myb_term_empirical_jacobian`,
  `myb_term_planet_radiation_jacobian`, `myb_term_yarkovsky_jacobian` (with
  `myb_force_yarkovsky_gradient` and `myb_force_yarkovsky_season_gradient`),
  `myb_term_surface_radiation_jacobian`, `myb_term_surface_drag_jacobian`,
  `myb_term_recoil_jacobian`.
- `myb_orbit_options.integrator`: DOP853 in `myb_orbit`; `myb_env_clone`.
- Geometry: `myb_ellipsoid_limb`, `myb_ellipsoid_terminator`, `myb_view_body`,
  `myb_view_intercept`; rays on a mesh go through a tree of bounding boxes.
- Files: `myb_spk_writer_add_raw`; `myb_pck_writer` (binary PCK); in `ext/sgp4`
  `myb_sgp4_spk10_constants`, `myb_sgp4_spk10_packet`, `myb_sgp4_spk10_size`,
  `myb_sgp4_spk10_words`; clocks, accuracy codes and flags of SP3 (`myb_sp3_clock`,
  `myb_sp3_sigma`, `myb_sp3_accuracy`, `myb_sp3_flags`).
- CCSDS: the XML notation (read automatically, written by `myb_oem_write_as`,
  `myb_opm_write_as`, `myb_omm_write_as`); covariances of OEM (`myb_oem_covariance`), OPM and
  OMM; the OCM message (`myb_ocm_*`).
- Packaging: ELF symbol versions (`cmake/myb.map`); pkg-config files of the libraries of
  `ext/`.
- Bindings (`bindings/`): C++17 — the header `myb/myb.hpp`, target `myb::cpp`; Python — the
  package `myb` over the shared library (cffi); Rust — the crates `myb-sys` and `myb`. Each
  has the whole C interface and classes for the usual path; results agree with the C calls
  to the bit. The example `examples/low_orbit.cpp`.

### Changed
- `myb_state_function` is declared in `myb/kin.h` (`myb/observe.h` includes it).
- `myb/force.h` includes `myb/frames.h`.
- The headers of `ext/` mark their functions with `MYB_EXT_API` instead of `MYB_API`.
- `ext/sgp4/myb_sgp4.h` includes `myb/ccsds.h`.
- New trailing fields (zero keeps the former behaviour): `myb_object.mass_function` and
  `mass_context`, `myb_radiation_term.moving` and `atmospheres`, `myb_drag_term.wind` and
  `wind_context`, `myb_yarkovsky.conductivity`, `myb_nongrav_term.delay`,
  `myb_sp3_data.clocks` … `flags`, `myb_oem_segment_data.covariance_*`, the covariances of
  `myb_opm` and `myb_omm`.
- The covariance section of an OEM is now read and checked: a matrix that does not have 21
  numbers is `MYB_E_FORMAT` (the section used to be passed over).
- The shadow of two or more bodies on the disk of the source in the radiation terms is the
  union of what they cover, not the product of the fractions; with one body the result is
  the same to the bit.
- `myb/geometry.h` includes `myb/env.h` and `myb/frames.h`; `myb/estim.h` includes
  `myb/surface.h`.
- `myb_term_drag_jacobian` takes the gradient of a wind into account (by differences of the
  wind function); without a wind the result is the same to the bit.

## 1.1.0.0 — 2026-11-10 (no tag)

### Added
- Compact binaries to 3PN and 3.5PN with spins: spin-orbit and spin-spin terms, spin-induced
  quadrupoles, spin precession (`myb_pn_binary_evaluate`).
- Electron column of power-law profiles and the dispersive delay (`myb/observe.h`).
- `ext/troposphere` (IERS Conventions 2010, ch. 9) and `ext/ionosphere` (IONEX, single layer).
- Charging of a sphere in a plasma: OML currents, photoemission, floating potential
  (`myb/force.h`).

## 1.0.0.0 — 2026-11-10

First release: numerical foundation (double-double, compensated sums, correctly rounded
elementary functions), reading of data (numbers, NAIF text kernels), time scales, the two-body
problem and Lambert's problem, integrators (DOP853, RKN12(10), IAS15, Gauss–Jackson,
symplectic compositions), the N-body problem (direct summation, interaction graph, threads,
Wisdom–Holman family), SPK ephemerides and reference frames (IAU 2006/2000A, EOP, WGCCRE, PCK,
equinox-based frames), relativity at the first post-Newtonian order, gravity fields of
extended bodies and tides, non-gravitational forces, force models with variational equations
and covariances, the restricted three-body problem, B-plane targeting, effect budgets,
attitude and rotational dynamics, regularisation, low thrust, compact binaries to 2.5PN,
observables, magnetic fields and the Lorentz force; `ext/nrlmsise00`, `ext/sgp4`.
