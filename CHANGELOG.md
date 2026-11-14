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
- `myb_budget_linear` (`myb/budget.h`): the budget of effects to first order — one
  propagation with a variation per term instead of a propagation per term; no floor of
  round-off.
- `myb/regular.h`: the time element of the KS equations (`MYB_REGULAR_KS_ELEMENT`,
  `myb_ks_element_from_cartesian`, `myb_ks_element_to_cartesian`, `myb_ks_element_rhs`); the
  Sundman transformation with any exponent α ≥ 0.
- The XML notation of the CCSDS messages: values in CDATA sections; a document type
  declaration with an internal subset is passed over.
- `myb/system.h`: a system of bodies with models (`myb_system`) — over the point masses of
  `myb_nbody`, levels of the graph: the field of a body on a point, the quadrupoles of two
  bodies, the tide by band with lags, the distortion of a body, the terms beyond 1PN of a
  pair (with spins), actions of any kind and light on the mesh of a body; the rotation of
  bodies (a mantle with a fluid core) integrated with the motion; encounters, collisions,
  rebound and merger.
- `myb_ias15_options.groups`, `group_count`: groups of coordinates of different kinds, each
  setting the step and judging the convergence by itself.
- `myb/iod.h`: initial orbit determination — `myb_iod_gibbs`, `myb_iod_herrick_gibbs`,
  `myb_iod_gauss` (from directions, refined by Newton's method).
- `myb/relative.h`: relative motion — the state in the frame of the chief, the linear map
  for any conic (`myb_relative_stm`), Clohessy–Wiltshire (`myb_cw_stm`), a two-impulse
  rendezvous, relative elements.
- `myb/conjunction.h`: the closest approaches of two objects (`myb_approach_find`,
  `myb_approach_orbits`) and the probability of collision (`myb_encounter_plane`,
  `myb_collision_probability`, `_chan`, `_max`).
- `myb/kepler.h`: `myb_mean_elements_average`; `ext/sgp4`: `myb_tle_fit`,
  `myb_tle_from_state`.
- `myb/cr3bp.h`: multiple shooting and pseudo-arclength continuation of families of periodic
  orbits (`myb_cr3bp_orbit_from_state`, `_correct`, `_continue`).
- `myb/thrust.h`: the Q-law of low thrust — `myb_qlaw`, `myb_qlaw_direction`,
  `myb_qlaw_arrived`, the model term `myb_term_qlaw` with an effectivity threshold; the same
  engine with the mass in the state — `myb_qlaw_flight_rhs`, `myb_qlaw_fly`.
- `myb/conjunction.h`: `myb_conjunction_find` — approaches with the covariances carried to
  them by the transition matrices of the orbits, and the probability of collision at each.
- `myb/optimal.h` (new): the indirect method for low thrust — `myb_indirect_solve` (the
  least propellant in a given time, the least time), `myb_indirect_sample`; the direct method
  of Sims and Flanagan — `myb_sims_flanagan_mismatch`, `myb_sims_flanagan_solve`; a solver of
  nonlinear programs — `myb_nlp_solve`.
- `myb/cr3bp.h`: quasi-periodic orbits as invariant circles of the flow —
  `myb_cr3bp_torus_from_orbit`, `myb_cr3bp_torus_continue` (families of one Jacobi constant,
  one period or one rotation), `myb_cr3bp_torus_state`.
- Files and frames: SPK type 15; a graph of frames (`myb_frame_graph_*`) and a state across
  frames (`myb_spk_state_in`); `myb/ck.h` (new) — the clock of a spacecraft, SCLK
  (`myb_sclk_*`), reading CK types 1–3 (`myb_ck_pointing_at`, `myb_orientation_ck`) and
  writing them (`myb_ck_writer_*`); a table of the nutation, `myb_cip_tabulate`; the
  difference of two scales by a function and by a table (`myb_timesys_set_difference`,
  `myb_timesys_tabulate`), TDB−TT from the ephemeris (`myb_spk_time_difference`).
- `myb/system.h`: bodies with a given motion (`myb_system_add_guided`), a mass that changes
  (`myb_system_set_mass_rate`), the Lense–Thirring effect (`myb_system_add_frame_dragging`).
- `myb/weather.h` (new): a table of space weather (`myb_space_weather_*`) and the input of
  the models of the upper atmosphere (`myb_thermosphere_point_at`); in `ext/nrlmsise00` —
  `myb_nrlmsise00_density`, the model as a `myb_density_function`.
- `myb/site.h` (new): the displacement of a site after IERS 2010 chapter 7 —
  `myb_site_displacement` (the tide of the solid Earth, the pole tide, the loads), ocean
  loading from BLQ tables (`myb_ocean_loading_*`), a grid of coefficients
  (`myb_site_grid_*`); the kernel `data/frames/iers2010-site.tk`. A station that moves:
  `myb_station.motion`.
- `myb/grav.h`: the secondary waves of the ocean tide by admittance (`MYB_TIDES_ADMITTANCE`,
  `myb_tides_secondary_waves`), waves added to those loaded (`myb_tides_add_waves`).
- `myb/estim.h`: measurements with the light time — `myb_light_measure`,
  `myb_measure_light_range`, `_range_rate`, `_radec`, `_azel`.
- `ext/troposphere`: `myb_tropo_vmf3`, GPT3 (`myb_gpt3_*`), `myb_tropo_zwd_askne`,
  `myb_tropo_gradient_delay_parts`. `ext/ionosphere`: integrals along a ray
  (`myb_iono_ray_integrals`), the density as a layer over the maps and as a grid
  (`myb_iono_layer_density`, `myb_iono_grid_density`), the higher orders
  (`myb_iono_orders_of`, `myb_iono_phase_delay`, `myb_iono_code_delay`, `myb_iono_bending`).
- `cmake/FetchData.cmake`: the GPT3 grid; large files with `-DLARGE=ON`.
- Phase 22. `myb/ck.h`: evaluating CK types 4, 5, 6 and writing them
  (`myb_ck_writer_add_chebyshev`, `myb_ck_writer_add_packets`,
  `myb_ck_writer_add_minisegments`).
- `myb/site.h`: the atmospheric S1/S2 loading from files — `myb_site_atmosphere_load`,
  `_load_file` (GGFC, TU Wien stations), `myb_site_atmosphere_grid_load`, `_grid_load_file`
  (the grids of TU Wien); the grid s1_s2_s3_cm_noib_grid.dat in `cmake/FetchData.cmake`.
- `ext/ionosphere`: three-dimensional IONEX maps — `myb_ionex_dimension`,
  `myb_ionex_density_grid`.
- `ext/jb2008` (new): JB2008 implemented from the papers — `myb_jb2008`, readers of the index
  files SOLFSMY and DTCFILE (`myb_jb2008_solar_*`, `myb_jb2008_storm_*`) and
  `myb_jb2008_density` as a `myb_density_function`.
- `myb/iod.h`: `myb_iod_laplace`, `myb_iod_double_r`, `myb_iod_gooding`.
- `myb/integ.h`: compositions of 7 and 15 stages (`MYB_SYMPLECTIC_6_S7`,
  `MYB_SYMPLECTIC_8_S15`); an interpolant of the steps of RKN12(10)
  (`myb_rkn12_options.interpolant`).
- `myb/orbit.h`: `MYB_ORBIT_RKN12` for force models of position-only terms
  (`myb_force_model_set_position_only`, `myb_force_model_position_only` in `myb/model.h`).
- `myb/conjunction.h`: the probability of a long encounter — `myb_collision_rate`,
  `myb_collision_probability_long`, `myb_conjunction_long`.
- `myb/force.h`: secondary emission and the thin sheath in charging —
  `myb_secondary_emission`, `myb_secondary_yield`, `myb_charging_model`,
  `myb_charging_current_model`, `myb_charging_equilibria`, `myb_debye_length`; atmospheric
  entry — `myb_air`, `myb_air_heat_ratio`, `myb_newtonian_cp_max`,
  `myb_heating_sutton_graves`, `myb_ablation`, `myb_ablation_rate`, `myb_ablation_rhs`,
  `myb_ablation_mass`.
- `myb/surface.h`: Newtonian continuum flow on plates and meshes
  (`myb_force_newtonian_plates`, `myb_force_newtonian_mesh`, `myb_continuum`); the light of
  the moment with thermal inertia (`myb_force_radiation_mesh_lagged`).
- `myb/system.h`: mergers of rotating bodies and of bodies of PN pairs; derivatives with
  respect to the state and the parameters (`myb_system_partials`); a search for neighbours in
  `myb_system_advance` (the same events, only the pairs whose paths over the step can meet
  are examined).
- `myb/contact.h` (new): contact of a body with a plane, an ellipsoid or a mesh — the law of
  Hunt and Crossley with smooth Coulomb friction (`myb_contact_force`,
  `myb_contact_surface_nearest`, `myb_contact_restitution`, `myb_contact_damping`), the term
  `myb_term_contact` and the gap `myb_event_contact`; in `myb/system.h` —
  `myb_system_add_contact`; in `myb/geometry.h` — `myb_polyhedron_nearest`.
- `myb/meteoroid.h` (new): random impacts of micrometeoroids after Grün et al. —
  `myb_meteoroid_flux`, `myb_meteoroid_rate`, `myb_meteoroid_mass`, `myb_meteoroid_next`,
  `myb_meteoroid_propagate`.
- `myb/cr3bp.h`: the stability and the manifolds of tori — `myb_cr3bp_torus_stability`,
  `myb_cr3bp_torus_manifold`; `myb_cr3bp_propagate` with a surface of section.
- `myb/dyn.h`: the fast multipole method as an approximate level of the graph —
  `myb_nbody_set_approximation`, `myb_nbody_approximation` (off by default).
- `myb/rel.h`: `MYB_PN_4` (4PN with the conservative tail), `MYB_PN_4RR` (the dissipative
  tail), `MYB_PN_SO25`, `MYB_PN_SO35` (spin-orbit terms at the next two orders),
  `MYB_PN_SORR` (spin-orbit radiation reaction); `myb_system_add_pn_pair` accepts them.
- `myb/optimal.h`: insertion into an orbit with a free phase in the indirect method
  (`myb_indirect_arrival`, `myb_indirect_problem.arrival`, `phase`); the field of a force model
  in Sims–Flanagan (`myb_sims_flanagan.field`); chains of legs with flybys — `myb_chain`,
  `myb_chain_size`, `myb_chain_guess`, `myb_chain_evaluate`, `myb_chain_solve`;
  `myb_nlp_options.weight`.
- Phase 21. Array forms of the calls (`*_many`, 26 functions: one call for a series, the
  bits of the single call in a loop, the index of the element that failed): in `myb/time.h`
  `myb_time_from_seconds_many`, `myb_time_seconds_many`, `myb_calendar_to_time_many`,
  `myb_time_to_calendar_many`, `myb_timesys_convert_many`, `myb_timesys_from_utc_many`,
  `myb_timesys_to_utc_many`; in `myb/kin.h` `myb_spk_state_many`, `myb_spk_state_in_many`; in
  `myb/frames.h` `myb_cip_terrestrial_many`, `myb_cip_xys_many`, `myb_eop_at_many`,
  `myb_equinox_rotation_many`, `myb_equinox_teme_pef_many`, `myb_wgccre_rotation_many`,
  `myb_pck_rotation_many`, `myb_orientation_many`, `myb_frame_graph_rotation_many`,
  `myb_rotation_apply_position_many`, `myb_rotation_apply_state_many`; in `myb/kepler.h`
  `myb_kepler_propagate_many`, `myb_kepler_stm_many`, `myb_pv_from_elements_many`,
  `myb_elements_from_pv_many`, `myb_pv_from_equinoctial_many`, `myb_equinoctial_from_pv_many`.
- Bindings for Go (`bindings/go`: cgo over the static library, the coarse layer — one call,
  one task; the package `myb/sgp4`) and Haskell (`bindings/haskell`: the whole C interface
  `Myb.C` generated from the headers, the upper layer `Myb` with `ForeignPtr` and `Either`;
  on Windows over the shared library). The pkg-config modules `myb-static` and
  `myb_sgp4-static`, installed with the library, name the archives by their paths.
- The upper layers of the C++, Python, Rust and Haskell bindings over the rest of the
  interface — N bodies and systems, estimation, the CR3BP, targeting and optimal transfers,
  attitude, observables, files, the array forms — and the same areas in the coarse layer of
  Go; in Rust also the frames of `myb/frames.h` and `OemSegment.covariance_frame`,
  `covariances`.
- `myb/contact.h`: friction that sticks — `myb_contact_friction`
  (`MYB_CONTACT_FRICTION_SMOOTH`, `MYB_CONTACT_FRICTION_STICK`), `myb_contact_force_stick`;
  in `myb/system.h` the deflections of its springs in the state — `myb_system_pack_contact`,
  `myb_system_unpack_contact`.
- `myb/integ.h`: the timescale of the step of IAS15 — `myb_ias15_timescale`
  (`MYB_IAS15_TIMESCALE_ACCELERATION`, `MYB_IAS15_TIMESCALE_JERK`); contact under IAS15 in
  `myb_orbit` with it.
- `myb/model.h`: the scheme of the breaks — `myb_break_scheme` (`MYB_BREAKS_RETAKE`, the
  default, and `MYB_BREAKS_PREDICT`), `myb_force_model_set_break_scheme`,
  `myb_force_model_break_scheme`.
- `myb/rel.h`: `MYB_PN_SS3` (spin-spin terms at next-to-leading order), `MYB_PN_SORR5`
  (spin-orbit radiation reaction at next-to-leading order), `MYB_PN_SOTAIL` (the dissipative
  spin-orbit tail on quasi-circular orbits), `MYB_PN_TAILX` (the 4PN tails on any osculating
  conic).
- `myb/optimal.h`: powered flybys in the chains — `myb_chain_flyby`,
  `myb_chain.impulse_exhaust`, `myb_chain_result.impulse`.
- `myb/cr3bp.h`: the normal rotation of elliptic tori — the fields `normal_rotation`,
  `normal_eigenvalue`, `normal_frequency`, `normal_residual`, `normal_re`, `normal_im` of
  `myb_cr3bp_torus_spectrum`; tori about vertical Lyapunov orbits.
- Fuzzing: 27 harnesses of the readers in `tests/fuzz` (`MYB_FUZZ`; libFuzzer with ASan and
  UBSan in the preset `linux-fuzz`, `tests/fuzz/run_fuzzers.cmake`; elsewhere a replay of the
  seeds and of damaged variants as ctest tests, label `fuzz`). Benchmarks: `bench/`
  (`MYB_BENCH`, the program `myb_bench`, the test `bench_smoke`). Coverage: `MYB_COVERAGE`, the
  presets `linux-coverage` and `windows-coverage`, `tests/tools/coverage.cmake`.
- macOS: the toolchain `cmake/toolchains/clang-macos.cmake` (Homebrew's LLVM), the presets
  `macos-dev`, `macos-release`, `macos-asan`. CI: Linux on ARM64, macOS, fuzzing, coverage,
  Go and GHC, the consumers built with MSVC.
- Tests: `test_msvc_consumer` (the headers as C17, C++17 and C++20 under `cl`, the example
  and the C++ face built by `cl` against `libmyb.dll`), `test_rust_binding_msvc`,
  `test_go_binding`, `test_haskell_binding`; on Linux under ASan `test_python_binding_asan`,
  `test_rust_binding_asan`.
- `docs/models/`: 45 model passports and an index, in English and Russian.
### Changed
- `MYB_TIDES_ALL` includes `MYB_TIDES_ADMITTANCE`: with a kernel that has a table of
  admittance and ocean waves loaded the result changes (the secondary waves are added).
- New fields at the end of structures (zero keeps the former behaviour):
  `myb_station.motion`, `myb_spk_body.graph`.
- `myb_measure_offset`: the derivatives include the turning of the offset after the orbit
  (they were those of the measurement at the point itself).
- `myb_ionex_vtec`: a value missing in the neighbouring map does not hinder a request
  exactly at the epoch of a map.
- `myb/kin.h` includes `myb/frames.h`; `myb/system.h` — `myb/kin.h`; `myb/estim.h` —
  `myb/site.h`; `myb/weather.h` — `myb/frames.h` and `myb/time.h`.
- `myb_state_function` is declared in `myb/kin.h` (`myb/observe.h` includes it).
- `myb/force.h` includes `myb/frames.h`.
- `myb/thrust.h` includes `myb/model.h`.
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
- `myb_sundman_rhs` no longer returns `MYB_E_DOMAIN` for an exponent other than 1, 3/2 and 2
  (only for a negative or non-finite one); for 1, 3/2 and 2 the result is the same to the
  bit. The description in `myb/regular.h` is corrected: steps uniform in the true anomaly
  belong to α = 2, not to 3/2.
- The C++ header: the functions that C calls moved to `myb::detail` and are pinned by
  `static_assert(noexcept(…))`.
- `myb_orbit` with DOP853 accepts `variational` (it was `MYB_E_DOMAIN`): the variations are
  integrated with the state.
- Phase 22. New trailing fields (zero keeps the former behaviour): `myb_offset_measure.gm`,
  `myb_system_light.theta`, `spin_axis`, `phases`, `myb_surface_drag_term.continuum`,
  `myb_rkn12_options.interpolant`; `myb_system_merger.spin`.
- `myb_system_merge` accepts bodies whose rotation is integrated and bodies of PN pairs (it
  was `MYB_E_DOMAIN`); `MYB_E_DOMAIN` remains for b with an integrated rotation and a without.
- CK: segments of types 4, 5, 6 are evaluated (they were loaded and passed over).
- `myb_ionex_load` reads three-dimensional maps (it was `MYB_E_DOMAIN`); `myb_iono_grid`:
  grids closed in longitude go round, a missing node is `MYB_E_RANGE`.
- `myb_ellipsoid_nearest` no longer returns `MYB_E_NOCONV` for points just off a sphere (its
  Newton steps fell below the resolution); the values it returned before are the same.
- `myb/myb.h` includes `myb/contact.h` and `myb/meteoroid.h`.
- `myb_wh_start` refuses (`MYB_E_DOMAIN`) a system with an approximation.
- `myb_pn_binary.r0` is also the gauge length of the 4PN logarithms (both constants r′₀ = r″₀).
- New trailing fields: `myb_indirect_problem.arrival`, `phase`, `myb_indirect_solution.phase`,
  `myb_indirect_point.costate`, `myb_sims_flanagan.field`, `myb_nlp_options.weight`.
- `myb_nlp_solve`: the penalty weight no longer grows once the constraints hold (the results of
  the earlier tests are the same to the bit).
- `myb/optimal.h` includes `myb/kin.h` and `myb/time.h`.
- Phase 21. New trailing fields (zero keeps the former behaviour): `myb_ias15_options.timescale`,
  `myb_orbit_options.timescale`, `myb_contact_law.kind`, `tangential_stiffness`,
  `tangential_damping`, `myb_chain.impulse_exhaust`, `myb_chain_result.impulse`, the
  `normal_*` fields of `myb_cr3bp_torus_spectrum`.
- `myb_system_partials`: the derivatives come from variational equations integrated with the
  state (the same signature; two bodies against the closed form 2.3·10⁻¹⁶ instead of
  1.8·10⁻¹¹ by differences); a system with contacts is still taken by differences of whole
  propagations.
- The FMM (`myb_nbody_set_approximation`) runs on the pool of `myb_nbody_set_pool` from 2048
  bodies on, with the bits of the evaluation without a pool.
- `myb_nlp_solve`: where the fall of the function is below 10⁻¹² of its values a step is
  accepted by the slope (approximate Wolfe conditions); problem 71 of Hock and Schittkowski
  moves by 10⁻¹⁴, a leg of Sims–Flanagan takes 436 evaluations instead of 435, and chains that
  stopped short of stationarity converge.
- `myb_system_add_contact` accepts the friction that sticks (3 coordinates a point after
  every other kind); `myb_contact_force` and `myb_term_contact` refuse it (`MYB_E_DOMAIN`).
- IAS15: a step whose sweeps diverge (a state at a node not finite) is rejected and tried four
  times shorter without calling the acceleration with such a state, under every timescale.
- `myb_sclk_decode` prints its fields as integers (a field of more than 18 digits is
  `MYB_E_RANGE`) and takes a NULL text of size 0 to give the length (`MYB_E_OVERFLOW`).
- The golden value of `test_repro` is `0xaa258e26488b41eb`; the prefixes 22, 20 and 19 are
  `0xe234728b95d060b4`, `0x7fc24de1c1d525be`, `0x4352d2d0bc185d9e` (the slope of the solver).
- Descriptions in the headers corrected after the audit of the model passports (the
  refusals of `myb_wh_start` and of the system, the XML readers, the light on a mesh, the
  ranges of JB2008 and of Mendes–Pavlis, the higher orders of the ionosphere, the errors of
  Lambert); no behaviour changes.

### Fixed
- Phase 21. Readers, found by fuzzing: `myb_sclk_decode` could overflow its text; times out of
  range from the kernels of time scales; WGCCRE epochs; SPK and PCK epochs beyond 2⁵¹ s and
  SPK type 10 packets; unchecked casts in the EOP reader; the length of CCSDS XML names was
  counted before the entities; TLE fields that the reader took, or the writer wrote, loosely
  (`ext/sgp4`).
- `myb_sclk_decode` printed large fields wrong with the C runtime of Windows (the `%.0f` of
  msvcrt gives 17 significant digits).
- `myb_fit_orbit` never reported convergence with exact observations (the weighted RMS falls
  to the rounding of the propagation); it also stops when the RMS is below the tolerance.
- The description of `myb_filter_start`: the covariance it reads is (6 + np + bias_count)²
  values.
- Python: the SP3 writer kept a pointer read back from a field of a structure, which did not
  keep the array of the satellites alive (it failed at random after the collector ran).
- Rust: with the MSVC toolchain the crates that depend on `myb-sys` did not find the import
  library; `build.rs` copies it into `OUT_DIR` as `myb.lib`.
- Go: the test keys the cache of cgo by a hash of the headers (a stale glue object read past
  the end of a structure that had grown); after installing another library run
  `go clean -cache`.

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
