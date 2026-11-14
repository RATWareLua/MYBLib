# myb for Python

The Python face of [myb](https://github.com/RATWareLua/MYBLib), a C library for the motion of
objects in space. The package loads the **shared library** of myb (it has no compiled part of
its own) and gives two layers:

- `myb.lib`, `myb.ffi` — every function and structure of the C interface, by the names of the
  C headers. The declarations are generated from the headers (`tools/make_cdef.py`).
- classes and functions over the whole interface, with ownership (an object of the library is
  released when Python drops it, or by `close()` and `with`), errors as exceptions
  (`myb.Error`) and Python callables in place of function pointers. Every name is a name of
  the package itself (`myb.System`, `myb.fit_orbit`); the modules hold them by subject.

No arithmetic is done in Python: a number has the bits of the same calls made in C. The tests
check exactly that, area by area, against the same calls made through `myb.lib`.

## What is there

| module | header(s) | classes and functions |
|---|---|---|
| `myb` | time.h, io.h (kernels), kin.h, grav.h, env.h, model.h, orbit.h, kepler.h, ccsds.h (OEM) | `KernelPool` (files, text, `describe`, `numbers`, `strings`), `TimeSystem` (`convert`, `from_utc`, `to_utc` and their `_many`), `Ephemeris` (`state`, `state_in`, `state_many`, `state_in_many`, segments and their words, images in memory, `body`), `GravityField`, `Polyhedron`, `Environment` (bodies with a state and an orientation from a provider), `ForceModel` (the terms of env.h, `add_cr3bp`, `add_thrust`, `add_qlaw`, Python terms, parameters), `Orbit`, the elements and Kepler (`state_from_elements`, `kepler_propagate`, `kepler_stm`, the equinoctial elements, and the `_many` of each), time (`time_from_seconds`, `seconds`, calendars, and their `_many`), `write_oem`, `read_oem`, `read_oem_file`, `oem_header` |
| `myb.nbody` | dyn.h, system.h | `NBody` (graph, relativity, threads through a `Pool`, the FMM approximation, accelerations, invariants), `WisdomHolman`, `System` (orientation, field, rotation, tide, mass rate, radius; the levels field, figures, tide, distortion, PN pair, action — a Python callable or the ready-made light on a mesh —, frame dragging, guided bodies, contact; pack/unpack of rotations and spins; accelerations, inertia, partials; `advance`, `merge`, `bounce`), `SystemIntegrator` |
| `myb.estim` | estim.h | `Station`, `Measurement` (the ready-made ones: `measure_position`, `measure_state`, `measure_range`, `measure_range_rate`, `measure_radec`, `measure_azel`, `measure_light_*`, `measure_offset`; your own with `Measurement.custom`), `Observation`, `fit_orbit` (parameters, shifts, priors, fixed unknowns), `Filter` (`start`, `predict`, `update`, `get`, shifts, `steps`, `smooth`, `step`) |
| `myb.cr3bp` | cr3bp.h | `cr3bp_lagrange`, `cr3bp_jacobi`, `cr3bp_acceleration`, `cr3bp_jacobian`, `cr3bp_correct`, `cr3bp_richardson`, `cr3bp_monodromy`, `cr3bp_stability`, `cr3bp_manifold`, `cr3bp_propagate` (sections), `cr3bp_from_ephemeris`, `cr3bp_to_ephemeris`, `Cr3bpOrbit` (`from_state`, `correct`, `continue_`), `Cr3bpTorus` (`from_orbit`, `continue_`, `state`, `stability`, `manifold`), `Cr3bpTorusSpectrum` (with the normal pair of an elliptic torus) |
| `myb.targeting` | transfer.h, lambert.h, optimal.h, thrust.h, iod.h, relative.h, conjunction.h, uncert.h | `bplane`, patched conics, `shoot`, `lambert`, `IndirectProblem` (`solve`, `sample`), `nlp_solve`, `SimsFlanagan` (`mismatch`, `solve`), `Chain` (powered flybys too: `size`, `guess`, `evaluate`, `solve`), `chain_flyby`, `Thrust`, `equinoctial_rates`, `edelbaum`, `edelbaum_state`, `QLaw`, `QLawTerm`, `QLawFlight` (`rhs`, `fly`), `iod_gibbs`, `iod_herrick_gibbs`, `iod_gauss`, `iod_laplace`, `iod_double_r`, `iod_gooding`, `relative_*`, `cw_stm`, `roe_*`, `approach_find`, `approach_orbits`, `encounter_plane`, `collision_probability*`, `conjunction_find`, `collision_rate`, `conjunction_long`, `Variational`, `covariance_transform`, `covariance_rtn`, `Rng`, `sample_statistics` |
| `myb.attitude` | attitude.h | quaternions (`quat_*`), `euler_to_matrix`, `euler_from_matrix`, `principal_axes`, `rigid_body_rates`, `core_mantle_rates`, the torques, `FreeRotation`, `rotation_split` (a torque in Python) |
| `myb.observe` | observe.h | `LightModel` (`light_time`, `round_trip`), `shapiro_delay`, `aberration`, `light_deflection`, `direction_radec`, `plasma_group_delay`, `plasma_column` |
| `myb.files` | io.h, kin.h, ccsds.h, ck.h, ext/sgp4 | `parse_double`, `format_*`, `Table`, `write_table`, `SpkWriter` (types 9 and 13 from states, 2 and 3 fitted to a provider, any type from its words), `Sp3`, `write_sp3`, `read_opm`, `write_opm`, `read_omm`, `write_omm`, `Ocm`, `write_ocm`, `write_ocm_line`, `Sclk`, `Ck` (with `orientation`, a provider), `CkWriter` (types 1–6), `tle_parse`, `tle_format`, `tle_epoch`, `tle_from_omm`, `omm_from_tle`, `tle_from_state`, `tle_fit` |
| `myb.frames` | frames.h | `Cip`, `Eop`, `Equinox`, `Wgccre`, `Pck`, `PckWriter`, `FrameGraph`, `orientation_fixed`, `orientation_many`, `rotation_apply_position_many`, `rotation_apply_state_many`: each call at one instant with its array form |

## Conventions

- A state is six floats (m, m/s); a vector a tuple of three. Matrices of six and more come out
  as lists of rows (`Orbit.partials`, covariances); the 3 x 3 matrices of attitude.h are nine
  numbers, rows one after another, as in C. Inputs take rows or flat numbers alike.
- Series — the `_many` array forms, the positions of many bodies, the epochs of a writer —
  come in as any sequence (of numbers, of tuples; an `array.array('d')` is taken without a
  copy) and go out as **lists**: of floats, of tuples, of instants. An array form is one
  crossing into C for the whole series.
- Instants are `myb_time` values; where instants are asked for, seconds past J2000 are taken
  too (made instants by `myb_time_from_seconds`, the same bits as `myb.time_from_seconds`).
  Calendar labels are `(year, month, day, hour, minute, second)`. Rotations are
  `myb.Rotation(m, dm)`, orientations `myb.Orientation(rotation, ut, xp, yp)`.
- A provider — of states (`myb_state_function`), of orientations, an attitude law — is a
  Python callable (`f(seconds) -> six numbers`, `f(seconds) -> Rotation`), a body of an
  `Ephemeris` (`eph.body(...)`), or a `myb.Provider(function, context, *kept)` of the C
  interface (`Pck.orientation`, `Ck.orientation`, `orientation_fixed`, ...).
- The structures of messages and elements are mappings with the names of the C members
  (`read_opm(text)["state"]`, `tle_parse(...)["mean_motion"]`).

## Errors and callables

A failed status raises `myb.Error`: `status` (`lib.MYB_E_DOMAIN`, ...), `name`, the `message`
of a parser with its line; `index`, for an array form, the element that failed; `partial`,
where the library leaves what it reached (`MYB_E_NOCONV` of the solvers, the state of a fit),
that. An exception raised in a Python callable that C code called — a term, a provider, an
event, a measurement, an action, the function of a solver — is kept, the callback reports a
failure, and the exception comes out of the call of the library that was going on. The
package catches nothing else. A broken contract of the C interface (an argument that the
headers say must not be NULL, a zero quaternion) stops the process with a message, as it does
in C.

## What stays on cffi only

Reached through `myb.lib` and `myb.ffi` (the two layers mix: every class has `pointer`, a
`Provider` takes any function of the C interface):

- the integrators of integ.h as such: an `NBody` integrated by RKN12, Gauss-Jackson or the
  symplectic compositions (`myb_nbody_acceleration` has their signature), a `System` with
  contacts integrated by DOP853 as a first-order system, `myb_qlaw_flight_rhs` under another
  integrator — the package has the IAS15 of a system (`SystemIntegrator`) and the integrators
  inside `Orbit`, `Variational`, `QLawFlight.fly`;
- what moves a station on its body (`myb_station_motion` of site.h: tides, loads), given to
  `Station` as cdata;
- allocators: every object is made with the default one (malloc);
- SGP4 itself and the SPK type 10 words of ext/sgp4 (`myb.extension("sgp4")`).

## Use

Build the library with `-DMYB_BUILD_SHARED=ON` (the `*-release` presets do), then

    pip install ./bindings/python
    export MYB_LIBRARY=/path/to/libmyb.so.1        # libmyb.dll on Windows

or copy the library into `myb/_lib/` before building a wheel, or install it where the system
looks for shared libraries. The libraries of `ext/` are found next to it
(`myb.extension("sgp4")`), or in the directories of `MYB_LIBRARY_PATH`.

```python
import math, myb

gm, radius = 3.986004418e14, 6378137.0
field = myb.GravityField(gm, radius, 2, [1, 0, 0, -1.08263e-3 / math.sqrt(5), 0, 0], [0] * 6)
env = myb.Environment()
planet = env.add_body(1, gm)
model = myb.ForceModel()
model.add_field("gravity", env, planet, field)
model.add_drag("drag", env, planet, mass=500.0, area=4.0, cd=2.2,
               layers=[(400e3, 3e-12, 58e3)], equatorial_radius=radius)
orbit = myb.Orbit(model)
y0 = myb.state_from_elements(gm, (6.778e6, 0.001, 0.9, 0.5, 1.0, 0.0))
orbit.start(0.0, y0, 86400.0)
node = orbit.next_event(lambda t, y: y[2], 86400.0, direction=1)   # (time, state) or None
print(orbit.state(86400.0))

# one crossing into C for a whole series
states = myb.kepler_propagate_many(gm, y0, [60.0 * k for k in range(1000)])

# an orbit of two bodies fitted to positions
two_body = myb.ForceModel()
two_body.add_point_mass("earth", gm)
fit = myb.fit_orbit(two_body, 0.0, states[0],
                    [myb.Observation(60.0 * k, myb.measure_position(), s[:3], (1.0,) * 3)
                     for k, s in enumerate(states[1:30], 1)])
print(fit.rms, fit.state)
```

Tests: `python -m unittest discover -s bindings/python/tests` with `MYB_LIBRARY` set (ctest:
`test_python_binding`) — 53 tests in 9 files, one per area.
