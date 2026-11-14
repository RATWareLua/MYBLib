# Frames and geodesy

**English** | [Русский](frames-geodesy.ru.md) · [All passports](index.en.md)

## What it is

The machinery of reference frames: rotations that carry their time derivative, so states pass
between frames that turn; the rotation models of bodies published by the IAU Working Group on
Cartographic Coordinates and Rotational Elements (WGCCRE) in NAIF text-PCK form; the fixed
frames of NAIF frame kernels (with SPICE's built-in inertial frames as data); a graph of frames
that composes the rotation between any two frames through their common ancestor; orientation
providers, the common form in which rotation enters a force model; geodetic coordinates on an
ellipsoid of revolution in closed form; the topocentric East-North-Up frame; and the orbital
frames RTN, VNC and LVLH.

## Source

- WGCCRE: B. A. Archinal et al., *Celest. Mech. Dyn. Astron.* 130, 22 (2018), as NAIF's text
  PCK `pck00011.tpc` writes it (WGCCRE 2015 rotation models and radii).
- TK frames: NAIF "Frames Required Reading"; the built-in inertial frames after Lieske (1977,
  1979), the FK4 equinox offset 0.525″, Galactic System II, the obliquities 84381.448″ (J2000)
  and 84404.836″ (B1950) of the Explanatory Supplement (1992).
- Geodetic coordinates: H. Vermeille, *J. Geodesy* 85, 105 (2011), the closed form valid
  everywhere including the interior; no iteration.
- Orbital frames: RTN (= RSW, RIC), VNC, and LVLH as CCSDS defines it.
- Implemented from these documents; no third-party code. CSPICE N0067 and ERFA served only as
  test oracles.

## Interface

Header [`myb/frames.h`](../../include/myb/frames.h):

- rotations: `myb_rotation`, `myb_rotation_apply_position`, `myb_rotation_apply_state`,
  `myb_rotation_compose`, `myb_rotation_invert`, `myb_rotation_euler313`;
- WGCCRE: `myb_wgccre`, `myb_wgccre_create`, `myb_wgccre_destroy`, `myb_wgccre_frame`,
  `myb_wgccre_euler`, `myb_wgccre_rotation`;
- TK frames: `myb_tkframe`;
- providers: `myb_orientation`, `myb_orientation_function`, `myb_wgccre_orientation`,
  `myb_orientation_wgccre`, `myb_orientation_fixed`, `myb_graph_orientation`,
  `myb_orientation_graph` (and, in their passports, `myb_orientation_cip`,
  `myb_orientation_pck`, `myb_orientation_ck`);
- the graph: `myb_frame_graph`, `MYB_FRAME_GRAPH_MAX_DEPTH`, `myb_frame_graph_options`,
  `myb_frame_graph_create`, `myb_frame_graph_destroy`, `myb_frame_graph_add`,
  `myb_frame_graph_add_fixed`, `myb_frame_graph_add_tkframe`, `myb_frame_graph_rotation`;
- geodesy: `myb_ellipsoid`, `myb_geodetic_from_cartesian`, `myb_geodetic_to_cartesian`,
  `myb_topocentric_enu`;
- orbital frames: `myb_orbital_frame` (`MYB_ORBITAL_RTN`, `MYB_ORBITAL_VNC`,
  `MYB_ORBITAL_LVLH`), `myb_orbital_rotation`.

See also [Earth orientation](earth-orientation.en.md) and
[PCK, CK and SCLK](pck-ck-sclk.en.md); SPK across frames is in
[SPK ephemerides](spk-ephemerides.en.md).

## Inputs and outputs

- A rotation maps coordinates of frame A to frame B: x_B = M x_A, v_B = M v_A + dM/dt x_A; M is
  proper orthogonal, stored by rows; dM/dt in s⁻¹. Composition and inversion keep the
  derivative exact to rounding. `myb_rotation_euler313`: M = R3(a₂) R1(a₁) R3(a₀), R1 and R3
  rotations of the frame, with rates in rad/s.
- WGCCRE: α = a₀ + a₁T + a₂T² + Σ Aᵢ sin θᵢ, δ = d₀ + d₁T + d₂T² + Σ Dᵢ cos θᵢ,
  W = w₀ + w₁d + w₂d² + Σ Wᵢ sin θᵢ (degrees), T in TDB Julian centuries and d in days from the
  epoch (J2000.0 unless `BODYn_CONSTANTS_JED_EPOCH`); the phase angles are polynomials of
  degree `BODYb_MAX_PHASE_DEGREE` (default 1), b = n/100 for 0 < n < 1000 (the system's
  barycentre) and b = n otherwise. M = R3(W) R1(90° − δ) R3(90° + α) from the reference frame
  (`BODYn_CONSTANTS_REF_FRAME`, default 1 = J2000) to the body-fixed frame. w₁d is kept in
  double-double and reduced modulo 360°.
- TK frames: `TKFRAME_<id>_SPEC` = 'ANGLES' (with `_AXES`, `_UNITS`), 'MATRIX' (by columns) or
  'QUATERNION' (SPICE style, scalar first); the rotation relative frame → frame `id`, constant.
- Providers: `myb_orientation` holds the rotation inertial axes → body-fixed frame, the body's
  rotation time `ut` (UT1 for the Earth; t otherwise) and the pole coordinates x_p, y_p (0 if
  not modelled); t is on the scale the provider is registered with.
- Graph: frames named by the caller's codes, each given in its base by a provider or a
  constant rotation; a frame without a definition that others refer to is a root (J2000, code
  1, for NAIF kernels). Instants are on the scale named at creation; a provider registered
  with another scale gets the same instant on its own (e.g. TT for the Earth's rotation in a
  graph on TDB).
- Geodesy: ellipsoid a (m) > 0, 0 ≤ f < 1 (from data, e.g. WGS-84 a = 6378137 m,
  f = 1/298.257223563); longitude in (−π, π], latitude in [−π/2, π/2], height in m; the centre
  gives the north pole at height −b. ENU: body-fixed → East-North-Up at a geodetic longitude
  and latitude, constant (dM/dt = 0).
- Orbital frames: the rotation from the frame of the state (r, v) to RTN, VNC (velocity,
  normal, co-normal = V × N) or LVLH (z nadir, y = −pole, x = y × z); the derivative follows
  from the acceleration given (a central force gives a constant pole).

## Domain of validity

- WGCCRE: `MYB_E_NOT_FOUND` (no `BODYn_POLE_RA`, missing phase angles), `MYB_E_FORMAT` (wrong
  counts, more terms than phase angles), `MYB_E_NOMEM`. The tests cover ±100 years.
- TK frames: `MYB_E_NOT_FOUND` (no definition), `MYB_E_FORMAT` (unknown SPEC or UNITS, wrong
  counts, a matrix that is not a rotation).
- Graph: `MYB_E_DOMAIN` (a code of 0, id == base, a second definition, a scale name without a
  time system), `MYB_E_NOT_FOUND` (an unknown scale; two frames in different trees),
  `MYB_E_FORMAT` (a chain longer than `MYB_FRAME_GRAPH_MAX_DEPTH` = 32: a cycle); the
  statuses of the providers and of the time system come back unchanged.
- Geodesy: valid everywhere, including the interior of the ellipsoid (tested from −6000 km to
  4·10⁵ km); `MYB_E_DOMAIN` for an invalid ellipsoid or latitude, or non-finite input.
- Orbital frames: `MYB_E_DOMAIN` for r or v zero, r parallel to v, or non-finite input.

## Accuracy

| What | Reference | Result |
|---|---|---|
| WGCCRE, 18 bodies (pck00011), 126 rotations | CSPICE `sxform` | ≤ 1.1·10⁻¹⁰ — growing with the rotation rate (Venus 10⁻¹³, the Moon 4·10⁻¹³, the Earth 2.5·10⁻¹¹, Vesta 1.1·10⁻¹⁰): CSPICE's W in one double |
| W of a fast rotator (Vesta), ±95 years | MPFR, 256 bits | 10⁻¹⁵ rad |
| SPICE's inertial frames from `naif-inertial.tk` | CSPICE `pxform` | 1.1·10⁻¹⁶ |
| Graph: 126 rotations between seven frames | CSPICE `sxform_c` | 2·10⁻¹⁶; the derivative 4·10⁻¹⁶ of the turning rate |
| Geodetic coordinates, −6000 km … 4·10⁵ km | the coordinates that generated the points | angles 3.3·10⁻¹⁶ rad, heights 5.6·10⁻¹⁶ of a + \|h\| |
| The same near the surface | ERFA `eraGc2gde` (Fukushima 2006) | 6·10⁻¹⁶ (far from the surface ERFA loses up to 2.6·10⁻¹¹ rad) |
| Orbital frames | orthonormality; central differences along a conic | 4·10⁻¹⁶; 5·10⁻¹¹ s⁻¹ |
| Euler 3-1-3 compose/invert | the identity | < 10⁻¹⁵, derivative < 10⁻¹⁸ s⁻¹ (tolerances) |

Tests: [`test_wgccre`](../../tests/unit/test_wgccre.c),
[`test_mpfr_frames`](../../tests/mpfr/test_mpfr_frames.c) (W),
[`test_tkframe`](../../tests/unit/test_tkframe.c),
[`test_frame_graph`](../../tests/unit/test_frame_graph.c),
[`test_erfa_frames`](../../tests/erfa/test_erfa_frames.c) (geodesy),
[`test_earth_frames`](../../tests/unit/test_earth_frames.c) (orbital frames).

## Cost

Not measured.

## Limitations

- Geodetic coordinates on an ellipsoid of revolution only; the ENU frame is constant.
- TK frames are constant rotations. `naif-inertial.tk` holds B1950, FK4, GALACTIC, ECLIPJ2000,
  ECLIPB1950; NAIF's DE-xxx frames (codes 4–12, 14, 15, 19–21) and MARSIAU (16) are not
  included.
- A WGCCRE model is as accurate as its published polynomials and series; against CSPICE the
  difference is CSPICE's single-double W, not the library's.
- The graph answers on one scale and needs each frame defined once; a provider's own status is
  returned as is. It serves several threads if its providers do.
- `myb_topocentric_enu` has no measurement of its own (the station tests of
  [site displacement](site-displacement.en.md) use it); `myb_orientation_wgccre` only wraps
  `myb_wgccre_rotation`.

## Data

In the repository: `data/frames/naif-inertial.tk` (SPICE's built-in inertial frames as TK
frames; facts) and `tests/data/wgccre/pck00011.tpc` (NAIF's text PCK, copied unchanged).
Downloaded: `pck00011.tpc` (NAIF, fixed release, SHA-256 in `cmake/FetchData.cmake`). The
ellipsoid constants are the caller's.
