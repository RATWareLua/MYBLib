# Magnetic fields, plasma and charging

**English** | [Русский](magnetic-field-charging.ru.md) · [All passports](index.en.md)

## What it is

The electromagnetic side of the environment of an object. Magnetic fields of bodies from Gauss
coefficients, read from tables in the IGRF format with their secular variation; the Lorentz
force; the electric field of a plasma that co-rotates with a planet; the momentum of a radial
plasma stream (the solar wind). And the charge a conducting sphere takes in a plasma: the
orbit-motion-limited currents of drifting species, photoemission, secondary electrons, the two
limits of the sheath, the floating potentials that balance the currents, the Debye length.

## Source

- Internal field: V = a Σₙ (a/r)ⁿ⁺¹ Σₘ (gₙₘ cos mλ + hₙₘ sin mλ) Pₙᵐ(sin φ), Schmidt
  semi-normalised (no Condon–Shortley phase), B = −∇V; evaluated by the spherical-harmonic code
  of [`myb/grav.h`](../../include/myb/grav.h) with Pₙᵐ = P̄ₙₘ / √(2n + 1).
- IGRF-14 (IAGA Division V Working Group V-MOD), the table of NOAA NCEI: linear in time between
  the epochs, the secular variation after the last (the IGRF prescription).
- Orbit-motion-limited currents: H. M. Mott-Smith, I. Langmuir, *Phys. Rev.* 28, 727 (1926); for
  drifting plasmas E. C. Whipple, *Rep. Prog. Phys.* 44, 1197 (1981) — closed forms in erf and
  erfc, a series in M² below M = 10⁻³. erf and erfc are those of fdlibm (FreeBSD msun `s_erf.c`,
  SunPro licence; see THIRD_PARTY_NOTICES.md).
- Secondary electrons: E. J. Sternglass (Westinghouse Sci. Pap. 1772, 1954) in the form of
  Whipple (1981) and NASCAP, δ(E) = e² δ_max (E/E_max) exp(−2 √(E/E_max)). Thin sheath: the
  flux of D. Bohm (1949) for attracted ions.
- ppigrf 2.1.0 (MIT licence) serves as an oracle only.
- Implemented from the papers; no third-party code besides erf.

## Interface

Header [`myb/force.h`](../../include/myb/force.h):

- fields: `myb_magnetic_field` (opaque), `myb_magnetic_field_create`,
  `myb_magnetic_field_load_table`, `myb_magnetic_field_destroy`, `myb_magnetic_field_degree`,
  `myb_magnetic_field_eval` (with a `myb_gravity_work` of
  [spherical harmonics](spherical-harmonics.en.md));
- forces: `myb_lorentz_acceleration`, `myb_corotation_field`, `myb_force_plasma_stream`;
- charging: `myb_plasma_species`, `myb_photoemission`, `myb_secondary_emission`,
  `myb_charging_model` (`MYB_SHEATH_THICK`, `MYB_SHEATH_THIN`), `myb_charging_current`,
  `myb_charging_current_model`, `myb_charging_equilibrium`, `myb_charging_equilibria`
  (`MYB_CHARGING_MAX_EQUILIBRIA`), `myb_secondary_yield`, `myb_debye_length`.

## Inputs and outputs

- Coefficients in tesla, indexed n (n + 1) / 2 + m (n = 0 unused); the reference radius a is the
  caller's (6371.2 km for the IGRF); r in the body-fixed frame (m), φ the geocentric latitude;
  B in tesla in the body-fixed frame.
- IGRF table: columns g/h, n, m, then the coefficients at the epochs of the header line
  "g/h n m 1900.0 … 2025.0 2025-30" in nT, the last column the secular variation in nT/year;
  the field is built for one decimal year.
- Lorentz: a = (q/m) (E + v × B), q/m in C/kg, E in V/m, B in T, all in one inertial frame.
  Corotation: E = −(Ω × r) × B.
- Plasma stream: ρ = ρ₀ (r₀/|r|)² flowing at `speed` away from the source; w = speed r̂ − v;
  a = (c_d/2) (A/m) ρ |w| w (c_d ≈ 2 for full momentum transfer).
- Species: charge (C), mass (kg), density (m⁻³), temperature as kT (J), drift relative to the
  body (m/s, ≥ 0). Photoemission: the emitted current density of the sunlit cross-section
  (A/m²), kT of the photoelectrons, e. The potential φ (V) is relative to the undisturbed
  plasma; a positive current charges the body positively. The charge of the floating sphere is
  Q = 4π ε₀ a φ, its q/m feeds `myb_lorentz_acceleration`.
- Secondaries: δ_max, E_max (J), kT of the secondaries (J); a fraction (1 + χ_s) e^(−χ_s) escapes
  a positive potential; for Maxwellian electrons the impact energies are weighted by their flux.
- Sheath: `MYB_SHEATH_THICK` (default) — a sphere small against the Debye length, OML currents;
  `MYB_SHEATH_THIN` — a ≫ λ_D: an attracted species brings no more than its current at zero
  potential (the drift included), attracted ions at least the flux of Bohm; a repelled species
  is held back by e^(−χ).
- `myb_charging_equilibrium`: the root of the current, bracketed outwards from 0 in steps of the
  largest kT/|q| and bisected to the resolution of the doubles. `myb_charging_equilibria`: roots
  over ±200 × the largest kT/|q| in steps of an eighth of it, bisected to adjacent doubles; the
  stable ones (the current falls with φ) in increasing order, each the last double with a
  positive current, at most 3.
- Debye length: √(ε₀ / Σ q² n / kT); 0 for no species.

## Domain of validity

- `myb_magnetic_field_create`: `MYB_E_DOMAIN` (degree < 1, radius ≤ 0, coefficients not
  finite), `MYB_E_NOMEM`. `myb_magnetic_field_load_table`: `MYB_E_FORMAT` (the line in the
  diagnostics), `MYB_E_RANGE` for a year before the first epoch, `MYB_E_NOMEM`. There is no upper
  limit: after the last epoch the secular variation goes on linearly (IGRF-14 gives it for
  2025–2030). `myb_magnetic_field_eval`: `MYB_E_DOMAIN` for r = 0 or not finite; the work must
  be made for at least the degree.
- Charging: `MYB_E_DOMAIN` for invalid species or radius, an unknown sheath, currents of one sign
  only (no balance), no root within ±200 kT/|q|; `MYB_E_NOCONV`.
- The thick limit needs a ≪ λ_D, the thin one a ≫ λ_D; between them (a ~ λ_D) neither holds.

## Accuracy

Tests [test_charge](../../tests/unit/test_charge.c) and
[test_charging](../../tests/unit/test_charging.c):

| What | Reference | Result |
|---|---|---|
| IGRF-14 from the NOAA table, 160 points 2000–2025, 1–4 Earth radii | ppigrf 2.1.0 | 1.4·10⁻¹⁵ of \|B\| |
| A dipole g₁₀ | its closed form | 8·10⁻¹⁶ |
| A charge in a uniform field, 50 gyrations (IAS15) | the exact helix | 5·10⁻¹⁵ |
| E × B drift | the straight drift line | exact |
| A proton in a dipole of Jupiter's scale, 2000 s | the speed and p_φ conserved | 3·10⁻¹⁶ |
| Corotation field | E ⊥ B, its drift the part of Ω × r across B | 10⁻¹⁵ |
| Solar wind, 5 protons cm⁻³ at 400 km/s, A/m = 0.01 m²/kg; at 2 au | ½ c_d (A/m) ρ v²; a quarter of it | 10⁻¹⁵ |
| OML current, M = 0…8, χ from −2 to 70 | the OML cross-section integrated over the drifting Maxwellian (mpmath, 30 digits) | 6·10⁻¹⁶ |
| Electron–proton plasma of equal temperatures | 1 − W(e √(m_p/m_e)) = −2.50385 kT/e | 10⁻¹⁵ |
| erf, erfc | mpmath | within 2 ulp |
| A grain in the solar wind at 1 au with photoemission | — | floats a few volts positive; negative in the shadow |
| Secondaries of Kapton (δ_max 2.1 at 150 eV), 40 eV electrons, −60…+60 V, both sheaths | a trapezoid of 4·10⁵ steps | 6.9·10⁻⁹ |
| A strong emitter (δ_max 3 at 300 eV) | — | positive from T_e = 14.2 eV |
| `myb_charging_equilibria` without secondaries | `myb_charging_equilibrium` | one root, 10⁻¹² |
| Thin sheath, a = 10 m | Bohm against Boltzmann, −3.338775 kT/e | 10⁻¹² |
| Debye length | by hand | 10⁻¹⁴ |

## Cost

Not measured.

## Limitations

- The field of internal sources only; its time dependence is the linear one of the table.
- No ready-made term of a force model: the caller puts together the field, the charge and the
  Lorentz acceleration.
- Charging gives floating potentials, not the charging transient; the sphere is a conductor far
  from others (Q = 4π ε₀ a φ).
- Secondary emission: Maxwellian electrons, their drift neglected. No bistability was found for
  Maxwellian plasmas (none of a grid of two-temperature plasmas gave two balances), and none is
  claimed.
- Only the two limits of the sheath; the plasma stream is radial with ρ ∝ r⁻² and one speed.

## Data

- `igrf14coeffs.txt` — IGRF-14 Gauss coefficients 1900–2025 with the secular variation (IAGA,
  NOAA NCEI), downloaded by `cmake/FetchData.cmake` and checked by SHA-256; not kept in the
  repository. Needed for the IGRF only; any field can be given by its coefficients.
- `tests/data/magnetic/igrf_reference.txt` — the ppigrf reference values of the test.
