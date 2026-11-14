# Post-Newtonian binaries

**English** | [Русский](pn-binary.ru.md) · [All passports](index.en.md)

## What it is

The relative motion of two compact bodies in harmonic coordinates of the centre-of-mass frame
beyond the first post-Newtonian order: the conservative 1PN, 2PN, 3PN and 4PN terms (4PN with
the conservative part of the gravitational-wave tail), the radiation reaction at 2.5PN, 3.5PN
and 4PN (the dissipative tail), the spin-orbit coupling at the leading order (1.5PN for maximal
spins) and at the next two orders, the spin-spin coupling, the quadrupoles induced by rotation,
the precession of the spins and the spin-orbit radiation reaction; with the Peters–Mathews
change of the period. Each term is a bit of a mask; a pair of bodies in a system of N bodies can
carry the same terms.

## Source

- L. Blanchet, "Post-Newtonian theory for gravitational waves", *Living Rev. Relativ.*,
  revision of 2024 (eqs. calAcoeff, calBcoeff = eqs. 203–204 of the 2014 edition); to 2.5PN
  also L. E. Kidder, *Phys. Rev. D* 52, 821 (1995), eq. 2.2. The 3.5PN radiation reaction at
  its first relative order: Iyer & Will 1993–1995, harmonic gauge.
- 4PN: L. Bernard, L. Blanchet, A. Bohé, G. Faye, S. Marsat, *Phys. Rev. D* 97, 044037 (2018):
  the instantaneous terms (eqs. 2.2–2.4), the tail (eq. 4.4) as the Fourier series of eq. 4.13
  in the mean anomaly, its dissipative part (eq. 6.2).
- Spins: Kidder 1995 (eqs. 2.2c, 2.2e, 2.4; covariant spin supplementary condition); the
  spin-induced quadrupoles of E. Poisson, *Phys. Rev. D* 57, 5287 (1998); the higher
  spin-orbit terms of A. Bohé, S. Marsat, G. Faye, L. Blanchet, *Class. Quantum Grav.* 30,
  075017 (2013), eqs. 3.5, 3.8–3.10, for spin variables of conserved norm; the spin-orbit
  radiation reaction of C. M. Will, *Phys. Rev. D* 71, 084027 (2005), eq. B3b; the geodetic
  precession of Barker & O'Connell (1975).
- P. C. Peters, J. Mathews, *Phys. Rev.* 131, 435 (1963).
- The 1PN … 3.5PN coefficients were checked against a machine parse of the LaTeX source of
  Blanchet's review (arXiv:1310.1528). The 4PN, higher spin-orbit and radiation-reaction blocks
  were transcribed by machine from the LaTeX sources of the papers, and the transcription was
  checked before the C: the 4PN energy and angular momentum are conserved identically through
  1/c⁸ along the paper's own equations. No third-party code.

## Interface

Header [`myb/rel.h`](../../include/myb/rel.h):

- the structure `myb_pn_binary`; `myb_pn_binary_evaluate`, `myb_pn_binary_acceleration`,
  `myb_peters_period_derivative`;
- the flags `MYB_PN_1`, `MYB_PN_2`, `MYB_PN_25`, `MYB_PN_3`, `MYB_PN_35`, `MYB_PN_SO`,
  `MYB_PN_SS`, `MYB_PN_ALL` (the seven above), and, asked for by name, `MYB_PN_4`,
  `MYB_PN_4RR`, `MYB_PN_SO25`, `MYB_PN_SO35`, `MYB_PN_SORR`.

A pair among N bodies (`myb_system_add_pn_pair`, `myb_system_pn_pair`, `myb_system_pack_spins`,
`myb_system_unpack_spins`): [system of bodies](system-of-bodies.en.md). The 1PN part is the
two-body limit of the EIH equations: [relativity at 1PN](relativity-1pn.en.md).

## Inputs and outputs

- `myb_pn_binary`: `gm1`, `gm2` (m³/s²), `c` (m/s), `orders`, `r0` (m, the gauge length r0′ of
  the logarithms, needed by `MYB_PN_3` and `MYB_PN_4`; at 4PN both gauge constants r′₀ and r″₀
  are r0), `kappa1`, `kappa2` (the spin-induced quadrupole constants, ≥ 0: 1 for a Kerr black
  hole, 2–10 for neutron stars; 0 omits the quadrupoles).
- x = x₁ − x₂ (m), v = v₁ − v₂ (m/s); a (m/s²) the acceleration of body 1 relative to body 2,
  a = −(Gm/r²) [(1 + A) n + B v], the terms of order k divided by c²ᵏ.
- Spins as dimensionless χ_A = c S_A / (G m_A²) (|χ| ≤ 1 for black holes); the outputs
  dχ_A/dt (1/s); |χ_A| stays constant. `chi1`, `chi2` NULL: no spin; `dchi1`, `dchi2` NULL: not
  wanted.
- r0 changes the coordinates, not the physics: gauge-invariant quantities such as Ω_φ/Ω_r do
  not depend on it.
- The tail is evaluated on the osculating Newtonian orbit through the current state, as a
  Fourier series in the mean anomaly with up to 1024 harmonics.
- `myb_pn_binary_acceleration`: the same without spins, 3PN and 4PN (r0 unset).
- `myb_peters_period_derivative`: dP/dt = −(192π/5)(2π G m_c/(P c³))^{5/3} f(e) (s/s), NaN for
  invalid input.
- General relativity only (β = γ = 1).

## Domain of validity

- `MYB_E_DOMAIN`: r = 0, gm1 or gm2 ≤ 0, c ≤ 0, input that is not finite, unknown bits,
  `MYB_PN_3` or `MYB_PN_4` without r0 > 0, `MYB_PN_SS` with a negative κ, `MYB_PN_4` or
  `MYB_PN_4RR` when the osculating Newtonian orbit is unbound or has e above about 0.83 (more
  than 1024 harmonics: e = 0.8 at periapsis is accepted, 0.9 refused); `MYB_PN_4` in
  `myb_pn_binary_acceleration`.
- In a system the tails need the pair's osculating orbit to stay bound; the integration fails
  otherwise.

## Accuracy

Tests: [`test_pn_binary`](../../tests/unit/test_pn_binary.c),
[`test_pn_higher`](../../tests/unit/test_pn_higher.c), [`test_pn4`](../../tests/unit/test_pn4.c),
[`test_system`](../../tests/unit/test_system.c); references in `tests/data/pn/`.

| What | Reference | Result |
|---|---|---|
| Coefficients 1PN … 3.5PN, 40 states | machine parse of Blanchet's LaTeX, mpmath 40 digits (`test_pn_higher`) | 4·10⁻¹⁶ |
| 2PN: (K − 1 − 3x)/x² on near-circular orbits, ν = 0.25 and 0.01 | 27/2 − 7ν (Damour, Jaranowski, Schäfer 2000) (`test_pn_binary`) | differs by 88·x, a remainder linear in x |
| Circular inspiral at x = 5·10⁻⁵ | Peters' dΩ/dt (`test_pn_binary`) | ratio 1.0006 |
| PSR B1913+16, 200 periastron passages (IAS15) | Peters–Mathews for the eccentricity of the integrated orbit (`test_pn_binary`) | Ṗ = −2.40234·10⁻¹², ratio 1.000031; periastron 4.2256°/yr (observed 4.2266) |
| Kepler law Ω²(r) through 3PN, ν = 0.25 and 0.1 | Blanchet, eq. keplerlaw (`test_pn_higher`) | 3PN coefficient to 10⁻⁶ (Richardson) |
| Radiation reaction on a circular orbit | the 1PN correction −743/336 − 11/4 ν (`test_pn_higher`) | 4·10⁻⁴ |
| Inspiral from 200 to 190 Gm/c² | eq. rOmdot with −1751/336 − 7/4 ν (`test_pn_higher`) | 98.7 % of the correction (the rest is 2PN) |
| K at 3PN, ν = 0.25 and 0.1 | 135/2 + (−649/4 + 123/32 π²) ν + 7ν² (`test_pn_higher`) | 0.06 %; independent of r0 (× 1000) to 10⁻⁴ |
| Spin-orbit, test-mass limit | `myb_rel_lense_thirring` (`test_pn_higher`) | 4·10⁻¹⁶ |
| Spin-spin with the quadrupoles | conservation of E and J (`test_pn_higher`) | 10⁻¹⁷ |
| Spin-orbit at the leading order | Kidder's E_N + E_SO (`test_pn_higher`) | conserved to the second order in the spins (spins × 1/10 → error × 1/100); aligned spins do not precess (exact zero) |
| PSR J0737−3039B, geodetic precession | Barker–O'Connell; Breton et al. 2008 (`test_pn_higher`) | 5.0735°/yr, 1.3·10⁻⁵ from the formula; the prediction is 5.0734°/yr |
| 4PN tail on eccentric orbits | its Fourier series evaluated in mpmath (`test_pn4`) | 10⁻¹³ of the tail (of its periapsis value near the apoapsis of a very eccentric orbit) by the header, 2·10⁻¹³ by the plan; the test bound is 10⁻⁹ |
| Tail on circular orbits | closed forms (eqs. omtail, eomcirc) (`test_pn4`) | 10⁻¹² (bound) |
| 4PN coefficient of the energy of circular orbits | eq. Ecirc, −54.92 (`test_pn4`) | to 0.035 (off by 50 without the tail energy) |
| 4PN coefficient of the periastron advance | eq. Kcirc, 237.5 (`test_pn4`) | to about 1 (the tail's share is 32) |
| 4PN radiation reaction on a circular orbit | the factor 1 + 4π x^{3/2} (`test_pn4`) | 4π to 3·10⁻⁴ (bound) |
| Energy lost to the tail in a period, e = 0.5 | the tail flux from an independent Fourier analysis (`test_pn4`) | 2.4·10⁻⁶ |
| Spin-orbit coefficients of E(x) at 1.5, 2.5, 3.5PN | arXiv:1212.5520, eq. Eofx (`test_pn4`) | 1.38111456 (exact), 2.967553 (2.967556), 5.0703 (5.0671) |
| Spin-orbit radiation reaction on circular orbits | Kidder's factor −0.2610081306 (`test_pn4`) | 2·10⁻⁹ |
| A pair among three bodies, with and without spins | the binary itself (`test_system`) | 9·10⁻¹⁷ |

## Cost

The tail (`MYB_PN_4`, `MYB_PN_4RR`), about the square of the number of harmonics in operations,
on one core of the Windows host (header `myb/rel.h`):

| e | Harmonics | Time per evaluation |
|---|---|---|
| 0.1 | 25 | 0.05 ms |
| 0.5 | 132 | 0.5 ms |
| 0.7 | 360 | 3 ms |
| 0.82 | 900 | 18 ms |

The other terms: not measured.

## Limitations

- General relativity only (β = γ = 1).
- The tail is taken on the osculating Newtonian orbit (an error of 5PN), on bound orbits with e
  below about 0.83.
- The 3.5PN radiation reaction is at its first relative order (harmonic gauge); at 4PN the
  radiation reaction is the dissipative tail.
- Spins: covariant supplementary condition, spin variables of conserved norm; the spin-orbit
  radiation reaction is at the leading order and leaves the spins alone (Will, eq. 1.6).
- In a system the relative acceleration is shared by the masses (the Newtonian momentum is
  kept); `MYB_PN_1` only if the point masses do not already carry the 1PN terms of the pair.
- Phase 21 adds: the spin-orbit terms of the tail, the spin-orbit radiation reaction beyond the
  leading order, spin-spin terms beyond the leading order, and the tail on unbound orbits.

## Data

None for the library. The tests read `tests/data/pn/blanchet_coefficients.txt`,
`pn4_coefficients.txt`, `so_terms.txt` and `tail_terms.txt`, computed by the authors of myb from
machine transcriptions of the papers' LaTeX sources (sympy, mpmath); only the evaluated numbers
are stored.
