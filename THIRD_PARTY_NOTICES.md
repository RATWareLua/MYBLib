# Third-party notices

myb is distributed under the BSD 3-Clause licence (`LICENSE`). It contains code and data
from the projects below, under their own licences.

## CORE-MATH

- Files: `third_party/core-math/binary64/<f>/` for sin, cos, tan, atan, atan2, asin, acos,
  exp, expm1, log, log1p, sinh, cosh, tanh, asinh, acosh, atanh, cbrt, hypot (modified
  mechanically by `third_party/core-math/import.cmake`, see
  `third_party/core-math/README.myb.md`); test data `tests/data/core-math/*.wc` (samples of
  the upstream worst cases).
- Upstream: https://gitlab.inria.fr/core-math/core-math,
  revision `660067dfd3e751c8f07acc31c899ff410aa7206e`.
- Copyright (c) Paul Zimmermann, Alexei Sibidanov and the CORE-MATH contributors.
- Licence: MIT — full text in `third_party/core-math/LICENSE`:

```
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

The worst-case test data originates from the CORE-MATH project, which credits
https://www.vinc17.net/research/testlibm/ (Vincent Lefèvre).

## ERFA (Essential Routines for Fundamental Astronomy)

- Code: `tests/third_party/erfa/` — unmodified, compiled **only into tests** as the IAU SOFA
  reference (`tests/third_party/erfa/README.myb.md`).
- Data derived from ERFA tables: `data/time/tdb-tt-fb1990.tk` (coefficients of `dtdb.c`),
  `data/time/utc-tai.tk` (pre-1972 rows of `dat.c`), `data/frames/iau1980.tk` (the 106 terms of
  the IAU 1980 nutation in `nut80.c`, the polynomials of `prec76.c`, `obl80.c`, `gmst82.c`).
  The test oracle also uses ERFA's `ab.c` and `ld.c` (aberration, light deflection).
- Upstream: https://github.com/liberfa/erfa, revision
  `1a8044cde5b7763295d472a6443387239127c6c8`.
- Copyright (C) 2013-2023, NumFOCUS Foundation. Derived, with permission, from the
  International Astronomical Union's "Standards of Fundamental Astronomy" library.
- Licence: BSD 3-Clause — full text in `tests/third_party/erfa/LICENSE`.

## SciPy (coefficients of DOP853)

- File: `src/integ/dop853_coefficients.h` — the numerical coefficients of the Dormand-Prince
  8(5,3) method and its dense output (E. Hairer, S. P. Norsett, G. Wanner, *Solving Ordinary
  Differential Equations I*, 1993), copied as text from
  `scipy/integrate/_ivp/dop853_coefficients.py`; `src/integ/dop853.c` follows the structure
  of the error estimate and of the dense output in `scipy/integrate/_ivp/rk.py` (no code is
  copied). The Gauss-Radau constants of IAS15 are computed by myb itself
  (`tests/mpfr/test_mpfr_radau.c`).
- Upstream: https://github.com/scipy/scipy, branch `main`, retrieved 2026-11-07.
- Licence: BSD 3-Clause:

```
Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## OrdinaryDiffEq.jl (coefficients of RKN12(10))

- Files: `src/integ/rkn12_coefficients.h` (the 17-stage Runge-Kutta-Nystrom pair of
  J. R. Dormand, M. E. A. El-Mikkawy and P. J. Prince, IMA J. Numer. Anal. 7, 423, 1987, rounded
  to double) and `tests/data/rkn12/dprkn12.txt` (the same numbers as exact fractions), taken from
  `DPRKN12Tableau` in `lib/OrdinaryDiffEqRKN/src/rkn_tableaus.jl`. No code is copied.
- Upstream: https://github.com/SciML/OrdinaryDiffEq.jl, revision
  `6a6d63e6429f0e02e8eac75fc19f7ed7bf34d060` (retrieved 2026-11-07).
- Licence: MIT ("Expat"):

```
The OrdinaryDiffEq.jl package is licensed under the MIT "Expat" License:

> Copyright (c) 2016-2020: ChrisRackauckas, Yingbo Ma, Julia Computing Inc, and
> other contributors:
> 
> https://github.com/SciML/OrdinaryDiffEq.jl/graphs/contributors
> 
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
> 
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
> 
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.
```

## IERS

`data/time/utc-tai.tk` (rows since 1972 and the expiry date) follows the IERS file
`Leap_Second.dat` (https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat, updated through
IERS Bulletin C 72, July 2026).

`data/frames/iau2006-cip.tk` holds the coefficients of Tables 5.2a, 5.2b and 5.2d of the IERS
Conventions (2010), IERS Technical Note 36 (G. Petit, B. Luzum, eds.), as published by the
IERS Conventions Centre (https://iers-conventions.obspm.fr). `data/frames/eop-tides.tk` holds
the coefficient tables (IERS Conventions 2010, Tables 8.2, 8.3, 5.1a, 5.1b) of the IERS
programs `interp.f` (Ch. Bizouard, after R. Eanes) and `PMSDNUT2.F`, `UTLIBR.F` of the IERS
Conventions software; only the numbers are used, no code. The source files and their SHA-256
are listed in each kernel.

Test data `tests/data/eop/` are excerpts of the IERS EOP 20 C04 series (IERS Earth
Orientation Centre, Paris Observatory) and of `finals2000A.all` (IERS Rapid Service /
Prediction Centre); `tests/data/eop/reference.txt` was computed from them by myb's authors.

`data/grav/iers2010-tides.tk` holds the Love numbers, the tables of frequency-dependent
corrections (Tables 6.3, 6.5a-c), the permanent-tide and pole-tide constants of chapter 6 and
the secular pole of chapter 7 of the IERS Conventions (2010, version of 2018) — numbers only,
with the SHA-256 of the source PDF in the kernel. `tests/data/grav/fes2004_excerpt.dat` is an
unchanged excerpt of `fes2004_Cnm-Snm.dat` (FES2004 ocean tide model of F. Lyard et al.,
converted to Stokes coefficients for the IERS Conventions).

## Gravity field models (ICGEM)

`tests/data/grav/EGM2008_n20.gfc` is EGM2008 (N. K. Pavlis, S. A. Holmes, S. C. Kenyon,
J. K. Factor, U.S. National Geospatial-Intelligence Agency) truncated to degree 20, from the
International Centre for Global Earth Models (ICGEM, GFZ Potsdam); coefficients unchanged.
GOCO06s (A. Kvas et al., 2019, doi:10.5880/ICGEM.2019.002, CC BY 4.0) is downloaded by
`scripts/fetch-data.cmd` and not distributed.

## ILRS

`tests/data/sp3/ilrsa.lageos1.250607.excerpt.sp3` is an excerpt (every third epoch) of the
combined LAGEOS-1 orbit `ilrsa.orb.lageos1.250607.v80.sp3` of the International Laser Ranging
Service (ILRS combination of the analysis centres ASI, BKG, CNES, DGFI, ESA, GFZ, JCET, NSGF),
obtained from the EUROLAS Data Center (EDC, DGFI-TUM) on 2026-11-09; data of the ILRS are
freely available, with acknowledgement of the ILRS (M. R. Pearlman et al., Adv. Space Res. 30,
135, 2002).

## NRLMSISE-00 (ext/nrlmsise00)

The optional module `ext/nrlmsise00` (target `myb::nrlmsise00`, not part of libmyb) is
derived from the C release of NRLMSISE-00 by Dominik Brodowski, placed by him in the public
domain, of the model of J. M. Picone, A. E. Hedin and D. P. Drob (Naval Research Laboratory).
Notice under 17 U.S.C. 403: the model, its algorithms and its coefficient tables are U.S.
Government material, not subject to copyright protection. This module is not endorsed by, nor
affiliated with, the Naval Research Laboratory. Details and the original statements:
`ext/nrlmsise00/NOTICE`.

## SGP4 (ext/sgp4)

The optional module `ext/sgp4` (target `myb::sgp4`, not part of libmyb) is a C translation of
the SGP4/SDP4 code of D. A. Vallado, P. Crawford, R. Hujsak and T. S. Kelso ("Revisiting
Spacetrack Report #3", AIAA 2006-6753; version of 12 March 2020, distributed by CelesTrak for
free use), made from its Python translation in python-sgp4 by Brandon Rhodes
(https://github.com/brandon-rhodes/python-sgp4, MIT licence). The module is under the same MIT
licence; its verification data `tests/data/sgp4/SGP4-VER.TLE` and `tcppver.out` are those of
Vallado et al. as distributed with python-sgp4. Full text: `ext/sgp4/NOTICE`.

## IGRF and ppigrf

The International Geomagnetic Reference Field (IGRF-14, IAGA Division V Working Group V-MOD)
is downloaded at build time from NOAA NCEI (`igrf14coeffs.txt`, checked by SHA-256) and is not
part of the repository. The reference values `tests/data/magnetic/igrf_reference.txt` were
computed by the authors of myb with ppigrf 2.1.0 (MIT licence; https://github.com/IAGA-VMOD/ppigrf),
an independent implementation used as an oracle only.

## mpmath

The reference values `tests/data/elliptic/reference.txt` were computed by the authors of myb
with mpmath 1.4.1 (https://mpmath.org, BSD 3-Clause licence) at 60 digits; mpmath itself is
not part of myb or of its tests.

## NAIF (JPL)

`tests/data/wgccre/pck00011.tpc` is NAIF's generic text PCK (WGCCRE 2015 rotation models and
radii, https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc), a NASA/JPL data
product, copied unchanged. `data/frames/naif-inertial.tk` restates the definitions of the
inertial frames built into SPICE as given in NAIF's Frames Required Reading.

The SPK test files `tests/data/spk/t*.bsp` (including the two-line-element segments of
`t10.bsp`), `chain_*.bsp` and the reference states `tests/data/spk/reference.txt`,
`t10_reference.txt`, `tests/data/wgccre/reference.txt`, `tests/data/observe/spice_lighttime.txt`
were generated by myb's authors with the
NAIF SPICE Toolkit (CSPICE N0067, through spiceypy), used as an external oracle only; no SPICE
code is part of myb or of its tests. The layout of SPK type 10 and its evaluation (blending of
two element sets, TEME as SPICE defines it) follow NAIF's documentation of SPKE10 and ZZTEME.

## JPL Horizons

`tests/data/spk/eros.bsp` (SPK of asteroid 433 Eros, solution JPL#659),
`tests/data/spk/eros_horizons.txt` and `tests/data/frames/horizons_moon_sites.txt` are
outputs of the JPL Horizons system (https://ssd.jpl.nasa.gov/horizons/, Solar System
Dynamics Group, JPL/Caltech), obtained through its API on 2026-11-09.

Spacecraft trajectories `tests/data/horizons/*.txt` (Tesla Roadster, OSIRIS-REx, JWST) and the
periodic orbits `tests/data/cr3bp/jpl_periodic_orbits.txt` (NASA/JPL Three-Body Periodic Orbits
API, https://ssd-api.jpl.nasa.gov/periodic_orbits.api) were obtained from the JPL Solar System
Dynamics Group on 2026-11-10; the queries are recorded in the files.

## Test-only dependency (not distributed)

GNU MPFR and GMP are used by `tests/mpfr/` as the arbitrary-precision reference. They are
linked only into test executables and are not part of the library.
