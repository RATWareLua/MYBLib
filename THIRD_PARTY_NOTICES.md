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
  `data/time/utc-tai.tk` (pre-1972 rows of `dat.c`).
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

## Test-only dependency (not distributed)

GNU MPFR and GMP are used by `tests/mpfr/` as the arbitrary-precision reference. They are
linked only into test executables and are not part of the library.
