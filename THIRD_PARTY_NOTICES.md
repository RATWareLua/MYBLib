# Third-party notices

myb is distributed under the BSD 3-Clause licence (`LICENSE`). It contains code and data
from the projects below, under their own licences.

## CORE-MATH

- Files: `third_party/core-math/binary64/sin/sin.c`, `third_party/core-math/binary64/cos/cos.c`
  (modified, see `third_party/core-math/README.myb.md`);
  test data `tests/data/core-math/sin.wc`, `tests/data/core-math/cos.wc` (samples).
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

## IERS

`data/time/utc-tai.tk` (rows since 1972 and the expiry date) follows the IERS file
`Leap_Second.dat` (https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat, updated through
IERS Bulletin C 72, July 2026).

## Test-only dependency (not distributed)

GNU MPFR and GMP are used by `tests/mpfr/` as the arbitrary-precision reference. They are
linked only into test executables and are not part of the library.
