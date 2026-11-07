# CORE-MATH in myb

Upstream: https://gitlab.inria.fr/core-math/core-math (MIT licence, see `LICENSE` here).
Imported revision: `660067dfd3e751c8f07acc31c899ff410aa7206e` (2026-11-01).

| File | Upstream path |
|---|---|
| `binary64/sin/sin.c` | `src/binary64/sin/sin.c` |
| `binary64/cos/cos.c` | `src/binary64/cos/cos.c` |

## Modifications

Exactly two, applied mechanically to each file:

1. One line added before the first `#include`:
   `#include "myb_coremath_prelude.h"`.
2. Every `__builtin_fma` replaced by `MYB_CM_FMA`.

`myb_coremath_prelude.h` maps `MYB_CM_FMA` to myb's exact integer `fma`
(`src/cmath/fma_soft.h`) and renames the entry points `cr_sin`/`cr_cos` to
`myb_cm_sin`/`myb_cm_cos`.

Why: CORE-MATH is correct only if `fma` is correctly rounded. On the x86-64-v2
baseline there is no hardware FMA, so `__builtin_fma` becomes a call to the C runtime,
and mingw-w64's `fma()` returns a wrong result for about 1 random input in 6
(measured on 2·10⁷ inputs against hardware FMA). The tests compare the results
with MPFR bit for bit.

## Updating

Copy the new upstream files and apply the two changes again; the check is:

```
diff <(sed 's/__builtin_fma\b/MYB_CM_FMA/g' upstream/sin.c) <(grep -v myb_coremath_prelude binary64/sin/sin.c)
```

which must print nothing.

The files are compiled with the library's numerical flags (`-ffp-contract=off`, no
fast-math), which CORE-MATH supports ("FMA contraction everywhere or nowhere").
