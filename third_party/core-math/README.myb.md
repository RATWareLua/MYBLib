# CORE-MATH in myb

Upstream: https://gitlab.inria.fr/core-math/core-math (MIT licence, see `LICENSE` here).
Imported revision: `660067dfd3e751c8f07acc31c899ff410aa7206e` (2026-11-01).

Functions (binary64, `src/binary64/<f>/<f>.c` and the local headers they include):
`sin cos tan atan atan2 asin acos exp expm1 log log1p sinh cosh tanh asinh acosh atanh cbrt
hypot` (`atan2/tint.h`, `log/dint.h`).

## Import and modifications

The files are produced by `import.cmake`, never edited by hand:

```
cmake -DUPSTREAM=<core-math checkout> -P third_party/core-math/import.cmake
```

It applies exactly three mechanical changes:

1. `#include "myb_coremath_prelude.h"` before the first preprocessor line of each `.c`.
2. `__builtin_fma` → `MYB_CM_FMA`.
3. `printf(` → `MYB_CM_PRINTF(`, `exit(` → `MYB_CM_EXIT(`.

`myb_coremath_prelude.h` maps `MYB_CM_FMA` to myb's exact integer `fma`
(`src/cmath/fma_soft.h`), drops debug output, turns the "unexpected worst case" path of
`atan2` (proved unreachable; upstream prints to stdout and calls `exit(1)`) into myb's
fail-fast check, and renames the entry points `cr_<f>` to `myb_cm_<f>`. Running the script
again on the same upstream revision leaves the tree unchanged.

Why the fma change: CORE-MATH is correct only if `fma` is correctly rounded. On the x86-64-v2
baseline there is no hardware FMA, so `__builtin_fma` becomes a call to the C runtime, and
mingw-w64's `fma()` returns a wrong result for about 1 random input in 6 (measured on 2·10⁷
inputs against hardware FMA).

## Run-time dependencies

Besides `memcpy`-class functions the vendored code calls: `fegetround`, `feholdexcept`,
`feupdateenv`, `fetestexcept`, `feclearexcept`, `fegetexceptflag`, `fesetexceptflag`,
`feraiseexcept` (flags are saved and restored; the environment on exit is the one on entry),
and `fmin`/`fmax` in `atan2` on absolute values only (no signed-zero or NaN ambiguity). No
libm function with platform-dependent rounding is used (`tests/tools/check_imports.cmake`).

## Verification

`tests/mpfr/test_mpfr_cmath.c` compares every function bit for bit with MPFR on random inputs
over the whole domain, special values, arguments near k·π/2 and a sample of the upstream
worst cases (`tests/data/core-math`). The files are compiled with myb's numerical flags
(`-ffp-contract=off`, no fast-math), which CORE-MATH supports.
