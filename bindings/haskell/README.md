# myb for Haskell

The Haskell face of [myb](https://github.com/RATWareLua/MYBLib), a C library for the motion of
objects in space. A cabal package with two layers:

- `Myb.C` — the whole C interface: every function, structure, enumeration and constant of
  `include/myb`, by the names of the C headers, generated from them (`tools/Generate.hs`). The
  libraries of `ext/` are the public sublibraries `myb:sgp4`, `myb:nrlmsise00`,
  `myb:troposphere`, `myb:ionosphere`, `myb:jb2008`.
- `Myb` — the path most programs take, in the manner of the language: kernels and time scales,
  Kepler, SPK ephemerides, gravity fields, the environment and the force model with its
  ready-made terms, terms written in Haskell, the orbit with states, events, impulses and the
  transition matrix, OEM messages written and read. An object of the library is released by its
  finalizer (`ForeignPtr`); a call that can fail returns `Result a`, that is
  `Either Error a` with the status of the library; terms, providers and events are Haskell
  functions.

No arithmetic is done in Haskell: a number has the bits of the same calls made in C. The test
builds the orbit of `examples/low_orbit.c` from `Myb` and from the C functions and compares the
15 ascending nodes, the final state and the 776 steps bit for bit.

The package needs GHC and cabal-install and nothing else: it depends on the boot packages of GHC
only (`base`, and `containers`, `directory`, `filepath`, `process` for the generator and the
test), and `cabal.project` names no repository, so neither `cabal update` nor a network is
needed. Checked with GHC 9.10.3 and cabal 3.16.1.0.

## Linking

The package links `libmyb` from the directory given with `--extra-lib-dirs`, the static
library or the import library of the shared one, whichever the linker finds there (with both in
one directory it takes the shared one; give a directory that holds `libmyb.a` alone to link
statically).

- **Linux and others**: the static library. A program needs nothing at run time.
- **Windows**: the shared library, `libmyb.dll` (the `windows-release` preset builds it). GHC on
  Windows brings its own mingw-w64 toolchain, built for the UCRT (`-D_UCRT`, the
  `api-ms-win-crt-*` DLLs), while the library is built by Clang against the MSYS2 MINGW64
  sysroot, for `msvcrt.dll`. The static archive does not link into a program of GHC: it needs
  the stack protector (`__stack_chk_guard`, `__stack_chk_fail`) of the newer mingw-w64 of MSYS2,
  which the toolchain of GHC does not have, and objects built for one C runtime are not to be
  mixed with another anyway. The DLL carries its own C runtime (it imports `KERNEL32.dll` and
  `msvcrt.dll` only), and nothing of a C runtime crosses the interface: the headers use
  `stddef.h` and `stdint.h` alone, no `FILE`, and what the library allocates it frees. Put the
  directory of `libmyb.dll` (and of `libmyb_sgp4.dll`, ...) on `PATH` to run a program.

## Use

Build the library (`cmake --preset linux-release`, `cmake --build --preset linux-release`;
on Windows the `windows-release` preset), then, from `bindings/haskell`:

    cabal build --extra-lib-dirs=/path/to/the/library
    cabal test  --extra-lib-dirs=/path/to/the/library     # with MYB_DATA_DIR, MYB_TEST_DATA_DIR

GHC and cabal are taken from `PATH`; another GHC is named with `-w /path/to/ghc` (its `hsc2hs`
is found next to it). A program of another project takes the package by its path in its
`cabal.project` (`packages: path/to/bindings/haskell`) and `build-depends: myb` (and
`myb:sgp4`, ... for the libraries of `ext/`, linked as `libmyb_sgp4`, ...). cabal warns that
the include directories of the package are outside it (`../../include`): the package is built
from the repository, as the other bindings are.

The example `examples/LowOrbit.hs` (the executable `myb-example`, built and run by the test):

```haskell
import Control.Monad.IO.Class (liftIO)
import Control.Monad.Trans.Except (ExceptT (..), runExceptT)
import Myb (State (..), Vec3 (..))
import qualified Myb
import System.Exit (exitFailure)

main :: IO ()
main = do
  let gm = 3.986004418e14
      radius = 6378137.0
  done <- runExceptT $ do
    field <- ExceptT (Myb.newGravityField gm radius 2 [1, 0, 0, -1.08263e-3 / sqrt 5, 0, 0] (replicate 6 0))
    env <- ExceptT (Myb.newEnvironment (Myb.timeFromSeconds 0))
    planet <- ExceptT (Myb.addOrigin env Myb.defaultBody {Myb.bodyId = 1, Myb.bodyGm = gm})
    model <- ExceptT Myb.newForceModel
    _ <- ExceptT (Myb.addField model "gravity" env planet field 0 Nothing)
    _ <- ExceptT (Myb.addDrag model "drag" env planet (Myb.ExponentialAtmosphere radius 0 [(400e3, 3e-12, 58e3)])
                    (Myb.Object 500 4 1.3 2.2))
    orbit <- ExceptT (Myb.newOrbit model Myb.defaultOrbitOptions)
    y0 <- ExceptT (pure (Myb.stateFromElements gm (Myb.Elements 6.778e6 0.001 0.9 0.5 1.0 0)))
    ExceptT (Myb.start orbit 0 y0 86400)
    node <- ExceptT (Myb.nextEvent orbit (\_ y -> pure (let Vec3 _ _ z = position y in z)) 86400 Myb.Rising)
    liftIO (print node)                                 -- Just (time, state), or Nothing
    ExceptT (Myb.stateAt orbit 86400)
  either (\e -> putStrLn (Myb.errorText e) >> exitFailure) print done
```

(`ExceptT` is from `transformers`, a boot package of GHC; the library of the package itself
needs only `base`. `ExceptT` only carries the `Either` of each call: nothing is thrown.)

## The modules

| module | |
|---|---|
| `Myb` | the upper layer, all of it |
| `Myb.Error` | `Status`, `Error`, `Result`, `statusName`, `errorText` |
| `Myb.Types` | `Vec3`, `State`, `Elements`, `Calendar`, `Time` and their structures of C |
| `Myb.Time` | `KernelPool`, `TimeSystem`: `convertTime`, `fromUtc`, `toUtc` |
| `Myb.Kepler` | `stateFromElements`, `elementsFromState`, `keplerPropagate` (pure) |
| `Myb.Ephemeris` | `Ephemeris`: `loadEphemerisFile`, `ephemerisState`, `segments` |
| `Myb.Gravity` | `GravityField`: `newGravityField`, `loadIcgemFile` |
| `Myb.Environment` | `Environment`, `Body`: `addOrigin`, `addFromEphemeris`, `addWith` (a Haskell provider), `bodyState` |
| `Myb.Force` | `ForceModel`: `addPointMass`, `addThirdBody`, `addField`, `addRadiation`, `addDrag`, `addSchwarzschild`, `addTerm` (a Haskell term), `addRaw`, `acceleration` |
| `Myb.Orbit` | `Orbit`: `start`, `stateAt`, `statesAt`, `nextEvent` (a Haskell event), `kick`, `partials` |
| `Myb.Message` | `writeOem`, `parseOem` and the reading of a message |
| `Myb.C` | the whole C interface; `Myb.C.<Header>` one module per header (`Myb.C.Orbit` is `myb/orbit.h`) |
| `Myb.C.Ext.<Name>` | a library of `ext/`, in the sublibrary `myb:<name>` |

The upper layers for the rest of the interface (N bodies, estimation, CR3BP, targeting,
attitude, observables, the writers of files) come in a second round for all the bindings at
once, as modules beside these (`Myb.NBody`, `Myb.Estimation`, ...). Until then they are taken
through `Myb.C`; the two mix freely: every object gives its `ForeignPtr` (`forceModelPtr`,
`environmentPtr`, ...), and `forceModelKeep`, `environmentKeep` tie the life of what a term or
a body added through `Myb.C` points to.

### The lower layer

- A function is a foreign import by its C name: `myb_orbit_state :: Ptr Myb_orbit -> CDouble
  -> Ptr CDouble -> IO Myb_status`; its address is `p'myb_orbit_state`. Every call is in IO and
  `safe`, since any of them may call back into Haskell.
- A structure `myb_x` is the record `Myb_x` with the fields `myb_x'<field>` (an array is a
  list), a `Storable` instance whose size, alignment and offsets hsc2hs takes from the header,
  and `zero`, the value of all bits zero: `zero {myb_point_mass_term'gm = 3.986004418e14}`. An
  object known by pointer only (`myb_orbit`) is an empty type.
- An enumeration is a newtype over the integer type of the header; its constants, those of
  anonymous enumerations and the integer `#define`s are pattern synonyms by their C names
  (`MYB_OK`, `MYB_ODM_XML`, `MYB_ENV_MAX_OCCULTERS`).
- A type of function pointer `myb_t` is `Myb_t = FunPtr (...)` with `mk'myb_t`, a Haskell
  function as one (free it with `freeHaskellFunPtr`), and `call'myb_t`.
- A foreign import passes no structure by value. The functions that take or return one
  (`myb_time_from_seconds`, `myb_timesys_convert`, ...) are Haskell functions of the same name
  over a C function of the package (`generated/cbits`) that takes the structure by pointer and
  calls the library; nothing is computed there. A callback that takes a structure by value
  (`myb_state_function`, ...) goes through a trampoline of the package: `mk'myb_state_function`
  returns the C function and the context to give it.

### Ownership and failures

An object is destroyed by its finalizer once it is unreachable; a force model keeps alive the
environments and fields its terms point to (the structures of the terms are its own), an
environment its ephemerides and the Haskell functions of its providers, an orbit its model, so
any of them may be dropped in any order. A finalizer of GHC is not certain to run when the
program ends; the system takes the memory back then. An object of the library has working
memory for one thread: use it from one thread at a time.

Nothing is thrown or caught by the package. A Haskell function the library calls (a term, a
provider, an event) reports a failure through its result (`Left Range`), which comes back as
the status of the call that was going on. An exception that leaves such a function is not
caught: it ends the program (GHC's handler of a call from C prints it and exits), as a broken
contract stops at once in the library (`MYB_CHECK`); it never unwinds through the frames of C.

## The generator

`tools/Generate.hs` (the executable `myb-generate` of the package) writes `generated/` and the
lists of modules and sublibraries in `myb.cabal` between the marks `-- BEGIN myb-generate` and
`-- END myb-generate`. The generated files are committed; to write them again after a change of
the headers, from `bindings/haskell`:

    cabal run myb-generate -- ../..                    # clang on PATH, or --clang <clang>, or MYB_CLANG
    cabal run myb-generate -- --check ../..           # exits with 1 if the files are out of date

The headers go through the preprocessor of Clang (`-E`, freestanding, with its line markers),
so each declaration is known by its header; a small parser of C declarations reads them into
structures, enumerations, function pointers and functions, and each header becomes a module of
hsc2hs. The text holds no size, offset or path: it is the same on every system, and hsc2hs
lays the structures out from the headers when the package is built. A declaration the reader
does not know stops it with a message; nothing is skipped in silence.

It is written in Haskell, not in Python as `make_cdef.py` is: the check that the files agree
with the headers then needs only the tools the package needs already (GHC, cabal, and the Clang
that builds the library), and the boot packages of GHC suffice for it.

## The test

`ctest -R test_haskell_binding` (registered when GHC and cabal are found: `MYB_GHC`,
`MYB_CABAL`, looked for on `PATH`, in `toolchains/` of the source tree or of the home directory,
and in `~/.ghcup/bin`) runs `tools/test.cmake`: the generator with `--check`, the test suite,
then the example. Everything cabal writes, its configuration, store and logs and the build, stays in
`build/<preset>/haskell` (`CABAL_DIR`). `MYB_HASKELL_LINK` (`static`, or `shared`, the default
on Windows) chooses the library; it is copied alone into `build/<preset>/haskell/lib`.

The test suite (160 checks) holds: the low orbit of the example bit for bit with the C calls;
failures as values with the status and the line of a parser; terms, providers and events as
Haskell functions, and their failures; the radiation and relativity terms; that a model keeps
its environment alive and both are finalized once dropped (watched through their finalizers);
time scales, an SPK set and a body of it, OEM messages in both notations; the lower layer (the
constants, a structure passed and returned by value, a callback that takes one, `myb:sgp4`);
and, run as a second process, that an exception in a term or an event ends the process.
