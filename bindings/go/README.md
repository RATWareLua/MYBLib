# myb for Go

The Go face of [myb](https://github.com/RATWareLua/MYBLib), a C library for the motion of
objects in space: the package `github.com/RATWareLua/MYBLib/bindings/go/myb`, plain cgo over
the **static library** `libmyb.a`, which cgo finds through pkg-config (the module
`myb-static`). The module depends on no other module.

## One call, one task

A crossing from Go into C costs tens of nanoseconds whatever the work inside, so the package
is the coarse layer of the library only — every call is a whole task:

- **a force model assembled from data**: `ModelSpec` lists the bodies (their GM, radii,
  luminosities; states from an `Ephemeris`; orientations: fixed, the Earth's with observed EOP,
  a binary PCK, a WGCCRE model) and the terms — `PointMass`, `ThirdBody`, `Field` (a
  `GravityField`), `Drag` (an exponential atmosphere), `Radiation` (with shadows),
  `Schwarzschild`, `LenseThirring`, `DeSitter`, `Empirical` — and the parameters a fit
  estimates; `NewModel` turns it into the environment and the model of the library in one call;
- **a propagation** (`Model.Propagate`: states at many times, the transition matrices and the
  sensitivities, impulsive manoeuvres), **a search for events** (`Model.Events`: nodes,
  apsides, distances, altitudes, the elevation over a station, eclipses), **a fit of an orbit**
  (`Model.Fit`: positions, states, range, range rate, angles from stations, with the light
  time, systematic shifts, a priori information) — one call each;
- **files**: OEM messages written and read in both notations (`WriteOEM`, `ParseOEM`,
  `ReadOEMFile`), SPK files written (`WriteSPK`, `EncodeSPK`) and read (`Ephemeris`);
- **everything "per epoch" as one call over a slice** — the array forms of the C library,
  which give the bits of the single calls in a loop: `Ephemeris.States`, `TimeSystem.ConvertMany`,
  `FromUTCMany`, `ToUTCMany`, `TimesFromSeconds`, `CIP.TerrestrialMany`, `CIP.XYSMany`,
  `EOP.AtMany`, `Equinox.RotationMany`, `WGCCRE.RotationMany`, `PCK.RotationMany`,
  `EarthOrientation.Sample` (and the other orientations), `ApplyRotations`,
  `KeplerPropagateMany`, `KeplerSTMMany`, `StatesFromElements`, `ElementsFromStates`, ...
  The single calls are there as well.

**No arithmetic is done in Go.** Numbers are copied between Go and the structures of the
library and its functions are called in the order a C program calls them, so a number has the
bits of the same calls made in C (the Go compiler may fuse `x*y + z` into one operation on some
processors: nothing in the package gives it the chance). The test checks it on the orbit of
`examples/low_orbit.c` against the C program `tests/bindings/low_orbit_reference.c`: the times
of the 15 ascending nodes, the states there, the state a day later and the count of steps agree
bit for bit.

**Go functions** may be terms (`GoTerm`) and events (`GoEvent`): the library calls them back
through cgo at every evaluation, a separate and plainly slow path (a crossing per call,
hundreds of thousands of calls per orbit). A term reports a failure by returning an error; the
call of the library stops with its status (that of a `*myb.Error`, `MYB_E_DOMAIN` for another
error) and the error it returns wraps the term's. A panic in such a function is not
recovered: it ends the program.

**Ownership.** Every object of the library is owned by a Go value with a `Close` method;
`runtime.AddCleanup` releases it should `Close` be forgotten, as a safety net only. A model
keeps what it uses (an ephemeris, a field, a time system) until it is closed itself: closing
such an object earlier makes it unusable from Go at once (`ErrClosed`) and releases it when the
model lets go. Values may be used from several goroutines; a model serves one call at a time.

**Errors** are values, `*myb.Error`: the status of the library (`errors.Is(err, myb.ErrRange)`),
the index of the element of a slice that failed (`Index`), the line a parser names (`Line`).
Nothing panics for a failure of the library; arguments that the library would take as a broken
contract (an instant out of range, slices of different lengths, a body that the model does not
have) are refused with `MYB_E_DOMAIN` before they reach it.

## Use

Install the library (static; `-DMYB_BUILD_SHARED=ON` may be on too) and point pkg-config at
it:

    cmake --preset linux-release && cmake --build --preset linux-release
    cmake --install build/linux-release --prefix /opt/myb
    export PKG_CONFIG_PATH=/opt/myb/lib/pkgconfig       # holds myb-static.pc

cgo needs a C compiler of the GNU kind; take the Clang the library is built with. On Linux:

    export CC=clang

On Windows the library is built for the mingw-w64 target with the C runtime of the MSYS2
MINGW64 sysroot; cgo must compile and link for the same:

    set CC=clang --target=x86_64-w64-windows-gnu --sysroot=C:/msys64/mingw64
    set CGO_LDFLAGS=-fuse-ld=lld
    set PKG_CONFIG_PATH=C:/myb/lib/pkgconfig

The module is not published (the library lives on GitHub only): refer to the checkout in the
`go.mod` of the program,

    require github.com/RATWareLua/MYBLib/bindings/go v0.0.0
    replace github.com/RATWareLua/MYBLib/bindings/go => /path/to/MYBLib/bindings/go

```go
package main

import (
	"fmt"
	"log"

	"github.com/RATWareLua/MYBLib/bindings/go/myb"
)

func main() {
	const gm, radius, day = 3.986004418e14, 6378137.0, 86400.0
	c20 := -4.841695e-4 // -J2 / sqrt(5), computed where the bits matter (Go may fuse)
	field, err := myb.NewGravityField(gm, radius, 2, []float64{1, 0, 0, c20, 0, 0}, make([]float64, 6))
	if err != nil {
		log.Fatal(err)
	}
	defer field.Close()
	model, err := myb.NewModel(myb.ModelSpec{
		Bodies: []myb.Body{{ID: 1, GM: gm}},
		Terms: []myb.Term{
			myb.Field{Name: "gravity", Body: 0, Field: field, NMax: -1},
			myb.Drag{Name: "drag", Body: 0, Object: myb.Object{Mass: 500, Area: 4, Cd: 2.2},
				EquatorialRadius: radius, Layers: []myb.Layer{{Base: 400e3, Density: 3e-12, Scale: 58e3}}},
		},
	})
	if err != nil {
		log.Fatal(err)
	}
	defer model.Close()
	y0, err := myb.StateFromElements(gm, myb.Elements{P: 6.778e6, E: 0.001, I: 0.9, RAAN: 0.5, ArgP: 1})
	if err != nil {
		log.Fatal(err)
	}
	nodes, err := model.Events(myb.EventSearch{T0: 0, Y0: y0, TEnd: day, Events: []myb.EventSpec{
		myb.PlaneCrossing{Body: myb.Origin, Normal: [3]float64{0, 0, 1}, Direction: myb.Rising},
	}})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(len(nodes.Events), "ascending nodes; a day later", nodes.Final)
}
```

## Test

`test_go_binding` of ctest (registered when CMake finds Go — the cache variable `MYB_GO`, else
`go` on `PATH` or in `toolchains/go*/bin` of the source tree — and pkg-config) builds the
package against the build tree: the pkg-config file
`build/<preset>/bindings/go/pkgconfig/myb-static.pc`, the compiler and the linker of the build
passed to cgo (`tests/tools/run_go_test.cmake`). Go keeps its caches, its work and the temporary
files of the tests in `build/<preset>/go`; it downloads nothing (`GOTOOLCHAIN=local`,
`GOPROXY=off`) and keeps no telemetry (`GOTELEMETRY=off`). By hand, from `bindings/go`:

    go test ./...       # with PKG_CONFIG_PATH, CC (and CGO_LDFLAGS on Windows) as above

The bit-for-bit test of the low orbit runs when `MYB_GO_REFERENCE` names the program
`go_low_orbit_reference` of the build tree, and is skipped otherwise.
