// SPDX-License-Identifier: BSD-3-Clause

// Package myb is the Go face of the C library myb (motion of bodies in space): plain cgo over
// the static library libmyb.a, found through pkg-config (the module myb-static).
//
// # One crossing, one task
//
// A call from Go into C costs tens of nanoseconds whatever the work inside, so the package
// is the coarse layer only: each call is a whole task. A force model is assembled from data —
// Go values that describe the bodies and the terms (point masses, third bodies from an
// ephemeris, a gravity field, drag, radiation, relativity, empirical accelerations) — and
// turned into the model of the library inside one call ([NewModel]); a propagation, a search
// for events, a fit of an orbit is one call more ([Model.Propagate], [Model.Events],
// [Model.Fit]); a file is read or written in one call ([ParseOEM], [WriteOEM], [WriteSPK]).
// Whatever goes "per epoch" or "per point" goes as one call over a slice: the ephemerides
// ([Ephemeris.States]), the rotations of frames ([CIP.TerrestrialMany], [EOP.AtMany],
// [EarthOrientation.Sample], ...), the time scales ([TimeSystem.ConvertMany]) and Kepler's
// problem ([KeplerPropagateMany], [StatesFromElements]) — the array forms of the C library,
// which give the bits of the single calls in a loop. The single calls are there as well.
//
// # No arithmetic
//
// The package computes nothing: numbers are copied between Go and the structures of the
// library, and its functions are called in the order a C program calls them. A number that
// comes back has the bits of the same calls made in C (the orbit of examples/low_orbit.c is
// tested bit for bit against the C program). The Go compiler may fuse x*y + z into one
// operation on some processors; nothing here gives it the chance.
//
// # Go functions in the library
//
// A term of a model ([GoTerm]) or an event ([GoEvent]) may be a Go function: the library
// calls it back through cgo at every evaluation — a separate, plainly slow path (a crossing
// per call, hundreds of thousands per orbit). A term reports a failure by returning an error:
// the call of the library stops with its status (that of a [*Error], MYB_E_DOMAIN for any
// other error), and the error returned wraps it. A panic in such a function is not recovered
// — it unwinds into the frames of C and ends the program.
//
// # Ownership and errors
//
// An object of the library is owned by a Go value with a Close method; a cleanup attached to
// the value (runtime.AddCleanup) releases it should Close be forgotten, as a safety net only.
// An object that a model uses (an ephemeris, a field, a time system) lives until the model is
// closed: closing it before makes it unusable from Go at once and releases it after. Values
// are safe for use by several goroutines; a model serves one call at a time.
//
// A failure of the library is an error value, [*Error], with the status of the library, the
// index of the element of a slice that failed and the line a parser names; nothing panics
// for a failure of the library. Arguments the library would take as a broken contract (a
// time out of its range, slices of different lengths) are refused with MYB_E_DOMAIN before
// they reach it.
package myb
