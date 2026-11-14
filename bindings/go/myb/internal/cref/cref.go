// SPDX-License-Identifier: BSD-3-Clause

// Package cref is the C reference of the tests of the Go binding: the tasks of the binding
// written as a C program writes them (cref.c), for the tests to compare bit for bit. Every
// function returns the status of the library as an int (0: MYB_OK).
package cref

/*
#cgo pkg-config: myb-static
#include <stdlib.h>
#include "cref.h"
*/
import "C"

import "unsafe"

func d(s []float64) *C.double {
	if len(s) == 0 {
		return nil
	}
	return (*C.double)(unsafe.Pointer(&s[0]))
}

// SystemRun: states (2 n each), events (7 each), the final state.
func SystemRun(gm, radius, y0, times []float64) (int, []float64, []float64, []float64) {
	n := len(y0)
	states := make([]float64, n*len(times))
	events := make([]float64, 7*64)
	final := make([]float64, n)
	var count C.size_t
	status := C.cref_system_run(d(gm), d(radius), C.size_t(len(gm)), d(y0), d(times),
		C.size_t(len(times)), d(states), d(events), &count, d(final))
	return int(status), states, events[:7*int(count)], final
}

// Invariants: the high words of the energies and the momenta (3 each) of count states.
func Invariants(gm, y []float64, count int) (int, []float64, []float64) {
	energy := make([]float64, count)
	momentum := make([]float64, 3*count)
	status := C.cref_invariants(d(gm), C.size_t(len(gm)), d(y), C.size_t(count), d(energy), d(momentum))
	return int(status), energy, momentum
}

// SystemPartials: y1, phi, the sensitivity to the gm of body 0.
func SystemPartials(gm, y0 []float64, t1 float64) (int, []float64, []float64, []float64) {
	n := len(y0)
	y1, phi, s := make([]float64, n), make([]float64, n*n), make([]float64, n)
	status := C.cref_system_partials(d(gm), C.size_t(len(gm)), d(y0), C.double(t1), d(y1), d(phi), d(s))
	return int(status), y1, phi, s
}

// WH: SABA3 without corrector; positions and velocities at the times.
func WH(gm, x, v []float64, h float64, times []float64) (int, []float64, []float64) {
	xs, vs := make([]float64, len(x)*len(times)), make([]float64, len(v)*len(times))
	status := C.cref_wh(d(gm), C.size_t(len(gm)), d(x), d(v), C.double(h), d(times), C.size_t(len(times)),
		d(xs), d(vs))
	return int(status), xs, vs
}

// Filter: innovations (3 each), the times and states of the steps (smoothed).
func Filter(gm float64, y0, covariance []float64, tEnd float64, times, positions []float64, sigma float64) (int, []float64, []float64, []float64) {
	n := len(times)
	innovations := make([]float64, 3*n)
	stepT := make([]float64, n+1)
	stepU := make([]float64, 6*(n+1))
	var steps C.size_t
	status := C.cref_filter(C.double(gm), d(y0), d(covariance), C.double(tEnd), d(times), d(positions),
		C.size_t(n), C.double(sigma), d(innovations), d(stepT), d(stepU), &steps)
	return int(status), innovations, stepT[:int(steps)], stepU[:6*int(steps)]
}

// Range: the range and its partials from a station on a sphere at the origin.
func Range(gm, longitude, latitude, radius float64, t, y []float64) (int, []float64, []float64) {
	values, partials := make([]float64, len(t)), make([]float64, 6*len(t))
	status := C.cref_range(C.double(gm), C.double(longitude), C.double(latitude), C.double(radius),
		d(t), d(y), C.size_t(len(t)), d(values), d(partials))
	return int(status), values, partials
}

// CR3BPFamily: the periods and the states (6 patches each) of count members.
func CR3BPFamily(mu float64, y0 []float64, period float64, patches int, step float64, count int) (int, []float64, []float64) {
	periods := make([]float64, count)
	states := make([]float64, 6*patches*count)
	status := C.cref_cr3bp_family(C.double(mu), d(y0), C.double(period), C.size_t(patches), C.double(step),
		C.size_t(count), d(periods), d(states))
	return int(status), periods, states
}

// Torus: the rotations of two members and of the second's spectrum the normal rotation, its
// eigenvalue, two components of the eigenfunction and the count of elliptic eigenvalues.
func Torus(mu float64, y0 []float64, period float64, points int, size, step float64) (int, []float64, []float64) {
	rotations, normal := make([]float64, 2), make([]float64, 6)
	status := C.cref_torus(C.double(mu), d(y0), C.double(period), C.size_t(points), C.double(size),
		C.double(step), d(rotations), d(normal))
	return int(status), rotations, normal
}

// Manifold: a state on the unstable manifold (side +1, epsilon 1e-6) at a phase.
func Manifold(mu float64, y0 []float64, period, phase float64) (int, []float64) {
	y := make([]float64, 6)
	status := C.cref_manifold(C.double(mu), d(y0), C.double(period), C.double(phase), d(y))
	return int(status), y
}

// LambertCell: the zero-revolution solution between two conics of mu at two times.
func LambertCell(mu float64, from, to []float64, departure, arrival float64) (int, []float64, []float64) {
	v1, v2 := make([]float64, 3), make([]float64, 3)
	status := C.cref_lambert_cell(C.double(mu), d(from), d(to), C.double(departure), C.double(arrival), d(v1), d(v2))
	return int(status), v1, v2
}

// Conjunctions: the count, the time and the probability of the first.
func Conjunctions(gm float64, y1, y2, covariance []float64, radius, t1, step float64) (int, int, float64, float64) {
	var count C.size_t
	var t, p C.double
	status := C.cref_conjunctions(C.double(gm), d(y1), d(y2), d(covariance), C.double(radius), C.double(t1),
		C.double(step), &t, &p, &count)
	return int(status), int(count), float64(t), float64(p)
}

// QLawFly: the state and the time of the arrival.
func QLawFly(mu float64, target, weight []float64, thrust, exhaust float64, y0 []float64, tLimit, step float64) (int, []float64, float64) {
	y := make([]float64, 7)
	var t C.double
	status := C.cref_qlaw_fly(C.double(mu), d(target), d(weight), C.double(thrust), C.double(exhaust), d(y0),
		C.double(tLimit), C.double(step), d(y), &t)
	return int(status), y, float64(t)
}

// Chain: the guess built in, its propellant and its constraints.
func Chain(mu float64, states, gm, periapsis, impulse, bounds []float64, unknowns, constraints int) (int, []float64, float64, []float64) {
	x, c := make([]float64, unknowns), make([]float64, constraints)
	var p C.double
	status := C.cref_chain(C.double(mu), d(states), d(gm), d(periapsis), d(impulse), d(bounds), d(x),
		C.size_t(unknowns), &p, d(c), C.size_t(constraints))
	return int(status), x, float64(p), c
}

// SimsFlanagan: the mismatch and the mass at the arrival.
func SimsFlanagan(mu float64, start, target []float64, t1 float64, throttles []float64) (int, []float64, float64) {
	m := make([]float64, 6)
	var mass C.double
	status := C.cref_sims_flanagan(C.double(mu), d(start), d(target), C.double(t1), C.size_t(len(throttles)/3),
		d(throttles), d(m), &mass)
	return int(status), m, float64(mass)
}

// SunNadir: the attitudes (9 each).
func SunNadir(gm float64, sun []float64, sunGM float64, t, y []float64) (int, []float64) {
	m := make([]float64, 9*len(t))
	status := C.cref_sun_nadir(C.double(gm), d(sun), C.double(sunGM), d(t), d(y), C.size_t(len(t)), d(m))
	return int(status), m
}

// LightTime: the light time and the range to a conic of mu from an observer.
func LightTime(c, mu float64, target, observer []float64, t float64) (int, float64, float64) {
	var lt, r C.double
	status := C.cref_light_time(C.double(c), C.double(mu), d(target), d(observer), C.double(t), &lt, &r)
	return int(status), float64(lt), float64(r)
}

// CKPointing: the quaternion of an instrument at the ticks.
func CKPointing(path string, instrument int, ticks float64) (int, []float64) {
	cpath := C.CString(path)
	defer C.free(unsafe.Pointer(cpath))
	q := make([]float64, 4)
	status := C.cref_ck_pointing(cpath, C.int(instrument), C.double(ticks), d(q))
	return int(status), q
}

// SCLKTime: the instant of encoded ticks of a clock.
func SCLKTime(kernel string, id int64, ticks float64) (int, int64, float64) {
	ck := C.CString(kernel)
	defer C.free(unsafe.Pointer(ck))
	var s C.longlong
	var f C.double
	status := C.cref_sclk_time(ck, C.longlong(id), C.double(ticks), &s, &f)
	return int(status), int64(s), float64(f)
}

// ShadowOrbit: the states (6 each) at the times of an orbit through the shadow of a planet
// with breaks, by a break scheme and a timescale, and the steps.
func ShadowOrbit(gm, radius, sunGM, sunRadius float64, sun []float64, sunMU, luminosity, area, cr float64,
	scheme, timescale int, y0 []float64, tEnd float64, times []float64) (int, []float64, int64) {
	states := make([]float64, 6*len(times))
	var steps C.longlong
	status := C.cref_shadow_orbit(C.double(gm), C.double(radius), C.double(sunGM),
		C.double(sunRadius), d(sun),
		C.double(sunMU), C.double(luminosity), C.double(area), C.double(cr), C.int(scheme),
		C.int(timescale), d(y0), C.double(tEnd), d(times), C.size_t(len(times)), d(states), &steps)
	return int(status), states, int64(steps)
}

func points(p []float64) []C.myb_contact_point {
	out := make([]C.myb_contact_point, len(p)/4)
	for k := range out {
		for i := range 3 {
			out[k].offset[i] = C.double(p[4*k+i])
		}
		out[k].radius = C.double(p[4*k+3])
	}
	return out
}

func pointsOf(p []C.myb_contact_point) *C.myb_contact_point {
	if len(p) == 0 {
		return nil
	}
	return &p[0]
}

// ContactSystem: the acceleration of the packed state (n), the states at the times (2 n each)
// and the deflections of the last; pts holds 4 values a point (offset, radius).
func ContactSystem(gm, point, normal, law, pts, bodies, deflections, times []float64) (int, []float64, []float64, []float64) {
	cp := points(pts)
	n := 6 + 3*len(cp)
	a := make([]float64, n)
	states := make([]float64, 2*n*len(times))
	final := make([]float64, 3*len(cp))
	var coordinates C.size_t
	status := C.cref_contact_system(d(gm), d(point), d(normal), d(law), pointsOf(cp),
		C.size_t(len(cp)), d(bodies), d(deflections), d(times), C.size_t(len(times)), d(a),
		d(states), d(final), &coordinates)
	if int(coordinates) != n {
		return -1, nil, nil, nil
	}
	return int(status), a, states, final
}

// ContactForces: 9 values a state (force, torque, energy, gap, touching), the rates of the
// deflections (with the friction that sticks), the state that failed.
func ContactForces(point, normal, law []float64, kind int, pts, x, v, attitude, omega, deflections []float64) (int, []float64, []float64, int) {
	cp := points(pts)
	count := len(x) / 3
	results := make([]float64, 9*count)
	rates := make([]float64, 3*len(cp)*count)
	var failed C.size_t
	status := C.cref_contact_forces(d(point), d(normal), d(law), C.int(kind), pointsOf(cp),
		C.size_t(len(cp)), d(x), d(v), d(attitude), d(omega), d(deflections), C.size_t(count),
		d(results), d(rates), &failed)
	return int(status), results, rates, int(failed)
}

// PNBinary: the accelerations (n each) of the states (12 values each) with the spins (3 each),
// and the state that failed.
func PNBinary(gm []float64, orders uint, c, r0, kappa1, kappa2 float64, bodies, chi1, chi2 []float64) (int, []float64, int) {
	count := len(bodies) / 12
	a := make([]float64, 12*count) // n = 6 + 6 with the spins
	var coordinates, failed C.size_t
	status := C.cref_pn_binary(d(gm), C.uint(orders), C.double(c), C.double(r0), C.double(kappa1),
		C.double(kappa2), d(bodies), d(chi1), d(chi2), C.size_t(count), d(a), &coordinates, &failed)
	if status == 0 && coordinates != 12 {
		return -1, nil, 0
	}
	return int(status), a, int(failed)
}
