// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include "glue.h"
*/
import "C"

import "unsafe"

// ---- Lambert's problem --------------------------------------------------------------------------

// LambertSolution is a solution of Lambert's problem (myb_lambert_solution).
type LambertSolution struct {
	V1, V2      [3]float64
	Revolutions int
	Branch      int
}

// LambertProblem is one problem: the two positions (m) and the time of flight (s).
type LambertProblem struct {
	R1, R2 [3]float64
	TOF    float64
}

func goLambert(s []C.myb_lambert_solution, n C.size_t) []LambertSolution {
	out := make([]LambertSolution, int(n))
	for k := range out {
		out[k] = LambertSolution{V1: goVec3(s[k].v1), V2: goVec3(s[k].v2),
			Revolutions: int(s[k].revolutions), Branch: int(s[k].branch)}
	}
	return out
}

// LambertMany solves the problems (myb_lambert), every solution with up to maxRevolutions
// whole turns, in one call.
func LambertMany(mu float64, problems []LambertProblem, retrograde bool, maxRevolutions int) ([][]LambertSolution, error) {
	const op = "LambertMany"
	if maxRevolutions < 0 || maxRevolutions > 1000 {
		return nil, refuse(op, -1, "revolutions out of range")
	}
	capacity := 1 + 2*maxRevolutions
	r1 := make([][3]float64, len(problems))
	r2 := make([][3]float64, len(problems))
	tof := make([]float64, len(problems))
	for k, p := range problems {
		r1[k], r2[k], tof[k] = p.R1, p.R2, p.TOF
	}
	solutions := make([]C.myb_lambert_solution, capacity*len(problems))
	counts := make([]C.size_t, len(problems))
	var failed C.size_t
	status := C.mybgo_lambert_many(C.double(mu), (*C.double)(first(r1)), (*C.double)(first(r2)),
		cDoubles(tof), C.size_t(len(problems)), flag(retrograde), C.int(maxRevolutions),
		C.size_t(capacity), (*C.myb_lambert_solution)(first(solutions)), (*C.size_t)(first(counts)),
		&failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	out := make([][]LambertSolution, len(problems))
	for k := range out {
		out[k] = goLambert(solutions[capacity*k:capacity*(k+1)], counts[k])
	}
	return out, nil
}

// LambertGrid is the grid of transfers between two bodies (a porkchop): the states of the
// departure body at the departures, of the arrival body at the arrivals, and for each pair
// with the arrival later the solutions of Lambert between the positions
// (Cells[i][j]: departure i, arrival j; empty where the arrival is not later).
type LambertGrid struct {
	Departures []State
	Arrivals   []State
	Cells      [][][]LambertSolution
}

// SolveLambertGrid computes the grid in one call.
func SolveLambertGrid(mu float64, departure, arrival StateSource, departures, arrivals []Time, retrograde bool, maxRevolutions int) (*LambertGrid, error) {
	const op = "SolveLambertGrid"
	if departure == nil || arrival == nil {
		return nil, refuse(op, -1, "no provider")
	}
	if maxRevolutions < 0 || maxRevolutions > 1000 {
		return nil, refuse(op, -1, "revolutions out of range")
	}
	if err := checkTimes(op, departures); err != nil {
		return nil, err
	}
	if err := checkTimes(op, arrivals); err != nil {
		return nil, err
	}
	k := newKeeper()
	defer k.done()
	var from, to C.mybgo_provider
	if err := departure.provider(&from, k); err != nil {
		return nil, err
	}
	if err := arrival.provider(&to, k); err != nil {
		return nil, err
	}
	capacity := 1 + 2*maxRevolutions
	cells := len(departures) * len(arrivals)
	out := &LambertGrid{Departures: make([]State, len(departures)), Arrivals: make([]State, len(arrivals))}
	solutions := make([]C.myb_lambert_solution, capacity*cells)
	counts := make([]C.size_t, cells)
	var failed C.size_t
	status := C.mybgo_lambert_grid(C.double(mu), &from, &to, cTimes(departures),
		C.size_t(len(departures)), cTimes(arrivals), C.size_t(len(arrivals)), flag(retrograde),
		C.int(maxRevolutions), C.size_t(capacity), cStates(out.Departures), cStates(out.Arrivals),
		(*C.myb_lambert_solution)(first(solutions)), (*C.size_t)(first(counts)), &failed)
	if err := k.failure(op, status, int(failed)); err != nil {
		return nil, err
	}
	out.Cells = make([][][]LambertSolution, len(departures))
	for i := range out.Cells {
		out.Cells[i] = make([][]LambertSolution, len(arrivals))
		for j := range out.Cells[i] {
			c := i*len(arrivals) + j
			out.Cells[i][j] = goLambert(solutions[capacity*c:capacity*(c+1)], counts[c])
		}
	}
	return out, nil
}

// ---- B-plane and flybys -------------------------------------------------------------------------

// BPlane is the B-plane of a hyperbola (myb_bplane).
type BPlane struct {
	VInfinity, C3      float64
	S, T, R, B         [3]float64
	BT, BR, BMagnitude float64
	Theta              float64
	Periapsis          float64
	Eccentricity       float64
	TimeToPeriapsis    float64
}

// BPlanes are the B-planes of the states relative to the target of gravitational parameter
// gm with the pole given, in one call.
func BPlanes(gm float64, states []State, pole [3]float64) ([]BPlane, error) {
	c := make([]C.myb_bplane, len(states))
	var failed C.size_t
	status := C.mybgo_bplanes(C.double(gm), cStates(states), C.size_t(len(states)),
		(*C.double)(unsafe.Pointer(&pole[0])), (*C.myb_bplane)(first(c)), &failed)
	if err := failAt("BPlanes", status, failed); err != nil {
		return nil, err
	}
	out := make([]BPlane, len(c))
	for k, b := range c {
		out[k] = BPlane{VInfinity: float64(b.v_infinity), C3: float64(b.c3), S: goVec3(b.s),
			T: goVec3(b.t), R: goVec3(b.r), B: goVec3(b.b), BT: float64(b.b_t), BR: float64(b.b_r),
			BMagnitude: float64(b.b_magnitude), Theta: float64(b.theta),
			Periapsis: float64(b.periapsis), Eccentricity: float64(b.eccentricity),
			TimeToPeriapsis: float64(b.time_to_periapsis)}
	}
	return out, nil
}

// SphereOfInfluence is a (m/M)^(2/5) (myb_sphere_of_influence).
func SphereOfInfluence(a, massRatio float64) float64 {
	return float64(C.myb_sphere_of_influence(C.double(a), C.double(massRatio)))
}

// HillRadius is a (1 - e) (m / 3M)^(1/3) (myb_hill_radius).
func HillRadius(a, e, massRatio float64) float64 {
	return float64(C.myb_hill_radius(C.double(a), C.double(e), C.double(massRatio)))
}

// FlybyTurn is the turn of the excess velocity of a flyby (myb_flyby_turn).
func FlybyTurn(gm, vInfinity, periapsis float64) float64 {
	return float64(C.myb_flyby_turn(C.double(gm), C.double(vInfinity), C.double(periapsis)))
}

// DepartureImpulse is the impulse from a circle of radius r to the excess speed vInfinity
// (myb_departure_impulse).
func DepartureImpulse(gm, radius, vInfinity float64) float64 {
	return float64(C.myb_departure_impulse(C.double(gm), C.double(radius), C.double(vInfinity)))
}

// ---- Shooting -----------------------------------------------------------------------------------

// ShootFunc fills f (m values) at the unknowns u (n values) and, when jacobian is not nil,
// the m x n Jacobian (row-major). The slices are valid during the call only.
type ShootFunc func(u, f, jacobian []float64) error

// Shoot solves f(u) = 0 by Newton's method (myb_shoot) with a Go residual function — the slow
// path of callbacks: u (n values) is updated in place; m equations; the iterations and the
// residual max |f_i|.
func Shoot(f ShootFunc, u []float64, m int, tolerance float64, maxIterations int) (int, float64, error) {
	const op = "Shoot"
	n := len(u)
	if f == nil || n < 1 || n > 64 || m < 1 || m > 64 {
		return 0, 0, refuse(op, -1, "a function and 1 .. 64 unknowns and equations")
	}
	k := newKeeper()
	defer k.done()
	h := k.handle(&shootCallback{f: f, n: n, m: m, state: k.state})
	options := C.myb_shoot_options{tolerance: C.double(tolerance), max_iterations: C.int(maxIterations)}
	var iterations C.int
	var residual C.double
	status := C.mybgo_shoot(h, C.size_t(n), C.size_t(m), cDoubles(u), &options, &iterations, &residual)
	return int(iterations), float64(residual), k.failure(op, status, -1)
}

// ---- Initial orbit determination ----------------------------------------------------------------

// IODSolution is an orbit from directions (myb_iod_solution).
type IODSolution struct {
	R, V       [3]float64
	Range      [3]float64
	Iterations int
}

func goIOD(c *C.myb_iod_solution) IODSolution {
	return IODSolution{R: goVec3(c.r), V: goVec3(c.v), Range: goVec3(c._range),
		Iterations: int(c.iterations)}
}

// IODGibbs is the velocity at the second of three positions of one orbit (myb_iod_gibbs) and
// the angle of the first out of the plane of the others.
func IODGibbs(mu float64, r1, r2, r3 [3]float64) ([3]float64, float64, error) {
	var v [3]float64
	var angle C.double
	status := C.myb_iod_gibbs(C.double(mu), (*C.double)(unsafe.Pointer(&r1[0])),
		(*C.double)(unsafe.Pointer(&r2[0])), (*C.double)(unsafe.Pointer(&r3[0])),
		(*C.double)(unsafe.Pointer(&v[0])), &angle)
	return v, float64(angle), fail("IODGibbs", status)
}

// IODHerrickGibbs is the velocity at the second of three close positions (myb_iod_herrick_gibbs).
func IODHerrickGibbs(mu float64, r1 [3]float64, t1 float64, r2 [3]float64, t2 float64, r3 [3]float64, t3 float64) ([3]float64, error) {
	var v [3]float64
	status := C.myb_iod_herrick_gibbs(C.double(mu), (*C.double)(unsafe.Pointer(&r1[0])), C.double(t1),
		(*C.double)(unsafe.Pointer(&r2[0])), C.double(t2), (*C.double)(unsafe.Pointer(&r3[0])),
		C.double(t3), (*C.double)(unsafe.Pointer(&v[0])))
	return v, fail("IODHerrickGibbs", status)
}

// IODSights are three directions (unit vectors) from three observers at three times.
type IODSights struct {
	T        [3]float64
	Observer [3][3]float64
	Sight    [3][3]float64
}

func (s *IODSights) c() (*C.double, *[3]C.double, *[3]C.double) {
	return (*C.double)(unsafe.Pointer(&s.T[0])), (*[3]C.double)(unsafe.Pointer(&s.Observer[0])),
		(*[3]C.double)(unsafe.Pointer(&s.Sight[0]))
}

func iodMany(op string, status C.myb_status, c []C.myb_iod_solution, n C.size_t) ([]IODSolution, error) {
	if err := fail(op, status); err != nil {
		return nil, err
	}
	out := make([]IODSolution, int(n))
	for k := range out {
		out[k] = goIOD(&c[k])
	}
	return out, nil
}

// IODGauss is Gauss's method with the refinement (myb_iod_gauss): every orbit found.
func IODGauss(mu float64, s IODSights) ([]IODSolution, error) {
	t, observer, sight := s.c()
	c := make([]C.myb_iod_solution, C.MYB_IOD_MAX_SOLUTIONS)
	var n C.size_t
	status := C.myb_iod_gauss(C.double(mu), t, observer, sight, &c[0], &n)
	return iodMany("IODGauss", status, c, n)
}

// IODLaplace is Laplace's method with the refinement (myb_iod_laplace).
func IODLaplace(mu float64, s IODSights) ([]IODSolution, error) {
	t, observer, sight := s.c()
	c := make([]C.myb_iod_solution, C.MYB_IOD_MAX_SOLUTIONS)
	var n C.size_t
	status := C.myb_iod_laplace(C.double(mu), t, observer, sight, &c[0], &n)
	return iodMany("IODLaplace", status, c, n)
}

// IODDoubleR is the double-r iteration from guesses of the distances at the first two times
// (myb_iod_double_r).
func IODDoubleR(mu float64, s IODSights, distance [2]float64) (IODSolution, error) {
	t, observer, sight := s.c()
	var c C.myb_iod_solution
	status := C.myb_iod_double_r(C.double(mu), t, observer, sight,
		(*C.double)(unsafe.Pointer(&distance[0])), &c)
	return goIOD(&c), fail("IODDoubleR", status)
}

// IODGooding is Gooding's method from guesses of the first and third ranges (myb_iod_gooding).
func IODGooding(mu float64, s IODSights, ranges [2]float64, revolutions, branch int) (IODSolution, error) {
	t, observer, sight := s.c()
	var c C.myb_iod_solution
	status := C.myb_iod_gooding(C.double(mu), t, observer, sight,
		(*C.double)(unsafe.Pointer(&ranges[0])), C.int(revolutions), C.int(branch), &c)
	return goIOD(&c), fail("IODGooding", status)
}

// ---- Relative motion ----------------------------------------------------------------------------

func relative(op string, kind C.int, mu float64, chief, other []State) ([]State, error) {
	if len(chief) != len(other) {
		return nil, refuse(op, -1, "slices of different lengths")
	}
	out := make([]State, len(chief))
	var failed C.size_t
	status := C.mybgo_relative(kind, C.double(mu), cStates(chief), cStates(other),
		C.size_t(len(chief)), cStates(out), &failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// RelativeFromInertial are the relative states of the deputies in the turning frames of the
// chiefs (myb_relative_from_inertial), in one call.
func RelativeFromInertial(chief, deputy []State) ([]State, error) {
	return relative("RelativeFromInertial", C.MYBGO_RELATIVE_FROM_INERTIAL, 0, chief, deputy)
}

// RelativeToInertial are the inertial states of the deputies (myb_relative_to_inertial).
func RelativeToInertial(chief, relativeStates []State) ([]State, error) {
	return relative("RelativeToInertial", C.MYBGO_RELATIVE_TO_INERTIAL, 0, chief, relativeStates)
}

// ROEFromStates are the relative orbital elements of the deputies (myb_roe_from_states).
func ROEFromStates(mu float64, chief, deputy []State) ([][6]float64, error) {
	s, err := relative("ROEFromStates", C.MYBGO_RELATIVE_ROE, mu, chief, deputy)
	if err != nil {
		return nil, err
	}
	out := make([][6]float64, len(s))
	for k := range s {
		out[k] = s[k]
	}
	return out, nil
}

// ROEPositions are the relative positions the elements give at the mean arguments of latitude
// u of a chief of semi-major axis a (myb_roe_to_position), in one call.
func ROEPositions(a float64, u []float64, roe [6]float64) [][3]float64 {
	out := make([][3]float64, len(u))
	C.mybgo_roe_positions(C.double(a), cDoubles(u), C.size_t(len(u)),
		(*C.double)(unsafe.Pointer(&roe[0])), (*C.double)(first(out)))
	return out
}

// RelativeSTM are the transition matrices of the linear relative motion about the chief (a
// two-body orbit of mu) for the offsets dt (myb_relative_stm), in one call.
func RelativeSTM(mu float64, chief State, dt []float64) ([][36]float64, error) {
	out := make([][36]float64, len(dt))
	var failed C.size_t
	status := C.mybgo_relative_stm(C.double(mu), (*C.double)(unsafe.Pointer(&chief)), 0,
		cDoubles(dt), C.size_t(len(dt)), (*C.double)(first(out)), &failed)
	if err := failAt("RelativeSTM", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// CWSTM are the matrices of Clohessy and Wiltshire for the mean motion n (myb_cw_stm).
func CWSTM(n float64, dt []float64) ([][36]float64, error) {
	out := make([][36]float64, len(dt))
	var failed C.size_t
	status := C.mybgo_relative_stm(0, nil, C.double(n), cDoubles(dt), C.size_t(len(dt)),
		(*C.double)(first(out)), &failed)
	if err := failAt("CWSTM", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// Rendezvous are the two impulses of a rendezvous in the linear motion
// (myb_relative_rendezvous); targetVelocity nil: at rest with the chief.
func Rendezvous(mu float64, chief State, dt float64, relativeState State, target [3]float64, targetVelocity *[3]float64) ([3]float64, [3]float64, error) {
	var c C.myb_pv
	for i := range 3 {
		c.r[i], c.v[i] = C.double(chief[i]), C.double(chief[3+i])
	}
	var dv0, dv1 [3]float64
	var tv *C.double
	if targetVelocity != nil {
		v := *targetVelocity
		tv = (*C.double)(unsafe.Pointer(&v[0]))
	}
	status := C.myb_relative_rendezvous(C.double(mu), &c, C.double(dt),
		(*C.double)(unsafe.Pointer(&relativeState)), (*C.double)(unsafe.Pointer(&target[0])), tv,
		(*C.double)(unsafe.Pointer(&dv0[0])), (*C.double)(unsafe.Pointer(&dv1[0])))
	return dv0, dv1, fail("Rendezvous", status)
}
