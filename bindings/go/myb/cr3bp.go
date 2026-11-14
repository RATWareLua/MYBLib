// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include "glue.h"
*/
import "C"

import "unsafe"

// The restricted problem of three bodies (myb/cr3bp.h): lengths in the distance of the
// primaries, times in 1 / their mean motion, mu = m2 / (m1 + m2).

// CR3BPLagrange are the five Lagrange points (myb_cr3bp_lagrange).
func CR3BPLagrange(mu float64) ([5][3]float64, error) {
	var points [5][3]float64
	status := C.myb_cr3bp_lagrange(C.double(mu), (*[3]C.double)(unsafe.Pointer(&points[0])))
	return points, fail("CR3BPLagrange", status)
}

// CR3BPJacobi are the Jacobi constants of the states, in one call.
func CR3BPJacobi(mu float64, y []State) []float64 {
	out := make([]float64, len(y))
	C.mybgo_cr3bp_jacobi(C.double(mu), cStates(y), C.size_t(len(y)), cDoubles(out))
	return out
}

// CR3BPSymmetry is the symmetry of a periodic orbit (myb_cr3bp_symmetry).
type CR3BPSymmetry int

// The symmetries the single shooting corrects.
const (
	SymmetryXZ CR3BPSymmetry = C.MYB_CR3BP_SYMMETRY_XZ
	SymmetryX  CR3BPSymmetry = C.MYB_CR3BP_SYMMETRY_X
)

// CR3BPCorrection are the settings of the single shooting (myb_cr3bp_correction).
type CR3BPCorrection struct {
	Symmetry      CR3BPSymmetry
	Fixed         int // the component of y0 held fixed
	Crossing      int // the half period ends at this crossing of y = 0 (1 = the first)
	Tolerance     float64
	MaxIterations int
}

// CR3BPRichardson is the third-order halo orbit of Richardson about L1 (1) or L2 (2): the
// state at the crossing of the xz-plane and the period.
func CR3BPRichardson(mu float64, point int, az float64, branch int) (State, float64, error) {
	var y State
	var period C.double
	status := C.myb_cr3bp_richardson(C.double(mu), C.int(point), C.double(az), C.int(branch),
		(*C.double)(unsafe.Pointer(&y)), &period)
	return y, float64(period), fail("CR3BPRichardson", status)
}

// CR3BPCorrect corrects y0 to a symmetric periodic orbit and gives its period.
func CR3BPCorrect(mu float64, c CR3BPCorrection, y0 State) (State, float64, error) {
	cc := C.myb_cr3bp_correction{symmetry: C.myb_cr3bp_symmetry(c.Symmetry), fixed: C.int(c.Fixed),
		crossing: C.int(c.Crossing), tolerance: C.double(c.Tolerance),
		max_iterations: C.int(c.MaxIterations)}
	var period C.double
	status := C.myb_cr3bp_correct(C.double(mu), &cc, (*C.double)(unsafe.Pointer(&y0)), &period)
	return y0, float64(period), fail("CR3BPCorrect", status)
}

// CR3BPMonodromy is the monodromy matrix of the periodic orbit (row-major).
func CR3BPMonodromy(mu float64, y0 State, period float64) ([36]float64, error) {
	var m [36]float64
	status := C.myb_cr3bp_monodromy(C.double(mu), (*C.double)(unsafe.Pointer(&y0)), C.double(period),
		(*C.double)(unsafe.Pointer(&m[0])))
	return m, fail("CR3BPMonodromy", status)
}

// CR3BPStability are the eigenvalues of a monodromy matrix and the stability index.
func CR3BPStability(monodromy [36]float64) (re, im [6]float64, index float64, err error) {
	var ci C.double
	status := C.myb_cr3bp_stability((*C.double)(unsafe.Pointer(&monodromy[0])),
		(*C.double)(unsafe.Pointer(&re[0])), (*C.double)(unsafe.Pointer(&im[0])), &ci)
	return re, im, float64(ci), fail("CR3BPStability", status)
}

// CR3BPManifold are states on the unstable (or stable) manifold of the periodic orbit at the
// phases (fractions of the period), displaced by epsilon on the side sign, in one call.
func CR3BPManifold(mu float64, y0 State, period float64, phases []float64, unstable bool, sign int, epsilon float64) ([]State, error) {
	out := make([]State, len(phases))
	var failed C.size_t
	status := C.mybgo_cr3bp_manifold(C.double(mu), (*C.double)(unsafe.Pointer(&y0)),
		C.double(period), cDoubles(phases), C.size_t(len(phases)), flag(unstable), C.int(sign),
		C.double(epsilon), cStates(out), &failed)
	if err := failAt("CR3BPManifold", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// CR3BPOrbit is a periodic orbit for multiple shooting (myb_cr3bp_orbit): the states at the
// cuts k period / patches, the period, and the tangent of its family (nil until corrected).
type CR3BPOrbit struct {
	States  []State
	Period  float64
	Tangent []float64
}

func (o CR3BPOrbit) c() (C.myb_cr3bp_orbit, error) {
	var out C.myb_cr3bp_orbit
	if len(o.States) == 0 || len(o.States) > C.MYB_CR3BP_MAX_PATCHES ||
		(o.Tangent != nil && len(o.Tangent) != 6*len(o.States)+1) {
		return out, refuse("CR3BPOrbit", -1, "1 .. 10 patches, a tangent of 6 patches + 1 values")
	}
	out.patches = C.size_t(len(o.States))
	for k, s := range o.States {
		for i := range 6 {
			out.states[6*k+i] = C.double(s[i])
		}
	}
	out.period = C.double(o.Period)
	if o.Tangent != nil {
		out.has_tangent = 1
		for i, v := range o.Tangent {
			out.tangent[i] = C.double(v)
		}
	}
	return out, nil
}

func goOrbit(c *C.myb_cr3bp_orbit) CR3BPOrbit {
	n := int(c.patches)
	out := CR3BPOrbit{States: make([]State, n), Period: float64(c.period)}
	for k := range n {
		for i := range 6 {
			out.States[k][i] = float64(c.states[6*k+i])
		}
	}
	if c.has_tangent != 0 {
		out.Tangent = make([]float64, 6*n+1)
		for i := range out.Tangent {
			out.Tangent[i] = float64(c.tangent[i])
		}
	}
	return out
}

// CR3BPOrbitFromState cuts the orbit through y0 with the period into patches.
func CR3BPOrbitFromState(mu float64, y0 State, period float64, patches int) (CR3BPOrbit, error) {
	if patches < 1 || patches > C.MYB_CR3BP_MAX_PATCHES {
		return CR3BPOrbit{}, refuse("CR3BPOrbitFromState", -1, "1 .. 10 patches")
	}
	var c C.myb_cr3bp_orbit
	status := C.myb_cr3bp_orbit_from_state(C.double(mu), (*C.double)(unsafe.Pointer(&y0)),
		C.double(period), C.size_t(patches), &c)
	if err := fail("CR3BPOrbitFromState", status); err != nil {
		return CR3BPOrbit{}, err
	}
	return goOrbit(&c), nil
}

// CR3BPFamily corrects the orbit and follows its family count - 1 steps of `step` along the
// tangent (pseudo-arclength), in one call: the members found. On an error the members found
// before it come back with it (Index: the member that failed; take a shorter step there).
func CR3BPFamily(mu float64, orbit CR3BPOrbit, step float64, count int, tolerance float64, maxIterations int) ([]CR3BPOrbit, error) {
	const op = "CR3BPFamily"
	c, err := orbit.c()
	if err != nil {
		return nil, err
	}
	if count < 1 {
		return nil, refuse(op, -1, "no member asked for")
	}
	members := make([]C.myb_cr3bp_orbit, count)
	var made C.size_t
	status := C.mybgo_cr3bp_family(C.double(mu), &c, C.double(step), C.size_t(count),
		C.double(tolerance), C.int(maxIterations), &members[0], &made)
	out := make([]CR3BPOrbit, int(made))
	for k := range out {
		out[k] = goOrbit(&members[k])
	}
	return out, failAt(op, status, made)
}

// CR3BPTorus is an invariant circle of a torus (myb_cr3bp_torus): its points, the time T of a
// turn, the rotation rho across it, the defect, and the tangent of its family (nil without).
type CR3BPTorus struct {
	States   []State
	Period   float64
	Rotation float64
	Defect   float64
	Tangent  []float64
}

func (t CR3BPTorus) c() (C.myb_cr3bp_torus, error) {
	var out C.myb_cr3bp_torus
	n := len(t.States)
	if n < 5 || n > C.MYB_CR3BP_MAX_CIRCLE || n%2 == 0 || (t.Tangent != nil && len(t.Tangent) != 6*n+2) {
		return out, refuse("CR3BPTorus", -1, "an odd number of points 5 .. 35, a tangent of 6 points + 2")
	}
	out.points = C.size_t(n)
	for k, s := range t.States {
		for i := range 6 {
			out.states[6*k+i] = C.double(s[i])
		}
	}
	out.period = C.double(t.Period)
	out.rotation = C.double(t.Rotation)
	out.defect = C.double(t.Defect)
	if t.Tangent != nil {
		out.has_tangent = 1
		for i, v := range t.Tangent {
			out.tangent[i] = C.double(v)
		}
	}
	return out, nil
}

func goTorus(c *C.myb_cr3bp_torus) CR3BPTorus {
	n := int(c.points)
	out := CR3BPTorus{States: make([]State, n), Period: float64(c.period),
		Rotation: float64(c.rotation), Defect: float64(c.defect)}
	for k := range n {
		for i := range 6 {
			out.States[k][i] = float64(c.states[6*k+i])
		}
	}
	if c.has_tangent != 0 {
		out.Tangent = make([]float64, 6*n+2)
		for i := range out.Tangent {
			out.Tangent[i] = float64(c.tangent[i])
		}
	}
	return out
}

// TorusFamily is the value a family of tori keeps (myb_cr3bp_torus_family).
type TorusFamily int

// The families of tori.
const (
	TorusJacobi   TorusFamily = C.MYB_CR3BP_TORUS_JACOBI
	TorusPeriod   TorusFamily = C.MYB_CR3BP_TORUS_PERIOD
	TorusRotation TorusFamily = C.MYB_CR3BP_TORUS_ROTATION
)

// CR3BPTorusFromOrbit is the smallest torus about the periodic orbit (myb_cr3bp_torus_from_orbit).
func CR3BPTorusFromOrbit(mu float64, y0 State, period float64, points, pair int, size float64) (CR3BPTorus, error) {
	var c C.myb_cr3bp_torus
	status := C.myb_cr3bp_torus_from_orbit(C.double(mu), (*C.double)(unsafe.Pointer(&y0)),
		C.double(period), C.size_t(points), C.int(pair), C.double(size), &c)
	if err := fail("CR3BPTorusFromOrbit", status); err != nil {
		return CR3BPTorus{}, err
	}
	return goTorus(&c), nil
}

// CR3BPTorusFamily corrects the torus in the family (a step of 0) and follows it count - 1
// steps, in one call: the members found (with the error, as CR3BPFamily).
func CR3BPTorusFamily(mu float64, torus CR3BPTorus, family TorusFamily, step float64, count int, tolerance float64, maxIterations int) ([]CR3BPTorus, error) {
	const op = "CR3BPTorusFamily"
	c, err := torus.c()
	if err != nil {
		return nil, err
	}
	if count < 1 {
		return nil, refuse(op, -1, "no member asked for")
	}
	members := make([]C.myb_cr3bp_torus, count)
	var made C.size_t
	status := C.mybgo_cr3bp_torus_family(C.double(mu), &c, C.myb_cr3bp_torus_family(family),
		C.double(step), C.size_t(count), C.double(tolerance), C.int(maxIterations), &members[0], &made)
	out := make([]CR3BPTorus, int(made))
	for k := range out {
		out[k] = goTorus(&members[k])
	}
	return out, failAt(op, status, made)
}

// StatesAt are the points of the circle at the angles (rad), by its trigonometric
// interpolant, in one call.
func (t CR3BPTorus) StatesAt(angles []float64) ([]State, error) {
	c, err := t.c()
	if err != nil {
		return nil, err
	}
	out := make([]State, len(angles))
	C.mybgo_cr3bp_torus_states(&c, cDoubles(angles), C.size_t(len(angles)), cStates(out))
	return out, nil
}

// TorusSpectrum is the stability of a torus (myb_cr3bp_torus_spectrum): the eigenvalues of the
// map, their classes, the hyperbolic pair with its eigenfunctions, and the normal rotation of
// an elliptic torus with its complex eigenfunction.
type TorusSpectrum struct {
	Re, Im            []float64
	Trivial           int
	Elliptic          int
	Hyperbolic        int
	Unresolved        int
	Index             float64
	Unstable, Stable  float64
	UnstableDirection []State
	StableDirection   []State
	Residual, Tail    float64
	NormalRotation    float64 // sigma (rad); 0: no elliptic pair
	NormalEigenvalue  [2]float64
	NormalFrequency   float64
	NormalResidual    float64
	NormalRe          []State // the eigenfunction psi = NormalRe + i NormalIm at the points
	NormalIm          []State
	c                 C.myb_cr3bp_torus_spectrum
}

// CR3BPTorusStability is the spectrum of the torus.
func CR3BPTorusStability(mu float64, torus CR3BPTorus, tolerance float64) (*TorusSpectrum, error) {
	c, err := torus.c()
	if err != nil {
		return nil, err
	}
	out := &TorusSpectrum{}
	status := C.myb_cr3bp_torus_stability(C.double(mu), &c, C.double(tolerance), &out.c)
	if err := fail("CR3BPTorusStability", status); err != nil {
		return nil, err
	}
	s := &out.c
	count := int(s.count)
	points := count / 6
	out.Re, out.Im = make([]float64, count), make([]float64, count)
	for i := range count {
		out.Re[i], out.Im[i] = float64(s.re[i]), float64(s.im[i])
	}
	out.Trivial, out.Elliptic = int(s.trivial), int(s.elliptic)
	out.Hyperbolic, out.Unresolved = int(s.hyperbolic), int(s.unresolved)
	out.Index, out.Unstable, out.Stable = float64(s.index), float64(s.unstable), float64(s.stable)
	out.Residual, out.Tail = float64(s.residual), float64(s.tail)
	out.NormalRotation = float64(s.normal_rotation)
	out.NormalEigenvalue = [2]float64{float64(s.normal_eigenvalue[0]), float64(s.normal_eigenvalue[1])}
	out.NormalFrequency, out.NormalResidual = float64(s.normal_frequency), float64(s.normal_residual)
	grid := func(v *[C.MYB_CR3BP_MAX_SPECTRUM]C.double) []State {
		states := make([]State, points)
		for k := range points {
			for i := range 6 {
				states[k][i] = float64(v[6*k+i])
			}
		}
		return states
	}
	out.UnstableDirection = grid(&s.unstable_direction)
	out.StableDirection = grid(&s.stable_direction)
	out.NormalRe = grid(&s.normal_re)
	out.NormalIm = grid(&s.normal_im)
	return out, nil
}

// CR3BPTorusManifold are states on the unstable (or stable) manifold of the torus at the
// angles, displaced by epsilon on the side sign, in one call.
func CR3BPTorusManifold(torus CR3BPTorus, spectrum *TorusSpectrum, angles []float64, unstable bool, sign int, epsilon float64) ([]State, error) {
	const op = "CR3BPTorusManifold"
	c, err := torus.c()
	if err != nil {
		return nil, err
	}
	if spectrum == nil {
		return nil, refuse(op, -1, "no spectrum")
	}
	out := make([]State, len(angles))
	var failed C.size_t
	status := C.mybgo_cr3bp_torus_manifold(&c, &spectrum.c, cDoubles(angles), C.size_t(len(angles)),
		flag(unstable), C.int(sign), C.double(epsilon), cStates(out), &failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// CR3BPSection is a plane of section: the component Axis (0 .. 5) of the state equal to
// Value, crossed in the Direction (+1, -1, 0).
type CR3BPSection struct {
	Axis      int
	Value     float64
	Direction int
}

// CR3BPCrossing is where a propagation ended: the state, the time from the start and whether
// the section was reached.
type CR3BPCrossing struct {
	State   State
	Elapsed float64
	Crossed bool
}

// CR3BPPropagate carries every state for the time duration, stopping at the first crossing of
// the section when one is given (myb_cr3bp_propagate), in one call.
func CR3BPPropagate(mu float64, y0 []State, duration float64, section *CR3BPSection) ([]CR3BPCrossing, error) {
	var cs *C.myb_cr3bp_section
	if section != nil {
		cs = &C.myb_cr3bp_section{axis: C.int(section.Axis), value: C.double(section.Value),
			direction: C.int(section.Direction)}
	}
	y := make([]State, len(y0))
	elapsed := make([]float64, len(y0))
	crossed := make([]C.int, len(y0))
	var failed C.size_t
	status := C.mybgo_cr3bp_propagate(C.double(mu), cStates(y0), C.size_t(len(y0)), C.double(duration),
		cs, cStates(y), cDoubles(elapsed), (*C.int)(first(crossed)), &failed)
	if err := failAt("CR3BPPropagate", status, failed); err != nil {
		return nil, err
	}
	out := make([]CR3BPCrossing, len(y0))
	for k := range out {
		out[k] = CR3BPCrossing{State: y[k], Elapsed: elapsed[k], Crossed: crossed[k] != 0}
	}
	return out, nil
}

// CR3BPFromEphemeris turns states (m, m/s) into the pulsating rotating frame of two bodies at
// each instant (their states body1[k], body2[k]), in one call; CR3BPToEphemeris back.
func CR3BPFromEphemeris(body1, body2 []State, gm1, gm2 float64, states []State) ([]State, error) {
	return cr3bpFrames("CR3BPFromEphemeris", 0, body1, body2, gm1, gm2, states)
}

// CR3BPToEphemeris turns rotating states back into states of the frame of the bodies.
func CR3BPToEphemeris(body1, body2 []State, gm1, gm2 float64, rotating []State) ([]State, error) {
	return cr3bpFrames("CR3BPToEphemeris", 1, body1, body2, gm1, gm2, rotating)
}

func cr3bpFrames(op string, to C.int, body1, body2 []State, gm1, gm2 float64, in []State) ([]State, error) {
	if len(body1) != len(in) || len(body2) != len(in) {
		return nil, refuse(op, -1, "slices of different lengths")
	}
	out := make([]State, len(in))
	var failed C.size_t
	status := C.mybgo_cr3bp_frames(to, cStates(body1), cStates(body2), C.double(gm1), C.double(gm2),
		cStates(in), C.size_t(len(in)), cStates(out), &failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	return out, nil
}
