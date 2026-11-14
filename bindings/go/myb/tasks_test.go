// SPDX-License-Identifier: BSD-3-Clause

// The coarse tasks of the binding against the same tasks written in C (package internal/cref:
// the functions of the library called as a C program calls them), bit for bit: systems of
// bodies with a merger, their invariants and derivatives, Wisdom-Holman runs; the filter and
// the smoother, measurements over slices; families of periodic orbits and of tori with the
// normal rotation; Lambert over a grid of dates, conjunctions, the Q-law, chains with a
// powered flyby, Sims and Flanagan; attitude laws and the light time over slices; CK and the
// clocks of SCLK. Messages (OPM, OMM, SP3, CK) written and read back, and the errors of bad
// input.
package myb_test

import (
	"errors"
	"math"
	"testing"

	"github.com/RATWareLua/MYBLib/bindings/go/myb"
	"github.com/RATWareLua/MYBLib/bindings/go/myb/internal/cref"
)

func sameFloats(t *testing.T, what string, got, want []float64) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("%s: %d values, the C calls %d", what, len(got), len(want))
	}
	for k := range got {
		if math.Float64bits(got[k]) != math.Float64bits(want[k]) {
			t.Errorf("%s[%d] = %v, the C calls %v", what, k, got[k], want[k])
			return
		}
	}
}

func flatten[T ~[3]float64 | ~[4]float64 | ~[6]float64 | ~[9]float64 | ~[36]float64](s []T) []float64 {
	var out []float64
	for _, v := range s {
		for i := range len(v) {
			out = append(out, v[i])
		}
	}
	return out
}

func okStatus(t *testing.T, what string, status int) {
	t.Helper()
	if status != 0 {
		t.Fatalf("%s: the C reference failed with status %d", what, status)
	}
}

// ---- Systems of bodies --------------------------------------------------------------------

func TestSystemRunWithMerger(t *testing.T) {
	gms := []float64{gm, 4.9e12, 1e9}
	radii := []float64{6.4e6, 1.7e6, 0}
	bodies := []myb.State{
		{0, 0, 0, 0, 0, 0},
		{1.2e7, 1e5, 0, -3000, 200, 0},
		{0, 4.2e7, 0, -3075, 0, 0},
	}
	spec := myb.SystemSpec{}
	for i := range gms {
		spec.Bodies = append(spec.Bodies, myb.SystemBody{GM: gms[i], Radius: radii[i]})
	}
	system, err := myb.NewSystem(spec)
	must(t, err)
	defer system.Close()
	y0, err := system.PackState(bodies, nil, nil)
	must(t, err)
	times := []float64{100, 400, 900, 1600, 2500}
	run, err := system.Run(myb.SystemRun{T0: 0, Y0: y0, Times: times, OnContact: myb.MergeAtContact})
	must(t, err)
	status, states, events, final := cref.SystemRun(gms, radii, y0, times)
	okStatus(t, "SystemRun", status)
	var flat []float64
	for _, s := range run.States {
		flat = append(flat, s...)
	}
	sameFloats(t, "states", flat, states)
	sameFloats(t, "final", run.Final, final)
	if len(run.Events) != 1 || len(events) != 7 {
		t.Fatalf("events: %d, the C calls %d", len(run.Events), len(events)/7)
	}
	e := run.Events[0]
	sameFloats(t, "event", []float64{float64(e.Kind), float64(e.A), float64(e.B), e.T, e.Distance,
		float64(e.Survivor), e.Merger.Energy}, events)
	if e.Kind != myb.Contact || e.Action != myb.MergeAtContact || e.Survivor != 0 {
		t.Errorf("event %+v", e)
	}
	if present, _, err := system.Present(1); err != nil || present {
		t.Errorf("the merged body is still there: %v", err)
	}
	if _, after, _ := system.Present(0); after == gms[0] {
		t.Errorf("the survivor did not take the mass")
	}

	// the invariants and the derivatives on a system of its own
	fresh, err := myb.NewSystem(spec)
	must(t, err)
	defer fresh.Close()
	sample := [][]float64{y0, run.States[0]}
	inv, err := fresh.Invariants(sample)
	must(t, err)
	status, energy, momentum := cref.Invariants(gms, append(append([]float64{}, y0...), run.States[0]...), 2)
	okStatus(t, "Invariants", status)
	sameFloats(t, "energy", []float64{inv[0].Energy.Hi, inv[1].Energy.Hi}, energy)
	sameFloats(t, "momentum", append(inv[0].Momentum[:], inv[1].Momentum[:]...), momentum)
	partials, err := fresh.Partials(myb.IAS15Options{}, 0, y0, 50, []int{0})
	must(t, err)
	status, y1, phi, sens := cref.SystemPartials(gms, y0, 50)
	okStatus(t, "SystemPartials", status)
	sameFloats(t, "y1", partials.Y1, y1)
	sameFloats(t, "phi", partials.Phi, phi)
	sameFloats(t, "sensitivity", partials.Sensitivity, sens)
	accelerations, err := fresh.Accelerations([]float64{0}, [][]float64{y0})
	must(t, err)
	if len(accelerations[0]) != 9 {
		t.Errorf("accelerations of %d coordinates", len(accelerations[0]))
	}
	back, _, _, err := fresh.UnpackState(y0, nil, nil)
	must(t, err)
	if back[1] != bodies[1] {
		t.Errorf("the state unpacked is not the one packed")
	}

	// errors: a link to a body that is not there, a state of another length
	_, err = myb.NewSystem(myb.SystemSpec{Bodies: spec.Bodies, Links: []myb.Link{myb.FieldLink{Body: 0, Point: 9}}})
	var me *myb.Error
	if !errors.As(err, &me) || me.Status != myb.Domain || me.Message != "link" || me.Index != 0 {
		t.Errorf("a link to no body: %v", err)
	}
	if _, err := fresh.Run(myb.SystemRun{Y0: y0[:5], Times: times}); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a short state: %v", err)
	}
}

func TestWisdomHolman(t *testing.T) {
	gm := []float64{1.32712440018e20, 1.26686534e17}
	x := [][3]float64{{0, 0, 0}, {7.78e11, 0, 0}}
	v := [][3]float64{{0, 0, 0}, {0, 13060, 100}}
	times := []float64{1e7, 3e7, 1e8}
	run, err := myb.RunWH(myb.WHRun{GM: gm, H: 2e5, Scheme: myb.WHSABA3, X: x, V: v, Times: times})
	must(t, err)
	status, xs, vs := cref.WH(gm, flatten(x), flatten(v), 2e5, times)
	okStatus(t, "WH", status)
	var gx, gv []float64
	for k := range times {
		gx = append(gx, flatten(run.X[k])...)
		gv = append(gv, flatten(run.V[k])...)
	}
	sameFloats(t, "positions", gx, xs)
	sameFloats(t, "velocities", gv, vs)
	if run.Stats.Steps == 0 {
		t.Errorf("no steps")
	}
	if _, err := myb.RunWH(myb.WHRun{GM: gm, H: 0, X: x, V: v, Times: times}); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a step of 0: %v", err)
	}
}

// ---- Estimation ------------------------------------------------------------------------------

func TestFilterAndMeasurements(t *testing.T) {
	model, err := myb.NewModel(myb.ModelSpec{Bodies: []myb.Body{{ID: 399, GM: gm}},
		Terms: []myb.Term{myb.PointMass{Name: "planet", GM: gm}}})
	must(t, err)
	defer model.Close()
	y0 := myb.State{7e6, 0, 0, 0, 6.5e3, 3.8e3}
	times := []float64{300, 600, 900, 1200, 1500, 1800}
	truth, err := model.Propagate(myb.Propagation{Y0: y0, TEnd: 2000, Times: times})
	must(t, err)
	var positions []float64
	var observations []myb.Observation
	for k, s := range truth.States {
		positions = append(positions, s[0], s[1], s[2])
		observations = append(observations, myb.Observation{T: times[k], Kind: myb.MeasurePosition,
			Value: []float64{s[0], s[1], s[2]}, Sigma: []float64{10, 10, 10}})
	}
	guess := y0
	guess[0] += 30
	covariance := make([]float64, 36)
	for i := range 6 {
		covariance[7*i] = 1e4
		if i >= 3 {
			covariance[7*i] = 1
		}
	}
	result, err := model.Filter(myb.FilterRequest{Y0: guess, Covariance: covariance, TEnd: 2000,
		Observations: observations, Smooth: true})
	must(t, err)
	status, innovations, stepT, stepU := cref.Filter(gm, guess[:], covariance, 2000, times, positions, 10)
	okStatus(t, "Filter", status)
	var gi, gt, gu []float64
	for _, in := range result.Innovations {
		gi = append(gi, in...)
	}
	for _, s := range result.Steps {
		gt = append(gt, s.T)
		gu = append(gu, s.Unknowns...)
	}
	sameFloats(t, "innovations", gi, innovations)
	sameFloats(t, "step times", gt, stepT)
	sameFloats(t, "smoothed unknowns", gu, stepU)
	if result.T != 1800 {
		t.Errorf("the filter stands at %v", result.T)
	}

	station := myb.Station{Body: 0, Ellipsoid: myb.Ellipsoid{A: radius}, Longitude: 0.3, Latitude: 0.5}
	values, partials, err := model.Measure(myb.MeasureRange, &station, nil, times, truth.States, true)
	must(t, err)
	status, ref, refPartials := cref.Range(gm, 0.3, 0.5, radius, times, flatten(truth.States))
	okStatus(t, "Range", status)
	var gv, gp []float64
	for k := range values {
		gv = append(gv, values[k]...)
		gp = append(gp, partials[k]...)
	}
	sameFloats(t, "ranges", gv, ref)
	sameFloats(t, "range partials", gp, refPartials)
	if _, _, err := model.Measure(myb.MeasureRange, nil, nil, times, truth.States, false); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a range without a station: %v", err)
	}
	if _, err := model.Filter(myb.FilterRequest{Y0: guess, Covariance: covariance[:6]}); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a short covariance: %v", err)
	}

	// the covariance carried along the orbit, in the radial axes
	var p0 [36]float64
	copy(p0[:], covariance)
	states, carried, err := model.PropagateCovariance(myb.OrbitStart{Y0: y0, TEnd: 2000, Covariance: p0}, times)
	must(t, err)
	if _, err := myb.CovariancesRTN(states, carried); err != nil {
		t.Errorf("CovariancesRTN: %v", err)
	}
	if carried[5][0] <= p0[0] {
		t.Errorf("the position covariance did not grow")
	}
}

// ---- The restricted problem ------------------------------------------------------------------

func TestCR3BPFamiliesAndTori(t *testing.T) {
	const mu = 1.215058560962404e-02
	lyapunov, period, err := myb.CR3BPCorrect(mu, myb.CR3BPCorrection{Symmetry: myb.SymmetryXZ},
		myb.State{8.2522726142834313e-01, 0, 0, 0, 1.0748883813645686e-01, 0})
	must(t, err)
	orbit, err := myb.CR3BPOrbitFromState(mu, lyapunov, period, 4)
	must(t, err)
	members, err := myb.CR3BPFamily(mu, orbit, 1e-3, 3, 0, 0)
	must(t, err)
	status, periods, states := cref.CR3BPFamily(mu, lyapunov[:], period, 4, 1e-3, 3)
	okStatus(t, "CR3BPFamily", status)
	var gp, gs []float64
	for _, m := range members {
		gp = append(gp, m.Period)
		gs = append(gs, flatten(m.States)...)
	}
	sameFloats(t, "periods", gp, periods)
	sameFloats(t, "family states", gs, states)
	phases := []float64{0, 0.25, 0.5}
	manifold, err := myb.CR3BPManifold(mu, lyapunov, period, phases, true, 1, 1e-6)
	must(t, err)
	for k, phase := range phases {
		status, y := cref.Manifold(mu, lyapunov[:], period, phase)
		okStatus(t, "Manifold", status)
		sameFloats(t, "manifold", manifold[k][:], y)
	}
	crossings, err := myb.CR3BPPropagate(mu, manifold, 5, &myb.CR3BPSection{Axis: 0, Value: 1 - mu})
	must(t, err)
	if len(crossings) != 3 {
		t.Errorf("crossings %d", len(crossings))
	}

	// tori about a distant retrograde orbit: the family, the spectrum and its normal rotation
	dro, droPeriod, err := myb.CR3BPCorrect(mu, myb.CR3BPCorrection{Symmetry: myb.SymmetryXZ},
		myb.State{8.2915209119468902e-01, 0, 0, 0, 4.9478973105640744e-01, 0})
	must(t, err)
	torus, err := myb.CR3BPTorusFromOrbit(mu, dro, droPeriod, 21, 0, 1e-3)
	must(t, err)
	tori, err := myb.CR3BPTorusFamily(mu, torus, myb.TorusJacobi, 1e-4, 2, 0, 0)
	must(t, err)
	spectrum, err := myb.CR3BPTorusStability(mu, tori[1], 0)
	must(t, err)
	status, rotations, normal := cref.Torus(mu, dro[:], droPeriod, 21, 1e-3, 1e-4)
	okStatus(t, "Torus", status)
	sameFloats(t, "rotations", []float64{tori[0].Rotation, tori[1].Rotation}, rotations)
	sameFloats(t, "normal rotation", []float64{spectrum.NormalRotation, spectrum.NormalEigenvalue[0],
		spectrum.NormalEigenvalue[1], spectrum.NormalRe[0][0], spectrum.NormalIm[0][5],
		float64(spectrum.Elliptic)}, normal)
	if spectrum.NormalRotation <= 0 {
		t.Errorf("an elliptic torus without its normal rotation")
	}
	if _, err := myb.CR3BPTorusManifold(tori[1], spectrum, []float64{0}, true, 1, 1e-6); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("the manifold of an elliptic torus: %v", err)
	}
	if _, err := myb.CR3BPOrbitFromState(mu, lyapunov, period, 11); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("11 patches: %v", err)
	}
}

// ---- Targeting and transfers ---------------------------------------------------------------

const sunGM = 1.32712440018e20

func circle(r, phase float64) myb.State {
	v := math.Sqrt(sunGM / r)
	return myb.State{r * math.Cos(phase), r * math.Sin(phase), 0, -v * math.Sin(phase), v * math.Cos(phase), 0}
}

func TestLambertGridAndTransfers(t *testing.T) {
	earth, mars := circle(1.496e11, 0), circle(2.279e11, 0.8)
	departures, err := myb.TimesFromSeconds([]float64{0, 2.592e6})
	must(t, err)
	arrivals, err := myb.TimesFromSeconds([]float64{1.728e7, 2.246e7})
	must(t, err)
	grid, err := myb.SolveLambertGrid(sunGM, myb.KeplerOrbit{MU: sunGM, State: earth},
		myb.KeplerOrbit{MU: sunGM, State: mars}, departures, arrivals, false, 0)
	must(t, err)
	for i, d := range []float64{0, 2.592e6} {
		for j, a := range []float64{1.728e7, 2.246e7} {
			status, v1, v2 := cref.LambertCell(sunGM, earth[:], mars[:], d, a)
			okStatus(t, "LambertCell", status)
			cell := grid.Cells[i][j]
			if len(cell) != 1 {
				t.Fatalf("cell %d %d: %d solutions", i, j, len(cell))
			}
			sameFloats(t, "v1", cell[0].V1[:], v1)
			sameFloats(t, "v2", cell[0].V2[:], v2)
		}
	}
	solutions, err := myb.LambertMany(sunGM, []myb.LambertProblem{{R1: [3]float64{1.496e11, 0, 0},
		R2: [3]float64{0, 2.279e11, 0}, TOF: 2e7}}, false, 1)
	must(t, err)
	if len(solutions[0]) == 0 {
		t.Errorf("no solution of Lambert")
	}

	// two orbits that cross: their conjunctions with the probability
	y1 := myb.State{7e6, 0, 0, 0, 7546, 0}
	y2 := myb.State{7e6, 0, 30, 0, 7540, 300}
	var covariance [36]float64
	for i := range 6 {
		covariance[7*i] = 100
		if i >= 3 {
			covariance[7*i] = 1e-2
		}
	}
	model, err := myb.NewModel(myb.ModelSpec{Terms: []myb.Term{myb.PointMass{Name: "planet", GM: gm}}})
	must(t, err)
	defer model.Close()
	start := func(y myb.State) myb.OrbitStart {
		return myb.OrbitStart{Options: myb.OrbitOptions{Variational: true}, Y0: y, TEnd: 12000, Covariance: covariance}
	}
	found, err := myb.FindConjunctions(model, start(y1), model, start(y2), 10, 6000, 10)
	must(t, err)
	status, count, firstT, firstP := cref.Conjunctions(gm, y1[:], y2[:], covariance[:], 10, 6000, 10)
	okStatus(t, "Conjunctions", status)
	if len(found) != count || count == 0 {
		t.Fatalf("conjunctions %d, the C calls %d", len(found), count)
	}
	sameFloats(t, "conjunction", []float64{found[0].Approach.T, found[0].Probability}, []float64{firstT, firstP})

	// the Q-law flown to a higher orbit
	target := []float64{7.1e6, 0.001, 0.9, 0, 0}
	weight := []float64{1, 0, 0, 0, 0}
	law := myb.QLaw{MU: gm, Target: [5]float64(target), Weight: [5]float64(weight)}
	y0 := [7]float64{7e6, 0, 0, 0, 4686, 5885, 1000}
	arrival, err := myb.QLawFly(law, myb.QLawEngine{Thrust: 1, Exhaust: 30000,
		Tolerance: [5]float64{1e3, 1e-3, 1e-3, 1e-3, 1e-3}}, nil, 0, 0, 0, y0, 2e5, 600)
	must(t, err)
	status, y, at := cref.QLawFly(gm, target, weight, 1, 30000, y0[:], 2e5, 600)
	okStatus(t, "QLawFly", status)
	sameFloats(t, "Q-law arrival", append(arrival.Y[:], arrival.T), append(y, at))
	if !arrival.Arrived {
		t.Errorf("the Q-law did not arrive")
	}
	elements, err := myb.ElementsFromStates(gm, []myb.State{y1, y2})
	must(t, err)
	if _, _, _, _, err := myb.QLawDirections(law, elements); err != nil {
		t.Errorf("QLawDirections: %v", err)
	}

	// Sims and Flanagan: the mismatch of the halves
	leg := myb.SimsFlanagan{MU: sunGM, Thrust: 0.5, Exhaust: 30000, Mass0: 1000, T1: 1.5e7,
		Start: earth, Target: circle(1.6e11, 0.5), Segments: 6}
	throttles := [][3]float64{{0.3, 0.1, 0}, {0.3, 0.1, 0}, {0, 0, 0}, {0, 0, 0.2}, {0.1, 0, 0}, {0.2, 0.2, 0}}
	mismatch, mass1, _, err := myb.SimsFlanaganMismatch(leg, throttles)
	must(t, err)
	status, refMismatch, refMass := cref.SimsFlanagan(sunGM, earth[:], leg.Target[:], 1.5e7, flatten(throttles))
	okStatus(t, "SimsFlanagan", status)
	sameFloats(t, "mismatch", append(mismatch[:], mass1), append(refMismatch, refMass))

	// relative motion and the B-plane over slices
	chief := []myb.State{y1, y1}
	deputy := []myb.State{y2, {7e6 + 100, 0, 0, 0, 7546, 1}}
	relative, err := myb.RelativeFromInertial(chief, deputy)
	must(t, err)
	inertial, err := myb.RelativeToInertial(chief, relative)
	must(t, err)
	for k := range inertial {
		for i := range 6 {
			if math.Abs(inertial[k][i]-deputy[k][i]) > 1e-6*math.Max(1, math.Abs(deputy[k][i])) {
				t.Errorf("the relative state does not come back")
			}
		}
	}
	if _, err := myb.BPlanes(gm, []myb.State{{7e6, 0, 0, 0, 12e3, 0}}, [3]float64{0, 0, 1}); err != nil {
		t.Errorf("BPlanes: %v", err)
	}
	if _, err := myb.BPlanes(gm, []myb.State{y1}, [3]float64{0, 0, 1}); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("the B-plane of an ellipse: %v", err)
	}
}

func TestChainWithPoweredFlyby(t *testing.T) {
	states := []myb.State{circle(1.496e11, 0), circle(1.082e11, 0.9), circle(2.279e11, 2.1)}
	gms := []float64{3.986004418e14, 3.24858592e14, 4.282837e13}
	periapsis := []float64{6678e3, 6351.8e3, 3696e3}
	impulse := []float64{0, 3000, 0}
	bounds := []float64{0, 1.728e6, 1.2e7, 1.4e7, 1.6e7, 2.0e7}
	chain := myb.Chain{MU: sunGM, Thrust: 1, Exhaust: 30000, Mass0: 1000, Launch: 2500,
		Departure: [2]float64{bounds[0], bounds[1]},
		Legs:      []myb.ChainLeg{{Segments: 6, Flight: [2]float64{bounds[2], bounds[3]}}, {Segments: 6, Flight: [2]float64{bounds[4], bounds[5]}}}}
	for k := range states {
		chain.Bodies = append(chain.Bodies, myb.ChainBody{State: myb.KeplerOrbit{MU: sunGM, State: states[k]},
			GM: gms[k], Periapsis: periapsis[k], ImpulseExhaust: impulse[k]})
	}
	unknowns, equalities, inequalities, err := chain.Size()
	must(t, err)
	if unknowns != 7*2+1+3*12+2 || equalities != 13 || inequalities != 4 {
		t.Errorf("the size of the chain: %d %d %d", unknowns, equalities, inequalities)
	}
	x, err := chain.Guess()
	must(t, err)
	evaluation, err := chain.Evaluate(x, true)
	must(t, err)
	status, refX, refP, refC := cref.Chain(sunGM, flatten(states), gms, periapsis, impulse, bounds,
		unknowns, equalities+inequalities)
	okStatus(t, "Chain", status)
	sameFloats(t, "guess", x, refX)
	sameFloats(t, "propellant and constraints", append([]float64{evaluation.Propellant}, evaluation.Constraints...),
		append([]float64{refP}, refC...))
	if len(evaluation.Jacobian) != unknowns*(equalities+inequalities) {
		t.Errorf("a Jacobian of %d values", len(evaluation.Jacobian))
	}
	dv, turn, _, err := myb.ChainFlyby(gms[1], 3000, 3500, 7e6)
	must(t, err)
	if dv <= 0 || turn <= 0 {
		t.Errorf("a powered flyby: dv %v turn %v", dv, turn)
	}
	chain.Bodies[0].ImpulseExhaust = 3000 // an impulse at the departure
	if _, _, _, err := chain.Size(); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("an impulse at the departure: %v", err)
	}
	if _, _, _, err := myb.ChainFlyby(gms[1], -1, 3500, 7e6); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a negative excess speed: %v", err)
	}
}

// ---- Attitude, observables ---------------------------------------------------------------

func TestAttitudeAndLightTime(t *testing.T) {
	sun := myb.State{1.496e11, 0, 0, 0, 29780, 0}
	model, err := myb.NewModel(myb.ModelSpec{Bodies: []myb.Body{{ID: 399, GM: gm},
		{ID: 10, GM: sunGM, State: myb.KeplerOrbit{MU: sunGM, State: sun}}}})
	must(t, err)
	defer model.Close()
	times := []float64{0, 1000, 5000}
	states := []myb.State{{7e6, 0, 1e5, 0, 7.5e3, 0}, {6e6, 3e6, 1e5, -3e3, 6.5e3, 100}, {-7e6, 1e6, 0, 0, -7.5e3, 0}}
	m, err := model.Attitudes(myb.AttitudeSunNadir{Body: 0, Source: 1}, times, states)
	must(t, err)
	status, ref := cref.SunNadir(gm, sun[:], sunGM, times, flatten(states))
	okStatus(t, "SunNadir", status)
	sameFloats(t, "attitudes", flatten(m), ref)
	if _, err := model.Attitudes(myb.AttitudeOrbital{Body: 5}, times, states); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a law of a body the model does not have: %v", err)
	}
	spin, err := myb.Attitudes(myb.AttitudeSpin{M0: [9]float64{1, 0, 0, 0, 1, 0, 0, 0, 1}, Rate: 1e-3}, times, states)
	must(t, err)
	if spin[0][0] != 1 {
		t.Errorf("a spin at its phase 0")
	}
	q, l, err := myb.FreeRotation([3]float64{1, 2, 3}, [4]float64{1, 0, 0, 0}, [3]float64{0.1, 0.2, 0.3}, times)
	must(t, err)
	if len(q) != 3 || l[2] == [3]float64{} {
		t.Errorf("a free rotation")
	}

	const c = 299792458.0
	target := myb.State{3.84e8, 0, 0, 0, 1022, 0}
	instants, err := myb.TimesFromSeconds(times)
	must(t, err)
	paths, err := myb.LightTimes(myb.LightModel{C: c, Gamma: 1}, myb.KeplerOrbit{MU: gm, State: target},
		instants, states, myb.Reception)
	must(t, err)
	for k := range times {
		status, lt, r := cref.LightTime(c, gm, target[:], states[k][:], times[k])
		okStatus(t, "LightTime", status)
		sameFloats(t, "light time", []float64{paths[k].LightTime, paths[k].Range}, []float64{lt, r})
	}
	ra, dec := myb.RaDec([][3]float64{{1, 0, 0}, {0, 0, 1}})
	if ra[0] != 0 || dec[1] != math.Pi/2 {
		t.Errorf("directions %v %v", ra, dec)
	}
	if _, err := myb.Aberrations([][3]float64{{1, 0, 0}}, nil, c); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("slices of different lengths: %v", err)
	}
}

// ---- Files ---------------------------------------------------------------------------------

func TestMessages(t *testing.T) {
	cal := myb.Calendar{Year: 2026, Month: 11, Day: 14, Hour: 12, Minute: 30, Second: 15.5}
	cov := [36]float64{}
	for i := range 6 {
		for j := range 6 {
			cov[6*i+j] = float64(1 + i + j)
		}
	}
	opm := &myb.OPM{Originator: "myb", ObjectName: "TEST", ObjectID: "2026-001A", CenterName: "EARTH",
		RefFrame: "EME2000", TimeSystem: "UTC", CreationDate: cal, Epoch: cal,
		State:      myb.State{7e6, 1e3, -2.5e3, 0.25, 7.5e3, -1e3},
		Elements:   &myb.OPMElements{SemiMajorAxis: 7e6, Eccentricity: 0.001, Inclination: 0.5, GM: gm},
		Spacecraft: &myb.OPMSpacecraft{Mass: 500, SolarRadArea: 4, SolarRadCoeff: 1.3, DragArea: 4, DragCoeff: 2.2},
		Maneuvers:  []myb.OPMManeuver{{EpochIgnition: cal, Duration: 60, DeltaMass: -1, RefFrame: "RTN", DV: [3]float64{1, 0, 0}}},
		Covariance: &cov}
	for _, notation := range []myb.Notation{myb.KVN, myb.XML} {
		text, err := myb.WriteOPM(opm, notation)
		must(t, err)
		back, err := myb.ParseOPM([]byte(text))
		must(t, err)
		if back.ObjectID != opm.ObjectID || back.State != opm.State || back.Epoch != opm.Epoch ||
			len(back.Maneuvers) != 1 || back.Maneuvers[0].DV != opm.Maneuvers[0].DV || back.Spacecraft == nil {
			t.Errorf("OPM %d read back: %+v", notation, back)
		}
	}
	omm := &myb.OMM{Originator: "myb", ObjectName: "ISS", ObjectID: "1998-067A", CenterName: "EARTH",
		RefFrame: "TEME", TimeSystem: "UTC", MeanElementTheory: "SGP4", CreationDate: cal, Epoch: cal,
		HasMeanMotion: true, MeanMotion: 15.5, Eccentricity: 0.0005, Inclination: 51.6, RAAN: 10,
		ArgPericenter: 20, MeanAnomaly: 30,
		TLE: &myb.OMMTLE{ClassificationType: 'U', NoradCatID: 25544, ElementSetNo: 999, RevAtEpoch: 1, BStar: 1e-4}}
	text, err := myb.WriteOMM(omm, myb.KVN)
	must(t, err)
	backOMM, err := myb.ParseOMM([]byte(text))
	must(t, err)
	if backOMM.MeanMotion != omm.MeanMotion || backOMM.TLE == nil || backOMM.TLE.NoradCatID != 25544 {
		t.Errorf("OMM read back: %+v", backOMM)
	}
	if _, err := myb.WriteOPM(&myb.OPM{ObjectName: string(make([]byte, 80))}, myb.KVN); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a name of 80 characters: %v", err)
	}
	if _, err := myb.ParseOMM([]byte("CCSDS_OMM_VERS = 2.0\n")); !errors.Is(err, myb.ErrFormat) {
		t.Errorf("a broken OMM: %v", err)
	}

	sp3 := &myb.SP3{Version: 'd', FileType: 'L', TimeSystem: "GPS", Frame: "SLR20", OrbitType: "FIT",
		Agency: "MYB", DataUsed: "SLR", Satellites: []string{"L51", "L52"},
		Epochs:   []myb.Calendar{{Year: 2026, Month: 1, Day: 1}, {Year: 2026, Month: 1, Day: 1, Minute: 2}},
		Interval: 120, Positions: [][3]float64{{7e6, 1e3, 2e3}, {-7e6, 5e3, 1e6}, {7.001e6, 1.1e3, 2e3}, {math.NaN(), 0, 0}}}
	text, err = myb.WriteSP3(sp3)
	must(t, err)
	backSP3, err := myb.ParseSP3([]byte(text))
	must(t, err)
	if len(backSP3.Satellites) != 2 || backSP3.Positions[0] != sp3.Positions[0] || !math.IsNaN(backSP3.Positions[3][0]) {
		t.Errorf("SP3 read back: %v", backSP3.Positions)
	}
}

func TestCKAndSCLK(t *testing.T) {
	pool, err := myb.NewKernelPool()
	must(t, err)
	defer pool.Close()
	kernel := testDataPath("ck/test.tsc")
	must(t, pool.LoadFile(kernel))
	clock, err := myb.NewSCLK(pool, 9099)
	must(t, err)
	defer clock.Close()
	ck, err := myb.NewCK()
	must(t, err)
	defer ck.Close()
	path := testDataPath("ck/test.bc")
	must(t, ck.LoadFile(path))
	segments, err := ck.Segments()
	must(t, err)
	var seg myb.CKSegmentInfo
	for _, s := range segments {
		if s.Instrument == -9099001 {
			seg = s
			break
		}
	}
	if seg.Instrument == 0 {
		t.Fatalf("no segment of the instrument -9099001")
	}
	ticks := []float64{seg.Start, seg.Start/2 + seg.Stop/2, seg.Stop}
	pointings, err := ck.Pointings(seg.Instrument, ticks, 0, false)
	must(t, err)
	times, err := clock.Times(ticks)
	must(t, err)
	for k, tick := range ticks {
		status, q := cref.CKPointing(path, int(seg.Instrument), tick)
		okStatus(t, "CKPointing", status)
		sameFloats(t, "quaternion", pointings[k].Q[:], q)
		status, s, f := cref.SCLKTime(kernel, 9099, tick)
		okStatus(t, "SCLKTime", status)
		if times[k] != (myb.Time{S: s, F: f}) {
			t.Errorf("the instant of tick %v: %+v, the C calls %d %v", tick, times[k], s, f)
		}
	}
	back, err := clock.Ticks(times)
	must(t, err)
	if len(back) != 3 {
		t.Errorf("ticks back")
	}
	if _, err := ck.Pointings(seg.Instrument, []float64{seg.Stop + 1e12}, 0, false); err == nil {
		t.Errorf("a pointing beyond the coverage: %v", err)
	}

	// a CK written and read back, its instrument through the clock as an orientation
	quaternions := [][4]float64{{1, 0, 0, 0}, {0.9, 0.1, 0, 0}, {0.8, 0.2, 0.1, 0}}
	norm := func(q [4]float64) [4]float64 {
		s := math.Sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3])
		return [4]float64{q[0] / s, q[1] / s, q[2] / s, q[3] / s}
	}
	image, err := myb.EncodeCK("myb test", []myb.CKSegment{{Instrument: -77001, Frame: 1, Name: "TEST",
		Type: 3, Ticks: []float64{1e6, 2e6, 3e6}, Quaternions: quaternions}})
	must(t, err)
	written, err := myb.NewCK()
	must(t, err)
	defer written.Close()
	must(t, written.LoadBytes(image))
	got, err := written.Pointings(-77001, []float64{1e6, 3e6}, 0, false)
	must(t, err)
	for k, i := range []int{0, 2} {
		want := norm(quaternions[i])
		for j := range 4 {
			if math.Abs(got[k].Q[j]-want[j]) > 1e-15 && math.Abs(got[k].Q[j]+want[j]) > 1e-15 {
				t.Errorf("the pointing %d read back: %v against %v", i, got[k].Q, want)
			}
		}
	}
	if _, err := myb.EncodeCK("myb", []myb.CKSegment{{Instrument: 1, Frame: 1, Name: "X", Type: 5,
		Ticks: []float64{1}, Quaternions: quaternions[:1]}}); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a segment of type 5: %v", err)
	}
	orientation := myb.CKOrientation{CK: ck, SCLK: clock, Instrument: -9099001, Base: 1}
	if _, err := orientation.Sample(times[1:2]); err != nil {
		t.Errorf("the CK as an orientation: %v", err)
	}
}
