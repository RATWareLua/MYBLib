// SPDX-License-Identifier: BSD-3-Clause

// The remaining paths: a provider in Go against the one in C (bit for bit), a system with a
// field, an orientation and an action in Go, shooting with a residual in Go, orbits from
// directions, the probabilities of collision, the long encounter, a leg of Sims and Flanagan
// solved, and the errors of the optimal transfers.
package myb_test

import (
	"errors"
	"math"
	"testing"

	"github.com/RATWareLua/MYBLib/bindings/go/myb"
)

func TestGoProviderMatchesC(t *testing.T) {
	sun := myb.State{1.496e11, 0, 0, 0, 29780, 0}
	calls := 0
	goSun := myb.GoState{Func: func(at myb.Time) (myb.State, error) {
		calls++
		seconds, err := myb.Seconds(at) // the epoch is 0: the seconds are the offset itself
		if err != nil {
			return myb.State{}, err
		}
		return myb.KeplerPropagate(sunGM, sun, seconds)
	}}
	instants, err := myb.TimesFromSeconds([]float64{0, 1e5, 3.3e6 + 0.25})
	must(t, err)
	byGo, err := myb.States(goSun, instants)
	must(t, err)
	byC, err := myb.States(myb.KeplerOrbit{MU: sunGM, State: sun}, instants)
	must(t, err)
	sameBytes(t, "a provider in Go and in C", byGo, byC)
	if calls != 3 {
		t.Errorf("the Go provider was called %d times", calls)
	}
	// in a model: the attitude that follows the Sun, both ways
	attitude := func(src myb.StateSource) [][9]float64 {
		model, err := myb.NewModel(myb.ModelSpec{Bodies: []myb.Body{{ID: 399, GM: gm}, {ID: 10, GM: sunGM, State: src}}})
		must(t, err)
		defer model.Close()
		m, err := model.Attitudes(myb.AttitudeSunNadir{Body: 0, Source: 1}, []float64{0, 5000},
			[]myb.State{{7e6, 0, 1e5, 0, 7.5e3, 0}, {-7e6, 1e6, 0, 0, -7.5e3, 0}})
		must(t, err)
		return m
	}
	sameBytes(t, "attitudes with the Sun in Go and in C", attitude(goSun), attitude(myb.KeplerOrbit{MU: sunGM, State: sun}))
	// a Go provider that fails stops the call with its error
	refused := errors.New("no ephemeris here")
	_, err = myb.States(myb.GoState{Func: func(myb.Time) (myb.State, error) { return myb.State{}, refused }}, instants)
	if !errors.Is(err, refused) {
		t.Errorf("a failing provider: %v", err)
	}
}

func TestSystemLevels(t *testing.T) {
	field, err := myb.NewGravityField(gm, radius, 2, []float64{1, 0, 0, -4.84e-4, 0, 0}, make([]float64, 6))
	must(t, err)
	defer field.Close()
	pushes := 0
	system, err := myb.NewSystem(myb.SystemSpec{
		Bodies: []myb.SystemBody{
			{GM: gm, Field: field, Orientation: myb.FixedOrientation{Rotation: myb.Rotation{M: [9]float64{1, 0, 0, 0, 1, 0, 0, 0, 1}}}},
			{GM: 1},
		},
		Links: []myb.Link{
			myb.FieldLink{Body: 0, Point: 1},
			myb.ActionLink{Body: 1, Func: func(_ float64, view myb.SystemView) ([3]float64, [3]float64, error) {
				pushes++
				if view.Body != 1 || len(view.X) != 2 {
					return [3]float64{}, [3]float64{}, errors.New("a wrong view")
				}
				return [3]float64{0, 1e-6, 0}, [3]float64{}, nil
			}},
		},
	})
	must(t, err)
	defer system.Close()
	n, err := system.Coordinates()
	must(t, err)
	y0, err := system.PackState([]myb.State{{}, {7e6, 0, 0, 0, 7546, 0}}, nil, nil)
	must(t, err)
	run, err := system.Run(myb.SystemRun{Y0: y0, Times: []float64{600, 1200}})
	must(t, err)
	if n != 6 || len(run.States) != 2 || pushes == 0 || len(run.Events) != 0 {
		t.Errorf("a system with levels: %d coordinates, %d states, %d pushes", n, len(run.States), pushes)
	}
}

func TestShootingIODAndProbabilities(t *testing.T) {
	// Newton's method with a residual in Go: u^2 = 2
	u := []float64{1}
	iterations, residual, err := myb.Shoot(func(u, f, jacobian []float64) error {
		f[0] = u[0]*u[0] - 2
		if jacobian != nil {
			jacobian[0] = 2 * u[0]
		}
		return nil
	}, u, 1, 1e-14, 0)
	must(t, err)
	if math.Abs(u[0]-math.Sqrt2) > 1e-14 || iterations == 0 || residual > 1e-14 {
		t.Errorf("shooting: %v after %d iterations", u[0], iterations)
	}
	if _, _, err := myb.Shoot(nil, u, 1, 1e-14, 0); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("shooting without a function: %v", err)
	}

	// an orbit from three directions seen from the centre
	y0 := myb.State{7e6, 1e5, -2e5, 0, 7.5e3, 1e3}
	states, err := myb.KeplerPropagateMany(gm, y0, []float64{0, 300, 600})
	must(t, err)
	var sights myb.IODSights
	for k, s := range states {
		r := math.Sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2])
		sights.T[k] = 300 * float64(k)
		sights.Sight[k] = [3]float64{s[0] / r, s[1] / r, s[2] / r}
		sights.Observer[k] = [3]float64{1e3 * float64(k+1), -2e3, 5e2} // observers off the centre
		sights.Sight[k] = direction(s, sights.Observer[k])
	}
	solutions, err := myb.IODGauss(gm, sights)
	must(t, err)
	best := math.Inf(1)
	for _, s := range solutions {
		best = math.Min(best, math.Abs(s.R[0]-states[1][0]))
	}
	if best > 1 {
		t.Errorf("Gauss's method: %d solutions, the nearest %v m off", len(solutions), best)
	}

	// the probabilities of an encounter
	encounter, err := myb.EncounterPlane([3]float64{20, 5, 0}, [3]float64{0, 0, 1e4},
		[9]float64{400, 0, 0, 0, 400, 0, 0, 0, 400}, 10)
	must(t, err)
	p, err := myb.CollisionProbability(encounter)
	must(t, err)
	chan_ := myb.CollisionProbabilityChan(encounter)
	pmax, factor, err := myb.CollisionProbabilityMax(encounter)
	must(t, err)
	if math.Abs(chan_/p-1) > 1e-10 || pmax < p || factor <= 0 {
		t.Errorf("probabilities %v %v %v (factor %v)", p, chan_, pmax, factor)
	}
	if _, err := myb.EncounterPlane([3]float64{20, 5, 0}, [3]float64{}, [9]float64{}, 10); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("an encounter without relative velocity: %v", err)
	}

	// the long encounter of two orbits
	var covariance [36]float64
	for i := range 6 {
		covariance[7*i] = 100
		if i >= 3 {
			covariance[7*i] = 1e-4
		}
	}
	model, err := myb.NewModel(myb.ModelSpec{Terms: []myb.Term{myb.PointMass{Name: "planet", GM: gm}}})
	must(t, err)
	defer model.Close()
	a := myb.OrbitStart{Y0: myb.State{7e6, 0, 0, 0, 7546, 0}, TEnd: 2000, Covariance: covariance}
	b := myb.OrbitStart{Y0: myb.State{7e6, 30, 0, 0, 7546, 0.05}, TEnd: 2000, Covariance: covariance}
	probability, _, err := myb.LongEncounter(model, a, model, b, 10, 600, 30)
	must(t, err)
	if !(probability >= 0) {
		t.Errorf("the long encounter: %v", probability)
	}
}

func direction(s myb.State, from [3]float64) [3]float64 {
	d := [3]float64{s[0] - from[0], s[1] - from[1], s[2] - from[2]}
	n := math.Sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
	return [3]float64{d[0] / n, d[1] / n, d[2] / n}
}

func TestOptimalTransfers(t *testing.T) {
	// the transfer of the tests of the library (7000 km to the circle of 8000 km), 10 segments
	low := func(r, angle float64) myb.State {
		v := math.Sqrt(gm / r)
		return myb.State{r * math.Cos(angle), r * math.Sin(angle), 0, -v * math.Sin(angle), v * math.Cos(angle), 0}
	}
	leg := myb.SimsFlanagan{MU: gm, Thrust: 300, Exhaust: 3000, T0: 100, Mass0: 1000, T1: 4600,
		Start: low(7e6, 0), Target: low(8e6, 4.4), Segments: 10}
	result, err := myb.SimsFlanaganSolve(leg, nil, nil)
	must(t, err)
	if !(result.Mass1 < 1000 && result.Mass1 > 500) || math.Abs(result.Mismatch[0]) > 1e3 {
		t.Errorf("Sims and Flanagan: %v kg left, mismatch %v", result.Mass1, result.Mismatch)
	}
	mismatch, mass1, _, err := myb.SimsFlanaganMismatch(leg, result.Throttles)
	must(t, err)
	if mass1 != result.Mass1 || mismatch != result.Mismatch {
		t.Errorf("the solution's mismatch differs from its evaluation")
	}
	if _, err := myb.IndirectSolve(myb.IndirectProblem{MU: gm, Thrust: -1, Exhaust: 3000, Mass0: 100, T1: 1000,
		Start: myb.State{7e6, 0, 0, 0, 7546, 0}, Target: myb.State{8e6, 0, 0, 0, 7059, 0}}, nil); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a negative thrust: %v", err)
	}
	if _, err := myb.SimsFlanaganSolve(leg, nil, make([][3]float64, 2)); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("throttles of another count: %v", err)
	}
}
