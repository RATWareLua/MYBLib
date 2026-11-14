// SPDX-License-Identifier: BSD-3-Clause

// The physics of the last parts of phase 21 through the binding against the same calls in C
// (package internal/cref), bit for bit: the break schemes of IAS15 and its timescales, the
// friction that sticks (the deflections of the springs in the state of a system, the force of a
// contact over slices of states), the PN pairs with the new spin and tail terms; and the errors
// of bad input, with the statuses and the indices that C gives.
package myb_test

import (
	"errors"
	"math"
	"testing"

	"github.com/RATWareLua/MYBLib/bindings/go/myb"
	"github.com/RATWareLua/MYBLib/bindings/go/myb/internal/cref"
)

// ---- Break schemes and timescales -----------------------------------------------------------

func TestBreakSchemesAndTimescales(t *testing.T) {
	const (
		sunGM      = 1.32712440018e20
		sunMU      = 1.32712838618e20
		sunRadius  = 6.957e8
		luminosity = 3.828e26
		area       = 0.2
		cr         = 1.5
	)
	sun := myb.State{1.496e11, 0, 0, 0, -29780, 0}
	y0 := myb.State{1.2e7, 0, 0, 0, 5763, 10}
	times := []float64{21600, 43200, day}
	model := func(scheme myb.BreakScheme) (*myb.Model, error) {
		return myb.NewModel(myb.ModelSpec{
			Bodies: []myb.Body{
				{ID: 399, GM: gm, Radius: radius},
				{ID: 10, GM: sunGM, Radius: sunRadius, Luminosity: luminosity,
					State: myb.KeplerOrbit{MU: sunMU, State: sun}},
			},
			Terms: []myb.Term{
				myb.PointMass{Name: "earth", GM: gm},
				myb.Radiation{Name: "srp", Source: 1, Object: myb.Object{Mass: 1, Area: area, Cr: cr},
					C: 299792458, Occulters: []int{0}, Breaks: true},
			},
			BreakScheme: scheme,
		})
	}
	cases := []struct {
		scheme    myb.BreakScheme
		timescale myb.IAS15Timescale
	}{
		{myb.BreaksRetake, myb.TimescaleAcceleration},
		{myb.BreaksPredict, myb.TimescaleAcceleration},
		{myb.BreaksPredict, myb.TimescaleJerk},
	}
	var steps []int64
	for _, c := range cases {
		m, err := model(c.scheme)
		must(t, err)
		defer m.Close()
		traj, err := m.Propagate(myb.Propagation{Options: myb.OrbitOptions{Timescale: c.timescale},
			Y0: y0, TEnd: day, Times: times})
		must(t, err)
		status, states, cSteps := cref.ShadowOrbit(gm, radius, sunGM, sunRadius, sun[:], sunMU, luminosity, area,
			cr, int(c.scheme), int(c.timescale), y0[:], day, times)
		okStatus(t, "ShadowOrbit", status)
		sameFloats(t, "states", flatten(traj.States), states)
		if traj.Steps != cSteps {
			t.Errorf("scheme %d, timescale %d: %d steps, the C calls %d", c.scheme, c.timescale,
				traj.Steps, cSteps)
		}
		steps = append(steps, traj.Steps)
	}
	if steps[0] == steps[1] && steps[1] == steps[2] {
		t.Errorf("the schemes and the timescales took the same steps: %v", steps)
	}

	// errors: an unknown scheme (a contract of the library, refused here), an unknown timescale
	var me *myb.Error
	if _, err := model(myb.BreakScheme(7)); !errors.As(err, &me) || me.Status != myb.Domain {
		t.Errorf("an unknown break scheme: %v", err)
	}
	m, err := model(myb.BreaksRetake)
	must(t, err)
	defer m.Close()
	_, err = m.Propagate(myb.Propagation{Options: myb.OrbitOptions{Timescale: 9}, Y0: y0,
		TEnd: day, Times: times})
	status, _, _ := cref.ShadowOrbit(gm, radius, sunGM, sunRadius, sun[:], sunMU, luminosity, area, cr, 0, 9,
		y0[:], day, times)
	if !errors.As(err, &me) || int(me.Status) != status || status == 0 || me.Index != 0 {
		t.Errorf("an unknown timescale: %v, the C status %d", err, status)
	}
}

// ---- The friction that sticks ---------------------------------------------------------------

const moonGM, moonRadius = 4.9028e12, 1737400.0

var (
	plane     = myb.ContactSurface{Shape: myb.ContactPlane, Point: [3]float64{0, 0, moonRadius}, Normal: [3]float64{0, 0, 1}}
	stickLaw  = myb.ContactLaw{Stiffness: 1e4, Exponent: 1, Damping: 20, Friction: 0.6, Kind: myb.FrictionStick, TangentialStiffness: 1e4, TangentialDamping: 20}
	smoothLaw = myb.ContactLaw{Stiffness: 1e4, Exponent: 1, Damping: 20, Friction: 0.6, Slip: 0.01}
)

func lawValues(l myb.ContactLaw) []float64 {
	return []float64{l.Stiffness, l.Exponent, l.Damping, l.Friction, l.Slip, l.TangentialStiffness,
		l.TangentialDamping}
}

func pointValues(points []myb.ContactPoint) []float64 {
	var out []float64
	for _, p := range points {
		out = append(out, p.Offset[0], p.Offset[1], p.Offset[2], p.Radius)
	}
	return out
}

func TestSystemWithStickingFriction(t *testing.T) {
	gms := []float64{moonGM, 1}
	points := []myb.ContactPoint{{Offset: [3]float64{0, 0, 0}, Radius: 1}}
	contact := myb.ContactLink{Body: 1, Ground: 0, Surface: plane, Law: stickLaw, Points: points}
	system, err := myb.NewSystem(myb.SystemSpec{
		Bodies: []myb.SystemBody{{GM: gms[0]}, {GM: gms[1]}},
		Links:  []myb.Link{contact},
	})
	must(t, err)
	defer system.Close()
	if n, err := system.Coordinates(); err != nil || n != 9 {
		t.Fatalf("coordinates %d (%v), not 6 and 3 for the point", n, err)
	}
	bodies := []myb.State{{0, 0, 0, 0, 0, 0}, {0, 0, moonRadius + 1 - 1.6e-4, 0.01, 0, 0}}
	deflections := [][3]float64{{1e-5, 0, 0}}
	y0, err := system.PackState(bodies, nil, nil)
	must(t, err)
	must(t, system.PackContacts(y0, []myb.ContactDeflections{{Contact: 0, Deflections: deflections}}))
	a, err := system.Accelerations([]float64{0}, [][]float64{y0})
	must(t, err)
	times := []float64{0.05, 0.1, 0.2}
	run, err := system.Run(myb.SystemRun{Options: myb.IAS15Options{Timescale: myb.TimescaleJerk},
		Y0: y0, Times: times})
	must(t, err)
	back, err := system.UnpackContacts(run.States[len(times)-1], []int{0})
	must(t, err)

	status, ca, states, final := cref.ContactSystem(gms, plane.Point[:], plane.Normal[:],
		lawValues(stickLaw), pointValues(points), flatten(bodies), flatten(deflections), times)
	okStatus(t, "ContactSystem", status)
	sameFloats(t, "accelerations", a[0], ca)
	var flat []float64
	for _, s := range run.States {
		flat = append(flat, s...)
	}
	sameFloats(t, "states", flat, states)
	if len(back) != 1 || back[0].Contact != 0 {
		t.Fatalf("unpacked %+v", back)
	}
	sameFloats(t, "deflections", flatten(back[0].Deflections), final)
	if back[0].Deflections[0] == deflections[0] {
		t.Errorf("the springs of the friction did not move")
	}
	again, err := system.UnpackContacts(y0, []int{0})
	must(t, err)
	if again[0].Deflections[0] != deflections[0] {
		t.Errorf("the deflections unpacked are not those packed: %v", again[0].Deflections)
	}

	// errors: no such contact, a count of points that is not the contact's, a contact whose
	// friction does not stick, an unknown friction, IAS15 with an unknown timescale
	var me *myb.Error
	if _, err := system.UnpackContacts(y0, []int{0, 1}); !errors.As(err, &me) || me.Status != myb.Domain || me.Index != 1 {
		t.Errorf("no such contact: %v", err)
	}
	two := [][3]float64{{0, 0, 0}, {0, 0, 0}}
	if err := system.PackContacts(y0, []myb.ContactDeflections{{Contact: 0, Deflections: two}}); !errors.As(err, &me) || me.Status != myb.Domain || me.Index != 0 {
		t.Errorf("two deflections for one point: %v", err)
	}
	smooth := contact
	smooth.Law = smoothLaw
	other, err := myb.NewSystem(myb.SystemSpec{Bodies: []myb.SystemBody{{GM: gms[0]}, {GM: gms[1]}},
		Links: []myb.Link{smooth}})
	must(t, err)
	defer other.Close()
	y1, err := other.PackState(bodies, nil, nil)
	must(t, err)
	if err := other.PackContacts(y1, []myb.ContactDeflections{{Contact: 0, Deflections: deflections}}); !errors.As(err, &me) || me.Status != myb.Domain {
		t.Errorf("deflections of a smooth friction: %v", err)
	}
	unknown := contact
	unknown.Law.Kind = 5
	_, err = myb.NewSystem(myb.SystemSpec{Bodies: []myb.SystemBody{{GM: gms[0]}, {GM: gms[1]}},
		Links: []myb.Link{unknown}})
	if !errors.As(err, &me) || me.Status != myb.Domain || me.Message != "link" || me.Index != 0 {
		t.Errorf("an unknown friction: %v", err)
	}
	if _, err := system.Run(myb.SystemRun{Options: myb.IAS15Options{Timescale: 3}, Y0: y0, Times: times}); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("an unknown timescale: %v", err)
	}
}

// The contacts are numbered among those in place: a merger that takes a body away drops its
// contacts, and those after them move up.
func TestStickingContactsThroughAMerger(t *testing.T) {
	lander := []myb.ContactPoint{{Radius: 1}}
	feet := []myb.ContactPoint{{Offset: [3]float64{0.5, 0, 0}, Radius: 0.5}, {Offset: [3]float64{-0.5, 0, 0}, Radius: 0.5}}
	side := myb.ContactSurface{Shape: myb.ContactPlane, Point: [3]float64{0, -moonRadius, 0}, Normal: [3]float64{0, -1, 0}}
	system, err := myb.NewSystem(myb.SystemSpec{
		Bodies: []myb.SystemBody{{GM: moonGM, Radius: moonRadius}, {GM: 1}, {GM: 2, Radius: 1}},
		Links: []myb.Link{
			myb.ContactLink{Body: 2, Ground: 0, Surface: side, Law: stickLaw, Points: feet},
			myb.ContactLink{Body: 1, Ground: 0, Surface: plane, Law: stickLaw, Points: lander},
		},
	})
	must(t, err)
	defer system.Close()
	if n, err := system.Coordinates(); err != nil || n != 18 {
		t.Fatalf("coordinates %d (%v), not 9, 6 and 3", n, err)
	}
	bodies := []myb.State{
		{0, 0, 0, 0, 0, 0},
		{0, 0, moonRadius + 1 - 1.6e-4, 0.01, 0, 0},
		{0, -moonRadius - 10, 0, 0, 100, 0},
	}
	y0, err := system.PackState(bodies, nil, nil)
	must(t, err)
	must(t, system.PackContacts(y0, []myb.ContactDeflections{
		{Contact: 1, Deflections: [][3]float64{{1e-5, 0, 0}}},
		{Contact: 0, Deflections: [][3]float64{{2e-6, 0, 0}, {0, 3e-6, 0}}},
	}))
	before, err := system.UnpackContacts(y0, []int{0, 1})
	must(t, err)
	if len(before[0].Deflections) != 2 || before[0].Deflections[1] != [3]float64{0, 3e-6, 0} ||
		len(before[1].Deflections) != 1 || before[1].Deflections[0] != [3]float64{1e-5, 0, 0} {
		t.Fatalf("unpacked %+v", before)
	}
	run, err := system.Run(myb.SystemRun{Options: myb.IAS15Options{Timescale: myb.TimescaleJerk},
		Y0: y0, Times: []float64{0.2}, OnContact: myb.MergeAtContact})
	must(t, err)
	if len(run.Events) != 1 || run.Events[0].Survivor != 0 {
		t.Fatalf("events %+v", run.Events)
	}
	if present, _, err := system.Present(2); err != nil || present {
		t.Fatalf("the body merged is still there (%v)", err)
	}
	after, err := system.UnpackContacts(run.Final, []int{0})
	must(t, err)
	if len(after[0].Deflections) != 1 {
		t.Errorf("contact 0 after the merger is not the lander's: %+v", after)
	}
	var me *myb.Error
	if _, err := system.UnpackContacts(run.Final, []int{1}); !errors.As(err, &me) || me.Status != myb.Domain || me.Index != 0 {
		t.Errorf("a contact dropped: %v", err)
	}
}

func TestContactForces(t *testing.T) {
	points := []myb.ContactPoint{{Offset: [3]float64{0.5, 0, 0}, Radius: 0.5}, {Offset: [3]float64{-0.5, 0, 0}, Radius: 0.5}}
	attitude := [9]float64{0.6, 0.8, 0, -0.8, 0.6, 0, 0, 0, 1}
	omega := [3]float64{0, 0.01, 0.02}
	samples := []myb.ContactSample{
		{X: [3]float64{0, 0, moonRadius + 0.5 - 1e-4}, V: [3]float64{0.01, 0, -0.001},
			Attitude: &attitude, Omega: &omega, Deflections: [][3]float64{{1e-5, 0, 0}, {0, 2e-6, 0}}},
		{X: [3]float64{0.3, 0, moonRadius + 0.5 - 3e-4}, V: [3]float64{-0.02, 0.01, 0.002},
			Attitude: &attitude, Omega: &omega, Deflections: [][3]float64{{-4e-5, 1e-5, 0}, {0, 0, 0}}},
		{X: [3]float64{0, 0, moonRadius + 2}, V: [3]float64{0, 0, -1},
			Attitude: &attitude, Omega: &omega, Deflections: [][3]float64{{0, 0, 0}, {0, 0, 0}}},
	}
	var x, v, z []float64
	for _, s := range samples {
		x = append(x, s.X[:]...)
		v = append(v, s.V[:]...)
		z = append(z, flatten(s.Deflections)...)
	}
	values := func(results []myb.ContactResult) ([]float64, []float64) {
		var out, rates []float64
		for _, r := range results {
			out = append(out, r.Force[:]...)
			out = append(out, r.Torque[:]...)
			out = append(out, r.Energy, r.Gap, float64(r.Touching))
			rates = append(rates, flatten(r.DeflectionRates)...)
		}
		return out, rates
	}
	for _, law := range []myb.ContactLaw{stickLaw, smoothLaw} {
		results, err := myb.ContactForces(plane, law, points, samples)
		must(t, err)
		got, gotRates := values(results)
		status, want, rates, _ := cref.ContactForces(plane.Point[:], plane.Normal[:], lawValues(law),
			int(law.Kind), pointValues(points), x, v, attitude[:], omega[:], z)
		okStatus(t, "ContactForces", status)
		sameFloats(t, "results", got, want)
		if law.Kind == myb.FrictionStick {
			sameFloats(t, "rates", gotRates, rates)
		} else if gotRates != nil {
			t.Errorf("rates of a smooth friction")
		}
		if results[0].Touching != 2 || results[2].Touching != 0 {
			t.Errorf("touching %d, %d", results[0].Touching, results[2].Touching)
		}
	}

	// errors: a deflection that is not finite (the library, at the sample of C), not one
	// deflection a point (refused here), a friction that sticks without its state in C
	bad := append([]myb.ContactSample{}, samples...)
	bad[1].Deflections = [][3]float64{{math.NaN(), 0, 0}, {0, 0, 0}}
	_, err := myb.ContactForces(plane, stickLaw, points, bad)
	z[6] = math.NaN()
	status, _, _, failed := cref.ContactForces(plane.Point[:], plane.Normal[:], lawValues(stickLaw),
		int(myb.FrictionStick), pointValues(points), x, v, attitude[:], omega[:], z)
	var me *myb.Error
	if !errors.As(err, &me) || int(me.Status) != status || status == 0 || me.Index != failed || failed != 1 {
		t.Errorf("a deflection not finite: %v, the C status %d at %d", err, status, failed)
	}
	bad[1].Deflections = bad[1].Deflections[:1]
	if _, err := myb.ContactForces(plane, stickLaw, points, bad); !errors.As(err, &me) || me.Status != myb.Domain || me.Index != 1 {
		t.Errorf("one deflection for two points: %v", err)
	}
	if _, err := myb.ContactForces(plane, stickLaw, nil, nil); err != nil {
		t.Errorf("no samples: %v", err)
	}
	alone := []myb.ContactSample{{X: samples[2].X, V: samples[2].V}}
	if _, err := myb.ContactForces(plane, stickLaw, nil, alone); !errors.As(err, &me) || me.Status != myb.Domain || me.Index != 0 || me.Message != "" {
		t.Errorf("no points: %v", err)
	}
}

// ---- PN pairs: the new spin and tail terms ---------------------------------------------------

func TestPNSpinsAndTails(t *testing.T) {
	gms := []float64{1.327e21, 6.6e20}
	const c, r0 = 299792458.0, 1e6
	all := myb.PN1 | myb.PN2 | myb.PN25 | myb.PN3 | myb.PN35 | myb.PNSO | myb.PNSS | myb.PN4 |
		myb.PN4RR | myb.PNSO25 | myb.PNSO35 | myb.PNSORR | myb.PNSS3 | myb.PNSORR5 | myb.PNSOTAIL |
		myb.PNTAILX
	// a circular orbit and a hyperbolic flyby: the tails of the second need PNTAILX
	states := [][]myb.State{
		{{0, 0, 0, 0, 0, 0}, {1e6, 0, 0, 0, 4.46e7, 1e6}},
		{{0, 0, 0, 0, 0, 0}, {1e6, 2e5, 0, -1e7, 8e7, 0}},
	}
	chi1 := [][3]float64{{0.1, 0.2, 0.5}, {0.3, -0.1, 0.2}}
	chi2 := [][3]float64{{-0.3, 0.1, 0.4}, {0, 0.5, -0.2}}
	accelerations := func(orders myb.PNOrders) ([]float64, error) {
		system, err := myb.NewSystem(myb.SystemSpec{
			Bodies: []myb.SystemBody{{GM: gms[0]}, {GM: gms[1]}},
			Links:  []myb.Link{myb.PNPair{A: 0, B: 1, Orders: orders, C: c, R0: r0, Kappa1: 1, Kappa2: 1.5}},
		})
		if err != nil {
			return nil, err
		}
		defer system.Close()
		var ys [][]float64
		for k, s := range states {
			y, err := system.PackState(s, nil, []myb.SpinState{{Pair: 0, Chi1: chi1[k], Chi2: chi2[k]}})
			if err != nil {
				return nil, err
			}
			ys = append(ys, y)
		}
		a, err := system.Accelerations([]float64{0, 0}, ys)
		if err != nil {
			return nil, err
		}
		var flat []float64
		for _, v := range a {
			flat = append(flat, v...)
		}
		return flat, nil
	}
	var bodies []float64
	for _, s := range states {
		bodies = append(bodies, flatten(s)...)
	}
	got, err := accelerations(all)
	must(t, err)
	status, want, _ := cref.PNBinary(gms, uint(all), c, r0, 1, 1.5, bodies, flatten(chi1), flatten(chi2))
	okStatus(t, "PNBinary", status)
	sameFloats(t, "accelerations", got, want)
	older, err := accelerations(all &^ (myb.PNSS3 | myb.PNSORR5 | myb.PNSOTAIL))
	must(t, err)
	if same := func() bool {
		for k := range got {
			if got[k] != older[k] {
				return false
			}
		}
		return true
	}(); same {
		t.Errorf("the new spin terms changed nothing")
	}

	// errors: the flyby without PNTAILX, at the state of C
	_, err = accelerations(all &^ myb.PNTAILX)
	status, _, failed := cref.PNBinary(gms, uint(all&^myb.PNTAILX), c, r0, 1, 1.5, bodies,
		flatten(chi1), flatten(chi2))
	var me *myb.Error
	if !errors.As(err, &me) || int(me.Status) != status || status == 0 || me.Index != failed || failed != 1 {
		t.Errorf("a flyby without the tails on any conic: %v, the C status %d at %d", err, status, failed)
	}
}
