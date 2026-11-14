// SPDX-License-Identifier: BSD-3-Clause

// The package of ext/sgp4 against the calls it makes: a set read and written back character
// for character, the theory over a slice of times with the bits of the single calls, the
// exchange with a mean-elements message, a fit, and the errors of bad input.
package sgp4_test

import (
	"errors"
	"math"
	"testing"

	"github.com/RATWareLua/MYBLib/bindings/go/myb"
	"github.com/RATWareLua/MYBLib/bindings/go/myb/sgp4"
)

const (
	line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
	line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"
)

func must(t *testing.T, err error) {
	t.Helper()
	if err != nil {
		t.Fatal(err)
	}
}

func TestTLE(t *testing.T) {
	tle, err := sgp4.Parse(line1, line2, false)
	must(t, err)
	a, b, err := tle.Format()
	must(t, err)
	if a != line1 || b != line2 {
		t.Errorf("the set written back:\n%s\n%s", a, b)
	}
	if epoch := tle.Epoch(); epoch.Year != 2000 || epoch.Month != 6 {
		t.Errorf("the epoch %+v", epoch)
	}
	record, err := sgp4.NewRecord(tle, sgp4.WGS72, sgp4.Improved)
	must(t, err)
	minutes := []float64{0, 360, 720, 1440, 4320}
	states, codes, err := record.PropagateMany(minutes)
	must(t, err)
	for k, m := range minutes {
		s, code, err := record.Propagate(m)
		must(t, err)
		for i := range 6 {
			if math.Float64bits(s[i]) != math.Float64bits(states[k][i]) || code != codes[k] {
				t.Errorf("the state at %v minutes differs from the single call", m)
			}
		}
	}
	// the state at the epoch of the verification output (tcppver.out: 7022.46529266 km)
	if math.Abs(states[0][0]-7022465.29266) > 1e-2 {
		t.Errorf("x at the epoch %v", states[0][0])
	}
	omm, err := tle.OMM()
	must(t, err)
	omm.Originator, omm.ObjectName = "myb", "TEST"
	text, err := myb.WriteOMM(omm, myb.KVN)
	must(t, err)
	read, err := myb.ParseOMM([]byte(text))
	must(t, err)
	back, err := sgp4.FromOMM(read)
	must(t, err)
	if c, d, _ := back.Format(); c != line1 || d != line2 {
		t.Errorf("the set through an OMM:\n%s\n%s", c, d)
	}

	// a fit to the positions of the theory itself comes back to the set
	positions := make([][3]float64, len(minutes))
	for k := range minutes {
		positions[k] = [3]float64{states[k][0], states[k][1], states[k][2]}
	}
	guess, err := sgp4.FromState(tle, states[0], sgp4.WGS72)
	must(t, err)
	fitted, rms, _, err := sgp4.Fit(guess, minutes, positions, sgp4.FitOptions{Gravity: sgp4.WGS72})
	must(t, err)
	if rms > 1 || math.Abs(fitted.MeanMotion-tle.MeanMotion) > 1e-6 {
		t.Errorf("the fit: rms %v m, mean motion %v against %v", rms, fitted.MeanMotion, tle.MeanMotion)
	}

	// errors
	_, err = sgp4.Parse(line1[:60], line2, false)
	var e *myb.Error
	if !errors.As(err, &e) || e.Status != myb.Format {
		t.Errorf("a short line: %v", err)
	}
	bad := tle
	bad.Eccentricity = 1.5
	if _, _, err := bad.Format(); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("an eccentricity of 1.5: %v", err)
	}
	if _, err := sgp4.NewRecord(bad, sgp4.WGS72, sgp4.Improved); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a record of e = 1.5: %v", err)
	}
}
