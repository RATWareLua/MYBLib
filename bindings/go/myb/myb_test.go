// SPDX-License-Identifier: BSD-3-Clause

// The Go binding against the C library it wraps:
//
//   - the low orbit of examples/low_orbit.c built from the types of the package gives the
//     nodes, the states there, the final state and the steps of the same calls made in C
//     (the program tests/bindings/low_orbit_reference.c, named by MYB_GO_REFERENCE) bit for
//     bit — the nodes found both by a ready-made event and by a Go function;
//   - the forms over slices give the bytes of the single calls made through the package:
//     time scales, an SPK set, Kepler's problem, the rotations of the Earth and of other
//     bodies, the accelerations of a model;
//   - events, a fit of an orbit with a parameter, OEM written and read in both notations, SPK
//     written and read back;
//   - failures are error values with the status of the library, the index of the element and
//     the line of a parser; objects closed and objects kept alive by a model; a Go term that
//     fails; a panic in a Go function that C called ends the process.
//
// MYB_DATA_DIR and MYB_TEST_DATA_DIR point at data/ and tests/data/ of the repository (found
// from this file when they are not set).
package myb_test

import (
	"bufio"
	"bytes"
	"errors"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
	"unsafe"

	"github.com/RATWareLua/MYBLib/bindings/go/myb"
)

const (
	gm     = 3.986004418e14
	radius = 6378137.0
	day    = 86400.0
)

func repository() string {
	return filepath.Join("..", "..", "..")
}

func dataPath(relative string) string {
	if dir := os.Getenv("MYB_DATA_DIR"); dir != "" {
		return filepath.Join(dir, relative)
	}
	return filepath.Join(repository(), "data", relative)
}

func testDataPath(relative string) string {
	if dir := os.Getenv("MYB_TEST_DATA_DIR"); dir != "" {
		return filepath.Join(dir, relative)
	}
	return filepath.Join(repository(), "tests", "data", relative)
}

// bytesOf is the memory of a slice: two slices with the same bytes hold the same bits.
func bytesOf[T any](s []T) []byte {
	if len(s) == 0 {
		return nil
	}
	var zero T
	return unsafe.Slice((*byte)(unsafe.Pointer(&s[0])), len(s)*int(unsafe.Sizeof(zero)))
}

func sameBytes[T any](t *testing.T, what string, a, b []T) {
	t.Helper()
	if len(a) != len(b) || !bytes.Equal(bytesOf(a), bytesOf(b)) {
		t.Errorf("%s: the slice form differs from the single calls", what)
	}
}

func must(t *testing.T, err error) {
	t.Helper()
	if err != nil {
		t.Fatal(err)
	}
}

// ---- The low orbit, bit for bit ------------------------------------------------------------

type reference struct {
	inputs map[string]float64
	nodes  []myb.Event
	final  myb.State
	steps  int64
}

func bitsOf(t *testing.T, word string) float64 {
	t.Helper()
	u, err := strconv.ParseUint(word, 16, 64)
	if err != nil {
		t.Fatalf("the reference: %v", err)
	}
	return math.Float64frombits(u)
}

// runReference runs the C program and reads the bits it prints.
func runReference(t *testing.T) reference {
	t.Helper()
	program := os.Getenv("MYB_GO_REFERENCE")
	if program == "" {
		t.Skip("MYB_GO_REFERENCE names no C reference program (ctest sets it)")
	}
	out, err := exec.Command(program).Output()
	if err != nil {
		t.Fatalf("%s: %v", program, err)
	}
	ref := reference{inputs: map[string]float64{}}
	scanner := bufio.NewScanner(bytes.NewReader(out))
	for scanner.Scan() {
		f := strings.Fields(scanner.Text())
		switch {
		case len(f) == 3 && f[0] == "input":
			ref.inputs[f[1]] = bitsOf(t, f[2])
		case len(f) == 8 && f[0] == "node":
			e := myb.Event{T: bitsOf(t, f[1])}
			for i := range 6 {
				e.State[i] = bitsOf(t, f[2+i])
			}
			ref.nodes = append(ref.nodes, e)
		case len(f) == 7 && f[0] == "final":
			for i := range 6 {
				ref.final[i] = bitsOf(t, f[1+i])
			}
		case len(f) == 2 && f[0] == "steps":
			ref.steps, err = strconv.ParseInt(f[1], 10, 64)
			must(t, err)
		default:
			t.Fatalf("the reference printed %q", scanner.Text())
		}
	}
	if len(ref.nodes) < 14 || ref.steps == 0 || len(ref.inputs) != 7 {
		t.Fatalf("the reference is incomplete: %d nodes, %d steps", len(ref.nodes), ref.steps)
	}
	return ref
}

func lowOrbitModel(t *testing.T, c20 float64) *myb.Model {
	t.Helper()
	field, err := myb.NewGravityField(gm, radius, 2, []float64{1, 0, 0, c20, 0, 0}, make([]float64, 6))
	must(t, err)
	model, err := myb.NewModel(myb.ModelSpec{
		Bodies: []myb.Body{{ID: 1, GM: gm}},
		Terms: []myb.Term{
			myb.Field{Name: "gravity", Body: 0, Field: field, NMax: -1},
			myb.Drag{Name: "drag", Body: 0, Object: myb.Object{Mass: 500, Area: 4, Cr: 1.3, Cd: 2.2},
				EquatorialRadius: radius, Layers: []myb.Layer{{Base: 400000, Density: 3e-12, Scale: 58000}}},
		},
	})
	must(t, err)
	must(t, field.Close()) // the model keeps the field
	return model
}

func sameEvents(t *testing.T, what string, got, want []myb.Event) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("%s: %d nodes, the C calls %d", what, len(got), len(want))
	}
	for k := range got {
		if math.Float64bits(got[k].T) != math.Float64bits(want[k].T) || got[k].State != want[k].State {
			t.Errorf("%s: node %d at %v, the C calls %v", what, k, got[k].T, want[k].T)
		}
	}
}

func TestLowOrbitBitForBit(t *testing.T) {
	ref := runReference(t)
	model := lowOrbitModel(t, ref.inputs["c20"])
	defer model.Close()
	if names := model.Terms(); len(names) != 2 || names[1] != "drag" {
		t.Errorf("terms %v", names)
	}
	y0, err := myb.StateFromElements(gm, myb.Elements{P: ref.inputs["p"], E: ref.inputs["e"],
		I: ref.inputs["i"], RAAN: ref.inputs["raan"], ArgP: ref.inputs["argp"], Nu: ref.inputs["nu"]})
	must(t, err)

	// the nodes as data: the plane of the equator through the origin
	nodes, err := model.Events(myb.EventSearch{T0: 0, Y0: y0, TEnd: day, Events: []myb.EventSpec{
		myb.PlaneCrossing{Body: myb.Origin, Normal: [3]float64{0, 0, 1}, Direction: myb.Rising, Limit: 32},
	}})
	must(t, err)
	sameEvents(t, "PlaneCrossing", nodes.Events, ref.nodes)
	if nodes.Final != ref.final || nodes.Steps != ref.steps {
		t.Errorf("final state or steps (%d) differ from the C calls (%d)", nodes.Steps, ref.steps)
	}

	// the nodes by a Go function, the event of the example itself
	byGo, err := model.Events(myb.EventSearch{T0: 0, Y0: y0, TEnd: day, Events: []myb.EventSpec{
		myb.GoEvent{Func: func(_ float64, y myb.State) float64 { return y[2] }, Direction: myb.Rising, Limit: 32},
	}})
	must(t, err)
	sameEvents(t, "GoEvent", byGo.Events, ref.nodes)
	if byGo.Final != ref.final || byGo.Steps != ref.steps {
		t.Errorf("final state or steps differ from the C calls (GoEvent)")
	}

	// the state a day later by a propagation, and the states at the nodes
	times := make([]float64, 0, len(ref.nodes)+1)
	for _, n := range ref.nodes {
		times = append(times, n.T)
	}
	times = append(times, day)
	trajectory, err := model.Propagate(myb.Propagation{T0: 0, Y0: y0, TEnd: day, Times: times})
	must(t, err)
	if trajectory.States[len(times)-1] != ref.final || trajectory.Steps != ref.steps {
		t.Errorf("the propagated state differs from the C calls")
	}
	for k, n := range ref.nodes {
		if trajectory.States[k] != n.State {
			t.Errorf("the state at node %d differs from the C calls", k)
		}
	}
}

// ---- Slices against single calls ---------------------------------------------------------

func timeSystem(t *testing.T, extra ...string) (*myb.KernelPool, *myb.TimeSystem) {
	t.Helper()
	pool, err := myb.NewKernelPool()
	must(t, err)
	for _, kernel := range append([]string{"time/time-scales.tk", "time/tdb-tt-fb1990.tk", "time/utc-tai.tk"}, extra...) {
		must(t, pool.LoadFile(dataPath(kernel)))
	}
	system, err := myb.NewTimeSystem(pool)
	must(t, err)
	return pool, system
}

func TestTimeSlices(t *testing.T) {
	pool, system := timeSystem(t)
	defer system.Close()
	must(t, pool.Close()) // the system keeps nothing of the pool
	seconds := []float64{-1.5e9, -3.25, 0, 0.5, 6.4e8 + 0.125, 7.1e8 + 1.0/3.0, 7.9e8}
	many, err := myb.TimesFromSeconds(seconds)
	must(t, err)
	single := make([]myb.Time, len(seconds))
	for k, s := range seconds {
		single[k], err = myb.TimeFromSeconds(s)
		must(t, err)
	}
	sameBytes(t, "TimesFromSeconds", many, single)
	back, err := myb.SecondsOf(many)
	must(t, err)
	for k := range many {
		s, err := myb.Seconds(many[k])
		must(t, err)
		if math.Float64bits(s) != math.Float64bits(back[k]) {
			t.Errorf("SecondsOf %d", k)
		}
	}
	labels, err := myb.TimesToCalendars(many)
	must(t, err)
	fromLabels, err := myb.CalendarsToTimes(labels)
	must(t, err)
	for k := range many {
		c, err := myb.TimeToCalendar(many[k])
		must(t, err)
		u, err := myb.CalendarToTime(c)
		must(t, err)
		if c != labels[k] || u != fromLabels[k] {
			t.Errorf("calendar %d", k)
		}
	}
	for _, pair := range [][2]string{{"TT", "TDB"}, {"TDB", "TT"}, {"TT", "TAI"}, {"TDB", "TCB"}} {
		converted, err := system.ConvertMany(many, pair[0], pair[1])
		must(t, err)
		for k := range many {
			one, err := system.Convert(many[k], pair[0], pair[1])
			must(t, err)
			if one != converted[k] {
				t.Errorf("Convert %s -> %s at %d", pair[0], pair[1], k)
			}
		}
	}
	recent := many[3:]
	utc, err := system.ToUTCMany(recent, "TT")
	must(t, err)
	tdb, err := system.FromUTCMany(utc, "TDB")
	must(t, err)
	for k := range recent {
		label, err := system.ToUTC(recent[k], "TT")
		must(t, err)
		instant, err := system.FromUTC(label, "TDB")
		must(t, err)
		if label != utc[k] || instant != tdb[k] {
			t.Errorf("UTC %d", k)
		}
	}
	// failures: the first label that fails names its index
	bad := append([]myb.Calendar(nil), utc...)
	bad[2].Month = 13
	_, err = system.FromUTCMany(bad, "TT")
	var e *myb.Error
	if !errors.As(err, &e) || e.Status != myb.Domain || e.Index != 2 {
		t.Errorf("FromUTCMany with a wrong label: %v", err)
	}
	if _, err := system.ConvertMany(many, "TT", "NO SUCH SCALE"); !errors.Is(err, myb.ErrNotFound) {
		t.Errorf("an unknown scale: %v", err)
	}
	if out, err := system.ConvertMany(nil, "TT", "TDB"); err != nil || len(out) != 0 {
		t.Errorf("an empty slice: %v", err)
	}
}

func TestEphemerisSlices(t *testing.T) {
	eph, err := myb.NewEphemeris()
	must(t, err)
	defer eph.Close()
	must(t, eph.LoadFile(testDataPath("spk/t02_03.bsp")))
	segments, err := eph.Segments()
	must(t, err)
	s := segments[0]
	seconds := make([]float64, 9)
	for k := range seconds {
		seconds[k] = s.Start + (s.Stop-s.Start)*(float64(k)+0.5)/float64(len(seconds))
	}
	epochs, err := myb.TimesFromSeconds(seconds)
	must(t, err)
	states, frames, err := eph.States(s.Target, s.Center, epochs)
	must(t, err)
	for k, epoch := range epochs {
		state, frame, err := eph.State(s.Target, s.Center, epoch)
		must(t, err)
		if state != states[k] || frame != frames[k] {
			t.Errorf("Ephemeris.States %d", k)
		}
	}
	// beyond the coverage at index 4
	late, err := myb.TimeFromSeconds(s.Stop + 1e9)
	must(t, err)
	epochs[4] = late
	_, _, err = eph.States(s.Target, s.Center, epochs)
	var e *myb.Error
	if !errors.As(err, &e) || e.Status != myb.Range || e.Index != 4 {
		t.Errorf("an epoch beyond the coverage: %v", err)
	}
}

func TestKeplerSlices(t *testing.T) {
	s0 := myb.State{6.9e6, 1.2e5, -3.4e5, -150, 7.3e3, 1.1e3}
	dt := []float64{-2e4, -1.5, 0, 600, 5400.25, 3e4, 1e6}
	many, err := myb.KeplerPropagateMany(gm, s0, dt)
	must(t, err)
	withSTM, stms, err := myb.KeplerSTMMany(gm, s0, dt)
	must(t, err)
	for k := range dt {
		one, err := myb.KeplerPropagate(gm, s0, dt[k])
		must(t, err)
		two, stm, err := myb.KeplerSTM(gm, s0, dt[k])
		must(t, err)
		if one != many[k] || two != withSTM[k] || stm != stms[k] {
			t.Errorf("Kepler %d", k)
		}
	}
	elements, err := myb.ElementsFromStates(gm, many)
	must(t, err)
	states, err := myb.StatesFromElements(gm, elements)
	must(t, err)
	equinoctial, err := myb.EquinoctialFromStates(gm, many)
	must(t, err)
	fromEquinoctial, err := myb.StatesFromEquinoctial(gm, equinoctial)
	must(t, err)
	for k := range many {
		e, err := myb.ElementsFromState(gm, many[k])
		must(t, err)
		s, err := myb.StateFromElements(gm, e)
		must(t, err)
		q, err := myb.EquinoctialFromState(gm, many[k])
		must(t, err)
		r, err := myb.StateFromEquinoctial(gm, q)
		must(t, err)
		if e != elements[k] || s != states[k] || q != equinoctial[k] || r != fromEquinoctial[k] {
			t.Errorf("elements %d", k)
		}
	}
	dt[2] = math.NaN()
	_, err = myb.KeplerPropagateMany(gm, s0, dt)
	var e *myb.Error
	if !errors.As(err, &e) || e.Status != myb.Domain || e.Index != 2 || e.Op != "KeplerPropagateMany" {
		t.Errorf("a time that is NaN: %v", err)
	}
}

func TestRotationSlices(t *testing.T) {
	pool, system := timeSystem(t, "frames/iau2006-cip.tk", "frames/iau1980.tk", "frames/eop-tides.tk")
	defer system.Close()
	defer pool.Close()
	cip, err := myb.NewCIP(pool, "IAU2006A")
	must(t, err)
	defer cip.Close()
	equinox, err := myb.NewEquinox(pool, "IAU1980")
	must(t, err)
	defer equinox.Close()
	eop, err := myb.NewEOP()
	must(t, err)
	defer eop.Close()
	must(t, eop.LoadC04File(testDataPath("eop/eopc04-excerpt.txt")))
	must(t, eop.SetTides(pool, "IERS10_EOPT"))
	seconds := make([]float64, 10)
	for k := range seconds {
		seconds[k] = 6.32e8 + 2.5e6*float64(k) + 0.375
	}
	tt, err := myb.TimesFromSeconds(seconds)
	must(t, err)

	xys, rates, err := cip.XYSMany(tt)
	must(t, err)
	values, ut1, err := eop.AtMany(system, tt, "TT")
	must(t, err)
	terrestrial, err := cip.TerrestrialMany(tt, ut1, values)
	must(t, err)
	tod, err := equinox.RotationMany(myb.TOD, tt)
	must(t, err)
	pef, err := equinox.TEMEPEFMany(ut1)
	must(t, err)
	earth := myb.EarthOrientation{CIP: cip, EOP: eop, TimeSystem: system, TT: "TT"}
	sampled, err := earth.Sample(tt)
	must(t, err)
	state := myb.State{6.5e6, 1.0e6, -2.0e6, 1.0e3, 7.0e3, -1.5e3}
	repeated := make([]myb.State, len(tt))
	for k := range repeated {
		repeated[k] = state
	}
	carried, err := myb.ApplyRotations(terrestrial, repeated)
	must(t, err)
	for k := range tt {
		x, r, err := cip.XYS(tt[k])
		must(t, err)
		v, u, err := eop.At(system, tt[k], "TT")
		must(t, err)
		m, err := cip.Terrestrial(tt[k], u, &v)
		must(t, err)
		q, err := equinox.Rotation(myb.TOD, tt[k])
		must(t, err)
		p, err := equinox.TEMEPEF(u)
		must(t, err)
		o, err := earth.Sample(tt[k : k+1])
		must(t, err)
		if x != xys[k] || r != rates[k] || v != values[k] || u != ut1[k] || m != terrestrial[k] ||
			q != tod[k] || p != pef[k] || o[0] != sampled[k] || m.Apply(state) != carried[k] {
			t.Errorf("rotations of the Earth at %d", k)
		}
	}
	if sampled[3].Rotation != terrestrial[3] || sampled[3].UT != ut1[3] {
		t.Errorf("the orientation of the Earth is not its terrestrial rotation")
	}
	// an instant before the series at index 6
	early, err := myb.TimeFromSeconds(3.0e8)
	must(t, err)
	tt[6] = early
	_, _, err = eop.AtMany(system, tt, "TT")
	var e *myb.Error
	if !errors.As(err, &e) || e.Status != myb.Range || e.Index != 6 {
		t.Errorf("an instant outside the EOP: %v", err)
	}
	if _, err := earth.Sample(tt); !errors.Is(err, myb.ErrRange) {
		t.Errorf("an orientation outside the EOP: %v", err)
	}

	// other bodies: Mars by the WGCCRE, a frame of a binary PCK
	planets, err := myb.NewKernelPool()
	must(t, err)
	defer planets.Close()
	must(t, planets.LoadFile(testDataPath("wgccre/pck00011.tpc")))
	mars, err := myb.NewWGCCRE(planets, 499)
	must(t, err)
	defer mars.Close()
	marsRotations, err := mars.RotationMany(tt[:5])
	must(t, err)
	marsSampled, err := myb.WGCCREOrientation{Model: mars}.Sample(tt[:5])
	must(t, err)
	for k := range marsRotations {
		r, err := mars.Rotation(tt[k])
		must(t, err)
		if r != marsRotations[k] || r != marsSampled[k].Rotation {
			t.Errorf("Mars at %d", k)
		}
	}
	pck, err := myb.NewPCK()
	must(t, err)
	defer pck.Close()
	must(t, pck.LoadFile(testDataPath("frames/graph.bpc")))
	pckSeconds := []float64{600002748, 600003000.25, 600004000.5, 600004677}
	pckTimes, err := myb.TimesFromSeconds(pckSeconds)
	must(t, err)
	frames, bases, err := pck.RotationMany(31099, pckTimes)
	must(t, err)
	for k := range frames {
		r, base, err := pck.Rotation(31099, pckTimes[k])
		must(t, err)
		if r != frames[k] || base != bases[k] {
			t.Errorf("the PCK at %d", k)
		}
	}
}

// ---- Models: accelerations, events, Go terms, a fit ------------------------------------------

func TestAccelerationsAndGoTerms(t *testing.T) {
	push := [3]float64{1e-6, -2e-6, 3e-6}
	calls := 0
	model, err := myb.NewModel(myb.ModelSpec{Terms: []myb.Term{
		myb.PointMass{Name: "planet", GM: gm},
		myb.GoTerm{Name: "push", Func: func(_ float64, _, _ [3]float64) ([3]float64, error) {
			calls++
			return push, nil
		}},
	}})
	must(t, err)
	defer model.Close()
	points := []myb.State{{7e6, 1e5, -2e5, 0, 7.5e3, 1e3}, {-6.8e6, 2e6, 1e5, -1e3, -7e3, 5e2}}
	times := []float64{0, 100}
	many, err := model.Accelerations(times, points)
	must(t, err)
	for k := range points {
		one, err := model.Acceleration(times[k], points[k])
		must(t, err)
		if one != many[k] {
			t.Errorf("Accelerations %d", k)
		}
	}
	if calls != 4 {
		t.Errorf("the Go term was called %d times", calls)
	}
	alone, err := myb.NewModel(myb.ModelSpec{Terms: []myb.Term{myb.PointMass{Name: "planet", GM: gm}}})
	must(t, err)
	defer alone.Close()
	a, err := alone.Acceleration(0, points[0])
	must(t, err)
	for i := range 3 {
		if math.Abs(many[0][i]-a[i]-push[i]) > 1e-14 {
			t.Errorf("the push is not added: %v %v", many[0], a)
		}
	}
	// an orbit with the Go term
	trajectory, err := model.Propagate(myb.Propagation{T0: 0, Y0: points[0], TEnd: 3000, Times: []float64{1000, 3000}})
	must(t, err)
	if len(trajectory.States) != 2 || calls < 100 {
		t.Errorf("a propagation with a Go term: %d calls", calls)
	}

	// a Go term that fails: its error comes back wrapped, with its status
	refused := errors.New("refused")
	failing, err := myb.NewModel(myb.ModelSpec{Terms: []myb.Term{
		myb.GoTerm{Name: "refuses", Func: func(float64, [3]float64, [3]float64) ([3]float64, error) {
			return [3]float64{}, refused
		}},
	}})
	must(t, err)
	defer failing.Close()
	_, err = failing.Acceleration(0, points[0])
	var e *myb.Error
	if !errors.Is(err, refused) || !errors.As(err, &e) || e.Status != myb.Domain {
		t.Errorf("a failing Go term: %v", err)
	}
	ranged, err := myb.NewModel(myb.ModelSpec{Terms: []myb.Term{
		myb.GoTerm{Name: "range", Func: func(float64, [3]float64, [3]float64) ([3]float64, error) {
			return [3]float64{}, myb.ErrRange
		}},
	}})
	must(t, err)
	defer ranged.Close()
	if _, err := ranged.Propagate(myb.Propagation{T0: 0, Y0: points[0], TEnd: 10, Times: []float64{5}}); !errors.Is(err, myb.ErrRange) {
		t.Errorf("a Go term failing with MYB_E_RANGE: %v", err)
	}
}

func TestEvents(t *testing.T) {
	model, err := myb.NewModel(myb.ModelSpec{Bodies: []myb.Body{{ID: 399, GM: gm}},
		Terms: []myb.Term{myb.PointMass{Name: "planet", GM: gm}}})
	must(t, err)
	defer model.Close()
	y0, err := myb.StateFromElements(gm, myb.Elements{P: 8e6, E: 0.1, I: 0.7, RAAN: 0.3, ArgP: 0.5, Nu: 0.2})
	must(t, err)
	result, err := model.Events(myb.EventSearch{T0: 0, Y0: y0, TEnd: 30000, Events: []myb.EventSpec{
		myb.Apsis{Body: myb.Origin, Direction: myb.Rising},
		myb.Apsis{Body: myb.Origin, Direction: myb.Falling},
		myb.PlaneCrossing{Body: 0, Normal: [3]float64{0, 0, 1}, Direction: myb.Any},
		myb.Distance{Body: myb.Origin, Range: 8.5e6, Direction: myb.Falling},
	}})
	must(t, err)
	count := map[int]int{}
	for k, e := range result.Events {
		count[e.Index]++
		if k > 0 && e.T < result.Events[k-1].T {
			t.Errorf("the events are not in the order of time")
		}
	}
	// a period of 7.4e3 s: four pericentres and apocentres, eight nodes in 30000 s
	if count[0] < 3 || count[1] < 3 || count[2] < 7 || count[3] < 3 {
		t.Errorf("events found: %v", count)
	}
	// the Kepler orbit agrees at a pericentre
	for _, e := range result.Events {
		if e.Index == 0 {
			kepler, err := myb.KeplerPropagate(gm, y0, e.T)
			must(t, err)
			for i := range 3 {
				if math.Abs(kepler[i]-e.State[i]) > 1e-3 {
					t.Errorf("the pericentre at %v is not on the orbit", e.T)
				}
			}
			break
		}
	}
	if _, err := model.Events(myb.EventSearch{T0: 0, Y0: y0, TEnd: 100, Events: []myb.EventSpec{
		myb.Altitude{Body: 7}}}); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("an event of a body the model does not have: %v", err)
	}
}

func TestFit(t *testing.T) {
	truth := 2e-7
	spec := func(constant float64) myb.ModelSpec {
		return myb.ModelSpec{
			Terms: []myb.Term{
				myb.PointMass{Name: "planet", GM: gm},
				myb.Empirical{Name: "empirical", Constant: [3]float64{0, constant, 0}},
			},
			Parameters: []myb.Parameter{{Name: "along", Term: "empirical", Kind: myb.ParameterConstant, Index: 1, Step: 1e-9}},
		}
	}
	exact, err := myb.NewModel(spec(truth))
	must(t, err)
	defer exact.Close()
	y0 := myb.State{7e6, 0, 0, 0, 6.5e3, 3.8e3}
	times := make([]float64, 40)
	for k := range times {
		times[k] = 300 * float64(k+1)
	}
	observed, err := exact.Propagate(myb.Propagation{T0: 0, Y0: y0, TEnd: 12001, Times: times})
	must(t, err)
	// positions with errors of a centimetre in a fixed pattern (noise-free data would never let
	// the weighted RMS settle, which is what the fit takes as converged)
	observations := make([]myb.Observation, len(times))
	for k := range times {
		value := make([]float64, 3)
		for i := range 3 {
			value[i] = observed.States[k][i] + 0.005*float64((7*k+3*i)%5-2)
		}
		observations[k] = myb.Observation{T: times[k], Kind: myb.MeasurePosition, Value: value,
			Sigma: []float64{0.01, 0.01, 0.01}}
	}
	model, err := myb.NewModel(spec(0))
	must(t, err)
	defer model.Close()
	guess := y0
	guess[0] += 200
	guess[4] += 0.2
	result, err := model.Fit(myb.FitRequest{T0: 0, Y0: guess, Observations: observations,
		Covariance: true, Residuals: true})
	must(t, err)
	if !result.Converged || result.RMS > 1.5 || result.RMSBefore < 1e3 || result.Values != 120 {
		t.Errorf("fit: converged %v, rms %g (before %g), %d values", result.Converged, result.RMS,
			result.RMSBefore, result.Values)
	}
	for i := range 3 {
		if math.Abs(result.Y0[i]-y0[i]) > 0.05 {
			t.Errorf("the position at the epoch: %v against %v", result.Y0, y0)
		}
	}
	if math.Abs(result.Parameters[0]-truth) > 2e-9 {
		t.Errorf("the empirical acceleration: %g against %g", result.Parameters[0], truth)
	}
	now, err := model.Parameters()
	must(t, err)
	if now[0] != result.Parameters[0] || len(result.Covariance) != 49 || len(result.Residuals) != 120 {
		t.Errorf("the model does not keep the estimate")
	}
	observations[3].Value = []float64{1, 2}
	if _, err := model.Fit(myb.FitRequest{Y0: guess, Observations: observations}); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("an observation with too few values: %v", err)
	}
}

// ---- Files ------------------------------------------------------------------------------

func TestOEM(t *testing.T) {
	segment := myb.OEMSegment{
		ObjectName: "A & B", ObjectID: "2026-001A", CenterName: "EARTH", RefFrame: "EME2000",
		TimeSystem: "UTC", Interpolation: "HERMITE", InterpolationDegree: 7,
		Epochs: []myb.Calendar{
			{Year: 2026, Month: 3, Day: 1, Hour: 12},
			{Year: 2026, Month: 3, Day: 1, Hour: 12, Minute: 1, Second: 0.5},
		},
		States: []myb.State{{7e6, 1, -2.5, 0.25, 7.5e3, -1e3}, {7e6 + 15, 4.5e5, -6e4, -480, 7.4e3, -990}},
	}
	created := myb.Calendar{Year: 2026, Month: 3, Day: 2}
	for _, notation := range []myb.Notation{myb.KVN, myb.XML} {
		text, err := myb.WriteOEM("myb", created, []myb.OEMSegment{segment}, notation)
		must(t, err)
		if strings.HasPrefix(text, "<") != (notation == myb.XML) {
			t.Errorf("notation %d", notation)
		}
		message, err := myb.ParseOEM([]byte(text))
		must(t, err)
		if len(message.Segments) != 1 || message.Version != "2.0" || message.Originator != "myb" {
			t.Fatalf("the message read: %+v", message)
		}
		s := message.Segments[0]
		if s.ObjectName != "A & B" || s.InterpolationDegree != 7 || s.Meta["REF_FRAME"] != "EME2000" ||
			s.Accelerations != nil {
			t.Errorf("metadata %+v", s.Meta)
		}
		sameBytes(t, "OEM states", s.States, segment.States)
		sameBytes(t, "OEM epochs", s.Epochs, segment.Epochs)
	}
	path := filepath.Join(t.TempDir(), "orbit.oem")
	must(t, myb.WriteOEMFile(path, "myb", created, []myb.OEMSegment{segment}, myb.KVN))
	fromFile, err := myb.ReadOEMFile(path)
	must(t, err)
	sameBytes(t, "OEM file", fromFile.Segments[0].States, segment.States)

	_, err = myb.ParseOEM([]byte("CCSDS_OEM_VERS = 2.0\nMETA_START\n"))
	var e *myb.Error
	if !errors.As(err, &e) || e.Status != myb.Format || e.Line != 2 {
		t.Errorf("a broken message: %v", err)
	}
	uneven := segment
	uneven.States = uneven.States[:1]
	if _, err := myb.WriteOEM("myb", created, []myb.OEMSegment{uneven}, myb.KVN); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("uneven epochs and states: %v", err)
	}
	if _, err := myb.WriteOEM("myb", created, nil, myb.KVN); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("no segments: %v", err)
	}
}

func TestSPK(t *testing.T) {
	s0 := myb.State{6.9e6, 1.2e5, -3.4e5, -150, 7.3e3, 1.1e3}
	epochs := make([]float64, 12)
	for k := range epochs {
		epochs[k] = 7e8 + 300*float64(k)
	}
	offsets := make([]float64, len(epochs))
	for k := range offsets {
		offsets[k] = 300 * float64(k)
	}
	states, err := myb.KeplerPropagateMany(gm, s0, offsets)
	must(t, err)
	segments := []myb.SPKSegment{{Target: -999, Center: 399, Frame: 1, Name: "TEST ORBIT", Type: 13,
		Degree: 5, Epochs: epochs, States: states}}
	path := filepath.Join(t.TempDir(), "orbit.bsp")
	must(t, myb.WriteSPK(path, "myb test", segments))
	image, err := myb.EncodeSPK("myb test", segments)
	must(t, err)
	fromFile, err := myb.NewEphemeris()
	must(t, err)
	defer fromFile.Close()
	must(t, fromFile.LoadFile(path))
	fromBytes, err := myb.NewEphemeris()
	must(t, err)
	defer fromBytes.Close()
	must(t, fromBytes.LoadBytes(image))
	times, err := myb.TimesFromSeconds(epochs)
	must(t, err)
	a, _, err := fromFile.States(-999, 399, times)
	must(t, err)
	b, _, err := fromBytes.States(-999, 399, times)
	must(t, err)
	sameBytes(t, "SPK from a file and from bytes", a, b)
	for k := range states {
		for i := range 6 {
			if math.Abs(a[k][i]-states[k][i]) > 1e-12*math.Abs(states[k][i])+1e-12 {
				t.Errorf("the state %d read back: %v against %v", k, a[k], states[k])
			}
		}
	}
	segments[0].Degree = 4 // Hermite needs an odd degree
	if err := myb.WriteSPK(path, "myb test", segments); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("an even Hermite degree: %v", err)
	} else if e := err.(*myb.Error); e.Index != 0 {
		t.Errorf("the segment refused: %d", e.Index)
	}
}

// ---- Failures and ownership ---------------------------------------------------------------

func TestErrors(t *testing.T) {
	_, err := myb.StateFromElements(gm, myb.Elements{P: -1})
	var e *myb.Error
	if !errors.As(err, &e) || e.Status != myb.Domain || e.Status.String() != "MYB_E_DOMAIN" ||
		!errors.Is(err, myb.ErrDomain) || errors.Is(err, myb.ErrRange) {
		t.Errorf("elements with p < 0: %v", err)
	}
	if _, err := myb.TimeFromSeconds(math.NaN()); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("seconds that are NaN: %v", err)
	}
	if _, err := myb.TimesFromSeconds([]float64{0, math.Inf(1)}); !errors.As(err, &e) || e.Index != 1 {
		t.Errorf("an infinite second: %v", err)
	}
	if _, err := myb.Seconds(myb.Time{S: 1 << 60}); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("an instant beyond 2^53 s: %v", err)
	}
	eph, err := myb.NewEphemeris()
	must(t, err)
	if err := eph.LoadFile(testDataPath("spk/missing.bsp")); !errors.Is(err, myb.ErrIO) {
		t.Errorf("a missing file: %v", err)
	}
	must(t, eph.LoadFile(testDataPath("spk/t02_03.bsp")))
	segments, err := eph.Segments()
	must(t, err)
	s := segments[0]
	middle, err := myb.TimeFromSeconds(s.Start/2 + s.Stop/2)
	must(t, err)
	there, _, err := eph.State(s.Target, s.Center, middle)
	must(t, err)

	// a model keeps the ephemeris it uses after its Close
	model, err := myb.NewModel(myb.ModelSpec{Epoch: middle, Bodies: []myb.Body{
		{ID: int64(s.Center), GM: gm},
		{ID: int64(s.Target), GM: 1, State: &myb.SPKBody{Ephemeris: eph, Target: s.Target, Center: s.Center}},
	}, Terms: []myb.Term{myb.PointMass{Name: "centre", GM: gm}, myb.ThirdBody{Name: "third", Body: 1}}})
	must(t, err)
	must(t, eph.Close())
	must(t, eph.Close()) // twice is harmless
	if _, _, err := eph.State(s.Target, s.Center, middle); !errors.Is(err, myb.ErrClosed) {
		t.Errorf("a closed ephemeris answered: %v", err)
	}
	far := myb.State{there[0] * 2, there[1] * 2, there[2]*2 + 1e7, 0, 0, 0}
	if _, err := model.Acceleration(0, far); err != nil {
		t.Errorf("the model lost its ephemeris: %v", err)
	}
	must(t, model.Close())
	if _, err := model.Acceleration(0, far); !errors.Is(err, myb.ErrClosed) {
		t.Errorf("a closed model answered: %v", err)
	}

	// descriptions the library refuses, and the record that failed
	_, err = myb.NewModel(myb.ModelSpec{Bodies: []myb.Body{{ID: 1, GM: gm}},
		Terms: []myb.Term{myb.PointMass{Name: "a", GM: gm}, myb.ThirdBody{Name: "nobody", Body: 3}}})
	if !errors.As(err, &e) || e.Status != myb.Domain || e.Index != 1 || e.Message != "term" {
		t.Errorf("a term of a body the model does not have: %v", err)
	}
	_, err = myb.NewModel(myb.ModelSpec{Terms: []myb.Term{myb.PointMass{Name: "same", GM: gm}, myb.PointMass{Name: "same", GM: gm}}})
	if !errors.As(err, &e) || e.Index != 1 {
		t.Errorf("two terms of one name: %v", err)
	}
	_, err = myb.NewModel(myb.ModelSpec{Terms: []myb.Term{myb.PointMass{Name: "a", GM: gm}},
		Parameters: []myb.Parameter{{Name: "cd", Term: "drag", Kind: myb.ParameterCd, Step: 0.01}}})
	if !errors.Is(err, myb.ErrNotFound) {
		t.Errorf("a parameter of a term the model does not have: %v", err)
	}
	_, err = myb.NewModel(myb.ModelSpec{Terms: []myb.Term{myb.PointMass{Name: "a", GM: gm}},
		Parameters: []myb.Parameter{{Name: "cd", Term: "a", Kind: myb.ParameterCd, Step: 0.01}}})
	if !errors.As(err, &e) || e.Message != "parameter" {
		t.Errorf("a drag coefficient of a point mass: %v", err)
	}
	if _, err := myb.NewGravityField(gm, radius, 4, []float64{1}, []float64{0}); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("a field without its coefficients: %v", err)
	}
	if _, err := myb.ApplyRotations(make([]myb.Rotation, 2), make([]myb.State, 3)); !errors.Is(err, myb.ErrDomain) {
		t.Errorf("slices of different lengths: %v", err)
	}
	if myb.Version() < 1020000 {
		t.Errorf("version %d", myb.Version())
	}
}

// A panic in a Go function that C calls is not recovered: the process ends (the test runs
// itself again as the process that panics).
func TestPanicInGoFunctionEndsTheProcess(t *testing.T) {
	if what := os.Getenv("MYB_GO_PANIC_IN"); what != "" {
		model, err := myb.NewModel(myb.ModelSpec{Terms: []myb.Term{
			myb.GoTerm{Name: "broken", Func: func(float64, [3]float64, [3]float64) ([3]float64, error) {
				panic("in a term")
			}},
		}})
		must(t, err)
		_, _ = model.Acceleration(0, myb.State{7e6, 0, 0, 0, 7.5e3, 0})
		os.Exit(0) // not reached: leaving normally is what the parent takes as a failure
	}
	command := exec.Command(os.Args[0], "-test.run=^TestPanicInGoFunctionEndsTheProcess$")
	command.Env = append(os.Environ(), "MYB_GO_PANIC_IN=term")
	output, err := command.CombinedOutput()
	if err == nil {
		t.Fatalf("the process went on after a panic in a term: %s", output)
	}
	if !strings.Contains(string(output), "panic: in a term") {
		t.Errorf("the process ended without the panic's message: %s", output)
	}
}
