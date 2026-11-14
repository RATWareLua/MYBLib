// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include <stdlib.h>
#include "glue.h"
*/
import "C"

import (
	"cmp"
	"runtime/cgo"
	"slices"
	"sync"
	"unsafe"
)

// ---- The description of a model ----------------------------------------------------------------

// ModelSpec describes a force model: the environment (the instant of t = 0 and its scale,
// the bodies) and the terms, with the parameters a fit estimates. NewModel turns it into the
// model of the library in one call; the objects it points to (ephemerides, fields, ...) are
// kept by the model, everything else is copied.
type ModelSpec struct {
	Epoch      Time        // the instant of the model time t = 0
	TimeSystem *TimeSystem // nil: one time scale, every provider gets Epoch + t
	Scale      string      // the scale of Epoch and t (needed with a time system), e.g. "TT"
	TableStep  float64     // s: the table of the differences of the scales; 0 = 3600
	Bodies     []Body
	Terms      []Term
	Parameters []Parameter
	// how IAS15 ends its steps at the breaks of the model (the edges of shadows of a Radiation
	// term with Breaks): BreaksRetake (the default) or BreaksPredict
	BreakScheme BreakScheme
}

// BreakScheme is how myb_orbit with IAS15 ends its steps at the breaks of a model
// (myb_break_scheme).
type BreakScheme int

// The schemes.
const (
	BreaksRetake  BreakScheme = C.MYB_BREAKS_RETAKE  // a step across a break is withdrawn and retaken
	BreaksPredict BreakScheme = C.MYB_BREAKS_PREDICT // the steps are aimed at the breaks predicted
)

// Body is a body of the environment, by its index in ModelSpec.Bodies. Everything is optional:
// zero means absent. The body without a State is the origin of the coordinates (one at most).
type Body struct {
	ID               int64
	GM               float64 // m^3/s^2
	Radius           float64 // m: the sphere that casts shadows
	Luminosity       float64 // W: the body as a source of light
	Limb             [2]float64
	Spin             [3]float64  // m^2/s, body-fixed
	State            StateSource // nil: the body is the origin
	StateScale       string      // the scale of the state's argument; "": the model's
	Orientation      OrientationSource
	OrientationScale string // "": the model's
}

// Object is the body whose motion is modelled, as the surface forces see it.
type Object struct {
	Mass float64 // kg
	Area float64 // m^2
	Cr   float64 // radiation pressure coefficient
	Cd   float64 // drag coefficient
}

// Ellipsoid is an oblate ellipsoid: equatorial radius A (m) and flattening F.
type Ellipsoid struct {
	A, F float64
}

// Layer is a layer of an exponential atmosphere: above the geodetic height Base (m) the
// density Density (kg/m^3) falls with the scale height Scale (m).
type Layer struct {
	Base, Density, Scale float64
}

// PPN are the speed of light (m/s) and the post-Newtonian parameters beta and gamma.
type PPN struct {
	C, Beta, Gamma float64
}

// Term is a term of a model; its name (at most 31 characters) is unique in the model.
type Term interface {
	describe(t *C.mybgo_term, k *keeper) error
}

// PointMass is the attraction of a point mass at the origin.
type PointMass struct {
	Name string
	GM   float64
}

// ThirdBody is the point mass of a body of the environment: with a body at the origin, the
// perturbation direct minus indirect.
type ThirdBody struct {
	Name string
	Body int
}

// Field is the gravity field of a body in its body-fixed frame, degrees NMin .. NMax (NMax
// -1: all). With Update the coefficients follow time (the model copies the field first).
type Field struct {
	Name       string
	Body       int
	Field      *GravityField
	NMin, NMax int
	Update     bool
}

// Drag is the drag of an exponential atmosphere of a body that rotates with it.
type Drag struct {
	Name             string
	Body             int
	Object           Object
	EquatorialRadius float64 // the ellipsoid of the heights
	Flattening       float64
	Layers           []Layer // 1 .. 32, by increasing Base
}

// Radiation is the pressure of the light of a source (a body with a luminosity) on the object
// as a sphere, with the shadows of the occulters (bodies with radii). C is the speed of light
// (m/s); Moving adds the terms of the velocity (Poynting-Robertson); Breaks declares the
// edges of the shadows to the integrator.
type Radiation struct {
	Name      string
	Source    int
	Object    Object
	C         float64
	Occulters []int
	Moving    bool
	Breaks    bool
}

// Schwarzschild is the relativistic correction of the field of the body Central.
type Schwarzschild struct {
	Name    string
	Central int
	PPN     PPN
}

// LenseThirring is the frame dragging of the spin of the body Central.
type LenseThirring struct {
	Name    string
	Central int
	PPN     PPN
}

// DeSitter is the geodetic precession of the motion about Central in the field of External.
type DeSitter struct {
	Name              string
	Central, External int
	PPN               PPN
}

// Empirical are accelerations in the radial, transverse and normal directions: a constant and
// a once-per-revolution part.
type Empirical struct {
	Name                   string
	Constant, Cosine, Sine [3]float64
}

// GoTerm is a term computed by a Go function, called back through cgo at every evaluation: a
// plainly slow path. Its Jacobian is taken by central differences. PositionOnly declares
// that it does not depend on the velocity.
type GoTerm struct {
	Name         string
	Func         TermFunc
	PositionOnly bool
}

// ParameterKind names the number of a term a parameter is.
type ParameterKind int

// The numbers of terms that can be estimated.
const (
	ParameterCd       ParameterKind = C.MYBGO_PARAMETER_CD       // the Cd of a Drag term
	ParameterCr       ParameterKind = C.MYBGO_PARAMETER_CR       // the Cr of a Radiation term
	ParameterConstant ParameterKind = C.MYBGO_PARAMETER_CONSTANT // Empirical.Constant[Index]
	ParameterCosine   ParameterKind = C.MYBGO_PARAMETER_COSINE   // Empirical.Cosine[Index]
	ParameterSine     ParameterKind = C.MYBGO_PARAMETER_SINE     // Empirical.Sine[Index]
)

// Parameter is a number of a term that the variational equations differentiate and a fit
// estimates: the Kind of the term named Term, with the step of its central differences.
type Parameter struct {
	Name  string
	Term  string
	Kind  ParameterKind
	Index int
	Step  float64
}

func cObject(o Object) C.myb_object {
	return C.myb_object{mass: C.double(o.Mass), area: C.double(o.Area), cr: C.double(o.Cr),
		cd: C.double(o.Cd)}
}

func cPPN(p PPN) C.myb_ppn {
	return C.myb_ppn{c: C.double(p.C), beta: C.double(p.Beta), gamma: C.double(p.Gamma)}
}

func flag(b bool) C.int {
	if b {
		return 1
	}
	return 0
}

func (p PointMass) describe(t *C.mybgo_term, k *keeper) error {
	t.kind = C.MYBGO_TERM_POINT_MASS
	t.name = k.cstring(p.Name)
	t.gm = C.double(p.GM)
	return nil
}

func (p ThirdBody) describe(t *C.mybgo_term, k *keeper) error {
	t.kind = C.MYBGO_TERM_THIRD_BODY
	t.name = k.cstring(p.Name)
	t.body = C.size_t(p.Body)
	return nil
}

func (f Field) describe(t *C.mybgo_term, k *keeper) error {
	t.kind = C.MYBGO_TERM_FIELD
	t.name = k.cstring(f.Name)
	t.body = C.size_t(f.Body)
	field, err := k.hold(f.Field.obj())
	if err != nil {
		return err
	}
	t.field = (*C.myb_gravity_field)(field)
	t.n_min = C.int(f.NMin)
	t.n_max = C.int(f.NMax)
	t.update = flag(f.Update)
	return nil
}

func (d Drag) describe(t *C.mybgo_term, k *keeper) error {
	if len(d.Layers) == 0 || len(d.Layers) > C.MYB_ATMOSPHERE_LAYERS {
		return refuse("NewModel", -1, "a drag term needs 1 .. 32 layers")
	}
	t.kind = C.MYBGO_TERM_DRAG
	t.name = k.cstring(d.Name)
	t.body = C.size_t(d.Body)
	t.object = cObject(d.Object)
	t.atmosphere.ellipsoid = C.myb_ellipsoid{a: C.double(d.EquatorialRadius), f: C.double(d.Flattening)}
	t.atmosphere.count = C.size_t(len(d.Layers))
	for i, l := range d.Layers {
		t.atmosphere.base[i] = C.double(l.Base)
		t.atmosphere.density[i] = C.double(l.Density)
		t.atmosphere.scale[i] = C.double(l.Scale)
	}
	return nil
}

func (r Radiation) describe(t *C.mybgo_term, k *keeper) error {
	if len(r.Occulters) > C.MYB_ENV_MAX_OCCULTERS {
		return refuse("NewModel", -1, "more than 8 occulters")
	}
	t.kind = C.MYBGO_TERM_RADIATION
	t.name = k.cstring(r.Name)
	t.body = C.size_t(r.Source)
	t.object = cObject(r.Object)
	t.c = C.double(r.C)
	for i, o := range r.Occulters {
		t.occulters[i] = C.size_t(o)
	}
	t.occulter_count = C.size_t(len(r.Occulters))
	t.moving = flag(r.Moving)
	t.breaks = flag(r.Breaks)
	return nil
}

func (s Schwarzschild) describe(t *C.mybgo_term, k *keeper) error {
	t.kind = C.MYBGO_TERM_SCHWARZSCHILD
	t.name = k.cstring(s.Name)
	t.body = C.size_t(s.Central)
	t.ppn = cPPN(s.PPN)
	return nil
}

func (l LenseThirring) describe(t *C.mybgo_term, k *keeper) error {
	t.kind = C.MYBGO_TERM_LENSE_THIRRING
	t.name = k.cstring(l.Name)
	t.body = C.size_t(l.Central)
	t.ppn = cPPN(l.PPN)
	return nil
}

func (d DeSitter) describe(t *C.mybgo_term, k *keeper) error {
	t.kind = C.MYBGO_TERM_DE_SITTER
	t.name = k.cstring(d.Name)
	t.body = C.size_t(d.Central)
	t.external = C.size_t(d.External)
	t.ppn = cPPN(d.PPN)
	return nil
}

func (g GoTerm) describe(t *C.mybgo_term, k *keeper) error {
	if g.Func == nil {
		return refuse("NewModel", -1, "a Go term without a function")
	}
	h := cgo.NewHandle(&termCallback{f: g.Func, state: k.state})
	k.handles = append(k.handles, h)
	t.kind = C.MYBGO_TERM_GO
	t.name = k.cstring(g.Name)
	t.handle = C.uintptr_t(h)
	t.position_only = flag(g.PositionOnly)
	return nil
}

func (e Empirical) describe(t *C.mybgo_term, k *keeper) error {
	t.kind = C.MYBGO_TERM_EMPIRICAL
	t.name = k.cstring(e.Name)
	for i := range 3 {
		t.constant[i] = C.double(e.Constant[i])
		t.cosine[i] = C.double(e.Cosine[i])
		t.sine[i] = C.double(e.Sine[i])
	}
	return nil
}

// ---- The model ----------------------------------------------------------------------------------

// Model is a force model of the library with its environment, made from a ModelSpec. It
// serves one call at a time (calls from several goroutines wait for each other).
type Model struct {
	o          *object
	mu         sync.Mutex
	state      *callState
	terms      []string
	parameters []string
}

// NewModel makes the environment and the model of a description (one call into C).
func NewModel(spec ModelSpec) (*Model, error) {
	const op = "NewModel"
	k := newKeeper()
	ok := false
	defer func() {
		if !ok {
			k.done()
		}
	}()
	if !validTime(spec.Epoch) {
		return nil, refuse(op, -1, "an epoch out of range")
	}
	if spec.BreakScheme != BreaksRetake && spec.BreakScheme != BreaksPredict {
		return nil, refuse(op, -1, "an unknown break scheme")
	}
	options := C.mybgo_model_options{epoch: cTime(spec.Epoch), scale: k.optional(spec.Scale),
		table_step: C.double(spec.TableStep), break_scheme: C.int(spec.BreakScheme)}
	sys, err := k.hold(spec.TimeSystem.obj())
	if err != nil {
		return nil, err
	}
	options.timesys = (*C.myb_timesys)(sys)

	bodies := make([]C.mybgo_body, len(spec.Bodies))
	for i, b := range spec.Bodies {
		d := &bodies[i]
		d.id = C.int64_t(b.ID)
		d.gm = C.double(b.GM)
		d.radius = C.double(b.Radius)
		d.luminosity = C.double(b.Luminosity)
		d.limb = [2]C.double{C.double(b.Limb[0]), C.double(b.Limb[1])}
		d.spin = [3]C.double{C.double(b.Spin[0]), C.double(b.Spin[1]), C.double(b.Spin[2])}
		if err := describeProvider(b.State, &d.state, k); err != nil {
			return nil, err
		}
		if b.State != nil {
			d.state_scale = k.optional(b.StateScale)
		}
		if b.Orientation != nil {
			if err := b.Orientation.describe(&d.orientation, k); err != nil {
				return nil, err
			}
			d.orientation_scale = k.optional(b.OrientationScale)
		}
	}

	terms := make([]C.mybgo_term, len(spec.Terms))
	names := make([]string, len(spec.Terms))
	for j, term := range spec.Terms {
		if term == nil {
			return nil, refuse(op, j, "a nil term")
		}
		if err := term.describe(&terms[j], k); err != nil {
			if e, ok := err.(*Error); ok && e.Index < 0 {
				e.Index = j
			}
			return nil, err
		}
		names[j] = C.GoString(terms[j].name)
	}

	parameters := make([]C.mybgo_parameter, len(spec.Parameters))
	parameterNames := make([]string, len(spec.Parameters))
	for p, par := range spec.Parameters {
		index := slices.Index(names, par.Term)
		if index < 0 {
			return nil, &Error{Status: NotFound, Op: op, Index: p, Message: "no term " + par.Term}
		}
		parameters[p] = C.mybgo_parameter{name: k.cstring(par.Name), term: C.size_t(index),
			what: C.int(par.Kind), index: C.int(par.Index), step: C.double(par.Step)}
		parameterNames[p] = par.Name
	}

	var model *C.mybgo_model
	var stage C.int
	var item C.size_t
	status := C.mybgo_model_create(&options, (*C.mybgo_body)(first(bodies)), C.size_t(len(bodies)),
		(*C.mybgo_term)(first(terms)), C.size_t(len(terms)),
		(*C.mybgo_parameter)(first(parameters)), C.size_t(len(parameters)), &model, &stage, &item)
	if status != C.MYB_OK {
		return nil, stageError(op, status, stage, item, "body", "term", "parameter")
	}
	held, handles := k.keep()
	ok = true
	m := &Model{state: k.state, terms: names, parameters: parameterNames}
	m.o = newObject(unsafe.Pointer(model), func(q unsafe.Pointer) {
		C.mybgo_model_destroy((*C.mybgo_model)(q))
		for _, h := range handles {
			h.Delete()
		}
	}, held...)
	return m, nil
}

// stageError is the error of a description that C refused: the record of the stage named.
func stageError(op string, status C.myb_status, stage C.int, item C.size_t, names ...string) error {
	e := &Error{Status: Status(status), Op: op, Index: -1, Message: "environment"}
	switch stage {
	case C.MYBGO_STAGE_BODY:
		e.Index, e.Message = int(item), names[0]
	case C.MYBGO_STAGE_TERM:
		e.Index, e.Message = int(item), names[1]
	case C.MYBGO_STAGE_PARAMETER:
		e.Index, e.Message = int(item), names[2]
	}
	return e
}

// Close releases the model, and with it its hold on the objects it uses.
func (m *Model) Close() error {
	if m != nil {
		m.o.close()
	}
	return nil
}

// begin holds the model for one call: the lock, the object, a clean record of the errors of
// the Go functions.
func (m *Model) begin() (*C.mybgo_model, error) {
	if m == nil {
		return nil, ErrClosed
	}
	m.mu.Lock()
	ptr, err := m.o.acquire()
	if err != nil {
		m.mu.Unlock()
		return nil, err
	}
	m.state.err = nil
	return (*C.mybgo_model)(ptr), nil
}

func (m *Model) end() {
	m.o.release()
	m.mu.Unlock()
}

// failure is the error of a call of the model: with the error of a Go function that caused it.
func (m *Model) failure(op string, status C.myb_status, index int) error {
	if status == C.MYB_OK {
		return nil
	}
	return &Error{Status: Status(status), Op: op, Index: index, Err: m.state.err}
}

// Terms are the names of the terms, in the order of the description.
func (m *Model) Terms() []string {
	return slices.Clone(m.terms)
}

// ParameterNames are the names of the parameters, in the order of the description.
func (m *Model) ParameterNames() []string {
	return slices.Clone(m.parameters)
}

// Parameters are the values of the parameters now (a fit changes them).
func (m *Model) Parameters() ([]float64, error) {
	model, err := m.begin()
	if err != nil {
		return nil, err
	}
	defer m.end()
	values := make([]float64, int(C.mybgo_model_parameter_count(model)))
	C.mybgo_model_get_parameters(model, cDoubles(values))
	return values, nil
}

// SetParameters gives the parameters new values (MYB_E_DOMAIN if one is not finite).
func (m *Model) SetParameters(values []float64) error {
	model, err := m.begin()
	if err != nil {
		return err
	}
	defer m.end()
	if len(values) != int(C.mybgo_model_parameter_count(model)) {
		return refuse("Model.SetParameters", -1, "not one value for every parameter")
	}
	return fail("Model.SetParameters", C.mybgo_model_set_parameters(model, cDoubles(values)))
}

// Acceleration is the sum of the terms at the model time t and the state s.
func (m *Model) Acceleration(t float64, s State) ([3]float64, error) {
	a, err := m.Accelerations([]float64{t}, []State{s})
	if err != nil {
		if e, ok := err.(*Error); ok {
			e.Op, e.Index = "Model.Acceleration", -1
		}
		return [3]float64{}, err
	}
	return a[0], nil
}

// Accelerations is the sum of the terms at every point (t[k], s[k]), in one call.
func (m *Model) Accelerations(t []float64, s []State) ([][3]float64, error) {
	const op = "Model.Accelerations"
	if len(t) != len(s) {
		return nil, refuse(op, -1, "times and states of different lengths")
	}
	model, err := m.begin()
	if err != nil {
		return nil, err
	}
	defer m.end()
	out := make([][3]float64, len(t))
	var failed C.size_t
	status := C.mybgo_model_accelerations(model, cDoubles(t), cStates(s), C.size_t(len(t)),
		(*C.double)(first(out)), &failed)
	return out, m.failure(op, status, int(failed))
}

// ---- Propagation ---------------------------------------------------------------------------------

// Integrator selects the integrator of a propagation (myb_orbit_integrator).
type Integrator int

// The integrators of the library.
const (
	IAS15  Integrator = C.MYB_ORBIT_IAS15
	DOP853 Integrator = C.MYB_ORBIT_DOP853
	RKN12  Integrator = C.MYB_ORBIT_RKN12
)

// OrbitOptions are the options of a propagation (myb_orbit_options); zero values are the
// defaults of the library.
type OrbitOptions struct {
	Epsilon     float64 // IAS15 precision parameter; 0 = 1e-10
	HMax        float64 // largest step (s); 0 = unlimited
	MaxSteps    int64   // budget of steps; 0 = 10^9
	Variational bool    // integrate the variational equations (set by Propagation.Partials)
	Integrator  Integrator
	RTol, ATol  float64 // DOP853 and RKN12 tolerances; 0 = 1e-12, 1e-9
	Timescale   IAS15Timescale
}

// IAS15Timescale is the timescale that sets the step of IAS15 (myb_ias15_timescale).
type IAS15Timescale int

// The timescales.
const (
	TimescaleAcceleration IAS15Timescale = C.MYB_IAS15_TIMESCALE_ACCELERATION // the default
	TimescaleJerk         IAS15Timescale = C.MYB_IAS15_TIMESCALE_JERK         // for a contact
)

func (o OrbitOptions) c() C.myb_orbit_options {
	return C.myb_orbit_options{epsilon: C.double(o.Epsilon), h_max: C.double(o.HMax),
		max_steps: C.int64_t(o.MaxSteps), variational: flag(o.Variational),
		integrator: C.myb_orbit_integrator(o.Integrator), rtol: C.double(o.RTol),
		atol: C.double(o.ATol), timescale: C.myb_ias15_timescale(o.Timescale)}
}

// Impulse is an impulsive manoeuvre: at the model time T the velocity changes by DV (m/s).
type Impulse struct {
	T  float64
	DV [3]float64
}

// Propagation asks for the states of an object from Y0 at T0 towards TEnd (never passed) at
// Times, in the order of the propagation; with Partials also the transition matrices and the
// sensitivities to the parameters; with the Impulses applied on the way (in the same order; a
// state asked at the time of an impulse is the one after it).
type Propagation struct {
	Options  OrbitOptions
	T0       float64
	Y0       State
	TEnd     float64
	Times    []float64
	Partials bool
	Impulses []Impulse
}

// Trajectory is what a propagation gives.
type Trajectory struct {
	States      []State
	Phi         [][36]float64 // with Partials: d y / d y0, row-major
	Sensitivity [][]float64   // with Partials: d y / d p, 6 rows of len(parameters)
	Steps       int64         // the accepted steps of the integrator
}

// Propagate integrates the object through the model (one call into C).
func (m *Model) Propagate(p Propagation) (*Trajectory, error) {
	const op = "Model.Propagate"
	model, err := m.begin()
	if err != nil {
		return nil, err
	}
	defer m.end()
	options := p.Options
	if p.Partials {
		options.Variational = true
	}
	co := options.c()
	kicks := make([]float64, 0, 4*len(p.Impulses))
	for _, imp := range p.Impulses {
		kicks = append(kicks, imp.T, imp.DV[0], imp.DV[1], imp.DV[2])
	}
	out := &Trajectory{States: make([]State, len(p.Times))}
	var sens []float64
	np := int(C.mybgo_model_parameter_count(model))
	if p.Partials {
		out.Phi = make([][36]float64, len(p.Times))
		sens = make([]float64, 6*np*len(p.Times))
	}
	var steps C.int64_t
	var failed C.size_t
	y0 := p.Y0
	status := C.mybgo_propagate(model, &co, C.double(p.T0), (*C.double)(unsafe.Pointer(&y0)), C.double(p.TEnd),
		cDoubles(p.Times), C.size_t(len(p.Times)), cStates(out.States), (*C.double)(first(out.Phi)),
		cDoubles(sens), cDoubles(kicks), C.size_t(len(p.Impulses)), &steps, &failed)
	out.Steps = int64(steps)
	if err := m.failure(op, status, int(failed)); err != nil {
		return nil, err
	}
	if p.Partials {
		out.Sensitivity = make([][]float64, len(p.Times))
		for k := range out.Sensitivity {
			out.Sensitivity[k] = sens[6*np*k : 6*np*(k+1) : 6*np*(k+1)]
		}
	}
	return out, nil
}

// ---- Events -------------------------------------------------------------------------------------

// Direction selects the roots of an event function: where it grows, where it falls, or both.
type Direction int

// The directions of a root.
const (
	Any     Direction = 0
	Rising  Direction = 1
	Falling Direction = -1
)

// Origin stands for the origin of the environment in place of a body (MYB_EVENT_ORIGIN).
const Origin = -1

// The most roots of an event a search reports when its Limit is 0.
const defaultLimit = 1 << 20

// EventSpec is an event to look for along an orbit.
type EventSpec interface {
	describe(e *C.mybgo_event, k *keeper) error
}

// PlaneCrossing is the crossing of the plane through a body (or the Origin) with the normal
// given: Rising at the ascending node (myb_event_plane).
type PlaneCrossing struct {
	Body      int
	Normal    [3]float64
	Direction Direction
	Limit     int
}

// Distance is the crossing of a sphere of radius Range about a body: Falling on the way in
// (myb_event_range).
type Distance struct {
	Body      int
	Range     float64
	Direction Direction
	Limit     int
}

// Apsis is a zero of the rate of the distance to a body: Rising at the pericentre, Falling at
// the apocentre (myb_event_range_rate).
type Apsis struct {
	Body      int
	Direction Direction
	Limit     int
}

// Altitude is the crossing of a geodetic height over the ellipsoid of a body with an
// orientation: Falling on the way down (myb_event_altitude).
type Altitude struct {
	Body      int
	Ellipsoid Ellipsoid
	Altitude  float64
	Direction Direction
	Limit     int
}

// Elevation is the crossing of an elevation over the horizon of a station on a body: Rising
// at the rise (myb_event_elevation).
type Elevation struct {
	Station   Station
	Elevation float64
	Direction Direction
	Limit     int
}

// Eclipse is an edge of the shadow of Occulter on the light of Source: Falling on the way in;
// Inner for the umbra (myb_event_eclipse).
type Eclipse struct {
	Source, Occulter int
	Inner            bool
	Direction        Direction
	Limit            int
}

// GoEvent is an event function in Go, called back through cgo at every evaluation: a plainly
// slow path.
type GoEvent struct {
	Func      EventFunc
	Direction Direction
	Limit     int
}

func bodyOrOrigin(body int) C.size_t {
	if body == Origin {
		return ^C.size_t(0)
	}
	return C.size_t(body)
}

func limit(n int) C.size_t {
	if n <= 0 {
		return defaultLimit
	}
	return C.size_t(n)
}

func (p PlaneCrossing) describe(e *C.mybgo_event, k *keeper) error {
	e.kind = C.MYBGO_EVENT_PLANE
	e.direction = C.int(p.Direction)
	e.limit = limit(p.Limit)
	e.body = bodyOrOrigin(p.Body)
	e.normal = [3]C.double{C.double(p.Normal[0]), C.double(p.Normal[1]), C.double(p.Normal[2])}
	return nil
}

func (d Distance) describe(e *C.mybgo_event, k *keeper) error {
	e.kind = C.MYBGO_EVENT_RANGE
	e.direction = C.int(d.Direction)
	e.limit = limit(d.Limit)
	e.body = bodyOrOrigin(d.Body)
	e._range = C.double(d.Range)
	return nil
}

func (a Apsis) describe(e *C.mybgo_event, k *keeper) error {
	e.kind = C.MYBGO_EVENT_RANGE_RATE
	e.direction = C.int(a.Direction)
	e.limit = limit(a.Limit)
	e.body = bodyOrOrigin(a.Body)
	return nil
}

func (a Altitude) describe(e *C.mybgo_event, k *keeper) error {
	e.kind = C.MYBGO_EVENT_ALTITUDE
	e.direction = C.int(a.Direction)
	e.limit = limit(a.Limit)
	e.body = C.size_t(a.Body)
	e.ellipsoid = C.myb_ellipsoid{a: C.double(a.Ellipsoid.A), f: C.double(a.Ellipsoid.F)}
	e.altitude = C.double(a.Altitude)
	return nil
}

func (s Elevation) describe(e *C.mybgo_event, k *keeper) error {
	e.kind = C.MYBGO_EVENT_ELEVATION
	e.direction = C.int(s.Direction)
	e.limit = limit(s.Limit)
	e.body = C.size_t(s.Station.Body)
	e.ellipsoid = C.myb_ellipsoid{a: C.double(s.Station.Ellipsoid.A), f: C.double(s.Station.Ellipsoid.F)}
	e.longitude = C.double(s.Station.Longitude)
	e.latitude = C.double(s.Station.Latitude)
	e.height = C.double(s.Station.Height)
	e.elevation = C.double(s.Elevation)
	return nil
}

func (s Eclipse) describe(e *C.mybgo_event, k *keeper) error {
	e.kind = C.MYBGO_EVENT_ECLIPSE
	e.direction = C.int(s.Direction)
	e.limit = limit(s.Limit)
	e.source = C.size_t(s.Source)
	e.occulter = C.size_t(s.Occulter)
	e.inner = flag(s.Inner)
	return nil
}

func (g GoEvent) describe(e *C.mybgo_event, k *keeper) error {
	if g.Func == nil {
		return refuse("Model.Events", -1, "a Go event without a function")
	}
	h := cgo.NewHandle(&eventCallback{f: g.Func})
	k.handles = append(k.handles, h)
	e.kind = C.MYBGO_EVENT_GO
	e.direction = C.int(g.Direction)
	e.limit = limit(g.Limit)
	e.handle = C.uintptr_t(h)
	return nil
}

// EventSearch asks for the roots of the Events along the orbit from Y0 at T0 to TEnd. Each
// event is searched on an orbit of its own — the same orbit, since the steps do not depend
// on what is asked —, so the work grows with the number of events.
type EventSearch struct {
	Options OrbitOptions
	T0      float64
	Y0      State
	TEnd    float64
	Events  []EventSpec
}

// Event is a root found: the index of its EventSpec, the time and the state there.
type Event struct {
	Index int
	T     float64
	State State
}

// EventResult are the roots in the order of the propagation (those at one time in the order
// of the events), the state at TEnd and the accepted steps of the propagation.
type EventResult struct {
	Events []Event
	Final  State
	Steps  int64
}

// Events finds the events of an orbit (one call into C).
func (m *Model) Events(s EventSearch) (*EventResult, error) {
	const op = "Model.Events"
	k := newKeeper()
	defer k.done()
	events := make([]C.mybgo_event, len(s.Events))
	for i, e := range s.Events {
		if e == nil {
			return nil, refuse(op, i, "a nil event")
		}
		if err := e.describe(&events[i], k); err != nil {
			return nil, err
		}
	}
	model, err := m.begin()
	if err != nil {
		return nil, err
	}
	defer m.end()
	co := s.Options.c()
	var found *C.mybgo_found
	var count, failed C.size_t
	var steps C.int64_t
	y0 := s.Y0
	var final State
	status := C.mybgo_events(model, &co, C.double(s.T0), (*C.double)(unsafe.Pointer(&y0)),
		C.double(s.TEnd), (*C.mybgo_event)(first(events)), C.size_t(len(events)), &found, &count,
		(*C.double)(unsafe.Pointer(&final)), &steps, &failed)
	if status != C.MYB_OK {
		return nil, m.failure(op, status, int(failed))
	}
	defer C.mybgo_free(unsafe.Pointer(found))
	out := &EventResult{Events: make([]Event, int(count)), Final: final, Steps: int64(steps)}
	for i, f := range unsafe.Slice(found, int(count)) {
		out.Events[i] = Event{Index: int(f.event), T: float64(f.t),
			State: *(*State)(unsafe.Pointer(&f.y[0]))}
	}
	sortEvents(out.Events, s.TEnd > s.T0)
	return out, nil
}

// ---- Estimation ---------------------------------------------------------------------------------

// Station is an observer at geodetic coordinates (rad, m) on the ellipsoid of a body.
type Station struct {
	Body                        int
	Ellipsoid                   Ellipsoid
	Longitude, Latitude, Height float64
}

// MeasureKind is what an observation measures.
type MeasureKind int

// The measurements of the library: of the object itself, and from a station.
const (
	MeasurePosition  MeasureKind = C.MYBGO_MEASURE_POSITION   // 3 values, m
	MeasureState     MeasureKind = C.MYBGO_MEASURE_STATE      // 6 values, m and m/s
	MeasureRange     MeasureKind = C.MYBGO_MEASURE_RANGE      // 1 value, m
	MeasureRangeRate MeasureKind = C.MYBGO_MEASURE_RANGE_RATE // 1 value, m/s
	MeasureRaDec     MeasureKind = C.MYBGO_MEASURE_RADEC      // 2 values, rad
	MeasureAzEl      MeasureKind = C.MYBGO_MEASURE_AZEL       // 2 values, rad
)

func (k MeasureKind) values() int {
	switch k {
	case MeasurePosition:
		return 3
	case MeasureState:
		return 6
	case MeasureRange, MeasureRangeRate:
		return 1
	case MeasureRaDec, MeasureAzEl:
		return 2
	}
	return 0
}

// Observation is a measurement at the model time T: Value and Sigma hold as many numbers as
// the Kind gives; Station indexes FitRequest.Stations; LightTime takes the light time into
// account (FitRequest.Light); Bias > 0 adds the systematic shifts Bias - 1 + i.
type Observation struct {
	T         float64
	Kind      MeasureKind
	Station   int
	LightTime bool
	Value     []float64
	Sigma     []float64
	Bias      int
}

// LightTime describes the measurements with the light time: the speed of light C (m/s), the GM
// of the body at the origin (0: straight lines), two-way, aberration.
type LightTime struct {
	C, GM              float64
	TwoWay, Aberration bool
}

// FitRequest is a fit of the state at T0 (first guess Y0) and of the parameters of the model
// to the observations, ordered in time away from T0 (myb_fit_orbit). Fixed and PriorSigma
// (6 + parameters + biases each, or nil) keep unknowns and give a priori information; Bias
// is the first guess of the systematic shifts.
type FitRequest struct {
	T0            float64
	Y0            State
	Observations  []Observation
	Stations      []Station
	Light         *LightTime
	MaxIterations int     // 0 = 10
	Tolerance     float64 // 0 = 1e-6
	Epsilon       float64 // 0 = 1e-10
	Fixed         []bool
	PriorSigma    []float64
	Bias          []float64
	Covariance    bool // ask for the covariance of the unknowns
	Residuals     bool // ask for the residuals at the estimate
}

// FitResult is the estimate.
type FitResult struct {
	Y0         State
	Parameters []float64
	Bias       []float64
	Iterations int
	Converged  bool
	Values     int     // scalar measurements
	RMSBefore  float64 // weighted RMS at the first guess
	RMS        float64 // weighted RMS at the estimate
	Covariance []float64
	Residuals  []float64
}

// Fit estimates the orbit and the parameters (one call into C); the parameters of the model
// keep the estimate. On an error the model keeps the values of the last iteration.
func (m *Model) Fit(f FitRequest) (*FitResult, error) {
	const op = "Model.Fit"
	observations := make([]C.mybgo_observation, len(f.Observations))
	values := 0
	for i, o := range f.Observations {
		n := o.Kind.values()
		if n == 0 || len(o.Value) != n || len(o.Sigma) != n {
			return nil, refuse(op, i, "an observation without the values its kind needs")
		}
		d := &observations[i]
		d.t = C.double(o.T)
		d.kind = C.int(o.Kind)
		if o.LightTime {
			d.kind |= C.MYBGO_MEASURE_LIGHT
		}
		d.station = C.size_t(o.Station)
		d.bias = C.int(o.Bias)
		for j := range n {
			d.value[j] = C.double(o.Value[j])
			d.sigma[j] = C.double(o.Sigma[j])
		}
		values += n
	}
	stations := make([]C.mybgo_station, len(f.Stations))
	for i, s := range f.Stations {
		stations[i] = C.mybgo_station{body: C.size_t(s.Body),
			ellipsoid: C.myb_ellipsoid{a: C.double(s.Ellipsoid.A), f: C.double(s.Ellipsoid.F)},
			longitude: C.double(s.Longitude), latitude: C.double(s.Latitude),
			height: C.double(s.Height)}
	}
	var light *C.mybgo_light
	if f.Light != nil {
		light = &C.mybgo_light{c: C.double(f.Light.C), gm: C.double(f.Light.GM),
			two_way: flag(f.Light.TwoWay), aberration: flag(f.Light.Aberration)}
	}
	model, err := m.begin()
	if err != nil {
		return nil, err
	}
	defer m.end()
	np := int(C.mybgo_model_parameter_count(model))
	unknowns := 6 + np + len(f.Bias)
	if (f.Fixed != nil && len(f.Fixed) != unknowns) ||
		(f.PriorSigma != nil && len(f.PriorSigma) != unknowns) {
		return nil, refuse(op, -1, "Fixed and PriorSigma need 6 + parameters + biases entries")
	}
	var fixed []C.uchar
	if f.Fixed != nil {
		fixed = make([]C.uchar, unknowns)
		for i, b := range f.Fixed {
			if b {
				fixed[i] = 1
			}
		}
	}
	options := C.mybgo_fit_options{max_iterations: C.int(f.MaxIterations),
		tolerance: C.double(f.Tolerance), epsilon: C.double(f.Epsilon),
		bias_count: C.size_t(len(f.Bias))}
	out := &FitResult{Y0: f.Y0, Bias: slices.Clone(f.Bias)}
	if f.Covariance {
		out.Covariance = make([]float64, unknowns*unknowns)
	}
	if f.Residuals {
		out.Residuals = make([]float64, values)
	}
	var result C.myb_fit_result
	var failed C.size_t
	y0 := f.Y0
	status := C.mybgo_fit(model, C.double(f.T0), (*C.double)(unsafe.Pointer(&y0)),
		(*C.mybgo_station)(first(stations)), C.size_t(len(stations)), light,
		(*C.mybgo_observation)(first(observations)), C.size_t(len(observations)), &options,
		(*C.uchar)(first(fixed)), cDoubles(f.PriorSigma), cDoubles(out.Bias), &result,
		cDoubles(out.Covariance), cDoubles(out.Residuals), &failed)
	if status != C.MYB_OK {
		index := -1
		if int(failed) < len(observations) {
			index = int(failed)
		}
		return nil, m.failure(op, status, index)
	}
	out.Y0 = y0
	out.Parameters = make([]float64, np)
	C.mybgo_model_get_parameters(model, cDoubles(out.Parameters))
	out.Iterations = int(result.iterations)
	out.Converged = result.converged != 0
	out.Values = int(result.values)
	out.RMSBefore = float64(result.rms_before)
	out.RMS = float64(result.rms)
	return out, nil
}

// sortEvents orders the roots by time in the direction of the propagation, those at one time
// by the order of the events: comparisons only.
func sortEvents(events []Event, forward bool) {
	slices.SortStableFunc(events, func(a, b Event) int {
		if forward {
			return cmp.Compare(a.T, b.T)
		}
		return cmp.Compare(b.T, a.T)
	})
}
