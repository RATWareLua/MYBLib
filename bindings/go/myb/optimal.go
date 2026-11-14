// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include "glue.h"
*/
import "C"

import "unsafe"

// fieldOf holds the force model of an optional *Model for a call (nil: none).
func fieldOf(m *Model) (*C.myb_force_model, func(), error) {
	if m == nil {
		return nil, func() {}, nil
	}
	model, err := m.begin()
	if err != nil {
		return nil, nil, err
	}
	return C.mybgo_model_force(model), m.end, nil
}

func fieldFailure(op string, status C.myb_status, m *Model) error {
	if status == C.MYB_OK {
		return nil
	}
	e := &Error{Status: Status(status), Op: op, Index: -1}
	if m != nil {
		e.Err = m.state.err
	}
	return e
}

// ---- The Q-law ----------------------------------------------------------------------------------

// QLaw is the feedback law of Petropoulos (myb_qlaw): the target a, e, i, raan, argp, their
// weights (0 leaves one free), the lowest periapsis and its weight.
type QLaw struct {
	MU       float64
	Target   [5]float64
	Weight   [5]float64
	RPMin    float64
	RPWeight float64
}

func (q QLaw) c() C.myb_qlaw {
	out := C.myb_qlaw{mu: C.double(q.MU), rp_min: C.double(q.RPMin), rp_weight: C.double(q.RPWeight)}
	for i := range 5 {
		out.target[i] = C.double(q.Target[i])
		out.weight[i] = C.double(q.Weight[i])
	}
	return out
}

// QLawDirections are the directions of the thrust (radial, transverse, normal), Q, its rate
// and the effectivity at every set of elements (myb_qlaw_direction), in one call.
func QLawDirections(law QLaw, elements []Elements) (direction [][3]float64, q, rate, effectivity []float64, err error) {
	cl := law.c()
	direction = make([][3]float64, len(elements))
	q = make([]float64, len(elements))
	rate = make([]float64, len(elements))
	effectivity = make([]float64, len(elements))
	var failed C.size_t
	status := C.mybgo_qlaw_directions(&cl, (*C.myb_elements)(first(elements)), C.size_t(len(elements)),
		(*C.double)(first(direction)), cDoubles(q), cDoubles(rate), cDoubles(effectivity), &failed)
	if err := failAt("QLawDirections", status, failed); err != nil {
		return nil, nil, nil, nil, err
	}
	return direction, q, rate, effectivity, nil
}

// QLawEngine is the engine that flies the law (myb_qlaw_term): thrust (N), exhaust velocity
// (m/s), the tolerances of the elements and the threshold of the effectivity.
type QLawEngine struct {
	Thrust, Exhaust float64
	Tolerance       [5]float64
	Effectivity     float64
}

// QLawArrival is where a flight ended (myb_qlaw_arrival): the time, x, v and the mass.
type QLawArrival struct {
	T       float64
	Y       [7]float64
	Arrived bool
}

// QLawFly flies the law by DOP853 from y0 (x, v, mass) at t0 until the orbit is within the
// tolerances, looked at every step, or until tLimit (myb_qlaw_fly). field: everything but the
// engine (nil: the point mass of the law).
func QLawFly(law QLaw, engine QLawEngine, field *Model, rtol, atol, t0 float64, y0 [7]float64, tLimit, step float64) (QLawArrival, error) {
	const op = "QLawFly"
	force, release, err := fieldOf(field)
	if err != nil {
		return QLawArrival{}, err
	}
	defer release()
	cl := law.c()
	term := C.myb_qlaw_term{thrust: C.double(engine.Thrust), exhaust: C.double(engine.Exhaust),
		effectivity: C.double(engine.Effectivity)}
	for i := range 5 {
		term.tolerance[i] = C.double(engine.Tolerance[i])
	}
	var a C.myb_qlaw_arrival
	status := C.mybgo_qlaw_fly(&cl, &term, force, C.double(rtol), C.double(atol), C.double(t0),
		(*C.double)(unsafe.Pointer(&y0[0])), C.double(tLimit), C.double(step), &a)
	out := QLawArrival{T: float64(a.t), Arrived: a.arrived != 0}
	for i := range 7 {
		out.Y[i] = float64(a.y[i])
	}
	return out, fieldFailure(op, status, field)
}

// ---- The indirect method ---------------------------------------------------------------------

// IndirectGoal is what the indirect method minimises (myb_indirect_goal).
type IndirectGoal int

// The goals.
const (
	LeastFuel IndirectGoal = C.MYB_INDIRECT_FUEL
	LeastTime IndirectGoal = C.MYB_INDIRECT_TIME
)

// IndirectArrival is where a transfer of the indirect method ends (myb_indirect_arrival).
type IndirectArrival int

// The arrivals.
const (
	ArriveRendezvous IndirectArrival = C.MYB_INDIRECT_RENDEZVOUS
	ArriveOnOrbit    IndirectArrival = C.MYB_INDIRECT_ORBIT // the phase free
)

// IndirectProblem is a transfer for the indirect method (myb_indirect_problem).
type IndirectProblem struct {
	MU        float64
	Field     *Model // everything but the engine; nil: the point mass MU
	Thrust    float64
	Exhaust   float64
	T0        float64
	Start     State
	Mass0     float64
	T1        float64
	Target    State
	Goal      IndirectGoal
	Smoothing float64
	Tolerance float64
	Tries     int
	Arrival   IndirectArrival
	Phase     float64
}

func (p IndirectProblem) c(field *C.myb_force_model) C.myb_indirect_problem {
	out := C.myb_indirect_problem{mu: C.double(p.MU), field: field, thrust: C.double(p.Thrust),
		exhaust: C.double(p.Exhaust), t0: C.double(p.T0), mass0: C.double(p.Mass0), t1: C.double(p.T1),
		goal: C.myb_indirect_goal(p.Goal), smoothing: C.double(p.Smoothing),
		tolerance: C.double(p.Tolerance), tries: C.int(p.Tries),
		arrival: C.myb_indirect_arrival(p.Arrival), phase: C.double(p.Phase)}
	for i := range 6 {
		out.start[i] = C.double(p.Start[i])
		out.target[i] = C.double(p.Target[i])
	}
	return out
}

// IndirectSolution is a solution of the indirect method (myb_indirect_solution).
type IndirectSolution struct {
	Costate   [8]float64
	T1        float64
	Mass1     float64
	Smoothing float64
	Residual  float64
	Tries     int
	Stages    int
	Phase     float64
}

func goIndirect(c *C.myb_indirect_solution) IndirectSolution {
	out := IndirectSolution{T1: float64(c.t1), Mass1: float64(c.mass1), Smoothing: float64(c.smoothing),
		Residual: float64(c.residual), Tries: int(c.tries), Stages: int(c.stages), Phase: float64(c.phase)}
	for i := range 8 {
		out.Costate[i] = float64(c.costate[i])
	}
	return out
}

func (s IndirectSolution) c() C.myb_indirect_solution {
	out := C.myb_indirect_solution{t1: C.double(s.T1), mass1: C.double(s.Mass1),
		smoothing: C.double(s.Smoothing), residual: C.double(s.Residual), tries: C.int(s.Tries),
		stages: C.int(s.Stages), phase: C.double(s.Phase)}
	for i := range 8 {
		out.costate[i] = C.double(s.Costate[i])
	}
	return out
}

// IndirectSolve solves the problem (myb_indirect_solve); guess: the costates of a solution
// nearby, or nil. On MYB_E_NOCONV the last solution found comes back with the error.
func IndirectSolve(p IndirectProblem, guess *[8]float64) (IndirectSolution, error) {
	const op = "IndirectSolve"
	force, release, err := fieldOf(p.Field)
	if err != nil {
		return IndirectSolution{}, err
	}
	defer release()
	cp := p.c(force)
	var g *C.double
	if guess != nil {
		copied := *guess
		g = (*C.double)(unsafe.Pointer(&copied[0]))
	}
	var s C.myb_indirect_solution
	status := C.myb_indirect_solve(&cp, g, nil, &s)
	return goIndirect(&s), fieldFailure(op, status, p.Field)
}

// IndirectPoint is a sample of a solution (myb_indirect_point).
type IndirectPoint struct {
	Y           [7]float64
	Direction   [3]float64
	Throttle    float64
	Switching   float64
	Hamiltonian float64
	Costate     [7]float64
}

// IndirectSample is the solution at the times (myb_indirect_sample), in one call.
func IndirectSample(p IndirectProblem, s IndirectSolution, times []float64) ([]IndirectPoint, error) {
	const op = "IndirectSample"
	force, release, err := fieldOf(p.Field)
	if err != nil {
		return nil, err
	}
	defer release()
	cp, cs := p.c(force), s.c()
	c := make([]C.myb_indirect_point, len(times))
	status := C.myb_indirect_sample(&cp, &cs, cDoubles(times), C.size_t(len(times)), nil,
		(*C.myb_indirect_point)(first(c)))
	if err := fieldFailure(op, status, p.Field); err != nil {
		return nil, err
	}
	out := make([]IndirectPoint, len(c))
	for k, v := range c {
		out[k] = IndirectPoint{Direction: goVec3(v.direction), Throttle: float64(v.throttle),
			Switching: float64(v.switching), Hamiltonian: float64(v.hamiltonian)}
		for i := range 7 {
			out[k].Y[i] = float64(v.y[i])
			out[k].Costate[i] = float64(v.costate[i])
		}
	}
	return out, nil
}

// ---- Sims and Flanagan ----------------------------------------------------------------------------

// NLPOptions are the options of the solver of nonlinear programs (myb_nlp_options).
type NLPOptions struct {
	Feasibility, Optimality float64
	MaxOuter, MaxInner      int
	Weight                  float64
}

func (o *NLPOptions) c() *C.myb_nlp_options {
	if o == nil {
		return nil
	}
	return &C.myb_nlp_options{feasibility: C.double(o.Feasibility), optimality: C.double(o.Optimality),
		max_outer: C.int(o.MaxOuter), max_inner: C.int(o.MaxInner), weight: C.double(o.Weight)}
}

// NLPResult is what the solver reached (myb_nlp_result).
type NLPResult struct {
	Objective, Violation, Stationarity float64
	Outer, Inner, Evaluations          int
}

func goNLP(c *C.myb_nlp_result) NLPResult {
	return NLPResult{Objective: float64(c.objective), Violation: float64(c.violation),
		Stationarity: float64(c.stationarity), Outer: int(c.outer), Inner: int(c.inner),
		Evaluations: int(c.evaluations)}
}

// SimsFlanagan is a leg of the direct method (myb_sims_flanagan).
type SimsFlanagan struct {
	MU       float64
	Thrust   float64
	Exhaust  float64
	T0       float64
	Start    State
	Mass0    float64
	T1       float64
	Target   State
	Segments int
	Field    *Model // nil: the conics of MU
}

func (l SimsFlanagan) c(field *C.myb_force_model) C.myb_sims_flanagan {
	out := C.myb_sims_flanagan{mu: C.double(l.MU), thrust: C.double(l.Thrust), exhaust: C.double(l.Exhaust),
		t0: C.double(l.T0), mass0: C.double(l.Mass0), t1: C.double(l.T1), segments: C.size_t(l.Segments),
		field: field}
	for i := range 6 {
		out.start[i] = C.double(l.Start[i])
		out.target[i] = C.double(l.Target[i])
	}
	return out
}

// SimsFlanaganMismatch is the mismatch of the two halves for the throttles (3 per segment),
// the mass at the arrival and the Jacobian of the mismatch (6 x 3 segments)
// (myb_sims_flanagan_mismatch).
func SimsFlanaganMismatch(leg SimsFlanagan, throttles [][3]float64) ([6]float64, float64, []float64, error) {
	const op = "SimsFlanaganMismatch"
	var mismatch [6]float64
	if len(throttles) != leg.Segments {
		return mismatch, 0, nil, refuse(op, -1, "not one throttle for every segment")
	}
	force, release, err := fieldOf(leg.Field)
	if err != nil {
		return mismatch, 0, nil, err
	}
	defer release()
	cl := leg.c(force)
	var mass1 C.double
	jacobian := make([]float64, 18*leg.Segments)
	status := C.myb_sims_flanagan_mismatch(&cl, (*C.double)(first(throttles)),
		(*C.double)(unsafe.Pointer(&mismatch[0])), &mass1, cDoubles(jacobian))
	return mismatch, float64(mass1), jacobian, fieldFailure(op, status, leg.Field)
}

// SimsFlanaganResult is the solution of a leg (myb_sims_flanagan_result).
type SimsFlanaganResult struct {
	Throttles [][3]float64
	Mass1     float64
	Mismatch  [6]float64
	NLP       NLPResult
}

// SimsFlanaganSolve finds the throttles with the least propellant (myb_sims_flanagan_solve);
// throttles: the first guess (nil: the one built in). On MYB_E_NOCONV what was reached comes
// back with the error.
func SimsFlanaganSolve(leg SimsFlanagan, options *NLPOptions, throttles [][3]float64) (*SimsFlanaganResult, error) {
	const op = "SimsFlanaganSolve"
	if throttles != nil && len(throttles) != leg.Segments {
		return nil, refuse(op, -1, "not one throttle for every segment")
	}
	if leg.Segments < 0 || leg.Segments > C.MYB_SIMS_FLANAGAN_MAX_SEGMENTS {
		return nil, refuse(op, -1, "segments out of range")
	}
	force, release, err := fieldOf(leg.Field)
	if err != nil {
		return nil, err
	}
	defer release()
	out := &SimsFlanaganResult{Throttles: make([][3]float64, leg.Segments)}
	copy(out.Throttles, throttles)
	cl := leg.c(force)
	var r C.myb_sims_flanagan_result
	status := C.myb_sims_flanagan_solve(&cl, options.c(), (*C.double)(first(out.Throttles)), nil, &r)
	out.Mass1 = float64(r.mass1)
	for i := range 6 {
		out.Mismatch[i] = float64(r.mismatch[i])
	}
	out.NLP = goNLP(&r.nlp)
	return out, fieldFailure(op, status, leg.Field)
}

// ---- Chains of legs with flybys -------------------------------------------------------------------

// ChainBody is a body of a chain: its provider (asked at Epoch + t), and for a body flown by
// its GM and the least distance of a flyby; ImpulseExhaust > 0 makes the flyby powered (an
// impulse at the periapsis by an engine of that exhaust velocity).
type ChainBody struct {
	State          StateSource
	GM             float64
	Periapsis      float64
	ImpulseExhaust float64
}

// ChainLeg is a leg: its segments and the bounds of its time of flight (s).
type ChainLeg struct {
	Segments int
	Flight   [2]float64
}

// Chain is a chain of legs of Sims and Flanagan joined by flybys (myb_chain): legs + 1 bodies.
type Chain struct {
	MU        float64
	Field     *Model
	Epoch     Time
	Thrust    float64
	Exhaust   float64
	Mass0     float64
	Legs      []ChainLeg
	Bodies    []ChainBody
	Departure [2]float64
	Launch    float64
	Arrival   float64
}

// chainCall describes the chain for C and holds what it needs for the call.
func (c Chain) begin(op string) (*C.myb_chain, []C.mybgo_provider, func(), *keeper, error) {
	if len(c.Legs) == 0 || len(c.Legs) > C.MYB_CHAIN_MAX_LEGS || len(c.Bodies) != len(c.Legs)+1 {
		return nil, nil, nil, nil, refuse(op, -1, "1 .. 8 legs and one body more")
	}
	if !validTime(c.Epoch) {
		return nil, nil, nil, nil, refuse(op, -1, "an epoch out of range")
	}
	k := newKeeper()
	providers := make([]C.mybgo_provider, len(c.Bodies))
	cc := &C.myb_chain{mu: C.double(c.MU), epoch: cTime(c.Epoch), thrust: C.double(c.Thrust),
		exhaust: C.double(c.Exhaust), mass0: C.double(c.Mass0), legs: C.size_t(len(c.Legs)),
		departure: [2]C.double{C.double(c.Departure[0]), C.double(c.Departure[1])},
		launch:    C.double(c.Launch), arrival: C.double(c.Arrival)}
	for i, l := range c.Legs {
		cc.segments[i] = C.size_t(l.Segments)
		cc.flight[i] = [2]C.double{C.double(l.Flight[0]), C.double(l.Flight[1])}
	}
	for i, b := range c.Bodies {
		if b.State == nil {
			k.done()
			return nil, nil, nil, nil, refuse(op, i, "a body without a provider")
		}
		if err := b.State.provider(&providers[i], k); err != nil {
			k.done()
			return nil, nil, nil, nil, err
		}
		cc.bodies[i].gm = C.double(b.GM)
		cc.bodies[i].periapsis = C.double(b.Periapsis)
		cc.impulse_exhaust[i] = C.double(b.ImpulseExhaust)
	}
	force, release, err := fieldOf(c.Field)
	if err != nil {
		k.done()
		return nil, nil, nil, nil, err
	}
	cc.field = force
	return cc, providers, func() { release(); k.done() }, k, nil
}

func (c Chain) failure(op string, status C.myb_status, k *keeper) error {
	if status == C.MYB_OK {
		return nil
	}
	e := &Error{Status: Status(status), Op: op, Index: -1, Err: k.state.err}
	if e.Err == nil && c.Field != nil {
		e.Err = c.Field.state.err
	}
	return e
}

// Size is the count of the unknowns and of the constraints of the chain (myb_chain_size).
func (c Chain) Size() (unknowns, equalities, inequalities int, err error) {
	const op = "Chain.Size"
	cc, providers, release, k, err := c.begin(op)
	if err != nil {
		return 0, 0, 0, err
	}
	defer release()
	var u, e, i C.size_t
	status := C.mybgo_chain_size(cc, &providers[0], &u, &e, &i)
	return int(u), int(e), int(i), c.failure(op, status, k)
}

// Guess is the first guess built in (myb_chain_guess).
func (c Chain) Guess() ([]float64, error) {
	const op = "Chain.Guess"
	u, _, _, err := c.Size()
	if err != nil {
		return nil, err
	}
	cc, providers, release, k, err := c.begin(op)
	if err != nil {
		return nil, err
	}
	defer release()
	x := make([]float64, u)
	status := C.mybgo_chain_guess(cc, &providers[0], cDoubles(x))
	if err := c.failure(op, status, k); err != nil {
		return nil, err
	}
	return x, nil
}

// ChainEvaluation is the propellant of the chain at the unknowns and its constraints, with the
// gradient and the Jacobian (constraints x unknowns) when asked (myb_chain_evaluate).
type ChainEvaluation struct {
	Propellant  float64
	Constraints []float64
	Gradient    []float64
	Jacobian    []float64
}

// Evaluate evaluates the chain at x.
func (c Chain) Evaluate(x []float64, derivatives bool) (*ChainEvaluation, error) {
	const op = "Chain.Evaluate"
	u, e, i, err := c.Size()
	if err != nil {
		return nil, err
	}
	if len(x) != u {
		return nil, refuse(op, -1, "not the count of unknowns of the chain")
	}
	cc, providers, release, k, err := c.begin(op)
	if err != nil {
		return nil, err
	}
	defer release()
	out := &ChainEvaluation{Constraints: make([]float64, e+i)}
	if derivatives {
		out.Gradient = make([]float64, u)
		out.Jacobian = make([]float64, (e+i)*u)
	}
	var propellant C.double
	status := C.mybgo_chain_evaluate(cc, &providers[0], cDoubles(x), &propellant,
		cDoubles(out.Constraints), cDoubles(out.Gradient), cDoubles(out.Jacobian))
	out.Propellant = float64(propellant)
	if err := c.failure(op, status, k); err != nil {
		return nil, err
	}
	return out, nil
}

// ChainResult is the solution of a chain (myb_chain_result), per body: the mass after it, the
// excess speeds coming in and going out, the turn and periapsis of its flyby, the impulse of a
// powered flyby; per leg the mismatch left; the two passes of the solver.
type ChainResult struct {
	X         []float64
	Mass      []float64
	ExcessIn  []float64
	ExcessOut []float64
	Turn      []float64
	Periapsis []float64
	Impulse   []float64
	Mismatch  [][6]float64
	First     NLPResult
	NLP       NLPResult
}

// Solve finds the chain with the least propellant (myb_chain_solve); x: the first guess (nil:
// the one built in). On MYB_E_NOCONV what was reached comes back with the error.
func (c Chain) Solve(options *NLPOptions, x []float64) (*ChainResult, error) {
	const op = "Chain.Solve"
	u, _, _, err := c.Size()
	if err != nil {
		return nil, err
	}
	if x != nil && len(x) != u {
		return nil, refuse(op, -1, "not the count of unknowns of the chain")
	}
	cc, providers, release, k, err := c.begin(op)
	if err != nil {
		return nil, err
	}
	defer release()
	out := &ChainResult{X: make([]float64, u)}
	copy(out.X, x)
	var r C.myb_chain_result
	status := C.mybgo_chain_solve(cc, &providers[0], options.c(), cDoubles(out.X), &r)
	n := len(c.Bodies)
	for i := range n {
		out.Mass = append(out.Mass, float64(r.mass[i]))
		out.ExcessIn = append(out.ExcessIn, float64(r.excess[i][0]))
		out.ExcessOut = append(out.ExcessOut, float64(r.excess[i][1]))
		out.Turn = append(out.Turn, float64(r.turn[i]))
		out.Periapsis = append(out.Periapsis, float64(r.periapsis[i]))
		out.Impulse = append(out.Impulse, float64(r.impulse[i]))
	}
	for l := range len(c.Legs) {
		var m [6]float64
		for i := range 6 {
			m[i] = float64(r.mismatch[l][i])
		}
		out.Mismatch = append(out.Mismatch, m)
	}
	out.First, out.NLP = goNLP(&r.first), goNLP(&r.nlp)
	return out, c.failure(op, status, k)
}

// ChainFlyby is the impulse and the turn of a powered flyby for the excess speeds and the
// periapsis, with their derivatives (myb_chain_flyby: d dv / d(in, out, rp), then d turn).
func ChainFlyby(gm, vIn, vOut, periapsis float64) (dv, turn float64, partials [6]float64, err error) {
	var cdv, cturn C.double
	status := C.myb_chain_flyby(C.double(gm), C.double(vIn), C.double(vOut), C.double(periapsis),
		&cdv, &cturn, (*C.double)(unsafe.Pointer(&partials[0])))
	return float64(cdv), float64(cturn), partials, fail("ChainFlyby", status)
}
