// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include "glue.h"
*/
import "C"

// StateSource is a provider of states as functions of time (myb_state_function of
// myb/kin.h), described by data: the state of a body of a model, a guided body of a system,
// a body of a chain of transfers, a participant of the light time, an object of a close
// approach. SPKBody and KeplerOrbit are evaluated in C; GoState is a Go function called back
// through cgo — the slow path.
type StateSource interface {
	provider(p *C.mybgo_provider, k *keeper) error
}

// SPKBody is the state of a body from an ephemeris: target relative to center; Frame != 0
// requires the segments to be in that frame (NAIF code, 1 = J2000).
type SPKBody struct {
	Ephemeris      *Ephemeris
	Target, Center int32
	Frame          int32
}

// KeplerOrbit is the two-body orbit about a centre of gravitational parameter MU through
// State at the instant Epoch: the state at t is myb_kepler_propagate over the seconds from
// Epoch to t (a planet on its orbit, as a C program would write it).
type KeplerOrbit struct {
	MU    float64
	Epoch Time
	State State
}

// StateFunc is the state (m, m/s) of a point at the instant t. A non-nil error stops the call
// of the library going on.
type StateFunc func(t Time) (State, error)

// GoState is a provider in Go, called back through cgo at every instant asked: the slow path.
type GoState struct {
	Func StateFunc
}

func (s SPKBody) provider(p *C.mybgo_provider, k *keeper) error {
	spk, err := k.hold(s.Ephemeris.obj())
	if err != nil {
		return err
	}
	p.kind = C.MYBGO_PROVIDER_SPK
	p.spk = (*C.myb_spk)(spk)
	p.target = C.int32_t(s.Target)
	p.center = C.int32_t(s.Center)
	p.frame = C.int32_t(s.Frame)
	return nil
}

func (o KeplerOrbit) provider(p *C.mybgo_provider, k *keeper) error {
	if !validTime(o.Epoch) {
		return refuse("StateSource", -1, "an epoch out of range")
	}
	p.kind = C.MYBGO_PROVIDER_KEPLER
	p.mu = C.double(o.MU)
	for i := range 6 {
		p.state[i] = C.double(o.State[i])
	}
	p.epoch = cTime(o.Epoch)
	return nil
}

func (g GoState) provider(p *C.mybgo_provider, k *keeper) error {
	if g.Func == nil {
		return refuse("StateSource", -1, "a Go provider without a function")
	}
	p.kind = C.MYBGO_PROVIDER_GO
	p.handle = k.handle(&stateCallback{f: g.Func, state: k.state})
	return nil
}

// describeProvider fills a provider record; nil is no provider (the origin of a model).
func describeProvider(src StateSource, p *C.mybgo_provider, k *keeper) error {
	if src == nil {
		p.kind = C.MYBGO_PROVIDER_NONE
		return nil
	}
	return src.provider(p, k)
}

// States are the states of a provider at the instants, in one call (a loop over the provider
// in C).
func States(src StateSource, t []Time) ([]State, error) {
	const op = "States"
	if src == nil {
		return nil, refuse(op, -1, "no provider")
	}
	if err := checkTimes(op, t); err != nil {
		return nil, err
	}
	k := newKeeper()
	defer k.done()
	var p C.mybgo_provider
	if err := src.provider(&p, k); err != nil {
		return nil, err
	}
	out := make([]State, len(t))
	var failed C.size_t
	status := C.mybgo_provider_states(&p, cTimes(t), C.size_t(len(t)), cStates(out), &failed)
	if err := k.failure(op, status, int(failed)); err != nil {
		return nil, err
	}
	return out, nil
}
