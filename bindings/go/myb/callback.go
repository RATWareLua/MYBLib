// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include "glue.h"
*/
import "C"

import (
	"runtime/cgo"
	"unsafe"
)

// TermFunc is the acceleration (m/s^2) of a term of a model at the model time t (s), the
// position x (m) and the velocity v (m/s) of the object. A non-nil error stops the call of
// the library going on.
type TermFunc func(t float64, x, v [3]float64) ([3]float64, error)

// EventFunc is an event function: its roots in time are the events. It has no error: a value
// that is NaN is no change of sign.
type EventFunc func(t float64, y State) float64

// callState is what the Go functions of one call record: the first error returned. A model
// serves one call at a time, and the library calls back on the goroutine of the call.
type callState struct {
	err error
}

type termCallback struct {
	f     TermFunc
	state *callState
}

type eventCallback struct {
	f EventFunc
}

//export mybgoTerm
func mybgoTerm(handle C.uintptr_t, t C.double, x, v, a *C.double) C.int {
	cb := cgo.Handle(handle).Value().(*termCallback)
	acceleration, err := cb.f(float64(t), *(*[3]float64)(unsafe.Pointer(x)),
		*(*[3]float64)(unsafe.Pointer(v)))
	if err != nil {
		if cb.state.err == nil {
			cb.state.err = err
		}
		return C.int(statusOf(err))
	}
	*(*[3]float64)(unsafe.Pointer(a)) = acceleration
	return C.MYB_OK
}

//export mybgoEvent
func mybgoEvent(handle C.uintptr_t, t C.double, y *C.double) C.double {
	cb := cgo.Handle(handle).Value().(*eventCallback)
	return C.double(cb.f(float64(t), *(*State)(unsafe.Pointer(y))))
}

// record keeps the first error of a call and gives its status.
func (s *callState) record(err error) C.int {
	if s.err == nil {
		s.err = err
	}
	return C.int(statusOf(err))
}

type stateCallback struct {
	f     StateFunc
	state *callState
}

//export mybgoState
func mybgoState(handle C.uintptr_t, s C.int64_t, f C.double, state *C.double) C.int {
	cb := cgo.Handle(handle).Value().(*stateCallback)
	y, err := cb.f(Time{S: int64(s), F: float64(f)})
	if err != nil {
		return cb.state.record(err)
	}
	*(*State)(unsafe.Pointer(state)) = y
	return C.MYB_OK
}

type actionCallback struct {
	f     ActionFunc
	state *callState
}

//export mybgoAction
func mybgoAction(handle C.uintptr_t, t C.double, view *C.myb_system_view, acceleration,
	torque *C.double) C.int {
	cb := cgo.Handle(handle).Value().(*actionCallback)
	n := int(view.bodies)
	v := SystemView{Body: int(view.body),
		X:  unsafe.Slice((*[3]float64)(unsafe.Pointer(view.x)), n),
		V:  unsafe.Slice((*[3]float64)(unsafe.Pointer(view.v)), n),
		GM: unsafe.Slice((*float64)(unsafe.Pointer(view.gm)), n)}
	if view.frame != nil {
		v.Frame = (*[9]float64)(unsafe.Pointer(view.frame))
	}
	a, w, err := cb.f(float64(t), v)
	if err != nil {
		return cb.state.record(err)
	}
	*(*[3]float64)(unsafe.Pointer(acceleration)) = a
	*(*[3]float64)(unsafe.Pointer(torque)) = w
	return C.MYB_OK
}

type shootCallback struct {
	f     ShootFunc
	n, m  int
	state *callState
}

//export mybgoShoot
func mybgoShoot(handle C.uintptr_t, u, f, jacobian *C.double) C.int {
	cb := cgo.Handle(handle).Value().(*shootCallback)
	var jac []float64
	if jacobian != nil {
		jac = unsafe.Slice((*float64)(unsafe.Pointer(jacobian)), cb.n*cb.m)
	}
	err := cb.f(unsafe.Slice((*float64)(unsafe.Pointer(u)), cb.n),
		unsafe.Slice((*float64)(unsafe.Pointer(f)), cb.m), jac)
	if err != nil {
		return cb.state.record(err)
	}
	return C.MYB_OK
}

type gaussianCallback struct {
	f     GaussianFunc
	state *callState
}

//export mybgoGaussian
func mybgoGaussian(handle C.uintptr_t, t C.double, mean, covariance *C.double) C.int {
	cb := cgo.Handle(handle).Value().(*gaussianCallback)
	m, p, err := cb.f(float64(t))
	if err != nil {
		return cb.state.record(err)
	}
	*(*State)(unsafe.Pointer(mean)) = m
	*(*[36]float64)(unsafe.Pointer(covariance)) = p
	return C.MYB_OK
}
