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
