// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include "glue.h"
*/
import "C"

import "unsafe"

// State is a position (m) and a velocity (m/s): x, y, z, vx, vy, vz (myb_pv).
type State [6]float64

// Elements are the classical elements (myb_elements): the semi-latus rectum P (m), the
// eccentricity, the inclination, the node, the argument of periapsis and the true anomaly
// (rad).
type Elements struct {
	P, E, I, RAAN, ArgP, Nu float64
}

// Equinoctial are the modified equinoctial elements (myb_equinoctial).
type Equinoctial struct {
	P, F, G, H, K, L float64
}

func cPV(s *State) *C.myb_pv {
	return (*C.myb_pv)(unsafe.Pointer(s))
}

func cPVs(s []State) *C.myb_pv {
	return (*C.myb_pv)(first(s))
}

// StateFromElements is the state of a set of elements about a body of gravitational
// parameter mu (myb_pv_from_elements).
func StateFromElements(mu float64, e Elements) (State, error) {
	var s State
	status := C.myb_pv_from_elements(C.double(mu), (*C.myb_elements)(unsafe.Pointer(&e)), cPV(&s))
	return s, fail("StateFromElements", status)
}

// StatesFromElements is StateFromElements over a slice, in one call.
func StatesFromElements(mu float64, e []Elements) ([]State, error) {
	out := make([]State, len(e))
	var failed C.size_t
	status := C.myb_pv_from_elements_many(C.double(mu), (*C.myb_elements)(first(e)),
		C.size_t(len(e)), cPVs(out), &failed)
	if err := failAt("StatesFromElements", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// ElementsFromState is the set of elements of a state (myb_elements_from_pv).
func ElementsFromState(mu float64, s State) (Elements, error) {
	var e Elements
	status := C.myb_elements_from_pv(C.double(mu), cPV(&s), (*C.myb_elements)(unsafe.Pointer(&e)))
	return e, fail("ElementsFromState", status)
}

// ElementsFromStates is ElementsFromState over a slice, in one call.
func ElementsFromStates(mu float64, s []State) ([]Elements, error) {
	out := make([]Elements, len(s))
	var failed C.size_t
	status := C.myb_elements_from_pv_many(C.double(mu), cPVs(s), C.size_t(len(s)),
		(*C.myb_elements)(first(out)), &failed)
	if err := failAt("ElementsFromStates", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// EquinoctialFromState is the set of modified equinoctial elements of a state
// (myb_equinoctial_from_pv).
func EquinoctialFromState(mu float64, s State) (Equinoctial, error) {
	var e Equinoctial
	status := C.myb_equinoctial_from_pv(C.double(mu), cPV(&s),
		(*C.myb_equinoctial)(unsafe.Pointer(&e)))
	return e, fail("EquinoctialFromState", status)
}

// EquinoctialFromStates is EquinoctialFromState over a slice, in one call.
func EquinoctialFromStates(mu float64, s []State) ([]Equinoctial, error) {
	out := make([]Equinoctial, len(s))
	var failed C.size_t
	status := C.myb_equinoctial_from_pv_many(C.double(mu), cPVs(s), C.size_t(len(s)),
		(*C.myb_equinoctial)(first(out)), &failed)
	if err := failAt("EquinoctialFromStates", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// StateFromEquinoctial is the state of a set of modified equinoctial elements
// (myb_pv_from_equinoctial).
func StateFromEquinoctial(mu float64, e Equinoctial) (State, error) {
	var s State
	status := C.myb_pv_from_equinoctial(C.double(mu), (*C.myb_equinoctial)(unsafe.Pointer(&e)),
		cPV(&s))
	return s, fail("StateFromEquinoctial", status)
}

// StatesFromEquinoctial is StateFromEquinoctial over a slice, in one call.
func StatesFromEquinoctial(mu float64, e []Equinoctial) ([]State, error) {
	out := make([]State, len(e))
	var failed C.size_t
	status := C.myb_pv_from_equinoctial_many(C.double(mu), (*C.myb_equinoctial)(first(e)),
		C.size_t(len(e)), cPVs(out), &failed)
	if err := failAt("StatesFromEquinoctial", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// KeplerPropagate is the state dt seconds after s0 on its two-body orbit
// (myb_kepler_propagate).
func KeplerPropagate(mu float64, s0 State, dt float64) (State, error) {
	var s State
	status := C.myb_kepler_propagate(C.double(mu), cPV(&s0), C.double(dt), cPV(&s))
	return s, fail("KeplerPropagate", status)
}

// KeplerPropagateMany is KeplerPropagate to every offset of a slice, in one call
// (myb_kepler_propagate_many).
func KeplerPropagateMany(mu float64, s0 State, dt []float64) ([]State, error) {
	out := make([]State, len(dt))
	var failed C.size_t
	status := C.myb_kepler_propagate_many(C.double(mu), cPV(&s0), cDoubles(dt), C.size_t(len(dt)),
		cPVs(out), &failed)
	if err := failAt("KeplerPropagateMany", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// KeplerSTM is KeplerPropagate with the state transition matrix, stm[6 i + j] = d y_i / d
// y0_j (myb_kepler_stm).
func KeplerSTM(mu float64, s0 State, dt float64) (State, [36]float64, error) {
	var s State
	var stm [36]float64
	status := C.myb_kepler_stm(C.double(mu), cPV(&s0), C.double(dt), cPV(&s),
		(*C.double)(unsafe.Pointer(&stm[0])))
	return s, stm, fail("KeplerSTM", status)
}

// KeplerSTMMany is KeplerSTM to every offset of a slice, in one call (myb_kepler_stm_many).
func KeplerSTMMany(mu float64, s0 State, dt []float64) ([]State, [][36]float64, error) {
	out := make([]State, len(dt))
	stm := make([][36]float64, len(dt))
	var failed C.size_t
	status := C.myb_kepler_stm_many(C.double(mu), cPV(&s0), cDoubles(dt), C.size_t(len(dt)),
		cPVs(out), (*C.double)(first(stm)), &failed)
	if err := failAt("KeplerSTMMany", status, failed); err != nil {
		return nil, nil, err
	}
	return out, stm, nil
}
