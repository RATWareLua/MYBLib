// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include "glue.h"
*/
import "C"

import "unsafe"

// ApproachInfo is a closest approach of two objects (myb_approach).
type ApproachInfo struct {
	T        float64
	Distance float64
	Position [3]float64 // of the second from the first
	Velocity [3]float64
}

func goApproach(a *C.myb_approach) ApproachInfo {
	return ApproachInfo{T: float64(a.t), Distance: float64(a.distance), Position: goVec3(a.position),
		Velocity: goVec3(a.velocity)}
}

// FindApproaches are the closest approaches of two providers between t0 and t1 (seconds, as
// myb_time_from_seconds), looked for in steps of `step` (myb_approach_find), in one call.
func FindApproaches(first, second StateSource, t0, t1, step float64) ([]ApproachInfo, error) {
	const op = "FindApproaches"
	if first == nil || second == nil {
		return nil, refuse(op, -1, "no provider")
	}
	k := newKeeper()
	defer k.done()
	var a, b C.mybgo_provider
	if err := first.provider(&a, k); err != nil {
		return nil, err
	}
	if err := second.provider(&b, k); err != nil {
		return nil, err
	}
	var found *C.myb_approach
	var count C.size_t
	status := C.mybgo_approach_find(&a, &b, C.double(t0), C.double(t1), C.double(step), &found, &count)
	defer C.mybgo_free(unsafe.Pointer(found))
	if err := k.failure(op, status, -1); err != nil {
		return nil, err
	}
	out := make([]ApproachInfo, int(count))
	for i, v := range unsafe.Slice(found, int(count)) {
		out[i] = goApproach(&v)
	}
	return out, nil
}

// Encounter is an encounter in the plane across the relative velocity (myb_encounter).
type Encounter struct {
	Miss   [2]float64
	Sigma  [2]float64
	Radius float64
}

func (e Encounter) c() C.myb_encounter {
	return C.myb_encounter{miss: [2]C.double{C.double(e.Miss[0]), C.double(e.Miss[1])},
		sigma: [2]C.double{C.double(e.Sigma[0]), C.double(e.Sigma[1])}, radius: C.double(e.Radius)}
}

func goEncounter(c *C.myb_encounter) Encounter {
	return Encounter{Miss: [2]float64{float64(c.miss[0]), float64(c.miss[1])},
		Sigma: [2]float64{float64(c.sigma[0]), float64(c.sigma[1])}, Radius: float64(c.radius)}
}

// EncounterPlane is the encounter of the relative position and velocity at the closest
// approach with the covariance of the relative position (3 x 3) (myb_encounter_plane).
func EncounterPlane(position, velocity [3]float64, covariance [9]float64, radius float64) (Encounter, error) {
	var c C.myb_encounter
	status := C.myb_encounter_plane((*C.double)(unsafe.Pointer(&position[0])),
		(*C.double)(unsafe.Pointer(&velocity[0])), (*C.double)(unsafe.Pointer(&covariance[0])),
		C.double(radius), &c)
	return goEncounter(&c), fail("EncounterPlane", status)
}

// CollisionProbability is the integral of the Gaussian over the disk (myb_collision_probability).
func CollisionProbability(e Encounter) (float64, error) {
	c := e.c()
	var p C.double
	status := C.myb_collision_probability(&c, &p)
	return float64(p), fail("CollisionProbability", status)
}

// CollisionProbabilityChan is the series of Chan (myb_collision_probability_chan; NaN for an
// invalid encounter).
func CollisionProbabilityChan(e Encounter) float64 {
	c := e.c()
	return float64(C.myb_collision_probability_chan(&c))
}

// CollisionProbabilityMax is the largest probability for a covariance known in shape only,
// and the factor at which it is reached (myb_collision_probability_max).
func CollisionProbabilityMax(e Encounter) (float64, float64, error) {
	c := e.c()
	var p, factor C.double
	status := C.myb_collision_probability_max(&c, &p, &factor)
	return float64(p), float64(factor), fail("CollisionProbabilityMax", status)
}

// CollisionRate is the rate of entry into the sphere of radius for the relative state
// (myb_collision_rate).
func CollisionRate(mean State, covariance [36]float64, radius float64) (float64, error) {
	var rate C.double
	status := C.myb_collision_rate((*C.double)(unsafe.Pointer(&mean)),
		(*C.double)(unsafe.Pointer(&covariance[0])), C.double(radius), &rate)
	return float64(rate), fail("CollisionRate", status)
}

// GaussianFunc is the relative state as a Gaussian at the time t: its mean and covariance.
type GaussianFunc func(t float64) (State, [36]float64, error)

// CollisionProbabilityLong is the expected number of entries over [t0, t1] of the relative
// Gaussian a Go function gives (myb_collision_probability_long) — the slow path of callbacks;
// and the estimate of its error.
func CollisionProbabilityLong(f GaussianFunc, t0, t1, step, radius float64) (float64, float64, error) {
	const op = "CollisionProbabilityLong"
	if f == nil {
		return 0, 0, refuse(op, -1, "no function")
	}
	k := newKeeper()
	defer k.done()
	h := k.handle(&gaussianCallback{f: f, state: k.state})
	var p, e C.double
	status := C.mybgo_collision_long(h, C.double(t0), C.double(t1), C.double(step), C.double(radius), &p, &e)
	return float64(p), float64(e), k.failure(op, status, -1)
}

// Conjunction is an approach with the uncertainty carried to it (myb_conjunction).
type Conjunction struct {
	Approach    ApproachInfo
	Covariance  [9]float64 // of the relative position
	Encounter   Encounter
	Probability float64
}

// FindConjunctions propagates two orbits (each of its model, started at one time, with the
// covariance of its state) to t1 and gives their approaches with the probability of collision
// for the sum of the radii (myb_conjunction_find), in one call.
func FindConjunctions(first *Model, a OrbitStart, second *Model, b OrbitStart, radius, t1, step float64) ([]Conjunction, error) {
	const op = "FindConjunctions"
	ma, mb, release, err := twoModels(first, second)
	if err != nil {
		return nil, err
	}
	defer release()
	ca, cb := a.c(), b.c()
	var found *C.myb_conjunction
	var count C.size_t
	status := C.mybgo_conjunctions(ma, &ca, mb, &cb, C.double(radius), C.double(t1), C.double(step),
		&found, &count, nil, nil)
	defer C.mybgo_free(unsafe.Pointer(found))
	if err := modelsFailure(op, status, first, second); err != nil {
		return nil, err
	}
	out := make([]Conjunction, int(count))
	for i, v := range unsafe.Slice(found, int(count)) {
		out[i] = Conjunction{Approach: goApproach(&v.approach), Encounter: goEncounter(&v.encounter),
			Probability: float64(v.probability)}
		for j := range 9 {
			out[i].Covariance[j] = float64(v.covariance[j])
		}
	}
	return out, nil
}

// LongEncounter is the expected number of entries of the second into the sphere about the
// first over the propagation of both to t1 (myb_conjunction_long), and its error estimate.
func LongEncounter(first *Model, a OrbitStart, second *Model, b OrbitStart, radius, t1, step float64) (float64, float64, error) {
	const op = "LongEncounter"
	ma, mb, release, err := twoModels(first, second)
	if err != nil {
		return 0, 0, err
	}
	defer release()
	ca, cb := a.c(), b.c()
	var count C.size_t
	var p, e C.double
	status := C.mybgo_conjunctions(ma, &ca, mb, &cb, C.double(radius), C.double(t1), C.double(step),
		nil, &count, &p, &e)
	return float64(p), float64(e), modelsFailure(op, status, first, second)
}

// twoModels holds two models (one when they are the same) for a call.
func twoModels(first, second *Model) (*C.mybgo_model, *C.mybgo_model, func(), error) {
	if second == first {
		ma, err := first.begin()
		if err != nil {
			return nil, nil, nil, err
		}
		return ma, ma, first.end, nil
	}
	// the locks in one order whatever the order of the arguments, so that two calls given the
	// same models the other way round cannot wait for each other (Go does not move objects of
	// the heap, so their addresses order them)
	lo, hi := first, second
	if uintptr(unsafe.Pointer(second)) < uintptr(unsafe.Pointer(first)) {
		lo, hi = second, first
	}
	pl, err := lo.begin()
	if err != nil {
		return nil, nil, nil, err
	}
	ph, err := hi.begin()
	if err != nil {
		lo.end()
		return nil, nil, nil, err
	}
	if lo == first {
		return pl, ph, func() { hi.end(); lo.end() }, nil
	}
	return ph, pl, func() { hi.end(); lo.end() }, nil
}

func modelsFailure(op string, status C.myb_status, first, second *Model) error {
	if status == C.MYB_OK {
		return nil
	}
	cause := first.state.err
	if cause == nil {
		cause = second.state.err
	}
	return &Error{Status: Status(status), Op: op, Index: -1, Err: cause}
}
