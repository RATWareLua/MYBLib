// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include "glue.h"
*/
import "C"

import "unsafe"

// ---- Attitude laws -------------------------------------------------------------------------------

// AttitudeLaw is a law of the attitude of an object (myb_attitude_function of myb/surface.h):
// the matrix (rows) from the inertial axes of a model to the body frame.
type AttitudeLaw interface {
	attitude(a *C.mybgo_attitude)
}

// AttitudeFixed is a constant attitude (myb_attitude_fixed).
type AttitudeFixed struct {
	M [9]float64
}

// AttitudeOrbital are the axes of the orbit about a body of the model (myb_attitude_orbital).
type AttitudeOrbital struct {
	Body int
}

// AttitudeSunNadir is the yaw steering towards Body with the panels across the plane of the
// source (myb_attitude_sun_nadir).
type AttitudeSunNadir struct {
	Body, Source int
}

// AttitudeSpin is a uniform rotation about the body z axis (myb_attitude_spin).
type AttitudeSpin struct {
	M0          [9]float64
	Rate, Phase float64
}

func (l AttitudeFixed) attitude(a *C.mybgo_attitude) {
	a.kind = C.MYBGO_ATTITUDE_FIXED
	for i := range 9 {
		a.m[i] = C.double(l.M[i])
	}
}

func (l AttitudeOrbital) attitude(a *C.mybgo_attitude) {
	a.kind, a.body = C.MYBGO_ATTITUDE_ORBITAL, C.size_t(l.Body)
}

func (l AttitudeSunNadir) attitude(a *C.mybgo_attitude) {
	a.kind, a.body, a.source = C.MYBGO_ATTITUDE_SUN_NADIR, C.size_t(l.Body), C.size_t(l.Source)
}

func (l AttitudeSpin) attitude(a *C.mybgo_attitude) {
	a.kind, a.rate, a.phase = C.MYBGO_ATTITUDE_SPIN, C.double(l.Rate), C.double(l.Phase)
	for i := range 9 {
		a.m[i] = C.double(l.M0[i])
	}
}

func attitudes(op string, model *C.mybgo_model, law AttitudeLaw, t []float64, y []State) ([][9]float64, error) {
	if law == nil {
		return nil, refuse(op, -1, "no law")
	}
	if len(t) != len(y) {
		return nil, refuse(op, -1, "times and states of different lengths")
	}
	var a C.mybgo_attitude
	law.attitude(&a)
	out := make([][9]float64, len(t))
	var failed C.size_t
	status := C.mybgo_attitudes(model, &a, cDoubles(t), cStates(y), C.size_t(len(t)),
		(*C.double)(first(out)), &failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// Attitudes samples a law at the times and states of a model (the laws that follow its bodies
// need its environment), in one call.
func (m *Model) Attitudes(law AttitudeLaw, t []float64, y []State) ([][9]float64, error) {
	model, err := m.begin()
	if err != nil {
		return nil, err
	}
	defer m.end()
	return attitudes("Model.Attitudes", model, law, t, y)
}

// Attitudes samples a law that needs no environment (AttitudeFixed, AttitudeSpin), in one call.
func Attitudes(law AttitudeLaw, t []float64, y []State) ([][9]float64, error) {
	return attitudes("Attitudes", nil, law, t, y)
}

// FreeRotation is the torque-free rotation of a body of the principal moments from the
// attitude q0 and the body angular momentum l0 at t = 0, at the times
// (myb_free_rotation_init, myb_free_rotation_state), in one call.
func FreeRotation(moments [3]float64, q0 [4]float64, l0 [3]float64, t []float64) ([][4]float64, [][3]float64, error) {
	q := make([][4]float64, len(t))
	l := make([][3]float64, len(t))
	status := C.mybgo_free_rotation((*C.double)(unsafe.Pointer(&moments[0])),
		(*C.double)(unsafe.Pointer(&q0[0])), (*C.double)(unsafe.Pointer(&l0[0])), cDoubles(t),
		C.size_t(len(t)), (*C.double)(first(q)), (*C.double)(first(l)))
	if err := fail("FreeRotation", status); err != nil {
		return nil, nil, err
	}
	return q, l, nil
}

// ---- Observables -----------------------------------------------------------------------------

// LightBody is a body that delays and deflects light.
type LightBody struct {
	GM    float64
	State StateSource
}

// LightModel is the model of the light time (myb_light_model): the speed of light, the PPN
// gamma and the bodies (none: Newtonian light time).
type LightModel struct {
	C      float64
	Gamma  float64
	Bodies []LightBody
}

// LightDirection is the sense of the light (myb_light_direction).
type LightDirection int

// The two senses.
const (
	Reception    LightDirection = C.MYB_LIGHT_RECEPTION
	Transmission LightDirection = C.MYB_LIGHT_TRANSMISSION
)

// LightPath is a solution of the light time (myb_light_path).
type LightPath struct {
	TTarget    Time
	Target     State
	LightTime  float64
	Range      float64
	Delay      float64
	RangeRate  float64
	Iterations int
}

func (l LightModel) c(k *keeper) (C.mybgo_light_model, []float64, []C.mybgo_provider, error) {
	gm := make([]float64, len(l.Bodies))
	providers := make([]C.mybgo_provider, len(l.Bodies))
	for i, b := range l.Bodies {
		if b.State == nil {
			return C.mybgo_light_model{}, nil, nil, refuse("LightModel", i, "a body without a provider")
		}
		if err := b.State.provider(&providers[i], k); err != nil {
			return C.mybgo_light_model{}, nil, nil, err
		}
		gm[i] = b.GM
	}
	return C.mybgo_light_model{c: C.double(l.C), gamma: C.double(l.Gamma)}, gm, providers, nil
}

// LightTimes solves the light time between the observers (state observers[k] at t[k]) and the
// target (myb_light_time), in one call.
func LightTimes(model LightModel, target StateSource, t []Time, observers []State, direction LightDirection) ([]LightPath, error) {
	const op = "LightTimes"
	if target == nil {
		return nil, refuse(op, -1, "no target")
	}
	if len(t) != len(observers) {
		return nil, refuse(op, -1, "times and observers of different lengths")
	}
	if err := checkTimes(op, t); err != nil {
		return nil, err
	}
	k := newKeeper()
	defer k.done()
	cm, gm, bodies, err := model.c(k)
	if err != nil {
		return nil, err
	}
	var aim C.mybgo_provider
	if err := target.provider(&aim, k); err != nil {
		return nil, err
	}
	paths := make([]C.myb_light_path, len(t))
	var failed C.size_t
	status := C.mybgo_light_times(&cm, cDoubles(gm), (*C.mybgo_provider)(first(bodies)),
		C.size_t(len(bodies)), &aim, cTimes(t), cStates(observers), C.size_t(len(t)), C.int(direction),
		(*C.myb_light_path)(first(paths)), &failed)
	if err := k.failure(op, status, int(failed)); err != nil {
		return nil, err
	}
	out := make([]LightPath, len(paths))
	for i, p := range paths {
		out[i] = LightPath{TTarget: Time{S: int64(p.t_target.s), F: float64(p.t_target.f)},
			LightTime: float64(p.light_time), Range: float64(p._range), Delay: float64(p.delay),
			RangeRate: float64(p.range_rate), Iterations: int(p.iterations)}
		for j := range 6 {
			out[i].Target[j] = float64(p.target[j])
		}
	}
	return out, nil
}

// RoundTrip is a round trip received at a station (myb_round_trip): the total time, the
// instants of the bounce and of the transmission.
type RoundTrip struct {
	Total    float64
	Bounce   Time
	Transmit Time
}

// RoundTrips solves the round trips received by the station at the instants (myb_round_trip),
// in one call.
func RoundTrips(model LightModel, station, target StateSource, t []Time) ([]RoundTrip, error) {
	const op = "RoundTrips"
	if station == nil || target == nil {
		return nil, refuse(op, -1, "no station or target")
	}
	if err := checkTimes(op, t); err != nil {
		return nil, err
	}
	k := newKeeper()
	defer k.done()
	cm, gm, bodies, err := model.c(k)
	if err != nil {
		return nil, err
	}
	var from, aim C.mybgo_provider
	if err := station.provider(&from, k); err != nil {
		return nil, err
	}
	if err := target.provider(&aim, k); err != nil {
		return nil, err
	}
	total := make([]float64, len(t))
	bounce := make([]Time, len(t))
	transmit := make([]Time, len(t))
	var failed C.size_t
	status := C.mybgo_round_trips(&cm, cDoubles(gm), (*C.mybgo_provider)(first(bodies)),
		C.size_t(len(bodies)), &from, &aim, cTimes(t), C.size_t(len(t)), cDoubles(total),
		cTimes(bounce), cTimes(transmit), &failed)
	if err := k.failure(op, status, int(failed)); err != nil {
		return nil, err
	}
	out := make([]RoundTrip, len(t))
	for i := range out {
		out[i] = RoundTrip{Total: total[i], Bounce: bounce[i], Transmit: transmit[i]}
	}
	return out, nil
}

// Aberrations are the directions p[k] as seen by observers of velocities v[k] (myb_aberration),
// in one call.
func Aberrations(p, v [][3]float64, c float64) ([][3]float64, error) {
	if len(p) != len(v) {
		return nil, refuse("Aberrations", -1, "slices of different lengths")
	}
	out := make([][3]float64, len(p))
	C.mybgo_aberrations((*C.double)(first(p)), (*C.double)(first(v)), C.size_t(len(p)), C.double(c),
		(*C.double)(first(out)))
	return out, nil
}

// Deflections are the directions p[k] deflected by a body (q: body -> source, e: body ->
// observer, em: the distances body - observer) (myb_light_deflection), in one call.
func Deflections(gm, c, gamma float64, p, q, e [][3]float64, em []float64, dlim float64) ([][3]float64, error) {
	if len(q) != len(p) || len(e) != len(p) || len(em) != len(p) {
		return nil, refuse("Deflections", -1, "slices of different lengths")
	}
	out := make([][3]float64, len(p))
	C.mybgo_deflections(C.double(gm), C.double(c), C.double(gamma), (*C.double)(first(p)),
		(*C.double)(first(q)), (*C.double)(first(e)), cDoubles(em), C.size_t(len(p)), C.double(dlim),
		(*C.double)(first(out)))
	return out, nil
}

// ShapiroDelays are the relativistic delays of the rays x1[k] -> x2[k] by a body at xb[k]
// (myb_shapiro_delay), in one call.
func ShapiroDelays(gm, c, gamma float64, x1, x2, xb [][3]float64) ([]float64, error) {
	if len(x2) != len(x1) || len(xb) != len(x1) {
		return nil, refuse("ShapiroDelays", -1, "slices of different lengths")
	}
	out := make([]float64, len(x1))
	C.mybgo_shapiro_delays(C.double(gm), C.double(c), C.double(gamma), (*C.double)(first(x1)),
		(*C.double)(first(x2)), (*C.double)(first(xb)), C.size_t(len(x1)), cDoubles(out))
	return out, nil
}

// RaDec are the right ascensions and declinations of the directions (myb_direction_radec).
func RaDec(p [][3]float64) ([]float64, []float64) {
	ra := make([]float64, len(p))
	dec := make([]float64, len(p))
	C.mybgo_radecs((*C.double)(first(p)), C.size_t(len(p)), cDoubles(ra), cDoubles(dec))
	return ra, dec
}

// PlasmaProfile is a power-law electron density about a body (myb_plasma_profile).
type PlasmaProfile struct {
	Radius  float64
	Density []float64 // 1 .. 8 terms
	Power   []float64
}

// PlasmaColumns are the electron columns of the profile along the segments p1[k] -> p2[k]
// (myb_plasma_column), in one call.
func PlasmaColumns(profile PlasmaProfile, p1, p2 [][3]float64) ([]float64, error) {
	const op = "PlasmaColumns"
	if len(p1) != len(p2) {
		return nil, refuse(op, -1, "slices of different lengths")
	}
	if len(profile.Density) == 0 || len(profile.Density) > C.MYB_PLASMA_TERMS ||
		len(profile.Power) != len(profile.Density) {
		return nil, refuse(op, -1, "1 .. 8 terms, a power for every density")
	}
	cp := C.myb_plasma_profile{radius: C.double(profile.Radius), count: C.size_t(len(profile.Density))}
	for i := range profile.Density {
		cp.density[i] = C.double(profile.Density[i])
		cp.power[i] = C.double(profile.Power[i])
	}
	out := make([]float64, len(p1))
	var failed C.size_t
	status := C.mybgo_plasma_columns(&cp, (*C.double)(first(p1)), (*C.double)(first(p2)),
		C.size_t(len(p1)), cDoubles(out), &failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	return out, nil
}
