// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include "glue.h"
*/
import "C"

import "slices"

// cObservations are the records of the observations, the count of their scalar values and the
// count of each.
func cObservations(op string, observations []Observation) ([]C.mybgo_observation, int, []int, error) {
	out := make([]C.mybgo_observation, len(observations))
	counts := make([]int, len(observations))
	values := 0
	for i, o := range observations {
		n := o.Kind.values()
		if n == 0 || len(o.Value) != n || len(o.Sigma) != n {
			return nil, 0, nil, refuse(op, i, "an observation without the values its kind needs")
		}
		d := &out[i]
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
		counts[i] = n
		values += n
	}
	return out, values, counts, nil
}

func cStation(s Station) C.mybgo_station {
	return C.mybgo_station{body: C.size_t(s.Body),
		ellipsoid: C.myb_ellipsoid{a: C.double(s.Ellipsoid.A), f: C.double(s.Ellipsoid.F)},
		longitude: C.double(s.Longitude), latitude: C.double(s.Latitude), height: C.double(s.Height)}
}

func cStations(stations []Station) []C.mybgo_station {
	out := make([]C.mybgo_station, len(stations))
	for i, s := range stations {
		out[i] = cStation(s)
	}
	return out
}

func cLight(l *LightTime) *C.mybgo_light {
	if l == nil {
		return nil
	}
	return &C.mybgo_light{c: C.double(l.C), gm: C.double(l.GM), two_way: flag(l.TwoWay),
		aberration: flag(l.Aberration)}
}

// FilterRequest is a run of the extended Kalman filter (myb_filter) from the state Y0 at T0
// with the covariance of the 6 + np unknowns (row-major) towards TEnd over the observations
// (in time order), with the white acceleration ProcessNoise (m^2 s^-3) between them; Bias are
// the first values of the systematic shifts (their count is that of the shifts); Smooth runs
// the smoother (Rauch, Tung, Striebel) back over every instant at the end.
type FilterRequest struct {
	T0           float64
	Y0           State
	Covariance   []float64
	TEnd         float64
	Epsilon      float64
	ProcessNoise float64
	Bias         []float64
	Observations []Observation
	Stations     []Station
	Light        *LightTime
	Smooth       bool
}

// FilterStep is an instant the filter kept: the 6 + np + biases unknowns and their covariance
// (the smoothed ones after the smoother).
type FilterStep struct {
	T          float64
	Unknowns   []float64
	Covariance []float64
}

// FilterResult: the innovations (observed - predicted before each correction), the instants
// kept (with Smooth: the start and every observation), the final state, its covariance and
// the biases.
type FilterResult struct {
	Innovations [][]float64
	Steps       []FilterStep
	T           float64
	Y           State
	Covariance  []float64
	Bias        []float64
}

// Filter runs the filter (and the smoother) over the observations in one call. The parameters
// of the model are corrected in place.
func (m *Model) Filter(r FilterRequest) (*FilterResult, error) {
	const op = "Model.Filter"
	observations, values, counts, err := cObservations(op, r.Observations)
	if err != nil {
		return nil, err
	}
	stations := cStations(r.Stations)
	model, err := m.begin()
	if err != nil {
		return nil, err
	}
	defer m.end()
	np := int(C.mybgo_model_parameter_count(model))
	base := 6 + np
	unknowns := base + len(r.Bias)
	if len(r.Covariance) != base*base {
		return nil, refuse(op, -1, "a covariance of (6 + parameters)^2 values is needed")
	}
	options := C.mybgo_filter_options{epsilon: C.double(r.Epsilon),
		process_noise: C.double(r.ProcessNoise), bias_count: C.size_t(len(r.Bias)),
		smooth: flag(r.Smooth)}
	room := len(r.Observations) + 1
	innovations := make([]float64, values)
	stepT := make([]float64, room)
	stepU := make([]float64, room*unknowns)
	stepP := make([]float64, room*unknowns*unknowns)
	out := &FilterResult{Covariance: make([]float64, unknowns*unknowns),
		Bias: slices.Clone(r.Bias)}
	finalBias := make([]float64, len(r.Bias))
	y0 := r.Y0
	var steps, failed C.size_t
	var tFinal C.double
	var final [1]State
	status := C.mybgo_filter_run(model, &options, C.double(r.T0), cStates([]State{y0}),
		cDoubles(r.Covariance), C.double(r.TEnd), cDoubles(r.Bias),
		(*C.mybgo_station)(first(stations)), C.size_t(len(stations)), cLight(r.Light),
		(*C.mybgo_observation)(first(observations)), C.size_t(len(observations)),
		cDoubles(innovations), &steps, cDoubles(stepT), cDoubles(stepU), cDoubles(stepP),
		&tFinal, cStates(final[:]), cDoubles(out.Covariance), cDoubles(finalBias), &failed)
	if err := m.failure(op, status, int(failed)); err != nil {
		return nil, err
	}
	out.T = float64(tFinal)
	out.Y = final[0]
	at := 0
	for _, n := range counts {
		out.Innovations = append(out.Innovations, innovations[at:at+n:at+n])
		at += n
	}
	for k := range int(steps) {
		out.Steps = append(out.Steps, FilterStep{T: stepT[k],
			Unknowns:   stepU[unknowns*k : unknowns*(k+1) : unknowns*(k+1)],
			Covariance: stepP[unknowns*unknowns*k : unknowns*unknowns*(k+1) : unknowns*unknowns*(k+1)]})
	}
	out.Bias = finalBias
	return out, nil
}

// Measure evaluates a measurement at every instant and state, in one call: the values (as many
// as the kind gives, each) and, with partials, their derivatives by the state (6 per value).
// station is needed by the measurements from a station; light (nil: geometric) takes the light
// time into account.
func (m *Model) Measure(kind MeasureKind, station *Station, light *LightTime, t []float64, y []State, partials bool) ([][]float64, [][]float64, error) {
	const op = "Model.Measure"
	n := kind.values()
	if n == 0 {
		return nil, nil, refuse(op, -1, "an unknown measurement")
	}
	if len(t) != len(y) {
		return nil, nil, refuse(op, -1, "times and states of different lengths")
	}
	ckind := C.int(kind)
	if light != nil {
		ckind |= C.MYBGO_MEASURE_LIGHT
	}
	var cs *C.mybgo_station
	if station != nil {
		s := cStation(*station)
		cs = &s
	}
	model, err := m.begin()
	if err != nil {
		return nil, nil, err
	}
	defer m.end()
	values := make([]float64, n*len(t))
	var parts []float64
	if partials {
		parts = make([]float64, 6*n*len(t))
	}
	var failed C.size_t
	status := C.mybgo_measure(model, ckind, cs, cLight(light), cDoubles(t), cStates(y),
		C.size_t(len(t)), cDoubles(values), cDoubles(parts), &failed)
	if err := m.failure(op, status, int(failed)); err != nil {
		return nil, nil, err
	}
	outV := make([][]float64, len(t))
	var outP [][]float64
	for k := range outV {
		outV[k] = values[n*k : n*(k+1) : n*(k+1)]
		if partials {
			outP = append(outP, parts[6*n*k:6*n*(k+1):6*n*(k+1)])
		}
	}
	return outV, outP, nil
}

// OrbitStart is an orbit of a model from Y0 at T0 towards TEnd with the covariance of its
// state then (6 x 6, row-major, m and m/s); the variational equations are integrated.
type OrbitStart struct {
	Options    OrbitOptions
	T0         float64
	Y0         State
	TEnd       float64
	Covariance [36]float64
}

func (s OrbitStart) c() C.mybgo_orbit_start {
	out := C.mybgo_orbit_start{options: s.Options.c(), t0: C.double(s.T0), t_end: C.double(s.TEnd)}
	for i := range 6 {
		out.y0[i] = C.double(s.Y0[i])
	}
	for i := range 36 {
		out.covariance[i] = C.double(s.Covariance[i])
	}
	return out
}

// PropagateCovariance carries the covariance along the orbit: the state and P(t) = Phi P0
// Phi^T at every time (myb_orbit_partials, myb_covariance_transform), in one call.
func (m *Model) PropagateCovariance(start OrbitStart, times []float64) ([]State, [][36]float64, error) {
	const op = "Model.PropagateCovariance"
	model, err := m.begin()
	if err != nil {
		return nil, nil, err
	}
	defer m.end()
	cs := start.c()
	states := make([]State, len(times))
	covariances := make([][36]float64, len(times))
	var failed C.size_t
	status := C.mybgo_covariance_propagate(model, &cs, cDoubles(times), C.size_t(len(times)),
		cStates(states), (*C.double)(first(covariances)), &failed)
	if err := m.failure(op, status, int(failed)); err != nil {
		return nil, nil, err
	}
	return states, covariances, nil
}

// CovariancesRTN turns state covariances into the radial, transverse, normal axes of their
// states (myb_covariance_rtn), in one call.
func CovariancesRTN(states []State, covariances [][36]float64) ([][36]float64, error) {
	const op = "CovariancesRTN"
	if len(states) != len(covariances) {
		return nil, refuse(op, -1, "states and covariances of different lengths")
	}
	out := make([][36]float64, len(states))
	var failed C.size_t
	status := C.mybgo_covariance_rtn(cStates(states), (*C.double)(first(covariances)),
		C.size_t(len(states)), (*C.double)(first(out)), &failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	return out, nil
}
