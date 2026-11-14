// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include <stdlib.h>
#include "glue.h"
*/
import "C"

import "unsafe"

// text copies a string into a text field of a message (MYB_ODM_TEXT with the NUL).
func setText(op string, dst *[C.MYB_ODM_TEXT]C.char, s string) error {
	if len(s) >= C.MYB_ODM_TEXT {
		return refuse(op, -1, "a text of 64 characters or more: "+s)
	}
	for i := range dst {
		dst[i] = 0
	}
	for i := range len(s) {
		dst[i] = C.char(s[i])
	}
	return nil
}

func getText(src *[C.MYB_ODM_TEXT]C.char) string {
	return C.GoString(&src[0])
}

func cCalendar(c Calendar) C.myb_calendar {
	return *(*C.myb_calendar)(unsafe.Pointer(&c))
}

func goCalendar(c C.myb_calendar) Calendar {
	return *(*Calendar)(unsafe.Pointer(&c))
}

func writeMessage(op string, status C.myb_status, text *C.char, length C.size_t) (string, error) {
	if err := fail(op, status); err != nil {
		return "", err
	}
	defer C.mybgo_free(unsafe.Pointer(text))
	return C.GoStringN(text, C.int(length)), nil
}

// ---- OPM ------------------------------------------------------------------------------------------

// OPMElements are the osculating elements of an OPM (m, rad, m^3/s^2).
type OPMElements struct {
	SemiMajorAxis, Eccentricity, Inclination, RAAN, ArgPericenter, Anomaly float64
	MeanAnomaly                                                            bool
	GM                                                                     float64
}

// OPMSpacecraft are the parameters of the spacecraft (kg, m^2).
type OPMSpacecraft struct {
	Mass, SolarRadArea, SolarRadCoeff, DragArea, DragCoeff float64
}

// OPMManeuver is a manoeuvre of an OPM.
type OPMManeuver struct {
	EpochIgnition Calendar
	Duration      float64
	DeltaMass     float64
	RefFrame      string
	DV            [3]float64
}

// OPM is an orbit parameter message (myb_opm): the state at one epoch (m, m/s) and what may
// come with it.
type OPM struct {
	Originator, ObjectName, ObjectID, CenterName, RefFrame, TimeSystem string
	CreationDate, Epoch                                                Calendar
	State                                                              State
	Elements                                                           *OPMElements
	Spacecraft                                                         *OPMSpacecraft
	Maneuvers                                                          []OPMManeuver
	Covariance                                                         *[36]float64
	CovRefFrame                                                        string
}

func (o *OPM) c(op string) (*C.myb_opm, error) {
	c := (*C.myb_opm)(C.calloc(1, C.sizeof_myb_opm))
	for _, t := range []struct {
		dst *[C.MYB_ODM_TEXT]C.char
		s   string
	}{{&c.originator, o.Originator}, {&c.object_name, o.ObjectName}, {&c.object_id, o.ObjectID},
		{&c.center_name, o.CenterName}, {&c.ref_frame, o.RefFrame}, {&c.time_system, o.TimeSystem},
		{&c.cov_ref_frame, o.CovRefFrame}} {
		if err := setText(op, t.dst, t.s); err != nil {
			C.free(unsafe.Pointer(c))
			return nil, err
		}
	}
	c.creation_date, c.epoch = cCalendar(o.CreationDate), cCalendar(o.Epoch)
	for i := range 6 {
		c.state[i] = C.double(o.State[i])
	}
	if e := o.Elements; e != nil {
		c.has_elements = 1
		c.semi_major_axis, c.eccentricity = C.double(e.SemiMajorAxis), C.double(e.Eccentricity)
		c.inclination, c.ra_of_asc_node = C.double(e.Inclination), C.double(e.RAAN)
		c.arg_of_pericenter, c.anomaly = C.double(e.ArgPericenter), C.double(e.Anomaly)
		c.anomaly_is_mean, c.gm = flag(e.MeanAnomaly), C.double(e.GM)
	}
	if s := o.Spacecraft; s != nil {
		c.has_spacecraft = 1
		c.mass, c.solar_rad_area, c.solar_rad_coeff = C.double(s.Mass), C.double(s.SolarRadArea), C.double(s.SolarRadCoeff)
		c.drag_area, c.drag_coeff = C.double(s.DragArea), C.double(s.DragCoeff)
	}
	if len(o.Maneuvers) > C.MYB_OPM_MAX_MANEUVERS {
		C.free(unsafe.Pointer(c))
		return nil, refuse(op, -1, "more than 16 manoeuvres")
	}
	c.maneuver_count = C.size_t(len(o.Maneuvers))
	for i, m := range o.Maneuvers {
		d := &c.maneuvers[i]
		d.epoch_ignition, d.duration, d.delta_mass = cCalendar(m.EpochIgnition), C.double(m.Duration), C.double(m.DeltaMass)
		d.dv = cVec3(m.DV)
		if err := setText(op, &d.ref_frame, m.RefFrame); err != nil {
			C.free(unsafe.Pointer(c))
			return nil, err
		}
	}
	if o.Covariance != nil {
		c.has_covariance = 1
		for i := range 36 {
			c.covariance[i] = C.double(o.Covariance[i])
		}
	}
	return c, nil
}

func goOPM(c *C.myb_opm) *OPM {
	o := &OPM{Originator: getText(&c.originator), ObjectName: getText(&c.object_name),
		ObjectID: getText(&c.object_id), CenterName: getText(&c.center_name),
		RefFrame: getText(&c.ref_frame), TimeSystem: getText(&c.time_system),
		CovRefFrame: getText(&c.cov_ref_frame), CreationDate: goCalendar(c.creation_date),
		Epoch: goCalendar(c.epoch)}
	for i := range 6 {
		o.State[i] = float64(c.state[i])
	}
	if c.has_elements != 0 {
		o.Elements = &OPMElements{SemiMajorAxis: float64(c.semi_major_axis),
			Eccentricity: float64(c.eccentricity), Inclination: float64(c.inclination),
			RAAN: float64(c.ra_of_asc_node), ArgPericenter: float64(c.arg_of_pericenter),
			Anomaly: float64(c.anomaly), MeanAnomaly: c.anomaly_is_mean != 0, GM: float64(c.gm)}
	}
	if c.has_spacecraft != 0 {
		o.Spacecraft = &OPMSpacecraft{Mass: float64(c.mass), SolarRadArea: float64(c.solar_rad_area),
			SolarRadCoeff: float64(c.solar_rad_coeff), DragArea: float64(c.drag_area),
			DragCoeff: float64(c.drag_coeff)}
	}
	for i := range int(c.maneuver_count) {
		m := &c.maneuvers[i]
		o.Maneuvers = append(o.Maneuvers, OPMManeuver{EpochIgnition: goCalendar(m.epoch_ignition),
			Duration: float64(m.duration), DeltaMass: float64(m.delta_mass),
			RefFrame: getText(&m.ref_frame), DV: goVec3(m.dv)})
	}
	if c.has_covariance != 0 {
		var p [36]float64
		for i := range 36 {
			p[i] = float64(c.covariance[i])
		}
		o.Covariance = &p
	}
	return o
}

// ParseOPM reads an orbit parameter message in either notation (myb_opm_parse).
func ParseOPM(text []byte) (*OPM, error) {
	ctext := C.CBytes(text)
	defer C.free(ctext)
	c := (*C.myb_opm)(C.calloc(1, C.sizeof_myb_opm))
	defer C.free(unsafe.Pointer(c))
	var diag C.myb_diag
	if err := failDiag("ParseOPM", C.myb_opm_parse((*C.char)(ctext), C.size_t(len(text)), c, &diag),
		&diag); err != nil {
		return nil, err
	}
	return goOPM(c), nil
}

// WriteOPM is the text of the message in the notation given (myb_opm_write_as).
func WriteOPM(o *OPM, notation Notation) (string, error) {
	const op = "WriteOPM"
	if o == nil {
		return "", refuse(op, -1, "no message")
	}
	c, err := o.c(op)
	if err != nil {
		return "", err
	}
	defer C.free(unsafe.Pointer(c))
	var text *C.char
	var length C.size_t
	status := C.mybgo_opm_write(c, C.int(notation), &text, &length)
	return writeMessage(op, status, text, length)
}

// ---- OMM ------------------------------------------------------------------------------------------

// OMMTLE are the parameters of an OMM that belong to a TLE.
type OMMTLE struct {
	EphemerisType                        int
	ClassificationType                   byte
	NoradCatID, ElementSetNo, RevAtEpoch int64
	BStar, MeanMotionDot, MeanMotionDDot float64
}

// OMM is a mean-elements message (myb_omm), its numbers in the units of the message (rev/day,
// degrees, km, km^3/s^2).
type OMM struct {
	Originator, ObjectName, ObjectID, CenterName, RefFrame, TimeSystem string
	MeanElementTheory                                                  string
	CreationDate, Epoch                                                Calendar
	HasMeanMotion                                                      bool // MEAN_MOTION, else SEMI_MAJOR_AXIS
	MeanMotion, SemiMajorAxis                                          float64
	Eccentricity, Inclination, RAAN, ArgPericenter, MeanAnomaly        float64
	GM                                                                 *float64
	TLE                                                                *OMMTLE
	Covariance                                                         *[36]float64
	CovRefFrame                                                        string
}

func (o *OMM) c(op string) (*C.myb_omm, error) {
	c := (*C.myb_omm)(C.calloc(1, C.sizeof_myb_omm))
	for _, t := range []struct {
		dst *[C.MYB_ODM_TEXT]C.char
		s   string
	}{{&c.originator, o.Originator}, {&c.object_name, o.ObjectName}, {&c.object_id, o.ObjectID},
		{&c.center_name, o.CenterName}, {&c.ref_frame, o.RefFrame}, {&c.time_system, o.TimeSystem},
		{&c.mean_element_theory, o.MeanElementTheory}, {&c.cov_ref_frame, o.CovRefFrame}} {
		if err := setText(op, t.dst, t.s); err != nil {
			C.free(unsafe.Pointer(c))
			return nil, err
		}
	}
	c.creation_date, c.epoch = cCalendar(o.CreationDate), cCalendar(o.Epoch)
	c.has_mean_motion = flag(o.HasMeanMotion)
	c.mean_motion, c.semi_major_axis = C.double(o.MeanMotion), C.double(o.SemiMajorAxis)
	c.eccentricity, c.inclination = C.double(o.Eccentricity), C.double(o.Inclination)
	c.ra_of_asc_node, c.arg_of_pericenter = C.double(o.RAAN), C.double(o.ArgPericenter)
	c.mean_anomaly = C.double(o.MeanAnomaly)
	if o.GM != nil {
		c.has_gm, c.gm = 1, C.double(*o.GM)
	}
	if t := o.TLE; t != nil {
		c.has_tle = 1
		c.ephemeris_type, c.classification_type = C.int(t.EphemerisType), C.char(t.ClassificationType)
		c.norad_cat_id, c.element_set_no, c.rev_at_epoch = C.long(t.NoradCatID), C.long(t.ElementSetNo), C.long(t.RevAtEpoch)
		c.bstar, c.mean_motion_dot, c.mean_motion_ddot = C.double(t.BStar), C.double(t.MeanMotionDot), C.double(t.MeanMotionDDot)
	}
	if o.Covariance != nil {
		c.has_covariance = 1
		for i := range 36 {
			c.covariance[i] = C.double(o.Covariance[i])
		}
	}
	return c, nil
}

func goOMM(c *C.myb_omm) *OMM {
	o := &OMM{Originator: getText(&c.originator), ObjectName: getText(&c.object_name),
		ObjectID: getText(&c.object_id), CenterName: getText(&c.center_name),
		RefFrame: getText(&c.ref_frame), TimeSystem: getText(&c.time_system),
		MeanElementTheory: getText(&c.mean_element_theory), CovRefFrame: getText(&c.cov_ref_frame),
		CreationDate: goCalendar(c.creation_date), Epoch: goCalendar(c.epoch),
		HasMeanMotion: c.has_mean_motion != 0, MeanMotion: float64(c.mean_motion),
		SemiMajorAxis: float64(c.semi_major_axis), Eccentricity: float64(c.eccentricity),
		Inclination: float64(c.inclination), RAAN: float64(c.ra_of_asc_node),
		ArgPericenter: float64(c.arg_of_pericenter), MeanAnomaly: float64(c.mean_anomaly)}
	if c.has_gm != 0 {
		gm := float64(c.gm)
		o.GM = &gm
	}
	if c.has_tle != 0 {
		o.TLE = &OMMTLE{EphemerisType: int(c.ephemeris_type), ClassificationType: byte(c.classification_type),
			NoradCatID: int64(c.norad_cat_id), ElementSetNo: int64(c.element_set_no),
			RevAtEpoch: int64(c.rev_at_epoch), BStar: float64(c.bstar),
			MeanMotionDot: float64(c.mean_motion_dot), MeanMotionDDot: float64(c.mean_motion_ddot)}
	}
	if c.has_covariance != 0 {
		var p [36]float64
		for i := range 36 {
			p[i] = float64(c.covariance[i])
		}
		o.Covariance = &p
	}
	return o
}

// ParseOMM reads a mean-elements message in either notation (myb_omm_parse).
func ParseOMM(text []byte) (*OMM, error) {
	ctext := C.CBytes(text)
	defer C.free(ctext)
	c := (*C.myb_omm)(C.calloc(1, C.sizeof_myb_omm))
	defer C.free(unsafe.Pointer(c))
	var diag C.myb_diag
	if err := failDiag("ParseOMM", C.myb_omm_parse((*C.char)(ctext), C.size_t(len(text)), c, &diag),
		&diag); err != nil {
		return nil, err
	}
	return goOMM(c), nil
}

// WriteOMM is the text of the message in the notation given (myb_omm_write_as).
func WriteOMM(o *OMM, notation Notation) (string, error) {
	const op = "WriteOMM"
	if o == nil {
		return "", refuse(op, -1, "no message")
	}
	c, err := o.c(op)
	if err != nil {
		return "", err
	}
	defer C.free(unsafe.Pointer(c))
	var text *C.char
	var length C.size_t
	status := C.mybgo_omm_write(c, C.int(notation), &text, &length)
	return writeMessage(op, status, text, length)
}

// OMMToC and OMMFromC let the package of ext/sgp4 pass a message to C: a copy in C memory
// (free it with FreeC), and a message read back from one.
func OMMToC(o *OMM) (unsafe.Pointer, error) {
	c, err := o.c("OMMToC")
	return unsafe.Pointer(c), err
}

// OMMFromC is the message of a myb_omm in C memory.
func OMMFromC(p unsafe.Pointer) *OMM {
	return goOMM((*C.myb_omm)(p))
}

// FreeC releases memory that OMMToC allocated.
func FreeC(p unsafe.Pointer) {
	C.free(p)
}
