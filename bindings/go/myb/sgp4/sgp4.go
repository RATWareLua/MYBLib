// SPDX-License-Identifier: BSD-3-Clause

// Package sgp4 is the Go face of ext/sgp4 of myb: two-line element sets read and written, the
// SGP4/SDP4 theory over slices of times, element sets fitted to positions, and their exchange
// with the mean-elements messages of the CCSDS. Plain cgo over the static libraries
// libmyb_sgp4.a and libmyb.a (the pkg-config module myb_sgp4-static), the coarse layer as in
// package myb: one call a task, nothing computed in Go.
package sgp4

/*
#cgo pkg-config: myb_sgp4-static
#include <stdlib.h>
#include <string.h>
#include "myb_sgp4.h"

// The theory at count times of a record: states 6 each (m, m/s), the codes of the theory;
// *failed: the time that failed (count when none did).
static myb_status mybgo_sgp4_many(const myb_sgp4* record, const double* minutes, size_t count,
                                  double* states, int* codes, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    myb_sgp4_error code = MYB_SGP4_OK;
    status = myb_sgp4_propagate(record, minutes[k], states + 6 * k, states + 6 * k + 3, &code);
    codes[k] = (int)code;
    if (status != MYB_OK) {
      break;
    }
  }
  *failed = k;
  return status;
}

// The two lines of a set into one buffer of 140 characters (two lines of 70 with their NULs).
static myb_status mybgo_tle_format(const myb_tle* tle, char* lines) {
  return myb_tle_format(tle, lines, lines + 70);
}
*/
import "C"

import (
	"unsafe"

	"github.com/RATWareLua/MYBLib/bindings/go/myb"
)

// TLE is a two-line element set (myb_tle), in the units of the catalogue.
type TLE struct {
	Satnum             string // as written (5 characters)
	CatalogNumber      int64
	Classification     byte
	Designator         string
	EpochYear          int
	EpochDay           float64
	NDot, NDDot, BStar float64
	EphemerisType      int
	ElementNumber      int
	Inclination        float64 // deg
	Node               float64
	Eccentricity       float64
	ArgumentOfPerigee  float64
	MeanAnomaly        float64
	MeanMotion         float64 // rev/day
	Revolution         int64
}

func fail(op string, status C.myb_status) error {
	if status == C.MYB_OK {
		return nil
	}
	return &myb.Error{Status: myb.Status(status), Op: op, Index: -1}
}

func setChars(dst []C.char, s string) {
	for i := range dst {
		dst[i] = 0
	}
	for i := 0; i < len(s) && i < len(dst)-1; i++ {
		dst[i] = C.char(s[i])
	}
}

func (t TLE) c() C.myb_tle {
	var c C.myb_tle
	setChars(c.satnum[:], t.Satnum)
	setChars(c.designator[:], t.Designator)
	c.catalog_number, c.classification = C.long(t.CatalogNumber), C.char(t.Classification)
	c.epoch_year, c.epoch_day = C.int(t.EpochYear), C.double(t.EpochDay)
	c.ndot, c.nddot, c.bstar = C.double(t.NDot), C.double(t.NDDot), C.double(t.BStar)
	c.ephemeris_type, c.element_number = C.int(t.EphemerisType), C.int(t.ElementNumber)
	c.inclination, c.node, c.eccentricity = C.double(t.Inclination), C.double(t.Node), C.double(t.Eccentricity)
	c.argument_of_perigee, c.mean_anomaly = C.double(t.ArgumentOfPerigee), C.double(t.MeanAnomaly)
	c.mean_motion, c.revolution = C.double(t.MeanMotion), C.long(t.Revolution)
	return c
}

func goTLE(c *C.myb_tle) TLE {
	return TLE{Satnum: C.GoString(&c.satnum[0]), CatalogNumber: int64(c.catalog_number),
		Classification: byte(c.classification), Designator: C.GoString(&c.designator[0]),
		EpochYear: int(c.epoch_year), EpochDay: float64(c.epoch_day), NDot: float64(c.ndot),
		NDDot: float64(c.nddot), BStar: float64(c.bstar), EphemerisType: int(c.ephemeris_type),
		ElementNumber: int(c.element_number), Inclination: float64(c.inclination),
		Node: float64(c.node), Eccentricity: float64(c.eccentricity),
		ArgumentOfPerigee: float64(c.argument_of_perigee), MeanAnomaly: float64(c.mean_anomaly),
		MeanMotion: float64(c.mean_motion), Revolution: int64(c.revolution)}
}

// Parse reads the two lines of a set (myb_tle_parse); noChecksum accepts lines whose checksum
// does not match.
func Parse(line1, line2 string, noChecksum bool) (TLE, error) {
	c1, c2 := C.CString(line1), C.CString(line2)
	defer C.free(unsafe.Pointer(c1))
	defer C.free(unsafe.Pointer(c2))
	var flags C.uint
	if noChecksum {
		flags = C.MYB_TLE_NO_CHECKSUM
	}
	var c C.myb_tle
	var diag C.myb_diag
	status := C.myb_tle_parse(c1, c2, flags, &c, &diag)
	if status != C.MYB_OK {
		return TLE{}, &myb.Error{Status: myb.Status(status), Op: "sgp4.Parse", Index: -1,
			Line: int(diag.line), Message: C.GoString(&diag.message[0])}
	}
	return goTLE(&c), nil
}

// Format is the two lines of the set with their checksums (myb_tle_format).
func (t TLE) Format() (string, string, error) {
	c := t.c()
	var lines [140]C.char
	if err := fail("sgp4.Format", C.mybgo_tle_format(&c, &lines[0])); err != nil {
		return "", "", err
	}
	return C.GoString(&lines[0]), C.GoString(&lines[70]), nil
}

// Epoch is the epoch of the set as a UTC calendar date (myb_tle_epoch).
func (t TLE) Epoch() myb.Calendar {
	c := t.c()
	var cal C.myb_calendar
	C.myb_tle_epoch(&c, &cal)
	return *(*myb.Calendar)(unsafe.Pointer(&cal))
}

// FromOMM is the set of a mean-elements message (myb_tle_from_omm).
func FromOMM(o *myb.OMM) (TLE, error) {
	p, err := myb.OMMToC(o)
	if err != nil {
		return TLE{}, err
	}
	defer myb.FreeC(p)
	var c C.myb_tle
	if err := fail("sgp4.FromOMM", C.myb_tle_from_omm((*C.myb_omm)(p), &c)); err != nil {
		return TLE{}, err
	}
	return goTLE(&c), nil
}

// OMM is the mean-elements message of the set (myb_omm_from_tle); the originator, the name
// and the creation date are the caller's to fill.
func (t TLE) OMM() (*myb.OMM, error) {
	c := t.c()
	p := C.calloc(1, C.sizeof_myb_omm)
	defer C.free(p)
	if err := fail("sgp4.OMM", C.myb_omm_from_tle(&c, (*C.myb_omm)(p))); err != nil {
		return nil, err
	}
	return myb.OMMFromC(unsafe.Pointer(p)), nil
}

// Gravity are the constants of the theory (myb_sgp4_gravity).
type Gravity int

// The standard constants.
const (
	WGS72Old Gravity = C.MYB_SGP4_WGS72_OLD
	WGS72    Gravity = C.MYB_SGP4_WGS72
	WGS84    Gravity = C.MYB_SGP4_WGS84
)

// Mode is the variant of the theory (myb_sgp4_mode).
type Mode int

// The variants.
const (
	Improved Mode = C.MYB_SGP4_IMPROVED
	AFSPC    Mode = C.MYB_SGP4_AFSPC
)

// Record is an initialised record of the theory (myb_sgp4): a value of numbers only.
type Record struct {
	c C.myb_sgp4
}

// NewRecord initialises the theory for a set (myb_sgp4_init).
func NewRecord(t TLE, gravity Gravity, mode Mode) (*Record, error) {
	c := t.c()
	r := &Record{}
	if err := fail("sgp4.NewRecord", C.myb_sgp4_init(&c, C.myb_sgp4_gravity(gravity),
		C.myb_sgp4_mode(mode), &r.c)); err != nil {
		return nil, err
	}
	return r, nil
}

// Propagate is the TEME state (m, m/s) at the minutes since the epoch, and the code of the
// theory (myb_sgp4_propagate: codes 1 .. 4 are errors, 6 a decayed orbit with its state).
func (r *Record) Propagate(minutes float64) (myb.State, int, error) {
	var s myb.State
	var code C.myb_sgp4_error
	status := C.myb_sgp4_propagate(&r.c, C.double(minutes), (*C.double)(unsafe.Pointer(&s[0])),
		(*C.double)(unsafe.Pointer(&s[3])), &code)
	return s, int(code), fail("sgp4.Propagate", status)
}

// PropagateMany is Propagate at every time of a slice, in one call: the states and the codes.
func (r *Record) PropagateMany(minutes []float64) ([]myb.State, []int, error) {
	states := make([]myb.State, len(minutes))
	codes := make([]C.int, len(minutes))
	var failed C.size_t
	var pm *C.double
	var ps *C.double
	var pc *C.int
	if len(minutes) > 0 {
		pm, ps, pc = (*C.double)(unsafe.Pointer(&minutes[0])), (*C.double)(unsafe.Pointer(&states[0])), &codes[0]
	}
	status := C.mybgo_sgp4_many(&r.c, pm, C.size_t(len(minutes)), ps, pc, &failed)
	out := make([]int, len(codes))
	for k, c := range codes {
		out[k] = int(c)
	}
	if status != C.MYB_OK {
		return nil, out, &myb.Error{Status: myb.Status(status), Op: "sgp4.PropagateMany", Index: int(failed)}
	}
	return states, out, nil
}

// FitOptions are the options of a fit (myb_tle_fit_options).
type FitOptions struct {
	BStar         bool
	MaxIterations int
	Tolerance     float64
	Gravity       Gravity
	Mode          Mode
}

// Fit fits the elements of the guess to TEME positions (m) at the minutes from its epoch
// (myb_tle_fit): the set, the root mean square of the residuals and the steps taken. On
// MYB_E_NOCONV the best set so far comes back with the error.
func Fit(guess TLE, minutes []float64, positions [][3]float64, options FitOptions) (TLE, float64, int, error) {
	if len(minutes) != len(positions) || len(minutes) == 0 {
		return guess, 0, 0, &myb.Error{Status: myb.Domain, Op: "sgp4.Fit", Index: -1,
			Message: "times and positions of different lengths"}
	}
	c := guess.c()
	co := C.myb_tle_fit_options{bstar: boolInt(options.BStar), max_iterations: C.int(options.MaxIterations),
		tolerance: C.double(options.Tolerance), gravity: C.myb_sgp4_gravity(options.Gravity),
		mode: C.myb_sgp4_mode(options.Mode)}
	var rms C.double
	var iterations C.int
	status := C.myb_tle_fit((*C.double)(unsafe.Pointer(&minutes[0])), (*C.double)(unsafe.Pointer(&positions[0])),
		C.size_t(len(minutes)), &co, &c, &rms, &iterations)
	return goTLE(&c), float64(rms), int(iterations), fail("sgp4.Fit", status)
}

// FromState writes the osculating elements of a TEME state as the mean ones of the set, for
// the constants of gravity (a first guess for Fit; myb_tle_from_state).
func FromState(t TLE, s myb.State, gravity Gravity) (TLE, error) {
	c := t.c()
	status := C.myb_tle_from_state((*C.double)(unsafe.Pointer(&s[0])), (*C.double)(unsafe.Pointer(&s[3])),
		C.myb_sgp4_gravity(gravity), &c)
	return goTLE(&c), fail("sgp4.FromState", status)
}

func boolInt(b bool) C.int {
	if b {
		return 1
	}
	return 0
}
