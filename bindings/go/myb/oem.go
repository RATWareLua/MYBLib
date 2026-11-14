// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include <stdlib.h>
#include "glue.h"
*/
import "C"

import (
	"os"
	"strconv"
	"unsafe"
)

// Notation is the notation of a CCSDS message.
type Notation int

// The two notations of the orbit data messages.
const (
	KVN Notation = C.MYB_ODM_KVN // KEYWORD = value lines
	XML Notation = C.MYB_ODM_XML
)

// OEMSegment is a segment of an ephemeris message: its metadata and its lines (epochs as
// labels on the segment's time system, states in m and m/s).
type OEMSegment struct {
	ObjectName          string
	ObjectID            string
	CenterName          string
	RefFrame            string
	TimeSystem          string
	Interpolation       string // "": none
	InterpolationDegree int
	Epochs              []Calendar
	States              []State
	Accelerations       [][3]float64      // read: nil when no line gives one
	Meta                map[string]string // read: the metadata keywords of the message (oemKeys)
}

// OEM is an ephemeris message as read.
type OEM struct {
	Version      string
	Originator   string
	CreationDate string
	Segments     []OEMSegment
}

// oemKeys are the metadata keywords a segment is read for.
var oemKeys = []string{"OBJECT_NAME", "OBJECT_ID", "CENTER_NAME", "REF_FRAME", "REF_FRAME_EPOCH",
	"TIME_SYSTEM", "START_TIME", "USEABLE_START_TIME", "USEABLE_STOP_TIME", "STOP_TIME",
	"INTERPOLATION", "INTERPOLATION_DEGREE"}

// WriteOEM is the text of an ephemeris message (version 2.0) of the segments, in the notation
// given (myb_oem_write_as): read back, it gives the states exactly but for the few values the
// conversion between km and m cannot reach, one unit of the last place away.
func WriteOEM(originator string, created Calendar, segments []OEMSegment, notation Notation) (string, error) {
	const op = "WriteOEM"
	k := &keeper{}
	defer k.freeStrings()
	names := make([]*C.char, 0, 6*len(segments))
	degrees := make([]C.int, len(segments))
	counts := make([]C.size_t, len(segments))
	var epochs []Calendar
	var states []State
	for j, s := range segments {
		if len(s.Epochs) != len(s.States) {
			return "", refuse(op, j, "epochs and states of different lengths")
		}
		names = append(names, k.cstring(s.ObjectName), k.cstring(s.ObjectID),
			k.cstring(s.CenterName), k.cstring(s.RefFrame), k.cstring(s.TimeSystem),
			k.optional(s.Interpolation))
		degrees[j] = C.int(s.InterpolationDegree)
		counts[j] = C.size_t(len(s.Epochs))
		epochs = append(epochs, s.Epochs...)
		states = append(states, s.States...)
	}
	var text *C.char
	var length C.size_t
	status := C.mybgo_oem_write(k.cstring(originator), (*C.myb_calendar)(unsafe.Pointer(&created)),
		C.int(notation), C.size_t(len(segments)), (**C.char)(first(names)), (*C.int)(first(degrees)),
		(*C.size_t)(first(counts)), (*C.myb_calendar)(first(epochs)), cStates(states), &text, &length)
	if err := fail(op, status); err != nil {
		return "", err
	}
	defer C.mybgo_free(unsafe.Pointer(text))
	return C.GoStringN(text, C.int(length)), nil
}

// WriteOEMFile writes the message to a file (created or replaced).
func WriteOEMFile(path, originator string, created Calendar, segments []OEMSegment, notation Notation) error {
	text, err := WriteOEM(originator, created, segments, notation)
	if err != nil {
		return err
	}
	if err := os.WriteFile(path, []byte(text), 0o644); err != nil {
		return &Error{Status: IO, Op: "WriteOEMFile", Index: -1, Err: err}
	}
	return nil
}

// ParseOEM reads an ephemeris message in either notation (myb_oem_load).
func ParseOEM(text []byte) (*OEM, error) {
	ctext := C.CBytes(text)
	defer C.free(ctext)
	var oem *C.myb_oem
	var diag C.myb_diag
	status := C.myb_oem_load((*C.char)(ctext), C.size_t(len(text)), nil, &oem, &diag)
	if err := failDiag("ParseOEM", status, &diag); err != nil {
		return nil, err
	}
	defer C.myb_oem_destroy(oem)
	return readOEM(oem), nil
}

// ReadOEMFile reads an ephemeris message from a file (myb_oem_load_file).
func ReadOEMFile(path string) (*OEM, error) {
	cpath := C.CString(path)
	defer C.free(unsafe.Pointer(cpath))
	var oem *C.myb_oem
	var diag C.myb_diag
	status := C.myb_oem_load_file(cpath, nil, &oem, &diag)
	if err := failDiag("ReadOEMFile", status, &diag); err != nil {
		return nil, err
	}
	defer C.myb_oem_destroy(oem)
	return readOEM(oem), nil
}

func header(oem *C.myb_oem, key string) string {
	ckey := C.CString(key)
	defer C.free(unsafe.Pointer(ckey))
	if v := C.myb_oem_header(oem, ckey); v != nil {
		return C.GoString(v)
	}
	return ""
}

// readOEM copies a message into Go values: a call per segment for its lines and one for its
// metadata.
func readOEM(oem *C.myb_oem) *OEM {
	out := &OEM{Version: header(oem, "CCSDS_OEM_VERS"), Originator: header(oem, "ORIGINATOR"),
		CreationDate: header(oem, "CREATION_DATE")}
	ckeys := make([]*C.char, len(oemKeys)) // C strings: no Go pointers for C to see
	for i, key := range oemKeys {
		ckeys[i] = C.CString(key)
	}
	defer func() {
		for _, key := range ckeys {
			C.free(unsafe.Pointer(key))
		}
	}()
	values := make([]*C.char, len(oemKeys))
	count := int(C.myb_oem_segment_count(oem))
	out.Segments = make([]OEMSegment, count)
	for j := range out.Segments {
		s := &out.Segments[j]
		C.mybgo_oem_meta(oem, C.size_t(j), (**C.char)(first(ckeys)), C.size_t(len(oemKeys)),
			(**C.char)(first(values)))
		s.Meta = make(map[string]string)
		for i, v := range values {
			if v != nil {
				s.Meta[oemKeys[i]] = C.GoString(v)
			}
		}
		s.ObjectName = s.Meta["OBJECT_NAME"]
		s.ObjectID = s.Meta["OBJECT_ID"]
		s.CenterName = s.Meta["CENTER_NAME"]
		s.RefFrame = s.Meta["REF_FRAME"]
		s.TimeSystem = s.Meta["TIME_SYSTEM"]
		s.Interpolation = s.Meta["INTERPOLATION"]
		if degree, err := strconv.Atoi(s.Meta["INTERPOLATION_DEGREE"]); err == nil {
			s.InterpolationDegree = degree
		}
		n := int(C.myb_oem_state_count(oem, C.size_t(j)))
		s.Epochs = make([]Calendar, n)
		s.States = make([]State, n)
		accelerations := make([][3]float64, n)
		with := C.mybgo_oem_lines(oem, C.size_t(j), (*C.myb_calendar)(first(s.Epochs)),
			cStates(s.States), (*C.double)(first(accelerations)))
		if with > 0 {
			s.Accelerations = accelerations
		}
	}
	return out
}
