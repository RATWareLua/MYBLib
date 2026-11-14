// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include <stdlib.h>
#include "glue.h"
*/
import "C"

import (
	"sync"
	"unsafe"
)

// Ephemeris is a set of SPK files (myb_spk): the states of bodies relative to others as
// functions of the ephemeris time (TDB for the JPL ephemerides), in m and m/s.
type Ephemeris struct {
	o      *object
	images *images
}

// images are the copies in C memory of the files given as bytes: the set reads them until
// it is destroyed.
type images struct {
	mu   sync.Mutex
	list []unsafe.Pointer
}

// Segment describes a segment of an SPK set (myb_spk_segment_info).
type Segment struct {
	Target, Center, Frame, Type int32
	Start, Stop                 float64 // the coverage, seconds past J2000 on the ephemeris' scale
	File                        int     // the order of the file in the set
	Name                        string
}

// NewEphemeris is an empty set.
func NewEphemeris() (*Ephemeris, error) {
	var s *C.myb_spk
	if err := fail("NewEphemeris", C.myb_spk_create(nil, &s)); err != nil {
		return nil, err
	}
	held := &images{}
	return &Ephemeris{images: held, o: newObject(unsafe.Pointer(s), func(q unsafe.Pointer) {
		C.myb_spk_destroy((*C.myb_spk)(q))
		held.mu.Lock()
		for _, p := range held.list {
			C.free(p)
		}
		held.list = nil
		held.mu.Unlock()
	})}, nil
}

func (e *Ephemeris) obj() *object {
	if e == nil {
		return nil
	}
	return e.o
}

// Close releases the set (after the models that use it).
func (e *Ephemeris) Close() error {
	e.obj().close()
	return nil
}

// LoadFile adds an SPK file; later files take precedence.
func (e *Ephemeris) LoadFile(path string) error {
	ptr, err := e.obj().acquire()
	if err != nil {
		return err
	}
	defer e.o.release()
	cpath := C.CString(path)
	defer C.free(unsafe.Pointer(cpath))
	var diag C.myb_diag
	return failDiag("Ephemeris.LoadFile", C.myb_spk_load_file((*C.myb_spk)(ptr), cpath, &diag),
		&diag)
}

// LoadBytes adds an SPK file given as its bytes (copied).
func (e *Ephemeris) LoadBytes(data []byte) error {
	ptr, err := e.obj().acquire()
	if err != nil {
		return err
	}
	defer e.o.release()
	if len(data) == 0 {
		return refuse("Ephemeris.LoadBytes", -1, "no bytes")
	}
	image := C.CBytes(data)
	var diag C.myb_diag
	status := C.myb_spk_load_memory((*C.myb_spk)(ptr), image, C.size_t(len(data)), &diag)
	if status != C.MYB_OK {
		C.free(image)
		return failDiag("Ephemeris.LoadBytes", status, &diag)
	}
	e.images.mu.Lock()
	e.images.list = append(e.images.list, image)
	e.images.mu.Unlock()
	return nil
}

// Segments lists the segments of the set in load order.
func (e *Ephemeris) Segments() ([]Segment, error) {
	ptr, err := e.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer e.o.release()
	spk := (*C.myb_spk)(ptr)
	count := int(C.myb_spk_segment_count(spk))
	out := make([]Segment, count)
	for k := range out {
		var info C.myb_spk_segment_info
		C.myb_spk_segment(spk, C.size_t(k), &info)
		out[k] = Segment{Target: int32(info.target), Center: int32(info.center),
			Frame: int32(info.frame), Type: int32(info._type), Start: float64(info.start),
			Stop: float64(info.stop), File: int(info.file), Name: C.GoString(&info.name[0])}
	}
	return out, nil
}

// State is the state of target relative to center at t, and the frame code of the segments
// (myb_spk_state).
func (e *Ephemeris) State(target, center int32, t Time) (State, int32, error) {
	if !validTime(t) {
		return State{}, 0, refuse("Ephemeris.State", -1, "an instant out of range")
	}
	ptr, err := e.obj().acquire()
	if err != nil {
		return State{}, 0, err
	}
	defer e.o.release()
	var s State
	var frame C.int32_t
	status := C.myb_spk_state((*C.myb_spk)(ptr), C.int32_t(target), C.int32_t(center),
		C.myb_time{s: C.int64_t(t.S), f: C.double(t.F)}, (*C.double)(unsafe.Pointer(&s)), &frame)
	return s, int32(frame), fail("Ephemeris.State", status)
}

// States is State at every instant of a slice, in one call (myb_spk_state_many).
func (e *Ephemeris) States(target, center int32, t []Time) ([]State, []int32, error) {
	if err := checkTimes("Ephemeris.States", t); err != nil {
		return nil, nil, err
	}
	ptr, err := e.obj().acquire()
	if err != nil {
		return nil, nil, err
	}
	defer e.o.release()
	states := make([]State, len(t))
	frames := make([]int32, len(t))
	var failed C.size_t
	status := C.myb_spk_state_many((*C.myb_spk)(ptr), C.int32_t(target), C.int32_t(center),
		cTimes(t), C.size_t(len(t)), cStates(states), (*C.int32_t)(first(frames)), &failed)
	if err := failAt("Ephemeris.States", status, failed); err != nil {
		return nil, nil, err
	}
	return states, frames, nil
}

// ---- Writing SPK files -------------------------------------------------------------------------

// SPKSegment is a segment to write: states at increasing epochs (seconds past J2000, TDB),
// interpolated by Lagrange polynomials (Type 9) or Hermite polynomials of odd degree (Type
// 13) of the given degree (1 .. 27).
type SPKSegment struct {
	Target, Center, Frame int32
	Name                  string // 1 .. 40 printable ASCII characters
	Type, Degree          int
	Epochs                []float64
	States                []State
}

// spkWrite lays the segments end to end for the writer: to path, or into bytes.
func spkWrite(op, path, internalName string, segments []SPKSegment) ([]byte, error) {
	headers := make([]C.myb_spk_segment_header, len(segments))
	types := make([]C.int, len(segments))
	degrees := make([]C.int, len(segments))
	counts := make([]C.size_t, len(segments))
	var epochs []float64
	var states []State
	names := make([]*C.char, len(segments))
	defer func() {
		for _, n := range names {
			C.free(unsafe.Pointer(n))
		}
	}()
	for j, s := range segments {
		if len(s.Epochs) != len(s.States) {
			return nil, refuse(op, j, "epochs and states of different lengths")
		}
		names[j] = C.CString(s.Name)
		headers[j] = C.myb_spk_segment_header{target: C.int32_t(s.Target),
			center: C.int32_t(s.Center), frame: C.int32_t(s.Frame), name: names[j]}
		types[j] = C.int(s.Type)
		degrees[j] = C.int(s.Degree)
		counts[j] = C.size_t(len(s.Epochs))
		epochs = append(epochs, s.Epochs...)
		states = append(states, s.States...)
	}
	cname := C.CString(internalName)
	defer C.free(unsafe.Pointer(cname))
	var cpath *C.char
	if path != "" {
		cpath = C.CString(path)
		defer C.free(unsafe.Pointer(cpath))
	}
	var data unsafe.Pointer
	var size, failed C.size_t
	status := C.mybgo_spk_write(cpath, cname, C.size_t(len(segments)),
		(*C.myb_spk_segment_header)(first(headers)), (*C.int)(first(types)),
		(*C.int)(first(degrees)), (*C.size_t)(first(counts)), cDoubles(epochs), cStates(states),
		&data, &size, &failed)
	if status != C.MYB_OK {
		e := &Error{Status: Status(status), Op: op, Index: -1}
		if int(failed) < len(segments) {
			e.Index = int(failed) // the segment the writer refused
		}
		return nil, e
	}
	if data == nil {
		return nil, nil
	}
	defer C.mybgo_free(data)
	return C.GoBytes(data, C.int(size)), nil
}

// WriteSPK writes the segments as one SPK file at path (created or replaced); internalName
// labels the file (at most 60 printable ASCII characters).
func WriteSPK(path, internalName string, segments []SPKSegment) error {
	if path == "" {
		return refuse("WriteSPK", -1, "no path")
	}
	_, err := spkWrite("WriteSPK", path, internalName, segments)
	return err
}

// EncodeSPK is the SPK file of the segments as bytes.
func EncodeSPK(internalName string, segments []SPKSegment) ([]byte, error) {
	return spkWrite("EncodeSPK", "", internalName, segments)
}
