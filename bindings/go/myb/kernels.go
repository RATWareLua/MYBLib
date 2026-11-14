// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include <stdlib.h>
#include "glue.h"
*/
import "C"

import "unsafe"

// ---- SP3 ------------------------------------------------------------------------------------------

// SP3 is a file of precise orbits (myb_sp3_data): records at the epochs (index epoch *
// satellites + satellite), positions in m (NaN: missing), the optional parts nil when absent.
type SP3 struct {
	Version    byte // 'c' or 'd'
	FileType   byte
	TimeSystem string
	Frame      string
	OrbitType  string
	Agency     string
	DataUsed   string
	Satellites []string
	Epochs     []Calendar
	Interval   float64
	Positions  [][3]float64
	Velocities [][3]float64
	Clocks     []float64 // s; NaN: none
	ClockRates []float64
	Accuracy   []float64 // m a satellite; 0: unknown
	Sigmas     [][8]float64
	Flags      []uint8
}

// WriteSP3 is the text of the file (myb_sp3_write).
func WriteSP3(s *SP3) (string, error) {
	const op = "WriteSP3"
	if s == nil {
		return "", refuse(op, -1, "no file")
	}
	records := len(s.Epochs) * len(s.Satellites)
	if len(s.Positions) != records || (s.Velocities != nil && len(s.Velocities) != records) ||
		(s.Clocks != nil && len(s.Clocks) != records) || (s.ClockRates != nil && len(s.ClockRates) != records) ||
		(s.Accuracy != nil && len(s.Accuracy) != len(s.Satellites)) ||
		(s.Sigmas != nil && len(s.Sigmas) != records) || (s.Flags != nil && len(s.Flags) != records) {
		return "", refuse(op, -1, "parts of other lengths than epochs x satellites")
	}
	k := newKeeper()
	defer k.done()
	head := C.myb_sp3_data{version: C.char(s.Version), file_type: C.char(s.FileType),
		time_system: k.cstring(s.TimeSystem), frame: k.cstring(s.Frame), orbit_type: k.cstring(s.OrbitType),
		agency: k.cstring(s.Agency), data_used: k.cstring(s.DataUsed),
		satellite_count: C.size_t(len(s.Satellites)), epoch_count: C.size_t(len(s.Epochs)),
		interval: C.double(s.Interval)}
	names := make([]*C.char, len(s.Satellites))
	for i, n := range s.Satellites {
		names[i] = k.cstring(n)
	}
	var text *C.char
	var length C.size_t
	status := C.mybgo_sp3_write(&head, (**C.char)(first(names)), (*C.myb_calendar)(first(s.Epochs)),
		(*C.double)(first(s.Positions)), (*C.double)(first(s.Velocities)), cDoubles(s.Clocks),
		cDoubles(s.ClockRates), cDoubles(s.Accuracy), (*C.double)(first(s.Sigmas)),
		(*C.uchar)(first(s.Flags)), &text, &length)
	return writeMessage(op, status, text, length)
}

// ParseSP3 reads an SP3 file (myb_sp3_load): the satellites, epochs, time system, frame,
// accuracies and every record (the version, type, orbit type, agency and data used are not
// read back).
func ParseSP3(text []byte) (*SP3, error) {
	ctext := C.CBytes(text)
	defer C.free(ctext)
	var sp3 *C.myb_sp3
	var diag C.myb_diag
	status := C.myb_sp3_load((*C.char)(ctext), C.size_t(len(text)), nil, &sp3, &diag)
	if err := failDiag("ParseSP3", status, &diag); err != nil {
		return nil, err
	}
	defer C.myb_sp3_destroy(sp3)
	ns, ne := int(C.myb_sp3_satellite_count(sp3)), int(C.myb_sp3_epoch_count(sp3))
	out := &SP3{TimeSystem: C.GoString(C.myb_sp3_time_system(sp3)), Frame: C.GoString(C.myb_sp3_frame(sp3))}
	for i := range ns {
		out.Satellites = append(out.Satellites, C.GoString(C.myb_sp3_satellite(sp3, C.size_t(i))))
		out.Accuracy = append(out.Accuracy, float64(C.myb_sp3_accuracy(sp3, C.size_t(i))))
	}
	out.Epochs = make([]Calendar, ne)
	for e := range ne {
		C.myb_sp3_epoch(sp3, C.size_t(e), (*C.myb_calendar)(unsafe.Pointer(&out.Epochs[e])))
	}
	n := ns * ne
	out.Positions, out.Velocities = make([][3]float64, n), make([][3]float64, n)
	out.Clocks, out.ClockRates = make([]float64, n), make([]float64, n)
	out.Sigmas, out.Flags = make([][8]float64, n), make([]uint8, n)
	given := make([]C.int, n)
	C.mybgo_sp3_records(sp3, (*C.double)(first(out.Positions)), (*C.double)(first(out.Velocities)),
		(*C.int)(first(given)), cDoubles(out.Clocks), cDoubles(out.ClockRates),
		(*C.double)(first(out.Sigmas)), (*C.uchar)(first(out.Flags)))
	return out, nil
}

// ---- Clocks of spacecraft --------------------------------------------------------------------------

// SCLK is the clock of a spacecraft (myb_sclk) from the variables of an SCLK kernel.
type SCLK struct {
	o *object
}

// NewSCLK reads the clock of spacecraft id from a pool (myb_sclk_create).
func NewSCLK(pool *KernelPool, id int64) (*SCLK, error) {
	pp, err := pool.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer pool.o.release()
	var s *C.myb_sclk
	var diag C.myb_diag
	if err := failDiag("NewSCLK", C.myb_sclk_create((*C.myb_kpool)(pp), C.int64_t(id), nil, &s, &diag),
		&diag); err != nil {
		return nil, err
	}
	return &SCLK{o: newObject(unsafe.Pointer(s), func(q unsafe.Pointer) {
		C.myb_sclk_destroy((*C.myb_sclk)(q))
	})}, nil
}

func (s *SCLK) obj() *object {
	if s == nil {
		return nil
	}
	return s.o
}

// Close releases the clock (after what uses it).
func (s *SCLK) Close() error {
	s.obj().close()
	return nil
}

// TimeSystem is the scale of the clock's table (1 TDB, 2 TT); Partitions their count.
func (s *SCLK) TimeSystem() (int, int, error) {
	ptr, err := s.obj().acquire()
	if err != nil {
		return 0, 0, err
	}
	defer s.o.release()
	return int(C.myb_sclk_time_system((*C.myb_sclk)(ptr))),
		int(C.myb_sclk_partitions((*C.myb_sclk)(ptr))), nil
}

// Times are the instants of encoded ticks (myb_sclk_to_time), in one call.
func (s *SCLK) Times(ticks []float64) ([]Time, error) {
	ptr, err := s.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer s.o.release()
	out := make([]Time, len(ticks))
	var failed C.size_t
	status := C.mybgo_sclk_times((*C.myb_sclk)(ptr), cDoubles(ticks), C.size_t(len(ticks)), cTimes(out), &failed)
	if err := failAt("SCLK.Times", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// Ticks are the encoded ticks of instants (myb_sclk_from_time), in one call.
func (s *SCLK) Ticks(t []Time) ([]float64, error) {
	if err := checkTimes("SCLK.Ticks", t); err != nil {
		return nil, err
	}
	ptr, err := s.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer s.o.release()
	out := make([]float64, len(t))
	var failed C.size_t
	status := C.mybgo_sclk_ticks((*C.myb_sclk)(ptr), cTimes(t), C.size_t(len(t)), cDoubles(out), &failed)
	if err := failAt("SCLK.Ticks", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// Encode is the encoded ticks of a clock string (myb_sclk_encode).
func (s *SCLK) Encode(text string) (float64, error) {
	ptr, err := s.obj().acquire()
	if err != nil {
		return 0, err
	}
	defer s.o.release()
	ctext := C.CString(text)
	defer C.free(unsafe.Pointer(ctext))
	var ticks C.double
	status := C.myb_sclk_encode((*C.myb_sclk)(ptr), ctext, C.size_t(len(text)), &ticks)
	return float64(ticks), fail("SCLK.Encode", status)
}

// Decode is the clock string of encoded ticks (myb_sclk_decode).
func (s *SCLK) Decode(ticks float64) (string, error) {
	ptr, err := s.obj().acquire()
	if err != nil {
		return "", err
	}
	defer s.o.release()
	var buffer [128]C.char
	var length C.size_t
	status := C.myb_sclk_decode((*C.myb_sclk)(ptr), C.double(ticks), &buffer[0], C.size_t(len(buffer)), &length)
	if err := fail("SCLK.Decode", status); err != nil {
		return "", err
	}
	return C.GoStringN(&buffer[0], C.int(length)), nil
}

// ---- CK --------------------------------------------------------------------------------------------

// CK is a set of CK files (myb_ck): the pointing of instruments in encoded ticks of a clock.
type CK struct {
	o      *object
	images *images
}

// CKSegmentInfo describes a segment of a CK set (myb_ck_segment_info).
type CKSegmentInfo struct {
	Instrument, Frame, Type int32
	HasRates                bool
	Start, Stop             float64 // encoded ticks
	File                    int
	Name                    string
}

// CKPointing is a pointing (myb_ck_pointing): the quaternion of the C-matrix, the angular
// velocity, the instant, the base frame and the rotation base -> instrument.
type CKPointing struct {
	Q        [4]float64
	Rate     [3]float64
	HasRate  bool
	Ticks    float64
	Frame    int32
	Rotation Rotation
}

// NewCK is an empty set.
func NewCK() (*CK, error) {
	var c *C.myb_ck
	if err := fail("NewCK", C.myb_ck_create(nil, &c)); err != nil {
		return nil, err
	}
	held := &images{}
	return &CK{images: held, o: newObject(unsafe.Pointer(c), func(q unsafe.Pointer) {
		C.myb_ck_destroy((*C.myb_ck)(q))
		held.mu.Lock()
		for _, p := range held.list {
			C.free(p)
		}
		held.list = nil
		held.mu.Unlock()
	})}, nil
}

func (c *CK) obj() *object {
	if c == nil {
		return nil
	}
	return c.o
}

// Close releases the set.
func (c *CK) Close() error {
	c.obj().close()
	return nil
}

// LoadFile adds a CK file.
func (c *CK) LoadFile(path string) error {
	ptr, err := c.obj().acquire()
	if err != nil {
		return err
	}
	defer c.o.release()
	cpath := C.CString(path)
	defer C.free(unsafe.Pointer(cpath))
	var diag C.myb_diag
	return failDiag("CK.LoadFile", C.myb_ck_load_file((*C.myb_ck)(ptr), cpath, &diag), &diag)
}

// LoadBytes adds a CK file given as its bytes (copied).
func (c *CK) LoadBytes(data []byte) error {
	ptr, err := c.obj().acquire()
	if err != nil {
		return err
	}
	defer c.o.release()
	if len(data) == 0 {
		return refuse("CK.LoadBytes", -1, "no bytes")
	}
	image := C.CBytes(data)
	var diag C.myb_diag
	status := C.myb_ck_load_memory((*C.myb_ck)(ptr), image, C.size_t(len(data)), &diag)
	if status != C.MYB_OK {
		C.free(image)
		return failDiag("CK.LoadBytes", status, &diag)
	}
	c.images.mu.Lock()
	c.images.list = append(c.images.list, image)
	c.images.mu.Unlock()
	return nil
}

// Segments lists the segments of the set.
func (c *CK) Segments() ([]CKSegmentInfo, error) {
	ptr, err := c.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer c.o.release()
	ck := (*C.myb_ck)(ptr)
	out := make([]CKSegmentInfo, int(C.myb_ck_segment_count(ck)))
	for k := range out {
		var info C.myb_ck_segment_info
		C.myb_ck_segment(ck, C.size_t(k), &info)
		out[k] = CKSegmentInfo{Instrument: int32(info.instrument), Frame: int32(info.frame),
			Type: int32(info._type), HasRates: info.has_rates != 0, Start: float64(info.start),
			Stop: float64(info.stop), File: int(info.file), Name: C.GoString(&info.name[0])}
	}
	return out, nil
}

// Pointings are the pointings of the instrument at the ticks (myb_ck_pointing_at), in one call.
func (c *CK) Pointings(instrument int32, ticks []float64, tolerance float64, needRate bool) ([]CKPointing, error) {
	ptr, err := c.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer c.o.release()
	p := make([]C.myb_ck_pointing, len(ticks))
	var failed C.size_t
	status := C.mybgo_ck_pointings((*C.myb_ck)(ptr), C.int32_t(instrument), cDoubles(ticks),
		C.size_t(len(ticks)), C.double(tolerance), flag(needRate), (*C.myb_ck_pointing)(first(p)), &failed)
	if err := failAt("CK.Pointings", status, failed); err != nil {
		return nil, err
	}
	out := make([]CKPointing, len(p))
	for k, v := range p {
		out[k] = CKPointing{Rate: goVec3(v.rate), HasRate: v.has_rate != 0, Ticks: float64(v.ticks),
			Frame: int32(v.frame), Rotation: *(*Rotation)(unsafe.Pointer(&v.rotation))}
		for i := range 4 {
			out[k].Q[i] = float64(v.q[i])
		}
	}
	return out, nil
}

// CKSegment is a segment to write: type 1 or 3 (pointings at Ticks with Quaternions and
// Rates — nil: none —, type 3 interpolating inside the intervals that start at IntervalStarts,
// nil: one) or type 2 (intervals Ticks[k] .. Stops[k], each with its quaternion, rate and
// SecondsPerTick).
type CKSegment struct {
	Instrument, Frame int32
	Name              string
	Type              int
	Ticks             []float64
	Stops             []float64
	Quaternions       [][4]float64
	Rates             [][3]float64
	SecondsPerTick    []float64
	IntervalStarts    []float64
}

func ckWrite(op, path, internalName string, segments []CKSegment) ([]byte, error) {
	k := newKeeper()
	defer k.done()
	cs := make([]C.mybgo_ck_segment, len(segments))
	var ticks, stops, spt, starts []float64
	var quaternions [][4]float64
	var rates [][3]float64
	for j, s := range segments {
		n := len(s.Ticks)
		if len(s.Quaternions) != n || (s.Rates != nil && len(s.Rates) != n) {
			return nil, refuse(op, j, "ticks, quaternions and rates of different lengths")
		}
		c := &cs[j]
		c.header = C.myb_ck_segment_header{instrument: C.int32_t(s.Instrument), frame: C.int32_t(s.Frame),
			name: k.cstring(s.Name)}
		c.count = C.size_t(n)
		c._type = C.int(s.Type)
		switch s.Type {
		case 1, 3:
			c.kind = C.MYBGO_CK_POINTINGS
			c.has_rates = flag(s.Rates != nil)
			c.interval_count = C.size_t(len(s.IntervalStarts))
			starts = append(starts, s.IntervalStarts...)
			stops = append(stops, make([]float64, n)...)
			spt = append(spt, make([]float64, n)...)
		case 2:
			if len(s.Stops) != n || len(s.SecondsPerTick) != n || len(s.Rates) != n {
				return nil, refuse(op, j, "a segment of type 2 needs stops, rates and seconds a tick")
			}
			c.kind = C.MYBGO_CK_INTERVALS
			stops = append(stops, s.Stops...)
			spt = append(spt, s.SecondsPerTick...)
		default:
			return nil, refuse(op, j, "a segment type other than 1, 2, 3")
		}
		ticks = append(ticks, s.Ticks...)
		quaternions = append(quaternions, s.Quaternions...)
		if s.Rates != nil {
			rates = append(rates, s.Rates...)
		} else {
			rates = append(rates, make([][3]float64, n)...)
		}
	}
	cname := k.cstring(internalName)
	var cpath *C.char
	if path != "" {
		cpath = k.cstring(path)
	}
	var data unsafe.Pointer
	var size, failed C.size_t
	status := C.mybgo_ck_write(cpath, cname, (*C.mybgo_ck_segment)(first(cs)), C.size_t(len(cs)),
		cDoubles(ticks), cDoubles(stops), (*C.double)(first(quaternions)), (*C.double)(first(rates)),
		cDoubles(spt), cDoubles(starts), &data, &size, &failed)
	if status != C.MYB_OK {
		e := &Error{Status: Status(status), Op: op, Index: -1}
		if int(failed) < len(segments) {
			e.Index = int(failed)
		}
		return nil, e
	}
	if data == nil {
		return nil, nil
	}
	defer C.mybgo_free(data)
	return C.GoBytes(data, C.int(size)), nil
}

// WriteCK writes the segments as one CK file at path; internalName labels the file.
func WriteCK(path, internalName string, segments []CKSegment) error {
	if path == "" {
		return refuse("WriteCK", -1, "no path")
	}
	_, err := ckWrite("WriteCK", path, internalName, segments)
	return err
}

// EncodeCK is the CK file of the segments as bytes.
func EncodeCK(internalName string, segments []CKSegment) ([]byte, error) {
	return ckWrite("EncodeCK", "", internalName, segments)
}
