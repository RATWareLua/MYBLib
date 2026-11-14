// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include <stdlib.h>
#include "glue.h"
*/
import "C"

import "unsafe"

// Time is an instant of a time scale (myb_time): whole seconds S since J2000.0 of the scale
// and a fraction F in [0, 1).
type Time struct {
	S int64
	F float64
}

// Calendar is a date and time of day on the proleptic Gregorian calendar (myb_calendar),
// with astronomical numbering of years.
type Calendar struct {
	Year   int64
	Month  int32
	Day    int32
	Hour   int32
	Minute int32
	Second float64
}

// The range of the instants and seconds the library takes (its contracts).
const (
	maxWhole   = int64(1) << 62
	maxSeconds = 0x1p62
	maxExact   = int64(1) << 53
)

// validSeconds is false for NaN and the infinities too (every comparison with NaN fails).
func validSeconds(s float64) bool {
	return s > -maxSeconds && s < maxSeconds
}

func validTime(t Time) bool {
	return t.S > -maxWhole && t.S < maxWhole && t.F >= 0 && t.F < 1
}

func checkTimes(op string, t []Time) error {
	for k, v := range t {
		if !validTime(v) {
			return refuse(op, k, "an instant out of range")
		}
	}
	return nil
}

// TimeFromSeconds is the instant of the seconds since J2000.0 (myb_time_from_seconds):
// MYB_E_DOMAIN for seconds not finite or beyond 2^62.
func TimeFromSeconds(seconds float64) (Time, error) {
	if !validSeconds(seconds) {
		return Time{}, refuse("TimeFromSeconds", -1, "seconds not finite or beyond 2^62")
	}
	t := C.myb_time_from_seconds(C.double(seconds))
	return Time{S: int64(t.s), F: float64(t.f)}, nil
}

// TimesFromSeconds is TimeFromSeconds over a slice, in one call.
func TimesFromSeconds(seconds []float64) ([]Time, error) {
	for k, s := range seconds {
		if !validSeconds(s) {
			return nil, refuse("TimesFromSeconds", k, "seconds not finite or beyond 2^62")
		}
	}
	out := make([]Time, len(seconds))
	C.myb_time_from_seconds_many(cDoubles(seconds), C.size_t(len(seconds)), cTimes(out))
	return out, nil
}

// Seconds is the instant as seconds since J2000.0 rounded to a double (myb_time_seconds).
func Seconds(t Time) (float64, error) {
	if !validTime(t) || t.S <= -maxExact || t.S >= maxExact {
		return 0, refuse("Seconds", -1, "an instant out of range")
	}
	return float64(C.myb_time_seconds(C.myb_time{s: C.int64_t(t.S), f: C.double(t.F)})), nil
}

// SecondsOf is Seconds over a slice, in one call.
func SecondsOf(t []Time) ([]float64, error) {
	for k, v := range t {
		if !validTime(v) || v.S <= -maxExact || v.S >= maxExact {
			return nil, refuse("SecondsOf", k, "an instant out of range")
		}
	}
	out := make([]float64, len(t))
	C.myb_time_seconds_many(cTimes(t), C.size_t(len(t)), cDoubles(out))
	return out, nil
}

// CalendarToTime is the instant of a date and time of a continuous scale
// (myb_calendar_to_time).
func CalendarToTime(c Calendar) (Time, error) {
	var t C.myb_time
	status := C.myb_calendar_to_time((*C.myb_calendar)(unsafe.Pointer(&c)), &t)
	return Time{S: int64(t.s), F: float64(t.f)}, fail("CalendarToTime", status)
}

// CalendarsToTimes is CalendarToTime over a slice, in one call.
func CalendarsToTimes(c []Calendar) ([]Time, error) {
	out := make([]Time, len(c))
	var failed C.size_t
	status := C.myb_calendar_to_time_many((*C.myb_calendar)(first(c)), C.size_t(len(c)),
		cTimes(out), &failed)
	if err := failAt("CalendarsToTimes", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// TimeToCalendar is the date and time of an instant of a continuous scale
// (myb_time_to_calendar).
func TimeToCalendar(t Time) (Calendar, error) {
	if !validTime(t) {
		return Calendar{}, refuse("TimeToCalendar", -1, "an instant out of range")
	}
	var c Calendar
	C.myb_time_to_calendar(C.myb_time{s: C.int64_t(t.S), f: C.double(t.F)},
		(*C.myb_calendar)(unsafe.Pointer(&c)))
	return c, nil
}

// TimesToCalendars is TimeToCalendar over a slice, in one call.
func TimesToCalendars(t []Time) ([]Calendar, error) {
	if err := checkTimes("TimesToCalendars", t); err != nil {
		return nil, err
	}
	out := make([]Calendar, len(t))
	C.myb_time_to_calendar_many(cTimes(t), C.size_t(len(t)), (*C.myb_calendar)(first(out)))
	return out, nil
}

// ---- Kernel pools --------------------------------------------------------------------------

// KernelPool holds the variables of NAIF text kernels (myb_kpool): the time scales, the
// models of the Earth's orientation and the rotation models are made from them.
type KernelPool struct {
	o *object
}

// NewKernelPool is an empty pool.
func NewKernelPool() (*KernelPool, error) {
	var p *C.myb_kpool
	if err := fail("NewKernelPool", C.myb_kpool_create(nil, &p)); err != nil {
		return nil, err
	}
	return &KernelPool{o: newObject(unsafe.Pointer(p), func(q unsafe.Pointer) {
		C.myb_kpool_destroy((*C.myb_kpool)(q))
	})}, nil
}

func (p *KernelPool) obj() *object {
	if p == nil {
		return nil
	}
	return p.o
}

// Close releases the pool.
func (p *KernelPool) Close() error {
	p.obj().close()
	return nil
}

// LoadFile adds the variables of a kernel file (path in UTF-8).
func (p *KernelPool) LoadFile(path string) error {
	ptr, err := p.obj().acquire()
	if err != nil {
		return err
	}
	defer p.o.release()
	cpath := C.CString(path)
	defer C.free(unsafe.Pointer(cpath))
	var diag C.myb_diag
	return failDiag("KernelPool.LoadFile", C.myb_kpool_load_file((*C.myb_kpool)(ptr), cpath, &diag),
		&diag)
}

// LoadText adds the variables of a kernel given as text.
func (p *KernelPool) LoadText(text string) error {
	ptr, err := p.obj().acquire()
	if err != nil {
		return err
	}
	defer p.o.release()
	ctext := C.CString(text)
	defer C.free(unsafe.Pointer(ctext))
	var diag C.myb_diag
	return failDiag("KernelPool.LoadText", C.myb_kpool_load_text((*C.myb_kpool)(ptr), ctext,
		C.size_t(len(text)), &diag), &diag)
}

// ---- Time systems --------------------------------------------------------------------------

// TimeSystem is a tree of time scales (myb_timesys) built from kernel variables
// (data/time/time-scales.tk, and data/time/utc-tai.tk for UTC). Scales are named as the
// kernels name them ("TT", "TDB", "TAI", ...).
type TimeSystem struct {
	o *object
}

// NewTimeSystem builds the scales from a pool; the pool may be closed afterwards.
func NewTimeSystem(pool *KernelPool) (*TimeSystem, error) {
	pp, err := pool.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer pool.o.release()
	var s *C.myb_timesys
	var diag C.myb_diag
	status := C.myb_timesys_create((*C.myb_kpool)(pp), nil, &s, &diag)
	if err := failDiag("NewTimeSystem", status, &diag); err != nil {
		return nil, err
	}
	return &TimeSystem{o: newObject(unsafe.Pointer(s), func(q unsafe.Pointer) {
		C.myb_timesys_destroy((*C.myb_timesys)(q))
	})}, nil
}

func (s *TimeSystem) obj() *object {
	if s == nil {
		return nil
	}
	return s.o
}

// Close releases the system (after the models that use it).
func (s *TimeSystem) Close() error {
	s.obj().close()
	return nil
}

func scaleOf(op string, s *C.myb_timesys, name string) (C.myb_scale, error) {
	cname := C.CString(name)
	defer C.free(unsafe.Pointer(cname))
	var scale C.myb_scale
	return scale, fail(op, C.myb_timesys_scale(s, cname, &scale))
}

// scales resolves the names of two scales of the system.
func (s *TimeSystem) scales(op string, from, to string) (*C.myb_timesys, C.myb_scale, C.myb_scale, error) {
	ptr, err := s.obj().acquire()
	if err != nil {
		return nil, C.myb_scale{}, C.myb_scale{}, err
	}
	sys := (*C.myb_timesys)(ptr)
	a, err := scaleOf(op, sys, from)
	if err == nil && to != "" {
		var b C.myb_scale
		b, err = scaleOf(op, sys, to)
		if err == nil {
			return sys, a, b, nil
		}
	}
	if err != nil {
		s.o.release()
		return nil, C.myb_scale{}, C.myb_scale{}, err
	}
	return sys, a, C.myb_scale{}, nil
}

// Convert is the same instant on another scale (myb_timesys_convert).
func (s *TimeSystem) Convert(t Time, from, to string) (Time, error) {
	if !validTime(t) {
		return Time{}, refuse("TimeSystem.Convert", -1, "an instant out of range")
	}
	sys, a, b, err := s.scales("TimeSystem.Convert", from, to)
	if err != nil {
		return Time{}, err
	}
	defer s.o.release()
	var out C.myb_time
	status := C.myb_timesys_convert(sys, C.myb_time{s: C.int64_t(t.S), f: C.double(t.F)}, a, b, &out)
	return Time{S: int64(out.s), F: float64(out.f)}, fail("TimeSystem.Convert", status)
}

// ConvertMany is Convert over a slice, in one call (myb_timesys_convert_many).
func (s *TimeSystem) ConvertMany(t []Time, from, to string) ([]Time, error) {
	if err := checkTimes("TimeSystem.ConvertMany", t); err != nil {
		return nil, err
	}
	sys, a, b, err := s.scales("TimeSystem.ConvertMany", from, to)
	if err != nil {
		return nil, err
	}
	defer s.o.release()
	out := make([]Time, len(t))
	var failed C.size_t
	status := C.myb_timesys_convert_many(sys, cTimes(t), C.size_t(len(t)), a, b, cTimes(out), &failed)
	if err := failAt("TimeSystem.ConvertMany", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// FromUTC is the instant on the scale `to` of a UTC label (myb_timesys_from_utc).
func (s *TimeSystem) FromUTC(utc Calendar, to string) (Time, error) {
	sys, scale, _, err := s.scales("TimeSystem.FromUTC", to, "")
	if err != nil {
		return Time{}, err
	}
	defer s.o.release()
	var out C.myb_time
	status := C.myb_timesys_from_utc(sys, (*C.myb_calendar)(unsafe.Pointer(&utc)), scale, &out)
	return Time{S: int64(out.s), F: float64(out.f)}, fail("TimeSystem.FromUTC", status)
}

// FromUTCMany is FromUTC over a slice, in one call (myb_timesys_from_utc_many).
func (s *TimeSystem) FromUTCMany(utc []Calendar, to string) ([]Time, error) {
	sys, scale, _, err := s.scales("TimeSystem.FromUTCMany", to, "")
	if err != nil {
		return nil, err
	}
	defer s.o.release()
	out := make([]Time, len(utc))
	var failed C.size_t
	status := C.myb_timesys_from_utc_many(sys, (*C.myb_calendar)(first(utc)), C.size_t(len(utc)),
		scale, cTimes(out), &failed)
	if err := failAt("TimeSystem.FromUTCMany", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// ToUTC is the UTC label of an instant of the scale `from` (myb_timesys_to_utc).
func (s *TimeSystem) ToUTC(t Time, from string) (Calendar, error) {
	if !validTime(t) {
		return Calendar{}, refuse("TimeSystem.ToUTC", -1, "an instant out of range")
	}
	sys, scale, _, err := s.scales("TimeSystem.ToUTC", from, "")
	if err != nil {
		return Calendar{}, err
	}
	defer s.o.release()
	var c Calendar
	status := C.myb_timesys_to_utc(sys, C.myb_time{s: C.int64_t(t.S), f: C.double(t.F)}, scale,
		(*C.myb_calendar)(unsafe.Pointer(&c)))
	return c, fail("TimeSystem.ToUTC", status)
}

// ToUTCMany is ToUTC over a slice, in one call (myb_timesys_to_utc_many).
func (s *TimeSystem) ToUTCMany(t []Time, from string) ([]Calendar, error) {
	if err := checkTimes("TimeSystem.ToUTCMany", t); err != nil {
		return nil, err
	}
	sys, scale, _, err := s.scales("TimeSystem.ToUTCMany", from, "")
	if err != nil {
		return nil, err
	}
	defer s.o.release()
	out := make([]Calendar, len(t))
	var failed C.size_t
	status := C.myb_timesys_to_utc_many(sys, cTimes(t), C.size_t(len(t)), scale,
		(*C.myb_calendar)(first(out)), &failed)
	if err := failAt("TimeSystem.ToUTCMany", status, failed); err != nil {
		return nil, err
	}
	return out, nil
}
