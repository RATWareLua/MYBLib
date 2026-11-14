// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include <stdlib.h>
#include "glue.h"
*/
import "C"

import "unsafe"

// Rotation maps coordinates of frame A to frame B with its time derivative (myb_rotation):
// the rows of M and of dM/dt (s^-1).
type Rotation struct {
	M  [9]float64
	DM [9]float64
}

// EOPValues are the Earth orientation parameters at one instant (myb_eop_values): polar
// motion XP, YP (rad), UT1 - UTC (s), the excess length of day (s), the celestial pole
// offsets DX, DY (rad).
type EOPValues struct {
	XP, YP, DUT1, LOD, DX, DY float64
}

// Orientation is what an orientation provider gives (myb_orientation): the rotation from the
// inertial axes to the body-fixed frame, the body's rotation time (UT1 for the Earth) and the
// coordinates of its pole.
type Orientation struct {
	Rotation Rotation
	UT       Time
	XP, YP   float64
}

func cTime(t Time) C.myb_time {
	return C.myb_time{s: C.int64_t(t.S), f: C.double(t.F)}
}

func cRotation(r *Rotation) *C.myb_rotation {
	return (*C.myb_rotation)(unsafe.Pointer(r))
}

func cRotations(r []Rotation) *C.myb_rotation {
	return (*C.myb_rotation)(first(r))
}

// Apply is the state s carried by the rotation: (M r, M v + dM/dt r)
// (myb_rotation_apply_state).
func (r Rotation) Apply(s State) State {
	var out State
	C.myb_rotation_apply_state(cRotation(&r), cStates([]State{s}), (*C.double)(unsafe.Pointer(&out)))
	return out
}

// ApplyRotations carries states[k] by rotations[k], in one call (myb_rotation_apply_state_many).
func ApplyRotations(rotations []Rotation, states []State) ([]State, error) {
	if len(rotations) != len(states) {
		return nil, refuse("ApplyRotations", -1, "rotations and states of different lengths")
	}
	out := make([]State, len(states))
	C.myb_rotation_apply_state_many(cRotations(rotations), cStates(states), C.size_t(len(states)),
		cStates(out))
	return out, nil
}

// ---- The celestial intermediate pole -----------------------------------------------------------

// CIP is the CIO-based model of the Earth's orientation (myb_cip, IERS Conventions 2010),
// read from kernel variables (data/frames/iau2006-cip.tk with the prefix "IAU2006A").
type CIP struct {
	o *object
}

// NewCIP reads the model <prefix>_* of a pool; the pool may be closed afterwards.
func NewCIP(pool *KernelPool, prefix string) (*CIP, error) {
	pp, err := pool.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer pool.o.release()
	cprefix := C.CString(prefix)
	defer C.free(unsafe.Pointer(cprefix))
	var c *C.myb_cip
	var diag C.myb_diag
	if err := failDiag("NewCIP", C.myb_cip_create((*C.myb_kpool)(pp), cprefix, nil, &c, &diag),
		&diag); err != nil {
		return nil, err
	}
	return &CIP{o: newObject(unsafe.Pointer(c), func(q unsafe.Pointer) {
		C.myb_cip_destroy((*C.myb_cip)(q))
	})}, nil
}

func (c *CIP) obj() *object {
	if c == nil {
		return nil
	}
	return c.o
}

// Close releases the model (after the models of forces that use it).
func (c *CIP) Close() error {
	c.obj().close()
	return nil
}

// Tabulate replaces the series by a table every step seconds of TT from tt0 to tt1
// (myb_cip_tabulate).
func (c *CIP) Tabulate(tt0, tt1 Time, step float64) error {
	if !validTime(tt0) || !validTime(tt1) {
		return refuse("CIP.Tabulate", -1, "an instant out of range")
	}
	ptr, err := c.obj().acquire()
	if err != nil {
		return err
	}
	defer c.o.release()
	return fail("CIP.Tabulate", C.myb_cip_tabulate((*C.myb_cip)(ptr), cTime(tt0), cTime(tt1),
		C.double(step)))
}

// XYS are the coordinates X, Y of the pole and the locator s (rad) at a TT instant, with
// their rates (myb_cip_xys).
func (c *CIP) XYS(tt Time) (xys, rates [3]float64, err error) {
	if !validTime(tt) {
		return xys, rates, refuse("CIP.XYS", -1, "an instant out of range")
	}
	ptr, err := c.obj().acquire()
	if err != nil {
		return xys, rates, err
	}
	defer c.o.release()
	C.myb_cip_xys((*C.myb_cip)(ptr), cTime(tt), (*C.double)(unsafe.Pointer(&xys[0])),
		(*C.double)(unsafe.Pointer(&rates[0])))
	return xys, rates, nil
}

// XYSMany is XYS at every instant of a slice, in one call (myb_cip_xys_many).
func (c *CIP) XYSMany(tt []Time) (xys, rates [][3]float64, err error) {
	if err := checkTimes("CIP.XYSMany", tt); err != nil {
		return nil, nil, err
	}
	ptr, err := c.obj().acquire()
	if err != nil {
		return nil, nil, err
	}
	defer c.o.release()
	xys = make([][3]float64, len(tt))
	rates = make([][3]float64, len(tt))
	C.myb_cip_xys_many((*C.myb_cip)(ptr), cTimes(tt), C.size_t(len(tt)), (*C.double)(first(xys)),
		(*C.double)(first(rates)))
	return xys, rates, nil
}

// Terrestrial is the rotation celestial -> terrestrial at the instant given as TT and UT1,
// with the EOP (nil: all zero) (myb_cip_terrestrial).
func (c *CIP) Terrestrial(tt, ut1 Time, eop *EOPValues) (Rotation, error) {
	if !validTime(tt) || !validTime(ut1) {
		return Rotation{}, refuse("CIP.Terrestrial", -1, "an instant out of range")
	}
	ptr, err := c.obj().acquire()
	if err != nil {
		return Rotation{}, err
	}
	defer c.o.release()
	var r Rotation
	C.myb_cip_terrestrial((*C.myb_cip)(ptr), cTime(tt), cTime(ut1),
		(*C.myb_eop_values)(unsafe.Pointer(eop)), cRotation(&r))
	return r, nil
}

// TerrestrialMany is Terrestrial at every instant (tt[k], ut1[k]) with eop[k] (eop may be
// nil: all zero), in one call (myb_cip_terrestrial_many).
func (c *CIP) TerrestrialMany(tt, ut1 []Time, eop []EOPValues) ([]Rotation, error) {
	if len(ut1) != len(tt) || (eop != nil && len(eop) != len(tt)) {
		return nil, refuse("CIP.TerrestrialMany", -1, "slices of different lengths")
	}
	if err := checkTimes("CIP.TerrestrialMany", tt); err != nil {
		return nil, err
	}
	if err := checkTimes("CIP.TerrestrialMany", ut1); err != nil {
		return nil, err
	}
	ptr, err := c.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer c.o.release()
	out := make([]Rotation, len(tt))
	C.myb_cip_terrestrial_many((*C.myb_cip)(ptr), cTimes(tt), cTimes(ut1),
		(*C.myb_eop_values)(first(eop)), C.size_t(len(tt)), cRotations(out))
	return out, nil
}

// ---- Earth orientation parameters ---------------------------------------------------------------

// EOP is a daily series of Earth orientation parameters with its interpolation (myb_eop).
type EOP struct {
	o *object
}

// NewEOP is an empty series.
func NewEOP() (*EOP, error) {
	var e *C.myb_eop
	if err := fail("NewEOP", C.myb_eop_create(nil, &e)); err != nil {
		return nil, err
	}
	return &EOP{o: newObject(unsafe.Pointer(e), func(q unsafe.Pointer) {
		C.myb_eop_destroy((*C.myb_eop)(q))
	})}, nil
}

func (e *EOP) obj() *object {
	if e == nil {
		return nil
	}
	return e.o
}

// Close releases the series (after the models that use it).
func (e *EOP) Close() error {
	e.obj().close()
	return nil
}

// LoadC04File replaces the table by an EOP 20 C04 file.
func (e *EOP) LoadC04File(path string) error {
	return e.load("EOP.LoadC04File", path, false)
}

// LoadFinalsFile replaces the table by a finals2000A file.
func (e *EOP) LoadFinalsFile(path string) error {
	return e.load("EOP.LoadFinalsFile", path, true)
}

func (e *EOP) load(op, path string, finals bool) error {
	ptr, err := e.obj().acquire()
	if err != nil {
		return err
	}
	defer e.o.release()
	cpath := C.CString(path)
	defer C.free(unsafe.Pointer(cpath))
	var diag C.myb_diag
	var status C.myb_status
	if finals {
		status = C.myb_eop_load_finals_file((*C.myb_eop)(ptr), cpath, &diag)
	} else {
		status = C.myb_eop_load_c04_file((*C.myb_eop)(ptr), cpath, &diag)
	}
	return failDiag(op, status, &diag)
}

// SetTides adds the diurnal and semidiurnal variations of kernel variables <prefix>_*
// (data/frames/eop-tides.tk, prefix "IERS10_EOPT").
func (e *EOP) SetTides(pool *KernelPool, prefix string) error {
	held, err := acquireAll(e.obj(), pool.obj())
	if err != nil {
		return err
	}
	defer releaseAll(held)
	if len(held) != 2 {
		return refuse("EOP.SetTides", -1, "no pool")
	}
	cprefix := C.CString(prefix)
	defer C.free(unsafe.Pointer(cprefix))
	var diag C.myb_diag
	return failDiag("EOP.SetTides", C.myb_eop_set_tides((*C.myb_eop)(e.o.ptr),
		(*C.myb_kpool)(pool.o.ptr), cprefix, &diag), &diag)
}

// At gives the EOP at the instant t of a scale of the system (e.g. "TT") and the same instant
// on UT1 (myb_eop_at).
func (e *EOP) At(system *TimeSystem, t Time, scale string) (EOPValues, Time, error) {
	if !validTime(t) {
		return EOPValues{}, Time{}, refuse("EOP.At", -1, "an instant out of range")
	}
	held, err := acquireAll(e.obj(), system.obj())
	if err != nil {
		return EOPValues{}, Time{}, err
	}
	defer releaseAll(held)
	if len(held) != 2 {
		return EOPValues{}, Time{}, refuse("EOP.At", -1, "no time system")
	}
	sys := (*C.myb_timesys)(system.o.ptr)
	sc, err := scaleOf("EOP.At", sys, scale)
	if err != nil {
		return EOPValues{}, Time{}, err
	}
	var v EOPValues
	var ut1 C.myb_time
	status := C.myb_eop_at((*C.myb_eop)(e.o.ptr), sys, cTime(t), sc,
		(*C.myb_eop_values)(unsafe.Pointer(&v)), &ut1)
	return v, Time{S: int64(ut1.s), F: float64(ut1.f)}, fail("EOP.At", status)
}

// AtMany is At at every instant of a slice, in one call (myb_eop_at_many).
func (e *EOP) AtMany(system *TimeSystem, t []Time, scale string) ([]EOPValues, []Time, error) {
	if err := checkTimes("EOP.AtMany", t); err != nil {
		return nil, nil, err
	}
	held, err := acquireAll(e.obj(), system.obj())
	if err != nil {
		return nil, nil, err
	}
	defer releaseAll(held)
	if len(held) != 2 {
		return nil, nil, refuse("EOP.AtMany", -1, "no time system")
	}
	sys := (*C.myb_timesys)(system.o.ptr)
	sc, err := scaleOf("EOP.AtMany", sys, scale)
	if err != nil {
		return nil, nil, err
	}
	values := make([]EOPValues, len(t))
	ut1 := make([]Time, len(t))
	var failed C.size_t
	status := C.myb_eop_at_many((*C.myb_eop)(e.o.ptr), sys, cTimes(t), C.size_t(len(t)), sc,
		(*C.myb_eop_values)(first(values)), cTimes(ut1), &failed)
	if err := failAt("EOP.AtMany", status, failed); err != nil {
		return nil, nil, err
	}
	return values, ut1, nil
}

// ---- Equinox-based systems ----------------------------------------------------------------------

// EquinoxFrame names a frame of the equinox-based systems (myb_equinox_frame).
type EquinoxFrame int

// The frames reached from J2000 by Equinox.Rotation.
const (
	MOD  EquinoxFrame = C.MYB_EQUINOX_MOD
	TOD  EquinoxFrame = C.MYB_EQUINOX_TOD
	TEME EquinoxFrame = C.MYB_EQUINOX_TEME
)

// Equinox is the classical model of precession and nutation (myb_equinox, IAU 1976/1980),
// read from kernel variables (data/frames/iau1980.tk, prefix "IAU1980").
type Equinox struct {
	o *object
}

// NewEquinox reads the model <prefix>_* of a pool.
func NewEquinox(pool *KernelPool, prefix string) (*Equinox, error) {
	pp, err := pool.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer pool.o.release()
	cprefix := C.CString(prefix)
	defer C.free(unsafe.Pointer(cprefix))
	var q *C.myb_equinox
	var diag C.myb_diag
	if err := failDiag("NewEquinox", C.myb_equinox_create((*C.myb_kpool)(pp), cprefix, nil, &q,
		&diag), &diag); err != nil {
		return nil, err
	}
	return &Equinox{o: newObject(unsafe.Pointer(q), func(p unsafe.Pointer) {
		C.myb_equinox_destroy((*C.myb_equinox)(p))
	})}, nil
}

func (q *Equinox) obj() *object {
	if q == nil {
		return nil
	}
	return q.o
}

// Close releases the model.
func (q *Equinox) Close() error {
	q.obj().close()
	return nil
}

func validFrame(frame EquinoxFrame) bool {
	return frame == MOD || frame == TOD || frame == TEME
}

// Rotation is the rotation J2000 -> frame at a TT instant (myb_equinox_rotation).
func (q *Equinox) Rotation(frame EquinoxFrame, tt Time) (Rotation, error) {
	if !validFrame(frame) || !validTime(tt) {
		return Rotation{}, refuse("Equinox.Rotation", -1, "an unknown frame or an instant out of range")
	}
	ptr, err := q.obj().acquire()
	if err != nil {
		return Rotation{}, err
	}
	defer q.o.release()
	var r Rotation
	C.myb_equinox_rotation((*C.myb_equinox)(ptr), C.myb_equinox_frame(frame), cTime(tt), cRotation(&r))
	return r, nil
}

// RotationMany is Rotation at every instant of a slice, in one call
// (myb_equinox_rotation_many).
func (q *Equinox) RotationMany(frame EquinoxFrame, tt []Time) ([]Rotation, error) {
	if !validFrame(frame) {
		return nil, refuse("Equinox.RotationMany", -1, "an unknown frame")
	}
	if err := checkTimes("Equinox.RotationMany", tt); err != nil {
		return nil, err
	}
	ptr, err := q.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer q.o.release()
	out := make([]Rotation, len(tt))
	C.myb_equinox_rotation_many((*C.myb_equinox)(ptr), C.myb_equinox_frame(frame), cTimes(tt),
		C.size_t(len(tt)), cRotations(out))
	return out, nil
}

// TEMEPEF is the rotation TEME -> PEF at a UT1 instant (myb_equinox_teme_pef).
func (q *Equinox) TEMEPEF(ut1 Time) (Rotation, error) {
	if !validTime(ut1) {
		return Rotation{}, refuse("Equinox.TEMEPEF", -1, "an instant out of range")
	}
	ptr, err := q.obj().acquire()
	if err != nil {
		return Rotation{}, err
	}
	defer q.o.release()
	var r Rotation
	C.myb_equinox_teme_pef((*C.myb_equinox)(ptr), cTime(ut1), cRotation(&r))
	return r, nil
}

// TEMEPEFMany is TEMEPEF at every instant of a slice, in one call
// (myb_equinox_teme_pef_many).
func (q *Equinox) TEMEPEFMany(ut1 []Time) ([]Rotation, error) {
	if err := checkTimes("Equinox.TEMEPEFMany", ut1); err != nil {
		return nil, err
	}
	ptr, err := q.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer q.o.release()
	out := make([]Rotation, len(ut1))
	C.myb_equinox_teme_pef_many((*C.myb_equinox)(ptr), cTimes(ut1), C.size_t(len(ut1)),
		cRotations(out))
	return out, nil
}

// ---- Rotation models of the WGCCRE --------------------------------------------------------------

// WGCCRE is the rotation model of a body as the IAU WGCCRE publishes it (myb_wgccre), read
// from a text PCK (pck00011.tpc) in a pool.
type WGCCRE struct {
	o *object
}

// NewWGCCRE reads the model of body (its NAIF code) from a pool.
func NewWGCCRE(pool *KernelPool, body int64) (*WGCCRE, error) {
	pp, err := pool.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer pool.o.release()
	var w *C.myb_wgccre
	var diag C.myb_diag
	if err := failDiag("NewWGCCRE", C.myb_wgccre_create((*C.myb_kpool)(pp), C.int64_t(body), nil,
		&w, &diag), &diag); err != nil {
		return nil, err
	}
	return &WGCCRE{o: newObject(unsafe.Pointer(w), func(q unsafe.Pointer) {
		C.myb_wgccre_destroy((*C.myb_wgccre)(q))
	})}, nil
}

func (w *WGCCRE) obj() *object {
	if w == nil {
		return nil
	}
	return w.o
}

// Close releases the model.
func (w *WGCCRE) Close() error {
	w.obj().close()
	return nil
}

// Rotation is the rotation reference frame -> body-fixed frame at a TDB instant
// (myb_wgccre_rotation).
func (w *WGCCRE) Rotation(t Time) (Rotation, error) {
	if !validTime(t) {
		return Rotation{}, refuse("WGCCRE.Rotation", -1, "an instant out of range")
	}
	ptr, err := w.obj().acquire()
	if err != nil {
		return Rotation{}, err
	}
	defer w.o.release()
	var r Rotation
	C.myb_wgccre_rotation((*C.myb_wgccre)(ptr), cTime(t), cRotation(&r))
	return r, nil
}

// RotationMany is Rotation at every instant of a slice, in one call
// (myb_wgccre_rotation_many).
func (w *WGCCRE) RotationMany(t []Time) ([]Rotation, error) {
	if err := checkTimes("WGCCRE.RotationMany", t); err != nil {
		return nil, err
	}
	ptr, err := w.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer w.o.release()
	out := make([]Rotation, len(t))
	C.myb_wgccre_rotation_many((*C.myb_wgccre)(ptr), cTimes(t), C.size_t(len(t)), cRotations(out))
	return out, nil
}

// ---- Binary PCK --------------------------------------------------------------------------------

// PCK is a set of binary PCK files (myb_pck): the orientation of frames as Euler angles.
type PCK struct {
	o *object
}

// NewPCK is an empty set.
func NewPCK() (*PCK, error) {
	var p *C.myb_pck
	if err := fail("NewPCK", C.myb_pck_create(nil, &p)); err != nil {
		return nil, err
	}
	return &PCK{o: newObject(unsafe.Pointer(p), func(q unsafe.Pointer) {
		C.myb_pck_destroy((*C.myb_pck)(q))
	})}, nil
}

func (p *PCK) obj() *object {
	if p == nil {
		return nil
	}
	return p.o
}

// Close releases the set.
func (p *PCK) Close() error {
	p.obj().close()
	return nil
}

// LoadFile adds a binary PCK file.
func (p *PCK) LoadFile(path string) error {
	ptr, err := p.obj().acquire()
	if err != nil {
		return err
	}
	defer p.o.release()
	cpath := C.CString(path)
	defer C.free(unsafe.Pointer(cpath))
	var diag C.myb_diag
	return failDiag("PCK.LoadFile", C.myb_pck_load_file((*C.myb_pck)(ptr), cpath, &diag), &diag)
}

// Rotation is the rotation from the base frame to frame class id at t, and the code of the
// base (myb_pck_rotation).
func (p *PCK) Rotation(id int32, t Time) (Rotation, int32, error) {
	if !validTime(t) {
		return Rotation{}, 0, refuse("PCK.Rotation", -1, "an instant out of range")
	}
	ptr, err := p.obj().acquire()
	if err != nil {
		return Rotation{}, 0, err
	}
	defer p.o.release()
	var r Rotation
	var base C.int32_t
	status := C.myb_pck_rotation((*C.myb_pck)(ptr), C.int32_t(id), cTime(t), cRotation(&r), &base)
	return r, int32(base), fail("PCK.Rotation", status)
}

// RotationMany is Rotation at every instant of a slice, in one call (myb_pck_rotation_many).
func (p *PCK) RotationMany(id int32, t []Time) ([]Rotation, []int32, error) {
	if err := checkTimes("PCK.RotationMany", t); err != nil {
		return nil, nil, err
	}
	ptr, err := p.obj().acquire()
	if err != nil {
		return nil, nil, err
	}
	defer p.o.release()
	out := make([]Rotation, len(t))
	bases := make([]int32, len(t))
	var failed C.size_t
	status := C.myb_pck_rotation_many((*C.myb_pck)(ptr), C.int32_t(id), cTimes(t), C.size_t(len(t)),
		cRotations(out), (*C.int32_t)(first(bases)), &failed)
	if err := failAt("PCK.RotationMany", status, failed); err != nil {
		return nil, nil, err
	}
	return out, bases, nil
}

// ---- Orientation providers ------------------------------------------------------------------------

// OrientationSource is the orientation of a body of a model, and a provider that can be
// sampled over a series of instants in one call (myb_orientation_many).
type OrientationSource interface {
	// Sample is the orientation at every instant of a slice.
	Sample(t []Time) ([]Orientation, error)
	describe(b *C.mybgo_body, k *keeper) error
}

// FixedOrientation is a constant rotation from the inertial axes to the body-fixed frame.
type FixedOrientation struct {
	Rotation Rotation
}

// EarthOrientation is the CIO-based model with observed EOP (myb_orientation_cip): instants
// on the scale TT of the time system, the rotation celestial -> terrestrial, UT1 and the
// polar motion.
type EarthOrientation struct {
	CIP        *CIP
	EOP        *EOP
	TimeSystem *TimeSystem
	TT         string // the name of TT in the time system
}

// PCKOrientation is a frame class of a binary PCK (myb_orientation_pck); Base (nil: none) is
// the rotation from the inertial axes in use to the base frame of the segments.
type PCKOrientation struct {
	PCK  *PCK
	ID   int32
	Base *Rotation
}

// WGCCREOrientation is a rotation model of the WGCCRE (myb_orientation_wgccre); Base as for
// PCKOrientation, relative to the model's reference frame.
type WGCCREOrientation struct {
	Model *WGCCRE
	Base  *Rotation
}

func (f FixedOrientation) describe(b *C.mybgo_body, k *keeper) error {
	b.orientation_kind = C.MYBGO_ORIENTATION_FIXED
	b.fixed = *cRotation(&f.Rotation)
	return nil
}

func (e EarthOrientation) describe(b *C.mybgo_body, k *keeper) error {
	b.orientation_kind = C.MYBGO_ORIENTATION_CIP
	cip, err := k.hold(e.CIP.obj())
	if err != nil {
		return err
	}
	eop, err := k.hold(e.EOP.obj())
	if err != nil {
		return err
	}
	sys, err := k.hold(e.TimeSystem.obj())
	if err != nil {
		return err
	}
	b.cip = (*C.myb_cip)(cip)
	b.eop = (*C.myb_eop)(eop)
	b.timesys = (*C.myb_timesys)(sys)
	b.tt_scale = k.cstring(e.TT)
	return nil
}

func (p PCKOrientation) describe(b *C.mybgo_body, k *keeper) error {
	b.orientation_kind = C.MYBGO_ORIENTATION_PCK
	pck, err := k.hold(p.PCK.obj())
	if err != nil {
		return err
	}
	b.pck = (*C.myb_pck)(pck)
	b.pck_id = C.int32_t(p.ID)
	if p.Base != nil {
		b.base = *cRotation(p.Base)
		b.has_base = 1
	}
	return nil
}

func (w WGCCREOrientation) describe(b *C.mybgo_body, k *keeper) error {
	b.orientation_kind = C.MYBGO_ORIENTATION_WGCCRE
	model, err := k.hold(w.Model.obj())
	if err != nil {
		return err
	}
	b.wgccre = (*C.myb_wgccre)(model)
	if w.Base != nil {
		b.base = *cRotation(w.Base)
		b.has_base = 1
	}
	return nil
}

// sample calls an orientation provider of the library over the instants; context holds C
// pointers only.
func sample(op string, provider C.myb_orientation_function, context unsafe.Pointer, t []Time) ([]Orientation, error) {
	out := make([]Orientation, len(t))
	var failed C.size_t
	status := C.myb_orientation_many(provider, context, cTimes(t), C.size_t(len(t)),
		(*C.myb_orientation)(first(out)), &failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	return out, nil
}

// Sample gives the rotation at every instant.
func (f FixedOrientation) Sample(t []Time) ([]Orientation, error) {
	if err := checkTimes("FixedOrientation.Sample", t); err != nil {
		return nil, err
	}
	r := f.Rotation
	return sample("FixedOrientation.Sample", C.myb_orientation_function(C.myb_orientation_fixed),
		unsafe.Pointer(&r), t)
}

// Sample gives the orientation of the Earth at every instant of TT.
func (e EarthOrientation) Sample(t []Time) ([]Orientation, error) {
	const op = "EarthOrientation.Sample"
	if err := checkTimes(op, t); err != nil {
		return nil, err
	}
	held, err := acquireAll(e.CIP.obj(), e.EOP.obj(), e.TimeSystem.obj())
	if err != nil {
		return nil, err
	}
	defer releaseAll(held)
	if len(held) != 3 {
		return nil, refuse(op, -1, "the model, the EOP and the time system are all needed")
	}
	sys := (*C.myb_timesys)(e.TimeSystem.o.ptr)
	tt, err := scaleOf(op, sys, e.TT)
	if err != nil {
		return nil, err
	}
	context := C.myb_cip_orientation{cip: (*C.myb_cip)(e.CIP.o.ptr), eop: (*C.myb_eop)(e.EOP.o.ptr),
		timesys: sys, tt: tt}
	return sample(op, C.myb_orientation_function(C.myb_orientation_cip), unsafe.Pointer(&context), t)
}

// Sample gives the orientation of the frame at every instant.
func (p PCKOrientation) Sample(t []Time) ([]Orientation, error) {
	const op = "PCKOrientation.Sample"
	if err := checkTimes(op, t); err != nil {
		return nil, err
	}
	ptr, err := p.PCK.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer p.PCK.o.release()
	var base *C.myb_rotation
	if p.Base != nil {
		base = (*C.myb_rotation)(C.malloc(C.sizeof_myb_rotation))
		defer C.free(unsafe.Pointer(base))
		*base = *cRotation(p.Base)
	}
	context := C.myb_pck_orientation{pck: (*C.myb_pck)(ptr), id: C.int32_t(p.ID), base: base}
	return sample(op, C.myb_orientation_function(C.myb_orientation_pck), unsafe.Pointer(&context), t)
}

// Sample gives the orientation of the body at every instant of TDB.
func (w WGCCREOrientation) Sample(t []Time) ([]Orientation, error) {
	const op = "WGCCREOrientation.Sample"
	if err := checkTimes(op, t); err != nil {
		return nil, err
	}
	ptr, err := w.Model.obj().acquire()
	if err != nil {
		return nil, err
	}
	defer w.Model.o.release()
	var base *C.myb_rotation
	if w.Base != nil {
		base = (*C.myb_rotation)(C.malloc(C.sizeof_myb_rotation))
		defer C.free(unsafe.Pointer(base))
		*base = *cRotation(w.Base)
	}
	context := C.myb_wgccre_orientation{model: (*C.myb_wgccre)(ptr), base: base}
	return sample(op, C.myb_orientation_function(C.myb_orientation_wgccre), unsafe.Pointer(&context), t)
}
