// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include <stdlib.h>
#include "glue.h"
*/
import "C"

import "unsafe"

// GravityField is the spherical-harmonic field of a body (myb_gravity_field).
type GravityField struct {
	o *object
}

// The largest degree the library takes.
const maxDegree = 2700

// NewGravityField is a field of degree 0 .. 2700 from fully normalised coefficients indexed
// n (n + 1) / 2 + m, 0 <= m <= n <= degree (myb_gravity_field_create); c and s are copied.
func NewGravityField(gm, radius float64, degree int, c, s []float64) (*GravityField, error) {
	if degree < 0 || degree > maxDegree {
		return nil, refuse("NewGravityField", -1, "a degree outside 0 .. 2700")
	}
	if need := (degree + 1) * (degree + 2) / 2; len(c) < need || len(s) < need {
		return nil, refuse("NewGravityField", -1, "fewer coefficients than the degree needs")
	}
	var f *C.myb_gravity_field
	status := C.myb_gravity_field_create(C.double(gm), C.double(radius), C.int(degree),
		cDoubles(c), cDoubles(s), nil, &f)
	if err := fail("NewGravityField", status); err != nil {
		return nil, err
	}
	return ownField(f), nil
}

// LoadICGEMFile reads a field from an ICGEM file, truncated to maxDegree (0: all)
// (myb_gravity_field_load_icgem_file).
func LoadICGEMFile(path string, maxDegree int) (*GravityField, error) {
	cpath := C.CString(path)
	defer C.free(unsafe.Pointer(cpath))
	var f *C.myb_gravity_field
	var diag C.myb_diag
	status := C.myb_gravity_field_load_icgem_file(cpath, C.int(maxDegree), nil, &f, &diag)
	if err := failDiag("LoadICGEMFile", status, &diag); err != nil {
		return nil, err
	}
	return ownField(f), nil
}

func ownField(f *C.myb_gravity_field) *GravityField {
	return &GravityField{o: newObject(unsafe.Pointer(f), func(q unsafe.Pointer) {
		C.myb_gravity_field_destroy((*C.myb_gravity_field)(q))
	})}
}

func (f *GravityField) obj() *object {
	if f == nil {
		return nil
	}
	return f.o
}

// Close releases the field (after the models that use it).
func (f *GravityField) Close() error {
	f.obj().close()
	return nil
}

// Degree is the degree of the field.
func (f *GravityField) Degree() (int, error) {
	ptr, err := f.obj().acquire()
	if err != nil {
		return 0, err
	}
	defer f.o.release()
	return int(C.myb_gravity_field_degree((*C.myb_gravity_field)(ptr))), nil
}
