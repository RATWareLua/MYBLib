// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#cgo pkg-config: myb-static
#include "glue.h"
*/
import "C"

import "unsafe"

// The Go records shared with C have the layout of the C structures (glue.h): an array of
// lengths that are negative unless the two agree fails to compile.
var (
	_ [unsafe.Sizeof(Time{}) - uintptr(C.MYBGO_SIZE_TIME)]struct{}
	_ [uintptr(C.MYBGO_SIZE_TIME) - unsafe.Sizeof(Time{})]struct{}
	_ [unsafe.Sizeof(Calendar{}) - uintptr(C.MYBGO_SIZE_CALENDAR)]struct{}
	_ [uintptr(C.MYBGO_SIZE_CALENDAR) - unsafe.Sizeof(Calendar{})]struct{}
	_ [unsafe.Offsetof(Calendar{}.Month) - uintptr(C.MYBGO_OFFSET_CALENDAR_MONTH)]struct{}
	_ [uintptr(C.MYBGO_OFFSET_CALENDAR_MONTH) - unsafe.Offsetof(Calendar{}.Month)]struct{}
	_ [unsafe.Offsetof(Calendar{}.Second) - uintptr(C.MYBGO_OFFSET_CALENDAR_SECOND)]struct{}
	_ [uintptr(C.MYBGO_OFFSET_CALENDAR_SECOND) - unsafe.Offsetof(Calendar{}.Second)]struct{}
	_ [unsafe.Sizeof(State{}) - uintptr(C.MYBGO_SIZE_PV)]struct{}
	_ [uintptr(C.MYBGO_SIZE_PV) - unsafe.Sizeof(State{})]struct{}
	_ [unsafe.Sizeof(Elements{}) - uintptr(C.MYBGO_SIZE_ELEMENTS)]struct{}
	_ [uintptr(C.MYBGO_SIZE_ELEMENTS) - unsafe.Sizeof(Elements{})]struct{}
	_ [unsafe.Sizeof(Equinoctial{}) - uintptr(C.MYBGO_SIZE_EQUINOCTIAL)]struct{}
	_ [uintptr(C.MYBGO_SIZE_EQUINOCTIAL) - unsafe.Sizeof(Equinoctial{})]struct{}
	_ [unsafe.Sizeof(Rotation{}) - uintptr(C.MYBGO_SIZE_ROTATION)]struct{}
	_ [uintptr(C.MYBGO_SIZE_ROTATION) - unsafe.Sizeof(Rotation{})]struct{}
	_ [unsafe.Sizeof(EOPValues{}) - uintptr(C.MYBGO_SIZE_EOP)]struct{}
	_ [uintptr(C.MYBGO_SIZE_EOP) - unsafe.Sizeof(EOPValues{})]struct{}
	_ [unsafe.Sizeof(Orientation{}) - uintptr(C.MYBGO_SIZE_ORIENTATION)]struct{}
	_ [uintptr(C.MYBGO_SIZE_ORIENTATION) - unsafe.Sizeof(Orientation{})]struct{}
	_ [unsafe.Offsetof(Orientation{}.UT) - uintptr(C.MYBGO_OFFSET_ORIENTATION_UT)]struct{}
	_ [uintptr(C.MYBGO_OFFSET_ORIENTATION_UT) - unsafe.Offsetof(Orientation{}.UT)]struct{}
	_ [unsafe.Offsetof(Orientation{}.XP) - uintptr(C.MYBGO_OFFSET_ORIENTATION_XP)]struct{}
	_ [uintptr(C.MYBGO_OFFSET_ORIENTATION_XP) - unsafe.Offsetof(Orientation{}.XP)]struct{}
)

// Version is the version of the linked library in the format of MYB_VERSION_NUMBER
// (major.minor.patch.tweak, two digits each: 1.2.2.0 is 1020200).
func Version() int {
	return int(C.myb_version_number())
}

// The pointer to the first element of a slice for C, or nil for an empty one (the library
// takes NULL with a count of 0). The slices given are of Go types laid out as the C ones and
// hold no Go pointers.
func first[T any](s []T) unsafe.Pointer {
	if len(s) == 0 {
		return nil
	}
	return unsafe.Pointer(&s[0])
}

func cTimes(s []Time) *C.myb_time {
	return (*C.myb_time)(first(s))
}

func cDoubles(s []float64) *C.double {
	return (*C.double)(first(s))
}

func cStates(s []State) *C.double {
	return (*C.double)(first(s))
}
