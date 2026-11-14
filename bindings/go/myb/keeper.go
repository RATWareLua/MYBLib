// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include <stdlib.h>
#include "glue.h"
*/
import "C"

import (
	"runtime/cgo"
	"unsafe"
)

// keeper holds what a description needs while it is passed to C: the objects of the library
// it points to, the C strings, the handles of the Go functions, and the record of their
// errors. A call that keeps nothing lets go of it all with done; an object made from the
// description (a model, a system) takes the objects and the handles with keep.
type keeper struct {
	held    []*object
	strings []*C.char
	handles []cgo.Handle
	state   *callState // the first error of the Go functions of the description
}

func newKeeper() *keeper {
	return &keeper{state: &callState{}}
}

func (k *keeper) hold(o *object) (unsafe.Pointer, error) {
	if o == nil {
		return nil, nil // the description lacks the object: the C side refuses it
	}
	p, err := o.acquire()
	if err != nil {
		return nil, err
	}
	k.held = append(k.held, o)
	return p, nil
}

func (k *keeper) cstring(s string) *C.char {
	c := C.CString(s)
	k.strings = append(k.strings, c)
	return c
}

func (k *keeper) optional(s string) *C.char {
	if s == "" {
		return nil
	}
	return k.cstring(s)
}

func (k *keeper) handle(value any) C.uintptr_t {
	h := cgo.NewHandle(value)
	k.handles = append(k.handles, h)
	return C.uintptr_t(h)
}

func (k *keeper) freeStrings() {
	for _, s := range k.strings {
		C.free(unsafe.Pointer(s))
	}
	k.strings = nil
}

// done lets go of everything: the end of a call that keeps nothing, or of a description that
// failed.
func (k *keeper) done() {
	k.freeStrings()
	releaseAll(k.held)
	k.held = nil
	for _, h := range k.handles {
		h.Delete()
	}
	k.handles = nil
}

// keep hands the objects and the handles over to the caller (to be released when what was
// made is destroyed) and frees the strings.
func (k *keeper) keep() ([]*object, []cgo.Handle) {
	k.freeStrings()
	held, handles := k.held, k.handles
	k.held, k.handles = nil, nil
	return held, handles
}

// failure is the error of a call whose Go functions may have caused it.
func (k *keeper) failure(op string, status C.myb_status, index int) error {
	if status == C.MYB_OK {
		return nil
	}
	return &Error{Status: Status(status), Op: op, Index: index, Err: k.state.err}
}
