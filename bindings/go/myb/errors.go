// SPDX-License-Identifier: BSD-3-Clause

package myb

/*
#include "glue.h"
*/
import "C"

import (
	"errors"
	"fmt"
)

// Status is a status code of the library (myb_status of myb/core.h).
type Status int

// The statuses of the library; their values are those of myb_status.
const (
	OK       Status = C.MYB_OK
	Domain   Status = C.MYB_E_DOMAIN    // an argument outside the domain of the model
	Range    Status = C.MYB_E_RANGE     // an epoch or argument outside the coverage of the data
	NoConv   Status = C.MYB_E_NOCONV    // an iteration budget exhausted
	StepMin  Status = C.MYB_E_STEP_MIN  // the step of an integrator below its minimum
	NoMem    Status = C.MYB_E_NOMEM     // an allocation failed
	Overflow Status = C.MYB_E_OVERFLOW  // a size or an index overflowed
	IO       Status = C.MYB_E_IO        // a file could not be read or written
	Format   Status = C.MYB_E_FORMAT    // malformed data
	NotFound Status = C.MYB_E_NOT_FOUND // a named item does not exist
)

// String is the name the library gives the status ("MYB_E_RANGE").
func (s Status) String() string {
	return C.GoString(C.myb_status_name(C.myb_status(s)))
}

// Error is a failure of the library.
type Error struct {
	Status  Status
	Op      string // the function of the package that failed
	Index   int    // the element of a slice (or the record of a description) that failed; -1: none
	Line    int    // the line of the text a parser names; 0: none
	Message string // what the parser or the description says
	Err     error  // the error a Go function returned to the library, if that was the cause
}

// The statuses as errors to compare with errors.Is: errors.Is(err, myb.ErrRange) holds for
// any *Error with the status MYB_E_RANGE.
var (
	ErrDomain   = &Error{Status: Domain, Index: -1}
	ErrRange    = &Error{Status: Range, Index: -1}
	ErrNoConv   = &Error{Status: NoConv, Index: -1}
	ErrStepMin  = &Error{Status: StepMin, Index: -1}
	ErrNoMem    = &Error{Status: NoMem, Index: -1}
	ErrOverflow = &Error{Status: Overflow, Index: -1}
	ErrIO       = &Error{Status: IO, Index: -1}
	ErrFormat   = &Error{Status: Format, Index: -1}
	ErrNotFound = &Error{Status: NotFound, Index: -1}
)

// ErrClosed is returned by the methods of an object after its Close.
var ErrClosed = errors.New("myb: the object is closed")

func (e *Error) Error() string {
	text := e.Status.String()
	if e.Op != "" {
		text = "myb." + e.Op + ": " + text
	}
	if e.Index >= 0 {
		text += fmt.Sprintf(" (element %d)", e.Index)
	}
	if e.Line > 0 {
		text += fmt.Sprintf(" at line %d", e.Line)
	}
	if e.Message != "" {
		text += ": " + e.Message
	}
	if e.Err != nil {
		text += ": " + e.Err.Error()
	}
	return text
}

// Unwrap gives the error of a Go function that stopped the library.
func (e *Error) Unwrap() error {
	return e.Err
}

// Is holds for a target that is an *Error of the same status and names no function — one of
// ErrDomain, ErrRange, ...
func (e *Error) Is(target error) bool {
	t, ok := target.(*Error)
	return ok && t.Op == "" && t.Status == e.Status
}

// statusOf gives the status a Go error stands for in the library.
func statusOf(err error) Status {
	var e *Error
	if errors.As(err, &e) && e.Status != OK {
		return e.Status
	}
	return Domain
}

func fail(op string, status C.myb_status) error {
	if status == C.MYB_OK {
		return nil
	}
	return &Error{Status: Status(status), Op: op, Index: -1}
}

// failAt is fail with the index of the element of a slice that failed.
func failAt(op string, status C.myb_status, index C.size_t) error {
	if status == C.MYB_OK {
		return nil
	}
	return &Error{Status: Status(status), Op: op, Index: int(index)}
}

// failDiag is fail with what a parser recorded.
func failDiag(op string, status C.myb_status, diag *C.myb_diag) error {
	if status == C.MYB_OK {
		return nil
	}
	return &Error{Status: Status(status), Op: op, Index: -1, Line: int(diag.line),
		Message: C.GoString(&diag.message[0])}
}

// refuse is the error for arguments that are not passed on to the library.
func refuse(op string, index int, message string) error {
	return &Error{Status: Domain, Op: op, Index: index, Message: message}
}
