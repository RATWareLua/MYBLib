// SPDX-License-Identifier: BSD-3-Clause

package myb

import (
	"runtime"
	"sync"
	"unsafe"
)

// object owns one object of the library. It is destroyed by close, or — as a safety net
// only — by a cleanup when the Go value holding it becomes unreachable. A call into C holds
// it (acquire ... release), and so does every object built on it for its whole life (an
// ephemeris used by a model): close then makes it unusable at once and destroys it when the
// last holder lets go.
type object struct {
	mu      sync.Mutex
	ptr     unsafe.Pointer
	users   int
	closed  bool
	destroy func(unsafe.Pointer)
	deps    []*object // held by this object until it is destroyed
	cleanup runtime.Cleanup
}

// teardown is what the cleanup does: it holds no reference to the object itself.
type teardown struct {
	ptr     unsafe.Pointer
	destroy func(unsafe.Pointer)
	deps    []*object
}

func (r teardown) run() {
	r.destroy(r.ptr)
	for _, d := range r.deps {
		d.release()
	}
}

// newObject takes ptr, and the holds already acquired on deps.
func newObject(ptr unsafe.Pointer, destroy func(unsafe.Pointer), deps ...*object) *object {
	o := &object{ptr: ptr, destroy: destroy, deps: deps}
	o.cleanup = runtime.AddCleanup(o, teardown.run, teardown{ptr: ptr, destroy: destroy, deps: deps})
	return o
}

// acquire holds the object for a call; ErrClosed after close.
func (o *object) acquire() (unsafe.Pointer, error) {
	if o == nil {
		return nil, ErrClosed
	}
	o.mu.Lock()
	defer o.mu.Unlock()
	if o.closed {
		return nil, ErrClosed
	}
	o.users++
	return o.ptr, nil
}

func (o *object) release() {
	o.mu.Lock()
	o.users--
	ptr := o.ptr
	free := o.closed && o.users == 0 && ptr != nil
	if free {
		o.ptr = nil
	}
	o.mu.Unlock()
	if free {
		o.free(ptr)
	}
}

func (o *object) close() {
	if o == nil {
		return
	}
	o.mu.Lock()
	if o.closed {
		o.mu.Unlock()
		return
	}
	o.closed = true
	ptr := o.ptr
	free := o.users == 0 && ptr != nil
	if free {
		o.ptr = nil
	}
	o.mu.Unlock()
	if free {
		o.free(ptr)
	}
}

func (o *object) free(ptr unsafe.Pointer) {
	o.cleanup.Stop()
	teardown{ptr: ptr, destroy: o.destroy, deps: o.deps}.run()
	o.deps = nil
}

// acquireAll holds every object given (nil ones are skipped) or none of them.
func acquireAll(objects ...*object) ([]*object, error) {
	held := make([]*object, 0, len(objects))
	for _, o := range objects {
		if o == nil {
			continue
		}
		if _, err := o.acquire(); err != nil {
			releaseAll(held)
			return nil, err
		}
		held = append(held, o)
	}
	return held, nil
}

func releaseAll(objects []*object) {
	for _, o := range objects {
		o.release()
	}
}
