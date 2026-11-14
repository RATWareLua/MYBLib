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

// ---- Shapes ---------------------------------------------------------------------------------

// Polyhedron is a closed triangular mesh (myb_polyhedron): the surface of a contact, the shape
// of a body.
type Polyhedron struct {
	o *object
}

// NewPolyhedron is the mesh of the vertices (m) and the faces (three vertex indices each,
// counter-clockwise seen from outside), of density so that G rho = gRho (s^-2)
// (myb_polyhedron_create).
func NewPolyhedron(vertices [][3]float64, faces [][3]int, gRho float64) (*Polyhedron, error) {
	cfaces := make([]C.size_t, 0, 3*len(faces))
	for k, f := range faces {
		for _, v := range f {
			if v < 0 || v >= len(vertices) {
				return nil, refuse("NewPolyhedron", k, "a vertex index out of range")
			}
			cfaces = append(cfaces, C.size_t(v))
		}
	}
	var p *C.myb_polyhedron
	status := C.myb_polyhedron_create((*C.double)(first(vertices)), C.size_t(len(vertices)),
		(*C.size_t)(first(cfaces)), C.size_t(len(faces)), C.double(gRho), nil, &p)
	if err := fail("NewPolyhedron", status); err != nil {
		return nil, err
	}
	return &Polyhedron{o: newObject(unsafe.Pointer(p), func(q unsafe.Pointer) {
		C.myb_polyhedron_destroy((*C.myb_polyhedron)(q))
	})}, nil
}

func (p *Polyhedron) obj() *object {
	if p == nil {
		return nil
	}
	return p.o
}

// Close releases the mesh (after the systems that use it).
func (p *Polyhedron) Close() error {
	p.obj().close()
	return nil
}

// ---- Systems of bodies ------------------------------------------------------------------------

// Interaction is the level at which a source acts on a receiver (myb_interaction).
type Interaction int

// The levels of the point masses.
const (
	NoInteraction Interaction = C.MYB_INTERACTION_NONE
	Newton        Interaction = C.MYB_INTERACTION_NEWTON
	EIH           Interaction = C.MYB_INTERACTION_EIH // with the 1PN terms
)

// FMM asks for the fast multipole method for the Newtonian sums (myb_nbody_set_approximation).
type FMM struct {
	Order int
	Theta float64
}

// SystemRotation is an integrated rotation (myb_system_rotation): moments per unit mass (m^2),
// the distortion (radius, k2, lag, mean motion) and a fluid core.
type SystemRotation struct {
	Moments      [3]float64
	Radius       float64
	K2           float64
	Lag          float64
	MeanMotion   float64
	Core         [2]float64
	CoreFriction float64
}

// SystemTide is the answer of a body to tides by band (myb_system_tide).
type SystemTide struct {
	Radius       float64
	K2           [3]float64
	OrbitalLag   [3]float64
	RotationLag  [3]float64
	RotationRate float64
}

// MassRate is a mass that changes: gm (1 + Rate (t - Epoch)) (myb_system_set_mass_rate).
type MassRate struct {
	Rate, Epoch float64
}

// SystemBody is a body of a system with its model; everything but GM is optional.
type SystemBody struct {
	GM          float64
	Radius      float64           // m: a sphere that may touch others; 0: none
	Orientation OrientationSource // the frame given from outside
	Field       *GravityField     // the figure, up to FieldDegree (0: all of it)
	FieldDegree int
	Rotation    *SystemRotation // an integrated rotation
	Tide        *SystemTide
	MassRate    *MassRate
}

// Link is a level of the graph between bodies of a system (indices in SystemSpec.Bodies).
type Link interface {
	link(l *C.mybgo_system_link, k *keeper, points *[]C.myb_contact_point) error
}

// InteractionLink sets the level at which Source acts on Receiver.
type InteractionLink struct {
	Receiver, Source int
	Level            Interaction
}

// FieldLink: the figure of Body acts on Point and back.
type FieldLink struct {
	Body, Point int
}

// FiguresLink: the quadrupoles of A and B act on each other.
type FiguresLink struct {
	A, B int
}

// TideLink: the tide that Raiser raises on Body acts on Point.
type TideLink struct {
	Body, Raiser, Point int
}

// DistortionLink: Raiser distorts Body, whose rotation is integrated.
type DistortionLink struct {
	Body, Raiser int
}

// PNOrders are the terms beyond 1PN of a pair (MYB_PN_* of myb/rel.h).
type PNOrders uint

// The orders of the pairs beyond 1PN.
const (
	PN1    PNOrders = C.MYB_PN_1
	PN2    PNOrders = C.MYB_PN_2
	PN25   PNOrders = C.MYB_PN_25
	PN3    PNOrders = C.MYB_PN_3
	PN35   PNOrders = C.MYB_PN_35
	PNSO   PNOrders = C.MYB_PN_SO
	PNSS   PNOrders = C.MYB_PN_SS
	PN4    PNOrders = C.MYB_PN_4
	PN4RR  PNOrders = C.MYB_PN_4RR
	PNSO25 PNOrders = C.MYB_PN_SO25
	PNSO35 PNOrders = C.MYB_PN_SO35
	PNSORR PNOrders = C.MYB_PN_SORR
	// spin-spin at next-to-leading order (1/c^6), with its precession (needs the kappas)
	PNSS3 PNOrders = C.MYB_PN_SS3
	// spin-orbit radiation reaction at next-to-leading order (1/c^10)
	PNSORR5 PNOrders = C.MYB_PN_SORR5
	// the spin-orbit part of the dissipative tail, quasi-circular orbits
	PNSOTAIL PNOrders = C.MYB_PN_SOTAIL
	// the tails of PN4 and PN4RR on any osculating conic (unbound, parabolic, e above 0.83)
	PNTAILX PNOrders = C.MYB_PN_TAILX
)

// PNPair adds the terms beyond 1PN of two bodies (myb_system_add_pn_pair).
type PNPair struct {
	A, B           int
	Orders         PNOrders
	C, R0          float64
	Kappa1, Kappa2 float64
}

// FrameDragging: the spin of Body (m^2/s, axes of the system) drags Point.
type FrameDragging struct {
	Body, Point int
	Spin        [3]float64
}

// ContactShape is the kind of a surface of contact.
type ContactShape int

// The surfaces of contact.
const (
	ContactPlane     ContactShape = C.MYB_CONTACT_PLANE
	ContactEllipsoid ContactShape = C.MYB_CONTACT_ELLIPSOID
	ContactMesh      ContactShape = C.MYB_CONTACT_MESH
)

// ContactSurface is the surface of a ground in its frame (myb_contact_surface).
type ContactSurface struct {
	Shape  ContactShape
	Point  [3]float64
	Normal [3]float64
	Radii  [3]float64
	Mesh   *Polyhedron
}

// ContactFriction is the friction of a law of contact (myb_contact_friction).
type ContactFriction int

// The frictions.
const (
	// Coulomb's made smooth, a function of the state: Slip
	FrictionSmooth ContactFriction = C.MYB_CONTACT_FRICTION_SMOOTH
	// a spring with a slip limit, a state of its own (3 values a point, ContactDeflections):
	// TangentialStiffness, TangentialDamping
	FrictionStick ContactFriction = C.MYB_CONTACT_FRICTION_STICK
)

// ContactLaw is the law of the force of contact (myb_contact_law): k, n, c of the normal
// springs, mu of the friction, and the friction's kind with its parameters (v_s of the smooth
// one, k_t and c_t of the one that sticks).
type ContactLaw struct {
	Stiffness, Exponent, Damping, Friction, Slip float64
	Kind                                         ContactFriction
	TangentialStiffness, TangentialDamping       float64
}

func (l ContactLaw) c() C.myb_contact_law {
	return C.myb_contact_law{stiffness: C.double(l.Stiffness), exponent: C.double(l.Exponent),
		damping: C.double(l.Damping), friction: C.double(l.Friction), slip: C.double(l.Slip),
		kind: C.myb_contact_friction(l.Kind), tangential_stiffness: C.double(l.TangentialStiffness),
		tangential_damping: C.double(l.TangentialDamping)}
}

// c is the surface of the library; a mesh is held by k.
func (s ContactSurface) c(k *keeper) (C.myb_contact_surface, error) {
	out := C.myb_contact_surface{shape: C.myb_contact_shape(s.Shape), point: cVec3(s.Point),
		normal: cVec3(s.Normal), radii: cVec3(s.Radii)}
	if s.Shape == ContactMesh {
		mesh, err := k.hold(s.Mesh.obj())
		if err != nil {
			return out, err
		}
		out.mesh = (*C.myb_polyhedron)(mesh)
	}
	return out, nil
}

func cContactPoints(points []ContactPoint) []C.myb_contact_point {
	out := make([]C.myb_contact_point, len(points))
	for k, q := range points {
		out[k] = C.myb_contact_point{offset: cVec3(q.Offset), radius: C.double(q.Radius)}
	}
	return out
}

// ContactPoint is a contact point of a body in its frame (myb_contact_point).
type ContactPoint struct {
	Offset [3]float64
	Radius float64
}

// ContactLink: Body touches the surface of Ground with its points.
type ContactLink struct {
	Body, Ground int
	Surface      ContactSurface
	Law          ContactLaw
	Points       []ContactPoint
}

// SystemView is what an action sees: the body acted on, the positions, velocities and GM of
// all bodies now, and the frame of the body (rows, nil without one). The slices are valid
// during the call only.
type SystemView struct {
	Body  int
	X, V  [][3]float64
	GM    []float64
	Frame *[9]float64
}

// ActionFunc is the acceleration (m/s^2) and the torque per unit mass (m^2/s^2, body frame) of
// an action at the time t.
type ActionFunc func(t float64, view SystemView) ([3]float64, [3]float64, error)

// ActionLink: a Go function acts on Body (the slow path of callbacks).
type ActionLink struct {
	Body int
	Func ActionFunc
}

// Guided is a body whose motion is given (myb_system_add_guided).
type Guided struct {
	GM    float64
	State StateSource
}

// SystemSpec describes a system of bodies with their models (myb/system.h).
type SystemSpec struct {
	Bodies        []SystemBody
	Relativity    *PPN // the 1PN terms of the point masses
	Approximation *FMM
	Links         []Link
	Guided        []Guided
}

func (l InteractionLink) link(c *C.mybgo_system_link, k *keeper, p *[]C.myb_contact_point) error {
	c.kind, c.a, c.b, c.level = C.MYBGO_LINK_INTERACTION, C.size_t(l.Receiver), C.size_t(l.Source), C.int(l.Level)
	return nil
}

func (l FieldLink) link(c *C.mybgo_system_link, k *keeper, p *[]C.myb_contact_point) error {
	c.kind, c.a, c.b = C.MYBGO_LINK_FIELD, C.size_t(l.Body), C.size_t(l.Point)
	return nil
}

func (l FiguresLink) link(c *C.mybgo_system_link, k *keeper, p *[]C.myb_contact_point) error {
	c.kind, c.a, c.b = C.MYBGO_LINK_FIGURES, C.size_t(l.A), C.size_t(l.B)
	return nil
}

func (l TideLink) link(c *C.mybgo_system_link, k *keeper, p *[]C.myb_contact_point) error {
	c.kind, c.a, c.b, c.c = C.MYBGO_LINK_TIDE, C.size_t(l.Body), C.size_t(l.Raiser), C.size_t(l.Point)
	return nil
}

func (l DistortionLink) link(c *C.mybgo_system_link, k *keeper, p *[]C.myb_contact_point) error {
	c.kind, c.a, c.b = C.MYBGO_LINK_DISTORTION, C.size_t(l.Body), C.size_t(l.Raiser)
	return nil
}

func (l PNPair) link(c *C.mybgo_system_link, k *keeper, p *[]C.myb_contact_point) error {
	c.kind, c.a, c.b = C.MYBGO_LINK_PN_PAIR, C.size_t(l.A), C.size_t(l.B)
	c.pn = C.myb_system_pn_pair{orders: C.uint(l.Orders), c: C.double(l.C), r0: C.double(l.R0),
		kappa1: C.double(l.Kappa1), kappa2: C.double(l.Kappa2)}
	return nil
}

func (l FrameDragging) link(c *C.mybgo_system_link, k *keeper, p *[]C.myb_contact_point) error {
	c.kind, c.a, c.b = C.MYBGO_LINK_DRAGGING, C.size_t(l.Body), C.size_t(l.Point)
	c.spin = [3]C.double{C.double(l.Spin[0]), C.double(l.Spin[1]), C.double(l.Spin[2])}
	return nil
}

func cVec3(v [3]float64) [3]C.double {
	return [3]C.double{C.double(v[0]), C.double(v[1]), C.double(v[2])}
}

func (l ContactLink) link(c *C.mybgo_system_link, k *keeper, p *[]C.myb_contact_point) error {
	c.kind, c.a, c.b = C.MYBGO_LINK_CONTACT, C.size_t(l.Body), C.size_t(l.Ground)
	surface, err := l.Surface.c(k)
	if err != nil {
		return err
	}
	c.surface = surface
	c.law = l.Law.c()
	c.first_point = C.size_t(len(*p))
	c.point_count = C.size_t(len(l.Points))
	*p = append(*p, cContactPoints(l.Points)...)
	return nil
}

func (l ActionLink) link(c *C.mybgo_system_link, k *keeper, p *[]C.myb_contact_point) error {
	if l.Func == nil {
		return refuse("NewSystem", -1, "an action without a function")
	}
	c.kind, c.a = C.MYBGO_LINK_ACTION, C.size_t(l.Body)
	c.handle = k.handle(&actionCallback{f: l.Func, state: k.state})
	return nil
}

// System is a system of bodies of the library (myb_system) made from a SystemSpec. Runs change
// it (mergers); it serves one call at a time.
type System struct {
	o      *object
	mu     sync.Mutex
	state  *callState
	bodies int
	points int // the most points of a contact
}

// NewSystem makes the system of a description (one call into C).
func NewSystem(spec SystemSpec) (*System, error) {
	const op = "NewSystem"
	k := newKeeper()
	ok := false
	defer func() {
		if !ok {
			k.done()
		}
	}()
	var options C.mybgo_system_options
	if spec.Relativity != nil {
		options.has_relativity = 1
		options.relativity = cPPN(*spec.Relativity)
	}
	if spec.Approximation != nil {
		options.fmm_order = C.int(spec.Approximation.Order)
		options.fmm_theta = C.double(spec.Approximation.Theta)
	}
	bodies := make([]C.mybgo_system_body, len(spec.Bodies))
	for i, b := range spec.Bodies {
		d := &bodies[i]
		d.gm = C.double(b.GM)
		d.radius = C.double(b.Radius)
		if b.Orientation != nil {
			if err := b.Orientation.describe(&d.orientation, k); err != nil {
				return nil, err
			}
		}
		if b.Field != nil {
			field, err := k.hold(b.Field.obj())
			if err != nil {
				return nil, err
			}
			d.field = (*C.myb_gravity_field)(field)
			d.degree = C.int(b.FieldDegree)
		}
		if r := b.Rotation; r != nil {
			d.has_rotation = 1
			d.rotation = C.myb_system_rotation{moments: cVec3(r.Moments), radius: C.double(r.Radius),
				k2: C.double(r.K2), lag: C.double(r.Lag), mean_motion: C.double(r.MeanMotion),
				core:          [2]C.double{C.double(r.Core[0]), C.double(r.Core[1])},
				core_friction: C.double(r.CoreFriction)}
		}
		if t := b.Tide; t != nil {
			d.has_tide = 1
			d.tide = C.myb_system_tide{radius: C.double(t.Radius), k2: cVec3(t.K2),
				orbital_lag: cVec3(t.OrbitalLag), rotation_lag: cVec3(t.RotationLag),
				rotation_rate: C.double(t.RotationRate)}
		}
		if m := b.MassRate; m != nil {
			d.has_mass_rate = 1
			d.mass_rate = C.double(m.Rate)
			d.mass_epoch = C.double(m.Epoch)
		}
	}
	links := make([]C.mybgo_system_link, len(spec.Links))
	var points []C.myb_contact_point
	most := 0
	for j, l := range spec.Links {
		if l == nil {
			return nil, refuse(op, j, "a nil link")
		}
		if err := l.link(&links[j], k, &points); err != nil {
			return nil, err
		}
		if links[j].kind == C.MYBGO_LINK_CONTACT {
			most = max(most, int(links[j].point_count))
		}
	}
	guided := make([]C.mybgo_provider, len(spec.Guided))
	guidedGM := make([]float64, len(spec.Guided))
	for g, b := range spec.Guided {
		if b.State == nil {
			return nil, refuse(op, g, "a guided body without a provider")
		}
		if err := b.State.provider(&guided[g], k); err != nil {
			return nil, err
		}
		guidedGM[g] = b.GM
	}
	var system *C.mybgo_system
	var stage C.int
	var item C.size_t
	status := C.mybgo_system_create(&options, (*C.mybgo_system_body)(first(bodies)),
		C.size_t(len(bodies)), (*C.mybgo_system_link)(first(links)), C.size_t(len(links)),
		(*C.myb_contact_point)(first(points)), C.size_t(len(points)),
		(*C.mybgo_provider)(first(guided)), cDoubles(guidedGM), C.size_t(len(guided)),
		&system, &stage, &item)
	if status != C.MYB_OK {
		return nil, stageError(op, status, stage, item, "body", "link", "guided body")
	}
	held, handles := k.keep()
	ok = true
	s := &System{state: k.state, bodies: len(spec.Bodies), points: most}
	s.o = newObject(unsafe.Pointer(system), func(q unsafe.Pointer) {
		C.mybgo_system_destroy((*C.mybgo_system)(q))
		for _, h := range handles {
			h.Delete()
		}
	}, held...)
	return s, nil
}

// Close releases the system.
func (s *System) Close() error {
	if s != nil {
		s.o.close()
	}
	return nil
}

func (s *System) begin() (*C.mybgo_system, error) {
	if s == nil {
		return nil, ErrClosed
	}
	s.mu.Lock()
	ptr, err := s.o.acquire()
	if err != nil {
		s.mu.Unlock()
		return nil, err
	}
	s.state.err = nil
	return (*C.mybgo_system)(ptr), nil
}

func (s *System) end() {
	s.o.release()
	s.mu.Unlock()
}

func (s *System) failure(op string, status C.myb_status, index int) error {
	if status == C.MYB_OK {
		return nil
	}
	return &Error{Status: Status(status), Op: op, Index: index, Err: s.state.err}
}

// Coordinates is n, the count of the coordinates of a state (3 a body, 9 an integrated
// rotation, 6 a pair with spins, 3 a point of a contact whose friction sticks); a state holds
// 2 n values.
func (s *System) Coordinates() (int, error) {
	sys, err := s.begin()
	if err != nil {
		return 0, err
	}
	defer s.end()
	return int(C.mybgo_system_coordinates(sys)), nil
}

// Present tells whether a body is still there (not merged into another), and its GM now.
func (s *System) Present(body int) (bool, float64, error) {
	sys, err := s.begin()
	if err != nil {
		return false, 0, err
	}
	defer s.end()
	if body < 0 || body >= s.bodies {
		return false, 0, refuse("System.Present", -1, "no such body")
	}
	return C.mybgo_system_present(sys, C.size_t(body)) != 0,
		float64(C.mybgo_system_gm(sys, C.size_t(body))), nil
}

// RotationState is the integrated rotation of a body in a state: the attitude quaternion,
// the angular velocity of the mantle and of its core (rad/s, mantle frame).
type RotationState struct {
	Body      int
	Q         [4]float64
	Omega     [3]float64
	OmegaCore *[3]float64 // nil: no core
}

// SpinState are the two dimensionless spins of a PN pair (by its order among the pairs).
type SpinState struct {
	Pair       int
	Chi1, Chi2 [3]float64
}

func cRotationStates(r []RotationState) []C.mybgo_rotation_state {
	out := make([]C.mybgo_rotation_state, len(r))
	for k, v := range r {
		out[k].body = C.size_t(v.Body)
		for i := range 4 {
			out[k].q[i] = C.double(v.Q[i])
		}
		out[k].omega = cVec3(v.Omega)
		if v.OmegaCore != nil {
			out[k].has_core = 1
			out[k].omega_core = cVec3(*v.OmegaCore)
		}
	}
	return out
}

func cSpinStates(s []SpinState) []C.mybgo_spin_state {
	out := make([]C.mybgo_spin_state, len(s))
	for k, v := range s {
		out[k] = C.mybgo_spin_state{pair: C.size_t(v.Pair), chi1: cVec3(v.Chi1), chi2: cVec3(v.Chi2)}
	}
	return out
}

// PackState is the state y (2 n values) of the bodies' states, the rotations and the spins.
func (s *System) PackState(bodies []State, rotations []RotationState, spins []SpinState) ([]float64, error) {
	const op = "System.PackState"
	sys, err := s.begin()
	if err != nil {
		return nil, err
	}
	defer s.end()
	if len(bodies) != s.bodies {
		return nil, refuse(op, -1, "not one state for every body")
	}
	y := make([]float64, 2*int(C.mybgo_system_coordinates(sys)))
	cr, cs := cRotationStates(rotations), cSpinStates(spins)
	status := C.mybgo_system_pack(sys, cStates(bodies), (*C.mybgo_rotation_state)(first(cr)),
		C.size_t(len(cr)), (*C.mybgo_spin_state)(first(cs)), C.size_t(len(cs)), cDoubles(y))
	if err := fail(op, status); err != nil {
		return nil, err
	}
	return y, nil
}

// UnpackState reads the bodies' states of y, and the rotations and spins named by the records
// given (their Body and Pair, and whether a core is there).
func (s *System) UnpackState(y []float64, rotations []RotationState, spins []SpinState) ([]State, []RotationState, []SpinState, error) {
	const op = "System.UnpackState"
	sys, err := s.begin()
	if err != nil {
		return nil, nil, nil, err
	}
	defer s.end()
	if len(y) != 2*int(C.mybgo_system_coordinates(sys)) {
		return nil, nil, nil, refuse(op, -1, "a state of another length")
	}
	bodies := make([]State, s.bodies)
	cr, cs := cRotationStates(rotations), cSpinStates(spins)
	status := C.mybgo_system_unpack(sys, cDoubles(y), cStates(bodies),
		(*C.mybgo_rotation_state)(first(cr)), C.size_t(len(cr)), (*C.mybgo_spin_state)(first(cs)),
		C.size_t(len(cs)))
	if err := fail(op, status); err != nil {
		return nil, nil, nil, err
	}
	outR := make([]RotationState, len(cr))
	for k, v := range cr {
		outR[k] = RotationState{Body: int(v.body)}
		for i := range 4 {
			outR[k].Q[i] = float64(v.q[i])
		}
		for i := range 3 {
			outR[k].Omega[i] = float64(v.omega[i])
		}
		if v.has_core != 0 {
			core := [3]float64{float64(v.omega_core[0]), float64(v.omega_core[1]), float64(v.omega_core[2])}
			outR[k].OmegaCore = &core
		}
	}
	outS := make([]SpinState, len(cs))
	for k, v := range cs {
		outS[k] = SpinState{Pair: int(v.pair)}
		for i := range 3 {
			outS[k].Chi1[i] = float64(v.chi1[i])
			outS[k].Chi2[i] = float64(v.chi2[i])
		}
	}
	return bodies, outR, outS, nil
}

// ContactDeflections are the deflections of the springs of a friction that sticks in a state
// (myb_system_pack_contact): Contact is the number of the contact among those in place (the
// ContactLinks in their order, less those of a body a merger took away), Deflections one a
// point of it (m, in the axes of its surface; zero: the springs at rest).
type ContactDeflections struct {
	Contact     int
	Deflections [][3]float64
}

// PackContacts writes the deflections of contacts into the state y (of PackState), in place.
// The deflections of a contact not given stay as they are in y (zero after PackState).
func (s *System) PackContacts(y []float64, contacts []ContactDeflections) error {
	const op = "System.PackContacts"
	sys, err := s.begin()
	if err != nil {
		return err
	}
	defer s.end()
	if len(y) != 2*int(C.mybgo_system_coordinates(sys)) {
		return refuse(op, -1, "a state of another length")
	}
	numbers := make([]C.size_t, len(contacts))
	counts := make([]C.size_t, len(contacts))
	var values [][3]float64
	for k, c := range contacts {
		if c.Contact < 0 {
			return refuse(op, k, "no such contact")
		}
		numbers[k], counts[k] = C.size_t(c.Contact), C.size_t(len(c.Deflections))
		values = append(values, c.Deflections...)
	}
	var failed C.size_t
	status := C.mybgo_system_contacts(sys, 1, cDoubles(y), (*C.size_t)(first(numbers)),
		(*C.size_t)(first(counts)), C.size_t(len(contacts)), (*C.double)(first(values)),
		C.size_t(3*len(values)), &failed)
	return failAt(op, status, failed)
}

// UnpackContacts reads the deflections of the contacts numbered (among those in place) from the
// state y.
func (s *System) UnpackContacts(y []float64, contacts []int) ([]ContactDeflections, error) {
	const op = "System.UnpackContacts"
	sys, err := s.begin()
	if err != nil {
		return nil, err
	}
	defer s.end()
	if len(y) != 2*int(C.mybgo_system_coordinates(sys)) {
		return nil, refuse(op, -1, "a state of another length")
	}
	numbers := make([]C.size_t, len(contacts))
	for k, c := range contacts {
		if c < 0 {
			return nil, refuse(op, k, "no such contact")
		}
		numbers[k] = C.size_t(c)
	}
	counts := make([]C.size_t, len(contacts))
	values := make([][3]float64, s.points*len(contacts))
	var failed C.size_t
	status := C.mybgo_system_contacts(sys, 0, cDoubles(y), (*C.size_t)(first(numbers)),
		(*C.size_t)(first(counts)), C.size_t(len(contacts)), (*C.double)(first(values)),
		C.size_t(3*len(values)), &failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	out := make([]ContactDeflections, len(contacts))
	at := 0
	for k, c := range contacts {
		n := int(counts[k])
		out[k] = ContactDeflections{Contact: c, Deflections: values[at : at+n : at+n]}
		at += n
	}
	return out, nil
}

// IAS15Options are the options of the IAS15 integrator of a system (myb_ias15_options: the
// precision, the step limits, the budget, the timescale of the step).
type IAS15Options struct {
	Epsilon, HInitial, HMax, HMin float64
	MaxSteps                      int64
	// TimescaleJerk for an acceleration that passes through zero: a body held by the springs
	// of a contact
	Timescale IAS15Timescale
}

func (o IAS15Options) c() C.myb_ias15_options {
	return C.myb_ias15_options{epsilon: C.double(o.Epsilon), h_initial: C.double(o.HInitial),
		h_max: C.double(o.HMax), h_min: C.double(o.HMin), max_steps: C.int64_t(o.MaxSteps),
		timescale: C.myb_ias15_timescale(o.Timescale)}
}

// ---- The force of a contact -------------------------------------------------------------------

// ContactSample is a state of a body over a surface (myb_contact_force): X and V of its centre
// of mass in the frame of the surface (m, m/s), Attitude the rotation from the axes of the
// surface to the body axes (rows: the body axes; nil: the same axes), Omega its angular
// velocity relative to the surface in the body axes (rad/s; nil: none), and with a friction
// that sticks the Deflections of its springs, one a point (myb_contact_force_stick).
type ContactSample struct {
	X, V        [3]float64
	Attitude    *[9]float64
	Omega       *[3]float64
	Deflections [][3]float64
}

// ContactResult is the force of a surface on a body (myb_contact_result): the force per unit
// mass in the axes of the surface (m/s^2), the torque per unit mass in the body axes
// (m^2/s^2), the energy of the springs (m^2/s^2), the least height of a sphere over the
// surface less its radius (m), the points that press in, and with a friction that sticks the
// rates of the deflections (m/s), one a point.
type ContactResult struct {
	Force, Torque   [3]float64
	Energy, Gap     float64
	Touching        int
	DeflectionRates [][3]float64
}

// ContactForces is the force of the surface on the points of a body (in its frame) at each
// sample, in one call (myb_contact_force_stick, which with the smooth friction is
// myb_contact_force).
func ContactForces(surface ContactSurface, law ContactLaw, points []ContactPoint, samples []ContactSample) ([]ContactResult, error) {
	const op = "ContactForces"
	sticks := law.Kind == FrictionStick
	cs := make([]C.mybgo_contact_sample, len(samples))
	var deflections []float64
	for k, q := range samples {
		cs[k].x, cs[k].v = cVec3(q.X), cVec3(q.V)
		if q.Attitude != nil {
			cs[k].has_attitude = 1
			for i := range 9 {
				cs[k].attitude[i] = C.double(q.Attitude[i])
			}
		}
		if q.Omega != nil {
			cs[k].has_omega = 1
			cs[k].omega = cVec3(*q.Omega)
		}
		if sticks {
			if len(q.Deflections) != len(points) {
				return nil, refuse(op, k, "not one deflection for every point")
			}
			for _, z := range q.Deflections {
				deflections = append(deflections, z[0], z[1], z[2])
			}
		}
	}
	k := newKeeper()
	defer k.done()
	cSurface, err := surface.c(k)
	if err != nil {
		return nil, err
	}
	cLaw := law.c()
	cPoints := cContactPoints(points)
	var rates []float64
	if sticks {
		rates = make([]float64, len(deflections))
	}
	results := make([]C.myb_contact_result, len(samples))
	var failed C.size_t
	status := C.mybgo_contact_forces(&cSurface, &cLaw, (*C.myb_contact_point)(first(cPoints)),
		C.size_t(len(cPoints)), (*C.mybgo_contact_sample)(first(cs)), C.size_t(len(cs)),
		cDoubles(deflections), cDoubles(rates), (*C.myb_contact_result)(first(results)), &failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	out := make([]ContactResult, len(results))
	for i, r := range results {
		out[i] = ContactResult{Energy: float64(r.energy), Gap: float64(r.gap), Touching: int(r.touching)}
		for j := range 3 {
			out[i].Force[j] = float64(r.force[j])
			out[i].Torque[j] = float64(r.torque[j])
		}
		if sticks {
			out[i].DeflectionRates = make([][3]float64, len(points))
			for p := range points {
				at := 3 * (len(points)*i + p)
				out[i].DeflectionRates[p] = [3]float64{rates[at], rates[at+1], rates[at+2]}
			}
		}
	}
	return out, nil
}

// ContactPolicy is what a run does at a contact of two bodies.
type ContactPolicy int

// The policies of a run.
const (
	StopAtContact     ContactPolicy = C.MYBGO_CONTACT_STOP
	MergeAtContact    ContactPolicy = C.MYBGO_CONTACT_MERGE  // the lighter into the heavier
	BounceAtContact   ContactPolicy = C.MYBGO_CONTACT_BOUNCE // with Restitution
	ContinueAtContact ContactPolicy = C.MYBGO_CONTACT_CONTINUE
)

// SystemEventKind is the kind of an event of a run (myb_system_event_kind).
type SystemEventKind int

// The events of a run.
const (
	Contact  SystemEventKind = C.MYB_SYSTEM_CONTACT
	Approach SystemEventKind = C.MYB_SYSTEM_APPROACH
)

// Merger is what a merger takes out of the motion (myb_system_merger, G as a factor).
type Merger struct {
	Energy          float64
	AngularMomentum [3]float64
	Spin            [3]float64
}

// SystemEvent is an event of a run: a contact acted on by the policy, or an approach.
type SystemEvent struct {
	Kind     SystemEventKind
	A, B     int
	T        float64
	Distance float64
	Action   ContactPolicy // what was done
	Survivor int           // MergeAtContact: the body that remains
	Merger   Merger
}

// SystemRun asks for the states of a system from Y0 at T0 at Times (in the order of the
// integration): the integrator ends a step at every time asked and at every event. Watch > 0
// reports the closest approaches within that distance. Limit (0: 1024) bounds the events.
type SystemRun struct {
	Options     IAS15Options
	T0          float64
	Y0          []float64
	Times       []float64
	Watch       float64
	OnContact   ContactPolicy
	Restitution float64
	Limit       int
}

// IAS15Stats are the counts of an IAS15 integration.
type IAS15Stats struct {
	Accepted, Rejected, Evaluations, Iterations, Unconverged int64
}

// SystemResult is what a run gives: the states at the times reached, the events, where it
// stopped.
type SystemResult struct {
	States [][]float64
	Events []SystemEvent
	T      float64
	Final  []float64
	Stats  IAS15Stats
}

// Run integrates the system (one call into C); mergers change the system for later calls.
func (s *System) Run(r SystemRun) (*SystemResult, error) {
	const op = "System.Run"
	sys, err := s.begin()
	if err != nil {
		return nil, err
	}
	defer s.end()
	n := 2 * int(C.mybgo_system_coordinates(sys))
	if len(r.Y0) != n {
		return nil, refuse(op, -1, "a first state of another length")
	}
	limit := r.Limit
	if limit <= 0 {
		limit = 1024
	}
	options := C.mybgo_system_run_options{ias15: r.Options.c(), watch: C.double(r.Watch),
		policy: C.int(r.OnContact), restitution: C.double(r.Restitution), limit: C.size_t(limit)}
	flat := make([]float64, n*len(r.Times))
	final := make([]float64, n)
	var events *C.mybgo_system_event
	var count, failed C.size_t
	var tEnd C.double
	var stats C.myb_ias15_stats
	status := C.mybgo_system_run(sys, &options, C.double(r.T0), cDoubles(r.Y0), cDoubles(r.Times),
		C.size_t(len(r.Times)), cDoubles(flat), &events, &count, &tEnd, cDoubles(final), &stats,
		&failed)
	defer C.mybgo_free(unsafe.Pointer(events))
	out := &SystemResult{T: float64(tEnd), Final: final, States: make([][]float64, int(failed)),
		Stats: IAS15Stats{Accepted: int64(stats.accepted), Rejected: int64(stats.rejected),
			Evaluations: int64(stats.evaluations), Iterations: int64(stats.iterations),
			Unconverged: int64(stats.unconverged)}}
	for k := range out.States {
		out.States[k] = flat[n*k : n*(k+1) : n*(k+1)]
	}
	for _, e := range unsafe.Slice(events, int(count)) {
		out.Events = append(out.Events, SystemEvent{Kind: SystemEventKind(e.kind), A: int(e.a),
			B: int(e.b), T: float64(e.t), Distance: float64(e.distance), Action: ContactPolicy(e.action),
			Survivor: int(e.survivor), Merger: Merger{Energy: float64(e.merger.energy),
				AngularMomentum: goVec3(e.merger.angular_momentum), Spin: goVec3(e.merger.spin)}})
	}
	return out, s.failure(op, status, int(failed))
}

func goVec3(v [3]C.double) [3]float64 {
	return [3]float64{float64(v[0]), float64(v[1]), float64(v[2])}
}

// flatStates checks that every state has 2 n values and lays them end to end.
func flatStates(op string, y [][]float64, n int) ([]float64, error) {
	flat := make([]float64, 0, n*len(y))
	for k, v := range y {
		if len(v) != n {
			return nil, refuse(op, k, "a state of another length")
		}
		flat = append(flat, v...)
	}
	return flat, nil
}

// Accelerations are the second derivatives of the coordinates (n each) at the times and
// states, in one call (myb_system_acceleration).
func (s *System) Accelerations(t []float64, y [][]float64) ([][]float64, error) {
	const op = "System.Accelerations"
	sys, err := s.begin()
	if err != nil {
		return nil, err
	}
	defer s.end()
	n := int(C.mybgo_system_coordinates(sys))
	if len(t) != len(y) {
		return nil, refuse(op, -1, "times and states of different lengths")
	}
	flat, err := flatStates(op, y, 2*n)
	if err != nil {
		return nil, err
	}
	a := make([]float64, n*len(y))
	var failed C.size_t
	status := C.mybgo_system_accelerations(sys, cDoubles(t), cDoubles(flat), C.size_t(len(y)),
		cDoubles(a), &failed)
	if err := s.failure(op, status, int(failed)); err != nil {
		return nil, err
	}
	out := make([][]float64, len(y))
	for k := range out {
		out[k] = a[n*k : n*(k+1) : n*(k+1)]
	}
	return out, nil
}

// DD is a double-word number: Hi + Lo.
type DD struct {
	Hi, Lo float64
}

// Invariants are the conserved quantities of the point masses (myb_nbody_invariants, G as a
// factor).
type Invariants struct {
	Energy          DD
	Momentum        [3]float64
	AngularMomentum [3]float64
	Center          [3]float64
	CenterVelocity  [3]float64
}

// PNInvariants are the 1PN conserved quantities (myb_nbody_pn_invariants).
type PNInvariants struct {
	Energy         DD
	Momentum       [3]float64
	Center         [3]float64
	CenterVelocity [3]float64
}

// Invariants of the point masses at every state, in one call.
func (s *System) Invariants(y [][]float64) ([]Invariants, error) {
	const op = "System.Invariants"
	sys, err := s.begin()
	if err != nil {
		return nil, err
	}
	defer s.end()
	flat, err := flatStates(op, y, 2*int(C.mybgo_system_coordinates(sys)))
	if err != nil {
		return nil, err
	}
	c := make([]C.myb_nbody_invariants, len(y))
	var failed C.size_t
	status := C.mybgo_system_invariants(sys, cDoubles(flat), C.size_t(len(y)),
		(*C.myb_nbody_invariants)(first(c)), &failed)
	if err := s.failure(op, status, int(failed)); err != nil {
		return nil, err
	}
	out := make([]Invariants, len(c))
	for k, v := range c {
		out[k] = Invariants{Energy: DD{Hi: float64(v.energy.hi), Lo: float64(v.energy.lo)},
			Momentum: goVec3(v.momentum), AngularMomentum: goVec3(v.angular_momentum),
			Center: goVec3(v.center), CenterVelocity: goVec3(v.center_velocity)}
	}
	return out, nil
}

// PNInvariants of the point masses at every state, in one call (needs the relativity).
func (s *System) PNInvariants(y [][]float64) ([]PNInvariants, error) {
	const op = "System.PNInvariants"
	sys, err := s.begin()
	if err != nil {
		return nil, err
	}
	defer s.end()
	flat, err := flatStates(op, y, 2*int(C.mybgo_system_coordinates(sys)))
	if err != nil {
		return nil, err
	}
	c := make([]C.myb_nbody_pn_invariants, len(y))
	var failed C.size_t
	status := C.mybgo_system_pn_invariants(sys, cDoubles(flat), C.size_t(len(y)),
		(*C.myb_nbody_pn_invariants)(first(c)), &failed)
	if err := s.failure(op, status, int(failed)); err != nil {
		return nil, err
	}
	out := make([]PNInvariants, len(c))
	for k, v := range c {
		out[k] = PNInvariants{Energy: DD{Hi: float64(v.energy.hi), Lo: float64(v.energy.lo)},
			Momentum: goVec3(v.momentum), Center: goVec3(v.center),
			CenterVelocity: goVec3(v.center_velocity)}
	}
	return out, nil
}

// SystemPartials are the derivatives of the state at T1: Phi (2 n x 2 n, row-major) and by
// the GM of the bodies asked (2 n x len(bodies)).
type SystemPartials struct {
	Y1          []float64
	Phi         []float64
	Sensitivity []float64
}

// Partials propagates the variational equations of the system from (t0, y0) to t1
// (myb_system_partials) with the derivatives by the GM of `bodies`.
func (s *System) Partials(options IAS15Options, t0 float64, y0 []float64, t1 float64, bodies []int) (*SystemPartials, error) {
	const op = "System.Partials"
	sys, err := s.begin()
	if err != nil {
		return nil, err
	}
	defer s.end()
	n := 2 * int(C.mybgo_system_coordinates(sys))
	if len(y0) != n {
		return nil, refuse(op, -1, "a first state of another length")
	}
	cb := make([]C.size_t, len(bodies))
	for k, b := range bodies {
		cb[k] = C.size_t(b)
	}
	out := &SystemPartials{Y1: make([]float64, n), Phi: make([]float64, n*n),
		Sensitivity: make([]float64, n*len(bodies))}
	co := options.c()
	status := C.mybgo_system_partials(sys, &co, C.double(t0), cDoubles(y0), C.double(t1),
		(*C.size_t)(first(cb)), C.size_t(len(cb)), cDoubles(out.Y1), cDoubles(out.Phi),
		cDoubles(out.Sensitivity))
	if err := s.failure(op, status, -1); err != nil {
		return nil, err
	}
	return out, nil
}

// ---- Wisdom-Holman integrations -----------------------------------------------------------

// WHScheme is a splitting scheme of the Wisdom-Holman family (myb_wh_scheme).
type WHScheme int

// The schemes.
const (
	WHDKD   WHScheme = C.MYB_WH_DKD
	WHSABA2 WHScheme = C.MYB_WH_SABA2
	WHSABA3 WHScheme = C.MYB_WH_SABA3
	WHSABA4 WHScheme = C.MYB_WH_SABA4
	WHSBAB1 WHScheme = C.MYB_WH_SBAB1
	WHSBAB2 WHScheme = C.MYB_WH_SBAB2
	WHSBAB3 WHScheme = C.MYB_WH_SBAB3
	WHSBAB4 WHScheme = C.MYB_WH_SBAB4
)

// WHRun asks for a Wisdom-Holman integration of point masses: GM of the massive bodies (the
// central one first), Test test bodies after them, the step, scheme and corrector, the
// positions and velocities at T0 (3 per body each) and the times.
type WHRun struct {
	GM        []float64
	Test      int
	H         float64
	Scheme    WHScheme
	Corrector int
	T0        float64
	X, V      [][3]float64
	Times     []float64
}

// WHStats are the counts of a Wisdom-Holman integration.
type WHStats struct {
	Steps, Evaluations, Drifts int64
}

// WHResult are the positions and velocities at the times.
type WHResult struct {
	X, V  [][][3]float64
	Stats WHStats
}

// RunWH integrates point masses by a Wisdom-Holman scheme (one call into C).
func RunWH(r WHRun) (*WHResult, error) {
	const op = "RunWH"
	bodies := len(r.GM) + r.Test
	if r.Test < 0 || len(r.X) != bodies || len(r.V) != bodies {
		return nil, refuse(op, -1, "not one position and velocity for every body")
	}
	options := C.myb_wh_options{h: C.double(r.H), scheme: C.myb_wh_scheme(r.Scheme),
		corrector: C.int(r.Corrector)}
	xs := make([][3]float64, bodies*len(r.Times))
	vs := make([][3]float64, bodies*len(r.Times))
	var stats C.myb_wh_stats
	var failed C.size_t
	status := C.mybgo_wh_run(cDoubles(r.GM), C.size_t(len(r.GM)), C.size_t(r.Test), &options,
		C.double(r.T0), (*C.double)(first(r.X)), (*C.double)(first(r.V)), cDoubles(r.Times),
		C.size_t(len(r.Times)), (*C.double)(first(xs)), (*C.double)(first(vs)), &stats, &failed)
	if err := failAt(op, status, failed); err != nil {
		return nil, err
	}
	out := &WHResult{Stats: WHStats{Steps: int64(stats.steps), Evaluations: int64(stats.evaluations),
		Drifts: int64(stats.drifts)}}
	for k := range r.Times {
		out.X = append(out.X, xs[bodies*k:bodies*(k+1):bodies*(k+1)])
		out.V = append(out.V, vs[bodies*k:bodies*(k+1):bodies*(k+1)])
	}
	return out, nil
}
