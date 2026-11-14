# SPDX-License-Identifier: BSD-3-Clause
"""Systems of bodies: the N-body problem of point masses with its graph of interactions, its
threads, the fast multipole approximation, its invariants and the Wisdom-Holman integrators
(myb/dyn.h); the system of bodies with their models — fields, tides, figures, rotations,
pairs beyond 1PN, actions, contact — and its encounters (myb/system.h).

Positions and velocities of many bodies come in as a sequence of triples (or flat, three
numbers a body) and go out as a list of triples. The state of a System is flat, as the C
interface lays it out: myb_system_coordinates(n) coordinates, then their rates (2 n
numbers); pack_rotation and its relatives fill and read its parts."""
from collections import namedtuple

from ._core import (Provider, _Calls, _check, _doubles, _flat, _optional, _orientation_provider,
                    _Owner, _rows, _settle, _state_provider, _tuple)
from ._native import ffi, lib

__all__ = [
    "NBody", "NBodyInvariants", "NBodyPNInvariants", "Pool", "System", "SystemEvent",
    "SystemIntegrator", "SystemMerger", "SystemView", "WisdomHolman", "WisdomHolmanStats",
    "Ias15Stats",
]

NBodyInvariants = namedtuple("NBodyInvariants",
                             "energy momentum angular_momentum center center_velocity")
NBodyPNInvariants = namedtuple("NBodyPNInvariants", "energy momentum center center_velocity")
WisdomHolmanStats = namedtuple("WisdomHolmanStats", "steps evaluations drifts")
Ias15Stats = namedtuple("Ias15Stats", "accepted rejected evaluations iterations unconverged")
SystemView = namedtuple("SystemView", "body bodies x v gm frame")
SystemEvent = namedtuple("SystemEvent", "kind a b t distance")
SystemMerger = namedtuple("SystemMerger", "energy angular_momentum spin")

_levels = {"none": lib.MYB_INTERACTION_NONE, "newton": lib.MYB_INTERACTION_NEWTON,
           "eih": lib.MYB_INTERACTION_EIH}
_schemes = {"dkd": lib.MYB_WH_DKD, "saba2": lib.MYB_WH_SABA2, "saba3": lib.MYB_WH_SABA3,
            "saba4": lib.MYB_WH_SABA4, "sbab1": lib.MYB_WH_SBAB1, "sbab2": lib.MYB_WH_SBAB2,
            "sbab3": lib.MYB_WH_SBAB3, "sbab4": lib.MYB_WH_SBAB4}
_shapes = {"plane": lib.MYB_CONTACT_PLANE, "ellipsoid": lib.MYB_CONTACT_ELLIPSOID,
           "mesh": lib.MYB_CONTACT_MESH}


def _nothing(pointer):
    pass


def _positions(values, count, what):
    return _doubles(values, 3 * count, what)


class Pool(_Owner):
    """Threads for the accelerations of an NBody (myb_pool): threads working on a job, the
    calling thread included (0: one per logical processor)."""

    _destroy = staticmethod(lib.myb_pool_destroy)

    def __init__(self, threads=0):
        out = ffi.new("myb_pool**")
        _check(lib.myb_pool_create(int(threads), ffi.NULL, out))
        super().__init__(out[0])

    @property
    def threads(self):
        return lib.myb_pool_threads(self.pointer)


class NBody(_Owner):
    """The Newtonian N-body problem (myb_nbody): len(gm) massive bodies, then `test` test
    bodies that feel them and attract nothing."""

    _destroy = staticmethod(lib.myb_nbody_destroy)

    def __init__(self, gm, test=0):
        gm = _flat(gm)
        out = ffi.new("myb_nbody**")
        _check(lib.myb_nbody_create(len(gm), int(test), ffi.new("double[]", gm or 1), ffi.NULL,
                                    out))
        super().__init__(out[0])
        self._owner = None

    @classmethod
    def _borrow(cls, pointer, owner):
        """The point masses inside a System: they live as long as the system does."""
        made = cls.__new__(cls)
        made._pointer = ffi.gc(pointer, _nothing)
        made._kept = [owner]
        made._owner = owner
        return made

    @property
    def pointer(self):
        if self._owner is not None and self._owner._pointer is None:
            raise ValueError("the System of this NBody is closed")
        return _Owner.pointer.fget(self)

    @property
    def count(self):
        return lib.myb_nbody_count(self.pointer)

    @property
    def massive(self):
        return lib.myb_nbody_massive(self.pointer)

    def gm(self, body):
        return lib.myb_nbody_gm(self.pointer, int(body))

    def set_pool(self, pool):
        """Threads for the accelerations (None: the calling thread only)."""
        lib.myb_nbody_set_pool(self.pointer, ffi.NULL if pool is None else pool.pointer)
        self._kept.append(pool)

    def set_interaction(self, receiver, source, level):
        """How source acts on receiver: "none", "newton", "eih" (or lib.MYB_INTERACTION_*)."""
        level = _levels[level] if isinstance(level, str) else int(level)
        _check(lib.myb_nbody_set_interaction(self.pointer, int(receiver), int(source), level))

    def set_relativity(self, c=None, beta=1.0, gamma=1.0):
        """The 1PN terms with the speed of light c and the PPN beta, gamma (c None: back to
        Newtonian)."""
        ppn = ffi.NULL if c is None else ffi.new("myb_ppn*", {"c": float(c), "beta": float(beta),
                                                              "gamma": float(gamma)})
        _check(lib.myb_nbody_set_relativity(self.pointer, ppn))

    def acceleration(self, t, x):
        """The accelerations at the positions x (myb_nbody_acceleration): a list of triples."""
        n = self.count
        a = ffi.new("double[]", 3 * n)
        _check(lib.myb_nbody_acceleration(self.pointer, float(t), _positions(x, n, "x"), a))
        return _rows(a, n, 3)

    def acceleration_xv(self, t, x, v):
        """The same with the velocities, which the 1PN terms need."""
        n = self.count
        a = ffi.new("double[]", 3 * n)
        _check(lib.myb_nbody_acceleration_xv(self.pointer, float(t), _positions(x, n, "x"),
                                             _positions(v, n, "v"), a))
        return _rows(a, n, 3)

    def invariants(self, x, v):
        """The conserved quantities of the massive bodies (energy a pair (hi, lo))."""
        n = self.count
        out = ffi.new("myb_nbody_invariants*")
        _check(lib.myb_nbody_invariants_of(self.pointer, _positions(x, n, "x"),
                                           _positions(v, n, "v"), out))
        return NBodyInvariants((out.energy.hi, out.energy.lo), tuple(out.momentum),
                               tuple(out.angular_momentum), tuple(out.center),
                               tuple(out.center_velocity))

    def pn_invariants(self, x, v):
        """The conserved quantities of the 1PN problem."""
        n = self.count
        out = ffi.new("myb_nbody_pn_invariants*")
        _check(lib.myb_nbody_pn_invariants_of(self.pointer, _positions(x, n, "x"),
                                              _positions(v, n, "v"), out))
        return NBodyPNInvariants((out.energy.hi, out.energy.lo), tuple(out.momentum),
                                 tuple(out.center), tuple(out.center_velocity))

    def set_approximation(self, order, theta=0.0):
        """The fast multipole method of the order (1 .. 20) and opening angle theta; order 0:
        back to exact direct sums."""
        _check(lib.myb_nbody_set_approximation(self.pointer, int(order), float(theta)))

    @property
    def approximation(self):
        """(order, theta) of the approximation, (0, 0.0) for exact sums."""
        theta = ffi.new("double*")
        order = lib.myb_nbody_approximation(self.pointer, theta)
        return order, theta[0]


class WisdomHolman(_Owner):
    """A Wisdom-Holman splitting integrator of an NBody (myb_wh): step h, scheme "dkd",
    "saba2" .. "saba4", "sbab1" .. "sbab4", a symplectic corrector of odd order 3 .. 17 for
    "dkd". The integrator keeps the system alive."""

    _destroy = staticmethod(lib.myb_wh_destroy)

    def __init__(self, nbody, h, scheme="dkd", corrector=0):
        options = ffi.new("myb_wh_options*", {"h": float(h), "corrector": int(corrector)})
        options.scheme = _schemes[scheme] if isinstance(scheme, str) else int(scheme)
        out = ffi.new("myb_wh**")
        _check(lib.myb_wh_create(nbody.pointer, options, ffi.NULL, out))
        super().__init__(out[0])
        self._kept.append(nbody)
        self._count = nbody.count

    def start(self, t0, x, v):
        n = self._count
        _check(lib.myb_wh_start(self.pointer, float(t0), _positions(x, n, "x"),
                                _positions(v, n, "v")))

    def step(self):
        _check(lib.myb_wh_step(self.pointer))

    def integrate(self, t_end):
        _check(lib.myb_wh_integrate(self.pointer, float(t_end)))

    @property
    def time(self):
        return lib.myb_wh_time(self.pointer)

    def state(self):
        """(positions, velocities) at the current node, corrected."""
        n = self._count
        x = ffi.new("double[]", 3 * n)
        v = ffi.new("double[]", 3 * n)
        _check(lib.myb_wh_state(self.pointer, x, v))
        return _rows(x, n, 3), _rows(v, n, 3)

    def dense(self, t):
        """(positions, velocities) at t within a step of the node."""
        n = self._count
        x = ffi.new("double[]", 3 * n)
        v = ffi.new("double[]", 3 * n)
        _check(lib.myb_wh_dense(self.pointer, float(t), x, v))
        return _rows(x, n, 3), _rows(v, n, 3)

    @property
    def statistics(self):
        s = lib.myb_wh_statistics(self.pointer)
        return WisdomHolmanStats(s.steps, s.evaluations, s.drifts)


# ---- The system of bodies ------------------------------------------------------------------------

def _ias15_options(epsilon=0.0, h_initial=0.0, h_max=0.0, h_min=0.0, max_steps=0):
    return ffi.new("myb_ias15_options*", {"epsilon": float(epsilon),
                                          "h_initial": float(h_initial), "h_max": float(h_max),
                                          "h_min": float(h_min), "max_steps": int(max_steps)})


class System(_Owner):
    """A system of bodies with their models (myb_system): len(gm) massive bodies. The point
    masses — their graph, the 1PN terms, the threads, the approximation — are `nbody`; every
    other level is a method here. Python callables (actions, orientations, guided bodies) may
    raise: the exception comes out of the call of the library that was evaluating them."""

    _destroy = staticmethod(lib.myb_system_destroy)

    def __init__(self, gm):
        gm = _flat(gm)
        out = ffi.new("myb_system**")
        _check(lib.myb_system_create(len(gm), ffi.new("double[]", gm or 1), ffi.NULL, out))
        super().__init__(out[0])
        self._calls = _Calls()
        self._nbody = None

    def _settle(self, status):
        _settle(status, self)

    @property
    def bodies(self):
        return lib.myb_system_bodies(self.pointer)

    @property
    def nbody(self):
        """The point masses (an NBody that lives as long as the system)."""
        if self._nbody is None:
            self._nbody = NBody._borrow(lib.myb_system_nbody(self.pointer), self)
        return self._nbody

    # ---- the model of a body

    def set_orientation(self, body, orientation):
        """The body-fixed frame from outside: a Provider of myb/frames.h or a callable
        f(seconds) -> Rotation or Orientation."""
        provider = _orientation_provider(orientation, self._calls)
        _check(lib.myb_system_set_orientation(self.pointer, int(body), provider.function,
                                              provider.context))
        self._kept.append(provider)

    def set_field(self, body, field, degree=0):
        """The figure of a body: a GravityField in its frame, up to `degree` (0: all)."""
        _check(lib.myb_system_set_field(self.pointer, int(body), field.pointer, int(degree)))
        self._kept.append(field)

    def set_rotation(self, body, moments, radius=0.0, k2=0.0, lag=0.0, mean_motion=0.0,
                     core=(0.0, 0.0), core_friction=0.0):
        """An integrated rotation: moments A <= B <= C per unit mass (m^2), the distortion
        (radius, k2, lag, mean_motion) and a fluid core (A_c, C_c) with its friction."""
        rotation = ffi.new("myb_system_rotation*")
        rotation.moments = list(_doubles(moments, 3, "moments"))
        rotation.radius = float(radius)
        rotation.k2 = float(k2)
        rotation.lag = float(lag)
        rotation.mean_motion = float(mean_motion)
        rotation.core = list(_doubles(core, 2, "core"))
        rotation.core_friction = float(core_friction)
        _check(lib.myb_system_set_rotation(self.pointer, int(body), rotation))

    def set_tide(self, body, radius, k2, orbital_lag=(0.0, 0.0, 0.0),
                 rotation_lag=(0.0, 0.0, 0.0), rotation_rate=0.0):
        """The answer of a body to tides by band: Love numbers k2 and lags of the three
        bands."""
        tide = ffi.new("myb_system_tide*")
        tide.radius = float(radius)
        tide.k2 = list(_doubles(k2, 3, "k2"))
        tide.orbital_lag = list(_doubles(orbital_lag, 3, "orbital_lag"))
        tide.rotation_lag = list(_doubles(rotation_lag, 3, "rotation_lag"))
        tide.rotation_rate = float(rotation_rate)
        _check(lib.myb_system_set_tide(self.pointer, int(body), tide))

    def set_mass_rate(self, body, rate, epoch=0.0):
        _check(lib.myb_system_set_mass_rate(self.pointer, int(body), float(rate), float(epoch)))

    def set_radius(self, body, radius):
        _check(lib.myb_system_set_radius(self.pointer, int(body), float(radius)))

    def radius(self, body):
        return lib.myb_system_radius(self.pointer, int(body))

    def present(self, body):
        """False for a body merged into another."""
        return lib.myb_system_present(self.pointer, int(body)) != 0

    # ---- the levels of the graph

    def add_field(self, body, point):
        _check(lib.myb_system_add_field(self.pointer, int(body), int(point)))

    def add_figures(self, a, b):
        _check(lib.myb_system_add_figures(self.pointer, int(a), int(b)))

    def add_tide(self, body, raiser, point):
        _check(lib.myb_system_add_tide(self.pointer, int(body), int(raiser), int(point)))

    def add_distortion(self, body, raiser):
        _check(lib.myb_system_add_distortion(self.pointer, int(body), int(raiser)))

    def add_pn_pair(self, a, b, orders, c=0.0, r0=0.0, kappa1=0.0, kappa2=0.0):
        """The terms beyond 1PN of two bodies: orders lib.MYB_PN_2 | ... (or a sequence of
        them)."""
        if not isinstance(orders, int):
            bits = 0
            for order in orders:
                bits |= int(order)
            orders = bits
        pair = ffi.new("myb_system_pn_pair*", {"orders": int(orders), "c": float(c),
                                               "r0": float(r0), "kappa1": float(kappa1),
                                               "kappa2": float(kappa2)})
        _check(lib.myb_system_add_pn_pair(self.pointer, int(a), int(b), pair))

    def add_frame_dragging(self, body, point, spin):
        _check(lib.myb_system_add_frame_dragging(self.pointer, int(body), int(point),
                                                 _doubles(spin, 3, "spin")))

    def add_action(self, body, action):
        """What pushes or turns a body otherwise: a callable action(t, view) -> (acceleration,
        torque), view a SystemView (body, bodies, x, v, gm, frame: nine numbers or None), the
        acceleration in inertial axes (m/s^2) and the torque per unit mass in the body's frame
        (m^2/s^2); or a Provider of a myb_system_action of the C interface."""
        if isinstance(action, Provider):
            provider = action
        else:
            calls = self._calls

            @ffi.callback("myb_status(void*, double, const myb_system_view*, double[3], "
                          "double[3])")
            def function(context, t, view, acceleration, torque):
                try:
                    n = view.bodies
                    seen = SystemView(view.body, n, _rows(view.x, n, 3), _rows(view.v, n, 3),
                                      _tuple(view.gm, n),
                                      None if view.frame == ffi.NULL else _tuple(view.frame, 9))
                    a, m = action(t, seen)
                    acceleration[0:3] = list(_doubles(a, 3, "acceleration"))
                    torque[0:3] = list(_doubles(m, 3, "torque"))
                    return lib.MYB_OK
                except BaseException as error:  # noqa: B902 - nothing may cross the C code
                    return calls.failed(error)

            provider = Provider(function, ffi.NULL, action)
        _check(lib.myb_system_add_action(self.pointer, int(body), provider.function,
                                         provider.context))
        self._kept.append(provider)

    def add_light_action(self, body, source, shape, flux, distance, c, mass, specular=0.0,
                         diffuse=0.0, thermal=0, shadows=0, theta=0.0, spin_axis=(0.0, 0.0, 1.0),
                         phases=0):
        """Light of the body `source` on the mesh `shape` (a Polyhedron) of the body: the
        ready-made action myb_system_light_action."""
        light = ffi.new("myb_system_light*")
        light.source = int(source)
        light.shape = shape.pointer
        light.surface.specular = float(specular)
        light.surface.diffuse = float(diffuse)
        light.surface.thermal = int(thermal)
        light.flux = float(flux)
        light.distance = float(distance)
        light.c = float(c)
        light.mass = float(mass)
        light.shadows = int(shadows)
        light.theta = float(theta)
        light.spin_axis = list(_doubles(spin_axis, 3, "spin_axis"))
        light.phases = int(phases)
        self.add_action(body, Provider(lib.myb_system_light_action, light, shape))

    def add_guided(self, gm, state):
        """A body whose motion is given — a body of an Ephemeris, a Provider or a callable
        f(seconds) -> six numbers — that attracts the system; its index among the guided
        bodies."""
        provider = _state_provider(state, self._calls)
        index = ffi.new("size_t*")
        _check(lib.myb_system_add_guided(self.pointer, float(gm), provider.function,
                                         provider.context, index))
        self._kept.append(provider)
        return index[0]

    @property
    def guided(self):
        return lib.myb_system_guided(self.pointer)

    def add_contact(self, body, ground, surface, law, points):
        """`body` touches the surface of `ground`. surface: a mapping with shape ("plane",
        "ellipsoid", "mesh"), point and normal (plane), radii (ellipsoid), mesh (a
        Polyhedron); law: stiffness, exponent, damping, friction, slip; points: pairs (offset,
        radius) in the frame of the body."""
        points = list(points)
        made = ffi.new("myb_contact_point[]", max(len(points), 1))
        for k, (offset, radius) in enumerate(points):
            made[k].offset = list(_doubles(offset, 3, "offset"))
            made[k].radius = float(radius)
        contact = ffi.new("myb_system_contact*")
        shape = surface.get("shape", "plane")
        contact.surface.shape = _shapes[shape] if isinstance(shape, str) else int(shape)
        contact.surface.point = list(_doubles(surface.get("point", (0.0, 0.0, 0.0)), 3, "point"))
        contact.surface.normal = list(_doubles(surface.get("normal", (0.0, 0.0, 0.0)), 3,
                                               "normal"))
        contact.surface.radii = list(_doubles(surface.get("radii", (0.0, 0.0, 0.0)), 3, "radii"))
        mesh = surface.get("mesh")
        if mesh is not None:
            contact.surface.mesh = mesh.pointer
        for key in ("stiffness", "exponent", "damping", "friction", "slip"):
            setattr(contact.law, key, float(law.get(key, 0.0)))
        contact.points = made
        contact.point_count = len(points)
        _check(lib.myb_system_add_contact(self.pointer, int(body), int(ground), contact))
        self._kept += [made, mesh]

    # ---- the state

    @property
    def coordinates(self):
        """n: 3 per body, 9 per integrated rotation, 6 per PN pair with spins."""
        return lib.myb_system_coordinates(self.pointer)

    def _state(self, y):
        return _doubles(y, 2 * self.coordinates, "a state of the system")

    def pack_rotation(self, y, body, q, omega, omega_core=None):
        """The state y with the rotation of `body` put in: a new tuple."""
        made = self._state(y)
        _check(lib.myb_system_pack_rotation(self.pointer, int(body), _doubles(q, 4, "q"),
                                            _doubles(omega, 3, "omega"),
                                            _optional(omega_core, 3, "omega_core"), made))
        return _tuple(made, 2 * self.coordinates)

    def unpack_rotation(self, y, body):
        """(q, omega, omega_core) of `body` in the state y."""
        q = ffi.new("double[4]")
        omega = ffi.new("double[3]")
        core = ffi.new("double[3]")
        _check(lib.myb_system_unpack_rotation(self.pointer, int(body), self._state(y), q, omega,
                                              core))
        return tuple(q), tuple(omega), tuple(core)

    def pack_spins(self, y, pair, chi1, chi2):
        made = self._state(y)
        _check(lib.myb_system_pack_spins(self.pointer, int(pair), _doubles(chi1, 3, "chi1"),
                                         _doubles(chi2, 3, "chi2"), made))
        return _tuple(made, 2 * self.coordinates)

    def unpack_spins(self, y, pair):
        chi1 = ffi.new("double[3]")
        chi2 = ffi.new("double[3]")
        _check(lib.myb_system_unpack_spins(self.pointer, int(pair), self._state(y), chi1, chi2))
        return tuple(chi1), tuple(chi2)

    def integrator(self, epsilon=0.0, h_initial=0.0, h_max=0.0, h_min=0.0, max_steps=0):
        """An IAS15 integrator of the system as it stands (myb_system_integrator)."""
        return SystemIntegrator(self, _ias15_options(epsilon, h_initial, h_max, h_min,
                                                     max_steps))

    def acceleration(self, t, x, v):
        """The second derivatives of the n coordinates (myb_system_acceleration): a tuple."""
        n = self.coordinates
        a = ffi.new("double[]", max(n, 1))
        self._settle(lib.myb_system_acceleration(self.pointer, float(t), _doubles(x, n, "x"),
                                                 _doubles(v, n, "v"), a))
        return _tuple(a, n)

    def inertia(self, body, t, y):
        """The inertia tensor (nine numbers, m^2) of an integrated rotation at (t, y)."""
        out = ffi.new("double[9]")
        self._settle(lib.myb_system_inertia(self.pointer, int(body), float(t), self._state(y),
                                            out))
        return tuple(out)

    def partials(self, t0, y0, t1, bodies=(), epsilon=0.0, h_initial=0.0, h_max=0.0, h_min=0.0,
                 max_steps=0):
        """(y1, phi, sensitivity): the state at t1, d y1 / d y0 (2 n rows of 2 n) and
        d y1 / d gm of the bodies given (2 n rows of len(bodies)), by central differences of
        whole propagations (myb_system_partials)."""
        n2 = 2 * self.coordinates
        bodies = [int(b) for b in bodies]
        chosen = ffi.new("size_t[]", bodies or 1)
        y1 = ffi.new("double[]", n2)
        phi = ffi.new("double[]", n2 * n2)
        sensitivity = ffi.new("double[]", max(n2 * len(bodies), 1))
        self._settle(lib.myb_system_partials(self.pointer, _ias15_options(epsilon, h_initial,
                                                                          h_max, h_min,
                                                                          max_steps),
                                             float(t0), self._state(y0), float(t1), chosen,
                                             len(bodies), y1, phi,
                                             sensitivity if bodies else ffi.NULL))
        return _tuple(y1, n2), _rows(phi, n2, n2), _rows(sensitivity, n2, len(bodies))

    # ---- encounters

    def advance(self, integrator, t_end, watch=0.0):
        """Steps the integrator towards t_end up to the first contact of two bodies with a
        radius or, with watch > 0, the closest approach of two within watch: a SystemEvent
        (kind lib.MYB_SYSTEM_ARRIVED, _CONTACT or _APPROACH)."""
        event = ffi.new("myb_system_event*")
        self._settle(lib.myb_system_advance(self.pointer, integrator.pointer, float(t_end),
                                            float(watch), event))
        return SystemEvent(event.kind, event.a, event.b, event.t, event.distance)

    def merge(self, integrator, a, b):
        """Body b merged into body a; what the merger took out of the motion."""
        merger = ffi.new("myb_system_merger*")
        self._settle(lib.myb_system_merge(self.pointer, integrator.pointer, int(a), int(b),
                                          merger))
        return SystemMerger(merger.energy, tuple(merger.angular_momentum), tuple(merger.spin))

    def bounce(self, integrator, a, b, restitution):
        self._settle(lib.myb_system_bounce(self.pointer, integrator.pointer, int(a), int(b),
                                           float(restitution)))


class SystemIntegrator(_Owner):
    """The IAS15 integrator of a System (myb_system_integrator, an myb_ias15 of myb/integ.h)
    that integrates myb_system_acceleration of that system; it keeps the system alive."""

    _destroy = staticmethod(lib.myb_ias15_destroy)

    def __init__(self, system, options):
        out = ffi.new("myb_ias15**")
        _check(lib.myb_system_integrator(system.pointer, options, ffi.NULL, out))
        super().__init__(out[0])
        self._system = system
        self._kept.append(system)
        self._n = system.coordinates

    def start(self, t0, y0):
        self._system._settle(lib.myb_ias15_start(self.pointer, lib.myb_system_acceleration,
                                                 self._system.pointer, float(t0),
                                                 _doubles(y0, 2 * self._n, "y0")))

    def step(self, t_limit):
        self._system._settle(lib.myb_ias15_step(self.pointer, float(t_limit)))

    def integrate(self, t_end):
        self._system._settle(lib.myb_ias15_integrate(self.pointer, float(t_end)))

    @property
    def time(self):
        return lib.myb_ias15_time(self.pointer)

    def state(self):
        """The 2 n numbers of the state now."""
        return _tuple(lib.myb_ias15_state(self.pointer), 2 * self._n)

    def dense(self, t):
        y = ffi.new("double[]", 2 * self._n)
        _check(lib.myb_ias15_dense(self.pointer, float(t), y))
        return _tuple(y, 2 * self._n)

    @property
    def statistics(self):
        s = lib.myb_ias15_statistics(self.pointer)
        return Ias15Stats(s.accepted, s.rejected, s.evaluations, s.iterations, s.unconverged)
