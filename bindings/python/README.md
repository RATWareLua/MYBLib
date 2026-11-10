# myb for Python

The Python face of [myb](https://github.com/RATWareLua/MYBLib), a C library for the motion of
objects in space. The package loads the **shared library** of myb (it has no compiled part of
its own) and gives two layers:

- `myb.lib`, `myb.ffi` — every function and structure of the C interface, by the names of the
  C headers. The declarations are generated from the headers (`tools/make_cdef.py`).
- classes and functions for the path most programs take: `KernelPool`, `TimeSystem`,
  `Ephemeris`, `GravityField`, `Environment`, `ForceModel`, `Orbit`, `state_from_elements`,
  `kepler_propagate`, `write_oem`, `read_oem`. Errors are exceptions (`myb.Error`), terms,
  providers and events are Python callables.

No arithmetic is done in Python: a number has the bits of the same calls made in C.

## Use

Build the library with `-DMYB_BUILD_SHARED=ON` (the `*-release` presets do), then

    pip install ./bindings/python
    export MYB_LIBRARY=/path/to/libmyb.so.1        # libmyb.dll on Windows

or copy the library into `myb/_lib/` before building a wheel, or install it where the system
looks for shared libraries. The libraries of `ext/` are found next to it
(`myb.extension("sgp4")`), or in the directories of `MYB_LIBRARY_PATH`.

```python
import math, myb

gm, radius = 3.986004418e14, 6378137.0
field = myb.GravityField(gm, radius, 2, [1, 0, 0, -1.08263e-3 / math.sqrt(5), 0, 0], [0] * 6)
env = myb.Environment()
planet = env.add_body(1, gm)
model = myb.ForceModel()
model.add_field("gravity", env, planet, field)
model.add_drag("drag", env, planet, mass=500.0, area=4.0, cd=2.2,
               layers=[(400e3, 3e-12, 58e3)], equatorial_radius=radius)
orbit = myb.Orbit(model)
orbit.start(0.0, myb.state_from_elements(gm, (6.778e6, 0.001, 0.9, 0.5, 1.0, 0.0)), 86400.0)
node = orbit.next_event(lambda t, y: y[2], 86400.0, direction=1)   # (time, state) or None
print(orbit.state(86400.0))
```

Tests: `python -m unittest discover -s bindings/python/tests` with `MYB_LIBRARY` set.
