# SPDX-License-Identifier: BSD-3-Clause
"""The shared library and its declarations.

`ffi` knows every type and function of the public C interface (from _cdef.h, which is
generated from the headers); `lib` is the loaded library. The library is looked for, in
this order:

  1. the file named by the environment variable MYB_LIBRARY (and nowhere else if it is set);
  2. the directory _lib/ of this package (a library shipped with it);
  3. where the system looks for shared libraries.

The libraries of ext/ are looked for next to the library of the core, in the directories
listed in MYB_LIBRARY_PATH (separated as PATH is), in _lib/, and where the system looks.
"""
import os
import sys

from cffi import FFI

ffi = FFI()
_here = os.path.dirname(os.path.abspath(__file__))


def _declarations(name):
    with open(os.path.join(_here, name), encoding="ascii") as f:
        return f.read()


ffi.cdef(_declarations("_cdef.h"))


def _names(stem):
    if sys.platform == "win32":
        return ["lib%s.dll" % stem, "%s.dll" % stem]
    if sys.platform == "darwin":
        return ["lib%s.1.dylib" % stem, "lib%s.dylib" % stem]
    return ["lib%s.so.1" % stem, "lib%s.so" % stem]


def _open(candidates, what):
    errors = []
    for path in candidates:
        try:
            # a path is made absolute (the loader of Windows refuses a relative one);
            # a bare name is left to the system's search
            target = os.path.abspath(path) if os.path.dirname(path) else path
            return ffi.dlopen(target), target
        except OSError as error:
            errors.append("  %s: %s" % (path, error))
    raise ImportError(
        "the shared library %s was not found. Build myb with -DMYB_BUILD_SHARED=ON and set the "
        "environment variable MYB_LIBRARY to the path of the library, or put it into %s.\n%s"
        % (what, os.path.join(_here, "_lib"), "\n".join(errors)))


def _candidates():
    named = os.environ.get("MYB_LIBRARY")
    if named:
        return [named]
    names = _names("myb")
    return [os.path.join(_here, "_lib", name) for name in names] + names


lib, library_path = _open(_candidates(), "libmyb")

# The declarations belong to one release: the layouts of the structures may change between
# minor versions.
if lib.myb_version_number() // 10000 != lib.MYB_VERSION_NUMBER // 10000:
    raise ImportError(
        "%s is myb %d, the Python package was made for %d: they must agree in the major and "
        "minor version" % (library_path, lib.myb_version_number(), lib.MYB_VERSION_NUMBER))

_extensions = {}


def extension(name):
    """The library of an optional module of ext/ ("sgp4", "nrlmsise00", "troposphere",
    "ionosphere", "jb2008"): its functions are attributes of what is returned, its types are known to
    `ffi`. It is looked for next to the library of the core, then as that one is."""
    if name in _extensions:
        return _extensions[name]
    file = "_cdef_%s.h" % name
    if not os.path.exists(os.path.join(_here, file)):
        raise ValueError("no such module of ext/: %r" % (name,))
    ffi.cdef(_declarations(file))
    names = _names("myb_" + name)
    directories = [os.path.dirname(library_path)]
    directories += [d for d in os.environ.get("MYB_LIBRARY_PATH", "").split(os.pathsep) if d]
    directories.append(os.path.join(_here, "_lib"))
    found = [os.path.join(d, n) for d in directories if d for n in names]
    opened, _ = _open(found + names, "libmyb_" + name)
    _extensions[name] = opened
    return opened
