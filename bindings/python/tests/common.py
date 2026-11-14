# SPDX-License-Identifier: BSD-3-Clause
"""What the tests of the package share: the directories of the data, constants, and the
comparison of numbers bit for bit (the bytes of each double, so that -0.0 differs from 0.0
and a NaN equals itself)."""
import math
import os
import struct
import unittest

import myb
from myb import ffi, lib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
DATA = os.environ.get("MYB_DATA_DIR", os.path.join(ROOT, "data"))
TEST_DATA = os.environ.get("MYB_TEST_DATA_DIR", os.path.join(ROOT, "tests", "data"))

GM = 3.986004418e14
RADIUS = 6378137.0
ELEMENTS = myb.Elements(p=(RADIUS + 700000.0) * (1.0 - 0.01 ** 2), e=0.01,
                        i=98.0 * math.pi / 180.0, raan=0.3, argp=0.7, nu=0.2)


def bits(value):
    """Numbers as bytes, at any depth of tuples and lists; instants as (s, bytes of f)."""
    if isinstance(value, float):
        return struct.pack("<d", value)
    if isinstance(value, (list, tuple)):
        return [bits(v) for v in value]
    if isinstance(value, dict):
        return {k: bits(v) for k, v in value.items()}
    if isinstance(value, ffi.CData) and ffi.typeof(value).cname == "myb_time":
        return (value.s, struct.pack("<d", value.f))
    return value


def doubles(buffer, count):
    return tuple(ffi.unpack(buffer, count))


def rows(buffer, n, m):
    flat = ffi.unpack(buffer, n * m)
    return [tuple(flat[k * m:k * m + m]) for k in range(n)]


class Case(unittest.TestCase):
    def assertSame(self, a, b):
        """Equal bit for bit."""
        self.assertEqual(bits(a), bits(b))

    def assertStatus(self, status, call, *args, **kwargs):
        """The call raises myb.Error with this status of the library."""
        with self.assertRaises(myb.Error) as caught:
            call(*args, **kwargs)
        self.assertEqual(caught.exception.status, status)
        self.assertEqual(caught.exception.name, ffi.string(lib.myb_status_name(status)).decode())
        return caught.exception


def ok(status):
    assert status == lib.MYB_OK, ffi.string(lib.myb_status_name(status)).decode()


def kernels(*names):
    return [os.path.join(DATA, *name.split("/")) for name in names]
