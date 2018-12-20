"""
vec3.py: Defines the Vec3 class used by OpenMM

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import absolute_import, division
__author__ = "Peter Eastman"
__version__ = "1.0"

from .. import unit
from collections import namedtuple

class Vec3(namedtuple('Vec3', ['x', 'y', 'z'])):
    """Vec3 is a 3-element tuple that supports many math operations."""

    def __new__(cls, x, y, z):
        """Create a new Vec3."""
        return tuple.__new__(cls, (x, y, z))

    def __getnewargs__(self):
        "Support for pickle protocol 2: http://docs.python.org/2/library/pickle.html#pickling-and-unpickling-normal-class-instances"
        return self[0], self[1], self[2]

    def __add__(self, other):
        """Add two Vec3s."""
        return Vec3(self.x+other[0], self.y+other[1], self.z+other[2])

    def __radd__(self, other):
        """Add two Vec3s."""
        return Vec3(self.x+other[0], self.y+other[1], self.z+other[2])

    def __sub__(self, other):
        """Add two Vec3s."""
        return Vec3(self.x-other[0], self.y-other[1], self.z-other[2])

    def __rsub__(self, other):
        """Add two Vec3s."""
        return Vec3(other[0]-self.x, other[1]-self.y, other[2]-self.z)

    def __mul__(self, other):
        """Multiply a Vec3 by a constant."""
        if unit.is_unit(other):
            return unit.Quantity(self, other)
        return Vec3(other*self.x, other*self.y, other*self.z)

    def __rmul__(self, other):
        """Multiply a Vec3 by a constant."""
        if unit.is_unit(other):
            return unit.Quantity(self, other)
        return Vec3(other*self.x, other*self.y, other*self.z)

    def __div__(self, other):
        """Divide a Vec3 by a constant."""
        return Vec3(self.x/other, self.y/other, self.z/other)
    __truediv__ = __div__

    def __deepcopy__(self, memo):
        return Vec3(self.x, self.y, self.z)

    def __neg__(self):
        return Vec3(-self.x, -self.y, -self.z)
