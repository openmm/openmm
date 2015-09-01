"""
unitcell.py: Routines for converting between different representations of the periodic unit cell.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2015 Stanford University and the Authors.
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
from __future__ import absolute_import
__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm import Vec3
from simtk.unit import nanometers, is_quantity, norm, dot, radians
import math


def computePeriodicBoxVectors(a_length, b_length, c_length, alpha, beta, gamma):
    """Convert lengths and angles to periodic box vectors.
    
    Lengths should be given in nanometers and angles in radians (or as Quantity
    instances)
    """

    if is_quantity(a_length): a_length = a_length.value_in_unit(nanometers)
    if is_quantity(b_length): b_length = b_length.value_in_unit(nanometers)
    if is_quantity(c_length): c_length = c_length.value_in_unit(nanometers)
    if is_quantity(alpha): alpha = alpha.value_in_unit(radians)
    if is_quantity(beta): beta = beta.value_in_unit(radians)
    if is_quantity(gamma): gamma = gamma.value_in_unit(radians)

    # Compute the vectors.

    a = [a_length, 0, 0]
    b = [b_length*math.cos(gamma), b_length*math.sin(gamma), 0]
    cx = c_length*math.cos(beta)
    cy = c_length*(math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma)
    cz = math.sqrt(c_length*c_length-cx*cx-cy*cy)
    c = [cx, cy, cz]

    # If any elements are very close to 0, set them to exactly 0.

    for i in range(3):
        if abs(a[i]) < 1e-6:
            a[i] = 0.0
        if abs(b[i]) < 1e-6:
            b[i] = 0.0
        if abs(c[i]) < 1e-6:
            c[i] = 0.0
    a = Vec3(*a)
    b = Vec3(*b)
    c = Vec3(*c)

    # Make sure they're in the reduced form required by OpenMM.

    c = c - b*round(c[1]/b[1])
    c = c - a*round(c[0]/a[0])
    b = b - a*round(b[0]/a[0])
    return (a, b, c)*nanometers

def reducePeriodicBoxVectors(periodicBoxVectors):
    """ Reduces the representation of the PBC. periodicBoxVectors is expected to
    be an unpackable iterable of length-3 iterables
    """
    if is_quantity(periodicBoxVectors):
        a, b, c = periodicBoxVectors.value_in_unit(nanometers)
    else:
        a, b, c = periodicBoxVectors
    a = Vec3(*a)
    b = Vec3(*b)
    c = Vec3(*c)

    c = c - b*round(c[1]/b[1])
    c = c - a*round(c[0]/a[0])
    b = b - a*round(b[0]/a[0])

    return (a, b, c) * nanometers

def computeLengthsAndAngles(periodicBoxVectors):
    """Convert periodic box vectors to lengths and angles.

    Lengths are returned in nanometers and angles in radians.
    """
    if is_quantity(periodicBoxVectors):
        (a, b, c) = periodicBoxVectors.value_in_unit(nanometers)
    else:
        a, b, c = periodicBoxVectors
    a_length = norm(a)
    b_length = norm(b)
    c_length = norm(c)
    alpha = math.acos(dot(b, c)/(b_length*c_length))
    beta = math.acos(dot(c, a)/(c_length*a_length))
    gamma = math.acos(dot(a, b)/(a_length*b_length))
    return (a_length, b_length, c_length, alpha, beta, gamma)
