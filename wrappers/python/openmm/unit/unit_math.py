#!/bin/env python
"""
Module openmm.unit.math

Arithmetic methods on Quantities and Units

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Christopher M. Bruns
Contributors: Peter Eastman

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
from __future__ import division, print_function, absolute_import

__author__ = "Christopher M. Bruns"
__version__ = "0.5"


import math
from .quantity import is_quantity
from .unit_definitions import *

####################
### TRIGONOMETRY ###
####################

def sin(angle):
    """
    Examples

    >>> sin(90*degrees)
    1.0
    """
    if is_quantity(angle):
        return math.sin(angle/radians)
    else:
        return math.sin(angle)

def sinh(angle):
    if is_quantity(angle):
        return math.sinh(angle/radians)
    else:
        return math.sinh(angle)

def cos(angle):
    """
    Examples

    >>> cos(180*degrees)
    -1.0
    """
    if is_quantity(angle):
        return math.cos(angle/radians)
    else:
        return math.cos(angle)

def cosh(angle):
    if is_quantity(angle):
        return math.cosh(angle/radians)
    else:
        return math.cosh(angle)

def tan(angle):
    if is_quantity(angle):
        return math.tan(angle/radians)
    else:
        return math.tan(angle)

def tanh(angle):
    if is_quantity(angle):
        return math.tanh(angle/radians)
    else:
        return math.tanh(angle)

def acos(x):
    """
    >>> acos(1.0)
    Quantity(value=0.0, unit=radian)
    >>> print(acos(1.0))
    0.0 rad
    """
    return math.acos(x) * radians

def acosh(x):
    return math.acosh(x) * radians

def asin(x):
    return math.asin(x) * radians

def asinh(x):
    return math.asinh(x) * radians

def atan(x):
    return math.atan(x) * radians

def atanh(x):
    return math.atanh(x) * radians

def atan2(x, y):
    return math.atan2(x, y) * radians

###################
### SQUARE ROOT ###
###################

def sqrt(val):
    """
    >>> sqrt(9.0)
    3.0
    >>> print(sqrt(meter*meter))
    meter
    >>> sqrt(9.0*meter*meter)
    Quantity(value=3.0, unit=meter)
    >>> sqrt(9.0*meter*meter*meter)
    Traceback (most recent call last):
    ...
    ArithmeticError: Exponents in Unit.sqrt() must be even.
    """
    try:
        return val.sqrt()
    except AttributeError:
        return math.sqrt(val)

###########
### SUM ###
###########

def sum(val):
    """
    >>> sum((1.0, 2.0))
    3.0
    >>> sum((2.0*meter, 3.0*meter))
    Quantity(value=5.0, unit=meter)
    >>> sum((2.0*meter, 30.0*centimeter))
    Quantity(value=2.3, unit=meter)
    """
    try:
        return val.sum()
    except AttributeError:
        pass
    if len(val) == 0:
        return 0
    result = val[0]
    for i in range(1, len(val)):
        result += val[i]
    return result

###################
### VECTOR MATH ###
###################

def dot(x, y):
    """
    >>> dot((2, 3)*meter, (4, 5)*meter)
    Quantity(value=23, unit=meter**2)
    """
    sum = x[0]*y[0]
    for i in range(1, len(x)):
        sum += x[i]*y[i]
    return sum

def norm(x):
    """
    >>> norm((3, 4)*meter)
    Quantity(value=5.0, unit=meter)
    """
    return sqrt(dot(x, x))

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
