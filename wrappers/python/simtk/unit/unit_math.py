#!/bin/env python
"""
Module simtk.unit.math

Arithmetic methods on Quantities and Units
"""

__author__ = "Christopher M. Bruns"
__version__ = "0.5"


import math
from quantity import is_quantity
from unit_definitions import *

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
    >>> print acos(1.0)
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
    >>> print sqrt(meter*meter)
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

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
