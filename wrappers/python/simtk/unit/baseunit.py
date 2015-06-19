#!/bin/env python


"""
Module simtk.unit.baseunit

Contains BaseUnit class, which is a component of the Unit class.

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
from __future__ import print_function, division, absolute_import

__author__ = "Christopher M. Bruns"
__version__ = "0.6"

class BaseUnit(object):
    '''
    Physical unit expressed in exactly one BaseDimension.

    For example, meter_base_unit could be a BaseUnit for the length dimension.
    The BaseUnit class is used internally in the more general Unit class.
    '''
    __array_priority__ = 100

    def __init__(self, base_dim, name, symbol):
        """Creates a new BaseUnit.

        Parameters
         - self: The newly created BaseUnit.
         - base_dim: (BaseDimension) The dimension of the new unit, e.g. 'mass'
         - name: (string) Name of the unit, e.g. "kilogram"
         - symbol: (string) Symbol for the unit, e.g. 'kg'.  This symbol will appear in
            Quantity string descriptions.
        """
        self.dimension = base_dim
        self.name = name
        self.symbol = symbol
        self._conversion_factor_to = {}
        self._conversion_factor_to[self] = 1.0
        self._conversion_factor_to_by_name = {}
        self._conversion_factor_to_by_name[self.name] = 1.0

    def __lt__(self, other):
        """
        Comparison function that sorts BaseUnits by BaseDimension
        """
        # First sort on dimension
        if self.dimension != other.dimension:
            return self.dimension < other.dimension
        # Second on conversion factor
        return self.conversion_factor_to(other) < 1.0

    def iter_base_dimensions(self):
        """
        Returns a dictionary of BaseDimension:exponent pairs, describing the dimension of this unit.
        """
        yield (self.dimension, 1)

    def iter_base_units(self):
        yield (self, 1)

    def get_dimension_tuple(self):
        """
        Returns a sorted tuple of (BaseDimension, exponent) pairs, that can be used as a dictionary key.
        """
        l = list(self.iter_base_dimensions())
        l.sort()
        return tuple(l)

    def __str__(self):
        """Returns a string with the name of this BaseUnit
        """
        return self.name

    def __repr__(self):
        return 'BaseUnit(base_dim=%s, name="%s", symbol="%s")' % (self.dimension, self.name, self.symbol)

    def define_conversion_factor_to(self, other, factor):
        """
        Defines a conversion factor between two BaseUnits.

        self * factor = other

        Parameters:
         - self: (BaseUnit) 'From' unit in conversion.
         - other: (BaseUnit) 'To' unit in conversion.
         - factor: (float) Conversion factor.

        After calling this method, both self and other will have stored
        conversion factors for one another, plus all other BaseUnits which
        self and other have previously defined.

        Both self and other must have the same dimension, otherwise a TypeError
        will be raised.

        Returns None.
        """
        if self.dimension != other.dimension:
            raise TypeError('Cannot define conversion for BaseUnits with different dimensions.')
        assert(factor != 0)
        assert(not self is other)
        # import all transitive conversions
        self._conversion_factor_to[other] = factor
        self._conversion_factor_to_by_name[other.name] = factor
        for (unit, cfac) in other._conversion_factor_to.items():
            if unit is self: continue
            if unit in self._conversion_factor_to: continue
            self._conversion_factor_to[unit] = factor * cfac
            unit._conversion_factor_to[self] = pow(factor * cfac, -1)
            self._conversion_factor_to_by_name[unit.name] = factor * cfac
            unit._conversion_factor_to_by_name[self.name] = pow(factor * cfac, -1)
        # and for the other guy
        invFac = pow(factor, -1.0)
        other._conversion_factor_to[self] = invFac
        other._conversion_factor_to_by_name[self.name] = invFac
        for (unit, cfac) in self._conversion_factor_to.items():
            if unit is other: continue
            if unit in other._conversion_factor_to: continue
            other._conversion_factor_to[unit] = invFac * cfac
            unit._conversion_factor_to[other] = pow(invFac * cfac, -1)
            other._conversion_factor_to_by_name[unit.name] = invFac * cfac
            unit._conversion_factor_to_by_name[other.name] = pow(invFac * cfac, -1)

    def conversion_factor_to(self, other):
        """Returns a conversion factor from this BaseUnit to another BaseUnit.

        It does not matter which existing BaseUnit you define the conversion factor to.
        Conversions for all other known BaseUnits will be computed at the same time.

        Raises TypeError if dimension does not match.
        Raises LookupError if no conversion has been defined. (see define_conversion_factor_to).

        """
        if self is other: return 1.0
        if self.dimension != other.dimension:
            raise TypeError('Cannot get conversion for BaseUnits with different dimensions.')
        if not other.name in self._conversion_factor_to_by_name:
            raise LookupError('No conversion defined from BaseUnit "%s" to "%s".' % (self, other))
        return self._conversion_factor_to_by_name[other.name]

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
