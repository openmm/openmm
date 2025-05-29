#!/bin/env python


"""
Module openmm.unit.baseunit

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

import collections

class BaseUnit(object):
    '''
    Physical unit expressed in exactly one BaseDimension.

    For example, meter_base_unit could be a BaseUnit for the length dimension.
    The BaseUnit class is used internally in the more general Unit class.
    '''
    __array_priority__ = 100

    # Global table of conversion factors between base units
    _conversion_factors = collections.defaultdict(lambda: collections.defaultdict(dict))

    def __init__(self, base_dim, name, symbol):
        """Creates a new BaseUnit.

        Parameters
         - self: The newly created BaseUnit.
         - base_dim: (BaseDimension) The dimension of the new unit, e.g. 'mass'
         - name: (string) Name of the unit, e.g. "kilogram".  This will be used
            to distinguish between BaseUnit objects with the same dimension: two
            BaseUnit objects with the same dimension that are given the same
            name will be treated as equal to each other.
         - symbol: (string) Symbol for the unit, e.g. 'kg'.  This symbol will appear in
            Quantity string descriptions.
        """
        self.dimension = base_dim
        self.name = name
        self.symbol = symbol
        BaseUnit._conversion_factors[self.dimension][self.name][self.name] = 1.0

    def __lt__(self, other):
        """
        Comparison function that sorts BaseUnits by BaseDimension
        """
        # First sort on dimension
        if self.dimension != other.dimension:
            return self.dimension < other.dimension
        # Second on conversion factor
        return self.conversion_factor_to(other) < 1.0

    def __eq__(self, other):
        if not isinstance(other, BaseUnit):
            return False
        return self.dimension == other.dimension and self.name == other.name

    def __hash__(self):
        return hash((self.dimension, self.name))

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
        assert(self != other)
        conversion_factors = BaseUnit._conversion_factors[self.dimension]
        conversion_factors_self = conversion_factors[self.name]
        conversion_factors_other = conversion_factors[other.name]
        # import all transitive conversions
        conversion_factors_self[other.name] = factor
        for (unit_name, cfac) in conversion_factors_other.items():
            if unit_name == self.name: continue
            if unit_name in conversion_factors_self: continue
            conversion_factors_self[unit_name] = factor * cfac
            conversion_factors[unit_name][self.name] = pow(factor * cfac, -1)
        # and for the other guy
        invFac = pow(factor, -1.0)
        conversion_factors_other[self.name] = invFac
        for (unit_name, cfac) in conversion_factors_self.items():
            if unit_name == other.name: continue
            if unit_name in conversion_factors_other: continue
            conversion_factors_other[unit_name] = invFac * cfac
            conversion_factors[unit_name][other.name] = pow(invFac * cfac, -1)

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
        if self.name == other.name: return 1.0
        conversion_factors_self = BaseUnit._conversion_factors[self.dimension][self.name]
        if not other.name in conversion_factors_self:
            raise LookupError('No conversion defined from BaseUnit "%s" to "%s".' % (self, other))
        return conversion_factors_self[other.name]

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
