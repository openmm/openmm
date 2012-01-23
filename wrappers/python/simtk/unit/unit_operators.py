#!/bin/env python
"""
Module simtk.unit.unit_operators

Physical quantities with units, intended to produce similar functionality
to Boost.Units package in C++ (but with a runtime cost).
Uses similar API as Scientific.Physics.PhysicalQuantities
but different internals to satisfy our local requirements.
In particular, there is no underlying set of 'canonical' base
units, whereas in Scientific.Physics.PhysicalQuantities all
units are secretly in terms of SI units.  Also, it is easier
to add new fundamental dimensions to simtk.dimensions.  You
might want to make new dimensions for, say, "currency" or 
"information".

Two possible enhancements that have not been implemented are
  1) Include uncertainties with propagation of errors
  2) Incorporate offsets for celsius <-> kelvin conversion
"""

__author__ = "Christopher M. Bruns"
__version__ = "0.5"

from unit import Unit, is_unit
from quantity import Quantity, is_quantity

# Attach methods of Unit class that return a Quantity to Unit class.
# I put them here to avoid circular dependence in imports.
# i.e. Quantity depends on Unit, but not vice versa

def _unit_class_rdiv(self, other):
    """
    Divide another object type by a Unit.
    
    Returns a new Quantity with a value of other and units
    of the inverse of self.
    """
    if is_unit(other):
        raise NotImplementedError('programmer is surprised __rtruediv__ was called instead of __truediv__')
    else:
        # print "R scalar / unit"
        unit = pow(self, -1.0)
        value = other
        return Quantity(value, unit).reduce_unit(self)

Unit.__rtruediv__ = _unit_class_rdiv


def _unit_class_mul(self, other):
    """Multiply a Unit by an object.
    
    If other is another Unit, returns a new composite Unit.  
    Exponents of similar dimensions are added.  If self and 
    other share similar BaseDimension, but
    with different BaseUnits, the resulting BaseUnit for that
    BaseDimension will be that used in self.
    
    If other is a not another Unit, this method returns a 
    new Quantity...  UNLESS other is a Quantity and the resulting
    unit is dimensionless, in which case the underlying value type
    of the Quantity is returned.
    """
    if is_unit(other):
        if self in Unit._multiplication_cache:
            if other in Unit._multiplication_cache[self]:
                return Unit._multiplication_cache[self][other]
        else:
            Unit._multiplication_cache[self] = {}
        # print "unit * unit"
        result1 = {} # dictionary of dimensionTuple: (BaseOrScaledUnit, exponent)
        for unit, exponent in self.iter_base_or_scaled_units():
            d = unit.get_dimension_tuple()
            if d not in result1:
                result1[d] = {}
            assert unit not in result1[d]
            result1[d][unit] = exponent
        for unit, exponent in other.iter_base_or_scaled_units():
            d = unit.get_dimension_tuple()
            if d not in result1:
                result1[d] = {}
            if unit not in result1[d]:
                result1[d][unit] = 0
            result1[d][unit] += exponent
        result2 = {} # stripped of zero exponents
        for d in result1:
            for unit in result1[d]:
                exponent = result1[d][unit]
                if exponent != 0:
                    assert unit not in result2
                    result2[unit] = exponent
        new_unit = Unit(result2)
        Unit._multiplication_cache[self][other] = new_unit
        return new_unit
    elif is_quantity(other):
        # print "unit * quantity"
        value = other._value
        unit = self * other.unit
        return Quantity(value, unit).reduce_unit(self)
    else:
        # print "scalar * unit"
        value = other
        unit = self
        # Is reduce_unit needed here?  I hope not, there is a performance issue...
        # return Quantity(other, self).reduce_unit(self)
        return Quantity(other, self)
        
Unit.__mul__ = _unit_class_mul
Unit.__rmul__ = Unit.__mul__
Unit._multiplication_cache = {}


# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
