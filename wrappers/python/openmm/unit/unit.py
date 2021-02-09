#!/bin/env python
"""
Module openmm.unit

Contains classes Unit and ScaledUnit.

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
import sys
from .mymatrix import MyMatrix, zeros
from .baseunit import BaseUnit
from .standard_dimensions import *

class Unit(object):
    """
    Physical unit such as meter or ampere.
    """

    __array_priority__ = 100

    def __init__(self, base_or_scaled_units):
        """Create a new Unit.

        Parameters
        ----------
        self : Unit
            The newly created Unit.
        base_or_scaled_units : dict
            Keys are BaseUnits or ScaledUnits.  Values are exponents (numbers).
        """
        # Unit contents are of two types: BaseUnits and ScaledUnits
        self._top_base_units = {}
        self._all_base_units = {}
        self._scaled_units = []
        for (base_or_scaled_unit, power) in base_or_scaled_units.items():
            if power == 0:
                continue
            if isinstance(base_or_scaled_unit, BaseUnit):
                bu = base_or_scaled_unit
                dim = bu.dimension
                if dim not in self._top_base_units:
                    self._top_base_units[dim] = {}
                if bu not in self._top_base_units[dim]:
                    self._top_base_units[dim][bu] = 0
                self._top_base_units[dim][bu] += power
            else:
                self._scaled_units.append((base_or_scaled_unit, power))
        # Populate self._all_base_units
        # first, deep copy of self._top_base_units
        self._all_base_units = {}
        for d in self._top_base_units:
            self._all_base_units[d] = {}
            for u in self._top_base_units[d]:
                self._all_base_units[d][u] = self._top_base_units[d][u]
        # second, BaseUnits from self._scaled_units
        for scaled_unit, exponent1 in self._scaled_units:
            for base_unit, exponent2 in scaled_unit.iter_base_units():
                dim = base_unit.dimension
                if dim not in self._all_base_units:
                    self._all_base_units[dim] = {}
                if base_unit not in self._all_base_units[dim]:
                    self._all_base_units[dim][base_unit] = 0
                self._all_base_units[dim][base_unit] += exponent1 * exponent2
        # What about heterogenous units that cancel? --> leave them
        self._scaled_units.sort()

    def create_unit(self, scale, name, symbol):
        """
        Convenience method for creating a new simple unit from another simple unit.
        Both units must consist of a single BaseUnit.
        """
        # TODO - also handle non-simple units, i.e. units with multiple BaseUnits/ScaledUnits
        assert len(self._top_base_units) == 1
        assert len(self._scaled_units) == 0
        dimension = next(iter(self._top_base_units))
        base_unit_dict = self._top_base_units[dimension]
        assert len(base_unit_dict) == 1
        parent_base_unit = next(iter(base_unit_dict))
        parent_exponent = base_unit_dict[parent_base_unit]
        new_base_unit = BaseUnit(parent_base_unit.dimension, name, symbol)
        # BaseUnit scale might be different depending on exponent
        true_scale = scale
        if parent_exponent != 1.0:
            true_scale = math.pow(scale, 1.0/parent_exponent)
        new_base_unit.define_conversion_factor_to(parent_base_unit, true_scale)
        new_unit = Unit({new_base_unit: 1.0})
        return new_unit

    def iter_base_dimensions(self):
        """
        Yields (BaseDimension, exponent) tuples comprising this unit.
        """
        # There might be two units with the same dimension? No.
        for dimension in sorted(self._all_base_units.keys()):
            exponent = sum(self._all_base_units[dimension].values())
            if exponent != 0:
                yield (dimension, exponent)

    def iter_all_base_units(self):
        """
        Yields (BaseUnit, exponent) tuples comprising this unit, including those BaseUnits
        found within ScaledUnits.

        There might be multiple BaseUnits with the same dimension.
        """
        for dimension in sorted(self._all_base_units.keys()):
            for base_unit in sorted(self._all_base_units[dimension].keys()):
                exponent = self._all_base_units[dimension][base_unit]
                yield (base_unit, exponent)

    def iter_top_base_units(self):
        """
        Yields (BaseUnit, exponent) tuples in this Unit, excluding those within BaseUnits.
        """
        for dimension in sorted(self._top_base_units.keys()):
            for unit in sorted(self._top_base_units[dimension].keys()):
                exponent = self._top_base_units[dimension][unit]
                yield (unit, exponent)

    def iter_scaled_units(self):
        for unit, exponent in self._scaled_units:
            yield (unit, exponent)

    def iter_base_or_scaled_units(self):
        for item in self.iter_top_base_units():
            yield item
        for item in self.iter_scaled_units():
            yield item

    def get_conversion_factor_to_base_units(self):
        """
        There may be ScaleUnit components to this Unit.
        Returns conversion factor to the set of BaseUnits returned by iter_all_base_units().

        Units comprised of only BaseUnits return 1.0
        """
        factor = 1.0
        for scaled_unit, exponent in self._scaled_units:
            # print scaled_unit.factor
            factor *= scaled_unit.factor ** exponent
        return factor

    def __eq__(self, other):
        if not is_unit(other):
            return False
        return self.get_name() == other.get_name()

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        """Compare two Units.

        Raises a TypeError if the units have different dimensions.

        Returns True if self < other, False otherwise.
        """
        if not self.is_compatible(other):
            raise TypeError('Unit "%s" is not compatible with Unit "%s".', (self, other))
        return self.conversion_factor_to(other) < 1.0

    def __le__(self, other):
        return self.__lt__(other) or self.__eq__(other)

    def __gt__(self, other):
        return other.__lt__(self) 

    def __ge__(self, other):
        return other.__lt__(self) or other.__eq__(self)

    def __hash__(self):
        """
        Compute a hash code for this object.
        """
        try:
            return self._hash
        except AttributeError:
            pass
        self._hash = hash(self.get_name())
        return self._hash

    # def __mul__(self, other):
    # See unit_operators.py for Unit.__mul__ operator

    def __truediv__(self, other):
        """Divide a Unit by another object.

        Returns a composite Unit if other is another Unit.

        Returns a Quantity otherwise.  UNLESS other is a Quantity AND
        the resulting unit type is dimensionless, in which case the underlying
        value type of the Quantity is returned.
        """
        return self * pow(other, -1)

    __div__ = __truediv__

    # def __rtruediv__(self, other):
    # Because rtruediv returns a Quantity, look in quantity.py for definition of Unit.__rtruediv__

    _pow_cache = {}

    def __pow__(self, exponent):
        """Raise a Unit to a power.

        Returns a new Unit with different exponents on the BaseUnits.
        """
        if self in Unit._pow_cache:
            if exponent in Unit._pow_cache[self]:
                return Unit._pow_cache[self][exponent]
        else:
            Unit._pow_cache[self] = {}
        result = {} # dictionary of unit: exponent
        for unit, exponent2 in self.iter_base_or_scaled_units():
            result[unit] = exponent2 * exponent
        new_unit = Unit(result)
        Unit._pow_cache[self][exponent] = new_unit
        return new_unit

    def sqrt(self):
        """
        Returns square root of a unit.

        Raises ArithmeticError if component exponents are not even.
        This behavior can be changed if you present a reasonable real life case to me.
        """
        new_units = {}
        # There might be odd exponents in base and scaled units that
        # boil down to even exponents in base dimensions.
        # But if ScaledUnits and BaseUnits have even exponents, we should use them.
        nice_and_even = True
        for u, exponent in self.iter_base_or_scaled_units():
            if exponent%2 != 0:
                # This isn't going to work, we need to bust apart the ScaledUnits
                nice_and_even = False
                break
            new_units[u] = exponent/2
        if not nice_and_even:
            # Create a new unit formed from inner BaseUnits
            new_units = {}
            base_units_by_dimension = {}
            # Choose the first BaseUnit for each dimension
            for base_unit, exponent in self.iter_all_base_units():
                d = base_unit.dimension
                if d not in base_units_by_dimension:
                    base_units_by_dimension[d] = base_unit
                    new_units[base_unit] = exponent
                else:
                    # Already assigned a BaseUnit to this dimension, just update exponent
                    bu = base_units_by_dimension[d]
                    new_units[bu] += exponent
            # If exponents are not even by now, they never will be even
            for u, exponent in new_units.items():
                if exponent%2 != 0:
                    raise ArithmeticError('Exponents in Unit.sqrt() must be even.')
                new_units[u] = exponent/2
        return Unit(new_units)

    def __str__(self):
        """Returns the human-readable name of this unit"""
        return self.get_name()

    def __repr__(self):
        """
        Returns a unit name (string) for this Unit, composed of its various
        BaseUnit symbols.  e.g. 'kilogram meter**2 second**-1'
        """
        units = {}
        for unit, power in self.iter_base_or_scaled_units():
            units[unit] = power
        return 'Unit(%s)' % repr(units)

    # Performance
    _is_compatible_cache = {}

    def is_compatible(self, other):
        """
        Returns True if two Units share the same dimension.
        Returns False otherwise.
        """
        if self in Unit._is_compatible_cache:
            if other in Unit._is_compatible_cache[self]:
                return Unit._is_compatible_cache[self][other]
        if not is_unit(other):
            if self.is_dimensionless():
                return True
            else:
                return False
        self_dims = {}
        for dimension, exponent in self.iter_base_dimensions():
            self_dims[dimension] = exponent
        other_dims = {}
        for dimension, exponent in other.iter_base_dimensions():
            other_dims[dimension] = exponent
        if len(self_dims) != len(other_dims):
            result = False
        else:
            result = (self_dims == other_dims)
        if not self in Unit._is_compatible_cache:
            Unit._is_compatible_cache[self] = {}
        Unit._is_compatible_cache[self][other] = result
        return result

    _is_dimensionless_cache = {}

    def is_dimensionless(self):
        """Returns True if this Unit has no dimensions.
        Returns False otherwise.
        """
        if self in Unit._is_dimensionless_cache:
            return Unit._is_dimensionless_cache[self]
        for dimension, exponent in self.iter_base_dimensions():
            if exponent != 0:
                Unit._is_dimensionless_cache[self] = False
                return False
        Unit._is_dimensionless_cache[self] = True
        return True

    # Performance
    _conversion_factor_cache = {}

    def conversion_factor_to(self, other):
        """
        Returns conversion factor for computing all of the common dimensions
        between self and other from self base units to other base units.

        The two units need not share all of the same dimensions.  In case they
        do not, the conversion factor applies only to the BaseUnits of self
        that correspond to different BaseUnits in other.

        This method requires strict compatibility between the two units.
        """
        factor = 1.0
        if (self is other):
            return factor
        if self in Unit._conversion_factor_cache:
            if other in Unit._conversion_factor_cache[self]:
                return Unit._conversion_factor_cache[self][other]
        assert self.is_compatible(other)
        factor *= self.get_conversion_factor_to_base_units()
        factor /= other.get_conversion_factor_to_base_units()

        # Organize both units' base units by dimension.  Since so many conversion factors
        # are powers of ten, we accumulate them separately as an integer power to reduce
        # numerical error.

        canonical_units = {} # dimension: BaseUnit
        powers_of_ten = 0
        for unit, power in self.iter_all_base_units():
            d = unit.dimension
            if d in canonical_units:
                if unit != canonical_units[d]:
                    conversion = unit.conversion_factor_to(canonical_units[d])
                    log_conversion = math.log10(conversion)
                    if log_conversion == int(log_conversion):
                        powers_of_ten += power*int(log_conversion)
                    else:
                        factor *= conversion**power
            else:
                canonical_units[d] = unit
        for unit, power in other.iter_all_base_units():
            d = unit.dimension
            if d in canonical_units:
                if unit != canonical_units[d]:
                    conversion = unit.conversion_factor_to(canonical_units[d])
                    log_conversion = math.log10(conversion)
                    if log_conversion == int(log_conversion):
                        powers_of_ten -= power*int(log_conversion)
                    else:
                        factor /= conversion**power
            else:
                canonical_units[d] = unit
        factor *= 10**powers_of_ten
        if not self in Unit._conversion_factor_cache:
            Unit._conversion_factor_cache[self] = {}
        Unit._conversion_factor_cache[self][other] = factor
        return factor

    def in_unit_system(self, system):
        """
        Returns a new Unit with the same dimensions as this one, expressed in a particular unit system.

        Strips off any ScaledUnits in the Unit, leaving only BaseUnits.

        Parameters
        ----------
        system : a dictionary of (BaseDimension, BaseUnit) pairs
        """
        return system.express_unit(self)

    def get_symbol(self):
        """
        Returns a unit symbol (string) for this Unit, composed of its various
        BaseUnit symbols.  e.g. 'kg m**2 s**-1'
        """
        symbol = ""
        # emit positive exponents first
        pos = ""
        pos_count = 0
        for unit, power in self.iter_base_or_scaled_units():
            if power > 0:
                pos_count += 1
                if pos_count > 1: pos += " "
                pos += unit.symbol
                if power != 1.0:
                    pos += "**%g" % power
        # emit negative exponents second
        neg = ""
        neg_count = 0
        simple_denominator = True
        for unit, power in self.iter_base_or_scaled_units():
            if power < 0:
                neg_count += 1
                if neg_count > 1: neg += " "
                neg += unit.symbol
                if power != -1.0:
                    neg += "**%g" % -power
                    simple_denominator = False
        # Format of denominator depends on number of terms
        if 0 == neg_count:
            neg_string = ""
        elif 1 == neg_count and simple_denominator:
            neg_string = "/%s" % neg
        else:
            neg_string = "/(%s)" % neg
        if 0 == pos_count:
            pos_string = ""
        else:
            pos_string = pos
        if 0 == pos_count == neg_count:
            symbol = "dimensionless"
        else:
            symbol = "%s%s" % (pos_string, neg_string)
        return symbol

    def get_name(self):
        """
        Returns a unit name (string) for this Unit, composed of its various
        BaseUnit symbols.  e.g. 'kilogram meter**2 secon**-1'.
        """
        try:
            return self._name
        except AttributeError:
            pass
        # emit positive exponents first
        pos = ""
        pos_count = 0
        for unit, power in self.iter_base_or_scaled_units():
            if power > 0:
                pos_count += 1
                if pos_count > 1: pos += "*"
                pos += unit.name
                if power != 1.0:
                    pos += "**%g" % power
        # emit negative exponents second
        neg = ""
        neg_count = 0
        simple_denominator = True
        for unit, power in self.iter_base_or_scaled_units():
            if power < 0:
                neg_count += 1
                if neg_count > 1: neg += "*"
                neg += unit.name
                if power != -1.0:
                    neg += "**%g" % -power
                    simple_denominator = False
        # Format of denominator depends on number of terms
        if 0 == neg_count:
            neg_string = ""
        elif 1 == neg_count and simple_denominator:
            neg_string = "/%s" % neg
        else:
            neg_string = "/(%s)" % neg
        if 0 == pos_count:
            pos_string = ""
        else:
            pos_string = pos
        if 0 == pos_count == neg_count:
            name = "dimensionless"
        else:
            name = "%s%s" % (pos_string, neg_string)
        self._name = name
        return name


class ScaledUnit(object):
    """
    ScaledUnit is like a BaseUnit, but it is based on another Unit.

    ScaledUnit and BaseUnit are both used in the internals of Unit.  They
    should only be used during the construction of Units.
    """
    __array_priority__ = 100

    def __init__(self, factor, master, name, symbol):
        self.factor = factor
        # Convert to one base_unit per dimension
        base_units = {}
        for bu, exponent in master.iter_all_base_units():
            dim = bu.dimension
            if dim not in base_units:
                base_units[dim] = [bu, exponent]
            else:
                base_units[dim][1] += exponent
                self.factor *= base_units[dim][0].conversion_factor_to(bu)
        for sbu, exponent in master.iter_scaled_units():
            self.factor *= sbu.factor**exponent
        self.base_units = base_units
        self.master = master
        self.name = name
        self.symbol = symbol

    def __iter__(self):
        for dim in sorted(self.base_units.keys()):
            yield self.base_units[dim]

    def iter_base_units(self):
        for base_unit, exponent in self:
            yield(base_unit, exponent)

    def iter_base_dimensions(self):
        """
        Returns a sorted tuple of (BaseDimension, exponent) pairs, describing the dimension of this unit.
        """
        for base_unit, exponent in self:
            if exponent != 0:
                yield (base_unit.dimension, exponent)

    def get_dimension_tuple(self):
        """
        Returns a sorted tuple of (BaseDimension, exponent) pairs, that can be used as a dictionary key.
        """
        l = list(self.iter_base_dimensions())
        l.sort()
        return tuple(l)

    def get_conversion_factor_to_base_units(self):
        return self.factor

    def conversion_factor_to(self, other):
        # Create fake unit based on base units
        if self is other:
            return 1.0
        u = {}
        for base_unit, exponent in self.iter_base_units():
            u[base_unit] = exponent
        if isinstance(other, Unit):
            other_u = other
        else:
            other_u = Unit({other: 1.0})
        return self.factor * Unit(u).conversion_factor_to(other_u)

    def __lt__(self, other):
        """Compare two ScaledUnits.
        """
        return hash(self) < hash(other)

    def __str__(self):
        """Returns a string with the name of this ScaledUnit
        """
        return self.name

    def __repr__(self):
        """
        """
        base_units = ""
        for base_unit, power in self.iter_base_units():
            if len(base_units) > 0:
                base_units += ", "
            base_units += "%s: %d" % (base_unit, power)
        return "ScaledUnit(factor=" + repr(self.factor) + \
                ", master="+str(self.master)+", name=" + repr(self.name)\
                + ", symbol=" + repr(self.symbol) + ")"

class UnitSystem(object):
    """
    A complete system of units defining the *base* unit in each dimension

    Parameters
    ----------
    units : list
        List of base units from which to construct the unit system
    """
    def __init__(self, units):
        self.units = units
        self._unit_conversion_cache = {}
        # Create a set of base units to be used for dimension conversion
        base_units = {}
        for unit in self.units:
            for base_unit, exponent in unit.iter_base_units():
                d = base_unit.dimension
                if d not in base_units:
                    base_units[d] = base_unit
        self.base_units = base_units
        if not len(self.base_units) == len(self.units):
            raise ArithmeticError("UnitSystem must have same number of units as base dimensions")
        # self.dimensions is a dict of {BaseDimension: index}
        dimensions = sorted(base_units.keys())
        self.dimensions = {}
        for d in range(len(dimensions)):
            self.dimensions[dimensions[d]] = d
        # Create units->base units exponent matrix
        to_base_units = zeros(len(self.units))
        for m in range(len(self.units)):
            unit = self.units[m]
            for dim, power in unit.iter_base_dimensions():
                n = self.dimensions[dim]
                to_base_units[m][n] = power
        try:
            self.from_base_units = ~to_base_units
        except ArithmeticError as e:
        # for compatibility between python 2.5 and python 3.0,
        # try replacing line above with the following two lines:
        # except ArithmeticError:
        #     e=sys.exc_info[1]
            raise ArithmeticError("UnitSystem is not a valid basis set.  " + str(e))

    def __iter__(self):
        for unit in self.units:
            yield unit

    def __str__(self):
        """
        """
        result = "UnitSystem(["
        sep = ""
        for unit in self:
            result += sep
            result += str(unit)
            sep = ", "
        result += "])"
        return result

    def express_unit(self, old_unit):
        """
        """
        if old_unit in self._unit_conversion_cache:
            return self._unit_conversion_cache[old_unit]
        # First express unit in terms of base dimensions found in this unit system
        # (plus other dimensions not found)
        m = len(self.dimensions)
        base_dims = [0] * m
        other_dims = {}
        for dim, exponent in old_unit.iter_base_dimensions():
            if dim in self.dimensions:
                base_dims[self.dimensions[dim]] = exponent
            else:
                other_dims[dim] = exponent
        # Multiply by self.from_base_units to convert to unit system units
        u = MyMatrix([base_dims,]) * self.from_base_units
        new_unit = dimensionless
        for i in range(m):
            exponent = u[0][i]
            if exponent != 0:
                new_unit *= Unit({self.units[i]: exponent})
        if len(other_dims) > 0:
            # Find one base unit for each dimension
            found_dims = {}
            for base_unit, useless_exponent in old_unit.iter_all_base_units():
                dim = base_unit.dimension
                if dim not in other_dims:
                    continue # this dimension is in the unit system
                if dim in found_dims:
                    continue # already got a BaseUnit for this dimension
                found_dims[dim] = base_unit
                exponent = other_dims[dim]
                new_unit *= Unit({base_unit: exponent})
        self._unit_conversion_cache[old_unit] = new_unit
        return new_unit

def is_unit(x):
    """
    Returns True if x is a Unit, False otherwise.

    Examples
    --------
    >>> is_unit(16)
    False
    """
    return isinstance(x, Unit)

dimensionless = Unit({})

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
