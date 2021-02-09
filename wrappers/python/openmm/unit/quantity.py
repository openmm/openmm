#!/bin/env python
"""
Module openmm.unit.quantity

Physical quantities with units, intended to produce similar functionality
to Boost.Units package in C++ (but with a runtime cost).
Uses similar API as Scientific.Physics.PhysicalQuantities
but different internals to satisfy our local requirements.
In particular, there is no underlying set of 'canonical' base
units, whereas in Scientific.Physics.PhysicalQuantities all
units are secretly in terms of SI units.  Also, it is easier
to add new fundamental dimensions to basedimension.  You
might want to make new dimensions for, say, "currency" or
"information".

Some features of this implementation:
  * Quantities are a combination of a value and a unit.  The value
    part can be any python type, including numbers, lists, numpy
    arrays, and anything else.  The unit part must be a openmm.unit.Unit.
  * Operations like adding incompatible units raises an error.
  * Multiplying or dividing units/quantities creates new units.
  * Users can create new Units and Dimensions, but most of the useful
    ones are predefined.
  * Conversion factors between units are applied transitively, so all
    possible conversions are available.
  * I want dimensioned Quantities that are compatible with numpy arrays,
    but do not necessarily require the python numpy package. In other
    words, Quantities can be based on either numpy arrays or on built in
    python types.
  * Units are NOT necessarily stored in terms of SI units internally.
    This is very important for me, because one important application
    area for us is at the molecular scale. Using SI units internally
    can lead to exponent overflow in commonly used molecular force
    calculations. Internally, all unit systems are equally fundamental
    in SimTK.

Two possible enhancements that have not been implemented are
  1) Include uncertainties with propagation of errors
  2) Incorporate offsets for celsius <-> kelvin conversion



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
import copy
from .standard_dimensions import *
from .unit import Unit, is_unit, dimensionless

class Quantity(object):
    """Physical quantity, such as 1.3 meters per second.

    Quantities contain both a value, such as 1.3; and a unit,
    such as 'meters per second'.

    Supported value types include:
      1 - numbers (float, int, long)
      2 - lists of numbers, e.g. [1,2,3]
      3 - tuples of numbers, e.g. (1,2,3)
            Note - unit conversions will cause tuples to be converted to lists
      4 - lists of tuples of numbers, lists of lists of ... etc. of numbers
      5 - numpy.arrays
    """
    __array_priority__ = 99

    def __init__(self, value=None, unit=None):
        """
        Create a new Quantity from a value and a unit.

        Parameters
         - value: (any type, usually a number) Measure of this quantity
         - unit: (Unit) the physical unit, e.g. openmm.unit.meters.
        """
        # When no unit is specified, bend over backwards to handle all one-argument possibilities
        if unit is None: # one argument version, copied from UList
            if is_unit(value):
                # Unit argument creates an empty list with that unit attached
                unit = value
                value = []
            elif is_quantity(value):
                # Ulist of a Quantity is just the Quantity itself
                unit = value.unit
                value = value._value
            elif _is_string(value):
                unit = dimensionless
            else:
                # Is value a container?
                is_container = True
                try:
                    _ = iter(value)
                except TypeError:
                    is_container = False
                if is_container:
                    if len(value) < 1:
                        unit = dimensionless
                    else:
                        first_item = next(iter(value))
                        # Avoid infinite recursion for string, because a one-character
                        # string is its own first element
                        try:
                            isstr = bool(value == first_item)
                        except ValueError:
                            # For numpy, value == first_item returns a numpy
                            # array of booleans, which cannot be evaluated for
                            # truthiness (a ValueError is raised). So in this
                            # case, we don't have a string
                            isstr = False
                        if isstr:
                            unit = dimensionless
                        else:
                            unit = Quantity(first_item).unit
                     # Notice that tuples, lists, and numpy.arrays can all be initialized with a list
                    new_container = Quantity([], unit)
                    for item in value:
                        new_container.append(Quantity(item)) # Strips off units into list new_container._value
                    # __class__ trick does not work for numpy.arrays
                    try:
                        import numpy
                        if isinstance(value, numpy.ndarray):
                            value = numpy.array(new_container._value)
                        else:
                            # delegate contruction to container class from list
                            value = value.__class__(new_container._value)
                    except ImportError:
                        # delegate contruction to container class from list
                        value = value.__class__(new_container._value)
                else:
                    # Non-Quantity, non container
                    # Wrap in a dimensionless Quantity
                    unit = dimensionless
        # Accept simple scalar quantities as units
        if is_quantity(unit):
            value = value * unit._value
            unit = unit.unit
        # Use empty list for unspecified values
        if value is None:
            value = []

        self._value = value
        self.unit = unit

    def __getstate__(self):
        state = dict()
        state['_value'] = self._value
        state['unit'] = self.unit
        return state

    def __setstate__(self, state):
        self._value = state['_value']
        self.unit = state['unit']
        return

    def __copy__(self):
        """
        Shallow copy produces a new Quantity with the shallow copy of value and the same unit.
        Because we want copy operations to work just the same way they would on the underlying value.
        """
        return Quantity(copy.copy(self._value), self.unit)

    def __deepcopy__(self, memo):
        """
        Deep copy produces a new Quantity with a deep copy of the value, and the same unit.
        Because we want copy operations to work just the same way they would on the underlying value.
        """
        return Quantity(copy.deepcopy(self._value, memo), self.unit)

    def __getattr__(self, attribute):
        """
        Delegate unrecognized attribute calls to the underlying value type.
        """
        ret_val = getattr(self._value, attribute)
        return ret_val

    def __str__(self):
        """Printable string version of this Quantity.

        Returns a string consisting of quantity number followed by unit abbreviation.
        """
        return str(self._value) + ' ' + str(self.unit.get_symbol())

    def __repr__(self):
        """
        """
        return (Quantity.__name__ + '(value=' + repr(self._value) + ', unit=' +
                str(self.unit) + ')')

    def format(self, format_spec):
        return format_spec % self._value + ' ' + str(self.unit.get_symbol())

    def __add__(self, other):
        """Add two Quantities.

        Only Quantities with the same dimensions (e.g. length)
        can be added.  Raises TypeError otherwise.

        Parameters
         - self: left hand member of sum
         - other: right hand member of sum

        Returns a new Quantity that is the sum of the two arguments.
        """
        # can only add using like units
        if not self.unit.is_compatible(other.unit):
            raise TypeError('Cannot add two quantities with incompatible units "%s" and "%s".' % (self.unit, other.unit))
        value = self._value + other.value_in_unit(self.unit)
        unit = self.unit
        return Quantity(value, unit)

    def __sub__(self, other):
        """Subtract two Quantities.

        Only Quantities with the same dimensions (e.g. length)
        can be subtracted.  Raises TypeError otherwise.

        Parameters
         - self: left hand member (a) of a - b.
         - other: right hand member (b) of a - b.

        Returns a new Quantity that is the difference of the two arguments.
        """
        if not self.unit.is_compatible(other.unit):
            raise TypeError('Cannot subtract two quantities with incompatible units "%s" and "%s".' % (self.unit, other.unit))
        value = self._value - other.value_in_unit(self.unit)
        unit = self.unit
        return Quantity(value, unit)

    def __eq__(self, other):
        """
        """
        if not is_quantity(other):
            return False
        if not self.unit.is_compatible(other.unit):
            return False
        return self.value_in_unit(other.unit) == other._value

    def __ne__(self, other):
        """
        """
        return not self == other

    def __lt__(self, other):
        """Compares two quantities.

        Raises TypeError if the Quantities are of different dimension (e.g. length vs. mass)

        Returns True if self < other, False otherwise.
        """
        return self._value < other.value_in_unit(self.unit)

    def __ge__(self, other):
        return self._value >= (other.value_in_unit(self.unit))
    def __gt__(self, other):
        return self._value > (other.value_in_unit(self.unit))
    def __le__(self, other):
        return self._value <= (other.value_in_unit(self.unit))
    def __lt__(self, other):
        return self._value < (other.value_in_unit(self.unit))

    __hash__ = None

    _reduce_cache = {}

    def reduce_unit(self, guide_unit=None):
        """
        Combine similar component units and scale, to form an
        equal Quantity in simpler units.

        Returns underlying value type if unit is dimensionless.
        """
        key = (self.unit, guide_unit)
        if key in Quantity._reduce_cache:
            (unit, value_factor) = Quantity._reduce_cache[key]
        else:
            value_factor = 1.0
            canonical_units = {} # dict of dimensionTuple: (Base/ScaledUnit, exponent)
            # Bias result toward guide units
            if guide_unit is not None:
                for u, exponent in guide_unit.iter_base_or_scaled_units():
                    d = u.get_dimension_tuple()
                    if d not in canonical_units:
                        canonical_units[d] = [u, 0]
            for u, exponent in self.unit.iter_base_or_scaled_units():
                d = u.get_dimension_tuple()
                # Take first unit found in a dimension as canonical
                if d not in canonical_units:
                    canonical_units[d] = [u, exponent]
                else:
                    value_factor *= (u.conversion_factor_to(canonical_units[d][0])**exponent)
                    canonical_units[d][1] += exponent
            new_base_units = {}
            for d in canonical_units:
                u, exponent = canonical_units[d]
                if exponent != 0:
                    assert u not in new_base_units
                    new_base_units[u] = exponent
            # Create new unit
            if len(new_base_units) == 0:
                unit = dimensionless
            else:
                unit = Unit(new_base_units)
            # There might be a factor due to unit conversion, even though unit is dimensionless
            # e.g. suppose unit is meter/centimeter
            if unit.is_dimensionless():
                unit_factor = unit.conversion_factor_to(dimensionless)
                if unit_factor != 1.0:
                    value_factor *= unit_factor
                    # print "value_factor = %s" % value_factor
                unit = dimensionless
            Quantity._reduce_cache[key] = (unit, value_factor)
        # Create Quantity, then scale (in case value is a container)
        # That's why we don't just scale the value.
        result = Quantity(self._value, unit)
        if value_factor != 1.0:
            # __mul__ strips off dimensionless, if appropriate
            result = result * value_factor
        if unit.is_dimensionless():
            assert unit is dimensionless # should have been set earlier in this method
            if is_quantity(result):
                result = copy.deepcopy(result._value)
        return result

    def __mul__(self, other):
        """Multiply a quantity by another object

        Returns a new Quantity that is the product of the self * other,
        unless the resulting unit is dimensionless, in which case the
        underlying value type is returned, instead of a Quantity.
        """
        if is_unit(other):
            # print "quantity * unit"
            # Many other mul/div operations delegate to here because I was debugging
            # a dimensionless unit conversion problem, which I ended up fixing within
            # the reduce_unit() method.
            unit = self.unit * other
            return Quantity(self._value, unit).reduce_unit(self.unit)
        elif is_quantity(other):
            # print "quantity * quantity"
            # Situations where the units cancel can result in scale factors from the unit cancellation.
            # To simplify things, delegate Quantity * Quantity to (Quantity * scalar) * unit
            return (self * other._value) * other.unit
        else:
            # print "quantity * scalar"
            return self._change_units_with_factor(self.unit, other, post_multiply=False)

    # value type might not be commutative for multiplication
    def __rmul__(self, other):
        """Multiply a scalar by a Quantity

        Returns a new Quantity with the same units as self, but with the value
        multiplied by other.
        """
        if is_unit(other):
            raise NotImplementedError('programmer is surprised __rmul__ was called instead of __mul__')
            # print "R unit * quantity"
        elif is_quantity(other):
            # print "R quantity * quantity"
            raise NotImplementedError('programmer is surprised __rmul__ was called instead of __mul__')
        else:
            # print "scalar * quantity"
            return self._change_units_with_factor(self.unit, other, post_multiply=True)
            # return Quantity(other * self._value, self.unit)

    def __truediv__(self, other):
        """Divide a Quantity by another object

        Returns a new Quantity, unless the resulting unit type is dimensionless,
        in which case the underlying value type is returned.
        """
        if is_unit(other):
            # print "quantity / unit"
            return self * pow(other, -1.0)
            # unit = self.unit / other
            # return Quantity(self._value, unit).reduce_unit(self.unit)
        elif is_quantity(other):
            # print "quantity / quantity"
            # Delegate quantity/quantity to (quantity/scalar)/unit
            return (self/other._value) / other.unit
        else:
            # print "quantity / scalar"
            return self * pow(other, -1.0)
            # return Quantity(self._value / other, self.unit)

    __div__ = __truediv__

    def __rtruediv__(self, other):
        """Divide a scalar by a quantity.

        Returns a new Quantity.  The resulting units are the inverse of the self argument units.
        """
        if is_unit(other):
            # print "R unit / quantity"
            raise NotImplementedError('programmer is surprised __rtruediv__ was called instead of __truediv__')
        elif is_quantity(other):
            raise NotImplementedError('programmer is surprised __rtruediv__ was called instead of __truediv__')
        else:
            # print "R scalar / quantity"
            return other * pow(self, -1.0)
            # return Quantity(other / self._value, pow(self.unit, -1.0))

    __rdiv__ = __rtruediv__

    def __pow__(self, exponent):
        """Raise a Quantity to a power.

        Generally both the value and the unit of the Quantity are affected by this operation.

        Returns a new Quantity equal to self**exponent.
        """
        return Quantity(pow(self._value, exponent), pow(self.unit, exponent))

    def sqrt(self):
        """
        Returns square root of a Quantity.

        Raises ArithmeticError if component exponents are not even.
        This behavior can be changed if you present a reasonable real life case to me.
        """
        # There might be a conversion factor from taking the square root of the unit
        new_value = math.sqrt(self._value)
        new_unit = self.unit.sqrt()
        unit_factor = self.unit.conversion_factor_to(new_unit*new_unit)
        if unit_factor != 1.0:
            new_value *= math.sqrt(unit_factor)
        return Quantity(value=new_value, unit=new_unit)

    def sum(self, *args, **kwargs):
        """
        Computes the sum of a sequence, with the result having the same unit as
        the current sequence.

        If the value is not iterable, it raises a TypeError (same behavior as if
        you tried to iterate over, for instance, an integer).

        This function can take as arguments any arguments recognized by
        `numpy.sum`. If arguments are passed to a non-numpy array, a TypeError
        is raised
        """
        try:
            # This will be much faster for numpy arrays
            mysum = self._value.sum(*args, **kwargs)
        except AttributeError:
            if args or kwargs:
                raise TypeError('Unsupported arguments for Quantity.sum')
            if len(self._value) == 0:
                mysum = 0
            else:
                mysum = self._value[0]
                for i in range(1, len(self._value)):
                    mysum += self._value[i]
        return Quantity(mysum, self.unit)

    def mean(self, *args, **kwargs):
        """
        Computes the mean of a sequence, with the result having the same unit as
        the current sequence.

        If the value is not iterable, it raises a TypeError

        This function can take as arguments any arguments recognized by
        `numpy.mean`. If arguments are passed to a non-numpy array, a TypeError
        is raised
        """
        try:
            # Faster for numpy arrays
            mean = self._value.mean(*args, **kwargs)
        except AttributeError:
            if args or kwargs:
                raise TypeError('Unsupported arguments for Quantity.mean')
            mean = (self.sum() / len(self._value))._value
        return Quantity(mean, self.unit)

    def std(self, *args, **kwargs):
        """
        Computes the square root of the variance of a sequence, with the result
        having the same unit as the current sequence.

        If the value is not iterable, it raises a TypeError

        This function can take as arguments any arguments recognized by
        `numpy.std`. If arguments are passed to a non-numpy array, a TypeError
        is raised
        """
        try:
            # Faster for numpy arrays
            std = self._value.std(*args, **kwargs)
        except AttributeError:
            if args or kwargs:
                raise TypeError('Unsupported arguments for Quantity.std')
            mean = self.mean()._value
            var = 0
            for val in self._value:
                res = mean - val
                var += res * res
            var /= len(self._value)
            std = math.sqrt(var)
        return Quantity(std, self.unit)

    def max(self, *args, **kwargs):
        """
        Computes the maximum value of the sequence, with the result having the
        same unit as the current sequence.

        If the value is not iterable, it raises a TypeError

        This function can take as arguments any arguments recognized by
        `numpy.max`. If arguments are passed to a non-numpy array, a TypeError
        is raised
        """
        try:
            # Faster for numpy arrays
            mymax = self._value.max(*args, **kwargs)
        except AttributeError:
            if args or kwargs:
                raise TypeError('Unsupported arguments for Quantity.max')
            mymax = max(self._value)
        return Quantity(mymax, self.unit)

    def min(self, *args, **kwargs):
        """
        Computes the minimum value of the sequence, with the result having the
        same unit as the current sequence.

        If the value is not iterable, it raises a TypeError

        This function can take as arguments any arguments recognized by
        `numpy.min`. If arguments are passed to a non-numpy array, a TypeError
        is raised
        """
        try:
            # Faster for numpy arrays
            mymin = self._value.min(*args, **kwargs)
        except AttributeError:
            if args or kwargs:
                raise TypeError('Unsupported arguments for Quantity.min')
            mymin = min(self._value)
        return Quantity(mymin, self.unit)

    def reshape(self, shape, order='C'):
        """
        Same as numpy.ndarray.reshape, except the result is a Quantity with the
        same units as the current object rather than a plain numpy.ndarray
        """
        try:
            return Quantity(self._value.reshape(shape, order=order), self.unit)
        except AttributeError:
            raise AttributeError('Only numpy array Quantity objects can be '
                                 'reshaped')

    def __abs__(self):
        """
        Return absolute value of a Quantity.

        The unit is unchanged.  A negative value of self will result in a positive value
        in the result.
        """
        return Quantity(abs(self._value), self.unit)

    def __pos__(self):
        """
        Returns a reference to self.
        """
        return Quantity(+(self._value), self.unit)

    def __neg__(self):
        """Negate a Quantity.

        Returns a new Quantity with a different sign on the value.
        """
        return Quantity(-(self._value), self.unit)

    def __nonzero__(self):
        """Returns True if value underlying Quantity is zero, False otherwise.
        """
        return bool(self._value)

    def __bool__(self):
        return bool(self._value)

    def __complex__(self):
        return Quantity(complex(self._value), self.unit)
    def __float__(self):
        return Quantity(float(self._value), self.unit)
    def __int__(self):
        return Quantity(int(self._value), self.unit)
    def __long__(self):
        return Quantity(int(self._value), self.unit)

    def value_in_unit(self, unit):
        """
        Returns underlying value, in the specified units.
        """
        val = self.in_units_of(unit)
        if is_quantity(val):
            return val._value
        else: # naked dimensionless
            return val

    def value_in_unit_system(self, system):
        """
        Returns the underlying value type, after conversion to a particular unit system.
        """
        result = self.in_unit_system(system)
        if is_quantity(result):
            return result._value
        else:
            return result # dimensionless

    def in_unit_system(self, system):
        """
        Returns a new Quantity equal to this one, expressed in a particular unit system.
        """
        new_units = system.express_unit(self.unit)
        f = self.unit.conversion_factor_to(new_units)
        return self._change_units_with_factor(new_units, f)

    def in_units_of(self, other_unit):
        """
        Returns an equal Quantity expressed in different units.

        If the units are the same as those in self, a reference to self is returned.
        Raises a TypeError if the new unit is not compatible with the original unit.

        The post_multiply argument is used in case the multiplication operation is not commutative.
          i.e. result = factor * value when post_multiply is False
          and  result = value * factor when post_multiply is True
        """
        if not self.unit.is_compatible(other_unit):
            raise TypeError('Unit "%s" is not compatible with Unit "%s".' % (self.unit, other_unit))
        f = self.unit.conversion_factor_to(other_unit)
        return self._change_units_with_factor(other_unit, f)

    def _change_units_with_factor(self, new_unit, factor, post_multiply=True):
        # numpy arrays cannot be compared with 1.0, so just "try"
        factor_is_identity = False
        try:
            if (factor == 1.0):
                factor_is_identity = True
        except ValueError:
            pass
        if factor_is_identity:
            # No multiplication required
            result = Quantity(copy.deepcopy(self._value), new_unit)
        else:
            try:
                # multiply operator, if it exists, is preferred
                if post_multiply:
                    value = self._value * factor # works for number, numpy.array, or vec3, e.g.
                else:
                    value = factor * self._value # works for number, numpy.array, or vec3, e.g.
                result = Quantity(value, new_unit)
            except TypeError:
                value = copy.deepcopy(self._value)
                result = Quantity(self._scale_sequence(value, factor, post_multiply), new_unit)
        if (new_unit.is_dimensionless()):
            return result._value
        else:
            return result

    def _scale_sequence(self, value, factor, post_multiply):
        try:
            if post_multiply:
                value = value*factor
            else:
                value = factor*value
        except TypeError:
            try:
                if post_multiply:
                    if isinstance(value, tuple):
                        value = tuple([x*factor for x in value])
                    else:
                        for i in range(len(value)):
                            value[i] = value[i]*factor
                else:
                    if isinstance(value, tuple):
                        value = tuple([factor*x for x in value])
                    else:
                        for i in range(len(value)):
                            value[i] = factor*value[i]
            except TypeError:
                if isinstance(value, tuple):
                    value = tuple([self._scale_sequence(x, factor, post_multiply) for x in value])
                else:
                    for i in range(len(value)):
                        value[i] = self._scale_sequence(value[i], factor, post_multiply)
        return value



    ####################################
    ### Sequence methods of Quantity ###
    ###  in case value is a sequence ###
    ####################################

    def __len__(self):
        """
        Return size of internal value type.
        """
        return len(self._value)

    def __getitem__(self, key):
        """
        Keep the same units on contained elements.
        """
        assert not is_quantity(self._value[key])
        return Quantity(self._value[key], self.unit)

    def __setitem__(self, key, value):
        # Delegate slices to one-at-a time ___setitem___
        if isinstance(key, slice): # slice
            indices = key.indices(len(self))
            for value_idx, self_idx in enumerate(range(*indices)):
                self[self_idx] = value[value_idx]
        else: # single index
            # Check unit compatibility
            if self.unit.is_dimensionless() and is_dimensionless(value):
                pass # OK
            elif not self.unit.is_compatible(value.unit):
                raise TypeError('Unit "%s" is not compatible with Unit "%s".' % (self.unit, value.unit))
            self._value[key] = value / self.unit
            assert not is_quantity(self._value[key])

    def __delitem__(self, key):
        del(self._value[key])

    def __contains__(self, item):
        return self._value.__contains__(item.value_in_unit(self.unit))

    def __iter__(self):
        for item in self._value:
            yield Quantity(item, self.unit)

    def count(self, item):
        return self._value.count(item.value_in_unit(self.unit))
    def index(self, item):
        return self._value.index(item.value_in_unit(self.unit))
    def append(self, item):
        if is_quantity(item):
            return self._value.append(item.value_in_unit(self.unit))
        elif is_dimensionless(self.unit):
            return self._value.append(item)
        else:
            raise TypeError("Cannot append item without units into list with units")
    def extend(self, rhs):
        self._value.extend(rhs.value_in_unit(self.unit))
    def insert(self, index, item):
        self._value.insert(index, item.value_in_unit(self.unit))
    def remove(self, item):
        self._value.remove(item)
    def pop(self, *args):
        return self._value.pop(*args) * self.unit
    # list.reverse will automatically delegate correctly
    # list.sort with no arguments will delegate correctly
    # list.sort with a comparison function cannot be done correctly


def is_quantity(x):
    """
    Returns True if x is a Quantity, False otherwise.
    """
    return isinstance(x, Quantity)

def is_dimensionless(x):
    """
    """
    if is_unit(x):
        return x.is_dimensionless()
    elif is_quantity(x):
        return x.unit.is_dimensionless()
    else:
        # everything else in the universe is dimensionless
        return True

# Strings can cause trouble
# as can any container that has infinite levels of containment
def _is_string(x):
     # step 1) String is always a container
     # and its contents are themselves containers.
     if isinstance(x, str):
         return True
     try:
         first_item = next(iter(x))
         inner_item = next(iter(first_item))
         if first_item is inner_item:
             return True
         else:
             return False
     except TypeError:
         return False
     except StopIteration:
         return False

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
