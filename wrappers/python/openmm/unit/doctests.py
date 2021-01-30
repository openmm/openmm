#!/bin/env python

"""
Module openmm.unit.doctests

Lots of in-place doctests would no longer work after I rearranged
so that specific unit definitions are defined late.  So those tests
are here.

Examples

>>> furlong = BaseUnit(length_dimension, "furlong", "fur")

Examples

>>> furlong_base_unit = BaseUnit(length_dimension, "furlong", "fur")
>>> furlong_base_unit.define_conversion_factor_to(meter_base_unit, 201.16800)
>>> furlong_base_unit.conversion_factor_to(angstrom_base_unit)
2011680000000.0

Examples

>>> furlong_base_unit = BaseUnit(length_dimension, "furlong", "fur")
>>> furlong_base_unit.define_conversion_factor_to(meter_base_unit, 201.16800)

Examples

Some of these example test methods from Unit and Quantity

from unit.is_unit
>>> is_unit(meter)
True
>>> is_unit(5*meter)
False

>>> c = 1.0*calories
>>> c
Quantity(value=1.0, unit=calorie)
>>> print(calorie.conversion_factor_to(joule))
4.184
>>> print(joule.conversion_factor_to(calorie))
0.239005736138
>>> c.in_units_of(joules)
Quantity(value=4.184, unit=joule)
>>> j = 1.0*joules
>>> j
Quantity(value=1.0, unit=joule)
>>> j.in_units_of(calories)
Quantity(value=0.2390057361376673, unit=calorie)
>>> j/joules
1.0
>>> print(j/calories)
0.239005736138
>>> print(c/joules)
4.184
>>> c/calories
1.0
>>> c**2
Quantity(value=1.0, unit=calorie**2)
>>> (c**2).in_units_of(joule*joule)
Quantity(value=17.505856, unit=joule**2)

>>> ScaledUnit(1000.0, kelvin, "kilokelvin", "kK")
ScaledUnit(factor=1000.0, master=kelvin, name='kilokelvin', symbol='kK')

>>> str(ScaledUnit(1000.0, kelvin, "kilokelvin", "kK"))
'kilokelvin'

Examples

>>> meters > centimeters
True
>>> angstroms > centimeters
False

Examples

>>> print(meter / second)
meter/second
>>> print(meter / meter)
dimensionless

Heterogeneous units are not reduced unless they are in a quantity.
>>> print(meter / centimeter)
meter/centimeter

Examples

>>> meters_per_second = Unit({meter_base_unit: 1.0, second_base_unit: -1.0})
>>> print(meters_per_second)
meter/second

>>> us = UnitSystem([ScaledUnit(1.0, coulomb/second, "ampere", "A"), second_base_unit])
>>> print(us.express_unit(second))
second
>>> print(us.express_unit(coulomb/second))
ampere
>>> print(us.express_unit(coulomb))
second*ampere
>>> print(us.express_unit(meter/second))
meter/second

>>> us = UnitSystem([ScaledUnit(1.0, coulomb/second, "ampere", "A"), second_base_unit])
>>> print(us)
UnitSystem([ampere, second])

Examples

>>> meter.is_dimensionless()
False
>>> (meter/meter).is_dimensionless()
True

>>> print((meter*meter).sqrt())
meter
>>> meter.sqrt()
Traceback (most recent call last):
  ...
ArithmeticError: Exponents in Unit.sqrt() must be even.
>>> (meter*meter*meter).sqrt()
Traceback (most recent call last):
  ...
ArithmeticError: Exponents in Unit.sqrt() must be even.
>>> print((meter*meter/second/second).sqrt())
meter/second

Mixture of BaseUnits and ScaledUnits should cause no trouble:
>>> print(sqrt(kilogram*joule))
kilogram*meter/second
>>> print(sqrt(kilogram*calorie))
kilogram*meter/second

Examples

>>> newton.get_name()
'newton'
>>> meter.get_name()
'meter'

Examples

>>> newton.get_symbol()
'N'
>>> meter.get_symbol()
'm'

Examples

>>> print(angstrom.in_unit_system(si_unit_system))
meter
>>> print(angstrom.in_unit_system(cgs_unit_system))
centimeter
>>> print(angstrom.in_unit_system(md_unit_system))
nanometer
>>> u = meter/second**2
>>> print(u)
meter/(second**2)
>>> print(u.in_unit_system(si_unit_system))
meter/(second**2)
>>> print(u.in_unit_system(cgs_unit_system))
centimeter/(second**2)
>>> print(u.in_unit_system(md_unit_system))
nanometer/(picosecond**2)

Examples

>>> meter.is_compatible(centimeter)
True
>>> meter.is_compatible(meter)
True
>>> meter.is_compatible(kelvin)
False
>>> meter.is_compatible(meter/second)
False
>>> joule.is_compatible(calorie)
True

Examples

>>> meter.conversion_factor_to(centimeter)
100.0

>>> print((md_kilocalorie/mole/angstrom).conversion_factor_to(md_kilojoule/mole/nanometer))
41.84

Examples

>>> print(meter)
meter

>>> print(meter * second * second * kilogram)
kilogram*meter*second**2
>>> print(meter / second / second / kilogram)
meter/(kilogram*second**2)

Examples

>>> print(meter**3)
meter**3
>>> print(meter**3)
meter**3

>>> meter.get_conversion_factor_to_base_units()
1.0

Simple ScaledUnit in calorie
>>> print(calorie.get_conversion_factor_to_base_units())
4.184

Compound ScaledUnit in md_kilocalorie
>>> print(md_kilocalorie.get_conversion_factor_to_base_units())
4.184

calorie in a more complex unit
>>> print((md_kilocalorie/mole/angstrom).get_conversion_factor_to_base_units())
4.184

Examples

Create simple Quantities with either the multiply operator or the Quantity constructor.
>>> print(5 * centimeters)
5 cm
>>> print(Quantity(value=5, unit=centimeter))
5 cm
>>> print(Quantity(5, centimeter))
5 cm

Extract the underlying value using either division or the value_in_unit() method.
>>> i = 5 * centimeters
>>> print(i / millimeters)
50.0
>>> print(i.value_in_unit(millimeters))
50.0

Collections of numbers can also be used as values.
>>> s = [1,2,3] * centimeters
>>> print(s)
[1, 2, 3] cm
>>> print(s / millimeters)
[10.0, 20.0, 30.0]
>>> s2 = [[1,2,3],[4,5,6]] * centimeters
>>> print(s2)
[[1, 2, 3], [4, 5, 6]] cm
>>> print(s2 / millimeters)
[[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]
>>> s3 = [(1,2,3),(4,5,6)] * centimeters
>>> print(s3)
[(1, 2, 3), (4, 5, 6)] cm
>>> print(s3 / millimeters)
[(10.0, 20.0, 30.0), (40.0, 50.0, 60.0)]
>>> s4 = ((1,2,3),(4,5,6)) * centimeters
>>> print(s4)
((1, 2, 3), (4, 5, 6)) cm
>>> print(s4 / millimeters)
((10.0, 20.0, 30.0), (40.0, 50.0, 60.0))
>>> t = (1,2,3) * centimeters
>>> print(t)
(1, 2, 3) cm
>>> print(t / millimeters)
(10.0, 20.0, 30.0)

Numpy examples are commented out because not all systems have numpy installed
# >>> import numpy
# >>>
# >>> a = Quantity(numpy.array([1,2,3]), centimeters)
# >>> print(a)
# [1 2 3] cm
# >>> print(a / millimeters)
# [ 10.  20.  30.]
# >>>
# >>> a2 = Quantity(numpy.array([[1,2,3],[4,5,6]]), centimeters)
# >>> print(a2)
# [[1 2 3]
#  [4 5 6]] cm
# >>> print(a2 / millimeters)
# [[ 10.  20.  30.]
#  [ 40.  50.  60.]]

Addition, subtraction, multiplication, division, and powers of Quantities
exhibit correct dimensional analysis and unit conversion.
>>> x = 1.3 * meters
>>> y = 75.2 * centimeters
>>> print(x + y)
2.052 m
>>> print(x - y)
0.548 m
>>> print(x/y)
1.72872340426
>>> print(x*y)
0.9776 m**2

The following examples are derived from the C++ Boost.Units examples at
http://www.boost.org/doc/libs/1_37_0/doc/html/boost_units/Examples.html
>>>
>>> l = 2.0 * meters
>>>
>>> print(l + 2.0 * nanometers)
2.000000002 m
>>> print(2.0 * nanometers + l)
2000000002.0 nm
>>>
>>> print(l)
2.0 m
>>> print(l+l)
4.0 m
>>> print(l-l)
0.0 m
>>> print(l*l)
4.0 m**2
>>> print(l/l)
1.0
>>> print(l * meter)
2.0 m**2
>>> print(kilograms * (l/seconds) * (l/seconds))
4.0 kg m**2/(s**2)
>>> print(kilograms * (l/seconds)**2)
4.0 kg m**2/(s**2)
>>> print(l ** 3)
8.0 m**3
>>> print(l ** (3.0/2.0))
2.82842712475 m**1.5
>>> print(l ** 0.5)
1.41421356237 m**0.5
>>> print(l ** (2.0/3.0))
1.58740105197 m**0.666667
>>> # complex example
>>> l = (3.0 + 4.0j) * meters
>>> print(l)
(3+4j) m
>>> print(l+l)
(6+8j) m
>>> print(l-l)
0j m
>>> print(l*l)
(-7+24j) m**2
>>> # Numerical error yields tiny imaginary component of l/l on linux CentOS5
>>> err = abs(l/l - 1)
>>> assert err < 1e-8
>>> print(l * meter)
(3+4j) m**2
>>> print(kilograms * (l/seconds) * (l/seconds))
(-7+24j) kg m**2/(s**2)
>>> print(kilograms * (l/seconds)**2)
(-7+24j) kg m**2/(s**2)
>>> print(l ** 3)
(-117+44j) m**3
>>> print(l ** (3.0/2.0))
(2+11j) m**1.5
>>> print(l ** 0.5)
(2+1j) m**0.5
>>> print(l ** (2.0/3.0))
(2.38285471252+1.69466313833j) m**0.666667
>>> # kitchen sink example
... s1 = 2.0
>>> x1 = 2
>>> x2 = 4.0/3.0
>>> u1 = kilogram * meter / second**2
>>> u2 = u1 * meter
>>> q1 = 1.0*u1
>>> q2 = 2.0*u2
>>> print(s1)
2.0
>>> print(x1)
2
>>> print(x2)
1.33333333333
>>> print(u1)
kilogram*meter/(second**2)
>>> print(u2)
kilogram*meter**2/(second**2)
>>> print(q1)
1.0 kg m/(s**2)
>>> print(q2)
2.0 kg m**2/(s**2)
>>> print(u1*s1)
2.0 kg m/(s**2)
>>> print(s1*u1)
2.0 kg m/(s**2)
>>> print(u1/s1)
0.5 kg m/(s**2)
>>> print(s1/u1)
2.0 s**2/(kg m)
>>> print(u1*u1)
kilogram**2*meter**2/(second**4)
>>> print(u1/u1)
dimensionless
>>> print(u1*u2)
kilogram**2*meter**3/(second**4)
>>> print(u1/u2)
/meter
>>> print(u1**x1)
kilogram**2*meter**2/(second**4)
>>> print(u1**(1.0/x1))
kilogram**0.5*meter**0.5/second
>>> print(u1**x2)
kilogram**1.33333*meter**1.33333/(second**2.66667)
>>> print(u1**(1.0/x2))
kilogram**0.75*meter**0.75/(second**1.5)
>>> l1 = 1.0*meters
>>> l2 = 2.0*meters
>>> print(l1 == l2)
False
>>> print(l1 != l2)
True
>>> print(l1 <= l2)
True
>>> print(l1 < l2)
True
>>> print(l1 >= l2)
False
>>> print(l1 > l2)
False
>>>
>>> def work(f, dx):
...   return f * dx
...
>>> F = 1.0 * kilogram * meter / second**2
>>> dx = 1.0 * meter
>>> E = work(F, dx)
>>>
>>> print("F = ", F)
F =  1.0 kg m/(s**2)
>>> print("dx = ", dx)
dx =  1.0 m
>>>
>>> def idealGasLaw(P, V, T):
...     R = MOLAR_GAS_CONSTANT_R
...     print("P * V = ", P * V)
...     print("R * T = ", R * T)
...     return (P * V / (R * T)).in_units_of(mole)
...
>>> T = (273.0 + 37.0) * kelvin
>>> P = 1.01325e5 * pascals
>>> r = 0.5e-6 * meters
>>> V = 4.0/3.0 * 3.14159 * r**3
>>> n = idealGasLaw(P, V, T)
P * V =  5.3053601125e-14 m**3 Pa
R * T =  2577.48646608 J/mol
>>> R = MOLAR_GAS_CONSTANT_R
>>>
>>> print("r = ", r)
r =  5e-07 m
>>> print("P = ", P)
P =  101325.0 Pa
>>> print("V = ", V)
V =  5.23598333333e-19 m**3
>>> print("T = ", T)
T =  310.0 K
>>> print("n = ", n)
n =  2.05834644811e-17 mol
>>> print("R = ", R)
R =  8.31447247122 J/(K mol)
>>> print("E = ", E)
E =  1.0 kg m**2/(s**2)
>>> print("is_quantity(V) = ", is_quantity(V))
is_quantity(V) =  True
>>> print((1.0*radians) / degrees)
57.2957795131
>>> print((1.0*radians).in_units_of(degrees))
57.2957795131 deg
>>> print((1.0*angstroms).in_units_of(nanometers))
0.1 nm
>>>
>>> print((90*degrees)/radians)
1.57079632679
>>> print(sin(90*degrees))
1.0
>>> x = 90 * degrees
>>> x += 0.3 * radians
>>> print(x)
107.188733854 deg
>>> print(1 * nanometers > 1 * angstroms)
True
>>> print(1 * nanometers > 1 * degrees)
Traceback (most recent call last):
   ...
TypeError: Unit "degree" is not compatible with Unit "nanometer".
>>>
>>> x = 1.5 * nanometers
>>> print(x / meters)
1.5e-09
>>> x = 1.5 * angstroms
>>> print(x / meters)
1.5e-10
>>> print(x / nanometers)
0.15

Examples

>>> print(is_quantity(meters))
False
>>> print(is_quantity(2.3*meters))
True
>>> print(is_quantity(2.3))
False

Examples

>>> x = 100.0 * millimeter
>>> print(x.value_in_unit_system(si_unit_system))
0.1
>>> print(x.value_in_unit_system(cgs_unit_system))
10.0
>>> print(x.value_in_unit_system(md_unit_system))
100000000.0
>>>
>>> y = 20 * millimeters / millisecond**2
>>> print(y.value_in_unit_system(si_unit_system))
20000.0
>>> print(y.value_in_unit_system(cgs_unit_system))
2000000.0
>>> print(y.value_in_unit_system(md_unit_system))
2e-11
>>> eps = Quantity(1.0, md_kilocalorie/mole)
>>> epsQ = eps.value_in_unit_system(md_unit_system)
>>> print(epsQ)
4.184

Dimensionless quantities return their unmodified values.
>>> Quantity(5, dimensionless).value_in_unit_system(md_unit_system)
5

Examples

>>> x = 2.3*meters
>>> print(x.value_in_unit(centimeters))
230.0

Examples

>>> print(bool(2.3*meters))
True
>>> print(bool(0*meters))
False


Examples

>>> print(-(2.3*meters))
-2.3 m
>>> print(-(-2.3*meters))
2.3 m

Examples

>>> print(+(2.3*meters))
2.3 m

Examples

>>> print(abs(-2.3*meters))
2.3 m

>>> (9.0*meter*meter).sqrt()
Quantity(value=3.0, unit=meter)
>>> (9.0*meter).sqrt()
Traceback (most recent call last):
  ...
ArithmeticError: Exponents in Unit.sqrt() must be even.
>>> (9.0*meter*meter*meter).sqrt()
Traceback (most recent call last):
  ...
ArithmeticError: Exponents in Unit.sqrt() must be even.
>>> (9.0*meter*meter/second/second).sqrt()
Quantity(value=3.0, unit=meter/second)

Mixture of BaseUnits and ScaledUnits should cause no trouble:
>>> sqrt(1.0 * kilogram * joule)
Quantity(value=1.0, unit=kilogram*meter/second)
>>> sqrt(1.0 * kilogram * calorie)
Quantity(value=2.0454828280872954, unit=kilogram*meter/second)

Examples

>>> print((2.3*meters)**2)
5.29 m**2

Examples

>>> x = 4.2 * centimeters
>>> print(8.4 / x)
2.0 /cm


        Examples

        >>> x = 4.3 * meters
        >>> print(x/centimeters)
        430.0
        >>> print(x/seconds)
        4.3 m/s
        >>> x = [1,2,3]*centimeter
        >>> x/millimeter
        [10.0, 20.0, 30.0]


        Examples

        >>> x = 1.2*meters
        >>> print(5*x)
        6.0 m

        Examples

        >>> x = 1.2*meters
        >>> y = 72*centimeters
        >>> print(x*y)
        0.864 m**2
        >>> x = [1,2,3]*centimeter
        >>> x
        Quantity(value=[1, 2, 3], unit=centimeter)
        >>> x * meter
        Quantity(value=[100.0, 200.0, 300.0], unit=centimeter**2)

        >>> u = nanometer**2/angstrom**2
        >>> print(u)
        nanometer**2/(angstrom**2)
        >>> q = Quantity(2.0, u)
        >>> q
        Quantity(value=2.0, unit=nanometer**2/(angstrom**2))
        >>> "%.1f" % q.reduce_unit()
        '200.0'

        Examples

        >>> 1.2*meters < 72*centimeters
        False
        >>> meter is not None
        True
        >>> meter is None
        False

        Examples

        >>> print(1.2 * meters - 72 * centimeters)
        0.48 m

        Examples

        >>> print(1.2 * meters + 72 * centimeters)
        1.92 m

        Examples

        >>> print(repr(1.2*meter))
        Quantity(value=1.2, unit=meter)


        Examples

        >>> print(5.0 * nanometers)
        5.0 nm

        Examples

        >>> Quantity(5.0, meters)
        Quantity(value=5.0, unit=meter)

        >>> Quantity([1*angstrom,2*nanometer,3*angstrom])
        Quantity(value=[1, 20.0, 3], unit=angstrom)
        >>> Quantity((1,2,3))
        Quantity(value=(1, 2, 3), unit=dimensionless)
        >>> Quantity([1*angstrom,2*nanometer,3*angstrom])
        Quantity(value=[1, 20.0, 3], unit=angstrom)
        >>> Quantity([1*angstrom,2*nanometer,3*second])
        Traceback (most recent call last):
          ...
        TypeError: Unit "second" is not compatible with Unit "angstrom".
        >>> Quantity(5)
        Quantity(value=5, unit=dimensionless)

        Passing a unit to the constructor yields a Quantity with an empty list value.
        >>> Quantity(angstrom)
        Quantity(value=[], unit=angstrom)

        >>> Quantity(5*angstrom)
        Quantity(value=5, unit=angstrom)
        >>> Quantity(([1*angstrom,2*nanometer,3*angstrom], [1*angstrom,4*nanometer,3*angstrom]))
        Quantity(value=([1, 20.0, 3], [1, 40.0, 3]), unit=angstrom)
        >>> Quantity([])
        Quantity(value=[], unit=dimensionless)

        A simple scalar Quantity can be used as the unit argument.
        >>> Quantity(value=5.0, unit=100.0*meters)
        Quantity(value=500.0, unit=meter)


        Examples

        >>> x = 2.3*meters
        >>> y = x.in_units_of(centimeters)
        >>> print(y)
        230.0 cm

        >>> x = 2.3*meters
        >>> print(x.in_units_of(centimeters))
        230.0 cm
        >>> print(x.in_units_of(seconds))
        Traceback (most recent call last):
           ...
        TypeError: Unit "meter" is not compatible with Unit "second".

        Examples

        >>> x = 100.0 * millimeter
        >>> print(x)
        100.0 mm
        >>> print(x.in_unit_system(si_unit_system))
        0.1 m
        >>> print(x.in_unit_system(cgs_unit_system))
        10.0 cm
        >>> print(x.in_unit_system(md_unit_system))
        100000000.0 nm
        >>> y = 20 * millimeters / millisecond**2
        >>> print(y)
        20 mm/(ms**2)
        >>> print(y.in_unit_system(si_unit_system))
        20000.0 m/(s**2)
        >>> print(y.in_unit_system(cgs_unit_system))
        2000000.0 cm/(s**2)
        >>> print(y.in_unit_system(md_unit_system))
        2e-11 nm/(ps**2)

        Sometimes mixed internal units have caused trouble:
        >>> q = 1.0 * md_kilocalorie/mole/angstrom
        >>> print(q.in_units_of(md_kilojoule/mole/nanometer))
        41.84 kJ/(nm mol)

        Examples

        >>> class Foo:
        ...     def bar(self):
        ...         print("bar")
        ...
        >>> x = Foo()
        >>> x.bar()
        bar
        >>> y = x * nanometers
        >>> y.bar()
        bar

    Examples

    >>> print(meters * centimeters)
    centimeter*meter
    >>> print(meters * meters)
    meter**2
    >>> print(meter * meter )
    meter**2


    Examples

    >>> print(meter / 2)
    0.5 m

    Examples

    >>> define_prefixed_units(kelvin_base_unit, sys.modules["__main__"])
    >>> from __main__ import millikelvin
    >>> print(5.0 * millikelvin)
    5.0 mK


        Creating a new BaseUnit:
        >>> ms = milli * second_base_unit
        >>> ms
        BaseUnit(base_dim=BaseDimension("time"), name="millisecond", symbol="ms")
        >>> ms.conversion_factor_to(second_base_unit)
        0.001

        Creating a new ScaledUnit:
        >>> mC = milli * ScaledUnit(4.184, joule, "calorie", "cal")
        >>> mC
        ScaledUnit(factor=0.004184, master=joule, name='millicalorie', symbol='mcal')

        Creating a new Unit:
        >>> ms = milli * second
        >>> ms
        Unit({BaseUnit(base_dim=BaseDimension("time"), name="millisecond", symbol="ms"): 1.0})

        Don't try a Quantity though:
        >>> ms = milli * (1.0 * second)
        Traceback (most recent call last):
           ...
        TypeError: Unit prefix "milli" can only be applied to a Unit, BaseUnit, or ScaledUnit.

    Comparison of dimensionless quantities issue (fixed in svn 513)
    >>> x = Quantity(1.0, dimensionless)
    >>> y = Quantity(1.0, dimensionless)
    >>> assert not x is y
    >>> assert x == y


    Formatting of Quantities
    >>> x = 5.439999999 * picosecond
    >>> x
    Quantity(value=5.439999999, unit=picosecond)
    >>> x.format("%.3f")
    '5.440 ps'

    # Bug report Dec 17, 2009 from John Chodera
    # deepcopy of Quantity containing numpy array wrongly strips units
    >>> try:
    ...     import numpy
    ...     import copy
    ...     x = Quantity(numpy.zeros([2,3]), nanometer)
    ...     y = copy.deepcopy(x)
    ...     assert x[0][0] == y[0][0]
    ... except ImportError:
    ...     pass

    # Passing a string through Quantity constructor should return a string/dimensionless
    >>> x = Quantity("string").value_in_unit_system(md_unit_system)
    >>> assert x == "string"

# Trouble with complicated unit conversion factors
# Jan 29 1010 email from John Chodera
>>> p1 = 1.0 * atmosphere
>>> p2 = (1.0 * atmosphere).in_units_of(joule/nanometer**3)
>>> V = 2.4 * nanometer**3
>>> beta = 4.e-4 * mole/joule
>>> x1 = beta*p1*V
>>> # print(x1)
... y1 = x1 * AVOGADRO_CONSTANT_NA
>>> print(y1)
0.0585785776197

# Wrong answer is 5.85785776197e+25

>>> x2 = beta*p2*V
>>> # print(x2)
... y2 = x2 * AVOGADRO_CONSTANT_NA
>>> print(y2)
0.0585785776197
>>> assert( abs(y1 - y2) < 0.01)

# division of numpy arrays error
# April 2010, thanks to John Chodera for reporting
    >>> try:
    ...     import numpy
    ...     x = Quantity(numpy.array([1.,2.]), nanometer)
    ...     y = Quantity(numpy.array([3.,4.]), picosecond)
    ...     assert str(x/y) == '[ 0.33333333  0.5       ] nm/ps'
    ... except ImportError:
    ...     pass


# another numpy problem from retarded implementation of == operator
# Thanks to Kyle Beauchamp July 2010
    >>> try:
    ...    import numpy
    ...    from openmm.unit.quantity import _is_string
    ...    a = numpy.array([[1,2,3],[4,5,6]])
    ...    assert isinstance("", str)
    ...    assert _is_string("")
    ...    assert _is_string("t")
    ...    assert _is_string("test")
    ...    assert not _is_string(3)
    ...    assert not _is_string(a)
    ... except ImportError:
    ...    pass

"""

from __future__ import print_function

__author__ = "Christopher M. Bruns"
__version__ = "0.5"

# This unit code might be found in different packages...
# So use local import
from baseunit import BaseUnit
from standard_dimensions import *
from unit import is_unit, dimensionless
from quantity import Quantity, is_quantity, is_dimensionless
from unit_definitions import *
from unit_math import *
from constants import *

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    (failed, passed) = doctest.testmod(sys.modules[__name__])
    # For use in automated testing, return number of failed tests as exit code
    exit(failed)

