"""
Tests the functionality in the simtk.unit package.
"""
from __future__ import division

from simtk import unit as u
import copy
import math
import unittest
try:
    import numpy as np
except ImportError:
    np = None
try:
    from itertools import izip as zip
except ImportError:
    pass # Python 3... zip _is_ izip

class QuantityTestCase(unittest.TestCase):

    def assertAlmostEqualQuantities(self, item1, item2, places=6):
        try:
            val1 = item1.value_in_unit(item1.unit)
            val2 = item2.value_in_unit(item1.unit)
        except TypeError:
            raise self.failureException('Incompatible units %s and %s' %
                                        (item1.unit, item2.unit))
        try:
            if len(val1) != len(val2):
                raise self.failureException('collections are different lengths')
            for x, y in zip(val1, val2):
                self.assertAlmostEqual(x, y, places=places)
        except TypeError:
            self.assertAlmostEqual(val1, val2, places=places)

class TestUnits(QuantityTestCase):

    def testBaseUnit(self):
        """ Tests the creation of a base unit furlong """
        furlong = u.BaseUnit(u.length_dimension, "furlong", "fur")
        furlong.define_conversion_factor_to(u.meter_base_unit, 201.168)
        self.assertEqual(furlong.conversion_factor_to(u.angstrom_base_unit),
                         201.168e10)

    def testIsUnit(self):
        """ Tests the unit.is_unit introspective function """
        self.assertTrue(u.is_unit(u.meters))
        self.assertFalse(u.is_unit(None))
        self.assertFalse(u.is_unit(10))
        # Distinguish between "units" and "quantities"
        self.assertFalse(u.is_unit(1*u.meters))

    def testCalorieConversion(self):
        """ Tests unit conversion with calories to joules """
        c = 1.0 * u.calories
        j = 1.0 * u.joules
        j2c = j.in_units_of(u.calories)
        self.assertNotEqual(c, 1.0) # Units do not equal scalars!
        self.assertIs(c.unit, u.calories)
        self.assertEqual(c._value, 1)
        self.assertEqual(u.calorie.conversion_factor_to(u.joule), 4.184)
        self.assertEqual(u.joule.conversion_factor_to(u.calorie), 1/4.184)
        self.assertIs(j2c.unit, u.calories)
        self.assertEqual(j2c._value, 1/4.184)
        self.assertEqual(c / u.calories, 1.0) # Dimensionless now
        self.assertEqual(j / u.calories, j2c._value)
        c2 = c**2
        cc = c * c
        self.assertEqual(cc, c2)
        self.assertTrue(c2.unit, u.calories**2)
        self.assertEqual(c2.in_units_of(u.joules**2), (4.184*u.joules)**2)

    def testScaledUnit(self):
        """ Tests ScaledUnit class with kilokelvins """
        kK = u.ScaledUnit(1000.0, u.kelvin, "kilokelvin", "kK")
        self.assertIs(kK.master, u.kelvin)
        self.assertEqual(kK.factor, 1000)
        self.assertEqual(str(kK), 'kilokelvin')
        self.assertEqual(kK.name, 'kilokelvin')
        self.assertEqual(kK.symbol, 'kK')

    def testComparisonOperators(self):
        """ Tests unit comparison operators """
        self.assertGreater(u.meters, u.centimeters)
        self.assertLess(u.angstroms, u.centimeters)

    def testUnitDivision(self):
        """ Tests the division of units and is_dimensionless """
        mps = u.meter / u.second
        mpm = u.meter / u.meter
        mpc = u.meter / u.centimeter
        self.assertTrue(u.is_unit(mps))
        self.assertTrue(u.is_unit(mpm))
        self.assertTrue(u.is_unit(mpc))
        self.assertFalse(mps.is_dimensionless())
        self.assertTrue(mpm.is_dimensionless())
        self.assertTrue(mpc.is_dimensionless())

    def testCompositeUnits(self):
        """ Tests the creation of a composite unit """
        mps = u.Unit({u.meter_base_unit : 1.0, u.second_base_unit : -1.0})
        self.assertTrue(u.is_unit(mps))
        self.assertEqual(str(mps), 'meter/second')

    def testUnitSystem(self):
        """ Tests the creation of a UnitSystem and its behavior """
        us = u.UnitSystem([u.ScaledUnit(1.0, u.coulomb/u.second, 'ampere', 'A'),
                           u.second_base_unit])
        self.assertEqual(us.express_unit(u.second), u.second)
        self.assertEqual(us.express_unit(u.hour), u.second)
        self.assertNotEqual(u.hour, u.second)
        self.assertEqual(us.express_unit(u.coulomb/u.second), u.ampere)
        self.assertEqual(us.express_unit(u.meter/u.second), u.meter/u.second)
        self.assertEqual(us.express_unit(u.kilometer/u.hour),
                         u.kilometer/u.second)
        x = 100 * u.millimeters
        self.assertEqual(x.value_in_unit_system(u.si_unit_system), 0.1)
        self.assertEqual(x.value_in_unit_system(u.cgs_unit_system), 10.0)
        self.assertAlmostEqual(x.value_in_unit_system(u.md_unit_system), 1e8)
        y = 20 * u.millimeters / u.millisecond**2
        self.assertAlmostEqual(y.value_in_unit_system(u.si_unit_system), 2e4)
        self.assertAlmostEqual(y.value_in_unit_system(u.cgs_unit_system), 2e6)
        self.assertAlmostEqual(y.value_in_unit_system(u.md_unit_system), 2e-11,
                               places=18)
        kcal = 1 * u.md_kilocalorie / u.mole
        self.assertEqual(kcal.value_in_unit_system(u.md_unit_system), 4.184)

    def testUnitSqrt(self):
        """ Tests taking the square root of units """
        self.assertEqual((u.meter * u.meter).sqrt(), u.meter)
        self.assertEqual((u.meter**4).sqrt(), u.meter**2)

    def testQuantitySqrt(self):
        """ Tests taking the square root of quantities """
        self.assertEqual((9.0*u.meters**2).sqrt(), 3.0*u.meters)
        self.assertEqual((9.0*u.meters**2/u.second**2).sqrt(),
                         3.0*u.meters/u.second)

    def testUnitBadSqrt(self):
        """ Tests that illegal sqrt calls on incompatible units fails """
        mps2 = u.meters/u.second**2
        self.assertRaises(ArithmeticError, lambda: u.meter.sqrt())
        self.assertRaises(ArithmeticError, lambda: (u.meters**3).sqrt())
        self.assertRaises(ArithmeticError, lambda: mps2.sqrt())

    def testQuantityBadSqrt(self):
        " Tests that taking sqrt of Quantities with incompatible units fails "
        self.assertRaises(ArithmeticError, lambda: (9.0*u.meters).sqrt())
        self.assertRaises(ArithmeticError, lambda: (9.0*u.meters**3).sqrt())

    def testBaseScaleMix(self):
        """ Test mixing of BaseUnit and ScaledUnit instances """
        kgj = u.kilogram * u.joule
        self.assertTrue(u.is_unit(kgj))
        self.assertEqual(str(kgj.sqrt()), 'kilogram*meter/second')

    def testGetUnitAttributes(self):
        """ Tests the unit attribute `getters' """
        self.assertEqual(u.newton.get_name(), 'newton')
        self.assertEqual(u.newtons.get_name(), 'newton')
        self.assertEqual(u.newton.get_symbol(), 'N')
        self.assertEqual(u.ampere.get_symbol(), 'A')
        self.assertEqual(u.meter.get_name(), 'meter')
        self.assertEqual(u.meter.get_symbol(), 'm')

    def testPresetUnitSystems(self):
        """ Tests some of the pre-set UnitSystem's """
        self.assertEqual(u.angstrom.in_unit_system(u.si_unit_system), u.meter)
        self.assertEqual(u.angstrom.in_unit_system(u.cgs_unit_system),
                         u.centimeter)
        self.assertEqual(u.angstrom.in_unit_system(u.md_unit_system),
                         u.nanometer)
        mps = u.meter / u.second**2
        self.assertEqual(str(mps), 'meter/(second**2)')
        self.assertEqual(mps.in_unit_system(u.si_unit_system), mps)
        self.assertEqual(mps.in_unit_system(u.cgs_unit_system),
                         u.centimeter/u.second**2)
        self.assertEqual(mps.in_unit_system(u.md_unit_system),
                         u.nanometer/u.picosecond**2)

    def testIsCompatible(self):
        """ Tests the is_compatible attribute of units """
        self.assertTrue(u.meter.is_compatible(u.centimeter))
        self.assertTrue(u.centimeter.is_compatible(u.meter))
        self.assertTrue(u.meter.is_compatible(u.meter))
        self.assertTrue(u.joule.is_compatible(u.calorie))
        self.assertFalse(u.meter.is_compatible(u.kelvin))
        self.assertFalse(u.meter.is_compatible(u.meter/u.second))
        self.assertFalse(u.meter.is_compatible(u.joule))

    def testConversionFactorTo(self):
        """ Tests the "conversion_factor_to" attribute of Units """
        self.assertEqual(u.meter.conversion_factor_to(u.centimeter), 100)
        self.assertEqual(u.kilocalorie.conversion_factor_to(u.joule), 4184)
        self.assertEqual(u.calorie.conversion_factor_to(u.kilojoule), 4.184e-3)
        self.assertEqual((u.kilocalorie/u.mole/u.angstrom).conversion_factor_to(
                                u.kilojoule/u.mole/u.nanometer), 41.84)

    def testUnitString(self):
        """ Test the Unit->str casting functionality """
        # Always alphabetical order
        self.assertEqual(str(u.meter * u.second * u.second * u.kilogram),
                         'kilogram*meter*second**2')
        self.assertEqual(str(u.meter / u.second / u.second / u.kilogram),
                         'meter/(kilogram*second**2)')
        self.assertEqual(str(u.meter**3), 'meter**3')

    def testToBaseUnit(self):
        """ Tests the "get_conversion_factor_to_base_units" method """
        self.assertEqual(u.meter.get_conversion_factor_to_base_units(), 1)
        self.assertEqual(u.calorie.get_conversion_factor_to_base_units(), 4.184)
        kcpma = u.md_kilocalorie/u.mole/u.angstrom
        self.assertEqual(kcpma.get_conversion_factor_to_base_units(), 4.184)

    def testScalarQuantityMultiplyDivide(self):
        """ Tests creating a scalar Quantity object by * or / by a Unit """
        self.assertTrue(u.is_quantity(5 * u.centimeters))
        self.assertTrue(u.is_quantity(1 / u.centimeters))
        self.assertTrue(u.is_quantity(10 * u.centimeters))
        self.assertTrue(u.is_quantity(9.81 * u.meters / u.second**2))

    def testScalarQuantityConstructor(self):
        """ Tests creating a Quantity using the Quantity constructor """
        self.assertTrue(u.is_quantity(u.Quantity(5, u.centimeters)))
        self.assertTrue(u.is_quantity(u.Quantity(5, u.centimeters**-1)))
        x = u.Quantity(value=5.0, unit=100.0*u.meters)
        self.assertTrue(u.is_quantity(x))
        self.assertEqual(x, 500*u.meters)

    def testValueInUnit(self):
        """ Tests the value_in_unit functionality for Quantity """
        i = 5 * u.centimeters
        self.assertEqual(i.value_in_unit(u.millimeters), 50)
        self.assertEqual(i / u.millimeters, 50)

    def testCollectionQuantities(self):
        """ Tests the use of collections as Quantity values """
        s = [1, 2, 3] * u.centimeters
        self.assertEqual(str(s), '[1, 2, 3] cm')
        self.assertTrue(u.is_quantity(s))
        s2 = s / u.millimeters
        self.assertEqual(s2, [10.0, 20.0, 30.0])
        self.assertEqual(s2, s.value_in_unit(u.millimeters))
        # Test 2-D list
        s = [[1, 2, 3], [4, 5, 6]]
        s *= u.centimeters
        self.assertTrue(u.is_quantity(s))
        s2 = s / u.millimeters
        self.assertEqual(s2, [[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]])
        self.assertEqual(s.value_in_unit(u.millimeters), s2)
        # Test tuples
        s = (1, 2, 3) * u.centimeters
        self.assertTrue(u.is_quantity(s))
        self.assertEqual(str(s), '(1, 2, 3) cm')
        s2 = s / u.millimeters
        self.assertEqual(s2, (10, 20, 30))
        self.assertIsInstance(s2, tuple)
        self.assertEqual(s.value_in_unit(u.millimeters), s2)
        self.assertIsInstance(s.value_in_unit(u.millimeters), tuple)
        x = [1, 2, 3] * u.centimeters
        x *= u.meters
        self.assertEqual(x, [100, 200, 300] * u.centimeters**2)

    def testReduceUnit(self):
        """ Tests the reduce_unit functionality """
        x = u.nanometer**2 / u.angstrom**2
        self.assertEqual(str(x), 'nanometer**2/(angstrom**2)')
        self.assertTrue(x.is_dimensionless())
        q = u.Quantity(2.0, x)
        self.assertEqual(str(q), '2.0 nm**2/(A**2)')
        self.assertEqual(q.reduce_unit(), 200)

    def testCollectionQuantityOperations(self):
        """ Tests that Quantity collections behave correctly """
        # Tests that __getitem__ returns a unit
        s = [1, 2, 3, 4] * u.angstroms
        self.assertTrue(u.is_quantity(s[0]))
        for i, val in enumerate(s):
            self.assertTrue(u.is_quantity(val))
            self.assertEqual(val, (i+1) * u.angstroms)
        # Tests that __setitem__ fails when an incompatible type is added
        def fail(s): s[0] = 5
        self.assertRaises(AttributeError, lambda: fail(s))
        def fail(s): s[0] = 5 * u.joules
        self.assertRaises(TypeError, lambda: fail(s))
        def fail(s): s[0] /= 10 * u.meters
        self.assertRaises(AttributeError, lambda: fail(s))
        # Tests that __setitem__ converts to the unit of the container
        s[0] = 1 * u.nanometers
        self.assertEqual(s[0]._value, 10)
        # Tests that __setitem__ handles slice assignment correctly
        x = [0, 1, 2, 3, 4] * u.kelvin
        x[2:4] = [-2, -3] * u.kelvin
        self.assertEqual(x._value, [0, 1, -2, -3, 4])
        # Tests standard unit conversions
        x = [1, 2, 3] * u.centimeters
        self.assertEqual(x / u.millimeters, [10, 20, 30])
        # Test the construction of a container in which each element is a
        # Quantity, passed to the Quantity constructor
        x = u.Quantity([1*u.angstrom, 2*u.nanometer, 3*u.angstrom])
        self.assertEqual(x._value, [1, 20, 3])
        self.assertEqual(x.unit, u.angstrom)
        x = u.Quantity((1, 2, 3))
        self.assertTrue(u.is_quantity(x))
        self.assertTrue(x.unit.is_dimensionless())
        x = u.Quantity(([1*u.angstrom, 2*u.nanometer, 3*u.angstrom],
                        [1*u.angstrom, 4*u.nanometer, 3*u.angstrom]))
        self.assertEqual(x._value, ([1, 20, 3], [1, 40, 3]))
        self.assertEqual(x.unit, u.angstrom)
        self.assertTrue(u.is_quantity(u.Quantity([])))

    def testMutableQuantityOperations(self):
        " Tests that mutable Quantity objects do not get unexpectedly changed "
        # This used to be a bug -- t and s._value were the same object, so
        # changing t would also change s silently
        s = [1, 2, 3, 4] * u.angstroms
        t = s / u.angstroms
        self.assertEqual(t, [1, 2, 3, 4])
        self.assertEqual(s, u.Quantity([1, 2, 3, 4], u.angstroms))
        t[0] = 2
        self.assertEqual(t, [2, 2, 3, 4])
        self.assertEqual(s, u.Quantity([1, 2, 3, 4], u.angstroms))
        t = s.value_in_unit(u.angstroms)
        self.assertEqual(t, [1, 2, 3, 4])
        self.assertEqual(s, u.Quantity([1, 2, 3, 4], u.angstroms))
        t[0] = 2
        self.assertEqual(t, [2, 2, 3, 4])
        self.assertEqual(s, u.Quantity([1, 2, 3, 4], u.angstroms))
        s = [1, 2, 3, 4] * u.nanometers
        self.assertEqual(s, u.Quantity([1, 2, 3, 4], u.nanometers))
        t = s.in_units_of(u.nanometers)
        self.assertEqual(s, t)
        t[0] = 1 * u.meters
        self.assertAlmostEqualQuantities(t,
                u.Quantity([1e9, 2, 3, 4], u.nanometers))
        self.assertEqual(s, u.Quantity([1, 2, 3, 4], u.nanometers))

    def testQuantityMaths(self):
        """ Tests dimensional analysis & maths on and b/w Quantity objects """
        x = 1.3 * u.meters
        y = 75.2 * u.centimeters
        self.assertEqual((x + y) / u.meters, 2.052)
        self.assertEqual((x - y) / u.meters, 0.548)
        self.assertEqual(x / y, 1.3 / 0.752)
        self.assertEqual(x * y, 1.3*0.752*u.meters**2)
        d1 = 2.0*u.meters
        d2 = 2.0*u.nanometers
        self.assertEqual(d1 + d2, (2+2e-9)*u.meters)
        self.assertAlmostEqual((d2+d1-(2e9+2)*u.nanometers)._value, 0, places=6)
        self.assertEqual(d1 + d1, 4.0*u.meters)
        self.assertEqual(d1 - d1, 0.0*u.meters)
        self.assertEqual(d1 / d1, 1.0)
        self.assertEqual(d1 * u.meters, 2.0*u.meters**2)
        self.assertEqual(u.kilograms*(d1/u.seconds)*(d1/u.seconds),
                         4*u.kilograms*u.meters**2/u.seconds**2)
        self.assertEqual(u.kilograms*(d1/u.seconds)**2,
                         4*u.kilograms*u.meters**2/u.seconds**2)
        self.assertEqual(d1**3, 8.0*u.meters**3)
        x = d1**(3/2)
        self.assertAlmostEqual(x._value, math.sqrt(2)**3)
        self.assertEqual(x.unit, u.meters**(3/2))
        self.assertAlmostEqual((d1**0.5)._value, math.sqrt(2))
        self.assertEqual((d1**0.5).unit, u.meters**0.5)
        comp = (3.0 + 4.0j) * u.meters
        self.assertTrue(u.is_quantity(comp))
        self.assertEqual(comp.unit, u.meters)
        self.assertEqual(str(comp), '(3+4j) m')
        self.assertEqual(comp + comp, (6.0 + 8.0j)*u.meters)
        self.assertEqual(comp - comp, 0*u.meters)
        self.assertEqual(comp * comp, (3.0 + 4.0j)**2 * u.meters**2)
        self.assertAlmostEqual(comp / comp, 1)
        self.assertAlmostEqual(1.5*u.nanometers / u.meters, 1.5e-9, places=15)
        self.assertEqual((2.3*u.meters)**2, 2.3**2*u.meters**2)
        x = 4.3 * u.meters
        self.assertEqual(x / u.centimeters, 430)
        self.assertEqual(str(x / u.seconds), '4.3 m/s')
        self.assertEqual(str(8.4 / (4.2*u.centimeters)), '2.0 /cm')
        x = 1.2 * u.meters
        self.assertEqual(x * 5, u.Quantity(6.0, u.meters))

    def testQuantityComplicatedMaths(self):
        """ Tests a complicated set of mathematical operations on a Quantity """
        s1 = 2.0
        x1 = 2
        x2 = 4.0 / 3.0
        u1 = u.kilogram * u.meter / u.second**2
        u2 = u1 * u.meter
        q1 = 1.0 * u1
        q2 = 2.0 * u2
        self.assertEqual(s1, 2.0)
        self.assertEqual(x1, 2)
        self.assertAlmostEqual(x2, 1.33333333333333)
        self.assertEqual(str(u1), 'kilogram*meter/(second**2)')
        self.assertEqual(str(u2), 'kilogram*meter**2/(second**2)')
        self.assertEqual(str(q1), '1.0 kg m/(s**2)')
        self.assertEqual(str(q2), '2.0 kg m**2/(s**2)')
        self.assertEqual(str(u1*s1), '2.0 kg m/(s**2)')
        self.assertEqual(str(s1*u1), '2.0 kg m/(s**2)')
        self.assertEqual(str(u1/s1), '0.5 kg m/(s**2)')
        self.assertEqual(str(s1/u1), '2.0 s**2/(kg m)')
        self.assertEqual(str(u1*u1), 'kilogram**2*meter**2/(second**4)')
        self.assertEqual(u1/u1, u.dimensionless)
        self.assertEqual(str(u1/u1), 'dimensionless')
        self.assertEqual(str(u1*u2), 'kilogram**2*meter**3/(second**4)')
        self.assertEqual(u1/u2, u.meters**-1)
        self.assertEqual(str(u1/u2), '/meter')
        self.assertEqual(u1**x1, u.kilogram**2*u.meter**2/(u.second**4))
        self.assertEqual(u1**(1/x1), u.kilogram**0.5*u.meter**0.5/u.second)
        self.assertEqual(u1**x2,
                         u.kilogram**(1+1/3)*u.meter**(1+1/3)/u.second**(2+2/3))
        q = 1.0 * u.md_kilocalorie/u.mole/u.angstrom
        self.assertEqual(str(q.in_units_of(u.md_kilojoule/u.mole/u.nanometer)),
                         '41.84 kJ/(nm mol)')

    def testQuantityComparisons(self):
        """ Tests binary comparison operators between Quantity """
        l1 = 1.0 * u.meters
        l2 = 2.0 * u.meters
        l3 = 1.0 * u.meters
        self.assertEqual(l1, l3)
        self.assertNotEqual(l1, l2)
        self.assertLessEqual(l1, l2)
        self.assertLess(l1, l2)
        self.assertGreater(l2, l1)
        self.assertGreaterEqual(l2, l1)
        self.assertLessEqual(l1, l3)
        self.assertGreaterEqual(l1, l3)
        self.assertGreater(1.2*u.meters, 72*u.centimeters)
        self.assertLess(1*u.meters, 200*u.centimeters)

    def testChemistryProblems(self):
        """ Tests some gen-chem applications with Quantity's """
        def work(f, dx):
            return f * dx

        F = 1.0 * u.kilogram * u.meter / u.second**2
        dx = 1.0 * u.meter
        self.assertEqual(work(F, dx), 1.0 * u.joule)
        self.assertEqual(F, 1.0 * u.newton)

        def ideal_gas_law(P, V, T):
            R = u.MOLAR_GAS_CONSTANT_R
            return (P * V / (R * T)).in_units_of(u.mole)

        T = (273.0 + 37.0) * u.kelvin
        P = (1.01325e5) * u.pascals
        r = 0.5e-6 * u.meters
        V = 4/3 * math.pi * r**3
        n = ideal_gas_law(P, V, T)
        val = 4/3*math.pi*0.5e-6**3*1
        self.assertAlmostEqualQuantities(P*V, val * u.atmospheres*u.meters**3)
        self.assertAlmostEqualQuantities(n, 2.05834818672e-17 * u.mole)
        self.assertAlmostEqualQuantities(V, 5.2359833333333e-19 * u.meters**3)
        self.assertEqual(str(T), '310.0 K')
        self.assertEqual(str(1*u.joules/u.kelvin/u.mole), '1 J/(K mol)')
        self.assertTrue(u.is_quantity(V))
        # Checks trouble with complicated unit conversion factors
        p1 = 1.0 * u.atmospheres
        p2 = p1.in_units_of(u.joules/u.nanometers**3)
        V = 2.4 * u.nanometers**3
        beta = 4.e-4 * u.mole/u.joule
        x1 = beta * p1 * V
        y1 = x1 * u.AVOGADRO_CONSTANT_NA
        self.assertAlmostEqual(y1, 0.0585785776197)
        x2 = beta * p2 * V
        y2 = x2 * u.AVOGADRO_CONSTANT_NA
        self.assertAlmostEqual(y1, y2)

    def testAngleQuantities(self):
        """ Tests angle measurements """
        self.assertEqual(1.0*u.radians / u.degrees, 180 / math.pi)
        self.assertTrue(u.is_quantity(1.0*u.radians))
        self.assertTrue(u.is_quantity(1.0*u.degrees))
        self.assertEqual((1.0*u.radians).in_units_of(u.degrees),
                         (180 / math.pi)*u.degrees)
        self.assertEqual(90*u.degrees/u.radians, math.pi/2)
        q = 90 * u.degrees + 0.3 * u.radians
        self.assertEqual(q._value, 90 + 180*0.3/math.pi)
        self.assertEqual(q.unit, u.degrees)

    def testIsQuantity(self):
        """ Tests if is_quantity can detect Quantities vs. scalars and units """
        self.assertTrue(u.is_quantity(1/u.second))
        self.assertFalse(u.is_quantity(u.second**-1))
        self.assertTrue(u.is_quantity(10*u.meters))
        self.assertFalse(u.is_quantity(10))

    def testBadQuantityMaths(self):
        """ Tests that Maths of incompatible units properly fails """
        self.assertRaises(TypeError, lambda:1*u.meters + 1*u.liters)
        self.assertRaises(AttributeError, lambda: 1*u.liters + 5)

    def testBadQuantityComparisons(self):
        """ Checks that comparisons of incompatible units fails """
        self.assertRaises(TypeError, lambda: 1*u.meters > 1*u.liters)
        self.assertRaises(TypeError, lambda: 1*u.nanometers > 1*u.degrees)

    def testDimensionless(self):
        """ Tests the properties of unit.dimensionless """
        x = 5 * u.dimensionless
        y = u.Quantity(5, u.dimensionless)
        self.assertTrue(u.is_quantity(x))
        self.assertTrue(u.is_quantity(y))
        self.assertNotEqual(x, 5)
        self.assertNotEqual(y, 5)
        self.assertEqual(x, y)
        self.assertEqual(x.value_in_unit_system(u.si_unit_system), 5)
        self.assertEqual(x.value_in_unit_system(u.cgs_unit_system), 5)
        self.assertEqual(x.value_in_unit_system(u.md_unit_system), 5)
        x = u.Quantity(1.0, u.dimensionless)
        y = u.Quantity(1.0, u.dimensionless)
        self.assertIsNot(x, y)
        self.assertEqual(x, y)

    def testTruthiness(self):
        """ Tests the truthiness of units """
        true = 2.3 * u.meters
        false = 0.0 * u.meters
        self.assertTrue(true)
        self.assertFalse(false)
        self.assertTrue(bool(true))
        self.assertFalse(bool(false))

    def testUnaryOperators(self):
        """ Tests unary operators on units """
        self.assertEqual(-(2.3*u.meters), u.Quantity(-2.3, u.meters))
        self.assertEqual(-(2.3*u.meters), -u.Quantity(2.3, u.meters))
        self.assertEqual(+(2.3*u.meters), u.Quantity(2.3, u.meters))
        self.assertEqual(2.3*u.meters, +u.Quantity(2.3, u.meters))

    def testUnitMathModule(self):
        """ Tests the unit_math functions on Quantity objects """
        self.assertEqual(u.sqrt(1.0*u.kilogram*u.joule),
                         1.0*u.kilogram*u.meter/u.second)
        self.assertEqual(u.sqrt(1.0*u.kilogram*u.calorie),
                         math.sqrt(4.184)*u.kilogram*u.meter/u.second)
        self.assertEqual(u.sqrt(9), 3) # Test on a scalar
        self.assertEqual(u.sin(90*u.degrees), 1)
        self.assertEqual(u.sin(math.pi/2*u.radians), 1)
        self.assertEqual(u.sin(math.pi/2), 1)
        self.assertEqual(u.cos(180*u.degrees), -1)
        self.assertEqual(u.cos(math.pi*u.radians), -1)
        self.assertEqual(u.cos(math.pi), -1)
        self.assertAlmostEqual(u.tan(45*u.degrees), 1)
        self.assertAlmostEqual(u.tan(math.pi/4*u.radians), 1)
        self.assertAlmostEqual(u.tan(math.pi/4), 1)
        acos = u.acos(1.0)
        asin = u.asin(1.0)
        atan = u.atan(1.0)
        self.assertTrue(u.is_quantity(acos))
        self.assertTrue(u.is_quantity(asin))
        self.assertTrue(u.is_quantity(atan))
        self.assertEqual(acos.unit, u.radians)
        self.assertEqual(asin.unit, u.radians)
        self.assertEqual(atan.unit, u.radians)
        self.assertEqual(acos.value_in_unit(u.degrees), 0)
        self.assertEqual(acos / u.radians, 0)
        self.assertEqual(asin.value_in_unit(u.degrees), 90)
        self.assertEqual(asin / u.radians, math.pi/2)
        self.assertAlmostEqual(atan.value_in_unit(u.degrees), 45)
        self.assertAlmostEqual(atan / u.radians, math.pi/4)
        # Check some sequence maths
        seq = [1, 2, 3, 4] * u.meters
        self.assertEqual(u.sum(seq), 10*u.meters)
        self.assertEqual(u.dot(seq, seq), (1+4+9+16)*u.meters**2)
        self.assertEqual(u.norm(seq), math.sqrt(30)*u.meters)

    def testUnitMathModuleBadInput(self):
        """ Tests that bad units to unit_math fails appropriately """
        self.assertRaises(ArithmeticError, lambda: u.sqrt(9*u.meters))
        self.assertRaises(TypeError, lambda: u.sin(1*u.meters))

    def testBadQuantityConversions(self):
        """ Tests that conversions to incompatible units fails """
        self.assertRaises(TypeError,
                lambda: (1.0*u.meters).in_units_of(u.seconds))
        self.assertRaises(TypeError,
                lambda: (1.0*u.meters).value_in_unit(u.seconds))

    def testCreateNewBaseUnit(self):
        """ Tests creating a new base unit """
        ms = u.milli * u.second_base_unit
        self.assertIsInstance(ms, u.BaseUnit)
        self.assertEqual(repr(ms),
                'BaseUnit(base_dim=BaseDimension("time"), name="millisecond", '
                'symbol="ms")')
        self.assertEqual(ms.conversion_factor_to(u.second_base_unit), 0.001)

    def testCreateNewScaledUnit(self):
        """ Tests creating a new ScaledUnit """
        mC = u.milli * u.ScaledUnit(4.184, u.joule, "calorie", "cal")
        self.assertIsInstance(mC, u.ScaledUnit)
        self.assertEqual(repr(mC),
                "ScaledUnit(factor=0.004184, master=joule, name='millicalorie',"
                " symbol='mcal')")
        self.assertFalse(u.is_unit(mC))

    def testCreateNewUnit(self):
        """ Tests creating a new Unit """
        ms = u.milli * u.second
        self.assertTrue(u.is_unit(ms))
        self.assertEqual(repr(ms),
                'Unit({BaseUnit(base_dim=BaseDimension("time"), '
                'name="millisecond", symbol="ms"): 1.0})')

    def testIllegalQuantityOps(self):
        """ Test that illegal operations on Quantity objects fails """
        self.assertRaises(TypeError, lambda: u.milli * (1.0 * u.second))

    def testQuantityFormat(self):
        """ Tests the format method on Quantity instances """
        x = 5.439999999 * u.picosecond
        self.assertEqual(x.format("%.3f"), '5.440 ps')

    def testQuantityCollectionMethods(self):
        """ Tests some sequence methods on Quantity objects (no numpy) """
        seq = [1, 2, 3, 4] * u.meters
        self.assertEqual(seq.sum(), 10*u.meters)
        self.assertEqual(seq.mean(), 2.5*u.meters)
        self.assertEqual(seq.max(), 4*u.meters)
        self.assertEqual(seq.min(), 1*u.meters)
        self.assertAlmostEqual(seq.std(), 1.1180339887498949*u.meters)

    def testString(self):
        """ Tests unit handling with strings, which should be dimensionless """
        s = u.Quantity("string")
        self.assertEqual(s.value_in_unit_system(u.md_unit_system), "string")

    def testMisc(self):
        """ Miscellaneous tests for the unit package """
        self.assertTrue(u.meter is not None)
        self.assertFalse(u.meter is None)
        self.assertEqual(repr(1.2*u.meters), 'Quantity(value=1.2, unit=meter)')
        class Foo(object):
            def bar(self):
                return 'bar'
        x = Foo()
        self.assertEqual(x.bar(), 'bar')
        y = x * u.nanometers
        self.assertEqual(y.bar(), 'bar')
        self.assertEqual(str(u.meters*u.centimeters), 'centimeter*meter')
        self.assertEqual(str(u.meters*u.meters), 'meter**2')
        self.assertEqual(str(u.meter*u.meter), 'meter**2')

@unittest.skipIf(np is None, 'Skipping numpy units tests')
class TestNumpyUnits(QuantityTestCase):

    def testNumpyQuantity(self):
        """ Tests that numpy arrays can form Quantity values """
        q = u.Quantity(np.array([1, 2, 3]), u.centimeters)
        self.assertTrue(u.is_quantity(q))
        self.assertIsInstance(q._value, np.ndarray)
        self.assertTrue(np.all(q / u.millimeters == np.array([1, 2, 3])*10))

    def testNumpyDeepCopy(self):
        """ Check that deepcopy on numpy array does not strip units """
        x = u.Quantity(np.zeros((2, 3)), u.nanometer)
        y = copy.deepcopy(x)
        self.assertTrue(np.all(x == y))
        self.assertTrue(u.is_quantity(x))
        self.assertTrue(u.is_quantity(y))

    def testNumpyDivision(self):
        """ Tests that division of numpy Quantities works correctly """
        x = u.Quantity(np.asarray([1., 2.]), u.nanometers)
        y = u.Quantity(np.asarray([3., 4.]), u.picoseconds)
        xy = x / y
        self.assertTrue(u.is_quantity(xy))
        self.assertEqual(xy.unit, u.nanometers/u.picoseconds)
        self.assertEqual(xy[0].value_in_unit(u.nanometers/u.picoseconds), 1/3)
        self.assertEqual(xy[1].value_in_unit(u.nanometers/u.picoseconds), 0.5)

    def testNumpyIsString(self):
        """ Tests the internal _is_string method with numpy Quantities """
        from simtk.unit.quantity import _is_string
        a = np.array([[1, 2, 3], [4, 5, 6]])
        self.assertIsInstance("", str)
        self.assertTrue(_is_string(""))
        self.assertTrue(_is_string("t"))
        self.assertTrue(_is_string("test"))
        self.assertFalse(_is_string(3))
        self.assertFalse(_is_string(a))

    def testNumpyFunctions(self):
        """ Tests various numpy attributes that they result in Quantities """
        a = u.Quantity(np.arange(10), u.seconds)
        self.assertEqual(a.max(), 9*u.seconds)
        self.assertEqual(a.min(), 0*u.seconds)
        self.assertEqual(a.mean(), 4.5*u.seconds)
        self.assertAlmostEqualQuantities(a.std(), 2.8722813232690143*u.seconds)
        b = a.reshape((5, 2))
        self.assertTrue(u.is_quantity(b))

    def testMultiplication(self):
        """ Tests that units override numpy.ndarray multiplication """
        self.assertIsInstance(np.arange(10)*u.angstroms, u.Quantity)
        # This only works with versions of numpy > 1.7 due to a bug in older
        # versions. Since Travis-CI installs Python 1.6.1 from aptitude, and we
        # don't want it to report test failures *every time*, just disable this
        # particular test for the numpy versions known to be bad.
        if np.version.version > '1.7':
            x = np.array([1]) * u.liters
            self.assertIsInstance(x, u.Quantity)
            self.assertIsInstance(np.arange(10) * x, u.Quantity)
